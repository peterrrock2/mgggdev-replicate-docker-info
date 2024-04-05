//! Main CLI for frcw.
use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use clap::{value_t, App, Arg};
use frcw::config::parse_region_weights_config;
use frcw::init::from_networkx;
use frcw::recom::run::multi_chain;
use frcw::recom::{RecomParams, RecomVariant};
use frcw::stats::{
    AssignmentsOnlyWriter, BenWriter, CanonicalWriter, JSONLWriter, PcompressWriter, StatsWriter,
    TSVWriter,
};
use serde_json::json;
use sha3::{Digest, Sha3_256};
use std::path::PathBuf;
use std::{fs, io};

fn main() {
    let mut cli = App::new("frcw")
        .version("0.1.0")
        .author("Parker J. Rule <parker.rule@tufts.edu>")
        .about("A minimal implementation of the ReCom Markov chain")
        .arg(
            Arg::with_name("graph_json")
                .long("graph-json")
                .takes_value(true)
                .required(true)
                .help("The path of the dual graph (in NetworkX format)."),
        )
        .arg(
            Arg::with_name("n_steps")
                .long("n-steps")
                .takes_value(true)
                .required(true)
                .help("The number of proposals to generate."),
        )
        .arg(
            Arg::with_name("tol")
                .long("tol")
                .takes_value(true)
                .required(true)
                .help("The relative population tolerance."),
        )
        .arg(
            Arg::with_name("pop_col")
                .long("pop-col")
                .takes_value(true)
                .required(true)
                .help("The name of the total population column in the graph metadata."),
        )
        .arg(
            Arg::with_name("assignment_col")
                .long("assignment-col")
                .takes_value(true)
                .required(true)
                .help("The name of the assignment column in the graph metadata."),
        )
        .arg(
            Arg::with_name("rng_seed")
                .long("rng-seed")
                .takes_value(true)
                .required(true)
                .help("The seed of the RNG used to draw proposals."),
        )
        .arg(
            Arg::with_name("balance_ub")
                .long("balance-ub")
                .short("M") // Variable used in RevReCom paper
                .takes_value(true)
                .default_value("0") // TODO: just use unwrap_or_default() instead?
                .help("The normalizing constant (reversible ReCom only)."),
        )
        .arg(
            Arg::with_name("n_threads")
                .long("n-threads")
                .takes_value(true)
                .required(true)
                .help("The number of threads to use."),
        )
        .arg(
            Arg::with_name("batch_size")
                .long("batch-size")
                .takes_value(true)
                .required(true)
                .help("The number of proposals per batch job."),
        )
        .arg(
            Arg::with_name("variant")
                .long("variant")
                .takes_value(true)
                .default_value("reversible"),
        ) // other options: cut_edges, district_pairs
        .arg(
            Arg::with_name("writer")
                .long("writer")
                .takes_value(true)
                .default_value("jsonl"),
        ) // other options: jsonl-full, tsv
        .arg(
            Arg::with_name("sum_cols")
                .long("sum-cols")
                .multiple(true)
                .takes_value(true),
        )
        .arg(
            Arg::with_name("region_weights")
                .long("region-weights")
                .takes_value(true)
                .help("Region columns with weights for region-aware ReCom."),
        )
        .arg(Arg::with_name("cut_edges_count").long("cut-edges-count"))
        .arg(
            Arg::with_name("output-file")
                .long("output-file")
                .short("o")
                .takes_value(true)
                .help("The path to write the output to."),
        );
    if cfg!(feature = "linalg") {
        cli = cli.arg(Arg::with_name("spanning_tree_counts").long("st-counts"));
    }
    let matches = cli.get_matches();
    let n_steps = value_t!(matches.value_of("n_steps"), u64).unwrap_or_else(|e| e.exit());
    let rng_seed = value_t!(matches.value_of("rng_seed"), u64).unwrap_or_else(|e| e.exit());
    let tol = value_t!(matches.value_of("tol"), f64).unwrap_or_else(|e| e.exit());
    let balance_ub = value_t!(matches.value_of("balance_ub"), u32).unwrap_or_else(|e| e.exit());
    let n_threads = value_t!(matches.value_of("n_threads"), usize).unwrap_or_else(|e| e.exit());
    let batch_size = value_t!(matches.value_of("batch_size"), usize).unwrap_or_else(|e| e.exit());
    let graph_json = fs::canonicalize(PathBuf::from(matches.value_of("graph_json").unwrap()))
        .unwrap()
        .into_os_string()
        .into_string()
        .unwrap();
    let pop_col = matches.value_of("pop_col").unwrap();
    let assignment_col = matches.value_of("assignment_col").unwrap();
    let variant_str = matches.value_of("variant").unwrap();
    let writer_str = matches.value_of("writer").unwrap();
    let st_counts = matches.is_present("spanning_tree_counts");
    let cut_edges_count = matches.is_present("cut_edges_count");
    let mut sum_cols: Vec<String> = matches
        .values_of("sum_cols")
        .unwrap_or_default()
        .map(|c| c.to_string())
        .collect();
    let region_weights_raw = matches.value_of("region_weights").unwrap_or_default();

    let variant = match variant_str {
        "reversible" => RecomVariant::Reversible,
        "cut-edges-ust" => RecomVariant::CutEdgesUST,
        "cut-edges-rmst" => RecomVariant::CutEdgesRMST,
        "cut-edges-region-aware" => RecomVariant::CutEdgesRegionAware,
        "district-pairs-ust" => RecomVariant::DistrictPairsUST,
        "district-pairs-rmst" => RecomVariant::DistrictPairsRMST,
        "district-pairs-region-aware" => RecomVariant::DistrictPairsRegionAware,
        bad => panic!("Parameter error: invalid variant '{}'", bad),
    };

    let output_buffer: Box<dyn io::Write + Send> = match matches.value_of("output-file") {
        Some(path) => {
            let path = std::path::Path::new(path);
            if path.exists() {
                panic!("Output file already exists.");
            };
            Box::new(io::BufWriter::new(fs::File::create(path).unwrap()))
        }
        None => Box::new(io::BufWriter::new(std::io::stdout())),
    };

    let writer: Box<dyn StatsWriter> = match writer_str {
        "tsv" => Box::new(TSVWriter::new(output_buffer)),
        "jsonl" => Box::new(JSONLWriter::new(
            false,
            st_counts,
            cut_edges_count,
            output_buffer,
        )),
        "pcompress" => Box::new(PcompressWriter::new(output_buffer)),
        "jsonl-full" => Box::new(JSONLWriter::new(
            true,
            st_counts,
            cut_edges_count,
            output_buffer,
        )),
        "assignments" => Box::new(AssignmentsOnlyWriter::new(false, output_buffer)),
        "canonicalized-assignments" => Box::new(AssignmentsOnlyWriter::new(true, output_buffer)),
        "canonical" => Box::new(CanonicalWriter::new(output_buffer)),
        "ben" => Box::new(BenWriter::new(output_buffer)),
        bad => panic!("Parameter error: invalid writer '{}'", bad),
    };
    if variant == RecomVariant::Reversible && balance_ub == 0 {
        panic!("For reversible ReCom, specify M > 0.");
    }

    assert!(tol >= 0.0 && tol <= 1.0);

    let region_weights = parse_region_weights_config(region_weights_raw);
    // Add the keys in the region weights to sum_cols if they are not there already
    // so that the user doesn't have to
    if let Some(weight_pairs_vec) = &region_weights {
        for (key, _) in weight_pairs_vec.iter() {
            if !sum_cols.contains(&key) {
                sum_cols.push(key.clone().to_string());
            }
        }
    }

    let (graph, partition) = from_networkx(&graph_json, pop_col, assignment_col, sum_cols).unwrap();
    let avg_pop = (graph.total_pop as f64) / (partition.num_dists as f64);

    let params = RecomParams {
        min_pop: ((1.0 - tol) * avg_pop as f64).floor() as u32,
        max_pop: ((1.0 + tol) * avg_pop as f64).ceil() as u32,
        num_steps: n_steps,
        rng_seed: rng_seed,
        balance_ub: balance_ub,
        variant: variant,
        region_weights: region_weights.clone(),
    };

    let mut graph_file = fs::File::open(&graph_json).unwrap();
    let mut graph_hasher = Sha3_256::new();
    io::copy(&mut graph_file, &mut graph_hasher).unwrap();
    let graph_hash = format!("{:x}", graph_hasher.finalize());
    let mut meta = json!({
        "assignment_col": assignment_col,
        "tol": tol,
        "pop_col": pop_col,
        "graph_path": graph_json,
        "graph_sha3": graph_hash,
        "batch_size": batch_size,
        "rng_seed": rng_seed,
        "num_threads": n_threads,
        "num_steps": n_steps,
        "parallel": true,
        "graph_json": graph_json,
        "chain_variant": variant_str,
    });
    if variant == RecomVariant::Reversible {
        meta.as_object_mut()
            .unwrap()
            .insert("balance_ub".to_string(), json!(balance_ub));
    }
    if region_weights.is_some() {
        meta.as_object_mut()
            .unwrap()
            .insert("region_weights".to_string(), json!(region_weights));
    }
    if writer_str == "jsonl" || writer_str == "jsonl-full" {
        // hotfix for pcompress writing
        // TODO: move this into init
        println!("{}", json!({ "meta": meta }).to_string());
    }
    multi_chain(&graph, &partition, writer, &params, n_threads, batch_size);
}
