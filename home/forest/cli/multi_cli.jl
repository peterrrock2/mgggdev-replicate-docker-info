using ArgParse

include("run_sampler.jl")

parser = ArgParseSettings()

@add_arg_table! parser begin
    "--input-file-name"
        help = "Name of input file (rest of path assumed)"
    "--output-file-name"
        help = "Name of output file (rest of path assumed)"
    "--subregion-name"
        help = "Label for the subregion column"
    "--region-name"
        help = "Label for the region column"
    "--pop-name"
        help = "Label for the population column"
    "--num-dists"
        help = "Number of districts"
        arg_type = Int
    "--rng-seed"
        help = "Seed for the rng"
        arg_type = Int
    "--pop-dev"
        help = "Allowable population deviance (between 0 and 1)"
        arg_type = Float64
    "--gamma"
        help = "Value for gamma in the multiscale"
        arg_type = Float64
    "--steps"
        help = "Number of steps allowed"
        arg_type = Int
end

args = parse_args(parser)


if args["output-file-name"] == nothing
    output_path = nothing
else
    output_path = Some(args["output-file-name"])
end

run_multiscale2(
    pctGraphPath=args["input-file-name"],
    subregion_name=args["subregion-name"],
    region_name=args["region-name"],
    population_col=args["pop-name"],
    output_path=output_path,
    num_dists=args["num-dists"],
    rng_seed=args["rng-seed"],
    pop_dev=args["pop-dev"],
    gamma=args["gamma"],
    steps=args["steps"]
)
