# MGGG Replicate Docker Information

This repository is meant to serve as an archive containing the information necessary to
rebuild the Docker image that we use as a part of our [gerrytools](https://github.com/mggg/gerrytools)
package. Specifically, the dockerfile in this repository can be used to reconstruct the
image [mgggdev/replicate](https://hub.docker.com/repository/docker/mgggdev/replicate/general) 
published on Docker Hub. We use this image as a part of the `ben` and `mgrp` submodules
in `gerrytools`.

In the event that you would like to build the image yourself, simply run the command


```
docker build -t <image-name> .
```

from the root of this cloned repository (this will take ~30 minutes). You may then specify 
your custom image in the `gerrytools.mgrp` and `gerrytools.ben` api by using


```python 
with ReplicatorContainer(
    configuration=my_configuration,
    docker_image_name="<image-name>", 
) as rep:
```

or 

```python 
with BenContainer(
    docker_image_name="<image-name>", 
) as bc:
```

## What's in the /home folder

The `/home/` folder of this repository is copied into the home folder of the docker image
that is generated by the dockerfile. Within the home folder, you will find the cli that
are used as a part of the `gerrytools.mgrp` module as well as the source code for some
packages that we are working with. Specifically, the code that is invoked by the 
`gerrytools.mgrp.recom.py` and `gerrytools.mgrp.forest.py` are contained in these folders,
and these are, in turn, just clones of the [frcw](https://github.com/peterrrock2/frcw) and
[Multi Scale Map Sampler](https://github.com/pterrrock2/Multi-Scale-Map-Sampler). 


### Last update

The versions of the code in the `/home` folder were last updated on 5 April 2024. If these
need to be updated for some reason, please open an issue.
