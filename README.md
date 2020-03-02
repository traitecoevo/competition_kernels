# Competition kernels

Code for the competition kernels.


Additional accessory files are also including:

- `DECRIPTION`: A machine-readable [compendium]() file containing key metadata and dependencies 
- `LICENSE`: License for the materials
- `Dockerfile` & `.binder/Dockerfile`: files used to generate docker containers for long-term reproducibility

## Running the code

All analyses were done in `R`, and the paper is written in LaTeX. All code needed to reproduce the submitted products is included in this repository. To reproduce this paper, run the code contained in the `analysis.R` file. Figures will be output to a directory called `output` and the paper and supplementary materials in the folder `ms`.


The paper was written in 2016 using a version of R available at the time. With some minor updates, the code has been updated and was last seen running wild and free on R 3.6.1. You can try running it on your current version and it may work. 

To ensure [computational reproducibility](https://www.britishecologicalsociety.org/wp-content/uploads/2017/12/guide-to-reproducible-code.pdf) into the future, we have also generated [Docker](http://dockerhub.com) and [Binder](https://mybinder.org) containers, enabling you to launch a compute environment built off R 3.6.1 with all the dependencies included.

### Running locally

If reproducing these results on your own machine, first download the code and then install the required packages, listed under `Depends` in the `DESCRIPTION` file. This can be achieved by opening the Rstudio project and running:

```{r}
#install.packages("devtools")
devtools::install_deps()
```

Then run `analysis.R`. 

### Running on Binder 

You can launch the analysis on the web in an interactive RStudio session with the required software pre-installed. This session is hosted by binder and can be accessed by clicking on the following:

[![Launch Rstudio Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/traitecoevo/competition_kernels/master?urlpath=rstudio)

### Running via Docker

If you have Docker installed, you can recreate the compute environment as follows. 

First fetch the container:

```
docker pull traitecoevo/competition_kernels
```

Then launch it via:

```
docker run --user root -v $(pwd):/home/rstudio/ -p 8787:8787 -e DISABLE_AUTH=true traitecoevo/competition_kernels
```

The code above initialises a docker container, which runs an rstudio session, which is accessed by pointing your browser to [localhost:8787](http://localhost:8787). For more instructions on running docker, see the info from [rocker](https://hub.docker.com/r/rocker/rstudio).

Note, this container does not contain the actual github repo, only the software environment. If you run the above command from within your downloaded repo, it will map the working directory as the current working directory inside the docker container.

### Building the docker images (optional)

For posterity, the docker image was built off [`rocker/verse:3.6.1` container](https://hub.docker.com/r/rocker/verse) via the following command, in a terminal contained within the downloaded repo:

```
docker build -t traitecoevo/competition_kernels .
```

and was then pushed to dockerhub ([here](https://cloud.docker.com/u/traitecoevo/repository/docker/traitecoevo/competition_kernels)). The image used by binder builds off this container, adding extra features needed bi binder, as described in [rocker/binder](https://hub.docker.com/r/rocker/binder/dockerfile).


Contributors
------------------------
Daniel Falster
Rich FitzJohn
Georges Kunstler
