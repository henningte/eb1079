FROM rocker/verse:4.3.0

# system dependencies
RUN export DEBIAN_FRONTEND=noninteractive; apt-get -y update \
  && apt-get install -y gdal-bin=3.4.1+dfsg-1build4 \
	libgdal-dev=3.4.1+dfsg-1build4 \
	libgeos-dev=3.10.2-1 \
	libgeos++-dev=3.10.2-1 \
	libudunits2-dev=2.2.28-3 \
	libnng-dev=1.5.2-1build1 \
	liblzma5=5.2.5-2ubuntu1 \
	xz-utils=5.2.5-2ubuntu1 \
	make=4.3-4.1build1 \
	pandoc=2.9.2.1-3ubuntu2 \
	pandoc-citeproc=0.17.0.1-1build5

# working directory name
WORKDIR /home/rstudio/eb1079

# install renv R package
ENV RENV_VERSION 0.16.0
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cran.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

# Restore renv cache
RUN mkdir -p renv
COPY renv.lock renv.lock
COPY .Rprofile .Rprofile
COPY renv/activate.R renv/activate.R
RUN chown -R rstudio . \
 && sudo -u rstudio Rscript -e 'renv::consent(provided = TRUE)' \
 && sudo -u rstudio Rscript -e 'Sys.setenv(R_INSTALL_STAGED = FALSE)' \
 && sudo -u rstudio R -e 'renv::restore(repos = c(RSPM = "https://packagemanager.rstudio.com/all/latest", CRAN = "https://cran.r-project.org/"))'

# labels
LABEL maintainer = "Henning Teickner <henning.teickner@uni-muenster.de>" \
  org.opencontainers.image.authors = "Henning Teickner <henning.teickner@uni-muenster.de>" \
  author.orcid = "0000-0002-3993-1182" \
  org.opencontainers.image.version = "0.0.0.9000" \
  org.opencontainers.image.licenses = "GPL-3"

# instructions

# to build the image, navigate to the directory with the Dockerfile and run:
# docker build -t eb1079:0.0.0.9000 .

# to run the image in a container, do:
# docker run --name eb1079_c -e PASSWORD=eb1079 --net r-db -p 8794:8787 -v $(pwd):/home/rstudio/eb1079 eb1079:0.0.0.9000

# to reproduce the analyses from all projects, run:
# docker run --rm -v $(pwd):/home/rstudio eb1079:0.0.0.9000 R -e 'setwd("/home/rstudio/eb1079");
#
# To reproduce the analysis, you will have to set up the pmird database in an own MariaDB server that is connected to the eb1079 container via the `r-db` network. Instructions are available from https://doi.org/10.5281/zenodo.17092587. The mid-infrared spectra from the pmird database also need to be downloaded from https://doi.org/10.5281/zenodo.17092587 and the folder `pmird_prepared_data` stored in the `data/raw_data/pmird` folder available from the repository where this Dockerfile here is also available from.
# After this setup you can execute 'run.R' in the folders in the targets folder (takes about 20 h and occupies additional 50 Gb disk space).

