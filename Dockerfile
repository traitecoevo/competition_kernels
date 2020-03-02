FROM rocker/verse:3.6.1
LABEL maintainer="Daniel Falster"
LABEL email="daniel.falster@unsw.edu.au"

# Install major libraries
RUN    apt-get update \
    && apt-get install -y --no-install-recommends \
        zip \
        unzip

# ---------------------------------------------

ENV NB_USER rstudio
ENV NB_UID 1000

# And set ENV for R! It doesn't read from the environment...
RUN echo "PATH=${PATH}" >> /usr/local/lib/R/etc/Renviron
RUN echo "export PATH=${PATH}" >> ${HOME}/.profile

# The `rsession` binary that is called by nbrsessionproxy to start R doesn't seem to start
# without this being explicitly set
ENV LD_LIBRARY_PATH /usr/local/lib/R/lib

ENV HOME /home/${NB_USER}
WORKDIR ${HOME}

# ---------------------------------------------

# Install extra latex style files 

## Saves downloading whenever tinytex runs
## Inspired by Yihui - https://github.com/yihui/tinytex/issues/135#issuecomment-514351695
## NB: tinytex extras installed at  /opt/TinyTeX/tlpkg/TeXLive/

RUN R --quiet -e 'tinytex::tlmgr_install(c("adjustbox", "a4wide", "algorithms", "appendix", "babel-english", "bbm-macros", "beamer", "breakurl", "catoptions", "charter", "cite", "cleveref", "collectbox", "colortbl", "comment", "courier", "eepic", "enumitem", "eso-pic", "eurosym", "extsizes", "fancyhdr", "floatrow", "fontaxes", "fpl", "hardwrap", "koma-script", "lastpage", "lettrine", "libertine", "lineno", "lipsum",  "ltxkeys", "ly1", "mathalpha", "mathpazo", "mathtools", "mdframed", "mdwtools", "microtype", "morefloats", "ms", "multirow", "mweights", "ncctools", "ncntrsbk", "needspace", "newtx", "ntgclass", "numname", "palatino", "pbox", "pdfpages", "pgf", "picinpar", "preprint", "preview", "psnfss", "refstyle", "roboto", "sectsty", "setspace", "siunitx", "srcltx", "standalone", "stmaryrd", "sttools", "subfig", "subfigure", "symbol", "tabu", "textcase", "threeparttable", "thumbpdf", "titlesec", "tufte-latex", "ucs", "ulem", "units", "varwidth", "vmargin", "wallpaper", "wrapfig", "xargs", "xcolor", "xstring", "xwatermark"))';

# ---------------------------------------------

# Add custom installations here

# install specific version of rootsolve
# install plant
# install other packages for analysis
# install remake

RUN R --quiet -e 'options(repos = list(CRAN = "http://mran.revolutionanalytics.com/snapshot/2019-08-26/"));         \
remotes::install_version("rootSolve", version="1.6.6");               \
install.packages(c("Rcpp", "R6", "crayon", "nleqslv", "BB" ,"BH", "numDeriv", "deSolve", "yaml", "digest", "crayon", "optparse"));   \
remotes::install_github(c("smbache/loggr", "traitecoevo/plant@v1.1.0", "richfitz/storr", "richfitz/remake"));'
