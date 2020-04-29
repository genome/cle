FROM rocker/r-ver:3.6.3
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL \
  description="Haplotect Report image"

RUN Rscript -e 'install.packages("ggplot2", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("dplyr", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("binom", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("gridExtra", repos="http://cran.us.r-project.org", dependencies=TRUE)'

COPY HaplotectReport.R /usr/local/bin/HaplotectReport.R  
