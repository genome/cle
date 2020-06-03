FROM rocker/r-ver:3.6.3
MAINTAINER Feiyu Du <fdu@wustl.edu>

LABEL \
  description="Somatic Concordance Report Graphic Docker Image"

RUN Rscript -e 'install.packages("ggplot2", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("dplyr", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("reshape2", repos="http://cran.us.r-project.org", dependencies=TRUE)'
RUN Rscript -e 'install.packages("gridExtra", repos="http://cran.us.r-project.org", dependencies=TRUE)'

COPY SomaticConcordanceReport.R /usr/local/bin/SomaticConcordanceReport.R  
