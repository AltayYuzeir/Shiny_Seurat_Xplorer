# get shiny server plus tidyverse packages image
FROM rocker/shiny-verse:4.2.1

# install system libraries of general use
# Add components if needed
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libpoppler-cpp-dev \
	libglpk-dev \
	libmagic-dev \
	libmagick++-dev\
	libgeos-dev\
    && apt-get clean && rm -rf /var/lib/apt/lists/*


# Install required R packages from install.R file
# Modify install.R for your needs
COPY --chown=shiny:shiny install.R /tmp/
RUN Rscript /tmp/install.R && \
    rm -rf /tmp/*


############################################################
### Do not modify this section for usage on FASTGenomics ###
############################################################
# Copy FASTGenomics config 
COPY fg_config/shiny-server.conf /etc/shiny-server/
COPY fg_config/run-shiny.sh /usr/bin
# Allow permissions for FASTGenomics config
RUN chown -R shiny:shiny /etc/shiny-server/
RUN chown -R shiny:shiny /usr/bin/run-shiny.sh
RUN chmod u+rwx /usr/bin/run-shiny.sh
#
# Copy APP data
COPY APP/ /app/
# Allow permissions for APP data
RUN chown -R shiny:shiny /app/
############################################################


EXPOSE 3838

USER shiny

# run app command
# If you start the container locally you can reach the app under: localhost:3838/proxy/shiny
# Don't forget to map port 3838 when testing locally
CMD ["/usr/bin/run-shiny.sh", "shiny"]
