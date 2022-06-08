FROM rstudio/r-base:4.2.0-focal

# R is already installed, but we need the dev package for Python to communicate
# with R. We also install a number of SSL related packages that we'll need later
# when adding the GPG key for NGINX Unit (or any other future packages)
RUN apt-get update --fix-missing
RUN apt-get install -y \
  dirmngr \
  gnupg \
  apt-transport-https \
  ca-certificates \
  software-properties-common \
  wget \
  curl

# Install R development library. https://cran.r-project.org/bin/linux/ubuntu/
# update indices
RUN apt-get update -qq
# add the signing key (by Michael Rutter) for these repos
# To verify key, run gpg --show-keys /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# Fingerprint: E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
# add the R 4.0 repo from CRAN -- adjust 'focal' to 'groovy' or 'bionic' as needed
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/"
#RUN apt-get update -qq
RUN apt-get install -y r-base-dev

# Setup NGINX Unit
RUN curl --output /usr/share/keyrings/nginx-keyring.gpg https://unit.nginx.org/keys/nginx-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/nginx-keyring.gpg] https://packages.nginx.org/unit/ubuntu/ focal unit \n\
          deb-src [signed-by=/usr/share/keyrings/nginx-keyring.gpg] https://packages.nginx.org/unit/ubuntu/ focal unit" > /etc/apt/sources.list.d/unit.list
RUN apt-get update
RUN apt-get install -y unit unit-python3.8

# Install rSalvador dependencies
RUN R --vanilla -e "install.packages('hypergeo', repos='http://cran.us.r-project.org')"
RUN R --vanilla -e "install.packages('gdata', repos='http://cran.us.r-project.org')"

# Install rSalvador
RUN mkdir rSalvador
WORKDIR rSalvador
ENV RSALVADOR_VERSION=1.9
RUN wget https://github.com/eeeeeric/rSalvador/releases/download/v1.9/rsalvador_1.9_R_x86_64-pc-linux-gnu.tar.gz
RUN R CMD INSTALL rsalvador_1.9_R_x86_64-pc-linux-gnu.tar.gz

# Config for NGINX Unit (the config is applied as part of the entrypoint)
COPY config.json .
EXPOSE 80

# Copy app
COPY app app
WORKDIR app

# Setup Python
RUN apt-get install -y python3-pip
RUN pip3 install -r requirements.txt

COPY entrypoint.sh .
ENTRYPOINT ["./entrypoint.sh"]
CMD ["/bin/bash"]