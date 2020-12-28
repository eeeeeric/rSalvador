# This Dockerfile originally used the ubuntu:20.04 image, but installing R
# ourselves proved problematic, as building on Digital Ocean did not play well
# with installing the GPG key (gpg: keyserver receive failed: Address already in use)
FROM rstudio/r-base:4.0-focal

# R is already installed, but we need the dev package for Python to communicate
# with R. We also install a number of SSL related packages that we'll need later
# when adding the GPG key for NGINX Unit (or any other future packages)
RUN apt-get update --fix-missing
RUN apt-get install -y \
  dirmngr \
  gnupg \
  apt-transport-https \
  ca-certificates \
  software-properties-common
RUN apt-get install -y r-base-dev

# Install rSalvador dependencies
RUN R --vanilla -e "install.packages('hypergeo', repos='http://cran.us.r-project.org')"
RUN R --vanilla -e "install.packages('gdata', repos='http://cran.us.r-project.org')"

# Install rSalvador
RUN mkdir rSalvador
WORKDIR rSalvador
ENV RSALVADOR_VERSION=1.8
RUN apt-get install -y wget
RUN wget https://github.com/eeeeeric/rSalvador/releases/download/v1.8/rsalvador_1.8_R_x86_64-pc-linux-gnu.tar.gz
RUN R CMD INSTALL rsalvador_1.8_R_x86_64-pc-linux-gnu.tar.gz

# Setup NGINX Unit
RUN apt-get install -y curl
RUN curl -sL https://nginx.org/keys/nginx_signing.key | apt-key add -
RUN echo "deb https://packages.nginx.org/unit/ubuntu/ focal unit \
deb-src https://packages.nginx.org/unit/ubuntu/ focal unit" > /etc/apt/sources.list.d/unit.list
RUN apt-get update
RUN apt-get install -y unit unit-python3.8

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
