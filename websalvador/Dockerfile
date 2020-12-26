FROM ubuntu:20.04

# Install R
RUN apt-get update --fix-missing
RUN apt-get install -y dirmngr gnupg apt-transport-https ca-certificates software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'
RUN apt-get install -y r-base r-base-dev
RUN apt-get install -y build-essential

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

# Copy app
COPY app app
WORKDIR app

# Setup Python
RUN apt-get install -y python3-pip
RUN pip3 install -r requirements.txt

# Setup NGINX Unit
RUN apt-get install curl
RUN curl -sL https://nginx.org/keys/nginx_signing.key | apt-key add -
COPY unit.list /etc/apt/sources.list.d/unit.list
RUN apt-get update
RUN apt-get install unit
RUN apt-get install unit-dev unit-python3.8

# Config
COPY config.json .
EXPOSE 80

RUN apt-get install -y less
COPY entrypoint.sh .
ENTRYPOINT ["./entrypoint.sh"]