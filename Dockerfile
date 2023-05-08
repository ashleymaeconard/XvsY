FROM ubuntu:latest

ARG DEBIAN_FRONTEND=noninteractive

# RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-cran-randomforest python3.10 python3-pip python3-setuptools python3-dev
RUN apt-get update && apt-get install -y --no-install-recommends \
 build-essential \
 r-base \
 r-base-dev \
 python3.10 \
 python3-pip \
 python3-setuptools \
 python3-dev

RUN apt-get install -y --no-install-recommends \
 libxml2-dev \
 libpq-dev \
 libssl-dev \
 libcurl4-openssl-dev \
 libfontconfig1-dev \
 libharfbuzz-dev \
 libfribidi-dev \
 libgit2-dev \
 libfreetype6-dev \
 libpng-dev \
 libtiff5-dev \
 libjpeg-dev \
 libcairo2-dev \
 libxt-dev \
 libproj-dev

# Copy Necessary Python and R Requirement Files and Install
COPY ./home ./home
RUN Rscript /home/requirements.r
RUN pip3 install -r /home/requirements.txt

# Set up the User Interface
RUN apt-get -y install sudo

RUN useradd -m -s /bin/bash xvsy-user && \
  echo "xvsy-user ALL=(ALL:ALL) NOPASSWD: ALL" > /etc/sudoers.d/xvsy-init

USER xvsy-user
RUN rm -f ~/.bash_logout

WORKDIR /home/analysis

CMD ["/bin/bash", "-l"]