FROM docker/dev-environments-default:stable-1
WORKDIR /
COPY . /ear-oracle
RUN requiredpackage='build-essential cmake libcgal-dev tmux' \
&&  apt-get update \
&&  apt-get -y upgrade \
&&  apt-get install -y $requiredpackage
ENTRYPOINT sleep infinity
