# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/make_prg:c61d9b

# This container has make_prg installed, based on Michael Hall's Singularity recipe.
# to build: sudo docker build . -t leandroishilima/make_prg:c61d9b

FROM ubuntu:focal

#ARGS
ARG URL="https://github.com/leoisl/make_prg.git"
ARG TAG="c61d9b5af33a580c90c60da250280460aa302fbc"


# INSTALL DEPS
RUN apt update && apt -y upgrade && apt install -y git build-essential python3 python3-pip python-is-python3


# INSTALL MAKE_PRG
RUN git clone "$URL" && \
    cd /make_prg || exit 1 && \
    git checkout "$TAG" && \
    python -m pip install --no-cache-dir . && \
    make_prg --help


# TESTS
RUN cd /make_prg && \
    python -m pip install --no-cache-dir nose hypothesis pandas && \
    nosetests tests/

