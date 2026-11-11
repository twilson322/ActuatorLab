# Reproducible build + test environment (Ubuntu 24.04, same toolchain as development).
#   docker build -t actuator-lab .
#   docker run --rm actuator-lab                       # runs the unit tests
#   docker run --rm -v $PWD/results:/src/results actuator-lab ./scripts/run_all.sh
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake ninja-build libeigen3-dev libgtest-dev \
      python3 python3-numpy python3-scipy python3-pandas python3-matplotlib \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release && cmake --build build
CMD ["bash", "-c", "cd build && ctest --output-on-failure"]
