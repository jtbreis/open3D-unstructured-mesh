FROM --platform=linux/amd64 julian-nvim-ubuntu

RUN apt-get update \
    && apt-get install -y build-essential git cmake

RUN git clone https://github.com/isl-org/Open3D
WORKDIR Open3D
RUN sed -i 's/SUDO=\${SUDO:=sudo}/SUDO=command/' util/install_deps_ubuntu.sh \
    && sed -i 's/apt-get install/apt-get install -y/g' util/install_deps_ubuntu.sh \
    && util/install_deps_ubuntu.sh [ assume-yes ]

RUN mkdir build && cd build && cmake -DBUILD_PYTHON_MODULE=OFF -DBUILD_GUI=OFF -DBUILD_EXAMPLES=OFF .. && make -j 4
RUN cd build && make install

WORKDIR /code/src

