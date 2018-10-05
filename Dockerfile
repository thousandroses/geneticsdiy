# Dockerfile for genomics app

FROM ubuntu


# Get g++ for compiling, wget to download Boost, git to clone source code repo,
# and make to automate program compilation with Makefile provided
RUN apt-get update \
  && apt-get install -y git \
                        g++ \
                        make \
                        wget

# Download boost, untar, setup install with bootstrap and only do the Program Options library,
# and then install
RUN cd /home && wget http://downloads.sourceforge.net/project/boost/boost/1.60.0/boost_1_60_0.tar.gz \
  && tar xfz boost_1_60_0.tar.gz \
  && rm boost_1_60_0.tar.gz \
  && cd boost_1_60_0 \
  && ./bootstrap.sh --prefix=/usr/local --with-libraries=program_options \
  && ./b2 install \
  && cd /home \
  && rm -rf boost_1_60_0
# Clone git repo for the ancestry program we need
RUN cd /home \
  && git clone https://github.com/pickleriiick/ancestry_mirror
RUN apt-get update && \
    apt-get install -y \
        zlib1g-dev
RUN apt-get update && \
    apt-get install -y \
        libgsl0-dev
RUN apt-get update && \
    apt-get install -y \
        libbz2-dev
RUN DEBIAN_FRONTEND=noninteractive apt-get install -y python3-tk
RUN apt-get install -y liblzma-dev
RUN apt-get install -y python-pip python-virtualenv python3-pip
RUN pip install --upgrade pip 
RUN pip install --upgrade virtualenv 
RUN cd /usr/bin && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2
RUN cd /usr/bin && tar -vxjf htslib-1.9.tar.bz2
RUN cd /usr/bin/htslib-1.9/ && make
RUN cd /usr/bin/htslib-1.9/ && make install
ENV LD_LIBRARY_PATH /usr/bin/htslib-1.9/
# Configure and make the ancestry program
RUN cd /home/ancestry/ && ./configure && make clean && make
# Add the folder files to the container's root directory
ADD . .
RUN chmod +x /home/ancestry/src/ancestry.cpp

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        libatlas-base-dev gfortran
RUN apt-get -y install graphviz libgraphviz-dev pkg-config
RUN apt-get -y install python-tk
RUN apt-get -y install libgeos-dev
# Install basemap
RUN pip install --user https://github.com/matplotlib/basemap/archive/master.zip
RUN pip install --upgrade pip
# Install the python modules specified in requirements.txt
RUN pip install --trusted-host pypi.python.org -r /requirements.txt

# Run controller.py when the container launches
ENTRYPOINT ["python", "/controller.py"]
