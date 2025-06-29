# Build image using the following command:
#
#   docker build -t birch:latest .
#
# Run the container from the built image using the following command:
#
#   docker run --rm -it birch:latest
#
# Once in the container, use the following to import BirchGenus:
#
#   from ternary_birch import BirchGenus
#

FROM sagemath/sagemath:latest
COPY ./src /var/src/birch
WORKDIR /var/src/birch
RUN sudo sage setup.py install
WORKDIR /home/sage
USER root
