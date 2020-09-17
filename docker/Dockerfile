FROM python:3.7
ARG version=latest

WORKDIR /usr/src/app

# TODO add labels
# LABEL

RUN echo \
    # Navigate to source directory
    && cd .. \
    # Clone MASSpy GitHub repository and navigate to it
    && git clone https://github.com/SBRG/MASSpy.git \
    && cd MASSpy \
    # Checkout a specific version of MASSpy if desired
    && if [ "${version}" != "latest" ] ; then git checkout "v${version}";fi \
    && pip install . \
    # Return to source directory
    && cd ..