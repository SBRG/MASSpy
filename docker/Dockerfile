FROM python:3.7-stretch

WORKDIR /usr/src/app

# Copy MASSpy contents
COPY ./ ./

RUN echo \
    # Install MASSpy
    && pip install .[all] 