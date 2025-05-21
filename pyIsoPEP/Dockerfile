FROM python:3.13-slim

RUN apt-get update \
 && apt-get install -y --no-install-recommends \
      build-essential \
      gfortran \
      libatlas-base-dev \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Copy the requirements file and install dependencies
COPY pyproject.toml LICENSE README.md /app/
COPY pyIsoPEP /app/pyIsoPEP

RUN pip install --upgrade pip hatchling build \
 && pip install --no-cache-dir .

COPY example/ /example/

# Set the entrypoint to run the main script.
ENTRYPOINT ["pyisopep"]
CMD ["--help"]
