FROM python:3.13-slim

WORKDIR /app

# Copy the requirements file and install dependencies
COPY pyproject.toml LICENSE /app/
COPY pyIsoPEP /app/pyIsoPEP
RUN pip install --no-cache-dir .
COPY example/ /example/

# Set the entrypoint to run the main script.
ENTRYPOINT ["pyisopep"]
CMD ["--help"]
