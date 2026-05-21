FROM python:3.13-slim

WORKDIR /app

COPY pyproject.toml LICENSE README.md /app/
COPY pyIsoPEP /app/pyIsoPEP

RUN pip install --upgrade pip hatchling build \
 && pip install --no-cache-dir .

COPY example/ /example/

ENTRYPOINT ["pyisopep"]
CMD ["--help"]
