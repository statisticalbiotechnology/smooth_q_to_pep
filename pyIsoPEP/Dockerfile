FROM python:3.13-slim

WORKDIR /app

# Copy the requirements file and install dependencies
COPY pyIsoPEP/requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt
COPY pyIsoPEP/ ./
COPY example/ /example/

# Set the entrypoint to run the main script.
# This will allow users to pass CLI arguments directly.
ENTRYPOINT ["python", "/app/main.py"]
