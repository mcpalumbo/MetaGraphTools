# Use Python 3.11 base image
FROM python:3.11-slim

# Set working directory
WORKDIR /app

# Copy the MetaGraphTools folder
COPY . /app

# Install the package
RUN pip install --no-cache-dir -e /app

# Create data directory
RUN mkdir -p /data

# Set environment variable to avoid Python buffering
ENV PYTHONUNBUFFERED=1

# Set the default working directory to /data for user convenience
WORKDIR /data

# Default command shows help
CMD ["MetaGraphTools", "--help"]
