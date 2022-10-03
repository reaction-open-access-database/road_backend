FROM python:3.10.7

# Install Rust
ENV HOME=/home/app
RUN curl https://sh.rustup.rs -sSf | sh -s -- -y
ENV PATH="/home/app/.cargo/bin:$PATH"
RUN pip install --upgrade pip
RUN cargo install cargo-chef

# Set up the working directory
ENV DOCKER_HOME=/home/app/road
RUN mkdir -p $DOCKER_HOME
WORKDIR $DOCKER_HOME

# Python environment variables
ENV PYTHONDONTWRITEBYTECODE 1
ENV PYTHONUNBUFFERED 1

# Install python dependencies
COPY ./requirements.txt ./
RUN pip install -r requirements.txt

# Build rust dependencies
RUN mkdir -p ./road/query_parser
COPY ./road/query_parser/Cargo.toml ./road/query_parser/
WORKDIR $DOCKER_HOME/road/query_parser
RUN cargo chef prepare
RUN cargo chef cook --release
WORKDIR $DOCKER_HOME

# Build rust project
COPY ./road/query_parser ./road/query_parser
WORKDIR $DOCKER_HOME/road/query_parser
RUN maturin build -r
RUN pip install target/wheels/*

# Entrypoint
WORKDIR $DOCKER_HOME
COPY ./docker-entrypoint.sh ./
RUN chmod +x ./docker-entrypoint.sh
ENTRYPOINT ["./docker-entrypoint.sh"]

# Copy the rest of the application
COPY . .