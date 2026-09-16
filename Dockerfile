FROM rust:1.96-slim-bookworm AS builder
WORKDIR /app
COPY Cargo.toml Cargo.lock ./
COPY src ./src
RUN cargo build --release --locked

FROM debian:bookworm-slim
RUN apt-get update && apt-get install -y --no-install-recommends ca-certificates && rm -rf /var/lib/apt/lists/*
COPY --from=builder /app/target/release/protein-tools /usr/local/bin/protein-tools
ENV CACHE_DB=/data/cache.sqlite CACHE_TTL_DAYS=30
EXPOSE 8080
CMD ["protein-tools"]
