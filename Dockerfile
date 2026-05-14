# Development image for SemiSim.
#
# Based on the official dolfinx image, which provides FEniCSx (dolfinx, ufl,
# basix, ffcx), PETSc, MPI, and a compatible Python. We layer the project's
# Python dev dependencies and install the package in editable mode on top.
#
# Build:   docker build -t semisim:dev .
# Run:     docker compose up -d dev   (preferred, see docker-compose.yml)

FROM ghcr.io/fenics/dolfinx/dolfinx:stable

ARG USERNAME=dev
ARG USER_UID=1000
ARG USER_GID=1000

ENV DEBIAN_FRONTEND=noninteractive \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

# System tools useful inside a dev container. The base image already has
# build-essential, git, and Python; we add interactive conveniences.
RUN apt-get update \
 && apt-get install -y --no-install-recommends \
        sudo \
        vim \
        less \
        curl \
        ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# Create a non-root user that matches the host UID/GID so bind-mounted files
# stay owned by the host user. If a user or group already exists at the
# requested UID/GID (the dolfinx base image ships an `ubuntu` user at 1000),
# rename it to ${USERNAME} rather than failing.
RUN set -eux; \
    if getent group ${USER_GID} >/dev/null; then \
        existing_group="$(getent group ${USER_GID} | cut -d: -f1)"; \
        if [ "${existing_group}" != "${USERNAME}" ]; then \
            groupmod -n ${USERNAME} "${existing_group}"; \
        fi; \
    else \
        groupadd --gid ${USER_GID} ${USERNAME}; \
    fi; \
    if getent passwd ${USER_UID} >/dev/null; then \
        existing_user="$(getent passwd ${USER_UID} | cut -d: -f1)"; \
        if [ "${existing_user}" != "${USERNAME}" ]; then \
            usermod -l ${USERNAME} -d /home/${USERNAME} -m "${existing_user}"; \
            usermod -g ${USER_GID} ${USERNAME}; \
        fi; \
    else \
        useradd --uid ${USER_UID} --gid ${USER_GID} -m -s /bin/bash ${USERNAME}; \
    fi; \
    echo "${USERNAME} ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/${USERNAME}; \
    chmod 0440 /etc/sudoers.d/${USERNAME}

WORKDIR /workspaces/SemiSim

# Two-stage install so the heavy dependency layer caches across source-only
# edits (Dockerfile change #1 in CI-speedup pass):
#
#   1. Copy *only* pyproject.toml + README.md plus minimal package stubs and
#      run `pip install -e ".[dev,server,test]"`. This layer's cache key
#      depends only on pyproject.toml, so a typical PR (which only edits
#      semi/ or tests/) hits the cache and skips re-resolving wheels.
#   2. `COPY . .` overlays the real source. The hatchling editable install
#      from step 1 emits a redirector for packages=["semi", "semisim_server"]
#      that resolves through PYTHONPATH/finders to whatever lives at those
#      paths at runtime, so we do *not* need a second `pip install -e`.
#
# We add jupyterlab + ipykernel here because they aren't pyproject deps but
# are wanted in the dev image. BuildKit's pip cache mount keeps wheels
# across builds without bloating the final image (PIP_NO_CACHE_DIR is
# unchanged for runtime; the mount supplies its own cache directory).
COPY pyproject.toml README.md ./
COPY schemas ./schemas
RUN mkdir -p semi semisim_server \
 && touch semi/__init__.py semisim_server/__init__.py
RUN --mount=type=cache,target=/root/.cache/pip,sharing=locked \
    PIP_NO_CACHE_DIR=0 pip install --break-system-packages \
        -e ".[dev,server,test]" \
        "jupyterlab>=4.0" \
        "ipykernel>=6.0"

# Overlay the real source. The bind mount in docker-compose / CI will
# overlay /workspaces/SemiSim again at runtime; this COPY guarantees
# `import semi` works when the image is run without a bind mount.
COPY . .

RUN chown -R ${USER_UID}:${USER_GID} /workspaces

USER ${USERNAME}

CMD ["bash"]
