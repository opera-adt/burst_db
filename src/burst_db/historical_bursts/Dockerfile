# This is the image we are going add micromaba to:
ARG BASE=ubuntu:22.04

# bring in the micromamba image so we can copy files from it
FROM mambaorg/micromamba:1.5.0 as micromamba

FROM $BASE

# Needed to redefine the arg so we can use it later
ARG BASE

# Install CA certificates if the base image has apt-get, otherwise we can't install conda packages
RUN if command -v apt-get > /dev/null; then \
    apt-get update && \
    apt-get install -y ca-certificates && \
    rm -rf /var/lib/apt/lists/*; \
    fi

ARG MAMBA_USER=mamba
ARG MAMBA_USER_ID=1000
ARG MAMBA_USER_GID=1000
ENV MAMBA_USER=$MAMBA_USER
ENV MAMBA_ROOT_PREFIX="/opt/conda"
ENV MAMBA_EXE="/bin/micromamba"

COPY --from=micromamba "$MAMBA_EXE" "$MAMBA_EXE"
COPY --from=micromamba /usr/local/bin/_activate_current_env.sh /usr/local/bin/_activate_current_env.sh
COPY --from=micromamba /usr/local/bin/_dockerfile_shell.sh /usr/local/bin/_dockerfile_shell.sh
COPY --from=micromamba /usr/local/bin/_entrypoint.sh /usr/local/bin/_entrypoint.sh
# The `RUN true` is a weird hack to prevent "layer does not exist" errors
# apparently arising from multi-stage builds? https://stackoverflow.com/a/62409523/4174466
RUN true
COPY --from=micromamba /usr/local/bin/_activate_current_env.sh /usr/local/bin/_activate_current_env.sh
RUN true
COPY --from=micromamba /usr/local/bin/_dockerfile_initialize_user_accounts.sh /usr/local/bin/_dockerfile_initialize_user_accounts.sh
RUN true
COPY --from=micromamba /usr/local/bin/_dockerfile_setup_root_prefix.sh /usr/local/bin/_dockerfile_setup_root_prefix.sh

RUN /usr/local/bin/_dockerfile_initialize_user_accounts.sh && \
    /usr/local/bin/_dockerfile_setup_root_prefix.sh


# Location to install source
WORKDIR /src
# Allow us to create new files in the workdir when pip installing
RUN chown $MAMBA_USER:$MAMBA_USER /src

USER $MAMBA_USER

SHELL ["/usr/local/bin/_dockerfile_shell.sh"]


# Install conda packages
# https://github.com/mamba-org/micromamba-docker#quick-start
COPY --chown=$MAMBA_USER:$MAMBA_USER specfile.txt /tmp/specfile.txt
RUN micromamba install --yes --channel conda-forge -n base -f /tmp/specfile.txt && \
    micromamba clean --all --yes

COPY --chown=$MAMBA_USER:$MAMBA_USER . .

# Activate, otherwise python will not be found
# https://github.com/mamba-org/micromamba-docker#running-commands-in-dockerfile-within-the-conda-environment
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN python -m pip install https://github.com/scottstanie/asfsmd/archive/s3fs-client.zip

# Run the entrypoint from the /work directory
# This means if they mount a volume to /work, we won't be mixing in the repository code
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
# CMD ["download.py"]
