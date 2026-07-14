# =============================================================================
# SimVascular full-GUI build environment
# Base: Ubuntu 22.04 (GCC 11) -- chosen so that the 2024.05 external set
# (Qt 6.6.2 / VTK 9.3.0 / ITK 5.4.0 / OpenCASCADE 7.6.0 / MITK Qt6) compiles
# with a compiler generation it was released against.
#
# This image contains ONLY the toolchain and system libraries.
# Externals and SimVascular itself are built inside the container by the
# scripts in Docker/scripts/, with all big outputs going to the mounted
# volume (/work), so nothing is lost when the container is removed.
# =============================================================================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# ---- Core toolchain --------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran cmake ninja-build \
      git wget curl ca-certificates unzip zip patch pkg-config file \
      python3 python3-pip python3-venv swig \
    && rm -rf /var/lib/apt/lists/*

# ---- Library dependencies --------------------------------------------------
# Sources of truth for this list:
#   - Externals/Make/2022.10/compilation-notes-linux.txt (in this repo)
#   - .github/workflows/build.yml (ubuntu job)
#   - SimVascular issue #1065 (runtime libs)
RUN apt-get update && apt-get install -y --no-install-recommends \
      # CPython build (ssl/sqlite/lzma/readline/ffi/uuid/zlib)
      libssl-dev libsqlite3-dev libbz2-dev liblzma-dev libreadline-dev \
      libffi-dev uuid-dev zlib1g-dev \
      # OpenGL / VTK
      libglu1-mesa-dev freeglut3-dev mesa-common-dev libglvnd-dev \
      libglew-dev libxt-dev \
      # OpenCASCADE
      libxmu-dev libxi-dev libfontconfig1-dev \
      # MITK
      libtiff-dev libboost-all-dev doxygen \
      # Qt6 xcb platform + runtime
      libx11-xcb-dev libxrender-dev libxkbcommon-dev libxkbcommon-x11-dev \
      libxcb1 libxcb-cursor0 libxcb-glx0 libxcb-icccm4 libxcb-image0 \
      libxcb-keysyms1 libxcb-randr0 libxcb-render-util0 libxcb-render0 \
      libxcb-shape0 libxcb-shm0 libxcb-sync1 libxcb-util1 libxcb-xfixes0 \
      libxcb-xinerama0 libxcb-xkb1 \
      # Qt WebEngine / GUI runtime
      libnss3 libasound2 libxss1 libxcursor1 libxcomposite1 libxdamage1 \
      libxrandr2 libxtst6 libgl1 libegl1 libdbus-1-3 libfontconfig1 \
      libxkbfile1 libgbm1 \
      # X11 dev headers commonly needed by Qt/VTK builds
      libx11-dev libxext-dev libxfixes-dev libxinerama-dev libxcursor-dev \
    && rm -rf /var/lib/apt/lists/*

# aqtinstall fetches official Qt binary releases from Qt mirrors.
RUN pip3 install --no-cache-dir aqtinstall

# Non-root build user so files created on the mounted volume are not root-owned.
ARG USERNAME=sv
ARG UID=1000
ARG GID=1000
RUN groupadd -g ${GID} ${USERNAME} && \
    useradd -m -u ${UID} -g ${GID} -s /bin/bash ${USERNAME}
USER ${USERNAME}
WORKDIR /work

CMD ["/bin/bash"]
