pkg_name=aptamer-scripts
pkg_origin=chrisortman
pkg_version="2.0.0"
pkg_deps=(
  core/groff
  core/python2
  core/glibc
  core/gcc-libs
  core/gcc
  core/zlib
  chrisortman/ViennaRNA
  chrisortman/mfold
)
pkg_build_deps=(
  core/cacerts
  core/virtualenv
  core/git
  core/tar
  core/coreutils
)
pkg_bin_dirs=(bin)

do_begin() {
  SRC_PATH=$CACHE_PATH
  return 0
}

# The default implementation is that the software specified in $pkg_source is
# downloaded, checksum-verified, and placed in $HAB_CACHE_SRC_PATH/$pkgfilename,
# which resolves to a path like /hab/cache/src/filename.tar.gz. You should
# override this behavior if you need to change how your binary source is
# downloaded, if you are not downloading any source code at all, or if your are
# cloning from git. If you do clone a repo from git, you must override
# do_verify() to return 0.
do_download() {
  return 0
}

do_verify() {
  return 0
}

do_unpack() {
  build_line "Cloning source files"
  cd $PLAN_CONTEXT/../
  { git ls-files; git ls-files --exclude-standard --others; } \
    | _tar_pipe_app_cp_to $HAB_CACHE_SRC_PATH/${pkg_dirname}

}

do_prepare() {
 virtualenv "$pkg_prefix"
 # shellcheck source=/dev/null
 source "$pkg_prefix/bin/activate"
}

do_build() {
  return 0
}

do_check() {

  bin/python -c "import numpy; import scipy.stats; from Bio import SeqIO; import Levenshtein; import scipy.special"
}

do_install() {
  pip install --upgrade pip
  pip install numpy
  pip install scipy
  pip install biopython
  pip install python-Levenshtein
  
  virtualenv --relocatable $pkg_prefix
  
  mkdir -p ${pkg_prefix}/scripts
  cp -a aptamer_functions.py ${pkg_prefix}/scripts
  cp -a create_graph.py ${pkg_prefix}/scripts
  cp -a find_families.py ${pkg_prefix}/scripts
  cp -a predict_structures.py ${pkg_prefix}/scripts

  cat << DO_SCRIPT > ${pkg_prefix}/bin/python-wrapper
#!/bin/sh

# This script is used to invoke python
# but ensure that gcc libraries are able to be
# found (that's why we set LD_LIBRARY_PATH)

# It also ensures that python output is not buffered.

export LD_LIBRARY_PATH=$(pkg_path_for core/gcc-libs)/lib:$(pkg_path_for core/gcc)/lib
export PATH=$(pkg_path_for chrisortman/ViennaRNA)/bin:$(pkg_path_for chrisortman/mfold)/bin\$PATH
${pkg_prefix}/bin/python -u "\$@"
DO_SCRIPT

  chmod +x ${pkg_prefix}/bin/python-wrapper

  cat << DO_SCRIPT > ${pkg_prefix}/bin/run-aptamer
#!/bin/sh

# This script ensures that we are invoked within
# our environment

hab pkg exec ${pkg_origin}/${pkg_name} python-wrapper "\$@"
DO_SCRIPT

  chmod +x ${pkg_prefix}/bin/run-aptamer
}

do_strip() {
  return 0
}

do_end() {
  return 0
}

_tar_pipe_app_cp_to() {
  local dst_path tar
  dst_path="$1"
  tar="$(pkg_path_for tar)/bin/tar"

  mkdir -p $dst_path

  "$tar" -cp \
      --owner=root:0 \
      --group=root:0 \
      --no-xattrs \
      --exclude-backups \
      --exclude-vcs \
      --exclude='habitat' \
      --exclude='node_modules' \
      --exclude='vendor/bundle' \
      --exclude='results' \
      --files-from=- \
      -f - \
  | "$tar" -x \
      -C "$dst_path" \
      -f -
}

