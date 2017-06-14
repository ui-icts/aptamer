# This file is the heart of your application's habitat.
# See full docs at https://www.habitat.sh/docs/reference/plan-syntax/

# Required.
# Sets the name of the package. This will be used in along with `pkg_origin`,
# and `pkg_version` to define the fully-qualified package name, which determines
# where the package is installed to on disk, how it is referred to in package
# metadata, and so on.
pkg_name=aptamer-scripts

# Required unless overridden by the `HAB_ORIGIN` environment variable.
# The origin is used to denote a particular upstream of a package.
pkg_origin=chrisortman

# Required.
# Sets the version of the package.
pkg_version="0.1.0"

# Optional.
# The name and email address of the package maintainer.
# pkg_maintainer="The Habitat Maintainers <humans@habitat.sh>"

# Optional.
# An array of valid software licenses that relate to this package.
# Please choose a license from http://spdx.org/licenses/
# pkg_license=('Apache-2.0')

# Required.
# A URL that specifies where to download the source from. Any valid wget url
# will work. Typically, the relative path for the URL is partially constructed
# from the pkg_name and pkg_version values; however, this convention is not
# required.
pkg_source="https://github.com/ui-icts/aptamer/archive/${pkg-name}-${pkg-version}.tar.bz2"

pkg_shasum="TODO"
pkg_deps=(
  core/groff
  core/python2
  core/gcc
)
pkg_build_deps=(
  core/cacerts
  core/virtualenv
  core/gcc
  core/git
)
pkg_bin_dirs=(bin)

do_begin() {
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
 export GIT_SSL_CAINFO="$(pkg_path_for core/cacerts)/ssl/certs/cacert.pem"

      # This is a way of getting the git code that I found in the chef plan
 build_line "Fake download! Creating archive of latest repository commit from $PLAN_CONTEXT"
          cd $PLAN_CONTEXT/..
  git archive --prefix=${pkg_name}-${pkg_version}/ --output=$HAB_CACHE_SRC_PATH/${pkg_filename} HEAD

  pkg_shasum=$(trim $(sha256sum $HAB_CACHE_SRC_PATH/${pkg_filename} | cut -d " " -f 1))

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

  bin/python -c "import numpy; import scipy.stats; from Bio import SeqIO; import Levenshtein;"
}

do_install() {
  pip install --upgrade pip
  pip install numpy
  pip install scipy
  pip install biopython
  pip install python-Levenshtein

  mkdir -p ${pkg_prefix}/scripts
  cp -a aptamer_functions.py ${pkg_prefix}/scripts
  cp -a create_graph.py ${pkg_prefix}/scripts
  cp -a find_families.py ${pkg_prefix}/scripts
  cp -a predict_structures.py ${pkg_prefix}/scripts

}

do_strip() {
  return 0
}

do_end() {
  return 0
}

