###
#
#  @file pastix6.rb
#  @copyright 2013-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
#                       Univ. Bordeaux. All rights reserved.
#
#  @brief Homebrew formula for PaStiX 6
#
#  @version 6.0.3
#  @author Pierre Ramet
#  @date 2019-11-12
#
###
class Pastix6 < Formula
  desc "Parallel solver for sparse linear systems based on direct methods"
  homepage "https://gitlab.inria.fr/solverstack/pastix"
  url "https://gitlab.inria.fr/solverstack/pastix//uploads/bd29252eb0390fa7bcb4e4d238c5b30a/pastix-6.0.3.tar.gz"
  sha256 "8d98799f8f90fdf4849340d9f3214d4b041bca0cdb2e9d8b410587707d569f8c"
  head "https://gitlab.inria.fr/solverstack/pastix.git"

  bottle :disable, "needs to be rebuilt"

  depends_on "openblas"
  depends_on "hwloc"              # Could be optinal but strongly recommanded
  depends_on "scotch5"            # Could be optinal but strongly recommanded
  depends_on "metis" => :optional # Use METIS ordering.
  depends_on "gcc"   => :build    # GNU Fortran is now provided as part of GCC
  depends_on "cmake" => :build

  def install
    args = ["-DCMAKE_INSTALL_PREFIX=#{prefix}",
            "-DBUILD_SHARED_LIBS=ON",
            "-DBUILD_DOCUMENTATION=OFF",
            "-DBUILD_64bits=OFF",
            "-DPASTIX_INT64=OFF",
            "-DPASTIX_ORDERING_SCOTCH=ON",
            "-DPASTIX_WITH_FORTRAN=ON",
            "-DPASTIX_WITH_MPI=OFF",
            "-DPASTIX_WITH_CUDA=OFF",
            "-DPASTIX_WITH_STARPU=OFF",
            "-DPASTIX_WITH_PARSEC=OFF"]
    args += ["-DPASTIX_ORDERING_METIS=ON"] if build.with? "metis"
    mkdir "build" do
      system "cmake", "..", *args
      system "make"
      system "make", "install"
    end
    pkgshare.install "example" # Contains all test programs.
  end

  def caveats; <<~EOS
    Set the PYTHONPATH environment variable:
      export PYTHONPATH=#{prefix}/lib/python:$PYTHONPATH
    Try python example with:
      python #{prefix}/examples/simple.py
    Or simple example with:
      #{prefix}/examples/simple -9 10:10:10
    EOS
  end

  test do
    system "#{prefix}/examples/simple", "-9", "10:10:10"
    ohai "All test output is in #{HOMEBREW_LOGS}/pastix. Please check."
  end
end
