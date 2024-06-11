{
  description = "OpenLB";

  inputs = {
    nixpkgs.url = github:NixOS/nixpkgs/nixos-22.11;
  };

  outputs = { self, nixpkgs, ... }: let
    system = "x86_64-linux";
    pkgs = import nixpkgs {
      inherit system;
      config.allowUnfree = true;
    };

    common-env = with pkgs;  [
    # common make dependencies
      gnumake

    # debugging
      gdb
      valgrind

    # result presentation
      gnuplot

    # external dependencies
      tinyxml
      zlib
      gtest
    ];

  in {
    packages.${system} = {
      userguide = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-userguide";
        src = ./.;
        buildInputs = with pkgs; let
          custom-texlive = pkgs.texlive.combine {
            inherit (pkgs.texlive) scheme-small collection-langgerman latexmk xpatch xstring siunitx biblatex logreq palatino courier mathpazo helvetic multirow elsarticle widetable makecell pgfplots spath3 placeins abstract tocloft;
          };

        in [
          gnumake
          biber
          custom-texlive
        ];
        buildPhase = ''
          make userguide
        '';
        installPhase = ''
          mkdir -p $out
          cp doc/userGuide/olb-ug.pdf $out/
        '';
      };

      doxygen = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-doxygen";
        src = ./.;
        buildInputs = with pkgs; [
          gnumake
          doxygen
          graphviz
        ];
        buildPhase = ''
          make doxygen
        '';
        installPhase = ''
          mkdir -p $out
          cp -r doc/doxygen/html $out/
        '';
      };

      env-gcc = pkgs.mkShell {
        name = "openlb-env-gcc";
        buildInputs = common-env ++ (with pkgs; [
          gcc11
        ]);
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++17"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-gcc-openmpi = pkgs.mkShell {
        name = "openlb-env-gcc-openmpi";
        buildInputs = common-env ++ (with pkgs; [
          gcc11
          openmpi
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++17"

          export PARALLEL_MODE=MPI

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-clang = pkgs.mkShell {
        name = "openlb-env-clang";
        buildInputs = common-env ++ (with pkgs; [
          clang_12
          llvmPackages_12.bintools-unwrapped
          llvmPackages_12.openmp
        ]);
        shellHook = ''
          export CXX=clang++
          export CC=clang

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++17"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD"
        '';
      };

      env-cuda = pkgs.mkShell {
        name = "openlb-env-cuda";
        buildInputs = common-env ++ (with pkgs; [
          cudatoolkit_11
        ]);
        shellHook = ''
          export CXX=nvcc
          export CC=nvcc

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      env-gcc-cuda = pkgs.mkShell {
        name = "openlb-env-gcc-cuda";
        buildInputs = common-env ++ (with pkgs; [
          gcc11
          cudatoolkit_11
        ]);
        shellHook = ''
          export CXX=g++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++17"

          export PARALLEL_MODE=NONE

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++17"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      env-gcc-openmpi-cuda = pkgs.mkShell {
        name = "openlb-env-gcc-openmpi-cuda";
        buildInputs = common-env ++ (with pkgs; [
          gcc11
          cudatoolkit_11
          (openmpi.override {
            cudaSupport = true;
            cudatoolkit = cudatoolkit_11;
          })
        ]);
        shellHook = ''
          export CXX=mpic++
          export CC=gcc

          export CXXFLAGS="-O3 -Wall -march=native -mtune=native -std=c++17"

          export PARALLEL_MODE=MPI

          export PLATFORMS="CPU_SISD GPU_CUDA"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++17"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      env-heterogeneity = pkgs.mkShell {
        name = "openlb-heterogeneity";
        buildInputs = common-env ++ (with pkgs; [
          gcc11
          cudatoolkit_11
          (openmpi.override {
            cudaSupport = true;
            cudatoolkit = cudatoolkit_11;
          })
        ]);
        shellHook = let
          arch = "tigerlake"; # match target CPU as closely as possible
        in ''
          export CXX=mpic++
          export CC=gcc

          export PLATFORMS="CPU_SISD CPU_SIMD GPU_CUDA"

          export CXXFLAGS="-O3 -Wall -march=${arch} -mtune=${arch} -std=c++17"

          export PARALLEL_MODE=HYBRID
          export OMPFLAGS="-fopenmp"

          export CUDA_CXX=nvcc
          export CUDA_CXXFLAGS="-O3 -std=c++17 --forward-unknown-to-host-compiler"
          export CUDA_LDFLAGS="-L/run/opengl-driver/lib"
          # try to auto-fill CUDA_ARCH for first GPU (override when in doubt)
          export CUDA_ARCH=$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader | head -n 1 | grep -o [0-9] | tr -d '\n')

          export FLOATING_POINT_TYPE=float

          export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/run/opengl-driver/lib
        '';
      };

      ## Bootstrap persistent conda environment for code generation (impure)
      # (to be replaced once cppyy can be consistently built using Nix's python facilities)
      env-generate-code = pkgs.stdenvNoCC.mkDerivation rec {
        name = "openlb-env-generate-code";
        buildInputs = [ pkgs.micromamba ];
        buildCommand = ''
          export MAMBA_ROOT_PREFIX=$out/
          mkdir -p $MAMBA_ROOT_PREFIX
          ## Create environment either from yml (fixed version) or explicit spec file (fixed hashes)
          micromamba create --cacert-path ${pkgs.cacert}/etc/ssl/certs/ca-bundle.crt --yes -f ${./script/codegen/environment.yml}
        '';
      };

      # Call "make cse" in code generation environment
      generate-code = pkgs.buildFHSUserEnv {
        name = "openlb-generate-code";
        targetPkgs = pkgs: [
          pkgs.micromamba
          (pkgs.writeScriptBin "activate-codegen-env" ''
            #!/usr/bin/env bash
            # Enable micromamba shell without completion (preventing error)
            eval "$(micromamba shell hook --shell=bash)"
            # Use pre-generated environment (unsandboxed)
            export MAMBA_ROOT_PREFIX=${self.packages.${system}.env-generate-code}
            export CLING_STANDARD_PCH=/tmp/cling-pch
            micromamba activate openlb-codegen
            make cse
          '')
        ];
        runScript = "/usr/bin/activate-codegen-env";
      };
    };
  };
}
