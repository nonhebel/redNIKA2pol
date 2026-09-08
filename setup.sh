
#!/bin/bash
set -e

SOURCE="W3IRS4" # enter source name here

# For safe rerunning
relink() {
    rm -f "$2"
    ln -s "$1" "$2"
}

# Copy only if destination doesn't already exist, so re-running setup.sh
# never silently overwrites edits made inside a target's scripts.
safe_copy() {
    local src="$1"
    local dest="$2"
    if [ -f "$dest" ]; then
        echo "Skipping $dest (already exists)"
    else
        cp "$src" "$dest"
    fi
}

setup_target() {
    local target="$1"
    local raw_dir="$2"

    # 01_general
    mkdir -p "reductions/$target/01_general/red"
    mkdir -p "reductions/$target/01_general/stat"
    relink "$(pwd)/$raw_dir" "reductions/$target/01_general/imbfitsDir"
    safe_copy "scripts/01_general/source_a13_MP.piic" "reductions/$target/01_general/${target}_a13_MP.piic"
    safe_copy "scripts/01_general/source_setMPpar.piic" "reductions/$target/01_general/${target}_setMPpar.piic"
    safe_copy "scripts/01_general/check_reduction_convergence.py" "reductions/$target/01_general/"
    
    # 02_individual
    mkdir -p "reductions/$target/02_individual/red"
    mkdir -p "reductions/$target/02_individual/stat"
    relink "$(pwd)/$raw_dir" "reductions/$target/02_individual/imbfitsDir"
    relink "$(pwd)/reductions/$target/01_general/source.pol" "reductions/$target/02_individual/source.pol"
    relink "$(pwd)/reductions/$target/01_general/base.pol" "reductions/$target/02_individual/base.pol"
    relink "$(pwd)/reductions/$target/01_general/${target}_a13.LIST" "reductions/$target/02_individual/${target}_a13.LIST"
    safe_copy "$(pwd)/scripts/02_individual/source_a13_save_ind.piic" "reductions/$target/02_individual/${target}_a13_save_ind.piic"

    # 03_corr
    mkdir -p "reductions/$target/03_corr/red"
    relink "$(pwd)/reductions/$target/02_individual/red" "reductions/$target/03_corr/red_uncorr"
    relink "$(pwd)/raw/calibrators" "reductions/$target/03_corr/red_calibrators"
    safe_copy "$(pwd)/scripts/03_corr/combine_maps.piic" "reductions/$target/03_corr/combine_maps.piic"
    safe_copy "$(pwd)/scripts/03_corr/corr_instr_pol.piic" "reductions/$target/03_corr/corr_instr_pol.piic"
    safe_copy "$(pwd)/scripts/03_corr/source_ip_corr_MP.py" "reductions/$target/03_corr/${target}_ip_corr_MP.py"

    # 04_bootstrap
    mkdir -p "reductions/$target/04_bootstrap"
    relink "$(pwd)/reductions/$target/03_corr/red" "reductions/$target/04_bootstrap/red"
    relink "$(pwd)/reductions/$target/03_corr/all_corr_maps.LIST" "reductions/$target/04_bootstrap/all_maps.LIST"
    relink "$(pwd)/reductions/$target/03_corr/all_rgw.LIST" "reductions/$target/04_bootstrap/all_rgw.LIST"
    safe_copy "$(pwd)/scripts/04_bootstrap/combine_maps.piic" "reductions/$target/04_bootstrap/combine_maps.piic"
    safe_copy "$(pwd)/scripts/04_bootstrap/bsm.piic" "reductions/$target/04_bootstrap/bsm.piic"
    safe_copy "$(pwd)/scripts/04_bootstrap/source_bootstrap_err.py" "reductions/$target/04_bootstrap/${target}_bootstrap_err.py"

    # 05_products
    mkdir -p "reductions/$target/05_products"
    safe_copy "$(pwd)/scripts/05_products/bsm.piic" "reductions/$target/05_products/bsm.piic"
    safe_copy "$(pwd)/scripts/05_products/creaRMS_IQU.piic" "reductions/$target/05_products/creaRMS_IQU.piic"
    safe_copy "$(pwd)/scripts/05_products/source_create_products.py" "reductions/$target/05_products/${target}_create_products.py"
    echo "Set up $target"
}

setup_target "$SOURCE" "raw/science"
