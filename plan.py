# Create a pretty, multipage PDF summarizing the final directory structure
# and where each provided routine is used in the new OOP layout.

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
from textwrap import wrap

# ---------- Content assembly ----------

title = "PIC Code OOP Refactor — Final Directory Structure & Routine Mapping"

directory_structure = r"""
./
├── CMakeLists.txt
└── Src
    ├── main.f90
    │
    ├── modules/                      # Interfaces (module headers)
    │   ├── core/
    │   │   ├── mod_kinds.f90
    │   │   ├── mod_constants.f90
    │   │   ├── mod_logger.f90
    │   │   ├── mod_timer.f90
    │   │   ├── mod_rng.f90
    │   │   └── mod_config.f90
    │   ├── domain/
    │   │   ├── mod_geometry.f90
    │   │   └── mod_boundary.f90
    │   ├── field/
    │   │   ├── mod_field.f90
    │   │   ├── mod_bfield.f90
    │   │   └── mod_poisson.f90
    │   ├── particles/
    │   │   ├── mod_species.f90
    │   │   ├── mod_injector.f90
    │   │   ├── mod_collisions.f90
    │   │   └── mod_reactions.f90
    │   ├── io/
    │   │   ├── mod_io.f90
    │   │   └── mod_restart.f90
    │   ├── diag/
    │   │   └── mod_diagnostics.f90
    │   └── driver/
    │       ├── mod_pic.f90
    │       └── mod_factory.f90
    │
    └── submodules/                   # Implementations (submodules)
        ├── core/
        │   ├── submod_kinds.f90
        │   ├── submod_constants.f90
        │   ├── submod_logger.f90
        │   ├── submod_timer.f90
        │   ├── submod_rng.f90
        │   └── submod_config.f90
        ├── domain/
        │   ├── submod_geometry.f90
        │   └── submod_boundary.f90
        ├── field/
        │   ├── submod_field.f90
        │   ├── submod_bfield_gaussian.f90
        │   ├── submod_bfield_map.f90
        │   ├── submod_poisson_mg.f90
        │   ├── submod_poisson_sor.f90
        │   └── submod_poisson_jacobi.f90
        ├── particles/
        │   ├── submod_species.f90
        │   ├── submod_injector.f90
        │   ├── submod_collisions.f90
        │   └── submod_reactions.f90
        ├── io/
        │   ├── submod_io.f90
        │   └── submod_restart.f90
        ├── diag/
        │   └── submod_diagnostics.f90
        └── driver/
            ├── submod_pic.f90
            └── submod_factory.f90
"""

mapping_sections = [
("Bfield.f90",
"""read_Bfield_map              → modules/field/mod_bfield.f90  → submodules/field/submod_bfield_map.f90
gaussian_Bfield              → modules/field/mod_bfield.f90  → submodules/field/submod_bfield_gaussian.f90
PG_current                   → modules/field/mod_bfield.f90  → submodules/field/submod_bfield_gaussian.f90
EE_magnets                   → modules/field/mod_bfield.f90  → submodules/field/submod_bfield_gaussian.f90
Bfield_Hall_thruster         → modules/field/mod_bfield.f90  → submodules/field/submod_bfield_gaussian.f90
Bx_magnet, By_magnet         → modules/field/mod_bfield.f90  → submodules/field/submod_bfield_gaussian.f90
find_B_dir                   → modules/field/mod_bfield.f90  → (utility) implemented in submodules/field/*"""),

("calc_rho.f90",
"""calc_rho                     → modules/field/mod_field.f90   → submodules/field/submod_field.f90
dens_red                      → modules/field/mod_field.f90   → submodules/field/submod_field.f90"""),

("collisions.f90",
"""collisions                   → modules/particles/mod_collisions.f90 → submodules/particles/submod_collisions.f90
collision_OMP                 → modules/particles/mod_collisions.f90 → submodules/particles/submod_collisions.f90
scatter                       → modules/particles/mod_collisions.f90 → submodules/particles/submod_collisions.f90"""),

("Efield.f90",
"""calc_Efield                  → modules/field/mod_field.f90   → submodules/field/submod_field.f90"""),

("generate_boundary.f90",
"""generate_boundary            → modules/domain/mod_boundary.f90 → submodules/domain/submod_boundary.f90"""),

("indexx.f90",
"""indexx                       → modules/particles/mod_reactions.f90 (or utils) → submodules/particles/submod_reactions.f90"""),

("load_part.f90",
"""load_part                    → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90
load_part_OMP                 → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90
load_gauss                    → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90"""),

("mg.f90",
"""mg                           → modules/field/mod_poisson.f90 → submodules/field/submod_poisson_mg.f90
restriction, getres, prolongation → modules/field/mod_poisson.f90 → submodules/field/submod_poisson_mg.f90"""),

("part_expmover.f90",
"""part_mover                   → modules/particles/mod_species.f90 → submodules/particles/submod_species.f90
eheating                      → modules/particles/mod_species.f90 → submodules/particles/submod_species.f90
charge_deposition             → modules/particles/mod_species.f90 → submodules/particles/submod_species.f90"""),

("part_flux_injection.f90",
"""part_flux_injection          → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90
load_flux_OMP                 → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90"""),

("part_injection.f90",
"""part_injection               → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90
shifted_maxwellian_flux       → modules/particles/mod_injector.f90 → submodules/particles/submod_injector.f90"""),

("pdesolver.f90",
"""pdesolver                    → modules/field/mod_poisson.f90 → submodules/field/submod_poisson_sor.f90 (or wrapper to select solver)"""),

("ran2.f90",
"""ran2                         → modules/core/mod_rng.f90 → submodules/core/submod_rng.f90"""),

("read_input.f90",
"""read_input                   → modules/core/mod_config.f90 → submodules/core/submod_config.f90"""),

("read_reactions.f90",
"""read_reactions               → modules/particles/mod_reactions.f90 → submodules/particles/submod_reactions.f90
ordering                      → modules/particles/mod_reactions.f90 → submodules/particles/submod_reactions.f90"""),

("restart.f90",
"""restart                      → modules/io/mod_restart.f90 → submodules/io/submod_restart.f90
restart_OMP                   → modules/io/mod_restart.f90 → submodules/io/submod_restart.f90"""),

("sors.f90",
"""sor_rb                       → modules/field/mod_poisson.f90 → submodules/field/submod_poisson_sor.f90
jacobi                        → modules/field/mod_poisson.f90 → submodules/field/submod_poisson_jacobi.f90"""),

("sorting.f90",
"""part_sorting_OMP             → modules/particles/mod_species.f90 → submodules/particles/submod_species.f90
part_sorting                  → modules/particles/mod_species.f90 → submodules/particles/submod_species.f90
part_sorting_dom              → modules/particles/mod_species.f90 → submodules/particles/submod_species.f90"""),

("utils.f90",
"""part_moments                 → modules/diag/mod_diagnostics.f90 → submodules/diag/submod_diagnostics.f90
calc_avg                      → modules/diag/mod_diagnostics.f90 → submodules/diag/submod_diagnostics.f90
MSTIMER                       → modules/core/mod_timer.f90 → submodules/core/submod_timer.f90
stop_calculation              → modules/core/mod_logger.f90 → submodules/core/submod_logger.f90"""),

("write_data.f90",
"""write_data                   → modules/io/mod_io.f90 → submodules/io/submod_io.f90"""),
]

# ---------- PDF generation ----------

pdf_path = "./PIC_OOP_Directory_and_Mapping.pdf"
pp = PdfPages(pdf_path)

def add_page(title_text, body_text, fontsize=11, title_size=16):
    fig = plt.figure(figsize=(8.5, 11))  # US Letter portrait
    # Title
    plt.text(0.5, 0.96, title_text, ha='center', va='top', fontsize=title_size, fontweight='bold')
    # Body text
    # Wrap lines to fit; use a monospaced look by manual formatting
    left_margin = 0.06
    top = 0.92
    line_height = 0.020  # tweak for spacing
    max_chars = 100      # rough wrap width
    y = top
    for line in body_text.strip("\n").split("\n"):
        if not line:
            y -= line_height
            continue
        wrapped = wrap(line, width=max_chars, break_long_words=False, replace_whitespace=False, drop_whitespace=False)
        for w in wrapped:
            plt.text(left_margin, y, w, ha='left', va='top', fontsize=fontsize, family='monospace')
            y -= line_height
            if y < 0.06:  # new page
                plt.axis('off')
                pp.savefig(fig, bbox_inches='tight')
                plt.close(fig)
                fig = plt.figure(figsize=(8.5, 11))
                plt.text(0.5, 0.96, title_text + " (cont.)", ha='center', va='top', fontsize=title_size, fontweight='bold')
                y = top
    plt.axis('off')
    pp.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# Page 1: Title + Overview
overview = """
This document presents the final directory layout for the modern, object-oriented PIC code and a
one-to-one mapping of all your original routines to their new module/submodule homes.

Conventions:
- modules/ contain headers & public interfaces (types + procedure declarations).
- submodules/ contain the actual implementations (bodies) for those interfaces.
- Naming: module mod_X  ↔  submodule (mod_X) submod_X
"""

add_page(title, overview)

# Page 2: Directory tree
add_page("Final Directory Structure", directory_structure, fontsize=10)

# Subsequent pages: mappings by original file
for section_title, content in mapping_sections:
    add_page(f"Mapping — {section_title}", content, fontsize=11)

pp.close()

pdf_path
