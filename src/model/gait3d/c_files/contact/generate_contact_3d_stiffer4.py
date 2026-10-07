"""Variant: quadruple Hertz stiffness coefficient (F = 4 * C_hertz * d^1.5)."""
from generate_contact_3d import generate_contact_3d

if __name__ == "__main__":
    generate_contact_3d(dissipation_val=2.0, stiffness_scale=4.0,
                         out_name='contact_3d_smoothsphere_stiffer4_al.c')
