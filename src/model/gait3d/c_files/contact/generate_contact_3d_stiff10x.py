"""Variant: 10x Hertz stiffness coefficient (F = 10 * C_hertz * d^1.5), default damping."""
from generate_contact_3d import generate_contact_3d

if __name__ == "__main__":
    generate_contact_3d(stiffness_scale=10.0, dissipation_val=2.0,
                         out_name='contact_3d_smoothsphere_stiff10x_al.c')
