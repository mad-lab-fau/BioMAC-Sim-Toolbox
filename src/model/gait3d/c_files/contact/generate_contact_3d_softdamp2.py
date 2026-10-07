"""Variant: half Hunt-Crossley damping (dissipation 2.0 -> 1.0, HC bracket coeff 3.0 -> 1.5)."""
from generate_contact_3d import generate_contact_3d

if __name__ == "__main__":
    generate_contact_3d(dissipation_val=1.0, stiffness_scale=1.0,
                         out_name='contact_3d_smoothsphere_softdamp2_al.c')
