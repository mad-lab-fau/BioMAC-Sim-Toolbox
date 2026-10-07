"""Variant: effectively zero Hunt-Crossley damping (dissipation_val~0; exactly 0.0 would divide-by-zero in the non-negativity guard, see generate_contact_3d.py)."""
from generate_contact_3d import generate_contact_3d

if __name__ == "__main__":
    generate_contact_3d(stiffness_scale=1.0, dissipation_val=1e-6,
                         out_name='contact_3d_smoothsphere_nodamp_al.c')
