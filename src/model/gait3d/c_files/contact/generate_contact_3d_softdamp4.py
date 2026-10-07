"""Variant: quarter Hunt-Crossley damping (dissipation 2.0 -> 0.5, HC bracket coeff 3.0 -> 0.75).

Note: bracket coeff 0.75 exactly matches the Nitschke contact model's damping
coefficient (see contact_3d.al: `(1 - 0.75*ycdot)`) -- this variant isolates
the stiffness-law difference (linear vs Hertzian) with damping held equal.
"""
from generate_contact_3d import generate_contact_3d

if __name__ == "__main__":
    generate_contact_3d(dissipation_val=0.5, stiffness_scale=1.0,
                         out_name='contact_3d_smoothsphere_softdamp4_al.c')
