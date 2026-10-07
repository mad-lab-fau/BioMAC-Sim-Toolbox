import sympy as sp
import os

def generate_contact_3d(dissipation_val=2.0, stiffness_scale=1.0, out_name='contact_3d_smoothsphere_al.c'):
    """Generate the smoothsphere contact model C code.

    dissipation_val: Hunt-Crossley dissipation coefficient c (s/m). Default 2.0
        matches the original model (HC bracket factor 1.5*c = 3.0).
    stiffness_scale: multiplier applied to the Hertz force coefficient C_hertz,
        i.e. F = stiffness_scale * C_hertz * d^1.5. Default 1.0 matches the
        original model.
    out_name: output C filename, written into this script's directory.
    """
    # Define symbolic variables
    xp, yp, zp = sp.symbols('contact->xp contact->yp contact->zp')

    # -----------------------------------------------------------------------
    # Physical parameters (matching PredSim's calcn geometry)
    # -----------------------------------------------------------------------
    # stiffness_val: plain-strain modulus E* (N/m^2) — as stored in the .osim.
    #   NOT a direct force coefficient. Must be converted via Hertz theory.
    stiffness_val     = 1000000.0   # N/m^2  (plain-strain modulus S = 1e6)
    # radius: sphere radius R = 0.032 m (CP marker must be placed at the sphere center y = -0.01 m)
    radius            = 0.032       # m
    dynamic_friction_val   = 0.8
    static_friction_val    = 0.8   # us = 0.8
    viscous_friction_val   = 0.5
    transition_velocity_val = 0.2  # m/s

    # Hertz coefficient C [N/m^1.5]:
    #   Simbody derivation: k = 0.5 * E*^(2/3),  C = (4/3) * k * sqrt(R * k)
    k_hertz = 0.5 * stiffness_val**(2.0 / 3.0)
    C_hertz_num = stiffness_scale * (4.0 / 3.0) * k_hertz * (radius * k_hertz)**0.5

    # Hunt-Crossley smoothing parameters (matching Simbody defaults)
    bv = 50.0   # hunt_crossley_smoothing = 50  (non-negativity tanh steepness)

    # -----------------------------------------------------------------------
    # Symbolic inputs
    # -----------------------------------------------------------------------
    fk    = sp.symbols('fk0:12')
    fkdot = sp.symbols('fkdot0:12')
    x     = sp.symbols('x0:6')
    xdot  = sp.symbols('xdot0:6')

    # rename for clarity matching contact_3d.al
    # Fx=x0, Fy=x1, Fz=x2, xc=x3, yc=x4, zc=x5
    Fx, Fy, Fz, xc, yc, zc = x[0], x[1], x[2], x[3], x[4], x[5]
    xcdot, ycdot, zcdot = xdot[3], xdot[4], xdot[5]

    # -----------------------------------------------------------------------
    # Kinematics: Ploc = R' * (Pglb - p)
    #   p = [fk0, fk1, fk2], R = [fk3,fk4,fk5; fk6,fk7,fk8; fk9,fk10,fk11]
    # -----------------------------------------------------------------------
    xcloc = fk[3]*(xc-fk[0]) + fk[6]*(yc-fk[1]) + fk[9]*(zc-fk[2])
    ycloc = fk[4]*(xc-fk[0]) + fk[7]*(yc-fk[1]) + fk[10]*(zc-fk[2])
    zcloc = fk[5]*(xc-fk[0]) + fk[8]*(yc-fk[1]) + fk[11]*(zc-fk[2])

    f1 = xcloc - xp
    f2 = ycloc - yp
    f3 = zcloc - zp

    # -----------------------------------------------------------------------
    # Normal force: Hertz + Hunt-Crossley, with non-negativity guard
    #
    # Fix 1 (stiffness): Use correct Hertz coefficient C_hertz (N/m^1.5)
    #   calculated from E* = 1e6 N/m^2 and R = 0.032 m, times stiffness_scale.
    #
    # Fix 2 (damping): Simbody Hunt-Crossley factor is (1 + 1.5*c*vn),
    #   not (1 + c*vn).
    #
    # Fix 3 (non-negativity): Without a guard, the force goes negative when
    #   ycdot > 1/(1.5*c) m/s, creating adhesive ground during push-off.
    #   Multiply by 0.5*(1+tanh(bv*(vn + 1/(1.5*c)))) to smoothly zero it at
    #   the point the damping bracket would cross zero (matches Simbody's
    #   hunt_crossley_smoothing parameter bv=50).
    # -----------------------------------------------------------------------
    eps_y = 0.001  # smoothing constant for penetration depth [m]

    # CP marker sits at the SPHERE CENTRE. f1-f3 pin (xc,yc,zc) to the centre,
    # and R - yc is exactly Simbody's indentation for a ground plane at y = 0.
    a  = radius - yc
    d  = sp.Rational(1, 2) * (a + sp.sqrt(a**2 + eps_y**2))   # smooth max(a, 0)
    vn = -ycdot  # indentation rate: positive = compressing


    # Bodyweight symbol (provided at runtime by get_bodyweight())
    bw = sp.symbols('bw')

    # Hertz force
    fh = C_hertz_num * d**sp.Rational(3, 2)

    # Hunt-Crossley damping: factor is 1.5*c (Fix 2)
    fhc = fh * (1.0 + 1.5 * dissipation_val * vn)

    # Non-negativity guard: zero force before damping bracket goes negative (Fix 3)
    #   The damping bracket (1 + 1.5*c*vn) crosses zero at vn = -1/(1.5*c) = -1/3 m/s,
    #   i.e. when ycdot = +1/3 m/s (foot lifting off).
    #   Guard = 0.5*(1 + tanh(bv*(vn + 1/(1.5*c)))) which → 1 during compression
    #   and → 0 during rapid liftoff, smoothly centred at vn = -1/(1.5*c).
    guard_zero_vn = -1.0 / (1.5 * dissipation_val)   # vn where bracket = 0
    fhc = fhc * sp.Rational(1, 2) * (1.0 + sp.tanh(bv * (vn - guard_zero_vn)))

    # Bodyweight-normalised residual (forces in BW, matching BioMAC convention)
    f4 = Fy - fhc / bw

    # -----------------------------------------------------------------------
    # Friction force: vectorial Stribeck friction model (matching Simbody/OpenSim)
    #   Includes rotational angular velocity contribution (omega x r) to slip velocity.
    # -----------------------------------------------------------------------
    # angular velocity of the segment in ground: omega_hat = skew(Rdot * R^T)
    Rm = sp.Matrix(3, 3, [fk[3], fk[4], fk[5], fk[6], fk[7], fk[8], fk[9], fk[10], fk[11]])
    Rd = sp.Matrix(3, 3, [fkdot[3], fkdot[4], fkdot[5], fkdot[6], fkdot[7], fkdot[8],
                          fkdot[9], fkdot[10], fkdot[11]])
    W  = sp.Rational(1, 2) * (Rd*Rm.T - (Rd*Rm.T).T)   # antisymmetrised, robust to drift in R
    wx, wz = W[2, 1], W[1, 0]

    # station sits at (0, -(R - d/2), 0) from the centre in GROUND coords, as in Simbody
    lever = radius - d/2
    vx = xcdot + lever*wz
    vz = zcdot - lever*wx

    v_slide = sp.sqrt(vx**2 + vz**2 + 1e-5)   # cf = constant_contact_force = 1e-5
    vrel    = v_slide / transition_velocity_val
    # gate: smooth min(vrel, 1.0) to scale down friction coefficient at low velocities
    gate    = 0.5 * (vrel + 1.0 - sp.sqrt((vrel - 1.0)**2 + 1e-6))

    # mu: friction coefficient containing dynamic, static (Stribeck), and viscous friction
    mu = gate * (dynamic_friction_val
                 + 2.0 * (static_friction_val - dynamic_friction_val) / (1.0 + vrel**2)) \
         + viscous_friction_val * v_slide

    f5 = Fx + Fy * mu * (vx / v_slide)
    f6 = Fz + Fy * mu * (vz / v_slide)

    f = [f1, f2, f3, f4, f5, f6]

    # -----------------------------------------------------------------------
    # Compute Jacobians
    # -----------------------------------------------------------------------
    df_dfk    = [[sp.diff(fi, fkj)    for fkj    in fk]   for fi in f]
    df_dfkdot = [[sp.diff(fi, fkj)    for fkj    in fkdot] for fi in f]
    df_dx     = [[sp.diff(fi, xj)     for xj     in x]    for fi in f]
    df_dxdot  = [[sp.diff(fi, xdotj)  for xdotj  in xdot] for fi in f]

    # -----------------------------------------------------------------------
    # Generate C code
    # -----------------------------------------------------------------------
    c_code = []
    c_code.append("// This file was generated by generate_contact_3d.py using SymPy\n")
    c_code.append("// Physics: Hertz+Hunt-Crossley with non-negativity guard (fixed vs original)")
    c_code.append(f"//   C_hertz = {C_hertz_num:.6e} N/m^1.5  (R={radius} m, S={stiffness_val}, stiffness_scale={stiffness_scale})")
    guard_zero_vn_val = -1.0 / (1.5 * dissipation_val)
    c_code.append(f"//   dissipation_val={dissipation_val}, HC factor: 1.5*c={1.5*dissipation_val}, guard zeroes at vn={guard_zero_vn_val:.4f} m/s (ycdot={-guard_zero_vn_val:.4f} m/s)\n")
    c_code.append("#include <math.h>")
    c_code.append('#include "gait3d_smoothsphere_contact.h"\n')
    c_code.append("double get_bodyweight();\n")
    c_code.append("void contact_al(contactprop* contact,")
    c_code.append("\tdouble fk[12], double fkdot[12], double x[NCVAR], double xdot[NCVAR],")
    c_code.append("\tdouble f[NCF], double df_dfk[NCF][12], double df_dfkdot[NCF][12],")
    c_code.append("\tdouble df_dx[NCF][NCVAR], double df_dxdot[NCF][NCVAR]) {\n")

    # Variable declarations
    for i in range(12):
        c_code.append(f"\tdouble fk{i} = fk[{i}];")
        c_code.append(f"\tdouble fkdot{i} = fkdot[{i}];")
    for i in range(6):
        c_code.append(f"\tdouble x{i} = x[{i}];")
        c_code.append(f"\tdouble xdot{i} = xdot[{i}];")

    c_code.append("\tdouble bw = get_bodyweight();")
    c_code.append("")

    # Residuals
    for i in range(6):
        c_code.append(f"\tf[{i}] = {sp.ccode(f[i])};")

    c_code.append("")
    # df/dfk
    for i in range(6):
        for j in range(12):
            expr = df_dfk[i][j]
            c_code.append(f"\tdf_dfk[{i}][{j}] = {sp.ccode(expr) if expr != 0 else '0.0'};")

    c_code.append("")
    # df/dfkdot
    for i in range(6):
        for j in range(12):
            expr = df_dfkdot[i][j]
            c_code.append(f"\tdf_dfkdot[{i}][{j}] = {sp.ccode(expr) if expr != 0 else '0.0'};")

    c_code.append("")
    # df/dx
    for i in range(6):
        for j in range(6):
            expr = df_dx[i][j]
            c_code.append(f"\tdf_dx[{i}][{j}] = {sp.ccode(expr) if expr != 0 else '0.0'};")

    c_code.append("")
    # df/dxdot
    for i in range(6):
        for j in range(6):
            expr = df_dxdot[i][j]
            c_code.append(f"\tdf_dxdot[{i}][{j}] = {sp.ccode(expr) if expr != 0 else '0.0'};")

    c_code.append("}")

    # Write output
    out_dir  = os.path.dirname(os.path.abspath(__file__))
    out_path = os.path.join(out_dir, out_name)
    with open(out_path, "w") as out:
        out.write("\n".join(c_code) + "\n")
    print(f"Generated {out_path}")

if __name__ == "__main__":
    generate_contact_3d()
