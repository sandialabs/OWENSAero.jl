_lb_wrap_to_pi(alpha) = mod(alpha + pi, 2*pi) - pi

function _lb_ideal_cl(alpha, alpha_zero_lift, lift_slope, reference_flag)
    alpha_id = _lb_wrap_to_pi(alpha - alpha_zero_lift)
    direction = 1.0
    if alpha_id > pi / 2
        alpha_id = pi - alpha_id
        direction = -1.0
    elseif alpha_id < -pi / 2
        alpha_id = -pi - alpha_id
        direction = -1.0
    end

    if reference_flag == 1
        return direction * lift_slope * alpha_id
    end

    cutoff_alpha = 30.0 * pi / 180.0
    d1 = 1.8
    ideal_slope = 2.0 * pi
    if abs(alpha_id) < cutoff_alpha
        return direction * ideal_slope * alpha_id
    end
    return direction * (
        ideal_slope * (cutoff_alpha - sin(d1 * cutoff_alpha) / d1) +
        ideal_slope * sin(d1 * alpha_id) / d1
    )
end

"""
    Leishman_Beddoes(af, alpha75, alpha50, ds, mach, reynolds, lift_slope,
                     alpha_zero_lift, cl_critical_positive,
                     cl_critical_negative, cl_ref_le_last, cl_ref_last,
                     fstat_last, cv_last, dp, df, dcnv, le_separation_state,
                     s_lev; Tp=1.7, family_factor=0.0)

Advance one section of the incompressible Leishman-Beddoes dynamic-stall model.
Angles are radians. `alpha75` is the three-quarter-chord angle of attack used
for the reference and static airfoil response, and `alpha50` is the half-chord
angle used to project leading-edge vortex lift and drag. `ds = 2 U dt / c` is
the nondimensional time step.

The airfoil callback must accept `af(alpha, reynolds, mach, family_factor)` and
return `(cl, cd, cm)`. The return tuple contains `(cl, cd, cm)` followed by the
updated state variables in the same order as the state inputs, including the
updated leading-edge vortex age `s_lev`.
"""
function Leishman_Beddoes(
    af,
    alpha75,
    alpha50,
    ds,
    mach,
    reynolds,
    lift_slope,
    alpha_zero_lift,
    cl_critical_positive,
    cl_critical_negative,
    cl_ref_le_last,
    cl_ref_last,
    fstat_last,
    cv_last,
    dp,
    df,
    dcnv,
    le_separation_state,
    s_lev;
    Tp = 1.7,
    family_factor = 0.0,
)
    drag_separation_factor = 0.1

    cl_static, cd_static, _ = af(alpha75, reynolds, mach, family_factor)

    cl_ref = _lb_ideal_cl(alpha75, alpha_zero_lift, lift_slope, 1)
    cl_id = _lb_ideal_cl(alpha75, alpha_zero_lift, lift_slope, 0)

    transition = cos(alpha75 - alpha_zero_lift)^2
    dcl_ref_le = transition * dp
    dalpha_ref_le = dcl_ref_le / lift_slope

    cl_ref_le = cl_ref - dcl_ref_le
    cl_rate_flag = cl_ref_le * (cl_ref_le - cl_ref_le_last) > 0 ? 1 : 0
    alpha_ref_le = _lb_wrap_to_pi(alpha75 - dalpha_ref_le)

    cl_static_f, _, cm = af(alpha_ref_le, reynolds, mach, family_factor)
    cl_id_f = _lb_ideal_cl(alpha_ref_le, alpha_zero_lift, lift_slope, 0)
    cl_ratio_f = abs(cl_id_f) < 0.001 ? 999.0 : cl_static_f / cl_id_f

    fstat = cl_ratio_f > 0.25 ? min((sqrt(4.0 * cl_ratio_f) - 1.0)^2, 1.0) : 0.0
    f_dynamic = clamp(fstat - df, 0.0, 1.0)

    cl_ratio = abs(cl_id) < 0.001 ? 999.0 : cl_static / cl_id
    cl_id_used = cl_ratio > 1.0 ? cl_static : cl_id
    cl_sep = cl_ratio > 0.25 ? cl_id_used / 4.0 : cl_static
    cl_f = cl_sep + cl_id_used * 0.25 * (f_dynamic + 2.0 * sqrt(f_dynamic))
    dcd_f = drag_separation_factor * (cl_static - cl_f) * sign(cl_static)

    dcl_v = dcnv * cos(alpha50)
    dcd_v = dcnv * sin(alpha50)

    cv = cl_id_used - cl_f
    dcv = cv - cv_last
    cutoff_alpha = 50.0 * pi / 180.0
    if sign(dcv * cl_ref_le) < 0 || abs(alpha75 - alpha_zero_lift) > cutoff_alpha ||
       cl_rate_flag == 0
        dcv = 0.0
    end

    cl = cl_f + dcl_v
    cd = cd_static + dcd_f + dcd_v

    tf_ref = 3.0
    tv_ref = 6.0
    tv_l = 11.0

    if le_separation_state == 0 &&
       (cl_ref_le > cl_critical_positive || cl_ref_le < cl_critical_negative)
        le_separation_state = 1
        s_lev = 0.0
    elseif le_separation_state == 1 &&
           (cl_critical_negative < cl_ref_le < cl_critical_positive)
        le_separation_state = 0
        s_lev = 0.0
    end

    if le_separation_state == 1
        if s_lev < tv_l
            if cl_rate_flag > 0
                tf = 4.0 * tf_ref
                tv = tv_ref
            else
                tf = tf_ref / 2.0
                tv = tv_ref / 2.0
            end
        elseif s_lev < 2.0 * tv_l
            if cl_rate_flag > 0
                tf = 2.0 * tf_ref
                tv = tv_ref
            else
                tf = tf_ref / 2.0
                tv = tv_ref / 2.0
            end
        else
            tf = tf_ref
            tv = tv_ref
        end
    else
        tf = f_dynamic > 0.7 ? tf_ref : 2.0 * tf_ref
        tv = tv_ref
    end

    if le_separation_state == 1
        s_lev += ds
    end

    dp_new = dp * exp(-ds / Tp) + (cl_ref - cl_ref_last) * exp(-ds / (2 * Tp))
    df_new = df * exp(-ds / tf) + (fstat - fstat_last) * exp(-ds / (2 * tf))
    dcnv_new = dcnv * exp(-ds / tv) + dcv * exp(-ds / (2 * tv))

    return cl,
    cd,
    cm,
    cl_ref_le,
    cl_ref,
    fstat,
    cv,
    dp_new,
    df_new,
    dcnv_new,
    le_separation_state,
    s_lev
end

