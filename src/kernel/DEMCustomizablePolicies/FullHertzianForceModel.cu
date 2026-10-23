// DEM force calculation strategies, modifiable

// No need to do any contact force calculation if no contact. And it can happen,
// when we added extra contact margins
if (overlapDepth > 0) {
    // Promote inputs before arithmetic; mixed float3/double3 operators can narrow results.
    // Keep force and rolling resistance in double until the existing float output boundary.
    const double3 normal = to_double3(B2A);
    double3 force_double = to_double3(force);
    double3 torque_double = to_double3(torque_only_force);
    double contact_time = static_cast<double>(delta_time);

    // Material properties
    double E_cnt, G_cnt, CoR_cnt, mu_cnt, Crr_cnt;
    {
        // E and nu are associated with each material, so obtain them this way
        double E_A = E[bodyAMatType];
        double nu_A = nu[bodyAMatType];
        double E_B = E[bodyBMatType];
        double nu_B = nu[bodyBMatType];
        matProxy2ContactParam<double>(E_cnt, G_cnt, E_A, nu_A, E_B, nu_B);
        // CoR, mu and Crr are pair-wise, so obtain them this way
        CoR_cnt = CoR[bodyAMatType][bodyBMatType];
        mu_cnt = mu[bodyAMatType][bodyBMatType];
        Crr_cnt = Crr[bodyAMatType][bodyBMatType];
    }

    // We also need the relative velocity between A and B in global frame to use in the damping terms
    // To get that, we need contact points' rotational velocity in GLOBAL frame
    // This is local rotational velocity (the portion of linear vel contributed by rotation)
    double3 rotVelCPA = cross(to_double3(ARotVel), to_double3(locCPA));
    double3 rotVelCPB = cross(to_double3(BRotVel), to_double3(locCPB));
    applyOriQToVector3(rotVelCPA, make_double4(AOriQ.x, AOriQ.y, AOriQ.z, AOriQ.w));
    applyOriQToVector3(rotVelCPB, make_double4(BOriQ.x, BOriQ.y, BOriQ.z, BOriQ.w));

    // A few re-usables
    double mass_eff, beta, cnt_rad;
    double3 vrel_tan;
    double3 delta_tan = make_double3(delta_tan_x, delta_tan_y, delta_tan_z);

    // Normal force part
    {
        // The (total) relative linear velocity of A relative to B
        const double3 velB2A = (to_double3(ALinVel) + rotVelCPA) - (to_double3(BLinVel) + rotVelCPB);
        const double projection = dot(velB2A, normal);
        vrel_tan = velB2A - projection * normal;

        // Now we already have sufficient info to update contact history
        {
            delta_tan += static_cast<double>(ts) * vrel_tan;
            const double disp_proj = dot(delta_tan, normal);
            delta_tan -= disp_proj * normal;
            contact_time += static_cast<double>(ts);
        }

        mass_eff = (static_cast<double>(AOwnerMass) * static_cast<double>(BOwnerMass)) /
                   (static_cast<double>(AOwnerMass) + static_cast<double>(BOwnerMass));

        // Contact radius (radial distance from contact center axis) called cnt_rad, computed from area
        cnt_rad = sqrt(overlapArea / deme::PI);
        const double Sn = 2. * E_cnt * cnt_rad;

        const double loge = (CoR_cnt < DEME_TINY_FLOAT) ? log(static_cast<double>(DEME_TINY_FLOAT)) : log(CoR_cnt);
        beta = loge / sqrt(loge * loge + deme::PI_SQUARED);

        const double k_n = (2. / 3.) * Sn;
        const double gamma_n = (2. * sqrt(5. / 6.)) * beta * sqrt(Sn * mass_eff);

        force_double += (k_n * overlapDepth + gamma_n * projection) * normal;
        // printf("normal force: %f, %f, %f\n", force_double.x, force_double.y, force_double.z);
    }

    // Rolling resistance part
    if (Crr_cnt > 0.) {
        // Figure out if we should apply rolling resistance force
        bool should_add_rolling_resistance = true;
        {
            // Use geometric equivalence R_eff ≈ cnt_rad^2 / overlapDepth
            const double R_eff = (cnt_rad * cnt_rad) / overlapDepth;

            const double kn_simple = (4. / 3.) * E_cnt * sqrt(R_eff);
            const double gn_simple = -2. * sqrt((5. / 3.) * mass_eff * E_cnt) * beta * pow(R_eff, 0.25);

            const double d_coeff = gn_simple / (2. * sqrt(kn_simple * mass_eff));

            if (d_coeff < 1.) {
                double t_collision = deme::PI * sqrt(mass_eff / (kn_simple * (1. - d_coeff * d_coeff)));
                if (contact_time <= t_collision) {
                    should_add_rolling_resistance = false;
                }
            }
        }
        // If should, then compute it (using Schwartz model)
        if (should_add_rolling_resistance) {
            // Tangential velocity (only rolling contribution) of B relative to A, at contact point, in global
            const double3 v_rot = rotVelCPB - rotVelCPA;
            // This v_rot is only used for identifying resistance direction
            const double v_rot_mag = length(v_rot);
            if (v_rot_mag > DEME_TINY_FLOAT) {
                // You should know that Crr_cnt * normal_force is the underlying formula, and in our model,
                // it is a `force' that produces torque only, instead of also cancelling out friction.
                // Its direction is that it `resists' rotation, see picture in
                // https://en.wikipedia.org/wiki/Rolling_resistance.
                torque_double = (v_rot / v_rot_mag) * (Crr_cnt * length(force_double));
                // printf("torque force: %f, %f, %f\n", torque_double.x, torque_double.y, torque_double.z);
            }
        }
    }

    // Tangential force part
    if (mu_cnt > 0.) {
        const double kt = 8. * G_cnt * cnt_rad;
        const double gt = -2. * sqrt(5. / 3.) * beta * sqrt(mass_eff * kt);
        double3 tangent_force = -kt * delta_tan - gt * vrel_tan;
        const double ft = length(tangent_force);
        if (ft > DEME_TINY_FLOAT) {
            // Reverse-engineer to get tangential displacement
            const double ft_max = length(force_double) * mu_cnt;
            if (ft > ft_max) {
                tangent_force = (ft_max / ft) * tangent_force;
                delta_tan = (tangent_force + gt * vrel_tan) / (-kt);
            }
        } else {
            tangent_force = make_double3(0, 0, 0);
        }
        // Use force to collect tangent_force
        force_double += tangent_force;
        // printf("tangent force: %f, %f, %f\n", tangent_force.x, tangent_force.y, tangent_force.z);
    }

    // Finally, make sure we update those wildcards (in this case, contact history)
    // The solver stores forces and contact wildcards as floats; narrow only at writeback.
    force = to_float3(force_double);
    torque_only_force = to_float3(torque_double);
    delta_time = static_cast<float>(contact_time);
    delta_tan_x = static_cast<float>(delta_tan.x);
    delta_tan_y = static_cast<float>(delta_tan.y);
    delta_tan_z = static_cast<float>(delta_tan.z);
} else {
    // This is to be more rigorous. If in fact no physical contact, then contact wildcards (such as contact history)
    // should be cleared (they will also be automatically cleared if the contact is no longer detected).
    delta_time = 0;
    delta_tan_x = 0;
    delta_tan_y = 0;
    delta_tan_z = 0;
}
