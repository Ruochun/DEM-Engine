// DEM force calculation strategies, modifiable

if (overlapDepth > 0) {
    // Promote before vector arithmetic and quaternion products, then narrow only at writeback.
    const double3 normal = to_double3(B2A);
    double3 force_double = to_double3(force);

    // Material properties
    double E_cnt, CoR_cnt;
    {
        matProxy2ContactParam<double>(E_cnt, E[bodyAMatType], nu[bodyAMatType], E[bodyBMatType], nu[bodyBMatType]);
        // CoR is pair-wise, so obtain it this way
        CoR_cnt = CoR[bodyAMatType][bodyBMatType];
    }

    // We also need the relative velocity between A and B in global frame to use in the damping terms
    // To get that, we need contact points' rotational velocity in GLOBAL frame
    // This is local rotational velocity (the portion of linear vel contributed by rotation)
    double3 rotVelCPA = cross(to_double3(ARotVel), to_double3(locCPA));
    double3 rotVelCPB = cross(to_double3(BRotVel), to_double3(locCPB));
    // This is mapping from local rotational velocity to global
    applyOriQToVector3(rotVelCPA, make_double4(AOriQ.x, AOriQ.y, AOriQ.z, AOriQ.w));
    applyOriQToVector3(rotVelCPB, make_double4(BOriQ.x, BOriQ.y, BOriQ.z, BOriQ.w));

    // The (total) relative linear velocity of A relative to B
    const double3 velB2A = (to_double3(ALinVel) + rotVelCPA) - (to_double3(BLinVel) + rotVelCPB);
    const double projection = dot(velB2A, normal);

    const double mass_eff = (static_cast<double>(AOwnerMass) * static_cast<double>(BOwnerMass)) /
                            (static_cast<double>(AOwnerMass) + static_cast<double>(BOwnerMass));

    // Contact radius (radial distance from contact center axis) called cnt_rad, computed from area
    const double cnt_rad = sqrt(overlapArea / deme::PI);
    const double Sn = 2. * E_cnt * cnt_rad;

    const double loge = (CoR_cnt < DEME_TINY_FLOAT) ? log(static_cast<double>(DEME_TINY_FLOAT)) : log(CoR_cnt);
    const double beta = loge / sqrt(loge * loge + deme::PI_SQUARED);

    const double k_n = (2. / 3.) * Sn;
    const double gamma_n = (2. * sqrt(5. / 6.)) * beta * sqrt(Sn * mass_eff);

    force_double += (k_n * overlapDepth + gamma_n * projection) * normal;
    force = to_float3(force_double);
}
