// Unlike the dedicated collection kernel, this policy consumes locCPA/locCPB directly from the force-calculation
// kernel. Those values are computed locally for active contacts and explicitly zeroed with the forces for inactive
// margin candidates, so this route cannot read an unwritten global contact-point slot. A separate zero-resultant guard
// is therefore unnecessary for the stochastic contact-ramp fix.
float3 forceA = force + A_only_force;
float3 forceB = -force + B_only_force;
float3 torqueA = A_only_torque;
float3 torqueB = B_only_torque;

// Take care of A
{
    const bool bad_cp = !isfinite3(locCPA);
    const bool bad_vec = !isfinite3(forceA) || !isfinite3(torqueA);
    if (bad_cp || bad_vec) {
        DEME_ABORT_KERNEL(
            "In-place force reduction found invalid force/torque/contact point for A: contact owner %llu, type "
            "%u, bad_vec=%d, bad_cp=%d.\n",
            static_cast<unsigned long long>(AOwner), static_cast<unsigned int>(ContactType), static_cast<int>(bad_vec),
            static_cast<int>(bad_cp));
    }

    atomicAdd(granData->aX + AOwner, forceA.x / AOwnerMass);
    atomicAdd(granData->aY + AOwner, forceA.y / AOwnerMass);
    atomicAdd(granData->aZ + AOwner, forceA.z / AOwnerMass);

    // Direct torque contributes angular acceleration independently of the contact force
    float3 myF = forceA;
    // F is in global frame, but it needs to be in local to coordinate with moi and cntPnt
    applyOriQToVector3(myF, make_float4(-AOriQ.x, -AOriQ.y, -AOriQ.z, AOriQ.w));
    applyOriQToVector3(torqueA, make_float4(-AOriQ.x, -AOriQ.y, -AOriQ.z, AOriQ.w));
    const float3 angAcc = (cross(locCPA, myF) + torqueA) / AOwnerMOI;
    atomicAdd(granData->alphaX + AOwner, angAcc.x);
    atomicAdd(granData->alphaY + AOwner, angAcc.y);
    atomicAdd(granData->alphaZ + AOwner, angAcc.z);
}

// Take care of B
{
    const bool bad_cp = !isfinite3(locCPB);
    const bool bad_vec = !isfinite3(forceB) || !isfinite3(torqueB);
    if (bad_cp || bad_vec) {
        DEME_ABORT_KERNEL(
            "In-place force reduction found invalid force/torque/contact point for B: contact owner %llu, type "
            "%u, bad_vec=%d, bad_cp=%d.\n",
            static_cast<unsigned long long>(BOwner), static_cast<unsigned int>(ContactType), static_cast<int>(bad_vec),
            static_cast<int>(bad_cp));
    }

    atomicAdd(granData->aX + BOwner, forceB.x / BOwnerMass);
    atomicAdd(granData->aY + BOwner, forceB.y / BOwnerMass);
    atomicAdd(granData->aZ + BOwner, forceB.z / BOwnerMass);

    // Direct torque contributes angular acceleration independently of the contact force
    float3 myF = forceB;
    // F is in global frame, but it needs to be in local to coordinate with moi and cntPnt
    applyOriQToVector3(myF, make_float4(-BOriQ.x, -BOriQ.y, -BOriQ.z, BOriQ.w));
    applyOriQToVector3(torqueB, make_float4(-BOriQ.x, -BOriQ.y, -BOriQ.z, BOriQ.w));
    const float3 angAcc = (cross(locCPB, myF) + torqueB) / BOwnerMOI;
    atomicAdd(granData->alphaX + BOwner, angAcc.x);
    atomicAdd(granData->alphaY + BOwner, angAcc.y);
    atomicAdd(granData->alphaZ + BOwner, angAcc.z);
}
