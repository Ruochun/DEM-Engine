// Double-precision local DFC calculations. DEME inputs/history and final force retain solver storage types.
if (overlapDepth > 0) {
    // Get material properties
    //
    double mortar_layer = mortar_layer_mat[bodyAMatType];
    double ENm = ENm_mat[bodyAMatType];
    double ENa = ENa_mat[bodyAMatType];
    double alpha = alpha_mat[bodyAMatType];
    double beta = beta_mat[bodyAMatType];
    double np = np_mat[bodyAMatType];
    double sgmTmax = sgmTmax_mat[bodyAMatType];
    double sgmTau0 = sgmTau0_mat[bodyAMatType];
    double kappa0 = kappa0_mat[bodyAMatType];
    double eta_inf = eta_inf_mat[bodyAMatType];
    double flocbeta = flocbeta_mat[bodyAMatType];
    double flocm = flocm_mat[bodyAMatType];
    double flocTcr = flocTcr_mat[bodyAMatType];

    //
    double h = mortar_layer;
    double R1 = 5000.0;
    double R2 = 5000.0;
    double Rmin, Rmax, delta_rmin;
    double Lij, L0ij;
    // Material parameters of viscous part which should be defined inside
    // contactinfo
    double eta0 = kappa0 * eta_inf;
    double Deps0 = sgmTau0 / eta0;

    //
    // Get state variables relevant to this contact
    //
    double3 statevar;
    statevar.x = contact_info_strain_x;
    statevar.y = contact_info_strain_y;
    statevar.z = contact_info_strain_z;
    double lambda0 = contact_info_lambda;
    //
    // ==== SIGMA-t influenced by Flocculation ======
    sgmTmax = sgmTmax * (1. + lambda0);
    //
    if (lambda0 == 0.) {
    lambda0 = lambda_init_mat[bodyAMatType];
    double val0 = ParticleThixoAvg[AOwner];
    double val1 = ParticleThixoAvg[BOwner];
    //
    if (val0 <= 0.) {
        val0 = lambda0;
    }
    if (val1 <= 0.) {
        val1 = lambda0;
    }
    lambda0 = DEME_MIN(val0, val1);
    contact_info_lambda = lambda0;
    // std::cout << "lambda-4 " << map_contact_info[mykey].lambda <<  "key "
    // << mykey << std::endl;
    contact_info_step_time = time;
    }
    //
    // initialize force values
    //
    force = make_float3(0, 0, 0);

    double3 Vx, Vy, Vz;
    double3 normal_dir = make_double3(-(double)B2A.x, -(double)B2A.y, -(double)B2A.z); // DEME uses B2A as normal direction, opposite
                            // of Chrono as it seems

    XdirToDxDyDz(normal_dir, make_double3(0, 1, 0), Vx, Vy, Vz);
    // std::cout<<" Vy: "<<Vy<<" Vz: "<<Vz<<"\t";

    double delta = overlapDepth;

    R1 = ARadius;
    R2 = BRadius;

    if (ContactType == deme::SPHERE_SPHERE_CONTACT) {
        // Concrete - Concrete Interaction
        // Contact of two DFC spheres
        // center to center distance between 2 objects
        L0ij = R1 + R2 - h;
        Lij = abs(R1 + R2 - delta);
        Rmin = DEME_MIN(R1, R2);
        Rmax = DEME_MAX(R1, R2);
        delta_rmin = delta / 2 * (2.0 * Rmax - delta) / Lij;

    } else if (ContactType == deme::SPHERE_MESH_CONTACT) { 
        // Concrete-container interaction
        R2 = 5.; //  0.;
        h = mortar_layer / 2.;
        // sgmTmax=this->material->Get_sgmTmax()/2;
        L0ij = R1 + R2 - h;
        Lij = abs(R1 + R2 - delta);
        Rmin = R1;
        // Rmin=std::min(R1, R2);
        Rmax = DEME_MAX(R1, R2);
        delta_rmin = delta;

    } else { 
        /// Concrete-analytical interaction
        R2 = 5.; // Match Chrono floor/plate box half-thickness.
        h = mortar_layer / 2.;
        // sgmTmax=this->material->Get_sgmTmax()/2;
        L0ij = R1 + R2 - h;
        Lij = abs(R1 + R2 - delta);
        Rmin = R1;
        Rmax = DEME_MAX(R1, R2);
        delta_rmin = delta;
    }
    double viscous_shear_multiplier = (ContactType == deme::SPHERE_MESH_CONTACT) ? 0.0 : 1.0;

    double ai = abs(Rmin - delta_rmin);
    double radius2 = Rmin * Rmin - ai * ai;
    double contact_area = deme::PI * radius2;

    //
    // modify penetration according to initial overlap
    //
    double delta_new = delta - h;
    //
    //
    // Relative velocities at contact at midplane (p0=(p1+p2)/2)
    //
    // In DEME, contact point in global is just called contactPnt

    // In DEME, local contact point is called make_double3(locCPA.x, locCPA.y, locCPA.z) and make_double3(locCPB.x, locCPB.y, locCPB.z)
    // DEME uses its own way to calculate rel vel
    double3 rotVelCPA, rotVelCPB;
    {
    // We also need the relative velocity between A and B in global frame to
    // use in the damping terms To get that, we need contact points'
    // rotational velocity in GLOBAL frame This is local rotational velocity
    // (the portion of linear vel contributed by rotation)
    rotVelCPA = cross(make_double3(ARotVel.x, ARotVel.y, ARotVel.z), make_double3(locCPA.x, locCPA.y, locCPA.z));
    rotVelCPB = cross(make_double3(BRotVel.x, BRotVel.y, BRotVel.z), make_double3(locCPB.x, locCPB.y, locCPB.z));
    // This is mapping from local rotational velocity to global
    applyOriQToVector3<double, double>(rotVelCPA.x, rotVelCPA.y,
                                            rotVelCPA.z, AOriQ.w, AOriQ.x,
                                            AOriQ.y, AOriQ.z);
    applyOriQToVector3<double, double>(rotVelCPB.x, rotVelCPB.y,
                                            rotVelCPB.z, BOriQ.w, BOriQ.x,
                                            BOriQ.y, BOriQ.z);
    }
    double time_scale_val = 1.;
    double3 relvel = (make_double3(BLinVel.x, BLinVel.y, BLinVel.z) + rotVelCPB) - (make_double3(ALinVel.x, ALinVel.y, ALinVel.z) + rotVelCPA); // A2B vel
    relvel = relvel / time_scale_val;

    double relvel_n_mag = dot(relvel, normal_dir);
    double3 relvel_n = relvel_n_mag * normal_dir;
    double3 relvel_t = relvel - relvel_n;
    double relvel_t_mag = length(relvel_t);
    //
    // Calculate displacement increment in normal and tangential direction
    //
    double dT = ts;
    double delta_t = relvel_t_mag * dT;
    double delta_n = relvel_n_mag * dT;
    double3 v_delta_t = relvel_t * dT;
    //
    //  Calculate the strain increment in each local direction
    //
    double depsN = delta_n / Lij;
    double delta_M = dot(relvel, Vy) * dT; // v_delta_t.Dot(Vy);
    double delta_L = dot(relvel, Vz) * dT; // v_delta_t.Dot(Vz);
    double depsM = delta_M / Lij;
    double depsL = delta_L / Lij;
    //
    //
    //
    double epsA = log(1 - h / (L0ij));
    double epsN = log(Lij / L0ij);


    
    // double epsN=statevar.x+depsN;
    double epsM = statevar.y + depsM;
    double epsL = statevar.z + depsL;
    //
    double epsT = pow((epsM * epsM + epsL * epsL), 0.5);
    double epsQ = pow(epsN * epsN + alpha * epsT * epsT, 0.5);
    //
    //
    //
    statevar.x = epsN;
    statevar.y = epsM;
    statevar.z = epsL;
    contact_info_strain_x = statevar.x;
    contact_info_strain_y = statevar.y;
    contact_info_strain_z = statevar.z;
    contact_info_step_time = time;
    //
    //
    //
    if (epsN < 0) {
    ////////////////////////////////////////////////////////////
    // Compressive contact;
    ////////////////////////////////////////////////////////////
    double sgmN = 0;
    double sgmM = 0;
    double sgmL = 0;
    double sgmT = 0;
    //

    double stot;
    if (epsN >= epsA) {
        // stot=epsQ*ENm;
        sgmN = epsN * ENm;
        contact_info_strain_y = 0;
        contact_info_strain_z = 0;
    } else {
        // Deep compression branch ported from FreshConcreteContact_V10.cpp.
        // Aggregate stiffness contribution is disabled there.
        sgmN = epsN * ENm;
        sgmM = alpha * ENm * epsM;
        sgmM = ((epsM > 0) - (epsM < 0)) * DEME_MIN(abs(sgmM), abs(sgmN)) * 0.;
        sgmL = alpha * ENm * epsL;
        sgmL = ((epsL > 0) - (epsL < 0)) * DEME_MIN(abs(sgmL), abs(sgmN)) * 0.;
    }

    // std::cout<<"delta_new: "<<delta_new<<" epsN: "<<epsN<< " epsQ:
    // "<<epsQ<<" stot: "<<stot<<" sgmN: "<<sgmN<<std::endl;
    //////////////////////////////////////////////////////////
    //
    // Viscous stresses
    //
    //////////////////////////////////////////////////////////
    //
    // depsN=relvel_n_mag/Lij*dT;
    // depsT=relvel_t_mag/Lij*dT;
    // calculate strain rates
    double vdepsN = depsN / dT;
    double vdepsM = depsM / dT;
    double vdepsL = depsL / dT;
    double vDeps = pow(
        (beta * vdepsN * vdepsN + vdepsM * vdepsM + vdepsL * vdepsL), 0.5);
    //
    // std::cout << "Lij : " << Lij << " vdeps_Comp : " << vDeps <<
    // std::endl;
    double eta;
    double lambda = lambda0;

    {
        // double lambda = rungeKutta4<double>(dfloc, lambda0, current_time,
        // flocbeta, vDeps, flocm, flocTcr, dT, 1); lambda=0;
        // map_contact_info[mykey].lambda=lambda;
        if (vDeps <= Deps0) {
        eta = eta0;
        } else {
        eta = eta_inf * pow(abs(vDeps), np - 1.) +
                sgmTau0 * (1.0 + lambda) / vDeps;
        }
    }

    double sgmN_vis = beta * eta * vdepsN;
    double sgmM_vis = viscous_shear_multiplier * eta * vdepsM;
    double sgmL_vis = viscous_shear_multiplier * eta * vdepsL;
    double sgmT_vis = pow(sgmM_vis * sgmM_vis + sgmL_vis * sgmL_vis, 0.5);
    //
    //
    //
    double forceN = contact_area * (sgmN + sgmN_vis);
    double forceM = contact_area * (sgmM + sgmM_vis);
    double forceL = contact_area * (sgmL + sgmL_vis);
    // std::cout<<"epsN "<<epsN<<" sgmN "<<sgmN<<std::endl;
    // force[0]=forceN;force[1]=contact_area*(sgmM-sgmM_vis);force[2]=contact_area*(sgmL-sgmL_vis);
    const double3 force_double = forceN * normal_dir + forceM * Vy + forceL * Vz;
    force = make_float3((float)force_double.x, (float)force_double.y, (float)force_double.z);
    // std::cout<<" force "<<force<<std::endl;
    /*
    force = -forceN * normal_dir;
    if (relvel_t_mag >= sys.GetSlipVelocityThreshold())
            force -= (forceT / relvel_t_mag) * relvel_t;
*/
    // return -force;

    } else {
    ////////////////////////////////////////////////////////////
    // Tensile contact
    ////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////
    //
    // Calculate Stress from material stiffness
    //
    //////////////////////////////////////////////////////////
    double sgmN = ENm * epsN;
    double sgmM = 0;
    double sgmL = 0;
    if (sgmN > sgmTmax)
        sgmN = sgmTmax;
    //
    //////////////////////////////////////////////////////////
    //
    // Viscous stresses
    //
    //////////////////////////////////////////////////////////
    //
    // calculate strain rates
    double vdepsN = depsN / dT;
    double vdepsM = depsM / dT;
    double vdepsL = depsL / dT;
    double vDeps = pow(
        (beta * vdepsN * vdepsN + vdepsM * vdepsM + vdepsL * vdepsL), 0.5);
    //
    // std::cout << "Lij : " << Lij << " vdeps_Tension : " << vDeps <<
    // std::endl;
    double eta;
    double lambda = lambda0;

    {
        // double lambda = rungeKutta4<double>(dfloc, lambda0, current_time,
        // flocbeta, vDeps, flocm, flocTcr, dT, 1); lambda=0;
        // map_contact_info[mykey].lambda=lambda;
        if (vDeps <= Deps0) {
        eta = eta0;
        } else {
        eta = eta_inf * pow(abs(vDeps), np - 1.) +
                sgmTau0 * (1.0 + lambda) / vDeps;
        }
    }
    // std::cout<<"relvel"<<relvel<<" Deps0: "<<Deps0<<" vDeps"<<vDeps<<"
    // eta: "<<eta<<std::endl;
    double sgmN_vis = beta * eta * vdepsN;
    double sgmM_vis = viscous_shear_multiplier * eta * vdepsM;
    double sgmL_vis = viscous_shear_multiplier * eta * vdepsL;
    // double sgmT=pow(sgmM*sgmM+sgmL*sgmL,0.5);
    //////////////////////////////////////////////////////////
    //
    // Combine Viscous stresses and stiffness stresses and calculate forces
    //
    //////////////////////////////////////////////////////////
    // exit(0);
    double forceN = contact_area * (sgmN + sgmN_vis);
    double forceM = contact_area * (sgmM + sgmM_vis);
    double forceL = contact_area * (sgmL + sgmL_vis);
    // std::cout<<"epsN "<<epsN<<" sgmN "<<sgmN<<std::endl;
    // double forceT = sgmT * contact_area;
    // force[0]=forceN;force[1]=contact_area*(sgmM-sgmM_vis);force[2]=contact_area*(sgmL-sgmL_vis);
    const double3 force_double = forceN * normal_dir + forceM * Vy + forceL * Vz;
    force = make_float3((float)force_double.x, (float)force_double.y, (float)force_double.z);
    /*
    force = -forceN * normal_dir;
    if (relvel_t_mag >= sys.GetSlipVelocityThreshold())
            force -= (forceT / relvel_t_mag) * relvel_t;
    */
    // return -force;
    }
    //   ChVector3d torque = Vcross((p1 - p0), -force);
    // return std::make_pair(-force, torque);
    //   return {-force, torque};

    // if (isnan(length(force))) {
    //     printf("Force NaN!!!!\n");
    // }
} else {
    // No contact
    contact_info_step_time = 0.f;
    contact_info_lambda = 0.f;
    contact_info_strain_x = 0.f;
    contact_info_strain_y = 0.f;
    contact_info_strain_z = 0.f;
    force = make_float3(0, 0, 0);
}
