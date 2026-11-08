
        void __device__ XdirToDxDyDz(const double3 &Vxdir,
                                    const double3 &Vsingular,
                                    double3 &Vx,
                                    double3 &Vy,
                                    double3 &Vz) {
            double zlen;
            if (length(Vxdir) < 1e-8) {
                Vx = make_double3(1, 0, 0);
            } else {
                Vx = Vxdir / length(Vxdir);
            }
            Vz = cross(Vx, Vsingular);
            zlen = length(Vz);

            if (zlen < 1e-6) {
                double3 mVsingualr;
                if (abs(Vsingular.z) < 0.9)
                    mVsingualr = make_double3(0, 0, 1);
                if (abs(Vsingular.y) < 0.9)
                    mVsingualr = make_double3(0, 1, 0);
                if (abs(Vsingular.x) < 0.9)
                    mVsingualr = make_double3(1, 0, 0);
                Vz = cross(Vx, mVsingualr);
                zlen = length(Vz);
            }

            Vz = Vz / length(Vz);
            Vy = cross(Vz, Vx);
        };
        