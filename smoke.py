'''
- current_time [sec] (double)
- prev_time [sec] (double)
- particles (list):
    - particle (dict):
        - X, Y, Z (lists): 
            - position [m] (double)
            - velocity [m/sec] (double)
            - acceleration [m/sec^2]) (double)
- params dict:
    - particles_per_sec (double)
    - ini_vel [m/sec] (doubles list - (x,y,z))
    - vel_std [m/sec] (doubles list - (x,y,z))
    - bouyancy [m/sec^2] (double)
    - buoyancy_const [sec] (double)
    - wind_vel [m/sec] (doubles list - (x,y))
    - drag (double) 
    - aperture [m] (double)
'''
def smoke(current_time, prev_time, particles, params):
    particles_updated = []

    _update_particles(particles, particles_updated, params)
    _add_particles(particles_updated, params)

    return particles_updated


def _update_particles(particles, particles_updated, params):
    pass


def _add_particles(particles, params):
    pass
