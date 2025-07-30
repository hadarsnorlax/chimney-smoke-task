import math
import copy

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
    - buoyancy [m/sec^2] (double)
    - buoyancy_const [sec] (double)
    - wind_vel [m/sec] (doubles list - (x,y))
    - drag (double) 
    - aperture [m] (double)
'''
def smoke(current_time, prev_time, particles, params):
    particles_updated = copy.deepcopy(particles)
    delta_t = current_time - prev_time

    _update_particles(particles_updated, delta_t, params)
    _add_particles(particles_updated, params)

    return particles_updated


def _update_particles(particles, delta_t, params):
    for particle in particles:
        _update_acceleration(particle, delta_t, params)
        initial_velocity = (particle["X"][1], particle["Y"][1], particle["Z"][1])
        _update_velocity(particle, delta_t)
        _update_location(particle, initial_velocity, delta_t)


def _add_particles(particles, params):
    pass


def _get_gravity_acceleration():
    # Option to use scipy.constants.g
    g = -9.81 # in [m/sec^2], z^ direction
    return g


# in [m/sec^2], (x^,y^) direction
def _get_drag_acceleration(particle, params):
    # v diff between particle and wind
    v_relative = (particle["X"][1] - params["wind_vel"][0], particle["Y"][1] - params["wind_vel"][1])
    # drag_a = -drag * v_relative
    drag_a = (-params["drag"] * v_relative[0], -params["drag"] * v_relative[1])
    return drag_a


# use time_alive for new particles
def _get_buoyancy_acceleration(particle, params, time_alive=None):
    buoyancy_a = None # in [m/sec^2], z^ direction
    if (time_alive):
        # buoyancy_a = buoyancy * e^(-t_alive / buoyancy_const)
        buoyancy_a = params["buoyancy"] * math.exp(-time_alive / params["buoyancy_const"])
        return buoyancy_a
    
    # estimated time_alive for old particles by their location
    # further in xy plane is colder than higher in z, estimate time by d/v
    # buoyancy_a = buoyancy * e^(-sqrt(x^2 + y^2) / τ) * e^(-0.5 * z / τ), τ = buoyancy_const * v
    v_xy_plane = math.sqrt(particle["X"][1]**2 + particle["Y"][1]**2)
    tau_xy_plane = params["buoyancy_const"] * v_xy_plane
    tau_z = params["buoyancy_const"] * particle["Z"][1]
    d_xy_plane = math.sqrt(particle["X"][0]**2 + particle["Y"][0]**2)
    buoyancy_a = params["buoyancy"] * math.exp(-(d_xy_plane / tau_xy_plane) - (particle["Z"][0] / tau_z))
    return buoyancy_a


def _update_acceleration(particle, delta_t, params):
    g = _get_gravity_acceleration()
    drag_a = _get_drag_acceleration(particle, params)
    buoyancy_a = _get_buoyancy_acceleration(particle, params)

    particle["Z"][2] += g
    particle["Z"][2] += buoyancy_a
    particle["X"][2] += drag_a[0]
    particle["Y"][2] += drag_a[1]


def _update_velocity(particle, delta_t):
    # V(t) = V0 + a*t
    particle["X"][1] += particle["X"][2] * delta_t
    particle["Y"][1] += particle["Y"][2] * delta_t
    particle["Z"][1] += particle["Z"][2] * delta_t


# initial_velocity (x,y,z) [m/sec]
def _update_location(particle, initial_v, delta_t):
    # X(t) = X0 + V0*t + 0.5*a*t^2
    particle["X"][0] += initial_v[0] * delta_t + 0.5 * particle["X"][2] * delta_t ** 2
    particle["Y"][0] += initial_v[0] * delta_t + 0.5 * particle["Y"][2] * delta_t ** 2
    particle["Z"][0] += initial_v[0] * delta_t + 0.5 * particle["Z"][2] * delta_t ** 2
