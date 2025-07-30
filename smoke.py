import math
import copy

import numpy as np

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

    for idx, particle in enumerate(particles_updated):
        _update_particle(particle, delta_t, False, params)

    new_particles_amount = int(delta_t * params["particles_per_sec"])
    for idx in range(len(particles), len(particles) + new_particles_amount):
        particle = _create_particle(current_time, idx, params)
        _update_particle(particle, delta_t, True, params)
        particles_updated.append(particle)

    return particles_updated


def _update_particle(particle, delta_t, is_new_particle, params):
    initial_velocity = (particle["X"][1], particle["Y"][1], particle["Z"][1])
    _update_acceleration(particle, delta_t, is_new_particle, params)
    _update_velocity(particle, delta_t)
    _update_location(particle, initial_velocity, delta_t)


def _update_acceleration(particle, delta_t, is_new_particle, params):
    g = _get_gravity_acceleration()
    drag_a = _get_drag_acceleration(particle, params)
    buoyancy_a = _get_buoyancy_acceleration(particle, delta_t, is_new_particle, params)

    particle["Z"][2] = g + buoyancy_a
    particle["X"][2] = drag_a[0]
    particle["Y"][2] = drag_a[1]


def _get_gravity_acceleration():
    # Option to use scipy.constants.g
    g = -9.81 # in [m/sec^2], z^ direction
    return g


# drag acceleration from wind in [m/sec^2], (x^,y^) direction
def _get_drag_acceleration(particle, params):
    # v diff between particle and wind
    v_relative = (particle["X"][1] - params["wind_vel"][0], particle["Y"][1] - params["wind_vel"][1])
    # drag_a = -drag * v_relative
    drag_a = (-params["drag"] * v_relative[0], -params["drag"] * v_relative[1])
    return drag_a


def _get_buoyancy_acceleration(particle, delta_t, is_new_particle, params):
    if is_new_particle:
        return params["buoyancy"]
    
    g = _get_gravity_acceleration()
    prev_buoyancy = particle["Z"][2] - g
    # buoyancy_a = buoyancy * e^(-t / buoyancy_const)
    decayed_buoyancy = prev_buoyancy * math.exp(-delta_t / params["buoyancy_const"])

    return decayed_buoyancy


def _update_velocity(particle, delta_t):
    axes = ["X", "Y", "Z"]
    for axis in axes:
        # V(t) = V0 + a*t
        particle[axis][1] += particle[axis][2] * delta_t


# initial_velocity (x,y,z) [m/sec]
def _update_location(particle, initial_v, delta_t):
    axes = ["X", "Y", "Z"]
    for idx, axis in enumerate(axes):
        # X(t) = X0 + V0*t + 0.5*a*t^2
        particle[axis][0] += initial_v[idx] * delta_t + 0.5 * particle[axis][2] * delta_t ** 2


def _create_particle(current_time, particle_idx, params):
    rng = _get_rng_particle_seed(current_time, particle_idx)
    initial_loc = _get_random_chimney_loc(rng, params)
    initial_v = _get_random_initial_velocity(rng, params)

    particle = {
        "X": [initial_loc[0], initial_v[0], 0],
        "Y": [initial_loc[1], initial_v[1], 0],
        "Z": [initial_loc[2], initial_v[2], 0]
    }

    return particle


# create unique seed to each particle, for consistent results with random authentic look
def _get_rng_particle_seed(current_time, particle_number):
    seed = int(current_time * 1000) + particle_number
    rng = np.random.default_rng(seed)
    return rng


# use rng_seed of specific particle from _get_rng_particle_seed
def _get_random_chimney_loc(rng_seed, params):
    chimney_radius = params["aperture"] / 2
    radius = chimney_radius * rng_seed.random()
    theta = 2 * np.pi * rng_seed.random()
    x = radius * np.cos(theta)
    y = radius * np.sin(theta)
    return [x, y, 0]


# use rng_seed of specific particle from _get_rng_particle_seed
def _get_random_initial_velocity(rng_seed, params):
    # random gaussian distribution
    v = rng_seed.normal(loc=params["ini_vel"], scale=params["vel_std"])
    return v