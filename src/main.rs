use rand::{distributions::Uniform, Rng};
use std::fs::File;
use std::io::Write;

fn main() {
    const PARTICLE_COUNT: usize = 5;
    const TIME_DELTA: f64 = 0.01;
    const TIME_STEPS: usize = 100;

    let mut pos_file = File::create("pos_out").expect("Unable to create file!");
    let mut vel_file = File::create("vel_out").expect("Unable to create file!");
    let mut settings_file = File::create("settings_out").expect("Unable to create file!");

    let mut rng = rand::thread_rng();

    // probability distributions
    let locality_range = Uniform::new(-1.0, 1.0);
    let mass_range = Uniform::new(0.0, 1.0);

    // current positions as uniform random distribution on 3D cube
    let mut pos: Vec<Vec<f64>> = Vec::new();
    for _id in 0..PARTICLE_COUNT {
        pos.push((0..3).map(|_| rng.sample(&locality_range)).collect())
    };

    // current velocities as uniform random distribution on 3D cube
    let mut vel: Vec<Vec<f64>> = Vec::new();
    for _id in 0..PARTICLE_COUNT {
        vel.push((0..3).map(|_| rng.sample(&locality_range)).collect())
    };

    // buffer for current force
    let mut fpm: Vec<Vec<f64>> = vec![vec![0.0; 3]; PARTICLE_COUNT];

    let masses: Vec<f64> = (0..PARTICLE_COUNT).map(|_| rng.sample(&mass_range)).collect();

    for _timestep in 0..TIME_STEPS {
        
        // compute force on each particle
        for id in 0..PARTICLE_COUNT {
            let b = Particle::new(pos[id][0], pos[id][1], pos[id][2], masses[id]);

            for id_other in 0..PARTICLE_COUNT {
                if id != id_other {
                    let a = Particle::new(pos[id_other][0], pos[id_other][1], pos[id_other][2], masses[id_other]);
                    // add contribution due to single particle a to total force per mass fpm on particle b
                    fpm[id] = fpm[id].iter().zip(force_per_mass(&a, &b).iter()).map(|(&f, &fid)| f + fid).collect();
                }
            }
        };
        // update position and velocity of each particle using Euler forward
        for id in 0..PARTICLE_COUNT {
            vel[id] = vel[id].iter().zip(fpm[id].iter()).map(|(&v, &dv)| v + TIME_DELTA * dv).collect();
            pos[id] = pos[id].iter().zip(vel[id].iter()).map(|(&x, &dx)| x + TIME_DELTA * dx).collect();

        };
        for id in 0..PARTICLE_COUNT {
            write!(pos_file, "{} {} {}\n", pos[id][0], pos[id][1], pos[id][2]).expect("Unable to write to file!");
            write!(vel_file, "{} {} {}\n", vel[id][0], vel[id][1], vel[id][2]).expect("Unable to write to file!");
        }
        fpm = vec![vec![0.0; 3]; PARTICLE_COUNT];
    };
    write!(settings_file, "PARTICLE_COUNT TIME_DELTA TIME_STEPS\n").expect("Unable to write to file!");
    write!(settings_file, "{PARTICLE_COUNT} {TIME_DELTA} {TIME_STEPS}\n").expect("Unable to write to file!");
}

struct Particle {
    x: f64,
    y: f64,
    z: f64,
    m: f64,
}

impl Particle {
    fn new(x: f64, y: f64, z: f64, m: f64) -> Particle {
        Particle { x, y, z, m }
    }
}
// force per mass due to particle a on particle b, where mass = b.m
fn force_per_mass(a: &Particle, b: &Particle) -> Vec<f64> {
    const GRAVITATIONAL_CONST: f64 = 1.0;

    let delta = vec![a.x - b.x, a.y - b.y, a.z - b.z];
    let norm = (&delta[0].powi(2) + &delta[1].powi(2) + &delta[2].powi(2)).sqrt();
    let scalar = GRAVITATIONAL_CONST * a.m / norm.powi(3);

    delta.iter().map(|&x| &scalar * x).collect()
}