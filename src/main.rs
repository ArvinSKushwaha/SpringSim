use std::{
    collections::VecDeque,
    f64::consts::PI,
    ops::{Add, AddAssign, Mul},
    time::Instant,
};

use eframe::{
    egui::{
        self,
        plot::{Line, Plot, PlotPoint, PlotUi, Polygon, Text},
        RichText,
    },
    emath::Align2,
};
use glam::{Vec2, vec2};

#[derive(Default)]
struct System {
    masses: Vec<f32>,
    positions: Vec<Vec2>,
    velocities: Vec<Vec2>,
    accelerations: Vec<Vec2>,
    interactions: Vec<Force>,
    time: f32,
}

const GRAVITATIONAL_ACCELERATION: f32 = 9.81;

#[derive(Copy, Clone)]
enum Force {
    Spring {
        m1: usize,
        m2: usize,
        k: f32,
        l: f32,
    },
    Gravity {
        m: usize,
    },
}

impl System {
    fn compute_forces(&mut self) {
        for acc in self.accelerations.iter_mut() {
            *acc *= 0.0;
        }

        for force in &self.interactions {
            match *force {
                Force::Spring { m1, m2, k, l } => {
                    let displacement = self.positions[m2] - self.positions[m1];
                    let magnitude = displacement.length();
                    self.accelerations[m1] +=
                        displacement * (magnitude - l) / (self.masses[m1] * magnitude);
                    self.accelerations[m2] +=
                        -displacement * (magnitude - l) / (self.masses[m2] * magnitude);
                }
                Force::Gravity { m } => {
                    self.accelerations[m] += GRAVITATIONAL_ACCELERATION * Vec2::NEG_Y;
                }
            }
        }
    }
}

fn main() {
    let mut system = System::default();
    system.masses.extend_from_slice(&[1e+7, 1.0, 1.0]);
    system.positions.extend_from_slice(&[vec2(0.0, 0.0), vec2(1.2, 0.0), vec2(2.2, 0.0)]);
    system.velocities.extend_from_slice(&[Vec2::ZERO; 3]);
    system.accelerations.extend_from_slice(&[Vec2::ZERO; 3]);
    system.interactions.push(Force::Spring { m1: 0, m2: 1, k: 100.0, l: 10.0 });
    system.interactions.push(Force::Spring { m1: 1, m2: 2, k: 200.0, l: 10.0 });
    system.interactions.push(Force::Gravity { m: 1 });
    system.interactions.push(Force::Gravity { m: 2 });

    eframe::run_native(
        "Simulation",
        eframe::NativeOptions::default(),
        Box::new(|_cc| {
            Box::new(MyApp {
                system,
                ..Default::default()
            })
        }),
    );
}

#[derive(Default)]
struct MyApp {
    system: System,
    state: SimState,
    elapsed_time: VecDeque<Instant>,
}

struct SimState {
    step_size: u32,
}

impl Default for SimState {
    fn default() -> Self {
        Self { step_size: 1000 }
    }
}

impl System {
    fn update(&mut self, dt: f32) {
        self.compute_forces();

        for ((pos, vel), acc) in self.positions.iter_mut().zip(self.velocities.iter_mut()).zip(self.accelerations.iter()) {
            *pos += *vel * dt + *acc * dt * dt / 2.0;
            *vel += *acc * dt;
        }

        self.time += dt;
    }

    fn draw(&self, plot_ui: &mut PlotUi) {
        for i in 0..self.masses.len() {
            let pos = self.positions[i];
            let (x, y) = (pos.x as f64, pos.y as f64);

            plot_ui.polygon(Polygon::new(vec![
                [x - 0.25, y - 0.25],
                [x + 0.25, y - 0.25],
                [x + 0.25, y + 0.25],
                [x - 0.25, y + 0.25],
            ]));
        }

        for interactions in &self.interactions {
            if let &Force::Spring { m1, m2, .. } = interactions {
                let start = self.positions[m1];
                let end = self.positions[m2];
                plot_ui.line(Line::new(vec![[start.x as f64, start.y as f64], [end.x as f64, end.y as f64]]));
            }
        }
    }
}

impl eframe::App for MyApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        if self.elapsed_time.len() >= 100 {
            self.elapsed_time.pop_front();
        }

        self.elapsed_time.push_back(Instant::now());

        egui::Window::new("Settings").show(ctx, |ui| {
            ui.add(
                egui::Slider::new(&mut self.state.step_size, 1..=100000)
                    .text("Steps")
                    .logarithmic(true),
            );
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            Plot::new("Sim")
                .show_axes([true, true])
                .data_aspect(1.0)
                .show(ui, |plot| {
                    for _ in 0..self.state.step_size {
                        self.system.update(0.00001);
                    }

                    self.system.draw(plot);

                    let [x, y] = plot.plot_bounds().min();
                    let w = plot.plot_bounds().width();
                    let h = plot.plot_bounds().height();
                    plot.text(
                        Text::new(
                            PlotPoint::new(x + w / 10., y + h / 10.),
                            RichText::new(format!("Time: {:.4}", self.system.time)).size(30.),
                        )
                        .anchor(Align2::LEFT_BOTTOM),
                    );
                });
        });

        ctx.request_repaint();

        let fps = 1.0 / self.elapsed_time[0].elapsed().as_secs_f64() * self.elapsed_time.len() as f64;
    }
}
