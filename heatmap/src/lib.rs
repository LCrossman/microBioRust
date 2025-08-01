//! # A Heatmap in Rust web assembly calling d3.js
//!
//! You will need to use wasm-pack to build instead of cargo
//!  wasm-pack build --target web
//! And some way of serving locally
//!  http-server .
//! It requires the index.html in the static directory
//! Currently working with fixed data
//! and a rusty colour theme
#![allow(non_snake_case)]
pub mod canvas;
pub mod heatmap_data;

// internal imports
use canvas::drawing::draw_responsive_heatmap;
use heatmap_data::HeatmapData;

// external imports
use std::rc::Rc;
use wasm_bindgen::prelude::*;
use wasm_bindgen::JsValue;
use web_sys::console;
use web_sys::{window, CanvasRenderingContext2d, HtmlCanvasElement};

//returns a JsValue to javascript
#[wasm_bindgen(start)]
pub fn start() -> Result<(), JsValue> {
    // Get the window and document
    console::log_1(&JsValue::from_str(&format!("literal start")));
    let window = window().ok_or(JsValue::from_str("should have a window in this context"))?;
    let window = Rc::new(window);
    let window_clone = Rc::clone(&window);
    let document = window.document().ok_or(JsValue::from_str("no document"))?;
    console::log_1(&JsValue::from_str(&format!(
        "up in the start of the function"
    )));
    // Get the canvas element
    let canvas = document
        .get_element_by_id("heatmap")
        .ok_or(JsValue::from_str("Canvas element not found"))?
        .dyn_into::<HtmlCanvasElement>()?;
    console::log_1(&JsValue::from_str(&format!("called the canvas")));
    let heatmap_values = vec![
        vec![2, 1, 0, 1, 0], // row 1
        vec![1, 2, 0, 0, 1], // row 2
        vec![2, 0, 1, 2, 1], // row 3
        vec![0, 0, 0, 2, 0], // row 4
        vec![1, 2, 0, 1, 1], // row 5
    ];
    console::log_1(&JsValue::from_str(&format!("called the heatmap vals")));
    let x_labels: Vec<String> = vec!["A", "B", "C", "D", "E"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let y_labels: Vec<String> = vec!["R1", "R2", "R3", "R4", "R5"]
        .iter()
        .map(|s| s.to_string())
        .collect();

    let num_rows = heatmap_values.len(); // Should be 5
    let num_cols = heatmap_values[0].len(); // Should be 5
    let mut heatmap_data = HeatmapData::new();
    heatmap_data.values = heatmap_values.clone();
    heatmap_data.x_labels = x_labels.clone();
    heatmap_data.y_labels = y_labels.clone();
    let box_size = 100.0;
    let device_pixel_ratio = window.device_pixel_ratio();
    console::log_1(&JsValue::from_str(&format!(
        "num rows are {:?} num cols are {:?}",
        &num_rows, &num_cols
    )));

    // Dynamically set canvas size based on number of rows and columns
    let canvas_width = num_cols as f64 * box_size; // 6 columns * 50px
    let canvas_height = num_rows as f64 * box_size; // 6 rows * 50px
    canvas.set_width(canvas_width as u32);
    canvas.set_height(canvas_height as u32);
    console::log_1(&JsValue::from_str(&format!(
        "Canvas width: {}, height: {}",
        canvas.width(),
        canvas.height()
    )));

    let context = canvas
        .get_context("2d")?
        .ok_or(JsValue::from_str("Context not found"))?
        .dyn_into::<CanvasRenderingContext2d>()?;

    // Define the heatmap matrix (3x3) with values representing different colors
    context.scale(device_pixel_ratio, device_pixel_ratio)?;

    draw_responsive_heatmap(
        &context,
        heatmap_values.clone(),
        x_labels.clone(),
        y_labels.clone(),
        canvas_width,
        canvas_height,
        device_pixel_ratio,
    )?;

    let handle_heatmap_resize = move || -> Result<(), JsValue> {
        let new_width = window_clone
            .inner_width()
            .map_err(|_| JsValue::from_str("error getting inner width"))?
            .as_f64()
            .ok_or(JsValue::from_str("error converting width to f64"))?;

        let new_height = window_clone
            .inner_height()
            .map_err(|_| JsValue::from_str("error getting inner height"))?
            .as_f64()
            .ok_or(JsValue::from_str("error converting height to f64"))?;

        let canvas_new_width = (num_cols as f64 * box_size).min(new_width);
        let canvas_new_height = (num_rows as f64 * box_size).min(new_height);

        canvas.set_width(canvas_new_width as u32);
        canvas.set_height(canvas_new_height as u32);

        context
            .set_transform(1.0, 0.0, 0.0, 1.0, 0.0, 0.0)
            .map_err(|_| JsValue::from_str("error setting transform"))?;
        context
            .scale(device_pixel_ratio, device_pixel_ratio)
            .map_err(|_| JsValue::from_str("error scaling context"))?;

        draw_responsive_heatmap(
            &context,
            heatmap_values.clone(),
            x_labels.clone(),
            y_labels.clone(),
            canvas_new_width,
            canvas_new_height,
            device_pixel_ratio,
        )?;
        Ok(())
    };

    // Wrap the closure_func to handle errors
    let error_handled_heatmap_resize = move || {
        if let Err(e) = handle_heatmap_resize() {
            console::error_1(&e);
        }
    };

    let closure = Closure::wrap(Box::new(error_handled_heatmap_resize) as Box<dyn FnMut()>);

    window.add_event_listener_with_callback("resize", closure.as_ref().unchecked_ref())?;
    closure.forget();

    Ok(())
}
