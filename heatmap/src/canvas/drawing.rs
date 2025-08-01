use wasm_bindgen::JsValue;
use web_sys::{console, CanvasRenderingContext2d};

pub fn draw_responsive_heatmap(
    context: &CanvasRenderingContext2d,
    values: Vec<Vec<i32>>,
    x_labels: Vec<String>,
    y_labels: Vec<String>,
    canvas_width: f64,
    canvas_height: f64,
    device_pixel_ratio: f64,
) -> Result<(), JsValue> {
    let rows = values.len();
    let cols = values[0].len();
    console::log_1(&JsValue::from_str(&format!("up in the draw function")));
    // Get canvas dimensions
    // Calculate dynamic padding and box size
    let adj_canvas_width = canvas_width * device_pixel_ratio;
    let adj_canvas_height = canvas_height * device_pixel_ratio;
    let padding_left = adj_canvas_width * 0.05;
    let padding_top = adj_canvas_height * 0.05;
    let padding_bottom = adj_canvas_height * 0.05;
    let _padding_right = adj_canvas_width * 0.05;

    //  let box_width = (adj_canvas_width - padding_left - padding_right) / (cols as f64 * 1.1);
    //  let box_height = (adj_canvas_height - padding_top - padding_bottom) / (rows as f64 * 1.1);

    let box_width = 30.0;
    let box_height = 30.0;
    // Clear the canvas
    console::log_1(&JsValue::from_str(&format!(
        "pad left {} pad bottom {}",
        &padding_left, &padding_bottom
    )));
    context.clear_rect(0.0, 0.0, adj_canvas_width, adj_canvas_height);
    println!("cleared rec");
    // Draw the heatmap
    for row in 0..rows {
        for col in 0..cols {
            let value = values[row][col];

            // Set color based on value
            let color = match value {
                0 => "#fee0d2",
                1 => "#fc9272",
                2 => "#de2d26",
                _ => "#FFFFFF",
            };
            //context.set_fill_style_str(&JsValue::from(color));
            context.set_fill_style_str(color);

            let x = padding_left + (col as f64 * box_width);
            let y = padding_top + (row as f64 * box_height);
            context.fill_rect(x, y, box_width, box_height);

            // Draw box borders
            //context.set_stroke_style(&JsValue::from("#FFFFFF"));
            context.set_stroke_style_str("#FFFFFF");
            context.set_line_width(2.0 / device_pixel_ratio);

            if row < rows - 1 {
                context.begin_path();
                context.move_to(x, y + box_height);
                context.line_to(x + box_width, y + box_height);
                context.stroke();
            }

            if col < cols - 1 {
                context.begin_path();
                context.move_to(x + box_width, y);
                context.line_to(x + box_width, y + box_height);
                context.stroke();
            }
        }
    }
    console::log_1(&JsValue::from_str(&format!(
        "after the rows and cols padding bottom: {}, height: {}",
        &padding_bottom,
        &(box_height * rows as f64),
    )));

    // Draw X-axis
    context.begin_path();
    //context.set_stroke_style_str(&JsValue::from("#000000"));
    context.set_stroke_style_str("#000000");
    context.move_to(padding_left, (box_height * rows as f64) + padding_bottom);
    context.line_to(
        (box_height * rows as f64) + padding_bottom,
        (box_height * rows as f64) + padding_left,
    );
    context.stroke();

    // Draw Y-axis
    context.begin_path();
    context.move_to(padding_left, padding_top);
    context.line_to(padding_left, (box_height * rows as f64) + padding_bottom);
    context.stroke();

    // Draw X-axis ticks and labels
    let label_font_size = (box_height * 0.3).min(box_width * 0.3).max(12.0);
    context.set_font(&format!("{}px Arial", label_font_size));
    context.set_text_align("center");
    context.set_text_baseline("top");

    for col in 0..cols {
        let x = padding_left + col as f64 * box_width + box_width / 2.0;
        let y = (box_height * rows as f64) + padding_bottom + 5.0; // Position below the heatmap
        context
            .fill_text(&x_labels[col], x, y)
            .map_err(|_| JsValue::from_str(&format!("Failed to draw text at column {}", col)))?;

        // Draw ticks
        context.begin_path();
        context.move_to(x, (box_height * rows as f64) + padding_bottom);
        context.line_to(x, (box_height * rows as f64) + padding_bottom + 5.0);
        context.stroke();
    }

    // Draw Y-axis ticks and labels
    context.set_text_align("right");
    context.set_text_baseline("middle");

    for row in 0..rows {
        let x = padding_left - 10.0; // Position to the left of the heatmap
        let y = padding_top + row as f64 * box_height + box_height / 2.0;
        context
            .fill_text(&y_labels[row], x, y)
            .map_err(|_| JsValue::from_str(&format!("Failed to draw text at row {}", row)))?;

        // Draw ticks
        context.begin_path();
        context.move_to(padding_left, y);
        context.line_to(padding_left - 5.0, y);
        context.stroke();
    }
    console::log_1(&JsValue::from_str(&format!(
        "at the end of draw funct Canvas width: {}, height: {}",
        &adj_canvas_width, &adj_canvas_height
    )));
    Ok(())
}
