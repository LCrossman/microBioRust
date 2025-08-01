//! This module contains the data structure for the heatmap

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone, Debug)]
pub struct HeatmapData {
    pub values: Vec<Vec<i32>>,
    pub x_labels: Vec<String>,
    pub y_labels: Vec<String>,
}

impl HeatmapData {
    // Constructor method
    pub fn new() -> Self {
        HeatmapData {
            values: vec![vec![0]],
            x_labels: Vec::new(),
            y_labels: Vec::new(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let heatmap_data = HeatmapData::new();
        assert_eq!(heatmap_data.values, vec![vec![0]]);
        assert_eq!(heatmap_data.x_labels, Vec::<String>::new());
        assert_eq!(heatmap_data.y_labels, Vec::<String>::new());
    }
}
