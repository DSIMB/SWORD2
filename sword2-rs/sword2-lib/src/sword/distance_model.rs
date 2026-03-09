//! Distance-to-model quality scoring.
//!
//! Port of `DistanceModel.pm`. Computes a signed distance from a decision
//! boundary in normalised (CR, CPD) space that separates good from bad
//! domain decompositions.

/// Compute the signed distance of a point (max_cr, mean_density) to the linear
/// decision boundary.
///
/// Returns a positive value when the point is inside the "good" region and a
/// negative value when it is outside.
///
/// If `_abs` is non-zero the absolute value semantics from the original Perl
/// code are preserved (the flag is present in the API for compatibility but
/// the computation is identical — the original code never actually changed
/// behaviour based on it).
pub fn distance_model(mx: f64, my: f64, _abs: i32) -> f64 {
    // Normalisation bounds (from training data)
    let min_x: f64 = 0.008501;
    let max_x: f64 = 1.099581;
    let min_y: f64 = 1.797981;
    let max_y: f64 = 4.124218;

    // Min-max normalisation
    let mx = (mx - min_x) / (max_x - min_x);
    let my = (my - min_y) / (max_y - min_y);

    // Model parameters
    let horizontal: f64 = 0.5884363;
    let vertical: f64 = 0.2046999;
    let diagonal_intercept: f64 = 0.4474396;
    let diagonal_slope: f64 = 1.680319;
    let diagonal_vertical: f64 = 0.2075470;
    let diagonal_horizontal: f64 = 0.7867766;

    // Distances to the three boundaries
    let distance_to_diag = (diagonal_slope * mx - my + diagonal_intercept).abs()
        / (1.0 + diagonal_slope * diagonal_slope).sqrt();
    let distance_to_horizontal = (my - horizontal).abs();
    let distance_to_vertical = (mx - vertical).abs();

    // Which side of the diagonal line is the point on?
    // Line through A=(0, diagonal_intercept) and B=(0.3288426, 1)
    let ax: f64 = 0.0;
    let ay: f64 = diagonal_intercept;
    let bx: f64 = 0.3288426;
    let by: f64 = 1.0;
    let d = (bx - ax) * (my - ay) - (by - ay) * (mx - ax);

    // Zone classification
    if my > horizontal && mx < vertical && d < 0.0 {
        // ZONE 1
        let dist = distance_to_diag
            .min(distance_to_horizontal)
            .min(distance_to_vertical);
        return dist * -1.0;
    }

    if (my < diagonal_horizontal && my > horizontal && d < 0.0)
        || (mx > diagonal_vertical && mx < vertical && d < 0.0)
        || (my < diagonal_horizontal && mx > diagonal_vertical && d < 0.0)
    {
        // ZONE 7, 8, 9
        return distance_to_diag * -1.0;
    }

    if my < horizontal {
        // ZONE 6
        return distance_to_horizontal * -1.0;
    }

    if mx > vertical {
        // ZONE 5
        return distance_to_vertical * -1.0;
    }

    // ZONE 2, 3, or 4 — inside
    let dist = distance_to_diag
        .min(distance_to_horizontal)
        .min(distance_to_vertical);
    dist // positive ⇒ inside
}

/// Map a real-valued distance into a 1–5 star rating.
///
/// Port of `step_function()` from SWORD Perl script.
pub fn step_function(input: f64) -> usize {
    if input >= 0.15 {
        5
    } else if input >= 0.05 {
        4
    } else if input >= -0.05 {
        3
    } else if input >= -0.15 {
        2
    } else {
        1
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_step_function() {
        assert_eq!(step_function(0.20), 5);
        assert_eq!(step_function(0.10), 4);
        assert_eq!(step_function(0.00), 3);
        assert_eq!(step_function(-0.10), 2);
        assert_eq!(step_function(-0.20), 1);
    }

    #[test]
    fn test_distance_model_returns_value() {
        // Verify the function returns a finite value for typical inputs
        let d = distance_model(0.05, 3.0, 0);
        assert!(d.is_finite(), "Expected finite distance, got {}", d);
    }
}
