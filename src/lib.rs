use std::f64::consts::{PI, TAU, FRAC_PI_4};

#[macro_use]
extern crate approx;

pub const MAX_ZOOM: u64 = 29;

type D = f64;
type M = f64;
type P = u64;
type Z = u64;
type T = u64;


#[derive(Debug)]
pub struct GlobalMercator {
    tile_size: u64,
    a: f64,
    initial_resolution: f64,
}

impl GlobalMercator {
    pub fn new(tile_size: u64) -> Self {
        let a = 6378137.0;
        let initial_resolution = TAU * a / (tile_size as f64);

        Self {
            tile_size,
            initial_resolution,
            a,
        }
    }

    /// Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:3857
    pub fn lat_lon_to_meters(&self, lon: D, lat: D) -> (M, M) {
        let mx = self.a * lon.to_radians();
        let my = self.a * f64::ln(f64::tan(FRAC_PI_4 + lat.to_radians() / 2.0));
        (mx, my)
    }

    /// Converts XY point from Spherical Mercator EPSG:3857 to lat/lon in WGS84 Datum
    pub fn meters_2_lat_lon(&self, mx: M, my: M) -> (D, D) {
        let lon = (mx / self.a).to_degrees();
        let lat = (f64::atan(f64::sinh(my / self.a))).to_degrees();
        (lon, lat)
    }

    /// Converts pixel coordinates in given Zoom level of pyramid to EPSG:3857
    pub fn pixels_to_meters(&self, px: P, py: P, zoom: Z) -> (M, M) {
        let res = self.resolution(zoom);
        let mx = (px as f64) * res - PI * self.a;
        let my = (py as f64) * res - PI * self.a;
        (mx, my)
    }

    /// Converts EPSG:3857 to pyramid pixel coordinates in given Zoom level
    pub fn meters_to_pixels(&self, mx: M, my: M, zoom: Z) -> (P, P) {
        let res = self.resolution(zoom);
        let px = M::floor((mx + PI * self.a) / res) as P;
        let py = M::floor((my + PI * self.a) / res) as P;
        (px, py)
    }

    /// Returns a Tile covering region in given pixel coordinates
    pub fn pixels_to_tile(&self, px: P, py: P) -> (T, T) {
        let tx = f64::ceil(px as f64 / self.tile_size as f64) as u64 - 1;
        let ty = f64::ceil(py as f64 / self.tile_size as f64) as u64 - 1;
        (tx, ty)
    }

    /// Move the origin of pixel coordinates to top-left corner
    pub fn pixels_to_raster(&self, px: P, py: P, zoom: Z) -> (P, P) {
        let map_size = (self.tile_size) << zoom;
        (px, map_size - py)
    }

    /// Returns tile for given mercator coordinates
    pub fn meters_to_tile(&self, mx: M, my: M, zoom: Z) -> (T, T) {
        let (px, py) = self.meters_to_pixels(mx, my, zoom);
        self.pixels_to_tile(px, py)
    }

    /// Returns bounds of the given tile in EPSG:3857 coordinates
    pub fn tile_bounds(&self, tx: T, ty: T, zoom: Z) -> (M, M, M, M) {
        let (minx, miny) = self.pixels_to_meters(tx * self.tile_size, ty * self.tile_size, zoom);
        let (maxx, maxy) = self.pixels_to_meters(
            (tx + 1) * self.tile_size,
            (ty + 1) * self.tile_size,
            zoom,
        );
        (minx, miny, maxx, maxy)
    }

    /// Returns bounds of the given tile in latitude/longitude using WGS84 datum
    pub fn tile_lat_lon_bounds(&self, tx: T, ty: T, zoom: Z) -> (D, D, D, D) {
        let (min_x, min_y, max_x, max_y) = self.tile_bounds(tx, ty, zoom);
        let (min_lon, min_lat) = self.meters_2_lat_lon(min_x, min_y);
        let (max_lon, max_lat) = self.meters_2_lat_lon(max_x, max_y);
        (min_lon, min_lat, max_lon, max_lat)
    }

    /// Resolution (meters/pixel) for given zoom level (measured at Equator)
    pub fn resolution(&self, zoom: Z) -> f64 {
        self.initial_resolution / (2.0f64.powi(zoom as i32))
    }

    /// Maximal scaledown zoom of the pyramid closest to the pixelSize.
    pub fn zoom_for_pixel_size(&self, pixel_size: f64) -> Z {
        for i in 0..=MAX_ZOOM {
            if pixel_size > self.resolution(i) {
                return std::cmp::max(0, i - 1);
            }
        }
        panic!("Invalid pixel_size: {}", pixel_size);
    }

    /// Converts TMS tile coordinates to Google Tile coordinates and vice versa
    pub fn google_tile(&self, tx: T, ty: T, zoom: Z) -> (T, T) {
        (tx, (2u64.pow(zoom as u32) - 1) - ty)
    }

    /// Converts TMS tile coordinates to Microsoft quad_tree
    pub fn quad_tree(&self, tx: T, ty: T, zoom: Z) -> String {
        let mut quad_key = String::new();
        let ty = (f64::powi(2.0, zoom as i32) - 1.0) as T - ty;
        for i in (1..(zoom + 1) as i32).rev() {
            let mut digit = 0;
            let mask = 1 << (i - 1);
            if (tx & mask) != 0 {
                digit += 1;
            }
            if (ty & mask) != 0 {
                digit += 2;
            }
            quad_key.push_str(format!("{}", digit).as_str());
        }
        quad_key
    }
}

// Add tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let transformer = GlobalMercator::new(256);
        assert_eq!(transformer.tile_size, 256);
    }

    #[test]
    fn lat_lon_to_meters() {
        let transformer = GlobalMercator::new(256);
        let (x, y) = transformer.lat_lon_to_meters(0.0, 0.0);
        assert_relative_eq!(x, 0.0, epsilon = 1e-9);
        assert_relative_eq!(y, 0.0, epsilon = 1e-9);

        let (x, y) = transformer.lat_lon_to_meters(10.0, 10.0);
        assert_relative_eq!(x, 1113194.91, epsilon = 1e-2);
        assert_relative_eq!(y, 1118889.97, epsilon = 1e-2);

        // These are the max lat lon for Google maps. More North or South and numerical errors occur
        let (x, y) = transformer.lat_lon_to_meters(179.999_996_920_671_83, 85.051_128_514_163);
        assert_relative_eq!(x, 20037508.0, epsilon = 1e-2);
        assert_relative_eq!(y, 20037508.0, epsilon = 1e-2);

        let (x, y) = transformer.lat_lon_to_meters(180.0, 85.0);
        assert_relative_eq!(x, 20037508.34, epsilon = 5e-3);
        assert_relative_eq!(y, 19971868.88, epsilon = 5e-3);
        //
        let (x, y) = transformer.lat_lon_to_meters(-180.0, -85.0);
        assert_relative_eq!(x, -20037508.34, epsilon = 5e-3);
        assert_relative_eq!(y, -19971868.88, epsilon = 5e-3);
    }

    #[test]
    fn meters_to_lat_lon() {
        let transformer = GlobalMercator::new(256);
        let (lon, lat) = transformer.meters_2_lat_lon(0.0, 0.0);
        assert_relative_eq!(lon, 0.0, epsilon = f64::EPSILON);
        assert_relative_eq!(lat, 0.0, epsilon = f64::EPSILON);

        let (lon, lat) = transformer.meters_2_lat_lon(-20037508.0, -20037508.0);
        assert_relative_eq!(lon, -179.999997, epsilon = 1e-6);
        assert_relative_eq!(lat, -85.051129, epsilon = 1e-6);

        let (lon, lat) = transformer.meters_2_lat_lon(20037508.0, 20037508.0);
        assert_relative_eq!(lon, 179.999997, epsilon = 1e-6);
        assert_relative_eq!(lat, 85.051129, epsilon = 1e-6);
    }

    #[test]
    fn tile_bounds() {
        let transformer = GlobalMercator::new(256);

        let c = 20_037_508.34;

        for z in 0..MAX_ZOOM {
            let num_tiles = 2u64.pow(z as u32);
            let m = c * 2.0 / (num_tiles as f64) - c;

            let bounds = transformer.tile_bounds(0, 0, z);
            assert_relative_eq!(bounds.0, -c, epsilon = 5e-3);
            assert_relative_eq!(bounds.1, -c, epsilon = 5e-3);
            assert_relative_eq!(bounds.2, m, epsilon = 5e-3);
            assert_relative_eq!(bounds.3, m, epsilon = 5e-3);

            let bounds = transformer.tile_bounds(num_tiles - 1, 0, z);
            assert_relative_eq!(bounds.0, -m, epsilon = 5e-3);
            assert_relative_eq!(bounds.1, -c, epsilon = 5e-3);
            assert_relative_eq!(bounds.2, c, epsilon = 5e-3);
            assert_relative_eq!(bounds.3, m, epsilon = 5e-3);
        }

    }
}