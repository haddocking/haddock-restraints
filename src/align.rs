#[derive(Debug, Clone, Copy)]
struct Point {
    x: f64,
    y: f64,
    z: f64,
}

impl Point {
    fn new(x: f64, y: f64, z: f64) -> Self {
        Point { x, y, z }
    }
}

use nalgebra::{DMatrix, Matrix3, Vector3, SVD};
use std::f64;
const D0: f64 = 3.0;
const D1: f64 = 4.0;
const MAX_KEPT: usize = 20;
const GAP_MAX: usize = 30;

// #[derive(Clone, Debug)]
// struct Afp {
//     first: i32,
//     second: i32,
// }

fn translate_to_origin(points: &[Point], centroid: Point) -> Vec<Point> {
    points
        .iter()
        .map(|p| Point::new(p.x - centroid.x, p.y - centroid.y, p.z - centroid.z))
        .collect()
}

fn centroid(points: &[Point]) -> Point {
    let mut x_sum = 0.0;
    let mut y_sum = 0.0;
    let mut z_sum = 0.0;

    for point in points {
        x_sum += point.x;
        y_sum += point.y;
        z_sum += point.z;
    }

    let len = points.len() as f64;
    Point::new(x_sum / len, y_sum / len, z_sum / len)
}

fn covariance_matrix(points1: &[Point], points2: &[Point]) -> [[f64; 3]; 3] {
    let mut cov_matrix = [[0.0; 3]; 3];

    for (p1, p2) in points1.iter().zip(points2.iter()) {
        cov_matrix[0][0] += p1.x * p2.x;
        cov_matrix[0][1] += p1.x * p2.y;
        cov_matrix[0][2] += p1.x * p2.z;

        cov_matrix[1][0] += p1.y * p2.x;
        cov_matrix[1][1] += p1.y * p2.y;
        cov_matrix[1][2] += p1.y * p2.z;

        cov_matrix[2][0] += p1.z * p2.x;
        cov_matrix[2][1] += p1.z * p2.y;
        cov_matrix[2][2] += p1.z * p2.z;
    }

    cov_matrix
}

fn optimal_rotation(cov_matrix: [[f64; 3]; 3]) -> Matrix3<f64> {
    let matrix = Matrix3::from(cov_matrix);
    let svd = SVD::new(matrix, true, true);
    let u = svd.u.unwrap();
    let v_t = svd.v_t.unwrap();

    let det = (u * v_t).determinant();
    let d = Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, det);
    u * d * v_t
}

fn apply_rotation(points: &[Point], rotation: Matrix3<f64>) -> Vec<Point> {
    points
        .iter()
        .map(|p| {
            let vec = Vector3::new(p.x, p.y, p.z);
            let rotated_vec = rotation * vec;
            Point::new(rotated_vec.x, rotated_vec.y, rotated_vec.z)
        })
        .collect()
}

fn get_coords(pdb_file: &str) -> Vec<Point> {
    let mut coords = Vec::new();

    let pdb = match pdbtbx::open_pdb(pdb_file, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    for atom in pdb.atoms() {
        coords.push(Point {
            x: atom.x(),
            y: atom.y(),
            z: atom.z(),
        });
    }

    coords
}

fn calc_distance_matrix(coords: &[Point]) -> Vec<Vec<f64>> {
    let len = coords.len();
    let mut dm = vec![vec![0.0; len]; len];

    for row in 0..len {
        for col in 0..len {
            dm[row][col] = ((coords[row].x - coords[col].x).powi(2)
                + (coords[row].y - coords[col].y).powi(2)
                + (coords[row].z - coords[col].z).powi(2))
            .sqrt();
        }
    }

    dm
}

fn calc_s(
    d1: &Vec<Vec<f64>>,
    d2: &Vec<Vec<f64>>,
    len_a: usize,
    len_b: usize,
    win_size: usize,
) -> Vec<Vec<f64>> {
    let sum_size = ((win_size - 1) * (win_size - 2)) as f64 / 2.0;

    // Initialize the 2D similarity matrix
    let mut s = vec![vec![-1.0; len_b]; len_a];

    for i_a in 0..len_a {
        for i_b in 0..len_b {
            if i_a > len_a - win_size || i_b > len_b - win_size {
                continue;
            }

            let mut score = 0.0;

            for row in 0..win_size - 2 {
                for col in row + 2..win_size {
                    score += (d1[i_a + row][i_a + col] - d2[i_b + row][i_b + col]).abs();
                }
            }

            s[i_a][i_b] = score / sum_size;
        }
    }

    s
}

fn kabsch_algorithm(points1: &[Point], points2: &[Point]) -> Vec<Point> {
    assert_eq!(
        points1.len(),
        points2.len(),
        "Point sets must have the same length"
    );

    let centroid1 = centroid(points1);
    let centroid2 = centroid(points2);

    let points1_centered = translate_to_origin(points1, centroid1);
    let points2_centered = translate_to_origin(points2, centroid2);

    let cov_matrix = covariance_matrix(&points1_centered, &points2_centered);

    let rotation = optimal_rotation(cov_matrix);

    let rotated = apply_rotation(&points1_centered, rotation);

    // Translate back to original position
    rotated
        .iter()
        .map(|p| Point::new(p.x + centroid2.x, p.y + centroid2.y, p.z + centroid2.z))
        .collect()
}
use std::ptr;

// // Define the path and afp structures
// struct Path {
//     afps: Vec<Afp>,
// }

struct Afp {
    first: i32,
    second: i32,
}

#[derive(Clone)] // Add this line to derive the Clone trait for Path
struct Path {
    afps: Vec<Afp>,
}

impl Clone for Afp {
    fn clone(&self) -> Self {
        Self {
            first: self.first,
            second: self.second,
        }
    }
}

fn find_path(
    s: &Vec<Vec<f64>>,
    da: &Vec<Vec<f64>>,
    db: &Vec<Vec<f64>>,
    len_a: usize,
    len_b: usize,
    win_size: usize,
) -> Vec<Vec<(i32, i32)>> {
    // CE-specific cutoffs
    const D0: f64 = 3.0;
    const D1: f64 = 4.0;
    const MAX_KEPT: usize = 20;
    const GAP_MAX: usize = 30;

    // the best Path's score
    let mut best_path_score = 1e6;
    let mut best_path_length = 0;

    // length of longest possible alignment
    let smaller = if len_a < len_b { len_a } else { len_b };
    let win_sum = (win_size - 1) * (win_size - 2) / 2;

    // BEST PATH
    let mut best_path = Path {
        afps: vec![
            Afp {
                first: -1,
                second: -1
            };
            smaller
        ],
    };

    // for storing the best 20 paths
    let mut buffer_index = 0;
    let mut buffer_size = 0;
    let mut len_buffer = vec![0; MAX_KEPT];
    let mut score_buffer = vec![1e6; MAX_KEPT];
    let mut path_buffer: Vec<Option<Path>> = vec![None; MAX_KEPT];

    // winCache
    let mut win_cache = Vec::with_capacity(smaller);
    for i in 0..smaller {
        win_cache.push((i + 1) * i * win_size / 2 + (i + 1) * win_sum);
    }

    // allScoreBuffer
    let mut all_score_buffer = Vec::with_capacity(smaller);
    for _ in 0..smaller {
        let mut buffer = vec![1e6; GAP_MAX * 2 + 1];
        all_score_buffer.push(buffer);
    }

    let mut t_index = vec![0; smaller];
    let mut gap_best_index = -1;

    // Start the search through the CE matrix.
    for i_a in 0..len_a {
        if i_a > len_a - win_size * (best_path_length - 1) {
            break;
        }

        for i_b in 0..len_b {
            if s[i_a][i_b] >= D0 {
                continue;
            }

            if s[i_a][i_b] == -1.0 {
                continue;
            }

            if i_b > len_b - win_size * (best_path_length - 1) {
                break;
            }

            let mut cur_path = Path {
                afps: vec![
                    Afp {
                        first: -1,
                        second: -1
                    };
                    smaller
                ],
            };

            cur_path.afps[0] = Afp {
                first: i_a as i32,
                second: i_b as i32,
            };
            let mut cur_path_length = 1;
            t_index[cur_path_length - 1] = 0;
            let mut cur_total_score = 0.0;
            let mut done = false;

            while !done {
                let mut gap_best_score = 1e6;
                gap_best_index = -1;

                for g in 0..(GAP_MAX * 2) + 1 {
                    let mut j_a = cur_path.afps[cur_path_length - 1].first + win_size as i32;
                    let mut j_b = cur_path.afps[cur_path_length - 1].second + win_size as i32;

                    if (g + 1) % 2 == 0 {
                        j_a += (g + 1) as i32 / 2;
                    } else {
                        j_b += (g + 1) as i32 / 2;
                    }

                    if j_a > len_a as i32 - win_size as i32 - 1
                        || j_b > len_b as i32 - win_size as i32 - 1
                    {
                        continue;
                    }

                    if s[j_a as usize][j_b as usize] > D0 {
                        continue;
                    }

                    if s[j_a as usize][j_b as usize] == -1.0 {
                        continue;
                    }

                    let mut cur_score = 0.0;
                    for s in 0..cur_path_length {
                        cur_score += (da[cur_path.afps[s].first as usize][j_a as usize]
                            - db[cur_path.afps[s].second as usize][j_b as usize])
                            .abs();
                        cur_score += (da[cur_path.afps[s].first as usize + win_size - 1]
                            [j_a as usize + win_size - 1]
                            - db[cur_path.afps[s].second as usize + win_size - 1]
                                [j_b as usize + win_size - 1])
                            .abs();
                        for k in 1..win_size - 1 {
                            cur_score += (da[cur_path.afps[s].first as usize + k]
                                [j_a as usize + win_size - 1 - k]
                                - db[cur_path.afps[s].second as usize + k]
                                    [j_b as usize + win_size - 1 - k])
                                .abs();
                        }
                    }

                    cur_score /= win_size as f64 * cur_path_length as f64;

                    if cur_score >= D1 {
                        continue;
                    }

                    if cur_score < gap_best_score {
                        cur_path.afps[cur_path_length].first = j_a;
                        cur_path.afps[cur_path_length].second = j_b;
                        gap_best_score = cur_score;
                        gap_best_index = g as i32;
                        all_score_buffer[cur_path_length - 1][g] = cur_score;
                    }
                }
                /// ROF -- END GAP SEARCHING
                if gap_best_index != -1 {
                    let j_gap = (gap_best_index + 1) / 2;
                    let g_a;
                    let g_b;
                    if (gap_best_index + 1) % 2 == 0 {
                        g_a = cur_path.afps[cur_path_length - 1].first + win_size as i32 + j_gap;
                        g_b = cur_path.afps[cur_path_length - 1].second + win_size as i32;
                    } else {
                        g_a = cur_path.afps[cur_path_length - 1].first + win_size as i32;
                        g_b = cur_path.afps[cur_path_length - 1].second + win_size as i32 + j_gap;
                    }

                    let score1 = (all_score_buffer[cur_path_length - 1][gap_best_index as usize]
                        * win_size as f64
                        * cur_path_length as f64
                        + s[g_a as usize][g_b as usize] * win_sum as f64)
                        / (win_size as f64 * cur_path_length as f64 + win_sum as f64);

                    let score2 = if cur_path_length > 1 {
                        (all_score_buffer[cur_path_length - 2][t_index[cur_path_length - 1]]
                            * win_cache[cur_path_length - 1] as f64
                            + score1
                                * (win_cache[cur_path_length] - win_cache[cur_path_length - 1])
                                    as f64)
                            / win_cache[cur_path_length] as f64
                    } else {
                        (s[i_a][i_b] * win_cache[cur_path_length - 1] as f64
                            + score1
                                * (win_cache[cur_path_length] - win_cache[cur_path_length - 1])
                                    as f64)
                            / win_cache[cur_path_length] as f64
                    };

                    cur_total_score = score2;

                    if cur_total_score > D1 {
                        done = true;
                        gap_best_index = -1;
                        break;
                    } else {
                        all_score_buffer[cur_path_length - 1][gap_best_index as usize] =
                            cur_total_score;
                        t_index[cur_path_length] = gap_best_index as usize;
                        cur_path_length += 1;
                    }
                } else {
                    done = true;
                    cur_path_length -= 1;
                    break;
                }

                if cur_path_length > best_path_length
                    || (cur_path_length == best_path_length && cur_total_score < best_path_score)
                {
                    best_path_length = cur_path_length;
                    best_path_score = cur_total_score;
                    best_path.afps = cur_path.afps.clone();
                }
            }
            /// END WHILE
            if best_path_length > len_buffer[buffer_index]
                || (best_path_length == len_buffer[buffer_index]
                    && best_path_score < score_buffer[buffer_index])
            {
                buffer_index = if buffer_index == MAX_KEPT - 1 {
                    0
                } else {
                    buffer_index + 1
                };
                buffer_size = if buffer_size < MAX_KEPT {
                    buffer_size + 1
                } else {
                    MAX_KEPT
                };
                let path_copy = Path {
                    afps: best_path.afps.clone(),
                };

                if buffer_index == 0 && buffer_size == MAX_KEPT {
                    path_buffer[MAX_KEPT - 1] = Some(path_copy);
                    score_buffer[MAX_KEPT - 1] = best_path_score;
                    len_buffer[MAX_KEPT - 1] = best_path_length;
                } else {
                    path_buffer[buffer_index - 1] = Some(path_copy);
                    score_buffer[buffer_index - 1] = best_path_score;
                    len_buffer[buffer_index - 1] = best_path_length;
                }
            }
        }
    }

    let mut r_val = Vec::new();
    let mut i = 0;

    while i < buffer_size {
        let mut cur_list = Vec::new();
        let mut j = 0;

        while j < smaller {
            if path_buffer[i].as_ref().unwrap().afps[j].first != -1 {
                let mut k = 0;
                while k < win_size {
                    cur_list.push((
                        path_buffer[i].as_ref().unwrap().afps[j].first + k as i32,
                        path_buffer[i].as_ref().unwrap().afps[j].second + k as i32,
                    ));
                    k += 1;
                }
            }
            j += 1;
        }
        r_val.push(cur_list);
        i += 1;
    }

    r_val
}

fn rmsd(points1: &[Point], points2: &[Point]) -> f64 {
    assert_eq!(
        points1.len(),
        points2.len(),
        "Point sets must have the same length"
    );

    let sum_of_squares: f64 = points1
        .iter()
        .zip(points2.iter())
        .map(|(p1, p2)| {
            let dx = p1.x - p2.x;
            let dy = p1.y - p2.y;
            let dz = p1.z - p2.z;
            dx * dx + dy * dy + dz * dz
        })
        .sum();

    (sum_of_squares / points1.len() as f64).sqrt()
}

fn simp_align(
    mat1: &[Point],
    mat2: &[Point],
    name1: &str,
    name2: &str,
    mol1: Option<&[Point]>,
    mol2: Option<&[Point]>,
    align: i32,
    l: usize,
) -> f64 {
    // check for consistency
    assert_eq!(mat1.len(), mat2.len());

    // must always center the two proteins to avoid affine transformations.
    // Center the two proteins to their selections.
    let com1 = mat1.iter().fold(Point::new(0.0, 0.0, 0.0), |acc, p| {
        Point::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
    });
    let com1 = Point::new(com1.x / l as f64, com1.y / l as f64, com1.z / l as f64);
    let com2 = mat2.iter().fold(Point::new(0.0, 0.0, 0.0), |acc, p| {
        Point::new(acc.x + p.x, acc.y + p.y, acc.z + p.z)
    });
    let com2 = Point::new(com2.x / l as f64, com2.y / l as f64, com2.z / l as f64);
    let mut mat1 = mat1
        .iter()
        .map(|p| Point::new(p.x - com1.x, p.y - com1.y, p.z - com1.z))
        .collect::<Vec<_>>();
    let mut mat2 = mat2
        .iter()
        .map(|p| Point::new(p.x - com2.x, p.y - com2.y, p.z - com2.z))
        .collect::<Vec<_>>();

    // Initial residual, see Kabsch.
    let e0 = mat1
        .iter()
        .fold(0.0, |acc, p| acc + p.x * p.x + p.y * p.y + p.z * p.z)
        + mat2
            .iter()
            .fold(0.0, |acc, p| acc + p.x * p.x + p.y * p.y + p.z * p.z);

    // This beautiful step provides the answer. V and Wt are the orthonormal
    // bases that when multiplied by each other give us the rotation matrix, U.
    // S, (Sigma, from SVD) provides us with the error! Isn't SVD great!
    let (v, s, wt) = {
        let m1 = nalgebra::DMatrix::from_iterator(
            3,
            mat1.len(),
            mat1.iter().flat_map(|p| vec![p.x, p.y, p.z]),
        );
        let m2 = nalgebra::DMatrix::from_iterator(
            3,
            mat2.len(),
            mat2.iter().flat_map(|p| vec![p.x, p.y, p.z]),
        );
        let m1t = m1.transpose();
        let svd = nalgebra::linalg::SVD::new(m1t * m2, true, true);
        (svd.u.unwrap(), svd.singular_values, svd.v_t.unwrap())
    };

    // we already have our solution, in the results from SVD.
    // we just need to check for reflections and then produce
    // the rotation. V and Wt are orthonormal, so their det's
    // are +/-1.

    let reflect = v.determinant() * wt.determinant();
    let s = if reflect == -1.0 {
        s.iter().map(|&x| -x).collect::<Vec<_>>()
    } else {
        s.iter().map(|&x| x).collect::<Vec<_>>()
    };
    println!("here!");

    let rmsd = e0 - 2.0 * s.iter().sum::<f64>();
    let rmsd = (rmsd.abs() / l as f64).sqrt();

    if align == 0 {
        return rmsd;
    }

    assert!(mol1.is_some());
    assert!(mol2.is_some());

    // U is simply V*Wt
    let u = v * wt;

    // rotate and translate the molecule
    let mol2 = mol2
        .unwrap()
        .iter()
        .map(|p| {
            let np = nalgebra::DVector::from_vec(vec![p.x - com2.x, p.y - com2.y, p.z - com2.z]);
            let np = u.clone() * np + nalgebra::DVector::from_vec(vec![com1.x, com1.y, com1.z]);
            Point::new(np[0], np[1], np[2])
        })
        .collect::<Vec<_>>();

    // let PyMol know about the changes to the coordinates
    // cmd.alter_state(1, name2, "(x,y,z)=stored.sel2.pop(0)");

    println!("NumAligned={}", l);
    println!("RMSD={}", rmsd);

    rmsd
}
pub fn align(input_a: &str, input_b: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("Reading input file: {}", input_a);
    println!("Reading input file: {}", input_b);

    let coords_a = get_coords(input_a);
    let coords_b = get_coords(input_b);

    let dm_a = calc_distance_matrix(&coords_a);
    let dm_b = calc_distance_matrix(&coords_b);

    let s = calc_s(&dm_a, &dm_b, dm_a.len(), dm_b.len(), 5);

    let paths = find_path(&s, &dm_a, &dm_b, dm_a.len(), dm_b.len(), 5);

    for cur_alignment in paths {
        let seq_count = cur_alignment.len();
        let mut mat_a: Option<Vec<Point>> = None;
        let mut mat_b: Option<Vec<Point>> = None;

        if seq_count == 0 {
            continue;
        }

        for afp in cur_alignment {
            let (first, second) = afp;
            if mat_a.is_none() && mat_b.is_none() {
                mat_a = Some(vec![coords_a[(first - 1) as usize]]);
                mat_b = Some(vec![coords_b[(second - 1) as usize]]);
            } else {
                mat_a.as_mut().unwrap().push(coords_a[(first - 1) as usize]);
                mat_b
                    .as_mut()
                    .unwrap()
                    .push(coords_b[(second - 1) as usize]);
            }
        }
        let score = simp_align(
            mat_a.as_ref().unwrap(),
            mat_b.as_ref().unwrap(),
            input_a,
            input_b,
            Some(&coords_a),
            Some(&coords_b),
            1,
            seq_count,
        );
    }

    // println!("{:?}", s.len());

    // for (i, p) in paths.iter().enumerate() {
    //     // println!("first: {:?} second: {:?}", p.first, p.second);
    //     println!("{} {:?}", i, p);
    // }

    // let points1 = vec![
    //     Point::new(1.0, 2.0, 3.0),
    //     Point::new(4.0, 5.0, 6.0),
    //     Point::new(7.0, 8.0, 9.0),
    // ];

    // let points2 = vec![
    //     Point::new(1.0, 2.0, 4.0),
    //     Point::new(4.0, 6.0, 5.0),
    //     Point::new(7.0, 9.0, 8.0),
    // ];

    // let centroid_p = centroid(&coords_a);

    // let ori_coords_a = translate_to_origin(&coords_a, centroid(&coords_a));

    // let aligned_points = kabsch_algorithm(&coords_a, &coords_b);

    // let rmsd_value = rmsd(&aligned_points, &coords_a);
    // println!("RMSD: {:.6}", rmsd_value);

    Ok(())
}
