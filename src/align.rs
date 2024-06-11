use nalgebra::{DMatrix, Matrix3, Vector3, SVD};
use rand::seq::index;
use std::{char::MAX, f64};

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

#[derive(Debug, Clone)]
struct Afp {
    first: i32,
    second: i32,
}

impl Afp {
    fn new(first: i32, second: i32) -> Self {
        Afp { first, second }
    }
}

#[derive(Clone, Debug)]
struct Path {
    afps: Vec<Afp>,
}

const D0: f64 = 3.0;
const D1: f64 = 4.0;
// const MAX_KEPT: usize = 20;
const GAP_MAX: usize = 30;

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

fn get_coords(pdb_file: &str) -> Vec<Point> {
    let mut coords = Vec::new();

    let pdb = match pdbtbx::open_pdb(pdb_file, pdbtbx::StrictnessLevel::Loose) {
        Ok((pdb, _warnings)) => pdb,
        Err(e) => {
            panic!("Error opening PDB file: {:?}", e);
        }
    };

    for atom in pdb.atoms() {
        // Only get CA atoms
        if atom.name() != "CA" {
            continue;
        }

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
    d1: &[Vec<f64>],
    d2: &[Vec<f64>],
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

fn find_path(
    da: &[Vec<f64>],
    db: &[Vec<f64>],
    s: &[Vec<f64>],
    win_size: usize,
) -> Vec<Vec<(i32, i32)>> {
    let mut best_path_len = 0;
    let mut best_path_score = f64::MAX;
    // Longest possible alignment path
    let longest_possible_path = if da.len() < db.len() {
        da.len()
    } else {
        db.len()
    };

    let win_sum = (win_size - 1) * (win_size - 2) / 2;

    // Create a vector to store the values
    let mut win_cache = vec![0; longest_possible_path];

    // Fill the vector with the computed values
    for i in 0..longest_possible_path {
        win_cache[i] = ((i + 1) * i * win_size / 2) + (((i + 1) as i32 * win_sum as i32) as usize);
    }

    let mut best_path = Path {
        afps: vec![Afp::new(-1, -1); longest_possible_path],
    };

    let mut all_score_buffer: Vec<Vec<f64>> =
        vec![vec![1e6; GAP_MAX * 2 + 1]; longest_possible_path];

    let mut t_index: Vec<i32> = vec![0; longest_possible_path];

    for (i_a, _element_a) in da.iter().enumerate() {
        // // Limit how deep this loop can go;
        // if i_a as i32 > da.len() as i32 - (win_size * best_path.afps.len() as i32 - 1) {
        //     continue;
        // }

        for (i_b, _element_b) in db.iter().enumerate() {
            let similarity_score = s[i_a][i_b];
            let mut current_path_len = 1;
            t_index[current_path_len - 1] = 0;
            // Check all possible paths starting from iA, iB
            let mut path = Path {
                afps: vec![Afp::new(-1, -1); longest_possible_path],
            };

            path.afps[0] = Afp::new(i_a as i32, i_b as i32);

            let mut done = false;
            while !done {
                let mut gap_best_score = f64::MAX;
                let mut gap_best_index = -1;

                for g in 0..(GAP_MAX * 2) + 1 {
                    // println!("g: {}", g);
                    let mut j_a = path.afps[current_path_len - 1].first + win_size as i32;
                    let mut j_b = path.afps[current_path_len - 1].second + win_size as i32;

                    if j_a < 0 || j_b < 0 {
                        continue;
                    }

                    if ((g + 1) % 2) == 0 {
                        j_a += (g as i32 + 1) / 2;
                    } else {
                        j_b += (g as i32 + 1) / 2;
                    }

                    if j_a > da.len() as i32 - win_size as i32 - 1
                        || j_b > db.len() as i32 - win_size as i32 - 1
                    {
                        continue;
                    }

                    if s[j_a as usize][j_b as usize] > D0 {
                        continue;
                    }

                    if s[j_a as usize][j_b as usize] == -1.0 {
                        continue;
                    }

                    // Calculate the score of this path
                    let mut score = 0.0;

                    let anchor_a = path.afps[current_path_len - 1].first;
                    let anchor_b = path.afps[current_path_len - 1].second;

                    // If any of the indexes are out of bounds, skip this iteration
                    if j_a >= da.len() as i32 || j_b >= db.len() as i32 {
                        continue;
                    }

                    // Calculate the score
                    score += (da[anchor_a as usize][j_a as usize]
                        - db[anchor_b as usize][j_b as usize])
                        .abs();

                    score += (da[anchor_a as usize + win_size - 1][j_a as usize + win_size - 1]
                        - db[anchor_b as usize + win_size - 1][j_b as usize + win_size - 1])
                        .abs();

                    for k in 0..win_size - 1 {
                        score += (da[anchor_a as usize + k][j_a as usize + win_size - 1 - k]
                            - db[anchor_b as usize + k][j_b as usize + win_size - 1 - k])
                            .abs();
                    }

                    score /= win_size as f64 * current_path_len as f64;

                    // println!("score: {}", score);

                    if score < gap_best_score {
                        path.afps[current_path_len].first = j_a;
                        path.afps[current_path_len].second = j_b;
                        gap_best_score = score;
                        gap_best_index = g as i32;
                        all_score_buffer[current_path_len - 1][g] = score;
                    }
                } // end of gap search

                // Calculate current_total_score
                let mut current_total_score = 0.0;

                let g_a: i32;
                let g_b: i32;
                if gap_best_index != -1 {
                    let j_gap = (gap_best_index + 1) / 2;
                    if ((gap_best_index + 1) % 2) == 0 {
                        g_a = path.afps[current_path_len - 1].first + win_size as i32 + j_gap;
                        g_b = path.afps[current_path_len - 1].second + win_size as i32;
                    } else {
                        g_a = path.afps[current_path_len - 1].first + win_size as i32;
                        g_b = path.afps[current_path_len - 1].second + win_size as i32 + j_gap;
                    }

                    let score1 = (all_score_buffer[current_path_len - 1][gap_best_index as usize]
                        * win_size as f64
                        * current_path_len as f64
                        + s[g_a as usize][g_b as usize] * win_size as f64)
                        / (win_sum as f64 * current_path_len as f64 + win_sum as f64);

                    let score2 = if current_path_len > 1 {
                        all_score_buffer[current_path_len - 2]
                            [t_index[current_path_len - 1] as usize]
                    } else {
                        s[i_a][i_b]
                    } * win_cache[current_path_len - 1] as f64
                        + score1
                            * (win_cache[current_path_len] - win_cache[current_path_len - 1])
                                as f64
                            / win_cache[current_path_len] as f64;

                    current_total_score = score2;

                    if current_total_score > D1 {
                        done = true;
                        gap_best_index -= 1;
                        continue;
                    } else {
                        all_score_buffer[current_path_len - 1][gap_best_index as usize] =
                            current_total_score;
                        t_index[current_path_len] = gap_best_index;
                        current_path_len += 1;
                    }
                } else {
                    done = true;
                    current_path_len -= 1;
                    continue;
                }

                if current_path_len > best_path_len
                    || (current_path_len == best_path_len && current_total_score < best_path_score)
                {
                    best_path_len = current_path_len;
                    best_path_score = current_total_score;
                    best_path = path.clone();

                    // for i in 0..current_path_len {
                    //     best_path.afps[i] = path.afps[i];
                    // }
                }

                done = true;
            } // end loop
        } // end ib
    } // end ia

    // Loop over the best paths and expand them into continuous lists
    let mut result: Vec<Vec<(i32, i32)>> = Vec::new();
    for p in best_path.afps.iter() {
        // println!("{:?}", p);
        let mut data: Vec<(i32, i32)> = Vec::new();
        if p.first == -1 || p.second == -1 {
            continue;
        }
        for k in 0..win_size {
            println!("{} {}", p.first + k as i32, p.second + k as i32);
            data.push((p.first + k as i32, p.second + k as i32));
        }
        result.push(data);
    }
    result
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

    let rmsd = e0 - 2.0 * s.iter().sum::<f64>();
    let rmsd = (rmsd.abs() / l as f64).sqrt();

    if align == 0 {
        return rmsd;
    }

    assert!(mol1.is_some());
    assert!(mol2.is_some());

    // U is simply V*Wt
    // let u = v * wt;

    // // rotate and translate the molecule
    // let mol2 = mol2
    //     .unwrap()
    //     .iter()
    //     .map(|p| {
    //         let np = nalgebra::DVector::from_vec(vec![p.x - com2.x, p.y - com2.y, p.z - com2.z]);
    //         let np = u.clone() * np + nalgebra::DVector::from_vec(vec![com1.x, com1.y, com1.z]);
    //         Point::new(np[0], np[1], np[2])
    //     })
    //     .collect::<Vec<_>>();

    // let PyMol know about the changes to the coordinates
    // cmd.alter_state(1, name2, "(x,y,z)=stored.sel2.pop(0)");

    println!("NumAligned={}", l);
    println!("RMSD={}", rmsd);

    rmsd
}
pub fn align(input_a: &str, input_b: &str) -> Result<(), Box<dyn std::error::Error>> {
    println!("Reading input file: {}", input_a);
    println!("Reading input file: {}", input_b);

    // 1. Get CA coordinates, these are added to a Vec<Point> (Python in cealign.py)
    let coords_a = get_coords(input_a);
    let coords_b = get_coords(input_b);

    // 2. Calculate the Distance Matrix for each protein
    let dm_a = calc_distance_matrix(&coords_a);
    let dm_b = calc_distance_matrix(&coords_b);

    // 3. Compute the Similarity Matrix;
    //  - The function iterates over all possible starting positions iA in d1 and iB in d2.
    //  - For each pair of starting positions, it initializes S[iA][iB] to -1.0 as a default value.
    //  - It checks if the current starting positions allow a complete segment of size winSize to fit within the bounds of d1 and d2. If not, it continues to the next pair.
    //  - For valid segments, it calculates a similarity score by comparing the distances between pairs of residues within the window. It skips distances between neighboring residues (a heuristic decision to ignore these commonly similar distances).
    //  - The score is computed as the sum of the absolute differences between the corresponding distances in d1 and d2, normalized by sumSize.
    //  - After computing the similarity scores for all pairs of segments, the function returns the 2D array S, which contains the normalized similarity scores.
    let s = calc_s(&dm_a, &dm_b, dm_a.len(), dm_b.len(), 8);

    // 4. Find the best path through the similarity matrix;
    //  - The function iterates over all possible starting positions iA and iB in the similarity matrix S.
    //  - For each pair (iA, iB), it checks if the similarity score S[iA][iB] is below the threshold D0 and not equal to -1.0.
    //  - If the starting positions are valid, it initializes a new path curPath starting from (iA, iB).
    //  - The function then searches for the best path starting from this position, considering all possible gaps up to gapMax.
    // let paths = find_path(&s, &dm_a, &dm_b, 8);
    let aln = find_path(&dm_a, &dm_b, &s, 8);

    println!("{:?}", aln);
    let bestPathScore = 100000;
    for c_aln in aln {
        // let seq_count = cur_alignment.len();
        let mut mat_a: Option<Vec<Point>> = None;
        let mut mat_b: Option<Vec<Point>> = None;

        // if seq_count == 0 {
        //     continue;
        // }
        for afp in &c_aln {
            // for afp in cur_alignment.iter() {
            // let (first, second) = afp;
            if mat_a.is_none() && mat_b.is_none() {
                mat_a = Some(vec![coords_a[(afp.0 - 1) as usize]]);
                mat_b = Some(vec![coords_b[(afp.1 - 1) as usize]]);
            } else {
                mat_a.as_mut().unwrap().push(coords_a[(afp.0 - 1) as usize]);
                mat_b.as_mut().unwrap().push(coords_b[(afp.1 - 1) as usize]);
            }
            // }
            let rmsd = simp_align(
                mat_a.as_ref().unwrap(),
                mat_b.as_ref().unwrap(),
                input_a,
                input_b,
                Some(&coords_a),
                Some(&coords_b),
                1,
                mat_a.as_ref().unwrap().len(),
            );

            //     internalGaps = 0.0;
            // for g in range(0, seqCount-1):
            // 	if (not curAlignment[g][0] + 1 == curAlignment[g+1][0]):
            // 		internalGaps += curAlignment[g+1][0]
            // 	if ( not curAlignment[g][1] + 1 == curAlignment[g+1][1] ):
            // 		internalGaps += curAlignment[g+1][1]

            // 	aliLen = float( len(curAlignment))
            // 	numGap = internalGaps;
            // 	curScore = float((curScore/aliLen)*(1.0+(numGap/aliLen)));
            // let mut internal_gaps = 0.0;
            // for g in 0..c_aln.len() - 1 {
            //     if c_aln[g].0 + 1 != c_aln[g + 1].0 {
            //         internal_gaps += afp[g + 1][0]
            //     }
            //     if c_aln[g].1 + 1 != c_aln[g + 1].1 {
            //         // internalGaps += curAlignment[g+1][1]
            //     }
            // }
        }
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
