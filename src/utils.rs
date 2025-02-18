use std::fs::File;
use std::io::Write;
use std::path::Path;

pub fn write_string_to_file(content: &str, file_path: &str) -> std::io::Result<()> {
    let path = Path::new(file_path);

    let mut file =
        File::create(path).unwrap_or_else(|e| panic!("Could not create {}: {}", file_path, e));

    file.write_all(content.as_bytes())
        .unwrap_or_else(|e| panic!("Could not write file {}: {}", file_path, e));

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::path::Path;

    #[test]
    fn test_write_string_to_file() {
        let temp_file = "test_file.txt";
        let test_content = "test";

        // Clean up any existing test file
        if Path::new(temp_file).exists() {
            fs::remove_file(temp_file).expect("Failed to remove existing test file");
        }

        // Test the function
        let result = write_string_to_file(test_content, temp_file);
        assert!(result.is_ok());

        // Verify content was written correctly
        let read_content = fs::read_to_string(temp_file).expect("Failed to read test file");
        assert_eq!(read_content, test_content);

        // Clean up
        fs::remove_file(temp_file).expect("Failed to clean up test file");
    }

    #[test]
    fn test_write_string_to_nonexistent_directory() {
        let invalid_path = "nonexistent_dir/test_file.txt";
        let test_content = "Hello, world!";

        // This should panic with our custom error message
        let result = std::panic::catch_unwind(|| write_string_to_file(test_content, invalid_path));

        assert!(result.is_err());
    }
}
