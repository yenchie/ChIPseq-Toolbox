# folder.structure.fun.R
# 20240515
# YCW

print_folder_structure <- function(path, prefix = "", level = 0, max_level = Inf, show_files = TRUE, show_folders = TRUE, max_files = 3) {
    # Stop if the current level exceeds max_level
    if (level > max_level) {
        return()
    }

    # Get the list of files and directories
    items <- list.files(path, full.names = TRUE)

    # Separate directories and files
    dirs <- items[dir.exists(items)]
    files <- items[!dir.exists(items)]

    # Print directories
    if (show_folders) {
        for (dir in dirs) {
            cat(prefix, "|--", basename(dir), "\n", sep = "")
            # Recursively print the contents of the directory with updated prefix and level
            print_folder_structure(dir, paste0(prefix, "   "), level + 1, max_level, show_files, show_folders, max_files)
        }
    }

    # Print files
    if (show_files && level < max_level) {
        num_files <- length(files)
        for (i in seq_len(min(num_files, max_files))) {
            cat(prefix, "|--", basename(files[i]), "\n", sep = "")
        }
        if (num_files > max_files) {
            cat(prefix, "...\n")
        }
    }
}
