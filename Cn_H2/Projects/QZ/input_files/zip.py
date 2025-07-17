import zipfile
import os

def zip_folder(folder_path, output_zip, compression=zipfile.ZIP_DEFLATED, level=6):
    """
    Compresses a folder (including subdirectories) into a zip file.

    Args:
        folder_path (str): Path to the folder to compress.
        output_zip (str): Output zip file (e.g. "archive.zip").
        compression: One of zipfile.ZIP_STORED, ZIP_DEFLATED, ZIP_BZIP2, ZIP_LZMA.
        level (int): Compression level:
                     - For DEFLATED: 0–9 (0=no compression, 9=max; default -1 ≈6)
                     - For BZIP2 or LZMA: 1–9 (will be ignored for STORED/LZMA has no effect)
    """
    with zipfile.ZipFile(output_zip, 'w',
                         compression=compression,
                         compresslevel=level) as zf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                abs_path = os.path.join(root, file)
                rel_path = os.path.relpath(abs_path, start=folder_path)
                zf.write(abs_path, arcname=rel_path)
    print(f"Created '{output_zip}' from '{folder_path}' with level={level}")
# Fastest (least compression):
zip_folder('Molpro_CP', 'Molpro_CP.zip', level=9)
