import hashlib
from pathlib import Path
from typing import Union


def calculate_checksum(filename: Union[str, Path]):
    assert Path(filename).exists(), f"file '{filename}' does not exist!"
    with open(filename, 'rb') as file_to_check:
        # read contents of the file
        data = file_to_check.read()
        # pipe contents of the file through
        md5_checksum = hashlib.md5(data).hexdigest()
        return md5_checksum


def write_checksums():
    videos_path = Path(__file__).parent.parent.joinpath('rnalysis/gui/videos')
    checksums_path = videos_path.joinpath('checksums')

    for item in videos_path.iterdir():
        if item.is_file() and item.suffix == '.webp':
            checksum = calculate_checksum(item)
            with open(checksums_path.joinpath(f'{item.stem}.txt'), 'w') as f:
                f.write(checksum)


if __name__ == '__main__':
    write_checksums()
