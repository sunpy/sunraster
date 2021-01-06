import os.path

import sunraster

__all__ = ["test_data_path"]

test_data_dir = os.path.join(os.path.dirname(sunraster.__file__), "tests", "data")
