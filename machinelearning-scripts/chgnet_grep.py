from chgnet.model.model import CHGNet
import numpy as np
from pymatgen.core import Structure
from chgnet.utils import parse_vasp_dir
dataset_dict = parse_vasp_dir(base_dir="./", save_path="./chgnet/chgnet_dataset.json"
)
