from chgnet.utils import read_json
from pymatgen.core import Structure
from chgnet.data.dataset import StructureData, get_train_val_test_loader
from chgnet.model import CHGNet
from chgnet.trainer import Trainer


dataset_dict = read_json("./chgnet_dataset.json")
structures = [Structure.from_dict(struct) for struct in dataset_dict["structure"]]
energies = dataset_dict["energy_per_atom"]
forces = dataset_dict["force"]
stresses = dataset_dict.get("stress") or None
magmoms = dataset_dict.get("magmom") or None


dataset = StructureData(
    structures=structures,
    energies=energies,
    forces=forces,
    stresses=stresses,  # can be None
    magmoms=magmoms,  # can be None
)
train_loader, val_loader, test_loader = get_train_val_test_loader(
    dataset, batch_size=8, train_ratio=0.9, val_ratio=0.05
)


# Load pretrained CHGNet
chgnet = CHGNet.load()


# Optionally fix the weights of some layers
for layer in [
    chgnet.atom_embedding,
    chgnet.bond_embedding,
    chgnet.angle_embedding,
    chgnet.bond_basis_expansion,
    chgnet.angle_basis_expansion,
    chgnet.atom_conv_layers[:-1],
    chgnet.bond_conv_layers,
    chgnet.angle_layers,
]:
    for param in layer.parameters():
        param.requires_grad = False

#         Define Trainer
trainer = Trainer(
    model=chgnet,
    targets="efs",
    optimizer="Adam",
    scheduler="CosLR",
    criterion="MSE",
    epochs=20,
    learning_rate=1e-2,
    use_device="cuda",
    print_freq=6,
)
trainer.train(train_loader, val_loader, test_loader)

