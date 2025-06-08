from debyecalculator import DebyeCalculator
import torch
import matplotlib.pyplot as plt


# Initialize the DebyeCalculator object
calc = DebyeCalculator(qmin=1.0, qmax=20.0, qstep=0.05, qdamp=0.0,biso=0.0,device='cpu')


# Define structure sources
xyz_file = "output.xyz"
q, iq= calc.iq(xyz_file)

å

# Plot
fig, ax = plt.subplots(figsize=(6,3))
ax.plot(q, iq, label=f"Simulated XRD for {xyz_file.split('/')[-1]}")
ax.set(xlabel='Q [$Å^{-1}$]', ylabel='I(Q) [a.u.]', yticks=[])
ax.grid(alpha=0.2)
ax.legend()
plt.show()
plt.savefig("IQ.png")
