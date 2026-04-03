"""
Sparse Autoencoder model.

Anthropic-style ReLU SAE with L1 sparsity penalty and unit-norm decoder
columns, following Bricken et al. (2023).
"""
import torch
import torch.nn as nn


class SparseAutoencoder(nn.Module):
    """ReLU sparse autoencoder for decomposing neural network activations.

    Args:
        d_input: Dimensionality of the input activations.
        d_hidden: Dictionary size (number of latent features).
    """

    def __init__(self, d_input: int, d_hidden: int):
        super().__init__()
        self.d_input = d_input
        self.d_hidden = d_hidden
        self.W_enc = nn.Linear(d_input, d_hidden, bias=True)
        self.W_dec = nn.Linear(d_hidden, d_input, bias=True)
        with torch.no_grad():
            self.W_dec.weight.data = nn.functional.normalize(
                self.W_dec.weight.data, dim=0
            )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        """Encode input activations into sparse feature activations."""
        return torch.relu(self.W_enc(x - self.W_dec.bias))

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """Decode sparse features back to activation space."""
        return self.W_dec(z)

    def forward(self, x: torch.Tensor):
        """Full forward pass: encode then decode.

        Returns:
            x_hat: Reconstructed activations.
            z: Sparse feature activations.
        """
        z = self.encode(x)
        x_hat = self.decode(z)
        return x_hat, z

    def loss(self, x: torch.Tensor, lam: float):
        """Compute SAE loss = MSE + lambda * L1.

        Returns:
            total_loss, mse, l1, z
        """
        x_hat, z = self.forward(x)
        mse = (x - x_hat).pow(2).mean()
        l1 = z.abs().mean()
        return mse + lam * l1, mse, l1, z

    def normalize_decoder(self):
        """Constrain decoder column norms to unity (call after each step)."""
        with torch.no_grad():
            self.W_dec.weight.data = nn.functional.normalize(
                self.W_dec.weight.data, dim=0
            )
