import numpy as np

from QuantyLF.QuantyLF import QuantyLF


def test_blend_spectra_normalizes_weights():
    quanty = QuantyLF.__new__(QuantyLF)
    spectra = [
        np.array([[0.0, 2.0], [1.0, 4.0]]),
        np.array([[0.0, 6.0], [1.0, 8.0]]),
    ]

    blended = quanty._blend_spectra(spectra, weights=[1.0, 3.0])

    np.testing.assert_allclose(blended, [[0.0, 5.0], [1.0, 7.0]])
