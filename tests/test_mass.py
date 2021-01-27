try:
    import mass
except ImportError:
    mass = None


def test_mass_import():
    return mass is not None
