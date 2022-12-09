try:
    import mass
except ImportError:
    mass = None


def test_mass_import():
    assert mass is not None
