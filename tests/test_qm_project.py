"""
Unit and regression test for the qm_project package.
"""

# Import package, test suite, and other packages as needed
import qm_project
import pytest
import sys

def test_qm_project_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "qm_project" in sys.modules
