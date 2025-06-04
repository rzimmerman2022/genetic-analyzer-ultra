import sys
import os

# Ensure the project root (c:/DNA) is in sys.path so modules can be imported by tests
project_root = os.path.abspath(os.path.dirname(__file__))
if project_root not in sys.path:
    sys.path.insert(0, project_root)
