#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Simple script to validate the test_noise.py test suite.
This checks that the test file is syntactically correct and can be imported.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def check_syntax():
    """Check that the test file has valid Python syntax."""
    test_file = os.path.join(os.path.dirname(__file__), 'test_noise.py')
    
    print("Checking syntax of test_noise.py...")
    with open(test_file, 'r') as f:
        code = f.read()
    
    try:
        compile(code, test_file, 'exec')
        print("✓ Syntax check passed")
        return True
    except SyntaxError as e:
        print(f"✗ Syntax error: {e}")
        return False

def check_imports():
    """Check which dependencies are available."""
    print("\nChecking available dependencies:")
    
    dependencies = {
        'numpy': False,
        'scipy': False,
        'torch': False,
        'lal': False,
        'lalsimulation': False,
        'gwpy': False,
        'bilby': False,
    }
    
    for dep in dependencies:
        try:
            __import__(dep)
            dependencies[dep] = True
            print(f"  ✓ {dep} available")
        except ImportError:
            print(f"  ✗ {dep} not available")
    
    return dependencies

def check_test_structure():
    """Verify test structure without importing (to avoid dependency errors)."""
    test_file = os.path.join(os.path.dirname(__file__), 'test_noise.py')
    
    print("\nChecking test structure...")
    with open(test_file, 'r') as f:
        content = f.read()
    
    # Check for test classes
    test_classes = [
        'TestLALSimulationPSD',
        'TestAdvancedLIGO',
        'TestKnownPSDs',
        'TestBilbyComparison'
    ]
    
    for test_class in test_classes:
        if f"class {test_class}" in content:
            print(f"  ✓ Found test class: {test_class}")
        else:
            print(f"  ✗ Missing test class: {test_class}")
    
    # Count test methods
    import re
    test_methods = re.findall(r'def (test_\w+)\(self\)', content)
    print(f"\n  Found {len(test_methods)} test methods:")
    
    # Key tests we expect
    key_tests = [
        'test_generated_noise_psd_matches_input_psd_long_duration',
        'test_generated_noise_psd_matches_input_psd_multiple_realizations',
        'test_noise_mean_comparison_with_bilby',
    ]
    
    for key_test in key_tests:
        if key_test in test_methods:
            print(f"    ✓ {key_test}")
        else:
            print(f"    ✗ {key_test} (missing)")
    
    print(f"\n  Total test methods: {len(test_methods)}")

def try_import_tests():
    """Try to import the test module (will fail if dependencies missing)."""
    print("\nAttempting to import test_noise module...")
    try:
        import test_noise
        print("✓ Successfully imported test_noise")
        
        # Try to list test methods
        import unittest
        loader = unittest.TestLoader()
        suite = loader.loadTestsFromModule(test_noise)
        test_count = suite.countTestCases()
        print(f"✓ Found {test_count} test cases")
        return True
    except ImportError as e:
        print(f"✗ Import failed: {e}")
        print("  (This is expected if dependencies are not installed)")
        return False

def main():
    """Run all validation checks."""
    print("=" * 60)
    print("Validating test_noise.py test suite")
    print("=" * 60)
    
    # Check syntax first (should always pass)
    if not check_syntax():
        print("\n❌ Syntax check failed - fix syntax errors first")
        return 1
    
    # Check dependencies
    deps = check_imports()
    
    # Check structure
    check_test_structure()
    
    # Try to import (might fail due to missing deps)
    import_ok = try_import_tests()
    
    print("\n" + "=" * 60)
    print("Validation Summary")
    print("=" * 60)
    print("✓ Test file syntax is valid")
    print(f"  {sum(deps.values())}/{len(deps)} dependencies available")
    
    if import_ok:
        print("✓ Tests can be imported and run")
        print("\nTo run tests:")
        print("  python -m unittest tests.test_noise -v")
    else:
        print("⚠ Tests cannot be imported due to missing dependencies")
        print("\nInstall dependencies with:")
        print("  pip install numpy scipy torch gwpy")
        print("  conda install -c conda-forge lalsuite")
        print("  pip install bilby  # optional")
    
    print("=" * 60)
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
