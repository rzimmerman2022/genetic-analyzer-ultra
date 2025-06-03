#!/usr/bin/env python3
"""
GenoInsight Pro - Setup and Launch Script
=========================================
This script checks your environment and helps you run the genetic analysis.
"""

import subprocess
import sys
import os

def check_python_version():
    """Ensure Python 3.7+ is installed."""
    version = sys.version_info
    if version.major < 3 or (version.major == 3 and version.minor < 7):
        print("❌ Python 3.7 or higher is required")
        print(f"   You have: Python {version.major}.{version.minor}")
        return False
    print(f"✅ Python {version.major}.{version.minor} detected")
    return True

def check_package(package_name):
    """Check if a package is installed."""
    try:
        __import__(package_name)
        return True
    except ImportError:
        return False

def install_requirements():
    """Check and install required packages."""
    packages = {
        'pandas': 'Data manipulation',
        'numpy': 'Numerical computing',
        'matplotlib': 'Plotting',
        'seaborn': 'Statistical visualization',
        'scipy': 'Scientific computing',
        'sklearn': 'Machine learning',
        'networkx': 'Network analysis',
        'requests': 'HTTP requests'
    }
    
    print("\nChecking required packages...")
    missing = []
    
    for package, description in packages.items():
        if check_package(package):
            print(f"✅ {package:<12} - {description}")
        else:
            print(f"❌ {package:<12} - {description} (not installed)")
            missing.append(package)
    
    if missing:
        print(f"\nMissing {len(missing)} packages: {', '.join(missing)}")
        response = input("Would you like to install them now? (y/n): ")
        
        if response.lower() == 'y':
            print("\nInstalling packages...")
            # sklearn is installed as scikit-learn
            if 'sklearn' in missing:
                missing.remove('sklearn')
                missing.append('scikit-learn')
            
            cmd = [sys.executable, '-m', 'pip', 'install'] + missing
            subprocess.check_call(cmd)
            print("✅ All packages installed successfully!")
        else:
            print("\nPlease install the missing packages manually:")
            print(f"pip install {' '.join(missing)}")
            return False
    else:
        print("\n✅ All required packages are installed!")
    
    return True

def check_data_file():
    """Check for 23andMe data file."""
    print("\nChecking for genetic data files...")
    
    # Look for common 23andMe file patterns
    potential_files = []
    for file in os.listdir('.'):
        if file.endswith('.txt') and any(pattern in file.lower() for pattern in ['23andme', 'genome', 'raw', 'dna', 'genetic', 'paste']):
            potential_files.append(file)
    
    if potential_files:
        print(f"Found {len(potential_files)} potential data file(s):")
        for i, file in enumerate(potential_files, 1):
            size = os.path.getsize(file) / (1024 * 1024)  # Size in MB
            print(f"  {i}. {file} ({size:.1f} MB)")
        return True
    else:
        print("❌ No genetic data files found")
        print("\nTo get your data:")
        print("1. Log into 23andMe.com")
        print("2. Go to Settings → Privacy & Data → Download Your Data")
        print("3. Download your raw genetic data")
        print("4. Place the .txt file in this directory")
        return False

def create_info_file():
    """Create an info file about the analysis."""
    info_content = """
GenoInsight Pro - Advanced Genetic Analysis Tool
================================================

WHAT THIS TOOL DOES:
- Analyzes your 23andMe raw genetic data using cutting-edge scientific research
- Calculates polygenic risk scores for complex diseases
- Provides pharmacogenomic insights for drug metabolism
- Screens for rare pathogenic variants
- Generates publication-quality visualizations
- Creates a comprehensive 15-20 page scientific report

KEY FEATURES:
1. Disease Risk Analysis - Based on latest GWAS studies (2023-2025)
2. Polygenic Risk Scores - For heart disease, diabetes, Alzheimer's, etc.
3. Pharmacogenomics - How you metabolize common medications
4. Rare Variant Screening - Check for known pathogenic mutations
5. Ancestry Analysis - Using multiple ancestry-informative markers
6. Trait Analysis - Athletic performance, caffeine metabolism, etc.

OUTPUT FILES:
- Comprehensive text report with all findings and citations
- Multiple visualization plots in the genetic_analysis_plots/ folder
- Console output with key findings highlighted

IMPORTANT NOTES:
- This is for educational and research purposes only
- Not a substitute for medical advice or clinical genetic testing
- Always consult healthcare professionals for medical decisions
- Consider genetic counseling for significant findings

HOW TO RUN:
python genoinsight_pro.py

TYPICAL ANALYSIS TIME: 2-5 minutes

For questions or issues, ensure all packages are installed:
pip install pandas numpy matplotlib seaborn scipy scikit-learn networkx requests
"""
    
    with open('GenoInsight_Pro_README.txt', 'w') as f:
        f.write(info_content)
    
    print("\n📄 Created GenoInsight_Pro_README.txt with detailed information")

def main():
    """Main setup and launch function."""
    print("="*60)
    print("GenoInsight Pro - Setup and Launch Assistant")
    print("="*60)
    
    # Check Python version
    if not check_python_version():
        sys.exit(1)
    
    # Check and install packages
    if not install_requirements():
        print("\n⚠️  Please install missing packages before running the analysis")
        sys.exit(1)
    
    # Check for data files
    data_available = check_data_file()
    
    # Create info file
    create_info_file()
    
    # Check if main script exists
    if not os.path.exists('genoinsight_pro.py'):
        print("\n❌ genoinsight_pro.py not found!")
        print("   Make sure you've saved the main analysis script as 'genoinsight_pro.py'")
        sys.exit(1)
    
    print("\n" + "="*60)
    print("SETUP COMPLETE!")
    print("="*60)
    
    if data_available:
        print("\n✅ Everything is ready for analysis!")
        response = input("\nWould you like to run the analysis now? (y/n): ")
        
        if response.lower() == 'y':
            print("\nLaunching GenoInsight Pro...")
            print("-"*60)
            subprocess.call([sys.executable, 'genoinsight_pro.py'])
    else:
        print("\n⚠️  Please add your 23andMe data file before running the analysis")
        print("   Once you have your data file, run: python genoinsight_pro.py")

if __name__ == "__main__":
    main()