#!/usr/bin/env python3
"""
download_db.py: Download transposase profile HMM reference database for pyISFinder
"""

import os
import sys
import urllib.request

HMM_URL = "https://raw.githubusercontent.com/xiezhq/ISEScan/master/pHMMs/clusters.faa.hmm"
REF_URL = "https://raw.githubusercontent.com/xiezhq/ISEScan/master/pHMMs/clusters.single.faa"

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    db_dir = os.path.join(script_dir, "db")
    os.makedirs(db_dir, exist_ok=True)
    
    hmm_path = os.path.join(db_dir, "clusters.faa.hmm")
    ref_path = os.path.join(db_dir, "clusters.single.faa")
    
    print("=" * 60)
    print("        pyISFinder - Reference Database Setup")
    print("=" * 60)
    print(f"Target Database Directory: {db_dir}\n")
    
    # Download HMM database
    if os.path.exists(hmm_path) and os.path.getsize(hmm_path) > 1000:
        print(f"[✓] Transposase profile HMM database already exists at: {hmm_path}")
    else:
        print("[↓] Downloading transposase profile HMM database (clusters.faa.hmm, ~40MB)...")
        try:
            urllib.request.urlretrieve(HMM_URL, hmm_path)
            print(f"[✓] Successfully downloaded: {hmm_path}")
        except Exception as e:
            print(f"[✗] Error downloading HMM database: {e}")
            sys.exit(1)
            
    # Download single reference FASTA
    if os.path.exists(ref_path) and os.path.getsize(ref_path) > 1000:
        print(f"[✓] Reference FASTA database already exists at: {ref_path}")
    else:
        print("[↓] Downloading reference transposase FASTA (clusters.single.faa, ~100KB)...")
        try:
            urllib.request.urlretrieve(REF_URL, ref_path)
            print(f"[✓] Successfully downloaded: {ref_path}")
        except Exception as e:
            print(f"[✗] Error downloading reference FASTA: {e}")
            sys.exit(1)
            
    print("\n" + "=" * 60)
    print("Database Setup Complete! You can now run pyISFinder.")
    print("=" * 60)

if __name__ == "__main__":
    main()
