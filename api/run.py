#!/usr/bin/env python3
"""
Convenience script to run the Shedding Hub API
"""

import subprocess
import sys
import os
from pathlib import Path

def main():
    """Run the API with proper setup."""
    
    # Get the directory containing this script
    api_dir = Path(__file__).parent
    project_root = api_dir.parent
    
    print("Starting Shedding Hub API...")
    
    # Check if we're in the right directory
    if not (project_root / "data").exists():
        print("Error: Cannot find data directory. Make sure you're running this from the shedding-hub project root.")
        sys.exit(1)
    
    # Check dependencies
    try:
        import fastapi
        import uvicorn
        print("Dependencies found")
    except ImportError:
        print("📦 Installing API dependencies...")
        subprocess.run([sys.executable, "-m", "pip", "install", "-r", "requirements.txt"], 
                      cwd=api_dir, check=True)
        print("Dependencies installed")
    
    # Start the server
    print("Starting server on http://localhost:8004")
    print("API documentation will be available at http://localhost:8004/docs")
    print("Use Ctrl+C to stop the server")
    print("")
    
    try:
        import uvicorn
        import os
        
        # Change to the API directory so uvicorn can find main.py
        os.chdir(api_dir)
        
        # Run uvicorn with the module string for reload support
        uvicorn.run("main:app", host="0.0.0.0", port=8004, reload=True)
        
    except KeyboardInterrupt:
        print("\nShutting down API server...")
    except Exception as e:
        print(f"Error starting server: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
