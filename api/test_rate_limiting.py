#!/usr/bin/env python3
"""
Test rate limiting functionality of the Shedding Hub API.
"""

import requests
import time
import concurrent.futures
from typing import List

API_BASE_URL = "http://localhost:8000"

def make_request(endpoint: str) -> tuple:
    """Make a single request and return status code and response time."""
    start_time = time.time()
    try:
        response = requests.get(f"{API_BASE_URL}{endpoint}", timeout=5)
        end_time = time.time()
        return response.status_code, end_time - start_time, None
    except Exception as e:
        end_time = time.time()
        return None, end_time - start_time, str(e)

def test_rate_limiting(endpoint: str, num_requests: int = 20, max_workers: int = 10):
    """Test rate limiting on a specific endpoint."""
    print(f"\n[TEST] Testing rate limiting on {endpoint}")
    print(f"Making {num_requests} concurrent requests with {max_workers} workers...")
    
    start_total = time.time()
    
    # Make concurrent requests
    with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(make_request, endpoint) for _ in range(num_requests)]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    
    end_total = time.time()
    
    # Analyze results
    success_count = sum(1 for status, _, _ in results if status == 200)
    rate_limited_count = sum(1 for status, _, _ in results if status == 429)
    error_count = sum(1 for status, _, _ in results if status is None)
    other_errors = sum(1 for status, _, _ in results if status not in [200, 429, None])
    
    avg_response_time = sum(response_time for _, response_time, _ in results) / len(results)
    
    print(f"[RESULT] Total time: {end_total - start_total:.2f}s")
    print(f"[RESULT] Successful requests: {success_count}")
    print(f"[RESULT] Rate limited requests (429): {rate_limited_count}")
    print(f"[RESULT] Network errors: {error_count}")
    print(f"[RESULT] Other errors: {other_errors}")
    print(f"[RESULT] Average response time: {avg_response_time:.3f}s")
    
    # Print some example rate limited responses
    rate_limited_responses = [(s, e) for s, _, e in results if s == 429]
    if rate_limited_responses:
        print(f"[INFO] Rate limiting is working! Got {len(rate_limited_responses)} 429 responses")
    else:
        print(f"[WARNING] No rate limiting detected. All requests succeeded.")
    
    return success_count, rate_limited_count, error_count

def test_different_endpoints():
    """Test rate limiting on different types of endpoints."""
    print("=" * 60)
    print("[TEST] Rate Limiting Test Suite")
    print("=" * 60)
    
    # Test default rate limit (100/minute)
    print(f"\n[INFO] Testing default endpoints (100/minute limit)")
    test_rate_limiting("/", 15, 5)
    test_rate_limiting("/health", 15, 5)
    
    # Small delay to avoid cross-endpoint interference
    time.sleep(2)
    
    # Test search endpoint (30/minute limit)
    print(f"\n[INFO] Testing search endpoint (30/minute limit)")
    test_rate_limiting("/search?q=test", 20, 8)
    
    # Small delay
    time.sleep(2)
    
    # Test expensive endpoint (10/minute limit)
    print(f"\n[INFO] Testing expensive endpoint (10/minute limit)")
    # Get a dataset ID first
    response = requests.get(f"{API_BASE_URL}/datasets?limit=1")
    if response.status_code == 200:
        datasets = response.json().get('datasets', [])
        if datasets:
            dataset_id = datasets[0]['dataset_id']
            test_rate_limiting(f"/datasets/{dataset_id}", 15, 8)
        else:
            print("[WARNING] No datasets found for testing expensive endpoint")
    else:
        print("[ERROR] Could not get dataset list for testing")

def test_rate_limit_recovery():
    """Test that rate limits recover after the time window."""
    print(f"\n[INFO] Testing rate limit recovery")
    
    endpoint = "/health"
    
    # Make enough requests to trigger rate limiting
    print("Making initial burst of requests...")
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
        futures = [executor.submit(make_request, endpoint) for _ in range(20)]
        results = [future.result() for future in concurrent.futures.as_completed(futures)]
    
    rate_limited = sum(1 for status, _, _ in results if status == 429)
    print(f"Rate limited requests: {rate_limited}")
    
    if rate_limited > 0:
        print("Waiting 10 seconds for rate limit to reset...")
        time.sleep(10)
        
        # Try a single request
        status, _, _ = make_request(endpoint)
        if status == 200:
            print("[SUCCESS] Rate limit recovered successfully!")
        else:
            print(f"[WARNING] Rate limit may not have recovered. Status: {status}")
    else:
        print("[INFO] No rate limiting triggered in initial burst")

def main():
    """Run all rate limiting tests."""
    # Check if API is running
    try:
        response = requests.get(API_BASE_URL, timeout=5)
        if response.status_code != 200:
            print(f"[ERROR] API not responding properly. Status: {response.status_code}")
            return
    except requests.exceptions.RequestException as e:
        print(f"[ERROR] Cannot connect to API: {e}")
        print("Make sure the API server is running at http://localhost:8000")
        return
    
    print("[OK] API server is responding")
    
    # Run tests
    test_different_endpoints()
    test_rate_limit_recovery()
    
    print("\n" + "=" * 60)
    print("[DONE] Rate limiting tests completed")
    print("=" * 60)

if __name__ == "__main__":
    main()