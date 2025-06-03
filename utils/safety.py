import functools
import json
import traceback
import datetime as dt
import pathlib as pl

def safeguard(stage: str):
    """
    A decorator to catch exceptions in major analysis stages,
    log them, and dump partial results for post-mortem analysis.
    """
    def wrap(fn):
        @functools.wraps(fn)
        def inner(analyzer_instance, *args, **kwargs):
            try:
                return fn(analyzer_instance, *args, **kwargs)
            except Exception as exc:  # noqa: BLE001 (broad exception catch is intended here)
                # Ensure the crash_dumps directory exists
                crash_dump_dir = pl.Path("crash_dumps")
                crash_dump_dir.mkdir(exist_ok=True)
                
                # Prepare data for dumping
                # Accessing 'self.results' from the passed 'analyzer_instance'
                results_to_dump = {}
                if hasattr(analyzer_instance, 'results') and isinstance(analyzer_instance.results, dict):
                    # Attempt to serialize, converting non-serializable items to string
                    try:
                        # A more robust serialization might be needed for complex objects
                        results_to_dump = json.loads(json.dumps(analyzer_instance.results, default=str))
                    except TypeError:
                        results_to_dump = {"error": "Could not serialize all results."}
                else:
                    results_to_dump = {"error": "Analyzer instance has no 'results' attribute or it's not a dict."}

                dump_data = {
                    "stage_failed": stage,
                    "timestamp_utc": dt.datetime.utcnow().isoformat(),
                    "error_type": type(exc).__name__,
                    "error_message": str(exc),
                    "traceback": traceback.format_exc(),
                    "partial_results_dumped": results_to_dump,
                }
                
                # Define dump file path
                dump_file_path = crash_dump_dir / f"crash_{stage}_{int(dt.datetime.utcnow().timestamp())}.json"
                
                # Write dump to file
                try:
                    dump_file_path.write_text(json.dumps(dump_data, indent=2))
                    print(f"\nCRITICAL ERROR in stage '{stage}'. Traceback logged and partial results dumped to: {dump_file_path}")
                except Exception as dump_exc:
                    print(f"\nCRITICAL ERROR in stage '{stage}'. Failed to write dump file: {dump_exc}")
                
                # Re-raise the original exception to halt further execution if desired,
                # or handle it (e.g., allow other stages to run)
                raise  # Or: return None / some error indicator
        return inner
    return wrap
