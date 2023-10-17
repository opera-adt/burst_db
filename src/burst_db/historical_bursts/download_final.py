from pathlib import Path

from tqdm.contrib.concurrent import thread_map

from burst_db.historical_bursts.download import download_safe_metadata

batch_size = 50
products = Path("../missing_products.txt").read_text().splitlines()
print(f"{len(products) = }")
remaining_products = [p for p in products if not Path(p + ".SAFE").exists()]
print(f"{len(remaining_products) = }")
batches = [
    remaining_products[i : i + batch_size]
    for i in range(0, len(remaining_products), batch_size)
]

thread_map(download_safe_metadata, batches, max_workers=5)
