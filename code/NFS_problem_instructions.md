# NFS failure

## Context Behind the NFS Failure

On compute nodes, the shared working directory `/scratch/lema` was not a direct mount to a local disk, but a virtual network mount running over **NFSv4**.

This setup suffered from two critical configuration vulnerabilities:

* **Constrained I/O buffers (`wsize=8192`):** Every generated file is fragmented into 8 KB packets.
* **Strict blocking semantics (`hard`):** The NFS client retries network requests indefinitely without allowing user-level interrupts.

Running parallel worker processes that simultaneously wrote Parquet files, competed to create identical subdirectory trees, and redirected standard logging into a single file on `$HOME` (also an NFS mount) overwhelmed the filesystem metadata layer and clogged the server's RPC queue. The Linux kernel transitioned these processes into **`D (disk sleep)`**—an uninterruptible sleep state where the operating system ignores all external signals (including `kill -9`). This deadlock stalled the entire filesystem, prevented resource cleanup, and ultimately rendered the node unresponsive.

---

## Purpose of the Diagnostic Routine

The audit routine was developed to prevent similar failures on other cluster nodes by meeting three operational objectives:

* **Unmask mount point abstractions:** Determine whether directories like `/scratch` map to genuine local disks (`ext4`/`xfs`) or network mounts (`autofs`/`nfs4`) that introduce latency.
* **Pinpoint optimal physical storage:** Locate bare partitions on local solid-state hardware (such as `/export/scratch` backed by NVMe drives), allowing workloads to write at full PCIe bus speeds without touching the operating system's network stack.
* **Perform write permission validation:** Query disk capacity and verify write permissions.

echo "=== 1. Physical disks ==="
lsblk -e7 -o NAME,ROTA,TYPE,TRAN,SIZE,FSTYPE,MOUNTPOINT

echo -e "\n=== 2. File system in scratch routes ==="
for dir in /tmp /scratch /scratch/$(hostname) /export/scratch; do
  if [ -e "$dir" ]; then
    findmnt -T "$dir" -o TARGET,SOURCE,FSTYPE,OPTIONS
  fi
done

echo -e "\n=== 3. Free space (no NFS) ==="
df -h -x nfs -x nfs4 -x cifs

echo -e "\n=== 4. Writing test in local scratch ==="
for candidate in "/export/scratch/tmp/$USER" "/tmp/$USER"; do
  if mkdir -p "$candidate" 2>/dev/null && touch "$candidate/.test" 2>/dev/null; then
    rm -f "$candidate/.test"
    echo "Ruta local escribible confirmada: $candidate"
    break
  fi
done