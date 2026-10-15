# Algorithm

Overburden selects templates in two passes over the input BAM file. Pass 1 chooses the
templates, and pass 2 writes them. The problem is stated in [PLAN.md](PLAN.md).

Template selection on spliced, multi-segment templates is a set multicover problem
(NP-hard). For contiguous single-segment reads, the sweep below is exactly optimal
(exchange argument). For everything else it is a greedy heuristic: it guarantees the depth
target, but not minimality.

## Pass 1: sweep (selection)

Each reference sequence is processed left to right. Positions are evaluated only at
events: block starts, and `end + 1` of blocks. Both depths are constant between events,
and selection only ever increases the selected depth, so checking at events is enough.

### State

* `selected::Set{UInt64}`: hashes of selected QNAMEs. The only state that grows with the
  file size (about 8 bytes per selected template).
* `pending`: min-heap of blocks not yet reached, keyed by start. A spliced record's later
  blocks start after the record's start, so they wait here.
* `active`: the blocks covering the current position, as a bucket queue keyed by end. It is
  a ring of linked lists, where bucket `end & (W - 1)` holds the blocks ending at `end`, for
  a ring size `W` that is a power of two. All active ends lie in `p:(p + W - 1)`, so each
  bucket holds blocks of a single end. When a block longer than the ring is activated,
  the ring is doubled until the block fits. A lower bound `lo` of the active ends is kept.
  The next end is found by scanning forward from `lo` to the first non-empty bucket. The
  number of active blocks is the full depth `d(p)`.
* `selected_depth`: the number of active blocks whose template is selected. This is a
  counter, not a second heap. An expiring block is selected iff its template is not a key
  of `unselected`, so expiry needs no second structure.
* `unselected::Dict{UInt64, Int}`: for each unselected template with active blocks,
  the number of its active blocks. Selecting a template moves all of its active blocks into
  `selected_depth` in one step.
* `candidates`: max-heap of `(end, template)` entries for the active blocks of unselected
  templates. It is maintained lazily: entries for ended blocks or for selected templates
  are discarded when popped. It is compacted when it holds more than `2n + 1024` entries,
  where `n` is the number of active blocks of unselected templates. Without compaction,
  stale entries pile up in regions where no selection is needed. Compaction is an in-place
  loop rather than `filter!`, because `filter!` shrinks the vector's memory, which then has
  to be reallocated as the heap grows again.

Records with QNAME `*` all hash to one template, since they have no name to tell them apart.

### Events

At each event position `p`:

1. Remove blocks from `active` that ended before `p`, by emptying the buckets from `lo`
   to `p - 1`. Decrement `selected_depth` for selected templates. Otherwise, decrement the
   template's count in `unselected`, and remove the template when the count reaches zero.
2. Activate pending blocks starting at or before `p`. If the block's template is selected,
   it counts towards `selected_depth`. Otherwise it is added to `unselected` and
   `candidates`.
3. Let `target = min(D, d(p))`. While `selected_depth < target`, pop the candidate with the
   furthest end. Discard it if it is stale. Otherwise add its template to `selected`, and
   move all of that template's active blocks into `selected_depth`.
4. Compact `candidates` if needed.

Every active block of an unselected template has an entry in `candidates`. So while
`selected_depth < target`, the heap cannot run empty.

### Reference boundaries

At the end of a reference, all blocks have ended. Every entry left in `candidates` is then
stale, so it must be discarded before the next reference starts. Otherwise an entry from
the previous reference, with a larger position than the current one, could win the
furthest-end comparison against blocks of the new reference. This would select a template
for a worse reason than the rule above. The same template name can occur on several
references (e.g. the two mates of a pair on different chromosomes), so this happens in
practice.

### Correctness

* The target depth is met at every position `p`. Selection at an event covers the
  positions up to the next event, because the depth only changes at events.
* Segments of a selected template that have not been reached yet are counted as selected
  once the sweep reaches them. They can't contribute before their start, so nothing is
  missed.
* Segments passed before their template was selected are still output in pass 2. This
  only adds surplus coverage.
* QNAMEs are tracked as 64-bit hashes. A collision merges two templates into one selection
  unit. This only costs optimality, never correctness, since both passes use the same hash.
* The hash is only compared within one process, so it does not need to be stable across
  Julia versions.

### Optimality

For contiguous single-segment reads on one reference, the rule "when short of depth,
select the candidate that extends the furthest" is exactly optimal, by an exchange argument:
any optimal solution can be changed to include the furthest-reaching candidate without
losing coverage. Tests check this against brute force (see below).

For spliced or multi-segment templates, one template covers several blocks, and choosing a
template can help several positions at once. The rule is then a heuristic.

### Complexity

Each block causes one push and one pop on `pending`, one push on `candidates` if
unselected, and constant-time insertion into and removal from `active`. So the heap work is
`O(B log B)` for `B` blocks.

Finding the next end in `active` also scans empty buckets. `lo` only moves back when a
new block ends before it, so typical data scans a few positions per event. In the worst
case, a long block among sparse short ones, the stretch up to the long block's end is
scanned again after each short block ends. Memory is `O(max depth + longest block)` plus
the selected set.

### Pipeline

Pass 1 runs in two concurrent tasks. The calling task reads records and checks their sort
order and CIGARs, computing their blocks and QNAME hashes. A spawned task runs the sweep.
The reader reuses one record (`copy = false`), so the blocks are copied into batches of
4096 records, as parallel vectors. Records without blocks are left out, since they don't
affect the sweep. Batches go to the sweep through one channel, and return for reuse through
another. Four batches are in use, which bounds memory, and `put!` never blocks.

All errors from bad input come from the reader, so they are thrown on the calling task
unwrapped. On an error, the reader closes the channel and waits for the sweep task to end.
If the sweep task fails, the channel of free batches is closed (`bind`), so the reader
doesn't wait forever.

On CDIFFICILE (D=20, 4 threads), each stage alone takes about 0.57 s, and pass 1 takes
0.69 s, down from 1.09 s with one task. Batch sizes from 1024 to 16384 perform the same.

## Pass 2: output

Seek back to the first record (`XAM2.seek_first_record`), stream all records again, and
write each record whose QNAME hash is in `selected`. Records are copied as raw bytes, without
re-encoding. The output header is the input header, with an added `@PG` line. Records keep
their input order.

## Testing

* `test/sweep.jl` checks the depth target against a per-base brute-force depth. It checks
  that no template is output partially (no singletons), and that templates covering nothing
  are never selected.
* It checks exact optimality on tiny contiguous single-segment instances against a brute
  force search.
* It runs the sweep directly on blocks up to 20 kb long spread over 300 kb, so the ring of
  `active` grows and wraps around many times, with gaps where nothing is active. It checks
  the depth target, and checks that the number of selected blocks matches the greedy rule
  evaluated naively at every event.
* It checks that the selection doesn't depend on the batch size or the number of batches,
  and that an error in a later batch is thrown unwrapped.
* `test/output.jl` checks that pass 2 writes exactly the records of the selected templates.

## Possible later improvements

* Tie-breaking and scoring: prefer templates whose mate lies ahead of `p` on the same
  contig (`PNEXT > p`), or higher MAPQ.
* A pruning pass that removes templates which turned out to be redundant.
* Report a lower bound (e.g. per-contig `max_p t(p)`), to estimate the optimality gap.
* Internal coordinate sorting, and stdin input via spooling to a temporary file.
* A bitmap of non-empty buckets in `active`, so the scan for the next end checks 64
  positions at a time. This only matters in the worst case of the scan described above.
* Pending heap: each record's blocks are already sorted, so they could be buffered and
  merged in order, instead of one heap push and pop per block.
* Faster lookups: map each template to its first-occurrence index and store a `BitVector`
  of kept templates. This makes the pass-2 lookup nearly free, at the cost of memory in
  pass 1 (see the TODO in `src/sweep.jl`).
