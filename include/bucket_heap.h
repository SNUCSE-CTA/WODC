#ifndef BUCKET_HEAP_H_
#define BUCKET_HEAP_H_

class BucketHeap {
 public:
  BucketHeap(ui capacity, ui depth);
  ~BucketHeap();

  void push(ui u, ui d);
  void pop(ui u);
  bool is_popped(ui u);

  void inc(ui u);
  void dec(ui u);

  void inc_all();
  void dec_all();

  inline ui cap(ui d);
  inline ui cap(ui d, ui cnt);
  inline ui get(ui d, ui i);
  inline ui lab(ui u);
  inline ui getpos(ui u);

  void make_capacity(ui d, ui n);

  ui offset;
  ui capacity, depth;

 private:
  struct Node {
    ui label, pos;
  }* nodes;
  ui* ns;
  ui* sub;
};

BucketHeap::BucketHeap(ui capacity, ui depth)
    : capacity(capacity), depth(depth) {
  nodes = new Node[capacity];
  ns = new ui[depth * 2 + 1]();
  sub = new ui[capacity * (depth * 2 + 1)];

  for (ui i = 0; i < capacity; ++i) {
    nodes[i].label = 0;
    nodes[i].pos = INVALID;
  }
  offset = depth;
}

BucketHeap::~BucketHeap() {
  delete[] nodes;
  delete[] ns;
  delete[] sub;
}

void BucketHeap::push(ui u, ui d) {
  assert(is_popped(u));
  assert(offset + d < depth * 2 + 1);
  this->nodes[u].pos = (offset + d) * capacity + this->ns[offset + d]++;
  this->sub[this->nodes[u].pos] = u;
  this->nodes[u].label = d - (depth - offset);
  assert(lab(u) == d);
}

void BucketHeap::pop(ui u) {
  assert(!is_popped(u));
  const ui ud = lab(u);
  const ui up = this->nodes[u].pos;
  this->nodes[u].pos = INVALID;
  assert(up != INVALID);
  assert(this->sub[up] == u);

  const ui v = this->sub[(offset + ud) * capacity + --this->ns[offset + ud]];
  if (u != v) {
    assert(lab(v) == ud);
    this->nodes[v].pos = up;
    this->sub[up] = v;
  }
}

void BucketHeap::inc(ui u) {
  const ui ud = lab(u);
  const ui up = this->nodes[u].pos;
  const ui v = this->sub[(offset + ud) * capacity + --this->ns[offset + ud]];
  if (u != v) {
    this->nodes[v].pos = up;
    this->sub[up] = v;
  }

  this->nodes[u].pos =
      (offset + ud + 1) * capacity + this->ns[offset + ud + 1]++;
  this->sub[this->nodes[u].pos] = u;
  this->nodes[u].label = offset + ud + 1 - depth;
}

void BucketHeap::dec(ui u) {
  const ui ud = lab(u);
  const ui up = this->nodes[u].pos;

  const ui v = this->sub[(offset + ud) * capacity + --this->ns[offset + ud]];
  if (u != v) {
    this->nodes[v].pos = up;
    this->sub[up] = v;
  }

  this->nodes[u].pos =
      (offset + ud - 1) * capacity + this->ns[offset + ud - 1]++;
  this->sub[this->nodes[u].pos] = u;
  this->nodes[u].label = ud - 1 - (depth - offset);
}

void BucketHeap::inc_all() { offset -= 1; }

void BucketHeap::dec_all() { offset += 1; }

bool BucketHeap::is_popped(ui u) { return this->nodes[u].pos == INVALID; }

inline ui BucketHeap::cap(ui d) { return this->ns[offset + d]; }

inline ui BucketHeap::cap(ui d, ui cnt) {
  ui res = 0;

  for (ui i = d + 1; i-- > 0;) {
    res += this->ns[offset + i];
    if (cnt-- == 0) {
      break;
    }
  }
  return res;
}

inline ui BucketHeap::get(ui d, ui i) {
  return this->sub[(offset + d) * capacity + i];
}

inline ui BucketHeap::lab(ui u) {
  return this->nodes[u].label + depth - offset;
}

inline ui BucketHeap::getpos(ui u) { return this->nodes[u].pos; }

void BucketHeap::make_capacity(ui d, ui n) { this->ns[offset + d] = n; }

#endif  // BUCKET_HEAP_H_
