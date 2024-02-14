#pragma once

#include <iostream>
#include <vector>
#include <memory>
#include <list>

#include <array>
#include <cmath>

template <size_t N>
class StackStorage {
 public:
  std::array<char, N> data_;
  size_t last = 0;
  StackStorage() = default;
  StackStorage(const StackStorage& other) = delete;
  StackStorage& operator=(const StackStorage& other) = delete;
};

template <typename T, size_t N>
class StackAllocator {
 public:
  StackStorage<N>* stack_storage_ = nullptr;
  using value_type = T;
  StackAllocator()
      : stack_storage_(nullptr) {
  }
  StackAllocator(StackStorage<N>& stack_storage)
      : stack_storage_(&stack_storage) {
  }
  template <typename U>
  StackAllocator(const StackAllocator<U, N>& other)
      : stack_storage_(other.stack_storage_) {
  }
  StackAllocator(StackAllocator<T, N>&& other)
      : stack_storage_(std::move(other.stack_storage_)) {
  }
  template <typename U>
  StackAllocator<T, N>& operator=(const StackAllocator<U, N>& other) {
    stack_storage_ = other.stack_storage_;
    return *this;
  }

  T* allocate(size_t num) {
    stack_storage_->last +=
        ((sizeof(T) - (stack_storage_->last % sizeof(T))) % sizeof(T));
    auto to_ret = stack_storage_->data_.data() + stack_storage_->last;
    stack_storage_->last += num * sizeof(T);
    return reinterpret_cast<T*>(to_ret);
  }

  void deallocate(T*, size_t) {
  }

  template <typename... Args>
  void construct(T* ptr, Args&... args) {
    new (ptr) T(args...);
  }

  void destroy(T* ptr) {
    ptr->~T();
  }

  template <typename U>
  struct rebind {
    using other = StackAllocator<U, N>;
  };

  ~StackAllocator() = default;
};

template <typename T, typename Alloc = std::allocator<T>>
class List {
 private:
  struct BaseNode {
    BaseNode* next_ = nullptr;
    BaseNode* prev_ = nullptr;
    BaseNode() = default;
    BaseNode(BaseNode* next_node, BaseNode* prev_node)
        : next_(next_node),
          prev_(prev_node) {
    }
    ~BaseNode() = default;
  };

  struct Node : public BaseNode {
    T value_;
    Node() = default;
    explicit Node(const T& value)
        : value_(value) {
    }
    explicit Node(T&& value)
        : value_(std::move(value)) {
    }
    ~Node() = default;
  };

  size_t size_;
  BaseNode* fake_node_;

  using AllocTraits = std::allocator_traits<Alloc>;
  [[no_unique_address]] typename AllocTraits::template rebind_alloc<Node>
      node_alloc_ = Alloc();
  [[no_unique_address]] typename AllocTraits::template rebind_alloc<BaseNode>
      basenode_alloc_ = Alloc();
  using node_alloc_type_ = typename AllocTraits::template rebind_alloc<Node>;
  using basenode_alloc_type_ =
      typename AllocTraits::template rebind_alloc<BaseNode>;

  void default_push_back() {
    Node* new_node =
        std::allocator_traits<node_alloc_type_>::allocate(node_alloc_, 1);
    try {
      std::allocator_traits<node_alloc_type_>::construct(node_alloc_, new_node);
    } catch (...) {
      std::allocator_traits<node_alloc_type_>::deallocate(node_alloc_, new_node,
                                                          1);
      throw;
    }

    auto before = end().base_node_ptr_->prev_;

    end().base_node_ptr_->prev_ = new_node;
    new_node->prev_ = before;
    new_node->next_ = end().base_node_ptr_;
    before->next_ = new_node;
    ++size_;
  }

 public:
  template <bool constant>
  class ListIterator {
   private:
    BaseNode* base_node_ptr_;

   public:
    using pointer = std::conditional_t<constant, const T*, T*>;
    using reference = std::conditional_t<constant, const T&, T&>;
    using value_type = std::conditional_t<constant, const T, T>;
    using iterator_category = std::bidirectional_iterator_tag;
    using difference_type = std::ptrdiff_t;

    ListIterator() = default;
    ListIterator(const ListIterator& other) = default;
    explicit ListIterator(BaseNode* base_node_ptr)
        : base_node_ptr_(base_node_ptr) {
    }

    reference operator*() const {
      return static_cast<Node*>(base_node_ptr_)->value_;
    }

    pointer operator->() const {
      return &(this->operator*());
    }

    ListIterator& operator++() {
      base_node_ptr_ = base_node_ptr_->next_;
      return *this;
    }

    ListIterator operator++(int) {
      ListIterator to_ret = *this;
      base_node_ptr_ = base_node_ptr_->next_;
      return to_ret;
    }

    ListIterator& operator--() {
      base_node_ptr_ = base_node_ptr_->prev_;
      return *this;
    }

    ListIterator operator--(int) {
      ListIterator to_ret = *this;
      base_node_ptr_ = base_node_ptr_->prev_;
      return to_ret;
    }

    operator ListIterator<true>() const {
      return ListIterator<true>(base_node_ptr_);
    }

    bool operator==(const ListIterator& other) const {
      return base_node_ptr_ == other.base_node_ptr_;
    }

    bool operator!=(const ListIterator& other) const {
      return base_node_ptr_ != other.base_node_ptr_;
    }

    friend ListIterator<true> List::insert(ListIterator<true> iter,
                                           const T& val);
    friend ListIterator<true> List::erase(ListIterator<true> iter);
    friend void List::default_push_back();
    friend ListIterator<true> List::list_insert_helper(ListIterator<true> iter,
                                                       Node* new_node);
    ~ListIterator() = default;
  };

  using iterator = ListIterator<false>;
  using const_iterator = ListIterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<ListIterator<true>>;
  iterator begin() {
    return iterator(fake_node_->next_);
  }

  iterator end() {
    return iterator(fake_node_);
  }

  const_iterator cbegin() {
    return const_iterator(fake_node_->next_);
  }

  const_iterator cend() {
    return const_iterator(fake_node_);
  }

  const_iterator begin() const {
    return const_iterator(fake_node_->next_);
  }

  const_iterator end() const {
    return const_iterator(fake_node_);
  }

  std::reverse_iterator<iterator> rbegin() {
    return std::reverse_iterator<iterator>(end());
  }

  std::reverse_iterator<iterator> rend() {
    return std::reverse_iterator<iterator>(begin());
  }

  std::reverse_iterator<const_iterator> rbegin() const {
    return std::reverse_iterator<const_iterator>(end());
  }

  std::reverse_iterator<const_iterator> rend() const {
    return std::reverse_iterator<const_iterator>(begin());
  }

  std::reverse_iterator<const_iterator> crbegin() const {
    return std::reverse_iterator<const_iterator>(cend());
  }

  std::reverse_iterator<const_iterator> crend() const {
    return std::reverse_iterator<const_iterator>(cbegin());
  }

  void clear() {
    while (!empty()) {
      pop_back();
    }
    std::allocator_traits<basenode_alloc_type_>::destroy(basenode_alloc_,
                                                         fake_node_);
    std::allocator_traits<basenode_alloc_type_>::deallocate(basenode_alloc_,
                                                            fake_node_, 1);
  }

  bool empty() const {
    return size_ == 0;
  }

  void push_back(T&& val) {
    try {
      insert(end(), std::forward(val));
    } catch (...) {
      throw;
    }
  }

  void push_front(T&& val) {
    try {
      insert(begin(), std::forward(val));
    } catch (...) {
      throw;
    }
  }

  void pop_back() {
    try {
      erase(--end());
    } catch (...) {
      throw;
    }
  }

  void pop_front() {
    try {
      erase(begin());
    } catch (...) {
      throw;
    }
  }

  size_t size() const {
    return size_;
  }
  node_alloc_type_ get_allocator() const {
    return node_alloc_;
  }

  List()
      : size_(0),
        fake_node_(std::allocator_traits<basenode_alloc_type_>::allocate(
            basenode_alloc_, 1)) {
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);
  }

  explicit List(size_t num)
      : size_(0),
        fake_node_(std::allocator_traits<basenode_alloc_type_>::allocate(
            basenode_alloc_, 1)) {
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);

    try {
      while (num > 0) {
        default_push_back();
        --num;
      }
    } catch (...) {
      clear();
      throw;
    }
  }

  List(size_t num, const T& val)
      : size_(0),
        fake_node_(std::allocator_traits<basenode_alloc_type_>::allocate(
            basenode_alloc_, 1)) {
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);

    try {
      while (num > 0) {
        push_back(val);
        --num;
      }
    } catch (...) {
      clear();
      throw;
    }
  }

 private:
  List get_correct_list(const node_alloc_type_& other_alloc) {
    if (std::allocator_traits<
            Alloc>::propagate_on_container_copy_assignment::value) {
      return List(other_alloc);
    } else {
      return List(node_alloc_);
    }
  }
  ListIterator<true> list_insert_helper(ListIterator<true> iter,
                                        Node* new_node) {
    auto before = iter.base_node_ptr_->prev_;
    iter.base_node_ptr_->prev_ = new_node;
    new_node->prev_ = before;
    new_node->next_ = iter.base_node_ptr_;
    before->next_ = new_node;
    ++size_;
    return ListIterator<true>(new_node--);
  }

 public:
  List& operator=(const List& other) {
    if (&other == this) {
      return *this;
    }
    List to_copy(get_correct_list(other.node_alloc_));

    for (auto& it : other) {
      to_copy.push_back(it);
    }

    std::swap(fake_node_, to_copy.fake_node_);
    std::swap(basenode_alloc_, to_copy.basenode_alloc_);
    std::swap(node_alloc_, to_copy.node_alloc_);
    std::swap(size_, to_copy.size_);
    return *this;
  }

  List& operator=(List&& other) noexcept {
    if (&other == this) {
      return *this;
    }

    fake_node_ = other.fake_node_;
    size_ = other.size_;
    node_alloc_ = other.node_alloc_;
    other.size_ = 0;
    other.fake_node_ = std::allocator_traits<basenode_alloc_type_>::allocate(
        other.basenode_alloc_, 1);
    other.fake_node_->prev_ = other.fake_node_;
    other.fake_node_->next_ = other.fake_node_;
    return *this;
  }

  List(const List& other)
      : size_(0) {
    node_alloc_ =
        (std::allocator_traits<Alloc>::select_on_container_copy_construction(
            other.node_alloc_));
    basenode_alloc_ =
        (std::allocator_traits<Alloc>::select_on_container_copy_construction(
            other.basenode_alloc_));

    fake_node_ = std::allocator_traits<basenode_alloc_type_>::allocate(
        basenode_alloc_, 1);
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);
    try {
      for (auto& it : other) {
        push_back(it);
      }
    } catch (...) {
      clear();
      throw;
    }
  }

  List(List&& other) noexcept
      : size_(other.size_),
        fake_node_(other.fake_node_),
        node_alloc_(other.node_alloc_) {
    other.size_ = 0;
    other.fake_node_ = std::allocator_traits<basenode_alloc_type_>::allocate(
        other.basenode_alloc_, 1);
    other.fake_node_->prev_ = other.fake_node_;
    other.fake_node_->next_ = other.fake_node_;
  }

  List(size_t num, const Alloc& alloc)
      : size_(0),
        node_alloc_(alloc),
        basenode_alloc_(alloc) {
    fake_node_ = std::allocator_traits<basenode_alloc_type_>::allocate(
        basenode_alloc_, 1);
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);
    try {
      while (num > 0) {
        default_push_back();
        --num;
      }
    } catch (...) {
      clear();
      throw;
    }
  }

  List(size_t num, const T& val, const Alloc& alloc)
      : size_(0),
        node_alloc_(alloc),
        basenode_alloc_(alloc) {
    fake_node_ = std::allocator_traits<basenode_alloc_type_>::allocate(
        basenode_alloc_, 1);
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);
    try {
      while (num > 0) {
        push_back(val);
        --num;
      }
    } catch (...) {
      clear();
      throw;
    }
  }

  List(const Alloc& alloc)
      : size_(0),
        node_alloc_(alloc),
        basenode_alloc_(alloc) {
    fake_node_ = std::allocator_traits<basenode_alloc_type_>::allocate(
        basenode_alloc_, 1);
    std::allocator_traits<basenode_alloc_type_>::construct(
        basenode_alloc_, fake_node_, fake_node_, fake_node_);
  }

  ListIterator<true> insert(ListIterator<true> iter, T&& val) {
    Node* new_node =
        std::allocator_traits<node_alloc_type_>::allocate(node_alloc_, 1);
    try {
      std::allocator_traits<node_alloc_type_>::construct(node_alloc_, new_node,
                                                         std::move(val));
    } catch (...) {
      std::allocator_traits<node_alloc_type_>::deallocate(node_alloc_, new_node,
                                                          1);
      throw;
    }
    return list_insert_helper(iter, new_node);
  }

  ListIterator<true> insert(ListIterator<true> iter, const T& val) {
    Node* new_node =
        std::allocator_traits<node_alloc_type_>::allocate(node_alloc_, 1);
    try {
      std::allocator_traits<node_alloc_type_>::construct(node_alloc_, new_node,
                                                         val);
    } catch (...) {
      std::allocator_traits<node_alloc_type_>::deallocate(node_alloc_, new_node,
                                                          1);
      throw;
    }
    return list_insert_helper(iter, new_node);
  }

  ListIterator<true> erase(ListIterator<true> iter) {
    BaseNode* cur_node = iter.base_node_ptr_;
    BaseNode* prev_node = cur_node->prev_;
    BaseNode* next_node = cur_node->next_;
    prev_node->next_ = next_node;
    next_node->prev_ = prev_node;
    std::allocator_traits<node_alloc_type_>::destroy(
        node_alloc_, static_cast<Node*>(cur_node));
    std::allocator_traits<node_alloc_type_>::deallocate(
        node_alloc_, static_cast<Node*>(cur_node), 1);
    --size_;
    return ListIterator<true>(next_node);
  }

  ~List() {
    clear();
  }
};

template <typename Key, typename Value, typename Hash = std::hash<Key>,
          typename Equal = std::equal_to<Key>,
          typename Alloc = std::allocator<std::pair<const Key, Value>>>
class UnorderedMap {
 private:
  using NodeType = std::pair<const Key, Value>;
  using NodeTypeAllocatorType =
      typename std::allocator_traits<Alloc>::template rebind_alloc<NodeType*>;
  using ListIterator =
      typename List<NodeType*, NodeTypeAllocatorType>::const_iterator;
  using ListIteratorAllocatorType = typename std::allocator_traits<
      Alloc>::template rebind_alloc<ListIterator>;

  Hash hashgen_{};
  Equal equal_to_{};
  Alloc allocator_{};
  List<NodeType*, NodeTypeAllocatorType> global_list_;
  std::vector<ListIterator, ListIteratorAllocatorType> buckets_;
  std::vector<ListIterator, ListIteratorAllocatorType> bucket_ends_;
  size_t size_;
  double max_load_factor_ = 0.6;

 public:
  UnorderedMap()
      : size_(0) {
    buckets_.resize(1, global_list_.end());
    bucket_ends_.resize(1, global_list_.end());
  }
  UnorderedMap(const UnorderedMap& other)
      : hashgen_(other.hashgen_),
        equal_to_(other.equal_to_),
        allocator_(other.allocator_),
        size_(0),
        max_load_factor_(other.max_load_factor_) {
    buckets_.resize(1, global_list_.end());
    bucket_ends_.resize(1, global_list_.end());
    for (auto& val : other) {
      insert(val);
    }
  }

  UnorderedMap(UnorderedMap&& other) noexcept
      : hashgen_(std::move(other.hashgen_)),
        equal_to_(std::move(other.equal_to_)),
        allocator_(other.allocator_),
        global_list_(std::move(other.global_list_)),
        buckets_(std::move(other.buckets_)),
        bucket_ends_(std::move(other.bucket_ends_)),
        size_(other.size_),
        max_load_factor_(other.max_load_factor_) {
    other.size_ = 0;
  }

  ~UnorderedMap() = default;

  template <bool constant>
  class MapIterator {
   private:
    ListIterator list_iterator_;

   public:
    using pointer = std::conditional_t<constant, const NodeType*, NodeType*>;
    using reference = std::conditional_t<constant, const NodeType&, NodeType&>;
    using value_type = std::conditional_t<constant, const NodeType, NodeType>;
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;

    MapIterator() = default;
    MapIterator(const MapIterator& other) = default;
    MapIterator(const ListIterator& list_iterator)
        : list_iterator_(list_iterator) {
    }

    ListIterator get_list_iter() const {
      return list_iterator_;
    }

    reference operator*() const {
      return *(*list_iterator_);
    }

    pointer operator->() const {
      return *list_iterator_;
    }

    MapIterator& operator++() {
      ++list_iterator_;
      return *this;
    }

    MapIterator operator++(int) {
      MapIterator to_ret = *this;
      ++list_iterator_;
      return to_ret;
    }

    MapIterator& operator--() {
      --list_iterator_;
      return *this;
    }

    ListIterator operator--(int) {
      MapIterator to_ret = *this;
      --list_iterator_;
      return to_ret;
    }

    operator MapIterator<true>() const {
      return MapIterator<true>(list_iterator_);
    }

    bool operator!=(const MapIterator& other) const {
      return list_iterator_ != other.list_iterator_;
    }

    bool operator==(const MapIterator& other) const {
      return list_iterator_ == other.list_iterator_;
    }

    ~MapIterator() = default;
  };

  using iterator = MapIterator<false>;
  using const_iterator = MapIterator<true>;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<MapIterator<true>>;

  iterator begin() {
    return iterator(global_list_.begin());
  }

  iterator end() {
    return iterator(global_list_.end());
  }

  const_iterator cbegin() {
    return const_iterator(global_list_.begin());
  }

  const_iterator cend() {
    return const_iterator(global_list_.end());
  }

  const_iterator begin() const {
    return const_iterator(global_list_.begin());
  }

  const_iterator end() const {
    return const_iterator(global_list_.end());
  }

  std::reverse_iterator<iterator> rbegin() {
    return std::reverse_iterator<iterator>(end());
  }

  std::reverse_iterator<iterator> rend() {
    return std::reverse_iterator<iterator>(begin());
  }

  std::reverse_iterator<const_iterator> rbegin() const {
    return std::reverse_iterator<const_iterator>(end());
  }

  std::reverse_iterator<const_iterator> rend() const {
    return std::reverse_iterator<const_iterator>(begin());
  }

  std::reverse_iterator<const_iterator> crbegin() const {
    return std::reverse_iterator<const_iterator>(cend());
  }

  std::reverse_iterator<const_iterator> crend() const {
    return std::reverse_iterator<const_iterator>(cbegin());
  }

 private:
  size_t get_bucket(const Key& key) const {
    return hashgen_(key) % buckets_.size();
  }

  ListIterator find_helper(const Key& key) const {
    size_t current_bucket = get_bucket(key);
    if (buckets_[current_bucket] == global_list_.end()) {
      return global_list_.end();
    }
    for (ListIterator current_node = buckets_[current_bucket];
         current_node != std::next(bucket_ends_[current_bucket]);
         current_node = std::next(current_node)) {
      auto& [iter_key, val] = *(*current_node);
      if (equal_to_(iter_key, key)) {
        return current_node;
      }
    }
    return global_list_.end();
  }

 public:
  size_t size() const {
    return size_;
  }

  UnorderedMap& operator=(const UnorderedMap& other) {
    UnorderedMap temp(other);
    *this = std::move(temp);
    return *this;
  }

  UnorderedMap& operator=(UnorderedMap&& other) noexcept {
    if (this == &other) {
      return *this;
    }
    size_ = other.size_;
    max_load_factor_ = other.max_load_factor_;
    allocator_ = other.allocator_;
    global_list_.clear();
    global_list_ = std::move(other.global_list_);
    hashgen_ = std::move(other.hashgen_);
    equal_to_ = std::move(other.equal_to_);
    buckets_ = std::move(other.buckets_);
    bucket_ends_ = std::move(other.bucket_ends_);
    return *this;
  }

  iterator find(const Key& key) {
    return find_helper(key);
  }

  const_iterator find(const Key& key) const {
    return find_helper(key);
  }

  void max_load_factor(double new_max_load_factor) {
    max_load_factor_ = new_max_load_factor;
  }

  double max_load_factor() const {
    return max_load_factor_;
  }

  double load_factor() const {
    return 1.0 * size() / buckets_.size();
  }

  std::pair<iterator, bool> insert(const NodeType& node) {
    return emplace(node);
  }

  template <typename InputIterator>
  void insert(InputIterator begin, InputIterator end) {
    for (auto it = begin; it != end; it = std::next(it)) {
      insert(*it);
    }
  }

  template <typename T>
  std::pair<iterator, bool> insert(T&& val) {
    return emplace(std::forward<T>(val));
  }

  void rehash(size_t new_size) {
    auto old_global_list = std::move(global_list_);
    buckets_.resize(new_size);
    buckets_.assign(new_size, global_list_.end());
    bucket_ends_.resize(new_size);
    bucket_ends_.assign(new_size, global_list_.end());
    auto current_end = global_list_.end();
    size_t inserted_num = 0;
    try {
      for (auto& it : old_global_list) {
        auto& [key, val] = *it;
        auto current_bucket = get_bucket(key);
        if (buckets_[current_bucket] == current_end) {
          bucket_ends_[current_bucket] = buckets_[current_bucket] =
              global_list_.insert(global_list_.end(), it);
        } else {
          bucket_ends_[current_bucket] =
              global_list_.insert(std::next(buckets_[current_bucket]), it);
        }
        ++inserted_num;
      }
    } catch (...) {
      for (int i = 0; i < inserted_num; ++i) {
        global_list_.erase(global_list_.begin());
      }
      throw;
    }
  }

  template <typename... Args>
  std::pair<iterator, bool> emplace(Args&&... args) {
    if (max_load_factor_ * size_ > static_cast<double>(buckets_.size())) {
      rehash(static_cast<size_t>(
          ceil((2 * (size_ + 1)) * 1.0 / max_load_factor_)));
    }
    auto new_node = std::allocator_traits<Alloc>::allocate(allocator_, 1);
    std::allocator_traits<Alloc>::construct(allocator_, new_node,
                                            std::forward<Args>(args)...);
    auto& [key, val] = *new_node;
    iterator finded = find(key);
    if (finded != end()) {
      return {finded, false};
    }

    size_t current_bucket = get_bucket(key);
    if (buckets_[current_bucket] == global_list_.end()) {
      bucket_ends_[current_bucket] = buckets_[current_bucket] =
          global_list_.insert(global_list_.end(), new_node);
    } else {
      bucket_ends_[current_bucket] = global_list_.insert(
          std::next(bucket_ends_[current_bucket]), new_node);
    }
    ++size_;
    return {iterator(bucket_ends_[current_bucket]), true};
  }

  iterator erase(const_iterator to_erase) {
    ListIterator current_list_iterator = to_erase.get_list_iter();
    ListIterator to_ret = std::next(current_list_iterator);
    auto& [key, value] = *to_erase;
    size_t current_bucket = get_bucket(key);

    auto& [first_elem_key, first_elem_val] = *(*buckets_[current_bucket]);
    auto& [last_elem_key, last_elem_val] = *(*bucket_ends_[current_bucket]);
    if (equal_to_(key, first_elem_key)) {
      if (bucket_ends_[current_bucket] == buckets_[current_bucket]) {
        buckets_[current_bucket] = bucket_ends_[current_bucket] =
            global_list_.end();
      } else {
        buckets_[current_bucket] = std::next(buckets_[current_bucket]);
      }
    } else if (equal_to_(key, last_elem_key)) {
      if (bucket_ends_[current_bucket] == buckets_[current_bucket]) {
        buckets_[current_bucket] = bucket_ends_[current_bucket] =
            global_list_.end();
      } else {
        bucket_ends_[current_bucket] = std::prev(bucket_ends_[current_bucket]);
      }
    }
    global_list_.erase(current_list_iterator);
    --size_;
    return to_ret;
  }

  template <typename EraseIterator>
  const_iterator erase(EraseIterator begin, EraseIterator end) {
    for (auto cur_iter = begin; cur_iter != end; cur_iter = erase(cur_iter)) {
      if (cur_iter == std::prev(end)) {
        return erase(cur_iter);
      }
    }
    return cend();
  }

  Value& operator[](const Key& key) {
    iterator finded = find(key);
    if (finded != end()) {
      return finded->second;
    }
    if constexpr (std::is_default_constructible<Value>()) {
      auto [iter, is_inserted] = emplace(key, Value());
      return iter->second;
    }
    // and if we here we're returning some bullshit,
    // cause it's a user problem
    return finded->second;
  }

  Value& at(const Key& key) {
    iterator finded = find(key);
    if (finded == end()) {
      throw std::out_of_range("loshara");
    }
    return finded->second;
  }

  void reserve(size_t count) {
    if (size() >= count) {
      return;
    }
    rehash(static_cast<size_t>(ceil(count * 1.0 / max_load_factor_)));
  }
};
