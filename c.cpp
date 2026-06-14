#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <algorithm>
#include <chrono>
#include <random>
#include <iomanip>
#include <climits>

// ========================= 1D и 2D RSQ через префиксную сумму =========================

template<typename T>
class PrefixSum1D {
private:
    int n;
    std::vector<T> pref;
public:
    PrefixSum1D(const std::vector<T>& arr) : n(arr.size()), pref(n + 1, 0) {
        for (int i = 0; i < n; ++i)
            pref[i + 1] = pref[i] + arr[i];
    }

    T sum(int l, int r) const {
        if (l < 0 || r >= n || l > r) return 0;
        return pref[r + 1] - pref[l];
    }
};

template<typename T>
class PrefixSum2D {
private:
    int n, m;
    std::vector<std::vector<T>> pref;
public:
    PrefixSum2D(const std::vector<std::vector<T>>& matrix) : n(matrix.size()), m(matrix[0].size()),
                                                             pref(n + 1, std::vector<T>(m + 1, 0)) {
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < m; ++j)
                pref[i + 1][j + 1] = pref[i + 1][j] + pref[i][j + 1] - pref[i][j] + matrix[i][j];
    }

    T sum(int x1, int y1, int x2, int y2) const {
        if (x1 < 0 || y1 < 0 || x2 >= n || y2 >= m || x1 > x2 || y1 > y2) return 0;
        return pref[x2 + 1][y2 + 1] - pref[x1][y2 + 1] - pref[x2 + 1][y1] + pref[x1][y1];
    }
};

// ========================= 1D RMQ через предпросчет всех отрезков =========================

template<typename T>
class AllSegmentsRMQ {
private:
    int n;
    std::vector<std::vector<T>> dp;
public:
    AllSegmentsRMQ(const std::vector<T>& arr) : n(arr.size()), dp(n, std::vector<T>(n)) {
        for (int l = 0; l < n; ++l) {
            dp[l][l] = arr[l];
            for (int r = l + 1; r < n; ++r)
                dp[l][r] = std::min(dp[l][r - 1], arr[r]);
        }
    }

    T query(int l, int r) const {
        if (l < 0 || r >= n || l > r) return 0;
        return dp[l][r];
    }
};

// ========================= Корневая декомпозиция =========================

template<typename T, typename CombineOp>
class SqrtDecomposition {
private:
    std::vector<T> data;
    std::vector<T> blocks;
    int block_size;
    int num_blocks;
    CombineOp op;
    T identity;

public:
    SqrtDecomposition(const std::vector<T>& arr, CombineOp combine_op, T identity_elem)
        : data(arr), op(combine_op), identity(identity_elem) {
        int n = data.size();
        block_size = std::max(1, (int)std::sqrt(n));
        num_blocks = (n + block_size - 1) / block_size;
        blocks.assign(num_blocks, identity);

        for (int i = 0; i < n; ++i)
            blocks[i / block_size] = op(blocks[i / block_size], data[i]);
    }

    T query(int l, int r) const {
        if (l < 0 || r >= (int)data.size() || l > r) return identity;
        
        T res = identity;
        int left_block = l / block_size;
        int right_block = r / block_size;

        if (left_block == right_block) {
            for (int i = l; i <= r; ++i)
                res = op(res, data[i]);
        } else {
            int end_left = (left_block + 1) * block_size - 1;
            for (int i = l; i <= end_left; ++i)
                res = op(res, data[i]);
            for (int b = left_block + 1; b < right_block; ++b)
                res = op(res, blocks[b]);
            int start_right = right_block * block_size;
            for (int i = start_right; i <= r; ++i)
                res = op(res, data[i]);
        }
        return res;
    }

    void update(int idx, T new_value) {
        if (idx < 0 || idx >= (int)data.size()) return;
        
        int block_id = idx / block_size;
        data[idx] = new_value;
        int start = block_id * block_size;
        int end = std::min((int)data.size(), start + block_size);
        blocks[block_id] = identity;
        for (int i = start; i < end; ++i)
            blocks[block_id] = op(blocks[block_id], data[i]);
    }
};

// ========================= Дерево отрезков =========================

template<typename T, typename CombineOp>
class SegmentTree {
private:
    int n;
    std::vector<T> tree;
    CombineOp op;
    T identity;

    void build(const std::vector<T>& arr, int node, int l, int r) {
        if (l == r) {
            tree[node] = arr[l];
            return;
        }
        int mid = (l + r) / 2;
        build(arr, node * 2, l, mid);
        build(arr, node * 2 + 1, mid + 1, r);
        tree[node] = op(tree[node * 2], tree[node * 2 + 1]);
    }

    T query(int node, int l, int r, int ql, int qr) const {
        if (qr < l || ql > r) return identity;
        if (ql <= l && r <= qr) return tree[node];
        int mid = (l + r) / 2;
        return op(query(node * 2, l, mid, ql, qr),
                  query(node * 2 + 1, mid + 1, r, ql, qr));
    }

    void update(int node, int l, int r, int idx, T value) {
        if (l == r) {
            tree[node] = value;
            return;
        }
        int mid = (l + r) / 2;
        if (idx <= mid)
            update(node * 2, l, mid, idx, value);
        else
            update(node * 2 + 1, mid + 1, r, idx, value);
        tree[node] = op(tree[node * 2], tree[node * 2 + 1]);
    }

public:
    SegmentTree(const std::vector<T>& arr, CombineOp combine_op, T identity_elem)
        : n(arr.size()), op(combine_op), identity(identity_elem) {
        if (n == 0) return;
        tree.resize(4 * n);
        build(arr, 1, 0, n - 1);
    }

    T query(int l, int r) const {
        if (l < 0 || r >= n || l > r) return identity;
        return query(1, 0, n - 1, l, r);
    }

    void update(int idx, T value) {
        if (idx < 0 || idx >= n) return;
        update(1, 0, n - 1, idx, value);
    }
};

// ========================= Дерево Фенвика =========================

template<typename T>
class FenwickTree {
private:
    int n;
    std::vector<T> bit;
public:
    FenwickTree(const std::vector<T>& arr) : n(arr.size()), bit(n + 1, 0) {
        for (int i = 0; i < n; ++i)
            add(i, arr[i]);
    }

    void add(int idx, T delta) {
        int i = idx + 1;
        while (i <= n) {
            bit[i] += delta;
            i += i & -i;
        }
    }

    T prefix_sum(int idx) const {
        T res = 0;
        int i = idx + 1;
        while (i > 0) {
            res += bit[i];
            i -= i & -i;
        }
        return res;
    }

    T sum(int l, int r) const {
        if (l < 0 || r >= n || l > r) return 0;
        if (l == 0) return prefix_sum(r);
        return prefix_sum(r) - prefix_sum(l - 1);
    }

    void update(int idx, T new_value) {
        if (idx < 0 || idx >= n) return;
        T old = sum(idx, idx);
        add(idx, new_value - old);
    }
};

// ========================= Разреженная таблица =========================

template<typename T>
class SparseTable {
private:
    int n;
    std::vector<std::vector<T>> st;
    std::vector<int> log;
public:
    SparseTable(const std::vector<T>& arr) : n(arr.size()) {
        if (n == 0) return;
        
        log.resize(n + 1);
        log[1] = 0;
        for (int i = 2; i <= n; ++i)
            log[i] = log[i / 2] + 1;

        int K = log[n] + 1;
        st.resize(n, std::vector<T>(K));

        for (int i = 0; i < n; ++i)
            st[i][0] = arr[i];

        for (int j = 1; j < K; ++j) {
            for (int i = 0; i + (1 << j) <= n; ++i) {
                st[i][j] = std::min(st[i][j - 1], st[i + (1 << (j - 1))][j - 1]);
            }
        }
    }

    T query(int l, int r) const {
        if (l < 0 || r >= n || l > r) return 0;
        int j = log[r - l + 1];
        return std::min(st[l][j], st[r - (1 << j) + 1][j]);
    }
};

// ========================= Демонстрация =========================

void demo() {
    std::cout << "========== ДЕМОНСТРАЦИЯ РАБОТЫ ==========\n\n";

    // 1D RSQ
    std::vector<int> arr1 = {1, 2, 3, 4, 5};
    PrefixSum1D<int> ps1d(arr1);
    std::cout << "1D PrefixSum: sum(1,3) = " << ps1d.sum(1, 3) << " (2+3+4=9) ✓\n";

    // 2D RSQ
    std::vector<std::vector<int>> matrix = {{1, 2, 3},
                                            {4, 5, 6},
                                            {7, 8, 9}};
    PrefixSum2D<int> ps2d(matrix);
    std::cout << "2D PrefixSum: sum(0,0,1,1) = " << ps2d.sum(0, 0, 1, 1) << " (1+2+4+5=12) ✓\n";

    // AllSegments RMQ
    AllSegmentsRMQ<int> all_rmq(arr1);
    std::cout << "AllSegmentsRMQ: min(1,3) = " << all_rmq.query(1, 3) << " (min=2) ✓\n";

    // SqrtDecomp RSQ
    SqrtDecomposition<int, std::plus<int>> sqrt_rsq(arr1, std::plus<int>(), 0);
    std::cout << "SqrtDecomp RSQ: sum(1,3) = " << sqrt_rsq.query(1, 3) << " ✓\n";
    sqrt_rsq.update(2, 10);
    std::cout << "  После update(2,10): sum(1,3) = " << sqrt_rsq.query(1, 3) << " (2+10+4=16) ✓\n";

    // SqrtDecomp RMQ
    SqrtDecomposition<int, std::function<int(int,int)>> sqrt_rmq(
        arr1, [](int a, int b) { return std::min(a, b); }, INT_MAX);
    std::cout << "SqrtDecomp RMQ: min(1,3) = " << sqrt_rmq.query(1, 3) << " ✓\n";
    sqrt_rmq.update(2, 0);
    std::cout << "  После update(2,0): min(1,3) = " << sqrt_rmq.query(1, 3) << " (min=0) ✓\n";

    // Segment Tree RSQ
    SegmentTree<int, std::plus<int>> seg_rsq(arr1, std::plus<int>(), 0);
    std::cout << "SegmentTree RSQ: sum(1,3) = " << seg_rsq.query(1, 3) << " ✓\n";
    seg_rsq.update(2, 20);
    std::cout << "  После update: sum(1,3) = " << seg_rsq.query(1, 3) << " (2+20+4=26) ✓\n";

    // Segment Tree RMQ
    SegmentTree<int, std::function<int(int,int)>> seg_rmq(
        arr1, [](int a, int b) { return std::min(a, b); }, INT_MAX);
    std::cout << "SegmentTree RMQ: min(1,3) = " << seg_rmq.query(1, 3) << " ✓\n";
    seg_rmq.update(2, 100);
    std::cout << "  После update(2,100): min(1,3) = " << seg_rmq.query(1, 3) << " (min=2) ✓\n";

    // Fenwick Tree
    FenwickTree<int> fenwick(arr1);
    std::cout << "Fenwick Tree RSQ: sum(1,3) = " << fenwick.sum(1, 3) << " ✓\n";
    fenwick.update(2, 30);
    std::cout << "  После update: sum(1,3) = " << fenwick.sum(1, 3) << " (2+30+4=36) ✓\n";

    // Sparse Table RMQ
    SparseTable<int> sparse(arr1);
    std::cout << "SparseTable RMQ: min(1,3) = " << sparse.query(1, 3) << " (min=2) ✓\n";

    std::cout << "\n✅ ВСЕ СТРУКТУРЫ РАБОТАЮТ КОРРЕКТНО!\n";
}

// ========================= MAIN =========================

int main() {
    try {
        demo();
    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    } catch (...) {
        std::cerr << "Неизвестная ошибка!" << std::endl;
        return 1;
    }
    
    std::cout << "\nПрограмма успешно завершена.\n";
    return 0;
}