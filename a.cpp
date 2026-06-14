#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <random>
#include <iomanip>
#include <climits>

long long op_count = 0;
void reset_counter() { op_count = 0; }

class PrefixSum1D {
private:
    std::vector<int> pref;
public:
    PrefixSum1D(const std::vector<int>& arr) {
        int n = arr.size();
        pref.assign(n + 1, 0);
        for (int i = 0; i < n; i++) {
            pref[i + 1] = pref[i] + arr[i];
        }
    }
    
    int sum(int l, int r) {
        op_count += 2;
        return pref[r + 1] - pref[l];
    }
};

class PrefixSum2D {
private:
    std::vector<std::vector<int>> pref;
    int n, m;
public:
    PrefixSum2D(const std::vector<std::vector<int>>& matrix) {
        n = matrix.size();
        m = matrix[0].size();
        pref.assign(n + 1, std::vector<int>(m + 1, 0));
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                pref[i + 1][j + 1] = pref[i + 1][j] + pref[i][j + 1] - pref[i][j] + matrix[i][j];
            }
        }
    }
    
    int sum(int x1, int y1, int x2, int y2) {
        op_count += 4;
        return pref[x2 + 1][y2 + 1] - pref[x1][y2 + 1] - pref[x2 + 1][y1] + pref[x1][y1];
    }
};

class AllSegmentsRMQ {
private:
    std::vector<std::vector<int>> dp;
public:
    AllSegmentsRMQ(const std::vector<int>& arr) {
        int n = arr.size();
        dp.assign(n, std::vector<int>(n));
        for (int l = 0; l < n; l++) {
            dp[l][l] = arr[l];
            for (int r = l + 1; r < n; r++) {
                dp[l][r] = std::min(dp[l][r-1], arr[r]);
            }
        }
    }
    
    int query(int l, int r) {
        op_count++;
        return dp[l][r];
    }
};

class SqrtDecompRSQ {
private:
    std::vector<int> data;
    std::vector<int> blocks;
    int block_size;
public:
    SqrtDecompRSQ(const std::vector<int>& arr) : data(arr) {
        int n = arr.size();
        block_size = std::max(1, (int)sqrt(n));
        int num_blocks = (n + block_size - 1) / block_size;
        blocks.assign(num_blocks, 0);
        
        for (int i = 0; i < n; i++) {
            blocks[i / block_size] += data[i];
        }
    }
    
    int query(int l, int r) {
        int res = 0;
        int left_block = l / block_size;
        int right_block = r / block_size;
        
        if (left_block == right_block) {
            for (int i = l; i <= r; i++) {
                res += data[i];
                op_count++;
            }
        } else {
            int end_left = (left_block + 1) * block_size - 1;
            for (int i = l; i <= end_left; i++) {
                res += data[i];
                op_count++;
            }
            for (int b = left_block + 1; b < right_block; b++) {
                res += blocks[b];
                op_count++;
            }
            int start_right = right_block * block_size;
            for (int i = start_right; i <= r; i++) {
                res += data[i];
                op_count++;
            }
        }
        return res;
    }
};

class SqrtDecompRMQ {
private:
    std::vector<int> data;
    std::vector<int> blocks;
    int block_size;
public:
    SqrtDecompRMQ(const std::vector<int>& arr) : data(arr) {
        int n = arr.size();
        block_size = std::max(1, (int)sqrt(n));
        int num_blocks = (n + block_size - 1) / block_size;
        blocks.assign(num_blocks, INT_MAX);
        
        for (int i = 0; i < n; i++) {
            blocks[i / block_size] = std::min(blocks[i / block_size], data[i]);
        }
    }
    
    int query(int l, int r) {
        int res = INT_MAX;
        int left_block = l / block_size;
        int right_block = r / block_size;
        
        if (left_block == right_block) {
            for (int i = l; i <= r; i++) {
                res = std::min(res, data[i]);
                op_count++;
            }
        } else {
            int end_left = (left_block + 1) * block_size - 1;
            for (int i = l; i <= end_left; i++) {
                res = std::min(res, data[i]);
                op_count++;
            }
            for (int b = left_block + 1; b < right_block; b++) {
                res = std::min(res, blocks[b]);
                op_count++;
            }
            int start_right = right_block * block_size;
            for (int i = start_right; i <= r; i++) {
                res = std::min(res, data[i]);
                op_count++;
            }
        }
        return res;
    }
};

class SegmentTreeRSQ {
private:
    std::vector<int> tree;
    int n;
    
    void build(const std::vector<int>& arr, int node, int l, int r) {
        if (l == r) {
            tree[node] = arr[l];
            return;
        }
        int mid = (l + r) / 2;
        build(arr, node*2, l, mid);
        build(arr, node*2+1, mid+1, r);
        tree[node] = tree[node*2] + tree[node*2+1];
    }
    
    int query(int node, int l, int r, int ql, int qr) {
        if (qr < l || ql > r) return 0;
        if (ql <= l && r <= qr) return tree[node];
        int mid = (l + r) / 2;
        int left_res = query(node*2, l, mid, ql, qr);
        int right_res = query(node*2+1, mid+1, r, ql, qr);
        op_count++;
        return left_res + right_res;
    }
    
public:
    SegmentTreeRSQ(const std::vector<int>& arr) {
        n = arr.size();
        tree.resize(4*n);
        build(arr, 1, 0, n-1);
    }
    
    int query(int l, int r) {
        return query(1, 0, n-1, l, r);
    }
};

class SegmentTreeRMQ {
private:
    std::vector<int> tree;
    int n;
    
    void build(const std::vector<int>& arr, int node, int l, int r) {
        if (l == r) {
            tree[node] = arr[l];
            return;
        }
        int mid = (l + r) / 2;
        build(arr, node*2, l, mid);
        build(arr, node*2+1, mid+1, r);
        tree[node] = std::min(tree[node*2], tree[node*2+1]);
    }
    
    int query(int node, int l, int r, int ql, int qr) {
        if (qr < l || ql > r) return INT_MAX;
        if (ql <= l && r <= qr) return tree[node];
        int mid = (l + r) / 2;
        int left_res = query(node*2, l, mid, ql, qr);
        int right_res = query(node*2+1, mid+1, r, ql, qr);
        op_count++;
        return std::min(left_res, right_res);
    }
    
public:
    SegmentTreeRMQ(const std::vector<int>& arr) {
        n = arr.size();
        tree.resize(4*n);
        build(arr, 1, 0, n-1);
    }
    
    int query(int l, int r) {
        return query(1, 0, n-1, l, r);
    }
};

class FenwickTree {
private:
    std::vector<int> bit;
    int n;
    
    int prefix_sum(int idx) {
        int res = 0;
        int i = idx + 1;
        while (i > 0) {
            res += bit[i];
            op_count++;
            i -= i & -i;
        }
        return res;
    }
    
public:
    FenwickTree(const std::vector<int>& arr) {
        n = arr.size();
        bit.assign(n + 1, 0);
        for (int i = 0; i < n; i++) {
            int idx = i + 1;
            while (idx <= n) {
                bit[idx] += arr[i];
                idx += idx & -idx;
            }
        }
    }
    
    int sum(int l, int r) {
        if (l == 0) {
            return prefix_sum(r);
        }
        int res = prefix_sum(r) - prefix_sum(l - 1);
        op_count++;
        return res;
    }
};

class SparseTable {
private:
    std::vector<std::vector<int>> st;
    std::vector<int> log;
    
public:
    SparseTable(const std::vector<int>& arr) {
        int n = arr.size();
        log.resize(n+1);
        log[1] = 0;
        for (int i = 2; i <= n; i++)
            log[i] = log[i/2] + 1;
        
        int k = log[n] + 1;
        st.assign(n, std::vector<int>(k));
        
        for (int i = 0; i < n; i++)
            st[i][0] = arr[i];
            
        for (int j = 1; j < k; j++) {
            for (int i = 0; i + (1<<j) <= n; i++) {
                st[i][j] = std::min(st[i][j-1], st[i + (1<<(j-1))][j-1]);
            }
        }
    }
    
    int query(int l, int r) {
        int j = log[r - l + 1];
        op_count++;
        return std::min(st[l][j], st[r - (1<<j) + 1][j]);
    }
};

void test_structure(const std::string& name, 
    std::function<void(const std::vector<int>&, int, int)> test_func,
    const std::vector<int>& sizes) {
    std::cout << "\n" << name << "\n";
    std::cout << "n\tкол-во опер.\n";
    std::random_device rd;
    std::mt19937 gen(42);
    std::uniform_int_distribution<int> val_dist(1, 100);
    
    for (int n : sizes) {
        std::vector<int> arr(n);
        for (int i = 0; i < n; i++)
            arr[i] = val_dist(gen);
        
        reset_counter();
        test_func(arr, n, 1000);
        double avg_ops = (double)op_count / 1000.0;
        std::cout << n << "\t" << std::fixed << std::setprecision(2) << avg_ops << "\n";
    }
}

void test_2d_prefix(const std::vector<int>& arr, int n, int queries) {
    std::vector<std::vector<int>> matrix(n, std::vector<int>(n));
    std::random_device rd;
    std::mt19937 gen(123);
    std::uniform_int_distribution<int> val_dist(1, 100);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = val_dist(gen);
    
    PrefixSum2D ps2d(matrix);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int x1 = idx_dist(gen);
        int y1 = idx_dist(gen);
        int x2 = idx_dist(gen);
        int y2 = idx_dist(gen);
        if (x1 > x2) std::swap(x1, x2);
        if (y1 > y2) std::swap(y1, y2);
        ps2d.sum(x1, y1, x2, y2);
    }
}

void test_prefix_sum(const std::vector<int>& arr, int n, int queries) {
    PrefixSum1D ps(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        ps.sum(l, r);
    }
}

void test_all_segments(const std::vector<int>& arr, int n, int queries) {
    if (n > 500) return;
    AllSegmentsRMQ rmq(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        rmq.query(l, r);
    }
}

void test_sqrt_rsq(const std::vector<int>& arr, int n, int queries) {
    SqrtDecompRSQ sqrt(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        sqrt.query(l, r);
    }
}

void test_sqrt_rmq(const std::vector<int>& arr, int n, int queries) {
    SqrtDecompRMQ sqrt(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        sqrt.query(l, r);
    }
}

void test_seg_rsq(const std::vector<int>& arr, int n, int queries) {
    SegmentTreeRSQ seg(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        seg.query(l, r);
    }
}

void test_seg_rmq(const std::vector<int>& arr, int n, int queries) {
    SegmentTreeRMQ seg(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        seg.query(l, r);
    }
}

void test_fenwick(const std::vector<int>& arr, int n, int queries) {
    FenwickTree fw(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        fw.sum(l, r);
    }
}

void test_sparse(const std::vector<int>& arr, int n, int queries) {
    SparseTable st(arr);
    std::uniform_int_distribution<int> idx_dist(0, n-1);
    std::mt19937 gen(123);
    
    reset_counter();
    
    for (int q = 0; q < queries; q++) {
        int l = idx_dist(gen);
        int r = idx_dist(gen);
        if (l > r) std::swap(l, r);
        st.query(l, r);
    }
}

int main() {
    std::vector<int> sizes = {100, 500, 1000, 5000, 10000};
    std::vector<int> small_sizes = {50, 100, 200, 300, 500};
    std::vector<int> matrix_sizes = {10, 20, 50, 100, 200};
    
    test_structure("1) 1D RSQ O(1)", test_prefix_sum, sizes);
    test_structure("2) 2D RSQ O(1)", test_2d_prefix, matrix_sizes);
    test_structure("3) All Segments RMQ O(1)", test_all_segments, small_sizes);
    test_structure("4) Sqrt Decomp RSQ O(√n)", test_sqrt_rsq, sizes);
    test_structure("5) Sqrt Decomp RMQ O(√n)", test_sqrt_rmq, sizes);
    test_structure("6) Segment Tree RSQ O(log n)", test_seg_rsq, sizes);
    test_structure("7) Segment Tree RMQ O(log n)", test_seg_rmq, sizes);
    test_structure("8) Fenwick Tree RSQ O(log n)", test_fenwick, sizes);
    test_structure("9) Sparse Table RMQ O(1)", test_sparse, sizes);
    
    return 0;
}