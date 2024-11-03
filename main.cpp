#include <algorithm>
#include <cstdint>
#include <iostream>
#include <utility>
#include <vector>
#include <unordered_set>

const int p_a = 2147079149;
const int p_b = 2147249981;

int BinPow(int a_, int64_t b_, int ring) {
    if (b_ == 0) {
        return 1;
    }
    int64_t result = static_cast<int64_t>(BinPow(a_, b_ / 2, ring));
    int64_t result_square = result * result;
    if (b_ % 2 == 0) {
        int result_final = static_cast<int>(result_square % ring);
        return result_final;
    } else {
        int64_t result_mod = a_ * (result_square % ring);
        int result_final = static_cast<int>(result_mod % ring);
        return result_final;
    }
}
 
int Inv(int a_, int field) {
    return BinPow(a_, field - 2, field);
}

struct Hash {
    int first;
    int second;
    bool operator==(const Hash& hash) const {
        return (first == hash.first and second == hash.second);
    }
    bool operator!=(const Hash& hash) const {
        return !(*this == hash);
    }
};

int N_inp;
int M_inp;
std::vector<std::vector<Hash>> prefix_hash;
std::vector<std::vector<int>> table;
std::vector<int> a_inv_coeff;
std::vector<int> b_inv_coeff;

Hash GetHash(int i_one, int j_one, int width) {
    int i_two = i_one + width - 1;
    int j_two = j_one + width - 1;

    int result_a = prefix_hash[i_two][j_two].first;
    int result_b = prefix_hash[i_two][j_two].second;

    if (j_one > 0) {
        result_a -= prefix_hash[i_two][j_one - 1].first;
        if (result_a < 0) {
            result_a = p_a + result_a;
        }
        result_b -= prefix_hash[i_two][j_one - 1].second;
        if (result_b < 0) {
            result_b = p_b + result_b;
        }
    }

    if (i_one > 0) {
        result_a -= prefix_hash[i_one - 1][j_two].first;
        if (result_a < 0) {
            result_a = p_a + result_a;
        }
        result_b -= prefix_hash[i_one - 1][j_two].second;
        if (result_b < 0) {
            result_b = p_b + result_b;
        }
    }

    if (i_one > 0 and j_one > 0) {
        result_a = static_cast<int>(
            (static_cast<int64_t>(result_a) + prefix_hash[i_one - 1][j_one - 1].first) % p_a);
        result_b = static_cast<int>(
            (static_cast<int64_t>(result_b) + prefix_hash[i_one - 1][j_one - 1].second) % p_b);
    }

    result_a = static_cast<int>(
        (static_cast<int64_t>(result_a) * a_inv_coeff[i_one * M_inp + j_one]) % p_a);
    result_b = static_cast<int>(
        (static_cast<int64_t>(result_b) * b_inv_coeff[i_one * M_inp + j_one]) % p_b);
    return {result_a, result_b};
}

struct Square {
    int i_start;
    int j_start;
    int height;
    Hash hash;
};

bool EqualSquares(const Square& one, const Square& two) {
    for (int i = 0; i < one.height; ++i) {
        for (int j = 0; j < one.height; ++j) {
            if (table[one.i_start + i][one.j_start + j] !=
                table[two.i_start + i][two.j_start + j]) {
                return false;
            }
        }
    }
    return true;
}

bool Exist(int width) {
    std::vector<Square> squares;
    squares.reserve((N_inp - width + 1) * (M_inp - width + 1));
    for (int i = width - 1; i < N_inp; ++i) {
        for (int j = width - 1; j < M_inp; ++j) {
            squares.push_back({i - width + 1, j - width + 1, width,
                               GetHash(i - width + 1, j - width + 1, width)});
        }
    }
    auto cmp = [](const Square& one, const Square& two) {
        if (one.hash.first == two.hash.first) {
            return one.hash.second < two.hash.second;
        }
        return one.hash.first < two.hash.first;
    };
    std::sort(squares.begin(), squares.end(), cmp);

    size_t now = 1;
    while (now < squares.size()) {
        if (squares[now].hash == squares[now - 1].hash) {
            return true;
        }
        ++now;
    }
    return false;
}

std::pair<std::pair<int, int>, std::pair<int, int>> FindExisting(int width) {
    std::vector<Square> squares;
    squares.reserve((N_inp - width + 1) * (M_inp - width + 1));
    for (int i = width - 1; i < N_inp; ++i) {
        for (int j = width - 1; j < M_inp; ++j) {
            squares.push_back({i - width + 1, j - width + 1, width,
                               GetHash(i - width + 1, j - width + 1, width)});
        }
    }
    auto cmp = [](const Square& one, const Square& two) {
        if (one.hash.first == two.hash.first) {
            return one.hash.second < two.hash.second;
        }
        return one.hash.first < two.hash.first;
    };
    std::sort(squares.begin(), squares.end(), cmp);

    size_t prev = 0;
    size_t now = 1;
    while (now < squares.size()) {
        if (squares[now].hash != squares[now - 1].hash) {
            if (now - prev > 1) {
                for (size_t i = prev; i < now; ++i) {
                    for (size_t j = i + 1; j < now; ++j) {
                        if (EqualSquares(squares[i], squares[j])) {
                            return {{squares[i].i_start, squares[i].j_start},
                                      {squares[j].i_start, squares[j].j_start}};
                        }
                    }
                }
            }
            prev = now;
        }
        ++now;
    }

    if (now - prev > 1) {
        for (size_t i = prev; i < now; ++i) {
            for (size_t j = i + 1; j < now; ++j) {
                if (EqualSquares(squares[i], squares[j])) {
                    return {{squares[i].i_start, squares[i].j_start},
                              {squares[j].i_start, squares[j].j_start}};
                }
            }
        }
    }

    return {{-1, -1}, {-1, -1}};
}

std::vector<int> GenerateCoeff(int coeff, int size, int field) {
    std::vector<int> a_coeff(size);
    a_coeff[0] = 1;
    for (int i = 1; i < size; ++i) {
        a_coeff[i] = static_cast<int>((static_cast<int64_t>(a_coeff[i - 1]) * coeff) % field);
    }
    return a_coeff;
}

int main() {
    std::cin >> N_inp >> M_inp;
    table.resize(N_inp);
    for (int i = 0; i < N_inp; ++i) {
        char symbol;
        table[i].resize(M_inp);
        for (int j = 0; j < M_inp; ++j) {
            std::cin >> symbol;
            table[i][j] = static_cast<int>(symbol);
        }
    }

    int a_hash = 2;
    int a_inv = Inv(a_hash, p_a);
    std::vector<int> a_coeff = GenerateCoeff(a_hash, N_inp * M_inp, p_a);
    a_inv_coeff = GenerateCoeff(a_inv, N_inp * M_inp, p_a);

    int b_hash = 3;
    int b_inv = Inv(b_hash, p_b);
    std::vector<int> b_coeff = GenerateCoeff(b_hash, N_inp * M_inp, p_b);
    b_inv_coeff = GenerateCoeff(b_inv, N_inp * M_inp, p_b);

    prefix_hash.resize(N_inp);
    for (int i = 0; i < N_inp; ++i) {
        prefix_hash[i].resize(M_inp);
        for (int j = 0; j < M_inp; ++j) {
            prefix_hash[i][j].first = static_cast<int>(
                (static_cast<int64_t>(a_coeff[i * M_inp + j]) * table[i][j]) % p_a);
            prefix_hash[i][j].second = static_cast<int>(
                (static_cast<int64_t>(b_coeff[i * M_inp + j]) * table[i][j]) % p_b);
        }
    }

    for (int i = 1; i < N_inp; ++i) {
        for (int j = 0; j < M_inp; ++j) {
            prefix_hash[i][j].first = static_cast<int>(
                (static_cast<int64_t>(prefix_hash[i][j].first) + prefix_hash[i - 1][j].first) %
                p_a);
            prefix_hash[i][j].second = static_cast<int>(
                (static_cast<int64_t>(prefix_hash[i][j].second) + prefix_hash[i - 1][j].second) %
                p_b);
        }
    }
    for (int j = 1; j < M_inp; ++j) {
        for (int i = 0; i < N_inp; ++i) {
            prefix_hash[i][j].first = static_cast<int>(
                (static_cast<int64_t>(prefix_hash[i][j].first) + prefix_hash[i][j - 1].first) %
                p_a);
            prefix_hash[i][j].second = static_cast<int>(
                (static_cast<int64_t>(prefix_hash[i][j].second) + prefix_hash[i][j - 1].second) %
                p_b);
        }
    }

    int left = 1;
    int right;
    if (N_inp != M_inp) {
        right = std::min(N_inp, M_inp);
    } else {
        right = N_inp - 1;
    }
    int ans = 0;
    while (left <= right) {
        int middle = (left + right) / 2;
        if (Exist(middle)) {
            ans = middle;
            left = middle + 1;
        } else {
            right = middle - 1;
        }
    }
    std::cout << ans << "\n";
    if (ans != 0) {
        std::cout << answer.first.first + 1 << " " << answer.first.second + 1 << "\n";
        std::cout << answer.second.first + 1 << " " << answer.second.second + 1 << "\n";
    }
}
