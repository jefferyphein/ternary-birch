#ifndef __SET_COVER_H_
#define __SET_COVER_H_

#include <vector>

class SetCover
{
public:
    SetCover() = default;
    SetCover(Z64 total, const std::vector<std::vector<W64>>& vectors, Z64 method=SetCover::METHOD_GREEDY);

    const std::vector<Z64>& positions(void) const;
    size_t num_positions(void) const;

    static constexpr Z64 METHOD_GREEDY = 0;
    static constexpr Z64 METHOD_KINDA_GREEDY = 1;
    static constexpr Z64 METHOD_BRUTE_FORCE = 2;

private:
    std::vector<Z64> positions_;
    Z64 total_;

    void greedy(const std::vector<std::vector<W64>>& vectors);
    void kinda_greedy(const std::vector<std::vector<W64>>& vectors);
    void brute_force(const std::vector<std::vector<W64>>& vectors);

    bool is_set_cover(const std::vector<std::vector<W64>>& vectors, const std::vector<Z64>& pivots) const;
};

#endif // __SET_COVER_H_
