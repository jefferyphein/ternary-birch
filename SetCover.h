#ifndef __SET_COVER_H_
#define __SET_COVER_H_

#include <vector>

class SetCover
{
public:
    SetCover() = default;
    SetCover(int64_t total, const std::vector<std::vector<uint64_t>>& vectors, int64_t method=SetCover::METHOD_GREEDY);

    const std::vector<int64_t>& positions(void) const;
    size_t num_positions(void) const;

    static constexpr int64_t METHOD_GREEDY = 0;
    static constexpr int64_t METHOD_KINDA_GREEDY = 1;
    static constexpr int64_t METHOD_BRUTE_FORCE = 2;

private:
    std::vector<int64_t> positions_;
    int64_t total_;

    void greedy(const std::vector<std::vector<uint64_t>>& vectors);
    void kinda_greedy(const std::vector<std::vector<uint64_t>>& vectors);
    void brute_force(const std::vector<std::vector<uint64_t>>& vectors);

    bool is_set_cover(const std::vector<std::vector<uint64_t>>& vectors, const std::vector<int64_t>& pivots) const;
};

#endif // __SET_COVER_H_
