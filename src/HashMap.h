#ifndef __HASHMAP_H_
#define __HASHMAP_H_

#include "birch.h"

template<typename Key>
class HashMap
{
public:
    HashMap(Z64 capacity = 1LL << DEFAULT_LG2_CAPACITY)
    {
        int lg2_capacity = (int)ceil(log2(capacity));
        if (lg2_capacity < DEFAULT_LG2_CAPACITY)
        {
            lg2_capacity = DEFAULT_LG2_CAPACITY;
        }

        this->capacity_ = 1LL << lg2_capacity;
        this->mask = (1LL << (lg2_capacity+1))-1;
        this->keys_.reserve(this->capacity_);
        this->vals.reserve(this->capacity_);
        this->keyptr.resize(this->capacity_ << 1, -1);
        this->num_stored = 0;
    }

    size_t capacity(void) const { return this->capacity_; }
    size_t size(void) const { return this->keys_.size(); }

    bool add(const Key& key)
    {
        W64 value = std::hash<Key>{}(key);
        return this->add(key, value, true);
    }

    bool add(const Key& key, W64 value)
    {
        return this->add(key, value, true);
    }

    const Key& get(size_t index) const
    {
        return this->keys_[index];
    }

    bool exists(const Key& key) const
    {
        W64 value = std::hash<Key>{}(key);
        size_t index = value & this->mask;
        while (1)
        {
            Z64 offset = this->keyptr[index];
            if (offset == -1)
            {
                return false;
            }

            if (key == this->keys_[offset])
            {
                return true;
            }

            index = (index + 1) & this->mask;
        }

        return false;
    }

    size_t indexof(const Key& key) const
    {
        W64 value = std::hash<Key>{}(key);
        size_t index = value & this->mask;
        while (1)
        {
            Z64 offset = this->keyptr[index];
            if (offset == -1)
            {
                throw std::invalid_argument("Key not found.");
            }

            if (key == this->keys_[offset])
            {
                return offset;
            }

            index = (index + 1) & this->mask;
        }

        return -1;
    }

    Key& at(size_t index)
    {
        return this->keys_[index];
    }

    const Key& last(void) const
    {
        return this->keys_.back();
    }

    const std::vector<Key>& keys(void) const
    {
        return this->keys_;
    }

private:
    bool add(const Key& key, W64 value, bool do_push_back)
    {
        Z64 index = value & this->mask;
        while (1)
        {
            Z64 offset = this->keyptr[index];
            if (offset == -1)
            {
                return this->insert(key, value, index, do_push_back);
            }

            if (key == this->keys_[offset])
            {
                return false;
            }

            index = (index + 1) & this->mask;
        }
        return false;
    }

    bool insert(const Key& key, W64 value, Z64 index, bool do_push_back)
    {
        if (this->num_stored == this->capacity_)
        {
            this->expand();
            return this->add(key);
        }

        Z64 offset = this->num_stored;
        ++this->num_stored;

        if (do_push_back)
        {
            this->keys_.emplace_back(key);
            this->vals.emplace_back(value);
        }
        this->keyptr[index] = offset;

        return true;
    }

    void expand(void)
    {
        this->capacity_ <<= 1;
        this->mask = (this->capacity_ << 1)-1;
        this->keys_.reserve(this->capacity_);
        this->vals.reserve(this->capacity_);

        this->keyptr.clear();
        this->keyptr.resize(this->capacity_ << 1, -1);

        Z64 stored = this->num_stored;
        this->num_stored = 0;

        for (Z64 n=0; n<stored; n++)
        {
            this->add(this->keys_[n], this->vals[n], false);
        }
    }

    size_t capacity_;
    W64 mask;
    size_t num_stored;
    std::vector<Key> keys_;
    std::vector<W64> vals;
    std::vector<Z64> keyptr;

    static constexpr int DEFAULT_LG2_CAPACITY = 4;
};

#endif // __HASHMAP_H_
