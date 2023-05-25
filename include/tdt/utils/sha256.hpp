#pragma once

#include <array>

#include <kassert/kassert.hpp>
#include <openssl/evp.h>

#include "tdt/assertion_levels.hpp"
#include "tdt/checking_casts.hpp"
#include "tdt/utils/concepts.hpp"

// using Sha256Digest = std::basic_string<unsigned char>;
using Sha256Digest = std::array<unsigned char, 32>;
class Sha256 {
public:
    Sha256() {
        _ssl_context = EVP_MD_CTX_new();
        KASSERT(_ssl_context != nullptr);
        reset();
    }

    ~Sha256() {
        EVP_MD_CTX_free(_ssl_context);
    }

    template <typename T>
    requires TriviallyCopyable<T>
    void update(T&& data) {
        [[maybe_unused]] int ret = EVP_DigestUpdate(_ssl_context, &data, sizeof(data));
        KASSERT(ret == 1, "EVP_DigestUpdate failed", tdt::assert::light);
    }

    template <typename T>
    requires requires(T t) {
        t.data();
        t.size();
    }
    void update(T& container) {
        KASSERT(_ssl_context, "Sha256 context is null", tdt::assert::light);
        const size_t          count = container.size() * sizeof(typename T::value_type);
        [[maybe_unused]] auto ret   = EVP_DigestUpdate(_ssl_context, container.data(), count);
        KASSERT(ret == 1, "EVP_DigestUpdate failed", tdt::assert::normal);
    }

    void digest(Sha256Digest& digest) {
        // TODO What's the difference between digest_len and digest_size?
        KASSERT(_ssl_context, "Sha256 context is null", tdt::assert::light);

        size_t const digest_size = asserting_cast<size_t>(EVP_MD_size(EVP_sha256()));
        unsigned int digest_len  = 0;

        [[maybe_unused]] auto ret = EVP_DigestFinal_ex(_ssl_context, digest.data(), &digest_len);
        KASSERT(ret == 1, "EVP_DigestFinal_ex failed", tdt::assert::light);
        KASSERT(digest_len == 32u, "Digest length is not 32", tdt::assert::light);
        KASSERT(digest_len == digest_size, "Digest length is not equal to digest size", tdt::assert::light);
        KASSERT(digest_size > 0ul, "Digest length is 0", tdt::assert::light);
    }

    void reset() {
        KASSERT(_ssl_context, "Sha256 context is null", tdt::assert::light);
        [[maybe_unused]] auto ret = EVP_DigestInit(_ssl_context, EVP_sha256());
        KASSERT(ret == 1, "EVP_DigestInit failed", tdt::assert::light);
    }

    Sha256Digest digest() {
        Sha256Digest digest;
        return digest;
    }

private:
    EVP_MD_CTX* _ssl_context = nullptr;
};
