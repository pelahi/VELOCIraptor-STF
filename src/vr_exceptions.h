#include <exception>
#include <string>

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

// Forward declaration
namespace NBody
{
class Particle;
}

namespace vr
{

/// Base class for all exceptions in velociraptor
class exception : public std::exception
{
public:
    exception(std::string what)
      : m_what(std::move(what))
    {}

    const char *what() const noexcept override
    {
        return m_what.c_str();
    }

private:
    std::string m_what;
};

/// Exception thrown when a particle with non-positive density is found
class non_positive_density : public exception
{
public:
    non_positive_density(const NBody::Particle &part, const std::string &function);
};

} // namespace vr

#endif
