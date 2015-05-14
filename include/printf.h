#pragma once

template<typename OS>
static void fprintf(OS &os, const char *s)
{
  while (*s) {
    if (*s == '%') {
      if (*(s + 1) == '%') {
        ++s;
      }
      else {
        throw "invalid format string: missing arguments";
      }
    }
    os << *s++;
  }
}

template<typename OS, typename T, typename... Args>
static void fprintf(OS &os, const char *s, T& value, Args... args)
{
  while (*s) {
    if (*s == '%') {
      if (*(s + 1) == '%') {
        ++s;
      }
      else {
        os << value;
        fprintf(os, s + 1, args...); // call even when *s == 0 to detect extra arguments
        return;
      }
    }
    os << *s++;
  }
  throw "extra arguments provided to printf";
}

