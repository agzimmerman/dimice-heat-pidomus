#include <iostream>
#include <fstream>

class LogOStream
{
public:
  LogOStream() : log_fstream("log.txt") {}; // check if opening file succeeded!!
  // for regular output of variables and stuff
  template<typename T> LogOStream& operator<<(const T& something)
  {
    std::cout << something;
    log_fstream << something;
    return *this;
  }
  // for manipulators like std::endl
  typedef std::ostream& (*stream_function)(std::ostream&);
  LogOStream& operator<<(stream_function func)
  {
    func(std::cout);
    func(log_fstream);
    return *this;
  }
private:
  std::ofstream log_fstream;
};