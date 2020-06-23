// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _AUX_IO_
#define _AUX_IO_

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>


// ansi color codes
#ifndef NCOLOR
	#define RESET__ "\033[0m"		//!< reset color
	
	#define BLACK__ "\033[1;30m"		//!< set black
	#define RED__ "\033[1;31m"		//!< set red
	#define GREEN__ "\033[1;32m"		//!< set green
	#define YELLOW__ "\033[1;33m"		//!< set yellow
	#define BLUE__ "\033[1;34m"		//!< set blue
	#define MAGENTA__ "\033[1;35m"		//!< set magenta
	#define CYAN__ "\033[1;36m"		//!< set cyan
	#define WHITE__ "\033[1;37m"		//!< set white
#else
	#define RESET__ ""
	
	#define BLACK__ ""
	#define RED__ ""
	#define GREEN__ ""
	#define YELLOW__ ""
	#define BLUE__ ""
	#define MAGENTA__ ""
	#define CYAN__ ""
	#define WHITE__ ""
#endif

namespace aux {
	template <typename T>
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& inp) noexcept;
}

namespace aux {
	
	//! returns timestamp string
	inline std::string timeStamp() noexcept {
		
		std::string res;
		{
			std::locale::global(std::locale("en_US.utf8"));
			std::time_t t = std::time(nullptr);
			char mbstr[100];
			std::strftime(mbstr, sizeof(mbstr), "%c", std::localtime(&t));
			res = mbstr;
		}
		
		std::replace(res.begin(),res.end(),' ','-');
		return res;
	}


	/** Open file helper function. \n
	 * Template arguments are: \n
	 * - ST: filestream, must be std::ifstream or std::ofstream
	 * @param fileName name of the file to be opened
	 * @param mode open mode
	 */
	template <class ST>
	inline ST openFile(const std::string& fileName,
			     const std::ios_base::openmode mode = std::ios_base::out) {
		// open file and complain
		ST file;
		file.open(fileName,mode);
		if (!file.good())
			throw(std::invalid_argument("file \'"+fileName+"\' not found"));

		// imbue to print commas every 6 digits
		class cnp_: public std::numpunct<char> {
		protected:
			virtual char do_thousands_sep() const { return ','; }
			virtual std::string do_grouping() const { return "\06"; }
		};
		file.imbue(std::locale(std::locale(), new cnp_));

		return file;
	}
	
	
	namespace detail {
		/** prints a range with whitespace
		 * @param os the output stream
		 * @param b iterator to the beginning of the range
		 * @param e iterator past the end of the range
		 * @param ws the whitespace character
		 */
		template<class ITR>
		std::ostream& pws(std::ostream& os, ITR b, ITR e, const char ws) {
			using std::distance;
			if (!distance(b,e)) return os;
			--e;
			while (b!=e) os << *b++ << ws;
			return (os << *e);
		};
	}


	//! vector printing
	template <typename T>
	std::ostream& operator<<(std::ostream& os, const std::vector<T>& inp) noexcept {
		return detail::pws(os,inp.cbegin(),inp.end(),'\n');
	}
	//! vector printing for size_t
	template<>
	inline std::ostream& operator<<<size_t>(std::ostream& os, const std::vector<size_t>& inp) noexcept {
		return detail::pws(os,inp.cbegin(),inp.cend(),' ');
	}
	//! vector printing for double
	template<>
	inline std::ostream& operator<<<double>(std::ostream& os, const std::vector<double>& inp) noexcept {
		return detail::pws(os,inp.cbegin(),inp.cend(),' ');
	}
	//! vector printing for bool
	template<>
	inline std::ostream& operator<<<bool>(std::ostream& os, const std::vector<bool>& inp) noexcept {
		return detail::pws(os,inp.cbegin(),inp.cend(),' ');
	}
	//! vector printing for std::string
	template<>
	inline std::ostream& operator<<<std::string>(std::ostream& os,
						const std::vector<std::string>& inp) noexcept {
		return detail::pws(os,inp.cbegin(),inp.cend(),' ');
	}
}



//! Tee buffer class.
template <class CT, class TT = std::char_traits<CT>>
class aux_b_teebuf: public std::basic_streambuf<CT,TT> {
public:
	/** @name types
	 */
	typedef typename TT::int_type int_type;		//!< integer type
	typedef typename TT::char_type char_type;	//!< char type

public:
	/** @name constructors
	 */
	//! constructor from streambuffers
	aux_b_teebuf(std::basic_streambuf<CT,TT>* sb1,
		     std::basic_streambuf<CT,TT>* sb2):
	sb1_(sb1), sb2_(sb2) {}

private:
	//! stream sync
	virtual int sync() noexcept {
		const int r1 = sb1_->pubsync();
		const int r2 = sb2_->pubsync();
		return !r1 && !r2 ? 0: -1;
	}
	//! overflow to two streams
	virtual int_type overflow(const int_type c) noexcept {
		
		const int_type eof = TT::eof();
		if (TT::eq_int_type(c,eof))
			return TT::not_eof(c);

		const char_type ch = TT::to_char_type(c);
		const int_type r1 = sb1_->sputc(ch);
		const int_type r2 = sb2_->sputc(ch);

		return TT::eq_int_type(r1,eof) ||
		       TT::eq_int_type(r2,eof) ? eof: c;
	}

private:
	/** @name member variables
	 */
	std::basic_streambuf<CT,TT>* sb1_;	//!< first streambuffer
	std::basic_streambuf<CT,TT>* sb2_;	//!< second streambuffer
};

//! teebuffer shorthand
typedef aux_b_teebuf<char> aux_teebuf;


//! Tee stream class.
class aux_tee: public std::ostream {
public:
	/** @name constructors
	 */
	//! constructor from two streams
	inline aux_tee(std::ostream& os1, std::ostream& os2):
		std::ostream(&tbuf_), tbuf_(os1.rdbuf(),os2.rdbuf()) {}

private:
	/** @name member variables
	 */
	aux_teebuf tbuf_;	//!< tee streambuffer
};

#endif // _AUX_IO_

/** @}
 */
