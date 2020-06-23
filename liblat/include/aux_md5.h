// 2014-2019, ETH Zurich, Integrated Systems Laboratory
// Authors: Christian Stieger

/** @addtogroup liblat
 * @{
 */

#ifndef _AUX_MD5_
#define _AUX_MD5_

#include "aux_io.h"
#include <openssl/md5.h>
#include <string>
#include <sstream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cstring>

namespace aux {

	/** compute md5 string from raw bytes
	 * @param buff array of raw bytes
	 * @param len length of buff
	 */
	inline std::string md5(const char* buff, const size_t len) noexcept {
		
		// MD5 magic
		unsigned char* digest[MD5_DIGEST_LENGTH];
		MD5((unsigned char*)buff,len,(unsigned char*)digest);

		// produce result
		char res[2*MD5_DIGEST_LENGTH+1];
		for (auto id = (unsigned char*)digest,
			  ie = (unsigned char*)digest+MD5_DIGEST_LENGTH,
			  j  = (unsigned char*)res; id!=ie; ++id, j+=2)
			sprintf((char*)j,"%02x",(unsigned int)(*id));

		return std::string(res);
	}

	//! get md5 for a file
	inline std::string md5(const std::string& fileName) {
		
		// open file
		auto file = aux::openFile<std::ifstream>(fileName);

		// read file contens into string
		std::string raw((std::istreambuf_iterator<char>(file)),
				 std::istreambuf_iterator<char>());

		return aux::md5(raw.c_str(),raw.size());
	}
}

#endif // _AUX_MD5_

/** @}
 */
