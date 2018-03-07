/*
 * Common.h
 *	A common file to declare variables which are globally used
 *  Created on: Feb 26, 2018
 *      Author: shiyongtan
 */

#ifndef INC_LIBSTB_COMMON_H_
#define INC_LIBSTB_COMMON_H_

// Data for debug
enum DebugMode {
		NO_SKIP = 0,
		SKIP_IPR_2D_POSITION,
		SKIP_IPR_TRIANGULATION,
		SKIP_IPR_SHAKING,
		SKIP_IPR,
		SKIP_PREDITIVE_FIELD,
		SKIP_INITIAL_PHASE,
		SKIP_PREVIOUS_TRACKS
};
extern DebugMode debug_mode;
extern int debug_frame_number;

//Data for error
enum ERROR {
	NONE = 0,
	NO_FILE, //Can't find required file
};

extern ERROR error;


#endif /* INC_LIBSTB_COMMON_H_ */
