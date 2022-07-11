#ifndef DIP_ARCH_H
#define DIP_ARCH_H

#define DIP_ARCH Linuxa64

/* The unsigned integer types (uint) */
typedef unsigned char dip_uint8;
typedef unsigned short dip_uint16;
typedef unsigned int dip_uint32;

/* The signed integer types (sint) */
typedef signed char dip_sint8;
typedef signed short dip_sint16;
typedef signed int dip_sint32;

/* This is a 64-bit machine */
#define DIP_64BITS
typedef unsigned long dip_uint64;
typedef signed long dip_sint64;

/* The generic integer types */
typedef dip_sint64 dip_int;
typedef dip_sint64 dip_sint;
typedef dip_uint64 dip_uint;


#define _POSIX_SOURCE 1

#define _GNU_SOURCE 1

#define DIP_PORT_HAS_POSIX_TIME
#define DIP_PORT_HAS_CHMOD
#define DIP_PORT_HAS_STRCASECMP

#define DIP_EXPORT
#define DIPIO_EXPORT
#define DML_EXPORT

#define DIP_ERROR DIP_EXPORT dip_Error
#define DIPIO_ERROR DIPIO_EXPORT dip_Error

#define DIPIO_DIRECTORY_SEPARATOR '/'
#define DIPIO_EXTENSION_SEPARATOR '.'

#endif
