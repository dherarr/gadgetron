#ifndef IMAGEWRITERGADGET_H
#define IMAGEWRITERGADGET_H

#include <complex>

#include "Gadget.h"
#include "NDArray.h"
#include "gadgetheaders.h"


class ImageWriterGadget :
public Gadget2<GadgetMessageImage, NDArray< std::complex<float> > >
{
 public:
  ImageWriterGadget()
    : calls_(0)
    {}

 protected:
  virtual int process( GadgetContainerMessage< GadgetMessageImage>* m1,
		       GadgetContainerMessage< NDArray< std::complex<float> > >* m2);

  long calls_;

};

#endif //IMAGEWRITERGADGET_H