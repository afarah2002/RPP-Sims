#include <boost/python.hpp>
#include "G4QuadrupoleMagField.hh"

using namespace boost::python


// ====================================================================
// module definition
// ====================================================================

void export_G4QuadrupoleMagField()
{
	class_<G4QuadrupoleMagField, G4QuadrupoleMagField*,
		bases<G4MagneticField, G4ElectroMagneticField, G4Field> >
		("G4QuadrupoleMagField", "quadrupole magnetic field",  no_init)
		// constructions
		// .def(init<const G4ThreeVector&>())
		// .def(init<const G4double, G4double, G4double>())
		// ---
		.def("GetFieldValue", &G4QuadrupoleMagField::GetFieldValue)
		;
}