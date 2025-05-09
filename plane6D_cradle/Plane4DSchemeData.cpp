#include "Plane4DSchemeData.hpp"
#include <algorithm>
#include <typeinfo>
#include <macro.hpp>
bool Plane4DSchemeData::operator==(const Plane4DSchemeData& rhs) const {
	if (this == &rhs) return true;
	else if (!SchemeData::operator==(rhs)) return false;
	else if (wMax != rhs.wMax) return false;
	else if ((aranges.size() != rhs.aranges.size()) || !std::equal(aranges.cbegin(), aranges.cend(), rhs.aranges.cbegin())) return false;
	else return true;
}
bool Plane4DSchemeData::operator==(const SchemeData& rhs) const {
	if (this == &rhs) return true;
	else if (typeid(*this) != typeid(rhs)) return false;
	else return operator==(dynamic_cast<const Plane4DSchemeData&>(rhs));
}

bool Plane4DSchemeData::hamFunc(
	beacls::UVec& hamValue_uvec,
	const FLOAT_TYPE,
	const beacls::UVec&,
	const std::vector<beacls::UVec>& derivs,
	const size_t begin_index,
	const size_t length
	) const {
	const levelset::HJI_Grid *hji_grid = get_grid();

	const beacls::FloatVec &xs0 = hji_grid->get_xs(0);
	const beacls::FloatVec &xs1 = hji_grid->get_xs(1);
	const beacls::FloatVec &xs2 = hji_grid->get_xs(2);
	const beacls::FloatVec &xs3 = hji_grid->get_xs(3);
	const beacls::FloatVec &xs4 = hji_grid->get_xs(4);
	const beacls::FloatVec &xs5 = hji_grid->get_xs(5);
	beacls::reallocateAsSrc(hamValue_uvec, derivs[0]);
	FLOAT_TYPE* hamValue = beacls::UVec_<FLOAT_TYPE>(hamValue_uvec).ptr();
	const FLOAT_TYPE* deriv0 = beacls::UVec_<FLOAT_TYPE>(derivs[0]).ptr();
	const FLOAT_TYPE* deriv1 = beacls::UVec_<FLOAT_TYPE>(derivs[1]).ptr();
	const FLOAT_TYPE* deriv2 = beacls::UVec_<FLOAT_TYPE>(derivs[2]).ptr();
	const FLOAT_TYPE* deriv3 = beacls::UVec_<FLOAT_TYPE>(derivs[3]).ptr();
	const FLOAT_TYPE* deriv4 = beacls::UVec_<FLOAT_TYPE>(derivs[4]).ptr();
	const FLOAT_TYPE* deriv5 = beacls::UVec_<FLOAT_TYPE>(derivs[5]).ptr();

	for (size_t i = 0; i<length; ++i) {
		FLOAT_TYPE deriv0_i = deriv0[i];
		FLOAT_TYPE deriv1_i = deriv1[i];
		FLOAT_TYPE deriv2_i = deriv2[i];
		FLOAT_TYPE deriv3_i = deriv3[i];
		FLOAT_TYPE deriv4_i = deriv4[i];
		FLOAT_TYPE deriv5_i = deriv5[i];
		FLOAT_TYPE xs0_i = xs0[begin_index + i];
		FLOAT_TYPE xs1_i = xs1[begin_index + i];
		FLOAT_TYPE xs2_i = xs2[begin_index + i];
		FLOAT_TYPE xs3_i = xs3[begin_index + i];
		FLOAT_TYPE xs4_i = xs4[begin_index + i];
		FLOAT_TYPE xs5_i = xs5[begin_index + i];
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.5;
		FLOAT_TYPE Vd = 1.5;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.05; 
		FLOAT_TYPE Yv = 15/37; 
			
		if (deriv0_i == 0 && deriv1_i == 0){
			hamValue[i] = -((-V0 + xs2_i ) * deriv0_i
			- Kp * xs0_i * deriv2_i
			- Kd * (-V0 + xs2_i) * deriv2_i
			+ Vd * HjiFabs(deriv0_i)
			+ Wd * HjiFabs(deriv2_i)
			+ xs3_i * deriv1_i
			- Kp * xs1_i * deriv3_i
			- Kd * (xs3_i) * deriv3_i
			+ Vd * HjiFabs(deriv1_i)
			+ Wd * HjiFabs(deriv3_i)
			+ xs5_i * deriv4_i
			+ (-Yv*(xs3_i*cos(xs4_i) - xs2_i*sin(xs4_i) - xs5_i) - 2*(xs2_i*cos(xs4_i) + xs3_i*sin(xs4_i))*xs5_i + (-(- Kp * xs0_i - Kd * (- V0 + xs2_i))*sin(xs4_i) + (- Kp * xs1_i - Kd * (xs3_i))*cos(xs4_i))) * deriv5_i
			+ Wr * HjiFabs(deriv5_i));


		}
		else if (deriv0_i == 0 && deriv1_i != 0){
			hamValue[i] = -((-V0 + xs2_i ) * deriv0_i
			- Kp * xs0_i * deriv2_i
			- Kd * (-V0 + xs2_i) * deriv2_i
			+ Vd * HjiFabs(deriv0_i)
			+ Wd * HjiFabs(deriv2_i)
			+ xs3_i * deriv1_i
			- Kp * xs1_i * deriv3_i
			- Kd * (xs3_i + Vd * deriv1_i/HjiFabs(deriv1_i)) * deriv3_i
			+ Vd * HjiFabs(deriv1_i)
			+ Wd * HjiFabs(deriv3_i)
			+ xs5_i * deriv4_i
			+ (-Yv*(xs3_i*cos(xs4_i) - xs2_i*sin(xs4_i) - xs5_i) - 2*(xs2_i*cos(xs4_i) + xs3_i*sin(xs4_i))*xs5_i + (-(- Kp * xs0_i - Kd * (- V0 + xs2_i))*sin(xs4_i) + (- Kp * xs1_i - Kd * (xs3_i + Vd*deriv1_i/HjiFabs(deriv1_i)))*cos(xs4_i))) * deriv5_i
			+ Wr * HjiFabs(deriv5_i));



		}
		else if (deriv0_i != 0 && deriv1_i == 0){
			hamValue[i] = -((-V0 + xs2_i ) * deriv0_i
			- Kp * xs0_i * deriv2_i
			- Kd * (-V0 + xs2_i + Vd * deriv0_i/HjiFabs(deriv0_i)) * deriv2_i
			+ Vd * HjiFabs(deriv0_i)
			+ Wd * HjiFabs(deriv2_i)
			+ xs3_i * deriv1_i
			- Kp * xs1_i * deriv3_i
			- Kd * (xs3_i) * deriv3_i
			+ Vd * HjiFabs(deriv1_i)
			+ Wd * HjiFabs(deriv3_i)
			+ xs5_i * deriv4_i
			+ (-Yv*(xs3_i*cos(xs4_i) - xs2_i*sin(xs4_i) - xs5_i) - 2*(xs2_i*cos(xs4_i) + xs3_i*sin(xs4_i))*xs5_i + (-(- Kp * xs0_i - Kd * (- V0 + xs2_i + Vd*deriv0_i/HjiFabs(deriv0_i)))*sin(xs4_i) + (- Kp * xs1_i - Kd * (xs3_i))*cos(xs4_i))) * deriv5_i
			+ Wr * HjiFabs(deriv5_i));

		}
		else{
			hamValue[i] = -((-V0 + xs2_i ) * deriv0_i
			- Kp * xs0_i * deriv2_i
			- Kd * (-V0 + xs2_i + Vd * deriv0_i/HjiFabs(deriv0_i)) * deriv2_i
			+ Vd * HjiFabs(deriv0_i)
			+ Wd * HjiFabs(deriv2_i)
			+ xs3_i * deriv1_i
			- Kp * xs1_i * deriv3_i
			- Kd * (xs3_i + Vd * deriv1_i/HjiFabs(deriv1_i)) * deriv3_i
			+ Vd * HjiFabs(deriv1_i)
			+ Wd * HjiFabs(deriv3_i)
			+ xs5_i * deriv4_i
			+ (-Yv*(xs3_i*cos(xs4_i) - xs2_i*sin(xs4_i) - xs5_i) - 2*(xs2_i*cos(xs4_i) + xs3_i*sin(xs4_i))*xs5_i + (-(- Kp * xs0_i - Kd * (- V0 + xs2_i + Vd*deriv0_i/HjiFabs(deriv0_i)))*sin(xs4_i) + (- Kp * xs1_i - Kd * (xs3_i + Vd*deriv1_i/HjiFabs(deriv1_i)))*cos(xs4_i))) * deriv5_i
			+ Wr * HjiFabs(deriv5_i));		

		
		}
	}
	return true;
}
bool Plane4DSchemeData::partialFunc(
	beacls::UVec& alphas_uvec,
	const FLOAT_TYPE,
	const beacls::UVec&,
	const std::vector<beacls::UVec>&,
	const std::vector<beacls::UVec>&,
	const size_t dim,
	const size_t begin_index,
	const size_t length
) const {
	if (alphas_uvec.type() != beacls::UVecType_Vector) alphas_uvec = beacls::UVec(beacls::type_to_depth<FLOAT_TYPE>(), beacls::UVecType_Vector, length);
	else alphas_uvec.resize(length);
	FLOAT_TYPE* alphas = beacls::UVec_<FLOAT_TYPE>(alphas_uvec).ptr();
	const levelset::HJI_Grid *hji_grid = get_grid();
	switch (dim) {
	case 0:
	{
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.3;
		FLOAT_TYPE Vd = 1.0;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.1; 
		FLOAT_TYPE Yv = 15/37; 
		const beacls::FloatVec &xs2 = hji_grid->get_xs(2);
		for (size_t i = 0; i<length; ++i) {
			// alphas[i] = HjiFabs(xs3[begin_index + i] * cos(xs2[begin_index + i]));
			alphas[i] = HjiFabs(V0 + Vd + HjiFabs(xs2[begin_index + i]));
		}
	}
	break;
	case 1:
	{
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.3;
		FLOAT_TYPE Vd = 1.0;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.1; 
		FLOAT_TYPE Yv = 15/37; 
		const beacls::FloatVec &xs3 = hji_grid->get_xs(3);
		for (size_t i = 0; i<length; ++i) {
			alphas[i] = HjiFabs(V0 + Vd + HjiFabs(xs3[begin_index + i]));
		}
	}
	break;
	case 2:
	{
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.3;
		FLOAT_TYPE Vd = 1.0;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.1; 
		FLOAT_TYPE Yv = 15/37; 
		const beacls::FloatVec &xs0 = hji_grid->get_xs(0);
		for (size_t i = 0; i<length; ++i) {
			alphas[i] = Kp*HjiFabs(xs0[begin_index + i]) + Wd;
		}
	}
	break;
	case 3:
	{
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.3;
		FLOAT_TYPE Vd = 1.0;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.1; 
		FLOAT_TYPE Yv = 15/37; 
		const beacls::FloatVec &xs1 = hji_grid->get_xs(1);
		for (size_t i = 0; i<length; ++i) {
			alphas[i] = Kp*HjiFabs(xs1[begin_index + i]) + Wd;
		}
	}
	break;
	case 4:
	{
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.3;
		FLOAT_TYPE Vd = 1.0;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.1; 
		FLOAT_TYPE Yv = 15/37; 
		const beacls::FloatVec &xs5 = hji_grid->get_xs(5);
		for (size_t i = 0; i<length; ++i) {
			alphas[i] = HjiFabs(xs5[begin_index + i]);
		}
	}
	break;
	case 5:
	{
		FLOAT_TYPE V0 = 2.5;
		FLOAT_TYPE Wd = 0.3;
		FLOAT_TYPE Vd = 1.0;
		FLOAT_TYPE Kp = 1.0;
		FLOAT_TYPE Kd = 1.0;
		FLOAT_TYPE Wr = 0.1; 
		FLOAT_TYPE Yv = 15/37; 
		const beacls::FloatVec &xs0 = hji_grid->get_xs(0);
		const beacls::FloatVec &xs1 = hji_grid->get_xs(1);
		const beacls::FloatVec &xs2 = hji_grid->get_xs(2);
		const beacls::FloatVec &xs3 = hji_grid->get_xs(3);
		const beacls::FloatVec &xs4 = hji_grid->get_xs(4);
		const beacls::FloatVec &xs5 = hji_grid->get_xs(5);
		for (size_t i = 0; i<length; ++i) {
			FLOAT_TYPE xs0_i = xs0[begin_index + i];
			FLOAT_TYPE xs1_i = xs1[begin_index + i];
			FLOAT_TYPE xs2_i = xs2[begin_index + i];
			FLOAT_TYPE xs3_i = xs3[begin_index + i];
			FLOAT_TYPE xs4_i = xs4[begin_index + i];
			FLOAT_TYPE xs5_i = xs5[begin_index + i];	
			alphas[i] = (HjiFabs(Yv*(xs3_i*cos(xs4_i) - xs2_i*sin(xs4_i) - xs5_i)) + HjiFabs(((- Kp * xs0_i - Kd * (- V0 + xs2_i + Vd))*cos(xs4_i) - (- Kp * xs1_i - Kd * (xs3_i + Vd))*sin(xs4_i))));

		}
	}
	break;		
	default:
		printf("error:%s,%dPartials for the game of two identical vehicles only exist in dimensions 1-3\n", __func__, __LINE__);
		return false;
	}
	return true;
}
