/* -*- mia-c++  -*-
 *
 * This file is part of MIA - a toolbox for medical image analysis
 * Copyright (c) Leipzig, Madrid 1999-2017 Gert Wollny
 *
 * MIA is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with MIA; if not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <mia/3d/cost/ncc.hh>

#include <mia/core/threadedmsg.hh>
#include <mia/core/nccsum.hh>
#include <mia/core/parallel.hh>


NS_BEGIN(NS)

using namespace mia;


CNCC3DImageCost::CNCC3DImageCost()
{
       m_copy_to_double = produce_3dimage_filter("convert:repn=double,map=copy");
}

struct FEvaluateNCCSum {
       FEvaluateNCCSum(const C3DDImage& mov, const C3DDImage& ref);
       NCCSums operator ()(const C1DParallelRange& range, const NCCSums& sumacc) const;
private:
       const C3DDImage& m_mov;
       const C3DDImage& m_ref;
};


FEvaluateNCCSum::FEvaluateNCCSum(const C3DDImage& mov, const C3DDImage& ref):
       m_mov(mov), m_ref(ref)
{
}

NCCSums FEvaluateNCCSum::operator ()(const C1DParallelRange& range, const NCCSums& sumacc) const
{
       CThreadMsgStream msks;
       NCCSums sum;

       for (auto z = range.begin(); z != range.end(); ++z) {
              auto ia = m_mov.begin_at(0, 0, z);
              auto ib = m_ref.begin_at(0, 0, z);
              auto eb = m_ref.begin_at(0, 0, z + 1);

              while (ib != eb) {
                     sum.add(*ia, *ib);
                     ++ia;
                     ++ib;
              }
       }

       return sum + sumacc;
};

double CNCC3DImageCost::do_value(const Data& a, const Data& b) const
{
       auto a_double_ptr = m_copy_to_double->filter(a);
       auto b_double_ptr = m_copy_to_double->filter(b);
       const C3DDImage& mov = static_cast<const C3DDImage&>(*a_double_ptr);
       const C3DDImage& ref = static_cast<const C3DDImage&>(*b_double_ptr);
       FEvaluateNCCSum ev(mov, ref);
       NCCSums sum;
       sum = preduce(C1DParallelRange(0, mov.get_size().z, 1), sum, ev,
       [](const NCCSums & x, const NCCSums & y) {
              return x + y;
       });
       return sum.value();
}


class FEvalCostForce : public TFilter<float>
{
       C3DFVectorfield& m_force;
public:
       FEvalCostForce(C3DFVectorfield& force):
              m_force(force)
       {}

       float operator () ( const C3DDImage& mov, const C3DDImage& ref) const
       {
              CThreadMsgStream msks;
              NCCSums sum;
              FEvaluateNCCSum ev(mov, ref);
              sum = preduce(C1DParallelRange(0, mov.get_size().z, 1), sum, ev,
              [](const NCCSums & x, const NCCSums & y) {
                     return x + y;
              });
              auto geval = sum.get_grad_helper();
              auto grad = get_gradient(mov);
              auto grad_eval = [this, &mov, &ref, &grad, &geval](const C1DParallelRange & range) {
                     for (auto z = range.begin(); z != range.end(); ++z) {
                            auto ig = grad.begin_at(0, 0, z);
                            auto iforce = m_force.begin_at(0, 0, z);
                            auto ia = mov.begin_at(0, 0, z);
                            auto ib = ref.begin_at(0, 0, z);
                            auto eb = ref.begin_at(0, 0, z + 1);

                            while (ib != eb) {
                                   *iforce = geval.second.get_gradient_scale(*ia, *ib) * *ig;
                                   ++ig;
                                   ++iforce;
                                   ++ia;
                                   ++ib;
                            }
                     };
              };
              pfor(C1DParallelRange(0, mov.get_size().z, 1), grad_eval);
              return geval.first;
       }

};

double CNCC3DImageCost::do_evaluate_force(const Data& a, const Data& b, Force& force) const
{
       FEvalCostForce ecostforce(force);
       auto a_double_ptr = m_copy_to_double->filter(a);
       auto b_double_ptr = m_copy_to_double->filter(b);
       const C3DDImage& mov = static_cast<const C3DDImage&>(*a_double_ptr);
       const C3DDImage& ref = static_cast<const C3DDImage&>(*b_double_ptr);
       return ecostforce(mov, ref);
}


CNCC3DImageCostPlugin::CNCC3DImageCostPlugin():
       C3DImageCostPlugin("ncc")
{
}

C3DImageCost *CNCC3DImageCostPlugin::do_create() const
{
       return new CNCC3DImageCost();
}

const std::string CNCC3DImageCostPlugin::do_get_descr() const
{
       return "normalized cross correlation.";
}

extern "C" EXPORT CPluginBase *get_plugin_interface()
{
       return new CNCC3DImageCostPlugin();
}

NS_END
