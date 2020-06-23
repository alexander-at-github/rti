#pragma once

#include "rti/trace/hit_accumulator.hpp"

namespace rti { namespace trace {
  template<typename numeric_type>
  class hit_accumulator_combination : public rti::trace::i_hit_accumulator<numeric_type> {
    //using internal_numeric_type = double;
    using internal_numeric_type = typename rti::trace::i_hit_accumulator<numeric_type>::internal_numeric_type;
  public:
    hit_accumulator_combination
    (i_hit_accumulator<numeric_type>& h1, i_hit_accumulator<numeric_type>& h2) :
      h1(h1),
      h2(h2) {
      if (h1.get_size() != h2.get_size()) {
        throw std::invalid_argument("arguments are not of equal size");
      }
      prepare_values();
    }

  private:
    void prepare_values() {
      size = h1.get_size();
      cnts.reserve(size);
      energyvalues.reserve(size);
      relativeerrors.reserve(size);
      vovs.reserve(size);
      exposedareas.reserve(size);
      cntsSum = 0;

      auto h1cnts = h1.get_cnts();
      auto h2cnts = h2.get_cnts();
      auto h1evs = h1.get_energy_values();
      auto h2evs = h2.get_energy_values();
      auto h1res = h1.get_relative_error();
      auto h2res = h2.get_relative_error();
      auto h1vovs = h1.get_vov();
      auto h2vovs = h2.get_vov();
      auto h1expa = h1.get_exposed_areas();
      auto h2expa = h2.get_exposed_areas();
      for (size_t idx; idx < size; ++idx) {
        auto vv = (internal_numeric_type) 0;
        if (h1res[idx] < h2res[idx]) {
          cnts.push_back(h1cnts[idx]);
          energyvalues.push_back(h1evs[idx]);
          relativeerrors.push_back(h1res[idx]);
          vovs.push_back(h1vovs[idx]);
          exposedareas.push_back(h1expa[idx]);
          //RLOG_TRACE << "F";
        } else {
          cnts.push_back(h2cnts[idx]);
          energyvalues.push_back(h2evs[idx]);
          relativeerrors.push_back(h2res[idx]);
          vovs.push_back(h2vovs[idx]);
          exposedareas.push_back(h2expa[idx]);
          //RLOG_TRACE << "S";
        }
        assert(cnts.size() == idx+1 && "Correctness Assertion");
        cntsSum += cnts[idx];
      }
      cnts.shrink_to_fit();
      energyvalues.shrink_to_fit();
      relativeerrors.shrink_to_fit();
      vovs.shrink_to_fit();
      exposedareas.shrink_to_fit();
    }

  public:
    virtual void use(unsigned int pPrimID, numeric_type value) override final
    {
      throw std::logic_error("Not implemented");
    }

    std::vector<internal_numeric_type> get_energy_values() override final
    {
      return energyvalues;
    }

    std::vector<size_t> get_cnts() override final
    {
      return cnts;
    }


    size_t get_cnts_sum() override final
    {
      return cntsSum;
    }

    std::vector<internal_numeric_type> get_relative_error() override final
    {
      return relativeerrors;
    }

    std::vector<internal_numeric_type> get_vov() override final
    {
      return vovs;
    }

    void set_exposed_areas(std::vector<internal_numeric_type>& areas) override final
    {
      exposedareas = areas;
    }

    std::vector<internal_numeric_type> get_exposed_areas() override final
    {
      return exposedareas;
    }

    void print(std::ostream& pOs) const override final
    {
      throw std::logic_error("Not implemented");
    }

    numeric_type get_total_energy_used() override final
    {
      throw std::logic_error("Not supported");
      // This is a reason to refactor: This class should not need to have this method.
    }

    size_t get_size() override final
    {
      return size;
    }

  private:
    i_hit_accumulator<numeric_type>& h1;
    i_hit_accumulator<numeric_type>& h2;
    size_t size;

    std::vector<size_t> cnts;
    std::vector<internal_numeric_type> energyvalues;
    std::vector<internal_numeric_type> relativeerrors;
    std::vector<internal_numeric_type> vovs;

    size_t cntsSum;

    std::vector<internal_numeric_type> exposedareas;
  };
}}
