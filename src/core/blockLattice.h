/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006-2008 Jonas Latt
 *                2008-2021 Mathias Krause
 *                2022      Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef BLOCK_LATTICE_H
#define BLOCK_LATTICE_H


#include "utilities/aliases.h"

#include "core/cell.h"
#include "core/stages.h"

#include "postProcessing.h"
#include "latticeStatistics.h"
#include "serializer.h"

#include "functors/analytical/analyticalF.h"

#include "geometry/blockGeometry.h"

#include "data.h"
#include "fieldArrayD.h"

#include "platform/cpu/sisd.h"

#ifdef PLATFORM_CPU_SIMD
#include "platform/cpu/simd.h"
#endif

#ifdef PLATFORM_GPU_CUDA
#include "platform/gpu/cuda.h"
#endif

#include "blockDynamicsMap.h"
#include "blockPostProcessorMap.h"

#include "communication/communicatable.h"

#include <memory>
#include <vector>
#include <map>
#include <functional>

namespace olb {


template<typename T, typename DESCRIPTOR> class SuperLattice;
template<typename T, typename DESCRIPTOR> struct Dynamics;
template<typename T, typename DESCRIPTOR, Platform PLATFORM> class ConcreteBlockLattice;

/// Curried ConcreteBlockLattice template for use in callUsingConcretePlatform
template<typename T, typename DESCRIPTOR>
struct ConcretizableBlockLattice {

using base_t = BlockLattice<T,DESCRIPTOR>;

template <Platform PLATFORM>
using type = ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>;

};

/// Platform-abstracted block lattice for external access and inter-block interaction
template<typename T, typename DESCRIPTOR>
class BlockLattice : public BlockStructure<DESCRIPTOR>
                   , public Serializable {
protected:
  const Platform _platform;

  bool _statisticsEnabled;
  LatticeStatistics<T>* _statistics;

public:
  BlockLattice(Vector<int,DESCRIPTOR::d> size, int padding, Platform platform);
  virtual ~BlockLattice();

  /// Execute the collide step on the non-overlapping block cells
  virtual void collide() = 0;
  /// Apply the streaming step to the entire block
  virtual void stream() = 0;

  /// Set processing context
  /**
   * This is currently used to trigger data transfers between host
   * and GPU data for Platform::GPU_CUDA.
   **/
  virtual void setProcessingContext(ProcessingContext) = 0;

  /// Returns pointers to host-side population locations of iCell
  /**
   * Used only for high-level olb::(Const)Cell interface
   **/
  virtual Vector<T*,DESCRIPTOR::q> getPopulationPointers(CellID iCell) = 0;

  /// Return platform used to process lattice
  Platform getPlatform() const {
    return _platform;
  }

  /// Return whether FIELD_TYPE is available / has been allocated
  template<typename FIELD_TYPE>
  bool hasData()
  {
    return callUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      _platform,
      this,
      [&](const auto* lattice) -> bool {
        return lattice->template hasData<FIELD_TYPE>();
      });
  }
  /// Return abstract interface for concrete FIELD_TYPE data
  template<typename FIELD_TYPE>
  const auto& getData(FIELD_TYPE field = FIELD_TYPE{}) const
  {
    return *callUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      _platform,
      this,
      [&](const auto* lattice) -> const auto* {
        return &(lattice->template getData<FIELD_TYPE>().asAbstract());
      });
  }

  /// Return abstract interface for concrete FIELD_TYPE data
  template<typename FIELD_TYPE>
  auto& getData(FIELD_TYPE field = FIELD_TYPE{})
  {
    return *callUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      _platform,
      this,
      [&](auto* lattice) -> auto* {
        return &(lattice->template getData<FIELD_TYPE>().asAbstract());
      });
  }

  /// Return abstract interface for FIELD array
  template<typename FIELD>
  const auto& getField(FIELD field = FIELD{}) const
  {
    return getData(Array<FIELD>{});
  }

  /// Return abstract interface for FIELD array
  template<typename FIELD>
  auto& getField(FIELD field = FIELD{})
  {
    return getData(Array<FIELD>{});
  }

  /// Get Cell interface for index iCell
  Cell<T,DESCRIPTOR> get(CellID iCell)
  {
    return Cell<T,DESCRIPTOR>(*this, iCell);
  }
  /// Get ConstCell interface for index iCell
  ConstCell<T,DESCRIPTOR> get(CellID iCell) const
  {
    return ConstCell<T,DESCRIPTOR>(*this, iCell);
  }

  /// Get Cell interface for location loc
  Cell<T,DESCRIPTOR> get(LatticeR<DESCRIPTOR::d> loc)
  {
    return get(this->getCellId(loc));
  }
  /// Get ConstCell interface for location loc
  ConstCell<T,DESCRIPTOR> get(LatticeR<DESCRIPTOR::d> loc) const
  {
    return get(this->getCellId(loc));
  }

  /// Get Cell interface for componentwise location latticeR
  template <typename... R>
  std::enable_if_t<sizeof...(R) == DESCRIPTOR::d, Cell<T,DESCRIPTOR>>
  get(R... latticeR)
  {
    return get(this->getCellId(latticeR...));
  }
  /// Get ConstCell interface for componentwise location latticeR
  template <typename... R>
  std::enable_if_t<sizeof...(R) == DESCRIPTOR::d, ConstCell<T,DESCRIPTOR>>
  get(R... latticeR) const
  {
    return get(this->getCellId(latticeR...));
  }

  /// Initialize the lattice cells to become ready for simulation
  void initialize();

  bool statisticsEnabled() const {
    return _statisticsEnabled;
  }
  void setStatisticsEnabled(bool state) {
    _statisticsEnabled = state;
  }

  /// Set dynamics at iCell to promised dynamics
  virtual void setDynamics(CellID iCell, DynamicsPromise<T,DESCRIPTOR>&&) = 0;
  // Legacy dynamics setter, to be deprecated
  virtual void setDynamics(CellID iCell, Dynamics<T,DESCRIPTOR>*) = 0;

  /// Return pointer to dynamics yielded by promise
  virtual Dynamics<T,DESCRIPTOR>* getDynamics(DynamicsPromise<T,DESCRIPTOR>&&) = 0;
  /// Return pointer to DYNAMICS (legacy)
  template<typename DYNAMICS>
  Dynamics<T,DESCRIPTOR>* getDynamics() {
    return getDynamics(DynamicsPromise(meta::id<DYNAMICS>{}));
  }

  /// Return pointer to dynamics at iCell
  virtual Dynamics<T,DESCRIPTOR>* getDynamics(CellID iCell) = 0;
  /// Return pointer to dynamics assigned to latticeR
  template <typename... R>
  std::enable_if_t<sizeof...(R) == DESCRIPTOR::d, Dynamics<T,DESCRIPTOR>*>
  getDynamics(R... latticeR)
  {
    return getDynamics(this->getCellId(latticeR...));
  }

  /// Assign promised DYNAMICS to latticeR
  void defineDynamics(LatticeR<DESCRIPTOR::d> latticeR,
                      DynamicsPromise<T,DESCRIPTOR>&& promise)
  {
    setDynamics(this->getCellId(latticeR),
                std::forward<DynamicsPromise<T,DESCRIPTOR>&&>(promise));
  }
  /// Assign DYNAMICS to latticeR
  template <template<typename,typename> typename DYNAMICS>
  void defineDynamics(LatticeR<DESCRIPTOR::d> latticeR)
  {
    setDynamics(this->getCellId(latticeR),
                DynamicsPromise(meta::id<DYNAMICS<T,DESCRIPTOR>>{}));
  }
  /// Assign dynamics to latticeR via pointer (legacy)
  void defineDynamics(LatticeR<DESCRIPTOR::d> latticeR,
                      Dynamics<T,DESCRIPTOR>* dynamics)
  {
    setDynamics(this->getCellId(latticeR), dynamics);
  }

  /// Define the dynamics on a domain described by an indicator
  void defineDynamics(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                      Dynamics<T,DESCRIPTOR>* dynamics);
  /// Define DYNAMICS on a domain described by an indicator
  template <typename DYNAMICS>
  void defineDynamics(BlockIndicatorF<T,DESCRIPTOR::d>& indicator);
  /// Define the same dynamics on all cells
  void defineDynamics(Dynamics<T,DESCRIPTOR>* dynamics);

  /// Set value of parameter FIELD for any dynamics that provide it
  /**
   * Most common way of defining parameters. E.g. to set the relaxation time
   * `descriptors::OMEGA` for all dynamics to `0.6`:
   * \code{.cpp}
   * blockLattice.template setParameter<descriptors::OMEGA>(0.6);
   * \endcode
   **/
  template <typename FIELD>
  void setParameter(FieldD<T,DESCRIPTOR,FIELD> value) {
    callUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      _platform,
      this,
      [&](auto* lattice) {
        lattice->template setParameter<FIELD>(value);
      });
  }

  template <typename PARAMETER, typename FIELD>
  void setParameter(AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& fieldArray) {
    callUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      _platform,
      this,
      [&](auto* lattice) {
        lattice->template setParameter<PARAMETER>(fieldArray);
      });
  }
  template <typename PARAMETER, Platform PLATFORM, typename FIELD>
  void setParameter(FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>& fieldArray) {
    callUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      _platform,
      this,
      [&](auto* lattice) {
        lattice->template setParameter<PARAMETER>(fieldArray);
      });
  }

  /// Returns true if stage contains post processor
  virtual bool hasPostProcessor(std::type_index stage,
                                PostProcessorPromise<T,DESCRIPTOR>&& promise) = 0;
  /// Schedule post processor for application to latticeR in stage
  virtual void addPostProcessor(std::type_index stage,
                                LatticeR<DESCRIPTOR::d> latticeR,
                                PostProcessorPromise<T,DESCRIPTOR>&& promise) = 0;
  /// Schedule post processor for application to entire block in stage
  virtual void addPostProcessor(std::type_index stage,
                                PostProcessorPromise<T,DESCRIPTOR>&& promise) = 0;
  /// Schedule post processor for application to indicated cells in stage
  virtual void addPostProcessor(std::type_index stage,
                                BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                PostProcessorPromise<T,DESCRIPTOR>&& promise) = 0;

  /// Schedule legacy post processor for application in stage
  virtual void addPostProcessor(std::type_index stage, PostProcessor<T,DESCRIPTOR>* postProcessor) = 0;
  /// Schedule legacy post processor in stage
  virtual void addPostProcessor(std::type_index stage, const PostProcessorGenerator<T,DESCRIPTOR>& ppGen) = 0;
  /// Schedule legacy post processor for application to indicated cells in stage
  virtual void addPostProcessor(std::type_index stage,
                                BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                                const PostProcessorGenerator<T,DESCRIPTOR>& ppGen) = 0;

  /// Schedule legacy post processor for application in STAGE
  template <typename STAGE=stage::PostStream>
  void addPostProcessor(PostProcessor<T,DESCRIPTOR>* postProcessor) {
    addPostProcessor(typeid(STAGE), postProcessor);
  }
  /// Schedule legacy post processor for application in STAGE
  template <typename STAGE=stage::PostStream>
  void addPostProcessor(const PostProcessorGenerator<T,DESCRIPTOR>& ppGen) {
    addPostProcessor(typeid(STAGE), ppGen);
  }
  /// Schedule legacy post processor for application to indicated cells in STAGE
  template <typename STAGE=stage::PostStream>
  void addPostProcessor(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                        const PostProcessorGenerator<T,DESCRIPTOR>& ppGen) {
    addPostProcessor(typeid(STAGE), indicator, ppGen);
  }

  /// Execute post processors of stage
  virtual void postProcess(std::type_index stage = typeid(stage::PostStream)) = 0;
  /// Execute post processors of STAGE
  template <typename STAGE>
  void postProcess() {
    postProcess(typeid(STAGE));
  }

  virtual bool hasCommunicatable(std::type_index) const = 0;
  virtual Communicatable& getCommunicatable(std::type_index) = 0;

  /// Define a field on a domain described by an indicator
  /**
   * \param indicator     Block indicator describing the target domain
   * \param field         Analytical functor (global)
   **/
  template <typename FIELD>
  void defineField(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                   AnalyticalF<DESCRIPTOR::d,T,T>& field);
  /// Define a field on a domain described by an indicator
  /**
   * \param indicator     Block indicator describing the target domain
   * \param field         Block functor
   **/
  template <typename FIELD>
  void defineField(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                   BlockF<T,DESCRIPTOR::d>& field);
  /// Define a field on a domain described by an analytical indicator
  /**
   * \param indicatorF Domain indicator to be reduced to BlockIndicatorFfromIndicatorF3D
   **/
  template <typename FIELD>
  void defineField(BlockGeometry<T,DESCRIPTOR::d>& blockGeometry,
                   IndicatorF<T,DESCRIPTOR::d>& indicatorF,
                   AnalyticalF<DESCRIPTOR::d,T,T>& field);


  /// Define rho on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   **/
  void defineRho(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                 AnalyticalF<DESCRIPTOR::d,T,T>& rho);

  /// Define u on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param u         Analytical functor (global)
   **/
  void defineU(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
               AnalyticalF<DESCRIPTOR::d,T,T>& u);

  /// Define rho and u on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  void defineRhoU(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                  AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u);

  /// Define a population on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param Pop       Analytical functor (global), target dimension DESCRIPTOR::q
   **/
  void definePopulations(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                         AnalyticalF<DESCRIPTOR::d,T,T>& Pop);
  /**
   * \param indicator Block indicator describing the target domain
   * \param Pop       Block functor, target dimension DESCRIPTOR::q
   **/
  void definePopulations(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                         BlockF<T,DESCRIPTOR::d>& Pop);

  /// Initialize by equilibrium on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  void iniEquilibrium(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u);
  void iniEquilibrium(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, BlockF<T,DESCRIPTOR::d>& u);

  /// Initialize by non- and equilibrium on a domain described by an indicator
  /**
   * \param indicator Block indicator describing the target domain
   * \param rho       Analytical functor (global)
   * \param u         Analytical functor (global)
   **/
  void iniRegularized(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                      AnalyticalF<DESCRIPTOR::d,T,T>& rho, AnalyticalF<DESCRIPTOR::d,T,T>& u, AnalyticalF<DESCRIPTOR::d,T,T>& pi);
  /// Subtract the given offset from all densities
  void stripeOffDensityOffset(T offset);

  /// Add a non-local post-processing step (legacy)
  void addLatticeCoupling(LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                          std::vector<BlockStructureD<DESCRIPTOR::d>*> partners);
  void addLatticeCoupling(BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                          LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                          std::vector<BlockStructureD<DESCRIPTOR::d>*> partners);
  /// Execute couplings steps (legacy)
  void executeCoupling();

  /// Return a handle to the LatticeStatistics object
  LatticeStatistics<T>& getStatistics();
  /// Return a constant handle to the LatticeStatistics object
  LatticeStatistics<T> const& getStatistics() const;

};

/// Implementation of BlockLattice on a concrete PLATFORM
template<typename T, typename DESCRIPTOR, Platform PLATFORM=Platform::CPU_SISD>
class ConcreteBlockLattice final : public BlockLattice<T,DESCRIPTOR> {
private:
  /// Field data
  Data<T,DESCRIPTOR,PLATFORM> _data;
  /// Index of DESCRIPTOR-declared field arrays
  utilities::FixedTypeIndexedMap<typename DESCRIPTOR::fields_t, ColumnVectorBase*> _descriptorFields;
  /// Pointers to Communicatable-casted FieldArrayD instances for overlap communication
  std::map<std::type_index, std::unique_ptr<Communicatable>> _communicatables;
  /// Assignments of dynamics instances to cell indices
  BlockDynamicsMap<T,DESCRIPTOR,PLATFORM> _dynamicsMap;
  /// Optional custom callable replacing default collision application
  std::optional<std::function<void(ConcreteBlockLattice&)>> _customCollisionO;
  /// Map of post processor stages
  std::map<std::type_index,
           std::map<int,
                    BlockPostProcessorMap<T,DESCRIPTOR,PLATFORM>,
                    std::less<int>>
           > _postProcessors;

public:
  using value_t = T;
  using descriptor_t = DESCRIPTOR;

#ifdef PLATFORM_CPU_SIMD
  static_assert(PLATFORM != Platform::CPU_SIMD || std::is_same_v<T,double> || std::is_same_v<T,float>,
                "SIMD blocks must use either single or double precision as fundamental type");
#endif

  static constexpr Platform platform = PLATFORM;

  static_assert(DESCRIPTOR::template provides<descriptors::POPULATION>(),
                "Lattice DESCRIPTOR must provide POPULATION field");

  ConcreteBlockLattice(Vector<int,DESCRIPTOR::d> size, int padding=0);

  void setProcessingContext(ProcessingContext context) override {
    _data.setProcessingContext(context);
  }

  template<typename FIELD_TYPE>
  bool hasData() const;
  template<typename FIELD_TYPE>
  const auto& getData() const;
  template<typename FIELD_TYPE>
  auto& getData();

  template<typename FIELD>
  const auto& getField(FIELD field = FIELD()) const;
  template<typename FIELD>
  auto& getField(FIELD field = FIELD());

  bool hasCommunicatable(std::type_index field) const override {
    return _communicatables.find(field) != _communicatables.end();
  }
  Communicatable& getCommunicatable(std::type_index field) override {
    return *_communicatables.at(field).get();
  }

  /// Apply collision step of non-overlap interior
  void collide() override;
  /// Perform propagation step on the whole block
  /**
   * Rotates the cyclic arrays storing the POPULATION field
   * to perform implicit propagation using the PS pattern.
   *
   * - Kummerländer, A., Dorn, M., Frank, M., and Krause, M. J. Implicit
   *   Propagation of Directly Addressed Grids in Lattice Boltzmann Methods.
   *   DOI: 10.13140/RG.2.2.35085.87523
   **/
  void stream() override;

  /// Replace default collision logic of BlockDynamicsMap
  /**
   * May be used to inject domain knowledge for improving performance
   * by e.g. reducing the need to use virtual dispatching for non-dominant
   * dynamics. This is necessarily platform specific.
   **/
  void setCollisionO(std::function<void(ConcreteBlockLattice&)>&& op) {
    _customCollisionO = op;
  }

  BlockDynamicsMap<T,DESCRIPTOR,PLATFORM>& getDynamicsMap() {
    return _dynamicsMap;
  }

  Dynamics<T,DESCRIPTOR>* getDynamics(DynamicsPromise<T,DESCRIPTOR>&& promise) override {
    return _dynamicsMap.get(std::forward<decltype(promise)>(promise));
  }

  /// Get reference to dynamics of cell by index
  Dynamics<T,DESCRIPTOR>* getDynamics(CellID iCell) override
  {
    return _dynamicsMap.get(iCell);
  }

  void setDynamics(CellID iCell, DynamicsPromise<T,DESCRIPTOR>&& promise) override
  {
    _dynamicsMap.set(iCell, std::forward<decltype(promise)>(promise));
    auto cell = this->get(iCell);
    getDynamics(iCell)->initialize(cell);
  }
  void setDynamics(CellID iCell, Dynamics<T,DESCRIPTOR>* dynamics) override
  {
    _dynamicsMap.set(iCell, dynamics);
  }

  template <typename FIELD>
  void setParameter(FieldD<T,DESCRIPTOR,FIELD> value)
  {
    _data.template forEachCastable<AbstractedConcreteParameters<T,DESCRIPTOR>>([&](auto* parameters) {
      auto& params = parameters->asAbstract();
      if (params.template provides<FIELD>()) {
        params.template set<FIELD>(value);
        parameters->setProcessingContext(ProcessingContext::Simulation);
      }
    });
  }

  template <typename PARAMETER, Platform _PLATFORM, typename FIELD>
  void setParameter(FieldArrayD<T,DESCRIPTOR,_PLATFORM,FIELD>& fieldArray)
  {
    if constexpr (PLATFORM == Platform::GPU_CUDA) {
      static_assert(PLATFORM == _PLATFORM, "FieldArrayD must be available on PLATFORM");
    }
    FieldD<T,DESCRIPTOR,PARAMETER> fieldArrayPointers;
    for (unsigned iD=0; iD < fieldArray.d; ++iD) {
      if constexpr (PLATFORM == Platform::GPU_CUDA) {
        fieldArrayPointers[iD] = fieldArray[iD].deviceData();
      } else {
        fieldArrayPointers[iD] = fieldArray[iD].data();
      }
    }
    setParameter<PARAMETER>(std::move(fieldArrayPointers));
  }

  template <typename PARAMETER, typename FIELD>
  void setParameter(AbstractFieldArrayD<T,DESCRIPTOR,FIELD>& abstractFieldArray)
  {
    auto& fieldArray = dynamic_cast<FieldArrayD<T,DESCRIPTOR,PLATFORM,FIELD>&>(abstractFieldArray);
    setParameter<PARAMETER>(fieldArray);
  }

  bool hasPostProcessor(std::type_index stage,
                        PostProcessorPromise<T,DESCRIPTOR>&& promise) override
  {
    auto [postProcessorsOfPriority, _] = _postProcessors[stage].try_emplace(promise.priority(), this);
    return std::get<1>(*postProcessorsOfPriority).contains(
      std::forward<decltype(promise)>(promise));
  }

  void addPostProcessor(std::type_index stage,
                        LatticeR<DESCRIPTOR::d> latticeR,
                        PostProcessorPromise<T,DESCRIPTOR>&& promise) override
  {
    if (!this->isPadding(latticeR)) {
      auto [postProcessorsOfPriority, _] = _postProcessors[stage].try_emplace(promise.priority(), this);
      std::get<1>(*postProcessorsOfPriority).add(this->getCellId(latticeR),
                                                 std::forward<decltype(promise)>(promise));
    }
  }

  void addPostProcessor(std::type_index stage,
                        BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                        PostProcessorPromise<T,DESCRIPTOR>&& promise) override
  {
    this->forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
      if (indicator(loc)) {
        addPostProcessor(stage, loc, std::forward<decltype(promise)>(promise));
      }
    });
  }

  void addPostProcessor(std::type_index stage,
                        PostProcessorPromise<T,DESCRIPTOR>&& promise) override
  {
    if (promise.scope() == OperatorScope::PerBlock) {
      auto [postProcessorsOfPriority, _] = _postProcessors[stage].try_emplace(promise.priority(), this);
      std::get<1>(*postProcessorsOfPriority).add(std::forward<decltype(promise)>(promise));
    } else {
      this->forCoreSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
        addPostProcessor(stage, loc, std::forward<decltype(promise)>(promise));
      });
    }
  }

  void addPostProcessor(std::type_index stage, PostProcessor<T,DESCRIPTOR>* postProcessor) override;
  void addPostProcessor(std::type_index stage, const PostProcessorGenerator<T,DESCRIPTOR>& ppGen) override;
  void addPostProcessor(std::type_index stage,
                        BlockIndicatorF<T,DESCRIPTOR::d>& indicator,
                        const PostProcessorGenerator<T,DESCRIPTOR>& ppGen) override;

  /// Execute post processors of stage
  void postProcess(std::type_index stage) override;

  /// Return pointers to population values of cell index iCell
  /**
   * Performance optimization for access via virtual Cell
   **/
  Vector<T*,DESCRIPTOR::q> getPopulationPointers(CellID iCell) override {
    auto& pops = getField<descriptors::POPULATION>();
    return Vector<T*,DESCRIPTOR::q>([&](unsigned iPop) -> T* {
      return &pops[iPop][iCell];
    });
  }

  /// Return reference to Data's FieldTypeRegistry
  /**
   * Only of interest for implementing specific device support
   **/
  auto& getDataRegistry() {
    return _data.getRegistry();
  }

  /// Number of data blocks for the serializable interface
  std::size_t getNblock() const override;
  /// Binary size for the serializer
  std::size_t getSerializableSize() const override;
  /// Return a pointer to the memory of the current block and its size for the serializable interface
  bool* getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode) override;
  /// Reinit population structure after deserialization
  void postLoad() override;


};


/// Wrapper for a local heterogeneous block communication request
/**
 * Specialized for Platform::GPU_CUDA as SOURCE resp. TARGET
 **/
template <typename T, typename DESCRIPTOR, Platform SOURCE, Platform TARGET>
class HeterogeneousCopyTask;

template <typename T, typename DESCRIPTOR, Platform PLATFORM>
class ConcreteBlockCommunicator<ConcreteBlockLattice<T,DESCRIPTOR,PLATFORM>>
  final : public BlockCommunicator {
private:
  const int _iC;
#ifdef PARALLEL_MODE_MPI
  MPI_Comm _mpiCommunicator;
#endif

#ifdef PARALLEL_MODE_MPI
  class SendTask;
  class RecvTask;

  std::vector<std::unique_ptr<SendTask>> _sendTasks;
  std::vector<std::unique_ptr<RecvTask>> _recvTasks;
#endif

public:
  struct CopyTask;

  ConcreteBlockCommunicator(SuperLattice<T,DESCRIPTOR>& super,
                            LoadBalancer<T>& loadBalancer,
#ifdef PARALLEL_MODE_MPI
                            SuperCommunicationTagCoordinator<T>& tagCoordinator,
                            MPI_Comm comm,
#endif
                            int iC,
                            const BlockCommunicationNeighborhood<T,DESCRIPTOR::d>& neighborhood);
  ~ConcreteBlockCommunicator();

#ifdef PARALLEL_MODE_MPI
  void receive() override;
  void send() override;
  void unpack() override;
  void wait() override;
#else
  void copy() override;
#endif

private:
  class HomogeneousCopyTask;

  std::vector<std::unique_ptr<CopyTask>> _copyTasks;

};


}

#endif
