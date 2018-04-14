#ifndef TEMPLATE_INSTANCE_TYPEDEFS_H
#define TEMPLATE_INSTANCE_TYPEDEFS_H

#include "mdt.h"
#include "defines.h" // for Key
#include "NodePayload.h"
#include "DataInterface.h"

// have to do this because charmxi doesn't like occurrence
// of keyword 'typename' in dependent type declaration
// in .ci file.
#define SPH_VISITOR_TYPE_REMOTE_NODE_TYPE_MACRO typename SphVisitorType::RemoteNodeType
#define SPH_VISITOR_TYPE_REMOTE_PARTICLE_TYPE_MACRO typename SphVisitorType::RemoteParticleType

// home tree node type
typedef TreeBuildDataInterface::LocalNodeType MyLocalNodeType;

// tree piece type
typedef MultiphaseDistree::Piece<TreeBuildDataInterface> MyPieceType;

// tree manager type
typedef MultiphaseDistree::Manager<TreeBuildDataInterface> MyTreeManagerType;

// tree handle
typedef MultiphaseDistree::Handle<TreeBuildDataInterface> MyTreeHandleType;

// traversals
typedef MultiphaseDistree::LeafNodeTopDownTraversalManager<GravityTraversalDataInterface> MyGravityTraversalType;
typedef MultiphaseDistree::LeafNodeBottomUpTraversalManager<SphDensityTraversalDataInterface> MySphDensityTraversalType;
typedef MultiphaseDistree::LeafNodeTopDownTraversalManager<SphBallTraversalDataInterface> MySphBallTraversalType;

// cache 
typedef MultiphaseDistree::CacheManager<GravityTraversalDataInterface> MyGravityCacheManagerType;
typedef MultiphaseDistree::CacheManager<SphBallTraversalDataInterface> MySphBallCacheManagerType;

// cache client traversal container 
typedef MultiphaseDistree::CacheClientTraversals CacheClientsType;
typedef MultiphaseDistree::NodeCacheClientTraversals NodeCacheClientsType;
typedef MultiphaseDistree::LeafCacheClientTraversals LeafCacheClientsType;

// traversal manager handles
typedef MultiphaseDistree::TraversalHandle<MyGravityTraversalType> MyGravityTraversalHandleType;
typedef MultiphaseDistree::TraversalHandle<MySphDensityTraversalType> MySphDensityTraversalHandleType;
typedef MultiphaseDistree::TraversalHandle<MySphBallTraversalType> MySphBallTraversalHandleType;

#endif // TEMPLATE_INSTANCE_TYPEDEFS_H
