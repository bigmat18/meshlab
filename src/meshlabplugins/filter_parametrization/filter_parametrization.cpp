/****************************************************************************
* MeshLab                                                           o o     *
* A versatile mesh processing toolbox                             o     o   *
*                                                                _   O  _   *
* Copyright(C) 2005-2021                                           \/)\/    *
* Visual Computing Lab                                            /\/|      *
* ISTI - Italian National Research Council                           |      *
*                                                                    \      *
* All rights reserved.                                                      *
*                                                                           *
* This program is free software; you can redistribute it and/or modify      *
* it under the terms of the GNU General Public License as published by      *
* the Free Software Foundation; either version 2 of the License, or         *
* (at your option) any later version.                                       *
*                                                                           *
* This program is distributed in the hope that it will be useful,           *
* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
* GNU General Public License (http://www.gnu.org/licenses/gpl.txt)          *
* for more details.                                                         *
*                                                                           *
****************************************************************************/

#include "filter_parametrization.h"

#include <common/utilities/eigen_mesh_conversions.h>

#include <igl/boundary_loop.h>
#include <igl/harmonic.h>
#include <igl/lscm.h>
#include <igl/map_vertices_to_circle.h>
#include <vcg/complex/algorithms/update/texture.h>
#include <vcg/complex/algorithms/parametrization/uv_utils.h>

#include <vcg/complex/algorithms/parametrization/distortion.h>
#include <vcg/complex/algorithms/geodesic.h>
#include <vcg/complex/algorithms/curve_on_manifold.h>
#include <vcg/complex/algorithms/crease_cut.h>
#include<vcg/complex/algorithms/cut_tree.h>

using namespace vcg;

class CutEdge;
class CutFace;
class CutVertex;
struct CutUsedTypes : public UsedTypes<Use<CutVertex>::AsVertexType, Use<CutEdge>::AsEdgeType, Use<CutFace>::AsFaceType>{};

class CutVertex  : public Vertex< CutUsedTypes, vertex::Coord3f, vertex::Normal3f, vertex::Qualityf, vertex::Color4b, vertex::VEAdj, vertex::VFAdj,vertex::BitFlags  >{};
class CutEdge    : public Edge<   CutUsedTypes, edge::VertexRef, edge::VEAdj,     edge::EEAdj, edge::BitFlags> {};
class CutFace    : public Face  < CutUsedTypes, face::VertexRef,  face::Normal3f, face::Qualityf, face::Color4b, face::VFAdj, face::FFAdj, face::Mark, face::Color4b, face::BitFlags > {};
class CutMesh : public tri::TriMesh< std::vector<CutVertex>, std::vector<CutEdge>, std::vector<CutFace> >{};

FilterParametrizationPlugin::FilterParametrizationPlugin()
{
	typeList = { FP_HARMONIC_PARAM, FP_LEAST_SQUARES_PARAM, FP_MAX_DISTORTION_CUT, FP_TOPOLOGICAL_CUT};

	for(const ActionIDType& tt : typeList)
		actionList.push_back(new QAction(filterName(tt), this));
}

QString FilterParametrizationPlugin::pluginName() const
{
	return "FilterParametrization";
}

QString FilterParametrizationPlugin::vendor() const
{
	return "CNR-ISTI-VCLab";
}

QString FilterParametrizationPlugin::filterName(ActionIDType filterId) const
{
	switch(filterId) {
	case FP_HARMONIC_PARAM :
		return "Parametrization: Harmonic";
	case FP_LEAST_SQUARES_PARAM:
		return "Parametrization: LSCM";
	case FP_MAX_DISTORTION_CUT:
		return "Cut mesh from max distortion point";
	case FP_TOPOLOGICAL_CUT:
		return "Topological cut";
	default :
		assert(0);
		return "";
	}
}

QString FilterParametrizationPlugin::pythonFilterName(ActionIDType filter) const
{
	switch(filter) {
	case FP_HARMONIC_PARAM :
		return "compute_texcoord_parametrization_harmonic";
	case FP_LEAST_SQUARES_PARAM:
		return "compute_texcoord_parametrization_least_squares_conformal_maps";
	case FP_MAX_DISTORTION_CUT:
		return "compute_cut_from_max_distortion_point";
	case FP_TOPOLOGICAL_CUT:
		return "compute_topological_cut";
	default :
		assert(0);
		return "";
	}
}

QString FilterParametrizationPlugin::filterInfo(ActionIDType filterId) const
{
	QString commonDescription =
		"The resulting parametrization is saved in the per vertex texture coordinates.<br>"
		"The filter uses the original code provided in the "
		"<a href=\"https://libigl.github.io/\">libigl library</a>.<br>";
	switch(filterId) {
	case FP_HARMONIC_PARAM :
		return "Computes a single patch, fixed boundary harmonic parametrization of a mesh. The "
			   "filter requires that the input mesh has a single fixed boundary." +
			   commonDescription;
	case FP_LEAST_SQUARES_PARAM:
		return "Computes a least squares conformal maps (LSCM) parametrization of a mesh. " +
			   commonDescription;
	case FP_MAX_DISTORTION_CUT:
		return "Compute a cut on mesh from point with max distortion";
	case FP_TOPOLOGICAL_CUT:
		return "Compute a topological cut on a mesh";
	default :
		assert(0);
		return "Unknown Filter";
	}
}

FilterParametrizationPlugin::FilterClass FilterParametrizationPlugin::getClass(const QAction *a) const
{
	switch(ID(a)) {
	case FP_HARMONIC_PARAM :
	case FP_LEAST_SQUARES_PARAM:
		return FilterPlugin::Texture;
	case FP_MAX_DISTORTION_CUT:
		return FilterPlugin::Texture;
	case FP_TOPOLOGICAL_CUT:
		return FilterPlugin::Texture;
	default :
		assert(0);
		return FilterPlugin::Generic;
	}
}

FilterPlugin::FilterArity FilterParametrizationPlugin::filterArity(const QAction*) const
{
	return SINGLE_MESH;
}

int FilterParametrizationPlugin::getPreConditions(const QAction*) const
{
	return MeshModel::MM_VERTCOORD | MeshModel::MM_FACEVERT;
}

int FilterParametrizationPlugin::getRequirements(const QAction*)
{
	return MeshModel::MM_VERTTEXCOORD;
}

int FilterParametrizationPlugin::postCondition(const QAction*) const
{
	return MeshModel::MM_VERTTEXCOORD;
}

RichParameterList FilterParametrizationPlugin::initParameterList(const QAction *action, const MeshModel &)
{
	RichParameterList parlst;
	switch(ID(action)) {
	case FP_HARMONIC_PARAM :
		parlst.addParam(RichInt("harm_function", 1,"N-Harmonic Function", "1 denotes harmonic function, 2 biharmonic, 3 triharmonic, etc."));			
		parlst.addParam(RichBool("harm_wedge",true,"Per Wedge UV","If true it generates per wedge texture coordinates, otherwise it generate per-vertex texcoords."));
		parlst.addParam(RichBool("harm_uv_fit",true,"UV fit","If true it rescale the generate texture coords so that it lies in the [0..1]x[0..1] UV space."));
		break;
	case FP_LEAST_SQUARES_PARAM:
		parlst.addParam(RichBool("lscm_wedge",true,"Per Wedge UV","If true it generates per wedge texture coordinates, otherwise it generate per-vertex texcoords."));
		parlst.addParam(RichBool("lscm_uv_fit",true,"UV fit","If true it rescale the generate texture coords so that it lies in the [0..1]x[0..1] UV space."));
		break;
	case FP_MAX_DISTORTION_CUT:
		parlst.addParam(RichInt("distortion_fun", 1, "Type of distortion (1 Area Distortion, 2 Edge Distortion, 3 Angle Distortion)", "1 Area Distortion, 2 Edge Distortion, 3 Angle Distortion" ));
		break;
	case FP_TOPOLOGICAL_CUT:
		break;
	default :
		assert(0);
	}
	return parlst;
}

std::map<std::string, QVariant> FilterParametrizationPlugin::applyFilter(
		const QAction * action,
		const RichParameterList & par,
		MeshDocument &md,
		unsigned int& /*postConditionMask*/,
		vcg::CallBackPos *)
{
	switch(ID(action)) {
	case FP_HARMONIC_PARAM : {
		int f = par.getInt("harm_function");
		if (f < 1)
			throw MLException("Invalid N-Harmonic Function value. Must be >= 1.");

		EigenMatrixX3m v = meshlab::vertexMatrix(md.mm()->cm);
		Eigen::MatrixX3d verts = v.cast<double>();
		Eigen::MatrixX3i faces = meshlab::faceMatrix(md.mm()->cm);

		Eigen::MatrixXd V_uv, bnd_uv;
		Eigen::VectorXi bnd;

		igl::boundary_loop(faces,bnd);
		if (bnd.size() == 0)
			throw MLException(
				"Harmonic Parametrization can be applied only on meshes that have a boundary.");

		md.mm()->updateDataMask(MeshModel::MM_FACEFACETOPO);
		if(tri::Clean<CMeshO>::CountConnectedComponents(md.mm()->cm) > 1) 
			throw MLException(
				"Harmonic Parametrization can be applied only on meshes that "
				"have only one connected components");


		igl::map_vertices_to_circle(verts, bnd, bnd_uv);
		igl::harmonic(verts,faces,bnd,bnd_uv,1,V_uv);

		unsigned int i = 0;
		for (auto& v : md.mm()->cm.vert){
			v.T().u() =V_uv(i, 0);
			v.T().v() =V_uv(i, 1);
			i++;
		}
		
		// if requested it rescale the generate texture coords so that it lies in the [0..1]x[0..1] UV space
		if(par.getBool("harm_uv_fit"))
		{
			// use the code in the static function RegularizeTexArea in the class voronoiTexture
			vcg::tri::UV_Utils<CMeshO>::PerVertScaleToUnitSpace(md.mm()->cm);
		}
		

		if(par.getBool("harm_wedge"))
		{
			md.mm()->updateDataMask(MeshModel::MM_WEDGTEXCOORD);
			tri::UpdateTexture<CMeshO>::WedgeTexFromVertexTex(md.mm()->cm);			
		}
	
		break;
	}
	case FP_LEAST_SQUARES_PARAM : {
		EigenMatrixX3m v = meshlab::vertexMatrix(md.mm()->cm);
		Eigen::MatrixX3d verts = v.cast<double>();
		Eigen::MatrixX3i faces = meshlab::faceMatrix(md.mm()->cm);
		Eigen::VectorXi bnd, boundaryPoints(2, 1);

		Eigen::MatrixXd V_uv;

		igl::boundary_loop(faces,bnd);
		if (bnd.size() == 0)
			throw MLException(
				"Least Squares Conformal Maps Parametrization can be applied only on meshes that "
				"have a boundary.");
		
		for(auto& f : md.mm()->cm.face) {
			double area = vcg::DoubleArea(f);
			if(area == 0)
				throw MLException("Least Squares Conformal Maps Parametrization can be applied only "
								  "on meshes that haven't faces with area value equals to 0");
		}

		if(tri::Clean<CMeshO>::RemoveUnreferencedVertex(md.mm()->cm, false))
			throw MLException(
				"Least Squares Conformal Maps Parametrization can be applied only on meshes that "
				"have no unreference vertex");


		md.mm()->updateDataMask(MeshModel::MM_FACEFACETOPO);
		if(tri::Clean<CMeshO>::CountConnectedComponents(md.mm()->cm) > 1) 
			throw MLException(
					"Least Squares Conformal Maps Parametrization can be applied only on meshes that "
					"have only one connected components");
		

		boundaryPoints(0) = bnd(0);
		boundaryPoints(1) = bnd(bnd.size()/2);

		Eigen::MatrixXd bc(2,2);
		bc<<0,0,1,1;

		// LSCM parametrization
		bool ret = igl::lscm(verts,faces,boundaryPoints,bc,V_uv);
		if(!ret)
			throw MLException("Least Squares Conformal Maps Parametrization failed.");
		
		unsigned int i = 0;
		for (auto& v : md.mm()->cm.vert){
			v.T().u() =V_uv(i, 0);
			v.T().v() =V_uv(i, 1);
			i++;
		}
		// if requested it rescale the generate texture coords so that it lies in the [0..1]x[0..1] UV space
		if(par.getBool("lscm_uv_fit"))
		{
			// use the code in the static function RegularizeTexArea in the class voronoiTexture
			vcg::tri::UV_Utils<CMeshO>::PerVertScaleToUnitSpace(md.mm()->cm);
		}
		

		if(par.getBool("lscm_wedge"))
		{
			md.mm()->updateDataMask(MeshModel::MM_WEDGTEXCOORD);
			tri::UpdateTexture<CMeshO>::WedgeTexFromVertexTex(md.mm()->cm);			
		}
			
		break;
	}
	case FP_MAX_DISTORTION_CUT : {
		
		MeshModel *m = md.mm();
		m->updateDataMask(
			MeshModel::MM_WEDGTEXCOORD | 
			MeshModel::MM_VERTTEXCOORD |
			MeshModel::MM_FACEQUALITY | 
			MeshModel::MM_VERTQUALITY | 
			MeshModel::MM_VERTFACETOPO | 
			MeshModel::MM_FACEFACETOPO | 
			MeshModel::MM_FACEMARK |
			MeshModel::MM_VERTCOLOR
		);

		// Choose witch distorion type use to calculate path
		vcg::tri::Distortion<CMeshO, true>::DistType type;
		switch (par.getInt("distortion_fun"))
		{
			case 1: 
				type = vcg::tri::Distortion<CMeshO, true>::DistType::AreaDist; break;
			case 2: 
				type = vcg::tri::Distortion<CMeshO, true>::DistType::EdgeDist; break;
			case 3: 
				type = vcg::tri::Distortion<CMeshO, true>::DistType::AngleDist; break;
			default:
				throw MLException("Distortion function options must be: 1 Area Distortion, 2 Edge Distortion, 3 Angle Distortion");
		}

		// Calculate distorion for each face and vertex
		vcg::tri::Distortion<CMeshO, true>::SetQasDistorsion(m->cm, type);
		vcg::tri::UpdateFlags<CMeshO>::VertexBorderFromNone(m->cm);
		vcg::tri::UpdateQuality<CMeshO>::VertexNormalize(m->cm);

		vcg::tri::UpdateColor<CMeshO>::PerVertexConstant(m->cm, Color4b(0, 0, 0, 0));

		// search the vertex with max distortion value that there is not boundary
		float maxDistortion = 0;
		int vertexIndex = 0;
		std::vector<CMeshO::VertexPointer> bnd;

		for (auto vi = m->cm.vert.begin(); vi != m->cm.vert.end(); ++vi) {
			vi->C() = vcg::Color4b(static_cast<unsigned char>(vi->Q() * 255), 0, 0, 255);
			if(vi->IsB()) {
				bnd.push_back(&(*vi));	
			} else if(vi->Q() > maxDistortion) {
				maxDistortion = vi->Q();
				vertexIndex = vi->Index();
			}
		}

		std::cout << vertexIndex << std::endl;

		if (bnd.empty())
			throw MLException(
				"Cut can be applied only on meshes that have a boundary.");

		m->cm.vert[vertexIndex].C() = vcg::Color4b(0, 255, 0, 255);

		// Calculate for each vertex of m the distance of border. Inside parents is stored the path fo boundary from each vertex.
	    CMeshO::PerVertexAttributeHandle<CMeshO::VertexPointer> parents;
		parents = vcg::tri::Allocator<CMeshO>::GetPerVertexAttribute<CMeshO::VertexPointer>(m->cm);

		vcg::tri::EuclideanDistance<CMeshO> dd;
    	tri::UpdateQuality<CMeshO>::VertexConstant(m->cm,0);
		vcg::tri::Geodesic<CMeshO>::Compute(m->cm, bnd, dd, std::numeric_limits<CMeshO::ScalarType>::max(), nullptr, nullptr, &parents);

		// store each edge route to boundary from max distortion vertex
		CMeshO polyline;
		while (parents[vertexIndex]->Index() != vertexIndex) {
			// std::cout << vertexIndex << " " << parents[vertexIndex]->Index() << std::endl;
			vcg::tri::Allocator<CMeshO>::AddEdge(polyline,m->cm.vert[vertexIndex].P(), parents[vertexIndex]->P());
			vertexIndex = parents[vertexIndex]->Index();
		}; 
		
		// generate a cut for polyline
		vcg::tri::CoM<CMeshO> cc(m->cm);
		cc.Init();
		bool ret = cc.TagFaceEdgeSelWithPolyLine(polyline);
		if(ret) {
			vcg::tri::UpdateTopology<CMeshO>::FaceFace(m->cm);
			vcg::tri::CutMeshAlongSelectedFaceEdges<CMeshO>(m->cm);
		}

		break;
	}
	case FP_TOPOLOGICAL_CUT : {
		MeshModel *m = md.mm();
		m->updateDataMask(
			MeshModel::MM_WEDGTEXCOORD | 
			MeshModel::MM_VERTTEXCOORD |
			MeshModel::MM_VERTFACETOPO | 
			MeshModel::MM_FACEFACETOPO | 
			MeshModel::MM_FACEMARK
		);

		CutMesh polyline, cm;
		vcg::tri::Append<CutMesh,CMeshO>::MeshCopy(cm,m->cm);

		srand(time(nullptr));
		vcg::tri::CutTree<CutMesh> ct(cm);
		ct.Build(polyline, rand() % cm.fn);

		vcg::tri::CoM<CutMesh> cc(cm);
		cc.Init();
		if(cc.TagFaceEdgeSelWithPolyLine(polyline)) {
			vcg::tri::UpdateTopology<CutMesh>::FaceFace(cm);
			vcg::tri::CutMeshAlongSelectedFaceEdges<CutMesh>(cm);
		}
		vcg::tri::Append<CMeshO,CutMesh>::MeshCopy(m->cm,cm);
		break;
	}
	default :
		wrongActionCalled(action);
	}
	return std::map<std::string, QVariant>();
}

MESHLAB_PLUGIN_NAME_EXPORTER(FilterSamplePlugin)
