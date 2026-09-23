#include "filter_xpbd_warp.h"

#include <map>
#include <set>
#include <vcg/complex/algorithms/update/bounding.h>
#include <vcg/complex/algorithms/update/normal.h>
#include <vcg/complex/algorithms/update/topology.h>
#include <vector>

using namespace vcg;
using namespace vcg::tri;

FilterXPBD::FilterXPBD()
{
	// Modern MeshLab uses an initializer list for typeList
	typeList = {FP_XPBD_WARP};

	// Populate the QActions for the UI menu
	for (ActionIDType tt : types()) {
		actionList.push_back(new QAction(filterName(tt), this));
	}
}

QString FilterXPBD::pluginName() const
{
	return "FilterXPBD";
}

QString FilterXPBD::filterName(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_XPBD_WARP: return "XPBD Mesh Inflation & Relaxation";
	default: assert(0); return "";
	}
}

QString FilterXPBD::filterInfo(ActionIDType filterId) const
{
	switch (filterId) {
	case FP_XPBD_WARP:
		return "Applies Extended Position-Based Dynamics (XPBD) to inflate or relax the mesh while "
			   "preserving edge-length constraints.";
	default: assert(0); return "";
	}
}

FilterPlugin::FilterClass FilterXPBD::getClass(const QAction* a) const
{
	switch (ID(a)) {
	// This maps the filter to the "Smoothing, Fairing and Deformation" menu
	case FP_XPBD_WARP: return FilterPlugin::Smoothing;
	default: assert(0); return FilterPlugin::Smoothing;
	}
}

FilterPlugin::FilterArity FilterXPBD::filterArity(const QAction* a) const
{
	switch (ID(a)) {
	case FP_XPBD_WARP: return SINGLE_MESH;
	default: assert(0); return NONE;
	}
}

RichParameterList FilterXPBD::initParameterList(const QAction* action, const MeshModel& /*m*/)
{
	RichParameterList par;
	switch (ID(action)) {
	case FP_XPBD_WARP:
		par.addParam(RichInt(
			"iterations",
			20,
			"Solver Iterations",
			"Number of XPBD projection iterations per step. Higher is stiffer but slower."));
		par.addParam(RichFloat(
			"compliance",
			0.001f,
			"Edge Compliance",
			"Inverse stiffness of the mesh edges. 0 is fully rigid."));
		par.addParam(RichFloat(
			"pressure",
			0.05f,
			"Inflation Pressure",
			"Outward force applied along vertex normals."));
		par.addParam(RichBool(
			"use_selection",
			false,
			"Lock Unselected",
			"If true, only selected vertices will be deformed; others act as rigid anchors."));
		break;
	default: assert(0);
	}
	return par;
}

int FilterXPBD::getPreConditions(const QAction* action) const
{
	switch (ID(action)) {
	case FP_XPBD_WARP: return MeshModel::MM_NONE;
	default: assert(0); return 0;
	}
}

std::map<std::string, QVariant> FilterXPBD::applyFilter(
	const QAction*           action,
	const RichParameterList& par,
	MeshDocument&            md,
	unsigned int& /*postConditionMask*/,
	vcg::CallBackPos* cb)
{
	switch (ID(action)) {
	case FP_XPBD_WARP: {
		MeshModel* mm = md.mm();
		if (!mm)
			return std::map<std::string, QVariant>();

		// Explicitly demand normals so VCG doesn't crash on unallocated components
		mm->updateDataMask(MeshModel::MM_VERTNORMAL | MeshModel::MM_FACENORMAL);

		CMeshO& m = mm->cm;

		int iterations = par.getInt("iterations");
		if (iterations <= 0)
			iterations = 1; // Prevent divide-by-zero crashes

		float compliance   = par.getFloat("compliance");
		float pressure     = par.getFloat("pressure");
		bool  useSelection = par.getBool("use_selection");

		// dt2 is still needed for XPBD compliance calculation (alpha)
		float dt  = 0.016f;
		float dt2 = dt * dt;

		vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(m);
		vcg::tri::UpdateNormal<CMeshO>::PerVertexNormalized(m);

		// 1. Allocate exactly the physical capacity of the VCG container
		size_t                    max_verts = m.vert.size();
		std::vector<vcg::Point3f> predictedPos(max_verts);
		std::vector<float>        invMass(max_verts, -1.0f); // -1.0f flags invalid/deleted holes

		// Calculate a scale-invariant multiplier based on the mesh size
		float scaleFactor = m.bbox.Diag() / 100.0f;

		int idx = 0;
		for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi, ++idx) {
			if (vi->IsD())
				continue;

			float im     = (useSelection && !vi->IsS()) ? 0.0f : 1.0f;
			invMass[idx] = im; // Safe: idx is strictly bound to max_verts

			// Apply pressure directly as a percentage of the mesh's bounding box
			vcg::Point3f externalForce = vi->N() * (pressure * scaleFactor);
			predictedPos[idx]          = vi->P() + (externalForce * im);
		}

		// 2. Extract edges safely using exact memory offsets
		struct EdgeConstraint
		{
			int   v1;
			int   v2;
			float rest_len;
			bool  operator<(const EdgeConstraint& o) const
			{
				if (v1 != o.v1)
					return v1 < o.v1;
				return v2 < o.v2;
			}
		};

		std::set<EdgeConstraint> unique_edges;
		for (auto fi = m.face.begin(); fi != m.face.end(); ++fi) {
			if (fi->IsD())
				continue;

			for (int j = 0; j < 3; ++j) {
				auto vp1 = fi->V(j);
				auto vp2 = fi->V((j + 1) % 3);

				if (vp1->IsD() || vp2->IsD())
					continue;

				// vcg::tri::Index gets the exact integer position in the array
				int id1 = vcg::tri::Index(m, vp1);
				int id2 = vcg::tri::Index(m, vp2);

				// Unbreakable hardware bounds check against corrupted faces
				if (id1 < 0 || id1 >= max_verts || id2 < 0 || id2 >= max_verts)
					continue;

				if (id1 > id2)
					std::swap(id1, id2);

				float rlen = (vp1->P() - vp2->P()).Norm();
				unique_edges.insert({id1, id2, rlen});
			}
		}

		// 3. XPBD Solver Iterations
		float alpha = compliance / dt2;
		for (int iter = 0; iter < iterations; ++iter) {
			if (cb && (iter % 5 == 0)) {
				bool canceled = (*cb)(100 * iter / iterations, "Solving XPBD...");
				if (canceled)
					break;
			}

			// Iterate the set directly (no vector conversion needed)
			for (const auto& c : unique_edges) {
				float w1 = invMass[c.v1];
				float w2 = invMass[c.v2];

				if (w1 < 0.0f || w2 < 0.0f)
					continue; // Skip invalid vertex holes

				float wSum = w1 + w2;
				if (wSum == 0.0f)
					continue;

				vcg::Point3f& p1 = predictedPos[c.v1];
				vcg::Point3f& p2 = predictedPos[c.v2];

				vcg::Point3f dir = p1 - p2;
				float        len = dir.Norm();
				if (len < 1e-6f)
					continue;

				float        C          = len - c.rest_len;
				float        dLambda    = -C / (wSum + alpha);
				vcg::Point3f correction = (dir / len) * dLambda;

				p1 += correction * w1;
				p2 -= correction * w2;
			}
		}

		// 4. Finalize positions
		idx = 0;
		for (auto vi = m.vert.begin(); vi != m.vert.end(); ++vi, ++idx) {
			if (vi->IsD())
				continue;
			if (idx >= max_verts)
				break; // Ultimate bounds check

			if (invMass[idx] > 0.0f) {
				vi->P() = predictedPos[idx];
			}
		}

		// 5. Cleanup
		vcg::tri::UpdateNormal<CMeshO>::PerFaceNormalized(m);
		vcg::tri::UpdateNormal<CMeshO>::PerVertexNormalized(m); // Safer than AngleWeighted
		vcg::tri::UpdateBounding<CMeshO>::Box(m);

		break;
	}
	default: wrongActionCalled(action); break;
	}

	return std::map<std::string, QVariant>();
}

int FilterXPBD::postCondition(const QAction* action) const
{
	switch (ID(action)) {
	case FP_XPBD_WARP:
		// MM_VERTCOORD tells MeshLab to visually update the vertex positions in the viewport
		return MeshModel::MM_VERTCOORD | MeshModel::MM_VERTNORMAL | MeshModel::MM_FACENORMAL;
	default: assert(0); return 0;
	}
}

// Add the missing version check implementation
std::pair<std::string, bool> FilterXPBD::getMLVersion() const
{
	// MESHLAB_VERSION is automatically passed by CMake during build.
	// false = single precision
	return std::make_pair(std::string(std::to_string(MESHLAB_VERSION)), false);
}

// Essential for modern MeshLab plugins to export the class properly
MESHLAB_PLUGIN_NAME_EXPORTER(FilterXPBD)
