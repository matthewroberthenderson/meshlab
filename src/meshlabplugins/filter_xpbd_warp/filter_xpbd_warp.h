#ifndef FILTER_XPBD_WARP_H
#define FILTER_XPBD_WARP_H

#include <QObject>
#include <common/plugins/interfaces/filter_plugin.h>
#include <string>  // Needed for std::string
#include <utility> // Needed for std::pair

class FilterXPBD : public QObject, public FilterPlugin
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID "vcg.meshlab.FilterPlugin/1.0")
	Q_INTERFACES(FilterPlugin)

public:
	enum { FP_XPBD_WARP };

	FilterXPBD();

	QString pluginName() const override;
	QString filterName(ActionIDType filterId) const override;
	QString filterInfo(ActionIDType filterId) const override;

	FilterClass getClass(const QAction* a) const override;
	FilterArity filterArity(const QAction* a) const override;

	RichParameterList initParameterList(const QAction* action, const MeshModel& m) override;
	int               getPreConditions(const QAction* action) const override;

	std::map<std::string, QVariant> applyFilter(
		const QAction*           action,
		const RichParameterList& parameters,
		MeshDocument&            md,
		unsigned int&            postConditionMask,
		vcg::CallBackPos*        cb) override;

	int postCondition(const QAction* action) const override;

	// REQUIRED by modern MeshLabPlugin base class
	std::pair<std::string, bool> getMLVersion() const override;
};

#endif // FILTER_XPBD_WARP_H
