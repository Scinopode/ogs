/**
 * \file VtkVisTabWidget.cpp
 * 18/2/2010 LB Initial implementation
 *
 * Implementation of VtkVisTabWidget
 */

// ** INCLUDES **
#include "VtkVisTabWidget.h"
#include "VtkVisPipelineItem.h"
#include "VtkCompositeColorByHeightFilter.h"
#include "VtkColorByHeightFilter.h"

#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformFilter.h>

#include "ColorTableModel.h"
#include "ColorTableView.h"

#include "VtkAlgorithmProperties.h"
#include "VtkAlgorithmPropertyLineEdit.h"
#include "VtkCompositeFilter.h"
#include "VtkAlgorithmPropertyCheckbox.h"
#include "VtkAlgorithmPropertyVectorEdit.h"

#include <vtkPointData.h>
#include <vtkCellData.h>

VtkVisTabWidget::VtkVisTabWidget( QWidget* parent /*= 0*/ )
: QWidget(parent)
{
	setupUi(this);

	connect(this->vtkVisPipelineView, SIGNAL(itemSelected(VtkVisPipelineItem*)),
		this, SLOT(setActiveItem(VtkVisPipelineItem*)));

	connect(this->activeScalarComboBox, SIGNAL(currentIndexChanged(int)),
		this, SLOT(SetActiveAttributeOnItem(int)));

}

void VtkVisTabWidget::setActiveItem( VtkVisPipelineItem* item )
{
	if (item)
	{
		_item = item;

		vtkActor* actor = dynamic_cast<vtkActor*>(_item->actor());
		if (actor)
		{
			actorPropertiesGroupBox->setEnabled(true);
			vtkProperty* vtkProps = actor->GetProperty();
			diffuseColorPickerButton->setColor(vtkProps->GetDiffuseColor());
			visibleEdgesCheckBox->setChecked(vtkProps->GetEdgeVisibility());
			edgeColorPickerButton->setColor(vtkProps->GetEdgeColor());
			opacitySlider->setValue((int)(vtkProps->GetOpacity() * 100.0));
			vtkTransform* transform = 
				static_cast<vtkTransform*>(_item->transformFilter()->GetTransform());
			double scale[3];
			transform->GetScale(scale);
			scaleZ->setText(QString::number(scale[2]));

			vtkDataSet* dataSet = vtkDataSet::SafeDownCast(item->algorithm()->GetOutputDataObject(0));
			QStringList dataSetAttributesList;
			if (dataSet)
			{
				vtkPointData* pointData = dataSet->GetPointData();
				std::cout << "  #point data arrays: " << pointData->GetNumberOfArrays() << std::endl;
				for (int i = 0; i < pointData->GetNumberOfArrays(); i++)
				{
					std::cout << "    Name: " << pointData->GetArrayName(i) << std::endl;
					dataSetAttributesList.push_back(pointData->GetArrayName(i));
				}

				vtkCellData* cellData = dataSet->GetCellData();
				std::cout << "  #cell data arrays: " << cellData->GetNumberOfArrays() << std::endl;
				for (int i = 0; i < cellData->GetNumberOfArrays(); i++)
				{
					std::cout << "    Name: " << cellData->GetArrayName(i) << std::endl;
					dataSetAttributesList.push_back(cellData->GetArrayName(i));
				}

			}
			this->activeScalarComboBox->clear();
			this->activeScalarComboBox->insertItems(0, dataSetAttributesList);
		}
		else
			actorPropertiesGroupBox->setEnabled(false);

		QFormLayout* layout = static_cast<QFormLayout*>(this->scrollAreaWidgetContents->layout());
		while(layout->count())
				delete layout->takeAt(0)->widget();

		QMap<QString, QVariant>* propMap = NULL;
		QMap<QString, QList<QVariant> >* propVecMap = NULL;
		VtkAlgorithmProperties* algProps = NULL;
		if (item->compositeFilter())
		{
			algProps = item->compositeFilter();
			propMap = item->compositeFilter()->GetAlgorithmUserProperties();
			propVecMap = item->compositeFilter()->GetAlgorithmUserVectorProperties();
		}
		else
		{
			algProps = dynamic_cast<VtkAlgorithmProperties*>(item->algorithm());
			if (algProps)
			{
				propMap = algProps->GetAlgorithmUserProperties();
				propVecMap = algProps->GetAlgorithmUserVectorProperties();
			}
		}

		if (propMap && algProps)
		{
			QMapIterator<QString, QVariant> i(*propMap);
			while (i.hasNext())
			{
				i.next();
				QString key = i.key();
				QVariant value = i.value();

				VtkAlgorithmPropertyLineEdit* lineEdit;
				VtkAlgorithmPropertyCheckbox* checkbox;
				switch (value.type())
				{
					case QVariant::Double:
						lineEdit = new VtkAlgorithmPropertyLineEdit(QString::number(value.toDouble()), key, QVariant::Double, algProps);
						connect(lineEdit, SIGNAL(editingFinished()), this, SIGNAL(requestViewUpdate()));
						layout->addRow(key, lineEdit);
						break;

					case QVariant::Int:
						lineEdit = new VtkAlgorithmPropertyLineEdit(QString::number(value.toInt()), key, QVariant::Int, algProps);
						connect(lineEdit, SIGNAL(editingFinished()), this, SIGNAL(requestViewUpdate()));
						layout->addRow(key, lineEdit);
						break;

					case QVariant::Bool:
						checkbox = new VtkAlgorithmPropertyCheckbox(value.toBool(), key, algProps);
						connect(checkbox, SIGNAL(stateChanged(int)), this, SIGNAL(requestViewUpdate()));
						layout->addRow(key, checkbox);
						break;


					default:
						break;
				}
			}
		}

		if (propVecMap && algProps)
		{
			QMapIterator<QString, QList<QVariant> > i(*propVecMap);
			while (i.hasNext())
			{
				i.next();
				QString key = i.key();
				QList<QVariant> values = i.value();

				VtkAlgorithmPropertyVectorEdit* vectorEdit;
				if (values.size() > 0)
				{
					QList<QString> valuesAsString;
					foreach (QVariant variant, values)
						valuesAsString.push_back(variant.toString());

					vectorEdit = new VtkAlgorithmPropertyVectorEdit(valuesAsString, key, values.front().type(), algProps);
					connect(vectorEdit, SIGNAL(editingFinished()), this, SIGNAL(requestViewUpdate()));
					layout->addRow(key, vectorEdit);
				}
			}
		}
		
		//
		///* Integrating colour tables into property-window (test!) */
		//VtkStationSource* test = dynamic_cast<VtkStationSource*>(_item->algorithm());
		//if (test)
		//{
		//	std::map<std::string, GEOLIB::Color> colors = test->getColorLookupTable();
		//	if (!colors.empty())
		//	{
		//		ColorTableModel* ctm = new ColorTableModel(colors);
		//		ColorTableView* ctv = new ColorTableView();
		//		ctv->setModel(ctm);
		//		ctv->setItemDelegate(new ColorTableViewDelegate);
		//		vbox->addWidget(ctv);
		//		ctv->resizeRowsToContents();
		//	}
		//}

		/**/

		emit requestViewUpdate();
	}
	else
	{
		actorPropertiesGroupBox->setEnabled(false);
		this->activeScalarComboBox->clear();
	}

}

void VtkVisTabWidget::on_diffuseColorPickerButton_colorPicked( QColor color )
{
	static_cast<vtkActor*>(_item->actor())->GetProperty()->SetDiffuseColor(
		color.redF(), color.greenF(), color.blueF());

	emit requestViewUpdate();
}

void VtkVisTabWidget::on_visibleEdgesCheckBox_stateChanged( int state )
{
	if (state == Qt::Checked)
	{
		static_cast<vtkActor*>(_item->actor())->GetProperty()->SetEdgeVisibility(1);
		edgeColorPickerButton->setEnabled(true);
	}
	else
	{
		static_cast<vtkActor*>(_item->actor())->GetProperty()->SetEdgeVisibility(0);
		edgeColorPickerButton->setEnabled(false);
	}

	emit requestViewUpdate();
}

void VtkVisTabWidget::on_edgeColorPickerButton_colorPicked( QColor color )
{
	static_cast<vtkActor*>(_item->actor())->GetProperty()->SetEdgeColor(
		color.redF(), color.greenF(), color.blueF());
	emit requestViewUpdate();
}

void VtkVisTabWidget::on_opacitySlider_sliderMoved( int value )
{
	static_cast<vtkActor*>(_item->actor())->GetProperty()->SetOpacity(value / 100.0);
	emit requestViewUpdate();
}

void VtkVisTabWidget::on_scaleZ_textChanged(const QString &text)
{
	bool ok=true;
	double scale = text.toDouble(&ok);

	if (ok)
	{
		vtkTransform* transform = 
			static_cast<vtkTransform*>(_item->transformFilter()->GetTransform());
		transform->Identity();
		transform->Scale(1.0, 1.0, scale);
		_item->transformFilter()->Modified();

		for (int i = 0; i < _item->childCount(); i++)
		{
			VtkVisPipelineItem* childItem = _item->child(i);
			if (childItem)
			{
				VtkCompositeColorByHeightFilter* colorFilter =
					dynamic_cast<VtkCompositeColorByHeightFilter*>
					(childItem->compositeFilter());
				if (colorFilter)
					VtkColorByHeightFilter::SafeDownCast(
					colorFilter->GetOutputAlgorithm())->SetTableRangeScaling(scale);
			}
		}		

		emit requestViewUpdate();
	}
}

void VtkVisTabWidget::addColorTable()
{

}

void VtkVisTabWidget::SetActiveAttributeOnItem( int idx )
{
	_item->SetActiveAttribute(idx, 0);
	emit requestViewUpdate();
}
