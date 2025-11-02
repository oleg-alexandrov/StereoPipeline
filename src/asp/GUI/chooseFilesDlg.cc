// __BEGIN_LICENSE__
//  Copyright (c) 2006-2013, United States Government as represented by the
//  Administrator of the National Aeronautics and Space Administration. All
//  rights reserved.
//
//  The NGT platform is licensed under the Apache License, Version 2.0 (the
//  "License"); you may not use this file except in compliance with the
//  License. You may obtain a copy of the License at
//  http://www.apache.org/licenses/LICENSE-2.0
//
//  Unless required by applicable law or agreed to in writing, software
//  distributed under the License is distributed on an "AS IS" BASIS,
//  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//  See the License for the specific language governing permissions and
//  limitations under the License.
// __END_LICENSE__

#include <asp/GUI/chooseFilesDlg.h>
#include <asp/GUI/GuiBase.h>
#include <asp/Core/StereoSettings.h>

#include <QWidget>
#include <QHeaderView>
#include <QVBoxLayout>
#include <QScrollBar>

namespace asp {

using namespace vw;

// Allow the user to choose which files to hide/show in the GUI.
// User's choice will be processed by MainWidget::showFilesChosenByUser().
chooseFilesDlg::chooseFilesDlg(QWidget * parent):
  QWidget(parent){

  setWindowModality(Qt::ApplicationModal);
  
  int spacing = 0;
  
  QVBoxLayout * vBoxLayout = new QVBoxLayout(this);
  vBoxLayout->setSpacing(spacing);
  vBoxLayout->setAlignment(Qt::AlignLeft);
  
  // The layout having the file names. It will be filled in
  // dynamically later.
  m_filesTable = new QTableWidget();
  
  //m_filesTable->horizontalHeader()->hide();
  m_filesTable->verticalHeader()->hide();
    
  vBoxLayout->addWidget(m_filesTable);
  
  return;
}
  
chooseFilesDlg::~chooseFilesDlg(){}

void chooseFilesDlg::chooseFiles(const std::vector<std::string> & image_files) {

  // See the top of this file for documentation.

  int numFiles = image_files.size();
  int numCols = 2;
  m_filesTable->setRowCount(numFiles);
  m_filesTable->setColumnCount(numCols);

  for (int fileIter = 0; fileIter < numFiles; fileIter++){

    // Checkbox
    QTableWidgetItem *item = new QTableWidgetItem(1);
    item->data(Qt::CheckStateRole);
    if (!asp::stereo_settings().hide_all)
      item->setCheckState(Qt::Checked);
    else
      item->setCheckState(Qt::Unchecked);
      
    m_filesTable->setItem(fileIter, 0, item);

    // Set the filename in the table
    std::string fileName = image_files[fileIter];
    item = new QTableWidgetItem(fileName.c_str());
    item->setFlags(Qt::NoItemFlags);
    item->setForeground(QColor::fromRgb(0, 0, 0));
    m_filesTable->setItem(fileIter, numCols - 1, item);

    // To be able to quickly look up an image
    image_to_row[fileName] = fileIter;
  }

  QStringList rowNamesList;
  for (int fileIter = 0; fileIter < numFiles; fileIter++) rowNamesList << "";
  m_filesTable->setVerticalHeaderLabels(rowNamesList);

  QStringList colNamesList;
  for (int colIter = 0; colIter < numCols; colIter++) colNamesList << "";
  m_filesTable->setHorizontalHeaderLabels(colNamesList);
  QTableWidgetItem * hs = m_filesTable->horizontalHeaderItem(0);
  hs->setBackground(QBrush(QColor("lightgray")));

  m_filesTable->setSelectionMode(QTableWidget::ExtendedSelection);
  std::string style = std::string("QTableWidget::indicator:unchecked ")
    + "{background-color:white; border: 1px solid black;}; " +
    "selection-background-color: rgba(128, 128, 128, 40);";

  m_filesTable->setSelectionMode(QTableWidget::NoSelection);
  m_filesTable->setStyleSheet(style.c_str());
  
  // Horizontal header caption
  QTableWidgetItem *item = new QTableWidgetItem("Hide/show all");
  item->setFlags(Qt::NoItemFlags);
  item->setForeground(QColor::fromRgb(0, 0, 0));
  m_filesTable->setHorizontalHeaderItem(1, item);
  
  m_filesTable->resizeColumnsToContents();
  m_filesTable->resizeRowsToContents();

  // The processing of user's choice happens in MainWidget::showFilesChosenByUser()

  return;
}

// Quickly find in what table row the current image is  
int chooseFilesDlg::imageRow(std::string const& image) const {
  auto it = image_to_row.find(image);
  if (it == image_to_row.end()) {
    popUp("Cannot find image in table.");
    return 0;
  }
  return it->second;
}
  
// Check if the given image is hidden (not shown) based on the table checkbox  
bool chooseFilesDlg::isHidden(std::string const& image) const {

  int row = imageRow(image);
  QTableWidgetItem *item = m_filesTable->item(row, 0);
  // TODO(oalexan1): Use below a function called image(int id).
  // There are more places like that.
  std::string curr_image = (m_filesTable->item(row, 1)->data(0)).toString().toStdString();
  if (image == curr_image)
    return (item->checkState() == Qt::Unchecked);
  return false;
}

// Hide the given image  
void chooseFilesDlg::hide(std::string const& image) {
  int image_id = imageRow(image);
  chooseFilesDlg::hide(image_id);
}
void chooseFilesDlg::hide(int image_id) {
  QTableWidgetItem *item = m_filesTable->item(image_id, 0);
  item->setCheckState(Qt::Unchecked);
}

// Show the given image by turning on the checkbox in the table
void chooseFilesDlg::unhide(std::string const& image) {
  int image_id = imageRow(image);
  chooseFilesDlg::unhide(image_id);
}
void chooseFilesDlg::unhide(int image_id) {
  QTableWidgetItem *item = m_filesTable->item(image_id, 0);
  item->setCheckState(Qt::Checked);
}

// Show this many of the first several input images
void chooseFilesDlg::setNumImagesToShow(int num) {

  int rows = m_filesTable->rowCount();
  for (int row = 0; row < std::min(num, rows); row++) {
    QTableWidgetItem *item = m_filesTable->item(row, 0);
    item->setCheckState(Qt::Checked);
  }
  for (int row = num; row < rows; row++) {
    QTableWidgetItem *item = m_filesTable->item(row, 0);
    item->setCheckState(Qt::Unchecked);
  }
}
  
// Show all images
void chooseFilesDlg::showAllImages() {
  int rows = m_filesTable->rowCount();
  for (int row = 0; row < rows; row++) {
    QTableWidgetItem *item = m_filesTable->item(row, 0);
    item->setCheckState(Qt::Checked);
  }
}

// Number of images being shown
int chooseFilesDlg::numShown() {
  int num = 0;
  int rows = m_filesTable->rowCount();
  for (int row = 0; row < rows; row++) {
    QTableWidgetItem *item = m_filesTable->item(row, 0);
    num += (item->checkState() == Qt::Checked);
  }

  return num;
}

// If some images are shown, hide all. Else, show all.
void chooseFilesDlg::hideShowAll() {

  int rows = m_filesTable->rowCount();

  // See if all files are hidden
  bool allOff = true;
  for (int rowIter = 0; rowIter < rows; rowIter++){
    QTableWidgetItem *item = m_filesTable->item(rowIter, 0);
    if (item->checkState() == Qt::Checked){
      allOff = false;
    }
  }
  
  // If all files are hidden, we will show all. Else hide all.
  for (int rowIter = 0; rowIter < rows; rowIter++){
    QTableWidgetItem *item = m_filesTable->item(rowIter, 0);
    // TODO(oalexan1): Use below a function called image(int id).
    std::string fileName = (m_filesTable->item(rowIter, 1)->data(0)).toString().toStdString();
    if (allOff)
      item->setCheckState(Qt::Checked);
    else
      item->setCheckState(Qt::Unchecked);
  }

  // Force the horizontal scrollbar in the table to go left, so one can see
  // the checkboxes.
  QScrollBar * hScrollBar = m_filesTable->horizontalScrollBar();
  hScrollBar->triggerAction(QScrollBar::SliderToMinimum);
}

void chooseFilesDlg::viewOtherImage(int delta) {
  if (delta != -1 && delta != 1) 
    return;
  
  int rows = m_filesTable->rowCount();
  
  if (rows == 0) 
    return;
  
  // First see how many images have a checkbox now, so are being shown
  std::set<int> shown;
  for (int rowIter = 0; rowIter < rows; rowIter++) {
    QTableWidgetItem *item = m_filesTable->item(rowIter, 0);
    if (item->checkState() == Qt::Checked)
      shown.insert(rowIter);
  }
  
  // If no images are being shown or more than one, show the first
  int shownRow = 0;
  if (shown.size() == 1) {
    // Else show the next or previous image. Note how we add 'rows'
    // before we find the remainder, as delta could be negative.
    shownRow = *shown.begin();
    shownRow = (shownRow + delta + rows) % rows;
  }
  
  // Show the next/previous one, and hide the rest 
  for (int rowIter = 0; rowIter < rows; rowIter++){
    QTableWidgetItem *item = m_filesTable->item(rowIter, 0);
    if (rowIter == shownRow)
      item->setCheckState(Qt::Checked);
    else
      item->setCheckState(Qt::Unchecked);
  }
  
  // Print count and image file (count starts from 1)
  // TODO(oalexan1): Implement a function called image(int id) to avoid this
  // lengthy text. It can be used in other places too.
  std::string fileName = (m_filesTable->item(shownRow, 1)->data(0)).toString().toStdString();
  vw_out() << "Image: " << shownRow + 1  << ' ' << fileName << "\n";
}
  
void chooseFilesDlg::keyPressEvent(QKeyEvent *event) {
  // std::cout << "Key was pressed " << event->key() << std::endl;
}

} // namespace asp
