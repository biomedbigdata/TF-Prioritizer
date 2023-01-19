import {Component} from '@angular/core';
import {DataManagerService} from "../../../services/data-manager.service";
import {Subject} from "rxjs";

@Component({
  selector: 'app-regression-coefficients',
  templateUrl: './regression-coefficients.component.html',
  styleUrls: ['./regression-coefficients.component.scss']
})
export class RegressionCoefficientsComponent {
  allData = this.dataService.regressionCoefficients;
  $selectedData: Subject<{ [tf: string]: number }> = new Subject<{ [tf: string]: number }>();
  selectedData: { 'name': string, 'value': number }[] = [];
  plotData: Subject<{ 'name': string, 'value': number }[]> = new Subject<{ 'name': string, 'value': number }[]>();

  threshold: number = 0.05;
  $threshold: Subject<number> = new Subject<number>();

  constructor(private dataService: DataManagerService) {
    this.$selectedData.subscribe(data => {
      this.selectedData = dataService.formatPlotData(data).sort((a, b) => b.value - a.value);
      this.updatePlotData();
    });

    this.$threshold.subscribe(threshold => {
      this.threshold = threshold;
      console.log(threshold);
      this.updatePlotData();
    });
  }

  private updatePlotData() {
    this.plotData.next(this.selectedData.filter(d => d.value > this.threshold || d.value < -this.threshold));
  }
}
