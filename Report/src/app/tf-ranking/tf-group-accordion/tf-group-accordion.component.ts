import {Component, Input, OnInit} from '@angular/core';
import {TranscriptionFactorGroup} from "../../types/types";

@Component({
  selector: 'tf-group-accordion',
  templateUrl: './tf-group-accordion.component.html',
  styleUrls: ['./tf-group-accordion.component.css']
})
export class TfGroupAccordionComponent implements OnInit {
  @Input()
    // @ts-ignore
  tfGroup: TranscriptionFactorGroup;

  @Input()
  index: number | undefined;

  visible: boolean = false
  title: string = "Undefined";

  constructor() {
  }

  toggleVisibility() {
    this.visible = !this.visible;
  }

  ngOnInit(): void {
    if (this.index) {
      this.title = this.index + ". " + this.tfGroup.name;
    } else {
      this.title = "Basic data";
    }
  }
}
