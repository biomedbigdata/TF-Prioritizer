import {ComponentFixture, TestBed} from '@angular/core/testing';

import {TfRankingComponent} from './tf-ranking.component';

describe('TfGroupsContainerComponent', () => {
  let component: TfRankingComponent;
  let fixture: ComponentFixture<TfRankingComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [TfRankingComponent]
    })
      .compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(TfRankingComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});