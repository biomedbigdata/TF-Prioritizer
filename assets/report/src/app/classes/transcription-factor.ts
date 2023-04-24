import {dcg, ranking} from "../../interfaces";

export class TranscriptionFactor {
  private dcgs: dcg = {};

  constructor(private _name: string, ranking: ranking) {
    for (let assay of Object.keys(ranking).filter(val => val != 'default')) {
      for (let pairing of Object.keys(ranking[assay])) {
        let ranked = ranking[assay][pairing];

        let rank = ranked.indexOf(_name);
        if (rank == -1) {
          continue;
        }

        let dcgValue = 1 - (rank / ranked.length);

        if (!(assay in this.dcgs)) {
          this.dcgs[assay] = {};
        }

        this.dcgs[assay][pairing] = dcgValue;
      }
    }
  }


  get name(): string {
    return this._name;
  }

  getDcg(assay: string, pairing: string) {
    if (!(assay in this.dcgs)) {
      return 0;
    }

    if (!(pairing in this.dcgs[assay])) {
      return 0;
    }

    return this.dcgs[assay][pairing];
  }

  getDcgSum(assays: Set<string>, pairings: Set<string>) {
    let sum = 0;

    for (let assay of assays) {
      for (let pairing of pairings) {
        sum += this.getDcg(assay, pairing);
      }
    }

    return sum;
  }
}
