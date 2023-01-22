export interface tf {
  "symbol": string;
  "groupMeanExpression": {
    [group: string]: number;
  }

  ensgs: string[];
  "groupTpm": {
    [group: string]: number;
  }
  "pairingLog2fc": {
    [pairing: string]: number;
  }
}

export interface tfGroup {
  "symbol": string;
  "hmDcg": {
    [hm: string]: number;
  }
  "transcriptionFactors": tf[];
  "biophysicalLogo"?: string;
  "tfSequence": {};
  "heatmaps": {};
  "igv": {};
  distributionPlots: {};
}

export interface dataInterface {
  "groups": tfGroup[];
  coOccurrence: {
    "name": string;
    "series": {
      "name": string;
      "value": number;
    }[];
  }[];
  "configs": {
    [mod: string]: {
      [config: string]: string | number | boolean;
    }
  }
  "regressionCoefficients": {
    [hm: string]: {
      [pairing: string]: {
        [tf: string]: number;
      }
    }
  }
}
