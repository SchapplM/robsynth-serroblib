% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR12_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR12_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:03
% EndTime: 2019-02-26 20:55:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (48->21), mult. (142->35), div. (0->0), fcn. (98->4), ass. (0->17)
t138 = sin(qJ(3));
t140 = cos(qJ(3));
t150 = r_i_i_C(3) + qJ(4);
t151 = pkin(3) - r_i_i_C(2);
t155 = -(t150 * t138 + t151 * t140) * qJD(3) + t140 * qJD(4);
t154 = t151 * qJD(3) - qJD(4);
t152 = t151 * t138 - t150 * t140 + qJ(2);
t139 = sin(qJ(1));
t149 = qJD(1) * t139;
t141 = cos(qJ(1));
t148 = qJD(1) * t141;
t147 = qJD(3) * t139;
t146 = qJD(3) * t141;
t144 = qJD(1) * t151;
t143 = qJD(1) * t150;
t142 = qJD(2) + (-pkin(1) - pkin(7) - r_i_i_C(1)) * qJD(1) - t155;
t1 = [t142 * t141 - t152 * t149, t148 (t141 * t144 + t150 * t147) * t140 + (-t154 * t139 + t141 * t143) * t138, t138 * t147 - t140 * t148, 0, 0; t142 * t139 + t152 * t148, t149 (t139 * t144 - t150 * t146) * t140 + (t139 * t143 + t154 * t141) * t138, -t138 * t146 - t140 * t149, 0, 0; 0, 0, t155, qJD(3) * t140, 0, 0;];
JaD_transl  = t1;
