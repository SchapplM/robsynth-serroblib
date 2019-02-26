% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:31:57
% EndTime: 2019-02-26 21:31:57
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->34), mult. (305->51), div. (0->0), fcn. (276->6), ass. (0->26)
t148 = sin(pkin(6));
t169 = t148 * (pkin(8) - r_i_i_C(3) - qJ(4));
t168 = r_i_i_C(1) + qJ(3);
t150 = sin(qJ(2));
t151 = sin(qJ(1));
t167 = t151 * t150;
t152 = cos(qJ(2));
t166 = t151 * t152;
t153 = cos(qJ(1));
t165 = t153 * t150;
t164 = t153 * t152;
t163 = qJD(1) * t148;
t162 = qJD(2) * t150;
t161 = t148 * qJD(4);
t160 = pkin(2) + pkin(3) - r_i_i_C(2);
t149 = cos(pkin(6));
t158 = t149 * t167;
t157 = t149 * t164;
t156 = qJD(2) * t149 + qJD(1);
t155 = t149 * t166 + t165;
t154 = t149 * t165 + t166;
t143 = -qJD(1) * t158 - t151 * t162 + t156 * t164;
t142 = t155 * qJD(1) + t154 * qJD(2);
t141 = t154 * qJD(1) + t155 * qJD(2);
t140 = -qJD(1) * t157 - qJD(2) * t164 + t156 * t167;
t1 = [-t153 * t161 - (-t157 + t167) * qJD(3) - t168 * t142 - t160 * t143 + (-t153 * pkin(1) - t151 * t169) * qJD(1) -(t158 - t164) * qJD(3) - t168 * t141 + t160 * t140, -t140, -t153 * t163, 0, 0; -t151 * t161 + t155 * qJD(3) - t168 * t140 - t160 * t141 + (-t151 * pkin(1) + t153 * t169) * qJD(1), t154 * qJD(3) - t160 * t142 + t168 * t143, t142, -t151 * t163, 0, 0; 0 (t150 * qJD(3) + (-t160 * t150 + t168 * t152) * qJD(2)) * t148, t148 * t162, 0, 0, 0;];
JaD_transl  = t1;
