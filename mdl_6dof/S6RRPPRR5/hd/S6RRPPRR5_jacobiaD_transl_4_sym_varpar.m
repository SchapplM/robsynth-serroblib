% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR5
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
% Datum: 2019-02-26 21:30
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:30:48
% EndTime: 2019-02-26 21:30:48
% DurationCPUTime: 0.12s
% Computational Cost: add. (107->34), mult. (305->51), div. (0->0), fcn. (276->6), ass. (0->26)
t147 = sin(pkin(6));
t168 = t147 * (pkin(8) - r_i_i_C(3) - qJ(4));
t167 = r_i_i_C(2) + qJ(3);
t149 = sin(qJ(2));
t150 = sin(qJ(1));
t166 = t150 * t149;
t151 = cos(qJ(2));
t165 = t150 * t151;
t152 = cos(qJ(1));
t164 = t152 * t149;
t163 = t152 * t151;
t162 = qJD(1) * t147;
t161 = qJD(2) * t149;
t160 = t147 * qJD(4);
t159 = pkin(2) + pkin(3) + r_i_i_C(1);
t148 = cos(pkin(6));
t157 = t148 * t166;
t156 = t148 * t163;
t155 = qJD(2) * t148 + qJD(1);
t154 = t148 * t165 + t164;
t153 = t148 * t164 + t165;
t142 = -qJD(1) * t157 - t150 * t161 + t155 * t163;
t141 = t154 * qJD(1) + t153 * qJD(2);
t140 = t153 * qJD(1) + t154 * qJD(2);
t139 = -qJD(1) * t156 - qJD(2) * t163 + t155 * t166;
t1 = [-t152 * t160 - (-t156 + t166) * qJD(3) - t167 * t141 - t159 * t142 + (-t152 * pkin(1) - t150 * t168) * qJD(1) -(t157 - t163) * qJD(3) - t167 * t140 + t159 * t139, -t139, -t152 * t162, 0, 0; -t150 * t160 + t154 * qJD(3) - t167 * t139 - t159 * t140 + (-t150 * pkin(1) + t152 * t168) * qJD(1), t153 * qJD(3) - t159 * t141 + t167 * t142, t141, -t150 * t162, 0, 0; 0 (t149 * qJD(3) + (-t159 * t149 + t167 * t151) * qJD(2)) * t147, t147 * t161, 0, 0, 0;];
JaD_transl  = t1;
