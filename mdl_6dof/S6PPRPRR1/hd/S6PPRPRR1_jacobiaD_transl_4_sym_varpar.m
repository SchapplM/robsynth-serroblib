% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PPRPRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PPRPRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:39:48
% EndTime: 2019-02-26 19:39:48
% DurationCPUTime: 0.10s
% Computational Cost: add. (54->30), mult. (201->65), div. (0->0), fcn. (198->12), ass. (0->28)
t150 = sin(pkin(13));
t155 = cos(pkin(13));
t160 = sin(qJ(3));
t161 = cos(qJ(3));
t169 = (t150 * t161 + t155 * t160) * qJD(3);
t168 = pkin(3) * qJD(3);
t152 = sin(pkin(11));
t154 = sin(pkin(6));
t167 = t152 * t154;
t159 = cos(pkin(6));
t166 = t152 * t159;
t157 = cos(pkin(11));
t165 = t154 * t157;
t164 = t157 * t159;
t148 = (t150 * t160 - t155 * t161) * qJD(3);
t158 = cos(pkin(7));
t156 = cos(pkin(12));
t153 = sin(pkin(7));
t151 = sin(pkin(12));
t147 = -t151 * t166 + t157 * t156;
t146 = -t157 * t151 - t156 * t166;
t145 = t151 * t164 + t152 * t156;
t144 = -t152 * t151 + t156 * t164;
t143 = t158 * t169;
t142 = t158 * t148;
t141 = t153 * t169;
t140 = t153 * t148;
t1 = [0, 0 (-t141 * t167 - t146 * t143 + t147 * t148) * r_i_i_C(1) + (t140 * t167 + t146 * t142 + t147 * t169) * r_i_i_C(2) + (-t147 * t161 + (-t146 * t158 - t153 * t167) * t160) * t168, 0, 0, 0; 0, 0 (t141 * t165 - t144 * t143 + t145 * t148) * r_i_i_C(1) + (-t140 * t165 + t144 * t142 + t145 * t169) * r_i_i_C(2) + (-t145 * t161 + (-t144 * t158 + t153 * t165) * t160) * t168, 0, 0, 0; 0, 0 (-t153 * t160 * t168 - t141 * r_i_i_C(1) + t140 * r_i_i_C(2)) * t159 + ((-t143 * t156 + t148 * t151) * r_i_i_C(1) + (t142 * t156 + t151 * t169) * r_i_i_C(2) + (-t156 * t158 * t160 - t151 * t161) * t168) * t154, 0, 0, 0;];
JaD_transl  = t1;
