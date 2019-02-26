% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRRR1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRRR1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:18:47
% EndTime: 2019-02-26 20:18:47
% DurationCPUTime: 0.16s
% Computational Cost: add. (180->39), mult. (333->75), div. (0->0), fcn. (304->10), ass. (0->38)
t224 = sin(pkin(12));
t226 = cos(pkin(12));
t229 = sin(qJ(2));
t227 = cos(pkin(6));
t231 = cos(qJ(2));
t247 = t227 * t231;
t256 = -t224 * t229 + t226 * t247;
t255 = r_i_i_C(3) + pkin(9) + pkin(8);
t223 = qJ(3) + qJ(4);
t220 = sin(t223);
t222 = qJD(3) + qJD(4);
t254 = t220 * t222;
t221 = cos(t223);
t253 = t221 * t222;
t225 = sin(pkin(6));
t252 = t222 * t225;
t228 = sin(qJ(3));
t250 = t225 * t228;
t249 = t225 * t229;
t248 = t227 * t229;
t215 = t224 * t231 + t226 * t248;
t210 = t256 * qJD(2);
t240 = t226 * t252 - t210;
t246 = (-t215 * t253 + t240 * t220) * r_i_i_C(1) + (t215 * t254 + t240 * t221) * r_i_i_C(2);
t236 = t224 * t248 - t226 * t231;
t237 = t224 * t247 + t226 * t229;
t212 = t237 * qJD(2);
t239 = -t224 * t252 + t212;
t245 = (t239 * t220 + t236 * t253) * r_i_i_C(1) + (t239 * t221 - t236 * t254) * r_i_i_C(2);
t241 = qJD(2) * t225 * t231;
t235 = -t222 * t227 - t241;
t243 = t222 * t249;
t244 = (t235 * t220 - t221 * t243) * r_i_i_C(1) + (t220 * t243 + t235 * t221) * r_i_i_C(2);
t230 = cos(qJ(3));
t238 = pkin(3) * t230 + r_i_i_C(1) * t221 - r_i_i_C(2) * t220 + pkin(2);
t234 = qJD(2) * t238;
t233 = -pkin(3) * qJD(3) * t228 + (-r_i_i_C(1) * t220 - r_i_i_C(2) * t221) * t222;
t1 = [0, -t255 * t212 - t233 * t237 + t236 * t234 (t212 * t228 + (-t224 * t250 + t230 * t236) * qJD(3)) * pkin(3) + t245, t245, 0, 0; 0, t255 * t210 - t215 * t234 + t233 * t256 (-t210 * t228 + (-t215 * t230 + t226 * t250) * qJD(3)) * pkin(3) + t246, t246, 0, 0; 0 (t233 * t231 + (-t238 * t229 + t255 * t231) * qJD(2)) * t225 (-t228 * t241 + (-t227 * t228 - t230 * t249) * qJD(3)) * pkin(3) + t244, t244, 0, 0;];
JaD_transl  = t1;
