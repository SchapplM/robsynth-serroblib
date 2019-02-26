% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:16
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:16:01
% EndTime: 2019-02-26 21:16:01
% DurationCPUTime: 0.25s
% Computational Cost: add. (381->55), mult. (424->84), div. (0->0), fcn. (334->10), ass. (0->51)
t241 = sin(qJ(3));
t242 = cos(qJ(4));
t232 = pkin(4) * t242 + pkin(3);
t239 = qJ(4) + qJ(5);
t235 = sin(t239);
t274 = r_i_i_C(2) * t235;
t236 = cos(t239);
t275 = r_i_i_C(1) * t236;
t250 = t232 - t274 + t275;
t243 = cos(qJ(3));
t272 = r_i_i_C(3) + pkin(9) + pkin(8);
t279 = t272 * t243;
t246 = -t250 * t241 + t279;
t285 = qJD(1) * t246;
t240 = sin(qJ(4));
t271 = pkin(4) * qJD(4);
t263 = t240 * t271;
t284 = (-t232 * t241 + t279) * qJD(3) - t243 * t263;
t237 = qJD(4) + qJD(5);
t266 = qJD(1) * t243;
t255 = -t237 + t266;
t282 = t236 * t255;
t281 = t240 * (-qJD(4) + t266);
t256 = t237 * t243 - qJD(1);
t265 = qJD(3) * t241;
t278 = -t235 * t265 + t256 * t236;
t273 = r_i_i_C(2) * t236;
t252 = r_i_i_C(1) * t235 + t273;
t277 = t252 * t237 + t263;
t276 = pkin(4) * t240;
t269 = t237 * t241;
t238 = qJ(1) + pkin(11);
t233 = sin(t238);
t234 = cos(t238);
t251 = t255 * t235;
t227 = t233 * t251 - t278 * t234;
t248 = t256 * t235 + t236 * t265;
t228 = t233 * t282 + t248 * t234;
t268 = t227 * r_i_i_C(1) + t228 * r_i_i_C(2);
t229 = t278 * t233 + t234 * t251;
t230 = t248 * t233 - t234 * t282;
t267 = -t229 * r_i_i_C(1) + t230 * r_i_i_C(2);
t264 = qJD(3) * t243;
t262 = t242 * t271;
t261 = pkin(7) + t276;
t259 = t272 * t241;
t249 = -t232 * t243 - pkin(2) - t259;
t247 = t240 * t265 + (-qJD(4) * t243 + qJD(1)) * t242;
t245 = t277 * t241 + (-t250 * t243 - t259) * qJD(3);
t231 = t269 * t274;
t1 = [t234 * t262 + t230 * r_i_i_C(1) + t229 * r_i_i_C(2) - t284 * t233 + (-cos(qJ(1)) * pkin(1) - t261 * t233 + t249 * t234) * qJD(1), 0, -t233 * t285 + t245 * t234 (t233 * t281 + t247 * t234) * pkin(4) + t268, t268, 0; t233 * t262 - t228 * r_i_i_C(1) + t227 * r_i_i_C(2) + t284 * t234 + (-sin(qJ(1)) * pkin(1) + t261 * t234 + t249 * t233) * qJD(1), 0, t245 * t233 + t234 * t285 (t247 * t233 - t234 * t281) * pkin(4) + t267, t267, 0; 0, 0, t246 * qJD(3) - t277 * t243, t231 + (-t237 * t275 - t262) * t241 + (-t252 - t276) * t264, -t264 * t273 + t231 + (-t235 * t264 - t236 * t269) * r_i_i_C(1), 0;];
JaD_transl  = t1;
