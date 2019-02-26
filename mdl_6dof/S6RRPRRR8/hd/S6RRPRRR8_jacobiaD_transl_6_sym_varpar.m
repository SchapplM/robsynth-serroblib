% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:58:17
% EndTime: 2019-02-26 21:58:17
% DurationCPUTime: 0.28s
% Computational Cost: add. (785->76), mult. (601->103), div. (0->0), fcn. (481->12), ass. (0->66)
t255 = pkin(11) + qJ(4);
t253 = qJ(5) + t255;
t248 = cos(t253);
t251 = cos(t255);
t242 = pkin(4) * t251 + pkin(5) * t248;
t236 = cos(pkin(11)) * pkin(3) + pkin(2) + t242;
t247 = sin(t253);
t250 = sin(t255);
t293 = pkin(4) * qJD(4);
t256 = qJD(4) + qJD(5);
t299 = pkin(5) * t256;
t237 = -t247 * t299 - t250 * t293;
t300 = pkin(5) * t247;
t241 = -pkin(4) * t250 - t300;
t257 = sin(qJ(2));
t259 = cos(qJ(2));
t282 = t257 * qJD(3);
t294 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8) + qJ(3);
t304 = t294 * t259;
t308 = (-t236 * t257 + t304) * qJD(2) + (pkin(7) + sin(pkin(11)) * pkin(3) - t241) * qJD(1) + t259 * t237 + t282;
t249 = qJ(6) + t253;
t244 = sin(t249);
t297 = r_i_i_C(2) * t244;
t245 = cos(t249);
t298 = r_i_i_C(1) * t245;
t266 = t236 - t297 + t298;
t262 = -t266 * t257 + t304;
t260 = cos(qJ(1));
t252 = qJD(6) + t256;
t274 = t252 * t259 - qJD(1);
t306 = t260 * t274;
t296 = r_i_i_C(2) * t245;
t271 = r_i_i_C(1) * t244 + t296;
t303 = t271 * t252 - t237;
t287 = qJD(1) * t259;
t273 = -t252 + t287;
t258 = sin(qJ(1));
t285 = qJD(2) * t257;
t279 = t258 * t285;
t302 = t273 * t260 - t279;
t291 = t252 * t257;
t283 = qJD(2) * t260;
t278 = t257 * t283;
t264 = t273 * t258 + t278;
t232 = t264 * t244 - t245 * t306;
t233 = t244 * t306 + t264 * t245;
t290 = t232 * r_i_i_C(1) + t233 * r_i_i_C(2);
t268 = t274 * t258;
t234 = t302 * t244 + t245 * t268;
t235 = t244 * t268 - t302 * t245;
t289 = -t234 * r_i_i_C(1) + t235 * r_i_i_C(2);
t288 = qJD(1) * t258;
t286 = qJD(1) * t260;
t284 = qJD(2) * t259;
t281 = t248 * t299;
t280 = t252 * t298;
t277 = t294 * t257;
t272 = -t256 + t287;
t270 = t241 * t287 - t237;
t269 = t248 * (-t256 * t259 + qJD(1));
t238 = t251 * t293 + t281;
t265 = qJD(1) * t242 - t238 * t259 - t241 * t285;
t263 = t238 + (-t236 * t259 - pkin(1) - t277) * qJD(1);
t261 = qJD(3) * t259 + t303 * t257 + (-t266 * t259 - t277) * qJD(2);
t240 = t291 * t297;
t1 = [t235 * r_i_i_C(1) + t234 * r_i_i_C(2) - t308 * t258 + t263 * t260, t261 * t260 - t262 * t288, -t257 * t288 + t259 * t283, -t270 * t258 + t265 * t260 + t290 (t260 * t269 + (t272 * t258 + t278) * t247) * pkin(5) + t290, t290; -t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t263 * t258 + t308 * t260, t261 * t258 + t262 * t286, t257 * t286 + t258 * t284, t265 * t258 + t270 * t260 + t289 (t258 * t269 + (-t272 * t260 + t279) * t247) * pkin(5) + t289, t289; 0, t262 * qJD(2) - t303 * t259 + t282, t285, t240 + (-t238 - t280) * t257 + (t241 - t271) * t284, t240 + (-t280 - t281) * t257 + (-t271 - t300) * t284, -t284 * t296 + t240 + (-t244 * t284 - t245 * t291) * r_i_i_C(1);];
JaD_transl  = t1;
