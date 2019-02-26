% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRR3
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
% Datum: 2019-02-26 21:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:55:11
% EndTime: 2019-02-26 21:55:11
% DurationCPUTime: 0.35s
% Computational Cost: add. (715->73), mult. (598->96), div. (0->0), fcn. (471->12), ass. (0->64)
t267 = qJ(4) + qJ(5);
t260 = sin(t267);
t269 = sin(qJ(4));
t306 = pkin(4) * qJD(4);
t264 = qJD(4) + qJD(5);
t312 = pkin(5) * t264;
t248 = -t260 * t312 - t269 * t306;
t261 = cos(t267);
t272 = cos(qJ(4));
t252 = t272 * pkin(4) + pkin(5) * t261;
t250 = pkin(3) + t252;
t313 = pkin(5) * t260;
t251 = t269 * pkin(4) + t313;
t265 = qJ(2) + pkin(11);
t257 = sin(t265);
t258 = cos(t265);
t307 = r_i_i_C(3) + pkin(10) + pkin(9) + pkin(8);
t281 = t307 * t258 - sin(qJ(2)) * pkin(2);
t322 = (-t250 * t257 + t281) * qJD(2) + (t251 + qJ(3) + pkin(7)) * qJD(1) + t258 * t248;
t262 = qJ(6) + t267;
t254 = sin(t262);
t310 = r_i_i_C(2) * t254;
t255 = cos(t262);
t311 = r_i_i_C(1) * t255;
t282 = t250 - t310 + t311;
t276 = -t282 * t257 + t281;
t274 = cos(qJ(1));
t259 = qJD(6) + t264;
t289 = t258 * t259 - qJD(1);
t320 = t274 * t289;
t318 = -t307 * t257 - cos(qJ(2)) * pkin(2);
t309 = r_i_i_C(2) * t255;
t286 = r_i_i_C(1) * t254 + t309;
t317 = t286 * t259 - t248;
t299 = qJD(1) * t258;
t288 = -t259 + t299;
t271 = sin(qJ(1));
t296 = qJD(2) * t257;
t292 = t271 * t296;
t316 = t288 * t274 - t292;
t304 = t257 * t259;
t291 = t274 * t296;
t278 = t288 * t271 + t291;
t243 = t278 * t254 - t255 * t320;
t244 = t254 * t320 + t278 * t255;
t302 = t243 * r_i_i_C(1) + t244 * r_i_i_C(2);
t283 = t289 * t271;
t245 = t316 * t254 + t255 * t283;
t246 = t254 * t283 - t316 * t255;
t301 = -t245 * r_i_i_C(1) + t246 * r_i_i_C(2);
t298 = qJD(1) * t271;
t297 = qJD(1) * t274;
t295 = qJD(2) * t258;
t294 = t261 * t312;
t293 = t259 * t311;
t287 = -t264 + t299;
t285 = t251 * t299 + t248;
t284 = t261 * (-t258 * t264 + qJD(1));
t249 = t272 * t306 + t294;
t279 = qJD(1) * t252 - t249 * t258 + t251 * t296;
t277 = qJD(3) + t249 + (-t250 * t258 - pkin(1) + t318) * qJD(1);
t275 = t317 * t257 + (-t282 * t258 + t318) * qJD(2);
t247 = t304 * t310;
t1 = [t246 * r_i_i_C(1) + t245 * r_i_i_C(2) - t322 * t271 + t277 * t274, t275 * t274 - t276 * t298, t297, t285 * t271 + t279 * t274 + t302 (t274 * t284 + (t287 * t271 + t291) * t260) * pkin(5) + t302, t302; -t244 * r_i_i_C(1) + t243 * r_i_i_C(2) + t277 * t271 + t322 * t274, t275 * t271 + t276 * t297, t298, t279 * t271 - t285 * t274 + t301 (t271 * t284 + (-t287 * t274 + t292) * t260) * pkin(5) + t301, t301; 0, t276 * qJD(2) - t317 * t258, 0, t247 + (-t249 - t293) * t257 + (-t251 - t286) * t295, t247 + (-t293 - t294) * t257 + (-t286 - t313) * t295, -t295 * t309 + t247 + (-t254 * t295 - t255 * t304) * r_i_i_C(1);];
JaD_transl  = t1;
