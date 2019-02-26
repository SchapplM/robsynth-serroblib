% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:04:32
% EndTime: 2019-02-26 21:04:32
% DurationCPUTime: 0.36s
% Computational Cost: add. (417->71), mult. (438->105), div. (0->0), fcn. (327->10), ass. (0->58)
t268 = qJ(3) + qJ(4);
t262 = pkin(10) + t268;
t261 = cos(t262);
t272 = cos(qJ(6));
t324 = r_i_i_C(1) * t272 + pkin(5);
t325 = t261 * t324;
t316 = pkin(9) + r_i_i_C(3);
t323 = qJD(1) * t325;
t270 = sin(qJ(3));
t298 = t316 * t261;
t260 = sin(t262);
t312 = pkin(5) * t260;
t264 = sin(t268);
t315 = pkin(4) * t264;
t321 = -pkin(3) * t270 - qJ(2) + t298 - t312 - t315;
t265 = cos(t268);
t269 = sin(qJ(6));
t318 = -r_i_i_C(2) * t261 * t269 + pkin(4) * t265;
t274 = cos(qJ(1));
t285 = qJD(1) * t260 + qJD(6);
t320 = t285 * t274;
t319 = t316 * t260;
t271 = sin(qJ(1));
t267 = qJD(3) + qJD(4);
t306 = t267 * t274;
t317 = -t261 * t306 + t285 * t271;
t313 = pkin(4) * t267;
t310 = -pkin(1) - qJ(5) - pkin(8) - pkin(7);
t309 = pkin(3) * qJD(3);
t308 = t260 * t269;
t307 = t267 * t272;
t304 = qJD(1) * t271;
t263 = qJD(1) * t274;
t303 = qJD(6) * t269;
t302 = qJD(6) * t272;
t273 = cos(qJ(3));
t299 = t273 * t309;
t297 = t267 * t308;
t295 = t261 * t267 * t271;
t290 = t261 * t303;
t289 = t261 * t302;
t287 = r_i_i_C(2) * t289;
t286 = -qJD(6) * t260 - qJD(1);
t283 = t286 * t274;
t282 = qJD(1) * (pkin(3) * t273 + t318);
t281 = -r_i_i_C(2) * t308 - t298;
t280 = t271 * r_i_i_C(2) * t297 + t323 * t274 + t316 * (t260 * t263 + t295);
t279 = qJD(1) * t318;
t278 = t260 * t307 + t290;
t277 = (r_i_i_C(1) * t290 + t287) * t274 + t323 * t271 + (t316 * t304 + t324 * t306) * t260;
t276 = qJD(2) + t265 * t313 + t299 + (pkin(5) * t261 + t319) * t267;
t275 = (-t318 - t319 - t325) * t267 + (r_i_i_C(1) * t303 + r_i_i_C(2) * t302) * t260;
t240 = -t264 * t313 - t270 * t309;
t237 = t272 * t320 + (t261 * t307 + t286 * t269) * t271;
t236 = t286 * t272 * t271 + (-t295 - t320) * t269;
t235 = t269 * t283 - t317 * t272;
t234 = t317 * t269 + t272 * t283;
t1 = [t235 * r_i_i_C(1) + t234 * r_i_i_C(2) - t271 * qJD(5) + t276 * t274 + (t321 * t271 + t310 * t274) * qJD(1), t263, t274 * t282 + (-t278 * r_i_i_C(1) - t267 * t312 + t240 - t287) * t271 + t280, t274 * t279 + ((-r_i_i_C(1) * t269 - r_i_i_C(2) * t272) * t261 * qJD(6) + (-t260 * t324 - t315) * t267) * t271 + t280, -t304, r_i_i_C(1) * t236 - r_i_i_C(2) * t237; t237 * r_i_i_C(1) + t236 * r_i_i_C(2) + qJD(5) * t274 + t276 * t271 + (t310 * t271 - t321 * t274) * qJD(1), t304, t271 * t282 + (t281 * t267 - t240) * t274 + t277, t271 * t279 + (t281 + t315) * t306 + t277, t263, -r_i_i_C(1) * t234 + r_i_i_C(2) * t235; 0, 0, t275 - t299, t275, 0, t278 * r_i_i_C(2) + (-t289 + t297) * r_i_i_C(1);];
JaD_transl  = t1;
