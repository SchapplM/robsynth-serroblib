% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:09:49
% EndTime: 2019-02-26 21:09:49
% DurationCPUTime: 0.39s
% Computational Cost: add. (513->74), mult. (521->105), div. (0->0), fcn. (400->9), ass. (0->57)
t324 = pkin(5) + r_i_i_C(1);
t268 = qJD(3) + qJD(4);
t270 = sin(qJ(5));
t317 = r_i_i_C(2) * t270;
t331 = t268 * t317 + qJD(6);
t267 = pkin(10) + qJ(3);
t265 = qJ(4) + t267;
t260 = sin(t265);
t261 = cos(t265);
t272 = cos(qJ(5));
t305 = qJD(5) * t272;
t306 = qJD(5) * t270;
t330 = t331 * t261 + (r_i_i_C(2) * t305 + t306 * t324) * t260;
t271 = sin(qJ(1));
t273 = cos(qJ(1));
t289 = qJD(1) * t261 - qJD(5);
t283 = t289 * t273;
t290 = qJD(5) * t261 - qJD(1);
t285 = t290 * t272;
t310 = t268 * t271;
t300 = t260 * t310;
t237 = t271 * t285 + (t283 - t300) * t270;
t262 = t272 * pkin(5) + pkin(4);
t329 = r_i_i_C(1) * t272 + t262;
t263 = sin(t267);
t314 = pkin(3) * qJD(3);
t302 = t263 * t314;
t269 = -qJ(6) - pkin(9);
t315 = r_i_i_C(3) - t269;
t325 = (-pkin(5) * t306 + t315 * t268) * t261 + (pkin(5) * t270 + pkin(7) + pkin(8) + qJ(2)) * qJD(1) - (t262 * t268 - qJD(6)) * t260 - t302;
t319 = pkin(3) * t263;
t316 = r_i_i_C(3) * t261;
t313 = t261 * t268;
t312 = t261 * t269;
t311 = t268 * t260;
t309 = t268 * t273;
t308 = qJD(1) * t271;
t307 = qJD(1) * t273;
t299 = t260 * t309;
t298 = t260 * t308;
t297 = t260 * t307;
t286 = t329 * t268;
t284 = t290 * t270;
t282 = -r_i_i_C(2) * t272 - t270 * t324;
t281 = -t260 * t329 - t312;
t280 = t269 * t300 + t330 * t271 + t297 * t317 + t307 * t316;
t279 = t269 * t299 + t330 * t273 + t329 * t298 + t308 * t312;
t278 = (-r_i_i_C(3) * t260 - t261 * t329) * t268;
t277 = t289 * t271 + t299;
t264 = cos(t267);
t276 = -t264 * t314 + t278;
t275 = pkin(5) * t305 + qJD(2) + (-t315 * t260 - t261 * t262 - pkin(3) * t264 - cos(pkin(10)) * pkin(2) - pkin(1)) * qJD(1);
t235 = t277 * t270 - t273 * t285;
t274 = r_i_i_C(3) * t313 + (t282 * qJD(5) - t268 * t269) * t261 + (-t286 + t331) * t260;
t238 = -t272 * t283 + (t272 * t311 + t284) * t271;
t236 = t277 * t272 + t273 * t284;
t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) - t325 * t271 + t275 * t273, t307 (-t260 * t317 - t316 + t319) * t308 + t276 * t273 + t279 (-r_i_i_C(3) * t309 - t308 * t317) * t260 + (-r_i_i_C(3) * t308 - t273 * t286) * t261 + t279, t236 * r_i_i_C(2) + t324 * t235, t261 * t309 - t298; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t275 * t271 + t325 * t273, t308, t276 * t271 + (t281 - t319) * t307 + t280, t271 * t278 + t281 * t307 + t280, t238 * r_i_i_C(2) - t324 * t237, t261 * t310 + t297; 0, 0, t274 - t302, t274, t282 * t313 + (-t272 * t324 + t317) * t260 * qJD(5), t311;];
JaD_transl  = t1;
