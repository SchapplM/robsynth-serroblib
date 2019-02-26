% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_jacobiaD_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:39:33
% EndTime: 2019-02-26 21:39:35
% DurationCPUTime: 1.92s
% Computational Cost: add. (12327->154), mult. (24081->307), div. (726->12), fcn. (31139->15), ass. (0->132)
t327 = sin(pkin(11));
t328 = cos(pkin(11));
t331 = sin(qJ(2));
t333 = cos(qJ(2));
t314 = t327 * t331 - t333 * t328;
t329 = cos(pkin(6));
t311 = t314 * t329;
t306 = qJD(2) * t311;
t355 = t327 * t333 + t328 * t331;
t313 = t355 * qJD(2);
t334 = cos(qJ(1));
t400 = sin(pkin(6));
t364 = t334 * t400;
t312 = t355 * t329;
t401 = sin(qJ(1));
t367 = t401 * t312;
t406 = -t401 * t313 - qJD(1) * t367 + (-qJD(1) * t314 - t306) * t334 - qJD(4) * t364;
t326 = qJ(4) + pkin(12);
t324 = sin(t326);
t325 = cos(t326);
t348 = -t334 * t312 + t401 * t314;
t284 = -t324 * t348 + t325 * t364;
t365 = t333 * t400;
t366 = t331 * t400;
t342 = -t327 * t365 - t328 * t366;
t300 = -t324 * t342 - t325 * t329;
t272 = atan2(-t284, t300);
t267 = sin(t272);
t268 = cos(t272);
t250 = -t267 * t284 + t268 * t300;
t248 = 0.1e1 / t250 ^ 2;
t349 = -t334 * t314 - t367;
t360 = t401 * t400;
t345 = -t324 * t349 + t325 * t360;
t283 = t345 ^ 2;
t246 = t248 * t283 + 0.1e1;
t290 = t324 * t360 + t325 * t349;
t341 = t348 * qJD(1) + t401 * t306 - t334 * t313;
t359 = qJD(1) * t364;
t254 = t290 * qJD(4) + t324 * t341 - t325 * t359;
t394 = t248 * t345;
t282 = t284 ^ 2;
t298 = 0.1e1 / t300 ^ 2;
t271 = t282 * t298 + 0.1e1;
t269 = 0.1e1 / t271;
t354 = qJD(1) * t360;
t378 = qJD(4) * t325;
t256 = t406 * t324 - t325 * t354 - t348 * t378;
t301 = t324 * t329 - t325 * t342;
t309 = -t327 * t366 + t328 * t365;
t305 = t309 * qJD(2);
t280 = t301 * qJD(4) + t305 * t324;
t297 = 0.1e1 / t300;
t385 = t284 * t298;
t353 = -t256 * t297 + t280 * t385;
t238 = t353 * t269;
t356 = -t267 * t300 - t268 * t284;
t233 = t356 * t238 - t267 * t256 + t268 * t280;
t247 = 0.1e1 / t250;
t249 = t247 * t248;
t398 = t233 * t249;
t376 = 0.2e1 * (-t254 * t394 - t283 * t398) / t246 ^ 2;
t405 = t280 * t298;
t292 = -t334 * t311 - t355 * t401;
t350 = -t292 * t297 + t309 * t385;
t404 = t324 * t350;
t257 = (qJD(4) * t348 + t354) * t324 + t406 * t325;
t332 = cos(qJ(6));
t295 = t401 * t311 - t334 * t355;
t330 = sin(qJ(6));
t383 = t295 * t330;
t266 = t290 * t332 - t383;
t260 = 0.1e1 / t266;
t261 = 0.1e1 / t266 ^ 2;
t403 = -0.2e1 * t284;
t402 = -0.2e1 * t345;
t387 = t297 * t405;
t397 = (t256 * t385 - t282 * t387) / t271 ^ 2;
t255 = t345 * qJD(4) + t324 * t359 + t325 * t341;
t307 = t329 * t313;
t347 = t314 * qJD(2);
t276 = t292 * qJD(1) - t401 * t307 - t334 * t347;
t382 = t295 * t332;
t265 = t290 * t330 + t382;
t377 = qJD(6) * t265;
t243 = t255 * t332 + t276 * t330 - t377;
t396 = t243 * t260 * t261;
t395 = t248 * t254;
t242 = t266 * qJD(6) + t255 * t330 - t276 * t332;
t259 = t265 ^ 2;
t253 = t259 * t261 + 0.1e1;
t391 = t261 * t265;
t393 = 0.1e1 / t253 ^ 2 * (t242 * t391 - t259 * t396);
t392 = t260 * t330;
t390 = t265 * t332;
t389 = t267 * t345;
t388 = t268 * t345;
t386 = t284 * t297;
t384 = t295 * t324;
t375 = -0.2e1 * t397;
t374 = t249 * t402;
t373 = -0.2e1 * t393;
t372 = 0.2e1 * t393;
t371 = t297 * t397;
t370 = t265 * t396;
t369 = t248 * t389;
t368 = t248 * t388;
t363 = 0.2e1 * t370;
t362 = t387 * t403;
t286 = -t324 * t364 - t325 * t348;
t357 = qJD(6) * t295 * t325 - t341;
t264 = -t286 * t332 + t292 * t330;
t263 = -t286 * t330 - t292 * t332;
t352 = t261 * t390 - t392;
t351 = -t286 * t297 + t301 * t385;
t346 = -t267 + (t268 * t386 + t267) * t269;
t343 = -qJD(4) * t384 + qJD(6) * t349 - t276 * t325;
t304 = t342 * qJD(2);
t281 = -t300 * qJD(4) + t305 * t325;
t278 = t295 * qJD(1) - t334 * t307 + t401 * t347;
t274 = t325 * t382 + t330 * t349;
t273 = t325 * t383 - t332 * t349;
t251 = 0.1e1 / t253;
t244 = 0.1e1 / t246;
t241 = t269 * t404;
t240 = t351 * t269;
t237 = t346 * t345;
t235 = (-t267 * t292 + t268 * t309) * t324 + t356 * t241;
t234 = t356 * t240 - t267 * t286 + t268 * t301;
t231 = t351 * t375 + (t301 * t362 - t257 * t297 + (t256 * t301 + t280 * t286 + t281 * t284) * t298) * t269;
t230 = t375 * t404 + (t350 * t378 + (t309 * t362 - t278 * t297 + (t256 * t309 + t280 * t292 + t284 * t304) * t298) * t324) * t269;
t1 = [t371 * t402 + (-t254 * t297 - t345 * t405) * t269, t230, 0, t231, 0, 0; t284 * t247 * t376 + (-t256 * t247 + (t233 * t284 + t237 * t254) * t248) * t244 - (-t237 * t248 * t376 + (-0.2e1 * t237 * t398 + (-t238 * t269 * t386 + t375) * t369 + (t371 * t403 - t238 + (t238 - t353) * t269) * t368 - t346 * t395) * t244) * t345 (-t235 * t394 - t247 * t384) * t376 + (-t235 * t395 + (-t276 * t324 + t295 * t378) * t247 + (t235 * t374 - t248 * t384) * t233 + (t309 * t378 - t230 * t284 - t241 * t256 + t304 * t324 + (-t241 * t300 - t292 * t324) * t238) * t368 + (-t292 * t378 - t230 * t300 - t241 * t280 - t278 * t324 + (t241 * t284 - t309 * t324) * t238) * t369) * t244, 0 (-t234 * t394 - t247 * t290) * t376 + (t234 * t233 * t374 + t255 * t247 + (-t290 * t233 - t234 * t254 + (-t231 * t284 - t240 * t256 + t281 + (-t240 * t300 - t286) * t238) * t388 + (-t231 * t300 - t240 * t280 - t257 + (t240 * t284 - t301) * t238) * t389) * t248) * t244, 0, 0; (-t260 * t263 + t264 * t391) * t372 + ((t264 * qJD(6) - t257 * t330 - t278 * t332) * t260 + t264 * t363 + (-t263 * t243 - (-t263 * qJD(6) - t257 * t332 + t278 * t330) * t265 - t264 * t242) * t261) * t251 (-t260 * t273 + t274 * t391) * t372 + (t274 * t363 + t357 * t260 * t332 + t343 * t392 + (t357 * t265 * t330 - t274 * t242 - t273 * t243 - t343 * t390) * t261) * t251, 0, -t352 * t345 * t373 + (t352 * t254 - ((-qJD(6) * t260 - 0.2e1 * t370) * t332 + (t242 * t332 + (t243 - t377) * t330) * t261) * t345) * t251, 0, t373 + 0.2e1 * (t242 * t261 * t251 + (-t251 * t396 - t261 * t393) * t265) * t265;];
JaD_rot  = t1;
