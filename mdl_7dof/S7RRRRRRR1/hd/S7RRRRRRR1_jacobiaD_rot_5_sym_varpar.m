% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JaD_rot [3x7]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S7RRRRRRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [7x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:54:21
% EndTime: 2019-02-26 22:54:23
% DurationCPUTime: 2.39s
% Computational Cost: add. (6638->201), mult. (19789->386), div. (992->12), fcn. (24542->13), ass. (0->158)
t307 = sin(qJ(4));
t308 = sin(qJ(3));
t309 = sin(qJ(2));
t375 = qJD(3) * t309;
t314 = cos(qJ(2));
t313 = cos(qJ(3));
t344 = qJD(2) * t313 - qJD(4);
t411 = t344 * t314;
t417 = (-t308 * t375 + t411) * t307;
t310 = sin(qJ(1));
t377 = qJD(2) * t314;
t315 = cos(qJ(1));
t382 = t309 * t315;
t416 = qJD(1) * t382 + t310 * t377;
t312 = cos(qJ(4));
t345 = qJD(4) * t313 - qJD(2);
t414 = (qJD(3) * t308 * t312 + t307 * t345) * t309 - t312 * t411;
t380 = t313 * t314;
t386 = t308 * t315;
t291 = t310 * t380 + t386;
t384 = t309 * t312;
t274 = t291 * t307 - t310 * t384;
t383 = t309 * t313;
t288 = t307 * t383 + t312 * t314;
t267 = atan2(-t274, t288);
t255 = sin(t267);
t256 = cos(t267);
t236 = -t255 * t274 + t256 * t288;
t234 = 0.1e1 / t236 ^ 2;
t379 = t315 * t313;
t381 = t310 * t308;
t294 = t314 * t379 - t381;
t279 = t294 * t307 - t312 * t382;
t273 = t279 ^ 2;
t232 = t234 * t273 + 0.1e1;
t346 = qJD(1) * t314 + qJD(3);
t376 = qJD(2) * t315;
t323 = t309 * t376 + t310 * t346;
t374 = qJD(3) * t314;
t347 = qJD(1) + t374;
t263 = t313 * t323 + t347 * t386;
t280 = t294 * t312 + t307 * t382;
t378 = qJD(1) * t310;
t329 = t309 * t378 - t314 * t376;
t240 = qJD(4) * t280 - t263 * t307 + t312 * t329;
t400 = t240 * t234;
t272 = t274 ^ 2;
t286 = 0.1e1 / t288 ^ 2;
t266 = t272 * t286 + 0.1e1;
t257 = 0.1e1 / t266;
t352 = t308 * t374;
t385 = t309 * t310;
t357 = qJD(2) * t385;
t265 = -t308 * t378 - t310 * t352 - t313 * t357 + t346 * t379;
t276 = t291 * t312 + t307 * t385;
t242 = qJD(4) * t276 + t265 * t307 - t312 * t416;
t339 = t345 * t312;
t260 = t309 * t339 + t417;
t285 = 0.1e1 / t288;
t390 = t274 * t286;
t335 = -t242 * t285 + t260 * t390;
t223 = t335 * t257;
t340 = -t255 * t288 - t256 * t274;
t217 = t223 * t340 - t242 * t255 + t256 * t260;
t233 = 0.1e1 / t236;
t235 = t233 * t234;
t405 = t217 * t235;
t370 = 0.2e1 * (-t273 * t405 + t279 * t400) / t232 ^ 2;
t413 = t260 * t286;
t290 = t314 * t381 - t379;
t387 = t308 * t309;
t361 = t274 * t387;
t330 = -t285 * t290 + t286 * t361;
t412 = t307 * t330;
t289 = -t307 * t314 + t312 * t383;
t283 = t289 * t315;
t351 = t313 * t375;
t409 = -qJD(5) * t283 + t308 * t329 - t315 * t351;
t373 = qJD(4) * t309;
t243 = (-qJD(4) * t291 + t416) * t307 + (t310 * t373 + t265) * t312;
t306 = sin(qJ(5));
t311 = cos(qJ(5));
t331 = t310 * t313 + t314 * t386;
t254 = t280 * t311 - t306 * t331;
t248 = 0.1e1 / t254;
t249 = 0.1e1 / t254 ^ 2;
t408 = -0.2e1 * t274;
t407 = 0.2e1 * t279;
t241 = (t315 * t373 - t263) * t312 + (-qJD(4) * t294 - t329) * t307;
t262 = t308 * t323 - t347 * t379;
t225 = qJD(5) * t254 + t241 * t306 - t262 * t311;
t253 = t280 * t306 + t311 * t331;
t247 = t253 ^ 2;
t239 = t247 * t249 + 0.1e1;
t397 = t249 * t253;
t371 = qJD(5) * t253;
t226 = t241 * t311 + t262 * t306 - t371;
t402 = t226 * t248 * t249;
t404 = (t225 * t397 - t247 * t402) / t239 ^ 2;
t392 = t285 * t413;
t403 = (t242 * t390 - t272 * t392) / t266 ^ 2;
t401 = t234 * t279;
t399 = t248 * t306;
t398 = t248 * t311;
t396 = t253 * t306;
t395 = t253 * t311;
t394 = t255 * t279;
t393 = t256 * t279;
t391 = t274 * t285;
t389 = t331 * t307;
t388 = t331 * t312;
t372 = qJD(4) * t312;
t369 = -0.2e1 * t404;
t368 = 0.2e1 * t404;
t367 = -0.2e1 * t403;
t366 = t235 * t407;
t365 = t285 * t403;
t364 = t253 * t402;
t363 = t234 * t394;
t362 = t234 * t393;
t360 = t308 * t382;
t350 = t217 * t366;
t349 = 0.2e1 * t364;
t348 = t392 * t408;
t341 = -qJD(5) * t388 - t263;
t252 = -t276 * t311 + t290 * t306;
t251 = -t276 * t306 - t290 * t311;
t336 = qJD(5) * t360 + t289 * t378 + t414 * t315;
t334 = t249 * t395 - t399;
t333 = -t276 * t285 + t289 * t390;
t281 = t288 * t310;
t292 = t307 * t380 - t384;
t332 = t281 * t285 + t292 * t390;
t328 = -t308 * t377 - t351;
t326 = -t255 + (t256 * t391 + t255) * t257;
t325 = qJD(1) * t288;
t324 = qJD(4) * t389 - qJD(5) * t294 + t262 * t312;
t282 = t288 * t315;
t271 = -t283 * t311 + t306 * t360;
t270 = -t283 * t306 - t311 * t360;
t269 = -t294 * t306 - t311 * t388;
t268 = t294 * t311 - t306 * t388;
t264 = t331 * qJD(1) + t291 * qJD(3) - t308 * t357;
t259 = t314 * t339 + (-t344 * t309 - t352) * t307;
t246 = t260 * t310 + t315 * t325;
t237 = 0.1e1 / t239;
t230 = 0.1e1 / t232;
t229 = t257 * t412;
t228 = t332 * t257;
t227 = t333 * t257;
t222 = t326 * t279;
t221 = (t255 * t290 - t256 * t387) * t307 - t340 * t229;
t219 = t228 * t340 + t255 * t281 + t256 * t292;
t218 = t227 * t340 - t255 * t276 + t256 * t289;
t216 = t332 * t367 + (t292 * t348 + t246 * t285 + (t242 * t292 + t259 * t274 - t260 * t281) * t286) * t257;
t214 = t333 * t367 + (t289 * t348 - t243 * t285 + (t242 * t289 + t260 * t276 - t274 * t414) * t286) * t257;
t213 = 0.2e1 * t403 * t412 + (-t330 * t372 + (0.2e1 * t361 * t392 + t264 * t285 + (-t242 * t387 - t260 * t290 + t274 * t328) * t286) * t307) * t257;
t1 = [t365 * t407 + (-t240 * t285 + t279 * t413) * t257, t216, t213, t214, 0, 0, 0; t274 * t233 * t370 + (-t242 * t233 + (t217 * t274 - t222 * t240) * t234) * t230 + (t222 * t234 * t370 + (0.2e1 * t222 * t405 - (-t223 * t257 * t391 + t367) * t363 - (t365 * t408 - t223 + (t223 - t335) * t257) * t362 - t326 * t400) * t230) * t279 (t219 * t401 + t233 * t282) * t370 + (t219 * t350 + (t282 * t217 - t219 * t240 - (-t216 * t274 - t228 * t242 + t259 + (-t228 * t288 + t281) * t223) * t393 - (-t216 * t288 - t228 * t260 + t246 + (t228 * t274 - t292) * t223) * t394) * t234 + (t310 * t325 + (-t345 * t384 - t417) * t315) * t233) * t230 (t221 * t401 + t233 * t389) * t370 + (-t221 * t400 + (t262 * t307 - t331 * t372) * t233 + (t221 * t366 + t234 * t389) * t217 - (t290 * t372 - t213 * t288 + t229 * t260 + t264 * t307 + (-t229 * t274 + t307 * t387) * t223) * t363 - (-t372 * t387 - t213 * t274 - (-t223 * t288 - t242) * t229 + (t223 * t290 + t328) * t307) * t362) * t230 (t218 * t401 - t233 * t280) * t370 + (t218 * t350 + t241 * t233 + (-t280 * t217 - t218 * t240 - (-t214 * t274 - t227 * t242 - t414 + (-t227 * t288 - t276) * t223) * t393 - (-t214 * t288 - t227 * t260 - t243 + (t227 * t274 - t289) * t223) * t394) * t234) * t230, 0, 0, 0; (-t248 * t251 + t252 * t397) * t368 + ((qJD(5) * t252 - t243 * t306 - t264 * t311) * t248 + t252 * t349 + (-t251 * t226 - (-t251 * qJD(5) - t243 * t311 + t264 * t306) * t253 - t252 * t225) * t249) * t237 (-t248 * t270 + t271 * t397) * t368 + (t271 * t349 + t336 * t399 + t409 * t398 + (-t271 * t225 - t270 * t226 - t336 * t395 + t396 * t409) * t249) * t237 (-t248 * t268 + t269 * t397) * t368 + (t269 * t349 + t341 * t398 + t324 * t399 + (-t269 * t225 - t268 * t226 - t324 * t395 + t341 * t396) * t249) * t237, t334 * t279 * t369 + (t334 * t240 + ((-qJD(5) * t248 - 0.2e1 * t364) * t311 + (t225 * t311 + (t226 - t371) * t306) * t249) * t279) * t237, t369 + 0.2e1 * (t225 * t237 * t249 + (-t237 * t402 - t249 * t404) * t253) * t253, 0, 0;];
JaD_rot  = t1;
