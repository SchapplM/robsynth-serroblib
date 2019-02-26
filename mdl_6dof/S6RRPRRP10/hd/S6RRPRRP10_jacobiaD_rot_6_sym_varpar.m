% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:51:04
% EndTime: 2019-02-26 21:51:07
% DurationCPUTime: 2.86s
% Computational Cost: add. (15414->179), mult. (28394->346), div. (951->12), fcn. (36086->13), ass. (0->148)
t313 = cos(pkin(6));
t316 = sin(qJ(1));
t315 = sin(qJ(2));
t374 = qJD(2) * t315;
t352 = t316 * t374;
t375 = qJD(1) * t316;
t355 = t315 * t375;
t318 = cos(qJ(2));
t319 = cos(qJ(1));
t376 = t318 * t319;
t281 = -t313 * t355 - t352 + (qJD(2) * t313 + qJD(1)) * t376;
t378 = t316 * t318;
t379 = t315 * t319;
t302 = t313 * t379 + t378;
t311 = pkin(11) + qJ(4);
t309 = sin(t311);
t310 = cos(t311);
t312 = sin(pkin(6));
t382 = t312 * t319;
t292 = t302 * t309 + t310 * t382;
t357 = t312 * t375;
t258 = qJD(4) * t292 - t281 * t310 - t309 * t357;
t293 = -t302 * t310 + t309 * t382;
t314 = sin(qJ(5));
t317 = cos(qJ(5));
t380 = t315 * t316;
t348 = -t313 * t376 + t380;
t270 = t293 * t314 + t317 * t348;
t336 = t313 * t378 + t379;
t326 = qJD(1) * t336 + qJD(2) * t302;
t243 = qJD(5) * t270 - t258 * t317 + t314 * t326;
t342 = t348 * t314;
t271 = t293 * t317 - t342;
t418 = qJD(5) * t271 + t258 * t314 + t317 * t326;
t262 = t270 ^ 2;
t384 = t312 * t315;
t298 = t309 * t313 + t310 * t384;
t377 = t317 * t318;
t288 = t298 * t314 + t312 * t377;
t285 = 0.1e1 / t288 ^ 2;
t251 = t262 * t285 + 0.1e1;
t247 = 0.1e1 / t251;
t297 = -t309 * t384 + t310 * t313;
t353 = qJD(2) * t312 * t318;
t283 = qJD(4) * t297 + t310 * t353;
t381 = t314 * t318;
t289 = t298 * t317 - t312 * t381;
t354 = t312 * t374;
t259 = qJD(5) * t289 + t283 * t314 - t317 * t354;
t284 = 0.1e1 / t288;
t386 = t270 * t285;
t340 = -t259 * t386 + t284 * t418;
t227 = t340 * t247;
t252 = atan2(t270, t288);
t245 = sin(t252);
t246 = cos(t252);
t341 = -t245 * t288 + t246 * t270;
t222 = t227 * t341 + t245 * t418 + t246 * t259;
t239 = t245 * t270 + t246 * t288;
t237 = 0.1e1 / t239 ^ 2;
t417 = t222 * t237;
t303 = -t313 * t380 + t376;
t383 = t312 * t316;
t295 = t303 * t310 + t309 * t383;
t272 = t295 * t314 - t317 * t336;
t402 = 0.2e1 * t272;
t236 = 0.1e1 / t239;
t412 = t236 * t417;
t349 = t402 * t412;
t280 = -qJD(1) * t302 - qJD(2) * t336;
t294 = -t303 * t309 + t310 * t383;
t356 = qJD(1) * t382;
t255 = qJD(4) * t294 + t280 * t310 + t309 * t356;
t273 = t295 * t317 + t314 * t336;
t279 = t313 * t352 + t355 + (-qJD(1) * t313 - qJD(2)) * t376;
t240 = qJD(5) * t273 + t255 * t314 + t279 * t317;
t396 = t240 * t237;
t416 = -t396 + t349;
t265 = 0.1e1 / t273 ^ 2;
t287 = t294 ^ 2;
t389 = t265 * t287;
t253 = 0.1e1 + t389;
t254 = -qJD(4) * t295 - t280 * t309 + t310 * t356;
t241 = -qJD(5) * t272 + t255 * t317 - t279 * t314;
t264 = 0.1e1 / t273;
t395 = t241 * t264 * t265;
t362 = t287 * t395;
t388 = t265 * t294;
t400 = (t254 * t388 - t362) / t253 ^ 2;
t413 = 0.2e1 * t400;
t411 = -0.2e1 * t272;
t410 = t259 * t285;
t337 = t284 * t292 - t297 * t386;
t409 = t314 * t337;
t249 = 0.1e1 / t253;
t391 = t249 * t265;
t407 = t241 * t391 + t264 * t413;
t263 = t272 ^ 2;
t235 = t237 * t263 + 0.1e1;
t233 = 0.1e1 / t235;
t401 = (-t263 * t412 + t272 * t396) / t235 ^ 2;
t406 = -t233 * t417 - 0.2e1 * t236 * t401;
t347 = t388 * t400;
t358 = t294 * t395;
t405 = 0.2e1 * t249 * t358 - t254 * t391 + 0.2e1 * t347;
t370 = 0.2e1 * t401;
t397 = t237 * t272;
t404 = t233 * t416 + t370 * t397;
t256 = qJD(4) * t293 - t281 * t309 + t310 * t357;
t403 = 0.2e1 * t270;
t390 = t284 * t410;
t399 = (-t262 * t390 + t386 * t418) / t251 ^ 2;
t398 = t233 * t236;
t394 = t245 * t272;
t393 = t246 * t272;
t392 = t249 * t264;
t387 = t270 * t284;
t385 = t294 * t314;
t373 = qJD(4) * t309;
t372 = qJD(5) * t310;
t371 = qJD(5) * t317;
t369 = -0.2e1 * t399;
t365 = t284 * t399;
t364 = t233 * t397;
t361 = t249 * t388;
t350 = t390 * t403;
t339 = t271 * t284 - t289 * t386;
t274 = -t302 * t317 - t310 * t342;
t296 = (t310 * t381 - t315 * t317) * t312;
t338 = -t274 * t284 - t296 * t386;
t332 = t310 * t336;
t331 = -t245 + (-t246 * t387 + t245) * t247;
t330 = qJD(4) * t336;
t329 = -qJD(5) * t332 - t280;
t328 = qJD(5) * t303 + t279 * t310 + t309 * t330;
t282 = -qJD(4) * t298 - t309 * t353;
t261 = ((-qJD(2) + t372) * t377 + (-t318 * t373 + (-qJD(2) * t310 + qJD(5)) * t315) * t314) * t312;
t260 = -qJD(5) * t288 + t283 * t317 + t314 * t354;
t244 = (-t348 * t372 - t281) * t317 + (t302 * qJD(5) - t310 * t326 + t348 * t373) * t314;
t232 = t247 * t409;
t231 = t338 * t247;
t230 = t339 * t247;
t225 = (t245 * t292 + t246 * t297) * t314 + t341 * t232;
t223 = t230 * t341 + t245 * t271 + t246 * t289;
t221 = t338 * t369 + (t296 * t350 - t244 * t284 + (t259 * t274 - t261 * t270 - t296 * t418) * t285) * t247;
t219 = t339 * t369 + (t289 * t350 - t243 * t284 + (-t259 * t271 - t260 * t270 - t289 * t418) * t285) * t247;
t218 = t369 * t409 + (t337 * t371 + (t297 * t350 - t256 * t284 + (-t259 * t292 - t270 * t282 - t297 * t418) * t285) * t314) * t247;
t1 = [t365 * t402 + (-t240 * t284 + t272 * t410) * t247, t221, 0, t218, t219, 0; t418 * t398 - (t331 * t240 + ((t227 * t247 * t387 + t369) * t245 + (t365 * t403 - t227 + (t227 - t340) * t247) * t246) * t272) * t364 + t406 * t270 + t404 * t331 * t272 (t314 * t328 + t317 * t329) * t398 - ((t221 * t270 + t231 * t418 + t261 + (-t231 * t288 - t274) * t227) * t246 + (-t221 * t288 - t231 * t259 - t244 + (-t231 * t270 - t296) * t227) * t245) * t364 + t406 * (-t303 * t317 - t314 * t332) + t404 * (t231 * t341 - t245 * t274 + t246 * t296) 0 (t225 * t397 - t236 * t385) * t370 + ((t254 * t314 + t294 * t371) * t236 + t416 * t225 + (-t385 * t222 - (t297 * t371 + t218 * t270 + t232 * t418 + t282 * t314 + (-t232 * t288 + t292 * t314) * t227) * t393 - (t292 * t371 - t218 * t288 - t232 * t259 - t256 * t314 + (-t232 * t270 - t297 * t314) * t227) * t394) * t237) * t233 (t223 * t397 - t236 * t273) * t370 + (t223 * t349 + t241 * t236 + (-t273 * t222 - t223 * t240 - (t219 * t270 + t230 * t418 + t260 + (-t230 * t288 + t271) * t227) * t393 - (-t219 * t288 - t230 * t259 - t243 + (-t230 * t270 - t289) * t227) * t394) * t237) * t233, 0; t243 * t361 - t256 * t392 + t271 * t405 - t292 * t407 -(-t314 * t329 + t317 * t328) * t361 + (-t279 * t309 + t310 * t330) * t392 - t407 * t309 * t336 + t405 * (t303 * t314 - t317 * t332) 0 (t264 * t295 + t317 * t389) * t413 + (0.2e1 * t317 * t362 - t255 * t264 + (qJD(5) * t287 * t314 - 0.2e1 * t254 * t294 * t317 + t241 * t295) * t265) * t249, t347 * t411 + (t358 * t411 + (t240 * t294 + t254 * t272) * t265) * t249, 0;];
JaD_rot  = t1;
