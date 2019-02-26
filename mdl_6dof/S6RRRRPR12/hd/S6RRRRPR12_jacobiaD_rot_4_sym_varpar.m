% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR12
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR12_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_jacobiaD_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:59
% EndTime: 2019-02-26 22:37:01
% DurationCPUTime: 1.98s
% Computational Cost: add. (8022->174), mult. (25165->321), div. (705->12), fcn. (31370->15), ass. (0->160)
t337 = cos(pkin(7));
t339 = sin(qJ(3));
t436 = cos(pkin(6));
t437 = sin(qJ(2));
t388 = t436 * t437;
t438 = sin(qJ(1));
t374 = t438 * t388;
t341 = cos(qJ(2));
t342 = cos(qJ(1));
t417 = t342 * t341;
t357 = t374 - t417;
t397 = t341 * t436;
t358 = t342 * t437 + t438 * t397;
t336 = sin(pkin(6));
t435 = sin(pkin(7));
t400 = t336 * t435;
t382 = t438 * t400;
t439 = cos(qJ(3));
t300 = -t357 * t439 + (-t358 * t337 + t382) * t339;
t420 = t336 * t337;
t392 = t438 * t420;
t319 = t358 * t435 + t392;
t338 = sin(qJ(4));
t340 = cos(qJ(4));
t280 = t300 * t338 - t319 * t340;
t447 = 0.2e1 * t280;
t391 = t437 * t439;
t418 = t339 * t341;
t364 = t337 * t391 + t418;
t367 = t337 * t418 + t391;
t384 = t436 * t435;
t381 = t339 * t384;
t290 = qJD(3) * t381 + (t364 * qJD(2) + t367 * qJD(3)) * t336;
t401 = t437 * t339;
t402 = t439 * t341;
t366 = t337 * t402 - t401;
t373 = t439 * t384;
t316 = -t366 * t336 - t373;
t314 = 0.1e1 / t316 ^ 2;
t446 = t290 * t314;
t335 = t438 * t437;
t380 = -t342 * t397 + t335;
t324 = t438 * t341 + t342 * t388;
t390 = t342 * t400;
t375 = t439 * t390;
t359 = -t324 * t339 - t375;
t363 = t380 * t439;
t361 = t337 * t363;
t294 = t361 - t359;
t292 = t294 ^ 2;
t286 = t292 * t314 + 0.1e1;
t284 = 0.1e1 / t286;
t371 = t380 * t339;
t404 = t324 * t439;
t353 = -t337 * t371 + t404;
t310 = t358 * qJD(1) + t324 * qJD(2);
t311 = -qJD(1) * t374 - qJD(2) * t335 + (qJD(2) * t436 + qJD(1)) * t417;
t331 = t339 * t390;
t369 = t439 * t382;
t403 = t337 * t439;
t362 = -qJD(1) * t369 - qJD(3) * t331 + t310 * t403 + t311 * t339;
t269 = t353 * qJD(3) + t362;
t313 = 0.1e1 / t316;
t422 = t294 * t314;
t379 = -t269 * t313 + t290 * t422;
t251 = t379 * t284;
t287 = atan2(-t294, t316);
t282 = sin(t287);
t283 = cos(t287);
t383 = -t282 * t316 - t283 * t294;
t246 = t383 * t251 - t282 * t269 + t283 * t290;
t263 = -t282 * t294 + t283 * t316;
t260 = 0.1e1 / t263;
t261 = 0.1e1 / t263 ^ 2;
t356 = t358 * t439;
t444 = -t337 * t356 + t339 * t357 + t369;
t293 = t444 ^ 2;
t257 = t293 * t261 + 0.1e1;
t255 = 0.1e1 / t257;
t429 = t255 * t261;
t309 = t324 * qJD(1) + t358 * qJD(2);
t351 = qJD(1) * t380 + t357 * qJD(2);
t349 = t351 * t439;
t267 = -qJD(1) * t375 + t300 * qJD(3) - t309 * t339 - t337 * t349;
t427 = t261 * t444;
t433 = t246 * t260 * t261;
t434 = (-t267 * t427 - t293 * t433) / t257 ^ 2;
t445 = -t246 * t429 - 0.2e1 * t260 * t434;
t440 = -0.2e1 * t444;
t395 = t433 * t440;
t415 = 0.2e1 * t434;
t443 = t255 * t395 - t267 * t429 - t415 * t427;
t442 = -(qJD(1) * t382 - t324 * qJD(3) - t310 * t337) * t339 + qJD(3) * t375 - t311 * t439;
t281 = t300 * t340 + t319 * t338;
t275 = 0.1e1 / t281;
t276 = 0.1e1 / t281 ^ 2;
t441 = -0.2e1 * t294;
t350 = t351 * t339;
t268 = qJD(1) * t331 + qJD(3) * t444 - t309 * t439 + t337 * t350;
t405 = t342 * t420;
t301 = qJD(1) * t405 - t351 * t435;
t258 = t281 * qJD(4) + t268 * t338 - t301 * t340;
t274 = t280 ^ 2;
t266 = t274 * t276 + 0.1e1;
t425 = t276 * t280;
t416 = qJD(4) * t280;
t259 = t268 * t340 + t301 * t338 - t416;
t428 = t259 * t275 * t276;
t432 = (t258 * t425 - t274 * t428) / t266 ^ 2;
t424 = t313 * t446;
t431 = (t269 * t422 - t292 * t424) / t286 ^ 2;
t430 = t255 * t260;
t264 = 0.1e1 / t266;
t426 = t264 * t276;
t423 = t294 * t313;
t419 = t337 * t339;
t414 = -0.2e1 * t432;
t413 = -0.2e1 * t431;
t411 = t276 * t432;
t410 = t313 * t431;
t408 = t255 * t427;
t407 = t258 * t426;
t406 = t280 * t428;
t399 = t338 * t435;
t398 = t340 * t435;
t394 = 0.2e1 * t406;
t393 = t424 * t441;
t372 = t337 * t380;
t354 = t339 * t372 - t404;
t298 = t331 + t354;
t318 = -t380 * t435 + t405;
t279 = t298 * t340 + t318 * t338;
t278 = t298 * t338 - t318 * t340;
t378 = -t338 * t275 + t340 * t425;
t296 = -t331 + t353;
t317 = t367 * t336 + t381;
t377 = -t296 * t313 + t317 * t422;
t306 = t324 * t403 - t371;
t323 = t364 * t336;
t376 = -t306 * t313 + t323 * t422;
t308 = t357 * t419 - t356;
t289 = t308 * t340 - t357 * t399;
t370 = -t308 * t338 - t357 * t398;
t368 = -t282 + (t283 * t423 + t282) * t284;
t365 = -t337 * t401 + t402;
t360 = t439 * t372;
t307 = -t358 * t339 - t357 * t403;
t303 = (t366 * qJD(2) + t365 * qJD(3)) * t336;
t302 = -qJD(1) * t392 - t310 * t435;
t291 = qJD(3) * t373 + (t365 * qJD(2) + t366 * qJD(3)) * t336;
t273 = t311 * t403 - t310 * t339 + (-t324 * t419 - t363) * qJD(3);
t272 = -t307 * qJD(3) + t309 * t419 + t349;
t271 = qJD(3) * t360 + t442;
t270 = -qJD(3) * t361 - t442;
t254 = t376 * t284;
t253 = t377 * t284;
t247 = t383 * t253 - t282 * t296 + t283 * t317;
t245 = t376 * t413 + (t323 * t393 - t273 * t313 + (t269 * t323 + t290 * t306 + t294 * t303) * t314) * t284;
t244 = t377 * t413 + (t317 * t393 - t270 * t313 + (t269 * t317 + t290 * t296 + t291 * t294) * t314) * t284;
t1 = [t410 * t440 + (-t267 * t313 - t444 * t446) * t284, t245, t244, 0, 0, 0; (t354 * qJD(3) - t362) * t430 + (t368 * t267 - ((-t251 * t284 * t423 + t413) * t282 + (t410 * t441 - t251 + (t251 - t379) * t284) * t283) * t444) * t408 + t445 * (-t360 + t359) - t443 * t368 * t444 (t308 * qJD(3) - t309 * t403 + t350) * t430 + ((-t245 * t294 - t254 * t269 + t303 + (-t254 * t316 - t306) * t251) * t283 + (-t245 * t316 - t254 * t290 - t273 + (t254 * t294 - t323) * t251) * t282) * t408 + t445 * t307 + t443 * (t383 * t254 - t282 * t306 + t283 * t323) (-t247 * t427 - t260 * t300) * t415 + (t247 * t395 + t268 * t260 + (-t300 * t246 - t247 * t267 - (-(-t244 * t294 - t253 * t269 + t291 + (-t253 * t316 - t296) * t251) * t283 - (-t244 * t316 - t253 * t290 - t270 + (t253 * t294 - t317) * t251) * t282) * t444) * t261) * t255, 0, 0, 0; 0.2e1 * (-t275 * t278 + t279 * t425) * t432 + ((t279 * qJD(4) + t271 * t338 - t302 * t340) * t275 + t279 * t394 + (-t278 * t259 - (-t278 * qJD(4) + t271 * t340 + t302 * t338) * t280 - t279 * t258) * t276) * t264 (t411 * t447 - t407) * t289 - (-t259 * t426 + t275 * t414) * t370 + ((t289 * qJD(4) + t272 * t338 + t309 * t398) * t275 - (t370 * qJD(4) + t272 * t340 - t309 * t399) * t425 + t289 * t394) * t264, -t378 * t444 * t414 + (t378 * t267 - ((-qJD(4) * t275 - 0.2e1 * t406) * t340 + (t258 * t340 + (t259 - t416) * t338) * t276) * t444) * t264, t414 + (t407 + (-t264 * t428 - t411) * t280) * t447, 0, 0;];
JaD_rot  = t1;
