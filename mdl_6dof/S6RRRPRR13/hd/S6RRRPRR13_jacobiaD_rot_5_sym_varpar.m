% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR13_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:23:05
% EndTime: 2019-02-26 22:23:07
% DurationCPUTime: 2.06s
% Computational Cost: add. (8441->175), mult. (25165->321), div. (705->12), fcn. (31370->15), ass. (0->161)
t348 = cos(pkin(7));
t349 = sin(qJ(3));
t445 = cos(pkin(6));
t446 = sin(qJ(2));
t397 = t445 * t446;
t447 = sin(qJ(1));
t383 = t447 * t397;
t350 = cos(qJ(2));
t351 = cos(qJ(1));
t426 = t351 * t350;
t366 = t383 - t426;
t406 = t350 * t445;
t367 = t351 * t446 + t447 * t406;
t347 = sin(pkin(6));
t444 = sin(pkin(7));
t407 = t347 * t444;
t391 = t447 * t407;
t448 = cos(qJ(3));
t308 = -t366 * t448 + (-t367 * t348 + t391) * t349;
t429 = t347 * t348;
t401 = t447 * t429;
t327 = t367 * t444 + t401;
t346 = pkin(13) + qJ(5);
t344 = sin(t346);
t345 = cos(t346);
t288 = t308 * t344 - t327 * t345;
t456 = 0.2e1 * t288;
t400 = t448 * t446;
t427 = t349 * t350;
t374 = t348 * t400 + t427;
t376 = t348 * t427 + t400;
t393 = t445 * t444;
t390 = t349 * t393;
t298 = qJD(3) * t390 + (t374 * qJD(2) + t376 * qJD(3)) * t347;
t410 = t446 * t349;
t411 = t448 * t350;
t375 = t348 * t411 - t410;
t382 = t448 * t393;
t324 = -t375 * t347 - t382;
t322 = 0.1e1 / t324 ^ 2;
t455 = t298 * t322;
t343 = t447 * t446;
t389 = -t351 * t406 + t343;
t332 = t447 * t350 + t351 * t397;
t399 = t351 * t407;
t384 = t448 * t399;
t368 = -t332 * t349 - t384;
t372 = t389 * t448;
t370 = t348 * t372;
t302 = t370 - t368;
t300 = t302 ^ 2;
t294 = t300 * t322 + 0.1e1;
t292 = 0.1e1 / t294;
t380 = t389 * t349;
t413 = t332 * t448;
t362 = -t348 * t380 + t413;
t318 = t367 * qJD(1) + t332 * qJD(2);
t319 = -qJD(1) * t383 - qJD(2) * t343 + (qJD(2) * t445 + qJD(1)) * t426;
t339 = t349 * t399;
t378 = t448 * t391;
t412 = t348 * t448;
t371 = -qJD(1) * t378 - qJD(3) * t339 + t318 * t412 + t319 * t349;
t277 = t362 * qJD(3) + t371;
t321 = 0.1e1 / t324;
t431 = t302 * t322;
t388 = -t277 * t321 + t298 * t431;
t259 = t388 * t292;
t295 = atan2(-t302, t324);
t290 = sin(t295);
t291 = cos(t295);
t392 = -t290 * t324 - t291 * t302;
t254 = t392 * t259 - t290 * t277 + t291 * t298;
t271 = -t290 * t302 + t291 * t324;
t268 = 0.1e1 / t271;
t269 = 0.1e1 / t271 ^ 2;
t365 = t367 * t448;
t453 = -t348 * t365 + t349 * t366 + t378;
t301 = t453 ^ 2;
t265 = t269 * t301 + 0.1e1;
t263 = 0.1e1 / t265;
t439 = t263 * t269;
t317 = t332 * qJD(1) + t367 * qJD(2);
t360 = qJD(1) * t389 + t366 * qJD(2);
t358 = t360 * t448;
t275 = -qJD(1) * t384 + t308 * qJD(3) - t317 * t349 - t348 * t358;
t437 = t269 * t453;
t442 = t254 * t268 * t269;
t443 = (-t275 * t437 - t301 * t442) / t265 ^ 2;
t454 = -t254 * t439 - 0.2e1 * t268 * t443;
t449 = -0.2e1 * t453;
t404 = t442 * t449;
t424 = 0.2e1 * t443;
t452 = t263 * t404 - t275 * t439 - t424 * t437;
t451 = -(qJD(1) * t391 - t332 * qJD(3) - t318 * t348) * t349 + qJD(3) * t384 - t319 * t448;
t289 = t308 * t345 + t327 * t344;
t283 = 0.1e1 / t289;
t284 = 0.1e1 / t289 ^ 2;
t450 = -0.2e1 * t302;
t433 = t321 * t455;
t441 = (t277 * t431 - t300 * t433) / t294 ^ 2;
t440 = t263 * t268;
t359 = t360 * t349;
t276 = qJD(1) * t339 + qJD(3) * t453 - t317 * t448 + t348 * t359;
t414 = t351 * t429;
t309 = qJD(1) * t414 - t360 * t444;
t425 = qJD(5) * t288;
t267 = t276 * t345 + t309 * t344 - t425;
t438 = t267 * t283 * t284;
t282 = t288 ^ 2;
t274 = t282 * t284 + 0.1e1;
t272 = 0.1e1 / t274;
t436 = t272 * t284;
t266 = t289 * qJD(5) + t276 * t344 - t309 * t345;
t434 = t284 * t288;
t435 = 0.1e1 / t274 ^ 2 * (t266 * t434 - t282 * t438);
t432 = t302 * t321;
t428 = t348 * t349;
t423 = -0.2e1 * t441;
t422 = -0.2e1 * t435;
t420 = t284 * t435;
t419 = t321 * t441;
t417 = t263 * t437;
t416 = t266 * t436;
t415 = t288 * t438;
t409 = t344 * t444;
t408 = t345 * t444;
t403 = 0.2e1 * t415;
t402 = t433 * t450;
t381 = t348 * t389;
t363 = t349 * t381 - t413;
t306 = t339 + t363;
t326 = -t389 * t444 + t414;
t287 = t306 * t345 + t326 * t344;
t286 = t306 * t344 - t326 * t345;
t387 = -t283 * t344 + t345 * t434;
t304 = -t339 + t362;
t325 = t376 * t347 + t390;
t386 = -t304 * t321 + t325 * t431;
t314 = t332 * t412 - t380;
t331 = t374 * t347;
t385 = -t314 * t321 + t331 * t431;
t316 = t366 * t428 - t365;
t297 = t316 * t345 - t366 * t409;
t379 = -t316 * t344 - t366 * t408;
t377 = -t290 + (t291 * t432 + t290) * t292;
t373 = -t348 * t410 + t411;
t369 = t448 * t381;
t315 = -t367 * t349 - t366 * t412;
t311 = (t375 * qJD(2) + t373 * qJD(3)) * t347;
t310 = -qJD(1) * t401 - t318 * t444;
t299 = qJD(3) * t382 + (t373 * qJD(2) + t375 * qJD(3)) * t347;
t281 = t319 * t412 - t318 * t349 + (-t332 * t428 - t372) * qJD(3);
t280 = -t315 * qJD(3) + t317 * t428 + t358;
t279 = qJD(3) * t369 + t451;
t278 = -qJD(3) * t370 - t451;
t262 = t385 * t292;
t261 = t386 * t292;
t255 = t392 * t261 - t290 * t304 + t291 * t325;
t253 = t385 * t423 + (t331 * t402 - t281 * t321 + (t277 * t331 + t298 * t314 + t302 * t311) * t322) * t292;
t252 = t386 * t423 + (t325 * t402 - t278 * t321 + (t277 * t325 + t298 * t304 + t299 * t302) * t322) * t292;
t1 = [t419 * t449 + (-t275 * t321 - t453 * t455) * t292, t253, t252, 0, 0, 0; (t363 * qJD(3) - t371) * t440 + (t377 * t275 - ((-t259 * t292 * t432 + t423) * t290 + (t419 * t450 - t259 + (t259 - t388) * t292) * t291) * t453) * t417 + t454 * (-t369 + t368) - t452 * t377 * t453 (t316 * qJD(3) - t317 * t412 + t359) * t440 + ((-t253 * t302 - t262 * t277 + t311 + (-t262 * t324 - t314) * t259) * t291 + (-t253 * t324 - t262 * t298 - t281 + (t262 * t302 - t331) * t259) * t290) * t417 + t454 * t315 + t452 * (t392 * t262 - t290 * t314 + t291 * t331) (-t255 * t437 - t268 * t308) * t424 + (t255 * t404 + t276 * t268 + (-t308 * t254 - t255 * t275 - (-(-t252 * t302 - t261 * t277 + t299 + (-t261 * t324 - t304) * t259) * t291 - (-t252 * t324 - t261 * t298 - t278 + (t261 * t302 - t325) * t259) * t290) * t453) * t269) * t263, 0, 0, 0; 0.2e1 * (-t283 * t286 + t287 * t434) * t435 + ((t287 * qJD(5) + t279 * t344 - t310 * t345) * t283 + t287 * t403 + (-t286 * t267 - (-t286 * qJD(5) + t279 * t345 + t310 * t344) * t288 - t287 * t266) * t284) * t272 (t420 * t456 - t416) * t297 - (-t267 * t436 + t283 * t422) * t379 + ((t297 * qJD(5) + t280 * t344 + t317 * t408) * t283 - (t379 * qJD(5) + t280 * t345 - t317 * t409) * t434 + t297 * t403) * t272, -t387 * t453 * t422 + (t387 * t275 - ((-qJD(5) * t283 - 0.2e1 * t415) * t345 + (t266 * t345 + (t267 - t425) * t344) * t284) * t453) * t272, 0, t422 + (t416 + (-t272 * t438 - t420) * t288) * t456, 0;];
JaD_rot  = t1;
