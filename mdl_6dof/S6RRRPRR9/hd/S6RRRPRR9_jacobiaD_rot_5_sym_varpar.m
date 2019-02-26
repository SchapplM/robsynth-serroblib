% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR9
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
% Datum: 2019-02-26 22:20
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR9_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:20:28
% EndTime: 2019-02-26 22:20:31
% DurationCPUTime: 2.72s
% Computational Cost: add. (13845->185), mult. (41275->344), div. (705->12), fcn. (52669->17), ass. (0->140)
t364 = sin(pkin(13));
t370 = sin(qJ(3));
t446 = cos(pkin(13));
t447 = cos(qJ(3));
t391 = -t370 * t364 + t447 * t446;
t452 = t391 * qJD(3);
t392 = t447 * t364 + t370 * t446;
t451 = qJD(3) * t392;
t367 = cos(pkin(7));
t345 = t391 * t367;
t368 = cos(pkin(6));
t374 = cos(qJ(2));
t375 = cos(qJ(1));
t421 = t374 * t375;
t371 = sin(qJ(2));
t372 = sin(qJ(1));
t424 = t371 * t372;
t349 = -t368 * t421 + t424;
t422 = t372 * t374;
t423 = t371 * t375;
t350 = t368 * t423 + t422;
t366 = sin(pkin(6));
t365 = sin(pkin(7));
t388 = t365 * t391;
t386 = t366 * t388;
t394 = -t345 * t349 - t350 * t392 - t375 * t386;
t308 = t394 ^ 2;
t427 = t366 * t374;
t429 = t366 * t371;
t324 = -t345 * t427 - t368 * t388 + t392 * t429;
t321 = 0.1e1 / t324 ^ 2;
t295 = t308 * t321 + 0.1e1;
t293 = 0.1e1 / t295;
t396 = t368 * t422 + t423;
t333 = t396 * qJD(1) + t350 * qJD(2);
t413 = t368 * t424;
t419 = qJD(2) * t371;
t334 = -qJD(1) * t413 - t372 * t419 + (qJD(2) * t368 + qJD(1)) * t421;
t343 = t367 * t451;
t385 = t365 * t451;
t383 = t366 * t385;
t384 = qJD(1) * t386;
t287 = -t333 * t345 - t334 * t392 + t349 * t343 - t350 * t452 + t372 * t384 + t375 * t383;
t306 = t366 * t345 * t419 + t452 * t429 + t368 * t385 + (qJD(2) * t392 + t343) * t427;
t320 = 0.1e1 / t324;
t432 = t394 * t321;
t401 = t287 * t320 - t306 * t432;
t268 = t401 * t293;
t296 = atan2(t394, t324);
t291 = sin(t296);
t292 = cos(t296);
t403 = -t291 * t324 + t292 * t394;
t263 = t403 * t268 + t287 * t291 + t292 * t306;
t278 = t291 * t394 + t292 * t324;
t276 = 0.1e1 / t278 ^ 2;
t450 = t263 * t276;
t449 = t306 * t321;
t275 = 0.1e1 / t278;
t428 = t366 * t372;
t337 = t365 * t396 + t367 * t428;
t369 = sin(qJ(5));
t373 = cos(qJ(5));
t344 = t392 * t365;
t346 = t392 * t367;
t395 = t413 - t421;
t389 = t344 * t428 - t346 * t396 - t391 * t395;
t304 = t337 * t369 + t373 * t389;
t298 = 0.1e1 / t304;
t299 = 0.1e1 / t304 ^ 2;
t317 = -t345 * t396 + t372 * t386 + t392 * t395;
t309 = t317 ^ 2;
t274 = t276 * t309 + 0.1e1;
t331 = t349 * qJD(1) + t395 * qJD(2);
t332 = t350 * qJD(1) + t396 * qJD(2);
t284 = t331 * t345 + t332 * t392 + t343 * t396 - t372 * t383 + t375 * t384 + t395 * t452;
t438 = t284 * t276;
t444 = t275 * t450;
t445 = (-t309 * t444 + t317 * t438) / t274 ^ 2;
t341 = t452 * t365;
t342 = t452 * t367;
t420 = qJD(1) * t344;
t285 = t331 * t346 - t332 * t391 - t396 * t342 + t395 * t451 + (t341 * t372 + t375 * t420) * t366;
t411 = qJD(1) * t366 * t367;
t328 = -t331 * t365 + t375 * t411;
t279 = t304 * qJD(5) + t285 * t369 - t328 * t373;
t303 = -t337 * t373 + t369 * t389;
t297 = t303 ^ 2;
t283 = t297 * t299 + 0.1e1;
t435 = t299 * t303;
t418 = qJD(5) * t303;
t280 = t285 * t373 + t328 * t369 - t418;
t439 = t280 * t298 * t299;
t443 = (t279 * t435 - t297 * t439) / t283 ^ 2;
t434 = t320 * t449;
t441 = (t287 * t432 - t308 * t434) / t295 ^ 2;
t440 = t276 * t317;
t437 = t291 * t317;
t436 = t292 * t317;
t433 = t394 * t320;
t431 = t365 * t369;
t430 = t365 * t373;
t426 = t366 * t375;
t416 = 0.2e1 * t445;
t415 = 0.2e1 * t443;
t414 = 0.2e1 * t441;
t408 = -0.2e1 * t320 * t441;
t407 = 0.2e1 * t303 * t439;
t406 = 0.2e1 * t394 * t434;
t405 = -0.2e1 * t317 * t444;
t336 = -t349 * t365 + t367 * t426;
t390 = t344 * t426 + t349 * t346 - t350 * t391;
t302 = t336 * t369 + t373 * t390;
t301 = -t336 * t373 + t369 * t390;
t400 = -t369 * t298 + t373 * t435;
t323 = t344 * t368 + (t346 * t374 + t371 * t391) * t366;
t399 = -t320 * t390 + t323 * t432;
t325 = -t345 * t350 + t349 * t392;
t330 = (t345 * t371 + t374 * t392) * t366;
t398 = -t320 * t325 + t330 * t432;
t327 = t346 * t395 - t391 * t396;
t397 = -t327 * t369 - t395 * t430;
t311 = t327 * t373 - t395 * t431;
t393 = t291 + (t292 * t433 - t291) * t293;
t382 = t333 * t346 - t334 * t391 + t349 * t342 + t350 * t451 + (t341 * t375 - t372 * t420) * t366;
t329 = -t333 * t365 - t372 * t411;
t326 = -t345 * t395 - t392 * t396;
t307 = (-t343 * t371 + t452 * t374 + (t345 * t374 - t371 * t392) * qJD(2)) * t366;
t305 = t341 * t368 + (t342 * t374 - t451 * t371 + (-t346 * t371 + t374 * t391) * qJD(2)) * t366;
t290 = t333 * t392 - t334 * t345 + t343 * t350 + t349 * t452;
t289 = t331 * t391 + t332 * t346 + t342 * t395 + t396 * t451;
t281 = 0.1e1 / t283;
t272 = 0.1e1 / t274;
t271 = t398 * t293;
t269 = t399 * t293;
t267 = t393 * t317;
t265 = -t403 * t271 + t291 * t325 + t292 * t330;
t264 = -t403 * t269 + t291 * t390 + t292 * t323;
t262 = t398 * t414 + (t330 * t406 + t290 * t320 + (-t287 * t330 - t306 * t325 - t307 * t394) * t321) * t293;
t260 = t399 * t414 + (t323 * t406 + t382 * t320 + (-t287 * t323 - t305 * t394 - t306 * t390) * t321) * t293;
t1 = [t317 * t408 + (t284 * t320 - t317 * t449) * t293, t262, t260, 0, 0, 0; -0.2e1 * (t267 * t440 + t275 * t394) * t445 + ((t438 + t405) * t267 + t287 * t275 - t394 * t450 + (t393 * t284 + ((-t268 * t293 * t433 + t414) * t291 + (t394 * t408 + t268 + (-t268 + t401) * t293) * t292) * t317) * t440) * t272 (-t265 * t440 - t275 * t326) * t416 + ((t331 * t392 - t332 * t345 + t343 * t395 - t396 * t452) * t275 + t265 * t405 + (-t326 * t263 + t265 * t284 + (t262 * t394 - t271 * t287 + t307 + (t271 * t324 + t325) * t268) * t436 + (-t262 * t324 + t271 * t306 + t290 + (t271 * t394 - t330) * t268) * t437) * t276) * t272 (-t264 * t440 - t275 * t389) * t416 + (t285 * t275 + t264 * t405 + (-t389 * t263 + t264 * t284 + (t260 * t394 - t269 * t287 + t305 + (t269 * t324 + t390) * t268) * t436 + (-t260 * t324 + t269 * t306 + t382 + (t269 * t394 - t323) * t268) * t437) * t276) * t272, 0, 0, 0; (-t298 * t301 + t302 * t435) * t415 + ((t302 * qJD(5) - t329 * t373 + t369 * t382) * t298 + t302 * t407 + (-t301 * t280 - (-t301 * qJD(5) + t329 * t369 + t373 * t382) * t303 - t302 * t279) * t299) * t281 (t298 * t397 + t311 * t435) * t415 + ((t311 * qJD(5) + t289 * t369 + t332 * t430) * t298 + t311 * t407 + (t397 * t280 - (t397 * qJD(5) + t289 * t373 - t332 * t431) * t303 - t311 * t279) * t299) * t281, t400 * t317 * t415 + (-t400 * t284 + ((qJD(5) * t298 + t407) * t373 + (-t279 * t373 + (-t280 + t418) * t369) * t299) * t317) * t281, 0, -0.2e1 * t443 + 0.2e1 * (t279 * t299 * t281 + (-t281 * t439 - t299 * t443) * t303) * t303, 0;];
JaD_rot  = t1;
