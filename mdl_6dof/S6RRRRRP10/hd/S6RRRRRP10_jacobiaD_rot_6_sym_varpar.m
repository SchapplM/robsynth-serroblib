% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRP10_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP10_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP10_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP10_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:45:01
% EndTime: 2019-02-26 22:45:04
% DurationCPUTime: 3.00s
% Computational Cost: add. (19416->183), mult. (35772->350), div. (1197->12), fcn. (45463->13), ass. (0->153)
t336 = qJ(4) + qJ(5);
t333 = sin(t336);
t338 = cos(pkin(6));
t341 = sin(qJ(1));
t340 = sin(qJ(2));
t402 = qJD(2) * t340;
t383 = t341 * t402;
t403 = qJD(1) * t341;
t385 = t340 * t403;
t343 = cos(qJ(2));
t344 = cos(qJ(1));
t404 = t343 * t344;
t305 = -t338 * t385 - t383 + (qJD(2) * t338 + qJD(1)) * t404;
t405 = t341 * t343;
t406 = t340 * t344;
t326 = t338 * t406 + t405;
t339 = sin(qJ(3));
t342 = cos(qJ(3));
t337 = sin(pkin(6));
t408 = t337 * t344;
t316 = t326 * t339 + t342 * t408;
t335 = qJD(4) + qJD(5);
t407 = t340 * t341;
t379 = -t338 * t404 + t407;
t371 = t379 * t335;
t387 = t337 * t403;
t363 = -qJD(3) * t316 + t305 * t342 + t339 * t387 + t371;
t449 = t333 * t363;
t334 = cos(t336);
t317 = -t326 * t342 + t339 * t408;
t365 = t338 * t405 + t406;
t351 = qJD(1) * t365 + qJD(2) * t326;
t447 = t317 * t335 + t351;
t267 = t447 * t333 + t363 * t334;
t443 = t317 * t333 + t334 * t379;
t286 = t443 ^ 2;
t410 = t337 * t342;
t325 = t338 * t339 + t340 * t410;
t409 = t337 * t343;
t309 = t325 * t333 + t334 * t409;
t307 = 0.1e1 / t309 ^ 2;
t273 = t286 * t307 + 0.1e1;
t271 = 0.1e1 / t273;
t412 = t334 * t335;
t266 = -t317 * t412 - t334 * t351 + t449;
t364 = -t325 * t335 + t337 * t402;
t411 = t337 * t339;
t324 = t338 * t342 - t340 * t411;
t384 = qJD(2) * t409;
t373 = qJD(3) * t324 - t335 * t409 + t342 * t384;
t277 = t333 * t373 - t334 * t364;
t306 = 0.1e1 / t309;
t415 = t443 * t307;
t369 = -t266 * t306 - t277 * t415;
t251 = t369 * t271;
t274 = atan2(t443, t309);
t269 = sin(t274);
t270 = cos(t274);
t370 = -t269 * t309 + t270 * t443;
t245 = t251 * t370 - t269 * t266 + t270 * t277;
t263 = t269 * t443 + t270 * t309;
t261 = 0.1e1 / t263 ^ 2;
t448 = t245 * t261;
t327 = -t338 * t407 + t404;
t319 = t327 * t342 + t341 * t411;
t296 = t319 * t333 - t334 * t365;
t431 = 0.2e1 * t296;
t260 = 0.1e1 / t263;
t444 = t260 * t448;
t380 = t431 * t444;
t303 = t338 * t383 + t385 + (-qJD(1) * t338 - qJD(2)) * t404;
t304 = -qJD(1) * t326 - qJD(2) * t365;
t318 = -t327 * t339 + t341 * t410;
t386 = qJD(1) * t408;
t281 = qJD(3) * t318 + t304 * t342 + t339 * t386;
t355 = t335 * t365 + t281;
t264 = t303 * t334 + t319 * t412 + t333 * t355;
t425 = t264 * t261;
t446 = -t425 + t380;
t297 = t319 * t334 + t333 * t365;
t289 = 0.1e1 / t297 ^ 2;
t311 = t318 ^ 2;
t418 = t289 * t311;
t279 = 0.1e1 + t418;
t280 = -qJD(3) * t319 - t304 * t339 + t342 * t386;
t265 = t355 * t334 + (-t319 * t335 - t303) * t333;
t288 = 0.1e1 / t297;
t424 = t265 * t288 * t289;
t392 = t311 * t424;
t417 = t289 * t318;
t428 = (t280 * t417 - t392) / t279 ^ 2;
t445 = 0.2e1 * t428;
t372 = t379 * t333;
t442 = t317 * t334 - t372;
t441 = -0.2e1 * t296;
t440 = t277 * t307;
t366 = t306 * t316 - t324 * t415;
t439 = t333 * t366;
t275 = 0.1e1 / t279;
t420 = t275 * t289;
t437 = t265 * t420 + t288 * t445;
t287 = t296 ^ 2;
t259 = t261 * t287 + 0.1e1;
t257 = 0.1e1 / t259;
t430 = (-t287 * t444 + t296 * t425) / t259 ^ 2;
t436 = -t257 * t448 - 0.2e1 * t260 * t430;
t378 = t417 * t428;
t388 = t318 * t424;
t435 = 0.2e1 * t275 * t388 - t280 * t420 + 0.2e1 * t378;
t400 = 0.2e1 * t430;
t426 = t261 * t296;
t434 = t446 * t257 + t400 * t426;
t282 = qJD(3) * t317 - t305 * t339 + t342 * t387;
t432 = 0.2e1 * t443;
t419 = t306 * t440;
t429 = (-t266 * t415 - t286 * t419) / t273 ^ 2;
t427 = t257 * t260;
t423 = t269 * t296;
t422 = t270 * t296;
t421 = t275 * t288;
t416 = t443 * t306;
t413 = t318 * t333;
t401 = qJD(3) * t339;
t399 = -0.2e1 * t429;
t396 = t306 * t429;
t394 = t257 * t426;
t391 = t275 * t417;
t381 = t419 * t432;
t310 = t325 * t334 - t333 * t409;
t368 = t306 * t442 - t310 * t415;
t298 = -t326 * t334 - t342 * t372;
t320 = (t333 * t342 * t343 - t334 * t340) * t337;
t367 = -t298 * t306 - t320 * t415;
t358 = t342 * t365;
t357 = -t269 + (-t270 * t416 + t269) * t271;
t356 = qJD(3) * t365;
t354 = -t335 * t358 - t304;
t353 = t303 * t342 + t327 * t335 + t339 * t356;
t312 = -qJD(3) * t325 - t339 * t384;
t285 = ((t335 * t342 - qJD(2)) * t343 * t334 + (-t343 * t401 + (-qJD(2) * t342 + t335) * t340) * t333) * t337;
t278 = t333 * t364 + t334 * t373;
t268 = (-t342 * t371 - t305) * t334 + (t326 * t335 - t342 * t351 + t379 * t401) * t333;
t256 = t271 * t439;
t255 = t367 * t271;
t254 = t368 * t271;
t249 = (t269 * t316 + t270 * t324) * t333 + t370 * t256;
t247 = t254 * t370 + t269 * t442 + t270 * t310;
t246 = t378 * t441 + (t388 * t441 + (t264 * t318 + t280 * t296) * t289) * t275;
t244 = t367 * t399 + (t320 * t381 - t268 * t306 + (t266 * t320 + t277 * t298 - t285 * t443) * t307) * t271;
t242 = t368 * t399 + (t310 * t381 - t267 * t306 + (t266 * t310 - t277 * t442 - t278 * t443) * t307) * t271;
t241 = t399 * t439 + (t366 * t412 + (t324 * t381 - t282 * t306 + (t266 * t324 - t277 * t316 - t312 * t443) * t307) * t333) * t271;
t240 = (t247 * t426 - t260 * t297) * t400 + (t247 * t380 + t265 * t260 + (-t297 * t245 - t247 * t264 - (t242 * t443 - t254 * t266 + t278 + (-t254 * t309 + t442) * t251) * t422 - (-t242 * t309 - t254 * t277 - t267 + (-t254 * t443 - t310) * t251) * t423) * t261) * t257;
t1 = [t396 * t431 + (-t264 * t306 + t296 * t440) * t271, t244, t241, t242, t242, 0; (t447 * t334 - t449) * t427 - (t357 * t264 + ((t251 * t271 * t416 + t399) * t269 + (t396 * t432 - t251 + (t251 - t369) * t271) * t270) * t296) * t394 + t436 * t443 + t434 * t357 * t296 (t333 * t353 + t334 * t354) * t427 - ((t244 * t443 - t255 * t266 + t285 + (-t255 * t309 - t298) * t251) * t270 + (-t244 * t309 - t255 * t277 - t268 + (-t255 * t443 - t320) * t251) * t269) * t394 + t436 * (-t327 * t334 - t333 * t358) + t434 * (t255 * t370 - t269 * t298 + t270 * t320) (t249 * t426 - t260 * t413) * t400 + ((t280 * t333 + t318 * t412) * t260 + t446 * t249 + (-t413 * t245 - (t324 * t412 + t241 * t443 - t256 * t266 + t312 * t333 + (-t256 * t309 + t316 * t333) * t251) * t422 - (t316 * t412 - t241 * t309 - t256 * t277 - t282 * t333 + (-t256 * t443 - t324 * t333) * t251) * t423) * t261) * t257, t240, t240, 0; t267 * t391 - t282 * t421 - t437 * t316 + t435 * t442 -(-t333 * t354 + t334 * t353) * t391 + (-t303 * t339 + t342 * t356) * t421 - t437 * t339 * t365 + t435 * (t327 * t333 - t334 * t358) (t288 * t319 + t334 * t418) * t445 + (0.2e1 * t334 * t392 - t281 * t288 + (-0.2e1 * t280 * t318 * t334 + t311 * t333 * t335 + t265 * t319) * t289) * t275, t246, t246, 0;];
JaD_rot  = t1;
