% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR11_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:40
% EndTime: 2019-02-26 21:06:43
% DurationCPUTime: 3.24s
% Computational Cost: add. (12717->155), mult. (38126->312), div. (702->12), fcn. (49058->17), ass. (0->144)
t433 = sin(pkin(12));
t438 = cos(pkin(6));
t394 = t438 * t433;
t436 = cos(pkin(12));
t439 = sin(qJ(1));
t441 = cos(qJ(1));
t337 = -t439 * t394 + t441 * t436;
t334 = t337 * qJD(1);
t347 = sin(qJ(3));
t348 = cos(qJ(3));
t336 = t441 * t394 + t439 * t436;
t434 = sin(pkin(7));
t435 = sin(pkin(6));
t392 = t435 * t434;
t381 = t441 * t392;
t437 = cos(pkin(7));
t395 = t438 * t436;
t448 = -t441 * t395 + t439 * t433;
t449 = t448 * t437;
t358 = t449 + t381;
t444 = t336 * t347 + t348 * t358;
t367 = t439 * t395 + t441 * t433;
t333 = t367 * qJD(1);
t378 = t439 * t392;
t455 = qJD(1) * t378 - t333 * t437;
t297 = t444 * qJD(3) - t334 * t348 - t455 * t347;
t346 = sin(qJ(4));
t393 = t437 * t435;
t379 = t439 * t393;
t368 = qJD(1) * t379 + t333 * t434;
t440 = cos(qJ(4));
t414 = t336 * t348;
t317 = t347 * t358 - t414;
t359 = -t441 * t393 + t434 * t448;
t463 = -t317 * t346 - t359 * t440;
t275 = qJD(4) * t463 + t297 * t440 - t368 * t346;
t306 = t317 * t440 - t359 * t346;
t465 = t306 * qJD(4) + t297 * t346 + t368 * t440;
t300 = t463 ^ 2;
t364 = t436 * t393 + t434 * t438;
t391 = t435 * t433;
t328 = t347 * t364 + t348 * t391;
t335 = -t436 * t392 + t438 * t437;
t384 = -t328 * t346 + t335 * t440;
t310 = 0.1e1 / t384 ^ 2;
t288 = t300 * t310 + 0.1e1;
t286 = 0.1e1 / t288;
t313 = t328 * t440 + t335 * t346;
t327 = -t347 * t391 + t348 * t364;
t322 = t327 * qJD(3);
t298 = qJD(4) * t313 + t322 * t346;
t309 = 0.1e1 / t384;
t419 = t463 * t310;
t388 = t298 * t419 - t309 * t465;
t255 = t388 * t286;
t289 = atan2(-t463, -t384);
t284 = sin(t289);
t285 = cos(t289);
t390 = t284 * t384 - t285 * t463;
t250 = t255 * t390 + t284 * t465 + t285 * t298;
t265 = -t284 * t463 - t285 * t384;
t263 = 0.1e1 / t265 ^ 2;
t462 = t250 * t263;
t459 = -t367 * t437 + t378;
t262 = 0.1e1 / t265;
t458 = t262 * t462;
t319 = t337 * t348 + t459 * t347;
t357 = t367 * t434 + t379;
t355 = t357 * t440;
t307 = t319 * t346 - t355;
t442 = 0.2e1 * t307;
t397 = t442 * t458;
t332 = t336 * qJD(1);
t411 = qJD(3) * t347;
t446 = t459 * t348;
t447 = t358 * qJD(1);
t293 = t446 * qJD(3) - t332 * t348 - t337 * t411 + t447 * t347;
t308 = t319 * t440 + t357 * t346;
t356 = t359 * qJD(1);
t271 = t308 * qJD(4) + t293 * t346 + t440 * t356;
t426 = t271 * t263;
t454 = -t426 + t397;
t453 = (t347 * t449 - t414) * qJD(3) - t334 * t347 + t455 * t348 + t381 * t411;
t301 = t307 ^ 2;
t261 = t301 * t263 + 0.1e1;
t409 = 0.2e1 * (-t301 * t458 + t307 * t426) / t261 ^ 2;
t410 = qJD(4) * t346;
t272 = qJD(4) * t355 + t293 * t440 - t319 * t410 - t346 * t356;
t292 = qJD(3) * t319 - t332 * t347 - t447 * t348;
t344 = sin(pkin(13));
t345 = cos(pkin(13));
t267 = t272 * t345 + t292 * t344;
t318 = t337 * t347 - t446;
t283 = t308 * t345 + t318 * t344;
t278 = 0.1e1 / t283 ^ 2;
t452 = t267 * t278;
t451 = t298 * t310;
t386 = -t309 * t444 + t327 * t419;
t450 = t346 * t386;
t401 = qJD(4) * t440;
t277 = 0.1e1 / t283;
t443 = -0.2e1 * t463;
t266 = t272 * t344 - t292 * t345;
t282 = t308 * t344 - t318 * t345;
t276 = t282 ^ 2;
t270 = t276 * t278 + 0.1e1;
t425 = t278 * t282;
t427 = t277 * t452;
t431 = (t266 * t425 - t276 * t427) / t270 ^ 2;
t421 = t309 * t451;
t429 = (t300 * t421 - t419 * t465) / t288 ^ 2;
t428 = t263 * t307;
t402 = t318 * t440;
t291 = t319 * t344 - t345 * t402;
t424 = t278 * t291;
t423 = t284 * t307;
t422 = t285 * t307;
t420 = t463 * t309;
t417 = t318 * t346;
t413 = t344 * t277;
t412 = t345 * t282;
t408 = 0.2e1 * t431;
t407 = -0.2e1 * t429;
t405 = t309 * t429;
t404 = t282 * t427;
t399 = 0.2e1 * t404;
t398 = t421 * t443;
t387 = -t306 * t309 + t313 * t419;
t377 = -t440 * t292 + t318 * t410;
t376 = -t284 + (-t285 * t420 + t284) * t286;
t323 = t328 * qJD(3);
t299 = t384 * qJD(4) + t322 * t440;
t290 = -t319 * t345 - t344 * t402;
t281 = t306 * t345 - t344 * t444;
t280 = t306 * t344 + t345 * t444;
t268 = 0.1e1 / t270;
t259 = 0.1e1 / t261;
t258 = t286 * t450;
t256 = t387 * t286;
t252 = (t284 * t444 + t285 * t327) * t346 + t390 * t258;
t251 = t256 * t390 + t284 * t306 + t285 * t313;
t248 = t387 * t407 + (-t313 * t398 - t275 * t309 + (-t298 * t306 + t299 * t463 - t313 * t465) * t310) * t286;
t247 = t407 * t450 + ((-t327 * t398 + t453 * t309 + (-t298 * t444 - t323 * t463 - t327 * t465) * t310) * t346 + t386 * t401) * t286;
t1 = [-t405 * t442 + (t271 * t309 + t307 * t451) * t286, 0, t247, t248, 0, 0; t463 * t262 * t409 + (t465 * t262 + t463 * t462 - (t376 * t271 + ((t255 * t286 * t420 + t407) * t284 + (-t405 * t443 - t255 + (t255 - t388) * t286) * t285) * t307) * t428) * t259 + (t454 * t259 + t428 * t409) * t376 * t307, 0 (t252 * t428 + t262 * t417) * t409 + ((-t292 * t346 - t318 * t401) * t262 + t454 * t252 + (t417 * t250 - (t327 * t401 - t247 * t463 + t258 * t465 - t323 * t346 + (t258 * t384 + t346 * t444) * t255) * t422 - (t444 * t401 + t247 * t384 - t258 * t298 - t453 * t346 + (t258 * t463 - t327 * t346) * t255) * t423) * t263) * t259 (t251 * t428 - t262 * t308) * t409 + (t251 * t397 + t272 * t262 + (-t308 * t250 - t251 * t271 - (-t248 * t463 + t256 * t465 + t299 + (t256 * t384 + t306) * t255) * t422 - (t248 * t384 - t256 * t298 + t275 + (t256 * t463 - t313) * t255) * t423) * t263) * t259, 0, 0; (-t277 * t280 + t281 * t425) * t408 + ((t275 * t344 - t345 * t453) * t277 + t281 * t399 + (-t280 * t267 - (t275 * t345 + t344 * t453) * t282 - t281 * t266) * t278) * t268, 0, -0.2e1 * t290 * t277 * t431 + t282 * t408 * t424 + ((-t293 * t345 + t344 * t377) * t277 - t290 * t452 - (t293 * t344 + t345 * t377) * t425 - t266 * t424 + t291 * t399) * t268 (-t278 * t412 + t413) * t307 * t408 + (-0.2e1 * t307 * t345 * t404 - t271 * t413 + (t271 * t412 + (t266 * t345 + t267 * t344) * t307) * t278) * t268, 0, 0;];
JaD_rot  = t1;
