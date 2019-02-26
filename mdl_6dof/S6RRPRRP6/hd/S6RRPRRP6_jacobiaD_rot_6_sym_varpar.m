% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP6
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
% Datum: 2019-02-26 21:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRRP6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:48:51
% EndTime: 2019-02-26 21:48:54
% DurationCPUTime: 3.62s
% Computational Cost: add. (18619->185), mult. (52466->359), div. (951->12), fcn. (68108->15), ass. (0->157)
t457 = sin(pkin(11));
t459 = cos(pkin(6));
t399 = t459 * t457;
t458 = cos(pkin(11));
t400 = t459 * t458;
t460 = sin(qJ(2));
t461 = cos(qJ(2));
t348 = t461 * t399 + t460 * t400;
t353 = t460 * t457 - t461 * t458;
t470 = t460 * t399 - t461 * t400;
t343 = t470 * qJD(2);
t384 = t461 * t457 + t460 * t458;
t352 = t384 * qJD(2);
t362 = sin(qJ(1));
t462 = cos(qJ(1));
t393 = -t462 * t343 - t362 * t352;
t416 = qJD(1) * t462;
t436 = qJD(1) * t362;
t314 = -t348 * t436 - t353 * t416 + t393;
t361 = sin(qJ(4));
t364 = cos(qJ(4));
t392 = -t462 * t348 + t362 * t353;
t359 = sin(pkin(6));
t418 = t359 * t462;
t325 = -t361 * t392 + t364 * t418;
t417 = t359 * t436;
t287 = t325 * qJD(4) - t314 * t364 - t361 * t417;
t324 = -t361 * t418 - t364 * t392;
t360 = sin(qJ(5));
t363 = cos(qJ(5));
t437 = -t362 * t384 - t462 * t470;
t414 = t437 * t363;
t299 = t324 * t360 + t414;
t378 = t362 * t470 - t462 * t384;
t379 = t348 * qJD(2);
t383 = t353 * qJD(2);
t374 = t378 * qJD(1) + t362 * t383 - t462 * t379;
t372 = t374 * t360;
t481 = t299 * qJD(5) + t287 * t363 + t372;
t415 = t437 * t360;
t301 = t324 * t363 - t415;
t274 = t301 * qJD(5) - t287 * t360 + t374 * t363;
t294 = t299 ^ 2;
t346 = t353 * t359;
t347 = t384 * t359;
t390 = -t347 * t364 - t459 * t361;
t320 = -t346 * t363 - t360 * t390;
t318 = 0.1e1 / t320 ^ 2;
t281 = t294 * t318 + 0.1e1;
t279 = 0.1e1 / t281;
t336 = -t347 * t361 + t459 * t364;
t342 = t359 * t383;
t316 = t336 * qJD(4) - t342 * t364;
t439 = t346 * t360;
t321 = -t363 * t390 + t439;
t341 = qJD(2) * t347;
t291 = t321 * qJD(5) + t316 * t360 - t341 * t363;
t317 = 0.1e1 / t320;
t441 = t299 * t318;
t397 = -t274 * t317 + t291 * t441;
t259 = t397 * t279;
t282 = atan2(-t299, t320);
t277 = sin(t282);
t278 = cos(t282);
t398 = -t277 * t320 - t278 * t299;
t254 = t398 * t259 - t274 * t277 + t278 * t291;
t271 = -t277 * t299 + t278 * t320;
t269 = 0.1e1 / t271 ^ 2;
t480 = t254 * t269;
t391 = -t362 * t348 - t462 * t353;
t438 = t359 * t362;
t328 = t361 * t438 + t364 * t391;
t304 = t328 * t360 + t363 * t378;
t463 = 0.2e1 * t304;
t268 = 0.1e1 / t271;
t475 = t268 * t480;
t411 = t463 * t475;
t327 = -t361 * t391 + t364 * t438;
t382 = t392 * qJD(1) + t362 * t343 - t462 * t352;
t409 = t359 * t416;
t284 = t327 * qJD(4) + t361 * t409 + t364 * t382;
t305 = t328 * t363 - t360 * t378;
t312 = qJD(1) * t437 - t362 * t379 - t462 * t383;
t272 = qJD(5) * t305 + t284 * t360 - t312 * t363;
t451 = t272 * t269;
t479 = -t451 + t411;
t283 = -t328 * qJD(4) - t361 * t382 + t364 * t409;
t297 = 0.1e1 / t305 ^ 2;
t322 = t327 ^ 2;
t444 = t297 * t322;
t290 = 0.1e1 + t444;
t273 = -qJD(5) * t304 + t284 * t363 + t312 * t360;
t296 = 0.1e1 / t305;
t450 = t273 * t296 * t297;
t423 = t322 * t450;
t443 = t297 * t327;
t455 = (t283 * t443 - t423) / t290 ^ 2;
t476 = 0.2e1 * t455;
t474 = -0.2e1 * t304;
t473 = t291 * t318;
t394 = t317 * t325 + t336 * t441;
t472 = t360 * t394;
t288 = 0.1e1 / t290;
t446 = t288 * t297;
t469 = t273 * t446 + t296 * t476;
t295 = t304 ^ 2;
t267 = t269 * t295 + 0.1e1;
t265 = 0.1e1 / t267;
t456 = (-t295 * t475 + t304 * t451) / t267 ^ 2;
t468 = -t265 * t480 - 0.2e1 * t268 * t456;
t410 = t443 * t455;
t419 = t327 * t450;
t467 = -t283 * t446 + 0.2e1 * t288 * t419 + 0.2e1 * t410;
t431 = 0.2e1 * t456;
t452 = t269 * t304;
t466 = t265 * t479 + t431 * t452;
t465 = t324 * qJD(4) + t314 * t361 - t364 * t417;
t464 = -0.2e1 * t299;
t445 = t317 * t473;
t454 = (t274 * t441 - t294 * t445) / t281 ^ 2;
t453 = t265 * t268;
t449 = t277 * t304;
t448 = t278 * t304;
t447 = t288 * t296;
t442 = t299 * t317;
t440 = t327 * t360;
t435 = qJD(4) * t361;
t434 = qJD(5) * t360;
t433 = qJD(5) * t363;
t432 = t364 * qJD(5);
t430 = -0.2e1 * t454;
t427 = t317 * t454;
t426 = t265 * t452;
t421 = t288 * t443;
t412 = t445 * t464;
t396 = -t301 * t317 + t321 * t441;
t308 = t392 * t363 + t364 * t415;
t329 = -t347 * t363 - t364 * t439;
t395 = -t308 * t317 + t329 * t441;
t385 = -t277 + (t278 * t442 + t277) * t279;
t376 = t364 * t378;
t375 = qJD(4) * t378;
t373 = qJD(5) * t376 - t382;
t371 = qJD(5) * t391 - t312 * t364 - t361 * t375;
t315 = t390 * qJD(4) + t342 * t361;
t293 = (-t346 * t432 + t342) * t363 + (qJD(5) * t347 - t341 * t364 + t346 * t435) * t360;
t292 = -t320 * qJD(5) + t316 * t363 + t341 * t360;
t276 = t364 * t372 - t415 * t435 + t414 * t432 - (t391 * qJD(1) + t393) * t363 - t392 * t434;
t264 = t279 * t472;
t263 = t395 * t279;
t262 = t396 * t279;
t257 = (t277 * t325 + t278 * t336) * t360 + t398 * t264;
t255 = t398 * t262 - t277 * t301 + t278 * t321;
t253 = t395 * t430 + (t329 * t412 - t276 * t317 + (t274 * t329 + t291 * t308 + t293 * t299) * t318) * t279;
t251 = t396 * t430 + (t321 * t412 + t481 * t317 + (t274 * t321 + t291 * t301 + t292 * t299) * t318) * t279;
t250 = t430 * t472 + (t394 * t433 + (t336 * t412 + t465 * t317 + (t274 * t336 - t291 * t325 + t299 * t315) * t318) * t360) * t279;
t1 = [t427 * t463 + (-t272 * t317 + t304 * t473) * t279, t253, 0, t250, t251, 0; -t274 * t453 - (t385 * t272 + ((-t259 * t279 * t442 + t430) * t277 + (t427 * t464 - t259 + (t259 - t397) * t279) * t278) * t304) * t426 - t468 * t299 + t466 * t385 * t304 (t360 * t371 + t363 * t373) * t453 - ((-t253 * t299 - t263 * t274 + t293 + (-t263 * t320 - t308) * t259) * t278 + (-t253 * t320 - t263 * t291 - t276 + (t263 * t299 - t329) * t259) * t277) * t426 + t468 * (t360 * t376 - t363 * t391) + t466 * (t398 * t263 - t277 * t308 + t278 * t329) 0 (t257 * t452 - t268 * t440) * t431 + ((t283 * t360 + t327 * t433) * t268 + t479 * t257 + (-t440 * t254 - (t336 * t433 - t250 * t299 - t264 * t274 + t315 * t360 + (-t264 * t320 + t325 * t360) * t259) * t448 - (t325 * t433 - t250 * t320 - t264 * t291 + t465 * t360 + (t264 * t299 - t336 * t360) * t259) * t449) * t269) * t265 (t255 * t452 - t268 * t305) * t431 + (t255 * t411 + t273 * t268 + (-t305 * t254 - t255 * t272 - (-t251 * t299 - t262 * t274 + t292 + (-t262 * t320 - t301) * t259) * t448 - (-t251 * t320 - t262 * t291 + t481 + (t262 * t299 - t321) * t259) * t449) * t269) * t265, 0; -t467 * t301 - t469 * t325 - t481 * t421 + t465 * t447 -(-t360 * t373 + t363 * t371) * t421 + (t312 * t361 - t364 * t375) * t447 + t469 * t361 * t378 + t467 * (t360 * t391 + t363 * t376) 0 (t296 * t328 + t363 * t444) * t476 + (0.2e1 * t363 * t423 - t284 * t296 + (-0.2e1 * t283 * t327 * t363 + t273 * t328 + t322 * t434) * t297) * t288, t410 * t474 + (t419 * t474 + (t272 * t327 + t283 * t304) * t297) * t288, 0;];
JaD_rot  = t1;
