% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR11_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR11_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR11_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:21:51
% EndTime: 2019-02-26 22:21:55
% DurationCPUTime: 3.66s
% Computational Cost: add. (14567->183), mult. (41112->340), div. (968->12), fcn. (52257->15), ass. (0->149)
t413 = sin(qJ(2));
t508 = cos(pkin(6));
t509 = sin(qJ(1));
t460 = t508 * t509;
t406 = t413 * t460;
t417 = cos(qJ(2));
t510 = cos(qJ(1));
t461 = t508 * t510;
t471 = t509 * t413;
t385 = -qJD(1) * t406 - qJD(2) * t471 + (t510 * qJD(1) + qJD(2) * t461) * t417;
t412 = sin(qJ(3));
t416 = cos(qJ(3));
t433 = t413 * t461 + t509 * t417;
t507 = sin(pkin(6));
t459 = t510 * t507;
t426 = t433 * t412 + t416 * t459;
t458 = t509 * t507;
t445 = qJD(1) * t458;
t357 = t426 * qJD(3) - t385 * t416 - t412 * t445;
t390 = -t412 * t459 + t433 * t416;
t411 = sin(qJ(5));
t415 = cos(qJ(5));
t366 = t390 * t415 + t426 * t411;
t427 = t390 * qJD(3) + t385 * t412 - t416 * t445;
t327 = t366 * qJD(5) - t357 * t411 - t427 * t415;
t470 = t412 * t507;
t397 = t413 * t470 - t508 * t416;
t469 = t416 * t507;
t398 = t508 * t412 + t413 * t469;
t381 = t397 * t411 + t398 * t415;
t457 = qJD(2) * t417 * t507;
t388 = t398 * qJD(3) + t412 * t457;
t389 = -t397 * qJD(3) + t416 * t457;
t338 = t381 * qJD(5) - t388 * t415 + t389 * t411;
t451 = t397 * t415 - t398 * t411;
t339 = t451 * qJD(5) + t388 * t411 + t389 * t415;
t538 = t390 * t411 - t426 * t415;
t361 = t538 ^ 2;
t377 = 0.1e1 / t451 ^ 2;
t342 = t361 * t377 + 0.1e1;
t340 = 0.1e1 / t342;
t376 = 0.1e1 / t451;
t488 = t538 * t377;
t441 = t366 * t376 + t381 * t488;
t524 = t338 * t377;
t494 = t376 * t524;
t511 = -0.2e1 * t538;
t464 = t494 * t511;
t503 = (t327 * t488 + t361 * t494) / t342 ^ 2;
t477 = -0.2e1 * t503;
t541 = qJD(5) * t538 + t357 * t415 - t427 * t411;
t299 = (t541 * t376 - (t327 * t381 + t338 * t366 + t339 * t538) * t377 + t381 * t464) * t340 - t441 * t477;
t444 = t327 * t376 + t338 * t488;
t309 = t444 * t340;
t543 = t441 * t340;
t549 = t299 * t451 - (t538 * t543 - t381) * t309 + t543 * t338 - t541;
t548 = t299 * t538 + (t451 * t543 - t366) * t309 - t543 * t327 + t339;
t343 = atan2(-t538, -t451);
t335 = sin(t343);
t336 = cos(t343);
t454 = t335 * t451 - t336 * t538;
t305 = -t335 * t366 + t336 * t381 + t454 * t543;
t432 = t510 * t413 + t417 * t460;
t383 = t433 * qJD(1) + t432 * qJD(2);
t455 = t510 * t417 - t406;
t393 = t412 * t458 + t455 * t416;
t446 = qJD(1) * t459;
t425 = t393 * qJD(3) - t383 * t412 - t416 * t446;
t436 = t455 * t412 - t416 * t458;
t354 = -t436 * qJD(3) - t383 * t416 + t412 * t446;
t491 = t354 * t415;
t513 = t393 * t411 - t436 * t415;
t325 = -qJD(5) * t513 + t425 * t411 + t491;
t372 = t393 * t415 + t436 * t411;
t410 = sin(qJ(6));
t414 = cos(qJ(6));
t351 = t372 * t414 - t410 * t432;
t399 = -t417 * t461 + t471;
t382 = t399 * qJD(1) - t455 * qJD(2);
t317 = t351 * qJD(6) + t325 * t410 - t382 * t414;
t350 = t372 * t410 + t414 * t432;
t481 = qJD(6) * t350;
t318 = t325 * t414 + t382 * t410 - t481;
t345 = 0.1e1 / t351;
t346 = 0.1e1 / t351 ^ 2;
t501 = t318 * t345 * t346;
t512 = 0.2e1 * t350;
t465 = t501 * t512;
t534 = (((-t318 + t481) * t410 - t317 * t414) * t346 + (qJD(6) * t345 + t465) * t414) * t513;
t493 = t346 * t350;
t443 = -t345 * t410 + t414 * t493;
t527 = t443 * t513;
t303 = t454 * t309 - t335 * t327 + t336 * t338;
t322 = -t335 * t538 - t336 * t451;
t320 = 0.1e1 / t322 ^ 2;
t525 = t303 * t320;
t449 = t411 * t412 + t415 * t416;
t521 = t449 * t432;
t518 = -t411 * t416 + t412 * t415;
t520 = t518 * t432;
t519 = t393 * qJD(5) - t425;
t516 = qJD(5) - qJD(3);
t319 = 0.1e1 / t322;
t362 = t513 ^ 2;
t316 = t320 * t362 + 0.1e1;
t492 = t354 * t411;
t324 = t372 * qJD(5) - t425 * t415 + t492;
t499 = t320 * t513;
t505 = t319 * t525;
t506 = (t324 * t499 - t362 * t505) / t316 ^ 2;
t500 = t320 * t324;
t344 = t350 ^ 2;
t332 = t344 * t346 + 0.1e1;
t330 = 0.1e1 / t332;
t498 = t330 * t346;
t497 = 0.1e1 / t332 ^ 2 * (t317 * t493 - t344 * t501);
t496 = t335 * t513;
t495 = t336 * t513;
t489 = t538 * t376;
t479 = 0.2e1 * t506;
t478 = 0.2e1 * t505;
t476 = -0.2e1 * t497;
t475 = 0.2e1 * t497;
t474 = t346 * t497;
t473 = t376 * t503;
t472 = t317 * t498;
t468 = -0.2e1 * t319 * t506;
t467 = t320 * t479;
t466 = t513 * t478;
t349 = -t366 * t414 + t399 * t410;
t348 = -t366 * t410 - t399 * t414;
t373 = t518 * t399;
t430 = t411 * t469 - t415 * t470;
t394 = t430 * t417;
t440 = t373 * t376 + t394 * t488;
t437 = -t335 + (-t336 * t489 + t335) * t340;
t359 = -t410 * t521 + t455 * t414;
t360 = -t455 * t410 - t414 * t521;
t435 = qJD(5) * t436;
t384 = t432 * qJD(1) + t433 * qJD(2);
t358 = -t430 * t413 * qJD(2) + t516 * t417 * (t411 * t470 + t415 * t469);
t334 = -t516 * t399 * t449 + t384 * t518;
t333 = t449 * t382 - t516 * t520;
t314 = 0.1e1 / t316;
t313 = t440 * t340;
t308 = t437 * t513;
t306 = t454 * t313 - t335 * t373 + t336 * t394;
t302 = t440 * t477 + (-t394 * t464 + t334 * t376 + (t327 * t394 + t338 * t373 + t358 * t538) * t377) * t340;
t1 = [-0.2e1 * t513 * t473 + (t324 * t376 + t513 * t524) * t340, t302, t299, 0, -t299, 0; -t538 * t468 + (-t327 * t319 + (t303 * t538 - t308 * t324) * t320) * t314 + (t308 * t467 + (t308 * t478 - t437 * t500 + (-(t309 * t340 * t489 + t477) * t496 - (-t473 * t511 - t309 + (t309 - t444) * t340) * t495) * t320) * t314) * t513 (t306 * t499 - t319 * t520) * t479 + (t306 * t466 + (-t520 * t303 - t306 * t324 - (-t302 * t538 - t313 * t327 + t358 + (t313 * t451 - t373) * t309) * t495 - (t302 * t451 - t313 * t338 - t334 + (t313 * t538 - t394) * t309) * t496) * t320 + (-t382 * t518 - t516 * t521) * t319) * t314, -t372 * t468 + ((t519 * t411 - t415 * t435 - t491) * t319 + t372 * t525 - (t549 * t335 - t548 * t336) * t499) * t314 - (t513 * t467 + (-t500 + t466) * t314) * t305, 0 (t305 * t499 - t319 * t372) * t479 + (t305 * t466 + t325 * t319 + (-t372 * t303 - t305 * t324 - t548 * t495 + t549 * t496) * t320) * t314, 0; (-t345 * t348 + t349 * t493) * t475 + ((t349 * qJD(6) - t384 * t414 + t410 * t541) * t345 + t349 * t465 + (-t348 * t318 - (-t348 * qJD(6) + t384 * t410 + t414 * t541) * t350 - t349 * t317) * t346) * t330 (t474 * t512 - t472) * t360 + (-t318 * t498 + t345 * t476) * t359 + ((t360 * qJD(6) + t333 * t410 - t383 * t414) * t345 - (-t359 * qJD(6) + t333 * t414 + t383 * t410) * t493 + t360 * t465) * t330, t475 * t527 + (-t443 * (t411 * t435 + t519 * t415 + t492) + t534) * t330, 0, t476 * t527 + (t443 * t324 - t534) * t330, t476 + (t472 + (-t330 * t501 - t474) * t350) * t512;];
JaD_rot  = t1;
