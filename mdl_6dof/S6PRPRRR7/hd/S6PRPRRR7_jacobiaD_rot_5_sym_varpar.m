% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR7
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRPRRR7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_jacobiaD_rot_5_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:57:31
% EndTime: 2019-02-26 19:57:35
% DurationCPUTime: 2.90s
% Computational Cost: add. (15417->170), mult. (47519->336), div. (538->12), fcn. (61232->19), ass. (0->149)
t411 = sin(pkin(14));
t412 = sin(pkin(8));
t413 = sin(pkin(7));
t415 = cos(pkin(14));
t416 = cos(pkin(8));
t417 = cos(pkin(7));
t498 = cos(pkin(13));
t499 = cos(pkin(6));
t457 = t499 * t498;
t497 = sin(pkin(13));
t500 = sin(qJ(2));
t501 = cos(qJ(2));
t443 = -t457 * t501 + t497 * t500;
t414 = sin(pkin(6));
t463 = t414 * t498;
t438 = -t413 * t463 - t417 * t443;
t444 = t457 * t500 + t497 * t501;
t511 = (-t411 * t444 + t415 * t438) * t416 + (t413 * t443 - t417 * t463) * t412;
t456 = t499 * t497;
t409 = -t456 * t500 + t498 * t501;
t445 = t456 * t501 + t498 * t500;
t462 = t414 * t497;
t439 = t413 * t462 - t417 * t445;
t433 = t409 * t411 - t415 * t439;
t440 = t413 * t445 + t417 * t462;
t509 = -t412 * t440 + t416 * t433;
t464 = t500 * t415;
t467 = t501 * t411;
t483 = t412 * t413;
t507 = (t500 * t483 + (-t417 * t464 - t467) * t416) * t414;
t442 = t417 * t444;
t510 = (t411 * t443 - t415 * t442) * t416 + t444 * t483;
t383 = t411 * t438 + t415 * t444;
t404 = t443 * qJD(2);
t405 = t444 * qJD(2);
t480 = t415 * t417;
t508 = (t404 * t411 - t405 * t480) * t416 - t383 * qJD(4);
t419 = sin(qJ(4));
t421 = cos(qJ(4));
t482 = t412 * t421;
t469 = t413 * t482;
t484 = t411 * t417;
t505 = qJD(4) * t511 - t404 * t415 - t405 * t484;
t337 = -t405 * t469 + t419 * t505 - t421 * t508;
t361 = t383 * t419 - t421 * t511;
t359 = t361 ^ 2;
t461 = t499 * t413;
t468 = t414 * t501;
t399 = t414 * t464 + (t417 * t468 + t461) * t411;
t465 = t500 * t411;
t466 = t501 * t415;
t447 = t417 * t466 - t465;
t458 = t413 * t468;
t453 = (t414 * t447 + t415 * t461) * t416 + (t417 * t499 - t458) * t412;
t502 = -t399 * t419 + t421 * t453;
t372 = 0.1e1 / t502 ^ 2;
t353 = t359 * t372 + 0.1e1;
t488 = t361 * t372;
t375 = t399 * t421 + t419 * t453;
t403 = (-t417 * t465 + t466) * t414;
t401 = qJD(2) * t403;
t446 = t507 * qJD(2);
t357 = qJD(4) * t375 + t401 * t419 - t421 * t446;
t371 = 0.1e1 / t502;
t489 = t357 * t371 * t372;
t506 = -0.2e1 * (t337 * t488 + t359 * t489) / t353 ^ 2;
t406 = t445 * qJD(2);
t407 = t409 * qJD(2);
t504 = -qJD(4) * t509 - t406 * t415 - t407 * t484;
t395 = -t409 * t484 - t415 * t445;
t394 = -t409 * t480 + t411 * t445;
t449 = t394 * t416 + t409 * t483;
t503 = -t395 * t419 + t421 * t449;
t354 = atan2(-t361, -t502);
t349 = sin(t354);
t350 = cos(t354);
t331 = -t349 * t361 - t350 * t502;
t328 = 0.1e1 / t331;
t384 = t409 * t415 + t411 * t439;
t365 = t384 * t421 - t419 * t509;
t376 = t412 * t433 + t416 * t440;
t418 = sin(qJ(5));
t420 = cos(qJ(5));
t348 = t365 * t420 + t376 * t418;
t344 = 0.1e1 / t348;
t329 = 0.1e1 / t331 ^ 2;
t345 = 0.1e1 / t348 ^ 2;
t351 = 0.1e1 / t353;
t321 = (t337 * t371 + t357 * t488) * t351;
t454 = t349 * t502 - t350 * t361;
t317 = t321 * t454 - t349 * t337 + t350 * t357;
t496 = t317 * t328 * t329;
t389 = t406 * t411 - t407 * t480;
t470 = t419 * t483;
t477 = qJD(4) * t419;
t340 = t389 * t416 * t419 - t384 * t477 + t407 * t470 + t421 * t504;
t481 = t413 * t416;
t377 = -t389 * t412 + t407 * t481;
t332 = qJD(5) * t348 + t340 * t418 - t377 * t420;
t347 = t365 * t418 - t376 * t420;
t343 = t347 ^ 2;
t336 = t343 * t345 + 0.1e1;
t492 = t345 * t347;
t475 = qJD(5) * t347;
t333 = t340 * t420 + t377 * t418 - t475;
t493 = t333 * t344 * t345;
t495 = (t332 * t492 - t343 * t493) / t336 ^ 2;
t364 = t384 * t419 + t421 * t509;
t494 = t329 * t364;
t491 = t349 * t364;
t490 = t350 * t364;
t479 = t416 * t421;
t476 = qJD(4) * t421;
t360 = t364 ^ 2;
t327 = t329 * t360 + 0.1e1;
t339 = t384 * t476 - t389 * t479 - t407 * t469 + t419 * t504;
t473 = 0.2e1 * (t339 * t494 - t360 * t496) / t327 ^ 2;
t472 = -0.2e1 * t495;
t471 = t347 * t493;
t460 = 0.2e1 * t364 * t496;
t459 = 0.2e1 * t361 * t489;
t369 = t395 * t421 + t419 * t449;
t379 = -t394 * t412 + t409 * t481;
t356 = t369 * t420 + t379 * t418;
t355 = t369 * t418 - t379 * t420;
t452 = -t344 * t418 + t420 * t492;
t363 = t383 * t421 + t419 * t511;
t451 = t363 * t371 + t375 * t488;
t393 = -t411 * t442 - t415 * t443;
t367 = t393 * t419 - t421 * t510;
t380 = t403 * t419 - t421 * t507;
t450 = t367 * t371 + t380 * t488;
t392 = t406 * t484 - t407 * t415;
t391 = t406 * t480 + t407 * t411;
t378 = -t391 * t412 - t406 * t481;
t366 = t403 * t476 + t507 * t477 + (-t458 * t482 + ((-t417 * t467 - t464) * t419 + t447 * t479) * t414) * qJD(2);
t358 = qJD(4) * t502 + t401 * t421 + t446 * t419;
t342 = t392 * t421 + (t391 * t416 - t406 * t483) * t419 + t503 * qJD(4);
t341 = (t404 * t484 - t405 * t415) * t419 + (-(t404 * t480 + t405 * t411) * t416 + t404 * t483) * t421 + (t393 * t421 + t419 * t510) * qJD(4);
t338 = t405 * t470 + t419 * t508 + t421 * t505;
t334 = 0.1e1 / t336;
t325 = 0.1e1 / t327;
t323 = t450 * t351;
t322 = t451 * t351;
t319 = t323 * t454 - t349 * t367 + t350 * t380;
t318 = t322 * t454 - t349 * t363 + t350 * t375;
t316 = t450 * t506 + (t380 * t459 + t341 * t371 + (t337 * t380 + t357 * t367 + t361 * t366) * t372) * t351;
t314 = t451 * t506 + (t375 * t459 + t338 * t371 + (t337 * t375 + t357 * t363 + t358 * t361) * t372) * t351;
t1 = [0, t316, 0, t314, 0, 0; 0 (t319 * t494 + t328 * t503) * t473 + (t319 * t460 + (t503 * t317 - t319 * t339 - (-t316 * t361 - t323 * t337 + t366 + (t323 * t502 - t367) * t321) * t490 - (t316 * t502 - t323 * t357 - t341 + (t323 * t361 - t380) * t321) * t491) * t329 + (qJD(4) * t369 - t391 * t479 + t392 * t419 + t406 * t469) * t328) * t325, 0 (t318 * t494 - t328 * t365) * t473 + (t318 * t460 + t340 * t328 + (-t365 * t317 - t318 * t339 - (-t314 * t361 - t322 * t337 + t358 + (t322 * t502 - t363) * t321) * t490 - (t314 * t502 - t322 * t357 - t338 + (t322 * t361 - t375) * t321) * t491) * t329) * t325, 0, 0; 0, 0.2e1 * (-t344 * t355 + t356 * t492) * t495 + ((qJD(5) * t356 + t342 * t418 - t378 * t420) * t344 + 0.2e1 * t356 * t471 + (-t355 * t333 - (-qJD(5) * t355 + t342 * t420 + t378 * t418) * t347 - t356 * t332) * t345) * t334, 0, t452 * t364 * t472 + (t452 * t339 + ((-qJD(5) * t344 - 0.2e1 * t471) * t420 + (t332 * t420 + (t333 - t475) * t418) * t345) * t364) * t334, t472 + 0.2e1 * (t332 * t345 * t334 + (-t334 * t493 - t345 * t495) * t347) * t347, 0;];
JaD_rot  = t1;
