% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR9
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
% Datum: 2019-02-26 21:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_jacobiaD_rot_6_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:05:30
% EndTime: 2019-02-26 21:05:33
% DurationCPUTime: 3.54s
% Computational Cost: add. (17692->166), mult. (41237->320), div. (726->12), fcn. (53008->17), ass. (0->146)
t493 = sin(pkin(12));
t498 = cos(pkin(6));
t453 = t498 * t493;
t496 = cos(pkin(12));
t499 = sin(qJ(1));
t500 = cos(qJ(1));
t392 = -t499 * t453 + t500 * t496;
t389 = t392 * qJD(1);
t401 = sin(qJ(3));
t403 = cos(qJ(3));
t391 = t500 * t453 + t499 * t496;
t494 = sin(pkin(7));
t495 = sin(pkin(6));
t451 = t495 * t494;
t439 = t500 * t451;
t497 = cos(pkin(7));
t454 = t498 * t496;
t507 = -t500 * t454 + t499 * t493;
t508 = t507 * t497;
t413 = t508 + t439;
t503 = t391 * t401 + t413 * t403;
t422 = t499 * t454 + t500 * t493;
t388 = t422 * qJD(1);
t436 = t499 * t451;
t513 = qJD(1) * t436 - t388 * t497;
t351 = t503 * qJD(3) - t389 * t403 - t513 * t401;
t470 = qJ(4) + pkin(13);
t399 = sin(t470);
t452 = t497 * t495;
t437 = t499 * t452;
t423 = qJD(1) * t437 + t388 * t494;
t459 = cos(t470);
t474 = t391 * t403;
t372 = t413 * t401 - t474;
t414 = -t500 * t452 + t494 * t507;
t522 = -t372 * t399 - t414 * t459;
t329 = qJD(4) * t522 + t351 * t459 - t423 * t399;
t360 = t372 * t459 - t414 * t399;
t524 = t360 * qJD(4) + t351 * t399 + t423 * t459;
t354 = t522 ^ 2;
t419 = t496 * t452 + t494 * t498;
t450 = t495 * t493;
t383 = t419 * t401 + t403 * t450;
t390 = -t496 * t451 + t498 * t497;
t429 = -t383 * t399 + t390 * t459;
t364 = 0.1e1 / t429 ^ 2;
t338 = t354 * t364 + 0.1e1;
t332 = 0.1e1 / t338;
t367 = t383 * t459 + t390 * t399;
t382 = -t401 * t450 + t419 * t403;
t377 = t382 * qJD(3);
t352 = t367 * qJD(4) + t377 * t399;
t363 = 0.1e1 / t429;
t479 = t522 * t364;
t445 = t352 * t479 - t363 * t524;
t309 = t445 * t332;
t339 = atan2(-t522, -t429);
t330 = sin(t339);
t331 = cos(t339);
t447 = t330 * t429 - t331 * t522;
t304 = t447 * t309 + t330 * t524 + t331 * t352;
t321 = -t330 * t522 - t331 * t429;
t319 = 0.1e1 / t321 ^ 2;
t521 = t304 * t319;
t518 = -t422 * t497 + t436;
t374 = t392 * t403 + t518 * t401;
t412 = t422 * t494 + t437;
t362 = t374 * t459 + t412 * t399;
t505 = t518 * t403;
t373 = t392 * t401 - t505;
t400 = sin(qJ(6));
t402 = cos(qJ(6));
t342 = t362 * t400 - t373 * t402;
t517 = 0.2e1 * t342;
t318 = 0.1e1 / t321;
t516 = t318 * t521;
t387 = t391 * qJD(1);
t473 = qJD(3) * t401;
t506 = t413 * qJD(1);
t347 = t505 * qJD(3) - t387 * t403 - t392 * t473 + t506 * t401;
t411 = qJD(1) * t414;
t325 = t362 * qJD(4) + t347 * t399 + t459 * t411;
t410 = t412 * t459;
t361 = t374 * t399 - t410;
t501 = 0.2e1 * t361;
t458 = t501 * t516;
t512 = -t319 * t325 + t458;
t511 = (t401 * t508 - t474) * qJD(3) - t389 * t401 + t513 * t403 + t439 * t473;
t355 = t361 ^ 2;
t315 = t319 * t355 + 0.1e1;
t487 = t319 * t361;
t469 = 0.2e1 * (t325 * t487 - t355 * t516) / t315 ^ 2;
t510 = t352 * t364;
t442 = -t363 * t503 + t382 * t479;
t509 = t399 * t442;
t448 = qJD(4) * t459;
t343 = t362 * t402 + t373 * t400;
t335 = 0.1e1 / t343;
t336 = 0.1e1 / t343 ^ 2;
t502 = -0.2e1 * t522;
t481 = t363 * t510;
t490 = (t354 * t481 - t479 * t524) / t338 ^ 2;
t472 = qJD(4) * t399;
t326 = qJD(4) * t410 + t347 * t459 - t374 * t472 - t399 * t411;
t346 = t374 * qJD(3) - t387 * t401 - t506 * t403;
t471 = qJD(6) * t342;
t317 = t326 * t402 + t346 * t400 - t471;
t489 = t317 * t335 * t336;
t334 = t342 ^ 2;
t324 = t334 * t336 + 0.1e1;
t322 = 0.1e1 / t324;
t486 = t322 * t336;
t316 = t343 * qJD(6) + t326 * t400 - t346 * t402;
t482 = t336 * t342;
t485 = 0.1e1 / t324 ^ 2 * (t316 * t482 - t334 * t489);
t484 = t330 * t361;
t483 = t331 * t361;
t480 = t522 * t363;
t477 = t373 * t399;
t468 = -0.2e1 * t490;
t466 = -0.2e1 * t485;
t465 = t336 * t485;
t464 = t363 * t490;
t463 = t316 * t486;
t462 = t342 * t489;
t457 = 0.2e1 * t462;
t456 = t481 * t502;
t449 = t373 * t459;
t341 = t360 * t402 - t400 * t503;
t340 = t360 * t400 + t402 * t503;
t444 = -t335 * t400 + t402 * t482;
t443 = -t360 * t363 + t367 * t479;
t435 = -t330 + (-t331 * t480 + t330) * t332;
t434 = qJD(6) * t449 + t347;
t424 = qJD(6) * t374 - t459 * t346 + t373 * t472;
t378 = t383 * qJD(3);
t353 = t429 * qJD(4) + t377 * t459;
t345 = t374 * t400 - t402 * t449;
t313 = 0.1e1 / t315;
t312 = t332 * t509;
t310 = t443 * t332;
t306 = (t330 * t503 + t331 * t382) * t399 + t447 * t312;
t305 = t447 * t310 + t330 * t360 + t331 * t367;
t302 = t443 * t468 + (-t367 * t456 - t329 * t363 + (-t352 * t360 + t353 * t522 - t367 * t524) * t364) * t332;
t301 = t468 * t509 + ((-t382 * t456 + t511 * t363 + (-t352 * t503 - t378 * t522 - t382 * t524) * t364) * t399 + t442 * t448) * t332;
t1 = [-t464 * t501 + (t325 * t363 + t361 * t510) * t332, 0, t301, t302, 0, 0; t522 * t318 * t469 + (t524 * t318 + t522 * t521 - (t435 * t325 + ((t309 * t332 * t480 + t468) * t330 + (-t464 * t502 - t309 + (t309 - t445) * t332) * t331) * t361) * t487) * t313 + (t512 * t313 + t487 * t469) * t435 * t361, 0 (t306 * t487 + t318 * t477) * t469 + ((-t346 * t399 - t373 * t448) * t318 + t512 * t306 + (t477 * t304 - (t382 * t448 - t301 * t522 + t312 * t524 - t378 * t399 + (t312 * t429 + t399 * t503) * t309) * t483 - (t503 * t448 + t301 * t429 - t312 * t352 - t511 * t399 + (t312 * t522 - t382 * t399) * t309) * t484) * t319) * t313 (t305 * t487 - t318 * t362) * t469 + (t305 * t458 + t326 * t318 + (-t362 * t304 - t305 * t325 - (-t302 * t522 + t310 * t524 + t353 + (t310 * t429 + t360) * t309) * t483 - (t302 * t429 - t310 * t352 + t329 + (t310 * t522 - t367) * t309) * t484) * t319) * t313, 0, 0; 0.2e1 * (-t335 * t340 + t341 * t482) * t485 + ((t341 * qJD(6) + t329 * t400 - t402 * t511) * t335 + t341 * t457 + (-t340 * t317 - (-t340 * qJD(6) + t329 * t402 + t400 * t511) * t342 - t341 * t316) * t336) * t322, 0 (t465 * t517 - t463) * t345 + (-t317 * t486 + t335 * t466) * (-t374 * t402 - t400 * t449) + ((t424 * t400 - t434 * t402) * t335 - (t434 * t400 + t424 * t402) * t482 + t345 * t457) * t322, t444 * t361 * t466 + (t444 * t325 + ((-qJD(6) * t335 - 0.2e1 * t462) * t402 + (t316 * t402 + (t317 - t471) * t400) * t336) * t361) * t322, 0, t466 + (t463 + (-t322 * t489 - t465) * t342) * t517;];
JaD_rot  = t1;
