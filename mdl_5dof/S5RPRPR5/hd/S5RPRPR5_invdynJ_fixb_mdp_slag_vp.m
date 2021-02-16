% Calculate vector of inverse dynamics joint torques for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% MDP [24x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S5RPRPR5_convert_par2_MPV_fixb.m
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:56
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR5_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(24,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [24 1]), ...
  'S5RPRPR5_invdynJ_fixb_mdp_slag_vp: MDP has to be [24x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:55:25
% EndTime: 2021-01-15 11:55:37
% DurationCPUTime: 4.42s
% Computational Cost: add. (2705->395), mult. (6665->533), div. (0->0), fcn. (4900->14), ass. (0->202)
t430 = sin(pkin(8));
t424 = t430 ^ 2;
t539 = 0.2e1 * t424;
t438 = cos(qJ(3));
t497 = qJD(3) * t438;
t551 = qJ(2) * t497;
t429 = sin(pkin(9));
t431 = cos(pkin(9));
t435 = sin(qJ(3));
t383 = t429 * t438 + t431 * t435;
t450 = qJD(1) * t383;
t351 = t430 * t450;
t437 = cos(qJ(5));
t514 = t437 * t351;
t502 = qJD(1) * t430;
t480 = t435 * t502;
t464 = t429 * t480;
t500 = qJD(1) * t438;
t479 = t430 * t500;
t355 = t431 * t479 - t464;
t434 = sin(qJ(5));
t526 = t355 * t434;
t310 = t514 + t526;
t432 = cos(pkin(8));
t501 = qJD(1) * t432;
t407 = -qJD(3) + t501;
t401 = -qJD(5) + t407;
t527 = t310 * t401;
t529 = qJDD(1) * pkin(1);
t436 = sin(qJ(1));
t439 = cos(qJ(1));
t544 = -g(2) * t439 - g(3) * t436;
t454 = -qJDD(2) + t529 + t544;
t488 = qJDD(1) * t435;
t475 = t430 * t488;
t490 = qJD(1) * qJD(3);
t477 = t438 * t490;
t550 = t430 * t477 + t475;
t491 = qJD(1) * qJD(2);
t492 = qJ(2) * qJDD(1);
t453 = t491 + t492;
t390 = -pkin(2) * t432 - pkin(6) * t430 - pkin(1);
t369 = t390 * qJDD(1) + qJDD(2);
t361 = t438 * t369;
t489 = qJDD(1) * t432;
t405 = -qJDD(3) + t489;
t530 = qJ(2) * t438;
t406 = t432 * t530;
t496 = qJD(4) * t430;
t499 = qJD(2) * t432;
t447 = -t435 * t499 - t438 * t496;
t522 = t430 * t438;
t484 = qJ(4) * t522;
t531 = qJ(2) * t435;
t486 = t432 * t531;
t448 = -t484 - t486;
t523 = t430 * t435;
t485 = qJ(4) * t523;
t370 = t390 * qJD(1) + qJD(2);
t498 = qJD(3) * t370;
t288 = -t435 * t498 - pkin(3) * t405 + t361 + t448 * qJDD(1) + ((-t406 + t485) * qJD(3) + t447) * qJD(1);
t442 = t448 * qJD(3) - t435 * t496;
t478 = t438 * t491;
t467 = qJDD(1) * t406 + t435 * t369 + t370 * t497 + t432 * t478;
t294 = -qJ(4) * t475 + t442 * qJD(1) + t467;
t277 = t431 * t288 - t294 * t429;
t487 = qJDD(1) * t438;
t474 = t430 * t487;
t389 = t431 * t474;
t449 = t383 * qJD(3);
t476 = t429 * t488;
t323 = -t389 + (qJD(1) * t449 + t476) * t430;
t275 = -pkin(4) * t405 + pkin(7) * t323 + t277;
t362 = t438 * t370;
t331 = t448 * qJD(1) + t362;
t321 = -pkin(3) * t407 + t331;
t332 = -qJ(4) * t480 + qJD(1) * t406 + t370 * t435;
t521 = t431 * t332;
t299 = t429 * t321 + t521;
t537 = pkin(7) * t351;
t285 = t299 - t537;
t495 = qJD(5) * t434;
t549 = t434 * t275 - t285 * t495;
t278 = t429 * t288 + t431 * t294;
t482 = t429 * t474 + t550 * t431;
t322 = qJD(3) * t464 - t482;
t276 = pkin(7) * t322 + t278;
t371 = pkin(3) * t480 + qJ(2) * t502 + qJD(4);
t330 = pkin(4) * t351 + t371;
t426 = qJ(3) + pkin(9);
t421 = qJ(5) + t426;
t413 = cos(t421);
t412 = sin(t421);
t512 = t439 * t412;
t520 = t432 * t436;
t347 = t413 * t520 - t512;
t519 = t432 * t439;
t349 = t412 * t436 + t413 * t519;
t535 = g(1) * t430;
t548 = g(2) * t347 - g(3) * t349 - t437 * t276 + t330 * t310 + t413 * t535 - t549;
t396 = -qJDD(5) + t405;
t386 = t396 * MDP(22);
t457 = -t351 * t434 + t437 * t355;
t547 = t310 * MDP(18) * t457 + (-t310 ^ 2 + t457 ^ 2) * MDP(19) - t386;
t545 = t453 * t432;
t528 = t457 * t401;
t513 = t438 * t439;
t518 = t435 * t436;
t374 = -t432 * t518 - t513;
t516 = t436 * t438;
t517 = t435 * t439;
t376 = t432 * t517 - t516;
t543 = -g(2) * t374 - g(3) * t376;
t524 = t429 * t435;
t456 = -t431 * t438 + t524;
t542 = t456 * qJD(3);
t346 = -t412 * t520 - t413 * t439;
t348 = -t413 * t436 + t432 * t512;
t471 = t437 * t275 - t434 * t276;
t541 = -g(2) * t346 - g(3) * t348 - t330 * t457 + t412 * t535 + t471;
t540 = (pkin(3) * t438 + pkin(2)) * t432 + t430 * (qJ(4) + pkin(6)) + pkin(1);
t470 = t437 * t322 + t323 * t434;
t280 = t457 * qJD(5) - t470;
t425 = t432 ^ 2;
t538 = pkin(3) * t429;
t536 = pkin(7) * t355;
t532 = MDP(9) * t430;
t440 = qJD(1) ^ 2;
t525 = t424 * t440;
t326 = t429 * t332;
t298 = t431 * t321 - t326;
t283 = -pkin(4) * t407 + t298 - t536;
t515 = t437 * t283;
t509 = t390 * t497 + t438 * t499;
t319 = t442 + t509;
t320 = (-t406 + (qJ(4) * t430 - t390) * t435) * qJD(3) + t447;
t292 = t431 * t319 + t429 * t320;
t302 = t431 * t331 - t326;
t381 = t438 * t390;
t336 = -t484 + t381 + (-pkin(3) - t531) * t432;
t508 = t435 * t390 + t406;
t341 = -t485 + t508;
t305 = t429 * t336 + t431 * t341;
t511 = -t432 * t450 + t449;
t510 = -t456 * t501 + t542;
t379 = (pkin(3) * t497 + qJD(2)) * t430;
t384 = pkin(3) * t523 + t430 * qJ(2);
t506 = t424 + t425;
t428 = t438 ^ 2;
t505 = t435 ^ 2 - t428;
t504 = MDP(10) * t430;
t494 = t405 * MDP(11);
t493 = qJD(3) + t407;
t483 = -qJD(5) * t514 + t434 * t322 - t437 * t323;
t472 = t506 * t440;
t291 = -t319 * t429 + t431 * t320;
t301 = -t331 * t429 - t521;
t304 = t431 * t336 - t341 * t429;
t469 = qJD(1) * t493;
t468 = t405 + t489;
t466 = 0.2e1 * t506;
t465 = qJD(3) * t486;
t461 = g(2) * t436 - g(3) * t439;
t460 = qJD(5) * t383 + t511;
t459 = -qJD(5) * t456 - t510;
t458 = -t434 * t283 - t437 * t285;
t367 = t383 * t430;
t368 = t456 * t430;
t324 = t437 * t367 - t368 * t434;
t325 = -t367 * t434 - t368 * t437;
t337 = t550 * pkin(3) + t453 * t430 + qJDD(4);
t455 = qJD(3) * (t407 + t501);
t415 = pkin(3) * t431 + pkin(4);
t452 = t415 * t434 + t437 * t538;
t451 = t415 * t437 - t434 * t538;
t279 = -t355 * t495 + t483;
t446 = -t407 ^ 2 - t525;
t444 = t466 * t491 - t461;
t419 = cos(t426);
t418 = sin(t426);
t414 = -pkin(3) * t435 - qJ(2);
t377 = t432 * t513 + t518;
t375 = t432 * t516 - t517;
t366 = t418 * t436 + t419 * t519;
t365 = t418 * t519 - t419 * t436;
t364 = -t418 * t439 + t419 * t520;
t363 = -t418 * t520 - t419 * t439;
t358 = t430 * t449;
t354 = t430 * t542;
t340 = pkin(3) * t479 + pkin(4) * t355;
t338 = pkin(4) * t367 + t384;
t333 = -pkin(4) * t354 + t379;
t303 = -pkin(4) * t322 + t337;
t300 = -pkin(7) * t367 + t305;
t297 = -pkin(4) * t432 + pkin(7) * t368 + t304;
t296 = t325 * qJD(5) - t437 * t354 - t358 * t434;
t295 = -t324 * qJD(5) + t354 * t434 - t358 * t437;
t290 = t302 - t536;
t289 = t301 + t537;
t282 = pkin(7) * t354 + t292;
t281 = pkin(7) * t358 + t291;
t1 = [qJDD(1) * MDP(1) + t544 * MDP(2) + t461 * MDP(3) + (t466 * t492 + t444) * MDP(5) + (t454 * pkin(1) + (t492 * t506 + t444) * qJ(2)) * MDP(6) + (qJDD(1) * t428 - 0.2e1 * t435 * t477) * t424 * MDP(7) + (-t435 * t487 + t505 * t490) * MDP(8) * t539 + (t435 * t455 - t468 * t438) * t532 + (t468 * t435 + t438 * t455) * t504 + (-g(2) * t377 - g(3) * t375 - t381 * t405 + (t539 + t425) * qJD(1) * t551 + (qJD(3) * t390 * t407 + t453 * t539) * t435) * MDP(12) + ((-t465 + t509) * t407 + t508 * t405 + g(2) * t376 - g(3) * t374 + (t478 + (-t435 * t490 + t487) * qJ(2)) * t539) * MDP(13) + (-g(2) * t366 - g(3) * t364 - t291 * t407 - t304 * t405 - t322 * t384 + t337 * t367 + t351 * t379 - t354 * t371) * MDP(14) + (g(2) * t365 - g(3) * t363 + t292 * t407 + t305 * t405 - t323 * t384 - t337 * t368 + t355 * t379 - t358 * t371) * MDP(15) + (t277 * t368 - t278 * t367 - t291 * t355 - t292 * t351 + t298 * t358 + t299 * t354 + t304 * t323 + t305 * t322 + t430 * t544) * MDP(16) + (t278 * t305 + t299 * t292 + t277 * t304 + t298 * t291 + t337 * t384 + t371 * t379 - g(2) * (-t414 * t436 + t439 * t540) - g(3) * (t414 * t439 + t436 * t540)) * MDP(17) + (t279 * t325 + t295 * t457) * MDP(18) + (-t279 * t324 - t280 * t325 - t295 * t310 - t296 * t457) * MDP(19) + (-t295 * t401 - t325 * t396) * MDP(20) + (t296 * t401 + t324 * t396) * MDP(21) + (-(t297 * t437 - t300 * t434) * t396 + t333 * t310 + t338 * t280 + t303 * t324 + t330 * t296 - g(2) * t349 - g(3) * t347 + (-t281 * t437 + t282 * t434 - (-t297 * t434 - t300 * t437) * qJD(5)) * t401) * MDP(23) + (g(2) * t348 - g(3) * t346 + t338 * t279 + t330 * t295 + t303 * t325 + t333 * t457 + ((-qJD(5) * t300 + t281) * t401 + t297 * t396) * t434 + ((qJD(5) * t297 + t282) * t401 + t300 * t396) * t437) * MDP(24) + ((t454 + t529) * MDP(4) + t494 + (-t361 + t407 * t551 + (qJ(2) * t405 + qJD(2) * t407 + t498 + t545) * t435) * MDP(12) + (-qJD(1) * t465 + t467) * MDP(13) - t277 * MDP(14) + t278 * MDP(15) - t279 * MDP(20) + t280 * MDP(21) + t386 + (-t458 * qJD(5) - t471) * MDP(23) + ((qJD(5) * t283 + t276) * t437 + t549) * MDP(24)) * t432; -MDP(4) * t489 - MDP(5) * t472 + (-qJ(2) * t472 - t454) * MDP(6) + (-t405 * t438 + t446 * t435) * MDP(12) + (t405 * t435 + t446 * t438) * MDP(13) + (-t351 * t502 + t405 * t456 + t511 * t407) * MDP(14) + (-t355 * t502 + t383 * t405 - t510 * t407) * MDP(15) + (t322 * t383 - t323 * t456 + t510 * t351 + t511 * t355) * MDP(16) + (-t277 * t456 + t278 * t383 - t511 * t298 - t510 * t299 - t371 * t502 - t544) * MDP(17) + (-(-t383 * t434 - t437 * t456) * t396 - t310 * t502 + (t459 * t434 + t460 * t437) * t401) * MDP(23) + ((t383 * t437 - t434 * t456) * t396 - t457 * t502 + (-t460 * t434 + t459 * t437) * t401) * MDP(24); t438 * t435 * MDP(7) * t525 - t505 * MDP(8) * t525 + (-t435 * t469 + t487) * t532 + (-t438 * t469 - t488) * t504 - t494 + (t361 + (-t432 * t469 - t525) * t530 + (-t493 * t370 + t535 - t545) * t435 + t543) * MDP(12) + (g(1) * t522 + g(2) * t375 - g(3) * t377 - t362 * t407 + (t493 * t501 + t525) * t531 - t467) * MDP(13) + (t418 * t535 - g(2) * t363 - g(3) * t365 + t301 * t407 - t355 * t371 + (-t351 * t479 - t405 * t431) * pkin(3) + t277) * MDP(14) + (t419 * t535 + g(2) * t364 - g(3) * t366 - t302 * t407 + t351 * t371 + (-t355 * t479 + t405 * t429) * pkin(3) - t278) * MDP(15) + ((t299 + t301) * t355 + (-t298 + t302) * t351 + (t322 * t429 + t323 * t431) * pkin(3)) * MDP(16) + (-t298 * t301 - t299 * t302 + (t278 * t429 + t277 * t431 + (g(1) * t435 - t371 * t500) * t430 + t543) * pkin(3)) * MDP(17) + (t279 - t527) * MDP(20) + (-t280 - t528) * MDP(21) + (-t451 * t396 + (t289 * t437 - t290 * t434) * t401 - t340 * t310 + (t452 * t401 + t458) * qJD(5) + t541) * MDP(23) + (t452 * t396 - (t289 * t434 + t290 * t437) * t401 - t340 * t457 + (t451 * t401 - t515) * qJD(5) + t548) * MDP(24) + t547; (-t355 * t407 + t482) * MDP(14) + (t351 * t407 + t389) * MDP(15) + (-t351 ^ 2 - t355 ^ 2) * MDP(16) + (g(1) * t432 + t298 * t355 + t299 * t351 + t337) * MDP(17) + (t280 - t528) * MDP(23) + (t279 + t527) * MDP(24) + (-MDP(15) * t476 - t461 * MDP(17) + (-MDP(14) * t524 - t383 * MDP(15)) * t490) * t430; (t483 - t527) * MDP(20) + (t470 - t528) * MDP(21) + (t458 * t401 + t541) * MDP(23) + (-(-t285 * t434 + t515) * t401 + t548) * MDP(24) + (-MDP(20) * t526 - t457 * MDP(21) + t458 * MDP(23) - MDP(24) * t515) * qJD(5) + t547;];
tau = t1;
