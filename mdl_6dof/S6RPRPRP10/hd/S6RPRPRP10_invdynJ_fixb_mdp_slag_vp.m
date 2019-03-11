% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% MDP [28x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP10_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP10_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(28,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP10_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP10_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_invdynJ_fixb_mdp_slag_vp: pkin has to be [8x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [28 1]), ...
  'S6RPRPRP10_invdynJ_fixb_mdp_slag_vp: MDP has to be [28x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:32:43
% EndTime: 2019-03-09 03:32:50
% DurationCPUTime: 5.74s
% Computational Cost: add. (3049->504), mult. (5364->590), div. (0->0), fcn. (2904->6), ass. (0->203)
t405 = sin(qJ(5));
t408 = cos(qJ(5));
t406 = sin(qJ(3));
t493 = qJD(1) * t406;
t466 = t408 * t493;
t485 = qJD(5) * t405;
t409 = cos(qJ(3));
t479 = qJD(1) * qJD(3);
t462 = t409 * t479;
t478 = qJDD(1) * t406;
t543 = t462 + t478;
t302 = qJD(3) * t485 - qJD(5) * t466 - t408 * qJDD(3) - t405 * t543;
t491 = qJD(3) * t405;
t344 = -t466 + t491;
t492 = qJD(1) * t409;
t373 = qJD(5) + t492;
t516 = t344 * t373;
t544 = t302 - t516;
t545 = t544 * MDP(26);
t542 = t373 ^ 2;
t412 = -pkin(1) - pkin(7);
t525 = pkin(4) - t412;
t410 = cos(qJ(1));
t395 = g(2) * t410;
t407 = sin(qJ(1));
t396 = g(1) * t407;
t540 = t396 - t395;
t541 = t540 * t406;
t389 = t409 * qJ(4);
t499 = t406 * pkin(3) - t389;
t527 = pkin(8) * t406;
t454 = -t499 - t527;
t339 = qJ(2) - t454;
t357 = t525 * t409;
t501 = t408 * t339 + t405 * t357;
t489 = qJD(3) * t408;
t346 = t405 * t493 + t489;
t486 = qJD(5) * t346;
t303 = qJDD(3) * t405 - t408 * t543 + t486;
t387 = t409 * qJDD(1);
t463 = t406 * t479;
t539 = -t463 + t387;
t343 = -qJDD(5) - t539;
t328 = t408 * t343;
t538 = -t373 * t485 - t328;
t398 = qJDD(1) * qJ(2);
t400 = qJD(1) * qJD(2);
t471 = 0.2e1 * t400;
t537 = 0.2e1 * t398 + t471;
t365 = qJD(1) * t412 + qJD(2);
t353 = t409 * t365;
t536 = qJD(4) - t353;
t535 = qJDD(1) * t412;
t474 = MDP(23) + MDP(25);
t473 = MDP(24) - MDP(27);
t352 = t406 * t365;
t324 = -pkin(4) * t493 + t352;
t401 = qJD(3) * qJ(4);
t317 = t324 + t401;
t291 = pkin(5) * t344 - qJ(6) * t346 + t317;
t331 = -t352 - t401;
t534 = -t331 * MDP(17) + t291 * MDP(28) + t346 * t473;
t506 = t409 * t410;
t334 = t405 * t506 + t407 * t408;
t508 = t407 * t409;
t336 = -t405 * t508 + t408 * t410;
t394 = g(3) * t406;
t441 = pkin(3) * t543 + qJ(4) * t463 + t398 + t400;
t459 = qJD(3) * pkin(8) - qJD(4);
t480 = qJ(4) * qJDD(1);
t290 = pkin(8) * t478 + (qJD(1) * t459 - t480) * t409 + t441;
t411 = -pkin(3) - pkin(8);
t490 = qJD(3) * t406;
t348 = t365 * t490;
t362 = qJDD(2) + t535;
t457 = -t362 * t409 + qJDD(4);
t442 = t348 + t457;
t296 = pkin(4) * t539 + t411 * qJDD(3) + t442;
t482 = pkin(4) * t492 + t536;
t312 = qJD(3) * t411 + t482;
t484 = qJD(5) * t408;
t469 = t408 * t290 + t405 * t296 + t312 * t484;
t444 = -t389 + t527;
t500 = pkin(3) * t493 + qJD(1) * qJ(2);
t315 = qJD(1) * t444 + t500;
t487 = qJD(5) * t315;
t533 = -g(1) * t336 - g(2) * t334 - (t487 + t394) * t405 + t469;
t351 = t406 * t362;
t397 = qJDD(3) * qJ(4);
t399 = qJD(3) * qJD(4);
t488 = qJD(3) * t409;
t305 = -t365 * t488 - t351 - t397 - t399;
t521 = qJDD(3) * pkin(3);
t306 = t442 - t521;
t460 = qJD(3) * pkin(3) - qJD(4);
t329 = -t353 - t460;
t532 = (t329 * t406 - t331 * t409) * qJD(3) - t305 * t406 - t306 * t409;
t330 = -qJ(4) * t492 + t500;
t355 = qJ(2) + t499;
t475 = qJDD(3) * t412;
t531 = (qJD(1) * t355 + t330) * qJD(3) + t475;
t530 = t346 ^ 2;
t528 = pkin(5) * t343;
t526 = g(3) * t409;
t524 = pkin(1) * qJDD(1);
t414 = qJD(1) ^ 2;
t523 = qJ(2) * t414;
t522 = qJ(6) * t343;
t289 = t312 * t405 + t315 * t408;
t520 = t289 * t373;
t519 = t302 * t408;
t518 = t303 * t405;
t517 = t344 * t346;
t515 = t344 * t405;
t514 = t344 * t408;
t513 = t346 * t373;
t512 = t346 * t405;
t511 = t346 * t408;
t510 = t406 * t407;
t509 = t406 * t410;
t507 = t408 * t409;
t505 = t411 * t343;
t350 = pkin(3) * t492 + qJ(4) * t493;
t323 = pkin(8) * t492 + t350;
t503 = t408 * t323 + t405 * t324;
t446 = pkin(5) * t408 + qJ(6) * t405;
t432 = -pkin(4) - t446;
t425 = t432 * t409;
t502 = -qJD(1) * t425 + qJD(5) * t446 - qJD(6) * t408 + t536;
t498 = t410 * pkin(1) + t407 * qJ(2);
t403 = t406 ^ 2;
t404 = t409 ^ 2;
t496 = -t403 - t404;
t495 = t403 - t404;
t413 = qJD(3) ^ 2;
t494 = t413 + t414;
t483 = qJD(5) * t411;
t288 = t312 * t408 - t315 * t405;
t481 = qJD(6) - t288;
t477 = qJDD(3) * t406;
t476 = qJDD(3) * t409;
t472 = g(1) * t510;
t470 = t406 * t409 * t414;
t468 = -g(1) * t508 + g(2) * t506 + t394;
t465 = pkin(3) * t488 + qJ(4) * t490 + qJD(2);
t464 = g(1) * (pkin(3) * t508 + qJ(4) * t510);
t461 = qJDD(2) - t540;
t456 = t405 * t290 - t408 * t296 + t312 * t485 + t315 * t484;
t455 = pkin(3) * t510 + t410 * pkin(7) + t498;
t333 = t405 * t407 - t408 * t506;
t335 = t405 * t410 + t407 * t507;
t453 = g(1) * t333 - g(2) * t335;
t452 = g(1) * t334 - g(2) * t336;
t451 = g(1) * t410 + g(2) * t407;
t277 = qJD(6) * t373 - t315 * t485 + t469 - t522;
t282 = -pkin(5) * t373 + t481;
t449 = -t282 * t492 - t277;
t278 = qJDD(6) + t456 + t528;
t283 = qJ(6) * t373 + t289;
t448 = -t283 * t492 + t278;
t447 = -pkin(3) * t409 - qJ(4) * t406;
t445 = pkin(5) * t405 - qJ(6) * t408;
t439 = -t277 * t405 + t278 * t408;
t438 = t282 * t408 - t283 * t405;
t437 = t282 * t405 + t283 * t408;
t390 = t410 * qJ(2);
t431 = pkin(3) * t509 - qJ(4) * t506 + t390;
t354 = qJ(4) + t445;
t430 = -t373 * t483 - t526;
t429 = t405 * t343 - t373 * t484;
t428 = 0.2e1 * qJ(2) * t479 + t475;
t314 = t409 * t459 + t465;
t341 = t525 * t490;
t427 = t408 * t314 - t339 * t485 - t405 * t341 + t357 * t484;
t297 = -pkin(4) * t543 - t305;
t279 = pkin(5) * t303 + qJ(6) * t302 - qJD(6) * t346 + t297;
t426 = -t279 - t430;
t424 = -t412 * t413 - t451;
t422 = t330 * t492 + t457 - t468;
t421 = t291 * t373 - t505;
t420 = -t405 * t473 + t408 * t474;
t419 = g(1) * t335 + g(2) * t333 - t394 * t408 - t456;
t418 = t424 + t537;
t298 = (-qJD(1) * qJD(4) - t480) * t409 + t441;
t321 = -qJD(4) * t409 + t465;
t417 = -qJD(1) * t321 - qJDD(1) * t355 - t298 - t424;
t416 = t291 * t346 + qJDD(6) - t419;
t415 = (-t512 + t514) * MDP(26) - t437 * MDP(28) + (t405 * t474 + t408 * t473) * t373;
t382 = t406 * t412;
t371 = t412 * t488;
t358 = g(2) * t405 * t509;
t356 = -pkin(4) * t406 + t382;
t342 = -pkin(4) * t488 + t371;
t313 = t406 * t432 + t382;
t307 = pkin(5) * t346 + qJ(6) * t344;
t301 = -pkin(5) * t409 + t339 * t405 - t357 * t408;
t300 = qJ(6) * t409 + t501;
t295 = pkin(5) * t493 + t323 * t405 - t324 * t408;
t294 = -qJ(6) * t493 + t503;
t287 = t371 + (qJD(5) * t445 - qJD(6) * t405) * t406 + qJD(3) * t425;
t281 = pkin(5) * t490 + qJD(5) * t501 + t314 * t405 + t341 * t408;
t280 = -qJ(6) * t490 + qJD(6) * t409 + t427;
t1 = [(-t456 * t409 + t342 * t344 + t356 * t303 + ((-qJD(5) * t357 - t314) * t373 + t339 * t343) * t405 + ((-qJD(5) * t339 - t341) * t373 - t357 * t343 - t317 * t488) * t408 + (-qJD(3) * t288 - t297 * t408 + t317 * t485) * t406 + t452) * MDP(23) + (t496 * t535 - t532 + t540) * MDP(14) + t540 * MDP(2) + (-t451 + t537) * MDP(5) + ((t373 * t489 - t303) * t409 + (qJD(3) * t344 + t538) * t406) * MDP(21) + (t298 * t355 + t330 * t321 - g(1) * (t407 * t412 + t431) - g(2) * (-t389 * t407 + t455) + t532 * t412) * MDP(17) + (-t406 * t428 + t409 * t418) * MDP(13) + (t406 * t418 + t409 * t428) * MDP(12) + (t406 * t531 + t417 * t409) * MDP(16) + (t417 * t406 - t409 * t531) * MDP(15) + t451 * MDP(3) + qJDD(1) * MDP(1) + (t277 * t300 + t283 * t280 + t279 * t313 + t291 * t287 + t278 * t301 + t282 * t281 - g(1) * (-pkin(5) * t334 + pkin(8) * t509 - qJ(6) * t333 + t431) - g(2) * (pkin(4) * t410 + pkin(5) * t336 + qJ(6) * t335 + t455) + (g(1) * t525 - g(2) * t444) * t407) * MDP(28) + (qJDD(1) * t404 - 0.2e1 * t406 * t462) * MDP(7) + (-t406 * t413 + t476) * MDP(9) + (-t409 * t413 - t477) * MDP(10) + (-t280 * t344 + t281 * t346 - t300 * t303 - t301 * t302 + t437 * t488 + (qJD(5) * t438 + t277 * t408 + t278 * t405 - t451) * t406) * MDP(26) + (-t302 * t405 * t406 + (t405 * t488 + t406 * t484) * t346) * MDP(18) + (-t281 * t373 + t287 * t344 + t301 * t343 + t303 * t313 + (-t291 * t489 - t278) * t409 + (qJD(3) * t282 - t279 * t408 + t291 * t485) * t406 + t452) * MDP(25) + (-t343 * t409 - t373 * t490) * MDP(22) + (t280 * t373 - t287 * t346 - t300 * t343 + t302 * t313 + (-t291 * t491 + t277) * t409 + (-qJD(3) * t283 - t279 * t405 - t291 * t484) * t406 + t453) * MDP(27) + ((t373 * t491 - t302) * t409 + (-qJD(3) * t346 - t429) * t406) * MDP(20) + 0.2e1 * (-t387 * t406 + t479 * t495) * MDP(8) + (-t427 * t373 + t501 * t343 + t342 * t346 - t356 * t302 + ((qJD(3) * t317 + t487) * t405 - t469) * t409 + (qJD(3) * t289 + t297 * t405 + t317 * t484) * t406 - t453) * MDP(24) + ((t511 - t515) * t488 + (-t519 - t518 + (-t512 - t514) * qJD(5)) * t406) * MDP(19) + (-(qJDD(2) - t524) * pkin(1) - g(1) * (-pkin(1) * t407 + t390) - g(2) * t498 + (t471 + t398) * qJ(2)) * MDP(6) + (t461 - 0.2e1 * t524) * MDP(4); -t414 * MDP(5) + (t461 - t523) * MDP(6) + (MDP(12) - MDP(15)) * (-t406 * t494 + t476) + (-MDP(13) + MDP(16)) * (t409 * t494 + t477) + (MDP(14) * t496 - pkin(1) * MDP(6) + MDP(4)) * qJDD(1) + (-t330 * MDP(17) + t415) * qJD(1) + (-t305 * MDP(17) + t279 * MDP(28) - t473 * t302 + (t329 * MDP(17) + (-t511 - t515) * MDP(26) - t438 * MDP(28) + t420 * t373) * qJD(3)) * t406 + (-t306 * MDP(17) + (t518 - t519) * MDP(26) + t439 * MDP(28) + t534 * qJD(3) + t420 * t343 + t415 * qJD(5)) * t409 + t474 * (t406 * t303 + t344 * t488) + (-MDP(17) - MDP(28)) * t540; MDP(7) * t470 - t495 * MDP(8) * t414 + MDP(9) * t387 - MDP(10) * t478 + qJDD(3) * MDP(11) + ((t362 - t523) * t409 + t468) * MDP(12) + (t526 - t351 + (t540 + t523) * t406) * MDP(13) + (t447 * qJDD(1) + ((-t331 - t401) * t409 + (t329 + t460) * t406) * qJD(1)) * MDP(14) + (t350 * t493 + t422 - 0.2e1 * t521) * MDP(15) + (t351 + 0.2e1 * t397 + 0.2e1 * t399 + (qJD(1) * t350 - g(3)) * t409 + (-qJD(1) * t330 - t540) * t406) * MDP(16) + (-t306 * pkin(3) + g(3) * t499 - t305 * qJ(4) - t329 * t352 - t330 * t350 - t331 * t536 - t447 * t395 - t464) * MDP(17) + (-t405 * t513 - t519) * MDP(18) + ((-t303 - t513) * t408 + (t302 + t516) * t405) * MDP(19) + ((-t373 * t405 * t409 + t346 * t406) * qJD(1) + t538) * MDP(20) + ((-t344 * t406 - t373 * t507) * qJD(1) + t429) * MDP(21) + t373 * MDP(22) * t493 + (t288 * t493 + qJ(4) * t303 + t358 + t482 * t344 + (-t505 + (t317 - t324) * t373) * t408 + (-t472 - t526 + t297 + (t323 - t483) * t373) * t405) * MDP(23) + (-qJ(4) * t302 + t503 * t373 - t289 * t493 + t482 * t346 + (-t317 * t373 + t505) * t405 + (t297 + t430 - t541) * t408) * MDP(24) + (-t282 * t493 + t295 * t373 + t303 * t354 + t358 + t502 * t344 + t421 * t408 + (-t426 - t472) * t405) * MDP(25) + (t294 * t344 - t295 * t346 + (t302 * t411 + (-t344 * t411 - t283) * qJD(5) + t448) * t408 + (-t303 * t411 + (t346 * t411 - t282) * qJD(5) + t449) * t405 + t468) * MDP(26) + (t283 * t493 - t294 * t373 + t302 * t354 - t502 * t346 + t421 * t405 + (t426 + t541) * t408) * MDP(27) + (t279 * t354 - t283 * t294 - t282 * t295 - t464 - g(3) * (t409 * t445 + t454) - (pkin(8) * t409 + t406 * t445) * t396 + t502 * t291 + (qJD(5) * t437 - t439) * t411 - (-t354 * t406 + t409 * t411) * t395) * MDP(28); MDP(14) * t387 + (qJDD(3) - t470) * MDP(15) + (-t404 * t414 - t413) * MDP(16) + (t348 + t422 - t521) * MDP(17) - t328 * MDP(23) - t468 * MDP(28) + (-t344 * t474 - t534) * qJD(3) + (-t343 * MDP(25) + t545 + (qJD(5) * t283 - t448) * MDP(28) - t473 * t542) * t408 + ((t346 * t492 - t303 + t486) * MDP(26) + (qJD(5) * t282 - t449) * MDP(28) + t473 * t343 - t474 * t542) * t405; MDP(18) * t517 + (-t344 ^ 2 + t530) * MDP(19) - t544 * MDP(20) + (t513 - t303) * MDP(21) - t343 * MDP(22) + (-t317 * t346 + t419 + t520) * MDP(23) + (t288 * t373 + t317 * t344 - t533) * MDP(24) + (-t307 * t344 - t416 + t520 - 0.2e1 * t528) * MDP(25) + (pkin(5) * t302 - qJ(6) * t303 + (t283 - t289) * t346 + (t282 - t481) * t344) * MDP(26) + (-0.2e1 * t522 - t291 * t344 + t307 * t346 + (0.2e1 * qJD(6) - t288) * t373 + t533) * MDP(27) + (t277 * qJ(6) - t278 * pkin(5) - t291 * t307 - t282 * t289 - g(1) * (-pkin(5) * t335 + qJ(6) * t336) - g(2) * (-pkin(5) * t333 + qJ(6) * t334) - t446 * t394 + t481 * t283) * MDP(28); (t343 + t517) * MDP(25) - t545 + (-t530 - t542) * MDP(27) + (-t283 * t373 + t416 + t528) * MDP(28);];
tau  = t1;
