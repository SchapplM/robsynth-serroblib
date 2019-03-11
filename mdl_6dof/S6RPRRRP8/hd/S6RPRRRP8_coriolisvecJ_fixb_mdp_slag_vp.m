% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% MDP [31x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP8_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP8_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(31,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [31 1]), ...
  'S6RPRRRP8_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [31x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:24:53
% EndTime: 2019-03-09 06:24:59
% DurationCPUTime: 4.35s
% Computational Cost: add. (5679->416), mult. (11679->539), div. (0->0), fcn. (7654->6), ass. (0->180)
t419 = sin(qJ(3));
t421 = cos(qJ(3));
t529 = sin(qJ(4));
t530 = cos(qJ(4));
t384 = t530 * t419 + t529 * t421;
t377 = t384 * qJD(1);
t539 = qJD(5) + t377;
t544 = t539 ^ 2;
t477 = t529 * t419;
t464 = qJD(1) * t477;
t478 = t530 * t421;
t376 = -qJD(1) * t478 + t464;
t418 = sin(qJ(5));
t420 = cos(qJ(5));
t485 = qJD(3) + qJD(4);
t468 = t420 * t485;
t354 = -t376 * t418 - t468;
t543 = t354 * t539;
t438 = t420 * t376 - t418 * t485;
t542 = t438 * t539;
t491 = qJD(5) * t420;
t511 = t377 * t420;
t541 = t491 + t511;
t492 = qJD(5) * t418;
t512 = t377 * t418;
t540 = -t492 - t512;
t475 = t530 * qJD(4);
t436 = t530 * qJD(3) + t475;
t435 = t436 * t421;
t499 = t485 * t464;
t346 = qJD(1) * t435 - t499;
t508 = t420 * t346;
t538 = (t492 * t539 - t508) * pkin(9);
t537 = MDP(8) * (t419 ^ 2 - t421 ^ 2);
t536 = t377 * t485;
t422 = -pkin(1) - pkin(7);
t394 = t422 * qJD(1) + qJD(2);
t495 = qJD(1) * t421;
t372 = -pkin(8) * t495 + t421 * t394;
t367 = qJD(3) * pkin(3) + t372;
t471 = pkin(8) * qJD(1) - t394;
t494 = qJD(3) * t419;
t368 = t471 * t494;
t493 = qJD(3) * t421;
t369 = t471 * t493;
t496 = qJD(1) * t419;
t371 = -pkin(8) * t496 + t394 * t419;
t473 = t529 * qJD(4);
t303 = t367 * t475 + t529 * t368 - t530 * t369 - t371 * t473;
t414 = qJD(1) * qJD(2);
t488 = qJD(1) * qJD(3);
t472 = t421 * t488;
t386 = pkin(3) * t472 + t414;
t307 = t346 * pkin(4) + pkin(9) * t536 + t386;
t366 = t530 * t371;
t336 = t529 * t367 + t366;
t329 = t485 * pkin(9) + t336;
t390 = pkin(3) * t496 + qJD(1) * qJ(2);
t337 = pkin(4) * t377 + pkin(9) * t376 + t390;
t440 = t420 * t303 + t418 * t307 - t329 * t492 + t337 * t491;
t523 = qJ(6) * t346;
t275 = qJD(6) * t539 + t440 + t523;
t465 = t418 * t303 - t420 * t307 + t329 * t491 + t337 * t492;
t528 = pkin(5) * t346;
t277 = t465 - t528;
t459 = t275 * t420 + t277 * t418;
t383 = t477 - t478;
t405 = t419 * pkin(3) + qJ(2);
t348 = pkin(4) * t384 + pkin(9) * t383 + t405;
t526 = pkin(8) - t422;
t387 = t526 * t419;
t388 = t526 * t421;
t353 = -t530 * t387 - t529 * t388;
t500 = t418 * t348 + t420 * t353;
t352 = -t529 * t387 + t530 * t388;
t535 = t540 * pkin(5) + t541 * qJ(6) + t418 * qJD(6);
t339 = t529 * t372 + t366;
t449 = pkin(3) * t473 - t339;
t534 = -t529 * qJD(3) - t473;
t533 = qJ(2) * MDP(6) + MDP(5);
t532 = MDP(26) + MDP(28);
t486 = MDP(27) - MDP(30);
t317 = -t438 * qJD(5) - t418 * t536;
t531 = t438 ^ 2;
t527 = pkin(5) * t376;
t525 = pkin(3) * qJD(4);
t304 = t367 * t473 - t530 * t368 - t529 * t369 + t371 * t475;
t316 = -qJD(5) * t468 - t376 * t492 + t420 * t536;
t279 = pkin(5) * t317 + qJ(6) * t316 + qJD(6) * t438 + t304;
t522 = t279 * t418;
t365 = t529 * t371;
t335 = t530 * t367 - t365;
t328 = -t485 * pkin(4) - t335;
t295 = t354 * pkin(5) + qJ(6) * t438 + t328;
t521 = t295 * t438;
t520 = t316 * t418;
t519 = t317 * t420;
t407 = t529 * pkin(3) + pkin(9);
t518 = t346 * t407;
t517 = t354 * t418;
t516 = t354 * t420;
t515 = t438 * t354;
t514 = t438 * t418;
t513 = t438 * t420;
t510 = t383 * t420;
t509 = t418 * t346;
t424 = qJD(1) ^ 2;
t507 = t421 * t424;
t423 = qJD(3) ^ 2;
t506 = t422 * t423;
t505 = t336 + t535;
t504 = -t449 + t535;
t347 = -pkin(4) * t376 + pkin(9) * t377;
t502 = t420 * t335 + t418 * t347;
t340 = t530 * t372 - t365;
t341 = pkin(3) * t495 + t347;
t501 = t420 * t340 + t418 * t341;
t395 = pkin(3) * t493 + qJD(2);
t298 = -t329 * t418 + t337 * t420;
t489 = qJD(6) - t298;
t484 = 0.2e1 * qJD(1);
t483 = t530 * pkin(3);
t481 = t418 * t530;
t480 = t420 * t530;
t470 = t420 * t539;
t350 = -t436 * t419 + t534 * t421;
t469 = t350 * t485;
t467 = pkin(3) * t475;
t299 = t329 * t420 + t337 * t418;
t462 = -t299 * t376 + t304 * t418 + t328 * t491;
t461 = -t420 * pkin(5) - t418 * qJ(6);
t460 = pkin(5) * t418 - qJ(6) * t420;
t289 = -pkin(5) * t539 + t489;
t290 = qJ(6) * t539 + t299;
t458 = t289 * t420 - t290 * t418;
t457 = t289 * t418 + t290 * t420;
t456 = -t519 - t520;
t454 = t328 * t377 - t518;
t453 = -t335 * t418 + t347 * t420;
t452 = -t513 + t517;
t389 = -pkin(4) + t461;
t448 = t290 * t376 - t295 * t511 - t522;
t447 = -t279 * t420 - t289 * t376 + t295 * t492;
t446 = t298 * t376 - t304 * t420 + t328 * t492;
t444 = -t350 * t418 + t383 * t491;
t443 = t350 * t420 + t383 * t492;
t442 = t299 * t539 - t465;
t441 = t390 * t376 - t304;
t351 = t534 * t419 + t435;
t314 = pkin(4) * t351 - pkin(9) * t350 + t395;
t381 = t526 * t494;
t382 = qJD(3) * t388;
t318 = -t352 * qJD(4) + t529 * t381 - t530 * t382;
t439 = t418 * t314 + t420 * t318 + t348 * t491 - t353 * t492;
t437 = (-t491 * t539 - t509) * pkin(9);
t434 = -t407 * t492 + t420 * t467;
t433 = t541 * t289 + t540 * t290 + t459;
t432 = t458 * qJD(5) + t459;
t431 = t452 * qJD(5) + t456;
t430 = -t436 * t495 + t499;
t429 = ((-t316 - t543) * t420 + (-t317 + t542) * t418) * MDP(22) + (-t438 * t470 - t520) * MDP(21) + (-t354 * t376 - t418 * t544 + t508) * MDP(24) + (-t376 * t438 + t470 * t539 + t509) * MDP(23) + (-t376 * t485 + t430) * MDP(17) + (t376 ^ 2 - t377 ^ 2) * MDP(15) + (-MDP(14) * t377 + MDP(25) * t539) * t376;
t428 = t452 * MDP(29) + t458 * MDP(31);
t427 = t390 * t377 - t303;
t319 = t353 * qJD(4) - t530 * t381 - t529 * t382;
t408 = -t483 - pkin(4);
t379 = -t483 + t389;
t373 = t376 * qJ(6);
t324 = -pkin(5) * t438 + qJ(6) * t354;
t320 = -t460 * t383 + t352;
t310 = -pkin(5) * t384 - t348 * t420 + t353 * t418;
t309 = qJ(6) * t384 + t500;
t297 = -t453 + t527;
t296 = -t373 + t502;
t294 = t340 * t418 - t341 * t420 + t527;
t293 = -t373 + t501;
t291 = -t316 + t543;
t282 = t460 * t350 + (t461 * qJD(5) + qJD(6) * t420) * t383 + t319;
t281 = -pkin(5) * t351 + t500 * qJD(5) - t314 * t420 + t318 * t418;
t280 = qJ(6) * t351 + qJD(6) * t384 + t439;
t1 = [(-t277 * t384 - t281 * t539 + t282 * t354 - t289 * t351 - t444 * t295 - t310 * t346 + t317 * t320 - t383 * t522) * MDP(28) + ((t514 - t516) * t350 + (-t520 + t519 + (-t513 - t517) * qJD(5)) * t383) * MDP(22) + (-t299 * t351 - t304 * t510 - t352 * t316 - t319 * t438 + t443 * t328 - t500 * t346 - t440 * t384 - t439 * t539) * MDP(27) + (t316 * t510 - t438 * t443) * MDP(21) + (t275 * t384 + t279 * t510 + t280 * t539 + t282 * t438 + t290 * t351 - t443 * t295 + t309 * t346 + t316 * t320) * MDP(30) + (-t316 * t384 - t351 * t438 - t383 * t508 + t443 * t539) * MDP(23) + (-t317 * t384 - t351 * t354 + t383 * t509 + t444 * t539) * MDP(24) + (-t421 * t506 + (-qJ(2) * t494 + qJD(2) * t421) * t484) * MDP(13) + qJ(2) * t493 * t484 * MDP(12) + (-t318 * t485 + t390 * t350 - t395 * t376 - t386 * t383 - t405 * t536) * MDP(20) + (t383 * t346 - t350 * t377 + t376 * t351 + t384 * t536) * MDP(15) + (-t376 * t350 + t383 * t536) * MDP(14) + (-t280 * t354 - t281 * t438 - t309 * t317 - t310 * t316 + t458 * t350 + (t457 * qJD(5) + t275 * t418 - t277 * t420) * t383) * MDP(29) - t423 * t421 * MDP(10) + (-t319 * t485 + t405 * t346 + t390 * t351 + t395 * t377 + t386 * t384) * MDP(19) + (-t465 * t384 + t298 * t351 + t319 * t354 + t352 * t317 + ((-qJD(5) * t353 + t314) * t539 + t348 * t346 - t328 * qJD(5) * t383) * t420 + ((-qJD(5) * t348 - t318) * t539 - t353 * t346 - t304 * t383 + t328 * t350) * t418) * MDP(26) + 0.2e1 * t488 * t537 - t351 * t485 * MDP(17) + MDP(16) * t469 + (t346 * t384 + t351 * t539) * MDP(25) + (t275 * t309 + t277 * t310 + t279 * t320 + t280 * t290 + t281 * t289 + t282 * t295) * MDP(31) + 0.2e1 * t533 * t414 + ((qJD(2) * t484 - t506) * MDP(12) - t423 * MDP(9) - 0.2e1 * MDP(7) * t472) * t419; MDP(19) * t469 + (t279 * t383 - t295 * t350) * MDP(31) - t533 * t424 + (-t485 * MDP(20) + (-t514 - t516) * MDP(29) + t457 * MDP(31)) * t351 + (-t377 * MDP(19) + t376 * MDP(20) + t428) * qJD(1) + (t532 * (-qJD(1) * t420 - t351 * t418) + t486 * (qJD(1) * t418 - t351 * t420)) * t539 + (t456 * MDP(29) + t459 * MDP(31) + (-t532 * t418 - t486 * t420) * t346 + ((t486 * t418 - t532 * t420) * t539 + t428) * qJD(5)) * t384 + t532 * (t383 * t317 - t350 * t354) + (MDP(12) * t419 + MDP(13) * t421) * (-t423 - t424) - t486 * (t316 * t383 - t350 * t438); (t339 * t485 + (-t377 * t495 - t485 * t473) * pkin(3) + t441) * MDP(19) + (t408 * t317 + t454 * t418 + t449 * t354 + ((-qJD(5) * t407 - t341) * t420 + (-t467 + t340) * t418) * t539 + t446) * MDP(26) + (-t408 * t316 + t454 * t420 - t449 * t438 + (-t434 + t501) * t539 + t462) * MDP(27) + (t379 * t317 + (t295 * t377 - t518) * t418 - t504 * t354 + (-t407 * t491 - t418 * t467 + t294) * t539 + t447) * MDP(28) + (t293 * t354 + t294 * t438 + (-t354 * t480 - t438 * t481) * t525 + t431 * t407 + t433) * MDP(29) + (t379 * t316 + (-qJD(5) * t295 + t518) * t420 - t504 * t438 + (-t293 + t434) * t539 + t448) * MDP(30) + (t340 * t485 + (t376 * t495 - t485 * t475) * pkin(3) + t427) * MDP(20) + t429 + (t279 * t379 - t289 * t294 - t290 * t293 - t504 * t295 + (t289 * t481 + t290 * t480) * t525 + t432 * t407) * MDP(31) - t424 * t537 + t419 * MDP(7) * t507 + (t424 * t419 * MDP(13) - MDP(12) * t507) * qJ(2); (-pkin(4) * t317 + t328 * t512 - t336 * t354 - t453 * t539 + t437 + t446) * MDP(26) + (pkin(4) * t316 + t328 * t511 + t336 * t438 + t502 * t539 + t462 + t538) * MDP(27) + (t295 * t512 + t297 * t539 + t317 * t389 - t505 * t354 + t437 + t447) * MDP(28) + (t431 * pkin(9) + t296 * t354 + t297 * t438 + t433) * MDP(29) + (-t295 * t491 - t296 * t539 + t316 * t389 - t438 * t505 + t448 - t538) * MDP(30) + (t432 * pkin(9) + t279 * t389 - t289 * t297 - t290 * t296 - t505 * t295) * MDP(31) + (t335 * t485 + t427) * MDP(20) + t429 + (t336 * t485 + t441) * MDP(19); -MDP(21) * t515 + (-t354 ^ 2 + t531) * MDP(22) + t291 * MDP(23) + (-t317 - t542) * MDP(24) + t346 * MDP(25) + (t328 * t438 + t442) * MDP(26) + (t298 * t539 + t328 * t354 - t440) * MDP(27) + (-t324 * t354 + t442 + t521 + 0.2e1 * t528) * MDP(28) + (pkin(5) * t316 - qJ(6) * t317 - (t290 - t299) * t438 + (t289 - t489) * t354) * MDP(29) + (0.2e1 * t523 - t295 * t354 - t324 * t438 + (0.2e1 * qJD(6) - t298) * t539 + t440) * MDP(30) + (-pkin(5) * t277 + qJ(6) * t275 - t289 * t299 + t489 * t290 - t295 * t324) * MDP(31); (t430 - t515) * MDP(28) + t291 * MDP(29) + (-t531 - t544) * MDP(30) + (-t290 * t539 + t277 - t521) * MDP(31);];
tauc  = t1;
