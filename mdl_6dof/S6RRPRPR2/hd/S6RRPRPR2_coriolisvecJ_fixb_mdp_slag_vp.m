% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPR2_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR2_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRPR2_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:14:05
% EndTime: 2019-03-09 10:14:14
% DurationCPUTime: 5.29s
% Computational Cost: add. (5508->405), mult. (14229->507), div. (0->0), fcn. (10896->8), ass. (0->192)
t470 = cos(qJ(6));
t469 = sin(qJ(2));
t471 = cos(qJ(2));
t557 = sin(pkin(10));
t558 = cos(pkin(10));
t489 = t469 * t557 - t471 * t558;
t438 = t489 * qJD(1);
t488 = -t469 * t558 - t471 * t557;
t439 = t488 * qJD(1);
t468 = sin(qJ(4));
t479 = qJD(2) * t439;
t519 = qJD(1) * qJD(2);
t492 = t489 * t519;
t564 = cos(qJ(4));
t515 = qJD(4) * t564;
t526 = qJD(4) * t468;
t494 = -t438 * t515 + t439 * t526 + t468 * t479 - t492 * t564;
t348 = t470 * t494;
t497 = t438 * t468 + t439 * t564;
t572 = qJD(6) - t497;
t467 = sin(qJ(6));
t586 = t572 * t467;
t573 = t572 * t586 - t348;
t541 = t467 * t494;
t585 = t572 * t470;
t490 = -t572 * t585 - t541;
t463 = qJD(2) + qJD(4);
t525 = qJD(6) * t467;
t394 = -t438 * t564 + t439 * t468;
t482 = qJD(2) * t488;
t477 = t564 * t482;
t475 = qJD(1) * t477 + qJD(4) * t497 + t468 * t492;
t524 = qJD(6) * t470;
t532 = -t394 * t524 - t467 * t475;
t329 = -t463 * t525 + t532;
t328 = t329 * t470;
t349 = t470 * t475;
t379 = -t394 * t467 + t463 * t470;
t330 = qJD(6) * t379 + t349;
t545 = t394 * t463;
t486 = t494 - t545;
t547 = t497 * t463;
t542 = t463 * t467;
t376 = t394 * t470 + t542;
t587 = t376 * t572;
t589 = (t475 - t547) * MDP(16) + t486 * MDP(15) + (-t379 * t586 + t328) * MDP(24) + (-t379 * t585 + (-t329 + t587) * t467 - t470 * t330) * MDP(25);
t560 = -qJ(3) - pkin(7);
t455 = t560 * t471;
t451 = qJD(1) * t455;
t441 = t557 * t451;
t454 = t560 * t469;
t450 = qJD(1) * t454;
t559 = qJD(2) * pkin(2);
t445 = t450 + t559;
t397 = t445 * t558 + t441;
t561 = t439 * pkin(8);
t373 = qJD(2) * pkin(3) + t397 + t561;
t512 = t558 * t451;
t398 = t445 * t557 - t512;
t562 = t438 * pkin(8);
t374 = t398 - t562;
t339 = t373 * t468 + t374 * t564;
t336 = -qJ(5) * t463 - t339;
t563 = t394 * pkin(5);
t320 = -t336 + t563;
t588 = t320 * t572;
t581 = t497 ^ 2;
t580 = pkin(4) * t497;
t579 = pkin(5) * t497;
t556 = qJ(5) * t394;
t578 = t320 * t497;
t518 = -t471 * pkin(2) - pkin(1);
t500 = t518 * qJD(1);
t453 = qJD(3) + t500;
t403 = pkin(3) * t438 + t453;
t478 = qJ(5) * t497 + t403;
t340 = -pkin(4) * t394 + t478;
t577 = t340 * t497;
t576 = t403 * t497;
t565 = pkin(4) + pkin(9);
t575 = t497 * t565;
t401 = -t450 * t557 + t512;
t380 = t401 + t562;
t402 = t450 * t558 + t441;
t381 = t402 + t561;
t517 = pkin(2) * t557;
t457 = t468 * t517;
t516 = t558 * pkin(2);
t459 = t516 + pkin(3);
t574 = qJD(4) * t457 + t380 * t468 + t381 * t564 - t459 * t515;
t338 = -t373 * t564 + t374 * t468;
t522 = qJD(5) + t338;
t571 = -0.2e1 * t519;
t570 = MDP(5) * (t469 ^ 2 - t471 ^ 2);
t487 = t459 * t468 + t517 * t564;
t534 = -qJD(4) * t487 - t380 * t564 + t381 * t468;
t533 = -qJD(5) + t574;
t523 = t522 - t579;
t513 = qJD(2) * t560;
t435 = qJD(3) * t471 + t469 * t513;
t416 = t435 * qJD(1);
t436 = -qJD(3) * t469 + t471 * t513;
t417 = t436 * qJD(1);
t375 = -t416 * t557 + t417 * t558;
t359 = pkin(8) * t492 + t375;
t378 = t416 * t558 + t417 * t557;
t360 = pkin(8) * t479 + t378;
t313 = -t359 * t564 + t360 * t468 + t373 * t526 + t374 * t515;
t568 = -t463 * t534 + t313;
t318 = -t463 * t565 + t523;
t323 = -t394 * t565 + t478;
t307 = t318 * t470 - t323 * t467;
t308 = t318 * t467 + t323 * t470;
t567 = -MDP(13) * t497 + MDP(19) * t403 - MDP(22) * t340 + MDP(28) * t572 + MDP(29) * t307 - MDP(30) * t308;
t386 = -t435 * t557 + t436 * t558;
t483 = qJD(2) * t489;
t365 = pkin(8) * t483 + t386;
t387 = t435 * t558 + t436 * t557;
t366 = pkin(8) * t482 + t387;
t405 = t454 * t558 + t455 * t557;
t388 = pkin(8) * t488 + t405;
t406 = t454 * t557 - t455 * t558;
t389 = -pkin(8) * t489 + t406;
t315 = -t365 * t468 - t366 * t564 - t388 * t515 + t389 * t526;
t555 = t315 * t463;
t498 = t388 * t468 + t389 * t564;
t316 = qJD(4) * t498 - t365 * t564 + t366 * t468;
t554 = t316 * t463;
t480 = t564 * t489;
t399 = -t468 * t488 + t480;
t484 = t468 * t489;
t400 = -t488 * t564 - t484;
t419 = pkin(3) * t489 + t518;
t476 = -t400 * qJ(5) + t419;
t331 = t399 * t565 + t476;
t553 = t331 * t494;
t552 = t339 * t463;
t551 = t376 * t394;
t550 = t379 * t394;
t546 = t394 ^ 2;
t544 = t399 * t467;
t433 = -t459 * t564 - pkin(4) + t457;
t429 = -pkin(9) + t433;
t543 = t429 * t494;
t473 = qJD(2) ^ 2;
t540 = t469 * t473;
t539 = t471 * t473;
t474 = qJD(1) ^ 2;
t538 = t471 * t474;
t537 = t565 * t494;
t502 = t359 * t468 + t360 * t564 + t373 * t515 - t374 * t526;
t312 = -qJD(5) * t463 - t502;
t303 = pkin(5) * t475 - t312;
t536 = t303 * t467 + t320 * t524;
t535 = t534 + t563;
t531 = -t533 - t579;
t527 = qJD(1) * t469;
t461 = t469 * t559;
t514 = t469 * t519;
t410 = pkin(2) * t527 - pkin(3) * t439;
t509 = pkin(1) * t571;
t499 = t410 - t556;
t504 = -qJD(6) * t429 + t499 - t575;
t503 = qJD(6) * t565 - t556 - t575;
t345 = -t388 * t564 + t389 * t468;
t354 = -qJD(4) * t484 - t468 * t483 - t488 * t515 - t477;
t496 = t354 * t467 + t399 * t524;
t495 = t313 - t577;
t332 = pkin(5) * t400 + t345;
t491 = t303 * t399 + t320 * t354 - t332 * t494;
t411 = -pkin(3) * t482 + t461;
t458 = pkin(2) * t514;
t404 = -pkin(3) * t479 + t458;
t353 = t463 * t480 - t468 * t482 - t488 * t526;
t317 = t354 * pkin(4) + t353 * qJ(5) - t400 * qJD(5) + t411;
t314 = -pkin(4) * t475 - qJ(5) * t494 + qJD(5) * t497 + t404;
t432 = qJ(5) + t487;
t352 = -t556 - t580;
t344 = t399 * pkin(4) + t476;
t343 = t499 - t580;
t337 = t494 * t400;
t335 = -pkin(4) * t463 + t522;
t333 = -pkin(5) * t399 + t498;
t322 = t339 + t563;
t311 = t354 * pkin(9) + t317;
t310 = -pkin(5) * t353 + t316;
t309 = -pkin(5) * t354 - t315;
t306 = -pkin(9) * t475 + t314;
t305 = pkin(5) * t494 + t313;
t304 = t470 * t305;
t302 = t303 * t470;
t1 = [MDP(6) * t539 + (t304 * t400 - t307 * t353 + t309 * t376 + t333 * t330 + (-t306 * t400 - t311 * t572 - t553) * t467 + (t310 * t572 - t491) * t470 + ((-t331 * t470 - t332 * t467) * t572 - t308 * t400 + t320 * t544) * qJD(6)) * MDP(29) + (t308 * t353 + t309 * t379 + t333 * t329 + (-(qJD(6) * t332 + t311) * t572 - t553 - (qJD(6) * t318 + t306) * t400 + t320 * qJD(6) * t399) * t470 + (-(-qJD(6) * t331 + t310) * t572 - (-qJD(6) * t323 + t305) * t400 + t491) * t467) * MDP(30) + ((-t376 * t467 + t379 * t470) * t354 + (t328 - t330 * t467 + (-t376 * t470 - t379 * t467) * qJD(6)) * t399) * MDP(25) + (t354 * t403 - t394 * t411 + t399 * t404 - t419 * t475 - t554) * MDP(18) + (-t314 * t399 + t317 * t394 - t340 * t354 + t344 * t475 + t554) * MDP(21) + (-t314 * t400 + t317 * t497 + t340 * t353 - t344 * t494 - t555) * MDP(22) + (-t353 * t403 + t400 * t404 - t411 * t497 + t419 * t494 + t555) * MDP(19) + (t329 * t544 + t379 * t496) * MDP(24) + (t399 * t348 - t330 * t400 + t353 * t376 + (t354 * t470 - t399 * t525) * t572) * MDP(27) - MDP(7) * t540 + (pkin(7) * t540 + t471 * t509) * MDP(10) + (t329 * t400 - t353 * t379 + t399 * t541 + t496 * t572) * MDP(26) + (-pkin(7) * t539 + t469 * t509) * MDP(9) + (t375 * t405 + t378 * t406 + t397 * t386 + t398 * t387 + (t453 + t500) * t461) * MDP(12) + (t375 * t488 + t386 * t439 - t387 * t438 + t398 * t482 + t405 * t492 + t406 * t479 + (qJD(2) * t397 - t378) * t489) * MDP(11) + (-t353 * t394 + t354 * t497 - t399 * t494 + t400 * t475) * MDP(14) + (t312 * t399 + t313 * t400 - t315 * t394 - t316 * t497 - t335 * t353 + t336 * t354 + t345 * t494 + t475 * t498) * MDP(20) + 0.2e1 * t471 * MDP(4) * t514 + t570 * t571 + (t353 * t497 + t337) * MDP(13) + (-t353 * t572 + t337) * MDP(28) + (-t312 * t498 + t313 * t345 + t314 * t344 + t315 * t336 + t316 * t335 + t317 * t340) * MDP(23) + (-MDP(15) * t353 - MDP(16) * t354) * t463; (-t343 * t497 - t463 * t533 - t312) * MDP(22) + (-MDP(20) * t335 - t567) * t394 + t589 + (t432 * t329 + t302 + t504 * t585 + t531 * t379 + (t535 * t572 - t543 - t588) * t467) * MDP(30) + (-t546 + t581) * MDP(14) + (t479 * t517 + t492 * t516 + (-t398 - t401) * t439 - (-t402 + t397) * t438) * MDP(11) + (t475 * t432 - t394 * t533 + t433 * t494 + (t336 + t534) * t497) * MDP(20) + (t432 * t330 + (t543 - t578) * t470 + t531 * t376 + (t467 * t504 - t470 * t535) * t572 + t536) * MDP(29) + (-t397 * t401 - t398 * t402 + (t375 * t558 + t378 * t557 - t453 * t527) * pkin(2)) * MDP(12) + (t490 + t551) * MDP(27) + (t394 * t410 - t568 + t576) * MDP(18) + (-t550 - t573) * MDP(26) + (-t312 * t432 + t313 * t433 - t335 * t534 + t336 * t533 - t340 * t343) * MDP(23) + (-t343 * t394 + t568 - t577) * MDP(21) + (MDP(9) * t469 * t474 + MDP(10) * t538) * pkin(1) + (t410 * t497 + t463 * t574 - t502) * MDP(19) - t469 * MDP(4) * t538 + t474 * t570; (-t438 ^ 2 - t439 ^ 2) * MDP(11) + (-t397 * t439 + t398 * t438 + t458) * MDP(12) + (-t546 - t581) * MDP(20) + (t335 * t497 + t336 * t394 + t314) * MDP(23) + (t490 - t551) * MDP(29) + (-t550 + t573) * MDP(30) + (MDP(19) - MDP(22)) * (t494 + t545) + (-MDP(18) + MDP(21)) * (t475 + t547); t581 * MDP(14) + (-t313 + t552 + t576) * MDP(18) + (-t338 * t463 - t502) * MDP(19) + (-pkin(4) * t494 + qJ(5) * t475 - (-t336 - t339) * t497) * MDP(20) + (t495 - t552) * MDP(21) + (-t352 * t497 + t463 * t522 - t312) * MDP(22) + (-pkin(4) * t313 - qJ(5) * t312 - t335 * t339 - t336 * t522 - t340 * t352) * MDP(23) - t573 * MDP(26) + t490 * MDP(27) + (qJ(5) * t330 + (-t537 - t578) * t470 + (-t322 * t470 + t467 * t503) * t572 + t523 * t376 + t536) * MDP(29) + (qJ(5) * t329 + t302 + t503 * t585 + t523 * t379 + (t322 * t572 + t537 - t588) * t467) * MDP(30) - ((t335 - t522) * MDP(20) + t352 * MDP(21) + t379 * MDP(26) - t376 * MDP(27) + MDP(14) * t394 + t567) * t394 + t589; t486 * MDP(20) - t497 * t394 * MDP(21) + (-t463 ^ 2 - t581) * MDP(22) + (t336 * t463 + t495) * MDP(23) + (-t376 * t463 + t348) * MDP(29) + (-t379 * t463 - t541) * MDP(30) + (-MDP(29) * t586 - MDP(30) * t585) * t572; t379 * t376 * MDP(24) + (-t376 ^ 2 + t379 ^ 2) * MDP(25) + (t532 + t587) * MDP(26) + (t379 * t572 - t349) * MDP(27) + t494 * MDP(28) + (-t306 * t467 + t308 * t572 - t320 * t379 + t304) * MDP(29) + (-t305 * t467 - t306 * t470 + t307 * t572 + t320 * t376) * MDP(30) + (-MDP(26) * t542 - MDP(27) * t379 - MDP(29) * t308 - MDP(30) * t307) * qJD(6);];
tauc  = t1;
