% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:16:52
% EndTime: 2019-03-09 06:17:04
% DurationCPUTime: 5.91s
% Computational Cost: add. (6905->422), mult. (17636->549), div. (0->0), fcn. (13591->8), ass. (0->180)
t487 = cos(pkin(10));
t492 = cos(qJ(3));
t541 = qJD(1) * t492;
t474 = t487 * t541;
t486 = sin(pkin(10));
t490 = sin(qJ(3));
t542 = qJD(1) * t490;
t526 = t486 * t542;
t442 = t474 - t526;
t454 = t486 * t492 + t487 * t490;
t443 = t454 * qJD(1);
t404 = pkin(3) * t443 - pkin(8) * t442;
t491 = cos(qJ(4));
t394 = t491 * t404;
t489 = sin(qJ(4));
t578 = pkin(8) + pkin(9);
t529 = qJD(4) * t578;
t576 = pkin(7) + qJ(2);
t465 = t576 * t486;
t455 = qJD(1) * t465;
t466 = t576 * t487;
t456 = qJD(1) * t466;
t587 = -t455 * t492 - t490 * t456;
t603 = pkin(4) * t443 - t489 * t587 + t394 + (-pkin(9) * t442 + t529) * t491;
t546 = t489 * t404 + t491 * t587;
t565 = t442 * t489;
t602 = -pkin(9) * t565 + t489 * t529 + t546;
t538 = qJD(4) * t489;
t601 = t538 - t565;
t471 = qJD(3) * t474;
t509 = qJD(3) * t526 - t471;
t595 = qJD(3) * qJD(4) - t509;
t376 = -t443 * t538 + t491 * t595;
t419 = qJD(3) * t491 - t489 * t443;
t420 = qJD(3) * t489 + t443 * t491;
t488 = sin(qJ(5));
t577 = cos(qJ(5));
t522 = t577 * qJD(5);
t537 = qJD(4) * t491;
t530 = t443 * t537 + t489 * t595;
t536 = qJD(5) * t488;
t328 = -t577 * t376 - t419 * t522 + t420 * t536 + t488 * t530;
t508 = t488 * t419 + t420 * t577;
t329 = qJD(5) * t508 + t488 * t376 + t577 * t530;
t364 = -t577 * t419 + t420 * t488;
t361 = t364 ^ 2;
t437 = qJD(4) - t442;
t432 = qJD(5) + t437;
t445 = t454 * qJD(3);
t433 = qJD(1) * t445;
t579 = t508 ^ 2;
t600 = t433 * MDP(26) + (t432 * t508 - t329) * MDP(25) + t364 * t508 * MDP(22) + (t364 * t432 - t328) * MDP(24) + (-t361 + t579) * MDP(23);
t599 = t364 * qJ(6);
t598 = (t486 ^ 2 + t487 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t408 = -t490 * t455 + t492 * t456;
t588 = t601 * pkin(4) - t408;
t528 = t577 * t489;
t458 = t488 * t491 + t528;
t584 = qJD(4) + qJD(5);
t413 = t584 * t458;
t548 = t458 * t442 - t413;
t556 = t488 * t489;
t507 = t577 * t491 - t556;
t585 = t577 * qJD(4) + t522;
t547 = t507 * t442 - t491 * t585 + t556 * t584;
t597 = t454 * qJD(2);
t467 = t578 * t489;
t468 = t578 * t491;
t544 = -t488 * t467 + t577 * t468;
t596 = -qJD(5) * t544 + t602 * t488 - t603 * t577;
t453 = t486 * t490 - t492 * t487;
t444 = t453 * qJD(3);
t524 = t454 * t537;
t594 = -t444 * t489 + t524;
t402 = -qJD(3) * pkin(3) - t587;
t360 = -pkin(4) * t419 + t402;
t478 = -pkin(2) * t487 - pkin(1);
t464 = qJD(1) * t478 + qJD(2);
t380 = -pkin(3) * t442 - pkin(8) * t443 + t464;
t403 = qJD(3) * pkin(8) + t408;
t356 = t380 * t489 + t403 * t491;
t499 = t453 * qJD(2);
t369 = -qJD(1) * t499 + qJD(3) * t587;
t388 = t433 * pkin(3) + pkin(8) * t509;
t385 = t491 * t388;
t498 = -qJD(4) * t356 - t369 * t489 + t385;
t314 = pkin(4) * t433 - pkin(9) * t376 + t498;
t503 = t491 * t369 + t380 * t537 + t489 * t388 - t403 * t538;
t320 = -pkin(9) * t530 + t503;
t355 = t491 * t380 - t403 * t489;
t342 = -pkin(9) * t420 + t355;
t332 = pkin(4) * t437 + t342;
t343 = pkin(9) * t419 + t356;
t515 = t488 * t314 + t577 * t320 + t332 * t522 - t343 * t536;
t593 = t360 * t364 - t515;
t591 = qJ(6) * t508;
t335 = pkin(5) * t364 + qJD(6) + t360;
t590 = t335 * t508;
t406 = pkin(3) * t453 - pkin(8) * t454 + t478;
t399 = t491 * t406;
t418 = -t465 * t490 + t466 * t492;
t561 = t454 * t491;
t351 = pkin(4) * t453 - pkin(9) * t561 - t418 * t489 + t399;
t411 = t491 * t418;
t545 = t489 * t406 + t411;
t562 = t454 * t489;
t357 = -pkin(9) * t562 + t545;
t549 = t488 * t351 + t577 * t357;
t417 = t465 * t492 + t490 * t466;
t586 = t467 * t522 + t468 * t536 + t603 * t488 + t602 * t577;
t339 = t577 * t343;
t317 = t488 * t332 + t339;
t496 = -qJD(5) * t317 + t577 * t314 - t488 * t320;
t583 = -t360 * t508 + t496;
t582 = -t432 * t547 + t433 * t458;
t581 = t328 * t507 - t508 * t548;
t573 = t364 * t443;
t572 = t508 * t443;
t571 = t376 * t489;
t570 = t419 * t437;
t569 = t419 * t443;
t568 = t420 * t437;
t567 = t420 * t443;
t563 = t444 * t491;
t337 = t488 * t343;
t555 = t489 * t433;
t423 = t491 * t433;
t316 = t577 * t332 - t337;
t306 = t316 - t591;
t305 = pkin(5) * t432 + t306;
t554 = t305 - t306;
t553 = t577 * t342 - t337;
t552 = qJ(6) * t548 + qJD(6) * t507 - t586;
t551 = -pkin(5) * t443 + qJ(6) * t547 - t458 * qJD(6) + t596;
t540 = qJD(3) * t490;
t539 = qJD(3) * t492;
t533 = qJD(1) * qJD(2);
t483 = -pkin(4) * t491 - pkin(3);
t525 = t454 * t538;
t521 = -t342 * t488 - t339;
t519 = t577 * t351 - t357 * t488;
t517 = -t577 * t467 - t468 * t488;
t516 = t437 * t491;
t370 = qJD(1) * t597 - t455 * t540 + t456 * t539;
t382 = -t465 * t540 + t466 * t539 + t597;
t514 = -t458 * t329 + t364 * t547;
t513 = t432 * t548 + t507 * t433;
t383 = pkin(4) * t562 + t417;
t510 = -t601 * t437 + t423;
t359 = pkin(4) * t594 + t382;
t505 = -t525 - t563;
t504 = -pkin(8) * t433 + t402 * t437;
t381 = -qJD(3) * t417 - t499;
t405 = pkin(3) * t445 + pkin(8) * t444;
t502 = t491 * t381 + t489 * t405 + t406 * t537 - t418 * t538;
t395 = t491 * t405;
t324 = pkin(9) * t563 + pkin(4) * t445 - t381 * t489 + t395 + (-t411 + (pkin(9) * t454 - t406) * t489) * qJD(4);
t326 = -pkin(9) * t594 + t502;
t501 = t488 * t324 + t577 * t326 + t351 * t522 - t357 * t536;
t348 = pkin(4) * t530 + t370;
t311 = t329 * pkin(5) + t348;
t495 = -qJD(5) * t549 + t577 * t324 - t488 * t326;
t482 = pkin(4) * t577 + pkin(5);
t409 = t433 * t453;
t397 = t507 * t454;
t396 = t458 * t454;
t391 = qJ(6) * t507 + t544;
t390 = -qJ(6) * t458 + t517;
t334 = -t444 * t528 - t488 * t525 - t536 * t562 + (-t444 * t488 + t454 * t585) * t491;
t333 = t413 * t454 + t444 * t507;
t321 = -qJ(6) * t396 + t549;
t318 = pkin(5) * t453 - qJ(6) * t397 + t519;
t309 = t553 - t591;
t308 = t521 + t599;
t307 = t317 - t599;
t304 = -qJ(6) * t334 - qJD(6) * t396 + t501;
t303 = t445 * pkin(5) + t333 * qJ(6) - t397 * qJD(6) + t495;
t302 = -qJ(6) * t329 - qJD(6) * t364 + t515;
t301 = t433 * pkin(5) + t328 * qJ(6) - qJD(6) * t508 + t496;
t1 = [(-t443 * t444 - t454 * t509) * MDP(8) + (-t454 * t433 - t444 * t442 - t443 * t445 + t453 * t509) * MDP(9) + (t433 * t478 + t445 * t464) * MDP(13) + (-t464 * t444 - t478 * t509) * MDP(14) + (t376 * t561 + t420 * t505) * MDP(15) + (-(t419 * t491 - t420 * t489) * t444 + (-t491 * t530 - t571 + (-t419 * t489 - t420 * t491) * qJD(4)) * t454) * MDP(16) + (t376 * t453 + t420 * t445 + t423 * t454 + t437 * t505) * MDP(17) + (t419 * t445 - t437 * t594 - t453 * t530 - t454 * t555) * MDP(18) + (t437 * t445 + t409) * MDP(19) + ((-t418 * t537 + t395) * t437 + t399 * t433 + (-t403 * t537 + t385) * t453 + t355 * t445 - t382 * t419 + t417 * t530 + t402 * t524 + ((-qJD(4) * t406 - t381) * t437 - t418 * t433 + (-qJD(4) * t380 - t369) * t453 + t370 * t454 - t402 * t444) * t489) * MDP(20) + (-t356 * t445 + t370 * t561 + t417 * t376 + t382 * t420 + t402 * t505 - t433 * t545 - t437 * t502 - t453 * t503) * MDP(21) + (-t328 * t397 - t333 * t508) * MDP(22) + (t328 * t396 - t329 * t397 + t333 * t364 - t334 * t508) * MDP(23) + (-t328 * t453 - t333 * t432 + t397 * t433 + t445 * t508) * MDP(24) + (-t329 * t453 - t334 * t432 - t364 * t445 - t396 * t433) * MDP(25) + (t432 * t445 + t409) * MDP(26) + (t316 * t445 + t383 * t329 + t360 * t334 + t348 * t396 + t359 * t364 + t432 * t495 + t433 * t519 + t453 * t496) * MDP(27) + (-t317 * t445 - t383 * t328 - t360 * t333 + t348 * t397 + t359 * t508 - t432 * t501 - t433 * t549 - t453 * t515) * MDP(28) + (-t301 * t397 - t302 * t396 - t303 * t508 - t304 * t364 + t305 * t333 - t307 * t334 + t318 * t328 - t321 * t329) * MDP(29) + (t302 * t321 + t307 * t304 + t301 * t318 + t305 * t303 + t311 * (pkin(5) * t396 + t383) + t335 * (pkin(5) * t334 + t359)) * MDP(30) + 0.2e1 * t533 * t598 + (-t444 * MDP(10) - t445 * MDP(11) - t382 * MDP(13) - MDP(14) * t381) * qJD(3); t471 * MDP(14) + (t510 + t569) * MDP(20) + (-t437 ^ 2 * t491 - t555 - t567) * MDP(21) + (t513 - t573) * MDP(27) + (-t572 - t582) * MDP(28) + (t514 + t581) * MDP(29) + (t301 * t507 + t302 * t458 + t305 * t548 - t307 * t547 - t335 * t443) * MDP(30) + ((t486 * t541 + t487 * t542 + t443) * MDP(13) + (t442 - t526) * MDP(14)) * qJD(3) - qJD(1) ^ 2 * t598; -t442 ^ 2 * MDP(9) + (t471 + (-t442 - t526) * qJD(3)) * MDP(10) + (qJD(3) * t408 - t370) * MDP(13) + (-t442 * t464 + t453 * t533) * MDP(14) + (t420 * t516 + t571) * MDP(15) + ((t376 + t570) * t491 + (-t530 - t568) * t489) * MDP(16) + (t437 * t516 + t555 - t567) * MDP(17) + (t510 - t569) * MDP(18) + (-pkin(3) * t530 - t370 * t491 + t408 * t419 + (-pkin(8) * t537 - t394) * t437 + (t437 * t587 + t504) * t489) * MDP(20) + (-pkin(3) * t376 + t370 * t489 - t408 * t420 + (pkin(8) * t538 + t546) * t437 + t504 * t491) * MDP(21) + (-t328 * t458 - t508 * t547) * MDP(22) + (t514 - t581) * MDP(23) + (-t572 + t582) * MDP(24) + (t513 + t573) * MDP(25) + (t483 * t329 - t348 * t507 - t360 * t548 + t364 * t588 + t432 * t596 + t517 * t433) * MDP(27) + (-t483 * t328 + t348 * t458 - t547 * t360 + t432 * t586 - t544 * t433 + t508 * t588) * MDP(28) + (-t301 * t458 + t302 * t507 + t305 * t547 + t307 * t548 + t328 * t390 - t329 * t391 - t364 * t552 - t508 * t551) * MDP(29) + (t302 * t391 + t301 * t390 + t311 * (-pkin(5) * t507 + t483) + (-pkin(5) * t548 + t588) * t335 + t552 * t307 + t551 * t305) * MDP(30) + (-MDP(13) * t464 - MDP(19) * t437 - MDP(20) * t355 + MDP(21) * t356 - MDP(26) * t432 - MDP(27) * t316 + MDP(28) * t317 - MDP(8) * t442 + MDP(9) * t443) * t443; -t420 * t419 * MDP(15) + (-t419 ^ 2 + t420 ^ 2) * MDP(16) + (t376 - t570) * MDP(17) + (-t530 + t568) * MDP(18) + t433 * MDP(19) + (t356 * t437 - t402 * t420 + t498) * MDP(20) + (t355 * t437 - t402 * t419 - t503) * MDP(21) + (-t521 * t432 + (-t420 * t364 - t432 * t536 + t433 * t577) * pkin(4) + t583) * MDP(27) + (t553 * t432 + (-t420 * t508 - t432 * t522 - t488 * t433) * pkin(4) + t593) * MDP(28) + (-t305 * t364 + t307 * t508 + t308 * t508 + t309 * t364 + t482 * t328 + (-t329 * t488 + (-t364 * t577 + t488 * t508) * qJD(5)) * pkin(4)) * MDP(29) + (-pkin(5) * t590 + t301 * t482 - t305 * t308 - t307 * t309 + (t302 * t488 - t335 * t420 + (-t305 * t488 + t307 * t577) * qJD(5)) * pkin(4)) * MDP(30) + t600; (t317 * t432 + t583) * MDP(27) + (t316 * t432 + t593) * MDP(28) + (pkin(5) * t328 - t364 * t554) * MDP(29) + (t554 * t307 + (t301 - t590) * pkin(5)) * MDP(30) + t600; (-t361 - t579) * MDP(29) + (t305 * t508 + t307 * t364 + t311) * MDP(30);];
tauc  = t1;
