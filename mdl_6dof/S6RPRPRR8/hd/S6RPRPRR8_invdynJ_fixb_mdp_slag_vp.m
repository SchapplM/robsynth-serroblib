% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
% MDP [29x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRR8_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR8_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(29,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR8_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR8_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [29 1]), ...
  'S6RPRPRR8_invdynJ_fixb_mdp_slag_vp: MDP has to be [29x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:59:43
% EndTime: 2019-03-09 03:59:53
% DurationCPUTime: 6.91s
% Computational Cost: add. (4740->501), mult. (9720->647), div. (0->0), fcn. (7128->14), ass. (0->232)
t510 = sin(qJ(3));
t514 = cos(qJ(3));
t624 = sin(pkin(10));
t625 = cos(pkin(10));
t452 = -t624 * t510 + t625 * t514;
t447 = t452 * qJD(1);
t509 = sin(qJ(5));
t513 = cos(qJ(5));
t576 = t513 * qJD(3);
t417 = t447 * t509 - t576;
t525 = t510 * t625 + t514 * t624;
t634 = t525 * qJD(1);
t642 = qJD(5) + t634;
t648 = t417 * t642;
t516 = -pkin(1) - pkin(7);
t465 = qJDD(1) * t516 + qJDD(2);
t456 = t514 * t465;
t466 = qJD(1) * t516 + qJD(2);
t570 = qJDD(1) * t514;
t572 = qJD(1) * qJD(4);
t573 = qJD(1) * qJD(3);
t583 = qJD(3) * t510;
t376 = -t514 * t572 - t466 * t583 + qJDD(3) * pkin(3) + t456 + (t510 * t573 - t570) * qJ(4);
t582 = qJD(3) * t514;
t386 = (-qJ(4) * qJD(1) + t466) * t582 + (-qJ(4) * qJDD(1) + t465 - t572) * t510;
t347 = t376 * t625 - t624 * t386;
t343 = -qJDD(3) * pkin(4) - t347;
t478 = pkin(3) * t624 + pkin(8);
t500 = qJ(3) + pkin(10);
t486 = sin(t500);
t487 = cos(t500);
t511 = sin(qJ(1));
t515 = cos(qJ(1));
t635 = g(1) * t511 - g(2) * t515;
t524 = -g(3) * t486 + t487 * t635;
t647 = qJD(5) * t478 * t642 + t343 + t524;
t512 = cos(qJ(6));
t419 = qJD(3) * t509 + t447 * t513;
t508 = sin(qJ(6));
t613 = t419 * t508;
t362 = t512 * t417 + t613;
t432 = qJD(6) + t642;
t646 = t362 * t432;
t537 = t417 * t508 - t512 * t419;
t645 = t432 * t537;
t581 = qJD(5) * t509;
t611 = t634 * t509;
t644 = t581 + t611;
t518 = qJD(1) ^ 2;
t526 = -qJ(2) * t518 - t635;
t551 = t642 * t513;
t552 = qJDD(1) * t624;
t553 = qJDD(1) * t625;
t407 = -qJD(3) * t447 - t510 * t553 - t514 * t552;
t405 = -qJDD(5) + t407;
t601 = t509 * t405;
t643 = t551 * t642 - t601;
t584 = qJD(1) * t514;
t436 = -qJ(4) * t584 + t514 * t466;
t429 = qJD(3) * pkin(3) + t436;
t585 = qJD(1) * t510;
t435 = -qJ(4) * t585 + t466 * t510;
t561 = t625 * t435;
t380 = t624 * t429 + t561;
t370 = qJD(3) * pkin(8) + t380;
t458 = pkin(3) * t585 + qJD(1) * qJ(2) + qJD(4);
t381 = pkin(4) * t634 - pkin(8) * t447 + t458;
t346 = t370 * t513 + t381 * t509;
t337 = -pkin(9) * t417 + t346;
t579 = qJD(6) * t508;
t334 = t337 * t579;
t422 = t624 * t435;
t379 = t429 * t625 - t422;
t369 = -qJD(3) * pkin(4) - t379;
t354 = t417 * pkin(5) + t369;
t506 = qJ(5) + qJ(6);
t493 = cos(t506);
t604 = t493 * t511;
t492 = sin(t506);
t605 = t492 * t515;
t426 = t486 * t604 + t605;
t603 = t493 * t515;
t606 = t492 * t511;
t428 = t486 * t603 - t606;
t628 = g(3) * t487;
t641 = g(1) * t426 - g(2) * t428 + t354 * t362 + t493 * t628 + t334;
t425 = -t486 * t606 + t603;
t427 = t486 * t605 + t604;
t348 = t624 * t376 + t625 * t386;
t344 = qJDD(3) * pkin(8) + t348;
t408 = qJD(3) * t634 + t510 * t552 - t514 * t553;
t501 = qJDD(1) * qJ(2);
t502 = qJD(1) * qJD(2);
t563 = t514 * t573;
t571 = qJDD(1) * t510;
t536 = qJDD(4) + t501 + t502 + (t563 + t571) * pkin(3);
t351 = -pkin(4) * t407 + pkin(8) * t408 + t536;
t350 = t513 * t351;
t358 = qJD(5) * t576 + t509 * qJDD(3) - t513 * t408 - t447 * t581;
t324 = -pkin(5) * t405 - pkin(9) * t358 - qJD(5) * t346 - t344 * t509 + t350;
t554 = -t513 * qJDD(3) - t408 * t509;
t359 = t419 * qJD(5) + t554;
t580 = qJD(5) * t513;
t529 = t513 * t344 + t509 * t351 - t370 * t581 + t381 * t580;
t325 = -pkin(9) * t359 + t529;
t556 = t512 * t324 - t508 * t325;
t640 = -g(1) * t425 - g(2) * t427 + t354 * t537 + t492 * t628 + t556;
t403 = -qJDD(6) + t405;
t639 = (-t362 ^ 2 + t537 ^ 2) * MDP(24) - t403 * MDP(27) - t362 * t537 * MDP(23);
t602 = t508 * t513;
t455 = t509 * t512 + t602;
t398 = t455 * t452;
t633 = qJD(5) + qJD(6);
t411 = t633 * t455;
t454 = t508 * t509 - t512 * t513;
t636 = t454 * t403 - t411 * t432;
t543 = g(1) * t515 + g(2) * t511;
t568 = 0.2e1 * t502;
t632 = 0.2e1 * t501 + t568 - t543;
t410 = t633 * t454;
t592 = -t454 * t634 - t410;
t616 = t403 * t455;
t631 = -t432 * t592 + t616;
t555 = t358 * t508 + t512 * t359;
t331 = -qJD(6) * t537 + t555;
t627 = g(3) * t510;
t496 = t510 * pkin(3);
t626 = pkin(9) + t478;
t623 = pkin(1) * qJDD(1);
t345 = -t370 * t509 + t513 * t381;
t336 = -pkin(9) * t419 + t345;
t333 = pkin(5) * t642 + t336;
t621 = t333 * t512;
t620 = t337 * t512;
t619 = t358 * t509;
t618 = t362 * t447;
t617 = t537 * t447;
t615 = t417 * t447;
t614 = t419 * t447;
t557 = qJD(3) * t624;
t558 = qJD(3) * t625;
t445 = t510 * t557 - t514 * t558;
t612 = t432 * t445;
t446 = -t510 * t558 - t514 * t557;
t610 = t446 * t509;
t609 = t446 * t513;
t608 = t452 * t509;
t607 = t452 * t513;
t600 = t509 * t511;
t599 = t509 * t515;
t598 = t511 * t513;
t397 = t513 * t405;
t595 = qJ(4) - t516;
t459 = t595 * t510;
t460 = t595 * t514;
t413 = -t459 * t625 - t460 * t624;
t409 = t513 * t413;
t597 = t513 * t515;
t596 = qJ(2) + t496;
t392 = t436 * t625 - t422;
t394 = pkin(3) * t584 + pkin(4) * t447 + pkin(8) * t634;
t594 = t513 * t392 + t509 * t394;
t406 = pkin(4) * t525 - pkin(8) * t452 + t596;
t593 = t509 * t406 + t409;
t389 = t455 * t634;
t591 = t411 + t389;
t590 = t515 * pkin(1) + t511 * qJ(2);
t505 = t514 ^ 2;
t588 = t510 ^ 2 - t505;
t517 = qJD(3) ^ 2;
t587 = -t517 - t518;
t586 = qJD(1) * t458;
t578 = qJD(6) * t512;
t575 = pkin(3) * t582 + qJD(2);
t569 = qJDD(3) * t510;
t567 = t512 * t358 - t508 * t359 - t417 * t578;
t566 = t452 * t581;
t565 = t452 * t580;
t562 = qJD(5) * t626;
t550 = -qJD(5) * t381 - t344;
t549 = qJD(6) * t333 + t325;
t548 = qJDD(2) - t623;
t547 = -qJD(5) * t525 - qJD(1);
t479 = -pkin(3) * t625 - pkin(4);
t391 = t436 * t624 + t561;
t546 = pkin(5) * t644 - t391;
t541 = -t389 * t432 + t636;
t540 = -t370 * t580 + t350;
t388 = t513 * t394;
t451 = t626 * t513;
t539 = pkin(5) * t447 + qJD(6) * t451 - t392 * t509 + t388 + (pkin(9) * t634 + t562) * t513;
t450 = t626 * t509;
t538 = pkin(9) * t611 + qJD(6) * t450 + t509 * t562 + t594;
t433 = -qJD(4) * t514 + t583 * t595;
t434 = -qJD(3) * t460 - qJD(4) * t510;
t384 = -t625 * t433 + t434 * t624;
t412 = -t459 * t624 + t625 * t460;
t327 = t333 * t508 + t620;
t535 = -t642 * t644 - t397;
t534 = t565 + t610;
t533 = -t566 + t609;
t531 = t432 * t455;
t530 = 0.2e1 * qJ(2) * t573 + qJDD(3) * t516;
t385 = t433 * t624 + t434 * t625;
t393 = -pkin(4) * t445 - pkin(8) * t446 + t575;
t528 = t513 * t385 + t509 * t393 + t406 * t580 - t413 * t581;
t330 = -t419 * t579 + t567;
t527 = t369 * t642 + t478 * t405;
t520 = t347 * t452 + t348 * t525 + t379 * t446 - t380 * t445 - t635;
t519 = -t516 * t517 + t632;
t507 = -qJ(4) - pkin(7);
t495 = t515 * qJ(2);
t490 = qJDD(3) * t514;
t463 = -t513 * pkin(5) + t479;
t441 = t486 * t597 - t600;
t440 = t486 * t599 + t598;
t439 = t486 * t598 + t599;
t438 = -t486 * t600 + t597;
t401 = t513 * t406;
t399 = t454 * t452;
t383 = t513 * t393;
t377 = pkin(5) * t608 + t412;
t353 = pkin(5) * t534 + t384;
t352 = -pkin(9) * t608 + t593;
t342 = pkin(5) * t525 - pkin(9) * t607 - t413 * t509 + t401;
t339 = t446 * t602 - t508 * t566 - t579 * t608 + (t607 * t633 + t610) * t512;
t338 = -t398 * t633 - t454 * t446;
t332 = t359 * pkin(5) + t343;
t329 = -pkin(9) * t534 + t528;
t328 = -pkin(9) * t609 - pkin(5) * t445 - t385 * t509 + t383 + (-t409 + (pkin(9) * t452 - t406) * t509) * qJD(5);
t326 = -t337 * t508 + t621;
t1 = [t632 * MDP(5) + (qJDD(1) * t505 - 0.2e1 * t510 * t563) * MDP(7) + t543 * MDP(3) + (t510 * t519 + t514 * t530) * MDP(12) + (-t510 * t530 + t514 * t519) * MDP(13) + qJDD(1) * MDP(1) + (-t330 * t399 - t338 * t537) * MDP(23) + (-t330 * t398 + t331 * t399 - t338 * t362 + t339 * t537) * MDP(24) + (t384 * t447 - t385 * t634 + t407 * t413 - t408 * t412 - t520) * MDP(14) + (g(1) * t427 - g(2) * t425 + t327 * t445 + t377 * t330 - t332 * t399 + t334 * t525 + t354 * t338 - t353 * t537 + (-(-qJD(6) * t352 + t328) * t432 + t342 * t403 - t324 * t525) * t508 + (-(qJD(6) * t342 + t329) * t432 + t352 * t403 - t549 * t525) * t512) * MDP(29) + (t330 * t525 + t338 * t432 + t399 * t403 + t445 * t537) * MDP(25) + (-t403 * t525 - t612) * MDP(27) + (-t331 * t525 - t339 * t432 + t362 * t445 + t398 * t403) * MDP(26) + ((t328 * t512 - t329 * t508) * t432 - (t342 * t512 - t352 * t508) * t403 + t556 * t525 - t326 * t445 + t353 * t362 + t377 * t331 + t332 * t398 + t354 * t339 - g(1) * t428 - g(2) * t426 + ((-t342 * t508 - t352 * t512) * t432 - t327 * t525) * qJD(6)) * MDP(28) + t635 * MDP(2) + (qJDD(2) - t635 - 0.2e1 * t623) * MDP(4) + (-t405 * t525 - t445 * t642) * MDP(20) + ((-t413 * t580 + t383) * t642 - t401 * t405 + t540 * t525 - t345 * t445 + t384 * t417 + t412 * t359 + t369 * t565 - g(1) * t441 - g(2) * t439 + ((-qJD(5) * t406 - t385) * t642 + t413 * t405 + t550 * t525 + t343 * t452 + t369 * t446) * t509) * MDP(21) + (t358 * t525 - t397 * t452 - t419 * t445 + t533 * t642) * MDP(18) + (-t359 * t525 + t417 * t445 + t452 * t601 - t534 * t642) * MDP(19) + (g(1) * t440 - g(2) * t438 + t343 * t607 + t346 * t445 + t412 * t358 + t369 * t533 + t384 * t419 + t405 * t593 - t525 * t529 - t528 * t642) * MDP(22) + (-t514 * t517 - t569) * MDP(10) + 0.2e1 * (-t510 * t570 + t573 * t588) * MDP(8) + (-t548 * pkin(1) - g(1) * (-pkin(1) * t511 + t495) - g(2) * t590 + (t568 + t501) * qJ(2)) * MDP(6) + (t358 * t607 + t419 * t533) * MDP(16) + ((-t417 * t513 - t419 * t509) * t446 + (-t619 - t359 * t513 + (t417 * t509 - t419 * t513) * qJD(5)) * t452) * MDP(17) + (t348 * t413 + t380 * t385 - t347 * t412 - t379 * t384 + t536 * t596 + t458 * t575 - g(1) * (t515 * t496 + t495 + (-pkin(1) + t507) * t511) - g(2) * (t496 * t511 - t507 * t515 + t590)) * MDP(15) + (-t510 * t517 + t490) * MDP(9); qJDD(1) * MDP(4) - t518 * MDP(5) + (t548 + t526) * MDP(6) + (t510 * t587 + t490) * MDP(12) + (t514 * t587 - t569) * MDP(13) + (t407 * t525 + t408 * t452 + t445 * t634 - t446 * t447) * MDP(14) + (t520 - t586) * MDP(15) + (t525 * t601 - t359 * t452 - t417 * t446 + (t445 * t509 + t513 * t547) * t642) * MDP(21) + (t525 * t397 - t358 * t452 - t419 * t446 + (t445 * t513 - t509 * t547) * t642) * MDP(22) + (-t452 * t331 - t446 * t362 + t445 * t531 + t454 * t432 * qJD(1) - (-t410 * t432 - t616) * t525) * MDP(28) + (qJD(1) * t531 - t452 * t330 + t446 * t537 - t454 * t612 - t525 * t636) * MDP(29); MDP(9) * t570 - MDP(10) * t571 + qJDD(3) * MDP(11) + (t514 * t526 + t456 + t627) * MDP(12) + (g(3) * t514 + (-t465 - t526) * t510) * MDP(13) + (-(t379 - t392) * t634 + (t407 * t624 + t408 * t625) * pkin(3)) * MDP(14) + (t379 * t391 - t380 * t392 + (t624 * t348 + t625 * t347 + t627 + (-t635 - t586) * t514) * pkin(3)) * MDP(15) + (t419 * t551 + t619) * MDP(16) + ((t358 - t648) * t513 + (-t419 * t642 - t359) * t509) * MDP(17) + (-t614 + t643) * MDP(18) + (t535 + t615) * MDP(19) + (t479 * t359 - t388 * t642 - t391 * t417 + (t392 * t642 + t527) * t509 - t647 * t513) * MDP(21) + (t479 * t358 - t391 * t419 + t509 * t647 + t527 * t513 + t594 * t642) * MDP(22) + (t330 * t455 - t537 * t592) * MDP(23) + (-t330 * t454 - t331 * t455 - t362 * t592 + t537 * t591) * MDP(24) + (t617 - t631) * MDP(25) + (t541 + t618) * MDP(26) + (-(-t450 * t512 - t451 * t508) * t403 + t463 * t331 + t332 * t454 + (t508 * t538 - t512 * t539) * t432 + t546 * t362 + t591 * t354 - t524 * t493) * MDP(28) + ((-t450 * t508 + t451 * t512) * t403 + t463 * t330 + t332 * t455 + (t508 * t539 + t512 * t538) * t432 - t546 * t537 + t592 * t354 + t524 * t492) * MDP(29) + (MDP(7) * t510 * t514 - MDP(8) * t588) * t518 + ((t380 - t391) * MDP(14) - t642 * MDP(20) - t345 * MDP(21) + t346 * MDP(22) - t432 * MDP(27) - t326 * MDP(28) + t327 * MDP(29)) * t447; (-t447 ^ 2 - t634 ^ 2) * MDP(14) + (t379 * t447 + t380 * t634 + t536 - t543) * MDP(15) + (t535 - t615) * MDP(21) + (-t614 - t643) * MDP(22) + (t541 - t618) * MDP(28) + (t617 + t631) * MDP(29); t419 * t417 * MDP(16) + (-t417 ^ 2 + t419 ^ 2) * MDP(17) + (t358 + t648) * MDP(18) + (-t554 + (-qJD(5) + t642) * t419) * MDP(19) - t405 * MDP(20) + (-g(1) * t438 - g(2) * t440 + t346 * t642 - t369 * t419 + (t550 + t628) * t509 + t540) * MDP(21) + (g(1) * t439 - g(2) * t441 + t345 * t642 + t369 * t417 + t513 * t628 - t529) * MDP(22) + (t330 + t646) * MDP(25) + (-t331 - t645) * MDP(26) + (-(-t336 * t508 - t620) * t432 - t327 * qJD(6) + (-t362 * t419 - t512 * t403 - t432 * t579) * pkin(5) + t640) * MDP(28) + ((-t337 * t432 - t324) * t508 + (t336 * t432 - t549) * t512 + (t508 * t403 + t419 * t537 - t432 * t578) * pkin(5) + t641) * MDP(29) + t639; (t567 + t646) * MDP(25) + (-t555 - t645) * MDP(26) + (t327 * t432 + t640) * MDP(28) + (-t508 * t324 - t512 * t325 + t326 * t432 + t641) * MDP(29) + (-MDP(25) * t613 + MDP(26) * t537 - MDP(28) * t327 - MDP(29) * t621) * qJD(6) + t639;];
tau  = t1;
