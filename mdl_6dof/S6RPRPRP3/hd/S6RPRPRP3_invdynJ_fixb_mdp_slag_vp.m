% Calculate vector of inverse dynamics joint torques for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% MDP [26x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRPRP3_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRP3_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(26,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [26 1]), ...
  'S6RPRPRP3_invdynJ_fixb_mdp_slag_vp: MDP has to be [26x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:50
% EndTime: 2019-03-09 03:09:59
% DurationCPUTime: 8.21s
% Computational Cost: add. (5308->524), mult. (10999->660), div. (0->0), fcn. (7484->14), ass. (0->218)
t620 = qJDD(3) * pkin(3);
t504 = sin(pkin(9));
t487 = pkin(1) * t504 + pkin(7);
t469 = t487 * qJD(1);
t509 = sin(qJ(3));
t511 = cos(qJ(3));
t586 = qJD(3) * t511;
t467 = t487 * qJDD(1);
t641 = -qJD(2) * qJD(3) - t467;
t545 = t511 * qJDD(2) - t469 * t586 + t509 * t641;
t380 = qJDD(4) - t545 - t620;
t498 = g(3) * t511;
t500 = qJ(1) + pkin(9);
t494 = cos(t500);
t492 = sin(t500);
t625 = g(2) * t492;
t557 = g(1) * t494 + t625;
t527 = t557 * t509 - t498;
t642 = t380 - t527;
t649 = -t642 + t620;
t503 = sin(pkin(10));
t622 = pkin(8) + qJ(4);
t465 = t622 * t503;
t505 = cos(pkin(10));
t466 = t622 * t505;
t508 = sin(qJ(5));
t631 = cos(qJ(5));
t407 = -t508 * t465 + t466 * t631;
t495 = t511 * qJDD(1);
t583 = qJD(1) * qJD(3);
t634 = t509 * t583 - t495;
t450 = qJDD(5) + t634;
t499 = pkin(10) + qJ(5);
t491 = sin(t499);
t648 = t407 * t450 + t491 * t527;
t647 = MDP(21) + MDP(23);
t646 = -MDP(22) + MDP(25);
t588 = qJD(3) * t505;
t591 = qJD(1) * t509;
t447 = -t503 * t591 + t588;
t575 = t505 * t591;
t448 = qJD(3) * t503 + t575;
t397 = -t631 * t447 + t448 * t508;
t645 = t397 ^ 2;
t539 = -t508 * t447 - t448 * t631;
t632 = t539 ^ 2;
t590 = qJD(1) * t511;
t481 = -qJD(5) + t590;
t644 = t397 * t481;
t643 = t481 * t539;
t537 = -t508 * t503 + t505 * t631;
t570 = qJD(5) * t631;
t585 = qJD(5) * t508;
t635 = -t503 * t585 + t505 * t570;
t595 = -t537 * t590 + t635;
t452 = t503 * t631 + t508 * t505;
t439 = t452 * qJD(5);
t531 = t511 * t452;
t594 = -qJD(1) * t531 + t439;
t572 = t505 * t586;
t573 = t503 * t586;
t385 = t439 * t509 + t508 * t573 - t572 * t631;
t429 = t537 * t509;
t600 = t385 * t481 + t429 * t450;
t386 = qJD(3) * t531 + t509 * t635;
t428 = t452 * t509;
t599 = t386 * t481 - t428 * t450;
t432 = qJD(2) * t511 - t509 * t469;
t553 = pkin(3) * t509 - qJ(4) * t511;
t457 = t553 * qJD(1);
t390 = -t432 * t503 + t505 * t457;
t606 = t505 * t511;
t548 = pkin(4) * t509 - pkin(8) * t606;
t361 = qJD(1) * t548 + t390;
t391 = t505 * t432 + t503 * t457;
t576 = t503 * t590;
t372 = -pkin(8) * t576 + t391;
t538 = -t465 * t631 - t508 * t466;
t638 = qJD(4) * t537 + qJD(5) * t538 - t508 * t361 - t631 * t372;
t637 = qJD(4) * t452 + qJD(5) * t407 + t361 * t631 - t508 * t372;
t554 = pkin(3) * t511 + qJ(4) * t509;
t547 = -pkin(2) - t554;
t506 = cos(pkin(9));
t630 = pkin(1) * t506;
t442 = t547 - t630;
t431 = t505 * t442;
t607 = t505 * t509;
t382 = -pkin(8) * t607 + t431 + (-t487 * t503 - pkin(4)) * t511;
t405 = t503 * t442 + t487 * t606;
t609 = t503 * t509;
t389 = -pkin(8) * t609 + t405;
t636 = t508 * t382 + t631 * t389;
t568 = t511 * t583;
t581 = qJDD(1) * t509;
t536 = t568 + t581;
t580 = qJDD(3) * t503;
t520 = t505 * t536 + t580;
t477 = t503 * t581;
t593 = t503 * t568 + t477;
t552 = qJDD(3) * t505 - t593;
t348 = -t447 * t570 + t448 * t585 - t508 * t552 - t631 * t520;
t587 = qJD(3) * t509;
t565 = t348 * t511 - t539 * t587;
t349 = -qJD(5) * t539 + t508 * t520 - t631 * t552;
t564 = -t349 * t511 + t397 * t587;
t436 = qJD(3) * t553 - qJD(4) * t509;
t574 = t487 * t587;
t394 = t505 * t436 + t503 * t574;
t370 = qJD(3) * t548 + t394;
t425 = t503 * t436;
t471 = t509 * t487;
t608 = t503 * t511;
t381 = t425 + (-pkin(8) * t608 - t471 * t505) * qJD(3);
t633 = -qJD(5) * t636 + t370 * t631 - t508 * t381;
t502 = t511 ^ 2;
t629 = pkin(4) * t503;
t628 = pkin(5) * t450;
t627 = g(1) * t492;
t624 = g(2) * t494;
t623 = g(3) * t509;
t621 = qJ(6) * t450;
t433 = t509 * qJD(2) + t511 * t469;
t421 = qJD(3) * qJ(4) + t433;
t424 = t442 * qJD(1);
t364 = -t421 * t503 + t505 * t424;
t355 = -pkin(4) * t590 - pkin(8) * t448 + t364;
t365 = t505 * t421 + t503 * t424;
t358 = pkin(8) * t447 + t365;
t331 = t508 * t355 + t358 * t631;
t619 = t331 * t481;
t616 = t397 * t539;
t614 = t491 * t492;
t493 = cos(t499);
t613 = t493 * t509;
t612 = t493 * t511;
t611 = t494 * t509;
t610 = t494 * t511;
t605 = t622 * t509;
t604 = qJDD(2) - g(3);
t603 = -t429 * t349 + t385 * t397;
t408 = pkin(4) * t576 + t433;
t602 = pkin(5) * t594 - qJ(6) * t595 - qJD(6) * t452 - t408;
t374 = qJDD(3) * qJ(4) + qJDD(2) * t509 + t467 * t511 + (qJD(4) + t432) * qJD(3);
t383 = qJD(1) * t436 + qJDD(1) * t442;
t345 = t505 * t374 + t503 * t383;
t598 = -qJ(6) * t591 + t638;
t597 = pkin(5) * t591 + t637;
t427 = pkin(4) * t573 + t487 * t586;
t435 = pkin(4) * t609 + t471;
t501 = t509 ^ 2;
t592 = t501 - t502;
t489 = -pkin(2) - t630;
t470 = qJD(1) * t489;
t330 = t355 * t631 - t508 * t358;
t584 = qJD(6) - t330;
t579 = qJDD(3) * t511;
t512 = cos(qJ(1));
t578 = t512 * pkin(1) + t494 * pkin(2) + t492 * pkin(7);
t488 = pkin(4) * t505 + pkin(3);
t577 = qJ(4) * t587;
t510 = sin(qJ(1));
t567 = -pkin(1) * t510 + t494 * pkin(7);
t566 = -qJD(3) * pkin(3) + qJD(4);
t344 = -t503 * t374 + t505 * t383;
t563 = qJD(1) * t592;
t334 = pkin(4) * t634 - t520 * pkin(8) + t344;
t339 = pkin(8) * t552 + t345;
t562 = -t631 * t334 + t508 * t339 + t355 * t585 + t358 * t570;
t474 = t509 * t627;
t560 = -g(2) * t611 + t474;
t415 = t493 * t494 + t511 * t614;
t417 = t491 * t610 - t492 * t493;
t559 = -g(1) * t415 + g(2) * t417;
t416 = -t494 * t491 + t492 * t612;
t418 = t493 * t610 + t614;
t558 = g(1) * t416 - g(2) * t418;
t556 = -t624 + t627;
t555 = g(1) * t510 - g(2) * t512;
t551 = -t344 * t503 + t345 * t505;
t550 = -t348 * t428 - t386 * t539;
t549 = -t364 * t503 + t365 * t505;
t546 = pkin(5) * t493 + qJ(6) * t491 + t488;
t542 = t382 * t631 - t508 * t389;
t535 = t505 * t581 + t580;
t534 = t508 * t334 + t631 * t339 + t355 * t570 - t358 * t585;
t533 = t508 * t370 + t631 * t381 + t382 * t570 - t389 * t585;
t532 = -qJD(1) * t470 + t557;
t419 = -t432 + t566;
t530 = g(1) * t493 * t611 - g(3) * t612 + t450 * t538 + t613 * t625;
t513 = qJD(3) ^ 2;
t529 = 0.2e1 * qJDD(1) * t489 + t487 * t513 + t624;
t528 = 0.2e1 * qJD(3) * t470 - qJDD(3) * t487;
t525 = -t511 * t557 - t623;
t392 = -pkin(4) * t447 + t419;
t523 = g(1) * t417 + g(2) * t415 + t491 * t623 - t562;
t342 = pkin(5) * t397 + qJ(6) * t539 + t392;
t519 = -t342 * t539 + qJDD(6) - t523;
t518 = -g(1) * t418 - g(2) * t416 - g(3) * t613 + t534;
t357 = -pkin(4) * t552 + t380;
t325 = t349 * pkin(5) + t348 * qJ(6) + qJD(6) * t539 + t357;
t462 = -t509 * t513 + t579;
t461 = qJDD(3) * t509 + t511 * t513;
t404 = -t487 * t608 + t431;
t395 = -t505 * t574 + t425;
t393 = -pkin(5) * t537 - qJ(6) * t452 - t488;
t359 = pkin(5) * t428 - qJ(6) * t429 + t435;
t354 = -pkin(5) * t539 + qJ(6) * t397;
t347 = t511 * pkin(5) - t542;
t346 = -qJ(6) * t511 + t636;
t336 = pkin(5) * t386 + qJ(6) * t385 - qJD(6) * t429 + t427;
t335 = -t348 - t644;
t329 = -t481 * qJ(6) + t331;
t328 = t481 * pkin(5) + t584;
t327 = -pkin(5) * t587 - t633;
t326 = qJ(6) * t587 - qJD(6) * t511 + t533;
t324 = qJDD(6) + t562 - t628;
t323 = -qJD(6) * t481 + t534 + t621;
t1 = [(t599 - t564) * MDP(19) + (t330 * t587 + t435 * t349 + t357 * t428 + t392 * t386 + t427 * t397 + t542 * t450 - t481 * t633 + t562 * t511 + t558) * MDP(21) + (g(1) * t512 + g(2) * t510) * MDP(3) + t461 * MDP(7) + t462 * MDP(8) + (-t348 * t429 + t385 * t539) * MDP(16) + (-t323 * t511 - t325 * t429 - t326 * t481 + t329 * t587 + t336 * t539 + t342 * t385 + t346 * t450 + t348 * t359 - t559) * MDP(25) + (-t323 * t428 + t324 * t429 - t326 * t397 - t327 * t539 - t328 * t385 - t329 * t386 - t346 * t349 - t347 * t348 + t560) * MDP(24) + (-t331 * t587 - t435 * t348 + t357 * t429 - t392 * t385 - t427 * t539 - t450 * t636 + t481 * t533 + t511 * t534 + t559) * MDP(22) + 0.2e1 * (-qJD(3) * t563 + t495 * t509) * MDP(6) + (t323 * t346 + t329 * t326 + t325 * t359 + t342 * t336 + t324 * t347 + t328 * t327 - g(1) * (-pkin(5) * t416 - qJ(6) * t415 + t494 * t629 + t567) - g(2) * (pkin(5) * t418 + qJ(6) * t417 + t488 * t610 + t494 * t605 + t578) + (-g(1) * (-t488 * t511 - pkin(2) - t605) - g(2) * t629) * t492) * MDP(26) + (-t557 * t505 + (t380 * t505 + t535 * t487 + (-qJD(1) * t405 - t365) * qJD(3)) * t509 + (t395 * qJD(1) + t405 * qJDD(1) + t345 - t556 * t503 + (t419 * t505 + (t448 + t575) * t487) * qJD(3)) * t511) * MDP(13) + (t565 + t600) * MDP(18) + (-t550 + t603) * MDP(17) + (t395 * t447 - t405 * t593 - t394 * t448 + (-t404 * qJDD(3) - t345 * t509 - t365 * t586) * t503 + (t405 * qJDD(3) - t344 * t509 - t364 * t586 - t404 * t536) * t505 + t560) * MDP(14) + (-t450 * t511 - t481 * t587) * MDP(20) + (t324 * t511 + t325 * t428 + t327 * t481 - t328 * t587 + t336 * t397 + t342 * t386 - t347 * t450 + t349 * t359 + t558) * MDP(23) + (qJDD(1) * t501 + 0.2e1 * t509 * t568) * MDP(5) + (t345 * t405 + t365 * t395 + t344 * t404 + t364 * t394 - g(1) * t567 - g(2) * (t494 * t554 + t578) - t547 * t627 + (t380 * t509 + t419 * t586) * t487) * MDP(15) + (t528 * t509 + (-t529 + t627) * t511) * MDP(10) + qJDD(1) * MDP(1) + (t509 * t529 + t511 * t528 - t474) * MDP(11) + (t555 + (t504 ^ 2 + t506 ^ 2) * qJDD(1) * pkin(1)) * pkin(1) * MDP(4) + t555 * MDP(2) + (-t557 * t503 + (-t487 * t552 + t380 * t503 + (t404 * qJD(1) + t364) * qJD(3)) * t509 + (-t394 * qJD(1) - t404 * qJDD(1) - t344 + t556 * t505 + (t419 * t503 - t447 * t487) * qJD(3)) * t511) * MDP(12); t604 * MDP(4) + t462 * MDP(10) - t461 * MDP(11) + ((t552 + t477) * t511 + (-t509 * t447 - t503 * t563) * qJD(3)) * MDP(12) + (-t503 * t579 + (t509 * t448 + (-t592 - t502) * t505 * qJD(1)) * qJD(3)) * MDP(13) + (t447 * t572 + t448 * t573 + t520 * t609 + t552 * t607) * MDP(14) + (-t380 * t511 - g(3) + t551 * t509 + (t419 * t509 + t511 * t549) * qJD(3)) * MDP(15) + (t550 + t603) * MDP(24) + (t323 * t429 + t324 * t428 - t325 * t511 + t328 * t386 - t329 * t385 + t342 * t587 - g(3)) * MDP(26) + t647 * (t564 + t599) + t646 * (t600 - t565); MDP(7) * t581 + MDP(8) * t495 + qJDD(3) * MDP(9) + (qJD(3) * t433 + t509 * t532 - t498 + t545) * MDP(10) + (qJD(3) * t432 + (qJD(3) * t469 - t604) * t509 + (t532 + t641) * t511) * MDP(11) + (t503 * qJ(4) * t495 - pkin(3) * t593 + t433 * t447 + t649 * t505 + (-t364 * t509 + t390 * t511 + (-t577 + (qJD(4) - t419) * t511) * t503) * qJD(1)) * MDP(12) + (-t433 * t448 - t553 * t505 * qJDD(1) - t649 * t503 + (t365 * t509 - t391 * t511 + (-t577 + (-t419 + t566) * t511) * t505) * qJD(1)) * MDP(13) + (t390 * t448 - t391 * t447 + (qJ(4) * t552 + qJD(4) * t447 + t364 * t590 + t345) * t505 + (qJ(4) * t520 + qJD(4) * t448 + t365 * t590 - t344) * t503 + t525) * MDP(14) + (-t364 * t390 - t365 * t391 - t419 * t433 + t549 * qJD(4) - t642 * pkin(3) + (t525 + t551) * qJ(4)) * MDP(15) + (-t348 * t452 - t539 * t595) * MDP(16) + (-t348 * t537 - t349 * t452 - t397 * t595 + t539 * t594) * MDP(17) + (t450 * t452 - t481 * t595 + t539 * t591) * MDP(18) + (t397 * t591 + t450 * t537 + t481 * t594) * MDP(19) + t481 * MDP(20) * t591 + (-t330 * t591 - t488 * t349 - t357 * t537 + t594 * t392 - t408 * t397 + t481 * t637 + t530) * MDP(21) + (t331 * t591 + t488 * t348 + t357 * t452 + t595 * t392 + t408 * t539 + t481 * t638 - t648) * MDP(22) + (-t325 * t537 + t328 * t591 + t342 * t594 + t349 * t393 + t397 * t602 + t481 * t597 + t530) * MDP(23) + (t323 * t537 + t324 * t452 + t328 * t595 - t329 * t594 + t348 * t538 - t349 * t407 - t397 * t598 - t539 * t597 + t525) * MDP(24) + (-t325 * t452 - t329 * t591 - t342 * t595 + t348 * t393 - t481 * t598 + t539 * t602 + t648) * MDP(25) + (-g(3) * t605 + t323 * t407 - t324 * t538 + t325 * t393 + t328 * t597 + t329 * t598 + t342 * t602 - t498 * t546 + t557 * (t509 * t546 - t511 * t622)) * MDP(26) + (-MDP(5) * t509 * t511 + MDP(6) * t592) * qJD(1) ^ 2; (-t448 * t590 - t552) * MDP(12) + ((-t447 + t588) * t590 + t535) * MDP(13) + (-t447 ^ 2 - t448 ^ 2) * MDP(14) + (t364 * t448 - t365 * t447 + t642) * MDP(15) + (-t632 - t645) * MDP(24) + (t328 * t539 + t329 * t397 + t325 - t527) * MDP(26) + t646 * (t348 - t644) + t647 * (t349 + t643); -MDP(16) * t616 + (t632 - t645) * MDP(17) + t335 * MDP(18) + (-t349 + t643) * MDP(19) + t450 * MDP(20) + (t392 * t539 + t523 - t619) * MDP(21) + (-t330 * t481 + t392 * t397 - t518) * MDP(22) + (-t354 * t397 - t519 - t619 + 0.2e1 * t628) * MDP(23) + (pkin(5) * t348 - qJ(6) * t349 - (t329 - t331) * t539 + (t328 - t584) * t397) * MDP(24) + (0.2e1 * t621 - t342 * t397 - t354 * t539 + (-0.2e1 * qJD(6) + t330) * t481 + t518) * MDP(25) + (t323 * qJ(6) - t324 * pkin(5) - t342 * t354 - t328 * t331 - g(1) * (-pkin(5) * t417 + qJ(6) * t418) - g(2) * (-pkin(5) * t415 + qJ(6) * t416) - (-pkin(5) * t491 + qJ(6) * t493) * t623 + t584 * t329) * MDP(26); (-t450 - t616) * MDP(23) + t335 * MDP(24) + (-t481 ^ 2 - t632) * MDP(25) + (t329 * t481 + t519 - t628) * MDP(26);];
tau  = t1;
