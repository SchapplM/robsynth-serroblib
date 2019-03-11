% Calculate vector of inverse dynamics joint torques for
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% MDP [25x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRPP1_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPRPP1_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(25,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [25 1]), ...
  'S6RRPRPP1_invdynJ_fixb_mdp_slag_vp: MDP has to be [25x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:48:07
% EndTime: 2019-03-09 09:48:18
% DurationCPUTime: 8.28s
% Computational Cost: add. (9849->576), mult. (22809->717), div. (0->0), fcn. (16735->14), ass. (0->243)
t570 = sin(pkin(9));
t572 = cos(pkin(9));
t576 = sin(qJ(2));
t579 = cos(qJ(2));
t525 = t570 * t579 + t572 * t576;
t569 = sin(pkin(10));
t571 = cos(pkin(10));
t575 = sin(qJ(4));
t578 = cos(qJ(4));
t691 = -t569 * t575 + t571 * t578;
t453 = t691 * t525;
t513 = t525 * qJD(1);
t636 = t578 * qJD(2);
t484 = t513 * t575 - t636;
t486 = qJD(2) * t575 + t513 * t578;
t429 = t571 * t484 + t486 * t569;
t639 = qJD(1) * t576;
t659 = t572 * t579;
t510 = qJD(1) * t659 - t570 * t639;
t501 = qJD(4) - t510;
t696 = t429 * t501;
t524 = t569 * t578 + t571 * t575;
t509 = t524 * qJD(4);
t645 = t524 * t510 - t509;
t637 = qJD(4) * t578;
t638 = qJD(4) * t575;
t644 = t691 * t510 + t569 * t638 - t571 * t637;
t523 = t570 * t576 - t659;
t515 = t523 * qJD(2);
t625 = t525 * t637;
t695 = -t515 * t575 + t625;
t602 = -t484 * t569 + t571 * t486;
t694 = t602 ^ 2;
t563 = t579 * pkin(2);
t555 = t563 + pkin(1);
t614 = t501 * t578;
t566 = qJ(2) + pkin(9);
t557 = sin(t566);
t577 = sin(qJ(1));
t580 = cos(qJ(1));
t611 = g(1) * t580 + g(2) * t577;
t693 = t557 * t611;
t448 = pkin(2) * t639 + pkin(3) * t513 - pkin(8) * t510;
t444 = t578 * t448;
t679 = qJ(3) + pkin(7);
t533 = t679 * t579;
t528 = qJD(1) * t533;
t516 = t570 * t528;
t532 = t679 * t576;
t527 = qJD(1) * t532;
t477 = -t527 * t572 - t516;
t401 = -qJ(5) * t510 * t578 + pkin(4) * t513 - t477 * t575 + t444;
t642 = t575 * t448 + t578 * t477;
t670 = t510 * t575;
t408 = -qJ(5) * t670 + t642;
t550 = pkin(2) * t570 + pkin(8);
t650 = qJ(5) + t550;
t615 = qJD(4) * t650;
t489 = qJD(5) * t578 - t575 * t615;
t591 = -qJD(5) * t575 - t578 * t615;
t648 = (-t401 + t591) * t571 + (t408 - t489) * t569;
t620 = g(1) * t577 - g(2) * t580;
t660 = t572 * t528;
t476 = -t527 * t570 + t660;
t692 = -t476 + (t638 - t670) * pkin(4);
t690 = MDP(20) + MDP(23);
t619 = qJD(2) * t679;
t508 = -qJD(3) * t576 - t579 * t619;
t467 = qJDD(2) * pkin(2) + qJD(1) * t508 - qJDD(1) * t532;
t507 = qJD(3) * t579 - t576 * t619;
t475 = qJD(1) * t507 + qJDD(1) * t533;
t421 = t467 * t572 - t570 * t475;
t419 = -qJDD(2) * pkin(3) - t421;
t559 = cos(t566);
t547 = g(3) * t559;
t590 = -t547 + t693;
t689 = -qJD(4) * t550 * t501 - t419 + t590;
t678 = qJD(2) * pkin(2);
t519 = -t527 + t678;
t471 = t519 * t572 - t516;
t458 = -qJD(2) * pkin(3) - t471;
t427 = pkin(4) * t484 + qJD(5) + t458;
t388 = pkin(5) * t429 - qJ(6) * t602 + t427;
t565 = qJ(4) + pkin(10);
t558 = cos(t565);
t556 = sin(t565);
t666 = t556 * t577;
t490 = t558 * t580 + t559 * t666;
t651 = t580 * t556;
t654 = t577 * t558;
t492 = t559 * t651 - t654;
t632 = qJDD(1) * t579;
t633 = qJDD(1) * t576;
t604 = -t570 * t633 + t572 * t632;
t473 = -qJD(2) * t513 + t604;
t634 = qJD(1) * qJD(2);
t622 = t579 * t634;
t623 = t576 * t634;
t474 = qJDD(1) * t525 - t570 * t623 + t572 * t622;
t592 = pkin(2) * t623 - qJDD(1) * t555 + qJDD(3);
t414 = -pkin(3) * t473 - pkin(8) * t474 + t592;
t410 = t578 * t414;
t531 = -qJD(1) * t555 + qJD(3);
t438 = -pkin(3) * t510 - pkin(8) * t513 + t531;
t472 = t570 * t519 + t660;
t459 = qJD(2) * pkin(8) + t472;
t412 = t438 * t575 + t459 * t578;
t422 = t570 * t467 + t572 * t475;
t420 = qJDD(2) * pkin(8) + t422;
t425 = qJD(4) * t636 + t575 * qJDD(2) + t578 * t474 - t513 * t638;
t511 = t525 * qJD(2);
t468 = qJD(1) * t511 + qJDD(4) - t604;
t364 = pkin(4) * t468 - qJ(5) * t425 - qJD(4) * t412 - qJD(5) * t486 - t420 * t575 + t410;
t426 = qJD(4) * t486 - t578 * qJDD(2) + t474 * t575;
t596 = t575 * t414 + t578 * t420 + t438 * t637 - t459 * t638;
t367 = -qJ(5) * t426 - qJD(5) * t484 + t596;
t355 = t571 * t364 - t569 * t367;
t621 = -qJDD(6) + t355;
t681 = g(3) * t557;
t688 = g(1) * t492 + g(2) * t490 - t388 * t602 + t556 * t681 + t621;
t687 = t559 * t611 + t681;
t686 = pkin(2) * t572;
t685 = pkin(5) * t468;
t680 = g(3) * t579;
t403 = -qJ(5) * t484 + t412;
t399 = t571 * t403;
t411 = t578 * t438 - t459 * t575;
t402 = -qJ(5) * t486 + t411;
t375 = t402 * t569 + t399;
t677 = t375 * t602;
t676 = t403 * t569;
t675 = t425 * t575;
t674 = t484 * t501;
t673 = t484 * t513;
t672 = t486 * t501;
t671 = t486 * t513;
t668 = t525 * t575;
t667 = t525 * t578;
t665 = t557 * t577;
t664 = t557 * t580;
t554 = pkin(4) * t578 + pkin(3);
t529 = t559 * t554;
t663 = t559 * t580;
t658 = t679 * t580;
t657 = t575 * t468;
t656 = t575 * t577;
t655 = t575 * t580;
t653 = t577 * t578;
t455 = t578 * t468;
t482 = -t532 * t570 + t533 * t572;
t478 = t578 * t482;
t652 = t578 * t580;
t356 = t569 * t364 + t571 * t367;
t631 = t576 * t678;
t449 = pkin(3) * t511 + pkin(8) * t515 + t631;
t445 = t578 * t449;
t447 = t507 * t572 + t508 * t570;
t470 = pkin(3) * t523 - pkin(8) * t525 - t555;
t601 = qJ(5) * t515 - qJD(5) * t525;
t378 = pkin(4) * t511 - t447 * t575 + t445 + t601 * t578 + (-t478 + (qJ(5) * t525 - t470) * t575) * qJD(4);
t627 = t578 * t447 + t575 * t449 + t470 * t637;
t382 = -qJ(5) * t625 + (-qJD(4) * t482 + t601) * t575 + t627;
t361 = t569 * t378 + t571 * t382;
t396 = pkin(4) * t501 + t402;
t372 = t569 * t396 + t399;
t380 = t569 * t401 + t571 * t408;
t457 = t578 * t470;
t407 = pkin(4) * t523 - qJ(5) * t667 - t482 * t575 + t457;
t641 = t575 * t470 + t478;
t415 = -qJ(5) * t668 + t641;
t386 = t569 * t407 + t571 * t415;
t649 = -pkin(5) * t645 + qJ(6) * t644 - qJD(6) * t524 + t692;
t647 = pkin(5) * t513 - t648;
t373 = qJ(6) * t513 + t380;
t435 = t571 * t489 + t569 * t591;
t646 = t435 - t373;
t643 = t501 * t670 + t455;
t567 = t576 ^ 2;
t640 = -t579 ^ 2 + t567;
t376 = t402 * t571 - t676;
t635 = qJD(6) - t376;
t629 = t559 * t655;
t628 = t468 * qJ(6) + t356;
t618 = t650 * t575;
t391 = t425 * t569 + t571 * t426;
t446 = t507 * t570 - t572 * t508;
t481 = t572 * t532 + t533 * t570;
t616 = t580 * t555 + t577 * t679;
t613 = -qJD(4) * t438 - t420;
t612 = g(1) * t665 - g(2) * t664;
t609 = pkin(4) * t668 + t481;
t573 = -qJ(5) - pkin(8);
t608 = -t557 * t573 + t529 + t563;
t607 = -t459 * t637 + t410;
t606 = -pkin(2) * t576 - t559 * t573;
t605 = pkin(5) * t558 + qJ(6) * t556;
t360 = t378 * t571 - t382 * t569;
t371 = t396 * t571 - t676;
t385 = t407 * t571 - t415 * t569;
t392 = t425 * t571 - t426 * t569;
t600 = -t554 - t686;
t599 = pkin(4) * t695 + t446;
t598 = -0.2e1 * pkin(1) * t634 - pkin(7) * qJDD(2);
t502 = t559 * t656 + t652;
t597 = -t515 * t578 - t525 * t638;
t595 = t458 * t501 - t550 * t468;
t589 = pkin(4) * t656 + t554 * t663 - t573 * t664 + t616;
t581 = qJD(2) ^ 2;
t588 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t581 + t620;
t582 = qJD(1) ^ 2;
t587 = pkin(1) * t582 - pkin(7) * qJDD(1) + t611;
t390 = pkin(4) * t426 + qJDD(5) + t419;
t586 = t658 + t573 * t665 + pkin(4) * t655 + (-t555 - t529) * t577;
t520 = t650 * t578;
t462 = t520 * t569 + t571 * t618;
t463 = t571 * t520 - t569 * t618;
t583 = -t463 * t391 + t392 * t462 - t435 * t429 - t687;
t357 = pkin(5) * t391 - qJ(6) * t392 - qJD(6) * t602 + t390;
t552 = -pkin(3) - t686;
t551 = -pkin(4) * t571 - pkin(5);
t546 = pkin(4) * t569 + qJ(6);
t544 = pkin(4) * t653;
t505 = t559 * t652 + t656;
t504 = -t629 + t653;
t503 = -t559 * t653 + t655;
t493 = t558 * t663 + t666;
t491 = t559 * t654 - t651;
t452 = t524 * t525;
t450 = -pkin(5) * t691 - qJ(6) * t524 + t600;
t417 = t509 * t525 + t515 * t691;
t416 = -qJD(4) * t453 + t515 * t524;
t397 = pkin(5) * t452 - qJ(6) * t453 + t609;
t394 = pkin(4) * t486 + pkin(5) * t602 + qJ(6) * t429;
t384 = -pkin(5) * t523 - t385;
t383 = qJ(6) * t523 + t386;
t370 = qJ(6) * t501 + t372;
t369 = -pkin(5) * t501 + qJD(6) - t371;
t368 = -pkin(5) * t416 + qJ(6) * t417 - qJD(6) * t453 + t599;
t359 = -pkin(5) * t511 - t360;
t358 = qJ(6) * t511 + qJD(6) * t523 + t361;
t354 = -t621 - t685;
t353 = qJD(6) * t501 + t628;
t1 = [qJDD(1) * MDP(1) + (qJDD(1) * t567 + 0.2e1 * t576 * t622) * MDP(4) + 0.2e1 * (t576 * t632 - t634 * t640) * MDP(5) + (qJDD(2) * t576 + t579 * t581) * MDP(6) + (qJDD(2) * t579 - t576 * t581) * MDP(7) + (t576 * t598 + t579 * t588) * MDP(9) + (-t576 * t588 + t579 * t598) * MDP(10) + (-t421 * t525 - t422 * t523 + t446 * t513 + t447 * t510 + t471 * t515 - t472 * t511 + t473 * t482 + t474 * t481 - t611) * MDP(11) + (t422 * t482 + t472 * t447 - t421 * t481 - t471 * t446 - t592 * t555 + t531 * t631 - g(1) * (-t555 * t577 + t658) - g(2) * t616) * MDP(12) + (t425 * t667 + t486 * t597) * MDP(13) + (-(-t484 * t578 - t486 * t575) * t515 + (-t675 - t426 * t578 + (t484 * t575 - t486 * t578) * qJD(4)) * t525) * MDP(14) + (t425 * t523 + t455 * t525 + t486 * t511 + t501 * t597) * MDP(15) + (-t426 * t523 - t484 * t511 - t501 * t695 - t525 * t657) * MDP(16) + (t468 * t523 + t501 * t511) * MDP(17) + ((-t482 * t637 + t445) * t501 + t457 * t468 + t607 * t523 + t411 * t511 + t446 * t484 + t481 * t426 + t458 * t625 - g(1) * t503 - g(2) * t505 + ((-qJD(4) * t470 - t447) * t501 - t482 * t468 + t613 * t523 + t419 * t525 - t458 * t515) * t575) * MDP(18) + (-(-t482 * t638 + t627) * t501 - t641 * t468 - t596 * t523 - t412 * t511 + t446 * t486 + t481 * t425 + t419 * t667 - g(1) * t502 - g(2) * t504 + t597 * t458) * MDP(19) + (-t355 * t453 - t356 * t452 - t360 * t602 - t361 * t429 + t371 * t417 + t372 * t416 - t385 * t392 - t386 * t391 + t612) * MDP(20) + (-g(1) * t586 - g(2) * t589 + t355 * t385 + t356 * t386 + t371 * t360 + t372 * t361 + t390 * t609 + t427 * t599) * MDP(21) + (g(1) * t491 - g(2) * t493 - t354 * t523 + t357 * t452 - t359 * t501 + t368 * t429 - t369 * t511 - t384 * t468 - t388 * t416 + t391 * t397) * MDP(22) + (-t353 * t452 + t354 * t453 - t358 * t429 + t359 * t602 - t369 * t417 + t370 * t416 - t383 * t391 + t384 * t392 + t612) * MDP(23) + (g(1) * t490 - g(2) * t492 + t353 * t523 - t357 * t453 + t358 * t501 - t368 * t602 + t370 * t511 + t383 * t468 + t388 * t417 - t392 * t397) * MDP(24) + (t353 * t383 + t370 * t358 + t357 * t397 + t388 * t368 + t354 * t384 + t369 * t359 - g(1) * (-pkin(5) * t491 - qJ(6) * t490 + t586) - g(2) * (pkin(5) * t493 + qJ(6) * t492 + t589)) * MDP(25) + t620 * MDP(2) + t611 * MDP(3); MDP(6) * t633 + MDP(7) * t632 + qJDD(2) * MDP(8) + (t576 * t587 - t680) * MDP(9) + (g(3) * t576 + t579 * t587) * MDP(10) + ((t472 - t476) * t513 + (t471 - t477) * t510 + (t473 * t570 - t474 * t572) * pkin(2)) * MDP(11) + (t471 * t476 - t472 * t477 + (-t680 + t421 * t572 + t422 * t570 + (-qJD(1) * t531 + t611) * t576) * pkin(2)) * MDP(12) + (t486 * t614 + t675) * MDP(13) + ((t425 - t674) * t578 + (-t426 - t672) * t575) * MDP(14) + (t501 * t614 + t657 - t671) * MDP(15) + (-t501 * t638 + t643 + t673) * MDP(16) - t501 * t513 * MDP(17) + (-t411 * t513 + t552 * t426 - t444 * t501 - t476 * t484 + (t477 * t501 + t595) * t575 + t689 * t578) * MDP(18) + (t412 * t513 + t552 * t425 - t476 * t486 + t642 * t501 - t575 * t689 + t595 * t578) * MDP(19) + (-t355 * t524 + t356 * t691 + t371 * t644 + t372 * t645 + t380 * t429 - t602 * t648 + t583) * MDP(20) + (t356 * t463 - t355 * t462 + t390 * t600 - g(3) * t608 + t692 * t427 + (t435 - t380) * t372 + t648 * t371 + t611 * (t554 * t557 - t606)) * MDP(21) + (-t357 * t691 + t369 * t513 - t388 * t645 + t391 * t450 + t429 * t649 - t462 * t468 - t501 * t647 + t558 * t590) * MDP(22) + (t353 * t691 + t354 * t524 - t369 * t644 + t370 * t645 + t373 * t429 + t602 * t647 + t583) * MDP(23) + (-t357 * t524 - t370 * t513 + t388 * t644 - t392 * t450 + t463 * t468 + t501 * t646 + t556 * t590 - t602 * t649) * MDP(24) + (t353 * t463 + t357 * t450 + t354 * t462 - g(3) * (t559 * t605 + t608) + t649 * t388 + t646 * t370 + t647 * t369 + t611 * (-(-t554 - t605) * t557 - t606)) * MDP(25) + (-MDP(4) * t576 * t579 + MDP(5) * t640) * t582; (-t510 ^ 2 - t513 ^ 2) * MDP(11) + (t471 * t513 - t472 * t510 + t592 - t620) * MDP(12) + (t643 - t673) * MDP(18) + (-t657 - t671) * MDP(19) + (t355 * t691 + t356 * t524 + t371 * t645 - t372 * t644 - t427 * t513 - t620) * MDP(21) + (-t429 * t513 + t468 * t691) * MDP(22) + t524 * t468 * MDP(24) + (t353 * t524 - t354 * t691 - t369 * t645 - t370 * t644 - t388 * t513 - t620) * MDP(25) + (t513 * MDP(24) - t645 * t690) * t602 + (-MDP(18) * t638 - MDP(19) * t614 + MDP(22) * t645 - MDP(24) * t644) * t501 + t690 * (-t524 * t391 - t392 * t691 + t429 * t644); t486 * t484 * MDP(13) + (-t484 ^ 2 + t486 ^ 2) * MDP(14) + (t425 + t674) * MDP(15) + (-t426 + t672) * MDP(16) + t468 * MDP(17) + (-g(1) * t504 + g(2) * t502 + t412 * t501 - t458 * t486 + (t613 + t681) * t575 + t607) * MDP(18) + (g(1) * t505 - g(2) * t503 + t411 * t501 + t458 * t484 + t578 * t681 - t596) * MDP(19) + (t372 * t602 - t677 + (-t391 * t569 - t392 * t571) * pkin(4) + (-t371 + t376) * t429) * MDP(20) + (-g(1) * t544 + t371 * t375 - t372 * t376 + (g(2) * t652 + t355 * t571 + t356 * t569 - t427 * t486 + t575 * t687) * pkin(4)) * MDP(21) + (t375 * t501 - t394 * t429 + (pkin(5) - t551) * t468 + t688) * MDP(22) + (t370 * t602 - t391 * t546 + t392 * t551 - t677 + (t369 - t635) * t429) * MDP(23) + (-t558 * t681 - g(1) * t493 - g(2) * t491 - t388 * t429 + t394 * t602 + t468 * t546 + (0.2e1 * qJD(6) - t376) * t501 + t628) * MDP(24) + (t353 * t546 + t354 * t551 - t388 * t394 - t369 * t375 - g(1) * (-pkin(4) * t629 - pkin(5) * t492 + qJ(6) * t493 + t544) - g(2) * (-pkin(4) * t502 - pkin(5) * t490 + qJ(6) * t491) + t635 * t370 - (-pkin(4) * t575 - pkin(5) * t556 + qJ(6) * t558) * t681) * MDP(25); (t371 * t602 + t372 * t429 + t390 + t547) * MDP(21) + (t501 * t602 + t391) * MDP(22) + (-t392 + t696) * MDP(24) + (-t369 * t602 + t370 * t429 + t357 + t547) * MDP(25) + (-MDP(21) - MDP(25)) * t693 + t690 * (-t429 ^ 2 - t694); (t429 * t602 - qJDD(4) + t473) * MDP(22) + (t392 + t696) * MDP(23) + (-t501 ^ 2 - t694) * MDP(24) + (-t370 * t501 - t685 - t688) * MDP(25);];
tau  = t1;
