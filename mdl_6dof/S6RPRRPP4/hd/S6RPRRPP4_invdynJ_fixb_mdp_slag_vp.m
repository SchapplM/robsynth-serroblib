% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% MDP [27x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRPP4_convert_par2_MPV_fixb.m
% 
% Output:
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational Coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP4_invdynJ_fixb_mdp_slag_vp(qJ, qJD, qJDD, g, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(27,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP4_invdynJ_fixb_mdp_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_invdynJ_fixb_mdp_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_invdynJ_fixb_mdp_slag_vp: pkin has to be [10x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [27 1]), ...
  'S6RPRRPP4_invdynJ_fixb_mdp_slag_vp: MDP has to be [27x1] (double)'); 

%% Symbolic Calculation
% From invdyn_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:41:01
% EndTime: 2019-03-09 04:41:11
% DurationCPUTime: 8.20s
% Computational Cost: add. (9334->564), mult. (21992->695), div. (0->0), fcn. (16718->14), ass. (0->229)
t568 = sin(pkin(9));
t570 = cos(pkin(9));
t574 = sin(qJ(3));
t577 = cos(qJ(3));
t520 = t568 * t577 + t570 * t574;
t567 = sin(pkin(10));
t569 = cos(pkin(10));
t573 = sin(qJ(4));
t576 = cos(qJ(4));
t693 = -t567 * t573 + t569 * t576;
t448 = t693 * t520;
t634 = t576 * qJD(3);
t692 = t520 * qJD(1);
t478 = t573 * t692 - t634;
t480 = qJD(3) * t573 + t576 * t692;
t431 = t569 * t478 + t480 * t567;
t639 = qJD(1) * t577;
t542 = t570 * t639;
t640 = qJD(1) * t574;
t622 = t568 * t640;
t508 = t542 - t622;
t497 = qJD(4) - t508;
t703 = t431 * t497;
t519 = t567 * t576 + t569 * t573;
t507 = t519 * qJD(4);
t649 = t519 * t508 - t507;
t635 = qJD(4) * t576;
t636 = qJD(4) * t573;
t648 = t693 * t508 + t567 * t636 - t569 * t635;
t553 = pkin(4) * t576 + pkin(3);
t565 = pkin(9) + qJ(3);
t557 = cos(t565);
t525 = t557 * t553;
t551 = pkin(2) * t570 + pkin(1);
t702 = -t551 - t525;
t632 = qJD(1) * qJD(2);
t682 = pkin(7) + qJ(2);
t688 = qJDD(1) * t682 + t632;
t491 = t688 * t568;
t492 = t688 * t570;
t528 = t682 * t568;
t522 = qJD(1) * t528;
t529 = t682 * t570;
t523 = qJD(1) * t529;
t637 = qJD(3) * t577;
t638 = qJD(3) * t574;
t598 = -t491 * t577 - t574 * t492 + t522 * t638 - t523 * t637;
t423 = -qJDD(3) * pkin(3) - t598;
t547 = g(3) * t557;
t555 = sin(t565);
t575 = sin(qJ(1));
t578 = cos(qJ(1));
t611 = g(1) * t578 + g(2) * t575;
t595 = t611 * t555;
t586 = t595 - t547;
t701 = -qJD(4) * pkin(8) * t497 - t423 + t586;
t518 = t568 * t574 - t577 * t570;
t511 = t518 * qJD(3);
t621 = t520 * t635;
t700 = -t511 * t573 + t621;
t695 = g(1) * t575 - g(2) * t578;
t699 = qJDD(2) - t695;
t604 = -t478 * t567 + t569 * t480;
t698 = t604 ^ 2;
t615 = t497 * t576;
t463 = pkin(3) * t692 - pkin(8) * t508;
t452 = t576 * t463;
t696 = -t522 * t577 - t574 * t523;
t407 = -qJ(5) * t508 * t576 + pkin(4) * t692 - t573 * t696 + t452;
t643 = t573 * t463 + t576 * t696;
t672 = t508 * t573;
t412 = -qJ(5) * t672 + t643;
t681 = qJ(5) + pkin(8);
t617 = qJD(4) * t681;
t503 = qJD(5) * t576 - t573 * t617;
t589 = -qJD(5) * t573 - t576 * t617;
t647 = (-t407 + t589) * t569 + (t412 - t503) * t567;
t474 = t528 * t577 + t574 * t529;
t470 = -t574 * t522 + t577 * t523;
t694 = -t470 + (t636 - t672) * pkin(4);
t691 = qJ(2) * qJDD(1);
t690 = MDP(22) + MDP(25);
t461 = -qJD(3) * pkin(3) - t696;
t429 = pkin(4) * t478 + qJD(5) + t461;
t391 = pkin(5) * t431 - qJ(6) * t604 + t429;
t566 = qJ(4) + pkin(10);
t558 = cos(t566);
t663 = t558 * t578;
t556 = sin(t566);
t664 = t556 * t575;
t484 = t557 * t664 + t663;
t651 = t578 * t556;
t654 = t575 * t558;
t486 = t557 * t651 - t654;
t526 = -qJD(1) * t551 + qJD(2);
t439 = -pkin(3) * t508 - pkin(8) * t692 + t526;
t462 = qJD(3) * pkin(8) + t470;
t415 = t439 * t573 + t462 * t576;
t630 = qJDD(1) * t577;
t631 = qJDD(1) * t574;
t623 = qJD(3) * t542 + t568 * t630 + t570 * t631;
t467 = -qJD(3) * t622 + t623;
t512 = t520 * qJD(3);
t538 = t570 * t630;
t606 = t568 * t631 - t538;
t468 = qJD(1) * t512 + t606;
t524 = -qJDD(1) * t551 + qJDD(2);
t424 = pkin(3) * t468 - pkin(8) * t467 + t524;
t421 = t576 * t424;
t603 = -t491 * t574 + t492 * t577;
t422 = qJDD(3) * pkin(8) + qJD(3) * t696 + t603;
t427 = qJD(4) * t634 + t573 * qJDD(3) + t576 * t467 - t636 * t692;
t460 = qJDD(4) + t468;
t368 = pkin(4) * t460 - qJ(5) * t427 - qJD(4) * t415 - qJD(5) * t480 - t422 * t573 + t421;
t428 = t480 * qJD(4) - t576 * qJDD(3) + t467 * t573;
t592 = t576 * t422 + t573 * t424 + t439 * t635 - t462 * t636;
t371 = -qJ(5) * t428 - qJD(5) * t478 + t592;
t359 = t569 * t368 - t567 * t371;
t619 = -qJDD(6) + t359;
t683 = g(3) * t555;
t689 = g(1) * t486 + g(2) * t484 - t391 * t604 + t556 * t683 + t619;
t587 = t557 * t611 + t683;
t687 = pkin(5) * t460;
t680 = qJDD(1) * pkin(1);
t406 = -qJ(5) * t478 + t415;
t402 = t569 * t406;
t414 = t576 * t439 - t462 * t573;
t405 = -qJ(5) * t480 + t414;
t377 = t405 * t567 + t402;
t679 = t377 * t604;
t678 = t406 * t567;
t677 = t427 * t573;
t676 = t478 * t497;
t675 = t478 * t692;
t674 = t480 * t497;
t673 = t480 * t692;
t670 = t520 * t573;
t669 = t520 * t576;
t666 = t555 * t575;
t665 = t555 * t578;
t659 = t570 * MDP(4);
t657 = t573 * t460;
t656 = t573 * t575;
t655 = t573 * t578;
t653 = t575 * t576;
t449 = t576 * t460;
t475 = -t528 * t574 + t529 * t577;
t471 = t576 * t475;
t652 = t576 * t578;
t360 = t567 * t368 + t569 * t371;
t440 = -t518 * qJD(2) - qJD(3) * t474;
t464 = pkin(3) * t512 + pkin(8) * t511;
t453 = t576 * t464;
t466 = pkin(3) * t518 - pkin(8) * t520 - t551;
t600 = qJ(5) * t511 - qJD(5) * t520;
t382 = pkin(4) * t512 - t440 * t573 + t453 + t600 * t576 + (-t471 + (qJ(5) * t520 - t466) * t573) * qJD(4);
t625 = t576 * t440 + t573 * t464 + t466 * t635;
t386 = -qJ(5) * t621 + (-qJD(4) * t475 + t600) * t573 + t625;
t365 = t567 * t382 + t569 * t386;
t400 = pkin(4) * t497 + t405;
t376 = t567 * t400 + t402;
t385 = t567 * t407 + t569 * t412;
t455 = t576 * t466;
t409 = pkin(4) * t518 - qJ(5) * t669 - t475 * t573 + t455;
t642 = t573 * t466 + t471;
t416 = -qJ(5) * t670 + t642;
t390 = t567 * t409 + t569 * t416;
t650 = -pkin(5) * t649 + qJ(6) * t648 - qJD(6) * t519 + t694;
t646 = pkin(5) * t692 - t647;
t380 = qJ(6) * t692 + t385;
t445 = t569 * t503 + t567 * t589;
t645 = t445 - t380;
t644 = t497 * t672 + t449;
t641 = t568 ^ 2 + t570 ^ 2;
t378 = t405 * t569 - t678;
t633 = qJD(6) - t378;
t627 = t557 * t655;
t626 = t460 * qJ(6) + t360;
t620 = t681 * t573;
t395 = t427 * t567 + t569 * t428;
t614 = -qJD(4) * t439 - t422;
t441 = qJD(2) * t520 - t528 * t638 + t529 * t637;
t613 = 0.2e1 * t641;
t612 = -g(1) * t666 + g(2) * t665;
t609 = pkin(4) * t670 + t474;
t608 = -t462 * t635 + t421;
t607 = pkin(5) * t558 + qJ(6) * t556;
t364 = t382 * t569 - t386 * t567;
t375 = t400 * t569 - t678;
t389 = t409 * t569 - t416 * t567;
t396 = t427 * t569 - t428 * t567;
t599 = t680 - t699;
t596 = pkin(4) * t700 + t441;
t498 = t557 * t656 + t652;
t594 = -t511 * t576 - t520 * t636;
t593 = -pkin(8) * t460 + t461 * t497;
t585 = pkin(4) * t656 + t682 * t575 - t578 * t702 + t681 * t665;
t584 = pkin(4) * t655 + t575 * t702 + t682 * t578 - t681 * t666;
t582 = t613 * t632 - t611;
t394 = pkin(4) * t428 + qJDD(5) + t423;
t530 = t681 * t576;
t476 = t530 * t567 + t569 * t620;
t477 = t569 * t530 - t567 * t620;
t581 = -t477 * t395 + t396 * t476 - t445 * t431 - t587;
t361 = pkin(5) * t395 - qJ(6) * t396 - qJD(6) * t604 + t394;
t550 = -pkin(4) * t569 - pkin(5);
t546 = pkin(4) * t567 + qJ(6);
t544 = pkin(4) * t653;
t501 = t557 * t652 + t656;
t500 = -t627 + t653;
t499 = -t557 * t653 + t655;
t487 = t557 * t663 + t664;
t485 = t557 * t654 - t651;
t465 = -pkin(5) * t693 - qJ(6) * t519 - t553;
t447 = t519 * t520;
t419 = t507 * t520 + t511 * t693;
t418 = -qJD(4) * t448 + t511 * t519;
t401 = pkin(5) * t447 - qJ(6) * t448 + t609;
t398 = pkin(4) * t480 + pkin(5) * t604 + qJ(6) * t431;
t388 = -pkin(5) * t518 - t389;
t387 = qJ(6) * t518 + t390;
t374 = qJ(6) * t497 + t376;
t373 = -pkin(5) * t497 + qJD(6) - t375;
t372 = -pkin(5) * t418 + qJ(6) * t419 - qJD(6) * t448 + t596;
t363 = -pkin(5) * t512 - t364;
t362 = qJ(6) * t512 + qJD(6) * t518 + t365;
t358 = -t619 - t687;
t357 = qJD(6) * t497 + t626;
t1 = [(pkin(1) * t599 + (t641 * t691 + t582) * qJ(2)) * MDP(7) + (t613 * t691 + t582) * MDP(6) + t611 * MDP(3) + (-g(1) * t584 - g(2) * t585 + t359 * t389 + t360 * t390 + t375 * t364 + t376 * t365 + t394 * t609 + t429 * t596) * MDP(23) + ((-t475 * t635 + t453) * t497 + t455 * t460 + t608 * t518 + t414 * t512 + t441 * t478 + t474 * t428 + t461 * t621 - g(1) * t499 - g(2) * t501 + ((-qJD(4) * t466 - t440) * t497 - t475 * t460 + t614 * t518 + t423 * t520 - t461 * t511) * t573) * MDP(20) + (-qJD(3) * t440 - qJDD(3) * t475 - t467 * t551 - t511 * t526 + t520 * t524 + t612) * MDP(14) + (-(-t478 * t576 - t480 * t573) * t511 + (-t677 - t428 * t576 + (t478 * t573 - t480 * t576) * qJD(4)) * t520) * MDP(16) + (-qJD(3) * t511 + qJDD(3) * t520) * MDP(10) + (-t467 * t518 - t468 * t520 - t508 * t511 - t512 * t692) * MDP(9) + (t467 * t520 - t511 * t692) * MDP(8) + (-qJD(3) * t441 - qJDD(3) * t474 - t468 * t551 + t512 * t526 + t518 * t524 + t557 * t695) * MDP(13) + t695 * MDP(2) + (-MDP(5) * t568 + t659) * (t599 + t680) + (-t359 * t448 - t360 * t447 - t364 * t604 - t365 * t431 + t375 * t419 + t376 * t418 - t389 * t396 - t390 * t395 - t612) * MDP(22) + (-t357 * t447 + t358 * t448 - t362 * t431 + t363 * t604 - t373 * t419 + t374 * t418 - t387 * t395 + t388 * t396 - t612) * MDP(25) + (g(1) * t484 - g(2) * t486 + t357 * t518 - t361 * t448 + t362 * t497 - t372 * t604 + t374 * t512 + t387 * t460 + t391 * t419 - t396 * t401) * MDP(26) + (t357 * t387 + t374 * t362 + t361 * t401 + t391 * t372 + t358 * t388 + t373 * t363 - g(1) * (-pkin(5) * t485 - qJ(6) * t484 + t584) - g(2) * (pkin(5) * t487 + qJ(6) * t486 + t585)) * MDP(27) + (-t428 * t518 - t478 * t512 - t700 * t497 - t520 * t657) * MDP(18) + (t427 * t518 + t449 * t520 + t480 * t512 + t497 * t594) * MDP(17) + (t427 * t669 + t480 * t594) * MDP(15) + (-(-t475 * t636 + t625) * t497 - t642 * t460 - t592 * t518 - t415 * t512 + t441 * t480 + t474 * t427 + t423 * t669 - g(1) * t498 - g(2) * t500 + t594 * t461) * MDP(21) + (t460 * t518 + t497 * t512) * MDP(19) + (g(1) * t485 - g(2) * t487 - t358 * t518 + t361 * t447 - t363 * t497 + t372 * t431 - t373 * t512 - t388 * t460 - t391 * t418 + t395 * t401) * MDP(24) + (-qJD(3) * t512 - qJDD(3) * t518) * MDP(11) + qJDD(1) * MDP(1); t699 * MDP(7) - t538 * MDP(13) + t623 * MDP(14) + (t644 - t675) * MDP(20) + (-t657 - t673) * MDP(21) + (t359 * t693 + t360 * t519 + t375 * t649 - t376 * t648 - t429 * t692 - t695) * MDP(23) + (-t431 * t692 + t460 * t693) * MDP(24) + t519 * t460 * MDP(26) + (t357 * t519 - t358 * t693 - t373 * t649 - t374 * t648 - t391 * t692 - t695) * MDP(27) + (MDP(26) * t692 - t649 * t690) * t604 + (-t659 - pkin(1) * MDP(7) + (MDP(13) * t574 + MDP(5)) * t568) * qJDD(1) + (-MDP(20) * t636 - MDP(21) * t615 + MDP(24) * t649 - MDP(26) * t648) * t497 + ((t568 * t639 + t570 * t640 + t692) * MDP(13) + (t508 - t622) * MDP(14)) * qJD(3) + (-qJ(2) * MDP(7) - MDP(6)) * qJD(1) ^ 2 * t641 + t690 * (-t519 * t395 - t396 * t693 + t431 * t648); -t508 ^ 2 * MDP(9) + ((-t508 - t622) * qJD(3) + t623) * MDP(10) - t606 * MDP(11) + qJDD(3) * MDP(12) + (qJD(3) * t470 + t586 + t598) * MDP(13) + (-t508 * t526 + t587 - t603) * MDP(14) + (t480 * t615 + t677) * MDP(15) + ((t427 - t676) * t576 + (-t428 - t674) * t573) * MDP(16) + (t497 * t615 + t657 - t673) * MDP(17) + (-t636 * t497 + t644 + t675) * MDP(18) + (-pkin(3) * t428 - t452 * t497 - t470 * t478 + (t497 * t696 + t593) * t573 + t701 * t576) * MDP(20) + (-pkin(3) * t427 - t470 * t480 + t643 * t497 - t573 * t701 + t593 * t576) * MDP(21) + (-t359 * t519 + t360 * t693 + t375 * t648 + t376 * t649 + t385 * t431 - t604 * t647 + t581) * MDP(22) + (t360 * t477 - t359 * t476 - t394 * t553 - g(3) * (t555 * t681 + t525) + t694 * t429 + (t445 - t385) * t376 + t647 * t375 + t611 * (t553 * t555 - t557 * t681)) * MDP(23) + (-t361 * t693 - t391 * t649 + t395 * t465 + t431 * t650 - t460 * t476 - t497 * t646 + t558 * t586) * MDP(24) + (t357 * t693 + t358 * t519 - t373 * t648 + t374 * t649 + t380 * t431 + t604 * t646 + t581) * MDP(25) + (-t361 * t519 + t391 * t648 - t396 * t465 + t460 * t477 + t497 * t645 + t556 * t586 - t604 * t650) * MDP(26) + (-g(3) * t525 + t357 * t477 + t358 * t476 + t361 * t465 + t650 * t391 + t645 * t374 + t646 * t373 + (-g(3) * t607 - t611 * t681) * t557 + (-g(3) * t681 + t611 * (t553 + t607)) * t555) * MDP(27) + (-t526 * MDP(13) - t497 * MDP(19) - t414 * MDP(20) + t415 * MDP(21) + t373 * MDP(24) - t374 * MDP(26) - MDP(8) * t508 + t692 * MDP(9)) * t692; t480 * t478 * MDP(15) + (-t478 ^ 2 + t480 ^ 2) * MDP(16) + (t427 + t676) * MDP(17) + (-t428 + t674) * MDP(18) + t460 * MDP(19) + (-g(1) * t500 + g(2) * t498 + t415 * t497 - t461 * t480 + (t614 + t683) * t573 + t608) * MDP(20) + (g(1) * t501 - g(2) * t499 + t414 * t497 + t461 * t478 + t576 * t683 - t592) * MDP(21) + (t376 * t604 - t679 + (-t395 * t567 - t396 * t569) * pkin(4) + (-t375 + t378) * t431) * MDP(22) + (-g(1) * t544 + t375 * t377 - t376 * t378 + (g(2) * t652 + t359 * t569 + t360 * t567 - t429 * t480 + t573 * t587) * pkin(4)) * MDP(23) + (t377 * t497 - t398 * t431 + (pkin(5) - t550) * t460 + t689) * MDP(24) + (t374 * t604 - t395 * t546 + t396 * t550 - t679 + (t373 - t633) * t431) * MDP(25) + (-t558 * t683 - g(1) * t487 - g(2) * t485 - t391 * t431 + t398 * t604 + t460 * t546 + (0.2e1 * qJD(6) - t378) * t497 + t626) * MDP(26) + (t357 * t546 + t358 * t550 - t391 * t398 - t373 * t377 - g(1) * (-pkin(4) * t627 - pkin(5) * t486 + qJ(6) * t487 + t544) - g(2) * (-pkin(4) * t498 - pkin(5) * t484 + qJ(6) * t485) + t633 * t374 - (-pkin(4) * t573 - pkin(5) * t556 + qJ(6) * t558) * t683) * MDP(27); (t375 * t604 + t376 * t431 + t394 + t547) * MDP(23) + (t497 * t604 + t395) * MDP(24) + (-t396 + t703) * MDP(26) + (-t373 * t604 + t374 * t431 + t361 + t547) * MDP(27) + (-MDP(23) - MDP(27)) * t595 + t690 * (-t431 ^ 2 - t698); (-qJD(3) * t692 + t431 * t604 - qJDD(4) - t606) * MDP(24) + (t396 + t703) * MDP(25) + (-t497 ^ 2 - t698) * MDP(26) + (-t374 * t497 - t687 - t689) * MDP(27);];
tau  = t1;
