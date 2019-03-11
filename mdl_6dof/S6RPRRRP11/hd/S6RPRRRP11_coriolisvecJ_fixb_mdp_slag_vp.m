% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRP11_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRP11_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [12x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RPRRRP11_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:41:31
% EndTime: 2019-03-09 06:41:50
% DurationCPUTime: 11.46s
% Computational Cost: add. (14557->555), mult. (47901->766), div. (0->0), fcn. (40859->12), ass. (0->231)
t551 = sin(pkin(12));
t557 = sin(qJ(3));
t553 = sin(pkin(6));
t643 = qJD(1) * t553;
t627 = t557 * t643;
t554 = cos(pkin(12));
t684 = cos(pkin(7));
t691 = cos(qJ(3));
t606 = t684 * t691;
t592 = t606 * t643;
t552 = sin(pkin(7));
t685 = cos(pkin(6));
t621 = t685 * t552;
t599 = t691 * t621;
t647 = -qJD(1) * t599 - t554 * t592;
t491 = t551 * t627 + t647;
t587 = qJD(4) + t491;
t619 = qJD(1) * t685;
t610 = pkin(1) * t619;
t541 = t554 * t610;
t670 = t551 * t553;
t570 = t685 * pkin(2) + (-pkin(9) * t684 - qJ(2)) * t670;
t490 = qJD(1) * t570 + t541;
t671 = t551 * t552;
t514 = (-pkin(2) * t554 - pkin(9) * t671 - pkin(1)) * t553;
t506 = qJD(1) * t514 + qJD(2);
t467 = -t490 * t552 + t684 * t506;
t622 = t557 * t684;
t502 = t553 * (t551 * t691 + t554 * t622) + t557 * t621;
t495 = qJD(1) * t502;
t417 = pkin(3) * t491 - pkin(10) * t495 + t467;
t620 = t684 * t490;
t478 = t557 * t620;
t628 = t554 * t643;
t517 = qJ(2) * t628 + t551 * t610;
t667 = t553 * t554;
t577 = (t667 * t684 + t621) * pkin(9);
t486 = qJD(1) * t577 + t517;
t480 = t691 * t486;
t668 = t552 * t557;
t441 = t506 * t668 + t478 + t480;
t602 = t685 * t684;
t533 = t552 * t628;
t634 = qJD(3) - t533;
t693 = qJD(1) * t602 + t634;
t421 = pkin(10) * t693 + t441;
t556 = sin(qJ(4));
t559 = cos(qJ(4));
t387 = t559 * t417 - t556 * t421;
t388 = t556 * t417 + t559 * t421;
t470 = t556 * t495 - t559 * t693;
t472 = t559 * t495 + t556 * t693;
t712 = t467 * MDP(13) + t472 * MDP(17) - t470 * MDP(18) + t387 * MDP(20) - t388 * MDP(21) - MDP(9) * t495;
t586 = t551 * t622 - t554 * t691;
t579 = t553 * t586;
t512 = qJD(1) * t579;
t624 = qJD(3) * t691;
t711 = -t552 * t624 - t512;
t710 = t441 - t587 * (pkin(4) * t556 - pkin(11) * t559);
t494 = t502 * qJD(3);
t709 = qJD(1) * t494;
t521 = t556 * t668 - t559 * t684;
t629 = t551 * t643;
t613 = t552 * t629;
t708 = qJD(4) * t521 + t556 * t613 + t711 * t559;
t584 = t551 * t606 + t554 * t557;
t578 = t553 * t584;
t626 = qJD(3) * t668;
t707 = qJD(1) * t578 - t626;
t581 = qJD(4) * t587;
t706 = -t556 * t709 - t559 * t581;
t485 = t491 * qJD(3);
t663 = t559 * t485;
t565 = -qJD(4) * t470 - t663;
t705 = qJD(5) * t587 + t565;
t687 = -qJ(6) - pkin(11);
t704 = -qJ(6) * t470 + qJD(5) * t687;
t703 = t551 * MDP(4) + t554 * MDP(5);
t700 = (t551 ^ 2 + t554 ^ 2) * MDP(6) * t553 ^ 2;
t572 = qJD(2) * t579;
t630 = t552 * t691;
t698 = -t557 * t486 + t490 * t606 + t506 * t630;
t411 = -qJD(1) * t572 + qJD(3) * t698;
t642 = qJD(2) * t553;
t611 = t642 * t671;
t600 = qJD(1) * t611;
t453 = pkin(3) * t709 + t485 * pkin(10) + t600;
t639 = qJD(4) * t559;
t641 = qJD(4) * t556;
t614 = t556 * t411 + t417 * t641 + t421 * t639 - t559 * t453;
t369 = -pkin(4) * t709 + t614;
t555 = sin(qJ(5));
t558 = cos(qJ(5));
t637 = qJD(5) * t558;
t396 = t472 * t637 + t555 * t705 - t558 * t709;
t359 = t396 * pkin(5) + t369;
t598 = t554 * t606;
t669 = t551 * t557;
t632 = t553 * t669;
t501 = -t553 * t598 - t599 + t632;
t631 = pkin(1) * t685;
t545 = t554 * t631;
t503 = t545 + t570;
t473 = -t503 * t552 + t684 * t514;
t432 = pkin(3) * t501 - pkin(10) * t502 + t473;
t520 = t552 * t667 - t602;
t646 = qJ(2) * t667 + t551 * t631;
t498 = t577 + t646;
t569 = t691 * t498 + (t503 * t684 + t514 * t552) * t557;
t439 = -t520 * pkin(10) + t569;
t658 = t556 * t432 + t559 * t439;
t391 = pkin(11) * t501 + t658;
t696 = -t557 * t498 + t503 * t606 + t514 * t630;
t438 = t520 * pkin(3) - t696;
t474 = t502 * t556 + t520 * t559;
t475 = t502 * t559 - t520 * t556;
t403 = t474 * pkin(4) - t475 * pkin(11) + t438;
t660 = t558 * t391 + t555 * t403;
t633 = pkin(10) * t641;
t699 = t555 * t633 - t710 * t558;
t522 = t556 * t684 + t559 * t668;
t697 = -qJD(4) * t522 + t711 * t556 - t559 * t613;
t462 = t495 * pkin(3) + pkin(10) * t491;
t657 = t556 * t462 + t559 * t698;
t400 = pkin(11) * t495 + t657;
t535 = -pkin(4) * t559 - pkin(11) * t556 - pkin(3);
t695 = t558 * t400 - t535 * t637 + t710 * t555;
t694 = t557 * (t506 * t552 + t620) + t480;
t447 = t558 * t472 + t555 * t587;
t692 = t447 ^ 2;
t690 = pkin(5) * t555;
t688 = t495 * pkin(4);
t686 = pkin(10) * qJD(4);
t681 = qJ(6) * t556;
t680 = t369 * t555;
t679 = t369 * t558;
t638 = qJD(5) * t555;
t395 = t472 * t638 - t555 * t709 - t558 * t705;
t678 = t395 * t555;
t666 = t556 * t485;
t443 = qJD(4) * t472 - t666;
t677 = t443 * t555;
t676 = t443 * t558;
t445 = t472 * t555 - t558 * t587;
t469 = qJD(5) + t470;
t675 = t445 * t469;
t674 = t447 * t469;
t673 = t491 * t556;
t672 = t491 * t559;
t665 = t556 * t558;
t664 = t558 * t559;
t382 = pkin(11) * t587 + t388;
t420 = -pkin(3) * t693 - t698;
t394 = t470 * pkin(4) - t472 * pkin(11) + t420;
t365 = -t382 * t555 + t558 * t394;
t361 = -qJ(6) * t447 + t365;
t360 = pkin(5) * t469 + t361;
t662 = t360 - t361;
t430 = pkin(4) * t472 + pkin(11) * t470;
t661 = t558 * t387 + t555 * t430;
t459 = -t491 * t664 + t495 * t555;
t546 = pkin(10) * t664;
t636 = qJD(6) * t558;
t656 = pkin(5) * t673 + qJ(6) * t459 + t400 * t555 - t556 * t636 + (pkin(5) * t556 - qJ(6) * t664) * qJD(4) + (-t546 + (-t535 + t681) * t555) * qJD(5) + t699;
t588 = -t558 * t522 + t555 * t630;
t655 = -qJD(5) * t588 - t555 * t708 + t707 * t558;
t507 = -t555 * t522 - t558 * t630;
t654 = -qJD(5) * t507 + t707 * t555 + t558 * t708;
t458 = -t558 * t495 - t555 * t672;
t653 = qJ(6) * t458 + (-qJ(6) * qJD(5) - t686) * t665 + (-qJD(6) * t556 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t559) * t555 - t695;
t651 = t555 * t704 + t636 - t661;
t425 = t558 * t430;
t650 = -pkin(5) * t472 - t425 + t704 * t558 + (-qJD(6) + t387) * t555;
t645 = t555 * t535 + t546;
t640 = qJD(4) * t558;
t635 = t420 * qJD(4);
t625 = t469 * t638;
t618 = -t391 * t555 + t558 * t403;
t617 = t432 * t559 - t556 * t439;
t436 = t556 * t698;
t616 = t462 * t559 - t436;
t615 = t469 * t558;
t604 = -t555 * t639 + t458;
t603 = t558 * t639 - t459;
t366 = t382 * t558 + t394 * t555;
t455 = t475 * t558 + t501 * t555;
t454 = t475 * t555 - t501 * t558;
t601 = (-qJ(2) * t629 + t541) * t551 - t517 * t554;
t426 = qJD(3) * t696 - t572;
t493 = (t599 + (t598 - t669) * t553) * qJD(3);
t457 = pkin(3) * t494 - pkin(10) * t493 + t611;
t597 = -t556 * t426 - t432 * t641 - t439 * t639 + t457 * t559;
t390 = -pkin(4) * t501 - t617;
t595 = -t469 * t637 - t677;
t381 = -pkin(4) * t587 - t387;
t594 = -pkin(11) * t443 + t381 * t469;
t591 = -t559 * t411 - t417 * t639 + t421 * t641 - t556 * t453;
t590 = t559 * t426 + t432 * t639 - t439 * t641 + t556 * t457;
t368 = pkin(11) * t709 - t591;
t380 = t443 * pkin(4) - pkin(11) * t565 + qJD(3) * t478 + t486 * t624 + t506 * t626 + (t551 * t592 + t554 * t627) * qJD(2);
t357 = t558 * t368 + t555 * t380 - t382 * t638 + t394 * t637;
t374 = pkin(11) * t494 + t590;
t571 = qJD(2) * t578;
t427 = qJD(3) * t569 + t571;
t451 = qJD(4) * t475 + t493 * t556;
t452 = -qJD(4) * t474 + t493 * t559;
t385 = t451 * pkin(4) - t452 * pkin(11) + t427;
t589 = t558 * t374 + t555 * t385 - t391 * t638 + t403 * t637;
t375 = -pkin(4) * t494 - t597;
t582 = t491 * t587;
t358 = -qJD(5) * t366 - t368 * t555 + t558 * t380;
t573 = -qJD(5) * t660 - t374 * t555 + t558 * t385;
t537 = t687 * t558;
t536 = t687 * t555;
t529 = t558 * t535;
t509 = -t555 * t681 + t645;
t500 = -qJ(6) * t665 + t529 + (-pkin(10) * t555 - pkin(5)) * t559;
t444 = t445 ^ 2;
t412 = qJD(1) * t571 + t694 * qJD(3);
t405 = -qJD(5) * t454 + t452 * t558 + t494 * t555;
t404 = qJD(5) * t455 + t452 * t555 - t494 * t558;
t399 = -t616 - t688;
t377 = t445 * pkin(5) + qJD(6) + t381;
t364 = -qJ(6) * t454 + t660;
t363 = pkin(5) * t474 - qJ(6) * t455 + t618;
t362 = -qJ(6) * t445 + t366;
t356 = -qJ(6) * t404 - qJD(6) * t454 + t589;
t355 = pkin(5) * t451 - qJ(6) * t405 - qJD(6) * t455 + t573;
t354 = -qJ(6) * t396 - qJD(6) * t445 + t357;
t353 = pkin(5) * t443 + qJ(6) * t395 - qJD(6) * t447 + t358;
t1 = [(t411 * t520 - t426 * t693 + t467 * t493 - t473 * t485 + 0.2e1 * t495 * t611) * MDP(14) + (t485 * t520 + t493 * t693) * MDP(10) + (t412 * t520 - t427 * t693 + t473 * t709 + t491 * t611 + t501 * t600) * MDP(13) + 0.2e1 * qJD(2) * qJD(1) * t700 + (t485 * t501 - t493 * t491 - t502 * t709) * MDP(9) + (t412 * t474 + t420 * t451 + t427 * t470 + t438 * t443 - t501 * t614 + t587 * t597 + t617 * t709) * MDP(20) + (t358 * t474 + t365 * t451 + t369 * t454 + t375 * t445 + t381 * t404 + t390 * t396 + t443 * t618 + t469 * t573) * MDP(27) + (-t475 * t443 - t472 * t451 - t452 * t470 - t474 * t565) * MDP(16) + (t472 * t452 + t475 * t565) * MDP(15) + (t354 * t364 + t362 * t356 + t353 * t363 + t360 * t355 + t359 * (pkin(5) * t454 + t390) + t377 * (pkin(5) * t404 + t375)) * MDP(30) + (-t443 * t501 - t451 * t587 - t474 * t709) * MDP(18) + (t452 * t587 + t475 * t709 + t501 * t565) * MDP(17) + (-t395 * t455 + t405 * t447) * MDP(22) + (t395 * t454 - t396 * t455 - t404 * t447 - t405 * t445) * MDP(23) + (-t353 * t455 - t354 * t454 - t355 * t447 - t356 * t445 - t360 * t405 - t362 * t404 + t363 * t395 - t364 * t396) * MDP(29) + (t412 * t475 + t420 * t452 + t427 * t472 + t438 * t565 + t501 * t591 - t587 * t590 - t658 * t709) * MDP(21) + (-t357 * t474 - t366 * t451 + t369 * t455 + t375 * t447 + t381 * t405 - t390 * t395 - t443 * t660 - t469 * t589) * MDP(28) + (-t485 * t502 + t493 * t495) * MDP(8) + (t443 * t474 + t451 * t469) * MDP(26) + (-t396 * t474 - t404 * t469 - t443 * t454 - t445 * t451) * MDP(25) + (-t395 * t474 + t405 * t469 + t443 * t455 + t447 * t451) * MDP(24) + (-0.2e1 * t703 * t619 + ((t554 * t646 + (qJ(2) * t670 - t545) * t551) * qJD(1) - t601) * MDP(7)) * t642 + ((-t634 + (t520 - t602) * qJD(1)) * MDP(11) + (qJD(4) + t647 + (t501 + t632) * qJD(1)) * MDP(19) + t712) * t494; t601 * MDP(7) * t643 + (-t491 * t613 + t684 * t709 + t693 * t707) * MDP(13) + (-t684 * t485 - t512 * t693 + (-t495 * t629 - t624 * t693) * t552) * MDP(14) + (-t443 * t630 - t470 * t707 - t521 * t709) * MDP(20) + (-t472 * t707 - t522 * t709 - t565 * t630) * MDP(21) + (t396 * t521 + t443 * t507 - t445 * t697 - t469 * t655) * MDP(27) + (-t395 * t521 + t443 * t588 - t447 * t697 + t469 * t654) * MDP(28) + (t395 * t507 + t396 * t588 + t445 * t654 + t447 * t655) * MDP(29) + (t353 * t507 - t354 * t588 + t359 * t521 - t360 * t655 - t362 * t654 - t377 * t697) * MDP(30) + (t697 * MDP(20) + MDP(21) * t708) * t587 + (t553 * t685 * t703 - t700) * qJD(1) ^ 2; -t491 ^ 2 * MDP(9) + (t491 * t693 - t485) * MDP(10) - MDP(11) * t709 + (-t441 * t533 + (t441 * t602 - t584 * t642) * qJD(1) + (t441 - t694) * qJD(3)) * MDP(13) + (-t698 * t533 + t467 * t491 + (t586 * t642 + t602 * t698) * qJD(1)) * MDP(14) + ((qJD(4) * t693 - t485) * t556 + t587 * t472) * t559 * MDP(15) + (-t556 * t443 + t559 * t565 + (-t641 - t673) * t472 + (-t639 - t672) * t470) * MDP(16) + (t559 * t582 - t706) * MDP(17) + (t559 * t709 + (-t581 - t582) * t556) * MDP(18) + (-pkin(3) * t443 + pkin(10) * t706 - t412 * t559 + t420 * t673 - t441 * t470 + t556 * t635 - t587 * t616) * MDP(20) + (-pkin(3) * t565 + t412 * t556 + t420 * t672 - t441 * t472 + (t633 + t657) * t587 + (-pkin(10) * t709 + t635) * t559) * MDP(21) + (-t395 * t665 + (-t556 * t638 + t603) * t447) * MDP(22) + (t445 * t459 + t447 * t458 + (-t445 * t558 - t447 * t555) * t639 + (t678 - t396 * t558 + (t445 * t555 - t447 * t558) * qJD(5)) * t556) * MDP(23) + (t395 * t559 + t603 * t469 + (t447 * t587 - t625 + t676) * t556) * MDP(24) + (t396 * t559 + t604 * t469 + (-t445 * t587 + t595) * t556) * MDP(25) + (t469 * t556 * t587 - t443 * t559) * MDP(26) + (-t381 * t458 - t399 * t445 + t529 * t443 + ((-qJD(5) * t535 + t400) * t555 + t699) * t469 + (t381 * t555 * qJD(4) - t358 + (qJD(4) * t445 + t595) * pkin(10)) * t559 + (pkin(10) * t396 + t365 * t587 + t381 * t637 + t680) * t556) * MDP(27) + (-t645 * t443 - t399 * t447 - t381 * t459 + t695 * t469 + (t381 * t640 + t357 + (qJD(4) * t447 + t625) * pkin(10)) * t559 + (-t381 * t638 + t679 - t587 * t366 + (t469 * t640 - t395) * pkin(10)) * t556) * MDP(28) + (t360 * t459 + t362 * t458 + t395 * t500 - t396 * t509 - t656 * t447 - t653 * t445 + (-t360 * t558 - t362 * t555) * t639 + (-t353 * t558 - t354 * t555 + (t360 * t555 - t362 * t558) * qJD(5)) * t556) * MDP(29) + (t353 * t500 + t354 * t509 + t359 * (pkin(10) + t690) * t556 + (t688 - t436 + (t462 + t686) * t559 + (t556 * t637 - t604) * pkin(5)) * t377 + t653 * t362 + t656 * t360) * MDP(30) + (-qJD(4) * t556 ^ 2 * MDP(15) + MDP(11) * t693 - t587 * MDP(19) + t491 * MDP(8) - t712) * t495; -t470 ^ 2 * MDP(16) + (t470 * t491 - t663) * MDP(17) + t666 * MDP(18) + MDP(19) * t709 + (t388 * t587 - t614) * MDP(20) + (t387 * t587 + t420 * t470 + t591) * MDP(21) + (t447 * t615 - t678) * MDP(22) + ((-t395 - t675) * t558 + (-t396 - t674) * t555) * MDP(23) + (t469 * t615 + t677) * MDP(24) + (-t469 ^ 2 * t555 + t676) * MDP(25) + (-pkin(4) * t396 - t679 - t388 * t445 + (-pkin(11) * t637 - t425) * t469 + (t387 * t469 + t594) * t555) * MDP(27) + (pkin(4) * t395 + t680 - t388 * t447 + (pkin(11) * t638 + t661) * t469 + t594 * t558) * MDP(28) + (t395 * t536 + t396 * t537 - t650 * t447 - t651 * t445 + (-t360 * t469 + t354) * t558 + (-t362 * t469 - t353) * t555) * MDP(29) + (-t354 * t537 + t353 * t536 + t359 * (-pkin(5) * t558 - pkin(4)) + (t469 * t690 - t388) * t377 + t651 * t362 + t650 * t360) * MDP(30) + (MDP(15) * t470 + MDP(16) * t472 + MDP(18) * t491 - t420 * MDP(20) - t447 * MDP(24) + t445 * MDP(25) - t469 * MDP(26) - t365 * MDP(27) + t366 * MDP(28)) * t472; t447 * t445 * MDP(22) + (-t444 + t692) * MDP(23) + (-t395 + t675) * MDP(24) + (-t396 + t674) * MDP(25) + t443 * MDP(26) + (t366 * t469 - t381 * t447 + t358) * MDP(27) + (t365 * t469 + t381 * t445 - t357) * MDP(28) + (pkin(5) * t395 - t445 * t662) * MDP(29) + (t662 * t362 + (-t377 * t447 + t353) * pkin(5)) * MDP(30); (-t444 - t692) * MDP(29) + (t360 * t447 + t362 * t445 + t359) * MDP(30);];
tauc  = t1;
