% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
% MDP [33x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPR1_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR1_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(33,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [33 1]), ...
  'S6RRRRPR1_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [33x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:55:10
% EndTime: 2019-03-09 21:55:22
% DurationCPUTime: 7.88s
% Computational Cost: add. (12753->463), mult. (33364->609), div. (0->0), fcn. (25733->10), ass. (0->234)
t568 = cos(qJ(6));
t642 = qJD(6) * t568;
t570 = cos(qJ(3));
t571 = cos(qJ(2));
t647 = qJD(1) * t571;
t636 = t570 * t647;
t566 = sin(qJ(3));
t567 = sin(qJ(2));
t648 = qJD(1) * t567;
t637 = t566 * t648;
t519 = -t636 + t637;
t521 = -t566 * t647 - t570 * t648;
t565 = sin(qJ(4));
t569 = cos(qJ(4));
t485 = -t569 * t519 + t521 * t565;
t563 = sin(pkin(11));
t598 = t519 * t565 + t569 * t521;
t685 = cos(pkin(11));
t701 = t485 * t685 + t563 * t598;
t715 = -t701 * t568 + t642;
t560 = qJD(2) + qJD(3);
t640 = qJD(1) * qJD(2);
t635 = t571 * t640;
t491 = qJD(3) * t636 - t560 * t637 + t570 * t635;
t690 = pkin(7) + pkin(8);
t542 = t690 * t571;
t536 = qJD(1) * t542;
t526 = t570 * t536;
t541 = t690 * t567;
t534 = qJD(1) * t541;
t686 = qJD(2) * pkin(2);
t528 = -t534 + t686;
t597 = -t528 * t566 - t526;
t638 = qJD(2) * t690;
t611 = qJD(1) * t638;
t529 = t567 * t611;
t530 = t571 * t611;
t622 = t566 * t529 - t570 * t530;
t581 = qJD(3) * t597 + t622;
t424 = -pkin(9) * t491 + t581;
t688 = pkin(9) * t519;
t463 = -t597 - t688;
t645 = qJD(4) * t565;
t455 = t463 * t645;
t556 = -pkin(2) * t571 - pkin(1);
t540 = t556 * qJD(1);
t498 = pkin(3) * t519 + t540;
t714 = -t565 * t424 - t498 * t485 + t455;
t564 = sin(qJ(6));
t643 = qJD(6) * t564;
t533 = t566 * t571 + t567 * t570;
t497 = t560 * t533;
t492 = t497 * qJD(1);
t644 = qJD(4) * t569;
t422 = t569 * t491 - t565 * t492 - t519 * t644 + t521 * t645;
t579 = -qJD(4) * t598 + t491 * t565 + t492 * t569;
t397 = t422 * t685 - t563 * t579;
t559 = qJD(4) + t560;
t659 = t568 * t397 + t559 * t642;
t700 = t485 * t563 - t598 * t685;
t377 = -t643 * t700 + t659;
t376 = t377 * t568;
t432 = t559 * t564 + t568 * t700;
t681 = t397 * t564;
t378 = t432 * qJD(6) + t681;
t674 = t700 * t564;
t430 = -t568 * t559 + t674;
t713 = -t564 * t378 - t715 * t430 + t376;
t375 = t377 * t564;
t396 = t422 * t563 + t579 * t685;
t393 = t564 * t396;
t439 = qJD(6) - t701;
t678 = t432 * t700;
t712 = t485 * MDP(18) * t598 + (-t485 ^ 2 + t598 ^ 2) * MDP(19) + (-t485 * t559 + t422) * MDP(20) + (-t559 * t598 - t579) * MDP(21) + (t715 * t432 + t375) * MDP(27) + (t393 - t678) * MDP(29) + (t715 * MDP(29) - t700 * MDP(31)) * t439;
t476 = t598 * qJ(5);
t515 = t521 * pkin(9);
t522 = t566 * t536;
t623 = t570 * t528 - t522;
t462 = t515 + t623;
t454 = pkin(3) * t560 + t462;
t456 = t565 * t463;
t627 = t569 * t454 - t456;
t413 = t627 + t476;
t409 = pkin(4) * t559 + t413;
t458 = t569 * t463;
t601 = -t454 * t565 - t458;
t684 = qJ(5) * t485;
t414 = -t601 + t684;
t669 = t563 * t414;
t383 = t409 * t685 - t669;
t381 = -t559 * pkin(5) - t383;
t683 = t381 * t701;
t679 = t430 * t700;
t620 = t534 * t566 - t526;
t465 = t620 + t688;
t652 = -t570 * t534 - t522;
t466 = t515 + t652;
t555 = pkin(2) * t570 + pkin(3);
t667 = t565 * t566;
t710 = -t555 * t644 - (-t566 * t645 + (t569 * t570 - t667) * qJD(3)) * pkin(2) + t565 * t465 + t569 * t466;
t666 = t566 * t569;
t709 = -t555 * t645 + (-t566 * t644 + (-t565 * t570 - t666) * qJD(3)) * pkin(2) - t569 * t465 + t466 * t565;
t646 = qJD(3) * t566;
t621 = -t566 * t530 - t536 * t646;
t695 = t570 * (qJD(3) * t528 - t529);
t423 = -pkin(9) * t492 + t621 + t695;
t602 = -t565 * t423 + t569 * t424;
t578 = qJD(4) * t601 + t602;
t673 = t498 * t598;
t708 = t578 + t673;
t617 = -qJD(4) * t454 - t423;
t707 = t617 * t569 + t714;
t371 = t485 * qJD(5) - t455 + (t424 + (qJD(4) * t519 - t491) * qJ(5)) * t565 + ((qJD(4) * t521 - t492) * qJ(5) - t617) * t569;
t575 = -qJ(5) * t422 + qJD(5) * t598 + t578;
t358 = t371 * t563 - t575 * t685;
t410 = t685 * t414;
t384 = t563 * t409 + t410;
t382 = pkin(10) * t559 + t384;
t450 = -pkin(4) * t485 + qJD(5) + t498;
t398 = -pkin(5) * t701 - pkin(10) * t700 + t450;
t368 = t382 * t568 + t398 * t564;
t609 = t358 * t564 + t368 * t700 + t381 * t642;
t604 = t382 * t564 - t398 * t568;
t591 = -t358 * t568 + t381 * t643 + t604 * t700;
t705 = pkin(5) * t700 - pkin(10) * t701;
t704 = pkin(4) * t598;
t703 = -t476 - t710;
t702 = t684 + t709;
t698 = -0.2e1 * t640;
t697 = MDP(4) * t567;
t696 = MDP(5) * (t567 ^ 2 - t571 ^ 2);
t694 = qJD(1) * t533;
t535 = t567 * t638;
t537 = t571 * t638;
t670 = t541 * t570;
t587 = -qJD(3) * t670 - t570 * t535 - t566 * t537 - t542 * t646;
t436 = -pkin(9) * t497 + t587;
t532 = t566 * t567 - t570 * t571;
t496 = t560 * t532;
t596 = t541 * t566 - t542 * t570;
t580 = qJD(3) * t596 + t535 * t566 - t570 * t537;
t437 = pkin(9) * t496 + t580;
t477 = -pkin(9) * t533 - t542 * t566 - t670;
t478 = -pkin(9) * t532 - t596;
t600 = -t477 * t565 - t478 * t569;
t692 = qJD(4) * t600 - t436 * t565 + t569 * t437;
t359 = t371 * t685 + t563 * t575;
t495 = -t532 * t565 + t533 * t569;
t427 = qJD(4) * t495 - t496 * t565 + t569 * t497;
t494 = t569 * t532 + t533 * t565;
t586 = t569 * t436 + t565 * t437 + t477 * t644 - t478 * t645;
t373 = -qJ(5) * t427 - qJD(5) * t494 + t586;
t426 = -qJD(4) * t494 - t496 * t569 - t497 * t565;
t577 = -qJ(5) * t426 - qJD(5) * t495 + t692;
t361 = t373 * t685 + t563 * t577;
t402 = t426 * t685 - t563 * t427;
t448 = t494 * t685 + t495 * t563;
t449 = -t563 * t494 + t495 * t685;
t503 = pkin(3) * t532 + t556;
t585 = pkin(4) * t494 + t503;
t404 = pkin(5) * t448 - pkin(10) * t449 + t585;
t421 = -qJ(5) * t494 - t600;
t584 = -qJ(5) * t495 + t477 * t569 - t478 * t565;
t392 = t421 * t685 + t563 * t584;
t608 = t358 * t449 - t392 * t396;
t691 = t381 * t402 - (qJD(6) * t404 + t361) * t439 - (qJD(6) * t398 + t359) * t448 + t608;
t689 = pkin(3) * t521;
t687 = pkin(3) * qJD(4);
t682 = t381 * t449;
t680 = t404 * t396;
t677 = t432 * t564;
t671 = t540 * t521;
t668 = t563 * t565;
t572 = qJD(2) ^ 2;
t665 = t567 * t572;
t394 = t568 * t396;
t664 = t571 * t572;
t573 = qJD(1) ^ 2;
t663 = t571 * t573;
t658 = t703 * t563 - t702 * t685;
t657 = t702 * t563 + t703 * t685;
t656 = t569 * t462 - t456;
t512 = -pkin(2) * t667 + t555 * t569 + pkin(4);
t516 = pkin(2) * t666 + t555 * t565;
t472 = t563 * t512 + t685 * t516;
t416 = t476 + t656;
t626 = -t462 * t565 - t458;
t593 = t626 - t684;
t632 = t685 * t565;
t654 = -t416 * t563 + t593 * t685 + (t563 * t569 + t632) * t687;
t653 = -t416 * t685 - t563 * t593 + (t569 * t685 - t668) * t687;
t554 = pkin(3) * t569 + pkin(4);
t514 = pkin(3) * t632 + t563 * t554;
t558 = t567 * t686;
t557 = pkin(2) * t648;
t634 = -pkin(2) * t560 - t528;
t633 = -pkin(3) * t559 - t454;
t467 = t492 * pkin(3) + qJD(2) * t557;
t487 = pkin(3) * t497 + t558;
t631 = pkin(1) * t698;
t619 = t439 * t564;
t610 = -t689 - t704;
t400 = t610 + t705;
t470 = pkin(10) + t472;
t614 = qJD(6) * t470 + t400 + t557;
t507 = pkin(10) + t514;
t613 = qJD(6) * t507 + t400;
t551 = pkin(4) * t563 + pkin(10);
t612 = qJD(6) * t551 - t704 + t705;
t607 = -t396 * t551 - t683;
t606 = -t396 * t470 - t683;
t605 = -t396 * t507 - t683;
t603 = t383 * t701 + t384 * t700;
t595 = pkin(4) * t427 + t487;
t594 = t394 + (t564 * t701 - t643) * t439;
t590 = t540 * t519 - t621;
t589 = t402 * t568 - t449 * t643;
t471 = t512 * t685 - t563 * t516;
t513 = -pkin(3) * t668 + t554 * t685;
t576 = pkin(4) * t579 + t467;
t574 = -t521 * t519 * MDP(11) + (-t432 * t619 + t713) * MDP(28) + (-t439 * t619 + t394 + t679) * MDP(30) + t491 * MDP(13) + (-t519 ^ 2 + t521 ^ 2) * MDP(12) + (t519 * MDP(13) + (-t521 - t694) * MDP(14)) * t560 + t712;
t552 = -pkin(4) * t685 - pkin(5);
t506 = -pkin(5) - t513;
t499 = t557 - t689;
t469 = -pkin(5) - t471;
t401 = t426 * t563 + t427 * t685;
t391 = t421 * t563 - t584 * t685;
t386 = t413 * t685 - t669;
t385 = t413 * t563 + t410;
t369 = pkin(5) * t401 - pkin(10) * t402 + t595;
t365 = t396 * pkin(5) - t397 * pkin(10) + t576;
t364 = t568 * t365;
t360 = t373 * t563 - t577 * t685;
t1 = [MDP(6) * t664 + (t498 * t427 + t467 * t494 - t487 * t485 + t579 * t503) * MDP(23) + (t556 * t492 + t540 * t497 + (qJD(1) * t532 + t519) * t558) * MDP(16) + (t376 * t449 + t432 * t589) * MDP(27) + ((-t430 * t568 - t677) * t402 + (-t375 - t378 * t568 + (t430 * t564 - t432 * t568) * qJD(6)) * t449) * MDP(28) + (-t449 * t393 - t378 * t448 - t401 * t430 + (-t402 * t564 - t449 * t642) * t439) * MDP(30) - MDP(7) * t665 + (pkin(7) * t665 + t571 * t631) * MDP(10) + (-pkin(7) * t664 + t567 * t631) * MDP(9) + (t377 * t448 + t394 * t449 + t401 * t432 + t439 * t589) * MDP(29) + t696 * t698 + (t358 * t391 + t359 * t392 - t383 * t360 + t384 * t361 + t450 * t595 + t576 * t585) * MDP(26) + (t426 * MDP(20) - t427 * MDP(21) + MDP(23) * t692 - t586 * MDP(24)) * t559 + 0.2e1 * t635 * t697 + (t396 * t448 + t401 * t439) * MDP(31) + (-t359 * t448 + t360 * t700 + t361 * t701 - t383 * t402 - t384 * t401 + t391 * t397 + t608) * MDP(25) + (t360 * t432 - t368 * t401 + t391 * t377 + (-(-qJD(6) * t392 + t369) * t439 - t680 - (-qJD(6) * t382 + t365) * t448 - qJD(6) * t682) * t564 + t691 * t568) * MDP(33) + (t360 * t430 + t364 * t448 - t604 * t401 + t391 * t378 + (t369 * t439 + t680 + (-t382 * t448 - t392 * t439 + t682) * qJD(6)) * t568 + t691 * t564) * MDP(32) + (t503 * t422 + t498 * t426 + t467 * t495 - t487 * t598) * MDP(24) + (t422 * t495 - t426 * t598) * MDP(18) + (-t422 * t494 + t426 * t485 + t427 * t598 - t495 * t579) * MDP(19) + (t556 * t491 - t540 * t496 + (-t521 + t694) * t558) * MDP(17) + (-t496 * MDP(13) - t497 * MDP(14) + t580 * MDP(16) - t587 * MDP(17)) * t560 + (-t491 * t532 - t492 * t533 + t496 * t519 + t497 * t521) * MDP(12) + (t491 * t533 + t496 * t521) * MDP(11); (t521 * t557 + t652 * t560 + (qJD(3) * t634 + t529) * t570 + t590) * MDP(17) + (-t519 * t557 + t671 - t620 * t560 + (t566 * t634 - t526) * qJD(3) + t622) * MDP(16) + (-t396 * t472 - t397 * t471 + t657 * t701 + t658 * t700 + t603) * MDP(25) + (t359 * t472 - t358 * t471 - t450 * (t557 + t610) + t657 * t384 - t658 * t383) * MDP(26) + t574 + (t469 * t377 + t606 * t568 + t658 * t432 + (t564 * t614 - t568 * t657) * t439 + t609) * MDP(33) + (t469 * t378 + t606 * t564 + t658 * t430 + (-t564 * t657 - t568 * t614) * t439 + t591) * MDP(32) - t663 * t697 + t573 * t696 + (t499 * t598 + t559 * t710 + t707) * MDP(24) + (t499 * t485 + t559 * t709 + t708) * MDP(23) + (MDP(9) * t567 * t573 + MDP(10) * t663) * pkin(1); (-t598 * t689 + t656 * t559 + (qJD(4) * t633 - t423) * t569 + t714) * MDP(24) + (-t485 * t689 + t673 - t626 * t559 + (t565 * t633 - t458) * qJD(4) + t602) * MDP(23) + (-t396 * t514 - t397 * t513 + t653 * t701 + t654 * t700 + t603) * MDP(25) + (t506 * t378 + t605 * t564 + t654 * t430 + (-t564 * t653 - t568 * t613) * t439 + t591) * MDP(32) + (t560 * t623 + t590 - t695) * MDP(17) + (-t358 * t513 + t359 * t514 - t383 * t654 + t384 * t653 - t450 * t610) * MDP(26) + (-t560 * t597 + t581 + t671) * MDP(16) + (t506 * t377 + t605 * t568 + t654 * t432 + (t564 * t613 - t568 * t653) * t439 + t609) * MDP(33) + t574; (-t559 * t601 + t708) * MDP(23) + (t559 * t627 + t707) * MDP(24) + ((-t396 * t563 - t397 * t685) * pkin(4) + (t383 - t386) * t701 + (t384 - t385) * t700) * MDP(25) + (t383 * t385 - t384 * t386 + (-t358 * t685 + t359 * t563 + t450 * t598) * pkin(4)) * MDP(26) + (-t439 * t677 + t713) * MDP(28) + (t594 + t679) * MDP(30) + (t552 * t378 - t385 * t430 + t607 * t564 + (t564 * t386 - t568 * t612) * t439 + t591) * MDP(32) + (t552 * t377 - t385 * t432 + t607 * t568 + (t568 * t386 + t564 * t612) * t439 + t609) * MDP(33) + t712; (-t700 ^ 2 - t701 ^ 2) * MDP(25) + (t383 * t700 - t384 * t701 + t576) * MDP(26) + (t594 - t679) * MDP(32) + (-t439 ^ 2 * t568 - t393 - t678) * MDP(33); t432 * t430 * MDP(27) + (-t430 ^ 2 + t432 ^ 2) * MDP(28) + (t430 * t439 + t659) * MDP(29) + (t432 * t439 - t681) * MDP(30) + t396 * MDP(31) + (-t359 * t564 + t368 * t439 - t381 * t432 + t364) * MDP(32) + (-t359 * t568 - t365 * t564 + t381 * t430 - t439 * t604) * MDP(33) + (-MDP(29) * t674 - MDP(30) * t432 - MDP(32) * t368 + MDP(33) * t604) * qJD(6);];
tauc  = t1;
