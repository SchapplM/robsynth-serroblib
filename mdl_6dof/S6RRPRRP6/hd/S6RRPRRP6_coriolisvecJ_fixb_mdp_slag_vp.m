% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% MDP [30x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRPRRP6_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRP6_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(30,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [30 1]), ...
  'S6RRPRRP6_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [30x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:14:29
% EndTime: 2019-03-09 12:14:47
% DurationCPUTime: 10.95s
% Computational Cost: add. (12328->600), mult. (36469->774), div. (0->0), fcn. (29083->10), ass. (0->249)
t547 = sin(pkin(6));
t553 = cos(qJ(2));
t689 = cos(pkin(11));
t624 = t689 * t553;
t611 = t547 * t624;
t532 = qJD(1) * t611;
t546 = sin(pkin(11));
t550 = sin(qJ(2));
t648 = qJD(1) * t547;
t633 = t550 * t648;
t591 = t546 * t633 - t532;
t711 = qJD(4) + t591;
t690 = cos(pkin(6));
t635 = pkin(1) * t690;
t614 = t553 * t635;
t538 = qJD(1) * t614;
t691 = pkin(8) + qJ(3);
t629 = t691 * t550;
t613 = t547 * t629;
t506 = -qJD(1) * t613 + t538;
t615 = t550 * t635;
t664 = t547 * t553;
t514 = t664 * t691 + t615;
t507 = t514 * qJD(1);
t625 = t689 * t507;
t461 = t506 * t546 + t625;
t549 = sin(qJ(4));
t552 = cos(qJ(4));
t715 = t461 - t711 * (pkin(4) * t549 - pkin(10) * t552);
t703 = qJD(4) * t711;
t568 = t552 * t703;
t586 = t546 * t553 + t550 * t689;
t519 = t586 * t648;
t623 = t690 * qJD(1);
t596 = t623 + qJD(2);
t487 = t549 * t519 - t552 * t596;
t639 = qJD(1) * qJD(2);
t627 = t547 * t639;
t612 = t550 * t627;
t578 = qJD(2) * t532 - t546 * t612;
t566 = t552 * t578;
t555 = -qJD(4) * t487 + t566;
t714 = qJD(5) * t711 + t555;
t579 = qJD(2) * t586;
t517 = t547 * t579;
t562 = qJD(1) * t517;
t713 = -t549 * t562 - t568;
t548 = sin(qJ(5));
t551 = cos(qJ(5));
t662 = t551 * t552;
t478 = t519 * t548 - t591 * t662;
t644 = qJD(5) * t548;
t604 = t549 * t644 + t478;
t645 = qJD(4) * t552;
t631 = t551 * t645;
t712 = t604 - t631;
t486 = qJD(5) + t487;
t710 = t715 * t551;
t543 = t547 ^ 2;
t709 = -0.2e1 * t543 * t639;
t647 = qJD(2) * t547;
t632 = t550 * t647;
t708 = pkin(2) * t632;
t707 = MDP(5) * (t550 ^ 2 - t553 ^ 2);
t706 = MDP(6) * t553;
t665 = t547 * t550;
t521 = t546 * t665 - t611;
t571 = pkin(2) * t690 - t613;
t505 = t614 + t571;
t470 = t546 * t505 + t689 * t514;
t460 = pkin(9) * t690 + t470;
t522 = t586 * t547;
t610 = (-pkin(2) * t553 - pkin(1)) * t547;
t476 = pkin(3) * t521 - pkin(9) * t522 + t610;
t653 = t552 * t460 + t549 * t476;
t411 = pkin(10) * t521 + t653;
t469 = t505 * t689 - t546 * t514;
t459 = -pkin(3) * t690 - t469;
t500 = t522 * t552 + t549 * t690;
t587 = -t522 * t549 + t552 * t690;
t415 = -pkin(4) * t587 - t500 * pkin(10) + t459;
t705 = t551 * t411 + t548 * t415;
t634 = t689 * pkin(2);
t542 = -t634 - pkin(3);
t528 = -t552 * pkin(4) - t549 * pkin(10) + t542;
t694 = pkin(2) * t546;
t541 = pkin(9) + t694;
t650 = t548 * t528 + t541 * t662;
t496 = t546 * t507;
t462 = t506 * t689 - t496;
t472 = pkin(2) * t633 + pkin(3) * t519 + pkin(9) * t591;
t652 = t552 * t462 + t549 * t472;
t406 = pkin(10) * t519 + t652;
t643 = qJD(5) * t551;
t704 = t551 * t406 - t528 * t643 + t548 * t715;
t564 = -t552 * t519 - t596 * t549;
t698 = qJD(4) * t564 - t549 * t578;
t702 = t698 * MDP(24);
t701 = t486 * MDP(24);
t490 = qJD(2) * pkin(2) + qJD(1) * t571 + t538;
t440 = t546 * t490 + t625;
t429 = pkin(9) * t596 + t440;
t592 = qJD(1) * t610;
t526 = qJD(3) + t592;
t457 = pkin(3) * t591 - pkin(9) * t519 + t526;
t403 = -t549 * t429 + t552 * t457;
t393 = -pkin(4) * t711 - t403;
t445 = -t548 * t564 - t551 * t711;
t447 = t548 * t711 - t551 * t564;
t373 = t445 * pkin(5) - t447 * qJ(6) + t393;
t692 = pkin(10) * t698;
t700 = t373 * t486 + t692;
t672 = t591 * t552;
t477 = -t551 * t519 - t548 * t672;
t605 = t548 * t645 - t477;
t576 = t549 * t643 + t605;
t396 = -t548 * t562 - t551 * t714 - t564 * t644;
t682 = t396 * t548;
t699 = -t447 * t576 + t549 * t682;
t539 = qJD(2) * t614;
t569 = (-qJD(2) * t629 + qJD(3) * t553) * t547;
t493 = t539 + t569;
t494 = -qJD(2) * t514 - qJD(3) * t665;
t431 = t493 * t689 + t546 * t494;
t666 = t546 * t550;
t518 = (t624 - t666) * t647;
t473 = pkin(3) * t517 - pkin(9) * t518 + t708;
t646 = qJD(4) * t549;
t583 = t552 * t431 - t460 * t646 + t549 * t473 + t476 * t645;
t379 = pkin(10) * t517 + t583;
t430 = t493 * t546 - t689 * t494;
t467 = qJD(4) * t500 + t518 * t549;
t468 = qJD(4) * t587 + t518 * t552;
t391 = pkin(4) * t467 - pkin(10) * t468 + t430;
t697 = -qJD(5) * t705 - t379 * t548 + t391 * t551;
t696 = t447 ^ 2;
t695 = t486 ^ 2;
t554 = qJD(1) ^ 2;
t693 = pkin(5) * t698;
t688 = qJ(6) * t698;
t606 = qJD(2) * t623;
t594 = pkin(1) * t606;
t536 = t553 * t594;
t483 = qJD(1) * t569 + t536;
t557 = qJD(1) * t494;
t421 = t483 * t689 + t546 * t557;
t535 = pkin(2) * t612;
t563 = qJD(2) * t519;
t458 = pkin(3) * t563 - pkin(9) * t578 + t535;
t616 = t549 * t421 + t429 * t645 + t457 * t646 - t552 * t458;
t372 = -pkin(4) * t562 + t616;
t397 = t548 * t714 - t551 * t562 - t564 * t643;
t363 = t397 * pkin(5) + t396 * qJ(6) - t447 * qJD(6) + t372;
t687 = t363 * t548;
t686 = t363 * t551;
t404 = t552 * t429 + t549 * t457;
t394 = pkin(10) * t711 + t404;
t439 = t490 * t689 - t496;
t428 = -pkin(3) * t596 - t439;
t400 = t487 * pkin(4) + pkin(10) * t564 + t428;
t370 = t394 * t551 + t400 * t548;
t366 = qJ(6) * t486 + t370;
t685 = t366 * t486;
t684 = t370 * t486;
t683 = t372 * t548;
t681 = t396 * t552;
t680 = t397 * t552;
t679 = t698 * t548;
t678 = t698 * t551;
t677 = t445 * t486;
t676 = t447 * t445;
t675 = t447 * t486;
t674 = t564 * t519;
t673 = t591 * t549;
t671 = t519 * t487;
t670 = t528 * t551;
t669 = t541 * t548;
t668 = t541 * t551;
t667 = t543 * t554;
t663 = t549 * t551;
t630 = t541 * t644;
t661 = -qJ(6) * t673 - (-qJD(6) - t630) * t552 - (qJ(6) - t668) * t646 + t704;
t626 = pkin(5) + t669;
t660 = pkin(5) * t673 - qJD(5) * t650 + t406 * t548 + t626 * t646 - t710;
t620 = -t549 * t462 + t552 * t472;
t405 = -pkin(4) * t519 - t620;
t602 = pkin(5) * t548 - qJ(6) * t551;
t590 = t541 + t602;
t603 = pkin(5) * t551 + qJ(6) * t548;
t659 = pkin(5) * t477 - qJ(6) * t478 + t405 - (qJD(5) * t603 - qJD(6) * t551) * t549 - t590 * t645;
t658 = qJD(6) * t548 - t486 * t602 + t404;
t432 = -pkin(4) * t564 + pkin(10) * t487;
t657 = t551 * t403 + t548 * t432;
t654 = t486 * t631 - t663 * t698;
t642 = t428 * qJD(4);
t369 = -t394 * t548 + t400 * t551;
t640 = qJD(6) - t369;
t638 = pkin(1) * t667;
t637 = pkin(10) * t644;
t622 = t447 * t646 + t681;
t621 = -t549 * t460 + t476 * t552;
t420 = t546 * t483 - t689 * t557;
t619 = t486 * t551;
t618 = t543 * t550 * t553 * MDP(4);
t584 = -t552 * t421 + t429 * t646 - t457 * t645 - t549 * t458;
t371 = pkin(10) * t562 - t584;
t389 = -pkin(4) * t698 - pkin(10) * t555 + t420;
t617 = t548 * t371 - t551 * t389 + t394 * t643 + t400 * t644;
t609 = t547 * t554 * t690;
t608 = pkin(1) * t709;
t365 = -pkin(5) * t486 + t640;
t601 = t365 * t551 - t366 * t548;
t600 = t365 * t548 + t366 * t551;
t597 = -t411 * t548 + t415 * t551;
t475 = t500 * t551 + t521 * t548;
t474 = t500 * t548 - t521 * t551;
t593 = -t549 * t431 - t460 * t645 + t473 * t552 - t476 * t646;
t410 = -pkin(4) * t521 - t621;
t589 = -t486 * t643 + t679;
t588 = t393 * t486 + t692;
t585 = t373 * t447 + t617;
t359 = t551 * t371 + t548 * t389 - t394 * t644 + t400 * t643;
t582 = t551 * t379 + t548 * t391 - t411 * t644 + t415 * t643;
t380 = -pkin(4) * t517 - t593;
t580 = -pkin(8) * t664 - t615;
t577 = t591 * t711;
t570 = -t397 * t663 + t445 * t712;
t560 = t552 * t562 + (-t703 - t577) * t549;
t357 = qJD(6) * t486 + t359 - t688;
t358 = t617 + t693;
t559 = qJD(5) * t601 + t357 * t551 + t358 * t548;
t558 = t596 * t580;
t534 = -pkin(4) - t603;
t513 = t590 * t549;
t492 = t552 * t626 - t670;
t491 = -qJ(6) * t552 + t650;
t412 = pkin(5) * t447 + qJ(6) * t445;
t409 = -qJD(5) * t474 + t468 * t551 + t517 * t548;
t408 = qJD(5) * t475 + t468 * t548 - t517 * t551;
t385 = pkin(5) * t474 - qJ(6) * t475 + t410;
t383 = pkin(5) * t564 + t403 * t548 - t432 * t551;
t382 = -qJ(6) * t564 + t657;
t381 = -t396 + t677;
t375 = pkin(5) * t587 - t597;
t374 = -qJ(6) * t587 + t705;
t364 = pkin(5) * t408 - qJ(6) * t409 - qJD(6) * t475 + t380;
t362 = -pkin(5) * t467 - t697;
t361 = qJ(6) * t467 - qJD(6) * t587 + t582;
t1 = [(t553 * t608 - (-pkin(8) * t632 + t539) * t596 - (-pkin(8) * t612 + t536) * t690) * MDP(10) + t707 * t709 + ((qJD(4) - t532) * t517 + (t517 * t666 + t521 * t579) * t648) * MDP(17) + 0.2e1 * t618 * t639 + (t357 * t374 + t358 * t375 + t361 * t366 + t362 * t365 + t363 * t385 + t364 * t373) * MDP(30) + (qJD(2) * t558 + t550 * t608 + t580 * t606) * MDP(9) + (t420 * t522 - t421 * t521 + t430 * t519 - t431 * t591 - t439 * t518 - t440 * t517 - t469 * t578 - t470 * t563) * MDP(11) + (-t420 * t469 + t421 * t470 - t439 * t430 + t440 * t431 + (t526 + t592) * t708) * MDP(12) + (-MDP(7) * t632 + t647 * t706) * (0.2e1 * t623 + qJD(2)) + (-t468 * t564 + t500 * t555) * MDP(13) + (-t404 * t517 + t420 * t500 + t428 * t468 - t430 * t564 + t459 * t555 + t521 * t584 - t563 * t653 - t583 * t711) * MDP(19) + (t468 * t711 + t500 * t562 - t517 * t564 + t521 * t555) * MDP(15) + (t403 * t517 - t420 * t587 + t428 * t467 + t430 * t487 - t459 * t698 - t521 * t616 + t562 * t621 + t593 * t711) * MDP(18) + (-t467 * t711 - t487 * t517 + t521 * t698 + t562 * t587) * MDP(16) + (t397 * t587 - t408 * t486 - t445 * t467 + t474 * t698) * MDP(23) + (t467 * t486 + t587 * t698) * MDP(24) + (t396 * t587 + t409 * t486 + t447 * t467 - t475 * t698) * MDP(22) + (-t357 * t587 + t361 * t486 - t363 * t475 - t364 * t447 + t366 * t467 - t373 * t409 - t374 * t698 + t385 * t396) * MDP(29) + (t358 * t587 - t362 * t486 + t363 * t474 + t364 * t445 - t365 * t467 + t373 * t408 + t375 * t698 + t385 * t397) * MDP(27) + (t369 * t467 + t372 * t474 + t380 * t445 + t393 * t408 + t410 * t397 + t486 * t697 + t587 * t617 - t597 * t698) * MDP(25) + (t467 * t564 - t468 * t487 + t500 * t698 + t555 * t587) * MDP(14) + (t359 * t587 - t370 * t467 + t372 * t475 + t380 * t447 + t393 * t409 - t410 * t396 - t486 * t582 + t698 * t705) * MDP(26) + (-t396 * t475 + t409 * t447) * MDP(20) + (-t357 * t474 + t358 * t475 - t361 * t445 + t362 * t447 + t365 * t409 - t366 * t408 - t374 * t397 - t375 * t396) * MDP(28) + (t396 * t474 - t397 * t475 - t408 * t447 - t409 * t445) * MDP(21); (-pkin(8) * t553 * t627 - qJD(1) * t558) * MDP(9) + ((-t594 + t638) * MDP(9) + MDP(7) * t609) * t550 - t609 * t706 + (t439 * t461 - t440 * t462 + (-t420 * t689 + t421 * t546 - t526 * t633) * pkin(2)) * MDP(12) + t667 * t707 + (t560 + t671) * MDP(16) + (t447 * t673 - t486 * t604 + t622 + t654) * MDP(22) + (-t396 * t663 - t447 * t712) * MDP(20) + (-t698 * t670 - t393 * t477 - t405 * t445 + (-t710 + (-qJD(5) * t528 + t406) * t548) * t486) * MDP(25) + (-t536 + t553 * t638 + (-pkin(8) * t633 + t538) * t623 + t538 * qJD(2)) * MDP(10) + (t674 - t713) * MDP(15) + (-t605 * t486 + t680) * MDP(23) + (-t365 * t478 + t366 * t477 - t396 * t492 - t397 * t491 + t661 * t445 - t660 * t447 + t601 * t645) * MDP(28) - t554 * t618 + (-t563 * t694 - t578 * t634 + (-t439 + t462) * t591) * MDP(11) + (t570 + t699) * MDP(21) + ((t646 + t673) * t564 + (-t645 - t672) * t487) * MDP(14) + (t428 * t672 + t461 * t564 + t542 * t555 + t652 * t711) * MDP(19) + (t428 * t673 - t461 * t487 + t541 * t713 - t542 * t698 - t620 * t711) * MDP(18) - (-t404 * MDP(19) + t403 * MDP(18) + (-t440 + t461) * MDP(11) + qJD(4) * t549 ^ 2 * MDP(13) + t711 * MDP(17)) * t519 + (t698 * MDP(14) + (t393 * t643 + t369 * t591 + t683 + t541 * t397 + (t486 * t669 + t369) * qJD(4)) * MDP(25) + (t366 * t711 + t373 * t644 - t686) * MDP(29) + (-t365 * t711 + t373 * t643 + t687) * MDP(27) + (-t393 * t644 - t370 * t591 + t372 * t551 - t541 * t396 + (t486 * t668 - t370) * qJD(4)) * MDP(26) + t711 * t701 + (-t445 * t711 + t589) * MDP(23) + (-qJD(5) * t600 - t357 * t548 + t358 * t551) * MDP(28) + (t541 * t703 + t420) * MDP(19) + t642 * MDP(18)) * t549 + ((-t541 * t563 + t642) * MDP(19) - t420 * MDP(18) + t555 * MDP(14) + t577 * MDP(15) + (t393 * t548 * qJD(4) + t617 + (qJD(4) * t445 + t589) * t541) * MDP(25) - t357 * MDP(29) + t358 * MDP(27) + (t486 * t630 + t359 + (t393 * t551 + t447 * t541) * qJD(4)) * MDP(26) + t702 + ((qJD(4) * t596 + t578) * t549 - t711 * t564) * MDP(13)) * t552 + (-t393 * t478 - t405 * t447 + t486 * t704 + t650 * t698) * MDP(26) + (t396 * t513 - t698 * t491 - t661 * t486 + t659 * t447 + (t478 - t631) * t373) * MDP(29) + (t605 * t373 + t397 * t513 - t659 * t445 + t660 * t486 + t492 * t698) * MDP(27) + (t357 * t491 + t358 * t492 + t363 * t513 - t365 * t660 - t366 * t661 - t373 * t659) * MDP(30); (-t519 ^ 2 - t591 ^ 2) * MDP(11) + (t439 * t519 + t440 * t591 + t535) * MDP(12) + (t560 - t671) * MDP(18) + (-t549 * t563 - t672 * t711 - t568 + t674) * MDP(19) + ((t447 * t591 + t678) * t549 + t712 * t486 + t622) * MDP(26) + (t570 - t699) * MDP(28) + (-t681 - t478 * t486 + (-t447 * t711 - t486 * t644) * t549 + t654) * MDP(29) + (-t365 * t477 - t366 * t478 + (qJD(4) * t600 - t363) * t552 + (t373 * t711 + t559) * t549) * MDP(30) + (MDP(25) + MDP(27)) * (-t680 + t445 * t646 + (t445 * t591 + t679) * t549 - t576 * t486); -t487 ^ 2 * MDP(14) + (t487 * t591 + t566) * MDP(15) + t698 * MDP(16) + MDP(17) * t562 + (t404 * t711 - t616) * MDP(18) + (t403 * t711 + t428 * t487 + t584) * MDP(19) + (t447 * t619 - t682) * MDP(20) + ((-t396 - t677) * t551 + (-t397 - t675) * t548) * MDP(21) + (t486 * t619 - t679) * MDP(22) + (-t548 * t695 - t678) * MDP(23) + (-pkin(4) * t397 - t404 * t445 + (-t372 + (-pkin(10) * qJD(5) - t432) * t486) * t551 + (t403 * t486 + t588) * t548) * MDP(25) + (pkin(4) * t396 + t683 - t404 * t447 + (t637 + t657) * t486 + t588 * t551) * MDP(26) + (-t686 + t397 * t534 + (-pkin(10) * t643 + t383) * t486 - t658 * t445 + t700 * t548) * MDP(27) + (t382 * t445 - t383 * t447 + (t357 + t486 * t365 + (qJD(5) * t447 - t397) * pkin(10)) * t551 + (t358 - t685 + (qJD(5) * t445 - t396) * pkin(10)) * t548) * MDP(28) + (-t687 + t396 * t534 + (-t382 - t637) * t486 + t658 * t447 - t700 * t551) * MDP(29) + (pkin(10) * t559 + t363 * t534 - t365 * t383 - t366 * t382 - t373 * t658) * MDP(30) - (t487 * MDP(13) - MDP(14) * t564 + MDP(16) * t711 - t428 * MDP(18) - t447 * MDP(22) + t445 * MDP(23) - t369 * MDP(25) + t370 * MDP(26) + t365 * MDP(27) - t366 * MDP(29) - t701) * t564; MDP(20) * t676 + (-t445 ^ 2 + t696) * MDP(21) + t381 * MDP(22) + (-t397 + t675) * MDP(23) - t702 + (-t393 * t447 - t617 + t684) * MDP(25) + (t369 * t486 + t393 * t445 - t359) * MDP(26) + (-t412 * t445 - t585 + t684 - 0.2e1 * t693) * MDP(27) + (pkin(5) * t396 - qJ(6) * t397 + (t366 - t370) * t447 + (t365 - t640) * t445) * MDP(28) + (-0.2e1 * t688 - t373 * t445 + t412 * t447 + (0.2e1 * qJD(6) - t369) * t486 + t359) * MDP(29) + (-pkin(5) * t358 + qJ(6) * t357 - t365 * t370 + t366 * t640 - t373 * t412) * MDP(30); (t676 + t698) * MDP(27) + t381 * MDP(28) + (-t695 - t696) * MDP(29) + (t585 - t685 + t693) * MDP(30);];
tauc  = t1;
