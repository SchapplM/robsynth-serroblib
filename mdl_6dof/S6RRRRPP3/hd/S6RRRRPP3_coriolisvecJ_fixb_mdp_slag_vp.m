% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% MDP [32x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRPP3_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP3_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(32,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [9x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [32 1]), ...
  'S6RRRRPP3_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [32x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:56:44
% EndTime: 2019-03-09 20:56:57
% DurationCPUTime: 7.65s
% Computational Cost: add. (7992->552), mult. (19043->652), div. (0->0), fcn. (13108->6), ass. (0->230)
t550 = cos(qJ(4));
t628 = qJD(4) * t550;
t549 = sin(qJ(2));
t551 = cos(qJ(2));
t673 = sin(qJ(3));
t614 = qJD(1) * t673;
t674 = cos(qJ(3));
t615 = qJD(1) * t674;
t499 = t549 * t614 - t551 * t615;
t648 = t499 * t550;
t689 = t628 + t648;
t548 = sin(qJ(4));
t629 = qJD(4) * t548;
t649 = t499 * t548;
t688 = t629 + t649;
t571 = -t549 * t674 - t551 * t673;
t623 = qJD(2) + qJD(3);
t466 = t623 * t571;
t691 = t466 * qJD(1);
t669 = pkin(4) + qJ(6);
t690 = t669 * t691;
t675 = -pkin(8) - pkin(7);
t529 = t675 * t551;
t520 = qJD(1) * t529;
t505 = t674 * t520;
t527 = t675 * t549;
t518 = qJD(1) * t527;
t462 = t673 * t518 - t505;
t612 = qJD(3) * t673;
t587 = pkin(2) * t612 - t462;
t515 = t549 * t673 - t551 * t674;
t569 = t515 * qJD(3);
t465 = -qJD(2) * t515 - t569;
t655 = t465 * t548;
t579 = -t571 * t628 + t655;
t501 = -t549 * t615 - t551 * t614;
t572 = t550 * t501 - t548 * t623;
t687 = qJD(4) * t572;
t539 = -pkin(2) * t551 - pkin(1);
t525 = t539 * qJD(1);
t440 = pkin(3) * t499 + pkin(9) * t501 + t525;
t668 = qJD(2) * pkin(2);
t507 = t518 + t668;
t460 = t673 * t507 - t505;
t444 = pkin(9) * t623 + t460;
t402 = -t550 * t440 + t444 * t548;
t584 = -pkin(5) * t572 + t402;
t627 = qJD(5) + t584;
t624 = qJD(1) * qJD(2);
t686 = -0.2e1 * t624;
t613 = qJD(3) * t674;
t602 = pkin(2) * t613;
t494 = qJD(4) + t499;
t656 = t691 * t550;
t685 = pkin(9) * (t494 * t629 + t656);
t657 = t691 * t548;
t684 = pkin(9) * (t494 * t628 - t657);
t683 = MDP(5) * (t549 ^ 2 - t551 ^ 2);
t597 = -pkin(4) * t649 + qJ(5) * t648;
t682 = t597 - t587;
t449 = t691 * qJ(5);
t485 = t494 * qJD(5);
t681 = -t449 + t485;
t604 = pkin(4) * t629 - qJD(5) * t548;
t680 = -(-qJ(5) * qJD(4) - qJD(6)) * t550 - t604 - t688 * qJ(6);
t475 = -t674 * t527 - t673 * t529;
t678 = -pkin(5) * t689 + t501 * t669;
t626 = qJD(5) + t402;
t392 = -pkin(4) * t494 + t626;
t403 = t548 * t440 + t550 * t444;
t393 = -qJ(5) * t494 - t403;
t677 = -t392 * t548 + t393 * t550;
t556 = t465 * qJD(1);
t421 = t548 * t556 - t687;
t471 = -t501 * t548 - t550 * t623;
t676 = t471 ^ 2;
t470 = t572 ^ 2;
t491 = t494 ^ 2;
t672 = pkin(4) * t691;
t671 = pkin(4) * t501;
t670 = pkin(5) * t471;
t667 = qJ(5) * t421;
t666 = qJ(5) * t471;
t407 = -t691 * pkin(3) + (pkin(9) * t569 + (t549 * pkin(2) + pkin(9) * t515) * qJD(2)) * qJD(1);
t619 = qJD(2) * t675;
t599 = qJD(1) * t619;
t508 = t549 * t599;
t509 = t551 * t599;
t415 = t507 * t613 + t508 * t674 + t509 * t673 + t520 * t612;
t601 = t548 * t407 + t550 * t415 + t440 * t628 - t444 * t629;
t359 = -t601 - t681;
t665 = t359 * t550;
t600 = t550 * t407 - t548 * t415 - t440 * t629 - t444 * t628;
t363 = -t600 + t672;
t362 = t363 * t548;
t504 = t673 * t520;
t459 = t507 * t674 + t504;
t443 = -pkin(3) * t623 - t459;
t560 = qJ(5) * t572 + t443;
t381 = t471 * t669 + t560;
t664 = t381 * t572;
t397 = t471 * pkin(4) + t560;
t661 = t397 * t572;
t603 = qJD(4) * t623;
t420 = -t501 * t629 + (-t556 - t603) * t550;
t660 = t420 * t548;
t659 = t421 * t550;
t537 = pkin(2) * t673 + pkin(9);
t658 = t691 * t537;
t654 = t471 * t572;
t653 = t471 * t494;
t652 = t471 * t548;
t651 = t572 * t494;
t650 = t572 * t550;
t647 = t571 * t548;
t646 = t571 * t550;
t552 = qJD(2) ^ 2;
t645 = t549 * t552;
t644 = t551 * t552;
t553 = qJD(1) ^ 2;
t643 = t551 * t553;
t593 = t548 * t602;
t568 = t537 * t628 + t593;
t454 = -pkin(3) * t501 + pkin(9) * t499;
t630 = qJD(1) * t549;
t442 = pkin(2) * t630 + t454;
t463 = t518 * t674 + t504;
t455 = t548 * t463;
t608 = -t442 * t550 + t455;
t642 = t608 - t568 + t678;
t607 = t454 * t550 - t548 * t459;
t641 = -pkin(9) * t628 - t607 + t678;
t592 = t550 * t602;
t492 = t501 * qJ(5);
t598 = pkin(5) * t649 - t492;
t633 = t548 * t442 + t550 * t463;
t640 = t598 + t633 - t592 - (-pkin(5) - t537) * t629;
t634 = t548 * t454 + t550 * t459;
t639 = t598 + t634 - (-pkin(5) - pkin(9)) * t629;
t414 = t597 + t460;
t638 = t414 + t680;
t637 = t680 + t682;
t493 = -qJ(5) * t628 + t604;
t636 = t414 - t493;
t635 = -t493 + t682;
t458 = pkin(3) * t515 + pkin(9) * t571 + t539;
t476 = t527 * t673 - t529 * t674;
t632 = t548 * t458 + t550 * t476;
t383 = t403 - t670;
t625 = -qJD(6) - t383;
t622 = t674 * pkin(2);
t621 = t549 * t668;
t620 = t392 * t628 + t393 * t629 + t362;
t617 = t571 * t629;
t611 = t549 * t624;
t610 = -t548 * qJ(5) - pkin(3);
t609 = pkin(1) * t686;
t467 = t548 * t476;
t606 = t458 * t550 - t467;
t605 = t494 * t550;
t416 = t507 * t612 + t673 * t508 - t674 * t509 - t520 * t613;
t408 = -qJ(5) * t515 - t632;
t566 = qJ(5) * t420 + qJD(5) * t572 + t416;
t367 = pkin(4) * t421 + t566;
t596 = t367 * t550 + t392 * t501 - t397 * t649;
t595 = -t403 * t501 + t416 * t548 + t443 * t628;
t591 = -qJ(5) * t550 + qJ(6) * t548;
t590 = t392 * t550 + t393 * t548;
t589 = t443 * t499 + t658;
t588 = -qJD(4) * t397 - t658;
t524 = -t550 * pkin(4) + t610;
t419 = -pkin(3) * t466 - pkin(9) * t465 + t621;
t519 = t549 * t619;
t521 = t551 * t619;
t425 = -qJD(3) * t475 + t674 * t519 + t673 * t521;
t586 = t419 * t550 - t548 * t425 - t458 * t629 - t476 * t628;
t585 = -pkin(4) * t647 + t475;
t583 = -t367 * t548 - t393 * t501 - t397 * t648;
t582 = -t402 * t501 - t416 * t550 + t443 * t629;
t578 = -t465 * t550 - t617;
t577 = -pkin(5) * t421 + t601;
t576 = -pkin(5) * t420 - t600;
t575 = t403 * t494 + t600;
t574 = t525 * t501 - t416;
t426 = t673 * t519 - t521 * t674 + t527 * t612 - t529 * t613;
t573 = t548 * t419 + t550 * t425 + t458 * t628 - t476 * t629;
t506 = -t550 * t669 + t610;
t391 = -t420 + t653;
t570 = -t381 * t471 + t577;
t567 = -t537 * t629 + t592;
t565 = pkin(4) * t579 - qJ(5) * t617 + t426;
t357 = qJD(6) * t471 + t421 * t669 + t566;
t372 = -t494 * t669 + t627;
t564 = -t357 * t550 - t372 * t501 + t381 * t688;
t377 = qJD(6) - t393 - t670;
t563 = -t357 * t548 + t377 * t501 - t381 * t689;
t562 = t576 + t690;
t354 = -qJD(6) * t494 + t562;
t356 = t577 + t681;
t561 = t354 * t548 + t356 * t550 + t372 * t689 - t377 * t688;
t365 = qJ(5) * t466 - qJD(5) * t515 - t573;
t559 = qJD(4) * t590 + t362 - t665;
t558 = ((-t420 - t653) * t550 + (-t421 + t651) * t548) * MDP(19) + (-t572 * t605 - t660) * MDP(18) + (-t471 * t501 - t491 * t548 - t656) * MDP(21) + (t494 * t605 - t501 * t572 - t657) * MDP(20) + (t499 * t623 + t556) * MDP(13) + (-t501 * t623 + t691) * MDP(14) + (-t499 ^ 2 + t501 ^ 2) * MDP(12) + (-MDP(11) * t499 + t494 * MDP(22)) * t501;
t557 = t525 * t499 - t415;
t544 = t550 * pkin(5);
t543 = t548 * pkin(5);
t538 = -t622 - pkin(3);
t528 = pkin(9) * t550 + t544;
t526 = pkin(9) * t548 + t543;
t512 = t537 * t550 + t544;
t511 = t537 * t548 + t543;
t510 = -t622 + t524;
t484 = -t622 + t506;
t429 = -pkin(4) * t572 + t666;
t424 = qJ(5) * t646 + t585;
t410 = -t571 * t591 + t585;
t409 = -pkin(4) * t515 - t606;
t406 = -t572 * t669 + t666;
t400 = -t607 + t671;
t399 = t492 - t634;
t396 = t608 + t671;
t395 = t492 - t633;
t394 = pkin(5) * t647 - t408;
t390 = t467 + (-pkin(5) * t571 - t458) * t550 - t669 * t515;
t369 = (-qJ(5) * t465 + qJD(5) * t571) * t550 + t565;
t368 = pkin(4) * t466 - t586;
t364 = t591 * t465 - (qJD(6) * t548 + (qJ(6) * qJD(4) - qJD(5)) * t550) * t571 + t565;
t360 = -pkin(5) * t579 - t365;
t358 = -pkin(5) * t578 - qJD(6) * t515 + t466 * t669 - t586;
t1 = [(-t421 * t515 + t466 * t471 - t494 * t579 - t647 * t691) * MDP(21) + (t363 * t515 + t367 * t647 + t368 * t494 - t369 * t471 - t392 * t466 - t397 * t579 - t409 * t691 - t421 * t424) * MDP(26) + (t402 * t466 - t416 * t647 + t475 * t421 + t426 * t471 + t443 * t579 + t494 * t586 + t515 * t600 - t606 * t691) * MDP(23) + (-t354 * t515 - t357 * t647 - t358 * t494 + t364 * t471 + t372 * t466 + t381 * t579 + t390 * t691 + t410 * t421) * MDP(31) + (-t359 * t515 - t365 * t494 + t367 * t646 + t369 * t572 + t393 * t466 + t397 * t578 + t408 * t691 + t420 * t424) * MDP(27) + (t356 * t515 + t357 * t646 + t360 * t494 + t364 * t572 - t377 * t466 + t381 * t578 - t394 * t691 + t410 * t420) * MDP(30) + (-t420 * t515 + t466 * t572 - t494 * t578 + t646 * t691) * MDP(20) + (t403 * t466 - t416 * t646 - t475 * t420 - t426 * t572 - t443 * t578 - t494 * t573 - t515 * t601 + t632 * t691) * MDP(24) + (-t465 * t499 - t501 * t466 - t515 * t556 - t571 * t691) * MDP(12) + (-t539 * t691 - t525 * t466 + (qJD(1) * t515 + t499) * t621) * MDP(16) + (-t466 * t494 - t515 * t691) * MDP(22) + (-t501 * t465 - t556 * t571) * MDP(11) + ((-t471 * t550 + t548 * t572) * t465 - (t660 - t659 + (t650 + t652) * qJD(4)) * t571) * MDP(19) + (-t358 * t572 - t360 * t471 - t390 * t420 - t394 * t421 + (t372 * t550 - t377 * t548) * t465 - (t354 * t550 - t356 * t548 + (-t372 * t548 - t377 * t550) * qJD(4)) * t571) * MDP(29) + (t365 * t471 - t368 * t572 + t408 * t421 - t409 * t420 + t590 * t465 - (qJD(4) * t677 + t359 * t548 + t363 * t550) * t571) * MDP(25) + (MDP(13) * t465 + t466 * MDP(14) - MDP(16) * t426 - MDP(17) * t425) * t623 - MDP(7) * t645 + (pkin(7) * t645 + t551 * t609) * MDP(10) + (-pkin(7) * t644 + t549 * t609) * MDP(9) + (t420 * t646 + t572 * t578) * MDP(18) + (-pkin(2) * t571 * t611 + t525 * t465 - t501 * t621 + t539 * t556) * MDP(17) + t683 * t686 + MDP(6) * t644 + (t359 * t408 + t363 * t409 + t365 * t393 + t367 * t424 + t368 * t392 + t369 * t397) * MDP(28) + (t354 * t390 + t356 * t394 + t357 * t410 + t358 * t372 + t360 * t377 + t364 * t381) * MDP(32) + 0.2e1 * t551 * MDP(4) * t611; t553 * t683 + (-t395 * t471 + t396 * t572 + (-t572 * t602 + t393 * t499 + (qJD(4) * t471 - t420) * t537) * t548 + (-t471 * t602 + t392 * t499 - t359 + (-t421 - t687) * t537) * t550 + t620) * MDP(25) + (-t420 * t511 - t421 * t512 + t471 * t640 + t572 * t642 + t561) * MDP(29) + (t420 * t484 - t494 * t640 - t512 * t691 - t572 * t637 + t563) * MDP(30) + (-t510 * t421 + t588 * t548 + t635 * t471 + (-t396 + t568) * t494 + t596) * MDP(26) + (t510 * t420 + t588 * t550 - t635 * t572 + (t395 + t567) * t494 + t583) * MDP(27) + (t462 * t623 + (-t499 * t630 - t612 * t623) * pkin(2) + t574) * MDP(16) + (t538 * t421 + t589 * t548 + t587 * t471 + (-t593 + t455 + (-qJD(4) * t537 - t442) * t550) * t494 + t582) * MDP(23) + (-t538 * t420 + t589 * t550 - t587 * t572 + (-t567 + t633) * t494 + t595) * MDP(24) + (t463 * t623 + (t501 * t630 - t613 * t623) * pkin(2) + t557) * MDP(17) - t549 * MDP(4) * t643 + t558 + (t367 * t510 - t392 * t396 - t393 * t395 - t635 * t397 + t559 * t537 - t602 * t677) * MDP(28) + (t354 * t511 + t356 * t512 + t357 * t484 - t372 * t642 - t377 * t640 - t381 * t637) * MDP(32) + (t421 * t484 - t471 * t637 + t494 * t642 + t511 * t691 + t564) * MDP(31) + (MDP(9) * t549 * t553 + MDP(10) * t643) * pkin(1); (pkin(9) * t559 + t367 * t524 - t392 * t400 - t393 * t399 - t397 * t636) * MDP(28) + (t354 * t526 + t356 * t528 + t357 * t506 - t372 * t641 - t377 * t639 - t381 * t638) * MDP(32) + (-t420 * t526 - t421 * t528 + t471 * t639 + t572 * t641 + t561) * MDP(29) + (t420 * t506 - t494 * t639 - t528 * t691 - t572 * t638 + t563) * MDP(30) + (-t397 * t629 - t400 * t494 - t421 * t524 + t471 * t636 + t596 + t684) * MDP(26) + (-t397 * t628 + t399 * t494 + t420 * t524 - t572 * t636 + t583 - t685) * MDP(27) + (-t665 - t399 * t471 + t400 * t572 + t590 * t499 + (-t660 - t659 + (-t650 + t652) * qJD(4)) * pkin(9) + t620) * MDP(25) + (t421 * t506 - t471 * t638 + t494 * t641 + t526 * t691 + t564) * MDP(31) + (t460 * t623 + t574) * MDP(16) + (t459 * t623 + t557) * MDP(17) + t558 + (pkin(3) * t420 + t443 * t648 + t460 * t572 + t494 * t634 + t595 + t685) * MDP(24) + (-pkin(3) * t421 + t443 * t649 - t460 * t471 - t494 * t607 + t582 - t684) * MDP(23); -MDP(18) * t654 + (t470 - t676) * MDP(19) + t391 * MDP(20) + (-t421 - t651) * MDP(21) - t691 * MDP(22) + (t443 * t572 + t575) * MDP(23) + (-t402 * t494 + t443 * t471 - t601) * MDP(24) + (pkin(4) * t420 - t667 - (-t393 - t403) * t572 + (t392 - t626) * t471) * MDP(25) + (t429 * t471 - t575 - t661 + 0.2e1 * t672) * MDP(26) + (-t397 * t471 - t429 * t572 + t494 * t626 - t359 - t449) * MDP(27) + (-pkin(4) * t363 - qJ(5) * t359 - t392 * t403 - t393 * t626 - t397 * t429) * MDP(28) + (-t667 + t420 * t669 - (t377 + t625) * t572 + (t372 - t627) * t471) * MDP(29) + (-t406 * t572 + t494 * t584 - 0.2e1 * t449 + 0.2e1 * t485 + t570) * MDP(30) + (t664 - t406 * t471 + (0.2e1 * qJD(6) + t383) * t494 - 0.2e1 * t690 - t576) * MDP(31) + (qJ(5) * t356 - t354 * t669 + t372 * t625 + t377 * t627 - t381 * t406) * MDP(32); (t393 * t494 + t363 - t661) * MDP(28) + (-t664 + (-qJD(6) - t377) * t494 + t562) * MDP(32) + (MDP(25) + MDP(29)) * t391 + (MDP(27) + MDP(30)) * (-t470 - t491) + (-MDP(26) + MDP(31)) * (t691 - t654); (t501 * t628 - t548 * t603 - t651) * MDP(29) - MDP(30) * t654 + (-t491 - t676) * MDP(31) + (t372 * t494 + t570 + t681) * MDP(32) + (-MDP(29) * t655 - t466 * MDP(30)) * qJD(1);];
tauc  = t1;
