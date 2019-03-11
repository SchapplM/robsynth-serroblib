% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RPRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% MDP [35x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RPRRRR5_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR5_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(35,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [11x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [35 1]), ...
  'S6RPRRRR5_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [35x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:10:05
% EndTime: 2019-03-09 07:10:19
% DurationCPUTime: 9.63s
% Computational Cost: add. (9228->490), mult. (24033->638), div. (0->0), fcn. (19856->10), ass. (0->217)
t572 = cos(pkin(11));
t579 = cos(qJ(3));
t661 = t572 * t579;
t571 = sin(pkin(11));
t576 = sin(qJ(3));
t662 = t571 * t576;
t596 = -t661 + t662;
t534 = t596 * qJD(1);
t542 = t571 * t579 + t572 * t576;
t535 = t542 * qJD(1);
t575 = sin(qJ(4));
t688 = cos(qJ(4));
t512 = t688 * t534 + t535 * t575;
t693 = qJD(5) + qJD(6);
t724 = t512 + t693;
t573 = sin(qJ(6));
t574 = sin(qJ(5));
t577 = cos(qJ(6));
t578 = cos(qJ(5));
t545 = t573 * t574 - t577 * t578;
t721 = t724 * t545;
t632 = qJD(1) * qJD(3);
t620 = t579 * t632;
t558 = t572 * t620;
t621 = t576 * t632;
t531 = -t571 * t621 + t558;
t593 = -t575 * t534 + t535 * t688;
t647 = t571 * t620 + t572 * t621;
t469 = qJD(4) * t593 + t575 * t531 + t688 * t647;
t703 = -qJD(5) - t512;
t507 = qJD(6) - t703;
t660 = t573 * t578;
t546 = t574 * t577 + t660;
t713 = t546 * t469 - t507 * t721;
t723 = t724 * t546;
t639 = qJD(5) * t578;
t707 = t512 * t578;
t722 = t639 + t707;
t708 = t512 * t574;
t720 = pkin(10) * t708;
t686 = pkin(7) + qJ(2);
t551 = t686 * t571;
t543 = qJD(1) * t551;
t552 = t686 * t572;
t544 = qJD(1) * t552;
t598 = t543 * t576 - t544 * t579;
t499 = -pkin(8) * t534 - t598;
t493 = t575 * t499;
t697 = -t579 * t543 - t544 * t576;
t498 = -pkin(8) * t535 + t697;
t442 = t498 * t688 - t493;
t622 = qJD(4) * t688;
t704 = pkin(3) * t622 - t442;
t640 = qJD(5) * t574;
t719 = (t640 + t708) * pkin(5);
t718 = pkin(5) * t593 + pkin(10) * t707;
t480 = -t647 * pkin(8) - qJD(2) * t534 + qJD(3) * t697;
t587 = t542 * qJD(2);
t586 = qJD(1) * t587;
t481 = -pkin(8) * t531 + qJD(3) * t598 - t586;
t495 = qJD(3) * pkin(3) + t498;
t641 = qJD(4) * t575;
t396 = t575 * t480 - t688 * t481 + t495 * t641 + t499 * t622;
t570 = qJD(3) + qJD(4);
t502 = t570 * t574 + t578 * t593;
t468 = t688 * t531 - t534 * t622 - t535 * t641 - t575 * t647;
t679 = t468 * t574;
t424 = qJD(5) * t502 + t679;
t381 = pkin(5) * t424 + t396;
t438 = t495 * t688 - t493;
t432 = -t570 * pkin(4) - t438;
t500 = -t578 * t570 + t574 * t593;
t414 = t500 * pkin(5) + t432;
t717 = t381 * t546 - t414 * t721;
t716 = t381 * t545 + t414 * t723;
t473 = pkin(4) * t593 + pkin(9) * t512;
t450 = pkin(3) * t535 + t473;
t715 = -t578 * t450 - t574 * t704;
t423 = t578 * t468 + t570 * t639 - t593 * t640;
t714 = t423 * t578 - t574 * t424 - t500 * t722;
t606 = -t545 * t469 - t507 * t723;
t464 = t574 * t469;
t649 = -t639 * t703 + t464;
t712 = -t703 * t707 + t649;
t637 = qJD(6) * t577;
t629 = t577 * t423 - t573 * t424 - t500 * t637;
t638 = qJD(6) * t573;
t382 = -t502 * t638 + t629;
t599 = t500 * t573 - t577 * t502;
t615 = t423 * t573 + t577 * t424;
t383 = -qJD(6) * t599 + t615;
t421 = t423 * t574;
t674 = t502 * t573;
t444 = t577 * t500 + t674;
t711 = (t502 * t722 + t421) * MDP(22) + t382 * t546 * MDP(29) + (-t382 * t545 - t546 * t383 + t444 * t721) * MDP(30) + (t721 * MDP(29) + MDP(30) * t723) * t599;
t710 = t444 * t507;
t709 = t507 * t599;
t494 = t688 * t499;
t439 = t575 * t495 + t494;
t433 = pkin(9) * t570 + t439;
t562 = -pkin(2) * t572 - pkin(1);
t550 = qJD(1) * t562 + qJD(2);
t520 = pkin(3) * t534 + t550;
t440 = pkin(4) * t512 - pkin(9) * t593 + t520;
t401 = t433 * t578 + t440 * t574;
t391 = -pkin(10) * t500 + t401;
t389 = t391 * t638;
t702 = t414 * t444 + t389;
t410 = pkin(3) * t647 + t469 * pkin(4) - t468 * pkin(9);
t407 = t578 * t410;
t582 = -t480 * t688 - t575 * t481 - t495 * t622 + t499 * t641;
t584 = -qJD(5) * t401 + t574 * t582 + t407;
t370 = pkin(5) * t469 - pkin(10) * t423 + t584;
t589 = t574 * t410 - t433 * t640 + t440 * t639 - t578 * t582;
t371 = -pkin(10) * t424 + t589;
t617 = t577 * t370 - t573 * t371;
t701 = t414 * t599 + t617;
t700 = t469 * MDP(33) + (-t444 ^ 2 + t599 ^ 2) * MDP(30) - t444 * MDP(29) * t599;
t517 = t542 * t688 - t575 * t596;
t474 = t546 * t517;
t441 = t575 * t498 + t494;
t608 = pkin(3) * t641 - t441;
t696 = t574 * t450 - t578 * t704;
t466 = t578 * t469;
t695 = -t640 * t703 - t466;
t663 = t551 * t579;
t505 = -pkin(8) * t542 - t552 * t576 - t663;
t597 = t551 * t576 - t552 * t579;
t506 = -pkin(8) * t596 - t597;
t694 = t688 * t505 - t575 * t506;
t692 = (t571 ^ 2 + t572 ^ 2) * (MDP(7) * qJ(2) + MDP(6));
t691 = -t593 * MDP(15) + MDP(16) * t512 - t520 * MDP(21);
t400 = -t433 * t574 + t578 * t440;
t390 = -pkin(10) * t502 + t400;
t385 = -pkin(5) * t703 + t390;
t683 = t385 * t577;
t374 = -t391 * t573 + t683;
t682 = t391 * t577;
t375 = t385 * t573 + t682;
t690 = MDP(16) * t593 - MDP(20) * t520 + MDP(26) * t703 - MDP(27) * t400 + MDP(28) * t401 - MDP(33) * t507 - MDP(34) * t374 + MDP(35) * t375;
t689 = -pkin(9) - pkin(10);
t687 = t578 * pkin(5);
t564 = pkin(3) * t575 + pkin(9);
t685 = -pkin(10) - t564;
t681 = t444 * t593;
t680 = t599 * t593;
t536 = t596 * qJD(3);
t537 = t542 * qJD(3);
t592 = -t575 * t542 - t596 * t688;
t478 = qJD(4) * t592 - t536 * t688 - t575 * t537;
t678 = t478 * t574;
t677 = t478 * t578;
t676 = t500 * t593;
t675 = t502 * t593;
t673 = t502 * t574;
t672 = t593 * t570;
t669 = t512 * t570;
t666 = t517 * t574;
t665 = t517 * t578;
t458 = t575 * t505 + t506 * t688;
t454 = t578 * t458;
t430 = t432 * t639;
t656 = t396 * t574 + t430;
t654 = t578 * t438 + t574 * t473;
t526 = pkin(3) * t596 + t562;
t459 = -pkin(4) * t592 - pkin(9) * t517 + t526;
t650 = t574 * t459 + t454;
t648 = t608 + t719;
t642 = qJD(3) * t535;
t633 = qJD(1) * qJD(2);
t627 = qJD(5) * t689;
t625 = qJD(1) * t662;
t623 = t517 * t640;
t619 = qJD(5) * t685;
t616 = -t396 * t578 + t432 * t640;
t614 = -t438 * t574 + t578 * t473;
t612 = t703 * t574;
t611 = qJD(6) * t385 + t371;
t565 = -pkin(3) * t688 - pkin(4);
t607 = -t439 + t719;
t567 = t578 * pkin(10);
t540 = t564 * t578 + t567;
t605 = qJD(6) * t540 - t578 * t619 - t715 + t718;
t555 = pkin(9) * t578 + t567;
t604 = qJD(6) * t555 - t578 * t627 + t614 + t718;
t539 = t685 * t574;
t603 = -qJD(6) * t539 - t574 * t619 + t696 + t720;
t554 = t689 * t574;
t602 = -qJD(6) * t554 - t574 * t627 + t654 + t720;
t600 = t432 * t512 - t469 * t564;
t594 = t703 * t708 - t695;
t591 = t517 * t639 + t678;
t590 = -t623 + t677;
t585 = -qJD(3) * t663 + qJD(2) * t661 + (-qJD(2) * t571 - qJD(3) * t552) * t576;
t484 = -pkin(8) * t537 + t585;
t581 = qJD(3) * t597 - t587;
t485 = pkin(8) * t536 + t581;
t404 = qJD(4) * t694 + t688 * t484 + t575 * t485;
t479 = qJD(4) * t517 - t575 * t536 + t537 * t688;
t413 = pkin(3) * t537 + pkin(4) * t479 - pkin(9) * t478;
t588 = t578 * t404 + t574 * t413 - t458 * t640 + t459 * t639;
t405 = qJD(4) * t458 + t575 * t484 - t485 * t688;
t566 = -pkin(4) - t687;
t553 = t565 - t687;
t475 = t545 * t517;
t456 = t578 * t459;
t431 = t469 * t592;
t428 = pkin(5) * t666 - t694;
t412 = t578 * t413;
t403 = -pkin(10) * t666 + t650;
t397 = -pkin(5) * t592 - pkin(10) * t665 - t458 * t574 + t456;
t387 = t478 * t660 - t573 * t623 - t638 * t666 + (t665 * t693 + t678) * t577;
t386 = -t474 * t693 - t545 * t478;
t384 = pkin(5) * t591 + t405;
t373 = -pkin(10) * t591 + t588;
t372 = -pkin(10) * t677 + pkin(5) * t479 - t404 * t574 + t412 + (-t454 + (pkin(10) * t517 - t459) * t574) * qJD(5);
t1 = [(MDP(17) * t478 - MDP(18) * t479 - MDP(20) * t405 - MDP(21) * t404) * t570 + (t526 * t468 + t520 * t478 + (t517 * t647 + t537 * t593) * pkin(3)) * MDP(21) + (t468 * t517 + t478 * t593) * MDP(15) + (t396 * t665 - t401 * t479 + t405 * t502 - t423 * t694 + t432 * t590 - t469 * t650 + t588 * t703 + t589 * t592) * MDP(28) + (-(-t458 * t639 + t412) * t703 + t456 * t469 - (-t433 * t639 + t407) * t592 + t400 * t479 + t405 * t500 - t694 * t424 + t517 * t430 + (-(-qJD(5) * t459 - t404) * t703 - t458 * t469 - (-qJD(5) * t440 + t582) * t592 + t396 * t517 + t432 * t478) * t574) * MDP(27) + (t424 * t592 - t464 * t517 - t479 * t500 + t591 * t703) * MDP(25) + (-t423 * t592 + t466 * t517 + t479 * t502 - t590 * t703) * MDP(24) + (-t479 * t703 - t431) * MDP(26) + 0.2e1 * t633 * t692 + ((-t500 * t578 - t673) * t478 + (-t421 - t424 * t578 + (t500 * t574 - t502 * t578) * qJD(5)) * t517) * MDP(23) + (-t382 * t475 - t386 * t599) * MDP(29) + (-t382 * t474 + t383 * t475 - t386 * t444 + t387 * t599) * MDP(30) + (t468 * t592 - t469 * t517 - t478 * t512 - t479 * t593) * MDP(16) + (t526 * t469 + t520 * t479 + (t537 * t512 - t592 * t647) * pkin(3)) * MDP(20) + ((t372 * t577 - t373 * t573) * t507 + (t397 * t577 - t403 * t573) * t469 - t617 * t592 + t374 * t479 + t384 * t444 + t428 * t383 + t381 * t474 + t414 * t387 + ((-t397 * t573 - t403 * t577) * t507 + t375 * t592) * qJD(6)) * MDP(34) + (t383 * t592 - t387 * t507 - t444 * t479 - t469 * t474) * MDP(32) + (-t375 * t479 - t381 * t475 + t428 * t382 - t384 * t599 + t414 * t386 - t389 * t592 + (-(-qJD(6) * t403 + t372) * t507 - t397 * t469 + t370 * t592) * t573 + (-(qJD(6) * t397 + t373) * t507 - t403 * t469 + t611 * t592) * t577) * MDP(35) + (-t382 * t592 + t386 * t507 - t469 * t475 - t479 * t599) * MDP(31) + (-t531 * t596 + t536 * t534 - t535 * t537 - t542 * t647) * MDP(9) + (t423 * t665 + t502 * t590) * MDP(22) + (t550 * t537 + t562 * t647) * MDP(13) + (t479 * t507 - t431) * MDP(33) + (-t536 * MDP(10) - t537 * MDP(11) + MDP(13) * t581 - MDP(14) * t585) * qJD(3) + (t531 * t542 - t535 * t536) * MDP(8) + (t562 * t531 - t550 * t536) * MDP(14); (t642 + t647) * MDP(13) + (t558 + (-t534 - t625) * qJD(3)) * MDP(14) + (t469 + t672) * MDP(20) + (t468 - t669) * MDP(21) + (t594 - t676) * MDP(27) + (-t578 * t703 ^ 2 - t464 - t675) * MDP(28) + (t606 - t681) * MDP(34) + (t680 - t713) * MDP(35) - qJD(1) ^ 2 * t692; t711 + (t442 * t570 + (-t535 * t593 - t570 * t622) * pkin(3) + t582) * MDP(21) + ((t539 * t577 - t540 * t573) * t469 + t553 * t383 + (t573 * t603 - t577 * t605) * t507 + t648 * t444 + t716) * MDP(34) + (-(t539 * t573 + t540 * t577) * t469 + t553 * t382 + (t573 * t605 + t577 * t603) * t507 - t648 * t599 + t717) * MDP(35) + (t680 + t713) * MDP(31) + (-t469 + t672) * MDP(18) + (t565 * t423 + t600 * t578 + t608 * t502 - (t564 * t640 + t696) * t703 + t656) * MDP(28) + (t565 * t424 + t600 * t574 + t608 * t500 - (-t564 * t639 + t715) * t703 + t616) * MDP(27) + (t606 + t681) * MDP(32) + (t673 * t703 + t714) * MDP(23) + (-t675 + t712) * MDP(24) + (t594 + t676) * MDP(25) + (t642 - t647) * MDP(11) + (-t550 * t535 - t586) * MDP(13) + t535 * t534 * MDP(8) + t690 * t593 - t691 * t512 + (-t534 ^ 2 + t535 ^ 2) * MDP(9) + (t550 * t534 + t596 * t633) * MDP(14) + (t558 + (t534 - t625) * qJD(3)) * MDP(10) + (t468 + t669) * MDP(17) + (t441 * t570 + (-t512 * t535 - t570 * t641) * pkin(3) - t396) * MDP(20); t468 * MDP(17) - t469 * MDP(18) + (t439 * t570 - t396) * MDP(20) + (t438 * t570 + t582) * MDP(21) + (t502 * t612 + t714) * MDP(23) + t712 * MDP(24) + (-t612 * t703 + t466) * MDP(25) + (-pkin(4) * t424 - pkin(9) * t649 + t432 * t708 - t439 * t500 + t614 * t703 + t616) * MDP(27) + (-pkin(4) * t423 + pkin(9) * t695 + t432 * t707 - t439 * t502 - t654 * t703 + t656) * MDP(28) + t713 * MDP(31) + t606 * MDP(32) + ((t554 * t577 - t555 * t573) * t469 + t566 * t383 + (t573 * t602 - t577 * t604) * t507 + t607 * t444 + t716) * MDP(34) + (-(t554 * t573 + t555 * t577) * t469 + t566 * t382 + (t573 * t604 + t577 * t602) * t507 - t607 * t599 + t717) * MDP(35) + (t570 * MDP(17) - t691) * t512 + (MDP(18) * t570 - MDP(24) * t502 + MDP(25) * t500 + MDP(31) * t599 + MDP(32) * t444 + t690) * t593 + t711; t502 * t500 * MDP(22) + (-t500 ^ 2 + t502 ^ 2) * MDP(23) + (-t500 * t703 + t423) * MDP(24) + (-t679 + (-qJD(5) - t703) * t502) * MDP(25) + t469 * MDP(26) + (-t401 * t703 - t432 * t502 + t584) * MDP(27) + (-t400 * t703 + t432 * t500 - t589) * MDP(28) + (t382 + t710) * MDP(31) + (-t383 - t709) * MDP(32) + (-(-t390 * t573 - t682) * t507 - t375 * qJD(6) + (-t444 * t502 + t469 * t577 - t507 * t638) * pkin(5) + t701) * MDP(34) + ((-t391 * t507 - t370) * t573 + (t390 * t507 - t611) * t577 + (-t469 * t573 + t502 * t599 - t507 * t637) * pkin(5) + t702) * MDP(35) + t700; (t629 + t710) * MDP(31) + (-t615 - t709) * MDP(32) + (t375 * t507 + t701) * MDP(34) + (-t573 * t370 - t577 * t371 + t374 * t507 + t702) * MDP(35) + (-MDP(31) * t674 + MDP(32) * t599 - MDP(34) * t375 - MDP(35) * t683) * qJD(6) + t700;];
tauc  = t1;
