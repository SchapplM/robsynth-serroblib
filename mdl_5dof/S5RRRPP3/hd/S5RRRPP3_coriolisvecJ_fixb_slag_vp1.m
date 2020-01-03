% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:36
% DurationCPUTime: 23.66s
% Computational Cost: add. (13578->600), mult. (15296->688), div. (0->0), fcn. (11815->6), ass. (0->361)
t781 = Icges(5,4) - Icges(4,5);
t780 = Icges(5,5) - Icges(4,6);
t384 = sin(qJ(3));
t386 = cos(qJ(3));
t779 = t780 * t384 - t781 * t386;
t778 = Icges(5,1) + Icges(4,3);
t777 = Icges(5,4) - Icges(6,5);
t379 = Icges(6,6) * t384;
t302 = -Icges(6,2) * t386 + t379;
t629 = Icges(5,6) * t384;
t476 = Icges(5,3) * t386 + t629;
t775 = t302 - t476;
t383 = qJ(1) + qJ(2);
t378 = cos(t383);
t774 = t779 * t378;
t627 = Icges(6,6) * t386;
t475 = -Icges(6,3) * t384 + t627;
t628 = Icges(5,6) * t386;
t479 = Icges(5,2) * t384 + t628;
t773 = t475 - t479;
t377 = sin(t383);
t604 = t377 * t386;
t605 = t377 * t384;
t765 = t778 * t378 + t781 * t604 - t605 * t780;
t772 = t778 * t377 + t774;
t351 = Icges(4,4) * t605;
t182 = Icges(4,1) * t604 - Icges(4,5) * t378 - t351;
t340 = Icges(6,6) * t605;
t185 = Icges(6,5) * t378 - Icges(6,3) * t604 - t340;
t771 = -t182 + t185;
t632 = Icges(4,4) * t384;
t314 = Icges(4,1) * t386 - t632;
t445 = t314 * t378;
t183 = Icges(4,5) * t377 + t445;
t474 = Icges(6,3) * t386 + t379;
t439 = t474 * t378;
t184 = Icges(6,5) * t377 + t439;
t760 = t183 + t184;
t477 = -Icges(5,3) * t384 + t628;
t440 = t477 * t378;
t186 = Icges(5,5) * t377 - t440;
t601 = t378 * t386;
t342 = Icges(6,6) * t601;
t602 = t378 * t384;
t188 = Icges(6,4) * t377 + Icges(6,2) * t602 + t342;
t759 = t186 + t188;
t187 = Icges(5,5) * t378 + Icges(5,6) * t604 - Icges(5,3) * t605;
t341 = Icges(6,6) * t604;
t189 = Icges(6,4) * t378 - Icges(6,2) * t605 - t341;
t770 = t187 + t189;
t769 = (Icges(6,4) + Icges(5,5)) * t386 + t777 * t384;
t180 = Icges(4,4) * t604 - Icges(4,2) * t605 - Icges(4,6) * t378;
t343 = Icges(5,6) * t605;
t191 = Icges(5,4) * t378 + Icges(5,2) * t604 - t343;
t768 = t180 * t384 - t191 * t386;
t301 = Icges(6,2) * t384 + t627;
t380 = Icges(4,4) * t386;
t484 = -Icges(4,2) * t384 + t380;
t739 = t301 - t477 - t484;
t305 = Icges(4,5) * t384 + Icges(4,6) * t386;
t751 = t305 - t769;
t480 = Icges(5,2) * t386 - t629;
t767 = t314 + t474 + t480;
t766 = -t183 * t604 - t186 * t605;
t763 = t182 * t386 - t187 * t384 - t768;
t313 = Icges(4,1) * t384 + t380;
t734 = t479 + t313;
t762 = t475 - t734;
t311 = Icges(4,2) * t386 + t632;
t747 = -t302 + t311;
t761 = -t476 - t747;
t463 = t311 * t384 - t313 * t386;
t750 = t775 * t384 - t773 * t386 - t463;
t307 = Icges(6,4) * t384 + Icges(6,5) * t386;
t442 = t307 * t378;
t192 = Icges(6,1) * t377 + t442;
t494 = -t184 * t604 - t188 * t605 + t192 * t378;
t444 = t484 * t378;
t181 = Icges(4,6) * t377 + t444;
t344 = Icges(5,6) * t602;
t190 = Icges(5,4) * t377 - Icges(5,2) * t601 + t344;
t757 = -t772 * t378 - t766;
t715 = -t181 * t605 - t190 * t604 + t757;
t758 = -t494 + t715;
t756 = t769 * t378;
t755 = t759 * t602 + t760 * t601 + (t192 + t772) * t377;
t193 = Icges(6,1) * t378 - Icges(6,4) * t605 - Icges(6,5) * t604;
t170 = t377 * t193;
t754 = t765 * t377 + t771 * t601 + t770 * t602 + t170;
t753 = t767 * qJD(3);
t752 = t739 * qJD(3);
t733 = t307 + t779;
t673 = t751 * t377;
t382 = qJD(1) + qJD(2);
t749 = t384 * ((Icges(4,5) - t777) * t382 + t762 * qJD(3)) + t386 * ((-Icges(6,4) - t780) * t382 + (-t311 + t775) * qJD(3));
t748 = t773 * t604 - t605 * t775 - t756;
t746 = t181 * t384 + t190 * t386;
t614 = t193 * t378;
t470 = t185 * t386 + t189 * t384;
t682 = t377 * t470;
t84 = t614 - t682;
t745 = t763 * t377 + t765 * t378 + t84;
t696 = -t180 * t602 + t191 * t601 - t754;
t695 = -t181 * t602 - t190 * t601 + t755;
t744 = t750 * t378 + t673;
t743 = t180 + t770;
t742 = t181 - t759;
t741 = t191 - t771;
t740 = -t190 + t760;
t646 = rSges(6,1) + pkin(4);
t335 = qJ(5) * t604;
t724 = -rSges(6,2) * t605 - rSges(6,3) * t604 - t335;
t707 = t646 * t378 + t724;
t735 = t707 * t382;
t732 = t753 * t386 + t752 * t384 + t751 * t382 + (t762 * t384 + t761 * t386) * qJD(3);
t731 = (Icges(6,1) + t778) * t382 - t751 * qJD(3);
t730 = t470 - t763;
t729 = t759 * t384 + t760 * t386 - t746;
t728 = -t733 * qJD(3) + t750 * t382;
t727 = t377 ^ 2;
t603 = t378 * t382;
t726 = ((t442 + t730 + t774) * t382 + t731 * t377) * t378;
t372 = t378 * pkin(7);
t265 = pkin(2) * t377 - t372;
t319 = pkin(7) * t603;
t725 = t382 * t265 + t319;
t723 = t744 * t382;
t722 = (t758 * t377 - t745 * t378) * qJD(3);
t721 = (t377 * t695 - t378 * t696) * qJD(3);
t720 = -t739 - t762;
t719 = t767 + t761;
t718 = Icges(5,3) * t601 + t378 * t747 + t344 - t740;
t717 = -Icges(6,3) * t602 - t378 * t734 + t342 - t742;
t610 = t305 * t378;
t145 = -t377 * t463 - t610;
t716 = (-t145 + t748) * t382;
t714 = rSges(6,3) + qJ(5);
t547 = qJD(5) * t386;
t713 = t377 * t547;
t600 = t382 * t384;
t540 = t377 * t600;
t376 = qJD(4) * t384;
t333 = t378 * t376;
t550 = qJD(3) * t386;
t525 = t378 * t550;
t570 = qJ(4) * t525 + t333;
t551 = qJD(3) * t384;
t526 = t378 * t551;
t599 = t382 * t386;
t705 = t377 * t599 + t526;
t140 = -pkin(3) * t705 - qJ(4) * t540 + t570;
t625 = qJ(4) * t384;
t644 = t386 * pkin(3);
t326 = t625 + t644;
t244 = t326 * t377;
t213 = t382 * t244;
t708 = t140 + t213;
t706 = t746 + t765;
t216 = t525 - t540;
t528 = t377 * t550;
t538 = t378 * t600;
t215 = t528 + t538;
t704 = t570 + t213 + t725;
t385 = sin(qJ(1));
t638 = pkin(1) * qJD(1);
t541 = t385 * t638;
t263 = rSges(3,1) * t377 + rSges(3,2) * t378;
t612 = t263 * t382;
t209 = -t541 - t612;
t642 = rSges(5,2) * t384;
t488 = rSges(5,3) * t386 + t642;
t703 = qJD(3) * t488;
t702 = -t716 + t722;
t701 = t721 + t723;
t606 = t377 * t382;
t700 = t729 * qJD(3) + (-t384 * t767 + t386 * t739) * t606 + t749 * t378;
t699 = t730 * qJD(3) + (t301 * t386 - t480 * t384) * t603 + ((-t440 - t444) * t386 + (-t439 - t445) * t384) * t382 - t749 * t377;
t698 = -t728 * t377 + t732 * t378;
t697 = t732 * t377 + t728 * t378;
t694 = t384 * t741 + t386 * t743;
t693 = t384 * t740 + t386 * t742;
t690 = -t614 + t755;
t689 = t731 * t378 - t729 * t382 - t733 * t606;
t687 = qJ(5) * t601 + t646 * t377;
t686 = 0.2e1 * qJD(3);
t685 = -rSges(6,2) - qJ(4);
t505 = rSges(5,1) * t603 + t705 * rSges(5,2) + rSges(5,3) * t525;
t139 = -rSges(5,3) * t540 + t505;
t389 = qJD(1) ^ 2;
t645 = pkin(1) * t385;
t544 = t389 * t645;
t498 = t382 * (-pkin(2) * t606 + t319) - t544;
t546 = qJD(3) * qJD(4);
t520 = t386 * t546;
t433 = t377 * t520 + t498 + (t140 + t333) * t382;
t321 = pkin(3) * t384 - qJ(4) * t386;
t561 = -t321 + t488;
t513 = t378 * t561;
t549 = qJD(4) * t386;
t256 = qJD(3) * t326 - t549;
t639 = rSges(5,3) * t384;
t641 = rSges(5,2) * t386;
t328 = t639 - t641;
t571 = -t328 * qJD(3) - t256;
t36 = t139 * t382 + (t571 * t377 + t382 * t513) * qJD(3) + t433;
t537 = t378 * t599;
t281 = rSges(5,2) * t537;
t137 = rSges(5,3) * t538 - t281 + (rSges(5,1) * t382 + t703) * t377;
t553 = qJD(3) * t377;
t257 = t321 * t553;
t387 = cos(qJ(1));
t381 = t387 * pkin(1);
t543 = t389 * t381;
t462 = t382 * t257 + t378 * t520 - t543;
t552 = qJD(3) * t378;
t529 = t377 * t551;
t296 = pkin(3) * t529;
t524 = t377 * t376;
t141 = pkin(3) * t537 + qJ(4) * t215 - t296 + t524;
t266 = t378 * pkin(2) + t377 * pkin(7);
t220 = t266 * t382;
t596 = -t141 - t220;
t37 = t571 * t552 + (-t137 + (-t376 - t703) * t377 + t596) * t382 + t462;
t460 = -pkin(2) - t326;
t634 = -rSges(5,3) - qJ(4);
t352 = rSges(5,3) * t605;
t207 = rSges(5,1) * t378 + rSges(5,2) * t604 - t352;
t454 = qJD(3) * t513 + t333;
t423 = t454 - t541;
t573 = -t244 - t265;
t73 = (t207 + t573) * t382 + t423;
t635 = t382 * t73;
t493 = -t257 + t524;
t458 = t488 * t553 + t493;
t555 = t377 * rSges(5,1) + rSges(5,3) * t602;
t205 = -rSges(5,2) * t601 + t555;
t249 = pkin(3) * t601 + qJ(4) * t602;
t503 = t249 + t266;
t533 = t205 + t503;
t542 = t387 * t638;
t74 = t382 * t533 + t458 + t542;
t684 = (t37 * (t460 + t641) + t74 * (t460 - t639) * t382 + (-t376 + (t386 * t634 - t642) * qJD(3) + (-rSges(5,1) - pkin(7)) * t382) * t73) * t377 + (-t74 * pkin(3) * t551 - t36 * t641 + t37 * rSges(5,1) + (t384 * t634 - pkin(2) - t644) * t635) * t378;
t492 = t244 * t553 + t249 * t552 - t549;
t548 = qJD(5) * t384;
t578 = rSges(6,2) * t602 + rSges(6,3) * t601 + t687;
t51 = t548 + (-t377 * t707 + t578 * t378) * qJD(3) + t492;
t683 = qJD(3) * t51;
t487 = rSges(6,2) * t386 - rSges(6,3) * t384;
t624 = qJ(5) * t384;
t512 = -t487 + t624;
t496 = -t321 - t512;
t455 = t378 * t496;
t288 = qJ(5) * t529;
t290 = rSges(6,3) * t529;
t680 = t288 + t290;
t679 = rSges(4,1) * t601 + t377 * rSges(4,3);
t554 = t378 ^ 2 + t727;
t388 = qJD(3) ^ 2;
t461 = t547 + t376;
t640 = rSges(6,2) * t384;
t322 = rSges(6,3) * t386 + t640;
t572 = -t322 * qJD(3) - t256;
t491 = -0.2e1 * t548 + t572;
t598 = rSges(6,2) * t528 + t713 - t680 + (t322 * t378 + t687) * t382;
t623 = qJ(5) * t386;
t22 = (qJD(3) * t491 - t388 * t623) * t378 + ((qJD(3) * t512 - t461) * t377 + t596 - t598) * t382 + t462;
t435 = qJD(3) * t455;
t545 = -pkin(3) - t714;
t502 = t545 * t386;
t514 = -pkin(2) - t625;
t331 = t378 * t547;
t559 = t331 + t333;
t64 = -t541 + t435 + (t573 + t707) * t382 + t559;
t636 = t382 * t64;
t260 = t487 * t553;
t534 = -t249 - t578;
t506 = t266 - t534;
t65 = t377 * t461 + t382 * t506 - t257 + t260 - t288 + t542;
t678 = (t22 * t646 + t65 * t545 * t551 + (t384 * t685 - pkin(2) + t502) * t636) * t378 + (t22 * t514 - t64 * t376 + (-t22 * pkin(3) + t64 * (qJD(3) * t685 - qJD(5))) * t386 + (t64 * (-pkin(7) - t646) + (t514 - t640 + t502) * t65) * t382) * t377 - t435 * t65;
t675 = rSges(6,2) * t525 + t646 * t603 + t331;
t677 = -t559 + t675 + t704 - t735;
t676 = -qJ(5) * t550 + t572;
t556 = rSges(4,2) * t605 + t378 * rSges(4,3);
t202 = rSges(4,1) * t604 - t556;
t173 = t382 * t202;
t446 = -t216 * rSges(4,2) + rSges(4,3) * t603;
t674 = -rSges(4,1) * t526 + t173 + t446 + t725;
t672 = -t610 + t756;
t671 = (-t720 * t384 + t719 * t386) * t382;
t670 = t718 * t384 + t717 * t386;
t669 = t340 - t343 - t351 + (-Icges(4,2) - Icges(6,2) - Icges(5,3)) * t604 + t741;
t668 = Icges(6,3) * t605 + t734 * t377 - t341 + t743;
t667 = t733 * t382;
t324 = rSges(4,1) * t384 + rSges(4,2) * t386;
t259 = t324 * t553;
t203 = -rSges(4,2) * t602 + t679;
t431 = t203 + t266;
t660 = t382 * t431 - t259;
t175 = t382 * t207;
t656 = t505 - t175 + t704;
t655 = t384 * t669 + t386 * t668;
t597 = -rSges(6,2) * t540 - t705 * t714 + t675;
t654 = (t534 * t382 + t598) * t377 + (t597 - t735) * t378;
t653 = 0.2e1 * t386;
t652 = m(5) / 0.2e1;
t651 = m(6) / 0.2e1;
t647 = t382 / 0.2e1;
t643 = rSges(4,1) * t386;
t530 = t324 * t552;
t436 = -t530 - t541;
t93 = (-t202 - t265) * t382 + t436;
t637 = t378 * t93;
t577 = -t205 - t249;
t575 = t377 * t244 + t378 * t249;
t246 = t321 * t378;
t574 = -t382 * t246 + t377 * t549;
t560 = -t326 - t328;
t558 = -t352 + t372;
t536 = t377 * t141 + t708 * t378;
t241 = t321 * t377;
t535 = -t241 * t553 - t246 * t552 + t376;
t532 = rSges(4,1) * t529 + t215 * rSges(4,2);
t521 = -pkin(2) - t643;
t519 = -t553 / 0.2e1;
t516 = t552 / 0.2e1;
t264 = t378 * rSges(3,1) - rSges(3,2) * t377;
t507 = t141 * t553 + t384 * t546 + t708 * t552;
t504 = t372 + t724;
t218 = rSges(3,1) * t603 - rSges(3,2) * t606;
t497 = t296 - t542;
t489 = -rSges(4,2) * t384 + t643;
t94 = t542 + t660;
t486 = -t377 * t94 - t637;
t456 = t503 + t555;
t98 = (t202 * t377 + t203 * t378) * qJD(3);
t21 = -t388 * t335 + (t331 + t597) * t382 + (t377 * t491 + t382 * t455) * qJD(3) + t433;
t432 = t21 * t377 + t22 * t378 + t683;
t427 = t377 * t521 + t372 + t556;
t425 = t503 + t578;
t393 = (t521 * t637 + (t93 * (-rSges(4,3) - pkin(7)) + t94 * t521) * t377) * t382;
t392 = (((t84 + t682 + t690) * t377 + ((t768 + t772) * t378 + t715 + t754 + t766) * t378) * qJD(3) + t723) * t516 + (((t706 * t378 - t690 + t695) * t378 + (t706 * t377 + t170 + t494 + t696 - t757) * t377) * qJD(3) + t702 + t716) * t519 + (t698 + t700) * t553 / 0.2e1 - (t697 - t699 + t701) * t552 / 0.2e1 + ((t145 + t694) * t377 + (t693 + t744) * t378) * qJD(3) * t647 + (t750 * qJD(3) + t753 * t384 - t752 * t386 + t748 * t519) * t382;
t334 = t378 * t549;
t284 = t489 * qJD(3);
t251 = t321 * t606;
t250 = t487 * t378;
t248 = t324 * t378;
t247 = t488 * t378;
t243 = t324 * t377;
t242 = t488 * t377;
t222 = t554 * t551;
t210 = t264 * t382 + t542;
t166 = -t218 * t382 - t543;
t165 = -t382 * t612 - t544;
t135 = t382 * t679 - t532;
t134 = -rSges(4,1) * t705 + t446;
t72 = (t205 * t378 - t207 * t377) * qJD(3) + t492;
t70 = -t543 - t284 * t552 + (-t135 - t220 + t259) * t382;
t69 = t134 * t382 + (-t284 * t377 - t324 * t603) * qJD(3) + t498;
t8 = ((t139 - t175) * t378 + (t577 * t382 + t137) * t377) * qJD(3) + t507;
t7 = (t547 + t654) * qJD(3) + t507;
t1 = [t392 + m(3) * (t166 * (-t263 - t645) + t165 * (t264 + t381) + (-t218 - t542 + t210) * t209) + (t22 * (t504 - t645) + t64 * (t497 + t680) + t21 * (t381 + t425) + (t64 + t677) * t65 + t678) * m(6) + (t37 * (t558 - t645) + t73 * (t281 + t497) + t36 * (t381 + t456) + (t73 - t423 - t541 + t656) * t74 + t684) * m(5) + (t70 * (t427 - t645) + t93 * (t532 - t542) + t69 * (t381 + t431) + t393 + (-t541 - t436 + t93 + t674) * t94) * m(4); t392 + (t21 * t425 + t22 * t504 + t506 * t636 + t677 * t65 + (t260 + t290 + t296 + t493 + t713) * t64 + t678) * m(6) + (t36 * t456 + t37 * t558 + t533 * t635 + (-t454 + t656) * t74 + (t458 + t281 + t296) * t73 + t684) * m(5) + (t70 * t427 + t69 * t431 + t393 + (t530 + t674) * t94 + (t532 + t660) * t93) * m(4) + (t165 * t264 - t166 * t263 - t209 * t218 - t210 * t612 - (-t209 * t264 - t210 * t263) * t382) * m(3); (t7 * t575 + (t21 * t496 - t7 * t707) * t377 + (t22 * t496 + t7 * t578) * t378 + (-t574 + (qJ(5) * t602 - t250 + t455) * t382 + t676 * t377) * t65 + (-t241 * t382 + t676 * t378 + t251 - t334) * t64 - (t65 * t377 + t64 * t378) * qJD(3) * (-t322 - t326 - t623) + (-t535 - t547 - (t250 * t378 + t727 * t487 - t554 * t624) * qJD(3) + t536 + t654) * t51) * m(6) + (-t73 * (t334 + (t241 - t242) * t382) - t74 * (t247 * t382 + t574) - t72 * t535 - ((t72 * t247 + t560 * t73) * t378 + (t72 * t242 + t560 * t74) * t377) * qJD(3) + t73 * t251 + t8 * t575 + t72 * t536 + (t37 * t561 + t73 * t571 + t8 * t205 + t72 * t139 + (-t72 * t207 + t561 * t74) * t382) * t378 + (t36 * t561 + t74 * t571 - t8 * t207 + t72 * t137 + (-t488 * t73 + t577 * t72) * t382) * t377) * m(5) + (-(t243 * t93 - t248 * t94) * t382 - (t98 * (-t243 * t377 - t248 * t378) + t486 * t489) * qJD(3) + 0.2e1 * t98 * ((t134 + t173) * t378 + (-t203 * t382 + t135) * t377) + t486 * t284 + ((-t382 * t94 - t70) * t378 + (t382 * t93 - t69) * t377) * t324) * m(4) - ((t719 * t384 + t720 * t386) * t382 + ((-t718 * t377 - t669 * t378) * t386 + (t717 * t377 + t668 * t378) * t384) * qJD(3)) * t382 / 0.2e1 + ((t693 * t382 + t699) * t378 + (t694 * t382 + t700) * t377) * t647 + ((t672 * t553 + t667) * t377 + ((t655 * t378 + (t670 + t673) * t377) * qJD(3) + t671) * t378) * t519 + ((-t673 * t552 - t667) * t378 + ((t670 * t377 + (t655 - t672) * t378) * qJD(3) + t671) * t377) * t516 + (t698 * t382 + (t695 * t603 + (t689 * t377 + t696 * t382 - t726) * t377) * t686) * t377 / 0.2e1 - (t697 * t382 + ((t758 * t382 + t726) * t378 + (-t689 * t378 + t382 * t745) * t377) * t686) * t378 / 0.2e1 + (t702 + t722) * t606 / 0.2e1 + (t701 + t721) * t603 / 0.2e1; -m(5) * (t215 * t74 + t216 * t73 + t222 * t72) - m(6) * (t215 * t65 + t216 * t64 + t222 * t51) + ((t552 * t73 + t553 * t74 - t8) * t652 + (t552 * t64 + t553 * t65 - t7) * t651) * t653 + 0.2e1 * ((qJD(3) * t72 + t36 * t377 + t37 * t378 + t603 * t74 - t606 * t73) * t652 + (t603 * t65 - t606 * t64 + t432) * t651) * t384; m(6) * t7 * t384 + (t432 * t651 - m(6) * t554 * t683 / 0.2e1) * t653;];
tauc = t1(:);
