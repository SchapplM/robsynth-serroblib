% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m [6x1]
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:28
% EndTime: 2022-01-20 11:42:48
% DurationCPUTime: 11.67s
% Computational Cost: add. (50836->504), mult. (33297->642), div. (0->0), fcn. (30164->10), ass. (0->333)
t480 = qJ(1) + qJ(2);
t475 = cos(t480);
t479 = qJ(3) + pkin(9);
t472 = cos(t479);
t474 = sin(t480);
t597 = t472 * t474;
t471 = sin(t479);
t601 = t471 * t474;
t509 = rSges(5,1) * t597 - rSges(5,2) * t601 - t475 * rSges(5,3);
t483 = cos(qJ(3));
t477 = t483 * pkin(3);
t468 = t477 + pkin(2);
t650 = -qJ(4) - pkin(7);
t549 = t474 * t468 + t475 * t650;
t275 = -t509 - t549;
t654 = sin(qJ(1)) * pkin(1);
t271 = t275 - t654;
t447 = t474 * t650;
t600 = t471 * t475;
t539 = -rSges(5,2) * t600 + t474 * rSges(5,3);
t647 = rSges(5,1) * t472;
t276 = -t447 + (t468 + t647) * t475 + t539;
t653 = cos(qJ(1)) * pkin(1);
t272 = t276 + t653;
t159 = t271 * t475 + t272 * t474;
t172 = t275 * t475 + t276 * t474;
t705 = m(6) / 0.2e1;
t706 = m(5) / 0.2e1;
t473 = qJ(5) + t479;
t466 = cos(t473);
t604 = t466 * t475;
t465 = sin(t473);
t608 = t465 * t475;
t543 = rSges(6,1) * t604 - rSges(6,2) * t608 + t474 * rSges(6,3);
t651 = pkin(4) * t472;
t422 = t468 + t651;
t478 = pkin(8) - t650;
t755 = t475 * t422 + t478 * t474;
t259 = t543 + t755;
t248 = t653 + t259;
t446 = t478 * t475;
t609 = t465 * t474;
t550 = rSges(6,2) * t609 + t475 * rSges(6,3);
t646 = rSges(6,1) * t466;
t721 = t446 + (-t422 - t646) * t474 + t550;
t737 = t721 - t654;
t727 = -t248 * t474 - t475 * t737;
t750 = -t259 * t474 - t475 * t721;
t640 = (t172 + t159) * t706 + (-t750 - t727) * t705;
t641 = (t159 - t172) * t706 + (t750 - t727) * t705;
t15 = t641 - t640;
t762 = t15 * qJD(1);
t761 = -t248 + t259;
t397 = rSges(6,1) * t465 + rSges(6,2) * t466;
t349 = t397 * t474;
t210 = t737 * t349;
t350 = t397 * t475;
t127 = -t248 * t350 + t210;
t760 = m(6) * t127;
t129 = -t259 * t350 + t349 * t721;
t759 = m(6) * t129;
t481 = sin(qJ(3));
t652 = pkin(3) * t481;
t423 = -pkin(4) * t471 - t652;
t508 = t397 - t423;
t285 = t508 * t475;
t729 = t508 * t474;
t745 = t474 * t729;
t725 = t285 * t475 + t745;
t499 = rSges(5,1) * t471 + rSges(5,2) * t472 + t652;
t730 = t499 * t475;
t731 = t499 * t474;
t736 = t474 * t731 + t475 * t730;
t577 = -m(6) * t725 / 0.2e1 - m(5) * t736 / 0.2e1;
t396 = t475 * t423;
t288 = t396 - t350;
t578 = (-t288 * t475 + t745) * t705 + t736 * t706;
t51 = t578 - t577;
t546 = qJD(1) + qJD(2);
t758 = t546 * t51;
t469 = t474 ^ 2;
t470 = t475 ^ 2;
t547 = t469 + t470;
t756 = t397 * t547;
t746 = t737 * t729;
t449 = Icges(6,4) * t466;
t393 = -Icges(6,2) * t465 + t449;
t394 = Icges(6,1) * t465 + t449;
t754 = t393 + t394;
t450 = Icges(5,4) * t472;
t409 = -Icges(5,2) * t471 + t450;
t410 = Icges(5,1) * t471 + t450;
t753 = t410 + t409;
t476 = Icges(4,4) * t483;
t432 = -Icges(4,2) * t481 + t476;
t433 = Icges(4,1) * t481 + t476;
t752 = t433 + t432;
t751 = -Icges(4,5) * t481 - Icges(5,5) * t471 - Icges(4,6) * t483 - Icges(5,6) * t472;
t106 = t248 * t721 - t259 * t737;
t707 = m(4) / 0.2e1;
t686 = t474 / 0.2e1;
t685 = -t475 / 0.2e1;
t684 = t475 / 0.2e1;
t682 = m(3) * (t653 * (-rSges(3,1) * t474 - rSges(3,2) * t475) + (t475 * rSges(3,1) - t474 * rSges(3,2)) * t654);
t742 = t721 * t729;
t632 = Icges(5,4) * t471;
t408 = Icges(5,2) * t472 + t632;
t411 = Icges(5,1) * t472 - t632;
t633 = Icges(4,4) * t481;
t431 = Icges(4,2) * t483 + t633;
t434 = Icges(4,1) * t483 - t633;
t741 = (-t431 + t434) * t483 - t752 * t481 + (-t408 + t411) * t472 - t753 * t471;
t687 = -t474 / 0.2e1;
t738 = t686 + t687;
t128 = t271 * t731 - t272 * t730;
t130 = t275 * t731 - t276 * t730;
t720 = t446 + t549;
t278 = t422 * t474 - t720;
t719 = t447 + t755;
t279 = -t475 * t468 + t719;
t735 = m(6) * ((t278 + t720) * t475 + (-t279 + (-t422 - t468) * t475 + t719) * t474);
t605 = t466 * t474;
t226 = t474 * (rSges(6,1) * t605 - t550) + t475 * t543;
t467 = t475 * pkin(7);
t566 = -t474 * (pkin(2) * t474 - t467 - t549) + t475 * (-t474 * pkin(7) - t447 + (-pkin(2) + t468) * t475);
t102 = t278 * t474 + t279 * t475 + t226 + t566;
t398 = -rSges(6,2) * t465 + t646;
t517 = Icges(6,5) * t465 + Icges(6,6) * t466;
t343 = t517 * t474;
t344 = t475 * t517;
t631 = Icges(6,4) * t465;
t395 = Icges(6,1) * t466 - t631;
t319 = Icges(6,5) * t474 + t395 * t475;
t392 = Icges(6,2) * t466 + t631;
t560 = -t392 * t475 + t319;
t418 = Icges(6,4) * t609;
t318 = Icges(6,1) * t605 - Icges(6,5) * t475 - t418;
t561 = -Icges(6,2) * t605 + t318 - t418;
t317 = Icges(6,6) * t474 + t393 * t475;
t562 = -t394 * t475 - t317;
t316 = Icges(6,4) * t605 - Icges(6,2) * t609 - Icges(6,6) * t475;
t563 = t394 * t474 + t316;
t715 = (-t560 * t474 + t561 * t475) * t465 + (t562 * t474 + t563 * t475) * t466;
t649 = (-t469 * t344 + (t474 * t343 + t715) * t475) * t686 + (-t470 * t343 + (t475 * t344 + t715) * t474) * t685;
t724 = t474 * t349 + t475 * t350;
t18 = t649 + m(6) * (-t102 * t724 + t398 * t725);
t732 = t18 * qJD(5);
t505 = t724 * t705;
t525 = m(6) * t756;
t165 = t505 + t525 / 0.2e1;
t728 = t546 * t165;
t124 = -t276 * t271 + t272 * t275;
t648 = rSges(4,1) * t483;
t542 = pkin(2) + t648;
t592 = t474 * t481;
t548 = rSges(4,2) * t592 + t475 * rSges(4,3);
t302 = -t542 * t474 + t467 + t548;
t289 = t302 - t654;
t588 = t475 * t481;
t443 = rSges(4,2) * t588;
t303 = -t443 + t542 * t475 + (rSges(4,3) + pkin(7)) * t474;
t290 = t303 + t653;
t131 = -t303 * t289 + t290 * t302;
t723 = t751 * t474;
t722 = t751 * t475;
t435 = rSges(4,1) * t481 + rSges(4,2) * t483;
t545 = ((-t290 + t303) * t475 + (t289 - t302) * t474) * t435 * t707 + (t761 * t285 - t742 + t746) * t705 + (t128 - t130) * t706;
t121 = t248 * t288 + t746;
t125 = t259 * t288 + t742;
t386 = t435 * t474;
t387 = t435 * t475;
t150 = t289 * t386 - t290 * t387;
t166 = t302 * t386 - t303 * t387;
t718 = (t166 + t150) * t707 + (t125 + t121) * t705 + (t130 + t128) * t706;
t335 = Icges(5,5) * t474 + t411 * t475;
t556 = -t408 * t475 + t335;
t333 = Icges(5,6) * t474 + t409 * t475;
t558 = -t410 * t475 - t333;
t717 = -t556 * t601 + t558 * t597;
t441 = Icges(4,4) * t592;
t591 = t474 * t483;
t355 = Icges(4,1) * t591 - Icges(4,5) * t475 - t441;
t553 = -Icges(4,2) * t591 + t355 - t441;
t353 = Icges(4,4) * t591 - Icges(4,2) * t592 - Icges(4,6) * t475;
t555 = t433 * t474 + t353;
t426 = Icges(5,4) * t601;
t334 = Icges(5,1) * t597 - Icges(5,5) * t475 - t426;
t557 = -Icges(5,2) * t597 + t334 - t426;
t332 = Icges(5,4) * t597 - Icges(5,2) * t601 - Icges(5,6) * t475;
t559 = t410 * t474 + t332;
t716 = t557 * t471 + t559 * t472 + t553 * t481 + t555 * t483;
t527 = t754 * t466 / 0.2e1 + (-t392 / 0.2e1 + t395 / 0.2e1) * t465;
t280 = t319 * t605;
t391 = Icges(6,5) * t466 - Icges(6,6) * t465;
t616 = t391 * t475;
t315 = Icges(6,3) * t474 + t616;
t533 = t315 * t475 - t280;
t145 = -t317 * t609 - t533;
t314 = Icges(6,5) * t605 - Icges(6,6) * t609 - Icges(6,3) * t475;
t571 = -t474 * t314 - t318 * t604;
t146 = -t316 * t608 - t571;
t570 = t474 * t315 + t319 * t604;
t147 = -t317 * t608 + t570;
t530 = t317 * t465 - t314;
t621 = t316 * t465;
t540 = ((t145 - t280 + (t315 + t621) * t475 + t571) * t475 + t570 * t474) * t685 + (-t146 * t475 + t147 * t474) * t684 + ((t530 * t474 + t145 + t146 + t533) * t474 + (t147 - t570 + (-t318 * t466 + t621) * t474 + (t530 + t314) * t475) * t475) * t686;
t714 = t753 * t472 / 0.2e1 + t752 * t483 / 0.2e1 + (t434 / 0.2e1 - t431 / 0.2e1) * t481 + (t411 / 0.2e1 - t408 / 0.2e1) * t471;
t712 = 0.4e1 * qJD(1);
t710 = 0.4e1 * qJD(2);
t709 = 2 * qJD(3);
t701 = t226 * t735;
t700 = t102 * t735;
t696 = m(6) * (t129 + t127);
t695 = m(6) * (t210 + (-t721 * t474 + t761 * t475) * t397);
t503 = t727 * t398;
t573 = -t285 * t349 + t350 * t729;
t693 = m(6) * (t503 + t573);
t502 = t750 * t398;
t692 = m(6) * (t502 + t573);
t501 = (-t288 * t474 - t475 * t729) * t397;
t691 = m(6) * (t503 + t501);
t690 = m(6) * (t502 + t501);
t291 = t335 * t597;
t407 = Icges(5,5) * t472 - Icges(5,6) * t471;
t614 = t407 * t475;
t331 = Icges(5,3) * t474 + t614;
t532 = t331 * t475 - t291;
t161 = -t333 * t601 - t532;
t330 = Icges(5,5) * t597 - Icges(5,6) * t601 - Icges(5,3) * t475;
t619 = t332 * t471;
t107 = -(-(-t334 * t472 + t619) * t474 - t330 * t475) * t475 + t161 * t474;
t596 = t472 * t475;
t569 = -t474 * t330 - t334 * t596;
t162 = -t332 * t600 - t569;
t568 = t474 * t331 + t335 * t596;
t163 = -t333 * t600 + t568;
t108 = -t162 * t475 + t163 * t474;
t354 = Icges(4,6) * t474 + t432 * t475;
t356 = Icges(4,5) * t474 + t434 * t475;
t311 = t356 * t591;
t430 = Icges(4,5) * t483 - Icges(4,6) * t481;
t612 = t430 * t475;
t352 = Icges(4,3) * t474 + t612;
t531 = t352 * t475 - t311;
t188 = -t354 * t592 - t531;
t351 = Icges(4,5) * t591 - Icges(4,6) * t592 - Icges(4,3) * t475;
t617 = t353 * t481;
t122 = -(-(-t355 * t483 + t617) * t474 - t351 * t475) * t475 + t188 * t474;
t587 = t475 * t483;
t565 = -t474 * t351 - t355 * t587;
t189 = -t353 * t588 - t565;
t564 = t474 * t352 + t356 * t587;
t190 = -t354 * t588 + t564;
t123 = -t189 * t475 + t190 * t474;
t529 = t333 * t471 - t330;
t38 = (t529 * t475 + t163 - t568) * t475 + (t529 * t474 + t162 + t532) * t474;
t39 = (t161 - t291 + (t331 + t619) * t475 + t569) * t475 + t568 * t474;
t528 = t354 * t481 - t351;
t45 = (t528 * t475 + t190 - t564) * t475 + (t528 * t474 + t189 + t531) * t474;
t46 = (t188 - t311 + (t352 + t617) * t475 + t565) * t475 + t564 * t474;
t2 = t700 + (-t46 / 0.2e1 - t39 / 0.2e1 + t123 / 0.2e1 + t108 / 0.2e1) * t475 + (t122 / 0.2e1 + t107 / 0.2e1 + t45 / 0.2e1 + t38 / 0.2e1) * t474 + t540;
t683 = t2 * qJD(3) - qJD(4) * t51;
t678 = m(4) * t131;
t676 = m(4) * t150;
t675 = m(4) * t166;
t672 = m(5) * t124;
t670 = m(5) * t128;
t669 = m(5) * t130;
t668 = m(5) * t159;
t667 = m(5) * t172;
t663 = m(6) * t106;
t661 = m(6) * t121;
t660 = m(6) * t125;
t659 = m(6) * t727;
t658 = m(6) * t750;
t639 = -t165 * qJD(4) + qJD(5) * t540;
t164 = t505 - t525 / 0.2e1;
t53 = t577 + t578;
t638 = t53 * qJD(3) + t164 * qJD(5);
t637 = t51 * qJD(3) + t165 * qJD(5);
t576 = (-t285 - t288) * t729;
t554 = -t433 * t475 - t354;
t552 = -t431 * t475 + t356;
t541 = rSges(5,2) * t471 - t477 - t647;
t524 = t701 / 0.2e1 + t540;
t523 = t696 / 0.2e1 + t527;
t507 = -t398 - t477 - t651;
t500 = (-t386 * t475 + t387 * t474) * t435;
t493 = (-t392 + t395) * t466 - t754 * t465;
t498 = -t540 + (t391 * t474 + t562 * t465 + t560 * t466 + t475 * t493) * t686 + (-t563 * t465 + t561 * t466 + t474 * t493 - t616) * t685;
t496 = -t527 + t738 * (t316 * t466 + t318 * t465);
t494 = -t552 * t481 + t554 * t483;
t490 = -t701 / 0.2e1 + t498;
t489 = (-t349 * t475 + t350 * t474) * t397;
t488 = t527 + t714;
t487 = t488 + t718;
t486 = t496 - t714 + (t332 * t472 + t334 * t471 + t353 * t483 + t355 * t481) * t738;
t485 = t53 * qJD(4) + (t498 - t700 + (t107 + t122 + t38 + t45) * t687 + (t558 * t471 + t556 * t472 + t554 * t481 + t552 * t483 + t741 * t475 + (t407 + t430) * t474) * t686 + (-t559 * t471 + t557 * t472 + t741 * t474 - t555 * t481 + t553 * t483 + t108 + t123 - t612 - t614) * t685 + (t39 + t46) * t684) * qJD(3);
t437 = -rSges(4,2) * t481 + t648;
t340 = t541 * t475;
t338 = t541 * t474;
t286 = t507 * t475;
t284 = t507 * t474;
t154 = t164 * qJD(4);
t136 = t475 * t396 + t423 * t469 - t724;
t105 = -t658 + t667;
t99 = -t659 + t668;
t88 = t527 + t759;
t87 = t527 + t760;
t85 = t690 / 0.2e1;
t82 = t691 / 0.2e1;
t81 = t692 / 0.2e1;
t79 = t693 / 0.2e1;
t69 = t695 / 0.2e1;
t40 = t663 + t672 + t678 + t682;
t35 = t488 + t660 + t669 + t675;
t34 = t488 + t661 + t670 + t676;
t22 = m(6) * (-t226 * t724 + t398 * t756) + t649;
t21 = t22 * qJD(5);
t20 = -t695 / 0.2e1 + t523;
t19 = t69 + t523;
t16 = t640 + t641;
t14 = t69 - t696 / 0.2e1 + t496;
t11 = t487 - t545;
t10 = t487 + t545;
t9 = t81 - t690 / 0.2e1 + t524;
t8 = t85 - t692 / 0.2e1 + t524;
t7 = t79 - t691 / 0.2e1 + t524;
t6 = t82 - t693 / 0.2e1 + t524;
t5 = t486 + t545 - t718;
t4 = t81 + t85 + t490;
t3 = t79 + t82 + t490;
t1 = [t40 * qJD(2) + t34 * qJD(3) + t99 * qJD(4) + t87 * qJD(5), t40 * qJD(1) + t10 * qJD(3) + t16 * qJD(4) + t19 * qJD(5) + 0.2e1 * (t106 * t705 + t124 * t706 + t131 * t707 + t682 / 0.2e1) * qJD(2), t34 * qJD(1) + t10 * qJD(2) + t3 * qJD(5) + ((t271 * t340 + t272 * t338) * t706 + ((-t289 * t475 - t290 * t474) * t437 + t500) * t707 + (t248 * t284 + t286 * t737 + t576) * t705) * t709 + t485, qJD(1) * t99 + qJD(2) * t16 + t638, t87 * qJD(1) + t19 * qJD(2) + t3 * qJD(3) + t154 + ((t503 + t489) * m(6) + t498) * qJD(5); t11 * qJD(3) - t15 * qJD(4) + t20 * qJD(5) + (-t663 / 0.4e1 - t672 / 0.4e1 - t678 / 0.4e1 - t682 / 0.4e1) * t712, qJD(3) * t35 + qJD(4) * t105 + qJD(5) * t88, t11 * qJD(1) + t35 * qJD(2) + t4 * qJD(5) + ((t275 * t340 + t276 * t338) * t706 + ((-t302 * t475 - t303 * t474) * t437 + t500) * t707 + (t259 * t284 + t286 * t721 + t576) * t705) * t709 + t485, qJD(2) * t105 + t638 - t762, t20 * qJD(1) + t88 * qJD(2) + t4 * qJD(3) + t154 + ((t502 + t489) * m(6) + t498) * qJD(5); t486 * qJD(1) + t5 * qJD(2) + t7 * qJD(5) + (-t661 / 0.4e1 - t670 / 0.4e1 - t676 / 0.4e1) * t712 + t683, t5 * qJD(1) + t486 * qJD(2) + t9 * qJD(5) + (-t660 / 0.4e1 - t669 / 0.4e1 - t675 / 0.4e1) * t710 + t683, (m(6) * (t102 * t136 - t284 * t729 - t285 * t286) + m(5) * (-t731 * t338 - t730 * t340 - (t474 * t509 + t475 * (rSges(5,1) * t596 + t539) + t566) * t499 * t547) + m(4) * ((t474 * (rSges(4,1) * t591 - t548) + t475 * (rSges(4,1) * t587 + t474 * rSges(4,3) - t443)) * (-t386 * t474 - t387 * t475) + t547 * t437 * t435) + t649 + ((t716 * t475 + (t494 - t723) * t474 + t717) * t475 + t722 * t469) * t686 + ((t494 * t474 + (t716 - t722) * t475 + t717) * t474 + t723 * t470) * t685) * qJD(3) + t732 + t546 * t2, -t758, t7 * qJD(1) + t9 * qJD(2) + t18 * qJD(3) + t732; t15 * qJD(2) + (t659 / 0.4e1 - t668 / 0.4e1) * t712 + t637, t762 + (t658 / 0.4e1 - t667 / 0.4e1) * t710 + t637, ((-t284 * t475 + t286 * t474) * t705 + (-t338 * t475 + t340 * t474) * t706) * t709 + t758, 0, t728; (t496 - t760) * qJD(1) + t14 * qJD(2) + t6 * qJD(3) + t639, t14 * qJD(1) + (t496 - t759) * qJD(2) + t8 * qJD(3) + t639, t6 * qJD(1) + t8 * qJD(2) + ((t136 * t226 + (-t284 * t474 - t286 * t475) * t397) * m(6) + t649) * qJD(3) + t21, -t728, qJD(3) * t22 + t540 * t546 + t21;];
Cq = t1;
