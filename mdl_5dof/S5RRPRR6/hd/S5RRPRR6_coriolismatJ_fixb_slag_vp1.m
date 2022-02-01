% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:16:51
% EndTime: 2022-01-20 11:17:12
% DurationCPUTime: 13.80s
% Computational Cost: add. (80557->516), mult. (74468->690), div. (0->0), fcn. (80559->10), ass. (0->329)
t486 = qJ(1) + qJ(2);
t483 = cos(t486);
t487 = sin(pkin(9));
t542 = t487 * (-pkin(8) - pkin(7));
t467 = t483 * t542;
t491 = cos(qJ(4));
t479 = pkin(4) * t491 + pkin(3);
t481 = sin(t486);
t488 = cos(pkin(9));
t485 = qJ(4) + qJ(5);
t480 = sin(t485);
t482 = cos(t485);
t590 = t483 * t488;
t434 = -t480 * t590 + t481 * t482;
t435 = t480 * t481 + t482 * t590;
t512 = t435 * rSges(6,1) + t434 * rSges(6,2);
t520 = -rSges(6,3) * t487 - pkin(2);
t489 = sin(qJ(4));
t622 = pkin(4) * t489;
t281 = -t467 + (qJ(3) + t622) * t481 + (t479 * t488 - t520) * t483 + t512;
t624 = cos(qJ(1)) * pkin(1);
t278 = t281 + t624;
t625 = sin(qJ(1)) * pkin(1);
t476 = t483 * qJ(3);
t589 = t483 * t489;
t595 = t481 * t488;
t535 = pkin(4) * t589 - t479 * t595 + t481 * t542;
t432 = t480 * t595 + t482 * t483;
t433 = -t483 * t480 + t482 * t595;
t695 = -t433 * rSges(6,1) + t432 * rSges(6,2);
t708 = t481 * t520 + t476 + t535 + t695;
t727 = t708 - t625;
t167 = t278 * t481 + t483 * t727;
t168 = t281 * t481 + t483 * t708;
t581 = t488 * t491;
t594 = t481 * t489;
t455 = t483 * t581 + t594;
t582 = t488 * t489;
t706 = -t481 * t491 + t483 * t582;
t514 = t455 * rSges(5,1) - rSges(5,2) * t706;
t623 = pkin(3) * t488;
t688 = (rSges(5,3) + pkin(7)) * t487 + pkin(2) + t623;
t314 = t481 * qJ(3) + t688 * t483 + t514;
t302 = t314 + t624;
t452 = t481 * t582 + t483 * t491;
t453 = t481 * t581 - t589;
t694 = -t453 * rSges(5,1) + t452 * rSges(5,2);
t709 = -t688 * t481 + t476 + t694;
t726 = t709 - t625;
t205 = t302 * t481 + t483 * t726;
t221 = t314 * t481 + t483 * t709;
t521 = rSges(4,1) * t488 + pkin(2);
t596 = t481 * t487;
t392 = rSges(4,2) * t596 + t483 * rSges(4,3) - t481 * t521 + t476;
t390 = t392 - t625;
t591 = t483 * t487;
t393 = -rSges(4,2) * t591 + t521 * t483 + (rSges(4,3) + qJ(3)) * t481;
t391 = t393 + t624;
t289 = t390 * t483 + t391 * t481;
t299 = t392 * t483 + t393 * t481;
t679 = m(6) / 0.2e1;
t680 = m(5) / 0.2e1;
t719 = m(4) / 0.2e1;
t537 = (t168 + t167) * t679 + (t221 + t205) * t680 + (t299 + t289) * t719;
t538 = (-t168 + t167) * t679 + (-t221 + t205) * t680 + (t289 - t299) * t719;
t18 = t538 - t537;
t744 = t18 * qJD(1);
t375 = -rSges(6,1) * t432 - rSges(6,2) * t433;
t376 = t434 * rSges(6,1) - t435 * rSges(6,2);
t141 = t281 * t376 - t375 * t708;
t743 = m(6) * t141;
t484 = t487 ^ 2;
t138 = t278 * t376 - t375 * t727;
t742 = m(6) * t138;
t349 = rSges(6,3) * t596 - t695;
t621 = pkin(7) * t487;
t364 = (t621 + t623) * t481 + t535;
t620 = -pkin(3) + t479;
t414 = -t488 * pkin(8) + t620 * t487;
t398 = t414 * t596;
t427 = -rSges(6,3) * t488 + (rSges(6,1) * t482 - rSges(6,2) * t480) * t487;
t399 = t427 * t596;
t201 = t398 + t399 + (t349 - t364) * t488;
t737 = t488 * t349 + t399;
t202 = -t364 * t488 + t398 + t737;
t741 = m(6) * (-t201 + t202);
t351 = rSges(6,3) * t591 + t512;
t545 = pkin(4) * t594;
t559 = t351 + t545 - t467 + (t620 * t488 - t621) * t483;
t203 = t559 * t488 + (t414 + t427) * t591;
t740 = t203 * t727;
t114 = t278 * t708 - t281 * t727;
t136 = t302 * t709 - t314 * t726;
t366 = rSges(5,3) * t596 - t694;
t448 = -rSges(5,3) * t488 + (rSges(5,1) * t491 - rSges(5,2) * t489) * t487;
t736 = t366 * t488 + t448 * t596;
t735 = -t487 / 0.2e1;
t734 = t203 * t708;
t610 = Icges(6,4) * t480;
t426 = -Icges(6,5) * t488 + (Icges(6,1) * t482 - t610) * t487;
t551 = t426 + (-Icges(6,2) * t482 - t610) * t487;
t729 = t480 * t551;
t613 = Icges(5,4) * t489;
t447 = -Icges(5,5) * t488 + (Icges(5,1) * t491 - t613) * t487;
t548 = t447 + (-Icges(5,2) * t491 - t613) * t487;
t728 = t489 * t548;
t416 = Icges(6,4) * t433;
t343 = -Icges(6,2) * t432 + Icges(6,6) * t596 + t416;
t415 = Icges(6,4) * t432;
t347 = -Icges(6,1) * t433 - Icges(6,5) * t596 + t415;
t724 = t434 * t343 - t347 * t435;
t437 = Icges(5,4) * t453;
t358 = -Icges(5,2) * t452 + Icges(5,6) * t596 + t437;
t436 = Icges(5,4) * t452;
t362 = -Icges(5,1) * t453 - Icges(5,5) * t596 + t436;
t722 = -t358 * t452 - t362 * t453;
t721 = -t358 * t706 - t362 * t455;
t659 = m(3) * (t624 * (-rSges(3,1) * t481 - rSges(3,2) * t483) + (t483 * rSges(3,1) - t481 * rSges(3,2)) * t625);
t655 = m(4) * (-t393 * t390 + t391 * t392);
t340 = Icges(6,5) * t433 - Icges(6,6) * t432 + Icges(6,3) * t596;
t714 = t340 * t591;
t163 = t714 + t724;
t697 = t163 - t714;
t718 = (t697 - t724) * t591;
t445 = -Icges(5,3) * t488 + (Icges(5,5) * t491 - Icges(5,6) * t489) * t487;
t612 = Icges(5,4) * t491;
t446 = -Icges(5,6) * t488 + (-Icges(5,2) * t489 + t612) * t487;
t716 = (t445 * t596 - t446 * t452 + t447 * t453) * t488;
t355 = Icges(5,5) * t453 - Icges(5,6) * t452 + Icges(5,3) * t596;
t712 = t355 * t591;
t194 = t203 * t481;
t122 = t202 * t483 - t194;
t368 = rSges(5,3) * t591 + t514;
t305 = t368 * t488 + t448 * t591;
t580 = t122 * t679 + (-t305 * t481 + t483 * t736) * t680;
t615 = (-t201 * t483 + t122 + t194) * t679;
t705 = t580 - t615;
t578 = -t202 * t281 - t734;
t509 = (t201 * t281 + t578 + t734) * t679;
t579 = t201 * t278 + t740;
t616 = (t578 + t579) * t679 + ((t302 - t314) * t736 + (-t709 + t726) * t305) * t680;
t704 = t616 - t509;
t459 = (-Icges(5,5) * t489 - Icges(5,6) * t491) * t487;
t583 = t488 * t459;
t703 = -t583 / 0.2e1 + t728 * t735;
t449 = (-Icges(6,5) * t480 - Icges(6,6) * t482) * t487;
t584 = t488 * t449;
t702 = -t584 / 0.2e1 + t729 * t735;
t439 = t452 * pkin(4);
t334 = t349 * t591;
t156 = t334 + (-t364 * t483 - t481 * t559) * t487;
t352 = t375 * t591;
t264 = -t376 * t596 + t352;
t456 = (-rSges(6,1) * t480 - rSges(6,2) * t482) * t487;
t315 = t488 * t375 + t456 * t596;
t316 = -t376 * t488 - t456 * t591;
t369 = -Icges(6,5) * t432 - Icges(6,6) * t433;
t561 = -Icges(6,2) * t433 - t347 - t415;
t563 = -Icges(6,1) * t432 - t343 - t416;
t142 = -t369 * t488 + (-t480 * t561 + t482 * t563) * t487;
t370 = Icges(6,5) * t434 - Icges(6,6) * t435;
t417 = Icges(6,4) * t434;
t348 = Icges(6,1) * t435 + Icges(6,5) * t591 + t417;
t560 = -Icges(6,2) * t435 + t348 + t417;
t611 = Icges(6,4) * t435;
t345 = Icges(6,2) * t434 + Icges(6,6) * t591 + t611;
t562 = Icges(6,1) * t434 - t345 - t611;
t143 = -t370 * t488 + (-t480 * t560 + t482 * t562) * t487;
t609 = Icges(6,4) * t482;
t425 = -Icges(6,6) * t488 + (-Icges(6,2) * t480 + t609) * t487;
t451 = (-Icges(6,1) * t480 - t609) * t487;
t552 = -t425 + t451;
t175 = -t432 * t551 + t433 * t552 + t449 * t596;
t176 = t434 * t551 + t435 * t552 + t449 * t591;
t526 = t591 / 0.2e1;
t530 = t596 / 0.2e1;
t606 = (-t584 + (t482 * t552 - t729) * t487) * t488;
t662 = -t488 / 0.2e1;
t546 = ((t370 * t591 + t434 * t560 + t435 * t562) * t591 + (t369 * t591 + t434 * t561 + t435 * t563) * t596 - t176 * t488) * t526 + ((t370 * t596 - t432 * t560 + t433 * t562) * t591 + (t369 * t596 - t432 * t561 + t433 * t563) * t596 - t175 * t488) * t530 + (-t606 + (t142 * t481 + t143 * t483) * t487) * t662;
t14 = t546 + m(6) * (t156 * t264 + t202 * t315 - t203 * t316);
t701 = t14 * qJD(5);
t182 = t712 + t721;
t696 = t182 - t712;
t547 = qJD(1) + qJD(2);
t357 = Icges(5,5) * t455 - Icges(5,6) * t706 + Icges(5,3) * t591;
t614 = Icges(5,4) * t455;
t360 = -Icges(5,2) * t706 + Icges(5,6) * t591 + t614;
t438 = Icges(5,4) * t706;
t363 = Icges(5,1) * t455 + Icges(5,5) * t591 - t438;
t691 = t357 * t596 - t452 * t360 + t453 * t363;
t327 = -t375 + t439;
t388 = -rSges(5,1) * t452 - rSges(5,2) * t453;
t389 = -rSges(5,1) * t706 - rSges(5,2) * t455;
t554 = t706 * pkin(4) - t376;
t567 = (t327 * t481 + t483 * t554) * t679 + (-t388 * t481 - t389 * t483) * t680;
t129 = -t278 * t554 + t327 * t727;
t134 = -t281 * t554 + t327 * t708;
t154 = t302 * t389 - t388 * t726;
t160 = t314 * t389 - t388 * t709;
t690 = (t160 + t154) * t680 + (t134 + t129) * t679;
t689 = t203 * t741;
t686 = 0.4e1 * qJD(1);
t685 = 0.2e1 * qJD(2);
t684 = 0.4e1 * qJD(2);
t683 = 2 * qJD(4);
t288 = t351 * t488 + t427 * t591;
t676 = t288 * t741;
t573 = t316 * t278 + t315 * t727;
t576 = -t202 * t375 - t203 * t376;
t672 = m(6) * (t573 + t576);
t572 = t316 * t281 + t315 * t708;
t671 = m(6) * (t572 + t576);
t670 = m(6) * ((t278 - t281) * t737 + (-t708 + t727) * t288);
t568 = t288 * t554 + t327 * t737;
t667 = m(6) * (t568 + t573);
t665 = m(6) * (t568 + t572);
t663 = m(6) * (t141 + t138);
t180 = t355 * t596 + t722;
t607 = t180 * t481;
t100 = -t716 + (t691 * t483 + t607) * t487;
t183 = t357 * t591 - t360 * t706 + t455 * t363;
t603 = (t445 * t591 - t446 * t706 + t447 * t455) * t488;
t101 = -t603 + (t182 * t481 + t183 * t483) * t487;
t41 = t716 + ((-t691 + t696 - t721) * t483 - t607) * t487;
t42 = -t603 + ((t180 + t183 - t722) * t483 + t696 * t481) * t487;
t164 = (Icges(6,5) * t435 + Icges(6,6) * t434 + Icges(6,3) * t591) * t591 + t434 * t345 + t435 * t348;
t605 = ((-Icges(6,3) * t488 + (Icges(6,5) * t482 - Icges(6,6) * t480) * t487) * t591 + t425 * t434 + t426 * t435) * t488;
t38 = -t605 + ((t340 * t596 + t164) * t483 + t697 * t481) * t487;
t531 = -t596 / 0.2e1;
t94 = -t605 + (t163 * t481 + t164 * t483) * t487;
t517 = t38 * t530 + t718 * t526 + t94 * t531;
t2 = ((t41 / 0.2e1 + t100 / 0.2e1) * t483 + (-t101 / 0.2e1 + t42 / 0.2e1) * t481) * t487 + t517 + t689;
t29 = t580 + t615 - t567;
t660 = t29 * qJD(3) + t2 * qJD(4);
t653 = m(4) * t289;
t652 = m(4) * t299;
t647 = m(5) * t136;
t645 = m(5) * t154;
t644 = m(5) * t160;
t643 = m(5) * t205;
t641 = m(5) * t221;
t637 = m(6) * t114;
t633 = m(6) * t129;
t632 = m(6) * t134;
t631 = m(6) * t167;
t630 = m(6) * t168;
t629 = m(6) * (-t288 * t481 + t483 * t737);
t628 = m(6) * (t315 * t481 - t316 * t483);
t626 = m(6) * (-t375 * t481 - t376 * t483);
t179 = t629 / 0.2e1;
t74 = t179 - t626 / 0.2e1;
t619 = t74 * qJD(3) + qJD(5) * t517;
t27 = t567 + t705;
t266 = t626 / 0.2e1;
t73 = t266 + t179;
t618 = t27 * qJD(4) + t73 * qJD(5);
t28 = t567 - t705;
t72 = t266 - t629 / 0.2e1;
t617 = t28 * qJD(4) + t72 * qJD(5);
t461 = (-Icges(5,1) * t489 - t612) * t487;
t549 = -t446 + t461;
t604 = (-t583 + (t491 * t549 - t728) * t487) * t488;
t592 = t482 * t487;
t587 = t487 * t491;
t577 = t202 * t327 + t203 * t554;
t566 = -t288 * t376 - t375 * t737;
t565 = -t305 * t389 - t388 * t736;
t558 = -Icges(5,1) * t452 - t358 - t437;
t557 = -Icges(5,1) * t706 - t360 - t614;
t556 = -Icges(5,2) * t453 - t362 - t436;
t555 = -Icges(5,2) * t455 + t363 - t438;
t529 = -t592 / 0.2e1;
t528 = t592 / 0.2e1;
t527 = -t591 / 0.2e1;
t523 = -t587 / 0.2e1;
t522 = t587 / 0.2e1;
t511 = t425 * t529 + t451 * t528 + t702;
t510 = t676 / 0.2e1 + t517;
t508 = t663 / 0.2e1 + t511;
t502 = t425 * t528 + t451 * t529 - t702;
t499 = t446 * t523 + t461 * t522 + t511 + t703;
t498 = t38 * t531 - t606 + t718 * t527 + (t143 + t176) * t526 + (t94 + t142 + t175) * t530;
t497 = t499 + t690;
t496 = -t676 / 0.2e1 + t498;
t494 = t446 * t522 + t461 * t523 + t502 - t703;
t379 = -Icges(5,5) * t452 - Icges(5,6) * t453;
t152 = -t379 * t488 + (-t489 * t556 + t491 * t558) * t487;
t380 = -Icges(5,5) * t706 - Icges(5,6) * t455;
t153 = -t380 * t488 + (-t489 * t555 + t491 * t557) * t487;
t216 = -t452 * t548 + t453 * t549 + t459 * t596;
t217 = t455 * t549 + t459 * t591 - t548 * t706;
t493 = t27 * qJD(3) + (t42 * t531 + t498 - t604 - t689 + (t101 + t152 + t216) * t530 + (t100 + t41) * t527 + (t153 + t217) * t526) * qJD(4);
t462 = (-rSges(5,1) * t489 - rSges(5,2) * t491) * t487;
t330 = -t389 * t488 - t462 * t591;
t329 = t388 * t488 + t462 * t596;
t263 = t554 * t488 + (-t456 * t487 + t484 * t622) * t483;
t262 = -t439 * t488 - t484 * t545 + t315;
t249 = -t351 * t596 + t334;
t220 = t352 + (-t439 * t483 + t481 * t554) * t487;
t210 = qJD(5) * t628;
t107 = t511 + t743;
t106 = t511 + t742;
t105 = t630 + t641 + t652;
t102 = t631 + t643 + t653;
t78 = t665 / 0.2e1;
t75 = t667 / 0.2e1;
t70 = t73 * qJD(3);
t61 = t670 / 0.2e1;
t59 = t671 / 0.2e1;
t57 = t672 / 0.2e1;
t56 = t499 + t632 + t644;
t55 = t499 + t633 + t645;
t49 = t637 + t647 + t655 + t659;
t22 = -t670 / 0.2e1 + t508;
t21 = t61 + t508;
t19 = t537 + t538;
t17 = t61 - t663 / 0.2e1 + t502;
t16 = m(6) * (t249 * t264 - t288 * t316 + t315 * t737) + t546;
t15 = t16 * qJD(5);
t13 = t497 - t704;
t12 = t497 + t704;
t9 = t494 + t509 + t616 - t690;
t8 = t78 - t671 / 0.2e1 + t510;
t7 = t59 - t665 / 0.2e1 + t510;
t6 = t75 - t672 / 0.2e1 + t510;
t5 = t57 - t667 / 0.2e1 + t510;
t4 = t496 + t59 + t78;
t3 = t496 + t57 + t75;
t1 = [t49 * qJD(2) + t102 * qJD(3) + t55 * qJD(4) + t106 * qJD(5), t49 * qJD(1) + t19 * qJD(3) + t12 * qJD(4) + t21 * qJD(5) + (t114 * t679 + t136 * t680 + t655 / 0.2e1 + t659 / 0.2e1) * t685, qJD(1) * t102 + qJD(2) * t19 + t618, t55 * qJD(1) + t12 * qJD(2) + t3 * qJD(5) + ((t262 * t727 + t263 * t278 + t577) * t679 + (t302 * t330 + t329 * t726 + t565) * t680) * t683 + t493, t106 * qJD(1) + t21 * qJD(2) + t70 + t3 * qJD(4) + (m(6) * (t566 + t573) + t498) * qJD(5); -t18 * qJD(3) + t13 * qJD(4) + t22 * qJD(5) + (-t637 / 0.4e1 - t647 / 0.4e1 - t655 / 0.4e1 - t659 / 0.4e1) * t686, qJD(3) * t105 + qJD(4) * t56 + qJD(5) * t107, qJD(2) * t105 + t618 - t744, t13 * qJD(1) + t56 * qJD(2) + t4 * qJD(5) + ((t262 * t708 + t263 * t281 + t577) * t679 + (t314 * t330 + t329 * t709 + t565) * t680) * t683 + t493, t22 * qJD(1) + t107 * qJD(2) + t70 + t4 * qJD(4) + (m(6) * (t566 + t572) + t498) * qJD(5); t18 * qJD(2) + (-t631 / 0.4e1 - t643 / 0.4e1 - t653 / 0.4e1) * t686 + t617, t744 + (-t630 / 0.4e1 - t641 / 0.4e1 - t652 / 0.4e1) * t684 + t617, 0, ((t329 * t481 - t330 * t483) * t680 + (t262 * t481 - t263 * t483) * t679) * t683 + t210 + t547 * t28, qJD(4) * t628 + t547 * t72 + t210; t9 * qJD(2) + t5 * qJD(5) + (-t633 / 0.4e1 - t645 / 0.4e1) * t686 + t660 + (t494 + 0.2e1 * (-t202 * t278 + t579 - t740) * t679) * qJD(1), t9 * qJD(1) + t494 * qJD(2) + t7 * qJD(5) + (-t632 / 0.4e1 - t644 / 0.4e1) * t684 + t509 * t685 + t660, t547 * t29, (m(5) * (t736 * t329 - t305 * t330 + (t366 * t483 - t368 * t481) * t484 * (t388 * t483 - t389 * t481)) + ((t380 * t591 + t455 * t557 - t555 * t706) * t591 + (t379 * t591 + t455 * t558 - t556 * t706) * t596 - t217 * t488) * t526 + ((t380 * t596 - t452 * t555 + t453 * t557) * t591 + (t379 * t596 - t452 * t556 + t453 * t558) * t596 - t216 * t488) * t530 + (-t604 + (t152 * t481 + t153 * t483) * t487) * t662 + m(6) * (t156 * t220 + t202 * t262 - t203 * t263) + t546) * qJD(4) + t701 + t547 * t2, t5 * qJD(1) + t7 * qJD(2) + t14 * qJD(4) + t701; (t502 - t742) * qJD(1) + t17 * qJD(2) + t6 * qJD(4) + t619, t17 * qJD(1) + (t502 - t743) * qJD(2) + t8 * qJD(4) + t619, t547 * t74, t6 * qJD(1) + t8 * qJD(2) + ((t220 * t249 + t262 * t737 - t263 * t288) * m(6) + t546) * qJD(4) + t15, qJD(4) * t16 + t517 * t547 + t15;];
Cq = t1;
