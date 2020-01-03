% Calculate vector of inverse dynamics joint torques for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:18
% EndTime: 2019-12-31 18:14:48
% DurationCPUTime: 25.10s
% Computational Cost: add. (8204->666), mult. (13037->775), div. (0->0), fcn. (9971->6), ass. (0->351)
t732 = Icges(6,4) + Icges(5,5);
t731 = Icges(5,6) - Icges(6,6);
t730 = Icges(4,3) + Icges(5,3);
t356 = sin(qJ(1));
t358 = cos(qJ(1));
t351 = qJ(3) + pkin(7);
t325 = cos(t351);
t324 = sin(t351);
t574 = Icges(5,4) * t324;
t419 = Icges(5,2) * t325 + t574;
t729 = t731 * t356 - t358 * t419;
t573 = Icges(5,4) * t325;
t422 = Icges(5,1) * t324 + t573;
t728 = -t732 * t356 + t358 * t422;
t355 = sin(qJ(3));
t357 = cos(qJ(3));
t725 = Icges(4,5) * t355 + Icges(5,5) * t324 + Icges(4,6) * t357 + Icges(5,6) * t325;
t555 = t325 * t358;
t557 = t324 * t358;
t722 = Icges(6,5) * t557 - Icges(6,3) * t555 + t729;
t418 = Icges(6,4) * t324 - Icges(6,6) * t325;
t147 = Icges(6,2) * t358 + t356 * t418;
t721 = t725 * t356 + t730 * t358;
t727 = t147 + t721;
t282 = Icges(6,5) * t555;
t726 = Icges(6,1) * t557 - t282 + t728;
t569 = Icges(6,5) * t324;
t233 = Icges(6,1) * t325 + t569;
t235 = Icges(5,1) * t325 - t574;
t716 = t233 + t235;
t568 = Icges(6,2) * t356;
t148 = Icges(6,4) * t557 - Icges(6,6) * t555 - t568;
t659 = -t730 * t356 + t725 * t358;
t713 = -t148 - t659;
t576 = Icges(4,4) * t355;
t420 = Icges(4,2) * t357 + t576;
t171 = -Icges(4,6) * t356 + t358 * t420;
t575 = Icges(4,4) * t357;
t423 = Icges(4,1) * t355 + t575;
t173 = -Icges(4,5) * t356 + t358 * t423;
t687 = t171 * t357 + t173 * t355 + t726 * t324 - t722 * t325;
t724 = t687 * t358;
t558 = t324 * t356;
t281 = Icges(6,5) * t558;
t556 = t325 * t356;
t566 = Icges(6,6) * t358;
t143 = -Icges(6,3) * t556 + t281 + t566;
t149 = Icges(5,6) * t358 + t356 * t419;
t723 = -t143 + t149;
t315 = Icges(6,5) * t325;
t421 = Icges(6,1) * t324 - t315;
t151 = Icges(6,4) * t358 + t356 * t421;
t283 = Icges(5,4) * t556;
t570 = Icges(5,5) * t358;
t153 = Icges(5,1) * t558 + t283 + t570;
t720 = -t151 - t153;
t718 = t716 * t358;
t231 = -Icges(5,2) * t324 + t573;
t655 = Icges(6,3) * t324 + t315;
t717 = -t231 + t655;
t170 = Icges(4,6) * t358 + t356 * t420;
t553 = t356 * t357;
t308 = Icges(4,4) * t553;
t554 = t355 * t356;
t571 = Icges(4,5) * t358;
t172 = Icges(4,1) * t554 + t308 + t571;
t693 = -t149 * t325 - t153 * t324 - t170 * t357 - t172 * t355;
t414 = -Icges(6,3) * t325 + t569;
t715 = t414 - t419;
t714 = -t421 - t422;
t635 = -t151 * t324 + t693;
t712 = Icges(4,5) * t357 - Icges(4,6) * t355 - t731 * t324 + t732 * t325;
t548 = t358 * t147 + t151 * t558;
t50 = -t143 * t556 + t548;
t670 = t149 * t556 + t153 * t558 + t170 * t553 + t172 * t554 + t358 * t721;
t644 = t50 + t670;
t266 = -Icges(4,2) * t355 + t575;
t268 = Icges(4,1) * t357 - t576;
t711 = t231 * t325 + t235 * t324 + t266 * t357 + t268 * t355;
t563 = t143 * t325;
t710 = t563 + t635;
t636 = -t143 * t555 - t727 * t356;
t701 = t233 * t324 - t325 * t655 + t711;
t181 = t231 * t358;
t497 = qJD(3) * t358;
t709 = qJD(3) * t181 - t655 * t497 + (t356 * t414 - t149 + t566) * qJD(1);
t174 = t655 * t356;
t498 = qJD(3) * t356;
t708 = qJD(3) * t174 - t231 * t498 + (t358 * t414 + t729) * qJD(1);
t707 = -t718 * qJD(3) + (t356 * t422 + t151 + t570) * qJD(1);
t184 = t235 * t356;
t706 = qJD(3) * t184 + t233 * t498 + (t358 * t421 + t728) * qJD(1);
t705 = t715 * qJD(3);
t704 = t714 * qJD(3);
t703 = (Icges(6,1) * t556 + t184 + t281 - t723) * t358 + (-t718 - t722) * t356;
t643 = -t171 * t553 - t173 * t554 + t358 * t713 + t556 * t722 - t558 * t726;
t642 = -t151 * t557 + t358 * t693 - t636;
t641 = t356 * t713 + t724;
t639 = t171 * t355 - t173 * t357 - t324 * t722 - t325 * t726;
t640 = t170 * t355 - t172 * t357 + t324 * t723 + t325 * t720;
t676 = t712 * t356;
t702 = -t266 * t355 + t268 * t357 + t324 * t717 + t325 * t716;
t700 = -t418 - t725;
t673 = t712 * t358;
t699 = (Icges(5,2) * t558 + t174 - t283 + t720) * t358 + (-Icges(6,3) * t557 + t181 - t282 + t726) * t356;
t698 = t659 * qJD(1);
t697 = -t714 - t717;
t696 = -t715 - t716;
t638 = t701 * t356 + t673;
t637 = t233 * t557 + t711 * t358 - t555 * t655 - t676;
t695 = t727 * qJD(1);
t694 = t358 ^ 2;
t610 = rSges(6,1) + pkin(4);
t578 = rSges(6,3) + qJ(5);
t314 = t325 * qJ(5);
t316 = t325 * rSges(6,3);
t436 = rSges(6,1) * t324 - t316;
t524 = pkin(4) * t324 - t314 + t436;
t606 = pkin(3) * t355;
t692 = -t524 - t606;
t521 = t266 + t423;
t522 = -t420 + t268;
t691 = (t697 * t324 + t696 * t325 + t355 * t521 - t357 * t522) * qJD(1);
t690 = t701 * qJD(1) + t700 * qJD(3);
t689 = t698 + t676 * qJD(3) + (t358 * t418 - t568 - t710) * qJD(1);
t688 = -t687 * qJD(1) - t673 * qJD(3) + t695;
t686 = t641 * t356 + t642 * t358;
t685 = t643 * t356 + t644 * t358;
t109 = qJD(1) * t171 + t266 * t498;
t220 = t268 * t356;
t111 = qJD(1) * t173 + qJD(3) * t220;
t684 = t640 * qJD(3) - t109 * t357 - t111 * t355 - t706 * t324 + t708 * t325 + t695;
t245 = t420 * qJD(3);
t246 = t423 * qJD(3);
t683 = -qJD(1) * t712 + t702 * qJD(3) - t245 * t357 - t246 * t355 + t704 * t324 + t705 * t325;
t219 = t266 * t358;
t108 = qJD(1) * t170 - qJD(3) * t219;
t221 = t268 * t358;
t110 = -qJD(3) * t221 + (t356 * t423 + t571) * qJD(1);
t682 = qJD(1) * t148 + t639 * qJD(3) + t108 * t357 + t110 * t355 + t707 * t324 - t709 * t325 + t698;
t368 = t356 * (t173 + t219) - t358 * (-Icges(4,2) * t554 + t172 + t308);
t369 = t356 * (t171 - t221) - t358 * (t170 - t220);
t681 = -t703 * t324 + t699 * t325 - t369 * t355 + t368 * t357;
t680 = rSges(4,2) * t355;
t240 = pkin(4) * t325 + qJ(5) * t324;
t241 = rSges(6,1) * t325 + rSges(6,3) * t324;
t523 = t240 + t241;
t588 = rSges(5,2) * t325;
t437 = rSges(5,1) * t324 + t588;
t389 = t437 + t606;
t501 = qJD(1) * t358;
t679 = t324 * t501 + t325 * t498;
t486 = qJD(3) ^ 2 * t606;
t678 = qJD(4) * qJD(1) + t486;
t677 = t638 * qJD(1);
t675 = t637 * qJD(1);
t674 = t700 * qJD(1);
t669 = rSges(4,2) * t357;
t347 = t358 * rSges(5,3);
t157 = rSges(5,1) * t558 + rSges(5,2) * t556 + t347;
t331 = qJD(2) * t358;
t465 = t357 * t497;
t496 = qJD(4) * t356;
t452 = -pkin(3) * t465 + t496;
t441 = -t331 + t452;
t312 = pkin(3) * t554;
t354 = -qJ(4) - pkin(6);
t595 = pkin(6) + t354;
t203 = -t358 * t595 + t312;
t274 = t358 * pkin(1) + t356 * qJ(2);
t602 = pkin(6) * t358;
t455 = t274 + t602;
t446 = t203 + t455;
t589 = rSges(5,2) * t324;
t592 = rSges(5,1) * t325;
t242 = -t589 + t592;
t471 = t242 * t497;
t49 = -t471 + (t157 + t446) * qJD(1) + t441;
t668 = t358 * t49;
t667 = t325 * t610;
t352 = t356 ^ 2;
t511 = t352 + t694;
t605 = pkin(3) * t357;
t661 = t511 * t605;
t202 = t356 * t595 + t358 * t606;
t333 = t358 * qJ(2);
t270 = pkin(1) * t356 - t333;
t253 = qJD(1) * t270;
t660 = qJD(1) * t202 - t253;
t348 = t358 * rSges(4,3);
t194 = rSges(4,1) * t554 + rSges(4,2) * t553 + t348;
t658 = t194 + t455;
t518 = -t356 * rSges(6,2) - rSges(6,3) * t555;
t349 = t358 * rSges(6,2);
t657 = t558 * t610 + t349;
t473 = t355 * t501;
t502 = qJD(1) * t356;
t656 = pkin(3) * t473 + t354 * t502;
t275 = -rSges(3,2) * t358 + t356 * rSges(3,3);
t470 = t324 * t498;
t654 = -t470 * t578 - t610 * t679;
t653 = -qJ(5) * t555 + t518;
t273 = rSges(4,1) * t357 - t680;
t223 = t273 * t358;
t438 = rSges(4,1) * t355 + t669;
t112 = -qJD(3) * t223 + (t356 * t438 + t348) * qJD(1);
t212 = qJD(1) * t274 - t331;
t247 = t438 * qJD(3);
t490 = qJD(1) * qJD(3);
t250 = qJDD(3) * t356 + t358 * t490;
t491 = qJD(1) * qJD(2);
t515 = qJDD(2) * t356 + t358 * t491;
t601 = pkin(6) * qJD(1) ^ 2;
t396 = -t358 * t601 + t515;
t195 = -t356 * rSges(4,3) + t358 * t438;
t603 = pkin(6) * t356;
t456 = -t270 - t603;
t448 = t195 + t456;
t33 = -t247 * t498 + t250 * t273 + (-t112 - t212) * qJD(1) + t448 * qJDD(1) + t396;
t466 = t357 * t498;
t476 = t501 * t669 + (t466 + t473) * rSges(4,1);
t499 = qJD(3) * t355;
t113 = (-rSges(4,2) * t499 - rSges(4,3) * qJD(1)) * t356 + t476;
t251 = qJDD(3) * t358 - t356 * t490;
t320 = qJ(2) * t501;
t330 = qJD(2) * t356;
t514 = t320 + t330;
t480 = qJD(1) * (-pkin(1) * t502 + t514) + qJDD(1) * t274 + t356 * t491;
t564 = qJDD(1) * pkin(6);
t388 = -t356 * t601 + t358 * t564 + t480;
t34 = qJD(1) * t113 + qJDD(1) * t194 - t251 * t273 + (qJD(3) * t247 - qJDD(2)) * t358 + t388;
t652 = t33 * t356 - t34 * t358;
t651 = qJD(1) * t523;
t650 = qJD(3) * t685 + t677;
t649 = qJD(3) * t686 - t675;
t648 = t356 * t690 - t358 * t683;
t647 = t356 * t683 + t358 * t690;
t646 = t710 * qJD(3) - t109 * t355 + t111 * t357 + t708 * t324 + t706 * t325;
t645 = t687 * qJD(3) - t108 * t355 + t110 * t357 + t709 * t324 + t707 * t325;
t634 = t689 * t694 + (t682 * t356 + (-t684 + t688) * t358) * t356;
t633 = t684 * t694 + (t688 * t356 + (-t682 + t689) * t358) * t356;
t311 = qJD(5) * t324;
t538 = -qJD(3) * t524 + t311;
t623 = -qJD(3) * (t311 + t538) + qJDD(5) * t325;
t517 = pkin(3) * t466 + qJD(4) * t358;
t450 = t517 + t656;
t117 = pkin(6) * t502 + t450;
t95 = -t471 + (t356 * t437 + t347) * qJD(1);
t477 = rSges(5,1) * t679 + t501 * t588;
t483 = qJD(3) * t589;
t97 = (-rSges(5,3) * qJD(1) - t483) * t356 + t477;
t622 = t358 * t95 + (-t117 - t97) * t356;
t617 = -m(5) - m(6);
t616 = -pkin(1) - pkin(6);
t614 = t250 / 0.2e1;
t613 = t251 / 0.2e1;
t612 = t356 / 0.2e1;
t609 = rSges(3,2) - pkin(1);
t608 = -rSges(6,2) - pkin(1);
t607 = -rSges(5,3) - pkin(1);
t600 = g(2) * t358;
t387 = qJDD(4) * t358 + t250 * t605 + t396;
t535 = t557 * t610 + t653;
t481 = t202 + t535;
t493 = qJD(5) * t358;
t307 = t354 * t501;
t118 = -t307 + (t312 - t602) * qJD(1) + t452;
t549 = -t118 - t212;
t192 = t241 * t358;
t594 = (pkin(4) * t502 - qJ(5) * t497) * t324 + (-qJ(5) * t502 + (-pkin(4) * qJD(3) + qJD(5)) * t358) * t325 - qJD(3) * t192 + (t356 * t436 + t349) * qJD(1);
t2 = t523 * t250 + (-t270 + t481) * qJDD(1) + (-t325 * t493 + t549 - t594) * qJD(1) + (-t564 - t623 - t678) * t356 + t387;
t599 = t2 * t356;
t361 = qJD(1) * t117 + qJDD(1) * t203 + qJDD(4) * t356 + t358 * t678 + t388;
t449 = -t523 - t605;
t494 = qJD(5) * t356;
t464 = t325 * t494;
t537 = -t556 * t578 + t657;
t593 = -(-qJ(5) * t501 - t494) * t325 - t518 * qJD(1) + t654;
t3 = (-qJDD(2) + t623) * t358 + t449 * t251 + t537 * qJDD(1) + (-t464 - t593) * qJD(1) + t361;
t598 = t3 * t358;
t587 = rSges(3,3) * t358;
t457 = -t242 - t605;
t214 = t437 * qJD(3);
t500 = qJD(3) * t214;
t13 = (-qJDD(2) + t500) * t358 + qJD(1) * t97 + qJDD(1) * t157 + t457 * t251 + t361;
t586 = t13 * t358;
t343 = t356 * rSges(5,3);
t159 = t358 * t437 - t343;
t447 = t202 + t456;
t395 = t159 + t447;
t475 = t330 + t517;
t451 = t242 * t498 + t475;
t48 = qJD(1) * t395 + t451;
t582 = t358 * t48;
t236 = t273 * t498;
t78 = qJD(1) * t448 + t236 + t330;
t581 = t358 * t78;
t550 = t118 * t497 - t251 * t202;
t536 = -t157 - t203;
t530 = t556 * t610 + t578 * t558;
t529 = -t240 * t358 - t192;
t271 = rSges(3,2) * t356 + t587;
t520 = -t270 + t271;
t201 = t274 + t275;
t513 = rSges(3,2) * t502 + rSges(3,3) * t501;
t512 = t330 - t253;
t495 = qJD(5) * t325;
t487 = -rSges(4,3) + t616;
t313 = pkin(3) * t553;
t485 = -t117 + t593;
t484 = rSges(5,1) * t555;
t482 = -t203 - t537;
t468 = t355 * t498;
t461 = -t498 / 0.2e1;
t460 = t498 / 0.2e1;
t459 = -t497 / 0.2e1;
t458 = t497 / 0.2e1;
t454 = t325 * t578;
t453 = t148 - t563;
t189 = rSges(5,1) * t556 - rSges(5,2) * t558;
t445 = g(1) * t356 - t600;
t439 = t333 + (-pkin(1) + t354) * t356;
t276 = rSges(2,1) * t358 - rSges(2,2) * t356;
t272 = rSges(2,1) * t356 + rSges(2,2) * t358;
t79 = qJD(1) * t658 - t273 * t497 - t331;
t424 = t356 * t78 - t358 * t79;
t413 = t112 * t358 - t113 * t356;
t403 = -t194 * t356 - t195 * t358;
t394 = -t354 * t358 + t274 + t312;
t393 = t307 - t441;
t391 = -t159 * t358 + t356 * t536;
t390 = -t464 + t475;
t362 = t498 * t523 + t390;
t290 = rSges(5,2) * t557;
t222 = t273 * t356;
t193 = t290 - t484;
t165 = t358 * t202;
t164 = t202 * t497;
t132 = qJD(1) * t201 - t331;
t131 = qJD(1) * t520 + t330;
t115 = t358 * t118;
t98 = t403 * qJD(3);
t63 = qJD(1) * t513 + qJDD(1) * t275 - qJDD(2) * t358 + t480;
t62 = t520 * qJDD(1) + (-qJD(1) * t275 - t212) * qJD(1) + t515;
t47 = qJD(3) * t391 - t164;
t39 = (-qJD(3) * t523 + t495) * t358 + (t446 + t537) * qJD(1) + t441;
t38 = (t447 + t535) * qJD(1) + t362;
t32 = -t164 + t311 + (t356 * t482 - t358 * t535) * qJD(3);
t12 = t242 * t250 + (-t486 - t500) * t356 + t395 * qJDD(1) + (-t95 - t496 + t549) * qJD(1) + t387;
t1 = qJDD(5) * t324 - t535 * t251 + t482 * t250 + (t485 * t356 + t594 * t358 + t495) * qJD(3) + t550;
t4 = [-m(2) * (-g(1) * t272 + g(2) * t276) - t637 * t250 / 0.2e1 + t639 * t614 + (((-t636 - t642 + t643) * t356 + (t548 - t724 + (t453 + t635 + t659) * t356 + t641 + t670) * t358) * qJD(3) + t677) * t461 + (-t701 * qJD(3) + t245 * t355 - t246 * t357 - t705 * t324 + t704 * t325) * qJD(1) + (-(-t38 + (t535 - t603) * qJD(1) + t362 + t660) * t39 + t38 * t393 + t39 * (t320 + t390 - t654 + t656) + (-t495 + (t324 * t578 + t667) * qJD(3)) * t358 * t38 + ((t38 * t608 - t39 * t454) * t358 + (t38 * (-qJ(2) + t692) + t39 * t608) * t356) * qJD(1) + (t3 - g(2)) * (-t356 * t454 + t394 + t657) + (t2 - g(1)) * ((t324 * t610 + t606) * t358 + t439 + t653)) * m(6) + (-(-t48 + (t159 - t603) * qJD(1) + t451 + t660) * t49 + t48 * (qJD(3) * t484 - t358 * t483 + t393) + t49 * (-rSges(5,2) * t470 + t450 + t477 + t514) + (t607 * t582 + (t48 * (-qJ(2) - t389) + t49 * t607) * t356) * qJD(1) + (t13 - g(2)) * (t394 + t157) + (t12 - g(1)) * (t358 * t389 - t343 + t439)) * m(5) + (-(t236 - t78 + (t195 - t603) * qJD(1) + t512) * t79 + t78 * (rSges(4,1) * t465 - t497 * t680 + t331) + t79 * (-rSges(4,2) * t468 + t476 + t514) + (t487 * t581 + (t78 * (-qJ(2) - t438) + t79 * t487) * t356) * qJD(1) + (t34 - g(2)) * t658 + (t33 - g(1)) * (t356 * t616 + t195 + t333)) * m(4) + (-(qJD(1) * t271 - t131 + t512) * t132 + t131 * t331 + t132 * (t513 + t514) + (t131 * t609 * t358 + (t131 * (-rSges(3,3) - qJ(2)) - t132 * pkin(1)) * t356) * qJD(1) + (-g(2) + t63) * t201 + (-g(1) + t62) * (t356 * t609 + t333 + t587)) * m(3) + (t638 - t640) * t613 + (t646 + t647) * t458 + (m(2) * (t272 ^ 2 + t276 ^ 2) + Icges(3,1) + Icges(2,3) + t702) * qJDD(1) + (((t356 * t453 - t50 + t548) * t356 + t659 * t352 + ((-t635 - t713) * t358 + t636 + t643) * t358) * qJD(3) + t649 + t675) * t459 + (t645 + t648 + t650) * t460; t617 * t445 + 0.2e1 * (t599 / 0.2e1 - t598 / 0.2e1) * m(6) + 0.2e1 * (t12 * t612 - t586 / 0.2e1) * m(5) + (-t358 * t63 + 0.2e1 * t62 * t612 - t445) * m(3) + (-t445 + t652) * m(4); t686 * t614 + t685 * t613 + (qJD(1) * t648 + qJD(3) * t633 - qJDD(1) * t637 + t250 * t641 + t251 * t642) * t612 + (t647 * qJD(1) + t634 * qJD(3) + t638 * qJDD(1) + t643 * t250 + t644 * t251) * t358 / 0.2e1 - ((t699 * t324 + t703 * t325 + t355 * t368 + t357 * t369) * qJD(3) + (t696 * t324 - t697 * t325 - t355 * t522 - t357 * t521) * qJD(1)) * qJD(1) / 0.2e1 + (t646 * t358 + t645 * t356 + (t640 * t356 + t639 * t358) * qJD(1)) * qJD(1) / 0.2e1 + (t639 * t356 - t640 * t358) * qJDD(1) / 0.2e1 - t650 * t502 / 0.2e1 + t649 * t501 / 0.2e1 + ((-t498 * t673 + t674) * t356 + ((t356 * t676 + t681) * qJD(3) + t691) * t358) * t461 + ((-t356 * t642 + t358 * t641) * qJD(1) + t633) * t460 + ((t497 * t676 + t674) * t358 + ((-t358 * t673 - t681) * qJD(3) - t691) * t356) * t459 + ((-t356 * t644 + t643 * t358) * qJD(1) + t634) * t458 + (t2 * t313 - t1 * t165 + (-t1 * t535 + t3 * t449 + t38 * t651 - t39 * t538) * t358 + (t2 * t523 + t38 * (-pkin(3) * t499 + t538) + t1 * t482 + t39 * t651) * t356 - (-t38 * t529 + t39 * t530) * qJD(1) - (t38 * t692 * t356 + t39 * t524 * t358) * qJD(3) - g(1) * (t313 + t530) - g(3) * (t314 + t316 - t606) - (-t605 - t667) * t600 + (g(3) * t610 - t38 * t494 + t39 * t493 + t578 * t600) * t324 + (t115 + (qJD(1) * t482 + t594) * t358 + (qJD(1) * t481 + t485) * t356 - t495 - (-t356 * t530 + t358 * t529 - t661) * qJD(3)) * t32) * m(6) + (t12 * (t242 * t356 + t313) + t48 * (-pkin(3) * t468 - t214 * t356) + t457 * t586 + t214 * t668 + (qJD(3) * t622 - t159 * t251 + t536 * t250 + t550) * (-t165 + t391) + t47 * (t115 + t622) - (t437 * t668 + t47 * (t193 * t358 - t661) + (-t47 * t189 - t389 * t48) * t356) * qJD(3) - g(1) * (t189 + t313) - g(2) * (t290 + (-t592 - t605) * t358) + g(3) * t389 + (t47 * (t536 * t358 + (t159 + t202) * t356) + (t356 * t49 + t582) * t242 + t48 * t193 - t49 * t189) * qJD(1)) * m(5) + (-(t222 * t79 + t223 * t78) * qJD(1) - (t98 * (-t222 * t356 - t223 * t358) - t424 * t438) * qJD(3) - g(1) * t222 + g(2) * t223 + g(3) * t438 + (qJD(3) * t413 - t194 * t250 - t195 * t251) * t403 + t98 * ((-t194 * t358 + t195 * t356) * qJD(1) + t413) - t424 * t247 + ((t356 * t79 + t581) * qJD(1) + t652) * t273) * m(4); t617 * (g(1) * t358 + g(2) * t356) + m(5) * (t12 * t358 + t13 * t356) + m(6) * (t2 * t358 + t3 * t356); ((t1 - g(3)) * t324 + (t445 + t598 - t599 + (-t511 + 0.1e1) * t32 * qJD(3)) * t325) * m(6);];
tau = t4;
