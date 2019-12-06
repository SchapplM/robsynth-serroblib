% Calculate vector of inverse dynamics joint torques for
% S5PRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:35
% EndTime: 2019-12-05 15:38:40
% DurationCPUTime: 52.45s
% Computational Cost: add. (16691->780), mult. (26649->1096), div. (0->0), fcn. (25508->8), ass. (0->394)
t712 = Icges(6,4) + Icges(5,5);
t711 = Icges(5,6) - Icges(6,6);
t686 = Icges(5,2) + Icges(6,3);
t740 = Icges(6,2) + Icges(5,3);
t364 = pkin(8) + qJ(4);
t354 = sin(t364);
t355 = cos(t364);
t370 = sin(qJ(2));
t739 = (-t712 * t354 - t711 * t355) * t370;
t366 = sin(pkin(7));
t368 = cos(pkin(7));
t371 = cos(qJ(2));
t537 = qJD(4) * t371;
t516 = t355 * t537;
t541 = qJD(2) * t370;
t520 = t366 * t541;
t176 = -t366 * t516 + (qJD(4) * t368 + t520) * t354;
t575 = t366 * t371;
t258 = t354 * t575 + t355 * t368;
t177 = -qJD(4) * t258 - t355 * t520;
t540 = qJD(2) * t371;
t519 = t366 * t540;
t738 = t711 * t176 + t712 * t177 + t740 * t519;
t518 = t368 * t541;
t583 = t354 * t366;
t178 = -qJD(4) * t583 + t354 * t518 - t368 * t516;
t572 = t368 * t371;
t260 = t354 * t572 - t366 * t355;
t179 = -qJD(4) * t260 - t355 * t518;
t517 = t368 * t540;
t737 = t711 * t178 + t712 * t179 + t740 * t517;
t736 = Icges(5,1) + Icges(6,1);
t735 = Icges(5,4) - Icges(6,5);
t259 = -t368 * t354 + t355 * t575;
t576 = t366 * t370;
t107 = Icges(5,5) * t259 - Icges(5,6) * t258 + Icges(5,3) * t576;
t109 = Icges(6,4) * t259 + Icges(6,2) * t576 + Icges(6,6) * t258;
t708 = t107 + t109;
t261 = t355 * t572 + t583;
t573 = t368 * t370;
t108 = Icges(5,5) * t261 - Icges(5,6) * t260 + Icges(5,3) * t573;
t110 = Icges(6,4) * t261 + Icges(6,2) * t573 + Icges(6,6) * t260;
t707 = t110 + t108;
t467 = Icges(5,5) * t355 - Icges(5,6) * t354;
t224 = Icges(5,3) * t370 + t371 * t467;
t470 = Icges(6,4) * t355 + Icges(6,6) * t354;
t226 = Icges(6,2) * t370 + t371 * t470;
t734 = -t739 * qJD(4) + (-t224 - t226) * qJD(2);
t590 = Icges(6,5) * t355;
t466 = Icges(6,3) * t354 + t590;
t222 = Icges(6,6) * t370 + t371 * t466;
t596 = Icges(5,4) * t355;
t471 = -Icges(5,2) * t354 + t596;
t228 = Icges(5,6) * t370 + t371 * t471;
t680 = -t222 + t228;
t223 = -Icges(5,3) * t371 + t370 * t467;
t225 = -Icges(6,2) * t371 + t370 * t470;
t723 = -t225 - t223;
t591 = Icges(6,5) * t354;
t474 = Icges(6,1) * t355 + t591;
t230 = Icges(6,4) * t370 + t371 * t474;
t597 = Icges(5,4) * t354;
t475 = Icges(5,1) * t355 - t597;
t232 = Icges(5,5) * t370 + t371 * t475;
t678 = -t230 - t232;
t733 = (-t686 * t355 + t591 - t597) * t370;
t732 = -t686 * t176 - t735 * t177 - t711 * t519;
t731 = -t686 * t178 - t735 * t179 - t711 * t517;
t730 = t735 * t176 + t736 * t177 + t712 * t519;
t729 = t735 * t178 + t736 * t179 + t712 * t517;
t240 = Icges(6,5) * t259;
t105 = Icges(6,6) * t576 + Icges(6,3) * t258 + t240;
t599 = Icges(5,4) * t259;
t111 = -Icges(5,2) * t258 + Icges(5,6) * t576 + t599;
t710 = t105 - t111;
t241 = Icges(6,5) * t261;
t106 = Icges(6,6) * t573 + Icges(6,3) * t260 + t241;
t598 = Icges(5,4) * t261;
t112 = -Icges(5,2) * t260 + Icges(5,6) * t573 + t598;
t709 = t106 - t112;
t593 = Icges(6,5) * t258;
t113 = Icges(6,1) * t259 + Icges(6,4) * t576 + t593;
t242 = Icges(5,4) * t258;
t115 = Icges(5,1) * t259 + Icges(5,5) * t576 - t242;
t728 = t113 + t115;
t592 = Icges(6,5) * t260;
t114 = Icges(6,1) * t261 + Icges(6,4) * t573 + t592;
t243 = Icges(5,4) * t260;
t116 = Icges(5,1) * t261 + Icges(5,5) * t573 - t243;
t727 = t114 + t116;
t726 = t680 * qJD(2) + t733 * qJD(4);
t275 = (-Icges(5,1) * t354 - t596) * t370;
t538 = qJD(4) * t370;
t725 = -(-Icges(6,1) * t354 + t590) * t538 - qJD(4) * t275 + t678 * qJD(2);
t580 = t355 * t370;
t326 = Icges(6,5) * t580;
t582 = t354 * t370;
t587 = Icges(6,6) * t371;
t221 = Icges(6,3) * t582 + t326 - t587;
t227 = -Icges(5,6) * t371 + t370 * t471;
t724 = t221 - t227;
t229 = -Icges(6,4) * t371 + t370 * t474;
t231 = -Icges(5,5) * t371 + t370 * t475;
t679 = t229 + t231;
t722 = t734 * t370 + t723 * t540;
t721 = t737 * t370 + t707 * t540;
t720 = t738 * t370 + t708 * t540;
t696 = -t710 * t176 + t728 * t177 + t732 * t258 + t730 * t259 + t720 * t366;
t695 = -t709 * t176 + t727 * t177 + t731 * t258 + t729 * t259 + t721 * t366;
t694 = -t710 * t178 + t728 * t179 + t732 * t260 + t730 * t261 + t720 * t368;
t693 = -t709 * t178 + t727 * t179 + t731 * t260 + t729 * t261 + t721 * t368;
t719 = t724 * t176 - t679 * t177 + t726 * t258 + t725 * t259 + t722 * t366;
t718 = t724 * t178 - t679 * t179 + t726 * t260 + t725 * t261 + t722 * t368;
t691 = t710 * t258 + t728 * t259 + t708 * t576;
t717 = t709 * t258 + t727 * t259 + t707 * t576;
t716 = t710 * t260 + t728 * t261 + t708 * t573;
t690 = t709 * t260 + t727 * t261 + t707 * t573;
t465 = t105 * t354 + t113 * t355;
t54 = -t109 * t371 + t370 * t465;
t461 = -t111 * t354 + t115 * t355;
t56 = -t107 * t371 + t370 * t461;
t689 = t54 + t56;
t464 = t106 * t354 + t114 * t355;
t55 = -t110 * t371 + t370 * t464;
t460 = -t112 * t354 + t116 * t355;
t57 = -t108 * t371 + t370 * t460;
t688 = t55 + t57;
t715 = t724 * t258 + t679 * t259 - t723 * t576;
t714 = t724 * t260 + t679 * t261 - t723 * t573;
t455 = t221 * t354 + t229 * t355;
t584 = t225 * t371;
t80 = t370 * t455 - t584;
t454 = -t227 * t354 + t231 * t355;
t585 = t223 * t371;
t81 = t370 * t454 - t585;
t653 = t80 + t81;
t469 = Icges(3,5) * t371 - Icges(3,6) * t370;
t713 = (-Icges(3,3) * t368 + t366 * t469) * t368;
t362 = t366 ^ 2;
t363 = t368 ^ 2;
t644 = t362 + t363;
t543 = qJD(2) * t366;
t308 = t368 * t538 + t543;
t542 = qJD(2) * t368;
t309 = t366 * t538 - t542;
t704 = (t308 * t717 + t309 * t691 - t537 * t715) * t366 + (t308 * t690 + t309 * t716 - t537 * t714) * t368;
t323 = rSges(3,1) * t371 - rSges(3,2) * t370;
t451 = t644 * t323;
t356 = qJD(3) * t370;
t339 = t366 * t356;
t320 = pkin(2) * t370 - qJ(3) * t371;
t430 = qJD(2) * t320;
t215 = -t366 * t430 + t339;
t341 = t368 * t356;
t216 = -t368 * t430 + t341;
t343 = qJ(3) * t575;
t344 = qJ(3) * t572;
t701 = t366 * t215 + t368 * t216 - (-pkin(2) * t576 + t343) * t543 - (-pkin(2) * t573 + t344) * t542 - t356;
t663 = rSges(6,3) + qJ(5);
t664 = rSges(6,1) + pkin(4);
t552 = t663 * t580 - t582 * t664;
t700 = 0.2e1 * qJD(2);
t699 = 2 * qJDD(2);
t531 = qJD(2) * qJD(4);
t421 = qJDD(4) * t370 + t371 * t531;
t530 = qJDD(2) * t366;
t238 = t368 * t421 + t530;
t529 = qJDD(2) * t368;
t239 = t366 * t421 - t529;
t315 = -qJDD(4) * t371 + t370 * t531;
t698 = t717 * t238 + t691 * t239 + t695 * t308 + t309 * t696 + t715 * t315 + t719 * t537;
t697 = t238 * t690 + t239 * t716 + t308 * t693 + t309 * t694 + t315 * t714 + t537 * t718;
t656 = ((t461 + t465) * qJD(2) - t738) * t371 + (t730 * t355 + t732 * t354 + (-t354 * t728 + t355 * t710) * qJD(4) + t708 * qJD(2)) * t370;
t655 = ((t460 + t464) * qJD(2) - t737) * t371 + (t729 * t355 + t731 * t354 + (-t354 * t727 + t355 * t709) * qJD(4) + t707 * qJD(2)) * t370;
t692 = t308 * t688 + t309 * t689 - t537 * t653;
t396 = -t370 * t466 + t587;
t186 = t396 * t366;
t192 = t227 * t366;
t685 = t186 + t192;
t187 = t396 * t368;
t193 = t227 * t368;
t684 = t187 + t193;
t194 = t229 * t366;
t196 = t231 * t366;
t683 = -t194 - t196;
t195 = t229 * t368;
t197 = t231 * t368;
t682 = -t195 - t197;
t439 = t224 - t454;
t440 = -t226 + t455;
t633 = (t225 * t368 + t464) * t308 + (t225 * t366 + t465) * t309;
t634 = -(-t223 * t368 - t460) * t308 - (-t223 * t366 - t461) * t309;
t677 = (-t633 - t634 + (-t439 + t440) * t537) * t370;
t676 = t739 * t537 + (t258 * t712 + t259 * t711) * t309 + (t260 * t712 + t261 * t711) * t308;
t675 = t308 * t707 + t309 * t708;
t367 = cos(pkin(8));
t352 = pkin(3) * t367 + pkin(2);
t616 = pkin(2) - t352;
t511 = t616 * t370;
t369 = -pkin(6) - qJ(3);
t570 = qJ(3) + t369;
t658 = -t371 * t570 + t511;
t403 = qJD(2) * t658;
t172 = t366 * t403;
t173 = t368 * t403;
t571 = t371 * t369;
t432 = t511 - t571;
t674 = t366 * t172 + t368 * t173 - (t366 * t432 - t343) * t543 - (t368 * t432 - t344) * t542 + t701;
t673 = -t584 - t585;
t672 = t714 * t370;
t671 = t715 * t370;
t670 = t717 * t368;
t669 = t716 * t366;
t365 = sin(pkin(8));
t611 = rSges(4,2) * t365;
t613 = rSges(4,1) * t367;
t489 = -t611 + t613;
t636 = rSges(4,3) * t371 - t370 * t489;
t666 = t644 * qJD(2) * t636;
t560 = -t320 + t658;
t499 = qJD(2) * t560;
t438 = t366 * t499 + t339;
t536 = qJD(5) * t258;
t487 = rSges(6,1) * t355 + rSges(6,3) * t354;
t660 = (pkin(4) * t355 + qJ(5) * t354 + t487) * t370;
t557 = -rSges(6,2) * t371 + t660;
t567 = rSges(6,2) * t573 + t260 * t663 + t261 * t664;
t44 = -t308 * t557 - t537 * t567 + t438 + t536;
t437 = t368 * t499 + t341;
t535 = qJD(5) * t260;
t568 = rSges(6,2) * t576 + t258 * t663 + t259 * t664;
t45 = t309 * t557 + t537 * t568 + t437 + t535;
t665 = t308 * t44 - t309 * t45 - g(3);
t321 = rSges(3,1) * t370 + rSges(3,2) * t371;
t662 = t321 * t644;
t532 = qJD(2) * qJD(3);
t659 = qJDD(3) * t370 + t371 * t532;
t654 = ((-t454 - t455) * qJD(2) - t734) * t371 + (t725 * t355 + t726 * t354 + (t354 * t679 - t355 * t724) * qJD(4) + t723 * qJD(2)) * t370;
t645 = t371 * pkin(2) + t370 * qJ(3);
t643 = (t679 + t733) * t537 + (t259 * t686 + t242 - t593 - t728) * t309 + (t261 * t686 + t243 - t592 - t727) * t308;
t642 = (Icges(6,1) * t582 - t275 - t326 - t724) * t537 + (-t258 * t736 + t240 - t599 + t710) * t309 + (-t260 * t736 + t241 - t598 + t709) * t308;
t641 = t676 * t370;
t479 = t56 * t366 + t57 * t368;
t480 = t54 * t366 + t55 * t368;
t640 = t479 + t480;
t639 = t368 * t690 + t669;
t638 = t366 * t691 + t670;
t637 = g(1) * t368 + g(2) * t366;
t534 = qJD(5) * t354;
t319 = t370 * t534;
t412 = -t370 * t570 - t371 * t616;
t577 = t365 * t368;
t170 = -pkin(3) * t577 + t366 * t412;
t578 = t365 * t366;
t171 = pkin(3) * t578 + t368 * t412;
t305 = t645 * t366;
t307 = t645 * t368;
t539 = qJD(3) * t371;
t449 = t305 * t543 + t307 * t542 + qJD(1) - t539;
t414 = t170 * t543 + t171 * t542 + t449;
t27 = t308 * t568 - t309 * t567 + t319 + t414;
t404 = -qJDD(3) * t371 + t215 * t543 + t216 * t542 + t305 * t530 + t307 * t529 + t370 * t532 + qJDD(1);
t383 = t170 * t530 + t171 * t529 + t172 * t543 + t173 * t542 + t404;
t416 = t354 * t540 + t355 * t538;
t614 = rSges(6,2) * t517 - t178 * t663 + t179 * t664 + t535;
t615 = rSges(6,2) * t519 - t176 * t663 + t177 * t664 + t536;
t5 = qJD(5) * t416 + qJDD(5) * t582 + t238 * t568 - t239 * t567 + t308 * t615 - t309 * t614 + t383;
t635 = t27 * t614 + t5 * t567;
t630 = t238 / 0.2e1;
t629 = t239 / 0.2e1;
t628 = -t644 * t541 / 0.2e1;
t627 = -t308 / 0.2e1;
t626 = t308 / 0.2e1;
t625 = -t309 / 0.2e1;
t624 = t309 / 0.2e1;
t623 = t315 / 0.2e1;
t620 = -t371 / 0.2e1;
t612 = rSges(5,1) * t355;
t360 = t370 * rSges(6,2);
t359 = t370 * rSges(4,3);
t358 = t370 * rSges(5,3);
t581 = t354 * t371;
t579 = t355 * t371;
t574 = t367 * t371;
t569 = t319 + t416 * qJ(5) + (-t354 * t538 + t355 * t540) * pkin(4) + (-rSges(6,1) * t354 + rSges(6,3) * t355) * t538 + (t371 * t487 + t360) * qJD(2);
t566 = -t258 * t664 + t259 * t663;
t565 = t260 * t664 - t261 * t663;
t347 = rSges(6,2) * t575;
t564 = -t366 * t660 + t347;
t350 = rSges(6,2) * t572;
t563 = t368 * t660 - t350;
t292 = qJD(2) * t645 - t539;
t562 = -t412 * qJD(2) - t292;
t558 = -(t371 * t489 + t359) * qJD(2) - t292;
t553 = -t320 + t636;
t551 = t366 * t305 + t368 * t307;
t524 = rSges(5,2) * t582;
t550 = rSges(5,3) * t575 + t366 * t524;
t549 = rSges(5,3) * t572 + t368 * t524;
t525 = t370 * t611;
t548 = rSges(4,3) * t575 + t366 * t525;
t547 = rSges(4,3) * t572 + t368 * t525;
t546 = t659 * t366;
t545 = t659 * t368;
t533 = -m(4) - m(5) - m(6);
t527 = t370 * t613;
t526 = rSges(5,1) * t580;
t279 = (-rSges(5,1) * t354 - rSges(5,2) * t355) * t370;
t488 = -rSges(5,2) * t354 + t612;
t133 = qJD(4) * t279 + (t371 * t488 + t358) * qJD(2);
t523 = -t133 + t562;
t235 = -rSges(5,3) * t371 + t370 * t488;
t522 = -t235 + t560;
t504 = -t537 / 0.2e1;
t503 = t537 / 0.2e1;
t327 = t371 * t352;
t500 = -t369 * t370 + t327;
t502 = t500 * t366;
t501 = t500 * t368;
t498 = qJD(2) * t553;
t497 = t562 - t569;
t496 = t366 * t170 + t368 * t171 + t551;
t494 = -t557 + t560;
t468 = -Icges(3,5) * t370 - Icges(3,6) * t371;
t118 = rSges(5,1) * t259 - rSges(5,2) * t258 + rSges(5,3) * t576;
t120 = rSges(5,1) * t261 - rSges(5,2) * t260 + rSges(5,3) * t573;
t459 = t118 * t368 - t120 * t366;
t293 = -t365 * t575 - t367 * t368;
t294 = t366 * t574 - t577;
t168 = rSges(4,1) * t294 + rSges(4,2) * t293 + rSges(4,3) * t576;
t295 = -t365 * t572 + t366 * t367;
t296 = t367 * t572 + t578;
t169 = rSges(4,1) * t296 + rSges(4,2) * t295 + rSges(4,3) * t573;
t456 = t168 * t366 + t169 * t368;
t450 = qJD(2) * t662;
t237 = rSges(5,1) * t579 - rSges(5,2) * t581 + t358;
t41 = t118 * t308 - t120 * t309 + t414;
t429 = t41 * t459;
t418 = qJD(2) * t468;
t415 = t27 * t615 + t5 * t568;
t413 = t44 * t567 - t45 * t568;
t409 = qJD(2) * t562 + qJDD(2) * t560;
t408 = qJD(2) * t558 + qJDD(2) * t553;
t401 = Icges(4,5) * t371 + (-Icges(4,1) * t367 + Icges(4,4) * t365) * t370;
t398 = Icges(4,6) * t371 + (-Icges(4,4) * t367 + Icges(4,2) * t365) * t370;
t390 = qJD(2) * t401;
t389 = qJD(2) * t398;
t385 = t366 * t409 + t546;
t384 = t368 * t409 + t545;
t376 = (t27 * t568 - t44 * t557) * t368 + (-t27 * t567 + t45 * t557) * t366;
t342 = t368 * t539;
t340 = t366 * t539;
t306 = t321 * t368;
t304 = t321 * t366;
t299 = t468 * t368;
t298 = t468 * t366;
t285 = t368 * t418;
t284 = t366 * t418;
t251 = Icges(3,3) * t366 + t368 * t469;
t214 = t401 * t368;
t213 = t401 * t366;
t212 = t398 * t368;
t211 = t398 * t366;
t203 = -t368 * t526 + t549;
t201 = -t366 * t526 + t550;
t185 = t368 * t390;
t184 = t366 * t390;
t183 = t368 * t389;
t182 = t366 * t389;
t166 = -rSges(5,1) * t260 - rSges(5,2) * t261;
t162 = -rSges(5,1) * t258 - rSges(5,2) * t259;
t159 = Icges(4,1) * t296 + Icges(4,4) * t295 + Icges(4,5) * t573;
t158 = Icges(4,1) * t294 + Icges(4,4) * t293 + Icges(4,5) * t576;
t157 = Icges(4,4) * t296 + Icges(4,2) * t295 + Icges(4,6) * t573;
t156 = Icges(4,4) * t294 + Icges(4,2) * t293 + Icges(4,6) * t576;
t155 = Icges(4,5) * t296 + Icges(4,6) * t295 + Icges(4,3) * t573;
t154 = Icges(4,5) * t294 + Icges(4,6) * t293 + Icges(4,3) * t576;
t137 = t368 * t498 + t341;
t136 = t366 * t498 + t339;
t97 = rSges(5,1) * t179 + rSges(5,2) * t178 + rSges(5,3) * t517;
t95 = rSges(5,1) * t177 + rSges(5,2) * t176 + rSges(5,3) * t519;
t77 = -qJD(2) * t450 + qJDD(2) * t451 + qJDD(1);
t76 = t368 * t408 + t545;
t75 = t366 * t408 + t546;
t66 = qJD(2) * t456 + t449;
t65 = t118 * t537 + t235 * t309 + t437;
t64 = -t120 * t537 - t235 * t308 + t438;
t40 = qJD(2) * t666 + t456 * qJDD(2) + t404;
t29 = -t118 * t315 + t133 * t309 + t235 * t239 + t537 * t95 + t384;
t28 = t120 * t315 - t133 * t308 - t235 * t238 - t537 * t97 + t385;
t8 = t118 * t238 - t120 * t239 + t308 * t95 - t309 * t97 + t383;
t7 = -qJD(5) * t178 + qJDD(5) * t260 + t239 * t557 + t309 * t569 - t315 * t568 + t537 * t615 + t384;
t6 = -qJD(5) * t176 + qJDD(5) * t258 - t238 * t557 - t308 * t569 + t315 * t567 - t537 * t614 + t385;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) + t533) * g(3) + m(3) * t77 + m(4) * t40 + m(5) * t8 + m(6) * t5; (t690 * t366 - t368 * t716) * t630 + (t366 * t717 - t691 * t368) * t629 + (((t260 * t680 + t261 * t678 + t669) * t371 + t672) * qJD(4) + (((t673 + t690) * qJD(4) + t675) * t371 + t677) * t368 + (t260 * t685 + t261 * t683) * t309 + (t260 * t684 + t261 * t682) * t308) * t627 + (t366 * t693 - t368 * t694) * t626 + (((t258 * t680 + t259 * t678 + t670) * t371 + t671) * qJD(4) + (((t673 + t691) * qJD(4) + t675) * t371 + t677) * t366 + (t258 * t685 + t259 * t683) * t309 + (t258 * t684 + t259 * t682) * t308) * t625 + (t366 * t695 - t368 * t696) * t624 + (t366 * t688 - t368 * t689) * t623 - (t299 * qJD(2) * t362 + (t212 * t295 + t214 * t296) * t543 + (-t295 * t211 - t296 * t213 - t298 * t366) * t542) * t543 / 0.2e1 + (t298 * qJD(2) * t363 - (t211 * t293 + t213 * t294) * t542 + (t293 * t212 + t294 * t214 - t299 * t368) * t543) * t542 / 0.2e1 - t692 * t538 / 0.2e1 + (((t187 * t354 - t195 * t355 + t110) * t308 + (t186 * t354 - t194 * t355 + t109) * t309 + t80 * qJD(4)) * t370 + ((-t440 * t371 + (-t222 * t354 - t230 * t355 - t225) * t370 + t480) * qJD(4) + t633) * t371 + ((t193 * t354 - t197 * t355 + t108) * t308 + (t192 * t354 - t196 * t355 + t107) * t309 + t81 * qJD(4)) * t370 + ((t439 * t371 + (t228 * t354 - t232 * t355 - t223) * t370 + t479) * qJD(4) + t634) * t371) * t503 + (t5 * t496 + (t45 * t497 + t494 * t7 + t635) * t368 + (t44 * t497 + t494 * t6 + t415) * t366 - g(1) * (-t368 * t571 + t350) - g(2) * (-t366 * t571 + t347) - g(3) * t327 - (-g(3) * t369 + t637 * (-t354 * t663 - t355 * t664 - t352)) * t370 - t45 * (-t319 * t368 + t342) - t44 * (-t319 * t366 + t340) - (-t44 * t502 - t45 * t501) * qJD(2) - (t413 * t370 + (t44 * t563 + t45 * t564 + t376) * t371) * qJD(4) + t665 * (t579 * t664 + t581 * t663 + t360) + (-t308 * t564 - t309 * t563 - t371 * t534 + t674) * t27) * m(6) + (-g(1) * t549 - g(2) * t550 - g(3) * (t237 + t500) - t637 * (-t571 + (-t352 - t612) * t370) + t8 * t496 + (t8 * t120 + t29 * t522 + t523 * t65) * t368 + (t8 * t118 + t28 * t522 + t523 * t64) * t366 - t65 * (t237 * t309 + t342) - t64 * (-t237 * t308 + t340) - (-t501 * t65 - t502 * t64) * qJD(2) - ((-t118 * t65 + t120 * t64) * t370 + (t65 * (t235 * t366 + t201) + t64 * (-t235 * t368 - t203) + t429) * t371) * qJD(4) + (-t201 * t308 + t203 * t309 + t366 * t95 + t368 * t97 + t674) * t41) * m(5) + (-g(1) * (t344 + t547) - g(2) * (t343 + t548) - t637 * t370 * (-pkin(2) - t613) - t137 * t342 - t136 * t340 + t40 * t551 + (t137 * t558 + t40 * t169 + t553 * t76) * t368 + (t136 * t558 + t40 * t168 + t553 * t75) * t366 + (-(-t136 * t366 - t137 * t368) * qJD(2) - g(3)) * (rSges(4,1) * t574 - t371 * t611 + t359 + t645) + (-(t366 * (-t366 * t527 + t548) + t368 * (-t368 * t527 + t547)) * qJD(2) + t666 + t701) * t66) * m(4) + (g(1) * t306 + g(2) * t304 - g(3) * t323 + t77 * t451 + (-t450 - (-t304 * t366 - t306 * t368) * qJD(2)) * (qJD(2) * t451 + qJD(1)) + (qJD(2) ^ 2 * t451 + qJDD(2) * t662) * t321) * m(3) + (t697 + ((t183 * t295 + t185 * t296 + t285 * t366) * t366 + (-t182 * t295 - t184 * t296 - t284 * t366) * t368) * t700 + ((-t154 * t573 - t156 * t295 - t158 * t296) * t368 + (t155 * t573 + t157 * t295 + t159 * t296 + t251 * t366 - t713) * t366) * t699) * t366 / 0.2e1 - (t698 + ((-t182 * t293 - t184 * t294 + t284 * t368) * t368 + (t183 * t293 + t185 * t294 - t285 * t368) * t366) * t700 + ((-t154 * t576 - t156 * t293 - t158 * t294 + t713) * t368 + (t155 * t576 + t157 * t293 + t159 * t294 - t251 * t368) * t366) * t699) * t368 / 0.2e1 + (t366 * t655 - t368 * t656 + t704) * t504; -t533 * g(3) * t371 + 0.2e1 * (t27 * t628 + t5 * t620) * m(6) + 0.2e1 * (t41 * t628 + t620 * t8) * m(5) + 0.2e1 * (t40 * t620 + t628 * t66) * m(4) + (t533 * t637 + m(4) * (qJD(2) * t66 + t366 * t75 + t368 * t76) + m(5) * (qJD(2) * t41 + t28 * t366 + t29 * t368) + m(6) * (qJD(2) * t27 + t366 * t6 + t368 * t7)) * t370; (t639 * t370 - t371 * t714) * t630 + (t638 * t370 - t371 * t715) * t629 + (t260 * t643 + t261 * t642 - t368 * t641) * t627 + (t718 * t371 + (t366 * t694 + t368 * t693) * t370 + (t371 * t639 + t672) * qJD(2)) * t626 + (t258 * t643 + t259 * t642 - t366 * t641) * t625 + (t719 * t371 + (t366 * t696 + t368 * t695) * t370 + (t371 * t638 + t671) * qJD(2)) * t624 + (t370 * t640 - t371 * t653) * t623 + (t238 * t688 + t239 * t689 + t655 * t308 + t656 * t309 + t653 * t315 + t654 * t537) * t620 + t698 * t576 / 0.2e1 + t697 * t573 / 0.2e1 + t692 * t541 / 0.2e1 + (t654 * t371 + (t366 * t656 + t368 * t655) * t370 + (t370 * t653 + t371 * t640) * qJD(2)) * t504 + (t676 * t371 + (t354 * t643 + t355 * t642) * t370) * t503 + (g(1) * t565 - g(2) * t566 - g(3) * t552 + (qJD(2) * t376 - t44 * t614 + t45 * t615 - t567 * t6 + t568 * t7) * t371 + (t413 * qJD(2) + (-t44 * t569 - t557 * t6 + t415) * t368 + (t45 * t569 + t557 * t7 - t635) * t366) * t370 - (t259 * t44 + t261 * t45 + t27 * t580) * qJD(5) - (t27 * t565 + t45 * t552) * t309 - (t27 * t566 - t44 * t552) * t308 - (t44 * t565 + t45 * t566) * t537) * m(6) + (-g(1) * t166 - g(2) * t162 - g(3) * t279 + (t29 * t118 - t28 * t120 - t64 * t97 + t65 * t95 + (t429 + (t366 * t65 - t368 * t64) * t235) * qJD(2)) * t371 + (t65 * (-qJD(2) * t118 + t133 * t366) + t64 * (qJD(2) * t120 - t133 * t368) + t8 * t459 + t41 * (-t366 * t97 + t368 * t95) + (-t28 * t368 + t29 * t366) * t235) * t370 - t65 * (t162 * t537 + t279 * t309) - t64 * (-t166 * t537 - t279 * t308) - t41 * (t162 * t308 - t166 * t309)) * m(5) + t704 * t540 / 0.2e1; (-t176 * t44 - t178 * t45 + t416 * t27 + (t5 + t665) * t582 + (t27 * t309 + t44 * t537 - g(1) + t7) * t260 + (-t27 * t308 - t45 * t537 - g(2) + t6) * t258) * m(6);];
tau = t1;
