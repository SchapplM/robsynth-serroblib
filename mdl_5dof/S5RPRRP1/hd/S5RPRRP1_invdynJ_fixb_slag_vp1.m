% Calculate vector of inverse dynamics joint torques for
% S5RPRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:59:14
% EndTime: 2019-12-05 17:59:56
% DurationCPUTime: 28.63s
% Computational Cost: add. (12022->736), mult. (17235->880), div. (0->0), fcn. (13424->6), ass. (0->387)
t781 = Icges(5,4) + Icges(6,4);
t780 = Icges(5,1) + Icges(6,1);
t779 = Icges(5,2) + Icges(6,2);
t395 = qJ(3) + qJ(4);
t371 = cos(t395);
t778 = t781 * t371;
t370 = sin(t395);
t777 = t781 * t370;
t776 = Icges(5,5) + Icges(6,5);
t775 = Icges(5,6) + Icges(6,6);
t774 = t779 * t371 + t777;
t768 = t780 * t370 + t778;
t772 = -t370 * t779 + t778;
t771 = t371 * t780 - t777;
t773 = Icges(5,3) + Icges(6,3);
t397 = sin(qJ(1));
t399 = cos(qJ(1));
t752 = t774 * t397 + t775 * t399;
t756 = -t775 * t397 + t774 * t399;
t766 = -t776 * t397 + t768 * t399;
t589 = t371 * t397;
t770 = t781 * t589;
t769 = t776 * t370 + t775 * t371;
t767 = t776 * t399;
t591 = t370 * t397;
t740 = t591 * t780 + t767 + t770;
t765 = -t775 * t370 + t776 * t371;
t763 = -t771 + t774;
t762 = t768 + t772;
t761 = t772 * t399 + t766;
t760 = -t771 * t399 + t756;
t759 = -t771 * t397 + t752;
t758 = t769 * t397 + t773 * t399;
t757 = -t773 * t397 + t769 * t399;
t754 = t771 * t370 + t772 * t371;
t674 = t766 * t370 + t756 * t371;
t708 = -t757 * t397 + t399 * t674;
t751 = t765 * t397;
t750 = t740 * t370 + t752 * t371;
t392 = qJD(3) + qJD(4);
t748 = t762 * t392;
t747 = t763 * t392;
t746 = -t766 * qJD(1) + t759 * t392;
t745 = t760 * t392 + (t768 * t397 + t767) * qJD(1);
t308 = t392 * t397;
t744 = -t756 * qJD(1) - t308 * t772 - t740 * t392;
t743 = -t752 * qJD(1) + t761 * t392;
t742 = t765 * t399;
t710 = -t757 * t399 - t756 * t589 - t591 * t766;
t707 = t754 * t397 + t742;
t705 = t754 * t399 - t751;
t741 = t758 * t397;
t396 = sin(qJ(3));
t634 = pkin(3) * t396;
t302 = pkin(4) * t370 + t634;
t622 = rSges(6,2) * t371;
t464 = rSges(6,1) * t370 + t622;
t739 = t302 + t464;
t738 = t757 * qJD(1);
t737 = t758 * qJD(1);
t711 = t758 * t399 + t752 * t589 + t740 * t591;
t709 = -t750 * t399 + t741;
t734 = t754 * qJD(1) - t392 * t769;
t733 = t750 * qJD(1) + t751 * t392 + t738;
t732 = -t674 * qJD(1) - t742 * t392 + t737;
t731 = t745 * t370 - t743 * t371 + t738;
t730 = t746 * t370 + t744 * t371 + t737;
t729 = t765 * qJD(1) + t748 * t370 + t747 * t371;
t309 = t392 * t399;
t693 = -t762 * qJD(1) + t760 * t308 - t759 * t309;
t694 = (t591 * t779 - t740 - t770) * t309 + t761 * t308 + t763 * qJD(1);
t728 = -t693 * t370 + t694 * t371;
t727 = rSges(5,2) * t370;
t400 = -pkin(7) - pkin(6);
t391 = -qJ(5) + t400;
t724 = t739 * t399 + (-rSges(6,3) + t391) * t397;
t723 = t705 * qJD(1);
t722 = t707 * qJD(1) + t710 * t308;
t721 = t733 * t397 + t730 * t399;
t720 = t732 * t397 - t731 * t399;
t719 = -t730 * t397 + t733 * t399;
t718 = t731 * t397 + t732 * t399;
t717 = t711 * t309 + t722;
t716 = t708 * t308 + t709 * t309 - t723;
t715 = t734 * t397 + t729 * t399;
t714 = -t729 * t397 + t734 * t399;
t713 = t744 * t370 - t746 * t371;
t712 = t743 * t370 + t745 * t371;
t706 = rSges(4,2) * t396;
t398 = cos(qJ(3));
t534 = qJD(3) * t398;
t520 = pkin(3) * t534;
t632 = pkin(4) * t371;
t257 = t392 * t632 + t520;
t624 = rSges(6,1) * t371;
t286 = -rSges(6,2) * t370 + t624;
t387 = t399 * rSges(6,3);
t532 = qJD(5) * t397;
t533 = qJD(3) * t399;
t501 = t398 * t533;
t344 = pkin(3) * t501;
t537 = qJD(1) * t400;
t553 = t399 * t537 + t344;
t587 = t391 * t399;
t609 = -t257 * t399 - t286 * t309 + t532 + t553 + (t387 - t587 + (-t634 + t739) * t397) * qJD(1);
t704 = -t752 * t370 + t740 * t371;
t703 = t756 * t370 - t371 * t766;
t526 = t399 * t634;
t580 = t397 * t400 + t526 - t724;
t692 = -qJD(1) * t769 - t742 * t308 + t751 * t309;
t538 = qJD(1) * t399;
t689 = t371 * t308 + t370 * t538;
t583 = t397 * t398;
t351 = Icges(4,4) * t583;
t586 = t396 * t397;
t601 = Icges(4,5) * t399;
t223 = Icges(4,1) * t586 + t351 + t601;
t606 = Icges(4,4) * t398;
t459 = Icges(4,1) * t396 + t606;
t224 = -Icges(4,5) * t397 + t399 * t459;
t315 = -Icges(4,2) * t396 + t606;
t263 = t315 * t399;
t420 = t397 * (t224 + t263) - t399 * (-Icges(4,2) * t586 + t223 + t351);
t607 = Icges(4,4) * t396;
t456 = Icges(4,2) * t398 + t607;
t221 = Icges(4,6) * t399 + t397 * t456;
t222 = -Icges(4,6) * t397 + t399 * t456;
t317 = Icges(4,1) * t398 - t607;
t264 = t317 * t397;
t265 = t317 * t399;
t421 = t397 * (t222 - t265) - t399 * (t221 - t264);
t688 = -t421 * t396 + t420 * t398;
t558 = t315 + t459;
t559 = -t456 + t317;
t687 = (t396 * t558 - t398 * t559) * qJD(1);
t528 = qJD(1) * qJD(3);
t298 = qJDD(3) * t397 + t399 * t528;
t202 = qJD(4) * t538 + qJDD(4) * t397 + t298;
t216 = t464 * t392;
t355 = pkin(3) * t586;
t401 = qJD(3) ^ 2;
t529 = qJD(1) * qJD(2);
t552 = qJDD(2) * t397 + t399 * t529;
t633 = pkin(3) * t398;
t427 = t298 * t633 - t355 * t401 + t552;
t626 = pkin(6) + t400;
t250 = t397 * t626 + t526;
t373 = t399 * qJ(2);
t319 = pkin(1) * t397 - t373;
t631 = pkin(6) * t397;
t493 = -t319 - t631;
t471 = t250 + t493;
t429 = t471 - t580;
t629 = pkin(6) * qJD(1) ^ 2;
t491 = qJDD(5) - t629;
t630 = pkin(6) * t399;
t162 = (t355 - t630) * qJD(1) - t553;
t323 = t399 * pkin(1) + t397 * qJ(2);
t369 = qJD(2) * t399;
t252 = qJD(1) * t323 - t369;
t578 = -t162 - t252;
t592 = t370 * t392;
t12 = t202 * t286 - t216 * t308 + t491 * t399 + (t202 * t371 - t308 * t592) * pkin(4) + t429 * qJDD(1) + (-t532 + t578 - t609) * qJD(1) + t427;
t686 = t12 - g(1);
t365 = qJDD(3) * t399;
t539 = qJD(1) * t397;
t203 = qJDD(4) * t399 - t392 * t539 + t365;
t367 = qJD(5) * t399;
t502 = t397 * t534;
t342 = pkin(3) * t502;
t506 = t396 * t538;
t509 = pkin(3) * t506 + t397 * t537 + t342;
t161 = pkin(6) * t539 + t509;
t251 = -t399 * t626 + t355;
t299 = -t397 * t528 + t365;
t361 = qJDD(1) * t630;
t368 = qJD(2) * t397;
t551 = qJ(2) * t538 + t368;
t514 = qJD(1) * (-pkin(1) * t539 + t551) + qJDD(1) * t323 + t397 * t529;
t431 = -qJDD(2) * t399 + t514;
t404 = qJD(1) * t161 + qJDD(1) * t251 - t299 * t633 + t401 * t526 + t361 + t431;
t677 = rSges(6,1) * t591 + rSges(6,2) * t589 + t397 * t302 + t387;
t579 = -t355 + (-t391 + t400) * t399 + t677;
t507 = t371 * t538;
t673 = rSges(6,1) * t689 + rSges(6,2) * t507 + t397 * t257 + t302 * t538 + t391 * t539 + t367;
t608 = (-rSges(6,2) * t592 - rSges(6,3) * qJD(1)) * t397 - t509 + t673;
t13 = t491 * t397 + t579 * qJDD(1) + (t367 + t608) * qJD(1) + (-t203 * t371 + t309 * t592) * pkin(4) + t309 * t216 - t203 * t286 + t404;
t685 = t13 - g(2);
t243 = rSges(5,1) * t589 - rSges(5,2) * t591;
t465 = rSges(5,1) * t370 + rSges(5,2) * t371;
t388 = t399 * rSges(5,3);
t196 = rSges(5,1) * t591 + rSges(5,2) * t589 + t388;
t492 = t323 + t630;
t470 = t251 + t492;
t625 = rSges(5,1) * t371;
t287 = t625 - t727;
t594 = t287 * t309;
t75 = -t594 - t344 - t369 + (t196 + t470) * qJD(1);
t683 = (qJD(1) * t243 + t309 * t465) * t75;
t442 = t222 * t398 + t224 * t396;
t678 = t442 * t399;
t301 = qJD(1) * t319;
t676 = qJD(1) * t250 - t301;
t389 = t399 * rSges(4,3);
t227 = rSges(4,1) * t586 + rSges(4,2) * t583 + t389;
t675 = t227 + t492;
t324 = -rSges(3,2) * t399 + t397 * rSges(3,3);
t242 = rSges(6,1) * t589 - rSges(6,2) * t591;
t337 = pkin(4) * t589;
t310 = qJD(1) * t337;
t672 = -qJD(1) * t242 + t286 * t539 - t309 * t464 + t310;
t311 = pkin(4) * t507;
t519 = t308 * t370;
t671 = pkin(4) * t519 - t397 * t216 + t286 * t538 + t308 * t464 + t311;
t590 = t370 * t399;
t335 = rSges(6,2) * t590;
t588 = t371 * t399;
t244 = -rSges(6,1) * t588 + t335;
t670 = -t309 * t244 + t399 * t609;
t322 = rSges(4,1) * t398 - t706;
t267 = t322 * t399;
t466 = rSges(4,1) * t396 + rSges(4,2) * t398;
t141 = -qJD(3) * t267 + (t397 * t466 + t389) * qJD(1);
t294 = t466 * qJD(3);
t228 = -t397 * rSges(4,3) + t399 * t466;
t472 = t228 + t493;
t524 = t399 * t629;
t535 = qJD(3) * t397;
t63 = -t524 - t294 * t535 + t298 * t322 + (-t141 - t252) * qJD(1) + t472 * qJDD(1) + t552;
t505 = t398 * t538;
t510 = rSges(4,2) * t505 + (t502 + t506) * rSges(4,1);
t536 = qJD(3) * t396;
t142 = (-rSges(4,2) * t536 - rSges(4,3) * qJD(1)) * t397 + t510;
t523 = t397 * t629;
t64 = -t523 + qJD(1) * t142 + qJDD(1) * t227 - t299 * t322 + t361 + (qJD(3) * t294 - qJDD(2)) * t399 + t514;
t669 = t63 * t397 - t64 * t399;
t124 = -t594 + (t397 * t465 + t388) * qJD(1);
t217 = t465 * t392;
t384 = t397 * rSges(5,3);
t198 = t399 * t465 - t384;
t436 = t198 + t471;
t29 = -t524 + t202 * t287 - t217 * t308 + (-t124 + t578) * qJD(1) + t436 * qJDD(1) + t427;
t511 = rSges(5,1) * t689 + rSges(5,2) * t507;
t126 = (-rSges(5,2) * t592 - rSges(5,3) * qJD(1)) * t397 + t511;
t30 = qJD(1) * t126 + qJDD(1) * t196 - t203 * t287 + t217 * t309 + t404 - t523;
t668 = t29 * t397 - t30 * t399;
t130 = t222 * t396 - t224 * t398;
t137 = qJD(1) * t221 - qJD(3) * t263;
t139 = -qJD(3) * t265 + (t397 * t459 + t601) * qJD(1);
t453 = Icges(4,5) * t396 + Icges(4,6) * t398;
t220 = -Icges(4,3) * t397 + t399 * t453;
t543 = qJD(1) * t220;
t667 = qJD(3) * t130 + t137 * t398 + t139 * t396 + t543;
t290 = t456 * qJD(3);
t291 = t459 * qJD(3);
t313 = Icges(4,5) * t398 - Icges(4,6) * t396;
t437 = t315 * t396 - t317 * t398;
t666 = qJD(1) * t313 + qJD(3) * t437 + t290 * t398 + t291 * t396;
t138 = qJD(1) * t222 + t315 * t535;
t140 = qJD(1) * t224 + qJD(3) * t264;
t443 = t221 * t396 - t223 * t398;
t219 = Icges(4,3) * t399 + t397 * t453;
t544 = qJD(1) * t219;
t665 = qJD(3) * t443 - t138 * t398 - t140 * t396 + t544;
t664 = g(1) * t397 - g(2) * t399;
t662 = t711 - t708;
t393 = t397 ^ 2;
t651 = -pkin(1) - pkin(6);
t649 = t202 / 0.2e1;
t648 = t203 / 0.2e1;
t647 = t298 / 0.2e1;
t646 = t299 / 0.2e1;
t645 = -t308 / 0.2e1;
t644 = t308 / 0.2e1;
t643 = -t309 / 0.2e1;
t642 = t309 / 0.2e1;
t641 = t397 / 0.2e1;
t640 = -t399 / 0.2e1;
t639 = t399 / 0.2e1;
t638 = -rSges(6,1) - pkin(4);
t637 = rSges(3,2) - pkin(1);
t636 = -rSges(5,3) - pkin(1);
t635 = -rSges(6,3) - pkin(1);
t628 = -qJD(1) / 0.2e1;
t627 = qJD(1) / 0.2e1;
t621 = rSges(3,3) * t399;
t494 = t286 + t632;
t555 = t342 + t368;
t409 = t308 * t494 + t367 + t555;
t61 = qJD(1) * t429 + t409;
t616 = t399 * t61;
t468 = t308 * t287 + t555;
t74 = qJD(1) * t436 + t468;
t615 = t399 * t74;
t475 = t369 - t532;
t62 = -t344 - t494 * t309 + (t470 + t579) * qJD(1) - t475;
t614 = t62 * t216;
t611 = qJDD(1) / 0.2e1;
t467 = -t250 * t533 - t251 * t535;
t49 = -t308 * t579 + t309 * t580 + t467;
t598 = qJD(1) * t49;
t271 = t322 * t535;
t105 = qJD(1) * t472 + t271 + t368;
t597 = t105 * t399;
t593 = t313 * t399;
t260 = t397 * t313;
t581 = t580 * t399;
t569 = -t196 - t251;
t320 = rSges(3,2) * t397 + t621;
t557 = -t319 + t320;
t247 = t323 + t324;
t556 = t397 * t286 + t337;
t550 = rSges(3,2) * t539 + rSges(3,3) * t538;
t549 = t368 - t301;
t540 = qJD(1) * t453;
t438 = t315 * t398 + t317 * t396;
t132 = t399 * t438 - t260;
t530 = t132 * qJD(1);
t527 = -rSges(4,3) + t651;
t525 = pkin(4) * t592;
t356 = pkin(3) * t583;
t522 = rSges(5,1) * t588;
t521 = pkin(3) * t536;
t516 = -t251 - t579;
t85 = t399 * t219 + t221 * t583 + t223 * t586;
t86 = -t399 * t220 - t222 * t583 - t224 * t586;
t504 = t396 * t535;
t500 = -t539 / 0.2e1;
t499 = t538 / 0.2e1;
t498 = -t535 / 0.2e1;
t497 = t535 / 0.2e1;
t496 = -t533 / 0.2e1;
t495 = t533 / 0.2e1;
t336 = rSges(5,2) * t590;
t245 = t336 - t522;
t482 = -t243 * t308 + t309 * t245;
t476 = -qJD(1) * t245 - t308 * t465;
t303 = t632 + t633;
t325 = rSges(2,1) * t399 - rSges(2,2) * t397;
t321 = rSges(2,1) * t397 + rSges(2,2) * t399;
t444 = t221 * t398 + t223 * t396;
t413 = qJD(1) * t444 + qJD(3) * t260 + t543;
t414 = -qJD(1) * t442 - qJD(3) * t593 + t544;
t463 = (t413 * t397 + t399 * t665) * t399 + (t414 * t397 - t399 * t667) * t397;
t462 = (-t397 * t665 + t413 * t399) * t399 + (t397 * t667 + t414 * t399) * t397;
t461 = t397 * t86 + t399 * t85;
t199 = t397 * t219;
t87 = -t399 * t444 + t199;
t88 = -t220 * t397 + t678;
t460 = t397 * t88 + t399 * t87;
t106 = qJD(1) * t675 - t322 * t533 - t369;
t450 = t105 * t397 - t106 * t399;
t449 = t141 * t399 - t142 * t397;
t441 = -t227 * t397 - t228 * t399;
t435 = (-t399 ^ 2 - t393) * t520;
t348 = pkin(3) * t505;
t434 = -pkin(3) * t504 + t348;
t430 = t370 * t638 - t622;
t428 = t465 + t634;
t419 = -t161 * t535 + t162 * t533 - t299 * t250 - t251 * t298;
t410 = t438 * qJD(1) - t453 * qJD(3);
t403 = (t397 * t708 + t399 * t709) * t649 + (t397 * t710 + t399 * t711) * t648 + (t692 * t397 + t728 * t399) * t645 + (t721 * t399 + t720 * t397 + (-t397 * t709 + t399 * t708) * qJD(1)) * t644 + (-t728 * t397 + t692 * t399) * t643 + (t719 * t399 + t718 * t397 + (-t397 * t711 + t399 * t710) * qJD(1)) * t642 + (t715 * qJD(1) - t705 * qJDD(1) + t708 * t202 + t709 * t203 + t720 * t308 + t309 * t721) * t641 + (qJD(1) * t714 + qJDD(1) * t707 + t202 * t710 + t203 * t711 + t308 * t718 + t309 * t719) * t639 + (t370 * t694 + t371 * t693) * t628 + (t713 * t399 + t712 * t397 + (-t397 * t704 + t399 * t703) * qJD(1)) * t627 + (t397 * t703 + t399 * t704) * t611 + t717 * t500 + t716 * t499;
t269 = t397 * t303;
t266 = t322 * t397;
t207 = t399 * t250;
t206 = (-t303 + t633) * t399;
t205 = -t356 + t269;
t175 = t399 * t198;
t167 = qJD(1) * t247 - t369;
t166 = qJD(1) * t557 + t368;
t149 = t399 * t162;
t131 = t397 * t438 + t593;
t128 = t131 * qJD(1);
t127 = t441 * qJD(3);
t108 = t399 * t124;
t92 = qJD(1) * t550 + qJDD(1) * t324 + t431;
t91 = t557 * qJDD(1) + (-qJD(1) * t324 - t252) * qJD(1) + t552;
t71 = -t196 * t308 - t198 * t309 + t467;
t68 = -t397 * t666 + t410 * t399;
t67 = t410 * t397 + t399 * t666;
t66 = qJD(3) * t442 - t137 * t396 + t139 * t398;
t65 = -qJD(3) * t444 - t138 * t396 + t140 * t398;
t48 = qJD(3) * t460 - t530;
t47 = qJD(3) * t461 + t128;
t18 = t124 * t309 - t126 * t308 - t196 * t202 - t198 * t203 + t419;
t9 = -t202 * t579 + t203 * t580 - t308 * t608 + t309 * t609 + t419;
t1 = [(t128 + ((-t87 + t199 + t86) * t397 + (t88 - t678 + (t220 - t444) * t397 + t85) * t399) * qJD(3)) * t498 - m(2) * (-g(1) * t321 + g(2) * t325) - t298 * t132 / 0.2e1 + t130 * t647 - t705 * t202 / 0.2e1 + t703 * t649 + (t131 - t443) * t646 + ((t662 + t708) * t309 + t722) * t645 + (t530 + (t220 * t393 + (-t199 + t86 + (t220 + t444) * t399) * t399) * qJD(3) + t48) * t496 + (t65 + t68) * t495 + (t66 + t67 + t47) * t497 + (-qJD(3) * t438 + t290 * t396 - t291 * t398 + t747 * t370 - t748 * t371) * qJD(1) + (-(-t61 + (-t580 - t631) * qJD(1) + t409 + t676) * t62 + t61 * t475 + t62 * (-rSges(6,2) * t519 + t551 + t673) + (t286 * t392 + t257) * t616 + ((t391 + t635) * t616 + (t61 * (-qJ(2) - t739) + t62 * t635) * t397) * qJD(1) + t685 * (t323 - t587 + t677) + t686 * (-t319 + t724)) * m(6) + (-(-t74 + (t198 - t631) * qJD(1) + t468 + t676) * t75 + t74 * (-t309 * t727 + t392 * t522 + t369 + t553) + t75 * (-rSges(5,2) * t519 + t509 + t511 + t551) + (t636 * t615 + (t74 * (-qJ(2) - t428) + t75 * t636) * t397) * qJD(1) + (t30 - g(2)) * (-t399 * t400 + t196 + t323 + t355) + (-g(1) + t29) * (t373 - t384 + (-pkin(1) + t400) * t397 + t428 * t399)) * m(5) + (-(-t105 + t271 + (t228 - t631) * qJD(1) + t549) * t106 + t105 * (rSges(4,1) * t501 - t533 * t706 + t369) + t106 * (-rSges(4,2) * t504 + t510 + t551) + (t527 * t597 + (t105 * (-qJ(2) - t466) + t106 * t527) * t397) * qJD(1) + (t64 - g(2)) * t675 + (t63 - g(1)) * (t397 * t651 + t228 + t373)) * m(4) + (-(qJD(1) * t320 - t166 + t549) * t167 + t166 * t369 + t167 * (t550 + t551) + (t166 * t637 * t399 + (t166 * (-rSges(3,3) - qJ(2)) - t167 * pkin(1)) * t397) * qJD(1) + (t92 - g(2)) * t247 + (t91 - g(1)) * (t397 * t637 + t373 + t621)) * m(3) + (t704 + t707) * t648 + (((t750 + t757) * t399 + t674 * t397 + t710 - t741) * t309 + (t662 - t711) * t308 + t716 + t723) * t643 + (t713 + t714) * t642 + (t712 + t715 + t717) * t644 + (m(2) * (t321 ^ 2 + t325 ^ 2) - t437 + Icges(2,3) + Icges(3,1) + t771 * t371 - t772 * t370) * qJDD(1); (-m(3) - m(6)) * t664 + 0.2e1 * (t12 * t641 + t13 * t640) * m(6) + 0.2e1 * (t640 * t92 + t641 * t91) * m(3) + (-t664 + t668) * m(5) + (-t664 + t669) * m(4); t403 + ((-t535 * t593 - t540) * t397 + (t687 + (t397 * t260 + t688) * qJD(3)) * t399) * t498 + t460 * t647 + t461 * t646 + (t130 * t397 - t399 * t443) * t611 + (t397 * t66 + t399 * t65 + (t130 * t399 + t397 * t443) * qJD(1)) * t627 + (qJD(1) * t67 + qJD(3) * t463 - qJDD(1) * t132 + t298 * t88 + t299 * t87) * t641 + (qJD(1) * t68 + qJD(3) * t462 + qJDD(1) * t131 + t298 * t86 + t299 * t85) * t639 + t48 * t499 + t47 * t500 + ((-t85 * t397 + t86 * t399) * qJD(1) + t462) * t495 + ((-t396 * t559 - t398 * t558) * qJD(1) + (t396 * t420 + t398 * t421) * qJD(3)) * t628 + ((t260 * t533 - t540) * t399 + (-t687 + (-t399 * t593 - t688) * qJD(3)) * t397) * t496 + ((-t87 * t397 + t88 * t399) * qJD(1) + t463) * t497 + (t12 * (t356 + t556) + t9 * (-t207 + t581) + (t13 * (-t286 - t303) + t614 + t516 * t598) * t399 + (t9 * t516 + (t250 - t580) * t598) * t397 - g(1) * (t242 + t269) - g(2) * (t335 + (-t303 - t624) * t399) - g(3) * (t430 - t634) + (-qJD(1) * t205 + t672) * t62 + (t348 + (-t521 - t525) * t397 - (-t206 - t244) * qJD(1) - t434 + t671) * t61 + (t149 + (-t161 - t608) * t397 - t206 * t309 - (-t205 - t242) * t308 - t435 + t670) * t49) * m(6) + (-g(1) * (t243 + t356) - g(2) * (t336 + (-t625 - t633) * t399) + g(3) * t428 - t74 * (t434 + t476) - t683 - t71 * (t435 + t482) + t29 * t356 + t74 * t348 + t18 * (-t175 - t207) + t71 * (t108 + t149) + (t30 * (-t287 - t633) + t75 * t217 + (t74 * t287 + t569 * t71) * qJD(1)) * t399 + (t29 * t287 + t74 * (-t217 - t521) + t18 * t569 + t71 * (-t126 - t161) + (t75 * t287 + t71 * (t198 + t250)) * qJD(1)) * t397) * m(5) + (-(t105 * t267 + t106 * t266) * qJD(1) - (t127 * (-t266 * t397 - t267 * t399) - t450 * t466) * qJD(3) + (qJD(3) * t449 - t227 * t298 - t228 * t299) * t441 + t127 * ((-t227 * t399 + t228 * t397) * qJD(1) + t449) - t450 * t294 + ((t106 * t397 + t597) * qJD(1) + t669) * t322 - g(1) * t266 + g(2) * t267 + g(3) * t466) * m(4); t403 + (t12 * t556 + (-t13 * t494 - t579 * t598 + t614) * t399 - g(1) * (t242 + t337) - g(2) * (t588 * t638 + t335) - g(3) * t430 + (-t397 * t579 + t581) * t9 + (-t310 + t672) * t62 + (qJD(1) * t244 - t397 * t525 - t311 + t671) * t61 + (t242 * t308 - (-t308 * t397 - t309 * t399) * t632 + (-qJD(1) * t580 - t608) * t397 + t670) * t49) * m(6) + (-t476 * t74 - t683 + t18 * (-t196 * t397 - t175) - (t397 * t74 - t399 * t75) * t217 + ((t397 * t75 + t615) * qJD(1) + t668) * t287 - g(1) * t243 - g(2) * t245 + g(3) * t465 + (-t482 - t126 * t397 + t108 + (-t196 * t399 + t198 * t397) * qJD(1)) * t71) * m(5); ((t308 * t49 + t686) * t399 + (-t309 * t49 + t685) * t397) * m(6);];
tau = t1;
