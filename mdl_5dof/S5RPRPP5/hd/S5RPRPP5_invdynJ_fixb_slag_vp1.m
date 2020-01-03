% Calculate vector of inverse dynamics joint torques for
% S5RPRPP5
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP5_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP5_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:08
% EndTime: 2019-12-31 18:16:41
% DurationCPUTime: 30.31s
% Computational Cost: add. (5376->674), mult. (13404->761), div. (0->0), fcn. (10269->4), ass. (0->355)
t770 = Icges(6,4) + Icges(5,5);
t769 = Icges(5,1) + Icges(6,1);
t763 = Icges(6,2) + Icges(5,3);
t768 = -Icges(5,6) + Icges(6,6);
t366 = sin(qJ(1));
t368 = cos(qJ(1));
t367 = cos(qJ(3));
t365 = sin(qJ(3));
t592 = Icges(4,4) * t365;
t428 = Icges(4,2) * t367 + t592;
t767 = -t368 * t428 + (Icges(4,6) - Icges(5,6)) * t366;
t766 = t770 * t367;
t591 = Icges(4,4) * t367;
t431 = Icges(4,1) * t365 + t591;
t765 = -t368 * t431 + (Icges(5,4) + Icges(4,5)) * t366;
t764 = Icges(5,4) - Icges(6,5);
t565 = t367 * t368;
t567 = t365 * t368;
t752 = Icges(5,5) * t567 - Icges(5,3) * t565 + t767;
t309 = Icges(5,5) * t565;
t762 = -Icges(5,1) * t567 + t309 + t765;
t262 = Icges(4,1) * t367 - t592;
t761 = t770 * t365;
t735 = t769 * t367 + t761;
t755 = t262 + t735;
t746 = t769 * t365 - t766;
t311 = Icges(6,4) * t565;
t586 = Icges(6,5) * t366;
t158 = Icges(6,1) * t567 - t311 + t586;
t740 = -t158 + t762;
t256 = -Icges(4,2) * t365 + t591;
t745 = t763 * t365 + t766;
t759 = -t256 + t745;
t568 = t365 * t366;
t758 = t770 * t568;
t757 = t768 * t368;
t566 = t366 * t367;
t753 = t763 * t566 + t757 - t758;
t424 = Icges(4,5) * t365 + Icges(4,6) * t367;
t149 = Icges(4,3) * t368 + t366 * t424;
t427 = Icges(5,4) * t365 - Icges(5,6) * t367;
t153 = Icges(5,2) * t368 + t366 * t427;
t756 = -t149 - t153;
t750 = t746 * t366 + t764 * t368;
t754 = t762 * t365 + t752 * t367;
t150 = -Icges(4,3) * t366 + t368 * t424;
t585 = Icges(5,2) * t366;
t154 = Icges(5,4) * t567 - Icges(5,6) * t565 - t585;
t751 = -t150 - t154;
t748 = t745 * t366;
t747 = -t367 * t763 + t761;
t744 = t751 * t368 + t752 * t566 + t740 * t568;
t155 = Icges(4,6) * t368 + t366 * t428;
t743 = -t155 - t753;
t581 = Icges(6,6) * t366;
t152 = Icges(6,4) * t567 - Icges(6,2) * t565 + t581;
t742 = -t152 - t752;
t312 = Icges(4,4) * t566;
t588 = Icges(4,5) * t368;
t161 = Icges(4,1) * t568 + t312 + t588;
t741 = t161 + t750;
t739 = t755 * t368;
t738 = (-Icges(4,5) - t764) * t367 + (Icges(4,6) + t768) * t365;
t737 = -t428 + t747;
t736 = -t431 - t746;
t578 = Icges(6,3) * t366;
t146 = Icges(6,5) * t567 - Icges(6,6) * t565 + t578;
t677 = t146 * t368 + t152 * t566 + t744;
t54 = t368 * t149 + t155 * t566 + t161 * t568;
t421 = Icges(6,5) * t365 - Icges(6,6) * t367;
t145 = -Icges(6,3) * t368 + t366 * t421;
t691 = t368 * t153 + t568 * t750;
t695 = -t145 * t368 + t753 * t566 + t691;
t678 = t54 + t695;
t734 = t753 * t367;
t732 = t754 * t368;
t670 = t756 * t366 + t565 * t753;
t405 = t256 * t367 + t262 * t365;
t720 = t365 * t735 - t367 * t745 + t405;
t413 = t155 * t367 + t161 * t365;
t668 = t365 * t750 + t413;
t731 = -t146 + t154;
t693 = rSges(6,3) + qJ(5);
t685 = pkin(4) * t568 - t368 * t693;
t730 = rSges(6,1) * t568 + t685;
t202 = t262 * t366;
t729 = (t769 * t566 + t202 + t743 + t758) * t368 + (-t739 + t742) * t366;
t728 = -t421 + t427;
t676 = -t145 * t366 - t413 * t368 - t750 * t567 - t670;
t574 = t152 * t367;
t415 = -t158 * t365 + t574;
t675 = -t368 * t415 - t732 + (t146 + t751) * t366;
t197 = t256 * t368;
t508 = qJD(3) * t368;
t727 = -qJD(3) * t197 + t745 * t508 + (-t747 * t366 + t155 + t757) * qJD(1);
t672 = t743 * t365 + t741 * t367;
t509 = qJD(3) * t366;
t726 = -t256 * t509 + t748 * qJD(3) + (t747 * t368 + t581 + t767) * qJD(1);
t671 = t742 * t365 + t740 * t367;
t725 = -t739 * qJD(3) + (t366 * t431 + t588 + t750) * qJD(1);
t724 = qJD(3) * t202 + t735 * t509 + (t746 * t368 + t586 - t765) * qJD(1);
t669 = t738 * t366;
t723 = t737 * qJD(3);
t722 = t736 * qJD(3);
t721 = t759 * t365 + t755 * t367;
t719 = -t415 - t754;
t718 = -t668 - t734;
t717 = -t424 - t728;
t664 = t738 * t368;
t674 = t720 * t366 - t664;
t673 = t368 * t405 - t565 * t745 + t735 * t567 + t669;
t716 = (-t145 - t756) * qJD(1);
t715 = t368 ^ 2;
t712 = (Icges(4,2) * t568 - t312 - t741 + t748) * t368 + (-t567 * t763 + t197 - t309 - t311 - t740) * t366;
t714 = -t729 * t365 + t712 * t367;
t713 = t720 * qJD(1) + t717 * qJD(3);
t517 = qJD(1) * t150;
t711 = t517 - t669 * qJD(3) + (t728 * t368 - t578 - t585 - t718) * qJD(1);
t710 = -t719 * qJD(1) + t664 * qJD(3) + t716;
t709 = t675 * t366 + t676 * t368;
t708 = t677 * t366 + t678 * t368;
t707 = t736 + t759;
t706 = t737 + t755;
t705 = -t672 * qJD(3) - t724 * t365 + t726 * t367 + t716;
t704 = t738 * qJD(1) + t721 * qJD(3) + t722 * t365 + t723 * t367;
t703 = t731 * qJD(1) + t671 * qJD(3) + t725 * t365 + t727 * t367 + t517;
t701 = (t707 * t365 + t706 * t367) * qJD(1);
t700 = t674 * qJD(1);
t699 = t673 * qJD(1);
t358 = t367 * rSges(6,2);
t264 = -rSges(6,1) * t365 + t358;
t234 = t264 * qJD(3);
t502 = qJD(1) * qJD(3);
t241 = qJDD(3) * t366 + t368 * t502;
t338 = qJD(5) * t366;
t303 = qJ(4) * t565;
t209 = pkin(3) * t567 - t303;
t345 = t368 * qJ(2);
t267 = pkin(1) * t366 - t345;
t616 = pkin(6) * t366;
t471 = -t267 - t616;
t454 = t209 + t471;
t322 = rSges(6,2) * t565;
t696 = rSges(6,1) + pkin(4);
t541 = t693 * t366 + t696 * t567 - t322;
t399 = t454 + t541;
t274 = t368 * pkin(1) + t366 * qJ(2);
t342 = qJD(2) * t368;
t177 = qJD(1) * t274 - t342;
t506 = qJD(4) * t368;
t482 = t367 * t506;
t324 = pkin(3) * t568;
t504 = t367 * qJD(3);
t484 = t368 * t504;
t505 = t365 * qJD(3);
t481 = t368 * t505;
t512 = qJD(1) * t366;
t687 = t367 * t512 + t481;
t492 = pkin(3) * t484 + t687 * qJ(4);
t91 = qJD(1) * t324 + t482 - t492;
t401 = -t177 - t91 - t482;
t370 = qJD(1) ^ 2;
t615 = pkin(6) * t370;
t467 = -qJDD(5) - t615;
t271 = rSges(6,1) * t367 + rSges(6,2) * t365;
t617 = pkin(4) * t367;
t472 = t271 + t617;
t618 = pkin(4) * t365;
t497 = qJD(3) ^ 2 * t618;
t501 = qJDD(4) * t367;
t503 = qJD(1) * qJD(2);
t525 = qJDD(2) * t366 + t368 * t503;
t211 = t271 * t368;
t528 = -pkin(4) * t484 - t338;
t561 = -qJD(3) * t211 + t528 + (-t264 * t366 + t685) * qJD(1);
t344 = t367 * qJ(4);
t619 = pkin(3) * t365;
t263 = t344 - t619;
t339 = qJD(4) * t365;
t176 = qJD(3) * t263 + t339;
t270 = pkin(3) * t367 + qJ(4) * t365;
t480 = t366 * t505;
t646 = qJD(4) * t480 + t176 * t509 + t241 * t270;
t2 = t467 * t368 + t472 * t241 + (qJD(3) * t234 - t497 - t501) * t366 + t399 * qJDD(1) + (t401 + t338 - t561) * qJD(1) + t525 + t646;
t698 = g(1) - t2;
t242 = qJDD(3) * t368 - t366 * t502;
t457 = -t270 - t472;
t361 = t368 * pkin(6);
t511 = qJD(1) * t368;
t331 = qJ(2) * t511;
t341 = qJD(2) * t366;
t524 = t331 + t341;
t494 = qJD(1) * (-pkin(1) * t512 + t524) + qJDD(1) * t274 + t366 * t503;
t461 = qJDD(1) * t361 + t494;
t465 = -t176 - t339;
t340 = qJD(4) * t367;
t483 = t366 * t340;
t543 = -rSges(6,2) * t566 + t730;
t485 = t366 * t504;
t291 = pkin(4) * t485;
t464 = -qJD(5) * t368 + t291;
t486 = t367 * t511;
t488 = t365 * t511;
t686 = t485 + t488;
t645 = t686 * rSges(6,1) + rSges(6,2) * t480 + pkin(4) * t488 + t693 * t512;
t560 = rSges(6,2) * t486 - t464 - t645;
t204 = -t344 * t366 + t324;
t493 = t686 * pkin(3) + qJ(4) * t480;
t507 = qJD(4) * t366;
t92 = (-qJ(4) * t511 - t507) * t367 + t493;
t647 = qJD(1) * t92 + qJDD(1) * t204 + t368 * t501;
t3 = t467 * t366 + t543 * qJDD(1) + t457 * t242 + (-t483 - t560) * qJD(1) + (t497 - qJD(5) * qJD(1) - qJDD(2) + (-t234 + t465) * qJD(3)) * t368 + t461 + t647;
t697 = t3 - g(2);
t272 = rSges(5,1) * t367 + rSges(5,3) * t365;
t692 = (-qJD(3) * t272 + t340) * t368;
t359 = t368 * rSges(4,3);
t167 = rSges(4,1) * t568 + rSges(4,2) * t566 + t359;
t649 = t361 + t274;
t690 = t167 + t649;
t688 = t154 + t734;
t684 = t708 * qJD(3) + t700;
t683 = t709 * qJD(3) - t699;
t682 = t718 * qJD(3) + t726 * t365 + t724 * t367;
t681 = t719 * qJD(3) - t727 * t365 + t725 * t367;
t680 = t713 * t366 - t704 * t368;
t679 = t704 * t366 + t713 * t368;
t667 = t711 * t715 + (t703 * t366 + (-t705 + t710) * t368) * t366;
t666 = t705 * t715 + (t710 * t366 + (-t703 + t711) * t368) * t366;
t665 = t717 * qJD(1);
t244 = qJD(1) * t267;
t655 = qJD(1) * t209 - t244;
t526 = -t366 * rSges(5,2) - rSges(5,3) * t565;
t275 = -rSges(3,2) * t368 + t366 * rSges(3,3);
t357 = t367 * rSges(5,3);
t265 = -rSges(5,1) * t365 + t357;
t500 = pkin(3) + t696;
t650 = -t500 * t365 + t358;
t273 = rSges(4,1) * t367 - rSges(4,2) * t365;
t213 = t273 * t368;
t445 = rSges(4,1) * t365 + rSges(4,2) * t367;
t113 = -qJD(3) * t213 + (t366 * t445 + t359) * qJD(1);
t236 = t445 * qJD(3);
t403 = -t361 * t370 + t525;
t170 = -t366 * rSges(4,3) + t368 * t445;
t456 = t170 + t471;
t25 = -t236 * t509 + t241 * t273 + (-t113 - t177) * qJD(1) + t456 * qJDD(1) + t403;
t490 = t686 * rSges(4,1) + rSges(4,2) * t486;
t116 = (-rSges(4,2) * t505 - rSges(4,3) * qJD(1)) * t366 + t490;
t397 = -t366 * t615 + t461;
t26 = qJD(1) * t116 + qJDD(1) * t167 - t242 * t273 + (qJD(3) * t236 - qJDD(2)) * t368 + t397;
t648 = t25 * t366 - t26 * t368;
t363 = t366 ^ 2;
t630 = m(5) / 0.2e1;
t629 = m(6) / 0.2e1;
t628 = -m(5) - m(6);
t627 = -pkin(1) - pkin(6);
t625 = t241 / 0.2e1;
t624 = t242 / 0.2e1;
t623 = t366 / 0.2e1;
t621 = rSges(5,1) + pkin(3);
t620 = rSges(3,2) - pkin(1);
t614 = g(2) * t368;
t613 = t2 * t366;
t612 = t3 * t368;
t212 = t272 * t368;
t360 = t368 * rSges(5,2);
t112 = -qJD(3) * t212 + (-t265 * t366 + t360) * qJD(1);
t235 = t265 * qJD(3);
t169 = rSges(5,1) * t567 + t526;
t402 = t169 + t454;
t5 = t241 * t272 + (qJD(3) * t235 - t501) * t366 + t402 * qJDD(1) + (-t112 + t401) * qJD(1) + t403 + t646;
t611 = t5 * t366;
t491 = t686 * rSges(5,1) + rSges(5,3) * t480;
t115 = qJD(1) * t526 + t491;
t527 = rSges(5,1) * t568 + t360;
t166 = -rSges(5,3) * t566 + t527;
t530 = -t270 - t272;
t6 = qJDD(1) * t166 + t530 * t242 + (t115 - t483) * qJD(1) + (-qJDD(2) + (-t235 + t465) * qJD(3)) * t368 + t397 + t647;
t610 = t6 * t368;
t604 = rSges(3,3) * t368;
t218 = t273 * t509;
t66 = qJD(1) * t456 + t218 + t341;
t600 = t368 * t66;
t598 = -rSges(6,2) - qJ(4);
t597 = -rSges(5,3) - qJ(4);
t595 = -t115 - t92;
t557 = t366 * t176 + t270 * t511;
t210 = t270 * t368;
t544 = -t210 * t508 + t340;
t542 = -t166 - t204;
t205 = pkin(3) * t566 + qJ(4) * t568;
t206 = rSges(6,1) * t566 + rSges(6,2) * t568;
t540 = -t205 - t206;
t207 = rSges(5,1) * t566 + rSges(5,3) * t568;
t539 = -t205 - t207;
t538 = -t270 * t508 - t342;
t268 = rSges(3,2) * t366 + t604;
t531 = -t267 + t268;
t174 = t274 + t275;
t523 = rSges(3,2) * t512 + rSges(3,3) * t511;
t522 = -t209 * t508 + t339;
t521 = t341 - t244;
t520 = t363 + t715;
t499 = -rSges(5,2) + t627;
t498 = -rSges(4,3) + t627;
t496 = -t92 + t560;
t495 = -t204 - t543;
t477 = -t509 / 0.2e1;
t476 = t509 / 0.2e1;
t475 = -t508 / 0.2e1;
t474 = t508 / 0.2e1;
t473 = t693 + t627;
t469 = t598 * t367;
t468 = t597 * t367;
t462 = qJD(4) * t504 + qJDD(4) * t365 - t242 * t209 + t91 * t508;
t460 = t342 + t492;
t458 = t324 + t649;
t455 = t204 + t649;
t452 = t627 * t366 + t345;
t451 = g(1) * t366 - t614;
t449 = qJD(1) * t205 - t365 * t506;
t448 = t341 - t483;
t446 = qJD(1) * t210 + t263 * t509 + t365 * t507;
t276 = rSges(2,1) * t368 - rSges(2,2) * t366;
t269 = rSges(2,1) * t366 + rSges(2,2) * t368;
t67 = t690 * qJD(1) - t273 * t508 - t342;
t432 = t366 * t66 - t368 * t67;
t420 = t113 * t368 - t116 * t366;
t410 = -t167 * t366 - t170 * t368;
t400 = t270 * t509 + t448;
t398 = t272 * t509 + t400;
t372 = t331 + t448 + t493;
t371 = t271 * t509 + t400 + t464;
t325 = pkin(4) * t566;
t224 = t366 * t270;
t220 = t270 * t512;
t208 = t273 * t366;
t184 = t520 * t504;
t183 = t480 - t486;
t175 = t368 * t209;
t122 = qJD(1) * t174 - t342;
t121 = qJD(1) * t531 + t341;
t89 = t368 * t91;
t70 = t410 * qJD(3);
t63 = qJD(1) * t523 + qJDD(1) * t275 - qJDD(2) * t368 + t494;
t62 = t531 * qJDD(1) + (-qJD(1) * t275 - t177) * qJD(1) + t525;
t49 = (-t169 * t368 + t366 * t542) * qJD(3) + t522;
t48 = t692 + (t166 + t455) * qJD(1) + t538;
t47 = qJD(1) * t402 + t398;
t44 = (-qJD(3) * t271 + t340) * t368 + (t455 + t543) * qJD(1) + t528 + t538;
t43 = qJD(1) * t399 + t371;
t42 = (t366 * t495 - t368 * t541) * qJD(3) + t522;
t4 = -t169 * t242 + t542 * t241 + (t112 * t368 + t366 * t595) * qJD(3) + t462;
t1 = -t541 * t242 + t495 * t241 + (t496 * t366 + t561 * t368) * qJD(3) + t462;
t7 = [-m(2) * (-g(1) * t269 + g(2) * t276) - t673 * t241 / 0.2e1 + t671 * t625 + (((t54 + (-t145 + t415) * t368 + t675 + t691 + t732) * t368 + ((-t145 + t574) * t366 + (t150 - t668 + t688) * t368 - t670 - t676 + t744) * t366) * qJD(3) + t700) * t477 + (-t720 * qJD(3) - t723 * t365 + t722 * t367) * qJD(1) + (-(-t43 + (t541 - t616) * qJD(1) + t371 + t655) * t44 + t43 * (t460 - t528) + t44 * (t291 + t372 + t645) + (t43 * (rSges(6,1) * t504 + rSges(6,2) * t505 - t340) - t44 * qJD(5)) * t368 + ((t43 * t473 + t44 * t469) * t368 + (t44 * t627 + (-qJ(2) + t650) * t43) * t366) * qJD(1) + t697 * (t366 * t469 + t458 + t730) - t698 * (t366 * t473 + t500 * t567 - t303 - t322 + t345)) * m(6) + (-(-t47 + (t169 - t616) * qJD(1) + t398 + t655) * t48 + t47 * t460 + t48 * (t372 + t491) - t47 * t692 + ((t468 * t48 + t47 * t499) * t368 + (t47 * (-qJ(2) + t265 - t619) + t48 * t499) * t366) * qJD(1) + (t6 - g(2)) * (t366 * t468 + t458 + t527) + (t5 - g(1)) * (t567 * t621 - t303 + t452 + t526)) * m(5) + (t66 * (rSges(4,1) * t484 - rSges(4,2) * t481 + t342) + t67 * (-rSges(4,2) * t480 + t490 + t524) + (t498 * t600 + (t66 * (-qJ(2) - t445) + t67 * t498) * t366) * qJD(1) - (t218 - t66 + (t170 - t616) * qJD(1) + t521) * t67 + (t26 - g(2)) * t690 + (t25 - g(1)) * (t170 + t452)) * m(4) + (t121 * t342 + t122 * (t523 + t524) + (t121 * t620 * t368 + (t121 * (-rSges(3,3) - qJ(2)) - t122 * pkin(1)) * t366) * qJD(1) - (qJD(1) * t268 - t121 + t521) * t122 + (t63 - g(2)) * t174 + (t62 - g(1)) * (t366 * t620 + t345 + t604)) * m(3) + (t672 + t674) * t624 + ((t150 * t363 + ((-t146 + t688) * t366 + t691 - t695) * t366 + ((t150 + t668 + t731) * t368 + t670 + t677) * t368) * qJD(3) + t683 + t699) * t475 + (t679 + t682) * t474 + (m(2) * (t269 ^ 2 + t276 ^ 2) + Icges(2,3) + Icges(3,1) + t721) * qJDD(1) + (t680 + t681 + t684) * t476; t628 * t451 + 0.2e1 * (t613 / 0.2e1 - t612 / 0.2e1) * m(6) + 0.2e1 * (t611 / 0.2e1 - t610 / 0.2e1) * m(5) + (-t368 * t63 + 0.2e1 * t62 * t623 - t451) * m(3) + (-t451 + t648) * m(4); t709 * t625 + t708 * t624 + (t680 * qJD(1) + t666 * qJD(3) - t673 * qJDD(1) + t675 * t241 + t676 * t242) * t623 + (t679 * qJD(1) + t667 * qJD(3) + t674 * qJDD(1) + t677 * t241 + t678 * t242) * t368 / 0.2e1 - ((t712 * t365 + t729 * t367) * qJD(3) + (-t706 * t365 + t707 * t367) * qJD(1)) * qJD(1) / 0.2e1 + (t682 * t368 + t681 * t366 + (-t672 * t366 + t671 * t368) * qJD(1)) * qJD(1) / 0.2e1 + (t671 * t366 + t672 * t368) * qJDD(1) / 0.2e1 - t684 * t512 / 0.2e1 + t683 * t511 / 0.2e1 + ((t664 * t509 + t665) * t366 + ((-t669 * t366 + t714) * qJD(3) - t701) * t368) * t477 + ((-t676 * t366 + t675 * t368) * qJD(1) + t666) * t476 + ((-t669 * t508 + t665) * t368 + ((t664 * t368 - t714) * qJD(3) + t701) * t366) * t475 + ((-t678 * t366 + t677 * t368) * qJD(1) + t667) * t474 + (-g(1) * (t325 - t540) - (t365 * t598 - t367 * t500) * t614 - (t344 + t650) * g(3) + t2 * (t224 + t325) - t1 * t175 + (t1 * t495 + t2 * t271) * t366 + (-t1 * t541 + t3 * t457) * t368 + (-t446 + t557 + (-(t264 - t618) * qJD(3) - pkin(4) * t505 + t234) * t366) * t43 + (-t449 + t220 + (-(-t263 - t264) * qJD(3) - t176 - t234) * t368 + (t271 * t366 - t206) * qJD(1)) * t44 + (-t544 - (-t211 * t368 + t366 * t540 - t520 * t617) * qJD(3) + t89 + (t496 + (t209 + t541) * qJD(1)) * t366 + (qJD(1) * t495 + t561) * t368) * t42) * m(6) + (g(1) * t539 - g(3) * (-t365 * t621 + t344 + t357) - (t365 * t597 - t367 * t621) * t614 - t47 * (qJD(1) * t212 + t446) - t48 * (qJD(1) * t207 + t449) - t49 * t544 - ((t48 * (-t263 - t265) - t49 * t212) * t368 + (t47 * t265 + t49 * t539) * t366) * qJD(3) + t5 * t224 + t47 * t557 + t48 * t220 - t4 * t175 + t49 * t89 + (t5 * t272 + t47 * t235 + t4 * t542 + t49 * t595 + (t48 * t272 + t49 * (t169 + t209)) * qJD(1)) * t366 + (t6 * t530 + t48 * (-t176 - t235) - t4 * t169 + t49 * t112 + (t47 * t272 + t49 * t542) * qJD(1)) * t368) * m(5) + ((qJD(3) * t420 - t167 * t241 - t170 * t242) * t410 + t70 * ((-t167 * t368 + t170 * t366) * qJD(1) + t420) - t432 * t236 + ((t366 * t67 + t600) * qJD(1) + t648) * t273 - (t208 * t67 + t213 * t66) * qJD(1) - (t70 * (-t208 * t366 - t213 * t368) - t432 * t445) * qJD(3) - g(1) * t208 + g(2) * t213 + g(3) * t445) * m(4); t628 * (g(3) * t365 - t367 * t451) - m(5) * (t183 * t47 + t184 * t49 - t48 * t687) - m(6) * (t183 * t43 + t184 * t42 - t44 * t687) + 0.2e1 * ((t47 * t509 - t48 * t508 + t4) * t630 + (t43 * t509 - t44 * t508 + t1) * t629) * t365 + 0.2e1 * ((qJD(3) * t49 - t47 * t511 - t48 * t512 + t610 - t611) * t630 + (qJD(3) * t42 - t43 * t511 - t44 * t512 + t612 - t613) * t629) * t367; (-t697 * t366 + t698 * t368) * m(6);];
tau = t7;
