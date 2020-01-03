% Calculate vector of inverse dynamics joint torques for
% S5RPRPP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:20
% EndTime: 2019-12-31 18:12:53
% DurationCPUTime: 27.46s
% Computational Cost: add. (10733->690), mult. (14440->793), div. (0->0), fcn. (11128->6), ass. (0->370)
t772 = Icges(5,4) - Icges(4,5);
t771 = Icges(5,5) - Icges(4,6);
t770 = Icges(5,1) + Icges(4,3);
t378 = pkin(7) + qJ(3);
t360 = sin(t378);
t361 = cos(t378);
t769 = t360 * t771 - t361 * t772;
t623 = Icges(4,4) * t360;
t260 = Icges(4,2) * t361 + t623;
t344 = Icges(6,6) * t360;
t251 = -Icges(6,2) * t361 + t344;
t612 = Icges(5,6) * t360;
t449 = Icges(5,3) * t361 + t612;
t763 = t251 - t449;
t768 = t260 - t763;
t256 = Icges(6,4) * t360 + Icges(6,5) * t361;
t384 = sin(qJ(1));
t385 = cos(qJ(1));
t169 = Icges(6,1) * t384 + t256 * t385;
t739 = t384 * t770 + t769 * t385;
t767 = t169 + t739;
t766 = t770 * t385;
t599 = t360 * t384;
t318 = Icges(4,4) * t599;
t597 = t361 * t384;
t618 = Icges(4,5) * t385;
t159 = Icges(4,1) * t597 - t318 - t618;
t307 = Icges(6,6) * t599;
t616 = Icges(6,5) * t385;
t162 = -Icges(6,3) * t597 - t307 + t616;
t765 = -t159 + t162;
t263 = Icges(4,1) * t361 - t623;
t160 = Icges(4,5) * t384 + t263 * t385;
t447 = Icges(6,3) * t361 + t344;
t161 = Icges(6,5) * t384 + t385 * t447;
t755 = t160 + t161;
t611 = Icges(5,6) * t361;
t450 = -Icges(5,3) * t360 + t611;
t163 = Icges(5,5) * t384 - t385 * t450;
t596 = t361 * t385;
t309 = Icges(6,6) * t596;
t598 = t360 * t385;
t620 = Icges(6,4) * t384;
t165 = Icges(6,2) * t598 + t309 + t620;
t754 = t163 + t165;
t617 = Icges(5,5) * t385;
t164 = Icges(5,6) * t597 - Icges(5,3) * t599 + t617;
t308 = Icges(6,6) * t597;
t619 = Icges(6,4) * t385;
t166 = -Icges(6,2) * t599 - t308 + t619;
t764 = t164 + t166;
t732 = t597 * t772 - t599 * t771 + t766;
t345 = Icges(4,4) * t361;
t262 = Icges(4,1) * t360 + t345;
t610 = Icges(6,6) * t361;
t753 = t610 - t611 + (-Icges(5,2) - Icges(6,3)) * t360;
t762 = -t262 + t753;
t613 = Icges(4,6) * t385;
t157 = Icges(4,4) * t597 - Icges(4,2) * t599 - t613;
t310 = Icges(5,6) * t599;
t621 = Icges(5,4) * t385;
t168 = Icges(5,2) * t597 - t310 + t621;
t761 = t157 * t360 - t168 * t361;
t760 = -t160 * t597 - t163 * t599;
t758 = t159 * t361 - t164 * t360 - t761;
t757 = (-Icges(6,4) - Icges(5,5)) * t361 + (-Icges(5,4) + Icges(6,5)) * t360;
t250 = Icges(6,2) * t360 + t610;
t457 = -Icges(4,2) * t360 + t345;
t745 = -t250 + t450 + t457;
t254 = Icges(4,5) * t360 + Icges(4,6) * t361;
t744 = t254 + t757;
t453 = Icges(5,2) * t361 - t612;
t756 = t263 + t447 + t453;
t483 = -t161 * t597 - t165 * t599 + t169 * t385;
t158 = Icges(4,6) * t384 + t385 * t457;
t311 = Icges(5,6) * t598;
t622 = Icges(5,4) * t384;
t167 = -Icges(5,2) * t596 + t311 + t622;
t752 = -t385 * t739 - t760;
t704 = -t158 * t599 - t167 * t597 + t752;
t672 = -t483 + t704;
t751 = t384 * t767 + t596 * t755 + t598 * t754;
t624 = Icges(6,1) * t385;
t170 = -Icges(6,4) * t599 - Icges(6,5) * t597 + t624;
t146 = t384 * t170;
t750 = t384 * t732 + t596 * t765 + t598 * t764 + t146;
t749 = t157 + t764;
t748 = t158 - t754;
t747 = t168 - t765;
t746 = -t167 + t755;
t742 = t762 * qJD(3);
t741 = t768 * qJD(3);
t431 = t260 * t360 - t262 * t361;
t734 = t360 * t763 - t361 * t753 - t431;
t737 = t757 * t385;
t703 = t597 * t753 - t599 * t763 + t737;
t600 = t254 * t385;
t72 = -t384 * t431 - t600;
t740 = t72 - t703;
t738 = t158 * t360 + t167 * t361;
t594 = t385 * t170;
t442 = t162 * t361 + t166 * t360;
t685 = t384 * t442;
t60 = t594 - t685;
t673 = t384 * t758 + t385 * t732 + t60;
t671 = -t157 * t598 + t168 * t596 - t750;
t670 = -t158 * t598 - t167 * t596 + t751;
t736 = t756 * qJD(3);
t735 = t745 * qJD(3);
t726 = t256 + t769;
t733 = t762 * t360 - t361 * t768;
t666 = t744 * t384;
t731 = -t251 + t260;
t687 = t385 * t734 + t666;
t669 = t360 * t747 + t361 * t749;
t668 = t360 * t746 + t361 * t748;
t730 = -t742 * t384 + (-t385 * t453 + t622 - t755) * qJD(1);
t729 = t742 * t385 + (-t384 * t756 + t616 + t618 - t621) * qJD(1);
t728 = t741 * t384 + (t250 * t385 - t158 + t163 + t620) * qJD(1);
t727 = t741 * t385 + (t384 * t745 - t613 + t617 + t619) * qJD(1);
t725 = t744 * qJD(3);
t724 = t442 - t758;
t723 = t360 * t754 + t361 * t755 - t738;
t640 = rSges(6,1) + pkin(4);
t722 = qJD(1) * t744 + qJD(3) * t733 - t360 * t735 + t361 * t736;
t721 = t384 * t670 - t385 * t671;
t720 = t384 * t672 - t385 * t673;
t719 = t767 * qJD(1);
t718 = t734 * qJD(1) - qJD(3) * t726;
t717 = (Icges(5,3) * t596 + t385 * t731 + t311 - t746) * t384 + (t307 - t310 - t318 + (-Icges(4,2) - Icges(6,2) - Icges(5,3)) * t597 + t747) * t385;
t716 = t385 ^ 2;
t702 = rSges(6,3) + qJ(5);
t715 = t687 * qJD(1);
t714 = t730 * t361 - t728 * t360 + t669 * qJD(3) + (t170 + t732) * qJD(1);
t713 = -qJD(3) * t668 + t360 * t727 + t361 * t729 + t719;
t712 = (Icges(6,3) * t599 - t308 + t749) * t385 + (-Icges(6,3) * t598 + t309 - t748) * t384;
t711 = t745 - t762;
t710 = -t449 - t731 + t756;
t707 = t740 * qJD(1);
t706 = -t725 * t385 + (-t384 * t726 + t624 - t723 + t766) * qJD(1);
t705 = qJD(1) * t724 - t384 * t725 + t719;
t569 = -rSges(6,2) * t599 + t385 * t640 - t597 * t702;
t697 = qJ(5) * t596 + t384 * t640;
t696 = t732 + t738;
t526 = qJD(3) * qJD(4);
t695 = qJDD(4) * t360 + t361 * t526;
t536 = qJD(3) * t385;
t511 = t360 * t536;
t539 = qJD(1) * t384;
t694 = t361 * t539 + t511;
t693 = qJD(3) * t721 + t715;
t692 = qJD(3) * t720 + t707;
t691 = -t384 * t718 + t385 * t722;
t690 = t384 * t722 + t385 * t718;
t689 = qJD(3) * t723 + t360 * t729 - t361 * t727;
t688 = qJD(3) * t724 + t360 * t730 + t361 * t728;
t686 = -t594 + t751;
t343 = t360 * qJ(4);
t681 = pkin(3) * t361 + t343;
t221 = t681 * t384;
t366 = t385 * qJ(2);
t296 = pkin(1) * t384 - t366;
t382 = cos(pkin(7));
t350 = pkin(2) * t382 + pkin(1);
t383 = -pkin(6) - qJ(2);
t353 = t385 * t383;
t553 = -t350 * t384 - t353;
t153 = t296 + t553;
t277 = qJD(1) * t296;
t684 = qJD(1) * t153 - t277;
t301 = qJ(4) * t597;
t218 = -pkin(3) * t599 + t301;
t219 = rSges(5,2) * t599 + rSges(5,3) * t597;
t683 = t218 + t219;
t346 = t360 * rSges(5,3);
t632 = rSges(5,2) * t361;
t477 = t346 - t632;
t347 = t360 * rSges(6,2);
t682 = rSges(6,3) * t361 + t347;
t510 = t361 * t536;
t538 = qJD(1) * t385;
t680 = rSges(6,2) * t510 + t538 * t640;
t365 = t384 * qJ(2);
t298 = pkin(1) * t385 + t365;
t264 = pkin(3) * t360 - qJ(4) * t361;
t475 = rSges(6,2) * t361 - rSges(6,3) * t360;
t607 = qJ(5) * t360;
t498 = -t475 + t607;
t485 = -t264 - t498;
t529 = qJD(5) * t385;
t292 = t361 * t529;
t533 = qJD(4) * t385;
t294 = t360 * t533;
t363 = qJD(2) * t384;
t555 = t294 + t363;
t515 = t292 + t555;
t678 = t485 * t536 + t515;
t524 = -pkin(3) - t702;
t677 = t361 * t524 - t343 - t347 - t350;
t675 = qJD(1) * t569;
t546 = t384 ^ 2 + t716;
t674 = qJD(3) * t546;
t667 = -qJD(1) * t221 + t684;
t665 = t600 + t737;
t606 = qJ(5) * t361;
t664 = t682 + t681 + t606;
t381 = sin(pkin(7));
t633 = rSges(3,2) * t381;
t635 = rSges(3,1) * t382;
t199 = t384 * rSges(3,3) + (-t633 + t635) * t385;
t663 = t360 * t717 + t361 * t712;
t662 = (-t360 * t711 + t361 * t710) * qJD(1);
t661 = t714 * t716 + (t706 * t384 + (-t705 + t713) * t385) * t384;
t660 = t705 * t716 + (t713 * t384 + (-t706 + t714) * t385) * t384;
t659 = t726 * qJD(1);
t532 = qJD(5) * t360;
t535 = qJD(4) * t361;
t180 = qJD(3) * t681 - t535;
t573 = -qJD(3) * t682 - t180;
t649 = -t532 + t573 + (-t606 + t664) * qJD(3);
t648 = m(5) / 0.2e1;
t647 = m(6) / 0.2e1;
t646 = -m(5) - m(6);
t527 = qJD(1) * qJD(3);
t274 = qJDD(3) * t384 + t385 * t527;
t645 = t274 / 0.2e1;
t275 = -qJDD(3) * t385 + t384 * t527;
t643 = t275 / 0.2e1;
t642 = t384 / 0.2e1;
t641 = -t385 / 0.2e1;
t639 = g(2) * t384;
t636 = pkin(1) - t350;
t634 = rSges(4,1) * t361;
t267 = rSges(4,1) * t360 + rSges(4,2) * t361;
t225 = t267 * t385;
t370 = t384 * rSges(4,3);
t182 = rSges(4,1) * t596 - rSges(4,2) * t598 + t370;
t364 = qJD(2) * t385;
t537 = qJD(3) * t384;
t493 = t350 * t385 - t383 * t384;
t154 = t493 - t298;
t581 = t154 + t298;
t66 = -t267 * t537 - t364 + (t182 + t581) * qJD(1);
t631 = t225 * t66;
t629 = t360 * rSges(5,2);
t373 = t384 * rSges(5,1);
t482 = -t267 * t536 + t363;
t181 = rSges(4,1) * t597 - rSges(4,2) * t599 - rSges(4,3) * t385;
t582 = t153 - t296;
t519 = -t181 + t582;
t65 = qJD(1) * t519 + t482;
t628 = t384 * t65;
t512 = t360 * t537;
t281 = qJ(5) * t512;
t530 = qJD(5) * t384;
t480 = t361 * t530 - t281;
t589 = t475 * t537 + t480 + (t385 * t682 + t697) * qJD(1);
t514 = t360 * t539;
t588 = -rSges(6,2) * t514 - t694 * t702 + t292 + t680;
t240 = qJD(1) * t298 - t364;
t339 = t383 * t539;
t583 = t339 - (-t385 * t636 - t365) * qJD(1) - t240;
t572 = -qJD(3) * t477 - t180;
t571 = rSges(6,2) * t598 + rSges(6,3) * t596 + t697;
t184 = -rSges(5,2) * t596 + rSges(5,3) * t598 + t373;
t336 = pkin(3) * t596;
t226 = qJ(4) * t598 + t336;
t570 = -t184 - t226;
t568 = t221 * t384 + t226 * t385;
t304 = qJ(4) * t596;
t223 = -pkin(3) * t598 + t304;
t534 = qJD(4) * t384;
t567 = qJD(1) * t223 + t361 * t534;
t142 = t199 + t298;
t566 = -t264 * t537 - t364;
t476 = rSges(5,3) * t361 + t629;
t559 = -t264 + t476;
t558 = -t681 - t477;
t556 = rSges(4,2) * t514 + rSges(4,3) * t538;
t522 = t384 * t635;
t340 = t384 * t633;
t549 = rSges(3,3) * t385 + t340;
t198 = t522 - t549;
t554 = -t296 - t198;
t224 = rSges(5,2) * t598 + rSges(5,3) * t596;
t551 = rSges(3,3) * t538 + qJD(1) * t340;
t550 = t339 + t364;
t528 = qJD(1) * qJD(2);
t548 = qJDD(2) * t384 + t385 * t528;
t354 = qJ(2) * t538;
t547 = t354 + t363;
t342 = qJD(4) * t360;
t531 = qJD(5) * t361;
t282 = qJ(4) * t510;
t92 = -pkin(3) * t694 - qJ(4) * t514 + t282 + t294;
t196 = t360 * t538 + t361 * t537;
t289 = pkin(3) * t512;
t508 = t360 * t534;
t481 = -t289 + t508;
t93 = qJ(4) * t196 + qJD(1) * t336 + t481;
t523 = t221 * t538 + t384 * t93 + t385 * t92;
t520 = -t93 + t583;
t518 = -t221 + t582;
t517 = t218 * t537 + t223 * t536 + t342;
t516 = -t226 - t571;
t507 = -pkin(1) - t635;
t503 = -t537 / 0.2e1;
t502 = t537 / 0.2e1;
t501 = -t536 / 0.2e1;
t500 = t536 / 0.2e1;
t499 = rSges(5,1) * t385 - rSges(5,3) * t599;
t492 = qJD(3) * t572;
t186 = rSges(5,2) * t597 + t499;
t491 = t186 + t518;
t490 = rSges(5,1) * t538 + rSges(5,2) * t694 + rSges(5,3) * t510;
t488 = t524 * t385;
t487 = g(1) * t385 + t639;
t478 = t221 * t537 + t226 * t536 - t535;
t299 = rSges(2,1) * t385 - rSges(2,2) * t384;
t297 = rSges(2,1) * t384 + rSges(2,2) * t385;
t271 = -rSges(4,2) * t360 + t634;
t461 = -t384 * t66 - t385 * t65;
t460 = t275 * t264 + t385 * t695 + t548;
t459 = t518 + t569;
t112 = -rSges(4,1) * t694 - rSges(4,2) * t510 + t556;
t220 = t267 * t384;
t113 = -qJD(3) * t220 + (t271 * t385 + t370) * qJD(1);
t446 = t112 * t385 + t113 * t384;
t437 = t181 * t384 + t182 * t385;
t430 = -t531 - t342;
t428 = t493 + t226;
t427 = -qJDD(2) * t385 + qJD(1) * (-pkin(1) * t539 + t547) + qJDD(1) * t298 + t384 * t528;
t424 = -qJDD(4) * t361 + t221 * t274 + t360 * t526 + t536 * t92 + t537 * t93;
t400 = (-qJ(5) * qJD(3) ^ 2 + qJDD(5)) * t361 + (-0.2e1 * t532 + t573) * qJD(3);
t410 = qJD(1) * (-t354 + (t384 * t636 - t353) * qJD(1)) + qJDD(1) * t154 + t427;
t401 = qJDD(1) * t226 + t410 + t695 * t384 + (t294 + t92) * qJD(1);
t2 = t401 + t485 * t274 + t571 * qJDD(1) + (t292 + t588) * qJD(1) + t400 * t384;
t3 = t498 * t275 + t400 * t385 + t459 * qJDD(1) + (t384 * t430 + t520 - t589) * qJD(1) + t460;
t43 = t532 + (-t384 * t569 + t385 * t571) * qJD(3) + t478;
t414 = qJD(3) * t43 + t2 * t384 + t3 * t385;
t412 = t536 * t559 + t555;
t406 = -t350 - t681 - t346;
t333 = rSges(6,2) * t596;
t327 = rSges(6,2) * t597;
t295 = t361 * t533;
t243 = t271 * qJD(3);
t230 = t264 * t539;
t227 = -rSges(6,3) * t598 + t333;
t222 = -rSges(6,3) * t599 + t327;
t197 = t510 - t514;
t195 = t360 * t674;
t123 = qJD(1) * t142 - t364;
t122 = qJD(1) * t554 + t363;
t117 = -rSges(5,3) * t514 + t490;
t115 = t476 * t537 + (t385 * t477 + t373) * qJD(1);
t88 = t437 * qJD(3);
t64 = qJDD(1) * t199 + qJD(1) * (-qJD(1) * t522 + t551) + t427;
t63 = t554 * qJDD(1) + (-qJD(1) * t199 - t240) * qJD(1) + t548;
t48 = (t184 * t385 - t186 * t384) * qJD(3) + t478;
t47 = (qJD(3) * t476 + t342) * t384 + (-t570 + t581) * qJD(1) + t566;
t46 = qJD(1) * t491 + t412;
t36 = -t281 + (qJD(3) * t475 - t430) * t384 + (-t516 + t581) * qJD(1) + t566;
t35 = qJD(1) * t459 + t678;
t26 = qJD(1) * t112 + qJDD(1) * t182 - t243 * t537 - t267 * t274 + t410;
t25 = -t243 * t536 + t267 * t275 + t519 * qJDD(1) + (-t113 + t583) * qJD(1) + t548;
t6 = -t186 * t274 + t570 * t275 + (t115 * t384 + t117 * t385) * qJD(3) + t424;
t5 = qJD(1) * t117 + qJDD(1) * t184 + t274 * t559 + t384 * t492 + t401;
t4 = -t476 * t275 + t385 * t492 + t491 * qJDD(1) + (-t115 - t508 + t520) * qJD(1) + t460;
t1 = qJDD(5) * t360 - t569 * t274 + t516 * t275 + (t384 * t589 + t385 * t588 + t531) * qJD(3) + t424;
t7 = [-m(2) * (-g(1) * t297 + g(2) * t299) - t703 * t275 / 0.2e1 + (((t60 + t685 + t686) * t384 + ((t739 + t761) * t385 + t704 + t750 + t760) * t385) * qJD(3) + t715) * t500 + (t734 * qJD(3) + t736 * t360 + t735 * t361) * qJD(1) + (-(-t35 + t667 + t675 + t678) * t36 + t35 * (-t480 - t481 + t550) + t36 * (t282 + t515 + t680) + (t35 * (-rSges(6,2) - qJ(4)) * t597 + (rSges(6,3) * t35 * t384 + t36 * t488) * t360) * qJD(3) + ((t35 * t677 - t36 * t383) * t385 + (-t35 * t640 + t36 * t677) * t384) * qJD(1) + (t2 - g(2)) * (t428 + t571) + (t3 - g(1)) * (-t221 + t553 + t569)) * m(6) + (-(qJD(1) * t186 + t412 - t46 + t667) * t47 + t46 * (t289 + t550) + t47 * (-pkin(3) * t511 + t282 + t490 + t555) + t46 * (-t342 + (-t629 + (-rSges(5,3) - qJ(4)) * t361) * qJD(3)) * t384 + ((-rSges(5,1) * t46 + t406 * t47) * t384 + (t46 * (t406 + t632) - t47 * t383) * t385) * qJD(1) + (t5 - g(2)) * (t184 + t428) + (t4 - g(1)) * ((-t343 + (rSges(5,2) - pkin(3)) * t361) * t384 + t499 + t553)) * m(5) + (t65 * t550 + t66 * (t363 + t556) + (t267 * t628 - t631) * qJD(3) + ((-t65 * rSges(4,3) + t66 * (-t350 - t634)) * t384 + (t65 * (-t271 - t350) - t66 * t383) * t385) * qJD(1) - (-qJD(1) * t181 + t482 - t65 + t684) * t66 + (t26 - g(2)) * (t182 + t493) + (t25 - g(1)) * (-t181 + t553)) * m(4) + (-(-qJD(1) * t198 - t122 - t277 + t363) * t123 + t122 * t364 + t123 * (t547 + t551) + (t122 * (t507 + t633) * t385 + (t122 * (-rSges(3,3) - qJ(2)) + t123 * t507) * t384) * qJD(1) + (t64 - g(2)) * t142 + (t63 - g(1)) * (t384 * t507 + t366 + t549)) * m(3) + (t72 + t669) * t643 + (t668 + t687) * t645 + (t689 + t691) * t502 + (m(2) * (t297 ^ 2 + t299 ^ 2) + Icges(3,2) * t382 ^ 2 + (Icges(3,1) * t381 + 0.2e1 * Icges(3,4) * t382) * t381 + Icges(2,3) - t733) * qJDD(1) + (((t385 * t696 + t670 - t686) * t385 + (t384 * t696 + t146 + t483 + t671 - t752) * t384) * qJD(3) + t692 - t707) * t503 + (-t688 + t690 + t693) * t501; (-m(3) - m(4) + t646) * (g(1) * t384 - g(2) * t385) + 0.2e1 * (t2 * t641 + t3 * t642) * m(6) + 0.2e1 * (t4 * t642 + t5 * t641) * m(5) + 0.2e1 * (t25 * t642 + t26 * t641) * m(4) + 0.2e1 * (t63 * t642 + t64 * t641) * m(3); t721 * t645 + t720 * t643 + (qJD(1) * t691 + t661 * qJD(3) + qJDD(1) * t687 + t670 * t274 + t671 * t275) * t642 + (t690 * qJD(1) + t660 * qJD(3) + qJDD(1) * t740 + t672 * t274 + t673 * t275) * t641 - ((t712 * t360 - t361 * t717) * qJD(3) + (t710 * t360 + t711 * t361) * qJD(1)) * qJD(1) / 0.2e1 + (t688 * t385 + t689 * t384 + (t669 * t384 + t668 * t385) * qJD(1)) * qJD(1) / 0.2e1 + (t668 * t384 - t669 * t385) * qJDD(1) / 0.2e1 + t692 * t539 / 0.2e1 + t693 * t538 / 0.2e1 + ((-t537 * t665 + t659) * t384 + ((t384 * t666 + t663) * qJD(3) + t662) * t385) * t503 + ((t384 * t671 + t385 * t670) * qJD(1) + t661) * t502 + ((t384 * t673 + t385 * t672) * qJD(1) + t660) * t501 + ((-t536 * t666 - t659) * t385 + ((t385 * t665 + t663) * qJD(3) + t662) * t384) * t500 + (t1 * t568 + (-t1 * t569 + t2 * t485) * t384 + (t1 * t571 + t3 * t485) * t385 - g(1) * (t304 + t333) - g(2) * (t301 + t327) - g(3) * t664 - (g(1) * t488 + t524 * t639) * t360 + (t360 * t530 - t567 + t649 * t384 + (qJ(5) * t598 + t385 * t485 - t227) * qJD(1)) * t36 + (t360 * t529 + t230 - t295 + t649 * t385 + (-t384 * t475 + t218 + t222) * qJD(1)) * t35 + (t523 + (qJD(1) * t516 + t589) * t384 + (t588 - t675) * t385 - t517 - t531 - (t222 * t384 + t227 * t385 - t546 * t607) * qJD(3)) * t43) * m(6) + (-g(1) * (t223 + t224) - g(2) * t683 + g(3) * t558 - t46 * (-qJD(1) * t683 + t295) - t47 * (qJD(1) * t224 + t567) - t48 * t517 - ((t224 * t48 + t46 * t558) * t385 + (t219 * t48 + t47 * t558) * t384) * qJD(3) + t46 * t230 + t6 * t568 + t48 * t523 + (t4 * t559 + t46 * t572 + t6 * t184 + t48 * t117 + (-t186 * t48 + t47 * t559) * qJD(1)) * t385 + (t5 * t559 + t47 * t572 - t6 * t186 + t48 * t115 + (-t46 * t476 + t48 * t570) * qJD(1)) * t384) * m(5) + (g(1) * t225 + g(2) * t220 - g(3) * t271 + (qJD(3) * t446 + t181 * t274 - t182 * t275) * t437 + t88 * ((t181 * t385 - t182 * t384) * qJD(1) + t446) + t461 * t243 + (-t25 * t385 - t26 * t384 + (-t385 * t66 + t628) * qJD(1)) * t267 - (t220 * t65 - t631) * qJD(1) - (t88 * (-t220 * t384 - t225 * t385) + t461 * t271) * qJD(3)) * m(4); t646 * (-g(3) * t361 + t360 * t487) - m(5) * (t195 * t48 + t196 * t47 + t197 * t46) - m(6) * (t195 * t43 + t196 * t36 + t197 * t35) + 0.2e1 * ((t46 * t536 + t47 * t537 - t6) * t648 + (t35 * t536 + t36 * t537 - t1) * t647) * t361 + 0.2e1 * ((qJD(3) * t48 + t384 * t5 + t385 * t4 - t46 * t539 + t47 * t538) * t648 + (-t35 * t539 + t36 * t538 + t414) * t647) * t360; ((t1 - g(3)) * t360 + (-t43 * t674 + t414 - t487) * t361) * m(6);];
tau = t7;
