% Calculate vector of inverse dynamics joint torques for
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:24
% EndTime: 2019-12-31 20:49:56
% DurationCPUTime: 26.99s
% Computational Cost: add. (18404->668), mult. (15659->786), div. (0->0), fcn. (12160->8), ass. (0->375)
t412 = qJ(3) + pkin(8);
t403 = sin(t412);
t404 = cos(t412);
t415 = sin(qJ(3));
t417 = cos(qJ(3));
t807 = Icges(4,5) * t417 + Icges(5,5) * t404 - Icges(4,6) * t415 - Icges(5,6) * t403;
t806 = Icges(4,3) + Icges(5,3);
t805 = Icges(6,4) + Icges(5,5);
t804 = Icges(5,6) - Icges(6,6);
t413 = qJ(1) + qJ(2);
t405 = sin(t413);
t292 = Icges(6,4) * t404 + Icges(6,6) * t403;
t406 = cos(t413);
t469 = t292 * t406;
t210 = Icges(6,2) * t405 + t469;
t802 = t807 * t406;
t776 = t806 * t405 + t802;
t803 = t210 + t776;
t650 = Icges(5,4) * t403;
t293 = Icges(5,2) * t404 + t650;
t371 = Icges(6,5) * t403;
t800 = Icges(6,3) * t404 + t293 - t371;
t649 = Icges(6,5) * t404;
t295 = Icges(6,1) * t403 - t649;
t372 = Icges(5,4) * t404;
t791 = -Icges(5,1) * t403 - t295 - t372;
t617 = t405 * t417;
t618 = t405 * t415;
t622 = t404 * t405;
t624 = t403 * t405;
t777 = -Icges(4,5) * t617 - Icges(5,5) * t622 + Icges(4,6) * t618 + Icges(5,6) * t624 + t806 * t406;
t504 = Icges(6,1) * t404 + t371;
t213 = -Icges(6,4) * t406 + t405 * t504;
t333 = Icges(5,4) * t624;
t215 = Icges(5,1) * t622 - Icges(5,5) * t406 - t333;
t795 = -t213 - t215;
t472 = t504 * t406;
t214 = Icges(6,4) * t405 + t472;
t298 = Icges(5,1) * t404 - t650;
t473 = t298 * t406;
t216 = Icges(5,5) * t405 + t473;
t801 = t214 + t216;
t211 = Icges(5,4) * t622 - Icges(5,2) * t624 - Icges(5,6) * t406;
t225 = Icges(4,4) * t617 - Icges(4,2) * t618 - Icges(4,6) * t406;
t799 = t211 * t403 + t225 * t415;
t651 = Icges(4,4) * t415;
t353 = Icges(4,1) * t417 - t651;
t474 = t353 * t406;
t228 = Icges(4,5) * t405 + t474;
t798 = -t216 * t622 - t228 * t617;
t288 = Icges(6,3) * t403 + t649;
t205 = -Icges(6,6) * t406 + t288 * t405;
t797 = t205 - t211;
t621 = t404 * t406;
t332 = Icges(6,5) * t621;
t623 = t403 * t406;
t206 = Icges(6,6) * t405 + Icges(6,3) * t623 + t332;
t502 = -Icges(5,2) * t403 + t372;
t470 = t502 * t406;
t212 = Icges(5,6) * t405 + t470;
t796 = t206 - t212;
t793 = t288 - t502;
t790 = t298 + t504;
t411 = qJD(1) + qJD(2);
t789 = t800 * qJD(3) - t804 * t411;
t788 = t791 * qJD(3) + t805 * t411;
t365 = Icges(4,4) * t618;
t227 = Icges(4,1) * t617 - Icges(4,5) * t406 - t365;
t787 = t215 * t404 + t227 * t417 - t799;
t786 = Icges(4,5) * t415 + Icges(4,6) * t417 + t805 * t403 + t804 * t404;
t517 = -t206 * t624 + t210 * t406 - t214 * t622;
t407 = Icges(4,4) * t417;
t503 = -Icges(4,2) * t415 + t407;
t471 = t503 * t406;
t226 = Icges(4,6) * t405 + t471;
t785 = -t776 * t406 - t798;
t742 = -t212 * t624 - t226 * t618 + t785;
t718 = -t517 + t742;
t784 = t212 * t403 + t226 * t415;
t614 = t406 * t417;
t783 = t206 * t623 + t228 * t614 + t803 * t405 + t801 * t621;
t209 = -Icges(6,2) * t406 + t292 * t405;
t193 = t405 * t209;
t782 = -t205 * t623 - t227 * t614 + t777 * t405 + t795 * t621 - t193;
t350 = Icges(4,2) * t417 + t651;
t352 = Icges(4,1) * t415 + t407;
t771 = t350 * t415 - t352 * t417 + t800 * t403 + t791 * t404;
t749 = rSges(6,1) + pkin(4);
t707 = t403 * t749;
t620 = t405 * t411;
t781 = t406 * t789 - t620 * t793;
t616 = t406 * t411;
t780 = t288 * t616 + t405 * t789 - t411 * t470;
t779 = t406 * t788 - t620 * t790;
t778 = (-t472 - t473) * t411 - t788 * t405;
t775 = t793 * qJD(3);
t774 = t790 * qJD(3);
t640 = t209 * t406;
t497 = t205 * t403 + t213 * t404;
t706 = t405 * t497;
t78 = -t640 + t706;
t721 = t405 * t787 + t406 * t777 + t78;
t615 = t406 * t415;
t720 = -t211 * t623 - t225 * t615 - t782;
t719 = -t212 * t623 - t226 * t615 + t783;
t773 = t292 + t807;
t772 = (Icges(6,2) + t806) * t411 - t786 * qJD(3);
t770 = -t350 * t417 - t352 * t415 + t403 * t791 - t404 * t800;
t769 = -t206 * t403 - t228 * t417 - t404 * t801 + t784;
t768 = -t497 - t787;
t767 = -t226 * t417 - t228 * t415 - t403 * t801 + t404 * t796;
t712 = -t225 * t417 - t227 * t415 + t403 * t795 + t404 * t797;
t710 = t786 * t406;
t709 = t786 * t405;
t741 = rSges(6,3) + qJ(5);
t717 = -t771 * t405 - t710;
t716 = -t771 * t406 + t709;
t699 = t741 * t622;
t766 = -t797 * t406 + (-Icges(6,1) * t623 + t295 * t406 + t332 + t796) * t405;
t765 = -t791 - t793;
t764 = t790 - t800;
t763 = (Icges(5,2) * t622 + t333 + t795) * t406 + (-t293 * t406 + t801) * t405;
t695 = t404 * rSges(6,1) + t403 * rSges(6,3);
t762 = t404 * pkin(4) + t403 * qJ(5) + t695;
t312 = t503 * qJD(3);
t313 = t353 * qJD(3);
t761 = t770 * qJD(3) - t312 * t415 + t313 * t417 + t775 * t403 + t774 * t404 + t411 * t786;
t447 = Icges(4,6) * t411 - qJD(3) * t350;
t152 = t405 * t447 + t411 * t471;
t450 = Icges(4,5) * t411 - qJD(3) * t352;
t154 = t405 * t450 + t411 * t474;
t760 = t152 * t415 - t154 * t417 + t778 * t404 - t780 * t403 + (-t209 + t777) * t411 - t712 * qJD(3);
t151 = t406 * t447 - t503 * t620;
t153 = -t353 * t620 + t406 * t450;
t759 = t767 * qJD(3) - t151 * t415 + t153 * t417 + t781 * t403 + t779 * t404 + t803 * t411;
t758 = t719 * t405 - t720 * t406;
t757 = t718 * t405 - t721 * t406;
t756 = t773 * qJD(3) + t771 * t411;
t755 = (t469 + t768 + t802) * t411 + t772 * t405;
t754 = t772 * t406 + t769 * t411 - t773 * t620;
t581 = -t741 * t404 + t707;
t753 = t716 * t411;
t752 = t717 * t411;
t570 = qJD(3) * t411;
t271 = -qJDD(3) * t406 + t405 * t570;
t410 = qJDD(1) + qJDD(2);
t613 = t417 * qJD(3) ^ 2;
t563 = pkin(3) * t613;
t565 = qJD(5) * t404;
t589 = -qJD(3) * t762 + t565;
t436 = -t563 + qJDD(5) * t403 + (t565 + t589) * qJD(3);
t416 = sin(qJ(1));
t418 = cos(qJ(1));
t420 = qJD(1) ^ 2;
t466 = (-qJDD(1) * t416 - t418 * t420) * pkin(1);
t567 = qJD(4) * t411;
t666 = pkin(3) * t415;
t442 = qJDD(4) * t405 + t271 * t666 + t406 * t567 + t466;
t393 = t406 * rSges(6,2);
t596 = t762 * t405 - t393;
t396 = t406 * pkin(7);
t307 = pkin(2) * t405 - t396;
t414 = -qJ(4) - pkin(7);
t368 = t406 * t414;
t408 = t417 * pkin(3);
t398 = t408 + pkin(2);
t578 = -t405 * t398 - t368;
t200 = t307 + t578;
t602 = t200 - t307;
t522 = -t596 + t602;
t566 = qJD(5) * t403;
t537 = t405 * t566;
t395 = t405 * pkin(7);
t619 = t405 * t414;
t325 = t411 * t619;
t569 = qJD(3) * t415;
t542 = t405 * t569;
t577 = pkin(3) * t542 + qJD(4) * t406;
t547 = t325 + t577;
t663 = pkin(2) - t398;
t146 = (-t406 * t663 - t395) * t411 - t547;
t308 = t406 * pkin(2) + t395;
t261 = t308 * t411;
t610 = -t146 - t261;
t572 = qJD(3) * t405;
t544 = t404 * t572;
t545 = t403 * t572;
t700 = t749 * t545;
t734 = t403 * t616 + t544;
t740 = t405 * rSges(6,2) + pkin(4) * t621;
t611 = t537 + t734 * qJ(5) + rSges(6,3) * t544 - t700 + (t406 * t695 + t740) * t411;
t8 = t581 * t271 + t522 * t410 + (-t537 + t610 - t611) * t411 + t436 * t406 + t442;
t751 = t8 - g(1);
t270 = qJDD(3) * t405 + t406 * t570;
t326 = t406 * t566;
t358 = pkin(7) * t616;
t376 = qJD(4) * t405;
t539 = t406 * t569;
t145 = -pkin(3) * t539 - t358 + t376 + (t405 * t663 - t368) * t411;
t524 = t406 * t398 - t619;
t201 = t524 - t308;
t409 = t418 * pkin(1);
t667 = pkin(1) * t416;
t520 = qJDD(1) * t409 - t420 * t667;
t477 = t411 * (-pkin(2) * t620 + t358) + t410 * t308 + t520;
t437 = -qJDD(4) * t406 + t411 * t145 + t410 * t201 + t405 * t567 + t477;
t519 = -t581 - t666;
t595 = rSges(6,1) * t621 + t741 * t623 + t740;
t557 = t404 * t620;
t571 = qJD(3) * t406;
t462 = -t403 * t571 - t557;
t559 = t403 * t620;
t543 = t404 * t571;
t701 = rSges(6,2) * t616 + t741 * t543;
t612 = t749 * t462 - t741 * t559 + t326 + t701;
t9 = t595 * t410 + (t326 + t612) * t411 + t519 * t270 + t436 * t405 + t437;
t750 = t9 - g(2);
t549 = -rSges(5,1) * t545 - t734 * rSges(5,2);
t697 = rSges(5,1) * t621 + t405 * rSges(5,3);
t143 = t411 * t697 + t549;
t375 = t404 * rSges(5,1);
t304 = -rSges(5,2) * t403 + t375;
t283 = t304 * qJD(3);
t301 = rSges(5,1) * t403 + rSges(5,2) * t404;
t218 = rSges(5,1) * t622 - rSges(5,2) * t624 - t406 * rSges(5,3);
t553 = -t218 + t602;
t23 = t271 * t301 + (-qJD(3) * t283 - t563) * t406 + (-t143 + t610) * t411 + t553 * t410 + t442;
t748 = -g(1) + t23;
t568 = qJD(3) * t417;
t541 = t405 * t568;
t562 = rSges(4,2) * t615;
t548 = rSges(4,1) * t542 + rSges(4,2) * t541 + t411 * t562;
t696 = rSges(4,1) * t614 + t405 * rSges(4,3);
t156 = t411 * t696 - t548;
t662 = rSges(4,1) * t417;
t361 = -rSges(4,2) * t415 + t662;
t321 = t361 * qJD(3);
t359 = rSges(4,1) * t415 + rSges(4,2) * t417;
t574 = rSges(4,2) * t618 + t406 * rSges(4,3);
t229 = rSges(4,1) * t617 - t574;
t590 = -t229 - t307;
t59 = -t321 * t571 + t271 * t359 + (-t156 - t261) * t411 + t590 * t410 + t466;
t747 = -g(1) + t59;
t586 = rSges(5,2) * t559 + rSges(5,3) * t616;
t141 = rSges(5,1) * t462 - rSges(5,2) * t543 + t586;
t220 = -rSges(5,2) * t623 + t697;
t24 = -t283 * t572 + t141 * t411 + t220 * t410 - t270 * t301 + (-t270 * t415 - t405 * t613) * pkin(3) + t437;
t746 = -g(2) + t24;
t538 = t406 * t568;
t475 = rSges(4,3) * t616 + (t411 * t618 - t538) * rSges(4,2);
t155 = (-t411 * t617 - t539) * rSges(4,1) + t475;
t230 = -t562 + t696;
t60 = t155 * t411 + t230 * t410 - t270 * t359 - t321 * t572 + t477;
t745 = -g(2) + t60;
t260 = rSges(3,1) * t616 - rSges(3,2) * t620;
t305 = rSges(3,1) * t405 + rSges(3,2) * t406;
t744 = -t260 * t411 - t305 * t410 - g(1) + t466;
t306 = t406 * rSges(3,1) - rSges(3,2) * t405;
t628 = t305 * t411;
t743 = t306 * t410 - t411 * t628 - g(2) + t520;
t441 = -t398 - t762;
t465 = t519 * t571;
t660 = pkin(1) * qJD(1);
t561 = t416 * t660;
t579 = t326 + t376;
t482 = -t561 + t579;
t57 = t411 * t522 + t465 + t482;
t560 = t418 * t660;
t457 = t537 - t577 + t560;
t552 = -t201 - t595;
t523 = t308 - t552;
t704 = t581 * t572;
t58 = t411 * t523 + t457 - t704;
t659 = t406 * t58;
t739 = ((-t58 * t414 + t441 * t57) * t406 + (-t57 * rSges(6,2) + t441 * t58) * t405) * t411 + (-t57 * t699 + (-t666 - t707) * t659) * qJD(3) - t465 * t58;
t738 = -t763 * t403 + t766 * t404;
t221 = t411 * t229;
t281 = t411 * t307;
t737 = -rSges(4,1) * t539 + t221 + t281 + t358 + t475;
t190 = t411 * t200;
t703 = t281 - t190;
t736 = -rSges(5,1) * t557 + t411 * t218 - t398 * t620 + t376 + t586 + t703;
t735 = t777 + t784;
t233 = -t561 - t628;
t732 = t596 * t411 - t579 + t701 + t703;
t731 = -t755 * t405 + t760 * t406;
t730 = t754 * t405 + t759 * t406;
t729 = t760 * t405 + t755 * t406;
t728 = t759 * t405 - t754 * t406;
t727 = t757 * qJD(3) + t752;
t726 = t758 * qJD(3) + t753;
t725 = t768 * qJD(3) - t152 * t417 - t154 * t415 + t778 * t403 + t780 * t404;
t724 = -t769 * qJD(3) + t151 * t417 + t153 * t415 + t779 * t403 - t781 * t404;
t723 = t761 * t405 - t756 * t406;
t722 = t756 * t405 + t761 * t406;
t575 = t352 + t503;
t576 = -t350 + t353;
t715 = (-t765 * t403 + t764 * t404 - t415 * t575 + t417 * t576) * t411;
t714 = t640 + t783;
t711 = t773 * t411;
t463 = -t301 - t666;
t601 = -t201 - t220;
t679 = t301 * t572 - t411 * (t308 - t601) + t577;
t76 = t560 - t679;
t708 = t463 * t76;
t702 = t304 + t408;
t698 = t741 * t621;
t693 = g(1) * t406 + g(2) * t405;
t691 = t325 + t700;
t690 = t408 + t762;
t188 = t230 + t308;
t683 = -t188 * t411 + t359 * t572;
t592 = -Icges(4,2) * t617 + t227 - t365;
t594 = t352 * t405 + t225;
t676 = -t415 * t592 - t417 * t594;
t675 = t270 / 0.2e1;
t674 = t271 / 0.2e1;
t673 = t405 / 0.2e1;
t672 = -t406 / 0.2e1;
t658 = t411 * t57;
t606 = -t200 * t572 + t201 * t571;
t46 = -t565 + (t405 * t596 + t406 * t595) * qJD(3) + t606;
t657 = t46 * t403;
t546 = t359 * t571;
t464 = -t546 - t561;
t120 = t411 * t590 + t464;
t643 = t120 * t406;
t605 = -t405 * t200 + t406 * t201;
t593 = -t352 * t406 - t226;
t591 = -t350 * t406 + t228;
t588 = -t749 * t624 + t699;
t587 = -t749 * t623 + t698;
t573 = t405 ^ 2 + t406 ^ 2;
t564 = pkin(3) * t615;
t555 = t145 * t571 + t146 * t572 - t270 * t200;
t554 = t405 * t146 + (t145 - t190) * t406;
t534 = -pkin(2) - t662;
t533 = -t572 / 0.2e1;
t532 = t572 / 0.2e1;
t531 = -t571 / 0.2e1;
t530 = t571 / 0.2e1;
t521 = t573 * t666;
t516 = t406 * t463;
t362 = rSges(2,1) * t418 - rSges(2,2) * t416;
t360 = rSges(2,1) * t416 + rSges(2,2) * t418;
t513 = t405 * t58 + t406 * t57;
t121 = t560 - t683;
t500 = -t121 * t405 - t643;
t499 = t141 * t406 + t143 * t405;
t492 = t218 * t405 + t220 * t406;
t489 = t229 * t405 + t230 * t406;
t481 = -pkin(3) * t568 + t589;
t476 = t513 * t404;
t165 = -t218 + t578;
t461 = qJD(3) * t516 + t376;
t458 = -t415 * t591 + t417 * t593;
t187 = t405 * t534 + t396 + t574;
t111 = t524 + t595;
t166 = t220 + t524;
t455 = -rSges(5,3) * t620 + t547 - t549;
t440 = t461 - t561;
t110 = t393 + (-t403 * t741 - t404 * t749) * t405 + t578;
t426 = (t534 * t643 + (t120 * (-rSges(4,3) - pkin(7)) + t121 * t534) * t405) * t411;
t75 = t411 * t553 + t440;
t423 = ((t75 * (-t398 - t375) - t76 * t414) * t411 + qJD(3) * t708) * t406;
t422 = (((t78 - t706 + t714) * t405 + ((t776 + t799) * t406 + t742 + t782 + t798) * t406) * qJD(3) + t753) * t530 + (-t771 * qJD(3) + t312 * t417 + t313 * t415 + t774 * t403 - t775 * t404) * t411 + (Icges(3,3) - t770) * t410 + (-t767 + t716) * t675 + (-t712 + t717) * t674 + (t722 + t724) * t532 + (((t735 * t406 - t714 + t719) * t406 + (t735 * t405 - t193 + t517 + t720 - t785) * t405) * qJD(3) + t727 - t752) * t533 + (t723 - t725 + t726) * t531;
t269 = t359 * t406;
t268 = t359 * t405;
t253 = t301 * t406;
t249 = t301 * t405;
t234 = t306 * t411 + t560;
t144 = t489 * qJD(3);
t71 = qJD(3) * t492 + t606;
t5 = -qJDD(5) * t404 + t596 * t270 + t552 * t271 + (t405 * t611 + t406 * t612 + t566) * qJD(3) + t555;
t1 = [Icges(2,3) * qJDD(1) + t422 + (t743 * (t306 + t409) + t744 * (-t305 - t667) + (-t260 - t560 + t234) * t233) * m(3) + (g(1) * t360 - g(2) * t362 + (t360 ^ 2 + t362 ^ 2) * qJDD(1)) * m(2) + (t57 * (-t457 + t691) + (t482 + t57 + t561 + t732) * t58 + t750 * (t111 + t409) + t751 * (t110 - t667) + t739) * m(6) + (t75 * (t455 - t560) + t423 + (-t561 + t75 - t440 + t736) * t76 + t746 * (t166 + t409) + t748 * (t165 - t667)) * m(5) + (t120 * (t548 - t560) + t426 + t745 * (t188 + t409) + t747 * (t187 - t667) + (t120 - t464 - t561 + t737) * t121) * m(4); t422 + (t523 * t658 + (t579 + t732) * t58 + t750 * t111 + t751 * t110 + (t691 - t704) * t57 + t739) * m(6) + (t423 + (-t461 + t736) * t76 + (t455 - t679) * t75 + t746 * t166 + t748 * t165) * m(5) + (t426 + t745 * t188 + t747 * t187 + (t546 + t737) * t121 + (t548 - t683) * t120) * m(4) + (-t233 * t260 - t234 * t628 + (t233 * t411 + t743) * t306 + (t234 * t411 - t744) * t305) * m(3); t758 * t675 + t757 * t674 + (t722 * t411 + t716 * t410 + t720 * t271 + t719 * t270 + (t730 * t405 + t731 * t406) * qJD(3)) * t673 + (t723 * t411 + t717 * t410 + t721 * t271 + t718 * t270 + (t728 * t405 + t729 * t406) * qJD(3)) * t672 + (-t405 * t767 + t712 * t406) * t410 / 0.2e1 - ((t764 * t403 + t765 * t404 + t415 * t576 + t417 * t575) * t411 + ((t405 * t591 - t406 * t592) * t417 + (t405 * t593 + t406 * t594) * t415 + t763 * t404 + t766 * t403) * qJD(3)) * t411 / 0.2e1 + ((-t411 * t767 + t725) * t406 + (-t712 * t411 + t724) * t405) * t411 / 0.2e1 + t727 * t620 / 0.2e1 + t726 * t616 / 0.2e1 + ((-t710 * t572 + t711) * t405 + ((-t676 * t406 + (t458 + t709) * t405 + t738) * qJD(3) + t715) * t406) * t533 + ((t719 * t411 + t731) * t406 + (t720 * t411 + t730) * t405) * t532 + ((t718 * t411 + t729) * t406 + (t721 * t411 + t728) * t405) * t531 + ((-t709 * t571 - t711) * t406 + ((t458 * t405 + (-t676 + t710) * t406 + t738) * qJD(3) + t715) * t405) * t530 + (-(-t71 * t521 + (-t71 * t253 - t702 * t75) * t406 + (-t71 * t249 - t702 * t76) * t405) * qJD(3) + t23 * t516 + t75 * (-pkin(3) * t538 - t283 * t406) + t76 * (-pkin(3) * t541 - t283 * t405) + (qJD(3) * t499 + t218 * t270 + t271 * t601 + t555) * (t492 + t605) + t71 * (t499 + t554) - g(3) * t702 + (t24 * t405 - t693) * t463 + (-t75 * t249 - t76 * (-t253 - t564) + (t71 * t218 + t708) * t406 + (t75 * t301 + t601 * t71) * t405) * t411) * m(5) + ((t229 * t270 - t230 * t271 + (t155 * t406 + t156 * t405) * qJD(3)) * t489 + t144 * ((t155 + t221) * t406 + (-t230 * t411 + t156) * t405) + t500 * t321 + ((-t121 * t411 - t59) * t406 + (t120 * t411 - t60) * t405) * t359 - (t120 * t268 - t121 * t269) * t411 - (t144 * (-t268 * t405 - t269 * t406) + t500 * t361) * qJD(3) + g(1) * t269 + g(2) * t268 - g(3) * t361) * m(4) + (-(t476 + t657) * qJD(5) - (-t57 * t588 + t58 * (-t564 + t587)) * t411 - (-t46 * t521 + (t46 * t587 - t57 * t690) * t406 + (t46 * t588 - t58 * t690) * t405) * qJD(3) + t5 * t605 + t46 * t554 + (t8 * t519 + t57 * t481 + t5 * t595 + t46 * t612 + (t46 * t596 + t519 * t58) * t411) * t406 + (t9 * t519 + t58 * t481 + t5 * t596 + t46 * t611 + (t46 * t552 + t57 * t581) * t411) * t405 - g(1) * (-t564 + t698) - g(2) * (-pkin(3) * t618 + t699) - g(3) * t690 + t693 * t707) * m(6); (-m(5) - m(6)) * (g(1) * t405 - g(2) * t406) + 0.2e1 * (t672 * t9 + t673 * t8) * m(6) + 0.2e1 * (t23 * t673 + t24 * t672) * m(5); (-(t573 * t657 + t476) * qJD(3) + (qJD(3) * t513 + g(3) - t5) * t404 + (qJD(3) * t46 + (t411 * t58 + t8) * t406 + (t9 - t658) * t405 - (-t405 * t57 + t659) * t411 - t693) * t403) * m(6);];
tau = t1;
