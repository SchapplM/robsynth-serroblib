% Calculate vector of inverse dynamics joint torques for
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRP9_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP9_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:48:36
% EndTime: 2019-12-31 18:49:10
% DurationCPUTime: 27.08s
% Computational Cost: add. (19247->783), mult. (18503->968), div. (0->0), fcn. (14536->8), ass. (0->425)
t416 = pkin(8) + qJ(3);
t391 = qJ(4) + t416;
t375 = sin(t391);
t376 = cos(t391);
t653 = Icges(6,5) * t376;
t296 = Icges(6,1) * t375 - t653;
t364 = Icges(5,4) * t376;
t298 = Icges(5,1) * t375 + t364;
t789 = t296 + t298;
t657 = Icges(5,4) * t375;
t294 = Icges(5,2) * t376 + t657;
t363 = Icges(6,5) * t375;
t482 = Icges(6,3) * t376 - t363;
t788 = t482 + t294;
t289 = Icges(6,3) * t375 + t653;
t423 = sin(qJ(1));
t424 = cos(qJ(1));
t207 = -Icges(6,6) * t424 + t289 * t423;
t628 = t376 * t423;
t630 = t375 * t423;
t648 = Icges(5,6) * t424;
t213 = Icges(5,4) * t628 - Icges(5,2) * t630 - t648;
t792 = -t207 + t213;
t627 = t376 * t424;
t335 = Icges(6,5) * t627;
t629 = t375 * t424;
t647 = Icges(6,6) * t423;
t208 = Icges(6,3) * t629 + t335 + t647;
t483 = -Icges(5,2) * t375 + t364;
t214 = Icges(5,6) * t423 + t424 * t483;
t791 = -t208 + t214;
t485 = Icges(6,1) * t376 + t363;
t215 = -Icges(6,4) * t424 + t423 * t485;
t336 = Icges(5,4) * t630;
t654 = Icges(5,5) * t424;
t217 = Icges(5,1) * t628 - t336 - t654;
t790 = t215 + t217;
t216 = Icges(6,4) * t423 + t424 * t485;
t299 = Icges(5,1) * t376 - t657;
t218 = Icges(5,5) * t423 + t299 * t424;
t778 = t216 + t218;
t787 = -t482 * t423 + t790;
t786 = (Icges(5,6) - Icges(6,6)) * t376 + (Icges(6,4) + Icges(5,5)) * t375;
t785 = t299 + t485 - t788;
t784 = t289 - t483 - t789;
t783 = t788 * t424 - t778;
t782 = -t298 * t424 - t791;
t781 = t789 * t423 + t792;
t291 = Icges(5,5) * t376 - Icges(5,6) * t375;
t210 = Icges(5,3) * t423 + t291 * t424;
t293 = Icges(6,4) * t376 + Icges(6,6) * t375;
t212 = Icges(6,2) * t423 + t293 * t424;
t780 = -t210 - t212;
t779 = -t788 * t375 + t789 * t376;
t777 = t291 + t293;
t478 = t213 * t375 - t217 * t376;
t480 = t207 * t375 + t215 * t376;
t776 = -t478 + t480;
t417 = qJD(3) + qJD(4);
t775 = t785 * t417;
t774 = t784 * t417;
t773 = -t778 * qJD(1) + t781 * t417;
t347 = t417 * t424;
t772 = -t296 * t347 + t782 * t417 + (-t299 * t423 - t215 + t654) * qJD(1);
t346 = t417 * t423;
t771 = -t294 * t346 + t787 * t417 + (-t289 * t424 + t214 - t647) * qJD(1);
t770 = -t783 * t417 + (-t423 * t483 + t207 + t648) * qJD(1);
t769 = t786 * t424;
t768 = t786 * t423;
t738 = t423 * t779 - t769;
t737 = t424 * t779 + t768;
t765 = rSges(6,3) + qJ(5);
t711 = t376 * rSges(6,1) + t375 * rSges(6,3);
t767 = -t375 * qJ(5) - t711;
t766 = t780 * qJD(1);
t681 = rSges(6,1) + pkin(4);
t610 = -t208 * t630 - t216 * t628;
t172 = t218 * t628;
t516 = t424 * t210 - t172;
t78 = -t214 * t630 - t516;
t741 = -t212 * t424 - t610 + t78;
t80 = t208 * t629 + t423 * t212 + t216 * t627;
t608 = t423 * t210 + t218 * t627;
t82 = -t214 * t629 + t608;
t739 = t80 + t82;
t764 = t786 * qJD(1) + t774 * t375 + t775 * t376;
t763 = -t770 * t375 + t772 * t376 - t766;
t644 = Icges(5,3) * t424;
t209 = Icges(5,5) * t628 - Icges(5,6) * t630 - t644;
t211 = -Icges(6,2) * t424 + t293 * t423;
t706 = qJD(1) * t211;
t762 = -qJD(1) * t209 + t771 * t375 + t773 * t376 - t706;
t761 = t781 * t347 + (-Icges(6,1) * t629 + t335 + t782) * t346 + t785 * qJD(1);
t760 = (-Icges(5,2) * t628 - t336 + t787) * t347 + t783 * t346 + t784 * qJD(1);
t759 = qJD(1) * t779 - t777 * t417;
t758 = t776 * qJD(1) + t768 * t417 + t766;
t640 = t214 * t375;
t757 = t706 + t769 * t417 + (t208 * t375 + t291 * t423 + t778 * t376 - t640 - t644) * qJD(1);
t756 = t376 * pkin(4) - t767;
t192 = t423 * t211;
t79 = t207 * t629 + t215 * t627 + t192;
t755 = t737 * qJD(1) - t347 * t79;
t332 = qJ(5) * t627;
t550 = t376 * t347;
t565 = qJD(1) * t424;
t754 = rSges(6,2) * t565 + rSges(6,3) * t550 + t417 * t332;
t580 = t765 * t628;
t579 = rSges(6,3) * t627 + t332;
t75 = -t424 * t211 + t423 * t480;
t753 = -t738 * qJD(1) + t347 * t75;
t752 = t758 * t423 + t762 * t424;
t751 = -t757 * t423 + t763 * t424;
t750 = t762 * t423 - t758 * t424;
t749 = t763 * t423 + t757 * t424;
t460 = t478 * t423;
t620 = t424 * t209;
t77 = -t460 - t620;
t748 = t741 * t346 - t347 * t77 - t753;
t609 = -t423 * t209 - t217 * t627;
t81 = -t213 * t629 - t609;
t747 = t739 * t346 - t347 * t81 + t755;
t746 = -t759 * t423 + t764 * t424;
t745 = t764 * t423 + t759 * t424;
t744 = t773 * t375 - t771 * t376;
t743 = t772 * t375 + t770 * t376;
t742 = t75 + t77;
t740 = t79 + t81;
t387 = sin(t416);
t736 = rSges(4,2) * t387;
t557 = qJD(1) * qJD(3);
t378 = t423 * t557;
t556 = qJD(1) * qJD(4);
t248 = t423 * t556 + t378 + (-qJDD(3) - qJDD(4)) * t424;
t388 = cos(t416);
t622 = t388 * qJD(3) ^ 2;
t554 = pkin(3) * t622;
t562 = qJD(5) * t417;
t449 = qJDD(5) * t375 + t376 * t562 - t554;
t421 = cos(pkin(8));
t377 = t421 * pkin(2) + pkin(1);
t374 = pkin(3) * t388;
t318 = t374 + t377;
t422 = -pkin(6) - qJ(2);
t415 = -pkin(7) + t422;
t587 = -t423 * t318 - t424 * t415;
t621 = t422 * t424;
t166 = t377 * t423 + t587 + t621;
t396 = t424 * qJ(2);
t675 = pkin(1) - t377;
t463 = t423 * t675 - t621;
t224 = -t396 + t463;
t348 = pkin(1) * t423 - t396;
t597 = t224 - t348;
t546 = t166 + t597;
t413 = t424 * rSges(6,2);
t598 = t756 * t423 - t413;
t488 = t546 - t598;
t561 = qJD(5) * t423;
t535 = t375 * t561;
t323 = -qJDD(3) * t424 + t378;
t558 = qJD(1) * qJD(2);
t573 = qJDD(2) * t423 + t424 * t558;
t680 = pkin(3) * t387;
t543 = t323 * t680 + t573;
t626 = t387 * t423;
t555 = pkin(3) * t626;
t337 = qJD(3) * t555;
t566 = qJD(1) * t423;
t360 = t415 * t566;
t369 = t422 * t566;
t583 = t318 - t377;
t141 = t565 * t583 - t337 - t360 + t369;
t394 = t423 * qJ(2);
t350 = t424 * pkin(1) + t394;
t393 = qJD(2) * t424;
t286 = qJD(1) * t350 - t393;
t604 = t369 - (-t424 * t675 - t394) * qJD(1) - t286;
t547 = -t141 + t604;
t301 = rSges(6,1) * t375 - rSges(6,3) * t376;
t589 = pkin(4) * t375 - qJ(5) * t376 + t301;
t526 = pkin(4) * t417 - qJD(5);
t605 = -t376 * t526 + t767 * t417;
t409 = t423 * rSges(6,2);
t614 = (pkin(4) * t565 + qJ(5) * t346) * t376 + (qJ(5) * t565 - t423 * t526) * t375 - t301 * t346 + (t424 * t711 + t409) * qJD(1);
t10 = t605 * t347 + t589 * t248 + t449 * t424 + t488 * qJDD(1) + (-t535 + t547 - t614) * qJD(1) + t543;
t578 = t337 + t393;
t466 = t535 - t578;
t307 = t424 * t318;
t355 = t424 * t377;
t571 = -t415 + t422;
t167 = t423 * t571 + t307 - t355;
t507 = -t422 * t423 + t355;
t225 = t507 - t350;
t596 = t225 + t350;
t545 = t167 + t596;
t593 = t681 * t627 + t765 * t629 + t409;
t59 = -t589 * t346 + (t545 + t593) * qJD(1) + t466;
t735 = t424 * (qJD(1) * t59 + t10);
t734 = t790 * t375 + t792 * t376;
t733 = t778 * t375 + t791 * t376;
t560 = qJD(5) * t424;
t327 = t375 * t560;
t455 = -t347 * t375 - t376 * t566;
t540 = t375 * t566;
t615 = t681 * t455 - t765 * t540 + t327 + t754;
t726 = -t681 * t630 + t580;
t725 = -t681 * t629 + t579;
t720 = t760 * t375 + t761 * t376;
t563 = qJD(3) * t424;
t564 = qJD(3) * t423;
t612 = -t166 * t564 + t167 * t563;
t49 = -qJD(5) * t376 + t346 * t598 + t347 * t593 + t612;
t541 = t424 * t681;
t678 = g(2) * t423;
t719 = -t49 * (qJD(5) * t375 + t726 * t346 + t725 * t347) - t59 * (t725 * qJD(1) - t346 * t756 + t376 * t561) - (-g(1) * t541 - t678 * t681) * t375;
t718 = t777 * qJD(1) - t769 * t346 + t768 * t347;
t717 = t726 * qJD(1) + t347 * t756 - t376 * t560 + t589 * t566;
t302 = rSges(5,1) * t375 + rSges(5,2) * t376;
t263 = t302 * t423;
t267 = t302 * t424;
t221 = rSges(5,1) * t628 - rSges(5,2) * t630 - t424 * rSges(5,3);
t406 = t423 * rSges(5,3);
t223 = rSges(5,1) * t627 - rSges(5,2) * t629 + t406;
t68 = t221 * t346 + t223 * t347 + t612;
t392 = qJD(2) * t423;
t537 = t387 * t563;
t505 = pkin(3) * t537;
t469 = t392 - t505;
t451 = -t302 * t347 + t469;
t504 = -t221 + t546;
t71 = qJD(1) * t504 + t451;
t367 = t376 * rSges(5,1);
t710 = -rSges(5,2) * t375 + t367;
t72 = -t302 * t346 + (t223 + t545) * qJD(1) - t578;
t716 = -(qJD(1) * t263 - t347 * t710) * t71 - t68 * (-t346 * t263 - t267 * t347) - t72 * (-qJD(1) * t267 - t346 * t710);
t322 = qJDD(3) * t423 + t424 * t557;
t247 = qJDD(4) * t423 + t424 * t556 + t322;
t140 = -t505 + (-t423 * t583 + t424 * t571) * qJD(1);
t382 = qJ(2) * t565;
t572 = t382 + t392;
t465 = -qJDD(2) * t424 + qJD(1) * (-pkin(1) * t566 + t572) + qJDD(1) * t350 + t423 * t558;
t452 = qJD(1) * (qJD(1) * t463 - t382) + qJDD(1) * t225 + t465;
t448 = qJD(1) * t140 + qJDD(1) * t167 + t452;
t631 = t322 * t387;
t11 = -pkin(3) * t631 + t605 * t346 - t589 * t247 + t593 * qJDD(1) + t449 * t423 + (t327 + t615) * qJD(1) + t448;
t715 = t11 * t423;
t325 = qJD(1) * t348;
t712 = qJD(1) * t224 - t325;
t373 = Icges(4,4) * t388;
t484 = -Icges(4,2) * t387 + t373;
t312 = Icges(4,1) * t387 + t373;
t707 = g(1) * t424 + t678;
t705 = qJD(1) * t166 + t712;
t420 = sin(pkin(8));
t672 = rSges(3,2) * t420;
t674 = rSges(3,1) * t421;
t272 = t423 * rSges(3,3) + (-t672 + t674) * t424;
t624 = t388 * t423;
t645 = Icges(4,3) * t424;
t228 = Icges(4,5) * t624 - Icges(4,6) * t626 - t645;
t354 = Icges(4,4) * t626;
t655 = Icges(4,5) * t424;
t232 = Icges(4,1) * t624 - t354 - t655;
t649 = Icges(4,6) * t424;
t230 = Icges(4,4) * t624 - Icges(4,2) * t626 - t649;
t639 = t230 * t387;
t476 = -t232 * t388 + t639;
t85 = -t424 * t228 - t423 * t476;
t309 = Icges(4,5) * t388 - Icges(4,6) * t387;
t308 = Icges(4,5) * t387 + Icges(4,6) * t388;
t457 = qJD(3) * t308;
t658 = Icges(4,4) * t387;
t313 = Icges(4,1) * t388 - t658;
t233 = Icges(4,5) * t423 + t313 * t424;
t231 = Icges(4,6) * t423 + t424 * t484;
t638 = t231 * t387;
t475 = -t233 * t388 + t638;
t700 = -t424 * t457 + (-t309 * t423 + t475 + t645) * qJD(1);
t229 = Icges(4,3) * t423 + t309 * t424;
t568 = qJD(1) * t229;
t699 = qJD(1) * t476 - t423 * t457 + t568;
t310 = Icges(4,2) * t388 + t658;
t470 = t310 * t387 - t312 * t388;
t696 = t470 * qJD(1) + t309 * qJD(3);
t695 = t423 * (-t310 * t424 + t233) - t424 * (-Icges(4,2) * t624 + t232 - t354);
t691 = t247 / 0.2e1;
t690 = t248 / 0.2e1;
t689 = t322 / 0.2e1;
t688 = t323 / 0.2e1;
t687 = -t346 / 0.2e1;
t686 = t346 / 0.2e1;
t685 = -t347 / 0.2e1;
t684 = t347 / 0.2e1;
t683 = t423 / 0.2e1;
t682 = -t424 / 0.2e1;
t677 = -qJD(1) / 0.2e1;
t676 = qJD(1) / 0.2e1;
t673 = rSges(4,1) * t388;
t671 = rSges(4,2) * t388;
t315 = rSges(4,1) * t387 + t671;
t280 = t315 * t424;
t407 = t423 * rSges(4,3);
t623 = t388 * t424;
t625 = t387 * t424;
t244 = rSges(4,1) * t623 - rSges(4,2) * t625 + t407;
t92 = -t315 * t564 - t393 + (t244 + t596) * qJD(1);
t669 = t280 * t92;
t499 = -t315 * t563 + t392;
t576 = rSges(4,2) * t626 + t424 * rSges(4,3);
t243 = rSges(4,1) * t624 - t576;
t544 = -t243 + t597;
t91 = qJD(1) * t544 + t499;
t665 = t423 * t91;
t664 = t71 * t302;
t663 = qJDD(1) / 0.2e1;
t642 = qJD(1) * t49;
t633 = t308 * t423;
t632 = t308 * t424;
t611 = -t423 * t166 + t424 * t167;
t607 = -t423 * t228 - t232 * t623;
t606 = t423 * t229 + t233 * t623;
t603 = t423 * t221 + t424 * t223;
t199 = t272 + t350;
t586 = -t310 + t313;
t585 = t312 + t484;
t582 = rSges(5,2) * t540 + rSges(5,3) * t565;
t581 = rSges(4,3) * t565 + t566 * t736;
t553 = t423 * t674;
t370 = t423 * t672;
t574 = t424 * rSges(3,3) + t370;
t271 = t553 - t574;
t577 = -t348 - t271;
t575 = rSges(3,3) * t565 + qJD(1) * t370;
t567 = qJD(1) * t309;
t108 = -t423 * t470 - t632;
t559 = t108 * qJD(1);
t551 = qJD(3) * t374;
t135 = rSges(5,1) * t455 - rSges(5,2) * t550 + t582;
t137 = -t417 * t263 + (t424 * t710 + t406) * qJD(1);
t549 = t424 * t135 + t423 * t137 + t221 * t565;
t548 = t424 * t140 + t423 * t141 - t166 * t565;
t538 = t387 * t565;
t534 = -pkin(1) - t674;
t533 = t566 / 0.2e1;
t532 = t565 / 0.2e1;
t531 = -t564 / 0.2e1;
t530 = t564 / 0.2e1;
t529 = -t563 / 0.2e1;
t528 = t563 / 0.2e1;
t456 = -t302 - t680;
t527 = -t377 - t673;
t525 = (-t423 ^ 2 - t424 ^ 2) * t387;
t183 = t233 * t624;
t515 = t424 * t229 - t183;
t514 = -t209 + t640;
t513 = -t228 + t638;
t508 = -t415 * t423 + t307;
t503 = t598 * t423 + t593 * t424;
t501 = -t589 - t680;
t246 = t710 * t417;
t500 = -t246 - t551;
t351 = rSges(2,1) * t424 - rSges(2,2) * t423;
t349 = rSges(2,1) * t423 + rSges(2,2) * t424;
t316 = t673 - t736;
t117 = t231 * t388 + t233 * t387;
t458 = qJD(3) * t310;
t145 = -t424 * t458 + (-t423 * t484 + t649) * qJD(1);
t459 = qJD(3) * t312;
t147 = -t424 * t459 + (-t313 * t423 + t655) * qJD(1);
t433 = -qJD(3) * t117 - t145 * t387 + t147 * t388 + t568;
t116 = t230 * t388 + t232 * t387;
t146 = qJD(1) * t231 - t423 * t458;
t148 = qJD(1) * t233 - t423 * t459;
t434 = qJD(1) * t228 - qJD(3) * t116 - t146 * t387 + t148 * t388;
t495 = -(t423 * t699 + t434 * t424) * t424 + (t423 * t700 + t433 * t424) * t423;
t494 = -(t434 * t423 - t424 * t699) * t424 + (t433 * t423 - t424 * t700) * t423;
t462 = t327 + t469;
t445 = -t347 * t589 + t462;
t58 = qJD(1) * t488 + t445;
t493 = t423 * t59 + t424 * t58;
t492 = -t423 * t72 - t424 * t71;
t86 = -t231 * t626 - t515;
t491 = t423 * t86 - t424 * t85;
t87 = -t230 * t625 - t607;
t88 = -t231 * t625 + t606;
t490 = t423 * t88 - t424 * t87;
t489 = -t423 * t92 - t424 * t91;
t149 = -t563 * t671 + (-t388 * t566 - t537) * rSges(4,1) + t581;
t279 = t315 * t423;
t150 = -qJD(3) * t279 + (t316 * t424 + t407) * qJD(1);
t481 = t149 * t424 + t150 * t423;
t474 = t243 * t423 + t244 * t424;
t471 = t310 * t388 + t312 * t387;
t468 = -t551 + t605;
t467 = t140 * t563 + t141 * t564 - t322 * t166 - t167 * t323;
t464 = t614 * t423 + t615 * t424 + t598 * t565;
t450 = t230 * t424 - t231 * t423;
t447 = (-t387 * t585 + t388 * t586) * qJD(1);
t444 = -t318 - t756;
t284 = t484 * qJD(3);
t285 = t313 * qJD(3);
t432 = qJD(1) * t308 - qJD(3) * t471 - t284 * t387 + t285 * t388;
t430 = -t387 * t695 + t450 * t388;
t429 = (t739 * t423 - t740 * t424) * t691 + (t741 * t423 - t742 * t424) * t690 + (t718 * t423 + t720 * t424) * t687 + (t752 * t424 + t751 * t423 + (t740 * t423 + t739 * t424) * qJD(1)) * t686 + (t750 * t424 + t749 * t423 + (t742 * t423 + t741 * t424) * qJD(1)) * t685 + (t720 * t423 - t718 * t424) * t684 + (t746 * qJD(1) + t737 * qJDD(1) + t739 * t247 + t740 * t248 + t751 * t346 + t752 * t347) * t683 + (t745 * qJD(1) + t738 * qJDD(1) + t741 * t247 + t742 * t248 + t749 * t346 + t750 * t347) * t682 + (t761 * t375 - t760 * t376) * t677 + (t744 * t424 + t743 * t423 + (t734 * t423 + t733 * t424) * qJD(1)) * t676 + (t733 * t423 - t734 * t424) * t663 + t748 * t533 + t747 * t532;
t287 = t316 * qJD(3);
t169 = qJD(1) * t199 - t393;
t168 = qJD(1) * t577 + t392;
t138 = t474 * qJD(3);
t109 = -t424 * t470 + t633;
t107 = t109 * qJD(1);
t90 = qJDD(1) * t272 + qJD(1) * (-qJD(1) * t553 + t575) + t465;
t89 = t577 * qJDD(1) + (-qJD(1) * t272 - t286) * qJD(1) + t573;
t67 = -qJD(3) * t475 + t145 * t388 + t147 * t387;
t66 = -t476 * qJD(3) + t146 * t388 + t148 * t387;
t65 = t432 * t423 - t424 * t696;
t64 = t423 * t696 + t432 * t424;
t55 = qJD(1) * t149 + qJDD(1) * t244 - t287 * t564 - t315 * t322 + t452;
t54 = -t287 * t563 + t315 * t323 + t544 * qJDD(1) + (-t150 + t604) * qJD(1) + t573;
t44 = qJD(3) * t490 + t107;
t43 = qJD(3) * t491 + t559;
t20 = qJD(1) * t135 + qJDD(1) * t223 - t246 * t346 - t247 * t302 + (-t423 * t622 - t631) * pkin(3) + t448;
t19 = -t424 * t554 - t246 * t347 + t248 * t302 + t504 * qJDD(1) + (-t137 + t547) * qJD(1) + t543;
t18 = t135 * t347 + t137 * t346 + t221 * t247 - t223 * t248 + t467;
t9 = -qJDD(5) * t376 + t247 * t598 - t248 * t593 + t346 * t614 + t347 * t615 + t375 * t562 + t467;
t1 = [(t107 + ((t86 - t183 + (t229 + t639) * t424 + t607) * t424 + t606 * t423) * qJD(3)) * t528 - m(2) * (-g(1) * t349 + g(2) * t351) + (t109 + t117) * t689 + (t108 + t116) * t688 + ((t78 + (t213 * t424 + t214 * t423) * t375 + t516 + t609) * t347 + (-t217 * t628 + t620 + t77 + (t213 * t423 - t214 * t424) * t375 + t608 + t80) * t346 + t755) * t684 + (t43 - t559 + ((t424 * t513 - t606 + t88) * t424 + (t423 * t513 + t515 + t87) * t423) * qJD(3)) * t531 + (t64 + t67) * t530 + (t66 + t65 + t44) * t529 + (-qJD(3) * t470 + t284 * t388 + t285 * t387 + t775 * t375 - t774 * t376) * qJD(1) + (-(-qJD(1) * t598 + t445 - t58 + t705) * t59 + t58 * (t360 - t466) + t59 * (t462 + t754) + (-t59 * t375 * t541 + (t375 * t681 - t376 * t765) * t423 * t58) * t417 + ((-t59 * t415 + t444 * t58) * t424 + (-t58 * rSges(6,2) + t444 * t59) * t423) * qJD(1) + (t10 - g(1)) * (t413 + (-t375 * t765 - t376 * t681) * t423 + t587) + (t11 - g(2)) * (t508 + t593)) * m(6) + (-(-qJD(1) * t221 + t451 + t705 - t71) * t72 + t71 * (t360 + t578) + t72 * (t469 + t582) + (-t267 * t72 + t423 * t664) * t417 + ((-t71 * rSges(5,3) + t72 * (-t318 - t367)) * t423 + (t71 * (-t318 - t710) - t72 * t415) * t424) * qJD(1) + (t20 - g(2)) * (t223 + t508) + (t19 - g(1)) * (-t221 + t587)) * m(5) + (-(-qJD(1) * t243 + t499 + t712 - t91) * t92 + t91 * (t369 + t393) + t92 * (t392 + t581) + (t315 * t665 - t669) * qJD(3) + ((-t91 * rSges(4,3) + t527 * t92) * t423 + (t91 * (-t316 - t377) - t92 * t422) * t424) * qJD(1) + (t55 - g(2)) * (t244 + t507) + (t54 - g(1)) * (t527 * t423 + t576 - t621)) * m(4) + (t168 * t393 + t169 * (t572 + t575) + (t168 * (t534 + t672) * t424 + (t168 * (-rSges(3,3) - qJ(2)) + t169 * t534) * t423) * qJD(1) - (-qJD(1) * t271 - t168 - t325 + t392) * t169 + (t90 - g(2)) * t199 + (t89 - g(1)) * (t534 * t423 + t396 + t574)) * m(3) + (t733 + t737) * t691 + (t734 + t738) * t690 + ((t514 * t424 - t460 - t608 + t82) * t347 + (t423 * t514 - t172 - t192 + t610 + (-t776 - t780) * t424 + t740) * t346 + t748 + t753) * t687 + (t743 + t746) * t686 + (-t744 + t745 + t747) * t685 + (m(2) * (t349 ^ 2 + t351 ^ 2) + t471 + Icges(2,3) + Icges(3,2) * t421 ^ 2 + (Icges(3,1) * t420 + 0.2e1 * Icges(3,4) * t421) * t420 + t788 * t376 + t789 * t375) * qJDD(1); (-m(3) - m(4) - m(5) - m(6)) * (g(1) * t423 - g(2) * t424) + 0.2e1 * (t10 * t683 + t11 * t682) * m(6) + 0.2e1 * (t19 * t683 + t20 * t682) * m(5) + 0.2e1 * (t54 * t683 + t55 * t682) * m(4) + 0.2e1 * (t682 * t90 + t683 * t89) * m(3); ((-t564 * t632 + t567) * t423 + (t447 + (t423 * t633 + t430) * qJD(3)) * t424) * t531 + (qJD(1) * t64 + qJD(3) * t495 + qJDD(1) * t109 + t322 * t88 + t323 * t87) * t683 + (qJD(1) * t65 + qJD(3) * t494 + qJDD(1) * t108 + t322 * t86 + t323 * t85) * t682 + ((t87 * t423 + t88 * t424) * qJD(1) + t495) * t530 + ((t85 * t423 + t86 * t424) * qJD(1) + t494) * t529 + ((-t563 * t633 - t567) * t424 + (t447 + (t424 * t632 + t430) * qJD(3)) * t423) * t528 + t43 * t533 + t429 + ((t387 * t586 + t388 * t585) * qJD(1) + (t450 * t387 + t388 * t695) * qJD(3)) * t677 + t44 * t532 + (t423 * t67 - t424 * t66 + (t116 * t423 + t117 * t424) * qJD(1)) * t676 + t491 * t688 + t490 * t689 + (-t116 * t424 + t117 * t423) * t663 + (-g(1) * (-pkin(3) * t625 + t579) - g(2) * (-t555 + t580) - g(3) * (t374 + t756) - (-t59 * t538 + (-t388 * t493 + t49 * t525) * qJD(3)) * pkin(3) + t9 * (t503 + t611) + t49 * (t464 + t548) + t501 * t735 + (t11 * t501 + t59 * t468 + (-t167 - t593) * t642) * t423 + (t424 * t468 + t717) * t58 + t719) * m(6) + (t18 * (t603 + t611) + t68 * (t548 + t549) + (t500 * t71 + (qJD(1) * t72 + t19) * t456) * t424 + (t20 * t456 + t72 * t500 + (t664 + t68 * (-t167 - t223)) * qJD(1)) * t423 - (-t72 * t538 + (t388 * t492 + t525 * t68) * qJD(3)) * pkin(3) - g(3) * (t710 + t374) - t707 * t456 + t716) * m(5) + (-(t279 * t91 - t669) * qJD(1) - (t138 * (-t279 * t423 - t280 * t424) + t489 * t316) * qJD(3) + (qJD(3) * t481 + t243 * t322 - t244 * t323) * t474 + t138 * ((t243 * t424 - t244 * t423) * qJD(1) + t481) + t489 * t287 + (-t55 * t423 - t54 * t424 + (-t424 * t92 + t665) * qJD(1)) * t315 + g(1) * t280 + g(2) * t279 - g(3) * t316) * m(4); t429 + (t9 * t503 + t49 * t464 + (t59 * t605 - t593 * t642) * t423 - g(1) * t579 - g(2) * t580 - g(3) * t756 + (-t715 - t735) * t589 + (t605 * t424 + t717) * t58 + t719) * m(6) + (t18 * t603 + t68 * (-t223 * t566 + t549) + t492 * t246 + (-t19 * t424 - t20 * t423 + (t423 * t71 - t424 * t72) * qJD(1)) * t302 + g(1) * t267 + g(2) * t263 - g(3) * t710 + t716) * m(5); ((-t346 * t59 - t347 * t58 + t417 * t493 + g(3) - t9) * t376 + (t10 * t424 + t715 + (-t346 * t423 - t347 * t424 + t417) * t49 - t707) * t375) * m(6);];
tau = t1;
