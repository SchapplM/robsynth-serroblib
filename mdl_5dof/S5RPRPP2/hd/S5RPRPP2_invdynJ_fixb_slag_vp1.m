% Calculate vector of inverse dynamics joint torques for
% S5RPRPP2
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPP2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPP2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:11:11
% DurationCPUTime: 26.84s
% Computational Cost: add. (10204->633), mult. (12584->744), div. (0->0), fcn. (9691->6), ass. (0->343)
t711 = Icges(5,6) - Icges(6,6);
t710 = Icges(5,1) + Icges(6,1);
t709 = Icges(5,4) - Icges(6,5);
t708 = Icges(6,2) + Icges(5,3);
t336 = sin(qJ(3));
t327 = Icges(6,4) * t336;
t338 = cos(qJ(3));
t403 = Icges(6,1) * t338 + t327;
t326 = Icges(5,5) * t336;
t404 = Icges(5,1) * t338 + t326;
t707 = t403 + t404;
t562 = Icges(4,4) * t336;
t265 = Icges(4,1) * t338 - t562;
t335 = qJ(1) + pkin(7);
t322 = sin(t335);
t323 = cos(t335);
t153 = Icges(4,5) * t322 + t265 * t323;
t702 = t322 * t709 + t323 * t707;
t692 = t153 + t702;
t529 = t323 * t338;
t706 = (Icges(6,4) + Icges(5,5)) * t529;
t705 = t711 * t322;
t530 = t323 * t336;
t695 = t530 * t708 + t705 + t706;
t253 = Icges(4,5) * t338 - Icges(4,6) * t336;
t141 = Icges(4,3) * t322 + t253 * t323;
t257 = Icges(5,4) * t338 + Icges(5,6) * t336;
t145 = Icges(5,2) * t322 + t257 * t323;
t704 = t141 + t145;
t659 = -t338 * t708 + t326 + t327;
t679 = -Icges(4,2) * t338 - t562 + t659;
t328 = Icges(4,4) * t338;
t402 = -Icges(4,2) * t336 + t328;
t558 = Icges(5,5) * t338;
t251 = Icges(5,3) * t336 + t558;
t560 = Icges(6,4) * t338;
t255 = Icges(6,2) * t336 + t560;
t693 = t251 + t255;
t703 = t402 - t693;
t264 = Icges(4,1) * t336 + t328;
t649 = t336 * t710 + t264 - t558 - t560;
t689 = (-Icges(4,6) + t711) * t338 + (-Icges(4,5) - t709) * t336;
t700 = t265 + t707;
t531 = t322 * t338;
t532 = t322 * t336;
t548 = Icges(4,3) * t323;
t140 = Icges(4,5) * t531 - Icges(4,6) * t532 - t548;
t288 = Icges(4,4) * t532;
t559 = Icges(4,5) * t323;
t152 = Icges(4,1) * t531 - t288 - t559;
t553 = Icges(4,6) * t323;
t146 = Icges(4,4) * t531 - Icges(4,2) * t532 - t553;
t540 = t146 * t336;
t390 = -t152 * t338 + t540;
t142 = Icges(6,6) * t323 + t255 * t322;
t148 = Icges(6,5) * t323 + t322 * t403;
t393 = t142 * t336 + t148 * t338;
t144 = -Icges(5,2) * t323 + t257 * t322;
t541 = t144 * t323;
t249 = Icges(6,5) * t338 + Icges(6,6) * t336;
t136 = Icges(6,3) * t323 + t249 * t322;
t543 = t136 * t323;
t138 = -Icges(5,6) * t323 + t251 * t322;
t150 = -Icges(5,4) * t323 + t322 * t404;
t397 = t138 * t336 + t150 * t338;
t631 = t322 * t397;
t661 = t322 * t393 - t541 + t543 + t631;
t620 = -t140 * t323 - t322 * t390 + t661;
t147 = Icges(4,6) * t322 + t323 * t402;
t120 = t153 * t531;
t439 = t141 * t323 - t120;
t51 = -t147 * t532 - t439;
t137 = -Icges(6,3) * t322 + t249 * t323;
t135 = t323 * t137;
t660 = -t145 * t323 + t531 * t702 + t532 * t695 + t135;
t619 = t51 + t660;
t699 = -t322 * t140 - t142 * t530 + (-t148 - t152) * t529;
t698 = t704 * t322 + t529 * t692 + t695 * t530;
t685 = t336 * t679 + t338 * t649;
t544 = t136 * t322;
t697 = t146 * t530 + t544 + t699;
t696 = t138 + t142;
t694 = -t148 - t150;
t617 = -t137 * t322 - t147 * t530 + t698;
t691 = t703 * qJD(3);
t690 = t700 * qJD(3);
t688 = t249 - t253 - t257;
t687 = t679 * qJD(3);
t686 = t649 * qJD(3);
t684 = -t336 * t649 + t338 * t679;
t613 = t689 * t323;
t614 = t689 * t322;
t651 = t322 * t685 + t613;
t650 = t323 * t685 - t614;
t683 = -t146 + t696;
t682 = -t147 + t695;
t681 = -t152 + t694;
t566 = -rSges(6,3) - qJ(5);
t633 = pkin(4) * t529 + t322 * t566;
t506 = rSges(6,1) * t529 + rSges(6,2) * t530 + t633;
t300 = pkin(3) * t529;
t325 = t336 * qJ(4);
t206 = t323 * t325 + t300;
t223 = pkin(2) * t323 + pkin(6) * t322;
t339 = cos(qJ(1));
t334 = t339 * pkin(1);
t626 = t334 + t223;
t635 = t206 + t626;
t678 = t635 + t506;
t677 = -t687 * t323 + (t322 * t402 - t553 - t696) * qJD(1);
t676 = -t687 * t322 + (t323 * t693 - t147 + t705) * qJD(1);
t675 = -t686 * t323 + (-t265 * t322 + t559 + t694) * qJD(1);
t674 = -qJD(1) * t692 + t322 * t686;
t673 = -qJD(1) * t689 + qJD(3) * t684 - t336 * t691 + t338 * t690;
t539 = t147 * t336;
t672 = t336 * t695 + t338 * t692 - t539;
t671 = t390 - t393 - t397;
t133 = t322 * t144;
t54 = t138 * t530 + t150 * t529 + t133;
t636 = t323 * t54;
t670 = t617 * t322 + t323 * t697 - t636;
t669 = t322 * t619 - t620 * t323;
t668 = t685 * qJD(1) + qJD(3) * t688;
t220 = rSges(3,1) * t322 + rSges(3,2) * t323;
t337 = sin(qJ(1));
t584 = pkin(1) * t337;
t208 = -t220 - t584;
t667 = t650 * qJD(1);
t666 = -t336 * t692 + t338 * t682;
t616 = t336 * t681 + t338 * t683;
t665 = t689 * qJD(3);
t664 = t651 * qJD(1);
t468 = qJD(1) * qJD(3);
t217 = -qJDD(3) * t323 + t322 * t468;
t473 = qJD(4) * t338;
t627 = pkin(3) * t338 + t325;
t210 = qJD(3) * t627 - t473;
t330 = t336 * rSges(6,2);
t628 = rSges(6,1) * t338 + t330;
t502 = -qJD(3) * t628 - t210;
t582 = pkin(4) * t338;
t360 = (-qJD(3) * t582 + t502) * qJD(3);
t324 = qJD(4) * t336;
t452 = t322 * t324;
t471 = qJD(5) * t323;
t372 = -t452 - t471;
t201 = t627 * t322;
t316 = t323 * pkin(6);
t222 = pkin(2) * t322 - t316;
t443 = -t222 - t584;
t430 = -t201 + t443;
t598 = -rSges(6,3) * t323 - t322 * t628;
t507 = pkin(4) * t531 + qJ(5) * t323 - t598;
t373 = t430 - t507;
t271 = pkin(3) * t336 - qJ(4) * t338;
t341 = qJD(1) ^ 2;
t464 = t341 * t334;
t467 = qJD(3) * qJD(4);
t634 = qJDD(4) * t336 + t338 * t467;
t375 = t217 * t271 + t323 * t634 - t464;
t576 = rSges(6,1) * t336;
t272 = -rSges(6,2) * t338 + t576;
t583 = pkin(4) * t336;
t440 = t272 + t583;
t475 = qJD(3) * t336;
t456 = t322 * t475;
t245 = pkin(4) * t456;
t436 = -t245 + t471;
t477 = qJD(3) * t322;
t525 = -t272 * t477 + t436 + (t323 * t628 + t633) * qJD(1);
t213 = t223 * qJD(1);
t474 = qJD(3) * t338;
t478 = qJD(1) * t336;
t177 = t322 * t474 + t323 * t478;
t246 = pkin(3) * t456;
t422 = -t246 + t452;
t80 = qJ(4) * t177 + qJD(1) * t300 + t422;
t632 = -t213 - t80;
t3 = -qJDD(5) * t322 + t440 * t217 + t360 * t323 + t373 * qJDD(1) + (t372 - t525 + t632) * qJD(1) + t375;
t663 = g(1) - t3;
t216 = qJDD(3) * t322 + t323 * t468;
t269 = t323 * t324;
t479 = qJD(1) * t323;
t303 = pkin(6) * t479;
t432 = qJDD(1) * t334 - t341 * t584;
t480 = qJD(1) * t322;
t374 = qJD(1) * (-pkin(2) * t480 + t303) + qJDD(1) * t223 + t432;
t454 = t323 * t475;
t361 = -t338 * t480 - t454;
t458 = t322 * t478;
t453 = t323 * t474;
t241 = qJ(4) * t453;
t500 = t241 + t269;
t79 = pkin(3) * t361 - qJ(4) * t458 + t500;
t351 = qJDD(1) * t206 + t374 + t634 * t322 + (t269 + t79) * qJD(1);
t428 = -t271 - t440;
t243 = rSges(6,2) * t453;
t472 = qJD(5) * t322;
t526 = -rSges(6,1) * t454 + pkin(4) * t361 - qJ(5) * t479 + qJD(1) * t598 + t243 - t472;
t2 = qJDD(5) * t323 + t506 * qJDD(1) + t526 * qJD(1) + t428 * t216 + (-qJD(1) * qJD(5) + t360) * t322 + t351;
t662 = -g(2) + t2;
t658 = (t136 - t144) * qJD(1);
t657 = qJD(3) * t669 + t664;
t656 = qJD(3) * t670 + t667;
t655 = qJD(3) * t671 + t336 * t674 + t338 * t676;
t654 = qJD(3) * t672 + t336 * t675 - t338 * t677;
t653 = -t322 * t668 + t323 * t673;
t652 = t322 * t673 + t323 * t668;
t618 = t54 - t697;
t648 = t541 + t698;
t647 = (-t137 + t704) * qJD(1);
t646 = t323 ^ 2;
t311 = t322 * rSges(5,2);
t166 = rSges(5,1) * t529 + rSges(5,3) * t530 + t311;
t64 = t635 + t166;
t645 = qJD(3) * t666 + t336 * t677 + t338 * t675 + t647;
t644 = -qJD(1) * t140 - qJD(3) * t616 - t336 * t676 + t338 * t674 + t658;
t643 = t649 + t703;
t642 = t679 + t700;
t641 = t323 * t679 + t692;
t640 = -t264 * t323 - t530 * t710 + t682 + t706;
t639 = qJD(1) * t671 + t322 * t665 + t647;
t638 = t658 + t665 * t323 + (-t253 * t322 + t548 - t672) * qJD(1);
t465 = -rSges(6,1) - pkin(3) - pkin(4);
t637 = t338 * t465 - pkin(2);
t310 = t322 * rSges(4,3);
t167 = rSges(4,1) * t529 - rSges(4,2) * t530 + t310;
t112 = t167 + t626;
t419 = rSges(5,1) * t338 + rSges(5,3) * t336;
t490 = -t627 - t419;
t219 = qJD(1) * t222;
t630 = -qJD(1) * t201 - t219;
t221 = rSges(3,1) * t323 - rSges(3,2) * t322;
t209 = t221 + t334;
t435 = t269 - t472;
t476 = qJD(3) * t323;
t625 = t428 * t476 + t435;
t624 = -t325 - t330 + t637;
t612 = t627 + t628 + t582;
t611 = (-t336 * t643 + t338 * t642) * qJD(1);
t610 = -t336 * t641 + t338 * t640;
t609 = t639 * t646 + (t645 * t322 + (-t638 + t644) * t323) * t322;
t608 = t644 * t646 + (t638 * t322 + (-t639 + t645) * t323) * t322;
t607 = t688 * qJD(1);
t606 = -Icges(4,2) * t531 + t322 * t659 - t288 - t681;
t605 = t322 * t649 - t683;
t594 = -pkin(4) * t474 + qJD(3) * t612 + t502;
t593 = t336 * t606 + t338 * t605;
t592 = m(5) / 0.2e1;
t591 = m(6) / 0.2e1;
t590 = -m(5) - m(6);
t589 = t216 / 0.2e1;
t588 = t217 / 0.2e1;
t585 = -rSges(5,1) - pkin(3);
t581 = g(2) * t322;
t578 = rSges(4,1) * t338;
t577 = rSges(5,1) * t336;
t274 = rSges(4,1) * t336 + rSges(4,2) * t338;
t205 = t274 * t323;
t62 = qJD(1) * t112 - t274 * t477;
t572 = t205 * t62;
t489 = rSges(4,2) * t532 + rSges(4,3) * t323;
t164 = rSges(4,1) * t531 - t489;
t431 = -t164 + t443;
t457 = t274 * t476;
t61 = qJD(1) * t431 - t457;
t571 = t322 * t61;
t570 = t323 * t61;
t568 = -rSges(6,2) - qJ(4);
t567 = -rSges(5,3) - qJ(4);
t505 = -t166 - t206;
t504 = t201 * t322 + t206 * t323;
t283 = qJ(4) * t529;
t202 = -pkin(3) * t530 + t283;
t503 = qJD(1) * t202 + t322 * t473;
t501 = -qJD(3) * t419 - t210;
t499 = rSges(5,2) * t479 + rSges(5,3) * t453;
t498 = rSges(4,2) * t458 + rSges(4,3) * t479;
t273 = -rSges(5,3) * t338 + t577;
t491 = -t271 - t273;
t487 = t322 ^ 2 + t646;
t463 = t201 * t479 + t322 * t80 + t323 * t79;
t281 = qJ(4) * t531;
t197 = -pkin(3) * t532 + t281;
t462 = t197 * t477 + t202 * t476 + t324;
t461 = -t206 - t506;
t459 = t323 * t585;
t451 = -pkin(2) - t578;
t447 = -t477 / 0.2e1;
t446 = t477 / 0.2e1;
t445 = -t476 / 0.2e1;
t444 = t476 / 0.2e1;
t441 = t316 - t584;
t438 = -t140 + t539;
t437 = qJD(3) * t501;
t434 = t323 * t465;
t280 = rSges(2,1) * t339 - rSges(2,2) * t337;
t275 = rSges(2,1) * t337 + rSges(2,2) * t339;
t279 = -rSges(4,2) * t336 + t578;
t407 = -t322 * t62 - t570;
t107 = rSges(4,1) * t361 - rSges(4,2) * t453 + t498;
t200 = t274 * t322;
t110 = -qJD(3) * t200 + (t279 * t323 + t310) * qJD(1);
t399 = t107 * t323 + t110 * t322;
t388 = t164 * t322 + t167 * t323;
t313 = t323 * rSges(5,2);
t163 = t322 * t419 - t313;
t381 = -t163 + t430;
t377 = t201 * t477 + t206 * t476 + qJD(2) - t473;
t376 = t476 * t491 + t269;
t359 = -qJDD(4) * t338 + t201 * t216 + t336 * t467 + t476 * t79 + t477 * t80 + qJDD(2);
t355 = t336 * t567 + t338 * t585 - pkin(2);
t295 = rSges(6,2) * t529;
t294 = rSges(5,3) * t529;
t291 = rSges(6,2) * t531;
t290 = rSges(5,3) * t531;
t270 = t323 * t473;
t235 = t279 * qJD(3);
t215 = t271 * t480;
t214 = t271 * t477;
t204 = -rSges(5,1) * t530 + t294;
t203 = -rSges(6,1) * t530 + t295;
t199 = -rSges(5,1) * t532 + t290;
t198 = -rSges(6,1) * t532 + t291;
t178 = t453 - t458;
t176 = t487 * t475;
t109 = -t273 * t477 + (t323 * t419 + t311) * qJD(1);
t106 = rSges(5,1) * t361 - rSges(5,3) * t458 + t499;
t60 = qJD(3) * t388 + qJD(2);
t45 = -t214 + (-qJD(3) * t273 + t324) * t322 + t64 * qJD(1);
t44 = qJD(1) * t381 + t376;
t43 = (t163 * t322 + t166 * t323) * qJD(3) + t377;
t42 = -t214 + (-qJD(3) * t272 + t324) * t322 + t678 * qJD(1) + t436;
t41 = qJD(1) * t373 + t625;
t34 = (t322 * t507 + t323 * t506) * qJD(3) + t377;
t33 = qJD(1) * t107 + qJDD(1) * t167 - t216 * t274 - t235 * t477 + t374;
t32 = -t464 - t235 * t476 + t217 * t274 + (-t110 - t213) * qJD(1) + t431 * qJDD(1);
t25 = qJD(3) * t399 + t164 * t216 - t167 * t217 + qJDD(2);
t18 = qJD(1) * t106 + qJDD(1) * t166 + t216 * t491 + t322 * t437 + t351;
t17 = t217 * t273 + t323 * t437 + t381 * qJDD(1) + (-t109 - t452 + t632) * qJD(1) + t375;
t4 = t163 * t216 + t505 * t217 + (t106 * t323 + t109 * t322) * qJD(3) + t359;
t1 = t507 * t216 + t461 * t217 + (t525 * t322 + t526 * t323) * qJD(3) + t359;
t5 = [-m(2) * (-g(1) * t275 + g(2) * t280) + ((-t220 * t341 - g(2) + t432) * t209 + (-t464 + (-0.2e1 * t221 - t334 + t209) * t341 - g(1)) * t208) * m(3) + ((-t636 + (t51 - t120 + (t141 + t540) * t323 + t699) * t323 + (-t631 + (-t137 - t393) * t322 + t648 + t661) * t322) * qJD(3) + t667) * t444 + (t685 * qJD(3) + t690 * t336 + t691 * t338) * qJD(1) + (-(-t41 + (-t507 - t584) * qJD(1) + t625 + t630) * t42 + t41 * (t245 + t246 + t372) + t42 * (t241 + t243 + t303 + t435) + (t42 * t336 * t434 + t41 * (t338 * t568 + t576) * t322) * qJD(3) + ((-t337 * t42 - t339 * t41) * pkin(1) + (t41 * t624 + t42 * t566) * t323 + (t41 * (-pkin(6) - t566) + t624 * t42) * t322) * qJD(1) + t662 * t678 - t663 * (t566 * t323 + (t336 * t568 + t637) * t322 + t441)) * m(6) + (-t44 * t422 + t45 * (t303 + t499 + t500) + (t45 * t336 * t459 + t44 * (t338 * t567 + t577) * t322) * qJD(3) + ((-t337 * t45 - t339 * t44) * pkin(1) + t44 * t355 * t323 + (t44 * (-rSges(5,2) - pkin(6)) + t45 * (-pkin(2) + t490)) * t322) * qJD(1) - (-t44 + (-t163 - t584) * qJD(1) + t376 + t630) * t45 + (-g(2) + t18) * t64 + (-g(1) + t17) * (t322 * t355 + t313 + t441)) * m(5) + (t62 * (t303 + t498) + (t274 * t571 - t572) * qJD(3) + ((-t337 * t62 - t339 * t61) * pkin(1) + (-pkin(2) - t279) * t570 + (t61 * (-rSges(4,3) - pkin(6)) + t62 * t451) * t322) * qJD(1) - (-t457 - t219 - t61 + (-t164 - t584) * qJD(1)) * t62 + (t33 - g(2)) * t112 + (-g(1) + t32) * (t322 * t451 + t441 + t489)) * m(4) + (m(2) * (t275 ^ 2 + t280 ^ 2) + m(3) * (t208 ^ 2 + t221 * t209) + Icges(2,3) + Icges(3,3) - t684) * qJDD(1) + (t650 - t666) * t589 + (-t616 + t651) * t588 + (((t323 * t438 + t543 + t617 - t648) * t323 + (t322 * t438 - t133 + t135 + t439 + t544 + t618 - t660) * t322) * qJD(3) + t657 - t664) * t447 + (t653 + t654) * t446 + (t652 - t655 + t656) * t445; m(3) * qJDD(2) + m(4) * t25 + m(5) * t4 + m(6) * t1 + (-m(3) - m(4) + t590) * g(3); t670 * t589 + t669 * t588 + (t653 * qJD(1) + t608 * qJD(3) + t650 * qJDD(1) + t617 * t216 + t618 * t217) * t322 / 0.2e1 - (t652 * qJD(1) + t609 * qJD(3) + t651 * qJDD(1) + t619 * t216 + t620 * t217) * t323 / 0.2e1 - (((t322 * t641 - t606 * t323) * t338 + (t322 * t640 + t605 * t323) * t336) * qJD(3) + (t642 * t336 + t643 * t338) * qJD(1)) * qJD(1) / 0.2e1 + (t655 * t323 + t654 * t322 + (-t322 * t616 - t323 * t666) * qJD(1)) * qJD(1) / 0.2e1 + (-t322 * t666 + t616 * t323) * qJDD(1) / 0.2e1 + t657 * t480 / 0.2e1 + t656 * t479 / 0.2e1 + ((t613 * t477 - t607) * t322 + ((t593 * t323 + (t610 - t614) * t322) * qJD(3) + t611) * t323) * t447 + ((t322 * t618 + t323 * t617) * qJD(1) + t608) * t446 + ((t322 * t620 + t619 * t323) * qJD(1) + t609) * t445 + ((t614 * t476 + t607) * t323 + ((t610 * t322 + (t593 - t613) * t323) * qJD(3) + t611) * t322) * t444 + (t1 * t504 + (t1 * t507 + t2 * t428) * t322 + (t1 * t506 + t3 * t428) * t323 - g(1) * (t283 + t295) - g(2) * (t281 + t291) - g(3) * t612 - (g(1) * t434 + t465 * t581) * t336 + (-t503 + t594 * t322 + (pkin(4) * t530 + t323 * t428 - t203) * qJD(1)) * t42 + (t215 - t270 + t594 * t323 + (t272 * t322 + t197 + t198) * qJD(1)) * t41 + (t463 + (qJD(1) * t461 + t525) * t322 + (qJD(1) * t507 + t526) * t323 - t462 - (t198 * t322 + t203 * t323 - t487 * t583) * qJD(3)) * t34) * m(6) + (-t44 * (t270 + (-t197 - t199) * qJD(1)) - t45 * (qJD(1) * t204 + t503) - t43 * t462 - ((t204 * t43 + t44 * t490) * t323 + (t199 * t43 + t45 * t490) * t322) * qJD(3) + t44 * t215 + t4 * t504 + t43 * t463 + (t17 * t491 + t44 * t501 + t4 * t166 + t43 * t106 + (t163 * t43 + t45 * t491) * qJD(1)) * t323 + (t18 * t491 + t45 * t501 + t4 * t163 + t43 * t109 + (t273 * t44 + t43 * t505) * qJD(1)) * t322 - g(1) * (t283 + t294) - g(2) * (t281 + t290) + g(3) * t490 - (g(1) * t459 + t581 * t585) * t336) * m(5) + (t25 * t388 + t60 * ((t164 * t323 - t167 * t322) * qJD(1) + t399) + t407 * t235 + (-t32 * t323 - t33 * t322 + (-t323 * t62 + t571) * qJD(1)) * t274 - (t200 * t61 - t572) * qJD(1) - (t60 * (-t200 * t322 - t205 * t323) + t407 * t279) * qJD(3) + g(1) * t205 + g(2) * t200 - g(3) * t279) * m(4); t590 * (-g(3) * t338 + (g(1) * t323 + t581) * t336) - m(5) * (t176 * t43 + t177 * t45 + t178 * t44) - m(6) * (t176 * t34 + t177 * t42 + t178 * t41) + 0.2e1 * ((t44 * t476 + t45 * t477 - t4) * t592 + (t41 * t476 + t42 * t477 - t1) * t591) * t338 + 0.2e1 * ((qJD(3) * t43 + t17 * t323 + t18 * t322 - t44 * t480 + t45 * t479) * t592 + (qJD(3) * t34 + t2 * t322 + t3 * t323 - t41 * t480 + t42 * t479) * t591) * t336; (t322 * t663 + t323 * t662) * m(6);];
tau = t5;
