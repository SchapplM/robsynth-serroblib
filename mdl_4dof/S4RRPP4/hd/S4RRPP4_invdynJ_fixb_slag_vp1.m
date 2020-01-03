% Calculate vector of inverse dynamics joint torques for
% S4RRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:59:17
% DurationCPUTime: 21.64s
% Computational Cost: add. (4662->583), mult. (11944->702), div. (0->0), fcn. (9343->4), ass. (0->318)
t676 = Icges(4,6) - Icges(5,6);
t675 = Icges(4,1) + Icges(5,1);
t674 = Icges(4,4) - Icges(5,5);
t673 = Icges(5,2) + Icges(4,3);
t323 = sin(qJ(2));
t306 = Icges(5,4) * t323;
t325 = cos(qJ(2));
t388 = Icges(5,1) * t325 + t306;
t304 = Icges(4,5) * t323;
t389 = Icges(4,1) * t325 + t304;
t672 = t388 + t389;
t324 = sin(qJ(1));
t506 = t323 * t324;
t281 = Icges(3,4) * t506;
t504 = t324 * t325;
t326 = cos(qJ(1));
t533 = Icges(3,5) * t326;
t152 = Icges(3,1) * t504 - t281 - t533;
t148 = Icges(5,5) * t326 + t324 * t388;
t150 = -Icges(4,4) * t326 + t324 * t389;
t660 = -t148 - t150;
t654 = -t152 + t660;
t536 = Icges(3,4) * t323;
t246 = Icges(3,1) * t325 - t536;
t153 = Icges(3,5) * t324 + t246 * t326;
t664 = t674 * t324 + t672 * t326;
t653 = t153 + t664;
t308 = Icges(3,4) * t325;
t532 = Icges(4,5) * t325;
t534 = Icges(5,4) * t325;
t644 = t675 * t323 - t532 - t534;
t665 = Icges(3,1) * t323 + t308 + t644;
t503 = t325 * t326;
t671 = (Icges(5,4) + Icges(4,5)) * t503;
t670 = t676 * t324;
t232 = Icges(4,3) * t323 + t532;
t138 = -Icges(4,6) * t326 + t232 * t324;
t236 = Icges(5,2) * t323 + t534;
t142 = Icges(5,6) * t326 + t236 * t324;
t662 = t138 + t142;
t505 = t323 * t326;
t661 = t673 * t505 + t670 + t671;
t234 = Icges(3,5) * t325 - Icges(3,6) * t323;
t141 = Icges(3,3) * t324 + t234 * t326;
t238 = Icges(4,4) * t325 + Icges(4,6) * t323;
t145 = Icges(4,2) * t324 + t238 * t326;
t669 = t141 + t145;
t643 = -t673 * t325 + t304 + t306;
t637 = -Icges(3,2) * t325 - t536 + t643;
t387 = -Icges(3,2) * t323 + t308;
t659 = t232 + t236;
t667 = -t387 + t659;
t652 = (-Icges(3,6) + t676) * t325 + (-Icges(3,5) - t674) * t323;
t663 = t246 + t672;
t522 = Icges(3,3) * t326;
t140 = Icges(3,5) * t504 - Icges(3,6) * t506 - t522;
t527 = Icges(3,6) * t326;
t146 = Icges(3,4) * t504 - Icges(3,2) * t506 - t527;
t514 = t146 * t323;
t375 = -t152 * t325 + t514;
t378 = t142 * t323 + t148 * t325;
t144 = -Icges(4,2) * t326 + t238 * t324;
t515 = t144 * t326;
t230 = Icges(5,5) * t325 + Icges(5,6) * t323;
t136 = Icges(5,3) * t326 + t230 * t324;
t517 = t136 * t326;
t382 = t138 * t323 + t150 * t325;
t597 = t324 * t382;
t611 = t324 * t378 - t515 + t517 + t597;
t587 = -t140 * t326 - t324 * t375 + t611;
t147 = Icges(3,6) * t324 + t326 * t387;
t119 = t153 * t504;
t421 = t141 * t326 - t119;
t50 = -t147 * t506 - t421;
t137 = -Icges(5,3) * t324 + t230 * t326;
t134 = t326 * t137;
t610 = -t145 * t326 + t664 * t504 + t661 * t506 + t134;
t586 = t50 + t610;
t658 = t669 * t324 + t653 * t503 + t661 * t505;
t132 = t324 * t144;
t657 = -t324 * t140 + t654 * t503 - t662 * t505 - t132;
t656 = -t146 + t662;
t655 = -t147 + t661;
t651 = t637 * qJD(2);
t650 = t665 * qJD(2);
t646 = t323 * t637 + t665 * t325;
t518 = t136 * t324;
t585 = -t146 * t505 - t518 - t657;
t584 = -t137 * t324 - t147 * t505 + t658;
t649 = t667 * qJD(2);
t648 = t663 * qJD(2);
t647 = t230 - t234 - t238;
t645 = -t323 * t665 + t637 * t325;
t580 = t652 * t326;
t581 = t652 * t324;
t540 = -rSges(5,3) - qJ(4);
t608 = pkin(3) * t503 + t324 * t540;
t484 = rSges(5,1) * t503 + rSges(5,2) * t505 + t608;
t601 = t324 * t646 + t580;
t600 = t326 * t646 - t581;
t642 = -t651 * t326 + (t324 * t387 - t527 - t662) * qJD(1);
t641 = -t651 * t324 + (t326 * t659 - t147 + t670) * qJD(1);
t640 = -t650 * t326 + (-t246 * t324 + t533 + t660) * qJD(1);
t639 = -qJD(1) * t653 + t324 * t650;
t638 = -t323 * t653 + t325 * t655;
t583 = t323 * t654 + t325 * t656;
t636 = t652 * qJD(2);
t513 = t147 * t323;
t635 = t323 * t661 + t325 * t653 - t513;
t634 = t375 - t378 - t382;
t633 = (t136 - t144) * qJD(1);
t632 = -qJD(1) * t652 + t645 * qJD(2) + t649 * t323 + t648 * t325;
t631 = t584 * t324 - t585 * t326;
t630 = t586 * t324 - t587 * t326;
t629 = (-t137 + t669) * qJD(1);
t628 = t646 * qJD(1) + t647 * qJD(2);
t627 = t324 * (t637 * t326 + t653) - t326 * (-Icges(3,2) * t504 + t643 * t324 - t281 - t654);
t626 = t326 ^ 2;
t625 = t600 * qJD(1);
t624 = t638 * qJD(2) + t642 * t323 + t640 * t325 + t629;
t623 = -qJD(1) * t140 - t583 * qJD(2) - t641 * t323 + t639 * t325 + t633;
t622 = -t656 * t326 + (t644 * t326 - t675 * t505 + t655 + t671) * t324;
t621 = -t667 + t665;
t620 = t637 + t663;
t617 = t601 * qJD(1);
t616 = t634 * qJD(1) + t636 * t324 + t629;
t615 = t633 + t636 * t326 + (-t234 * t324 + t522 - t635) * qJD(1);
t447 = qJD(1) * qJD(2);
t226 = -qJDD(2) * t326 + t324 * t447;
t451 = qJD(3) * t325;
t299 = t323 * qJ(3);
t593 = t325 * pkin(2) + t299;
t170 = qJD(2) * t593 - t451;
t310 = t323 * rSges(5,2);
t594 = t325 * rSges(5,1) + t310;
t481 = -qJD(2) * t594 - t170;
t554 = pkin(3) * t325;
t348 = (-t554 * qJD(2) + t481) * qJD(2);
t298 = qJD(3) * t323;
t431 = t324 * t298;
t448 = qJD(4) * t326;
t362 = -t431 - t448;
t199 = t593 * t324;
t319 = t326 * pkin(5);
t257 = pkin(1) * t324 - t319;
t478 = -t199 - t257;
t567 = -rSges(5,3) * t326 - t324 * t594;
t486 = pkin(3) * t504 + qJ(4) * t326 - t567;
t414 = t478 - t486;
t248 = rSges(5,1) * t323 - rSges(5,2) * t325;
t555 = pkin(3) * t323;
t422 = t248 + t555;
t247 = pkin(2) * t323 - qJ(3) * t325;
t446 = qJD(2) * qJD(3);
t609 = qJDD(3) * t323 + t325 * t446;
t442 = t226 * t247 + t326 * t609;
t454 = qJD(2) * t324;
t435 = t323 * t454;
t266 = pkin(3) * t435;
t418 = -t266 + t448;
t498 = -t248 * t454 + t418 + (t326 * t594 + t608) * qJD(1);
t258 = t326 * pkin(1) + t324 * pkin(5);
t224 = t258 * qJD(1);
t453 = qJD(2) * t325;
t455 = qJD(1) * t326;
t175 = t323 * t455 + t324 * t453;
t293 = pkin(2) * t503;
t267 = pkin(2) * t435;
t406 = -t267 + t431;
t83 = qJ(3) * t175 + qJD(1) * t293 + t406;
t598 = -t224 - t83;
t2 = -qJDD(4) * t324 + t422 * t226 + t348 * t326 + t414 * qJDD(1) + (t362 - t498 + t598) * qJD(1) + t442;
t614 = g(1) - t2;
t225 = qJDD(2) * t324 + t326 * t447;
t204 = qJ(3) * t505 + t293;
t450 = qJD(3) * t326;
t272 = t323 * t450;
t297 = pkin(5) * t455;
t456 = qJD(1) * t324;
t477 = qJD(1) * (-pkin(1) * t456 + t297) + qJDD(1) * t258;
t452 = qJD(2) * t326;
t434 = t323 * t452;
t349 = -t325 * t456 - t434;
t437 = t323 * t456;
t433 = t325 * t452;
t262 = qJ(3) * t433;
t468 = t262 + t272;
t82 = pkin(2) * t349 - qJ(3) * t437 + t468;
t351 = qJDD(1) * t204 + t477 + t609 * t324 + (t272 + t82) * qJD(1);
t412 = -t247 - t422;
t264 = rSges(5,2) * t433;
t449 = qJD(4) * t324;
t499 = -rSges(5,1) * t434 + pkin(3) * t349 - qJ(4) * t455 + qJD(1) * t567 + t264 - t449;
t3 = qJDD(4) * t326 + t484 * qJDD(1) + t499 * qJD(1) + t412 * t225 + (-qJD(1) * qJD(4) + t348) * t324 + t351;
t613 = -g(2) + t3;
t444 = -rSges(5,1) - pkin(2) - pkin(3);
t612 = t444 * t325 - pkin(1);
t403 = t325 * rSges(4,1) + t323 * rSges(4,3);
t469 = -t593 - t403;
t607 = qJD(2) * t630 + t617;
t606 = qJD(2) * t631 + t625;
t605 = t634 * qJD(2) + t639 * t323 + t641 * t325;
t604 = t635 * qJD(2) + t640 * t323 - t642 * t325;
t603 = -t324 * t628 + t326 * t632;
t602 = t324 * t632 + t326 * t628;
t599 = t515 + t658;
t228 = qJD(1) * t257;
t596 = -qJD(1) * t199 - t228;
t417 = t272 - t449;
t592 = t412 * t452 + t417;
t591 = -t299 - t310 + t612;
t312 = t324 * rSges(4,2);
t165 = rSges(4,1) * t503 + rSges(4,3) * t505 + t312;
t413 = t204 + t258;
t85 = t413 + t165;
t588 = qJD(1) * t486;
t579 = t593 + t594 + t554;
t578 = -t323 * t627 + t622 * t325;
t577 = (-t323 * t621 + t325 * t620) * qJD(1);
t576 = t616 * t626 + (t624 * t324 + (-t615 + t623) * t326) * t324;
t575 = t623 * t626 + (t615 * t324 + (-t616 + t624) * t326) * t324;
t574 = t647 * qJD(1);
t563 = -pkin(3) * t453 + qJD(2) * t579 + t481;
t562 = m(4) / 0.2e1;
t561 = m(5) / 0.2e1;
t560 = t225 / 0.2e1;
t559 = t226 / 0.2e1;
t556 = -rSges(4,1) - pkin(2);
t553 = g(2) * t324;
t550 = rSges(3,1) * t325;
t250 = rSges(3,1) * t323 + rSges(3,2) * t325;
t203 = t250 * t326;
t311 = t324 * rSges(3,3);
t166 = rSges(3,1) * t503 - rSges(3,2) * t505 + t311;
t113 = t166 + t258;
t64 = qJD(1) * t113 - t250 * t454;
t547 = t203 * t64;
t546 = t323 * rSges(4,1);
t436 = t250 * t452;
t465 = rSges(3,2) * t506 + t326 * rSges(3,3);
t163 = rSges(3,1) * t504 - t465;
t485 = -t163 - t257;
t63 = qJD(1) * t485 - t436;
t545 = t324 * t63;
t544 = t326 * t63;
t542 = -rSges(5,2) - qJ(3);
t541 = -rSges(4,3) - qJ(3);
t483 = -t165 - t204;
t482 = t324 * t199 + t326 * t204;
t480 = -t403 * qJD(2) - t170;
t276 = qJ(3) * t503;
t200 = -pkin(2) * t505 + t276;
t479 = qJD(1) * t200 + t324 * t451;
t249 = -rSges(4,3) * t325 + t546;
t470 = -t247 - t249;
t467 = rSges(4,2) * t455 + rSges(4,3) * t433;
t466 = rSges(3,2) * t437 + rSges(3,3) * t455;
t463 = t324 ^ 2 + t626;
t443 = t199 * t455 + t324 * t83 + t326 * t82;
t274 = qJ(3) * t504;
t195 = -pkin(2) * t506 + t274;
t441 = t195 * t454 + t200 * t452 + t298;
t316 = t326 * rSges(4,2);
t162 = t324 * t403 - t316;
t440 = -t162 + t478;
t439 = -t204 - t484;
t438 = t556 * t326;
t430 = -pkin(1) - t550;
t426 = -t454 / 0.2e1;
t425 = t454 / 0.2e1;
t424 = -t452 / 0.2e1;
t423 = t452 / 0.2e1;
t420 = -t140 + t513;
t419 = qJD(2) * t480;
t415 = t444 * t326;
t404 = t199 * t454 + t204 * t452 - t451;
t256 = rSges(2,1) * t326 - rSges(2,2) * t324;
t251 = rSges(2,1) * t324 + rSges(2,2) * t326;
t255 = -rSges(3,2) * t323 + t550;
t391 = -t324 * t64 - t544;
t106 = rSges(3,1) * t349 - rSges(3,2) * t433 + t466;
t198 = t250 * t324;
t109 = -qJD(2) * t198 + (t255 * t326 + t311) * qJD(1);
t384 = t106 * t326 + t109 * t324;
t373 = t163 * t324 + t166 * t326;
t363 = t452 * t470 + t272;
t361 = -qJDD(3) * t325 + t225 * t199 + t323 * t446 + t82 * t452 + t83 * t454;
t344 = t323 * t541 + t325 * t556 - pkin(1);
t288 = rSges(5,2) * t503;
t287 = rSges(4,3) * t503;
t284 = rSges(5,2) * t504;
t283 = rSges(4,3) * t504;
t273 = t325 * t450;
t219 = t255 * qJD(2);
t206 = t247 * t456;
t205 = t247 * t454;
t202 = -rSges(4,1) * t505 + t287;
t201 = -rSges(5,1) * t505 + t288;
t197 = -rSges(4,1) * t506 + t283;
t196 = -rSges(5,1) * t506 + t284;
t176 = t433 - t437;
t174 = t463 * t323 * qJD(2);
t108 = -t249 * t454 + (t326 * t403 + t312) * qJD(1);
t105 = rSges(4,1) * t349 - rSges(4,3) * t437 + t467;
t59 = t373 * qJD(2);
t44 = -t205 + (-qJD(2) * t249 + t298) * t324 + t85 * qJD(1);
t43 = qJD(1) * t440 + t363;
t42 = (t162 * t324 + t165 * t326) * qJD(2) + t404;
t41 = -t205 + (-qJD(2) * t248 + t298) * t324 + (t258 - t439) * qJD(1) + t418;
t40 = qJD(1) * t414 + t592;
t39 = (t324 * t486 + t326 * t484) * qJD(2) + t404;
t38 = qJD(1) * t106 + qJDD(1) * t166 - t219 * t454 - t225 * t250 + t477;
t37 = -t219 * t452 + t226 * t250 + t485 * qJDD(1) + (-t109 - t224) * qJD(1);
t18 = qJD(1) * t105 + qJDD(1) * t165 + t225 * t470 + t324 * t419 + t351;
t17 = t226 * t249 + t326 * t419 + t440 * qJDD(1) + (-t108 - t431 + t598) * qJD(1) + t442;
t4 = t162 * t225 + t483 * t226 + (t105 * t326 + t108 * t324) * qJD(2) + t361;
t1 = t486 * t225 + t439 * t226 + (t498 * t324 + t499 * t326) * qJD(2) + t361;
t5 = [-m(2) * (-g(1) * t251 + g(2) * t256) + (((t50 - t119 + (t141 + t514) * t326 + t657) * t326 + (-t597 + (-t137 - t378) * t324 + t599 + t611) * t324) * qJD(2) + t625) * t423 + (t646 * qJD(2) + t648 * t323 - t649 * t325) * qJD(1) + (-(-t40 - t588 + t592 + t596) * t41 + t40 * (t266 + t267 + t362) + t41 * (t262 + t264 + t297 + t417) + (t40 * t542 * t504 + (t40 * rSges(5,1) * t324 + t41 * t415) * t323) * qJD(2) + ((t40 * t591 + t41 * t540) * t326 + (t40 * (-pkin(5) - t540) + t591 * t41) * t324) * qJD(1) + t613 * (t413 + t484) - t614 * (t319 + t540 * t326 + (t323 * t542 + t612) * t324)) * m(5) + (-t43 * t406 + t44 * (t297 + t467 + t468) + (t44 * t323 * t438 + t43 * (t325 * t541 + t546) * t324) * qJD(2) + (t43 * t344 * t326 + (t43 * (-rSges(4,2) - pkin(5)) + t44 * (-pkin(1) + t469)) * t324) * qJD(1) - (-qJD(1) * t162 + t363 - t43 + t596) * t44 + (t18 - g(2)) * t85 + (t17 - g(1)) * (t324 * t344 + t316 + t319)) * m(4) + (-(-qJD(1) * t163 - t228 - t436 - t63) * t64 + t64 * (t297 + t466) + (t250 * t545 - t547) * qJD(2) + ((-pkin(1) - t255) * t544 + (t63 * (-rSges(3,3) - pkin(5)) + t64 * t430) * t324) * qJD(1) + (t38 - g(2)) * t113 + (t37 - g(1)) * (t430 * t324 + t319 + t465)) * m(3) + (m(2) * (t251 ^ 2 + t256 ^ 2) + Icges(2,3) - t645) * qJDD(1) + (-t638 + t600) * t560 + (-t583 + t601) * t559 + (((t326 * t420 + t517 + t584 - t599) * t326 + (t324 * t420 - t132 + t134 + t421 + t518 + t585 - t610) * t324) * qJD(2) + t607 - t617) * t426 + (t603 + t604) * t425 + (t602 - t605 + t606) * t424; t631 * t560 + t630 * t559 + (qJD(1) * t603 + t575 * qJD(2) + qJDD(1) * t600 + t584 * t225 + t585 * t226) * t324 / 0.2e1 - (qJD(1) * t602 + t576 * qJD(2) + qJDD(1) * t601 + t586 * t225 + t587 * t226) * t326 / 0.2e1 - ((t622 * t323 + t325 * t627) * qJD(2) + (t620 * t323 + t621 * t325) * qJD(1)) * qJD(1) / 0.2e1 + (t605 * t326 + t604 * t324 + (-t324 * t583 - t326 * t638) * qJD(1)) * qJD(1) / 0.2e1 + (-t324 * t638 + t326 * t583) * qJDD(1) / 0.2e1 + t607 * t456 / 0.2e1 + t606 * t455 / 0.2e1 + ((t454 * t580 - t574) * t324 + ((-t324 * t581 + t578) * qJD(2) + t577) * t326) * t426 + ((t324 * t585 + t326 * t584) * qJD(1) + t575) * t425 + ((t324 * t587 + t326 * t586) * qJD(1) + t576) * t424 + ((t452 * t581 + t574) * t326 + ((-t326 * t580 + t578) * qJD(2) + t577) * t324) * t423 + (-g(1) * (t276 + t288) - g(2) * (t274 + t284) - g(3) * t579 - (g(1) * t415 + t444 * t553) * t323 + t1 * t482 + (t1 * t486 + t3 * t412) * t324 + (t1 * t484 + t2 * t412) * t326 + (-t479 + t563 * t324 + (pkin(3) * t505 + t326 * t412 - t201) * qJD(1)) * t41 + (t206 - t273 + t563 * t326 + (t248 * t324 + t195 + t196) * qJD(1)) * t40 + (-t441 - (t196 * t324 + t201 * t326 - t463 * t555) * qJD(2) + t443 + (qJD(1) * t439 + t498) * t324 + (t499 + t588) * t326) * t39) * m(5) + (-g(1) * (t276 + t287) - g(2) * (t274 + t283) + g(3) * t469 - (g(1) * t438 + t553 * t556) * t323 - t43 * (t273 + (-t195 - t197) * qJD(1)) - t44 * (qJD(1) * t202 + t479) - t42 * t441 - ((t42 * t202 + t43 * t469) * t326 + (t42 * t197 + t44 * t469) * t324) * qJD(2) + t43 * t206 + t4 * t482 + t42 * t443 + (t17 * t470 + t43 * t480 + t4 * t165 + t42 * t105 + (t42 * t162 + t44 * t470) * qJD(1)) * t326 + (t18 * t470 + t44 * t480 + t4 * t162 + t42 * t108 + (t43 * t249 + t42 * t483) * qJD(1)) * t324) * m(4) + (-(t198 * t63 - t547) * qJD(1) - (t59 * (-t198 * t324 - t203 * t326) + t391 * t255) * qJD(2) + (qJD(2) * t384 + t163 * t225 - t166 * t226) * t373 + t59 * ((t163 * t326 - t166 * t324) * qJD(1) + t384) + t391 * t219 + (-t38 * t324 - t37 * t326 + (-t326 * t64 + t545) * qJD(1)) * t250 + g(1) * t203 + g(2) * t198 - g(3) * t255) * m(3); (-m(4) - m(5)) * (-g(3) * t325 + (g(1) * t326 + t553) * t323) - m(4) * (t174 * t42 + t175 * t44 + t176 * t43) - m(5) * (t174 * t39 + t175 * t41 + t176 * t40) + 0.2e1 * ((t43 * t452 + t44 * t454 - t4) * t562 + (t40 * t452 + t41 * t454 - t1) * t561) * t325 + 0.2e1 * ((qJD(2) * t42 + t17 * t326 + t18 * t324 - t43 * t456 + t44 * t455) * t562 + (qJD(2) * t39 + t2 * t326 + t3 * t324 - t40 * t456 + t41 * t455) * t561) * t323; (t324 * t614 + t613 * t326) * m(5);];
tau = t5;
