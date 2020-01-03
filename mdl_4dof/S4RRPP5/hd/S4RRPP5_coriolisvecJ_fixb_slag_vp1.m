% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:10
% EndTime: 2019-12-31 17:00:35
% DurationCPUTime: 21.46s
% Computational Cost: add. (4362->531), mult. (11594->631), div. (0->0), fcn. (9002->4), ass. (0->306)
t666 = Icges(4,4) - Icges(3,5);
t665 = Icges(4,5) - Icges(3,6);
t664 = Icges(4,1) + Icges(3,3);
t312 = sin(qJ(2));
t314 = cos(qJ(2));
t663 = t665 * t312 - t314 * t666;
t228 = Icges(5,4) * t312 + Icges(5,5) * t314;
t313 = sin(qJ(1));
t315 = cos(qJ(1));
t141 = Icges(5,1) * t313 + t228 * t315;
t635 = t664 * t313 + t663 * t315;
t662 = t141 + t635;
t298 = Icges(5,6) * t312;
t223 = -Icges(5,2) * t314 + t298;
t519 = Icges(4,6) * t312;
t383 = Icges(4,3) * t314 + t519;
t661 = t223 - t383;
t517 = Icges(5,6) * t314;
t518 = Icges(4,6) * t314;
t652 = t517 - t518 + (-Icges(4,2) - Icges(5,3)) * t312;
t660 = t664 * t315;
t500 = t313 * t314;
t502 = t312 * t313;
t654 = t500 * t666 - t665 * t502 + t660;
t530 = Icges(3,4) * t312;
t235 = Icges(3,1) * t314 - t530;
t132 = Icges(3,5) * t313 + t235 * t315;
t381 = Icges(5,3) * t314 + t298;
t133 = Icges(5,5) * t313 + t315 * t381;
t634 = t132 + t133;
t384 = -Icges(4,3) * t312 + t518;
t135 = Icges(4,5) * t313 - t315 * t384;
t499 = t314 * t315;
t271 = Icges(5,6) * t499;
t501 = t312 * t315;
t527 = Icges(5,4) * t313;
t137 = Icges(5,2) * t501 + t271 + t527;
t659 = t135 + t137;
t658 = (-Icges(5,4) - Icges(4,5)) * t314 + (-Icges(4,4) + Icges(5,5)) * t312;
t520 = Icges(3,6) * t315;
t129 = Icges(3,4) * t500 - Icges(3,2) * t502 - t520;
t272 = Icges(4,6) * t502;
t528 = Icges(4,4) * t315;
t140 = Icges(4,2) * t500 - t272 + t528;
t657 = t129 * t312 - t140 * t314;
t222 = Icges(5,2) * t312 + t517;
t301 = Icges(3,4) * t314;
t391 = -Icges(3,2) * t312 + t301;
t627 = t222 - t384 - t391;
t387 = Icges(4,2) * t314 - t519;
t656 = t235 + t381 + t387;
t226 = Icges(3,5) * t312 + Icges(3,6) * t314;
t641 = t226 + t658;
t655 = -t132 * t500 - t135 * t502;
t280 = Icges(3,4) * t502;
t525 = Icges(3,5) * t315;
t131 = Icges(3,1) * t500 - t280 - t525;
t524 = Icges(4,5) * t315;
t136 = Icges(4,6) * t500 - Icges(4,3) * t502 + t524;
t653 = t131 * t314 - t136 * t312 - t657;
t232 = Icges(3,2) * t314 + t530;
t632 = -t223 + t232;
t644 = -t383 - t632;
t234 = Icges(3,1) * t312 + t301;
t370 = t232 * t312 - t234 * t314;
t642 = t312 * t661 - t314 * t652 - t370;
t399 = -t133 * t500 - t137 * t502 + t141 * t315;
t130 = Icges(3,6) * t313 + t315 * t391;
t273 = Icges(4,6) * t501;
t529 = Icges(4,4) * t313;
t139 = -Icges(4,2) * t499 + t273 + t529;
t650 = -t635 * t315 - t655;
t612 = -t130 * t502 - t139 * t500 + t650;
t651 = -t399 + t612;
t649 = -t131 * t499 + t136 * t501 + t654 * t313;
t648 = t658 * t315;
t647 = t662 * t313 + t634 * t499 + t659 * t501;
t646 = t656 * qJD(2);
t645 = t627 * qJD(2);
t618 = t228 + t663;
t643 = t234 - t652;
t568 = t641 * t313;
t640 = t129 * t501 - t140 * t499 + t649;
t639 = t652 * t500 - t502 * t661 + t648;
t638 = t130 * t312 + t139 * t314;
t531 = Icges(5,1) * t315;
t142 = -Icges(5,4) * t502 - Icges(5,5) * t500 + t531;
t506 = t142 * t315;
t269 = Icges(5,6) * t502;
t523 = Icges(5,5) * t315;
t134 = -Icges(5,3) * t500 - t269 + t523;
t270 = Icges(5,6) * t500;
t526 = Icges(5,4) * t315;
t138 = -Icges(5,2) * t502 - t270 + t526;
t377 = t134 * t314 + t138 * t312;
t577 = t313 * t377;
t54 = t506 - t577;
t637 = t653 * t313 + t654 * t315 + t54;
t601 = -t130 * t501 - t139 * t499 + t647;
t636 = t642 * t315 + t568;
t633 = t135 - t130;
t631 = t129 + t136 + t138;
t630 = -t137 - t633;
t629 = t131 - t134 + t140;
t628 = -t139 + t634;
t625 = t646 * t314 + t645 * t312 + (-t312 * t643 + t314 * t644) * qJD(2) + t641 * qJD(1);
t624 = t644 * qJD(2);
t623 = t643 * qJD(2);
t622 = t377 - t653;
t621 = t659 * t312 + t634 * t314 - t638;
t620 = t642 * qJD(1) - qJD(2) * t618;
t548 = rSges(5,1) + pkin(3);
t619 = t636 * qJD(1);
t617 = t641 * qJD(2);
t616 = (t313 * t651 - t637 * t315) * qJD(2);
t123 = t313 * t142;
t50 = -t134 * t499 - t138 * t501 - t123;
t592 = t315 * t50;
t615 = (t601 * t313 + t640 * t315 - t592) * qJD(2);
t503 = t226 * t315;
t71 = -t313 * t370 - t503;
t614 = (t71 - t639) * qJD(1);
t613 = t313 * (Icges(4,3) * t499 + t632 * t315 + t273 - t628) + t315 * (t269 - t272 - t280 + (-Icges(3,2) - Icges(5,2) - Icges(4,3)) * t500 + t629);
t591 = rSges(5,3) + qJ(4);
t611 = (-t617 * t313 + (t622 + t662) * qJD(1)) * t315;
t610 = t638 + t654;
t514 = qJ(3) * t312;
t241 = pkin(2) * t314 + t514;
t191 = t241 * t313;
t161 = qJD(1) * t191;
t308 = t315 * pkin(5);
t246 = pkin(1) * t313 - t308;
t217 = qJD(1) * t246;
t449 = qJD(1) * t315;
t294 = pkin(5) * t449;
t443 = qJD(3) * t315;
t262 = t312 * t443;
t445 = qJD(2) * t315;
t427 = t314 * t445;
t464 = qJ(3) * t427 + t262;
t609 = t294 + t464 + t161 + t217;
t450 = qJD(1) * t313;
t608 = t614 + t616;
t607 = t615 + t619;
t606 = (t624 * t315 + (t627 * t313 + t520 - t524 - t526) * qJD(1)) * t314 + (-t623 * t315 + (-t313 * t656 + t523 + t525 - t528) * qJD(1)) * t312 + t621 * qJD(2);
t605 = (-t624 * t313 + (t222 * t315 + t527 + t633) * qJD(1)) * t314 + (t623 * t313 + (-t315 * t387 + t529 - t634) * qJD(1)) * t312 + t622 * qJD(2);
t604 = -t620 * t313 + t625 * t315;
t603 = t625 * t313 + t620 * t315;
t602 = t50 - t640;
t600 = t629 * t312 + t631 * t314;
t599 = t628 * t312 + t630 * t314;
t598 = (Icges(5,3) * t502 - t270 + t631) * t315 + (-Icges(5,3) * t501 + t271 - t630) * t313;
t597 = -t506 + t647;
t596 = -t627 + t643;
t595 = t656 + t644;
t539 = t312 * rSges(5,3);
t394 = rSges(5,2) * t314 - t539;
t192 = t394 * t313;
t586 = qJ(4) * t499 + t313 * t548;
t428 = t312 * t445;
t585 = t314 * t450 + t428;
t582 = -t617 * t315 + (-t313 * t618 + t531 - t621 + t660) * qJD(1);
t580 = 0.2e1 * qJD(2);
t438 = -pkin(2) - t591;
t579 = t438 * t314 - pkin(1);
t289 = pkin(2) * t499;
t196 = qJ(3) * t501 + t289;
t444 = qJD(3) * t314;
t447 = qJD(2) * t313;
t398 = t191 * t447 + t196 * t445 - t444;
t442 = qJD(4) * t312;
t573 = -rSges(5,2) * t502 - t500 * t591;
t479 = t548 * t315 + t573;
t481 = rSges(5,2) * t501 + rSges(5,3) * t499 + t586;
t37 = t442 + (-t313 * t479 + t315 * t481) * qJD(2) + t398;
t578 = qJD(2) * t37;
t296 = qJD(3) * t312;
t540 = t312 * rSges(4,2);
t395 = rSges(4,3) * t314 + t540;
t576 = (-qJD(2) * t395 - t296) * t313;
t304 = t313 * rSges(4,1);
t155 = -rSges(4,2) * t499 + rSges(4,3) * t501 + t304;
t247 = t315 * pkin(1) + t313 * pkin(5);
t406 = t196 + t247;
t575 = t155 + t406;
t457 = t313 ^ 2 + t315 ^ 2;
t440 = qJD(4) * t315;
t260 = t314 * t440;
t236 = pkin(2) * t312 - qJ(3) * t314;
t513 = qJ(4) * t312;
t414 = -t394 + t513;
t401 = -t236 - t414;
t365 = t401 * qJD(2);
t571 = t315 * t365 + t260 + t262;
t570 = rSges(5,2) * t427 + t449 * t548 + t260;
t569 = qJD(1) * t479;
t567 = t503 + t648;
t566 = t312 * t613 + t598 * t314;
t565 = (-t312 * t596 + t314 * t595) * qJD(1);
t564 = t618 * qJD(1);
t543 = rSges(5,2) * t312;
t237 = rSges(5,3) * t314 + t543;
t446 = qJD(2) * t314;
t160 = qJD(2) * t241 - t444;
t477 = -t237 * qJD(2) - t160;
t512 = qJ(4) * t314;
t554 = -qJ(4) * t446 - qJD(2) * (-t237 - t241 - t512) - t442 + t477;
t553 = 0.2e1 * t314;
t552 = m(4) / 0.2e1;
t551 = m(5) / 0.2e1;
t546 = qJD(1) / 0.2e1;
t545 = rSges(3,1) * t314;
t544 = rSges(4,2) * t314;
t542 = rSges(4,3) * t312;
t239 = rSges(3,1) * t312 + rSges(3,2) * t314;
t195 = t239 * t315;
t431 = t239 * t447;
t302 = t313 * rSges(3,3);
t153 = rSges(3,1) * t499 - rSges(3,2) * t501 + t302;
t482 = t153 + t247;
t62 = qJD(1) * t482 - t431;
t541 = t195 * t62;
t458 = rSges(3,2) * t502 + t315 * rSges(3,3);
t152 = rSges(3,1) * t500 - t458;
t430 = t239 * t445;
t61 = -t430 + (-t152 - t246) * qJD(1);
t538 = t313 * t61;
t537 = t315 * t61;
t536 = -rSges(5,2) - qJ(3);
t535 = -rSges(4,3) - qJ(3);
t216 = t247 * qJD(1);
t168 = t312 * t449 + t313 * t446;
t429 = t312 * t447;
t259 = pkin(2) * t429;
t81 = qJ(3) * t168 + qJD(1) * t289 + t296 * t313 - t259;
t534 = -t216 - t81;
t251 = qJ(4) * t429;
t441 = qJD(4) * t314;
t495 = qJD(2) * t192 + t313 * t441 - t251 + (t237 * t315 + t586) * qJD(1);
t433 = t312 * t450;
t494 = -rSges(5,2) * t433 - t585 * t591 + t570;
t480 = -t155 - t196;
t478 = t313 * t191 + t315 * t196;
t243 = t542 - t544;
t476 = -t243 * qJD(2) - t160;
t193 = t236 * t315;
t475 = -qJD(1) * t193 + t313 * t444;
t198 = t236 * t447;
t439 = qJD(2) * qJD(3);
t424 = t314 * t439;
t474 = qJD(1) * t198 + t315 * t424;
t473 = -t191 - t246;
t466 = -t236 + t395;
t465 = -t241 - t243;
t462 = rSges(3,2) * t433 + rSges(3,3) * t449;
t448 = qJD(2) * t312;
t80 = -pkin(2) * t585 - qJ(3) * t433 + t464;
t437 = t191 * t449 + t313 * t81 + t315 * t80;
t188 = t236 * t313;
t436 = -t188 * t447 - t193 * t445 + t296;
t435 = -t196 - t481;
t425 = -pkin(1) - t545;
t421 = -t447 / 0.2e1;
t418 = t445 / 0.2e1;
t416 = rSges(4,1) * t315 - rSges(4,3) * t502;
t415 = t315 * t466;
t409 = t312 * t439 + t81 * t447 + (t161 + t80) * t445;
t200 = qJD(1) * (-pkin(1) * t450 + t294);
t408 = t313 * t424 + t200 + (t262 + t80) * qJD(1);
t407 = rSges(4,1) * t449 + rSges(4,2) * t585 + rSges(4,3) * t427;
t396 = -rSges(3,2) * t312 + t545;
t393 = -t313 * t62 - t537;
t368 = -t441 - t296;
t367 = -pkin(1) - t241;
t364 = qJD(2) * t415 + t262;
t190 = t239 * t313;
t189 = t395 * t313;
t57 = (t152 * t313 + t153 * t315) * qJD(2);
t337 = (-qJD(2) * t512 - 0.2e1 * t442 + t477) * qJD(2);
t15 = t337 * t313 + ((t365 + t441) * t315 + t494) * qJD(1) + t408;
t16 = t337 * t315 + ((qJD(2) * t414 + t368) * t313 - t495 + t534) * qJD(1) + t474;
t345 = t15 * t313 + t16 * t315 + t578;
t263 = t314 * t443;
t212 = t396 * qJD(2);
t199 = t236 * t450;
t197 = t394 * t315;
t194 = t395 * t315;
t169 = t427 - t433;
t167 = t457 * t448;
t157 = rSges(4,2) * t500 + t416;
t105 = -rSges(4,3) * t433 + t407;
t103 = qJD(2) * t189 + (t243 * t315 + t304) * qJD(1);
t101 = -qJD(2) * t190 + (t315 * t396 + t302) * qJD(1);
t100 = -rSges(3,1) * t585 - rSges(3,2) * t427 + t462;
t44 = qJD(1) * t575 - t198 - t576;
t43 = (t157 + t473) * qJD(1) + t364;
t42 = -t212 * t445 + (-t101 - t216 + t431) * qJD(1);
t41 = -t212 * t447 + t200 + (t100 - t430) * qJD(1);
t40 = (t155 * t315 - t157 * t313) * qJD(2) + t398;
t39 = -t198 - t251 + (qJD(2) * t394 - t368) * t313 + (t247 - t435) * qJD(1);
t38 = (t473 + t479) * qJD(1) + t571;
t24 = t476 * t445 + (-t103 + t534 + t576) * qJD(1) + t474;
t23 = qJD(1) * t105 + (qJD(1) * t415 + t313 * t476) * qJD(2) + t408;
t2 = (t103 * t313 + t105 * t315 + (-t157 * t315 + t313 * t480) * qJD(1)) * qJD(2) + t409;
t1 = (t441 + t494 * t315 + t495 * t313 + (t313 * t435 - t315 * t479) * qJD(1)) * qJD(2) + t409;
t3 = [(t642 * qJD(2) + t646 * t312 - t645 * t314) * qJD(1) + (t16 * (t308 + t573) + t38 * (t251 + t259) + t15 * (t406 + t481) + (t16 * t548 + t38 * (t312 * t536 + t579) * qJD(1)) * t315 + (t16 * t367 + (t368 + (t314 * t536 + t539) * qJD(2) + (-pkin(5) - t548) * qJD(1)) * t38) * t313 + (t38 - t569 - t571 + t570 + t438 * t448 * t315 + (-t514 - t543 + t579) * t450 + t609) * t39) * m(5) + (t23 * t575 + (t259 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t314 + t535 * t312) * t449 + (-t296 + (t314 * t535 - t540) * qJD(2) + (-rSges(4,1) - pkin(5)) * qJD(1)) * t313) * t43 + (t308 + t416 + (t367 + t544) * t313) * t24 + (-pkin(2) * t428 - t364 + t407 + t43 + (-t157 + (t367 - t542) * t313) * qJD(1) + t609) * t44) * m(4) + (t42 * (t313 * t425 + t308 + t458) + t41 * t482 + t62 * (t294 + t462) + (t239 * t538 - t541) * qJD(2) + ((-pkin(1) - t396) * t537 + (t61 * (-rSges(3,3) - pkin(5)) + t62 * t425) * t313) * qJD(1) - (-qJD(1) * t152 - t217 - t430 - t61) * t62) * m(3) + ((-t592 + ((t635 + t657) * t315 + t612 + t649 + t655) * t315 + (t54 + t577 + t597) * t313) * qJD(2) + t619) * t418 + (t604 + t606) * t447 / 0.2e1 - (t603 - t605 + t607) * t445 / 0.2e1 + (t639 * qJD(1) + ((t315 * t610 - t597 + t601) * t315 + (t313 * t610 + t123 + t399 + t602 - t650) * t313) * qJD(2) + t608 - t614) * t421 + ((t71 + t600) * t313 + (t599 + t636) * t315) * qJD(2) * t546; (t1 * t478 + (-t1 * t479 + t15 * t401) * t313 + (t1 * t481 + t16 * t401) * t315 + (-t475 + (qJ(4) * t501 + t315 * t401 - t197) * qJD(1) + (t442 + t554) * t313) * t39 + (-qJD(1) * t188 + t312 * t440 + t554 * t315 + t199 - t263) * t38 + (t437 + (qJD(1) * t435 + t495) * t313 + (t494 - t569) * t315 - t436 - t441 - (t192 * t313 + t197 * t315 - t457 * t513) * qJD(2)) * t37) * m(5) + (t43 * t199 + t2 * t478 + t40 * t437 + (t24 * t466 + t43 * t476 + t2 * t155 + t40 * t105 + (-t40 * t157 + t44 * t466) * qJD(1)) * t315 + (t23 * t466 + t44 * t476 - t2 * t157 + t40 * t103 + (-t395 * t43 + t40 * t480) * qJD(1)) * t313 - t43 * (t263 + (t188 - t189) * qJD(1)) - t44 * (qJD(1) * t194 + t475) - t40 * t436 - ((t40 * t194 + t43 * t465) * t315 + (t40 * t189 + t44 * t465) * t313) * qJD(2)) * m(4) + (0.2e1 * t57 * (t100 * t315 + t101 * t313 + (t152 * t315 - t153 * t313) * qJD(1)) + t393 * t212 + (-t41 * t313 - t42 * t315 + (-t315 * t62 + t538) * qJD(1)) * t239 - (t190 * t61 - t541) * qJD(1) - (t57 * (-t190 * t313 - t195 * t315) + t393 * t396) * qJD(2)) * m(3) - ((t598 * t312 - t314 * t613) * qJD(2) + (t595 * t312 + t596 * t314) * qJD(1)) * qJD(1) / 0.2e1 + (t605 * t315 + t606 * t313 + (t313 * t600 + t315 * t599) * qJD(1)) * t546 + ((-t447 * t567 + t564) * t313 + ((t313 * t568 + t566) * qJD(2) + t565) * t315) * t421 + ((-t445 * t568 - t564) * t315 + ((t315 * t567 + t566) * qJD(2) + t565) * t313) * t418 + (t604 * qJD(1) + (t601 * t449 + (qJD(1) * t602 + t582 * t313 - t611) * t313) * t580) * t313 / 0.2e1 - (t603 * qJD(1) + ((qJD(1) * t651 + t611) * t315 + (qJD(1) * t637 - t582 * t315) * t313) * t580) * t315 / 0.2e1 + (t608 + t616) * t450 / 0.2e1 + (t607 + t615) * t449 / 0.2e1; -m(4) * (t167 * t40 + t168 * t44 + t169 * t43) - m(5) * (t167 * t37 + t168 * t39 + t169 * t38) + ((t43 * t445 + t44 * t447 - t2) * t552 + (t38 * t445 + t39 * t447 - t1) * t551) * t553 + 0.2e1 * ((qJD(2) * t40 + t23 * t313 + t24 * t315 - t43 * t450 + t44 * t449) * t552 + (-t38 * t450 + t39 * t449 + t345) * t551) * t312; m(5) * t1 * t312 + (t345 * t551 - m(5) * t457 * t578 / 0.2e1) * t553;];
tauc = t3(:);
