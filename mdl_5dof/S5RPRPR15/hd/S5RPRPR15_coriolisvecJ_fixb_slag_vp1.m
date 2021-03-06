% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR15_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR15_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:38
% EndTime: 2019-12-31 18:37:16
% DurationCPUTime: 30.79s
% Computational Cost: add. (15100->1010), mult. (27177->1332), div. (0->0), fcn. (25200->8), ass. (0->465)
t380 = cos(pkin(8));
t382 = sin(qJ(3));
t385 = cos(qJ(1));
t379 = sin(pkin(8));
t383 = sin(qJ(1));
t573 = t383 * t379;
t274 = t380 * t385 - t382 * t573;
t577 = t382 * t383;
t579 = t379 * t385;
t275 = t380 * t577 + t579;
t384 = cos(qJ(3));
t572 = t383 * t384;
t156 = Icges(5,4) * t275 + Icges(5,2) * t274 - Icges(5,6) * t572;
t159 = Icges(5,1) * t275 + Icges(5,4) * t274 - Icges(5,5) * t572;
t447 = Icges(4,5) * t382 + Icges(4,6) * t384;
t252 = Icges(4,3) * t385 + t383 * t447;
t597 = Icges(4,4) * t382;
t450 = Icges(4,2) * t384 + t597;
t254 = Icges(4,6) * t385 + t383 * t450;
t349 = Icges(4,4) * t572;
t592 = Icges(4,5) * t385;
t256 = Icges(4,1) * t577 + t349 + t592;
t700 = t274 * t156 + t275 * t159 + t385 * t252 + t254 * t572 + t256 * t577;
t153 = Icges(5,5) * t275 + Icges(5,6) * t274 - Icges(5,3) * t572;
t699 = -t153 * t572 + t700;
t576 = t382 * t385;
t276 = t379 * t576 + t380 * t383;
t277 = t380 * t576 - t573;
t571 = t384 * t385;
t155 = -Icges(5,5) * t277 + Icges(5,6) * t276 + Icges(5,3) * t571;
t158 = -Icges(5,4) * t277 + Icges(5,2) * t276 + Icges(5,6) * t571;
t161 = -Icges(5,1) * t277 + Icges(5,4) * t276 + Icges(5,5) * t571;
t568 = t274 * t158 + t275 * t161;
t253 = -Icges(4,3) * t383 + t385 * t447;
t255 = -Icges(4,6) * t383 + t385 * t450;
t596 = Icges(4,4) * t384;
t453 = Icges(4,1) * t382 + t596;
t257 = -Icges(4,5) * t383 + t385 * t453;
t93 = -t385 * t253 - t255 * t572 - t257 * t577;
t693 = -t155 * t572 + t568 + t93;
t51 = t153 * t571 + t276 * t156 - t277 * t159;
t230 = t383 * t252;
t433 = t254 * t384 + t256 * t382;
t94 = -t433 * t385 + t230;
t698 = t51 + t94;
t437 = -t158 * t276 + t161 * t277;
t431 = t255 * t384 + t257 * t382;
t655 = t431 * t385;
t95 = -t253 * t383 + t655;
t697 = t155 * t571 - t437 + t95;
t446 = Icges(5,5) * t380 - Icges(5,6) * t379;
t240 = Icges(5,3) * t382 + t384 * t446;
t449 = Icges(5,4) * t380 - Icges(5,2) * t379;
t242 = Icges(5,6) * t382 + t384 * t449;
t452 = Icges(5,1) * t380 - Icges(5,4) * t379;
t244 = Icges(5,5) * t382 + t384 * t452;
t319 = -Icges(4,2) * t382 + t596;
t321 = Icges(4,1) * t384 - t597;
t430 = t319 * t384 + t321 * t382;
t317 = Icges(4,5) * t384 - Icges(4,6) * t382;
t581 = t317 * t385;
t696 = -t240 * t572 + t242 * t274 + t244 * t275 + t383 * t430 + t581;
t378 = pkin(8) + qJ(5);
t363 = cos(t378);
t362 = sin(t378);
t574 = t383 * t362;
t248 = t363 * t385 - t382 * t574;
t249 = t362 * t385 + t363 * t577;
t119 = Icges(6,5) * t249 + Icges(6,6) * t248 - Icges(6,3) * t572;
t250 = t362 * t576 + t363 * t383;
t251 = t363 * t576 - t574;
t121 = -Icges(6,5) * t251 + Icges(6,6) * t250 + Icges(6,3) * t571;
t236 = Icges(6,4) * t251;
t124 = Icges(6,2) * t250 + Icges(6,6) * t571 - t236;
t235 = Icges(6,4) * t250;
t126 = Icges(6,1) * t251 - Icges(6,5) * t571 - t235;
t442 = -t124 * t250 - t126 * t251;
t595 = Icges(6,4) * t249;
t122 = Icges(6,2) * t248 - Icges(6,6) * t572 + t595;
t234 = Icges(6,4) * t248;
t125 = Icges(6,1) * t249 - Icges(6,5) * t572 + t234;
t617 = t248 * t122 + t249 * t125;
t695 = t442 + t617 + (-t119 * t383 - t121 * t385) * t384;
t130 = -rSges(6,1) * t251 + rSges(6,2) * t250 + rSges(6,3) * t571;
t522 = qJD(5) * t384;
t526 = qJD(3) * t385;
t305 = -t383 * t522 + t526;
t523 = qJD(5) * t382;
t352 = qJD(1) + t523;
t45 = t119 * t571 + t250 * t122 - t251 * t125;
t445 = Icges(6,5) * t363 - Icges(6,6) * t362;
t222 = Icges(6,3) * t382 + t384 * t445;
t593 = Icges(6,4) * t363;
t448 = -Icges(6,2) * t362 + t593;
t224 = Icges(6,6) * t382 + t384 * t448;
t594 = Icges(6,4) * t362;
t451 = Icges(6,1) * t363 - t594;
t226 = Icges(6,5) * t382 + t384 * t451;
t74 = t222 * t571 + t224 * t250 - t226 * t251;
t694 = -t305 * t45 - t352 * t74;
t638 = -pkin(1) - pkin(6);
t479 = -pkin(4) * t379 + t638;
t690 = t383 * t479;
t689 = -t230 - t568;
t348 = qJ(4) * t571;
t292 = pkin(3) * t576 - t348;
t366 = qJD(2) * t383;
t530 = qJD(1) * t385;
t540 = qJ(2) * t530 + t366;
t688 = -qJD(1) * t292 + t540;
t460 = rSges(6,1) * t363 - rSges(6,2) * t362;
t229 = rSges(6,3) * t382 + t384 * t460;
t521 = qJD(5) * t385;
t528 = qJD(3) * t383;
t304 = t384 * t521 + t528;
t327 = pkin(3) * t384 + qJ(4) * t382;
t296 = t327 * t528;
t365 = qJD(4) * t384;
t498 = t383 * t365;
t423 = t296 + t366 - t498;
t381 = -pkin(7) - qJ(4);
t570 = qJ(4) + t381;
t355 = pkin(4) * t380 + pkin(3);
t619 = pkin(3) - t355;
t644 = t382 * t570 + t384 * t619;
t687 = -t130 * t352 + t229 * t304 - t528 * t644 + t423;
t441 = t124 * t362 + t126 * t363;
t48 = t121 * t382 - t384 * t441;
t44 = -t121 * t572 + t248 * t124 - t126 * t249;
t686 = t696 * qJD(1);
t375 = t385 * rSges(4,3);
t264 = rSges(4,1) * t577 + rSges(4,2) * t572 + t375;
t329 = t385 * pkin(1) + t383 * qJ(2);
t376 = t385 * pkin(6);
t652 = t376 + t329;
t681 = t264 + t652;
t680 = (t383 * t697 + t385 * t698) * qJD(3);
t679 = (t693 * t383 + t385 * t699) * qJD(3);
t527 = qJD(3) * t384;
t496 = t383 * t527;
t678 = t382 * t530 + t496;
t531 = qJD(1) * t384;
t677 = t382 * t526 + t383 * t531;
t283 = t383 * t317;
t142 = t385 * t430 - t283;
t79 = t240 * t571 + t242 * t276 - t244 * t277;
t676 = (-t142 + t79) * qJD(1);
t286 = t319 * t385;
t400 = t383 * (t257 + t286) - t385 * (-Icges(4,2) * t577 + t256 + t349);
t287 = t321 * t383;
t288 = t321 * t385;
t401 = t383 * (t255 - t288) - t385 * (t254 - t287);
t439 = t153 * t385 + t155 * t383;
t675 = t400 * t384 + (-t401 - t439) * t382;
t533 = qJD(1) * t382;
t541 = t319 + t453;
t542 = -t450 + t321;
t673 = -(t382 * t541 - t384 * t542) * qJD(1) + t240 * t533;
t436 = t158 * t379 - t161 * t380;
t671 = 0.2e1 * qJD(3);
t670 = t679 + t686;
t669 = t676 + t680;
t172 = qJD(1) * t255 + t319 * t528;
t174 = qJD(1) * t257 + qJD(3) * t287;
t438 = t156 * t379 - t159 * t380;
t197 = -qJD(1) * t276 - t379 * t496;
t198 = qJD(1) * t277 + t380 * t496;
t501 = t382 * t528;
t503 = t384 * t530;
t280 = t501 - t503;
t85 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t280;
t87 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t280;
t89 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t280;
t668 = (-t379 * t87 + t380 * t89 + t174) * t384 + (t85 - t172) * t382 + (t153 * t384 + t382 * t438 - t433) * qJD(3);
t171 = qJD(1) * t254 - qJD(3) * t286;
t173 = -qJD(3) * t288 + (t383 * t453 + t592) * qJD(1);
t499 = t384 * t526;
t195 = qJD(1) * t274 + t379 * t499;
t196 = qJD(1) * t275 - t380 * t499;
t84 = Icges(5,5) * t196 + Icges(5,6) * t195 - Icges(5,3) * t677;
t86 = Icges(5,4) * t196 + Icges(5,2) * t195 - Icges(5,6) * t677;
t88 = Icges(5,1) * t196 + Icges(5,4) * t195 - Icges(5,5) * t677;
t667 = (-t379 * t86 + t380 * t88 + t173) * t384 + (t84 - t171) * t382 + (t155 * t384 + t382 * t436 + t431) * qJD(3);
t239 = Icges(5,3) * t384 - t382 * t446;
t215 = t239 * qJD(3);
t241 = Icges(5,6) * t384 - t382 * t449;
t216 = t241 * qJD(3);
t243 = Icges(5,5) * t384 - t382 * t452;
t217 = t243 * qJD(3);
t397 = t430 * qJD(1) - t447 * qJD(3);
t308 = t450 * qJD(3);
t309 = t453 * qJD(3);
t647 = qJD(1) * t317 + qJD(3) * (t319 * t382 - t321 * t384) + t308 * t384 + t309 * t382;
t666 = t195 * t242 + t196 * t244 + t215 * t571 + t216 * t276 - t217 * t277 - t240 * t677 + t397 * t383 + t385 * t647;
t665 = t197 * t242 + t198 * t244 - t215 * t572 + t216 * t274 + t217 * t275 + t240 * t280 - t383 * t647 + t397 * t385;
t432 = t254 * t382 - t256 * t384;
t664 = -t153 * t382 + t384 * t438 + t432;
t139 = t255 * t382 - t257 * t384;
t663 = t155 * t382 - t384 * t436 + t139;
t73 = -t222 * t572 + t224 * t248 + t226 * t249;
t662 = t304 * t44 + t73 * t352;
t661 = t382 * t619;
t221 = Icges(6,3) * t384 - t382 * t445;
t435 = t224 * t362 - t226 * t363;
t443 = t122 * t362 - t125 * t363;
t388 = t304 * (-t222 * t385 + t441) + t305 * (t222 * t383 + t443) + t352 * (t221 + t435);
t660 = t388 * t384;
t463 = rSges(5,1) * t380 - rSges(5,2) * t379;
t247 = rSges(5,3) * t382 + t384 * t463;
t657 = (qJD(3) * t247 - t365) * t385;
t656 = (-qJD(3) * t644 - t365) * t385;
t369 = t385 * qJ(2);
t324 = pkin(1) * t383 - t369;
t624 = pkin(6) * t383;
t484 = -t324 - t624;
t481 = -rSges(3,2) * t385 + t383 * rSges(3,3);
t653 = t329 + t481;
t434 = t242 * t379 - t244 * t380;
t651 = (t239 + t434) * qJD(1);
t478 = qJD(5) + t533;
t650 = t383 * t478 - t499;
t649 = t385 * t478 + t496;
t535 = qJD(1) * t253;
t648 = qJD(3) * t139 + t171 * t384 + t173 * t382 + t535;
t536 = qJD(1) * t252;
t646 = qJD(3) * t432 - t172 * t384 - t174 * t382 + t536;
t643 = t383 * t436 + t385 * t438;
t259 = (-Icges(6,2) * t363 - t594) * t384;
t395 = t304 * (Icges(6,2) * t251 - t126 + t235) + t305 * (-Icges(6,2) * t249 + t125 + t234) + t352 * (t226 + t259);
t260 = (-Icges(6,1) * t362 - t593) * t384;
t642 = t304 * (-Icges(6,1) * t250 + t124 - t236) + t305 * (-Icges(6,1) * t248 + t122 + t595) + t352 * (t224 - t260);
t641 = t383 ^ 2;
t640 = m(5) / 0.2e1;
t639 = m(6) / 0.2e1;
t517 = qJD(3) * qJD(5);
t494 = t382 * t517;
t219 = -qJD(1) * t304 + t383 * t494;
t637 = t219 / 0.2e1;
t220 = qJD(1) * t305 - t385 * t494;
t636 = t220 / 0.2e1;
t635 = -t304 / 0.2e1;
t634 = t304 / 0.2e1;
t633 = -t305 / 0.2e1;
t632 = t305 / 0.2e1;
t631 = -t352 / 0.2e1;
t630 = t352 / 0.2e1;
t629 = t382 / 0.2e1;
t628 = t383 / 0.2e1;
t626 = rSges(3,2) - pkin(1);
t625 = pkin(3) * t382;
t427 = t383 * t352;
t117 = -t362 * t649 - t363 * t427;
t118 = -t362 * t427 + t363 * t649;
t63 = Icges(6,5) * t118 + Icges(6,6) * t117 + Icges(6,3) * t280;
t65 = Icges(6,4) * t118 + Icges(6,2) * t117 + Icges(6,6) * t280;
t67 = Icges(6,1) * t118 + Icges(6,4) * t117 + Icges(6,5) * t280;
t8 = (qJD(3) * t443 + t63) * t382 + (qJD(3) * t119 - t362 * t65 + t363 * t67 + (-t122 * t363 - t125 * t362) * qJD(5)) * t384;
t623 = t8 * t305;
t426 = t385 * t352;
t115 = -t362 * t650 + t363 * t426;
t116 = t362 * t426 + t363 * t650;
t62 = Icges(6,5) * t116 + Icges(6,6) * t115 - Icges(6,3) * t677;
t64 = Icges(6,4) * t116 + Icges(6,2) * t115 - Icges(6,6) * t677;
t66 = Icges(6,1) * t116 + Icges(6,4) * t115 - Icges(6,5) * t677;
t9 = (qJD(3) * t441 + t62) * t382 + (qJD(3) * t121 - t362 * t64 + t363 * t66 + (-t124 * t363 + t126 * t362) * qJD(5)) * t384;
t622 = t9 * t304;
t258 = (-Icges(6,5) * t362 - Icges(6,6) * t363) * t384;
t131 = qJD(3) * t221 + qJD(5) * t258;
t223 = Icges(6,6) * t384 - t382 * t448;
t132 = qJD(3) * t223 + qJD(5) * t259;
t225 = Icges(6,5) * t384 - t382 * t451;
t133 = qJD(3) * t225 + qJD(5) * t260;
t27 = (qJD(3) * t435 + t131) * t382 + (qJD(3) * t222 - t132 * t362 + t133 * t363 + (-t224 * t363 - t226 * t362) * qJD(5)) * t384;
t493 = t384 * t517;
t81 = t222 * t382 - t384 * t435;
t618 = t27 * t352 + t81 * t493;
t615 = rSges(3,3) * t385;
t613 = rSges(5,3) * t384;
t611 = rSges(6,3) * t384;
t553 = t249 * rSges(6,1) + t248 * rSges(6,2);
t128 = -rSges(6,3) * t572 + t553;
t228 = -t382 * t460 + t611;
t263 = (-rSges(6,1) * t362 - rSges(6,2) * t363) * t384;
t135 = qJD(3) * t228 + qJD(5) * t263;
t213 = -t384 * t570 + t661;
t193 = t213 * qJD(3);
t509 = pkin(3) * t678 + qJ(4) * t501;
t525 = qJD(4) * t383;
t168 = (-qJ(4) * t530 - t525) * t384 + t509;
t386 = qJD(1) ^ 2;
t519 = qJD(1) * qJD(2);
t532 = qJD(1) * t383;
t546 = qJD(1) * (-pkin(1) * t532 + t540) + t383 * t519;
t428 = -t386 * t624 + t546;
t518 = qJD(1) * qJD(3);
t495 = t327 * t518;
t408 = qJD(1) * t168 + t383 * t495 + t428;
t322 = qJ(4) * t384 - t625;
t364 = qJD(4) * t382;
t270 = qJD(3) * t322 + t364;
t480 = -t270 - t364;
t514 = t118 * rSges(6,1) + t117 * rSges(6,2) + rSges(6,3) * t501;
t69 = -rSges(6,3) * t503 + t514;
t511 = t355 * t678 + t381 * t503;
t515 = pkin(4) * t573;
t83 = -t381 * t501 + (t348 - t515) * qJD(1) - t509 + t511;
t10 = -t135 * t305 - t219 * t229 + t352 * t69 + (t83 - t498) * qJD(1) + (-t644 * t532 + t128 * t522 + (-t193 + t480) * t385) * qJD(3) + t408;
t610 = t10 * t385;
t346 = t382 * t525;
t358 = t385 * t519;
t477 = -t376 * t386 + t358;
t407 = qJD(3) * t346 + t270 * t528 + t385 * t495 + t477;
t353 = pkin(3) * t577;
t508 = pkin(3) * t499 + qJ(4) * t677;
t524 = qJD(4) * t385;
t167 = qJD(1) * t353 + t384 * t524 - t508;
t367 = qJD(2) * t385;
t271 = qJD(1) * t329 - t367;
t557 = -t167 - t271;
t462 = rSges(6,1) * t116 + rSges(6,2) * t115;
t68 = -rSges(6,3) * t677 + t462;
t347 = pkin(4) * t579;
t578 = t381 * t384;
t580 = t355 * t384;
t82 = (t381 * t382 - t580) * t526 + (t347 + (t578 - t661) * t383) * qJD(1) + t508;
t11 = t135 * t304 + t220 * t229 - t352 * t68 + (-t130 * t522 + t193 * t383) * qJD(3) + (t557 - t82 + t656) * qJD(1) + t407;
t609 = t11 * t383;
t246 = -t382 * t463 + t613;
t227 = t246 * qJD(3);
t512 = t198 * rSges(5,1) + t197 * rSges(5,2) + rSges(5,3) * t501;
t91 = -rSges(5,3) * t503 + t512;
t35 = (t91 - t498) * qJD(1) + (t247 * t532 + (-t227 + t480) * t385) * qJD(3) + t408;
t607 = t35 * t385;
t465 = rSges(5,1) * t196 + rSges(5,2) * t195;
t90 = -rSges(5,3) * t677 + t465;
t36 = t227 * t528 + (t557 - t90 + t657) * qJD(1) + t407;
t605 = t36 * t383;
t507 = rSges(4,1) * t678 + rSges(4,2) * t503;
t529 = qJD(3) * t382;
t178 = (-rSges(4,2) * t529 - rSges(4,3) * qJD(1)) * t383 + t507;
t328 = rSges(4,1) * t384 - rSges(4,2) * t382;
t297 = t328 * t528;
t466 = rSges(4,1) * t382 + rSges(4,2) * t384;
t310 = t466 * qJD(3);
t75 = t310 * t526 + (t178 + t297) * qJD(1) + t428;
t604 = t385 * t75;
t47 = t119 * t382 - t384 * t443;
t603 = t47 * t219;
t602 = t48 * t220;
t294 = t328 * t385;
t177 = -qJD(3) * t294 + (t383 * t466 + t375) * qJD(1);
t502 = t328 * t526;
t76 = -t310 * t528 + (-t177 - t271 + t502) * qJD(1) + t477;
t601 = t76 * t383;
t600 = -rSges(5,3) - qJ(4);
t599 = -t168 - t83;
t598 = -t168 - t91;
t373 = t383 * rSges(4,3);
t265 = t466 * t385 - t373;
t112 = t297 + t366 + (t265 + t484) * qJD(1);
t583 = t112 * t385;
t575 = t383 * t327;
t563 = t135 + t193;
t562 = qJD(3) * t365 + t167 * t526;
t545 = t275 * rSges(5,1) + t274 * rSges(5,2);
t162 = -rSges(5,3) * t572 + t545;
t289 = -qJ(4) * t572 + t353;
t561 = -t162 - t289;
t464 = rSges(5,1) * t277 - rSges(5,2) * t276;
t164 = rSges(5,3) * t571 - t464;
t560 = -t164 + t292;
t510 = t355 * t577 + t381 * t572 + t347;
t165 = -t289 + t510;
t559 = -t165 - t289;
t543 = t355 * t576 + t381 * t571;
t166 = t292 + t515 - t543;
t558 = -t166 + t292;
t556 = -t644 + t229;
t552 = t383 * t270 + t327 * t530;
t293 = t327 * t385;
t547 = -t293 * t526 + t365;
t544 = -t327 * t526 - t367;
t539 = rSges(3,2) * t532 + rSges(3,3) * t530;
t538 = -t292 * t526 + t364;
t314 = qJD(1) * t324;
t537 = t366 - t314;
t534 = qJD(1) * t447;
t516 = -rSges(4,3) + t638;
t513 = -t128 + t559;
t491 = -t531 / 0.2e1;
t489 = -t528 / 0.2e1;
t488 = t528 / 0.2e1;
t487 = t527 / 0.2e1;
t486 = -t526 / 0.2e1;
t476 = t289 + t652;
t475 = t292 + t484;
t471 = qJD(5) * t487;
t470 = -t613 + t625;
t469 = qJD(1) * t575 - t382 * t524;
t468 = qJD(1) * t293 + t322 * t528 + t346;
t43 = -t119 * t572 + t617;
t459 = t383 * t44 + t385 * t43;
t458 = t383 * t43 - t385 * t44;
t46 = t121 * t571 - t442;
t457 = t383 * t46 + t385 * t45;
t456 = t383 * t45 - t385 * t46;
t455 = t383 * t48 + t385 * t47;
t454 = t383 * t47 - t385 * t48;
t113 = qJD(1) * t681 - t367 - t502;
t444 = t112 * t383 - t113 * t385;
t440 = t128 * t385 + t130 * t383;
t136 = (-t264 * t383 - t265 * t385) * qJD(3);
t406 = t119 * t305 + t121 * t304 + t222 * t352;
t405 = (Icges(6,5) * t248 - Icges(6,6) * t249) * t305 + (Icges(6,5) * t250 + Icges(6,6) * t251) * t304 + t258 * t352;
t399 = -qJD(1) * t431 - qJD(3) * t581 + t536;
t398 = qJD(1) * t433 + qJD(3) * t283 + t535;
t34 = -t128 * t304 + t130 * t305 + (t166 * t385 + t383 * t559) * qJD(3) + t538;
t39 = (-t166 + t475) * qJD(1) + t687;
t40 = t128 * t352 - t229 * t305 - t656 + (t165 + t476) * qJD(1) + t544;
t389 = t34 * t440 + (-t383 * t40 - t385 * t39) * t229;
t387 = qJD(3) * t643 + t651;
t325 = rSges(3,2) * t383 + t615;
t299 = t327 * t532;
t291 = t328 * t383;
t281 = (t385 ^ 2 + t641) * t527;
t267 = t385 * t292;
t218 = t247 * t528;
t208 = t247 * t385;
t207 = t247 * t383;
t206 = t244 * t385;
t205 = t244 * t383;
t204 = t242 * t385;
t203 = t242 * t383;
t200 = qJD(1) * t653 - t367;
t199 = t366 + (-t324 + t325) * qJD(1);
t190 = t229 * t385;
t189 = t229 * t383;
t188 = t226 * t385;
t187 = t226 * t383;
t186 = t224 * t385;
t185 = t224 * t383;
t180 = t644 * t385;
t179 = t644 * t383;
t176 = t358 + (-qJD(1) * t481 - t271) * qJD(1);
t175 = qJD(1) * t539 + t546;
t151 = rSges(6,1) * t250 + rSges(6,2) * t251;
t150 = rSges(6,1) * t248 - rSges(6,2) * t249;
t143 = t385 * t167;
t71 = -t657 + (t162 + t476) * qJD(1) + t544;
t70 = t218 + (-t164 + t475) * qJD(1) + t423;
t55 = (t164 * t385 + t383 * t561) * qJD(3) + t538;
t23 = (t385 * t90 + t598 * t383 + (t383 * t560 + t385 * t561) * qJD(1)) * qJD(3) + t562;
t22 = t117 * t224 + t118 * t226 - t131 * t572 + t132 * t248 + t133 * t249 + t222 * t280;
t21 = t115 * t224 + t116 * t226 + t131 * t571 + t132 * t250 - t133 * t251 - t222 * t677;
t18 = t304 * t48 + t305 * t47 + t352 * t81;
t13 = t304 * t46 - t694;
t12 = t305 * t43 + t662;
t7 = t117 * t124 - t118 * t126 + t121 * t280 + t248 * t64 + t249 * t66 - t572 * t62;
t6 = t117 * t122 + t118 * t125 + t119 * t280 + t248 * t65 + t249 * t67 - t572 * t63;
t5 = t115 * t124 - t116 * t126 - t121 * t677 + t250 * t64 - t251 * t66 + t571 * t62;
t4 = t115 * t122 + t116 * t125 - t119 * t677 + t250 * t65 - t251 * t67 + t571 * t63;
t3 = -t128 * t220 + t130 * t219 - t304 * t69 + t305 * t68 + (t385 * t82 + t599 * t383 + (t383 * t558 + t385 * t559) * qJD(1)) * qJD(3) + t562;
t2 = t219 * t43 + t22 * t352 + t220 * t44 + t304 * t7 + t305 * t6 + t493 * t73;
t1 = t21 * t352 + t219 * t45 + t220 * t46 + t304 * t5 + t305 * t4 + t493 * t74;
t14 = [t622 / 0.2e1 + t623 / 0.2e1 + t618 + t602 / 0.2e1 + t603 / 0.2e1 + t22 * t632 + ((t46 + t695) * t305 + t662) * t635 + t74 * t636 + t73 * t637 + (t21 + t12) * t634 + ((-t43 + t695) * t304 + t13 + t694) * t633 + (((-t655 + t95 + t700) * t385 + (-t384 * t439 - t94 + t93 + (t253 - t433) * t385 - t689) * t383) * qJD(3) + t686) * t489 + ((-t216 * t379 + t217 * t380 - t309) * t384 + (t308 + t215) * t382 + (t240 * t384 + t382 * t434 - t430) * qJD(3)) * qJD(1) + ((-t383 * t611 + t510 + t553 + t652) * t10 + (t369 + t543 + t690 - t130) * t11 + (t367 - t462 + (-t365 + (t580 + (rSges(6,3) - t381) * t382) * qJD(3)) * t385 + (t479 * t385 + (-t355 * t382 - qJ(2) - t578 + t611) * t383) * qJD(1)) * t39 + (t314 + t39 + t511 + t514 + (-t381 * t529 - t365) * t383 + (-t385 * t611 + t166 + t624 + t690) * qJD(1) - t687 + t688) * t40) * m(6) + (t36 * (-t348 + t369 + t464) + t70 * (t367 - t465 + t508) + t35 * (t353 + t652 + t545) + (qJD(1) * t164 - t218 - t296 + t509 + t512 - t537 + t688 + t70) * t71 + (t36 * t470 + t70 * (rSges(5,3) * t529 - t365) + (t384 * t600 * t71 + t638 * t70) * qJD(1)) * t385 + (-t71 * (-pkin(6) * qJD(1) - t365) + t36 * t638 + (-t71 * qJD(4) + t35 * t600) * t384 + (t70 * (-qJ(2) - t470) + t71 * t638) * qJD(1)) * t383) * m(5) + (t76 * (-t373 + t484) + t112 * t367 + t75 * t681 + t113 * (-rSges(4,2) * t501 + t507 + t540) + (qJD(3) * t112 * t328 + t466 * t76) * t385 + (t516 * t583 + (t112 * (-qJ(2) - t466) + t113 * t516) * t383) * qJD(1) - (-t112 + t297 + (t265 - t624) * qJD(1) + t537) * t113) * m(4) + (t176 * (t383 * t626 + t369 + t615) + t199 * t367 + t175 * t653 + t200 * (t539 + t540) + (t199 * t626 * t385 + (t199 * (-rSges(3,3) - qJ(2)) - t200 * pkin(1)) * t383) * qJD(1) - (qJD(1) * t325 - t199 + t537) * t200) * m(3) + ((t437 * t383 + t253 * t641 + (-t51 + (t253 + t433) * t385 + t689 + t693) * t385) * qJD(3) + t669 - t676) * t486 + (t666 + t667 + t670) * t488 + ((t79 + t663) * qJD(1) + t665 + t668) * t526 / 0.2e1 - (t142 * t385 + (-t664 + t696) * t383) * t518 / 0.2e1; 0.2e1 * (-t610 / 0.2e1 + t609 / 0.2e1) * m(6) + 0.2e1 * (-t607 / 0.2e1 + t605 / 0.2e1) * m(5) + 0.2e1 * (t601 / 0.2e1 - t604 / 0.2e1) * m(4) + 0.2e1 * (-t175 * t385 / 0.2e1 + t176 * t628) * m(3); -t18 * t522 / 0.2e1 + t455 * t471 - t383 * t12 * t523 / 0.2e1 + t13 * t521 * t629 + (-qJD(1) * t454 + t383 * t9 + t385 * t8) * t630 + (((-t185 * t362 + t187 * t363 + t119) * t305 + (t186 * t362 - t188 * t363 + t121) * t304 + (-t223 * t362 + t225 * t363 + t222) * t352 + t81 * qJD(5)) * t384 + (qJD(5) * t454 + t388) * t382) * t631 + (-qJD(1) * t458 + t383 * t7 + t385 * t6) * t632 + ((t185 * t248 + t187 * t249) * t305 + (-t186 * t248 - t188 * t249) * t304 + (t223 * t248 + t225 * t249) * t352 + (t384 * t73 - t44 * t576) * qJD(5) + ((qJD(5) * t43 + t406) * t382 - t660) * t383) * t633 + (-qJD(1) * t456 + t383 * t5 + t385 * t4) * t634 + ((t185 * t250 - t187 * t251) * t305 + (-t186 * t250 + t188 * t251) * t304 + (t223 * t250 - t225 * t251) * t352 + (t384 * t74 + t45 * t577) * qJD(5) + ((-qJD(5) * t46 - t406) * t382 + t660) * t385) * t635 + t457 * t636 + t459 * t637 - (((-t241 * t379 + t243 * t380 + t240) * qJD(1) + ((-t203 * t379 + t205 * t380 + t153) * t385 + (t204 * t379 - t206 * t380 + t155) * t383) * qJD(3)) * t384 + t387 * t382 + (-t382 * t542 - t384 * t541) * qJD(1) + (t382 * t400 + t384 * t401) * qJD(3)) * qJD(1) / 0.2e1 + (t668 * t385 + t667 * t383 + (t383 * t664 + t385 * t663) * qJD(1)) * qJD(1) / 0.2e1 + ((-t528 * t581 - t534) * t383 + (-t204 * t276 + t206 * t277) * t528 + (t241 * t276 - t243 * t277) * qJD(1) + (t387 * t384 + (t203 * t276 - t205 * t277 + t383 * t283 + t675) * qJD(3) - t673) * t385) * t489 + ((t283 * t526 - t534) * t385 + (t203 * t274 + t205 * t275) * t526 + (t241 * t274 + t243 * t275) * qJD(1) + (-t651 * t384 + (-t274 * t204 - t275 * t206 - t384 * t643 - t385 * t581 - t675) * qJD(3) + t673) * t383) * t486 + (-t39 * (-qJD(1) * t180 + t190 * t352 + t228 * t304 + t468) - t40 * (-qJD(1) * t179 + t189 * t352 - t228 * t305 + t469) - t34 * (-t189 * t304 - t190 * t305 + t547) - ((t128 * t40 - t130 * t39) * t384 + t389 * t382) * qJD(5) - ((t40 * (-t213 - t322) + t34 * t180) * t385 + (t39 * t213 + t34 * (t179 - t575)) * t383) * qJD(3) + t11 * t575 + t39 * t552 + t40 * t299 - t3 * t267 + t34 * t143 + (t11 * t556 + t39 * t563 + t3 * t513 + t34 * (-t69 + t599) + (t40 * t556 + t34 * (-t130 + t558)) * qJD(1)) * t383 + (t10 * (-t327 - t556) + t40 * (-t270 - t563) + t3 * (t130 + t166) + t34 * (t68 + t82) + (t34 * t513 + t39 * t556) * qJD(1)) * t385) * m(6) + (-t70 * (qJD(1) * t208 + t468) - t71 * (qJD(1) * t207 + t469) - t55 * t547 - ((t71 * (-t246 - t322) - t55 * t208) * t385 + (t70 * t246 + t55 * (-t207 - t575)) * t383) * qJD(3) + t36 * t575 + t70 * t552 + t71 * t299 - t23 * t267 + t55 * t143 + (t36 * t247 + t70 * t227 + t23 * t561 + t55 * t598 + (t71 * t247 + t55 * t560) * qJD(1)) * t383 + (t35 * (-t247 - t327) + t71 * (-t227 - t270) + t23 * t164 + t55 * t90 + (t70 * t247 + t55 * t561) * qJD(1)) * t385) * m(5) + (-(t112 * t294 + t113 * t291) * qJD(1) - (t136 * (-t291 * t383 - t294 * t385) - t444 * t466) * qJD(3) + 0.2e1 * t136 * (t177 * t385 - t178 * t383 + (-t264 * t385 + t265 * t383) * qJD(1)) - t444 * t310 + (t601 - t604 + (t113 * t383 + t583) * qJD(1)) * t328) * m(4) + (t666 * qJD(1) + t1 + ((-t153 * t677 + t156 * t195 + t159 * t196 + t276 * t87 - t277 * t89 + t385 * t646 + t571 * t85) * t385 + (-t155 * t677 + t158 * t195 + t161 * t196 + t276 * t86 - t277 * t88 + t399 * t383 + t571 * t84 + (t398 - t648) * t385) * t383 + (-t383 * t698 + t385 * t697) * qJD(1)) * t671) * t628 + (t665 * qJD(1) + t2 + ((t153 * t280 + t156 * t197 + t159 * t198 + t274 * t87 + t275 * t89 + t398 * t385 - t572 * t85) * t385 + (t155 * t280 + t158 * t197 + t161 * t198 + t274 * t86 + t275 * t88 + t648 * t383 - t572 * t84 + (t399 - t646) * t385) * t383 + (-t383 * t699 + t693 * t385) * qJD(1)) * t671) * t385 / 0.2e1 - (t12 + t670 + t679) * t532 / 0.2e1 + (t13 + t669 + t680) * t530 / 0.2e1; -m(5) * (t280 * t70 + t281 * t55 - t677 * t71) - m(6) * (t280 * t39 + t281 * t34 - t40 * t677) + 0.2e1 * ((-t526 * t71 + t528 * t70 + t23) * t640 + (t39 * t528 - t40 * t526 + t3) * t639) * t382 + 0.2e1 * ((qJD(3) * t55 - t530 * t70 - t532 * t71 - t605 + t607) * t640 + (qJD(3) * t34 - t39 * t530 - t40 * t532 - t609 + t610) * t639) * t384; -t2 * t572 / 0.2e1 + (t382 * t73 - t384 * t458) * t637 + ((qJD(3) * t458 + t22) * t382 + (-qJD(1) * t459 + qJD(3) * t73 - t383 * t6 + t385 * t7) * t384) * t632 + t1 * t571 / 0.2e1 + (t382 * t74 - t384 * t456) * t636 + ((qJD(3) * t456 + t21) * t382 + (-qJD(1) * t457 + qJD(3) * t74 - t383 * t4 + t385 * t5) * t384) * t634 + t18 * t487 + (t602 + t603 + t618 + t622 + t623) * t629 + (t382 * t81 - t384 * t454) * t471 + ((qJD(3) * t454 + t27) * t382 + (-qJD(1) * t455 + qJD(3) * t81 - t383 * t8 + t385 * t9) * t384) * t630 + (t248 * t395 - t249 * t642 - t405 * t572) * t633 + (t395 * t250 + t251 * t642 + t405 * t571) * t635 + (t405 * t382 + (-t362 * t395 - t363 * t642) * t384) * t631 + (t382 * t486 + t383 * t491) * t13 + (t382 * t488 + t385 * t491) * t12 + ((qJD(3) * t389 + t10 * t128 - t11 * t130 - t39 * t68 + t40 * t69) * t382 + (t39 * (-qJD(3) * t130 + t135 * t385) + t40 * (qJD(3) * t128 + t135 * t383) - t3 * t440 + t34 * (t128 * t532 - t130 * t530 - t383 * t68 - t385 * t69) + (t10 * t383 + t11 * t385 + (-t383 * t39 + t385 * t40) * qJD(1)) * t229) * t384 - t39 * (-t151 * t352 + t263 * t304) - t40 * (t150 * t352 - t263 * t305) - t34 * (-t150 * t304 + t151 * t305)) * m(6);];
tauc = t14(:);
