% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:54
% EndTime: 2019-12-31 18:19:40
% DurationCPUTime: 39.06s
% Computational Cost: add. (22637->851), mult. (21957->1140), div. (0->0), fcn. (19897->10), ass. (0->434)
t698 = Icges(4,3) + Icges(5,3);
t345 = qJ(3) + pkin(9);
t339 = sin(t345);
t341 = cos(t345);
t260 = Icges(5,5) * t341 - Icges(5,6) * t339;
t349 = sin(qJ(3));
t352 = cos(qJ(3));
t296 = Icges(4,5) * t352 - Icges(4,6) * t349;
t694 = t260 + t296;
t346 = qJ(1) + pkin(8);
t342 = cos(t346);
t697 = t698 * t342;
t340 = sin(t346);
t555 = t340 * t352;
t557 = t340 * t349;
t559 = t340 * t341;
t561 = t339 * t340;
t678 = -Icges(4,5) * t555 - Icges(5,5) * t559 + Icges(4,6) * t557 + Icges(5,6) * t561 + t697;
t677 = t698 * t340 + t694 * t342;
t576 = Icges(5,6) * t342;
t192 = Icges(5,4) * t559 - Icges(5,2) * t561 - t576;
t577 = Icges(4,6) * t342;
t204 = Icges(4,4) * t555 - Icges(4,2) * t557 - t577;
t696 = t192 * t339 + t204 * t349;
t586 = Icges(5,4) * t339;
t264 = Icges(5,1) * t341 - t586;
t195 = Icges(5,5) * t340 + t264 * t342;
t587 = Icges(4,4) * t349;
t300 = Icges(4,1) * t352 - t587;
t209 = Icges(4,5) * t340 + t300 * t342;
t695 = -t195 * t559 - t209 * t555;
t261 = Icges(5,2) * t341 + t586;
t327 = Icges(5,4) * t341;
t263 = Icges(5,1) * t339 + t327;
t297 = Icges(4,2) * t352 + t587;
t343 = Icges(4,4) * t352;
t299 = Icges(4,1) * t349 + t343;
t692 = t261 * t339 - t263 * t341 + t297 * t349 - t299 * t352;
t284 = Icges(5,4) * t561;
t581 = Icges(5,5) * t342;
t194 = Icges(5,1) * t559 - t284 - t581;
t308 = Icges(4,4) * t557;
t582 = Icges(4,5) * t342;
t208 = Icges(4,1) * t555 - t308 - t582;
t681 = -t194 * t341 - t208 * t352 + t696;
t259 = Icges(5,5) * t339 + Icges(5,6) * t341;
t295 = Icges(4,5) * t349 + Icges(4,6) * t352;
t693 = t295 + t259;
t691 = t677 * t342 + t695;
t550 = t342 * t352;
t554 = t341 * t342;
t660 = t195 * t554 + t209 * t550 + t677 * t340;
t690 = -t194 * t554 - t208 * t550 + t678 * t340;
t689 = -t681 * t340 + t678 * t342;
t428 = -Icges(5,2) * t339 + t327;
t193 = Icges(5,6) * t340 + t342 * t428;
t429 = -Icges(4,2) * t349 + t343;
t205 = Icges(4,6) * t340 + t342 * t429;
t667 = -t193 * t561 - t205 * t557 - t691;
t552 = t342 * t349;
t560 = t339 * t342;
t666 = -t192 * t560 - t204 * t552 - t690;
t665 = -t193 * t560 - t205 * t552 + t660;
t562 = t295 * t342;
t564 = t259 * t342;
t688 = -t692 * t340 - t562 - t564;
t563 = t295 * t340;
t565 = t259 * t340;
t687 = -t692 * t342 + t563 + t565;
t347 = -qJ(4) - pkin(6);
t314 = t342 * t347;
t615 = pkin(3) * t352;
t336 = pkin(2) + t615;
t525 = -t340 * t336 - t314;
t350 = sin(qJ(1));
t617 = pkin(1) * t350;
t686 = t525 - t617;
t685 = t193 * t339 + t205 * t349;
t253 = t428 * qJD(3);
t254 = t264 * qJD(3);
t276 = t429 * qJD(3);
t277 = t300 * qJD(3);
t684 = -t253 * t339 + t254 * t341 - t276 * t349 + t277 * t352 + (-t261 * t341 - t263 * t339 - t297 * t352 - t299 * t349) * qJD(3) + t693 * qJD(1);
t498 = rSges(5,1) * t559;
t683 = -t498 + t686;
t682 = t195 * t341 + t209 * t352 - t685;
t680 = t692 * qJD(1) + t694 * qJD(3);
t679 = t687 * qJD(1);
t676 = (t665 * t340 - t666 * t342) * qJD(3);
t675 = (t667 * t340 - t689 * t342) * qJD(3);
t674 = t688 * qJD(1);
t351 = cos(qJ(5));
t551 = t342 * t351;
t348 = sin(qJ(5));
t558 = t340 * t348;
t233 = t341 * t558 + t551;
t553 = t342 * t348;
t556 = t340 * t351;
t234 = t341 * t556 - t553;
t117 = Icges(6,5) * t234 - Icges(6,6) * t233 + Icges(6,3) * t561;
t221 = Icges(6,4) * t234;
t120 = -Icges(6,2) * t233 + Icges(6,6) * t561 + t221;
t220 = Icges(6,4) * t233;
t124 = -Icges(6,1) * t234 - Icges(6,5) * t561 + t220;
t658 = t120 * t348 + t124 * t351;
t47 = -t117 * t341 - t339 * t658;
t673 = t674 + t675;
t672 = t676 + t679;
t394 = qJD(3) * t261;
t112 = qJD(1) * t193 - t340 * t394;
t396 = qJD(3) * t263;
t114 = qJD(1) * t195 - t340 * t396;
t395 = qJD(3) * t297;
t140 = qJD(1) * t205 - t340 * t395;
t397 = qJD(3) * t299;
t143 = qJD(1) * t209 - t340 * t397;
t671 = t681 * qJD(3) - t112 * t341 - t114 * t339 - t140 * t352 - t143 * t349;
t111 = -t342 * t394 + (-t340 * t428 + t576) * qJD(1);
t113 = -t342 * t396 + (-t264 * t340 + t581) * qJD(1);
t139 = -t342 * t395 + (-t340 * t429 + t577) * qJD(1);
t142 = -t342 * t397 + (-t300 * t340 + t582) * qJD(1);
t670 = t682 * qJD(3) + t111 * t341 + t113 * t339 + t139 * t352 + t142 * t349;
t669 = t680 * t340 + t684 * t342;
t668 = t684 * t340 - t680 * t342;
t664 = t192 * t341 + t194 * t339 + t204 * t352 + t208 * t349;
t663 = t193 * t341 + t195 * t339 + t205 * t352 + t209 * t349;
t662 = t693 * qJD(3);
t661 = t678 + t685;
t659 = t677 * qJD(1);
t507 = qJD(5) * t339;
t511 = qJD(3) * t340;
t246 = t342 * t507 + t511;
t510 = qJD(3) * t342;
t247 = -t340 * t507 + t510;
t506 = qJD(5) * t341;
t311 = qJD(1) - t506;
t35 = t117 * t561 - t120 * t233 - t124 * t234;
t235 = -t341 * t553 + t556;
t236 = t341 * t551 + t558;
t119 = Icges(6,5) * t236 + Icges(6,6) * t235 + Icges(6,3) * t560;
t585 = Icges(6,4) * t236;
t122 = Icges(6,2) * t235 + Icges(6,6) * t560 + t585;
t222 = Icges(6,4) * t235;
t125 = Icges(6,1) * t236 + Icges(6,5) * t560 + t222;
t36 = t119 * t561 - t233 * t122 + t234 * t125;
t426 = Icges(6,5) * t351 - Icges(6,6) * t348;
t198 = -Icges(6,3) * t341 + t339 * t426;
t583 = Icges(6,4) * t351;
t427 = -Icges(6,2) * t348 + t583;
t202 = -Icges(6,6) * t341 + t339 * t427;
t584 = Icges(6,4) * t348;
t430 = Icges(6,1) * t351 - t584;
t206 = -Icges(6,5) * t341 + t339 * t430;
t64 = t198 * t561 - t202 * t233 + t206 * t234;
t10 = t246 * t36 - t247 * t35 + t311 * t64;
t37 = t117 * t560 + t235 * t120 - t124 * t236;
t38 = t119 * t560 + t235 * t122 + t236 * t125;
t65 = t198 * t560 + t202 * t235 + t206 * t236;
t11 = t246 * t38 - t247 * t37 + t65 * t311;
t355 = qJD(1) ^ 2;
t657 = rSges(4,2) * t349;
t653 = -qJD(3) * t663 - t111 * t339 + t113 * t341 - t139 * t349 + t142 * t352 + t659;
t652 = t678 * qJD(1) + t664 * qJD(3) + t112 * t339 - t114 * t341 + t140 * t349 - t143 * t352;
t613 = pkin(4) * t341;
t270 = pkin(7) * t339 + t613;
t228 = t270 * t340;
t334 = t342 * pkin(6);
t269 = pkin(2) * t340 - t334;
t186 = t269 + t525;
t470 = -t269 - t617;
t454 = t186 + t470;
t441 = rSges(6,1) * t234 - rSges(6,2) * t233;
t128 = rSges(6,3) * t561 + t441;
t440 = rSges(6,1) * t351 - rSges(6,2) * t348;
t210 = -rSges(6,3) * t341 + t339 * t440;
t319 = qJD(4) * t340;
t614 = pkin(4) * t339;
t268 = -pkin(7) * t341 + t614;
t616 = pkin(3) * t349;
t466 = -t268 - t616;
t446 = t466 * t342;
t648 = qJD(3) * t446 - t128 * t311 - t247 * t210 + t319;
t39 = (-t228 + t454) * qJD(1) + t648;
t130 = t236 * rSges(6,1) + t235 * rSges(6,2) + rSges(6,3) * t560;
t230 = pkin(4) * t554 + pkin(7) * t560;
t333 = t340 * pkin(6);
t271 = t342 * pkin(2) + t333;
t291 = t342 * t336;
t458 = -t340 * t347 + t291;
t187 = t458 - t271;
t353 = cos(qJ(1));
t344 = t353 * pkin(1);
t469 = t271 + t344;
t453 = t187 + t469;
t489 = t268 * t511;
t509 = qJD(3) * t349;
t293 = t340 * pkin(3) * t509;
t523 = qJD(4) * t342 + t293;
t40 = -t489 + t130 * t311 - t210 * t246 + (t230 + t453) * qJD(1) - t523;
t596 = t340 * t40;
t651 = t342 * t39 + t596;
t650 = t681 * qJD(1) - t662 * t340 + t659;
t649 = -t662 * t342 + (-t694 * t340 - t682 + t697) * qJD(1);
t647 = 0.2e1 * qJD(3);
t182 = qJD(1) * t186;
t258 = qJD(1) * t269;
t644 = t182 - t258;
t329 = t340 * rSges(4,3);
t213 = rSges(4,1) * t550 - rSges(4,2) * t552 + t329;
t643 = t213 + t469;
t642 = -rSges(5,2) * t561 - t342 * rSges(5,3);
t463 = t342 * rSges(3,1) - rSges(3,2) * t340;
t641 = t344 + t463;
t632 = t340 * (-t261 * t342 + t195) - t342 * (-Icges(5,2) * t559 + t194 - t284);
t530 = -Icges(4,2) * t555 + t208 - t308;
t533 = t299 * t340 + t204;
t631 = -t349 * t530 - t352 * t533;
t199 = Icges(6,3) * t339 + t341 * t426;
t420 = -t202 * t348 + t206 * t351;
t424 = -t122 * t348 + t125 * t351;
t630 = t246 * (-t198 * t342 - t424) - t247 * (-t198 * t340 + t658) + t311 * (t199 - t420);
t240 = (-Icges(6,2) * t351 - t584) * t339;
t629 = t246 * (-Icges(6,2) * t236 + t125 + t222) - t247 * (-Icges(6,2) * t234 - t124 - t220) + t311 * (t206 + t240);
t503 = qJD(3) * qJD(5);
t479 = t341 * t503;
t180 = qJD(1) * t246 + t340 * t479;
t628 = t180 / 0.2e1;
t181 = qJD(1) * t247 + t342 * t479;
t627 = t181 / 0.2e1;
t626 = -t246 / 0.2e1;
t625 = t246 / 0.2e1;
t624 = -t247 / 0.2e1;
t623 = t247 / 0.2e1;
t622 = -t311 / 0.2e1;
t621 = t311 / 0.2e1;
t620 = t340 / 0.2e1;
t619 = -t342 / 0.2e1;
t618 = -rSges(6,3) - pkin(7);
t512 = qJD(3) * t339;
t379 = t311 * t351 + t348 * t512;
t514 = qJD(1) * t341;
t457 = -qJD(5) + t514;
t103 = t340 * t379 - t457 * t553;
t378 = t311 * t348 - t351 * t512;
t104 = t340 * t378 + t457 * t551;
t513 = qJD(1) * t342;
t388 = t339 * t513 + t341 * t511;
t54 = Icges(6,5) * t104 + Icges(6,6) * t103 + Icges(6,3) * t388;
t56 = Icges(6,4) * t104 + Icges(6,2) * t103 + Icges(6,6) * t388;
t58 = Icges(6,1) * t104 + Icges(6,4) * t103 + Icges(6,5) * t388;
t8 = (-qJD(3) * t658 - t54) * t341 + (qJD(3) * t117 - t348 * t56 + t351 * t58 + (-t120 * t351 + t124 * t348) * qJD(5)) * t339;
t612 = t8 * t247;
t101 = t342 * t379 + t457 * t558;
t102 = t342 * t378 - t457 * t556;
t486 = t341 * t510;
t515 = qJD(1) * t340;
t492 = t339 * t515;
t387 = t486 - t492;
t53 = Icges(6,5) * t102 + Icges(6,6) * t101 + Icges(6,3) * t387;
t55 = Icges(6,4) * t102 + Icges(6,2) * t101 + Icges(6,6) * t387;
t57 = Icges(6,1) * t102 + Icges(6,4) * t101 + Icges(6,5) * t387;
t9 = (qJD(3) * t424 - t53) * t341 + (qJD(3) * t119 - t348 * t55 + t351 * t57 + (-t122 * t351 - t125 * t348) * qJD(5)) * t339;
t611 = t9 * t246;
t609 = qJD(1) / 0.2e1;
t608 = pkin(2) - t336;
t237 = (-Icges(6,5) * t348 - Icges(6,6) * t351) * t339;
t135 = qJD(3) * t199 + qJD(5) * t237;
t203 = Icges(6,6) * t339 + t341 * t427;
t138 = qJD(3) * t203 + qJD(5) * t240;
t207 = Icges(6,5) * t339 + t341 * t430;
t243 = (-Icges(6,1) * t348 - t583) * t339;
t141 = qJD(3) * t207 + qJD(5) * t243;
t23 = (qJD(3) * t420 - t135) * t341 + (qJD(3) * t198 - t138 * t348 + t141 * t351 + (-t202 * t351 - t206 * t348) * qJD(5)) * t339;
t480 = t339 * t503;
t77 = -t198 * t341 + t339 * t420;
t607 = t23 * t311 + t77 * t480;
t606 = rSges(4,1) * t352;
t604 = rSges(6,3) * t339;
t274 = pkin(7) * t486;
t487 = t339 * t510;
t389 = -t340 * t514 - t487;
t146 = pkin(4) * t389 - pkin(7) * t492 + t274;
t212 = t341 * t440 + t604;
t248 = (-rSges(6,1) * t348 - rSges(6,2) * t351) * t339;
t148 = qJD(3) * t212 + qJD(5) * t248;
t256 = t270 * qJD(3);
t318 = pkin(6) * t513;
t483 = t342 * t509;
t413 = -pkin(3) * t483 + t319;
t132 = -t318 + (t340 * t608 - t314) * qJD(1) + t413;
t500 = t355 * t617;
t455 = qJD(1) * (-pkin(2) * t515 + t318) - t500;
t504 = qJD(1) * qJD(4);
t403 = qJD(1) * t132 + t340 * t504 + t455;
t501 = qJD(3) ^ 2 * t615;
t496 = t102 * rSges(6,1) + t101 * rSges(6,2) + rSges(6,3) * t486;
t59 = -rSges(6,3) * t492 + t496;
t12 = -t340 * t501 + qJD(1) * t146 - t148 * t246 - t181 * t210 + t311 * t59 + (qJD(1) * t446 + t130 * t507 - t256 * t340) * qJD(3) + t403;
t602 = t12 * t342;
t147 = t388 * pkin(7) + (-t339 * t511 + t341 * t513) * pkin(4);
t499 = t355 * t344;
t414 = qJD(1) * t293 + t342 * t504 - t499;
t493 = t347 * t515 + t523;
t133 = (-t342 * t608 - t333) * qJD(1) - t493;
t257 = t271 * qJD(1);
t544 = -t133 - t257;
t442 = rSges(6,1) * t104 + rSges(6,2) * t103;
t60 = rSges(6,3) * t388 + t442;
t13 = -t342 * t501 - t148 * t247 + t180 * t210 - t311 * t60 + (-t128 * t507 - t256 * t342) * qJD(3) + (-t147 + t489 + t544) * qJD(1) + t414;
t601 = t13 * t340;
t302 = rSges(4,1) * t349 + rSges(4,2) * t352;
t250 = t302 * t342;
t488 = t302 * t511;
t98 = qJD(1) * t643 - t488;
t598 = t250 * t98;
t328 = t340 * rSges(5,3);
t520 = rSges(4,2) * t557 + t342 * rSges(4,3);
t211 = rSges(4,1) * t555 - t520;
t484 = t302 * t510;
t97 = -t484 + (-t211 + t470) * qJD(1);
t595 = t340 * t97;
t594 = t342 * t97;
t593 = t47 * t180;
t48 = -t119 * t341 + t339 * t424;
t592 = t48 * t181;
t265 = rSges(5,1) * t339 + rSges(5,2) * t341;
t196 = t498 + t642;
t468 = -t265 - t616;
t410 = t468 * t510;
t386 = t319 + t410;
t68 = (-t196 + t454) * qJD(1) + t386;
t591 = t68 * t265;
t546 = t128 + t228;
t545 = t130 + t230;
t540 = -t340 * t186 + t342 * t187;
t197 = rSges(5,1) * t554 - rSges(5,2) * t560 + t328;
t537 = -t187 - t197;
t536 = -t187 - t230;
t532 = -t299 * t342 - t205;
t529 = -t297 * t342 + t209;
t528 = -t261 + t264;
t527 = t263 + t428;
t526 = rSges(5,2) * t492 + rSges(5,3) * t513;
t524 = rSges(4,3) * t513 + t515 * t657;
t522 = -t297 + t300;
t521 = t299 + t429;
t517 = qJD(1) * t260;
t516 = qJD(1) * t296;
t508 = qJD(3) * t352;
t502 = pkin(3) * t552;
t497 = pkin(3) * t508;
t495 = t133 * t511 + (t132 - t182) * t510;
t494 = t342 * t132 + t340 * t133 - t186 * t513;
t490 = t265 * t511;
t482 = -t186 * t511 + t187 * t510 + qJD(2);
t481 = -pkin(2) - t606;
t477 = t513 / 0.2e1;
t476 = t512 / 0.2e1;
t475 = -t511 / 0.2e1;
t474 = t511 / 0.2e1;
t472 = t510 / 0.2e1;
t267 = rSges(5,1) * t341 - rSges(5,2) * t339;
t467 = -t267 - t615;
t456 = (-t340 ^ 2 - t342 ^ 2) * t616;
t452 = -t210 + t466;
t449 = qJD(5) * t476;
t255 = t267 * qJD(3);
t448 = -t255 - t497;
t69 = -t490 + (t197 + t453) * qJD(1) - t523;
t447 = t69 * t468;
t266 = rSges(3,1) * t340 + rSges(3,2) * t342;
t443 = t606 - t657;
t439 = t340 * t36 - t342 * t35;
t438 = t340 * t35 + t342 * t36;
t437 = t340 * t38 - t342 * t37;
t436 = t340 * t37 + t342 * t38;
t435 = t340 * t48 - t342 * t47;
t434 = t340 * t47 + t342 * t48;
t433 = -t340 * t98 - t594;
t423 = t128 * t342 - t130 * t340;
t417 = t211 * t340 + t213 * t342;
t412 = -t148 - t256 - t497;
t409 = -qJD(3) * t255 - t501;
t249 = t302 * t340;
t223 = t265 * t340;
t391 = -t270 - t604;
t385 = -t117 * t247 + t119 * t246 + t198 * t311;
t384 = (-Icges(6,5) * t233 - Icges(6,6) * t234) * t247 - (Icges(6,5) * t235 - Icges(6,6) * t236) * t246 - t237 * t311;
t383 = t192 * t342 - t193 * t340;
t382 = -t349 * t529 + t352 * t532;
t381 = t339 * t384;
t374 = (-t339 * t527 + t341 * t528) * qJD(1);
t373 = (-t349 * t521 + t352 * t522) * qJD(1);
t371 = (Icges(6,1) * t235 - t122 - t585) * t246 - (-Icges(6,1) * t233 - t120 - t221) * t247 + (-t202 + t243) * t311;
t149 = -rSges(4,2) * t342 * t508 + (-t352 * t515 - t483) * rSges(4,1) + t524;
t150 = -qJD(3) * t249 + (t342 * t443 + t329) * qJD(1);
t369 = t149 * t342 + t150 * t340 + (t211 * t342 - t213 * t340) * qJD(1);
t30 = t128 * t246 + t130 * t247 + (t228 * t340 + t230 * t342) * qJD(3) + t482;
t364 = t30 * t423 + (t340 * t39 - t342 * t40) * t210;
t357 = -t339 * t632 + t383 * t341;
t356 = t630 * t339;
t278 = t443 * qJD(3);
t229 = t268 * t342;
t227 = t268 * t340;
t224 = t265 * t342;
t170 = t210 * t342;
t169 = t210 * t340;
t168 = t206 * t342;
t167 = t206 * t340;
t166 = t202 * t342;
t165 = t202 * t340;
t158 = rSges(6,1) * t235 - rSges(6,2) * t236;
t157 = -rSges(6,1) * t233 - rSges(6,2) * t234;
t116 = -qJD(3) * t223 + (t267 * t342 + t328) * qJD(1);
t115 = rSges(5,1) * t389 - rSges(5,2) * t486 + t526;
t96 = qJD(3) * t417 + qJD(2);
t75 = -t499 - t278 * t510 + (-t150 - t257 + t488) * qJD(1);
t74 = -t278 * t511 + (t149 - t484) * qJD(1) + t455;
t52 = (t196 * t340 + t197 * t342) * qJD(3) + t482;
t49 = t369 * qJD(3);
t42 = t409 * t342 + (-t116 + t490 + t544) * qJD(1) + t414;
t41 = t409 * t340 + (t115 + t410) * qJD(1) + t403;
t17 = (t115 * t342 + t116 * t340 + (t196 * t342 + t340 * t537) * qJD(1)) * qJD(3) + t495;
t16 = t103 * t202 + t104 * t206 + t135 * t561 - t138 * t233 + t141 * t234 + t198 * t388;
t15 = t101 * t202 + t102 * t206 + t135 * t560 + t138 * t235 + t141 * t236 + t198 * t387;
t14 = t246 * t48 - t247 * t47 + t311 * t77;
t7 = t103 * t122 + t104 * t125 + t119 * t388 - t233 * t55 + t234 * t57 + t53 * t561;
t6 = t103 * t120 - t104 * t124 + t117 * t388 - t233 * t56 + t234 * t58 + t54 * t561;
t5 = t101 * t122 + t102 * t125 + t119 * t387 + t235 * t55 + t236 * t57 + t53 * t560;
t4 = t101 * t120 - t102 * t124 + t117 * t387 + t235 * t56 + t236 * t58 + t54 * t560;
t3 = t128 * t181 - t130 * t180 + t246 * t60 + t247 * t59 + (t146 * t342 + t147 * t340 + (t228 * t342 + t340 * t536) * qJD(1)) * qJD(3) + t495;
t2 = t16 * t311 + t180 * t35 + t181 * t36 + t246 * t7 - t247 * t6 + t480 * t64;
t1 = t15 * t311 + t180 * t37 + t181 * t38 + t246 * t5 - t247 * t4 + t480 * t65;
t18 = [t16 * t624 + t15 * t625 + m(3) * ((-t266 * t355 - t500) * t641 + (-t499 + (-0.2e1 * t463 - t344 + t641) * t355) * (-t266 - t617)) + t607 + t592 / 0.2e1 + t593 / 0.2e1 + t611 / 0.2e1 - t612 / 0.2e1 + t65 * t627 + t64 * t628 + ((t660 * t340 + ((t677 + t696) * t342 + t667 + t690 + t695) * t342) * qJD(3) + t679) * t472 + (t624 + t623) * t11 + (-t692 * qJD(3) + t253 * t341 + t254 * t339 + t276 * t352 + t277 * t349) * qJD(1) + (-(-t39 + (-t228 - t617) * qJD(1) + t644 + t648) * t40 + t13 * (-t441 + t686) + t39 * (-t442 + t493) + t12 * (t291 + t344 + t545) + t40 * (-pkin(4) * t487 + t274 + t413 + t496) + (t13 * t391 - t12 * t347 + t39 * (t341 * t618 + t614) * qJD(3)) * t340 + ((-t350 * t40 - t353 * t39) * pkin(1) + (t39 * (-t336 + t391) - t40 * t347) * t342 + (t339 * t618 - t336 - t613) * t596) * qJD(1)) * m(6) + (t42 * (-t642 + t683) + t41 * (t197 + t344 + t458) + (t340 * t591 + t342 * t447) * qJD(3) + (t493 + (-t328 - t344 + (-t267 - t336) * t342) * qJD(1)) * t68 + (t319 + t526 + t68 - t386 - t644 + (t196 + t617 + t683) * qJD(1)) * t69) * m(5) + (-(-t484 - t258 - t97 + (-t211 - t617) * qJD(1)) * t98 + t75 * (t340 * t481 + t334 + t520 - t617) + t74 * t643 + t98 * (t318 + t524) + (t302 * t595 - t598) * qJD(3) + ((-t350 * t98 - t353 * t97) * pkin(1) + (-pkin(2) - t443) * t594 + (t97 * (-rSges(4,3) - pkin(6)) + t98 * t481) * t340) * qJD(1)) * m(4) + (((t342 * t661 - t660 + t665) * t342 + (t340 * t661 + t666 + t691) * t340) * qJD(3) + t673 - t674) * t475 + (t669 + t670) * t474 - (t668 - t671 + t672) * t510 / 0.2e1 + ((t664 + t688) * t340 + (t663 + t687) * t342) * qJD(3) * t609; m(4) * t49 + m(5) * t17 + m(6) * t3; (qJD(1) * t434 + t340 * t9 - t342 * t8) * t621 + (((t166 * t348 - t168 * t351 + t119) * t246 - (t165 * t348 - t167 * t351 + t117) * t247 + (-t203 * t348 + t207 * t351 + t198) * t311 + t77 * qJD(5)) * t339 + (qJD(5) * t434 - t630) * t341) * t622 + ((t166 * t233 - t168 * t234) * t246 - (t165 * t233 - t167 * t234) * t247 + (-t203 * t233 + t207 * t234) * t311 + (t339 * t64 + t36 * t554) * qJD(5) + ((qJD(5) * t35 + t385) * t341 + t356) * t340) * t623 + (qJD(1) * t438 + t340 * t7 - t342 * t6) * t624 + (qJD(1) * t436 + t340 * t5 - t342 * t4) * t625 + ((-t166 * t235 - t168 * t236) * t246 - (-t165 * t235 - t167 * t236) * t247 + (t203 * t235 + t207 * t236) * t311 + (t339 * t65 + t37 * t559) * qJD(5) + ((qJD(5) * t38 + t385) * t341 + t356) * t342) * t626 - t14 * t507 / 0.2e1 + t437 * t627 + t439 * t628 + t435 * t449 - (((t340 * t529 - t342 * t530) * t352 + (t340 * t532 + t342 * t533) * t349 + t632 * t341 + t383 * t339) * qJD(3) + (t339 * t528 + t341 * t527 + t349 * t522 + t352 * t521) * qJD(1)) * qJD(1) / 0.2e1 + (t671 * t342 + t670 * t340 + (t340 * t664 + t342 * t663) * qJD(1)) * t609 + ((-t511 * t564 + t517) * t340 + (t374 + (t340 * t565 + t357) * qJD(3)) * t342 + (-t511 * t562 + t516) * t340 + (t373 + (-t631 * t342 + (t563 + t382) * t340) * qJD(3)) * t342) * t475 + ((-t510 * t565 - t517) * t342 + (t374 + (t342 * t564 + t357) * qJD(3)) * t340 + (-t510 * t563 - t516) * t342 + (t373 + (t382 * t340 + (t562 - t631) * t342) * qJD(3)) * t340) * t472 - (t10 * t340 + t11 * t342) * t506 / 0.2e1 + (-t39 * (qJD(1) * t227 + t169 * t311 - t212 * t247) - t40 * (-t170 * t311 - t212 * t246 + (-t229 - t502) * qJD(1)) - ((-t128 * t39 + t130 * t40) * t339 + t364 * t341) * qJD(5) + t3 * t540 + (t3 * t545 + t39 * t412 + (qJD(1) * t40 + t13) * t452) * t342 + (t12 * t452 + t40 * t412 + t3 * t546 + t39 * (t210 + t268) * qJD(1)) * t340 - t651 * qJD(3) * (-t270 - t615) + (t169 * t246 + t170 * t247 - (-t227 * t340 - t229 * t342 + t456) * qJD(3) + t494 + (qJD(1) * t546 + t146 + t59) * t342 + (t147 + t60 + (-t130 + t536) * qJD(1)) * t340) * t30) * m(6) + (-(t249 * t97 - t598) * qJD(1) - (t96 * (-t249 * t340 - t250 * t342) + t433 * t443) * qJD(3) + t49 * t417 + t96 * t369 + t433 * t278 + (-t74 * t340 - t75 * t342 + (-t342 * t98 + t595) * qJD(1)) * t302) * m(4) + (-(t68 * t223 + t69 * (-t224 - t502)) * qJD(1) - (t52 * t456 + (-t52 * t224 + t467 * t68) * t342 + (-t52 * t223 + t467 * t69) * t340) * qJD(3) + t17 * t540 + t52 * t494 + (t42 * t468 + t68 * t448 + t17 * t197 + t52 * t115 + (t52 * t196 + t447) * qJD(1)) * t342 + (t41 * t468 + t69 * t448 + t17 * t196 + t52 * t116 + (t52 * t537 + t591) * qJD(1)) * t340) * m(5) + (t1 + t669 * qJD(1) + ((t665 * qJD(1) + t652 * t342) * t342 + (t649 * t340 + t666 * qJD(1) + (-t650 + t653) * t342) * t340) * t647) * t620 + (t2 + t668 * qJD(1) + ((t667 * qJD(1) + t650 * t342) * t342 + (t653 * t340 + t689 * qJD(1) + (-t649 + t652) * t342) * t340) * t647) * t619 + (t10 + t673 + t675) * t515 / 0.2e1 + (t11 + t672 + t676) * t477; 0.2e1 * (-t602 / 0.2e1 + t601 / 0.2e1) * m(6) + 0.2e1 * (t41 * t619 + t42 * t620) * m(5); t1 * t560 / 0.2e1 + (t339 * t436 - t341 * t65) * t627 + ((qJD(3) * t436 - t15) * t341 + (-qJD(1) * t437 + qJD(3) * t65 + t340 * t4 + t342 * t5) * t339) * t625 + t2 * t561 / 0.2e1 + (t339 * t438 - t341 * t64) * t628 + ((qJD(3) * t438 - t16) * t341 + (-qJD(1) * t439 + qJD(3) * t64 + t340 * t6 + t342 * t7) * t339) * t624 + t14 * t476 - t341 * (t592 + t593 + t607 + t611 - t612) / 0.2e1 + (t339 * t434 - t341 * t77) * t449 + ((qJD(3) * t434 - t23) * t341 + (-qJD(1) * t435 + qJD(3) * t77 + t340 * t8 + t342 * t9) * t339) * t621 + (t235 * t629 + t371 * t236 - t342 * t381) * t626 + (-t233 * t629 + t234 * t371 - t340 * t381) * t623 + (t384 * t341 + (-t348 * t629 + t351 * t371) * t339) * t622 + (-t492 / 0.2e1 + t341 * t472) * t11 + (t339 * t477 + t341 * t474) * t10 + ((qJD(3) * t364 - t12 * t130 + t13 * t128 + t39 * t60 - t40 * t59) * t341 + (t39 * (-qJD(3) * t128 + t148 * t340) + t40 * (qJD(3) * t130 - t148 * t342) + t3 * t423 + t30 * (-t128 * t515 - t130 * t513 - t340 * t59 + t342 * t60) + (qJD(1) * t651 + t601 - t602) * t210) * t339 - t39 * (-t157 * t311 - t247 * t248) - t40 * (t158 * t311 - t246 * t248) - t30 * (t157 * t246 + t158 * t247)) * m(6);];
tauc = t18(:);
