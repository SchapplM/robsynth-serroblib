% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR7_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:05:58
% EndTime: 2019-12-31 17:06:36
% DurationCPUTime: 32.98s
% Computational Cost: add. (13823->804), mult. (21343->1099), div. (0->0), fcn. (19535->8), ass. (0->403)
t680 = Icges(3,3) + Icges(4,3);
t334 = qJ(2) + pkin(7);
t316 = sin(t334);
t317 = cos(t334);
t339 = sin(qJ(2));
t342 = cos(qJ(2));
t673 = Icges(3,5) * t342 + Icges(4,5) * t317 - Icges(3,6) * t339 - Icges(4,6) * t316;
t343 = cos(qJ(1));
t679 = t680 * t343;
t340 = sin(qJ(1));
t532 = t340 * t342;
t535 = t339 * t340;
t538 = t317 * t340;
t540 = t316 * t340;
t665 = -Icges(3,5) * t532 - Icges(4,5) * t538 + Icges(3,6) * t535 + Icges(4,6) * t540 + t679;
t674 = t680 * t340 + t673 * t343;
t675 = Icges(3,5) * t339 + Icges(4,5) * t316 + Icges(3,6) * t342 + Icges(4,6) * t317;
t558 = Icges(4,6) * t343;
t196 = Icges(4,4) * t538 - Icges(4,2) * t540 - t558;
t559 = Icges(3,6) * t343;
t206 = Icges(3,4) * t532 - Icges(3,2) * t535 - t559;
t678 = t196 * t316 + t206 * t339;
t287 = Icges(4,4) * t540;
t563 = Icges(4,5) * t343;
t198 = Icges(4,1) * t538 - t287 - t563;
t304 = Icges(3,4) * t535;
t564 = Icges(3,5) * t343;
t208 = Icges(3,1) * t532 - t304 - t564;
t660 = -t198 * t317 - t208 * t342 + t678;
t657 = -t660 * t340 + t665 * t343;
t568 = Icges(4,4) * t316;
t257 = Icges(4,1) * t317 - t568;
t199 = Icges(4,5) * t340 + t257 * t343;
t569 = Icges(3,4) * t339;
t278 = Icges(3,1) * t342 - t569;
t209 = Icges(3,5) * t340 + t278 * t343;
t676 = -t199 * t538 - t209 * t532;
t254 = Icges(4,2) * t317 + t568;
t307 = Icges(4,4) * t317;
t256 = Icges(4,1) * t316 + t307;
t275 = Icges(3,2) * t342 + t569;
t326 = Icges(3,4) * t342;
t277 = Icges(3,1) * t339 + t326;
t672 = t254 * t316 - t256 * t317 + t275 * t339 - t277 * t342;
t671 = t343 * t674 + t676;
t530 = t342 * t343;
t537 = t317 * t343;
t619 = -t199 * t537 - t209 * t530 - t340 * t674;
t670 = -t198 * t537 - t208 * t530 + t340 * t665;
t625 = t675 * t343;
t624 = t675 * t340;
t417 = -Icges(4,2) * t316 + t307;
t197 = Icges(4,6) * t340 + t343 * t417;
t418 = -Icges(3,2) * t339 + t326;
t207 = Icges(3,6) * t340 + t343 * t418;
t669 = t197 * t316 + t207 * t339;
t656 = -t197 * t540 - t207 * t535 - t671;
t534 = t339 * t343;
t539 = t316 * t343;
t655 = -t196 * t539 - t206 * t534 - t670;
t654 = -t197 * t539 - t207 * t534 - t619;
t634 = t196 * t317 + t198 * t316 + t206 * t342 + t208 * t339;
t633 = t197 * t317 + t199 * t316 + t207 * t342 + t209 * t339;
t667 = -t672 * t340 - t625;
t666 = -t672 * t343 + t624;
t663 = t675 * qJD(2);
t662 = t199 * t317 + t209 * t342 - t669;
t341 = cos(qJ(4));
t531 = t341 * t343;
t338 = sin(qJ(4));
t536 = t338 * t340;
t233 = t317 * t536 + t531;
t529 = t343 * t338;
t533 = t340 * t341;
t234 = t317 * t533 - t529;
t112 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t540;
t220 = Icges(5,4) * t234;
t115 = -Icges(5,2) * t233 + Icges(5,6) * t540 + t220;
t219 = Icges(5,4) * t233;
t119 = -Icges(5,1) * t234 - Icges(5,5) * t540 + t219;
t644 = t115 * t338 + t119 * t341;
t46 = -t112 * t317 - t644 * t316;
t238 = t417 * qJD(2);
t239 = t257 * qJD(2);
t264 = t418 * qJD(2);
t265 = t278 * qJD(2);
t661 = -t238 * t316 + t239 * t317 - t264 * t339 + t265 * t342 + (-t254 * t317 - t256 * t316 - t275 * t342 - t277 * t339) * qJD(2) + t675 * qJD(1);
t659 = t674 * qJD(1);
t658 = t672 * qJD(1) + t673 * qJD(2);
t653 = t666 * qJD(1);
t383 = qJD(2) * t254;
t106 = -t343 * t383 + (-t340 * t417 + t558) * qJD(1);
t385 = qJD(2) * t256;
t108 = -t343 * t385 + (-t257 * t340 + t563) * qJD(1);
t384 = qJD(2) * t275;
t143 = -t343 * t384 + (-t340 * t418 + t559) * qJD(1);
t386 = qJD(2) * t277;
t145 = -t343 * t386 + (-t278 * t340 + t564) * qJD(1);
t652 = -t633 * qJD(2) - t106 * t316 + t108 * t317 - t143 * t339 + t145 * t342 + t659;
t107 = qJD(1) * t197 - t340 * t383;
t109 = qJD(1) * t199 - t340 * t385;
t144 = qJD(1) * t207 - t340 * t384;
t146 = qJD(1) * t209 - t340 * t386;
t651 = t665 * qJD(1) + t634 * qJD(2) + t107 * t316 - t109 * t317 + t144 * t339 - t146 * t342;
t650 = (t654 * t340 - t655 * t343) * qJD(2);
t649 = (t656 * t340 - t657 * t343) * qJD(2);
t648 = t667 * qJD(1);
t591 = pkin(3) * t317;
t262 = pkin(6) * t316 + t591;
t230 = t262 * t340;
t332 = t343 * pkin(5);
t291 = pkin(1) * t340 - t332;
t337 = -qJ(3) - pkin(5);
t310 = t343 * t337;
t593 = pkin(2) * t342;
t312 = pkin(1) + t593;
t499 = -t340 * t312 - t310;
t192 = t291 + t499;
t514 = t192 - t291;
t429 = rSges(5,1) * t234 - rSges(5,2) * t233;
t121 = rSges(5,3) * t540 + t429;
t428 = rSges(5,1) * t341 - rSges(5,2) * t338;
t187 = -rSges(5,3) * t317 + t316 * t428;
t486 = qJD(4) * t316;
t487 = qJD(2) * t343;
t249 = -t340 * t486 + t487;
t485 = qJD(4) * t317;
t293 = qJD(1) - t485;
t318 = qJD(3) * t340;
t592 = pkin(3) * t316;
t261 = -pkin(6) * t317 + t592;
t594 = pkin(2) * t339;
t448 = -t261 - t594;
t433 = t448 * t343;
t623 = qJD(2) * t433 - t121 * t293 - t249 * t187 + t318;
t38 = (-t230 + t514) * qJD(1) + t623;
t235 = -t317 * t529 + t533;
t236 = t317 * t531 + t536;
t123 = t236 * rSges(5,1) + t235 * rSges(5,2) + rSges(5,3) * t539;
t232 = pkin(3) * t537 + pkin(6) * t539;
t488 = qJD(2) * t340;
t248 = t343 * t486 + t488;
t470 = t261 * t488;
t299 = t488 * t594;
t497 = qJD(3) * t343 + t299;
t331 = t340 * pkin(5);
t292 = t343 * pkin(1) + t331;
t297 = t343 * t312;
t441 = -t337 * t340 + t297;
t193 = t441 - t292;
t511 = t193 + t292;
t39 = -t470 + t123 * t293 - t187 * t248 + (t232 + t511) * qJD(1) - t497;
t575 = t340 * t39;
t647 = t343 * t38 + t575;
t646 = t660 * qJD(1) - t663 * t340 + t659;
t645 = -t663 * t343 + (-t340 * t673 - t662 + t679) * qJD(1);
t34 = t112 * t540 - t115 * t233 - t119 * t234;
t114 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t539;
t567 = Icges(5,4) * t236;
t117 = Icges(5,2) * t235 + Icges(5,6) * t539 + t567;
t221 = Icges(5,4) * t235;
t120 = Icges(5,1) * t236 + Icges(5,5) * t539 + t221;
t35 = t114 * t540 - t233 * t117 + t234 * t120;
t415 = Icges(5,5) * t341 - Icges(5,6) * t338;
t180 = -Icges(5,3) * t317 + t316 * t415;
t565 = Icges(5,4) * t341;
t416 = -Icges(5,2) * t338 + t565;
t182 = -Icges(5,6) * t317 + t316 * t416;
t566 = Icges(5,4) * t338;
t419 = Icges(5,1) * t341 - t566;
t184 = -Icges(5,5) * t317 + t316 * t419;
t61 = t180 * t540 - t182 * t233 + t184 * t234;
t10 = t248 * t35 - t249 * t34 + t293 * t61;
t36 = t112 * t539 + t235 * t115 - t119 * t236;
t37 = t114 * t539 + t235 * t117 + t236 * t120;
t62 = t180 * t539 + t182 * t235 + t184 * t236;
t11 = t248 * t37 - t249 * t36 + t62 * t293;
t643 = 0.2e1 * qJD(2);
t642 = t648 + t649;
t641 = t650 + t653;
t640 = t658 * t340 + t661 * t343;
t639 = t661 * t340 - t658 * t343;
t638 = t660 * qJD(2) - t107 * t317 - t109 * t316 - t144 * t342 - t146 * t339;
t637 = t662 * qJD(2) + t106 * t317 + t108 * t316 + t143 * t342 + t145 * t339;
t636 = rSges(3,2) * t339;
t200 = rSges(4,1) * t538 - rSges(4,2) * t540 - t343 * rSges(4,3);
t327 = t340 * rSges(4,3);
t201 = rSges(4,1) * t537 - rSges(4,2) * t539 + t327;
t407 = t200 * t340 + t201 * t343;
t518 = -t340 * t192 + t343 * t193;
t379 = t407 + t518;
t519 = -t192 * t488 + t193 * t487;
t63 = qJD(2) * t407 + t519;
t632 = qJD(2) * t379 + t63;
t371 = t206 * t343 - t207 * t340;
t372 = t196 * t343 - t197 * t340;
t608 = t340 * (-t254 * t343 + t199) - t343 * (-Icges(4,2) * t538 + t198 - t287);
t609 = t340 * (-t275 * t343 + t209) - t343 * (-Icges(3,2) * t532 + t208 - t304);
t628 = -t316 * t608 + t372 * t317 - t339 * t609 + t371 * t342;
t501 = t277 + t418;
t502 = -t275 + t278;
t504 = t256 + t417;
t505 = -t254 + t257;
t627 = (-t316 * t504 + t317 * t505 - t339 * t501 + t342 * t502) * qJD(1);
t626 = t673 * qJD(1);
t186 = qJD(1) * t192;
t269 = qJD(1) * t291;
t621 = t186 - t269;
t620 = t665 + t669;
t440 = qJD(1) * t317 - qJD(4);
t462 = t316 * t487;
t616 = t340 * t440 + t462;
t181 = Icges(5,3) * t316 + t317 * t415;
t410 = -t182 * t338 + t184 * t341;
t413 = -t117 * t338 + t120 * t341;
t607 = t248 * (-t180 * t343 - t413) - t249 * (-t180 * t340 + t644) + t293 * (t181 - t410);
t213 = (-Icges(5,2) * t341 - t566) * t316;
t606 = t248 * (-Icges(5,2) * t236 + t120 + t221) - t249 * (-Icges(5,2) * t234 - t119 - t219) + t293 * (t184 + t213);
t482 = qJD(2) * qJD(4);
t459 = t317 * t482;
t174 = qJD(1) * t248 + t340 * t459;
t605 = t174 / 0.2e1;
t175 = qJD(1) * t249 + t343 * t459;
t604 = t175 / 0.2e1;
t603 = -t248 / 0.2e1;
t602 = t248 / 0.2e1;
t601 = -t249 / 0.2e1;
t600 = t249 / 0.2e1;
t599 = -t293 / 0.2e1;
t598 = t293 / 0.2e1;
t597 = t340 / 0.2e1;
t596 = -t343 / 0.2e1;
t595 = -rSges(5,3) - pkin(6);
t401 = t293 * t341;
t463 = t316 * t488;
t101 = t340 * t401 + (-t343 * t440 + t463) * t338;
t402 = t293 * t338;
t489 = qJD(2) * t316;
t102 = t440 * t531 + (-t341 * t489 + t402) * t340;
t490 = qJD(1) * t343;
t377 = t316 * t490 + t317 * t488;
t54 = Icges(5,5) * t102 + Icges(5,6) * t101 + Icges(5,3) * t377;
t56 = Icges(5,4) * t102 + Icges(5,2) * t101 + Icges(5,6) * t377;
t58 = Icges(5,1) * t102 + Icges(5,4) * t101 + Icges(5,5) * t377;
t8 = (-qJD(2) * t644 - t54) * t317 + (qJD(2) * t112 - t338 * t56 + t341 * t58 + (-t115 * t341 + t119 * t338) * qJD(4)) * t316;
t590 = t8 * t249;
t100 = -t341 * t616 + t343 * t402;
t467 = t317 * t487;
t491 = qJD(1) * t340;
t473 = t316 * t491;
t376 = t467 - t473;
t99 = t338 * t616 + t343 * t401;
t53 = Icges(5,5) * t100 + Icges(5,6) * t99 + Icges(5,3) * t376;
t55 = Icges(5,4) * t100 + Icges(5,2) * t99 + Icges(5,6) * t376;
t57 = Icges(5,1) * t100 + Icges(5,4) * t99 + Icges(5,5) * t376;
t9 = (qJD(2) * t413 - t53) * t317 + (qJD(2) * t114 - t338 * t55 + t341 * t57 + (-t117 * t341 - t120 * t338) * qJD(4)) * t316;
t589 = t9 * t248;
t587 = qJD(1) / 0.2e1;
t586 = pkin(1) - t312;
t210 = (-Icges(5,5) * t338 - Icges(5,6) * t341) * t316;
t96 = qJD(2) * t181 + qJD(4) * t210;
t183 = Icges(5,6) * t316 + t317 * t416;
t97 = qJD(2) * t183 + qJD(4) * t213;
t185 = Icges(5,5) * t316 + t317 * t419;
t216 = (-Icges(5,1) * t338 - t565) * t316;
t98 = qJD(2) * t185 + qJD(4) * t216;
t18 = (qJD(2) * t410 - t96) * t317 + (qJD(2) * t180 - t338 * t97 + t341 * t98 + (-t182 * t341 - t184 * t338) * qJD(4)) * t316;
t460 = t316 * t482;
t65 = -t180 * t317 + t316 * t410;
t585 = t18 * t293 + t65 * t460;
t584 = rSges(3,1) * t342;
t583 = rSges(4,1) * t317;
t582 = rSges(5,3) * t316;
t188 = t317 * t428 + t582;
t222 = (-rSges(5,1) * t338 - rSges(5,2) * t341) * t316;
t103 = qJD(2) * t188 + qJD(4) * t222;
t272 = pkin(6) * t467;
t378 = -t317 * t491 - t462;
t139 = pkin(3) * t378 - pkin(6) * t473 + t272;
t241 = t262 * qJD(2);
t315 = pkin(5) * t490;
t466 = t339 * t487;
t400 = -pkin(2) * t466 + t318;
t137 = -t315 + (t340 * t586 - t310) * qJD(1) + t400;
t260 = qJD(1) * (-pkin(1) * t491 + t315);
t483 = qJD(1) * qJD(3);
t475 = qJD(1) * t137 + t340 * t483 + t260;
t480 = qJD(2) ^ 2 * t593;
t479 = t100 * rSges(5,1) + t99 * rSges(5,2) + rSges(5,3) * t467;
t59 = -rSges(5,3) * t473 + t479;
t12 = -t340 * t480 + qJD(1) * t139 - t103 * t248 - t175 * t187 + t293 * t59 + (qJD(1) * t433 + t123 * t486 - t241 * t340) * qJD(2) + t475;
t580 = t12 * t343;
t140 = t377 * pkin(6) + (t317 * t490 - t463) * pkin(3);
t500 = qJD(1) * t299 + t343 * t483;
t474 = t337 * t491 + t497;
t138 = (-t343 * t586 - t331) * qJD(1) - t474;
t268 = t292 * qJD(1);
t523 = -t138 - t268;
t430 = rSges(5,1) * t102 + rSges(5,2) * t101;
t60 = rSges(5,3) * t377 + t430;
t13 = -t343 * t480 - t103 * t249 + t174 * t187 - t293 * t60 + (-t121 * t486 - t241 * t343) * qJD(2) + (-t140 + t470 + t523) * qJD(1) + t500;
t579 = t13 * t340;
t328 = t340 * rSges(3,3);
t574 = t46 * t174;
t47 = -t114 * t317 + t316 * t413;
t573 = t47 * t175;
t258 = rSges(4,1) * t316 + rSges(4,2) * t317;
t450 = -t258 - t594;
t434 = t343 * t450;
t398 = qJD(2) * t434;
t375 = t318 + t398;
t68 = (-t200 + t514) * qJD(1) + t375;
t572 = t68 * t258;
t496 = rSges(3,2) * t535 + t343 * rSges(3,3);
t227 = rSges(3,1) * t532 - t496;
t281 = rSges(3,1) * t339 + rSges(3,2) * t342;
t468 = t281 * t487;
t128 = -t468 + (-t227 - t291) * qJD(1);
t553 = t128 * t340;
t552 = t128 * t343;
t469 = t281 * t488;
t228 = rSges(3,1) * t530 - rSges(3,2) * t534 + t328;
t506 = t228 + t292;
t129 = qJD(1) * t506 - t469;
t251 = t281 * t343;
t551 = t129 * t251;
t525 = t121 + t230;
t524 = t123 + t232;
t512 = -t193 - t232;
t503 = rSges(4,2) * t473 + rSges(4,3) * t490;
t498 = rSges(3,3) * t490 + t491 * t636;
t481 = pkin(2) * t534;
t478 = qJD(2) * t593;
t477 = t138 * t488 + (t137 - t186) * t487;
t476 = t343 * t137 + t340 * t138 - t192 * t490;
t471 = t258 * t488;
t465 = t342 * t487;
t461 = -pkin(1) - t584;
t457 = t490 / 0.2e1;
t456 = t489 / 0.2e1;
t455 = -t488 / 0.2e1;
t454 = t488 / 0.2e1;
t452 = t487 / 0.2e1;
t259 = -rSges(4,2) * t316 + t583;
t449 = -t259 - t593;
t439 = (-t340 ^ 2 - t343 ^ 2) * t594;
t438 = -t187 + t448;
t435 = qJD(4) * t456;
t431 = t584 - t636;
t427 = t34 * t343 - t340 * t35;
t426 = t34 * t340 + t343 * t35;
t425 = t340 * t37 - t343 * t36;
t424 = t340 * t36 + t343 * t37;
t423 = t340 * t47 - t343 * t46;
t422 = t340 * t46 + t343 * t47;
t412 = t121 * t343 - t123 * t340;
t411 = -t129 * t340 - t552;
t399 = -t103 - t241 - t478;
t240 = t259 * qJD(2);
t397 = -qJD(2) * t240 - t480;
t250 = t281 * t340;
t223 = t258 * t340;
t124 = (t227 * t340 + t228 * t343) * qJD(2);
t380 = -t262 - t582;
t374 = -t112 * t249 + t114 * t248 + t180 * t293;
t373 = (-Icges(5,5) * t233 - Icges(5,6) * t234) * t249 - (Icges(5,5) * t235 - Icges(5,6) * t236) * t248 - t210 * t293;
t370 = t316 * t373;
t362 = (Icges(5,1) * t235 - t117 - t567) * t248 - (-Icges(5,1) * t233 - t115 - t220) * t249 + (-t182 + t216) * t293;
t29 = t121 * t248 + t123 * t249 + (t230 * t340 + t232 * t343) * qJD(2) + t519;
t355 = t29 * t412 + (t340 * t38 - t343 * t39) * t187;
t346 = t607 * t316;
t266 = t431 * qJD(2);
t231 = t261 * t343;
t229 = t261 * t340;
t224 = t258 * t343;
t164 = t187 * t343;
t163 = t187 * t340;
t162 = t184 * t343;
t161 = t184 * t340;
t160 = t182 * t343;
t159 = t182 * t340;
t156 = rSges(5,1) * t235 - rSges(5,2) * t236;
t155 = -rSges(5,1) * t233 - rSges(5,2) * t234;
t148 = -qJD(2) * t250 + (t343 * t431 + t328) * qJD(1);
t147 = -rSges(3,2) * t465 + (-t342 * t491 - t466) * rSges(3,1) + t498;
t111 = -qJD(2) * t223 + (t259 * t343 + t327) * qJD(1);
t110 = rSges(4,1) * t378 - rSges(4,2) * t467 + t503;
t75 = -t266 * t487 + (-t148 - t268 + t469) * qJD(1);
t74 = -t266 * t488 + t260 + (t147 - t468) * qJD(1);
t69 = -t471 + (t201 + t511) * qJD(1) - t497;
t43 = t397 * t343 + (-t111 + t471 + t523) * qJD(1) + t500;
t42 = t397 * t340 + (t110 + t398) * qJD(1) + t475;
t16 = t101 * t182 + t102 * t184 + t180 * t377 - t233 * t97 + t234 * t98 + t540 * t96;
t15 = t100 * t184 + t180 * t376 + t182 * t99 + t235 * t97 + t236 * t98 + t539 * t96;
t14 = t248 * t47 - t249 * t46 + t293 * t65;
t7 = t121 * t175 - t123 * t174 + t248 * t60 + t249 * t59 + (t139 * t343 + t140 * t340 + (t230 * t343 + t340 * t512) * qJD(1)) * qJD(2) + t477;
t6 = t101 * t117 + t102 * t120 + t114 * t377 - t233 * t55 + t234 * t57 + t53 * t540;
t5 = t101 * t115 - t102 * t119 + t112 * t377 - t233 * t56 + t234 * t58 + t54 * t540;
t4 = t100 * t120 + t114 * t376 + t117 * t99 + t235 * t55 + t236 * t57 + t53 * t539;
t3 = -t100 * t119 + t112 * t376 + t115 * t99 + t235 * t56 + t236 * t58 + t539 * t54;
t2 = t16 * t293 + t174 * t34 + t175 * t35 + t248 * t6 - t249 * t5 + t460 * t61;
t1 = t15 * t293 + t174 * t36 + t175 * t37 + t248 * t4 - t249 * t3 + t460 * t62;
t17 = [t573 / 0.2e1 + t574 / 0.2e1 + t11 * t600 + t15 * t602 + t62 * t604 + t61 * t605 + t585 + t589 / 0.2e1 - t590 / 0.2e1 + (t16 + t11) * t601 + ((((t674 + t678) * t343 + t656 + t670 + t676) * t343 - t619 * t340) * qJD(2) + t653) * t452 + (-t672 * qJD(2) + t238 * t317 + t239 * t316 + t264 * t342 + t265 * t339) * qJD(1) + (-(-qJD(1) * t230 - t38 + t621 + t623) * t39 + t13 * (-t429 + t499) + t38 * (-t430 + t474) + t12 * (t297 + t524) + t39 * (-pkin(3) * t462 + t272 + t400 + t479) + (t13 * t380 - t12 * t337 + t38 * (t317 * t595 + t592) * qJD(2)) * t340 + ((t316 * t595 - t312 - t591) * t575 + (t38 * (-t312 + t380) - t39 * t337) * t343) * qJD(1)) * m(5) + (t43 * (-t200 + t499) + t68 * t474 + t42 * (t201 + t441) + t69 * (t318 + t503) + (t340 * t572 + t434 * t69) * qJD(2) + ((-t68 * rSges(4,3) + t69 * (-t312 - t583)) * t340 + (t68 * (-t259 - t312) - t69 * t337) * t343) * qJD(1) - (-qJD(1) * t200 + t375 + t621 - t68) * t69) * m(4) + (-(-qJD(1) * t227 - t128 - t269 - t468) * t129 + t75 * (t340 * t461 + t332 + t496) + t74 * t506 + t129 * (t315 + t498) + (t281 * t553 - t551) * qJD(2) + ((-pkin(1) - t431) * t552 + (t128 * (-rSges(3,3) - pkin(5)) + t129 * t461) * t340) * qJD(1)) * m(3) + (((t343 * t620 + t619 + t654) * t343 + (t340 * t620 + t655 + t671) * t340) * qJD(2) + t642 - t648) * t455 + (t637 + t640) * t454 - (-t638 + t639 + t641) * t487 / 0.2e1 + ((t634 + t667) * t340 + (t633 + t666) * t343) * qJD(2) * t587; (((t160 * t338 - t162 * t341 + t114) * t248 - (t159 * t338 - t161 * t341 + t112) * t249 + (-t183 * t338 + t185 * t341 + t180) * t293 + t65 * qJD(4)) * t316 + (qJD(4) * t422 - t607) * t317) * t599 - t14 * t486 / 0.2e1 + ((t160 * t233 - t162 * t234) * t248 - (t159 * t233 - t161 * t234) * t249 + (-t183 * t233 + t185 * t234) * t293 + (t316 * t61 + t35 * t537) * qJD(4) + ((qJD(4) * t34 + t374) * t317 + t346) * t340) * t600 + (qJD(1) * t426 + t340 * t6 - t343 * t5) * t601 + (qJD(1) * t424 - t3 * t343 + t340 * t4) * t602 + ((-t160 * t235 - t162 * t236) * t248 - (-t159 * t235 - t161 * t236) * t249 + (t183 * t235 + t185 * t236) * t293 + (t316 * t62 + t36 * t538) * qJD(4) + ((qJD(4) * t37 + t374) * t317 + t346) * t343) * t603 + t425 * t604 + (qJD(1) * t422 + t340 * t9 - t343 * t8) * t598 - t174 * t427 / 0.2e1 + t423 * t435 - ((t372 * t316 + t317 * t608 + t371 * t339 + t342 * t609) * qJD(2) + (t316 * t505 + t317 * t504 + t339 * t502 + t342 * t501) * qJD(1)) * qJD(1) / 0.2e1 + (t638 * t343 + t637 * t340 + (t340 * t634 + t343 * t633) * qJD(1)) * t587 + ((-t488 * t625 + t626) * t340 + ((t340 * t624 + t628) * qJD(2) + t627) * t343) * t455 + ((-t487 * t624 - t626) * t343 + ((t343 * t625 + t628) * qJD(2) + t627) * t340) * t452 - (t10 * t340 + t11 * t343) * t485 / 0.2e1 + (t7 * t518 + (t38 * t399 + t7 * t524 + (qJD(1) * t39 + t13) * t438) * t343 + (t12 * t438 + t39 * t399 + t7 * t525 + t38 * (t187 + t261) * qJD(1)) * t340 - t38 * (qJD(1) * t229 + t163 * t293 - t188 * t249) - t39 * (-t164 * t293 - t188 * t248 + (-t231 - t481) * qJD(1)) - ((-t121 * t38 + t123 * t39) * t316 + t355 * t317) * qJD(4) - t647 * qJD(2) * (-t262 - t593) + (t476 + (t525 * qJD(1) + t139 + t59) * t343 + (t140 + t60 + (-t123 + t512) * qJD(1)) * t340 + t163 * t248 + t164 * t249 - (-t229 * t340 - t231 * t343 + t439) * qJD(2)) * t29) * m(5) + (t43 * t434 - t68 * pkin(2) * t465 + (t110 * t487 + t111 * t488 + t477) * t379 + t63 * t476 + (-t68 * t240 + t63 * t110 + (t200 * t632 + t450 * t69) * qJD(1)) * t343 + (t42 * t450 + t69 * (-t240 - t478) + t63 * t111 + (t572 + t632 * (-t193 - t201)) * qJD(1)) * t340 - (t68 * t223 + t69 * (-t224 - t481)) * qJD(1) - (t63 * t439 + (-t63 * t224 + t449 * t68) * t343 + (-t63 * t223 + t449 * t69) * t340) * qJD(2)) * m(4) + (-(t128 * t250 - t551) * qJD(1) - (t124 * (-t250 * t340 - t251 * t343) + t411 * t431) * qJD(2) + 0.2e1 * t124 * (t147 * t343 + t148 * t340 + (t227 * t343 - t228 * t340) * qJD(1)) + t411 * t266 + (-t74 * t340 - t75 * t343 + (-t129 * t343 + t553) * qJD(1)) * t281) * m(3) + (t640 * qJD(1) + t1 + ((t654 * qJD(1) + t651 * t343) * t343 + (t645 * t340 + t655 * qJD(1) + (-t646 + t652) * t343) * t340) * t643) * t597 + (t639 * qJD(1) + t2 + ((t656 * qJD(1) + t646 * t343) * t343 + (t652 * t340 + t657 * qJD(1) + (-t645 + t651) * t343) * t340) * t643) * t596 + (t10 + t642 + t649) * t491 / 0.2e1 + (t11 + t641 + t650) * t457; 0.2e1 * (-t580 / 0.2e1 + t579 / 0.2e1) * m(5) + 0.2e1 * (t42 * t596 + t43 * t597) * m(4); t1 * t539 / 0.2e1 + (t316 * t424 - t317 * t62) * t604 + ((qJD(2) * t424 - t15) * t317 + (-qJD(1) * t425 + qJD(2) * t62 + t3 * t340 + t343 * t4) * t316) * t602 + t2 * t540 / 0.2e1 + (t316 * t426 - t317 * t61) * t605 + ((qJD(2) * t426 - t16) * t317 + (qJD(1) * t427 + qJD(2) * t61 + t340 * t5 + t343 * t6) * t316) * t601 + t14 * t456 - t317 * (t573 + t574 + t585 + t589 - t590) / 0.2e1 + (t316 * t422 - t317 * t65) * t435 + ((qJD(2) * t422 - t18) * t317 + (-qJD(1) * t423 + qJD(2) * t65 + t340 * t8 + t343 * t9) * t316) * t598 + (t235 * t606 + t362 * t236 - t343 * t370) * t603 + (-t233 * t606 + t234 * t362 - t340 * t370) * t600 + (t373 * t317 + (-t338 * t606 + t341 * t362) * t316) * t599 + (-t473 / 0.2e1 + t317 * t452) * t11 + (t316 * t457 + t317 * t454) * t10 + ((qJD(2) * t355 - t12 * t123 + t13 * t121 + t38 * t60 - t39 * t59) * t317 + (t38 * (-qJD(2) * t121 + t103 * t340) + t39 * (qJD(2) * t123 - t103 * t343) + t7 * t412 + t29 * (-t121 * t491 - t123 * t490 - t340 * t59 + t343 * t60) + (t647 * qJD(1) + t579 - t580) * t187) * t316 - t38 * (-t155 * t293 - t222 * t249) - t39 * (t156 * t293 - t222 * t248) - t29 * (t155 * t248 + t156 * t249)) * m(5);];
tauc = t17(:);
