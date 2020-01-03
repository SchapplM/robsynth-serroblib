% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR6_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR6_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:23
% EndTime: 2019-12-31 17:30:43
% DurationCPUTime: 70.42s
% Computational Cost: add. (37057->1324), mult. (100840->1874), div. (0->0), fcn. (114776->10), ass. (0->473)
t451 = sin(qJ(1));
t447 = sin(pkin(4));
t571 = qJD(2) * t447;
t435 = t451 * t571;
t450 = sin(qJ(2));
t453 = cos(qJ(2));
t454 = cos(qJ(1));
t615 = cos(pkin(4));
t540 = t451 * t615;
t479 = -t454 * t450 - t453 * t540;
t370 = -qJD(3) * t479 + t435;
t410 = -t450 * t540 + t453 * t454;
t449 = sin(qJ(3));
t633 = cos(qJ(3));
t557 = t447 * t633;
t484 = -t410 * t449 + t451 * t557;
t262 = -qJD(4) * t484 + t370;
t539 = t454 * t615;
t478 = -t450 * t539 - t451 * t453;
t363 = -t449 * t478 + t454 * t557;
t569 = qJD(2) * t454;
t436 = t447 * t569;
t525 = t453 * t539;
t589 = t450 * t451;
t477 = t525 - t589;
t371 = qJD(3) * t477 + t436;
t263 = -qJD(4) * t363 + t371;
t593 = t447 * t450;
t402 = t449 * t593 - t615 * t633;
t441 = qJD(2) * t615 + qJD(1);
t591 = t447 * t453;
t414 = -qJD(3) * t591 + t441;
t348 = qJD(4) * t402 + t414;
t590 = t447 * t454;
t364 = -t449 * t590 - t478 * t633;
t448 = sin(qJ(4));
t452 = cos(qJ(4));
t271 = t364 * t452 - t477 * t448;
t272 = t364 * t448 + t477 * t452;
t122 = Icges(5,5) * t271 - Icges(5,6) * t272 + Icges(5,3) * t363;
t609 = Icges(5,4) * t271;
t126 = Icges(5,2) * t272 - Icges(5,6) * t363 - t609;
t266 = Icges(5,4) * t272;
t129 = -Icges(5,1) * t271 - Icges(5,5) * t363 + t266;
t38 = t122 * t363 + t126 * t272 - t129 * t271;
t592 = t447 * t451;
t368 = t410 * t633 + t449 * t592;
t274 = -t368 * t448 - t452 * t479;
t275 = t368 * t452 - t448 * t479;
t124 = Icges(5,5) * t275 + Icges(5,6) * t274 - Icges(5,3) * t484;
t608 = Icges(5,4) * t275;
t127 = Icges(5,2) * t274 - Icges(5,6) * t484 + t608;
t267 = Icges(5,4) * t274;
t130 = Icges(5,1) * t275 - Icges(5,5) * t484 + t267;
t39 = t363 * t124 - t127 * t272 + t271 * t130;
t403 = t449 * t615 + t450 * t557;
t358 = -t403 * t448 - t452 * t591;
t498 = -t403 * t452 + t448 * t591;
t198 = -Icges(5,5) * t498 + Icges(5,6) * t358 + Icges(5,3) * t402;
t607 = Icges(5,4) * t498;
t199 = Icges(5,2) * t358 + Icges(5,6) * t402 - t607;
t349 = Icges(5,4) * t358;
t200 = -Icges(5,1) * t498 + Icges(5,5) * t402 + t349;
t697 = -t198 * t363 + t199 * t272 - t200 * t271;
t16 = t262 * t39 - t263 * t38 - t697 * t348;
t40 = -t122 * t484 - t126 * t274 - t129 * t275;
t41 = -t124 * t484 + t274 * t127 + t275 * t130;
t66 = -t198 * t484 + t199 * t274 + t200 * t275;
t17 = t262 * t41 - t263 * t40 + t66 * t348;
t48 = t122 * t402 - t126 * t358 + t129 * t498;
t519 = -rSges(5,1) * t271 + rSges(5,2) * t272;
t131 = rSges(5,3) * t363 - t519;
t201 = -rSges(5,1) * t498 + rSges(5,2) * t358 + rSges(5,3) * t402;
t630 = pkin(3) * t364;
t238 = -pkin(8) * t363 - t630;
t334 = pkin(3) * t403 + pkin(8) * t402;
t709 = -t131 * t348 - t201 * t263 + t238 * t414 - t334 * t371;
t281 = Icges(4,5) * t403 - Icges(4,6) * t402 - Icges(4,3) * t591;
t610 = Icges(4,4) * t403;
t282 = -Icges(4,2) * t402 - Icges(4,6) * t591 + t610;
t389 = Icges(4,4) * t402;
t283 = Icges(4,1) * t403 - Icges(4,5) * t591 - t389;
t668 = t281 * t477 + t282 * t363 - t283 * t364;
t203 = Icges(4,5) * t364 - Icges(4,6) * t363 - Icges(4,3) * t477;
t612 = Icges(4,4) * t364;
t207 = Icges(4,2) * t363 + Icges(4,6) * t477 - t612;
t350 = Icges(4,4) * t363;
t210 = -Icges(4,1) * t364 + Icges(4,5) * t477 + t350;
t69 = -t203 * t477 + t207 * t363 - t210 * t364;
t205 = Icges(4,5) * t368 + Icges(4,6) * t484 - Icges(4,3) * t479;
t611 = Icges(4,4) * t368;
t208 = Icges(4,2) * t484 - Icges(4,6) * t479 + t611;
t351 = Icges(4,4) * t484;
t211 = Icges(4,1) * t368 - Icges(4,5) * t479 + t351;
t70 = -t477 * t205 - t363 * t208 + t364 * t211;
t30 = t370 * t70 - t371 * t69 - t668 * t414;
t71 = -t203 * t479 - t207 * t484 - t210 * t368;
t72 = -t205 * t479 + t208 * t484 + t368 * t211;
t82 = -t281 * t479 + t282 * t484 + t283 * t368;
t31 = t370 * t72 - t371 * t71 + t82 * t414;
t74 = -t203 * t591 + t207 * t402 - t210 * t403;
t378 = Icges(3,3) * t615 + (Icges(3,5) * t450 + Icges(3,6) * t453) * t447;
t613 = Icges(3,4) * t450;
t379 = Icges(3,6) * t615 + (Icges(3,2) * t453 + t613) * t447;
t437 = Icges(3,4) * t591;
t380 = Icges(3,1) * t593 + Icges(3,5) * t615 + t437;
t469 = t378 * t590 - t379 * t477 + t380 * t478;
t703 = t469 * t441;
t572 = qJD(1) * t454;
t552 = t447 * t572;
t573 = qJD(1) * t451;
t523 = -pkin(1) * t573 + pkin(6) * t552;
t212 = rSges(4,1) * t364 - rSges(4,2) * t363 - rSges(4,3) * t477;
t284 = rSges(4,1) * t403 - rSges(4,2) * t402 - rSges(4,3) * t591;
t699 = -t212 * t414 - t284 * t371;
t286 = Icges(3,5) * t478 - Icges(3,6) * t477 + Icges(3,3) * t590;
t394 = Icges(3,4) * t478;
t288 = Icges(3,2) * t477 - Icges(3,6) * t590 - t394;
t393 = Icges(3,4) * t477;
t292 = Icges(3,1) * t478 + Icges(3,5) * t590 - t393;
t120 = -t615 * t286 + t447 * (t288 * t453 - t292 * t450);
t287 = Icges(3,5) * t410 + Icges(3,6) * t479 + Icges(3,3) * t592;
t614 = Icges(3,4) * t410;
t290 = Icges(3,2) * t479 + Icges(3,6) * t592 + t614;
t395 = Icges(3,4) * t479;
t293 = Icges(3,1) * t410 + Icges(3,5) * t592 + t395;
t112 = t287 * t592 + t290 * t479 + t410 * t293;
t512 = t288 * t477 + t292 * t478;
t693 = -t512 + t112;
t296 = -rSges(3,1) * t478 + rSges(3,2) * t477 - rSges(3,3) * t590;
t382 = t615 * rSges(3,3) + (rSges(3,1) * t450 + rSges(3,2) * t453) * t447;
t418 = pkin(1) * t451 - pkin(6) * t590;
t416 = qJD(1) * t418;
t175 = -t296 * t441 - t382 * t436 - t416;
t661 = t454 * pkin(1) + pkin(6) * t592;
t673 = t661 * qJD(1);
t692 = t382 * t435 - t673;
t344 = -pkin(2) * t478 - pkin(7) * t477;
t412 = (pkin(2) * t450 - pkin(7) * t453) * t447;
t467 = -t344 * t441 - t412 * t436 - t416;
t691 = -t288 * t479 + t292 * t410;
t331 = qJD(1) * t410 + qJD(2) * t477;
t534 = qJD(1) * t557;
t196 = qJD(3) * t364 + t331 * t449 - t451 * t534;
t330 = -qJD(1) * t479 - qJD(2) * t478;
t545 = qJD(1) * t571;
t426 = t451 * t545;
t278 = qJD(3) * t330 + t426;
t135 = qJD(4) * t196 + t278;
t329 = qJD(1) * t478 + qJD(2) * t479;
t194 = qJD(3) * t368 + t329 * t449 - t454 * t534;
t328 = -qJD(1) * t525 + t441 * t589 - t453 * t569;
t427 = t454 * t545;
t279 = -qJD(3) * t328 + t427;
t136 = qJD(4) * t194 + t279;
t195 = qJD(3) * t484 + t329 * t633 + t449 * t552;
t104 = -qJD(4) * t275 - t195 * t448 - t328 * t452;
t105 = qJD(4) * t274 + t195 * t452 - t328 * t448;
t550 = t453 * t571;
t357 = -qJD(3) * t402 + t550 * t633;
t570 = qJD(2) * t450;
t551 = t447 * t570;
t217 = qJD(4) * t498 - t357 * t448 + t452 * t551;
t218 = qJD(4) * t358 + t357 * t452 + t448 * t551;
t356 = qJD(3) * t403 + t449 * t550;
t97 = Icges(5,5) * t218 + Icges(5,6) * t217 + Icges(5,3) * t356;
t98 = Icges(5,4) * t218 + Icges(5,2) * t217 + Icges(5,6) * t356;
t99 = Icges(5,1) * t218 + Icges(5,4) * t217 + Icges(5,5) * t356;
t19 = t104 * t199 + t105 * t200 + t194 * t198 + t274 * t98 + t275 * t99 - t484 * t97;
t565 = qJD(3) * t450;
t549 = t447 * t565;
t533 = qJD(2) * t549;
t312 = qJD(4) * t356 + t533;
t553 = t447 * t573;
t197 = -qJD(3) * t363 + t331 * t633 + t449 * t553;
t106 = -qJD(4) * t271 - t197 * t448 + t330 * t452;
t107 = -qJD(4) * t272 + t197 * t452 + t330 * t448;
t55 = Icges(5,5) * t107 + Icges(5,6) * t106 + Icges(5,3) * t196;
t57 = Icges(5,4) * t107 + Icges(5,2) * t106 + Icges(5,6) * t196;
t59 = Icges(5,1) * t107 + Icges(5,4) * t106 + Icges(5,5) * t196;
t8 = -t104 * t126 - t105 * t129 + t122 * t194 + t274 * t57 + t275 * t59 - t484 * t55;
t54 = Icges(5,5) * t105 + Icges(5,6) * t104 + Icges(5,3) * t194;
t56 = Icges(5,4) * t105 + Icges(5,2) * t104 + Icges(5,6) * t194;
t58 = Icges(5,1) * t105 + Icges(5,4) * t104 + Icges(5,5) * t194;
t9 = t104 * t127 + t105 * t130 + t124 * t194 + t274 * t56 + t275 * t58 - t484 * t54;
t1 = t135 * t40 + t136 * t41 + t19 * t348 + t262 * t9 - t263 * t8 + t312 * t66;
t84 = Icges(4,5) * t197 - Icges(4,6) * t196 + Icges(4,3) * t330;
t86 = Icges(4,4) * t197 - Icges(4,2) * t196 + Icges(4,6) * t330;
t88 = Icges(4,1) * t197 - Icges(4,4) * t196 + Icges(4,5) * t330;
t21 = t194 * t207 - t195 * t210 - t203 * t328 + t368 * t88 - t479 * t84 + t484 * t86;
t83 = Icges(4,5) * t195 - Icges(4,6) * t194 - Icges(4,3) * t328;
t85 = Icges(4,4) * t195 - Icges(4,2) * t194 - Icges(4,6) * t328;
t87 = Icges(4,1) * t195 - Icges(4,4) * t194 - Icges(4,5) * t328;
t22 = -t194 * t208 + t195 * t211 - t205 * t328 + t368 * t87 - t479 * t83 + t484 * t85;
t219 = Icges(4,5) * t357 - Icges(4,6) * t356 + Icges(4,3) * t551;
t220 = Icges(4,4) * t357 - Icges(4,2) * t356 + Icges(4,6) * t551;
t221 = Icges(4,1) * t357 - Icges(4,4) * t356 + Icges(4,5) * t551;
t36 = -t194 * t282 + t195 * t283 - t219 * t479 + t220 * t484 + t221 * t368 - t281 * t328;
t689 = -t21 * t371 + t370 * t22 + t278 * t71 + t279 * t72 + t36 * t414 + t533 * t82 + t1;
t10 = -t106 * t126 - t107 * t129 + t122 * t196 + t271 * t59 - t272 * t57 + t363 * t55;
t11 = t106 * t127 + t107 * t130 + t124 * t196 + t271 * t58 - t272 * t56 + t363 * t54;
t20 = t106 * t199 + t107 * t200 + t196 * t198 + t271 * t99 - t272 * t98 + t363 * t97;
t2 = -t10 * t263 + t11 * t262 + t135 * t38 + t136 * t39 + t20 * t348 - t312 * t697;
t23 = t196 * t207 - t197 * t210 + t203 * t330 - t363 * t86 + t364 * t88 - t477 * t84;
t24 = -t196 * t208 + t197 * t211 + t205 * t330 - t363 * t85 + t364 * t87 - t477 * t83;
t37 = -t196 * t282 + t197 * t283 - t219 * t477 - t220 * t363 + t221 * t364 + t281 * t330;
t688 = -t23 * t371 + t24 * t370 + t278 * t69 + t279 * t70 + t37 * t414 - t533 * t668 + t2;
t49 = t124 * t402 + t127 * t358 - t130 * t498;
t621 = t49 * t136;
t622 = t48 * t135;
t13 = t124 * t356 + t127 * t217 + t130 * t218 + t358 * t56 + t402 * t54 - t498 * t58;
t625 = t13 * t262;
t12 = t122 * t356 - t126 * t217 - t129 * t218 + t358 * t57 + t402 * t55 - t498 * t59;
t626 = t12 * t263;
t26 = t198 * t356 + t199 * t217 + t200 * t218 + t358 * t98 + t402 * t97 - t498 * t99;
t68 = t198 * t402 + t199 * t358 - t200 * t498;
t628 = t26 * t348 + t68 * t312;
t3 = t621 + t622 + t625 - t626 + t628;
t108 = -t281 * t591 - t282 * t402 + t283 * t403;
t53 = -t220 * t402 + t221 * t403 - t282 * t356 + t283 * t357 + (-t219 * t453 + t281 * t570) * t447;
t618 = t108 * t533 + t53 * t414;
t75 = -t205 * t591 - t208 * t402 + t211 * t403;
t619 = t75 * t279;
t620 = t74 * t278;
t29 = -t208 * t356 + t211 * t357 - t402 * t85 + t403 * t87 + (t205 * t570 - t453 * t83) * t447;
t623 = t29 * t370;
t28 = t207 * t356 - t210 * t357 - t402 * t86 + t403 * t88 + (t203 * t570 - t453 * t84) * t447;
t624 = t28 * t371;
t687 = t618 + t619 + t620 + t623 - t624 + t3;
t685 = t30 + t16;
t684 = t31 + t17;
t100 = rSges(5,1) * t218 + rSges(5,2) * t217 + rSges(5,3) * t356;
t116 = t195 * pkin(3) + pkin(8) * t194;
t133 = t275 * rSges(5,1) + t274 * rSges(5,2) - rSges(5,3) * t484;
t226 = pkin(3) * t357 + pkin(8) * t356;
t240 = t368 * pkin(3) - pkin(8) * t484;
t413 = (pkin(2) * t453 + pkin(7) * t450) * t447;
t388 = qJD(2) * t413;
t494 = -t388 * t451 - t412 * t572;
t215 = t329 * pkin(2) - pkin(7) * t328;
t401 = qJD(1) * t523;
t585 = t441 * t215 + t401;
t60 = t105 * rSges(5,1) + t104 * rSges(5,2) + t194 * rSges(5,3);
t14 = -t100 * t262 + t116 * t414 + t133 * t312 - t136 * t201 - t226 * t370 - t279 * t334 + t348 * t60 + (t240 * t565 + t494) * t571 + t585;
t631 = pkin(3) * t197;
t117 = pkin(8) * t196 + t631;
t216 = pkin(2) * t331 + t330 * pkin(7);
t493 = -qJD(1) * t673 - t216 * t441 + t412 * t426;
t598 = t388 * t454;
t520 = -rSges(5,1) * t107 - rSges(5,2) * t106;
t61 = rSges(5,3) * t196 - t520;
t15 = -t100 * t263 - t117 * t414 - t131 * t312 + t135 * t201 - t226 * t371 + t278 * t334 - t348 * t61 + (t238 * t565 - t598) * t571 + t493;
t42 = t467 + t709;
t347 = t410 * pkin(2) - pkin(7) * t479;
t492 = t347 * t441 - t412 * t435 + t673;
t43 = t133 * t348 - t201 * t262 + t240 * t414 - t334 * t370 + t492;
t586 = t133 + t240;
t587 = t131 - t238;
t616 = t117 + t61;
t617 = t116 + t60;
t669 = t14 * t586 - t15 * t587 - t42 * t616 + t43 * t617;
t579 = t344 * t435 + t347 * t436;
t33 = t131 * t262 + t133 * t263 - t238 * t370 + t240 * t371 + t579;
t665 = t33 * t587;
t298 = t410 * rSges(3,1) + rSges(3,2) * t479 + rSges(3,3) * t592;
t163 = (t296 * t451 + t298 * t454) * t571;
t582 = -t290 * t477 + t293 * t478;
t660 = t215 - t467 + t523;
t468 = t215 * t436 + t216 * t435 + t344 * t427 - t347 * t426;
t7 = t116 * t371 + t117 * t370 + t131 * t136 - t133 * t135 - t238 * t279 - t240 * t278 + t262 * t61 + t263 * t60 + t468;
t659 = t33 * t617 + t7 * t586;
t495 = Icges(4,5) * t633 - Icges(4,6) * t449;
t658 = t370 * (-Icges(4,3) * t410 - t208 * t449 + t211 * t633 - t479 * t495) - t371 * (Icges(4,3) * t478 + t207 * t449 - t210 * t633 - t477 * t495) + t414 * (-t282 * t449 + t283 * t633 - (Icges(4,3) * t450 + t453 * t495) * t447);
t657 = -t17 / 0.2e1;
t656 = -t31 / 0.2e1;
t655 = t135 / 0.2e1;
t654 = t136 / 0.2e1;
t653 = -t262 / 0.2e1;
t652 = t262 / 0.2e1;
t651 = t263 / 0.2e1;
t650 = -t263 / 0.2e1;
t649 = t278 / 0.2e1;
t648 = t279 / 0.2e1;
t647 = t312 / 0.2e1;
t644 = -t348 / 0.2e1;
t643 = t348 / 0.2e1;
t642 = -t370 / 0.2e1;
t641 = t370 / 0.2e1;
t640 = t371 / 0.2e1;
t639 = -t371 / 0.2e1;
t636 = -t414 / 0.2e1;
t635 = t414 / 0.2e1;
t634 = -rSges(5,3) - pkin(8);
t603 = t175 * t454;
t600 = t286 * t451;
t599 = t287 * t454;
t596 = t477 * t449;
t594 = t479 * t449;
t588 = t100 + t226;
t584 = t201 + t334;
t222 = rSges(4,1) * t357 - rSges(4,2) * t356 + rSges(4,3) * t551;
t583 = -t222 - t388;
t581 = -t284 - t412;
t343 = pkin(2) * t477 - pkin(7) * t478;
t346 = pkin(2) * t479 + pkin(7) * t410;
t580 = t343 * t435 + t346 * t436;
t577 = t344 * t592 + t347 * t590;
t406 = (Icges(3,1) * t453 - t613) * t447;
t576 = -t379 + t406;
t575 = -Icges(3,2) * t593 + t380 + t437;
t567 = qJD(3) * t478;
t566 = qJD(3) * t410;
t564 = qJD(4) * t449;
t563 = t449 * t591;
t561 = -t388 - t588;
t89 = t195 * rSges(4,1) - t194 * rSges(4,2) - t328 * rSges(4,3);
t560 = t215 * t590 + t216 * t592 + t344 * t552;
t559 = -t412 - t584;
t111 = -t286 * t592 - t691;
t183 = t329 * rSges(3,1) + t328 * rSges(3,2) + rSges(3,3) * t552;
t214 = t368 * rSges(4,1) + rSges(4,2) * t484 - rSges(4,3) * t479;
t556 = t448 * t633;
t555 = t452 * t633;
t554 = t453 * t633;
t543 = -t571 / 0.2e1;
t542 = t571 / 0.2e1;
t541 = t42 * t584;
t538 = t615 * t344;
t537 = t454 * t581;
t530 = t450 * t542;
t529 = t451 * t543;
t528 = t451 * t542;
t527 = t454 * t543;
t526 = t454 * t542;
t524 = qJD(1) * t542;
t521 = -rSges(3,1) * t331 + rSges(3,2) * t330;
t518 = -rSges(5,1) * t452 + rSges(5,2) * t448;
t517 = -t216 * t615 + t412 * t553;
t516 = -Icges(5,1) * t452 + Icges(5,4) * t448;
t515 = -Icges(5,4) * t452 + Icges(5,2) * t448;
t514 = -Icges(5,5) * t452 + Icges(5,6) * t448;
t513 = t111 * t451 + t112 * t454;
t511 = -t288 * t450 - t292 * t453;
t510 = -t290 * t450 + t293 * t453;
t335 = Icges(3,5) * t477 + Icges(3,6) * t478;
t336 = Icges(3,5) * t479 - Icges(3,6) * t410;
t509 = t335 * t454 - t336 * t451;
t508 = t347 + t661;
t507 = t451 * t524;
t506 = t454 * t524;
t505 = qJD(3) * t530;
t503 = pkin(3) * t633 + pkin(8) * t449;
t502 = t441 * t346 - t413 * t435;
t501 = rSges(4,1) * t633 - rSges(4,2) * t449;
t497 = Icges(4,1) * t633 - Icges(4,4) * t449;
t496 = Icges(4,4) * t633 - Icges(4,2) * t449;
t491 = (rSges(3,1) * t453 - rSges(3,2) * t450) * t447;
t489 = -t344 - t418;
t404 = (Icges(3,5) * t453 - Icges(3,6) * t450) * t447;
t109 = t286 * t590 + t512;
t110 = -t287 * t590 - t582;
t488 = (-t109 * t454 + t110 * t451) * t447;
t487 = (-t111 * t454 + t112 * t451) * t447;
t480 = -t343 * t441 - t413 * t436;
t90 = rSges(4,1) * t197 - rSges(4,2) * t196 + rSges(4,3) * t330;
t476 = t33 * t616 + t587 * t7;
t475 = -(-Icges(5,5) * t272 - Icges(5,6) * t271) * t263 + (Icges(5,5) * t274 - Icges(5,6) * t275) * t262 + (Icges(5,5) * t358 + Icges(5,6) * t498) * t348;
t474 = -(-Icges(4,5) * t363 - Icges(4,6) * t364) * t371 + (Icges(4,5) * t484 - Icges(4,6) * t368) * t370 + (-Icges(4,5) * t402 - Icges(4,6) * t403) * t414;
t387 = qJD(2) * t491;
t472 = (-t382 * t572 - t387 * t451) * t447;
t339 = Icges(3,1) * t477 + t394;
t340 = Icges(3,1) * t479 - t614;
t471 = (-t290 + t340) * t451 - (-t288 + t339) * t454;
t337 = Icges(3,2) * t478 + t393;
t338 = -Icges(3,2) * t410 + t395;
t470 = (-t293 - t338) * t451 - (t292 - t337) * t454;
t466 = -t216 - t673;
t465 = (Icges(5,1) * t274 - t127 - t608) * t262 - (-Icges(5,1) * t272 + t126 - t609) * t263 + (Icges(5,1) * t358 - t199 + t607) * t348;
t464 = (-Icges(5,2) * t275 + t130 + t267) * t262 - (-Icges(5,2) * t271 - t129 - t266) * t263 + (Icges(5,2) * t498 + t200 + t349) * t348;
t463 = (Icges(4,1) * t484 - t208 - t611) * t370 - (-Icges(4,1) * t363 + t207 - t612) * t371 + (-Icges(4,1) * t402 - t282 - t610) * t414;
t462 = (Icges(4,2) * t368 - t211 - t351) * t370 - (Icges(4,2) * t364 + t210 + t350) * t371 + (Icges(4,2) * t403 - t283 + t389) * t414;
t177 = Icges(3,5) * t329 + Icges(3,6) * t328 + Icges(3,3) * t552;
t178 = Icges(3,5) * t331 - Icges(3,6) * t330 + Icges(3,3) * t553;
t179 = Icges(3,4) * t329 + Icges(3,2) * t328 + Icges(3,6) * t552;
t180 = Icges(3,4) * t331 - Icges(3,2) * t330 + Icges(3,6) * t553;
t181 = Icges(3,1) * t329 + Icges(3,4) * t328 + Icges(3,5) * t552;
t182 = Icges(3,1) * t331 - Icges(3,4) * t330 + Icges(3,5) * t553;
t460 = (qJD(1) * t513 - (t180 * t479 + t182 * t410 + t288 * t328 - t292 * t329 + (t178 * t451 - t286 * t572) * t447) * t454 + (t179 * t479 + t181 * t410 + t290 * t328 + t293 * t329 + (t177 * t451 + t287 * t572) * t447) * t451) * t447;
t459 = (t451 * (t179 * t477 - t181 * t478 - t290 * t330 + t293 * t331 + (-t177 * t454 + t287 * t573) * t447) - t454 * (t180 * t477 - t182 * t478 - t288 * t330 - t292 * t331 + (-t178 * t454 - t286 * t573) * t447) + (t109 * t451 + t110 * t454) * qJD(1)) * t447;
t121 = t615 * t287 + (t290 * t453 + t293 * t450) * t447;
t63 = t615 * t178 + (qJD(2) * t511 + t180 * t453 + t182 * t450) * t447;
t64 = t615 * t177 + (qJD(2) * t510 + t179 * t453 + t181 * t450) * t447;
t458 = (t451 * t64 - t454 * t63 + (t120 * t451 + t121 * t454) * qJD(1)) * t447;
t457 = (Icges(5,3) * t368 + t127 * t448 - t130 * t452 - t484 * t514) * t262 - (Icges(5,3) * t364 - t126 * t448 + t129 * t452 + t363 * t514) * t263 + (Icges(5,3) * t403 + t199 * t448 - t200 * t452 + t402 * t514) * t348;
t386 = qJD(2) * t406;
t385 = (Icges(3,4) * t453 - Icges(3,2) * t450) * t571;
t384 = qJD(2) * t404;
t383 = (t453 * t564 + t565) * t447;
t381 = t503 * t591;
t377 = (t448 * t450 + t452 * t554) * t447;
t376 = (-t448 * t554 + t450 * t452) * t447;
t369 = (rSges(4,3) * t450 + t453 * t501) * t447;
t362 = (Icges(4,5) * t450 + t453 * t497) * t447;
t361 = (Icges(4,6) * t450 + t453 * t496) * t447;
t342 = rSges(3,1) * t479 - rSges(3,2) * t410;
t341 = rSges(3,1) * t477 + rSges(3,2) * t478;
t333 = -pkin(3) * t402 + pkin(8) * t403;
t332 = -rSges(4,1) * t402 - rSges(4,2) * t403;
t317 = t615 * t347;
t316 = t479 * t564 + t566;
t315 = t477 * t564 - t567;
t314 = t503 * t479;
t313 = t503 * t477;
t311 = t410 * t448 + t479 * t555;
t310 = t410 * t452 - t479 * t556;
t309 = -t448 * t478 + t477 * t555;
t308 = -t452 * t478 - t477 * t556;
t265 = rSges(5,1) * t377 + rSges(5,2) * t376 + rSges(5,3) * t563;
t261 = Icges(5,1) * t377 + Icges(5,4) * t376 + Icges(5,5) * t563;
t260 = Icges(5,4) * t377 + Icges(5,2) * t376 + Icges(5,6) * t563;
t259 = Icges(5,5) * t377 + Icges(5,6) * t376 + Icges(5,3) * t563;
t258 = rSges(4,3) * t410 + t479 * t501;
t257 = -rSges(4,3) * t478 + t477 * t501;
t256 = Icges(4,5) * t410 + t479 * t497;
t255 = -Icges(4,5) * t478 + t477 * t497;
t254 = Icges(4,6) * t410 + t479 * t496;
t253 = -Icges(4,6) * t478 + t477 * t496;
t250 = rSges(5,3) * t403 + t402 * t518;
t243 = Icges(5,5) * t403 + t402 * t516;
t242 = Icges(5,6) * t403 + t402 * t515;
t239 = pkin(3) * t484 + pkin(8) * t368;
t236 = -pkin(3) * t363 + pkin(8) * t364;
t235 = rSges(4,1) * t484 - rSges(4,2) * t368;
t234 = -rSges(4,1) * t363 - rSges(4,2) * t364;
t227 = rSges(5,1) * t358 + rSges(5,2) * t498;
t202 = t615 * t215;
t184 = rSges(3,3) * t553 - t521;
t176 = t298 * t441 - t692;
t171 = rSges(5,3) * t368 - t484 * t518;
t170 = rSges(5,3) * t364 + t363 * t518;
t169 = Icges(5,5) * t368 - t484 * t516;
t168 = Icges(5,5) * t364 + t363 * t516;
t167 = Icges(5,6) * t368 - t484 * t515;
t166 = Icges(5,6) * t364 + t363 * t515;
t162 = t378 * t592 + t379 * t479 + t380 * t410;
t160 = rSges(5,1) * t311 + rSges(5,2) * t310 + rSges(5,3) * t594;
t159 = rSges(5,1) * t309 + rSges(5,2) * t308 + rSges(5,3) * t596;
t158 = Icges(5,1) * t311 + Icges(5,4) * t310 + Icges(5,5) * t594;
t157 = Icges(5,1) * t309 + Icges(5,4) * t308 + Icges(5,5) * t596;
t156 = Icges(5,4) * t311 + Icges(5,2) * t310 + Icges(5,6) * t594;
t155 = Icges(5,4) * t309 + Icges(5,2) * t308 + Icges(5,6) * t596;
t154 = Icges(5,5) * t311 + Icges(5,6) * t310 + Icges(5,3) * t594;
t153 = Icges(5,5) * t309 + Icges(5,6) * t308 + Icges(5,3) * t596;
t152 = rSges(5,1) * t274 - rSges(5,2) * t275;
t151 = -rSges(5,1) * t272 - rSges(5,2) * t271;
t138 = t162 * t441;
t137 = t615 * t384 + (t385 * t453 + t386 * t450 + (-t379 * t450 + t380 * t453) * qJD(2)) * t447;
t134 = t137 * t441;
t119 = qJD(1) * t692 - t184 * t441 - t387 * t436;
t118 = qJD(2) * t472 + t183 * t441 + t401;
t79 = -t330 * t379 + t331 * t380 + t385 * t477 - t386 * t478 + (t378 * t573 - t384 * t454) * t447;
t78 = t328 * t379 + t329 * t380 + t385 * t479 + t386 * t410 + (t378 * t572 + t384 * t451) * t447;
t77 = t214 * t414 - t284 * t370 + t492;
t76 = t467 + t699;
t73 = t212 * t370 + t214 * t371 + t579;
t52 = qJD(2) * t487 + t138;
t51 = qJD(2) * t488 - t703;
t35 = -t222 * t371 + t278 * t284 - t414 * t90 + (-t212 * t565 - t598) * t571 + t493;
t34 = -t222 * t370 - t279 * t284 + t414 * t89 + (t214 * t565 + t494) * t571 + t585;
t32 = t108 * t414 + t370 * t75 - t371 * t74;
t27 = t212 * t279 - t214 * t278 + t370 * t90 + t371 * t89 + t468;
t18 = t262 * t49 - t263 * t48 + t348 * t68;
t4 = [(t64 + t78) * t528 + (t63 + t79 + t52) * t527 - t668 * t649 + (t121 + t162) * t506 + t17 * t651 + (t15 * (t363 * t634 + t489 + t519 - t630) + t42 * (t196 * t634 + t466 + t520 - t631) + t14 * (t508 + t586) + (t42 + t617 + t660 - t709) * t43) * m(5) + (t138 + ((t109 + t693) * t451 + (t110 + (t599 - t600) * t447 - t111 + t582) * t454) * t571) * t526 + t625 / 0.2e1 - t626 / 0.2e1 + (t35 * (-t212 + t489) + t76 * (t466 - t90) + t34 * (t508 + t214) + (t660 - t699 + t76 + t89) * t77) * m(4) - t697 * t655 + t621 / 0.2e1 + t622 / 0.2e1 + (((t582 + t691) * t451 - t693 * t454 + ((t599 + t600) * t451 + t286 * t454 ^ 2) * t447 + t513) * t571 + t51 + t703) * t529 + t134 + t618 + t619 / 0.2e1 + t620 / 0.2e1 + (t119 * (-t296 - t418) + t175 * t521 + t118 * (t298 + t661) + (-pkin(1) * t603 + t175 * (-rSges(3,3) - pkin(6)) * t592) * qJD(1) + (t183 + t523) * t176) * m(3) + t628 + t31 * t640 + (t120 - t469) * t507 + t37 * t639 + t36 * t641 + t82 * t648 + t20 * t650 + t19 * t652 + t66 * t654 + t371 * t656 + t263 * t657 + t623 / 0.2e1 - t624 / 0.2e1; (-t42 * (-t131 * t383 - t159 * t348 + t201 * t315 - t263 * t265 - t313 * t414 - t371 * t381 + t480) - t43 * (t133 * t383 + t160 * t348 - t201 * t316 - t262 * t265 + t314 * t414 - t370 * t381 + t502) - t33 * (t131 * t316 - t133 * t315 + t159 * t262 + t160 * t263 + t313 * t370 + t314 * t371 + t580) - (t42 * (t238 * t593 - t334 * t478) + t43 * (t240 * t593 - t334 * t410) + t33 * (-t238 * t410 + t240 * t478)) * qJD(3) - t15 * t538 + t42 * t517 + t14 * t317 + t43 * t202 + t7 * t577 + t33 * t560 + ((t15 * t559 + t42 * t561 + (t43 * t559 + t665) * qJD(1) + t659) * t454 + (t14 * t559 + t43 * t561 + (t541 + t33 * (-t347 - t586)) * qJD(1) + t476) * t451) * t447 + t669 * t615) * m(5) + (-t668 * t615 + (t451 * t70 - t454 * t69) * t447) * t649 + ((-t205 * t478 - t254 * t363 + t256 * t364) * t370 - (-t203 * t478 - t253 * t363 + t255 * t364) * t371 + (-t281 * t478 - t361 * t363 + t362 * t364) * t414 + (t410 * t70 - t478 * t69 - t593 * t668) * qJD(3) + t658 * t477) * t640 + ((t124 * t563 + t127 * t376 + t130 * t377 + t154 * t402 + t156 * t358 - t158 * t498) * t262 + t49 * t316 - (t122 * t563 - t126 * t376 - t129 * t377 + t153 * t402 + t155 * t358 - t157 * t498) * t263 + t48 * t315 + (t198 * t563 + t199 * t376 + t200 * t377 + t259 * t402 + t260 * t358 - t261 * t498) * t348 + t68 * t383) * t644 + ((t124 * t594 + t127 * t310 + t130 * t311 - t154 * t484 + t156 * t274 + t158 * t275) * t262 + t41 * t316 - (t122 * t594 - t126 * t310 - t129 * t311 - t153 * t484 + t155 * t274 + t157 * t275) * t263 + t40 * t315 + (t198 * t594 + t199 * t310 + t200 * t311 - t259 * t484 + t260 * t274 + t261 * t275) * t348 + t66 * t383) * t653 + ((t124 * t596 + t127 * t308 + t130 * t309 + t154 * t363 - t156 * t272 + t158 * t271) * t262 + t39 * t316 - (t122 * t596 - t126 * t308 - t129 * t309 + t153 * t363 - t155 * t272 + t157 * t271) * t263 + t38 * t315 + (t198 * t596 + t199 * t308 + t200 * t309 + t259 * t363 - t260 * t272 + t261 * t271) * t348 - t697 * t383) * t651 + ((t404 * t592 + t410 * t576 + t479 * t575) * t441 + (t410 * t471 - t470 * t479 - t509 * t592) * t571) * t529 + (-t697 * t615 + (-t38 * t454 + t39 * t451) * t447) * t655 - t383 * t18 / 0.2e1 + (t615 * t79 + t459) * t527 + (t615 * t78 + t460) * t528 + ((-t404 * t590 + t477 * t575 - t478 * t576) * t441 + (-t470 * t477 - t471 * t478 + t509 * t590) * t571) * t526 - t32 * t549 / 0.2e1 + t30 * t567 / 0.2e1 + (t108 * t615 + (t451 * t75 - t454 * t74) * t447) * t505 + (t162 * t615 + t487) * t506 - t441 * (t615 * t404 * t441 + ((t450 * t576 + t453 * t575) * t441 + (t336 * t540 - t335 * t539 + ((t338 * t453 + t340 * t450 + t510) * t451 - (t337 * t453 + t339 * t450 + t511) * t454) * t447) * qJD(2)) * t447) / 0.2e1 + t441 * (t137 * t615 + t458) / 0.2e1 + (-t469 * t615 + t488) * t507 + (qJD(2) * t458 + t134 + t687) * t615 / 0.2e1 - (qJD(2) * t459 + t441 * t79 + t688) * t590 / 0.2e1 + (qJD(2) * t460 + t441 * t78 + t689) * t592 / 0.2e1 + ((t52 + t684) * t454 + (t51 + t685) * t451) * qJD(1) * t447 / 0.2e1 + (t53 * t615 + (-t28 * t454 + t29 * t451 + (t451 * t74 + t454 * t75) * qJD(1)) * t447) * t635 + (t37 * t615 + (-t23 * t454 + t24 * t451 + (t451 * t69 + t454 * t70) * qJD(1)) * t447) * t639 + (t36 * t615 + (-t21 * t454 + t22 * t451 + (t451 * t71 + t454 * t72) * qJD(1)) * t447) * t641 + (t26 * t615 + (-t12 * t454 + t13 * t451 + (t451 * t48 + t454 * t49) * qJD(1)) * t447) * t643 + (t68 * t615 + (t451 * t49 - t454 * t48) * t447) * t647 + (t82 * t615 + (t451 * t72 - t454 * t71) * t447) * t648 + (t20 * t615 + (-t10 * t454 + t11 * t451 + (t38 * t451 + t39 * t454) * qJD(1)) * t447) * t650 + (t19 * t615 + (t451 * t9 - t454 * t8 + (t40 * t451 + t41 * t454) * qJD(1)) * t447) * t652 + (t66 * t615 + (-t40 * t454 + t41 * t451) * t447) * t654 + t566 * t656 + t316 * t657 + (t119 * (-t296 * t615 - t382 * t590) - t175 * t615 * t184 + t118 * (t298 * t615 - t382 * t592) + t176 * (t183 * t615 + t472) + (t175 * (t382 * t573 - t387 * t454) + 0.2e1 * t163 * (t183 * t454 + t184 * t451 + (t296 * t454 - t298 * t451) * qJD(1))) * t447 - (-t175 * t341 + t176 * t342) * t441 - (t163 * (t341 * t451 + t342 * t454) + (-t176 * t451 - t603) * t491) * t571) * m(3) + ((-t254 * t402 + t256 * t403) * t370 - (-t253 * t402 + t255 * t403) * t371 + (-t361 * t402 + t362 * t403) * t414 + (t410 * t75 - t478 * t74) * qJD(3) + ((qJD(3) * t108 - t203 * t371 + t205 * t370 + t281 * t414) * t450 + t658 * t453) * t447) * t636 + (t35 * (-t212 * t615 + t447 * t537 - t538) + t34 * (t214 * t615 + t581 * t592 + t317) + t27 * ((t212 * t451 + t214 * t454) * t447 + t577) + (-t258 * t414 + t369 * t370 - t502 - (t214 * t593 - t284 * t410) * qJD(3) + t615 * t89 + t202 + (qJD(1) * t537 + t451 * t583) * t447) * t77 + (t257 * t414 + t369 * t371 - t480 - (-t212 * t593 - t284 * t478) * qJD(3) - t615 * t90 + (t284 * t573 + t454 * t583) * t447 + t517) * t76 + (-t257 * t370 - t258 * t371 - t580 - (t212 * t410 + t214 * t478) * qJD(3) + (t451 * t90 + t454 * t89 + (t212 * t454 + (-t214 - t347) * t451) * qJD(1)) * t447 + t560) * t73) * m(4) + ((t205 * t410 + t254 * t484 + t256 * t368) * t370 - (t203 * t410 + t253 * t484 + t255 * t368) * t371 + (t281 * t410 + t361 * t484 + t362 * t368) * t414 + (t410 * t72 - t478 * t71 + t593 * t82) * qJD(3) + t658 * t479) * t642 - t315 * t16 / 0.2e1; ((t124 * t403 + t167 * t358 - t169 * t498) * t262 - (t122 * t403 + t166 * t358 - t168 * t498) * t263 + (t198 * t403 + t242 * t358 - t243 * t498) * t348 + (t364 * t48 + t368 * t49 + t403 * t68) * qJD(4) + t457 * t402) * t644 + ((t124 * t368 + t167 * t274 + t169 * t275) * t262 - (t122 * t368 + t166 * t274 + t168 * t275) * t263 + (t198 * t368 + t242 * t274 + t243 * t275) * t348 + (t364 * t40 + t368 * t41 + t403 * t66) * qJD(4) - t457 * t484) * t653 + ((-t33 * t586 + t541) * t330 + (t43 * t584 - t665) * t328 - (-t14 * t584 - t43 * t588 + t476) * t479 - (t15 * t584 + t42 * t588 - t659) * t477 + ((-t42 * t587 + t43 * t586) * t570 - t669 * t453) * t447 - t42 * (-t170 * t348 - t236 * t414 - t250 * t263 - t333 * t371) - t43 * (t171 * t348 + t239 * t414 - t250 * t262 - t333 * t370) - t33 * (t170 * t262 + t171 * t263 + t236 * t370 + t239 * t371) - (t42 * (-t131 * t403 + t201 * t364) + t43 * (t133 * t403 - t201 * t368) + t33 * (t131 * t368 - t133 * t364)) * qJD(4)) * m(5) - (t16 * t364 + t17 * t368 + t18 * t403) * qJD(4) / 0.2e1 + (t32 + t18) * t530 + (-t23 * t477 - t24 * t479 - t328 * t70 + t330 * t69 + (-t37 * t453 - t570 * t668) * t447) * t639 + (-t477 * t69 - t479 * t70 + t591 * t668) * t649 + (t363 * t462 + t364 * t463 - t474 * t477) * t640 + (-t21 * t477 - t22 * t479 - t328 * t72 + t330 * t71 + (-t36 * t453 + t570 * t82) * t447) * t641 + (-t12 * t477 - t13 * t479 - t328 * t49 + t330 * t48 + (-t26 * t453 + t570 * t68) * t447) * t643 + (-t477 * t48 - t479 * t49 - t591 * t68) * t647 + (-t477 * t71 - t479 * t72 - t591 * t82) * t648 + (-t328 * t41 + t330 * t40 - t477 * t8 - t479 * t9 + (-t19 * t453 + t570 * t66) * t447) * t652 + (-t40 * t477 - t41 * t479 - t591 * t66) * t654 + (t368 * t463 - t462 * t484 - t474 * t479) * t642 + ((t124 * t364 - t167 * t272 + t169 * t271) * t262 - (t122 * t364 - t166 * t272 + t168 * t271) * t263 + (t198 * t364 - t242 * t272 + t243 * t271) * t348 + (t364 * t38 + t368 * t39 - t403 * t697) * qJD(4) + t457 * t363) * t651 + (-t10 * t477 - t11 * t479 - t328 * t39 + t330 * t38 + (-t20 * t453 - t570 * t697) * t447) * t650 + (-t38 * t477 - t39 * t479 + t591 * t697) * t655 + (-t108 * t591 - t477 * t74 - t479 * t75) * t505 + (-t28 * t477 - t29 * t479 - t328 * t75 + t330 * t74 + (t108 * t570 - t453 * t53) * t447) * t635 - t684 * t328 / 0.2e1 + t685 * t330 / 0.2e1 - t687 * t591 / 0.2e1 - t688 * t477 / 0.2e1 - t689 * t479 / 0.2e1 + (t402 * t462 + t403 * t463 - t474 * t591) * t636 + (t27 * (-t212 * t479 + t214 * t477) + (-t477 * t76 + t479 * t77) * t222 + (t328 * t77 + t330 * t76 + t34 * t479 - t35 * t477) * t284 + ((-t212 * t76 + t214 * t77) * t570 + (t212 * t35 - t214 * t34 + t76 * t90 - t77 * t89) * t453) * t447 - t76 * (-t234 * t414 - t332 * t371) - t77 * (t235 * t414 - t332 * t370) + (-t212 * t328 - t214 * t330 - t234 * t370 - t235 * t371 + t477 * t89 - t479 * t90) * t73) * m(4); t194 * t17 / 0.2e1 - t484 * t1 / 0.2e1 + (t363 * t40 + t402 * t66 - t41 * t484) * t654 + (t19 * t402 + t194 * t41 + t196 * t40 + t356 * t66 + t363 * t8 - t484 * t9) * t652 + t196 * t16 / 0.2e1 + t363 * t2 / 0.2e1 + (t363 * t38 - t39 * t484 - t402 * t697) * t655 + (t10 * t363 - t11 * t484 + t194 * t39 + t196 * t38 + t20 * t402 - t356 * t697) * t650 + t356 * t18 / 0.2e1 + t402 * t3 / 0.2e1 + (t363 * t48 + t402 * t68 - t484 * t49) * t647 + (t12 * t363 - t13 * t484 + t194 * t49 + t196 * t48 + t26 * t402 + t356 * t68) * t643 + (t274 * t464 + t275 * t465 - t475 * t484) * t653 + (t271 * t465 - t272 * t464 + t363 * t475) * t651 + (t358 * t464 + t402 * t475 - t465 * t498) * t644 + (t15 * (-t131 * t402 + t201 * t363) + t14 * (t133 * t402 + t201 * t484) + t7 * (-t131 * t484 - t133 * t363) + (t100 * t484 + t133 * t356 - t152 * t348 - t194 * t201 + t227 * t262 + t402 * t60) * t43 + (t100 * t363 - t131 * t356 + t151 * t348 + t196 * t201 + t227 * t263 - t402 * t61) * t42 + (t131 * t194 - t133 * t196 - t151 * t262 - t152 * t263 - t363 * t60 - t484 * t61) * t33) * m(5);];
tauc = t4(:);