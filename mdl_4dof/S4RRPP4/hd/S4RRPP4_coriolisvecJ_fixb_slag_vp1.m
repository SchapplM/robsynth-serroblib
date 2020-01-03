% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP4
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
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:59:15
% DurationCPUTime: 18.89s
% Computational Cost: add. (4231->505), mult. (11193->616), div. (0->0), fcn. (8700->4), ass. (0->284)
t630 = Icges(4,1) + Icges(5,1);
t629 = Icges(5,2) + Icges(4,3);
t628 = Icges(4,6) - Icges(5,6);
t627 = Icges(4,4) - Icges(5,5);
t299 = cos(qJ(2));
t287 = Icges(3,4) * t299;
t297 = sin(qJ(2));
t500 = Icges(4,5) * t299;
t502 = Icges(5,4) * t299;
t601 = t630 * t297 - t500 - t502;
t621 = Icges(3,1) * t297 + t287 + t601;
t285 = Icges(5,4) * t297;
t369 = Icges(5,1) * t299 + t285;
t283 = Icges(4,5) * t297;
t370 = Icges(4,1) * t299 + t283;
t626 = t369 + t370;
t298 = sin(qJ(1));
t473 = t297 * t298;
t265 = Icges(3,4) * t473;
t471 = t298 * t299;
t300 = cos(qJ(1));
t501 = Icges(3,5) * t300;
t145 = Icges(3,1) * t471 - t265 - t501;
t141 = Icges(5,5) * t300 + t298 * t369;
t143 = -Icges(4,4) * t300 + t298 * t370;
t602 = t141 + t143;
t595 = t145 + t602;
t504 = Icges(3,4) * t297;
t234 = Icges(3,1) * t299 - t504;
t146 = Icges(3,5) * t298 + t234 * t300;
t619 = t627 * t298 + t626 * t300;
t594 = t146 + t619;
t600 = -t629 * t299 + t283 + t285;
t593 = -Icges(3,2) * t299 - t504 + t600;
t470 = t299 * t300;
t625 = (Icges(5,4) + Icges(4,5)) * t470;
t368 = -Icges(3,2) * t297 + t287;
t220 = Icges(4,3) * t297 + t500;
t224 = Icges(5,2) * t297 + t502;
t612 = t220 + t224;
t624 = -t368 + t612;
t623 = t628 * t298;
t131 = -Icges(4,6) * t300 + t220 * t298;
t135 = Icges(5,6) * t300 + t224 * t298;
t603 = -t131 - t135;
t472 = t297 * t300;
t613 = t629 * t472 + t623 + t625;
t222 = Icges(3,5) * t299 - Icges(3,6) * t297;
t134 = Icges(3,3) * t298 + t222 * t300;
t226 = Icges(4,4) * t299 + Icges(4,6) * t297;
t138 = Icges(4,2) * t298 + t226 * t300;
t622 = t134 + t138;
t606 = (Icges(3,6) - t628) * t299 + (Icges(3,5) + t627) * t297;
t620 = t234 + t626;
t421 = qJD(2) * t297;
t419 = qJD(2) * t299;
t617 = -t593 * t297 - t621 * t299;
t490 = Icges(3,3) * t300;
t133 = Icges(3,5) * t471 - Icges(3,6) * t473 - t490;
t495 = Icges(3,6) * t300;
t139 = Icges(3,4) * t471 - Icges(3,2) * t473 - t495;
t481 = t139 * t297;
t357 = -t145 * t299 + t481;
t360 = t135 * t297 + t141 * t299;
t137 = -Icges(4,2) * t300 + t226 * t298;
t482 = t137 * t300;
t218 = Icges(5,5) * t299 + Icges(5,6) * t297;
t129 = Icges(5,3) * t300 + t218 * t298;
t484 = t129 * t300;
t364 = t131 * t297 + t143 * t299;
t553 = t298 * t364;
t574 = t298 * t360 - t482 + t484 + t553;
t616 = -t133 * t300 - t298 * t357 + t574;
t140 = Icges(3,6) * t298 + t300 * t368;
t112 = t146 * t471;
t389 = t134 * t300 - t112;
t50 = -t140 * t473 - t389;
t130 = -Icges(5,3) * t298 + t218 * t300;
t127 = t300 * t130;
t573 = -t138 * t300 + t619 * t471 + t613 * t473 + t127;
t615 = t50 + t573;
t614 = t594 * t297;
t611 = t622 * t298 + t594 * t470 + t613 * t472;
t125 = t298 * t137;
t610 = -t298 * t133 - t595 * t470 + t603 * t472 - t125;
t609 = t624 * qJD(2);
t608 = t620 * qJD(2);
t607 = t218 - t222 - t226;
t543 = t606 * t300;
t544 = t606 * t298;
t604 = -t593 * t419 + t621 * t421;
t485 = t129 * t298;
t564 = -t139 * t472 - t485 - t610;
t563 = -t130 * t298 - t140 * t472 + t611;
t599 = t617 * t298 + t543;
t598 = -t617 * t300 + t544;
t597 = t139 + t603;
t596 = t140 - t613;
t509 = -rSges(5,3) - qJ(4);
t590 = t608 * t299 + t609 * t297 + (-t297 * t621 + t593 * t299) * qJD(2) + t606 * qJD(1);
t589 = t606 * qJD(2);
t480 = t140 * t297;
t588 = t297 * t613 + t299 * t594 - t480;
t587 = t357 - t360 - t364;
t586 = -qJD(1) * t617 + qJD(2) * t607;
t585 = (t300 * t593 + t594) * t298 - (-Icges(3,2) * t471 + t298 * t600 - t265 + t595) * t300;
t584 = (-t298 * t589 + (-t130 + t587 + t622) * qJD(1)) * t300;
t422 = qJD(1) * t300;
t583 = t598 * qJD(1);
t582 = t597 * t300 + (t601 * t300 - t630 * t472 - t596 + t625) * t298;
t581 = (t298 * t563 - t300 * t564) * qJD(2);
t580 = (t615 * t298 - t300 * t616) * qJD(2);
t579 = t621 - t624;
t578 = t593 + t620;
t575 = t599 * qJD(1);
t416 = qJD(3) * t300;
t258 = t297 * t416;
t235 = pkin(2) * t297 - qJ(3) * t299;
t236 = rSges(5,1) * t297 - rSges(5,2) * t299;
t523 = pkin(3) * t297;
t392 = t236 + t523;
t382 = -t235 - t392;
t415 = qJD(4) * t298;
t418 = qJD(2) * t300;
t322 = t382 * t418 - t415;
t572 = t258 + t322;
t272 = pkin(3) * t470;
t571 = t298 * t509 + t272;
t570 = -t575 + t580;
t569 = t581 + t583;
t568 = t587 * qJD(2) + t604 * t298 + ((t300 * t612 - t140 + t623) * t299 - t614) * qJD(1);
t567 = t588 * qJD(2) - t604 * t300 + ((-t298 * t368 + t495 - t603) * t299 + (-t234 * t298 + t501 - t602) * t297) * qJD(1);
t566 = -t298 * t586 + t300 * t590;
t565 = t298 * t590 + t300 * t586;
t562 = t299 * t596 + t614;
t561 = t297 * t595 + t299 * t597;
t558 = t482 + t611;
t556 = -t589 * t300 + (-t222 * t298 + t129 - t137 + t490 - t588) * qJD(1);
t412 = -rSges(5,1) - pkin(2) - pkin(3);
t487 = qJ(3) * t297;
t518 = rSges(5,2) * t297;
t555 = qJD(1) * (t412 * t299 - pkin(1) - t487 - t518) - qJD(4);
t554 = 0.2e1 * qJD(2);
t188 = t236 * t298;
t514 = t297 * rSges(4,1);
t237 = -rSges(4,3) * t299 + t514;
t278 = qJD(3) * t297;
t552 = (qJD(2) * t237 - t278) * t298;
t240 = pkin(2) * t299 + t487;
t191 = t240 * t298;
t163 = qJD(1) * t191;
t293 = t300 * pkin(5);
t245 = pkin(1) * t298 - t293;
t216 = qJD(1) * t245;
t551 = -t163 - t216;
t246 = t300 * pkin(1) + t298 * pkin(5);
t289 = t298 * rSges(4,2);
t158 = rSges(4,1) * t470 + rSges(4,3) * t472 + t289;
t273 = pkin(2) * t470;
t196 = qJ(3) * t472 + t273;
t383 = t196 + t246;
t548 = t383 + t158;
t241 = rSges(5,1) * t299 + t518;
t533 = -rSges(5,3) * t300 - t241 * t298;
t453 = pkin(3) * t471 + qJ(4) * t300 - t533;
t545 = qJD(1) * t453;
t542 = -t297 * t585 + t582 * t299;
t541 = (-t579 * t297 + t578 * t299) * qJD(1);
t540 = t607 * qJD(1);
t417 = qJD(3) * t299;
t162 = qJD(2) * t240 - t417;
t448 = -t241 * qJD(2) - t162;
t522 = pkin(3) * t299;
t529 = -pkin(3) * t419 - qJD(2) * (-t240 - t241 - t522) + t448;
t528 = m(4) / 0.2e1;
t527 = m(5) / 0.2e1;
t524 = -rSges(4,1) - pkin(2);
t520 = qJD(1) / 0.2e1;
t519 = rSges(3,1) * t299;
t238 = rSges(3,1) * t297 + rSges(3,2) * t299;
t195 = t238 * t300;
t420 = qJD(2) * t298;
t406 = t238 * t420;
t288 = t298 * rSges(3,3);
t159 = rSges(3,1) * t470 - rSges(3,2) * t472 + t288;
t450 = t159 + t246;
t62 = qJD(1) * t450 - t406;
t515 = t195 * t62;
t432 = rSges(3,2) * t473 + t300 * rSges(3,3);
t156 = rSges(3,1) * t471 - t432;
t405 = t238 * t418;
t61 = -t405 + (-t156 - t245) * qJD(1);
t513 = t298 * t61;
t512 = t300 * t61;
t511 = -rSges(5,2) - qJ(3);
t510 = -rSges(4,3) - qJ(3);
t215 = t246 * qJD(1);
t167 = t297 * t422 + t298 * t419;
t404 = t297 * t420;
t255 = pkin(2) * t404;
t81 = qJ(3) * t167 + qJD(1) * t273 + t278 * t298 - t255;
t508 = -t215 - t81;
t402 = t299 * t418;
t252 = rSges(5,2) * t402;
t403 = t297 * t418;
t423 = qJD(1) * t298;
t330 = -t299 * t423 - t403;
t466 = -rSges(5,1) * t403 + pkin(3) * t330 - qJ(4) * t422 + qJD(1) * t533 + t252 - t415;
t254 = pkin(3) * t404;
t414 = qJD(4) * t300;
t386 = -t254 + t414;
t465 = -qJD(2) * t188 + t386 + (t241 * t300 + t571) * qJD(1);
t431 = rSges(5,1) * t470 + rSges(5,2) * t472;
t452 = t431 + t571;
t451 = -t158 - t196;
t449 = t298 * t191 + t300 * t196;
t242 = rSges(4,1) * t299 + rSges(4,3) * t297;
t447 = -t242 * qJD(2) - t162;
t192 = t235 * t300;
t446 = -qJD(1) * t192 + t298 * t417;
t197 = t235 * t420;
t413 = qJD(2) * qJD(3);
t399 = t299 * t413;
t445 = qJD(1) * t197 + t300 * t399;
t444 = -t191 - t245;
t437 = -t235 - t237;
t436 = -t240 - t242;
t435 = qJ(3) * t402 + t258;
t434 = rSges(4,2) * t422 + rSges(4,3) * t402;
t407 = t297 * t423;
t433 = rSges(3,2) * t407 + rSges(3,3) * t422;
t430 = t298 ^ 2 + t300 ^ 2;
t80 = pkin(2) * t330 - qJ(3) * t407 + t435;
t411 = t191 * t422 + t298 * t81 + t300 * t80;
t187 = t235 * t298;
t410 = -t187 * t420 - t192 * t418 + t278;
t409 = -t196 - t452;
t277 = pkin(5) * t422;
t408 = t277 + t435;
t400 = -pkin(1) - t519;
t396 = -t420 / 0.2e1;
t393 = t418 / 0.2e1;
t390 = t300 * t437;
t388 = -t133 + t480;
t387 = t297 * t413 + t81 * t420 + (t163 + t80) * t418;
t199 = qJD(1) * (-pkin(1) * t423 + t277);
t384 = t298 * t399 + t199 + (t258 + t80) * qJD(1);
t377 = t299 * t524 - pkin(1);
t375 = t191 * t420 + t196 * t418 - t417;
t373 = -rSges(3,2) * t297 + t519;
t372 = -t298 * t62 - t512;
t348 = qJD(2) * t390 + t258;
t190 = t238 * t298;
t189 = t237 * t298;
t57 = (t156 * t298 + t159 * t300) * qJD(2);
t329 = (-t522 * qJD(2) + t448) * qJD(2);
t291 = t300 * rSges(4,2);
t259 = t299 * t416;
t211 = t373 * qJD(2);
t198 = t235 * t423;
t194 = t237 * t300;
t193 = t236 * t300;
t168 = t402 - t407;
t166 = t430 * t421;
t155 = t242 * t298 - t291;
t105 = -qJD(2) * t190 + (t300 * t373 + t288) * qJD(1);
t104 = -qJD(2) * t189 + (t242 * t300 + t289) * qJD(1);
t102 = rSges(3,1) * t330 - rSges(3,2) * t402 + t433;
t101 = rSges(4,1) * t330 - rSges(4,3) * t407 + t434;
t44 = qJD(1) * t548 - t197 - t552;
t43 = (-t155 + t444) * qJD(1) + t348;
t42 = -t211 * t418 + (-t105 - t215 + t406) * qJD(1);
t41 = -t211 * t420 + t199 + (t102 - t405) * qJD(1);
t40 = (t155 * t298 + t158 * t300) * qJD(2) + t375;
t39 = -t197 + (-qJD(2) * t236 + t278) * t298 + (t246 - t409) * qJD(1) + t386;
t38 = (t444 - t453) * qJD(1) + t572;
t37 = (t298 * t453 + t300 * t452) * qJD(2) + t375;
t24 = t447 * t418 + (-t104 + t508 + t552) * qJD(1) + t445;
t23 = qJD(1) * t101 + (qJD(1) * t390 + t298 * t447) * qJD(2) + t384;
t16 = t329 * t300 + (-t414 + (qJD(2) * t392 - t278) * t298 - t465 + t508) * qJD(1) + t445;
t15 = t329 * t298 + (t322 + t466) * qJD(1) + t384;
t2 = (t101 * t300 + t104 * t298 + (t155 * t300 + t298 * t451) * qJD(1)) * qJD(2) + t387;
t1 = (t466 * t300 + t465 * t298 + (t298 * t409 + t300 * t453) * qJD(1)) * qJD(2) + t387;
t3 = [(-qJD(2) * t617 + t608 * t297 - t609 * t299) * qJD(1) + (t16 * t293 + t38 * (t254 + t255) + t15 * (t272 + t383 + t431) + (t16 * t509 + t555 * t38) * t300 + (-t16 * pkin(1) + t15 * t509 + (t16 * t511 + t38 * (rSges(5,1) * qJD(2) - qJD(3))) * t297 + (qJD(2) * t38 * t511 + t16 * t412) * t299 + t38 * (-pkin(5) - t509) * qJD(1)) * t298 + (t38 + t545 - t551 + t252 + t408 + (qJD(1) * t509 + t412 * t421) * t300 + t555 * t298 - t572) * t39) * m(5) + (-(-qJD(1) * t155 + t348 - t43 + t551) * t44 + t24 * (t291 + t293) + t43 * t255 + t23 * t548 + t44 * (t408 + t434) + (t44 * t524 * t421 + t43 * (t297 * t510 + t377) * qJD(1)) * t300 + (t24 * t377 + (-t43 * qJD(3) + t24 * t510) * t297 + t43 * (t299 * t510 + t514) * qJD(2) + (t43 * (-rSges(4,2) - pkin(5)) + t44 * (-pkin(1) + t436)) * qJD(1)) * t298) * m(4) + (-(-qJD(1) * t156 - t216 - t405 - t61) * t62 + t42 * (t298 * t400 + t293 + t432) + t41 * t450 + t62 * (t277 + t433) + (t238 * t513 - t515) * qJD(2) + ((-pkin(1) - t373) * t512 + (t61 * (-rSges(3,3) - pkin(5)) + t62 * t400) * t298) * qJD(1)) * m(3) + (((t50 - t112 + (t134 + t481) * t300 + t610) * t300 + (-t553 + (-t130 - t360) * t298 + t558 + t574) * t298) * qJD(2) + t583) * t393 + (((t300 * t388 + t484 - t558 + t563) * t300 + (t298 * t388 - t125 + t127 + t389 + t485 + t564 - t573) * t298) * qJD(2) + t570 + t575) * t396 + (t566 + t567) * t420 / 0.2e1 - (t565 - t568 + t569) * t418 / 0.2e1 + ((t561 - t599) * t298 + (t562 + t598) * t300) * qJD(2) * t520; (t1 * t449 + (t1 * t453 + t15 * t382) * t298 + (t1 * t452 + t16 * t382) * t300 + (-t446 + t529 * t298 + (pkin(3) * t472 + t300 * t382 + t193) * qJD(1)) * t39 + (-qJD(1) * t187 + t300 * t529 + t198 - t259) * t38 + (-t410 - (-t188 * t298 - t193 * t300 - t430 * t523) * qJD(2) + t411 + (qJD(1) * t409 + t465) * t298 + (t466 + t545) * t300) * t37) * m(5) + (t43 * t198 + t2 * t449 + t40 * t411 + (t24 * t437 + t43 * t447 + t2 * t158 + t40 * t101 + (t40 * t155 + t437 * t44) * qJD(1)) * t300 + (t23 * t437 + t44 * t447 + t2 * t155 + t40 * t104 + (t43 * t237 + t40 * t451) * qJD(1)) * t298 - t43 * (t259 + (t187 + t189) * qJD(1)) - t44 * (-qJD(1) * t194 + t446) - t40 * t410 - ((-t40 * t194 + t43 * t436) * t300 + (-t40 * t189 + t436 * t44) * t298) * qJD(2)) * m(4) + (-(t190 * t61 - t515) * qJD(1) - (t57 * (-t190 * t298 - t195 * t300) + t372 * t373) * qJD(2) + 0.2e1 * t57 * (t102 * t300 + t105 * t298 + (t156 * t300 - t159 * t298) * qJD(1)) + t372 * t211 + (-t41 * t298 - t42 * t300 + (-t300 * t62 + t513) * qJD(1)) * t238) * m(3) - ((t582 * t297 + t299 * t585) * qJD(2) + (t578 * t297 + t579 * t299) * qJD(1)) * qJD(1) / 0.2e1 + (t568 * t300 + t567 * t298 + (t298 * t561 + t300 * t562) * qJD(1)) * t520 + ((-t420 * t543 - t540) * t298 + ((t298 * t544 + t542) * qJD(2) + t541) * t300) * t396 + ((-t418 * t544 + t540) * t300 + ((t300 * t543 + t542) * qJD(2) + t541) * t298) * t393 + (t566 * qJD(1) + (t563 * t422 + (t564 * qJD(1) + t556 * t298 - t584) * t298) * t554) * t298 / 0.2e1 - (t565 * qJD(1) + ((qJD(1) * t615 + t584) * t300 + (qJD(1) * t616 - t556 * t300) * t298) * t554) * t300 / 0.2e1 + (t570 + t580) * t423 / 0.2e1 + (t569 + t581) * t422 / 0.2e1; -m(4) * (t166 * t40 + t167 * t44 + t168 * t43) - m(5) * (t166 * t37 + t167 * t39 + t168 * t38) + 0.2e1 * ((t418 * t43 + t420 * t44 - t2) * t528 + (t38 * t418 + t39 * t420 - t1) * t527) * t299 + 0.2e1 * ((qJD(2) * t40 + t23 * t298 + t24 * t300 + t422 * t44 - t423 * t43) * t528 + (qJD(2) * t37 + t15 * t298 + t16 * t300 - t38 * t423 + t39 * t422) * t527) * t297; m(5) * (t15 * t300 - t16 * t298);];
tauc = t3(:);
