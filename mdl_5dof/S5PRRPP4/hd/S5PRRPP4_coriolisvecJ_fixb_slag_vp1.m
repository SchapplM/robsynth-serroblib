% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:47
% EndTime: 2019-12-31 17:41:09
% DurationCPUTime: 19.03s
% Computational Cost: add. (9199->516), mult. (11499->628), div. (0->0), fcn. (8884->4), ass. (0->291)
t634 = Icges(5,1) + Icges(6,1);
t633 = Icges(6,2) + Icges(5,3);
t632 = Icges(5,6) - Icges(6,6);
t631 = Icges(5,4) - Icges(6,5);
t301 = sin(qJ(3));
t298 = Icges(6,4) * t301;
t302 = cos(qJ(3));
t370 = Icges(6,1) * t302 + t298;
t297 = Icges(5,5) * t301;
t371 = Icges(5,1) * t302 + t297;
t630 = t370 + t371;
t300 = pkin(7) + qJ(2);
t294 = sin(t300);
t479 = t294 * t301;
t266 = Icges(4,4) * t479;
t478 = t294 * t302;
t295 = cos(t300);
t507 = Icges(4,5) * t295;
t146 = Icges(4,1) * t478 - t266 - t507;
t142 = Icges(6,5) * t295 + t294 * t370;
t144 = -Icges(5,4) * t295 + t294 * t371;
t607 = t142 + t144;
t602 = t146 + t607;
t510 = Icges(4,4) * t301;
t249 = Icges(4,1) * t302 - t510;
t147 = Icges(4,5) * t294 + t249 * t295;
t624 = t631 * t294 + t630 * t295;
t601 = t147 + t624;
t598 = -t633 * t302 + t297 + t298;
t617 = -Icges(4,2) * t302 - t510 + t598;
t476 = t295 * t302;
t629 = (Icges(6,4) + Icges(5,5)) * t476;
t299 = Icges(4,4) * t302;
t369 = -Icges(4,2) * t301 + t299;
t506 = Icges(5,5) * t302;
t235 = Icges(5,3) * t301 + t506;
t508 = Icges(6,4) * t302;
t239 = Icges(6,2) * t301 + t508;
t618 = t235 + t239;
t628 = -t369 + t618;
t627 = t632 * t294;
t248 = Icges(4,1) * t301 + t299;
t597 = t634 * t301 + t248 - t506 - t508;
t132 = -Icges(5,6) * t295 + t235 * t294;
t136 = Icges(6,6) * t295 + t239 * t294;
t608 = -t132 - t136;
t477 = t295 * t301;
t619 = t633 * t477 + t627 + t629;
t237 = Icges(4,5) * t302 - Icges(4,6) * t301;
t135 = Icges(4,3) * t294 + t237 * t295;
t241 = Icges(5,4) * t302 + Icges(5,6) * t301;
t139 = Icges(5,2) * t294 + t241 * t295;
t626 = t135 + t139;
t611 = (Icges(4,6) - t632) * t302 + (Icges(4,5) + t631) * t301;
t625 = t249 + t630;
t418 = qJD(3) * t301;
t417 = qJD(3) * t302;
t623 = -t301 * t617 - t302 * t597;
t496 = Icges(4,3) * t295;
t134 = Icges(4,5) * t478 - Icges(4,6) * t479 - t496;
t501 = Icges(4,6) * t295;
t140 = Icges(4,4) * t478 - Icges(4,2) * t479 - t501;
t487 = t140 * t301;
t358 = -t146 * t302 + t487;
t361 = t136 * t301 + t142 * t302;
t138 = -Icges(5,2) * t295 + t241 * t294;
t488 = t138 * t295;
t233 = Icges(6,5) * t302 + Icges(6,6) * t301;
t130 = Icges(6,3) * t295 + t233 * t294;
t490 = t130 * t295;
t365 = t132 * t301 + t144 * t302;
t561 = t294 * t365;
t582 = t294 * t361 - t488 + t490 + t561;
t622 = -t134 * t295 - t294 * t358 + t582;
t141 = Icges(4,6) * t294 + t295 * t369;
t113 = t147 * t478;
t389 = t135 * t295 - t113;
t51 = -t141 * t479 - t389;
t131 = -Icges(6,3) * t294 + t233 * t295;
t128 = t295 * t131;
t581 = -t139 * t295 + t624 * t478 + t479 * t619 + t128;
t621 = t51 + t581;
t620 = t601 * t301;
t616 = t626 * t294 + t476 * t601 + t619 * t477;
t126 = t294 * t138;
t615 = -t294 * t134 - t476 * t602 + t608 * t477 - t126;
t614 = t628 * qJD(3);
t613 = t625 * qJD(3);
t612 = t233 - t237 - t241;
t551 = t611 * t295;
t552 = t611 * t294;
t609 = -t417 * t617 + t418 * t597;
t491 = t130 * t294;
t572 = -t140 * t477 - t491 - t615;
t571 = -t131 * t294 - t141 * t477 + t616;
t606 = t623 * t294 + t551;
t605 = -t295 * t623 + t552;
t604 = t140 + t608;
t603 = t141 - t619;
t515 = -rSges(6,3) - qJ(5);
t596 = t613 * t302 + t614 * t301 + (-t597 * t301 + t617 * t302) * qJD(3) + t611 * qJD(2);
t595 = t611 * qJD(3);
t486 = t141 * t301;
t594 = t619 * t301 + t601 * t302 - t486;
t593 = t358 - t361 - t365;
t592 = -qJD(2) * t623 + t612 * qJD(3);
t591 = (-t595 * t294 + (-t131 + t593 + t626) * qJD(2)) * t295;
t422 = qJD(2) * t295;
t590 = t605 * qJD(2);
t589 = (t571 * t294 - t572 * t295) * qJD(3);
t588 = (t621 * t294 - t622 * t295) * qJD(3);
t587 = t597 - t628;
t586 = t617 + t625;
t585 = t617 * t295 + t601;
t584 = -t248 * t295 - t634 * t477 - t603 + t629;
t583 = t606 * qJD(2);
t296 = qJD(4) * t301;
t251 = t295 * t296;
t253 = pkin(3) * t301 - qJ(4) * t302;
t254 = rSges(6,1) * t301 - rSges(6,2) * t302;
t529 = pkin(4) * t301;
t392 = t254 + t529;
t382 = -t253 - t392;
t415 = qJD(5) * t294;
t419 = qJD(3) * t295;
t321 = t382 * t419 - t415;
t580 = t251 + t321;
t273 = pkin(4) * t476;
t579 = t515 * t294 + t273;
t578 = -t583 + t588;
t577 = t589 + t590;
t576 = t593 * qJD(3) + t609 * t294 + ((t618 * t295 - t141 + t627) * t302 - t620) * qJD(2);
t575 = t594 * qJD(3) - t609 * t295 + ((-t294 * t369 + t501 - t608) * t302 + (-t249 * t294 + t507 - t607) * t301) * qJD(2);
t574 = -t592 * t294 + t596 * t295;
t573 = t596 * t294 + t592 * t295;
t570 = t603 * t302 + t620;
t569 = t602 * t301 + t604 * t302;
t566 = t488 + t616;
t564 = -t595 * t295 + (-t237 * t294 + t130 - t138 + t496 - t594) * qJD(2);
t412 = -rSges(6,1) - pkin(3) - pkin(4);
t493 = qJ(4) * t301;
t523 = rSges(6,2) * t301;
t563 = qJD(2) * (t412 * t302 - pkin(2) - t493 - t523) - qJD(5);
t562 = 0.2e1 * qJD(3);
t187 = t254 * t294;
t524 = rSges(5,1) * t301;
t255 = -rSges(5,3) * t302 + t524;
t560 = (qJD(3) * t255 - t296) * t294;
t257 = pkin(3) * t302 + t493;
t190 = t257 * t294;
t163 = qJD(2) * t190;
t289 = t295 * pkin(6);
t208 = pkin(2) * t294 - t289;
t205 = qJD(2) * t208;
t559 = -t163 - t205;
t209 = t295 * pkin(2) + t294 * pkin(6);
t285 = t294 * rSges(5,2);
t159 = rSges(5,1) * t476 + rSges(5,3) * t477 + t285;
t274 = pkin(3) * t476;
t195 = qJ(4) * t477 + t274;
t383 = t195 + t209;
t556 = t383 + t159;
t258 = rSges(6,1) * t302 + t523;
t539 = -rSges(6,3) * t295 - t258 * t294;
t453 = pkin(4) * t478 + qJ(5) * t295 - t539;
t553 = qJD(2) * t453;
t550 = (-t301 * t587 + t302 * t586) * qJD(2);
t549 = -t301 * t585 + t302 * t584;
t548 = t612 * qJD(2);
t547 = -Icges(4,2) * t478 + t598 * t294 - t266 + t602;
t546 = t597 * t294 + t604;
t535 = t301 * t547 + t302 * t546;
t534 = m(5) / 0.2e1;
t533 = m(6) / 0.2e1;
t530 = -rSges(5,1) - pkin(3);
t528 = pkin(4) * t302;
t526 = qJD(2) / 0.2e1;
t525 = rSges(4,1) * t302;
t256 = rSges(4,1) * t301 + rSges(4,2) * t302;
t194 = t256 * t295;
t420 = qJD(3) * t294;
t406 = t256 * t420;
t284 = t294 * rSges(4,3);
t160 = rSges(4,1) * t476 - rSges(4,2) * t477 + t284;
t450 = t160 + t209;
t60 = qJD(2) * t450 - t406;
t520 = t194 * t60;
t432 = rSges(4,2) * t479 + t295 * rSges(4,3);
t157 = rSges(4,1) * t478 - t432;
t405 = t256 * t419;
t59 = -t405 + (-t157 - t208) * qJD(2);
t519 = t294 * t59;
t518 = t295 * t59;
t517 = -rSges(6,2) - qJ(4);
t516 = -rSges(5,3) - qJ(4);
t202 = t209 * qJD(2);
t421 = qJD(2) * t301;
t166 = t294 * t417 + t295 * t421;
t404 = t294 * t418;
t230 = pkin(3) * t404;
t76 = qJ(4) * t166 + qJD(2) * t274 + t294 * t296 - t230;
t514 = -t202 - t76;
t401 = t295 * t417;
t227 = rSges(6,2) * t401;
t402 = t295 * t418;
t423 = qJD(2) * t294;
t329 = -t302 * t423 - t402;
t472 = -rSges(6,1) * t402 + pkin(4) * t329 - qJ(5) * t422 + qJD(2) * t539 + t227 - t415;
t229 = pkin(4) * t404;
t414 = qJD(5) * t295;
t386 = -t229 + t414;
t471 = -qJD(3) * t187 + t386 + (t258 * t295 + t579) * qJD(2);
t431 = rSges(6,1) * t476 + rSges(6,2) * t477;
t452 = t431 + t579;
t451 = -t159 - t195;
t449 = t294 * t190 + t295 * t195;
t191 = t253 * t295;
t416 = qJD(4) * t302;
t448 = -qJD(2) * t191 + t294 * t416;
t447 = -t190 - t208;
t203 = t253 * t420;
t413 = qJD(3) * qJD(4);
t399 = t302 * t413;
t446 = qJD(2) * t203 + t295 * t399;
t198 = qJD(3) * t257 - t416;
t445 = -t258 * qJD(3) - t198;
t259 = rSges(5,1) * t302 + rSges(5,3) * t301;
t444 = -t259 * qJD(3) - t198;
t443 = qJ(4) * t401 + t251;
t442 = rSges(5,2) * t422 + rSges(5,3) * t401;
t407 = t294 * t421;
t441 = rSges(4,2) * t407 + rSges(4,3) * t422;
t434 = -t253 - t255;
t433 = -t257 - t259;
t430 = t294 ^ 2 + t295 ^ 2;
t75 = pkin(3) * t329 - qJ(4) * t407 + t443;
t411 = t190 * t422 + t294 * t76 + t295 * t75;
t186 = t253 * t294;
t410 = -t186 * t420 - t191 * t419 + t296;
t409 = -t195 - t452;
t277 = pkin(6) * t422;
t408 = t277 + t443;
t400 = -pkin(2) - t525;
t396 = -t420 / 0.2e1;
t393 = t419 / 0.2e1;
t390 = t434 * t295;
t388 = -t134 + t486;
t387 = t301 * t413 + t76 * t420 + (t163 + t75) * t419;
t196 = qJD(2) * (-pkin(2) * t423 + t277);
t384 = t294 * t399 + t196 + (t251 + t75) * qJD(2);
t377 = t302 * t530 - pkin(2);
t374 = -rSges(4,2) * t301 + t525;
t373 = -t294 * t60 - t518;
t356 = t157 * t294 + t160 * t295;
t351 = -pkin(4) * t417 + t445;
t348 = t190 * t420 + t195 * t419 + qJD(1) - t416;
t347 = qJD(3) * t390 + t251;
t189 = t256 * t294;
t188 = t255 * t294;
t328 = (-t528 * qJD(3) + t445) * qJD(3);
t103 = rSges(4,1) * t329 - rSges(4,2) * t401 + t441;
t106 = -qJD(3) * t189 + (t295 * t374 + t284) * qJD(2);
t320 = t103 * t295 + t106 * t294 + (t157 * t295 - t160 * t294) * qJD(2);
t287 = t295 * rSges(5,2);
t252 = t295 * t416;
t221 = t374 * qJD(3);
t204 = t253 * t423;
t193 = t255 * t295;
t192 = t254 * t295;
t167 = t401 - t407;
t165 = t430 * t418;
t156 = t259 * t294 - t287;
t105 = -qJD(3) * t188 + (t259 * t295 + t285) * qJD(2);
t102 = rSges(5,1) * t329 - rSges(5,3) * t407 + t442;
t58 = qJD(3) * t356 + qJD(1);
t45 = qJD(2) * t556 - t203 - t560;
t44 = (-t156 + t447) * qJD(2) + t347;
t43 = -t221 * t419 + (-t106 - t202 + t406) * qJD(2);
t42 = -t221 * t420 + t196 + (t103 - t405) * qJD(2);
t41 = (t156 * t294 + t159 * t295) * qJD(3) + t348;
t40 = -t203 + (-qJD(3) * t254 + t296) * t294 + (t209 - t409) * qJD(2) + t386;
t39 = (t447 - t453) * qJD(2) + t580;
t32 = (t294 * t453 + t295 * t452) * qJD(3) + t348;
t25 = t320 * qJD(3);
t24 = t444 * t419 + (-t105 + t514 + t560) * qJD(2) + t446;
t23 = qJD(2) * t102 + (qJD(2) * t390 + t294 * t444) * qJD(3) + t384;
t16 = t328 * t295 + (-t414 + (qJD(3) * t392 - t296) * t294 - t471 + t514) * qJD(2) + t446;
t15 = t328 * t294 + (t321 + t472) * qJD(2) + t384;
t2 = (t102 * t295 + t105 * t294 + (t156 * t295 + t294 * t451) * qJD(2)) * qJD(3) + t387;
t1 = (t472 * t295 + t471 * t294 + (t294 * t409 + t295 * t453) * qJD(2)) * qJD(3) + t387;
t3 = [m(4) * t25 + m(5) * t2 + m(6) * t1; (-qJD(3) * t623 + t613 * t301 - t614 * t302) * qJD(2) + (t16 * t289 + t39 * (t229 + t230) + t15 * (t273 + t383 + t431) + (t16 * t515 + t563 * t39) * t295 + (-t16 * pkin(2) + t15 * t515 + (t16 * t517 + t39 * (rSges(6,1) * qJD(3) - qJD(4))) * t301 + (qJD(3) * t39 * t517 + t16 * t412) * t302 + t39 * (-pkin(6) - t515) * qJD(2)) * t294 + (t39 + t553 - t559 + t227 + t408 + (qJD(2) * t515 + t412 * t418) * t295 + t563 * t294 - t580) * t40) * m(6) + (-(-qJD(2) * t156 + t347 - t44 + t559) * t45 + t24 * (t287 + t289) + t44 * t230 + t23 * t556 + t45 * (t408 + t442) + (t45 * t530 * t418 + t44 * (t301 * t516 + t377) * qJD(2)) * t295 + (t24 * t377 + (-t44 * qJD(4) + t24 * t516) * t301 + t44 * (t302 * t516 + t524) * qJD(3) + (t44 * (-rSges(5,2) - pkin(6)) + t45 * (-pkin(2) + t433)) * qJD(2)) * t294) * m(5) + (t43 * (t294 * t400 + t289 + t432) + t42 * t450 + t60 * (t277 + t441) + (t256 * t519 - t520) * qJD(3) + ((-pkin(2) - t374) * t518 + (t59 * (-rSges(4,3) - pkin(6)) + t60 * t400) * t294) * qJD(2) - (-qJD(2) * t157 - t205 - t405 - t59) * t60) * m(4) + (((t51 - t113 + (t135 + t487) * t295 + t615) * t295 + (-t561 + (-t131 - t361) * t294 + t566 + t582) * t294) * qJD(3) + t590) * t393 + (((t295 * t388 + t490 - t566 + t571) * t295 + (t294 * t388 - t126 + t128 + t389 + t491 + t572 - t581) * t294) * qJD(3) + t578 + t583) * t396 + (t574 + t575) * t420 / 0.2e1 - (t573 - t576 + t577) * t419 / 0.2e1 + ((t569 - t606) * t294 + (t570 + t605) * t295) * qJD(3) * t526; (t1 * t449 + (t1 * t453 + t15 * t382) * t294 + (t1 * t452 + t16 * t382) * t295 + (t351 * t294 - t448 + (pkin(4) * t477 + t295 * t382 + t192) * qJD(2)) * t40 + (-t186 * qJD(2) + t351 * t295 + t204 - t252) * t39 - (t294 * t40 + t295 * t39) * qJD(3) * (-t257 - t258 - t528) + (-t410 - (-t187 * t294 - t192 * t295 - t430 * t529) * qJD(3) + t411 + (qJD(2) * t409 + t471) * t294 + (t472 + t553) * t295) * t32) * m(6) + (-t44 * (t252 + (t186 + t188) * qJD(2)) - t45 * (-qJD(2) * t193 + t448) - t41 * t410 - ((-t41 * t193 + t433 * t44) * t295 + (-t41 * t188 + t433 * t45) * t294) * qJD(3) + t44 * t204 + t2 * t449 + t41 * t411 + (t24 * t434 + t44 * t444 + t2 * t159 + t41 * t102 + (t41 * t156 + t434 * t45) * qJD(2)) * t295 + (t23 * t434 + t45 * t444 + t2 * t156 + t41 * t105 + (t44 * t255 + t41 * t451) * qJD(2)) * t294) * m(5) + (t25 * t356 + t58 * t320 + t373 * t221 + (-t42 * t294 - t43 * t295 + (-t295 * t60 + t519) * qJD(2)) * t256 - (t189 * t59 - t520) * qJD(2) - (t58 * (-t189 * t294 - t194 * t295) + t373 * t374) * qJD(3)) * m(4) - (((t294 * t585 - t547 * t295) * t302 + (t294 * t584 + t546 * t295) * t301) * qJD(3) + (t586 * t301 + t587 * t302) * qJD(2)) * qJD(2) / 0.2e1 + (t576 * t295 + t575 * t294 + (t294 * t569 + t295 * t570) * qJD(2)) * t526 + ((-t420 * t551 - t548) * t294 + ((t535 * t295 + (t549 + t552) * t294) * qJD(3) + t550) * t295) * t396 + ((-t419 * t552 + t548) * t295 + ((t549 * t294 + (t535 + t551) * t295) * qJD(3) + t550) * t294) * t393 + (t574 * qJD(2) + (t571 * t422 + (t572 * qJD(2) + t564 * t294 - t591) * t294) * t562) * t294 / 0.2e1 - (t573 * qJD(2) + ((t621 * qJD(2) + t591) * t295 + (t622 * qJD(2) - t564 * t295) * t294) * t562) * t295 / 0.2e1 + (t578 + t588) * t423 / 0.2e1 + (t577 + t589) * t422 / 0.2e1; -m(5) * (t165 * t41 + t166 * t45 + t167 * t44) - m(6) * (t165 * t32 + t166 * t40 + t167 * t39) + 0.2e1 * ((t419 * t44 + t420 * t45 - t2) * t534 + (t39 * t419 + t40 * t420 - t1) * t533) * t302 + 0.2e1 * ((qJD(3) * t41 + t23 * t294 + t24 * t295 + t422 * t45 - t423 * t44) * t534 + (qJD(3) * t32 + t15 * t294 + t16 * t295 - t39 * t423 + t40 * t422) * t533) * t301; m(6) * (t15 * t295 - t16 * t294);];
tauc = t3(:);
