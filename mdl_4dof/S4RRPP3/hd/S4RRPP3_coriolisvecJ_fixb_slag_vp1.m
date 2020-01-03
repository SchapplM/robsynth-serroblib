% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:23
% EndTime: 2019-12-31 16:57:47
% DurationCPUTime: 20.07s
% Computational Cost: add. (6934->523), mult. (10976->663), div. (0->0), fcn. (8528->6), ass. (0->296)
t635 = Icges(3,3) + Icges(4,3);
t298 = qJ(2) + pkin(6);
t275 = sin(t298);
t276 = cos(t298);
t302 = sin(qJ(2));
t304 = cos(qJ(2));
t634 = Icges(3,5) * t304 + Icges(4,5) * t276 - Icges(3,6) * t302 - Icges(4,6) * t275;
t206 = Icges(5,4) * t276 + Icges(5,6) * t275;
t303 = sin(qJ(1));
t305 = cos(qJ(1));
t143 = Icges(5,2) * t303 + t206 * t305;
t616 = t635 * t303 + t634 * t305;
t633 = t143 + t616;
t508 = Icges(4,4) * t275;
t207 = Icges(4,2) * t276 + t508;
t265 = Icges(5,5) * t275;
t632 = Icges(5,3) * t276 + t207 - t265;
t504 = Icges(5,5) * t276;
t209 = Icges(5,1) * t275 - t504;
t266 = Icges(4,4) * t276;
t631 = Icges(4,1) * t275 + t209 + t266;
t630 = t635 * t305;
t474 = t303 * t304;
t476 = t302 * t303;
t478 = t276 * t303;
t480 = t275 * t303;
t601 = -Icges(3,5) * t474 - Icges(4,5) * t478 + Icges(3,6) * t476 + Icges(4,6) * t480 + t630;
t375 = Icges(5,1) * t276 + t265;
t146 = -Icges(5,4) * t305 + t303 * t375;
t246 = Icges(4,4) * t480;
t505 = Icges(4,5) * t305;
t148 = Icges(4,1) * t478 - t246 - t505;
t623 = t146 + t148;
t147 = Icges(5,4) * t303 + t305 * t375;
t212 = Icges(4,1) * t276 - t508;
t149 = Icges(4,5) * t303 + t212 * t305;
t622 = t147 + t149;
t202 = Icges(5,3) * t275 + t504;
t373 = -Icges(4,2) * t275 + t266;
t629 = t202 - t373;
t628 = t212 + t375;
t499 = Icges(4,6) * t305;
t144 = Icges(4,4) * t478 - Icges(4,2) * t480 - t499;
t500 = Icges(3,6) * t305;
t159 = Icges(3,4) * t474 - Icges(3,2) * t476 - t500;
t627 = t144 * t275 + t159 * t302;
t617 = Icges(3,5) * t302 + Icges(3,6) * t304 + (Icges(4,6) - Icges(5,6)) * t276 + (Icges(5,4) + Icges(4,5)) * t275;
t509 = Icges(3,4) * t302;
t235 = Icges(3,1) * t304 - t509;
t162 = Icges(3,5) * t303 + t235 * t305;
t626 = -t149 * t478 - t162 * t474;
t138 = -Icges(5,6) * t305 + t202 * t303;
t625 = -t138 + t144;
t477 = t276 * t305;
t245 = Icges(5,5) * t477;
t479 = t275 * t305;
t498 = Icges(5,6) * t303;
t139 = Icges(5,3) * t479 + t245 + t498;
t145 = Icges(4,6) * t303 + t305 * t373;
t624 = t139 - t145;
t621 = t206 + t634;
t620 = t632 * qJD(2);
t619 = t631 * qJD(2);
t262 = Icges(3,4) * t476;
t506 = Icges(3,5) * t305;
t161 = Icges(3,1) * t474 - t262 - t506;
t618 = t148 * t276 + t161 * t304 - t627;
t232 = Icges(3,2) * t304 + t509;
t288 = Icges(3,4) * t304;
t234 = Icges(3,1) * t302 + t288;
t606 = t302 * t232 - t304 * t234 + t632 * t275 - t631 * t276;
t615 = t629 * qJD(2);
t614 = t628 * qJD(2);
t382 = -t139 * t480 + t143 * t305 - t147 * t478;
t374 = -Icges(3,2) * t302 + t288;
t160 = Icges(3,6) * t303 + t305 * t374;
t610 = -t616 * t305 - t626;
t585 = -t145 * t480 - t160 * t476 + t610;
t611 = -t382 + t585;
t609 = t145 * t275 + t160 * t302;
t473 = t304 * t305;
t608 = t139 * t479 + t162 * t473 + t633 * t303 + t622 * t477;
t142 = -Icges(5,2) * t305 + t206 * t303;
t130 = t303 * t142;
t607 = -t138 * t479 - t161 * t473 + t601 * t303 - t623 * t477 - t130;
t555 = t617 * t305;
t554 = t617 * t303;
t605 = t620 * t305 + (t303 * t373 - t138 - t499) * qJD(1);
t604 = t620 * t303 + (t202 * t305 - t145 + t498) * qJD(1);
t603 = -t619 * t305 + (-t212 * t303 - t146 + t505) * qJD(1);
t602 = -t622 * qJD(1) + t619 * t303;
t492 = t142 * t305;
t370 = t138 * t275 + t146 * t276;
t551 = t303 * t370;
t46 = -t492 + t551;
t599 = t618 * t303 + t601 * t305 + t46;
t475 = t302 * t305;
t583 = -t144 * t479 - t159 * t475 - t607;
t582 = -t145 * t479 - t160 * t475 + t608;
t598 = t606 * t303 + t555;
t597 = -t606 * t305 + t554;
t596 = -t160 * t304 - t162 * t302 - t622 * t275 + t624 * t276;
t560 = t159 * t304 + t161 * t302 + t623 * t275 + t625 * t276;
t595 = t617 * qJD(2);
t594 = t139 * t275 + t162 * t304 + t622 * t276 - t609;
t593 = -t370 - t618;
t592 = t625 * t305 + (-Icges(5,1) * t479 + t209 * t305 + t245 + t624) * t303;
t591 = t631 - t629;
t590 = -t632 + t628;
t589 = (Icges(4,2) * t478 + t246 - t623) * t305 + (-t207 * t305 + t622) * t303;
t221 = t374 * qJD(2);
t222 = t235 * qJD(2);
t588 = -t221 * t302 + t222 * t304 + t614 * t276 + t615 * t275 + (-t232 * t304 - t234 * t302 - t275 * t631 - t276 * t632) * qJD(2) + t617 * qJD(1);
t587 = t633 * qJD(1);
t586 = t606 * qJD(1) + t621 * qJD(2);
t584 = rSges(3,2) * t302;
t572 = rSges(5,3) + qJ(4);
t581 = t597 * qJD(1);
t344 = qJD(2) * t232;
t105 = -t305 * t344 + (-t303 * t374 + t500) * qJD(1);
t347 = qJD(2) * t234;
t107 = -t305 * t347 + (-t235 * t303 + t506) * qJD(1);
t580 = t596 * qJD(2) - t105 * t302 + t107 * t304 + t605 * t275 + t603 * t276 + t587;
t106 = qJD(1) * t160 - t303 * t344;
t108 = qJD(1) * t162 - t303 * t347;
t545 = qJD(1) * t142;
t579 = t601 * qJD(1) + t560 * qJD(2) + t106 * t302 - t108 * t304 - t604 * t275 + t602 * t276 - t545;
t578 = (t303 * t582 - t305 * t583) * qJD(2);
t577 = (t303 * t611 - t599 * t305) * qJD(2);
t576 = t598 * qJD(1);
t575 = qJD(1) * t593 - t303 * t595 + t587;
t574 = -t545 - t595 * t305 + (-t303 * t634 - t594 + t630) * qJD(1);
t573 = 0.2e1 * qJD(2);
t528 = rSges(5,1) + pkin(3);
t152 = rSges(4,1) * t478 - rSges(4,2) * t480 - t305 * rSges(4,3);
t289 = t303 * rSges(4,3);
t154 = rSges(4,1) * t477 - rSges(4,2) * t479 + t289;
t365 = t152 * t303 + t154 * t305;
t296 = t305 * pkin(5);
t251 = pkin(1) * t303 - t296;
t301 = -qJ(3) - pkin(5);
t269 = t305 * t301;
t526 = pkin(2) * t304;
t270 = pkin(1) + t526;
t436 = -t303 * t270 - t269;
t136 = t251 + t436;
t295 = t303 * pkin(5);
t252 = t305 * pkin(1) + t295;
t255 = t305 * t270;
t388 = -t301 * t303 + t255;
t137 = t388 - t252;
t464 = -t303 * t136 + t305 * t137;
t337 = t365 + t464;
t423 = qJD(2) * t305;
t424 = qJD(2) * t303;
t465 = -t136 * t424 + t137 * t423;
t41 = qJD(2) * t365 + t465;
t571 = qJD(2) * t337 + t41;
t422 = qJD(4) * t275;
t213 = pkin(3) * t275 - qJ(4) * t276;
t214 = rSges(5,1) * t275 - rSges(5,3) * t276;
t443 = t213 + t214;
t570 = (qJD(2) * t443 - t422) * t303;
t217 = rSges(5,1) * t276 + rSges(5,3) * t275;
t569 = pkin(3) * t276 + qJ(4) * t275 + t217;
t568 = t601 + t609;
t567 = -t576 + t577;
t566 = t578 + t581;
t565 = t303 * t586 + t305 * t588;
t564 = t303 * t588 - t305 * t586;
t563 = t593 * qJD(2) - t106 * t304 - t108 * t302 + t602 * t275 + t604 * t276;
t562 = t594 * qJD(2) + t105 * t304 + t107 * t302 + t603 * t275 - t605 * t276;
t332 = t159 * t305 - t160 * t303;
t533 = t303 * (-t232 * t305 + t162) - t305 * (-Icges(3,2) * t474 + t161 - t262);
t559 = -t589 * t275 + t276 * t592 - t302 * t533 + t332 * t304;
t439 = t234 + t374;
t440 = -t232 + t235;
t558 = (-t275 * t591 + t276 * t590 - t302 * t439 + t304 * t440) * qJD(1);
t557 = t492 + t608;
t556 = t621 * qJD(1);
t530 = t303 / 0.2e1;
t552 = t275 * t528;
t132 = qJD(1) * t136;
t226 = qJD(1) * t251;
t550 = t132 - t226;
t406 = t276 * t423;
t425 = qJD(1) * t305;
t548 = rSges(5,2) * t425 + t572 * t406;
t547 = t303 ^ 2 + t305 ^ 2;
t527 = pkin(2) * t302;
t386 = -t443 - t527;
t357 = t386 * qJD(2);
t237 = t305 * t422;
t277 = qJD(3) * t303;
t437 = t237 + t277;
t546 = t305 * t357 + t437;
t529 = -t305 / 0.2e1;
t524 = qJD(1) / 0.2e1;
t523 = pkin(1) - t270;
t426 = qJD(1) * t303;
t336 = -t275 * t423 - t276 * t426;
t411 = t275 * t426;
t522 = t336 * t528 - t411 * t572 + t237 + t548;
t176 = t214 * t303;
t291 = t303 * rSges(5,2);
t521 = (pkin(3) * t425 + qJ(4) * t424) * t276 + (qJ(4) * t425 + (-pkin(3) * qJD(2) + qJD(4)) * t303) * t275 - qJD(2) * t176 + (t217 * t305 + t291) * qJD(1);
t520 = rSges(3,1) * t304;
t519 = rSges(4,1) * t276;
t238 = rSges(3,1) * t302 + rSges(3,2) * t304;
t200 = t238 * t305;
t408 = t238 * t424;
t290 = t303 * rSges(3,3);
t184 = rSges(3,1) * t473 - rSges(3,2) * t475 + t290;
t448 = t184 + t252;
t93 = qJD(1) * t448 - t408;
t518 = t200 * t93;
t433 = rSges(3,2) * t476 + t305 * rSges(3,3);
t183 = rSges(3,1) * t474 - t433;
t407 = t238 * t423;
t92 = -t407 + (-t183 - t251) * qJD(1);
t517 = t303 * t92;
t516 = t305 * t92;
t215 = rSges(4,1) * t275 + rSges(4,2) * t276;
t395 = -t215 - t527;
t381 = t305 * t395;
t356 = qJD(2) * t381;
t335 = t277 + t356;
t461 = t136 - t251;
t44 = (-t152 + t461) * qJD(1) + t335;
t514 = t44 * t215;
t257 = t424 * t527;
t434 = qJD(3) * t305 + t257;
t412 = t301 * t426 + t434;
t102 = (-t305 * t523 - t295) * qJD(1) - t412;
t225 = t252 * qJD(1);
t469 = -t102 - t225;
t460 = -t137 - t154;
t421 = qJD(4) * t276;
t455 = -qJD(2) * t569 + t421;
t294 = t305 * rSges(5,2);
t454 = t303 * t569 - t294;
t453 = t477 * t528 + t479 * t572 + t291;
t450 = -t213 * t303 - t176;
t449 = t443 * t305;
t441 = rSges(4,2) * t411 + rSges(4,3) * t425;
t420 = qJD(1) * qJD(3);
t438 = qJD(1) * t257 + t305 * t420;
t435 = rSges(3,3) * t425 + t426 * t584;
t419 = pkin(2) * t475;
t418 = qJD(2) ^ 2 * t526;
t274 = pkin(5) * t425;
t405 = t302 * t423;
t101 = -pkin(2) * t405 - t274 + t277 + (t303 * t523 - t269) * qJD(1);
t417 = t305 * t101 + t303 * t102 - t136 * t425;
t416 = t102 * t424 + (t101 - t132) * t423;
t415 = qJD(2) * t526;
t219 = qJD(1) * (-pkin(1) * t426 + t274);
t414 = qJD(1) * t101 + t303 * t420 + t219;
t413 = -t137 - t453;
t409 = t215 * t424;
t404 = t304 * t423;
t402 = -pkin(1) - t520;
t399 = -t424 / 0.2e1;
t396 = t423 / 0.2e1;
t218 = -rSges(4,2) * t275 + t519;
t394 = -t218 - t526;
t387 = t547 * t527;
t385 = -t526 - t569;
t379 = t520 - t584;
t378 = -t303 * t93 - t516;
t359 = -t415 + t455;
t192 = t218 * qJD(2);
t355 = -qJD(2) * t192 - t418;
t199 = t238 * t303;
t177 = t215 * t303;
t88 = (t183 * t303 + t184 * t305) * qJD(2);
t328 = -t418 + (t421 + t455) * qJD(2);
t327 = -t270 - t569;
t223 = t379 * qJD(2);
t181 = t215 * t305;
t110 = -qJD(2) * t199 + (t305 * t379 + t290) * qJD(1);
t109 = -rSges(3,2) * t404 + (-t304 * t426 - t405) * rSges(3,1) + t435;
t87 = -qJD(2) * t177 + (t218 * t305 + t289) * qJD(1);
t85 = rSges(4,1) * t336 - rSges(4,2) * t406 + t441;
t55 = -t223 * t423 + (-t110 - t225 + t408) * qJD(1);
t54 = -t223 * t424 + t219 + (t109 - t407) * qJD(1);
t45 = -t409 + (t252 - t460) * qJD(1) - t434;
t36 = -t570 + (t252 - t413) * qJD(1) - t434;
t35 = (-t454 + t461) * qJD(1) + t546;
t32 = -t421 + (t303 * t454 + t305 * t453) * qJD(2) + t465;
t27 = t355 * t305 + (-t87 + t409 + t469) * qJD(1) + t438;
t26 = t355 * t303 + (t85 + t356) * qJD(1) + t414;
t11 = t328 * t305 + (t469 - t521 + t570) * qJD(1) + t438;
t10 = t328 * t303 + ((t357 + t422) * t305 + t522) * qJD(1) + t414;
t1 = (t422 + t522 * t305 + t521 * t303 + (t303 * t413 + t305 * t454) * qJD(1)) * qJD(2) + t416;
t2 = [(-t606 * qJD(2) + t221 * t304 + t222 * t302 + t614 * t275 - t615 * t276) * qJD(1) + (t11 * (t294 + t436) + t35 * t412 + t10 * (t255 + t453) + t36 * (t437 + t548) + (t36 * (-t527 - t552) * qJD(2) + (-t36 * t301 + t327 * t35) * qJD(1)) * t305 + (-t10 * t301 - t11 * t528 * t276 + (-t35 * qJD(4) - t11 * t572) * t275 + t35 * (-t276 * t572 + t552) * qJD(2) + (-t35 * rSges(5,2) + t327 * t36) * qJD(1)) * t303 - (-qJD(1) * t454 - t35 + t546 + t550) * t36) * m(5) + (-(-qJD(1) * t152 + t335 - t44 + t550) * t45 + t27 * (-t152 + t436) + t44 * t412 + t26 * (t154 + t388) + t45 * (t277 + t441) + (t514 * t303 + t381 * t45) * qJD(2) + ((-t44 * rSges(4,3) + t45 * (-t270 - t519)) * t303 + (t44 * (-t218 - t270) - t45 * t301) * t305) * qJD(1)) * m(4) + (t55 * (t303 * t402 + t296 + t433) + t54 * t448 + t93 * (t274 + t435) + (t238 * t517 - t518) * qJD(2) + ((-pkin(1) - t379) * t516 + (t92 * (-rSges(3,3) - pkin(5)) + t93 * t402) * t303) * qJD(1) - (-qJD(1) * t183 - t226 - t407 - t92) * t93) * m(3) + (((t46 - t551 + t557) * t303 + ((t616 + t627) * t305 + t585 + t607 + t626) * t305) * qJD(2) + t581) * t396 + (t562 + t565) * t424 / 0.2e1 + (((t305 * t568 - t557 + t582) * t305 + (t303 * t568 - t130 + t382 + t583 - t610) * t303) * qJD(2) + t567 + t576) * t399 - (-t563 + t564 + t566) * t423 / 0.2e1 + ((t560 - t598) * t303 + (-t596 + t597) * t305) * qJD(2) * t524; (-(t44 * t177 + t45 * (-t181 - t419)) * qJD(1) - (-t41 * t387 + (-t41 * t181 + t394 * t44) * t305 + (-t41 * t177 + t394 * t45) * t303) * qJD(2) + t27 * t381 - t44 * pkin(2) * t404 + (t423 * t85 + t424 * t87 + t416) * t337 + t41 * t417 + (-t44 * t192 + t41 * t85 + (t152 * t571 + t395 * t45) * qJD(1)) * t305 + (t26 * t395 + t45 * (-t192 - t415) + t41 * t87 + (t460 * t571 + t514) * qJD(1)) * t303) * m(4) + (-(t199 * t92 - t518) * qJD(1) - (t88 * (-t199 * t303 - t200 * t305) + t378 * t379) * qJD(2) + 0.2e1 * t88 * (t109 * t305 + t110 * t303 + (t183 * t305 - t184 * t303) * qJD(1)) + t378 * t223 + (-t54 * t303 - t55 * t305 + (-t305 * t93 + t517) * qJD(1)) * t238) * m(3) - ((t275 * t592 + t589 * t276 + t332 * t302 + t304 * t533) * qJD(2) + (t590 * t275 + t591 * t276 + t302 * t440 + t304 * t439) * qJD(1)) * qJD(1) / 0.2e1 + (t563 * t305 + t562 * t303 + (t303 * t560 - t305 * t596) * qJD(1)) * t524 + ((-t424 * t555 + t556) * t303 + ((t303 * t554 + t559) * qJD(2) + t558) * t305) * t399 + ((-t423 * t554 - t556) * t305 + ((t305 * t555 + t559) * qJD(2) + t558) * t303) * t396 + (-(t275 * t32 + (t303 * t36 + t305 * t35) * t276) * qJD(4) - (-t35 * t450 + t36 * (-t419 - t449)) * qJD(1) - (-t32 * t387 + (-t32 * t449 + t35 * t385) * t305 + (t32 * t450 + t36 * t385) * t303) * qJD(2) + t1 * t464 + t32 * t417 + (t11 * t386 + t35 * t359 + t1 * t453 + t32 * t522 + (t32 * t454 + t36 * t386) * qJD(1)) * t305 + (t10 * t386 + t36 * t359 + t1 * t454 + t32 * t521 + (t32 * t413 + t35 * t443) * qJD(1)) * t303) * m(5) + (t565 * qJD(1) + ((t582 * qJD(1) + t579 * t305) * t305 + (t574 * t303 + t583 * qJD(1) + (-t575 + t580) * t305) * t303) * t573) * t530 + (t564 * qJD(1) + ((t611 * qJD(1) + t575 * t305) * t305 + (t580 * t303 + t599 * qJD(1) + (-t574 + t579) * t305) * t303) * t573) * t529 + (t567 + t577) * t426 / 0.2e1 + (t566 + t578) * t425 / 0.2e1; 0.2e1 * (t10 * t529 + t11 * t530) * m(5) + 0.2e1 * (t26 * t529 + t27 * t530) * m(4); (-t1 * t276 + 0.2e1 * (t10 * t530 + t11 * t305 / 0.2e1 + (0.1e1 / 0.2e1 - t547 / 0.2e1) * qJD(2) * t32) * t275) * m(5);];
tauc = t2(:);
