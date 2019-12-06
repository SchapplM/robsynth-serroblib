% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP1_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:11
% EndTime: 2019-12-05 16:06:45
% DurationCPUTime: 22.26s
% Computational Cost: add. (12032->568), mult. (11291->712), div. (0->0), fcn. (8717->6), ass. (0->330)
t628 = Icges(4,3) + Icges(5,3);
t305 = qJ(3) + pkin(8);
t300 = sin(t305);
t302 = cos(t305);
t210 = Icges(5,5) * t302 - Icges(5,6) * t300;
t307 = sin(qJ(3));
t308 = cos(qJ(3));
t254 = Icges(4,5) * t308 - Icges(4,6) * t307;
t627 = t210 + t254;
t212 = Icges(6,4) * t302 + Icges(6,6) * t300;
t304 = pkin(7) + qJ(2);
t299 = sin(t304);
t301 = cos(t304);
t147 = Icges(6,2) * t299 + t212 * t301;
t604 = t628 * t299 + t627 * t301;
t626 = t147 + t604;
t515 = Icges(5,4) * t300;
t213 = Icges(5,2) * t302 + t515;
t282 = Icges(6,5) * t300;
t371 = Icges(6,3) * t302 - t282;
t623 = t213 + t371;
t511 = Icges(6,5) * t302;
t215 = Icges(6,1) * t300 - t511;
t286 = Icges(5,4) * t302;
t217 = Icges(5,1) * t300 + t286;
t622 = t215 + t217;
t625 = t628 * t301;
t484 = t299 * t308;
t485 = t299 * t307;
t486 = t299 * t302;
t487 = t299 * t300;
t595 = -Icges(4,5) * t484 - Icges(5,5) * t486 + Icges(4,6) * t485 + Icges(5,6) * t487 + t625;
t374 = Icges(6,1) * t302 + t282;
t150 = -Icges(6,4) * t301 + t299 * t374;
t242 = Icges(5,4) * t487;
t512 = Icges(5,5) * t301;
t152 = Icges(5,1) * t486 - t242 - t512;
t624 = t150 + t152;
t151 = Icges(6,4) * t299 + t301 * t374;
t218 = Icges(5,1) * t302 - t515;
t153 = Icges(5,5) * t299 + t218 * t301;
t618 = t151 + t153;
t506 = Icges(5,6) * t301;
t148 = Icges(5,4) * t486 - Icges(5,2) * t487 - t506;
t507 = Icges(4,6) * t301;
t160 = Icges(4,4) * t484 - Icges(4,2) * t485 - t507;
t621 = t148 * t300 + t160 * t307;
t516 = Icges(4,4) * t307;
t258 = Icges(4,1) * t308 - t516;
t163 = Icges(4,5) * t299 + t258 * t301;
t620 = -t153 * t486 - t163 * t484;
t264 = Icges(4,4) * t485;
t513 = Icges(4,5) * t301;
t162 = Icges(4,1) * t484 - t264 - t513;
t619 = t152 * t302 + t162 * t308 - t621;
t255 = Icges(4,2) * t308 + t516;
t303 = Icges(4,4) * t308;
t257 = Icges(4,1) * t307 + t303;
t607 = t255 * t307 - t257 * t308 + t623 * t300 - t622 * t302;
t208 = Icges(6,3) * t300 + t511;
t372 = -Icges(5,2) * t300 + t286;
t617 = (t208 - t372) * qJD(3);
t616 = (t218 + t374) * qJD(3);
t614 = t623 * qJD(3);
t613 = t622 * qJD(3);
t482 = t301 * t302;
t241 = Icges(6,5) * t482;
t483 = t300 * t301;
t505 = Icges(6,6) * t299;
t143 = Icges(6,3) * t483 + t241 + t505;
t381 = -t143 * t487 + t147 * t301 - t151 * t486;
t149 = Icges(5,6) * t299 + t301 * t372;
t373 = -Icges(4,2) * t307 + t303;
t161 = Icges(4,6) * t299 + t301 * t373;
t611 = -t604 * t301 - t620;
t583 = -t149 * t487 - t161 * t485 + t611;
t612 = -t381 + t583;
t480 = t301 * t308;
t610 = t143 * t483 + t163 * t480 + t626 * t299 + t618 * t482;
t146 = -Icges(6,2) * t301 + t212 * t299;
t132 = t299 * t146;
t142 = -Icges(6,6) * t301 + t208 * t299;
t609 = -t142 * t483 - t162 * t480 + t595 * t299 - t624 * t482 - t132;
t209 = Icges(5,5) * t300 + Icges(5,6) * t302;
t211 = Icges(6,4) * t300 - Icges(6,6) * t302;
t253 = Icges(4,5) * t307 + Icges(4,6) * t308;
t608 = t253 + t211 + t209;
t489 = t253 * t299;
t491 = t211 * t299;
t493 = t209 * t299;
t606 = t489 + t491 + t493;
t488 = t253 * t301;
t490 = t211 * t301;
t492 = t209 * t301;
t605 = -t492 - t490 - t488;
t603 = t149 * t300 + t161 * t307;
t499 = t146 * t301;
t369 = t142 * t300 + t150 * t302;
t560 = t299 * t369;
t48 = -t499 + t560;
t602 = t619 * t299 + t595 * t301 + t48;
t481 = t301 * t307;
t571 = -t148 * t483 - t160 * t481 - t609;
t570 = -t149 * t483 - t161 * t481 + t610;
t601 = t607 * t299 - t605;
t600 = -t607 * t301 + t606;
t599 = t614 * t301 + (t299 * t372 - t142 - t506) * qJD(2);
t598 = t614 * t299 + (t208 * t301 - t149 + t505) * qJD(2);
t597 = -t613 * t301 + (-t218 * t299 - t150 + t512) * qJD(2);
t596 = -t618 * qJD(2) + t613 * t299;
t594 = -t161 * t308 - t163 * t307 + (t143 - t149) * t302 - t618 * t300;
t568 = t160 * t308 + t162 * t307 + (-t142 + t148) * t302 + t624 * t300;
t233 = t373 * qJD(3);
t234 = t258 * qJD(3);
t593 = -t233 * t307 + t234 * t308 + t616 * t302 + t617 * t300 + (-t255 * t308 - t257 * t307 - t622 * t300 - t623 * t302) * qJD(3) + t608 * qJD(2);
t592 = t608 * qJD(3);
t591 = t143 * t300 + t163 * t308 + t618 * t302 - t603;
t590 = -t369 - t619;
t589 = (t212 + t627) * qJD(3) + t607 * qJD(2);
t581 = rSges(6,3) + qJ(5);
t588 = t600 * qJD(2);
t587 = (t299 * t570 - t301 * t571) * qJD(3);
t586 = (t612 * t299 - t602 * t301) * qJD(3);
t585 = t601 * qJD(2);
t584 = t626 * qJD(2);
t536 = rSges(6,1) + pkin(4);
t582 = rSges(4,2) * t307;
t423 = qJD(5) * t300;
t220 = pkin(4) * t300 - qJ(5) * t302;
t221 = rSges(6,1) * t300 - rSges(6,3) * t302;
t446 = t220 + t221;
t580 = (qJD(3) * t446 - t423) * t299;
t224 = pkin(4) * t302 + qJ(5) * t300;
t225 = rSges(6,1) * t302 + rSges(6,3) * t300;
t579 = t224 + t225;
t578 = t595 + t603;
t577 = -t585 + t586;
t576 = t587 + t588;
t343 = qJD(3) * t255;
t106 = qJD(2) * t161 - t299 * t343;
t346 = qJD(3) * t257;
t108 = qJD(2) * t163 - t299 * t346;
t575 = t590 * qJD(3) - t106 * t308 - t108 * t307 + t596 * t300 + t598 * t302;
t105 = -t301 * t343 + (-t299 * t373 + t507) * qJD(2);
t107 = -t301 * t346 + (-t258 * t299 + t513) * qJD(2);
t574 = t591 * qJD(3) + t105 * t308 + t107 * t307 + t597 * t300 - t599 * t302;
t573 = t589 * t299 + t301 * t593;
t572 = t593 * t299 - t301 * t589;
t567 = t594 * qJD(3) - t105 * t307 + t107 * t308 + t599 * t300 + t597 * t302 + t584;
t551 = qJD(2) * t146;
t566 = t595 * qJD(2) + t568 * qJD(3) + t106 * t307 - t108 * t308 - t598 * t300 + t596 * t302 - t551;
t565 = t499 + t610;
t564 = qJD(2) * t590 - t592 * t299 + t584;
t563 = -t551 - t592 * t301 + (-t299 * t627 - t591 + t625) * qJD(2);
t562 = 0.2e1 * qJD(3);
t294 = t301 * pkin(6);
t227 = pkin(2) * t299 - t294;
t306 = -qJ(4) - pkin(6);
t269 = t301 * t306;
t534 = pkin(3) * t308;
t296 = pkin(2) + t534;
t441 = -t299 * t296 - t269;
t138 = t227 + t441;
t293 = t299 * pkin(6);
t228 = t301 * pkin(2) + t293;
t249 = t301 * t296;
t388 = -t299 * t306 + t249;
t139 = t388 - t228;
t426 = qJD(3) * t301;
t427 = qJD(3) * t299;
t403 = -t138 * t427 + t139 * t426 + qJD(1);
t422 = qJD(5) * t302;
t289 = t299 * rSges(6,2);
t459 = t482 * t536 + t581 * t483 + t289;
t292 = t301 * rSges(6,2);
t460 = t299 * t579 - t292;
t23 = -t422 + (t299 * t460 + t301 * t459) * qJD(3) + t403;
t561 = qJD(3) * t23;
t559 = t300 * t536;
t134 = qJD(2) * t138;
t206 = qJD(2) * t227;
t558 = t134 - t206;
t405 = t302 * t426;
t428 = qJD(2) * t301;
t556 = rSges(6,2) * t428 + t405 * t581;
t555 = t299 ^ 2 + t301 ^ 2;
t535 = pkin(3) * t307;
t386 = -t446 - t535;
t356 = t386 * qJD(3);
t237 = t301 * t423;
t274 = qJD(4) * t299;
t442 = t237 + t274;
t554 = t301 * t356 + t442;
t461 = -t213 * t301 + t153;
t463 = -t371 * t301 + t151;
t553 = t461 + t463;
t462 = -Icges(5,2) * t486 + t152 - t242;
t464 = -t371 * t299 + t150;
t552 = t462 + t464;
t456 = -Icges(4,2) * t484 + t162 - t264;
t458 = t257 * t299 + t160;
t540 = t307 * t456 + t308 * t458;
t413 = -t139 - t459;
t273 = pkin(6) * t428;
t425 = qJD(3) * t307;
t404 = t301 * t425;
t531 = pkin(2) - t296;
t100 = -pkin(3) * t404 - t273 + t274 + (t299 * t531 - t269) * qJD(2);
t429 = qJD(2) * t299;
t251 = t299 * pkin(3) * t425;
t439 = qJD(4) * t301 + t251;
t412 = t306 * t429 + t439;
t101 = (-t301 * t531 - t293) * qJD(2) - t412;
t416 = t101 * t427 + (t100 - t134) * t426;
t180 = t221 * t299;
t529 = t224 * t428 + (-qJD(3) * t220 + t423) * t299 - qJD(3) * t180 + (t225 * t301 + t289) * qJD(2);
t336 = -t300 * t426 - t302 * t429;
t411 = t300 * t429;
t530 = t336 * t536 - t581 * t411 + t237 + t556;
t1 = (t423 + t530 * t301 + t529 * t299 + (t299 * t413 + t301 * t460) * qJD(2)) * qJD(3) + t416;
t539 = m(6) * t1;
t538 = t299 / 0.2e1;
t537 = -t301 / 0.2e1;
t532 = qJD(2) / 0.2e1;
t528 = rSges(4,1) * t308;
t527 = rSges(5,1) * t302;
t260 = rSges(4,1) * t307 + rSges(4,2) * t308;
t194 = t260 * t301;
t408 = t260 * t427;
t288 = t299 * rSges(4,3);
t165 = rSges(4,1) * t480 - rSges(4,2) * t481 + t288;
t454 = t165 + t228;
t92 = qJD(2) * t454 - t408;
t526 = t194 * t92;
t287 = t299 * rSges(5,3);
t436 = rSges(4,2) * t485 + t301 * rSges(4,3);
t164 = rSges(4,1) * t484 - t436;
t407 = t260 * t426;
t91 = -t407 + (-t164 - t227) * qJD(2);
t525 = t299 * t91;
t524 = t301 * t91;
t222 = rSges(5,1) * t300 + rSges(5,2) * t302;
t155 = rSges(5,1) * t486 - rSges(5,2) * t487 - t301 * rSges(5,3);
t395 = -t222 - t535;
t355 = t395 * t426;
t335 = t274 + t355;
t470 = t138 - t227;
t46 = (-t155 + t470) * qJD(2) + t335;
t522 = t46 * t222;
t205 = t228 * qJD(2);
t477 = -t101 - t205;
t473 = -t299 * t138 + t301 * t139;
t157 = rSges(5,1) * t482 - rSges(5,2) * t483 + t287;
t469 = -t139 - t157;
t468 = -t215 * t299 + t142;
t467 = -Icges(6,1) * t483 + t143 + t241;
t466 = -t217 * t299 - t148;
t465 = -t217 * t301 - t149;
t457 = -t257 * t301 - t161;
t455 = -t255 * t301 + t163;
t453 = -qJD(3) * t579 + t422;
t452 = -t220 * t299 - t180;
t451 = t446 * t301;
t450 = -t371 + t374;
t449 = t208 - t215;
t448 = -t213 + t218;
t447 = t217 + t372;
t444 = rSges(5,2) * t411 + rSges(5,3) * t428;
t420 = qJD(2) * qJD(4);
t443 = qJD(2) * t251 + t301 * t420;
t440 = rSges(4,3) * t428 + t429 * t582;
t438 = -t255 + t258;
t437 = t257 + t373;
t432 = qJD(2) * t210;
t431 = qJD(2) * t212;
t430 = qJD(2) * t254;
t424 = qJD(3) * t308;
t419 = pkin(3) * t481;
t418 = qJD(3) ^ 2 * t534;
t417 = t301 * t100 + t299 * t101 - t138 * t428;
t415 = pkin(3) * t424;
t195 = qJD(2) * (-pkin(2) * t429 + t273);
t414 = qJD(2) * t100 + t299 * t420 + t195;
t409 = t222 * t427;
t402 = -pkin(2) - t528;
t399 = -t427 / 0.2e1;
t396 = t426 / 0.2e1;
t226 = -rSges(5,2) * t300 + t527;
t394 = -t226 - t534;
t387 = t555 * t535;
t385 = -t534 - t579;
t203 = t226 * qJD(3);
t382 = -t203 - t415;
t47 = -t409 + (t228 - t469) * qJD(2) - t439;
t380 = t47 * t395;
t378 = t528 - t582;
t377 = -t299 * t92 - t524;
t362 = t164 * t299 + t165 * t301;
t358 = -t415 + t453;
t354 = -qJD(3) * t203 - t418;
t193 = t260 * t299;
t181 = t222 * t299;
t334 = t467 * t299 - t301 * t468;
t333 = t465 * t299 - t301 * t466;
t332 = -t307 * t455 + t308 * t457;
t331 = (t300 * t449 + t302 * t450) * qJD(2);
t330 = (-t300 * t447 + t302 * t448) * qJD(2);
t329 = (-t307 * t437 + t308 * t438) * qJD(2);
t328 = -t418 + (t422 + t453) * qJD(3);
t327 = -t296 - t579;
t111 = -rSges(4,2) * t301 * t424 + (-t308 * t429 - t404) * rSges(4,1) + t440;
t112 = -qJD(3) * t193 + (t301 * t378 + t288) * qJD(2);
t326 = t111 * t301 + t112 * t299 + (t164 * t301 - t165 * t299) * qJD(2);
t235 = t378 * qJD(3);
t185 = t222 * t301;
t96 = -qJD(3) * t181 + (t226 * t301 + t287) * qJD(2);
t94 = rSges(5,1) * t336 - rSges(5,2) * t405 + t444;
t72 = qJD(3) * t362 + qJD(1);
t57 = -t235 * t426 + (-t112 - t205 + t408) * qJD(2);
t56 = -t235 * t427 + t195 + (t111 - t407) * qJD(2);
t41 = (t155 * t299 + t157 * t301) * qJD(3) + t403;
t38 = -t580 + (t228 - t413) * qJD(2) - t439;
t37 = (-t460 + t470) * qJD(2) + t554;
t36 = t326 * qJD(3);
t25 = t354 * t301 + (-t96 + t409 + t477) * qJD(2) + t443;
t24 = t354 * t299 + (t94 + t355) * qJD(2) + t414;
t12 = t328 * t301 + (t477 - t529 + t580) * qJD(2) + t443;
t11 = t328 * t299 + ((t356 + t423) * t301 + t530) * qJD(2) + t414;
t2 = (t299 * t96 + t301 * t94 + (t155 * t301 + t299 * t469) * qJD(2)) * qJD(3) + t416;
t3 = [m(4) * t36 + m(5) * t2 + t539; (-t607 * qJD(3) + t233 * t308 + t234 * t307 + t616 * t300 - t617 * t302) * qJD(2) + (-(-qJD(2) * t460 - t37 + t554 + t558) * t38 + t12 * (t292 + t441) + t37 * t412 + t11 * (t249 + t459) + t38 * (t442 + t556) + (t38 * (-t535 - t559) * qJD(3) + (-t38 * t306 + t327 * t37) * qJD(2)) * t301 + (-t11 * t306 - t12 * t536 * t302 + (-t37 * qJD(5) - t12 * t581) * t300 + t37 * (-t302 * t581 + t559) * qJD(3) + (-t37 * rSges(6,2) + t327 * t38) * qJD(2)) * t299) * m(6) + (-(-qJD(2) * t155 + t335 - t46 + t558) * t47 + t25 * (-t155 + t441) + t46 * t412 + t24 * (t157 + t388) + t47 * (t274 + t444) + (t299 * t522 + t301 * t380) * qJD(3) + ((-t46 * rSges(5,3) + t47 * (-t296 - t527)) * t299 + (t46 * (-t226 - t296) - t47 * t306) * t301) * qJD(2)) * m(5) + (-(-qJD(2) * t164 - t206 - t407 - t91) * t92 + t57 * (t299 * t402 + t294 + t436) + t56 * t454 + t92 * (t273 + t440) + (t260 * t525 - t526) * qJD(3) + ((-pkin(2) - t378) * t524 + (t91 * (-rSges(4,3) - pkin(6)) + t92 * t402) * t299) * qJD(2)) * m(4) + (((t48 - t560 + t565) * t299 + ((t604 + t621) * t301 + t583 + t609 + t620) * t301) * qJD(3) + t588) * t396 + (t573 + t574) * t427 / 0.2e1 + (((t301 * t578 - t565 + t570) * t301 + (t299 * t578 - t132 + t381 + t571 - t611) * t299) * qJD(3) + t577 + t585) * t399 - (t572 - t575 + t576) * t426 / 0.2e1 + ((t568 - t601) * t299 + (-t594 + t600) * t301) * qJD(3) * t532; (-(t193 * t91 - t526) * qJD(2) - (t72 * (-t193 * t299 - t194 * t301) + t377 * t378) * qJD(3) + t36 * t362 + t72 * t326 + t377 * t235 + (-t56 * t299 - t57 * t301 + (-t301 * t92 + t525) * qJD(2)) * t260) * m(4) - (((t299 * t455 - t301 * t456) * t308 + (t299 * t457 + t301 * t458) * t307 + (t299 * t553 - t301 * t552) * t302 + (t334 + t333) * t300) * qJD(3) + (t307 * t438 + t308 * t437 + (t447 - t449) * t302 + (t448 + t450) * t300) * qJD(2)) * qJD(2) / 0.2e1 + (t575 * t301 + t574 * t299 + (t299 * t568 - t301 * t594) * qJD(2)) * t532 + ((t605 * t427 + t430 + t431 + t432) * t299 + (t331 + t330 + t329 + (((-t466 - t468) * t302 + t552 * t300 + t540) * t301 + (t332 + (t465 + t467) * t302 - t553 * t300 + t606) * t299) * qJD(3)) * t301) * t399 + ((-t426 * t491 - t431) * t301 + (t331 + (t490 * t301 + t334 * t302 + (-t299 * t463 + t301 * t464) * t300) * qJD(3)) * t299 + (-t426 * t493 - t432) * t301 + (t330 + (t492 * t301 + t333 * t302 + (-t299 * t461 + t301 * t462) * t300) * qJD(3)) * t299 + (-t426 * t489 - t430) * t301 + (t329 + (t332 * t299 + (t488 + t540) * t301) * qJD(3)) * t299) * t396 + (-(t23 * t300 + (t299 * t38 + t301 * t37) * t302) * qJD(5) - (-t37 * t452 + t38 * (-t419 - t451)) * qJD(2) - (-t23 * t387 + (-t23 * t451 + t37 * t385) * t301 + (t23 * t452 + t38 * t385) * t299) * qJD(3) + t1 * t473 + t23 * t417 + (t12 * t386 + t37 * t358 + t1 * t459 + t23 * t530 + (t23 * t460 + t38 * t386) * qJD(2)) * t301 + (t11 * t386 + t38 * t358 + t1 * t460 + t23 * t529 + (t23 * t413 + t37 * t446) * qJD(2)) * t299) * m(6) + (-(t46 * t181 + t47 * (-t185 - t419)) * qJD(2) - (-t41 * t387 + (-t41 * t185 + t394 * t46) * t301 + (-t41 * t181 + t394 * t47) * t299) * qJD(3) + t2 * t473 + t41 * t417 + (t25 * t395 + t46 * t382 + t2 * t157 + t41 * t94 + (t41 * t155 + t380) * qJD(2)) * t301 + (t24 * t395 + t47 * t382 + t2 * t155 + t41 * t96 + (t41 * t469 + t522) * qJD(2)) * t299) * m(5) + (t573 * qJD(2) + ((t570 * qJD(2) + t566 * t301) * t301 + (t563 * t299 + t571 * qJD(2) + (-t564 + t567) * t301) * t299) * t562) * t538 + (t572 * qJD(2) + ((t612 * qJD(2) + t564 * t301) * t301 + (t567 * t299 + t602 * qJD(2) + (-t563 + t566) * t301) * t299) * t562) * t537 + (t577 + t586) * t429 / 0.2e1 + (t576 + t587) * t428 / 0.2e1; 0.2e1 * (t11 * t537 + t12 * t538) * m(6) + 0.2e1 * (t24 * t537 + t25 * t538) * m(5); -t302 * t539 + 0.2e1 * (m(6) * (t11 * t299 + t12 * t301 + t561) / 0.2e1 - m(6) * t555 * t561 / 0.2e1) * t300;];
tauc = t3(:);
