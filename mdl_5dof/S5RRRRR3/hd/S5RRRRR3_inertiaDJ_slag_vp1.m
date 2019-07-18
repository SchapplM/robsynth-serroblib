% Calculate time derivative of joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR3_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR3_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR3_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:17:11
% EndTime: 2019-07-18 17:17:47
% DurationCPUTime: 16.15s
% Computational Cost: add. (54507->913), mult. (58419->1274), div. (0->0), fcn. (56582->10), ass. (0->480)
t411 = qJ(4) + qJ(5);
t403 = cos(t411);
t412 = qJ(2) + qJ(3);
t404 = cos(t412);
t418 = cos(qJ(1));
t401 = sin(t411);
t415 = sin(qJ(1));
t580 = t415 * t401;
t336 = -t403 * t418 - t404 * t580;
t579 = t415 * t403;
t337 = -t401 * t418 + t404 * t579;
t402 = sin(t412);
t589 = t402 * t415;
t243 = Icges(6,5) * t337 + Icges(6,6) * t336 + Icges(6,3) * t589;
t245 = Icges(6,4) * t337 + Icges(6,2) * t336 + Icges(6,6) * t589;
t247 = Icges(6,1) * t337 + Icges(6,4) * t336 + Icges(6,5) * t589;
t462 = -t245 * t401 + t247 * t403;
t122 = -t243 * t404 + t402 * t462;
t474 = Icges(6,5) * t403 - Icges(6,6) * t401;
t309 = -Icges(6,3) * t404 + t402 * t474;
t609 = Icges(6,4) * t403;
t478 = -Icges(6,2) * t401 + t609;
t310 = -Icges(6,6) * t404 + t402 * t478;
t610 = Icges(6,4) * t401;
t483 = Icges(6,1) * t403 - t610;
t311 = -Icges(6,5) * t404 + t402 * t483;
t162 = t309 * t589 + t310 * t336 + t311 * t337;
t663 = -t162 - t122;
t586 = t404 * t418;
t338 = -t401 * t586 + t579;
t339 = t403 * t586 + t580;
t588 = t402 * t418;
t244 = Icges(6,5) * t339 + Icges(6,6) * t338 + Icges(6,3) * t588;
t246 = Icges(6,4) * t339 + Icges(6,2) * t338 + Icges(6,6) * t588;
t248 = Icges(6,1) * t339 + Icges(6,4) * t338 + Icges(6,5) * t588;
t461 = -t246 * t401 + t248 * t403;
t123 = -t244 * t404 + t402 * t461;
t163 = t309 * t588 + t338 * t310 + t339 * t311;
t662 = -t163 - t123;
t554 = qJD(1) * t415;
t529 = t402 * t554;
t407 = qJD(4) + qJD(5);
t508 = -t404 * t407 + qJD(1);
t449 = t508 * t418;
t555 = qJD(1) * t404;
t507 = -t407 + t555;
t408 = qJD(2) + qJD(3);
t583 = t408 * t418;
t539 = t402 * t583;
t644 = t415 * t507 + t539;
t220 = t401 * t644 + t403 * t449;
t221 = t401 * t449 - t403 * t644;
t537 = t404 * t583;
t536 = t221 * rSges(6,1) + t220 * rSges(6,2) + rSges(6,3) * t537;
t148 = -rSges(6,3) * t529 + t536;
t413 = sin(qJ(4));
t544 = pkin(3) * qJD(4) * t418;
t582 = t413 * t418;
t548 = pkin(3) * t582;
t416 = cos(qJ(4));
t549 = qJD(4) * t416;
t560 = t415 * pkin(3) * t549 + qJD(1) * t548;
t399 = pkin(3) * t416 + pkin(2);
t627 = pkin(2) - t399;
t661 = t148 + t627 * t539 + (-t413 * t544 + t554 * t627) * t404 + t560;
t489 = -t337 * rSges(6,1) - t336 * rSges(6,2);
t251 = rSges(6,3) * t589 - t489;
t518 = t404 * t627;
t304 = -t415 * t518 - t548;
t660 = t251 + t304;
t409 = t415 ^ 2;
t410 = t418 ^ 2;
t558 = t409 + t410;
t587 = t404 * t408;
t521 = t587 / 0.2e1;
t659 = -t418 * t521 + t529 / 0.2e1;
t553 = qJD(1) * t418;
t516 = t553 / 0.2e1;
t658 = -t402 * t516 - t415 * t521;
t585 = t408 * t415;
t538 = t404 * t585;
t433 = t402 * t553 + t538;
t414 = sin(qJ(2));
t417 = cos(qJ(2));
t495 = rSges(3,1) * t417 - rSges(3,2) * t414;
t657 = t418 * rSges(3,3) - t415 * t495;
t615 = Icges(3,4) * t417;
t482 = -Icges(3,2) * t414 + t615;
t348 = Icges(3,6) * t415 + t418 * t482;
t616 = Icges(3,4) * t414;
t487 = Icges(3,1) * t417 - t616;
t350 = Icges(3,5) * t415 + t418 * t487;
t452 = t348 * t414 - t350 * t417;
t656 = t415 * t452;
t613 = Icges(4,4) * t404;
t480 = -Icges(4,2) * t402 + t613;
t325 = Icges(4,6) * t415 + t418 * t480;
t614 = Icges(4,4) * t402;
t485 = Icges(4,1) * t404 - t614;
t327 = Icges(4,5) * t415 + t418 * t485;
t454 = t325 * t402 - t327 * t404;
t655 = t415 * t454;
t347 = -Icges(3,6) * t418 + t415 * t482;
t349 = -Icges(3,5) * t418 + t415 * t487;
t453 = t347 * t414 - t349 * t417;
t653 = t418 * t453;
t324 = -Icges(4,6) * t418 + t415 * t480;
t326 = -Icges(4,5) * t418 + t415 * t485;
t455 = t324 * t402 - t326 * t404;
t652 = t418 * t455;
t252 = t339 * rSges(6,1) + t338 * rSges(6,2) + rSges(6,3) * t588;
t651 = -t415 * t251 - t418 * t252;
t475 = Icges(5,5) * t416 - Icges(5,6) * t413;
t239 = t475 * t587 + (Icges(5,3) * t408 + (-Icges(5,5) * t413 - Icges(5,6) * t416) * qJD(4)) * t402;
t611 = Icges(5,4) * t416;
t479 = -Icges(5,2) * t413 + t611;
t317 = -Icges(5,6) * t404 + t402 * t479;
t600 = t317 * t413;
t650 = -t408 * t600 - t239;
t576 = t416 * t418;
t578 = t415 * t413;
t359 = -t404 * t578 - t576;
t577 = t415 * t416;
t360 = t404 * t577 - t582;
t492 = -rSges(5,1) * t360 - rSges(5,2) * t359;
t279 = rSges(5,3) * t589 - t492;
t361 = -t404 * t582 + t577;
t362 = t404 * t576 + t578;
t280 = t362 * rSges(5,1) + t361 * rSges(5,2) + rSges(5,3) * t588;
t649 = -t415 * t279 - t418 * t280;
t551 = qJD(2) * t417;
t648 = -t414 * t553 - t415 * t551;
t620 = t415 * rSges(3,3);
t352 = t418 * t495 + t620;
t647 = -t352 * t415 - t418 * t657;
t476 = Icges(4,5) * t404 - Icges(4,6) * t402;
t322 = -Icges(4,3) * t418 + t415 * t476;
t646 = qJD(1) * t322;
t477 = Icges(3,5) * t417 - Icges(3,6) * t414;
t345 = -Icges(3,3) * t418 + t415 * t477;
t505 = -qJD(4) + t555;
t645 = t415 * t505 + t539;
t540 = t402 * t585;
t643 = t418 * t507 - t540;
t364 = Icges(4,2) * t404 + t614;
t365 = Icges(4,1) * t402 + t613;
t451 = t364 * t402 - t365 * t404;
t642 = qJD(1) * t451 + t476 * t408;
t641 = 2 * m(3);
t640 = 2 * m(4);
t639 = 2 * m(5);
t638 = 2 * m(6);
t637 = -t404 / 0.2e1;
t636 = t415 / 0.2e1;
t635 = -t418 / 0.2e1;
t634 = -rSges(5,3) - pkin(5);
t633 = -rSges(6,3) - pkin(5);
t381 = rSges(3,1) * t414 + rSges(3,2) * t417;
t632 = m(3) * t381;
t367 = rSges(4,1) * t402 + rSges(4,2) * t404;
t631 = m(4) * t367;
t630 = pkin(1) * t414;
t629 = pkin(1) * t417;
t628 = pkin(2) * t404;
t625 = rSges(4,1) * t404;
t623 = rSges(4,3) * t418;
t450 = t508 * t415;
t222 = -t401 * t643 + t403 * t450;
t223 = t401 * t450 + t403 * t643;
t143 = Icges(6,5) * t223 + Icges(6,6) * t222 + Icges(6,3) * t433;
t145 = Icges(6,4) * t223 + Icges(6,2) * t222 + Icges(6,6) * t433;
t147 = Icges(6,1) * t223 + Icges(6,4) * t222 + Icges(6,5) * t433;
t36 = (t408 * t462 - t143) * t404 + (t243 * t408 + (-t245 * t407 + t147) * t403 + (-t247 * t407 - t145) * t401) * t402;
t622 = t36 * t418;
t432 = -t529 + t537;
t142 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t432;
t144 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t432;
t146 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t432;
t37 = (t408 * t461 - t142) * t404 + (t244 * t408 + (-t246 * t407 + t146) * t403 + (-t248 * t407 - t144) * t401) * t402;
t621 = t37 * t415;
t405 = t415 * rSges(4,3);
t506 = -qJD(4) * t404 + qJD(1);
t257 = t506 * t577 + (-t418 * t505 + t540) * t413;
t584 = t408 * t416;
t258 = t505 * t576 + (-t402 * t584 + t413 * t506) * t415;
t167 = Icges(5,5) * t258 + Icges(5,6) * t257 + Icges(5,3) * t433;
t169 = Icges(5,4) * t258 + Icges(5,2) * t257 + Icges(5,6) * t433;
t171 = Icges(5,1) * t258 + Icges(5,4) * t257 + Icges(5,5) * t433;
t273 = Icges(5,5) * t360 + Icges(5,6) * t359 + Icges(5,3) * t589;
t275 = Icges(5,4) * t360 + Icges(5,2) * t359 + Icges(5,6) * t589;
t277 = Icges(5,1) * t360 + Icges(5,4) * t359 + Icges(5,5) * t589;
t460 = -t275 * t413 + t277 * t416;
t45 = (t408 * t460 - t167) * t404 + (-t169 * t413 + t171 * t416 + t273 * t408 + (-t275 * t416 - t277 * t413) * qJD(4)) * t402;
t618 = t45 * t418;
t447 = t506 * t418;
t255 = t413 * t645 + t416 * t447;
t256 = t413 * t447 - t416 * t645;
t166 = Icges(5,5) * t256 + Icges(5,6) * t255 + Icges(5,3) * t432;
t168 = Icges(5,4) * t256 + Icges(5,2) * t255 + Icges(5,6) * t432;
t170 = Icges(5,1) * t256 + Icges(5,4) * t255 + Icges(5,5) * t432;
t274 = Icges(5,5) * t362 + Icges(5,6) * t361 + Icges(5,3) * t588;
t276 = Icges(5,4) * t362 + Icges(5,2) * t361 + Icges(5,6) * t588;
t278 = Icges(5,1) * t362 + Icges(5,4) * t361 + Icges(5,5) * t588;
t459 = -t276 * t413 + t278 * t416;
t46 = (t408 * t459 - t166) * t404 + (-t168 * t413 + t170 * t416 + t274 * t408 + (-t276 * t416 - t278 * t413) * qJD(4)) * t402;
t617 = t46 * t415;
t612 = Icges(5,4) * t413;
t240 = t479 * t587 + (Icges(5,6) * t408 + (-Icges(5,2) * t416 - t612) * qJD(4)) * t402;
t602 = t240 * t413;
t601 = t311 * t403;
t599 = t347 * t417;
t598 = t348 * t417;
t597 = t349 * t414;
t596 = t350 * t414;
t593 = t364 * t408;
t592 = t365 * t408;
t591 = t402 * t407;
t590 = t402 * t408;
t494 = -rSges(4,2) * t402 + t625;
t343 = t494 * t408;
t581 = t415 * t343;
t575 = t417 * t418;
t490 = t223 * rSges(6,1) + t222 * rSges(6,2);
t149 = rSges(6,3) * t433 + t490;
t574 = t149 * t588 + t251 * t537;
t488 = rSges(6,1) * t403 - rSges(6,2) * t401;
t219 = t488 * t587 + (rSges(6,3) * t408 + (-rSges(6,1) * t401 - rSges(6,2) * t403) * t407) * t402;
t550 = qJD(4) * t413;
t524 = t402 * t550;
t298 = -pkin(3) * t524 - t408 * t518;
t572 = -t219 - t298;
t312 = -rSges(6,3) * t404 + t402 * t488;
t301 = t312 * t554;
t571 = t252 * t590 + t402 * t301;
t491 = rSges(5,1) * t416 - rSges(5,2) * t413;
t242 = t491 * t587 + (rSges(5,3) * t408 + (-rSges(5,1) * t413 - rSges(5,2) * t416) * qJD(4)) * t402;
t497 = pkin(5) * t402 + t628;
t344 = t497 * t408;
t570 = -t242 - t344;
t386 = pkin(2) * t586;
t396 = pkin(3) * t578;
t562 = t399 * t586 + t396;
t305 = -t386 + t562;
t568 = -t252 - t305;
t385 = pkin(5) * t588;
t356 = t386 + t385;
t567 = -t280 - t356;
t186 = t404 * t251 + t312 * t589;
t319 = -rSges(5,3) * t404 + t402 * t491;
t308 = t319 * t554;
t368 = pkin(2) * t402 - pkin(5) * t404;
t358 = t368 * t554;
t566 = t308 + t358;
t357 = t627 * t402;
t565 = -t312 + t357;
t328 = t415 * t494 - t623;
t329 = rSges(4,1) * t586 - rSges(4,2) * t588 + t405;
t259 = t415 * t328 + t418 * t329;
t564 = -t319 - t368;
t355 = t497 * t415;
t563 = t415 * t355 + t418 * t356;
t561 = rSges(4,2) * t529 + rSges(4,3) * t553;
t559 = t558 * t629;
t323 = Icges(4,3) * t415 + t418 * t476;
t557 = qJD(1) * t323;
t346 = Icges(3,3) * t415 + t418 * t477;
t556 = qJD(1) * t346;
t552 = qJD(2) * t414;
t547 = pkin(1) * t552;
t546 = pkin(1) * t551;
t545 = pkin(3) * t550;
t543 = t310 * t587;
t541 = t399 * t590;
t535 = -t344 + t572;
t484 = Icges(5,1) * t416 - t612;
t241 = t484 * t587 + (Icges(5,5) * t408 + (-Icges(5,1) * t413 - t611) * qJD(4)) * t402;
t316 = -Icges(5,3) * t404 + t402 * t475;
t318 = -Icges(5,5) * t404 + t402 * t484;
t534 = t402 * t416 * t241 + t404 * t318 * t584 + t316 * t590;
t533 = t256 * rSges(5,1) + t255 * rSges(5,2) + rSges(5,3) * t537;
t375 = pkin(2) * t540;
t376 = pkin(5) * t537;
t434 = -t404 * t554 - t539;
t532 = t415 * (pkin(5) * t433 + qJD(1) * t386 - t375) + t418 * (pkin(2) * t434 - pkin(5) * t529 + t376) + t355 * t553;
t333 = t357 * t554;
t531 = t301 - t333 + t358;
t530 = -t368 + t565;
t526 = t415 * t552;
t523 = t589 / 0.2e1;
t522 = t588 / 0.2e1;
t138 = -t273 * t404 + t402 * t460;
t180 = t316 * t589 + t317 * t359 + t318 * t360;
t520 = t138 / 0.2e1 + t180 / 0.2e1;
t139 = -t274 * t404 + t402 * t459;
t181 = t316 * t588 + t361 * t317 + t362 * t318;
t519 = t139 / 0.2e1 + t181 / 0.2e1;
t517 = t554 / 0.2e1;
t515 = -t367 - t630;
t514 = -t368 - t630;
t513 = t565 * t418;
t282 = t564 * t418;
t262 = -qJD(1) * t324 - t418 * t593;
t512 = t327 * t408 + t262;
t263 = qJD(1) * t325 - t415 * t593;
t511 = t326 * t408 + t263;
t264 = -qJD(1) * t326 - t418 * t592;
t510 = -t325 * t408 + t264;
t265 = qJD(1) * t327 - t415 * t592;
t509 = t324 * t408 - t265;
t504 = t404 * t149 + t219 * t589 + t312 * t433;
t161 = t563 - t649;
t499 = -t319 + t514;
t498 = -t344 - t546;
t211 = t530 * t418;
t182 = -t309 * t404 + (-t310 * t401 + t601) * t402;
t468 = t122 * t415 + t123 * t418;
t25 = t143 * t588 + t338 * t145 + t339 * t147 + t220 * t245 + t221 * t247 + t243 * t432;
t26 = t142 * t588 + t338 * t144 + t339 * t146 + t220 * t246 + t221 * t248 + t244 * t432;
t107 = t243 * t588 + t338 * t245 + t339 * t247;
t108 = t244 * t588 + t338 * t246 + t339 * t248;
t470 = t107 * t415 + t108 * t418;
t471 = t107 * t418 - t108 * t415;
t212 = t474 * t587 + (Icges(6,3) * t408 + (-Icges(6,5) * t401 - Icges(6,6) * t403) * t407) * t402;
t213 = t478 * t587 + (Icges(6,6) * t408 + (-Icges(6,2) * t403 - t610) * t407) * t402;
t214 = t483 * t587 + (Icges(6,5) * t408 + (-Icges(6,1) * t401 - t609) * t407) * t402;
t62 = t212 * t588 + t338 * t213 + t339 * t214 + t220 * t310 + t221 * t311 + t309 * t432;
t5 = (t408 * t470 - t62) * t404 + (qJD(1) * t471 + t163 * t408 + t25 * t415 + t26 * t418) * t402;
t105 = t243 * t589 + t245 * t336 + t247 * t337;
t106 = t244 * t589 + t246 * t336 + t248 * t337;
t472 = t105 * t415 + t106 * t418;
t55 = -t162 * t404 + t402 * t472;
t56 = -t163 * t404 + t402 * t470;
t27 = t143 * t589 + t336 * t145 + t337 * t147 + t222 * t245 + t223 * t247 + t243 * t433;
t28 = t142 * t589 + t336 * t144 + t337 * t146 + t222 * t246 + t223 * t248 + t244 * t433;
t473 = t105 * t418 - t106 * t415;
t63 = t212 * t589 + t336 * t213 + t337 * t214 + t222 * t310 + t223 * t311 + t309 * t433;
t6 = (t408 * t472 - t63) * t404 + (qJD(1) * t473 + t162 * t408 + t27 * t415 + t28 * t418) * t402;
t496 = t5 * t588 + t56 * t537 + t6 * t589 + (-t182 * t404 + t402 * t468) * t590 + t433 * t55;
t493 = t258 * rSges(5,1) + t257 * rSges(5,2);
t486 = Icges(3,1) * t414 + t615;
t481 = Icges(3,2) * t417 + t616;
t363 = Icges(4,5) * t402 + Icges(4,6) * t404;
t469 = t122 * t418 - t123 * t415;
t126 = t273 * t589 + t275 * t359 + t277 * t360;
t127 = t274 * t589 + t276 * t359 + t278 * t360;
t467 = t126 * t418 - t127 * t415;
t466 = t126 * t415 + t127 * t418;
t128 = t273 * t588 + t361 * t275 + t362 * t277;
t129 = t274 * t588 + t361 * t276 + t362 * t278;
t465 = t128 * t418 - t129 * t415;
t464 = t128 * t415 + t129 * t418;
t463 = t138 * t415 + t139 * t418;
t458 = t279 * t418 - t280 * t415;
t448 = t514 + t565;
t446 = t558 * t547;
t445 = -t242 + t498;
t268 = t499 * t418;
t102 = t415 * t304 + t418 * t305 + t563 - t651;
t444 = t367 * t408;
t443 = t498 + t572;
t440 = t408 * t363;
t439 = t304 * t418 + t415 * t568;
t136 = -t329 * t554 + t415 * (-t415 * t444 + (t418 * t494 + t405) * qJD(1)) + t418 * (rSges(4,1) * t434 - rSges(4,2) * t537 + t561) + t328 * t553;
t438 = qJD(2) * t486;
t437 = qJD(2) * t481;
t436 = qJD(2) * (-Icges(3,5) * t414 - Icges(3,6) * t417);
t202 = t448 * t418;
t435 = -t494 - t629;
t15 = qJD(1) * t470 - t25 * t418 + t26 * t415;
t190 = -t322 * t418 - t415 * t455;
t191 = -t323 * t418 - t655;
t192 = t415 * t322 - t652;
t193 = t415 * t323 - t418 * t454;
t39 = t167 * t588 + t361 * t169 + t362 * t171 + t255 * t275 + t256 * t277 + t273 * t432;
t40 = t166 * t588 + t361 * t168 + t362 * t170 + t255 * t276 + t256 * t278 + t274 * t432;
t21 = qJD(1) * t464 - t39 * t418 + t40 * t415;
t260 = -t418 * t440 - t646;
t261 = -t415 * t440 + t557;
t431 = (-t190 * t418 - t467 - t473) * t554 + (-t192 * t418 - t465 - t471) * t553 + (t15 + t21 + t191 * t554 + t193 * t553 + (t193 * qJD(1) + (t263 * t402 - t265 * t404 + t324 * t587 + t326 * t590 - t646) * t418) * t418 + ((t192 + t655) * qJD(1) + (-t261 + t510 * t404 - t512 * t402 + (t323 - t455) * qJD(1)) * t418 + t415 * t260) * t415) * t415;
t430 = t402 * t634 - t628 - t629;
t429 = -t399 * t404 + t402 * t633 - t629;
t428 = qJD(1) * t430;
t16 = qJD(1) * t472 - t27 * t418 + t28 * t415;
t427 = t15 * t522 + t16 * t523 + t5 * t636 + t6 * t635 + (qJD(1) * t468 + t621 - t622) * t637 + t55 * t517 + t56 * t516 - t469 * t590 / 0.2e1 + t659 * t471 + t658 * t473;
t426 = t429 * t415;
t425 = -t404 * t212 + t309 * t590 + t587 * t601 + (t214 * t402 - t310 * t591) * t403;
t172 = -rSges(5,3) * t529 + t533;
t173 = rSges(5,3) * t433 + t493;
t67 = t418 * t172 + t415 * t173 + t279 * t553 + t554 * t567 + t532;
t179 = t182 * t590;
t83 = (-t543 + (-t311 * t407 - t213) * t402) * t401 + t425;
t7 = t179 + (t408 * t468 - t83) * t404 + (qJD(1) * t469 + t36 * t415 + t37 * t418) * t402;
t424 = -t404 * t7 - t529 * t56 + t496;
t423 = t179 + (t36 + t63) * t523 + (t37 + t62) * t522 + t662 * t659 + t663 * t658;
t41 = t167 * t589 + t359 * t169 + t360 * t171 + t257 * t275 + t258 * t277 + t273 * t433;
t42 = t166 * t589 + t359 * t168 + t360 * t170 + t257 * t276 + t258 * t278 + t274 * t433;
t22 = qJD(1) * t466 - t41 * t418 + t42 * t415;
t30 = (t418 * t261 + (t191 + t652) * qJD(1)) * t418 + (t190 * qJD(1) + (-t262 * t402 + t264 * t404 - t325 * t587 - t327 * t590 + t557) * t415 + (-t260 + t509 * t404 + t511 * t402 + (-t322 - t454) * qJD(1)) * t418) * t415;
t422 = (-t16 - t22 - t30) * t418 + t431;
t341 = t480 * t408;
t342 = t485 * t408;
t421 = qJD(1) * t363 + (t342 - t593) * t404 + (-t341 - t592) * t402;
t199 = -t399 * t540 + t375 - t518 * t553 + (qJD(4) * t359 + t413 * t554) * pkin(3);
t38 = (-t356 + t568) * t554 + t532 + t660 * t553 + t661 * t418 + (t149 + t199) * t415;
t81 = t239 * t588 + t361 * t240 + t362 * t241 + t255 * t317 + t256 * t318 + t316 * t432;
t10 = (t408 * t464 - t81) * t404 + (qJD(1) * t465 + t408 * t181 + t39 * t415 + t40 * t418) * t402;
t82 = t239 * t589 + t359 * t240 + t360 * t241 + t257 * t317 + t258 * t318 + t316 * t433;
t11 = (t408 * t466 - t82) * t404 + (qJD(1) * t467 + t408 * t180 + t41 * t415 + t418 * t42) * t402;
t68 = -t180 * t404 + t402 * t466;
t69 = -t181 * t404 + t402 * t464;
t420 = t10 * t636 + t11 * t635 + t21 * t522 + t22 * t523 + (qJD(1) * t463 + t617 - t618) * t637 + t427 + t68 * t517 + t69 * t516 + (-t138 * t418 + t139 * t415) * t590 / 0.2e1 + t659 * t465 + t658 * t467;
t419 = -t622 / 0.2e1 + t621 / 0.2e1 - t618 / 0.2e1 + t617 / 0.2e1 + (t402 * t510 + t404 * t512 + t415 * t642 + t421 * t418 + t62 + t81) * t636 + (-t402 * t509 + t404 * t511 + t421 * t415 - t418 * t642 + t63 + t82) * t635 + (t324 * t404 + t326 * t402 - t363 * t418 - t415 * t451 + t138 + t180 - t663) * t517 + (t325 * t404 + t327 * t402 + t415 * t363 - t418 * t451 + t139 + t181 - t662) * t516;
t398 = pkin(1) * t575;
t392 = t554 * t630;
t391 = pkin(1) * t526;
t374 = t495 * qJD(2);
t321 = t515 * t418;
t320 = t515 * t415;
t307 = t329 + t398;
t306 = t415 * t435 + t623;
t291 = -rSges(3,1) * t526 + (rSges(3,1) * t575 + t620) * qJD(1) + t648 * rSges(3,2);
t290 = -t381 * t418 * qJD(2) + qJD(1) * t657;
t285 = t415 * t436 + t556;
t284 = -qJD(1) * t345 + t418 * t436;
t281 = t564 * t415;
t267 = t499 * t415;
t236 = pkin(1) * t648 - t367 * t553 - t581;
t235 = t367 * t554 + t392 + (-t343 - t546) * t418;
t232 = t251 * t588;
t224 = t559 + t259;
t218 = t391 + t367 * t585 + (t418 * t435 - t405) * qJD(1);
t217 = (-t625 - t629) * t554 + (-t444 - t547) * t418 + t561;
t210 = t530 * t415;
t209 = t398 - t567;
t208 = t415 * t430 + t492;
t207 = t415 * t346 - t418 * t452;
t206 = t415 * t345 - t653;
t205 = -t346 * t418 - t656;
t204 = -t345 * t418 - t415 * t453;
t201 = t448 * t415;
t197 = -t404 * t280 - t319 * t588;
t196 = t279 * t404 + t319 * t589;
t189 = t398 + t385 + t252 + t562;
t188 = t426 + t489 + t548;
t187 = -t404 * t252 - t312 * t588;
t185 = -t316 * t404 + (t318 * t416 - t600) * t402;
t184 = t185 * t590;
t183 = t458 * t402;
t176 = -t252 * t589 + t232;
t175 = qJD(1) * t282 + t415 * t570;
t174 = t418 * t570 + t566;
t158 = qJD(1) * t268 + t415 * t445;
t157 = t418 * t445 + t392 + t566;
t152 = t402 * t513 + t404 * t568;
t151 = t304 * t404 - t357 * t589 + t186;
t150 = t161 + t559;
t132 = -t446 + t136;
t121 = t418 * t428 + t538 * t634 + t375 + t391 - t493;
t120 = t376 + (-pkin(2) * t590 - t547) * t418 + t415 * t428 + t533;
t115 = t402 * t439 + t232;
t110 = qJD(1) * t211 + t415 * t535;
t109 = t418 * t535 + t531;
t104 = qJD(1) * t202 + t415 * t443;
t103 = t418 * t443 + t392 + t531;
t101 = t102 + t559;
t100 = t416 * t544 + t391 + (t541 + (t408 * t633 + t545) * t404) * t415 + (t418 * t429 - t396) * qJD(1) - t490;
t99 = t376 + (-t404 * t545 - t541 - t547) * t418 + qJD(1) * t426 + t536 + t560;
t98 = (t319 * t585 + t173) * t404 + (t415 * t242 - t279 * t408 + t319 * t553) * t402;
t97 = (-t319 * t583 - t172) * t404 + (-t242 * t418 + t280 * t408 + t308) * t402;
t95 = t650 * t404 + (-t602 + (-t317 * t416 - t318 * t413) * qJD(4)) * t402 + t534;
t92 = -t251 * t590 + t504;
t91 = -t219 * t588 + (-t312 * t583 - t148) * t404 + t571;
t70 = t458 * t587 + (qJD(1) * t649 - t172 * t415 + t173 * t418) * t402;
t59 = -t446 + t67;
t54 = -t252 * t538 + (qJD(1) * t651 - t148 * t415) * t402 + t574;
t48 = (-t357 * t585 + t199) * t404 + (t415 * t298 - t357 * t553 - t408 * t660) * t402 + t504;
t47 = (t408 * t513 - t661) * t404 + (t305 * t408 + t418 * t572 - t333) * t402 + t571;
t35 = -t446 + t38;
t24 = t439 * t587 + (t199 * t418 - t661 * t415 + (-t415 * t660 + t418 * t568) * qJD(1)) * t402 + t574;
t1 = [t534 + t425 + (t100 * t188 + t189 * t99) * t638 + (t120 * t209 + t121 * t208) * t639 + (t217 * t307 + t218 * t306) * t640 + (t290 * t352 - t291 * t657) * t641 - t318 * t524 + t365 * t587 - t364 * t590 + (t487 - t481) * t552 + (t486 + t482) * t551 + (-t311 * t591 - t543) * t401 + (t341 + t650) * t404 + (-t213 * t401 - t317 * t549 + t342 - t602) * t402; m(3) * ((-t290 * t415 + t291 * t418) * t381 + t647 * t374) + (t410 / 0.2e1 + t409 / 0.2e1) * t477 * qJD(2) + m(6) * (t100 * t202 + t103 * t188 + t104 * t189 + t201 * t99) + m(5) * (t120 * t267 + t121 * t268 + t157 * t208 + t158 * t209) + m(4) * (t217 * t320 + t218 * t321 + t235 * t306 + t236 * t307) + ((t598 / 0.2e1 + t596 / 0.2e1 - t352 * t632) * t418 + (t657 * t632 + t599 / 0.2e1 + t597 / 0.2e1) * t415) * qJD(1) + (-qJD(2) * t453 + (qJD(1) * t348 - t415 * t437) * t417 + (qJD(1) * t350 - t415 * t438) * t414) * t635 + (-qJD(2) * t452 + (-qJD(1) * t347 - t418 * t437) * t417 + (-qJD(1) * t349 - t418 * t438) * t414) * t636 + t419; t431 - t418 * t16 - t418 * t22 - t418 * t30 + (t101 * t35 + t103 * t202 + t104 * t201) * t638 + (t150 * t59 + t157 * t268 + t158 * t267) * t639 + (t132 * t224 + t235 * t321 + t236 * t320) * t640 + ((t352 * t418 - t415 * t657) * (qJD(1) * t647 + t418 * t290 + t415 * t291) + t558 * t381 * t374) * t641 + (-t206 * t418 + t207 * t415) * t553 + (-t204 * t418 + t205 * t415) * t554 - t418 * ((t418 * t285 + (t205 + t653) * qJD(1)) * t418 + (t204 * qJD(1) + (-t348 * t551 - t350 * t552 + t556) * t415 + (-t284 + (t597 + t599) * qJD(2) - t452 * qJD(1)) * t418) * t415) + t415 * ((t415 * t284 + (t206 + t656) * qJD(1)) * t415 + (t207 * qJD(1) + (t347 * t551 + t349 * t552) * t418 + (-t285 + (-t596 - t598) * qJD(2) + (t346 - t453) * qJD(1)) * t415) * t418); m(6) * (t100 * t211 + t109 * t188 + t110 * t189 + t210 * t99) + m(5) * (t120 * t281 + t121 * t282 + t174 * t208 + t175 * t209) + (-t217 * t415 - t218 * t418 + (t306 * t415 - t307 * t418) * qJD(1)) * t631 + m(4) * (-t306 * t418 - t307 * t415) * t343 + t419; t422 + m(4) * (-t321 * t343 * t418 + t259 * t132 + t136 * t224 - t320 * t581) + (-t235 * t418 - t236 * t415 + (-t320 * t418 + t321 * t415) * qJD(1)) * t631 + m(6) * (t101 * t38 + t102 * t35 + t103 * t211 + t104 * t210 + t109 * t202 + t110 * t201) + m(5) * (t150 * t67 + t157 * t282 + t158 * t281 + t161 * t59 + t174 * t268 + t175 * t267); (t102 * t38 + t109 * t211 + t110 * t210) * t638 + (t161 * t67 + t174 * t282 + t175 * t281) * t639 + (t343 * t367 * t558 + t136 * t259) * t640 + t422; (-t83 - t95 + (t415 * t520 + t418 * t519) * t408) * t404 + t423 + m(6) * (t100 * t151 + t152 * t99 + t188 * t48 + t189 * t47) + m(5) * (t120 * t197 + t121 * t196 + t208 * t98 + t209 * t97) + t184 + ((t81 / 0.2e1 + t46 / 0.2e1) * t418 + (t45 / 0.2e1 + t82 / 0.2e1) * t415 + (-t415 * t519 + t418 * t520) * qJD(1)) * t402; t420 + m(6) * (t101 * t24 + t103 * t151 + t104 * t152 + t115 * t35 + t201 * t47 + t202 * t48) + m(5) * (t150 * t70 + t157 * t196 + t158 * t197 + t183 * t59 + t267 * t97 + t268 * t98); t420 + m(6) * (t102 * t24 + t109 * t151 + t110 * t152 + t115 * t38 + t210 * t47 + t211 * t48) + m(5) * (t161 * t70 + t174 * t196 + t175 * t197 + t183 * t67 + t281 * t97 + t282 * t98); (t115 * t24 + t151 * t48 + t152 * t47) * t638 + (t183 * t70 + t196 * t98 + t197 * t97) * t639 + (t95 * t404 - t184 - t7 + (-t404 * t463 + t415 * t68 + t418 * t69) * t408) * t404 + (t418 * t10 + t415 * t11 + t463 * t590 + (-t185 * t408 - t45 * t415 - t46 * t418) * t404 + ((-t138 * t404 + t68) * t418 + (t139 * t404 - t56 - t69) * t415) * qJD(1)) * t402 + t496; t423 - t83 * t404 + m(6) * (t100 * t186 + t187 * t99 + t188 * t92 + t189 * t91); t427 + m(6) * (t101 * t54 + t103 * t186 + t104 * t187 + t176 * t35 + t201 * t91 + t202 * t92); t427 + m(6) * (t102 * t54 + t109 * t186 + t110 * t187 + t176 * t38 + t210 * t91 + t211 * t92); m(6) * (t115 * t54 + t151 * t92 + t152 * t91 + t176 * t24 + t186 * t48 + t187 * t47) + t424; (t176 * t54 + t186 * t92 + t187 * t91) * t638 + t424;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq  = res;
