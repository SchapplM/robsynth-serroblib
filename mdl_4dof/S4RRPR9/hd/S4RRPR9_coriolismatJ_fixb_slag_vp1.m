% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRPR9
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRPR9_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR9_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:22
% EndTime: 2019-12-31 17:09:44
% DurationCPUTime: 17.65s
% Computational Cost: add. (31208->656), mult. (45186->942), div. (0->0), fcn. (49206->8), ass. (0->372)
t387 = sin(qJ(2));
t388 = sin(qJ(1));
t514 = t387 * t388;
t389 = cos(qJ(2));
t381 = pkin(7) + qJ(4);
t373 = sin(t381);
t374 = cos(t381);
t556 = rSges(5,1) * t374;
t439 = -rSges(5,2) * t373 + t556;
t290 = -rSges(5,3) * t389 + t387 * t439;
t528 = qJ(3) * t389;
t561 = pkin(2) * t387;
t352 = -t528 + t561;
t386 = -pkin(6) - qJ(3);
t504 = qJ(3) + t386;
t385 = cos(pkin(7));
t372 = pkin(3) * t385 + pkin(2);
t559 = -pkin(2) + t372;
t460 = t387 * t559 + t389 * t504 + t290 + t352;
t195 = t460 * t388;
t390 = cos(qJ(1));
t197 = t460 * t390;
t384 = sin(pkin(7));
t557 = rSges(4,1) * t385;
t442 = -rSges(4,2) * t384 + t557;
t408 = t442 * t387;
t629 = t389 * rSges(4,3) - t408;
t484 = t352 - t629;
t256 = t484 * t388;
t258 = t484 * t390;
t380 = t390 * pkin(5);
t549 = rSges(4,3) + qJ(3);
t560 = pkin(2) * t389;
t608 = t387 * t549 + pkin(1) + t560;
t510 = t388 * t389;
t330 = t384 * t510 + t385 * t390;
t505 = t390 * t384;
t331 = t385 * t510 - t505;
t616 = -t331 * rSges(4,1) + t330 * rSges(4,2);
t200 = -t388 * t608 + t380 + t616;
t332 = t388 * t385 - t389 * t505;
t508 = t389 * t390;
t511 = t388 * t384;
t333 = t385 * t508 + t511;
t443 = t333 * rSges(4,1) + t332 * rSges(4,2);
t201 = t388 * pkin(5) + t390 * t608 + t443;
t497 = t200 * t508 + t201 * t510;
t446 = -pkin(3) * t505 + t372 * t510;
t548 = -rSges(5,3) + t386;
t307 = t373 * t510 + t374 * t390;
t506 = t390 * t373;
t308 = t374 * t510 - t506;
t617 = -t308 * rSges(5,1) + t307 * rSges(5,2);
t171 = t380 + (t387 * t548 - pkin(1)) * t388 - t446 + t617;
t512 = t387 * t390;
t358 = t386 * t512;
t309 = t388 * t374 - t389 * t506;
t310 = t388 * t373 + t374 * t508;
t440 = t310 * rSges(5,1) + t309 * rSges(5,2);
t521 = t372 * t389;
t555 = rSges(5,3) * t387;
t172 = -t358 + (pkin(3) * t384 + pkin(5)) * t388 + (pkin(1) + t521 + t555) * t390 + t440;
t502 = t171 * t508 + t172 * t510;
t599 = m(5) / 0.2e1;
t601 = m(4) / 0.2e1;
t546 = (-t256 * t512 + t258 * t514 + t497) * t601 + (-t195 * t512 + t197 * t514 + t502) * t599;
t212 = (t548 * t389 + (t372 + t439) * t387) * t388;
t474 = t387 * rSges(5,2) * t506 + rSges(5,3) * t508;
t516 = t386 * t389;
t213 = (-t516 + (-t372 - t556) * t387) * t390 + t474;
t370 = pkin(2) * t514;
t254 = t370 + (-t389 * t549 + t408) * t388;
t362 = qJ(3) * t508;
t473 = t387 * rSges(4,2) * t505 + rSges(4,3) * t508;
t255 = t362 + (-pkin(2) - t557) * t512 + t473;
t547 = ((t254 * t390 + t255 * t388) * t387 + t497) * t601 + ((t212 * t390 + t213 * t388) * t387 + t502) * t599;
t630 = t546 - t547;
t647 = t630 * qJD(1);
t429 = Icges(5,5) * t374 - Icges(5,6) * t373;
t283 = -Icges(5,3) * t389 + t387 * t429;
t536 = Icges(5,4) * t374;
t433 = -Icges(5,2) * t373 + t536;
t285 = -Icges(5,6) * t389 + t387 * t433;
t537 = Icges(5,4) * t373;
t435 = Icges(5,1) * t374 - t537;
t287 = -Icges(5,5) * t389 + t387 * t435;
t152 = t283 * t514 - t285 * t307 + t287 * t308;
t214 = Icges(5,5) * t308 - Icges(5,6) * t307 + Icges(5,3) * t514;
t295 = Icges(5,4) * t308;
t217 = -Icges(5,2) * t307 + Icges(5,6) * t514 + t295;
t294 = Icges(5,4) * t307;
t221 = -Icges(5,1) * t308 - Icges(5,5) * t514 + t294;
t116 = t214 * t514 - t217 * t307 - t221 * t308;
t216 = Icges(5,5) * t310 + Icges(5,6) * t309 + Icges(5,3) * t512;
t538 = Icges(5,4) * t310;
t219 = Icges(5,2) * t309 + Icges(5,6) * t512 + t538;
t296 = Icges(5,4) * t309;
t222 = Icges(5,1) * t310 + Icges(5,5) * t512 + t296;
t117 = t216 * t514 - t307 * t219 + t308 * t222;
t428 = t116 * t388 + t117 * t390;
t15 = t152 * t389 - t387 * t428;
t154 = t283 * t512 + t309 * t285 + t310 * t287;
t118 = t214 * t512 + t309 * t217 - t310 * t221;
t119 = t216 * t512 + t309 * t219 + t310 * t222;
t427 = t388 * t118 + t119 * t390;
t645 = -t154 * t389 + t387 * t427;
t240 = Icges(4,5) * t331 - Icges(4,6) * t330 + Icges(4,3) * t514;
t243 = Icges(4,4) * t331 - Icges(4,2) * t330 + Icges(4,6) * t514;
t246 = Icges(4,1) * t331 - Icges(4,4) * t330 + Icges(4,5) * t514;
t627 = t332 * t243 + t333 * t246;
t138 = t240 * t512 + t627;
t242 = Icges(4,5) * t333 + Icges(4,6) * t332 + Icges(4,3) * t512;
t245 = Icges(4,4) * t333 + Icges(4,2) * t332 + Icges(4,6) * t512;
t248 = Icges(4,1) * t333 + Icges(4,4) * t332 + Icges(4,5) * t512;
t642 = -t118 * t390 + t119 * t388;
t644 = (t242 * t512 + t332 * t245 + t333 * t248) * t388 - t138 * t390 + t642;
t62 = -t116 * t390 + t117 * t388;
t525 = t214 * t389;
t641 = t217 * t373 + t221 * t374;
t134 = t641 * t387 + t525;
t223 = rSges(5,3) * t514 - t617;
t174 = t223 * t389 + t290 * t514;
t382 = t388 ^ 2;
t383 = t390 ^ 2;
t470 = t382 + t383;
t576 = -t389 / 0.2e1;
t286 = Icges(5,6) * t387 + t389 * t433;
t288 = Icges(5,5) * t387 + t389 * t435;
t284 = Icges(5,3) * t387 + t389 * t429;
t518 = t374 * t287;
t520 = t373 * t285;
t420 = t518 - t520;
t410 = t284 - t420;
t523 = t283 * t389;
t394 = t387 * t410 + t523;
t112 = -t286 * t307 + t288 * t308 + t388 * t394;
t263 = t285 * t388;
t265 = t287 * t388;
t414 = -t283 * t388 + t641;
t94 = -t414 * t389 + (t263 * t373 - t265 * t374 + t214) * t387;
t637 = t112 + t94;
t113 = t309 * t286 + t310 * t288 + t390 * t394;
t264 = t285 * t390;
t266 = t287 * t390;
t424 = -t219 * t373 + t222 * t374;
t413 = -t283 * t390 - t424;
t95 = -t413 * t389 + (t264 * t373 - t266 * t374 + t216) * t387;
t636 = t113 + t95;
t421 = -t243 * t330 + t246 * t331;
t622 = t290 * t388;
t430 = Icges(4,5) * t385 - Icges(4,6) * t384;
t621 = t387 * t430;
t620 = t387 * t504;
t378 = Icges(3,4) * t389;
t533 = Icges(3,2) * t387;
t540 = Icges(3,1) * t387;
t436 = Icges(4,1) * t385 - Icges(4,4) * t384;
t609 = -Icges(4,5) * t389 + t387 * t436;
t434 = Icges(4,4) * t385 - Icges(4,2) * t384;
t610 = -Icges(4,6) * t389 + t387 * t434;
t619 = (t518 / 0.2e1 - t520 / 0.2e1 - t284 / 0.2e1 + t378 + t540 / 0.2e1 - t533 / 0.2e1 + t385 * t609 / 0.2e1 - t384 * t610 / 0.2e1 - Icges(4,3) * t387 / 0.2e1 + t430 * t576) * t389;
t238 = -rSges(5,1) * t307 - rSges(5,2) * t308;
t239 = rSges(5,1) * t309 - rSges(5,2) * t310;
t158 = (t238 * t390 - t239 * t388) * t387;
t513 = t387 * t389;
t472 = t470 * t513;
t618 = (m(4) / 0.4e1 + m(5) / 0.4e1) * (t472 - t513);
t343 = t470 * t387;
t530 = Icges(4,3) * t389;
t400 = t530 - t621;
t412 = t243 * t384 - t246 * t385 + t400 * t388;
t365 = Icges(3,4) * t514;
t315 = Icges(3,1) * t510 - Icges(3,5) * t390 - t365;
t480 = -Icges(3,2) * t510 + t315 - t365;
t615 = -t412 + t480;
t575 = -t390 / 0.2e1;
t614 = qJD(2) * t575;
t578 = t388 / 0.2e1;
t613 = qJD(2) * t578;
t612 = t242 * t514 - t330 * t245 + t331 * t248;
t314 = Icges(3,6) * t388 + (t378 - t533) * t390;
t437 = -t378 - t540;
t481 = t437 * t390 - t314;
t313 = Icges(3,4) * t510 - Icges(3,2) * t514 - Icges(3,6) * t390;
t482 = -t437 * t388 + t313;
t607 = (t481 * t388 + t390 * t482) * t389;
t318 = (-Icges(5,2) * t374 - t537) * t387;
t319 = (-Icges(5,1) * t373 - t536) * t387;
t606 = -t373 * (t287 / 0.2e1 + t318 / 0.2e1) + t374 * (t319 / 0.2e1 - t285 / 0.2e1);
t604 = 0.4e1 * qJD(1);
t603 = 0.2e1 * qJD(2);
t597 = -t645 / 0.2e1;
t596 = t62 / 0.2e1;
t595 = t642 / 0.2e1;
t232 = -Icges(5,5) * t307 - Icges(5,6) * t308;
t493 = -Icges(5,2) * t308 - t221 - t294;
t495 = -Icges(5,1) * t307 - t217 - t295;
t99 = -t232 * t389 + (-t373 * t493 + t374 * t495) * t387;
t594 = t99 / 0.2e1;
t268 = -t512 * t556 + t474;
t225 = rSges(5,3) * t512 + t440;
t423 = t223 * t390 - t225 * t388;
t125 = t423 * t389 + (-t268 * t388 - t390 * t622) * t387;
t291 = t389 * t439 + t555;
t144 = (t291 * t388 - t223) * t387;
t145 = (-t290 * t390 - t268) * t389 + (-t291 * t390 + t225) * t387;
t156 = t423 * t387;
t176 = t389 * t225 + t290 * t512;
t501 = t174 * t508 - t176 * t510;
t591 = m(5) * (-t125 * t389 + (t144 * t390 + t145 * t388 + t156) * t387 + t501);
t590 = m(5) * (t125 * t156 + t144 * t174 - t145 * t176);
t588 = m(5) * (t144 * t171 + t145 * t172 + t174 * t212 - t176 * t213);
t320 = (-rSges(5,1) * t373 - rSges(5,2) * t374) * t387;
t586 = m(5) * (-t195 * t239 + t197 * t238 + (-t171 * t390 - t172 * t388) * t320);
t355 = qJ(3) * t387 + t560;
t477 = t470 * t355;
t104 = (t223 - (t560 + t620) * t388 + t446) * t388 + (pkin(3) * t511 + t225 - t358 + (-t355 + t521) * t390) * t390 + t477;
t498 = -t195 * t510 - t197 * t508;
t585 = m(5) * (t104 * t343 + t498);
t582 = m(5) * (t156 * t343 + t501);
t580 = m(5) * (t171 * t212 + t172 * t213);
t577 = t388 / 0.4e1;
t574 = -t390 / 0.4e1;
t573 = t390 / 0.2e1;
t558 = rSges(3,1) * t389;
t451 = pkin(1) + t558;
t471 = rSges(3,2) * t514 + t390 * rSges(3,3);
t276 = -t388 * t451 + t380 + t471;
t367 = rSges(3,2) * t512;
t277 = -t367 + t451 * t390 + (rSges(3,3) + pkin(5)) * t388;
t353 = rSges(3,1) * t387 + rSges(3,2) * t389;
t340 = t353 * t388;
t342 = t353 * t390;
t572 = m(3) * (t276 * t340 - t277 * t342);
t148 = t388 * (rSges(4,3) * t514 - t616) + t390 * (rSges(4,3) * t512 + t443) + t477;
t490 = -t256 * t510 - t258 * t508;
t569 = m(4) * (t148 * t343 + t490);
t567 = m(4) * (t200 * t254 + t201 * t255);
t194 = t201 * t512;
t566 = m(4) * (-t200 * t514 + t194);
t165 = t172 * t512;
t565 = m(5) * (-t171 * t514 + t165);
t564 = m(5) * (-t174 * t514 - t176 * t512);
t162 = t388 * t238 + t239 * t390;
t563 = m(5) * (-t162 * t389 - t320 * t343);
t562 = m(5) * t158;
t553 = t388 * t15;
t550 = t390 * t645;
t539 = Icges(3,4) * t387;
t524 = t216 * t389;
t519 = t373 * t286;
t517 = t374 * t288;
t515 = t387 * t313;
t317 = (-Icges(5,5) * t373 - Icges(5,6) * t374) * t387;
t509 = t389 * t317;
t494 = Icges(5,1) * t309 - t219 - t538;
t492 = -Icges(5,2) * t310 + t222 + t296;
t311 = Icges(3,5) * t510 - Icges(3,6) * t514 - Icges(3,3) * t390;
t489 = -t388 * t311 - t315 * t508;
t432 = Icges(3,5) * t389 - Icges(3,6) * t387;
t312 = Icges(3,3) * t388 + t390 * t432;
t351 = Icges(3,1) * t389 - t539;
t316 = Icges(3,5) * t388 + t351 * t390;
t488 = t388 * t312 + t316 * t508;
t487 = -t285 + t319;
t486 = t287 + t318;
t483 = -rSges(4,3) * t387 - t389 * t442 - t355;
t348 = Icges(3,2) * t389 + t539;
t479 = -t348 * t390 + t316;
t478 = t388 * (qJ(3) * t510 - t370) + t390 * (-pkin(2) * t512 + t362);
t469 = qJD(1) * t387;
t468 = qJD(4) * t387;
t396 = t387 * t414 + t525;
t80 = t263 * t307 - t265 * t308 + t388 * t396;
t395 = t387 * t413 + t524;
t81 = t264 * t307 - t266 * t308 + t388 * t395;
t11 = (-t112 + t428) * t389 + (t388 * t80 + t390 * t81 + t152) * t387;
t84 = t232 * t514 - t307 * t493 + t308 * t495;
t233 = Icges(5,5) * t309 - Icges(5,6) * t310;
t85 = t233 * t514 - t307 * t492 + t308 * t494;
t43 = t85 * t388 - t390 * t84;
t463 = t43 / 0.2e1 - t11 / 0.2e1;
t82 = -t309 * t263 - t310 * t265 + t390 * t396;
t83 = -t309 * t264 - t310 * t266 + t390 * t395;
t12 = (-t113 + t427) * t389 + (t388 * t82 + t390 * t83 + t154) * t387;
t86 = t232 * t512 + t309 * t493 + t310 * t495;
t87 = t233 * t512 + t309 * t492 + t310 * t494;
t44 = t87 * t388 - t390 * t86;
t462 = t44 / 0.2e1 - t12 / 0.2e1;
t461 = t597 + t645 / 0.2e1;
t459 = -t389 * t559 - t291 - t355 + t620;
t457 = t514 / 0.4e1;
t449 = t479 * t388;
t278 = t316 * t510;
t448 = t312 * t390 - t278;
t447 = t387 * t314 - t311;
t100 = -t233 * t389 + (-t373 * t492 + t374 * t494) * t387;
t122 = -t307 * t486 + t308 * t487 + t317 * t514;
t123 = t309 * t486 + t310 * t487 + t317 * t512;
t445 = t586 / 0.2e1 + (t123 + t100) * t577 + (t122 + t99) * t574;
t431 = -Icges(3,5) * t387 - Icges(3,6) * t389;
t135 = t387 * t424 - t524;
t426 = -t134 * t388 + t135 * t390;
t418 = -t372 * t387 - t516;
t417 = m(5) * (-t171 * t238 + t172 * t239) - t509 / 0.2e1;
t411 = t245 * t384 - t248 * t385 + t400 * t390;
t407 = (-t240 * t390 + t242 * t388) * t389;
t406 = t411 * t388;
t405 = t15 * t577 + t645 * t574 - t553 / 0.4e1 + t550 / 0.4e1 + (t457 - t514 / 0.4e1) * t642;
t128 = -t410 * t389 + (t283 + t517 - t519) * t387;
t166 = t387 * t420 - t523;
t397 = t128 * t576 + t166 * t387 / 0.2e1 + t588 / 0.2e1 + t637 * t457 + t636 * t512 / 0.4e1 + (-t134 + t152) * t510 / 0.4e1 + (t135 + t154) * t508 / 0.4e1;
t302 = Icges(4,6) * t387 + t389 * t434;
t304 = Icges(4,5) * t387 + t389 * t436;
t391 = -t385 * t304 / 0.2e1 + t384 * t302 / 0.2e1 - t517 / 0.2e1 + t519 / 0.2e1 - t283 / 0.2e1 - t351 / 0.2e1 + t348 / 0.2e1 - t621 / 0.2e1 + t530 / 0.2e1;
t356 = -rSges(3,2) * t387 + t558;
t335 = t431 * t390;
t334 = t431 * t388;
t275 = t609 * t390;
t274 = t609 * t388;
t273 = t610 * t390;
t272 = t610 * t388;
t259 = t483 * t390;
t257 = t483 * t388;
t211 = 0.4e1 * t618;
t198 = t459 * t390;
t196 = t459 * t388;
t185 = -t389 * t239 - t320 * t512;
t184 = t238 * t389 + t320 * t514;
t180 = -t314 * t512 + t488;
t179 = -t313 * t512 - t489;
t178 = -t314 * t514 - t448;
t160 = t390 * (-t512 * t557 + t473) + t629 * t382 + t478;
t157 = -t562 / 0.2e1;
t146 = -t509 + (-t373 * t486 + t374 * t487) * t387;
t142 = t563 / 0.2e1;
t141 = -t179 * t390 + t180 * t388;
t140 = -(-t388 * (-t389 * t315 + t515) - t311 * t390) * t390 + t178 * t388;
t129 = (-t362 + t268 + (t418 + t561) * t390) * t390 + (t370 - t622 + (t418 - t528) * t388) * t388 + t478;
t127 = t564 / 0.2e1;
t88 = t582 / 0.2e1;
t78 = t565 + t566;
t75 = -(t240 * t514 + t421) * t390 + t388 * t612;
t61 = t387 * t606 + t417;
t56 = -t166 * t389 + t387 * t426;
t55 = (t178 - t278 + (t312 + t515) * t390 + t489) * t390 + t488 * t388;
t54 = (t390 * t447 + t180 - t488) * t390 + (t388 * t447 + t179 + t448) * t388;
t51 = t104 * t162 + (t195 * t388 + t197 * t390) * t320;
t45 = t569 + t585;
t42 = t83 * t388 - t390 * t82;
t41 = t81 * t388 - t390 * t80;
t36 = t127 + t562 / 0.2e1;
t35 = t157 + t127;
t34 = t157 - t564 / 0.2e1;
t32 = t591 / 0.2e1;
t31 = -t123 * t389 + (t388 * t86 + t390 * t87) * t387;
t30 = -t122 * t389 + (t388 * t84 + t390 * t85) * t387;
t29 = -t387 * t391 + t567 + t572 + t580 + t619;
t27 = t421 * t390 + (t138 - t612 - t627) * t388;
t22 = (-t128 + t426) * t389 + (t94 * t388 + t95 * t390 + t166) * t387;
t21 = t88 + t32 - t563 / 0.2e1;
t20 = t142 + t88 - t591 / 0.2e1;
t19 = t142 + t32 - t582 / 0.2e1;
t9 = t546 + t547;
t7 = m(5) * t51 + t43 * t575 + t44 * t578;
t6 = t461 * t514;
t5 = t590 + (t550 / 0.2e1 - t553 / 0.2e1 - t22 / 0.2e1) * t389 + (t12 * t573 + t11 * t578 + t56 / 0.2e1) * t387;
t4 = (t595 - t642 / 0.2e1 - t55 / 0.2e1 + t141 / 0.2e1) * t390 + (-t62 / 0.2e1 + t596 + t140 / 0.2e1 + t27 / 0.2e1 + t75 / 0.2e1 + t54 / 0.2e1) * t388;
t3 = t397 + t445;
t2 = -t588 / 0.2e1 + t405 + (-t166 / 0.2e1 + (-t113 / 0.4e1 - t95 / 0.4e1) * t390 + (-t112 / 0.4e1 - t94 / 0.4e1) * t388) * t387 + (t128 / 0.2e1 + (-t154 / 0.4e1 - t135 / 0.4e1) * t390 + (-t152 / 0.4e1 + t134 / 0.4e1) * t388) * t389 + t445;
t1 = -t586 / 0.2e1 + (-t123 / 0.4e1 - t100 / 0.4e1) * t388 + t397 + t405 + (t122 / 0.4e1 + t99 / 0.4e1) * t390;
t8 = [t29 * qJD(2) + t78 * qJD(3) + t61 * qJD(4), t29 * qJD(1) + t9 * qJD(3) + t3 * qJD(4) + (m(3) * ((-t276 * t390 - t277 * t388) * t356 + (-t340 * t390 + t342 * t388) * t353) / 0.2e1 + (t171 * t198 + t172 * t196 - t195 * t213 - t197 * t212) * t599 + (t200 * t259 + t201 * t257 - t254 * t258 - t255 * t256) * t601) * t603 + (t332 * t302 + t333 * t304 + (t479 - t411) * t389 + (t273 * t384 - t275 * t385 + t242 + t481) * t387 + t636) * t613 + (-t330 * t302 + t331 * t304 + t141 + t615 * t389 + (t272 * t384 - t274 * t385 + t240 - t482) * t387 + t637 + t644) * t614 + ((t382 / 0.2e1 + t383 / 0.2e1) * t432 + (t55 + t644) * t573 - (t140 + t27 + t75 + t54) * t388 / 0.2e1) * qJD(2), qJD(1) * t78 + qJD(2) * t9 + qJD(4) * t35, t61 * qJD(1) + t3 * qJD(2) + t35 * qJD(3) + (m(5) * (t171 * t184 + t172 * t185 - t174 * t238 - t176 * t239) - t146 * t389) * qJD(4) + ((t100 / 0.2e1 + t123 / 0.2e1) * t390 + (t594 + t122 / 0.2e1 - t461) * t388) * t468; t4 * qJD(2) + t630 * qJD(3) + t2 * qJD(4) + (-t580 / 0.4e1 - t567 / 0.4e1 - t572 / 0.4e1) * t604 + t391 * t469 - t619 * qJD(1), t4 * qJD(1) + (m(5) * (t104 * t129 - t195 * t196 - t197 * t198) + m(4) * (t148 * t160 - t256 * t257 - t258 * t259) + m(3) * ((t388 * (rSges(3,1) * t510 - t471) + t390 * (rSges(3,1) * t508 + t388 * rSges(3,3) - t367)) * (-t388 * t340 - t342 * t390) + t470 * t356 * t353)) * qJD(2) + t45 * qJD(3) + t7 * qJD(4) + (t42 + (-t332 * t273 - t333 * t275) * t388 + (t332 * t272 + t333 * t274 + t407 + (-t390 * t412 + t406) * t387) * t390 + t382 * t335 + (-t388 * t334 + t607 + (t390 * t480 - t449) * t387) * t390) * t613 + (t41 + t383 * t334 - (t272 * t330 - t274 * t331) * t390 + (-t390 * t335 + t607 + t330 * t273 - t331 * t275 + t407 + (t390 * t615 + t406 - t449) * t387) * t388) * t614, t647 + t45 * qJD(2) + t20 * qJD(4) + (-0.4e1 * t618 + 0.2e1 * (t599 + t601) * (-t343 * t389 + t472)) * qJD(3), t2 * qJD(1) + t7 * qJD(2) + t20 * qJD(3) + (-t56 / 0.2e1 + t462 * t390 + t463 * t388) * t468 + (m(5) * (t158 * t104 + t156 * t162 - t184 * t197 - t185 * t195 + (-t174 * t390 + t176 * t388) * t320) + t31 * t578 + t30 * t575 - t590 + (t22 / 0.2e1 + (t594 + t597) * t390 + (-t100 / 0.2e1 + t15 / 0.2e1) * t388) * t389) * qJD(4); -t630 * qJD(2) + t34 * qJD(4) + (-t565 / 0.4e1 - t566 / 0.4e1) * t604 + 0.2e1 * (t165 * t599 + t194 * t601 + (-t172 * t599 - t201 * t601) * t512) * qJD(1), -t647 + t211 * qJD(3) + t19 * qJD(4) + 0.4e1 * (-t585 / 0.4e1 - t569 / 0.4e1) * qJD(2) + ((-t389 * t129 + t498) * t599 + (-t389 * t160 + t490) * t601 + ((t196 * t388 + t198 * t390 + t104) * t599 + (t257 * t388 + t259 * t390 + t148) * t601) * t387) * t603, t211 * qJD(2), t34 * qJD(1) + t19 * qJD(2) + m(5) * (-t158 * t389 + (t184 * t390 + t185 * t388) * t387) * qJD(4); -t417 * qJD(1) + t1 * qJD(2) + t36 * qJD(3) + t6 * qJD(4) - t606 * t469, t1 * qJD(1) + (((t595 + t94 / 0.2e1) * t389 + t463) * t390 + ((t596 - t95 / 0.2e1) * t389 - t462) * t388 + ((t42 / 0.2e1 + t134 / 0.2e1) * t390 + (t41 / 0.2e1 + t135 / 0.2e1) * t388) * t387 + (t104 * t125 + t129 * t156 - t144 * t197 - t145 * t195 + t174 * t198 - t176 * t196 - t51) * m(5)) * qJD(2) + t21 * qJD(3) + t5 * qJD(4), qJD(1) * t36 + qJD(2) * t21, t6 * qJD(1) + t5 * qJD(2) + (m(5) * (t156 * t158 + t174 * t184 - t176 * t185) + t389 ^ 2 * t146 / 0.2e1 + (t31 * t573 + t30 * t578 + (t100 * t390 + t99 * t388) * t576) * t387) * qJD(4);];
Cq = t8;
