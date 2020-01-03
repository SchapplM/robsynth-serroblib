% Calculate matrix of centrifugal and coriolis load on the joints for
% S4RRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S4RRRR3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR3_coriolismatJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:24:17
% EndTime: 2019-12-31 17:24:30
% DurationCPUTime: 8.09s
% Computational Cost: add. (36438->483), mult. (34384->630), div. (0->0), fcn. (32184->8), ass. (0->299)
t368 = qJ(2) + qJ(3);
t360 = qJ(4) + t368;
t345 = cos(t360);
t370 = sin(qJ(1));
t485 = t345 * t370;
t344 = sin(t360);
t489 = t344 * t370;
t263 = rSges(5,1) * t489 + rSges(5,2) * t485;
t293 = rSges(5,1) * t344 + rSges(5,2) * t345;
t372 = cos(qJ(1));
t264 = t293 * t372;
t163 = -t370 * t263 - t372 * t264;
t366 = t370 ^ 2;
t367 = t372 ^ 2;
t434 = t366 + t367;
t348 = sin(t368);
t524 = pkin(3) * t348;
t148 = -t434 * t524 + t163;
t369 = sin(qJ(2));
t525 = pkin(2) * t369;
t313 = -t524 - t525;
t389 = t293 - t313;
t177 = t389 * t370;
t179 = t389 * t372;
t516 = rSges(5,1) * t345;
t294 = -rSges(5,2) * t344 + t516;
t349 = cos(t368);
t523 = pkin(3) * t349;
t426 = -t294 - t523;
t202 = t426 * t370;
t204 = t426 * t372;
t459 = -t177 * t202 - t179 * t204;
t227 = rSges(5,1) * t485 - rSges(5,2) * t489 - t372 * rSges(5,3);
t488 = t344 * t372;
t421 = -rSges(5,2) * t488 + rSges(5,3) * t370;
t484 = t345 * t372;
t157 = t370 * t227 + t372 * (rSges(5,1) * t484 + t421);
t371 = cos(qJ(2));
t364 = t371 * pkin(2);
t346 = t364 + pkin(1);
t312 = t346 + t523;
t373 = -pkin(6) - pkin(5);
t433 = -pkin(7) + t373;
t341 = t370 * t433;
t347 = t370 * t373;
t437 = -t370 * t312 - t372 * t433;
t461 = t372 * t373;
t108 = -t370 * (t370 * t346 + t437 + t461) + t157 + t372 * (-t341 + t347 + (t312 - t346) * t372);
t365 = t372 * pkin(5);
t522 = pkin(1) - t346;
t452 = -t370 * (t370 * t522 - t365 - t461) + t372 * (-pkin(5) * t370 - t372 * t522 - t347);
t89 = t108 + t452;
t54 = t89 * t148 + t459;
t480 = t348 * t372;
t422 = -rSges(4,2) * t480 + rSges(4,3) * t370;
t481 = t348 * t370;
t436 = rSges(4,2) * t481 + t372 * rSges(4,3);
t476 = t349 * t372;
t477 = t349 * t370;
t161 = t370 * (rSges(4,1) * t477 - t436) + t372 * (rSges(4,1) * t476 + t422);
t124 = t161 + t452;
t310 = rSges(4,1) * t348 + rSges(4,2) * t349;
t286 = t310 * t370;
t287 = t310 * t372;
t169 = -t370 * t286 - t372 * t287;
t517 = rSges(4,1) * t349;
t311 = -rSges(4,2) * t348 + t517;
t383 = t310 + t525;
t570 = t383 * t372;
t571 = t383 * t370;
t76 = t124 * t169 + (t370 * t571 + t372 * t570) * t311;
t582 = -m(4) * t76 - m(5) * t54;
t335 = Icges(5,4) * t345;
t290 = -Icges(5,2) * t344 + t335;
t291 = Icges(5,1) * t344 + t335;
t581 = t290 + t291;
t343 = Icges(4,4) * t349;
t307 = -Icges(4,2) * t348 + t343;
t308 = Icges(4,1) * t348 + t343;
t580 = t308 + t307;
t558 = m(4) / 0.2e1;
t557 = m(5) / 0.2e1;
t535 = t370 / 0.2e1;
t533 = -t372 / 0.2e1;
t577 = t372 / 0.2e1;
t504 = Icges(3,4) * t369;
t322 = Icges(3,2) * t371 + t504;
t325 = Icges(3,1) * t371 - t504;
t576 = (t325 / 0.2e1 - t322 / 0.2e1) * t369;
t181 = -t313 * t370 + t263;
t304 = t372 * t313;
t182 = t304 - t264;
t427 = t293 + t524;
t201 = t427 * t370;
t423 = t346 + t517;
t174 = -t370 * t423 + t436 - t461;
t175 = t372 * t423 - t347 + t422;
t386 = (-t174 * t372 - t175 * t370) * t311;
t403 = t427 * t372;
t164 = -t227 + t437;
t165 = -t341 + (t312 + t516) * t372 + t421;
t460 = t204 * t164 + t202 * t165;
t519 = (-t181 * t403 - t182 * t201 + t460) * t557 + (t386 + (t370 * t570 - t372 * t571) * t310) * t558;
t193 = pkin(3) * t481 + t263;
t520 = (t177 * t403 - t179 * t193 + t460) * t557 + (-t286 * t570 + t287 * t571 + t386) * t558;
t575 = t519 - t520;
t536 = -t370 / 0.2e1;
t574 = t535 + t536;
t359 = Icges(3,4) * t371;
t323 = -Icges(3,2) * t369 + t359;
t324 = Icges(3,1) * t369 + t359;
t503 = Icges(4,4) * t348;
t306 = Icges(4,2) * t349 + t503;
t309 = Icges(4,1) * t349 - t503;
t567 = t580 * t349 / 0.2e1 + (t309 / 0.2e1 - t306 / 0.2e1) * t348;
t502 = Icges(5,4) * t344;
t289 = Icges(5,2) * t345 + t502;
t292 = Icges(5,1) * t345 - t502;
t408 = t581 * t345 / 0.2e1 + (-t289 / 0.2e1 + t292 / 0.2e1) * t344;
t249 = Icges(4,6) * t370 + t307 * t372;
t251 = Icges(4,5) * t370 + t309 * t372;
t188 = t251 * t477;
t305 = Icges(4,5) * t349 - Icges(4,6) * t348;
t493 = t305 * t372;
t247 = Icges(4,3) * t370 + t493;
t413 = t372 * t247 - t188;
t133 = -t249 * t481 - t413;
t248 = Icges(4,4) * t477 - Icges(4,2) * t481 - Icges(4,6) * t372;
t246 = Icges(4,5) * t477 - Icges(4,6) * t481 - Icges(4,3) * t372;
t332 = Icges(4,4) * t481;
t250 = Icges(4,1) * t477 - Icges(4,5) * t372 - t332;
t455 = -t370 * t246 - t250 * t476;
t134 = -t248 * t480 - t455;
t454 = t370 * t247 + t251 * t476;
t135 = -t249 * t480 + t454;
t410 = t249 * t348 - t246;
t496 = t248 * t348;
t566 = (-t134 * t372 + t135 * t370) * t577 + ((t133 - t188 + (t247 + t496) * t372 + t455) * t372 + t454 * t370) * t533 + ((t370 * t410 + t133 + t134 + t413) * t370 + (t370 * (-t250 * t349 + t496) + t135 - t454 + (t246 + t410) * t372) * t372) * t535;
t220 = Icges(5,6) * t370 + t290 * t372;
t222 = Icges(5,5) * t370 + t292 * t372;
t183 = t222 * t485;
t288 = Icges(5,5) * t345 - Icges(5,6) * t344;
t495 = t288 * t372;
t218 = Icges(5,3) * t370 + t495;
t414 = t372 * t218 - t183;
t127 = -t220 * t489 - t414;
t219 = Icges(5,4) * t485 - Icges(5,2) * t489 - Icges(5,6) * t372;
t217 = Icges(5,5) * t485 - Icges(5,6) * t489 - Icges(5,3) * t372;
t316 = Icges(5,4) * t489;
t221 = Icges(5,1) * t485 - Icges(5,5) * t372 - t316;
t457 = -t370 * t217 - t221 * t484;
t128 = -t219 * t488 - t457;
t456 = t370 * t218 + t222 * t484;
t129 = -t220 * t488 + t456;
t412 = t220 * t344 - t217;
t497 = t219 * t344;
t425 = ((t127 - t183 + (t218 + t497) * t372 + t457) * t372 + t456 * t370) * t533 + (-t128 * t372 + t129 * t370) * t577 + ((t370 * t412 + t127 + t128 + t414) * t370 + (t129 - t456 + t370 * (-t221 * t345 + t497) + (t412 + t217) * t372) * t372) * t535;
t278 = Icges(3,5) * t370 + t325 * t372;
t438 = -t322 * t372 + t278;
t474 = t369 * t370;
t338 = Icges(3,4) * t474;
t470 = t370 * t371;
t277 = Icges(3,1) * t470 - Icges(3,5) * t372 - t338;
t439 = -Icges(3,2) * t470 + t277 - t338;
t276 = Icges(3,6) * t370 + t323 * t372;
t440 = -t324 * t372 - t276;
t275 = Icges(3,4) * t470 - Icges(3,2) * t474 - Icges(3,6) * t372;
t441 = t324 * t370 + t275;
t565 = (-t438 * t370 + t372 * t439) * t369 + (t440 * t370 + t372 * t441) * t371;
t442 = -t306 * t372 + t251;
t443 = -Icges(4,2) * t477 + t250 - t332;
t444 = -t308 * t372 - t249;
t445 = t308 * t370 + t248;
t564 = (-t442 * t370 + t372 * t443) * t348 + (t444 * t370 + t372 * t445) * t349;
t448 = -t289 * t372 + t222;
t449 = -Icges(5,2) * t485 + t221 - t316;
t450 = -t291 * t372 - t220;
t451 = t291 * t370 + t219;
t563 = (-t448 * t370 + t372 * t449) * t344 + (t450 * t370 + t372 * t451) * t345;
t562 = 4 * qJD(1);
t561 = 2 * qJD(2);
t559 = 2 * qJD(3);
t71 = t89 * t163;
t96 = t108 * t163;
t550 = m(5) * (t71 + t96 + ((t179 + t403) * t372 + (t177 + t201) * t370) * t294);
t116 = t157 * t148;
t385 = (-t202 * t370 - t204 * t372) * t293;
t56 = t71 + (t177 * t370 + t179 * t372) * t294;
t548 = m(5) * (t116 + t385 + t56);
t473 = t369 * t372;
t130 = t366 * (t313 + t525) + t372 * (pkin(2) * t473 + t304) + t163;
t64 = t96 + (t201 * t370 + t372 * t403) * t294;
t375 = t385 + t64;
t547 = m(5) * (t130 * t157 + t375);
t458 = -t201 * t202 - t204 * t403;
t543 = m(5) * (t108 * t130 + t458);
t387 = (-t164 * t372 - t165 * t370) * t294;
t540 = m(5) * (t177 * t264 - t179 * t263 + t387);
t539 = m(5) * (t387 + (-t181 * t372 - t182 * t370) * t293);
t538 = m(5) * (t201 * t264 - t263 * t403 + t387);
t537 = m(5) * (t387 + (-t193 * t372 + t370 * t403) * t293);
t518 = rSges(3,1) * t371;
t429 = pkin(1) + t518;
t435 = rSges(3,2) * t474 + t372 * rSges(3,3);
t205 = -t370 * t429 + t365 + t435;
t340 = rSges(3,2) * t473;
t206 = -t340 + t429 * t372 + (rSges(3,3) + pkin(5)) * t370;
t326 = rSges(3,1) * t369 + rSges(3,2) * t371;
t302 = t326 * t370;
t303 = t326 * t372;
t532 = m(3) * (t205 * t302 - t206 * t303);
t103 = t434 * t310 * t311 + t161 * t169;
t102 = m(4) * t103;
t531 = m(4) * (t174 * t571 - t175 * t570);
t530 = m(4) * (t174 * t286 - t175 * t287);
t100 = t434 * t293 * t294 + t157 * t163;
t529 = m(5) * t100;
t528 = m(5) * (t164 * t181 + t165 * t182);
t527 = m(5) * (t164 * t193 - t165 * t403);
t526 = m(5) * (t164 * t263 - t165 * t264);
t396 = Icges(5,5) * t344 + Icges(5,6) * t345;
t257 = t396 * t370;
t258 = t372 * t396;
t521 = (-t366 * t258 + (t370 * t257 + t563) * t372) * t535 + (-t367 * t257 + (t372 * t258 + t563) * t370) * t533;
t475 = t369 * t275;
t469 = t371 * t372;
t273 = Icges(3,5) * t470 - Icges(3,6) * t474 - Icges(3,3) * t372;
t447 = -t370 * t273 - t277 * t469;
t399 = Icges(3,5) * t371 - Icges(3,6) * t369;
t274 = Icges(3,3) * t370 + t372 * t399;
t446 = t370 * t274 + t278 * t469;
t432 = t547 / 0.2e1 + t521;
t431 = -t529 + t521;
t430 = t108 * t148 + t458;
t428 = -t311 - t364;
t397 = Icges(4,5) * t348 + Icges(4,6) * t349;
t280 = t397 * t370;
t281 = t372 * t397;
t424 = (-t366 * t281 + (t370 * t280 + t564) * t372) * t535 + (-t367 * t280 + (t372 * t281 + t564) * t370) * t533 + t521;
t223 = t278 * t470;
t411 = t372 * t274 - t223;
t409 = t369 * t276 - t273;
t407 = t434 * t525;
t404 = t102 + t424;
t398 = -Icges(3,5) * t369 - Icges(3,6) * t371;
t388 = t426 - t364;
t384 = t425 + t566;
t379 = (-t289 + t292) * t345 - t581 * t344;
t382 = -t425 + (t288 * t370 + t344 * t450 + t345 * t448 + t372 * t379) * t535 + (-t344 * t451 + t345 * t449 + t370 * t379 - t495) * t533;
t381 = -t408 + t574 * (t219 * t345 + t221 * t344);
t380 = t408 + t567;
t378 = (-t306 + t309) * t349 - t580 * t348;
t376 = t382 - t566 + (t305 * t370 + t348 * t444 + t349 * t442 + t372 * t378) * t535 + (-t348 * t445 + t349 * t443 + t370 * t378 - t493) * t533;
t374 = t381 - t567 + t574 * (t248 * t349 + t250 * t348);
t328 = -rSges(3,2) * t369 + t518;
t297 = t398 * t372;
t296 = t398 * t370;
t236 = t428 * t372;
t234 = t428 * t370;
t180 = t388 * t372;
t178 = t388 * t370;
t156 = -t407 + t169;
t144 = -t276 * t473 + t446;
t143 = -t275 * t473 - t447;
t142 = -t276 * t474 - t411;
t118 = -t407 + t130;
t105 = -t143 * t372 + t144 * t370;
t104 = -(-t370 * (-t371 * t277 + t475) - t372 * t273) * t372 + t142 * t370;
t75 = t537 / 0.2e1;
t74 = t408 + t526;
t73 = t538 / 0.2e1;
t70 = t539 / 0.2e1;
t68 = t540 / 0.2e1;
t50 = t380 + t527 + t530;
t49 = (t142 - t223 + (t274 + t475) * t372 + t447) * t372 + t446 * t370;
t48 = (t372 * t409 + t144 - t446) * t372 + (t370 * t409 + t143 + t411) * t370;
t35 = t548 / 0.2e1;
t26 = (t324 / 0.2e1 + t323 / 0.2e1) * t371 + t576 + t532 + t531 + t528 + t380;
t24 = t550 / 0.2e1;
t23 = t521 + t529;
t22 = t23 * qJD(4);
t19 = m(5) * t64 + t521;
t18 = m(5) * t56 + t521;
t17 = t404 + t543;
t15 = t424 - t582;
t14 = t24 - t548 / 0.2e1 + t432;
t13 = t24 + t35 - t547 / 0.2e1 + t521;
t12 = t35 - t550 / 0.2e1 + t432;
t11 = t73 - t537 / 0.2e1 + t425;
t10 = t75 - t538 / 0.2e1 + t425;
t9 = t68 - t539 / 0.2e1 + t425;
t8 = t70 - t540 / 0.2e1 + t425;
t7 = t73 + t75 + t382;
t6 = t68 + t70 + t382;
t4 = (t105 / 0.2e1 - t49 / 0.2e1) * t372 + (t48 / 0.2e1 + t104 / 0.2e1) * t370 + t384;
t3 = t384 + t575;
t2 = t384 - t575;
t1 = t376 + t519 + t520;
t5 = [t26 * qJD(2) + t50 * qJD(3) + qJD(4) * t74, t26 * qJD(1) + t1 * qJD(3) + t6 * qJD(4) + (m(3) * ((-t205 * t372 - t206 * t370) * t328 + (-t302 * t372 + t303 * t370) * t326) / 0.2e1 + (t164 * t180 + t165 * t178 - t177 * t182 - t179 * t181) * t557 + (t174 * t236 + t175 * t234) * t558) * t561 + (t376 + (t369 * t440 + t371 * t438) * t535 + t49 * t577 + (t48 + t104) * t536 + (-t369 * t441 + t371 * t439 + t105) * t533 + (t366 / 0.2e1 + t367 / 0.2e1) * t399) * qJD(2), t50 * qJD(1) + t1 * qJD(2) + t376 * qJD(3) + t7 * qJD(4) + ((t386 + (-t286 * t372 + t287 * t370) * t310) * t558 + (t460 + (-t193 + t201) * t403) * t557) * t559, t74 * qJD(1) + t6 * qJD(2) + t7 * qJD(3) + ((t387 + (-t263 * t372 + t264 * t370) * t293) * m(5) + t382) * qJD(4); (t374 - (t324 + t323) * t371 / 0.2e1 - t576) * qJD(1) + t4 * qJD(2) + t2 * qJD(3) + t9 * qJD(4) + (-t528 / 0.4e1 - t531 / 0.4e1 - t532 / 0.4e1) * t562, t4 * qJD(1) + (m(5) * (t118 * t89 - t177 * t178 - t179 * t180) + m(4) * (t124 * t156 - t234 * t571 - t236 * t570) + (t366 * t297 + (-t370 * t296 + t565) * t372) * t535 + m(3) * ((t370 * (rSges(3,1) * t470 - t435) + t372 * (rSges(3,1) * t469 + rSges(3,3) * t370 - t340)) * (-t302 * t370 - t303 * t372) + t434 * t328 * t326) + (t367 * t296 + (-t372 * t297 + t565) * t370) * t533 + t424) * qJD(2) + t15 * qJD(3) + t18 * qJD(4), t2 * qJD(1) + t15 * qJD(2) + t13 * qJD(4) + ((t430 + t54) * t557 + (t76 + t103) * t558) * t559 + (t424 - t102 - t543) * qJD(3), t9 * qJD(1) + t18 * qJD(2) + t13 * qJD(3) + (m(5) * (t56 + t100) + t431) * qJD(4); t374 * qJD(1) + t3 * qJD(2) + t384 * qJD(3) + t11 * qJD(4) + (-t530 / 0.4e1 - t527 / 0.4e1) * t562, t3 * qJD(1) + t17 * qJD(3) + t14 * qJD(4) + ((t108 * t118 + t130 * t89 - t178 * t201 - t180 * t403 + t459) * t557 + (t156 * t161 + (-t234 * t370 - t236 * t372) * t310 + t76) * t558) * t561 + (t424 + t582) * qJD(2), t384 * qJD(1) + t17 * qJD(2) + (m(5) * t430 + t404) * qJD(3) + t19 * qJD(4), t11 * qJD(1) + t14 * qJD(2) + t19 * qJD(3) + (m(5) * (t64 + t100) + t431) * qJD(4); (t381 - t526) * qJD(1) + t8 * qJD(2) + t10 * qJD(3) + t425 * qJD(4), t8 * qJD(1) + ((t118 * t157 + (-t178 * t370 - t180 * t372) * t293) * m(5) + t521) * qJD(2) + t12 * qJD(3) + t22, t10 * qJD(1) + t12 * qJD(2) + ((t116 - t64 + t375) * m(5) + t521) * qJD(3) + t22, qJD(1) * t425 + t22 + (qJD(2) + qJD(3)) * t23;];
Cq = t5;
