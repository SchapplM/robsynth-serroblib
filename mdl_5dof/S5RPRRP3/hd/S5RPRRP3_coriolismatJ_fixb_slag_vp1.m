% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPRRP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:46:56
% EndTime: 2020-01-03 11:47:14
% DurationCPUTime: 9.35s
% Computational Cost: add. (29351->386), mult. (22022->502), div. (0->0), fcn. (20324->8), ass. (0->254)
t579 = Icges(5,4) + Icges(6,4);
t578 = Icges(5,1) + Icges(6,1);
t572 = Icges(5,5) + Icges(6,5);
t576 = Icges(5,2) + Icges(6,2);
t571 = Icges(5,6) + Icges(6,6);
t340 = qJ(3) + qJ(4);
t333 = sin(t340);
t577 = t579 * t333;
t339 = qJ(1) + pkin(8);
t331 = sin(t339);
t332 = cos(t339);
t334 = cos(t340);
t420 = t332 * t334;
t421 = t332 * t333;
t564 = t571 * t331 + t579 * t420 - t576 * t421;
t575 = -t571 * t333 + t572 * t334;
t556 = t578 * t334 - t577;
t574 = t579 * t421;
t573 = t579 * t334;
t570 = Icges(5,3) + Icges(6,3);
t562 = t572 * t331 + t420 * t578 - t574;
t569 = -t333 * t576 + t573;
t563 = t556 * t331 - t572 * t332;
t557 = t576 * t334 + t577;
t568 = t575 * t331;
t567 = t564 * t333;
t566 = t333 * t578 + t573;
t546 = -t570 * t332 + t568;
t545 = t570 * t331 + t572 * t420 - t571 * t421;
t565 = t569 * t331 - t571 * t332;
t561 = -t562 * t334 + t567;
t464 = rSges(6,2) * t334;
t289 = rSges(6,1) * t333 + t464;
t329 = t331 ^ 2;
t225 = t329 * t289;
t248 = t289 * t332;
t330 = t332 ^ 2;
t381 = t329 + t330;
t473 = pkin(4) * t333;
t126 = -t248 * t332 - t381 * t473 - t225;
t341 = sin(qJ(3));
t474 = pkin(3) * t341;
t305 = -t473 - t474;
t353 = t289 - t305;
t177 = t353 * t331;
t179 = t353 * t332;
t465 = rSges(6,1) * t334;
t291 = -rSges(6,2) * t333 + t465;
t472 = pkin(4) * t334;
t374 = -t291 - t472;
t204 = t374 * t331;
t206 = pkin(4) * t420 + t332 * t291;
t405 = -t177 * t204 + t179 * t206;
t325 = t332 * pkin(6);
t343 = cos(qJ(3));
t336 = t343 * pkin(3);
t328 = t336 + pkin(2);
t345 = -pkin(7) - pkin(6);
t417 = t332 * t345;
t183 = t331 * (t417 + t325 + (-pkin(2) + t328) * t331);
t300 = t332 * t328;
t197 = t332 * pkin(2) - t300 + (pkin(6) + t345) * t331;
t338 = -qJ(5) + t345;
t380 = t338 - t345;
t293 = t328 + t472;
t529 = rSges(6,1) * t420 - rSges(6,2) * t421 + t332 * t293;
t402 = -t300 + t529 + (rSges(6,3) - t380) * t331;
t427 = t331 * t333;
t384 = -rSges(6,2) * t427 - t332 * rSges(6,3);
t426 = t331 * t334;
t403 = (rSges(6,1) * t426 + t380 * t332 + t384 + (t293 - t328) * t331) * t331;
t55 = t183 + (-t197 + t402) * t332 + t403;
t32 = t55 * t126 + t405;
t290 = rSges(5,1) * t333 + rSges(5,2) * t334;
t247 = t290 * t331;
t249 = t290 * t332;
t159 = -t331 * t247 - t249 * t332;
t466 = rSges(5,1) * t334;
t292 = -rSges(5,2) * t333 + t466;
t433 = t292 * t332;
t434 = t292 * t331;
t376 = t290 + t474;
t533 = t376 * t332;
t534 = t376 * t331;
t383 = -rSges(5,2) * t427 - t332 * rSges(5,3);
t187 = t331 * (rSges(5,1) * t426 + t383);
t370 = rSges(5,1) * t420 - rSges(5,2) * t421;
t220 = t331 * rSges(5,3) + t370;
t89 = t183 + t187 + (-t197 + t220) * t332;
t44 = t89 * t159 + t433 * t533 + t434 * t534;
t560 = -m(5) * t44 - m(6) * t32;
t559 = t563 * t426;
t558 = t562 * t426;
t555 = t572 * t333 + t571 * t334;
t554 = t576 * t420 - t562 + t574;
t553 = -t557 * t331 + t563;
t552 = t566 * t332 + t564;
t551 = -t566 * t331 - t565;
t550 = t566 + t569;
t549 = t179 * t332;
t548 = t546 * t332 + t565 * t427 - t559;
t547 = -t545 * t332 - t564 * t427 + t558;
t542 = t563 * t334;
t541 = t565 * t333;
t540 = -t545 * t331 + t561 * t332 - t559;
t539 = t546 * t331 + t563 * t420;
t503 = m(5) / 0.2e1;
t502 = m(6) / 0.2e1;
t484 = -t331 / 0.2e1;
t538 = t331 / 0.2e1;
t483 = -t332 / 0.2e1;
t452 = Icges(4,4) * t341;
t308 = Icges(4,2) * t343 + t452;
t311 = Icges(4,1) * t343 - t452;
t532 = (t311 / 0.2e1 - t308 / 0.2e1) * t341;
t531 = t555 * t332;
t530 = t555 * t331;
t278 = t331 * t305;
t181 = -t289 * t331 + t278;
t375 = t289 + t473;
t203 = t375 * t331;
t205 = t375 * t332;
t471 = sin(qJ(1)) * pkin(1);
t160 = t471 + t417 + (t328 + t466) * t331 + t383;
t337 = cos(qJ(1)) * pkin(1);
t161 = t300 + t337 + (rSges(5,3) - t345) * t331 + t370;
t368 = t160 * t433 - t161 * t434;
t157 = t471 + t332 * t338 + (t293 + t465) * t331 + t384;
t158 = t337 + (rSges(6,3) - t338) * t331 + t529;
t406 = t206 * t157 + t204 * t158;
t469 = (t179 * t203 + t181 * t205 + t406) * t502 + ((t331 * t533 - t332 * t534) * t290 + t368) * t503;
t352 = -t464 + (-rSges(6,1) - pkin(4)) * t333;
t195 = t352 * t331;
t196 = t352 * t332;
t470 = (-t177 * t196 + t179 * t195 + t406) * t502 + (-t247 * t533 + t249 * t534 + t368) * t503;
t528 = t469 - t470;
t526 = -t553 * t333 + t551 * t334;
t525 = t554 * t333 - t552 * t334;
t524 = (-t556 + t557) * t334 + t550 * t333;
t335 = Icges(4,4) * t343;
t309 = -Icges(4,2) * t341 + t335;
t519 = Icges(4,1) * t341 + t335;
t419 = t332 * t341;
t318 = Icges(4,4) * t419;
t418 = t332 * t343;
t232 = Icges(4,1) * t418 + Icges(4,5) * t331 - t318;
t389 = -Icges(4,2) * t418 + t232 - t318;
t230 = Icges(4,4) * t418 - Icges(4,2) * t419 + Icges(4,6) * t331;
t391 = t332 * t519 + t230;
t510 = -t341 * t389 - t343 * t391;
t231 = -Icges(4,5) * t332 + t311 * t331;
t390 = -t308 * t331 + t231;
t229 = -Icges(4,6) * t332 + t309 * t331;
t392 = t331 * t519 + t229;
t509 = -t341 * t390 - t343 * t392;
t349 = t550 * t334 / 0.2e1 + (-t557 / 0.2e1 + t556 / 0.2e1) * t333;
t351 = (t547 * t331 + t548 * t332) * t538 + (((t546 - t561) * t332 + t540) * t332 + ((t546 - t567) * t331 + (t541 + t542) * t332 - t539 + t558) * t331) * t484 + (((-t542 + t545) * t332 + t539 + t547) * t332 + ((t541 + t545) * t331 + t540 - t548) * t331) * t483;
t508 = 4 * qJD(1);
t507 = 2 * qJD(3);
t505 = 2 * qJD(4);
t504 = m(4) / 0.2e1;
t467 = rSges(4,1) * t343;
t377 = pkin(2) + t467;
t425 = t331 * t341;
t382 = -rSges(4,2) * t425 - t332 * rSges(4,3);
t169 = t331 * t377 - t325 + t382 + t471;
t320 = rSges(4,2) * t419;
t170 = -t320 + t337 + t377 * t332 + (rSges(4,3) + pkin(6)) * t331;
t312 = rSges(4,1) * t341 + rSges(4,2) * t343;
t272 = t312 * t331;
t273 = t312 * t332;
t500 = m(4) * (-t169 * t272 - t170 * t273);
t146 = t220 * t332 + t187;
t80 = t290 * t292 * t381 + t146 * t159;
t79 = m(5) * t80;
t494 = m(5) * (-t160 * t534 - t161 * t533);
t493 = m(5) * (-t160 * t247 - t161 * t249);
t404 = -t203 * t204 + t205 * t206;
t84 = t332 * t402 + t403;
t99 = t331 * (pkin(3) * t425 + t278) - t225 + (-(-t305 - t474) * t332 - t248) * t332;
t490 = m(6) * (t84 * t99 + t404);
t487 = m(6) * (t157 * t181 - t158 * t179);
t486 = m(6) * (t157 * t195 + t158 * t196);
t485 = m(6) * (-t157 * t332 + t158 * t331);
t482 = t332 / 0.2e1;
t479 = m(6) * (-t181 * t331 + t549);
t478 = m(6) * (-t177 * t331 - t549);
t477 = m(6) * t126;
t476 = m(6) * (-t195 * t331 - t196 * t332);
t475 = m(6) * (-t203 * t331 - t205 * t332);
t468 = m(6) * qJD(3);
t440 = t229 * t341;
t439 = t230 * t341;
t438 = t231 * t343;
t424 = t331 * t343;
t76 = 0.2e1 * (t126 / 0.4e1 - t99 / 0.4e1) * m(6);
t407 = t76 * qJD(2);
t379 = qJD(1) * t485;
t378 = t84 * t126 + t404;
t373 = (t531 * t329 + (t526 * t332 + (-t525 - t530) * t331) * t332) * t484 + (-t530 * t330 + (t525 * t331 + (-t526 + t531) * t332) * t331) * t483;
t372 = t381 * t474;
t371 = t79 + t373;
t364 = Icges(4,5) * t343 - Icges(4,6) * t341;
t363 = Icges(4,5) * t341 + Icges(4,6) * t343;
t355 = t232 * t343 - t439;
t347 = -t351 + (t524 * t332 + t552 * t333 + t554 * t334 - t568) * t484 + (-t524 * t331 - t575 * t332 + t551 * t333 + t553 * t334) * t483;
t346 = -t349 + (t562 * t333 + t564 * t334) * (t482 + t483);
t321 = pkin(3) * t418;
t314 = -rSges(4,2) * t341 + t467;
t267 = t363 * t332;
t266 = t363 * t331;
t228 = Icges(4,5) * t418 - Icges(4,6) * t419 + Icges(4,3) * t331;
t227 = -Icges(4,3) * t332 + t331 * t364;
t224 = t321 + t433;
t222 = (-t292 - t336) * t331;
t194 = t229 * t419;
t193 = t232 * t424;
t192 = t231 * t424;
t180 = t321 + t206;
t178 = (t374 - t336) * t331;
t166 = -t272 * t331 - t273 * t332;
t154 = m(5) * t159;
t141 = -t204 * t332 - t206 * t331;
t135 = t159 - t372;
t134 = t475 / 0.2e1;
t128 = m(6) * t141 * qJD(4);
t125 = t476 / 0.2e1;
t118 = t478 / 0.2e1;
t117 = t479 / 0.2e1;
t115 = t228 * t331 + t332 * t355;
t114 = -t227 * t331 - t231 * t418 + t194;
t113 = t228 * t332 + t230 * t425 - t193;
t112 = -t227 * t332 - t229 * t425 + t192;
t90 = -t372 + t99;
t75 = -t114 * t332 - t115 * t331;
t74 = -t112 * t332 - t113 * t331;
t73 = t134 - t476 / 0.2e1;
t72 = t134 + t125;
t71 = t125 - t475 / 0.2e1;
t48 = t154 + t477 / 0.2e1 + t99 * t502;
t47 = t118 + t117;
t46 = t118 - t479 / 0.2e1;
t45 = t117 - t478 / 0.2e1;
t28 = t349 + t486 + t493;
t25 = (-t113 + t194 + (t228 - t438) * t332) * t332 + (t112 - t192 + (t228 + t440) * t331) * t331;
t24 = (t114 + t193 - t194 + (t227 - t439) * t331) * t331 + (-t192 - t115 + (t227 + t355) * t332 + (t438 + t440) * t331) * t332;
t18 = (t519 / 0.2e1 + t309 / 0.2e1) * t343 + t532 + t500 + t494 + t487 + t349;
t7 = t371 + t490;
t6 = t373 - t560;
t4 = t351 + t528;
t3 = t351 - t528;
t2 = (-t75 / 0.2e1 - t25 / 0.2e1) * t332 + (-t24 / 0.2e1 + t74 / 0.2e1) * t331 + t351;
t1 = t347 + t469 + t470;
t5 = [t18 * qJD(3) + t28 * qJD(4) + qJD(5) * t485, 0, t18 * qJD(1) + t1 * qJD(4) + t47 * qJD(5) + (((t169 * t332 - t170 * t331) * t314 + (-t272 * t332 + t273 * t331) * t312) * t504 + (t160 * t224 + t161 * t222) * t503 + (t157 * t180 + t158 * t178 + (t177 + t181) * t179) * t502) * t507 + (t347 + (-t341 * t392 + t343 * t390) * t483 + t24 * t538 + (t341 * t391 - t343 * t389 + t74) * t484 + (t75 + t25) * t482 + (t330 / 0.2e1 + t329 / 0.2e1) * t364) * qJD(3), t28 * qJD(1) + t1 * qJD(3) + t347 * qJD(4) + t72 * qJD(5) + (((-t247 * t332 + t249 * t331) * t290 + t368) * t503 + (t195 * t205 - t196 * t203 + t406) * t502) * t505, t47 * qJD(3) + t72 * qJD(4) + t379; 0, 0, t48 * qJD(4) + (t135 * t503 + t166 * t504 + t90 * t502) * t507, t48 * qJD(3) + (t154 + t477) * qJD(4), 0; (t346 - (t309 + t519) * t343 / 0.2e1 - t532) * qJD(1) + t2 * qJD(3) + t3 * qJD(4) + t46 * qJD(5) + (-t500 / 0.4e1 - t487 / 0.4e1 - t494 / 0.4e1) * t508, qJD(4) * t76, t2 * qJD(1) + (m(6) * (-t177 * t178 + t179 * t180 + t55 * t90) + m(5) * (t135 * t89 - t222 * t534 + t224 * t533) + (t329 * t267 + (t509 * t332 + (-t266 - t510) * t331) * t332) * t484 + (-t330 * t266 + (t510 * t331 + (t267 - t509) * t332) * t331) * t483 + m(4) * (t312 * t314 * t381 + (t332 * (rSges(4,1) * t418 + rSges(4,3) * t331 - t320) + t331 * (rSges(4,1) * t424 + t382)) * t166) + t373) * qJD(3) + t6 * qJD(4), t3 * qJD(1) + t407 + t6 * qJD(3) + ((t378 + t32) * t502 + (t44 + t80) * t503) * t505 + (t373 - t79 - t490) * qJD(4), t46 * qJD(1); t346 * qJD(1) + t4 * qJD(3) + t351 * qJD(4) + t73 * qJD(5) + (-t493 / 0.4e1 - t486 / 0.4e1) * t508, -qJD(3) * t76, t4 * qJD(1) - t407 + t7 * qJD(4) + ((-t178 * t203 + t180 * t205 + t55 * t99 + t84 * t90 + t405) * t502 + (t135 * t146 + (-t222 * t331 + t224 * t332) * t290 + t44) * t503) * t507 + (t373 + t560) * qJD(3), t351 * qJD(1) + t7 * qJD(3) + (m(6) * t378 + t371) * qJD(4), t73 * qJD(1); t45 * qJD(3) + t71 * qJD(4) - t379, 0, t45 * qJD(1) + (-t178 * t332 - t180 * t331) * t468 + t128, t71 * qJD(1) + t141 * t468 + t128, 0;];
Cq = t5;
