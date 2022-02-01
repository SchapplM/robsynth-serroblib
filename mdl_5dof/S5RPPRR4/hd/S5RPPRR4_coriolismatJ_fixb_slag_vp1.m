% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RPPRR4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:18
% EndTime: 2022-01-23 09:16:35
% DurationCPUTime: 10.02s
% Computational Cost: add. (41539->464), mult. (40405->672), div. (0->0), fcn. (43710->10), ass. (0->288)
t370 = sin(pkin(9));
t484 = pkin(3) * t370;
t358 = qJ(2) + t484;
t375 = sin(qJ(1));
t350 = t358 * t375;
t374 = qJ(3) + pkin(6);
t367 = -pkin(7) - t374;
t371 = sin(pkin(8));
t376 = cos(qJ(1));
t459 = t371 * t376;
t349 = t367 * t459;
t372 = cos(pkin(9));
t361 = t372 * pkin(3) + pkin(2);
t373 = cos(pkin(8));
t340 = t361 * t373 + t371 * t374 + pkin(1);
t369 = pkin(9) + qJ(4);
t363 = cos(t369);
t345 = pkin(4) * t363 + t361;
t401 = t345 * t373 + pkin(1);
t397 = -t340 + t401;
t383 = t397 * t376 - t349;
t362 = sin(t369);
t483 = pkin(4) * t362;
t346 = t483 + t484;
t448 = qJ(2) + t346;
t399 = t448 * t375;
t196 = -t350 + t399 + t383;
t364 = qJ(5) + t369;
t359 = sin(t364);
t450 = t376 * t359;
t360 = cos(t364);
t453 = t375 * t360;
t314 = -t373 * t450 + t453;
t454 = t375 * t359;
t455 = t373 * t376;
t315 = t360 * t455 + t454;
t394 = t315 * rSges(6,1) + t314 * rSges(6,2);
t245 = rSges(6,3) * t459 + t394;
t277 = (t367 + t374) * t373 + (t345 - t361) * t371;
t297 = -rSges(6,3) * t373 + (rSges(6,1) * t360 - rSges(6,2) * t359) * t371;
t539 = (t277 + t297) * t459;
t100 = (t196 + t245) * t373 + t539;
t195 = (-t358 + t448) * t375 + t383;
t441 = -t195 - t245;
t101 = t441 * t373 - t539;
t569 = t100 + t101;
t472 = Icges(5,4) * t363;
t299 = -Icges(5,6) * t373 + (-Icges(5,2) * t362 + t472) * t371;
t473 = Icges(5,4) * t362;
t300 = -Icges(5,5) * t373 + (Icges(5,1) * t363 - t473) * t371;
t327 = (-Icges(5,2) * t363 - t473) * t371;
t328 = (-Icges(5,1) * t362 - t472) * t371;
t326 = (-Icges(5,5) * t362 - Icges(5,6) * t363) * t371;
t457 = t373 * t326;
t568 = t457 / 0.2e1 - (-(t300 / 0.2e1 + t327 / 0.2e1) * t362 + (t328 / 0.2e1 - t299 / 0.2e1) * t363) * t371;
t388 = rSges(6,3) * t371 + t401;
t169 = t388 * t376 - t349 + t394 + t399;
t312 = t360 * t376 + t373 * t454;
t313 = t373 * t453 - t450;
t266 = -t312 * rSges(6,1) - t313 * rSges(6,2);
t267 = t314 * rSges(6,1) - t315 * rSges(6,2);
t366 = t376 * qJ(2);
t460 = t371 * t375;
t413 = t376 * t346 + t367 * t460 + t366;
t535 = -t313 * rSges(6,1) + t312 * rSges(6,2);
t548 = -t388 * t375 + t413 + t535;
t567 = m(6) * (t169 * t267 - t266 * t548);
t243 = rSges(6,3) * t460 - t535;
t286 = t297 * t460;
t565 = t373 * t243 + t286;
t452 = t375 * t362;
t329 = t363 * t376 + t373 * t452;
t449 = t376 * t362;
t451 = t375 * t363;
t330 = t373 * t451 - t449;
t534 = -t330 * rSges(5,1) + t329 * rSges(5,2);
t257 = rSges(5,3) * t460 - t534;
t301 = -rSges(5,3) * t373 + (rSges(5,1) * t363 - rSges(5,2) * t362) * t371;
t564 = t257 * t373 + t301 * t460;
t563 = -t371 / 0.2e1;
t470 = Icges(6,4) * t359;
t296 = -Icges(6,5) * t373 + (Icges(6,1) * t360 - t470) * t371;
t428 = t296 + (-Icges(6,2) * t360 - t470) * t371;
t562 = t428 * t359;
t303 = Icges(6,4) * t313;
t237 = -Icges(6,2) * t312 + Icges(6,6) * t460 + t303;
t302 = Icges(6,4) * t312;
t241 = -Icges(6,1) * t313 - Icges(6,5) * t460 + t302;
t560 = t314 * t237 - t315 * t241;
t318 = Icges(5,4) * t330;
t251 = -Icges(5,2) * t329 + Icges(5,6) * t460 + t318;
t317 = Icges(5,4) * t329;
t255 = -Icges(5,1) * t330 - Icges(5,5) * t460 + t317;
t332 = t363 * t455 + t452;
t546 = t373 * t449 - t451;
t558 = t251 * t546 + t332 * t255;
t523 = m(5) / 0.2e1;
t521 = m(6) / 0.2e1;
t234 = Icges(6,5) * t313 - Icges(6,6) * t312 + Icges(6,3) * t460;
t553 = t234 * t459;
t105 = t553 + t560;
t538 = t105 - t553;
t556 = (t538 - t560) * t459;
t248 = Icges(5,5) * t330 - Icges(5,6) * t329 + Icges(5,3) * t460;
t551 = t248 * t459;
t431 = t546 * pkin(4) - t267;
t550 = t431 * t375;
t463 = t358 * t376;
t547 = t463 + (-rSges(5,3) * t371 - t340) * t375 + t534;
t259 = t332 * rSges(5,1) - rSges(5,2) * t546 + rSges(5,3) * t459;
t190 = t373 * t259 + t301 * t459;
t194 = t397 * t375 - t413 + t463;
t268 = t277 * t460;
t99 = t194 * t373 + t268 + t565;
t67 = t101 * t375 + t99 * t376;
t475 = (-t190 * t375 + t376 * t564) * t523 + t67 * t521;
t98 = t268 + t286 + (t194 + t243) * t373;
t481 = (t100 * t375 - t376 * t98 + t67) * t521;
t545 = t475 - t481;
t92 = t101 * t459;
t476 = (-t190 * t459 - t460 * t564) * t523 + (-t99 * t460 + t92) * t521;
t480 = t98 - t99;
t63 = (-t195 + t196) * t459;
t482 = (-t63 * t373 + t92 + (t100 * t376 + t480 * t375) * t371) * t521;
t544 = t476 - t482;
t309 = (-Icges(6,5) * t359 - Icges(6,6) * t360) * t371;
t458 = t373 * t309;
t543 = -t458 / 0.2e1 + t562 * t563;
t542 = t371 / 0.2e1;
t501 = -t373 / 0.2e1;
t541 = m(6) * t371;
t320 = t329 * pkin(4);
t246 = t266 * t459;
t164 = -t267 * t460 + t246;
t316 = (-rSges(6,1) * t359 - rSges(6,2) * t360) * t371;
t201 = t373 * t266 + t316 * t460;
t202 = -t373 * t267 - t316 * t459;
t469 = Icges(6,4) * t360;
t295 = -Icges(6,6) * t373 + (-Icges(6,2) * t359 + t469) * t371;
t311 = (-Icges(6,1) * t359 - t469) * t371;
t429 = -t295 + t311;
t107 = t309 * t460 - t428 * t312 + t429 * t313;
t108 = t309 * t459 + t428 * t314 + t429 * t315;
t260 = -Icges(6,5) * t312 - Icges(6,6) * t313;
t261 = Icges(6,5) * t314 - Icges(6,6) * t315;
t405 = t459 / 0.2e1;
t407 = t460 / 0.2e1;
t304 = Icges(6,4) * t314;
t242 = Icges(6,1) * t315 + Icges(6,5) * t459 + t304;
t437 = -Icges(6,2) * t315 + t242 + t304;
t438 = -Icges(6,2) * t313 - t241 - t302;
t471 = Icges(6,4) * t315;
t239 = Icges(6,2) * t314 + Icges(6,6) * t459 + t471;
t439 = Icges(6,1) * t314 - t239 - t471;
t440 = -Icges(6,1) * t312 - t237 - t303;
t138 = -t458 + (t429 * t360 - t562) * t371;
t466 = t138 * t373;
t84 = -t260 * t373 + (-t438 * t359 + t440 * t360) * t371;
t85 = -t261 * t373 + (-t437 * t359 + t439 * t360) * t371;
t422 = ((t261 * t459 + t437 * t314 + t439 * t315) * t459 + (t260 * t459 + t438 * t314 + t440 * t315) * t460 - t108 * t373) * t405 + ((t261 * t460 - t437 * t312 + t439 * t313) * t459 + (t260 * t460 - t438 * t312 + t440 * t313) * t460 - t107 * t373) * t407 + (-t466 + (t375 * t84 + t376 * t85) * t371) * t501;
t229 = t243 * t459;
t83 = t229 + (t194 * t376 + t441 * t375) * t371;
t6 = t422 + m(6) * (t101 * t202 + t83 * t164 + t99 * t201);
t540 = t6 * qJD(5);
t275 = -rSges(5,1) * t329 - rSges(5,2) * t330;
t276 = -rSges(5,1) * t546 - rSges(5,2) * t332;
t174 = (t275 * t376 - t276 * t375) * t371;
t112 = t551 - t558;
t537 = t112 - t551;
t199 = t340 * t376 + t259 + t350;
t217 = -t266 + t320;
t445 = (t375 * t217 + t376 * t431) * t521 + (-t375 * t275 - t276 * t376) * t523;
t446 = (t217 * t376 - t550) * t541 / 0.2e1 - m(5) * t174 / 0.2e1;
t344 = (-t375 ^ 2 - t376 ^ 2) * t371;
t530 = 0.2e1 * t344;
t529 = 2 * qJD(1);
t528 = 4 * qJD(1);
t527 = 2 * qJD(4);
t526 = 4 * qJD(4);
t525 = m(4) / 0.2e1;
t522 = m(5) / 0.4e1;
t520 = m(6) / 0.4e1;
t516 = m(5) * (t199 * t276 - t275 * t547);
t152 = -t245 * t460 + t229;
t183 = t373 * t245 + t297 * t459;
t515 = m(6) * (t152 * t63 - t480 * t183 + t569 * t565);
t514 = m(6) * (t100 * t99 + t101 * t98 + t63 * t83);
t447 = t202 * t169 + t201 * t548;
t511 = m(6) * (t101 * t267 - t266 * t99 + t447);
t509 = m(6) * (t183 * t431 + t217 * t565 + t447);
t504 = m(6) * (-t169 * t431 + t217 * t548);
t161 = t169 * t459;
t503 = m(6) * (-t460 * t548 + t161);
t500 = m(3) * ((rSges(3,2) * t460 + rSges(3,3) * t376 + t366) * t376 + (-rSges(3,2) * t459 + (rSges(3,3) + qJ(2)) * t375) * t375);
t400 = -pkin(2) * t373 - pkin(1) + (-rSges(4,3) - qJ(3)) * t371;
t226 = (t370 * rSges(4,1) + t372 * rSges(4,2) + qJ(2)) * t375 + ((rSges(4,1) * t372 - rSges(4,2) * t370) * t373 - t400) * t376;
t219 = t226 * t459;
t456 = t373 * t375;
t377 = -rSges(4,1) * (-t370 * t376 + t372 * t456) + rSges(4,2) * (t370 * t456 + t372 * t376) + t400 * t375 + t366;
t499 = m(4) * (-t377 * t460 + t219);
t498 = m(4) * (t226 * t375 + t376 * t377);
t187 = t199 * t459;
t495 = m(5) * (-t460 * t547 + t187);
t494 = m(5) * (t199 * t375 + t376 * t547);
t491 = m(6) * (t169 * t375 + t376 * t548);
t490 = m(6) * (-t183 * t459 - t460 * t565);
t489 = m(6) * (-t183 * t375 + t376 * t565);
t486 = (-t266 * t376 + t267 * t375) * t541;
t485 = m(6) * (-t375 * t266 - t267 * t376);
t479 = m(6) * qJD(4);
t478 = m(6) * qJD(5);
t474 = Icges(5,4) * t332;
t465 = ((-Icges(6,3) * t373 + (Icges(6,5) * t360 - Icges(6,6) * t359) * t371) * t459 + t314 * t295 + t315 * t296) * t373;
t461 = t360 * t371;
t435 = -Icges(5,1) * t329 - t251 - t318;
t253 = -Icges(5,2) * t546 + Icges(5,6) * t459 + t474;
t434 = -Icges(5,1) * t546 - t253 - t474;
t433 = -Icges(5,2) * t330 - t255 - t317;
t319 = Icges(5,4) * t546;
t256 = Icges(5,1) * t332 + Icges(5,5) * t459 - t319;
t432 = -Icges(5,2) * t332 + t256 - t319;
t427 = -t299 + t328;
t426 = t300 + t327;
t173 = (t521 + t523 + t525) * t530;
t423 = t173 * qJD(1);
t421 = t371 ^ 2 * t483;
t418 = (t537 + t558) * t542 * t376 + (t373 / 0.2e1 + t501) * ((-Icges(5,3) * t373 + (Icges(5,5) * t363 - Icges(5,6) * t362) * t371) * t460 - t299 * t329 + t300 * t330);
t113 = (Icges(5,5) * t332 - Icges(5,6) * t546 + Icges(5,3) * t459) * t459 - t546 * t253 + t332 * t256;
t417 = (t112 * t375 + t113 * t376) * t563 + ((t248 * t460 + t113) * t376 + t537 * t375) * t542;
t106 = (Icges(6,5) * t315 + Icges(6,6) * t314 + Icges(6,3) * t459) * t459 + t314 * t239 + t315 * t242;
t410 = -t461 / 0.2e1;
t409 = t461 / 0.2e1;
t408 = -t460 / 0.2e1;
t23 = -t465 + ((t234 * t460 + t106) * t376 + t538 * t375) * t371;
t52 = -t465 + (t105 * t375 + t106 * t376) * t371;
t398 = t23 * t407 + t556 * t405 + t52 * t408;
t393 = t295 * t410 + t311 * t409 + t543;
t392 = t515 / 0.2e1 + t398;
t382 = t23 * t408 - t556 * t459 / 0.2e1 + (t85 + t108) * t405 + (t52 + t84 + t107) * t407;
t380 = t295 * t409 + t311 * t410 - t543;
t379 = t382 - t466;
t333 = (-rSges(5,1) * t362 - rSges(5,2) * t363) * t371;
t270 = -Icges(5,5) * t546 - Icges(5,6) * t332;
t269 = -Icges(5,5) * t329 - Icges(5,6) * t330;
t210 = -t373 * t276 - t333 * t459;
t209 = t275 * t373 + t333 * t460;
t172 = (t520 + t522 + m(4) / 0.4e1) * t530 - (m(6) + m(5) + m(4)) * t344 / 0.2e1;
t166 = t485 / 0.2e1;
t160 = t486 / 0.2e1;
t157 = (-t316 * t371 + t421) * t376 + t431 * t373;
t156 = -t320 * t373 - t375 * t421 + t201;
t142 = -t457 + (-t426 * t362 + t427 * t363) * t371;
t137 = t246 + (-t320 * t376 + t550) * t371;
t136 = t201 * t375 - t202 * t376;
t132 = t136 * t478;
t119 = t489 / 0.2e1;
t118 = t326 * t459 + t427 * t332 - t426 * t546;
t117 = t326 * t460 - t426 * t329 + t427 * t330;
t109 = t490 / 0.2e1;
t89 = -t270 * t373 + (-t432 * t362 + t434 * t363) * t371;
t88 = -t269 * t373 + (-t433 * t362 + t435 * t363) * t371;
t87 = -t164 * t373 + (t201 * t376 + t202 * t375) * t371;
t86 = t87 * t478;
t62 = t393 + t567;
t55 = t495 + t499 + t503;
t46 = t491 + t494 + t498 + t500;
t44 = t119 - t485 / 0.2e1;
t43 = t166 + t119;
t42 = t166 - t489 / 0.2e1;
t40 = t509 / 0.2e1;
t38 = t109 - t486 / 0.2e1;
t37 = t160 + t109;
t36 = t160 - t490 / 0.2e1;
t32 = t511 / 0.2e1;
t31 = t516 + t504 + t393 - t568;
t16 = t475 + t481 - t445;
t15 = t445 - t545;
t14 = t445 + t545;
t11 = t476 + t482 - t446;
t10 = t446 - t544;
t9 = t446 + t544;
t8 = m(6) * (t152 * t164 - t183 * t202 + t201 * t565) + t422;
t7 = t8 * qJD(5);
t4 = t40 - t511 / 0.2e1 + t392;
t3 = t32 - t509 / 0.2e1 + t392;
t2 = t32 + t40 - t515 / 0.2e1 + t379;
t1 = t514 + (t417 * t375 + t418 * t376) * t371 + t398;
t5 = [t46 * qJD(2) + t55 * qJD(3) + t31 * qJD(4) + t62 * qJD(5), qJD(1) * t46 + qJD(3) * t172 + qJD(4) * t14 + qJD(5) * t43, qJD(1) * t55 + qJD(2) * t172 + qJD(4) * t9 + qJD(5) * t37, t31 * qJD(1) + t14 * qJD(2) + t9 * qJD(3) + t2 * qJD(5) - t514 * t526 / 0.4e1 + ((-t101 * t431 + t156 * t548 + t157 * t169 + t217 * t99) * t521 + (-t190 * t276 + t199 * t210 + t209 * t547 - t275 * t564) * t523) * t527 + (t382 + (-t138 - t142) * t373 + ((t89 / 0.2e1 + t118 / 0.2e1 - t418) * t376 + (t88 / 0.2e1 + t117 / 0.2e1 - t417) * t375) * t371) * qJD(4), t62 * qJD(1) + t43 * qJD(2) + t37 * qJD(3) + t2 * qJD(4) + ((-t183 * t267 - t266 * t565 + t447) * m(6) + t379) * qJD(5); t173 * qJD(3) + t15 * qJD(4) + t42 * qJD(5) + (-t491 / 0.4e1 - t494 / 0.4e1 - t500 / 0.4e1 - t498 / 0.4e1) * t528, 0, t423, t15 * qJD(1) + ((t209 * t375 - t210 * t376) * t523 + (t156 * t375 - t157 * t376) * t521) * t527 + t132, t42 * qJD(1) + t136 * t479 + t132; -t173 * qJD(2) + t10 * qJD(4) + t36 * qJD(5) + (-t503 / 0.4e1 - t495 / 0.4e1 - t499 / 0.4e1) * t528 + (t161 * t521 + t187 * t523 + t219 * t525 + (-t169 * t521 - t199 * t523 - t226 * t525) * t459) * t529, -t423, 0, t10 * qJD(1) + ((-t174 * t373 + (t209 * t376 + t210 * t375) * t371) * t523 + (-t137 * t373 + (t156 * t376 + t157 * t375) * t371) * t521) * t527 + t86, t36 * qJD(1) + t479 * t87 + t86; t16 * qJD(2) + t11 * qJD(3) + t1 * qJD(4) + t3 * qJD(5) + (-t504 / 0.4e1 - t516 / 0.4e1) * t528 + (t480 * t169 + t569 * t548) * t521 * t529 + (t380 + t568) * qJD(1), t16 * qJD(1), t11 * qJD(1), t1 * qJD(1) + (((t270 * t459 + t434 * t332 - t432 * t546) * t459 + (t269 * t459 + t435 * t332 - t433 * t546) * t460 - t118 * t373) * t405 + ((t270 * t460 - t432 * t329 + t434 * t330) * t459 + (t269 * t460 - t433 * t329 + t435 * t330) * t460 - t117 * t373) * t407 + (-t142 * t373 + (t88 * t375 + t89 * t376) * t371) * t501 + t422) * qJD(4) + t540 + (((t257 * t376 - t259 * t375) * t371 * t174 - t190 * t210 + t209 * t564) * t522 + (t101 * t157 + t137 * t83 + t156 * t99) * t520) * t526, t3 * qJD(1) + t6 * qJD(4) + t540; (t380 - t567) * qJD(1) + t44 * qJD(2) + t38 * qJD(3) + t4 * qJD(4) + t398 * qJD(5), t44 * qJD(1), t38 * qJD(1), t4 * qJD(1) + ((t137 * t152 + t156 * t565 - t157 * t183) * m(6) + t422) * qJD(4) + t7, qJD(1) * t398 + qJD(4) * t8 + t7;];
Cq = t5;
