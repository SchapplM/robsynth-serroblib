% Calculate time derivative of joint inertia matrix for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR15_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_inertiaDJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR15_inertiaDJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR15_inertiaDJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:12
% EndTime: 2019-12-31 20:41:48
% DurationCPUTime: 22.31s
% Computational Cost: add. (23780->922), mult. (40365->1285), div. (0->0), fcn. (38620->8), ass. (0->440)
t350 = cos(qJ(2));
t347 = sin(qJ(2));
t530 = Icges(4,6) * t347;
t540 = Icges(3,4) * t347;
t593 = -t530 - t540 + (-Icges(3,2) - Icges(4,3)) * t350;
t529 = Icges(4,6) * t350;
t539 = Icges(3,4) * t350;
t592 = -t529 - t539 + (-Icges(3,1) - Icges(4,2)) * t347;
t591 = t593 * qJD(2);
t590 = t592 * qJD(2);
t348 = sin(qJ(1));
t351 = cos(qJ(1));
t410 = -Icges(3,2) * t347 + t539;
t249 = -Icges(3,6) * t351 + t348 * t410;
t400 = -Icges(4,3) * t347 + t529;
t572 = Icges(4,5) * t351 + t348 * t400;
t589 = -t249 - t572;
t250 = Icges(3,6) * t348 + t351 * t410;
t254 = Icges(4,5) * t348 - t351 * t400;
t588 = t250 - t254;
t414 = Icges(3,1) * t350 - t540;
t252 = -Icges(3,5) * t351 + t348 * t414;
t402 = Icges(4,2) * t350 - t530;
t571 = Icges(4,4) * t351 + t348 * t402;
t587 = t252 + t571;
t253 = Icges(3,5) * t348 + t351 * t414;
t256 = Icges(4,4) * t348 - t351 * t402;
t586 = t253 - t256;
t486 = qJD(2) * t347;
t454 = -t486 / 0.2e1;
t488 = qJD(1) * t350;
t466 = t348 * t488;
t585 = t351 * t454 - t466 / 0.2e1;
t487 = qJD(1) * t351;
t455 = t487 / 0.2e1;
t584 = t348 * t454 + t350 * t455;
t483 = qJD(2) * t351;
t463 = t347 * t483;
t583 = t463 + t466;
t349 = cos(qJ(4));
t346 = sin(qJ(4));
t537 = Icges(5,4) * t346;
t407 = Icges(5,2) * t349 + t537;
t248 = Icges(5,6) * t347 - t350 * t407;
t536 = Icges(5,4) * t349;
t412 = Icges(5,1) * t346 + t536;
t251 = Icges(5,5) * t347 - t350 * t412;
t582 = t248 * t349 + t251 * t346;
t345 = qJ(4) + qJ(5);
t334 = cos(t345);
t333 = sin(t345);
t535 = Icges(6,4) * t333;
t406 = Icges(6,2) * t334 + t535;
t240 = Icges(6,6) * t347 - t350 * t406;
t534 = Icges(6,4) * t334;
t411 = Icges(6,1) * t333 + t534;
t241 = Icges(6,5) * t347 - t350 * t411;
t581 = t240 * t334 + t241 * t333;
t580 = qJD(2) / 0.2e1;
t386 = t254 * t347 - t256 * t350;
t579 = t348 * t386;
t387 = t250 * t347 - t253 * t350;
t578 = t348 * t387;
t385 = -t347 * t572 + t350 * t571;
t577 = t351 * t385;
t388 = t249 * t347 - t252 * t350;
t576 = t388 * t351;
t517 = t347 * t348;
t268 = t333 * t351 + t334 * t517;
t269 = t333 * t517 - t334 * t351;
t428 = -t269 * rSges(6,1) - t268 * rSges(6,2);
t514 = t348 * t350;
t179 = rSges(6,3) * t514 - t428;
t340 = t351 * pkin(3);
t352 = -pkin(8) - pkin(7);
t552 = -pkin(7) - t352;
t554 = pkin(4) * t346;
t365 = t347 * t554 + t350 * t552;
t328 = pkin(4) * t349 + pkin(3);
t511 = t351 * t328;
t218 = t348 * t365 + t340 - t511;
t504 = t179 + t218;
t575 = t504 * t351;
t512 = t350 * t351;
t574 = t348 * rSges(4,1) - rSges(4,2) * t512;
t294 = t348 * pkin(3) + pkin(7) * t512;
t516 = t347 * t351;
t573 = -rSges(3,2) * t516 + t348 * rSges(3,3);
t405 = Icges(3,5) * t350 - Icges(3,6) * t347;
t246 = -Icges(3,3) * t351 + t348 * t405;
t408 = Icges(4,4) * t350 - Icges(4,5) * t347;
t570 = Icges(4,1) * t351 + t348 * t408;
t266 = -t333 * t348 + t334 * t516;
t267 = t333 * t516 + t334 * t348;
t178 = t267 * rSges(6,1) + t266 * rSges(6,2) + rSges(6,3) * t512;
t474 = t346 * t516;
t376 = pkin(4) * t474 + t348 * t328 - t352 * t512;
t217 = t376 - t294;
t505 = t178 + t217;
t569 = -t348 * t504 - t351 * t505;
t568 = 2 * m(3);
t567 = 2 * m(4);
t566 = 2 * m(5);
t565 = 2 * m(6);
t564 = m(4) / 0.2e1;
t563 = m(5) / 0.2e1;
t562 = m(6) / 0.2e1;
t561 = t347 / 0.2e1;
t560 = t348 / 0.2e1;
t559 = -t351 / 0.2e1;
t558 = t351 / 0.2e1;
t557 = -rSges(6,3) - pkin(2);
t309 = rSges(3,1) * t347 + rSges(3,2) * t350;
t556 = m(3) * t309;
t555 = pkin(2) * t350;
t553 = qJD(1) / 0.2e1;
t551 = rSges(4,1) * t351;
t550 = rSges(4,2) * t347;
t549 = rSges(3,3) * t351;
t548 = rSges(5,3) * t347;
t547 = pkin(4) * qJD(4);
t546 = -rSges(4,3) - qJ(3);
t545 = rSges(6,3) - t352;
t403 = Icges(6,5) * t333 + Icges(6,6) * t334;
t239 = Icges(6,3) * t347 - t350 * t403;
t132 = t239 * t347 - t581 * t350;
t342 = qJD(4) + qJD(5);
t519 = t342 * t350;
t166 = (Icges(6,2) * t333 - t534) * t519 + (Icges(6,6) * t350 + t347 * t406) * qJD(2);
t165 = (-Icges(6,5) * t334 + Icges(6,6) * t333) * t519 + (Icges(6,3) * t350 + t347 * t403) * qJD(2);
t484 = qJD(2) * t350;
t417 = t333 * t240 * t519 + t347 * t165 + t239 * t484 + t581 * t486;
t167 = (-Icges(6,1) * t334 + t535) * t519 + (Icges(6,5) * t350 + t347 * t411) * qJD(2);
t525 = t167 * t333;
t544 = t132 * t484 + ((-t525 + (-t241 * t342 - t166) * t334) * t350 + t417) * t347;
t404 = Icges(5,5) * t346 + Icges(5,6) * t349;
t245 = Icges(5,3) * t347 - t350 * t404;
t141 = t245 * t347 - t582 * t350;
t481 = qJD(4) * t350;
t204 = (Icges(5,2) * t346 - t536) * t481 + (Icges(5,6) * t350 + t347 * t407) * qJD(2);
t207 = (-Icges(5,1) * t349 + t537) * t481 + (Icges(5,5) * t350 + t347 * t412) * qJD(2);
t201 = (-Icges(5,5) * t349 + Icges(5,6) * t346) * t481 + (Icges(5,3) * t350 + t347 * t404) * qJD(2);
t416 = t346 * t248 * t481 + t347 * t201 + t245 * t484 + t582 * t486;
t543 = t141 * t484 + ((-t207 * t346 + (-qJD(4) * t251 - t204) * t349) * t350 + t416) * t347;
t485 = qJD(2) * t348;
t464 = t347 * t485;
t465 = t350 * t487;
t364 = -t464 + t465;
t490 = qJD(1) * t347;
t449 = t342 + t490;
t461 = t348 * t484;
t359 = t351 * t449 + t461;
t450 = t342 * t347 + qJD(1);
t384 = t333 * t450;
t156 = t334 * t359 - t348 * t384;
t383 = t334 * t450;
t157 = t333 * t359 + t348 * t383;
t429 = t157 * rSges(6,1) + t156 * rSges(6,2);
t103 = rSges(6,3) * t364 + t429;
t542 = t103 * t512 + t178 * t464;
t527 = qJ(3) * t347;
t526 = qJ(3) * t350;
t515 = t348 * t349;
t279 = t346 * t351 + t347 * t515;
t513 = t349 * t351;
t280 = t346 * t517 - t513;
t431 = -rSges(5,1) * t280 - rSges(5,2) * t279;
t198 = rSges(5,3) * t514 - t431;
t524 = t198 * t351;
t518 = t346 * t350;
t317 = pkin(7) * t464;
t462 = t352 * t486;
t128 = t348 * t462 + t317 + (qJD(4) * t279 + t346 * t461) * pkin(4) + ((-pkin(3) + t328) * t348 + t365 * t351) * qJD(1);
t510 = t103 + t128;
t460 = t350 * t483;
t358 = -t348 * t449 + t460;
t158 = t334 * t358 - t351 * t384;
t159 = t333 * t358 + t351 * t383;
t508 = t159 * rSges(6,1) + t158 * rSges(6,2);
t104 = -rSges(6,3) * t583 + t508;
t332 = pkin(3) * t487;
t473 = t347 * t513;
t415 = t328 * t487 + t351 * t462 + t352 * t466 + t460 * t554 + t473 * t547;
t445 = qJD(4) + t490;
t129 = pkin(7) * t463 - t332 + (pkin(7) * t488 - t445 * t554) * t348 + t415;
t509 = t104 + t129;
t427 = rSges(6,1) * t333 + rSges(6,2) * t334;
t169 = (-rSges(6,1) * t334 + rSges(6,2) * t333) * t519 + (rSges(6,3) * t350 + t347 * t427) * qJD(2);
t242 = rSges(6,3) * t347 - t350 * t427;
t507 = t169 * t514 + t242 * t465;
t459 = t349 * t481;
t225 = -pkin(4) * t459 + qJD(2) * t365;
t506 = -t169 - t225;
t357 = -t348 * t445 + t460;
t446 = qJD(4) * t347 + qJD(1);
t380 = t446 * t346;
t184 = t349 * t357 - t351 * t380;
t381 = t349 * t446;
t185 = t346 * t357 + t351 * t381;
t503 = t185 * rSges(5,1) + t184 * rSges(5,2);
t272 = -pkin(4) * t518 + t347 * t552;
t502 = t242 + t272;
t424 = t527 + t555;
t283 = t424 * t348;
t284 = pkin(2) * t512 + qJ(3) * t516;
t501 = t348 * t283 + t351 * t284;
t273 = qJD(2) * t424 - qJD(3) * t350;
t426 = -rSges(4,2) * t350 + rSges(4,3) * t347;
t500 = -t426 * qJD(2) - t273;
t499 = -t284 - t294;
t307 = pkin(2) * t347 - t526;
t489 = qJD(1) * t348;
t285 = t307 * t489;
t467 = t347 * t489;
t498 = pkin(7) * t467 + t285;
t425 = rSges(4,3) * t350 + t550;
t497 = -t307 + t425;
t482 = qJD(3) * t347;
t496 = qJ(3) * t460 + t351 * t482;
t495 = rSges(3,2) * t467 + rSges(3,3) * t487;
t494 = t351 * pkin(1) + t348 * pkin(6);
t493 = t348 ^ 2 + t351 ^ 2;
t247 = Icges(3,3) * t348 + t351 * t405;
t492 = qJD(1) * t247;
t258 = Icges(4,1) * t348 - t351 * t408;
t491 = qJD(1) * t258;
t480 = -rSges(5,3) - pkin(2) - pkin(7);
t477 = t346 * t547;
t114 = Icges(5,5) * t185 + Icges(5,6) * t184 - Icges(5,3) * t583;
t116 = Icges(5,4) * t185 + Icges(5,2) * t184 - Icges(5,6) * t583;
t118 = Icges(5,1) * t185 + Icges(5,4) * t184 - Icges(5,5) * t583;
t277 = -t346 * t348 + t473;
t278 = t474 + t515;
t191 = Icges(5,5) * t278 + Icges(5,6) * t277 + Icges(5,3) * t512;
t193 = Icges(5,4) * t278 + Icges(5,2) * t277 + Icges(5,6) * t512;
t195 = Icges(5,1) * t278 + Icges(5,4) * t277 + Icges(5,5) * t512;
t396 = t193 * t349 + t195 * t346;
t31 = (qJD(2) * t396 + t114) * t347 + (qJD(2) * t191 - t116 * t349 - t118 * t346 + (t193 * t346 - t195 * t349) * qJD(4)) * t350;
t58 = t184 * t248 + t185 * t251 + t201 * t512 + t204 * t277 + t207 * t278 - t245 * t583;
t476 = t31 / 0.2e1 + t58 / 0.2e1;
t379 = t445 * t351;
t182 = t349 * t379 + (t349 * t484 - t380) * t348;
t183 = t348 * t381 + (t379 + t461) * t346;
t113 = Icges(5,5) * t183 + Icges(5,6) * t182 + Icges(5,3) * t364;
t115 = Icges(5,4) * t183 + Icges(5,2) * t182 + Icges(5,6) * t364;
t117 = Icges(5,1) * t183 + Icges(5,4) * t182 + Icges(5,5) * t364;
t192 = Icges(5,5) * t280 + Icges(5,6) * t279 + Icges(5,3) * t514;
t194 = Icges(5,4) * t280 + Icges(5,2) * t279 + Icges(5,6) * t514;
t196 = Icges(5,1) * t280 + Icges(5,4) * t279 + Icges(5,5) * t514;
t395 = t194 * t349 + t196 * t346;
t32 = (qJD(2) * t395 + t113) * t347 + (qJD(2) * t192 - t115 * t349 - t117 * t346 + (t194 * t346 - t196 * t349) * qJD(4)) * t350;
t57 = t182 * t248 + t183 * t251 + t201 * t514 + t204 * t279 + t207 * t280 + t245 * t364;
t475 = t32 / 0.2e1 + t57 / 0.2e1;
t124 = t245 * t512 + t248 * t277 + t251 * t278;
t95 = t191 * t347 - t350 * t396;
t472 = -t95 / 0.2e1 - t124 / 0.2e1;
t125 = t245 * t514 + t248 * t279 + t251 * t280;
t96 = t192 * t347 - t350 * t395;
t471 = t96 / 0.2e1 + t125 / 0.2e1;
t318 = pkin(2) * t464;
t470 = t348 * (pkin(2) * t465 + t348 * t482 - t318 + (t347 * t487 + t461) * qJ(3)) + t351 * (-pkin(2) * t583 - qJ(3) * t467 + t496) + t283 * t487;
t197 = t278 * rSges(5,1) + t277 * rSges(5,2) + rSges(5,3) * t512;
t331 = pkin(6) * t487;
t469 = t331 + t496;
t430 = rSges(5,1) * t346 + rSges(5,2) * t349;
t262 = -t350 * t430 + t548;
t468 = t262 * t489;
t458 = t514 / 0.2e1;
t457 = t512 / 0.2e1;
t456 = t405 * t580 - t408 * qJD(2) / 0.2e1;
t453 = -qJ(3) - t554;
t452 = -t347 * pkin(7) - t307;
t88 = t191 * t514 + t193 * t279 + t195 * t280;
t89 = t192 * t514 + t194 * t279 + t196 * t280;
t418 = t348 * t89 + t351 * t88;
t44 = t347 * t125 + t350 * t418;
t451 = t347 * t96 + t44;
t235 = t497 * t351;
t448 = t347 * t104 + t178 * t484 + t583 * t242;
t100 = Icges(6,4) * t159 + Icges(6,2) * t158 - Icges(6,6) * t583;
t102 = Icges(6,1) * t159 + Icges(6,4) * t158 - Icges(6,5) * t583;
t172 = Icges(6,5) * t267 + Icges(6,6) * t266 + Icges(6,3) * t512;
t174 = Icges(6,4) * t267 + Icges(6,2) * t266 + Icges(6,6) * t512;
t176 = Icges(6,1) * t267 + Icges(6,4) * t266 + Icges(6,5) * t512;
t398 = t174 * t334 + t176 * t333;
t98 = Icges(6,5) * t159 + Icges(6,6) * t158 - Icges(6,3) * t583;
t25 = (qJD(2) * t398 + t98) * t347 + (qJD(2) * t172 + (-t176 * t342 - t100) * t334 + (t174 * t342 - t102) * t333) * t350;
t101 = Icges(6,1) * t157 + Icges(6,4) * t156 + Icges(6,5) * t364;
t173 = Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t514;
t175 = Icges(6,4) * t269 + Icges(6,2) * t268 + Icges(6,6) * t514;
t177 = Icges(6,1) * t269 + Icges(6,4) * t268 + Icges(6,5) * t514;
t397 = t175 * t334 + t177 * t333;
t97 = Icges(6,5) * t157 + Icges(6,6) * t156 + Icges(6,3) * t364;
t99 = Icges(6,4) * t157 + Icges(6,2) * t156 + Icges(6,6) * t364;
t26 = (qJD(2) * t397 + t97) * t347 + (qJD(2) * t173 + (-t177 * t342 - t99) * t334 + (t175 * t342 - t101) * t333) * t350;
t112 = t239 * t514 + t240 * t268 + t241 * t269;
t74 = t172 * t514 + t174 * t268 + t176 * t269;
t75 = t173 * t514 + t175 * t268 + t177 * t269;
t422 = t348 * t75 + t351 * t74;
t40 = t112 * t347 + t350 * t422;
t84 = t172 * t347 - t350 * t398;
t85 = t173 * t347 - t350 * t397;
t420 = t348 * t84 - t351 * t85;
t421 = t348 * t85 + t351 * t84;
t19 = t100 * t268 + t102 * t269 + t156 * t174 + t157 * t176 + t172 * t364 + t514 * t98;
t20 = t101 * t269 + t156 * t175 + t157 * t177 + t173 * t364 + t268 * t99 + t514 * t97;
t47 = t156 * t240 + t157 * t241 + t165 * t514 + t166 * t268 + t167 * t269 + t239 * t364;
t55 = t348 * t74 - t351 * t75;
t5 = (-qJD(2) * t422 + t47) * t347 + (-qJD(1) * t55 + qJD(2) * t112 + t19 * t351 + t20 * t348) * t350;
t111 = t239 * t512 + t240 * t266 + t241 * t267;
t21 = t100 * t266 + t102 * t267 + t158 * t174 + t159 * t176 - t172 * t583 + t512 * t98;
t22 = t101 * t267 + t158 * t175 + t159 * t177 - t173 * t583 + t266 * t99 + t512 * t97;
t72 = t172 * t512 + t174 * t266 + t176 * t267;
t73 = t173 * t512 + t175 * t266 + t177 * t267;
t423 = t348 * t73 + t351 * t72;
t48 = t158 * t240 + t159 * t241 + t165 * t512 + t166 * t266 + t167 * t267 - t239 * t583;
t54 = t348 * t72 - t351 * t73;
t6 = (-qJD(2) * t423 + t48) * t347 + (-qJD(1) * t54 + qJD(2) * t111 + t21 * t351 + t22 * t348) * t350;
t447 = t5 * t514 + t6 * t512 + t40 * t465 + (t132 * t347 + t350 * t421) * t484 + t347 * (-t421 * t486 + (-qJD(1) * t420 + t25 * t351 + t26 * t348) * t350 + t544);
t295 = pkin(7) * t514 - t340;
t444 = t351 * t294 + t348 * t295 + t501;
t443 = rSges(4,1) * t487 + t583 * rSges(4,2) + rSges(4,3) * t460;
t442 = t494 + t284;
t39 = t111 * t347 + t350 * t423;
t86 = t191 * t512 + t193 * t277 + t195 * t278;
t87 = t192 * t512 + t194 * t277 + t196 * t278;
t419 = t348 * t87 + t351 * t86;
t43 = t347 * t124 + t350 * t419;
t441 = -t347 * t95 - t39 - t43;
t436 = -t262 + t452;
t435 = -pkin(7) * t484 - t273;
t434 = t347 * t546 - pkin(1);
t433 = rSges(3,1) * t350 - rSges(3,2) * t347;
t432 = rSges(5,1) * t183 + rSges(5,2) * t182;
t60 = t348 * t86 - t351 * t87;
t61 = t348 * t88 - t351 * t89;
t394 = -t197 * t351 - t198 * t348;
t393 = t197 * t348 - t524;
t382 = t452 - t502;
t263 = rSges(3,1) * t512 + t573;
t264 = rSges(4,3) * t516 + t574;
t216 = (-rSges(5,1) * t349 + rSges(5,2) * t346) * t481 + (rSges(5,3) * t350 + t347 * t430) * qJD(2);
t378 = -t216 + t435;
t377 = -pkin(1) - t433;
t190 = t436 * t351;
t375 = t348 * (t294 * qJD(1) - t317) + t351 * (-pkin(7) * t583 + t332) + t295 * t487 + t470;
t374 = t347 * t453 - pkin(1);
t373 = t435 + t506;
t372 = qJD(2) * t309;
t369 = qJD(2) * (Icges(4,4) * t347 + Icges(4,5) * t350);
t368 = qJD(2) * (-Icges(3,5) * t347 - Icges(3,6) * t350);
t151 = t382 * t351;
t362 = t350 * t480 - pkin(1) - t527;
t361 = (rSges(4,2) - pkin(2)) * t350 + t434;
t360 = t362 * t348;
t12 = qJD(1) * t422 + t19 * t348 - t20 * t351;
t13 = qJD(1) * t423 + t21 * t348 - t22 * t351;
t356 = t12 * t458 + t13 * t457 + t5 * t559 + t6 * t560 + (qJD(1) * t421 + t25 * t348 - t26 * t351) * t561 + t40 * t489 / 0.2e1 + t39 * t455 + t420 * t484 / 0.2e1 + t584 * t55 + t585 * t54;
t355 = (-pkin(2) - t545) * t350 + t374;
t354 = t544 + (t26 + t47) * t458 + (t25 + t48) * t457 + (t111 + t84) * t585 + (t112 + t85) * t584;
t353 = (-t348 * t40 - t351 * t39) * t486 - t39 * t466 + t447;
t339 = t351 * pkin(6);
t293 = t433 * qJD(2);
t265 = t348 * t426 - t551;
t261 = t348 * t433 - t549;
t234 = t497 * t348;
t230 = t263 + t494;
t229 = t348 * t377 + t339 + t549;
t226 = t242 * t514;
t215 = t570 * qJD(1) + t351 * t369;
t214 = t348 * t369 + t491;
t203 = t348 * t368 + t492;
t202 = -qJD(1) * t246 + t351 * t368;
t200 = t264 + t442;
t199 = t348 * t361 + t339 + t551;
t189 = t436 * t348;
t171 = t309 * t485 + ((-rSges(3,3) - pkin(6)) * t348 + t377 * t351) * qJD(1);
t170 = -rSges(3,1) * t583 - rSges(3,2) * t460 - pkin(1) * t489 + t331 + t495;
t168 = t347 * t178;
t163 = t179 * t512;
t153 = qJD(1) * t235 + t348 * t500;
t152 = t351 * t500 - t425 * t489 + t285;
t150 = t382 * t348;
t149 = t197 * t347 - t262 * t512;
t148 = -t198 * t347 + t262 * t514;
t147 = t348 * t385 + t351 * t570;
t146 = -t258 * t351 + t579;
t145 = -t348 * t570 + t577;
t144 = t258 * t348 + t386 * t351;
t143 = t247 * t348 - t387 * t351;
t142 = t246 * t348 - t576;
t140 = -t247 * t351 - t578;
t139 = -t246 * t351 - t348 * t388;
t138 = t442 + t197 + t294;
t137 = t339 + t340 + t360 + t431;
t135 = t264 * t351 + t265 * t348 + t501;
t134 = -t242 * t512 + t168;
t133 = -t179 * t347 + t226;
t131 = t318 + (-t482 + (t350 * t546 - t550) * qJD(2)) * t348 + ((-rSges(4,1) - pkin(6)) * t348 + t361 * t351) * qJD(1);
t130 = -pkin(2) * t463 + (t434 - t555) * t489 + t443 + t469;
t126 = t393 * t350;
t123 = t376 + t442 + t178;
t122 = t348 * t355 + t339 + t428 + t511;
t121 = -t178 * t514 + t163;
t120 = -rSges(5,3) * t583 + t503;
t119 = rSges(5,3) * t364 + t432;
t110 = qJD(1) * t190 + t348 * t378;
t109 = t351 * t378 + t468 + t498;
t92 = t217 * t347 - t502 * t512 + t168;
t91 = t272 * t514 - t347 * t504 + t226;
t90 = -t394 + t444;
t79 = t317 + t318 + (-t482 + (-t526 + t548) * qJD(2)) * t348 + ((-pkin(3) - pkin(6)) * t348 + t362 * t351) * qJD(1) - t432;
t78 = qJD(1) * t360 + t463 * t480 + t332 + t469 + t503;
t77 = qJD(1) * t151 + t348 * t373;
t76 = t351 * t373 + t489 * t502 + t498;
t71 = t163 + (t218 * t351 - t348 * t505) * t350;
t70 = t444 - t569;
t69 = (-t262 * t485 - t119) * t347 + (-qJD(2) * t198 + t216 * t348 + t262 * t487) * t350;
t68 = (t262 * t483 + t120) * t347 + (qJD(2) * t197 - t216 * t351 + t468) * t350;
t67 = -t351 * t477 + t318 + ((-t349 * t547 - qJD(3)) * t347 + (t347 * t545 + t350 * t453) * qJD(2)) * t348 + ((-pkin(6) - t328) * t348 + t355 * t351) * qJD(1) - t429;
t66 = t557 * t463 + (-t477 + (t350 * t557 + t374) * qJD(1)) * t348 + t415 + t469 + t508;
t65 = (qJD(1) * t265 + t443) * t351 + (t425 * t485 + (-t264 - t284 + t574) * qJD(1)) * t348 + t470;
t63 = -t103 * t347 + (-t179 * t350 - t242 * t517) * qJD(2) + t507;
t62 = -t169 * t512 + t448;
t49 = t393 * t486 + (qJD(1) * t394 + t119 * t351 - t120 * t348) * t350;
t41 = -t179 * t463 + (-t104 * t348 + (-t178 * t351 - t179 * t348) * qJD(1)) * t350 + t542;
t35 = t119 * t348 + t120 * t351 + (t524 + (-t197 + t499) * t348) * qJD(1) + t375;
t34 = (t225 * t348 + t272 * t487) * t350 - t510 * t347 + (-t350 * t504 - t502 * t517) * qJD(2) + t507;
t33 = (t272 * t483 + t129) * t347 + (qJD(2) * t217 + t272 * t489 + t351 * t506) * t350 + t448;
t30 = t113 * t512 + t115 * t277 + t117 * t278 + t184 * t194 + t185 * t196 - t192 * t583;
t29 = t114 * t512 + t116 * t277 + t118 * t278 + t184 * t193 + t185 * t195 - t191 * t583;
t28 = t113 * t514 + t115 * t279 + t117 * t280 + t182 * t194 + t183 * t196 + t192 * t364;
t27 = t114 * t514 + t116 * t279 + t118 * t280 + t182 * t193 + t183 * t195 + t191 * t364;
t18 = t509 * t351 + t510 * t348 + (t575 + (t499 - t505) * t348) * qJD(1) + t375;
t17 = (t217 * t348 - t575) * t486 + (t569 * qJD(1) + t128 * t351 - t509 * t348) * t350 + t542;
t16 = qJD(1) * t419 + t29 * t348 - t30 * t351;
t15 = qJD(1) * t418 + t27 * t348 - t28 * t351;
t9 = (-qJD(2) * t419 + t58) * t347 + (-t60 * qJD(1) + qJD(2) * t124 + t29 * t351 + t30 * t348) * t350;
t8 = (-qJD(2) * t418 + t57) * t347 + (-qJD(1) * t61 + qJD(2) * t125 + t27 * t351 + t28 * t348) * t350;
t1 = [-t251 * t459 - t334 * t241 * t519 + (t137 * t79 + t138 * t78) * t566 + (t122 * t67 + t123 * t66) * t565 + (t130 * t200 + t131 * t199) * t567 + (t170 * t230 + t171 * t229) * t568 + t416 - t207 * t518 + t417 + (-t166 * t334 - t204 * t349 - t525) * t350 + (t414 + t402 + t593) * t486 + (t410 + t400 - t592) * t484; m(5) * (t109 * t137 + t110 * t138 + t189 * t78 + t190 * t79) + m(6) * (t122 * t76 + t123 * t77 + t150 * t66 + t151 * t67) + m(4) * (t130 * t234 + t131 * t235 + t152 * t199 + t153 * t200) + (-t26 / 0.2e1 - t47 / 0.2e1 + m(3) * (-t171 * t309 - t229 * t293) + t456 * t351 - t475) * t351 + (t25 / 0.2e1 + t48 / 0.2e1 + m(3) * (-t170 * t309 - t230 * t293) + t456 * t348 + t476) * t348 + ((-t588 * qJD(2) + t590 * t351) * t560 + (t589 * qJD(2) + t590 * t348) * t559 + (t586 * t559 - t587 * t560) * qJD(1)) * t347 + ((t586 * qJD(2) + t591 * t351) * t560 + (t587 * qJD(2) + t591 * t348) * t559 + (t588 * t559 + t589 * t560) * qJD(1)) * t350 + ((t84 / 0.2e1 + t111 / 0.2e1 - t230 * t556 + (t250 / 0.2e1 - t254 / 0.2e1) * t350 + (t253 / 0.2e1 - t256 / 0.2e1) * t347 - t472) * t351 + (t85 / 0.2e1 + t112 / 0.2e1 + t229 * t556 + (t249 / 0.2e1 + t572 / 0.2e1) * t350 + (t252 / 0.2e1 + t571 / 0.2e1) * t347 + t471) * t348) * qJD(1); t348 * t13 - t351 * t12 + (t150 * t77 + t151 * t76 + t70 * t18) * t565 + t348 * t16 - t351 * t15 + (t109 * t190 + t110 * t189 + t35 * t90) * t566 + (t135 * t65 + t152 * t235 + t153 * t234) * t567 - t351 * ((t214 * t351 + (t146 - t577) * qJD(1)) * t351 + (t147 * qJD(1) + (t254 * t484 + t256 * t486 + t491) * t348 + (-t215 + (t347 * t571 + t350 * t572) * qJD(2) + t386 * qJD(1)) * t351) * t348) + ((t261 * t348 + t263 * t351) * ((qJD(1) * t261 - t351 * t372 + t495) * t351 + (-t348 * t372 + (-t263 + t573) * qJD(1)) * t348) + t493 * t309 * t293) * t568 + t348 * ((t215 * t348 + (t145 - t579) * qJD(1)) * t348 + (t144 * qJD(1) + (t484 * t572 + t486 * t571) * t351 + (-t214 + (t254 * t350 + t256 * t347) * qJD(2) + (t258 + t385) * qJD(1)) * t348) * t351) - t351 * ((t351 * t203 + (t140 + t576) * qJD(1)) * t351 + (t139 * qJD(1) + (-t250 * t484 - t253 * t486 + t492) * t348 + (-t202 + (t249 * t350 + t252 * t347) * qJD(2) - t387 * qJD(1)) * t351) * t348) + t348 * ((t348 * t202 + (t142 + t578) * qJD(1)) * t348 + (t143 * qJD(1) + (t249 * t484 + t252 * t486) * t351 + (-t203 + (-t250 * t350 - t253 * t347) * qJD(2) + (t247 - t388) * qJD(1)) * t348) * t351) + (t55 + t61 + (-t139 - t147) * t351 + (t140 + t146) * t348) * t489 + (t54 + t60 + (-t142 - t145) * t351 + (t143 + t144) * t348) * t487; 0.2e1 * ((t122 * t351 + t123 * t348) * t562 + (t137 * t351 + t138 * t348) * t563 + (t199 * t351 + t200 * t348) * t564) * t484 + 0.2e1 * ((-t122 * t489 + t123 * t487 + t348 * t66 + t351 * t67) * t562 + (-t137 * t489 + t138 * t487 + t348 * t78 + t351 * t79) * t563 + (t130 * t348 + t131 * t351 - t199 * t489 + t200 * t487) * t564) * t347; 0.2e1 * ((t150 * t485 + t151 * t483 - t18) * t562 + (t189 * t485 + t190 * t483 - t35) * t563 + (t234 * t485 + t235 * t483 - t65) * t564) * t350 + 0.2e1 * ((qJD(2) * t70 + t150 * t487 - t151 * t489 + t348 * t77 + t351 * t76) * t562 + (qJD(2) * t90 + t109 * t351 + t110 * t348 + t189 * t487 - t190 * t489) * t563 + (qJD(2) * t135 + t152 * t351 + t153 * t348 + t234 * t487 - t235 * t489) * t564) * t347; 0.4e1 * (t564 + t563 + t562) * (-0.1e1 + t493) * t347 * t484; (-t348 * t471 + t351 * t472) * t486 + (t476 * t351 + t475 * t348 + (t348 * t472 + t351 * t471) * qJD(1)) * t350 + m(5) * (t137 * t69 + t138 * t68 + t148 * t79 + t149 * t78) + m(6) * (t122 * t34 + t123 * t33 + t66 * t92 + t67 * t91) + t354 + t543; t356 + (t15 * t560 + t16 * t558 + (t348 * t95 - t351 * t96) * t580 + (t61 * t558 - t348 * t60 / 0.2e1) * qJD(1)) * t350 + (t60 * t454 + (qJD(1) * t95 - t32) * t561 - t8 / 0.2e1 + t43 * t553) * t351 + ((qJD(1) * t96 + t31) * t561 + t61 * t454 + t9 / 0.2e1 + t44 * t553) * t348 + m(5) * (t109 * t148 + t110 * t149 - t126 * t35 + t189 * t68 + t190 * t69 + t49 * t90) + m(6) * (t150 * t33 + t151 * t34 + t17 * t70 + t71 * t18 + t76 * t91 + t77 * t92); 0.2e1 * ((t148 * t483 + t149 * t485 - t49) * t563 + (t483 * t91 + t485 * t92 - t17) * t562) * t350 + 0.2e1 * ((-qJD(2) * t126 - t148 * t489 + t149 * t487 + t348 * t68 + t351 * t69) * t563 + (qJD(2) * t71 + t33 * t348 + t34 * t351 + t487 * t92 - t489 * t91) * t562) * t347; (t71 * t17 + t33 * t92 + t34 * t91) * t565 + (-t126 * t49 + t148 * t69 + t149 * t68) * t566 + ((t441 * t351 + (-t40 - t451) * t348) * qJD(2) + t543) * t347 + (t351 * t9 + t348 * t8 + t347 * (t31 * t351 + t32 * t348) + (t141 * t347 + (t348 * t96 + t351 * t95) * t350) * qJD(2) + (t348 * t441 + t351 * t451) * qJD(1)) * t350 + t447; m(6) * (t122 * t63 + t123 * t62 + t133 * t67 + t134 * t66) + t354; t356 + m(6) * (t121 * t18 + t133 * t76 + t134 * t77 + t150 * t62 + t151 * t63 + t41 * t70); m(6) * ((-t41 + (t133 * t351 + t134 * t348) * qJD(2)) * t350 + (qJD(2) * t121 + t348 * t62 + t351 * t63 + (-t133 * t348 + t134 * t351) * qJD(1)) * t347); m(6) * (t121 * t17 + t133 * t34 + t134 * t33 + t41 * t71 + t62 * t92 + t63 * t91) + t353; (t121 * t41 + t133 * t63 + t134 * t62) * t565 + t353;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
