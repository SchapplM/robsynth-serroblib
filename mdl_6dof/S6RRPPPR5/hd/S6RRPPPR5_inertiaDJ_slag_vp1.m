% Calculate time derivative of joint inertia matrix for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPPR5_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPPR5_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:22
% EndTime: 2019-03-09 08:22:05
% DurationCPUTime: 30.86s
% Computational Cost: add. (14215->989), mult. (39174->1322), div. (0->0), fcn. (40416->8), ass. (0->418)
t554 = Icges(5,4) + Icges(6,6);
t562 = Icges(4,5) - t554;
t553 = Icges(6,4) + Icges(4,6);
t561 = Icges(5,5) - t553;
t335 = sin(qJ(1));
t337 = cos(qJ(2));
t493 = t335 * t337;
t334 = sin(qJ(2));
t331 = sin(pkin(9));
t332 = cos(pkin(9));
t381 = -Icges(5,6) * t332 + Icges(5,3) * t331;
t388 = Icges(4,4) * t332 - Icges(4,2) * t331;
t391 = Icges(6,1) * t331 + Icges(6,5) * t332;
t560 = ((t381 - t388 + t391) * t337 + t561 * t334) * qJD(2);
t382 = -Icges(5,2) * t332 + Icges(5,6) * t331;
t383 = Icges(6,5) * t331 + Icges(6,3) * t332;
t392 = Icges(4,1) * t332 - Icges(4,4) * t331;
t559 = ((-t382 + t383 + t392) * t337 + t562 * t334) * qJD(2);
t240 = -Icges(4,6) * t337 + t334 * t388;
t241 = Icges(6,4) * t337 + t334 * t391;
t243 = -Icges(5,5) * t337 + t334 * t381;
t558 = t243 - t240 + t241;
t237 = Icges(6,6) * t337 + t334 * t383;
t242 = -Icges(4,5) * t337 + t334 * t392;
t244 = -Icges(5,4) * t337 + t334 * t382;
t557 = t244 - t242 - t237;
t556 = Icges(5,5) / 0.2e1 - t553 / 0.2e1;
t555 = Icges(4,5) / 0.2e1 - t554 / 0.2e1;
t525 = t335 / 0.2e1;
t338 = cos(qJ(1));
t522 = t338 / 0.2e1;
t552 = Icges(5,1) + Icges(6,2) + Icges(4,3);
t543 = -qJD(1) / 0.2e1;
t526 = m(7) / 0.2e1;
t527 = m(6) / 0.2e1;
t464 = t527 + t526;
t542 = 0.2e1 * t464;
t333 = sin(qJ(6));
t336 = cos(qJ(6));
t375 = t331 * t336 + t332 * t333;
t260 = t375 * t334;
t503 = Icges(3,4) * t337;
t390 = -Icges(3,2) * t334 + t503;
t253 = Icges(3,6) * t335 + t338 * t390;
t504 = Icges(3,4) * t334;
t394 = Icges(3,1) * t337 - t504;
t255 = Icges(3,5) * t335 + t338 * t394;
t376 = t253 * t334 - t255 * t337;
t360 = t376 * t335;
t252 = -Icges(3,6) * t338 + t335 * t390;
t254 = -Icges(3,5) * t338 + t335 * t394;
t377 = t252 * t334 - t254 * t337;
t361 = t377 * t338;
t494 = t334 * t338;
t541 = -rSges(3,2) * t494 + t335 * rSges(3,3);
t277 = t331 * t493 + t332 * t338;
t278 = -t331 * t338 + t332 * t493;
t495 = t334 * t335;
t154 = Icges(5,5) * t495 - Icges(5,6) * t278 + Icges(5,3) * t277;
t166 = Icges(4,4) * t278 - Icges(4,2) * t277 + Icges(4,6) * t495;
t168 = Icges(6,1) * t277 - Icges(6,4) * t495 + Icges(6,5) * t278;
t539 = t154 - t166 + t168;
t492 = t337 * t338;
t279 = t331 * t492 - t335 * t332;
t280 = t335 * t331 + t332 * t492;
t155 = Icges(5,5) * t494 - Icges(5,6) * t280 + Icges(5,3) * t279;
t167 = Icges(4,4) * t280 - Icges(4,2) * t279 + Icges(4,6) * t494;
t169 = Icges(6,1) * t279 - Icges(6,4) * t494 + Icges(6,5) * t280;
t538 = t155 - t167 + t169;
t156 = Icges(5,4) * t495 - Icges(5,2) * t278 + Icges(5,6) * t277;
t160 = Icges(6,5) * t277 - Icges(6,6) * t495 + Icges(6,3) * t278;
t170 = Icges(4,1) * t278 - Icges(4,4) * t277 + Icges(4,5) * t495;
t537 = -t156 + t160 + t170;
t157 = Icges(5,4) * t494 - Icges(5,2) * t280 + Icges(5,6) * t279;
t161 = Icges(6,5) * t279 - Icges(6,6) * t494 + Icges(6,3) * t280;
t171 = Icges(4,1) * t280 - Icges(4,4) * t279 + Icges(4,5) * t494;
t536 = -t157 + t161 + t171;
t535 = t561 * t331 + t332 * t562;
t385 = Icges(3,5) * t337 - Icges(3,6) * t334;
t250 = -Icges(3,3) * t338 + t335 * t385;
t534 = 2 * m(3);
t533 = 2 * m(4);
t532 = 2 * m(5);
t531 = 0.2e1 * m(6);
t530 = 0.2e1 * m(7);
t328 = t335 ^ 2;
t330 = t338 ^ 2;
t529 = m(4) / 0.2e1;
t528 = m(5) / 0.2e1;
t524 = -t337 / 0.2e1;
t521 = rSges(5,2) - pkin(3);
t520 = -rSges(6,3) - pkin(3);
t519 = pkin(2) * t334;
t518 = pkin(2) * t337;
t517 = pkin(5) * t277;
t468 = qJD(2) * t335;
t428 = t334 * t468;
t470 = qJD(1) * t338;
t431 = t331 * t470;
t471 = qJD(1) * t335;
t212 = -t331 * t428 - t332 * t471 + t337 * t431;
t516 = t212 * pkin(5);
t213 = qJD(1) * t280 - t332 * t428;
t515 = t213 * pkin(3);
t324 = t335 * pkin(7);
t514 = qJD(1) / 0.2e1;
t513 = -pkin(4) - qJ(3);
t512 = rSges(6,1) * t277;
t511 = rSges(6,2) * t337;
t510 = rSges(3,3) * t338;
t509 = rSges(5,3) * t277;
t508 = t212 * rSges(6,1);
t507 = t212 * rSges(5,3);
t506 = -rSges(5,1) - qJ(3);
t505 = -rSges(4,3) - qJ(3);
t499 = t331 * t334;
t498 = t331 * t337;
t497 = t332 * t334;
t496 = t332 * t337;
t191 = -t277 * t333 + t278 * t336;
t192 = t277 * t336 + t278 * t333;
t407 = -rSges(7,1) * t192 - rSges(7,2) * t191;
t107 = rSges(7,3) * t495 - t407;
t491 = pkin(8) * t495 + t107 + t517;
t193 = -t279 * t333 + t280 * t336;
t194 = t279 * t336 + t280 * t333;
t108 = t194 * rSges(7,1) + t193 * rSges(7,2) + rSges(7,3) * t494;
t490 = t279 * pkin(5) + pkin(8) * t494 + t108;
t374 = -t331 * t333 + t332 * t336;
t259 = t374 * t334;
t151 = rSges(7,1) * t260 + rSges(7,2) * t259 - rSges(7,3) * t337;
t489 = pkin(5) * t499 - pkin(8) * t337 + t151;
t198 = t280 * pkin(3) + t279 * qJ(4);
t318 = pkin(2) * t492;
t284 = qJ(3) * t494 + t318;
t488 = -t198 - t284;
t487 = -t212 * qJ(4) - t277 * qJD(4);
t486 = t213 * qJ(5) + t278 * qJD(5);
t405 = qJ(3) * t334 + t518;
t276 = qJD(2) * t405 - qJD(3) * t337;
t404 = pkin(3) * t332 + qJ(4) * t331;
t467 = qJD(2) * t337;
t485 = -qJD(4) * t499 - t404 * t467 - t276;
t409 = rSges(4,1) * t332 - rSges(4,2) * t331;
t484 = -(rSges(4,3) * t334 + t337 * t409) * qJD(2) - t276;
t281 = t404 * t334;
t297 = -qJ(3) * t337 + t519;
t285 = t297 * t471;
t483 = t281 * t471 + t285;
t248 = -rSges(4,3) * t337 + t334 * t409;
t482 = -t248 - t297;
t283 = t405 * t335;
t481 = t335 * t283 + t338 * t284;
t480 = -t281 - t297;
t466 = qJD(2) * t338;
t424 = t337 * t466;
t465 = qJD(3) * t334;
t479 = qJ(3) * t424 + t338 * t465;
t432 = t334 * t471;
t478 = rSges(3,2) * t432 + rSges(3,3) * t470;
t426 = t334 * t466;
t210 = qJD(1) * t277 + t331 * t426;
t477 = -t210 * pkin(5) + pkin(8) * t424;
t476 = t338 * pkin(1) + t324;
t265 = t277 * qJ(4);
t325 = t338 * pkin(7);
t475 = t325 - t265;
t474 = t334 ^ 2 - t337 ^ 2;
t473 = t328 + t330;
t251 = Icges(3,3) * t335 + t338 * t385;
t472 = qJD(1) * t251;
t469 = qJD(2) * t334;
t463 = rSges(6,2) + t513;
t354 = t337 * t471 + t426;
t211 = t332 * t354 - t431;
t89 = -qJD(6) * t194 + t210 * t333 - t211 * t336;
t90 = qJD(6) * t193 - t210 * t336 - t211 * t333;
t462 = t90 * rSges(7,1) + t89 * rSges(7,2) + rSges(7,3) * t424;
t460 = rSges(6,2) * t494;
t101 = Icges(7,5) * t192 + Icges(7,6) * t191 + Icges(7,3) * t495;
t103 = Icges(7,4) * t192 + Icges(7,2) * t191 + Icges(7,6) * t495;
t105 = Icges(7,1) * t192 + Icges(7,4) * t191 + Icges(7,5) * t495;
t32 = -t101 * t337 + t103 * t259 + t105 * t260;
t148 = Icges(7,5) * t260 + Icges(7,6) * t259 - Icges(7,3) * t337;
t149 = Icges(7,4) * t260 + Icges(7,2) * t259 - Icges(7,6) * t337;
t150 = Icges(7,1) * t260 + Icges(7,4) * t259 - Icges(7,5) * t337;
t47 = t148 * t495 + t149 * t191 + t150 * t192;
t459 = t32 / 0.2e1 + t47 / 0.2e1;
t102 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t494;
t104 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t494;
t106 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t494;
t33 = -t102 * t337 + t104 * t259 + t106 * t260;
t48 = t148 * t494 + t193 * t149 + t194 * t150;
t458 = t33 / 0.2e1 + t48 / 0.2e1;
t158 = Icges(5,1) * t495 - Icges(5,4) * t278 + Icges(5,5) * t277;
t457 = t158 * t495;
t456 = t158 * t494;
t159 = Icges(5,1) * t494 - Icges(5,4) * t280 + Icges(5,5) * t279;
t455 = t159 * t495;
t454 = t159 * t494;
t162 = Icges(4,5) * t278 - Icges(4,6) * t277 + Icges(4,3) * t495;
t453 = t162 * t495;
t452 = t162 * t494;
t163 = Icges(4,5) * t280 - Icges(4,6) * t279 + Icges(4,3) * t494;
t451 = t163 * t495;
t450 = t163 * t494;
t164 = Icges(6,4) * t277 - Icges(6,2) * t495 + Icges(6,6) * t278;
t449 = t164 * t495;
t448 = t164 * t494;
t165 = Icges(6,4) * t279 - Icges(6,2) * t494 + Icges(6,6) * t280;
t447 = t165 * t495;
t446 = t165 * t494;
t305 = pkin(2) * t428;
t425 = t335 * t467;
t353 = t334 * t470 + t425;
t445 = t335 * (qJ(3) * t353 + qJD(1) * t318 + t335 * t465 - t305) + t338 * (-pkin(2) * t354 - qJ(3) * t432 + t479) + t283 * t470;
t221 = pkin(4) * t494 + t280 * qJ(5);
t444 = -t221 + t488;
t443 = -t211 * rSges(4,1) + t210 * rSges(4,2) + rSges(4,3) * t424;
t442 = rSges(5,1) * t424 + t211 * rSges(5,2) - t210 * rSges(5,3);
t406 = -rSges(5,2) * t332 + rSges(5,3) * t331;
t441 = -(rSges(5,1) * t334 + t337 * t406) * qJD(2) + t485;
t440 = -qJD(5) * t497 - (pkin(4) * t334 + qJ(5) * t496) * qJD(2) + t485;
t289 = -pkin(4) * t337 + qJ(5) * t497;
t439 = t289 * t471 + t483;
t249 = -rSges(5,1) * t337 + t334 * t406;
t438 = -t249 + t480;
t177 = t280 * rSges(4,1) - t279 * rSges(4,2) + rSges(4,3) * t494;
t437 = -t289 + t480;
t322 = pkin(7) * t470;
t436 = t322 + t479;
t435 = -t210 * rSges(6,1) + rSges(6,2) * t432 - t211 * rSges(6,3);
t434 = t305 + t487;
t173 = rSges(5,1) * t494 - t280 * rSges(5,2) + t279 * rSges(5,3);
t433 = -pkin(1) - t518;
t427 = t334 * t467;
t423 = t467 / 0.2e1;
t422 = -rSges(7,3) - pkin(8) + t513;
t197 = pkin(3) * t278 + t265;
t196 = t482 * t338;
t421 = pkin(2) * t426;
t420 = t528 + t464;
t419 = t335 * t197 + t338 * t198 + t481;
t408 = rSges(6,1) * t331 + rSges(6,3) * t332;
t418 = -(-rSges(6,2) * t334 + t337 * t408) * qJD(2) + t440;
t247 = t334 * t408 + t511;
t417 = -t247 + t437;
t416 = t476 + t284;
t144 = t438 * t338;
t91 = -qJD(6) * t192 - t212 * t333 + t213 * t336;
t92 = qJD(6) * t191 + t212 * t336 + t213 * t333;
t415 = t92 * rSges(7,1) + t91 * rSges(7,2);
t414 = t237 / 0.2e1 + t242 / 0.2e1 - t244 / 0.2e1;
t413 = t241 / 0.2e1 - t240 / 0.2e1 + t243 / 0.2e1;
t412 = rSges(3,1) * t337 - rSges(3,2) * t334;
t298 = rSges(3,1) * t334 + rSges(3,2) * t337;
t411 = -t213 * rSges(4,1) + t212 * rSges(4,2);
t410 = -rSges(4,1) * t278 + rSges(4,2) * t277;
t28 = t101 * t495 + t103 * t191 + t105 * t192;
t29 = t102 * t495 + t104 * t191 + t106 * t192;
t403 = t28 * t338 - t29 * t335;
t402 = t28 * t335 + t29 * t338;
t30 = t101 * t494 + t193 * t103 + t194 * t105;
t31 = t102 * t494 + t193 * t104 + t194 * t106;
t401 = t30 * t338 - t31 * t335;
t400 = t30 * t335 + t31 * t338;
t399 = t32 * t338 - t33 * t335;
t398 = t32 * t335 + t33 * t338;
t178 = -qJD(6) * t260 + t374 * t467;
t179 = qJD(6) * t259 + t375 * t467;
t100 = rSges(7,1) * t179 + rSges(7,2) * t178 + rSges(7,3) * t469;
t397 = -t100 - (pkin(5) * t498 + pkin(8) * t334) * qJD(2) + t440;
t396 = t437 - t489;
t395 = t434 - t486;
t393 = Icges(3,1) * t334 + t503;
t389 = Icges(3,2) * t337 + t504;
t380 = t107 * t338 - t108 * t335;
t258 = rSges(3,1) * t492 + t541;
t176 = t279 * rSges(6,1) + t280 * rSges(6,3) - t460;
t138 = t417 * t338;
t373 = -pkin(1) - t412;
t266 = t278 * qJ(5);
t220 = pkin(4) * t495 + t266;
t370 = t335 * t220 + t338 * t221 + t419;
t26 = t335 * t491 + t338 * t490 + t370;
t363 = pkin(4) * t424 - t211 * qJ(5) + t280 * qJD(5);
t364 = -t211 * pkin(3) - t210 * qJ(4) + t279 * qJD(4);
t371 = t335 * (-t487 + t515) + t338 * t364 + t197 * t470 + t445;
t343 = t335 * (pkin(4) * t353 + t486) + t338 * (-pkin(4) * t432 + t363) + t220 * t470 + t371;
t45 = -rSges(7,3) * t432 + t462;
t46 = rSges(7,3) * t353 + t415;
t5 = (t45 + t477) * t338 + (pkin(8) * t425 + t46 + t516) * t335 + (t491 * t338 + (t444 - t490) * t335) * qJD(1) + t343;
t372 = t26 * t467 + t334 * t5;
t14 = t380 * t467 + (-t335 * t45 + t338 * t46 + (-t335 * t107 - t108 * t338) * qJD(1)) * t334;
t55 = t380 * t334;
t369 = t14 * t334 + t467 * t55;
t174 = -rSges(6,2) * t495 + rSges(6,3) * t278 + t512;
t15 = t335 * (-rSges(6,2) * t425 + t213 * rSges(6,3) + t508) + t338 * (-rSges(6,2) * t424 + t435) + (t338 * t174 + (-t176 + t444 - t460) * t335) * qJD(1) + t343;
t38 = t335 * t174 + t176 * t338 + t370;
t368 = t15 * t334 + t38 * t467;
t362 = qJD(2) * t298;
t86 = t396 * t338;
t359 = qJD(2) * t393;
t358 = qJD(2) * t389;
t357 = qJD(2) * (-Icges(3,5) * t334 - Icges(3,6) * t337);
t356 = t334 * t506 + t433;
t355 = t334 * t505 + t433;
t97 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t469;
t98 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t469;
t99 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t469;
t21 = t148 * t469 + t178 * t149 + t179 * t150 + t259 * t98 + t260 * t99 - t337 * t97;
t352 = -t424 + t432;
t351 = t198 + t416;
t350 = t334 * t463 + t433;
t346 = t334 * t422 + t433;
t345 = t356 * t335;
t344 = t355 * t335;
t342 = t364 + t436;
t341 = t346 * t335;
t340 = t221 + t351;
t339 = t342 + t363;
t290 = t412 * qJD(2);
t257 = t335 * t412 - t510;
t218 = t258 + t476;
t217 = t335 * t373 + t325 + t510;
t195 = t482 * t335;
t181 = t335 * t357 + t472;
t180 = -qJD(1) * t250 + t338 * t357;
t175 = rSges(4,3) * t495 - t410;
t172 = rSges(5,1) * t495 - rSges(5,2) * t278 + t509;
t146 = t298 * t468 + ((-rSges(3,3) - pkin(7)) * t335 + t373 * t338) * qJD(1);
t145 = -rSges(3,1) * t354 - rSges(3,2) * t424 - pkin(1) * t471 + t322 + t478;
t143 = t438 * t335;
t142 = t416 + t177;
t141 = t325 + t344 + t410;
t137 = t417 * t335;
t136 = t335 * t251 - t338 * t376;
t135 = t335 * t250 - t361;
t134 = -t251 * t338 - t360;
t133 = -t250 * t338 - t335 * t377;
t132 = Icges(4,1) * t213 - Icges(4,4) * t212 + Icges(4,5) * t353;
t131 = -Icges(4,1) * t211 + Icges(4,4) * t210 - Icges(4,5) * t352;
t130 = Icges(6,1) * t212 - Icges(6,4) * t353 + Icges(6,5) * t213;
t129 = -Icges(6,1) * t210 + Icges(6,4) * t352 - Icges(6,5) * t211;
t128 = Icges(4,4) * t213 - Icges(4,2) * t212 + Icges(4,6) * t353;
t127 = -Icges(4,4) * t211 + Icges(4,2) * t210 - Icges(4,6) * t352;
t122 = Icges(6,5) * t212 - Icges(6,6) * t353 + Icges(6,3) * t213;
t121 = -Icges(6,5) * t210 + Icges(6,6) * t352 - Icges(6,3) * t211;
t118 = Icges(5,4) * t353 - Icges(5,2) * t213 + Icges(5,6) * t212;
t117 = -Icges(5,4) * t352 + Icges(5,2) * t211 - Icges(5,6) * t210;
t116 = Icges(5,5) * t353 - Icges(5,6) * t213 + Icges(5,3) * t212;
t115 = -Icges(5,5) * t352 + Icges(5,6) * t211 - Icges(5,3) * t210;
t112 = qJD(1) * t196 + t335 * t484;
t111 = t248 * t471 + t338 * t484 + t285;
t96 = t351 + t173;
t95 = t278 * t521 + t345 + t475 - t509;
t85 = t396 * t335;
t82 = t176 + t340;
t81 = t278 * t520 + t335 * t350 - t266 + t475 - t512;
t80 = t335 * t175 + t177 * t338 + t481;
t79 = t305 + (t467 * t505 - t465) * t335 + (t338 * t355 - t324) * qJD(1) + t411;
t78 = qJD(1) * t344 - t421 + t436 + t443;
t77 = qJD(1) * t144 + t335 * t441;
t76 = t249 * t471 + t338 * t441 + t483;
t75 = -t337 * t108 - t151 * t494;
t74 = t107 * t337 + t151 * t495;
t73 = -t279 * t167 + t280 * t171 + t450;
t72 = -t279 * t166 + t280 * t170 + t452;
t71 = t280 * t161 + t279 * t169 - t446;
t70 = t280 * t160 + t279 * t168 - t448;
t69 = -t167 * t277 + t171 * t278 + t451;
t68 = -t166 * t277 + t170 * t278 + t453;
t67 = t161 * t278 + t169 * t277 - t447;
t66 = t160 * t278 + t168 * t277 - t449;
t65 = t279 * t155 - t280 * t157 + t454;
t64 = t279 * t154 - t280 * t156 + t456;
t63 = t155 * t277 - t157 * t278 + t455;
t62 = t154 * t277 - t156 * t278 + t457;
t59 = qJD(1) * t138 + t335 * t418;
t58 = t247 * t471 + t338 * t418 + t439;
t57 = -t148 * t337 + t149 * t259 + t150 * t260;
t56 = t57 * t469;
t54 = t340 + t490;
t53 = -t197 - t266 + t325 + t341 + t407 - t517;
t51 = t335 * t172 + t173 * t338 + t419;
t50 = -t507 + t521 * t213 + (t467 * t506 - t465) * t335 + (t338 * t356 - t324) * qJD(1) + t434;
t49 = qJD(1) * t345 + t342 - t421 + t442;
t44 = Icges(7,1) * t92 + Icges(7,4) * t91 + Icges(7,5) * t353;
t43 = Icges(7,1) * t90 + Icges(7,4) * t89 - Icges(7,5) * t352;
t42 = Icges(7,4) * t92 + Icges(7,2) * t91 + Icges(7,6) * t353;
t41 = Icges(7,4) * t90 + Icges(7,2) * t89 - Icges(7,6) * t352;
t40 = Icges(7,5) * t92 + Icges(7,6) * t91 + Icges(7,3) * t353;
t39 = Icges(7,5) * t90 + Icges(7,6) * t89 - Icges(7,3) * t352;
t37 = -t508 + t520 * t213 + (t463 * t467 - t465) * t335 + (t338 * t350 - t324) * qJD(1) + t395;
t36 = (-t511 - t519) * t466 + (t334 * t513 + t433) * t471 + t339 + t435;
t35 = qJD(1) * t86 + t335 * t397;
t34 = t338 * t397 + t471 * t489 + t439;
t27 = t335 * (rSges(4,3) * t425 - t411) + t338 * t443 + (t338 * t175 + (-t177 - t284) * t335) * qJD(1) + t445;
t25 = -t515 - t516 + (t422 * t467 - t465) * t335 + (t338 * t346 - t324) * qJD(1) + t395 - t415;
t24 = qJD(1) * t341 + t339 - t421 + t462 + t477;
t23 = (t151 * t468 + t46) * t337 + (-qJD(2) * t107 + t335 * t100 + t151 * t470) * t334;
t22 = (-t151 * t466 - t45) * t337 + (qJD(2) * t108 - t100 * t338 + t151 * t471) * t334;
t20 = t335 * (rSges(5,1) * t425 - t213 * rSges(5,2) + t507) + t338 * t442 + (t338 * t172 + (-t173 + t488) * t335) * qJD(1) + t371;
t17 = t148 * t353 + t91 * t149 + t92 * t150 + t191 * t98 + t192 * t99 + t495 * t97;
t16 = -t148 * t352 + t89 * t149 + t90 * t150 + t193 * t98 + t194 * t99 + t494 * t97;
t13 = t334 * t400 - t48 * t337;
t12 = t334 * t402 - t47 * t337;
t11 = t102 * t469 + t104 * t178 + t106 * t179 + t259 * t41 + t260 * t43 - t337 * t39;
t10 = t101 * t469 + t103 * t178 + t105 * t179 + t259 * t42 + t260 * t44 - t337 * t40;
t9 = t102 * t353 + t91 * t104 + t92 * t106 + t191 * t41 + t192 * t43 + t39 * t495;
t8 = t101 * t353 + t91 * t103 + t92 * t105 + t191 * t42 + t192 * t44 + t40 * t495;
t7 = -t102 * t352 + t89 * t104 + t90 * t106 + t193 * t41 + t194 * t43 + t39 * t494;
t6 = -t101 * t352 + t89 * t103 + t90 * t105 + t193 * t42 + t194 * t44 + t40 * t494;
t4 = qJD(1) * t402 + t9 * t335 - t338 * t8;
t3 = qJD(1) * t400 + t7 * t335 - t338 * t6;
t2 = (qJD(2) * t402 - t17) * t337 + (qJD(1) * t403 + qJD(2) * t47 + t335 * t8 + t338 * t9) * t334;
t1 = (qJD(2) * t400 - t16) * t337 + (qJD(1) * t401 + qJD(2) * t48 + t335 * t6 + t338 * t7) * t334;
t18 = [(t24 * t54 + t25 * t53) * t530 + (t36 * t82 + t37 * t81) * t531 + (t49 * t96 + t50 * t95) * t532 + (t141 * t79 + t142 * t78) * t533 + (t145 * t218 + t146 * t217) * t534 + t21 + t560 * t499 + t559 * t497 + (t334 * t535 - t337 * t552 - t389 + t394) * t469 + (t331 * t558 - t332 * t557 - t552 * t334 - t535 * t337 + t390 + t393) * t467; m(4) * (t111 * t141 + t112 * t142 + t195 * t78 + t196 * t79) + m(6) * (t137 * t36 + t138 * t37 + t58 * t81 + t59 * t82) + m(5) * (t143 * t49 + t144 * t50 + t76 * t95 + t77 * t96) + m(7) * (t24 * t85 + t25 * t86 + t34 * t53 + t35 * t54) + m(3) * ((-t145 * t335 - t146 * t338) * t298 + (-t217 * t338 - t218 * t335) * t290) + ((t556 * t212 + t555 * t213 + t253 * t543) * t338 + (t556 * t210 + t555 * t211 + t252 * t543) * t335 + t552 * t353 * t522) * t337 + ((t279 * t413 + t280 * t414 + t458) * t338 + (t277 * t413 + t278 * t414 + t459) * t335 + m(3) * (t217 * t335 - t218 * t338) * t298 + ((t253 / 0.2e1 + t165 / 0.2e1 - t163 / 0.2e1 - t159 / 0.2e1) * t338 + (t252 / 0.2e1 + t164 / 0.2e1 - t162 / 0.2e1 - t158 / 0.2e1) * t335) * t337 + (t331 * t538 + t332 * t536 + t255) * t334 * t522) * qJD(1) + ((t330 / 0.2e1 + t328 / 0.2e1) * t385 + t361 / 0.2e1 - t360 / 0.2e1) * qJD(2) + (t16 + t11 + (t358 * t338 + t352 * t552) * t337 + (t536 * t496 + t538 * t498) * qJD(2) + (-t338 * t359 + (-t117 + t121 + t131) * t332 + (t115 - t127 + t129) * t331 + (t159 + t163 - t165) * qJD(2) + (t331 * t539 + t332 * t537) * qJD(1)) * t334 + t559 * t280 + t560 * t279 + t557 * t211 - t558 * t210) * t525 - (t17 + t10 + t358 * t493 + (qJD(1) * t255 - t335 * t359 + (-t118 + t122 + t132) * t332 + (t116 - t128 + t130) * t331) * t334 + (t539 * t498 + t537 * t496 + (t158 + t162 - t164) * t334) * qJD(2) + t559 * t278 + t560 * t277 - t557 * t213 + t558 * t212) * t338 / 0.2e1; (t26 * t5 + t34 * t86 + t35 * t85) * t530 + (t137 * t59 + t138 * t58 + t15 * t38) * t531 + (t143 * t77 + t144 * t76 + t20 * t51) * t532 + (t111 * t196 + t112 * t195 + t27 * t80) * t533 + ((t335 * t257 + t258 * t338) * ((qJD(1) * t257 - t338 * t362 + t478) * t338 + (-t335 * t362 + (-t258 + t541) * qJD(1)) * t335) + t473 * t298 * t290) * t534 - t338 * ((t338 * t181 + (t134 + t361) * qJD(1)) * t338 + (t133 * qJD(1) + (-t253 * t467 - t255 * t469 + t472) * t335 + (-t180 + (t252 * t337 + t254 * t334) * qJD(2) - t376 * qJD(1)) * t338) * t335) - t338 * ((-t278 * t122 - t277 * t130 - t213 * t160 - t212 * t168 + (t67 + t448) * qJD(1)) * t338 + (t278 * t121 + t277 * t129 + t213 * t161 + t212 * t169 + (t66 - t446) * qJD(1)) * t335) + t335 * ((t280 * t121 + t279 * t129 - t211 * t161 - t210 * t169 + (t70 + t447) * qJD(1)) * t335 + (-t280 * t122 - t279 * t130 + t211 * t160 + t210 * t168 + (t71 - t449) * qJD(1)) * t338) - t338 * ((t277 * t128 - t278 * t132 + t212 * t166 - t213 * t170 + (t69 - t452) * qJD(1)) * t338 + (-t277 * t127 + t278 * t131 - t212 * t167 + t213 * t171 + (t68 + t450) * qJD(1)) * t335) + t335 * ((-t279 * t127 + t280 * t131 + t210 * t167 - t211 * t171 + (t72 - t451) * qJD(1)) * t335 + (t279 * t128 - t280 * t132 - t210 * t166 + t211 * t170 + (t73 + t453) * qJD(1)) * t338) - t338 * ((-t277 * t116 + t278 * t118 - t212 * t154 + t213 * t156 + (t63 - t456) * qJD(1)) * t338 + (t277 * t115 - t278 * t117 + t212 * t155 - t213 * t157 + (t62 + t454) * qJD(1)) * t335) + t335 * ((t279 * t115 - t280 * t117 - t210 * t155 + t211 * t157 + (t64 - t455) * qJD(1)) * t335 + (-t279 * t116 + t280 * t118 + t210 * t154 - t211 * t156 + (t65 + t457) * qJD(1)) * t338) - t338 * t4 + t335 * t3 + t335 * ((t335 * t180 + (t135 + t360) * qJD(1)) * t335 + (t136 * qJD(1) + (t252 * t467 + t254 * t469) * t338 + (-t181 + (-t253 * t337 - t255 * t334) * qJD(2) + (t251 - t377) * qJD(1)) * t335) * t338) + (-t403 + (-t133 - t62 - t66 - t68) * t338 + (t134 + t63 + t67 + t69) * t335) * t471 + (-t401 + (-t135 - t64 - t70 - t72) * t338 + (t136 + t65 + t71 + t73) * t335) * t470; 0.2e1 * ((t335 * t54 + t338 * t53) * t526 + (t335 * t96 + t338 * t95) * t528 + (t335 * t82 + t338 * t81) * t527 + (t141 * t338 + t142 * t335) * t529) * t467 + 0.2e1 * ((t24 * t335 + t25 * t338 + t470 * t54 - t471 * t53) * t526 + (t335 * t49 + t338 * t50 + t470 * t96 - t471 * t95) * t528 + (t335 * t36 + t338 * t37 + t470 * t82 - t471 * t81) * t527 + (-t141 * t471 + t142 * t470 + t335 * t78 + t338 * t79) * t529) * t334; 0.2e1 * ((t466 * t86 + t468 * t85 - t5) * t526 + (t137 * t468 + t138 * t466 - t15) * t527 + (t143 * t468 + t144 * t466 - t20) * t528 + (t195 * t468 + t196 * t466 - t27) * t529) * t337 + 0.2e1 * ((qJD(2) * t26 + t335 * t35 + t338 * t34 + t470 * t85 - t471 * t86) * t526 + (qJD(2) * t38 + t137 * t470 - t138 * t471 + t335 * t59 + t338 * t58) * t527 + (qJD(2) * t51 + t143 * t470 - t144 * t471 + t335 * t77 + t338 * t76) * t528 + (qJD(2) * t80 + t111 * t338 + t112 * t335 + t195 * t470 - t196 * t471) * t529) * t334; 0.4e1 * (t529 + t420) * (-0.1e1 + t473) * t427; m(7) * (-t210 * t53 + t212 * t54 + t24 * t277 + t25 * t279) + m(5) * (-t210 * t95 + t212 * t96 + t277 * t49 + t279 * t50) + m(6) * (-t210 * t81 + t212 * t82 + t277 * t36 + t279 * t37); m(7) * (-t210 * t86 + t212 * t85 + t277 * t35 + t279 * t34 + t331 * t372) + m(6) * (t137 * t212 - t138 * t210 + t277 * t59 + t279 * t58 + t331 * t368) + m(5) * (t143 * t212 - t144 * t210 + t277 * t77 + t279 * t76 + (t20 * t334 + t467 * t51) * t331); 0.2e1 * t420 * ((-t210 * t338 + t212 * t335 + (t277 * t338 - t279 * t335) * qJD(1)) * t334 + ((t277 * t335 + t279 * t338) * t337 + t474 * t331) * qJD(2)); 0.4e1 * t420 * (t331 ^ 2 * t427 - t210 * t279 + t212 * t277); m(7) * (-t211 * t53 + t213 * t54 + t24 * t278 + t25 * t280) + m(6) * (-t211 * t81 + t213 * t82 + t278 * t36 + t280 * t37); m(7) * (-t211 * t86 + t213 * t85 + t278 * t35 + t280 * t34 + t332 * t372) + m(6) * (t137 * t213 - t138 * t211 + t278 * t59 + t280 * t58 + t332 * t368); ((-t211 * t338 + t213 * t335 + (t278 * t338 - t280 * t335) * qJD(1)) * t334 + ((t278 * t335 + t280 * t338) * t337 + t474 * t332) * qJD(2)) * t542; (0.2e1 * t331 * t332 * t427 - t210 * t280 - t211 * t279 + t212 * t278 + t213 * t277) * t542; 0.4e1 * t464 * (t332 ^ 2 * t427 - t211 * t280 + t213 * t278); m(7) * (t22 * t54 + t23 * t53 + t24 * t75 + t25 * t74) + t56 + (-t21 + (t335 * t459 + t338 * t458) * qJD(2)) * t337 + ((t11 / 0.2e1 + t16 / 0.2e1) * t338 + (t10 / 0.2e1 + t17 / 0.2e1) * t335 + (-t335 * t458 + t338 * t459) * qJD(1)) * t334; m(7) * (t14 * t26 + t22 * t85 + t23 * t86 + t34 * t74 + t35 * t75 + t5 * t55) + (-t2 / 0.2e1 - t401 * t423 + t13 * t514 + (qJD(1) * t33 - t10) * t524) * t338 + (t12 * t514 - t403 * t423 + t1 / 0.2e1 + (qJD(1) * t32 + t11) * t524) * t335 + (t4 * t525 + t3 * t522 - qJD(2) * t399 / 0.2e1 + (t401 * t525 - t403 * t522) * qJD(1)) * t334; m(7) * ((-t14 + (t335 * t75 + t338 * t74) * qJD(2)) * t337 + (qJD(2) * t55 + t22 * t335 + t23 * t338 + (-t335 * t74 + t338 * t75) * qJD(1)) * t334); m(7) * (-t210 * t74 + t212 * t75 + t22 * t277 + t23 * t279 + t331 * t369); m(7) * (-t211 * t74 + t213 * t75 + t22 * t278 + t23 * t280 + t332 * t369); (t14 * t55 + t22 * t75 + t23 * t74) * t530 + (t21 * t337 - t56 + (t335 * t12 + t338 * t13 - t337 * t398) * qJD(2)) * t337 + (t338 * t1 + t335 * t2 - t337 * (t10 * t335 + t11 * t338) + (t334 * t398 - t57 * t337) * qJD(2) + (t338 * t12 - t335 * t13 - t337 * t399) * qJD(1)) * t334;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
