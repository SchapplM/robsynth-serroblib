% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR6_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:35
% EndTime: 2019-03-09 02:30:43
% DurationCPUTime: 6.37s
% Computational Cost: add. (11820->529), mult. (23297->707), div. (0->0), fcn. (21982->6), ass. (0->290)
t337 = sin(qJ(6));
t338 = sin(qJ(5));
t340 = cos(qJ(5));
t527 = cos(qJ(6));
t407 = t527 * t340;
t285 = t337 * t338 - t407;
t537 = -t285 / 0.2e1;
t286 = t337 * t340 + t338 * t527;
t584 = -t286 / 0.2e1;
t341 = cos(qJ(4));
t528 = t341 / 0.2e1;
t334 = pkin(1) + qJ(3);
t339 = sin(qJ(4));
t526 = pkin(4) * t339;
t296 = -pkin(8) * t341 + t334 + t526;
t274 = t340 * t296;
t333 = qJ(2) - pkin(7);
t393 = -t333 * t338 + pkin(5);
t438 = t340 * t341;
t424 = pkin(9) * t438;
t165 = t339 * t393 + t274 - t424;
t442 = t339 * t340;
t213 = t296 * t338 + t333 * t442;
t444 = t338 * t341;
t179 = -pkin(9) * t444 + t213;
t449 = t337 * t179;
t80 = t165 * t527 - t449;
t445 = t338 * t339;
t212 = -t333 * t445 + t274;
t178 = t212 - t424;
t90 = t178 * t527 - t449;
t513 = t80 - t90;
t331 = t338 ^ 2;
t332 = t340 ^ 2;
t429 = t331 + t332;
t583 = mrSges(6,3) * t429;
t562 = t338 * mrSges(6,1) + t340 * mrSges(6,2);
t574 = t562 * t341;
t582 = t333 * t574;
t188 = -Ifges(7,5) * t285 - Ifges(7,6) * t286;
t546 = -pkin(9) - pkin(8);
t305 = t546 * t338;
t307 = t546 * t340;
t203 = t337 * t305 - t307 * t527;
t387 = t527 * t305 + t307 * t337;
t38 = -t203 * mrSges(7,1) - t387 * mrSges(7,2) + t188;
t581 = t38 * qJD(6);
t252 = t286 * t341;
t232 = Ifges(7,4) * t252;
t250 = t337 * t444 - t341 * t407;
t143 = -Ifges(7,1) * t250 + t339 * Ifges(7,5) - t232;
t153 = Ifges(7,2) * t250 - t232;
t580 = t153 + t143;
t253 = -t337 * t445 + t339 * t407;
t538 = -t253 / 0.2e1;
t577 = t339 / 0.4e1;
t529 = -t341 / 0.2e1;
t409 = t527 * t179;
t81 = t337 * t165 + t409;
t88 = -t337 * t178 - t409;
t512 = t81 + t88;
t251 = t286 * t339;
t391 = -t253 * mrSges(7,1) + t251 * mrSges(7,2);
t572 = qJD(6) * t391;
t563 = t286 * mrSges(7,1) - t285 * mrSges(7,2);
t571 = qJD(6) * t563;
t492 = t252 * mrSges(7,3);
t207 = -mrSges(7,2) * t339 - t492;
t495 = t250 * mrSges(7,3);
t209 = mrSges(7,1) * t339 + t495;
t570 = t207 * t537 + t209 * t584;
t543 = -t209 / 0.2e1;
t545 = t207 / 0.2e1;
t569 = t203 * t543 + t387 * t545;
t540 = t251 / 0.2e1;
t371 = Ifges(7,5) * t538 + Ifges(7,6) * t540;
t330 = pkin(4) * t341;
t306 = pkin(8) * t339 + t330;
t289 = t340 * t306;
t176 = pkin(9) * t442 + t341 * t393 + t289;
t221 = t338 * t306 + t333 * t438;
t193 = pkin(9) * t445 + t221;
t89 = t176 * t527 - t337 * t193;
t91 = t337 * t176 + t193 * t527;
t556 = Ifges(7,3) * t528 + t89 * mrSges(7,1) / 0.2e1 - t91 * mrSges(7,2) / 0.2e1 + t371;
t483 = t340 * mrSges(6,1);
t484 = t338 * mrSges(6,2);
t560 = t484 - t483;
t568 = t560 / 0.2e1;
t486 = t286 * mrSges(7,3);
t566 = t339 * t341;
t493 = t252 * mrSges(7,1);
t496 = t250 * mrSges(7,2);
t152 = t493 - t496;
t564 = t152 + t574;
t329 = Ifges(6,4) * t340;
t504 = Ifges(6,2) * t338;
t561 = t329 - t504;
t302 = Ifges(6,1) * t338 + t329;
t220 = -t333 * t444 + t289;
t373 = -t220 * t338 + t221 * t340;
t295 = t339 * mrSges(6,1) - mrSges(6,3) * t438;
t439 = t340 * t295;
t559 = -t439 / 0.2e1 + t529 * t583;
t487 = t286 * mrSges(7,2);
t490 = t285 * mrSges(7,1);
t431 = -t490 / 0.2e1 - t487 / 0.2e1;
t433 = -t493 / 0.2e1 + t496 / 0.2e1;
t413 = -t486 / 0.2e1;
t488 = t285 * mrSges(7,3);
t414 = -t488 / 0.2e1;
t558 = t250 * t413 - t252 * t414 - t570;
t150 = -mrSges(7,1) * t250 - mrSges(7,2) * t252;
t557 = t150 * t529 + t209 * t538 - t251 * t207 / 0.2e1;
t453 = t286 * t250;
t457 = t285 * t252;
t555 = (-t457 / 0.2e1 + t453 / 0.2e1) * mrSges(7,3) + t570;
t508 = Ifges(7,4) * t250;
t141 = -Ifges(7,2) * t252 + t339 * Ifges(7,6) - t508;
t154 = -Ifges(7,1) * t252 + t508;
t278 = Ifges(7,4) * t285;
t189 = -Ifges(7,2) * t286 - t278;
t507 = Ifges(7,4) * t286;
t190 = -Ifges(7,2) * t285 + t507;
t191 = -Ifges(7,1) * t285 - t507;
t192 = Ifges(7,1) * t286 - t278;
t524 = pkin(5) * t338;
t392 = -t333 + t524;
t283 = t392 * t341;
t523 = pkin(5) * t340;
t321 = -pkin(4) - t523;
t554 = t188 * t577 + t321 * t150 / 0.2e1 + t283 * t563 / 0.2e1 - (t192 + t189) * t252 / 0.4e1 - t580 * t285 / 0.4e1 + (t154 / 0.4e1 - t141 / 0.4e1) * t286 + (-t191 / 0.4e1 + t190 / 0.4e1) * t250;
t553 = t339 ^ 2;
t552 = t341 ^ 2;
t551 = -m(6) / 0.2e1;
t550 = m(6) / 0.2e1;
t549 = -m(7) / 0.2e1;
t548 = m(7) / 0.2e1;
t547 = m(7) * pkin(5);
t208 = mrSges(7,1) * t341 + t253 * mrSges(7,3);
t544 = t208 / 0.2e1;
t542 = -t250 / 0.2e1;
t541 = t250 / 0.2e1;
t539 = -t252 / 0.2e1;
t536 = t286 / 0.2e1;
t535 = -t338 / 0.2e1;
t534 = t338 / 0.2e1;
t533 = t339 / 0.2e1;
t532 = -t340 / 0.2e1;
t531 = t340 / 0.2e1;
t525 = pkin(5) * t337;
t522 = t80 * mrSges(7,2);
t521 = t81 * mrSges(7,1);
t520 = t88 * mrSges(7,1);
t518 = t90 * mrSges(7,2);
t511 = m(7) * qJD(2);
t509 = Ifges(6,4) * t338;
t506 = Ifges(6,5) * t339;
t505 = Ifges(6,5) * t340;
t503 = Ifges(6,6) * t338;
t502 = Ifges(6,6) * t339;
t494 = t251 * mrSges(7,1);
t491 = t253 * mrSges(7,2);
t140 = -Ifges(7,4) * t253 + Ifges(7,2) * t251 + Ifges(7,6) * t341;
t142 = -Ifges(7,1) * t253 + Ifges(7,4) * t251 + Ifges(7,5) * t341;
t151 = -t491 - t494;
t206 = -mrSges(7,2) * t341 + t251 * mrSges(7,3);
t247 = Ifges(6,6) * t341 - t339 * t561;
t379 = Ifges(6,1) * t340 - t509;
t248 = Ifges(6,5) * t341 - t339 * t379;
t272 = t562 * t339;
t282 = t392 * t339;
t292 = -mrSges(6,2) * t341 + mrSges(6,3) * t445;
t421 = mrSges(6,3) * t444;
t293 = -t339 * mrSges(6,2) - t421;
t294 = t341 * mrSges(6,1) + mrSges(6,3) * t442;
t327 = t339 * mrSges(5,2);
t377 = -t503 + t505;
t3 = -t334 * t327 + t213 * t292 + t221 * t293 + t212 * t294 + t220 * t295 - t282 * t152 + t283 * t151 + t140 * t539 + t143 * t538 + t142 * t542 + t141 * t540 + t81 * t206 + t91 * t207 + t80 * t208 + t89 * t209 + m(6) * (t212 * t220 + t213 * t221) + m(7) * (-t282 * t283 + t80 * t89 + t81 * t91) + (Ifges(7,5) * t542 + Ifges(7,6) * t539 + t334 * mrSges(5,1) + t248 * t531 + t247 * t535 + t333 * t272 + (t505 / 0.2e1 - t503 / 0.2e1 - Ifges(5,4)) * t341) * t341 + (t582 + (Ifges(5,4) - t377) * t339 + (Ifges(6,3) - Ifges(5,1) + Ifges(5,2) + Ifges(7,3) - Ifges(6,1) * t332 / 0.2e1 - m(6) * t333 ^ 2 + (t329 - t504 / 0.2e1) * t338) * t341 + t371) * t339;
t485 = t3 * qJD(1);
t432 = -Ifges(7,5) * t252 + Ifges(7,6) * t250;
t479 = t81 * t250;
t364 = mrSges(7,3) * t479 + t283 * t150 + t432 * t533;
t6 = -pkin(5) * t152 * t438 - t88 * t209 - m(7) * (t80 * t88 + t81 * t90) - t90 * t207 + t213 * t295 + (t154 / 0.2e1 - t141 / 0.2e1) * t250 - (-t143 / 0.2e1 - t153 / 0.2e1 + t80 * mrSges(7,3)) * t252 + ((t506 + (-mrSges(6,2) * t333 - t509) * t341) * t338 + (t502 - t283 * t547 + t213 * mrSges(6,3) + (t333 * mrSges(6,1) + t329 + (Ifges(6,1) - Ifges(6,2)) * t338) * t341) * t340) * t341 - t364 + (-t293 - t421) * t212;
t482 = t6 * qJD(1);
t7 = -t209 * t81 + t154 * t542 + t141 * t541 + (t207 + t492) * t80 + t364 + t580 * t539;
t481 = t7 * qJD(1);
t480 = t80 * t285;
t476 = t560 - mrSges(5,1);
t382 = m(5) * (t552 + t553);
t351 = t382 / 0.2e1 + (t429 * t553 + t552) * t550 + (t251 ^ 2 + t253 ^ 2 + t552) * t548;
t358 = -m(5) / 0.2e1 + t429 * t551 + (-t285 ^ 2 - t286 ^ 2) * t548;
t52 = -m(4) - t351 + t358;
t475 = qJD(1) * t52;
t402 = -t445 / 0.2e1;
t346 = (-t251 * t512 - t253 * t513 - t523 * t552) * t548 + t552 * t568 + t293 * t402 + t559 * t339;
t426 = t547 / 0.2e1;
t355 = -t484 / 0.2e1 + t483 / 0.2e1 + (-t285 * t527 + t286 * t337) * t426 + t431;
t416 = t492 / 0.2e1;
t418 = -t495 / 0.2e1;
t10 = t251 * t416 + t253 * t418 - t346 + t355 - t557;
t474 = t10 * qJD(1);
t415 = -t492 / 0.2e1;
t417 = t495 / 0.2e1;
t376 = t251 * t415 + t253 * t417 + t557;
t14 = t376 - t431;
t473 = t14 * qJD(1);
t370 = -t494 / 0.2e1 - t491 / 0.2e1;
t21 = t370 + t555;
t472 = t21 * qJD(1);
t374 = -t212 * t338 + t213 * t340;
t440 = t340 * t293;
t446 = t338 * t295;
t23 = -mrSges(4,2) - mrSges(3,3) + t251 * t209 - t253 * t207 + (mrSges(5,3) * t341 + t564) * t341 - m(7) * (-t251 * t80 + t253 * t81 - t283 * t341) + (-m(6) * t552 - t382) * t333 + (-m(6) * t374 + mrSges(5,3) * t339 - t440 + t446) * t339 + (-m(4) - m(3)) * qJ(2);
t469 = t23 * qJD(1);
t466 = t250 * t337;
t448 = t338 * t293;
t32 = t339 * mrSges(5,1) + t341 * mrSges(5,2) + t286 * t207 - t285 * t209 + t448 + t439 + mrSges(4,3) + m(7) * (t81 * t286 - t480) + m(6) * (t212 * t340 + t213 * t338) + 0.4e1 * (m(5) / 0.4e1 + m(4) / 0.4e1) * t334;
t452 = t32 * qJD(1);
t450 = t333 * t341;
t447 = t338 * t294;
t441 = t340 * t292;
t436 = t341 * t563;
t428 = qJD(4) * t339;
t427 = mrSges(7,3) * t525;
t97 = t453 - t457;
t423 = t97 * t548;
t425 = qJD(4) * t423;
t422 = pkin(5) * t527;
t412 = t486 / 0.2e1;
t410 = -t436 / 0.2e1 + (t413 + t412) * t253;
t408 = t527 * t252;
t401 = -t444 / 0.2e1;
t399 = -t442 / 0.2e1;
t397 = -t438 / 0.2e1;
t396 = t438 / 0.2e1;
t187 = t487 + t490;
t395 = -t560 / 0.2e1 - t187 / 0.2e1;
t394 = pkin(8) * t429;
t386 = mrSges(7,3) * t422;
t385 = t81 * t412;
t301 = Ifges(6,2) * t340 + t509;
t350 = (-t251 * t527 + t253 * t337) * t426 + mrSges(6,1) * t402 + mrSges(6,2) * t399 + t370;
t369 = t446 / 0.2e1 - t440 / 0.2e1;
t356 = (t285 * t512 + t286 * t513) * t548 + t369;
t12 = t350 - t356 + t555;
t375 = t12 * qJD(1);
t372 = t251 * t285 + t253 * t286;
t368 = t301 * t534 + t302 * t532;
t347 = (t251 * t536 + t253 * t537) * mrSges(7,3) + (t332 / 0.2e1 + t331 / 0.2e1) * t339 * mrSges(6,3) + (t339 * t394 + t330) * t550 + (t203 * t253 - t251 * t387 - t321 * t341) * t548;
t348 = (-t340 * t220 - t338 * t221) * t550 + (t285 * t89 - t286 * t91) * t548 + t285 * t544 + t206 * t584 + t292 * t535 + t294 * t532;
t17 = -t327 + (mrSges(5,1) + t395) * t341 + t347 - t348;
t366 = -t17 * qJD(1) + qJD(3) * t423;
t363 = (t203 * t286 - t285 * t387) * t548;
t28 = t321 * t563 + (t191 / 0.2e1 - t190 / 0.2e1) * t286 + (-t192 / 0.2e1 - t189 / 0.2e1) * t285;
t35 = t410 - t433;
t357 = t203 * t417 + t387 * t416 + t81 * t413 + t554;
t349 = t357 + t385 + t569;
t5 = t349 - t556;
t362 = t5 * qJD(1) + t35 * qJD(3) + t28 * qJD(4);
t354 = (t337 * t543 + t527 * t545 + (t408 / 0.2e1 + t466 / 0.2e1) * mrSges(7,3)) * pkin(5);
t20 = (-t80 / 0.2e1 + t90 / 0.2e1) * mrSges(7,2) + (-t81 / 0.2e1 - t88 / 0.2e1) * mrSges(7,1) + t354;
t297 = (mrSges(7,1) * t337 + mrSges(7,2) * t527) * pkin(5);
t361 = -t20 * qJD(1) + t297 * qJD(5);
t62 = m(6) * (-0.1e1 + t429) * t566 + m(7) * (-t253 * t250 + t251 * t252 - t566);
t344 = (t374 * t341 + (t373 - 0.2e1 * t450) * t339) * t551 + (-t251 * t89 - t80 * t252 + t253 * t91 + t341 * t282 + t283 * t339 - t479) * t549 + t207 * t541 + t208 * t540 - t252 * t543 + t206 * t538;
t8 = (t151 / 0.2e1 - t272 / 0.2e1 + t369) * t341 + (-t152 / 0.2e1 - t574 / 0.2e1 + t447 / 0.2e1 - t441 / 0.2e1) * t339 + t363 + t344;
t360 = t8 * qJD(1) - t62 * qJD(3) - t97 * t511 / 0.2e1;
t19 = -pkin(4) * t562 + t379 * t534 + t531 * t561 + t28 - t368 + (m(7) * t321 + t187) * t524;
t343 = ((t283 * t338 + t321 * t438) * pkin(5) + t512 * t387 - t513 * t203) * t548 - t338 * (t341 * t561 + t502) / 0.4e1 + t377 * t577 + t330 * t568 + t152 * t524 / 0.2e1 + t301 * t397 - t582 / 0.2e1 + mrSges(7,3) * t480 / 0.2e1 + t88 * t413 + t90 * t414 + pkin(5) * t187 * t396 + (0.2e1 * t379 * t341 + t506) * t340 / 0.4e1 + (-t448 / 0.2e1 + t559) * pkin(8) + t569 - (0.2e1 * t302 + t561) * t444 / 0.4e1;
t345 = Ifges(6,3) * t528 + t220 * mrSges(6,1) / 0.2e1 - t221 * mrSges(6,2) / 0.2e1 + (t337 * t91 + t527 * t89) * t426 + Ifges(6,5) * t399 + Ifges(6,6) * t445 / 0.2e1 + t206 * t525 / 0.2e1 + t422 * t544 + t556;
t2 = t203 * t418 + t387 * t415 - t343 + t345 + t385 - t554;
t352 = (-t408 - t466) * t426 + mrSges(6,1) * t401 + mrSges(6,2) * t397 + t433;
t353 = -pkin(5) * t444 * t549 + t574 / 0.2e1;
t25 = t436 / 0.2e1 + t352 + t353;
t359 = -t2 * qJD(1) - t25 * qJD(3) + t19 * qJD(4);
t284 = t297 * qJD(6);
t53 = t351 + t358;
t36 = t410 + t433;
t26 = t352 - t353 + t410;
t22 = t370 + t558;
t18 = t341 * t395 + t347 + t348;
t16 = -t522 / 0.2e1 - t521 / 0.2e1 - t518 / 0.2e1 + t520 / 0.2e1 + t354 + t432;
t15 = t376 + t431;
t13 = t350 + t356 + t558;
t11 = t346 + t355 + t376;
t9 = t293 * t396 + t294 * t402 + t295 * t401 - t344 + t363 + (t151 - t272) * t529 + (t441 + t564) * t533;
t4 = t349 + t556;
t1 = t345 + t343 + t357;
t24 = [-qJD(2) * t23 + qJD(3) * t32 + qJD(4) * t3 - qJD(5) * t6 + qJD(6) * t7, t53 * qJD(3) + t18 * qJD(4) + t13 * qJD(5) + t22 * qJD(6) - t372 * t511 - t469, m(7) * qJD(3) * t372 + t53 * qJD(2) + t9 * qJD(4) + t11 * qJD(5) + t15 * qJD(6) + t452, t485 + t18 * qJD(2) + t9 * qJD(3) + t1 * qJD(5) + t4 * qJD(6) + (-Ifges(5,5) + (-m(6) * pkin(4) + t476) * t333 + t368) * t428 + (-Ifges(5,6) * t341 + t247 * t531 + t248 * t534 + t321 * t151 + t142 * t536 + t140 * t537 - t282 * t187 + pkin(4) * t272 - mrSges(5,2) * t450 + m(7) * (t203 * t91 - t282 * t321 + t387 * t89) - t91 * t488 - t89 * t486 + t192 * t538 + t190 * t540 + t203 * t206 + t387 * t208 + (m(6) * t373 + t441 - t447) * pkin(8) + (Ifges(6,5) * t338 + Ifges(7,5) * t286 + Ifges(6,6) * t340 - Ifges(7,6) * t285) * t528 + t373 * mrSges(6,3)) * qJD(4), -t482 + t13 * qJD(2) + t11 * qJD(3) + t1 * qJD(4) + (-Ifges(6,5) * t444 - Ifges(6,6) * t438 + t252 * t386 + (t337 * t90 + t527 * t88) * t547 + t250 * t427 - t518 + t520 - t212 * mrSges(6,2) - t213 * mrSges(6,1) + t432) * qJD(5) + t16 * qJD(6), t481 + t22 * qJD(2) + t15 * qJD(3) + t4 * qJD(4) + t16 * qJD(5) + (t432 - t521 - t522) * qJD(6); qJD(3) * t52 - qJD(4) * t17 - qJD(5) * t12 - qJD(6) * t21 + t469, 0, t425 + t475, t366 (m(7) * (t285 * t525 + t286 * t422) + t563 + t562) * qJD(5) + t571 - t375, qJD(5) * t563 - t472 + t571; -qJD(2) * t52 - qJD(4) * t8 - qJD(5) * t10 + qJD(6) * t14 - t452, t425 - t475, t62 * qJD(4), t26 * qJD(5) + t36 * qJD(6) + (t187 + t476) * t428 - t360 + ((t250 * t285 + t252 * t286) * mrSges(7,3) + (-mrSges(5,2) + t583) * t341 + 0.2e1 * (-t203 * t250 - t252 * t387 + t321 * t339) * t548 + 0.2e1 * (t341 * t394 - t526) * t550) * qJD(4), -t474 + t26 * qJD(4) + (m(7) * (-t251 * t525 - t253 * t422) + t560 * t339 + t391) * qJD(5) + t572, t36 * qJD(4) + qJD(5) * t391 + t473 + t572; qJD(2) * t17 + qJD(3) * t8 - qJD(5) * t2 + qJD(6) * t5 - t485, -t366, -t25 * qJD(5) + t35 * qJD(6) + t360, qJD(5) * t19 + qJD(6) * t28 (-t286 * t427 + t285 * t386 + (-t203 * t527 + t337 * t387) * t547 + t377 + t560 * pkin(8) + t38) * qJD(5) + t581 + t359, t38 * qJD(5) + t362 + t581; qJD(2) * t12 + qJD(3) * t10 + qJD(4) * t2 + qJD(6) * t20 + t482, t375, t25 * qJD(4) + t474, -t359, -t284, -t284 - t361; qJD(2) * t21 - qJD(3) * t14 - qJD(4) * t5 - qJD(5) * t20 - t481, t472, -t35 * qJD(4) - t473, -t362, t361, 0;];
Cq  = t24;
