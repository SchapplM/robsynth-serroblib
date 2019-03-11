% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:24
% EndTime: 2019-03-09 02:20:37
% DurationCPUTime: 7.11s
% Computational Cost: add. (22939->468), mult. (44487->638), div. (0->0), fcn. (49778->10), ass. (0->272)
t318 = sin(pkin(11));
t319 = cos(pkin(11));
t323 = sin(qJ(4));
t500 = cos(qJ(4));
t296 = t318 * t323 - t319 * t500;
t321 = sin(qJ(6));
t322 = sin(qJ(5));
t324 = cos(qJ(5));
t499 = cos(qJ(6));
t357 = t321 * t324 + t322 * t499;
t226 = t357 * t296;
t356 = t321 * t322 - t499 * t324;
t228 = t356 * t296;
t99 = -t226 * t356 + t228 * t357;
t560 = m(7) * t99 * qJD(3);
t298 = t318 * t500 + t323 * t319;
t227 = t357 * t298;
t469 = t227 * mrSges(7,3);
t168 = -mrSges(7,2) * t296 - t469;
t520 = t168 / 0.2e1;
t507 = t298 / 0.2e1;
t505 = t357 / 0.2e1;
t526 = m(7) / 0.2e1;
t559 = t526 * t99;
t367 = -cos(pkin(10)) * pkin(1) - pkin(3) * t319 - pkin(2);
t494 = pkin(8) * t298;
t224 = pkin(4) * t296 + t367 - t494;
t381 = sin(pkin(10)) * pkin(1) + qJ(3);
t366 = pkin(7) + t381;
t289 = t366 * t319;
t350 = t318 * t366;
t235 = t289 * t500 - t323 * t350;
t129 = t324 * t224 - t322 * t235;
t433 = t298 * t324;
t101 = -pkin(9) * t433 + t129;
t130 = t322 * t224 + t235 * t324;
t434 = t298 * t322;
t102 = -pkin(9) * t434 + t130;
t425 = t321 * t102;
t358 = t101 * t499 - t425;
t94 = t296 * pkin(5) + t101;
t55 = t499 * t94 - t425;
t349 = t55 - t358;
t250 = -Ifges(7,5) * t356 - Ifges(7,6) * t357;
t522 = -pkin(9) - pkin(8);
t307 = t522 * t322;
t308 = t522 * t324;
t267 = t321 * t307 - t308 * t499;
t387 = t499 * t307 + t308 * t321;
t50 = -t267 * mrSges(7,1) - t387 * mrSges(7,2) + t250;
t558 = t50 * qJD(6);
t221 = Ifges(7,4) * t227;
t225 = t356 * t298;
t117 = -Ifges(7,1) * t225 + Ifges(7,5) * t296 - t221;
t138 = Ifges(7,2) * t225 - t221;
t557 = t117 + t138;
t555 = t267 / 0.2e1;
t554 = t296 / 0.4e1;
t502 = t324 / 0.2e1;
t510 = t387 / 0.2e1;
t398 = t499 * t102;
t360 = -t321 * t101 - t398;
t56 = t321 * t94 + t398;
t348 = t56 + t360;
t135 = -t225 * mrSges(7,1) - t227 * mrSges(7,2);
t549 = qJD(6) * t135;
t248 = mrSges(7,1) * t357 - mrSges(7,2) * t356;
t548 = t248 * qJD(6);
t471 = t225 * mrSges(7,3);
t170 = mrSges(7,1) * t296 + t471;
t536 = t170 * t505 + t356 * t520;
t515 = t228 / 0.2e1;
t517 = t226 / 0.2e1;
t362 = Ifges(7,5) * t515 + Ifges(7,6) * t517;
t234 = t289 * t323 + t500 * t350;
t497 = pkin(4) * t298;
t247 = pkin(8) * t296 + t497;
t144 = -t324 * t234 + t322 * t247;
t436 = t296 * t322;
t120 = pkin(9) * t436 + t144;
t143 = t322 * t234 + t324 * t247;
t435 = t296 * t324;
t98 = t298 * pkin(5) + pkin(9) * t435 + t143;
t71 = -t321 * t120 + t499 * t98;
t72 = t120 * t499 + t321 * t98;
t531 = Ifges(7,3) * t507 + t71 * mrSges(7,1) / 0.2e1 - t72 * mrSges(7,2) / 0.2e1 + t362;
t547 = -t170 / 0.2e1;
t317 = t324 ^ 2;
t546 = -t317 / 0.2e1;
t462 = t357 * mrSges(7,3);
t534 = -mrSges(6,1) * t324 + t322 * mrSges(6,2);
t543 = t534 - mrSges(5,1);
t355 = t534 * t298;
t313 = Ifges(6,4) * t324;
t477 = Ifges(6,2) * t322;
t541 = t313 - t477;
t316 = t322 ^ 2;
t539 = t317 + t316;
t394 = t436 / 0.2e1;
t458 = t324 * mrSges(6,2);
t399 = t458 / 0.2e1;
t524 = m(7) * pkin(5);
t414 = t524 / 0.2e1;
t468 = t228 * mrSges(7,2);
t470 = t226 * mrSges(7,1);
t420 = t470 / 0.2e1 - t468 / 0.2e1;
t334 = (t226 * t499 + t228 * t321) * t414 + mrSges(6,1) * t394 + t296 * t399 + t420;
t518 = t225 / 0.2e1;
t396 = t357 * t518;
t382 = mrSges(7,3) * t396;
t331 = -t382 + t334;
t401 = t462 / 0.2e1;
t383 = t225 * t401;
t538 = t331 + t383;
t409 = mrSges(6,3) * t434;
t363 = -t296 * mrSges(6,2) - t409;
t364 = t296 * mrSges(6,1) - mrSges(6,3) * t433;
t537 = t322 * t363 + t324 * t364;
t368 = -t143 * t322 + t144 * t324;
t535 = -t322 * Ifges(6,1) - t313;
t533 = mrSges(6,3) * t539;
t463 = t356 * mrSges(7,3);
t402 = t463 / 0.2e1;
t532 = -t227 * t402 + t383 - t536;
t481 = Ifges(7,4) * t225;
t115 = -Ifges(7,2) * t227 + t296 * Ifges(7,6) - t481;
t139 = -Ifges(7,1) * t227 + t481;
t173 = pkin(5) * t434 + t234;
t295 = Ifges(7,4) * t356;
t251 = -Ifges(7,2) * t357 - t295;
t480 = Ifges(7,4) * t357;
t252 = -Ifges(7,2) * t356 + t480;
t253 = -Ifges(7,1) * t356 - t480;
t254 = Ifges(7,1) * t357 - t295;
t312 = -pkin(5) * t324 - pkin(4);
t530 = t471 * t555 + t469 * t510 + t312 * t135 / 0.2e1 + t250 * t554 + t173 * t248 / 0.2e1 - (t254 + t251) * t227 / 0.4e1 - t557 * t356 / 0.4e1 + (t139 / 0.4e1 - t115 / 0.4e1) * t357 + (-t253 / 0.4e1 + t252 / 0.4e1) * t225;
t528 = 2 * qJD(4);
t527 = m(6) / 0.2e1;
t137 = mrSges(7,1) * t227 - mrSges(7,2) * t225;
t521 = t137 / 0.2e1;
t519 = -t225 / 0.2e1;
t516 = -t227 / 0.2e1;
t249 = mrSges(7,1) * t356 + mrSges(7,2) * t357;
t513 = t249 / 0.2e1;
t512 = t252 / 0.2e1;
t511 = t254 / 0.2e1;
t509 = -t267 / 0.2e1;
t508 = t296 / 0.2e1;
t506 = -t356 / 0.2e1;
t504 = -t322 / 0.2e1;
t503 = t322 / 0.2e1;
t498 = m(7) * t173;
t496 = pkin(5) * t321;
t495 = pkin(5) * t322;
t493 = t55 * mrSges(7,2);
t492 = t56 * mrSges(7,1);
t483 = mrSges(6,3) * t296;
t482 = Ifges(6,4) * t322;
t479 = Ifges(6,5) * t296;
t478 = Ifges(6,5) * t324;
t476 = Ifges(6,6) * t296;
t475 = Ifges(6,6) * t322;
t473 = t143 * mrSges(6,1);
t472 = t144 * mrSges(6,2);
t461 = t322 * mrSges(6,1);
t263 = t296 * t298;
t346 = m(7) * (-t225 * t228 - t226 * t227 + t263);
t388 = t539 * t296;
t347 = m(6) * (-t298 * t388 + t263);
t41 = t346 / 0.2e1 + t347 / 0.2e1;
t455 = qJD(1) * t41;
t16 = t135 * t508 + t168 * t516 + t170 * t518 - (t225 ^ 2 + t227 ^ 2) * mrSges(7,3) / 0.2e1;
t343 = t355 / 0.2e1;
t412 = pkin(5) * t433;
t10 = -m(7) * (t225 * t349 - t227 * t348 + t296 * t412) / 0.2e1 + t296 * t343 - t16 + t537 * t507 + t298 ^ 2 * t533 / 0.2e1;
t454 = t10 * qJD(1);
t327 = (-t348 * t356 - t349 * t357) * t526 + t364 * t504 + t363 * t502;
t395 = t227 * t506;
t11 = -mrSges(7,3) * t395 - t327 + t331 + t536;
t453 = t11 * qJD(1);
t450 = t16 * qJD(1);
t17 = (-t395 - t396) * mrSges(7,3) + t420 + t536;
t449 = t17 * qJD(1);
t445 = t226 * t170;
t444 = t226 * t357;
t441 = t228 * t168;
t440 = t228 * t356;
t439 = t234 * t298;
t438 = t296 * t248;
t246 = t298 * mrSges(6,1) + mrSges(6,3) * t435;
t423 = t322 * t246;
t245 = -mrSges(6,2) * t298 + mrSges(6,3) * t436;
t422 = t324 * t245;
t419 = -Ifges(7,5) * t227 + Ifges(7,6) * t225;
t415 = mrSges(7,3) * t496;
t413 = qJD(4) * t559;
t411 = pkin(5) * t436;
t410 = pkin(5) * t499;
t400 = -t462 / 0.2e1;
t390 = t534 / 0.2e1 + t513;
t386 = mrSges(7,3) * t410;
t384 = t56 * t401;
t379 = t458 + t461;
t378 = Ifges(6,1) * t324 - t482;
t305 = t324 * Ifges(6,2) + t482;
t376 = -t475 + t478;
t114 = Ifges(7,4) * t228 + Ifges(7,2) * t226 + t298 * Ifges(7,6);
t116 = Ifges(7,1) * t228 + Ifges(7,4) * t226 + t298 * Ifges(7,5);
t136 = t468 - t470;
t167 = -mrSges(7,2) * t298 + t226 * mrSges(7,3);
t169 = mrSges(7,1) * t298 - t228 * mrSges(7,3);
t174 = t235 - t411;
t194 = Ifges(6,6) * t298 - t296 * t541;
t195 = Ifges(6,5) * t298 - t296 * t378;
t240 = t296 * t379;
t241 = t379 * t298;
t290 = t296 * mrSges(5,2);
t369 = t143 * t324 + t322 * t144;
t4 = -t367 * t290 + t130 * t245 + t129 * t246 - t234 * t240 + t235 * t241 + t114 * t516 + t117 * t515 + t116 * t519 + t115 * t517 + t56 * t167 + t72 * t168 + t55 * t169 + t71 * t170 + t173 * t136 + t174 * t137 + m(6) * (t129 * t143 + t130 * t144 + t234 * t235) + m(7) * (t173 * t174 + t55 * t71 + t56 * t72) + (t195 * t502 + t194 * t504 + t367 * mrSges(5,1) + Ifges(7,5) * t519 + Ifges(7,6) * t516 + (t478 / 0.2e1 - t475 / 0.2e1 - Ifges(5,4)) * t298 - t369 * mrSges(6,3)) * t298 + (t473 - t472 + (Ifges(6,1) * t546 + Ifges(6,3) + Ifges(5,2) + Ifges(7,3) - Ifges(5,1) + (t313 - t477 / 0.2e1) * t322) * t298 + t362 + (Ifges(5,4) - t376) * t296) * t296;
t370 = t129 * t322 - t130 * t324;
t374 = t226 * t55 + t228 * t56;
t9 = (-t225 * t72 - t227 * t71 + t374) * t526 + t441 / 0.2e1 + t167 * t519 + t445 / 0.2e1 + t169 * t516 + (t422 / 0.2e1 - t423 / 0.2e1 + t241 / 0.2e1 + t521 + t498 / 0.2e1 + (t234 + t368) * t527) * t298 + (-t240 / 0.2e1 + t136 / 0.2e1 + (t399 + t461 / 0.2e1) * t296 + t174 * t526 + (t235 + t370) * t527) * t296;
t375 = t4 * qJD(1) + t9 * qJD(2);
t352 = t173 * t135 + t419 * t508 + t56 * t471;
t5 = -t358 * t168 - t360 * t170 - m(7) * (t358 * t56 + t360 * t55) - t137 * t412 + t234 * t355 + t130 * t364 + (-t115 / 0.2e1 + t139 / 0.2e1) * t225 - (-t138 / 0.2e1 - t117 / 0.2e1 + t55 * mrSges(7,3)) * t227 + ((-Ifges(6,4) * t434 + t479) * t322 + (t476 - pkin(5) * t498 + t130 * mrSges(6,3) + (t313 + (Ifges(6,1) - Ifges(6,2)) * t322) * t298) * t324) * t298 - t352 + (-t363 - t409) * t129;
t373 = -t5 * qJD(1) - t10 * qJD(2);
t8 = -t170 * t56 + t115 * t518 + t139 * t519 + (t168 + t469) * t55 + t352 + t557 * t516;
t372 = t8 * qJD(1) + t16 * qJD(2);
t15 = t441 + t445 + (mrSges(5,3) * t298 + t137 + t241) * t298 + (mrSges(5,3) + t379) * t296 ^ 2 + m(7) * (t173 * t298 + t374) + m(6) * (t296 * t370 + t439) + m(5) * (-t235 * t296 + t439) + (m(4) * t381 + mrSges(4,3)) * (t318 ^ 2 + t319 ^ 2);
t371 = -t15 * qJD(1) - t41 * qJD(2);
t365 = (-pkin(8) * t388 - t497) * t527 + t526 * (t226 * t387 + t228 * t267 + t298 * t312);
t201 = t438 / 0.2e1;
t361 = t383 + t201 - t382;
t359 = t305 * t503 + t502 * t535;
t328 = t369 * t527 + (-t356 * t71 + t357 * t72) * t526 + t169 * t506 + t167 * t505 + t245 * t503 + t246 * t502;
t333 = (-t440 / 0.2e1 - t444 / 0.2e1) * mrSges(7,3) + (t546 - t316 / 0.2e1) * t483 + t365;
t19 = t290 + (-mrSges(5,1) + t390) * t298 - t328 + t333;
t353 = -t19 * qJD(1) + qJD(2) * t559;
t342 = t360 * mrSges(7,1);
t341 = t358 * mrSges(7,2);
t329 = (t321 * t547 + t499 * t520 + (t499 * t227 / 0.2e1 + t321 * t518) * mrSges(7,3)) * pkin(5) - t493 / 0.2e1 - t492 / 0.2e1;
t335 = -t342 / 0.2e1 + t341 / 0.2e1;
t14 = t329 + t335;
t303 = (mrSges(7,1) * t321 + mrSges(7,2) * t499) * pkin(5);
t65 = (-t387 / 0.2e1 + t510) * mrSges(7,2) + (t509 + t555) * mrSges(7,1);
t340 = -t14 * qJD(1) - t65 * qJD(4) + t303 * qJD(5);
t40 = t346 + t347;
t339 = -t9 * qJD(1) - t40 * qJD(2) - t560 / 0.2e1;
t332 = t379 * t508 + t411 * t526;
t29 = -t438 / 0.2e1 - t332 + t538;
t325 = ((t173 * t322 + t312 * t433) * pkin(5) + t348 * t387 - t349 * t267) * t526 + t234 * t379 / 0.2e1 + t267 * t547 + t168 * t510 + t376 * t554 - t322 * (t298 * t541 + t476) / 0.4e1 + pkin(4) * t343 + t495 * t521 - t305 * t433 / 0.2e1 + t55 * t402 + t360 * t400 - t358 * t463 / 0.2e1 + t412 * t513 + (0.2e1 * t378 * t298 + t479) * t324 / 0.4e1 - t537 * pkin(8) / 0.2e1 - t494 * t533 / 0.2e1 + (t535 / 0.2e1 - t541 / 0.4e1) * t434;
t326 = Ifges(6,3) * t507 + t473 / 0.2e1 - t472 / 0.2e1 + (t321 * t72 + t499 * t71) * t414 - Ifges(6,5) * t435 / 0.2e1 + Ifges(6,6) * t394 + t167 * t496 / 0.2e1 + t169 * t410 / 0.2e1 + t531;
t3 = -t530 + t384 + t326 - t325;
t38 = t312 * t248 - (t512 - t253 / 0.2e1) * t357 - (t251 / 0.2e1 + t511) * t356;
t31 = -pkin(4) * t379 + t378 * t503 + t502 * t541 - t359 + t38 + (m(7) * t312 + t249) * t495;
t338 = t3 * qJD(1) + t29 * qJD(2) - t31 * qJD(4);
t32 = t361 - t420;
t336 = t56 * t400 + t530;
t330 = t170 * t509 + t387 * t520 + t336 + t384;
t6 = t330 - t531;
t337 = -t6 * qJD(1) - t32 * qJD(2) - t38 * qJD(4);
t299 = t303 * qJD(6);
t33 = t361 + t420;
t30 = t201 + t332 + t538;
t20 = t298 * t390 + t328 + t333;
t18 = t420 + t532;
t13 = t329 - t335 + t419;
t12 = t327 + t334 + t532;
t7 = t330 + t531;
t2 = t336 + t326 + t325;
t1 = qJD(3) * t41 + qJD(4) * t9 - qJD(5) * t10 + qJD(6) * t16;
t21 = [qJD(3) * t15 + qJD(4) * t4 - qJD(5) * t5 + qJD(6) * t8, t1, t20 * qJD(4) + t12 * qJD(5) + t18 * qJD(6) - t371 + t560, t20 * qJD(3) + t2 * qJD(5) + t7 * qJD(6) + ((-pkin(4) * t235 + pkin(8) * t368) * t527 + (t174 * t312 + t267 * t72 + t387 * t71) * t526) * t528 + t375 + (-t71 * t462 - t72 * t463 + t194 * t502 + t195 * t503 + t312 * t136 + t116 * t505 - Ifges(5,6) * t298 + t114 * t506 + t267 * t167 + t387 * t169 + t174 * t249 + t226 * t512 + t228 * t511 + t234 * mrSges(5,2) + pkin(4) * t240 + (-Ifges(5,5) + t359) * t296 + (Ifges(6,5) * t322 + Ifges(7,5) * t357 + Ifges(6,6) * t324 - Ifges(7,6) * t356) * t507 + t543 * t235 + (t422 - t423) * pkin(8) + t368 * mrSges(6,3)) * qJD(4), t12 * qJD(3) + t2 * qJD(4) + (-Ifges(6,5) * t434 - Ifges(6,6) * t433 - t341 + t342 + t227 * t386 + (-t321 ^ 2 - t499 ^ 2) * t102 * t524 + t225 * t415 - t129 * mrSges(6,2) - t130 * mrSges(6,1) + t419) * qJD(5) + t13 * qJD(6) + t373, t18 * qJD(3) + t7 * qJD(4) + t13 * qJD(5) + (t419 - t492 - t493) * qJD(6) + t372; t1, t40 * qJD(4), t413 + t455, t30 * qJD(5) + t33 * qJD(6) + t365 * t528 - t339 + (t290 + (t249 + t543) * t298 - t539 * t483 + (-t440 - t444) * mrSges(7,3)) * qJD(4), -t454 + t30 * qJD(4) + (m(7) * (t225 * t410 - t227 * t496) + t355 - t135) * qJD(5) - t549, qJD(4) * t33 - qJD(5) * t135 + t450 - t549; -qJD(4) * t19 - qJD(5) * t11 - qJD(6) * t17 + t371, t413 - t455, 0, t353, -t453 - t548 + (-t248 - t379 + (-t356 * t496 - t357 * t410) * m(7)) * qJD(5), -qJD(5) * t248 - t449 - t548; qJD(3) * t19 - qJD(5) * t3 + qJD(6) * t6 - t375, -t29 * qJD(5) + t32 * qJD(6) + t339, -t353, qJD(5) * t31 + qJD(6) * t38 (-t357 * t415 + t356 * t386 + (-t267 * t499 + t321 * t387) * t524 + t376 + t534 * pkin(8) + t50) * qJD(5) + t558 - t338, t50 * qJD(5) - t337 + t558; qJD(3) * t11 + qJD(4) * t3 + qJD(6) * t14 - t373, t29 * qJD(4) + t454, t453, t65 * qJD(6) + t338, -t299, -t299 - t340; qJD(3) * t17 - qJD(4) * t6 - qJD(5) * t14 - t372, -t32 * qJD(4) - t450, t449, -t65 * qJD(5) + t337, t340, 0;];
Cq  = t21;
