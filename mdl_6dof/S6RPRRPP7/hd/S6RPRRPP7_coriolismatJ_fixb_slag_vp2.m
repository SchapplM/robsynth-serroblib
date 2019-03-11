% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:15
% EndTime: 2019-03-09 04:50:28
% DurationCPUTime: 6.91s
% Computational Cost: add. (6530->639), mult. (13400->816), div. (0->0), fcn. (10217->4), ass. (0->308)
t531 = Ifges(6,4) + Ifges(5,5);
t330 = sin(qJ(4));
t328 = t330 ^ 2;
t332 = cos(qJ(4));
t329 = t332 ^ 2;
t427 = t328 + t329;
t439 = t332 * qJ(5);
t499 = pkin(4) + pkin(5);
t536 = t499 * t330;
t538 = t536 - t439;
t537 = mrSges(6,2) + mrSges(5,3);
t530 = Ifges(5,6) + Ifges(7,6);
t252 = t330 * mrSges(6,1) - t332 * mrSges(6,3);
t333 = cos(qJ(3));
t202 = t333 * t252;
t320 = t332 * mrSges(7,2);
t459 = t330 * mrSges(7,1);
t253 = t320 - t459;
t203 = t253 * t333;
t535 = -t203 + t202;
t507 = m(6) / 0.4e1;
t424 = m(7) / 0.4e1 + t507;
t534 = 0.2e1 * t424;
t500 = m(6) + m(7);
t331 = sin(qJ(3));
t533 = -t331 / 0.2e1;
t488 = t331 / 0.2e1;
t485 = -t333 / 0.2e1;
t484 = t333 / 0.2e1;
t478 = m(7) * t331;
t417 = -t478 / 0.2e1;
t529 = pkin(8) - qJ(6);
t246 = t529 * t330;
t251 = t529 * t332;
t532 = m(7) * (t246 * t330 + t251 * t332);
t528 = mrSges(7,3) * t427;
t248 = t332 * mrSges(7,1) + t330 * mrSges(7,2);
t197 = t333 * t248;
t444 = t330 * t331;
t227 = -t333 * mrSges(5,2) + mrSges(5,3) * t444;
t237 = mrSges(6,2) * t444 + t333 * mrSges(6,3);
t527 = t227 + t237;
t443 = t330 * t333;
t293 = mrSges(7,3) * t443;
t319 = t331 * mrSges(7,2);
t228 = t293 + t319;
t456 = t331 * mrSges(5,2);
t229 = -mrSges(5,3) * t443 - t456;
t318 = t331 * mrSges(6,3);
t420 = mrSges(6,2) * t443;
t236 = t318 - t420;
t433 = t229 + t236;
t514 = t228 + t433;
t438 = t332 * t333;
t457 = t331 * mrSges(7,1);
t233 = -mrSges(7,3) * t438 - t457;
t235 = -mrSges(6,1) * t331 + mrSges(6,2) * t438;
t526 = t235 + t233;
t380 = t332 * mrSges(6,1) + t330 * mrSges(6,3);
t525 = -t380 - t248;
t476 = pkin(4) * t330;
t250 = -t439 + t476;
t482 = m(6) * t250;
t524 = t252 + t482;
t315 = t330 * qJ(5);
t429 = t332 * pkin(4) + t315;
t322 = Ifges(6,5) * t330;
t523 = Ifges(6,1) * t332 + t322;
t522 = Ifges(6,6) * t330 + t531 * t332;
t324 = Ifges(7,4) * t330;
t521 = Ifges(7,1) * t332 + t324;
t326 = Ifges(5,4) * t332;
t520 = -Ifges(5,2) * t330 + t326;
t265 = Ifges(5,1) * t330 + t326;
t437 = t332 * t499;
t518 = -t437 - t315;
t517 = Ifges(6,3) * t332 - t322;
t516 = Ifges(7,2) * t332 - t324;
t515 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t464 = Ifges(6,5) * t332;
t256 = Ifges(6,3) * t330 + t464;
t465 = Ifges(7,4) * t332;
t258 = Ifges(7,2) * t330 + t465;
t414 = -Ifges(7,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t503 = Ifges(6,6) / 0.2e1;
t390 = t503 + t414;
t513 = t390 * t333 + Ifges(6,6) * t484 + t520 * t488 + (t256 + t258) * t533 + t530 * t485;
t512 = 0.2e1 * m(7);
t511 = 2 * qJD(3);
t510 = m(5) / 0.2e1;
t509 = -m(6) / 0.2e1;
t508 = m(6) / 0.2e1;
t506 = -m(7) / 0.2e1;
t505 = m(7) / 0.2e1;
t504 = -mrSges(5,1) / 0.2e1;
t475 = pkin(8) * t331;
t271 = pkin(3) * t333 + t475;
t335 = -pkin(1) - pkin(7);
t436 = t333 * t335;
t288 = t330 * t436;
t412 = t499 * t333;
t51 = t288 - t412 + (qJ(6) * t331 - t271) * t332;
t502 = t51 / 0.2e1;
t105 = t271 * t332 - t288;
t88 = -pkin(4) * t333 - t105;
t501 = t88 / 0.2e1;
t474 = pkin(8) * t333;
t477 = pkin(3) * t331;
t244 = qJ(2) - t474 + t477;
t440 = t331 * t335;
t101 = t330 * t244 + t332 * t440;
t497 = t101 / 0.2e1;
t196 = t380 * t333;
t496 = t196 / 0.2e1;
t495 = t203 / 0.2e1;
t494 = t228 / 0.2e1;
t442 = t331 * t332;
t294 = mrSges(7,3) * t442;
t453 = t333 * mrSges(7,1);
t230 = t294 - t453;
t493 = -t230 / 0.2e1;
t492 = -t233 / 0.2e1;
t491 = t248 / 0.2e1;
t490 = t251 / 0.2e1;
t489 = t252 / 0.2e1;
t486 = -t332 / 0.2e1;
t245 = -pkin(3) - t429;
t483 = m(6) * t245;
t317 = t333 * qJ(5);
t291 = t330 * t317;
t104 = -t332 * t412 - t291;
t481 = m(7) * t104;
t480 = m(7) * t538;
t479 = m(7) * t330;
t473 = mrSges(6,2) - mrSges(7,3);
t472 = -mrSges(7,2) - mrSges(6,3);
t431 = -t332 * t244 + t330 * t440;
t69 = qJ(6) * t438 - t431;
t50 = -t331 * t499 - t69;
t471 = t50 + t69;
t290 = qJ(6) * t443;
t316 = t331 * qJ(5);
t81 = t101 + t316;
t56 = t290 + t81;
t70 = t101 + t290;
t470 = -t56 + t70;
t466 = Ifges(5,4) * t330;
t106 = t330 * t271 + t332 * t436;
t109 = (-t250 + t335) * t331;
t292 = t332 * t317;
t110 = -t292 + (-t335 + t476) * t333;
t199 = t252 * t331;
t200 = t253 * t331;
t254 = t330 * mrSges(5,1) + t332 * mrSges(5,2);
t201 = t254 * t331;
t226 = t333 * mrSges(7,2) - mrSges(7,3) * t444;
t231 = t333 * mrSges(5,1) + mrSges(5,3) * t442;
t419 = mrSges(6,2) * t442;
t454 = t333 * mrSges(6,1);
t232 = -t419 - t454;
t458 = t331 * mrSges(5,1);
t234 = -mrSges(5,3) * t438 + t458;
t455 = t331 * Ifges(7,5);
t154 = t333 * t521 - t455;
t156 = t331 * Ifges(6,4) + t333 * t523;
t266 = Ifges(5,1) * t332 - t466;
t158 = t331 * Ifges(5,5) + t333 * t266;
t415 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t391 = -Ifges(7,5) / 0.2e1 + t415;
t347 = -t154 / 0.2e1 - t156 / 0.2e1 - t158 / 0.2e1 - t391 * t331;
t348 = Ifges(7,5) * t485 + t391 * t333 + t531 * t484 + (t266 + t521 + t523) * t533;
t296 = Ifges(6,5) * t438;
t148 = Ifges(6,6) * t331 + Ifges(6,3) * t443 + t296;
t297 = Ifges(7,4) * t438;
t150 = Ifges(7,2) * t443 - t331 * Ifges(7,6) + t297;
t152 = t331 * Ifges(5,6) + t333 * t520;
t386 = -t148 / 0.2e1 - t150 / 0.2e1 + t152 / 0.2e1;
t87 = t317 + t106;
t58 = -qJ(6) * t444 + t87;
t79 = (-t335 + t538) * t331;
t80 = t292 + (t335 - t536) * t333;
t83 = -pkin(4) * t331 + t431;
t3 = -t110 * t199 - t80 * t200 + t109 * t202 + t79 * t203 + t56 * t226 + t101 * t227 + t58 * t228 + t106 * t229 + t50 * t230 - t431 * t231 + t83 * t232 + t51 * t233 + t105 * t234 + t88 * t235 + t87 * t236 + t81 * t237 + m(6) * (t109 * t110 + t81 * t87 + t83 * t88) + m(7) * (t50 * t51 + t56 * t58 + t79 * t80) + m(5) * (t101 * t106 - t105 * t431) + (-qJ(2) * mrSges(4,2) + Ifges(4,4) * t331 + t347 * t332 + (-t331 * t390 + t386) * t330) * t331 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t333 + t335 * t201 + t348 * t332 + t513 * t330 + (-Ifges(4,1) + Ifges(4,2) + (-m(5) * t335 + t254) * t335 + t515) * t331) * t333;
t460 = t3 * qJD(1);
t195 = pkin(4) * t438 + t291;
t381 = t332 * mrSges(5,1) - t330 * mrSges(5,2);
t198 = t381 * t333;
t204 = t517 * t333;
t205 = t516 * t333;
t259 = Ifges(5,2) * t332 + t466;
t206 = t333 * t259;
t207 = -Ifges(7,1) * t443 + t297;
t208 = -Ifges(6,1) * t443 + t296;
t209 = t333 * t265;
t295 = Ifges(6,6) * t438;
t432 = -t234 + t235;
t4 = t110 * t196 - t80 * t197 + t195 * t202 + t104 * t203 + t69 * t228 + t70 * t233 + t295 * t488 + t432 * t101 - t433 * t431 + m(6) * (t101 * t83 + t110 * t195 - t431 * t81) + m(7) * (t104 * t80 + t50 * t70 + t56 * t69) + (-t335 * t198 + (-t101 * mrSges(5,3) - t81 * mrSges(6,2) + t56 * mrSges(7,3) + t207 / 0.2e1 + t208 / 0.2e1 - t209 / 0.2e1 + t414 * t331 - t386) * t332 + (-t431 * mrSges(5,3) - t83 * mrSges(6,2) + t50 * mrSges(7,3) + t204 / 0.2e1 + t205 / 0.2e1 + t206 / 0.2e1 + t347) * t330) * t333;
t452 = t4 * qJD(1);
t392 = mrSges(6,3) / 0.2e1 + mrSges(7,2) / 0.2e1 - mrSges(5,2) / 0.2e1;
t416 = mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1;
t394 = t504 - t416;
t341 = t392 * t330 - t394 * t332 + t429 * t508 - t505 * t518;
t369 = t101 * t330 - t332 * t431;
t374 = t81 * t330 - t83 * t332;
t375 = t56 * t330 - t50 * t332;
t400 = t229 / 0.2e1 + t236 / 0.2e1;
t384 = t494 + t400;
t399 = t234 / 0.2e1 - t235 / 0.2e1;
t385 = t492 + t399;
t7 = (t385 * t332 + t384 * t330 + (t369 - t374) * t509 + (t70 * t330 + t69 * t332 - t375) * t506) * t331 + t341 + (t198 / 0.2e1 + t197 / 0.2e1 + t496 + 0.2e1 * t195 * t507 - t481 / 0.2e1 + t427 * t331 * (mrSges(6,2) / 0.2e1 + mrSges(5,3) / 0.2e1 - mrSges(7,3) / 0.2e1)) * t333;
t451 = t7 * qJD(1);
t450 = -t431 + t83;
t449 = t101 - t81;
t448 = -t381 - mrSges(4,1);
t354 = m(7) * t375;
t10 = t331 * mrSges(4,1) + t333 * mrSges(4,2) + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + (-t233 - t432) * t332 + t514 * t330 + t354 + m(6) * t374 + m(5) * t369;
t447 = qJD(1) * t10;
t15 = -t535 * t438 + (t228 + t236) * t331 + m(7) * (t331 * t56 + t438 * t80) + m(6) * (-t110 * t438 + t331 * t81);
t446 = qJD(1) * t15;
t26 = (t330 * t228 - t332 * t233 + t354) * t333;
t445 = qJD(1) * t26;
t441 = t331 * t333;
t388 = (t331 ^ 2 + t333 ^ 2) * t332 * t534;
t40 = (t505 + t508) * t332 + t388;
t434 = t40 * qJD(1);
t430 = t427 * t474;
t426 = qJD(3) * t330;
t425 = qJD(3) * t331;
t422 = m(7) * t438;
t421 = t79 * t505;
t418 = -t480 / 0.2e1;
t411 = t333 * t437;
t410 = -t444 / 0.2e1;
t409 = t444 / 0.2e1;
t407 = -t442 / 0.2e1;
t406 = t442 / 0.2e1;
t404 = t441 / 0.2e1;
t402 = t495 - t202 / 0.2e1;
t401 = -t226 / 0.2e1 - t237 / 0.2e1;
t398 = -t329 / 0.2e1 - t328 / 0.2e1;
t389 = (t380 / 0.2e1 + t491) * t332;
t383 = t517 / 0.2e1 + t516 / 0.2e1 + t259 / 0.2e1;
t261 = Ifges(7,1) * t330 - t465;
t263 = Ifges(6,1) * t330 - t464;
t382 = -t261 / 0.2e1 - t263 / 0.2e1 - t265 / 0.2e1;
t373 = t330 * t88 + t332 * t87;
t29 = 0.4e1 * (m(5) / 0.4e1 + t424) * (-0.1e1 + t427) * t441;
t350 = (-t246 * t332 + t251 * t330) * t505;
t367 = -t105 * t330 + t106 * t332;
t353 = m(5) * t367;
t356 = t110 + t373;
t357 = t330 * t83 + t332 * t81 - t109;
t359 = t330 * t51 + t332 * t58 - t80;
t360 = t330 * t50 + t332 * t56 + t79;
t368 = t101 * t332 + t330 * t431;
t5 = t350 + ((-t227 / 0.2e1 + t401) * t332 + (t493 - t232 / 0.2e1 + t231 / 0.2e1) * t330 + t359 * t506 + t356 * t509 - t353 / 0.2e1 + t402) * t331 + (t200 / 0.2e1 - t201 / 0.2e1 - t199 / 0.2e1 + (-t456 / 0.2e1 - t384) * t332 + (-t458 / 0.2e1 + t385) * t330 + t360 * t506 + t357 * t509 - m(5) * (t368 - 0.2e1 * t440) / 0.2e1) * t333;
t372 = -t5 * qJD(1) + t29 * qJD(2);
t210 = pkin(5) * t332 - t245;
t340 = t402 * t330 + (-t110 * t330 + (-t245 * t333 + t475) * t332) * t508 + (t210 * t438 + t251 * t331 + t330 * t80) * t505;
t364 = m(6) * t501 + m(7) * t502;
t11 = t419 - t294 + (t389 + t416) * t333 + t340 - t364;
t35 = -t210 * t479 + (t483 + t525) * t330;
t371 = qJD(1) * t11 - qJD(3) * t35;
t33 = t197 + (t411 / 0.4e1 + t291 / 0.4e1 - t104 / 0.4e1) * t512;
t43 = (-t439 / 0.4e1 + t536 / 0.4e1 + t538 / 0.4e1) * t512 - t253;
t370 = qJD(1) * t33 + qJD(3) * t43;
t352 = t101 + 0.2e1 * t316;
t342 = -t318 - t319 + t352 * t509 + (t290 + t352) * t506;
t362 = m(6) * t497 + t505 * t70;
t25 = t342 + t362;
t272 = qJ(5) * t500 - t472;
t365 = qJD(1) * t25 - qJD(4) * t272;
t363 = -mrSges(6,2) * pkin(4) + mrSges(7,3) * t499 - Ifges(7,5);
t361 = -qJ(5) * t473 - t530;
t358 = t392 * t332;
t355 = m(6) * t429;
t336 = (t110 * t250 + t195 * t245) * t508 + (t104 * t210 + t246 * t470 + t251 * t471 - t538 * t80) * t505 - pkin(3) * t198 / 0.2e1 + t104 * t491 + t110 * t489 - t195 * t380 / 0.2e1 - t210 * t197 / 0.2e1 - t538 * t495 + t245 * t496 - t246 * t228 / 0.2e1 + t250 * t202 / 0.2e1 + t233 * t490 + t80 * t253 / 0.2e1 + t522 * t331 / 0.4e1;
t337 = (-pkin(4) * t88 + qJ(5) * t87) * t509 + (qJ(5) * t58 - t499 * t51) * t506 + pkin(4) * t232 / 0.2e1 + t105 * t504 + t106 * mrSges(5,2) / 0.2e1 - t499 * t493 + mrSges(7,1) * t502 - t58 * mrSges(7,2) / 0.2e1 - t87 * mrSges(6,3) / 0.2e1 + mrSges(6,1) * t501;
t338 = (t497 - t81 / 0.2e1) * mrSges(6,2) + (-t70 / 0.2e1 + t56 / 0.2e1) * mrSges(7,3) + (t449 * t508 - t400) * pkin(8) + t148 / 0.4e1 + t150 / 0.4e1 - t152 / 0.4e1 + t207 / 0.4e1 + t208 / 0.4e1 - t209 / 0.4e1;
t339 = (-t431 / 0.2e1 + t83 / 0.2e1) * mrSges(6,2) + (-t69 / 0.2e1 - t50 / 0.2e1) * mrSges(7,3) + (t450 * t508 - t399) * pkin(8) + t154 / 0.4e1 + t156 / 0.4e1 + t158 / 0.4e1 - t204 / 0.4e1 - t205 / 0.4e1 - t206 / 0.4e1;
t343 = -t335 * t254 / 0.2e1 + t537 * pkin(8) * t398;
t344 = mrSges(7,3) * t490 + t266 / 0.4e1 + t523 / 0.4e1 + t521 / 0.4e1 - t259 / 0.4e1 - t516 / 0.4e1 - t517 / 0.4e1;
t345 = t246 * mrSges(7,3) / 0.2e1 - t265 / 0.4e1 - t263 / 0.4e1 - t261 / 0.4e1 - t520 / 0.4e1 + t258 / 0.4e1 + t256 / 0.4e1;
t1 = ((-0.3e1 / 0.4e1 * Ifges(7,5) + t415) * t331 + t344 * t333 + t339) * t332 + ((-0.3e1 / 0.4e1 * Ifges(7,6) - 0.3e1 / 0.4e1 * Ifges(5,6) + t503) * t331 + t345 * t333 + t338) * t330 + t336 + (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t343) * t333 + t337 + t401 * qJ(5);
t18 = t292 * t534 + (t482 / 0.2e1 + t480 / 0.2e1 + t254 / 0.2e1 + t489 - t253 / 0.2e1 + t358 + (pkin(4) * t509 - t499 * t505 + t394) * t330) * t333;
t9 = -pkin(3) * t254 - t538 * t248 - t250 * t380 + t524 * t245 + (t253 - t480) * t210 + (-t256 / 0.2e1 + t520 / 0.2e1 - t258 / 0.2e1 - t382) * t332 + (t266 / 0.2e1 + t521 / 0.2e1 + t523 / 0.2e1 - t383) * t330;
t351 = t1 * qJD(1) - t18 * qJD(2) + t9 * qJD(3);
t102 = (0.1e1 / 0.2e1 + t398) * t478;
t346 = m(7) * ((-t246 * t333 - t56) * t332 + (t251 * t333 - t50) * t330);
t16 = (-t319 / 0.2e1 + t494) * t332 + (t457 / 0.2e1 + t233 / 0.2e1) * t330 + t421 - t346 / 0.2e1;
t47 = t528 - t532;
t349 = -qJD(1) * t16 + qJD(2) * t102 + qJD(3) * t47;
t213 = t500 * t442;
t212 = t500 * t443;
t211 = (qJD(1) * t438 + t426) * m(7);
t103 = (t427 + 0.1e1) * t417;
t77 = m(7) * t251 + (m(6) * pkin(8) + t473) * t332;
t57 = t538 * t505 + t418;
t41 = t486 * t500 + t388;
t39 = (t291 + t411) * t505 + t481 / 0.2e1;
t23 = t293 - t342 + t362 - t420;
t19 = t253 * t484 + t394 * t443 + t292 * t505 + (-pkin(4) * t443 + t292) * t508 + (-t505 * t536 + t358 + t418) * t333 + (t254 + t524) * t485;
t17 = t346 / 0.2e1 + t330 * t492 + t228 * t486 + t421 + (-t320 / 0.2e1 + t459 / 0.2e1) * t331;
t12 = -t453 / 0.2e1 - t454 / 0.2e1 + t333 * t389 + t340 + t364;
t8 = (-t195 * t333 + (t330 * t449 + t332 * t450) * t331) * t508 + (t104 * t333 + (t330 * t470 + t332 * t471) * t331) * t505 + t234 * t407 + t341 - t197 * t484 + (t198 + t196) * t485 + t526 * t406 + t404 * t528 + t514 * t410 - t537 * t427 * t441 / 0.2e1;
t6 = -t200 * t484 + (t331 * t359 + t333 * t360) * t505 + t254 * t404 + (t331 * t356 + t333 * t357) * t508 + (t368 * t333 + (t367 - 0.2e1 * t436) * t331) * t510 + t231 * t410 + t350 + (-t201 - t199) * t485 + (t230 + t232) * t409 + (t226 + t527) * t406 + t514 * t438 / 0.2e1 + t535 * t488 + (-t234 / 0.2e1 + t526 / 0.2e1) * t443;
t2 = (t330 * t345 + t332 * t344 + t343) * t333 + (-t455 / 0.4e1 + t339) * t332 + ((-Ifges(5,6) / 0.4e1 - Ifges(7,6) / 0.4e1) * t331 + t338) * t330 + t336 - t337 + Ifges(7,5) * t406 + Ifges(6,6) * t410 + (t226 + t237) * qJ(5) / 0.2e1 + t530 * t409 + t531 * t407 + t515 * t484;
t13 = [qJD(2) * t10 + qJD(3) * t3 + qJD(4) * t4 + qJD(5) * t15 + qJD(6) * t26, qJD(3) * t6 + qJD(4) * t8 + qJD(5) * t41 + t447, t460 + t6 * qJD(2) + t2 * qJD(4) + t12 * qJD(5) + t17 * qJD(6) + (t109 * t483 / 0.2e1 + (t210 * t79 + t246 * t51 + t251 * t58) * t505) * t511 + (t88 * mrSges(6,2) - t105 * mrSges(5,3) - t51 * mrSges(7,3) + t348) * t426 + (-Ifges(4,5) + (-m(5) * pkin(3) + t448) * t335 + t382 * t332 + t383 * t330) * t425 + (-mrSges(4,2) * t436 - Ifges(4,6) * t333 + pkin(3) * t201 - t109 * t380 - t245 * t199 - t210 * t200 + t251 * t226 + t246 * t230 + t79 * t248 + (t87 * mrSges(6,2) + t106 * mrSges(5,3) - t58 * mrSges(7,3) - t513) * t332 + (t527 * t332 + (-t231 + t232) * t330 + t353 + m(6) * t373) * pkin(8)) * qJD(3), t8 * qJD(2) + t2 * qJD(3) + t23 * qJD(5) + t39 * qJD(6) + t452 + (-t101 * mrSges(5,1) - t101 * mrSges(6,1) - t70 * mrSges(7,1) + t431 * mrSges(5,2) + t69 * mrSges(7,2) - t431 * mrSges(6,3) + t295 + 0.2e1 * (qJ(5) * t69 - t499 * t70) * t505 + 0.2e1 * (-pkin(4) * t101 - qJ(5) * t431) * t508 + (t361 * t332 + (-t363 - t531) * t330) * t333) * qJD(4), qJD(2) * t41 + qJD(3) * t12 + qJD(4) * t23 + t446, qJD(3) * t17 + qJD(4) * t39 + t445; -qJD(3) * t5 - qJD(4) * t7 + qJD(5) * t40 - t447, t29 * qJD(3), t19 * qJD(4) + t212 * qJD(5) + t103 * qJD(6) + (t448 + t525) * t425 + (t210 * t417 + (t430 - t477) * t510 + (t245 * t331 + t430) * t508) * t511 + (t532 - mrSges(4,2) + t427 * (mrSges(5,3) + t473)) * qJD(3) * t333 + t372, -t451 + t19 * qJD(3) + t213 * qJD(5) + ((-mrSges(5,1) - mrSges(6,1) - mrSges(7,1)) * t332 + (mrSges(5,2) + t472) * t330 - t355 + m(7) * t518) * qJD(4) * t331, qJD(3) * t212 + qJD(4) * t213 + t434, t103 * qJD(3); qJD(2) * t5 + qJD(4) * t1 + qJD(5) * t11 - qJD(6) * t16 - t460, -qJD(4) * t18 + qJD(6) * t102 - t372, qJD(4) * t9 - qJD(5) * t35 + qJD(6) * t47, t77 * qJD(5) + t57 * qJD(6) + t351 + (m(7) * (-qJ(5) * t246 - t251 * t499) - t251 * mrSges(7,1) - t246 * mrSges(7,2) + t363 * t332 + t361 * t330 + (-t355 - t380 - t381) * pkin(8) + t522) * qJD(4), qJD(4) * t77 + t371, qJD(4) * t57 + t349; qJD(2) * t7 - qJD(3) * t1 - qJD(5) * t25 + qJD(6) * t33 - t452, qJD(3) * t18 + t451, qJD(6) * t43 - t351, t272 * qJD(5), -t365, t370; -qJD(2) * t40 - qJD(3) * t11 + qJD(4) * t25 - qJD(6) * t422 - t446, -t434, -qJD(6) * t479 - t371, t365, 0, -t211; qJD(3) * t16 - qJD(4) * t33 + qJD(5) * t422 - t445, -t102 * qJD(3), -qJD(4) * t43 + qJD(5) * t479 - t349, -t370, t211, 0;];
Cq  = t13;
