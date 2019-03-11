% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRP5_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP5_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP5_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP5_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:45:39
% EndTime: 2019-03-08 21:45:50
% DurationCPUTime: 6.35s
% Computational Cost: add. (6258->524), mult. (14584->708), div. (0->0), fcn. (12945->8), ass. (0->268)
t291 = cos(qJ(3));
t290 = cos(qJ(5));
t275 = t290 * mrSges(7,3);
t287 = sin(qJ(5));
t427 = t287 * mrSges(7,1);
t490 = t275 - t427;
t512 = t490 * t291;
t515 = t512 / 0.2e1;
t288 = sin(qJ(3));
t484 = -t288 * pkin(3) + qJ(4) * t291;
t514 = m(5) * t484;
t286 = sin(pkin(6));
t289 = sin(qJ(2));
t403 = t286 * t289;
t292 = cos(qJ(2));
t401 = t286 * t292;
t366 = -t401 / 0.2e1;
t513 = t366 * t514;
t282 = t287 ^ 2;
t284 = t290 ^ 2;
t384 = t282 + t284;
t508 = -mrSges(6,3) - mrSges(7,2);
t511 = t508 * t384;
t422 = t290 * mrSges(6,2);
t428 = t287 * mrSges(6,1);
t242 = t422 + t428;
t510 = t242 - t490;
t470 = m(7) / 0.2e1;
t509 = 0.2e1 * t470;
t465 = pkin(4) + pkin(8);
t507 = Ifges(7,4) + Ifges(6,5);
t293 = -pkin(3) - pkin(9);
t353 = -qJ(4) * t288 - pkin(2);
t212 = t291 * t293 + t353;
t251 = t465 * t288;
t114 = -t212 * t287 + t251 * t290;
t93 = -pkin(5) * t288 - t114;
t418 = t114 + t93;
t396 = t287 * t291;
t220 = t288 * mrSges(6,1) + mrSges(6,3) * t396;
t221 = -t288 * mrSges(7,1) - mrSges(7,2) * t396;
t506 = -t220 + t221;
t388 = t290 * t291;
t223 = -t288 * mrSges(6,2) - mrSges(6,3) * t388;
t423 = t288 * mrSges(7,3);
t224 = -mrSges(7,2) * t388 + t423;
t493 = t223 + t224;
t505 = t288 ^ 2 + t291 ^ 2;
t504 = t221 / 0.2e1 - t220 / 0.2e1;
t503 = m(6) + m(7);
t502 = -t288 / 0.2e1;
t449 = t288 / 0.2e1;
t501 = -t291 / 0.2e1;
t443 = t291 / 0.2e1;
t500 = mrSges(6,1) + mrSges(7,1);
t497 = -Ifges(5,6) - Ifges(4,4);
t216 = pkin(9) * t288 - t484;
t252 = t465 * t291;
t116 = -t216 * t287 + t252 * t290;
t100 = -pkin(5) * t291 - t116;
t496 = t100 * t470;
t389 = t290 * t212;
t398 = t287 * t251;
t115 = t389 + t398;
t393 = t288 * qJ(6);
t92 = t115 + t393;
t276 = Ifges(7,5) * t290;
t352 = -Ifges(7,1) * t287 + t276;
t177 = t288 * Ifges(7,4) + t291 * t352;
t430 = Ifges(6,4) * t290;
t339 = Ifges(6,1) * t287 + t430;
t312 = t339 * t291;
t179 = t288 * Ifges(6,5) - t312;
t495 = t177 + t179;
t391 = t288 * t290;
t222 = -t291 * mrSges(6,2) + mrSges(6,3) * t391;
t225 = mrSges(7,2) * t391 + t291 * mrSges(7,3);
t494 = t222 + t225;
t429 = Ifges(7,5) * t287;
t249 = t290 * Ifges(7,1) + t429;
t431 = Ifges(6,4) * t287;
t250 = t290 * Ifges(6,1) - t431;
t492 = t249 + t250;
t334 = pkin(5) * t287 - qJ(6) * t290;
t491 = -m(7) * t334 - t510;
t433 = mrSges(7,3) * t287;
t435 = mrSges(7,1) * t290;
t340 = t433 + t435;
t414 = qJ(6) * t287;
t438 = pkin(5) * t290;
t246 = t414 + t438;
t439 = m(7) * t246;
t489 = t340 + t439;
t442 = t293 / 0.2e1;
t488 = t384 * t291 * t442;
t336 = -t287 * Ifges(7,3) - t276;
t487 = t352 - t336;
t245 = t288 * mrSges(4,1) + t291 * mrSges(4,2);
t483 = -t225 / 0.2e1 - t222 / 0.2e1;
t381 = m(6) / 0.4e1 + m(7) / 0.4e1;
t482 = mrSges(5,3) + t510;
t481 = -t250 / 0.2e1 - t249 / 0.2e1;
t274 = m(7) * qJ(6) + mrSges(7,3);
t480 = t493 * t290;
t416 = cos(pkin(6));
t343 = t288 * t403 - t416 * t291;
t479 = (t288 * t343 - t403) * t286;
t373 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t466 = -Ifges(6,5) / 0.2e1;
t375 = t466 - Ifges(7,4) / 0.2e1;
t478 = -t287 * t375 - t290 * t373;
t477 = 0.2e1 * m(7);
t402 = t286 * t291;
t202 = t416 * t288 + t289 * t402;
t476 = 0.2e1 * t202;
t475 = 2 * qJD(3);
t474 = m(5) / 0.2e1;
t473 = -m(6) / 0.2e1;
t472 = m(6) / 0.2e1;
t471 = -m(7) / 0.2e1;
t469 = mrSges(6,1) / 0.2e1;
t468 = -mrSges(6,2) / 0.2e1;
t467 = mrSges(7,3) / 0.2e1;
t143 = t246 * t291 + t252;
t464 = -t143 / 0.2e1;
t463 = t143 / 0.2e1;
t461 = t202 / 0.2e1;
t207 = t340 * t291;
t460 = -t207 / 0.2e1;
t397 = t287 * t288;
t421 = t291 * mrSges(7,1);
t219 = mrSges(7,2) * t397 - t421;
t459 = -t219 / 0.2e1;
t248 = -t287 * Ifges(6,2) + t430;
t454 = t248 / 0.2e1;
t453 = t252 / 0.2e1;
t452 = -t287 / 0.2e1;
t450 = t287 / 0.2e1;
t447 = -t290 / 0.2e1;
t445 = t290 / 0.2e1;
t441 = m(5) * qJ(4);
t236 = qJ(4) + t334;
t440 = m(7) * t236;
t437 = m(7) * qJD(6);
t436 = mrSges(6,1) * t290;
t434 = mrSges(6,2) * t287;
t424 = t288 * mrSges(5,2);
t419 = t92 * t290;
t417 = -t115 + t92;
t395 = t287 * t292;
t118 = t286 * t395 + t290 * t343;
t372 = t290 * t401;
t119 = -t287 * t343 + t372;
t408 = t202 * t290;
t409 = t202 * t287;
t11 = t503 * (t118 * t408 - t119 * t409 - t202 * t343);
t413 = t11 * qJD(1);
t412 = t115 * t290;
t371 = t291 * t401;
t133 = t202 * t371;
t154 = t287 * t403 - t288 * t372;
t155 = (t288 * t395 + t289 * t290) * t286;
t12 = m(5) * (t202 * t402 + t479) * t292 + m(4) * (t479 * t292 + t133) + t503 * (-t118 * t154 - t119 * t155 + t133);
t411 = t12 * qJD(1);
t306 = t418 * t287 + t417 * t290;
t399 = t287 * t221;
t298 = t306 * t471 + t220 * t450 - t399 / 0.2e1 + t493 * t447;
t299 = (t334 * t470 - t275 / 0.2e1 + t427 / 0.2e1 + t422 / 0.2e1 + t428 / 0.2e1) * t288;
t305 = t508 * (-t284 / 0.2e1 - t282 / 0.2e1);
t13 = -t291 * t305 + t298 + t299;
t410 = t13 * qJD(2);
t400 = t287 * t155;
t392 = t288 * t119;
t390 = t290 * t154;
t117 = t290 * t216 + t287 * t252;
t385 = t505 * pkin(8) * t401;
t382 = qJD(3) * t290;
t380 = m(7) * t397;
t378 = t287 * t437;
t377 = t469 + mrSges(7,1) / 0.2e1;
t376 = mrSges(6,2) / 0.2e1 - mrSges(7,3) / 0.2e1;
t374 = Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t363 = -t396 / 0.4e1;
t335 = -Ifges(7,3) * t290 + t429;
t337 = Ifges(6,2) * t290 + t431;
t359 = Ifges(6,6) * t501 + Ifges(7,6) * t443 + t335 * t449 + t337 * t502;
t358 = (-t352 + t339) * t502 + t507 * t501;
t356 = t223 / 0.2e1 + t224 / 0.2e1;
t354 = mrSges(7,2) * qJ(6) - Ifges(7,6);
t350 = -m(5) * pkin(3) - mrSges(4,1) + mrSges(5,2);
t348 = mrSges(7,2) * pkin(5) - t507;
t341 = -t434 + t436;
t142 = (-t246 - t465) * t288;
t205 = t340 * t288;
t206 = t341 * t288;
t208 = t341 * t291;
t218 = t291 * mrSges(6,1) - mrSges(6,3) * t397;
t244 = -t291 * mrSges(5,3) - t424;
t327 = t114 * t290 + t115 * t287;
t331 = t287 * t92 - t290 * t93;
t94 = qJ(6) * t291 + t117;
t294 = -t513 + (-t205 - t206) * t461 - t504 * t408 + t493 * t409 / 0.2e1 + (t244 + t245) * t366 + ((-t251 + t327) * t472 + (t142 + t331) * t470) * t202 + (-t252 * t472 - t143 * t470 + t460 - t208 / 0.2e1) * t343 + (-t117 * t472 - t94 * t470 + t483) * t119 + (t116 * t472 - t496 + t218 / 0.2e1 + t459) * t118;
t323 = -t390 + t400;
t309 = t323 * t293;
t295 = (qJ(4) * t371 + t309) * t473 + (t236 * t371 + t309) * t471 + t513 + t245 * t401 / 0.2e1 - t508 * (-t390 / 0.2e1 + t400 / 0.2e1) + (t482 * t291 + t424) * t366;
t1 = t295 + t294;
t310 = t335 * t291;
t173 = t288 * Ifges(7,6) - t310;
t311 = t337 * t291;
t175 = t288 * Ifges(6,6) - t311;
t237 = -pkin(3) * t291 + t353;
t240 = t291 * mrSges(5,2) - t288 * mrSges(5,3);
t5 = -pkin(2) * t245 + t100 * t221 + t114 * t218 + t115 * t222 + t116 * t220 + t117 * t223 + t142 * t207 - t143 * t205 - t252 * t206 - t251 * t208 + t93 * t219 + t94 * t224 + t92 * t225 - t484 * t240 + (t244 - t514) * t237 + m(6) * (t114 * t116 + t115 * t117 - t251 * t252) + m(7) * (t100 * t93 + t142 * t143 + t92 * t94) + (t358 * t287 + t359 * t290 + (-t478 - t497) * t291) * t291 + (t497 * t288 + (-t173 / 0.2e1 + t175 / 0.2e1 - t373 * t288) * t290 + (t177 / 0.2e1 + t179 / 0.2e1 - t375 * t288) * t287 + (Ifges(5,2) + Ifges(4,1) - Ifges(4,2) - Ifges(5,3) + Ifges(6,3) + Ifges(7,2)) * t291) * t288;
t333 = t1 * qJD(1) + t5 * qJD(2);
t209 = t334 * t291;
t211 = t291 * t242;
t296 = (-t211 / 0.2e1 + t515) * t202 - t504 * t119 + t356 * t118 + (t417 * t118 - t418 * t119 - t202 * t209) * t470 - t508 * t291 * (t118 * t445 + t119 * t452);
t316 = m(7) * (-pkin(5) * t154 + qJ(6) * t155);
t6 = t376 * t155 + t377 * t154 - t316 / 0.2e1 + t296;
t269 = Ifges(6,6) * t396;
t8 = -t252 * t211 - t209 * t207 + t269 * t449 + (t391 * t466 + (-Ifges(7,4) * t290 - Ifges(7,6) * t287) * t449 + t173 * t452 + t327 * mrSges(6,3) + t331 * mrSges(7,2) + (t248 + t336) * t291 * t445 + (t492 * t291 + t175) * t450 + t495 * t447) * t291 + (-m(7) * t209 + t512) * t143 + (m(7) * t93 + t506) * t115 + (m(7) * t92 + t493) * t114;
t332 = t6 * qJD(1) + t8 * qJD(2);
t330 = -t100 * t290 + t287 * t94;
t21 = -t220 * t397 + (m(5) * t237 + t399 + t240 + t480 - m(7) * (-t93 * t287 - t419) - m(6) * (t114 * t287 - t412)) * t288;
t301 = (t473 + t471) * t323;
t325 = t118 * t287 + t119 * t290;
t313 = m(7) * t325;
t314 = m(6) * t325;
t22 = 0.2e1 * (t314 / 0.4e1 + t313 / 0.4e1) * t288 + t301;
t329 = qJD(1) * t22 - qJD(2) * t21;
t34 = t207 * t396 + m(7) * (t143 * t396 + t288 * t92) + t288 * t224;
t38 = (t154 / 0.4e1 + t202 * t363 + t392 / 0.4e1) * t477;
t328 = -qJD(1) * t38 + qJD(2) * t34;
t326 = t116 * t290 + t117 * t287;
t324 = t246 * t143 - t236 * t209;
t120 = (-t490 + t440) * t290;
t300 = (-t290 * t143 + (t236 * t291 + t288 * t293) * t287) * t471 + t207 * t445;
t317 = t496 - t421 / 0.2e1;
t29 = (t288 * mrSges(7,2) + t515) * t287 + t300 + t317;
t322 = qJD(2) * t29 + qJD(3) * t120;
t41 = t423 + (t393 / 0.2e1 + t389 / 0.4e1 + t398 / 0.4e1 - t115 / 0.4e1) * t477;
t321 = qJD(2) * t41 + qJD(5) * t274;
t319 = m(6) * t453 + m(7) * t463;
t318 = -Ifges(6,2) / 0.4e1 - Ifges(7,3) / 0.4e1 + Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1;
t26 = qJ(4) * t341 - t246 * t490 + t337 * t450 + t489 * t236 + (t339 + t248) * t447 + t487 * t445 + (t335 + t492) * t452;
t297 = (-pkin(5) * t100 + qJ(6) * t94) * t470 + pkin(5) * t459 + qJ(6) * t225 / 0.2e1 - t100 * mrSges(7,1) / 0.2e1 + t116 * t469 + t117 * t468 + t94 * t467;
t304 = qJ(4) * t211 / 0.2e1 - t209 * t490 / 0.2e1 - t236 * t512 / 0.2e1 + t246 * t460;
t4 = t324 * t471 + (-t293 * t305 + t374) * t291 + (mrSges(7,3) * t464 + mrSges(6,2) * t453 + t179 / 0.4e1 + t177 / 0.4e1 + (0.3e1 / 0.4e1 * Ifges(7,4) + 0.3e1 / 0.4e1 * Ifges(6,5)) * t288 + (t114 / 0.2e1 + t93 / 0.2e1) * mrSges(7,2) + (t418 * t471 - t504) * t293 + (t276 / 0.4e1 - t248 / 0.4e1 - t336 / 0.4e1 - t318 * t287) * t291) * t287 + (-t252 * mrSges(6,1) / 0.2e1 + t175 / 0.4e1 - t173 / 0.4e1 + mrSges(7,1) * t464 + (0.3e1 / 0.4e1 * Ifges(6,6) - 0.3e1 / 0.4e1 * Ifges(7,6)) * t288 + (-t115 / 0.2e1 + t92 / 0.2e1) * mrSges(7,2) + (t417 * t471 - t356) * t293 + (t250 / 0.4e1 + t249 / 0.4e1 + t318 * t290 + (-Ifges(6,4) + 0.3e1 / 0.4e1 * Ifges(7,5)) * t287) * t291) * t290 + t297 + t304;
t9 = (t438 / 0.4e1 + t414 / 0.4e1 - t246 / 0.4e1) * m(7) * t476;
t308 = t9 * qJD(1) + t4 * qJD(2) - t26 * qJD(3);
t302 = t326 * t473 + t330 * t471;
t17 = (t219 / 0.2e1 - t218 / 0.2e1 + t377 * t291) * t290 + (-t376 * t291 + t483) * t287 + t302 + t319;
t33 = (-t384 + 0.1e1) * t381 * t476;
t79 = t422 - t275 + t440 + mrSges(5,3) + t500 * t287 + (m(6) + m(5)) * qJ(4);
t307 = qJD(1) * t33 + qJD(2) * t17 + qJD(3) * t79;
t226 = (m(7) * t293 - mrSges(7,2)) * t287;
t39 = (t202 * t396 + t154 - t392) * t470;
t37 = t92 * t509 + t224;
t32 = -t450 * t512 - t300 + t317;
t28 = (t474 + t381) * t476 + t503 * t384 * t461;
t20 = (m(5) * t401 + t314 / 0.2e1 + t313 / 0.2e1) * t288 - t301;
t16 = t219 * t447 + t218 * t445 + (m(5) * pkin(8) - t287 * t376 + t290 * t377 + mrSges(5,1)) * t291 - t302 + t319 + t494 * t450;
t14 = -t443 * t511 - t298 + t299;
t10 = (-t434 / 0.2e1 + t436 / 0.2e1 + t435 / 0.2e1 + t433 / 0.2e1 + t439 / 0.2e1) * t202 + (t341 + t489) * t461;
t7 = t316 / 0.2e1 + t296 - t500 * t154 / 0.2e1 + (t468 + t467) * t155;
t3 = -t304 + t340 * t463 + (t293 * t306 + t324) * t470 + t341 * t453 + t374 * t291 + t478 * t288 + t297 + ((-Ifges(6,6) + Ifges(7,6)) * t290 - t507 * t287) * t288 / 0.4e1 - (-t311 + t175) * t290 / 0.4e1 + (-t310 + t173) * t290 / 0.4e1 + (t454 + t336 / 0.4e1) * t396 + t504 * t287 * t293 + t481 * t388 + t487 * t363 + t442 * t480 + t488 * mrSges(6,3) - (-t312 + t495) * t287 / 0.4e1 + (t418 * t452 + t412 / 0.2e1 - t419 / 0.2e1 + t488) * mrSges(7,2);
t2 = -t295 + t294;
t15 = [qJD(2) * t12 + qJD(3) * t11, t2 * qJD(3) + t20 * qJD(4) + t7 * qJD(5) + t39 * qJD(6) + t411 + (t493 * t155 + t506 * t154 + ((-t291 * mrSges(4,1) + t288 * mrSges(4,2) - mrSges(3,1) + t240) * t289 + (-mrSges(3,2) + (t207 + t208) * t291 + (mrSges(4,3) + mrSges(5,1)) * t505) * t292) * t286 + 0.2e1 * (-t114 * t154 + t115 * t155 + t252 * t371) * t472 + (t143 * t371 + t154 * t93 + t155 * t92) * t509 + 0.2e1 * (t237 * t403 + t385) * t474 + m(4) * (-pkin(2) * t403 + t385)) * qJD(2), t413 + t2 * qJD(2) + t28 * qJD(4) + t10 * qJD(5) + ((mrSges(4,2) - t482) * qJD(3) + (-qJ(4) * t472 - t236 * t470 - t441 / 0.2e1) * t475) * t343 + ((t350 + t511) * qJD(3) - t290 * t437 + (t472 + t470) * t475 * t384 * t293) * t202, qJD(2) * t20 + qJD(3) * t28, t7 * qJD(2) + t10 * qJD(3) + (t500 * t119 + (-mrSges(6,2) + mrSges(7,3)) * t118) * qJD(5) + ((pkin(5) * t119 + qJ(6) * t118) * qJD(5) / 0.2e1 - t119 * qJD(6) / 0.2e1) * t477, t39 * qJD(2) + (-t119 * qJD(5) - t202 * t382) * m(7); qJD(3) * t1 + qJD(4) * t22 + qJD(5) * t6 - qJD(6) * t38 - t411, qJD(3) * t5 - qJD(4) * t21 + qJD(5) * t8 + qJD(6) * t34, t16 * qJD(4) + t3 * qJD(5) + t32 * qJD(6) + (-t116 * mrSges(6,3) + t100 * mrSges(7,2) + (t218 - t219) * t293 - t358) * t382 + ((t236 * t142 + t293 * t330) * t470 + (-qJ(4) * t251 + t293 * t326) * t472) * t475 + t333 + (-qJ(4) * t206 - t142 * t490 - t236 * t205 - t251 * t242 + (-t94 * mrSges(7,2) - t117 * mrSges(6,3) + t494 * t293 + t359) * t287 + (-pkin(3) * mrSges(5,1) + pkin(8) * t350 + t287 * t373 - t290 * t375 - Ifges(5,4) + Ifges(4,5)) * t291 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + (t336 / 0.2e1 + t454) * t290 - t481 * t287 + (mrSges(4,2) - mrSges(5,3) - t441) * pkin(8)) * t288) * qJD(3), qJD(3) * t16 + qJD(5) * t14 + t329, t3 * qJD(3) + t14 * qJD(4) + t37 * qJD(6) + t332 + (t269 + (-m(7) * pkin(5) - t500) * t115 + (-mrSges(6,2) + t274) * t114 + (t287 * t354 + t290 * t348) * t291) * qJD(5), qJD(3) * t32 + qJD(5) * t37 + t328; -qJD(2) * t1 + qJD(4) * t33 - qJD(5) * t9 - t413, qJD(4) * t17 - qJD(5) * t4 - qJD(6) * t29 - t333, qJD(4) * t79 + qJD(5) * t26 - qJD(6) * t120, t307, t226 * qJD(6) - t308 + ((-Ifges(6,6) - t354) * t290 + t348 * t287 + t491 * t293) * qJD(5), qJD(5) * t226 - t322; -qJD(2) * t22 - qJD(3) * t33, -qJD(3) * t17 - qJD(5) * t13 + t288 * t378 - t329, -t307, 0, t491 * qJD(5) + t378 - t410 (qJD(2) * t288 + qJD(5)) * t287 * m(7); -qJD(2) * t6 + qJD(3) * t9, qJD(3) * t4 + qJD(4) * t13 + qJD(6) * t41 - t332, t308, t410, t274 * qJD(6), t321; t38 * qJD(2), qJD(3) * t29 - qJD(4) * t380 - qJD(5) * t41 - t328, t322, -qJD(2) * t380, -t321, 0;];
Cq  = t15;
