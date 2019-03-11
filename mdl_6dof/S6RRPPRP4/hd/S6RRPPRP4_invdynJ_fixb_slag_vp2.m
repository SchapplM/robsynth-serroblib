% Calculate vector of inverse dynamics joint torques for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:02
% EndTime: 2019-03-09 08:37:44
% DurationCPUTime: 29.12s
% Computational Cost: add. (6608->716), mult. (14808->887), div. (0->0), fcn. (10003->8), ass. (0->325)
t530 = Ifges(5,1) + Ifges(4,1);
t273 = cos(qJ(2));
t256 = t273 * qJDD(1);
t270 = sin(qJ(2));
t361 = qJD(1) * qJD(2);
t346 = t270 * t361;
t212 = -t256 + t346;
t428 = t212 / 0.2e1;
t213 = qJDD(1) * t270 + t273 * t361;
t266 = sin(pkin(9));
t267 = cos(pkin(9));
t164 = qJDD(2) * t266 + t213 * t267;
t433 = t164 / 0.2e1;
t528 = Ifges(4,5) + Ifges(5,4);
t536 = t428 * t528 + t530 * t433;
t499 = -Ifges(6,4) + Ifges(7,5);
t535 = t499 + Ifges(7,5);
t373 = qJD(1) * t270;
t353 = t266 * t373;
t362 = t267 * qJD(2);
t203 = t353 - t362;
t351 = t267 * t373;
t204 = qJD(2) * t266 + t351;
t269 = sin(qJ(5));
t272 = cos(qJ(5));
t123 = t203 * t269 + t204 * t272;
t305 = t272 * t203 - t204 * t269;
t250 = pkin(7) * t373;
t215 = -qJD(2) * pkin(2) + qJD(3) + t250;
t90 = t203 * pkin(3) - t204 * qJ(4) + t215;
t69 = -pkin(4) * t203 - t90;
t19 = -pkin(5) * t305 - qJ(6) * t123 + t69;
t118 = Ifges(6,4) * t305;
t372 = qJD(1) * t273;
t241 = qJD(5) + t372;
t402 = Ifges(7,5) * t305;
t498 = Ifges(7,4) + Ifges(6,5);
t500 = Ifges(6,1) + Ifges(7,1);
t488 = t123 * t500 + t241 * t498 + t118 - t402;
t533 = -t69 * mrSges(6,2) + t19 * mrSges(7,3) - t488 / 0.2e1;
t427 = -t241 / 0.2e1;
t438 = -t123 / 0.2e1;
t493 = Ifges(6,3) + Ifges(7,2);
t532 = t493 * t427 + t498 * t438;
t163 = -t267 * qJDD(2) + t213 * t266;
t40 = qJD(5) * t305 + t163 * t269 + t164 * t272;
t447 = t40 / 0.2e1;
t41 = qJD(5) * t123 - t272 * t163 + t164 * t269;
t445 = t41 / 0.2e1;
t531 = m(6) + m(7);
t434 = t163 / 0.2e1;
t206 = qJDD(5) - t212;
t430 = t206 / 0.2e1;
t439 = -t305 / 0.2e1;
t529 = -Ifges(4,4) + Ifges(5,5);
t496 = Ifges(4,6) - Ifges(5,6);
t527 = Ifges(6,6) - Ifges(7,6);
t494 = -Ifges(4,3) - Ifges(5,2);
t350 = t267 * t372;
t352 = t266 * t372;
t161 = t269 * t350 - t272 * t352;
t363 = qJD(5) * t272;
t364 = qJD(5) * t269;
t184 = t266 * t363 - t267 * t364;
t480 = t161 - t184;
t207 = t266 * t269 + t267 * t272;
t175 = t207 * t273;
t162 = qJD(1) * t175;
t183 = t207 * qJD(5);
t479 = t162 + t183;
t404 = Ifges(5,5) * t266;
t407 = Ifges(4,4) * t266;
t526 = t267 * t530 + t404 - t407;
t323 = t267 * mrSges(5,1) + t266 * mrSges(5,3);
t325 = mrSges(4,1) * t267 - mrSges(4,2) * t266;
t525 = -t325 - t323;
t426 = t241 / 0.2e1;
t437 = t123 / 0.2e1;
t440 = t305 / 0.2e1;
t524 = -Ifges(6,4) * t440 - Ifges(7,5) * t439 - t498 * t426 - t500 * t437 + t533;
t271 = sin(qJ(1));
t274 = cos(qJ(1));
t472 = g(1) * t274 + g(2) * t271;
t523 = -mrSges(4,3) - mrSges(5,2) + mrSges(6,3) + mrSges(7,2);
t507 = -t212 / 0.2e1;
t519 = m(4) * qJD(3);
t339 = -t266 * qJ(4) - pkin(2);
t214 = -t267 * pkin(3) + t339;
t366 = qJD(4) * t266;
t251 = pkin(7) * t372;
t377 = qJ(4) * t350 - t251;
t441 = -pkin(3) - pkin(4);
t475 = -t352 * t441 + t366 - t377;
t415 = -mrSges(6,2) + mrSges(7,3);
t462 = m(7) * qJ(6) + t415;
t416 = -mrSges(7,1) - mrSges(6,1);
t518 = m(7) * pkin(5) - t416;
t249 = Ifges(3,4) * t372;
t517 = t266 * (t204 * Ifges(5,5) - Ifges(5,6) * t372 + t203 * Ifges(5,3)) + Ifges(3,1) * t373 + Ifges(3,5) * qJD(2) + t249 + t267 * (t203 * t529 + t204 * t530 - t372 * t528);
t367 = qJD(3) * t270;
t392 = qJDD(1) * pkin(1);
t109 = pkin(2) * t212 - qJ(3) * t213 - qJD(1) * t367 - t392;
t248 = pkin(7) * t256;
t166 = qJDD(2) * qJ(3) + t248 + (qJD(3) - t250) * qJD(2);
t60 = t109 * t267 - t266 * t166;
t308 = qJDD(4) - t60;
t51 = -pkin(3) * t212 + t308;
t516 = mrSges(5,2) * t51 + t434 * t529 + t536;
t117 = Ifges(7,5) * t123;
t45 = t241 * Ifges(7,6) - Ifges(7,3) * t305 + t117;
t405 = Ifges(6,4) * t123;
t48 = Ifges(6,2) * t305 + t241 * Ifges(6,6) + t405;
t514 = t45 / 0.2e1 + t69 * mrSges(6,1) + t19 * mrSges(7,1) - t48 / 0.2e1;
t513 = -Ifges(6,2) * t440 + Ifges(7,3) * t439 - t426 * t527 + t437 * t499 + t514;
t375 = t273 * pkin(2) + t270 * qJ(3);
t492 = -pkin(1) - t375;
t195 = t492 * qJD(1);
t222 = qJD(2) * qJ(3) + t251;
t125 = t195 * t267 - t266 * t222;
t95 = pkin(3) * t372 + qJD(4) - t125;
t63 = pkin(4) * t372 - pkin(8) * t204 + t95;
t126 = t266 * t195 + t267 * t222;
t103 = -qJ(4) * t372 + t126;
t70 = pkin(8) * t203 + t103;
t21 = t269 * t63 + t272 * t70;
t23 = -pkin(8) * t164 + t212 * t441 + t308;
t365 = qJD(4) * t273;
t61 = t266 * t109 + t267 * t166;
t42 = t212 * qJ(4) - qJD(1) * t365 + t61;
t25 = pkin(8) * t163 + t42;
t4 = -qJD(5) * t21 + t23 * t272 - t25 * t269;
t2 = -pkin(5) * t206 + qJDD(6) - t4;
t202 = t213 * pkin(7);
t176 = -qJDD(2) * pkin(2) + qJDD(3) + t202;
t52 = t163 * pkin(3) - t164 * qJ(4) - t204 * qJD(4) + t176;
t30 = pkin(4) * t163 + t52;
t446 = -t41 / 0.2e1;
t5 = -pkin(5) * t41 + qJ(6) * t40 + qJD(6) * t123 + t30;
t512 = -mrSges(6,2) * t30 + mrSges(7,2) * t2 - mrSges(6,3) * t4 + mrSges(7,3) * t5 + Ifges(6,4) * t446 + 0.2e1 * t430 * t498 + t445 * t535 + 0.2e1 * t447 * t500;
t20 = -t269 * t70 + t272 * t63;
t482 = qJD(6) - t20;
t17 = -pkin(5) * t241 + t482;
t18 = qJ(6) * t241 + t21;
t409 = Ifges(3,4) * t270;
t485 = t273 * Ifges(3,2);
t319 = t409 + t485;
t510 = t17 * mrSges(7,1) + t21 * mrSges(6,2) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t319 / 0.2e1 + t527 * t439 + t494 * t372 / 0.2e1 + t528 * t204 / 0.2e1 - t496 * t203 / 0.2e1 - t18 * mrSges(7,3) - t20 * mrSges(6,1) + t532;
t504 = qJD(2) / 0.2e1;
t503 = mrSges(4,1) + mrSges(5,1);
t502 = mrSges(4,2) - mrSges(5,3);
t501 = -mrSges(3,3) + mrSges(2,2);
t26 = mrSges(6,1) * t206 - mrSges(6,3) * t40;
t27 = -t206 * mrSges(7,1) + t40 * mrSges(7,2);
t490 = t27 - t26;
t28 = -mrSges(6,2) * t206 - mrSges(6,3) * t41;
t29 = -mrSges(7,2) * t41 + mrSges(7,3) * t206;
t489 = t28 + t29;
t385 = t266 * t272;
t208 = -t267 * t269 + t385;
t486 = -pkin(5) * t480 + qJ(6) * t479 - qJD(6) * t208 + t475;
t411 = mrSges(6,3) * t305;
t81 = -mrSges(6,2) * t241 + t411;
t84 = mrSges(7,2) * t305 + mrSges(7,3) * t241;
t413 = t81 + t84;
t410 = mrSges(6,3) * t123;
t82 = mrSges(6,1) * t241 - t410;
t83 = -mrSges(7,1) * t241 + mrSges(7,2) * t123;
t412 = -t83 + t82;
t327 = mrSges(3,1) * t273 - mrSges(3,2) * t270;
t484 = -t327 - mrSges(2,1);
t483 = t207 * t270;
t298 = t270 * t208;
t478 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t203 - mrSges(4,2) * t204 - mrSges(3,3) * t373;
t371 = qJD(2) * t270;
t477 = qJ(4) * t371 - t365;
t476 = t206 * t493 + t40 * t498 - t41 * t527;
t201 = -pkin(7) * t346 + t248;
t474 = t201 * t273 + t202 * t270;
t473 = -t266 * t60 + t267 * t61;
t469 = m(4) + m(5) + t531;
t403 = Ifges(5,5) * t267;
t314 = Ifges(5,3) * t266 + t403;
t468 = t203 * (Ifges(5,6) * t270 + t273 * t314) + (t270 * t528 + t273 * t526) * t204;
t384 = t266 * t273;
t234 = pkin(7) * t384;
t261 = t273 * pkin(3);
t421 = pkin(8) * t270;
t104 = pkin(4) * t273 + t234 + t261 + (-t492 - t421) * t267;
t382 = t267 * t273;
t156 = pkin(7) * t382 + t266 * t492;
t138 = -qJ(4) * t273 + t156;
t386 = t266 * t270;
t124 = pkin(8) * t386 + t138;
t395 = t269 * t104 + t272 * t124;
t355 = -pkin(7) * t266 - pkin(3);
t286 = -pkin(8) * t382 + (-pkin(4) + t355) * t270;
t313 = pkin(2) * t270 - qJ(3) * t273;
t181 = qJD(2) * t313 - t367;
t390 = t181 * t267;
t73 = qJD(2) * t286 - t390;
t165 = t266 * t181;
t383 = t267 * t270;
t300 = -pkin(7) * t383 + pkin(8) * t384;
t74 = qJD(2) * t300 + t165 + t477;
t13 = -qJD(5) * t395 - t269 * t74 + t272 * t73;
t463 = t523 * t270;
t458 = Ifges(5,5) * t433 + Ifges(5,6) * t428 - t164 * Ifges(4,4) / 0.2e1 + Ifges(4,6) * t507 - t42 * mrSges(5,2) + (Ifges(5,3) + Ifges(4,2)) * t434;
t406 = Ifges(4,4) * t267;
t318 = -Ifges(4,2) * t266 + t406;
t322 = mrSges(5,1) * t266 - mrSges(5,3) * t267;
t324 = mrSges(4,1) * t266 + mrSges(4,2) * t267;
t455 = -t95 * (-mrSges(5,1) * t270 + mrSges(5,2) * t382) - t273 * t90 * t322 - t215 * t273 * t324 - t126 * (-mrSges(4,2) * t270 - mrSges(4,3) * t384) - t125 * (mrSges(4,1) * t270 - mrSges(4,3) * t382) - t103 * (-mrSges(5,2) * t384 + mrSges(5,3) * t270) + t203 * (Ifges(4,6) * t270 + t273 * t318) / 0.2e1;
t3 = t269 * t23 + t272 * t25 + t63 * t363 - t364 * t70;
t1 = qJ(6) * t206 + qJD(6) * t241 + t3;
t454 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t451 = -mrSges(6,1) * t30 - mrSges(7,1) * t5 + 0.2e1 * Ifges(7,3) * t445 - t40 * Ifges(6,4) / 0.2e1 - t206 * Ifges(6,6) / 0.2e1 + t535 * t447 + (-t527 + Ifges(7,6)) * t430 + (-t446 + t445) * Ifges(6,2);
t435 = -t163 / 0.2e1;
t422 = pkin(7) * t270;
t414 = -pkin(8) + qJ(3);
t211 = t313 * qJD(1);
t387 = t211 * t267;
t85 = qJD(1) * t286 - t387;
t189 = t266 * t211;
t246 = qJ(4) * t373;
t96 = qJD(1) * t300 + t189 + t246;
t34 = t269 * t85 + t272 * t96;
t408 = Ifges(3,4) * t273;
t391 = t176 * t270;
t381 = t270 * t274;
t380 = t271 * t273;
t379 = t273 * t274;
t378 = t274 * t266;
t348 = t273 * t362;
t376 = -qJ(4) * t348 - qJD(4) * t383;
t374 = t274 * pkin(1) + t271 * pkin(7);
t370 = qJD(2) * t273;
t369 = qJD(3) * t266;
t368 = qJD(3) * t267;
t358 = pkin(7) * t371;
t354 = pkin(3) * t266 + pkin(7);
t349 = t266 * t370;
t343 = -t372 / 0.2e1;
t338 = -t361 / 0.2e1;
t337 = t361 / 0.2e1;
t77 = t163 * mrSges(4,1) + t164 * mrSges(4,2);
t108 = -t212 * mrSges(5,1) + t164 * mrSges(5,2);
t76 = t163 * mrSges(5,1) - t164 * mrSges(5,3);
t155 = t267 * t492 - t234;
t190 = t267 * pkin(4) - t214;
t335 = pkin(3) * t382 + qJ(4) * t384 + t375;
t334 = pkin(2) * t379 + qJ(3) * t381 + t374;
t330 = t266 * t441 - pkin(7);
t329 = t355 * t270;
t185 = t266 * t380 + t267 * t274;
t186 = t267 * t380 - t378;
t263 = t274 * pkin(7);
t328 = -t186 * pkin(3) - qJ(4) * t185 + t263;
t326 = mrSges(3,1) * t270 + mrSges(3,2) * t273;
t317 = Ifges(5,4) * t267 + Ifges(5,6) * t266;
t316 = Ifges(3,5) * t273 - Ifges(3,6) * t270;
t315 = Ifges(4,5) * t267 - Ifges(4,6) * t266;
t33 = -t269 * t96 + t272 * t85;
t54 = t104 * t272 - t124 * t269;
t99 = t185 * t272 - t186 * t269;
t306 = t185 * t269 + t186 * t272;
t219 = t414 * t266;
t220 = t414 * t267;
t304 = t272 * t219 - t220 * t269;
t136 = t219 * t269 + t220 * t272;
t141 = -pkin(7) * t351 + t189;
t134 = -t267 * t358 + t165;
t302 = pkin(1) * t326;
t12 = t104 * t363 - t124 * t364 + t269 * t73 + t272 * t74;
t301 = t270 * (Ifges(3,1) * t273 - t409);
t232 = qJ(4) * t383;
t137 = t270 * t330 + t232;
t187 = -t271 * t267 + t273 * t378;
t188 = t271 * t266 + t267 * t379;
t289 = t188 * pkin(3) + t187 * qJ(4) + t334;
t288 = -g(1) * t187 - g(2) * t185 - g(3) * t386;
t94 = t330 * t370 - t376;
t281 = t273 * (Ifges(5,2) * t270 + t273 * t317);
t280 = t273 * (Ifges(4,3) * t270 + t273 * t315);
t223 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t372;
t172 = t269 * t383 - t270 * t385;
t167 = t270 * t354 - t232;
t160 = mrSges(5,1) * t372 + mrSges(5,2) * t204;
t159 = -mrSges(4,1) * t372 - mrSges(4,3) * t204;
t158 = mrSges(4,2) * t372 - mrSges(4,3) * t203;
t157 = -mrSges(5,2) * t203 - mrSges(5,3) * t372;
t145 = pkin(3) * t352 - t377;
t140 = pkin(7) * t353 + t387;
t139 = -t155 + t261;
t133 = t266 * t358 + t390;
t132 = qJD(1) * t329 - t387;
t130 = mrSges(5,1) * t203 - mrSges(5,3) * t204;
t129 = t141 + t246;
t128 = t354 * t370 + t376;
t119 = qJD(2) * t329 - t390;
t113 = t204 * Ifges(4,4) - t203 * Ifges(4,2) - Ifges(4,6) * t372;
t107 = mrSges(4,1) * t212 - mrSges(4,3) * t164;
t106 = -mrSges(4,2) * t212 - mrSges(4,3) * t163;
t105 = -mrSges(5,2) * t163 + mrSges(5,3) * t212;
t102 = t187 * t269 + t188 * t272;
t101 = -t187 * t272 + t188 * t269;
t92 = t134 + t477;
t89 = qJD(2) * t175 + qJD(5) * t298;
t88 = qJD(5) * t483 + t269 * t348 - t272 * t349;
t80 = pkin(5) * t207 - qJ(6) * t208 + t190;
t79 = qJD(5) * t136 + t269 * t368 - t272 * t369;
t78 = qJD(3) * t207 + qJD(5) * t304;
t59 = pkin(5) * t172 - qJ(6) * t483 + t137;
t58 = -mrSges(6,1) * t305 + mrSges(6,2) * t123;
t57 = -mrSges(7,1) * t305 - mrSges(7,3) * t123;
t56 = pkin(5) * t123 - qJ(6) * t305;
t44 = -pkin(5) * t273 - t54;
t43 = qJ(6) * t273 + t395;
t32 = pkin(5) * t373 - t33;
t31 = -qJ(6) * t373 + t34;
t16 = pkin(5) * t88 - qJ(6) * t89 - qJD(6) * t483 + t94;
t15 = t41 * mrSges(6,1) + t40 * mrSges(6,2);
t14 = t41 * mrSges(7,1) - t40 * mrSges(7,3);
t7 = pkin(5) * t371 - t13;
t6 = -qJ(6) * t371 + qJD(6) * t273 + t12;
t8 = [m(7) * (t1 * t43 + t16 * t19 + t17 * t7 + t18 * t6 + t2 * t44 - t5 * t59) + m(5) * (t103 * t92 + t119 * t95 + t128 * t90 + t138 * t42 + t139 * t51 + t167 * t52) + (t17 * mrSges(7,2) - t20 * mrSges(6,3) - t524) * t89 + (t52 * t322 - t337 * t485 + t314 * t434 + t318 * t435 + Ifges(3,5) * qJDD(2) + Ifges(3,1) * t213 + Ifges(3,4) * t507 + t526 * t433 + (t315 + t317) * t428) * t270 + (t316 * t504 - t455) * qJD(2) + (-qJDD(2) * mrSges(3,1) + t77) * t422 + t324 * t391 + qJDD(2) * Ifges(3,6) * t273 + t213 * t408 / 0.2e1 + (t213 * t422 + t474) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t474) + (Ifges(3,4) * t213 - Ifges(3,2) * t212 + t476) * t273 / 0.2e1 + (-t61 * mrSges(4,3) + t458) * t386 - t478 * pkin(7) * t370 + (-mrSges(7,2) * t1 - mrSges(6,3) * t3 + t451) * t172 - t302 * t361 + t327 * t392 - t113 * t349 / 0.2e1 + Ifges(2,3) * qJDD(1) + (-m(3) * t374 - m(4) * t334 - m(5) * t289 - t531 * (t188 * pkin(4) - pkin(8) * t381 + t289) + t484 * t274 + t501 * t271 - t503 * t188 + t502 * t187 - t518 * t102 - t462 * t101 + t523 * t381) * g(2) + (-m(5) * t328 - t531 * (-t186 * pkin(4) + t271 * t421 + t328) - t462 * t99 + t501 * t274 + (-m(3) - m(4)) * t263 + t503 * t186 - t502 * t185 + t518 * t306 + (m(3) * pkin(1) - t469 * t492 - t463 - t484) * t271) * g(1) - (-t163 * t496 + t164 * t528 - t212 * t494) * t273 / 0.2e1 + (t408 * t337 - Ifges(5,6) * t434 - Ifges(4,6) * t435 + t61 * mrSges(4,2) - t42 * mrSges(5,3) - t60 * mrSges(4,1) + t51 * mrSges(5,1) + Ifges(7,6) * t445 + Ifges(6,6) * t446 + pkin(7) * (-qJDD(2) * mrSges(3,2) - mrSges(3,3) * t212) + t498 * t447 - t528 * t433 + t493 * t430 + t494 * t428 + t454) * t273 + t517 * t370 / 0.2e1 + (-mrSges(4,3) * t60 + t516) * t383 - pkin(1) * (mrSges(3,1) * t212 + mrSges(3,2) * t213) + t319 * t507 + t468 * t504 + (-mrSges(7,2) * t18 - mrSges(6,3) * t21 + t513) * t88 + t512 * t483 + m(6) * (t12 * t21 + t13 * t20 - t137 * t30 + t3 * t395 + t4 * t54 + t69 * t94) + t395 * t28 + (-Ifges(6,6) * t440 - Ifges(7,6) * t439 - t493 * t426 - t498 * t437 + t510) * t371 + t167 * t76 + t155 * t107 + t156 * t106 + t92 * t157 + t134 * t158 + t133 * t159 + t119 * t160 + t128 * t130 + t137 * t15 + t138 * t105 + t139 * t108 + t94 * t58 + t12 * t81 + t13 * t82 + t7 * t83 + t6 * t84 + t16 * t57 + t59 * t14 + t44 * t27 + t54 * t26 + t43 * t29 + t301 * t337 - t223 * t358 + m(4) * (t125 * t133 + t126 * t134 + t155 * t60 + t156 * t61 + (t215 * t370 + t391) * pkin(7)) + (t281 + t280) * t338; (-t369 - t140) * t159 + t524 * t183 + (-t468 / 0.2e1 + t455 + (t280 / 0.2e1 + t281 / 0.2e1 + t302 - t301 / 0.2e1) * qJD(1)) * qJD(1) + (Ifges(6,4) * t439 + Ifges(7,5) * t440 + t498 * t427 + t500 * t438 + t533) * t162 + (-t403 + t406) * t433 - (g(1) * t379 + g(2) * t380) * qJ(3) * t469 + t223 * t250 - t412 * t79 + t413 * t78 + t486 * t57 + t489 * t136 + (t368 - t129) * t157 + (-t298 * t462 + t483 * t518 + (pkin(8) * t531 + t523) * t273 + t326 + (-t531 * (t267 * t441 + t339) + m(4) * pkin(2) - t214 * m(5) - t525) * t270) * t472 + (-t103 * t129 - t132 * t95 - t145 * t90 + t214 * t52) * m(5) + (t20 * t479 - t207 * t3 + t21 * t480) * mrSges(6,3) + (t369 - t132) * t160 + t478 * t251 + (-t1 * t207 - t17 * t479 + t18 * t480) * mrSges(7,2) + t473 * mrSges(4,3) + t475 * t58 + (-t366 - t145) * t130 + (t368 - t141) * t158 + t451 * t207 + (m(5) * (qJ(3) * t51 + qJD(3) * t95 - qJD(4) * t90) + (t108 - t107) * qJ(3) - t125 * t519 + t516 + t536) * t266 + (t1 * t136 - t304 * t2 - t5 * t80 + t486 * t19 + (t78 - t31) * t18 + (t79 - t32) * t17) * m(7) - t490 * t304 + (t304 * t4 + t136 * t3 - t190 * t30 + t475 * t69 + (t78 - t34) * t21 + (-t79 - t33) * t20) * m(6) + (-pkin(2) * t176 + t473 * qJ(3) - t125 * t140 - t126 * t141 - t215 * t251) * m(4) + t404 * t434 + t407 * t435 + t113 * t352 / 0.2e1 + (-m(4) * t375 - m(5) * t335 - t327 - t531 * (pkin(4) * t382 + t335 - t421) + t525 * t273 - t518 * t175 - t462 * (t269 * t382 - t272 * t384) + t463) * g(3) + (-Ifges(6,2) * t439 + Ifges(7,3) * t440 - t427 * t527 + t499 * t438 - t514) * t161 + (m(5) * (qJ(3) * t42 + qJD(3) * t103) - Ifges(5,3) * t434 + Ifges(4,2) * t435 + t496 * t428 - t458 + t126 * t519 + (t106 + t105) * qJ(3)) * t267 + (t249 + t517) * t343 + Ifges(3,3) * qJDD(2) - Ifges(3,6) * t212 + Ifges(3,5) * t213 + t214 * t76 + t513 * t184 + t512 * t208 - t201 * mrSges(3,2) - t202 * mrSges(3,1) + t190 * t15 + (-Ifges(3,2) * t343 - Ifges(6,6) * t439 - Ifges(7,6) * t440 - t510 - t532) * t373 - pkin(2) * t77 + t80 * t14 - t34 * t81 - t33 * t82 - t32 * t83 - t31 * t84 - t52 * t323 - t176 * t325 + t316 * t338; t416 * t41 + t415 * t40 + (-t160 + t159) * t204 + (t157 + t158) * t203 + t413 * t305 - t412 * t123 + t76 + t77 + (t123 * t17 + t18 * t305 + t5) * m(7) + (-t123 * t20 + t21 * t305 + t30) * m(6) + (t103 * t203 - t204 * t95 + t52) * m(5) + (t125 * t204 + t126 * t203 + t176) * m(4) + (t273 * g(3) - t270 * t472) * t469; t157 * t372 + (-t57 - t58 + t130) * t204 + (t241 * t413 - t490) * t272 + (-t241 * t412 + t489) * t269 + t108 + (t1 * t269 - t19 * t204 - t2 * t272 + t288 + t241 * (t17 * t269 + t18 * t272)) * m(7) + (-t204 * t69 + t269 * t3 + t272 * t4 + t288 - t241 * (t20 * t269 - t21 * t272)) * m(6) + (t103 * t372 + t204 * t90 + t288 + t51) * m(5); t454 + (-Ifges(6,2) * t123 + t118 + t488) * t439 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t21 + t18 * t482 - t19 * t56) * m(7) + (t410 + t412) * t21 + t48 * t437 + (t101 * t518 - t102 * t462) * g(1) + (t172 * t518 - t462 * t483) * g(3) + (t123 * t18 - t17 * t305) * mrSges(7,2) - t19 * (mrSges(7,1) * t123 - mrSges(7,3) * t305) - t69 * (mrSges(6,1) * t123 + mrSges(6,2) * t305) + (-t123 * t527 + t305 * t498) * t427 + (t305 * t500 + t117 - t405 + t45) * t438 + (Ifges(7,3) * t123 + t402) * t440 + (t411 - t413) * t20 + (-t306 * t462 - t518 * t99) * g(2) + qJD(6) * t84 - t56 * t57 - pkin(5) * t27 + qJ(6) * t29 + t476; t123 * t57 - t241 * t84 + (-g(1) * t101 + g(2) * t99 - g(3) * t172 + t123 * t19 - t18 * t241 + t2) * m(7) + t27;];
tau  = t8;
