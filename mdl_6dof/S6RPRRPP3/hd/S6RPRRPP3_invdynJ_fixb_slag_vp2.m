% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:08
% EndTime: 2019-03-09 04:35:40
% DurationCPUTime: 20.46s
% Computational Cost: add. (4792->680), mult. (9755->810), div. (0->0), fcn. (5730->10), ass. (0->297)
t462 = -Ifges(7,5) - Ifges(5,5);
t482 = Ifges(6,4) + t462;
t463 = -Ifges(6,5) - Ifges(7,4);
t481 = -Ifges(5,6) - t463;
t485 = Ifges(5,1) + Ifges(7,3);
t484 = Ifges(7,2) + Ifges(6,3);
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t205 = sin(qJ(3));
t208 = cos(qJ(3));
t336 = qJD(1) * qJD(3);
t307 = t208 * t336;
t232 = qJDD(1) * t205 + t307;
t339 = qJD(4) * t205;
t312 = t204 * t339;
t337 = t207 * qJD(3);
t76 = qJD(1) * t312 - qJD(4) * t337 - t204 * qJDD(3) - t207 * t232;
t417 = -t76 / 0.2e1;
t345 = qJD(3) * t204;
t347 = qJD(1) * t205;
t144 = t207 * t347 + t345;
t341 = qJD(4) * t144;
t77 = -t207 * qJDD(3) + t204 * t232 + t341;
t414 = t77 / 0.2e1;
t143 = t204 * t347 - t337;
t411 = -t143 / 0.2e1;
t408 = t144 / 0.2e1;
t346 = qJD(1) * t208;
t176 = -qJD(4) + t346;
t406 = t176 / 0.2e1;
t487 = qJD(3) / 0.2e1;
t486 = -mrSges(6,1) - mrSges(5,3);
t373 = Ifges(7,6) * t204;
t378 = Ifges(5,4) * t204;
t443 = t207 * t485 + t373 - t378;
t372 = Ifges(7,6) * t207;
t375 = Ifges(6,6) * t207;
t441 = t204 * t484 + t372 - t375;
t483 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t480 = -m(7) * pkin(5) - mrSges(7,1) + t486;
t426 = t481 * t204 - t482 * t207;
t266 = t207 * mrSges(6,2) - t204 * mrSges(6,3);
t268 = mrSges(5,1) * t207 - mrSges(5,2) * t204;
t369 = t204 * mrSges(7,2);
t479 = t268 + t369 - t266;
t198 = t208 * pkin(3);
t196 = t205 * pkin(8);
t315 = -pkin(2) - t196;
t241 = t315 - t198;
t202 = cos(pkin(9));
t402 = pkin(1) * t202;
t225 = (t241 - t402) * qJD(1);
t201 = sin(pkin(9));
t182 = pkin(1) * t201 + pkin(7);
t156 = t182 * qJD(1);
t344 = qJD(3) * t205;
t477 = qJD(2) * qJD(3) + t182 * qJDD(1);
t55 = t205 * qJDD(2) - t156 * t344 + t208 * t477;
t478 = qJDD(3) * pkin(8) + qJD(4) * t225 + t55;
t115 = t205 * qJD(2) + t208 * t156;
t103 = qJD(3) * pkin(8) + t115;
t340 = qJD(4) * t204;
t151 = qJDD(1) * t208 - t205 * t336;
t75 = -pkin(8) * t307 - t151 * pkin(3) + (t315 - t402) * qJDD(1);
t6 = -t103 * t340 + t204 * t75 + t207 * t478;
t338 = qJD(4) * t207;
t7 = -t103 * t338 - t204 * t478 + t207 * t75;
t272 = -t204 * t7 + t207 * t6;
t42 = t103 * t204 - t207 * t225;
t43 = t207 * t103 + t204 * t225;
t476 = t42 * t338 - t43 * t340 + t272;
t138 = qJDD(4) - t151;
t4 = -qJ(5) * t138 + qJD(5) * t176 - t6;
t237 = qJDD(5) - t7;
t5 = -pkin(4) * t138 + t237;
t273 = t204 * t5 - t207 * t4;
t30 = pkin(4) * t176 + qJD(5) + t42;
t31 = t176 * qJ(5) - t43;
t475 = t30 * t338 + t31 * t340 + t273;
t474 = -qJD(5) * t204 - t115 + (-t204 * t346 + t340) * pkin(4);
t343 = qJD(3) * t208;
t473 = t204 * t343 + t205 * t338;
t472 = -t143 * pkin(5) + qJD(6);
t240 = pkin(5) * t144 + t42;
t471 = t240 + qJD(5);
t415 = -t77 / 0.2e1;
t470 = (-Ifges(5,4) + Ifges(7,6)) * t414 + Ifges(6,6) * t415 + (t485 + Ifges(6,2)) * t417 + (-t462 / 0.2e1 - Ifges(6,4) / 0.2e1) * t138;
t187 = Ifges(4,4) * t346;
t135 = Ifges(5,4) * t143;
t374 = Ifges(7,6) * t143;
t455 = t144 * t485 + t176 * t462 - t135 + t374;
t133 = Ifges(7,6) * t144;
t370 = t144 * Ifges(6,6);
t456 = t143 * t484 + t176 * t463 + t133 - t370;
t469 = Ifges(4,1) * t347 + Ifges(4,5) * qJD(3) + t204 * t456 + t207 * t455 + t187;
t380 = Ifges(4,4) * t205;
t452 = t208 * Ifges(4,2);
t261 = t380 + t452;
t468 = t483 * t406 + t482 * t408 + t481 * t411 + Ifges(4,6) * t487 + qJD(1) * t261 / 0.2e1;
t467 = -m(6) - m(7);
t466 = t151 / 0.2e1;
t465 = t232 / 0.2e1;
t464 = -mrSges(4,3) + mrSges(3,2);
t461 = t484 * t77 + (Ifges(6,6) - Ifges(7,6)) * t76 - t463 * t138;
t36 = mrSges(6,1) * t77 - mrSges(6,3) * t138;
t37 = -t77 * mrSges(7,1) + t138 * mrSges(7,2);
t459 = -t36 + t37;
t33 = mrSges(5,1) * t138 + mrSges(5,3) * t76;
t38 = -t76 * mrSges(6,1) + t138 * mrSges(6,2);
t458 = t38 - t33;
t359 = qJ(5) * t207;
t242 = qJ(6) * t204 - t359;
t233 = t242 * t208;
t457 = -qJD(1) * t233 + qJD(4) * t242 - qJD(6) * t207 + t474;
t386 = mrSges(6,1) * t143;
t97 = mrSges(6,3) * t176 + t386;
t384 = mrSges(7,1) * t143;
t98 = -mrSges(7,2) * t176 - t384;
t387 = -t97 + t98;
t382 = mrSges(5,3) * t143;
t94 = mrSges(5,2) * t176 - t382;
t454 = t97 - t94;
t381 = mrSges(5,3) * t144;
t95 = -mrSges(5,1) * t176 - t381;
t385 = mrSges(6,1) * t144;
t99 = -mrSges(6,2) * t176 + t385;
t453 = t99 - t95;
t20 = mrSges(5,1) * t77 - mrSges(5,2) * t76;
t451 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t232 + t20;
t270 = mrSges(4,1) * t208 - mrSges(4,2) * t205;
t450 = -t270 - mrSges(3,1);
t413 = pkin(5) + pkin(8);
t163 = t413 * t207;
t353 = t207 * t208;
t331 = pkin(5) * t353;
t389 = pkin(4) + qJ(6);
t114 = qJD(2) * t208 - t205 * t156;
t397 = pkin(8) * t208;
t274 = pkin(3) * t205 - t397;
t147 = t274 * qJD(1);
t66 = -t204 * t114 + t147 * t207;
t449 = -(-t205 * t389 + t331) * qJD(1) + t66 + qJD(4) * t163;
t355 = t204 * t208;
t332 = pkin(5) * t355;
t67 = t207 * t114 + t204 * t147;
t448 = -(qJ(5) * t205 - t332) * qJD(1) - t67 - t413 * t340;
t447 = -qJ(5) * t338 + t346 * t359 + t474;
t325 = mrSges(4,3) * t347;
t446 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t143 + mrSges(5,2) * t144 + t325;
t152 = t182 * t343;
t445 = -t205 * t462 + t208 * t443;
t444 = -t205 * t463 + t208 * t441;
t376 = Ifges(6,6) * t204;
t442 = t207 * t484 - t373 + t376;
t440 = -t43 - t472;
t134 = Ifges(6,6) * t143;
t61 = -Ifges(6,4) * t176 - Ifges(6,2) * t144 + t134;
t371 = t144 * Ifges(5,4);
t62 = -Ifges(5,2) * t143 - Ifges(5,6) * t176 + t371;
t437 = t204 * t62 + t207 * t61;
t56 = qJDD(2) * t208 - t156 * t343 - t205 * t477;
t436 = -t205 * t56 + t208 * t55;
t435 = t486 * t205;
t434 = qJD(3) * t114;
t433 = qJD(3) * t115;
t432 = -m(5) + t467;
t431 = -mrSges(7,2) - mrSges(6,3) + mrSges(5,2);
t428 = t205 * t483 + t208 * t426;
t102 = -qJD(3) * pkin(3) - t114;
t226 = -qJ(5) * t144 + t102;
t25 = t143 * t389 + t226;
t364 = t207 * mrSges(7,2);
t264 = mrSges(7,3) * t204 - t364;
t265 = mrSges(6,2) * t204 + mrSges(6,3) * t207;
t267 = mrSges(5,1) * t204 + mrSges(5,2) * t207;
t41 = pkin(4) * t143 + t226;
t427 = -t102 * t267 - t25 * t264 + t41 * t265;
t425 = t138 * t483 + t481 * t77 + t482 * t76;
t308 = m(7) * qJ(6) + mrSges(7,3);
t289 = -qJ(5) * t204 - pkin(3);
t153 = -pkin(4) * t207 + t289;
t269 = mrSges(4,1) * t205 + mrSges(4,2) * t208;
t422 = t269 + t480 * t208 + (-m(7) * t289 - (-m(7) * t389 - mrSges(7,3)) * t207 - m(6) * t153 + m(5) * pkin(3) + t479) * t205;
t1 = -pkin(5) * t76 + qJD(6) * t176 - t138 * t389 + t237;
t17 = t176 * t389 + t471;
t2 = -pkin(5) * t77 + qJDD(6) - t4;
t22 = -t31 + t472;
t421 = t1 * t204 + t17 * t338 + t2 * t207 - t22 * t340;
t420 = mrSges(5,1) - mrSges(6,2) + t308;
t419 = -t7 * mrSges(5,1) + t6 * mrSges(5,2) - t5 * mrSges(6,2) - t2 * mrSges(7,2) + t4 * mrSges(6,3) + t1 * mrSges(7,3);
t416 = t76 / 0.2e1;
t412 = t138 / 0.2e1;
t410 = t143 / 0.2e1;
t409 = -t144 / 0.2e1;
t407 = -t176 / 0.2e1;
t206 = sin(qJ(1));
t401 = pkin(1) * t206;
t400 = pkin(5) * t205;
t209 = cos(qJ(1));
t199 = t209 * pkin(1);
t84 = -mrSges(7,2) * t144 + mrSges(7,3) * t143;
t87 = -mrSges(6,2) * t143 - mrSges(6,3) * t144;
t388 = t84 + t87;
t383 = mrSges(7,1) * t144;
t379 = Ifges(4,4) * t208;
t377 = Ifges(5,4) * t207;
t360 = qJ(5) * t143;
t200 = qJ(1) + pkin(9);
t190 = cos(t200);
t358 = t190 * t205;
t357 = t190 * t208;
t356 = t204 * t205;
t161 = t205 * t182;
t354 = t205 * t207;
t183 = -pkin(2) - t402;
t350 = t198 + t196;
t132 = t183 - t350;
t150 = t274 * qJD(3);
t352 = t132 * t338 + t204 * t150;
t146 = t182 * t353;
t83 = t204 * t132 + t146;
t351 = -pkin(4) * t356 + qJ(5) * t354;
t334 = t94 + t387;
t96 = mrSges(7,3) * t176 + t383;
t333 = t96 + t453;
t324 = mrSges(4,3) * t346;
t145 = t182 * t355;
t92 = t161 - t351;
t189 = sin(t200);
t316 = t190 * pkin(2) + t189 * pkin(7) + t199;
t314 = t182 * t344;
t294 = -t340 / 0.2e1;
t293 = t340 / 0.2e1;
t292 = -t338 / 0.2e1;
t291 = t338 / 0.2e1;
t35 = -t76 * mrSges(7,1) - t138 * mrSges(7,3);
t290 = t190 * pkin(7) - t401;
t288 = -t182 * t204 - pkin(4);
t287 = t336 / 0.2e1;
t284 = t182 * t207 - qJ(5);
t82 = t132 * t207 - t145;
t282 = pkin(4) * t473 + qJ(5) * t312 + t152;
t281 = pkin(4) * t353 + qJ(5) * t355 + t350;
t79 = qJ(5) * t208 - t83;
t271 = qJD(4) * t146 + t132 * t340 - t150 * t207;
t262 = Ifges(5,1) * t204 + t377;
t260 = -Ifges(5,2) * t204 + t377;
t259 = Ifges(5,2) * t207 + t378;
t257 = Ifges(6,4) * t204 + Ifges(6,5) * t207;
t256 = Ifges(7,4) * t207 - Ifges(7,5) * t204;
t254 = Ifges(4,5) * t208 - Ifges(4,6) * t205;
t252 = Ifges(5,5) * t204 + Ifges(5,6) * t207;
t251 = -Ifges(6,2) * t207 + t376;
t250 = Ifges(6,2) * t204 + t375;
t245 = -Ifges(7,3) * t204 + t372;
t243 = pkin(3) * t357 + pkin(8) * t358 + t316;
t235 = t183 * qJD(1) * t269;
t234 = t205 * (Ifges(4,1) * t208 - t380);
t108 = t189 * t355 + t190 * t207;
t110 = -t189 * t207 + t190 * t355;
t228 = -g(1) * t110 - g(2) * t108 - g(3) * t356;
t52 = -qJDD(3) * pkin(3) - t56;
t222 = Ifges(6,4) * t205 + t208 * t251;
t217 = Ifges(5,6) * t205 + t208 * t260;
t213 = qJ(5) * t76 - qJD(5) * t144 + t52;
t197 = t208 * pkin(4);
t162 = t413 * t204;
t159 = -qJD(3) * mrSges(4,2) + t324;
t155 = t183 * qJDD(1);
t131 = t267 * t205;
t130 = -t207 * t389 + t289;
t117 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t151;
t111 = t189 * t204 + t190 * t353;
t109 = t189 * t353 - t190 * t204;
t85 = pkin(4) * t144 + t360;
t81 = qJ(6) * t356 + t92;
t80 = t197 - t82;
t53 = -pkin(5) * t356 - t79;
t50 = -pkin(4) * t347 - t66;
t49 = -qJ(5) * t347 - t67;
t47 = t144 * t389 + t360;
t44 = qJ(6) * t208 + t145 + t197 + (-t132 + t400) * t207;
t39 = (-qJ(5) * t343 - qJD(5) * t205) * t207 + t282;
t34 = -mrSges(5,2) * t138 - mrSges(5,3) * t77;
t29 = t204 * t314 - t271;
t28 = (-t205 * t337 - t208 * t340) * t182 + t352;
t24 = t288 * t344 + t271;
t23 = (t182 * t340 + qJD(5)) * t208 + t284 * t344 - t352;
t21 = mrSges(7,2) * t76 + mrSges(7,3) * t77;
t19 = -mrSges(6,2) * t77 + mrSges(6,3) * t76;
t18 = qJD(3) * t233 + (qJD(6) * t204 + (qJ(6) * qJD(4) - qJD(5)) * t207) * t205 + t282;
t16 = -qJD(5) * t208 + (-pkin(5) * t354 - t145) * qJD(4) + (-t205 * t284 - t332) * qJD(3) + t352;
t14 = -t76 * Ifges(5,4) - t77 * Ifges(5,2) + t138 * Ifges(5,6);
t9 = -pkin(5) * t312 + qJD(6) * t208 + (t331 + (-qJ(6) + t288) * t205) * qJD(3) + t271;
t8 = pkin(4) * t77 + t213;
t3 = qJD(6) * t143 + t389 * t77 + t213;
t10 = [(-t42 * mrSges(5,1) - t43 * mrSges(5,2) + t30 * mrSges(6,2) + t22 * mrSges(7,2) - t115 * mrSges(4,3) - t31 * mrSges(6,3) - t17 * mrSges(7,3) - t468) * t344 + ((-t252 + t256 + t257) * t339 + t428 * qJD(3)) * t407 + (t461 / 0.2e1 - t14 / 0.2e1 - t6 * mrSges(5,3) - t2 * mrSges(7,1) + t4 * mrSges(6,1)) * t356 + (t217 * t411 + t222 * t409 + t254 * t487 + t445 * t408 + t444 * t410 + t235) * qJD(3) + t379 * t465 + t261 * t466 + (-m(3) * t199 - m(4) * t316 - m(5) * t243 - mrSges(2,1) * t209 + mrSges(2,2) * t206 + t467 * (t111 * pkin(4) + qJ(5) * t110 + t243) + t450 * t190 + t464 * t189 + t480 * t358 - t420 * t111 + t431 * t110) * g(2) + (-t425 / 0.2e1 - t483 * t412 + t379 * t287 - Ifges(5,6) * t415 - Ifges(6,4) * t416 + Ifges(4,4) * t465 + Ifges(4,2) * t466 + t462 * t417 + t463 * t414 + (-m(4) * t434 + t117) * t182 + t419 + Ifges(4,6) * qJDD(3)) * t208 + (mrSges(6,1) * t30 + mrSges(7,1) * t17 + mrSges(5,2) * t102 - mrSges(7,2) * t25 + mrSges(5,3) * t42 - mrSges(6,3) * t41) * (t208 * t337 - t312) - t155 * t270 + (Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(3,1) * t202 - 0.2e1 * mrSges(3,2) * t201 + m(3) * (t201 ^ 2 + t202 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + t436 * mrSges(4,3) - t159 * t314 + (m(3) * t401 + mrSges(2,1) * t206 + mrSges(2,2) * t209 + (-m(5) - m(4)) * t290 + t467 * (-t109 * pkin(4) - qJ(5) * t108 + t290) + t464 * t190 + t420 * t109 - t431 * t108 + (m(4) * pkin(2) - m(7) * (-pkin(2) - t198) + (-m(5) - m(6)) * t241 - t435 - t450) * t189) * g(1) + (t426 * t412 + t455 * t294 + t456 * t291 - t287 * t452 + t3 * t264 - t8 * t265 + Ifges(4,1) * t232 + t62 * t292 + t61 * t293 + Ifges(4,4) * t466 + t260 * t415 + t251 * t416 + t443 * t417 + t441 * t414 + (-m(4) * t433 + m(5) * t52) * t182 + (m(7) * t413 + mrSges(7,1)) * t189 * g(1) + Ifges(4,5) * qJDD(3)) * t205 + t53 * t37 + t44 * t35 + (-t114 * mrSges(4,3) - t437 / 0.2e1 + t469 / 0.2e1) * t343 + m(7) * (t1 * t44 + t16 * t22 + t17 * t9 + t18 * t25 + t2 * t53 + t3 * t81) + m(6) * (t23 * t31 + t24 * t30 + t39 * t41 + t4 * t79 + t5 * t80 + t8 * t92) + (mrSges(5,1) * t102 + mrSges(6,1) * t31 - mrSges(7,1) * t22 - mrSges(6,2) * t41 - mrSges(5,3) * t43 + mrSges(7,3) * t25) * t473 + t183 * (-t151 * mrSges(4,1) + mrSges(4,2) * t232) + m(5) * (t102 * t152 + t28 * t43 - t29 * t42 + t6 * t83 + t7 * t82) + m(4) * (t155 * t183 + t436 * t182) + (t250 * t409 - t259 * t411 + t442 * t410 + (t245 - t262) * t408) * t339 + (mrSges(6,1) * t5 + mrSges(7,1) * t1 - mrSges(5,3) * t7 + t470) * t354 + t79 * t36 + t80 * t38 + t81 * t21 + t82 * t33 + t83 * t34 + t18 * t84 + t39 * t87 + t92 * t19 + t28 * t94 + t29 * t95 + t9 * t96 + t23 * t97 + t16 * t98 + t24 * t99 + t446 * t152 + t451 * t161 + t52 * t131 + t234 * t287; m(3) * qJDD(2) + (-m(3) - m(4) + t432) * g(3) + (-t19 - t21 + (t204 * t333 + t207 * t334 + t159) * qJD(3) + m(4) * (t56 + t433) + m(6) * (t30 * t345 - t31 * t337 - t8) + m(7) * (t17 * t345 + t22 * t337 - t3) + m(5) * (t337 * t43 + t345 * t42 - t52) - t451) * t208 + (t117 + (t34 + t459) * t207 + (t35 + t458) * t204 + (t388 + t446) * qJD(3) + (-t204 * t334 + t207 * t333) * qJD(4) + m(4) * (t55 - t434) + m(6) * (qJD(3) * t41 + t475) + m(7) * (qJD(3) * t25 + t421) + m(5) * (qJD(3) * t102 + t476)) * t205; t468 * t347 + (t324 - t159) * t114 + (-m(7) * (t281 + t400) - m(5) * t350 - m(6) * t281 - t270 + (-t207 * t308 - t479) * t208 + t435) * g(3) + (-t251 / 0.2e1 + t443 / 0.2e1) * t341 + ((t217 / 0.2e1 - t444 / 0.2e1) * t143 - t43 * (-mrSges(5,2) * t205 - mrSges(5,3) * t355) - t22 * (-mrSges(7,1) * t355 + mrSges(7,2) * t205) - t31 * (mrSges(6,1) * t355 - mrSges(6,3) * t205) + t42 * (mrSges(5,1) * t205 - mrSges(5,3) * t353) - t17 * (mrSges(7,1) * t353 - mrSges(7,3) * t205) - t30 * (mrSges(6,1) * t353 + mrSges(6,2) * t205) - t235 + t428 * t406 + (t222 / 0.2e1 - t445 / 0.2e1) * t144 - t234 * qJD(1) / 0.2e1) * qJD(1) + (t426 * t407 - t427) * qJD(4) + (t262 + t250) * t417 + t3 * (-t207 * mrSges(7,3) - t369) - t55 * mrSges(4,2) + t56 * mrSges(4,1) - pkin(3) * t20 + (-t260 / 0.2e1 + t441 / 0.2e1) * qJD(4) * t143 - t254 * t336 / 0.2e1 - t461 * t207 / 0.2e1 + t455 * t291 + t456 * t293 + t457 * t84 + (t1 * t162 + t130 * t3 + t163 * t2 + t17 * t449 + t22 * t448 + t25 * t457) * m(7) + t475 * mrSges(6,1) + t476 * mrSges(5,3) + (t427 + t437 / 0.2e1) * t346 + Ifges(4,5) * t232 + (-t256 / 0.2e1 - t257 / 0.2e1) * t138 + (-g(3) * t205 + t421) * mrSges(7,1) + t8 * t266 - t52 * t268 + (-pkin(3) * t52 + t42 * t66 - t43 * t67) * m(5) + (t153 * t8 - t30 * t50 - t31 * t49 + t447 * t41) * m(6) + (((-t204 * t43 + t207 * t42) * qJD(4) + t272) * m(5) + t357 * t432 * g(1) + t453 * t338 + t454 * t340 + t458 * t204 + (t34 - t36) * t207 + ((t204 * t31 + t207 * t30) * qJD(4) + t273) * m(6)) * pkin(8) + (t397 * t432 + t422) * g(2) * t189 + (t259 + t442) * t415 + Ifges(4,3) * qJDD(3) - (-Ifges(4,2) * t347 + t187 + t469) * t346 / 0.2e1 + t470 * t204 - t67 * t94 - t66 * t95 - t49 * t97 - t50 * t99 + (-m(5) * t102 + t325 - t446) * t115 + t447 * t87 + t448 * t98 + t449 * t96 + t130 * t21 + t61 * t292 + t62 * t294 + t252 * t412 + t245 * t416 + Ifges(4,6) * t151 + t153 * t19 + t162 * t35 + t163 * t37 + t207 * t14 / 0.2e1 + t190 * t422 * g(1); (-t143 * t485 + t133 - t371 + t456) * t409 + t240 * t98 - t389 * t35 - t419 - pkin(4) * t38 + (t144 * t484 + t134 - t374 + t61) * t411 + t425 + (-m(6) * t30 + t381 - t453) * t43 + (-m(6) * t31 + t382 - t454) * t42 + t387 * qJD(5) + (-Ifges(5,2) * t144 - t135 + t455) * t410 + t459 * qJ(5) + (Ifges(6,2) * t143 + t370 + t62) * t408 + (t131 + t467 * t351 + (t204 * t308 - t265 - t364) * t205) * g(3) + (t467 * (-t110 * pkin(4) + qJ(5) * t111) + t431 * t111 + t420 * t110) * g(1) + (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t31 - t41 * t85) * m(6) + (t482 * t143 + t481 * t144) * t406 + (t467 * (-t108 * pkin(4) + qJ(5) * t109) + t431 * t109 + t420 * t108) * g(2) - t31 * t385 + t440 * t96 + t22 * t383 + t17 * t384 + (qJ(5) * t2 - t1 * t389 + t440 * t17 + t22 * t471 - t25 * t47) * m(7) - t47 * t84 - t85 * t87 - t41 * (-mrSges(6,2) * t144 + mrSges(6,3) * t143) - t25 * (mrSges(7,2) * t143 + mrSges(7,3) * t144) - t102 * (mrSges(5,1) * t144 - mrSges(5,2) * t143) + t30 * t386; t388 * t144 + t387 * t176 + t35 + t38 + (t144 * t25 + t176 * t22 + t1 + t228) * m(7) + (t144 * t41 - t176 * t31 + t228 + t5) * m(6); -t143 * t84 - t176 * t96 + (-g(1) * t111 - g(2) * t109 - g(3) * t354 - t25 * t143 - t17 * t176 + t2) * m(7) + t37;];
tau  = t10;
