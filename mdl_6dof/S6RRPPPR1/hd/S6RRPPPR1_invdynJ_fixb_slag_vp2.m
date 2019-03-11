% Calculate vector of inverse dynamics joint torques for
% S6RRPPPR1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRPPPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:06:05
% EndTime: 2019-03-09 08:06:43
% DurationCPUTime: 23.59s
% Computational Cost: add. (8517->740), mult. (19524->942), div. (0->0), fcn. (14244->12), ass. (0->339)
t461 = -m(7) - m(6);
t332 = m(5) - t461;
t249 = sin(pkin(9));
t253 = sin(qJ(2));
t360 = cos(pkin(9));
t387 = cos(qJ(2));
t210 = t249 * t387 + t253 * t360;
t190 = t210 * qJD(1);
t248 = sin(pkin(10));
t250 = cos(pkin(10));
t162 = -t250 * qJD(2) + t190 * t248;
t401 = t162 / 0.2e1;
t302 = t360 * t387;
t338 = qJD(1) * t253;
t188 = -qJD(1) * t302 + t249 * t338;
t395 = t188 / 0.2e1;
t280 = qJD(2) * t248 + t250 * t190;
t399 = t280 / 0.2e1;
t474 = mrSges(6,2) + mrSges(5,3);
t473 = Ifges(5,1) + Ifges(6,1);
t472 = -Ifges(5,4) + Ifges(6,5);
t458 = Ifges(6,4) + Ifges(5,5);
t457 = Ifges(6,2) + Ifges(5,3);
t456 = Ifges(5,6) - Ifges(6,6);
t295 = t250 * mrSges(6,1) + t248 * mrSges(6,3);
t297 = mrSges(5,1) * t250 - mrSges(5,2) * t248;
t252 = sin(qJ(6));
t255 = cos(qJ(6));
t211 = t248 * t255 - t250 * t252;
t281 = t248 * t252 + t250 * t255;
t463 = t281 * mrSges(7,1) + mrSges(7,2) * t211;
t471 = t295 + t297 + t463;
t470 = m(7) * pkin(8) + mrSges(7,3) - t474;
t362 = qJDD(2) / 0.2e1;
t430 = 0.2e1 * t362;
t309 = t387 * qJDD(1);
t331 = qJD(1) * qJD(2);
t311 = t253 * t331;
t215 = t309 - t311;
t313 = qJD(2) * t387;
t216 = qJD(1) * t313 + t253 * qJDD(1);
t156 = t249 * t215 + t216 * t360;
t126 = -t250 * qJDD(2) + t156 * t248;
t406 = t126 / 0.2e1;
t127 = qJDD(2) * t248 + t156 * t250;
t469 = -t127 / 0.2e1;
t405 = t127 / 0.2e1;
t155 = -t360 * t215 + t216 * t249;
t468 = -t155 / 0.2e1;
t403 = t155 / 0.2e1;
t467 = t395 * t458 + t399 * t473 + t401 * t472;
t254 = sin(qJ(1));
t466 = g(2) * t254;
t222 = -mrSges(3,1) * t387 + t253 * mrSges(3,2);
t465 = -m(3) * pkin(1) - mrSges(2,1) + t222;
t112 = t211 * t188;
t193 = t211 * qJD(6);
t439 = t112 - t193;
t113 = t281 * t188;
t192 = t281 * qJD(6);
t438 = t113 - t192;
t346 = t248 * qJ(5);
t378 = t250 * pkin(4);
t464 = t346 + t378;
t462 = t162 * t252 + t255 * t280;
t94 = t162 * t255 - t252 * t280;
t29 = qJD(6) * t94 + t126 * t252 + t127 * t255;
t416 = t29 / 0.2e1;
t30 = -qJD(6) * t462 + t126 * t255 - t127 * t252;
t415 = t30 / 0.2e1;
t185 = qJD(6) - t188;
t384 = Ifges(7,4) * t462;
t33 = Ifges(7,2) * t94 + Ifges(7,6) * t185 + t384;
t414 = t33 / 0.2e1;
t93 = Ifges(7,4) * t94;
t34 = Ifges(7,1) * t462 + Ifges(7,5) * t185 + t93;
t413 = t34 / 0.2e1;
t152 = qJDD(6) - t155;
t404 = t152 / 0.2e1;
t460 = t215 / 0.2e1;
t459 = mrSges(6,3) - mrSges(5,2);
t383 = pkin(2) * t249;
t232 = qJ(4) + t383;
t376 = -pkin(8) + t232;
t200 = t376 * t248;
t201 = t376 * t250;
t137 = t200 * t255 - t201 * t252;
t328 = pkin(2) * t338;
t123 = pkin(3) * t190 + qJ(4) * t188 + t328;
t329 = t387 * pkin(7);
t223 = qJ(3) * t387 + t329;
t214 = t223 * qJD(1);
t194 = t249 * t214;
t251 = -qJ(3) - pkin(7);
t221 = t251 * t253;
t213 = qJD(1) * t221;
t154 = t213 * t360 - t194;
t143 = t248 * t154;
t408 = pkin(4) + pkin(5);
t35 = t143 + (pkin(8) * t188 - t123) * t250 - t408 * t190;
t357 = t188 * t248;
t71 = t248 * t123 + t250 * t154;
t55 = t190 * qJ(5) + t71;
t40 = -pkin(8) * t357 + t55;
t455 = qJD(4) * t281 + qJD(6) * t137 - t252 * t35 - t255 * t40;
t454 = t126 * t472 + t127 * t473 + t155 * t458;
t78 = -mrSges(5,2) * t155 - mrSges(5,3) * t126;
t79 = -mrSges(6,2) * t126 + mrSges(6,3) * t155;
t453 = t78 + t79;
t80 = mrSges(5,1) * t155 - mrSges(5,3) * t127;
t81 = -t155 * mrSges(6,1) + t127 * mrSges(6,2);
t452 = t81 - t80;
t451 = -t162 * t456 + t188 * t457 + t280 * t458;
t138 = t200 * t252 + t201 * t255;
t447 = qJD(4) * t211 - qJD(6) * t138 + t252 * t40 - t255 * t35;
t366 = t190 * Ifges(4,4);
t443 = Ifges(4,6) * qJD(2);
t446 = Ifges(7,5) * t462 - t188 * Ifges(4,2) + t94 * Ifges(7,6) + t185 * Ifges(7,3) + t366 + t443;
t308 = t360 * t214;
t153 = t213 * t249 + t308;
t334 = qJD(5) * t248;
t359 = qJ(5) * t250;
t425 = t248 * t408 - t359;
t445 = -t188 * t425 + t153 + t334;
t444 = Ifges(4,5) * qJD(2);
t440 = qJD(2) * mrSges(4,1) - mrSges(5,1) * t162 - mrSges(5,2) * t280 - mrSges(4,3) * t190;
t437 = -t248 * t456 + t250 * t458;
t371 = Ifges(6,5) * t248;
t373 = Ifges(5,4) * t248;
t436 = t250 * t473 + t371 - t373;
t239 = pkin(7) * t309;
t206 = -pkin(7) * t311 + t239;
t207 = t216 * pkin(7);
t435 = t387 * t206 + t207 * t253;
t358 = qJDD(1) * pkin(1);
t175 = -pkin(2) * t215 + qJDD(3) - t358;
t53 = pkin(3) * t155 - qJ(4) * t156 - qJD(4) * t190 + t175;
t312 = t387 * qJD(3);
t337 = qJD(2) * t253;
t326 = pkin(7) * t337;
t150 = t215 * qJ(3) + t239 + (t312 - t326) * qJD(1);
t333 = t253 * qJD(3);
t259 = qJDD(2) * pkin(2) - t216 * qJ(3) - qJD(1) * t333 - t207;
t77 = t360 * t150 + t249 * t259;
t69 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t77;
t19 = -t248 * t69 + t250 * t53;
t20 = t248 * t53 + t250 * t69;
t434 = -t19 * t248 + t20 * t250;
t256 = cos(qJ(1));
t433 = g(1) * t256 + t466;
t275 = mrSges(3,1) * t253 + mrSges(3,2) * t387;
t374 = Ifges(3,4) * t253;
t432 = pkin(1) * t275 - t253 * (Ifges(3,1) * t387 - t374) / 0.2e1;
t431 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t429 = Ifges(5,4) * t469 + Ifges(6,5) * t405 + Ifges(5,6) * t468 + Ifges(6,6) * t403 + (Ifges(5,2) + Ifges(6,3)) * t406;
t246 = qJ(2) + pkin(9);
t243 = sin(t246);
t244 = cos(t246);
t279 = -pkin(3) - t464;
t377 = t250 * pkin(5);
t382 = pkin(2) * t253;
t428 = t332 * t382 + t470 * t244 + (-m(7) * (t279 - t377) - m(6) * t279 + m(5) * pkin(3) + t471) * t243;
t298 = t244 * mrSges(4,1) - t243 * mrSges(4,2);
t427 = t243 * t474 + t298;
t424 = qJ(4) * t332;
t423 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1);
t14 = t155 * qJ(5) + t188 * qJD(5) + t20;
t11 = pkin(8) * t126 + t14;
t301 = qJDD(5) - t19;
t7 = -pkin(8) * t127 - t155 * t408 + t301;
t245 = t387 * pkin(2);
t238 = t245 + pkin(1);
t217 = -qJD(1) * t238 + qJD(3);
t104 = t188 * pkin(3) - t190 * qJ(4) + t217;
t204 = qJD(2) * pkin(2) + t213;
t147 = t249 * t204 + t308;
t136 = qJD(2) * qJ(4) + t147;
t59 = t104 * t250 - t248 * t136;
t299 = qJD(5) - t59;
t26 = -pkin(8) * t280 - t188 * t408 + t299;
t60 = t248 * t104 + t250 * t136;
t48 = t188 * qJ(5) + t60;
t36 = pkin(8) * t162 + t48;
t9 = -t252 * t36 + t255 * t26;
t1 = qJD(6) * t9 + t11 * t255 + t252 * t7;
t10 = t252 * t26 + t255 * t36;
t2 = -qJD(6) * t10 - t11 * t252 + t255 * t7;
t422 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t146 = t360 * t204 - t194;
t278 = qJD(2) * pkin(3) - qJD(4) + t146;
t294 = mrSges(6,1) * t248 - mrSges(6,3) * t250;
t296 = mrSges(5,1) * t248 + mrSges(5,2) * t250;
t388 = t248 / 0.2e1;
t262 = qJ(5) * t280 + t278;
t58 = pkin(4) * t162 - t262;
t421 = -t248 * (Ifges(5,4) * t280 - t162 * Ifges(5,2) + Ifges(5,6) * t188) / 0.2e1 + (Ifges(6,5) * t280 + Ifges(6,6) * t188 + t162 * Ifges(6,3)) * t388 + t217 * mrSges(4,2) + t58 * t294 - t146 * mrSges(4,3) + t444 / 0.2e1 - t278 * t296;
t47 = -pkin(4) * t188 + t299;
t420 = t147 * mrSges(4,3) + t47 * mrSges(6,1) + t60 * mrSges(5,2) + t9 * mrSges(7,1) + t443 / 0.2e1 - t10 * mrSges(7,2) - t217 * mrSges(4,1) - t48 * mrSges(6,3) - t59 * mrSges(5,1);
t418 = Ifges(7,4) * t416 + Ifges(7,2) * t415 + Ifges(7,6) * t404;
t417 = Ifges(7,1) * t416 + Ifges(7,4) * t415 + Ifges(7,5) * t404;
t412 = -t94 / 0.2e1;
t411 = t94 / 0.2e1;
t410 = -t462 / 0.2e1;
t409 = t462 / 0.2e1;
t407 = -t126 / 0.2e1;
t402 = -t162 / 0.2e1;
t400 = -t280 / 0.2e1;
t398 = -t185 / 0.2e1;
t397 = t185 / 0.2e1;
t396 = -t188 / 0.2e1;
t392 = -t190 / 0.2e1;
t391 = t190 / 0.2e1;
t386 = mrSges(7,3) * t94;
t385 = mrSges(7,3) * t462;
t379 = g(3) * t243;
t233 = t244 * pkin(3);
t372 = Ifges(5,4) * t250;
t370 = Ifges(6,5) * t250;
t101 = mrSges(6,1) * t162 - mrSges(6,3) * t280;
t41 = -mrSges(7,1) * t94 + mrSges(7,2) * t462;
t361 = t101 - t41;
t186 = qJD(2) * t221 + t312;
t187 = -qJD(2) * t223 - t333;
t125 = t186 * t360 + t249 * t187;
t189 = t210 * qJD(2);
t266 = -t249 * t253 + t302;
t191 = t266 * qJD(2);
t327 = pkin(2) * t337;
t98 = pkin(3) * t189 - qJ(4) * t191 - qJD(4) * t210 + t327;
t50 = t250 * t125 + t248 * t98;
t356 = t188 * t250;
t355 = t191 * t248;
t354 = t191 * t250;
t352 = t210 * t248;
t351 = t210 * t250;
t231 = t243 * qJ(4);
t348 = t243 * t256;
t347 = t244 * t256;
t345 = t250 * t256;
t344 = t251 * t256;
t343 = t254 * t248;
t342 = t254 * t250;
t341 = t256 * t248;
t105 = -mrSges(6,2) * t162 + mrSges(6,3) * t188;
t106 = -mrSges(5,2) * t188 - mrSges(5,3) * t162;
t340 = t105 + t106;
t107 = mrSges(5,1) * t188 - mrSges(5,3) * t280;
t108 = -mrSges(6,1) * t188 + mrSges(6,2) * t280;
t339 = t107 - t108;
t145 = -pkin(3) * t266 - qJ(4) * t210 - t238;
t161 = t249 * t221 + t223 * t360;
t90 = t248 * t145 + t250 * t161;
t336 = qJD(4) * t248;
t335 = qJD(4) * t250;
t330 = Ifges(7,5) * t29 + Ifges(7,6) * t30 + Ifges(7,3) * t152;
t323 = Ifges(3,4) * t387;
t67 = -qJ(5) * t266 + t90;
t320 = t233 + t231 + t245;
t318 = t360 * pkin(2);
t314 = qJD(1) * t387;
t310 = -t238 - t233;
t307 = t155 * mrSges(4,1) + t156 * mrSges(4,2);
t66 = t126 * mrSges(5,1) + t127 * mrSges(5,2);
t65 = t126 * mrSges(6,1) - t127 * mrSges(6,3);
t110 = t248 * t125;
t49 = t250 * t98 - t110;
t70 = t123 * t250 - t143;
t76 = -t150 * t249 + t360 * t259;
t157 = t248 * t161;
t89 = t145 * t250 - t157;
t124 = t186 * t249 - t360 * t187;
t160 = -t360 * t221 + t223 * t249;
t306 = t256 * t238 - t254 * t251;
t305 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t338) * t329;
t31 = t189 * qJ(5) - qJD(5) * t266 + t50;
t237 = -t318 - pkin(3);
t8 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t291 = -Ifges(5,2) * t248 + t372;
t288 = Ifges(6,3) * t248 + t370;
t287 = pkin(4) * t248 - t359;
t286 = t248 * t59 - t250 * t60;
t46 = t157 + (-pkin(8) * t210 - t145) * t250 + t408 * t266;
t52 = pkin(8) * t352 + t67;
t17 = -t252 * t52 + t255 * t46;
t18 = t252 * t46 + t255 * t52;
t74 = -mrSges(7,2) * t185 + t386;
t75 = mrSges(7,1) * t185 - t385;
t285 = -t252 * t75 + t255 * t74;
t284 = t244 * t464 + t320;
t176 = t244 * t343 + t345;
t177 = t244 * t342 - t341;
t283 = t176 * t255 - t177 * t252;
t282 = -t176 * t252 - t177 * t255;
t277 = pkin(3) * t347 + qJ(4) * t348 + t306;
t72 = -qJDD(2) * pkin(3) + qJDD(4) - t76;
t274 = Ifges(3,2) * t387 + t374;
t273 = Ifges(3,5) * t387 - Ifges(3,6) * t253;
t131 = t211 * t210;
t270 = t237 - t346;
t268 = -qJD(5) * t351 + t124;
t178 = t244 * t341 - t342;
t264 = -g(1) * t178 - g(2) * t176 - t248 * t379;
t21 = t126 * pkin(4) - t127 * qJ(5) - qJD(5) * t280 + t72;
t240 = Ifges(3,4) * t314;
t220 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t314;
t199 = Ifges(3,1) * t338 + Ifges(3,5) * qJD(2) + t240;
t198 = Ifges(3,6) * qJD(2) + qJD(1) * t274;
t197 = t270 - t378;
t184 = Ifges(4,4) * t188;
t179 = t244 * t345 + t343;
t169 = t250 * t408 - t270;
t166 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t188;
t142 = mrSges(4,1) * t188 + mrSges(4,2) * t190;
t140 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t156;
t139 = -qJDD(2) * mrSges(4,2) - mrSges(4,3) * t155;
t134 = t190 * Ifges(4,1) - t184 + t444;
t132 = t281 * t210;
t115 = t178 * t252 + t179 * t255;
t114 = t178 * t255 - t179 * t252;
t91 = t210 * t287 + t160;
t82 = -t188 * t287 + t153;
t73 = -t210 * t425 - t160;
t68 = pkin(4) * t266 - t89;
t64 = t191 * t211 - t192 * t210;
t63 = qJD(6) * t131 + t191 * t281;
t56 = -pkin(4) * t190 - t70;
t54 = t191 * t287 + t268;
t39 = -t162 * t408 + t262;
t38 = -t191 * t425 - t268;
t37 = -pkin(4) * t189 - t49;
t25 = pkin(8) * t355 + t31;
t24 = t110 + (-pkin(8) * t191 - t98) * t250 - t408 * t189;
t23 = -mrSges(7,2) * t152 + mrSges(7,3) * t30;
t22 = mrSges(7,1) * t152 - mrSges(7,3) * t29;
t16 = -pkin(4) * t155 + t301;
t15 = pkin(5) * t126 + t21;
t4 = -qJD(6) * t18 + t24 * t255 - t25 * t252;
t3 = qJD(6) * t17 + t24 * t252 + t25 * t255;
t5 = [t354 * t467 + (t1 * t131 + t10 * t64 - t132 * t2 - t63 * t9) * mrSges(7,3) + (Ifges(7,1) * t63 + Ifges(7,4) * t64) * t409 + (Ifges(7,1) * t132 + Ifges(7,4) * t131) * t416 + (-t305 + t273 * qJD(2) / 0.2e1) * qJD(2) - qJDD(2) * mrSges(3,2) * t329 + (-m(4) * t76 + m(5) * t72 - t140 + t66) * t160 + m(5) * (t19 * t89 + t20 * t90 + t49 * t59 + t50 * t60) + m(4) * (t125 * t147 + t161 * t77 - t175 * t238 + t217 * t327) + (t215 * t329 + t435) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t435) + t429 * t352 - t432 * t331 + (-t14 * t352 + t16 * t351 + t354 * t47 - t355 * t48) * mrSges(6,2) + (-m(4) * t146 - m(5) * t278 - t440) * t124 + Ifges(3,6) * t387 * t362 + t216 * t323 / 0.2e1 + (-t446 / 0.2e1 - Ifges(4,4) * t391 - Ifges(7,5) * t409 - Ifges(4,2) * t396 + Ifges(5,6) * t402 + Ifges(6,6) * t401 - Ifges(7,6) * t411 - Ifges(7,3) * t397 + t457 * t395 + t458 * t399 - t420 + t451 / 0.2e1) * t189 + t274 * t460 + (Ifges(3,1) * t216 + Ifges(3,4) * t460 - pkin(7) * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t216) + t430 * Ifges(3,5)) * t253 + t454 * t351 / 0.2e1 + (Ifges(4,1) * t391 + Ifges(4,4) * t396 + t288 * t401 + t291 * t402 + t134 / 0.2e1 + t436 * t399 + t437 * t395 + t421) * t191 + (t175 * mrSges(4,2) - t76 * mrSges(4,3) + Ifges(4,1) * t156 - Ifges(4,4) * t155 + Ifges(4,5) * t430 + t21 * t294 + t288 * t406 + t291 * t407 + t72 * t296 + t403 * t437 + t405 * t436) * t210 + (qJD(1) * (-Ifges(3,2) * t253 + t323) + t199) * t313 / 0.2e1 - t198 * t337 / 0.2e1 + (Ifges(7,5) * t63 + Ifges(7,6) * t64) * t397 + (Ifges(7,5) * t132 + Ifges(7,6) * t131) * t404 + m(6) * (t14 * t67 + t16 * t68 + t21 * t91 + t31 * t48 + t37 * t47 + t54 * t58) + m(7) * (t1 * t18 + t10 * t3 - t15 * t73 + t17 * t2 + t38 * t39 + t4 * t9) + (-m(4) * t306 - m(5) * t277 - t115 * mrSges(7,1) - t114 * mrSges(7,2) + t461 * (t179 * pkin(4) + t178 * qJ(5) + t277) - t423 * t179 - t459 * t178 + t470 * t348 + (-t298 + t465) * t256 + t431 * t254) * g(2) - t220 * t326 + t387 * (Ifges(3,4) * t216 + Ifges(3,2) * t215 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (-t19 * t351 - t20 * t352 - t354 * t59 - t355 * t60) * mrSges(5,3) + (t330 / 0.2e1 - Ifges(4,2) * t155 + Ifges(4,4) * t156 + t77 * mrSges(4,3) - t19 * mrSges(5,1) + Ifges(7,3) * t404 - Ifges(5,6) * t407 + Ifges(7,6) * t415 + Ifges(7,5) * t416 + t16 * mrSges(6,1) + t20 * mrSges(5,2) - t14 * mrSges(6,3) - t175 * mrSges(4,1) + t422 + (-t405 + t469) * t458 + (-t403 + t468) * t457 + Ifges(4,6) * t430 + (-Ifges(6,6) + t456) * t406) * t266 + Ifges(2,3) * qJDD(1) - t238 * t307 + (Ifges(7,4) * t63 + Ifges(7,2) * t64) * t411 + (Ifges(7,4) * t132 + Ifges(7,2) * t131) * t415 + t38 * t41 + t39 * (-mrSges(7,1) * t64 + mrSges(7,2) * t63) + t142 * t327 + t73 * t8 + t3 * t74 + t4 * t75 + t67 * t79 + t68 * t81 + t89 * t80 + t90 * t78 + (m(5) * t344 - t282 * mrSges(7,1) + t283 * mrSges(7,2) + t461 * (-t177 * pkin(4) - t176 * qJ(5) - t344) + t423 * t177 + t459 * t176 + (m(4) * t251 + t431) * t256 + (-m(7) * t310 - (m(7) * (pkin(8) - qJ(4)) + mrSges(7,3)) * t243 + m(4) * t238 + (-m(6) - m(5)) * (t310 - t231) + t427 - t465) * t254) * g(1) + t91 * t65 + t54 * t101 + t31 * t105 + t50 * t106 + t49 * t107 + t37 * t108 + t63 * t413 + t64 * t414 + t132 * t417 + t131 * t418 - t15 * (-mrSges(7,1) * t131 + mrSges(7,2) * t132) + t161 * t139 + t125 * t166 - t222 * t358 + t17 * t22 + t18 * t23 - pkin(1) * (-mrSges(3,1) * t215 + mrSges(3,2) * t216); (-t244 * t424 + t428) * t466 + t356 * t467 + (-t334 - t82) * t101 + t445 * t41 + t446 * t391 + t447 * t75 + (t335 - t71) * t106 + (m(4) * t382 + mrSges(4,1) * t243 + mrSges(4,2) * t244 + t275) * t433 + (t146 * t153 - t147 * t154 - t217 * t328 + (t249 * t77 + t360 * t76) * pkin(2)) * m(4) + (-t356 * t59 - t357 * t60 + t434) * mrSges(5,3) + (t197 * t21 - t47 * t56 - t48 * t55 - t58 * t82) * m(6) + (qJD(1) * t432 + t305) * qJD(1) + (-m(7) * (-pkin(8) * t243 + t284) + t243 * mrSges(7,3) - m(6) * t284 - m(5) * t320 - m(4) * t245 + t222 + (-m(7) * t377 - t471) * t244 - t427) * g(3) + (t336 - t56) * t108 + (-t286 * qJD(4) + t153 * t278 + t232 * t434 + t237 * t72 - t59 * t70 - t60 * t71) * m(5) + (Ifges(7,5) * t211 - Ifges(7,6) * t281) * t404 + (Ifges(7,4) * t211 - Ifges(7,2) * t281) * t415 + (Ifges(7,1) * t211 - Ifges(7,4) * t281) * t416 + (-t1 * t281 + t10 * t439 - t2 * t211 - t438 * t9) * mrSges(7,3) - t281 * t418 + (t256 * t428 - t347 * t424) * g(1) + t140 * t318 + (-t336 - t70) * t107 + (t453 * t232 + Ifges(5,2) * t407 - Ifges(6,3) * t406 + t403 * t456 - t429 + (qJD(4) * t48 + t14 * t232) * m(6) + t14 * mrSges(6,2)) * t250 + (t1 * t138 + t10 * t455 + t137 * t2 - t15 * t169 + t39 * t445 + t447 * t9) * m(7) + t455 * t74 + (-Ifges(7,5) * t410 - Ifges(4,2) * t395 + Ifges(5,6) * t401 + Ifges(6,6) * t402 - Ifges(7,6) * t412 - Ifges(7,3) * t398 + t396 * t457 + t400 * t458 + t420) * t190 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-t366 + t451) * t392 + t454 * t388 + t371 * t406 + t373 * t407 + (pkin(7) * t220 + t198 / 0.2e1) * t338 + (t356 * t47 - t357 * t48) * mrSges(6,2) + (-t184 + t134) * t395 + (-Ifges(7,5) * t192 - Ifges(7,6) * t193) * t397 + (-Ifges(7,1) * t192 - Ifges(7,4) * t193) * t409 + (-Ifges(7,4) * t192 - Ifges(7,2) * t193) * t411 + (-Ifges(7,1) * t113 - Ifges(7,4) * t112) * t410 + (-Ifges(7,5) * t113 - Ifges(7,6) * t112) * t398 + (-Ifges(7,4) * t113 - Ifges(7,2) * t112) * t412 + (-Ifges(4,1) * t392 - t288 * t402 - t291 * t401 - t396 * t437 - t400 * t436 + t421) * t188 + (-mrSges(7,1) * t439 + mrSges(7,2) * t438) * t39 + t440 * t153 - (-Ifges(3,2) * t338 + t199 + t240) * t314 / 0.2e1 - t273 * t331 / 0.2e1 + (t335 - t55) * t105 - t142 * t328 - t21 * t295 - t72 * t297 + t237 * t66 + (-t370 + t372) * t405 - t15 * t463 + t438 * t413 + t439 * t414 + (t452 * t232 + t473 * t405 + t458 * t403 + (qJD(4) * t47 - qJD(5) * t58 + t16 * t232) * m(6) + t16 * mrSges(6,2)) * t248 + t76 * mrSges(4,1) - t77 * mrSges(4,2) + t139 * t383 + t211 * t417 + t137 * t22 + t138 * t23 - Ifges(4,6) * t155 + Ifges(4,5) * t156 - t154 * t166 + t169 * t8 + t197 * t65 - t206 * mrSges(3,2) - t207 * mrSges(3,1) + Ifges(3,6) * t215 + Ifges(3,5) * t216; -t281 * t22 + t211 * t23 + t439 * t75 + t438 * t74 - t452 * t250 + t453 * t248 + (-t361 + t440) * t190 + (-t248 * t339 + t250 * t340 + t166) * t188 + t307 + (-g(1) * t254 + g(2) * t256) * (m(4) + t332) + (t1 * t211 + t10 * t438 + t190 * t39 - t2 * t281 + t439 * t9) * m(7) + (t14 * t248 - t16 * t250 - t190 * t58 - (-t248 * t47 - t250 * t48) * t188) * m(6) + (-t188 * t286 + t19 * t250 + t190 * t278 + t20 * t248) * m(5) + (t146 * t190 + t147 * t188 + t175) * m(4); t339 * t280 + t340 * t162 + t94 * t74 - t462 * t75 + t65 + t66 - t8 + (t10 * t94 - t462 * t9 + t15) * m(7) + (t162 * t48 - t280 * t47 + t21) * m(6) + (t162 * t60 + t280 * t59 + t72) * m(5) + (t244 * g(3) - t243 * t433) * t332; t255 * t22 + t252 * t23 + t361 * t280 + t285 * qJD(6) + (-t105 - t285) * t188 + t81 + (t1 * t252 - t280 * t39 + t2 * t255 + t264 + t185 * (t10 * t255 - t252 * t9)) * m(7) + (-t188 * t48 + t280 * t58 + t16 + t264) * m(6); -t39 * (mrSges(7,1) * t462 + mrSges(7,2) * t94) + (Ifges(7,1) * t94 - t384) * t410 + t33 * t409 + (Ifges(7,5) * t94 - Ifges(7,6) * t462) * t398 - g(1) * (mrSges(7,1) * t114 - mrSges(7,2) * t115) - g(2) * (mrSges(7,1) * t283 + mrSges(7,2) * t282) - (mrSges(7,1) * t211 - mrSges(7,2) * t281) * t379 + t330 + (-t74 + t386) * t9 + (-Ifges(7,2) * t462 + t34 + t93) * t412 + (t385 + t75) * t10 + t422;];
tau  = t5;
