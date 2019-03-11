% Calculate vector of inverse dynamics joint torques for
% S6PRRPRP4
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:17
% EndTime: 2019-03-08 21:40:52
% DurationCPUTime: 22.26s
% Computational Cost: add. (4483->644), mult. (10011->843), div. (0->0), fcn. (6807->10), ass. (0->291)
t479 = Ifges(6,4) + Ifges(7,4);
t484 = -mrSges(4,1) + mrSges(5,2);
t480 = Ifges(6,1) + Ifges(7,1);
t457 = Ifges(6,5) + Ifges(7,5);
t478 = Ifges(6,2) + Ifges(7,2);
t456 = Ifges(6,6) + Ifges(7,6);
t216 = cos(qJ(5));
t483 = t479 * t216;
t213 = sin(qJ(5));
t482 = t479 * t213;
t217 = cos(qJ(3));
t337 = qJD(2) * t217;
t154 = -qJD(3) * t213 - t216 * t337;
t214 = sin(qJ(3));
t326 = qJD(2) * qJD(3);
t164 = qJDD(2) * t214 + t217 * t326;
t215 = sin(qJ(2));
t210 = sin(pkin(6));
t342 = qJD(1) * t210;
t311 = t215 * t342;
t165 = qJD(2) * pkin(8) + t311;
t211 = cos(pkin(6));
t324 = qJDD(1) * t211;
t333 = qJD(3) * t217;
t218 = cos(qJ(2));
t340 = qJD(2) * t210;
t304 = qJD(1) * t340;
t184 = t218 * t304;
t325 = qJDD(1) * t210;
t125 = t215 * t325 + t184;
t336 = qJD(3) * t211;
t471 = qJDD(2) * pkin(8) + qJD(1) * t336 + t125;
t34 = -t165 * t333 - t214 * t471 + t217 * t324;
t231 = qJDD(4) - t34;
t402 = pkin(3) + pkin(9);
t17 = pkin(4) * t164 - qJDD(3) * t402 + t231;
t163 = -t217 * qJDD(2) + t214 * t326;
t183 = t215 * t304;
t124 = t218 * t325 - t183;
t111 = -qJDD(2) * pkin(2) - t124;
t332 = qJD(4) * t214;
t221 = -qJ(4) * t164 - qJD(2) * t332 + t111;
t29 = t163 * t402 + t221;
t330 = qJD(5) * t216;
t331 = qJD(5) * t213;
t341 = qJD(1) * t211;
t191 = t217 * t341;
t87 = -t214 * (pkin(4) * qJD(2) + t165) + t191;
t445 = qJD(4) - t87;
t67 = -qJD(3) * t402 + t445;
t358 = qJ(4) * t214;
t293 = -pkin(2) - t358;
t149 = -t217 * t402 + t293;
t310 = t218 * t342;
t89 = qJD(2) * t149 - t310;
t3 = t213 * t17 + t216 * t29 + t67 * t330 - t331 * t89;
t334 = qJD(3) * t216;
t236 = t213 * t337 - t334;
t69 = qJD(5) * t236 - qJDD(3) * t213 + t163 * t216;
t2 = qJ(6) * t69 + qJD(6) * t154 + t3;
t481 = t2 * mrSges(7,2);
t459 = Ifges(4,5) - Ifges(5,4);
t458 = Ifges(5,5) - Ifges(4,6);
t477 = Ifges(6,3) + Ifges(7,3);
t339 = qJD(2) * t214;
t314 = mrSges(4,3) * t339;
t316 = mrSges(5,1) * t339;
t440 = qJD(3) * t484 + t314 + t316;
t433 = t213 * t457 + t216 * t456;
t431 = t216 * t478 + t482;
t428 = t213 * t480 + t483;
t153 = qJDD(5) + t164;
t26 = t213 * t67 + t216 * t89;
t4 = -qJD(5) * t26 + t216 * t17 - t213 * t29;
t68 = qJD(5) * t154 + qJDD(3) * t216 + t163 * t213;
t1 = pkin(5) * t153 - qJ(6) * t68 + qJD(6) * t236 + t4;
t355 = t210 * t215;
t312 = t214 * t355;
t138 = -t211 * t217 + t312;
t359 = sin(pkin(10));
t285 = t359 * t218;
t360 = cos(pkin(10));
t288 = t360 * t215;
t133 = t211 * t288 + t285;
t290 = t210 * t360;
t79 = t133 * t214 + t217 * t290;
t286 = t359 * t215;
t287 = t360 * t218;
t135 = -t211 * t286 + t287;
t289 = t210 * t359;
t81 = t135 * t214 - t217 * t289;
t407 = g(1) * t81 + g(2) * t79 + g(3) * t138;
t476 = t1 * t216 + t2 * t213 - t407;
t370 = Ifges(5,6) * t214;
t247 = -t217 * Ifges(5,3) - t370;
t195 = qJD(5) + t339;
t470 = t479 * t154;
t451 = t457 * t195 - t480 * t236 + t470;
t468 = t479 * t236;
t452 = t478 * t154 + t456 * t195 - t468;
t475 = Ifges(5,5) * qJD(3) + qJD(2) * t247 + t213 * t451 + t216 * t452;
t474 = -t478 * t69 / 0.2e1 - t479 * t68 / 0.2e1 - t456 * t153 / 0.2e1;
t465 = -m(7) - m(5);
t420 = -m(6) + t465;
t473 = -m(4) + t420;
t246 = pkin(9) * t214 - qJ(4) * t217;
t335 = qJD(3) * t214;
t280 = pkin(3) * t335 - t332;
t107 = qJD(3) * t246 + t280;
t350 = t214 * t218;
t109 = (t213 * t350 + t215 * t216) * t210;
t401 = pkin(4) + pkin(8);
t162 = t401 * t333;
t181 = t401 * t214;
t454 = -qJD(1) * t109 + t216 * t107 - t149 * t331 + t213 * t162 + t181 * t330;
t348 = t216 * t218;
t108 = (-t213 * t215 + t214 * t348) * t210;
t472 = -qJD(1) * t108 + t216 * t162;
t212 = -qJ(6) - pkin(9);
t264 = mrSges(5,2) * t217 - mrSges(5,3) * t214;
t269 = mrSges(4,1) * t217 - mrSges(4,2) * t214;
t353 = t213 * t214;
t469 = -m(7) * (pkin(5) * t353 - t212 * t217) - t217 * mrSges(7,3) - mrSges(3,1) - t269 + t264;
t466 = m(6) * pkin(9);
t467 = -m(7) * t212 - pkin(3) * t420 + t466 - t484;
t281 = qJ(6) * t217 - t149;
t327 = qJD(6) * t217;
t463 = pkin(5) * t333 + t281 * t330 + (-qJ(6) * t335 - qJD(5) * t181 - t107 + t327) * t213 + t472;
t329 = qJD(5) * t217;
t307 = t213 * t329;
t233 = t214 * t334 + t307;
t462 = qJ(6) * t233 - t216 * t327 + t454;
t461 = -mrSges(6,1) - mrSges(7,1);
t460 = -mrSges(6,2) - mrSges(7,2);
t455 = t153 * t457 + t479 * t69 + t480 * t68;
t76 = t216 * t149 + t213 * t181;
t453 = -qJD(5) * t76 - t107 * t213 + t472;
t450 = m(7) * pkin(5) + mrSges(7,1);
t346 = qJ(6) + t402;
t200 = pkin(3) * t339;
t120 = qJD(2) * t246 + t200;
t105 = t217 * t165 + t214 * t341;
t88 = pkin(4) * t337 + t105;
t36 = -t120 * t213 + t216 * t88;
t449 = -(pkin(5) * t217 - qJ(6) * t353) * qJD(2) - t36 - qJD(6) * t216 + t331 * t346;
t168 = t346 * t216;
t37 = t216 * t120 + t213 * t88;
t448 = -qJ(6) * t216 * t339 - qJD(5) * t168 - qJD(6) * t213 - t37;
t198 = pkin(5) * t216 + pkin(4);
t447 = -t191 - (-qJD(2) * t198 - t165) * t214 + pkin(5) * t330 + qJD(4);
t315 = mrSges(5,1) * t337;
t175 = -qJD(3) * mrSges(5,3) - t315;
t78 = -mrSges(6,1) * t154 - mrSges(6,2) * t236;
t446 = t78 - t175;
t132 = -t211 * t287 + t286;
t357 = t132 * t217;
t444 = -pkin(3) * t357 - t132 * t358;
t134 = t211 * t285 + t288;
t356 = t134 * t217;
t443 = -pkin(3) * t356 - t134 * t358;
t128 = mrSges(5,1) * t163 - qJDD(3) * mrSges(5,3);
t442 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t163 - t128;
t129 = t164 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t441 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t164 + t129;
t313 = mrSges(4,3) * t337;
t174 = -qJD(3) * mrSges(4,2) + t313;
t439 = t174 - t175;
t438 = t214 * t433 + t217 * t477;
t437 = t214 * t431 + t217 * t456;
t436 = t214 * t428 + t217 * t457;
t369 = Ifges(5,6) * t217;
t435 = t214 * (-Ifges(5,2) * t217 + t370) + t217 * (Ifges(5,3) * t214 - t369);
t169 = -pkin(3) * t217 + t293;
t106 = qJD(2) * t169 - t310;
t166 = -qJD(2) * pkin(2) - t310;
t434 = t166 * (mrSges(4,1) * t214 + mrSges(4,2) * t217) + t106 * (-mrSges(5,2) * t214 - mrSges(5,3) * t217);
t432 = -t213 * t456 + t216 * t457;
t430 = -t213 * t478 + t483;
t429 = t214 * t458 + t217 * t459;
t427 = t216 * t480 - t482;
t425 = t153 * t477 + t456 * t69 + t457 * t68;
t33 = -t165 * t335 + t214 * t324 + t217 * t471;
t424 = -t214 * t34 + t217 * t33;
t30 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t33;
t31 = -qJDD(3) * pkin(3) + t231;
t423 = t214 * t31 - t217 * t30;
t422 = t213 * t3 + t216 * t4;
t419 = -mrSges(6,1) - t450;
t199 = Ifges(4,4) * t337;
t418 = Ifges(4,1) * t339 + Ifges(4,5) * qJD(3) + t154 * t456 + t195 * t477 - t236 * t457 + t199;
t209 = qJD(3) * qJ(4);
t70 = t209 + t88;
t46 = -pkin(5) * t154 + qJD(6) + t70;
t77 = -mrSges(7,1) * t154 - mrSges(7,2) * t236;
t416 = -m(7) * t46 - t77;
t266 = mrSges(7,1) * t216 - mrSges(7,2) * t213;
t268 = mrSges(6,1) * t216 - mrSges(6,2) * t213;
t411 = t46 * t266 + t70 * t268;
t410 = -m(7) * t198 - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t197 = pkin(5) * t213 + qJ(4);
t265 = t213 * mrSges(7,1) + t216 * mrSges(7,2);
t267 = mrSges(6,1) * t213 + mrSges(6,2) * t216;
t408 = -m(7) * t197 + mrSges(4,2) - mrSges(5,3) - t265 - t267 + (-m(5) - m(6)) * qJ(4);
t406 = -m(6) * t401 + t410;
t220 = qJD(2) ^ 2;
t404 = t68 / 0.2e1;
t403 = t69 / 0.2e1;
t400 = t153 / 0.2e1;
t398 = t154 / 0.2e1;
t396 = -t236 / 0.2e1;
t394 = t195 / 0.2e1;
t204 = t217 * pkin(8);
t42 = mrSges(7,1) * t153 - mrSges(7,3) * t68;
t43 = mrSges(6,1) * t153 - mrSges(6,3) * t68;
t382 = t42 + t43;
t44 = -mrSges(7,2) * t153 + mrSges(7,3) * t69;
t45 = -mrSges(6,2) * t153 + mrSges(6,3) * t69;
t381 = -t44 - t45;
t380 = mrSges(6,3) * t154;
t379 = mrSges(6,3) * t236;
t378 = mrSges(7,3) * t154;
t377 = mrSges(7,3) * t236;
t376 = Ifges(4,4) * t214;
t375 = Ifges(4,4) * t217;
t354 = t210 * t218;
t352 = t213 * t217;
t351 = t214 * t216;
t349 = t216 * t217;
t347 = t217 * t218;
t100 = -mrSges(7,2) * t195 + t378;
t101 = -mrSges(6,2) * t195 + t380;
t345 = t100 + t101;
t102 = mrSges(7,1) * t195 + t377;
t103 = mrSges(6,1) * t195 + t379;
t344 = t102 + t103;
t182 = t217 * pkin(4) + t204;
t338 = qJD(2) * t215;
t328 = qJD(5) * t402;
t321 = -t77 - t446;
t309 = t210 * t338;
t308 = t218 * t340;
t19 = -t69 * mrSges(7,1) + t68 * mrSges(7,2);
t294 = -t330 / 0.2e1;
t122 = t132 * pkin(2);
t292 = pkin(8) * t133 - t122;
t123 = t134 * pkin(2);
t291 = pkin(8) * t135 - t123;
t25 = -t213 * t89 + t216 * t67;
t284 = -t326 / 0.2e1;
t104 = t165 * t214 - t191;
t259 = t217 * Ifges(4,2) + t376;
t248 = -Ifges(5,2) * t214 - t369;
t245 = t25 * t213 - t26 * t216;
t13 = qJ(6) * t236 + t25;
t86 = -t138 * t213 + t210 * t348;
t85 = t138 * t216 + t213 * t354;
t139 = t211 * t214 + t217 * t355;
t239 = t214 * (Ifges(4,1) * t217 - t376);
t232 = t213 * t335 - t216 * t329;
t18 = -pkin(4) * t163 - t30;
t222 = -qJD(5) * t245 + t422;
t167 = t346 * t213;
t161 = t401 * t335;
t160 = -qJ(4) * t337 + t200;
t159 = t269 * qJD(2);
t158 = t264 * qJD(2);
t157 = t216 * t181;
t143 = Ifges(5,4) * qJD(3) + qJD(2) * t248;
t140 = Ifges(4,6) * qJD(3) + qJD(2) * t259;
t137 = pkin(5) * t349 + t182;
t131 = -qJ(4) * t333 + t280;
t94 = -t209 - t105;
t93 = -qJD(3) * pkin(3) + qJD(4) + t104;
t92 = -mrSges(5,2) * t163 - mrSges(5,3) * t164;
t91 = mrSges(4,1) * t163 + mrSges(4,2) * t164;
t90 = -pkin(5) * t307 + (-pkin(8) - t198) * t335;
t84 = -qJD(3) * t312 + (t308 + t336) * t217;
t83 = qJD(3) * t139 + t214 * t308;
t82 = t135 * t217 + t214 * t289;
t80 = t133 * t217 - t214 * t290;
t75 = -t149 * t213 + t157;
t60 = -qJ(6) * t349 + t76;
t51 = pkin(5) * t214 + t213 * t281 + t157;
t35 = pkin(3) * t163 + t221;
t24 = qJD(5) * t85 + t213 * t83 + t216 * t309;
t23 = qJD(5) * t86 - t213 * t309 + t216 * t83;
t20 = -mrSges(6,1) * t69 + mrSges(6,2) * t68;
t14 = qJ(6) * t154 + t26;
t8 = pkin(5) * t195 + t13;
t5 = -pkin(5) * t69 + qJDD(6) + t18;
t6 = [m(2) * qJDD(1) + t381 * t86 + t382 * t85 + t440 * t83 + t345 * t24 + t344 * t23 + t441 * t138 + (t174 - t321) * t84 + (t19 + t20 + t442) * t139 + (-m(2) - m(3) + t473) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t220 - t91 - t92) * t218 + (-mrSges(3,1) * t220 - mrSges(3,2) * qJDD(2) + (t158 - t159) * qJD(2)) * t215) * t210 + m(4) * (t104 * t83 + t105 * t84 - t138 * t34 + t139 * t33 + (-t111 * t218 + t166 * t338) * t210) + m(5) * (t138 * t31 - t139 * t30 + t83 * t93 - t84 * t94 + (t106 * t338 - t218 * t35) * t210) + m(3) * (qJDD(1) * t211 ^ 2 + (t124 * t218 + t125 * t215) * t210) + m(7) * (t1 * t85 + t139 * t5 + t14 * t24 - t2 * t86 + t23 * t8 + t46 * t84) + m(6) * (t139 * t18 + t23 * t25 + t24 * t26 - t3 * t86 + t4 * t85 + t70 * t84); (-t14 * mrSges(7,2) - t26 * mrSges(6,2) + t8 * mrSges(7,1) + t25 * mrSges(6,1) - t143 / 0.2e1 + t93 * mrSges(5,1) + t104 * mrSges(4,3) + t418 / 0.2e1) * t333 + t3 * (-mrSges(6,2) * t214 - mrSges(6,3) * t349) + t440 * (pkin(8) * t333 - t214 * t310) + t164 * t375 / 0.2e1 + t46 * (-mrSges(7,1) * t233 + mrSges(7,2) * t232) + t70 * (-mrSges(6,1) * t233 + mrSges(6,2) * t232) + (t217 * (-Ifges(4,2) * t214 + t375) + t239) * t326 / 0.2e1 + t164 * Ifges(4,1) * t214 + t1 * mrSges(7,1) * t214 + (t124 + t183) * mrSges(3,1) + t5 * t266 * t217 + t18 * t268 * t217 + t159 * t311 + t423 * mrSges(5,1) + (t1 * t352 + t14 * t233 - t2 * t349 - t232 * t8) * mrSges(7,3) + t51 * t42 + t60 * t44 + t4 * (mrSges(6,1) * t214 + mrSges(6,3) * t352) + (-t125 + t184) * mrSges(3,2) + t349 * t474 + t424 * mrSges(4,3) + (-t140 / 0.2e1 + t94 * mrSges(5,1) - t105 * mrSges(4,3) - t439 * pkin(8) + t475 / 0.2e1) * t335 + (-m(6) * (-pkin(9) * t356 - t123 + t443) + mrSges(6,3) * t356 - m(4) * t291 + t461 * (-t134 * t353 + t135 * t216) + t460 * (-t134 * t351 - t135 * t213) + t465 * (t291 + t443) + t406 * t135 - t469 * t134) * g(1) + (-m(6) * (-pkin(9) * t357 - t122 + t444) + mrSges(6,3) * t357 - m(4) * t292 + t461 * (-t132 * t353 + t133 * t216) + t460 * (-t132 * t351 - t133 * t213) + t465 * (t292 + t444) + t406 * t133 - t469 * t132) * g(2) + (t460 * t108 + t461 * t109 + t473 * (pkin(2) * t354 + pkin(8) * t355) + (t420 * (pkin(3) * t347 + qJ(4) * t350) + (-mrSges(6,3) - t466) * t347 + t469 * t218 + (-m(6) * pkin(4) + t410) * t215) * t210) * g(3) + (-t311 + t131) * t158 - t25 * mrSges(6,3) * t232 + t26 * mrSges(6,3) * t233 + Ifges(3,3) * qJDD(2) + (t1 * t51 + t137 * t5 + t14 * t462 + t2 * t60 + t46 * t90 + t463 * t8) * m(7) + t462 * t100 + t463 * t102 + (t214 * t457 - t217 * t428) * t404 + (t214 * t459 - t217 * t458) * qJDD(3) / 0.2e1 + t452 * t307 / 0.2e1 + t453 * t103 + t454 * t101 + (-t161 * t70 + t18 * t182 + t453 * t25 + t454 * t26 + t3 * t76 + t4 * t75) * m(6) - t455 * t352 / 0.2e1 + (t214 * t456 - t217 * t431) * t403 + t442 * t204 + t435 * t284 + (t394 * t438 + t396 * t436 + t398 * t437 + t434) * qJD(3) + (-(t106 * t215 + (t214 * t93 - t217 * t94) * t218) * t342 + t106 * t131 + t169 * t35 + ((t214 * t94 + t217 * t93) * qJD(3) + t423) * pkin(8)) * m(5) + (-(t166 * t215 + (t104 * t214 + t105 * t217) * t218) * t342 - pkin(2) * t111 + ((t104 * t217 - t105 * t214) * qJD(3) + t424) * pkin(8)) * m(4) + t429 * qJD(3) ^ 2 / 0.2e1 + (-t394 * t432 - t396 * t427 - t398 * t430) * t329 + (-Ifges(4,4) * t163 + Ifges(4,5) * qJDD(3) + t425) * t214 / 0.2e1 + t217 * (Ifges(4,4) * t164 - Ifges(4,2) * t163 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t217 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t164 + Ifges(5,3) * t163) / 0.2e1 + t75 * t43 + t76 * t45 + t90 * t77 - pkin(2) * t91 + t163 * t247 / 0.2e1 - t164 * t248 / 0.2e1 + (t214 * t477 - t433 * t217) * t400 - t163 * t259 / 0.2e1 + t35 * t264 - t111 * t269 - t214 * t481 + t137 * t19 - t161 * t78 + t169 * t92 + t182 * t20 + (-m(6) * t70 + t416 - t439 - t78) * t217 * t310 + t441 * pkin(8) * t214 + t451 * t217 * t294 - t214 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t164 + Ifges(5,6) * t163) / 0.2e1; (-t14 * t330 + t331 * t8 - t476) * mrSges(7,3) + (-t239 / 0.2e1 + t435 / 0.2e1) * t220 + t143 * t337 / 0.2e1 + (t138 * t467 + t139 * t408) * g(3) + (t408 * t80 + t467 * t79) * g(2) + (t408 * t82 + t467 * t81) * g(1) - t93 * t315 - t94 * t316 + (t20 - t128) * qJ(4) + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + (-pkin(3) * t31 - qJ(4) * t30 - qJD(4) * t94 - t106 * t160) * m(5) - (t154 * t437 + t195 * t438 - t236 * t436) * qJD(2) / 0.2e1 - (t431 * t154 + t433 * t195 - t236 * t428) * qJD(5) / 0.2e1 - (t213 * t45 + t216 * t43) * t402 + (qJ(4) * t18 - t222 * t402 - t25 * t36 - t26 * t37 + t445 * t70) * m(6) - t30 * mrSges(5,3) + t31 * mrSges(5,2) - t33 * mrSges(4,2) + t34 * mrSges(4,1) + (t213 * t328 - t36) * t103 + (-t216 * t328 - t37) * t101 + t213 * t474 - t475 * t339 / 0.2e1 + t458 * t163 + t459 * t164 - t451 * t331 / 0.2e1 + t452 * t294 + t455 * t216 / 0.2e1 + t446 * qJD(4) + t447 * t77 + t448 * t100 + t449 * t102 + (-t1 * t168 + t14 * t448 - t167 * t2 + t197 * t5 + t447 * t46 + t449 * t8) * m(7) + (-m(5) * t94 - t313 + t439) * t104 + (-m(5) * t93 + t314 - t440) * t105 + (-t14 * (-mrSges(7,2) * t217 + mrSges(7,3) * t351) - t26 * (-mrSges(6,2) * t217 + mrSges(6,3) * t351) - t8 * (mrSges(7,1) * t217 - mrSges(7,3) * t353) - t25 * (mrSges(6,1) * t217 - mrSges(6,3) * t353) - t434) * qJD(2) + t427 * t404 + t429 * t284 + t430 * t403 + t432 * t400 + (t25 * t331 - t26 * t330 + t407 - t422) * mrSges(6,3) - (-Ifges(4,2) * t339 + t199 + t418) * t337 / 0.2e1 + (t140 / 0.2e1 + t411) * t339 + t411 * qJD(5) - t87 * t78 + t5 * t265 + t18 * t267 - pkin(3) * t129 - t160 * t158 - t167 * t44 - t168 * t42 + t197 * t19; t158 * t339 + t321 * qJD(3) + (t195 * t345 + t382) * t216 + (-t195 * t344 - t381) * t213 + t129 + (-qJD(3) * t46 + t195 * (t14 * t216 - t8 * t213) + t476) * m(7) + (-qJD(3) * t70 - t245 * t339 + t222 - t407) * m(6) + (qJD(3) * t94 + t106 * t339 + t31 - t407) * m(5); t450 * t1 + t452 * t396 + (t380 - t101) * t25 + t4 * mrSges(6,1) - t481 - t3 * mrSges(6,2) + t8 * t378 + (t480 * t154 + t468) * t236 / 0.2e1 + (-t377 - m(7) * (t13 - t8) + t102) * t14 + (-t379 + t103) * t26 - t46 * (-mrSges(7,1) * t236 + mrSges(7,2) * t154) - t70 * (-mrSges(6,1) * t236 + mrSges(6,2) * t154) - (t154 * t457 + t236 * t456) * t195 / 0.2e1 - (t478 * t236 + t451 + t470) * t154 / 0.2e1 + t425 + (t419 * t85 + t460 * t86) * g(3) + (t460 * (-t132 * t216 - t213 * t79) + t419 * (-t132 * t213 + t216 * t79)) * g(2) + (t460 * (-t134 * t216 - t213 * t81) + t419 * (-t134 * t213 + t216 * t81)) * g(1) - t13 * t100 + (-t236 * t416 + t42) * pkin(5); -t154 * t100 - t236 * t102 + (-g(1) * t82 - g(2) * t80 - g(3) * t139 - t14 * t154 - t236 * t8 + t5) * m(7) + t19;];
tau  = t6;
