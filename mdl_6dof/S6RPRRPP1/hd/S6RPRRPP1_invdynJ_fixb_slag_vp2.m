% Calculate vector of inverse dynamics joint torques for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRRPP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:54
% EndTime: 2019-03-09 04:28:33
% DurationCPUTime: 25.71s
% Computational Cost: add. (7331->676), mult. (15488->870), div. (0->0), fcn. (10039->14), ass. (0->297)
t241 = sin(pkin(9));
t224 = pkin(1) * t241 + pkin(7);
t203 = t224 * qJD(1);
t245 = sin(qJ(3));
t248 = cos(qJ(3));
t152 = qJD(2) * t248 - t245 * t203;
t139 = -qJD(3) * pkin(3) - t152;
t244 = sin(qJ(4));
t247 = cos(qJ(4));
t328 = qJD(3) * t247;
t331 = qJD(1) * t245;
t187 = -t244 * t331 + t328;
t114 = -pkin(4) * t187 + qJD(5) + t139;
t330 = qJD(1) * t248;
t217 = qJD(4) - t330;
t350 = cos(pkin(10));
t240 = sin(pkin(10));
t153 = t245 * qJD(2) + t248 * t203;
t140 = qJD(3) * pkin(8) + t153;
t291 = pkin(3) * t248 + pkin(8) * t245;
t272 = -pkin(2) - t291;
t242 = cos(pkin(9));
t376 = pkin(1) * t242;
t174 = t272 - t376;
t143 = t174 * qJD(1);
t85 = t140 * t247 + t143 * t244;
t74 = qJ(5) * t187 + t85;
t357 = t240 * t74;
t188 = qJD(3) * t244 + t247 * t331;
t84 = -t140 * t244 + t247 * t143;
t73 = -qJ(5) * t188 + t84;
t68 = pkin(4) * t217 + t73;
t22 = t350 * t68 - t357;
t20 = -t217 * pkin(5) + qJD(6) - t22;
t120 = -t350 * t187 + t188 * t240;
t266 = t240 * t187 + t188 * t350;
t31 = pkin(5) * t120 - qJ(6) * t266 + t114;
t395 = -t120 / 0.2e1;
t444 = Ifges(6,4) - Ifges(7,5);
t382 = t217 / 0.2e1;
t391 = t266 / 0.2e1;
t443 = Ifges(7,4) + Ifges(6,5);
t445 = Ifges(6,1) + Ifges(7,1);
t481 = t443 * t382 + t445 * t391;
t484 = t114 * mrSges(6,2) + mrSges(7,2) * t20 - mrSges(6,3) * t22 - t31 * mrSges(7,3) + t395 * t444 + t481;
t383 = -t217 / 0.2e1;
t392 = -t266 / 0.2e1;
t394 = t120 / 0.2e1;
t483 = Ifges(6,4) * t394 + Ifges(7,5) * t395 + t383 * t443 + t392 * t445 - t484;
t482 = Ifges(6,4) * t395 + Ifges(7,5) * t394 + t481 + t484;
t480 = mrSges(6,1) + mrSges(7,1);
t479 = -mrSges(6,2) + mrSges(7,3);
t442 = Ifges(7,2) + Ifges(6,3);
t465 = Ifges(5,3) + t442;
t290 = pkin(3) * t245 - pkin(8) * t248;
t190 = t290 * qJD(1);
t111 = t247 * t152 + t244 * t190;
t243 = -qJ(5) - pkin(8);
t296 = qJD(4) * t243;
t310 = t244 * t330;
t323 = qJD(5) * t247;
t478 = qJ(5) * t310 + t244 * t296 - t111 + t323;
t110 = -t152 * t244 + t247 * t190;
t334 = t247 * t248;
t270 = pkin(4) * t245 - qJ(5) * t334;
t477 = -qJD(1) * t270 - qJD(5) * t244 + t247 * t296 - t110;
t229 = pkin(4) * t247 + pkin(3);
t213 = t248 * t229;
t339 = t243 * t245;
t398 = m(6) + m(7);
t469 = -mrSges(6,3) - mrSges(7,2);
t422 = t469 * t245;
t475 = t422 - t398 * (t213 - t339);
t71 = t350 * t74;
t23 = t240 * t68 + t71;
t21 = qJ(6) * t217 + t23;
t474 = Ifges(6,4) * t391 + Ifges(7,5) * t392 + Ifges(6,6) * t382 + Ifges(7,6) * t383 + (Ifges(6,2) + Ifges(7,3)) * t395 - mrSges(6,1) * t114 - mrSges(7,1) * t31 + mrSges(7,2) * t21 + mrSges(6,3) * t23;
t322 = qJD(1) * qJD(3);
t192 = qJDD(1) * t248 - t245 * t322;
t179 = qJDD(4) - t192;
t387 = t179 / 0.2e1;
t193 = qJDD(1) * t245 + t248 * t322;
t117 = qJD(4) * t187 + qJDD(3) * t244 + t193 * t247;
t118 = -qJD(4) * t188 + qJDD(3) * t247 - t193 * t244;
t64 = t117 * t350 + t240 * t118;
t400 = t64 / 0.2e1;
t472 = t443 * t387 + t445 * t400;
t63 = t117 * t240 - t118 * t350;
t402 = -t63 / 0.2e1;
t226 = -pkin(2) - t376;
t202 = t226 * qJDD(1);
t115 = -pkin(3) * t192 - pkin(8) * t193 + t202;
t324 = qJD(4) * t247;
t326 = qJD(4) * t244;
t329 = qJD(3) * t245;
t458 = qJD(2) * qJD(3) + t224 * qJDD(1);
t103 = t245 * qJDD(2) - t203 * t329 + t248 * t458;
t93 = qJDD(3) * pkin(8) + t103;
t24 = t244 * t115 - t140 * t326 + t143 * t324 + t247 * t93;
t471 = t24 * mrSges(5,2);
t25 = -qJD(4) * t85 + t247 * t115 - t244 * t93;
t470 = t25 * mrSges(5,1);
t441 = Ifges(6,6) - Ifges(7,6);
t468 = pkin(4) * t398 + mrSges(5,1);
t466 = Ifges(6,2) * t395 - Ifges(7,3) * t394 + t444 * t391 + t474;
t238 = qJ(4) + pkin(10);
t231 = sin(t238);
t233 = cos(t238);
t286 = -mrSges(5,1) * t247 + mrSges(5,2) * t244;
t464 = m(5) * pkin(3) + t479 * t231 + t233 * t480 - t286;
t363 = Ifges(4,4) * t245;
t280 = t248 * Ifges(4,2) + t363;
t462 = t20 * mrSges(7,1) + t23 * mrSges(6,2) + Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t280 / 0.2e1 - t21 * mrSges(7,3) - t22 * mrSges(6,1);
t460 = t402 * t444 + t472;
t434 = -t240 * t478 + t350 * t477;
t431 = t240 * t477 + t350 * t478;
t327 = qJD(3) * t248;
t261 = t244 * t327 + t245 * t324;
t239 = qJ(1) + pkin(9);
t232 = sin(t239);
t234 = cos(t239);
t337 = t244 * t248;
t147 = t232 * t247 - t234 * t337;
t425 = -t153 + (-t310 + t326) * pkin(4);
t457 = g(1) * t234 + g(2) * t232;
t453 = -Ifges(6,2) * t394 + Ifges(7,3) * t395 - t383 * t441 - t392 * t444 + t474;
t401 = t63 / 0.2e1;
t449 = -m(5) - m(4);
t397 = t117 / 0.2e1;
t396 = t118 / 0.2e1;
t447 = -qJD(1) / 0.2e1;
t446 = mrSges(3,2) - mrSges(4,3);
t18 = t63 * mrSges(7,1) - t64 * mrSges(7,3);
t19 = t63 * mrSges(6,1) + t64 * mrSges(6,2);
t439 = -t18 - t19;
t34 = -mrSges(7,2) * t63 + mrSges(7,3) * t179;
t35 = -mrSges(6,2) * t179 - mrSges(6,3) * t63;
t438 = t34 + t35;
t36 = mrSges(6,1) * t179 - mrSges(6,3) * t64;
t37 = -t179 * mrSges(7,1) + t64 * mrSges(7,2);
t437 = -t36 + t37;
t295 = t350 * t244;
t181 = t240 * t247 + t295;
t141 = t181 * t330;
t294 = t350 * t247;
t340 = t240 * t244;
t265 = t294 - t340;
t259 = t248 * t265;
t142 = qJD(1) * t259;
t162 = t181 * qJD(4);
t163 = t265 * qJD(4);
t435 = -qJD(6) * t181 + t425 + (t142 - t163) * qJ(6) + (-t141 + t162) * pkin(5);
t89 = -mrSges(6,2) * t217 - mrSges(6,3) * t120;
t92 = -mrSges(7,2) * t120 + mrSges(7,3) * t217;
t367 = t89 + t92;
t90 = mrSges(6,1) * t217 - mrSges(6,3) * t266;
t91 = -mrSges(7,1) * t217 + mrSges(7,2) * t266;
t366 = t90 - t91;
t433 = pkin(5) * t331 - t434;
t432 = -qJ(6) * t331 + t431;
t66 = -mrSges(5,1) * t118 + mrSges(5,2) * t117;
t430 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t193 + t66;
t429 = t245 * t457;
t316 = mrSges(4,3) * t331;
t428 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t187 + mrSges(5,2) * t188 + t316;
t175 = Ifges(5,4) * t187;
t108 = t188 * Ifges(5,1) + t217 * Ifges(5,5) + t175;
t230 = Ifges(4,4) * t330;
t427 = Ifges(4,1) * t331 + Ifges(4,5) * qJD(3) + t247 * t108 + t230;
t189 = t224 * t334;
t124 = t244 * t174 + t189;
t325 = qJD(4) * t245;
t426 = t244 * t325 - t247 * t327;
t104 = qJDD(2) * t248 - t203 * t327 - t245 * t458;
t424 = t103 * t248 - t104 * t245;
t80 = mrSges(5,1) * t179 - mrSges(5,3) * t117;
t81 = -mrSges(5,2) * t179 + mrSges(5,3) * t118;
t423 = -t244 * t80 + t247 * t81;
t421 = t24 * t247 - t244 * t25;
t420 = t188 * Ifges(5,5) + t187 * Ifges(5,6) - t120 * t441 + t217 * t465 + t266 * t443;
t288 = mrSges(4,1) * t248 - mrSges(4,2) * t245;
t419 = t245 * mrSges(5,3) + mrSges(3,1) + t288;
t418 = Ifges(5,5) * t117 + Ifges(5,6) * t118 + t179 * t465 - t441 * t63 + t443 * t64;
t417 = m(7) * pkin(5) + t480;
t415 = m(7) * qJ(6) + t479;
t11 = qJ(5) * t118 + qJD(5) * t187 + t24;
t7 = pkin(4) * t179 - qJ(5) * t117 - qJD(5) * t188 + t25;
t4 = t350 * t11 + t240 * t7;
t1 = qJ(6) * t179 + qJD(6) * t217 + t4;
t94 = -qJDD(3) * pkin(3) - t104;
t49 = -pkin(4) * t118 + qJDD(5) + t94;
t5 = pkin(5) * t63 - qJ(6) * t64 - qJD(6) * t266 + t49;
t409 = mrSges(6,1) * t49 + mrSges(7,1) * t5 - mrSges(7,2) * t1 - mrSges(6,3) * t4 + 0.2e1 * Ifges(7,3) * t401 - t64 * Ifges(6,4) / 0.2e1 - t179 * Ifges(6,6) / 0.2e1 + Ifges(7,6) * t387 + (-t444 + Ifges(7,5)) * t400 + (-t402 + t401) * Ifges(6,2);
t408 = m(6) * pkin(4);
t405 = Ifges(5,1) * t397 + Ifges(5,4) * t396 + Ifges(5,5) * t387;
t384 = t188 / 0.2e1;
t246 = sin(qJ(1));
t375 = pkin(1) * t246;
t374 = pkin(4) * t188;
t373 = pkin(4) * t240;
t370 = g(3) * t245;
t249 = cos(qJ(1));
t237 = t249 * pkin(1);
t191 = t290 * qJD(3);
t309 = t224 * t329;
t332 = t247 * t191 + t244 * t309;
t40 = -t245 * t323 + t270 * qJD(3) + (-t189 + (qJ(5) * t245 - t174) * t244) * qJD(4) + t332;
t333 = t174 * t324 + t244 * t191;
t336 = t245 * t247;
t47 = (-qJ(5) * qJD(4) - qJD(3) * t224) * t336 + (-qJD(5) * t245 + (-qJ(5) * qJD(3) - qJD(4) * t224) * t248) * t244 + t333;
t13 = t240 * t40 + t350 * t47;
t365 = mrSges(5,3) * t187;
t364 = mrSges(5,3) * t188;
t362 = Ifges(4,4) * t248;
t361 = Ifges(5,4) * t244;
t360 = Ifges(5,4) * t247;
t359 = t188 * Ifges(5,4);
t338 = t244 * t245;
t109 = -qJ(5) * t338 + t124;
t155 = t247 * t174;
t97 = -qJ(5) * t336 + t155 + (-t224 * t244 - pkin(4)) * t248;
t46 = t350 * t109 + t240 * t97;
t346 = t232 * t244;
t344 = t232 * t248;
t343 = t234 * t244;
t341 = t234 * t248;
t209 = t245 * t224;
t218 = pkin(4) * t338;
t159 = t209 + t218;
t317 = m(5) * pkin(8) + mrSges(5,3);
t315 = mrSges(4,3) * t330;
t195 = t224 * t327;
t126 = pkin(4) * t261 + t195;
t312 = t234 * pkin(2) + t232 * pkin(7) + t237;
t311 = t350 * pkin(4);
t107 = t187 * Ifges(5,2) + t217 * Ifges(5,6) + t359;
t304 = -t244 * t107 / 0.2e1;
t298 = t234 * pkin(7) - t375;
t287 = mrSges(4,1) * t245 + mrSges(4,2) * t248;
t285 = mrSges(5,1) * t244 + mrSges(5,2) * t247;
t282 = Ifges(5,1) * t247 - t361;
t281 = Ifges(5,1) * t244 + t360;
t279 = -Ifges(5,2) * t244 + t360;
t278 = Ifges(5,2) * t247 + t361;
t277 = Ifges(4,5) * t248 - Ifges(4,6) * t245;
t276 = Ifges(5,5) * t247 - Ifges(5,6) * t244;
t275 = Ifges(5,5) * t244 + Ifges(5,6) * t247;
t274 = pkin(5) * t233 + qJ(6) * t231;
t3 = -t240 * t11 + t350 * t7;
t145 = t232 * t337 + t234 * t247;
t12 = -t240 * t47 + t350 * t40;
t269 = t139 * t285;
t268 = t226 * qJD(1) * t287;
t267 = t245 * (Ifges(4,1) * t248 - t363);
t45 = -t240 * t109 + t350 * t97;
t257 = Ifges(5,5) * t245 + t248 * t282;
t256 = Ifges(5,6) * t245 + t248 * t279;
t255 = Ifges(5,3) * t245 + t248 * t276;
t254 = (-t244 * t85 - t247 * t84) * qJD(4) + t421;
t225 = -t311 - pkin(5);
t219 = qJ(6) + t373;
t208 = t243 * t247;
t206 = -qJD(3) * mrSges(4,2) + t315;
t173 = t285 * t245;
t156 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t192;
t151 = t265 * t245;
t150 = t181 * t245;
t148 = t234 * t334 + t346;
t146 = -t232 * t334 + t343;
t137 = mrSges(5,1) * t217 - t364;
t136 = -mrSges(5,2) * t217 + t365;
t134 = t231 * t232 + t233 * t341;
t133 = t231 * t341 - t232 * t233;
t132 = -t234 * t231 + t233 * t344;
t131 = t231 * t344 + t233 * t234;
t128 = -t208 * t350 + t243 * t340;
t127 = -t208 * t240 - t243 * t295;
t123 = -t224 * t337 + t155;
t116 = -pkin(5) * t265 - qJ(6) * t181 - t229;
t96 = qJD(3) * t259 - t181 * t325;
t95 = t240 * t426 - t294 * t325 - t295 * t327;
t77 = pkin(5) * t150 - qJ(6) * t151 + t159;
t76 = -qJD(4) * t124 + t332;
t75 = (-t245 * t328 - t248 * t326) * t224 + t333;
t70 = mrSges(6,1) * t120 + mrSges(6,2) * t266;
t69 = mrSges(7,1) * t120 - mrSges(7,3) * t266;
t48 = pkin(5) * t266 + qJ(6) * t120 + t374;
t42 = t117 * Ifges(5,4) + t118 * Ifges(5,2) + t179 * Ifges(5,6);
t41 = t248 * pkin(5) - t45;
t39 = -qJ(6) * t248 + t46;
t28 = -pkin(5) * t95 - qJ(6) * t96 - qJD(6) * t151 + t126;
t27 = t350 * t73 - t357;
t26 = t240 * t73 + t71;
t10 = -pkin(5) * t329 - t12;
t8 = qJ(6) * t329 - qJD(6) * t248 + t13;
t2 = -t179 * pkin(5) + qJDD(6) - t3;
t6 = [t430 * t209 - t202 * t288 - t418 * t248 / 0.2e1 + (qJD(3) * t255 - t275 * t325 + t441 * t95) * t382 + (-t1 * t248 - t151 * t5) * mrSges(7,3) + m(7) * (t1 * t39 + t10 * t20 + t2 * t41 + t21 * t8 + t28 * t31 + t5 * t77) + m(6) * (t114 * t126 + t12 * t22 + t13 * t23 + t159 * t49 + t3 * t45 + t4 * t46) + m(4) * (t202 * t226 + ((-t152 * t248 - t153 * t245) * qJD(3) + t424) * t224) + t427 * t327 / 0.2e1 + t428 * t195 + m(5) * (t123 * t25 + t124 * t24 + t75 * t85 + t76 * t84 + (t139 * t327 + t245 * t94) * t224) + t151 * t460 + t187 * (qJD(3) * t256 - t278 * t325) / 0.2e1 + (t151 * t49 + t248 * t4) * mrSges(6,2) + (Ifges(7,5) * t151 - Ifges(7,6) * t248) * t401 + (Ifges(6,4) * t151 - Ifges(6,6) * t248) * t402 + (0.2e1 * (mrSges(3,1) * t242 - mrSges(3,2) * t241) * pkin(1) + m(3) * (t241 ^ 2 + t242 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + t482 * t96 + t409 * t150 + (t363 + t280) * t192 / 0.2e1 + (-t24 * t338 - t25 * t336 - t261 * t85 + t426 * t84) * mrSges(5,3) + t139 * (mrSges(5,1) * t261 - mrSges(5,2) * t426) + t248 * t471 + t248 * (Ifges(4,4) * t193 + Ifges(4,2) * t192) / 0.2e1 + (-m(3) * t237 - mrSges(2,1) * t249 - t148 * mrSges(5,1) + mrSges(2,2) * t246 - t147 * mrSges(5,2) + t449 * t312 - t398 * (pkin(4) * t346 + t312) + t446 * t232 - t417 * t134 - t415 * t133 + (-m(5) * t291 - t419 + t475) * t234) * g(2) + t336 * t405 + (-Ifges(5,6) * t248 + t245 * t279) * t396 + (-Ifges(5,5) * t248 + t245 * t282) * t397 + (qJD(3) * t257 - t281 * t325) * t384 + (-t441 * t150 + t443 * t151 + t245 * t276 - t248 * t465) * t387 + qJD(3) ^ 2 * t277 / 0.2e1 - t42 * t338 / 0.2e1 + t304 * t327 + (t267 + t248 * (-Ifges(4,2) * t245 + t362)) * t322 / 0.2e1 - (t247 * t107 + t244 * t108) * t325 / 0.2e1 + (-t152 * t327 + t424) * mrSges(4,3) + t466 * t95 + qJD(3) * t268 + (t442 * t382 + t420 / 0.2e1 - t153 * mrSges(4,3) - t85 * mrSges(5,2) + Ifges(7,6) * t394 + Ifges(6,6) * t395 + t84 * mrSges(5,1) + t443 * t391 - t462) * t329 - t248 * t470 + t2 * (mrSges(7,1) * t248 + mrSges(7,2) * t151) + t3 * (-mrSges(6,1) * t248 - mrSges(6,3) * t151) + qJDD(3) * (Ifges(4,5) * t245 + Ifges(4,6) * t248) + t193 * t362 / 0.2e1 - t206 * t309 + t226 * (-mrSges(4,1) * t192 + mrSges(4,2) * t193) + t94 * t173 + t159 * t19 + t75 * t136 + t76 * t137 + t126 * t70 + t123 * t80 + t124 * t81 + t13 * t89 + t12 * t90 + t10 * t91 + t8 * t92 + t193 * t245 * Ifges(4,1) + (t151 * t445 - t248 * t443) * t400 + (m(3) * t375 + mrSges(2,1) * t246 - t146 * mrSges(5,1) + mrSges(2,2) * t249 - t145 * mrSges(5,2) + t449 * t298 - t398 * (pkin(4) * t343 + t232 * t339 + t298) + t446 * t234 + t417 * t132 + t415 * t131 + (m(4) * pkin(2) - m(5) * t272 - t398 * (-pkin(2) - t213) + t419 - t422) * t232) * g(1) + t39 * t34 + t41 * t37 + t45 * t36 + t46 * t35 + t28 * t69 + t77 * t18 + t248 * t224 * t156; m(3) * qJDD(2) + t367 * t96 + t366 * t95 + t438 * t151 + t437 * t150 + (-m(3) - t398 + t449) * g(3) + ((t136 * t247 - t137 * t244 + t206) * qJD(3) - t430 + t439) * t248 + (t156 + (-t244 * t136 - t247 * t137) * qJD(4) + (t69 + t70 + t428) * qJD(3) + t423) * t245 + m(7) * (t1 * t151 + t150 * t2 - t20 * t95 + t21 * t96 - t248 * t5 + t31 * t329) + m(6) * (t114 * t329 - t150 * t3 + t151 * t4 + t22 * t95 + t23 * t96 - t248 * t49) + m(4) * (t103 * t245 + t104 * t248 + (-t152 * t245 + t153 * t248) * qJD(3)) + m(5) * ((-t94 + (-t244 * t84 + t247 * t85) * qJD(3)) * t248 + (qJD(3) * t139 + t254) * t245); (t304 + t269) * qJD(4) + t431 * t89 + t432 * t92 + t433 * t91 + (t114 * t425 - t127 * t3 + t128 * t4 + t22 * t434 - t229 * t49 + t23 * t431) * m(6) + t434 * t90 + (t1 * t128 + t116 * t5 + t127 * t2 + t20 * t433 + t21 * t432 + t31 * t435) * m(7) + t435 * t69 + t437 * t127 + t438 * t128 - (-t441 * t387 + t409) * t265 + (-t324 * t84 - t326 * t85 + t421) * mrSges(5,3) + (m(5) * t254 - t136 * t326 - t137 * t324 + t423) * pkin(8) - t420 * t331 / 0.2e1 + t425 * t70 - (-Ifges(4,2) * t331 + t230 + t427) * t330 / 0.2e1 + (-m(5) * t139 + t316 - t428) * t153 + t108 * t324 / 0.2e1 - t277 * t322 / 0.2e1 + t482 * t163 + t483 * t142 + (-pkin(3) * t94 - t110 * t84 - t111 * t85) * m(5) + (t187 * t256 + t188 * t257 + t217 * t255) * t447 + t453 * t141 + (-t245 * t317 - t288 + (-m(7) * t274 - t464) * t248 + t475) * g(3) + t244 * t405 + t275 * t387 + t278 * t396 + t281 * t397 - t269 * t330 + (t187 * t279 + t188 * t282 + t217 * t276) * qJD(4) / 0.2e1 + (-t441 * t382 - t466) * t162 + (t315 - t206) * t152 + (Ifges(6,6) * t394 + Ifges(7,6) * t395 + t442 * t383 + t443 * t392 + t462) * t331 + t457 * (t287 + (t243 * t398 - t317 + t469) * t248 + (m(6) * t229 - m(7) * (-t229 - t274) + t464) * t245) + t247 * t42 / 0.2e1 + t107 * t310 / 0.2e1 - t229 * t19 + Ifges(4,6) * t192 + Ifges(4,5) * t193 - t111 * t136 - t110 * t137 + t116 * t18 - t103 * mrSges(4,2) + t104 * mrSges(4,1) + (mrSges(6,2) * t49 + mrSges(7,2) * t2 - mrSges(6,3) * t3 - mrSges(7,3) * t5 + Ifges(6,4) * t402 + Ifges(7,5) * t401 + t460 + t472) * t181 + (-t85 * (-mrSges(5,2) * t245 - mrSges(5,3) * t337) - t84 * (mrSges(5,1) * t245 - mrSges(5,3) * t334) - t268 + t267 * t447) * qJD(1) + t94 * t286 + Ifges(4,3) * qJDD(3) - pkin(3) * t66; (m(6) * t22 + t366) * t26 + (-m(6) * t23 - t367) * t27 + t418 + (-m(6) * t114 - t70) * t374 - t483 * t120 + (-mrSges(5,2) * t146 + t131 * t417 - t132 * t415 + t145 * t468) * g(2) + (mrSges(5,2) * t148 + t417 * t133 - t415 * t134 - t147 * t468) * g(1) - (-mrSges(6,1) * t231 - mrSges(6,2) * t233 - t244 * t408) * t370 + t453 * t266 - t188 * (Ifges(5,1) * t187 - t359) / 0.2e1 + (t240 * t4 + t3 * t350) * t408 + t35 * t373 + (Ifges(5,5) * t187 - Ifges(5,6) * t188) * t383 + t107 * t384 + t36 * t311 - (-Ifges(5,2) * t188 + t108 + t175) * t187 / 0.2e1 + (-(-t231 * mrSges(7,1) + t233 * mrSges(7,3)) * t245 + t173) * g(3) + (-g(3) * (-t218 + (-pkin(5) * t231 + qJ(6) * t233) * t245) - t20 * t26 - t31 * t48 + t1 * t219 + t2 * t225 + (qJD(6) - t27) * t21) * m(7) - t471 + t470 + t225 * t37 + t219 * t34 - t139 * (mrSges(5,1) * t188 + mrSges(5,2) * t187) + qJD(6) * t92 + (t365 - t136) * t84 + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2) + (t364 + t137) * t85 - t48 * t69; t398 * t248 * g(3) + t366 * t266 + t367 * t120 + (t120 * t21 - t20 * t266 - t429 + t5) * m(7) + (t120 * t23 + t22 * t266 - t429 + t49) * m(6) - t439; t266 * t69 - t217 * t92 + (-g(1) * t133 - g(2) * t131 - t21 * t217 - t231 * t370 + t266 * t31 + t2) * m(7) + t37;];
tau  = t6;
