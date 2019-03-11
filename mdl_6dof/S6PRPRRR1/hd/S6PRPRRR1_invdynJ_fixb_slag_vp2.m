% Calculate vector of inverse dynamics joint torques for
% S6PRPRRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:25
% EndTime: 2019-03-08 20:22:51
% DurationCPUTime: 13.90s
% Computational Cost: add. (8171->608), mult. (18134->856), div. (0->0), fcn. (14236->16), ass. (0->287)
t446 = mrSges(6,2) - mrSges(7,3);
t249 = sin(qJ(6));
t253 = cos(qJ(6));
t296 = -mrSges(7,1) * t253 + mrSges(7,2) * t249;
t445 = -mrSges(6,1) + t296;
t250 = sin(qJ(5));
t251 = sin(qJ(4));
t254 = cos(qJ(5));
t255 = cos(qJ(4));
t205 = t250 * t255 + t251 * t254;
t198 = t205 * qJD(2);
t241 = qJD(4) + qJD(5);
t167 = -t198 * t249 + t241 * t253;
t168 = t198 * t253 + t241 * t249;
t385 = mrSges(6,3) * t198;
t373 = -mrSges(6,1) * t241 - mrSges(7,1) * t167 + mrSges(7,2) * t168 + t385;
t243 = sin(pkin(12));
t245 = sin(pkin(6));
t252 = sin(qJ(2));
t350 = t245 * t252;
t313 = qJD(1) * t350;
t214 = t243 * t313;
t246 = cos(pkin(12));
t256 = cos(qJ(2));
t338 = qJD(1) * t256;
t312 = t245 * t338;
t180 = t246 * t312 - t214;
t230 = pkin(2) * t243 + pkin(8);
t392 = pkin(9) + t230;
t307 = qJD(4) * t392;
t192 = t251 * t307;
t204 = t250 * t251 - t254 * t255;
t200 = t392 * t251;
t201 = t392 * t255;
t283 = -t200 * t254 - t201 * t250;
t300 = t255 * t307;
t432 = qJD(5) * t283 + t180 * t204 - t192 * t254 - t250 * t300;
t269 = t204 * qJD(5);
t157 = -qJD(4) * t204 - t269;
t270 = t205 * qJD(5);
t158 = qJD(4) * t205 + t270;
t281 = t243 * t256 + t246 * t252;
t186 = t281 * t245;
t177 = qJD(1) * t186;
t335 = qJD(4) * t251;
t326 = pkin(4) * t335;
t444 = pkin(5) * t158 - pkin(10) * t157 - t177 + t326;
t437 = -m(7) - m(6);
t438 = -m(4) - m(5);
t443 = t437 + t438;
t242 = qJ(4) + qJ(5);
t238 = sin(t242);
t239 = cos(t242);
t298 = -mrSges(5,1) * t255 + mrSges(5,2) * t251;
t412 = m(7) * pkin(10);
t413 = m(7) * pkin(5);
t442 = m(5) * pkin(3) - t298 + (t413 - t445) * t239 + (t412 - t446) * t238;
t329 = qJD(2) * qJD(4);
t208 = qJDD(2) * t255 - t251 * t329;
t209 = qJDD(2) * t251 + t255 * t329;
t115 = -qJD(2) * t270 + t208 * t254 - t209 * t250;
t244 = sin(pkin(11));
t247 = cos(pkin(11));
t248 = cos(pkin(6));
t343 = t248 * t256;
t441 = -t244 * t252 + t247 * t343;
t341 = t256 * t246;
t202 = t243 * t252 - t341;
t345 = t248 * t252;
t339 = -t243 * t343 - t246 * t345;
t131 = -t202 * t247 + t244 * t339;
t349 = t245 * t255;
t440 = -t131 * t251 + t244 * t349;
t159 = -t186 * t251 + t248 * t255;
t114 = -qJD(2) * t269 + t208 * t250 + t209 * t254;
t240 = qJDD(4) + qJDD(5);
t54 = qJD(6) * t167 + t114 * t253 + t240 * t249;
t410 = t54 / 0.2e1;
t55 = -qJD(6) * t168 - t114 * t249 + t240 * t253;
t409 = t55 / 0.2e1;
t197 = t204 * qJD(2);
t194 = qJD(6) + t197;
t380 = t168 * Ifges(7,4);
t68 = t167 * Ifges(7,2) + t194 * Ifges(7,6) + t380;
t439 = -t68 / 0.2e1;
t109 = qJDD(6) - t115;
t407 = t109 / 0.2e1;
t436 = t208 / 0.2e1;
t22 = -mrSges(7,1) * t55 + mrSges(7,2) * t54;
t99 = mrSges(6,1) * t240 - mrSges(6,3) * t114;
t391 = t22 - t99;
t235 = pkin(4) * t255 + pkin(3);
t398 = pkin(2) * t246;
t216 = -t235 - t398;
t139 = pkin(5) * t204 - pkin(10) * t205 + t216;
t144 = -t200 * t250 + t201 * t254;
t61 = t139 * t249 + t144 * t253;
t435 = -qJD(6) * t61 - t249 * t432 + t253 * t444;
t60 = t139 * t253 - t144 * t249;
t434 = qJD(6) * t60 + t249 * t444 + t253 * t432;
t295 = mrSges(7,1) * t249 + mrSges(7,2) * t253;
t226 = qJD(1) * t248 + qJD(3);
t215 = t255 * t226;
t210 = qJD(2) * pkin(2) + t312;
t166 = t210 * t243 + t246 * t313;
t162 = qJD(2) * pkin(8) + t166;
t305 = pkin(9) * qJD(2) + t162;
t112 = -t251 * t305 + t215;
t103 = qJD(4) * pkin(4) + t112;
t356 = t226 * t251;
t113 = t255 * t305 + t356;
t369 = t113 * t250;
t40 = t103 * t254 - t369;
t36 = -pkin(5) * t241 - t40;
t431 = t295 * t36;
t185 = t243 * t350 - t245 * t341;
t175 = t185 * t253;
t160 = t186 * t255 + t248 * t251;
t73 = t159 * t250 + t160 * t254;
t49 = -t249 * t73 + t175;
t428 = (mrSges(3,1) * t256 - mrSges(3,2) * t252) * t245 - t185 * mrSges(4,1) - t186 * mrSges(4,2);
t153 = -t186 * t238 + t239 * t248;
t154 = t186 * t239 + t238 * t248;
t427 = t153 * t445 + t154 * t446;
t354 = t244 * t245;
t91 = -t131 * t238 + t239 * t354;
t92 = t131 * t239 + t238 * t354;
t426 = t445 * t91 + t446 * t92;
t126 = t202 * t244 + t247 * t339;
t352 = t245 * t247;
t89 = t126 * t238 - t239 * t352;
t90 = -t126 * t239 - t238 * t352;
t425 = t445 * t89 + t446 * t90;
t330 = qJD(6) * t253;
t277 = t157 * t249 + t205 * t330;
t337 = qJD(2) * t251;
t321 = mrSges(5,3) * t337;
t219 = qJD(4) * mrSges(5,1) - t321;
t336 = qJD(2) * t255;
t320 = mrSges(5,3) * t336;
t220 = -qJD(4) * mrSges(5,2) + t320;
t424 = -t219 * t251 + t220 * t255;
t304 = qJD(2) * t313;
t348 = t245 * t256;
t187 = qJDD(1) * t348 - t304;
t370 = qJDD(2) * pkin(2);
t183 = t187 + t370;
t310 = qJD(2) * t338;
t188 = (qJDD(1) * t252 + t310) * t245;
t121 = t183 * t243 + t188 * t246;
t111 = qJDD(2) * pkin(8) + t121;
t225 = qJDD(1) * t248 + qJDD(3);
t334 = qJD(4) * t255;
t38 = t111 * t255 - t162 * t335 + t225 * t251 + t226 * t334;
t123 = t162 * t255 + t356;
t39 = -qJD(4) * t123 - t111 * t251 + t225 * t255;
t423 = -t251 * t39 + t255 * t38;
t28 = mrSges(7,1) * t109 - mrSges(7,3) * t54;
t29 = -mrSges(7,2) * t109 + mrSges(7,3) * t55;
t422 = -t249 * t28 + t253 * t29;
t150 = mrSges(6,1) * t197 + mrSges(6,2) * t198;
t420 = t150 + (-mrSges(4,1) + t298) * qJD(2);
t116 = -mrSges(7,2) * t194 + mrSges(7,3) * t167;
t117 = mrSges(7,1) * t194 - mrSges(7,3) * t168;
t386 = mrSges(6,3) * t197;
t171 = -mrSges(6,2) * t241 - t386;
t419 = t116 * t253 - t117 * t249 + t171;
t33 = qJDD(4) * pkin(4) - pkin(9) * t209 + t39;
t35 = pkin(9) * t208 + t38;
t342 = t254 * t113;
t41 = t103 * t250 + t342;
t9 = -qJD(5) * t41 - t250 * t35 + t254 * t33;
t418 = -mrSges(4,1) - t442;
t37 = pkin(10) * t241 + t41;
t165 = t210 * t246 - t214;
t148 = -qJD(2) * t235 - t165;
t71 = pkin(5) * t197 - pkin(10) * t198 + t148;
t20 = -t249 * t37 + t253 * t71;
t120 = t183 * t246 - t188 * t243;
t110 = -qJDD(2) * pkin(3) - t120;
t74 = -pkin(4) * t208 + t110;
t23 = -pkin(5) * t115 - pkin(10) * t114 + t74;
t332 = qJD(5) * t254;
t333 = qJD(5) * t250;
t8 = t103 * t332 - t113 * t333 + t250 * t33 + t254 * t35;
t5 = pkin(10) * t240 + t8;
t2 = qJD(6) * t20 + t23 * t249 + t253 * t5;
t21 = t249 * t71 + t253 * t37;
t3 = -qJD(6) * t21 + t23 * t253 - t249 * t5;
t417 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t416 = -m(5) * pkin(8) - mrSges(5,3) - mrSges(6,3) - t295;
t415 = -mrSges(4,2) - t416;
t414 = qJD(2) ^ 2;
t411 = Ifges(7,1) * t410 + Ifges(7,4) * t409 + Ifges(7,5) * t407;
t406 = -t167 / 0.2e1;
t405 = -t168 / 0.2e1;
t404 = t168 / 0.2e1;
t403 = -t194 / 0.2e1;
t400 = t198 / 0.2e1;
t399 = t253 / 0.2e1;
t397 = pkin(4) * t250;
t396 = pkin(4) * t254;
t395 = t2 * t253;
t394 = t249 * t3;
t384 = Ifges(5,4) * t251;
t383 = Ifges(5,4) * t255;
t382 = Ifges(7,4) * t249;
t381 = Ifges(7,4) * t253;
t379 = t198 * Ifges(6,4);
t371 = mrSges(4,2) * qJD(2);
t364 = t185 * t249;
t362 = t197 * t249;
t361 = t197 * t253;
t360 = t205 * t249;
t359 = t205 * t253;
t351 = t245 * t251;
t331 = qJD(6) * t249;
t227 = pkin(2) * t348;
t328 = Ifges(7,5) * t54 + Ifges(7,6) * t55 + Ifges(7,3) * t109;
t327 = pkin(4) * t337;
t163 = Ifges(7,4) * t167;
t69 = Ifges(7,1) * t168 + t194 * Ifges(7,5) + t163;
t316 = t69 * t399;
t314 = pkin(5) * t91 + pkin(10) * t92;
t308 = -t331 / 0.2e1;
t303 = t441 * pkin(2);
t302 = t440 * pkin(4);
t301 = t159 * pkin(4);
t152 = pkin(5) * t198 + pkin(10) * t197;
t294 = Ifges(7,1) * t253 - t382;
t293 = t255 * Ifges(5,2) + t384;
t292 = -Ifges(7,2) * t249 + t381;
t291 = Ifges(5,5) * t255 - Ifges(5,6) * t251;
t290 = Ifges(7,5) * t253 - Ifges(7,6) * t249;
t289 = -t20 * t249 + t21 * t253;
t50 = t253 * t73 + t364;
t122 = -t162 * t251 + t215;
t285 = -t122 * t251 + t123 * t255;
t284 = t159 * t254 - t160 * t250;
t279 = t126 * t251 - t247 * t349;
t276 = -t157 * t253 + t205 * t331;
t161 = -qJD(2) * pkin(3) - t165;
t274 = t161 * (mrSges(5,1) * t251 + mrSges(5,2) * t255);
t273 = t251 * (Ifges(5,1) * t255 - t384);
t271 = t202 * t248;
t266 = t279 * pkin(4);
t264 = -t394 + (-t20 * t253 - t21 * t249) * qJD(6);
t262 = t264 + t395;
t261 = (-t116 * t249 - t117 * t253) * qJD(6) + t422;
t136 = -t197 * Ifges(6,2) + t241 * Ifges(6,6) + t379;
t193 = Ifges(6,4) * t197;
t137 = t198 * Ifges(6,1) + t241 * Ifges(6,5) - t193;
t16 = Ifges(7,4) * t54 + Ifges(7,2) * t55 + Ifges(7,6) * t109;
t6 = -pkin(5) * t240 - t9;
t67 = Ifges(7,5) * t168 + t167 * Ifges(7,6) + t194 * Ifges(7,3);
t259 = Ifges(6,5) * t114 - t148 * (mrSges(6,1) * t198 - mrSges(6,2) * t197) + (Ifges(7,3) * t198 - t197 * t290) * t403 + (Ifges(7,5) * t198 - t197 * t294) * t405 + (Ifges(7,6) * t198 - t197 * t292) * t406 - t241 * (-Ifges(6,5) * t197 - Ifges(6,6) * t198) / 0.2e1 - t20 * (mrSges(7,1) * t198 + mrSges(7,3) * t361) - t21 * (-mrSges(7,2) * t198 + mrSges(7,3) * t362) + (t316 + t431) * qJD(6) + t68 * t308 - t8 * mrSges(6,2) + t9 * mrSges(6,1) + (-Ifges(6,2) * t198 + t137 - t193) * t197 / 0.2e1 + (t167 * t292 + t168 * t294 + t194 * t290) * qJD(6) / 0.2e1 - (-Ifges(6,1) * t197 - t379 + t67) * t198 / 0.2e1 + t69 * t361 / 0.2e1 + (Ifges(7,5) * t249 + Ifges(7,6) * t253) * t407 + (Ifges(7,2) * t253 + t382) * t409 + (Ifges(7,1) * t249 + t381) * t410 + t249 * t411 + t362 * t439 + t16 * t399 + t136 * t400 + t197 * t431 - t40 * t386 + mrSges(7,3) * t395 + t6 * t296 + Ifges(6,3) * t240 + Ifges(6,6) * t115;
t257 = -pkin(9) - pkin(8);
t236 = Ifges(5,4) * t336;
t234 = -pkin(5) - t396;
t231 = -pkin(3) - t398;
t196 = Ifges(5,1) * t337 + Ifges(5,5) * qJD(4) + t236;
t195 = Ifges(5,6) * qJD(4) + qJD(2) * t293;
t191 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t209;
t190 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t208;
t179 = t202 * t245 * qJD(2);
t178 = qJD(2) * t186;
t164 = -mrSges(5,1) * t208 + mrSges(5,2) * t209;
t151 = t153 * pkin(5);
t135 = t152 + t327;
t130 = t244 * t271 - t247 * t281;
t127 = -t244 * t281 - t247 * t271;
t100 = -mrSges(6,2) * t240 + mrSges(6,3) * t115;
t88 = qJD(4) * t159 - t179 * t255;
t87 = -qJD(4) * t160 + t179 * t251;
t85 = t89 * pkin(5);
t44 = -mrSges(6,1) * t115 + mrSges(6,2) * t114;
t43 = t112 * t254 - t369;
t42 = t112 * t250 + t342;
t27 = t152 * t249 + t253 * t40;
t26 = t152 * t253 - t249 * t40;
t25 = t135 * t249 + t253 * t43;
t24 = t135 * t253 - t249 * t43;
t19 = qJD(5) * t73 + t250 * t88 - t254 * t87;
t18 = qJD(5) * t284 + t250 * t87 + t254 * t88;
t11 = -qJD(6) * t50 + t178 * t253 - t18 * t249;
t10 = qJD(6) * t49 + t178 * t249 + t18 * t253;
t1 = [t179 * t371 + m(2) * qJDD(1) + t10 * t116 + t73 * t100 + t11 * t117 + t159 * t191 + t160 * t190 + t18 * t171 + t87 * t219 + t88 * t220 + t49 * t28 + t50 * t29 - t391 * t284 + (-mrSges(3,1) * t252 - mrSges(3,2) * t256) * t414 * t245 + t373 * t19 + (t44 + t164) * t185 + t420 * t178 + t428 * qJDD(2) + (-m(2) - m(3) + t443) * g(3) + m(3) * (qJDD(1) * t248 ^ 2 + (t187 * t256 + t188 * t252) * t245) + m(4) * (-t120 * t185 + t121 * t186 - t165 * t178 - t166 * t179 + t225 * t248) + m(6) * (t148 * t178 + t18 * t41 + t185 * t74 - t19 * t40 + t284 * t9 + t73 * t8) + m(5) * (t110 * t185 + t122 * t87 + t123 * t88 + t159 * t39 + t160 * t38 + t161 * t178) + m(7) * (t10 * t21 + t11 * t20 + t19 * t36 + t2 * t50 - t284 * t6 + t3 * t49); (t144 * t8 + t148 * t326 + t216 * t74 + t283 * t9 + t41 * t432) * m(6) + (t273 + t255 * (-Ifges(5,2) * t251 + t383)) * t329 / 0.2e1 + (t2 * t61 + t20 * t435 + t21 * t434 - t283 * t6 + t3 * t60) * m(7) + (t245 * t310 - t188) * mrSges(3,2) + t255 * (Ifges(5,4) * t209 + Ifges(5,2) * t208) / 0.2e1 + (m(5) * ((-t122 * t255 - t123 * t251) * qJD(4) + t423) + t255 * t190 - t251 * t191 - t219 * t334 - t220 * t335) * t230 + t209 * t383 / 0.2e1 + m(4) * (t120 * t246 + t121 * t243) * pkin(2) + (-t441 * mrSges(3,1) - (-t244 * t256 - t247 * t345) * mrSges(3,2) + t438 * t303 + t437 * (t126 * t257 + t127 * t235 + t303) + t418 * t127 + t415 * t126) * g(2) + (t437 * (-t185 * t235 - t186 * t257 + t227) + t438 * t227 + t416 * t186 + t442 * t185 - t428) * g(3) + (-t2 * t360 + t20 * t276 - t21 * t277 - t3 * t359) * mrSges(7,3) + (t274 + t291 * qJD(4) / 0.2e1) * qJD(4) + (-t122 * t334 - t123 * t335 + t423) * mrSges(5,3) + (-m(4) * t166 - m(5) * t285 + t371 - t424) * t180 + t150 * t326 + (-t157 * t40 - t158 * t41) * mrSges(6,3) + (t304 + t187) * mrSges(3,1) + (-t243 * t370 - t121) * mrSges(4,2) + (m(6) * t40 - m(7) * t36 - t373) * (-qJD(5) * t144 + t180 * t205 + t192 * t250 - t254 * t300) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-(t244 * t345 - t247 * t256) * mrSges(3,2) + t437 * (t130 * t235 - t131 * t257) + t418 * t130 - t415 * t131 + (pkin(2) * t443 - mrSges(3,1)) * (-t244 * t343 - t247 * t252)) * g(1) + t20 * mrSges(7,1) * t158 - t21 * mrSges(7,2) * t158 + (t246 * t370 + t120) * mrSges(4,1) + (m(5) * t231 + t298) * t110 + (t328 / 0.2e1 - t8 * mrSges(6,3) + Ifges(7,6) * t409 + Ifges(7,5) * t410 + Ifges(7,3) * t407 + t74 * mrSges(6,1) - Ifges(6,4) * t114 - Ifges(6,2) * t115 - Ifges(6,6) * t240 + t417) * t204 + t157 * t316 + t60 * t28 + t61 * t29 + (t74 * mrSges(6,2) - t9 * mrSges(6,3) + Ifges(6,1) * t114 + Ifges(6,4) * t115 + Ifges(6,5) * t240 + t290 * t407 + t292 * t409 + t294 * t410 + t295 * t6 + t308 * t69) * t205 - t391 * t283 - t16 * t360 / 0.2e1 + t196 * t334 / 0.2e1 - t195 * t335 / 0.2e1 + t359 * t411 + t277 * t439 + (Ifges(6,1) * t157 - Ifges(6,4) * t158) * t400 + (-Ifges(7,1) * t276 - Ifges(7,4) * t277 + Ifges(7,5) * t158) * t404 + t293 * t436 + t144 * t100 + t157 * t137 / 0.2e1 + t158 * t67 / 0.2e1 - t158 * t136 / 0.2e1 + t148 * (mrSges(6,1) * t158 + mrSges(6,2) * t157) + t432 * t171 + t434 * t116 + t435 * t117 + (Ifges(5,1) * t209 + Ifges(5,4) * t436) * t251 - t197 * (Ifges(6,4) * t157 - Ifges(6,2) * t158) / 0.2e1 + t216 * t44 + (m(4) * t165 - m(5) * t161 - m(6) * t148 - t420) * t177 + t231 * t164 + t241 * (Ifges(6,5) * t157 - Ifges(6,6) * t158) / 0.2e1 + t36 * (mrSges(7,1) * t277 - mrSges(7,2) * t276) + t194 * (-Ifges(7,5) * t276 - Ifges(7,6) * t277 + Ifges(7,3) * t158) / 0.2e1 + t167 * (-Ifges(7,4) * t276 - Ifges(7,2) * t277 + Ifges(7,6) * t158) / 0.2e1 + qJDD(4) * (Ifges(5,5) * t251 + Ifges(5,6) * t255); t251 * t190 + t255 * t191 + t391 * t204 + t373 * t158 + t424 * qJD(4) + t419 * t157 + (t100 + t261) * t205 + m(5) * (qJD(4) * t285 + t251 * t38 + t255 * t39) + m(7) * (t157 * t289 + t158 * t36 + t204 * t6 + t205 * t262) + m(4) * t225 + m(6) * (t157 * t41 - t158 * t40 - t204 * t9 + t205 * t8) - (-t248 * g(3) + (-g(1) * t244 + g(2) * t247) * t245) * t443; -(-Ifges(5,2) * t337 + t196 + t236) * t336 / 0.2e1 + (t320 - t220) * t122 + (m(7) * t262 - t116 * t331 - t117 * t330 + t422) * (pkin(10) + t397) + (-m(6) * t302 - t440 * mrSges(5,1) - (-t131 * t255 - t244 * t351) * mrSges(5,2) - m(7) * (t302 + t314) + t426) * g(1) + (t321 + t219) * t123 + (-m(6) * t301 - m(7) * (pkin(10) * t154 + t151 + t301) - mrSges(5,1) * t159 + mrSges(5,2) * t160 + t427) * g(3) + (-t20 * t330 - t21 * t331 - t394) * mrSges(7,3) + t259 - t38 * mrSges(5,2) + t39 * mrSges(5,1) - t414 * t273 / 0.2e1 + Ifges(5,3) * qJDD(4) + (-m(6) * t266 - t279 * mrSges(5,1) - (t126 * t255 + t247 * t351) * mrSges(5,2) - m(7) * (pkin(10) * t90 + t266 + t85) + t425) * g(2) + t195 * t337 / 0.2e1 - t291 * t329 / 0.2e1 - t150 * t327 + t100 * t397 + t41 * t385 + t99 * t396 + t373 * (pkin(4) * t333 - t42) - t43 * t171 + Ifges(5,6) * t208 + Ifges(5,5) * t209 + t234 * t22 + ((t250 * t8 + t254 * t9 + (-t250 * t40 + t254 * t41) * qJD(5)) * pkin(4) - t148 * t327 + t40 * t42 - t41 * t43) * m(6) - qJD(2) * t274 - t25 * t116 - t24 * t117 + t419 * pkin(4) * t332 + (-t20 * t24 - t21 * t25 - t36 * t42 + t234 * t6 + (t250 * t36 + t254 * t289) * qJD(5) * pkin(4)) * m(7); ((-g(2) * t90 - g(3) * t154) * m(7) + t261) * pkin(10) - m(7) * (t20 * t26 + t21 * t27 + t36 * t41) + t264 * mrSges(7,3) + t259 - pkin(5) * t22 + t262 * t412 - t6 * t413 + (-t373 + t385) * t41 + (-m(7) * t85 + t425) * g(2) + (-m(7) * t314 + t426) * g(1) + (-m(7) * t151 + t427) * g(3) - t40 * t171 - t27 * t116 - t26 * t117; -t36 * (mrSges(7,1) * t168 + mrSges(7,2) * t167) + (Ifges(7,1) * t167 - t380) * t405 + t68 * t404 + (Ifges(7,5) * t167 - Ifges(7,6) * t168) * t403 + t21 * t117 - t20 * t116 - g(1) * ((-t130 * t253 - t249 * t92) * mrSges(7,1) + (t130 * t249 - t253 * t92) * mrSges(7,2)) - g(2) * ((-t127 * t253 - t249 * t90) * mrSges(7,1) + (t127 * t249 - t253 * t90) * mrSges(7,2)) - g(3) * ((-t154 * t249 + t175) * mrSges(7,1) + (-t154 * t253 - t364) * mrSges(7,2)) + (t167 * t20 + t168 * t21) * mrSges(7,3) + t328 + (-Ifges(7,2) * t168 + t163 + t69) * t406 + t417;];
tau  = t1;
