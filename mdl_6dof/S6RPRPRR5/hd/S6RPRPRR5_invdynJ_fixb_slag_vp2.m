% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:59
% EndTime: 2019-03-09 03:48:31
% DurationCPUTime: 20.51s
% Computational Cost: add. (10691->648), mult. (25383->824), div. (0->0), fcn. (19051->12), ass. (0->286)
t429 = -m(6) - m(7);
t428 = mrSges(4,1) + mrSges(5,1);
t427 = -mrSges(5,2) - mrSges(4,3);
t419 = Ifges(4,1) + Ifges(5,1);
t398 = Ifges(4,5) + Ifges(5,4);
t226 = sin(pkin(10));
t349 = pkin(7) + qJ(2);
t189 = t349 * t226;
t173 = qJD(1) * t189;
t227 = cos(pkin(10));
t190 = t349 * t227;
t174 = qJD(1) * t190;
t232 = sin(qJ(3));
t359 = cos(qJ(3));
t116 = -t359 * t173 - t232 * t174;
t172 = t226 * t359 + t232 * t227;
t161 = t172 * qJD(1);
t81 = pkin(8) * t161 + t116;
t426 = qJD(4) - t81;
t309 = t226 ^ 2 + t227 ^ 2;
t301 = qJD(1) * qJD(2);
t203 = qJ(2) * qJDD(1) + t301;
t418 = Ifges(5,5) - Ifges(4,4);
t230 = sin(qJ(6));
t234 = cos(qJ(6));
t191 = -t234 * mrSges(7,1) + t230 * mrSges(7,2);
t425 = -m(7) * pkin(5) - mrSges(6,1) + t191;
t222 = -qJD(3) + qJD(5);
t289 = t359 * t227;
t247 = -t232 * t226 + t289;
t160 = t247 * qJD(1);
t231 = sin(qJ(5));
t235 = cos(qJ(5));
t393 = -t160 * t231 + t235 * t161;
t345 = mrSges(6,3) * t393;
t82 = t222 * t234 - t230 * t393;
t83 = t222 * t230 + t234 * t393;
t347 = -mrSges(6,1) * t222 - mrSges(7,1) * t82 + mrSges(7,2) * t83 + t345;
t237 = -pkin(3) - pkin(4);
t74 = qJD(3) * t237 + t426;
t225 = qJD(3) * qJ(4);
t117 = -t232 * t173 + t359 * t174;
t275 = -pkin(8) * t160 + t117;
t75 = t225 + t275;
t47 = -t231 * t75 + t235 * t74;
t40 = -pkin(5) * t222 - t47;
t424 = -m(7) * t40 - t347;
t221 = pkin(10) + qJ(3);
t216 = cos(t221);
t233 = sin(qJ(1));
t315 = t216 * t233;
t195 = qJ(4) * t315;
t215 = sin(t221);
t296 = t215 * t237;
t423 = t233 * t296 + t195;
t362 = t172 / 0.2e1;
t211 = pkin(2) * t227 + pkin(1);
t421 = m(4) * t211;
t182 = -qJD(1) * t211 + qJD(2);
t86 = -t160 * pkin(3) - t161 * qJ(4) + t182;
t71 = pkin(4) * t160 - t86;
t420 = t71 * mrSges(6,2);
t417 = Ifges(4,6) - Ifges(5,6);
t218 = -qJDD(3) + qJDD(5);
t101 = -t160 * t235 - t231 * t161;
t162 = t247 * qJD(3);
t112 = qJD(1) * t162 + qJDD(1) * t172;
t163 = t172 * qJD(3);
t300 = qJDD(1) * t226;
t113 = qJD(1) * t163 - qJDD(1) * t289 + t232 * t300;
t49 = qJD(5) * t101 + t112 * t235 + t113 * t231;
t28 = qJD(6) * t82 + t218 * t230 + t234 * t49;
t29 = -qJD(6) * t83 + t218 * t234 - t230 * t49;
t12 = -mrSges(7,1) * t29 + mrSges(7,2) * t28;
t37 = mrSges(6,1) * t218 - mrSges(6,3) * t49;
t416 = t12 - t37;
t146 = Ifges(4,4) * t160;
t338 = Ifges(5,5) * t160;
t415 = qJD(3) * t398 + t419 * t161 + t146 - t338;
t92 = Ifges(6,4) * t101;
t414 = t101 * Ifges(6,2);
t413 = t222 * Ifges(6,5);
t412 = t222 * Ifges(6,6);
t184 = t235 * qJ(4) + t231 * t237;
t396 = qJD(5) * t184 + t231 * t426 + t235 * t275;
t395 = qJD(3) * t428 + t161 * t427;
t131 = mrSges(5,2) * t160 + qJD(3) * mrSges(5,3);
t394 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t160 + t131;
t181 = -qJDD(1) * t211 + qJDD(2);
t411 = -qJ(4) * t112 - qJD(4) * t161 + t181;
t410 = t428 * t216 + (-mrSges(4,2) + mrSges(5,3)) * t215;
t115 = t172 * t235 - t231 * t247;
t304 = qJD(6) * t234;
t252 = -t172 * t231 - t235 * t247;
t67 = qJD(5) * t252 + t162 * t235 + t163 * t231;
t250 = t115 * t304 + t230 * t67;
t93 = qJD(6) - t101;
t60 = Ifges(6,1) * t393 + t413 + t92;
t409 = t47 * mrSges(6,3) - t60 / 0.2e1;
t297 = m(7) * pkin(9) + mrSges(7,3);
t406 = -mrSges(6,2) + t297;
t30 = t83 * Ifges(7,5) + t82 * Ifges(7,6) + t93 * Ifges(7,3);
t48 = t231 * t74 + t235 * t75;
t341 = Ifges(6,4) * t393;
t59 = t341 + t412 + t414;
t405 = t48 * mrSges(6,3) - t30 / 0.2e1 + t59 / 0.2e1;
t271 = -mrSges(3,1) * t227 + mrSges(3,2) * t226;
t404 = -m(3) * pkin(1) - mrSges(2,1) + t271 - t410;
t33 = -pkin(5) * t101 - pkin(9) * t393 + t71;
t41 = pkin(9) * t222 + t48;
t15 = -t230 * t41 + t234 * t33;
t279 = pkin(7) * qJDD(1) + t203;
t140 = t279 * t226;
t141 = t279 * t227;
t284 = qJD(3) * t359;
t308 = qJD(3) * t232;
t66 = -t140 * t359 - t232 * t141 + t173 * t308 - t174 * t284;
t243 = qJDD(4) - t66;
t42 = -t112 * pkin(8) + qJDD(3) * t237 + t243;
t65 = -t232 * t140 + t359 * t141 - t173 * t284 - t174 * t308;
t61 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t65;
t43 = pkin(8) * t113 + t61;
t10 = qJD(5) * t47 + t231 * t42 + t235 * t43;
t7 = pkin(9) * t218 + t10;
t34 = t113 * t237 - t411;
t50 = -qJD(5) * t393 - t112 * t231 + t113 * t235;
t9 = -pkin(5) * t50 - pkin(9) * t49 + t34;
t1 = qJD(6) * t15 + t230 * t9 + t234 * t7;
t16 = t230 * t33 + t234 * t41;
t2 = -qJD(6) * t16 - t230 * t7 + t234 * t9;
t403 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t64 = pkin(5) * t393 - pkin(9) * t101;
t105 = t225 + t117;
t402 = t182 * mrSges(4,1) + t86 * mrSges(5,1) - t105 * mrSges(5,2) - t117 * mrSges(4,3);
t388 = qJD(4) - t116;
t104 = -qJD(3) * pkin(3) + t388;
t401 = t182 * mrSges(4,2) + t104 * mrSges(5,2) - t116 * mrSges(4,3) - t86 * mrSges(5,3);
t400 = t71 * mrSges(6,1) + t15 * mrSges(7,1) - t16 * mrSges(7,2);
t285 = m(3) * qJ(2) + mrSges(3,3);
t399 = mrSges(2,2) + mrSges(6,3) - t285 + t429 * (-pkin(8) + t349) + t427;
t264 = t230 * mrSges(7,1) + t234 * mrSges(7,2);
t248 = t40 * t264;
t55 = -mrSges(7,2) * t93 + mrSges(7,3) * t82;
t56 = mrSges(7,1) * t93 - mrSges(7,3) * t83;
t259 = t230 * t56 - t234 * t55;
t346 = mrSges(6,3) * t101;
t84 = -mrSges(6,2) * t222 + t346;
t251 = t259 - t84;
t121 = -t232 * t189 + t359 * t190;
t236 = cos(qJ(1));
t314 = t216 * t236;
t316 = t215 * t236;
t392 = pkin(3) * t314 + qJ(4) * t316;
t46 = qJDD(6) - t50;
t13 = mrSges(7,1) * t46 - mrSges(7,3) * t28;
t14 = -mrSges(7,2) * t46 + mrSges(7,3) * t29;
t390 = -t230 * t13 + t234 * t14;
t389 = g(1) * t236 + g(2) * t233;
t11 = -qJD(5) * t48 - t231 * t43 + t235 * t42;
t317 = t215 * t235;
t124 = t231 * t315 - t233 * t317;
t153 = t215 * t231 + t216 * t235;
t125 = t153 * t233;
t387 = -t125 * mrSges(6,2) + t124 * t425;
t311 = t235 * t236;
t313 = t231 * t236;
t126 = -t215 * t313 - t216 * t311;
t127 = -t215 * t311 + t216 * t313;
t386 = t126 * mrSges(6,2) + t127 * t425;
t154 = -t216 * t231 + t317;
t385 = -t154 * mrSges(6,2) + t153 * t425;
t241 = t1 * t234 - t2 * t230 + (-t15 * t234 - t16 * t230) * qJD(6);
t305 = qJD(6) * t230;
t384 = m(7) * t241 - t56 * t304 - t55 * t305 + t390;
t361 = -t222 / 0.2e1;
t371 = -t101 / 0.2e1;
t373 = -t93 / 0.2e1;
t375 = -t83 / 0.2e1;
t377 = -t82 / 0.2e1;
t383 = Ifges(7,5) * t375 - Ifges(6,2) * t371 - Ifges(6,6) * t361 + Ifges(7,6) * t377 + Ifges(7,3) * t373 - t400;
t261 = Ifges(7,5) * t234 - Ifges(7,6) * t230;
t339 = Ifges(7,4) * t234;
t262 = -Ifges(7,2) * t230 + t339;
t340 = Ifges(7,4) * t230;
t263 = Ifges(7,1) * t234 - t340;
t352 = t83 * Ifges(7,4);
t31 = t82 * Ifges(7,2) + t93 * Ifges(7,6) + t352;
t78 = Ifges(7,4) * t82;
t32 = t83 * Ifges(7,1) + t93 * Ifges(7,5) + t78;
t327 = t234 * t32;
t343 = mrSges(7,3) * t234;
t344 = mrSges(7,3) * t230;
t360 = t230 / 0.2e1;
t370 = -t393 / 0.2e1;
t382 = -t420 + Ifges(6,1) * t370 + Ifges(6,5) * t361 + t15 * t343 + t16 * t344 + t261 * t373 + t262 * t377 + t263 * t375 - t248 - t327 / 0.2e1 + t31 * t360;
t381 = t28 / 0.2e1;
t380 = t29 / 0.2e1;
t378 = t46 / 0.2e1;
t376 = t82 / 0.2e1;
t374 = t83 / 0.2e1;
t372 = t93 / 0.2e1;
t369 = t393 / 0.2e1;
t367 = t160 / 0.2e1;
t366 = -t160 / 0.2e1;
t364 = t161 / 0.2e1;
t208 = t216 * pkin(3);
t351 = -qJD(3) / 0.2e1;
t350 = qJD(3) / 0.2e1;
t342 = Ifges(4,4) * t161;
t322 = t115 * t230;
t321 = t115 * t234;
t206 = t215 * qJ(4);
t106 = t161 * pkin(3) - t160 * qJ(4);
t310 = t208 + t206;
t307 = qJD(3) * t235;
t303 = m(5) - t429;
t299 = qJDD(1) * t227;
t298 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t46;
t291 = t327 / 0.2e1;
t79 = t163 * pkin(3) - t162 * qJ(4) - t172 * qJD(4);
t290 = t216 * pkin(4) + t310;
t283 = -t305 / 0.2e1;
t282 = t113 * mrSges(4,1) + t112 * mrSges(4,2);
t281 = t113 * mrSges(5,1) - t112 * mrSges(5,3);
t278 = -t211 - t206;
t120 = t359 * t189 + t190 * t232;
t198 = t236 * t211;
t277 = t233 * t349 + t198;
t94 = -qJDD(3) * mrSges(5,1) + t112 * mrSges(5,2);
t276 = pkin(4) * t314 + t198 + t392;
t77 = -pkin(4) * t161 - t106;
t273 = -t50 * mrSges(6,1) + t49 * mrSges(6,2);
t272 = -mrSges(3,1) * t299 + mrSges(3,2) * t300;
t260 = t15 * t230 - t16 * t234;
t111 = -pkin(3) * t247 - t172 * qJ(4) - t211;
t80 = pkin(4) * t247 - t111;
t44 = -pkin(5) * t252 - pkin(9) * t115 + t80;
t90 = -pkin(8) * t172 + t120;
t91 = -pkin(8) * t247 + t121;
t58 = t231 * t90 + t235 * t91;
t22 = -t230 * t58 + t234 * t44;
t23 = t230 * t44 + t234 * t58;
t57 = t231 * t91 - t235 * t90;
t183 = -qJ(4) * t231 + t235 * t237;
t255 = -t125 * t234 - t230 * t236;
t254 = t125 * t230 - t234 * t236;
t69 = -pkin(4) * t163 - t79;
t249 = t115 * t305 - t234 * t67;
t244 = t216 * mrSges(5,3) + (-m(5) * pkin(3) - mrSges(5,1)) * t215;
t242 = (t216 * t237 + t278) * t233;
t87 = -t189 * t284 + qJD(2) * t289 + (-qJD(2) * t226 - qJD(3) * t190) * t232;
t88 = t172 * qJD(2) + qJD(3) * t121;
t240 = -t162 * pkin(8) + t88;
t5 = t28 * Ifges(7,4) + t29 * Ifges(7,2) + t46 * Ifges(7,6);
t6 = t28 * Ifges(7,1) + t29 * Ifges(7,4) + t46 * Ifges(7,5);
t8 = -pkin(5) * t218 - t11;
t239 = Ifges(6,5) * t49 + Ifges(6,6) * t50 + Ifges(6,3) * t218 + t11 * mrSges(6,1) + t8 * t191 + t1 * t343 - t10 * mrSges(6,2) + t6 * t360 + t234 * t5 / 0.2e1 + (Ifges(7,1) * t230 + t339) * t381 + (Ifges(7,2) * t234 + t340) * t380 + (Ifges(7,5) * t230 + Ifges(7,6) * t234) * t378 + t31 * t283 - t2 * t344 + (-t15 * t304 - t16 * t305) * mrSges(7,3) + (t261 * t372 + t262 * t376 + t263 * t374 + t248 + t291) * qJD(6);
t214 = -qJDD(1) * pkin(1) + qJDD(2);
t197 = qJ(4) * t314;
t179 = pkin(5) - t183;
t145 = Ifges(5,5) * t161;
t123 = t161 * t230 + t234 * t307;
t122 = t161 * t234 - t230 * t307;
t119 = -t126 * t234 - t230 * t233;
t118 = t126 * t230 - t233 * t234;
t107 = -mrSges(5,1) * t160 - mrSges(5,3) * t161;
t97 = t160 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t342;
t96 = Ifges(5,6) * qJD(3) - t160 * Ifges(5,3) + t145;
t95 = -mrSges(5,2) * t113 + qJDD(3) * mrSges(5,3);
t72 = pkin(8) * t163 + t87;
t68 = qJD(5) * t115 + t162 * t231 - t235 * t163;
t63 = -mrSges(6,1) * t101 + mrSges(6,2) * t393;
t62 = -qJDD(3) * pkin(3) + t243;
t54 = pkin(3) * t113 + t411;
t52 = t231 * t275 + t235 * t81;
t38 = -mrSges(6,2) * t218 + mrSges(6,3) * t50;
t35 = -t64 + t77;
t25 = t230 * t64 + t234 * t47;
t24 = -t230 * t47 + t234 * t64;
t21 = pkin(5) * t68 - pkin(9) * t67 + t69;
t19 = -qJD(5) * t57 + t231 * t240 + t235 * t72;
t18 = t230 * t35 + t234 * t52;
t17 = -t230 * t52 + t234 * t35;
t4 = -qJD(6) * t23 - t19 * t230 + t21 * t234;
t3 = qJD(6) * t22 + t19 * t234 + t21 * t230;
t20 = [0.2e1 * t309 * t203 * mrSges(3,3) + (t172 * t62 + t247 * t61) * mrSges(5,2) + t54 * (-mrSges(5,1) * t247 - mrSges(5,3) * t172) + t247 * (Ifges(4,4) * t112 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t247 * (Ifges(5,5) * t112 + Ifges(5,6) * qJDD(3)) / 0.2e1 + (t112 * t120 - t172 * t66 + t247 * t65) * mrSges(4,3) + (t34 * mrSges(6,2) - t11 * mrSges(6,3) + Ifges(6,1) * t49 + Ifges(6,4) * t50 + Ifges(6,5) * t218 + t261 * t378 + t262 * t380 + t263 * t381 + t264 * t8 + t283 * t32) * t115 + m(5) * (t111 * t54 + t79 * t86) + (-t412 / 0.2e1 - t414 / 0.2e1 - Ifges(6,4) * t369 + Ifges(7,3) * t372 + Ifges(7,5) * t374 + Ifges(7,6) * t376 + t400 - t405) * t68 + (-m(4) * t116 + m(5) * t104 - t395) * t88 - t250 * t31 / 0.2e1 + (m(4) * t117 + m(5) * t105 + t394) * t87 + (-m(4) * t66 + m(5) * t62 - qJDD(3) * mrSges(4,1) + t94) * t120 + (-Ifges(7,4) * t249 - Ifges(7,2) * t250) * t376 + (t415 / 0.2e1 + Ifges(4,4) * t367 + Ifges(5,5) * t366 + t398 * t350 + t419 * t364 + t401) * t162 - t5 * t322 / 0.2e1 + (m(4) * t65 + m(5) * t61 - qJDD(3) * mrSges(4,2) + t95) * t121 + (-t121 * mrSges(4,3) + (-Ifges(5,3) - Ifges(4,2)) * t247 + 0.2e1 * t418 * t362) * t113 + (-Ifges(7,5) * t249 - Ifges(7,6) * t250) * t372 + (-Ifges(7,1) * t249 - Ifges(7,4) * t250) * t374 + (-m(6) * t47 - t424) * (qJD(5) * t58 + t231 * t72 - t235 * t240) + t111 * t281 - t211 * t282 + (-m(6) * t11 + m(7) * t8 + t416) * t57 + (t96 / 0.2e1 - t97 / 0.2e1 + Ifges(5,3) * t366 - Ifges(4,2) * t367 + t418 * t364 - t417 * t350 + t402) * t163 + (t398 * t172 + t417 * t247) * qJDD(3) / 0.2e1 + (t419 * t172 - t418 * t247) * t112 / 0.2e1 + (t398 * qJDD(3) + t419 * t112) * t362 + m(7) * (t1 * t23 + t15 * t4 + t16 * t3 + t2 * t22) + m(6) * (t10 * t58 + t19 * t48 + t34 * t80 + t69 * t71) + t6 * t321 / 0.2e1 + (-t1 * t322 + t15 * t249 - t16 * t250 - t2 * t321) * mrSges(7,3) - pkin(1) * t272 + t80 * t273 + (t413 / 0.2e1 + t92 / 0.2e1 + t420 + t291 + Ifges(6,1) * t369 - t409) * t67 + (Ifges(3,4) * t226 + Ifges(3,2) * t227) * t299 + (Ifges(3,1) * t226 + Ifges(3,4) * t227) * t300 + t40 * (mrSges(7,1) * t250 - mrSges(7,2) * t249) + t214 * t271 + (-mrSges(4,1) * t247 + mrSges(4,2) * t172 - t421) * t181 + (-t255 * mrSges(7,1) - t254 * mrSges(7,2) - (-pkin(5) * t125 + t242) * m(7) - m(6) * t242 + t125 * mrSges(6,1) + t406 * t124 + (t421 - m(5) * (t278 - t208) - t404) * t233 + ((-m(4) - m(5)) * t349 + t399) * t236) * g(1) + Ifges(2,3) * qJDD(1) + t79 * t107 + t19 * t84 + t69 * t63 + t3 * t55 + t4 * t56 + t58 * t38 + t22 * t13 + t23 * t14 + (Ifges(6,6) * t218 + Ifges(6,4) * t49 + Ifges(6,2) * t50 - t34 * mrSges(6,1) - Ifges(7,3) * t378 - Ifges(7,6) * t380 - Ifges(7,5) * t381 - t298 / 0.2e1 + t10 * mrSges(6,3) + t403) * t252 + (-m(4) * t277 - m(6) * t276 + t126 * mrSges(6,1) - m(5) * (t277 + t392) - m(7) * (-pkin(5) * t126 + t276) - t119 * mrSges(7,1) - t118 * mrSges(7,2) - t406 * t127 + t404 * t236 + t399 * t233) * g(2) + m(3) * (-pkin(1) * t214 + (t203 + t301) * qJ(2) * t309); t281 + t282 - t273 + t259 * qJD(6) - t251 * t101 + t347 * t393 + m(3) * t214 + t272 + t395 * t161 - t234 * t13 - t230 * t14 - t394 * t160 + (-g(1) * t233 + g(2) * t236) * (m(3) + m(4) + t303) - t285 * t309 * qJD(1) ^ 2 + (-t1 * t230 - t2 * t234 + t260 * t93 + t393 * t40) * m(7) + (t101 * t48 - t393 * t47 - t34) * m(6) + (-t104 * t161 - t105 * t160 + t54) * m(5) + (t116 * t161 - t117 * t160 + t181) * m(4); t396 * t347 + (-m(7) * (-pkin(9) * t154 + t290) + t154 * mrSges(7,3) - m(6) * t290 - m(5) * t310 + t385 - t410) * g(3) + t384 * (-pkin(9) + t184) + t398 * t112 + (-Ifges(4,2) * t161 + t146 + t415) * t366 - t417 * t113 + (Ifges(5,3) * t367 - t351 * t417 - t402) * t161 - (-t351 * t398 + t401) * t160 + (-pkin(3) * t62 + qJ(4) * t61 - t104 * t117 + t105 * t388 - t106 * t86) * m(5) + (-m(5) * t197 - t236 * t244 - m(6) * (t236 * t296 + t197) - m(7) * (pkin(9) * t126 + t237 * t316 + t197) - t126 * mrSges(7,3) + t386) * g(1) - (t419 * t160 + t145 - t342 + t96) * t161 / 0.2e1 + t338 * t367 + (-m(5) * t195 - t233 * t244 - m(6) * t423 - m(7) * (-pkin(9) * t125 + t423) + t125 * mrSges(7,3) + t387) * g(2) + (m(6) * t48 - m(7) * t260 - t251) * (qJD(4) * t235 + qJD(5) * t183) + (-t15 * t17 - t16 * t18 + t179 * t8 + t396 * t40) * m(7) + (t10 * t184 + t11 * t183 - t396 * t47 - t48 * t52 - t71 * t77) * m(6) - t239 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) + t389 * (mrSges(4,1) * t215 + mrSges(4,2) * t216) - t394 * t116 + t395 * t117 + t179 * t12 + t183 * t37 + t184 * t38 + qJD(4) * t131 - t106 * t107 + qJ(4) * t95 - pkin(3) * t94 - t52 * t84 - t77 * t63 + t61 * mrSges(5,3) - t62 * mrSges(5,1) - t65 * mrSges(4,2) + t66 * mrSges(4,1) - t18 * t55 - t17 * t56 + t97 * t364 - (-Ifges(6,4) * t370 + t383 + t405) * t393 - (Ifges(6,4) * t371 + t382 + t409) * t101; -qJD(3) * t131 - t122 * t56 - t123 * t55 + (-t63 + t107) * t161 + (-qJD(3) * t84 - qJD(5) * t251 - t416) * t235 + (t38 + (-t230 * t55 - t234 * t56) * qJD(6) + t222 * t347 + t390) * t231 + t94 + (g(3) * t216 - t215 * t389) * t303 + (-t122 * t15 - t123 * t16 + (-qJD(5) * t260 - t8) * t235 + (t222 * t40 + t241) * t231) * m(7) + (t10 * t231 + t11 * t235 - t161 * t71 + t222 * (-t231 * t47 + t235 * t48)) * m(6) + (-qJD(3) * t105 + t161 * t86 + t62) * m(5); -pkin(5) * t12 - t24 * t56 - t25 * t55 + t59 * t369 + t239 + (t346 - t84) * t47 + (t92 + t60) * t371 + (-t341 + t30) * t370 + (-t154 * t297 - t385) * g(3) + (-t125 * t297 - t387) * g(2) + (t126 * t297 - t386) * g(1) + (-pkin(5) * t8 - t15 * t24 - t16 * t25) * m(7) + (t345 + t424) * t48 + t383 * t393 + t382 * t101 + t384 * pkin(9); -t40 * (mrSges(7,1) * t83 + mrSges(7,2) * t82) + (Ifges(7,1) * t82 - t352) * t375 + t31 * t374 + (Ifges(7,5) * t82 - Ifges(7,6) * t83) * t373 - t15 * t55 + t16 * t56 - g(1) * (mrSges(7,1) * t118 - mrSges(7,2) * t119) - g(2) * (-mrSges(7,1) * t254 + mrSges(7,2) * t255) + g(3) * t264 * t154 + (t15 * t82 + t16 * t83) * mrSges(7,3) + t298 + (-Ifges(7,2) * t83 + t32 + t78) * t377 - t403;];
tau  = t20;
