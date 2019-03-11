% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR4
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:28
% EndTime: 2019-03-09 03:44:57
% DurationCPUTime: 19.40s
% Computational Cost: add. (7094->696), mult. (13942->924), div. (0->0), fcn. (8373->14), ass. (0->323)
t232 = sin(qJ(3));
t214 = t232 * qJD(1);
t234 = cos(qJ(6));
t235 = cos(qJ(5));
t306 = t235 * t214;
t230 = sin(qJ(6));
t231 = sin(qJ(5));
t344 = t230 * t231;
t118 = -t214 * t344 + t234 * t306;
t317 = qJD(5) + qJD(6);
t323 = qJD(6) * t230;
t327 = qJD(5) * t231;
t339 = t234 * t235;
t89 = -t230 * t327 - t231 * t323 + t317 * t339;
t443 = t89 + t118;
t264 = t230 * t235 + t234 * t231;
t253 = t264 * t232;
t119 = qJD(1) * t253;
t90 = t317 * t264;
t474 = t90 + t119;
t455 = -Ifges(5,4) + Ifges(4,5);
t454 = Ifges(5,5) - Ifges(4,6);
t236 = cos(qJ(3));
t342 = t231 * t232;
t262 = pkin(5) * t236 - pkin(9) * t342;
t405 = pkin(3) + pkin(8);
t380 = pkin(9) + t405;
t228 = sin(pkin(10));
t197 = pkin(1) * t228 + pkin(7);
t175 = t197 * qJD(1);
t129 = t232 * qJD(2) + t236 * t175;
t333 = qJD(1) * t236;
t104 = pkin(4) * t333 + t129;
t206 = pkin(3) * t214;
t352 = qJ(4) * t236;
t270 = pkin(8) * t232 - t352;
t130 = qJD(1) * t270 + t206;
t59 = t235 * t104 - t130 * t231;
t473 = -qJD(1) * t262 + t380 * t327 - t59;
t171 = t380 * t235;
t60 = t231 * t104 + t235 * t130;
t472 = pkin(9) * t306 + qJD(5) * t171 + t60;
t227 = qJ(5) + qJ(6);
t216 = sin(t227);
t217 = cos(t227);
t283 = t231 * mrSges(6,1) + mrSges(6,2) * t235;
t393 = pkin(5) * t231;
t471 = -m(7) * t393 - t216 * mrSges(7,1) - t217 * mrSges(7,2) - t283;
t238 = -pkin(9) - pkin(8);
t470 = -m(7) * (-pkin(3) + t238) + mrSges(7,3) + m(6) * t405 + mrSges(6,3);
t193 = t214 + qJD(5);
t329 = qJD(3) * t235;
t254 = t231 * t333 - t329;
t215 = t236 * qJD(2);
t103 = -t232 * (pkin(4) * qJD(1) + t175) + t215;
t431 = qJD(4) - t103;
t88 = -qJD(3) * t405 + t431;
t218 = t232 * qJ(4);
t296 = -pkin(2) - t218;
t229 = cos(pkin(10));
t396 = pkin(1) * t229;
t97 = (-t236 * t405 + t296 - t396) * qJD(1);
t42 = -t231 * t97 + t235 * t88;
t38 = pkin(9) * t254 + t42;
t34 = pkin(5) * t193 + t38;
t331 = qJD(3) * t231;
t162 = -t235 * t333 - t331;
t43 = t231 * t88 + t235 * t97;
t39 = pkin(9) * t162 + t43;
t367 = t230 * t39;
t13 = t234 * t34 - t367;
t362 = t234 * t39;
t14 = t230 * t34 + t362;
t290 = t234 * t162 + t230 * t254;
t321 = qJD(1) * qJD(3);
t169 = qJDD(1) * t232 + t236 * t321;
t158 = qJDD(5) + t169;
t149 = qJDD(6) + t158;
t168 = -t236 * qJDD(1) + t232 * t321;
t78 = qJD(5) * t162 + qJDD(3) * t235 + t168 * t231;
t79 = qJD(5) * t254 - qJDD(3) * t231 + t168 * t235;
t23 = qJD(6) * t290 + t230 * t79 + t234 * t78;
t86 = t162 * t230 - t234 * t254;
t24 = -qJD(6) * t86 - t230 * t78 + t234 * t79;
t315 = Ifges(7,5) * t23 + Ifges(7,6) * t24 + Ifges(7,3) * t149;
t397 = Ifges(7,4) * t86;
t188 = t214 + t317;
t400 = -t188 / 0.2e1;
t407 = -t86 / 0.2e1;
t326 = qJD(5) * t235;
t328 = qJD(3) * t236;
t463 = qJD(2) * qJD(3) + t197 * qJDD(1);
t71 = qJDD(2) * t236 - t175 * t328 - t232 * t463;
t252 = qJDD(4) - t71;
t49 = pkin(4) * t169 - qJDD(3) * t405 + t252;
t198 = -pkin(2) - t396;
t173 = t198 * qJDD(1);
t244 = -qJ(4) * t169 - qJD(4) * t214 + t173;
t55 = t168 * t405 + t244;
t11 = t231 * t49 + t235 * t55 + t88 * t326 - t327 * t97;
t10 = pkin(9) * t79 + t11;
t12 = -qJD(5) * t43 - t231 * t55 + t235 * t49;
t8 = pkin(5) * t158 - pkin(9) * t78 + t12;
t3 = -qJD(6) * t14 - t10 * t230 + t234 * t8;
t467 = t3 * mrSges(7,1);
t2 = qJD(6) * t13 + t10 * t234 + t230 * t8;
t468 = t2 * mrSges(7,2);
t226 = qJD(3) * qJ(4);
t94 = t226 + t104;
t67 = -pkin(5) * t162 + t94;
t469 = t315 + t467 - t468 + (Ifges(7,5) * t290 - Ifges(7,6) * t86) * t400 + (t13 * t290 + t14 * t86) * mrSges(7,3) - t67 * (mrSges(7,1) * t86 + mrSges(7,2) * t290) + (Ifges(7,1) * t290 - t397) * t407;
t222 = qJ(1) + pkin(10);
t209 = sin(t222);
t466 = g(2) * t209;
t322 = m(5) + m(6) + m(7);
t465 = -m(4) - t322;
t114 = -t226 - t129;
t65 = -qJDD(3) * pkin(3) + t252;
t464 = -qJD(3) * t114 - t65;
t269 = t11 * t231 + t12 * t235;
t462 = -t43 * t326 + t42 * t327 - t269;
t282 = t236 * mrSges(5,2) - t232 * mrSges(5,3);
t286 = mrSges(4,1) * t236 - mrSges(4,2) * t232;
t461 = t282 - t286;
t263 = -t339 + t344;
t460 = t474 * t13 - t443 * t14 - t2 * t264 + t263 * t3;
t81 = Ifges(7,4) * t290;
t459 = -Ifges(7,2) * t86 + t81;
t378 = mrSges(6,3) * t162;
t105 = -mrSges(6,2) * t193 + t378;
t377 = mrSges(6,3) * t254;
t106 = mrSges(6,1) * t193 + t377;
t265 = t235 * t105 - t231 * t106;
t56 = mrSges(6,1) * t158 - mrSges(6,3) * t78;
t57 = -mrSges(6,2) * t158 + mrSges(6,3) * t79;
t458 = -t265 * qJD(5) - t231 * t57 - t235 * t56;
t371 = Ifges(5,6) * t236;
t272 = -t232 * Ifges(5,2) - t371;
t457 = t14 * mrSges(7,2) + Ifges(5,4) * qJD(3) / 0.2e1 + qJD(1) * t272 / 0.2e1 - t13 * mrSges(7,1);
t417 = t23 / 0.2e1;
t416 = t24 / 0.2e1;
t404 = t149 / 0.2e1;
t170 = t380 * t231;
t93 = -t170 * t234 - t171 * t230;
t453 = -qJD(6) * t93 + t230 * t472 + t234 * t473;
t92 = t170 * t230 - t171 * t234;
t452 = qJD(6) * t92 + t230 * t473 - t234 * t472;
t284 = mrSges(6,1) * t235 - mrSges(6,2) * t231;
t447 = t94 * t284;
t418 = m(7) * pkin(5);
t446 = -mrSges(6,1) - t418;
t312 = mrSges(5,1) * t333;
t179 = -qJD(3) * mrSges(5,3) - t312;
t91 = -mrSges(6,1) * t162 - mrSges(6,2) * t254;
t445 = -t179 + t91;
t391 = pkin(5) * t235;
t204 = pkin(4) + t391;
t444 = pkin(5) * t326 + qJD(4) - t215 - (-qJD(1) * t204 - t175) * t232;
t220 = t236 * pkin(3);
t334 = t220 + t218;
t150 = t198 - t334;
t390 = pkin(8) * t236;
t126 = t150 - t390;
t381 = pkin(4) + t197;
t153 = t381 * t232;
t133 = t231 * t153;
t69 = t235 * t126 + t133;
t137 = mrSges(5,1) * t168 - qJDD(3) * mrSges(5,3);
t441 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t168 - t137;
t311 = mrSges(4,3) * t214;
t313 = mrSges(5,1) * t214;
t440 = -t311 - t313 + (mrSges(4,1) - mrSges(5,2)) * qJD(3);
t310 = mrSges(4,3) * t333;
t178 = -qJD(3) * mrSges(4,2) + t310;
t439 = -t178 + t179;
t372 = Ifges(5,6) * t232;
t438 = t232 * (-Ifges(5,2) * t236 + t372) + t236 * (Ifges(5,3) * t232 - t371);
t437 = t454 * t232 + t236 * t455;
t330 = qJD(3) * t232;
t70 = t232 * qJDD(2) - t175 * t330 + t236 * t463;
t435 = -t232 * t71 + t236 * t70;
t62 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t70;
t433 = t232 * t65 - t236 * t62;
t210 = cos(t222);
t432 = g(1) * t210 + t466;
t205 = Ifges(4,4) * t333;
t430 = Ifges(4,1) * t214 + Ifges(4,5) * qJD(3) - Ifges(6,5) * t254 + t86 * Ifges(7,5) + t162 * Ifges(6,6) + Ifges(7,6) * t290 + t193 * Ifges(6,3) + t188 * Ifges(7,3) + t205;
t429 = -mrSges(3,1) + t461;
t271 = -t236 * Ifges(5,3) - t372;
t368 = t254 * Ifges(6,4);
t73 = t162 * Ifges(6,2) + t193 * Ifges(6,6) - t368;
t152 = Ifges(6,4) * t162;
t74 = -Ifges(6,1) * t254 + t193 * Ifges(6,5) + t152;
t428 = Ifges(5,5) * qJD(3) + qJD(1) * t271 + t231 * t74 + t235 * t73;
t359 = t236 * mrSges(5,3);
t427 = -t359 + t471 * t236 + (m(5) * pkin(3) - mrSges(5,2) + t470) * t232;
t426 = t236 * t317;
t138 = t169 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t424 = -t138 + t458;
t422 = -m(6) * pkin(4) - m(7) * t204 - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t420 = Ifges(7,4) * t417 + Ifges(7,2) * t416 + Ifges(7,6) * t404;
t419 = Ifges(7,1) * t417 + Ifges(7,4) * t416 + Ifges(7,5) * t404;
t415 = -t78 * Ifges(6,4) / 0.2e1 - t79 * Ifges(6,2) / 0.2e1 - t158 * Ifges(6,6) / 0.2e1;
t36 = Ifges(7,2) * t290 + Ifges(7,6) * t188 + t397;
t414 = -t36 / 0.2e1;
t37 = Ifges(7,1) * t86 + Ifges(7,5) * t188 + t81;
t413 = -t37 / 0.2e1;
t412 = t37 / 0.2e1;
t411 = t78 / 0.2e1;
t410 = t79 / 0.2e1;
t409 = -t290 / 0.2e1;
t408 = t290 / 0.2e1;
t406 = t86 / 0.2e1;
t403 = t158 / 0.2e1;
t401 = -t254 / 0.2e1;
t399 = t188 / 0.2e1;
t233 = sin(qJ(1));
t395 = pkin(1) * t233;
t394 = pkin(5) * t254;
t237 = cos(qJ(1));
t221 = t237 * pkin(1);
t379 = mrSges(7,1) * t217;
t376 = Ifges(4,4) * t232;
t375 = Ifges(4,4) * t236;
t374 = Ifges(6,4) * t231;
t373 = Ifges(6,4) * t235;
t358 = t236 * mrSges(7,3);
t345 = t217 * t232;
t107 = -t209 * t216 + t210 * t345;
t346 = t216 * t232;
t108 = t209 * t217 + t210 * t346;
t354 = t107 * mrSges(7,1) - t108 * mrSges(7,2);
t351 = qJD(3) * t67;
t350 = qJD(3) * t94;
t347 = t210 * t236;
t341 = t231 * t236;
t340 = t232 * t235;
t337 = t235 * t236;
t185 = t236 * t197;
t336 = t236 * t238;
t109 = t209 * t345 + t210 * t216;
t110 = -t209 * t346 + t210 * t217;
t335 = t109 * mrSges(7,1) + t110 * mrSges(7,2);
t154 = t236 * pkin(4) + t185;
t325 = qJD(5) * t236;
t41 = -mrSges(7,1) * t290 + mrSges(7,2) * t86;
t316 = -t41 - t445;
t314 = Ifges(6,5) * t78 + Ifges(6,6) * t79 + Ifges(6,3) * t158;
t307 = t210 * pkin(2) + t209 * pkin(7) + t221;
t303 = t231 * t325;
t299 = -t326 / 0.2e1;
t199 = qJ(4) + t393;
t295 = pkin(9) * t236 - t126;
t294 = -t321 / 0.2e1;
t289 = pkin(3) * t330 - qJD(4) * t232;
t112 = qJD(3) * t270 + t289;
t140 = t381 * t328;
t291 = -t112 * t231 + t235 * t140;
t128 = t175 * t232 - t215;
t285 = mrSges(4,1) * t232 + mrSges(4,2) * t236;
t281 = Ifges(6,1) * t235 - t374;
t280 = Ifges(6,1) * t231 + t373;
t279 = t236 * Ifges(4,2) + t376;
t277 = -Ifges(6,2) * t231 + t373;
t276 = Ifges(6,2) * t235 + t374;
t274 = Ifges(6,5) * t235 - Ifges(6,6) * t231;
t273 = Ifges(6,5) * t231 + Ifges(6,6) * t235;
t134 = t235 * t153;
t58 = pkin(5) * t232 + t231 * t295 + t134;
t61 = -pkin(9) * t337 + t69;
t26 = -t230 * t61 + t234 * t58;
t27 = t230 * t58 + t234 * t61;
t268 = t42 * t231 - t43 * t235;
t261 = t296 - t220;
t122 = -t209 * t231 + t210 * t340;
t124 = t209 * t340 + t210 * t231;
t115 = (t261 - t396) * qJD(1);
t259 = t115 * (-mrSges(5,2) * t232 - t359);
t258 = t198 * qJD(1) * t285;
t257 = t232 * (Ifges(4,1) * t236 - t376);
t32 = t235 * t112 - t126 * t327 + t231 * t140 + t153 * t326;
t250 = t232 * t329 + t303;
t249 = t231 * t330 - t235 * t325;
t247 = Ifges(6,5) * t236 + t232 * t280;
t246 = Ifges(6,6) * t236 + t232 * t276;
t245 = Ifges(6,3) * t236 + t232 * t273;
t52 = -pkin(4) * t168 - t62;
t243 = -qJD(5) * t268 + t269;
t187 = t236 * t216 * mrSges(7,2);
t167 = -qJ(4) * t333 + t206;
t166 = t282 * qJD(1);
t151 = t284 * t236;
t144 = Ifges(4,6) * qJD(3) + qJD(1) * t279;
t141 = -qJ(4) * t328 + t289;
t139 = t381 * t330;
t136 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t169;
t132 = t264 * t236;
t131 = t263 * t236;
t125 = -t209 * t342 + t210 * t235;
t123 = t209 * t235 + t210 * t342;
t116 = pkin(5) * t337 + t154;
t111 = -qJD(3) * pkin(3) + qJD(4) + t128;
t82 = -pkin(5) * t303 + (-t197 - t204) * t330;
t68 = -t126 * t231 + t134;
t66 = pkin(3) * t168 + t244;
t64 = mrSges(7,1) * t188 - mrSges(7,3) * t86;
t63 = -mrSges(7,2) * t188 + mrSges(7,3) * t290;
t51 = -t263 * t330 + t264 * t426;
t50 = qJD(3) * t253 + t263 * t426;
t40 = -mrSges(6,1) * t79 + mrSges(6,2) * t78;
t33 = -qJD(5) * t69 + t291;
t31 = t78 * Ifges(6,1) + t79 * Ifges(6,4) + t158 * Ifges(6,5);
t29 = -pkin(5) * t79 + t52;
t28 = pkin(9) * t250 + t32;
t25 = t262 * qJD(3) + (t235 * t295 - t133) * qJD(5) + t291;
t18 = -mrSges(7,2) * t149 + mrSges(7,3) * t24;
t17 = mrSges(7,1) * t149 - mrSges(7,3) * t23;
t16 = t234 * t38 - t367;
t15 = -t230 * t38 - t362;
t9 = -mrSges(7,1) * t24 + mrSges(7,2) * t23;
t5 = -qJD(6) * t27 - t230 * t28 + t234 * t25;
t4 = qJD(6) * t26 + t230 * t25 + t234 * t28;
t1 = [(Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t229 - 0.2e1 * mrSges(3,2) * t228 + m(3) * (t228 ^ 2 + t229 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (m(4) * t198 - t286) * t173 + t232 * t467 + t437 * qJD(3) ^ 2 / 0.2e1 + t438 * t294 + t441 * t185 + (-m(3) * t221 - m(4) * t307 - mrSges(2,1) * t237 - t123 * mrSges(6,1) - t108 * mrSges(7,1) + mrSges(2,2) * t233 - t122 * mrSges(6,2) - t107 * mrSges(7,2) + (-m(6) * pkin(8) - mrSges(6,3)) * t347 - t322 * (pkin(3) * t347 + t210 * t218 + t307) + t422 * t209 + (-m(7) * (pkin(5) * t342 - t336) - t358 + t429) * t210) * g(2) + (t111 * mrSges(5,1) + t128 * mrSges(4,3) + t430 / 0.2e1 + Ifges(7,3) * t399 - t43 * mrSges(6,2) + t42 * mrSges(6,1) + Ifges(7,5) * t406 + Ifges(7,6) * t408 - t457) * t328 + m(7) * (t116 * t29 + t13 * t5 + t14 * t4 + t2 * t27 + t26 * t3 + t67 * t82) + t168 * t271 / 0.2e1 - t169 * t272 / 0.2e1 + t94 * (-mrSges(6,1) * t250 + mrSges(6,2) * t249) + (Ifges(7,5) * t50 + Ifges(7,6) * t51) * t399 + (t114 * mrSges(5,1) - t129 * mrSges(4,3) + t428 / 0.2e1 - t144 / 0.2e1) * t330 + (Ifges(7,4) * t50 + Ifges(7,2) * t51) * t408 - t168 * t279 / 0.2e1 + t66 * t282 + (t236 * (-Ifges(4,2) * t232 + t375) + t257) * t321 / 0.2e1 + (Ifges(7,1) * t50 + Ifges(7,4) * t51) * t406 + m(5) * (t115 * t141 + t150 * t66) + (m(5) * ((t111 * t236 + t114 * t232) * qJD(3) + t433) + m(4) * ((t128 * t236 - t129 * t232) * qJD(3) + t435) + (t138 - t136) * t232 + t439 * t330 - t440 * t328) * t197 - t42 * mrSges(6,3) * t249 + (-Ifges(7,5) * t132 + Ifges(7,6) * t131 + Ifges(7,3) * t232) * t404 + (-Ifges(7,4) * t132 + Ifges(7,2) * t131 + Ifges(7,6) * t232) * t416 + (-Ifges(7,1) * t132 + Ifges(7,4) * t131 + Ifges(7,5) * t232) * t417 + t29 * (-mrSges(7,1) * t131 - mrSges(7,2) * t132) + m(6) * (t11 * t69 + t12 * t68 - t139 * t94 + t154 * t52 + t32 * t43 + t33 * t42) + qJD(3) * t258 + qJD(3) * t259 + t435 * mrSges(4,3) + t433 * mrSges(5,1) + t43 * mrSges(6,3) * t250 + t12 * (mrSges(6,1) * t232 + mrSges(6,3) * t341) - t31 * t341 / 0.2e1 + t11 * (-mrSges(6,2) * t232 - mrSges(6,3) * t337) + t193 * (qJD(3) * t245 - t274 * t325) / 0.2e1 + t162 * (qJD(3) * t246 - t277 * t325) / 0.2e1 + (-Ifges(4,4) * t168 + Ifges(4,5) * qJDD(3) + t314 + t315) * t232 / 0.2e1 + t169 * t375 / 0.2e1 + t236 * (Ifges(4,4) * t169 - Ifges(4,2) * t168 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t236 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t169 + Ifges(5,3) * t168) / 0.2e1 + (-t13 * t50 + t131 * t2 + t132 * t3 + t14 * t51) * mrSges(7,3) - t232 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t169 + Ifges(5,6) * t168) / 0.2e1 + t73 * t303 / 0.2e1 + t198 * (mrSges(4,1) * t168 + mrSges(4,2) * t169) + t150 * (-mrSges(5,2) * t168 - mrSges(5,3) * t169) + t141 * t166 + t154 * t40 - t232 * t468 + (t232 * t455 - t236 * t454) * qJDD(3) / 0.2e1 + t169 * t232 * Ifges(4,1) + t236 * t74 * t299 + (qJD(3) * t247 - t281 * t325) * t401 + (Ifges(6,3) * t232 - t236 * t273) * t403 + (Ifges(6,6) * t232 - t236 * t276) * t410 + (Ifges(6,5) * t232 - t236 * t280) * t411 + t50 * t412 + t337 * t415 - t132 * t419 + t131 * t420 + (m(3) * t395 + mrSges(2,1) * t233 - t125 * mrSges(6,1) - t110 * mrSges(7,1) + mrSges(2,2) * t237 + t124 * mrSges(6,2) + t109 * mrSges(7,2) + t465 * (t210 * pkin(7) - t395) + t422 * t210 + (-m(7) * (-t199 * t232 - pkin(2)) - m(6) * t296 - m(5) * t261 + m(4) * pkin(2) + t470 * t236 - t429) * t209) * g(1) + t26 * t17 + t27 * t18 + t51 * t36 / 0.2e1 + t4 * t63 + t5 * t64 + t67 * (-mrSges(7,1) * t51 + mrSges(7,2) * t50) + t68 * t56 + t69 * t57 + t82 * t41 + t32 * t105 + t33 * t106 + t116 * t9 - t139 * t91 + t52 * t151; m(3) * qJDD(2) + m(7) * (t13 * t51 + t131 * t3 - t132 * t2 + t14 * t50) + t50 * t63 + t51 * t64 + t131 * t17 - t132 * t18 + (-m(3) + t465) * g(3) + (t40 + t9 + (t105 * t231 + t106 * t235 - t440) * qJD(3) + m(5) * (qJD(3) * t111 - t62) + m(4) * (qJD(3) * t128 + t70) + m(7) * t29 + m(6) * (t329 * t42 + t331 * t43 + t52) + t441) * t232 + (t136 + (t178 - t316) * qJD(3) + m(5) * t464 + m(4) * (qJD(3) * t129 + t71) + m(7) * t351 + m(6) * (t350 + t462) + t424) * t236; (mrSges(7,1) * t443 - mrSges(7,2) * t474) * t67 + (Ifges(7,5) * t119 + Ifges(7,6) * t118) * t400 + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + (-t258 - t259 - t43 * (-mrSges(6,2) * t236 + mrSges(6,3) * t340) - t42 * (mrSges(6,1) * t236 - mrSges(6,3) * t342)) * qJD(1) + (t40 - t137) * qJ(4) + (t438 / 0.2e1 - t257 / 0.2e1) * qJD(1) ^ 2 - (t162 * t276 + t193 * t273 - t254 * t280) * qJD(5) / 0.2e1 - (t162 * t246 + t193 * t245 - t247 * t254) * qJD(1) / 0.2e1 - t264 * t420 + t29 * (mrSges(7,1) * t264 - mrSges(7,2) * t263) + (-Ifges(7,5) * t263 - Ifges(7,6) * t264) * t404 + (-Ifges(7,4) * t263 - Ifges(7,2) * t264) * t416 + (-Ifges(7,1) * t263 - Ifges(7,4) * t264) * t417 - t263 * t419 + (-Ifges(7,5) * t90 - Ifges(7,6) * t89) * t399 + (-Ifges(7,1) * t90 - Ifges(7,4) * t89) * t406 + (-Ifges(7,4) * t90 - Ifges(7,2) * t89) * t408 + (-t322 * t352 + t427) * t466 + (Ifges(7,4) * t119 + Ifges(7,2) * t118) * t409 + (Ifges(7,1) * t119 + Ifges(7,4) * t118) * t407 + (-pkin(3) * t65 - qJ(4) * t62 - qJD(4) * t114 - t115 * t167) * m(5) + t443 * t414 + t444 * t41 + t445 * qJD(4) + (t447 + t144 / 0.2e1) * t214 + t437 * t294 + (-m(5) * t114 - t310 - t439) * t128 + (-m(5) * t111 + t311 + t440) * t129 + t432 * t285 - (-Ifges(4,2) * t214 + t205 + t430) * t333 / 0.2e1 + (-qJ(4) * t322 * t347 + t210 * t427) * g(1) - t428 * t214 / 0.2e1 + (qJ(4) * t52 - t42 * t59 - t43 * t60 + t431 * t94) * m(6) + (-m(6) * t243 + t458) * t405 + t460 * mrSges(7,3) + t462 * mrSges(6,3) + (Ifges(7,5) * t407 + Ifges(7,6) * t409 + Ifges(7,3) * t400 + t457) * t333 + t73 * t299 + qJD(5) * t447 + t52 * t283 - t74 * t327 / 0.2e1 - t111 * t312 - t114 * t313 + t235 * t31 / 0.2e1 + t199 * t9 - t167 * t166 + t452 * t63 + t453 * t64 + (t13 * t453 + t14 * t452 + t199 * t29 + t2 * t93 + t3 * t92 + t444 * t67) * m(7) + t454 * t168 + t455 * t169 + t274 * t403 + t277 * t410 + t281 * t411 - t90 * t412 + t119 * t413 + t231 * t415 + (-m(6) * (t334 + t390) - t236 * mrSges(6,3) - m(7) * (t334 - t336) - t358 - m(5) * t334 + t471 * t232 + t461) * g(3) - t62 * mrSges(5,3) + t65 * mrSges(5,2) - t70 * mrSges(4,2) + t71 * mrSges(4,1) + t92 * t17 + t93 * t18 - t103 * t91 - t60 * t105 - t59 * t106 - pkin(3) * t138; t264 * t18 - t263 * t17 - t474 * t64 + t443 * t63 + t316 * qJD(3) + t322 * t236 * g(3) + ((t166 + t265) * qJD(1) - t432 * t322) * t232 - t424 + (-t351 - t460) * m(7) + (-t214 * t268 + t243 - t350) * m(6) + (t115 * t214 - t464) * m(5); t469 - t86 * t414 + t254 * (Ifges(6,1) * t162 + t368) / 0.2e1 - (Ifges(6,2) * t254 + t152 + t74) * t162 / 0.2e1 - t193 * (Ifges(6,5) * t162 + Ifges(6,6) * t254) / 0.2e1 - t94 * (-mrSges(6,1) * t254 + mrSges(6,2) * t162) + (mrSges(6,2) * t123 + t122 * t446 - t354) * g(1) + (-mrSges(6,2) * t125 + t124 * t446 - t335) * g(2) + t459 * t409 + t314 + (-t377 + t106) * t43 - m(7) * (t13 * t15 + t14 * t16 - t394 * t67) + (t378 - t105) * t42 + (-t187 - (-m(7) * t391 - t379) * t236 + t151) * g(3) + t290 * t413 + t41 * t394 + (t18 * t230 - t323 * t64 + (qJD(6) * t63 + t17) * t234) * pkin(5) - t11 * mrSges(6,2) + t12 * mrSges(6,1) + t73 * t401 + (t2 * t230 + t234 * t3 + (-t13 * t230 + t14 * t234) * qJD(6)) * t418 - t16 * t63 - t15 * t64; t36 * t406 - t13 * t63 + t14 * t64 - g(1) * t354 - g(2) * t335 - g(3) * (-t236 * t379 + t187) + (t37 + t459) * t409 + t469;];
tau  = t1;
