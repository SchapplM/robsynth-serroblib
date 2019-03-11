% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR5
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:12
% EndTime: 2019-03-08 21:18:43
% DurationCPUTime: 21.60s
% Computational Cost: add. (5500->683), mult. (12226->926), div. (0->0), fcn. (8696->14), ass. (0->313)
t464 = mrSges(4,1) - mrSges(5,2);
t463 = qJD(3) / 0.2e1;
t243 = sin(qJ(3));
t339 = qJD(3) * t243;
t224 = pkin(3) * t339;
t246 = cos(qJ(3));
t337 = qJD(4) * t243;
t273 = -qJD(5) * t246 - t337;
t276 = -qJ(4) * t246 + qJ(5) * t243;
t110 = qJD(3) * t276 + t224 + t273;
t338 = qJD(3) * t246;
t396 = pkin(4) + pkin(8);
t182 = t396 * t338;
t236 = sin(pkin(11));
t238 = cos(pkin(11));
t244 = sin(qJ(2));
t247 = cos(qJ(2));
t349 = t243 * t247;
t271 = -t236 * t244 + t238 * t349;
t237 = sin(pkin(6));
t346 = qJD(1) * t237;
t438 = -t110 * t236 + t238 * t182 - t271 * t346;
t270 = t236 * t349 + t238 * t244;
t437 = t238 * t110 + t236 * t182 - t270 * t346;
t462 = -m(5) - m(6);
t335 = qJD(2) * qJD(3);
t183 = -t246 * qJDD(2) + t243 * t335;
t128 = -qJDD(3) * t236 + t183 * t238;
t184 = qJDD(2) * t243 + t246 * t335;
t319 = t244 * t346;
t187 = qJD(2) * pkin(8) + t319;
t239 = cos(pkin(6));
t333 = qJDD(1) * t239;
t344 = qJD(2) * t237;
t313 = qJD(1) * t344;
t205 = t247 * t313;
t334 = qJDD(1) * t237;
t144 = t244 * t334 + t205;
t340 = qJD(3) * t239;
t454 = qJDD(2) * pkin(8) + qJD(1) * t340 + t144;
t52 = -t187 * t338 - t243 * t454 + t246 * t333;
t256 = qJDD(4) - t52;
t381 = pkin(3) + qJ(5);
t32 = pkin(4) * t184 - qJD(3) * qJD(5) - qJDD(3) * t381 + t256;
t204 = t244 * t313;
t143 = t247 * t334 - t204;
t130 = -qJDD(2) * pkin(2) - t143;
t253 = -qJ(4) * t184 + t130;
t38 = qJD(2) * t273 + t183 * t381 + t253;
t13 = t236 * t32 + t238 * t38;
t10 = pkin(9) * t128 + t13;
t242 = sin(qJ(6));
t245 = cos(qJ(6));
t12 = -t236 * t38 + t238 * t32;
t129 = qJDD(3) * t238 + t183 * t236;
t5 = pkin(5) * t184 - pkin(9) * t129 + t12;
t336 = t238 * qJD(3);
t341 = qJD(2) * t246;
t171 = t236 * t341 - t336;
t343 = qJD(2) * t243;
t361 = qJ(4) * t243;
t307 = -pkin(2) - t361;
t169 = -t246 * t381 + t307;
t318 = t247 * t346;
t107 = qJD(2) * t169 - t318;
t345 = qJD(1) * t239;
t211 = t246 * t345;
t108 = -t243 * (pkin(4) * qJD(2) + t187) + t211;
t81 = -qJD(3) * t381 + qJD(4) - t108;
t39 = -t107 * t236 + t238 * t81;
t21 = pkin(5) * t343 + pkin(9) * t171 + t39;
t172 = -qJD(3) * t236 - t238 * t341;
t40 = t238 * t107 + t236 * t81;
t24 = pkin(9) * t172 + t40;
t8 = t21 * t245 - t24 * t242;
t1 = qJD(6) * t8 + t10 * t245 + t242 * t5;
t461 = t1 * mrSges(7,2);
t9 = t21 * t242 + t24 * t245;
t2 = -qJD(6) * t9 - t10 * t242 + t245 * t5;
t460 = t2 * mrSges(7,1);
t444 = Ifges(4,1) + Ifges(6,3);
t443 = -Ifges(5,4) + Ifges(4,5);
t442 = Ifges(5,5) - Ifges(4,6);
t356 = t236 * t243;
t272 = pkin(5) * t246 - pkin(9) * t356;
t459 = qJD(3) * t272 + t438;
t458 = -pkin(9) * t243 * t336 - t437;
t315 = t238 * t343;
t357 = t236 * t242;
t126 = t245 * t315 - t343 * t357;
t425 = t238 * t245 - t357;
t157 = t425 * qJD(6);
t432 = t157 + t126;
t275 = t245 * t236 + t242 * t238;
t264 = t275 * t243;
t127 = qJD(2) * t264;
t158 = t275 * qJD(6);
t431 = -t158 - t127;
t324 = mrSges(4,3) * t343;
t326 = mrSges(5,1) * t343;
t457 = qJD(3) * t464 - t324 - t326;
t373 = Ifges(5,6) * t246;
t279 = -t243 * Ifges(5,2) - t373;
t456 = t9 * mrSges(7,2) + Ifges(5,4) * t463 + qJD(2) * t279 / 0.2e1 - t8 * mrSges(7,1);
t421 = -m(7) + t462;
t455 = -m(4) + t421;
t286 = mrSges(5,2) * t246 - mrSges(5,3) * t243;
t231 = pkin(11) + qJ(6);
t226 = sin(t231);
t227 = cos(t231);
t287 = t226 * mrSges(7,1) + t227 * mrSges(7,2);
t291 = mrSges(4,1) * t246 - mrSges(4,2) * t243;
t387 = pkin(5) * t236;
t445 = m(7) * (-pkin(9) - qJ(5));
t453 = -t243 * (m(7) * t387 + t287) + t246 * (-mrSges(7,3) + t445) + t286 - mrSges(3,1) - t291;
t362 = sin(pkin(10));
t299 = t362 * t247;
t363 = cos(pkin(10));
t302 = t363 * t244;
t154 = t239 * t302 + t299;
t304 = t237 * t363;
t101 = t154 * t243 + t246 * t304;
t300 = t362 * t244;
t301 = t363 * t247;
t156 = -t239 * t300 + t301;
t303 = t237 * t362;
t103 = t156 * t243 - t246 * t303;
t354 = t237 * t244;
t322 = t243 * t354;
t160 = -t239 * t246 + t322;
t263 = -g(1) * t103 - g(2) * t101 - g(3) * t160;
t452 = -t1 * t275 - t2 * t425 - t431 * t8 - t432 * t9 - t263;
t417 = -m(6) * qJ(5) - mrSges(6,3);
t451 = -pkin(3) * t421 - t417 - t445 + t464;
t450 = -t227 * mrSges(7,1) + t226 * mrSges(7,2) - mrSges(5,1) + mrSges(3,2) - mrSges(4,3);
t295 = t171 * t242 + t245 * t172;
t30 = qJD(6) * t295 + t128 * t242 + t129 * t245;
t405 = t30 / 0.2e1;
t83 = t171 * t245 - t172 * t242;
t31 = qJD(6) * t83 + t128 * t245 - t129 * t242;
t404 = t31 / 0.2e1;
t173 = qJDD(6) + t184;
t393 = t173 / 0.2e1;
t449 = -t184 / 0.2e1;
t202 = t396 * t243;
t177 = t238 * t202;
t69 = pkin(5) * t243 + t177 + (pkin(9) * t246 - t169) * t236;
t350 = t238 * t246;
t92 = t238 * t169 + t236 * t202;
t74 = -pkin(9) * t350 + t92;
t23 = t242 * t69 + t245 * t74;
t448 = -qJD(6) * t23 + t242 * t458 + t245 * t459;
t22 = -t242 * t74 + t245 * t69;
t447 = qJD(6) * t22 + t242 * t459 - t245 * t458;
t120 = t246 * t187 + t243 * t345;
t235 = qJD(3) * qJ(4);
t114 = -t235 - t120;
t446 = m(5) * t114;
t11 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t64 = -t128 * mrSges(6,1) + t129 * mrSges(6,2);
t441 = t11 + t64;
t109 = pkin(4) * t341 + t120;
t222 = pkin(3) * t343;
t140 = qJD(2) * t276 + t222;
t56 = t238 * t109 - t140 * t236;
t42 = qJD(2) * t272 + t56;
t57 = t236 * t109 + t238 * t140;
t46 = pkin(9) * t315 + t57;
t379 = -pkin(9) - t381;
t185 = t379 * t236;
t186 = t379 * t238;
t97 = -t185 * t242 + t186 * t245;
t440 = -qJD(5) * t275 + qJD(6) * t97 - t242 * t42 - t245 * t46;
t98 = t185 * t245 + t186 * t242;
t439 = -qJD(5) * t425 - qJD(6) * t98 + t242 * t46 - t245 * t42;
t153 = -t239 * t301 + t300;
t360 = t153 * t246;
t436 = -pkin(3) * t360 - t153 * t361;
t155 = t239 * t299 + t302;
t359 = t155 * t246;
t435 = -pkin(3) * t359 - t155 * t361;
t147 = mrSges(5,1) * t183 - qJDD(3) * mrSges(5,3);
t434 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t183 - t147;
t148 = t184 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t433 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t184 + t148;
t323 = mrSges(4,3) * t341;
t195 = -qJD(3) * mrSges(4,2) + t323;
t325 = mrSges(5,1) * t341;
t196 = -qJD(3) * mrSges(5,3) - t325;
t429 = t195 - t196;
t374 = Ifges(5,6) * t243;
t428 = t243 * (-Ifges(5,2) * t246 + t374) + t246 * (Ifges(5,3) * t243 - t373);
t427 = t243 * t442 + t246 * t443;
t51 = -t187 * t339 + t243 * t333 + t246 * t454;
t424 = -t243 * t52 + t246 * t51;
t45 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t51;
t47 = -qJDD(3) * pkin(3) + t256;
t423 = t243 * t47 - t246 * t45;
t76 = -mrSges(6,2) * t184 + mrSges(6,3) * t128;
t77 = mrSges(6,1) * t184 - mrSges(6,3) * t129;
t422 = t236 * t76 + t238 * t77;
t277 = t12 * t238 + t13 * t236;
t215 = qJD(6) + t343;
t221 = Ifges(4,4) * t341;
t419 = Ifges(4,5) * qJD(3) - t171 * Ifges(6,5) - Ifges(7,5) * t83 + t172 * Ifges(6,6) + Ifges(7,6) * t295 + t215 * Ifges(7,3) + t343 * t444 + t221;
t41 = -mrSges(7,1) * t295 - mrSges(7,2) * t83;
t94 = -mrSges(6,1) * t172 - mrSges(6,2) * t171;
t330 = t196 - t94 - t41;
t290 = t238 * mrSges(6,1) - t236 * mrSges(6,2);
t278 = -t246 * Ifges(5,3) - t374;
t415 = t236 * (-t171 * Ifges(6,1) + t172 * Ifges(6,4) + Ifges(6,5) * t343) + t238 * (-t171 * Ifges(6,4) + t172 * Ifges(6,2) + Ifges(6,6) * t343) + Ifges(5,5) * qJD(3) + qJD(2) * t278;
t216 = qJ(4) + t387;
t289 = mrSges(6,1) * t236 + mrSges(6,2) * t238;
t414 = -m(7) * t216 + qJ(4) * t462 + mrSges(4,2) - mrSges(5,3) - t287 - t289;
t87 = qJD(5) + t235 + t109;
t61 = -pkin(5) * t172 + t87;
t413 = -m(6) * t87 - m(7) * t61 + t330;
t189 = -pkin(3) * t246 + t307;
t122 = qJD(2) * t189 - t318;
t188 = -qJD(2) * pkin(2) - t318;
t376 = Ifges(6,4) * t236;
t282 = Ifges(6,2) * t238 + t376;
t375 = Ifges(6,4) * t238;
t285 = Ifges(6,1) * t236 + t375;
t352 = t238 * t243;
t411 = -t122 * (-mrSges(5,2) * t243 - mrSges(5,3) * t246) - t188 * (mrSges(4,1) * t243 + mrSges(4,2) * t246) - t39 * (mrSges(6,1) * t246 - mrSges(6,3) * t356) - t40 * (-mrSges(6,2) * t246 + mrSges(6,3) * t352) + t243 * t87 * t290 - t172 * (Ifges(6,6) * t246 + t243 * t282) / 0.2e1 + t171 * (Ifges(6,5) * t246 + t243 * t285) / 0.2e1;
t409 = mrSges(6,1) * t356 + mrSges(6,2) * t352 - t453;
t220 = pkin(5) * t238 + pkin(4);
t380 = pkin(8) + t220;
t408 = -m(6) * t396 - m(7) * t380 - t290 + t450;
t248 = qJD(2) ^ 2;
t407 = Ifges(7,4) * t405 + Ifges(7,2) * t404 + Ifges(7,6) * t393;
t406 = Ifges(7,1) * t405 + Ifges(7,4) * t404 + Ifges(7,5) * t393;
t82 = Ifges(7,4) * t295;
t35 = -Ifges(7,1) * t83 + Ifges(7,5) * t215 + t82;
t402 = t35 / 0.2e1;
t401 = -t129 * Ifges(6,4) / 0.2e1 - t128 * Ifges(6,2) / 0.2e1 + Ifges(6,6) * t449;
t400 = -t295 / 0.2e1;
t399 = t295 / 0.2e1;
t398 = t83 / 0.2e1;
t397 = -t83 / 0.2e1;
t395 = t128 / 0.2e1;
t394 = t129 / 0.2e1;
t392 = t184 / 0.2e1;
t391 = -t215 / 0.2e1;
t390 = t215 / 0.2e1;
t388 = Ifges(7,4) * t83;
t229 = t246 * pkin(8);
t378 = Ifges(4,4) * t243;
t377 = Ifges(4,4) * t246;
t355 = t236 * t246;
t353 = t237 * t247;
t348 = t246 * t247;
t203 = t246 * pkin(4) + t229;
t342 = qJD(2) * t244;
t329 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t173;
t141 = t153 * pkin(2);
t321 = -t141 + t436;
t142 = t155 * pkin(2);
t320 = -t142 + t435;
t317 = t237 * t342;
t316 = t247 * t344;
t306 = pkin(8) * t154 - t141;
t305 = pkin(8) * t156 - t142;
t298 = -t335 / 0.2e1;
t119 = t187 * t243 - t211;
t284 = t246 * Ifges(4,2) + t378;
t280 = Ifges(6,5) * t236 + Ifges(6,6) * t238;
t100 = t160 * t236 - t238 * t353;
t99 = t160 * t238 + t236 * t353;
t44 = t100 * t245 + t242 * t99;
t43 = -t100 * t242 + t245 * t99;
t161 = t239 * t243 + t246 * t354;
t267 = t243 * (Ifges(4,1) * t246 - t378);
t137 = t425 * t246;
t102 = t154 * t246 - t243 * t304;
t104 = t156 * t246 + t243 * t303;
t262 = -g(1) * t104 - g(2) * t102 - g(3) * t161;
t250 = t243 * (Ifges(6,3) * t246 + t243 * t280);
t36 = -pkin(4) * t183 + qJDD(5) - t45;
t181 = t396 * t339;
t180 = -qJ(4) * t341 + t222;
t179 = t291 * qJD(2);
t178 = t286 * qJD(2);
t164 = Ifges(4,6) * qJD(3) + qJD(2) * t284;
t152 = pkin(5) * t350 + t203;
t150 = -qJ(4) * t338 + t224 - t337;
t138 = t275 * t246;
t134 = t380 * t339;
t125 = mrSges(6,1) * t343 + mrSges(6,3) * t171;
t124 = -mrSges(6,2) * t343 + mrSges(6,3) * t172;
t113 = -qJD(3) * pkin(3) + qJD(4) + t119;
t112 = -mrSges(5,2) * t183 - mrSges(5,3) * t184;
t111 = mrSges(4,1) * t183 + mrSges(4,2) * t184;
t106 = -qJD(3) * t322 + (t316 + t340) * t246;
t105 = qJD(3) * t161 + t243 * t316;
t91 = -t169 * t236 + t177;
t75 = t211 + (-qJD(2) * t220 - t187) * t243;
t71 = t158 * t246 + t339 * t425;
t70 = qJD(3) * t264 - qJD(6) * t137;
t68 = t105 * t236 + t238 * t317;
t67 = t105 * t238 - t236 * t317;
t66 = mrSges(7,1) * t215 + mrSges(7,3) * t83;
t65 = -mrSges(7,2) * t215 + mrSges(7,3) * t295;
t58 = pkin(3) * t183 - qJD(2) * t337 + t253;
t50 = t129 * Ifges(6,1) + t128 * Ifges(6,4) + t184 * Ifges(6,5);
t34 = Ifges(7,2) * t295 + Ifges(7,6) * t215 - t388;
t20 = -pkin(5) * t128 + t36;
t19 = -mrSges(7,2) * t173 + mrSges(7,3) * t31;
t18 = mrSges(7,1) * t173 - mrSges(7,3) * t30;
t15 = -qJD(6) * t44 - t242 * t68 + t245 * t67;
t14 = qJD(6) * t43 + t242 * t67 + t245 * t68;
t3 = [m(2) * qJDD(1) + t100 * t76 + t68 * t124 + t67 * t125 + t14 * t65 + t15 * t66 + t43 * t18 + t44 * t19 + t99 * t77 + t433 * t160 - t457 * t105 + (t434 + t441) * t161 + (t195 - t330) * t106 + (-m(2) - m(3) + t455) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t248 - t111 - t112) * t247 + (-mrSges(3,1) * t248 - mrSges(3,2) * qJDD(2) + (t178 - t179) * qJD(2)) * t244) * t237 + m(3) * (qJDD(1) * t239 ^ 2 + (t143 * t247 + t144 * t244) * t237) + m(5) * (t105 * t113 - t106 * t114 + t160 * t47 - t161 * t45 + (t122 * t342 - t247 * t58) * t237) + m(4) * (t105 * t119 + t106 * t120 - t160 * t52 + t161 * t51 + (-t130 * t247 + t188 * t342) * t237) + m(7) * (t1 * t44 + t106 * t61 + t14 * t9 + t15 * t8 + t161 * t20 + t2 * t43) + m(6) * (t100 * t13 + t106 * t87 + t12 * t99 + t161 * t36 + t39 * t67 + t40 * t68); (Ifges(7,1) * t70 + Ifges(7,4) * t71) * t397 - t50 * t355 / 0.2e1 + t12 * (mrSges(6,1) * t243 + mrSges(6,3) * t355) + t423 * mrSges(5,1) + (Ifges(7,5) * t70 + Ifges(7,6) * t71) * t390 + ((m(4) * (t119 * t246 - t120 * t243) + m(5) * (t113 * t246 + t114 * t243)) * pkin(8) - t411 + t427 * t463) * qJD(3) + t13 * (-mrSges(6,2) * t243 - mrSges(6,3) * t350) + t424 * mrSges(4,3) + (t204 + t143) * mrSges(3,1) + (t415 / 0.2e1 + t114 * mrSges(5,1) - t120 * mrSges(4,3) - t429 * pkin(8) - t164 / 0.2e1) * t339 + (-Ifges(7,5) * t138 - Ifges(7,6) * t137 + Ifges(7,3) * t243) * t393 + (-Ifges(7,4) * t138 - Ifges(7,2) * t137 + Ifges(7,6) * t243) * t404 + (-Ifges(7,1) * t138 - Ifges(7,4) * t137 + Ifges(7,5) * t243) * t405 + t20 * (mrSges(7,1) * t137 - mrSges(7,2) * t138) + (-t1 * t137 + t138 * t2 - t70 * t8 + t71 * t9) * mrSges(7,3) + t457 * (-pkin(8) * t338 + t243 * t318) + t179 * t319 + (t205 - t144) * mrSges(3,2) + (t246 * (-Ifges(4,2) * t243 + t377) + t267 + t250) * t335 / 0.2e1 + (t122 * t150 + t189 * t58 + t423 * pkin(8) - (t122 * t244 + (t113 * t243 - t114 * t246) * t247) * t346) * m(5) + (-pkin(2) * t130 + t424 * pkin(8) - (t188 * t244 + (t119 * t243 + t120 * t246) * t247) * t346) * m(4) + t428 * t298 + (t419 / 0.2e1 + t113 * mrSges(5,1) + t119 * mrSges(4,3) + Ifges(7,5) * t397 + Ifges(7,6) * t399 + Ifges(7,3) * t390 - t456) * t338 + (Ifges(6,5) * t243 - t246 * t285) * t394 + (Ifges(6,6) * t243 - t246 * t282) * t395 + t350 * t401 + t70 * t402 - t138 * t406 - t137 * t407 + (Ifges(7,4) * t70 + Ifges(7,2) * t71) * t399 + Ifges(3,3) * qJDD(2) - t130 * t291 - t183 * t284 / 0.2e1 + t58 * t286 + t183 * t278 / 0.2e1 + t246 * (Ifges(4,4) * t184 - Ifges(4,2) * t183 + Ifges(4,6) * qJDD(3)) / 0.2e1 - t246 * (Ifges(5,5) * qJDD(3) - Ifges(5,6) * t184 + Ifges(5,3) * t183) / 0.2e1 - t243 * (Ifges(5,4) * qJDD(3) - Ifges(5,2) * t184 + Ifges(5,6) * t183) / 0.2e1 + t203 * t64 + t189 * t112 - t181 * t94 + t152 * t11 - t134 * t41 - pkin(2) * t111 + t91 * t77 + t92 * t76 + t243 * t460 - t243 * t461 + t61 * (-mrSges(7,1) * t71 + mrSges(7,2) * t70) + t71 * t34 / 0.2e1 + t23 * t19 + t22 * t18 + (-t195 + t413) * t246 * t318 + t433 * pkin(8) * t243 + t434 * t229 + (-m(5) * (t305 + t435) - m(4) * t305 - m(7) * t320 - m(6) * (-qJ(5) * t359 + t320) + mrSges(6,3) * t359 + t408 * t156 + t409 * t155) * g(1) + (-m(5) * (t306 + t436) - m(4) * t306 - m(7) * t321 - m(6) * (-qJ(5) * t360 + t321) + mrSges(6,3) * t360 + t408 * t154 + t409 * t153) * g(2) + t437 * t124 + t438 * t125 + (t12 * t91 + t13 * t92 - t181 * t87 + t203 * t36 + t39 * t438 + t40 * t437) * m(6) + (t455 * (pkin(2) * t353 + pkin(8) * t354) + (t421 * (pkin(3) * t348 + qJ(4) * t349) - t270 * mrSges(6,1) - t271 * mrSges(6,2) + t417 * t348 + t453 * t247 + (-m(6) * pkin(4) - m(7) * t220 + t450) * t244) * t237) * g(3) + (t243 * t443 - t246 * t442) * qJDD(3) / 0.2e1 + (t243 * t444 - t246 * t280 + t377) * t392 + (-Ifges(4,4) * t183 + Ifges(4,5) * qJDD(3) + Ifges(6,5) * t129 + Ifges(6,6) * t128 + t184 * t444 + t329) * t243 / 0.2e1 + t447 * t65 + (t1 * t23 - t134 * t61 + t152 * t20 + t2 * t22 + t447 * t9 + t448 * t8) * m(7) + t448 * t66 + t36 * t290 * t246 + (-t319 + t150) * t178 + t279 * t449; t411 * qJD(2) - t415 * t343 / 0.2e1 + t164 * t343 / 0.2e1 - t113 * t325 - t114 * t326 - (-Ifges(4,2) * t343 + t221 + t419) * t341 / 0.2e1 + (Ifges(7,4) * t425 - Ifges(7,2) * t275) * t404 + (Ifges(7,1) * t425 - Ifges(7,4) * t275) * t405 + t20 * (mrSges(7,1) * t275 + mrSges(7,2) * t425) + (Ifges(7,5) * t425 - Ifges(7,6) * t275) * t393 + t425 * t406 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t400 + (-t108 * t87 - t39 * t56 - t40 * t57 + qJ(4) * t36 - t277 * t381 + (-t236 * t40 - t238 * t39) * qJD(5)) * m(6) - t422 * t381 + (t160 * t451 + t161 * t414) * g(3) + (t101 * t451 + t102 * t414) * g(2) + (t103 * t451 + t104 * t414) * g(1) + t452 * mrSges(7,3) + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t398 + (Ifges(7,5) * t127 + Ifges(7,6) * t126) * t391 - t275 * t407 + (-t147 + t64) * qJ(4) - t277 * mrSges(6,3) + t427 * t298 + (Ifges(7,5) * t398 + Ifges(7,6) * t400 + Ifges(7,3) * t391 + t456) * t341 + (-t124 * t236 - t125 * t238) * qJD(5) + (-Ifges(7,1) * t158 - Ifges(7,4) * t157) * t397 + (-Ifges(7,4) * t158 - Ifges(7,2) * t157) * t399 + (-Ifges(7,5) * t158 - Ifges(7,6) * t157) * t390 + (Ifges(6,1) * t238 - t376) * t394 + (-Ifges(6,2) * t236 + t375) * t395 + t236 * t401 - t158 * t402 + (Ifges(6,5) * t238 - Ifges(6,6) * t236) * t392 + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + t36 * t289 + t238 * t50 / 0.2e1 + t216 * t11 - t180 * t178 - pkin(3) * t148 - t127 * t35 / 0.2e1 - t57 * t124 - t56 * t125 - t108 * t94 + t97 * t18 + t98 * t19 - t75 * t41 + (-m(5) * t113 + t324 + t457) * t120 - t51 * mrSges(4,2) + t52 * mrSges(4,1) - t45 * mrSges(5,3) + t47 * mrSges(5,2) + (-t267 / 0.2e1 - t250 / 0.2e1 + t428 / 0.2e1) * t248 - t432 * t34 / 0.2e1 + (mrSges(7,1) * t432 + mrSges(7,2) * t431) * t61 + t439 * t66 + t440 * t65 + (t1 * t98 + t2 * t97 + t20 * t216 + t439 * t8 + t440 * t9 - t61 * t75) * m(7) + t442 * t183 + t443 * t184 + (-pkin(3) * t47 - qJ(4) * t45 - t122 * t180) * m(5) + (-t323 + t429 - t446) * t119 + (-t413 - t446) * qJD(4); t275 * t19 + t425 * t18 + t431 * t66 + t432 * t65 + t330 * qJD(3) + (t124 * t238 - t125 * t236 + t178) * t343 + t148 + (-qJD(3) * t61 - t452) * m(7) + (t263 - qJD(3) * t87 - (t236 * t39 - t238 * t40) * t343 + t277) * m(6) + (qJD(3) * t114 + t122 * t343 + t263 + t47) * m(5) + t422; -t172 * t124 - t171 * t125 - t295 * t65 - t83 * t66 + (-t295 * t9 - t8 * t83 + t20 + t262) * m(7) + (-t171 * t39 - t172 * t40 + t262 + t36) * m(6) + t441; -t461 + t460 - t61 * (-mrSges(7,1) * t83 + mrSges(7,2) * t295) + (Ifges(7,1) * t295 + t388) * t398 + t34 * t397 + (Ifges(7,5) * t295 + Ifges(7,6) * t83) * t391 - t8 * t65 + t9 * t66 - g(1) * ((t103 * t227 - t155 * t226) * mrSges(7,1) + (-t103 * t226 - t155 * t227) * mrSges(7,2)) - g(2) * ((t101 * t227 - t153 * t226) * mrSges(7,1) + (-t101 * t226 - t153 * t227) * mrSges(7,2)) - g(3) * ((t160 * t227 + t226 * t353) * mrSges(7,1) + (-t160 * t226 + t227 * t353) * mrSges(7,2)) + (t295 * t8 - t83 * t9) * mrSges(7,3) + t329 + (Ifges(7,2) * t83 + t35 + t82) * t400;];
tau  = t3;
