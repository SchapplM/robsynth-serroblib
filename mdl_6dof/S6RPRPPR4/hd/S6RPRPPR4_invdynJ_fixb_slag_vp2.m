% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:46:47
% EndTime: 2019-03-09 02:47:17
% DurationCPUTime: 21.38s
% Computational Cost: add. (8076->677), mult. (19590->854), div. (0->0), fcn. (14780->12), ass. (0->298)
t445 = -mrSges(6,2) - mrSges(5,3);
t444 = Ifges(5,1) + Ifges(6,1);
t443 = -Ifges(5,4) + Ifges(6,5);
t432 = -Ifges(5,5) - Ifges(6,4);
t431 = Ifges(5,6) - Ifges(6,6);
t430 = Ifges(5,3) + Ifges(6,2);
t300 = qJD(1) * qJD(2);
t210 = qJ(2) * qJDD(1) + t300;
t229 = sin(pkin(9));
t234 = sin(qJ(3));
t231 = cos(pkin(9));
t354 = cos(qJ(3));
t293 = t354 * t231;
t246 = -t234 * t229 + t293;
t183 = t246 * qJD(3);
t192 = t229 * t354 + t234 * t231;
t142 = qJD(1) * t183 + qJDD(1) * t192;
t228 = sin(pkin(10));
t230 = cos(pkin(10));
t120 = -t230 * qJDD(3) + t142 * t228;
t373 = t120 / 0.2e1;
t121 = qJDD(3) * t228 + t142 * t230;
t372 = t121 / 0.2e1;
t184 = t192 * qJD(3);
t299 = qJDD(1) * t229;
t143 = qJD(1) * t184 - qJDD(1) * t293 + t234 * t299;
t442 = -t143 / 0.2e1;
t370 = t143 / 0.2e1;
t273 = -mrSges(3,1) * t231 + mrSges(3,2) * t229;
t441 = -m(3) * pkin(1) - mrSges(2,1) + t273;
t306 = t229 ^ 2 + t231 ^ 2;
t344 = pkin(7) + qJ(2);
t201 = t344 * t231;
t194 = qJD(1) * t201;
t177 = t234 * t194;
t199 = t344 * t229;
t193 = qJD(1) * t199;
t178 = t354 * t193;
t144 = -t177 - t178;
t180 = t192 * qJD(1);
t153 = -t230 * qJD(3) + t180 * t228;
t179 = t246 * qJD(1);
t218 = pkin(2) * t231 + pkin(1);
t196 = -qJD(1) * t218 + qJD(2);
t251 = qJD(3) * pkin(3) - qJD(4) + t144;
t252 = qJD(3) * t228 + t230 * t180;
t267 = mrSges(6,1) * t228 - mrSges(6,3) * t230;
t269 = mrSges(5,1) * t228 + mrSges(5,2) * t230;
t355 = t228 / 0.2e1;
t400 = Ifges(4,5) * qJD(3);
t243 = qJ(5) * t252 + t251;
t60 = pkin(4) * t153 - t243;
t440 = t196 * mrSges(4,2) - t228 * (Ifges(5,4) * t252 - t153 * Ifges(5,2) - t179 * Ifges(5,6)) / 0.2e1 + (Ifges(6,5) * t252 - t179 * Ifges(6,6) + t153 * Ifges(6,3)) * t355 - t144 * mrSges(4,3) - t251 * t269 + t60 * t267 + t400 / 0.2e1;
t233 = sin(qJ(6));
t236 = cos(qJ(6));
t104 = -pkin(3) * t179 - qJ(4) * t180 + t196;
t145 = -t234 * t193 + t194 * t354;
t136 = qJD(3) * qJ(4) + t145;
t61 = t104 * t230 - t228 * t136;
t274 = qJD(5) - t61;
t375 = pkin(4) + pkin(5);
t30 = -pkin(8) * t252 + t179 * t375 + t274;
t62 = t228 * t104 + t230 * t136;
t48 = -t179 * qJ(5) + t62;
t35 = pkin(8) * t153 + t48;
t10 = t233 * t30 + t236 * t35;
t399 = Ifges(4,6) * qJD(3);
t47 = pkin(4) * t179 + t274;
t9 = -t233 * t35 + t236 * t30;
t439 = -t196 * mrSges(4,1) - t61 * mrSges(5,1) + t47 * mrSges(6,1) + t9 * mrSges(7,1) + t62 * mrSges(5,2) - t10 * mrSges(7,2) + t145 * mrSges(4,3) - t48 * mrSges(6,3) + t399 / 0.2e1;
t438 = m(7) * pkin(8) + mrSges(7,3) + t445;
t268 = t230 * mrSges(6,1) + t228 * mrSges(6,3);
t270 = mrSges(5,1) * t230 - mrSges(5,2) * t228;
t191 = t228 * t236 - t230 * t233;
t253 = t228 * t233 + t230 * t236;
t414 = t253 * mrSges(7,1) + mrSges(7,2) * t191;
t437 = t268 + t270 + t414;
t436 = -m(6) - m(7);
t235 = sin(qJ(1));
t434 = g(2) * t235;
t433 = mrSges(6,3) - mrSges(5,2);
t429 = t120 * t443 + t121 * t444 - t143 * t432;
t428 = -t153 * t431 - t179 * t430 - t252 * t432;
t427 = t153 * t443 + t179 * t432 + t252 * t444;
t176 = qJD(6) + t179;
t339 = Ifges(4,4) * t180;
t413 = t153 * t233 + t236 * t252;
t93 = t153 * t236 - t233 * t252;
t424 = Ifges(7,5) * t413 + t179 * Ifges(4,2) + t93 * Ifges(7,6) + t176 * Ifges(7,3) + t339 + t399;
t401 = qJD(3) * mrSges(4,1) - mrSges(5,1) * t153 - mrSges(5,2) * t252 - mrSges(4,3) * t180;
t423 = -t228 * t431 - t230 * t432;
t336 = Ifges(6,5) * t228;
t338 = Ifges(5,4) * t228;
t422 = t230 * t444 + t336 - t338;
t151 = -t234 * t199 + t201 * t354;
t112 = qJD(2) * t192 + qJD(3) * t151;
t319 = t192 * t230;
t421 = qJD(5) * t319 - t112;
t286 = m(3) * qJ(2) + mrSges(3,3);
t419 = -t286 - mrSges(4,3) + mrSges(2,2);
t418 = Ifges(6,5) * t372 + Ifges(6,6) * t370 - t121 * Ifges(5,4) / 0.2e1 + Ifges(5,6) * t442 + (Ifges(6,3) + Ifges(5,2)) * t373;
t226 = pkin(9) + qJ(3);
t222 = sin(t226);
t223 = cos(t226);
t272 = t223 * mrSges(4,1) - t222 * mrSges(4,2);
t416 = t222 * t445 - t272;
t412 = m(7) * pkin(5) + mrSges(5,1) + mrSges(6,1);
t195 = -qJDD(1) * t218 + qJDD(2);
t54 = pkin(3) * t143 - qJ(4) * t142 - qJD(4) * t180 + t195;
t281 = pkin(7) * qJDD(1) + t210;
t169 = t281 * t231;
t256 = t281 * t229;
t294 = -qJD(3) * t178 + t354 * t169 - t234 * t256;
t69 = qJDD(3) * qJ(4) + (qJD(4) - t177) * qJD(3) + t294;
t23 = t228 * t54 + t230 * t69;
t14 = t143 * qJ(5) - t179 * qJD(5) + t23;
t11 = pkin(8) * t120 + t14;
t22 = -t228 * t69 + t230 * t54;
t277 = qJDD(5) - t22;
t8 = -pkin(8) * t121 - t143 * t375 + t277;
t1 = qJD(6) * t9 + t11 * t236 + t233 * t8;
t2 = -qJD(6) * t10 - t11 * t233 + t236 * t8;
t411 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t28 = qJD(6) * t93 + t120 * t233 + t121 * t236;
t383 = t28 / 0.2e1;
t29 = -qJD(6) * t413 + t120 * t236 - t121 * t233;
t382 = t29 / 0.2e1;
t351 = Ifges(7,4) * t413;
t32 = Ifges(7,2) * t93 + Ifges(7,6) * t176 + t351;
t381 = t32 / 0.2e1;
t92 = Ifges(7,4) * t93;
t33 = Ifges(7,1) * t413 + Ifges(7,5) * t176 + t92;
t380 = t33 / 0.2e1;
t140 = qJDD(6) - t143;
t371 = t140 / 0.2e1;
t71 = -mrSges(6,2) * t120 + mrSges(6,3) * t143;
t72 = -mrSges(5,2) * t143 - mrSges(5,3) * t120;
t407 = t71 + t72;
t73 = mrSges(5,1) * t143 - mrSges(5,3) * t121;
t74 = -t143 * mrSges(6,1) + t121 * mrSges(6,2);
t406 = t73 - t74;
t43 = -mrSges(7,1) * t93 + mrSges(7,2) * t413;
t98 = mrSges(6,1) * t153 - mrSges(6,3) * t252;
t405 = -t98 + t43;
t343 = pkin(8) - qJ(4);
t198 = t343 * t228;
t200 = t343 * t230;
t150 = -t198 * t233 - t200 * t236;
t134 = t228 * t144;
t138 = pkin(3) * t180 - qJ(4) * t179;
t36 = t134 + (-pkin(8) * t179 - t138) * t230 - t375 * t180;
t324 = t179 * t228;
t81 = t228 * t138 + t230 * t144;
t55 = t180 * qJ(5) + t81;
t46 = pkin(8) * t324 + t55;
t403 = qJD(4) * t191 - qJD(6) * t150 + t233 * t46 - t236 * t36;
t148 = -t198 * t236 + t200 * t233;
t402 = qJD(4) * t253 + qJD(6) * t148 - t233 * t36 - t236 * t46;
t118 = t191 * t179;
t182 = t191 * qJD(6);
t395 = -t118 - t182;
t119 = t253 * t179;
t181 = t253 * qJD(6);
t394 = -t119 - t181;
t393 = -t354 * t199 - t234 * t201;
t392 = -t22 * t228 + t23 * t230;
t237 = cos(qJ(1));
t391 = g(1) * t237 + t434;
t302 = m(5) - t436;
t326 = qJ(5) * t228;
t285 = pkin(3) + t326;
t350 = pkin(4) * t230;
t197 = -t285 - t350;
t349 = pkin(5) * t230;
t389 = t438 * t223 + (-m(7) * (t197 - t349) - m(6) * t197 + m(5) * pkin(3) + t437) * t222;
t325 = qJ(5) * t230;
t387 = t228 * t375 - t325;
t386 = qJ(4) * t302;
t79 = -qJD(3) * t145 - t234 * t169 - t354 * t256;
t385 = Ifges(7,4) * t383 + Ifges(7,2) * t382 + Ifges(7,6) * t371;
t384 = Ifges(7,1) * t383 + Ifges(7,4) * t382 + Ifges(7,5) * t371;
t379 = -t93 / 0.2e1;
t378 = t93 / 0.2e1;
t377 = -t413 / 0.2e1;
t376 = t413 / 0.2e1;
t374 = -t120 / 0.2e1;
t369 = -t153 / 0.2e1;
t368 = t153 / 0.2e1;
t367 = -t252 / 0.2e1;
t366 = t252 / 0.2e1;
t365 = -t176 / 0.2e1;
t364 = t176 / 0.2e1;
t363 = t179 / 0.2e1;
t362 = -t179 / 0.2e1;
t361 = -t180 / 0.2e1;
t360 = t180 / 0.2e1;
t352 = mrSges(7,3) * t93;
t346 = g(3) * t222;
t345 = t10 * mrSges(7,3);
t213 = t223 * pkin(3);
t337 = Ifges(5,4) * t230;
t335 = Ifges(6,5) * t230;
t100 = pkin(3) * t184 - qJ(4) * t183 - qJD(4) * t192;
t111 = t246 * qJD(2) + qJD(3) * t393;
t52 = t228 * t100 + t230 * t111;
t323 = t179 * t230;
t322 = t183 * t228;
t321 = t183 * t230;
t320 = t192 * t228;
t212 = t222 * qJ(4);
t318 = t222 * t237;
t317 = t223 * t237;
t316 = t228 * t235;
t315 = t230 * t237;
t314 = t344 * t237;
t311 = t235 * t230;
t310 = t237 * t228;
t113 = -mrSges(6,2) * t153 - mrSges(6,3) * t179;
t114 = mrSges(5,2) * t179 - mrSges(5,3) * t153;
t309 = t113 + t114;
t115 = -mrSges(5,1) * t179 - mrSges(5,3) * t252;
t116 = mrSges(6,1) * t179 + mrSges(6,2) * t252;
t308 = t115 - t116;
t141 = -pkin(3) * t246 - qJ(4) * t192 - t218;
t90 = t228 * t141 + t230 * t151;
t307 = t213 + t212;
t305 = qJD(4) * t228;
t304 = qJD(4) * t230;
t303 = qJD(5) * t228;
t298 = qJDD(1) * t231;
t297 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t140;
t67 = -qJ(5) * t246 + t90;
t284 = -t218 - t213;
t283 = t143 * mrSges(4,1) + t142 * mrSges(4,2);
t59 = t120 * mrSges(5,1) + t121 * mrSges(5,2);
t58 = t120 * mrSges(6,1) - t121 * mrSges(6,3);
t102 = t228 * t111;
t51 = t100 * t230 - t102;
t80 = t138 * t230 - t134;
t146 = t228 * t151;
t89 = t141 * t230 - t146;
t279 = t237 * t218 + t235 * t344;
t34 = t184 * qJ(5) - qJD(5) * t246 + t52;
t278 = t307 + (t326 + t350) * t223;
t7 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t276 = -mrSges(3,1) * t298 + mrSges(3,2) * t299;
t264 = -Ifges(5,2) * t228 + t337;
t261 = Ifges(6,3) * t228 + t335;
t260 = pkin(4) * t228 - t325;
t258 = t228 * t61 - t230 * t62;
t44 = t146 + (-pkin(8) * t192 - t141) * t230 + t375 * t246;
t53 = pkin(8) * t320 + t67;
t17 = -t233 * t53 + t236 * t44;
t18 = t233 * t44 + t236 * t53;
t76 = -mrSges(7,2) * t176 + t352;
t77 = mrSges(7,1) * t176 - mrSges(7,3) * t413;
t257 = -t233 * t77 + t236 * t76;
t165 = t223 * t316 + t315;
t166 = t223 * t311 - t310;
t255 = t165 * t236 - t166 * t233;
t254 = -t165 * t233 - t166 * t236;
t250 = pkin(3) * t317 + qJ(4) * t318 + t279;
t127 = t191 * t192;
t167 = t223 * t310 - t311;
t244 = -g(1) * t167 - g(2) * t165 - t228 * t346;
t70 = -qJDD(3) * pkin(3) + qJDD(4) - t79;
t19 = t120 * pkin(4) - t121 * qJ(5) - qJD(5) * t252 + t70;
t219 = -qJDD(1) * pkin(1) + qJDD(2);
t185 = t230 * t375 + t285;
t174 = Ifges(4,4) * t179;
t168 = t223 * t315 + t316;
t157 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t179;
t131 = t180 * Ifges(4,1) + t174 + t400;
t128 = t253 * t192;
t123 = t167 * t233 + t168 * t236;
t122 = t167 * t236 - t168 * t233;
t91 = t192 * t260 - t393;
t82 = t179 * t260 + t145;
t78 = -qJD(3) * t177 + t294;
t75 = -t192 * t387 + t393;
t68 = pkin(4) * t246 - t89;
t66 = -t181 * t192 + t183 * t191;
t65 = qJD(6) * t127 + t183 * t253;
t57 = -pkin(4) * t180 - t80;
t56 = -t179 * t387 - t145;
t50 = t183 * t260 - t421;
t45 = -t153 * t375 + t243;
t42 = -pkin(4) * t184 - t51;
t41 = -t183 * t387 + t421;
t25 = pkin(8) * t322 + t34;
t24 = t102 + (-pkin(8) * t183 - t100) * t230 - t375 * t184;
t21 = -mrSges(7,2) * t140 + mrSges(7,3) * t29;
t20 = mrSges(7,1) * t140 - mrSges(7,3) * t28;
t16 = -pkin(4) * t143 + t277;
t15 = t120 * pkin(5) + t19;
t4 = -qJD(6) * t18 - t233 * t25 + t236 * t24;
t3 = qJD(6) * t17 + t233 * t24 + t236 * t25;
t5 = [(-m(4) * t144 - m(5) * t251 - t401) * t112 + (t195 * mrSges(4,2) - t79 * mrSges(4,3) + Ifges(4,1) * t142 - Ifges(4,4) * t143 + Ifges(4,5) * qJDD(3) + t19 * t267 + t261 * t373 + t264 * t374 + t269 * t70 + t370 * t423 + t372 * t422) * t192 + t418 * t320 + (t1 * t127 + t10 * t66 - t128 * t2 - t65 * t9) * mrSges(7,3) + (-t22 * t319 - t23 * t320 - t321 * t61 - t322 * t62) * mrSges(5,3) + (-t14 * t320 + t16 * t319 + t321 * t47 - t322 * t48) * mrSges(6,2) + (-m(5) * t314 - t254 * mrSges(7,1) + t255 * mrSges(7,2) + t436 * (-t166 * pkin(4) - qJ(5) * t165 + t314) + t412 * t166 + t433 * t165 + (-m(4) * t344 + t419) * t237 + (-m(7) * t284 - (m(7) * t343 + mrSges(7,3)) * t222 + m(4) * t218 + (-m(5) - m(6)) * (t284 - t212) - t416 - t441) * t235) * g(1) + (-m(4) * t279 - m(5) * t250 - t123 * mrSges(7,1) - t122 * mrSges(7,2) + t436 * (t168 * pkin(4) + qJ(5) * t167 + t250) - t412 * t168 - t433 * t167 + t438 * t318 + (-t272 + t441) * t237 + t419 * t235) * g(2) + (t131 / 0.2e1 + Ifges(4,1) * t360 + Ifges(4,4) * t363 + t261 * t368 + t264 * t369 + t422 * t366 + t423 * t362 + t440) * t183 + (-t424 / 0.2e1 + t428 / 0.2e1 + t430 * t362 - t432 * t366 - Ifges(4,4) * t360 - Ifges(4,2) * t363 - Ifges(7,3) * t364 + Ifges(6,6) * t368 + Ifges(5,6) * t369 - Ifges(7,5) * t376 - Ifges(7,6) * t378 - t439) * t184 + m(6) * (t14 * t67 + t16 * t68 + t19 * t91 + t34 * t48 + t42 * t47 + t50 * t60) + m(7) * (t1 * t18 + t10 * t3 - t15 * t75 + t17 * t2 + t4 * t9 + t41 * t45) + m(5) * (t22 * t89 + t23 * t90 + t51 * t61 + t52 * t62) + m(4) * (t111 * t145 + t151 * t78 - t195 * t218) + (Ifges(7,1) * t65 + Ifges(7,4) * t66) * t376 + (Ifges(7,1) * t128 + Ifges(7,4) * t127) * t383 + (Ifges(7,4) * t65 + Ifges(7,2) * t66) * t378 + (Ifges(7,4) * t128 + Ifges(7,2) * t127) * t382 - (-m(4) * t79 + m(5) * t70 - qJDD(3) * mrSges(4,1) + mrSges(4,3) * t142 + t59) * t393 + (t297 / 0.2e1 + Ifges(4,6) * qJDD(3) - Ifges(4,2) * t143 + Ifges(4,4) * t142 + t78 * mrSges(4,3) + Ifges(7,3) * t371 - Ifges(5,6) * t374 + Ifges(7,6) * t382 + Ifges(7,5) * t383 + t23 * mrSges(5,2) - t14 * mrSges(6,3) - t22 * mrSges(5,1) + t16 * mrSges(6,1) - t195 * mrSges(4,1) + t411 + 0.2e1 * t432 * t372 + (-t370 + t442) * t430 + (-Ifges(6,6) + t431) * t373) * t246 + m(3) * (-pkin(1) * t219 + (t210 + t300) * qJ(2) * t306) + Ifges(2,3) * qJDD(1) - t218 * t283 + t427 * t321 / 0.2e1 + t429 * t319 / 0.2e1 + t45 * (-mrSges(7,1) * t66 + mrSges(7,2) * t65) + t67 * t71 + t68 * t74 + t75 * t7 + t3 * t76 + t4 * t77 + t89 * t73 + t90 * t72 + t91 * t58 + 0.2e1 * t306 * t210 * mrSges(3,3) + (Ifges(7,5) * t65 + Ifges(7,6) * t66) * t364 + (Ifges(7,5) * t128 + Ifges(7,6) * t127) * t371 + t50 * t98 + t34 * t113 + t52 * t114 + t51 * t115 + t42 * t116 - t15 * (-mrSges(7,1) * t127 + mrSges(7,2) * t128) + (Ifges(3,4) * t229 + Ifges(3,2) * t231) * t298 + (Ifges(3,1) * t229 + Ifges(3,4) * t231) * t299 + t151 * (-qJDD(3) * mrSges(4,2) - mrSges(4,3) * t143) + t17 * t20 + t18 * t21 + t41 * t43 + t111 * t157 + t219 * t273 - pkin(1) * t276 + t65 * t380 + t66 * t381 + t128 * t384 + t127 * t385; t283 + t406 * t230 + t407 * t228 + (t405 + t401) * t180 + m(3) * t219 + t395 * t77 + t394 * t76 + t276 - (-t228 * t308 + t230 * t309 + t157) * t179 - t253 * t20 + t191 * t21 + (-g(1) * t235 + g(2) * t237) * (m(3) + m(4) + t302) - t286 * t306 * qJD(1) ^ 2 + (t1 * t191 + t10 * t394 + t180 * t45 - t2 * t253 + t395 * t9) * m(7) + (t14 * t228 - t16 * t230 - t180 * t60 + (-t228 * t47 - t230 * t48) * t179) * m(6) + (t179 * t258 + t180 * t251 + t22 * t230 + t23 * t228) * m(5) + (t144 * t180 - t145 * t179 + t195) * m(4); t424 * t360 - t15 * t414 + (t19 * t197 - t47 * t57 - t48 * t55 - t60 * t82) * m(6) + (-pkin(3) * t70 + qJ(4) * t392 - t258 * qJD(4) + t145 * t251 - t61 * t80 - t62 * t81) * m(5) + (Ifges(7,5) * t191 - Ifges(7,6) * t253) * t371 + (Ifges(7,4) * t191 - Ifges(7,2) * t253) * t382 + (Ifges(7,1) * t191 - Ifges(7,4) * t253) * t383 - t253 * t385 + (t237 * t389 - t317 * t386) * g(1) - (-Ifges(4,1) * t361 - t261 * t369 - t264 * t368 - t363 * t423 - t367 * t422 + t440) * t179 + (-m(7) * (-pkin(8) * t222 + t278) + t222 * mrSges(7,3) - m(5) * t307 - m(6) * t278 + (-m(7) * t349 - t437) * t223 + t416) * g(3) + (Ifges(7,4) * t119 + Ifges(7,2) * t118) * t379 + (-t335 + t337) * t372 + t338 * t374 - (mrSges(7,2) * t45 + Ifges(7,1) * t376 + Ifges(7,4) * t378 + Ifges(7,5) * t364 + t380) * t181 + (-t55 + t304) * t113 + (-t223 * t386 + t389) * t434 + (-t57 + t305) * t116 + t402 * t76 + (t1 * t150 + t148 * t2 - t15 * t185 + t403 * t9 + (t303 - t56) * t45 + t402 * t10) * m(7) + t403 * t77 + (-t81 + t304) * t114 + t336 * t373 + t405 * t303 + (-Ifges(7,5) * t377 - Ifges(4,2) * t362 + Ifges(5,6) * t368 + Ifges(6,6) * t369 - Ifges(7,6) * t379 - Ifges(7,3) * t365 + t430 * t363 - t367 * t432 + t439) * t180 + (-t1 * t253 - t10 * t118 - t191 * t2 - t394 * t9) * mrSges(7,3) + t401 * t145 + (Ifges(7,1) * t119 + Ifges(7,4) * t118) * t377 + t391 * (mrSges(4,1) * t222 + mrSges(4,2) * t223) + (-t323 * t47 + t324 * t48) * mrSges(6,2) + (t323 * t61 + t324 * t62 + t392) * mrSges(5,3) + (t174 + t131) * t362 + Ifges(4,3) * qJDD(3) + ((qJ(4) * t16 + qJD(4) * t47 - qJD(5) * t60) * m(6) - t406 * qJ(4) + t444 * t372 - t432 * t370 + t16 * mrSges(6,2)) * t228 - t427 * t323 / 0.2e1 + (-t339 + t428) * t361 + t429 * t355 + ((qJ(4) * t14 + qJD(4) * t48) * m(6) + t407 * qJ(4) - Ifges(6,3) * t373 + Ifges(5,2) * t374 + t431 * t370 + t14 * mrSges(6,2) - t418) * t230 - t78 * mrSges(4,2) + t79 * mrSges(4,1) + (Ifges(7,5) * t119 + Ifges(7,6) * t118) * t365 - t82 * t98 - t45 * (-mrSges(7,1) * t118 + mrSges(7,2) * t119) + Ifges(4,5) * t142 - Ifges(4,6) * t143 - t119 * t380 - t118 * t381 + t148 * t20 + t150 * t21 - t144 * t157 - t56 * t43 - pkin(3) * t59 + (-t305 - t80) * t115 + t185 * t7 + t197 * t58 - t19 * t268 - t70 * t270 - (-mrSges(7,1) * t45 + Ifges(7,4) * t376 + Ifges(7,2) * t378 + Ifges(7,6) * t364 + t345 + t381) * t182 + t191 * t384; t308 * t252 + t309 * t153 + t93 * t76 - t413 * t77 + t58 + t59 - t7 + (t10 * t93 - t413 * t9 + t15) * m(7) + (t153 * t48 - t252 * t47 + t19) * m(6) + (t153 * t62 + t252 * t61 + t70) * m(5) + (t223 * g(3) - t222 * t391) * t302; t236 * t20 + t233 * t21 - t405 * t252 + t257 * qJD(6) - (-t113 - t257) * t179 + t74 + (t1 * t233 - t252 * t45 + t2 * t236 + t244 + t176 * (t10 * t236 - t233 * t9)) * m(7) + (t179 * t48 + t252 * t60 + t16 + t244) * m(6); -t45 * (mrSges(7,1) * t413 + mrSges(7,2) * t93) + (Ifges(7,1) * t93 - t351) * t377 + t32 * t376 + (Ifges(7,5) * t93 - Ifges(7,6) * t413) * t365 + t413 * t345 + t10 * t77 - g(1) * (mrSges(7,1) * t122 - mrSges(7,2) * t123) - g(2) * (mrSges(7,1) * t255 + mrSges(7,2) * t254) - (mrSges(7,1) * t191 - mrSges(7,2) * t253) * t346 + t297 + (-t76 + t352) * t9 + (-Ifges(7,2) * t413 + t33 + t92) * t379 + t411;];
tau  = t5;
