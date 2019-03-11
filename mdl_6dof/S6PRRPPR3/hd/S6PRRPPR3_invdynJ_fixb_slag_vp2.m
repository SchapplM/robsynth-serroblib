% Calculate vector of inverse dynamics joint torques for
% S6PRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
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
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRPPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:08:11
% EndTime: 2019-03-08 21:08:35
% DurationCPUTime: 15.26s
% Computational Cost: add. (3503->623), mult. (7705->812), div. (0->0), fcn. (5058->10), ass. (0->290)
t393 = -m(7) - m(6);
t406 = mrSges(4,1) + mrSges(5,1);
t405 = Ifges(6,1) + Ifges(5,3);
t387 = Ifges(5,4) + Ifges(4,5);
t404 = -Ifges(4,6) + Ifges(5,6);
t403 = Ifges(6,6) + t387;
t194 = cos(qJ(3));
t300 = qJD(2) * t194;
t281 = mrSges(5,2) * t300;
t143 = qJD(3) * mrSges(5,3) + t281;
t190 = sin(qJ(6));
t193 = cos(qJ(6));
t297 = qJD(3) * t193;
t123 = t190 * t300 - t297;
t124 = qJD(3) * t190 + t193 * t300;
t277 = mrSges(6,3) * t300;
t378 = -qJD(3) * mrSges(6,1) + mrSges(7,1) * t123 + mrSges(7,2) * t124 + t277;
t283 = -t143 + t378;
t191 = sin(qJ(3));
t243 = pkin(5) * t191 + pkin(9) * t194;
t195 = cos(qJ(2));
t185 = sin(pkin(6));
t305 = qJD(1) * t185;
t157 = t195 * t305;
t132 = -qJD(2) * pkin(2) - t157;
t302 = qJD(2) * t191;
t73 = -pkin(3) * t300 - qJ(4) * t302 + t132;
t46 = pkin(4) * t300 + qJD(5) - t73;
t27 = qJD(2) * t243 + t46;
t357 = -pkin(4) - pkin(9);
t181 = -pkin(3) + t357;
t164 = qJ(5) * t302;
t192 = sin(qJ(2));
t273 = t192 * t305;
t131 = qJD(2) * pkin(8) + t273;
t187 = cos(pkin(6));
t304 = qJD(1) * t187;
t309 = t131 * t191 - t194 * t304;
t270 = qJD(4) + t309;
t245 = -t164 + t270;
t30 = qJD(3) * t181 + t245;
t10 = -t190 * t30 + t193 * t27;
t292 = qJD(2) * qJD(3);
t129 = -qJDD(2) * t194 + t191 * t292;
t130 = qJDD(2) * t191 + t194 * t292;
t174 = t191 * qJD(4);
t303 = qJD(2) * t185;
t267 = qJD(1) * t303;
t149 = t192 * t267;
t290 = qJDD(1) * t185;
t86 = t195 * t290 - t149;
t75 = -qJDD(2) * pkin(2) - t86;
t24 = pkin(3) * t129 - qJ(4) * t130 - qJD(2) * t174 + t75;
t210 = qJDD(5) - t24;
t5 = pkin(5) * t130 + t129 * t357 + t210;
t289 = qJDD(1) * t187;
t296 = qJD(3) * t194;
t299 = qJD(3) * t187;
t150 = t195 * t267;
t87 = t192 * t290 + t150;
t399 = qJDD(2) * pkin(8) + qJD(1) * t299 + t87;
t23 = -t131 * t296 - t191 * t399 + t194 * t289;
t206 = qJDD(4) - t23;
t291 = qJD(2) * qJD(5);
t198 = -qJ(5) * t130 - t191 * t291 + t206;
t8 = qJDD(3) * t181 + t198;
t1 = qJD(6) * t10 + t190 * t5 + t193 * t8;
t402 = t1 * mrSges(7,2);
t11 = t190 * t27 + t193 * t30;
t2 = -qJD(6) * t11 - t190 * t8 + t193 * t5;
t401 = t2 * mrSges(7,1);
t400 = Ifges(6,5) - t404;
t242 = mrSges(4,1) * t194 - mrSges(4,2) * t191;
t241 = mrSges(5,1) * t194 + mrSges(5,3) * t191;
t125 = t241 * qJD(2);
t240 = mrSges(6,1) * t191 - mrSges(6,2) * t194;
t128 = t240 * qJD(2);
t308 = -t125 - t128;
t398 = qJD(2) * t242 - t308;
t238 = t190 * mrSges(7,1) + t193 * mrSges(7,2);
t397 = mrSges(6,3) + t238 + mrSges(3,2) - mrSges(5,2) - mrSges(4,3);
t396 = m(7) * pkin(9);
t37 = qJD(6) * t123 - qJDD(3) * t190 + t129 * t193;
t395 = -t37 / 0.2e1;
t38 = qJD(6) * t124 - qJDD(3) * t193 - t129 * t190;
t394 = -t38 / 0.2e1;
t183 = qJD(3) * qJ(4);
t71 = t131 * t194 + t191 * t304;
t63 = t183 + t71;
t392 = m(5) * t63;
t56 = -qJ(5) * t300 + t71;
t41 = -t183 - t56;
t391 = m(6) * t41;
t121 = qJDD(6) + t130;
t390 = -t121 / 0.2e1;
t343 = pkin(8) - qJ(5);
t144 = t343 * t191;
t175 = t191 * qJ(4);
t133 = -pkin(3) * t194 - pkin(2) - t175;
t115 = pkin(4) * t194 - t133;
t72 = t243 + t115;
t29 = t144 * t193 + t190 * t72;
t312 = t191 * t195;
t212 = pkin(5) * t194 + t181 * t191;
t306 = qJ(4) * t296 + t174;
t40 = qJD(3) * t212 + t306;
t98 = -qJD(5) * t191 + t296 * t343;
t389 = -qJD(6) * t29 - t190 * t98 + t193 * t40 - (-t190 * t312 - t192 * t193) * t305;
t28 = -t144 * t190 + t193 * t72;
t310 = t193 * t195;
t388 = qJD(6) * t28 + t190 * t40 + t193 * t98 - (-t190 * t192 + t191 * t310) * t305;
t14 = -mrSges(7,1) * t38 + mrSges(7,2) * t37;
t88 = -qJDD(3) * mrSges(6,1) - mrSges(6,3) * t129;
t385 = t14 - t88;
t322 = sin(pkin(10));
t255 = t322 * t192;
t186 = cos(pkin(10));
t316 = t186 * t195;
t99 = t187 * t316 - t255;
t325 = t194 * t99;
t384 = pkin(3) * t325 + t175 * t99;
t254 = t322 * t195;
t317 = t186 * t192;
t101 = -t187 * t254 - t317;
t321 = t101 * t194;
t383 = pkin(3) * t321 + t101 * t175;
t91 = -qJDD(3) * mrSges(5,1) + mrSges(5,2) * t130;
t92 = qJDD(3) * mrSges(6,2) - mrSges(6,3) * t130;
t382 = t91 + t92;
t93 = -mrSges(5,2) * t129 + qJDD(3) * mrSges(5,3);
t381 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t129 + t93;
t39 = qJD(3) * pkin(5) - t41;
t380 = t238 * t39;
t342 = Ifges(4,4) * t191;
t233 = t194 * Ifges(4,2) + t342;
t161 = qJD(6) + t302;
t338 = Ifges(7,4) * t124;
t32 = Ifges(7,2) * t123 + Ifges(7,6) * t161 - t338;
t377 = Ifges(4,6) * qJD(3) + qJD(2) * t233 + t190 * t32;
t280 = mrSges(4,3) * t302;
t282 = mrSges(5,2) * t302;
t376 = qJD(3) * t406 - t280 - t282;
t279 = mrSges(4,3) * t300;
t142 = -qJD(3) * mrSges(4,2) + t279;
t375 = -t143 - t142;
t374 = t191 * t404 + t194 * t387;
t25 = mrSges(7,1) * t121 - mrSges(7,3) * t37;
t26 = -mrSges(7,2) * t121 + mrSges(7,3) * t38;
t373 = -t190 * t25 + t193 * t26;
t298 = qJD(3) * t191;
t22 = -t131 * t298 + t191 * t289 + t194 * t399;
t372 = -t191 * t23 + t194 * t22;
t20 = qJDD(3) * qJ(4) + qJD(3) * qJD(4) + t22;
t21 = -qJDD(3) * pkin(3) + t206;
t371 = t191 * t21 + t194 * t20;
t370 = t1 * t193 - t190 * t2;
t170 = Ifges(4,4) * t300;
t334 = Ifges(5,5) * t194;
t237 = Ifges(5,1) * t191 - t334;
t368 = Ifges(4,1) * t302 - t124 * Ifges(7,5) + t123 * Ifges(7,6) + t161 * Ifges(7,3) + qJD(2) * t237 + qJD(3) * t387 + t170;
t169 = Ifges(5,5) * t302;
t340 = Ifges(6,4) * t191;
t236 = -Ifges(6,1) * t194 - t340;
t116 = Ifges(7,4) * t123;
t33 = -Ifges(7,1) * t124 + Ifges(7,5) * t161 + t116;
t367 = -Ifges(5,3) * t300 + qJD(2) * t236 + t193 * t33 + t169 + (-Ifges(6,5) + Ifges(5,6)) * qJD(3);
t278 = mrSges(6,3) * t302;
t138 = qJD(3) * mrSges(6,2) - t278;
t366 = -t138 + t376;
t339 = Ifges(6,4) * t194;
t365 = t191 * t340 + (t334 - t339 + (-Ifges(6,2) + t405) * t191) * t194;
t189 = qJ(4) + pkin(5);
t239 = mrSges(7,1) * t193 - mrSges(7,2) * t190;
t364 = -m(7) * t189 - mrSges(6,1) + mrSges(4,2) - mrSges(5,3) - t239;
t363 = t343 * t393 + t397;
t362 = -mrSges(6,2) + t396 + t406;
t361 = -mrSges(3,1) - (mrSges(7,3) + t396) * t194 - (m(7) * pkin(5) + t239) * t191 - t240 - t242 - t241;
t360 = -m(7) * t39 + t283;
t320 = t185 * t192;
t276 = t191 * t320;
t104 = -t187 * t194 + t276;
t100 = t187 * t317 + t254;
t319 = t185 * t194;
t48 = t100 * t191 + t186 * t319;
t102 = -t187 * t255 + t316;
t256 = t185 * t322;
t50 = t102 * t191 - t194 * t256;
t211 = -g(1) * t50 - g(2) * t48 - g(3) * t104;
t197 = qJD(2) ^ 2;
t358 = Ifges(7,1) * t395 + Ifges(7,4) * t394 + Ifges(7,5) * t390;
t196 = -pkin(3) - pkin(4);
t355 = -t124 / 0.2e1;
t354 = t124 / 0.2e1;
t348 = pkin(8) * t194;
t345 = -qJD(2) / 0.2e1;
t344 = -qJD(6) / 0.2e1;
t341 = Ifges(4,4) * t194;
t337 = Ifges(7,4) * t190;
t336 = Ifges(7,4) * t193;
t335 = Ifges(5,5) * t191;
t318 = t185 * t195;
t315 = t190 * t191;
t314 = t190 * t194;
t313 = t191 * t193;
t311 = t193 * t194;
t307 = pkin(2) * t318 + pkin(8) * t320;
t301 = qJD(2) * t192;
t295 = qJD(6) * t190;
t294 = qJD(6) * t193;
t293 = qJD(6) * t194;
t286 = Ifges(7,5) * t37 + Ifges(7,6) * t38 + Ifges(7,3) * t121;
t285 = pkin(8) * t298;
t275 = t194 * t318;
t274 = t196 * qJD(3);
t272 = t185 * t301;
t271 = t195 * t303;
t84 = t99 * pkin(2);
t269 = pkin(8) * t100 + t84;
t85 = t101 * pkin(2);
t268 = pkin(8) * t102 + t85;
t43 = t48 * pkin(3);
t49 = -t185 * t186 * t191 + t100 * t194;
t258 = qJ(4) * t49 - t43;
t45 = t50 * pkin(3);
t51 = t102 * t194 + t191 * t256;
t257 = qJ(4) * t51 - t45;
t253 = -t292 / 0.2e1;
t252 = t292 / 0.2e1;
t59 = mrSges(6,1) * t130 + mrSges(6,2) * t129;
t105 = t187 * t191 + t192 * t319;
t95 = t104 * pkin(3);
t249 = qJ(4) * t105 - t95;
t248 = qJ(4) * t185 * t312 + pkin(3) * t275 + t307;
t235 = Ifges(7,1) * t193 - t337;
t234 = Ifges(7,1) * t190 + t336;
t231 = -t191 * Ifges(6,2) - t339;
t230 = -Ifges(7,2) * t190 + t336;
t229 = Ifges(7,2) * t193 + t337;
t227 = Ifges(6,5) * t191 - Ifges(6,6) * t194;
t226 = Ifges(7,5) * t193 - Ifges(7,6) * t190;
t225 = Ifges(7,5) * t190 + Ifges(7,6) * t193;
t224 = t10 * t193 + t11 * t190;
t68 = -mrSges(7,2) * t161 + mrSges(7,3) * t123;
t69 = mrSges(7,1) * t161 + mrSges(7,3) * t124;
t223 = -t190 * t69 + t193 * t68;
t222 = -t190 * t68 - t193 * t69;
t36 = t274 + t245;
t221 = t191 * t36 - t194 * t41;
t57 = -t104 * t190 + t185 * t310;
t58 = t104 * t193 + t190 * t318;
t219 = t46 * (mrSges(6,1) * t194 + mrSges(6,2) * t191);
t218 = t73 * (mrSges(5,1) * t191 - mrSges(5,3) * t194);
t217 = t132 * (mrSges(4,1) * t191 + mrSges(4,2) * t194);
t216 = t191 * (Ifges(4,1) * t194 - t342);
t209 = t190 * t293 + t191 * t297;
t208 = t190 * t298 - t193 * t293;
t207 = -g(1) * t101 - g(2) * t99 - g(3) * t318;
t203 = Ifges(7,5) * t194 + t191 * t235;
t202 = Ifges(7,6) * t194 + t191 * t230;
t201 = Ifges(7,3) * t194 + t191 * t226;
t200 = -qJD(6) * t224 + t370;
t13 = -qJ(5) * t129 + t194 * t291 - t20;
t166 = qJ(4) * t300;
t145 = -qJ(5) * t194 + t348;
t127 = pkin(3) * t302 - t166;
t107 = -Ifges(6,6) * qJD(3) + qJD(2) * t231;
t97 = pkin(3) * t298 - t306;
t96 = -qJ(5) * t298 + qJD(5) * t194 + t285;
t94 = t104 * pkin(4);
t90 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t130;
t83 = t196 * t302 + t166;
t74 = t191 * t274 + t306;
t62 = -qJD(3) * pkin(3) + t270;
t61 = mrSges(4,1) * t129 + mrSges(4,2) * t130;
t60 = mrSges(5,1) * t129 - mrSges(5,3) * t130;
t55 = -t164 + t309;
t54 = -qJD(3) * t276 + (t271 + t299) * t194;
t53 = qJD(3) * t105 + t191 * t271;
t52 = qJD(2) * t212 + t166;
t44 = t50 * pkin(4);
t42 = t48 * pkin(4);
t19 = t190 * t52 + t193 * t56;
t18 = -t190 * t56 + t193 * t52;
t17 = -pkin(4) * t129 + t210;
t16 = qJD(6) * t57 - t190 * t272 + t193 * t53;
t15 = -qJD(6) * t58 - t190 * t53 - t193 * t272;
t12 = qJDD(3) * t196 + t198;
t9 = qJDD(3) * pkin(5) - t13;
t6 = Ifges(7,4) * t37 + Ifges(7,2) * t38 + Ifges(7,6) * t121;
t3 = [m(2) * qJDD(1) + t15 * t69 + t16 * t68 + t57 * t25 + t58 * t26 - t366 * t53 + (-t90 + t382) * t104 + (t142 - t283) * t54 + (t381 + t385) * t105 + (-m(2) - m(3) - m(4) - m(5) + t393) * g(3) + ((mrSges(3,1) * qJDD(2) - mrSges(3,2) * t197 + t59 - t60 - t61) * t195 + (-mrSges(3,1) * t197 - mrSges(3,2) * qJDD(2) - qJD(2) * t398) * t192) * t185 + m(3) * (qJDD(1) * t187 ^ 2 + (t192 * t87 + t195 * t86) * t185) + m(4) * (-t104 * t23 + t105 * t22 + t53 * t309 + t54 * t71 + (t132 * t301 - t195 * t75) * t185) + m(5) * (t104 * t21 + t105 * t20 + t53 * t62 + t54 * t63 + (-t195 * t24 + t301 * t73) * t185) + m(6) * (t104 * t12 - t105 * t13 + t36 * t53 - t41 * t54 + (t17 * t195 - t301 * t46) * t185) + m(7) * (t1 * t58 + t10 * t15 + t105 * t9 + t11 * t16 + t2 * t57 + t39 * t54); t98 * t138 + t144 * t92 + m(6) * (t115 * t17 + t12 * t144 + t36 * t98 + t46 * t74) + t133 * t60 - t97 * t125 + t74 * t128 + t115 * t59 - pkin(2) * t61 + t28 * t25 + t29 * t26 + (t366 * t191 + (-t142 + t360) * t194) * t157 + t17 * t240 + (-m(5) * (t268 + t383) - m(4) * t268 + t393 * (pkin(4) * t321 + t383 + t85) + t363 * t102 + t361 * t101) * g(1) + (-m(5) * (t269 + t384) - m(4) * t269 + t393 * (pkin(4) * t325 + t384 + t84) + t361 * t99 + t363 * t100) * g(2) + t388 * t68 + t389 * t69 + (t1 * t29 + t10 * t389 + t11 * t388 + t145 * t9 + t2 * t28 - t39 * t96) * m(7) + (t378 + t391) * t96 + m(5) * (t133 * t24 + t73 * t97 + ((-t191 * t63 + t194 * t62) * qJD(3) + t371) * pkin(8)) + t375 * t285 + ((Ifges(4,6) / 0.2e1 + t400 / 0.2e1 + Ifges(6,5) / 0.2e1 - Ifges(5,6) / 0.2e1) * t194 + t403 * t191) * qJDD(3) + (-t87 + t150) * mrSges(3,2) + t39 * (mrSges(7,1) * t208 + mrSges(7,2) * t209) + (-t90 + t91) * pkin(8) * t191 - t9 * t238 * t194 - t24 * t241 - t75 * t242 + (t1 * t314 - t10 * t209 - t11 * t208 + t2 * t311) * mrSges(7,3) + qJD(3) * t217 + qJD(3) * t218 + qJD(3) * t219 + t191 * t401 + t121 * (Ifges(7,3) * t191 - t194 * t226) / 0.2e1 + t38 * (Ifges(7,6) * t191 - t194 * t230) / 0.2e1 - t130 * t231 / 0.2e1 - t129 * t233 / 0.2e1 + t37 * (Ifges(7,5) * t191 - t194 * t235) / 0.2e1 + (qJD(3) * t203 + t234 * t293) * t355 + t311 * t358 + t365 * t253 + (-m(6) * t13 + t385) * t145 + t381 * t348 + t6 * t314 / 0.2e1 + Ifges(3,3) * qJDD(2) + t123 * (qJD(3) * t202 + t229 * t293) / 0.2e1 + t161 * (qJD(3) * t201 + t225 * t293) / 0.2e1 + (t286 + (Ifges(4,1) + Ifges(5,1)) * t130 + (-Ifges(4,4) + Ifges(5,5)) * t129) * t191 / 0.2e1 - t191 * t402 + t398 * t273 + (-m(4) * t307 - m(5) * t248 + t393 * (pkin(4) * t275 + t248) + (t361 * t195 + (-qJ(5) * t393 + t397) * t192) * t185) * g(3) + (t190 * t33 + t193 * t32) * t293 / 0.2e1 + (t191 * Ifges(4,1) + t237 + t341) * t130 / 0.2e1 + (-t194 * Ifges(5,3) + t236 + t335) * t129 / 0.2e1 - t191 * (Ifges(6,4) * t129 - Ifges(6,2) * t130) / 0.2e1 + (-t63 * mrSges(5,2) - t71 * mrSges(4,3) - t377 / 0.2e1 - t41 * mrSges(6,3) + t367 / 0.2e1) * t298 + (t216 + t191 * (Ifges(5,1) * t194 + t335) + t194 * (-Ifges(4,2) * t191 + t341)) * t252 + (-t227 / 0.2e1 + t374 / 0.2e1) * qJD(3) ^ 2 + m(4) * (-pkin(2) * t75 + ((-t191 * t71 + t194 * t309) * qJD(3) + t372) * pkin(8)) + (-m(4) * (t132 * t192 + (t191 * t309 + t194 * t71) * t195) - m(5) * (t192 * t73 + (t191 * t62 + t194 * t63) * t195) - m(6) * (-t192 * t46 + t195 * t221)) * t305 + (t62 * mrSges(5,2) + t309 * mrSges(4,3) - t36 * mrSges(6,3) - t376 * pkin(8) + t368 / 0.2e1 + t10 * mrSges(7,1) - t107 / 0.2e1 - t11 * mrSges(7,2)) * t296 + t194 * (Ifges(4,4) * t130 - Ifges(4,2) * t129) / 0.2e1 - ((-Ifges(6,4) + Ifges(5,5)) * t130 + t405 * t129) * t194 / 0.2e1 + (t86 + t149) * mrSges(3,1) + t371 * mrSges(5,2) + t372 * mrSges(4,3) + (-t12 * t191 + t13 * t194) * mrSges(6,3); t227 * t252 + t196 * t92 - t193 * t6 / 0.2e1 + t189 * t14 - t56 * t138 + t127 * t125 - t83 * t128 - pkin(3) * t91 - t19 * t68 - t18 * t69 + t20 * mrSges(5,3) - t21 * mrSges(5,1) - t22 * mrSges(4,2) + t23 * mrSges(4,1) + t12 * mrSges(6,2) - t13 * mrSges(6,1) + t9 * t239 + (-t360 - t391 + t392) * qJD(4) + (m(7) * t200 - t294 * t69 - t295 * t68 + t373) * t181 + t374 * t253 + (-m(5) * t62 + t280 + t376) * t71 + t377 * t302 / 0.2e1 + t403 * t130 + (-m(7) * (-t42 - t43) - m(6) * (t258 - t42) - m(5) * t258 + t362 * t48 + t364 * t49) * g(2) + (-m(7) * (-t44 - t45) - m(6) * (t257 - t44) - m(5) * t257 + t362 * t50 + t364 * t51) * g(1) + (-m(7) * (-t94 - t95) - m(6) * (t249 - t94) - m(5) * t249 + t362 * t104 + t364 * t105) * g(3) + (t93 - t88) * qJ(4) - t302 * t380 + (t235 * t354 - t380) * qJD(6) + (-t10 * t18 - t11 * t19 + t189 * t9 + t39 * t55) * m(7) + (-pkin(3) * t21 + qJ(4) * t20 - t127 * t73) * m(5) + t190 * t358 + t225 * t390 + t229 * t394 + t234 * t395 - (Ifges(5,1) * t300 + t169 + t367) * t302 / 0.2e1 - (-Ifges(4,2) * t302 + t170 + t368) * t300 / 0.2e1 - t378 * t55 + (t10 * t294 + t11 * t295 - t211 - t370) * mrSges(7,3) + (-qJ(4) * t13 + t12 * t196 - t36 * t56 - t41 * t55 - t46 * t83) * m(6) + (t202 * t345 + t230 * t344) * t123 + t107 * t300 / 0.2e1 + t32 * t295 / 0.2e1 - t33 * t294 / 0.2e1 - t62 * t281 + (-t216 / 0.2e1 + t365 / 0.2e1) * t197 - t400 * t129 + (t201 * t345 + t226 * t344) * t161 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) - (t279 + t375 - t392) * t309 + t36 * t277 + t41 * t278 + t63 * t282 + (t203 * t354 - t10 * (mrSges(7,1) * t194 - mrSges(7,3) * t313) - t11 * (-mrSges(7,2) * t194 - mrSges(7,3) * t315) - t217 - t218 - t219) * qJD(2); t222 * qJD(6) + t283 * qJD(3) + (t222 + t308) * t302 + (-qJD(3) * t39 - t224 * t302 + t200 + t211) * m(7) + (qJD(3) * t41 - t302 * t46 + t12 + t211) * m(6) + (-qJD(3) * t63 + t302 * t73 + t21 + t211) * m(5) + t373 + t382; t190 * t26 + t193 * t25 + t223 * qJD(6) + (t1 * t190 + t193 * t2 + (-t10 * t190 + t11 * t193) * qJD(6) + t207) * m(7) + (t17 + t207) * m(6) + (-t378 * t194 + (t138 + t223) * t191 + m(6) * t221 - m(7) * (t10 * t315 - t11 * t313 - t194 * t39)) * qJD(2) + t59; -t402 + t401 - t39 * (-mrSges(7,1) * t124 + mrSges(7,2) * t123) + (Ifges(7,1) * t123 + t338) * t354 + t32 * t355 - t161 * (Ifges(7,5) * t123 + Ifges(7,6) * t124) / 0.2e1 - t10 * t68 + t11 * t69 - g(1) * ((t101 * t193 - t190 * t50) * mrSges(7,1) + (-t101 * t190 - t193 * t50) * mrSges(7,2)) - g(2) * ((-t190 * t48 + t193 * t99) * mrSges(7,1) + (-t190 * t99 - t193 * t48) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t57 - mrSges(7,2) * t58) + (t10 * t123 - t11 * t124) * mrSges(7,3) + t286 - (Ifges(7,2) * t124 + t116 + t33) * t123 / 0.2e1;];
tau  = t3;
