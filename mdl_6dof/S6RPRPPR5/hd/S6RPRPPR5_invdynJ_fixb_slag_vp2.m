% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:16
% EndTime: 2019-03-09 02:50:39
% DurationCPUTime: 19.48s
% Computational Cost: add. (8108->654), mult. (19181->821), div. (0->0), fcn. (14268->14), ass. (0->284)
t411 = -mrSges(5,2) + mrSges(4,1);
t412 = -mrSges(5,1) - mrSges(4,3);
t401 = Ifges(6,3) + Ifges(4,1);
t226 = sin(pkin(10));
t228 = cos(pkin(10));
t233 = sin(qJ(6));
t236 = cos(qJ(6));
t377 = -t226 * t233 + t228 * t236;
t171 = t377 * qJD(6);
t227 = sin(pkin(9));
t229 = cos(pkin(9));
t234 = sin(qJ(3));
t346 = cos(qJ(3));
t180 = t227 * t346 + t234 * t229;
t170 = t180 * qJD(1);
t382 = t377 * t170;
t380 = t171 + t382;
t293 = t227 ^ 2 + t229 ^ 2;
t289 = qJD(1) * qJD(2);
t198 = qJ(2) * qJDD(1) + t289;
t388 = Ifges(5,5) - Ifges(4,6);
t267 = -mrSges(3,1) * t229 + mrSges(3,2) * t227;
t410 = m(3) * pkin(1) + mrSges(2,1) - t267;
t220 = pkin(10) + qJ(6);
t214 = sin(t220);
t216 = cos(t220);
t263 = mrSges(6,1) * t226 + mrSges(6,2) * t228;
t343 = pkin(5) * t226;
t409 = -m(7) * t343 - mrSges(7,1) * t214 - mrSges(7,2) * t216 - t263;
t256 = t236 * t226 + t233 * t228;
t172 = t256 * qJD(6);
t248 = t256 * t170;
t378 = -t248 - t172;
t221 = pkin(9) + qJ(3);
t215 = sin(t221);
t217 = cos(t221);
t408 = -t411 * t217 + (mrSges(4,2) - mrSges(5,3)) * t215;
t210 = pkin(5) * t228 + pkin(4);
t264 = -t228 * mrSges(6,1) + t226 * mrSges(6,2);
t279 = m(3) * qJ(2) + mrSges(3,3);
t331 = pkin(7) + qJ(2);
t407 = t264 - m(6) * (pkin(4) + t331) - m(7) * t210 + mrSges(2,2) - t279 + t412;
t230 = -pkin(8) - qJ(5);
t333 = pkin(3) + qJ(5);
t406 = -m(7) * (-pkin(3) + t230) + mrSges(7,3) + m(6) * t333 + mrSges(6,3);
t283 = t346 * t229;
t272 = qJD(1) * t283;
t297 = t227 * t234;
t169 = qJD(1) * t297 - t272;
t211 = pkin(2) * t229 + pkin(1);
t186 = -qJD(1) * t211 + qJD(2);
t243 = -qJ(4) * t170 + t186;
t102 = pkin(3) * t169 + t243;
t187 = t331 * t227;
t181 = qJD(1) * t187;
t188 = t331 * t229;
t182 = qJD(1) * t188;
t137 = t346 * t181 + t182 * t234;
t127 = -qJD(3) * pkin(3) + qJD(4) + t137;
t145 = -qJD(3) * t226 + t169 * t228;
t386 = t145 * Ifges(6,6);
t144 = -t228 * qJD(3) - t226 * t169;
t387 = t144 * Ifges(6,5);
t77 = t169 * t333 + t243;
t250 = pkin(4) * t170 + t137;
t84 = -qJD(3) * t333 + qJD(4) + t250;
t40 = -t226 * t77 + t228 * t84;
t41 = t226 * t84 + t228 * t77;
t22 = pkin(5) * t170 + pkin(8) * t144 + t40;
t25 = pkin(8) * t145 + t41;
t8 = t22 * t236 - t233 * t25;
t9 = t22 * t233 + t236 * t25;
t405 = t127 * mrSges(5,1) + t40 * mrSges(6,1) + t8 * mrSges(7,1) + t186 * mrSges(4,2) - t41 * mrSges(6,2) - t9 * mrSges(7,2) + t137 * mrSges(4,3) - t102 * mrSges(5,3) + t386 / 0.2e1 - t387 / 0.2e1;
t404 = -m(5) - m(7);
t138 = -t234 * t181 + t346 * t182;
t225 = qJD(3) * qJ(4);
t128 = -t225 - t138;
t403 = m(5) * t128;
t235 = sin(qJ(1));
t402 = g(2) * t235;
t359 = -m(6) - m(7);
t291 = m(5) - t359;
t389 = Ifges(5,4) - Ifges(4,5);
t315 = qJDD(3) / 0.2e1;
t150 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t169;
t152 = mrSges(5,1) * t169 - qJD(3) * mrSges(5,3);
t400 = t150 - t152;
t381 = t411 * qJD(3) + t170 * t412;
t399 = -m(4) + t404;
t161 = Ifges(4,4) * t169;
t165 = qJD(6) + t170;
t274 = t144 * t233 + t236 * t145;
t85 = t144 * t236 - t145 * t233;
t398 = Ifges(4,5) * qJD(3) - Ifges(7,5) * t85 + Ifges(7,6) * t274 + t165 * Ifges(7,3) + t170 * t401 - t161 + t386 - t387;
t397 = 0.2e1 * t315;
t174 = t180 * qJD(3);
t277 = qJDD(1) * t346;
t286 = qJDD(1) * t234;
t136 = qJD(1) * t174 + t227 * t286 - t229 * t277;
t114 = -qJDD(3) * t226 + t136 * t228;
t292 = qJD(3) * t234;
t280 = t227 * t292;
t135 = qJD(1) * t280 - qJD(3) * t272 - t227 * t277 - t229 * t286;
t185 = -qJDD(1) * t211 + qJDD(2);
t241 = qJ(4) * t135 - qJD(4) * t170 + t185;
t29 = qJD(5) * t169 + t136 * t333 + t241;
t275 = pkin(7) * qJDD(1) + t198;
t158 = t275 * t227;
t159 = t275 * t229;
t278 = qJD(3) * t346;
t69 = -t158 * t346 - t234 * t159 + t181 * t292 - t182 * t278;
t244 = qJDD(4) - t69;
t46 = -t135 * pkin(4) - qJD(3) * qJD(5) - qJDD(3) * t333 + t244;
t15 = t226 * t46 + t228 * t29;
t10 = pkin(8) * t114 + t15;
t115 = qJDD(3) * t228 + t136 * t226;
t14 = -t226 * t29 + t228 * t46;
t5 = -pkin(5) * t135 - pkin(8) * t115 + t14;
t1 = qJD(6) * t8 + t10 * t236 + t233 * t5;
t2 = -qJD(6) * t9 - t10 * t233 + t236 * t5;
t395 = -t1 * t256 - t2 * t377 - t378 * t8 - t380 * t9;
t394 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t294 = t217 * pkin(3) + t215 * qJ(4);
t393 = -m(5) * t294 + t408;
t328 = Ifges(6,4) * t226;
t260 = Ifges(6,2) * t228 + t328;
t327 = Ifges(6,4) * t228;
t261 = Ifges(6,1) * t226 + t327;
t99 = -pkin(4) * t169 + t138;
t90 = qJD(5) + t225 + t99;
t391 = t138 * mrSges(4,3) + t102 * mrSges(5,2) - t90 * t264 - t128 * mrSges(5,1) + t144 * t261 / 0.2e1 - t145 * t260 / 0.2e1 - t186 * mrSges(4,1);
t32 = qJD(6) * t274 + t114 * t233 + t115 * t236;
t368 = t32 / 0.2e1;
t33 = qJD(6) * t85 + t114 * t236 - t115 * t233;
t367 = t33 / 0.2e1;
t358 = t114 / 0.2e1;
t357 = t115 / 0.2e1;
t133 = qJDD(6) - t135;
t356 = t133 / 0.2e1;
t355 = -t135 / 0.2e1;
t330 = -pkin(8) - t333;
t183 = t330 * t226;
t184 = t330 * t228;
t140 = t183 * t236 + t184 * t233;
t94 = t228 * t99;
t314 = qJ(4) * t169;
t95 = t170 * t333 + t314;
t34 = -pkin(5) * t169 + t94 + (-pkin(8) * t170 - t95) * t226;
t309 = t170 * t228;
t51 = t226 * t99 + t228 * t95;
t39 = pkin(8) * t309 + t51;
t385 = -qJD(5) * t377 - qJD(6) * t140 + t233 * t39 - t236 * t34;
t139 = -t183 * t233 + t184 * t236;
t384 = -qJD(5) * t256 + qJD(6) * t139 - t233 * t34 - t236 * t39;
t237 = cos(qJ(1));
t374 = g(1) * t237 + t402;
t383 = t217 * t374;
t299 = t217 * t237;
t303 = t215 * t237;
t379 = pkin(3) * t299 + qJ(4) * t303;
t108 = -mrSges(6,2) * t170 + mrSges(6,3) * t145;
t109 = mrSges(6,1) * t170 + mrSges(6,3) * t144;
t376 = -t108 * t226 - t109 * t228;
t63 = mrSges(6,2) * t135 + mrSges(6,3) * t114;
t64 = -mrSges(6,1) * t135 - mrSges(6,3) * t115;
t375 = t226 * t63 + t228 * t64;
t258 = t14 * t228 + t15 * t226;
t48 = -mrSges(7,1) * t274 - mrSges(7,2) * t85;
t91 = -mrSges(6,1) * t145 - mrSges(6,2) * t144;
t284 = -t48 - t91 + t152;
t373 = (-mrSges(5,3) + t409) * t217 + (m(5) * pkin(3) - mrSges(5,2) + t406) * t215;
t106 = t234 * (qJD(2) * t227 + qJD(3) * t188) - qJD(2) * t283 + t187 * t278;
t372 = qJ(4) * t291;
t370 = Ifges(7,4) * t368 + Ifges(7,2) * t367 + Ifges(7,6) * t356;
t369 = Ifges(7,1) * t368 + Ifges(7,4) * t367 + Ifges(7,5) * t356;
t83 = Ifges(7,4) * t274;
t37 = -Ifges(7,1) * t85 + Ifges(7,5) * t165 + t83;
t365 = t37 / 0.2e1;
t364 = Ifges(6,1) * t357 + Ifges(6,4) * t358 + Ifges(6,5) * t355;
t363 = -t274 / 0.2e1;
t362 = t274 / 0.2e1;
t361 = t85 / 0.2e1;
t360 = -t85 / 0.2e1;
t354 = -t165 / 0.2e1;
t353 = t165 / 0.2e1;
t352 = -t169 / 0.2e1;
t351 = t169 / 0.2e1;
t350 = -t170 / 0.2e1;
t349 = t170 / 0.2e1;
t344 = Ifges(7,4) * t85;
t340 = g(3) * t217;
t335 = -qJD(3) / 0.2e1;
t334 = qJD(3) / 0.2e1;
t178 = -t283 + t297;
t173 = -t229 * t278 + t280;
t254 = qJ(4) * t173 - qJD(4) * t180;
t60 = qJD(5) * t178 + t174 * t333 + t254;
t143 = -t234 * t187 + t188 * t346;
t107 = qJD(2) * t180 + qJD(3) * t143;
t79 = -t173 * pkin(4) + t107;
t31 = t226 * t79 + t228 * t60;
t320 = t170 * Ifges(4,4);
t319 = t170 * Ifges(5,6);
t318 = t217 * mrSges(7,3);
t142 = t346 * t187 + t188 * t234;
t116 = pkin(4) * t180 + t142;
t252 = -qJ(4) * t180 - t211;
t97 = t178 * t333 + t252;
t54 = t226 * t116 + t228 * t97;
t310 = t170 * t226;
t308 = t174 * t226;
t307 = t174 * t228;
t306 = t178 * t226;
t305 = t178 * t228;
t304 = t214 * t235;
t302 = t216 * t235;
t301 = t216 * t237;
t300 = t217 * t230;
t288 = qJDD(1) * t227;
t287 = qJDD(1) * t229;
t285 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t133;
t11 = -t33 * mrSges(7,1) + t32 * mrSges(7,2);
t204 = qJ(4) + t343;
t119 = -t135 * mrSges(5,1) + qJDD(3) * mrSges(5,2);
t55 = -t114 * mrSges(6,1) + t115 * mrSges(6,2);
t192 = t237 * t211;
t273 = t235 * t331 + t192;
t268 = -mrSges(3,1) * t287 + mrSges(3,2) * t288;
t259 = Ifges(6,5) * t226 + Ifges(6,6) * t228;
t257 = t226 * t41 + t228 * t40;
t104 = t228 * t116;
t38 = pkin(5) * t180 + t104 + (-pkin(8) * t178 - t97) * t226;
t49 = pkin(8) * t305 + t54;
t16 = -t233 * t49 + t236 * t38;
t17 = t233 * t38 + t236 * t49;
t120 = t377 * t178;
t68 = -t234 * t158 + t346 * t159 - t181 * t278 - t182 * t292;
t59 = -qJDD(3) * qJ(4) - qJD(3) * qJD(4) - t68;
t47 = -pkin(4) * t136 + qJDD(5) - t59;
t212 = -qJDD(1) * pkin(1) + qJDD(2);
t160 = Ifges(5,6) * t169;
t149 = -t215 * t304 + t301;
t148 = t214 * t237 + t215 * t302;
t147 = t214 * t303 + t302;
t146 = t215 * t301 - t304;
t134 = pkin(3) * t178 + t252;
t132 = t135 * mrSges(4,2);
t131 = t135 * mrSges(5,3);
t130 = -mrSges(5,2) * t169 - mrSges(5,3) * t170;
t129 = pkin(3) * t170 + t314;
t124 = -t169 * Ifges(4,2) + Ifges(4,6) * qJD(3) + t320;
t123 = Ifges(5,4) * qJD(3) - t170 * Ifges(5,2) + t160;
t122 = Ifges(5,5) * qJD(3) + t169 * Ifges(5,3) - t319;
t121 = t256 * t178;
t118 = mrSges(5,1) * t136 - qJDD(3) * mrSges(5,3);
t117 = -t178 * pkin(4) + t143;
t96 = pkin(3) * t174 + t254;
t82 = -t178 * t210 + t143;
t78 = -pkin(4) * t174 - t106;
t76 = t228 * t79;
t73 = -t170 * t210 - t137;
t72 = -Ifges(6,1) * t144 + Ifges(6,4) * t145 + Ifges(6,5) * t170;
t71 = -Ifges(6,4) * t144 + Ifges(6,2) * t145 + Ifges(6,6) * t170;
t67 = mrSges(7,1) * t165 + mrSges(7,3) * t85;
t66 = -mrSges(7,2) * t165 + mrSges(7,3) * t274;
t65 = -pkin(5) * t145 + t90;
t62 = -qJDD(3) * pkin(3) + t244;
t61 = -t174 * t210 - t106;
t58 = -t172 * t178 + t174 * t377;
t57 = qJD(6) * t120 + t174 * t256;
t53 = -t226 * t97 + t104;
t52 = pkin(3) * t136 + t241;
t50 = -t226 * t95 + t94;
t44 = t115 * Ifges(6,4) + t114 * Ifges(6,2) - t135 * Ifges(6,6);
t36 = Ifges(7,2) * t274 + Ifges(7,6) * t165 - t344;
t30 = -t226 * t60 + t76;
t23 = -pkin(5) * t114 + t47;
t21 = pkin(8) * t307 + t31;
t20 = -pkin(5) * t173 + t76 + (-pkin(8) * t174 - t60) * t226;
t19 = -mrSges(7,2) * t133 + mrSges(7,3) * t33;
t18 = mrSges(7,1) * t133 - mrSges(7,3) * t32;
t4 = -qJD(6) * t17 + t20 * t236 - t21 * t233;
t3 = qJD(6) * t16 + t20 * t233 + t21 * t236;
t6 = [(t185 * mrSges(4,1) + t59 * mrSges(5,1) - t52 * mrSges(5,2) - t68 * mrSges(4,3) + Ifges(5,6) * t135 + t259 * t355 + t260 * t358 + t261 * t357 + t47 * t264 + (Ifges(5,3) + Ifges(4,2)) * t136 + (-t355 + t135 / 0.2e1) * Ifges(4,4) + t388 * t397) * t178 + (-t149 * mrSges(7,1) + t148 * mrSges(7,2) + (t331 * t399 + t407) * t237 + (t406 * t217 + (m(6) * qJ(4) + m(7) * t204 + t263) * t215 + (m(6) - t399) * t211 - t393 + t410) * t235) * g(1) + (-m(6) * (qJ(5) * t299 + t192 + t379) - mrSges(6,3) * t299 - t147 * mrSges(7,1) - t146 * mrSges(7,2) - m(4) * t273 - t263 * t303 + t404 * (t273 + t379) + t407 * t235 + (-m(7) * (t215 * t343 - t300) - t318 + t408 - t410) * t237) * g(2) + (t1 * t120 - t121 * t2 - t57 * t8 + t58 * t9) * mrSges(7,3) + (-m(4) * t69 + m(5) * t62 - qJDD(3) * mrSges(4,1) - mrSges(4,3) * t135 + t119) * t142 + (-t14 * t306 + t15 * t305 + t307 * t41 - t308 * t40) * mrSges(6,3) + (m(4) * t68 - m(5) * t59 - qJDD(3) * mrSges(4,2) - mrSges(4,3) * t136 - t118) * t143 + (Ifges(7,4) * t57 + Ifges(7,2) * t58) * t362 + (Ifges(7,4) * t121 + Ifges(7,2) * t120) * t367 + m(5) * (t102 * t96 + t134 * t52) + m(6) * (t117 * t47 + t14 * t53 + t15 * t54 + t30 * t40 + t31 * t41 + t78 * t90) + m(7) * (t1 * t17 + t16 * t2 + t23 * t82 + t3 * t9 + t4 * t8 + t61 * t65) + 0.2e1 * t293 * t198 * mrSges(3,3) + t82 * t11 + t53 * t64 + t65 * (-mrSges(7,1) * t58 + mrSges(7,2) * t57) + t3 * t66 + t4 * t67 + t58 * t36 / 0.2e1 + t61 * t48 + t54 * t63 + t16 * t18 + t17 * t19 + (-m(4) * t185 - t136 * mrSges(4,1) + t132) * t211 + (Ifges(3,4) * t227 + Ifges(3,2) * t229) * t287 + (Ifges(3,1) * t227 + Ifges(3,4) * t229) * t288 + (Ifges(7,1) * t57 + Ifges(7,4) * t58) * t360 + (Ifges(7,1) * t121 + Ifges(7,4) * t120) * t368 + Ifges(2,3) * qJDD(1) + t212 * t267 - pkin(1) * t268 + m(3) * (-pkin(1) * t212 + (t198 + t289) * qJ(2) * t293) + (-t401 * t349 + t123 / 0.2e1 + Ifges(5,2) * t350 + Ifges(5,6) * t351 - Ifges(4,4) * t352 - Ifges(7,3) * t353 - Ifges(7,5) * t360 - Ifges(7,6) * t362 + t389 * t334 - t405 - t398 / 0.2e1) * t173 + t72 * t308 / 0.2e1 + t44 * t305 / 0.2e1 + t71 * t307 / 0.2e1 + t78 * t91 + t31 * t108 + t30 * t109 + t117 * t55 + t23 * (-mrSges(7,1) * t120 + mrSges(7,2) * t121) + (m(4) * t137 + m(5) * t127 - t381) * t107 + t96 * t130 + t134 * (-t136 * mrSges(5,2) + t131) + t306 * t364 + t57 * t365 + t121 * t369 + t120 * t370 + (t62 * mrSges(5,1) + t14 * mrSges(6,1) + t185 * mrSges(4,2) - t15 * mrSges(6,2) - t69 * mrSges(4,3) - t52 * mrSges(5,3) + Ifges(4,5) * t315 + Ifges(6,5) * t357 + Ifges(7,5) * t368 - Ifges(5,2) * t135 + Ifges(6,6) * t358 + Ifges(7,6) * t367 + Ifges(7,3) * t356 + t401 * t355 + (-Ifges(4,4) / 0.2e1 - Ifges(5,6)) * t136 - t397 * Ifges(5,4) + t394) * t180 + (-Ifges(4,4) * t136 + Ifges(4,5) * qJDD(3) + Ifges(6,5) * t115 + Ifges(6,6) * t114 - t401 * t135 + t285) * t180 / 0.2e1 + (Ifges(5,6) * t350 + Ifges(5,3) * t351 - Ifges(4,2) * t352 + t122 / 0.2e1 - t124 / 0.2e1 + (t259 - Ifges(4,4)) * t349 + t388 * t334 - t391) * t174 + (-m(4) * t138 - t400 + t403) * t106 + (Ifges(7,5) * t57 + Ifges(7,6) * t58) * t353 + (Ifges(7,5) * t121 + Ifges(7,6) * t120) * t356; t411 * t136 + m(3) * t212 - t380 * t67 + t131 - t132 + t228 * t63 + (t376 + t381) * t170 - t226 * t64 + t268 + (t150 - t284) * t169 - t256 * t18 + t377 * t19 + t378 * t66 + (-g(1) * t235 + g(2) * t237) * (m(3) + m(4) + t291) - t279 * t293 * qJD(1) ^ 2 + (t1 * t377 + t169 * t65 - t2 * t256 + t378 * t9 - t380 * t8) * m(7) + (-t14 * t226 + t15 * t228 + t169 * t90 - t170 * t257) * m(6) + (-t127 * t170 - t128 * t169 + t52) * m(5) + (-t137 * t170 + t138 * t169 + t185) * m(4); (-m(7) * (t294 - t300) - t318 - m(6) * (qJ(5) * t217 + t294) - t217 * mrSges(6,3) + t409 * t215 + t393) * g(3) + (-Ifges(7,5) * t361 + Ifges(5,2) * t349 - Ifges(7,6) * t363 - Ifges(7,3) * t354 + t389 * t335 + t405) * t169 + (-t217 * t372 + t373) * t402 + (t237 * t373 - t299 * t372) * g(1) + t377 * t369 + (Ifges(7,5) * t377 - Ifges(7,6) * t256) * t356 + (Ifges(7,4) * t377 - Ifges(7,2) * t256) * t367 + (Ifges(7,1) * t377 - Ifges(7,4) * t256) * t368 + t23 * (mrSges(7,1) * t256 + mrSges(7,2) * t377) + (Ifges(5,3) * t352 + t388 * t335 + t391) * t170 + (-pkin(3) * t62 - qJ(4) * t59 - t102 * t129 - t127 * t138 - t128 * t137) * m(5) + (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + (Ifges(7,1) * t248 + Ifges(7,4) * t382) * t361 + (Ifges(7,5) * t248 + Ifges(7,6) * t382) * t354 + (Ifges(7,4) * t248 + Ifges(7,2) * t382) * t363 + (t160 + t123) * t352 + (t319 + t124) * t349 + t374 * (mrSges(4,1) * t215 + mrSges(4,2) * t217) + (-t309 * t41 + t310 * t40 - t258) * mrSges(6,3) + t376 * qJD(5) + (mrSges(7,1) * t380 + mrSges(7,2) * t378) * t65 - t380 * t36 / 0.2e1 - t68 * mrSges(4,2) + t69 * mrSges(4,1) - t73 * t48 - t59 * mrSges(5,3) + t62 * mrSges(5,2) + (t55 - t118) * qJ(4) + t47 * t263 + (-Ifges(7,4) * t172 - Ifges(7,2) * t171) * t362 + (-Ifges(7,5) * t172 - Ifges(7,6) * t171) * t353 + (-Ifges(7,1) * t172 - Ifges(7,4) * t171) * t360 - t248 * t37 / 0.2e1 - t226 * t44 / 0.2e1 - t71 * t309 / 0.2e1 - t72 * t310 / 0.2e1 - t256 * t370 - t375 * t333 + (qJ(4) * t47 - qJD(5) * t257 + t250 * t90 - t258 * t333 - t40 * t50 - t41 * t51) * m(6) + t381 * t138 + t395 * mrSges(7,3) - t51 * t108 - t50 * t109 + t384 * t66 + (t1 * t140 + t139 * t2 + t204 * t23 + t384 * t9 + t385 * t8 - t65 * t73) * m(7) + t385 * t67 + (-Ifges(4,2) * t170 - t161 + t398) * t351 - pkin(3) * t119 + t400 * t137 - t129 * t130 + t388 * t136 + t389 * t135 + t139 * t18 + t140 * t19 + (Ifges(6,5) * t228 - Ifges(6,6) * t226) * t355 + (Ifges(6,1) * t228 - t328) * t357 + (-Ifges(6,2) * t226 + t327) * t358 + t228 * t364 - t172 * t365 + (-t169 * t401 + t170 * t259 + t122 - t320) * t350 + (m(6) * t90 + m(7) * t65 - t284 - t403) * qJD(4) + t204 * t11 + t250 * t91; t256 * t19 + t377 * t18 + t378 * t67 + t380 * t66 + (t108 * t228 - t109 * t226 + t130) * t170 + t284 * qJD(3) + t119 + (-qJD(3) * t65 - t395) * m(7) + (-qJD(3) * t90 - (t226 * t40 - t228 * t41) * t170 + t258) * m(6) + (qJD(3) * t128 + t102 * t170 + t62) * m(5) + (-t215 * t374 + t340) * t291 + t375; t359 * t215 * g(3) - t145 * t108 - t144 * t109 - t274 * t66 - t85 * t67 + t11 + t55 + (-t274 * t9 - t8 * t85 + t23 - t383) * m(7) + (-t144 * t40 - t145 * t41 - t383 + t47) * m(6); -t65 * (-mrSges(7,1) * t85 + mrSges(7,2) * t274) + (Ifges(7,1) * t274 + t344) * t361 + t36 * t360 + (Ifges(7,5) * t274 + Ifges(7,6) * t85) * t354 - t8 * t66 + t9 * t67 - g(1) * (mrSges(7,1) * t146 - mrSges(7,2) * t147) - g(2) * (mrSges(7,1) * t148 + mrSges(7,2) * t149) - (-mrSges(7,1) * t216 + mrSges(7,2) * t214) * t340 + (t274 * t8 - t85 * t9) * mrSges(7,3) + t285 + (Ifges(7,2) * t85 + t37 + t83) * t363 + t394;];
tau  = t6;
