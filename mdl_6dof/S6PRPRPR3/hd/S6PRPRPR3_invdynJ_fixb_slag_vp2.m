% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:31
% EndTime: 2019-03-08 19:34:52
% DurationCPUTime: 12.63s
% Computational Cost: add. (3864->555), mult. (8813->751), div. (0->0), fcn. (6669->12), ass. (0->260)
t381 = -mrSges(5,1) + mrSges(6,2);
t374 = m(6) + m(7);
t373 = Ifges(5,5) - Ifges(6,4);
t372 = Ifges(6,5) - Ifges(5,6);
t173 = sin(qJ(4));
t176 = cos(qJ(4));
t210 = pkin(9) * t173 - qJ(5) * t176;
t271 = qJD(4) * t173;
t232 = pkin(4) * t271 - qJD(5) * t173;
t168 = sin(pkin(6));
t166 = sin(pkin(11));
t169 = cos(pkin(11));
t174 = sin(qJ(2));
t177 = cos(qJ(2));
t207 = t166 * t177 + t169 * t174;
t102 = t207 * t168;
t93 = qJD(1) * t102;
t380 = qJD(4) * t210 + t232 - t93;
t172 = sin(qJ(6));
t175 = cos(qJ(6));
t263 = qJD(2) * qJD(4);
t135 = qJDD(2) * t173 + t176 * t263;
t124 = qJDD(6) + t135;
t272 = qJD(4) * t172;
t273 = qJD(2) * t176;
t127 = -t175 * t273 - t272;
t134 = -t176 * qJDD(2) + t173 * t263;
t71 = qJD(6) * t127 + qJDD(4) * t175 + t134 * t172;
t35 = mrSges(7,1) * t124 - mrSges(7,3) * t71;
t264 = t173 * qJD(2);
t155 = qJD(6) + t264;
t86 = -mrSges(7,2) * t155 + mrSges(7,3) * t127;
t270 = qJD(4) * t175;
t196 = t172 * t273 - t270;
t87 = mrSges(7,1) * t155 + mrSges(7,3) * t196;
t355 = -t172 * t87 + t175 * t86;
t72 = qJD(6) * t196 - qJDD(4) * t172 + t134 * t175;
t36 = -mrSges(7,2) * t124 + mrSges(7,3) * t72;
t348 = qJD(6) * t355 + t172 * t36 + t175 * t35;
t274 = qJD(1) * t177;
t243 = qJD(2) * t274;
t106 = (qJDD(1) * t174 + t243) * t168;
t289 = t168 * t174;
t248 = qJD(1) * t289;
t231 = qJD(2) * t248;
t287 = t168 * t177;
t105 = qJDD(1) * t287 - t231;
t296 = qJDD(2) * pkin(2);
t99 = t105 + t296;
t47 = -t166 * t106 + t169 * t99;
t40 = -qJDD(2) * pkin(3) - t47;
t180 = -qJ(5) * t135 - qJD(5) * t264 + t40;
t337 = pkin(4) + pkin(9);
t13 = t134 * t337 + t180;
t171 = cos(pkin(6));
t151 = qJDD(1) * t171 + qJDD(3);
t153 = qJD(1) * t171 + qJD(3);
t269 = qJD(4) * t176;
t48 = t169 * t106 + t166 * t99;
t41 = qJDD(2) * pkin(8) + t48;
t247 = t168 * t274;
t136 = qJD(2) * pkin(2) + t247;
t83 = t166 * t136 + t169 * t248;
t79 = qJD(2) * pkin(8) + t83;
t12 = t151 * t176 - t153 * t271 - t173 * t41 - t79 * t269;
t194 = qJDD(5) - t12;
t5 = pkin(5) * t135 - qJDD(4) * t337 + t194;
t139 = t176 * t153;
t42 = -t173 * (pkin(5) * qJD(2) + t79) + t139;
t352 = -t42 + qJD(5);
t24 = -qJD(4) * t337 + t352;
t300 = qJ(5) * t173;
t237 = -pkin(3) - t300;
t192 = -t176 * t337 + t237;
t138 = t166 * t248;
t82 = t136 * t169 - t138;
t46 = qJD(2) * t192 - t82;
t9 = -t172 * t46 + t175 * t24;
t1 = qJD(6) * t9 + t13 * t175 + t172 * t5;
t379 = t1 * mrSges(7,2);
t10 = t172 * t24 + t175 * t46;
t2 = -qJD(6) * t10 - t13 * t172 + t175 * t5;
t378 = t2 * mrSges(7,1);
t252 = mrSges(5,3) * t264;
t254 = mrSges(6,1) * t264;
t277 = t381 * qJD(4) + t252 + t254;
t167 = sin(pkin(10));
t170 = cos(pkin(10));
t284 = t171 * t177;
t377 = -t167 * t174 + t170 * t284;
t376 = m(7) * pkin(9) + pkin(4) * t374 - t381;
t228 = t1 * t172 + t175 * t2;
t267 = qJD(6) * t175;
t268 = qJD(6) * t172;
t375 = -t10 * t267 + t9 * t268 - t228;
t158 = pkin(2) * t166 + pkin(8);
t324 = pkin(5) + t158;
t123 = t324 * t176;
t112 = qJD(4) * t123;
t281 = t173 * t175;
t331 = pkin(2) * t169;
t103 = t192 - t331;
t122 = t324 * t173;
t59 = t103 * t175 + t122 * t172;
t96 = t169 * t247 - t138;
t371 = -qJD(6) * t59 + t112 * t175 - t172 * t380 - t281 * t96;
t283 = t172 * t173;
t58 = -t103 * t172 + t122 * t175;
t370 = qJD(6) * t58 + t112 * t172 + t175 * t380 - t283 * t96;
t279 = t177 * t169;
t125 = t166 * t174 - t279;
t193 = t125 * t171;
t61 = -t167 * t207 - t170 * t193;
t306 = t176 * t61;
t369 = pkin(4) * t306 + t61 * t300;
t64 = t167 * t193 - t170 * t207;
t305 = t176 * t64;
t368 = pkin(4) * t305 + t64 * t300;
t54 = t153 * t173 + t176 * t79;
t45 = -qJD(4) * qJ(5) - t54;
t8 = -qJDD(4) * pkin(4) + t194;
t367 = -qJD(4) * t45 - t8;
t101 = t166 * t289 - t168 * t279;
t295 = t101 * t176;
t366 = -pkin(4) * t295 - t101 * t300;
t160 = Ifges(5,4) * t273;
t365 = Ifges(5,1) * t264 + Ifges(5,5) * qJD(4) - Ifges(7,5) * t196 + t127 * Ifges(7,6) + t155 * Ifges(7,3) + t160;
t253 = mrSges(6,1) * t273;
t143 = -qJD(4) * mrSges(6,3) - t253;
t75 = -mrSges(7,1) * t127 - mrSges(7,2) * t196;
t364 = t75 - t143;
t110 = t135 * mrSges(6,1) + qJDD(4) * mrSges(6,2);
t278 = -qJDD(4) * mrSges(5,1) + mrSges(5,3) * t135 + t110;
t109 = mrSges(6,1) * t134 - qJDD(4) * mrSges(6,3);
t363 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t134 - t109;
t251 = mrSges(5,3) * t273;
t142 = -qJD(4) * mrSges(5,2) + t251;
t362 = t142 - t143;
t361 = (mrSges(3,1) * t177 - mrSges(3,2) * t174) * t168 - t101 * mrSges(4,1) - t102 * mrSges(4,2);
t318 = Ifges(6,6) * t176;
t319 = Ifges(6,6) * t173;
t360 = t173 * (-Ifges(6,2) * t176 + t319) + t176 * (Ifges(6,3) * t173 - t318);
t359 = t173 * t372 + t176 * t373;
t257 = t173 * t151 + t153 * t269 + t176 * t41;
t11 = -t271 * t79 + t257;
t354 = t11 * t176 - t12 * t173;
t311 = t173 * t79;
t7 = -qJDD(4) * qJ(5) + qJD(4) * (-qJD(5) + t311) - t257;
t353 = t173 * t8 - t176 * t7;
t211 = -t176 * Ifges(6,3) - t319;
t315 = t196 * Ifges(7,4);
t67 = t127 * Ifges(7,2) + t155 * Ifges(7,6) - t315;
t121 = Ifges(7,4) * t127;
t68 = -Ifges(7,1) * t196 + t155 * Ifges(7,5) + t121;
t350 = Ifges(6,5) * qJD(4) + qJD(2) * t211 + t172 * t68 + t175 * t67;
t222 = mrSges(6,2) * t176 - mrSges(6,3) * t173;
t129 = t222 * qJD(2);
t225 = mrSges(5,1) * t176 - mrSges(5,2) * t173;
t349 = t129 + (-mrSges(4,1) - t225) * qJD(2);
t285 = t171 * t174;
t276 = -t166 * t284 - t169 * t285;
t65 = -t170 * t125 + t167 * t276;
t60 = t167 * t125 + t170 * t276;
t224 = t175 * mrSges(7,1) - t172 * mrSges(7,2);
t223 = mrSges(7,1) * t172 + mrSges(7,2) * t175;
t346 = -qJ(5) * t374 + mrSges(5,2) - mrSges(6,3) - t223;
t288 = t168 * t176;
t31 = t170 * t288 - t173 * t60;
t33 = -t167 * t288 + t173 * t65;
t76 = t102 * t173 - t171 * t176;
t195 = -g(1) * t33 - g(2) * t31 - g(3) * t76;
t345 = mrSges(7,1) * t283 + mrSges(7,2) * t281 - t222 + t225;
t344 = -mrSges(4,1) - t345;
t343 = -m(7) * (pkin(5) + pkin(8)) - mrSges(6,1) - mrSges(5,3) - t224;
t342 = -mrSges(4,2) - t343;
t341 = qJD(2) ^ 2;
t340 = -t71 * Ifges(7,4) / 0.2e1 - t72 * Ifges(7,2) / 0.2e1 - t124 * Ifges(7,6) / 0.2e1;
t339 = t71 / 0.2e1;
t338 = t72 / 0.2e1;
t335 = t124 / 0.2e1;
t333 = -t196 / 0.2e1;
t323 = Ifges(5,4) * t173;
t322 = Ifges(5,4) * t176;
t321 = Ifges(7,4) * t172;
t320 = Ifges(7,4) * t175;
t301 = mrSges(4,2) * qJD(2);
t162 = pkin(5) * t273;
t39 = t162 - t45;
t298 = qJD(4) * t39;
t290 = t168 * t173;
t282 = t172 * t176;
t280 = t175 * t176;
t154 = pkin(2) * t287;
t275 = -t101 * pkin(3) + t154;
t266 = qJD(6) * t176;
t260 = Ifges(7,5) * t71 + Ifges(7,6) * t72 + Ifges(7,3) * t124;
t19 = -mrSges(7,1) * t72 + mrSges(7,2) * t71;
t258 = t19 + t363;
t256 = t142 + t364;
t255 = m(4) + m(5) + t374;
t244 = t172 * t266;
t239 = -t267 / 0.2e1;
t235 = -t263 / 0.2e1;
t53 = -t139 + t311;
t230 = t377 * pkin(2);
t229 = pkin(8) * t102 + t275;
t227 = t10 * t175 - t9 * t172;
t221 = Ifges(7,1) * t175 - t321;
t220 = Ifges(7,1) * t172 + t320;
t219 = t176 * Ifges(5,2) + t323;
t217 = -Ifges(7,2) * t172 + t320;
t216 = Ifges(7,2) * t175 + t321;
t214 = Ifges(7,5) * t175 - Ifges(7,6) * t172;
t213 = Ifges(7,5) * t172 + Ifges(7,6) * t175;
t212 = -t173 * Ifges(6,2) - t318;
t26 = t101 * t175 + t172 * t76;
t25 = -t101 * t172 + t175 * t76;
t77 = t102 * t176 + t171 * t173;
t206 = t61 * pkin(3) + t230;
t205 = -pkin(4) * t176 + t237;
t203 = -t167 * t284 - t170 * t174;
t57 = qJD(2) * t205 - t82;
t202 = t57 * (-mrSges(6,2) * t173 - mrSges(6,3) * t176);
t78 = -qJD(2) * pkin(3) - t82;
t201 = t78 * (mrSges(5,1) * t173 + mrSges(5,2) * t176);
t200 = t173 * (Ifges(5,1) * t176 - t323);
t190 = -pkin(8) * t60 + t206;
t189 = t203 * pkin(2);
t188 = t173 * t270 + t244;
t187 = t172 * t271 - t175 * t266;
t186 = t64 * pkin(3) + t189;
t185 = Ifges(7,5) * t176 + t173 * t220;
t184 = Ifges(7,6) * t176 + t173 * t216;
t183 = Ifges(7,3) * t176 + t173 * t213;
t182 = pkin(8) * t65 + t186;
t181 = qJD(6) * t227 + t228;
t161 = pkin(4) * t264;
t131 = -qJ(5) * t273 + t161;
t120 = t205 - t331;
t117 = Ifges(6,4) * qJD(4) + qJD(2) * t212;
t114 = Ifges(5,6) * qJD(4) + qJD(2) * t219;
t113 = -qJ(5) * t269 + t232;
t111 = t324 * t271;
t104 = qJD(2) * t210 + t161;
t95 = t125 * t168 * qJD(2);
t94 = qJD(2) * t102;
t81 = -mrSges(6,2) * t134 - mrSges(6,3) * t135;
t80 = mrSges(5,1) * t134 + mrSges(5,2) * t135;
t44 = -qJD(4) * pkin(4) + qJD(5) + t53;
t43 = t162 + t54;
t23 = -t102 * t271 + (qJD(4) * t171 - t95) * t176;
t22 = qJD(4) * t77 - t173 * t95;
t21 = t104 * t175 + t172 * t43;
t20 = -t104 * t172 + t175 * t43;
t18 = pkin(4) * t134 + t180;
t15 = t71 * Ifges(7,1) + t72 * Ifges(7,4) + t124 * Ifges(7,5);
t6 = -pkin(5) * t134 - t7;
t4 = qJD(6) * t25 + t172 * t22 + t175 * t94;
t3 = -qJD(6) * t26 - t172 * t94 + t175 * t22;
t14 = [t95 * t301 + m(2) * qJDD(1) + t25 * t35 + t26 * t36 + t3 * t87 + t4 * t86 + t278 * t76 + t277 * t22 + (-mrSges(3,1) * t174 - mrSges(3,2) * t177) * t341 * t168 + (t80 + t81) * t101 + t349 * t94 + t258 * t77 + t256 * t23 + t361 * qJDD(2) + (-m(2) - m(3) - t255) * g(3) + m(7) * (t1 * t26 + t10 * t4 + t2 * t25 + t23 * t39 + t3 * t9 + t6 * t77) + m(6) * (t101 * t18 + t22 * t44 - t23 * t45 + t57 * t94 - t7 * t77 + t76 * t8) + m(5) * (t101 * t40 + t11 * t77 - t12 * t76 + t22 * t53 + t23 * t54 + t78 * t94) + m(4) * (-t101 * t47 + t102 * t48 + t151 * t171 - t82 * t94 - t83 * t95) + m(3) * (qJDD(1) * t171 ^ 2 + (t105 * t177 + t106 * t174) * t168); (-Ifges(5,4) * t134 + Ifges(5,5) * qJDD(4) + t260) * t173 / 0.2e1 + (-m(4) * t154 - m(6) * (t229 + t366) - m(5) * t229 - m(7) * (-pkin(9) * t295 + t275 + t366) + t343 * t102 + t345 * t101 - t361) * g(3) + (-g(1) * t305 - g(2) * t306 + g(3) * t295 - t1 * t280 + t10 * t188 - t187 * t9 + t2 * t282) * mrSges(7,3) + (m(5) * t40 + t80) * (-pkin(3) - t331) + (t231 + t105) * mrSges(3,1) + t58 * t35 + t59 * t36 - t15 * t282 / 0.2e1 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + m(6) * (t113 * t57 + t120 * t18) + (-t203 * mrSges(3,1) - (t167 * t285 - t170 * t177) * mrSges(3,2) - m(5) * t182 - m(6) * (t182 + t368) - m(4) * t189 - m(7) * (pkin(9) * t305 + t186 + t368) + t344 * t64 - t342 * t65) * g(1) + t353 * mrSges(6,1) + t354 * mrSges(5,3) + t67 * t244 / 0.2e1 + t127 * (qJD(4) * t184 - t217 * t266) / 0.2e1 + t155 * (qJD(4) * t183 - t214 * t266) / 0.2e1 + (t169 * t296 + t47) * mrSges(4,1) + t135 * t322 / 0.2e1 + (m(4) * t82 - m(5) * t78 - m(6) * t57 - t349) * t93 + (-t166 * t296 - t48) * mrSges(4,2) + t359 * qJD(4) ^ 2 / 0.2e1 + t360 * t235 + ((-m(7) * t39 - t362 - t75) * t176 - m(5) * (t173 * t53 + t176 * t54) - m(6) * (t173 * t44 - t176 * t45) - m(4) * t83 + t301) * t96 + (-t114 / 0.2e1 + t350 / 0.2e1 + t45 * mrSges(6,1) - t54 * mrSges(5,3)) * t271 + (-t117 / 0.2e1 - t10 * mrSges(7,2) + t9 * mrSges(7,1) + t44 * mrSges(6,1) + t53 * mrSges(5,3) + t365 / 0.2e1) * t269 + (-t362 * t271 + t363 * t176 + t278 * t173 + m(6) * ((t45 * t173 + t44 * t176) * qJD(4) + t353) + m(5) * ((-t54 * t173 + t53 * t176) * qJD(4) + t354)) * t158 + (t176 * (-Ifges(5,2) * t173 + t322) + t200) * t263 / 0.2e1 - t40 * t225 - t134 * t219 / 0.2e1 + t18 * t222 - t135 * t212 / 0.2e1 + t134 * t211 / 0.2e1 + m(4) * (t166 * t48 + t169 * t47) * pkin(2) + (t168 * t243 - t106) * mrSges(3,2) + t39 * (-mrSges(7,1) * t188 + mrSges(7,2) * t187) - t173 * t379 + (-t377 * mrSges(3,1) - (-t167 * t177 - t170 * t285) * mrSges(3,2) - m(4) * t230 - m(5) * t190 - m(6) * (t190 + t369) - m(7) * (pkin(9) * t306 + t206 + t369) + t344 * t61 + t342 * t60) * g(2) + t277 * (t158 * t269 - t173 * t96) + t173 * t378 + t6 * t224 * t176 + t135 * t173 * Ifges(5,1) - t111 * t75 + t120 * t81 + t123 * t19 + t113 * t129 + t370 * t86 + t176 * t68 * t239 + t371 * t87 + (t1 * t59 + t10 * t370 - t111 * t39 + t123 * t6 + t2 * t58 + t371 * t9) * m(7) + (t173 * t373 - t176 * t372) * qJDD(4) / 0.2e1 + (qJD(4) * t185 - t221 * t266) * t333 + (Ifges(7,3) * t173 - t176 * t213) * t335 + (Ifges(7,6) * t173 - t176 * t216) * t338 + (Ifges(7,5) * t173 - t176 * t220) * t339 + t280 * t340 + qJD(4) * t201 + qJD(4) * t202 - t173 * (Ifges(6,4) * qJDD(4) - Ifges(6,2) * t135 + Ifges(6,6) * t134) / 0.2e1 + t176 * (Ifges(5,4) * t135 - Ifges(5,2) * t134 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t176 * (Ifges(6,5) * qJDD(4) - Ifges(6,6) * t135 + Ifges(6,3) * t134) / 0.2e1; m(4) * t151 + ((t172 * t86 + t175 * t87 + t277) * qJD(4) + m(6) * (qJD(4) * t44 - t7) + m(5) * (qJD(4) * t53 + t11) + m(7) * (t10 * t272 + t270 * t9 + t6) + t258) * t173 + (t256 * qJD(4) + m(6) * t367 + m(5) * (qJD(4) * t54 + t12) + m(7) * (t298 + t375) - t278 - t348) * t176 + (-t171 * g(3) + (-g(1) * t167 + g(2) * t170) * t168) * t255; (-pkin(4) * t8 - qJ(5) * t7 - qJD(5) * t45 - t131 * t57) * m(6) - (t127 * t216 + t155 * t213 - t196 * t220) * qJD(6) / 0.2e1 - (t127 * t184 + t155 * t183 - t185 * t196) * qJD(2) / 0.2e1 - t44 * t253 - t45 * t254 + (t19 - t109) * qJ(5) + (-m(7) * t181 - t348) * t337 + t67 * t239 - t21 * t86 - t20 * t87 - t42 * t75 + t8 * mrSges(6,2) - t11 * mrSges(5,2) + t12 * mrSges(5,1) - t7 * mrSges(6,3) + (t360 / 0.2e1 - t200 / 0.2e1) * t341 + (t346 * (-t170 * t290 - t176 * t60) + t376 * t31) * g(2) + t117 * t273 / 0.2e1 + (qJ(5) * t6 - t10 * t21 - t20 * t9 + t352 * t39) * m(7) - t68 * t268 / 0.2e1 + t155 * t39 * t224 + (-m(6) * t44 + t252 - t277) * t54 + t114 * t264 / 0.2e1 - t350 * t264 / 0.2e1 + t359 * t235 + (-m(6) * t45 - t251 + t362) * t53 + t6 * t223 + (Ifges(5,3) + Ifges(6,1)) * qJDD(4) + (t346 * t77 + t376 * t76) * g(3) + (t346 * (t167 * t290 + t176 * t65) + t376 * t33) * g(1) + (-t195 + t375) * mrSges(7,3) + (-t201 - t202 - t10 * (-mrSges(7,2) * t176 + mrSges(7,3) * t281) - t9 * (mrSges(7,1) * t176 - mrSges(7,3) * t283)) * qJD(2) + t364 * qJD(5) - (-Ifges(5,2) * t264 + t160 + t365) * t273 / 0.2e1 - pkin(4) * t110 - t131 * t129 + t372 * t134 + t373 * t135 + t214 * t335 + t217 * t338 + t221 * t339 + t172 * t340 + t175 * t15 / 0.2e1; -t364 * qJD(4) + (t129 + t355) * t264 + t110 + (t227 * t264 + t181 + t195 - t298) * m(7) + (t264 * t57 + t195 - t367) * m(6) + t348; -t379 + t378 - t39 * (-mrSges(7,1) * t196 + mrSges(7,2) * t127) + t196 * (Ifges(7,1) * t127 + t315) / 0.2e1 + t67 * t333 - t155 * (Ifges(7,5) * t127 + Ifges(7,6) * t196) / 0.2e1 - t9 * t86 + t10 * t87 - g(1) * ((t172 * t64 + t175 * t33) * mrSges(7,1) + (-t172 * t33 + t175 * t64) * mrSges(7,2)) - g(2) * ((t172 * t61 + t175 * t31) * mrSges(7,1) + (-t172 * t31 + t175 * t61) * mrSges(7,2)) - g(3) * (mrSges(7,1) * t25 - mrSges(7,2) * t26) + (-t10 * t196 + t127 * t9) * mrSges(7,3) + t260 - (Ifges(7,2) * t196 + t121 + t68) * t127 / 0.2e1;];
tau  = t14;
