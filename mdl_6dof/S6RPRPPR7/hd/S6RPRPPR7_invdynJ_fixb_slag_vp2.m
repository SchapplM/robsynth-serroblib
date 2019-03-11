% Calculate vector of inverse dynamics joint torques for
% S6RPRPPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:19
% EndTime: 2019-03-09 02:56:33
% DurationCPUTime: 11.88s
% Computational Cost: add. (5087->580), mult. (9919->718), div. (0->0), fcn. (6284->10), ass. (0->272)
t353 = -mrSges(6,2) + mrSges(5,1);
t176 = cos(qJ(3));
t352 = t176 / 0.2e1;
t347 = Ifges(5,6) - Ifges(6,5);
t165 = qJ(3) + pkin(9);
t159 = cos(t165);
t303 = pkin(3) * t176;
t351 = pkin(4) * t159 + t303;
t158 = sin(t165);
t173 = sin(qJ(3));
t212 = mrSges(4,1) * t173 + mrSges(4,2) * t176;
t350 = -mrSges(5,2) * t159 - t353 * t158 - t212;
t172 = sin(qJ(6));
t175 = cos(qJ(6));
t237 = qJD(1) * qJD(3);
t124 = qJDD(1) * t176 - t173 * t237;
t125 = -qJDD(1) * t173 - t176 * t237;
t170 = cos(pkin(9));
t262 = sin(pkin(9));
t77 = t124 * t262 - t170 * t125;
t118 = t170 * t173 + t176 * t262;
t110 = t118 * qJD(1);
t84 = -qJD(3) * t172 + t110 * t175;
t32 = qJD(6) * t84 + qJDD(3) * t175 + t172 * t77;
t321 = t32 / 0.2e1;
t85 = qJD(3) * t175 + t110 * t172;
t33 = -qJD(6) * t85 - qJDD(3) * t172 + t175 * t77;
t320 = t33 / 0.2e1;
t78 = t170 * t124 + t125 * t262;
t76 = qJDD(6) + t78;
t319 = t76 / 0.2e1;
t315 = -m(3) - m(4);
t314 = -m(6) - m(7);
t349 = t84 * Ifges(7,6);
t348 = Ifges(5,5) - Ifges(6,4);
t41 = -mrSges(7,1) * t84 + mrSges(7,2) * t85;
t290 = mrSges(6,1) * t110;
t91 = -qJD(3) * mrSges(6,3) + t290;
t346 = t41 - t91;
t102 = Ifges(5,4) * t110;
t226 = t262 * t173;
t245 = qJD(1) * t176;
t112 = -qJD(1) * t226 + t170 * t245;
t103 = qJD(6) + t112;
t341 = t103 * Ifges(7,3);
t345 = t112 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t85 * Ifges(7,5) - t102 + t341 + t349;
t66 = -qJDD(3) * mrSges(5,2) - mrSges(5,3) * t77;
t68 = mrSges(6,1) * t77 - qJDD(3) * mrSges(6,3);
t344 = t66 - t68;
t288 = mrSges(5,3) * t110;
t89 = -qJD(3) * mrSges(5,2) - t288;
t343 = -t89 + t91;
t287 = mrSges(5,3) * t112;
t289 = mrSges(6,1) * t112;
t342 = -qJD(3) * t353 + t287 + t289;
t210 = mrSges(7,1) * t175 - mrSges(7,2) * t172;
t297 = t110 * pkin(5);
t246 = qJD(1) * t173;
t178 = -pkin(1) - pkin(7);
t138 = qJD(1) * t178 + qJD(2);
t257 = t138 * t173;
t99 = -qJ(4) * t246 + t257;
t268 = t170 * t99;
t126 = t176 * t138;
t100 = -qJ(4) * t245 + t126;
t96 = qJD(3) * pkin(3) + t100;
t50 = t262 * t96 + t268;
t45 = -qJD(3) * qJ(5) - t50;
t31 = -t45 - t297;
t340 = t31 * t210;
t263 = qJDD(3) / 0.2e1;
t177 = cos(qJ(1));
t299 = g(2) * t177;
t174 = sin(qJ(1));
t300 = g(1) * t174;
t337 = t299 - t300;
t339 = t159 * t337;
t167 = qJD(1) * qJD(2);
t139 = qJDD(1) * qJ(2) + t167;
t213 = mrSges(4,1) * t176 - mrSges(4,2) * t173;
t338 = -m(5) * t303 - mrSges(5,1) * t159 + mrSges(5,2) * t158 - t213;
t137 = qJDD(1) * t178 + qJDD(2);
t244 = qJD(3) * t173;
t82 = t176 * t137 - t138 * t244;
t243 = qJD(3) * t176;
t83 = t173 * t137 + t138 * t243;
t198 = t173 * t83 + t176 * t82;
t18 = mrSges(7,1) * t76 - mrSges(7,3) * t32;
t19 = -mrSges(7,2) * t76 + mrSges(7,3) * t33;
t201 = -t172 * t19 - t175 * t18;
t283 = Ifges(4,4) * t176;
t336 = (-Ifges(4,1) * t173 - t283) * t352 + qJ(2) * t213;
t335 = 0.2e1 * t263;
t130 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t246;
t131 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t245;
t334 = (t130 * t176 - t131 * t173) * qJD(3);
t87 = -pkin(3) * t125 + qJDD(4) + t139;
t181 = -qJ(5) * t78 - qJD(5) * t112 + t87;
t313 = pkin(4) + pkin(8);
t11 = t313 * t77 + t181;
t93 = t262 * t99;
t49 = t170 * t96 - t93;
t220 = qJD(5) - t49;
t296 = t112 * pkin(5);
t28 = -qJD(3) * t313 + t220 + t296;
t127 = pkin(3) * t246 + qJD(1) * qJ(2) + qJD(4);
t194 = -qJ(5) * t112 + t127;
t34 = t110 * t313 + t194;
t7 = -t172 * t34 + t175 * t28;
t236 = qJD(1) * qJD(4);
t48 = qJDD(3) * pkin(3) - qJ(4) * t124 - t176 * t236 + t82;
t54 = qJ(4) * t125 - t173 * t236 + t83;
t22 = t170 * t48 - t262 * t54;
t214 = qJDD(5) - t22;
t9 = pkin(5) * t78 - qJDD(3) * t313 + t214;
t1 = qJD(6) * t7 + t11 * t175 + t172 * t9;
t8 = t172 * t28 + t175 * t34;
t2 = -qJD(6) * t8 - t11 * t172 + t175 * t9;
t216 = t172 * t7 - t175 * t8;
t333 = qJD(6) * t216 - t1 * t172 - t175 * t2;
t111 = qJD(3) * t226 - t170 * t243;
t113 = t118 * qJD(3);
t119 = t170 * t176 - t226;
t23 = t170 * t54 + t262 * t48;
t15 = -qJDD(3) * qJ(5) - qJD(3) * qJD(5) - t23;
t17 = -qJDD(3) * pkin(4) + t214;
t43 = -qJD(3) * pkin(4) + t220;
t332 = -t111 * t45 - t113 * t43 + t118 * t15 + t119 * t17;
t331 = t111 * t50 + t113 * t49 - t118 * t23 - t119 * t22;
t330 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t329 = m(4) * t198 + t176 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t124) + t173 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t125);
t52 = -mrSges(7,2) * t103 + mrSges(7,3) * t84;
t53 = mrSges(7,1) * t103 - mrSges(7,3) * t85;
t200 = t172 * t53 - t175 * t52;
t69 = t78 * mrSges(6,1) + qJDD(3) * mrSges(6,2);
t328 = -t200 * qJD(6) - t201 + t69;
t147 = t159 * qJ(5);
t327 = -mrSges(3,3) - m(7) * (pkin(8) * t158 - t147) - t158 * mrSges(7,3) - (-m(6) * qJ(5) - mrSges(6,3)) * t159 + mrSges(2,2) + t350;
t326 = -m(7) * pkin(5) - mrSges(2,1) - mrSges(6,1) + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t295 = t85 * Ifges(7,4);
t26 = t84 * Ifges(7,2) + t103 * Ifges(7,6) + t295;
t232 = -t175 * t26 / 0.2e1;
t81 = Ifges(7,4) * t84;
t27 = t85 * Ifges(7,1) + t103 * Ifges(7,5) + t81;
t305 = -t172 / 0.2e1;
t51 = pkin(4) * t110 + t194;
t325 = -t127 * mrSges(5,1) + t51 * mrSges(6,2) + t27 * t305 + t232;
t324 = t7 * mrSges(7,1) + t127 * mrSges(5,2) - t8 * mrSges(7,2) - t51 * mrSges(6,3);
t323 = qJD(1) ^ 2;
t322 = Ifges(7,1) * t321 + Ifges(7,4) * t320 + Ifges(7,5) * t319;
t318 = -t84 / 0.2e1;
t317 = -t85 / 0.2e1;
t316 = t85 / 0.2e1;
t312 = -t103 / 0.2e1;
t311 = -t110 / 0.2e1;
t310 = t110 / 0.2e1;
t309 = -t112 / 0.2e1;
t308 = t112 / 0.2e1;
t304 = pkin(3) * t170;
t302 = pkin(4) * t158;
t298 = g(3) * t158;
t163 = t173 * pkin(3);
t294 = -qJD(3) / 0.2e1;
t293 = qJD(3) / 0.2e1;
t286 = mrSges(7,3) * t172;
t285 = mrSges(7,3) * t175;
t284 = Ifges(4,4) * t173;
t282 = Ifges(7,4) * t172;
t281 = Ifges(7,4) * t175;
t10 = -pkin(5) * t77 - t15;
t280 = t10 * t118;
t277 = t112 * Ifges(5,4);
t276 = t112 * Ifges(6,6);
t261 = t118 * t172;
t260 = t118 * t175;
t256 = t158 * t174;
t255 = t172 * t174;
t254 = t172 * t177;
t252 = t174 * t175;
t251 = t175 * t177;
t151 = qJ(2) + t163;
t249 = qJ(4) - t178;
t248 = t147 - t163;
t247 = t177 * pkin(1) + t174 * qJ(2);
t241 = qJD(6) * t172;
t240 = qJD(6) * t175;
t239 = qJDD(1) * mrSges(3,2);
t140 = pkin(3) * t243 + qJD(2);
t238 = -m(5) + t314;
t235 = Ifges(7,5) * t32 + Ifges(7,6) * t33 + Ifges(7,3) * t76;
t234 = -t89 - t346;
t155 = pkin(3) * t245;
t231 = qJ(5) * t256 + t174 * t351;
t150 = -pkin(4) - t304;
t230 = t262 * pkin(3);
t228 = -t241 / 0.2e1;
t162 = t177 * qJ(2);
t227 = -pkin(1) * t174 + t162;
t225 = -t237 / 0.2e1;
t129 = t249 * t176;
t223 = (t139 + t167) * qJ(2);
t222 = qJ(5) * t110 + t155;
t221 = -m(7) * t313 - mrSges(7,3);
t217 = t172 * t8 + t175 * t7;
t215 = -qJ(5) * t119 + t151;
t97 = -qJD(4) * t176 + t244 * t249;
t98 = -qJD(3) * t129 - qJD(4) * t173;
t55 = -t170 * t97 + t262 * t98;
t59 = t100 * t262 + t268;
t209 = mrSges(7,1) * t172 + mrSges(7,2) * t175;
t208 = -t159 * mrSges(6,2) + t158 * mrSges(6,3);
t207 = t176 * Ifges(4,1) - t284;
t206 = Ifges(7,1) * t172 + t281;
t205 = -t173 * Ifges(4,2) + t283;
t204 = Ifges(7,2) * t175 + t282;
t203 = -Ifges(4,5) * t173 - Ifges(4,6) * t176;
t202 = Ifges(7,5) * t172 + Ifges(7,6) * t175;
t44 = t118 * t313 + t215;
t128 = t249 * t173;
t79 = -t128 * t262 + t170 * t129;
t57 = t119 * pkin(5) + t79;
t21 = t172 * t57 + t175 * t44;
t20 = -t172 * t44 + t175 * t57;
t199 = -t172 * t52 - t175 * t53;
t171 = -qJ(4) - pkin(7);
t197 = t177 * t163 + t174 * t171 + t227;
t196 = t174 * t163 - t171 * t177 + t247;
t71 = -mrSges(6,2) * t110 - mrSges(6,3) * t112;
t195 = t200 - t71;
t56 = t170 * t98 + t262 * t97;
t193 = -t111 * t172 + t118 * t240;
t192 = t111 * t175 + t118 * t241;
t189 = t173 * (-Ifges(4,2) * t176 - t284);
t60 = t170 * t100 - t93;
t80 = -t170 * t128 - t129 * t262;
t187 = -t199 + t342;
t185 = qJ(5) * t113 - qJD(5) * t119 + t140;
t157 = -qJDD(1) * pkin(1) + qJDD(2);
t148 = t230 + qJ(5);
t121 = t212 * qJD(1);
t115 = Ifges(4,5) * qJD(3) + qJD(1) * t207;
t114 = Ifges(4,6) * qJD(3) + qJD(1) * t205;
t107 = -t159 * t255 + t251;
t106 = -t159 * t252 - t254;
t105 = -t159 * t254 - t252;
t104 = -t159 * t251 + t255;
t101 = Ifges(6,6) * t110;
t75 = t78 * mrSges(6,3);
t74 = t78 * mrSges(5,2);
t72 = pkin(4) * t118 + t215;
t70 = mrSges(5,1) * t110 + mrSges(5,2) * t112;
t67 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t78;
t64 = -t110 * Ifges(5,2) + Ifges(5,6) * qJD(3) + t277;
t63 = Ifges(6,4) * qJD(3) - t112 * Ifges(6,2) + t101;
t62 = Ifges(6,5) * qJD(3) + t110 * Ifges(6,3) - t276;
t61 = pkin(4) * t112 + t222;
t58 = -t118 * pkin(5) + t80;
t40 = -pkin(4) * t111 + t185;
t39 = t112 * t313 + t222;
t38 = t60 - t296;
t37 = t59 - t297;
t36 = t111 * pkin(5) + t56;
t35 = -t113 * pkin(5) + t55;
t24 = -t111 * t313 + t185;
t16 = pkin(4) * t77 + t181;
t14 = t172 * t37 + t175 * t39;
t13 = -t172 * t39 + t175 * t37;
t12 = -mrSges(7,1) * t33 + mrSges(7,2) * t32;
t5 = t32 * Ifges(7,4) + t33 * Ifges(7,2) + t76 * Ifges(7,6);
t4 = -qJD(6) * t21 - t172 * t24 + t175 * t35;
t3 = qJD(6) * t20 + t172 * t35 + t175 * t24;
t6 = [(-m(3) * t227 - m(4) * t162 - m(5) * t197 - t105 * mrSges(7,1) - t104 * mrSges(7,2) + t314 * (t177 * t302 + t197) + (-m(4) * t178 - t326) * t174 + t327 * t177) * g(1) + (-m(5) * t196 - t107 * mrSges(7,1) - t106 * mrSges(7,2) + t315 * t247 + t314 * (pkin(4) * t256 + t196) + (-m(4) * pkin(7) + t326) * t177 + t327 * t174) * g(2) + m(6) * (t16 * t72 + t40 * t51) + m(5) * (t127 * t140 + t151 * t87) - t198 * mrSges(4,3) + t336 * t237 + (t212 + 0.2e1 * mrSges(3,3)) * t139 + (-m(5) * t49 + m(6) * t43 + t342) * t55 + (m(5) * t50 - m(6) * t45 - t343) * t56 + (m(5) * t23 - m(6) * t15 + t344) * t80 + (Ifges(5,4) * t308 - Ifges(6,6) * t309 - Ifges(6,3) * t310 + Ifges(5,2) * t311 + t64 / 0.2e1 - t62 / 0.2e1 + t347 * t293 + t325) * t111 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (-t341 / 0.2e1 - t349 / 0.2e1 - Ifges(7,5) * t316 - Ifges(5,1) * t308 + Ifges(6,2) * t309 + Ifges(6,6) * t310 - Ifges(5,4) * t311 + t63 / 0.2e1 - t348 * t293 - t324 - t345 / 0.2e1) * t113 + m(4) * t223 + t103 * (Ifges(7,5) * t193 - Ifges(7,6) * t192) / 0.2e1 + t84 * (Ifges(7,4) * t193 - Ifges(7,2) * t192) / 0.2e1 + t3 * t52 + t4 * t53 + t36 * t41 + t20 * t18 + t21 * t19 - t114 * t243 / 0.2e1 - t115 * t244 / 0.2e1 - pkin(1) * t239 + (t87 * mrSges(5,1) - t16 * mrSges(6,2) + t202 * t319 + t204 * t320 + t206 * t321 + t26 * t228 + (-Ifges(6,6) - Ifges(5,4)) * t78 + (Ifges(6,3) + Ifges(5,2)) * t77 - t347 * t335) * t118 + (Ifges(4,5) * t176 - Ifges(4,6) * t173) * t263 + m(7) * (t1 * t21 + t10 * t58 + t2 * t20 + t3 * t8 + t31 * t36 + t4 * t7) - t210 * t280 + t329 * t178 + t331 * mrSges(5,3) + t332 * mrSges(6,1) + t31 * (mrSges(7,1) * t192 + mrSges(7,2) * t193) + (t87 * mrSges(5,2) - t16 * mrSges(6,3) + Ifges(5,5) * t263 + Ifges(7,5) * t321 + Ifges(7,6) * t320 + Ifges(7,3) * t319 + Ifges(6,2) * t78 + (-Ifges(5,4) / 0.2e1 - Ifges(6,6)) * t77 - t335 * Ifges(6,4) + t330) * t119 + (-Ifges(5,4) * t77 + Ifges(5,5) * qJDD(3) + t235) * t119 / 0.2e1 + (Ifges(7,1) * t193 - Ifges(7,4) * t192) * t316 + t72 * (-t77 * mrSges(6,2) - t75) + (t1 * t260 - t192 * t8 - t193 * t7 - t2 * t261) * mrSges(7,3) + t78 * Ifges(5,1) * t119 + (Ifges(4,1) * t124 + Ifges(4,4) * t125 + Ifges(4,5) * qJDD(3)) * t352 + t261 * t322 + t125 * t205 / 0.2e1 + t124 * t207 / 0.2e1 + t58 * t12 + t40 * t71 + (qJD(6) * t27 + t5) * t260 / 0.2e1 + t189 * t225 + t178 * t334 + qJD(2) * t121 + qJ(2) * (-mrSges(4,1) * t125 + mrSges(4,2) * t124) + t140 * t70 + t151 * (t77 * mrSges(5,1) + t74) + t157 * mrSges(3,2) - t173 * (Ifges(4,4) * t124 + Ifges(4,2) * t125 + Ifges(4,6) * qJDD(3)) / 0.2e1 + m(3) * (-pkin(1) * t157 + t223) + qJD(3) ^ 2 * t203 / 0.2e1 + (-m(5) * t22 + m(6) * t17 - t67 + t69) * t79; t239 + t334 + (qJ(2) * t315 - mrSges(3,3)) * t323 + (t12 + t344) * t118 + t234 * t111 + t187 * t113 + (-t328 + t67) * t119 - m(5) * t331 - m(6) * t332 + m(3) * t157 + m(7) * (-t111 * t31 + t217 * t113 + t119 * t333 + t280) + (-m(5) * t127 - m(6) * t51 + m(7) * t216 - t121 + t195 - t70) * qJD(1) + t337 * (-t238 - t315) + t329; (-g(1) * t231 - t148 * t15 + t150 * t17 - t43 * t59 - t51 * t61 + (-qJD(5) + t60) * t45) * m(6) + (t232 + t340) * qJD(6) + t338 * t300 + (t189 / 0.2e1 - t336) * t323 + (-t240 * t8 + t241 * t7) * mrSges(7,3) - t342 * t59 + t343 * t60 + (-t102 + t345) * t310 + t346 * qJD(5) + (-Ifges(5,2) * t310 + Ifges(6,3) * t311 + t202 * t312 + t204 * t318 + t206 * t317 - t285 * t8 + t286 * t7 - t294 * t347 + t325 + t340) * t112 - t347 * t77 + t348 * t78 + (-Ifges(5,1) * t309 - Ifges(7,5) * t317 + Ifges(6,2) * t308 - Ifges(7,6) * t318 - Ifges(7,3) * t312 - t294 * t348 + t324) * t110 - t14 * t52 - t13 * t53 - t38 * t41 - t15 * mrSges(6,3) + t17 * mrSges(6,2) + t22 * mrSges(5,1) - t23 * mrSges(5,2) - t2 * t285 - t1 * t286 - t49 * t288 + (t276 + t64) * t308 + t131 * t257 - t130 * t126 + (t101 + t63) * t311 + (-t277 + t62) * t309 + (Ifges(4,3) + Ifges(6,1) + Ifges(5,3)) * qJDD(3) + t114 * t245 / 0.2e1 + t115 * t246 / 0.2e1 + (t10 * t148 - t13 * t7 - t14 * t8 + (-t38 + qJD(5)) * t31) * m(7) + (-t127 * t155 + t49 * t59 - t50 * t60 + (t170 * t22 + t23 * t262) * pkin(3)) * m(5) + t66 * t230 + (Ifges(7,5) * t175 - Ifges(7,6) * t172) * t319 + (-Ifges(7,2) * t172 + t281) * t320 + (Ifges(7,1) * t175 - t282) * t321 + t175 * t322 + t67 * t304 + t5 * t305 - t70 * t155 + (-m(7) * t231 + (-t208 - (m(7) * pkin(8) + mrSges(7,3)) * t159 - t209 * t158) * t174) * g(1) + t50 * t287 + t43 * t290 + t10 * t209 + (m(7) * t303 - t221 * t159 - (-m(7) * qJ(5) - t209) * t158 - m(6) * (-qJ(5) * t158 - t351) + t208 - t338) * t299 + (t12 - t68) * t148 - t61 * t71 + t82 * mrSges(4,1) - t83 * mrSges(4,2) - (t103 * t202 + t204 * t84 + t206 * t85) * qJD(6) / 0.2e1 + t203 * t225 + t27 * t228 + Ifges(4,5) * t124 + Ifges(4,6) * t125 + t150 * t69 - t45 * t289 + (-m(7) * t333 + t240 * t52 - t241 * t53 - t201) * (-pkin(8) + t150) + (-m(7) * t248 - t221 * t158 + m(5) * t163 - m(6) * (t248 - t302) + (-t209 - mrSges(6,3)) * t159 - t350) * g(3); -t172 * t18 + t175 * t19 + t74 - t75 + t353 * t77 + t199 * qJD(6) - t234 * t110 - t187 * t112 + (g(1) * t177 + g(2) * t174) * t238 + (t1 * t175 - t103 * t217 + t110 * t31 - t2 * t172) * m(7) + (-t110 * t45 - t112 * t43 + t16) * m(6) + (t110 * t50 + t112 * t49 + t87) * m(5); -t346 * qJD(3) - t195 * t112 + t314 * t298 + (-qJD(3) * t31 - t112 * t216 - t333 - t339) * m(7) + (qJD(3) * t45 + t112 * t51 + t17 - t339) * m(6) + t328; -t31 * (mrSges(7,1) * t85 + mrSges(7,2) * t84) + (Ifges(7,1) * t84 - t295) * t317 + t26 * t316 + (Ifges(7,5) * t84 - Ifges(7,6) * t85) * t312 - t7 * t52 + t8 * t53 - g(1) * (mrSges(7,1) * t106 - mrSges(7,2) * t107) - g(2) * (-mrSges(7,1) * t104 + mrSges(7,2) * t105) - t210 * t298 + (t7 * t84 + t8 * t85) * mrSges(7,3) + t235 + (-Ifges(7,2) * t85 + t27 + t81) * t318 + t330;];
tau  = t6;
