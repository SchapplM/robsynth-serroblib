% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:16
% EndTime: 2019-03-09 01:46:31
% DurationCPUTime: 10.00s
% Computational Cost: add. (5575->541), mult. (10395->720), div. (0->0), fcn. (6444->12), ass. (0->248)
t334 = -m(7) - m(6);
t173 = sin(pkin(10));
t175 = cos(pkin(10));
t179 = sin(qJ(4));
t181 = cos(qJ(4));
t122 = t173 * t181 + t175 * t179;
t116 = t122 * qJD(4);
t174 = sin(pkin(9));
t266 = cos(pkin(9));
t224 = t181 * t266;
t225 = t179 * t266;
t341 = (-t173 * t225 + t175 * t224) * qJD(1) + t174 * t116;
t253 = qJD(3) * t179;
t182 = -pkin(1) - pkin(2);
t141 = qJD(1) * t182 + qJD(2);
t163 = t266 * qJ(2);
t106 = qJD(1) * t163 + t174 * t141;
t97 = -qJD(1) * pkin(7) + t106;
t68 = t253 + (-qJ(5) * qJD(1) + t97) * t181;
t272 = t173 * t68;
t255 = qJD(1) * t179;
t78 = t181 * qJD(3) - t179 * t97;
t67 = qJ(5) * t255 + t78;
t61 = qJD(4) * pkin(4) + t67;
t27 = t175 * t61 - t272;
t25 = -qJD(4) * pkin(5) - t27;
t115 = t122 * qJD(1);
t281 = mrSges(6,3) * t115;
t178 = sin(qJ(6));
t180 = cos(qJ(6));
t87 = qJD(4) * t180 + t115 * t178;
t88 = qJD(4) * t178 - t115 * t180;
t284 = -qJD(4) * mrSges(6,1) - mrSges(7,1) * t87 + mrSges(7,2) * t88 - t281;
t340 = m(6) * t27 - m(7) * t25 - t284;
t339 = m(7) * pkin(8) + mrSges(7,3);
t139 = mrSges(5,1) * t181 - mrSges(5,2) * t179;
t338 = -m(5) * pkin(3) - mrSges(4,1) - t139;
t140 = qJDD(1) * t182 + qJDD(2);
t246 = qJD(1) * qJD(2);
t142 = qJDD(1) * qJ(2) + t246;
t91 = t174 * t140 + t266 * t142;
t85 = -qJDD(1) * pkin(7) + t91;
t337 = qJD(3) * qJD(4) + t85;
t120 = t173 * t179 - t175 * t181;
t114 = t120 * qJD(1);
t108 = qJD(6) - t114;
t335 = -m(4) - m(5);
t333 = t87 * Ifges(7,6);
t332 = mrSges(3,1) + mrSges(2,1);
t331 = -mrSges(3,3) + mrSges(2,2);
t245 = qJD(1) * qJD(4);
t128 = -qJDD(1) * t181 + t179 * t245;
t129 = -qJDD(1) * t179 - t181 * t245;
t76 = t128 * t173 + t129 * t175;
t40 = qJD(6) * t87 + qJDD(4) * t178 + t180 * t76;
t41 = -qJD(6) * t88 + qJDD(4) * t180 - t178 * t76;
t13 = -mrSges(7,1) * t41 + mrSges(7,2) * t40;
t65 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t76;
t330 = t13 - t65;
t256 = qJD(1) * t174;
t104 = t120 * t174;
t81 = -t180 * t104 - t178 * t266;
t329 = -qJD(6) * t81 + t178 * t341 - t180 * t256;
t80 = t178 * t104 - t180 * t266;
t328 = qJD(6) * t80 - t178 * t256 - t180 * t341;
t107 = Ifges(6,4) * t114;
t325 = t108 * Ifges(7,3);
t324 = t114 * Ifges(6,2);
t267 = qJDD(4) / 0.2e1;
t323 = Ifges(6,5) * qJD(4);
t322 = Ifges(6,6) * qJD(4);
t277 = Ifges(5,4) * t181;
t278 = Ifges(5,4) * t179;
t321 = t179 * (-Ifges(5,1) * t181 + t278) + t181 * (Ifges(5,2) * t179 - t277);
t251 = qJD(4) * t181;
t252 = qJD(4) * t179;
t320 = t173 * t252 - t175 * t251;
t38 = t179 * qJDD(3) + t181 * t337 - t252 * t97;
t164 = t181 * qJDD(3);
t79 = t181 * t97 + t253;
t39 = -qJD(4) * t79 - t179 * t85 + t164;
t319 = -t179 * t39 + t181 * t38;
t75 = t128 * t175 - t129 * t173;
t74 = qJDD(6) - t75;
t17 = mrSges(7,1) * t74 - mrSges(7,3) * t40;
t18 = -mrSges(7,2) * t74 + mrSges(7,3) * t41;
t318 = -t178 * t17 + t180 * t18;
t317 = 0.2e1 * t267;
t316 = -t334 - t335;
t90 = t140 * t266 - t174 * t142;
t84 = qJDD(1) * pkin(3) - t90;
t314 = m(5) * t84 - mrSges(5,1) * t128 + mrSges(5,2) * t129;
t171 = qJ(4) + pkin(10);
t160 = cos(t171);
t213 = -mrSges(7,1) * t180 + mrSges(7,2) * t178;
t192 = m(7) * pkin(5) - t213;
t159 = sin(t171);
t214 = mrSges(6,1) * t160 - mrSges(6,2) * t159;
t312 = -t192 * t160 - t214;
t57 = t175 * t68;
t28 = t173 * t61 + t57;
t26 = qJD(4) * pkin(8) + t28;
t254 = qJD(1) * t181;
t105 = -qJ(2) * t256 + t141 * t266;
t96 = qJD(1) * pkin(3) - t105;
t82 = pkin(4) * t254 + qJD(5) + t96;
t43 = -t114 * pkin(5) + t115 * pkin(8) + t82;
t11 = -t178 * t26 + t180 * t43;
t54 = -t128 * pkin(4) + qJDD(5) + t84;
t16 = -t75 * pkin(5) - t76 * pkin(8) + t54;
t244 = qJD(1) * qJD(5);
t23 = -t97 * t251 + qJDD(4) * pkin(4) - qJ(5) * t129 + t164 + (t244 - t337) * t179;
t24 = qJ(5) * t128 - t181 * t244 + t38;
t8 = t173 * t23 + t175 * t24;
t6 = qJDD(4) * pkin(8) + t8;
t1 = qJD(6) * t11 + t16 * t178 + t180 * t6;
t12 = t178 * t43 + t180 * t26;
t2 = -qJD(6) * t12 + t16 * t180 - t178 * t6;
t311 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t295 = t178 / 0.2e1;
t292 = Ifges(7,4) * t88;
t34 = Ifges(7,2) * t87 + Ifges(7,6) * t108 + t292;
t86 = Ifges(7,4) * t87;
t35 = Ifges(7,1) * t88 + Ifges(7,5) * t108 + t86;
t310 = -t82 * mrSges(6,2) - t180 * t35 / 0.2e1 + t34 * t295;
t309 = -m(5) * pkin(7) + mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t109 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t128;
t110 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t129;
t236 = mrSges(5,3) * t255;
t137 = qJD(4) * mrSges(5,1) + t236;
t235 = mrSges(5,3) * t254;
t138 = -qJD(4) * mrSges(5,2) - t235;
t308 = m(5) * ((-t179 * t79 - t181 * t78) * qJD(4) + t319) - t138 * t252 - t137 * t251 - t110 * t179 + t109 * t181;
t307 = m(4) * (-t105 * t174 + t106 * t266) + m(5) * (t174 * t96 + t224 * t79 - t225 * t78) + t138 * t224 - t137 * t225;
t306 = t82 * mrSges(6,1) + t11 * mrSges(7,1) - t12 * mrSges(7,2);
t183 = qJD(1) ^ 2;
t305 = t40 / 0.2e1;
t304 = t41 / 0.2e1;
t303 = t74 / 0.2e1;
t302 = -t87 / 0.2e1;
t301 = -t88 / 0.2e1;
t300 = t88 / 0.2e1;
t299 = -t108 / 0.2e1;
t298 = -t114 / 0.2e1;
t297 = -t115 / 0.2e1;
t296 = t115 / 0.2e1;
t294 = cos(qJ(1));
t293 = sin(qJ(1));
t291 = pkin(4) * t173;
t290 = pkin(4) * t175;
t288 = pkin(8) * t159;
t168 = t181 * pkin(4);
t283 = mrSges(4,1) * t174;
t282 = mrSges(6,3) * t114;
t280 = mrSges(7,3) * t178;
t279 = mrSges(7,3) * t180;
t276 = Ifges(6,4) * t115;
t275 = Ifges(7,4) * t178;
t274 = Ifges(7,4) * t180;
t273 = t159 * mrSges(7,3);
t263 = t122 * t178;
t262 = t122 * t180;
t261 = t160 * t178;
t260 = t160 * t180;
t131 = t174 * t182 + t163;
t125 = -pkin(7) + t131;
t259 = qJ(5) - t125;
t258 = t294 * pkin(1) + t293 * qJ(2);
t257 = mrSges(4,1) * qJDD(1);
t250 = qJD(6) * t178;
t249 = qJD(6) * t180;
t248 = qJDD(1) * mrSges(3,1);
t247 = qJDD(1) * mrSges(4,2);
t158 = t174 * qJD(2);
t240 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t74;
t237 = pkin(4) * t255;
t232 = t294 * pkin(2) + t258;
t42 = -t75 * mrSges(6,1) + t76 * mrSges(6,2);
t226 = t249 / 0.2e1;
t223 = -t245 / 0.2e1;
t222 = t259 * t179;
t221 = t266 * qJD(2);
t220 = qJD(4) * t259;
t217 = -pkin(1) * t293 + t294 * qJ(2);
t133 = -pkin(4) * t252 + t158;
t121 = -t174 * t293 - t266 * t294;
t123 = t174 * t294 - t266 * t293;
t216 = -g(1) * t123 + g(2) * t121;
t215 = mrSges(5,1) * t179 + mrSges(5,2) * t181;
t212 = mrSges(7,1) * t178 + mrSges(7,2) * t180;
t211 = -t179 * Ifges(5,1) - t277;
t210 = Ifges(7,1) * t180 - t275;
t209 = -t181 * Ifges(5,2) - t278;
t208 = -Ifges(7,2) * t178 + t274;
t207 = -Ifges(5,5) * t181 + Ifges(5,6) * t179;
t206 = Ifges(7,5) * t180 - Ifges(7,6) * t178;
t205 = -t11 * t178 + t12 * t180;
t7 = -t173 * t24 + t175 * t23;
t100 = t259 * t181;
t48 = -t175 * t100 + t173 * t222;
t130 = -t174 * qJ(2) + t182 * t266;
t124 = pkin(3) - t130;
t111 = t168 + t124;
t50 = -t120 * pkin(5) + t122 * pkin(8) + t111;
t20 = t178 * t50 + t180 * t48;
t19 = -t178 * t48 + t180 * t50;
t51 = -mrSges(7,2) * t108 + mrSges(7,3) * t87;
t52 = mrSges(7,1) * t108 - mrSges(7,3) * t88;
t204 = -t178 * t52 + t180 * t51;
t202 = t221 - qJD(5);
t98 = -qJD(4) * mrSges(6,2) + t282;
t201 = -t204 - t98;
t200 = t25 * t212;
t199 = t96 * t215;
t198 = t122 * t249 - t178 * t320;
t197 = t122 * t250 + t180 * t320;
t191 = -pkin(2) * t293 + t217;
t188 = t1 * t180 - t178 * t2 + (-t11 * t180 - t12 * t178) * qJD(6);
t185 = -t179 * t202 + t181 * t220;
t177 = -qJ(5) - pkin(7);
t157 = -qJDD(1) * pkin(1) + qJDD(2);
t155 = t168 + pkin(3);
t151 = -pkin(5) - t290;
t126 = t139 * qJD(1);
t119 = Ifges(5,5) * qJD(4) + qJD(1) * t211;
t118 = Ifges(5,6) * qJD(4) + qJD(1) * t209;
t103 = t122 * t174;
t69 = -mrSges(6,1) * t114 - mrSges(6,2) * t115;
t66 = t179 * t220 + t181 * t202;
t64 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t75;
t63 = -t121 * t260 + t123 * t178;
t62 = t121 * t261 + t123 * t180;
t60 = -t115 * Ifges(6,1) + t107 + t323;
t59 = -t276 + t322 + t324;
t55 = -pkin(5) * t115 - pkin(8) * t114 - t237;
t49 = -pkin(5) * t116 - pkin(8) * t320 + t133;
t33 = t88 * Ifges(7,5) + t325 + t333;
t32 = t175 * t67 - t272;
t31 = t173 * t185 + t175 * t66;
t30 = t173 * t67 + t57;
t15 = t178 * t55 + t180 * t32;
t14 = -t178 * t32 + t180 * t55;
t10 = t40 * Ifges(7,1) + t41 * Ifges(7,4) + t74 * Ifges(7,5);
t9 = t40 * Ifges(7,4) + t41 * Ifges(7,2) + t74 * Ifges(7,6);
t5 = -qJDD(4) * pkin(5) - t7;
t4 = -qJD(6) * t20 - t178 * t31 + t180 * t49;
t3 = qJD(6) * t19 + t178 * t49 + t180 * t31;
t21 = [(Ifges(4,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) + t307 * qJD(2) + t308 * t125 - t10 * t262 / 0.2e1 - t130 * t257 + t118 * t252 / 0.2e1 - t119 * t251 / 0.2e1 + m(3) * (-pkin(1) * t157 + (t142 + t246) * qJ(2)) + m(7) * (t1 * t20 + t11 * t4 + t12 * t3 + t19 * t2) + m(6) * (t111 * t54 + t133 * t82 + t28 * t31 + t48 * t8) + (qJD(6) * t35 + t9) * t263 / 0.2e1 + (t1 * t263 - t11 * t197 + t12 * t198 + t2 * t262) * mrSges(7,3) + t87 * (Ifges(7,4) * t197 + Ifges(7,2) * t198) / 0.2e1 + m(4) * (t90 * t130 + t91 * t131) + t126 * t158 + (t251 * t78 + t252 * t79 - t319) * mrSges(5,3) + t321 * t223 + (qJD(1) * t221 + t91) * mrSges(4,2) + (-t199 + Ifges(6,5) * t320 / 0.2e1 + Ifges(6,6) * t116 / 0.2e1 + t207 * qJD(4) / 0.2e1) * qJD(4) + (t116 * t28 - t27 * t320) * mrSges(6,3) - (-Ifges(6,1) * t297 - t60 / 0.2e1 - t107 / 0.2e1 + t310) * t320 + t314 * t124 + (-Ifges(7,3) * t303 - t240 / 0.2e1 - Ifges(7,6) * t304 - Ifges(7,5) * t305 + t8 * mrSges(6,3) + Ifges(6,4) * t76 + Ifges(6,2) * t75 - t54 * mrSges(6,1) + t317 * Ifges(6,6) + t311) * t120 + (-t54 * mrSges(6,2) + t7 * mrSges(6,3) - Ifges(6,1) * t76 - Ifges(6,4) * t75 - Ifges(6,5) * t317 - t206 * t303 - t208 * t304 - t210 * t305 - t5 * t212 + t34 * t226) * t122 + (Ifges(7,1) * t197 + Ifges(7,4) * t198) * t300 - t340 * (t173 * t66 - t175 * t185) + (-m(3) * t217 + t331 * t294 + t332 * t293 + t335 * t191 + t334 * (-t121 * t177 + t123 * t155 + t191) + (-t159 * t339 + t312 + t338) * t123 + (-t212 + t309) * t121) * g(1) + t4 * t52 + t48 * t64 + (-m(6) * t7 + m(7) * t5 + t330) * (-t100 * t173 - t175 * t222) + (Ifges(6,4) * t297 - Ifges(7,5) * t300 + t59 / 0.2e1 - t33 / 0.2e1 + t324 / 0.2e1 - t333 / 0.2e1 - t325 / 0.2e1 - t306) * t116 - t90 * mrSges(4,1) + t31 * t98 + (-m(3) * t258 - t63 * mrSges(7,1) - t62 * mrSges(7,2) - t332 * t294 + t331 * t293 + t335 * t232 + t334 * (-t121 * t155 - t123 * t177 + t232) + t309 * t123 + (-m(7) * (-pkin(5) * t160 - t288) + t273 + t214 - t338) * t121) * g(2) + t111 * t42 + t133 * t69 + t84 * t139 + 0.2e1 * t142 * mrSges(3,3) + t19 * t17 + t20 * t18 + t3 * t51 - t157 * mrSges(3,1) - t179 * (Ifges(5,1) * t129 + Ifges(5,4) * t128 + Ifges(5,5) * qJDD(4)) / 0.2e1 - t181 * (Ifges(5,4) * t129 + Ifges(5,2) * t128 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t108 * (Ifges(7,5) * t197 + Ifges(7,6) * t198) / 0.2e1 + t25 * (-mrSges(7,1) * t198 + mrSges(7,2) * t197) + t131 * t247 + pkin(1) * t248 + (-Ifges(5,5) * t179 - Ifges(5,6) * t181) * t267 + t246 * t283 + t128 * t209 / 0.2e1 + t129 * t211 / 0.2e1; m(3) * t157 - t104 * t64 + t80 * t17 + t81 * t18 - t248 - t341 * t98 + t329 * t52 + t328 * t51 + (-t126 - t69) * t256 + t330 * t103 + (t1 * t81 + t103 * t5 + t11 * t329 + t12 * t328 + t2 * t80) * m(7) + (-t7 * t103 - t8 * t104 - t256 * t82 - t28 * t341) * m(6) + (-m(3) * qJ(2) - mrSges(3,3) - t283) * t183 - t307 * qJD(1) + (m(4) * t90 - m(6) * t54 - mrSges(4,2) * t183 - t257 - t314 - t42) * t266 + (m(4) * t91 + t247 + t308) * t174 + t340 * ((t173 * t224 + t175 * t225) * qJD(1) + t320 * t174) + (-g(1) * t293 + g(2) * t294) * (m(3) + t316); m(4) * qJDD(3) + t179 * t109 + t181 * t110 + t330 * t120 + t284 * t116 + (-t137 * t179 + t138 * t181) * qJD(4) + t201 * t320 + (t64 + (-t178 * t51 - t180 * t52) * qJD(6) + t318) * t122 + m(5) * (t179 * t38 + t181 * t39 + (-t179 * t78 + t181 * t79) * qJD(4)) + m(6) * (-t116 * t27 - t120 * t7 + t122 * t8 - t28 * t320) + m(7) * (t116 * t25 + t120 * t5 + t122 * t188 - t205 * t320) + t316 * g(3); (-t235 - t138) * t78 + (t60 + t107) * t298 + t207 * t223 + (-t236 + t137) * t79 - t118 * t255 / 0.2e1 + t119 * t254 / 0.2e1 - t34 * t250 / 0.2e1 + t69 * t237 + (t108 * t206 + t208 * t87 + t210 * t88) * qJD(6) / 0.2e1 + qJD(6) * t200 + (-m(7) * (-t168 - t288) + t273 + t139 + m(6) * t168 - t312) * g(3) + t65 * t290 + t64 * t291 + t10 * t295 + t59 * t297 + (m(7) * t188 - t249 * t52 - t250 * t51 + t318) * (pkin(8) + t291) + t321 * t183 / 0.2e1 + (-t11 * t249 - t12 * t250) * mrSges(7,3) + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t276 + t33) * t296 + (g(1) * t121 + g(2) * t123) * (-t215 + t334 * pkin(4) * t179 + (-mrSges(6,2) + t339) * t160 + (-mrSges(6,1) - t192) * t159) + (-t11 * t14 - t12 * t15 + t151 * t5 - t25 * t30) * m(7) + (Ifges(7,5) * t178 + Ifges(7,6) * t180) * t303 + (Ifges(7,2) * t180 + t275) * t304 + (Ifges(7,1) * t178 + t274) * t305 + t35 * t226 - t14 * t52 - t284 * t30 + Ifges(6,6) * t75 + Ifges(6,5) * t76 - t32 * t98 + Ifges(5,6) * t128 + Ifges(5,5) * t129 + t7 * mrSges(6,1) - t8 * mrSges(6,2) - t2 * t280 - t28 * t281 - t38 * mrSges(5,2) + t39 * mrSges(5,1) - t15 * t51 + t151 * t13 + t180 * t9 / 0.2e1 + qJD(1) * t199 + (Ifges(6,2) * t298 - Ifges(7,3) * t299 - Ifges(7,5) * t301 - Ifges(7,6) * t302 - t322 / 0.2e1 + t306) * t115 + (Ifges(6,1) * t296 + t206 * t299 + t210 * t301 + t208 * t302 - t323 / 0.2e1 + t11 * t279 + t12 * t280 - t200 + t310) * t114 + t1 * t279 + t27 * t282 + ((t173 * t8 + t175 * t7) * pkin(4) + t237 * t82 + t27 * t30 - t28 * t32) * m(6) + t5 * t213; t204 * qJD(6) + t201 * t114 + t284 * t115 + t180 * t17 + t178 * t18 + t42 + (t1 * t178 + t108 * t205 + t115 * t25 + t180 * t2 + t216) * m(7) + (-t114 * t28 - t115 * t27 + t216 + t54) * m(6); -t25 * (mrSges(7,1) * t88 + mrSges(7,2) * t87) + (Ifges(7,1) * t87 - t292) * t301 + t34 * t300 + (Ifges(7,5) * t87 - Ifges(7,6) * t88) * t299 - t11 * t51 + t12 * t52 - g(1) * (mrSges(7,1) * t62 - mrSges(7,2) * t63) - g(2) * ((-t121 * t180 + t123 * t261) * mrSges(7,1) + (t121 * t178 + t123 * t260) * mrSges(7,2)) - g(3) * t212 * t159 + (t11 * t87 + t12 * t88) * mrSges(7,3) + t240 + (-Ifges(7,2) * t88 + t35 + t86) * t302 - t311;];
tau  = t21;
