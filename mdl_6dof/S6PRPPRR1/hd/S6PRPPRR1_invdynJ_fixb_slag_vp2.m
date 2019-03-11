% Calculate vector of inverse dynamics joint torques for
% S6PRPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
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
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:14:00
% EndTime: 2019-03-08 19:14:19
% DurationCPUTime: 10.95s
% Computational Cost: add. (5444->493), mult. (12755->683), div. (0->0), fcn. (10352->16), ass. (0->226)
t183 = sin(qJ(5));
t186 = cos(qJ(5));
t176 = sin(pkin(6));
t187 = cos(qJ(2));
t244 = qJD(1) * t187;
t226 = t176 * t244;
t147 = qJD(2) * pkin(2) + t226;
t174 = sin(pkin(11));
t178 = cos(pkin(11));
t184 = sin(qJ(2));
t253 = t176 * t184;
t227 = qJD(1) * t253;
t108 = t174 * t147 + t178 * t227;
t105 = qJD(2) * qJ(4) + t108;
t180 = cos(pkin(6));
t160 = qJD(1) * t180 + qJD(3);
t177 = cos(pkin(12));
t151 = t177 * t160;
t173 = sin(pkin(12));
t275 = pkin(8) * qJD(2);
t64 = t151 + (-t105 - t275) * t173;
t73 = t177 * t105 + t173 * t160;
t65 = t177 * t275 + t73;
t25 = -t183 * t65 + t186 * t64;
t21 = -qJD(5) * pkin(5) - t25;
t143 = t173 * t186 + t177 * t183;
t134 = t143 * qJD(2);
t182 = sin(qJ(6));
t185 = cos(qJ(6));
t110 = qJD(5) * t185 - t134 * t182;
t111 = qJD(5) * t182 + t134 * t185;
t279 = mrSges(6,3) * t134;
t265 = qJD(5) * mrSges(6,1) + mrSges(7,1) * t110 - mrSges(7,2) * t111 - t279;
t331 = -m(7) * t21 + t265;
t153 = t174 * t227;
t120 = t178 * t226 - t153;
t326 = -t120 + qJD(4);
t215 = -mrSges(7,1) * t185 + mrSges(7,2) * t182;
t329 = m(7) * pkin(5) + mrSges(6,1) - t215;
t313 = -m(7) * pkin(9) + mrSges(6,2) - mrSges(7,3);
t315 = mrSges(5,3) * (t173 ^ 2 + t177 ^ 2);
t141 = t173 * t183 - t186 * t177;
t162 = pkin(2) * t174 + qJ(4);
t282 = pkin(8) + t162;
t137 = t282 * t173;
t138 = t282 * t177;
t205 = -t186 * t137 - t138 * t183;
t316 = qJD(5) * t205 - t326 * t141;
t204 = t174 * t187 + t178 * t184;
t127 = t204 * t176;
t117 = qJD(1) * t127;
t135 = t141 * qJD(5);
t136 = t143 * qJD(5);
t328 = pkin(5) * t136 + pkin(9) * t135 - t117;
t322 = -m(6) - m(7);
t323 = -m(4) - m(5);
t327 = t322 + t323;
t172 = pkin(12) + qJ(5);
t168 = sin(t172);
t169 = cos(t172);
t217 = -mrSges(5,1) * t177 + mrSges(5,2) * t173;
t325 = m(5) * pkin(3) - t168 * t313 + t169 * t329 - t217;
t175 = sin(pkin(10));
t179 = cos(pkin(10));
t248 = t180 * t187;
t324 = -t175 * t184 + t179 * t248;
t97 = -qJD(2) * t136 - qJDD(2) * t141;
t133 = t141 * qJD(2);
t132 = qJD(6) + t133;
t96 = -qJD(2) * t135 + qJDD(2) * t143;
t44 = qJD(6) * t110 + qJDD(5) * t182 + t185 * t96;
t297 = t44 / 0.2e1;
t45 = -qJD(6) * t111 + qJDD(5) * t185 - t182 * t96;
t296 = t45 / 0.2e1;
t93 = qJDD(6) - t97;
t295 = t93 / 0.2e1;
t26 = t183 * t64 + t186 * t65;
t22 = qJD(5) * pkin(9) + t26;
t166 = pkin(4) * t177 + pkin(3);
t107 = t147 * t178 - t153;
t219 = qJD(4) - t107;
t90 = -qJD(2) * t166 + t219;
t43 = pkin(5) * t133 - pkin(9) * t134 + t90;
t9 = -t182 * t22 + t185 * t43;
t321 = t9 * mrSges(7,1);
t10 = t182 * t43 + t185 * t22;
t320 = t10 * mrSges(7,2);
t286 = pkin(2) * t178;
t152 = -t166 - t286;
t82 = pkin(5) * t141 - pkin(9) * t143 + t152;
t88 = -t137 * t183 + t138 * t186;
t31 = -t182 * t88 + t185 * t82;
t319 = qJD(6) * t31 + t328 * t182 + t185 * t316;
t32 = t182 * t82 + t185 * t88;
t318 = -qJD(6) * t32 - t182 * t316 + t328 * t185;
t17 = -mrSges(7,1) * t45 + mrSges(7,2) * t44;
t281 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t96 + t17;
t312 = mrSges(6,1) * t133 + mrSges(6,2) * t134 + t217 * qJD(2);
t247 = t187 * t178;
t126 = t174 * t253 - t176 * t247;
t116 = t126 * t185;
t103 = -t127 * t173 + t177 * t180;
t104 = t127 * t177 + t173 * t180;
t49 = t103 * t183 + t104 * t186;
t27 = -t182 * t49 + t116;
t310 = (mrSges(3,1) * t187 - mrSges(3,2) * t184) * t176 - t126 * mrSges(4,1) - t127 * mrSges(4,2);
t241 = qJD(6) * t185;
t199 = -t135 * t182 + t143 * t241;
t23 = mrSges(7,1) * t93 - mrSges(7,3) * t44;
t24 = -mrSges(7,2) * t93 + mrSges(7,3) * t45;
t309 = -t182 * t23 + t185 * t24;
t159 = qJDD(1) * t180 + qJDD(3);
t149 = t177 * t159;
t240 = qJD(2) * qJD(4);
t222 = qJD(2) * t227;
t252 = t176 * t187;
t130 = qJDD(1) * t252 - t222;
t264 = qJDD(2) * pkin(2);
t125 = t130 + t264;
t224 = qJD(2) * t244;
t131 = (qJDD(1) * t184 + t224) * t176;
t69 = t174 * t125 + t178 * t131;
t62 = qJDD(2) * qJ(4) + t240 + t69;
t46 = -t173 * t62 + t149;
t47 = t173 * t159 + t177 * t62;
t308 = -t173 * t46 + t177 * t47;
t68 = t125 * t178 - t174 * t131;
t207 = qJDD(4) - t68;
t54 = -qJDD(2) * t166 + t207;
t18 = -pkin(5) * t97 - pkin(9) * t96 + t54;
t35 = t149 + (-pkin(8) * qJDD(2) - t62) * t173;
t238 = qJDD(2) * t177;
t36 = pkin(8) * t238 + t47;
t5 = qJD(5) * t25 + t183 * t35 + t186 * t36;
t3 = qJDD(5) * pkin(9) + t5;
t1 = qJD(6) * t9 + t18 * t182 + t185 * t3;
t2 = -qJD(6) * t10 + t18 * t185 - t182 * t3;
t307 = t1 * t185 - t182 * t2;
t6 = -qJD(5) * t26 - t183 * t36 + t186 * t35;
t142 = t174 * t184 - t247;
t249 = t180 * t184;
t246 = -t174 * t248 - t178 * t249;
t81 = -t179 * t142 + t175 * t246;
t76 = t175 * t142 + t179 * t246;
t303 = -mrSges(4,1) - t325;
t302 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t214 = mrSges(7,1) * t182 + mrSges(7,2) * t185;
t301 = -m(5) * qJ(4) - mrSges(5,3) - mrSges(6,3) - t214;
t300 = -mrSges(4,2) - t301;
t299 = Ifges(7,1) * t297 + Ifges(7,4) * t296 + Ifges(7,5) * t295;
t294 = -t110 / 0.2e1;
t293 = -t111 / 0.2e1;
t292 = t111 / 0.2e1;
t291 = -t132 / 0.2e1;
t288 = t134 / 0.2e1;
t287 = t185 / 0.2e1;
t280 = mrSges(6,3) * t133;
t278 = Ifges(6,4) * t134;
t277 = Ifges(7,4) * t182;
t276 = Ifges(7,4) * t185;
t274 = t111 * Ifges(7,4);
t239 = qJDD(2) * t173;
t140 = -mrSges(5,1) * t238 + mrSges(5,2) * t239;
t37 = -t97 * mrSges(6,1) + t96 * mrSges(6,2);
t266 = t140 + t37;
t263 = t126 * t182;
t262 = t133 * t182;
t261 = t133 * t185;
t259 = t143 * t182;
t258 = t143 * t185;
t256 = t175 * t176;
t254 = t176 * t179;
t243 = qJD(2) * t120;
t242 = qJD(6) * t182;
t237 = Ifges(7,5) * t44 + Ifges(7,6) * t45 + Ifges(7,3) * t93;
t161 = pkin(2) * t252;
t106 = Ifges(7,4) * t110;
t40 = Ifges(7,1) * t111 + Ifges(7,5) * t132 + t106;
t231 = t40 * t287;
t223 = -t242 / 0.2e1;
t221 = t324 * pkin(2);
t220 = t10 * t185 - t9 * t182;
t213 = Ifges(7,1) * t185 - t277;
t212 = -Ifges(7,2) * t182 + t276;
t211 = Ifges(7,5) * t185 - Ifges(7,6) * t182;
t210 = -t173 * (-t105 * t173 + t151) + t177 * t73;
t70 = -mrSges(7,2) * t132 + mrSges(7,3) * t110;
t71 = mrSges(7,1) * t132 - mrSges(7,3) * t111;
t209 = -t182 * t71 + t185 * t70;
t28 = t185 * t49 + t263;
t206 = t186 * t103 - t104 * t183;
t121 = -qJD(5) * mrSges(6,2) - t280;
t203 = -t121 - t209;
t200 = t21 * t214;
t198 = t135 * t185 + t143 * t242;
t195 = t142 * t180;
t77 = -t175 * t204 - t179 * t195;
t80 = t175 * t195 - t179 * t204;
t196 = g(1) * t80 + g(2) * t77 - g(3) * t126;
t190 = (-t10 * t182 - t185 * t9) * qJD(6) + t307;
t188 = qJD(2) ^ 2;
t181 = -pkin(8) - qJ(4);
t167 = -pkin(3) - t286;
t129 = Ifges(6,4) * t133;
t119 = t142 * t176 * qJD(2);
t118 = qJD(2) * t127;
t102 = -qJD(2) * pkin(3) + t219;
t99 = t127 * t169 + t168 * t180;
t94 = pkin(5) * t134 + pkin(9) * t133;
t86 = t134 * Ifges(6,1) + Ifges(6,5) * qJD(5) - t129;
t85 = -t133 * Ifges(6,2) + Ifges(6,6) * qJD(5) + t278;
t84 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t97;
t63 = -qJDD(2) * pkin(3) + t207;
t53 = t168 * t256 + t169 * t81;
t51 = -t168 * t254 - t169 * t76;
t39 = t110 * Ifges(7,2) + t132 * Ifges(7,6) + t274;
t38 = t111 * Ifges(7,5) + t110 * Ifges(7,6) + t132 * Ifges(7,3);
t20 = qJD(5) * t49 - t119 * t143;
t19 = qJD(5) * t206 + t119 * t141;
t16 = t182 * t94 + t185 * t25;
t15 = -t182 * t25 + t185 * t94;
t13 = t44 * Ifges(7,4) + t45 * Ifges(7,2) + t93 * Ifges(7,6);
t8 = -qJD(6) * t28 + t118 * t185 - t182 * t19;
t7 = qJD(6) * t27 + t118 * t182 + t185 * t19;
t4 = -qJDD(5) * pkin(5) - t6;
t11 = [m(2) * qJDD(1) + t19 * t121 + t27 * t23 + t28 * t24 + t49 * t84 + t7 * t70 + t8 * t71 - t281 * t206 - t265 * t20 + (-mrSges(3,1) * t184 - mrSges(3,2) * t187) * t188 * t176 + t266 * t126 + t312 * t118 + (-t118 * mrSges(4,1) - (-mrSges(4,2) + t315) * t119) * qJD(2) + ((-t103 * t173 + t104 * t177) * mrSges(5,3) + t310) * qJDD(2) + (-m(2) - m(3) + t327) * g(3) + m(7) * (t1 * t28 + t10 * t7 + t2 * t27 + t20 * t21 - t206 * t4 + t8 * t9) + m(6) * (t118 * t90 + t126 * t54 + t19 * t26 - t20 * t25 + t206 * t6 + t49 * t5) + m(5) * (t102 * t118 + t103 * t46 + t104 * t47 - t119 * t210 + t126 * t63) + m(4) * (-t107 * t118 - t108 * t119 - t126 * t68 + t127 * t69 + t159 * t180) + m(3) * (qJDD(1) * t180 ^ 2 + (t130 * t187 + t131 * t184) * t176); (-t324 * mrSges(3,1) - (-t175 * t187 - t179 * t249) * mrSges(3,2) + t323 * t221 + t322 * (t77 * t166 + t76 * t181 + t221) + t303 * t77 + t300 * t76) * g(2) + (t222 + t130) * mrSges(3,1) - t281 * t205 + t308 * mrSges(5,3) - t199 * t39 / 0.2e1 + (t25 * t135 - t26 * t136) * mrSges(6,3) + qJD(5) * (-Ifges(6,5) * t135 - Ifges(6,6) * t136) / 0.2e1 + t21 * (mrSges(7,1) * t199 - mrSges(7,2) * t198) + t132 * (-Ifges(7,5) * t198 - Ifges(7,6) * t199 + Ifges(7,3) * t136) / 0.2e1 + t110 * (-Ifges(7,4) * t198 - Ifges(7,2) * t199 + Ifges(7,6) * t136) / 0.2e1 + (qJD(2) * t117 + t178 * t264 + t68) * mrSges(4,1) + (m(6) * t25 + t331) * (-qJD(5) * t88 - t143 * t326) + t316 * t121 + t318 * t71 + t319 * t70 + (t176 * t224 - t131) * mrSges(3,2) + t63 * t217 - t13 * t259 / 0.2e1 + (t237 / 0.2e1 - Ifges(6,4) * t96 - Ifges(6,2) * t97 - Ifges(6,6) * qJDD(5) + t54 * mrSges(6,1) - t5 * mrSges(6,3) + Ifges(7,3) * t295 + Ifges(7,6) * t296 + Ifges(7,5) * t297 + t302) * t141 - t312 * t117 + (-t1 * t259 - t10 * t199 + t9 * t198 - t2 * t258) * mrSges(7,3) - t136 * t320 + (t54 * mrSges(6,2) - t6 * mrSges(6,3) + Ifges(6,1) * t96 + Ifges(6,4) * t97 + Ifges(6,5) * qJDD(5) + t211 * t295 + t212 * t296 + t213 * t297 + t214 * t4 + t223 * t40) * t143 + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t322 * (-t126 * t166 - t127 * t181 + t161) + t323 * t161 + t301 * t127 + t325 * t126 - t310) * g(3) + (t1 * t32 + t10 * t319 + t2 * t31 - t205 * t4 + t318 * t9) * m(7) + ((t174 * t69 + t178 * t68) * pkin(2) + t107 * t117 - t108 * t120) * m(4) + (-t117 * t90 + t152 * t54 + t205 * t6 + t26 * t316 + t5 * t88) * m(6) - t133 * (-Ifges(6,4) * t135 - Ifges(6,2) * t136) / 0.2e1 + t90 * (mrSges(6,1) * t136 - mrSges(6,2) * t135) + (-Ifges(6,1) * t135 - Ifges(6,4) * t136) * t288 - t135 * t231 + (-t174 * t264 + t243 - t69) * mrSges(4,2) + t136 * t321 + t167 * t140 + t152 * t37 + t136 * t38 / 0.2e1 - t136 * t85 / 0.2e1 - t135 * t86 / 0.2e1 + t88 * t84 + (-t102 * t117 + t308 * t162 + t167 * t63 + t210 * t326) * m(5) + (-(t175 * t249 - t179 * t187) * mrSges(3,2) + t322 * (t80 * t166 - t181 * t81) + t303 * t80 - t300 * t81 + (pkin(2) * t327 - mrSges(3,1)) * (-t175 * t248 - t179 * t184)) * g(1) + (Ifges(5,4) * t173 + Ifges(5,2) * t177) * t238 + (Ifges(5,1) * t173 + Ifges(5,4) * t177) * t239 + (-Ifges(7,1) * t198 - Ifges(7,4) * t199 + Ifges(7,5) * t136) * t292 + t258 * t299 + (qJDD(2) * t162 + t240 - t243) * t315 + t31 * t23 + t32 * t24; t281 * t141 - t265 * t136 + t203 * t135 + (t84 + (-t182 * t70 - t185 * t71) * qJD(6) + t309) * t143 + m(6) * (-t135 * t26 - t136 * t25 - t141 * t6 + t143 * t5) + m(4) * t159 + m(5) * (t173 * t47 + t177 * t46) + m(7) * (-t135 * t220 + t136 * t21 + t141 * t4 + t143 * t190) - (-t180 * g(3) + (-g(1) * t175 + g(2) * t179) * t176) * t327; -t188 * t315 + t209 * qJD(6) - t203 * t133 + t265 * t134 + t182 * t24 + t185 * t23 + t266 + (t1 * t182 + t132 * t220 - t134 * t21 + t2 * t185 + t196) * m(7) + (t133 * t26 + t134 * t25 + t196 + t54) * m(6) + (-qJD(2) * t210 + t196 + t63) * m(5); (-pkin(5) * t4 - t10 * t16 - t15 * t9) * m(7) + (-Ifges(6,2) * t134 - t129 + t86) * t133 / 0.2e1 + (t231 + t200) * qJD(6) + ((-t241 - t261) * t9 + (-t242 - t262) * t10 + t307) * mrSges(7,3) + (m(7) * t190 - t241 * t71 - t242 * t70 + t309) * pkin(9) - (-Ifges(6,1) * t133 - t278 + t38) * t134 / 0.2e1 - qJD(5) * (-Ifges(6,5) * t133 - Ifges(6,6) * t134) / 0.2e1 - t90 * (mrSges(6,1) * t134 - mrSges(6,2) * t133) + (Ifges(7,3) * t134 - t133 * t211) * t291 + (Ifges(7,5) * t134 - t133 * t213) * t293 + (Ifges(7,6) * t134 - t133 * t212) * t294 + t133 * t200 + (t279 + t331) * t26 + t4 * t215 + t40 * t261 / 0.2e1 - t134 * t321 + (t110 * t212 + t111 * t213 + t132 * t211) * qJD(6) / 0.2e1 + (-t262 / 0.2e1 + t223) * t39 + (-t121 - t280) * t25 + (t313 * t53 - t329 * (-t168 * t81 + t169 * t256)) * g(1) + (t313 * t99 - t329 * (-t127 * t168 + t169 * t180)) * g(3) + (t313 * t51 - t329 * (t168 * t76 - t169 * t254)) * g(2) + Ifges(6,3) * qJDD(5) + Ifges(6,5) * t96 + Ifges(6,6) * t97 + t13 * t287 + t85 * t288 + (Ifges(7,5) * t182 + Ifges(7,6) * t185) * t295 + (Ifges(7,2) * t185 + t277) * t296 + (Ifges(7,1) * t182 + t276) * t297 + t182 * t299 + t134 * t320 - t5 * mrSges(6,2) + t6 * mrSges(6,1) - pkin(5) * t17 - t16 * t70 - t15 * t71; -t21 * (mrSges(7,1) * t111 + mrSges(7,2) * t110) + (Ifges(7,1) * t110 - t274) * t293 + t39 * t292 + (Ifges(7,5) * t110 - Ifges(7,6) * t111) * t291 - t9 * t70 + t10 * t71 - g(1) * ((-t182 * t53 - t185 * t80) * mrSges(7,1) + (t182 * t80 - t185 * t53) * mrSges(7,2)) - g(2) * ((-t182 * t51 - t185 * t77) * mrSges(7,1) + (t182 * t77 - t185 * t51) * mrSges(7,2)) - g(3) * ((-t182 * t99 + t116) * mrSges(7,1) + (-t185 * t99 - t263) * mrSges(7,2)) + (t10 * t111 + t110 * t9) * mrSges(7,3) + t237 + (-Ifges(7,2) * t111 + t106 + t40) * t294 + t302;];
tau  = t11;
