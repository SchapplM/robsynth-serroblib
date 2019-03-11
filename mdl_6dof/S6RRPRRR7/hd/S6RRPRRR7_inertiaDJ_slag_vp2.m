% Calculate time derivative of joint inertia matrix for
% S6RRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:46
% EndTime: 2019-03-09 13:55:59
% DurationCPUTime: 5.50s
% Computational Cost: add. (7097->507), mult. (14904->710), div. (0->0), fcn. (13412->8), ass. (0->205)
t152 = sin(qJ(5));
t156 = cos(qJ(5));
t174 = mrSges(6,1) * t152 + mrSges(6,2) * t156;
t111 = t174 * qJD(5);
t151 = sin(qJ(6));
t155 = cos(qJ(6));
t109 = t151 * t156 + t152 * t155;
t257 = qJD(5) + qJD(6);
t269 = t257 * t109;
t195 = qJD(6) * t155;
t197 = qJD(5) * t156;
t211 = t151 * t152;
t78 = -t155 * t197 - t156 * t195 + t211 * t257;
t39 = mrSges(7,1) * t269 - mrSges(7,2) * t78;
t276 = t111 + t39;
t153 = sin(qJ(4));
t202 = qJD(4) * t153;
t154 = sin(qJ(2));
t240 = pkin(7) - pkin(8);
t128 = t240 * t154;
t158 = cos(qJ(2));
t130 = t240 * t158;
t157 = cos(qJ(4));
t259 = t157 * t128 - t130 * t153;
t216 = qJD(4) * t259;
t106 = -t155 * t156 + t211;
t107 = t153 * t154 + t157 * t158;
t108 = -t153 * t158 + t154 * t157;
t213 = t108 * t156;
t258 = -t158 * pkin(2) - t154 * qJ(3);
t122 = -pkin(1) + t258;
t104 = t158 * pkin(3) - t122;
t64 = pkin(4) * t107 - pkin(9) * t108 + t104;
t92 = t128 * t153 + t130 * t157;
t36 = -t152 * t92 + t156 * t64;
t25 = pkin(5) * t107 - pkin(10) * t213 + t36;
t214 = t108 * t152;
t88 = t156 * t92;
t37 = t152 * t64 + t88;
t29 = -pkin(10) * t214 + t37;
t11 = -t151 * t29 + t155 * t25;
t12 = t151 * t25 + t155 * t29;
t124 = Ifges(6,5) * t152 + Ifges(6,6) * t156;
t198 = qJD(5) * t152;
t138 = Ifges(6,6) * t198;
t201 = qJD(4) * t157;
t204 = qJD(2) * t158;
t205 = qJD(2) * t154;
t80 = t153 * t204 + t154 * t201 - t157 * t205 - t158 * t202;
t81 = (qJD(2) - qJD(4)) * t107;
t159 = -pkin(2) - pkin(3);
t207 = qJ(3) * t204 + t154 * qJD(3);
t94 = t159 * t205 + t207;
t32 = pkin(4) * t80 - pkin(9) * t81 + t94;
t119 = t240 * t205;
t179 = qJD(2) * t130;
t51 = -t157 * t119 + t153 * t179 + t216;
t185 = -t152 * t51 + t156 * t32;
t219 = t156 * t81;
t6 = -pkin(10) * t219 + pkin(5) * t80 + (-t88 + (pkin(10) * t108 - t64) * t152) * qJD(5) + t185;
t220 = t152 * t81;
t167 = t108 * t197 + t220;
t9 = t152 * t32 + t156 * t51 + t64 * t197 - t198 * t92;
t8 = -pkin(10) * t167 + t9;
t2 = qJD(6) * t11 + t151 * t6 + t155 * t8;
t225 = Ifges(7,5) * t78 + Ifges(7,6) * t269;
t215 = qJD(4) * t92;
t52 = -t119 * t153 - t157 * t179 + t215;
t24 = pkin(5) * t167 + t52;
t66 = t109 * t108;
t67 = t106 * t108;
t27 = -Ifges(7,4) * t67 - Ifges(7,2) * t66 + Ifges(7,6) * t107;
t28 = -Ifges(7,1) * t67 - Ifges(7,4) * t66 + Ifges(7,5) * t107;
t3 = -qJD(6) * t12 - t151 * t8 + t155 * t6;
t17 = -t106 * t81 - t108 * t269;
t161 = t257 * t106;
t18 = t108 * t161 - t109 * t81;
t4 = Ifges(7,4) * t17 + Ifges(7,2) * t18 + t80 * Ifges(7,6);
t40 = Ifges(7,4) * t78 + Ifges(7,2) * t269;
t42 = Ifges(7,1) * t78 + Ifges(7,4) * t269;
t5 = Ifges(7,1) * t17 + Ifges(7,4) * t18 + t80 * Ifges(7,5);
t65 = pkin(5) * t214 - t259;
t82 = mrSges(7,1) * t106 + mrSges(7,2) * t109;
t275 = t24 * t82 + t67 * t42 / 0.2e1 + t66 * t40 / 0.2e1 + t65 * t39 - (t106 * t2 + t109 * t3 - t11 * t78 + t12 * t269) * mrSges(7,3) - (Ifges(5,6) - Ifges(7,5) * t109 / 0.2e1 + Ifges(7,6) * t106 / 0.2e1 - t124 / 0.2e1) * t80 + Ifges(5,5) * t81 - t259 * t111 - t106 * t4 / 0.2e1 + t109 * t5 / 0.2e1 - t51 * mrSges(5,2) - t78 * t28 / 0.2e1 - t269 * t27 / 0.2e1 + (Ifges(6,5) * t197 - t138 - t225) * t107 / 0.2e1;
t84 = Ifges(7,4) * t109 - Ifges(7,2) * t106;
t242 = t84 / 0.2e1;
t86 = Ifges(7,1) * t109 - Ifges(7,4) * t106;
t241 = t86 / 0.2e1;
t191 = pkin(5) * t198;
t121 = t157 * qJ(3) + t153 * t159;
t98 = t153 * qJD(3) + qJD(4) * t121;
t93 = t98 - t191;
t274 = t93 * t82;
t120 = -t153 * qJ(3) + t157 * t159;
t115 = pkin(4) - t120;
t231 = t156 * pkin(5);
t102 = t115 + t231;
t273 = t102 * t39;
t137 = -pkin(4) - t231;
t272 = t137 * t39;
t271 = (t151 * t269 - t155 * t78 + (t106 * t155 - t109 * t151) * qJD(6)) * mrSges(7,3);
t206 = t152 ^ 2 + t156 ^ 2;
t183 = t206 * mrSges(6,3);
t270 = mrSges(5,2) - t183;
t268 = t106 * t40 - t109 * t42;
t267 = -t269 * t84 - t78 * t86 + t268;
t265 = Ifges(6,5) * t219 + Ifges(6,3) * t80;
t100 = t106 * t153;
t46 = -t106 * t201 - t153 * t269;
t47 = -t109 * t201 + t153 * t161;
t99 = t109 * t153;
t264 = (-t100 * t269 + t106 * t46 + t109 * t47 + t78 * t99) * mrSges(7,3);
t123 = -t156 * mrSges(6,1) + mrSges(6,2) * t152;
t217 = mrSges(5,1) - t123;
t263 = t217 * t98;
t262 = m(4) * pkin(7) + mrSges(4,2);
t35 = mrSges(7,1) * t66 - mrSges(7,2) * t67;
t261 = m(7) * t65 + t35;
t260 = qJD(5) * t37;
t256 = -m(6) * pkin(4) - t217;
t253 = 2 * m(5);
t252 = 0.2e1 * m(6);
t251 = 0.2e1 * m(7);
t250 = -0.2e1 * pkin(1);
t249 = 0.2e1 * t52;
t248 = 0.2e1 * t94;
t247 = -0.2e1 * t111;
t246 = 0.2e1 * t122;
t239 = -pkin(10) - pkin(9);
t236 = -t108 / 0.2e1;
t223 = Ifges(6,4) * t152;
t125 = Ifges(6,2) * t156 + t223;
t235 = -t125 / 0.2e1;
t230 = t52 * t259;
t229 = t259 * t98;
t228 = qJD(5) / 0.2e1;
t227 = Ifges(4,5) - Ifges(3,4);
t116 = -pkin(9) + t121;
t226 = pkin(10) - t116;
t222 = Ifges(6,4) * t156;
t221 = Ifges(6,5) * t156;
t218 = t157 * t98;
t114 = Ifges(6,1) * t197 - Ifges(6,4) * t198;
t210 = t152 * t114;
t126 = Ifges(6,1) * t152 + t222;
t209 = t156 * t126;
t200 = qJD(5) * t116;
t199 = qJD(5) * t126;
t196 = qJD(6) * t151;
t194 = 0.2e1 * mrSges(7,3);
t193 = 0.2e1 * t158;
t192 = Ifges(7,5) * t17 + Ifges(7,6) * t18 + Ifges(7,3) * t80;
t190 = qJD(5) * t239;
t189 = t108 * t198;
t188 = t47 * mrSges(7,1) - t46 * mrSges(7,2);
t187 = -Ifges(6,6) * t152 - (2 * Ifges(5,4));
t186 = qJD(5) * t36 - t9;
t184 = qJD(5) * t226;
t97 = t157 * qJD(3) + qJD(4) * t120;
t182 = t206 * t97;
t10 = t185 - t260;
t181 = t10 + t260;
t180 = t116 * t206;
t166 = t189 - t219;
t22 = -Ifges(6,4) * t166 - Ifges(6,2) * t167 + Ifges(6,6) * t80;
t178 = t22 / 0.2e1 + t81 * t126 / 0.2e1;
t23 = -Ifges(6,1) * t166 - Ifges(6,4) * t167 + Ifges(6,5) * t80;
t177 = t23 / 0.2e1 + t81 * t235;
t175 = -t158 * mrSges(4,1) - t154 * mrSges(4,3);
t95 = t226 * t152;
t96 = t226 * t156;
t57 = t151 * t96 + t155 * t95;
t58 = t151 * t95 - t155 * t96;
t62 = t152 * t184 + t156 * t97;
t63 = -t152 * t97 + t156 * t184;
t20 = qJD(6) * t57 + t151 * t63 + t155 * t62;
t21 = -qJD(6) * t58 - t151 * t62 + t155 * t63;
t172 = t21 * mrSges(7,1) - t20 * mrSges(7,2) + t225;
t117 = t152 * t190;
t118 = t156 * t190;
t127 = t239 * t152;
t129 = t239 * t156;
t89 = t127 * t155 + t129 * t151;
t49 = qJD(6) * t89 + t117 * t155 + t118 * t151;
t91 = t127 * t151 - t129 * t155;
t50 = -qJD(6) * t91 - t117 * t151 + t118 * t155;
t171 = t50 * mrSges(7,1) - t49 * mrSges(7,2) - t225;
t113 = Ifges(6,4) * t197 - Ifges(6,2) * t198;
t169 = -t156 * t113 - t210;
t168 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t192;
t105 = (-mrSges(7,1) * t151 - mrSges(7,2) * t155) * qJD(6) * pkin(5);
t71 = mrSges(6,1) * t107 - mrSges(6,3) * t213;
t70 = -mrSges(6,2) * t107 - mrSges(6,3) * t214;
t68 = t174 * t108;
t56 = Ifges(6,5) * t107 + (Ifges(6,1) * t156 - t223) * t108;
t55 = Ifges(6,6) * t107 + (-Ifges(6,2) * t152 + t222) * t108;
t54 = mrSges(7,1) * t107 + mrSges(7,3) * t67;
t53 = -mrSges(7,2) * t107 - mrSges(7,3) * t66;
t34 = -mrSges(6,2) * t80 - mrSges(6,3) * t167;
t33 = mrSges(6,1) * t80 + mrSges(6,3) * t166;
t26 = mrSges(6,1) * t167 - mrSges(6,2) * t166;
t14 = -mrSges(7,2) * t80 + mrSges(7,3) * t18;
t13 = mrSges(7,1) * t80 - mrSges(7,3) * t17;
t7 = -mrSges(7,1) * t18 + mrSges(7,2) * t17;
t1 = [(mrSges(5,1) * t248 - 0.2e1 * t51 * mrSges(5,3) + t187 * t81 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t80 + t192 + t265) * t107 + 0.2e1 * (-t259 * t81 - t80 * t92) * mrSges(5,3) - 0.2e1 * t259 * t26 + ((mrSges(3,2) * t250 - 0.2e1 * t122 * mrSges(4,3) - t193 * t227) * t158 + (mrSges(3,1) * t250 + mrSges(4,1) * t246 + 0.2e1 * t227 * t154 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t193) * t154) * qJD(2) + t80 * (-Ifges(7,5) * t67 - Ifges(7,6) * t66) + (m(4) * t246 + 0.2e1 * t175) * (pkin(2) * t205 - t207) + (t10 * t36 + t37 * t9 - t230) * t252 + (t104 * t94 + t51 * t92 - t230) * t253 + (t11 * t3 + t12 * t2 + t24 * t65) * t251 + t68 * t249 + t56 * t219 + 0.2e1 * t11 * t13 + 0.2e1 * t12 * t14 + t18 * t27 + t17 * t28 - t55 * t220 + 0.2e1 * t24 * t35 + 0.2e1 * t36 * t33 + 0.2e1 * t37 * t34 + 0.2e1 * t2 * t53 + 0.2e1 * t3 * t54 + 0.2e1 * t65 * t7 - t66 * t4 - t67 * t5 + 0.2e1 * t9 * t70 + 0.2e1 * t10 * t71 + (mrSges(5,2) * t248 + mrSges(5,3) * t249 + 0.2e1 * Ifges(5,1) * t81 - t152 * t22 + t156 * t23 + (t187 + t221) * t80 + (-t107 * t124 - t152 * t56 - t156 * t55) * qJD(5)) * t108 + 0.2e1 * t104 * (mrSges(5,1) * t80 + mrSges(5,2) * t81); ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t158 + (-mrSges(4,2) * qJ(3) - Ifges(3,6) + Ifges(4,6)) * t154 + (m(4) * t258 - t158 * mrSges(3,1) + t154 * mrSges(3,2) + t175) * pkin(7)) * qJD(2) + t262 * qJD(3) * t158 + m(6) * (t115 * t52 - t229) + m(5) * (-t120 * t52 + t121 * t51 + t92 * t97 - t229) + (-t107 * t97 + t108 * t98 - t120 * t81 - t121 * t80) * mrSges(5,3) + m(7) * (t102 * t24 + t11 * t21 + t12 * t20 + t2 * t58 + t3 * t57 + t65 * t93) + t217 * t52 + (-qJD(5) * t56 / 0.2e1 + m(6) * (t116 * t9 - t200 * t36 + t37 * t97) + t116 * t34 + t97 * t70 - t71 * t200 + (-t114 / 0.2e1 + t125 * t228) * t108 + t186 * mrSges(6,3) - t178) * t156 + (m(6) * (-t10 * t116 - t200 * t37 - t36 * t97) - t116 * t33 - t97 * t71 + t55 * t228 - t70 * t200 + (t113 / 0.2e1 + t199 / 0.2e1) * t108 + t181 * mrSges(6,3) - t177) * t152 + t20 * t53 + t21 * t54 + t57 * t13 + t58 * t14 - t18 * t84 / 0.2e1 - t17 * t86 / 0.2e1 + t93 * t35 + t98 * t68 + t102 * t7 + t115 * t26 - t275; -0.2e1 * t273 - 0.2e1 * t274 + (-t152 * t125 + t209) * qJD(5) - t120 * t98 * t253 + (t102 * t93 + t20 * t58 + t21 * t57) * t251 + (t106 * t20 + t109 * t21 + t269 * t58 - t57 * t78) * t194 - t169 + 0.2e1 * t263 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) + (t252 * t98 + t247) * t115 + (t121 * t253 + t180 * t252 + 0.2e1 * mrSges(5,2) - 0.2e1 * t183) * t97 + t267; m(7) * (-t100 * t2 + t11 * t47 + t12 * t46 - t3 * t99) + t46 * t53 - t100 * t14 + t47 * t54 - t99 * t13 + t262 * t204 + (-t81 * mrSges(5,3) - t26 - t7 - m(7) * t24 - m(6) * t52 + m(5) * (-t52 + t215) + (-t107 * mrSges(5,3) - t152 * t71 + t156 * t70 + m(6) * (-t152 * t36 + t156 * t37)) * qJD(4)) * t157 + (-t80 * mrSges(5,3) - t152 * t33 + t156 * t34 + (-t152 * t70 - t156 * t71) * qJD(5) + m(6) * (-t10 * t152 + t9 * t156 - t197 * t36 - t198 * t37 - t216) + m(5) * (t51 - t216) + (t108 * mrSges(5,3) + t261 + t68) * qJD(4)) * t153; t276 * t157 + (t270 * t157 + (-t82 + t217) * t153) * qJD(4) + m(7) * (-t100 * t20 + t102 * t202 - t157 * t93 - t21 * t99 + t46 * t58 + t47 * t57) + m(6) * (-t218 + t153 * t182 + (t115 * t153 + t157 * t180) * qJD(4)) + m(5) * (t153 * t97 - t218 + (-t120 * t153 + t121 * t157) * qJD(4)) + t264; (-t100 * t46 - t47 * t99) * t251 + 0.4e1 * (m(6) * (-0.1e1 + t206) / 0.2e1 - m(7) / 0.2e1) * t153 * t201; m(7) * (t11 * t50 + t12 * t49 + t137 * t24 + t2 * t91 + t3 * t89) + t256 * t52 + (t113 * t236 - t10 * mrSges(6,3) + (-t55 / 0.2e1 - t37 * mrSges(6,3) + t126 * t236 + t261 * pkin(5)) * qJD(5) + (-m(6) * t181 - qJD(5) * t70 - t33) * pkin(9) + t177) * t152 + (t108 * t114 / 0.2e1 + t9 * mrSges(6,3) + (t56 / 0.2e1 - t36 * mrSges(6,3) + t108 * t235) * qJD(5) + (-m(6) * t186 - qJD(5) * t71 + t34) * pkin(9) + t178) * t156 - pkin(4) * t26 + t49 * t53 + t50 * t54 + t18 * t242 + t17 * t241 + t89 * t13 + t91 * t14 + t137 * t7 + t275; t273 - t272 + t274 - t263 + 0.2e1 * t242 * t269 + 0.2e1 * t241 * t78 + (pkin(4) + t115) * t111 - t270 * t97 + (-t209 + (-pkin(5) * t82 + t125) * t152) * qJD(5) + m(7) * (t102 * t191 + t137 * t93 + t20 * t91 + t21 * t89 + t49 * t58 + t50 * t57) + m(6) * (-pkin(4) * t98 + pkin(9) * t182) + ((-t58 + t91) * t269 + (t57 - t89) * t78 + (-t21 + t50) * t109 + (-t20 + t49) * t106) * mrSges(7,3) + t169 - t268; m(7) * (-t100 * t49 + t46 * t91 + t47 * t89 - t50 * t99) - t264 + (m(7) * t137 + t256 + t82) * t202 + (-m(7) * t191 + (m(6) * pkin(9) * t206 - t270) * qJD(4) - t276) * t157; t210 - t125 * t198 + (t137 * t191 + t49 * t91 + t50 * t89) * t251 + 0.2e1 * t82 * t191 + 0.2e1 * t272 + pkin(4) * t247 + (t113 + t199) * t156 + (-t106 * t49 - t109 * t50 - t269 * t91 + t78 * t89) * t194 + t267; -Ifges(6,5) * t189 + t10 * mrSges(6,1) - t9 * mrSges(6,2) - t167 * Ifges(6,6) + (m(7) * (-t11 * t196 + t12 * t195 + t151 * t2 + t155 * t3) + t53 * t195 + t151 * t14 - t54 * t196 + t155 * t13) * pkin(5) + t168 + t265; t138 - t174 * t97 + (t116 * t123 - t221) * qJD(5) + (m(7) * (t151 * t20 + t155 * t21 + (-t151 * t57 + t155 * t58) * qJD(6)) + t271) * pkin(5) + t172; (t153 * t198 - t156 * t201) * mrSges(6,2) + (-t152 * t201 - t153 * t197) * mrSges(6,1) + m(7) * (t151 * t46 + t155 * t47 + (-t100 * t155 + t151 * t99) * qJD(6)) * pkin(5) + t188; -t138 + (pkin(9) * t123 + t221) * qJD(5) + (m(7) * (t151 * t49 + t155 * t50 + (-t151 * t89 + t155 * t91) * qJD(6)) - t271) * pkin(5) + t171; 0.2e1 * t105; t168; t172; t188; t171; t105; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
