% Calculate time derivative of joint inertia matrix for
% S6RRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:52:46
% EndTime: 2019-03-09 11:52:58
% DurationCPUTime: 5.89s
% Computational Cost: add. (6624->440), mult. (14589->626), div. (0->0), fcn. (13840->8), ass. (0->179)
t238 = Ifges(7,4) + Ifges(6,5);
t236 = -Ifges(6,6) + Ifges(7,6);
t152 = sin(pkin(10));
t153 = cos(pkin(10));
t156 = sin(qJ(2));
t159 = cos(qJ(2));
t130 = t152 * t159 + t153 * t156;
t126 = t130 * qJD(2);
t154 = sin(qJ(5));
t155 = sin(qJ(4));
t157 = cos(qJ(5));
t158 = cos(qJ(4));
t133 = t154 * t158 + t155 * t157;
t232 = qJD(4) + qJD(5);
t102 = t232 * t133;
t129 = t152 * t156 - t153 * t159;
t127 = t129 * qJD(2);
t132 = t154 * t155 - t157 * t158;
t33 = -t102 * t130 + t132 * t127;
t195 = qJD(5) * t154;
t201 = t155 * t127;
t202 = t130 * t158;
t203 = t130 * t155;
t198 = qJD(4) * t155;
t187 = t130 * t198;
t199 = t158 * t127;
t233 = -t199 - t187;
t34 = -t195 * t203 + (t202 * t232 - t201) * t157 + t233 * t154;
t241 = (-Ifges(6,4) + Ifges(7,5)) * t34 + (Ifges(6,1) + Ifges(7,1)) * t33 + t238 * t126;
t101 = t232 * t132;
t240 = -t238 * t101 + t236 * t102;
t239 = -mrSges(6,1) - mrSges(7,1);
t237 = Ifges(7,2) + Ifges(6,3);
t197 = qJD(4) * t158;
t169 = t130 * t197 - t201;
t235 = m(7) + m(6);
t218 = -qJ(3) - pkin(7);
t138 = t218 * t156;
t139 = t218 * t159;
t110 = t138 * t152 - t139 * t153;
t100 = t158 * t110;
t149 = -pkin(2) * t159 - pkin(1);
t91 = t129 * pkin(3) - t130 * pkin(8) + t149;
t59 = t155 * t91 + t100;
t234 = -Ifges(5,5) * t199 + Ifges(5,3) * t126;
t175 = mrSges(5,1) * t155 + mrSges(5,2) * t158;
t134 = t175 * qJD(4);
t180 = qJD(2) * t218;
t125 = qJD(3) * t159 + t156 * t180;
t164 = -t156 * qJD(3) + t159 * t180;
t81 = t153 * t125 + t152 * t164;
t192 = pkin(2) * qJD(2) * t156;
t82 = pkin(3) * t126 + pkin(8) * t127 + t192;
t181 = -t155 * t81 + t158 * t82;
t15 = pkin(9) * t199 + pkin(4) * t126 + (-t100 + (pkin(9) * t130 - t91) * t155) * qJD(4) + t181;
t22 = -t110 * t198 + t155 * t82 + t158 * t81 + t91 * t197;
t19 = -pkin(9) * t169 + t22;
t58 = -t110 * t155 + t158 * t91;
t39 = pkin(4) * t129 - pkin(9) * t202 + t58;
t48 = -pkin(9) * t203 + t59;
t215 = t154 * t39 + t157 * t48;
t6 = -qJD(5) * t215 + t15 * t157 - t154 * t19;
t231 = 2 * m(6);
t230 = 2 * m(7);
t229 = 0.2e1 * pkin(4);
t228 = -2 * mrSges(4,3);
t109 = -t153 * t138 - t139 * t152;
t226 = 0.2e1 * t109;
t225 = 0.2e1 * t149;
t224 = t126 / 0.2e1;
t222 = -t130 / 0.2e1;
t210 = Ifges(5,4) * t155;
t140 = Ifges(5,2) * t158 + t210;
t220 = -t140 / 0.2e1;
t144 = pkin(2) * t152 + pkin(8);
t219 = pkin(9) + t144;
t24 = mrSges(6,1) * t126 - mrSges(6,3) * t33;
t25 = -t126 * mrSges(7,1) + t33 * mrSges(7,2);
t217 = -t24 + t25;
t26 = -mrSges(6,2) * t126 - mrSges(6,3) * t34;
t27 = -mrSges(7,2) * t34 + mrSges(7,3) * t126;
t216 = t26 + t27;
t84 = t133 * t130;
t69 = -mrSges(6,2) * t129 - mrSges(6,3) * t84;
t72 = -mrSges(7,2) * t84 + mrSges(7,3) * t129;
t214 = t69 + t72;
t85 = t132 * t130;
t70 = mrSges(6,1) * t129 + mrSges(6,3) * t85;
t71 = -mrSges(7,1) * t129 - mrSges(7,2) * t85;
t213 = -t70 + t71;
t209 = Ifges(5,4) * t158;
t208 = Ifges(5,6) * t155;
t80 = t125 * t152 - t153 * t164;
t207 = t109 * t80;
t206 = t126 * Ifges(5,5);
t205 = t126 * Ifges(5,6);
t204 = t129 * Ifges(5,6);
t194 = qJD(5) * t157;
t191 = pkin(4) * t198;
t190 = pkin(4) * t195;
t189 = pkin(4) * t194;
t128 = t219 * t158;
t183 = t219 * t155;
t89 = t154 * t128 + t157 * t183;
t188 = t89 * t195;
t146 = -pkin(2) * t153 - pkin(3);
t185 = t132 * t195;
t179 = qJD(4) * t219;
t124 = t155 * t179;
t176 = t158 * t179;
t56 = -qJD(5) * t89 - t157 * t124 - t154 * t176;
t90 = t157 * t128 - t154 * t183;
t57 = qJD(5) * t90 - t154 * t124 + t157 * t176;
t184 = t90 * t56 + t57 * t89;
t61 = t102 * mrSges(7,1) + t101 * mrSges(7,3);
t62 = t102 * mrSges(6,1) - t101 * mrSges(6,2);
t182 = -(2 * Ifges(4,4)) - t208;
t178 = 0.2e1 * t192;
t78 = pkin(4) * t203 + t109;
t174 = Ifges(5,1) * t158 - t210;
t173 = -Ifges(5,2) * t155 + t209;
t20 = -t154 * t48 + t157 * t39;
t137 = -pkin(4) * t158 + t146;
t170 = t126 * t237 + t236 * t34 + t238 * t33;
t54 = pkin(4) * t169 + t80;
t5 = t154 * t15 + t157 * t19 + t39 * t194 - t195 * t48;
t167 = -t61 - t62;
t165 = -pkin(5) * t102 - qJ(6) * t101 + qJD(6) * t133;
t163 = t239 * t57 + (-mrSges(6,2) + mrSges(7,3)) * t56 + t240;
t2 = qJ(6) * t126 + qJD(6) * t129 + t5;
t3 = -pkin(5) * t126 - t6;
t162 = t6 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t2 * mrSges(7,3) + t170;
t142 = qJD(6) + t189;
t160 = -mrSges(6,2) * t189 + t142 * mrSges(7,3) + t190 * t239;
t151 = qJD(6) * mrSges(7,3);
t150 = Ifges(5,5) * t197;
t148 = -pkin(4) * t157 - pkin(5);
t145 = pkin(4) * t154 + qJ(6);
t141 = Ifges(5,1) * t155 + t209;
t136 = t174 * qJD(4);
t135 = t173 * qJD(4);
t118 = t127 * mrSges(4,2);
t108 = Ifges(6,1) * t133 - Ifges(6,4) * t132;
t107 = Ifges(7,1) * t133 + Ifges(7,5) * t132;
t106 = Ifges(6,4) * t133 - Ifges(6,2) * t132;
t105 = Ifges(7,5) * t133 + Ifges(7,3) * t132;
t104 = mrSges(6,1) * t132 + mrSges(6,2) * t133;
t103 = mrSges(7,1) * t132 - mrSges(7,3) * t133;
t93 = mrSges(5,1) * t129 - mrSges(5,3) * t202;
t92 = -mrSges(5,2) * t129 - mrSges(5,3) * t203;
t86 = pkin(5) * t132 - qJ(6) * t133 + t137;
t74 = Ifges(5,5) * t129 + t130 * t174;
t73 = t130 * t173 + t204;
t68 = -mrSges(5,2) * t126 - mrSges(5,3) * t169;
t67 = mrSges(5,1) * t126 - mrSges(5,3) * t233;
t66 = -Ifges(6,1) * t101 - Ifges(6,4) * t102;
t65 = -Ifges(7,1) * t101 + Ifges(7,5) * t102;
t64 = -Ifges(6,4) * t101 - Ifges(6,2) * t102;
t63 = -Ifges(7,5) * t101 + Ifges(7,3) * t102;
t55 = mrSges(5,1) * t169 + mrSges(5,2) * t233;
t51 = mrSges(6,1) * t84 - mrSges(6,2) * t85;
t50 = mrSges(7,1) * t84 + mrSges(7,3) * t85;
t49 = -t165 + t191;
t45 = -Ifges(6,1) * t85 - Ifges(6,4) * t84 + Ifges(6,5) * t129;
t44 = -Ifges(7,1) * t85 + Ifges(7,4) * t129 + Ifges(7,5) * t84;
t43 = -Ifges(6,4) * t85 - Ifges(6,2) * t84 + Ifges(6,6) * t129;
t42 = -Ifges(7,5) * t85 + Ifges(7,6) * t129 + Ifges(7,3) * t84;
t41 = Ifges(5,1) * t233 - Ifges(5,4) * t169 + t206;
t40 = Ifges(5,4) * t233 - Ifges(5,2) * t169 + t205;
t35 = pkin(5) * t84 + qJ(6) * t85 + t78;
t23 = -qJD(4) * t59 + t181;
t17 = -pkin(5) * t129 - t20;
t16 = qJ(6) * t129 + t215;
t13 = mrSges(6,1) * t34 + mrSges(6,2) * t33;
t12 = mrSges(7,1) * t34 - mrSges(7,3) * t33;
t9 = Ifges(6,4) * t33 - Ifges(6,2) * t34 + t126 * Ifges(6,6);
t8 = Ifges(7,5) * t33 + t126 * Ifges(7,6) + Ifges(7,3) * t34;
t7 = pkin(5) * t34 - qJ(6) * t33 + qJD(6) * t85 + t54;
t1 = [0.2e1 * (-pkin(1) * (mrSges(3,1) * t156 + mrSges(3,2) * t159) + (-Ifges(3,2) + Ifges(3,1)) * t156 * t159 + (-t156 ^ 2 + t159 ^ 2) * Ifges(3,4)) * qJD(2) + (mrSges(4,2) * t178 + t158 * t41 - t155 * t40 - 0.2e1 * Ifges(4,1) * t127 + (Ifges(5,5) * t158 + t182) * t126 + (-t158 * t73 - t155 * t74 + t129 * (-Ifges(5,5) * t155 - Ifges(5,6) * t158)) * qJD(4) + 0.2e1 * (mrSges(4,3) + t175) * t80) * t130 + (t20 * t6 + t215 * t5 + t54 * t78) * t231 + 0.2e1 * t215 * t26 - t241 * t85 + (t8 - t9) * t84 - (mrSges(4,3) * t226 - t155 * t73 + t158 * t74) * t127 + 0.2e1 * m(4) * (t110 * t81 + t149 * t192 + t207) + 0.2e1 * m(5) * (t22 * t59 + t23 * t58 + t207) - t118 * t225 + t55 * t226 + (t16 * t2 + t17 * t3 + t35 * t7) * t230 + (mrSges(4,1) * t178 + t81 * t228 - t182 * t127 + ((2 * Ifges(4,2)) + Ifges(5,3) + t237) * t126 + t170 + t234) * t129 + (mrSges(4,1) * t225 + t110 * t228 + t236 * t84 - t238 * t85) * t126 + (t42 - t43) * t34 + (t44 + t45) * t33 + 0.2e1 * t20 * t24 + 0.2e1 * t17 * t25 + 0.2e1 * t16 * t27 + 0.2e1 * t35 * t12 + 0.2e1 * t7 * t50 + 0.2e1 * t54 * t51 + 0.2e1 * t58 * t67 + 0.2e1 * t59 * t68 + 0.2e1 * t5 * t69 + 0.2e1 * t6 * t70 + 0.2e1 * t3 * t71 + 0.2e1 * t2 * t72 + 0.2e1 * t78 * t13 + 0.2e1 * t22 * t92 + 0.2e1 * t23 * t93; (t130 * t136 / 0.2e1 - t127 * t141 / 0.2e1 + t22 * mrSges(5,3) + t40 / 0.2e1 - t80 * mrSges(5,1) + t205 / 0.2e1 + (t130 * t220 - t58 * mrSges(5,3) + t74 / 0.2e1) * qJD(4) + (-qJD(4) * t93 + m(5) * (-qJD(4) * t58 + t22) + t68) * t144) * t158 + m(6) * (t137 * t54 - t20 * t57 + t215 * t56 + t5 * t90 - t6 * t89) + (t241 / 0.2e1 + t3 * mrSges(7,2) - t6 * mrSges(6,3) + t238 * t224) * t133 + (m(5) * t80 + t55) * t146 + (-t2 * mrSges(7,2) - t5 * mrSges(6,3) + t8 / 0.2e1 - t9 / 0.2e1 + t236 * t224) * t132 + (t150 + t240) * t129 / 0.2e1 + ((-t126 * t152 + t127 * t153) * mrSges(4,3) + m(4) * (t152 * t81 - t153 * t80)) * pkin(2) + (Ifges(3,5) * t159 - Ifges(3,6) * t156 + (-mrSges(3,1) * t159 + mrSges(3,2) * t156) * pkin(7)) * qJD(2) + t217 * t89 + t213 * t57 + t214 * t56 + t216 * t90 + m(7) * (t16 * t56 + t17 * t57 + t2 * t90 + t3 * t89 + t35 * t49 + t7 * t86) - (t65 / 0.2e1 + t66 / 0.2e1) * t85 + (t63 / 0.2e1 - t64 / 0.2e1) * t84 + (t105 / 0.2e1 - t106 / 0.2e1) * t34 + (t107 / 0.2e1 + t108 / 0.2e1) * t33 + (t42 / 0.2e1 - t43 / 0.2e1) * t102 - (t44 / 0.2e1 + t45 / 0.2e1) * t101 + (t20 * t101 - t102 * t215) * mrSges(6,3) + (-t101 * t17 - t102 * t16) * mrSges(7,2) + (t135 * t222 - t127 * t220 - t23 * mrSges(5,3) + t41 / 0.2e1 + t80 * mrSges(5,2) + t206 / 0.2e1 + (t141 * t222 - t59 * mrSges(5,3) - t73 / 0.2e1 - t204 / 0.2e1 + (m(6) * t78 + t51) * pkin(4)) * qJD(4) + (-m(5) * t23 - t67 + (-m(5) * t59 - t92) * qJD(4)) * t144) * t155 + t49 * t50 + t35 * t61 + t78 * t62 - t80 * mrSges(4,1) - t81 * mrSges(4,2) + t86 * t12 + t7 * t103 + t54 * t104 - Ifges(4,6) * t126 - Ifges(4,5) * t127 + t109 * t134 + t137 * t13; 0.2e1 * t49 * t103 + 0.2e1 * t146 * t134 + t158 * t135 + t155 * t136 + 0.2e1 * t137 * t62 + 0.2e1 * t86 * t61 + (t65 + t66) * t133 + (t63 - t64) * t132 + (t105 - t106) * t102 - (t107 + t108) * t101 + (t158 * t141 + (t104 * t229 - t140) * t155) * qJD(4) + (t137 * t191 + t184) * t231 + (t49 * t86 + t184) * t230 + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * (-t101 * t89 - t102 * t90 - t132 * t56 + t133 * t57); m(4) * t192 + t126 * mrSges(4,1) + t155 * t68 + t158 * t67 - t118 + t216 * t133 + t217 * t132 + t213 * t102 - t214 * t101 + (-t155 * t93 + t158 * t92) * qJD(4) + m(7) * (-t101 * t16 + t102 * t17 + t132 * t3 + t133 * t2) + m(6) * (-t101 * t215 - t102 * t20 - t132 * t6 + t133 * t5) + m(5) * (t155 * t22 + t158 * t23 + (-t155 * t58 + t158 * t59) * qJD(4)); t235 * (-t101 * t90 + t102 * t89 + t132 * t57 + t133 * t56); 0.2e1 * t235 * (-t133 * t101 + t102 * t132); t162 + m(7) * (t142 * t16 + t145 * t2 + t148 * t3) + ((m(6) * t6 + t24 + (m(6) * t215 + t69) * qJD(5)) * t157 + (m(6) * t5 + t26 + (-m(6) * t20 + m(7) * t17 + t213) * qJD(5)) * t154) * pkin(4) - Ifges(5,5) * t187 - t169 * Ifges(5,6) - t22 * mrSges(5,2) + t23 * mrSges(5,1) + t142 * t72 + t145 * t27 + t148 * t25 + t234; m(7) * (t142 * t90 + t145 * t56 + t148 * t57) + t150 + (-t208 + (-mrSges(5,1) * t158 + mrSges(5,2) * t155) * t144) * qJD(4) + (-t101 * t148 - t102 * t145 - t132 * t142) * mrSges(7,2) + (t133 * mrSges(7,2) * t195 + m(7) * t188 + m(6) * (t154 * t56 - t157 * t57 + t194 * t90 + t188) + (t157 * t101 - t154 * t102 + (-t132 * t157 + t133 * t154) * qJD(5)) * mrSges(6,3)) * pkin(4) + t163; -t134 + m(7) * (-t101 * t145 + t102 * t148 + t133 * t142) + (m(6) * (-t101 * t154 - t102 * t157 + t133 * t194 + t185) / 0.2e1 + m(7) * t185 / 0.2e1) * t229 + t167; 0.2e1 * m(7) * (t142 * t145 + t148 * t190) + 0.2e1 * t160; t162 + m(7) * (-pkin(5) * t3 + qJ(6) * t2 + qJD(6) * t16) - pkin(5) * t25 + qJ(6) * t27 + qJD(6) * t72; m(7) * (-pkin(5) * t57 + qJ(6) * t56 + qJD(6) * t90) + (pkin(5) * t101 - qJ(6) * t102 - qJD(6) * t132) * mrSges(7,2) + t163; m(7) * t165 + t167; m(7) * (-pkin(5) * t190 + qJ(6) * t142 + qJD(6) * t145) + t151 + t160; 0.2e1 * m(7) * qJ(6) * qJD(6) + 0.2e1 * t151; m(7) * t3 + t25; m(7) * t57 - t101 * mrSges(7,2); m(7) * t102; m(7) * t190; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
