% Calculate time derivative of joint inertia matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:48
% EndTime: 2019-03-09 05:04:54
% DurationCPUTime: 3.27s
% Computational Cost: add. (2874->467), mult. (6504->651), div. (0->0), fcn. (5052->8), ass. (0->180)
t195 = Ifges(6,4) + Ifges(5,5);
t210 = -Ifges(6,2) - Ifges(5,3);
t128 = sin(qJ(4));
t131 = cos(qJ(4));
t171 = t128 ^ 2 + t131 ^ 2;
t129 = sin(qJ(3));
t132 = cos(qJ(3));
t169 = qJD(3) * t132;
t156 = t128 * t169;
t166 = qJD(4) * t131;
t136 = t129 * t166 + t156;
t161 = -cos(pkin(10)) * pkin(1) - pkin(2);
t209 = 0.2e1 * t161;
t127 = sin(qJ(6));
t130 = cos(qJ(6));
t141 = t127 * t131 - t128 * t130;
t69 = t141 * t129;
t208 = -t127 * mrSges(7,1) - t130 * mrSges(7,2);
t123 = t132 * pkin(4);
t76 = -pkin(3) * t132 - pkin(8) * t129 + t161;
t116 = sin(pkin(10)) * pkin(1) + pkin(7);
t177 = t116 * t132;
t91 = t128 * t177;
t27 = pkin(5) * t132 + t123 + t91 + (-pkin(9) * t129 - t76) * t131;
t176 = t128 * t129;
t174 = t131 * t132;
t92 = t116 * t174;
t44 = t128 * t76 + t92;
t36 = -qJ(5) * t132 + t44;
t30 = pkin(9) * t176 + t36;
t6 = -t127 * t30 + t130 * t27;
t168 = qJD(4) * t128;
t95 = (pkin(3) * t129 - pkin(8) * t132) * qJD(3);
t150 = qJD(4) * t92 - t131 * t95 + t76 * t168;
t152 = -t116 * t128 - pkin(4);
t167 = qJD(4) * t129;
t158 = t128 * t167;
t8 = pkin(9) * t158 + (-pkin(9) * t174 + (-pkin(5) + t152) * t129) * qJD(3) + t150;
t170 = qJD(3) * t129;
t117 = qJ(5) * t170;
t175 = t129 * t131;
t191 = t128 * t95 + t76 * t166;
t9 = t117 + (pkin(9) * qJD(4) - qJD(3) * t116) * t175 + (-qJD(5) + (pkin(9) * qJD(3) - qJD(4) * t116) * t128) * t132 + t191;
t1 = qJD(6) * t6 + t127 * t8 + t130 * t9;
t7 = t127 * t27 + t130 * t30;
t2 = -qJD(6) * t7 - t127 * t9 + t130 * t8;
t207 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t206 = qJD(4) - qJD(6);
t179 = qJ(5) * t128;
t201 = pkin(4) + pkin(5);
t205 = -t131 * t201 - t179;
t204 = 2 * m(5);
t203 = 2 * m(7);
t202 = 0.2e1 * t116;
t200 = pkin(8) - pkin(9);
t96 = -t127 * qJ(5) - t130 * t201;
t67 = t130 * qJD(5) + qJD(6) * t96;
t197 = t67 * mrSges(7,2);
t97 = t130 * qJ(5) - t127 * t201;
t68 = -t127 * qJD(5) - qJD(6) * t97;
t196 = t68 * mrSges(7,1);
t140 = t127 * t128 + t130 * t131;
t45 = t206 * t140;
t22 = t129 * t45 - t141 * t169;
t23 = t140 * t169 + t206 * t69;
t194 = Ifges(7,5) * t23 + Ifges(7,6) * t22;
t155 = t131 * t169;
t137 = t155 - t158;
t56 = mrSges(5,1) * t170 - mrSges(5,3) * t137;
t57 = mrSges(6,2) * t155 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t168) * t129;
t193 = -t56 + t57;
t58 = -mrSges(5,2) * t170 - mrSges(5,3) * t136;
t59 = -mrSges(6,2) * t136 + mrSges(6,3) * t170;
t192 = t58 + t59;
t185 = Ifges(6,5) * t128;
t146 = Ifges(6,1) * t131 + t185;
t65 = -Ifges(6,4) * t132 + t129 * t146;
t187 = Ifges(5,4) * t128;
t147 = Ifges(5,1) * t131 - t187;
t66 = -Ifges(5,5) * t132 + t129 * t147;
t190 = t65 + t66;
t87 = mrSges(5,2) * t132 - mrSges(5,3) * t176;
t90 = -mrSges(6,2) * t176 - mrSges(6,3) * t132;
t189 = t87 + t90;
t88 = -mrSges(5,1) * t132 - mrSges(5,3) * t175;
t89 = mrSges(6,1) * t132 + mrSges(6,2) * t175;
t188 = -t88 + t89;
t186 = Ifges(5,4) * t131;
t184 = Ifges(6,5) * t131;
t181 = t132 * Ifges(5,6);
t100 = -t131 * mrSges(5,1) + t128 * mrSges(5,2);
t180 = t100 - mrSges(4,1);
t178 = qJ(5) * t131;
t164 = qJD(5) * t131;
t173 = qJ(5) * t155 + t129 * t164;
t172 = t171 * pkin(8) * t169;
t165 = qJD(5) * t128;
t162 = t201 * t128;
t106 = t200 * t131;
t160 = t116 * t170;
t159 = t129 * t169;
t101 = -Ifges(6,3) * t131 + t185;
t102 = Ifges(5,2) * t131 + t187;
t154 = -t102 / 0.2e1 + t101 / 0.2e1;
t103 = Ifges(6,1) * t128 - t184;
t104 = Ifges(5,1) * t128 + t186;
t153 = t103 / 0.2e1 + t104 / 0.2e1;
t43 = t131 * t76 - t91;
t144 = Ifges(6,3) * t128 + t184;
t63 = -Ifges(6,6) * t132 + t129 * t144;
t145 = -Ifges(5,2) * t128 + t186;
t64 = t129 * t145 - t181;
t151 = t63 - t64 + t181;
t5 = -t22 * mrSges(7,1) + t23 * mrSges(7,2);
t149 = mrSges(5,1) * t128 + mrSges(5,2) * t131;
t99 = -t131 * mrSges(6,1) - t128 * mrSges(6,3);
t148 = mrSges(6,1) * t128 - mrSges(6,3) * t131;
t143 = -pkin(4) * t131 - t179;
t142 = pkin(4) * t128 - t178;
t105 = t200 * t128;
t51 = t105 * t130 - t106 * t127;
t52 = t105 * t127 + t106 * t130;
t139 = -t162 + t178;
t93 = t200 * t168;
t94 = qJD(4) * t106;
t24 = qJD(6) * t51 + t127 * t94 - t130 * t93;
t25 = -qJD(6) * t52 + t127 * t93 + t130 * t94;
t46 = t206 * t141;
t40 = Ifges(7,6) * t46;
t41 = Ifges(7,5) * t45;
t138 = t25 * mrSges(7,1) - t24 * mrSges(7,2) - t40 + t41;
t135 = -t136 * Ifges(6,6) - t195 * t155 + t210 * t170 + t194;
t20 = (-t131 * t170 - t132 * t168) * t116 + t191;
t134 = (m(6) * t143 + t100 + t99) * qJD(4);
t122 = Ifges(6,4) * t166;
t121 = Ifges(5,5) * t166;
t119 = Ifges(6,6) * t168;
t98 = -pkin(3) + t143;
t86 = t147 * qJD(4);
t85 = t146 * qJD(4);
t84 = t145 * qJD(4);
t83 = t144 * qJD(4);
t82 = t149 * qJD(4);
t81 = t148 * qJD(4);
t77 = pkin(3) - t205;
t75 = t149 * t129;
t74 = t148 * t129;
t72 = qJD(4) * t142 - t165;
t70 = t140 * t129;
t60 = qJD(4) * t139 + t165;
t55 = (t116 + t142) * t129;
t54 = mrSges(7,1) * t132 - mrSges(7,3) * t70;
t53 = -mrSges(7,2) * t132 - mrSges(7,3) * t69;
t49 = -Ifges(7,1) * t141 - Ifges(7,4) * t140;
t48 = -Ifges(7,4) * t141 - Ifges(7,2) * t140;
t47 = mrSges(7,1) * t140 - mrSges(7,2) * t141;
t42 = (-t116 + t139) * t129;
t39 = mrSges(5,1) * t136 + mrSges(5,2) * t137;
t38 = mrSges(6,1) * t136 - mrSges(6,3) * t137;
t37 = t123 - t43;
t35 = mrSges(7,1) * t69 + mrSges(7,2) * t70;
t34 = -t104 * t167 + (Ifges(5,5) * t129 + t132 * t147) * qJD(3);
t33 = -t103 * t167 + (Ifges(6,4) * t129 + t132 * t146) * qJD(3);
t32 = -t102 * t167 + (Ifges(5,6) * t129 + t132 * t145) * qJD(3);
t31 = -t101 * t167 + (Ifges(6,6) * t129 + t132 * t144) * qJD(3);
t29 = Ifges(7,1) * t70 - Ifges(7,4) * t69 + Ifges(7,5) * t132;
t28 = Ifges(7,4) * t70 - Ifges(7,2) * t69 + Ifges(7,6) * t132;
t26 = pkin(4) * t136 + qJ(5) * t158 + t116 * t169 - t173;
t21 = t128 * t160 - t150;
t17 = Ifges(7,1) * t45 - Ifges(7,4) * t46;
t16 = Ifges(7,4) * t45 - Ifges(7,2) * t46;
t15 = mrSges(7,1) * t46 + mrSges(7,2) * t45;
t14 = t152 * t170 + t150;
t13 = -mrSges(7,1) * t170 - mrSges(7,3) * t23;
t12 = mrSges(7,2) * t170 + mrSges(7,3) * t22;
t11 = t205 * t167 + (-t116 - t162) * t169 + t173;
t10 = -qJD(5) * t132 + t117 + t20;
t4 = Ifges(7,1) * t23 + Ifges(7,4) * t22 - Ifges(7,5) * t170;
t3 = Ifges(7,4) * t23 + Ifges(7,2) * t22 - Ifges(7,6) * t170;
t18 = [((mrSges(4,2) * t209 + 0.2e1 * Ifges(4,4) * t132 + t128 * t151 + t131 * t190 + t202 * t75) * qJD(3) + t135) * t132 + (t39 * t202 + (t33 + t34) * t131 + (t31 - t32) * t128 + (t151 * t131 + (t132 * t195 - t190) * t128) * qJD(4) + (-Ifges(7,5) * t70 + Ifges(7,6) * t69 + mrSges(4,1) * t209 + (-0.2e1 * Ifges(4,4) + t195 * t131 + (-Ifges(5,6) + Ifges(6,6)) * t128) * t129 + (t116 ^ 2 * t204 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(7,3)) + t210) * t132) * qJD(3)) * t129 + 0.2e1 * t7 * t12 + 0.2e1 * t6 * t13 + t22 * t28 + t23 * t29 + 0.2e1 * t11 * t35 + 0.2e1 * t42 * t5 + 0.2e1 * t1 * t53 + 0.2e1 * t2 * t54 + 0.2e1 * t55 * t38 + 0.2e1 * t43 * t56 + 0.2e1 * t37 * t57 + 0.2e1 * t44 * t58 + 0.2e1 * t36 * t59 - t69 * t3 + t70 * t4 + 0.2e1 * t26 * t74 + 0.2e1 * t20 * t87 + 0.2e1 * t21 * t88 + 0.2e1 * t14 * t89 + 0.2e1 * t10 * t90 + (t1 * t7 + t11 * t42 + t2 * t6) * t203 + (t20 * t44 + t21 * t43) * t204 + 0.2e1 * m(6) * (t10 * t36 + t14 * t37 + t26 * t55); m(7) * (t1 * t70 - t2 * t69 + t22 * t6 + t23 * t7) + t23 * t53 + t22 * t54 - t69 * t13 + t70 * t12 + (t5 - t38 - t39 - m(6) * t26 + m(7) * t11 + (t189 * t131 + t188 * t128 + m(6) * (t128 * t37 + t131 * t36) + (-t128 * t43 + t131 * t44 - t177) * m(5)) * qJD(3)) * t132 + (t192 * t131 + t193 * t128 + (-t128 * t189 + t131 * t188) * qJD(4) + m(6) * (t10 * t131 + t128 * t14 + t166 * t37 - t168 * t36) + (m(6) * t55 - m(7) * t42 - t35 + t74 + t75) * qJD(3) + (-t128 * t21 + t131 * t20 - t166 * t43 - t168 * t44 + t160) * m(5)) * t129; (-t22 * t69 + t23 * t70 - t159) * t203 + 0.4e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t171) * t159; (-t121 / 0.2e1 + t41 / 0.2e1 - t40 / 0.2e1 - t122 / 0.2e1 - t119 / 0.2e1 + (Ifges(4,5) + (-m(5) * pkin(3) + t180) * t116) * qJD(3)) * t132 + m(7) * (t1 * t52 + t11 * t77 + t2 * t51 + t24 * t7 + t25 * t6 + t42 * t60) + (t10 * mrSges(6,2) + t20 * mrSges(5,3) - t31 / 0.2e1 + t32 / 0.2e1 + t153 * t169 + (t37 * mrSges(6,2) - t43 * mrSges(5,3) + t65 / 0.2e1 + t66 / 0.2e1) * qJD(4) + (t188 * qJD(4) + m(6) * (t37 * qJD(4) + t10) + m(5) * (-t43 * qJD(4) + t20) + t192) * pkin(8)) * t131 + (-t21 * mrSges(5,3) + t14 * mrSges(6,2) + t33 / 0.2e1 + t34 / 0.2e1 + t154 * t169 + (-t36 * mrSges(6,2) - t44 * mrSges(5,3) + t63 / 0.2e1 - t64 / 0.2e1 + t181 / 0.2e1) * qJD(4) + (-t189 * qJD(4) + m(6) * (-t36 * qJD(4) + t14) + m(5) * (-t44 * qJD(4) - t21) + t193) * pkin(8)) * t128 - t140 * t3 / 0.2e1 + (t116 * t82 + (t85 / 0.2e1 + t86 / 0.2e1) * t131 + (t83 / 0.2e1 - t84 / 0.2e1) * t128 + (t116 * mrSges(4,2) + Ifges(7,5) * t141 / 0.2e1 + Ifges(7,6) * t140 / 0.2e1 - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t131 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t128) * qJD(3) + (-t128 * t153 + t131 * t154) * qJD(4)) * t129 - t141 * t4 / 0.2e1 + (-t1 * t140 + t141 * t2 - t6 * t45 - t7 * t46) * mrSges(7,3) + m(6) * (t26 * t98 + t55 * t72) - pkin(3) * t39 + t42 * t15 + t45 * t29 / 0.2e1 - t46 * t28 / 0.2e1 + t11 * t47 + t22 * t48 / 0.2e1 + t23 * t49 / 0.2e1 + t51 * t13 + t52 * t12 + t24 * t53 + t25 * t54 + t60 * t35 - t69 * t16 / 0.2e1 + t70 * t17 / 0.2e1 + t72 * t74 + t77 * t5 + t55 * t81 + t98 * t38 + t26 * t99; (t15 - t81 - t82) * t132 + ((-t47 + t99 + t180) * t129 + (-mrSges(4,2) + (mrSges(6,2) + mrSges(5,3)) * t171) * t132) * qJD(3) + m(6) * (-t132 * t72 + t170 * t98 + t172) + m(7) * (t132 * t60 - t170 * t77 + t22 * t51 + t23 * t52 + t24 * t70 - t25 * t69) + m(5) * (-pkin(3) * t170 + t172) + (-t140 * t23 + t141 * t22 + t45 * t69 - t46 * t70) * mrSges(7,3); 0.2e1 * t98 * t81 - 0.2e1 * pkin(3) * t82 - t46 * t48 - t140 * t16 + t45 * t49 - t141 * t17 + 0.2e1 * t60 * t47 + 0.2e1 * t77 * t15 + (t24 * t52 + t25 * t51 + t60 * t77) * t203 + (-t83 + t84) * t131 + (t85 + t86) * t128 + ((t103 + t104) * t131 + (t101 - t102) * t128) * qJD(4) + 0.2e1 * (m(6) * t98 + t99) * t72 + 0.2e1 * (-t140 * t24 + t141 * t25 - t45 * t51 - t46 * t52) * mrSges(7,3); -Ifges(5,6) * t156 + (Ifges(7,3) * qJD(3) + (-Ifges(5,6) * t131 - t128 * t195) * qJD(4)) * t129 + m(6) * (-pkin(4) * t14 + qJ(5) * t10 + qJD(5) * t36) + m(7) * (t1 * t97 + t2 * t96 + t6 * t68 + t67 * t7) - t135 + t10 * mrSges(6,3) - t14 * mrSges(6,1) - t20 * mrSges(5,2) + t21 * mrSges(5,1) - pkin(4) * t57 + qJ(5) * t59 + t67 * t53 + t68 * t54 + qJD(5) * t90 + t96 * t13 + t97 * t12 + t207; ((-mrSges(5,2) + mrSges(6,3)) * t131 + (-mrSges(5,1) - mrSges(6,1)) * t128) * t169 + m(6) * (-pkin(4) * t156 + t173) + m(7) * (t22 * t96 + t23 * t97 + t67 * t70 - t68 * t69) + t129 * t134 + t5; m(7) * (t24 * t97 + t25 * t96 + t51 * t68 + t52 * t67) + t122 + t119 + t121 - Ifges(5,6) * t168 + (qJD(4) * t143 + t164) * mrSges(6,2) + (-t140 * t67 + t141 * t68 - t45 * t96 - t46 * t97) * mrSges(7,3) + (m(6) * t164 + t134) * pkin(8) - t138; -0.2e1 * t196 + (t67 * t97 + t68 * t96) * t203 + 0.2e1 * t197 + 0.2e1 * (m(6) * qJ(5) + mrSges(6,3)) * qJD(5); t127 * t12 + t130 * t13 + (-t127 * t54 + t130 * t53) * qJD(6) + m(7) * (t1 * t127 + t130 * t2 + (-t127 * t6 + t130 * t7) * qJD(6)) + m(6) * t14 + t57; m(7) * (t127 * t23 + t130 * t22 + (t127 * t69 + t130 * t70) * qJD(6)) + t136 * m(6); m(7) * (t127 * t24 + t130 * t25 + (-t127 * t51 + t130 * t52) * qJD(6)) + (m(6) * pkin(8) + mrSges(6,2)) * t166 + (-t127 * t46 - t130 * t45 + (-t127 * t141 - t130 * t140) * qJD(6)) * mrSges(7,3); m(7) * (t127 * t67 + t130 * t68) + (m(7) * (-t127 * t96 + t130 * t97) - t208) * qJD(6); 0; -Ifges(7,3) * t170 + t194 - t207; -t5; t138; t196 - t197; t208 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t18(1) t18(2) t18(4) t18(7) t18(11) t18(16); t18(2) t18(3) t18(5) t18(8) t18(12) t18(17); t18(4) t18(5) t18(6) t18(9) t18(13) t18(18); t18(7) t18(8) t18(9) t18(10) t18(14) t18(19); t18(11) t18(12) t18(13) t18(14) t18(15) t18(20); t18(16) t18(17) t18(18) t18(19) t18(20) t18(21);];
Mq  = res;
