% Calculate time derivative of joint inertia matrix for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR16_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:45
% EndTime: 2019-12-31 20:44:52
% DurationCPUTime: 2.85s
% Computational Cost: add. (2693->403), mult. (7021->603), div. (0->0), fcn. (6044->8), ass. (0->178)
t123 = sin(qJ(5));
t126 = cos(qJ(5));
t127 = cos(qJ(4));
t154 = qJD(5) * t127;
t124 = sin(qJ(4));
t159 = qJD(4) * t124;
t131 = t123 * t159 - t126 * t154;
t207 = t123 / 0.2e1;
t190 = t126 / 0.2e1;
t163 = t123 ^ 2 + t126 ^ 2;
t121 = sin(pkin(5));
t125 = sin(qJ(2));
t168 = t121 * t125;
t111 = pkin(7) * t168;
t122 = cos(pkin(5));
t128 = cos(qJ(2));
t187 = pkin(1) * t128;
t150 = -pkin(2) - t187;
t50 = pkin(3) * t168 + t111 + (-pkin(8) + t150) * t122;
t129 = -pkin(2) - pkin(8);
t171 = qJ(3) * t125;
t61 = (t129 * t128 - pkin(1) - t171) * t121;
t181 = t124 * t50 + t127 * t61;
t206 = qJD(4) * t181;
t162 = qJD(2) * t121;
t148 = t125 * t162;
t108 = pkin(2) * t148;
t160 = qJD(3) * t125;
t47 = t108 + (-t160 + (pkin(8) * t125 - qJ(3) * t128) * qJD(2)) * t121;
t114 = t122 * t125 * pkin(1);
t167 = t121 * t128;
t196 = pkin(3) + pkin(7);
t62 = (t196 * t167 + t114) * qJD(2);
t12 = -t124 * t47 + t127 * t62 - t206;
t205 = 2 * m(6);
t204 = -0.2e1 * pkin(1);
t203 = 2 * mrSges(4,1);
t202 = -2 * mrSges(3,3);
t137 = -pkin(2) * t128 - t171;
t72 = (-pkin(1) + t137) * t121;
t201 = -0.2e1 * t72;
t151 = t124 * t167;
t80 = t122 * t127 - t151;
t58 = -t123 * t80 + t126 * t168;
t200 = t58 / 0.2e1;
t59 = t123 * t168 + t126 * t80;
t199 = t59 / 0.2e1;
t178 = Ifges(6,4) * t123;
t139 = Ifges(6,1) * t126 - t178;
t76 = Ifges(6,5) * t124 + t139 * t127;
t198 = t76 / 0.2e1;
t155 = qJD(5) * t126;
t116 = Ifges(6,5) * t155;
t156 = qJD(5) * t123;
t197 = -Ifges(6,6) * t156 / 0.2e1 + t116 / 0.2e1;
t195 = Ifges(6,5) * t207 + Ifges(6,6) * t190;
t101 = Ifges(6,2) * t126 + t178;
t194 = t101 / 0.2e1;
t193 = -t123 / 0.2e1;
t191 = -t126 / 0.2e1;
t161 = qJD(2) * t128;
t147 = t121 * t161;
t79 = t122 * t124 + t127 * t167;
t56 = -t79 * qJD(4) + t124 * t148;
t38 = mrSges(5,1) * t147 - mrSges(5,3) * t56;
t21 = -t59 * qJD(5) - t123 * t56 + t126 * t147;
t22 = t58 * qJD(5) + t123 * t147 + t126 * t56;
t6 = -mrSges(6,1) * t21 + mrSges(6,2) * t22;
t189 = t38 - t6;
t188 = m(6) * t124;
t186 = pkin(4) * t124;
t185 = pkin(9) * t127;
t153 = t122 * t187;
t109 = qJD(2) * t153;
t77 = -pkin(7) * t148 + t109;
t184 = t77 * mrSges(3,2);
t183 = Ifges(4,6) + Ifges(3,4);
t31 = -mrSges(6,1) * t58 + mrSges(6,2) * t59;
t65 = mrSges(5,1) * t168 - mrSges(5,3) * t80;
t182 = t31 - t65;
t180 = Ifges(5,4) * t124;
t179 = Ifges(5,4) * t127;
t177 = Ifges(6,4) * t126;
t176 = Ifges(6,6) * t123;
t17 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t79;
t175 = t123 * t17;
t82 = pkin(7) * t167 + t114;
t78 = t82 * qJD(2);
t174 = t125 * t78;
t18 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t79;
t173 = t126 * t18;
t172 = t127 * mrSges(6,3);
t166 = t124 * t129;
t95 = qJ(3) - t185 + t186;
t69 = -t123 * t166 + t126 * t95;
t170 = qJD(5) * t69;
t70 = t123 * t95 + t126 * t166;
t169 = qJD(5) * t70;
t165 = t127 * t129;
t164 = Ifges(3,5) * t147 + Ifges(4,5) * t148;
t158 = qJD(4) * t126;
t157 = qJD(4) * t127;
t57 = -qJD(4) * t151 + t122 * t157 - t127 * t148;
t3 = Ifges(6,5) * t22 + Ifges(6,6) * t21 + Ifges(6,3) * t57;
t152 = Ifges(5,5) * t56 - Ifges(5,6) * t57 + Ifges(5,3) * t147;
t71 = -t122 * qJ(3) - t82;
t149 = (mrSges(4,2) - mrSges(3,1)) * t78;
t146 = t129 * t157;
t84 = qJD(3) + (pkin(4) * t127 + pkin(9) * t124) * qJD(4);
t32 = t123 * t84 + t126 * t146 + t170;
t143 = t32 - t170;
t33 = -t123 * t146 + t126 * t84 - t169;
t142 = -t33 - t169;
t60 = pkin(3) * t167 - t71;
t117 = t122 * qJD(3);
t49 = -t196 * t148 + t109 + t117;
t15 = pkin(4) * t57 - pkin(9) * t56 + t49;
t11 = t124 * t62 + t127 * t47 + t50 * t157 - t61 * t159;
t7 = pkin(9) * t147 + t11;
t24 = pkin(9) * t168 + t181;
t30 = pkin(4) * t79 - pkin(9) * t80 + t60;
t9 = -t123 * t24 + t126 * t30;
t1 = t9 * qJD(5) + t123 * t15 + t126 * t7;
t10 = t123 * t30 + t126 * t24;
t2 = -t10 * qJD(5) - t123 * t7 + t126 * t15;
t141 = t1 * t126 - t123 * t2;
t99 = t124 * mrSges(5,1) + t127 * mrSges(5,2);
t98 = -mrSges(6,1) * t126 + mrSges(6,2) * t123;
t140 = mrSges(6,1) * t123 + mrSges(6,2) * t126;
t103 = Ifges(6,1) * t123 + t177;
t138 = -Ifges(6,2) * t123 + t177;
t13 = -mrSges(6,2) * t57 + mrSges(6,3) * t21;
t14 = mrSges(6,1) * t57 - mrSges(6,3) * t22;
t136 = -t123 * t14 + t126 * t13;
t27 = -t124 * t61 + t127 * t50;
t135 = t27 * t124 - t127 * t181;
t132 = t123 * t154 + t124 * t158;
t130 = -t10 * t156 - t9 * t155 + t141;
t40 = -t132 * Ifges(6,5) + t131 * Ifges(6,6) + Ifges(6,3) * t157;
t104 = Ifges(5,1) * t127 - t180;
t102 = -Ifges(5,2) * t124 + t179;
t94 = mrSges(6,1) * t124 - t126 * t172;
t93 = -mrSges(6,2) * t124 - t123 * t172;
t92 = (-Ifges(5,1) * t124 - t179) * qJD(4);
t91 = t139 * qJD(5);
t90 = (-Ifges(5,2) * t127 - t180) * qJD(4);
t89 = t138 * qJD(5);
t87 = (mrSges(5,1) * t127 - mrSges(5,2) * t124) * qJD(4);
t86 = t140 * qJD(5);
t85 = -mrSges(4,1) * t167 - mrSges(4,3) * t122;
t83 = t140 * t127;
t81 = -t111 + t153;
t75 = Ifges(6,6) * t124 + t138 * t127;
t74 = Ifges(6,3) * t124 + (Ifges(6,5) * t126 - t176) * t127;
t73 = t150 * t122 + t111;
t68 = -t117 - t77;
t67 = -mrSges(6,2) * t157 + t131 * mrSges(6,3);
t66 = mrSges(6,1) * t157 + t132 * mrSges(6,3);
t64 = -mrSges(5,2) * t168 - mrSges(5,3) * t79;
t63 = t108 + (-qJ(3) * t161 - t160) * t121;
t48 = -t131 * mrSges(6,1) - t132 * mrSges(6,2);
t44 = mrSges(5,1) * t79 + mrSges(5,2) * t80;
t42 = -t103 * t154 + (Ifges(6,5) * t127 - t139 * t124) * qJD(4);
t41 = -t101 * t154 + (Ifges(6,6) * t127 - t138 * t124) * qJD(4);
t39 = -mrSges(5,2) * t147 - mrSges(5,3) * t57;
t37 = Ifges(5,1) * t80 - Ifges(5,4) * t79 + Ifges(5,5) * t168;
t36 = Ifges(5,4) * t80 - Ifges(5,2) * t79 + Ifges(5,6) * t168;
t35 = mrSges(6,1) * t79 - mrSges(6,3) * t59;
t34 = -mrSges(6,2) * t79 + mrSges(6,3) * t58;
t29 = mrSges(5,1) * t57 + mrSges(5,2) * t56;
t26 = Ifges(5,1) * t56 - Ifges(5,4) * t57 + Ifges(5,5) * t147;
t25 = Ifges(5,4) * t56 - Ifges(5,2) * t57 + Ifges(5,6) * t147;
t23 = -pkin(4) * t168 - t27;
t16 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + Ifges(6,3) * t79;
t8 = -pkin(4) * t147 - t12;
t5 = Ifges(6,1) * t22 + Ifges(6,4) * t21 + Ifges(6,5) * t57;
t4 = Ifges(6,4) * t22 + Ifges(6,2) * t21 + Ifges(6,6) * t57;
t19 = [(t16 - t36) * t57 + (t3 - t25) * t79 + 0.2e1 * t1 * t34 + 0.2e1 * t2 * t35 + 0.2e1 * t27 * t38 + 0.2e1 * t8 * t31 + t21 * t17 + t22 * t18 + 0.2e1 * t23 * t6 + 0.2e1 * t10 * t13 + 0.2e1 * t9 * t14 + 0.2e1 * m(3) * (t77 * t82 - t78 * t81) + 0.2e1 * m(4) * (t63 * t72 + t68 * t71 + t73 * t78) + (t1 * t10 + t2 * t9 + t23 * t8) * t205 + (t125 * t152 + t174 * t203 + 0.2e1 * t63 * (mrSges(4,2) * t128 - mrSges(4,3) * t125) + 0.2e1 * (t128 * t77 + t174) * mrSges(3,3) + ((t71 * t203 + mrSges(4,2) * t201 + t82 * t202 + (Ifges(4,5) - (2 * Ifges(3,6))) * t122 + (mrSges(3,1) * t204 - 0.2e1 * t125 * t183) * t121) * t125 + (t73 * t203 + t81 * t202 + mrSges(4,3) * t201 + Ifges(5,5) * t80 - Ifges(5,6) * t79 + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t122 + (mrSges(3,2) * t204 + 0.2e1 * t128 * t183) * t121 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) + (2 * Ifges(4,2)) - (2 * Ifges(4,3)) + Ifges(5,3)) * t168) * t128) * qJD(2)) * t121 + (0.2e1 * t149 + t164 - 0.2e1 * t184) * t122 + 0.2e1 * t181 * t39 + 0.2e1 * m(5) * (t11 * t181 + t12 * t27 + t49 * t60) + 0.2e1 * t49 * t44 + t56 * t37 + t58 * t4 + t59 * t5 + 0.2e1 * t60 * t29 + 0.2e1 * t11 * t64 + 0.2e1 * t12 * t65 + t80 * t26 + 0.2e1 * t68 * t85; t32 * t34 + t33 * t35 + qJ(3) * t29 + t149 + (-t11 * mrSges(5,3) + t129 * t39 - t25 / 0.2e1 + t3 / 0.2e1) * t124 + (t44 - t85) * qJD(3) + t164 + t22 * t198 + t42 * t199 + t41 * t200 + (t128 * (Ifges(5,5) * t127 - Ifges(5,6) * t124) / 0.2e1 - t128 * Ifges(4,4) - t125 * Ifges(3,6) + t137 * mrSges(4,1)) * t162 - t184 + (t5 * t190 + t4 * t193 - t12 * mrSges(5,3) + t26 / 0.2e1 + t189 * t129 + (t17 * t191 + t18 * t193) * qJD(5)) * t127 + m(6) * (t70 * t1 + t32 * t10 - t8 * t165 + t69 * t2 + t33 * t9) + m(5) * (qJ(3) * t49 + qJD(3) * t60 + t11 * t166 + t12 * t165) + m(4) * (-pkin(2) * t78 - qJ(3) * t68 - qJD(3) * t71) + (t40 / 0.2e1 - t90 / 0.2e1) * t79 + (t74 / 0.2e1 - t102 / 0.2e1) * t57 + (t124 * t175 / 0.2e1 + (-Ifges(5,5) * t124 - Ifges(5,6) * t127) * t168 / 0.2e1 - t127 * t36 / 0.2e1 + t127 * t16 / 0.2e1 + t135 * mrSges(5,3) + (-m(5) * t135 + t182 * t124 + t127 * t64 + t23 * t188) * t129 - (t37 + t173) * t124 / 0.2e1) * qJD(4) + t23 * t48 + t9 * t66 + t10 * t67 - t68 * mrSges(4,3) + t69 * t14 + t70 * t13 + t21 * t75 / 0.2e1 + t8 * t83 + t60 * t87 + t80 * t92 / 0.2e1 + t1 * t93 + t2 * t94 + t49 * t99 + t56 * t104 / 0.2e1; (t70 * t32 + t69 * t33) * t205 + 0.2e1 * t32 * t93 + 0.2e1 * t70 * t67 + 0.2e1 * t33 * t94 + 0.2e1 * t69 * t66 + 0.2e1 * qJ(3) * t87 + 0.2e1 * (mrSges(4,3) + t99 + (m(4) + m(5)) * qJ(3)) * qJD(3) + (t40 - t90 + (t123 * t75 - t126 * t76 + 0.2e1 * t129 * t83 - t104) * qJD(4)) * t124 + (-t123 * t41 + t126 * t42 - 0.2e1 * t129 * t48 + t92 + (-t123 * t76 - t126 * t75) * qJD(5) + (-0.2e1 * t129 ^ 2 * t188 - t102 + t74) * qJD(4)) * t127; m(4) * t78 + mrSges(4,1) * t147 + ((-t123 * t35 + t126 * t34 + t64) * qJD(4) + m(6) * (-qJD(4) * t123 * t9 + t10 * t158 - t8) + m(5) * (t12 + t206) + t189) * t127 + (t39 + (-t123 * t34 - t126 * t35) * qJD(5) + t182 * qJD(4) + m(6) * (qJD(4) * t23 + t130) + m(5) * (-qJD(4) * t27 + t11) + t136) * t124; (-t48 + (m(6) * (-t123 * t69 + t126 * t70) + t126 * t93 - t123 * t94) * qJD(4)) * t127 + (m(6) * (-t123 * t33 + t126 * t32 - t69 * t155 - t70 * t156 - 0.2e1 * t146) - t93 * t156 + t126 * t67 - t94 * t155 - t123 * t66 + qJD(4) * t83) * t124; 0.2e1 * (-0.1e1 + t163) * t157 * t188; -t11 * mrSges(5,2) + t12 * mrSges(5,1) + t23 * t86 + t79 * t197 + t89 * t200 + t91 * t199 + t8 * t98 + t57 * t195 + t21 * t194 + t22 * t103 / 0.2e1 + t5 * t207 + t4 * t190 + (t173 / 0.2e1 - t175 / 0.2e1) * qJD(5) + (-m(6) * t8 - t6) * pkin(4) + ((-t10 * t123 - t126 * t9) * qJD(5) + t141) * mrSges(6,3) + (m(6) * t130 - t35 * t155 - t34 * t156 + t136) * pkin(9) + t152; -pkin(4) * t48 + (t197 + (-Ifges(5,5) + (-m(6) * pkin(4) - mrSges(5,1) + t98) * t129) * qJD(4)) * t124 + (qJD(5) * t198 - t103 * t159 / 0.2e1 + t41 / 0.2e1 + t143 * mrSges(6,3) + (m(6) * t143 - qJD(5) * t94 + t67) * pkin(9)) * t126 + (-qJD(5) * t75 / 0.2e1 + t159 * t194 + t42 / 0.2e1 + t142 * mrSges(6,3) + (m(6) * t142 - qJD(5) * t93 - t66) * pkin(9)) * t123 + (t91 * t190 + t89 * t193 - t129 * t86 + (t101 * t191 + t103 * t193) * qJD(5) + (-t129 * mrSges(5,2) - Ifges(5,6) + t195) * qJD(4)) * t127; -t127 * t86 + (t124 * t98 + m(6) * (t163 * t185 - t186) + t163 * t172 - t99) * qJD(4); -0.2e1 * pkin(4) * t86 + t123 * t91 + t126 * t89 + (-t101 * t123 + t103 * t126) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 + t3; mrSges(6,1) * t33 - mrSges(6,2) * t32 + t40; (t124 * t156 - t126 * t157) * mrSges(6,2) + (-t123 * t157 - t124 * t155) * mrSges(6,1); t116 + (pkin(9) * t98 - t176) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
