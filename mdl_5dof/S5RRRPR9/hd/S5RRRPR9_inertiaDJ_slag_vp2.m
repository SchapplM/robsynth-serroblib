% Calculate time derivative of joint inertia matrix for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR9_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR9_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR9_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR9_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:22:16
% EndTime: 2019-12-31 21:22:24
% DurationCPUTime: 2.95s
% Computational Cost: add. (3695->377), mult. (9079->578), div. (0->0), fcn. (7927->8), ass. (0->171)
t218 = -Ifges(4,3) - Ifges(5,3);
t156 = sin(pkin(9));
t157 = cos(pkin(9));
t159 = sin(qJ(3));
t162 = cos(qJ(3));
t124 = t156 * t162 + t157 * t159;
t118 = t124 * qJD(3);
t167 = t156 * t159 - t157 * t162;
t119 = t167 * qJD(3);
t185 = qJD(3) * t162;
t217 = Ifges(4,5) * t185 - Ifges(5,5) * t119 - Ifges(5,6) * t118;
t160 = sin(qJ(2));
t163 = cos(qJ(2));
t188 = qJD(2) * t163;
t164 = t159 * t188 + t160 * t185;
t137 = -pkin(2) * t163 - t160 * pkin(7) - pkin(1);
t193 = t162 * t163;
t147 = pkin(6) * t193;
t184 = qJD(4) * t162;
t135 = (pkin(2) * t160 - pkin(7) * t163) * qJD(2);
t189 = qJD(2) * t160;
t205 = pkin(6) * t159;
t190 = t162 * t135 + t189 * t205;
t41 = -t160 * t184 + (pkin(3) * t160 - qJ(4) * t193) * qJD(2) + (-t147 + (qJ(4) * t160 - t137) * t159) * qJD(3) + t190;
t191 = t159 * t135 + t137 * t185;
t196 = t160 * t162;
t48 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t196 + (-qJD(4) * t160 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t163) * t159 + t191;
t14 = -t156 * t48 + t157 * t41;
t68 = -t118 * t160 - t167 * t188;
t10 = pkin(4) * t189 - pkin(8) * t68 + t14;
t15 = t156 * t41 + t157 * t48;
t186 = qJD(3) * t160;
t67 = -t124 * t188 + t167 * t186;
t11 = pkin(8) * t67 + t15;
t158 = sin(qJ(5));
t161 = cos(qJ(5));
t110 = t167 * t160;
t126 = t162 * t137;
t84 = -qJ(4) * t196 + t126 + (-pkin(3) - t205) * t163;
t100 = t159 * t137 + t147;
t197 = t159 * t160;
t90 = -qJ(4) * t197 + t100;
t49 = -t156 * t90 + t157 * t84;
t31 = -pkin(4) * t163 + t110 * pkin(8) + t49;
t109 = t124 * t160;
t50 = t156 * t84 + t157 * t90;
t34 = -pkin(8) * t109 + t50;
t12 = -t158 * t34 + t161 * t31;
t2 = qJD(5) * t12 + t10 * t158 + t11 * t161;
t13 = t158 * t31 + t161 * t34;
t3 = -qJD(5) * t13 + t10 * t161 - t11 * t158;
t216 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t215 = 2 * m(4);
t214 = 2 * m(5);
t213 = 2 * m(6);
t212 = -0.2e1 * pkin(1);
t211 = 0.2e1 * pkin(6);
t209 = -t159 / 0.2e1;
t208 = t162 / 0.2e1;
t206 = pkin(3) * t156;
t203 = -qJ(4) - pkin(7);
t82 = -t124 * t158 - t161 * t167;
t42 = qJD(5) * t82 - t118 * t158 - t119 * t161;
t83 = t124 * t161 - t158 * t167;
t43 = -qJD(5) * t83 - t118 * t161 + t119 * t158;
t202 = Ifges(6,5) * t42 + Ifges(6,6) * t43;
t201 = Ifges(4,4) * t159;
t200 = Ifges(4,4) * t162;
t199 = Ifges(4,6) * t159;
t198 = t163 * Ifges(4,6);
t171 = Ifges(4,1) * t162 - t201;
t108 = -Ifges(4,5) * t163 + t160 * t171;
t195 = t162 * t108;
t141 = Ifges(4,1) * t159 + t200;
t194 = t162 * t141;
t176 = qJD(3) * t203;
t116 = t159 * t176 + t184;
t117 = -qJD(4) * t159 + t162 * t176;
t70 = t157 * t116 + t156 * t117;
t138 = t203 * t159;
t139 = t203 * t162;
t92 = t156 * t138 - t157 * t139;
t136 = pkin(3) * t197 + t160 * pkin(6);
t187 = qJD(3) * t159;
t62 = -t109 * t161 + t110 * t158;
t24 = qJD(5) * t62 + t158 * t67 + t161 * t68;
t63 = -t109 * t158 - t110 * t161;
t25 = -qJD(5) * t63 - t158 * t68 + t161 * t67;
t183 = -Ifges(6,5) * t24 - Ifges(6,6) * t25 - Ifges(6,3) * t189;
t182 = pkin(3) * t187;
t98 = t164 * pkin(3) + pkin(6) * t188;
t149 = -pkin(3) * t162 - pkin(2);
t181 = t159 * t186;
t178 = t162 * t188;
t35 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t6 = -t25 * mrSges(6,1) + t24 * mrSges(6,2);
t16 = -t43 * mrSges(6,1) + t42 * mrSges(6,2);
t177 = (2 * Ifges(3,4)) + t199;
t76 = t118 * mrSges(5,1) - t119 * mrSges(5,2);
t69 = -t116 * t156 + t157 * t117;
t91 = t157 * t138 + t139 * t156;
t71 = -pkin(8) * t124 + t91;
t72 = -pkin(8) * t167 + t92;
t32 = -t158 * t72 + t161 * t71;
t52 = pkin(8) * t119 + t69;
t53 = -pkin(8) * t118 + t70;
t8 = qJD(5) * t32 + t158 * t52 + t161 * t53;
t33 = t158 * t71 + t161 * t72;
t9 = -qJD(5) * t33 - t158 * t53 + t161 * t52;
t175 = t9 * mrSges(6,1) - t8 * mrSges(6,2) + t202;
t174 = -mrSges(4,1) * t162 + mrSges(4,2) * t159;
t173 = mrSges(4,1) * t159 + mrSges(4,2) * t162;
t148 = pkin(3) * t157 + pkin(4);
t114 = t148 * t161 - t158 * t206;
t105 = t114 * qJD(5);
t115 = t148 * t158 + t161 * t206;
t106 = t115 * qJD(5);
t172 = -t106 * mrSges(6,1) - t105 * mrSges(6,2);
t170 = -Ifges(4,2) * t159 + t200;
t140 = Ifges(4,2) * t162 + t201;
t169 = Ifges(4,5) * t159 + Ifges(4,6) * t162;
t58 = (-t162 * t189 - t163 * t187) * pkin(6) + t191;
t59 = -t100 * qJD(3) + t190;
t168 = -t159 * t59 + t162 * t58;
t166 = -Ifges(4,5) * t178 - Ifges(5,5) * t68 - Ifges(5,6) * t67 + t218 * t189 + t183;
t165 = t178 - t181;
t134 = -mrSges(4,1) * t163 - mrSges(4,3) * t196;
t133 = mrSges(4,2) * t163 - mrSges(4,3) * t197;
t132 = t171 * qJD(3);
t131 = t170 * qJD(3);
t130 = t173 * qJD(3);
t107 = t160 * t170 - t198;
t101 = pkin(4) * t167 + t149;
t99 = -t163 * t205 + t126;
t97 = -mrSges(4,2) * t189 - mrSges(4,3) * t164;
t96 = mrSges(4,1) * t189 - mrSges(4,3) * t165;
t95 = pkin(4) * t118 + t182;
t94 = -mrSges(5,1) * t163 + t110 * mrSges(5,3);
t93 = mrSges(5,2) * t163 - t109 * mrSges(5,3);
t89 = Ifges(5,1) * t124 - Ifges(5,4) * t167;
t88 = Ifges(5,4) * t124 - Ifges(5,2) * t167;
t87 = mrSges(5,1) * t167 + mrSges(5,2) * t124;
t85 = pkin(4) * t109 + t136;
t81 = mrSges(4,1) * t164 + mrSges(4,2) * t165;
t78 = -Ifges(5,1) * t119 - Ifges(5,4) * t118;
t77 = -Ifges(5,4) * t119 - Ifges(5,2) * t118;
t75 = -t141 * t186 + (Ifges(4,5) * t160 + t163 * t171) * qJD(2);
t74 = -t140 * t186 + (Ifges(4,6) * t160 + t163 * t170) * qJD(2);
t73 = mrSges(5,1) * t109 - mrSges(5,2) * t110;
t61 = -Ifges(5,1) * t110 - Ifges(5,4) * t109 - Ifges(5,5) * t163;
t60 = -Ifges(5,4) * t110 - Ifges(5,2) * t109 - Ifges(5,6) * t163;
t57 = mrSges(5,1) * t189 - mrSges(5,3) * t68;
t56 = -mrSges(5,2) * t189 + mrSges(5,3) * t67;
t55 = -mrSges(6,1) * t163 - t63 * mrSges(6,3);
t54 = mrSges(6,2) * t163 + t62 * mrSges(6,3);
t51 = -pkin(4) * t67 + t98;
t47 = Ifges(6,1) * t83 + Ifges(6,4) * t82;
t46 = Ifges(6,4) * t83 + Ifges(6,2) * t82;
t45 = -mrSges(6,1) * t82 + mrSges(6,2) * t83;
t30 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t29 = Ifges(5,1) * t68 + Ifges(5,4) * t67 + Ifges(5,5) * t189;
t28 = Ifges(5,4) * t68 + Ifges(5,2) * t67 + Ifges(5,6) * t189;
t27 = Ifges(6,1) * t63 + Ifges(6,4) * t62 - Ifges(6,5) * t163;
t26 = Ifges(6,4) * t63 + Ifges(6,2) * t62 - Ifges(6,6) * t163;
t20 = -mrSges(6,2) * t189 + mrSges(6,3) * t25;
t19 = mrSges(6,1) * t189 - mrSges(6,3) * t24;
t18 = Ifges(6,1) * t42 + Ifges(6,4) * t43;
t17 = Ifges(6,4) * t42 + Ifges(6,2) * t43;
t5 = Ifges(6,1) * t24 + Ifges(6,4) * t25 + Ifges(6,5) * t189;
t4 = Ifges(6,4) * t24 + Ifges(6,2) * t25 + Ifges(6,6) * t189;
t1 = [(t100 * t58 + t99 * t59) * t215 + (t12 * t3 + t13 * t2 + t51 * t85) * t213 + (t136 * t98 + t14 * t49 + t15 * t50) * t214 + 0.2e1 * t12 * t19 + 0.2e1 * t13 * t20 + t25 * t26 + t24 * t27 + 0.2e1 * t51 * t30 + 0.2e1 * t2 * t54 + 0.2e1 * t3 * t55 + 0.2e1 * t50 * t56 + 0.2e1 * t49 * t57 + t62 * t4 + t63 * t5 + t67 * t60 + t68 * t61 + 0.2e1 * t85 * t6 + 0.2e1 * t15 * t93 + 0.2e1 * t14 * t94 + 0.2e1 * t98 * t73 + 0.2e1 * t99 * t96 + 0.2e1 * t100 * t97 - t109 * t28 - t110 * t29 + 0.2e1 * t58 * t133 + 0.2e1 * t59 * t134 + 0.2e1 * t136 * t35 + ((mrSges(3,2) * t212 - t159 * t107 + t163 * t177 + t195) * qJD(2) + t166) * t163 + (t81 * t211 - t159 * t74 + t162 * t75 + (-t162 * t107 - t159 * t108 + t163 * t169) * qJD(3) + (Ifges(6,5) * t63 + Ifges(6,6) * t62 - Ifges(5,5) * t110 - Ifges(5,6) * t109 + mrSges(3,1) * t212 + (Ifges(4,5) * t162 - t177) * t160 + (pkin(6) ^ 2 * t215 + t173 * t211 + (2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(6,3) + t218) * t163) * qJD(2)) * t160; (t131 * t209 - Ifges(3,6) * qJD(2) + t132 * t208 + (-t162 * t140 / 0.2e1 + t141 * t209) * qJD(3) + (qJD(2) * mrSges(3,2) + t130) * pkin(6) + (Ifges(5,5) * t124 + Ifges(6,5) * t83 - Ifges(5,6) * t167 + Ifges(6,6) * t82 + t169) * qJD(2) / 0.2e1) * t160 + (-t50 * t118 + t49 * t119 - t14 * t124 - t15 * t167) * mrSges(5,3) - t167 * t28 / 0.2e1 + m(6) * (t101 * t51 + t12 * t9 + t13 * t8 + t2 * t33 + t3 * t32 + t85 * t95) + t74 * t208 + ((-t100 * t159 - t162 * t99) * qJD(3) + t168) * mrSges(4,3) + t32 * t19 + t33 * t20 + t42 * t27 / 0.2e1 + t43 * t26 / 0.2e1 + t25 * t46 / 0.2e1 + t24 * t47 / 0.2e1 + t51 * t45 + t8 * t54 + t9 * t55 + t62 * t17 / 0.2e1 + t63 * t18 / 0.2e1 - pkin(2) * t81 + t82 * t4 / 0.2e1 + t83 * t5 / 0.2e1 + t85 * t16 + t67 * t88 / 0.2e1 + t68 * t89 / 0.2e1 + t91 * t57 + t92 * t56 + t70 * t93 + t69 * t94 + t95 * t30 + t98 * t87 + m(5) * (t136 * t182 + t14 * t91 + t149 * t98 + t15 * t92 + t49 * t69 + t50 * t70) + t101 * t6 - t109 * t77 / 0.2e1 - t110 * t78 / 0.2e1 - t118 * t60 / 0.2e1 + (-t134 * t185 - t133 * t187 - t159 * t96 + t162 * t97 + m(4) * (-t100 * t187 - t185 * t99 + t168)) * pkin(7) - t119 * t61 / 0.2e1 + t124 * t29 / 0.2e1 + t136 * t76 + t149 * t35 + (t195 / 0.2e1 + (pkin(3) * t73 - t107 / 0.2e1 + t198 / 0.2e1) * t159) * qJD(3) + t159 * t75 / 0.2e1 + (t194 / 0.2e1 + t140 * t209 + Ifges(3,5) + (-m(4) * pkin(2) - mrSges(3,1) + t174) * pkin(6)) * t188 + (-t12 * t42 + t13 * t43 + t2 * t82 - t3 * t83) * mrSges(6,3) - (t202 + t217) * t163 / 0.2e1; -0.2e1 * pkin(2) * t130 + 0.2e1 * t101 * t16 - t118 * t88 - t119 * t89 - t167 * t77 + t124 * t78 + t162 * t131 + t159 * t132 + 0.2e1 * t149 * t76 + t82 * t17 + t83 * t18 + t42 * t47 + t43 * t46 + 0.2e1 * t95 * t45 + (t194 + (0.2e1 * pkin(3) * t87 - t140) * t159) * qJD(3) + (t101 * t95 + t32 * t9 + t33 * t8) * t213 + (t149 * t182 + t69 * t91 + t70 * t92) * t214 + 0.2e1 * (-t32 * t42 + t33 * t43 + t8 * t82 - t83 * t9) * mrSges(6,3) + 0.2e1 * (-t118 * t92 + t119 * t91 - t124 * t69 - t167 * t70) * mrSges(5,3); -Ifges(4,5) * t181 - t166 + (t156 * t56 + t157 * t57 + m(5) * (t14 * t157 + t15 * t156)) * pkin(3) + m(6) * (t105 * t13 - t106 * t12 + t114 * t3 + t115 * t2) - t164 * Ifges(4,6) + t14 * mrSges(5,1) - t15 * mrSges(5,2) - t58 * mrSges(4,2) + t59 * mrSges(4,1) + t105 * t54 - t106 * t55 + t114 * t19 + t115 * t20 + t216; m(6) * (t105 * t33 - t106 * t32 + t114 * t9 + t115 * t8) - t70 * mrSges(5,2) + t69 * mrSges(5,1) + (pkin(7) * t174 - t199) * qJD(3) + (m(5) * (t156 * t70 + t157 * t69) + (-t118 * t156 + t119 * t157) * mrSges(5,3)) * pkin(3) + (t105 * t82 + t106 * t83 - t114 * t42 + t115 * t43) * mrSges(6,3) + t175 + t217; 0.2e1 * m(6) * (t105 * t115 - t106 * t114) + 0.2e1 * t172; m(5) * t98 + m(6) * t51 + t35 + t6; m(5) * t182 + m(6) * t95 + t16 + t76; 0; 0; -t183 + t216; t175; t172; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
