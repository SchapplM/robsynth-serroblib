% Calculate time derivative of joint inertia matrix for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:30
% EndTime: 2019-03-09 02:25:37
% DurationCPUTime: 3.04s
% Computational Cost: add. (3276->393), mult. (6692->589), div. (0->0), fcn. (5677->8), ass. (0->173)
t200 = qJD(5) + qJD(6);
t116 = sin(qJ(6));
t117 = sin(qJ(5));
t119 = cos(qJ(6));
t120 = cos(qJ(5));
t173 = t119 * t120;
t78 = -t116 * t117 + t173;
t44 = t200 * t78;
t79 = t116 * t120 + t117 * t119;
t45 = t200 * t79;
t18 = mrSges(7,1) * t45 + mrSges(7,2) * t44;
t131 = mrSges(6,1) * t117 + mrSges(6,2) * t120;
t80 = t131 * qJD(5);
t208 = t18 + t80;
t121 = cos(qJ(4));
t165 = qJD(4) * t121;
t141 = t120 * t165;
t147 = t117 * t165;
t118 = sin(qJ(4));
t69 = t79 * t118;
t24 = t116 * t147 - t119 * t141 + t200 * t69;
t176 = t117 * t118;
t154 = t116 * t176;
t174 = t118 * t120;
t163 = qJD(5) * t118;
t144 = t117 * t163;
t201 = t144 - t141;
t25 = -qJD(6) * t154 + (t174 * t200 + t147) * t119 - t201 * t116;
t203 = -t25 * mrSges(7,1) + t24 * mrSges(7,2);
t162 = qJD(5) * t120;
t125 = t118 * t162 + t147;
t41 = -mrSges(6,1) * t125 + mrSges(6,2) * t201;
t207 = t41 + t203;
t206 = Ifges(5,1) - Ifges(5,2);
t205 = -Ifges(6,3) - Ifges(7,3);
t166 = qJD(4) * t118;
t148 = t117 * t166;
t115 = cos(pkin(10));
t168 = qJD(2) * t115;
t149 = t121 * t168;
t114 = sin(pkin(10));
t169 = qJD(2) * t114;
t189 = pkin(8) * t121;
t190 = pkin(4) * t118;
t71 = t169 + (t189 - t190) * qJD(4);
t122 = -pkin(1) - pkin(2);
t171 = t115 * qJ(2) + t114 * t122;
t77 = -pkin(7) + t171;
t126 = -t117 * t149 + t120 * t71 + t77 * t148;
t172 = t120 * t121;
t134 = -t114 * qJ(2) + t115 * t122;
t76 = pkin(3) - t134;
t61 = t121 * pkin(4) + t118 * pkin(8) + t76;
t68 = t77 * t172;
t7 = (-pkin(5) * t118 + pkin(9) * t172) * qJD(4) + (-t68 + (-pkin(9) * t118 - t61) * t117) * qJD(5) + t126;
t164 = qJD(5) * t117;
t140 = t121 * t164;
t146 = t120 * t166;
t14 = t61 * t162 + t117 * t71 + t120 * t149 + (-t140 - t146) * t77;
t8 = pkin(9) * t125 + t14;
t57 = t120 * t61;
t29 = pkin(9) * t174 + t57 + (-t117 * t77 + pkin(5)) * t121;
t34 = t117 * t61 + t68;
t30 = pkin(9) * t176 + t34;
t9 = -t116 * t30 + t119 * t29;
t2 = qJD(6) * t9 + t116 * t7 + t119 * t8;
t10 = t116 * t29 + t119 * t30;
t3 = -qJD(6) * t10 - t116 * t8 + t119 * t7;
t204 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t202 = qJD(5) * t34;
t170 = t117 ^ 2 + t120 ^ 2;
t199 = 2 * m(7);
t197 = m(6) * pkin(4);
t196 = m(7) * pkin(5);
t183 = Ifges(6,4) * t117;
t95 = Ifges(6,2) * t120 + t183;
t195 = t95 / 0.2e1;
t194 = -pkin(9) - pkin(8);
t192 = t117 / 0.2e1;
t191 = -t120 / 0.2e1;
t187 = Ifges(7,5) * t24 + Ifges(7,6) * t25;
t94 = -mrSges(6,1) * t120 + mrSges(6,2) * t117;
t186 = t94 - mrSges(5,1);
t182 = Ifges(6,4) * t120;
t181 = Ifges(6,6) * t117;
t111 = t118 ^ 2;
t180 = t111 * t77;
t113 = t121 ^ 2;
t179 = t113 * t77;
t178 = t118 * mrSges(5,2);
t177 = t121 * Ifges(6,6);
t175 = t117 * t121;
t161 = qJD(6) * t116;
t160 = qJD(6) * t119;
t46 = -mrSges(7,1) * t78 + mrSges(7,2) * t79;
t159 = t46 + t186;
t158 = 0.2e1 * t169;
t157 = Ifges(6,5) * t172;
t156 = pkin(5) * t164;
t155 = t77 * t165;
t153 = Ifges(6,5) * t144 + Ifges(6,6) * t125;
t70 = t118 * t173 - t154;
t38 = -mrSges(7,1) * t69 - mrSges(7,2) * t70;
t58 = (-pkin(5) * t117 + t77) * t118;
t152 = m(7) * t58 + t38;
t151 = qJD(5) * t194;
t150 = t114 * t168;
t145 = t118 * t165;
t142 = t118 * t168;
t72 = -t114 * t175 - t115 * t120;
t73 = t114 * t172 - t115 * t117;
t39 = -t116 * t73 + t119 * t72;
t49 = qJD(5) * t72 - t114 * t146;
t50 = -qJD(5) * t73 + t114 * t148;
t12 = qJD(6) * t39 + t116 * t50 + t119 * t49;
t40 = t116 * t72 + t119 * t73;
t13 = -qJD(6) * t40 - t116 * t49 + t119 * t50;
t139 = t13 * mrSges(7,1) - t12 * mrSges(7,2);
t33 = -t175 * t77 + t57;
t137 = -qJD(5) * t33 + t14;
t135 = 0.2e1 * t155;
t133 = t186 - t197;
t132 = (t111 - t113) * qJD(4) * t114;
t130 = Ifges(6,1) * t120 - t183;
t96 = Ifges(6,1) * t117 + t182;
t129 = -Ifges(6,2) * t117 + t182;
t97 = t194 * t117;
t98 = t194 * t120;
t52 = t116 * t98 + t119 * t97;
t53 = t116 * t97 - t119 * t98;
t86 = t117 * t151;
t87 = t120 * t151;
t27 = qJD(6) * t52 + t116 * t87 + t119 * t86;
t28 = -qJD(6) * t53 - t116 * t86 + t119 * t87;
t42 = Ifges(7,6) * t45;
t43 = Ifges(7,5) * t44;
t128 = t28 * mrSges(7,1) - t27 * mrSges(7,2) - t42 + t43;
t127 = -Ifges(7,3) * t166 + t187;
t124 = -t117 * t50 + t120 * t49 + (-t117 * t73 - t120 * t72) * qJD(5);
t123 = m(6) * t124;
t107 = Ifges(6,5) * t162;
t106 = -pkin(5) * t120 - pkin(4);
t89 = t114 ^ 2 * t145;
t88 = t111 * t150;
t85 = mrSges(6,1) * t121 + mrSges(6,3) * t174;
t84 = -mrSges(6,2) * t121 + mrSges(6,3) * t176;
t83 = t130 * qJD(5);
t82 = t129 * qJD(5);
t81 = (-t118 * mrSges(5,1) - t121 * mrSges(5,2)) * qJD(4);
t75 = (-mrSges(7,1) * t116 - mrSges(7,2) * t119) * qJD(6) * pkin(5);
t74 = t131 * t118;
t67 = Ifges(6,5) * t121 - t118 * t130;
t66 = -t118 * t129 + t177;
t62 = t168 * t180;
t60 = mrSges(6,2) * t166 + mrSges(6,3) * t125;
t59 = -mrSges(6,1) * t166 - mrSges(6,3) * t201;
t55 = mrSges(7,1) * t121 + mrSges(7,3) * t70;
t54 = -mrSges(7,2) * t121 + mrSges(7,3) * t69;
t48 = Ifges(7,1) * t79 + Ifges(7,4) * t78;
t47 = Ifges(7,4) * t79 + Ifges(7,2) * t78;
t37 = t96 * t163 + (-Ifges(6,5) * t118 - t121 * t130) * qJD(4);
t36 = t95 * t163 + (-Ifges(6,6) * t118 - t121 * t129) * qJD(4);
t35 = -pkin(5) * t125 + t142 + t155;
t32 = -Ifges(7,1) * t70 + Ifges(7,4) * t69 + Ifges(7,5) * t121;
t31 = -Ifges(7,4) * t70 + Ifges(7,2) * t69 + Ifges(7,6) * t121;
t20 = Ifges(7,1) * t44 - Ifges(7,4) * t45;
t19 = Ifges(7,4) * t44 - Ifges(7,2) * t45;
t17 = mrSges(7,2) * t166 + mrSges(7,3) * t25;
t16 = -mrSges(7,1) * t166 - mrSges(7,3) * t24;
t15 = t126 - t202;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t25 - Ifges(7,5) * t166;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t25 - Ifges(7,6) * t166;
t1 = [(t10 * t2 + t3 * t9 + t35 * t58) * t199 + t36 * t176 + 0.2e1 * m(5) * (t62 + (t114 * t76 + t115 * t179) * qJD(2)) - t37 * t174 - (0.2e1 * t142 + t135) * t74 + (-0.2e1 * (t113 + t111) * mrSges(5,3) + (2 * mrSges(4,2))) * t168 + (-t178 + mrSges(4,1)) * t158 + 0.2e1 * m(6) * (t145 * t77 ^ 2 + t14 * t34 + t15 * t33 + t62) + (0.2e1 * Ifges(5,4) * t165 + mrSges(5,1) * t158 + (-Ifges(6,3) * t118 - t157) * qJD(4) + t153 + t127 + (t205 + t206) * t166) * t121 + 0.2e1 * (mrSges(3,3) + m(4) * (-t114 * t134 + t115 * t171) + m(3) * qJ(2)) * qJD(2) + (Ifges(7,5) * t70 - Ifges(7,6) * t69 + (Ifges(6,5) * t120 - 0.2e1 * Ifges(5,4) - t181) * t118) * t166 + (t206 * t165 + 0.2e1 * t77 * t41) * t118 + t125 * t66 + t201 * t67 + 0.2e1 * t9 * t16 + 0.2e1 * t10 * t17 + t25 * t31 + t24 * t32 + 0.2e1 * t35 * t38 + 0.2e1 * t2 * t54 + 0.2e1 * t3 * t55 + 0.2e1 * t33 * t59 + 0.2e1 * t34 * t60 + t69 * t4 - t70 * t5 + 0.2e1 * t76 * t81 + 0.2e1 * t14 * t84 + 0.2e1 * t15 * t85 + 0.2e1 * t58 * t203; -t115 * t81 + t12 * t54 + t13 * t55 + t39 * t16 + t40 * t17 + t49 * t84 + t50 * t85 + t72 * t59 + t73 * t60 + (t207 * t118 + (t38 - t74) * t165) * t114 + m(7) * (t10 * t12 + t13 * t9 + t2 * t40 + t3 * t39 + (t118 * t35 + t165 * t58) * t114) + m(6) * (t114 * t118 * t135 + t14 * t73 + t15 * t72 + t33 * t50 + t34 * t49 + t88) + m(5) * (t88 + (t113 - 0.1e1) * t150); 0.2e1 * m(7) * (t12 * t40 + t13 * t39 + t89) + 0.2e1 * m(6) * (t49 * t73 + t50 * t72 + t89); m(7) * (-t10 * t24 + t2 * t70 - t25 * t9 - t3 * t69) - t24 * t54 + t70 * t17 - t25 * t55 - t69 * t16 + (m(6) * (t172 * t34 - t175 * t33 - t179 + t180) + t84 * t172 - t85 * t175) * qJD(4) + (m(6) * (-t117 * t15 + t120 * t14 - t162 * t33 - t164 * t34 - t149) - t84 * t164 + t120 * t60 - t85 * t162 - t117 * t59 + (t152 - t74) * qJD(4)) * t118 + (-m(7) * t35 - t207) * t121; m(7) * (t12 * t70 - t13 * t69 - t24 * t40 - t25 * t39 + t132) + m(6) * (t141 * t73 - t147 * t72 + t132) + t118 * t123; (-t24 * t70 + t25 * t69) * t199 + 0.4e1 * (m(6) * (-0.1e1 + t170) / 0.2e1 - m(7) / 0.2e1) * t145; m(7) * (t10 * t27 + t106 * t35 + t2 * t53 + t28 * t9 + t3 * t52) - pkin(4) * t41 + t44 * t32 / 0.2e1 - t45 * t31 / 0.2e1 + t35 * t46 + t25 * t47 / 0.2e1 + t24 * t48 / 0.2e1 + t52 * t16 + t53 * t17 + t27 * t54 + t28 * t55 + t58 * t18 + t69 * t19 / 0.2e1 - t70 * t20 / 0.2e1 + t78 * t4 / 0.2e1 + t79 * t5 / 0.2e1 + t106 * t203 + (t43 / 0.2e1 - t42 / 0.2e1 + t107 / 0.2e1 - mrSges(5,2) * t168 + (t133 * t77 - Ifges(5,5)) * qJD(4)) * t121 + (-t15 * mrSges(6,3) + t37 / 0.2e1 + t165 * t195 + (-t66 / 0.2e1 - t177 / 0.2e1 - t34 * mrSges(6,3) + t152 * pkin(5)) * qJD(5) + (m(6) * (-t15 - t202) - qJD(5) * t84 - t59) * pkin(8)) * t117 + (-t10 * t45 + t2 * t78 - t3 * t79 - t44 * t9) * mrSges(7,3) + (qJD(5) * t67 / 0.2e1 + t36 / 0.2e1 - t96 * t165 / 0.2e1 + t137 * mrSges(6,3) + (m(6) * t137 - qJD(5) * t85 + t60) * pkin(8)) * t120 + (t77 * t80 + t82 * t192 + t83 * t191 + (t120 * t195 + t192 * t96) * qJD(5) + (-Ifges(7,5) * t79 / 0.2e1 - Ifges(7,6) * t78 / 0.2e1 + Ifges(5,6) - Ifges(6,5) * t117 / 0.2e1 + Ifges(6,6) * t191 + t77 * mrSges(5,2)) * qJD(4) + t133 * t168) * t118; m(7) * (t12 * t53 + t13 * t52 + t27 * t40 + t28 * t39) + pkin(8) * t123 + (t12 * t78 - t13 * t79 - t39 * t44 - t40 * t45) * mrSges(7,3) + t124 * mrSges(6,3) + (t208 * t118 + (t121 * t159 + t178) * qJD(4) + m(7) * (pkin(5) * t144 + t106 * t165) - t165 * t197) * t114; m(7) * (-pkin(5) * t140 - t24 * t53 - t25 * t52 + t27 * t70 - t28 * t69) + (-t24 * t78 + t25 * t79 + t44 * t69 - t45 * t70) * mrSges(7,3) + (m(6) * (t170 * t189 - t190) + (m(7) * t106 + t159) * t118) * qJD(4) + ((mrSges(6,3) * t170 - mrSges(5,2)) * qJD(4) - t208) * t121; (t106 * t156 + t27 * t53 + t28 * t52) * t199 + t44 * t48 + t79 * t20 - t45 * t47 + t78 * t19 + 0.2e1 * t46 * t156 + 0.2e1 * t106 * t18 - t95 * t164 + t117 * t83 - 0.2e1 * pkin(4) * t80 + (qJD(5) * t96 + t82) * t120 + 0.2e1 * (t27 * t78 - t28 * t79 - t44 * t52 - t45 * t53) * mrSges(7,3); t15 * mrSges(6,1) - t14 * mrSges(6,2) + (t118 * t205 - t157) * qJD(4) + (m(7) * (t10 * t160 + t116 * t2 + t119 * t3 - t161 * t9) + t54 * t160 + t116 * t17 - t55 * t161 + t119 * t16) * pkin(5) + t153 + t187 + t204; t50 * mrSges(6,1) - t49 * mrSges(6,2) + (t116 * t12 + t119 * t13 + (-t116 * t39 + t119 * t40) * qJD(6)) * t196 + t139; (-t116 * t24 - t119 * t25 + (t116 * t69 + t119 * t70) * qJD(6)) * t196 + t207; t107 + (pkin(8) * t94 - t181) * qJD(5) + (m(7) * (t116 * t27 + t119 * t28 + (-t116 * t52 + t119 * t53) * qJD(6)) + (-t116 * t45 - t119 * t44 + (t116 * t79 + t119 * t78) * qJD(6)) * mrSges(7,3)) * pkin(5) + t128; 0.2e1 * t75; t127 + t204; t139; t203; t128; t75; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
