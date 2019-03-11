% Calculate time derivative of joint inertia matrix for
% S6PRPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:15
% EndTime: 2019-03-08 20:45:20
% DurationCPUTime: 2.75s
% Computational Cost: add. (2766->402), mult. (6838->614), div. (0->0), fcn. (5979->10), ass. (0->175)
t114 = sin(qJ(6));
t115 = sin(qJ(5));
t118 = cos(qJ(6));
t119 = cos(qJ(5));
t128 = t114 * t115 - t118 * t119;
t196 = qJD(5) + qJD(6);
t199 = t196 * t128;
t83 = t114 * t119 + t115 * t118;
t46 = t196 * t83;
t20 = mrSges(7,1) * t46 - mrSges(7,2) * t199;
t135 = mrSges(6,1) * t115 + mrSges(6,2) * t119;
t84 = t135 * qJD(5);
t201 = t20 + t84;
t116 = sin(qJ(4));
t122 = -pkin(2) - pkin(8);
t167 = t116 * t122;
t104 = t119 * t167;
t120 = cos(qJ(4));
t180 = pkin(9) * t120;
t181 = pkin(4) * t116;
t92 = qJ(3) - t180 + t181;
t63 = t115 * t92 + t104;
t200 = qJD(5) * t63;
t198 = t116 * mrSges(5,1) + t120 * mrSges(5,2) + mrSges(4,3);
t162 = qJD(4) * t116;
t141 = t115 * t162;
t160 = qJD(4) * t120;
t197 = Ifges(6,6) * t141 + Ifges(6,3) * t160;
t164 = t115 ^ 2 + t119 ^ 2;
t106 = -pkin(5) * t119 - pkin(4);
t96 = -mrSges(6,1) * t119 + mrSges(6,2) * t115;
t179 = t96 - mrSges(5,1);
t48 = mrSges(7,1) * t128 + mrSges(7,2) * t83;
t195 = m(7) * t106 + t179 + t48;
t113 = cos(pkin(6));
t112 = sin(pkin(6));
t117 = sin(qJ(2));
t170 = t112 * t117;
t146 = qJD(2) * t170;
t121 = cos(qJ(2));
t169 = t112 * t121;
t149 = t116 * t169;
t52 = -qJD(4) * t149 + t113 * t160 - t120 * t146;
t172 = t120 * t52;
t73 = t113 * t116 + t120 * t169;
t171 = qJD(4) * t73;
t51 = t116 * t146 - t171;
t74 = t113 * t120 - t149;
t194 = (t116 * t73 + t120 * t74) * qJD(4) + t116 * t51 - t172;
t193 = 0.2e1 * m(7);
t192 = -0.2e1 * t122;
t191 = 0.2e1 * t122;
t190 = m(5) / 0.2e1;
t189 = m(6) / 0.2e1;
t188 = -m(7) / 0.2e1;
t187 = m(6) * pkin(4);
t186 = m(7) * pkin(5);
t185 = -t128 / 0.2e1;
t184 = t83 / 0.2e1;
t183 = -pkin(10) - pkin(9);
t182 = -t115 / 0.2e1;
t35 = t73 * t52;
t178 = Ifges(6,4) * t115;
t177 = Ifges(6,4) * t119;
t176 = Ifges(6,5) * t115;
t175 = Ifges(6,6) * t115;
t174 = Ifges(6,6) * t119;
t173 = t116 * Ifges(6,6);
t168 = t115 * t120;
t166 = t119 * t120;
t163 = qJD(2) * t121;
t145 = t112 * t163;
t165 = qJ(3) * t145 + qJD(3) * t170;
t161 = qJD(4) * t119;
t159 = qJD(4) * t122;
t158 = qJD(5) * t115;
t157 = qJD(5) * t119;
t156 = qJD(5) * t120;
t155 = qJD(6) * t114;
t154 = qJD(6) * t118;
t27 = -t46 * t120 + t128 * t162;
t29 = t120 * t199 + t83 * t162;
t152 = Ifges(7,5) * t27 + Ifges(7,6) * t29 + Ifges(7,3) * t160;
t151 = pkin(5) * t158;
t54 = t115 * t170 + t119 * t74;
t10 = -t54 * qJD(5) - t115 * t51 + t119 * t145;
t53 = -t115 * t74 + t119 * t170;
t11 = t53 * qJD(5) + t115 * t145 + t119 * t51;
t18 = -t114 * t54 + t118 * t53;
t3 = t18 * qJD(6) + t10 * t114 + t11 * t118;
t19 = t114 * t53 + t118 * t54;
t4 = -t19 * qJD(6) + t10 * t118 - t11 * t114;
t150 = t4 * mrSges(7,1) - t3 * mrSges(7,2);
t64 = t73 * t162;
t68 = t83 * t120;
t70 = t128 * t120;
t40 = mrSges(7,1) * t68 - mrSges(7,2) * t70;
t81 = (pkin(5) * t115 - t122) * t120;
t148 = m(7) * t81 + t40;
t147 = qJD(5) * t183;
t144 = t120 * t159;
t143 = t115 * t156;
t142 = t119 * t156;
t140 = t116 * t158;
t26 = -qJD(4) * t70 - t46 * t116;
t28 = -qJD(4) * t68 + t116 * t199;
t139 = t28 * mrSges(7,1) - t26 * mrSges(7,2);
t138 = -Ifges(6,5) * t119 + (2 * Ifges(5,4));
t137 = -t115 * t122 + pkin(5);
t80 = qJD(3) + (pkin(4) * t120 + pkin(9) * t116) * qJD(4);
t33 = t115 * t80 + t119 * t144 - t122 * t140 + t92 * t157;
t79 = t119 * t92;
t62 = -t115 * t167 + t79;
t136 = -qJD(5) * t62 + t33;
t134 = Ifges(6,1) * t119 - t178;
t99 = Ifges(6,1) * t115 + t177;
t133 = -Ifges(6,2) * t115 + t177;
t98 = Ifges(6,2) * t119 + t178;
t132 = -t10 * t115 + t11 * t119;
t42 = -pkin(10) * t166 + t137 * t116 + t79;
t47 = -pkin(10) * t168 + t63;
t13 = -t114 * t47 + t118 * t42;
t14 = t114 * t42 + t118 * t47;
t100 = t183 * t115;
t101 = t183 * t119;
t55 = t100 * t118 + t101 * t114;
t90 = t115 * t147;
t91 = t119 * t147;
t31 = t55 * qJD(6) + t114 * t91 + t118 * t90;
t56 = t100 * t114 - t101 * t118;
t32 = -t56 * qJD(6) - t114 * t90 + t118 * t91;
t43 = Ifges(7,6) * t46;
t44 = Ifges(7,5) * t199;
t129 = t32 * mrSges(7,1) - t31 * mrSges(7,2) - t43 - t44;
t72 = t119 * t80;
t12 = t72 + (-t104 + (pkin(10) * t120 - t92) * t115) * qJD(5) + (pkin(10) * t116 * t119 + t137 * t120) * qJD(4);
t125 = t141 - t142;
t17 = t125 * pkin(10) + t33;
t5 = t13 * qJD(6) + t114 * t12 + t118 * t17;
t6 = -t14 * qJD(6) - t114 * t17 + t118 * t12;
t127 = t6 * mrSges(7,1) - t5 * mrSges(7,2) + t152;
t126 = t116 * t161 + t143;
t109 = Ifges(6,5) * t157;
t89 = mrSges(6,1) * t116 - mrSges(6,3) * t166;
t88 = -mrSges(6,2) * t116 - mrSges(6,3) * t168;
t87 = t134 * qJD(5);
t86 = t133 * qJD(5);
t85 = (mrSges(5,1) * t120 - mrSges(5,2) * t116) * qJD(4);
t77 = (-mrSges(7,1) * t114 - mrSges(7,2) * t118) * qJD(6) * pkin(5);
t76 = t135 * t120;
t69 = t128 * t116;
t67 = t83 * t116;
t66 = Ifges(6,5) * t116 + t134 * t120;
t65 = t133 * t120 + t173;
t61 = -mrSges(6,2) * t160 + t125 * mrSges(6,3);
t60 = mrSges(6,1) * t160 + t126 * mrSges(6,3);
t59 = -t125 * pkin(5) + t116 * t159;
t58 = mrSges(7,1) * t116 + t70 * mrSges(7,3);
t57 = -mrSges(7,2) * t116 - t68 * mrSges(7,3);
t50 = Ifges(7,1) * t83 - Ifges(7,4) * t128;
t49 = Ifges(7,4) * t83 - Ifges(7,2) * t128;
t41 = -t125 * mrSges(6,1) - t126 * mrSges(6,2);
t39 = -t99 * t156 + (Ifges(6,5) * t120 - t134 * t116) * qJD(4);
t38 = -t98 * t156 + (Ifges(6,6) * t120 - t133 * t116) * qJD(4);
t37 = -Ifges(7,1) * t70 - Ifges(7,4) * t68 + Ifges(7,5) * t116;
t36 = -Ifges(7,4) * t70 - Ifges(7,2) * t68 + Ifges(7,6) * t116;
t34 = -t115 * t144 - t200 + t72;
t22 = -Ifges(7,1) * t199 - Ifges(7,4) * t46;
t21 = -Ifges(7,4) * t199 - Ifges(7,2) * t46;
t16 = -mrSges(7,2) * t160 + t29 * mrSges(7,3);
t15 = mrSges(7,1) * t160 - t27 * mrSges(7,3);
t9 = -mrSges(7,1) * t29 + mrSges(7,2) * t27;
t8 = Ifges(7,1) * t27 + Ifges(7,4) * t29 + Ifges(7,5) * t160;
t7 = Ifges(7,4) * t27 + Ifges(7,2) * t29 + Ifges(7,6) * t160;
t1 = [0.2e1 * m(6) * (t10 * t53 + t11 * t54 + t35) + 0.2e1 * m(7) * (t18 * t4 + t19 * t3 + t35) + 0.2e1 * m(5) * (t112 ^ 2 * t117 * t163 + t74 * t51 + t35); t10 * t89 + t11 * t88 + t18 * t15 + t19 * t16 + t3 * t57 + t4 * t58 + t53 * t60 + t54 * t61 + (t41 + t9) * t73 + (t40 + t76) * t52 - t194 * mrSges(5,3) + (t117 * t85 + ((-mrSges(3,1) + mrSges(4,2)) * t117 + (-mrSges(3,2) + t198) * t121) * qJD(2)) * t112 + m(7) * (t13 * t4 + t14 * t3 + t6 * t18 + t5 * t19 + t52 * t81 + t59 * t73) + m(4) * (-pkin(2) * t146 + t165) + m(6) * (t62 * t10 + t63 * t11 + t33 * t54 + t34 * t53) + m(5) * t165 + ((t64 - t172) * t189 + t194 * t190) * t191; 0.2e1 * qJ(3) * t85 + 0.2e1 * t13 * t15 + 0.2e1 * t14 * t16 + t27 * t37 + t29 * t36 + 0.2e1 * t33 * t88 + 0.2e1 * t34 * t89 + 0.2e1 * t59 * t40 + 0.2e1 * t5 * t57 + 0.2e1 * t6 * t58 + 0.2e1 * t62 * t60 + 0.2e1 * t63 * t61 - t68 * t7 - t70 * t8 + 0.2e1 * t81 * t9 + 0.2e1 * m(6) * (t63 * t33 + t62 * t34) + (t13 * t6 + t14 * t5 + t59 * t81) * t193 + 0.2e1 * ((m(4) + m(5)) * qJ(3) + t198) * qJD(3) + ((t115 * t65 + t138 * t116 - t119 * t66 + t76 * t191) * qJD(4) + t152 + t197) * t116 + (-t115 * t38 + t119 * t39 + t41 * t192 + (-t119 * t65 - t115 * t66 + t116 * (-t174 - t176)) * qJD(5) + (-Ifges(7,5) * t70 - Ifges(7,6) * t68 + (-t138 - t175) * t120 + (-0.2e1 * m(6) * t122 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t116) * qJD(4)) * t120; m(4) * t146 + m(6) * t64 + m(7) * (t28 * t18 + t26 * t19 - t69 * t3 - t67 * t4 + t64) + 0.2e1 * ((-qJD(4) * t115 * t53 + t54 * t161 - t52) * t189 + t52 * t188 + (qJD(4) * t74 - t52) * t190) * t120 + 0.2e1 * ((-t53 * t157 - t54 * t158 + t132) * t189 + (t51 + t171) * t190) * t116; m(7) * (-t120 * t59 + t28 * t13 + t26 * t14 - t69 * t5 - t67 * t6) + t26 * t57 - t69 * t16 + t28 * t58 - t67 * t15 - t120 * t9 - t120 * t41 + (-t115 * t89 + m(6) * (-t115 * t62 + t119 * t63) + t119 * t88) * t160 + (-t89 * t157 - t115 * t60 + m(6) * (-t34 * t115 + t33 * t119 - t62 * t157 - t63 * t158) - t88 * t158 + t119 * t61 + (m(6) * t120 * t192 + t148 + t76) * qJD(4)) * t116; (-t69 * t26 - t67 * t28) * t193 + 0.4e1 * ((-0.1e1 + t164) * t189 + t188) * t116 * t160; -t51 * mrSges(5,2) + t201 * t73 + m(7) * (t73 * t151 + t32 * t18 + t31 * t19 + t56 * t3 + t55 * t4) + (-t128 * t3 + t18 * t199 - t19 * t46 - t4 * t83) * mrSges(7,3) + (-t187 + t195) * t52 + (m(6) * pkin(9) + mrSges(6,3)) * ((-t115 * t54 - t119 * t53) * qJD(5) + t132); -pkin(4) * t41 - t199 * t37 / 0.2e1 - t46 * t36 / 0.2e1 + t29 * t49 / 0.2e1 + t27 * t50 / 0.2e1 + t55 * t15 + t56 * t16 + t31 * t57 + t32 * t58 + t59 * t48 - t68 * t21 / 0.2e1 - t70 * t22 / 0.2e1 + t81 * t20 + t7 * t185 + t8 * t184 + m(7) * (t106 * t59 + t32 * t13 + t31 * t14 + t56 * t5 + t55 * t6) + t106 * t9 + (-t44 / 0.2e1 - t43 / 0.2e1 + t109 / 0.2e1 + (-Ifges(5,5) + (t179 - t187) * t122) * qJD(4)) * t116 + (t98 * t162 / 0.2e1 - t34 * mrSges(6,3) + t39 / 0.2e1 + (-t63 * mrSges(6,3) - t65 / 0.2e1 - t173 / 0.2e1 + t148 * pkin(5)) * qJD(5) + (-t60 + m(6) * (-t34 - t200) - qJD(5) * t88) * pkin(9)) * t115 + (-t128 * t5 + t13 * t199 - t14 * t46 - t6 * t83) * mrSges(7,3) + (-t99 * t162 / 0.2e1 + qJD(5) * t66 / 0.2e1 + t38 / 0.2e1 + t136 * mrSges(6,3) + (m(6) * t136 - qJD(5) * t89 + t61) * pkin(9)) * t119 + (t86 * t182 + t119 * t87 / 0.2e1 - t122 * t84 + (t99 * t182 - t119 * t98 / 0.2e1) * qJD(5) + (-t122 * mrSges(5,2) + Ifges(7,5) * t184 + Ifges(7,6) * t185 - Ifges(5,6) + t176 / 0.2e1 + t174 / 0.2e1) * qJD(4)) * t120; m(7) * (-pkin(5) * t143 + t56 * t26 + t55 * t28 - t31 * t69 - t32 * t67) + (-t128 * t26 - t199 * t67 - t28 * t83 + t46 * t69) * mrSges(7,3) + (m(6) * (t164 * t180 - t181) + t195 * t116) * qJD(4) + ((t164 * mrSges(6,3) - mrSges(5,2)) * qJD(4) - t201) * t120; 0.2e1 * t48 * t151 + 0.2e1 * t106 * t20 + (t106 * t151 + t56 * t31 + t55 * t32) * t193 - t199 * t50 + t83 * t22 - t46 * t49 - t128 * t21 - 0.2e1 * pkin(4) * t84 + t115 * t87 - t98 * t158 + (qJD(5) * t99 + t86) * t119 + 0.2e1 * (-t128 * t31 + t199 * t55 - t32 * t83 - t46 * t56) * mrSges(7,3); t10 * mrSges(6,1) - t11 * mrSges(6,2) + (t114 * t3 + t118 * t4 + (-t114 * t18 + t118 * t19) * qJD(6)) * t186 + t150; -Ifges(6,6) * t142 + t34 * mrSges(6,1) - t33 * mrSges(6,2) - t126 * Ifges(6,5) + (m(7) * (t114 * t5 + t118 * t6 - t13 * t155 + t14 * t154) + t57 * t154 + t114 * t16 - t58 * t155 + t118 * t15) * pkin(5) + t127 + t197; (-t119 * t160 + t140) * mrSges(6,2) + (-t115 * t160 - t116 * t157) * mrSges(6,1) + (t114 * t26 + t118 * t28 + (t114 * t67 - t118 * t69) * qJD(6)) * t186 + t139; t109 + (t96 * pkin(9) - t175) * qJD(5) + (m(7) * (t114 * t31 + t118 * t32 + (-t114 * t55 + t118 * t56) * qJD(6)) + (-t114 * t46 + t118 * t199 + (t114 * t83 - t118 * t128) * qJD(6)) * mrSges(7,3)) * pkin(5) + t129; 0.2e1 * t77; t150; t127; t139; t129; t77; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
