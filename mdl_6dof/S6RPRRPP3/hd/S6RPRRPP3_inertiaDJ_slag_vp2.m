% Calculate time derivative of joint inertia matrix for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
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
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:10
% EndTime: 2019-03-09 04:35:16
% DurationCPUTime: 3.06s
% Computational Cost: add. (1561->417), mult. (3703->554), div. (0->0), fcn. (2459->6), ass. (0->165)
t209 = -Ifges(7,4) - Ifges(6,5);
t193 = Ifges(5,5) + Ifges(7,5);
t118 = sin(qJ(4));
t120 = cos(qJ(4));
t204 = t118 ^ 2 + t120 ^ 2;
t121 = cos(qJ(3));
t167 = qJD(3) * t121;
t149 = t120 * t167;
t119 = sin(qJ(3));
t165 = qJD(4) * t119;
t151 = t118 * t165;
t124 = t149 - t151;
t164 = qJD(4) * t120;
t123 = t118 * t167 + t119 * t164;
t208 = -Ifges(6,1) - Ifges(7,1) - Ifges(5,3);
t154 = -cos(pkin(9)) * pkin(1) - pkin(2);
t207 = 0.2e1 * t154;
t206 = mrSges(6,3) + mrSges(7,2);
t175 = qJ(5) * t118;
t128 = -pkin(4) * t120 - t175;
t72 = -pkin(3) + t128;
t73 = t120 * mrSges(6,2) - t118 * mrSges(6,3);
t205 = m(6) * t72 + t73;
t74 = -t120 * mrSges(5,1) + t118 * mrSges(5,2);
t203 = -m(5) * pkin(3) - mrSges(4,1) + t74;
t117 = -pkin(4) - qJ(6);
t202 = m(7) * t117 - mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t201 = 0.2e1 * m(5);
t200 = 0.2e1 * m(7);
t199 = 2 * mrSges(7,1);
t102 = sin(pkin(9)) * pkin(1) + pkin(7);
t198 = 0.2e1 * t102;
t196 = pkin(5) + pkin(8);
t194 = mrSges(6,1) + mrSges(7,1);
t168 = qJD(3) * t119;
t31 = mrSges(6,1) * t123 - mrSges(6,3) * t168;
t32 = -mrSges(7,1) * t123 + mrSges(7,2) * t168;
t192 = -t31 + t32;
t28 = mrSges(5,1) * t168 - mrSges(5,3) * t124;
t176 = mrSges(6,1) * t149 + mrSges(6,2) * t168;
t33 = -mrSges(6,1) * t151 + t176;
t191 = t33 - t28;
t49 = -pkin(3) * t121 - t119 * pkin(8) + t154;
t70 = (pkin(3) * t119 - pkin(8) * t121) * qJD(3);
t190 = t118 * t70 + t49 * t164;
t171 = t119 * t120;
t61 = -mrSges(5,1) * t121 - mrSges(5,3) * t171;
t65 = mrSges(6,1) * t171 - mrSges(6,2) * t121;
t189 = -t61 + t65;
t173 = t118 * t119;
t63 = mrSges(6,1) * t173 + mrSges(6,3) * t121;
t64 = -mrSges(7,1) * t173 - mrSges(7,2) * t121;
t188 = -t63 + t64;
t170 = t120 * t121;
t67 = t102 * t170;
t24 = t118 * t49 + t67;
t185 = Ifges(5,4) * t118;
t184 = Ifges(5,4) * t120;
t183 = Ifges(6,6) * t118;
t182 = Ifges(6,6) * t120;
t181 = Ifges(7,6) * t118;
t180 = Ifges(7,6) * t120;
t179 = t121 * Ifges(6,4);
t178 = t121 * Ifges(5,6);
t177 = pkin(4) * t173 + t119 * t102;
t174 = qJD(3) * mrSges(7,3);
t172 = t118 * t121;
t169 = qJ(5) * qJD(5);
t166 = qJD(4) * t118;
t163 = qJD(6) * t119;
t162 = mrSges(5,2) - t206;
t138 = Ifges(5,1) * t120 - t185;
t37 = -Ifges(5,5) * t121 + t119 * t138;
t129 = Ifges(7,3) * t120 + t181;
t38 = -Ifges(7,5) * t121 + t119 * t129;
t136 = -Ifges(6,2) * t120 + t183;
t41 = t119 * t136 - t179;
t160 = t41 - t37 - t38;
t60 = mrSges(5,2) * t121 - mrSges(5,3) * t173;
t159 = t60 + t188;
t62 = mrSges(7,1) * t171 + mrSges(7,3) * t121;
t158 = t62 + t189;
t157 = qJ(5) * t124 + qJD(5) * t171;
t84 = t196 * t120;
t155 = mrSges(7,1) * t166;
t66 = t102 * t172;
t153 = t102 * t168;
t148 = -t102 * t118 - pkin(4);
t23 = t120 * t49 - t66;
t147 = t102 * t120 - qJ(5);
t146 = pkin(4) * t166 - qJD(5) * t118;
t145 = qJD(4) * t67 - t120 * t70 + t49 * t166;
t130 = -Ifges(7,3) * t118 + t180;
t135 = Ifges(6,2) * t118 + t182;
t81 = Ifges(5,1) * t118 + t184;
t144 = -t130 / 0.2e1 + t135 / 0.2e1 + t81 / 0.2e1;
t131 = Ifges(6,3) * t120 + t183;
t134 = Ifges(7,2) * t120 - t181;
t80 = Ifges(5,2) * t120 + t185;
t143 = -t131 / 0.2e1 - t134 / 0.2e1 - t80 / 0.2e1;
t17 = qJ(5) * t121 - t24;
t142 = -pkin(4) * t123 + t157;
t141 = mrSges(5,1) * t118 + mrSges(5,2) * t120;
t140 = -mrSges(6,2) * t118 - mrSges(6,3) * t120;
t139 = -mrSges(7,2) * t120 + mrSges(7,3) * t118;
t137 = -Ifges(5,2) * t118 + t184;
t133 = Ifges(7,2) * t118 + t180;
t132 = Ifges(6,3) * t118 - t182;
t36 = t119 * t137 - t178;
t39 = -Ifges(6,5) * t121 + t119 * t132;
t40 = -Ifges(7,4) * t121 + t119 * t133;
t127 = -t36 + t39 + t40 + t178;
t126 = -qJ(5) * t120 + qJ(6) * t118;
t125 = t120 * (m(6) * pkin(8) + t194);
t8 = t102 * t167 - t142;
t122 = -Ifges(6,4) * t151 + t123 * t209 - t193 * t149 + t208 * t168;
t113 = t121 * pkin(4);
t109 = Ifges(7,4) * t166;
t108 = Ifges(5,5) * t164;
t107 = Ifges(6,5) * t166;
t106 = Ifges(7,5) * t164;
t87 = mrSges(7,1) * t149;
t83 = t196 * t118;
t75 = -mrSges(7,2) * t118 - mrSges(7,3) * t120;
t69 = qJD(4) * t84;
t68 = t196 * t166;
t59 = t138 * qJD(4);
t58 = t137 * qJD(4);
t57 = t136 * qJD(4);
t56 = t133 * qJD(4);
t55 = t132 * qJD(4);
t54 = t129 * qJD(4);
t53 = t141 * qJD(4);
t52 = t140 * qJD(4);
t51 = t139 * qJD(4);
t48 = t141 * t119;
t47 = t140 * t119;
t46 = t139 * t119;
t45 = t117 * t120 - pkin(3) - t175;
t43 = -qJ(5) * t164 + t146;
t30 = t87 + (-t155 - t174) * t119;
t29 = -mrSges(5,2) * t168 - mrSges(5,3) * t123;
t27 = -qJ(5) * t171 + t177;
t26 = qJD(4) * t126 - qJD(6) * t120 + t146;
t22 = t119 * t126 + t177;
t21 = -mrSges(7,2) * t124 + mrSges(7,3) * t123;
t20 = mrSges(5,1) * t123 + mrSges(5,2) * t124;
t19 = -mrSges(6,2) * t123 - mrSges(6,3) * t124;
t18 = t113 - t23;
t16 = t135 * t165 + (Ifges(6,4) * t119 + t121 * t136) * qJD(3);
t15 = t134 * t165 + (Ifges(7,4) * t119 + t121 * t133) * qJD(3);
t14 = t131 * t165 + (Ifges(6,5) * t119 + t121 * t132) * qJD(3);
t13 = t130 * t165 + (Ifges(7,5) * t119 + t121 * t129) * qJD(3);
t12 = -t81 * t165 + (Ifges(5,5) * t119 + t121 * t138) * qJD(3);
t11 = -t80 * t165 + (Ifges(5,6) * t119 + t121 * t137) * qJD(3);
t10 = -pkin(5) * t173 - t17;
t9 = qJ(6) * t121 + t113 + t66 + (pkin(5) * t119 - t49) * t120;
t7 = t118 * t153 - t145;
t6 = (-t120 * t168 - t121 * t166) * t102 + t190;
t5 = t148 * t168 + t145;
t4 = (t102 * t166 + qJD(5)) * t121 + t147 * t168 - t190;
t3 = qJ(6) * t123 + t118 * t163 + t8;
t2 = -qJD(5) * t121 + (-pkin(5) * t171 - t66) * qJD(4) + (-pkin(5) * t172 - t119 * t147) * qJD(3) + t190;
t1 = -pkin(5) * t151 + qJD(6) * t121 + (pkin(5) * t170 + (-qJ(6) + t148) * t119) * qJD(3) + t145;
t25 = [0.2e1 * t6 * t60 + 0.2e1 * t7 * t61 + 0.2e1 * t1 * t62 + 0.2e1 * t4 * t63 + 0.2e1 * t2 * t64 + 0.2e1 * t5 * t65 + 0.2e1 * t3 * t46 + 0.2e1 * t8 * t47 + 0.2e1 * t27 * t19 + 0.2e1 * t23 * t28 + 0.2e1 * t24 * t29 + 0.2e1 * t9 * t30 + 0.2e1 * t17 * t31 + 0.2e1 * t10 * t32 + 0.2e1 * t18 * t33 + 0.2e1 * t22 * t21 + (t1 * t9 + t10 * t2 + t22 * t3) * t200 + 0.2e1 * m(6) * (t17 * t4 + t18 * t5 + t27 * t8) + (t23 * t7 + t24 * t6) * t201 + ((t48 * t198 + mrSges(4,2) * t207 + 0.2e1 * Ifges(4,4) * t121 + (-t160 + t179) * t120 + t127 * t118) * qJD(3) + t122) * t121 + (t20 * t198 + (t12 + t13 - t16) * t120 + (-t11 + t14 + t15) * t118 + (t127 * t120 + (t121 * t193 + t160) * t118) * qJD(4) + (mrSges(4,1) * t207 - 0.2e1 * Ifges(4,4) * t119 + (-Ifges(6,4) + t193) * t171 + (-Ifges(5,6) - t209) * t173 + (t102 ^ 2 * t201 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) + t208) * t121) * qJD(3)) * t119; (-t19 - t20 - t21 - m(6) * t8 - m(7) * t3 + (t159 * t120 + t158 * t118 + m(7) * (t10 * t120 + t118 * t9) + m(6) * (t118 * t18 - t120 * t17) + (-t102 * t121 - t118 * t23 + t120 * t24) * m(5)) * qJD(3)) * t121 + ((t29 + t192) * t120 + (t30 + t191) * t118 + (t46 + t47 + t48) * qJD(3) + (-t118 * t159 + t120 * t158) * qJD(4) + m(7) * (qJD(3) * t22 + t1 * t118 - t10 * t166 + t120 * t2 + t164 * t9) + m(6) * (qJD(3) * t27 + t118 * t5 - t120 * t4 + t164 * t18 + t166 * t17) + (-t118 * t7 + t120 * t6 - t164 * t23 - t166 * t24 + t153) * m(5)) * t119; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t204) * t119 * t167; -pkin(3) * t20 + t72 * t19 + t45 * t21 + t22 * t51 + t26 * t46 + t27 * t52 + t3 * t75 + t83 * t30 + t84 * t32 + t43 * t47 + t69 * t62 - t68 * t64 + t8 * t73 + m(7) * (t1 * t83 - t10 * t68 + t2 * t84 + t22 * t26 + t3 * t45 + t69 * t9) + m(6) * (t27 * t43 + t72 * t8) + (-t107 / 0.2e1 - t109 / 0.2e1 - t106 / 0.2e1 - t108 / 0.2e1 + (t203 * t102 + Ifges(4,5)) * qJD(3)) * t121 + (t11 / 0.2e1 - t14 / 0.2e1 - t15 / 0.2e1 + t6 * mrSges(5,3) + t2 * mrSges(7,1) - t4 * mrSges(6,1) + t144 * t167 + (t179 / 0.2e1 + t37 / 0.2e1 + t38 / 0.2e1 - t41 / 0.2e1 + t18 * mrSges(6,1) - t23 * mrSges(5,3) + t9 * mrSges(7,1)) * qJD(4) + (t29 - t31 + t189 * qJD(4) + m(5) * (-qJD(4) * t23 + t6) + m(6) * (qJD(4) * t18 - t4)) * pkin(8)) * t120 + (t5 * mrSges(6,1) - t7 * mrSges(5,3) + t1 * mrSges(7,1) + t12 / 0.2e1 + t13 / 0.2e1 - t16 / 0.2e1 + t143 * t167 + (t178 / 0.2e1 - t36 / 0.2e1 + t39 / 0.2e1 + t40 / 0.2e1 - t10 * mrSges(7,1) + t17 * mrSges(6,1) - t24 * mrSges(5,3)) * qJD(4) + ((-t60 + t63) * qJD(4) + m(5) * (-qJD(4) * t24 - t7) + m(6) * (qJD(4) * t17 + t5) + t191) * pkin(8)) * t118 + (t102 * t53 + (t54 / 0.2e1 - t57 / 0.2e1 + t59 / 0.2e1) * t120 + (t55 / 0.2e1 + t56 / 0.2e1 - t58 / 0.2e1) * t118 + (t102 * mrSges(4,2) - Ifges(4,6) + (-Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1 + Ifges(5,6) / 0.2e1) * t120 + (Ifges(7,5) / 0.2e1 - Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t118) * qJD(3) + (-t118 * t144 + t120 * t143) * qJD(4)) * t119; m(7) * (t118 * t69 - t120 * t68 + t164 * t83 - t166 * t84) * t119 + (-m(6) * t43 - m(7) * t26 - t51 - t52 - t53) * t121 + ((m(7) * t45 + t203 + t205 + t75) * t119 + (m(7) * (t118 * t83 + t120 * t84) - mrSges(4,2) + t204 * (mrSges(5,3) + t194)) * t121) * qJD(3) + (m(6) + m(5)) * t204 * pkin(8) * t167; (t26 * t45 - t68 * t84 + t69 * t83) * t200 - 0.2e1 * pkin(3) * t53 + 0.2e1 * t72 * t52 + 0.2e1 * t26 * t75 + 0.2e1 * t45 * t51 + 0.2e1 * t205 * t43 + (-t68 * t199 - t55 - t56 + t58) * t120 + (t69 * t199 + t54 - t57 + t59) * t118 + ((t199 * t83 - t130 + t135 + t81) * t120 + (-0.2e1 * mrSges(7,1) * t84 - t131 - t134 - t80) * t118) * qJD(4); -t122 + t188 * qJD(5) + t192 * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t10 - qJD(6) * t9 + t1 * t117) + m(6) * (-pkin(4) * t5 - qJ(5) * t4 - qJD(5) * t17) + t117 * t30 - qJD(6) * t62 - pkin(4) * t33 + t5 * mrSges(6,2) - t6 * mrSges(5,2) + t7 * mrSges(5,1) - t1 * mrSges(7,3) + t2 * mrSges(7,2) - t4 * mrSges(6,3) + (-Ifges(6,4) * t120 - Ifges(5,6) * t118) * t167 + (-Ifges(5,6) * t120 - t118 * t193) * t165; m(6) * t142 + m(7) * t157 + t202 * t165 * t120 + (-m(7) * t163 + t162 * t165) * t118 + (t202 * t118 - t162 * t120) * t167; t109 + t106 + m(7) * (-qJ(5) * t68 + qJD(5) * t84 - qJD(6) * t83 + t117 * t69) + t107 + t108 - t68 * mrSges(7,2) - t69 * mrSges(7,3) - qJD(6) * t118 * mrSges(7,1) + qJD(5) * t125 + ((-mrSges(6,1) * pkin(4) + mrSges(7,1) * t117 - Ifges(6,4)) * t120 + (-qJ(5) * t194 - Ifges(5,6)) * t118 + (m(6) * t128 + t73 + t74) * pkin(8)) * qJD(4); 0.2e1 * m(6) * t169 + 0.2e1 * qJD(6) * mrSges(7,3) + 0.2e1 * m(7) * (-qJD(6) * t117 + t169) + 0.2e1 * t206 * qJD(5); t87 + m(6) * t5 + m(7) * t1 + (-t166 * t194 - t174) * t119 + t176; (m(6) + m(7)) * t123; m(7) * t69 + qJD(4) * t125; -m(7) * qJD(6); 0; m(7) * t2 + t32; t124 * m(7); -m(7) * t68 - t155; m(7) * qJD(5); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
