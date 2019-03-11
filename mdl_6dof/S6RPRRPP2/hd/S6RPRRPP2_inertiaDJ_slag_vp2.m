% Calculate time derivative of joint inertia matrix for
% S6RPRRPP2
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
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:30
% EndTime: 2019-03-09 04:31:34
% DurationCPUTime: 2.62s
% Computational Cost: add. (1553->418), mult. (3696->546), div. (0->0), fcn. (2462->6), ass. (0->174)
t204 = -Ifges(6,4) - Ifges(5,5);
t203 = -Ifges(6,2) - Ifges(5,3);
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t197 = t106 ^ 2 + t108 ^ 2;
t107 = sin(qJ(3));
t109 = cos(qJ(3));
t157 = qJD(3) * t109;
t135 = t106 * t157;
t154 = qJD(4) * t108;
t111 = t107 * t154 + t135;
t163 = qJ(5) * t106;
t188 = pkin(4) + pkin(5);
t195 = -t108 * t188 - t163;
t139 = -cos(pkin(9)) * pkin(1) - pkin(2);
t202 = 0.2e1 * t139;
t201 = m(6) + m(7);
t155 = qJD(4) * t107;
t137 = t106 * t155;
t134 = t108 * t157;
t84 = mrSges(7,2) * t134;
t20 = -t111 * mrSges(7,1) - mrSges(7,2) * t137 + t84;
t141 = t188 * t106;
t160 = t107 * t108;
t173 = qJ(5) * t134 + qJD(5) * t160;
t94 = sin(pkin(9)) * pkin(1) + pkin(7);
t4 = t195 * t155 + (-t94 - t141) * t157 + t173;
t200 = m(7) * t4 + t20;
t162 = qJ(5) * t108;
t114 = -t141 + t162;
t153 = qJD(5) * t106;
t33 = qJD(4) * t114 + t153;
t156 = qJD(4) * t106;
t54 = -mrSges(7,1) * t156 + mrSges(7,2) * t154;
t199 = m(7) * t33 + t54;
t119 = -pkin(4) * t108 - t163;
t71 = -pkin(3) + t119;
t73 = -t108 * mrSges(6,1) - t106 * mrSges(6,3);
t198 = m(6) * t71 + t73;
t75 = -t108 * mrSges(5,1) + t106 * mrSges(5,2);
t196 = -m(5) * pkin(3) - mrSges(4,1) + t75;
t194 = 0.2e1 * m(5);
t193 = 0.2e1 * m(7);
t192 = -2 * mrSges(7,3);
t191 = 0.2e1 * t94;
t185 = -mrSges(6,1) - mrSges(7,1);
t184 = mrSges(6,2) - mrSges(7,3);
t183 = -mrSges(7,2) - mrSges(6,3);
t182 = Ifges(5,6) + Ifges(7,6);
t181 = pkin(8) - qJ(6);
t112 = t134 - t137;
t158 = qJD(3) * t107;
t28 = mrSges(5,1) * t158 - mrSges(5,3) * t112;
t85 = mrSges(6,2) * t134;
t29 = t85 + (-mrSges(6,1) * qJD(3) - mrSges(6,2) * t156) * t107;
t180 = -t28 + t29;
t31 = -mrSges(5,2) * t158 - mrSges(5,3) * t111;
t32 = -mrSges(6,2) * t111 + mrSges(6,3) * t158;
t179 = t31 + t32;
t50 = -pkin(3) * t109 - pkin(8) * t107 + t139;
t159 = t108 * t109;
t69 = t94 * t159;
t178 = qJD(4) * t69 + t50 * t156;
t70 = (pkin(3) * t107 - pkin(8) * t109) * qJD(3);
t177 = t106 * t70 + t50 * t154;
t161 = t106 * t107;
t63 = mrSges(5,2) * t109 - mrSges(5,3) * t161;
t67 = -mrSges(6,2) * t161 - mrSges(6,3) * t109;
t176 = -t63 - t67;
t65 = -mrSges(5,1) * t109 - mrSges(5,3) * t160;
t66 = mrSges(6,1) * t109 + mrSges(6,2) * t160;
t175 = -t65 + t66;
t24 = t106 * t50 + t69;
t171 = Ifges(5,4) * t106;
t170 = Ifges(5,4) * t108;
t169 = Ifges(7,4) * t106;
t168 = Ifges(7,4) * t108;
t167 = Ifges(6,5) * t106;
t166 = Ifges(6,5) * t108;
t165 = t109 * Ifges(7,5);
t164 = t94 * t109;
t152 = qJD(6) * t108;
t151 = -mrSges(5,1) + t185;
t149 = -Ifges(7,5) - t204;
t148 = qJ(5) * t158 + t177;
t123 = Ifges(7,1) * t108 + t169;
t39 = t107 * t123 + t165;
t124 = Ifges(6,1) * t108 + t167;
t40 = -Ifges(6,4) * t109 + t107 * t124;
t125 = Ifges(5,1) * t108 - t171;
t41 = -Ifges(5,5) * t109 + t107 * t125;
t147 = -t39 - t40 - t41;
t62 = -mrSges(7,2) * t109 + mrSges(7,3) * t161;
t146 = t62 - t176;
t64 = mrSges(7,1) * t109 - mrSges(7,3) * t160;
t145 = t64 + t175;
t143 = mrSges(7,3) * t159;
t142 = Ifges(5,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t140 = t94 * t158;
t138 = t108 * t158;
t133 = -t106 * t94 - pkin(4);
t76 = t181 * t108;
t68 = t106 * t164;
t23 = t108 * t50 - t68;
t132 = -t108 * t70 + t178;
t77 = -Ifges(6,3) * t108 + t167;
t78 = -Ifges(7,2) * t108 + t169;
t79 = Ifges(5,2) * t108 + t171;
t131 = t77 / 0.2e1 + t78 / 0.2e1 - t79 / 0.2e1;
t80 = Ifges(7,1) * t106 - t168;
t81 = Ifges(6,1) * t106 - t166;
t82 = Ifges(5,1) * t106 + t170;
t130 = t80 / 0.2e1 + t81 / 0.2e1 + t82 / 0.2e1;
t17 = -qJ(5) * t109 + t24;
t3 = -t94 * t138 + (-t156 * t94 - qJD(5)) * t109 + t148;
t5 = t133 * t158 + t132;
t129 = t106 * t5 + t108 * t3;
t6 = (-t109 * t156 - t138) * t94 + t177;
t7 = t106 * t140 - t132;
t128 = -t106 * t7 + t108 * t6;
t127 = mrSges(5,1) * t106 + mrSges(5,2) * t108;
t126 = mrSges(6,1) * t106 - mrSges(6,3) * t108;
t122 = -Ifges(5,2) * t106 + t170;
t121 = Ifges(7,2) * t106 + t168;
t120 = Ifges(6,3) * t106 + t166;
t118 = pkin(4) * t106 - t162;
t117 = -t111 * Ifges(6,6) + t204 * t134 + t203 * t158;
t116 = t108 * (m(6) * pkin(8) + t184);
t115 = m(6) * t119;
t36 = -Ifges(6,6) * t109 + t107 * t120;
t37 = t109 * Ifges(7,6) + t107 * t121;
t38 = -t109 * Ifges(5,6) + t107 * t122;
t113 = t182 * t109 + t36 + t37 - t38;
t102 = t109 * pkin(4);
t101 = Ifges(6,4) * t154;
t100 = Ifges(5,5) * t154;
t98 = Ifges(6,6) * t156;
t86 = mrSges(7,3) * t137;
t74 = t108 * mrSges(7,1) + t106 * mrSges(7,2);
t72 = t181 * t106;
t61 = t125 * qJD(4);
t60 = t124 * qJD(4);
t59 = t123 * qJD(4);
t58 = t122 * qJD(4);
t57 = t121 * qJD(4);
t56 = t120 * qJD(4);
t55 = t127 * qJD(4);
t53 = t126 * qJD(4);
t51 = pkin(3) - t195;
t49 = t127 * t107;
t48 = (-mrSges(7,1) * t106 + mrSges(7,2) * t108) * t107;
t47 = t126 * t107;
t45 = qJD(4) * t76 - qJD(6) * t106;
t44 = qJD(4) * t118 - t153;
t43 = -t156 * t181 - t152;
t30 = mrSges(7,2) * t158 + mrSges(7,3) * t111;
t27 = t86 + (-mrSges(7,1) * t107 - t143) * qJD(3);
t26 = (t118 + t94) * t107;
t22 = (-t94 + t114) * t107;
t21 = mrSges(5,1) * t111 + mrSges(5,2) * t112;
t19 = mrSges(6,1) * t111 - mrSges(6,3) * t112;
t18 = t102 - t23;
t16 = -t82 * t155 + (Ifges(5,5) * t107 + t109 * t125) * qJD(3);
t15 = -t81 * t155 + (Ifges(6,4) * t107 + t109 * t124) * qJD(3);
t14 = -t80 * t155 + (-Ifges(7,5) * t107 + t109 * t123) * qJD(3);
t13 = -t79 * t155 + (Ifges(5,6) * t107 + t109 * t122) * qJD(3);
t12 = -t78 * t155 + (-Ifges(7,6) * t107 + t109 * t121) * qJD(3);
t11 = -t77 * t155 + (Ifges(6,6) * t107 + t109 * t120) * qJD(3);
t10 = qJ(6) * t161 + t17;
t9 = pkin(5) * t109 + t102 + t68 + (-qJ(6) * t107 - t50) * t108;
t8 = pkin(4) * t111 + qJ(5) * t137 + t157 * t94 - t173;
t2 = -qJD(5) * t109 + (qJ(6) * qJD(4) - qJD(3) * t94) * t160 + (qJD(6) * t107 + (qJ(6) * qJD(3) - qJD(4) * t94) * t109) * t106 + t148;
t1 = (-qJ(6) * t157 - t70) * t108 + (qJ(6) * t156 - t152 + (-pkin(5) + t133) * qJD(3)) * t107 + t178;
t25 = [0.2e1 * t2 * t62 + 0.2e1 * t6 * t63 + 0.2e1 * t1 * t64 + 0.2e1 * t7 * t65 + 0.2e1 * t5 * t66 + 0.2e1 * t3 * t67 + 0.2e1 * t8 * t47 + 0.2e1 * t4 * t48 + 0.2e1 * t26 * t19 + 0.2e1 * t9 * t27 + 0.2e1 * t23 * t28 + 0.2e1 * t18 * t29 + 0.2e1 * t10 * t30 + 0.2e1 * t24 * t31 + 0.2e1 * t17 * t32 + 0.2e1 * t22 * t20 + (t1 * t9 + t10 * t2 + t22 * t4) * t193 + 0.2e1 * m(6) * (t17 * t3 + t18 * t5 + t26 * t8) + (t23 * t7 + t24 * t6) * t194 + ((t49 * t191 + mrSges(4,2) * t202 + 0.2e1 * Ifges(4,4) * t109 + (-t147 + t165) * t108 + t113 * t106) * qJD(3) + t117) * t109 + (t21 * t191 + (t14 + t15 + t16) * t108 + (t11 + t12 - t13) * t106 + (t113 * t108 + (t109 * t149 + t147) * t106) * qJD(4) + (mrSges(4,1) * t202 - 0.2e1 * Ifges(4,4) * t107 + t149 * t160 + (Ifges(6,6) - t182) * t161 + (t94 ^ 2 * t194 + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - (2 * Ifges(7,3)) + t203) * t109) * qJD(3)) * t107; (-t19 - t21 - m(6) * t8 + (t146 * t108 + t145 * t106 + m(7) * (t10 * t108 + t106 * t9) + m(6) * (t106 * t18 + t108 * t17) + (-t106 * t23 + t108 * t24 - t164) * m(5)) * qJD(3) + t200) * t109 + ((t30 + t179) * t108 + (t27 + t180) * t106 + (t47 - t48 + t49) * qJD(3) + (-t106 * t146 + t108 * t145) * qJD(4) + m(7) * (-qJD(3) * t22 + t1 * t106 - t10 * t156 + t108 * t2 + t154 * t9) + m(6) * (qJD(3) * t26 + t154 * t18 - t156 * t17 + t129) + (-t154 * t23 - t156 * t24 + t128 + t140) * m(5)) * t107; 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (-0.1e1 + t197) * t107 * t157; -pkin(3) * t21 + t71 * t19 + t51 * t20 + t22 * t54 + t26 * t53 + t72 * t27 + t76 * t30 + t33 * t48 + t4 * t74 + t43 * t62 + t44 * t47 + t45 * t64 + t8 * t73 + m(7) * (t1 * t72 + t10 * t43 + t2 * t76 + t22 * t33 + t4 * t51 + t45 * t9) + m(6) * (t26 * t44 + t71 * t8) + (t6 * mrSges(5,3) - t2 * mrSges(7,3) + t3 * mrSges(6,2) - t11 / 0.2e1 - t12 / 0.2e1 + t13 / 0.2e1) * t108 + (t14 / 0.2e1 + t15 / 0.2e1 + t16 / 0.2e1 - t7 * mrSges(5,3) - t1 * mrSges(7,3) + t5 * mrSges(6,2)) * t106 + (m(5) * t128 + m(6) * t129 + t106 * t180 + t108 * t179) * pkin(8) + (-t101 / 0.2e1 - t98 / 0.2e1 - t100 / 0.2e1 + (t131 * t106 + t130 * t108 + t196 * t94 + Ifges(4,5)) * qJD(3)) * t109 + (t94 * t55 + (t59 / 0.2e1 + t60 / 0.2e1 + t61 / 0.2e1) * t108 + (t56 / 0.2e1 + t57 / 0.2e1 - t58 / 0.2e1) * t106 + (t94 * mrSges(4,2) - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 + t142) * t108 + (-Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t106) * qJD(3)) * t107 + ((t39 / 0.2e1 + t40 / 0.2e1 + t41 / 0.2e1 + t18 * mrSges(6,2) - t23 * mrSges(5,3) - t9 * mrSges(7,3) + t165 / 0.2e1 + t131 * t107) * t108 + (t36 / 0.2e1 + t37 / 0.2e1 - t38 / 0.2e1 - t17 * mrSges(6,2) - t24 * mrSges(5,3) + t10 * mrSges(7,3) + t142 * t109 - t130 * t107) * t106 + (t175 * t108 + t176 * t106 + m(5) * (-t106 * t24 - t108 * t23) + m(6) * (-t106 * t17 + t108 * t18)) * pkin(8)) * qJD(4); m(7) * (t106 * t45 + t108 * t43 + t154 * t72 - t156 * t76) * t107 + (-m(6) * t44 + t199 - t53 - t55) * t109 + ((-m(7) * t51 + t196 + t198 - t74) * t107 + (m(7) * (t106 * t72 + t108 * t76) - mrSges(4,2) + t197 * (mrSges(5,3) + t184)) * t109) * qJD(3) + (m(5) + m(6)) * t197 * pkin(8) * t157; 0.2e1 * t71 * t53 + 0.2e1 * t33 * t74 + 0.2e1 * t51 * t54 + (t33 * t51 + t43 * t76 + t45 * t72) * t193 - 0.2e1 * pkin(3) * t55 + 0.2e1 * t198 * t44 + (t43 * t192 - t56 - t57 + t58) * t108 + (t45 * t192 + t59 + t60 + t61) * t106 + ((t72 * t192 + t80 + t81 + t82) * t108 + (0.2e1 * mrSges(7,3) * t76 + t77 + t78 - t79) * t106) * qJD(4); (t62 + t67) * qJD(5) + (t30 + t32) * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t10 - t1 * t188) + m(6) * (-pkin(4) * t5 + qJ(5) * t3 + qJD(5) * t17) - t188 * t27 - pkin(4) * t29 - t5 * mrSges(6,1) - t6 * mrSges(5,2) + t7 * mrSges(5,1) - t1 * mrSges(7,1) + t2 * mrSges(7,2) + t3 * mrSges(6,3) - t117 + (Ifges(7,3) * qJD(3) + (-t106 * t149 - t108 * t182) * qJD(4)) * t107 + (-Ifges(7,5) * t108 - t106 * t182) * t157; t84 + ((-mrSges(5,2) + mrSges(6,3)) * t108 + t151 * t106) * t157 + m(6) * (-pkin(4) * t135 + t173) + m(7) * (-t135 * t188 + t173) + (t151 * t108 + (mrSges(5,2) + t183) * t106 + t115 + m(7) * t195) * t155; m(7) * (qJ(5) * t43 + qJD(5) * t76 - t188 * t45) - t45 * mrSges(7,1) + t43 * mrSges(7,2) + t101 + t98 + t100 + qJD(5) * t116 + ((-mrSges(6,2) * pkin(4) + mrSges(7,3) * t188 - Ifges(7,5)) * t108 + (-qJ(5) * t184 - t182) * t106 + (t115 + t73 + t75) * pkin(8)) * qJD(4); 0.2e1 * (qJ(5) * t201 - t183) * qJD(5); -mrSges(6,2) * t137 + t85 + t86 + m(6) * t5 + m(7) * t1 + (t107 * t185 - t143) * qJD(3); t201 * t111; m(7) * t45 + qJD(4) * t116; 0; 0; t200; -m(7) * t158; t199; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
