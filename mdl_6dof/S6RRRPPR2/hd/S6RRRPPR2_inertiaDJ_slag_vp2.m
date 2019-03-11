% Calculate time derivative of joint inertia matrix for
% S6RRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 15:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:24:59
% EndTime: 2019-03-09 15:25:05
% DurationCPUTime: 2.42s
% Computational Cost: add. (5450->305), mult. (12002->442), div. (0->0), fcn. (11531->8), ass. (0->139)
t116 = sin(qJ(6));
t191 = -t116 / 0.2e1;
t119 = cos(qJ(6));
t190 = t119 / 0.2e1;
t168 = -mrSges(6,2) + mrSges(5,1);
t150 = qJD(6) * t119;
t114 = sin(pkin(10));
t115 = cos(pkin(10));
t186 = qJD(2) + qJD(3);
t117 = sin(qJ(3));
t118 = sin(qJ(2));
t120 = cos(qJ(3));
t121 = cos(qJ(2));
t88 = -t117 * t118 + t120 * t121;
t70 = t186 * t88;
t89 = t117 * t121 + t120 * t118;
t71 = t186 * t89;
t50 = t114 * t70 + t115 * t71;
t163 = t116 * t50;
t66 = t114 * t89 - t115 * t88;
t131 = t66 * t150 + t163;
t167 = Ifges(7,4) * t116;
t100 = Ifges(7,1) * t119 - t167;
t166 = Ifges(7,4) * t119;
t99 = -Ifges(7,2) * t116 + t166;
t189 = t100 * t116 + t119 * t99;
t98 = mrSges(7,1) * t116 + mrSges(7,2) * t119;
t187 = mrSges(6,3) + t98;
t67 = t114 * t88 + t115 * t89;
t109 = -pkin(2) * t121 - pkin(1);
t74 = -t88 * pkin(3) + t109;
t127 = -t67 * qJ(5) + t74;
t177 = pkin(4) + pkin(9);
t23 = t177 * t66 + t127;
t176 = -pkin(8) - pkin(7);
t101 = t176 * t118;
t102 = t176 * t121;
t157 = t102 * t117;
t72 = t120 * t101 + t157;
t61 = -qJ(4) * t89 + t72;
t73 = t117 * t101 - t120 * t102;
t62 = qJ(4) * t88 + t73;
t33 = t114 * t62 - t115 * t61;
t24 = pkin(5) * t67 + t33;
t13 = t116 * t24 + t119 * t23;
t155 = t13 * qJD(6);
t51 = -t114 * t71 + t115 * t70;
t64 = qJD(2) * t118 * pkin(2) + pkin(3) * t71;
t126 = -qJ(5) * t51 - qJD(5) * t67 + t64;
t10 = t177 * t50 + t126;
t143 = qJD(2) * t176;
t96 = t118 * t143;
t97 = t121 * t143;
t54 = -qJD(3) * t73 - t117 * t96 + t120 * t97;
t124 = -qJ(4) * t70 - qJD(4) * t89 + t54;
t152 = qJD(3) * t120;
t53 = qJD(3) * t157 + t101 * t152 + t117 * t97 + t120 * t96;
t31 = -qJ(4) * t71 + qJD(4) * t88 + t53;
t17 = t114 * t31 - t115 * t124;
t8 = pkin(5) * t51 + t17;
t2 = -t10 * t116 + t119 * t8 - t155;
t141 = -t2 - t155;
t185 = 2 * m(5);
t184 = 2 * m(6);
t183 = 2 * m(7);
t19 = pkin(4) * t50 + t126;
t182 = -0.2e1 * t19;
t181 = 0.2e1 * t64;
t151 = qJD(6) * t116;
t93 = -mrSges(7,1) * t150 + mrSges(7,2) * t151;
t180 = -0.2e1 * t93;
t179 = 0.2e1 * t109;
t175 = Ifges(7,6) * t67;
t12 = -t116 * t23 + t119 * t24;
t1 = qJD(6) * t12 + t10 * t119 + t116 * t8;
t174 = t1 * t116;
t156 = t115 * t117;
t165 = pkin(2) * qJD(3);
t79 = (t114 * t120 + t156) * t165;
t173 = t33 * t79;
t172 = t67 * t79;
t104 = t114 * t117 * pkin(2);
t80 = pkin(2) * t115 * t152 - qJD(3) * t104;
t75 = qJD(5) + t80;
t108 = pkin(2) * t120 + pkin(3);
t82 = pkin(2) * t156 + t108 * t114;
t77 = qJ(5) + t82;
t171 = t75 * t77;
t170 = t80 * mrSges(5,2);
t164 = t116 * t13;
t162 = t116 * t66;
t161 = t119 * t50;
t160 = t119 * t66;
t154 = t116 ^ 2 + t119 ^ 2;
t149 = 0.2e1 * t50;
t148 = 0.2e1 * t121;
t147 = t12 * t151;
t146 = t66 * t151;
t144 = t160 / 0.2e1;
t107 = -pkin(3) * t115 - pkin(4);
t140 = t154 * mrSges(7,3);
t139 = t154 * t79;
t138 = t131 * Ifges(7,5) + Ifges(7,6) * t161 + Ifges(7,3) * t51;
t81 = t108 * t115 - t104;
t78 = -pkin(4) - t81;
t18 = t114 * t124 + t115 * t31;
t34 = t114 * t61 + t115 * t62;
t137 = t17 * t33 + t18 * t34;
t136 = mrSges(7,1) * t119 - mrSges(7,2) * t116;
t135 = Ifges(7,1) * t116 + t166;
t134 = Ifges(7,2) * t119 + t167;
t133 = -Ifges(7,5) * t116 - Ifges(7,6) * t119;
t106 = pkin(3) * t114 + qJ(5);
t132 = qJD(5) * t77 + t106 * t75;
t130 = t146 - t161;
t128 = (-mrSges(4,1) * t117 - mrSges(4,2) * t120) * t165;
t94 = t134 * qJD(6);
t95 = t135 * qJD(6);
t125 = -t189 * qJD(6) + t116 * t94 - t119 * t95;
t20 = -mrSges(7,2) * t51 - mrSges(7,3) * t130;
t21 = mrSges(7,1) * t51 - mrSges(7,3) * t131;
t123 = m(7) * (t174 + t119 * t2 + (-t116 * t12 + t119 * t13) * qJD(6)) + t116 * t20 + t119 * t21;
t25 = -pkin(5) * t66 + t34;
t28 = t134 * t66 + t175;
t29 = Ifges(7,5) * t67 + t135 * t66;
t6 = Ifges(7,4) * t131 - Ifges(7,2) * t130 + Ifges(7,6) * t51;
t7 = Ifges(7,1) * t131 - Ifges(7,4) * t130 + Ifges(7,5) * t51;
t9 = -t50 * pkin(5) + t18;
t122 = -t95 * t162 / 0.2e1 - t28 * t150 / 0.2e1 + t7 * t190 + t6 * t191 - t25 * t93 + t9 * t98 + Ifges(4,5) * t70 - Ifges(4,6) * t71 - t53 * mrSges(4,2) + t54 * mrSges(4,1) + mrSges(7,3) * t147 + (-mrSges(5,2) + mrSges(6,3)) * t18 - t94 * t144 + (Ifges(7,5) * t190 + Ifges(7,6) * t191 - Ifges(6,4) + Ifges(5,5)) * t51 - t168 * t17 - (t66 * t99 + t29) * t151 / 0.2e1 + (t67 * t133 / 0.2e1 + t100 * t144) * qJD(6) + (-Ifges(5,6) + Ifges(6,5) + t189 / 0.2e1) * t50;
t105 = -pkin(9) + t107;
t76 = -pkin(9) + t78;
t44 = t51 * mrSges(6,3);
t43 = t51 * mrSges(5,2);
t40 = -mrSges(7,2) * t67 + mrSges(7,3) * t160;
t39 = mrSges(7,1) * t67 - mrSges(7,3) * t162;
t38 = t136 * t66;
t32 = t66 * pkin(4) + t127;
t14 = mrSges(7,1) * t130 + mrSges(7,2) * t131;
t3 = [(mrSges(4,1) * t71 + mrSges(4,2) * t70) * t179 - 0.2e1 * t88 * Ifges(4,2) * t71 + 0.2e1 * t70 * t89 * Ifges(4,1) + 0.2e1 * t74 * (t50 * mrSges(5,1) + t43) + 0.2e1 * t32 * (-t50 * mrSges(6,2) - t44) - 0.2e1 * t9 * t38 + 0.2e1 * t2 * t39 + 0.2e1 * t1 * t40 + 0.2e1 * t25 * t14 + 0.2e1 * t13 * t20 + 0.2e1 * t12 * t21 + t28 * t161 + t29 * t163 + 0.2e1 * m(4) * (t53 * t73 + t54 * t72) + (t64 * t74 + t137) * t185 + (t19 * t32 + t137) * t184 + (t1 * t13 + t12 * t2 + t25 * t9) * t183 + ((-mrSges(3,2) * pkin(1) + Ifges(3,4) * t121) * t148 + (0.2e1 * pkin(2) * (-mrSges(4,1) * t88 + mrSges(4,2) * t89) - 0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t179 - 0.2e1 * Ifges(3,4) * t118 + (-Ifges(3,2) + Ifges(3,1)) * t148) * t118) * qJD(2) + (mrSges(5,2) * t181 + mrSges(6,3) * t182 + (-Ifges(5,4) - Ifges(6,6)) * t149 + ((2 * Ifges(5,1)) + (2 * Ifges(6,2)) + Ifges(7,3)) * t51 + t138) * t67 + (mrSges(5,1) * t181 + mrSges(6,2) * t182 + t116 * t7 + t119 * t6 + (Ifges(6,3) + Ifges(5,2)) * t149 + (-0.2e1 * Ifges(5,4) - 0.2e1 * Ifges(6,6) - t133) * t51 + (t119 * t29 + (-t28 - t175) * t116) * qJD(6)) * t66 + 0.2e1 * (t70 * t88 - t71 * t89) * Ifges(4,4) + 0.2e1 * (t53 * t88 - t54 * t89 - t70 * t72 - t71 * t73) * mrSges(4,3) + (t17 * t67 - t18 * t66 + t33 * t51 - t34 * t50) * (2 * mrSges(6,1) + 2 * mrSges(5,3)); t122 + m(5) * (-t17 * t81 + t18 * t82 + t34 * t80 + t173) + m(6) * (t17 * t78 + t18 * t77 + t34 * t75 + t173) + m(7) * (t79 * t164 + t25 * t75 + t77 * t9 + (-t147 + t174) * t76) + t77 * t14 - t75 * t38 + (-t50 * t77 + t51 * t78 - t66 * t75 + t172) * mrSges(6,1) + (t141 * mrSges(7,3) + (m(7) * t12 + t39) * t79 + (-m(7) * t141 + qJD(6) * t40 + t21) * t76) * t119 + (m(4) * (t117 * t53 + t120 * t54 + (-t117 * t72 + t120 * t73) * qJD(3)) + (-t117 * t71 - t120 * t70 + (t117 * t89 + t120 * t88) * qJD(3)) * mrSges(4,3)) * pkin(2) + (-t50 * t82 - t51 * t81 - t66 * t80 + t172) * mrSges(5,3) + (-t1 * mrSges(7,3) + t79 * t40 + (-qJD(6) * t39 + t20) * t76) * t116 + (Ifges(3,5) * t121 - Ifges(3,6) * t118 + (-mrSges(3,1) * t121 + mrSges(3,2) * t118) * pkin(7)) * qJD(2); -0.2e1 * t170 + t77 * t180 + (t139 * t76 + t171) * t183 + t171 * t184 + t80 * t82 * t185 + t125 + 0.2e1 * t187 * t75 + 0.2e1 * t128 + (t184 * t78 - t185 * t81 - 0.2e1 * t140 - 0.2e1 * t168) * t79; t122 + t106 * t14 - qJD(5) * t38 + (t119 * t141 - t174) * mrSges(7,3) + m(7) * (qJD(5) * t25 + t106 * t9) + (m(5) * (t114 * t18 - t115 * t17) + (-t114 * t50 - t115 * t51) * mrSges(5,3)) * pkin(3) + (t150 * t40 - t151 * t39 + t123) * t105 + m(6) * (qJD(5) * t34 + t106 * t18 + t107 * t17) + (-qJD(5) * t66 - t106 * t50 + t107 * t51) * mrSges(6,1); -t170 + (-t77 - t106) * t93 + t128 + (-t140 - t168) * t79 + m(7) * (t105 * t139 + t132) + m(6) * (t107 * t79 + t132) + m(5) * (t114 * t80 - t115 * t79) * pkin(3) + t125 + t187 * (t75 + qJD(5)); t106 * t180 + 0.2e1 * ((m(6) + m(7)) * t106 + t187) * qJD(5) + t125; -t116 * t21 + t119 * t20 + t43 - t44 + t168 * t50 + (-t116 * t40 - t119 * t39) * qJD(6) + m(7) * (t1 * t119 - t116 * t2 + (-t119 * t12 - t164) * qJD(6)) + m(6) * t19 + m(5) * t64; 0; 0; 0; t51 * mrSges(6,1) + (-t116 * t39 + t119 * t40) * qJD(6) + m(6) * t17 + t123; 0.2e1 * (m(7) * t154 / 0.2e1 + m(6) / 0.2e1) * t79; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t146 + t138; t136 * t79 + ((-mrSges(7,2) * t76 - Ifges(7,6)) * t119 + (-mrSges(7,1) * t76 - Ifges(7,5)) * t116) * qJD(6); (-t105 * t98 + t133) * qJD(6); t93; -t98 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
