% Calculate time derivative of joint inertia matrix for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:13
% EndTime: 2019-03-09 04:50:20
% DurationCPUTime: 2.80s
% Computational Cost: add. (1517->409), mult. (3470->548), div. (0->0), fcn. (2230->4), ass. (0->169)
t192 = 2 * qJD(2);
t204 = Ifges(6,2) + Ifges(5,3);
t105 = sin(qJ(4));
t107 = cos(qJ(4));
t198 = t105 ^ 2 + t107 ^ 2;
t108 = cos(qJ(3));
t156 = qJD(4) * t108;
t138 = t105 * t156;
t106 = sin(qJ(3));
t160 = qJD(3) * t106;
t113 = t107 * t160 + t138;
t136 = t105 * t160;
t137 = t107 * t156;
t111 = t136 - t137;
t165 = qJ(5) * t105;
t189 = pkin(4) + pkin(5);
t195 = -t189 * t107 - t165;
t203 = m(6) + m(7);
t101 = t106 * qJD(5);
t159 = qJD(3) * t108;
t202 = qJ(5) * t159 + t101;
t18 = mrSges(7,1) * t111 - t113 * mrSges(7,2);
t154 = qJD(5) * t107;
t110 = -pkin(1) - pkin(7);
t164 = qJ(5) * t107;
t115 = -t189 * t105 + t164;
t197 = t110 + t115;
t3 = (t195 * qJD(4) + t154) * t108 - t197 * t160;
t201 = m(7) * t3 + t18;
t155 = qJD(5) * t105;
t32 = t115 * qJD(4) + t155;
t157 = qJD(4) * t107;
t158 = qJD(4) * t105;
t53 = -mrSges(7,1) * t158 + mrSges(7,2) * t157;
t200 = m(7) * t32 + t53;
t120 = pkin(4) * t107 + t165;
t68 = -pkin(3) - t120;
t73 = -t107 * mrSges(6,1) - t105 * mrSges(6,3);
t199 = m(6) * t68 + t73;
t75 = -t107 * mrSges(5,1) + t105 * mrSges(5,2);
t196 = -m(5) * pkin(3) - mrSges(4,1) + t75;
t194 = 0.2e1 * m(7);
t193 = -2 * mrSges(7,3);
t186 = -mrSges(6,1) - mrSges(7,1);
t185 = mrSges(6,2) - mrSges(7,3);
t184 = mrSges(7,2) + mrSges(6,3);
t183 = -Ifges(5,6) - Ifges(7,6);
t182 = -Ifges(6,6) + Ifges(7,6);
t181 = pkin(8) - qJ(6);
t25 = mrSges(5,1) * t159 + t113 * mrSges(5,3);
t26 = -mrSges(6,1) * t159 - t113 * mrSges(6,2);
t180 = -t25 + t26;
t28 = -mrSges(5,2) * t159 + t111 * mrSges(5,3);
t29 = t111 * mrSges(6,2) + mrSges(6,3) * t159;
t179 = t28 + t29;
t163 = t105 * t108;
t62 = -mrSges(5,2) * t106 - mrSges(5,3) * t163;
t66 = -mrSges(6,2) * t163 + mrSges(6,3) * t106;
t178 = -t62 - t66;
t161 = t107 * t108;
t64 = t106 * mrSges(5,1) - mrSges(5,3) * t161;
t65 = -t106 * mrSges(6,1) + mrSges(6,2) * t161;
t177 = -t64 + t65;
t140 = t107 * t159;
t175 = qJ(5) * t140 + t107 * t101;
t174 = t113 * mrSges(7,3);
t162 = t106 * t110;
t67 = pkin(3) * t106 - pkin(8) * t108 + qJ(2);
t31 = t105 * t67 + t107 * t162;
t172 = Ifges(5,4) * t105;
t171 = Ifges(5,4) * t107;
t170 = Ifges(7,4) * t105;
t169 = Ifges(7,4) * t107;
t168 = Ifges(6,5) * t105;
t167 = Ifges(6,5) * t107;
t166 = t106 * Ifges(7,5);
t153 = qJD(6) * t107;
t152 = -mrSges(5,1) + t186;
t151 = mrSges(5,2) - t184;
t149 = Ifges(6,4) + Ifges(5,5) - Ifges(7,5);
t121 = Ifges(6,3) * t105 + t167;
t34 = Ifges(6,6) * t106 + t121 * t108;
t122 = Ifges(7,2) * t105 + t169;
t35 = -t106 * Ifges(7,6) + t122 * t108;
t123 = -Ifges(5,2) * t105 + t171;
t36 = t106 * Ifges(5,6) + t123 * t108;
t148 = t36 - t34 - t35;
t139 = t110 * t159;
t51 = qJD(2) + (pkin(3) * t108 + pkin(8) * t106) * qJD(3);
t147 = t105 * t51 + t107 * t139 + t67 * t157;
t135 = t106 * t157;
t146 = t105 * t139 + t110 * t135 + t67 * t158;
t61 = mrSges(7,2) * t106 + mrSges(7,3) * t163;
t145 = t61 - t178;
t63 = -t106 * mrSges(7,1) - mrSges(7,3) * t161;
t144 = t63 + t177;
t22 = t106 * qJ(5) + t31;
t142 = -Ifges(5,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t89 = t105 * t162;
t141 = t105 * t159;
t76 = t181 * t107;
t30 = t107 * t67 - t89;
t133 = Ifges(5,6) * t136 + Ifges(6,6) * t137 + t159 * t204;
t77 = -Ifges(6,3) * t107 + t168;
t78 = -Ifges(7,2) * t107 + t170;
t79 = Ifges(5,2) * t107 + t172;
t132 = t78 / 0.2e1 - t79 / 0.2e1 + t77 / 0.2e1;
t80 = Ifges(7,1) * t105 - t169;
t81 = Ifges(6,1) * t105 - t167;
t82 = Ifges(5,1) * t105 + t171;
t131 = -t80 / 0.2e1 - t81 / 0.2e1 - t82 / 0.2e1;
t6 = -qJD(4) * t89 + t147;
t4 = t6 + t202;
t7 = t107 * t51 - t146;
t5 = -pkin(4) * t159 - t7;
t130 = t105 * t5 + t107 * t4;
t129 = -t105 * t7 + t107 * t6;
t128 = mrSges(5,1) * t105 + mrSges(5,2) * t107;
t127 = mrSges(6,1) * t105 - mrSges(6,3) * t107;
t126 = Ifges(5,1) * t107 - t172;
t125 = Ifges(6,1) * t107 + t168;
t124 = Ifges(7,1) * t107 + t170;
t119 = pkin(4) * t105 - t164;
t118 = t107 * (m(6) * pkin(8) + t185);
t117 = -t110 + t119;
t116 = m(6) * t120;
t37 = t124 * t108 - t166;
t38 = Ifges(6,4) * t106 + t125 * t108;
t39 = Ifges(5,5) * t106 + t126 * t108;
t114 = -t149 * t106 - t37 - t38 - t39;
t100 = Ifges(6,4) * t157;
t99 = Ifges(5,5) * t157;
t97 = Ifges(6,6) * t158;
t74 = t107 * mrSges(7,1) + t105 * mrSges(7,2);
t72 = t181 * t105;
t60 = t126 * qJD(4);
t59 = t125 * qJD(4);
t58 = t124 * qJD(4);
t57 = t123 * qJD(4);
t56 = t122 * qJD(4);
t55 = t121 * qJD(4);
t54 = t128 * qJD(4);
t52 = t127 * qJD(4);
t49 = pkin(3) - t195;
t48 = t128 * t108;
t47 = (-mrSges(7,1) * t105 + mrSges(7,2) * t107) * t108;
t46 = t127 * t108;
t42 = qJD(4) * t76 - qJD(6) * t105;
t41 = t119 * qJD(4) - t155;
t40 = -t181 * t158 - t153;
t33 = t117 * t108;
t27 = mrSges(7,2) * t159 - t111 * mrSges(7,3);
t24 = -mrSges(7,1) * t159 + t174;
t23 = -pkin(4) * t106 - t30;
t21 = t197 * t108;
t19 = -t111 * mrSges(5,1) - t113 * mrSges(5,2);
t17 = -t111 * mrSges(6,1) + t113 * mrSges(6,3);
t16 = qJ(6) * t163 + t22;
t15 = t89 + (-qJ(6) * t108 - t67) * t107 - t189 * t106;
t14 = -t82 * t156 + (Ifges(5,5) * t108 - t126 * t106) * qJD(3);
t13 = -t81 * t156 + (Ifges(6,4) * t108 - t125 * t106) * qJD(3);
t12 = -t80 * t156 + (-Ifges(7,5) * t108 - t124 * t106) * qJD(3);
t11 = -t79 * t156 + (Ifges(5,6) * t108 - t123 * t106) * qJD(3);
t10 = -t78 * t156 + (-Ifges(7,6) * t108 - t122 * t106) * qJD(3);
t9 = -t77 * t156 + (Ifges(6,6) * t108 - t121 * t106) * qJD(3);
t8 = (t120 * qJD(4) - t154) * t108 - t117 * t160;
t2 = qJ(6) * t137 + (qJD(6) * t108 + (-qJ(6) * qJD(3) - qJD(4) * t110) * t106) * t105 + t147 + t202;
t1 = (qJ(6) * t160 - t51) * t107 + (qJ(6) * t158 - t189 * qJD(3) - t153) * t108 + t146;
t20 = [0.2e1 * t2 * t61 + 0.2e1 * t6 * t62 + 0.2e1 * t1 * t63 + 0.2e1 * t7 * t64 + 0.2e1 * t5 * t65 + 0.2e1 * t4 * t66 + 0.2e1 * t8 * t46 + 0.2e1 * t3 * t47 + 0.2e1 * t15 * t24 + 0.2e1 * t23 * t26 + 0.2e1 * t16 * t27 + 0.2e1 * t22 * t29 + 0.2e1 * t30 * t25 + 0.2e1 * t31 * t28 + 0.2e1 * t33 * t17 + 0.2e1 * t21 * t18 + (mrSges(3,3) + (m(3) + m(4)) * qJ(2)) * t192 + 0.2e1 * m(6) * (t22 * t4 + t23 * t5 + t33 * t8) + (t1 * t15 + t16 * t2 + t21 * t3) * t194 + 0.2e1 * m(5) * (t30 * t7 + t31 * t6) + (mrSges(4,1) * t192 + (-0.2e1 * qJ(2) * mrSges(4,2) + 0.2e1 * Ifges(4,4) * t106 + 0.2e1 * t110 * t48 + (t182 * t106 + t148) * t105 + t114 * t107) * qJD(3) + t133) * t106 + (mrSges(4,2) * t192 - 0.2e1 * t110 * t19 + (t12 + t13 + t14) * t107 + (t10 - t11 + t9) * t105 + ((t183 * t106 - t148) * t107 + t114 * t105) * qJD(4) + (0.2e1 * qJ(2) * mrSges(4,1) - 0.2e1 * Ifges(4,4) * t108 + t149 * t161 + (-Ifges(5,6) - t182) * t163 + (-0.2e1 * m(5) * t110 ^ 2 - (2 * Ifges(4,1)) + (2 * Ifges(4,2)) + (2 * Ifges(7,3)) + t204) * t106) * qJD(3)) * t108; (-t17 - t19 - m(6) * t8 + (t145 * t107 + t144 * t105 + m(6) * (t105 * t23 + t107 * t22) + m(7) * (t105 * t15 + t107 * t16) + m(5) * (-t105 * t30 + t107 * t31)) * qJD(3) + t201) * t108 + ((t27 + t179) * t107 + (t24 + t180) * t105 + (t46 - t47 + t48) * qJD(3) + (-t145 * t105 + t144 * t107) * qJD(4) + m(6) * (qJD(3) * t33 + t23 * t157 - t22 * t158 + t130) + m(7) * (-qJD(3) * t21 + t1 * t105 + t107 * t2 + t15 * t157 - t16 * t158) + m(5) * (-t30 * t157 - t31 * t158 + t129 - 0.2e1 * t139)) * t106; 0.4e1 * (m(6) / 0.2e1 + m(7) / 0.2e1 + m(5) / 0.2e1) * (-0.1e1 + t198) * t106 * t159; -pkin(3) * t19 + t68 * t17 + t49 * t18 + t21 * t53 + t72 * t24 + t76 * t27 + t3 * t74 + t32 * t47 + t33 * t52 + t40 * t61 + t41 * t46 + t42 * t63 + t8 * t73 + m(7) * (t1 * t72 + t15 * t42 + t16 * t40 + t2 * t76 + t21 * t32 + t3 * t49) + m(6) * (t33 * t41 + t68 * t8) + (t6 * mrSges(5,3) - t2 * mrSges(7,3) + t4 * mrSges(6,2) - t9 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1) * t107 + (t5 * mrSges(6,2) - t7 * mrSges(5,3) - t1 * mrSges(7,3) + t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1) * t105 + (m(5) * t129 + m(6) * t130 + t180 * t105 + t179 * t107) * pkin(8) + (t100 / 0.2e1 + t97 / 0.2e1 + t99 / 0.2e1 + (-t132 * t105 + t131 * t107 + t196 * t110 - Ifges(4,5)) * qJD(3)) * t106 + (-t110 * t54 + (t58 / 0.2e1 + t59 / 0.2e1 + t60 / 0.2e1) * t107 + (t55 / 0.2e1 + t56 / 0.2e1 - t57 / 0.2e1) * t105 + (-t110 * mrSges(4,2) - Ifges(4,6) + (-Ifges(6,6) / 0.2e1 - t142) * t107 + (-Ifges(7,5) / 0.2e1 + Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t105) * qJD(3)) * t108 + ((t37 / 0.2e1 + t38 / 0.2e1 + t39 / 0.2e1 - t30 * mrSges(5,3) - t15 * mrSges(7,3) + t23 * mrSges(6,2) - t166 / 0.2e1 + t132 * t108) * t107 + (t34 / 0.2e1 + t35 / 0.2e1 - t36 / 0.2e1 - t31 * mrSges(5,3) + t16 * mrSges(7,3) - t22 * mrSges(6,2) + t142 * t106 + t131 * t108) * t105 + (t177 * t107 + t178 * t105 + m(5) * (-t31 * t105 - t30 * t107) + m(6) * (-t22 * t105 + t23 * t107)) * pkin(8)) * qJD(4); m(7) * (t105 * t42 + t107 * t40 + t157 * t72 - t158 * t76) * t106 + (-m(6) * t41 + t200 - t52 - t54) * t108 + ((-m(7) * t49 + t196 + t199 - t74) * t106 + (m(7) * (t105 * t72 + t107 * t76) - mrSges(4,2) + t198 * (mrSges(5,3) + t185)) * t108) * qJD(3) + (m(6) + m(5)) * t198 * pkin(8) * t159; 0.2e1 * t68 * t52 + 0.2e1 * t32 * t74 + 0.2e1 * t49 * t53 + (t32 * t49 + t40 * t76 + t42 * t72) * t194 - 0.2e1 * pkin(3) * t54 + 0.2e1 * t199 * t41 + (t40 * t193 - t55 - t56 + t57) * t107 + (t42 * t193 + t58 + t59 + t60) * t105 + ((t193 * t72 + t80 + t81 + t82) * t107 + (0.2e1 * mrSges(7,3) * t76 + t77 + t78 - t79) * t105) * qJD(4); t7 * mrSges(5,1) - t5 * mrSges(6,1) - t1 * mrSges(7,1) - t6 * mrSges(5,2) + t2 * mrSges(7,2) + t4 * mrSges(6,3) - pkin(4) * t26 - t189 * t24 + (t61 + t66) * qJD(5) + (t27 + t29) * qJ(5) + m(7) * (qJ(5) * t2 + qJD(5) * t16 - t1 * t189) + m(6) * (-pkin(4) * t5 + qJ(5) * t4 + qJD(5) * t22) + (-t105 * t149 + t107 * t183) * t156 + (Ifges(7,3) * t108 + (t105 * t182 - t107 * t149) * t106) * qJD(3) + t133; t152 * t141 + m(6) * (-pkin(4) * t141 + t175) + m(7) * (-t141 * t189 + t175) - t151 * t140 + (m(7) * t195 + t151 * t105 + t152 * t107 - t116) * t106 * qJD(4); m(7) * (qJ(5) * t40 + qJD(5) * t76 - t189 * t42) - t42 * mrSges(7,1) + t40 * mrSges(7,2) + t99 + t100 + t97 + qJD(5) * t118 + ((-pkin(4) * mrSges(6,2) + mrSges(7,3) * t189 - Ifges(7,5)) * t107 + (-qJ(5) * t185 + t183) * t105 + (-t116 + t73 + t75) * pkin(8)) * qJD(4); 0.2e1 * (t203 * qJ(5) + t184) * qJD(5); -mrSges(6,2) * t138 + m(6) * t5 + m(7) * t1 + (-t107 * t106 * mrSges(6,2) + t108 * t186) * qJD(3) + t174; t203 * (t135 + t141); m(7) * t42 + qJD(4) * t118; 0; 0; t201; -m(7) * t160; t200; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
