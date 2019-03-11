% Calculate time derivative of joint inertia matrix for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
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
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:30:33
% EndTime: 2019-03-09 02:30:38
% DurationCPUTime: 2.39s
% Computational Cost: add. (2350->351), mult. (5051->534), div. (0->0), fcn. (4025->6), ass. (0->148)
t171 = qJD(5) + qJD(6);
t104 = sin(qJ(6));
t105 = sin(qJ(5));
t107 = cos(qJ(6));
t108 = cos(qJ(5));
t69 = t104 * t108 + t105 * t107;
t176 = t171 * t69;
t137 = qJD(6) * t107;
t140 = qJD(5) * t108;
t149 = t104 * t105;
t39 = -t107 * t140 - t108 * t137 + t149 * t171;
t172 = mrSges(7,1) * t176 - t39 * mrSges(7,2);
t141 = qJD(5) * t105;
t70 = mrSges(6,1) * t141 + mrSges(6,2) * t140;
t175 = -t172 - t70;
t106 = sin(qJ(4));
t109 = cos(qJ(4));
t142 = qJD(4) * t109;
t126 = t108 * t142;
t174 = -t106 * t141 + t126;
t143 = qJD(4) * t106;
t127 = t105 * t143;
t173 = Ifges(6,6) * t127 + Ifges(6,3) * t142;
t152 = t105 ^ 2 + t108 ^ 2;
t170 = 2 * m(7);
t102 = (qJ(2) - pkin(7));
t169 = -2 * t102;
t168 = 2 * qJD(3);
t167 = m(7) * pkin(5);
t68 = -t107 * t108 + t149;
t166 = -t68 / 0.2e1;
t165 = t69 / 0.2e1;
t164 = -pkin(9) - pkin(8);
t163 = -t105 / 0.2e1;
t162 = pkin(8) * t109;
t161 = t106 * pkin(4);
t103 = pkin(1) + qJ(3);
t79 = -mrSges(6,1) * t108 + mrSges(6,2) * t105;
t159 = t79 - mrSges(5,1);
t147 = t106 * t108;
t77 = t103 + t161 - t162;
t50 = t102 * t147 + t105 * t77;
t158 = Ifges(6,4) * t105;
t157 = Ifges(6,4) * t108;
t156 = Ifges(6,5) * t105;
t155 = Ifges(6,6) * t105;
t154 = Ifges(6,6) * t106;
t153 = Ifges(6,6) * t108;
t99 = t106 ^ 2;
t151 = qJD(2) * t99;
t150 = t102 * t105;
t148 = t105 * t109;
t146 = t108 * t109;
t145 = qJD(2) * t106;
t144 = qJD(3) * t103;
t139 = qJD(5) * t109;
t138 = qJD(6) * t104;
t101 = t109 ^ 2;
t97 = t101 * qJD(2);
t22 = -t109 * t176 + t68 * t143;
t110 = t171 * t68;
t24 = t110 * t109 + t69 * t143;
t136 = Ifges(7,5) * t22 + Ifges(7,6) * t24 + Ifges(7,3) * t142;
t135 = pkin(5) * t141;
t56 = t69 * t109;
t58 = t68 * t109;
t32 = mrSges(7,1) * t56 - mrSges(7,2) * t58;
t66 = (pkin(5) * t105 - t102) * t109;
t134 = m(7) * t66 + t32;
t133 = qJD(5) * t164;
t132 = t105 * t142;
t130 = t105 * t139;
t129 = t106 * t140;
t128 = t108 * t139;
t21 = -qJD(4) * t58 - t106 * t176;
t23 = -qJD(4) * t56 + t110 * t106;
t125 = t23 * mrSges(7,1) - t21 * mrSges(7,2);
t123 = -Ifges(6,5) * t108 + (2 * Ifges(5,4));
t67 = qJD(3) + (pkin(4) * t109 + pkin(8) * t106) * qJD(4);
t19 = t174 * t102 + t105 * t67 + t108 * t145 + t77 * t140;
t64 = t108 * t77;
t49 = -t106 * t150 + t64;
t122 = -t49 * qJD(5) + t19;
t121 = m(6) * pkin(4) - t159;
t120 = Ifges(6,1) * t108 - t158;
t81 = Ifges(6,1) * t105 + t157;
t119 = -Ifges(6,2) * t105 + t157;
t80 = Ifges(6,2) * t108 + t158;
t33 = -pkin(9) * t146 + t64 + (pkin(5) - t150) * t106;
t38 = -pkin(9) * t148 + t50;
t9 = -t104 * t38 + t107 * t33;
t10 = t104 * t33 + t107 * t38;
t82 = t164 * t105;
t83 = t164 * t108;
t45 = t104 * t83 + t107 * t82;
t46 = t104 * t82 - t107 * t83;
t73 = -t106 * mrSges(6,2) - mrSges(6,3) * t148;
t74 = t106 * mrSges(6,1) - mrSges(6,3) * t146;
t118 = t105 * t74 - t108 * t73;
t75 = t105 * t133;
t76 = t108 * t133;
t26 = t45 * qJD(6) + t104 * t76 + t107 * t75;
t27 = -t46 * qJD(6) - t104 * t75 + t107 * t76;
t36 = Ifges(7,6) * t176;
t37 = Ifges(7,5) * t39;
t117 = t27 * mrSges(7,1) - t26 * mrSges(7,2) - t36 - t37;
t114 = -t102 * t142 - t145;
t115 = -t102 * t129 + t108 * t67;
t7 = (pkin(5) * t109 + pkin(9) * t147) * qJD(4) + ((pkin(9) * t109 - t77) * qJD(5) + t114) * t105 + t115;
t112 = t127 - t128;
t8 = t112 * pkin(9) + t19;
t2 = t9 * qJD(6) + t104 * t7 + t107 * t8;
t3 = -t10 * qJD(6) - t104 * t8 + t107 * t7;
t116 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + t136;
t113 = t108 * t143 + t130;
t96 = Ifges(6,5) * t140;
t91 = -pkin(5) * t108 - pkin(4);
t86 = t102 * t97;
t72 = t120 * qJD(5);
t71 = t119 * qJD(5);
t65 = (-mrSges(7,1) * t104 - mrSges(7,2) * t107) * qJD(6) * pkin(5);
t62 = (mrSges(6,1) * t105 + mrSges(6,2) * t108) * t109;
t57 = t68 * t106;
t55 = t69 * t106;
t54 = Ifges(6,5) * t106 + t120 * t109;
t53 = t119 * t109 + t154;
t52 = -mrSges(6,2) * t142 + t112 * mrSges(6,3);
t51 = mrSges(6,1) * t142 + t113 * mrSges(6,3);
t48 = mrSges(7,1) * t106 + mrSges(7,3) * t58;
t47 = -mrSges(7,2) * t106 - mrSges(7,3) * t56;
t44 = Ifges(7,1) * t69 - Ifges(7,4) * t68;
t43 = Ifges(7,4) * t69 - Ifges(7,2) * t68;
t42 = mrSges(7,1) * t68 + mrSges(7,2) * t69;
t41 = -t112 * pkin(5) - qJD(2) * t109 + t102 * t143;
t34 = -t112 * mrSges(6,1) - t113 * mrSges(6,2);
t31 = -t81 * t139 + (Ifges(6,5) * t109 - t120 * t106) * qJD(4);
t30 = -t80 * t139 + (Ifges(6,6) * t109 - t119 * t106) * qJD(4);
t29 = -Ifges(7,1) * t58 - Ifges(7,4) * t56 + Ifges(7,5) * t106;
t28 = -Ifges(7,4) * t58 - Ifges(7,2) * t56 + Ifges(7,6) * t106;
t20 = (-qJD(5) * t77 + t114) * t105 + t115;
t15 = -Ifges(7,1) * t39 - Ifges(7,4) * t176;
t14 = -Ifges(7,4) * t39 - Ifges(7,2) * t176;
t12 = -mrSges(7,2) * t142 + t24 * mrSges(7,3);
t11 = mrSges(7,1) * t142 - t22 * mrSges(7,3);
t6 = -mrSges(7,1) * t24 + mrSges(7,2) * t22;
t5 = Ifges(7,1) * t22 + Ifges(7,4) * t24 + Ifges(7,5) * t142;
t4 = Ifges(7,4) * t22 + Ifges(7,2) * t24 + Ifges(7,6) * t142;
t1 = [(mrSges(4,3) * t168) + 0.2e1 * t10 * t12 + 0.2e1 * t9 * t11 + 0.2e1 * t19 * t73 + 0.2e1 * t2 * t47 + 0.2e1 * t20 * t74 + t22 * t29 + t24 * t28 + 0.2e1 * t3 * t48 + 0.2e1 * t41 * t32 - t56 * t4 + 0.2e1 * t49 * t51 - t58 * t5 + 0.2e1 * t50 * t52 + 0.2e1 * t66 * t6 + 0.2e1 * (m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3) + (-t101 - t99) * mrSges(5,3)) * qJD(2) + 0.2e1 * m(5) * (t102 * t151 + t144 + t86) + 0.2e1 * m(4) * (qJ(2) * qJD(2) + t144) + (t10 * t2 + t3 * t9 + t41 * t66) * t170 + 0.2e1 * m(6) * (t50 * t19 + t49 * t20 + t86) + (mrSges(5,1) * t168 + (-0.2e1 * t103 * mrSges(5,2) + 0.2e1 * t102 * t62 + t105 * t53 + t123 * t106 - t108 * t54) * qJD(4) + t136 + t173) * t106 + (mrSges(5,2) * t168 - 0.2e1 * qJD(2) * t62 + t34 * t169 - t105 * t30 + t108 * t31 + (t106 * (-t153 - t156) - t108 * t53 - t105 * t54) * qJD(5) + (0.2e1 * t103 * mrSges(5,1) - Ifges(7,5) * t58 - Ifges(7,6) * t56 + (-t123 - t155) * t109 + (-0.2e1 * m(6) * (t102 ^ 2) - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t106) * qJD(4)) * t109; -t105 * t52 - t108 * t51 + t68 * t11 - t69 * t12 + t39 * t47 + t176 * t48 + t118 * qJD(5) + (-t109 * mrSges(5,1) + t106 * mrSges(5,2)) * qJD(4) + (-m(5) - m(4)) * qJD(3) + m(7) * (t10 * t39 + t176 * t9 - t2 * t69 + t3 * t68) + m(6) * (-t105 * t19 - t108 * t20 + (t105 * t49 - t108 * t50) * qJD(5)); (t176 * t68 - t39 * t69) * t170; m(4) * qJD(2) - t55 * t11 - t57 * t12 + t21 * t47 + t23 * t48 + (-t118 * qJD(4) - t34 - t6) * t109 + m(7) * (t21 * t10 - t109 * t41 - t57 * t2 + t23 * t9 - t55 * t3) + m(6) * (t50 * t126 - t49 * t132 + t97) + m(5) * (t97 + t151) + (-t74 * t140 - t105 * t51 + m(6) * (-t20 * t105 + t108 * t19 - t49 * t140 - t50 * t141) - t73 * t141 + t108 * t52 + (m(6) * t109 * t169 + t134 + t62) * qJD(4)) * t106; m(7) * (-t176 * t55 - t21 * t69 + t23 * t68 - t39 * t57); (-t57 * t21 - t55 * t23) * t170 + 0.4e1 * (m(6) * (-0.1e1 + t152) / 0.2e1 - m(7) / 0.2e1) * t106 * t142; t91 * t6 + t4 * t166 + t5 * t165 - t56 * t14 / 0.2e1 - t58 * t15 / 0.2e1 + t66 * t172 + t22 * t44 / 0.2e1 + t45 * t11 + t46 * t12 + t26 * t47 + t27 * t48 - t39 * t29 / 0.2e1 - t176 * t28 / 0.2e1 + t41 * t42 + t24 * t43 / 0.2e1 - pkin(4) * t34 + m(7) * (t10 * t26 + t2 * t46 + t27 * t9 + t3 * t45 + t41 * t91) + (t96 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 - qJD(2) * mrSges(5,2) + (-t121 * t102 - Ifges(5,5)) * qJD(4)) * t106 + (t31 / 0.2e1 + t80 * t143 / 0.2e1 - t20 * mrSges(6,3) + (-t154 / 0.2e1 - t50 * mrSges(6,3) - t53 / 0.2e1 + t134 * pkin(5)) * qJD(5) + (-t51 + m(6) * (-t50 * qJD(5) - t20) - qJD(5) * t73) * pkin(8)) * t105 + (-t10 * t176 - t2 * t68 - t3 * t69 + t9 * t39) * mrSges(7,3) + (t30 / 0.2e1 - t81 * t143 / 0.2e1 + qJD(5) * t54 / 0.2e1 + t122 * mrSges(6,3) + (m(6) * t122 - qJD(5) * t74 + t52) * pkin(8)) * t108 + (-t102 * t70 + t71 * t163 + t108 * t72 / 0.2e1 + (-t108 * t80 / 0.2e1 + t81 * t163) * qJD(5) + (-t102 * mrSges(5,2) + Ifges(7,5) * t165 + Ifges(7,6) * t166 - Ifges(5,6) + t156 / 0.2e1 + t153 / 0.2e1) * qJD(4) + t121 * qJD(2)) * t109; m(7) * (t176 * t45 - t26 * t69 + t27 * t68 + t39 * t46); m(7) * (-pkin(5) * t130 + t46 * t21 + t45 * t23 - t26 * t57 - t27 * t55) + (t176 * t57 - t21 * t68 - t23 * t69 - t39 * t55) * mrSges(7,3) + (m(6) * (t152 * t162 - t161) + (m(7) * t91 + t159 + t42) * t106) * qJD(4) + ((t152 * mrSges(6,3) - mrSges(5,2)) * qJD(4) + t175) * t109; -t39 * t44 + t69 * t15 - t176 * t43 - t68 * t14 + (t91 * t135 + t26 * t46 + t27 * t45) * t170 + 0.2e1 * t42 * t135 + 0.2e1 * t91 * t172 - 0.2e1 * pkin(4) * t70 - t80 * t141 + t105 * t72 + (qJD(5) * t81 + t71) * t108 + 0.2e1 * (-t176 * t46 - t26 * t68 - t27 * t69 + t39 * t45) * mrSges(7,3); -Ifges(6,6) * t128 + t20 * mrSges(6,1) - t19 * mrSges(6,2) - t113 * Ifges(6,5) + (-t48 * t138 + t107 * t11 + m(7) * (t10 * t137 + t104 * t2 + t107 * t3 - t9 * t138) + t47 * t137 + t104 * t12) * pkin(5) + t116 + t173; (t104 * t39 + t107 * t176 + (-t104 * t68 - t107 * t69) * qJD(6)) * t167 - t175; -t174 * mrSges(6,2) + (-t129 - t132) * mrSges(6,1) + (t104 * t21 + t107 * t23 + (t104 * t55 - t107 * t57) * qJD(6)) * t167 + t125; t96 + (t79 * pkin(8) - t155) * qJD(5) + (m(7) * (t104 * t26 + t107 * t27 + (-t104 * t45 + t107 * t46) * qJD(6)) + (-t104 * t176 + t107 * t39 + (t104 * t69 - t107 * t68) * qJD(6)) * mrSges(7,3)) * pkin(5) + t117; 0.2e1 * t65; t116; t172; t125; t117; t65; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
