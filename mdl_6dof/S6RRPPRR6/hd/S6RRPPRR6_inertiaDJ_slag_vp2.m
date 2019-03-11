% Calculate time derivative of joint inertia matrix for
% S6RRPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:13:07
% EndTime: 2019-03-09 09:13:13
% DurationCPUTime: 2.44s
% Computational Cost: add. (4593->323), mult. (9525->471), div. (0->0), fcn. (9248->8), ass. (0->136)
t94 = sin(qJ(6));
t97 = cos(qJ(6));
t73 = -t97 * mrSges(7,1) + mrSges(7,2) * t94;
t146 = mrSges(6,1) - t73;
t172 = -m(7) * pkin(5) - t146;
t144 = pkin(7) - qJ(4);
t96 = sin(qJ(2));
t177 = t144 * t96;
t98 = cos(qJ(2));
t75 = t144 * t98;
t92 = sin(pkin(10));
t93 = cos(pkin(10));
t43 = t177 * t92 + t93 * t75;
t158 = cos(qJ(5));
t95 = sin(qJ(5));
t104 = t158 * t93 - t95 * t92;
t99 = -pkin(2) - pkin(3);
t70 = -qJ(3) * t92 + t93 * t99;
t65 = -pkin(4) + t70;
t71 = t93 * qJ(3) + t92 * t99;
t39 = t158 * t65 - t95 * t71;
t28 = t104 * qJD(3) + t39 * qJD(5);
t179 = t28 * mrSges(6,2);
t40 = t158 * t71 + t95 * t65;
t178 = -t98 * pkin(2) - t96 * qJ(3);
t176 = m(4) * pkin(7) + mrSges(4,2);
t110 = t92 * t98 - t93 * t96;
t42 = t177 * t93 - t75 * t92;
t108 = pkin(8) * t110 + t42;
t109 = t92 * t96 + t93 * t98;
t33 = -pkin(8) * t109 + t43;
t175 = t158 * t108 - t33 * t95;
t174 = m(7) * pkin(9) + mrSges(7,3);
t105 = -t109 * t158 + t110 * t95;
t133 = qJD(6) * t97;
t53 = t110 * qJD(2);
t54 = t109 * qJD(2);
t22 = t105 * qJD(5) + t158 * t54 - t95 * t53;
t36 = -t109 * t95 - t110 * t158;
t23 = t36 * qJD(5) + t158 * t53 + t95 * t54;
t31 = t43 * qJD(2) + t110 * qJD(4);
t100 = -t54 * pkin(8) + t31;
t32 = t93 * (-qJD(2) * t177 - qJD(4) * t98) + t92 * (qJD(2) * t75 - qJD(4) * t96);
t27 = -pkin(8) * t53 + t32;
t6 = t175 * qJD(5) + t95 * t100 + t158 * t27;
t114 = mrSges(7,1) * t94 + mrSges(7,2) * t97;
t66 = t114 * qJD(6);
t76 = Ifges(7,5) * t94 + Ifges(7,6) * t97;
t134 = qJD(6) * t94;
t82 = Ifges(7,6) * t134;
t173 = (t76 / 0.2e1 - Ifges(6,6)) * t23 + Ifges(6,5) * t22 - t175 * t66 - t105 * (Ifges(7,5) * t133 - t82) / 0.2e1 - t6 * mrSges(6,2);
t68 = Ifges(7,4) * t133 - Ifges(7,2) * t134;
t69 = Ifges(7,1) * t133 - Ifges(7,4) * t134;
t156 = Ifges(7,4) * t94;
t77 = Ifges(7,2) * t97 + t156;
t155 = Ifges(7,4) * t97;
t78 = Ifges(7,1) * t94 + t155;
t101 = -(t94 * t77 - t97 * t78) * qJD(6) + t97 * t68 + t94 * t69;
t171 = 2 * m(6);
t170 = 0.2e1 * m(7);
t169 = -2 * pkin(1);
t17 = t95 * t108 + t158 * t33;
t7 = t17 * qJD(5) - t158 * t100 + t95 * t27;
t168 = 0.2e1 * t7;
t167 = -2 * mrSges(6,3);
t166 = -0.2e1 * t175;
t138 = qJD(2) * t96;
t137 = qJD(2) * t98;
t142 = qJ(3) * t137 + t96 * qJD(3);
t44 = t99 * t138 + t142;
t34 = pkin(4) * t53 + t44;
t165 = 0.2e1 * t34;
t164 = -0.2e1 * t66;
t72 = -pkin(1) + t178;
t163 = 0.2e1 * t72;
t161 = -t36 / 0.2e1;
t160 = -t77 / 0.2e1;
t159 = t175 * t7;
t157 = mrSges(7,3) * t36;
t154 = Ifges(7,5) * t97;
t61 = t158 * t92 + t95 * t93;
t29 = t61 * qJD(3) + t40 * qJD(5);
t153 = t175 * t29;
t152 = t29 * t104;
t52 = t61 * qJD(5);
t150 = t52 * t104;
t149 = t104 * t66;
t148 = t97 * t22;
t147 = qJD(6) / 0.2e1;
t145 = Ifges(3,4) - Ifges(4,5);
t143 = Ifges(7,5) * t148 + Ifges(7,3) * t23;
t141 = t94 ^ 2 + t97 ^ 2;
t57 = t98 * pkin(3) - t72;
t41 = pkin(4) * t109 + t57;
t15 = -pkin(5) * t105 - pkin(9) * t36 + t41;
t8 = t15 * t97 - t17 * t94;
t140 = t8 * qJD(6);
t9 = t15 * t94 + t17 * t97;
t139 = t9 * qJD(6);
t38 = -pkin(9) + t40;
t135 = qJD(6) * t38;
t131 = 0.2e1 * t98;
t130 = t36 * t134;
t128 = t141 * mrSges(7,3);
t127 = t141 * t28;
t51 = t104 * qJD(5);
t126 = t141 * t51;
t125 = t141 * t61;
t124 = t53 * mrSges(5,1) + t54 * mrSges(5,2);
t123 = t23 * mrSges(6,1) + t22 * mrSges(6,2);
t122 = -Ifges(7,6) * t94 - (2 * Ifges(6,4));
t12 = pkin(5) * t23 - pkin(9) * t22 + t34;
t1 = t12 * t94 + t6 * t97 + t140;
t121 = -t1 + t140;
t2 = t12 * t97 - t6 * t94 - t139;
t120 = t2 + t139;
t106 = t130 - t148;
t107 = t36 * t133 + t94 * t22;
t3 = -t106 * Ifges(7,4) - t107 * Ifges(7,2) + Ifges(7,6) * t23;
t119 = t3 / 0.2e1 + t22 * t78 / 0.2e1;
t4 = -t106 * Ifges(7,1) - t107 * Ifges(7,4) + Ifges(7,5) * t23;
t118 = t4 / 0.2e1 + t22 * t160;
t117 = -t8 * t94 + t9 * t97;
t116 = -t104 * t7 - t175 * t52;
t115 = -t98 * mrSges(4,1) - t96 * mrSges(4,3);
t24 = mrSges(7,2) * t105 - t94 * t157;
t25 = -mrSges(7,1) * t105 - t97 * t157;
t113 = t97 * t24 - t94 * t25;
t37 = pkin(5) - t39;
t19 = t114 * t36;
t14 = -Ifges(7,5) * t105 + (Ifges(7,1) * t97 - t156) * t36;
t13 = -Ifges(7,6) * t105 + (-Ifges(7,2) * t94 + t155) * t36;
t11 = -mrSges(7,2) * t23 - t107 * mrSges(7,3);
t10 = mrSges(7,1) * t23 + t106 * mrSges(7,3);
t5 = t107 * mrSges(7,1) - t106 * mrSges(7,2);
t16 = [0.2e1 * t44 * (mrSges(5,1) * t109 - mrSges(5,2) * t110) - 0.2e1 * t110 * t54 * Ifges(5,1) + 0.2e1 * t53 * Ifges(5,2) * t109 + 0.2e1 * t57 * t124 + 0.2e1 * t41 * t123 + 0.2e1 * t1 * t24 + 0.2e1 * t2 * t25 + t5 * t166 + t19 * t168 + 0.2e1 * t8 * t10 + 0.2e1 * t9 * t11 + t17 * t23 * t167 + (mrSges(6,3) * t166 - t94 * t13 + t97 * t14) * t22 + (t1 * t9 + t2 * t8 - t159) * t170 + (t17 * t6 + t34 * t41 - t159) * t171 + 0.2e1 * m(5) * (t31 * t42 + t32 * t43 + t44 * t57) - (mrSges(6,1) * t165 + t6 * t167 + ((2 * Ifges(6,2)) + Ifges(7,3)) * t23 + t122 * t22 + t143) * t105 + (mrSges(6,2) * t165 + mrSges(6,3) * t168 + 0.2e1 * Ifges(6,1) * t22 - t94 * t3 + t97 * t4 + (t122 + t154) * t23 + (t105 * t76 - t97 * t13 - t94 * t14) * qJD(6)) * t36 + 0.2e1 * (-t109 * t54 + t110 * t53) * Ifges(5,4) + 0.2e1 * (-t109 * t32 + t110 * t31 - t42 * t54 - t43 * t53) * mrSges(5,3) + (m(4) * t163 + 0.2e1 * t115) * (pkin(2) * t138 - t142) + (((mrSges(3,2) * t169) - 0.2e1 * t72 * mrSges(4,3) + t145 * t131) * t98 + ((mrSges(3,1) * t169) + mrSges(4,1) * t163 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t131 - 0.2e1 * t145 * t96) * t96) * qJD(2); t146 * t7 + m(6) * (t17 * t28 - t39 * t7 + t40 * t6 - t153) + m(7) * (t37 * t7 - t153) + (t105 * t28 - t22 * t39 - t23 * t40 + t29 * t36) * mrSges(6,3) + (-t53 * t71 - t54 * t70 + (-t109 * t93 - t110 * t92) * qJD(3)) * mrSges(5,3) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t98 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t96 + (m(4) * t178 - t98 * mrSges(3,1) + t96 * mrSges(3,2) + t115) * pkin(7)) * qJD(2) + (-t25 * t135 - qJD(6) * t14 / 0.2e1 + t28 * t24 + t38 * t11 + m(7) * (t1 * t38 - t8 * t135 + t28 * t9) + (t77 * t147 - t69 / 0.2e1) * t36 + t121 * mrSges(7,3) - t119) * t97 + (-t24 * t135 + t13 * t147 - t28 * t25 - t38 * t10 + m(7) * (-t9 * t135 - t2 * t38 - t28 * t8) + (t78 * t147 + t68 / 0.2e1) * t36 + t120 * mrSges(7,3) - t118) * t94 - Ifges(5,5) * t54 + Ifges(5,6) * t53 + t37 * t5 + t32 * mrSges(5,2) + t29 * t19 - t31 * mrSges(5,1) + m(5) * (t31 * t70 + t32 * t71 + (-t42 * t92 + t43 * t93) * qJD(3)) + t176 * qJD(3) * t98 - t173; t37 * t164 + 0.2e1 * t179 + (t38 * t127 + t29 * t37) * t170 + (t28 * t40 - t29 * t39) * t171 + t101 + 0.2e1 * t146 * t29 - 0.2e1 * t128 * t28 + 0.2e1 * (mrSges(4,3) + t93 * mrSges(5,2) + t92 * mrSges(5,1) + m(5) * (-t70 * t92 + t71 * t93) + m(4) * qJ(3)) * qJD(3); t52 * t19 - t104 * t5 + t113 * t51 + t176 * t137 + (-t53 * t92 - t54 * t93) * mrSges(5,3) + (-t94 * t10 + t97 * t11 + (-t24 * t94 - t25 * t97) * qJD(6)) * t61 + m(7) * (t117 * t51 + (t1 * t97 - t2 * t94 + (-t8 * t97 - t9 * t94) * qJD(6)) * t61 + t116) + m(6) * (t17 * t51 + t6 * t61 + t116) + m(5) * (t31 * t93 + t32 * t92) + (-t104 * t22 + t105 * t51 - t23 * t61 + t36 * t52) * mrSges(6,3); t149 + t146 * t52 + (mrSges(6,2) - t128) * t51 + m(7) * (t28 * t125 + t38 * t126 + t37 * t52 - t152) + m(6) * (t28 * t61 - t39 * t52 + t40 * t51 - t152); 0.2e1 * m(6) * (t51 * t61 - t150) + 0.2e1 * m(7) * (t51 * t125 - t150); t97 * t10 + t94 * t11 + t113 * qJD(6) + m(7) * (t117 * qJD(6) + t1 * t94 + t2 * t97) + m(6) * t34 + m(5) * t44 + t123 + t124; 0; 0; 0; -pkin(5) * t5 + t172 * t7 + (t1 * mrSges(7,3) + t36 * t69 / 0.2e1 + (t36 * t160 - t8 * mrSges(7,3) + t14 / 0.2e1) * qJD(6) + (-m(7) * t121 - qJD(6) * t25 + t11) * pkin(9) + t119) * t97 + (-t2 * mrSges(7,3) + t68 * t161 + (t78 * t161 - t9 * mrSges(7,3) - t13 / 0.2e1) * qJD(6) + (-m(7) * t120 - qJD(6) * t24 - t10) * pkin(9) + t118) * t94 + t173; -t179 + (t37 + pkin(5)) * t66 + t174 * t127 + t172 * t29 - t101; -t51 * mrSges(6,2) + t174 * t126 + t172 * t52 - t149; 0; pkin(5) * t164 + t101; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,5) * t130 - t107 * Ifges(7,6) + t143; t82 - t114 * t28 + (t73 * t38 - t154) * qJD(6); (t61 * t134 - t51 * t97) * mrSges(7,2) + (-t61 * t133 - t51 * t94) * mrSges(7,1); -t66; -t82 + (t73 * pkin(9) + t154) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
