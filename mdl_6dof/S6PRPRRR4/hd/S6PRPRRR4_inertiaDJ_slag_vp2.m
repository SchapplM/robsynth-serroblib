% Calculate time derivative of joint inertia matrix for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR4_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR4_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:36:05
% EndTime: 2019-03-08 20:36:11
% DurationCPUTime: 2.74s
% Computational Cost: add. (4758->381), mult. (11578->580), div. (0->0), fcn. (11787->12), ass. (0->161)
t129 = cos(qJ(5));
t120 = sin(pkin(12));
t122 = cos(pkin(12));
t126 = sin(qJ(4));
t130 = cos(qJ(4));
t97 = t120 * t126 - t122 * t130;
t93 = t97 * qJD(4);
t166 = t129 * t93;
t98 = t120 * t130 + t122 * t126;
t94 = t98 * qJD(4);
t196 = -Ifges(6,5) * t166 + Ifges(6,3) * t94;
t124 = sin(qJ(6));
t125 = sin(qJ(5));
t128 = cos(qJ(6));
t137 = t124 * t125 - t128 * t129;
t64 = t137 * t98;
t141 = mrSges(6,1) * t125 + mrSges(6,2) * t129;
t101 = t141 * qJD(5);
t190 = qJD(5) + qJD(6);
t77 = t190 * t137;
t100 = t124 * t129 + t125 * t128;
t78 = t190 * t100;
t40 = mrSges(7,1) * t78 - t77 * mrSges(7,2);
t195 = t101 + t40;
t194 = m(6) * pkin(9) + mrSges(6,3);
t174 = pkin(8) + qJ(3);
t106 = t174 * t120;
t107 = t174 * t122;
t193 = -t106 * t130 - t107 * t126;
t108 = -mrSges(6,1) * t129 + mrSges(6,2) * t125;
t192 = -m(6) * pkin(4) - mrSges(5,1) + t108;
t123 = cos(pkin(6));
t121 = sin(pkin(6));
t127 = sin(qJ(2));
t161 = t121 * t127;
t91 = -t120 * t161 + t122 * t123;
t92 = t120 * t123 + t122 * t161;
t191 = -t120 * t91 + t122 * t92;
t189 = 2 * m(7);
t188 = -2 * mrSges(5,3);
t83 = -t106 * t126 + t107 * t130;
t56 = qJD(3) * t98 + qJD(4) * t83;
t187 = 0.2e1 * t56;
t186 = -0.2e1 * t193;
t183 = m(7) * pkin(5);
t181 = -t98 / 0.2e1;
t180 = -pkin(10) - pkin(9);
t172 = Ifges(6,4) * t125;
t109 = Ifges(6,2) * t129 + t172;
t179 = -t109 / 0.2e1;
t178 = Ifges(6,5) * t94;
t177 = Ifges(6,6) * t94;
t176 = t56 * t193;
t131 = cos(qJ(2));
t157 = qJD(2) * t131;
t146 = t121 * t157;
t59 = t126 * t91 + t130 * t92;
t39 = qJD(4) * t59 + t146 * t98;
t58 = t126 * t92 - t130 * t91;
t26 = t58 * t39;
t175 = t97 * Ifges(6,6);
t173 = -Ifges(7,5) * t77 - Ifges(7,6) * t78;
t114 = -pkin(3) * t122 - pkin(2);
t70 = pkin(4) * t97 - pkin(9) * t98 + t114;
t76 = t129 * t83;
t37 = t125 * t70 + t76;
t171 = Ifges(6,4) * t129;
t170 = Ifges(6,6) * t125;
t167 = t125 * t98;
t165 = t129 * t98;
t162 = t121 ^ 2 * t127;
t160 = t121 * t131;
t159 = t120 ^ 2 + t122 ^ 2;
t158 = qJD(2) * t121;
t156 = qJD(5) * t125;
t155 = qJD(5) * t129;
t154 = qJD(6) * t124;
t153 = qJD(6) * t128;
t22 = t137 * t93 - t78 * t98;
t23 = t100 * t93 + t190 * t64;
t152 = Ifges(7,5) * t22 + Ifges(7,6) * t23 + Ifges(7,3) * t94;
t151 = pkin(5) * t156;
t147 = t127 * t158;
t38 = -qJD(4) * t58 - t146 * t97;
t50 = -t125 * t59 - t129 * t160;
t18 = qJD(5) * t50 + t125 * t147 + t129 * t38;
t133 = t125 * t160 - t129 * t59;
t19 = qJD(5) * t133 - t125 * t38 + t129 * t147;
t24 = t124 * t133 + t128 * t50;
t5 = qJD(6) * t24 + t124 * t19 + t128 * t18;
t25 = t124 * t50 - t128 * t133;
t6 = -qJD(6) * t25 - t124 * t18 + t128 * t19;
t150 = mrSges(7,1) * t6 - t5 * mrSges(7,2);
t149 = t98 * t156;
t148 = qJD(5) * t180;
t67 = t94 * mrSges(5,1) - mrSges(5,2) * t93;
t145 = -(2 * Ifges(5,4)) - t170;
t55 = -qJD(3) * t97 + qJD(4) * t193;
t69 = pkin(4) * t94 + pkin(9) * t93;
t144 = -t125 * t55 + t129 * t69;
t36 = -t125 * t83 + t129 * t70;
t142 = -t193 * t39 + t56 * t58;
t140 = Ifges(6,1) * t129 - t172;
t139 = -Ifges(6,2) * t125 + t171;
t27 = pkin(5) * t97 - pkin(10) * t165 + t36;
t32 = -pkin(10) * t167 + t37;
t12 = -t124 * t32 + t128 * t27;
t13 = t124 * t27 + t128 * t32;
t104 = t125 * t148;
t105 = t129 * t148;
t111 = t180 * t125;
t112 = t180 * t129;
t84 = t111 * t128 + t112 * t124;
t44 = qJD(6) * t84 + t104 * t128 + t105 * t124;
t85 = t111 * t124 - t112 * t128;
t45 = -qJD(6) * t85 - t104 * t124 + t105 * t128;
t138 = mrSges(7,1) * t45 - t44 * mrSges(7,2) + t173;
t10 = pkin(10) * t166 + pkin(5) * t94 + (-t76 + (pkin(10) * t98 - t70) * t125) * qJD(5) + t144;
t135 = -t125 * t93 + t155 * t98;
t14 = t125 * t69 + t129 * t55 + t155 * t70 - t156 * t83;
t11 = -pkin(10) * t135 + t14;
t2 = qJD(6) * t12 + t10 * t124 + t11 * t128;
t3 = -qJD(6) * t13 + t10 * t128 - t11 * t124;
t136 = mrSges(7,1) * t3 - mrSges(7,2) * t2 + t152;
t134 = t149 + t166;
t116 = Ifges(6,5) * t155;
t115 = -pkin(5) * t129 - pkin(4);
t110 = Ifges(6,1) * t125 + t171;
t103 = t140 * qJD(5);
t102 = t139 * qJD(5);
t95 = (-mrSges(7,1) * t124 - mrSges(7,2) * t128) * qJD(6) * pkin(5);
t81 = Ifges(7,1) * t100 - Ifges(7,4) * t137;
t80 = Ifges(7,4) * t100 - Ifges(7,2) * t137;
t79 = mrSges(7,1) * t137 + mrSges(7,2) * t100;
t72 = mrSges(6,1) * t97 - mrSges(6,3) * t165;
t71 = -mrSges(6,2) * t97 - mrSges(6,3) * t167;
t68 = t141 * t98;
t63 = t100 * t98;
t57 = pkin(5) * t167 - t193;
t53 = Ifges(6,5) * t97 + t140 * t98;
t52 = t139 * t98 + t175;
t49 = mrSges(7,1) * t97 + mrSges(7,3) * t64;
t48 = -mrSges(7,2) * t97 - mrSges(7,3) * t63;
t47 = -mrSges(6,2) * t94 - mrSges(6,3) * t135;
t46 = mrSges(6,1) * t94 + mrSges(6,3) * t134;
t42 = -Ifges(7,1) * t77 - Ifges(7,4) * t78;
t41 = -Ifges(7,4) * t77 - Ifges(7,2) * t78;
t35 = mrSges(6,1) * t135 - mrSges(6,2) * t134;
t34 = mrSges(7,1) * t63 - mrSges(7,2) * t64;
t33 = pkin(5) * t135 + t56;
t31 = -Ifges(7,1) * t64 - Ifges(7,4) * t63 + Ifges(7,5) * t97;
t30 = -Ifges(7,4) * t64 - Ifges(7,2) * t63 + Ifges(7,6) * t97;
t29 = -Ifges(6,1) * t134 - Ifges(6,4) * t135 + t178;
t28 = -Ifges(6,4) * t134 - Ifges(6,2) * t135 + t177;
t17 = -mrSges(7,2) * t94 + mrSges(7,3) * t23;
t16 = mrSges(7,1) * t94 - mrSges(7,3) * t22;
t15 = -qJD(5) * t37 + t144;
t9 = -mrSges(7,1) * t23 + mrSges(7,2) * t22;
t8 = Ifges(7,1) * t22 + Ifges(7,4) * t23 + Ifges(7,5) * t94;
t7 = Ifges(7,4) * t22 + Ifges(7,2) * t23 + t94 * Ifges(7,6);
t1 = [0.2e1 * m(7) * (t24 * t6 + t25 * t5 + t26) + 0.2e1 * m(6) * (-t133 * t18 + t19 * t50 + t26) + 0.2e1 * m(4) * (t121 * t191 - t162) * t157 + 0.2e1 * (-t157 * t162 + t38 * t59 + t26) * m(5); t24 * t16 + t25 * t17 + t18 * t71 + t19 * t72 + t50 * t46 - t133 * t47 + t5 * t48 + t6 * t49 + (t35 + t9) * t58 + (t34 + t68) * t39 + (-t38 * t97 + t39 * t98 - t58 * t93 - t59 * t94) * mrSges(5,3) + (-t131 * t67 + ((mrSges(4,3) * t159 - mrSges(3,2)) * t131 + (-mrSges(4,1) * t122 + mrSges(5,1) * t97 + mrSges(4,2) * t120 + mrSges(5,2) * t98 - mrSges(3,1)) * t127) * qJD(2)) * t121 + m(4) * (t191 * qJD(3) + (qJ(3) * t131 * t159 - pkin(2) * t127) * t158) + m(5) * (t114 * t147 + t38 * t83 + t55 * t59 + t142) + m(6) * (-t133 * t14 + t15 * t50 + t18 * t37 + t19 * t36 + t142) + m(7) * (t12 * t6 + t13 * t5 + t2 * t25 + t24 * t3 + t33 * t58 + t39 * t57); -(mrSges(5,3) * t186 - t125 * t52 + t129 * t53) * t93 + (-Ifges(7,5) * t64 - Ifges(7,6) * t63 + t83 * t188) * t94 + (mrSges(5,3) * t187 - 0.2e1 * Ifges(5,1) * t93 - t125 * t28 + t129 * t29 + (Ifges(6,5) * t129 + t145) * t94 + (t97 * (-Ifges(6,5) * t125 - Ifges(6,6) * t129) - t129 * t52 - t125 * t53) * qJD(5)) * t98 + (t55 * t188 - t145 * t93 + ((2 * Ifges(5,2)) + Ifges(6,3) + Ifges(7,3)) * t94 + t152 + t196) * t97 + 0.2e1 * t114 * t67 - t63 * t7 - t64 * t8 + 0.2e1 * t14 * t71 + 0.2e1 * t15 * t72 + 0.2e1 * t57 * t9 + 0.2e1 * t36 * t46 + 0.2e1 * t37 * t47 + 0.2e1 * t2 * t48 + 0.2e1 * t3 * t49 + t23 * t30 + t22 * t31 + 0.2e1 * t33 * t34 + 0.2e1 * t12 * t16 + 0.2e1 * t13 * t17 + 0.2e1 * m(6) * (t14 * t37 + t15 * t36 - t176) + 0.2e1 * m(5) * (t55 * t83 - t176) + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * t159 * qJD(3) + t35 * t186 + t68 * t187 + (t12 * t3 + t13 * t2 + t33 * t57) * t189; (m(4) + m(5)) * t147 + m(7) * (t100 * t5 - t137 * t6 - t24 * t78 - t25 * t77) + m(6) * (t125 * t18 + t129 * t19 + (-t125 * t50 - t129 * t133) * qJD(5)); t100 * t17 + t125 * t47 + t129 * t46 - t137 * t16 - t77 * t48 - t78 * t49 + (-t125 * t72 + t129 * t71) * qJD(5) + m(7) * (t100 * t2 - t12 * t78 - t13 * t77 - t137 * t3) + m(6) * (t125 * t14 + t129 * t15 + (-t125 * t36 + t129 * t37) * qJD(5)) + t67; (-t100 * t77 + t137 * t78) * t189; -t38 * mrSges(5,2) + t195 * t58 + m(7) * (t151 * t58 + t24 * t45 + t25 * t44 + t5 * t85 + t6 * t84) + (-t100 * t6 - t137 * t5 + t24 * t77 - t25 * t78) * mrSges(7,3) + (m(7) * t115 + t192 + t79) * t39 + t194 * (-t19 * t125 + t18 * t129 + (t125 * t133 - t129 * t50) * qJD(5)); (t177 / 0.2e1 + t28 / 0.2e1 + t98 * t103 / 0.2e1 - t93 * t110 / 0.2e1 + t14 * mrSges(6,3) + (t53 / 0.2e1 + t98 * t179 - t36 * mrSges(6,3)) * qJD(5) + (m(6) * (-t36 * qJD(5) + t14) + t47 - qJD(5) * t72) * pkin(9)) * t129 + (t178 / 0.2e1 + t29 / 0.2e1 + t102 * t181 - t93 * t179 - t15 * mrSges(6,3) + (-m(6) * t15 - t46) * pkin(9) + (-t175 / 0.2e1 - t52 / 0.2e1 + t110 * t181 - pkin(9) * t71 + pkin(5) * t34 + t57 * t183 - t194 * t37) * qJD(5)) * t125 + t192 * t56 + t115 * t9 + t100 * t8 / 0.2e1 - Ifges(5,5) * t93 - Ifges(5,6) * t94 + t84 * t16 + t85 * t17 - t77 * t31 / 0.2e1 - t78 * t30 / 0.2e1 + t33 * t79 + t23 * t80 / 0.2e1 + t22 * t81 / 0.2e1 - t63 * t41 / 0.2e1 - t64 * t42 / 0.2e1 - t55 * mrSges(5,2) + t57 * t40 + t44 * t48 + t45 * t49 - pkin(4) * t35 + (-t3 * t100 + t12 * t77 - t13 * t78 - t137 * t2) * mrSges(7,3) + t94 * (Ifges(7,5) * t100 - Ifges(7,6) * t137) / 0.2e1 - t137 * t7 / 0.2e1 + (t173 + t116) * t97 / 0.2e1 - t193 * t101 + m(7) * (t115 * t33 + t12 * t45 + t13 * t44 + t2 * t85 + t3 * t84); m(7) * (t100 * t44 - t137 * t45 - t77 * t85 - t78 * t84); -t77 * t81 + t100 * t42 + 0.2e1 * t79 * t151 + 0.2e1 * t115 * t40 - t78 * t80 - t137 * t41 + (t115 * t151 + t44 * t85 + t45 * t84) * t189 - 0.2e1 * pkin(4) * t101 - t109 * t156 + t125 * t103 + (qJD(5) * t110 + t102) * t129 + 0.2e1 * (-t100 * t45 - t137 * t44 + t77 * t84 - t78 * t85) * mrSges(7,3); t19 * mrSges(6,1) - t18 * mrSges(6,2) + (t124 * t5 + t128 * t6 + (-t124 * t24 + t128 * t25) * qJD(6)) * t183 + t150; -Ifges(6,5) * t149 + t15 * mrSges(6,1) - t14 * mrSges(6,2) - t135 * Ifges(6,6) + (m(7) * (-t12 * t154 + t124 * t2 + t128 * t3 + t13 * t153) + t48 * t153 + t124 * t17 - t49 * t154 + t128 * t16) * pkin(5) + t136 + t196; (-t124 * t77 - t128 * t78 + (t100 * t128 + t124 * t137) * qJD(6)) * t183 - t195; t116 + (pkin(9) * t108 - t170) * qJD(5) + (m(7) * (t124 * t44 + t128 * t45 + (-t124 * t84 + t128 * t85) * qJD(6)) + (-t124 * t78 + t128 * t77 + (t100 * t124 - t128 * t137) * qJD(6)) * mrSges(7,3)) * pkin(5) + t138; 0.2e1 * t95; t150; t136; -t40; t138; t95; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
