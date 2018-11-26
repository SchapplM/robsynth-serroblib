% Calculate joint inertia matrix for
% S6RPRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRPRR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:08:12
% EndTime: 2018-11-23 16:08:13
% DurationCPUTime: 1.72s
% Computational Cost: add. (5551->336), mult. (15062->507), div. (0->0), fcn. (17335->14), ass. (0->137)
t180 = Ifges(4,3) + Ifges(5,3);
t121 = sin(pkin(13));
t109 = pkin(3) * t121 + pkin(10);
t179 = 0.2e1 * t109;
t129 = sin(qJ(6));
t132 = cos(qJ(6));
t130 = sin(qJ(5));
t133 = cos(qJ(5));
t125 = cos(pkin(13));
t122 = sin(pkin(12));
t124 = sin(pkin(6));
t126 = cos(pkin(12));
t128 = cos(pkin(6));
t131 = sin(qJ(3));
t127 = cos(pkin(7));
t134 = cos(qJ(3));
t151 = t127 * t134;
t123 = sin(pkin(7));
t154 = t123 * t134;
t65 = t128 * t154 + (-t122 * t131 + t126 * t151) * t124;
t152 = t127 * t131;
t155 = t123 * t131;
t66 = t128 * t155 + (t122 * t134 + t126 * t152) * t124;
t51 = t121 * t65 + t125 * t66;
t153 = t124 * t126;
t86 = -t123 * t153 + t127 * t128;
t43 = t130 * t86 + t133 * t51;
t50 = t121 * t66 - t125 * t65;
t26 = -t129 * t43 + t132 * t50;
t42 = t130 * t51 - t133 * t86;
t17 = -mrSges(7,2) * t42 + mrSges(7,3) * t26;
t27 = t129 * t50 + t132 * t43;
t18 = mrSges(7,1) * t42 - mrSges(7,3) * t27;
t178 = -t129 * t18 + t132 * t17;
t177 = Ifges(4,5) * t66 + Ifges(5,5) * t51 + Ifges(4,6) * t65;
t93 = -mrSges(7,1) * t132 + mrSges(7,2) * t129;
t176 = -m(7) * pkin(5) - mrSges(6,1) + t93;
t80 = (t121 * t134 + t125 * t131) * t123;
t70 = -t133 * t127 + t130 * t80;
t175 = t70 ^ 2;
t78 = t121 * t155 - t125 * t154;
t174 = t78 ^ 2;
t173 = 2 * mrSges(3,1);
t172 = 0.2e1 * t128;
t171 = t26 / 0.2e1;
t170 = t27 / 0.2e1;
t169 = -t129 / 0.2e1;
t168 = t129 / 0.2e1;
t167 = t132 / 0.2e1;
t166 = pkin(1) * t128;
t165 = pkin(5) * t133;
t11 = -mrSges(7,1) * t26 + mrSges(7,2) * t27;
t29 = mrSges(6,1) * t50 - mrSges(6,3) * t43;
t164 = t11 - t29;
t88 = qJ(2) * t153 + t122 * t166;
t62 = (t123 * t128 + t127 * t153) * pkin(9) + t88;
t105 = t126 * t166;
t156 = t122 * t124;
t67 = pkin(2) * t128 + t105 + (-pkin(9) * t127 - qJ(2)) * t156;
t76 = (-pkin(9) * t122 * t123 - pkin(2) * t126 - pkin(1)) * t124;
t36 = -t131 * t62 + t67 * t151 + t76 * t154;
t31 = pkin(3) * t86 - qJ(4) * t66 + t36;
t37 = t134 * t62 + t67 * t152 + t76 * t155;
t35 = qJ(4) * t65 + t37;
t16 = t121 * t31 + t125 * t35;
t14 = pkin(10) * t86 + t16;
t53 = -t123 * t67 + t127 * t76;
t41 = -pkin(3) * t65 + t53;
t22 = pkin(4) * t50 - pkin(10) * t51 + t41;
t6 = t130 * t22 + t133 * t14;
t163 = Ifges(7,4) * t129;
t162 = Ifges(7,4) * t132;
t159 = t133 * t70;
t158 = t109 * t130;
t157 = t109 * t133;
t150 = t129 * t130;
t149 = t130 * t132;
t95 = Ifges(7,5) * t129 + Ifges(7,6) * t132;
t148 = Ifges(6,5) * t130 + Ifges(6,6) * t133;
t147 = t129 ^ 2 + t132 ^ 2;
t118 = t130 ^ 2;
t120 = t133 ^ 2;
t146 = t118 + t120;
t8 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t42;
t145 = Ifges(6,5) * t43 - Ifges(6,6) * t42 + Ifges(6,3) * t50;
t110 = -pkin(3) * t125 - pkin(4);
t32 = t50 * mrSges(5,1) + t51 * mrSges(5,2);
t15 = -t121 * t35 + t125 * t31;
t144 = t147 * t130;
t4 = pkin(11) * t50 + t6;
t13 = -pkin(4) * t86 - t15;
t7 = pkin(5) * t42 - pkin(11) * t43 + t13;
t1 = -t129 * t4 + t132 * t7;
t2 = t129 * t7 + t132 * t4;
t143 = -t1 * t129 + t132 * t2;
t94 = -t133 * mrSges(6,1) + t130 * mrSges(6,2);
t142 = mrSges(7,1) * t129 + mrSges(7,2) * t132;
t72 = t127 * t130 + t133 * t80;
t54 = -t129 * t72 + t132 * t78;
t55 = t129 * t78 + t132 * t72;
t141 = -t129 * t54 + t132 * t55;
t90 = -pkin(11) * t130 + t110 - t165;
t68 = -t129 * t157 + t132 * t90;
t69 = t129 * t90 + t132 * t157;
t140 = -t129 * t68 + t132 * t69;
t91 = mrSges(7,2) * t133 - mrSges(7,3) * t150;
t92 = -mrSges(7,1) * t133 - mrSges(7,3) * t149;
t139 = -t129 * t92 + t132 * t91;
t5 = -t130 * t14 + t133 * t22;
t138 = t130 * t70 + t133 * t72;
t137 = -Ifges(5,6) * t50 + t180 * t86 + t177;
t81 = Ifges(7,5) * t149 - Ifges(7,6) * t150 - Ifges(7,3) * t133;
t116 = t127 ^ 2;
t107 = t109 ^ 2;
t103 = mrSges(3,2) * t156;
t100 = t118 * t107;
t99 = Ifges(6,1) * t130 + Ifges(6,4) * t133;
t98 = Ifges(7,1) * t129 + t162;
t97 = Ifges(6,4) * t130 + Ifges(6,2) * t133;
t96 = Ifges(7,2) * t132 + t163;
t89 = t142 * t130;
t87 = -qJ(2) * t156 + t105;
t83 = -Ifges(7,5) * t133 + (Ifges(7,1) * t132 - t163) * t130;
t82 = -Ifges(7,6) * t133 + (-Ifges(7,2) * t129 + t162) * t130;
t57 = mrSges(4,1) * t86 - mrSges(4,3) * t66;
t56 = -mrSges(4,2) * t86 + mrSges(4,3) * t65;
t52 = -mrSges(4,1) * t65 + mrSges(4,2) * t66;
t45 = mrSges(5,1) * t86 - mrSges(5,3) * t51;
t44 = -mrSges(5,2) * t86 - mrSges(5,3) * t50;
t28 = -mrSges(6,2) * t50 - mrSges(6,3) * t42;
t23 = mrSges(6,1) * t42 + mrSges(6,2) * t43;
t20 = Ifges(6,1) * t43 - Ifges(6,4) * t42 + Ifges(6,5) * t50;
t19 = Ifges(6,4) * t43 - Ifges(6,2) * t42 + Ifges(6,6) * t50;
t10 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t42;
t9 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t42;
t3 = -pkin(5) * t50 - t5;
t12 = [(-0.2e1 * t88 * mrSges(3,2) + Ifges(3,3) * t128 + t173 * t87) * t128 + (-0.2e1 * pkin(1) * t103 + (-0.2e1 * t87 * mrSges(3,3) + Ifges(3,1) * t156 + Ifges(3,5) * t172) * t122 + (0.2e1 * t88 * mrSges(3,3) + Ifges(3,6) * t172 + (0.2e1 * Ifges(3,4) * t122 + Ifges(3,2) * t126 + pkin(1) * t173) * t124) * t126) * t124 + m(6) * (t13 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(3) * (pkin(1) ^ 2 * t124 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t36 ^ 2 + t37 ^ 2 + t53 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t41 ^ 2) + (-0.2e1 * Ifges(5,4) * t51 + Ifges(5,2) * t50 - Ifges(5,6) * t86 + t145) * t50 + (t8 - t19) * t42 + (0.2e1 * Ifges(4,4) * t66 + Ifges(4,2) * t65) * t65 + Ifges(4,1) * t66 ^ 2 + (t137 + t177) * t86 + Ifges(5,1) * t51 ^ 2 + Ifges(2,3) + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t13 * t23 + t26 * t9 + t27 * t10 + 0.2e1 * t6 * t28 + 0.2e1 * t5 * t29 + 0.2e1 * t41 * t32 + t43 * t20 + 0.2e1 * t16 * t44 + 0.2e1 * t15 * t45 + 0.2e1 * t53 * t52 + 0.2e1 * t37 * t56 + 0.2e1 * t36 * t57; t55 * t17 + t54 * t18 + t72 * t28 + t80 * t44 + t103 + (t23 - t45) * t78 + t164 * t70 + (t32 + t52) * t127 + (-m(3) * pkin(1) - t126 * mrSges(3,1)) * t124 + (t131 * t56 + t134 * t57) * t123 + m(7) * (t1 * t54 + t2 * t55 + t3 * t70) + m(6) * (t13 * t78 - t5 * t70 + t6 * t72) + m(5) * (t127 * t41 - t15 * t78 + t16 * t80) + m(4) * (t127 * t53 + (t131 * t37 + t134 * t36) * t123); m(3) + m(6) * (t72 ^ 2 + t174 + t175) + m(7) * (t54 ^ 2 + t55 ^ 2 + t175) + m(5) * (t80 ^ 2 + t116 + t174) + m(4) * (t116 + (t131 ^ 2 + t134 ^ 2) * t123 ^ 2); (t81 / 0.2e1 - t97 / 0.2e1) * t42 + t137 + (t121 * t44 + t125 * t45 + m(5) * (t121 * t16 + t125 * t15)) * pkin(3) + m(6) * (t110 * t13 + (-t5 * t130 + t6 * t133) * t109) + m(7) * (t1 * t68 + t158 * t3 + t2 * t69) + (-t5 * mrSges(6,3) + t9 * t169 + t10 * t167 + t20 / 0.2e1 + t164 * t109) * t130 + (t6 * mrSges(6,3) + t109 * t28 - t8 / 0.2e1 + t19 / 0.2e1) * t133 + t15 * mrSges(5,1) - t16 * mrSges(5,2) + t36 * mrSges(4,1) - t37 * mrSges(4,2) + t68 * t18 + t69 * t17 + t82 * t171 + t83 * t170 + t3 * t89 + t2 * t91 + t1 * t92 + t13 * t94 + t43 * t99 / 0.2e1 + t110 * t23 + t50 * t148 / 0.2e1; -t80 * mrSges(5,2) + t54 * t92 + t55 * t91 + t70 * t89 + (-mrSges(5,1) + t94) * t78 + (mrSges(4,1) * t134 - mrSges(4,2) * t131) * t123 + t138 * mrSges(6,3) + m(6) * (t109 * t138 + t110 * t78) + m(7) * (t158 * t70 + t54 * t68 + t55 * t69) + m(5) * (t121 * t80 - t125 * t78) * pkin(3); 0.2e1 * t110 * t94 + 0.2e1 * t68 * t92 + 0.2e1 * t69 * t91 + (-t81 + t97) * t133 + m(7) * (t68 ^ 2 + t69 ^ 2 + t100) + m(6) * (t107 * t120 + t110 ^ 2 + t100) + m(5) * (t121 ^ 2 + t125 ^ 2) * pkin(3) ^ 2 + (-t129 * t82 + t132 * t83 + t179 * t89 + t99) * t130 + 0.2e1 * (mrSges(5,1) * t125 - mrSges(5,2) * t121) * pkin(3) + t146 * mrSges(6,3) * t179 + t180; -t164 * t133 + (t28 + t178) * t130 + m(7) * (t130 * t143 - t133 * t3) + m(6) * (t130 * t6 + t133 * t5) + m(5) * t41 + t32; m(5) * t127 + m(6) * (t130 * t72 - t159) + m(7) * (t130 * t141 - t159); -t133 * t89 + (m(7) * (t140 - t157) + t139) * t130; m(5) + m(7) * (t118 * t147 + t120) + m(6) * t146; t9 * t167 + t3 * t93 + t98 * t170 + t96 * t171 + t42 * t95 / 0.2e1 + t10 * t168 - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t143 * mrSges(7,3) + t145 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t143 + t178) * pkin(11); -t72 * mrSges(6,2) + (m(7) * pkin(11) + mrSges(7,3)) * t141 + t176 * t70; t83 * t168 + t82 * t167 - pkin(5) * t89 + (m(7) * t140 + t139) * pkin(11) + (t109 * t176 + t98 * t167 + t96 * t169) * t130 + (-t95 / 0.2e1 - t109 * mrSges(6,2)) * t133 + t140 * mrSges(7,3) + t148; -t133 * t93 + m(7) * (pkin(11) * t144 + t165) + mrSges(7,3) * t144 - t94; Ifges(6,3) - 0.2e1 * pkin(5) * t93 + m(7) * (pkin(11) ^ 2 * t147 + pkin(5) ^ 2) + t129 * t98 + t132 * t96 + 0.2e1 * t147 * pkin(11) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t8; mrSges(7,1) * t54 - mrSges(7,2) * t55; mrSges(7,1) * t68 - mrSges(7,2) * t69 + t81; -t89; -pkin(11) * t142 + t95; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
