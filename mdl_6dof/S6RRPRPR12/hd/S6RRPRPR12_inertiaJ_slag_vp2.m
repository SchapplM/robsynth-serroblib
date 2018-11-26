% Calculate joint inertia matrix for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:08:19
% EndTime: 2018-11-23 17:08:20
% DurationCPUTime: 1.44s
% Computational Cost: add. (2535->341), mult. (5360->489), div. (0->0), fcn. (5702->10), ass. (0->128)
t175 = Ifges(4,1) + Ifges(3,3);
t174 = Ifges(5,3) + Ifges(6,3);
t128 = cos(qJ(4));
t120 = sin(pkin(11));
t122 = cos(pkin(11));
t125 = sin(qJ(4));
t80 = t120 * t125 - t122 * t128;
t81 = t120 * t128 + t122 * t125;
t173 = Ifges(5,5) * t128 - Ifges(6,5) * t80 - Ifges(6,6) * t81;
t130 = -pkin(2) - pkin(9);
t148 = -qJ(5) + t130;
t140 = t148 * t128;
t85 = t148 * t125;
t55 = t120 * t85 - t122 * t140;
t172 = t55 ^ 2;
t79 = t80 ^ 2;
t171 = -2 * mrSges(6,3);
t124 = sin(qJ(6));
t127 = cos(qJ(6));
t121 = sin(pkin(6));
t126 = sin(qJ(2));
t150 = t121 * t126;
t123 = cos(pkin(6));
t129 = cos(qJ(2));
t149 = t121 * t129;
t69 = -t123 * t125 - t128 * t149;
t70 = t123 * t128 - t125 * t149;
t42 = t120 * t69 + t122 * t70;
t27 = -t124 * t42 + t127 * t150;
t170 = t27 / 0.2e1;
t28 = t124 * t150 + t127 * t42;
t169 = t28 / 0.2e1;
t88 = Ifges(7,5) * t124 + Ifges(7,6) * t127;
t168 = t88 / 0.2e1;
t167 = t124 / 0.2e1;
t165 = -t127 / 0.2e1;
t164 = t127 / 0.2e1;
t163 = pkin(1) * t129;
t162 = t55 * t80;
t98 = pkin(8) * t150;
t71 = t123 * t163 - t98;
t161 = t71 * mrSges(3,1);
t72 = pkin(1) * t123 * t126 + pkin(8) * t149;
t160 = t72 * mrSges(3,2);
t11 = -mrSges(7,1) * t27 + mrSges(7,2) * t28;
t33 = mrSges(6,1) * t150 - mrSges(6,3) * t42;
t159 = t11 - t33;
t145 = -pkin(2) - t163;
t50 = pkin(3) * t150 + t98 + (-pkin(9) + t145) * t123;
t143 = -qJ(3) * t126 - pkin(1);
t59 = (t129 * t130 + t143) * t121;
t23 = -t125 * t59 + t128 * t50;
t17 = pkin(4) * t150 - qJ(5) * t70 + t23;
t24 = t125 * t50 + t128 * t59;
t20 = qJ(5) * t69 + t24;
t6 = t120 * t17 + t122 * t20;
t157 = Ifges(7,4) * t124;
t156 = Ifges(7,4) * t127;
t155 = t124 * t80;
t154 = t127 * t80;
t153 = t128 * mrSges(5,1);
t84 = mrSges(4,1) * t150 + mrSges(4,2) * t123;
t104 = pkin(4) * t120 + pkin(10);
t152 = t104 * t124;
t151 = t104 * t127;
t106 = pkin(4) * t125 + qJ(3);
t147 = t124 ^ 2 + t127 ^ 2;
t146 = t125 ^ 2 + t128 ^ 2;
t41 = t120 * t70 - t122 * t69;
t7 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t41;
t62 = -qJ(3) * t123 - t72;
t19 = t41 * mrSges(6,1) + mrSges(6,2) * t42;
t51 = t81 * mrSges(6,1) - mrSges(6,2) * t80;
t144 = m(5) * t146;
t142 = t146 * mrSges(5,3);
t141 = t104 * t147;
t139 = Ifges(3,5) * t150 + Ifges(3,6) * t149 + t123 * t175;
t58 = pkin(3) * t149 - t62;
t34 = -pkin(4) * t69 + t58;
t10 = pkin(5) * t41 - pkin(10) * t42 + t34;
t4 = pkin(10) * t150 + t6;
t1 = t10 * t127 - t124 * t4;
t2 = t10 * t124 + t127 * t4;
t138 = -t1 * t124 + t127 * t2;
t137 = -mrSges(7,1) * t124 - mrSges(7,2) * t127;
t5 = -t120 * t20 + t122 * t17;
t136 = t120 * t81 - t122 * t80;
t45 = pkin(5) * t81 + pkin(10) * t80 + t106;
t57 = t120 * t140 + t122 * t85;
t21 = -t124 * t57 + t127 * t45;
t22 = t124 * t45 + t127 * t57;
t135 = -t124 * t21 + t127 * t22;
t134 = t125 * t24 + t128 * t23;
t29 = -Ifges(7,5) * t154 + Ifges(7,6) * t155 + Ifges(7,3) * t81;
t133 = Ifges(5,5) * t70 + Ifges(6,5) * t42 + Ifges(5,6) * t69 - Ifges(6,6) * t41 + t150 * t174;
t131 = qJ(3) ^ 2;
t105 = -pkin(4) * t122 - pkin(5);
t92 = Ifges(5,1) * t128 - Ifges(5,4) * t125;
t91 = Ifges(7,1) * t124 + t156;
t90 = Ifges(5,4) * t128 - Ifges(5,2) * t125;
t89 = Ifges(7,2) * t127 + t157;
t87 = mrSges(5,1) * t125 + mrSges(5,2) * t128;
t86 = -mrSges(7,1) * t127 + mrSges(7,2) * t124;
t83 = -mrSges(4,1) * t149 - mrSges(4,3) * t123;
t78 = t81 ^ 2;
t64 = t123 * t145 + t98;
t63 = (-pkin(2) * t129 + t143) * t121;
t61 = mrSges(5,1) * t150 - mrSges(5,3) * t70;
t60 = -mrSges(5,2) * t150 + mrSges(5,3) * t69;
t53 = -Ifges(6,1) * t80 - Ifges(6,4) * t81;
t52 = -Ifges(6,4) * t80 - Ifges(6,2) * t81;
t49 = mrSges(7,1) * t81 + mrSges(7,3) * t154;
t48 = -mrSges(7,2) * t81 + mrSges(7,3) * t155;
t44 = -mrSges(5,1) * t69 + mrSges(5,2) * t70;
t43 = t137 * t80;
t36 = Ifges(5,1) * t70 + Ifges(5,4) * t69 + Ifges(5,5) * t150;
t35 = Ifges(5,4) * t70 + Ifges(5,2) * t69 + Ifges(5,6) * t150;
t32 = -mrSges(6,2) * t150 - mrSges(6,3) * t41;
t31 = Ifges(7,5) * t81 + (-Ifges(7,1) * t127 + t157) * t80;
t30 = Ifges(7,6) * t81 + (Ifges(7,2) * t124 - t156) * t80;
t16 = Ifges(6,1) * t42 - Ifges(6,4) * t41 + Ifges(6,5) * t150;
t15 = Ifges(6,4) * t42 - Ifges(6,2) * t41 + Ifges(6,6) * t150;
t13 = mrSges(7,1) * t41 - mrSges(7,3) * t28;
t12 = -mrSges(7,2) * t41 + mrSges(7,3) * t27;
t9 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t41;
t8 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t41;
t3 = -pkin(5) * t150 - t5;
t14 = [0.2e1 * t62 * t83 + 0.2e1 * t64 * t84 + t69 * t35 + t70 * t36 + 0.2e1 * t58 * t44 + 0.2e1 * t24 * t60 + 0.2e1 * t23 * t61 + t42 * t16 + 0.2e1 * t34 * t19 + t27 * t8 + t28 * t9 + 0.2e1 * t6 * t32 + 0.2e1 * t5 * t33 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + (m(3) * pkin(1) ^ 2 * t121 + (0.2e1 * t63 * mrSges(4,2) + 0.2e1 * t72 * mrSges(3,3) + (-(2 * Ifges(4,5)) + Ifges(3,6)) * t123 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + Ifges(4,3)) * t129) * t121) * t129 + (-0.2e1 * t71 * mrSges(3,3) - 0.2e1 * t63 * mrSges(4,3) + (-(2 * Ifges(4,4)) + Ifges(3,5)) * t123 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(4,2) + Ifges(3,1)) * t126 + 0.2e1 * (Ifges(3,4) + Ifges(4,6)) * t129) * t121 + t133) * t126) * t121 + m(3) * (t71 ^ 2 + t72 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2 + t58 ^ 2) + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t34 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + Ifges(2,3) + (t139 - 0.2e1 * t160 + 0.2e1 * t161) * t123 + (t7 - t15) * t41; (t130 * t60 - t35 / 0.2e1 - t24 * mrSges(5,3)) * t125 + (-t23 * mrSges(5,3) + t130 * t61 + t36 / 0.2e1) * t128 + (-t83 + t44) * qJ(3) + t106 * t19 + m(6) * (t106 * t34 - t5 * t55 + t57 * t6) + m(4) * (-pkin(2) * t64 - qJ(3) * t62) + m(7) * (t1 * t21 + t2 * t22 + t3 * t55) + (t7 / 0.2e1 - t15 / 0.2e1 - t6 * mrSges(6,3)) * t81 - pkin(2) * t84 + t58 * t87 + t69 * t90 / 0.2e1 + t70 * t92 / 0.2e1 + t57 * t32 - t62 * mrSges(4,3) + t64 * mrSges(4,2) + t3 * t43 + t2 * t48 + t1 * t49 + t34 * t51 + t42 * t53 / 0.2e1 + t21 * t13 + t22 * t12 + (-t52 / 0.2e1 + t29 / 0.2e1) * t41 + m(5) * (qJ(3) * t58 + t130 * t134) - t160 + t161 + (-t16 / 0.2e1 + t9 * t165 + t8 * t167 + t5 * mrSges(6,3)) * t80 + t139 + (-Ifges(4,4) * t126 - Ifges(4,5) * t129 + (-Ifges(5,6) * t125 + t173) * t126 / 0.2e1) * t121 + t31 * t169 + t30 * t170 + t159 * t55; -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t106 * t51 - t125 * t90 + t128 * t92 + 0.2e1 * t21 * t49 + 0.2e1 * t22 * t48 + 0.2e1 * t55 * t43 + (t171 * t57 + t29 - t52) * t81 + (t124 * t30 - t127 * t31 + t171 * t55 - t53) * t80 + m(7) * (t21 ^ 2 + t22 ^ 2 + t172) + m(6) * (t106 ^ 2 + t57 ^ 2 + t172) + m(5) * (t130 ^ 2 * t146 + t131) + m(4) * (pkin(2) ^ 2 + t131) + 0.2e1 * (t87 + mrSges(4,3)) * qJ(3) - 0.2e1 * t130 * t142 + t175; t125 * t60 + t128 * t61 + t159 * t80 + (t127 * t12 - t124 * t13 + t32) * t81 + m(7) * (t138 * t81 + t3 * t80) + m(6) * (-t5 * t80 + t6 * t81) + m(5) * t134 + m(4) * t64 + t84; -m(4) * pkin(2) - t79 * mrSges(6,3) + t80 * t43 + mrSges(4,2) - t142 + (-mrSges(6,3) * t81 - t124 * t49 + t127 * t48) * t81 + m(7) * (t135 * t81 + t162) + m(6) * (t57 * t81 + t162) + t130 * t144; m(4) + t144 + m(6) * (t78 + t79) + m(7) * (t147 * t78 + t79); t9 * t167 + t8 * t164 + t105 * t11 + t3 * t86 + t41 * t168 + t89 * t170 + t91 * t169 + t23 * mrSges(5,1) - t24 * mrSges(5,2) + t5 * mrSges(6,1) - t6 * mrSges(6,2) + t12 * t151 - t13 * t152 + m(7) * (t104 * t138 + t105 * t3) + (m(6) * (t120 * t6 + t122 * t5) + t122 * t33 + t120 * t32) * pkin(4) + t133 + t138 * mrSges(7,3); -t49 * t152 + t48 * t151 + t130 * t153 + t55 * t86 + t31 * t167 + t30 * t164 + t105 * t43 + t81 * t168 + m(7) * (t104 * t135 + t105 * t55) - t57 * mrSges(6,2) - t55 * mrSges(6,1) + (t165 * t91 + t167 * t89) * t80 + (-mrSges(5,2) * t130 - Ifges(5,6)) * t125 + t135 * mrSges(7,3) + (m(6) * (t120 * t57 - t122 * t55) - t136 * mrSges(6,3)) * pkin(4) + t173; t153 - t125 * mrSges(5,2) + (-mrSges(6,1) + t86) * t80 + (mrSges(7,3) * t147 - mrSges(6,2)) * t81 + m(7) * (t105 * t80 + t141 * t81) + m(6) * t136 * pkin(4); 0.2e1 * t105 * t86 + t124 * t91 + t127 * t89 + m(7) * (t104 ^ 2 * t147 + t105 ^ 2) + m(6) * (t120 ^ 2 + t122 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (mrSges(6,1) * t122 - mrSges(6,2) * t120) * pkin(4) + 0.2e1 * mrSges(7,3) * t141 + t174; t124 * t12 + t127 * t13 + m(7) * (t1 * t127 + t124 * t2) + m(6) * t34 + t19; t124 * t48 + t127 * t49 + m(7) * (t124 * t22 + t127 * t21) + m(6) * t106 + t51; 0; 0; m(7) * t147 + m(6); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t21 - mrSges(7,2) * t22 + t29; t137 * t81; t104 * t137 + t88; -t86; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
