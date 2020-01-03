% Calculate joint inertia matrix for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRPR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR12_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR12_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:06
% EndTime: 2019-12-31 21:37:10
% DurationCPUTime: 1.16s
% Computational Cost: add. (1814->308), mult. (4178->454), div. (0->0), fcn. (4435->10), ass. (0->124)
t160 = 2 * pkin(8);
t121 = sin(pkin(10));
t123 = cos(pkin(10));
t122 = sin(pkin(5));
t130 = cos(qJ(2));
t138 = t122 * t130;
t124 = cos(pkin(5));
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t127 = sin(qJ(2));
t139 = t122 * t127;
t77 = t124 * t126 + t129 * t139;
t48 = -t121 * t77 - t123 * t138;
t49 = -t121 * t138 + t123 * t77;
t76 = -t124 * t129 + t126 * t139;
t16 = Ifges(5,1) * t49 + Ifges(5,4) * t48 + Ifges(5,5) * t76;
t159 = t16 / 0.2e1;
t125 = sin(qJ(5));
t128 = cos(qJ(5));
t22 = -t125 * t49 + t128 * t48;
t158 = t22 / 0.2e1;
t23 = t125 * t48 + t128 * t49;
t157 = t23 / 0.2e1;
t141 = Ifges(5,4) * t123;
t66 = -Ifges(5,6) * t129 + (-Ifges(5,2) * t121 + t141) * t126;
t156 = t66 / 0.2e1;
t142 = Ifges(5,4) * t121;
t67 = -Ifges(5,5) * t129 + (Ifges(5,1) * t123 - t142) * t126;
t155 = t67 / 0.2e1;
t87 = t121 * t128 + t123 * t125;
t71 = t87 * t126;
t154 = -t71 / 0.2e1;
t86 = -t121 * t125 + t123 * t128;
t72 = t86 * t126;
t153 = t72 / 0.2e1;
t152 = t86 / 0.2e1;
t151 = t87 / 0.2e1;
t97 = Ifges(5,1) * t121 + t141;
t150 = t97 / 0.2e1;
t149 = pkin(1) * t130;
t148 = pkin(8) * t126;
t147 = pkin(8) * t129;
t105 = pkin(7) * t139;
t79 = t124 * t149 - t105;
t146 = t79 * mrSges(3,1);
t80 = t124 * t127 * pkin(1) + pkin(7) * t138;
t145 = t80 * mrSges(3,2);
t144 = pkin(9) + qJ(4);
t68 = t105 + (-pkin(2) - t149) * t124;
t28 = pkin(3) * t76 - qJ(4) * t77 + t68;
t69 = pkin(8) * t124 + t80;
t70 = (-pkin(2) * t130 - pkin(8) * t127 - pkin(1)) * t122;
t34 = t126 * t70 + t129 * t69;
t29 = -qJ(4) * t138 + t34;
t10 = t121 * t28 + t123 * t29;
t143 = -Ifges(4,5) * t77 + Ifges(4,6) * t76;
t45 = Ifges(6,5) * t87 + Ifges(6,6) * t86;
t91 = -pkin(3) * t129 - qJ(4) * t126 - pkin(2);
t61 = t121 * t91 + t123 * t147;
t140 = t121 * t126;
t137 = t123 * t126;
t78 = mrSges(5,1) * t140 + mrSges(5,2) * t137;
t136 = Ifges(4,5) * t126 + Ifges(4,6) * t129;
t135 = t121 ^ 2 + t123 ^ 2;
t4 = Ifges(6,5) * t23 + Ifges(6,6) * t22 + Ifges(6,3) * t76;
t134 = Ifges(5,5) * t121 / 0.2e1 + Ifges(5,6) * t123 / 0.2e1 + t45 / 0.2e1;
t133 = Ifges(3,5) * t139 + Ifges(3,6) * t138 + Ifges(3,3) * t124;
t24 = -t48 * mrSges(5,1) + t49 * mrSges(5,2);
t8 = -t22 * mrSges(6,1) + t23 * mrSges(6,2);
t40 = t71 * mrSges(6,1) + t72 * mrSges(6,2);
t44 = -t86 * mrSges(6,1) + t87 * mrSges(6,2);
t9 = -t121 * t29 + t123 * t28;
t33 = -t126 * t69 + t129 * t70;
t93 = -t123 * mrSges(5,1) + t121 * mrSges(5,2);
t35 = Ifges(6,5) * t72 - Ifges(6,6) * t71 - Ifges(6,3) * t129;
t30 = pkin(3) * t138 - t33;
t132 = pkin(8) ^ 2;
t120 = t129 ^ 2;
t119 = t126 ^ 2;
t116 = t119 * t132;
t111 = -pkin(4) * t123 - pkin(3);
t100 = Ifges(4,1) * t126 + Ifges(4,4) * t129;
t99 = Ifges(4,4) * t126 + Ifges(4,2) * t129;
t98 = -mrSges(4,1) * t129 + mrSges(4,2) * t126;
t96 = Ifges(5,2) * t123 + t142;
t94 = t144 * t123;
t92 = t144 * t121;
t90 = (pkin(4) * t121 + pkin(8)) * t126;
t89 = -mrSges(5,1) * t129 - mrSges(5,3) * t137;
t88 = mrSges(5,2) * t129 - mrSges(5,3) * t140;
t85 = t123 * t91;
t65 = -Ifges(5,3) * t129 + (Ifges(5,5) * t123 - Ifges(5,6) * t121) * t126;
t60 = -t121 * t147 + t85;
t56 = -mrSges(6,1) * t129 - mrSges(6,3) * t72;
t55 = mrSges(6,2) * t129 - mrSges(6,3) * t71;
t54 = -mrSges(4,1) * t138 - mrSges(4,3) * t77;
t53 = mrSges(4,2) * t138 - mrSges(4,3) * t76;
t52 = -t125 * t92 + t128 * t94;
t51 = -t125 * t94 - t128 * t92;
t50 = -pkin(9) * t140 + t61;
t47 = Ifges(6,1) * t87 + Ifges(6,4) * t86;
t46 = Ifges(6,4) * t87 + Ifges(6,2) * t86;
t42 = -pkin(9) * t137 + t85 + (-pkin(8) * t121 - pkin(4)) * t129;
t41 = mrSges(4,1) * t76 + mrSges(4,2) * t77;
t39 = Ifges(4,1) * t77 - Ifges(4,4) * t76 - Ifges(4,5) * t138;
t38 = Ifges(4,4) * t77 - Ifges(4,2) * t76 - Ifges(4,6) * t138;
t37 = Ifges(6,1) * t72 - Ifges(6,4) * t71 - Ifges(6,5) * t129;
t36 = Ifges(6,4) * t72 - Ifges(6,2) * t71 - Ifges(6,6) * t129;
t32 = mrSges(5,1) * t76 - mrSges(5,3) * t49;
t31 = -mrSges(5,2) * t76 + mrSges(5,3) * t48;
t18 = t125 * t42 + t128 * t50;
t17 = -t125 * t50 + t128 * t42;
t15 = Ifges(5,4) * t49 + Ifges(5,2) * t48 + Ifges(5,6) * t76;
t14 = Ifges(5,5) * t49 + Ifges(5,6) * t48 + Ifges(5,3) * t76;
t13 = -pkin(4) * t48 + t30;
t12 = mrSges(6,1) * t76 - mrSges(6,3) * t23;
t11 = -mrSges(6,2) * t76 + mrSges(6,3) * t22;
t7 = pkin(9) * t48 + t10;
t6 = Ifges(6,1) * t23 + Ifges(6,4) * t22 + Ifges(6,5) * t76;
t5 = Ifges(6,4) * t23 + Ifges(6,2) * t22 + Ifges(6,6) * t76;
t3 = pkin(4) * t76 - pkin(9) * t49 + t9;
t2 = t125 * t3 + t128 * t7;
t1 = -t125 * t7 + t128 * t3;
t19 = [0.2e1 * t1 * t12 + 0.2e1 * t10 * t31 + 0.2e1 * t2 * t11 + 0.2e1 * t13 * t8 + t48 * t15 + t49 * t16 + t22 * t5 + t23 * t6 + 0.2e1 * t30 * t24 + 0.2e1 * t9 * t32 + 0.2e1 * t33 * t54 + 0.2e1 * t34 * t53 + t77 * t39 + 0.2e1 * t68 * t41 + Ifges(2,3) + (t14 + t4 - t38) * t76 + (t133 - 0.2e1 * t145 + 0.2e1 * t146) * t124 + ((-0.2e1 * t79 * mrSges(3,3) + Ifges(3,5) * t124 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t127) * t122) * t127 + (0.2e1 * t80 * mrSges(3,3) + Ifges(3,6) * t124 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t127 + (Ifges(3,2) + Ifges(4,3)) * t130) * t122 + t143) * t130) * t122 + m(3) * (pkin(1) ^ 2 * t122 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2 + t68 ^ 2) + m(5) * (t10 ^ 2 + t30 ^ 2 + t9 ^ 2) + m(6) * (t1 ^ 2 + t13 ^ 2 + t2 ^ 2); -t136 * t138 / 0.2e1 - t145 + t146 + m(4) * (-pkin(2) * t68 + (-t33 * t126 + t34 * t129) * pkin(8)) + (t38 / 0.2e1 - t14 / 0.2e1 - t4 / 0.2e1 + pkin(8) * t53 + t34 * mrSges(4,3)) * t129 + m(6) * (t1 * t17 + t13 * t90 + t18 * t2) + (t39 / 0.2e1 + t123 * t159 - t121 * t15 / 0.2e1 - t33 * mrSges(4,3) + (-t54 + t24) * pkin(8)) * t126 + m(5) * (t10 * t61 + t30 * t148 + t60 * t9) + (-t99 / 0.2e1 + t65 / 0.2e1 + t35 / 0.2e1) * t76 + t10 * t88 + t9 * t89 + t90 * t8 + t68 * t98 + t77 * t100 / 0.2e1 + t30 * t78 + t2 * t55 + t1 * t56 + t60 * t32 + t61 * t31 - pkin(2) * t41 + t13 * t40 + t17 * t12 + t18 * t11 + t6 * t153 + t5 * t154 + t49 * t155 + t48 * t156 + t37 * t157 + t36 * t158 + t133; -0.2e1 * pkin(2) * t98 + 0.2e1 * t17 * t56 + 0.2e1 * t18 * t55 - t71 * t36 + t72 * t37 + 0.2e1 * t90 * t40 + 0.2e1 * t60 * t89 + 0.2e1 * t61 * t88 + Ifges(3,3) + (t119 + t120) * mrSges(4,3) * t160 + (-t35 + t99 - t65) * t129 + m(6) * (t17 ^ 2 + t18 ^ 2 + t90 ^ 2) + m(5) * (t60 ^ 2 + t61 ^ 2 + t116) + m(4) * (pkin(2) ^ 2 + t120 * t132 + t116) + (-t121 * t66 + t123 * t67 + t78 * t160 + t100) * t126; (t15 / 0.2e1 + qJ(4) * t31 + t10 * mrSges(5,3)) * t123 + (-t9 * mrSges(5,3) - qJ(4) * t32 + t159) * t121 - Ifges(4,3) * t138 + t134 * t76 + m(6) * (t1 * t51 + t111 * t13 + t2 * t52) - t143 + (-t1 * t87 + t2 * t86) * mrSges(6,3) + m(5) * (-pkin(3) * t30 + (t10 * t123 - t121 * t9) * qJ(4)) + t111 * t8 + t5 * t152 + t6 * t151 + t30 * t93 + t48 * t96 / 0.2e1 + t49 * t150 + t51 * t12 + t52 * t11 + t13 * t44 + t46 * t158 + t47 * t157 + t33 * mrSges(4,1) - t34 * mrSges(4,2) - pkin(3) * t24; t111 * t40 + t36 * t152 + t37 * t151 + t90 * t44 - pkin(3) * t78 + t46 * t154 + t47 * t153 + t52 * t55 + t51 * t56 + (-mrSges(4,1) + t93) * t148 + (-t17 * t87 + t18 * t86) * mrSges(6,3) + (t61 * mrSges(5,3) + qJ(4) * t88 + t126 * t150 + t156) * t123 + (t155 - t126 * t96 / 0.2e1 - qJ(4) * t89 - t60 * mrSges(5,3)) * t121 + m(5) * (-pkin(3) * t148 + (-t121 * t60 + t123 * t61) * qJ(4)) + m(6) * (t111 * t90 + t17 * t51 + t18 * t52) + (-pkin(8) * mrSges(4,2) - t134) * t129 + t136; -0.2e1 * pkin(3) * t93 + 0.2e1 * t111 * t44 + t121 * t97 + t123 * t96 + t86 * t46 + t87 * t47 + Ifges(4,3) + m(6) * (t111 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t135 * qJ(4) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t51 * t87 + t52 * t86) * mrSges(6,3) + 0.2e1 * t135 * qJ(4) * mrSges(5,3); m(5) * t30 + m(6) * t13 + t24 + t8; m(5) * t148 + m(6) * t90 + t40 + t78; -m(5) * pkin(3) + m(6) * t111 + t44 + t93; m(5) + m(6); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t4; mrSges(6,1) * t17 - mrSges(6,2) * t18 + t35; mrSges(6,1) * t51 - mrSges(6,2) * t52 + t45; 0; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t19(1), t19(2), t19(4), t19(7), t19(11); t19(2), t19(3), t19(5), t19(8), t19(12); t19(4), t19(5), t19(6), t19(9), t19(13); t19(7), t19(8), t19(9), t19(10), t19(14); t19(11), t19(12), t19(13), t19(14), t19(15);];
Mq = res;
