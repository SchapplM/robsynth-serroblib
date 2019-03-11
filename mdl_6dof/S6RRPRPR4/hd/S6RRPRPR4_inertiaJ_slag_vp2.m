% Calculate joint inertia matrix for
% S6RRPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:23
% EndTime: 2019-03-09 10:21:27
% DurationCPUTime: 1.64s
% Computational Cost: add. (3677->341), mult. (8662->505), div. (0->0), fcn. (9731->12), ass. (0->131)
t181 = Ifges(3,3) + Ifges(4,3);
t180 = Ifges(5,3) + Ifges(6,3);
t129 = sin(pkin(11));
t132 = cos(pkin(11));
t130 = sin(pkin(6));
t139 = cos(qJ(2));
t153 = t130 * t139;
t136 = sin(qJ(2));
t154 = t130 * t136;
t81 = t129 * t154 - t132 * t153;
t82 = (t129 * t139 + t132 * t136) * t130;
t179 = Ifges(4,5) * t82 - Ifges(4,6) * t81;
t135 = sin(qJ(4));
t138 = cos(qJ(4));
t128 = sin(pkin(12));
t131 = cos(pkin(12));
t93 = t128 * t135 - t131 * t138;
t95 = t128 * t138 + t131 * t135;
t178 = Ifges(5,5) * t135 + Ifges(6,5) * t95 + Ifges(5,6) * t138 - Ifges(6,6) * t93;
t115 = pkin(2) * t129 + pkin(9);
t152 = qJ(5) + t115;
t147 = t152 * t135;
t91 = t152 * t138;
t55 = t128 * t91 + t131 * t147;
t177 = t55 ^ 2;
t176 = t93 ^ 2;
t175 = 0.2e1 * t55;
t96 = (-pkin(2) * t139 - pkin(1)) * t130;
t174 = 0.2e1 * t96;
t134 = sin(qJ(6));
t137 = cos(qJ(6));
t133 = cos(pkin(6));
t68 = t133 * t138 - t135 * t82;
t69 = t133 * t135 + t138 * t82;
t39 = t128 * t68 + t131 * t69;
t25 = -t134 * t39 + t137 * t81;
t173 = t25 / 0.2e1;
t26 = t134 * t81 + t137 * t39;
t172 = t26 / 0.2e1;
t100 = Ifges(7,5) * t134 + Ifges(7,6) * t137;
t170 = t100 / 0.2e1;
t169 = -t134 / 0.2e1;
t168 = t134 / 0.2e1;
t167 = t137 / 0.2e1;
t166 = pkin(1) * t133;
t165 = t55 * t93;
t110 = t139 * t166;
t85 = -pkin(8) * t154 + t110;
t164 = t85 * mrSges(3,1);
t86 = pkin(8) * t153 + t136 * t166;
t163 = t86 * mrSges(3,2);
t11 = -mrSges(7,1) * t25 + mrSges(7,2) * t26;
t29 = mrSges(6,1) * t81 - mrSges(6,3) * t39;
t162 = t11 - t29;
t70 = pkin(2) * t133 + t110 + (-pkin(8) - qJ(3)) * t154;
t75 = qJ(3) * t153 + t86;
t45 = t129 * t70 + t132 * t75;
t43 = pkin(9) * t133 + t45;
t48 = pkin(3) * t81 - pkin(9) * t82 + t96;
t21 = -t135 * t43 + t138 * t48;
t13 = pkin(4) * t81 - qJ(5) * t69 + t21;
t22 = t135 * t48 + t138 * t43;
t19 = qJ(5) * t68 + t22;
t6 = t128 * t13 + t131 * t19;
t160 = Ifges(7,4) * t134;
t159 = Ifges(7,4) * t137;
t158 = t134 * t95;
t157 = t137 * t95;
t114 = pkin(4) * t128 + pkin(10);
t156 = t114 * t134;
t155 = t114 * t137;
t150 = t134 ^ 2 + t137 ^ 2;
t149 = t135 ^ 2 + t138 ^ 2;
t38 = t128 * t69 - t131 * t68;
t7 = Ifges(7,5) * t26 + Ifges(7,6) * t25 + Ifges(7,3) * t38;
t117 = -pkin(2) * t132 - pkin(3);
t20 = t38 * mrSges(6,1) + t39 * mrSges(6,2);
t87 = t95 * mrSges(6,2);
t64 = t93 * mrSges(6,1) + t87;
t44 = -t129 * t75 + t132 * t70;
t148 = t114 * t150;
t42 = -pkin(3) * t133 - t44;
t27 = -pkin(4) * t68 + t42;
t10 = pkin(5) * t38 - pkin(10) * t39 + t27;
t4 = pkin(10) * t81 + t6;
t1 = t10 * t137 - t134 * t4;
t2 = t10 * t134 + t137 * t4;
t146 = -t1 * t134 + t137 * t2;
t99 = -t138 * mrSges(5,1) + t135 * mrSges(5,2);
t145 = mrSges(7,1) * t134 + mrSges(7,2) * t137;
t5 = -t128 * t19 + t13 * t131;
t97 = -pkin(4) * t138 + t117;
t54 = pkin(5) * t93 - pkin(10) * t95 + t97;
t57 = -t128 * t147 + t131 * t91;
t30 = -t134 * t57 + t137 * t54;
t31 = t134 * t54 + t137 * t57;
t144 = -t134 * t30 + t137 * t31;
t51 = Ifges(7,5) * t157 - Ifges(7,6) * t158 + Ifges(7,3) * t93;
t143 = Ifges(5,5) * t69 + Ifges(6,5) * t39 + Ifges(5,6) * t68 - Ifges(6,6) * t38 + t180 * t81;
t142 = Ifges(3,5) * t154 + Ifges(3,6) * t153 + t181 * t133 + t179;
t116 = -pkin(4) * t131 - pkin(5);
t104 = Ifges(5,1) * t135 + Ifges(5,4) * t138;
t103 = Ifges(7,1) * t134 + t159;
t102 = Ifges(5,4) * t135 + Ifges(5,2) * t138;
t101 = Ifges(7,2) * t137 + t160;
t98 = -mrSges(7,1) * t137 + mrSges(7,2) * t134;
t92 = t95 ^ 2;
t76 = t82 * mrSges(4,2);
t72 = mrSges(4,1) * t133 - mrSges(4,3) * t82;
t71 = -mrSges(4,2) * t133 - mrSges(4,3) * t81;
t66 = Ifges(6,1) * t95 - Ifges(6,4) * t93;
t65 = Ifges(6,4) * t95 - Ifges(6,2) * t93;
t60 = mrSges(7,1) * t93 - mrSges(7,3) * t157;
t59 = -mrSges(7,2) * t93 - mrSges(7,3) * t158;
t58 = t145 * t95;
t53 = Ifges(7,5) * t93 + (Ifges(7,1) * t137 - t160) * t95;
t52 = Ifges(7,6) * t93 + (-Ifges(7,2) * t134 + t159) * t95;
t50 = mrSges(5,1) * t81 - mrSges(5,3) * t69;
t49 = -mrSges(5,2) * t81 + mrSges(5,3) * t68;
t41 = -mrSges(5,1) * t68 + mrSges(5,2) * t69;
t33 = Ifges(5,1) * t69 + Ifges(5,4) * t68 + Ifges(5,5) * t81;
t32 = Ifges(5,4) * t69 + Ifges(5,2) * t68 + Ifges(5,6) * t81;
t28 = -mrSges(6,2) * t81 - mrSges(6,3) * t38;
t17 = Ifges(6,1) * t39 - Ifges(6,4) * t38 + Ifges(6,5) * t81;
t16 = Ifges(6,4) * t39 - Ifges(6,2) * t38 + Ifges(6,6) * t81;
t15 = mrSges(7,1) * t38 - mrSges(7,3) * t26;
t14 = -mrSges(7,2) * t38 + mrSges(7,3) * t25;
t9 = Ifges(7,1) * t26 + Ifges(7,4) * t25 + Ifges(7,5) * t38;
t8 = Ifges(7,4) * t26 + Ifges(7,2) * t25 + Ifges(7,6) * t38;
t3 = -pkin(5) * t81 - t5;
t12 = [(t7 - t16) * t38 + (t142 - 0.2e1 * t163 + 0.2e1 * t164 + t179) * t133 + ((Ifges(3,5) * t136 + Ifges(3,6) * t139) * t133 + 0.2e1 * (-t136 * t85 + t139 * t86) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t139 + mrSges(3,2) * t136) + t136 * (Ifges(3,1) * t136 + Ifges(3,4) * t139) + t139 * (Ifges(3,4) * t136 + Ifges(3,2) * t139)) * t130) * t130 + t76 * t174 + Ifges(4,1) * t82 ^ 2 + m(3) * (t85 ^ 2 + t86 ^ 2) + (mrSges(4,1) * t174 - 0.2e1 * Ifges(4,4) * t82 + Ifges(4,2) * t81 + t143) * t81 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2 + t96 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2 + t42 ^ 2) + m(6) * (t27 ^ 2 + t5 ^ 2 + t6 ^ 2) + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t14 + 0.2e1 * t1 * t15 + t25 * t8 + t26 * t9 + 0.2e1 * t27 * t20 + 0.2e1 * t6 * t28 + 0.2e1 * t5 * t29 + t39 * t17 + 0.2e1 * t42 * t41 + 0.2e1 * t22 * t49 + Ifges(2,3) + 0.2e1 * t21 * t50 + t68 * t32 + t69 * t33 + 0.2e1 * t45 * t71 + 0.2e1 * t44 * t72; (t22 * mrSges(5,3) + t115 * t49 + t32 / 0.2e1) * t138 + (-t21 * mrSges(5,3) - t115 * t50 + t33 / 0.2e1) * t135 + t178 * t81 / 0.2e1 + m(7) * (t1 * t30 + t2 * t31 + t3 * t55) + m(6) * (t27 * t97 - t5 * t55 + t57 * t6) + (-t6 * mrSges(6,3) + t7 / 0.2e1 - t16 / 0.2e1) * t93 + (t51 / 0.2e1 - t65 / 0.2e1) * t38 + (t129 * t71 + t132 * t72 + m(4) * (t129 * t45 + t132 * t44)) * pkin(2) + t142 + t164 - t163 + t53 * t172 + t52 * t173 + (-t5 * mrSges(6,3) + t8 * t169 + t9 * t167 + t17 / 0.2e1) * t95 + t162 * t55 + m(5) * (t117 * t42 + (-t135 * t21 + t138 * t22) * t115) + t30 * t15 + t31 * t14 + t44 * mrSges(4,1) - t45 * mrSges(4,2) + t57 * t28 + t3 * t58 + t2 * t59 + t1 * t60 + t27 * t64 + t39 * t66 / 0.2e1 + t97 * t20 + t42 * t99 + t68 * t102 / 0.2e1 + t69 * t104 / 0.2e1 + t117 * t41; t138 * t102 + t135 * t104 + 0.2e1 * t117 * t99 + 0.2e1 * t30 * t60 + 0.2e1 * t31 * t59 + t58 * t175 + 0.2e1 * t97 * t64 + (-0.2e1 * mrSges(6,3) * t57 + t51 - t65) * t93 + (mrSges(6,3) * t175 - t134 * t52 + t137 * t53 + t66) * t95 + m(7) * (t30 ^ 2 + t31 ^ 2 + t177) + m(6) * (t57 ^ 2 + t97 ^ 2 + t177) + m(5) * (t149 * t115 ^ 2 + t117 ^ 2) + m(4) * (t129 ^ 2 + t132 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t132 - mrSges(4,2) * t129) * pkin(2) + 0.2e1 * t149 * t115 * mrSges(5,3) + t181; t81 * mrSges(4,1) + t135 * t49 + t138 * t50 + t76 + t162 * t93 + (-t134 * t15 + t137 * t14 + t28) * t95 + m(7) * (t146 * t95 + t3 * t93) + m(6) * (-t5 * t93 + t6 * t95) + m(5) * (t135 * t22 + t138 * t21) + m(4) * t96; t93 * t58 + (-t134 * t60 + t137 * t59) * t95 + m(7) * (t144 * t95 + t165) + m(6) * (t57 * t95 + t165); m(4) + m(5) * t149 + m(6) * (t92 + t176) + m(7) * (t150 * t92 + t176); t5 * mrSges(6,1) - t6 * mrSges(6,2) + t143 + t146 * mrSges(7,3) - t15 * t156 + t14 * t155 + t21 * mrSges(5,1) - t22 * mrSges(5,2) + t3 * t98 + t38 * t170 + t101 * t173 + t103 * t172 + t116 * t11 + t9 * t168 + t8 * t167 + m(7) * (t146 * t114 + t116 * t3) + (t128 * t28 + t131 * t29 + m(6) * (t128 * t6 + t131 * t5)) * pkin(4); t55 * t98 + t93 * t170 + m(7) * (t144 * t114 + t116 * t55) + t53 * t168 + t52 * t167 + t116 * t58 - t55 * mrSges(6,1) - t57 * mrSges(6,2) + t59 * t155 - t60 * t156 + (t101 * t169 + t103 * t167) * t95 + (-mrSges(5,1) * t135 - mrSges(5,2) * t138) * t115 + t144 * mrSges(7,3) + (m(6) * (t128 * t57 - t131 * t55) + (-t128 * t93 - t131 * t95) * mrSges(6,3)) * pkin(4) + t178; -t87 + (-mrSges(6,1) + t98) * t93 + t150 * t95 * mrSges(7,3) + m(7) * (t116 * t93 + t95 * t148) + m(6) * (t128 * t95 - t131 * t93) * pkin(4) - t99; t137 * t101 + t134 * t103 + 0.2e1 * t116 * t98 + m(7) * (t150 * t114 ^ 2 + t116 ^ 2) + m(6) * (t128 ^ 2 + t131 ^ 2) * pkin(4) ^ 2 + 0.2e1 * (mrSges(6,1) * t131 - mrSges(6,2) * t128) * pkin(4) + 0.2e1 * mrSges(7,3) * t148 + t180; t134 * t14 + t137 * t15 + m(7) * (t1 * t137 + t134 * t2) + m(6) * t27 + t20; t134 * t59 + t137 * t60 + m(7) * (t134 * t31 + t137 * t30) + m(6) * t97 + t64; 0; 0; m(7) * t150 + m(6); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t30 - mrSges(7,2) * t31 + t51; -t58; -t145 * t114 + t100; -t98; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
