% Calculate joint inertia matrix for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRRR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:56
% EndTime: 2019-12-31 19:16:59
% DurationCPUTime: 1.04s
% Computational Cost: add. (2424->262), mult. (6477->403), div. (0->0), fcn. (7197->12), ass. (0->113)
t145 = 2 * pkin(9);
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t101 = cos(pkin(6));
t102 = cos(pkin(5));
t97 = sin(pkin(11));
t99 = sin(pkin(5));
t133 = t99 * t97;
t100 = cos(pkin(11));
t135 = pkin(1) * t102;
t83 = t100 * t135;
t48 = pkin(2) * t102 + t83 + (-pkin(8) * t101 - qJ(2)) * t133;
t127 = t101 * t48;
t123 = t99 * t100;
t118 = t101 * t123;
t98 = sin(pkin(6));
t111 = t102 * t98 + t118;
t63 = qJ(2) * t123 + t97 * t135;
t42 = t111 * pkin(8) + t63;
t53 = (-pkin(8) * t97 * t98 - pkin(2) * t100 - pkin(1)) * t99;
t22 = -t105 * t42 + (t53 * t98 + t127) * t108;
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t72 = -mrSges(6,1) * t106 + mrSges(6,2) * t103;
t144 = -m(6) * pkin(4) - mrSges(5,1) + t72;
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t125 = t105 * t98;
t64 = -t107 * t101 + t104 * t125;
t143 = t64 ^ 2;
t142 = 2 * mrSges(3,1);
t141 = 0.2e1 * t102;
t47 = t111 * t105 + t108 * t133;
t61 = t101 * t102 - t98 * t123;
t36 = t104 * t61 + t107 * t47;
t124 = t108 * t98;
t46 = -t102 * t124 + t105 * t133 - t108 * t118;
t26 = -t103 * t36 + t106 * t46;
t140 = t26 / 0.2e1;
t27 = t103 * t46 + t106 * t36;
t139 = t27 / 0.2e1;
t138 = -t103 / 0.2e1;
t137 = t103 / 0.2e1;
t136 = t106 / 0.2e1;
t134 = pkin(9) * t107;
t11 = -mrSges(6,1) * t26 + mrSges(6,2) * t27;
t29 = mrSges(5,1) * t46 - mrSges(5,3) * t36;
t132 = t11 - t29;
t31 = t101 * t53 - t48 * t98;
t17 = pkin(3) * t46 - pkin(9) * t47 + t31;
t23 = t105 * t127 + t108 * t42 + t53 * t125;
t20 = pkin(9) * t61 + t23;
t6 = t104 * t17 + t107 * t20;
t74 = Ifges(6,5) * t103 + Ifges(6,6) * t106;
t131 = Ifges(5,5) * t104 + Ifges(5,6) * t107;
t130 = t103 ^ 2 + t106 ^ 2;
t129 = Ifges(6,4) * t103;
t128 = Ifges(6,4) * t106;
t126 = t104 * t64;
t122 = t103 * t104;
t121 = t104 * t106;
t35 = t104 * t47 - t61 * t107;
t7 = Ifges(6,5) * t27 + Ifges(6,6) * t26 + Ifges(6,3) * t35;
t120 = Ifges(5,5) * t36 - Ifges(5,6) * t35 + Ifges(5,3) * t46;
t119 = Ifges(4,5) * t47 - Ifges(4,6) * t46 + Ifges(4,3) * t61;
t19 = -pkin(3) * t61 - t22;
t10 = pkin(4) * t35 - pkin(10) * t36 + t19;
t4 = pkin(10) * t46 + t6;
t1 = t10 * t106 - t103 * t4;
t2 = t10 * t103 + t106 * t4;
t117 = -t1 * t103 + t106 * t2;
t115 = mrSges(6,1) * t103 + mrSges(6,2) * t106;
t71 = -pkin(4) * t107 - pkin(10) * t104 - pkin(3);
t54 = -t103 * t134 + t106 * t71;
t55 = t103 * t71 + t106 * t134;
t113 = -t103 * t54 + t106 * t55;
t5 = -t104 * t20 + t107 * t17;
t66 = t101 * t104 + t107 * t125;
t112 = t107 * t66 + t126;
t57 = Ifges(6,5) * t121 - Ifges(6,6) * t122 - Ifges(6,3) * t107;
t110 = pkin(9) ^ 2;
t96 = t107 ^ 2;
t94 = t104 ^ 2;
t92 = t98 ^ 2;
t91 = t94 * t110;
t85 = t92 * t108 ^ 2;
t81 = mrSges(3,2) * t133;
t78 = Ifges(5,1) * t104 + Ifges(5,4) * t107;
t77 = Ifges(6,1) * t103 + t128;
t76 = Ifges(5,4) * t104 + Ifges(5,2) * t107;
t75 = Ifges(6,2) * t106 + t129;
t73 = -mrSges(5,1) * t107 + mrSges(5,2) * t104;
t70 = -mrSges(6,1) * t107 - mrSges(6,3) * t121;
t69 = mrSges(6,2) * t107 - mrSges(6,3) * t122;
t67 = t115 * t104;
t62 = -qJ(2) * t133 + t83;
t59 = -Ifges(6,5) * t107 + (Ifges(6,1) * t106 - t129) * t104;
t58 = -Ifges(6,6) * t107 + (-Ifges(6,2) * t103 + t128) * t104;
t51 = -t103 * t124 + t106 * t66;
t50 = -t103 * t66 - t106 * t124;
t38 = mrSges(4,1) * t61 - mrSges(4,3) * t47;
t37 = -mrSges(4,2) * t61 - mrSges(4,3) * t46;
t30 = mrSges(4,1) * t46 + mrSges(4,2) * t47;
t28 = -mrSges(5,2) * t46 - mrSges(5,3) * t35;
t21 = mrSges(5,1) * t35 + mrSges(5,2) * t36;
t15 = Ifges(5,1) * t36 - Ifges(5,4) * t35 + Ifges(5,5) * t46;
t14 = Ifges(5,4) * t36 - Ifges(5,2) * t35 + Ifges(5,6) * t46;
t13 = mrSges(6,1) * t35 - mrSges(6,3) * t27;
t12 = -mrSges(6,2) * t35 + mrSges(6,3) * t26;
t9 = Ifges(6,1) * t27 + Ifges(6,4) * t26 + Ifges(6,5) * t35;
t8 = Ifges(6,4) * t27 + Ifges(6,2) * t26 + Ifges(6,6) * t35;
t3 = -pkin(4) * t46 - t5;
t16 = [t61 * t119 + t47 * (Ifges(4,1) * t47 + Ifges(4,5) * t61) + t36 * t15 + 0.2e1 * t23 * t37 + 0.2e1 * t22 * t38 + t26 * t8 + t27 * t9 + 0.2e1 * t6 * t28 + 0.2e1 * t5 * t29 + 0.2e1 * t31 * t30 + 0.2e1 * t19 * t21 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + 0.2e1 * t1 * t13 + Ifges(2,3) + (-0.2e1 * t63 * mrSges(3,2) + Ifges(3,3) * t102 + t62 * t142) * t102 + (-0.2e1 * pkin(1) * t81 + (-0.2e1 * t62 * mrSges(3,3) + Ifges(3,1) * t133 + Ifges(3,5) * t141) * t97 + (0.2e1 * t63 * mrSges(3,3) + Ifges(3,6) * t141 + (0.2e1 * Ifges(3,4) * t97 + Ifges(3,2) * t100 + pkin(1) * t142) * t99) * t100) * t99 + (t7 - t14) * t35 + (-0.2e1 * Ifges(4,4) * t47 + Ifges(4,2) * t46 - Ifges(4,6) * t61 + t120) * t46 + m(3) * (pkin(1) ^ 2 * t99 ^ 2 + t62 ^ 2 + t63 ^ 2) + m(4) * (t22 ^ 2 + t23 ^ 2 + t31 ^ 2) + m(5) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2); t101 * t30 + t51 * t12 + t50 * t13 + t66 * t28 + t81 + (-m(3) * pkin(1) - mrSges(3,1) * t100) * t99 + t132 * t64 + (t105 * t37 + (-t21 + t38) * t108) * t98 + m(6) * (t1 * t50 + t2 * t51 + t3 * t64) + m(5) * (-t19 * t124 - t5 * t64 + t6 * t66) + m(4) * (t101 * t31 + (t105 * t23 + t108 * t22) * t98); m(3) + m(4) * (t105 ^ 2 * t92 + t101 ^ 2 + t85) + m(5) * (t66 ^ 2 + t143 + t85) + m(6) * (t50 ^ 2 + t51 ^ 2 + t143); (t15 / 0.2e1 + t9 * t136 + t8 * t138 - t5 * mrSges(5,3) + t132 * pkin(9)) * t104 + (t14 / 0.2e1 - t7 / 0.2e1 + pkin(9) * t28 + t6 * mrSges(5,3)) * t107 + m(6) * (pkin(9) * t104 * t3 + t1 * t54 + t2 * t55) + (-t76 / 0.2e1 + t57 / 0.2e1) * t35 + t119 + t36 * t78 / 0.2e1 + t46 * t131 / 0.2e1 + t3 * t67 + t2 * t69 + t1 * t70 + t19 * t73 + t54 * t13 + t55 * t12 + t58 * t140 + t59 * t139 - pkin(3) * t21 + t22 * mrSges(4,1) - t23 * mrSges(4,2) + m(5) * (-pkin(3) * t19 + (-t5 * t104 + t6 * t107) * pkin(9)); t50 * t70 + t51 * t69 + t64 * t67 + t112 * mrSges(5,3) + (-t105 * mrSges(4,2) + (mrSges(4,1) - t73) * t108) * t98 + m(5) * (pkin(3) * t124 + t112 * pkin(9)) + m(6) * (pkin(9) * t126 + t50 * t54 + t51 * t55); -0.2e1 * pkin(3) * t73 + 0.2e1 * t54 * t70 + 0.2e1 * t55 * t69 + Ifges(4,3) + (t76 - t57) * t107 + (t94 + t96) * mrSges(5,3) * t145 + m(5) * (pkin(3) ^ 2 + t110 * t96 + t91) + m(6) * (t54 ^ 2 + t55 ^ 2 + t91) + (-t103 * t58 + t106 * t59 + t67 * t145 + t78) * t104; t9 * t137 + t8 * t136 + t3 * t72 + t75 * t140 + t35 * t74 / 0.2e1 + t77 * t139 - t6 * mrSges(5,2) + t5 * mrSges(5,1) + t117 * mrSges(6,3) + t120 + (-m(6) * t3 - t11) * pkin(4) + (m(6) * t117 - t103 * t13 + t106 * t12) * pkin(10); -t66 * mrSges(5,2) + (m(6) * pkin(10) + mrSges(6,3)) * (-t103 * t50 + t106 * t51) + t144 * t64; t59 * t137 + t58 * t136 - pkin(4) * t67 + (m(6) * t113 - t103 * t70 + t106 * t69) * pkin(10) + (t144 * pkin(9) + t77 * t136 + t75 * t138) * t104 + (-pkin(9) * mrSges(5,2) - t74 / 0.2e1) * t107 + t113 * mrSges(6,3) + t131; Ifges(5,3) + t103 * t77 + t106 * t75 - 0.2e1 * pkin(4) * t72 + m(6) * (t130 * pkin(10) ^ 2 + pkin(4) ^ 2) + 0.2e1 * t130 * pkin(10) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t7; mrSges(6,1) * t50 - mrSges(6,2) * t51; mrSges(6,1) * t54 - mrSges(6,2) * t55 + t57; -t115 * pkin(10) + t74; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t16(1), t16(2), t16(4), t16(7), t16(11); t16(2), t16(3), t16(5), t16(8), t16(12); t16(4), t16(5), t16(6), t16(9), t16(13); t16(7), t16(8), t16(9), t16(10), t16(14); t16(11), t16(12), t16(13), t16(14), t16(15);];
Mq = res;
