% Calculate joint inertia matrix for
% S6RRPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta5]';
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
% Datum: 2019-03-09 10:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:58:53
% EndTime: 2019-03-09 09:58:56
% DurationCPUTime: 1.25s
% Computational Cost: add. (1197->266), mult. (2162->351), div. (0->0), fcn. (1972->6), ass. (0->96)
t131 = m(6) + m(7);
t136 = (mrSges(7,2) + mrSges(6,3));
t135 = pkin(3) + pkin(7);
t95 = sin(qJ(2));
t97 = cos(qJ(2));
t134 = t95 ^ 2 + t97 ^ 2;
t133 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t132 = m(6) * pkin(4);
t92 = sin(pkin(9));
t93 = cos(pkin(9));
t94 = sin(qJ(4));
t96 = cos(qJ(4));
t58 = t92 * t94 - t93 * t96;
t59 = t92 * t96 + t93 * t94;
t130 = -t58 ^ 2 - t59 ^ 2;
t129 = -m(4) * pkin(2) + mrSges(4,2);
t128 = -t94 / 0.2e1;
t98 = -pkin(2) - pkin(8);
t118 = qJ(5) * t97;
t112 = -qJ(3) * t95 - pkin(1);
t55 = t98 * t97 + t112;
t70 = t135 * t95;
t62 = t96 * t70;
t12 = pkin(4) * t95 + t62 + (-t55 + t118) * t94;
t19 = t96 * t55 + t94 * t70;
t16 = -t96 * t118 + t19;
t4 = t92 * t12 + t93 * t16;
t127 = Ifges(5,4) * t94;
t126 = Ifges(5,4) * t96;
t125 = t94 * t97;
t124 = t96 * mrSges(5,1);
t123 = t96 * t97;
t43 = -t93 * t123 + t92 * t125;
t30 = mrSges(7,2) * t43 + mrSges(7,3) * t95;
t31 = -mrSges(6,2) * t95 + mrSges(6,3) * t43;
t122 = t30 + t31;
t44 = t59 * t97;
t32 = mrSges(6,1) * t95 + mrSges(6,3) * t44;
t33 = -t95 * mrSges(7,1) - t44 * mrSges(7,2);
t121 = -t32 + t33;
t120 = t134 * pkin(7) ^ 2;
t71 = t135 * t97;
t119 = -t94 ^ 2 - t96 ^ 2;
t77 = t94 * pkin(4) + qJ(3);
t117 = -qJ(5) + t98;
t110 = t117 * t96;
t65 = t117 * t94;
t27 = -t93 * t110 + t65 * t92;
t29 = t92 * t110 + t93 * t65;
t116 = t27 ^ 2 + t29 ^ 2;
t46 = pkin(4) * t123 + t71;
t114 = m(5) * t119;
t113 = t119 * mrSges(5,3);
t15 = -t43 * mrSges(6,1) - t44 * mrSges(6,2);
t21 = t59 * mrSges(6,1) - t58 * mrSges(6,2);
t14 = -t43 * mrSges(7,1) + t44 * mrSges(7,3);
t20 = t59 * mrSges(7,1) + t58 * mrSges(7,3);
t109 = 2 * t136;
t107 = -t94 * mrSges(5,2) + t124;
t106 = -Ifges(5,5) * t94 - Ifges(5,6) * t96;
t3 = t12 * t93 - t16 * t92;
t18 = -t55 * t94 + t62;
t105 = t18 * t96 + t19 * t94;
t74 = pkin(4) * t92 + qJ(6);
t76 = -pkin(4) * t93 - pkin(5);
t104 = t58 * t76 + t59 * t74;
t103 = t58 * t93 - t59 * t92;
t102 = (-Ifges(7,4) - Ifges(6,5)) * t44 + (Ifges(6,6) - Ifges(7,6)) * t43 + t133 * t95;
t99 = qJ(3) ^ 2;
t81 = Ifges(5,5) * t96;
t69 = Ifges(5,1) * t96 - t127;
t68 = -Ifges(5,2) * t94 + t126;
t67 = mrSges(5,1) * t94 + mrSges(5,2) * t96;
t66 = -pkin(2) * t97 + t112;
t64 = -mrSges(5,2) * t95 - mrSges(5,3) * t123;
t63 = mrSges(5,1) * t95 + mrSges(5,3) * t125;
t53 = t107 * t97;
t52 = Ifges(7,4) * t58;
t51 = Ifges(6,5) * t58;
t50 = Ifges(6,6) * t59;
t49 = Ifges(7,6) * t59;
t42 = Ifges(5,5) * t95 + (-Ifges(5,1) * t94 - t126) * t97;
t41 = Ifges(5,6) * t95 + (-Ifges(5,2) * t96 - t127) * t97;
t25 = -Ifges(6,1) * t58 - Ifges(6,4) * t59;
t24 = -Ifges(7,1) * t58 + Ifges(7,5) * t59;
t23 = -Ifges(6,4) * t58 - Ifges(6,2) * t59;
t22 = -Ifges(7,5) * t58 + Ifges(7,3) * t59;
t17 = pkin(5) * t59 + qJ(6) * t58 + t77;
t10 = -Ifges(6,1) * t44 + Ifges(6,4) * t43 + Ifges(6,5) * t95;
t9 = -Ifges(7,1) * t44 + Ifges(7,4) * t95 - Ifges(7,5) * t43;
t8 = -Ifges(6,4) * t44 + Ifges(6,2) * t43 + Ifges(6,6) * t95;
t7 = -Ifges(7,5) * t44 + Ifges(7,6) * t95 - Ifges(7,3) * t43;
t5 = -pkin(5) * t43 + qJ(6) * t44 + t46;
t2 = -pkin(5) * t95 - t3;
t1 = qJ(6) * t95 + t4;
t6 = [0.2e1 * t1 * t30 + 0.2e1 * t5 * t14 + 0.2e1 * t46 * t15 + 0.2e1 * t18 * t63 + 0.2e1 * t19 * t64 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t32 + 0.2e1 * t4 * t31 + 0.2e1 * t71 * t53 + Ifges(2,3) - (t9 + t10) * t44 + (t8 - t7) * t43 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t66 * mrSges(4,2) - t96 * t41 - t94 * t42 + (Ifges(4,3) + Ifges(3,2)) * t97) * t97 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t66 * mrSges(4,3) + (Ifges(4,2) + Ifges(3,1)) * t95 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t106) * t97 + t102) * t95 + m(4) * (t66 ^ 2 + t120) + m(3) * (pkin(1) ^ 2 + t120) + m(5) * (t18 ^ 2 + t19 ^ 2 + t71 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t46 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t134; qJ(3) * t53 + t17 * t14 + t77 * t15 + t5 * t20 + t46 * t21 + t71 * t67 - (t24 / 0.2e1 + t25 / 0.2e1) * t44 + (t23 / 0.2e1 - t22 / 0.2e1) * t43 + t122 * t29 + t121 * t27 + (t42 / 0.2e1 + t98 * t63 - t18 * mrSges(5,3)) * t96 + (-t41 / 0.2e1 + t98 * t64 - t19 * mrSges(5,3)) * t94 + (Ifges(5,6) * t128 + t81 / 0.2e1 - t51 / 0.2e1 - t50 / 0.2e1 - t52 / 0.2e1 + t49 / 0.2e1 - Ifges(4,4) + Ifges(3,5) - pkin(2) * mrSges(4,1)) * t95 + m(5) * (qJ(3) * t71 + t105 * t98) + m(6) * (-t27 * t3 + t29 * t4 + t77 * t46) + m(7) * (t1 * t29 + t17 * t5 + t2 * t27) + (-t4 * mrSges(6,3) - t1 * mrSges(7,2) + t7 / 0.2e1 - t8 / 0.2e1) * t59 + (-t2 * mrSges(7,2) + t3 * mrSges(6,3) - t9 / 0.2e1 - t10 / 0.2e1) * t58 + (-Ifges(4,5) + Ifges(3,6) - t96 * t68 / 0.2e1 + qJ(3) * mrSges(4,1) + t69 * t128) * t97 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t97 + (-mrSges(3,1) + t129) * t95) * pkin(7); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t17 * t20 + 0.2e1 * t77 * t21 - t94 * t68 + t96 * t69 + Ifges(4,1) + Ifges(3,3) + (-t29 * t109 + t22 - t23) * t59 + (-t27 * t109 - t24 - t25) * t58 + m(6) * (t77 ^ 2 + t116) + m(7) * (t17 ^ 2 + t116) + m(5) * (-t119 * t98 ^ 2 + t99) + m(4) * (pkin(2) ^ 2 + t99) + 0.2e1 * (t67 + mrSges(4,3)) * qJ(3) + 0.2e1 * t98 * t113; t96 * t63 + t94 * t64 + (m(4) * pkin(7) + mrSges(4,1)) * t95 + t122 * t59 + t121 * t58 + m(7) * (t1 * t59 + t2 * t58) + m(6) * (-t3 * t58 + t4 * t59) + m(5) * t105; t113 - t98 * t114 + t131 * (t27 * t58 + t59 * t29) + t129 + t136 * t130; -t131 * t130 + m(4) - t114; m(7) * (t1 * t74 + t2 * t76) + t106 * t97 + (t92 * t31 + t93 * t32 + m(6) * (t3 * t93 + t4 * t92)) * pkin(4) + t74 * t30 + t76 * t33 + t18 * mrSges(5,1) - t19 * mrSges(5,2) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2) + t102; -t27 * mrSges(6,1) + t29 * mrSges(7,3) - t29 * mrSges(6,2) - t27 * mrSges(7,1) + m(7) * (t27 * t76 + t29 * t74) + t81 - t52 + t49 - t51 - t50 + t98 * t124 + (-mrSges(5,2) * t98 - Ifges(5,6)) * t94 - t104 * mrSges(7,2) + (m(6) * (-t27 * t93 + t29 * t92) + t103 * mrSges(6,3)) * pkin(4); (-mrSges(6,2) + mrSges(7,3)) * t59 + (-mrSges(6,1) - mrSges(7,1)) * t58 + m(7) * t104 - t103 * t132 + t107; -0.2e1 * t76 * mrSges(7,1) + 0.2e1 * t74 * mrSges(7,3) + m(7) * (t74 ^ 2 + t76 ^ 2) + (0.2e1 * mrSges(6,1) * t93 - 0.2e1 * mrSges(6,2) * t92 + (t92 ^ 2 + t93 ^ 2) * t132) * pkin(4) + t133; m(6) * t46 + m(7) * t5 + t14 + t15; m(6) * t77 + m(7) * t17 + t20 + t21; 0; 0; t131; m(7) * t2 + t33; m(7) * t27 - t58 * mrSges(7,2); m(7) * t58; m(7) * t76 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
