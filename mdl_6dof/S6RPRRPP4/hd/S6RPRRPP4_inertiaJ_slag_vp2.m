% Calculate joint inertia matrix for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:47
% EndTime: 2019-03-09 04:38:50
% DurationCPUTime: 1.11s
% Computational Cost: add. (1769->249), mult. (3432->352), div. (0->0), fcn. (3748->8), ass. (0->89)
t138 = Ifges(7,4) + Ifges(6,5);
t137 = Ifges(6,6) - Ifges(7,6);
t136 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t135 = m(6) * pkin(4);
t134 = m(6) + m(7);
t101 = cos(qJ(4));
t95 = sin(pkin(10));
t97 = cos(pkin(10));
t99 = sin(qJ(4));
t107 = t97 * t101 - t95 * t99;
t73 = t101 * t95 + t97 * t99;
t133 = Ifges(5,5) * t99 + Ifges(5,6) * t101 + t137 * t107 + t138 * t73;
t100 = sin(qJ(3));
t128 = cos(qJ(3));
t124 = pkin(7) + qJ(2);
t96 = sin(pkin(9));
t76 = t124 * t96;
t98 = cos(pkin(9));
t77 = t124 * t98;
t49 = t100 * t77 + t128 * t76;
t132 = t49 ^ 2;
t92 = t98 ^ 2;
t131 = 0.2e1 * t49;
t86 = -pkin(2) * t98 - pkin(1);
t130 = 0.2e1 * t86;
t74 = t100 * t98 + t128 * t96;
t126 = t74 * t99;
t72 = t100 * t96 - t128 * t98;
t39 = pkin(3) * t72 - pkin(8) * t74 + t86;
t51 = -t100 * t76 + t128 * t77;
t17 = t101 * t51 + t99 * t39;
t13 = -qJ(5) * t126 + t17;
t114 = t101 * t74;
t16 = t101 * t39 - t51 * t99;
t7 = pkin(4) * t72 - qJ(5) * t114 + t16;
t4 = t97 * t13 + t95 * t7;
t127 = Ifges(5,4) * t99;
t123 = -qJ(5) - pkin(8);
t33 = t73 * t74;
t18 = -mrSges(6,2) * t72 - mrSges(6,3) * t33;
t21 = -mrSges(7,2) * t33 + mrSges(7,3) * t72;
t122 = t18 + t21;
t34 = t107 * t74;
t19 = mrSges(6,1) * t72 - mrSges(6,3) * t34;
t20 = -t72 * mrSges(7,1) + t34 * mrSges(7,2);
t121 = -t19 + t20;
t117 = t96 ^ 2 + t92;
t116 = t101 ^ 2 + t99 ^ 2;
t115 = Ifges(5,4) * t101;
t111 = t123 * t99;
t79 = t123 * t101;
t53 = -t97 * t111 - t79 * t95;
t55 = t95 * t111 - t97 * t79;
t113 = t53 ^ 2 + t55 ^ 2;
t87 = -pkin(4) * t101 - pkin(3);
t110 = -t98 * mrSges(3,1) + t96 * mrSges(3,2);
t15 = t33 * mrSges(6,1) + t34 * mrSges(6,2);
t44 = -mrSges(6,1) * t107 + t73 * mrSges(6,2);
t14 = t33 * mrSges(7,1) - t34 * mrSges(7,3);
t43 = -mrSges(7,1) * t107 - t73 * mrSges(7,3);
t24 = pkin(4) * t126 + t49;
t3 = -t13 * t95 + t7 * t97;
t108 = mrSges(5,1) * t99 + mrSges(5,2) * t101;
t78 = -t101 * mrSges(5,1) + t99 * mrSges(5,2);
t106 = -t43 - t44;
t105 = Ifges(5,5) * t114 + t136 * t72 - t137 * t33 + t138 * t34;
t85 = -pkin(4) * t97 - pkin(5);
t82 = pkin(4) * t95 + qJ(6);
t81 = Ifges(5,1) * t99 + t115;
t80 = Ifges(5,2) * t101 + t127;
t60 = t74 * mrSges(4,2);
t48 = Ifges(6,1) * t73 + Ifges(6,4) * t107;
t47 = Ifges(7,1) * t73 - Ifges(7,5) * t107;
t46 = Ifges(6,4) * t73 + Ifges(6,2) * t107;
t45 = Ifges(7,5) * t73 - Ifges(7,3) * t107;
t41 = mrSges(5,1) * t72 - mrSges(5,3) * t114;
t40 = -mrSges(5,2) * t72 - mrSges(5,3) * t126;
t38 = -pkin(5) * t107 - qJ(6) * t73 + t87;
t37 = t108 * t74;
t23 = Ifges(5,5) * t72 + (Ifges(5,1) * t101 - t127) * t74;
t22 = Ifges(5,6) * t72 + (-Ifges(5,2) * t99 + t115) * t74;
t11 = Ifges(6,1) * t34 - Ifges(6,4) * t33 + Ifges(6,5) * t72;
t10 = Ifges(7,1) * t34 + Ifges(7,4) * t72 + Ifges(7,5) * t33;
t9 = Ifges(6,4) * t34 - Ifges(6,2) * t33 + Ifges(6,6) * t72;
t8 = Ifges(7,5) * t34 + Ifges(7,6) * t72 + Ifges(7,3) * t33;
t5 = pkin(5) * t33 - qJ(6) * t34 + t24;
t2 = -pkin(5) * t72 - t3;
t1 = qJ(6) * t72 + t4;
t6 = [Ifges(3,2) * t92 - 0.2e1 * pkin(1) * t110 + t60 * t130 + 0.2e1 * t17 * t40 + 0.2e1 * t16 * t41 + t37 * t131 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t1 * t21 + 0.2e1 * t24 * t15 + 0.2e1 * t5 * t14 + Ifges(2,3) + (Ifges(3,1) * t96 + 0.2e1 * Ifges(3,4) * t98) * t96 + (t10 + t11) * t34 + (t8 - t9) * t33 + 0.2e1 * t117 * qJ(2) * mrSges(3,3) + (mrSges(4,3) * t131 + Ifges(4,1) * t74 + t101 * t23 - t22 * t99) * t74 + (mrSges(4,1) * t130 - 0.2e1 * mrSges(4,3) * t51 + Ifges(4,2) * t72 + (-Ifges(5,6) * t99 - (2 * Ifges(4,4))) * t74 + t105) * t72 + m(3) * (t117 * qJ(2) ^ 2 + pkin(1) ^ 2) + m(4) * (t51 ^ 2 + t86 ^ 2 + t132) + m(5) * (t16 ^ 2 + t17 ^ 2 + t132) + m(6) * (t24 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); -m(3) * pkin(1) + t72 * mrSges(4,1) + t101 * t41 + t99 * t40 + t60 + t122 * t73 - t121 * t107 + m(6) * (t107 * t3 + t4 * t73) + m(7) * (t1 * t73 - t107 * t2) + m(5) * (t101 * t16 + t17 * t99) + m(4) * t86 + t110; m(5) * t116 + m(3) + m(4) + t134 * (t107 ^ 2 + t73 ^ 2); t87 * t15 + Ifges(4,5) * t74 + t5 * t43 + t24 * t44 - t51 * mrSges(4,2) - pkin(3) * t37 + t38 * t14 + t122 * t55 + t121 * t53 + (t78 - mrSges(4,1)) * t49 + (t47 / 0.2e1 + t48 / 0.2e1) * t34 + (t45 / 0.2e1 - t46 / 0.2e1) * t33 + (t23 / 0.2e1 - t16 * mrSges(5,3) - t74 * t80 / 0.2e1 - pkin(8) * t41) * t99 + (t22 / 0.2e1 + t74 * t81 / 0.2e1 + pkin(8) * t40 + t17 * mrSges(5,3)) * t101 + m(5) * (-pkin(3) * t49 + (t17 * t101 - t16 * t99) * pkin(8)) + m(6) * (t24 * t87 - t3 * t53 + t4 * t55) + m(7) * (t1 * t55 + t2 * t53 + t38 * t5) + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1) * t73 - (-t1 * mrSges(7,2) - t4 * mrSges(6,3) + t8 / 0.2e1 - t9 / 0.2e1) * t107 + (-Ifges(4,6) + t133 / 0.2e1) * t72; t134 * (-t107 * t53 + t55 * t73); -0.2e1 * pkin(3) * t78 + t101 * t80 + 0.2e1 * t38 * t43 + 0.2e1 * t87 * t44 + t99 * t81 + Ifges(4,3) + m(7) * (t38 ^ 2 + t113) + m(6) * (t87 ^ 2 + t113) + m(5) * (t116 * pkin(8) ^ 2 + pkin(3) ^ 2) + (t47 + t48) * t73 - (t45 - t46) * t107 + 0.2e1 * t116 * pkin(8) * mrSges(5,3) + 0.2e1 * (t107 * t55 + t53 * t73) * (mrSges(7,2) + mrSges(6,3)); t105 + m(7) * (t1 * t82 + t2 * t85) - Ifges(5,6) * t126 + t82 * t21 + t85 * t20 + t16 * mrSges(5,1) - t17 * mrSges(5,2) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2) + t1 * mrSges(7,3) + (m(6) * (t3 * t97 + t4 * t95) + t95 * t18 + t97 * t19) * pkin(4); m(7) * (-t107 * t85 + t73 * t82) + (t107 * t97 + t73 * t95) * t135 + t106 - t78; m(7) * (t53 * t85 + t55 * t82) - t53 * mrSges(6,1) - t53 * mrSges(7,1) + t55 * mrSges(7,3) - t55 * mrSges(6,2) - t108 * pkin(8) + (t107 * t82 + t73 * t85) * mrSges(7,2) + (m(6) * (-t53 * t97 + t55 * t95) + (t107 * t95 - t73 * t97) * mrSges(6,3)) * pkin(4) + t133; -0.2e1 * t85 * mrSges(7,1) + 0.2e1 * t82 * mrSges(7,3) + m(7) * (t82 ^ 2 + t85 ^ 2) + (0.2e1 * mrSges(6,1) * t97 - 0.2e1 * mrSges(6,2) * t95 + (t95 ^ 2 + t97 ^ 2) * t135) * pkin(4) + t136; m(6) * t24 + m(7) * t5 + t14 + t15; 0; m(6) * t87 + m(7) * t38 - t106; 0; t134; m(7) * t2 + t20; -m(7) * t107; m(7) * t53 + t73 * mrSges(7,2); m(7) * t85 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
