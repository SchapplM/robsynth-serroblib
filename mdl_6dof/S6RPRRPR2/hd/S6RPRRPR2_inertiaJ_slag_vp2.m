% Calculate joint inertia matrix for
% S6RPRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2018-11-23 16:16
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:16:10
% EndTime: 2018-11-23 16:16:11
% DurationCPUTime: 1.05s
% Computational Cost: add. (1779->281), mult. (3513->405), div. (0->0), fcn. (3588->10), ass. (0->111)
t101 = sin(qJ(3));
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t95 = sin(pkin(11));
t97 = cos(pkin(11));
t68 = t100 * t97 + t103 * t95;
t55 = t68 * t101;
t67 = -t100 * t95 + t103 * t97;
t56 = t67 * t101;
t32 = t55 * mrSges(6,1) + t56 * mrSges(6,2);
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t26 = -t102 * t55 - t56 * t99;
t27 = t102 * t56 - t55 * t99;
t8 = t26 * mrSges(7,1) - t27 * mrSges(7,2);
t108 = -t32 + t8;
t110 = mrSges(5,1) * t100 + mrSges(5,2) * t103;
t65 = t110 * t101;
t135 = -t65 + t108;
t96 = sin(pkin(10));
t85 = pkin(1) * t96 + pkin(7);
t134 = 0.2e1 * t85;
t133 = m(6) * pkin(4);
t38 = t102 * t67 - t68 * t99;
t39 = t102 * t68 + t67 * t99;
t14 = -t38 * mrSges(7,1) + t39 * mrSges(7,2);
t40 = -t67 * mrSges(6,1) + t68 * mrSges(6,2);
t132 = -t14 - t40;
t115 = t101 * t103;
t131 = -Ifges(5,5) * t115 - Ifges(6,5) * t56 + Ifges(6,6) * t55;
t130 = m(6) + m(7);
t129 = pkin(4) * t95;
t104 = cos(qJ(3));
t128 = pkin(3) * t104;
t127 = pkin(8) * t101;
t86 = pkin(4) * t97 + pkin(5);
t59 = t102 * t86 - t129 * t99;
t126 = t59 * mrSges(7,1);
t60 = t102 * t129 + t86 * t99;
t125 = t60 * mrSges(7,2);
t124 = -qJ(5) - pkin(8);
t123 = -Ifges(7,5) * t27 - Ifges(7,6) * t26;
t98 = cos(pkin(10));
t87 = -pkin(1) * t98 - pkin(2);
t66 = -t127 + t87 - t128;
t58 = t103 * t66;
t28 = -qJ(5) * t115 + t58 + (-t100 * t85 - pkin(4)) * t104;
t116 = t100 * t101;
t117 = t104 * t85;
t44 = t100 * t66 + t103 * t117;
t33 = -qJ(5) * t116 + t44;
t13 = t95 * t28 + t97 * t33;
t75 = t124 * t100;
t77 = t124 * t103;
t46 = t95 * t75 - t97 * t77;
t76 = -mrSges(5,1) * t103 + mrSges(5,2) * t100;
t122 = t76 - mrSges(4,1);
t81 = t101 * t85;
t61 = pkin(4) * t116 + t81;
t121 = t100 ^ 2 + t103 ^ 2;
t92 = t101 ^ 2;
t94 = t104 ^ 2;
t120 = t92 + t94;
t119 = Ifges(5,4) * t100;
t118 = Ifges(5,4) * t103;
t114 = -Ifges(7,3) - Ifges(5,3) - Ifges(6,3);
t88 = -pkin(4) * t103 - pkin(3);
t113 = t121 * mrSges(5,3);
t12 = t97 * t28 - t33 * t95;
t45 = t97 * t75 + t77 * t95;
t4 = -pkin(5) * t104 - pkin(9) * t56 + t12;
t5 = -pkin(9) * t55 + t13;
t2 = t102 * t4 - t5 * t99;
t3 = t102 * t5 + t4 * t99;
t112 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t123;
t29 = -pkin(9) * t68 + t45;
t30 = pkin(9) * t67 + t46;
t10 = t102 * t29 - t30 * t99;
t11 = t102 * t30 + t29 * t99;
t36 = Ifges(7,6) * t38;
t37 = Ifges(7,5) * t39;
t111 = t10 * mrSges(7,1) - t11 * mrSges(7,2) + t36 + t37;
t43 = -t100 * t117 + t58;
t109 = -t43 * t100 + t44 * t103;
t90 = Ifges(5,5) * t100;
t89 = Ifges(5,6) * t103;
t84 = t85 ^ 2;
t80 = t92 * t84;
t79 = Ifges(5,1) * t100 + t118;
t78 = Ifges(5,2) * t103 + t119;
t73 = -mrSges(5,1) * t104 - mrSges(5,3) * t115;
t72 = mrSges(5,2) * t104 - mrSges(5,3) * t116;
t64 = Ifges(6,5) * t68;
t63 = Ifges(6,6) * t67;
t54 = -Ifges(5,5) * t104 + (Ifges(5,1) * t103 - t119) * t101;
t53 = -Ifges(5,6) * t104 + (-Ifges(5,2) * t100 + t118) * t101;
t49 = -pkin(5) * t67 + t88;
t48 = -mrSges(6,1) * t104 - mrSges(6,3) * t56;
t47 = mrSges(6,2) * t104 - mrSges(6,3) * t55;
t42 = Ifges(6,1) * t68 + Ifges(6,4) * t67;
t41 = Ifges(6,4) * t68 + Ifges(6,2) * t67;
t34 = pkin(5) * t55 + t61;
t23 = Ifges(6,1) * t56 - Ifges(6,4) * t55 - Ifges(6,5) * t104;
t22 = Ifges(6,4) * t56 - Ifges(6,2) * t55 - Ifges(6,6) * t104;
t18 = -mrSges(7,1) * t104 - mrSges(7,3) * t27;
t17 = mrSges(7,2) * t104 + mrSges(7,3) * t26;
t16 = Ifges(7,1) * t39 + Ifges(7,4) * t38;
t15 = Ifges(7,4) * t39 + Ifges(7,2) * t38;
t7 = Ifges(7,1) * t27 + Ifges(7,4) * t26 - Ifges(7,5) * t104;
t6 = Ifges(7,4) * t27 + Ifges(7,2) * t26 - Ifges(7,6) * t104;
t1 = [0.2e1 * t12 * t48 + 0.2e1 * t13 * t47 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t55 * t22 + t56 * t23 + t26 * t6 + t27 * t7 + 0.2e1 * t61 * t32 - 0.2e1 * t34 * t8 + 0.2e1 * t43 * t73 + 0.2e1 * t44 * t72 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * mrSges(4,2) * t87 + Ifges(4,1) * t101 - t100 * t53 + t103 * t54 + t65 * t134) * t101 + (-0.2e1 * t87 * mrSges(4,1) + (Ifges(4,2) - t114) * t104 + (Ifges(5,6) * t100 + (2 * Ifges(4,4))) * t101 + t123 + t131) * t104 + m(4) * (t84 * t94 + t87 ^ 2 + t80) + m(5) * (t43 ^ 2 + t44 ^ 2 + t80) + m(6) * (t12 ^ 2 + t13 ^ 2 + t61 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t34 ^ 2) + m(3) * (t96 ^ 2 + t98 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t98 - mrSges(3,2) * t96) * pkin(1) + t120 * mrSges(4,3) * t134; t27 * t17 + t26 * t18 + t56 * t47 - t55 * t48 + t135 * t104 + m(6) * (-t104 * t61 - t12 * t55 + t13 * t56) + m(7) * (-t104 * t34 + t2 * t26 + t27 * t3) + (t103 * t72 - t100 * t73 + m(5) * (t109 - t117)) * t101; m(3) + m(6) * (t55 ^ 2 + t56 ^ 2 + t94) + m(7) * (t26 ^ 2 + t27 ^ 2 + t94) + m(5) * (t121 * t92 + t94) + m(4) * t120; (t103 * t79 / 0.2e1 - t100 * t78 / 0.2e1 + Ifges(4,5) + t122 * t85) * t101 + m(5) * (-pkin(3) * t81 + pkin(8) * t109) + t88 * t32 + t68 * t23 / 0.2e1 + t61 * t40 - pkin(3) * t65 + t67 * t22 / 0.2e1 - t49 * t8 - t55 * t41 / 0.2e1 + t56 * t42 / 0.2e1 + t38 * t6 / 0.2e1 + t39 * t7 / 0.2e1 + t46 * t47 + t45 * t48 + t26 * t15 / 0.2e1 + t27 * t16 / 0.2e1 + t34 * t14 + t11 * t17 + t10 * t18 + (pkin(8) * t72 + t44 * mrSges(5,3) + t53 / 0.2e1) * t103 + (-pkin(8) * t73 - t43 * mrSges(5,3) + t54 / 0.2e1) * t100 + (-t12 * t68 + t13 * t67) * mrSges(6,3) + (-t2 * t39 + t3 * t38) * mrSges(7,3) + m(6) * (t12 * t45 + t13 * t46 + t61 * t88) + m(7) * (t10 * t2 + t11 * t3 + t34 * t49) + (-t85 * mrSges(4,2) - t90 / 0.2e1 - t89 / 0.2e1 - t64 / 0.2e1 - t63 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 + Ifges(4,6)) * t104; (-t26 * t39 + t27 * t38) * mrSges(7,3) + (t55 * t68 + t56 * t67) * mrSges(6,3) + (-mrSges(4,2) + t113) * t101 + (-t122 + t132) * t104 + m(6) * (-t104 * t88 - t45 * t55 + t46 * t56) + m(7) * (t10 * t26 - t104 * t49 + t11 * t27) + m(5) * (t121 * t127 + t128); -0.2e1 * pkin(3) * t76 + t100 * t79 + t103 * t78 + 0.2e1 * t49 * t14 + t38 * t15 + t39 * t16 + 0.2e1 * t88 * t40 + t67 * t41 + t68 * t42 + Ifges(4,3) + m(7) * (t10 ^ 2 + t11 ^ 2 + t49 ^ 2) + m(6) * (t45 ^ 2 + t46 ^ 2 + t88 ^ 2) + m(5) * (pkin(8) ^ 2 * t121 + pkin(3) ^ 2) + 0.2e1 * (-t10 * t39 + t11 * t38) * mrSges(7,3) + 0.2e1 * (-t45 * t68 + t46 * t67) * mrSges(6,3) + 0.2e1 * pkin(8) * t113; -Ifges(5,6) * t116 + m(7) * (t2 * t59 + t3 * t60) + t59 * t18 + t60 * t17 - t13 * mrSges(6,2) + t12 * mrSges(6,1) - t44 * mrSges(5,2) + t43 * mrSges(5,1) + t114 * t104 + (t97 * t48 + t95 * t47 + m(6) * (t12 * t97 + t13 * t95)) * pkin(4) + t112 - t131; m(7) * (t26 * t59 + t27 * t60) + (-t55 * t97 + t56 * t95) * t133 + t135; m(7) * (t10 * t59 + t11 * t60) + t45 * mrSges(6,1) - t46 * mrSges(6,2) + t64 + t63 + t90 + t89 - t110 * pkin(8) + (t38 * t60 - t39 * t59) * mrSges(7,3) + (m(6) * (t45 * t97 + t46 * t95) + (t67 * t95 - t68 * t97) * mrSges(6,3)) * pkin(4) + t111; 0.2e1 * t126 - 0.2e1 * t125 + m(7) * (t59 ^ 2 + t60 ^ 2) - t114 + (0.2e1 * mrSges(6,1) * t97 - 0.2e1 * mrSges(6,2) * t95 + (t95 ^ 2 + t97 ^ 2) * t133) * pkin(4); m(6) * t61 + m(7) * t34 - t108; -t130 * t104; m(6) * t88 + m(7) * t49 - t132; 0; t130; -Ifges(7,3) * t104 + t112; t8; t111; Ifges(7,3) - t125 + t126; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
