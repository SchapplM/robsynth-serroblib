% Calculate joint inertia matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR8_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR8_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:24
% EndTime: 2019-12-31 22:24:27
% DurationCPUTime: 0.75s
% Computational Cost: add. (1416->220), mult. (2721->324), div. (0->0), fcn. (2834->8), ass. (0->92)
t100 = cos(qJ(2));
t95 = sin(qJ(3));
t96 = sin(qJ(2));
t99 = cos(qJ(3));
t71 = t95 * t100 + t99 * t96;
t98 = cos(qJ(4));
t126 = t71 * t98;
t69 = -t99 * t100 + t95 * t96;
t141 = Ifges(5,5) * t126 + Ifges(5,3) * t69;
t94 = sin(qJ(4));
t123 = Ifges(5,5) * t94 + Ifges(5,6) * t98;
t93 = sin(qJ(5));
t97 = cos(qJ(5));
t68 = -t93 * t94 + t97 * t98;
t129 = t68 * mrSges(6,3);
t140 = t93 * pkin(4) * t129 + t123;
t134 = -pkin(7) - pkin(6);
t117 = t134 * t96;
t79 = t134 * t100;
t51 = -t99 * t117 - t95 * t79;
t139 = t51 ^ 2;
t70 = t93 * t98 + t97 * t94;
t42 = -t68 * mrSges(6,1) + t70 * mrSges(6,2);
t138 = 0.2e1 * t42;
t137 = 0.2e1 * t51;
t85 = -t100 * pkin(2) - pkin(1);
t136 = 0.2e1 * t85;
t133 = t99 * pkin(2);
t132 = Ifges(5,4) * t94;
t131 = Ifges(5,4) * t98;
t39 = t69 * pkin(3) - t71 * pkin(8) + t85;
t54 = t95 * t117 - t99 * t79;
t17 = t94 * t39 + t98 * t54;
t130 = t17 * t98;
t128 = t70 * mrSges(6,3);
t127 = t71 * t94;
t125 = t94 * mrSges(5,3);
t124 = Ifges(6,5) * t70 + Ifges(6,6) * t68;
t122 = t94 ^ 2 + t98 ^ 2;
t121 = t100 ^ 2 + t96 ^ 2;
t120 = 0.2e1 * mrSges(6,3);
t119 = t97 * t128;
t30 = t70 * t71;
t31 = t68 * t71;
t118 = Ifges(6,5) * t31 - Ifges(6,6) * t30 + Ifges(6,3) * t69;
t84 = -t98 * pkin(4) - pkin(3);
t82 = t95 * pkin(2) + pkin(8);
t116 = t122 * t82;
t16 = t98 * t39 - t94 * t54;
t43 = Ifges(6,4) * t70 + Ifges(6,2) * t68;
t44 = Ifges(6,1) * t70 + Ifges(6,4) * t68;
t75 = Ifges(5,2) * t98 + t132;
t76 = Ifges(5,1) * t94 + t131;
t115 = t68 * t43 + t70 * t44 + t98 * t75 + t94 * t76 + Ifges(4,3);
t114 = t94 * mrSges(5,1) + t98 * mrSges(5,2);
t113 = -t16 * t94 + t130;
t40 = -t69 * mrSges(5,2) - t71 * t125;
t41 = t69 * mrSges(5,1) - mrSges(5,3) * t126;
t112 = t98 * t40 - t94 * t41;
t111 = 0.2e1 * t122 * mrSges(5,3);
t64 = (-pkin(9) - t82) * t94;
t88 = t98 * pkin(9);
t65 = t98 * t82 + t88;
t36 = t97 * t64 - t93 * t65;
t37 = t93 * t64 + t97 * t65;
t110 = t36 * mrSges(6,1) - t37 * mrSges(6,2) + t124;
t77 = (-pkin(9) - pkin(8)) * t94;
t78 = t98 * pkin(8) + t88;
t50 = t97 * t77 - t93 * t78;
t53 = t93 * t77 + t97 * t78;
t109 = t50 * mrSges(6,1) - t53 * mrSges(6,2) + t124;
t10 = -pkin(9) * t127 + t17;
t7 = t69 * pkin(4) - pkin(9) * t126 + t16;
t3 = -t93 * t10 + t97 * t7;
t4 = t97 * t10 + t93 * t7;
t108 = t3 * mrSges(6,1) - t4 * mrSges(6,2) + t118;
t107 = (t99 * mrSges(4,1) - t95 * mrSges(4,2)) * pkin(2);
t106 = (t97 * mrSges(6,1) - t93 * mrSges(6,2)) * pkin(4);
t22 = Ifges(5,6) * t69 + (-Ifges(5,2) * t94 + t131) * t71;
t23 = Ifges(5,5) * t69 + (Ifges(5,1) * t98 - t132) * t71;
t27 = pkin(4) * t127 + t51;
t74 = -t98 * mrSges(5,1) + t94 * mrSges(5,2);
t8 = Ifges(6,4) * t31 - Ifges(6,2) * t30 + Ifges(6,6) * t69;
t9 = Ifges(6,1) * t31 - Ifges(6,4) * t30 + Ifges(6,5) * t69;
t105 = -t54 * mrSges(4,2) - t16 * t125 - t3 * t128 + t4 * t129 + t68 * t8 / 0.2e1 + t27 * t42 + mrSges(5,3) * t130 - t30 * t43 / 0.2e1 + t31 * t44 / 0.2e1 + t94 * t23 / 0.2e1 + t98 * t22 / 0.2e1 - t75 * t127 / 0.2e1 + t76 * t126 / 0.2e1 + t70 * t9 / 0.2e1 - Ifges(4,6) * t69 + Ifges(4,5) * t71 + (t74 - mrSges(4,1)) * t51 + (t124 + t123) * t69 / 0.2e1;
t83 = -pkin(3) - t133;
t73 = t84 - t133;
t38 = t114 * t71;
t19 = t69 * mrSges(6,1) - t31 * mrSges(6,3);
t18 = -t69 * mrSges(6,2) - t30 * mrSges(6,3);
t11 = t30 * mrSges(6,1) + t31 * mrSges(6,2);
t1 = [t96 * (Ifges(3,1) * t96 + Ifges(3,4) * t100) + t100 * (Ifges(3,4) * t96 + Ifges(3,2) * t100) - 0.2e1 * pkin(1) * (-t100 * mrSges(3,1) + t96 * mrSges(3,2)) + 0.2e1 * t17 * t40 + 0.2e1 * t16 * t41 + t38 * t137 + 0.2e1 * t3 * t19 + 0.2e1 * t27 * t11 - t30 * t8 + t31 * t9 + 0.2e1 * t4 * t18 + Ifges(2,3) + 0.2e1 * t121 * pkin(6) * mrSges(3,3) + (mrSges(4,2) * t136 + mrSges(4,3) * t137 + Ifges(4,1) * t71 - t94 * t22 + t98 * t23) * t71 + (mrSges(4,1) * t136 - 0.2e1 * t54 * mrSges(4,3) + Ifges(4,2) * t69 + (-Ifges(5,6) * t94 - (2 * Ifges(4,4))) * t71 + t118 + t141) * t69 + m(3) * (t121 * pkin(6) ^ 2 + pkin(1) ^ 2) + m(4) * (t54 ^ 2 + t85 ^ 2 + t139) + m(5) * (t16 ^ 2 + t17 ^ 2 + t139) + m(6) * (t27 ^ 2 + t3 ^ 2 + t4 ^ 2); (m(4) * (-t51 * t99 + t54 * t95) + (-t95 * t69 - t99 * t71) * mrSges(4,3)) * pkin(2) + m(5) * (t113 * t82 + t83 * t51) + t112 * t82 + m(6) * (t73 * t27 + t36 * t3 + t37 * t4) + (-t96 * mrSges(3,1) - t100 * mrSges(3,2)) * pkin(6) + Ifges(3,5) * t96 + Ifges(3,6) * t100 + t73 * t11 + t83 * t38 + t36 * t19 + t37 * t18 + t105; t73 * t138 + 0.2e1 * t83 * t74 + Ifges(3,3) + 0.2e1 * t107 + (-t36 * t70 + t37 * t68) * t120 + t82 * t111 + m(6) * (t36 ^ 2 + t37 ^ 2 + t73 ^ 2) + m(5) * (t122 * t82 ^ 2 + t83 ^ 2) + m(4) * (t95 ^ 2 + t99 ^ 2) * pkin(2) ^ 2 + t115; m(5) * (-pkin(3) * t51 + t113 * pkin(8)) + t112 * pkin(8) + m(6) * (t84 * t27 + t50 * t3 + t53 * t4) + t84 * t11 + t53 * t18 - pkin(3) * t38 + t50 * t19 + t105; (t83 - pkin(3)) * t74 + (t73 + t84) * t42 + t107 + m(6) * (t50 * t36 + t53 * t37 + t84 * t73) + m(5) * (-pkin(3) * t83 + pkin(8) * t116) + ((-t36 - t50) * t70 + (t37 + t53) * t68) * mrSges(6,3) + (t122 * pkin(8) + t116) * mrSges(5,3) + t115; -0.2e1 * pkin(3) * t74 + t84 * t138 + (-t50 * t70 + t53 * t68) * t120 + pkin(8) * t111 + m(6) * (t50 ^ 2 + t53 ^ 2 + t84 ^ 2) + m(5) * (t122 * pkin(8) ^ 2 + pkin(3) ^ 2) + t115; -Ifges(5,6) * t127 + t16 * mrSges(5,1) - t17 * mrSges(5,2) + (m(6) * (t3 * t97 + t4 * t93) + t93 * t18 + t97 * t19) * pkin(4) + t108 + t141; -t114 * t82 + (-t119 + m(6) * (t36 * t97 + t37 * t93)) * pkin(4) + t110 + t140; -t114 * pkin(8) + (-t119 + m(6) * (t50 * t97 + t53 * t93)) * pkin(4) + t109 + t140; Ifges(5,3) + Ifges(6,3) + m(6) * (t93 ^ 2 + t97 ^ 2) * pkin(4) ^ 2 + 0.2e1 * t106; t108; t110; t109; Ifges(6,3) + t106; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
