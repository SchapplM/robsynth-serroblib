% Calculate joint inertia matrix for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:20:04
% EndTime: 2019-03-09 05:20:06
% DurationCPUTime: 0.90s
% Computational Cost: add. (1667->207), mult. (2857->293), div. (0->0), fcn. (3118->8), ass. (0->87)
t81 = sin(qJ(6));
t84 = cos(qJ(6));
t112 = t81 ^ 2 + t84 ^ 2;
t142 = mrSges(7,3) * t112;
t59 = -mrSges(7,1) * t84 + mrSges(7,2) * t81;
t141 = t59 - mrSges(6,1);
t82 = sin(qJ(4));
t83 = sin(qJ(3));
t85 = cos(qJ(4));
t86 = cos(qJ(3));
t56 = -t82 * t86 - t85 * t83;
t58 = -t82 * t83 + t85 * t86;
t79 = sin(pkin(10));
t80 = cos(pkin(10));
t36 = t56 * t79 + t80 * t58;
t97 = -t80 * t56 + t58 * t79;
t140 = t36 * t80 + t79 * t97;
t128 = pkin(3) * t82;
t70 = pkin(3) * t85 + pkin(4);
t44 = -t128 * t79 + t70 * t80;
t45 = t80 * t128 + t79 * t70;
t139 = t36 * t44 + t45 * t97;
t138 = t86 ^ 2;
t137 = m(6) * pkin(4);
t117 = t81 * mrSges(7,3);
t15 = -mrSges(7,2) * t97 - t117 * t36;
t120 = t36 * t84;
t16 = mrSges(7,1) * t97 - mrSges(7,3) * t120;
t99 = t84 * t15 - t81 * t16;
t87 = -pkin(1) - pkin(7);
t126 = -pkin(8) + t87;
t107 = t126 * t86;
t108 = t126 * t83;
t39 = t82 * t107 + t85 * t108;
t22 = qJ(5) * t56 + t39;
t38 = t107 * t85 - t108 * t82;
t92 = -t58 * qJ(5) + t38;
t10 = t22 * t79 - t80 * t92;
t134 = t10 ^ 2;
t133 = t36 ^ 2;
t132 = 0.2e1 * t10;
t67 = t83 * pkin(3) + qJ(2);
t40 = -pkin(4) * t56 + t67;
t131 = 0.2e1 * t40;
t130 = 0.2e1 * t59;
t12 = t80 * t22 + t79 * t92;
t13 = pkin(5) * t97 - pkin(9) * t36 + t40;
t3 = t12 * t84 + t13 * t81;
t127 = t3 * t84;
t124 = Ifges(7,4) * t81;
t123 = Ifges(7,4) * t84;
t122 = t10 * t36;
t121 = t36 * t81;
t119 = t44 * mrSges(6,1);
t118 = t45 * mrSges(6,2);
t114 = Ifges(7,5) * t120 + Ifges(7,3) * t97;
t113 = Ifges(7,5) * t81 + Ifges(7,6) * t84;
t111 = t83 ^ 2 + t138;
t110 = t56 ^ 2 + t58 ^ 2;
t109 = m(4) * t111;
t43 = pkin(9) + t45;
t106 = t112 * t43;
t65 = pkin(4) * t79 + pkin(9);
t105 = t112 * t65;
t104 = t111 * mrSges(4,3);
t60 = Ifges(7,2) * t84 + t124;
t61 = Ifges(7,1) * t81 + t123;
t103 = t84 * t60 + t81 * t61 + Ifges(5,3) + Ifges(6,3);
t2 = -t12 * t81 + t13 * t84;
t102 = -t2 * t81 + t127;
t101 = t80 * mrSges(6,1) - t79 * mrSges(6,2);
t100 = mrSges(7,1) * t81 + mrSges(7,2) * t84;
t98 = t38 * t58 - t39 * t56;
t96 = t56 * t82 - t58 * t85;
t95 = 0.2e1 * t142;
t94 = (mrSges(5,1) * t85 - mrSges(5,2) * t82) * pkin(3);
t93 = t58 * mrSges(5,1) + t56 * mrSges(5,2) - t141 * t36 + (-mrSges(6,2) + t142) * t97;
t6 = Ifges(7,6) * t97 + (-Ifges(7,2) * t81 + t123) * t36;
t7 = Ifges(7,5) * t97 + (Ifges(7,1) * t84 - t124) * t36;
t91 = -t39 * mrSges(5,2) - t12 * mrSges(6,2) + mrSges(7,3) * t127 - t117 * t2 - t60 * t121 / 0.2e1 + t61 * t120 / 0.2e1 + Ifges(6,5) * t36 + t38 * mrSges(5,1) + t81 * t7 / 0.2e1 + Ifges(5,6) * t56 + t84 * t6 / 0.2e1 + Ifges(5,5) * t58 + (t113 / 0.2e1 - Ifges(6,6)) * t97 + t141 * t10;
t88 = qJ(2) ^ 2;
t66 = -pkin(4) * t80 - pkin(5);
t42 = -pkin(5) - t44;
t31 = t97 ^ 2;
t27 = t36 * mrSges(6,2);
t14 = t100 * t36;
t1 = [Ifges(4,1) * t138 + 0.2e1 * t67 * (-mrSges(5,1) * t56 + mrSges(5,2) * t58) + t58 * (Ifges(5,1) * t58 + Ifges(5,4) * t56) + t56 * (Ifges(5,4) * t58 + Ifges(5,2) * t56) + t27 * t131 + t14 * t132 + 0.2e1 * t3 * t15 + 0.2e1 * t2 * t16 + Ifges(3,1) + Ifges(2,3) - (2 * pkin(1) * mrSges(3,2)) + (mrSges(6,1) * t131 - 0.2e1 * t12 * mrSges(6,3) + Ifges(6,2) * t97 + t114) * t97 + (mrSges(6,3) * t132 + Ifges(6,1) * t36 - t81 * t6 + t84 * t7 + (-Ifges(7,6) * t81 - (2 * Ifges(6,4))) * t97) * t36 + m(4) * (t111 * t87 ^ 2 + t88) + m(3) * ((pkin(1) ^ 2) + t88) + m(5) * (t38 ^ 2 + t39 ^ 2 + t67 ^ 2) + m(6) * (t12 ^ 2 + t40 ^ 2 + t134) + m(7) * (t2 ^ 2 + t3 ^ 2 + t134) + (-0.2e1 * Ifges(4,4) * t86 + Ifges(4,2) * t83) * t83 + 0.2e1 * (mrSges(4,1) * t83 + mrSges(4,2) * t86 + mrSges(3,3)) * qJ(2) - 0.2e1 * mrSges(5,3) * t98 - 0.2e1 * t104 * t87; -m(3) * pkin(1) + mrSges(3,2) - (t36 * mrSges(6,3) + t14) * t36 - t110 * mrSges(5,3) - t104 + (-mrSges(6,3) * t97 + t99) * t97 + m(7) * (t102 * t97 - t122) + m(6) * (t12 * t97 - t122) + m(5) * t98 + t87 * t109; m(3) + m(7) * (t112 * t31 + t133) + m(6) * (t31 + t133) + m(5) * t110 + t109; t99 * t43 + m(6) * (-t10 * t44 + t12 * t45) + (mrSges(4,1) * t87 + Ifges(4,5)) * t86 + (-mrSges(4,2) * t87 - Ifges(4,6)) * t83 - t139 * mrSges(6,3) + (m(5) * (t38 * t85 + t39 * t82) + t96 * mrSges(5,3)) * pkin(3) + m(7) * (t10 * t42 + t102 * t43) + t91 + t42 * t14; t86 * mrSges(4,1) - t83 * mrSges(4,2) + m(7) * (t106 * t97 - t36 * t42) + m(6) * t139 - m(5) * t96 * pkin(3) + t93; 0.2e1 * t119 - 0.2e1 * t118 + t42 * t130 + Ifges(4,3) + 0.2e1 * t94 + t43 * t95 + m(6) * (t44 ^ 2 + t45 ^ 2) + m(7) * (t112 * t43 ^ 2 + t42 ^ 2) + m(5) * (t82 ^ 2 + t85 ^ 2) * pkin(3) ^ 2 + t103; (m(6) * (-t10 * t80 + t12 * t79) - t140 * mrSges(6,3)) * pkin(4) + t91 + (m(7) * t10 + t14) * t66 + (m(7) * t102 + t99) * t65; m(7) * (t105 * t97 - t36 * t66) + t140 * t137 + t93; m(7) * (t105 * t43 + t42 * t66) - t118 + t119 + (t42 + t66) * t59 + t94 + (m(6) * (t44 * t80 + t45 * t79) + t101) * pkin(4) + (t105 + t106) * mrSges(7,3) + t103; t66 * t130 + t65 * t95 + m(7) * (t112 * t65 ^ 2 + t66 ^ 2) + t103 + (0.2e1 * t101 + (t79 ^ 2 + t80 ^ 2) * t137) * pkin(4); t97 * mrSges(6,1) + t81 * t15 + t84 * t16 + t27 + m(7) * (t2 * t84 + t3 * t81) + m(6) * t40; 0; 0; 0; m(7) * t112 + m(6); mrSges(7,1) * t2 - mrSges(7,2) * t3 - Ifges(7,6) * t121 + t114; -t100 * t97; -t100 * t43 + t113; -t100 * t65 + t113; -t59; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
