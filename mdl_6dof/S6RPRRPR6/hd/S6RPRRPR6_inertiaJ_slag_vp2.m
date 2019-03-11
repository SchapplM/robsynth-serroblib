% Calculate joint inertia matrix for
% S6RPRRPR6
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:27
% EndTime: 2019-03-09 05:15:29
% DurationCPUTime: 1.08s
% Computational Cost: add. (2685->272), mult. (5204->400), div. (0->0), fcn. (5981->10), ass. (0->103)
t146 = Ifges(5,3) + Ifges(6,3);
t145 = m(6) * pkin(4);
t109 = sin(qJ(6));
t112 = cos(qJ(6));
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t80 = -t105 * t110 + t107 * t113;
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t111 = sin(qJ(3));
t137 = cos(qJ(3));
t83 = t137 * t106 + t111 * t108;
t46 = t80 * t83;
t126 = t113 * t83;
t81 = t106 * t111 - t137 * t108;
t96 = -pkin(2) * t108 - pkin(1);
t50 = pkin(3) * t81 - pkin(8) * t83 + t96;
t133 = pkin(7) + qJ(2);
t87 = t133 * t106;
t88 = t133 * t108;
t64 = -t111 * t87 + t137 * t88;
t33 = -t110 * t64 + t113 * t50;
t15 = pkin(4) * t81 - qJ(5) * t126 + t33;
t127 = t110 * t83;
t34 = t110 * t50 + t113 * t64;
t25 = -qJ(5) * t127 + t34;
t6 = -t105 * t25 + t107 * t15;
t4 = pkin(5) * t81 - pkin(9) * t46 + t6;
t82 = t105 * t113 + t107 * t110;
t45 = t82 * t83;
t7 = t105 * t15 + t107 * t25;
t5 = -pkin(9) * t45 + t7;
t2 = -t109 * t5 + t112 * t4;
t3 = t109 * t4 + t112 * t5;
t144 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t143 = Ifges(5,5) * t110 + Ifges(6,5) * t82 + Ifges(5,6) * t113 + Ifges(6,6) * t80;
t62 = t111 * t88 + t137 * t87;
t142 = t62 ^ 2;
t141 = 0.2e1 * t62;
t140 = 0.2e1 * t96;
t102 = t108 ^ 2;
t136 = pkin(4) * t105;
t95 = pkin(4) * t107 + pkin(5);
t70 = -t109 * t136 + t112 * t95;
t135 = t70 * mrSges(7,1);
t71 = t109 * t95 + t112 * t136;
t134 = t71 * mrSges(7,2);
t132 = -qJ(5) - pkin(8);
t54 = -t109 * t82 + t112 * t80;
t55 = t109 * t80 + t112 * t82;
t131 = Ifges(7,5) * t55 + Ifges(7,6) * t54;
t89 = t132 * t110;
t91 = t132 * t113;
t66 = t105 * t89 - t107 * t91;
t129 = Ifges(5,4) * t110;
t128 = Ifges(5,4) * t113;
t124 = t106 ^ 2 + t102;
t123 = t110 ^ 2 + t113 ^ 2;
t26 = -t109 * t46 - t112 * t45;
t27 = -t109 * t45 + t112 * t46;
t122 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t81;
t97 = -pkin(4) * t113 - pkin(3);
t29 = t45 * mrSges(6,1) + t46 * mrSges(6,2);
t59 = -t80 * mrSges(6,1) + t82 * mrSges(6,2);
t10 = -t26 * mrSges(7,1) + t27 * mrSges(7,2);
t30 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t121 = -t108 * mrSges(3,1) + t106 * mrSges(3,2);
t65 = t105 * t91 + t107 * t89;
t39 = pkin(4) * t127 + t62;
t90 = -t113 * mrSges(5,1) + t110 * mrSges(5,2);
t120 = mrSges(5,1) * t110 + mrSges(5,2) * t113;
t40 = -pkin(9) * t82 + t65;
t41 = pkin(9) * t80 + t66;
t20 = -t109 * t41 + t112 * t40;
t21 = t109 * t40 + t112 * t41;
t119 = t20 * mrSges(7,1) - t21 * mrSges(7,2) + t131;
t118 = -t30 - t59;
t117 = Ifges(5,5) * t126 + Ifges(6,5) * t46 - Ifges(6,6) * t45 + t146 * t81 + t122;
t93 = Ifges(5,1) * t110 + t128;
t92 = Ifges(5,2) * t113 + t129;
t73 = t83 * mrSges(4,2);
t67 = -pkin(5) * t80 + t97;
t61 = Ifges(6,1) * t82 + Ifges(6,4) * t80;
t60 = Ifges(6,4) * t82 + Ifges(6,2) * t80;
t57 = mrSges(5,1) * t81 - mrSges(5,3) * t126;
t56 = -mrSges(5,2) * t81 - mrSges(5,3) * t127;
t49 = t120 * t83;
t38 = Ifges(5,5) * t81 + (Ifges(5,1) * t113 - t129) * t83;
t37 = Ifges(5,6) * t81 + (-Ifges(5,2) * t110 + t128) * t83;
t36 = mrSges(6,1) * t81 - mrSges(6,3) * t46;
t35 = -mrSges(6,2) * t81 - mrSges(6,3) * t45;
t32 = Ifges(7,1) * t55 + Ifges(7,4) * t54;
t31 = Ifges(7,4) * t55 + Ifges(7,2) * t54;
t28 = pkin(5) * t45 + t39;
t17 = Ifges(6,1) * t46 - Ifges(6,4) * t45 + Ifges(6,5) * t81;
t16 = Ifges(6,4) * t46 - Ifges(6,2) * t45 + Ifges(6,6) * t81;
t14 = mrSges(7,1) * t81 - mrSges(7,3) * t27;
t13 = -mrSges(7,2) * t81 + mrSges(7,3) * t26;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t81;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t81;
t1 = [Ifges(3,2) * t102 + (mrSges(4,1) * t140 - 0.2e1 * t64 * mrSges(4,3) + Ifges(4,2) * t81 + (-Ifges(5,6) * t110 - (2 * Ifges(4,4))) * t83 + t117) * t81 + (mrSges(4,3) * t141 + Ifges(4,1) * t83 - t110 * t37 + t113 * t38) * t83 + m(4) * (t64 ^ 2 + t96 ^ 2 + t142) + m(5) * (t33 ^ 2 + t34 ^ 2 + t142) + 0.2e1 * t34 * t56 + 0.2e1 * t33 * t57 + 0.2e1 * t39 * t29 - t45 * t16 + t46 * t17 + t27 * t9 + 0.2e1 * t28 * t10 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t36 + t26 * t8 + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t14 + m(3) * (qJ(2) ^ 2 * t124 + pkin(1) ^ 2) - 0.2e1 * pkin(1) * t121 + 0.2e1 * t124 * qJ(2) * mrSges(3,3) + t73 * t140 + t49 * t141 + Ifges(2,3) + (Ifges(3,1) * t106 + 0.2e1 * Ifges(3,4) * t108) * t106 + m(6) * (t39 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(7) * (t2 ^ 2 + t28 ^ 2 + t3 ^ 2); -m(3) * pkin(1) + t81 * mrSges(4,1) + t110 * t56 + t113 * t57 + t55 * t13 + t54 * t14 + t82 * t35 + t80 * t36 + t73 + m(7) * (t2 * t54 + t3 * t55) + m(6) * (t6 * t80 + t7 * t82) + m(5) * (t110 * t34 + t113 * t33) + m(4) * t96 + t121; m(3) + m(4) + m(5) * t123 + m(6) * (t80 ^ 2 + t82 ^ 2) + m(7) * (t54 ^ 2 + t55 ^ 2); (t131 + t143) * t81 / 0.2e1 + (-t2 * t55 + t3 * t54) * mrSges(7,3) + (-t6 * t82 + t7 * t80) * mrSges(6,3) + t97 * t29 - Ifges(4,6) * t81 + t82 * t17 / 0.2e1 + Ifges(4,5) * t83 + t80 * t16 / 0.2e1 - t45 * t60 / 0.2e1 + t46 * t61 / 0.2e1 - t64 * mrSges(4,2) + t65 * t36 + t66 * t35 + t67 * t10 - pkin(3) * t49 + t54 * t8 / 0.2e1 + t55 * t9 / 0.2e1 + t39 * t59 + t28 * t30 + t26 * t31 / 0.2e1 + t27 * t32 / 0.2e1 + t20 * t14 + t21 * t13 + m(5) * (-pkin(3) * t62 + (-t33 * t110 + t34 * t113) * pkin(8)) + (t90 - mrSges(4,1)) * t62 + (t83 * t93 / 0.2e1 + pkin(8) * t56 + t34 * mrSges(5,3) + t37 / 0.2e1) * t113 + (-t33 * mrSges(5,3) - t83 * t92 / 0.2e1 - pkin(8) * t57 + t38 / 0.2e1) * t110 + m(6) * (t39 * t97 + t6 * t65 + t66 * t7) + m(7) * (t2 * t20 + t21 * t3 + t28 * t67); m(6) * (t65 * t80 + t66 * t82) + m(7) * (t20 * t54 + t21 * t55); -0.2e1 * pkin(3) * t90 + t110 * t93 + t113 * t92 + 0.2e1 * t67 * t30 + t54 * t31 + t55 * t32 + 0.2e1 * t97 * t59 + t80 * t60 + t82 * t61 + Ifges(4,3) + m(7) * (t20 ^ 2 + t21 ^ 2 + t67 ^ 2) + m(6) * (t65 ^ 2 + t66 ^ 2 + t97 ^ 2) + m(5) * (pkin(8) ^ 2 * t123 + pkin(3) ^ 2) + 0.2e1 * (-t20 * t55 + t21 * t54) * mrSges(7,3) + 0.2e1 * (-t65 * t82 + t66 * t80) * mrSges(6,3) + 0.2e1 * t123 * pkin(8) * mrSges(5,3); m(7) * (t2 * t70 + t3 * t71) - Ifges(5,6) * t127 + t70 * t14 + t71 * t13 + t33 * mrSges(5,1) - t34 * mrSges(5,2) + t6 * mrSges(6,1) - t7 * mrSges(6,2) + t117 + (m(6) * (t105 * t7 + t107 * t6) + t107 * t36 + t105 * t35) * pkin(4) + t144; m(7) * (t54 * t70 + t55 * t71) + (t105 * t82 + t107 * t80) * t145 + t118 - t90; m(7) * (t20 * t70 + t21 * t71) + t65 * mrSges(6,1) - t66 * mrSges(6,2) - t120 * pkin(8) + (t54 * t71 - t55 * t70) * mrSges(7,3) + (m(6) * (t105 * t66 + t107 * t65) + (t105 * t80 - t107 * t82) * mrSges(6,3)) * pkin(4) + t119 + t143; 0.2e1 * t135 - 0.2e1 * t134 + Ifges(7,3) + m(7) * (t70 ^ 2 + t71 ^ 2) + (0.2e1 * mrSges(6,1) * t107 - 0.2e1 * mrSges(6,2) * t105 + (t105 ^ 2 + t107 ^ 2) * t145) * pkin(4) + t146; m(6) * t39 + m(7) * t28 + t10 + t29; 0; m(6) * t97 + m(7) * t67 - t118; 0; m(6) + m(7); t122 + t144; -t30; t119; Ifges(7,3) - t134 + t135; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
