% Calculate joint inertia matrix for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:50:30
% EndTime: 2019-03-09 12:50:33
% DurationCPUTime: 1.12s
% Computational Cost: add. (1350->281), mult. (2375->370), div. (0->0), fcn. (2162->6), ass. (0->102)
t142 = m(7) + m(6);
t141 = pkin(3) + pkin(7);
t140 = -mrSges(6,1) - mrSges(7,1);
t139 = -mrSges(6,2) + mrSges(7,3);
t138 = (-mrSges(7,2) - mrSges(6,3));
t129 = Ifges(7,2) + Ifges(6,3);
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t137 = t96 ^ 2 + t99 ^ 2;
t94 = sin(qJ(5));
t95 = sin(qJ(4));
t97 = cos(qJ(5));
t98 = cos(qJ(4));
t60 = t94 * t95 - t97 * t98;
t61 = t94 * t98 + t95 * t97;
t122 = t60 ^ 2 + t61 ^ 2;
t136 = -m(4) * pkin(2) + mrSges(4,2);
t135 = 2 * mrSges(7,3);
t100 = -pkin(2) - pkin(8);
t118 = -qJ(3) * t96 - pkin(1);
t50 = t100 * t99 + t118;
t72 = t141 * t96;
t64 = t98 * t72;
t13 = pkin(4) * t96 + t64 + (pkin(9) * t99 - t50) * t95;
t130 = t98 * t99;
t21 = t98 * t50 + t95 * t72;
t18 = -pkin(9) * t130 + t21;
t6 = t94 * t13 + t97 * t18;
t134 = Ifges(5,4) * t95;
t133 = Ifges(5,4) * t98;
t132 = Ifges(5,6) * t95;
t131 = t95 * t99;
t128 = -pkin(9) + t100;
t45 = -t97 * t130 + t94 * t131;
t34 = mrSges(7,2) * t45 + mrSges(7,3) * t96;
t35 = -mrSges(6,2) * t96 + mrSges(6,3) * t45;
t127 = t34 + t35;
t46 = t61 * t99;
t36 = mrSges(6,1) * t96 + mrSges(6,3) * t46;
t37 = -t96 * mrSges(7,1) - t46 * mrSges(7,2);
t126 = -t36 + t37;
t125 = t137 * pkin(7) ^ 2;
t73 = t141 * t99;
t124 = -t95 ^ 2 - t98 ^ 2;
t77 = t95 * pkin(4) + qJ(3);
t120 = t128 * t98;
t67 = t128 * t95;
t31 = -t97 * t120 + t67 * t94;
t33 = t94 * t120 + t97 * t67;
t123 = t31 ^ 2 + t33 ^ 2;
t48 = pkin(4) * t130 + t73;
t121 = m(5) * t124;
t119 = t124 * mrSges(5,3);
t116 = 2 * t138;
t114 = t98 * mrSges(5,1) - t95 * mrSges(5,2);
t113 = -Ifges(5,5) * t95 - Ifges(5,6) * t98;
t112 = pkin(5) * t60 - qJ(6) * t61;
t5 = t13 * t97 - t18 * t94;
t20 = -t50 * t95 + t64;
t111 = t20 * t98 + t21 * t95;
t76 = pkin(4) * t94 + qJ(6);
t79 = -pkin(4) * t97 - pkin(5);
t110 = t60 * t79 + t61 * t76;
t109 = t60 * t97 - t61 * t94;
t108 = t129 * t96 + (-Ifges(7,4) - Ifges(6,5)) * t46 + (Ifges(6,6) - Ifges(7,6)) * t45;
t107 = (mrSges(6,1) * t97 - mrSges(6,2) * t94) * pkin(4);
t106 = t139 * t61 + t140 * t60;
t53 = Ifges(7,6) * t61;
t54 = Ifges(6,6) * t61;
t55 = Ifges(6,5) * t60;
t56 = Ifges(7,4) * t60;
t105 = t139 * t33 + t140 * t31 + t53 - t54 - t55 - t56;
t2 = qJ(6) * t96 + t6;
t3 = -pkin(5) * t96 - t5;
t104 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) + t108;
t101 = qJ(3) ^ 2;
t83 = Ifges(5,5) * t98;
t81 = Ifges(5,3) * t96;
t71 = Ifges(5,1) * t98 - t134;
t70 = -Ifges(5,2) * t95 + t133;
t69 = mrSges(5,1) * t95 + mrSges(5,2) * t98;
t68 = -pkin(2) * t99 + t118;
t66 = -mrSges(5,2) * t96 - mrSges(5,3) * t130;
t65 = mrSges(5,1) * t96 + mrSges(5,3) * t131;
t49 = t114 * t99;
t44 = Ifges(5,5) * t96 + (-Ifges(5,1) * t95 - t133) * t99;
t43 = t96 * Ifges(5,6) + (-Ifges(5,2) * t98 - t134) * t99;
t29 = -Ifges(6,1) * t60 - Ifges(6,4) * t61;
t28 = -Ifges(7,1) * t60 + Ifges(7,5) * t61;
t27 = -Ifges(6,4) * t60 - Ifges(6,2) * t61;
t26 = -Ifges(7,5) * t60 + Ifges(7,3) * t61;
t25 = mrSges(6,1) * t61 - mrSges(6,2) * t60;
t24 = mrSges(7,1) * t61 + mrSges(7,3) * t60;
t19 = pkin(5) * t61 + qJ(6) * t60 + t77;
t17 = -mrSges(6,1) * t45 - mrSges(6,2) * t46;
t16 = -mrSges(7,1) * t45 + mrSges(7,3) * t46;
t12 = -Ifges(6,1) * t46 + Ifges(6,4) * t45 + Ifges(6,5) * t96;
t11 = -Ifges(7,1) * t46 + Ifges(7,4) * t96 - Ifges(7,5) * t45;
t10 = -Ifges(6,4) * t46 + Ifges(6,2) * t45 + Ifges(6,6) * t96;
t9 = -Ifges(7,5) * t46 + Ifges(7,6) * t96 - Ifges(7,3) * t45;
t7 = -pkin(5) * t45 + qJ(6) * t46 + t48;
t1 = [0.2e1 * t7 * t16 + 0.2e1 * t48 * t17 + 0.2e1 * t2 * t34 + 0.2e1 * t20 * t65 + 0.2e1 * t21 * t66 + 0.2e1 * t3 * t37 + 0.2e1 * t6 * t35 + 0.2e1 * t5 * t36 + 0.2e1 * t73 * t49 + Ifges(2,3) - (t11 + t12) * t46 + (t10 - t9) * t45 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * t68 * mrSges(4,2) - t98 * t43 - t95 * t44 + (Ifges(4,3) + Ifges(3,2)) * t99) * t99 + (-0.2e1 * pkin(1) * mrSges(3,2) - 0.2e1 * t68 * mrSges(4,3) + t81 + (Ifges(3,1) + Ifges(4,2)) * t96 + ((2 * Ifges(3,4)) + (2 * Ifges(4,6)) + t113) * t99 + t108) * t96 + m(4) * (t68 ^ 2 + t125) + m(3) * (pkin(1) ^ 2 + t125) + m(5) * (t20 ^ 2 + t21 ^ 2 + t73 ^ 2) + m(6) * (t48 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + 0.2e1 * (mrSges(4,1) + mrSges(3,3)) * pkin(7) * t137; qJ(3) * t49 + t19 * t16 + t77 * t17 + t7 * t24 + t48 * t25 + t73 * t69 - (t29 / 0.2e1 + t28 / 0.2e1) * t46 + (t27 / 0.2e1 - t26 / 0.2e1) * t45 + t127 * t33 + t126 * t31 + (t100 * t65 - t20 * mrSges(5,3) + t44 / 0.2e1) * t98 + (t100 * t66 - t21 * mrSges(5,3) - t43 / 0.2e1) * t95 + (-pkin(2) * mrSges(4,1) - t132 / 0.2e1 + t83 / 0.2e1 - t55 / 0.2e1 - t54 / 0.2e1 - t56 / 0.2e1 + t53 / 0.2e1 - Ifges(4,4) + Ifges(3,5)) * t96 + m(5) * (qJ(3) * t73 + t111 * t100) + m(6) * (-t31 * t5 + t33 * t6 + t48 * t77) + m(7) * (t19 * t7 + t2 * t33 + t3 * t31) + (-t10 / 0.2e1 + t9 / 0.2e1 - t2 * mrSges(7,2) - t6 * mrSges(6,3)) * t61 + (-t11 / 0.2e1 - t12 / 0.2e1 - t3 * mrSges(7,2) + t5 * mrSges(6,3)) * t60 + (-t98 * t70 / 0.2e1 - t95 * t71 / 0.2e1 + qJ(3) * mrSges(4,1) - Ifges(4,5) + Ifges(3,6)) * t99 + ((m(4) * qJ(3) - mrSges(3,2) + mrSges(4,3)) * t99 + (-mrSges(3,1) + t136) * t96) * pkin(7); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t19 * t24 + 0.2e1 * t77 * t25 - t95 * t70 + t98 * t71 + Ifges(4,1) + Ifges(3,3) + (t33 * t116 + t26 - t27) * t61 + (t31 * t116 - t28 - t29) * t60 + m(7) * (t19 ^ 2 + t123) + m(6) * (t77 ^ 2 + t123) + m(5) * (-t124 * t100 ^ 2 + t101) + m(4) * (pkin(2) ^ 2 + t101) + 0.2e1 * (t69 + mrSges(4,3)) * qJ(3) + 0.2e1 * t100 * t119; t98 * t65 + t95 * t66 + (m(4) * pkin(7) + mrSges(4,1)) * t96 + t127 * t61 + t126 * t60 + m(7) * (t2 * t61 + t3 * t60) + m(6) * (-t5 * t60 + t6 * t61) + m(5) * t111; t119 - t100 * t121 + t142 * (t31 * t60 + t61 * t33) + t136 + t138 * t122; t142 * t122 + m(4) - t121; (m(6) * (t5 * t97 + t6 * t94) + t94 * t35 + t97 * t36) * pkin(4) + m(7) * (t2 * t76 + t3 * t79) + t113 * t99 + t81 + t79 * t37 + t76 * t34 + t104 + t20 * mrSges(5,1) - t21 * mrSges(5,2); m(7) * (t31 * t79 + t33 * t76) + t83 - t132 + t114 * t100 - t110 * mrSges(7,2) + (m(6) * (-t31 * t97 + t33 * t94) + t109 * mrSges(6,3)) * pkin(4) + t105; -m(6) * t109 * pkin(4) + m(7) * t110 + t106 + t114; -0.2e1 * t79 * mrSges(7,1) + t76 * t135 + Ifges(5,3) + 0.2e1 * t107 + m(7) * (t76 ^ 2 + t79 ^ 2) + m(6) * (t94 ^ 2 + t97 ^ 2) * pkin(4) ^ 2 + t129; -pkin(5) * t37 + qJ(6) * t34 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) + t104; m(7) * (-pkin(5) * t31 + qJ(6) * t33) + t112 * mrSges(7,2) + t105; -m(7) * t112 + t106; m(7) * (-pkin(5) * t79 + qJ(6) * t76) + t107 + (t76 + qJ(6)) * mrSges(7,3) + (-t79 + pkin(5)) * mrSges(7,1) + t129; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t135 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t129; m(7) * t3 + t37; m(7) * t31 - t60 * mrSges(7,2); m(7) * t60; m(7) * t79 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
