% Calculate joint inertia matrix for
% S6RPRRPR8
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:42
% EndTime: 2019-03-09 05:22:44
% DurationCPUTime: 1.06s
% Computational Cost: add. (1681->292), mult. (3287->416), div. (0->0), fcn. (3354->8), ass. (0->114)
t142 = Ifges(5,3) + Ifges(6,3);
t109 = (-pkin(1) - pkin(7));
t141 = -2 * t109;
t140 = m(6) * pkin(4);
t108 = cos(qJ(3));
t101 = sin(pkin(10));
t102 = cos(pkin(10));
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t71 = t101 * t107 + t102 * t104;
t59 = t71 * t108;
t70 = -t101 * t104 + t102 * t107;
t61 = t70 * t108;
t31 = t59 * mrSges(6,1) + t61 * mrSges(6,2);
t103 = sin(qJ(6));
t106 = cos(qJ(6));
t26 = -t103 * t61 - t106 * t59;
t28 = -t103 * t59 + t106 * t61;
t7 = -t26 * mrSges(7,1) + t28 * mrSges(7,2);
t139 = -t31 - t7;
t105 = sin(qJ(3));
t122 = t107 * t108;
t78 = pkin(3) * t105 - pkin(8) * t108 + qJ(2);
t69 = t107 * t78;
t37 = -qJ(5) * t122 + t69 + (-t104 * t109 + pkin(4)) * t105;
t124 = t104 * t108;
t123 = t105 * t109;
t51 = t104 * t78 + t107 * t123;
t45 = -qJ(5) * t124 + t51;
t12 = -t101 * t45 + t102 * t37;
t6 = pkin(5) * t105 - pkin(9) * t61 + t12;
t13 = t101 * t37 + t102 * t45;
t9 = -pkin(9) * t59 + t13;
t2 = -t103 * t9 + t106 * t6;
t3 = t103 * t6 + t106 * t9;
t138 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t38 = -t103 * t71 + t106 * t70;
t39 = t103 * t70 + t106 * t71;
t14 = -t38 * mrSges(7,1) + t39 * mrSges(7,2);
t42 = -t70 * mrSges(6,1) + t71 * mrSges(6,2);
t137 = -t14 - t42;
t136 = 0.2e1 * qJ(2);
t135 = m(6) + m(7);
t133 = pkin(4) * t101;
t87 = pkin(4) * t102 + pkin(5);
t62 = -t103 * t133 + t106 * t87;
t132 = t62 * mrSges(7,1);
t63 = t103 * t87 + t106 * t133;
t131 = t63 * mrSges(7,2);
t130 = -qJ(5) - pkin(8);
t79 = t130 * t104;
t81 = t130 * t107;
t47 = t101 * t79 - t102 * t81;
t80 = -mrSges(5,1) * t107 + mrSges(5,2) * t104;
t129 = -t80 + mrSges(4,1);
t128 = t104 ^ 2 + t107 ^ 2;
t97 = t105 ^ 2;
t99 = t108 ^ 2;
t127 = t99 + t97;
t126 = Ifges(5,4) * t104;
t125 = Ifges(5,4) * t107;
t121 = t108 * t109;
t120 = Ifges(7,5) * t28 + Ifges(7,6) * t26 + Ifges(7,3) * t105;
t88 = -pkin(4) * t107 - pkin(3);
t119 = t128 * mrSges(5,3);
t118 = t127 * mrSges(4,3);
t58 = t71 * t105;
t60 = t70 * t105;
t25 = -t103 * t60 - t106 * t58;
t27 = -t103 * t58 + t106 * t60;
t117 = t25 * mrSges(7,1) - t27 * mrSges(7,2);
t46 = t101 * t81 + t102 * t79;
t72 = pkin(4) * t124 - t121;
t29 = -pkin(9) * t71 + t46;
t30 = pkin(9) * t70 + t47;
t10 = -t103 * t30 + t106 * t29;
t11 = t103 * t29 + t106 * t30;
t35 = Ifges(7,6) * t38;
t36 = Ifges(7,5) * t39;
t116 = t10 * mrSges(7,1) - t11 * mrSges(7,2) + t35 + t36;
t115 = mrSges(5,1) * t104 + mrSges(5,2) * t107;
t50 = -t104 * t123 + t69;
t114 = -t50 * t104 + t51 * t107;
t113 = Ifges(5,5) * t122 + Ifges(6,5) * t61 - Ifges(6,6) * t59 + t142 * t105 + t120;
t110 = qJ(2) ^ 2;
t100 = t109 ^ 2;
t95 = Ifges(5,5) * t104;
t94 = Ifges(5,6) * t107;
t90 = t99 * t109;
t89 = t99 * t100;
t83 = Ifges(5,1) * t104 + t125;
t82 = Ifges(5,2) * t107 + t126;
t77 = mrSges(5,1) * t105 - mrSges(5,3) * t122;
t76 = -mrSges(5,2) * t105 - mrSges(5,3) * t124;
t67 = t115 * t108;
t66 = Ifges(6,5) * t71;
t65 = Ifges(6,6) * t70;
t57 = Ifges(5,5) * t105 + (Ifges(5,1) * t107 - t126) * t108;
t56 = Ifges(5,6) * t105 + (-Ifges(5,2) * t104 + t125) * t108;
t52 = -pkin(5) * t70 + t88;
t49 = mrSges(6,1) * t105 - mrSges(6,3) * t61;
t48 = -mrSges(6,2) * t105 - mrSges(6,3) * t59;
t44 = Ifges(6,1) * t71 + Ifges(6,4) * t70;
t43 = Ifges(6,4) * t71 + Ifges(6,2) * t70;
t41 = pkin(5) * t59 + t72;
t24 = Ifges(6,1) * t61 - Ifges(6,4) * t59 + Ifges(6,5) * t105;
t23 = Ifges(6,4) * t61 - Ifges(6,2) * t59 + Ifges(6,6) * t105;
t18 = mrSges(7,1) * t105 - mrSges(7,3) * t28;
t17 = -mrSges(7,2) * t105 + mrSges(7,3) * t26;
t16 = Ifges(7,1) * t39 + Ifges(7,4) * t38;
t15 = Ifges(7,4) * t39 + Ifges(7,2) * t38;
t5 = Ifges(7,1) * t28 + Ifges(7,4) * t26 + Ifges(7,5) * t105;
t4 = Ifges(7,4) * t28 + Ifges(7,2) * t26 + Ifges(7,6) * t105;
t1 = [-(2 * pkin(1) * mrSges(3,2)) + mrSges(3,3) * t136 + 0.2e1 * t12 * t49 + 0.2e1 * t13 * t48 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 - t59 * t23 + t61 * t24 + t26 * t4 + t28 * t5 + 0.2e1 * t72 * t31 + 0.2e1 * t41 * t7 + 0.2e1 * t50 * t77 + 0.2e1 * t51 * t76 + Ifges(3,1) + Ifges(2,3) + t118 * t141 + (mrSges(4,2) * t136 + Ifges(4,1) * t108 - t104 * t56 + t107 * t57 + t67 * t141) * t108 + (mrSges(4,1) * t136 + Ifges(4,2) * t105 + (-Ifges(5,6) * t104 - (2 * Ifges(4,4))) * t108 + t113) * t105 + m(4) * (t100 * t97 + t110 + t89) + m(3) * ((pkin(1) ^ 2) + t110) + m(5) * (t50 ^ 2 + t51 ^ 2 + t89) + m(6) * (t12 ^ 2 + t13 ^ 2 + t72 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t41 ^ 2); -m(3) * pkin(1) + t27 * t17 + t25 * t18 + t60 * t48 - t58 * t49 + mrSges(3,2) + (-t104 * t77 + t107 * t76) * t105 - t118 + (-t67 + t139) * t108 + m(7) * (-t108 * t41 + t2 * t25 + t27 * t3) + m(6) * (-t108 * t72 - t12 * t58 + t13 * t60) + m(5) * (t105 * t114 + t90) + m(4) * (t109 * t97 + t90); m(3) + m(7) * (t25 ^ 2 + t27 ^ 2 + t99) + m(6) * (t58 ^ 2 + t60 ^ 2 + t99) + m(5) * (t128 * t97 + t99) + m(4) * t127; (t95 / 0.2e1 + t94 / 0.2e1 + t66 / 0.2e1 + t65 / 0.2e1 + t36 / 0.2e1 + t35 / 0.2e1 - Ifges(4,6) - t109 * mrSges(4,2)) * t105 + (-t2 * t39 + t3 * t38) * mrSges(7,3) + (-t12 * t71 + t13 * t70) * mrSges(6,3) + (Ifges(4,5) + t107 * t83 / 0.2e1 - t104 * t82 / 0.2e1 + t129 * t109) * t108 + m(5) * (pkin(3) * t121 + pkin(8) * t114) + t88 * t31 - pkin(3) * t67 + t70 * t23 / 0.2e1 + t71 * t24 / 0.2e1 + t72 * t42 - t59 * t43 / 0.2e1 + t61 * t44 / 0.2e1 + t47 * t48 + t46 * t49 + t52 * t7 + t38 * t4 / 0.2e1 + t39 * t5 / 0.2e1 + t41 * t14 + t26 * t15 / 0.2e1 + t28 * t16 / 0.2e1 + t11 * t17 + t10 * t18 + m(6) * (t12 * t46 + t13 * t47 + t72 * t88) + m(7) * (t10 * t2 + t11 * t3 + t41 * t52) + (t57 / 0.2e1 - pkin(8) * t77 - t50 * mrSges(5,3)) * t104 + (t56 / 0.2e1 + pkin(8) * t76 + t51 * mrSges(5,3)) * t107; (-t25 * t39 + t27 * t38) * mrSges(7,3) + (t58 * t71 + t60 * t70) * mrSges(6,3) + (-mrSges(4,2) + t119) * t105 + (t129 + t137) * t108 + m(7) * (t10 * t25 - t108 * t52 + t11 * t27) + m(6) * (-t108 * t88 - t46 * t58 + t47 * t60) + m(5) * (pkin(8) * t105 * t128 + pkin(3) * t108); -0.2e1 * pkin(3) * t80 + t104 * t83 + t107 * t82 + 0.2e1 * t52 * t14 + t38 * t15 + t39 * t16 + 0.2e1 * t88 * t42 + t70 * t43 + t71 * t44 + Ifges(4,3) + m(7) * (t10 ^ 2 + t11 ^ 2 + t52 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t88 ^ 2) + m(5) * (pkin(8) ^ 2 * t128 + pkin(3) ^ 2) + 0.2e1 * (-t10 * t39 + t11 * t38) * mrSges(7,3) + 0.2e1 * (-t46 * t71 + t47 * t70) * mrSges(6,3) + 0.2e1 * pkin(8) * t119; m(7) * (t2 * t62 + t3 * t63) + t113 + (m(6) * (t101 * t13 + t102 * t12) + t101 * t48 + t102 * t49) * pkin(4) + t63 * t17 + t62 * t18 + t50 * mrSges(5,1) - t51 * mrSges(5,2) + t12 * mrSges(6,1) - t13 * mrSges(6,2) - Ifges(5,6) * t124 + t138; -t58 * mrSges(6,1) - t60 * mrSges(6,2) - t115 * t105 + m(7) * (t25 * t62 + t27 * t63) + (t101 * t60 - t102 * t58) * t140 + t117; m(7) * (t10 * t62 + t11 * t63) - t47 * mrSges(6,2) + t46 * mrSges(6,1) + t95 + t94 + t66 + t65 - t115 * pkin(8) + (t38 * t63 - t39 * t62) * mrSges(7,3) + (m(6) * (t101 * t47 + t102 * t46) + (t101 * t70 - t102 * t71) * mrSges(6,3)) * pkin(4) + t116; 0.2e1 * t132 - 0.2e1 * t131 + Ifges(7,3) + m(7) * (t62 ^ 2 + t63 ^ 2) + (0.2e1 * mrSges(6,1) * t102 - 0.2e1 * mrSges(6,2) * t101 + (t101 ^ 2 + t102 ^ 2) * t140) * pkin(4) + t142; m(6) * t72 + m(7) * t41 - t139; -t135 * t108; m(6) * t88 + m(7) * t52 - t137; 0; t135; t120 + t138; t117; t116; Ifges(7,3) - t131 + t132; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
