% Calculate joint inertia matrix for
% S6RRRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:56:45
% EndTime: 2019-03-09 15:56:48
% DurationCPUTime: 1.33s
% Computational Cost: add. (1802->333), mult. (3418->463), div. (0->0), fcn. (3378->8), ass. (0->122)
t117 = sin(qJ(3));
t120 = cos(qJ(3));
t154 = t117 ^ 2 + t120 ^ 2;
t153 = 2 * pkin(7);
t116 = sin(qJ(6));
t119 = cos(qJ(6));
t114 = sin(pkin(10));
t115 = cos(pkin(10));
t121 = cos(qJ(2));
t109 = t121 * pkin(3);
t118 = sin(qJ(2));
t81 = -pkin(2) * t121 - pkin(8) * t118 - pkin(1);
t146 = pkin(7) * t121;
t95 = t117 * t146;
t30 = pkin(4) * t121 + t109 + t95 + (-qJ(5) * t118 - t81) * t120;
t137 = t117 * t118;
t52 = t117 * t81 + t120 * t146;
t47 = -qJ(4) * t121 + t52;
t37 = qJ(5) * t137 + t47;
t10 = -t114 * t37 + t115 * t30;
t127 = t114 * t117 + t115 * t120;
t60 = t127 * t118;
t3 = pkin(5) * t121 - pkin(9) * t60 + t10;
t11 = t114 * t30 + t115 * t37;
t67 = -t114 * t120 + t115 * t117;
t59 = t67 * t118;
t6 = pkin(9) * t59 + t11;
t1 = -t116 * t6 + t119 * t3;
t2 = t116 * t3 + t119 * t6;
t152 = t1 * mrSges(7,1) - t2 * mrSges(7,2);
t122 = -pkin(3) - pkin(4);
t149 = Ifges(6,5) * t67;
t148 = Ifges(6,6) * t127;
t147 = pkin(3) * t117;
t78 = -qJ(4) * t114 + t115 * t122;
t73 = -pkin(5) + t78;
t79 = qJ(4) * t115 + t114 * t122;
t38 = -t116 * t79 + t119 * t73;
t145 = t38 * mrSges(7,1);
t39 = t116 * t73 + t119 * t79;
t144 = t39 * mrSges(7,2);
t143 = pkin(8) - qJ(5);
t82 = t143 * t117;
t85 = t143 * t120;
t45 = t114 * t82 + t115 * t85;
t136 = t118 * t120;
t76 = t121 * mrSges(5,1) + mrSges(5,2) * t136;
t142 = Ifges(4,4) * t117;
t141 = Ifges(4,4) * t120;
t140 = Ifges(5,5) * t117;
t139 = Ifges(5,5) * t120;
t138 = t121 * Ifges(5,6);
t135 = t154 * pkin(8) ^ 2;
t134 = Ifges(5,2) + Ifges(4,3) + Ifges(6,3);
t80 = -t120 * pkin(3) - t117 * qJ(4) - pkin(2);
t23 = -t116 * t60 + t119 * t59;
t24 = t116 * t59 + t119 * t60;
t133 = Ifges(7,5) * t24 + Ifges(7,6) * t23 + Ifges(7,3) * t121;
t27 = -t59 * mrSges(6,1) + t60 * mrSges(6,2);
t40 = mrSges(6,1) * t127 + t67 * mrSges(6,2);
t7 = -t23 * mrSges(7,1) + t24 * mrSges(7,2);
t35 = -t116 * t67 - t119 * t127;
t36 = -t116 * t127 + t119 * t67;
t12 = -t35 * mrSges(7,1) + t36 * mrSges(7,2);
t44 = -t114 * t85 + t115 * t82;
t51 = t120 * t81 - t95;
t64 = t120 * pkin(4) - t80;
t66 = -t114 * t116 + t115 * t119;
t68 = t114 * t119 + t115 * t116;
t131 = t66 * mrSges(7,1) - t68 * mrSges(7,2);
t130 = t117 * mrSges(4,1) + t120 * mrSges(4,2);
t129 = t117 * mrSges(5,1) - t120 * mrSges(5,3);
t128 = qJ(4) * t120 - t147;
t33 = Ifges(7,6) * t35;
t34 = Ifges(7,5) * t36;
t25 = -pkin(9) * t67 + t44;
t26 = -pkin(9) * t127 + t45;
t8 = -t116 * t26 + t119 * t25;
t9 = t116 * t25 + t119 * t26;
t126 = t8 * mrSges(7,1) - t9 * mrSges(7,2) + t33 + t34;
t90 = qJ(4) * t136;
t46 = t90 + (t122 * t117 - pkin(7)) * t118;
t125 = Ifges(6,5) * t60 - Ifges(5,6) * t137 + Ifges(6,6) * t59 + t133 + (-Ifges(5,4) - Ifges(4,5)) * t136;
t124 = pkin(7) ^ 2;
t113 = t121 ^ 2;
t111 = t118 ^ 2;
t105 = t111 * t124;
t103 = Ifges(5,4) * t117;
t102 = Ifges(4,5) * t117;
t101 = Ifges(4,6) * t120;
t89 = Ifges(4,1) * t117 + t141;
t88 = Ifges(5,1) * t117 - t139;
t87 = Ifges(4,2) * t120 + t142;
t86 = -Ifges(5,3) * t120 + t140;
t84 = -mrSges(4,1) * t120 + mrSges(4,2) * t117;
t83 = -mrSges(5,1) * t120 - mrSges(5,3) * t117;
t77 = -mrSges(5,2) * t137 - mrSges(5,3) * t121;
t75 = -mrSges(4,1) * t121 - mrSges(4,3) * t136;
t74 = mrSges(4,2) * t121 - mrSges(4,3) * t137;
t63 = t130 * t118;
t62 = t129 * t118;
t58 = -t90 + (pkin(7) + t147) * t118;
t57 = -Ifges(4,5) * t121 + (Ifges(4,1) * t120 - t142) * t118;
t56 = -Ifges(5,4) * t121 + (Ifges(5,1) * t120 + t140) * t118;
t55 = -Ifges(4,6) * t121 + (-Ifges(4,2) * t117 + t141) * t118;
t54 = -t138 + (Ifges(5,3) * t117 + t139) * t118;
t50 = t109 - t51;
t49 = mrSges(6,1) * t121 - mrSges(6,3) * t60;
t48 = -mrSges(6,2) * t121 + mrSges(6,3) * t59;
t43 = pkin(5) * t127 + t64;
t42 = Ifges(6,1) * t67 - Ifges(6,4) * t127;
t41 = Ifges(6,4) * t67 - Ifges(6,2) * t127;
t22 = Ifges(6,1) * t60 + Ifges(6,4) * t59 + t121 * Ifges(6,5);
t21 = Ifges(6,4) * t60 + Ifges(6,2) * t59 + t121 * Ifges(6,6);
t20 = -pkin(5) * t59 + t46;
t16 = mrSges(7,1) * t121 - mrSges(7,3) * t24;
t15 = -mrSges(7,2) * t121 + mrSges(7,3) * t23;
t14 = Ifges(7,1) * t36 + Ifges(7,4) * t35;
t13 = Ifges(7,4) * t36 + Ifges(7,2) * t35;
t5 = Ifges(7,1) * t24 + Ifges(7,4) * t23 + Ifges(7,5) * t121;
t4 = Ifges(7,4) * t24 + Ifges(7,2) * t23 + Ifges(7,6) * t121;
t17 = [0.2e1 * t1 * t16 + 0.2e1 * t10 * t49 + 0.2e1 * t11 * t48 + 0.2e1 * t2 * t15 + 0.2e1 * t20 * t7 + t59 * t21 + t60 * t22 + t23 * t4 + t24 * t5 + 0.2e1 * t46 * t27 + 0.2e1 * t47 * t77 + 0.2e1 * t50 * t76 + 0.2e1 * t51 * t75 + 0.2e1 * t52 * t74 + 0.2e1 * t58 * t62 + Ifges(2,3) + (t111 + t113) * mrSges(3,3) * t153 + m(3) * (pkin(1) ^ 2 + t113 * t124 + t105) + m(4) * (t51 ^ 2 + t52 ^ 2 + t105) + m(5) * (t47 ^ 2 + t50 ^ 2 + t58 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t46 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t20 ^ 2) + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t118 + t63 * t153 + (t56 + t57) * t120 + (t54 - t55) * t117) * t118 + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t134) * t121 + (Ifges(4,6) * t117 + (2 * Ifges(3,4))) * t118 + t125) * t121; m(6) * (t10 * t44 + t11 * t45 + t46 * t64) + m(7) * (t1 * t8 + t2 * t9 + t20 * t43) + (-t103 / 0.2e1 - t102 / 0.2e1 - t101 / 0.2e1 + t149 / 0.2e1 - t148 / 0.2e1 + t34 / 0.2e1 + t33 / 0.2e1 + Ifges(3,6)) * t121 + m(4) * (-pkin(2) * pkin(7) * t118 + (-t51 * t117 + t52 * t120) * pkin(8)) + m(5) * (t58 * t80 + (t50 * t117 + t47 * t120) * pkin(8)) + (-t10 * t67 - t11 * t127) * mrSges(6,3) - t127 * t21 / 0.2e1 + (-t1 * t36 + t2 * t35) * mrSges(7,3) + (t50 * mrSges(5,2) - t51 * mrSges(4,3) + t56 / 0.2e1 + t57 / 0.2e1 + (t86 / 0.2e1 - t87 / 0.2e1) * t118 + (-t75 + t76) * pkin(8)) * t117 + (-t121 * mrSges(3,2) + (-mrSges(3,1) + t84) * t118) * pkin(7) + Ifges(3,5) * t118 + t80 * t62 + t58 * t83 - pkin(2) * t63 + t64 * t27 + t67 * t22 / 0.2e1 + t59 * t41 / 0.2e1 + t60 * t42 / 0.2e1 + t43 * t7 + t46 * t40 + t45 * t48 + t44 * t49 + t35 * t4 / 0.2e1 + t36 * t5 / 0.2e1 + t8 * t16 + t20 * t12 + t23 * t13 / 0.2e1 + t24 * t14 / 0.2e1 + t9 * t15 + (t47 * mrSges(5,2) + t52 * mrSges(4,3) - t54 / 0.2e1 + t55 / 0.2e1 + t138 / 0.2e1 + (t88 / 0.2e1 + t89 / 0.2e1) * t118 + (t74 + t77) * pkin(8)) * t120; -0.2e1 * pkin(2) * t84 + 0.2e1 * t43 * t12 + t35 * t13 + t36 * t14 + 0.2e1 * t64 * t40 - t127 * t41 + t67 * t42 + 0.2e1 * t80 * t83 + Ifges(3,3) + (t87 - t86) * t120 + (t88 + t89) * t117 + m(5) * (t80 ^ 2 + t135) + m(4) * (pkin(2) ^ 2 + t135) + m(6) * (t44 ^ 2 + t45 ^ 2 + t64 ^ 2) + m(7) * (t43 ^ 2 + t8 ^ 2 + t9 ^ 2) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * pkin(8) * t154 + 0.2e1 * (t35 * t9 - t36 * t8) * mrSges(7,3) + 0.2e1 * (-t127 * t45 - t44 * t67) * mrSges(6,3); m(6) * (t10 * t78 + t11 * t79) + m(5) * (-pkin(3) * t50 + qJ(4) * t47) + m(7) * (t1 * t38 + t2 * t39) - t134 * t121 - Ifges(4,6) * t137 - pkin(3) * t76 + qJ(4) * t77 + t78 * t49 + t79 * t48 + t51 * mrSges(4,1) - t52 * mrSges(4,2) + t47 * mrSges(5,3) - t50 * mrSges(5,1) + t38 * t16 + t39 * t15 - t10 * mrSges(6,1) + t11 * mrSges(6,2) - t125 - t152; -t44 * mrSges(6,1) + t45 * mrSges(6,2) - t149 - Ifges(5,6) * t120 + t148 + t101 + t102 + t103 + (t35 * t39 - t36 * t38) * mrSges(7,3) + (-t127 * t79 - t67 * t78) * mrSges(6,3) + t128 * mrSges(5,2) + m(6) * (t44 * t78 + t45 * t79) + m(7) * (t38 * t8 + t39 * t9) + (m(5) * t128 - t129 - t130) * pkin(8) - t126; 0.2e1 * pkin(3) * mrSges(5,1) - 0.2e1 * t78 * mrSges(6,1) - 0.2e1 * t145 + 0.2e1 * t79 * mrSges(6,2) + 0.2e1 * t144 + 0.2e1 * qJ(4) * mrSges(5,3) + Ifges(7,3) + m(6) * (t78 ^ 2 + t79 ^ 2) + m(7) * (t38 ^ 2 + t39 ^ 2) + m(5) * (pkin(3) ^ 2 + qJ(4) ^ 2) + t134; t114 * t48 + t115 * t49 + t68 * t15 + t66 * t16 + m(7) * (t1 * t66 + t2 * t68) + m(6) * (t10 * t115 + t11 * t114) + m(5) * t50 + t76; (m(5) * pkin(8) + mrSges(5,2)) * t117 + (t35 * t68 - t36 * t66) * mrSges(7,3) + (-t114 * t127 - t115 * t67) * mrSges(6,3) + m(7) * (t66 * t8 + t68 * t9) + m(6) * (t114 * t45 + t115 * t44); -m(5) * pkin(3) - t115 * mrSges(6,1) + t114 * mrSges(6,2) - mrSges(5,1) + m(6) * (t114 * t79 + t115 * t78) + m(7) * (t38 * t66 + t39 * t68) - t131; m(5) + m(6) * (t114 ^ 2 + t115 ^ 2) + m(7) * (t66 ^ 2 + t68 ^ 2); m(6) * t46 + m(7) * t20 + t27 + t7; m(6) * t64 + m(7) * t43 + t12 + t40; 0; 0; m(6) + m(7); t133 + t152; t126; -Ifges(7,3) - t144 + t145; t131; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t17(1) t17(2) t17(4) t17(7) t17(11) t17(16); t17(2) t17(3) t17(5) t17(8) t17(12) t17(17); t17(4) t17(5) t17(6) t17(9) t17(13) t17(18); t17(7) t17(8) t17(9) t17(10) t17(14) t17(19); t17(11) t17(12) t17(13) t17(14) t17(15) t17(20); t17(16) t17(17) t17(18) t17(19) t17(20) t17(21);];
Mq  = res;
