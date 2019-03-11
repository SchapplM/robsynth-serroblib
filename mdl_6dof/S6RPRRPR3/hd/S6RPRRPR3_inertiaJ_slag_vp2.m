% Calculate joint inertia matrix for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:43
% EndTime: 2019-03-09 05:04:45
% DurationCPUTime: 1.06s
% Computational Cost: add. (1024->278), mult. (1936->370), div. (0->0), fcn. (1676->8), ass. (0->107)
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t104 = t82 ^ 2 + t85 ^ 2;
t79 = sin(pkin(10));
t65 = pkin(1) * t79 + pkin(7);
t129 = 0.2e1 * t65;
t121 = pkin(4) * t82;
t102 = qJ(5) * t85;
t83 = sin(qJ(3));
t57 = t83 * t102;
t20 = -t57 + (t65 + t121) * t83;
t97 = t82 * mrSges(6,1) - t85 * mrSges(6,3);
t31 = t97 * t83;
t128 = m(6) * t20 + t31;
t109 = t83 * t85;
t110 = t82 * t83;
t127 = -Ifges(6,6) * t110 + (-Ifges(6,4) - Ifges(5,5)) * t109;
t126 = (mrSges(6,2) + mrSges(5,3)) * t104;
t125 = pkin(4) + pkin(5);
t124 = pkin(8) - pkin(9);
t86 = cos(qJ(3));
t122 = pkin(3) * t86;
t120 = pkin(8) * t83;
t119 = Ifges(5,4) * t82;
t118 = Ifges(5,4) * t85;
t117 = Ifges(6,5) * t82;
t116 = Ifges(6,5) * t85;
t115 = Ifges(5,6) * t86;
t114 = Ifges(6,6) * t86;
t81 = sin(qJ(6));
t84 = cos(qJ(6));
t45 = -qJ(5) * t81 - t125 * t84;
t113 = t45 * mrSges(7,1);
t46 = qJ(5) * t84 - t125 * t81;
t112 = t46 * mrSges(7,2);
t111 = t65 * t86;
t49 = -mrSges(5,1) * t85 + mrSges(5,2) * t82;
t108 = mrSges(4,1) - t49;
t107 = Ifges(6,2) + Ifges(5,3);
t80 = cos(pkin(10));
t66 = -pkin(1) * t80 - pkin(2);
t33 = -t120 + t66 - t122;
t12 = t85 * t111 + t82 * t33;
t41 = t86 * mrSges(6,1) + mrSges(6,2) * t109;
t106 = t104 * t120;
t105 = t104 * pkin(8) ^ 2;
t76 = t83 ^ 2;
t78 = t86 ^ 2;
t103 = t76 + t78;
t38 = -t81 * t85 + t82 * t84;
t27 = t38 * t83;
t93 = t81 * t82 + t84 * t85;
t28 = t93 * t83;
t101 = Ifges(7,5) * t28 + Ifges(7,6) * t27 + Ifges(7,3) * t86;
t100 = qJ(5) * t82 + pkin(3);
t43 = t82 * t111;
t11 = t33 * t85 - t43;
t8 = -qJ(5) * t86 + t12;
t98 = t82 * mrSges(5,1) + t85 * mrSges(5,2);
t7 = -t27 * mrSges(7,1) + t28 * mrSges(7,2);
t96 = t84 * mrSges(7,1) - t81 * mrSges(7,2);
t95 = t102 - t121;
t94 = -t11 * t82 + t12 * t85;
t55 = t124 * t82;
t56 = t124 * t85;
t16 = t55 * t84 - t56 * t81;
t17 = t55 * t81 + t56 * t84;
t35 = Ifges(7,6) * t93;
t36 = Ifges(7,5) * t38;
t92 = t16 * mrSges(7,1) - t17 * mrSges(7,2) - t35 + t36;
t74 = t86 * pkin(4);
t3 = pkin(5) * t86 + t43 + t74 + (-pkin(9) * t83 - t33) * t85;
t6 = pkin(9) * t110 + t8;
t1 = t3 * t84 - t6 * t81;
t2 = t3 * t81 + t6 * t84;
t91 = t1 * mrSges(7,1) - t2 * mrSges(7,2) + t101;
t39 = mrSges(5,2) * t86 - mrSges(5,3) * t110;
t40 = -mrSges(5,1) * t86 - mrSges(5,3) * t109;
t42 = -mrSges(6,2) * t110 - mrSges(6,3) * t86;
t9 = -t11 + t74;
t90 = m(6) * (t8 * t85 + t82 * t9) + (t39 + t42) * t85 + (-t40 + t41) * t82;
t71 = Ifges(6,4) * t82;
t70 = Ifges(5,5) * t82;
t69 = Ifges(5,6) * t85;
t64 = t65 ^ 2;
t54 = t76 * t64;
t53 = Ifges(5,1) * t82 + t118;
t52 = Ifges(6,1) * t82 - t116;
t51 = Ifges(5,2) * t85 + t119;
t50 = -Ifges(6,3) * t85 + t117;
t48 = -mrSges(6,1) * t85 - mrSges(6,3) * t82;
t47 = -pkin(4) * t85 - t100;
t34 = t125 * t85 + t100;
t32 = t98 * t83;
t26 = -Ifges(5,5) * t86 + (Ifges(5,1) * t85 - t119) * t83;
t25 = -Ifges(6,4) * t86 + (Ifges(6,1) * t85 + t117) * t83;
t24 = -t115 + (-Ifges(5,2) * t82 + t118) * t83;
t23 = -t114 + (Ifges(6,3) * t82 + t116) * t83;
t19 = mrSges(7,1) * t86 - mrSges(7,3) * t28;
t18 = -mrSges(7,2) * t86 + mrSges(7,3) * t27;
t15 = Ifges(7,1) * t38 - Ifges(7,4) * t93;
t14 = Ifges(7,4) * t38 - Ifges(7,2) * t93;
t13 = mrSges(7,1) * t93 + mrSges(7,2) * t38;
t10 = t57 + (-t125 * t82 - t65) * t83;
t5 = Ifges(7,1) * t28 + Ifges(7,4) * t27 + Ifges(7,5) * t86;
t4 = Ifges(7,4) * t28 + Ifges(7,2) * t27 + Ifges(7,6) * t86;
t21 = [0.2e1 * t1 * t19 + 0.2e1 * t10 * t7 + 0.2e1 * t11 * t40 + 0.2e1 * t12 * t39 + 0.2e1 * t2 * t18 + 0.2e1 * t20 * t31 + t27 * t4 + t28 * t5 + 0.2e1 * t9 * t41 + 0.2e1 * t8 * t42 + Ifges(2,3) + Ifges(3,3) + (-0.2e1 * t66 * mrSges(4,1) + (Ifges(4,2) + t107) * t86 + t101 + t127) * t86 + (0.2e1 * t66 * mrSges(4,2) + Ifges(4,1) * t83 + 0.2e1 * Ifges(4,4) * t86 + t32 * t129 + (t25 + t26) * t85 + (t23 - t24 + t115) * t82) * t83 + m(4) * (t64 * t78 + t66 ^ 2 + t54) + m(5) * (t11 ^ 2 + t12 ^ 2 + t54) + m(6) * (t20 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(7) * (t1 ^ 2 + t10 ^ 2 + t2 ^ 2) + m(3) * (t79 ^ 2 + t80 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t80 - mrSges(3,2) * t79) * pkin(1) + t103 * mrSges(4,3) * t129; t28 * t18 + t27 * t19 + m(7) * (t1 * t27 + t2 * t28) + (m(5) * (t94 - t111) + t90) * t83 + (m(7) * t10 - t128 - t32 + t7) * t86; m(3) + m(7) * (t27 ^ 2 + t28 ^ 2 + t78) + m(4) * t103 + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (t104 * t76 + t78); -t93 * t4 / 0.2e1 + t38 * t5 / 0.2e1 + t20 * t48 + t28 * t15 / 0.2e1 - pkin(3) * t32 + t34 * t7 + t17 * t18 + t16 * t19 + t27 * t14 / 0.2e1 + t10 * t13 + (-t1 * t38 - t2 * t93) * mrSges(7,3) + m(7) * (t1 * t16 + t10 * t34 + t17 * t2) + (-t71 / 0.2e1 - t70 / 0.2e1 - t69 / 0.2e1 + t36 / 0.2e1 - t35 / 0.2e1 + Ifges(4,6) - t65 * mrSges(4,2)) * t86 + (-t23 / 0.2e1 + t24 / 0.2e1 + t114 / 0.2e1 + t12 * mrSges(5,3) + t8 * mrSges(6,2)) * t85 + (t25 / 0.2e1 + t26 / 0.2e1 + t9 * mrSges(6,2) - t11 * mrSges(5,3)) * t82 + (m(5) * t94 + t90) * pkin(8) + (Ifges(4,5) + (t52 / 0.2e1 + t53 / 0.2e1) * t85 + (t50 / 0.2e1 - t51 / 0.2e1) * t82 + (-m(5) * pkin(3) - t108) * t65) * t83 + t128 * t47; (-t27 * t38 - t28 * t93) * mrSges(7,3) + (t13 - t48 + t108) * t86 + m(6) * (-t47 * t86 + t106) + m(7) * (t16 * t27 + t17 * t28 + t34 * t86) + m(5) * (t106 + t122) + (-mrSges(4,2) + t126) * t83; -0.2e1 * pkin(3) * t49 + 0.2e1 * t34 * t13 - t93 * t14 + t38 * t15 + 0.2e1 * t47 * t48 + Ifges(4,3) + (-t50 + t51) * t85 + (t52 + t53) * t82 + 0.2e1 * (-t16 * t38 - t17 * t93) * mrSges(7,3) + m(7) * (t16 ^ 2 + t17 ^ 2 + t34 ^ 2) + m(6) * (t47 ^ 2 + t105) + m(5) * (pkin(3) ^ 2 + t105) + 0.2e1 * pkin(8) * t126; m(7) * (t1 * t45 + t2 * t46) + m(6) * (-pkin(4) * t9 + qJ(5) * t8) - t107 * t86 - t91 - pkin(4) * t41 + qJ(5) * t42 + t45 * t19 + t46 * t18 + t8 * mrSges(6,3) - t9 * mrSges(6,1) + t11 * mrSges(5,1) - t12 * mrSges(5,2) - Ifges(5,6) * t110 - t127; ((-mrSges(5,2) + mrSges(6,3)) * t85 + (-mrSges(5,1) - mrSges(6,1)) * t82) * t83 + m(6) * (-pkin(4) * t110 + t57) + m(7) * (t27 * t45 + t28 * t46) + t7; m(7) * (t16 * t45 + t17 * t46) + t71 - Ifges(6,6) * t85 + t70 + t69 + (-t38 * t45 - t46 * t93) * mrSges(7,3) + t95 * mrSges(6,2) + (m(6) * t95 - t97 - t98) * pkin(8) - t92; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t113 + 0.2e1 * t112 + 0.2e1 * qJ(5) * mrSges(6,3) + Ifges(7,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + m(7) * (t45 ^ 2 + t46 ^ 2) + t107; t81 * t18 + t84 * t19 + m(7) * (t1 * t84 + t2 * t81) + m(6) * t9 + t41; m(6) * t110 + m(7) * (t27 * t84 + t28 * t81); m(7) * (t16 * t84 + t17 * t81) + (m(6) * pkin(8) + mrSges(6,2)) * t82 + (-t38 * t84 - t81 * t93) * mrSges(7,3); -mrSges(6,1) - m(6) * pkin(4) + m(7) * (t45 * t84 + t46 * t81) - t96; m(6) + m(7) * (t81 ^ 2 + t84 ^ 2); t91; -t7; t92; -Ifges(7,3) - t112 + t113; t96; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
