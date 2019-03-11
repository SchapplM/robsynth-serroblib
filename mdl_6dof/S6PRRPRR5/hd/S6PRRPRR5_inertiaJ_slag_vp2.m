% Calculate joint inertia matrix for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:26
% EndTime: 2019-03-08 22:16:28
% DurationCPUTime: 0.97s
% Computational Cost: add. (1686->282), mult. (3695->414), div. (0->0), fcn. (4013->12), ass. (0->116)
t144 = 2 * pkin(8);
t108 = sin(qJ(6));
t112 = cos(qJ(6));
t104 = sin(pkin(12));
t106 = cos(pkin(12));
t109 = sin(qJ(5));
t113 = cos(qJ(5));
t77 = -t104 * t109 + t106 * t113;
t78 = t104 * t113 + t106 * t109;
t39 = -t108 * t78 + t112 * t77;
t40 = t108 * t77 + t112 * t78;
t15 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t44 = -t77 * mrSges(6,1) + t78 * mrSges(6,2);
t143 = t15 + t44;
t110 = sin(qJ(3));
t62 = t78 * t110;
t63 = t77 * t110;
t142 = -Ifges(6,5) * t63 + Ifges(6,6) * t62;
t29 = -t108 * t63 - t112 * t62;
t30 = -t108 * t62 + t112 * t63;
t11 = -t29 * mrSges(7,1) + t30 * mrSges(7,2);
t33 = t62 * mrSges(6,1) + t63 * mrSges(6,2);
t125 = t106 * t110;
t128 = t104 * t110;
t68 = mrSges(5,1) * t128 + mrSges(5,2) * t125;
t141 = t11 + t33 + t68;
t107 = cos(pkin(6));
t114 = cos(qJ(3));
t105 = sin(pkin(6));
t111 = sin(qJ(2));
t127 = t105 * t111;
t65 = -t107 * t114 + t110 * t127;
t64 = t65 ^ 2;
t140 = -t104 / 0.2e1;
t139 = t106 / 0.2e1;
t138 = pkin(8) * t114;
t98 = t110 * pkin(8);
t137 = -Ifges(7,3) - Ifges(6,3);
t136 = pkin(9) + qJ(4);
t135 = -Ifges(7,5) * t30 - Ifges(7,6) * t29;
t82 = -pkin(3) * t114 - qJ(4) * t110 - pkin(2);
t73 = t106 * t82;
t41 = -pkin(9) * t125 + t73 + (-pkin(8) * t104 - pkin(4)) * t114;
t55 = t104 * t82 + t106 * t138;
t49 = -pkin(9) * t128 + t55;
t19 = t109 * t41 + t113 * t49;
t83 = t136 * t104;
t85 = t136 * t106;
t51 = -t109 * t83 + t113 * t85;
t84 = -t106 * mrSges(5,1) + t104 * mrSges(5,2);
t134 = t84 - mrSges(4,1);
t81 = pkin(4) * t128 + t98;
t133 = Ifges(5,4) * t104;
t132 = Ifges(5,4) * t106;
t131 = t110 * t65;
t67 = t107 * t110 + t114 * t127;
t130 = t114 * t67;
t129 = t104 ^ 2 + t106 ^ 2;
t115 = cos(qJ(2));
t126 = t105 * t115;
t47 = -t104 * t67 - t106 * t126;
t48 = -t104 * t126 + t106 * t67;
t20 = -t109 * t48 + t113 * t47;
t21 = t109 * t47 + t113 * t48;
t5 = -t108 * t21 + t112 * t20;
t6 = t108 * t20 + t112 * t21;
t124 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t95 = -pkin(4) * t106 - pkin(3);
t18 = -t109 * t49 + t113 * t41;
t50 = -t109 * t85 - t113 * t83;
t10 = -pkin(5) * t114 - pkin(10) * t63 + t18;
t14 = -pkin(10) * t62 + t19;
t2 = t10 * t112 - t108 * t14;
t3 = t10 * t108 + t112 * t14;
t123 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t135;
t31 = -pkin(10) * t78 + t50;
t32 = pkin(10) * t77 + t51;
t12 = -t108 * t32 + t112 * t31;
t13 = t108 * t31 + t112 * t32;
t37 = Ifges(7,6) * t39;
t38 = Ifges(7,5) * t40;
t122 = t12 * mrSges(7,1) - t13 * mrSges(7,2) + t37 + t38;
t121 = -t104 * t47 + t106 * t48;
t54 = -t104 * t138 + t73;
t120 = -t54 * t104 + t55 * t106;
t119 = (mrSges(7,1) * t112 - mrSges(7,2) * t108) * pkin(5);
t117 = pkin(8) ^ 2;
t103 = t114 ^ 2;
t102 = t110 ^ 2;
t100 = t105 ^ 2;
t97 = t102 * t117;
t93 = t100 * t115 ^ 2;
t88 = -mrSges(4,1) * t114 + mrSges(4,2) * t110;
t87 = Ifges(5,1) * t104 + t132;
t86 = Ifges(5,2) * t106 + t133;
t80 = -mrSges(5,1) * t114 - mrSges(5,3) * t125;
t79 = mrSges(5,2) * t114 - mrSges(5,3) * t128;
t71 = Ifges(6,5) * t78;
t70 = Ifges(6,6) * t77;
t61 = -Ifges(5,5) * t114 + (Ifges(5,1) * t106 - t133) * t110;
t60 = -Ifges(5,6) * t114 + (-Ifges(5,2) * t104 + t132) * t110;
t56 = -pkin(5) * t77 + t95;
t53 = -mrSges(6,1) * t114 - mrSges(6,3) * t63;
t52 = mrSges(6,2) * t114 - mrSges(6,3) * t62;
t46 = Ifges(6,1) * t78 + Ifges(6,4) * t77;
t45 = Ifges(6,4) * t78 + Ifges(6,2) * t77;
t43 = pkin(5) * t62 + t81;
t28 = Ifges(6,1) * t63 - Ifges(6,4) * t62 - Ifges(6,5) * t114;
t27 = Ifges(6,4) * t63 - Ifges(6,2) * t62 - Ifges(6,6) * t114;
t23 = -mrSges(7,1) * t114 - mrSges(7,3) * t30;
t22 = mrSges(7,2) * t114 + mrSges(7,3) * t29;
t17 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t16 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t29 - Ifges(7,5) * t114;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t29 - Ifges(7,6) * t114;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t64) + m(6) * (t20 ^ 2 + t21 ^ 2 + t64) + m(5) * (t47 ^ 2 + t48 ^ 2 + t64) + m(4) * (t67 ^ 2 + t64 + t93) + m(3) * (t100 * t111 ^ 2 + t107 ^ 2 + t93); mrSges(4,3) * t130 + t20 * t53 + t21 * t52 + t6 * t22 + t5 * t23 + t47 * t80 + t48 * t79 + (-t111 * mrSges(3,2) + (mrSges(3,1) - t88) * t115) * t105 + (t110 * mrSges(4,3) + t141) * t65 + m(7) * (t2 * t5 + t3 * t6 + t43 * t65) + m(6) * (t18 * t20 + t19 * t21 + t65 * t81) + m(5) * (pkin(8) * t131 + t47 * t54 + t48 * t55) + m(4) * (pkin(2) * t126 + (t130 + t131) * pkin(8)); -0.2e1 * pkin(2) * t88 + 0.2e1 * t43 * t11 + 0.2e1 * t18 * t53 + 0.2e1 * t19 * t52 + 0.2e1 * t2 * t23 + 0.2e1 * t3 * t22 - t62 * t27 + t63 * t28 + t29 * t7 + t30 * t8 + 0.2e1 * t81 * t33 + 0.2e1 * t54 * t80 + 0.2e1 * t55 * t79 + Ifges(3,3) + (t102 + t103) * mrSges(4,3) * t144 + (Ifges(4,1) * t110 - t104 * t60 + t106 * t61 + t68 * t144) * t110 + m(4) * (pkin(2) ^ 2 + t103 * t117 + t97) + m(5) * (t54 ^ 2 + t55 ^ 2 + t97) + m(6) * (t18 ^ 2 + t19 ^ 2 + t81 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t43 ^ 2) + ((Ifges(5,3) + Ifges(4,2) - t137) * t114 + (-Ifges(5,5) * t106 + Ifges(5,6) * t104 + (2 * Ifges(4,4))) * t110 + t135 + t142) * t114; -t67 * mrSges(4,2) + (t39 * t6 - t40 * t5) * mrSges(7,3) + (-t20 * t78 + t21 * t77) * mrSges(6,3) + t121 * mrSges(5,3) + (t134 + t143) * t65 + m(7) * (t12 * t5 + t13 * t6 + t56 * t65) + m(6) * (t20 * t50 + t21 * t51 + t65 * t95) + m(5) * (-pkin(3) * t65 + qJ(4) * t121); (-t2 * t40 + t3 * t39) * mrSges(7,3) + (-t18 * t78 + t19 * t77) * mrSges(6,3) + m(5) * (-pkin(3) * t98 + qJ(4) * t120) + (pkin(8) * t134 + t139 * t87 + t140 * t86 + Ifges(4,5)) * t110 + (Ifges(5,5) * t140 - Ifges(5,6) * t106 / 0.2e1 - t71 / 0.2e1 - t70 / 0.2e1 - t38 / 0.2e1 - t37 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t114 + t120 * mrSges(5,3) + t104 * t61 / 0.2e1 + t77 * t27 / 0.2e1 + t78 * t28 / 0.2e1 + t81 * t44 + t95 * t33 - pkin(3) * t68 + t50 * t53 + t56 * t11 - t62 * t45 / 0.2e1 + t63 * t46 / 0.2e1 + t39 * t7 / 0.2e1 + t40 * t8 / 0.2e1 + t43 * t15 + t51 * t52 + t29 * t16 / 0.2e1 + t30 * t17 / 0.2e1 + t13 * t22 + t12 * t23 + t60 * t139 + m(6) * (t18 * t50 + t19 * t51 + t81 * t95) + m(7) * (t12 * t2 + t13 * t3 + t43 * t56) + (-t104 * t80 + t106 * t79) * qJ(4); -0.2e1 * pkin(3) * t84 + t104 * t87 + t106 * t86 + 0.2e1 * t56 * t15 + t39 * t16 + t40 * t17 + 0.2e1 * t95 * t44 + t77 * t45 + t78 * t46 + Ifges(4,3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t56 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2 + t95 ^ 2) + m(5) * (qJ(4) ^ 2 * t129 + pkin(3) ^ 2) + 0.2e1 * (-t12 * t40 + t13 * t39) * mrSges(7,3) + 0.2e1 * (-t50 * t78 + t51 * t77) * mrSges(6,3) + 0.2e1 * t129 * qJ(4) * mrSges(5,3); 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1 + m(5) / 0.2e1) * t65; m(5) * t98 + m(6) * t81 + m(7) * t43 + t141; -m(5) * pkin(3) + m(6) * t95 + m(7) * t56 + t143 + t84; m(5) + m(6) + m(7); t20 * mrSges(6,1) - t21 * mrSges(6,2) + m(7) * (t108 * t6 + t112 * t5) * pkin(5) + t124; t18 * mrSges(6,1) - t19 * mrSges(6,2) + t137 * t114 + (t108 * t22 + m(7) * (t108 * t3 + t112 * t2) + t112 * t23) * pkin(5) + t123 - t142; t50 * mrSges(6,1) - t51 * mrSges(6,2) + t70 + t71 + (m(7) * (t108 * t13 + t112 * t12) + (t108 * t39 - t112 * t40) * mrSges(7,3)) * pkin(5) + t122; 0; m(7) * (t108 ^ 2 + t112 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t119 - t137; t124; -Ifges(7,3) * t114 + t123; t122; 0; Ifges(7,3) + t119; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
