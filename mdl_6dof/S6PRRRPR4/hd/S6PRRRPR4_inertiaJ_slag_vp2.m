% Calculate joint inertia matrix for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:43
% EndTime: 2019-03-08 23:15:45
% DurationCPUTime: 1.05s
% Computational Cost: add. (1787->304), mult. (3902->446), div. (0->0), fcn. (4208->12), ass. (0->122)
t148 = 2 * pkin(8);
t147 = m(6) * pkin(4);
t113 = sin(qJ(3));
t107 = sin(pkin(12));
t109 = cos(pkin(12));
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t75 = t107 * t116 + t109 * t112;
t62 = t75 * t113;
t74 = -t107 * t112 + t109 * t116;
t63 = t74 * t113;
t33 = t62 * mrSges(6,1) + mrSges(6,2) * t63;
t111 = sin(qJ(6));
t115 = cos(qJ(6));
t29 = -t111 * t63 - t115 * t62;
t30 = -t111 * t62 + t115 * t63;
t9 = -t29 * mrSges(7,1) + mrSges(7,2) * t30;
t146 = t33 + t9;
t39 = -t111 * t75 + t115 * t74;
t40 = t111 * t74 + t115 * t75;
t15 = -t39 * mrSges(7,1) + mrSges(7,2) * t40;
t44 = -t74 * mrSges(6,1) + mrSges(6,2) * t75;
t145 = t15 + t44;
t129 = t113 * t116;
t144 = -Ifges(5,5) * t129 - Ifges(6,5) * t63 + Ifges(6,6) * t62;
t110 = cos(pkin(6));
t117 = cos(qJ(3));
t108 = sin(pkin(6));
t114 = sin(qJ(2));
t132 = t108 * t114;
t67 = -t110 * t117 + t113 * t132;
t66 = t67 ^ 2;
t143 = pkin(4) * t107;
t142 = pkin(8) * t117;
t101 = t113 * pkin(8);
t96 = pkin(4) * t109 + pkin(5);
t64 = -t111 * t143 + t115 * t96;
t141 = t64 * mrSges(7,1);
t65 = t111 * t96 + t115 * t143;
t140 = t65 * mrSges(7,2);
t139 = -qJ(5) - pkin(9);
t138 = -Ifges(7,5) * t30 - Ifges(7,6) * t29;
t84 = -pkin(3) * t117 - pkin(9) * t113 - pkin(2);
t77 = t116 * t84;
t41 = -qJ(5) * t129 + t77 + (-pkin(8) * t112 - pkin(4)) * t117;
t130 = t112 * t113;
t55 = t112 * t84 + t116 * t142;
t49 = -qJ(5) * t130 + t55;
t19 = t107 * t41 + t109 * t49;
t85 = t139 * t112;
t88 = t139 * t116;
t51 = t107 * t85 - t109 * t88;
t86 = -mrSges(5,1) * t116 + mrSges(5,2) * t112;
t137 = t86 - mrSges(4,1);
t136 = Ifges(5,4) * t112;
t135 = Ifges(5,4) * t116;
t134 = t113 * t67;
t69 = t110 * t113 + t117 * t132;
t133 = t117 * t69;
t83 = pkin(4) * t130 + t101;
t118 = cos(qJ(2));
t131 = t108 * t118;
t128 = t112 ^ 2 + t116 ^ 2;
t127 = -Ifges(7,3) - Ifges(6,3) - Ifges(5,3);
t47 = -t112 * t69 - t116 * t131;
t48 = -t112 * t131 + t116 * t69;
t20 = -t107 * t48 + t109 * t47;
t21 = t107 * t47 + t109 * t48;
t5 = -t111 * t21 + t115 * t20;
t6 = t111 * t20 + t115 * t21;
t126 = mrSges(7,1) * t5 - t6 * mrSges(7,2);
t97 = -pkin(4) * t116 - pkin(3);
t18 = -t107 * t49 + t109 * t41;
t50 = t107 * t88 + t109 * t85;
t11 = -pkin(5) * t117 - pkin(10) * t63 + t18;
t14 = -pkin(10) * t62 + t19;
t2 = t11 * t115 - t111 * t14;
t3 = t11 * t111 + t115 * t14;
t125 = mrSges(7,1) * t2 - t3 * mrSges(7,2) - t138;
t124 = mrSges(5,1) * t112 + mrSges(5,2) * t116;
t123 = -t112 * t47 + t116 * t48;
t31 = -pkin(10) * t75 + t50;
t32 = pkin(10) * t74 + t51;
t12 = -t111 * t32 + t115 * t31;
t13 = t111 * t31 + t115 * t32;
t37 = Ifges(7,6) * t39;
t38 = Ifges(7,5) * t40;
t122 = mrSges(7,1) * t12 - t13 * mrSges(7,2) + t37 + t38;
t120 = pkin(8) ^ 2;
t106 = t117 ^ 2;
t104 = t113 ^ 2;
t102 = t108 ^ 2;
t100 = t104 * t120;
t99 = Ifges(5,5) * t112;
t98 = Ifges(5,6) * t116;
t94 = t102 * t118 ^ 2;
t90 = Ifges(5,1) * t112 + t135;
t89 = Ifges(5,2) * t116 + t136;
t87 = -mrSges(4,1) * t117 + mrSges(4,2) * t113;
t82 = -mrSges(5,1) * t117 - mrSges(5,3) * t129;
t81 = mrSges(5,2) * t117 - mrSges(5,3) * t130;
t73 = t124 * t113;
t72 = Ifges(6,5) * t75;
t71 = Ifges(6,6) * t74;
t61 = -Ifges(5,5) * t117 + (Ifges(5,1) * t116 - t136) * t113;
t60 = -Ifges(5,6) * t117 + (-Ifges(5,2) * t112 + t135) * t113;
t56 = -pkin(5) * t74 + t97;
t54 = -t112 * t142 + t77;
t53 = -mrSges(6,1) * t117 - mrSges(6,3) * t63;
t52 = mrSges(6,2) * t117 - mrSges(6,3) * t62;
t46 = Ifges(6,1) * t75 + Ifges(6,4) * t74;
t45 = Ifges(6,4) * t75 + Ifges(6,2) * t74;
t42 = pkin(5) * t62 + t83;
t28 = Ifges(6,1) * t63 - Ifges(6,4) * t62 - Ifges(6,5) * t117;
t27 = Ifges(6,4) * t63 - Ifges(6,2) * t62 - Ifges(6,6) * t117;
t23 = -mrSges(7,1) * t117 - mrSges(7,3) * t30;
t22 = mrSges(7,2) * t117 + mrSges(7,3) * t29;
t17 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t16 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t8 = Ifges(7,1) * t30 + Ifges(7,4) * t29 - Ifges(7,5) * t117;
t7 = Ifges(7,4) * t30 + Ifges(7,2) * t29 - Ifges(7,6) * t117;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t66) + m(6) * (t20 ^ 2 + t21 ^ 2 + t66) + m(5) * (t47 ^ 2 + t48 ^ 2 + t66) + m(4) * (t69 ^ 2 + t66 + t94) + m(3) * (t102 * t114 ^ 2 + t110 ^ 2 + t94); mrSges(4,3) * t133 + t20 * t53 + t21 * t52 + t6 * t22 + t5 * t23 + t47 * t82 + t48 * t81 + (-t114 * mrSges(3,2) + (mrSges(3,1) - t87) * t118) * t108 + (t113 * mrSges(4,3) + t146 + t73) * t67 + m(7) * (t2 * t5 + t3 * t6 + t42 * t67) + m(6) * (t18 * t20 + t19 * t21 + t67 * t83) + m(5) * (pkin(8) * t134 + t47 * t54 + t48 * t55) + m(4) * (pkin(2) * t131 + (t133 + t134) * pkin(8)); -0.2e1 * pkin(2) * t87 + 0.2e1 * t18 * t53 + 0.2e1 * t19 * t52 + 0.2e1 * t2 * t23 + 0.2e1 * t3 * t22 - t62 * t27 + t63 * t28 + t29 * t7 + t30 * t8 + 0.2e1 * t83 * t33 + 0.2e1 * t42 * t9 + 0.2e1 * t54 * t82 + 0.2e1 * t55 * t81 + Ifges(3,3) + (t104 + t106) * mrSges(4,3) * t148 + (Ifges(4,1) * t113 - t112 * t60 + t116 * t61 + t148 * t73) * t113 + m(4) * (pkin(2) ^ 2 + t106 * t120 + t100) + m(5) * (t54 ^ 2 + t55 ^ 2 + t100) + m(6) * (t18 ^ 2 + t19 ^ 2 + t83 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t42 ^ 2) + ((Ifges(4,2) - t127) * t117 + (Ifges(5,6) * t112 + (2 * Ifges(4,4))) * t113 + t138 + t144) * t117; -t69 * mrSges(4,2) + (t39 * t6 - t40 * t5) * mrSges(7,3) + (-t20 * t75 + t21 * t74) * mrSges(6,3) + t123 * mrSges(5,3) + (t137 + t145) * t67 + m(7) * (t12 * t5 + t13 * t6 + t56 * t67) + m(6) * (t20 * t50 + t21 * t51 + t67 * t97) + m(5) * (-pkin(3) * t67 + pkin(9) * t123); (t60 / 0.2e1 + pkin(9) * t81 + t55 * mrSges(5,3)) * t116 + (t61 / 0.2e1 - pkin(9) * t82 - t54 * mrSges(5,3)) * t112 + t97 * t33 + t74 * t27 / 0.2e1 + t75 * t28 / 0.2e1 + t83 * t44 + t63 * t46 / 0.2e1 - pkin(3) * t73 + t56 * t9 - t62 * t45 / 0.2e1 + t39 * t7 / 0.2e1 + t40 * t8 / 0.2e1 + t42 * t15 + t51 * t52 + t50 * t53 + t29 * t16 / 0.2e1 + t30 * t17 / 0.2e1 + t13 * t22 + t12 * t23 + (-t18 * t75 + t19 * t74) * mrSges(6,3) + (Ifges(4,5) - t112 * t89 / 0.2e1 + t116 * t90 / 0.2e1 + t137 * pkin(8)) * t113 + m(5) * (-pkin(3) * t101 + (-t112 * t54 + t116 * t55) * pkin(9)) + (-t2 * t40 + t3 * t39) * mrSges(7,3) + (-t99 / 0.2e1 - t98 / 0.2e1 - t72 / 0.2e1 - t71 / 0.2e1 - t38 / 0.2e1 - t37 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t117 + m(6) * (t18 * t50 + t19 * t51 + t83 * t97) + m(7) * (t12 * t2 + t13 * t3 + t42 * t56); -0.2e1 * pkin(3) * t86 + t112 * t90 + t116 * t89 + 0.2e1 * t56 * t15 + t39 * t16 + t40 * t17 + 0.2e1 * t97 * t44 + t74 * t45 + t75 * t46 + Ifges(4,3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t56 ^ 2) + m(6) * (t50 ^ 2 + t51 ^ 2 + t97 ^ 2) + m(5) * (pkin(9) ^ 2 * t128 + pkin(3) ^ 2) + 0.2e1 * (-t12 * t40 + t13 * t39) * mrSges(7,3) + 0.2e1 * (-t50 * t75 + t51 * t74) * mrSges(6,3) + 0.2e1 * t128 * pkin(9) * mrSges(5,3); t47 * mrSges(5,1) + t20 * mrSges(6,1) - t48 * mrSges(5,2) - t21 * mrSges(6,2) + m(7) * (t5 * t64 + t6 * t65) + (t107 * t21 + t109 * t20) * t147 + t126; -Ifges(5,6) * t130 + m(7) * (t2 * t64 + t3 * t65) + t65 * t22 + t64 * t23 - t19 * mrSges(6,2) + t18 * mrSges(6,1) - t55 * mrSges(5,2) + t54 * mrSges(5,1) + t127 * t117 + (m(6) * (t107 * t19 + t109 * t18) + t107 * t52 + t109 * t53) * pkin(4) + t125 - t144; m(7) * (t12 * t64 + t13 * t65) + t50 * mrSges(6,1) - t51 * mrSges(6,2) + t72 + t71 + t99 + t98 - t124 * pkin(9) + (t39 * t65 - t40 * t64) * mrSges(7,3) + (m(6) * (t107 * t51 + t109 * t50) + (t107 * t74 - t109 * t75) * mrSges(6,3)) * pkin(4) + t122; 0.2e1 * t141 - 0.2e1 * t140 + m(7) * (t64 ^ 2 + t65 ^ 2) - t127 + (0.2e1 * mrSges(6,1) * t109 - 0.2e1 * mrSges(6,2) * t107 + (t107 ^ 2 + t109 ^ 2) * t147) * pkin(4); 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t67; m(6) * t83 + m(7) * t42 + t146; m(6) * t97 + m(7) * t56 + t145; 0; m(6) + m(7); t126; -Ifges(7,3) * t117 + t125; t122; Ifges(7,3) - t140 + t141; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
