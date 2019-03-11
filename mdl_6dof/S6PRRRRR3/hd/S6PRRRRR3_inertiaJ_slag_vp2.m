% Calculate joint inertia matrix for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:48:11
% EndTime: 2019-03-09 00:48:13
% DurationCPUTime: 1.14s
% Computational Cost: add. (2060->322), mult. (4454->471), div. (0->0), fcn. (4808->12), ass. (0->129)
t154 = 2 * pkin(8);
t113 = sin(qJ(3));
t111 = sin(qJ(5));
t112 = sin(qJ(4));
t116 = cos(qJ(5));
t117 = cos(qJ(4));
t77 = t111 * t117 + t112 * t116;
t63 = t77 * t113;
t76 = -t111 * t112 + t116 * t117;
t64 = t76 * t113;
t153 = -Ifges(6,5) * t64 + Ifges(6,6) * t63;
t109 = cos(pkin(6));
t118 = cos(qJ(3));
t108 = sin(pkin(6));
t114 = sin(qJ(2));
t140 = t108 * t114;
t68 = -t109 * t118 + t113 * t140;
t67 = t68 ^ 2;
t152 = -pkin(10) - pkin(9);
t151 = pkin(4) * t111;
t150 = pkin(8) * t118;
t102 = t113 * pkin(8);
t110 = sin(qJ(6));
t115 = cos(qJ(6));
t96 = pkin(4) * t116 + pkin(5);
t66 = t110 * t96 + t115 * t151;
t149 = t66 * mrSges(7,2);
t148 = -Ifges(7,3) - Ifges(6,3);
t30 = -t110 * t64 - t115 * t63;
t31 = -t110 * t63 + t115 * t64;
t147 = -Ifges(7,5) * t31 - Ifges(7,6) * t30;
t137 = t113 * t117;
t84 = -pkin(3) * t118 - pkin(9) * t113 - pkin(2);
t75 = t117 * t84;
t41 = -pkin(10) * t137 + t75 + (-pkin(8) * t112 - pkin(4)) * t118;
t138 = t112 * t113;
t56 = t112 * t84 + t117 * t150;
t49 = -pkin(10) * t138 + t56;
t20 = t111 * t41 + t116 * t49;
t89 = t152 * t112;
t90 = t152 * t117;
t52 = t111 * t89 - t116 * t90;
t85 = -mrSges(5,1) * t117 + mrSges(5,2) * t112;
t146 = t85 - mrSges(4,1);
t145 = Ifges(5,4) * t112;
t144 = Ifges(5,4) * t117;
t143 = t110 * mrSges(7,2);
t142 = t113 * t68;
t70 = t109 * t113 + t118 * t140;
t141 = t118 * t70;
t83 = pkin(4) * t138 + t102;
t119 = cos(qJ(2));
t139 = t108 * t119;
t136 = t112 ^ 2 + t117 ^ 2;
t135 = pkin(5) * t143;
t134 = -Ifges(5,3) + t148;
t47 = -t112 * t70 - t117 * t139;
t48 = -t112 * t139 + t117 * t70;
t22 = -t111 * t48 + t116 * t47;
t23 = t111 * t47 + t116 * t48;
t5 = -t110 * t23 + t115 * t22;
t6 = t110 * t22 + t115 * t23;
t133 = t5 * mrSges(7,1) - t6 * mrSges(7,2);
t97 = -pkin(4) * t117 - pkin(3);
t19 = -t111 * t49 + t116 * t41;
t51 = t111 * t90 + t116 * t89;
t65 = -t110 * t151 + t115 * t96;
t62 = t65 * mrSges(7,1);
t132 = Ifges(7,3) + t62 - t149;
t14 = -pkin(11) * t63 + t20;
t9 = -pkin(5) * t118 - pkin(11) * t64 + t19;
t2 = -t110 * t14 + t115 * t9;
t3 = t110 * t9 + t115 * t14;
t131 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t147;
t130 = mrSges(5,1) * t112 + mrSges(5,2) * t117;
t129 = -t112 * t47 + t117 * t48;
t33 = -pkin(11) * t77 + t51;
t34 = pkin(11) * t76 + t52;
t12 = -t110 * t34 + t115 * t33;
t13 = t110 * t33 + t115 * t34;
t39 = -t110 * t77 + t115 * t76;
t37 = Ifges(7,6) * t39;
t40 = t110 * t76 + t115 * t77;
t38 = Ifges(7,5) * t40;
t128 = t12 * mrSges(7,1) - t13 * mrSges(7,2) + t37 + t38;
t127 = t22 * mrSges(6,1) - t23 * mrSges(6,2) + t133;
t126 = (mrSges(6,1) * t116 - mrSges(6,2) * t111) * pkin(4);
t125 = t19 * mrSges(6,1) - t20 * mrSges(6,2) + t131 - t153;
t72 = Ifges(6,6) * t76;
t73 = Ifges(6,5) * t77;
t124 = t51 * mrSges(6,1) - t52 * mrSges(6,2) + t128 + t72 + t73;
t121 = pkin(8) ^ 2;
t107 = t118 ^ 2;
t105 = t113 ^ 2;
t103 = t108 ^ 2;
t101 = t105 * t121;
t100 = Ifges(5,5) * t112;
t99 = Ifges(5,6) * t117;
t98 = t115 * pkin(5) * mrSges(7,1);
t94 = t103 * t119 ^ 2;
t91 = Ifges(5,5) * t137;
t88 = Ifges(5,1) * t112 + t144;
t87 = Ifges(5,2) * t117 + t145;
t86 = -mrSges(4,1) * t118 + mrSges(4,2) * t113;
t82 = -mrSges(5,1) * t118 - mrSges(5,3) * t137;
t81 = mrSges(5,2) * t118 - mrSges(5,3) * t138;
t71 = t130 * t113;
t61 = -Ifges(5,5) * t118 + (Ifges(5,1) * t117 - t145) * t113;
t60 = -Ifges(5,6) * t118 + (-Ifges(5,2) * t112 + t144) * t113;
t57 = -pkin(5) * t76 + t97;
t55 = -t112 * t150 + t75;
t54 = -mrSges(6,1) * t118 - mrSges(6,3) * t64;
t53 = mrSges(6,2) * t118 - mrSges(6,3) * t63;
t46 = Ifges(6,1) * t77 + Ifges(6,4) * t76;
t45 = Ifges(6,4) * t77 + Ifges(6,2) * t76;
t44 = -mrSges(6,1) * t76 + mrSges(6,2) * t77;
t42 = pkin(5) * t63 + t83;
t32 = mrSges(6,1) * t63 + mrSges(6,2) * t64;
t29 = Ifges(6,1) * t64 - Ifges(6,4) * t63 - Ifges(6,5) * t118;
t28 = Ifges(6,4) * t64 - Ifges(6,2) * t63 - Ifges(6,6) * t118;
t25 = -mrSges(7,1) * t118 - mrSges(7,3) * t31;
t24 = mrSges(7,2) * t118 + mrSges(7,3) * t30;
t17 = Ifges(7,1) * t40 + Ifges(7,4) * t39;
t16 = Ifges(7,4) * t40 + Ifges(7,2) * t39;
t15 = -mrSges(7,1) * t39 + mrSges(7,2) * t40;
t10 = -mrSges(7,1) * t30 + mrSges(7,2) * t31;
t8 = Ifges(7,1) * t31 + Ifges(7,4) * t30 - Ifges(7,5) * t118;
t7 = Ifges(7,4) * t31 + Ifges(7,2) * t30 - Ifges(7,6) * t118;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t67) + m(6) * (t22 ^ 2 + t23 ^ 2 + t67) + m(5) * (t47 ^ 2 + t48 ^ 2 + t67) + m(4) * (t70 ^ 2 + t67 + t94) + m(3) * (t103 * t114 ^ 2 + t109 ^ 2 + t94); mrSges(4,3) * t141 + t22 * t54 + t23 * t53 + t6 * t24 + t5 * t25 + t47 * t82 + t48 * t81 + (-t114 * mrSges(3,2) + (mrSges(3,1) - t86) * t119) * t108 + (t113 * mrSges(4,3) + t10 + t32 + t71) * t68 + m(7) * (t2 * t5 + t3 * t6 + t42 * t68) + m(6) * (t19 * t22 + t20 * t23 + t68 * t83) + m(5) * (pkin(8) * t142 + t47 * t55 + t48 * t56) + m(4) * (pkin(2) * t139 + (t141 + t142) * pkin(8)); -0.2e1 * pkin(2) * t86 + 0.2e1 * t42 * t10 + 0.2e1 * t19 * t54 + 0.2e1 * t2 * t25 + 0.2e1 * t20 * t53 + 0.2e1 * t3 * t24 - t63 * t28 + t64 * t29 + t30 * t7 + t31 * t8 + 0.2e1 * t83 * t32 + 0.2e1 * t55 * t82 + 0.2e1 * t56 * t81 + Ifges(3,3) + (t105 + t107) * mrSges(4,3) * t154 + (Ifges(4,1) * t113 - t112 * t60 + t117 * t61 + t71 * t154) * t113 + m(4) * (pkin(2) ^ 2 + t107 * t121 + t101) + m(5) * (t55 ^ 2 + t56 ^ 2 + t101) + m(6) * (t19 ^ 2 + t20 ^ 2 + t83 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t42 ^ 2) + (-t91 + (Ifges(4,2) - t134) * t118 + (Ifges(5,6) * t112 + (2 * Ifges(4,4))) * t113 + t147 + t153) * t118; -t70 * mrSges(4,2) + (t39 * t6 - t40 * t5) * mrSges(7,3) + (-t22 * t77 + t23 * t76) * mrSges(6,3) + t129 * mrSges(5,3) + (t15 + t44 + t146) * t68 + m(7) * (t12 * t5 + t13 * t6 + t57 * t68) + m(6) * (t22 * t51 + t23 * t52 + t68 * t97) + m(5) * (-pkin(3) * t68 + t129 * pkin(9)); (t60 / 0.2e1 + pkin(9) * t81 + t56 * mrSges(5,3)) * t117 + (-pkin(9) * t82 - t55 * mrSges(5,3) + t61 / 0.2e1) * t112 + m(6) * (t19 * t51 + t20 * t52 + t83 * t97) + m(7) * (t12 * t2 + t13 * t3 + t42 * t57) + (-t100 / 0.2e1 - t99 / 0.2e1 - t73 / 0.2e1 - t72 / 0.2e1 - t38 / 0.2e1 - t37 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t118 + m(5) * (-pkin(3) * t102 + (-t55 * t112 + t56 * t117) * pkin(9)) + (t117 * t88 / 0.2e1 - t112 * t87 / 0.2e1 + Ifges(4,5) + t146 * pkin(8)) * t113 + t83 * t44 + t97 * t32 - pkin(3) * t71 + t76 * t28 / 0.2e1 + t77 * t29 / 0.2e1 - t63 * t45 / 0.2e1 + t64 * t46 / 0.2e1 + t52 * t53 + t51 * t54 + t57 * t10 + (-t2 * t40 + t3 * t39) * mrSges(7,3) + t39 * t7 / 0.2e1 + t40 * t8 / 0.2e1 + t42 * t15 + t30 * t16 / 0.2e1 + t31 * t17 / 0.2e1 + t13 * t24 + t12 * t25 + (-t19 * t77 + t20 * t76) * mrSges(6,3); -0.2e1 * pkin(3) * t85 + t112 * t88 + t117 * t87 + 0.2e1 * t57 * t15 + t39 * t16 + t40 * t17 + 0.2e1 * t97 * t44 + t76 * t45 + t77 * t46 + Ifges(4,3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t57 ^ 2) + m(6) * (t51 ^ 2 + t52 ^ 2 + t97 ^ 2) + m(5) * (t136 * pkin(9) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t12 * t40 + t13 * t39) * mrSges(7,3) + 0.2e1 * (-t51 * t77 + t52 * t76) * mrSges(6,3) + 0.2e1 * t136 * pkin(9) * mrSges(5,3); t47 * mrSges(5,1) - t48 * mrSges(5,2) + m(7) * (t5 * t65 + t6 * t66) + m(6) * (t111 * t23 + t116 * t22) * pkin(4) + t127; m(7) * (t2 * t65 + t3 * t66) + t134 * t118 - Ifges(5,6) * t138 + t91 + t125 + t65 * t25 + t66 * t24 + t55 * mrSges(5,1) - t56 * mrSges(5,2) + (m(6) * (t111 * t20 + t116 * t19) + t111 * t53 + t116 * t54) * pkin(4); m(7) * (t12 * t65 + t13 * t66) + t100 + t99 - t130 * pkin(9) + (t39 * t66 - t40 * t65) * mrSges(7,3) + (m(6) * (t111 * t52 + t116 * t51) + (t111 * t76 - t116 * t77) * mrSges(6,3)) * pkin(4) + t124; -0.2e1 * t149 + 0.2e1 * t62 + 0.2e1 * t126 + m(7) * (t65 ^ 2 + t66 ^ 2) + m(6) * (t111 ^ 2 + t116 ^ 2) * pkin(4) ^ 2 - t134; m(7) * (t110 * t6 + t115 * t5) * pkin(5) + t127; t148 * t118 + (m(7) * (t110 * t3 + t115 * t2) + t110 * t24 + t115 * t25) * pkin(5) + t125; (m(7) * (t110 * t13 + t115 * t12) + (t110 * t39 - t115 * t40) * mrSges(7,3)) * pkin(5) + t124; Ifges(6,3) + t98 + t126 + (m(7) * (t110 * t66 + t115 * t65) - t143) * pkin(5) + t132; -0.2e1 * t135 + 0.2e1 * t98 + m(7) * (t110 ^ 2 + t115 ^ 2) * pkin(5) ^ 2 - t148; t133; -Ifges(7,3) * t118 + t131; t128; t132; Ifges(7,3) + t98 - t135; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
