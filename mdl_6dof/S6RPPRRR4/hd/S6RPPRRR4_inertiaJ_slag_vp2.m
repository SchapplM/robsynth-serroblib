% Calculate joint inertia matrix for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:30
% EndTime: 2019-03-09 02:25:31
% DurationCPUTime: 0.83s
% Computational Cost: add. (1112->242), mult. (1926->345), div. (0->0), fcn. (1793->8), ass. (0->98)
t78 = sin(qJ(6));
t79 = sin(qJ(5));
t81 = cos(qJ(6));
t82 = cos(qJ(5));
t46 = t78 * t82 + t81 * t79;
t80 = sin(qJ(4));
t33 = t46 * t80;
t103 = t81 * t82;
t106 = t79 * t80;
t35 = t80 * t103 - t78 * t106;
t119 = -t33 * mrSges(7,1) - t35 * mrSges(7,2);
t92 = -t79 * mrSges(6,1) - t82 * mrSges(6,2);
t39 = t92 * t80;
t121 = t119 + t39;
t83 = cos(qJ(4));
t120 = Ifges(6,6) * t106 + Ifges(6,3) * t83;
t118 = m(7) * pkin(5);
t117 = t79 / 0.2e1;
t116 = -pkin(9) - pkin(8);
t115 = t80 * pkin(8);
t114 = Ifges(6,4) * t79;
t113 = Ifges(6,4) * t82;
t112 = Ifges(6,5) * t83;
t76 = sin(pkin(10));
t77 = cos(pkin(10));
t84 = -pkin(1) - pkin(2);
t50 = t77 * qJ(2) + t76 * t84;
t44 = -pkin(7) + t50;
t109 = t44 * t76;
t108 = t76 * t80;
t107 = t76 * t83;
t105 = t79 * t83;
t104 = t80 * t82;
t102 = t82 * t83;
t49 = -t76 * qJ(2) + t77 * t84;
t43 = pkin(3) - t49;
t69 = t83 * pkin(4);
t25 = t43 + t69 + t115;
t9 = t44 * t102 + t79 * t25;
t51 = -t82 * mrSges(6,1) + t79 * mrSges(6,2);
t101 = t51 - mrSges(5,1);
t100 = t79 ^ 2 + t82 ^ 2;
t73 = t80 ^ 2;
t75 = t83 ^ 2;
t99 = t73 + t75;
t98 = -Ifges(7,5) * t35 + Ifges(7,6) * t33 + Ifges(7,3) * t83;
t97 = t100 * mrSges(6,3);
t96 = t99 * mrSges(5,3);
t36 = -t76 * t105 - t82 * t77;
t37 = t76 * t102 - t79 * t77;
t12 = t81 * t36 - t78 * t37;
t13 = t78 * t36 + t81 * t37;
t95 = t12 * mrSges(7,1) - t13 * mrSges(7,2);
t23 = t82 * t25;
t8 = -t44 * t105 + t23;
t93 = -t79 * t8 + t82 * t9;
t91 = -t36 * t79 + t37 * t82;
t47 = -t83 * mrSges(6,2) + mrSges(6,3) * t106;
t48 = t83 * mrSges(6,1) + mrSges(6,3) * t104;
t90 = t82 * t47 - t79 * t48;
t55 = t116 * t79;
t56 = t116 * t82;
t18 = t81 * t55 + t78 * t56;
t19 = t78 * t55 - t81 * t56;
t45 = -t78 * t79 + t103;
t40 = Ifges(7,6) * t45;
t41 = Ifges(7,5) * t46;
t89 = t18 * mrSges(7,1) - t19 * mrSges(7,2) + t40 + t41;
t4 = pkin(9) * t104 + t23 + (-t44 * t79 + pkin(5)) * t83;
t5 = pkin(9) * t106 + t9;
t2 = t81 * t4 - t78 * t5;
t3 = t78 * t4 + t81 * t5;
t88 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t98;
t87 = (t81 * mrSges(7,1) - t78 * mrSges(7,2)) * pkin(5);
t71 = t77 ^ 2;
t70 = t76 ^ 2;
t68 = Ifges(6,5) * t79;
t67 = Ifges(6,6) * t82;
t64 = t83 * mrSges(5,1);
t61 = -t82 * pkin(5) - pkin(4);
t59 = t73 * t70;
t54 = Ifges(6,1) * t79 + t113;
t53 = Ifges(6,2) * t82 + t114;
t52 = -t80 * mrSges(5,2) + t64;
t42 = t44 ^ 2;
t38 = t73 * t42;
t31 = t73 * t109;
t30 = t112 + (-Ifges(6,1) * t82 + t114) * t80;
t29 = Ifges(6,6) * t83 + (Ifges(6,2) * t79 - t113) * t80;
t24 = (-pkin(5) * t79 + t44) * t80;
t21 = t83 * mrSges(7,1) + t35 * mrSges(7,3);
t20 = -t83 * mrSges(7,2) + t33 * mrSges(7,3);
t16 = Ifges(7,1) * t46 + Ifges(7,4) * t45;
t15 = Ifges(7,4) * t46 + Ifges(7,2) * t45;
t14 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t7 = -Ifges(7,1) * t35 + Ifges(7,4) * t33 + Ifges(7,5) * t83;
t6 = -Ifges(7,4) * t35 + Ifges(7,2) * t33 + Ifges(7,6) * t83;
t1 = [(2 * pkin(1) * mrSges(3,1)) - 0.2e1 * t49 * mrSges(4,1) + 0.2e1 * t50 * mrSges(4,2) + 0.2e1 * qJ(2) * mrSges(3,3) + 0.2e1 * t24 * t119 + 0.2e1 * t2 * t21 + 0.2e1 * t3 * t20 + t33 * t6 - t35 * t7 + 0.2e1 * t43 * t52 + 0.2e1 * t9 * t47 + 0.2e1 * t8 * t48 + Ifges(3,2) + Ifges(2,3) + Ifges(4,3) - 0.2e1 * t44 * t96 + (Ifges(5,2) * t83 + t120 + t98) * t83 + (Ifges(5,1) * t80 + 0.2e1 * Ifges(5,4) * t83 + t79 * t29 + 0.2e1 * t44 * t39 + (-t30 - t112) * t82) * t80 + m(5) * (t75 * t42 + t43 ^ 2 + t38) + m(4) * (t49 ^ 2 + t50 ^ 2) + m(6) * (t8 ^ 2 + t9 ^ 2 + t38) + m(7) * (t2 ^ 2 + t24 ^ 2 + t3 ^ 2) + m(3) * ((pkin(1) ^ 2) + qJ(2) ^ 2); -m(3) * pkin(1) + t12 * t21 + t13 * t20 + t36 * t48 + t37 * t47 - mrSges(3,1) + (-mrSges(4,1) - t52) * t77 + (t121 * t80 + mrSges(4,2) - t96) * t76 + m(7) * (t24 * t108 + t12 * t2 + t13 * t3) + m(6) * (t36 * t8 + t37 * t9 + t31) + m(5) * (t75 * t109 - t77 * t43 + t31) + m(4) * (t77 * t49 + t76 * t50); m(3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t59) + m(6) * (t36 ^ 2 + t37 ^ 2 + t59) + m(5) * (t75 * t70 + t59 + t71) + m(4) * (t70 + t71); m(7) * (-t33 * t2 - t24 * t83 + t35 * t3) + t35 * t20 - t33 * t21 - t83 * t119 - t83 * t39 + (m(6) * (-t44 * t83 + t93) + t90) * t80; m(6) * (t91 - t107) * t80 + (-t107 * t80 - t33 * t12 + t35 * t13) * m(7); m(4) + m(5) * t99 + m(6) * (t100 * t73 + t75) + m(7) * (t33 ^ 2 + t35 ^ 2 + t75); t82 * t29 / 0.2e1 + t30 * t117 + t46 * t7 / 0.2e1 + t61 * t119 + m(7) * (t18 * t2 + t19 * t3 + t61 * t24) - pkin(4) * t39 + t45 * t6 / 0.2e1 + t33 * t15 / 0.2e1 - t35 * t16 / 0.2e1 + t19 * t20 + t18 * t21 + t24 * t14 + (t68 / 0.2e1 + t67 / 0.2e1 + t41 / 0.2e1 + t40 / 0.2e1 - Ifges(5,6) - t44 * mrSges(5,2)) * t83 + (-t2 * t46 + t3 * t45) * mrSges(7,3) + t93 * mrSges(6,3) + (m(6) * t93 + t90) * pkin(8) + (-Ifges(5,5) - t82 * t54 / 0.2e1 + t53 * t117 + (-m(6) * pkin(4) + t101) * t44) * t80; (-t12 * t46 + t13 * t45) * mrSges(7,3) + t91 * mrSges(6,3) + (-t83 * mrSges(5,2) + (t14 + t101) * t80) * t76 + m(7) * (t61 * t108 + t18 * t12 + t19 * t13) + m(6) * (-pkin(4) * t108 + t91 * pkin(8)); t64 + (-t14 - t51) * t83 + (t33 * t46 + t35 * t45) * mrSges(7,3) + (-mrSges(5,2) + t97) * t80 + m(6) * (t100 * t115 + t69) + m(7) * (-t18 * t33 + t19 * t35 - t61 * t83); -0.2e1 * pkin(4) * t51 + 0.2e1 * t61 * t14 + t45 * t15 + t46 * t16 + t82 * t53 + t79 * t54 + Ifges(5,3) + m(7) * (t18 ^ 2 + t19 ^ 2 + t61 ^ 2) + m(6) * (t100 * pkin(8) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t18 * t46 + t19 * t45) * mrSges(7,3) + 0.2e1 * pkin(8) * t97; -Ifges(6,5) * t104 + t8 * mrSges(6,1) - t9 * mrSges(6,2) + (m(7) * (t2 * t81 + t3 * t78) + t78 * t20 + t81 * t21) * pkin(5) + t88 + t120; t36 * mrSges(6,1) - t37 * mrSges(6,2) + (t12 * t81 + t13 * t78) * t118 + t95; (-t33 * t81 + t35 * t78) * t118 + t121; t67 + t68 + t92 * pkin(8) + (m(7) * (t18 * t81 + t19 * t78) + (t78 * t45 - t81 * t46) * mrSges(7,3)) * pkin(5) + t89; Ifges(6,3) + Ifges(7,3) + m(7) * (t78 ^ 2 + t81 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t87; t88; t95; t119; t89; Ifges(7,3) + t87; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
