% Calculate joint inertia matrix for
% S6PRPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:26:58
% EndTime: 2019-03-08 20:27:00
% DurationCPUTime: 0.88s
% Computational Cost: add. (1059->244), mult. (2407->356), div. (0->0), fcn. (2485->12), ass. (0->100)
t76 = sin(qJ(6));
t77 = sin(qJ(5));
t80 = cos(qJ(6));
t81 = cos(qJ(5));
t46 = t76 * t81 + t77 * t80;
t78 = sin(qJ(4));
t36 = t46 * t78;
t45 = -t76 * t77 + t80 * t81;
t37 = t45 * t78;
t12 = t36 * mrSges(7,1) + t37 * mrSges(7,2);
t92 = mrSges(6,1) * t77 + mrSges(6,2) * t81;
t41 = t92 * t78;
t116 = -t12 - t41;
t72 = sin(pkin(12));
t61 = pkin(2) * t72 + pkin(8);
t117 = 0.2e1 * t61;
t73 = sin(pkin(6));
t74 = cos(pkin(12));
t79 = sin(qJ(2));
t83 = cos(qJ(2));
t33 = (t72 * t83 + t74 * t79) * t73;
t75 = cos(pkin(6));
t82 = cos(qJ(4));
t17 = t33 * t78 - t75 * t82;
t16 = t17 ^ 2;
t31 = (t72 * t79 - t74 * t83) * t73;
t115 = t31 ^ 2;
t114 = m(7) * pkin(5);
t113 = t81 / 0.2e1;
t112 = -pkin(10) - pkin(9);
t111 = pkin(9) * t78;
t110 = t82 * pkin(4);
t109 = Ifges(6,4) * t77;
t108 = Ifges(6,4) * t81;
t107 = t17 * t78;
t19 = t33 * t82 + t75 * t78;
t106 = t19 * t82;
t105 = t61 * t82;
t104 = t77 * t78;
t103 = t78 * t81;
t102 = t82 * t17;
t101 = -Ifges(7,3) - Ifges(6,3);
t100 = -Ifges(7,5) * t37 + Ifges(7,6) * t36;
t62 = -pkin(2) * t74 - pkin(3);
t42 = -t110 + t62 - t111;
t15 = t81 * t105 + t77 * t42;
t50 = -mrSges(6,1) * t81 + mrSges(6,2) * t77;
t99 = t50 - mrSges(5,1);
t98 = t77 ^ 2 + t81 ^ 2;
t69 = t78 ^ 2;
t71 = t82 ^ 2;
t97 = t69 + t71;
t20 = -mrSges(7,1) * t45 + mrSges(7,2) * t46;
t96 = -t20 - t99;
t7 = -t19 * t77 + t31 * t81;
t8 = t19 * t81 + t31 * t77;
t2 = t7 * t80 - t76 * t8;
t3 = t7 * t76 + t8 * t80;
t95 = t2 * mrSges(7,1) - t3 * mrSges(7,2);
t94 = t98 * mrSges(6,3);
t93 = -t7 * t77 + t8 * t81;
t39 = t81 * t42;
t14 = -t77 * t105 + t39;
t91 = -t14 * t77 + t15 * t81;
t47 = mrSges(6,2) * t82 - mrSges(6,3) * t104;
t48 = -mrSges(6,1) * t82 - mrSges(6,3) * t103;
t90 = t81 * t47 - t77 * t48;
t11 = -pkin(10) * t103 + t39 + (-t61 * t77 - pkin(5)) * t82;
t13 = -pkin(10) * t104 + t15;
t5 = t11 * t80 - t13 * t76;
t6 = t11 * t76 + t13 * t80;
t89 = t5 * mrSges(7,1) - t6 * mrSges(7,2) - t100;
t56 = t112 * t77;
t57 = t112 * t81;
t24 = t56 * t80 + t57 * t76;
t25 = t56 * t76 - t57 * t80;
t43 = Ifges(7,6) * t45;
t44 = Ifges(7,5) * t46;
t88 = t24 * mrSges(7,1) - t25 * mrSges(7,2) + t43 + t44;
t87 = (mrSges(7,1) * t80 - mrSges(7,2) * t76) * pkin(5);
t67 = t75 ^ 2;
t65 = Ifges(6,5) * t77;
t64 = Ifges(6,6) * t81;
t63 = -pkin(5) * t81 - pkin(4);
t59 = t61 ^ 2;
t58 = Ifges(6,5) * t103;
t54 = t69 * t59;
t53 = Ifges(6,1) * t77 + t108;
t52 = Ifges(6,2) * t81 + t109;
t51 = -t82 * mrSges(5,1) + t78 * mrSges(5,2);
t40 = (pkin(5) * t77 + t61) * t78;
t35 = -Ifges(6,5) * t82 + (Ifges(6,1) * t81 - t109) * t78;
t34 = -Ifges(6,6) * t82 + (-Ifges(6,2) * t77 + t108) * t78;
t27 = -mrSges(7,1) * t82 - mrSges(7,3) * t37;
t26 = mrSges(7,2) * t82 - mrSges(7,3) * t36;
t22 = Ifges(7,1) * t46 + Ifges(7,4) * t45;
t21 = Ifges(7,4) * t46 + Ifges(7,2) * t45;
t10 = Ifges(7,1) * t37 - Ifges(7,4) * t36 - Ifges(7,5) * t82;
t9 = Ifges(7,4) * t37 - Ifges(7,2) * t36 - Ifges(7,6) * t82;
t1 = [m(2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t16) + m(6) * (t7 ^ 2 + t8 ^ 2 + t16) + m(5) * (t19 ^ 2 + t115 + t16) + m(4) * (t33 ^ 2 + t115 + t67) + m(3) * (t67 + (t79 ^ 2 + t83 ^ 2) * t73 ^ 2); mrSges(5,3) * t106 - t33 * mrSges(4,2) + t2 * t27 + t3 * t26 + t8 * t47 + t7 * t48 + (mrSges(3,1) * t83 - mrSges(3,2) * t79) * t73 + (-mrSges(4,1) + t51) * t31 + (t78 * mrSges(5,3) - t116) * t17 + m(7) * (t17 * t40 + t2 * t5 + t3 * t6) + m(6) * (t61 * t107 + t14 * t7 + t15 * t8) + m(5) * (t62 * t31 + (t106 + t107) * t61) + m(4) * (-t31 * t74 + t33 * t72) * pkin(2); t37 * t10 + 0.2e1 * t40 * t12 + 0.2e1 * t14 * t48 + 0.2e1 * t15 * t47 + 0.2e1 * t6 * t26 + 0.2e1 * t5 * t27 - t36 * t9 + 0.2e1 * t62 * t51 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t78 + t41 * t117 - t34 * t77 + t35 * t81) * t78 + m(5) * (t59 * t71 + t62 ^ 2 + t54) + m(6) * (t14 ^ 2 + t15 ^ 2 + t54) + m(7) * (t40 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(4) * (t72 ^ 2 + t74 ^ 2) * pkin(2) ^ 2 + (-t58 + (Ifges(5,2) - t101) * t82 + (Ifges(6,6) * t77 + (2 * Ifges(5,4))) * t78 + t100) * t82 + 0.2e1 * (mrSges(4,1) * t74 - mrSges(4,2) * t72) * pkin(2) + t97 * mrSges(5,3) * t117; m(4) * t75 + m(7) * (-t2 * t36 + t3 * t37 - t102) + m(6) * (t93 * t78 - t102) + m(5) * (t19 * t78 - t102); m(7) * (-t36 * t5 + t37 * t6) + t37 * t26 - t36 * t27 + (m(6) * (t91 - t105) + t90) * t78 + (-m(7) * t40 + t116) * t82; m(4) + m(5) * t97 + m(6) * (t98 * t69 + t71) + m(7) * (t36 ^ 2 + t37 ^ 2 + t71); -t19 * mrSges(5,2) + (-t2 * t46 + t3 * t45) * mrSges(7,3) + t93 * mrSges(6,3) - t96 * t17 + m(7) * (t17 * t63 + t2 * t24 + t25 * t3) + m(6) * (-pkin(4) * t17 + t93 * pkin(9)); t34 * t113 + t77 * t35 / 0.2e1 + t63 * t12 + t45 * t9 / 0.2e1 + t46 * t10 / 0.2e1 + m(7) * (t24 * t5 + t25 * t6 + t40 * t63) - t36 * t21 / 0.2e1 + t37 * t22 / 0.2e1 + t40 * t20 - pkin(4) * t41 + t25 * t26 + t24 * t27 + (-t65 / 0.2e1 - t64 / 0.2e1 - t44 / 0.2e1 - t43 / 0.2e1 + Ifges(5,6) - t61 * mrSges(5,2)) * t82 + (t45 * t6 - t46 * t5) * mrSges(7,3) + t91 * mrSges(6,3) + (m(6) * t91 + t90) * pkin(9) + (Ifges(5,5) + t53 * t113 - t77 * t52 / 0.2e1 + (-m(6) * pkin(4) + t99) * t61) * t78; (t36 * t46 + t37 * t45) * mrSges(7,3) + t96 * t82 + (-mrSges(5,2) + t94) * t78 + m(6) * (t98 * t111 + t110) + m(7) * (-t24 * t36 + t25 * t37 - t63 * t82); -0.2e1 * pkin(4) * t50 + 0.2e1 * t63 * t20 + t45 * t21 + t46 * t22 + t81 * t52 + t77 * t53 + Ifges(5,3) + m(7) * (t24 ^ 2 + t25 ^ 2 + t63 ^ 2) + m(6) * (t98 * pkin(9) ^ 2 + pkin(4) ^ 2) + 0.2e1 * (-t24 * t46 + t25 * t45) * mrSges(7,3) + 0.2e1 * pkin(9) * t94; t7 * mrSges(6,1) - t8 * mrSges(6,2) + (t2 * t80 + t3 * t76) * t114 + t95; -Ifges(6,6) * t104 + t14 * mrSges(6,1) - t15 * mrSges(6,2) + t58 + t101 * t82 + (m(7) * (t5 * t80 + t6 * t76) + t76 * t26 + t80 * t27) * pkin(5) + t89; (-t36 * t80 + t37 * t76) * t114 + t116; t64 + t65 - t92 * pkin(9) + (m(7) * (t24 * t80 + t25 * t76) + (t45 * t76 - t46 * t80) * mrSges(7,3)) * pkin(5) + t88; m(7) * (t76 ^ 2 + t80 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t87 - t101; t95; -Ifges(7,3) * t82 + t89; -t12; t88; Ifges(7,3) + t87; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
