% Calculate joint inertia matrix for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:15
% EndTime: 2019-03-08 20:18:17
% DurationCPUTime: 0.68s
% Computational Cost: add. (557->209), mult. (1221->279), div. (0->0), fcn. (999->8), ass. (0->93)
t122 = Ifges(7,2) + Ifges(6,3);
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t91 = t65 ^ 2 + t68 ^ 2;
t63 = sin(pkin(6));
t67 = sin(qJ(2));
t107 = t63 * t67;
t70 = cos(qJ(2));
t106 = t63 * t70;
t64 = cos(pkin(6));
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t21 = -t66 * t106 + t64 * t69;
t6 = t65 * t107 + t21 * t68;
t114 = t6 * t68;
t4 = -t68 * t107 + t21 * t65;
t121 = t4 * t65 + t114;
t120 = m(7) + m(6);
t119 = mrSges(7,2) + mrSges(6,3);
t118 = -m(7) * pkin(5) - mrSges(7,1);
t117 = t119 * t91;
t19 = t69 * t106 + t64 * t66;
t17 = t19 ^ 2;
t113 = Ifges(6,4) * t65;
t112 = Ifges(6,4) * t68;
t111 = Ifges(7,5) * t65;
t110 = Ifges(7,5) * t68;
t109 = Ifges(6,6) * t66;
t108 = Ifges(7,6) * t66;
t105 = t65 * t66;
t104 = t65 * t69;
t103 = t66 * t21;
t71 = -pkin(2) - pkin(8);
t102 = t66 * t71;
t29 = t66 * pkin(4) - t69 * pkin(9) + qJ(3);
t101 = t68 * t29;
t100 = t68 * t69;
t99 = t69 * t19;
t79 = t65 * mrSges(7,1) - t68 * mrSges(7,3);
t22 = t79 * t69;
t80 = t65 * mrSges(6,1) + t68 * mrSges(6,2);
t23 = t80 * t69;
t98 = -t22 - t23;
t25 = -t66 * mrSges(6,2) - mrSges(6,3) * t104;
t28 = -mrSges(7,2) * t104 + t66 * mrSges(7,3);
t97 = t25 + t28;
t26 = t66 * mrSges(6,1) - mrSges(6,3) * t100;
t27 = -t66 * mrSges(7,1) + mrSges(7,2) * t100;
t96 = -t26 + t27;
t32 = -t68 * mrSges(6,1) + t65 * mrSges(6,2);
t95 = t32 - mrSges(5,1);
t94 = t66 * mrSges(5,1) + t69 * mrSges(5,2) + mrSges(4,3);
t10 = t68 * t102 + t65 * t29;
t93 = t91 * pkin(9) * t66;
t92 = t91 * pkin(9) ^ 2;
t59 = t66 ^ 2;
t61 = t69 ^ 2;
t90 = -t61 - t59;
t31 = -t68 * mrSges(7,1) - t65 * mrSges(7,3);
t89 = -t31 - t95;
t87 = t90 * mrSges(5,3);
t86 = t121 * pkin(9);
t84 = Ifges(7,6) * t104 + t122 * t66 + (Ifges(7,4) + Ifges(6,5)) * t100;
t7 = t66 * qJ(6) + t10;
t8 = -t101 + (t65 * t71 - pkin(5)) * t66;
t82 = t65 * t8 + t68 * t7;
t9 = -t65 * t102 + t101;
t81 = t10 * t68 - t65 * t9;
t78 = -pkin(5) * t65 + qJ(6) * t68;
t77 = -t99 + t103;
t75 = t96 * t65 + t97 * t68;
t74 = m(7) * t78 - t79 - t80;
t72 = qJ(3) ^ 2;
t62 = t71 ^ 2;
t57 = t63 ^ 2;
t54 = Ifges(7,4) * t65;
t53 = Ifges(6,5) * t65;
t51 = Ifges(6,6) * t68;
t49 = t61 * t71;
t48 = t61 * t62;
t47 = t57 * t67 ^ 2;
t40 = qJ(3) * t107;
t37 = Ifges(6,1) * t65 + t112;
t36 = Ifges(7,1) * t65 - t110;
t35 = Ifges(6,2) * t68 + t113;
t34 = -Ifges(7,3) * t68 + t111;
t30 = -t68 * pkin(5) - t65 * qJ(6) - pkin(4);
t15 = Ifges(6,5) * t66 + (Ifges(6,1) * t68 - t113) * t69;
t14 = Ifges(7,4) * t66 + (Ifges(7,1) * t68 + t111) * t69;
t13 = t109 + (-Ifges(6,2) * t65 + t112) * t69;
t12 = t108 + (Ifges(7,3) * t65 + t110) * t69;
t11 = (-t71 - t78) * t69;
t1 = [m(2) + m(5) * (t21 ^ 2 + t17 + t47) + 0.2e1 * (m(4) / 0.2e1 + m(3) / 0.2e1) * (t57 * t70 ^ 2 + t64 ^ 2 + t47) + t120 * (t4 ^ 2 + t6 ^ 2 + t17); -mrSges(5,3) * t103 + t97 * t6 + t96 * t4 + (t69 * mrSges(5,3) - t98) * t19 + ((mrSges(3,1) - mrSges(4,2)) * t70 + (-mrSges(3,2) + t94) * t67) * t63 + m(6) * (t10 * t6 - t9 * t4 - t71 * t99) + m(7) * (t11 * t19 + t8 * t4 + t7 * t6) + m(5) * (t77 * t71 + t40) + m(4) * (pkin(2) * t106 + t40); -0.2e1 * pkin(2) * mrSges(4,2) + 0.2e1 * t10 * t25 + 0.2e1 * t11 * t22 + 0.2e1 * t9 * t26 + 0.2e1 * t8 * t27 + 0.2e1 * t7 * t28 + Ifges(4,1) + Ifges(3,3) + (Ifges(5,2) * t66 + t84) * t66 + m(5) * (t59 * t62 + t48 + t72) + m(4) * (pkin(2) ^ 2 + t72) + m(6) * (t10 ^ 2 + t9 ^ 2 + t48) + m(7) * (t11 ^ 2 + t7 ^ 2 + t8 ^ 2) + (Ifges(5,1) * t69 - 0.2e1 * Ifges(5,4) * t66 - 0.2e1 * t71 * t23 + (t14 + t15) * t68 + (t12 - t13 - t109) * t65) * t69 + 0.2e1 * t94 * qJ(3) + 0.2e1 * t71 * t87; -m(4) * t106 + m(5) * t77 + t120 * (t4 * t105 + t66 * t114 - t99); -m(4) * pkin(2) + mrSges(4,2) + t98 * t69 + t87 + t75 * t66 + m(7) * (-t69 * t11 + t82 * t66) + m(6) * (t81 * t66 + t49) + m(5) * (t59 * t71 + t49); m(4) - m(5) * t90 + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t91 * t59 + t61); -t21 * mrSges(5,2) - t89 * t19 + m(6) * (-pkin(4) * t19 + t86) + m(7) * (t30 * t19 + t86) + t119 * t121; -pkin(4) * t23 + t30 * t22 + (m(7) * t30 + t31) * t11 + (-t71 * mrSges(5,2) + t53 / 0.2e1 + t51 / 0.2e1 - Ifges(5,6) + t54 / 0.2e1) * t66 + (t10 * mrSges(6,3) + t7 * mrSges(7,2) - t108 / 0.2e1 - t12 / 0.2e1 + t13 / 0.2e1) * t68 + (t8 * mrSges(7,2) - t9 * mrSges(6,3) + t14 / 0.2e1 + t15 / 0.2e1) * t65 + (m(6) * t81 + m(7) * t82 + t75) * pkin(9) + (Ifges(5,5) + (t36 / 0.2e1 + t37 / 0.2e1) * t68 + (t34 / 0.2e1 - t35 / 0.2e1) * t65 + (m(6) * pkin(4) - t95) * t71) * t69; t89 * t69 + m(7) * (-t30 * t69 + t93) + m(6) * (pkin(4) * t69 + t93) + (-mrSges(5,2) + t117) * t66; -0.2e1 * pkin(4) * t32 + 0.2e1 * t30 * t31 + Ifges(5,3) + (t35 - t34) * t68 + (t37 + t36) * t65 + m(7) * (t30 ^ 2 + t92) + m(6) * (pkin(4) ^ 2 + t92) + 0.2e1 * pkin(9) * t117; (m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3)) * t6 + (-mrSges(6,1) + t118) * t4; -Ifges(6,6) * t104 - pkin(5) * t27 + m(7) * (-pkin(5) * t8 + qJ(6) * t7) + qJ(6) * t28 + t7 * mrSges(7,3) - t8 * mrSges(7,1) - t10 * mrSges(6,2) + t9 * mrSges(6,1) + t84; t74 * t66; t78 * mrSges(7,2) - Ifges(7,6) * t68 + t74 * pkin(9) + t51 + t53 + t54; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t122; m(7) * t4; m(7) * t8 + t27; m(7) * t105; (m(7) * pkin(9) + mrSges(7,2)) * t65; t118; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
