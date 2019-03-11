% Calculate joint inertia matrix for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:10
% EndTime: 2019-03-09 09:45:13
% DurationCPUTime: 1.20s
% Computational Cost: add. (1868->261), mult. (3565->376), div. (0->0), fcn. (3860->8), ass. (0->89)
t138 = Ifges(7,4) + Ifges(6,5);
t137 = -Ifges(6,6) + Ifges(7,6);
t136 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t135 = m(6) * pkin(4);
t134 = m(6) + m(7);
t100 = sin(qJ(4));
t102 = cos(qJ(4));
t96 = sin(pkin(10));
t98 = cos(pkin(10));
t71 = t100 * t96 - t102 * t98;
t74 = t100 * t98 + t102 * t96;
t133 = Ifges(5,5) * t100 + Ifges(5,6) * t102 + t137 * t71 + t138 * t74;
t101 = sin(qJ(2));
t127 = -qJ(3) - pkin(7);
t111 = t127 * t101;
t103 = cos(qJ(2));
t79 = t127 * t103;
t97 = sin(pkin(9));
t99 = cos(pkin(9));
t53 = -t111 * t99 - t79 * t97;
t132 = t53 ^ 2;
t131 = 0.2e1 * t53;
t89 = -pkin(2) * t103 - pkin(1);
t130 = 0.2e1 * t89;
t75 = t101 * t99 + t103 * t97;
t117 = t100 * t75;
t73 = t101 * t97 - t103 * t99;
t43 = pkin(3) * t73 - pkin(8) * t75 + t89;
t55 = t111 * t97 - t99 * t79;
t17 = t100 * t43 + t102 * t55;
t13 = -qJ(5) * t117 + t17;
t116 = t102 * t75;
t16 = -t100 * t55 + t102 * t43;
t7 = pkin(4) * t73 - qJ(5) * t116 + t16;
t4 = t13 * t98 + t7 * t96;
t34 = t74 * t75;
t18 = -mrSges(6,2) * t73 - mrSges(6,3) * t34;
t21 = -mrSges(7,2) * t34 + mrSges(7,3) * t73;
t126 = t18 + t21;
t35 = t71 * t75;
t19 = mrSges(6,1) * t73 + mrSges(6,3) * t35;
t20 = -mrSges(7,1) * t73 - mrSges(7,2) * t35;
t125 = -t19 + t20;
t121 = t100 ^ 2 + t102 ^ 2;
t120 = t101 ^ 2 + t103 ^ 2;
t119 = Ifges(5,4) * t100;
t118 = Ifges(5,4) * t102;
t86 = pkin(2) * t97 + pkin(8);
t115 = qJ(5) + t86;
t110 = t115 * t100;
t69 = t115 * t102;
t39 = t110 * t98 + t69 * t96;
t41 = -t110 * t96 + t98 * t69;
t114 = t39 ^ 2 + t41 ^ 2;
t88 = -pkin(2) * t99 - pkin(3);
t15 = t34 * mrSges(6,1) - mrSges(6,2) * t35;
t47 = t71 * mrSges(6,1) + mrSges(6,2) * t74;
t14 = mrSges(7,1) * t34 + t35 * mrSges(7,3);
t46 = mrSges(7,1) * t71 - t74 * mrSges(7,3);
t25 = pkin(4) * t117 + t53;
t3 = -t13 * t96 + t7 * t98;
t78 = -t102 * mrSges(5,1) + t100 * mrSges(5,2);
t109 = mrSges(5,1) * t100 + mrSges(5,2) * t102;
t77 = -pkin(4) * t102 + t88;
t108 = -t46 - t47;
t107 = Ifges(5,5) * t116 + t136 * t73 + t137 * t34 - t138 * t35;
t87 = -pkin(4) * t98 - pkin(5);
t83 = pkin(4) * t96 + qJ(6);
t81 = Ifges(5,1) * t100 + t118;
t80 = Ifges(5,2) * t102 + t119;
t61 = t75 * mrSges(4,2);
t51 = Ifges(6,1) * t74 - Ifges(6,4) * t71;
t50 = Ifges(7,1) * t74 + Ifges(7,5) * t71;
t49 = Ifges(6,4) * t74 - Ifges(6,2) * t71;
t48 = Ifges(7,5) * t74 + Ifges(7,3) * t71;
t45 = mrSges(5,1) * t73 - mrSges(5,3) * t116;
t44 = -mrSges(5,2) * t73 - mrSges(5,3) * t117;
t42 = t109 * t75;
t33 = pkin(5) * t71 - qJ(6) * t74 + t77;
t24 = Ifges(5,5) * t73 + (Ifges(5,1) * t102 - t119) * t75;
t23 = Ifges(5,6) * t73 + (-Ifges(5,2) * t100 + t118) * t75;
t11 = -Ifges(6,1) * t35 - Ifges(6,4) * t34 + Ifges(6,5) * t73;
t10 = -Ifges(7,1) * t35 + Ifges(7,4) * t73 + Ifges(7,5) * t34;
t9 = -Ifges(6,4) * t35 - Ifges(6,2) * t34 + Ifges(6,6) * t73;
t8 = -Ifges(7,5) * t35 + Ifges(7,6) * t73 + Ifges(7,3) * t34;
t5 = pkin(5) * t34 + qJ(6) * t35 + t25;
t2 = -pkin(5) * t73 - t3;
t1 = qJ(6) * t73 + t4;
t6 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t103 + mrSges(3,2) * t101) + t101 * (Ifges(3,1) * t101 + Ifges(3,4) * t103) + t103 * (Ifges(3,4) * t101 + Ifges(3,2) * t103) + t61 * t130 + 0.2e1 * t17 * t44 + 0.2e1 * t16 * t45 + t42 * t131 + 0.2e1 * t25 * t15 + 0.2e1 * t5 * t14 + 0.2e1 * t4 * t18 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t1 * t21 + Ifges(2,3) - (t10 + t11) * t35 + (-t9 + t8) * t34 + 0.2e1 * t120 * pkin(7) * mrSges(3,3) + (mrSges(4,3) * t131 + Ifges(4,1) * t75 - t100 * t23 + t102 * t24) * t75 + (mrSges(4,1) * t130 - 0.2e1 * mrSges(4,3) * t55 + Ifges(4,2) * t73 + (-Ifges(5,6) * t100 - (2 * Ifges(4,4))) * t75 + t107) * t73 + m(3) * (pkin(7) ^ 2 * t120 + pkin(1) ^ 2) + m(4) * (t55 ^ 2 + t89 ^ 2 + t132) + m(5) * (t16 ^ 2 + t17 ^ 2 + t132) + m(6) * (t25 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2); (t78 - mrSges(4,1)) * t53 + m(5) * (t53 * t88 + (-t100 * t16 + t102 * t17) * t86) + (t75 * t81 / 0.2e1 + t86 * t44 + t17 * mrSges(5,3) + t23 / 0.2e1) * t102 + ((-t73 * t97 - t75 * t99) * mrSges(4,3) + m(4) * (-t53 * t99 + t55 * t97)) * pkin(2) + (-t75 * t80 / 0.2e1 - t86 * t45 - t16 * mrSges(5,3) + t24 / 0.2e1) * t100 + m(6) * (t25 * t77 - t3 * t39 + t4 * t41) + m(7) * (t1 * t41 + t2 * t39 + t33 * t5) + (t2 * mrSges(7,2) - t3 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1) * t74 + (-t1 * mrSges(7,2) - t4 * mrSges(6,3) - t9 / 0.2e1 + t8 / 0.2e1) * t71 + Ifges(3,5) * t101 - (t50 / 0.2e1 + t51 / 0.2e1) * t35 + Ifges(3,6) * t103 + t88 * t42 + Ifges(4,5) * t75 + t77 * t15 - Ifges(4,6) * t73 + (t48 / 0.2e1 - t49 / 0.2e1) * t34 + t5 * t46 + t25 * t47 - t55 * mrSges(4,2) + t33 * t14 + t125 * t39 + t126 * t41 + t133 * t73 / 0.2e1 + (-mrSges(3,1) * t101 - mrSges(3,2) * t103) * pkin(7); t100 * t81 + t102 * t80 + 0.2e1 * t33 * t46 + 0.2e1 * t77 * t47 + 0.2e1 * t88 * t78 + Ifges(3,3) + Ifges(4,3) + (t50 + t51) * t74 + (t48 - t49) * t71 + m(7) * (t33 ^ 2 + t114) + m(6) * (t77 ^ 2 + t114) + m(5) * (t121 * t86 ^ 2 + t88 ^ 2) + m(4) * (t97 ^ 2 + t99 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t99 - mrSges(4,2) * t97) * pkin(2) + 0.2e1 * t121 * t86 * mrSges(5,3) + 0.2e1 * (t39 * t74 - t41 * t71) * (mrSges(7,2) + mrSges(6,3)); t73 * mrSges(4,1) + t100 * t44 + t102 * t45 + t61 + t126 * t74 + t125 * t71 + m(6) * (-t3 * t71 + t4 * t74) + m(7) * (t1 * t74 + t2 * t71) + m(5) * (t100 * t17 + t102 * t16) + m(4) * t89; t134 * (t39 * t71 + t41 * t74); m(5) * t121 + m(4) + t134 * (t71 ^ 2 + t74 ^ 2); (m(6) * (t3 * t98 + t4 * t96) + t96 * t18 + t98 * t19) * pkin(4) + m(7) * (t1 * t83 + t2 * t87) - Ifges(5,6) * t117 + t107 + t87 * t20 + t83 * t21 + t16 * mrSges(5,1) - t17 * mrSges(5,2) + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t3 * mrSges(6,1) - t4 * mrSges(6,2); m(7) * (t39 * t87 + t41 * t83) + t41 * mrSges(7,3) - t41 * mrSges(6,2) - t39 * mrSges(6,1) - t39 * mrSges(7,1) - t109 * t86 + (-t71 * t83 + t74 * t87) * mrSges(7,2) + (m(6) * (-t39 * t98 + t41 * t96) + (-t71 * t96 - t74 * t98) * mrSges(6,3)) * pkin(4) + t133; m(7) * (t71 * t87 + t74 * t83) + (-t71 * t98 + t74 * t96) * t135 + t108 - t78; -0.2e1 * t87 * mrSges(7,1) + 0.2e1 * t83 * mrSges(7,3) + m(7) * (t83 ^ 2 + t87 ^ 2) + (0.2e1 * mrSges(6,1) * t98 - 0.2e1 * mrSges(6,2) * t96 + (t96 ^ 2 + t98 ^ 2) * t135) * pkin(4) + t136; m(6) * t25 + m(7) * t5 + t14 + t15; m(6) * t77 + m(7) * t33 - t108; 0; 0; t134; m(7) * t2 + t20; m(7) * t39 + t74 * mrSges(7,2); m(7) * t71; m(7) * t87 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
