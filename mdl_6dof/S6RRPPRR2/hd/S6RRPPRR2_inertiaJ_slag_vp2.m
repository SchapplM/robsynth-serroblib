% Calculate joint inertia matrix for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:52
% EndTime: 2019-03-09 08:49:54
% DurationCPUTime: 0.99s
% Computational Cost: add. (2676->268), mult. (5113->395), div. (0->0), fcn. (5864->10), ass. (0->101)
t103 = sin(pkin(10));
t105 = cos(pkin(10));
t108 = sin(qJ(2));
t111 = cos(qJ(2));
t84 = t103 * t111 + t105 * t108;
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t85 = t102 * t110 + t104 * t107;
t45 = t85 * t84;
t83 = -t102 * t107 + t104 * t110;
t46 = t83 * t84;
t82 = t103 * t108 - t105 * t111;
t136 = Ifges(6,5) * t46 - Ifges(6,6) * t45 + Ifges(6,3) * t82;
t128 = -qJ(3) - pkin(7);
t120 = t128 * t108;
t91 = t128 * t111;
t64 = -t103 * t91 - t105 * t120;
t135 = t64 ^ 2;
t134 = 0.2e1 * t64;
t96 = -pkin(2) * t111 - pkin(1);
t133 = 0.2e1 * t96;
t132 = t104 / 0.2e1;
t130 = pkin(2) * t103;
t93 = qJ(4) + t130;
t131 = pkin(8) + t93;
t129 = pkin(2) * t105;
t123 = t104 * t84;
t52 = pkin(3) * t82 - qJ(4) * t84 + t96;
t66 = t103 * t120 - t105 * t91;
t33 = -t102 * t66 + t104 * t52;
t18 = pkin(4) * t82 - pkin(8) * t123 + t33;
t124 = t102 * t84;
t34 = t102 * t52 + t104 * t66;
t25 = -pkin(8) * t124 + t34;
t7 = t107 * t18 + t110 * t25;
t74 = t131 * t102;
t75 = t131 * t104;
t50 = -t107 * t74 + t110 * t75;
t51 = mrSges(5,1) * t124 + mrSges(5,2) * t123;
t127 = t102 ^ 2 + t104 ^ 2;
t126 = Ifges(5,4) * t102;
t125 = Ifges(5,4) * t104;
t122 = t108 ^ 2 + t111 ^ 2;
t106 = sin(qJ(6));
t109 = cos(qJ(6));
t26 = -t106 * t46 - t109 * t45;
t27 = -t106 * t45 + t109 * t46;
t121 = Ifges(7,5) * t27 + Ifges(7,6) * t26 + Ifges(7,3) * t82;
t95 = -pkin(3) - t129;
t29 = t45 * mrSges(6,1) + t46 * mrSges(6,2);
t60 = -t83 * mrSges(6,1) + t85 * mrSges(6,2);
t10 = -t26 * mrSges(7,1) + t27 * mrSges(7,2);
t56 = -t106 * t85 + t109 * t83;
t57 = t106 * t83 + t109 * t85;
t30 = -t56 * mrSges(7,1) + t57 * mrSges(7,2);
t88 = -t104 * mrSges(5,1) + t102 * mrSges(5,2);
t6 = -t107 * t25 + t110 * t18;
t49 = -t107 * t75 - t110 * t74;
t41 = pkin(4) * t124 + t64;
t119 = -t33 * t102 + t34 * t104;
t37 = -pkin(9) * t85 + t49;
t38 = pkin(9) * t83 + t50;
t12 = -t106 * t38 + t109 * t37;
t13 = t106 * t37 + t109 * t38;
t54 = Ifges(7,6) * t56;
t55 = Ifges(7,5) * t57;
t118 = t12 * mrSges(7,1) - t13 * mrSges(7,2) + t54 + t55;
t87 = -pkin(4) * t104 + t95;
t4 = pkin(5) * t82 - pkin(9) * t46 + t6;
t5 = -pkin(9) * t45 + t7;
t2 = -t106 * t5 + t109 * t4;
t3 = t106 * t4 + t109 * t5;
t117 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t121;
t116 = (mrSges(7,1) * t109 - mrSges(7,2) * t106) * pkin(5);
t115 = -t30 - t60;
t90 = Ifges(5,1) * t102 + t125;
t89 = Ifges(5,2) * t104 + t126;
t81 = Ifges(6,5) * t85;
t80 = Ifges(6,6) * t83;
t76 = t84 * mrSges(4,2);
t67 = -pkin(5) * t83 + t87;
t62 = Ifges(6,1) * t85 + Ifges(6,4) * t83;
t61 = Ifges(6,4) * t85 + Ifges(6,2) * t83;
t59 = mrSges(5,1) * t82 - mrSges(5,3) * t123;
t58 = -mrSges(5,2) * t82 - mrSges(5,3) * t124;
t40 = t82 * Ifges(5,5) + (Ifges(5,1) * t104 - t126) * t84;
t39 = t82 * Ifges(5,6) + (-Ifges(5,2) * t102 + t125) * t84;
t36 = mrSges(6,1) * t82 - mrSges(6,3) * t46;
t35 = -mrSges(6,2) * t82 - mrSges(6,3) * t45;
t32 = Ifges(7,1) * t57 + Ifges(7,4) * t56;
t31 = Ifges(7,4) * t57 + Ifges(7,2) * t56;
t28 = pkin(5) * t45 + t41;
t20 = Ifges(6,1) * t46 - Ifges(6,4) * t45 + Ifges(6,5) * t82;
t19 = Ifges(6,4) * t46 - Ifges(6,2) * t45 + Ifges(6,6) * t82;
t17 = mrSges(7,1) * t82 - mrSges(7,3) * t27;
t16 = -mrSges(7,2) * t82 + mrSges(7,3) * t26;
t9 = Ifges(7,1) * t27 + Ifges(7,4) * t26 + Ifges(7,5) * t82;
t8 = Ifges(7,4) * t27 + Ifges(7,2) * t26 + Ifges(7,6) * t82;
t1 = [-0.2e1 * pkin(1) * (-t111 * mrSges(3,1) + t108 * mrSges(3,2)) + t108 * (Ifges(3,1) * t108 + Ifges(3,4) * t111) + t111 * (Ifges(3,4) * t108 + Ifges(3,2) * t111) + 0.2e1 * t34 * t58 + 0.2e1 * t33 * t59 - t45 * t19 + t46 * t20 + 0.2e1 * t7 * t35 + 0.2e1 * t6 * t36 + 0.2e1 * t41 * t29 + 0.2e1 * t122 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t133 - 0.2e1 * t66 * mrSges(4,3) + (Ifges(5,3) + Ifges(4,2)) * t82 + (Ifges(5,5) * t104 - Ifges(5,6) * t102 - (2 * Ifges(4,4))) * t84 + t121 + t136) * t82 + t26 * t8 + t27 * t9 + 0.2e1 * t28 * t10 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + t76 * t133 + t51 * t134 + m(3) * (t122 * pkin(7) ^ 2 + pkin(1) ^ 2) + Ifges(2,3) + (mrSges(4,3) * t134 + Ifges(4,1) * t84 - t102 * t39 + t104 * t40) * t84 + m(4) * (t66 ^ 2 + t96 ^ 2 + t135) + m(5) * (t33 ^ 2 + t34 ^ 2 + t135) + m(7) * (t2 ^ 2 + t28 ^ 2 + t3 ^ 2) + m(6) * (t41 ^ 2 + t6 ^ 2 + t7 ^ 2); (-t108 * mrSges(3,1) - t111 * mrSges(3,2)) * pkin(7) + (t93 * t58 + t39 / 0.2e1) * t104 + (-t93 * t59 + t40 / 0.2e1) * t102 + (-t2 * t57 + t3 * t56) * mrSges(7,3) + (-t6 * t85 + t7 * t83) * mrSges(6,3) + m(4) * (t103 * t66 - t105 * t64) * pkin(2) + Ifges(3,5) * t108 + Ifges(3,6) * t111 + t85 * t20 / 0.2e1 + t87 * t29 + t95 * t51 + t83 * t19 / 0.2e1 + t56 * t8 / 0.2e1 + t57 * t9 / 0.2e1 + t41 * t60 - t45 * t61 / 0.2e1 + t46 * t62 / 0.2e1 - t66 * mrSges(4,2) + t67 * t10 + t49 * t36 + t50 * t35 + m(5) * (t119 * t93 + t64 * t95) + t119 * mrSges(5,3) + t28 * t30 + t26 * t31 / 0.2e1 + t27 * t32 / 0.2e1 + t13 * t16 + t12 * t17 + (t88 - mrSges(4,1)) * t64 + m(6) * (t41 * t87 + t49 * t6 + t50 * t7) + m(7) * (t12 * t2 + t13 * t3 + t28 * t67) + (-t102 * t89 / 0.2e1 + t90 * t132 + Ifges(4,5) - mrSges(4,3) * t129) * t84 + (Ifges(5,5) * t102 / 0.2e1 + Ifges(5,6) * t132 + t81 / 0.2e1 + t80 / 0.2e1 + t55 / 0.2e1 + t54 / 0.2e1 - Ifges(4,6) - mrSges(4,3) * t130) * t82; t102 * t90 + t104 * t89 + 0.2e1 * t67 * t30 + t56 * t31 + t57 * t32 + 0.2e1 * t87 * t60 + t83 * t61 + t85 * t62 + 0.2e1 * t95 * t88 + Ifges(3,3) + Ifges(4,3) + m(7) * (t12 ^ 2 + t13 ^ 2 + t67 ^ 2) + m(6) * (t49 ^ 2 + t50 ^ 2 + t87 ^ 2) + m(5) * (t127 * t93 ^ 2 + t95 ^ 2) + m(4) * (t103 ^ 2 + t105 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t105 - mrSges(4,2) * t103) * pkin(2) + 0.2e1 * (-t12 * t57 + t13 * t56) * mrSges(7,3) + 0.2e1 * (-t49 * t85 + t50 * t83) * mrSges(6,3) + 0.2e1 * t127 * t93 * mrSges(5,3); t82 * mrSges(4,1) + t102 * t58 + t104 * t59 + t57 * t16 + t56 * t17 + t85 * t35 + t83 * t36 + t76 + m(7) * (t2 * t56 + t3 * t57) + m(6) * (t6 * t83 + t7 * t85) + m(5) * (t102 * t34 + t104 * t33) + m(4) * t96; m(7) * (t12 * t56 + t13 * t57) + m(6) * (t49 * t83 + t50 * t85); m(4) + m(5) * t127 + m(6) * (t83 ^ 2 + t85 ^ 2) + m(7) * (t56 ^ 2 + t57 ^ 2); m(5) * t64 + m(6) * t41 + m(7) * t28 + t10 + t29 + t51; m(5) * t95 + m(6) * t87 + m(7) * t67 - t115 + t88; 0; m(5) + m(6) + m(7); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (t109 * t17 + m(7) * (t106 * t3 + t109 * t2) + t106 * t16) * pkin(5) + t117 + t136; t49 * mrSges(6,1) - t50 * mrSges(6,2) + t80 + t81 + (m(7) * (t106 * t13 + t109 * t12) + (t106 * t56 - t109 * t57) * mrSges(7,3)) * pkin(5) + t118; m(7) * (t106 * t57 + t109 * t56) * pkin(5) + t115; 0; Ifges(6,3) + Ifges(7,3) + m(7) * (t106 ^ 2 + t109 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t116; t117; t118; -t30; 0; Ifges(7,3) + t116; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
