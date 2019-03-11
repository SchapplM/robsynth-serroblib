% Calculate joint inertia matrix for
% S6PRPRRR1
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPRRR1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:27
% EndTime: 2019-03-08 20:22:28
% DurationCPUTime: 0.76s
% Computational Cost: add. (1030->199), mult. (2257->300), div. (0->0), fcn. (2384->12), ass. (0->85)
t79 = sin(qJ(6));
t83 = cos(qJ(6));
t104 = t79 ^ 2 + t83 ^ 2;
t129 = mrSges(7,3) * t104;
t55 = -mrSges(7,1) * t83 + mrSges(7,2) * t79;
t128 = t55 - mrSges(6,1);
t109 = t79 * mrSges(7,3);
t80 = sin(qJ(5));
t81 = sin(qJ(4));
t84 = cos(qJ(5));
t85 = cos(qJ(4));
t51 = t80 * t81 - t84 * t85;
t53 = t80 * t85 + t81 * t84;
t24 = -mrSges(7,2) * t51 - t109 * t53;
t110 = t53 * t83;
t25 = mrSges(7,1) * t51 - mrSges(7,3) * t110;
t127 = t83 * t24 - t79 * t25;
t75 = sin(pkin(12));
t62 = pkin(2) * t75 + pkin(8);
t117 = pkin(9) + t62;
t101 = t117 * t81;
t43 = t117 * t85;
t19 = t101 * t84 + t43 * t80;
t111 = t53 * t79;
t22 = -mrSges(7,1) * t111 - mrSges(7,2) * t110;
t126 = m(7) * t19 - t22;
t120 = pkin(5) * t51;
t77 = cos(pkin(12));
t63 = -pkin(2) * t77 - pkin(3);
t54 = -pkin(4) * t85 + t63;
t17 = -pkin(10) * t53 + t120 + t54;
t21 = -t101 * t80 + t84 * t43;
t5 = t17 * t83 - t21 * t79;
t6 = t17 * t79 + t21 * t83;
t97 = -t5 * t79 + t6 * t83;
t125 = m(7) * t97 + t127;
t76 = sin(pkin(6));
t82 = sin(qJ(2));
t86 = cos(qJ(2));
t36 = (t75 * t86 + t77 * t82) * t76;
t78 = cos(pkin(6));
t26 = -t36 * t81 + t78 * t85;
t27 = t36 * t85 + t78 * t81;
t9 = -t26 * t84 + t27 * t80;
t124 = t9 ^ 2;
t123 = t19 ^ 2;
t34 = (t75 * t82 - t77 * t86) * t76;
t33 = t34 ^ 2;
t122 = t51 ^ 2;
t74 = t85 ^ 2;
t121 = m(6) * pkin(4);
t119 = t19 * t9;
t118 = t51 * t9;
t115 = mrSges(7,3) * t83;
t114 = Ifges(7,4) * t79;
t113 = Ifges(7,4) * t83;
t112 = t51 * t19;
t106 = Ifges(7,5) * t110 + Ifges(7,3) * t51;
t105 = Ifges(7,5) * t79 + Ifges(7,6) * t83;
t103 = t81 ^ 2 + t74;
t57 = Ifges(7,2) * t83 + t114;
t58 = Ifges(7,1) * t79 + t113;
t102 = t57 * t83 + t58 * t79 + Ifges(6,3);
t100 = t104 * pkin(10);
t65 = pkin(4) * t80 + pkin(10);
t99 = t104 * t65;
t28 = mrSges(6,1) * t51 + t53 * mrSges(6,2);
t11 = t26 * t80 + t27 * t84;
t2 = -t11 * t79 + t34 * t83;
t3 = t11 * t83 + t34 * t79;
t98 = -t2 * t79 + t3 * t83;
t56 = -t85 * mrSges(5,1) + t81 * mrSges(5,2);
t96 = -mrSges(7,1) * t79 - mrSges(7,2) * t83;
t95 = -t26 * t81 + t27 * t85;
t94 = 0.2e1 * t129;
t93 = t129 * t53 + t51 * t55 - t28;
t92 = (mrSges(6,1) * t84 - mrSges(6,2) * t80) * pkin(4);
t91 = -t11 * mrSges(6,2) - t2 * t109 + t3 * t115 + t128 * t9;
t14 = Ifges(7,6) * t51 + (-Ifges(7,2) * t79 + t113) * t53;
t15 = Ifges(7,5) * t51 + (Ifges(7,1) * t83 - t114) * t53;
t90 = -t21 * mrSges(6,2) - t5 * t109 + t6 * t115 + t79 * t15 / 0.2e1 + t83 * t14 / 0.2e1 - t57 * t111 / 0.2e1 + t58 * t110 / 0.2e1 + Ifges(6,5) * t53 + (t105 / 0.2e1 - Ifges(6,6)) * t51 + t128 * t19;
t70 = t78 ^ 2;
t66 = -pkin(4) * t84 - pkin(5);
t48 = t53 ^ 2;
t1 = [m(2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t124) + m(6) * (t11 ^ 2 + t124 + t33) + m(5) * (t26 ^ 2 + t27 ^ 2 + t33) + m(4) * (t36 ^ 2 + t33 + t70) + m(3) * (t70 + (t82 ^ 2 + t86 ^ 2) * t76 ^ 2); -t36 * mrSges(4,2) + t2 * t25 - t9 * t22 + t3 * t24 + (mrSges(3,1) * t86 - mrSges(3,2) * t82) * t76 + (-t11 * t51 + t53 * t9) * mrSges(6,3) + t95 * mrSges(5,3) + (-mrSges(4,1) + t28 + t56) * t34 + m(7) * (t2 * t5 + t3 * t6 + t119) + m(6) * (t11 * t21 + t34 * t54 + t119) + m(5) * (t63 * t34 + t62 * t95) + m(4) * (-t34 * t77 + t36 * t75) * pkin(2); Ifges(5,2) * t74 - 0.2e1 * t19 * t22 + 0.2e1 * t6 * t24 + 0.2e1 * t5 * t25 + 0.2e1 * t54 * t28 + 0.2e1 * t63 * t56 + Ifges(3,3) + Ifges(4,3) + (Ifges(5,1) * t81 + 0.2e1 * Ifges(5,4) * t85) * t81 + (-0.2e1 * t21 * mrSges(6,3) + Ifges(6,2) * t51 + t106) * t51 + (0.2e1 * t19 * mrSges(6,3) + Ifges(6,1) * t53 - t79 * t14 + t83 * t15 + (-Ifges(7,6) * t79 - (2 * Ifges(6,4))) * t51) * t53 + m(7) * (t5 ^ 2 + t6 ^ 2 + t123) + m(6) * (t21 ^ 2 + t54 ^ 2 + t123) + m(5) * (t103 * t62 ^ 2 + t63 ^ 2) + m(4) * (t75 ^ 2 + t77 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (t77 * mrSges(4,1) - t75 * mrSges(4,2)) * pkin(2) + 0.2e1 * t103 * t62 * mrSges(5,3); m(4) * t78 + m(7) * (t53 * t98 + t118) + m(6) * (t11 * t53 + t118) + m(5) * (t26 * t85 + t27 * t81); -t51 * t22 + t127 * t53 + m(7) * (t53 * t97 + t112) + m(6) * (t21 * t53 + t112); m(4) + m(5) * t103 + m(6) * (t48 + t122) + m(7) * (t104 * t48 + t122); t26 * mrSges(5,1) - t27 * mrSges(5,2) + m(7) * (t65 * t98 + t66 * t9) + (t11 * t80 - t84 * t9) * t121 + t91; t90 + (m(6) * (-t19 * t84 + t21 * t80) + (-t51 * t80 - t53 * t84) * mrSges(6,3)) * pkin(4) + (-mrSges(5,1) * t81 - mrSges(5,2) * t85) * t62 + Ifges(5,6) * t85 + Ifges(5,5) * t81 + t126 * t66 + t125 * t65; m(7) * (t66 * t51 + t53 * t99) + (-t51 * t84 + t53 * t80) * t121 + t93 - t56; 0.2e1 * t66 * t55 + Ifges(5,3) + 0.2e1 * t92 + t65 * t94 + m(7) * (t104 * t65 ^ 2 + t66 ^ 2) + m(6) * (t80 ^ 2 + t84 ^ 2) * pkin(4) ^ 2 + t102; m(7) * (-pkin(5) * t9 + pkin(10) * t98) + t91; -pkin(5) * t126 + pkin(10) * t125 + t90; m(7) * (t100 * t53 - t120) + t93; m(7) * (-pkin(5) * t66 + pkin(10) * t99) + (t66 - pkin(5)) * t55 + t92 + (t99 + t100) * mrSges(7,3) + t102; -0.2e1 * pkin(5) * t55 + m(7) * (pkin(10) ^ 2 * t104 + pkin(5) ^ 2) + pkin(10) * t94 + t102; mrSges(7,1) * t2 - mrSges(7,2) * t3; mrSges(7,1) * t5 - mrSges(7,2) * t6 - Ifges(7,6) * t111 + t106; t22; t65 * t96 + t105; pkin(10) * t96 + t105; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
