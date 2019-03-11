% Calculate joint inertia matrix for
% S6PRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
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
% Datum: 2019-03-08 22:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRPP1_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP1_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP1_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP1_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:42:54
% EndTime: 2019-03-08 22:42:56
% DurationCPUTime: 0.96s
% Computational Cost: add. (1167->281), mult. (2582->393), div. (0->0), fcn. (2585->10), ass. (0->106)
t138 = 2 * pkin(8);
t137 = m(6) * pkin(4);
t136 = m(6) + m(7);
t135 = mrSges(7,2) + mrSges(6,3);
t101 = sin(qJ(3));
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t96 = sin(pkin(11));
t98 = cos(pkin(11));
t64 = t100 * t98 + t103 * t96;
t50 = t64 * t101;
t63 = t100 * t96 - t98 * t103;
t51 = t63 * t101;
t15 = t50 * mrSges(7,1) + t51 * mrSges(7,3);
t16 = t50 * mrSges(6,1) - t51 * mrSges(6,2);
t134 = t15 + t16;
t22 = t63 * mrSges(7,1) - t64 * mrSges(7,3);
t23 = t63 * mrSges(6,1) + t64 * mrSges(6,2);
t133 = t22 + t23;
t104 = cos(qJ(3));
t102 = sin(qJ(2));
t97 = sin(pkin(6));
t122 = t102 * t97;
t99 = cos(pkin(6));
t53 = t101 * t122 - t99 * t104;
t52 = t53 ^ 2;
t132 = pkin(8) * t104;
t90 = t101 * pkin(8);
t130 = -qJ(5) - pkin(9);
t118 = t101 * t103;
t71 = -pkin(3) * t104 - pkin(9) * t101 - pkin(2);
t66 = t103 * t71;
t20 = -qJ(5) * t118 + t66 + (-pkin(8) * t100 - pkin(4)) * t104;
t119 = t100 * t101;
t40 = t100 * t71 + t103 * t132;
t30 = -qJ(5) * t119 + t40;
t5 = t96 * t20 + t98 * t30;
t35 = -mrSges(7,2) * t50 - mrSges(7,3) * t104;
t36 = mrSges(6,2) * t104 - mrSges(6,3) * t50;
t129 = t35 + t36;
t37 = -mrSges(6,1) * t104 + mrSges(6,3) * t51;
t38 = t104 * mrSges(7,1) - t51 * mrSges(7,2);
t128 = -t37 + t38;
t72 = -mrSges(5,1) * t103 + mrSges(5,2) * t100;
t127 = t72 - mrSges(4,1);
t70 = pkin(4) * t119 + t90;
t126 = t100 ^ 2 + t103 ^ 2;
t125 = Ifges(5,4) * t100;
t124 = Ifges(5,4) * t103;
t123 = t101 * t53;
t55 = t101 * t99 + t104 * t122;
t121 = t104 * t55;
t105 = cos(qJ(2));
t120 = t105 * t97;
t113 = t130 * t100;
t74 = t130 * t103;
t32 = -t113 * t98 - t74 * t96;
t34 = t113 * t96 - t98 * t74;
t117 = t32 ^ 2 + t34 ^ 2;
t116 = -Ifges(5,3) - Ifges(7,2) - Ifges(6,3);
t28 = -t100 * t55 - t103 * t120;
t29 = -t100 * t120 + t103 * t55;
t7 = -t98 * t28 + t29 * t96;
t9 = t28 * t96 + t29 * t98;
t115 = t32 * t7 + t34 * t9;
t85 = -pkin(4) * t103 - pkin(3);
t112 = -Ifges(5,5) * t118 + (Ifges(7,4) + Ifges(6,5)) * t51 + (Ifges(6,6) - Ifges(7,6)) * t50;
t4 = t20 * t98 - t30 * t96;
t110 = mrSges(5,1) * t100 + mrSges(5,2) * t103;
t109 = -t100 * t28 + t103 * t29;
t107 = pkin(8) ^ 2;
t95 = t104 ^ 2;
t93 = t101 ^ 2;
t91 = t97 ^ 2;
t89 = t93 * t107;
t88 = Ifges(5,5) * t100;
t87 = Ifges(5,6) * t103;
t84 = -pkin(4) * t98 - pkin(5);
t81 = t91 * t105 ^ 2;
t80 = pkin(4) * t96 + qJ(6);
t76 = Ifges(5,1) * t100 + t124;
t75 = Ifges(5,2) * t103 + t125;
t73 = -mrSges(4,1) * t104 + mrSges(4,2) * t101;
t69 = -mrSges(5,1) * t104 - mrSges(5,3) * t118;
t68 = mrSges(5,2) * t104 - mrSges(5,3) * t119;
t62 = t110 * t101;
t61 = Ifges(7,4) * t64;
t60 = Ifges(6,5) * t64;
t59 = Ifges(6,6) * t63;
t58 = Ifges(7,6) * t63;
t49 = -Ifges(5,5) * t104 + (Ifges(5,1) * t103 - t125) * t101;
t48 = -Ifges(5,6) * t104 + (-Ifges(5,2) * t100 + t124) * t101;
t39 = -t100 * t132 + t66;
t27 = Ifges(6,1) * t64 - Ifges(6,4) * t63;
t26 = Ifges(7,1) * t64 + Ifges(7,5) * t63;
t25 = Ifges(6,4) * t64 - Ifges(6,2) * t63;
t24 = Ifges(7,5) * t64 + Ifges(7,3) * t63;
t17 = pkin(5) * t63 - qJ(6) * t64 + t85;
t14 = -Ifges(6,1) * t51 - Ifges(6,4) * t50 - Ifges(6,5) * t104;
t13 = -Ifges(7,1) * t51 - Ifges(7,4) * t104 + Ifges(7,5) * t50;
t12 = -Ifges(6,4) * t51 - Ifges(6,2) * t50 - Ifges(6,6) * t104;
t11 = -Ifges(7,5) * t51 - Ifges(7,6) * t104 + Ifges(7,3) * t50;
t10 = pkin(5) * t50 + qJ(6) * t51 + t70;
t3 = pkin(5) * t104 - t4;
t2 = -qJ(6) * t104 + t5;
t1 = [m(2) + m(5) * (t28 ^ 2 + t29 ^ 2 + t52) + m(4) * (t55 ^ 2 + t52 + t81) + m(3) * (t102 ^ 2 * t91 + t99 ^ 2 + t81) + t136 * (t7 ^ 2 + t9 ^ 2 + t52); mrSges(4,3) * t121 + t28 * t69 + t29 * t68 + t129 * t9 + t128 * t7 + (-t102 * mrSges(3,2) + (mrSges(3,1) - t73) * t105) * t97 + (t101 * mrSges(4,3) + t134 + t62) * t53 + m(6) * (-t4 * t7 + t5 * t9 + t53 * t70) + m(7) * (t10 * t53 + t2 * t9 + t3 * t7) + m(5) * (t123 * pkin(8) + t28 * t39 + t29 * t40) + m(4) * (pkin(2) * t120 + (t121 + t123) * pkin(8)); -0.2e1 * pkin(2) * t73 + 0.2e1 * t10 * t15 + 0.2e1 * t70 * t16 + 0.2e1 * t2 * t35 + 0.2e1 * t3 * t38 + 0.2e1 * t5 * t36 + 0.2e1 * t4 * t37 + 0.2e1 * t39 * t69 + 0.2e1 * t40 * t68 + Ifges(3,3) - (t13 + t14) * t51 + (t11 - t12) * t50 + (t93 + t95) * mrSges(4,3) * t138 + (Ifges(4,1) * t101 - t100 * t48 + t103 * t49 + t62 * t138) * t101 + m(4) * (pkin(2) ^ 2 + t107 * t95 + t89) + m(5) * (t39 ^ 2 + t40 ^ 2 + t89) + m(6) * (t4 ^ 2 + t5 ^ 2 + t70 ^ 2) + m(7) * (t10 ^ 2 + t2 ^ 2 + t3 ^ 2) + ((Ifges(4,2) - t116) * t104 + (Ifges(5,6) * t100 + (2 * Ifges(4,4))) * t101 + t112) * t104; -t55 * mrSges(4,2) + t109 * mrSges(5,3) + (t127 + t133) * t53 + m(6) * (t53 * t85 + t115) + m(7) * (t17 * t53 + t115) + m(5) * (-pkin(3) * t53 + pkin(9) * t109) + t135 * (-t63 * t9 + t64 * t7); -pkin(3) * t62 + t10 * t22 + t17 * t15 + t85 * t16 + t70 * t23 - (t26 / 0.2e1 + t27 / 0.2e1) * t51 + (t24 / 0.2e1 - t25 / 0.2e1) * t50 + t129 * t34 + t128 * t32 + (-t88 / 0.2e1 - t87 / 0.2e1 - t60 / 0.2e1 + t59 / 0.2e1 - t61 / 0.2e1 - t58 / 0.2e1 + Ifges(4,6) - pkin(8) * mrSges(4,2)) * t104 + (t48 / 0.2e1 + pkin(9) * t68 + t40 * mrSges(5,3)) * t103 + (t49 / 0.2e1 - pkin(9) * t69 - t39 * mrSges(5,3)) * t100 + (Ifges(4,5) + t103 * t76 / 0.2e1 - t100 * t75 / 0.2e1 + t127 * pkin(8)) * t101 + m(5) * (-pkin(3) * t90 + (-t39 * t100 + t40 * t103) * pkin(9)) + m(6) * (-t32 * t4 + t34 * t5 + t70 * t85) + m(7) * (t10 * t17 + t2 * t34 + t3 * t32) + (t13 / 0.2e1 + t14 / 0.2e1 + t3 * mrSges(7,2) - t4 * mrSges(6,3)) * t64 + (t11 / 0.2e1 - t12 / 0.2e1 - t2 * mrSges(7,2) - t5 * mrSges(6,3)) * t63; -0.2e1 * pkin(3) * t72 + t100 * t76 + t103 * t75 + 0.2e1 * t17 * t22 + 0.2e1 * t85 * t23 + Ifges(4,3) + m(7) * (t17 ^ 2 + t117) + m(6) * (t85 ^ 2 + t117) + m(5) * (pkin(9) ^ 2 * t126 + pkin(3) ^ 2) + (t26 + t27) * t64 + (t24 - t25) * t63 + 0.2e1 * t126 * pkin(9) * mrSges(5,3) + 0.2e1 * (t32 * t64 - t34 * t63) * t135; t28 * mrSges(5,1) - t29 * mrSges(5,2) + (-mrSges(6,2) + mrSges(7,3)) * t9 + (-mrSges(6,1) - mrSges(7,1)) * t7 + m(7) * (t7 * t84 + t80 * t9) + (-t7 * t98 + t9 * t96) * t137; -Ifges(5,6) * t119 - t3 * mrSges(7,1) - t5 * mrSges(6,2) + t4 * mrSges(6,1) + t84 * t38 + t2 * mrSges(7,3) + t80 * t35 + m(7) * (t2 * t80 + t3 * t84) - t40 * mrSges(5,2) + t39 * mrSges(5,1) + t116 * t104 + (m(6) * (t4 * t98 + t5 * t96) + t96 * t36 + t98 * t37) * pkin(4) - t112; m(7) * (t32 * t84 + t34 * t80) - t32 * mrSges(7,1) + t34 * mrSges(7,3) - t34 * mrSges(6,2) - t32 * mrSges(6,1) + t60 + t58 - t59 + t88 + t87 + t61 - t110 * pkin(9) + (-t63 * t80 + t64 * t84) * mrSges(7,2) + (m(6) * (-t32 * t98 + t34 * t96) + (-t63 * t96 - t64 * t98) * mrSges(6,3)) * pkin(4); -0.2e1 * t84 * mrSges(7,1) + 0.2e1 * t80 * mrSges(7,3) + m(7) * (t80 ^ 2 + t84 ^ 2) - t116 + (0.2e1 * mrSges(6,1) * t98 - 0.2e1 * mrSges(6,2) * t96 + (t96 ^ 2 + t98 ^ 2) * t137) * pkin(4); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * t53; m(6) * t70 + m(7) * t10 + t134; m(6) * t85 + m(7) * t17 + t133; 0; t136; m(7) * t7; m(7) * t3 + t38; m(7) * t32 + t64 * mrSges(7,2); m(7) * t84 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
