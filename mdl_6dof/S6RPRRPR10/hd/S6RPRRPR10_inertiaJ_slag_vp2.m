% Calculate joint inertia matrix for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:29
% EndTime: 2019-03-09 05:34:32
% DurationCPUTime: 1.16s
% Computational Cost: add. (967->278), mult. (1794->372), div. (0->0), fcn. (1526->6), ass. (0->107)
t129 = Ifges(6,2) + Ifges(5,3);
t80 = sin(qJ(4));
t83 = cos(qJ(4));
t106 = t80 ^ 2 + t83 ^ 2;
t86 = (-pkin(1) - pkin(7));
t128 = -2 * t86;
t127 = (mrSges(6,2) + mrSges(5,3)) * t106;
t84 = cos(qJ(3));
t111 = t83 * t84;
t113 = t80 * t84;
t81 = sin(qJ(3));
t126 = Ifges(6,6) * t113 + t129 * t81 + (Ifges(6,4) + Ifges(5,5)) * t111;
t125 = 2 * qJ(2);
t124 = pkin(4) + pkin(5);
t123 = pkin(8) - pkin(9);
t121 = Ifges(5,4) * t80;
t120 = Ifges(5,4) * t83;
t119 = Ifges(6,5) * t80;
t118 = Ifges(6,5) * t83;
t117 = Ifges(5,6) * t81;
t116 = Ifges(6,6) * t81;
t79 = sin(qJ(6));
t82 = cos(qJ(6));
t44 = -qJ(5) * t79 - t124 * t82;
t115 = t44 * mrSges(7,1);
t45 = qJ(5) * t82 - t124 * t79;
t114 = t45 * mrSges(7,2);
t112 = t81 * t86;
t49 = -mrSges(5,1) * t83 + mrSges(5,2) * t80;
t110 = mrSges(4,1) - t49;
t39 = -t79 * t83 + t80 * t82;
t29 = t39 * t84;
t93 = t79 * t80 + t82 * t83;
t30 = t93 * t84;
t109 = -Ifges(7,5) * t30 - Ifges(7,6) * t29;
t46 = pkin(3) * t81 - pkin(8) * t84 + qJ(2);
t19 = t83 * t112 + t80 * t46;
t108 = t106 * pkin(8) * t81;
t107 = t106 * pkin(8) ^ 2;
t75 = t81 ^ 2;
t77 = t84 ^ 2;
t105 = t77 + t75;
t104 = qJ(5) * t83;
t12 = t81 * qJ(5) + t19;
t103 = t105 * mrSges(4,3);
t42 = -t81 * mrSges(6,1) + mrSges(6,2) * t111;
t102 = qJ(5) * t80 + pkin(3);
t56 = t80 * t112;
t18 = t46 * t83 - t56;
t100 = t80 * mrSges(5,1) + t83 * mrSges(5,2);
t99 = t80 * mrSges(6,1) - t83 * mrSges(6,3);
t27 = t39 * t81;
t28 = t93 * t81;
t98 = t27 * mrSges(7,1) - t28 * mrSges(7,2);
t97 = t82 * mrSges(7,1) - t79 * mrSges(7,2);
t96 = -pkin(4) * t80 + t104;
t13 = -pkin(4) * t81 - t18;
t95 = t12 * t83 + t13 * t80;
t94 = -t18 * t80 + t19 * t83;
t54 = t123 * t80;
t55 = t123 * t83;
t14 = t54 * t82 - t55 * t79;
t15 = t54 * t79 + t55 * t82;
t35 = Ifges(7,6) * t93;
t36 = Ifges(7,5) * t39;
t92 = t14 * mrSges(7,1) - t15 * mrSges(7,2) - t35 + t36;
t40 = -t81 * mrSges(5,2) - mrSges(5,3) * t113;
t41 = mrSges(5,1) * t81 - mrSges(5,3) * t111;
t43 = -mrSges(6,2) * t113 + t81 * mrSges(6,3);
t91 = (t40 + t43) * t83 + (-t41 + t42) * t80;
t6 = t56 + (-pkin(9) * t84 - t46) * t83 - t124 * t81;
t7 = pkin(9) * t113 + t12;
t1 = t6 * t82 - t7 * t79;
t2 = t6 * t79 + t7 * t82;
t90 = t1 * mrSges(7,1) - t2 * mrSges(7,2) - Ifges(7,3) * t81 - t109;
t89 = m(6) * t96 - t100 - t99;
t87 = qJ(2) ^ 2;
t78 = t86 ^ 2;
t71 = Ifges(6,4) * t80;
t70 = Ifges(5,5) * t80;
t68 = Ifges(5,6) * t83;
t65 = t77 * t86;
t64 = t77 * t78;
t53 = Ifges(5,1) * t80 + t120;
t52 = Ifges(6,1) * t80 - t118;
t51 = Ifges(5,2) * t83 + t121;
t50 = -Ifges(6,3) * t83 + t119;
t48 = -mrSges(6,1) * t83 - mrSges(6,3) * t80;
t47 = -pkin(4) * t83 - t102;
t34 = t124 * t83 + t102;
t33 = t100 * t84;
t32 = t99 * t84;
t26 = Ifges(5,5) * t81 + (Ifges(5,1) * t83 - t121) * t84;
t25 = Ifges(6,4) * t81 + (Ifges(6,1) * t83 + t119) * t84;
t24 = t117 + (-Ifges(5,2) * t80 + t120) * t84;
t23 = t116 + (Ifges(6,3) * t80 + t118) * t84;
t20 = (-t86 - t96) * t84;
t17 = -mrSges(7,1) * t81 - mrSges(7,3) * t30;
t16 = mrSges(7,2) * t81 + mrSges(7,3) * t29;
t11 = (-t124 * t80 + t104 + t86) * t84;
t10 = Ifges(7,1) * t39 - Ifges(7,4) * t93;
t9 = Ifges(7,4) * t39 - Ifges(7,2) * t93;
t8 = mrSges(7,1) * t93 + mrSges(7,2) * t39;
t5 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t4 = Ifges(7,1) * t30 + Ifges(7,4) * t29 - Ifges(7,5) * t81;
t3 = Ifges(7,4) * t30 + Ifges(7,2) * t29 - Ifges(7,6) * t81;
t21 = [-(2 * pkin(1) * mrSges(3,2)) + (mrSges(3,3) * t125) + 0.2e1 * t1 * t17 + 0.2e1 * t11 * t5 + 0.2e1 * t12 * t43 + 0.2e1 * t13 * t42 + 0.2e1 * t2 * t16 + 0.2e1 * t18 * t41 + 0.2e1 * t19 * t40 + 0.2e1 * t20 * t32 + t29 * t3 + t30 * t4 + Ifges(3,1) + Ifges(2,3) + t103 * t128 + (mrSges(4,1) * t125 + (Ifges(4,2) + Ifges(7,3)) * t81 + t109 + t126) * t81 + ((mrSges(4,2) * t125) + Ifges(4,1) * t84 - 0.2e1 * Ifges(4,4) * t81 + t33 * t128 + (t25 + t26) * t83 + (t23 - t24 - t117) * t80) * t84 + m(4) * (t75 * t78 + t64 + t87) + (m(3) * (pkin(1) ^ 2 + t87)) + m(5) * (t18 ^ 2 + t19 ^ 2 + t64) + m(6) * (t12 ^ 2 + t13 ^ 2 + t20 ^ 2) + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2); -(m(3) * pkin(1)) + t28 * t16 + t27 * t17 + mrSges(3,2) - t103 + (-t32 - t33 + t5) * t84 + t91 * t81 + m(7) * (t1 * t27 + t11 * t84 + t2 * t28) + m(6) * (-t20 * t84 + t95 * t81) + m(5) * (t94 * t81 + t65) + m(4) * (t75 * t86 + t65); m(3) + m(7) * (t27 ^ 2 + t28 ^ 2 + t77) + m(4) * t105 + 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * (t106 * t75 + t77); t47 * t32 + t30 * t10 / 0.2e1 - pkin(3) * t33 + t34 * t5 - t93 * t3 / 0.2e1 + t39 * t4 / 0.2e1 + t29 * t9 / 0.2e1 + t11 * t8 + t15 * t16 + t14 * t17 + (-t1 * t39 - t2 * t93) * mrSges(7,3) + m(7) * (t1 * t14 + t11 * t34 + t15 * t2) + (t71 / 0.2e1 + t70 / 0.2e1 + t68 / 0.2e1 - t36 / 0.2e1 + t35 / 0.2e1 - Ifges(4,6) - (t86 * mrSges(4,2))) * t81 + (-t116 / 0.2e1 - t23 / 0.2e1 + t24 / 0.2e1 + t19 * mrSges(5,3) + t12 * mrSges(6,2)) * t83 + (t25 / 0.2e1 + t26 / 0.2e1 + t13 * mrSges(6,2) - t18 * mrSges(5,3)) * t80 + (m(5) * t94 + m(6) * t95 + t91) * pkin(8) + (Ifges(4,5) + (t52 / 0.2e1 + t53 / 0.2e1) * t83 + (t50 / 0.2e1 - t51 / 0.2e1) * t80 + (m(5) * pkin(3) + t110) * t86) * t84 + (m(6) * t47 + t48) * t20; (-t27 * t39 - t28 * t93) * mrSges(7,3) + (-t48 + t8 + t110) * t84 + m(6) * (-t47 * t84 + t108) + m(7) * (t14 * t27 + t15 * t28 + t34 * t84) + m(5) * (pkin(3) * t84 + t108) + (-mrSges(4,2) + t127) * t81; -0.2e1 * pkin(3) * t49 + t39 * t10 + 0.2e1 * t34 * t8 - t93 * t9 + 0.2e1 * t47 * t48 + Ifges(4,3) + (t51 - t50) * t83 + (t52 + t53) * t80 + 0.2e1 * (-t14 * t39 - t15 * t93) * mrSges(7,3) + m(7) * (t14 ^ 2 + t15 ^ 2 + t34 ^ 2) + m(6) * (t47 ^ 2 + t107) + m(5) * (pkin(3) ^ 2 + t107) + 0.2e1 * pkin(8) * t127; m(7) * (t1 * t44 + t2 * t45) + m(6) * (-pkin(4) * t13 + qJ(5) * t12) - t90 - pkin(4) * t42 + qJ(5) * t43 + t44 * t17 + t45 * t16 + t18 * mrSges(5,1) - t19 * mrSges(5,2) + t12 * mrSges(6,3) - t13 * mrSges(6,1) - Ifges(5,6) * t113 + t126; m(7) * (t27 * t44 + t28 * t45) + t89 * t81 - t98; m(7) * (t14 * t44 + t15 * t45) + t68 + t71 - Ifges(6,6) * t83 + t70 + (-t39 * t44 - t45 * t93) * mrSges(7,3) + t96 * mrSges(6,2) + t89 * pkin(8) - t92; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t115 + 0.2e1 * t114 + 0.2e1 * qJ(5) * mrSges(6,3) + Ifges(7,3) + m(6) * (pkin(4) ^ 2 + qJ(5) ^ 2) + m(7) * (t44 ^ 2 + t45 ^ 2) + t129; t79 * t16 + t82 * t17 + m(7) * (t1 * t82 + t2 * t79) + m(6) * t13 + t42; m(6) * t80 * t81 + m(7) * (t27 * t82 + t28 * t79); m(7) * (t14 * t82 + t15 * t79) + (m(6) * pkin(8) + mrSges(6,2)) * t80 + (-t39 * t82 - t79 * t93) * mrSges(7,3); -mrSges(6,1) - m(6) * pkin(4) + m(7) * (t44 * t82 + t45 * t79) - t97; m(6) + m(7) * (t79 ^ 2 + t82 ^ 2); t90; t98; t92; -Ifges(7,3) - t114 + t115; t97; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
