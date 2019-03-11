% Calculate joint inertia matrix for
% S6RRPPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
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
% Datum: 2019-03-09 08:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP4_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:37:03
% EndTime: 2019-03-09 08:37:06
% DurationCPUTime: 1.04s
% Computational Cost: add. (1041->278), mult. (2068->363), div. (0->0), fcn. (1871->6), ass. (0->106)
t136 = Ifges(7,2) + Ifges(6,3);
t96 = sin(pkin(9));
t97 = cos(pkin(9));
t135 = t96 ^ 2 + t97 ^ 2;
t134 = 2 * pkin(7);
t133 = m(7) + m(6);
t119 = -mrSges(6,2) + mrSges(7,3);
t118 = mrSges(7,2) + mrSges(6,3);
t132 = -m(7) * pkin(5) - mrSges(7,1);
t131 = mrSges(6,1) - t132;
t130 = m(7) * qJ(6) + t119;
t129 = pkin(3) + pkin(4);
t128 = Ifges(4,4) * t96;
t127 = Ifges(4,4) * t97;
t126 = Ifges(5,5) * t96;
t125 = Ifges(5,5) * t97;
t101 = cos(qJ(2));
t124 = pkin(7) * t101;
t123 = t96 * mrSges(5,3);
t99 = sin(qJ(2));
t122 = t96 * t99;
t121 = t97 * t99;
t120 = -mrSges(6,1) - mrSges(7,1);
t117 = -pkin(8) + qJ(3);
t100 = cos(qJ(5));
t67 = -pkin(2) * t101 - qJ(3) * t99 - pkin(1);
t82 = t96 * t124;
t91 = t101 * pkin(3);
t15 = pkin(4) * t101 + t82 + t91 + (-pkin(8) * t99 - t67) * t97;
t38 = t124 * t97 + t67 * t96;
t31 = -qJ(4) * t101 + t38;
t17 = pkin(8) * t122 + t31;
t98 = sin(qJ(5));
t4 = t100 * t17 + t15 * t98;
t113 = t100 * t96;
t49 = -t113 * t99 + t121 * t98;
t33 = -mrSges(7,2) * t49 + mrSges(7,3) * t101;
t34 = -mrSges(6,2) * t101 - mrSges(6,3) * t49;
t116 = t33 + t34;
t60 = t100 * t97 + t96 * t98;
t50 = t60 * t99;
t35 = mrSges(6,1) * t101 - mrSges(6,3) * t50;
t36 = -mrSges(7,1) * t101 + mrSges(7,2) * t50;
t115 = t35 - t36;
t64 = mrSges(5,1) * t101 + mrSges(5,2) * t121;
t53 = mrSges(4,1) * t122 + mrSges(4,2) * t121;
t114 = t135 * qJ(3) ^ 2;
t111 = t117 * t96;
t70 = t117 * t97;
t26 = -t100 * t111 + t70 * t98;
t28 = t100 * t70 + t111 * t98;
t112 = t26 ^ 2 + t28 ^ 2;
t110 = qJ(4) * t96 + pkin(2);
t37 = t67 * t97 - t82;
t32 = -t37 + t91;
t107 = t31 * t97 + t32 * t96;
t106 = -t37 * t96 + t38 * t97;
t3 = t100 * t15 - t17 * t98;
t104 = (Ifges(7,4) + Ifges(6,5)) * t50 + (-Ifges(6,6) + Ifges(7,6)) * t49 + t136 * t101;
t51 = t129 * t97 + t110;
t76 = qJ(4) * t121;
t29 = -t76 - (-t129 * t96 - pkin(7)) * t99;
t103 = pkin(7) ^ 2;
t95 = t101 ^ 2;
t94 = t99 ^ 2;
t90 = t94 * t103;
t86 = t96 * mrSges(4,2);
t78 = mrSges(5,1) * t122;
t74 = Ifges(4,1) * t96 + t127;
t73 = Ifges(5,1) * t96 - t125;
t72 = Ifges(4,2) * t97 + t128;
t71 = -Ifges(5,3) * t97 + t126;
t69 = -t97 * mrSges(4,1) + t86;
t68 = -t97 * mrSges(5,1) - t123;
t66 = -pkin(3) * t97 - t110;
t65 = -mrSges(5,2) * t122 - mrSges(5,3) * t101;
t63 = -mrSges(4,1) * t101 - mrSges(4,3) * t121;
t62 = mrSges(4,2) * t101 - mrSges(4,3) * t122;
t61 = -t97 * t98 + t113;
t57 = Ifges(7,4) * t61;
t56 = Ifges(6,5) * t61;
t55 = Ifges(6,6) * t60;
t54 = Ifges(7,6) * t60;
t52 = -mrSges(5,3) * t121 + t78;
t48 = -t76 + (pkin(3) * t96 + pkin(7)) * t99;
t47 = -Ifges(4,5) * t101 + (Ifges(4,1) * t97 - t128) * t99;
t46 = -Ifges(5,4) * t101 + (Ifges(5,1) * t97 + t126) * t99;
t45 = -Ifges(4,6) * t101 + (-Ifges(4,2) * t96 + t127) * t99;
t44 = -Ifges(5,6) * t101 + (Ifges(5,3) * t96 + t125) * t99;
t24 = Ifges(6,1) * t61 - Ifges(6,4) * t60;
t23 = Ifges(7,1) * t61 + Ifges(7,5) * t60;
t22 = Ifges(6,4) * t61 - Ifges(6,2) * t60;
t21 = Ifges(7,5) * t61 + Ifges(7,3) * t60;
t20 = t60 * mrSges(6,1) + t61 * mrSges(6,2);
t19 = t60 * mrSges(7,1) - t61 * mrSges(7,3);
t13 = t49 * mrSges(6,1) + t50 * mrSges(6,2);
t12 = t49 * mrSges(7,1) - t50 * mrSges(7,3);
t11 = Ifges(6,1) * t50 - Ifges(6,4) * t49 + Ifges(6,5) * t101;
t10 = Ifges(7,1) * t50 + Ifges(7,4) * t101 + Ifges(7,5) * t49;
t9 = Ifges(6,4) * t50 - Ifges(6,2) * t49 + Ifges(6,6) * t101;
t8 = Ifges(7,5) * t50 + Ifges(7,6) * t101 + Ifges(7,3) * t49;
t7 = pkin(5) * t60 - qJ(6) * t61 + t51;
t5 = -pkin(5) * t49 + qJ(6) * t50 + t29;
t2 = -pkin(5) * t101 - t3;
t1 = qJ(6) * t101 + t4;
t6 = [0.2e1 * t1 * t33 - 0.2e1 * t5 * t12 - 0.2e1 * t29 * t13 + 0.2e1 * t2 * t36 + 0.2e1 * t3 * t35 + 0.2e1 * t31 * t65 + 0.2e1 * t32 * t64 + 0.2e1 * t4 * t34 + 0.2e1 * t37 * t63 + 0.2e1 * t38 * t62 + 0.2e1 * t48 * t52 + Ifges(2,3) + (t10 + t11) * t50 + (t8 - t9) * t49 + (t94 + t95) * mrSges(3,3) * t134 + m(3) * (pkin(1) ^ 2 + t103 * t95 + t90) + m(4) * (t37 ^ 2 + t38 ^ 2 + t90) + m(5) * (t31 ^ 2 + t32 ^ 2 + t48 ^ 2) + m(6) * (t29 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(5,2) + Ifges(4,3) + Ifges(3,2)) * t101 + t104) * t101 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t99 + t53 * t134 + (t46 + t47) * t97 + (t44 - t45) * t96 + ((2 * Ifges(3,4)) + (-Ifges(5,4) - Ifges(4,5)) * t97 + (Ifges(4,6) - Ifges(5,6)) * t96) * t101) * t99; -pkin(2) * t53 + t7 * t12 + t51 * t13 - t5 * t19 - t29 * t20 + t48 * t68 + t66 * t52 + (-t44 / 0.2e1 + t45 / 0.2e1) * t97 + (t47 / 0.2e1 + t46 / 0.2e1) * t96 + (t23 / 0.2e1 + t24 / 0.2e1) * t50 + (t21 / 0.2e1 - t22 / 0.2e1) * t49 + t116 * t28 - t115 * t26 + (t56 / 0.2e1 - t55 / 0.2e1 + t57 / 0.2e1 + t54 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2) + (-Ifges(4,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * t97 + (-Ifges(4,5) / 0.2e1 - Ifges(5,4) / 0.2e1) * t96) * t101 + t106 * mrSges(4,3) + t107 * mrSges(5,2) + (t10 / 0.2e1 + t11 / 0.2e1 + t2 * mrSges(7,2) - t3 * mrSges(6,3)) * t61 + (t8 / 0.2e1 - t9 / 0.2e1 - t1 * mrSges(7,2) - t4 * mrSges(6,3)) * t60 + ((t62 + t65) * t97 + (-t63 + t64) * t96) * qJ(3) + (Ifges(3,5) + (t73 / 0.2e1 + t74 / 0.2e1) * t97 + (t71 / 0.2e1 - t72 / 0.2e1) * t96 + (t69 - mrSges(3,1)) * pkin(7)) * t99 + m(4) * (-pkin(2) * pkin(7) * t99 + qJ(3) * t106) + m(5) * (qJ(3) * t107 + t48 * t66) + m(6) * (-t26 * t3 + t28 * t4 - t29 * t51) + m(7) * (t1 * t28 + t2 * t26 - t5 * t7); -0.2e1 * pkin(2) * t69 + 0.2e1 * t7 * t19 + 0.2e1 * t51 * t20 + 0.2e1 * t66 * t68 + Ifges(3,3) + (-t71 + t72) * t97 + (t73 + t74) * t96 + (0.2e1 * t118 * t26 + t23 + t24) * t61 + (-0.2e1 * t118 * t28 + t21 - t22) * t60 + m(5) * (t66 ^ 2 + t114) + m(4) * (pkin(2) ^ 2 + t114) + m(6) * (t51 ^ 2 + t112) + m(7) * (t7 ^ 2 + t112) + 0.2e1 * (mrSges(5,2) + mrSges(4,3)) * qJ(3) * t135; t78 + (m(4) * pkin(7) - mrSges(5,3) * t97) * t99 + t119 * t50 + t120 * t49 + m(5) * t48 + m(6) * t29 + m(7) * t5 + t53; -m(4) * pkin(2) - t123 + t86 + (-mrSges(5,1) - mrSges(4,1)) * t97 + t119 * t61 + t120 * t60 + m(5) * t66 - m(6) * t51 - m(7) * t7; m(4) + m(5) + t133; t116 * t98 + t115 * t100 + m(7) * (t1 * t98 - t100 * t2) + m(6) * (t100 * t3 + t4 * t98) + m(5) * t32 + t64; (m(5) * qJ(3) + mrSges(5,2)) * t96 + t133 * (-t100 * t26 + t28 * t98) + t118 * (-t100 * t61 - t60 * t98); 0; m(5) + 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t100 ^ 2 + t98 ^ 2); t3 * mrSges(6,1) - t4 * mrSges(6,2) - pkin(5) * t36 + m(7) * (-pkin(5) * t2 + qJ(6) * t1) + qJ(6) * t33 + t1 * mrSges(7,3) - t2 * mrSges(7,1) + t104; t57 + t54 + t56 - t55 + (-pkin(5) * t61 - qJ(6) * t60) * mrSges(7,2) + t130 * t28 - t131 * t26; 0; t100 * t131 + t130 * t98; 0.2e1 * pkin(5) * mrSges(7,1) + 0.2e1 * qJ(6) * mrSges(7,3) + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t136; m(7) * t2 + t36; m(7) * t26 + t61 * mrSges(7,2); 0; -m(7) * t100; t132; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t6(1) t6(2) t6(4) t6(7) t6(11) t6(16); t6(2) t6(3) t6(5) t6(8) t6(12) t6(17); t6(4) t6(5) t6(6) t6(9) t6(13) t6(18); t6(7) t6(8) t6(9) t6(10) t6(14) t6(19); t6(11) t6(12) t6(13) t6(14) t6(15) t6(20); t6(16) t6(17) t6(18) t6(19) t6(20) t6(21);];
Mq  = res;
