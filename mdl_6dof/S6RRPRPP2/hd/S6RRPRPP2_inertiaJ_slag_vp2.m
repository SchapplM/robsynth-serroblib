% Calculate joint inertia matrix for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 16:57
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:57:24
% EndTime: 2018-11-23 16:57:25
% DurationCPUTime: 0.93s
% Computational Cost: add. (1151->269), mult. (2165->344), div. (0->0), fcn. (2130->6), ass. (0->96)
t124 = Ifges(6,2) + Ifges(5,3);
t79 = sin(qJ(4));
t81 = cos(qJ(4));
t123 = t79 ^ 2 + t81 ^ 2;
t122 = 0.2e1 * t123;
t121 = m(4) * pkin(2);
t102 = -qJ(3) - pkin(7);
t82 = cos(qJ(2));
t51 = t102 * t82;
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t80 = sin(qJ(2));
t92 = t102 * t80;
t28 = -t51 * t77 - t78 * t92;
t120 = t28 ^ 2;
t119 = -2 * mrSges(7,3);
t118 = 0.2e1 * t28;
t65 = -pkin(2) * t82 - pkin(1);
t117 = 0.2e1 * t65;
t116 = m(6) + m(7);
t83 = -pkin(4) - pkin(5);
t115 = pkin(2) * t77;
t114 = pkin(2) * t78;
t113 = Ifges(5,4) * t79;
t112 = Ifges(5,4) * t81;
t111 = Ifges(7,4) * t79;
t110 = Ifges(7,4) * t81;
t109 = Ifges(6,5) * t79;
t108 = Ifges(6,5) * t81;
t45 = t77 * t80 - t78 * t82;
t107 = Ifges(7,5) * t45;
t46 = t77 * t82 + t78 * t80;
t106 = t46 * t79;
t105 = t46 * t81;
t104 = mrSges(6,2) - mrSges(7,3);
t103 = -Ifges(5,6) - Ifges(7,6);
t21 = -mrSges(5,2) * t45 - mrSges(5,3) * t106;
t25 = -mrSges(6,2) * t106 + mrSges(6,3) * t45;
t101 = t21 + t25;
t23 = mrSges(5,1) * t45 - mrSges(5,3) * t105;
t24 = -t45 * mrSges(6,1) + mrSges(6,2) * t105;
t100 = -t23 + t24;
t19 = pkin(3) * t45 - pkin(8) * t46 + t65;
t30 = -t78 * t51 + t77 * t92;
t7 = t79 * t19 + t81 * t30;
t63 = pkin(8) + t115;
t99 = t123 * t63 ^ 2;
t66 = t79 * qJ(5);
t98 = t81 * pkin(4) + t66;
t49 = t81 * mrSges(7,1) + t79 * mrSges(7,2);
t96 = t80 ^ 2 + t82 ^ 2;
t95 = qJ(5) * t81;
t94 = qJ(6) * t46;
t93 = -qJ(6) + t63;
t3 = t45 * qJ(5) + t7;
t64 = -pkin(3) - t114;
t26 = t79 * t30;
t6 = t19 * t81 - t26;
t17 = -mrSges(7,1) * t106 + mrSges(7,2) * t105;
t22 = -t45 * mrSges(7,1) - mrSges(7,3) * t105;
t90 = Ifges(6,6) * t106 + t124 * t45 + (Ifges(6,4) + Ifges(5,5)) * t105;
t89 = t79 * mrSges(5,1) + t81 * mrSges(5,2);
t88 = t79 * mrSges(6,1) - t81 * mrSges(6,3);
t87 = pkin(4) * t79 - t95;
t42 = t64 - t98;
t84 = qJ(5) ^ 2;
t71 = Ifges(6,4) * t79;
t70 = Ifges(5,5) * t79;
t69 = Ifges(5,6) * t81;
t57 = Ifges(5,1) * t79 + t112;
t56 = Ifges(6,1) * t79 - t108;
t55 = Ifges(7,1) * t79 - t110;
t54 = Ifges(5,2) * t81 + t113;
t53 = -Ifges(7,2) * t81 + t111;
t52 = -Ifges(6,3) * t81 + t109;
t50 = -t81 * mrSges(5,1) + t79 * mrSges(5,2);
t48 = -t81 * mrSges(6,1) - t79 * mrSges(6,3);
t44 = t93 * t81;
t43 = t93 * t79;
t39 = t46 * mrSges(4,2);
t31 = pkin(5) * t81 - t42;
t20 = mrSges(7,2) * t45 + mrSges(7,3) * t106;
t18 = t89 * t46;
t16 = t88 * t46;
t14 = Ifges(5,5) * t45 + (Ifges(5,1) * t81 - t113) * t46;
t13 = Ifges(6,4) * t45 + (Ifges(6,1) * t81 + t109) * t46;
t12 = -t107 + (Ifges(7,1) * t81 + t111) * t46;
t11 = Ifges(5,6) * t45 + (-Ifges(5,2) * t79 + t112) * t46;
t10 = -Ifges(7,6) * t45 + (Ifges(7,2) * t79 + t110) * t46;
t9 = Ifges(6,6) * t45 + (Ifges(6,3) * t79 + t108) * t46;
t8 = t46 * t87 + t28;
t5 = (t79 * t83 + t95) * t46 - t28;
t4 = -pkin(4) * t45 - t6;
t2 = t79 * t94 + t3;
t1 = t26 + (-t19 - t94) * t81 + t83 * t45;
t15 = [-0.2e1 * pkin(1) * (-t82 * mrSges(3,1) + t80 * mrSges(3,2)) + t80 * (Ifges(3,1) * t80 + Ifges(3,4) * t82) + t82 * (Ifges(3,4) * t80 + Ifges(3,2) * t82) + t39 * t117 + 0.2e1 * t2 * t20 + 0.2e1 * t7 * t21 + 0.2e1 * t1 * t22 + 0.2e1 * t6 * t23 + 0.2e1 * t4 * t24 + 0.2e1 * t3 * t25 + t18 * t118 + 0.2e1 * t8 * t16 + 0.2e1 * t5 * t17 + Ifges(2,3) + 0.2e1 * t96 * pkin(7) * mrSges(3,3) + (mrSges(4,1) * t117 - 0.2e1 * t30 * mrSges(4,3) + (Ifges(7,3) + Ifges(4,2)) * t45 + t90) * t45 + m(3) * (pkin(7) ^ 2 * t96 + pkin(1) ^ 2) + m(4) * (t30 ^ 2 + t65 ^ 2 + t120) + m(5) * (t6 ^ 2 + t7 ^ 2 + t120) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (mrSges(4,3) * t118 + Ifges(4,1) * t46 - 0.2e1 * Ifges(4,4) * t45 + (t12 + t13 + t14 - t107) * t81 + (t103 * t45 + t10 - t11 + t9) * t79) * t46; -t30 * mrSges(4,2) + Ifges(3,5) * t80 + Ifges(3,6) * t82 + t42 * t16 + t31 * t17 + t64 * t18 + t44 * t20 + t43 * t22 + t8 * t48 + t5 * t49 + (t50 - mrSges(4,1)) * t28 + (-t80 * mrSges(3,1) - t82 * mrSges(3,2)) * pkin(7) + (t71 / 0.2e1 + t70 / 0.2e1 + t69 / 0.2e1 - Ifges(4,6) - mrSges(4,3) * t115) * t45 + (t7 * mrSges(5,3) - t2 * mrSges(7,3) + t3 * mrSges(6,2) - t9 / 0.2e1 - t10 / 0.2e1 + t11 / 0.2e1 + t101 * t63 + (Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t45) * t81 + (-t6 * mrSges(5,3) - t1 * mrSges(7,3) + t4 * mrSges(6,2) - t107 / 0.2e1 + t12 / 0.2e1 + t13 / 0.2e1 + t14 / 0.2e1 + t100 * t63) * t79 + m(5) * (t28 * t64 + (-t6 * t79 + t7 * t81) * t63) + m(6) * (t42 * t8 + (t3 * t81 + t4 * t79) * t63) + m(7) * (t1 * t43 + t2 * t44 + t31 * t5) + (-t28 * t78 + t30 * t77) * t121 + (-mrSges(4,3) * t114 + Ifges(4,5) + (t55 / 0.2e1 + t56 / 0.2e1 + t57 / 0.2e1) * t81 + (t52 / 0.2e1 + t53 / 0.2e1 - t54 / 0.2e1) * t79) * t46; 0.2e1 * t31 * t49 + 0.2e1 * t42 * t48 + 0.2e1 * t64 * t50 + Ifges(3,3) + Ifges(4,3) + (t119 * t44 - t52 - t53 + t54) * t81 + (t119 * t43 + t55 + t56 + t57) * t79 + m(6) * (t42 ^ 2 + t99) + m(7) * (t31 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t64 ^ 2 + t99) + (mrSges(6,2) + mrSges(5,3)) * t63 * t122 + (0.2e1 * mrSges(4,1) * t78 - 0.2e1 * mrSges(4,2) * t77 + (t77 ^ 2 + t78 ^ 2) * t121) * pkin(2); t45 * mrSges(4,1) + t39 + (-t22 - t100) * t81 + (t20 + t101) * t79 + m(6) * (t3 * t79 - t4 * t81) + m(7) * (-t1 * t81 + t2 * t79) + m(5) * (t6 * t81 + t7 * t79) + m(4) * t65; m(7) * (-t43 * t81 + t44 * t79); m(4) + (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t122; t6 * mrSges(5,1) - t4 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t3 * mrSges(6,3) + Ifges(7,3) * t45 - pkin(4) * t24 + t83 * t22 + (t20 + t25) * qJ(5) + m(6) * (-pkin(4) * t4 + qJ(5) * t3) + m(7) * (qJ(5) * t2 + t1 * t83) + (-Ifges(7,5) * t81 + t103 * t79) * t46 + t90; m(7) * (qJ(5) * t44 + t43 * t83) - t43 * mrSges(7,1) + t44 * mrSges(7,2) + t71 + t70 + t69 + (-pkin(4) * mrSges(6,2) - t83 * mrSges(7,3) - Ifges(7,5)) * t79 + (qJ(5) * t104 - Ifges(6,6) + Ifges(7,6)) * t81 + (-m(6) * t87 - t88 - t89) * t63; (mrSges(5,1) + mrSges(6,1)) * t81 + (-mrSges(5,2) + mrSges(6,3)) * t79 + m(6) * t98 + m(7) * (-t81 * t83 + t66) + t49; 0.2e1 * pkin(4) * mrSges(6,1) - 0.2e1 * t83 * mrSges(7,1) + Ifges(7,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t84) + m(7) * (t83 ^ 2 + t84) + t124; m(6) * t4 + m(7) * t1 + t22 + t24; m(7) * t43 + (m(6) * t63 + t104) * t79; -t116 * t81; -m(6) * pkin(4) + m(7) * t83 - mrSges(6,1) - mrSges(7,1); t116; m(7) * t5 + t17; m(7) * t31 + t49; 0; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
