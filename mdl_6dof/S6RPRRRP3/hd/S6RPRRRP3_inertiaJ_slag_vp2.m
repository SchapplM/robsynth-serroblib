% Calculate joint inertia matrix for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2018-11-23 16:25
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:25:12
% EndTime: 2018-11-23 16:25:13
% DurationCPUTime: 1.06s
% Computational Cost: add. (1299->279), mult. (2520->377), div. (0->0), fcn. (2335->8), ass. (0->101)
t90 = sin(qJ(5));
t91 = sin(qJ(4));
t93 = cos(qJ(5));
t94 = cos(qJ(4));
t61 = t90 * t94 + t91 * t93;
t92 = sin(qJ(3));
t48 = t61 * t92;
t60 = t90 * t91 - t93 * t94;
t50 = t60 * t92;
t16 = t48 * mrSges(7,1) + t50 * mrSges(7,3);
t17 = t48 * mrSges(6,1) - t50 * mrSges(6,2);
t101 = -t16 - t17;
t106 = mrSges(5,1) * t91 + mrSges(5,2) * t94;
t54 = t106 * t92;
t132 = -t54 + t101;
t88 = sin(pkin(10));
t75 = pkin(1) * t88 + pkin(7);
t131 = 0.2e1 * t75;
t130 = mrSges(7,2) + mrSges(6,3);
t129 = 0.2e1 * mrSges(7,3);
t128 = -pkin(9) - pkin(8);
t95 = cos(qJ(3));
t127 = pkin(3) * t95;
t126 = pkin(8) * t92;
t121 = t92 * t94;
t89 = cos(pkin(10));
t77 = -pkin(1) * t89 - pkin(2);
t55 = -t126 + t77 - t127;
t52 = t94 * t55;
t14 = -pkin(9) * t121 + t52 + (-t75 * t91 - pkin(4)) * t95;
t122 = t91 * t92;
t123 = t75 * t95;
t21 = t94 * t123 + t91 * t55;
t18 = -pkin(9) * t122 + t21;
t6 = t90 * t14 + t93 * t18;
t125 = Ifges(5,4) * t91;
t124 = Ifges(5,4) * t94;
t70 = t92 * t75;
t119 = -Ifges(6,3) - Ifges(7,2);
t34 = -mrSges(7,2) * t48 - mrSges(7,3) * t95;
t35 = mrSges(6,2) * t95 - mrSges(6,3) * t48;
t118 = t34 + t35;
t36 = -mrSges(6,1) * t95 + mrSges(6,3) * t50;
t37 = t95 * mrSges(7,1) - t50 * mrSges(7,2);
t117 = -t36 + t37;
t66 = -mrSges(5,1) * t94 + mrSges(5,2) * t91;
t116 = t66 - mrSges(4,1);
t53 = pkin(4) * t122 + t70;
t115 = t91 ^ 2 + t94 ^ 2;
t85 = t92 ^ 2;
t87 = t95 ^ 2;
t114 = t85 + t87;
t111 = t128 * t91;
t71 = t128 * t94;
t31 = -t93 * t111 - t71 * t90;
t33 = t90 * t111 - t93 * t71;
t113 = t31 ^ 2 + t33 ^ 2;
t112 = -Ifges(5,3) + t119;
t80 = -pkin(4) * t94 - pkin(3);
t110 = t115 * mrSges(5,3);
t109 = t31 * t48 - t33 * t50;
t107 = (Ifges(7,4) + Ifges(6,5)) * t50 + (Ifges(6,6) - Ifges(7,6)) * t48;
t105 = -pkin(5) * t48 - qJ(6) * t50;
t5 = t14 * t93 - t18 * t90;
t20 = -t91 * t123 + t52;
t104 = -t20 * t91 + t21 * t94;
t102 = (mrSges(6,1) * t93 - mrSges(6,2) * t90) * pkin(4);
t2 = -qJ(6) * t95 + t6;
t3 = pkin(5) * t95 - t5;
t100 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) - t107;
t56 = Ifges(7,6) * t60;
t57 = Ifges(6,6) * t60;
t58 = Ifges(6,5) * t61;
t59 = Ifges(7,4) * t61;
t99 = t56 - t57 + t58 + t59 + (-mrSges(6,2) + mrSges(7,3)) * t33 + (-mrSges(6,1) - mrSges(7,1)) * t31;
t83 = Ifges(5,5) * t91;
t82 = Ifges(5,6) * t94;
t79 = -pkin(4) * t93 - pkin(5);
t76 = pkin(4) * t90 + qJ(6);
t74 = t75 ^ 2;
t72 = Ifges(5,5) * t121;
t69 = t85 * t74;
t68 = Ifges(5,1) * t91 + t124;
t67 = Ifges(5,2) * t94 + t125;
t64 = -mrSges(5,1) * t95 - mrSges(5,3) * t121;
t63 = mrSges(5,2) * t95 - mrSges(5,3) * t122;
t46 = -Ifges(5,5) * t95 + (Ifges(5,1) * t94 - t125) * t92;
t45 = -Ifges(5,6) * t95 + (-Ifges(5,2) * t91 + t124) * t92;
t27 = Ifges(6,1) * t61 - Ifges(6,4) * t60;
t26 = Ifges(7,1) * t61 + Ifges(7,5) * t60;
t25 = Ifges(6,4) * t61 - Ifges(6,2) * t60;
t24 = Ifges(7,5) * t61 + Ifges(7,3) * t60;
t23 = mrSges(6,1) * t60 + mrSges(6,2) * t61;
t22 = mrSges(7,1) * t60 - mrSges(7,3) * t61;
t19 = pkin(5) * t60 - qJ(6) * t61 + t80;
t12 = -Ifges(6,1) * t50 - Ifges(6,4) * t48 - Ifges(6,5) * t95;
t11 = -Ifges(7,1) * t50 - Ifges(7,4) * t95 + Ifges(7,5) * t48;
t10 = -Ifges(6,4) * t50 - Ifges(6,2) * t48 - Ifges(6,6) * t95;
t9 = -Ifges(7,5) * t50 - Ifges(7,6) * t95 + Ifges(7,3) * t48;
t7 = -t105 + t53;
t1 = [0.2e1 * t7 * t16 + 0.2e1 * t53 * t17 + 0.2e1 * t2 * t34 + 0.2e1 * t20 * t64 + 0.2e1 * t21 * t63 + 0.2e1 * t3 * t37 + 0.2e1 * t6 * t35 + 0.2e1 * t5 * t36 + Ifges(2,3) + Ifges(3,3) - (t11 + t12) * t50 + (-t10 + t9) * t48 + (0.2e1 * mrSges(4,2) * t77 + Ifges(4,1) * t92 + t54 * t131 - t45 * t91 + t46 * t94) * t92 + (-0.2e1 * t77 * mrSges(4,1) - t72 + (Ifges(4,2) - t112) * t95 + (Ifges(5,6) * t91 + (2 * Ifges(4,4))) * t92 + t107) * t95 + m(4) * (t74 * t87 + t77 ^ 2 + t69) + m(5) * (t20 ^ 2 + t21 ^ 2 + t69) + m(6) * (t5 ^ 2 + t53 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + m(3) * (t88 ^ 2 + t89 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (mrSges(3,1) * t89 - mrSges(3,2) * t88) * pkin(1) + t114 * mrSges(4,3) * t131; -t118 * t50 + t117 * t48 + t132 * t95 + m(7) * (-t2 * t50 + t3 * t48 - t7 * t95) + m(6) * (-t48 * t5 - t50 * t6 - t53 * t95) + (t94 * t63 - t91 * t64 + m(5) * (t104 - t123)) * t92; m(3) + m(5) * (t115 * t85 + t87) + m(4) * t114 + (m(7) + m(6)) * (t48 ^ 2 + t50 ^ 2 + t87); -pkin(3) * t54 + t19 * t16 + t80 * t17 + t7 * t22 + t53 * t23 + (-t83 / 0.2e1 - t82 / 0.2e1 - t58 / 0.2e1 + t57 / 0.2e1 - t59 / 0.2e1 - t56 / 0.2e1 + Ifges(4,6) - t75 * mrSges(4,2)) * t95 - (t26 / 0.2e1 + t27 / 0.2e1) * t50 + (-t25 / 0.2e1 + t24 / 0.2e1) * t48 + t118 * t33 + t117 * t31 + (t45 / 0.2e1 + pkin(8) * t63 + t21 * mrSges(5,3)) * t94 + (t46 / 0.2e1 - pkin(8) * t64 - t20 * mrSges(5,3)) * t91 + (Ifges(4,5) + t94 * t68 / 0.2e1 - t91 * t67 / 0.2e1 + t116 * t75) * t92 + m(5) * (-pkin(3) * t70 + t104 * pkin(8)) + m(6) * (-t31 * t5 + t33 * t6 + t53 * t80) + m(7) * (t19 * t7 + t2 * t33 + t3 * t31) + (t11 / 0.2e1 + t12 / 0.2e1 - t5 * mrSges(6,3) + t3 * mrSges(7,2)) * t61 + (t9 / 0.2e1 - t10 / 0.2e1 - t2 * mrSges(7,2) - t6 * mrSges(6,3)) * t60; (-mrSges(4,2) + t110) * t92 + (-t22 - t23 - t116) * t95 + m(7) * (-t19 * t95 + t109) + m(6) * (-t80 * t95 + t109) + m(5) * (t115 * t126 + t127) + t130 * (t48 * t61 + t50 * t60); -0.2e1 * pkin(3) * t66 + 0.2e1 * t19 * t22 + 0.2e1 * t80 * t23 + t94 * t67 + t91 * t68 + Ifges(4,3) + m(7) * (t19 ^ 2 + t113) + m(6) * (t80 ^ 2 + t113) + m(5) * (t115 * pkin(8) ^ 2 + pkin(3) ^ 2) + (t26 + t27) * t61 + (t24 - t25) * t60 + 0.2e1 * pkin(8) * t110 + 0.2e1 * (t31 * t61 - t33 * t60) * t130; t100 + m(7) * (t2 * t76 + t3 * t79) + t112 * t95 + (t93 * t36 + t90 * t35 + m(6) * (t5 * t93 + t6 * t90)) * pkin(4) + t72 + t79 * t37 + t76 * t34 - t21 * mrSges(5,2) + t20 * mrSges(5,1) - Ifges(5,6) * t122; m(7) * (t48 * t79 - t50 * t76) + m(6) * (-t48 * t93 - t50 * t90) * pkin(4) + t132; m(7) * (t31 * t79 + t33 * t76) + t83 + t82 - t106 * pkin(8) + (-t60 * t76 + t61 * t79) * mrSges(7,2) + (m(6) * (-t31 * t93 + t33 * t90) + (-t60 * t90 - t61 * t93) * mrSges(6,3)) * pkin(4) + t99; -0.2e1 * t79 * mrSges(7,1) + t76 * t129 + 0.2e1 * t102 + m(7) * (t76 ^ 2 + t79 ^ 2) + m(6) * (t90 ^ 2 + t93 ^ 2) * pkin(4) ^ 2 - t112; -pkin(5) * t37 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) + qJ(6) * t34 + t119 * t95 + t100; m(7) * t105 + t101; m(7) * (-pkin(5) * t31 + qJ(6) * t33) + (-pkin(5) * t61 - qJ(6) * t60) * mrSges(7,2) + t99; m(7) * (-pkin(5) * t79 + qJ(6) * t76) + t102 + (t76 + qJ(6)) * mrSges(7,3) + (-t79 + pkin(5)) * mrSges(7,1) - t119; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t129 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t119; m(7) * t3 + t37; m(7) * t48; m(7) * t31 + t61 * mrSges(7,2); m(7) * t79 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
