% Calculate joint inertia matrix for
% S6RRPRPP3
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPP3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_inertiaJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:16
% EndTime: 2019-03-09 09:54:18
% DurationCPUTime: 0.93s
% Computational Cost: add. (1228->289), mult. (2481->366), div. (0->0), fcn. (2378->6), ass. (0->104)
t133 = 2 * pkin(7);
t132 = m(6) + m(7);
t101 = sin(pkin(9));
t131 = -t101 / 0.2e1;
t102 = cos(pkin(9));
t130 = t102 / 0.2e1;
t129 = cos(qJ(4));
t106 = cos(qJ(2));
t128 = pkin(7) * t106;
t105 = sin(qJ(2));
t95 = t105 * pkin(7);
t127 = mrSges(5,1) - mrSges(6,2);
t125 = pkin(4) + qJ(6);
t124 = pkin(8) + qJ(3);
t104 = sin(qJ(4));
t118 = t102 * t105;
t80 = -pkin(2) * t106 - qJ(3) * t105 - pkin(1);
t73 = t102 * t80;
t23 = -pkin(8) * t118 + t73 + (-pkin(7) * t101 - pkin(3)) * t106;
t119 = t101 * t105;
t47 = t101 * t80 + t102 * t128;
t35 = -pkin(8) * t119 + t47;
t7 = t104 * t23 + t129 * t35;
t62 = mrSges(4,1) * t119 + mrSges(4,2) * t118;
t115 = t129 * t102;
t61 = -t104 * t119 + t105 * t115;
t40 = t61 * mrSges(7,1) + t106 * mrSges(7,3);
t79 = pkin(3) * t119 + t95;
t123 = t101 ^ 2 + t102 ^ 2;
t122 = Ifges(4,4) * t101;
t121 = Ifges(4,4) * t102;
t75 = t101 * t104 - t115;
t120 = qJ(5) * t75;
t117 = -Ifges(5,3) - Ifges(6,1) - Ifges(7,1);
t114 = t124 * t101;
t82 = t124 * t102;
t36 = t104 * t82 + t129 * t114;
t38 = -t104 * t114 + t129 * t82;
t116 = t36 ^ 2 + t38 ^ 2;
t91 = -pkin(3) * t102 - pkin(2);
t76 = t129 * t101 + t104 * t102;
t60 = t76 * t105;
t18 = -t61 * mrSges(7,2) + t60 * mrSges(7,3);
t26 = -t76 * mrSges(7,2) + t75 * mrSges(7,3);
t81 = -t102 * mrSges(4,1) + t101 * mrSges(4,2);
t113 = -qJ(5) * t61 + t79;
t3 = qJ(5) * t106 - t7;
t6 = -t104 * t35 + t129 * t23;
t43 = t61 * mrSges(6,1) - t106 * mrSges(6,2);
t42 = -t60 * mrSges(7,1) - t106 * mrSges(7,2);
t46 = -t101 * t128 + t73;
t112 = -t46 * t101 + t47 * t102;
t111 = (Ifges(6,4) - Ifges(5,5) - Ifges(7,5)) * t61 + (-Ifges(7,4) - Ifges(6,5) + Ifges(5,6)) * t60;
t4 = t106 * pkin(4) - t6;
t110 = -qJ(5) * t76 + t91;
t109 = pkin(7) ^ 2;
t107 = qJ(5) ^ 2;
t100 = t106 ^ 2;
t99 = t105 ^ 2;
t94 = t99 * t109;
t84 = Ifges(4,1) * t101 + t121;
t83 = Ifges(4,2) * t102 + t122;
t78 = -mrSges(4,1) * t106 - mrSges(4,3) * t118;
t77 = mrSges(4,2) * t106 - mrSges(4,3) * t119;
t71 = Ifges(6,4) * t76;
t70 = Ifges(7,4) * t75;
t69 = Ifges(5,5) * t76;
t68 = Ifges(6,5) * t75;
t67 = Ifges(7,5) * t76;
t66 = Ifges(5,6) * t75;
t65 = t76 * mrSges(6,3);
t64 = t76 * mrSges(5,2);
t59 = -Ifges(4,5) * t106 + (Ifges(4,1) * t102 - t122) * t105;
t58 = -Ifges(4,6) * t106 + (-Ifges(4,2) * t101 + t121) * t105;
t51 = t61 * mrSges(6,3);
t49 = t61 * mrSges(5,2);
t45 = -mrSges(5,1) * t106 - mrSges(5,3) * t61;
t44 = mrSges(5,2) * t106 - mrSges(5,3) * t60;
t41 = mrSges(6,1) * t60 + mrSges(6,3) * t106;
t34 = Ifges(5,1) * t76 - Ifges(5,4) * t75;
t33 = Ifges(5,4) * t76 - Ifges(5,2) * t75;
t32 = -Ifges(6,2) * t76 + Ifges(6,6) * t75;
t31 = Ifges(7,2) * t75 + Ifges(7,6) * t76;
t30 = -Ifges(6,6) * t76 + Ifges(6,3) * t75;
t29 = Ifges(7,6) * t75 + Ifges(7,3) * t76;
t28 = -t75 * mrSges(6,2) - t65;
t27 = t75 * mrSges(5,1) + t64;
t22 = pkin(4) * t75 + t110;
t20 = -t60 * mrSges(6,2) - t51;
t19 = t60 * mrSges(5,1) + t49;
t17 = -t75 * pkin(5) + t38;
t16 = pkin(5) * t76 + t36;
t15 = Ifges(5,1) * t61 - Ifges(5,4) * t60 - Ifges(5,5) * t106;
t14 = Ifges(5,4) * t61 - Ifges(5,2) * t60 - Ifges(5,6) * t106;
t13 = -Ifges(6,4) * t106 - Ifges(6,2) * t61 + Ifges(6,6) * t60;
t12 = -Ifges(7,4) * t106 + Ifges(7,2) * t60 + Ifges(7,6) * t61;
t11 = -Ifges(6,5) * t106 - Ifges(6,6) * t61 + Ifges(6,3) * t60;
t10 = -Ifges(7,5) * t106 + Ifges(7,6) * t60 + Ifges(7,3) * t61;
t9 = t125 * t75 + t110;
t8 = pkin(4) * t60 + t113;
t5 = t125 * t60 + t113;
t2 = -pkin(5) * t60 - t3;
t1 = t61 * pkin(5) + t106 * qJ(6) + t4;
t21 = [0.2e1 * t1 * t40 + 0.2e1 * t5 * t18 + 0.2e1 * t79 * t19 + 0.2e1 * t2 * t42 + 0.2e1 * t8 * t20 + 0.2e1 * t3 * t41 + 0.2e1 * t4 * t43 + 0.2e1 * t7 * t44 + 0.2e1 * t6 * t45 + 0.2e1 * t46 * t78 + 0.2e1 * t47 * t77 + Ifges(2,3) + (t100 + t99) * mrSges(3,3) * t133 + (t15 + t10 - t13) * t61 + (t11 + t12 - t14) * t60 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t105 - t101 * t58 + t102 * t59 + t133 * t62) * t105 + m(3) * (pkin(1) ^ 2 + t100 * t109 + t94) + m(5) * (t6 ^ 2 + t7 ^ 2 + t79 ^ 2) + m(4) * (t46 ^ 2 + t47 ^ 2 + t94) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + m(6) * (t3 ^ 2 + t4 ^ 2 + t8 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(4,3) + Ifges(3,2) - t117) * t106 + (-Ifges(4,5) * t102 + Ifges(4,6) * t101 + (2 * Ifges(3,4))) * t105 + t111) * t106; (-t41 + t44) * t38 + (-t45 + t43) * t36 + (-t101 * t78 + t102 * t77) * qJ(3) + m(4) * (-pkin(2) * t95 + t112 * qJ(3)) + t101 * t59 / 0.2e1 + t91 * t19 + t79 * t27 - pkin(2) * t62 + t8 * t28 + t16 * t40 + t17 * t42 + t9 * t18 + t22 * t20 + t5 * t26 + (t83 * t131 + t84 * t130 + Ifges(3,5) + (t81 - mrSges(3,1)) * pkin(7)) * t105 + (-pkin(7) * mrSges(3,2) + Ifges(4,5) * t131 - Ifges(4,6) * t102 / 0.2e1 - t69 / 0.2e1 + t66 / 0.2e1 - t70 / 0.2e1 - t67 / 0.2e1 + t71 / 0.2e1 - t68 / 0.2e1 + Ifges(3,6)) * t106 + t112 * mrSges(4,3) + (t4 * mrSges(6,1) - t6 * mrSges(5,3) + t1 * mrSges(7,1) + t15 / 0.2e1 + t10 / 0.2e1 - t13 / 0.2e1) * t76 + (t29 / 0.2e1 + t34 / 0.2e1 - t32 / 0.2e1) * t61 + (-t2 * mrSges(7,1) + t3 * mrSges(6,1) - t7 * mrSges(5,3) + t11 / 0.2e1 + t12 / 0.2e1 - t14 / 0.2e1) * t75 + t58 * t130 + m(5) * (-t36 * t6 + t38 * t7 + t79 * t91) + m(6) * (t22 * t8 - t3 * t38 + t36 * t4) + m(7) * (t1 * t16 + t17 * t2 + t5 * t9) + (t30 / 0.2e1 + t31 / 0.2e1 - t33 / 0.2e1) * t60; -0.2e1 * pkin(2) * t81 + t101 * t84 + t102 * t83 + 0.2e1 * t22 * t28 + 0.2e1 * t9 * t26 + 0.2e1 * t91 * t27 + Ifges(3,3) + m(4) * (t123 * qJ(3) ^ 2 + pkin(2) ^ 2) + m(5) * (t91 ^ 2 + t116) + m(6) * (t22 ^ 2 + t116) + m(7) * (t16 ^ 2 + t17 ^ 2 + t9 ^ 2) + (0.2e1 * t16 * mrSges(7,1) + t29 - t32 + t34) * t76 + (-0.2e1 * t17 * mrSges(7,1) + t30 + t31 - t33) * t75 + 0.2e1 * t123 * qJ(3) * mrSges(4,3) + 0.2e1 * (t36 * t76 - t38 * t75) * (mrSges(6,1) + mrSges(5,3)); m(4) * t95 + m(5) * t79 + m(6) * t8 + m(7) * t5 + t127 * t60 + t18 + t49 - t51 + t62; -m(4) * pkin(2) + m(5) * t91 + m(6) * t22 + m(7) * t9 + t127 * t75 + t26 + t64 - t65 + t81; m(4) + m(5) + t132; (-t41 + t42) * qJ(5) + m(7) * (qJ(5) * t2 - t1 * t125) + m(6) * (-pkin(4) * t4 - qJ(5) * t3) + t117 * t106 - t111 - t125 * t40 - pkin(4) * t43 - t7 * mrSges(5,2) - t1 * mrSges(7,3) + t2 * mrSges(7,2) - t3 * mrSges(6,3) + t4 * mrSges(6,2) + t6 * mrSges(5,1); t17 * mrSges(7,2) - t16 * mrSges(7,3) - t66 + t67 + t68 + t69 + t70 - t71 + (-mrSges(5,2) + mrSges(6,3)) * t38 - t127 * t36 + (-t125 * t76 - t120) * mrSges(7,1) + (-pkin(4) * t76 - t120) * mrSges(6,1) + m(6) * (-pkin(4) * t36 + qJ(5) * t38) + m(7) * (qJ(5) * t17 - t125 * t16); 0; -0.2e1 * pkin(4) * mrSges(6,2) + 0.2e1 * t125 * mrSges(7,3) + 0.2e1 * (mrSges(6,3) + mrSges(7,2)) * qJ(5) + m(6) * (pkin(4) ^ 2 + t107) + m(7) * (t125 ^ 2 + t107) - t117; m(6) * t4 + m(7) * t1 + t40 + t43; (mrSges(7,1) + mrSges(6,1)) * t76 + m(6) * t36 + m(7) * t16; 0; -m(6) * pkin(4) - m(7) * t125 + mrSges(6,2) - mrSges(7,3); t132; m(7) * t2 + t42; m(7) * t17 - t75 * mrSges(7,1); 0; m(7) * qJ(5) + mrSges(7,2); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t21(1) t21(2) t21(4) t21(7) t21(11) t21(16); t21(2) t21(3) t21(5) t21(8) t21(12) t21(17); t21(4) t21(5) t21(6) t21(9) t21(13) t21(18); t21(7) t21(8) t21(9) t21(10) t21(14) t21(19); t21(11) t21(12) t21(13) t21(14) t21(15) t21(20); t21(16) t21(17) t21(18) t21(19) t21(20) t21(21);];
Mq  = res;
