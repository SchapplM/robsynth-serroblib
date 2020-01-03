% Calculate joint inertia matrix for
% S5RRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR14_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR14_inertiaJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR14_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR14_inertiaJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR14_inertiaJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR14_inertiaJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:35:49
% EndTime: 2019-12-31 20:35:51
% DurationCPUTime: 0.90s
% Computational Cost: add. (1717->257), mult. (3877->380), div. (0->0), fcn. (4231->10), ass. (0->105)
t100 = sin(pkin(5));
t107 = cos(qJ(2));
t116 = t100 * t107;
t104 = sin(qJ(4));
t126 = pkin(8) + qJ(3);
t99 = sin(pkin(10));
t114 = t126 * t99;
t132 = cos(qJ(4));
t101 = cos(pkin(10));
t77 = t126 * t101;
t51 = t104 * t77 + t114 * t132;
t141 = t51 ^ 2;
t140 = 0.2e1 * t51;
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t102 = cos(pkin(5));
t105 = sin(qJ(2));
t117 = t100 * t105;
t65 = t101 * t102 - t117 * t99;
t66 = t101 * t117 + t102 * t99;
t41 = t104 * t65 + t132 * t66;
t25 = -t103 * t41 - t106 * t116;
t139 = t25 / 0.2e1;
t26 = -t103 * t116 + t106 * t41;
t138 = t26 / 0.2e1;
t81 = Ifges(6,5) * t103 + Ifges(6,6) * t106;
t137 = t81 / 0.2e1;
t136 = -t103 / 0.2e1;
t135 = t103 / 0.2e1;
t134 = t106 / 0.2e1;
t131 = pkin(1) * t107;
t130 = pkin(9) * t103;
t129 = pkin(9) * t106;
t86 = pkin(7) * t117;
t67 = t102 * t131 - t86;
t128 = t67 * mrSges(3,1);
t68 = t102 * t105 * pkin(1) + pkin(7) * t116;
t127 = t68 * mrSges(3,2);
t59 = qJ(3) * t102 + t68;
t60 = (-pkin(2) * t107 - qJ(3) * t105 - pkin(1)) * t100;
t32 = t101 * t60 - t59 * t99;
t18 = -pkin(3) * t116 - pkin(8) * t66 + t32;
t33 = t101 * t59 + t99 * t60;
t22 = pkin(8) * t65 + t33;
t9 = t104 * t18 + t132 * t22;
t40 = t104 * t66 - t132 * t65;
t125 = -Ifges(5,5) * t41 + Ifges(5,6) * t40;
t74 = -t101 * t132 + t104 * t99;
t75 = t104 * t101 + t132 * t99;
t124 = Ifges(5,5) * t75 - Ifges(5,6) * t74;
t123 = t101 ^ 2 + t99 ^ 2;
t122 = t103 ^ 2 + t106 ^ 2;
t121 = Ifges(6,4) * t103;
t120 = Ifges(6,4) * t106;
t119 = t103 * t75;
t118 = t106 * t75;
t5 = Ifges(6,5) * t26 + Ifges(6,6) * t25 + Ifges(6,3) * t40;
t115 = Ifges(3,5) * t117 + Ifges(3,6) * t116 + Ifges(3,3) * t102;
t90 = -pkin(3) * t101 - pkin(2);
t43 = -t65 * mrSges(4,1) + t66 * mrSges(4,2);
t16 = t40 * mrSges(5,1) + t41 * mrSges(5,2);
t48 = t74 * mrSges(5,1) + t75 * mrSges(5,2);
t76 = -t101 * mrSges(4,1) + t99 * mrSges(4,2);
t62 = t86 + (-pkin(2) - t131) * t102;
t42 = -pkin(3) * t65 + t62;
t10 = pkin(4) * t40 - pkin(9) * t41 + t42;
t4 = -pkin(9) * t116 + t9;
t1 = t10 * t106 - t103 * t4;
t2 = t10 * t103 + t106 * t4;
t113 = -t1 * t103 + t106 * t2;
t112 = t101 * t33 - t32 * t99;
t111 = mrSges(6,1) * t103 + mrSges(6,2) * t106;
t27 = Ifges(6,5) * t118 - Ifges(6,6) * t119 + Ifges(6,3) * t74;
t8 = -t104 * t22 + t132 * t18;
t83 = Ifges(6,1) * t103 + t120;
t82 = Ifges(6,2) * t106 + t121;
t80 = -mrSges(6,1) * t106 + mrSges(6,2) * t103;
t79 = Ifges(4,1) * t99 + Ifges(4,4) * t101;
t78 = Ifges(4,4) * t99 + Ifges(4,2) * t101;
t55 = -mrSges(4,1) * t116 - mrSges(4,3) * t66;
t54 = mrSges(4,2) * t116 + mrSges(4,3) * t65;
t53 = -t104 * t114 + t132 * t77;
t50 = Ifges(5,1) * t75 - Ifges(5,4) * t74;
t49 = Ifges(5,4) * t75 - Ifges(5,2) * t74;
t47 = mrSges(6,1) * t74 - mrSges(6,3) * t118;
t46 = -mrSges(6,2) * t74 - mrSges(6,3) * t119;
t45 = pkin(4) * t74 - pkin(9) * t75 + t90;
t44 = t111 * t75;
t35 = Ifges(4,1) * t66 + Ifges(4,4) * t65 - Ifges(4,5) * t116;
t34 = Ifges(4,4) * t66 + Ifges(4,2) * t65 - Ifges(4,6) * t116;
t31 = -mrSges(5,1) * t116 - mrSges(5,3) * t41;
t30 = mrSges(5,2) * t116 - mrSges(5,3) * t40;
t29 = Ifges(6,5) * t74 + (Ifges(6,1) * t106 - t121) * t75;
t28 = Ifges(6,6) * t74 + (-Ifges(6,2) * t103 + t120) * t75;
t20 = t103 * t45 + t106 * t53;
t19 = -t103 * t53 + t106 * t45;
t15 = Ifges(5,1) * t41 - Ifges(5,4) * t40 - Ifges(5,5) * t116;
t14 = Ifges(5,4) * t41 - Ifges(5,2) * t40 - Ifges(5,6) * t116;
t13 = mrSges(6,1) * t40 - mrSges(6,3) * t26;
t12 = -mrSges(6,2) * t40 + mrSges(6,3) * t25;
t11 = -mrSges(6,1) * t25 + mrSges(6,2) * t26;
t7 = Ifges(6,1) * t26 + Ifges(6,4) * t25 + Ifges(6,5) * t40;
t6 = Ifges(6,4) * t26 + Ifges(6,2) * t25 + Ifges(6,6) * t40;
t3 = pkin(4) * t116 - t8;
t17 = [0.2e1 * t1 * t13 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + t41 * t15 + 0.2e1 * t42 * t16 + t25 * t6 + t26 * t7 + 0.2e1 * t9 * t30 + 0.2e1 * t8 * t31 + 0.2e1 * t32 * t55 + 0.2e1 * t33 * t54 + t65 * t34 + t66 * t35 + 0.2e1 * t62 * t43 + Ifges(2,3) + (t5 - t14) * t40 + (t115 - 0.2e1 * t127 + 0.2e1 * t128) * t102 + ((-0.2e1 * t67 * mrSges(3,3) + Ifges(3,5) * t102 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t105) * t100) * t105 + (0.2e1 * t68 * mrSges(3,3) - Ifges(4,5) * t66 + Ifges(3,6) * t102 - Ifges(4,6) * t65 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t105 + (Ifges(3,2) + Ifges(4,3) + Ifges(5,3)) * t107) * t100 + t125) * t107) * t100 + m(3) * (pkin(1) ^ 2 * t100 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(4) * (t32 ^ 2 + t33 ^ 2 + t62 ^ 2) + m(5) * (t42 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2); m(5) * (t42 * t90 - t51 * t8 + t53 * t9) + m(6) * (t1 * t19 + t2 * t20 + t3 * t51) + (t5 / 0.2e1 - t14 / 0.2e1 - t9 * mrSges(5,3)) * t74 + t115 + (-t49 / 0.2e1 + t27 / 0.2e1) * t40 + t29 * t138 + t28 * t139 + t112 * mrSges(4,3) + m(4) * (-pkin(2) * t62 + qJ(3) * t112) - (Ifges(4,5) * t99 + Ifges(4,6) * t101 + t124) * t116 / 0.2e1 + t99 * t35 / 0.2e1 + t101 * t34 / 0.2e1 + t65 * t78 / 0.2e1 + t66 * t79 / 0.2e1 + t90 * t16 + t62 * t76 - pkin(2) * t43 + t3 * t44 + t2 * t46 + t1 * t47 + t42 * t48 + t41 * t50 / 0.2e1 + t53 * t30 + t19 * t13 + t20 * t12 + (t15 / 0.2e1 + t6 * t136 + t7 * t134 - t8 * mrSges(5,3)) * t75 + t128 - t127 + (t101 * t54 - t99 * t55) * qJ(3) + (t11 - t31) * t51; -0.2e1 * pkin(2) * t76 + t101 * t78 + 0.2e1 * t19 * t47 + 0.2e1 * t20 * t46 + t44 * t140 + 0.2e1 * t90 * t48 + t99 * t79 + Ifges(3,3) + 0.2e1 * t123 * qJ(3) * mrSges(4,3) + (-0.2e1 * mrSges(5,3) * t53 + t27 - t49) * t74 + m(6) * (t19 ^ 2 + t20 ^ 2 + t141) + m(5) * (t53 ^ 2 + t90 ^ 2 + t141) + m(4) * (qJ(3) ^ 2 * t123 + pkin(2) ^ 2) + (mrSges(5,3) * t140 - t103 * t28 + t106 * t29 + t50) * t75; t103 * t12 + t106 * t13 + m(6) * (t1 * t106 + t103 * t2) + m(5) * t42 + m(4) * t62 + t16 + t43; -m(4) * pkin(2) + t103 * t46 + t106 * t47 + m(6) * (t103 * t20 + t106 * t19) + m(5) * t90 + t76 + t48; m(6) * t122 + m(4) + m(5); -Ifges(5,3) * t116 + m(6) * (-pkin(4) * t3 + pkin(9) * t113) + t12 * t129 - t13 * t130 + t3 * t80 + t6 * t134 + t7 * t135 - pkin(4) * t11 + t83 * t138 + t82 * t139 + t40 * t137 - t9 * mrSges(5,2) + t8 * mrSges(5,1) + t113 * mrSges(6,3) - t125; t74 * t137 + t29 * t135 + t28 * t134 - pkin(4) * t44 + t46 * t129 - t47 * t130 - t53 * mrSges(5,2) + (t134 * t83 + t136 * t82) * t75 + t124 + (m(6) * pkin(9) + mrSges(6,3)) * (-t103 * t19 + t106 * t20) + (-m(6) * pkin(4) - mrSges(5,1) + t80) * t51; 0; Ifges(5,3) + t103 * t83 + t106 * t82 + m(6) * (pkin(9) ^ 2 * t122 + pkin(4) ^ 2) - 0.2e1 * pkin(4) * t80 + 0.2e1 * t122 * pkin(9) * mrSges(6,3); mrSges(6,1) * t1 - mrSges(6,2) * t2 + t5; mrSges(6,1) * t19 - mrSges(6,2) * t20 + t27; -t80; -pkin(9) * t111 + t81; Ifges(6,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
