% Calculate joint inertia matrix for
% S6RRPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
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
% Datum: 2019-03-09 09:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR7_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:17:12
% EndTime: 2019-03-09 09:17:14
% DurationCPUTime: 1.13s
% Computational Cost: add. (1392->293), mult. (2877->400), div. (0->0), fcn. (2696->8), ass. (0->112)
t147 = Ifges(4,4) + Ifges(3,5);
t146 = Ifges(4,2) + Ifges(3,3) + Ifges(5,3);
t101 = cos(qJ(6));
t100 = sin(qJ(2));
t95 = sin(pkin(6));
t121 = t95 * t100;
t102 = cos(qJ(5));
t103 = cos(qJ(2));
t122 = t103 * t95;
t96 = cos(pkin(6));
t99 = sin(qJ(5));
t43 = -t102 * t122 - t96 * t99;
t98 = sin(qJ(6));
t25 = t101 * t121 - t43 * t98;
t42 = -t96 * t102 + t99 * t122;
t12 = mrSges(7,2) * t42 + mrSges(7,3) * t25;
t26 = t101 * t43 + t98 * t121;
t13 = -mrSges(7,1) * t42 - mrSges(7,3) * t26;
t145 = t101 * t12 - t98 * t13;
t55 = -mrSges(7,1) * t101 + mrSges(7,2) * t98;
t144 = -m(7) * pkin(5) - mrSges(6,1) + t55;
t34 = -t95 * pkin(1) - pkin(2) * t122 - qJ(3) * t121;
t143 = -0.2e1 * t34;
t142 = t25 / 0.2e1;
t141 = t26 / 0.2e1;
t140 = t98 / 0.2e1;
t104 = -pkin(2) - pkin(3);
t28 = pkin(3) * t122 - t34;
t17 = (pkin(4) * t100 + pkin(9) * t103) * t95 + t28;
t69 = pkin(8) * t121;
t114 = -qJ(4) * t121 + t69;
t134 = pkin(1) * t103;
t119 = -pkin(2) - t134;
t115 = -pkin(3) + t119;
t19 = (-pkin(9) + t115) * t96 + t114;
t5 = t102 * t17 - t19 * t99;
t3 = -pkin(5) * t121 - t5;
t139 = t3 * t99;
t138 = -t101 / 0.2e1;
t137 = t101 / 0.2e1;
t92 = t99 ^ 2;
t94 = t102 ^ 2;
t127 = t94 + t92;
t136 = m(6) * t127 + m(5);
t135 = Ifges(7,4) * t98;
t44 = t96 * t134 - t69;
t133 = t44 * mrSges(3,1);
t45 = t96 * t100 * pkin(1) + pkin(8) * t122;
t132 = t45 * mrSges(3,2);
t130 = t98 * t99;
t97 = qJ(3) + pkin(4);
t11 = -mrSges(7,1) * t25 + mrSges(7,2) * t26;
t30 = mrSges(6,1) * t121 - mrSges(6,3) * t43;
t129 = t11 - t30;
t6 = t102 * t19 + t99 * t17;
t57 = Ifges(7,5) * t98 + Ifges(7,6) * t101;
t128 = t101 ^ 2 + t98 ^ 2;
t126 = Ifges(7,4) * t101;
t124 = t101 * t99;
t90 = -pkin(9) + t104;
t123 = t102 * t90;
t7 = Ifges(7,5) * t26 + Ifges(7,6) * t25 - Ifges(7,3) * t42;
t120 = Ifges(6,5) * t43 + Ifges(6,6) * t42 + Ifges(6,3) * t121;
t33 = t96 * qJ(3) + t45;
t118 = t128 * t99;
t117 = t127 * mrSges(6,3);
t48 = -t96 * mrSges(4,1) + mrSges(4,2) * t121;
t56 = t102 * mrSges(6,1) - t99 * mrSges(6,2);
t47 = t96 * mrSges(5,2) - mrSges(5,3) * t121;
t27 = qJ(4) * t122 - t33;
t21 = t96 * pkin(4) - t27;
t10 = -pkin(5) * t42 - pkin(10) * t43 + t21;
t4 = pkin(10) * t121 + t6;
t1 = t10 * t101 - t4 * t98;
t2 = t10 * t98 + t101 * t4;
t113 = -t1 * t98 + t101 * t2;
t112 = t102 * t6 - t5 * t99;
t111 = -mrSges(7,1) * t98 - mrSges(7,2) * t101;
t88 = t102 * pkin(5);
t54 = pkin(10) * t99 + t88 + t97;
t31 = t101 * t54 - t98 * t123;
t32 = t101 * t123 + t54 * t98;
t110 = t101 * t32 - t31 * t98;
t52 = -mrSges(7,2) * t102 + mrSges(7,3) * t130;
t53 = mrSges(7,1) * t102 + mrSges(7,3) * t124;
t109 = t101 * t52 - t98 * t53;
t36 = -Ifges(7,5) * t124 + Ifges(7,6) * t130 + Ifges(7,3) * t102;
t29 = -mrSges(6,2) * t121 + mrSges(6,3) * t42;
t108 = t29 + t145;
t107 = Ifges(3,6) * t122 + t147 * t121 + t146 * t96;
t105 = qJ(3) ^ 2;
t89 = t90 ^ 2;
t76 = t92 * t90;
t75 = t92 * t89;
t64 = mrSges(5,1) * t121;
t61 = -Ifges(6,1) * t99 - Ifges(6,4) * t102;
t60 = Ifges(7,1) * t98 + t126;
t59 = -Ifges(6,4) * t99 - Ifges(6,2) * t102;
t58 = Ifges(7,2) * t101 + t135;
t50 = mrSges(4,2) * t122 + mrSges(4,3) * t96;
t49 = -mrSges(5,1) * t96 + mrSges(5,3) * t122;
t46 = t111 * t99;
t38 = Ifges(7,5) * t102 + (-Ifges(7,1) * t101 + t135) * t99;
t37 = Ifges(7,6) * t102 + (Ifges(7,2) * t98 - t126) * t99;
t35 = t119 * t96 + t69;
t22 = t115 * t96 + t114;
t20 = -mrSges(6,1) * t42 + mrSges(6,2) * t43;
t16 = Ifges(6,1) * t43 + Ifges(6,4) * t42 + Ifges(6,5) * t121;
t15 = Ifges(6,4) * t43 + Ifges(6,2) * t42 + Ifges(6,6) * t121;
t9 = Ifges(7,1) * t26 + Ifges(7,4) * t25 - Ifges(7,5) * t42;
t8 = Ifges(7,4) * t26 + Ifges(7,2) * t25 - Ifges(7,6) * t42;
t14 = [0.2e1 * t1 * t13 + 0.2e1 * t3 * t11 + 0.2e1 * t2 * t12 + t43 * t16 + 0.2e1 * t21 * t20 + 0.2e1 * t22 * t47 + t25 * t8 + t26 * t9 + 0.2e1 * t27 * t49 + 0.2e1 * t28 * t64 + 0.2e1 * t6 * t29 + 0.2e1 * t5 * t30 + 0.2e1 * t33 * t50 + 0.2e1 * t35 * t48 + Ifges(2,3) + (t15 - t7) * t42 + (t107 - 0.2e1 * t132 + 0.2e1 * t133) * t96 + m(3) * (t44 ^ 2 + t45 ^ 2) + m(4) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) + m(6) * (t21 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t22 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + (m(3) * pkin(1) ^ 2 * t95 + (mrSges(4,1) * t143 - 0.2e1 * t28 * mrSges(5,2) + 0.2e1 * t45 * mrSges(3,3) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(5,1) + Ifges(4,3) + Ifges(3,2)) * t103) * t95 + ((2 * Ifges(5,5)) + Ifges(3,6) - (2 * Ifges(4,6))) * t96) * t103 + (-0.2e1 * t44 * mrSges(3,3) + mrSges(4,3) * t143 + (-0.2e1 * pkin(1) * mrSges(3,2) + (Ifges(5,2) + Ifges(4,1) + Ifges(3,1)) * t100) * t95 + 0.2e1 * (Ifges(3,4) + Ifges(5,4) - Ifges(4,5)) * t122 + ((2 * Ifges(5,6)) + t147) * t96 + t120) * t100) * t95; m(6) * (t112 * t90 + t21 * t97) + t107 + (t9 * t138 + t5 * mrSges(6,3) + t8 * t140 - t16 / 0.2e1 + t129 * t90) * t99 + m(7) * (t1 * t31 + t90 * t139 + t2 * t32) + t38 * t141 + t37 * t142 - t132 + t133 + t104 * t47 + t97 * t20 + t2 * t52 + t1 * t53 + t21 * t56 + t43 * t61 / 0.2e1 + t3 * t46 - pkin(2) * t48 + t31 * t13 + t32 * t12 + t33 * mrSges(4,3) - t35 * mrSges(4,1) + t22 * mrSges(5,2) - t27 * mrSges(5,1) + ((-Ifges(6,5) * t99 / 0.2e1 - Ifges(6,6) * t102 / 0.2e1 + Ifges(5,6)) * t100 + (-Ifges(4,6) + Ifges(5,5)) * t103) * t95 + (t59 / 0.2e1 - t36 / 0.2e1) * t42 + m(5) * (-qJ(3) * t27 + t104 * t22) + m(4) * (-pkin(2) * t35 + qJ(3) * t33) + (t50 - t49) * qJ(3) + (t90 * t29 - t6 * mrSges(6,3) - t15 / 0.2e1 + t7 / 0.2e1) * t102; 0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * t104 * mrSges(5,2) + 0.2e1 * t31 * t53 + 0.2e1 * t32 * t52 + 0.2e1 * t97 * t56 + (t36 - t59) * t102 + (-t101 * t38 + t37 * t98 + 0.2e1 * t46 * t90 - t61) * t99 + m(7) * (t31 ^ 2 + t32 ^ 2 + t75) + m(6) * (t89 * t94 + t97 ^ 2 + t75) + m(4) * (pkin(2) ^ 2 + t105) + m(5) * (t104 ^ 2 + t105) + 0.2e1 * (mrSges(5,1) + mrSges(4,3)) * qJ(3) - 0.2e1 * t90 * t117 + t146; t129 * t99 + t108 * t102 + m(7) * (t113 * t102 + t139) + m(6) * t112 + m(5) * t22 + m(4) * t35 + t47 + t48; -m(4) * pkin(2) + t99 * t46 - mrSges(4,1) + mrSges(5,2) + t109 * t102 - t117 + m(7) * (t110 * t102 + t76) + m(6) * (t90 * t94 + t76) + m(5) * t104; m(4) + m(7) * (t128 * t94 + t92) + t136; -mrSges(5,2) * t122 + t64 - t129 * t102 + t108 * t99 + m(7) * (-t102 * t3 + t113 * t99) + m(6) * (t102 * t5 + t6 * t99) + m(5) * t28; -t102 * t46 + (m(7) * (t110 - t123) + t109) * t99; m(7) * (-0.1e1 + t128) * t99 * t102; m(7) * (t128 * t92 + t94) + t136; t58 * t142 - t42 * t57 / 0.2e1 + t60 * t141 + t3 * t55 + t9 * t140 + t8 * t137 - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t113 * mrSges(7,3) + t120 + (-m(7) * t3 - t11) * pkin(5) + (m(7) * t113 + t145) * pkin(10); t37 * t137 - pkin(5) * t46 + t38 * t140 + t110 * mrSges(7,3) + (m(7) * t110 + t109) * pkin(10) + (t60 * t138 + t58 * t140 + t144 * t90 - Ifges(6,5)) * t99 + (t57 / 0.2e1 - Ifges(6,6) - t90 * mrSges(6,2)) * t102; t144 * t99 + (-mrSges(6,2) + (m(7) * pkin(10) + mrSges(7,3)) * t128) * t102; -t102 * t55 + m(7) * (pkin(10) * t118 + t88) + mrSges(7,3) * t118 + t56; Ifges(6,3) + m(7) * (t128 * pkin(10) ^ 2 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t55 + t98 * t60 + t101 * t58 + 0.2e1 * t128 * pkin(10) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t31 - mrSges(7,2) * t32 + t36; t111 * t102; t46; t111 * pkin(10) + t57; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t14(1) t14(2) t14(4) t14(7) t14(11) t14(16); t14(2) t14(3) t14(5) t14(8) t14(12) t14(17); t14(4) t14(5) t14(6) t14(9) t14(13) t14(18); t14(7) t14(8) t14(9) t14(10) t14(14) t14(19); t14(11) t14(12) t14(13) t14(14) t14(15) t14(20); t14(16) t14(17) t14(18) t14(19) t14(20) t14(21);];
Mq  = res;
