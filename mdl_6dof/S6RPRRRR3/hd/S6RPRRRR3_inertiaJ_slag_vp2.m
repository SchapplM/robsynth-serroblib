% Calculate joint inertia matrix for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:24
% EndTime: 2019-03-09 07:00:27
% DurationCPUTime: 1.21s
% Computational Cost: add. (2052->299), mult. (4015->433), div. (0->0), fcn. (4094->10), ass. (0->117)
t102 = sin(qJ(3));
t100 = sin(qJ(5));
t101 = sin(qJ(4));
t104 = cos(qJ(5));
t105 = cos(qJ(4));
t69 = t100 * t105 + t101 * t104;
t57 = t69 * t102;
t68 = -t100 * t101 + t104 * t105;
t58 = t68 * t102;
t31 = t57 * mrSges(6,1) + t58 * mrSges(6,2);
t103 = cos(qJ(6));
t99 = sin(qJ(6));
t27 = -t103 * t57 - t58 * t99;
t28 = t103 * t58 - t57 * t99;
t8 = t27 * mrSges(7,1) - t28 * mrSges(7,2);
t113 = -t31 + t8;
t116 = mrSges(5,1) * t101 + mrSges(5,2) * t105;
t64 = t116 * t102;
t141 = -t64 + t113;
t97 = sin(pkin(11));
t86 = pkin(1) * t97 + pkin(7);
t140 = 0.2e1 * t86;
t139 = -Ifges(6,5) * t58 + Ifges(6,6) * t57;
t138 = -pkin(9) - pkin(8);
t106 = cos(qJ(3));
t137 = pkin(3) * t106;
t136 = pkin(4) * t100;
t135 = pkin(8) * t102;
t88 = pkin(4) * t104 + pkin(5);
t63 = t103 * t136 + t88 * t99;
t134 = t63 * mrSges(7,2);
t133 = t99 * mrSges(7,2);
t132 = -Ifges(7,3) - Ifges(6,3);
t131 = -Ifges(7,5) * t28 - Ifges(7,6) * t27;
t123 = t102 * t105;
t98 = cos(pkin(11));
t87 = -pkin(1) * t98 - pkin(2);
t65 = -t135 + t87 - t137;
t60 = t105 * t65;
t29 = -pkin(9) * t123 + t60 + (-t101 * t86 - pkin(4)) * t106;
t124 = t101 * t102;
t125 = t106 * t86;
t41 = t101 * t65 + t105 * t125;
t34 = -pkin(9) * t124 + t41;
t14 = t100 * t29 + t104 * t34;
t81 = t138 * t101;
t82 = t138 * t105;
t47 = t100 * t81 - t104 * t82;
t76 = -mrSges(5,1) * t105 + mrSges(5,2) * t101;
t130 = t76 - mrSges(4,1);
t80 = t102 * t86;
t61 = pkin(4) * t124 + t80;
t129 = t101 ^ 2 + t105 ^ 2;
t94 = t102 ^ 2;
t96 = t106 ^ 2;
t128 = t94 + t96;
t127 = Ifges(5,4) * t101;
t126 = Ifges(5,4) * t105;
t122 = pkin(5) * t133;
t121 = -Ifges(5,3) + t132;
t89 = -pkin(4) * t105 - pkin(3);
t120 = t129 * mrSges(5,3);
t13 = -t100 * t34 + t104 * t29;
t46 = t100 * t82 + t104 * t81;
t62 = t103 * t88 - t99 * t136;
t56 = t62 * mrSges(7,1);
t119 = Ifges(7,3) + t56 - t134;
t4 = -pkin(5) * t106 - pkin(10) * t58 + t13;
t5 = -pkin(10) * t57 + t14;
t2 = t103 * t4 - t5 * t99;
t3 = t103 * t5 + t4 * t99;
t118 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t131;
t32 = -pkin(10) * t69 + t46;
t33 = pkin(10) * t68 + t47;
t11 = t103 * t32 - t33 * t99;
t12 = t103 * t33 + t32 * t99;
t38 = t103 * t68 - t69 * t99;
t36 = Ifges(7,6) * t38;
t39 = t103 * t69 + t68 * t99;
t37 = Ifges(7,5) * t39;
t117 = t11 * mrSges(7,1) - t12 * mrSges(7,2) + t36 + t37;
t40 = -t101 * t125 + t60;
t115 = -t40 * t101 + t41 * t105;
t114 = (mrSges(6,1) * t104 - mrSges(6,2) * t100) * pkin(4);
t112 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t118 - t139;
t66 = Ifges(6,6) * t68;
t67 = Ifges(6,5) * t69;
t111 = t46 * mrSges(6,1) - t47 * mrSges(6,2) + t117 + t66 + t67;
t92 = Ifges(5,5) * t101;
t91 = Ifges(5,6) * t105;
t90 = t103 * pkin(5) * mrSges(7,1);
t85 = t86 ^ 2;
t83 = Ifges(5,5) * t123;
t79 = t94 * t85;
t78 = Ifges(5,1) * t101 + t126;
t77 = Ifges(5,2) * t105 + t127;
t74 = -mrSges(5,1) * t106 - mrSges(5,3) * t123;
t73 = mrSges(5,2) * t106 - mrSges(5,3) * t124;
t55 = -Ifges(5,5) * t106 + (Ifges(5,1) * t105 - t127) * t102;
t54 = -Ifges(5,6) * t106 + (-Ifges(5,2) * t101 + t126) * t102;
t50 = -pkin(5) * t68 + t89;
t49 = -mrSges(6,1) * t106 - mrSges(6,3) * t58;
t48 = mrSges(6,2) * t106 - mrSges(6,3) * t57;
t44 = Ifges(6,1) * t69 + Ifges(6,4) * t68;
t43 = Ifges(6,4) * t69 + Ifges(6,2) * t68;
t42 = -mrSges(6,1) * t68 + mrSges(6,2) * t69;
t35 = pkin(5) * t57 + t61;
t26 = Ifges(6,1) * t58 - Ifges(6,4) * t57 - Ifges(6,5) * t106;
t25 = Ifges(6,4) * t58 - Ifges(6,2) * t57 - Ifges(6,6) * t106;
t19 = -mrSges(7,1) * t106 - mrSges(7,3) * t28;
t18 = mrSges(7,2) * t106 + mrSges(7,3) * t27;
t17 = Ifges(7,1) * t39 + Ifges(7,4) * t38;
t16 = Ifges(7,4) * t39 + Ifges(7,2) * t38;
t15 = -mrSges(7,1) * t38 + mrSges(7,2) * t39;
t7 = Ifges(7,1) * t28 + Ifges(7,4) * t27 - Ifges(7,5) * t106;
t6 = Ifges(7,4) * t28 + Ifges(7,2) * t27 - Ifges(7,6) * t106;
t1 = [0.2e1 * t13 * t49 + 0.2e1 * t14 * t48 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 - t57 * t25 + t58 * t26 + t27 * t6 + t28 * t7 + 0.2e1 * t61 * t31 - 0.2e1 * t35 * t8 + 0.2e1 * t40 * t74 + 0.2e1 * t41 * t73 + Ifges(2,3) + Ifges(3,3) + (0.2e1 * t87 * mrSges(4,2) + Ifges(4,1) * t102 - t101 * t54 + t105 * t55 + t140 * t64) * t102 + (-0.2e1 * t87 * mrSges(4,1) - t83 + (Ifges(4,2) - t121) * t106 + (Ifges(5,6) * t101 + (2 * Ifges(4,4))) * t102 + t131 + t139) * t106 + m(4) * (t85 * t96 + t87 ^ 2 + t79) + m(5) * (t40 ^ 2 + t41 ^ 2 + t79) + m(6) * (t13 ^ 2 + t14 ^ 2 + t61 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t35 ^ 2) + m(3) * (t97 ^ 2 + t98 ^ 2) * pkin(1) ^ 2 + 0.2e1 * (t98 * mrSges(3,1) - t97 * mrSges(3,2)) * pkin(1) + t128 * mrSges(4,3) * t140; t28 * t18 + t27 * t19 + t58 * t48 - t57 * t49 + t141 * t106 + m(7) * (-t106 * t35 + t2 * t27 + t28 * t3) + m(6) * (-t106 * t61 - t13 * t57 + t14 * t58) + (t105 * t73 - t101 * t74 + m(5) * (t115 - t125)) * t102; m(3) + m(7) * (t27 ^ 2 + t28 ^ 2 + t96) + m(6) * (t57 ^ 2 + t58 ^ 2 + t96) + m(5) * (t129 * t94 + t96) + m(4) * t128; (t54 / 0.2e1 + pkin(8) * t73 + t41 * mrSges(5,3)) * t105 + (t55 / 0.2e1 - pkin(8) * t74 - t40 * mrSges(5,3)) * t101 + m(6) * (t13 * t46 + t14 * t47 + t61 * t89) + m(7) * (t11 * t2 + t12 * t3 + t35 * t50) + t89 * t31 + t68 * t25 / 0.2e1 + t69 * t26 / 0.2e1 + t58 * t44 / 0.2e1 + t61 * t42 - pkin(3) * t64 + t47 * t48 + t46 * t49 - t50 * t8 - t57 * t43 / 0.2e1 + t35 * t15 + t38 * t6 / 0.2e1 + t39 * t7 / 0.2e1 + t27 * t16 / 0.2e1 + t28 * t17 / 0.2e1 + (Ifges(4,5) + t105 * t78 / 0.2e1 - t101 * t77 / 0.2e1 + t130 * t86) * t102 + m(5) * (-pkin(3) * t80 + t115 * pkin(8)) + t12 * t18 + t11 * t19 + (-t13 * t69 + t14 * t68) * mrSges(6,3) + (-t2 * t39 + t3 * t38) * mrSges(7,3) + (-t86 * mrSges(4,2) - t92 / 0.2e1 - t91 / 0.2e1 - t67 / 0.2e1 - t66 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 + Ifges(4,6)) * t106; (-t27 * t39 + t28 * t38) * mrSges(7,3) + (t57 * t69 + t58 * t68) * mrSges(6,3) + (-mrSges(4,2) + t120) * t102 + (-t15 - t42 - t130) * t106 + m(7) * (-t106 * t50 + t11 * t27 + t12 * t28) + m(6) * (-t106 * t89 - t46 * t57 + t47 * t58) + m(5) * (t129 * t135 + t137); -0.2e1 * pkin(3) * t76 + t101 * t78 + t105 * t77 + 0.2e1 * t50 * t15 + t38 * t16 + t39 * t17 + 0.2e1 * t89 * t42 + t68 * t43 + t69 * t44 + Ifges(4,3) + m(7) * (t11 ^ 2 + t12 ^ 2 + t50 ^ 2) + m(6) * (t46 ^ 2 + t47 ^ 2 + t89 ^ 2) + m(5) * (t129 * pkin(8) ^ 2 + pkin(3) ^ 2) + 0.2e1 * (-t11 * t39 + t12 * t38) * mrSges(7,3) + 0.2e1 * (-t46 * t69 + t47 * t68) * mrSges(6,3) + 0.2e1 * pkin(8) * t120; (t100 * t48 + t104 * t49 + m(6) * (t100 * t14 + t104 * t13)) * pkin(4) + m(7) * (t2 * t62 + t3 * t63) + t121 * t106 + t83 + t62 * t19 + t63 * t18 + t40 * mrSges(5,1) - t41 * mrSges(5,2) - Ifges(5,6) * t124 + t112; m(7) * (t27 * t62 + t28 * t63) + m(6) * (t100 * t58 - t104 * t57) * pkin(4) + t141; m(7) * (t11 * t62 + t12 * t63) + t92 + t91 - t116 * pkin(8) + (t38 * t63 - t39 * t62) * mrSges(7,3) + (m(6) * (t100 * t47 + t104 * t46) + (t100 * t68 - t104 * t69) * mrSges(6,3)) * pkin(4) + t111; -0.2e1 * t134 + 0.2e1 * t56 + 0.2e1 * t114 + m(7) * (t62 ^ 2 + t63 ^ 2) + m(6) * (t100 ^ 2 + t104 ^ 2) * pkin(4) ^ 2 - t121; t132 * t106 + (t99 * t18 + t103 * t19 + m(7) * (t103 * t2 + t3 * t99)) * pkin(5) + t112; m(7) * (t103 * t27 + t28 * t99) * pkin(5) + t113; (m(7) * (t103 * t11 + t12 * t99) + (-t103 * t39 + t38 * t99) * mrSges(7,3)) * pkin(5) + t111; Ifges(6,3) + t90 + t114 + (m(7) * (t103 * t62 + t63 * t99) - t133) * pkin(5) + t119; -0.2e1 * t122 + 0.2e1 * t90 + m(7) * (t103 ^ 2 + t99 ^ 2) * pkin(5) ^ 2 - t132; -Ifges(7,3) * t106 + t118; t8; t117; t119; Ifges(7,3) + t90 - t122; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
