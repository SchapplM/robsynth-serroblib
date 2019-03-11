% Calculate joint inertia matrix for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2019-03-09 22:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR3_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:01:57
% EndTime: 2019-03-09 22:01:59
% DurationCPUTime: 0.91s
% Computational Cost: add. (2138->247), mult. (3867->331), div. (0->0), fcn. (4165->8), ass. (0->103)
t81 = sin(qJ(6));
t85 = cos(qJ(6));
t120 = -t81 ^ 2 - t85 ^ 2;
t86 = cos(qJ(4));
t133 = pkin(3) * t86;
t70 = -pkin(4) - t133;
t67 = -pkin(10) + t70;
t145 = t120 * t67;
t137 = pkin(4) + pkin(10);
t144 = t120 * t137;
t83 = sin(qJ(3));
t135 = pkin(2) * t83;
t87 = cos(qJ(3));
t71 = pkin(2) * t87 + pkin(3);
t82 = sin(qJ(4));
t48 = t86 * t135 + t71 * t82;
t44 = qJ(5) + t48;
t59 = mrSges(7,1) * t81 + mrSges(7,2) * t85;
t35 = t44 * t59;
t47 = -t82 * t135 + t71 * t86;
t46 = -pkin(4) - t47;
t40 = t46 * mrSges(6,2);
t43 = t47 * mrSges(5,1);
t143 = t35 + t40 + t43;
t134 = pkin(3) * t82;
t68 = qJ(5) + t134;
t49 = t68 * t59;
t64 = t70 * mrSges(6,2);
t73 = mrSges(5,1) * t133;
t142 = t49 + t64 + t73;
t141 = t44 ^ 2;
t140 = t68 ^ 2;
t84 = sin(qJ(2));
t88 = cos(qJ(2));
t54 = -t83 * t84 + t87 * t88;
t55 = t83 * t88 + t84 * t87;
t33 = -t86 * t54 + t55 * t82;
t34 = t54 * t82 + t55 * t86;
t72 = -pkin(2) * t88 - pkin(1);
t39 = -pkin(3) * t54 + t72;
t97 = -qJ(5) * t34 + t39;
t11 = pkin(4) * t33 + t97;
t139 = -0.2e1 * t11;
t138 = 0.2e1 * t39;
t136 = -pkin(8) - pkin(7);
t132 = pkin(4) * mrSges(6,2);
t131 = mrSges(7,1) * t85;
t130 = Ifges(7,4) * t81;
t129 = Ifges(7,4) * t85;
t128 = t33 * t81;
t127 = t33 * t85;
t126 = t44 * mrSges(6,3);
t125 = t44 * t68;
t124 = t48 * mrSges(5,2);
t123 = t68 * mrSges(6,3);
t122 = mrSges(6,1) + mrSges(5,3);
t62 = t136 * t84;
t63 = t136 * t88;
t38 = t83 * t62 - t87 * t63;
t121 = t84 ^ 2 + t88 ^ 2;
t119 = qJ(5) * t44;
t118 = qJ(5) * t68;
t117 = mrSges(5,2) * t134;
t116 = Ifges(7,5) * t128 + Ifges(7,6) * t127 + Ifges(7,3) * t34;
t25 = pkin(9) * t54 + t38;
t37 = t62 * t87 + t63 * t83;
t98 = -pkin(9) * t55 + t37;
t14 = t25 * t82 - t86 * t98;
t16 = t86 * t25 + t82 * t98;
t115 = t14 ^ 2 + t16 ^ 2;
t114 = m(7) * t120;
t111 = t120 * mrSges(7,3);
t42 = -pkin(10) + t46;
t110 = t120 * t42;
t74 = Ifges(7,5) * t85;
t107 = -Ifges(7,6) * t81 + t74;
t4 = t137 * t33 + t97;
t5 = pkin(5) * t34 + t14;
t1 = -t4 * t81 + t5 * t85;
t2 = t4 * t85 + t5 * t81;
t106 = t1 * t85 + t2 * t81;
t105 = mrSges(6,2) + t111;
t104 = -t81 * mrSges(7,2) + t131;
t19 = mrSges(7,1) * t34 - mrSges(7,3) * t128;
t20 = -mrSges(7,2) * t34 + mrSges(7,3) * t127;
t103 = t85 * t19 + t81 * t20;
t102 = 0.2e1 * t111;
t60 = -Ifges(7,2) * t81 + t129;
t61 = Ifges(7,1) * t85 - t130;
t101 = -t81 * t60 + t85 * t61 + Ifges(6,1) + Ifges(5,3);
t100 = Ifges(4,3) + t101;
t99 = (mrSges(4,1) * t87 - mrSges(4,2) * t83) * pkin(2);
t58 = qJ(5) * t59;
t80 = qJ(5) * mrSges(6,3);
t96 = t101 + t58 + t80 - t132;
t10 = Ifges(7,5) * t34 + (Ifges(7,1) * t81 + t129) * t33;
t6 = -t33 * pkin(5) + t16;
t9 = Ifges(7,6) * t34 + (Ifges(7,2) * t85 + t130) * t33;
t95 = t61 * t128 / 0.2e1 + t60 * t127 / 0.2e1 + t6 * t59 - t81 * t9 / 0.2e1 + t85 * t10 / 0.2e1 + (-mrSges(5,2) + mrSges(6,3)) * t16 - t106 * mrSges(7,3) + (Ifges(6,5) - Ifges(5,6)) * t33 + (mrSges(6,2) - mrSges(5,1)) * t14 + (t107 / 0.2e1 + Ifges(5,5) - Ifges(6,4)) * t34;
t94 = t37 * mrSges(4,1) - t38 * mrSges(4,2) + Ifges(4,5) * t55 + Ifges(4,6) * t54 + t95;
t90 = qJ(5) ^ 2;
t18 = t104 * t33;
t3 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t88 + mrSges(3,2) * t84) + t88 * (Ifges(3,4) * t84 + Ifges(3,2) * t88) + t84 * (Ifges(3,1) * t84 + Ifges(3,4) * t88) + 0.2e1 * t72 * (-mrSges(4,1) * t54 + mrSges(4,2) * t55) + t54 * (Ifges(4,4) * t55 + Ifges(4,2) * t54) + t55 * (Ifges(4,1) * t55 + Ifges(4,4) * t54) + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 - 0.2e1 * t6 * t18 + Ifges(2,3) + 0.2e1 * (-t37 * t55 + t38 * t54) * mrSges(4,3) + 0.2e1 * t121 * pkin(7) * mrSges(3,3) + (mrSges(5,2) * t138 + mrSges(6,3) * t139 + (Ifges(5,1) + Ifges(6,2)) * t34 + 0.2e1 * t122 * t14 + t116) * t34 + (mrSges(5,1) * t138 + mrSges(6,2) * t139 + t81 * t10 + t85 * t9 + (Ifges(5,2) + Ifges(6,3)) * t33 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t34 - 0.2e1 * t122 * t16) * t33 + m(3) * (t121 * pkin(7) ^ 2 + pkin(1) ^ 2) + m(4) * (t37 ^ 2 + t38 ^ 2 + t72 ^ 2) + m(5) * (t39 ^ 2 + t115) + m(6) * (t11 ^ 2 + t115) + m(7) * (t1 ^ 2 + t2 ^ 2 + t6 ^ 2); m(7) * (t106 * t42 + t44 * t6) + t94 + Ifges(3,6) * t88 + Ifges(3,5) * t84 - t44 * t18 + t103 * t42 + m(6) * (t14 * t46 + t16 * t44) + m(5) * (-t14 * t47 + t16 * t48) + (-t84 * mrSges(3,1) - t88 * mrSges(3,2)) * pkin(7) + (-t48 * t33 - t47 * t34) * mrSges(5,3) + (-t44 * t33 + t46 * t34) * mrSges(6,1) + (m(4) * (t37 * t87 + t38 * t83) + (t54 * t83 - t55 * t87) * mrSges(4,3)) * pkin(2); -0.2e1 * t124 + 0.2e1 * t126 + Ifges(3,3) + 0.2e1 * t35 + 0.2e1 * t40 + 0.2e1 * t43 + 0.2e1 * t99 + t42 * t102 + m(7) * (-t120 * t42 ^ 2 + t141) + m(6) * (t46 ^ 2 + t141) + m(5) * (t47 ^ 2 + t48 ^ 2) + m(4) * (t83 ^ 2 + t87 ^ 2) * pkin(2) ^ 2 + t100; (-t68 * t33 + t70 * t34) * mrSges(6,1) + (m(5) * (-t14 * t86 + t16 * t82) + (-t82 * t33 - t86 * t34) * mrSges(5,3)) * pkin(3) + m(7) * (t106 * t67 + t6 * t68) + t94 - t68 * t18 + t103 * t67 + m(6) * (t14 * t70 + t16 * t68); t99 + (t44 + t68) * mrSges(6,3) + (-t48 - t134) * mrSges(5,2) + m(7) * (-t145 * t42 + t125) + m(6) * (t46 * t70 + t125) + m(5) * (t47 * t86 + t48 * t82) * pkin(3) + (t145 + t110) * mrSges(7,3) + t100 + t142 + t143; -0.2e1 * t117 + 0.2e1 * t123 + 0.2e1 * t49 + 0.2e1 * t64 + 0.2e1 * t73 + t67 * t102 + m(7) * (-t120 * t67 ^ 2 + t140) + m(6) * (t70 ^ 2 + t140) + m(5) * (t82 ^ 2 + t86 ^ 2) * pkin(3) ^ 2 + t100; -t103 * t137 + m(6) * (-pkin(4) * t14 + qJ(5) * t16) + (-pkin(4) * t34 - qJ(5) * t33) * mrSges(6,1) + m(7) * (qJ(5) * t6 - t106 * t137) + t95 - qJ(5) * t18; -t124 + t126 + m(7) * (t144 * t42 + t119) + m(6) * (-pkin(4) * t46 + t119) + (-t144 + t110) * mrSges(7,3) + t96 + t143; -t117 + t123 + m(7) * (t144 * t67 + t118) + m(6) * (-pkin(4) * t70 + t118) + (-t144 + t145) * mrSges(7,3) + t96 + t142; -0.2e1 * t132 + 0.2e1 * t58 + 0.2e1 * t80 - t137 * t102 + m(7) * (-t120 * t137 ^ 2 + t90) + m(6) * (pkin(4) ^ 2 + t90) + t101; m(6) * t14 + m(7) * t106 + t34 * mrSges(6,1) + t103; m(6) * t46 - t42 * t114 + t105; m(6) * t70 - m(7) * t145 + t105; -m(6) * pkin(4) + m(7) * t144 + t105; m(6) - t114; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t116; t104 * t42 + t107; t104 * t67 + t107; -t137 * t131 + t74 + (mrSges(7,2) * t137 - Ifges(7,6)) * t81; t104; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
