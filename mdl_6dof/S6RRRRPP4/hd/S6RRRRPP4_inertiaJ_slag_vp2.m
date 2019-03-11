% Calculate joint inertia matrix for
% S6RRRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 21:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:58:31
% EndTime: 2019-03-09 20:58:34
% DurationCPUTime: 1.28s
% Computational Cost: add. (2540->339), mult. (5090->469), div. (0->0), fcn. (5309->8), ass. (0->120)
t160 = 2 * pkin(7);
t159 = -2 * mrSges(7,1);
t158 = 2 * mrSges(7,3);
t157 = -pkin(9) - pkin(8);
t125 = sin(qJ(4));
t156 = pkin(3) * t125;
t130 = cos(qJ(2));
t155 = pkin(7) * t130;
t127 = sin(qJ(2));
t118 = t127 * pkin(7);
t128 = cos(qJ(4));
t112 = pkin(3) * t128 + pkin(4);
t123 = sin(pkin(10));
t124 = cos(pkin(10));
t86 = t112 * t124 - t123 * t156;
t154 = t86 * mrSges(6,1);
t87 = t112 * t123 + t124 * t156;
t153 = t87 * mrSges(6,2);
t126 = sin(qJ(3));
t129 = cos(qJ(3));
t148 = t127 * t129;
t101 = -pkin(2) * t130 - pkin(8) * t127 - pkin(1);
t93 = t129 * t101;
t62 = -pkin(9) * t148 + t93 + (-pkin(7) * t126 - pkin(3)) * t130;
t149 = t126 * t127;
t76 = t101 * t126 + t129 * t155;
t68 = -pkin(9) * t149 + t76;
t32 = -t125 * t68 + t128 * t62;
t94 = -t125 * t126 + t128 * t129;
t85 = t94 * t127;
t15 = -pkin(4) * t130 - qJ(5) * t85 + t32;
t33 = t125 * t62 + t128 * t68;
t95 = t125 * t129 + t126 * t128;
t84 = t95 * t127;
t24 = -qJ(5) * t84 + t33;
t6 = t123 * t15 + t124 * t24;
t143 = t157 * t126;
t144 = t157 * t129;
t71 = t125 * t143 - t128 * t144;
t151 = Ifges(4,4) * t126;
t150 = Ifges(4,4) * t129;
t48 = -t123 * t84 + t124 * t85;
t37 = mrSges(7,1) * t130 + mrSges(7,2) * t48;
t100 = pkin(3) * t149 + t118;
t147 = t126 ^ 2 + t129 ^ 2;
t70 = t125 * t144 + t128 * t143;
t137 = -t95 * qJ(5) + t70;
t50 = qJ(5) * t94 + t71;
t21 = t123 * t50 - t124 * t137;
t23 = t123 * t137 + t124 * t50;
t146 = t21 ^ 2 + t23 ^ 2;
t145 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t113 = -pkin(3) * t129 - pkin(2);
t47 = t123 * t85 + t124 * t84;
t14 = t47 * mrSges(6,1) + mrSges(6,2) * t48;
t60 = t123 * t95 - t124 * t94;
t61 = t123 * t94 + t124 * t95;
t26 = t60 * mrSges(6,1) + mrSges(6,2) * t61;
t13 = mrSges(7,1) * t47 - t48 * mrSges(7,3);
t25 = mrSges(7,1) * t60 - t61 * mrSges(7,3);
t142 = Ifges(4,3) + t145;
t63 = pkin(4) * t84 + t100;
t141 = mrSges(4,1) * t126 + mrSges(4,2) * t129;
t140 = t124 * mrSges(6,1) - t123 * mrSges(6,2);
t5 = -t123 * t24 + t124 * t15;
t139 = -Ifges(5,5) * t85 + Ifges(5,6) * t84 + (-Ifges(7,4) - Ifges(6,5)) * t48 + (Ifges(6,6) - Ifges(7,6)) * t47;
t77 = -pkin(4) * t94 + t113;
t138 = (mrSges(5,1) * t128 - mrSges(5,2) * t125) * pkin(3);
t2 = -qJ(6) * t130 + t6;
t3 = pkin(5) * t130 - t5;
t136 = mrSges(5,1) * t32 + mrSges(6,1) * t5 - t3 * mrSges(7,1) - t33 * mrSges(5,2) - t6 * mrSges(6,2) + mrSges(7,3) * t2 - t139;
t56 = Ifges(7,6) * t60;
t57 = Ifges(6,6) * t60;
t58 = Ifges(6,5) * t61;
t59 = Ifges(7,4) * t61;
t90 = Ifges(5,6) * t94;
t91 = Ifges(5,5) * t95;
t135 = t70 * mrSges(5,1) - t71 * mrSges(5,2) + t56 - t57 + t58 + t59 + t90 + t91 + (-mrSges(6,2) + mrSges(7,3)) * t23 + (-mrSges(6,1) - mrSges(7,1)) * t21;
t132 = pkin(7) ^ 2;
t122 = t130 ^ 2;
t120 = t127 ^ 2;
t117 = t120 * t132;
t116 = Ifges(4,5) * t126;
t115 = Ifges(4,6) * t129;
t111 = -pkin(4) * t124 - pkin(5);
t110 = pkin(4) * t123 + qJ(6);
t107 = Ifges(4,5) * t148;
t104 = Ifges(4,1) * t126 + t150;
t103 = Ifges(4,2) * t129 + t151;
t102 = -mrSges(4,1) * t129 + mrSges(4,2) * t126;
t99 = -mrSges(4,1) * t130 - mrSges(4,3) * t148;
t98 = mrSges(4,2) * t130 - mrSges(4,3) * t149;
t89 = t141 * t127;
t83 = -pkin(5) - t86;
t82 = qJ(6) + t87;
t81 = -Ifges(4,5) * t130 + (Ifges(4,1) * t129 - t151) * t127;
t80 = -Ifges(4,6) * t130 + (-Ifges(4,2) * t126 + t150) * t127;
t75 = -t126 * t155 + t93;
t73 = -mrSges(5,1) * t130 - mrSges(5,3) * t85;
t72 = mrSges(5,2) * t130 - mrSges(5,3) * t84;
t67 = Ifges(5,1) * t95 + Ifges(5,4) * t94;
t66 = Ifges(5,4) * t95 + Ifges(5,2) * t94;
t65 = -mrSges(5,1) * t94 + mrSges(5,2) * t95;
t51 = mrSges(5,1) * t84 + mrSges(5,2) * t85;
t46 = Ifges(5,1) * t85 - Ifges(5,4) * t84 - Ifges(5,5) * t130;
t45 = Ifges(5,4) * t85 - Ifges(5,2) * t84 - Ifges(5,6) * t130;
t36 = -mrSges(6,1) * t130 - mrSges(6,3) * t48;
t35 = mrSges(6,2) * t130 - mrSges(6,3) * t47;
t34 = -mrSges(7,2) * t47 - mrSges(7,3) * t130;
t30 = Ifges(6,1) * t61 - Ifges(6,4) * t60;
t29 = Ifges(7,1) * t61 + Ifges(7,5) * t60;
t28 = Ifges(6,4) * t61 - Ifges(6,2) * t60;
t27 = Ifges(7,5) * t61 + Ifges(7,3) * t60;
t17 = pkin(5) * t60 - qJ(6) * t61 + t77;
t11 = Ifges(6,1) * t48 - Ifges(6,4) * t47 - Ifges(6,5) * t130;
t10 = Ifges(7,1) * t48 - Ifges(7,4) * t130 + Ifges(7,5) * t47;
t9 = Ifges(6,4) * t48 - Ifges(6,2) * t47 - Ifges(6,6) * t130;
t8 = Ifges(7,5) * t48 - Ifges(7,6) * t130 + Ifges(7,3) * t47;
t7 = pkin(5) * t47 - qJ(6) * t48 + t63;
t1 = [0.2e1 * t100 * t51 + 0.2e1 * t7 * t13 + 0.2e1 * t63 * t14 + 0.2e1 * t2 * t34 + 0.2e1 * t3 * t37 + 0.2e1 * t32 * t73 + 0.2e1 * t33 * t72 + 0.2e1 * t6 * t35 + 0.2e1 * t5 * t36 - t84 * t45 + t85 * t46 + 0.2e1 * t75 * t99 + 0.2e1 * t76 * t98 + Ifges(2,3) + (t10 + t11) * t48 + (t8 - t9) * t47 + (t120 + t122) * mrSges(3,3) * t160 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t127 - t126 * t80 + t129 * t81 + t89 * t160) * t127 + m(3) * (pkin(1) ^ 2 + t122 * t132 + t117) + m(4) * (t75 ^ 2 + t76 ^ 2 + t117) + m(5) * (t100 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t63 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t107 + (Ifges(3,2) + t142) * t130 + (Ifges(4,6) * t126 + (2 * Ifges(3,4))) * t127 + t139) * t130; (t34 + t35) * t23 + (-t36 + t37) * t21 + m(5) * (t100 * t113 + t32 * t70 + t33 * t71) + m(6) * (-t21 * t5 + t23 * t6 + t63 * t77) + m(7) * (t17 * t7 + t2 * t23 + t21 * t3) + (t10 / 0.2e1 + t11 / 0.2e1 + t3 * mrSges(7,2) - t5 * mrSges(6,3)) * t61 + (t8 / 0.2e1 - t9 / 0.2e1 - t6 * mrSges(6,3) - t2 * mrSges(7,2)) * t60 + t100 * t65 + t113 * t51 - pkin(2) * t89 + t94 * t45 / 0.2e1 + t95 * t46 / 0.2e1 - t84 * t66 / 0.2e1 + t85 * t67 / 0.2e1 + t71 * t72 + t70 * t73 + t77 * t14 + t63 * t26 + (Ifges(3,5) + t129 * t104 / 0.2e1 - t126 * t103 / 0.2e1 + (t102 - mrSges(3,1)) * pkin(7)) * t127 + (t29 / 0.2e1 + t30 / 0.2e1) * t48 + (t27 / 0.2e1 - t28 / 0.2e1) * t47 + t7 * t25 + t17 * t13 + (-t116 / 0.2e1 - t115 / 0.2e1 - t91 / 0.2e1 - t90 / 0.2e1 - t58 / 0.2e1 + t57 / 0.2e1 - t59 / 0.2e1 - t56 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2)) * t130 + m(4) * (-pkin(2) * t118 + (-t126 * t75 + t129 * t76) * pkin(8)) + (-t32 * t95 + t33 * t94) * mrSges(5,3) + (t80 / 0.2e1 + pkin(8) * t98 + t76 * mrSges(4,3)) * t129 + (t81 / 0.2e1 - pkin(8) * t99 - t75 * mrSges(4,3)) * t126; -0.2e1 * pkin(2) * t102 + t129 * t103 + t126 * t104 + 0.2e1 * t113 * t65 + 0.2e1 * t17 * t25 + 0.2e1 * t77 * t26 + t94 * t66 + t95 * t67 + Ifges(3,3) + (t29 + t30) * t61 + (t27 - t28) * t60 + m(4) * (pkin(8) ^ 2 * t147 + pkin(2) ^ 2) + m(5) * (t113 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(6) * (t77 ^ 2 + t146) + m(7) * (t17 ^ 2 + t146) + 0.2e1 * (-t70 * t95 + t71 * t94) * mrSges(5,3) + 0.2e1 * t147 * pkin(8) * mrSges(4,3) + 0.2e1 * (t21 * t61 - t23 * t60) * (mrSges(7,2) + mrSges(6,3)); m(7) * (t2 * t82 + t3 * t83) + m(6) * (t5 * t86 + t6 * t87) + t107 - t142 * t130 + t87 * t35 + t82 * t34 + t83 * t37 + t86 * t36 + t75 * mrSges(4,1) - t76 * mrSges(4,2) - Ifges(4,6) * t149 + (m(5) * (t125 * t33 + t128 * t32) + t128 * t73 + t125 * t72) * pkin(3) + t136; -t141 * pkin(8) + (-t60 * t87 - t61 * t86) * mrSges(6,3) + (-t60 * t82 + t61 * t83) * mrSges(7,2) + t115 + t116 + (m(5) * (t125 * t71 + t128 * t70) + (t125 * t94 - t128 * t95) * mrSges(5,3)) * pkin(3) + t135 + m(7) * (t21 * t83 + t23 * t82) + m(6) * (-t21 * t86 + t23 * t87); 0.2e1 * t154 + t83 * t159 - 0.2e1 * t153 + t82 * t158 + 0.2e1 * t138 + m(6) * (t86 ^ 2 + t87 ^ 2) + m(7) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t125 ^ 2 + t128 ^ 2) * pkin(3) ^ 2 + t142; m(7) * (t110 * t2 + t111 * t3) - t145 * t130 + t110 * t34 + t111 * t37 + (t124 * t36 + t123 * t35 + m(6) * (t123 * t6 + t124 * t5)) * pkin(4) + t136; (-t110 * t60 + t111 * t61) * mrSges(7,2) + m(7) * (t110 * t23 + t111 * t21) + (m(6) * (t123 * t23 - t124 * t21) + (-t123 * t60 - t124 * t61) * mrSges(6,3)) * pkin(4) + t135; m(7) * (t110 * t82 + t111 * t83) - t153 + t154 + t138 + (t110 + t82) * mrSges(7,3) + (-t111 - t83) * mrSges(7,1) + (m(6) * (t123 * t87 + t124 * t86) + t140) * pkin(4) + t145; t111 * t159 + t110 * t158 + m(7) * (t110 ^ 2 + t111 ^ 2) + t145 + (0.2e1 * t140 + m(6) * (t123 ^ 2 + t124 ^ 2) * pkin(4)) * pkin(4); m(6) * t63 + m(7) * t7 + t13 + t14; m(6) * t77 + m(7) * t17 + t25 + t26; 0; 0; m(6) + m(7); m(7) * t3 + t37; m(7) * t21 + t61 * mrSges(7,2); m(7) * t83 - mrSges(7,1); m(7) * t111 - mrSges(7,1); 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
