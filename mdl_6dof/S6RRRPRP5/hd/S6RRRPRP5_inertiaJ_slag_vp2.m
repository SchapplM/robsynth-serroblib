% Calculate joint inertia matrix for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:47:47
% EndTime: 2019-03-09 16:47:50
% DurationCPUTime: 1.22s
% Computational Cost: add. (2385->326), mult. (4795->445), div. (0->0), fcn. (5031->8), ass. (0->115)
t155 = 2 * pkin(7);
t125 = sin(qJ(2));
t126 = cos(qJ(3));
t141 = t125 * t126;
t121 = sin(pkin(10));
t122 = cos(pkin(10));
t124 = sin(qJ(3));
t93 = t121 * t126 + t122 * t124;
t82 = t93 * t125;
t92 = -t121 * t124 + t122 * t126;
t83 = t92 * t125;
t154 = -Ifges(4,5) * t141 - Ifges(5,5) * t83 + Ifges(5,6) * t82;
t153 = 2 * mrSges(7,3);
t152 = cos(qJ(5));
t151 = pkin(3) * t121;
t127 = cos(qJ(2));
t150 = pkin(7) * t127;
t116 = t125 * pkin(7);
t110 = pkin(3) * t122 + pkin(4);
t123 = sin(qJ(5));
t85 = t110 * t152 - t123 * t151;
t149 = t85 * mrSges(6,1);
t86 = t123 * t110 + t152 * t151;
t148 = t86 * mrSges(6,2);
t146 = Ifges(7,2) + Ifges(6,3);
t145 = -qJ(4) - pkin(8);
t101 = -pkin(2) * t127 - pkin(8) * t125 - pkin(1);
t95 = t126 * t101;
t61 = -qJ(4) * t141 + t95 + (-pkin(7) * t124 - pkin(3)) * t127;
t142 = t124 * t125;
t74 = t124 * t101 + t126 * t150;
t67 = -qJ(4) * t142 + t74;
t31 = -t121 * t67 + t122 * t61;
t17 = -pkin(4) * t127 - pkin(9) * t83 + t31;
t32 = t121 * t61 + t122 * t67;
t24 = -pkin(9) * t82 + t32;
t6 = t123 * t17 + t152 * t24;
t137 = t145 * t124;
t138 = t145 * t126;
t69 = t121 * t137 - t122 * t138;
t144 = Ifges(4,4) * t124;
t143 = Ifges(4,4) * t126;
t48 = -t123 * t82 + t152 * t83;
t36 = t127 * mrSges(7,1) + t48 * mrSges(7,2);
t100 = pkin(3) * t142 + t116;
t140 = t124 ^ 2 + t126 ^ 2;
t68 = t121 * t138 + t122 * t137;
t131 = -t93 * pkin(9) + t68;
t49 = pkin(9) * t92 + t69;
t21 = t123 * t49 - t131 * t152;
t23 = t123 * t131 + t152 * t49;
t139 = t21 ^ 2 + t23 ^ 2;
t111 = -pkin(3) * t126 - pkin(2);
t50 = t82 * mrSges(5,1) + t83 * mrSges(5,2);
t64 = -t92 * mrSges(5,1) + t93 * mrSges(5,2);
t47 = t123 * t83 + t152 * t82;
t14 = t47 * mrSges(6,1) + t48 * mrSges(6,2);
t59 = t123 * t93 - t152 * t92;
t60 = t123 * t92 + t152 * t93;
t26 = t59 * mrSges(6,1) + t60 * mrSges(6,2);
t13 = t47 * mrSges(7,1) - t48 * mrSges(7,3);
t25 = t59 * mrSges(7,1) - t60 * mrSges(7,3);
t136 = Ifges(5,3) + Ifges(4,3) + t146;
t135 = (-Ifges(7,4) - Ifges(6,5)) * t48 + (Ifges(6,6) - Ifges(7,6)) * t47;
t62 = pkin(4) * t82 + t100;
t134 = mrSges(4,1) * t124 + mrSges(4,2) * t126;
t75 = -pkin(4) * t92 + t111;
t5 = -t123 * t24 + t152 * t17;
t2 = -qJ(6) * t127 + t6;
t3 = t127 * pkin(5) - t5;
t133 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) - t135;
t55 = Ifges(7,6) * t59;
t56 = Ifges(6,6) * t59;
t57 = Ifges(6,5) * t60;
t58 = Ifges(7,4) * t60;
t132 = t55 - t56 + t57 + t58 + (-mrSges(6,2) + mrSges(7,3)) * t23 + (-mrSges(6,1) - mrSges(7,1)) * t21;
t129 = pkin(7) ^ 2;
t120 = t127 ^ 2;
t118 = t125 ^ 2;
t115 = t118 * t129;
t114 = Ifges(4,5) * t124;
t113 = Ifges(4,6) * t126;
t104 = Ifges(4,1) * t124 + t143;
t103 = Ifges(4,2) * t126 + t144;
t102 = -mrSges(4,1) * t126 + mrSges(4,2) * t124;
t99 = -mrSges(4,1) * t127 - mrSges(4,3) * t141;
t98 = mrSges(4,2) * t127 - mrSges(4,3) * t142;
t91 = t134 * t125;
t90 = Ifges(5,5) * t93;
t89 = Ifges(5,6) * t92;
t84 = -pkin(5) - t85;
t81 = qJ(6) + t86;
t80 = -Ifges(4,5) * t127 + (Ifges(4,1) * t126 - t144) * t125;
t79 = -Ifges(4,6) * t127 + (-Ifges(4,2) * t124 + t143) * t125;
t73 = -t124 * t150 + t95;
t71 = -mrSges(5,1) * t127 - mrSges(5,3) * t83;
t70 = mrSges(5,2) * t127 - mrSges(5,3) * t82;
t66 = Ifges(5,1) * t93 + Ifges(5,4) * t92;
t65 = Ifges(5,4) * t93 + Ifges(5,2) * t92;
t45 = Ifges(5,1) * t83 - Ifges(5,4) * t82 - Ifges(5,5) * t127;
t44 = Ifges(5,4) * t83 - Ifges(5,2) * t82 - Ifges(5,6) * t127;
t35 = -mrSges(6,1) * t127 - mrSges(6,3) * t48;
t34 = mrSges(6,2) * t127 - mrSges(6,3) * t47;
t33 = -mrSges(7,2) * t47 - mrSges(7,3) * t127;
t30 = Ifges(6,1) * t60 - Ifges(6,4) * t59;
t29 = Ifges(7,1) * t60 + Ifges(7,5) * t59;
t28 = Ifges(6,4) * t60 - Ifges(6,2) * t59;
t27 = Ifges(7,5) * t60 + Ifges(7,3) * t59;
t18 = pkin(5) * t59 - qJ(6) * t60 + t75;
t11 = Ifges(6,1) * t48 - Ifges(6,4) * t47 - Ifges(6,5) * t127;
t10 = Ifges(7,1) * t48 - Ifges(7,4) * t127 + Ifges(7,5) * t47;
t9 = Ifges(6,4) * t48 - Ifges(6,2) * t47 - Ifges(6,6) * t127;
t8 = Ifges(7,5) * t48 - Ifges(7,6) * t127 + Ifges(7,3) * t47;
t7 = pkin(5) * t47 - qJ(6) * t48 + t62;
t1 = [0.2e1 * t100 * t50 + 0.2e1 * t7 * t13 + 0.2e1 * t62 * t14 + 0.2e1 * t2 * t33 + 0.2e1 * t3 * t36 + 0.2e1 * t31 * t71 + 0.2e1 * t32 * t70 + 0.2e1 * t6 * t34 + 0.2e1 * t5 * t35 - t82 * t44 + t83 * t45 + 0.2e1 * t73 * t99 + 0.2e1 * t74 * t98 + Ifges(2,3) + (t10 + t11) * t48 + (t8 - t9) * t47 + (t118 + t120) * mrSges(3,3) * t155 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t125 - t124 * t79 + t126 * t80 + t91 * t155) * t125 + m(3) * (pkin(1) ^ 2 + t120 * t129 + t115) + m(4) * (t73 ^ 2 + t74 ^ 2 + t115) + m(5) * (t100 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t5 ^ 2 + t6 ^ 2 + t62 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) + (Ifges(3,2) + t136) * t127 + (Ifges(4,6) * t124 + (2 * Ifges(3,4))) * t125 + t135 + t154) * t127; (t79 / 0.2e1 + pkin(8) * t98 + t74 * mrSges(4,3)) * t126 + (t80 / 0.2e1 - pkin(8) * t99 - t73 * mrSges(4,3)) * t124 + (t33 + t34) * t23 + (-t35 + t36) * t21 + t100 * t64 + t111 * t50 - pkin(2) * t91 + t92 * t44 / 0.2e1 + t93 * t45 / 0.2e1 - t82 * t65 / 0.2e1 + t83 * t66 / 0.2e1 + t69 * t70 + t68 * t71 + t75 * t14 + t62 * t26 + t7 * t25 + t18 * t13 + (-t114 / 0.2e1 - t113 / 0.2e1 - t90 / 0.2e1 - t89 / 0.2e1 - t57 / 0.2e1 + t56 / 0.2e1 - t58 / 0.2e1 - t55 / 0.2e1 + Ifges(3,6) - pkin(7) * mrSges(3,2)) * t127 + m(5) * (t100 * t111 + t31 * t68 + t32 * t69) + m(6) * (-t21 * t5 + t23 * t6 + t62 * t75) + m(7) * (t18 * t7 + t2 * t23 + t21 * t3) + m(4) * (-pkin(2) * t116 + (-t73 * t124 + t74 * t126) * pkin(8)) + (Ifges(3,5) - t124 * t103 / 0.2e1 + t126 * t104 / 0.2e1 + (t102 - mrSges(3,1)) * pkin(7)) * t125 + (t3 * mrSges(7,2) - t5 * mrSges(6,3) + t10 / 0.2e1 + t11 / 0.2e1) * t60 + (-t2 * mrSges(7,2) - t6 * mrSges(6,3) + t8 / 0.2e1 - t9 / 0.2e1) * t59 + (t29 / 0.2e1 + t30 / 0.2e1) * t48 + (t27 / 0.2e1 - t28 / 0.2e1) * t47 + (-t31 * t93 + t32 * t92) * mrSges(5,3); -0.2e1 * pkin(2) * t102 + t126 * t103 + t124 * t104 + 0.2e1 * t111 * t64 + 0.2e1 * t18 * t25 + 0.2e1 * t75 * t26 + t92 * t65 + t93 * t66 + Ifges(3,3) + (t29 + t30) * t60 + (t27 - t28) * t59 + m(4) * (pkin(8) ^ 2 * t140 + pkin(2) ^ 2) + m(5) * (t111 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t75 ^ 2 + t139) + m(7) * (t18 ^ 2 + t139) + 0.2e1 * (-t68 * t93 + t69 * t92) * mrSges(5,3) + 0.2e1 * t140 * pkin(8) * mrSges(4,3) + 0.2e1 * (t21 * t60 - t23 * t59) * (mrSges(7,2) + mrSges(6,3)); t86 * t34 - Ifges(4,6) * t142 + t133 + t81 * t33 + t84 * t36 + t85 * t35 + t73 * mrSges(4,1) - t74 * mrSges(4,2) + t31 * mrSges(5,1) - t32 * mrSges(5,2) + (t121 * t70 + m(5) * (t121 * t32 + t122 * t31) + t122 * t71) * pkin(3) - t136 * t127 + m(7) * (t2 * t81 + t3 * t84) + m(6) * (t5 * t85 + t6 * t86) - t154; m(6) * (-t21 * t85 + t23 * t86) + m(7) * (t21 * t84 + t23 * t81) - t134 * pkin(8) + (-t59 * t86 - t60 * t85) * mrSges(6,3) + (-t59 * t81 + t60 * t84) * mrSges(7,2) + t113 + t114 + t90 + t89 + t132 + t68 * mrSges(5,1) - t69 * mrSges(5,2) + (m(5) * (t121 * t69 + t122 * t68) + (t121 * t92 - t122 * t93) * mrSges(5,3)) * pkin(3); 0.2e1 * t149 - 0.2e1 * t84 * mrSges(7,1) - 0.2e1 * t148 + t81 * t153 + m(7) * (t81 ^ 2 + t84 ^ 2) + m(6) * (t85 ^ 2 + t86 ^ 2) + t136 + (0.2e1 * mrSges(5,1) * t122 - 0.2e1 * mrSges(5,2) * t121 + m(5) * (t121 ^ 2 + t122 ^ 2) * pkin(3)) * pkin(3); m(5) * t100 + m(6) * t62 + m(7) * t7 + t13 + t14 + t50; m(5) * t111 + m(6) * t75 + m(7) * t18 + t25 + t26 + t64; 0; m(5) + m(6) + m(7); -pkin(5) * t36 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) + qJ(6) * t33 - t146 * t127 + t133; m(7) * (-pkin(5) * t21 + qJ(6) * t23) + (-pkin(5) * t60 - qJ(6) * t59) * mrSges(7,2) + t132; t149 + m(7) * (-pkin(5) * t84 + qJ(6) * t81) - t148 + (t81 + qJ(6)) * mrSges(7,3) + (-t84 + pkin(5)) * mrSges(7,1) + t146; 0; 0.2e1 * pkin(5) * mrSges(7,1) + qJ(6) * t153 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) + t146; m(7) * t3 + t36; m(7) * t21 + t60 * mrSges(7,2); m(7) * t84 - mrSges(7,1); 0; -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
