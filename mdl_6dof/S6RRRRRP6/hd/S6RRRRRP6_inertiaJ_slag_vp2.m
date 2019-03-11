% Calculate joint inertia matrix for
% S6RRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP6_inertiaJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:26:37
% EndTime: 2019-03-10 01:26:41
% DurationCPUTime: 1.42s
% Computational Cost: add. (2809->359), mult. (5541->494), div. (0->0), fcn. (5745->8), ass. (0->134)
t170 = 2 * pkin(7);
t127 = sin(qJ(2));
t125 = sin(qJ(4));
t126 = sin(qJ(3));
t129 = cos(qJ(4));
t130 = cos(qJ(3));
t93 = t125 * t130 + t126 * t129;
t81 = t93 * t127;
t92 = -t125 * t126 + t129 * t130;
t82 = t92 * t127;
t169 = -Ifges(5,5) * t82 + Ifges(5,6) * t81;
t124 = sin(qJ(5));
t109 = pkin(4) * t124 + qJ(6);
t104 = t109 * mrSges(7,3);
t128 = cos(qJ(5));
t164 = pkin(4) * t128;
t113 = mrSges(6,1) * t164;
t168 = t104 + t113;
t123 = qJ(6) * mrSges(7,3);
t132 = pkin(5) * mrSges(7,1);
t167 = t123 + t132;
t166 = -pkin(9) - pkin(8);
t165 = pkin(3) * t125;
t131 = cos(qJ(2));
t163 = pkin(7) * t131;
t118 = t127 * pkin(7);
t111 = pkin(3) * t129 + pkin(4);
t84 = t111 * t128 - t124 * t165;
t83 = -pkin(5) - t84;
t162 = t83 * mrSges(7,1);
t85 = t124 * t111 + t128 * t165;
t161 = t85 * mrSges(6,2);
t159 = -Ifges(6,3) - Ifges(7,2);
t153 = t127 * t130;
t99 = -pkin(2) * t131 - pkin(8) * t127 - pkin(1);
t91 = t130 * t99;
t58 = -pkin(9) * t153 + t91 + (-pkin(7) * t126 - pkin(3)) * t131;
t154 = t126 * t127;
t72 = t126 * t99 + t130 * t163;
t64 = -pkin(9) * t154 + t72;
t32 = -t125 * t64 + t129 * t58;
t13 = -pkin(4) * t131 - pkin(10) * t82 + t32;
t33 = t125 * t58 + t129 * t64;
t24 = -pkin(10) * t81 + t33;
t6 = t124 * t13 + t128 * t24;
t147 = t166 * t126;
t148 = t166 * t130;
t67 = t125 * t147 - t129 * t148;
t158 = Ifges(4,4) * t126;
t157 = Ifges(4,4) * t130;
t110 = -pkin(5) - t164;
t156 = t110 * mrSges(7,1);
t155 = t124 * mrSges(6,2);
t46 = -t124 * t81 + t128 * t82;
t37 = t131 * mrSges(7,1) + t46 * mrSges(7,2);
t98 = pkin(3) * t154 + t118;
t152 = t126 ^ 2 + t130 ^ 2;
t66 = t125 * t148 + t129 * t147;
t139 = -t93 * pkin(10) + t66;
t49 = pkin(10) * t92 + t67;
t21 = t124 * t49 - t128 * t139;
t23 = t124 * t139 + t128 * t49;
t151 = t21 ^ 2 + t23 ^ 2;
t150 = pkin(4) * t155;
t149 = -Ifges(5,3) + t159;
t112 = -pkin(3) * t130 - pkin(2);
t146 = Ifges(4,3) - t149;
t45 = t124 * t82 + t128 * t81;
t145 = (-Ifges(7,4) - Ifges(6,5)) * t46 + (Ifges(6,6) - Ifges(7,6)) * t45;
t59 = pkin(4) * t81 + t98;
t144 = mrSges(4,1) * t126 + mrSges(4,2) * t130;
t5 = -t124 * t24 + t128 * t13;
t73 = -pkin(4) * t92 + t112;
t80 = qJ(6) + t85;
t74 = t80 * mrSges(7,3);
t79 = t84 * mrSges(6,1);
t143 = -t159 + t74 + t79 - t161;
t142 = (mrSges(5,1) * t129 - mrSges(5,2) * t125) * pkin(3);
t2 = -qJ(6) * t131 + t6;
t3 = pkin(5) * t131 - t5;
t141 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t6 * mrSges(6,2) + t2 * mrSges(7,3) - t145;
t56 = t124 * t93 - t128 * t92;
t52 = Ifges(7,6) * t56;
t53 = Ifges(6,6) * t56;
t57 = t124 * t92 + t128 * t93;
t54 = Ifges(6,5) * t57;
t55 = Ifges(7,4) * t57;
t140 = t52 - t53 + t54 + t55 + (-mrSges(6,2) + mrSges(7,3)) * t23 + (-mrSges(6,1) - mrSges(7,1)) * t21;
t138 = t32 * mrSges(5,1) - t33 * mrSges(5,2) + t141 - t169;
t88 = Ifges(5,6) * t92;
t89 = Ifges(5,5) * t93;
t137 = t66 * mrSges(5,1) - t67 * mrSges(5,2) + t140 + t88 + t89;
t134 = pkin(7) ^ 2;
t122 = t131 ^ 2;
t120 = t127 ^ 2;
t117 = t120 * t134;
t116 = Ifges(4,5) * t126;
t115 = Ifges(4,6) * t130;
t105 = Ifges(4,5) * t153;
t102 = Ifges(4,1) * t126 + t157;
t101 = Ifges(4,2) * t130 + t158;
t100 = -mrSges(4,1) * t130 + mrSges(4,2) * t126;
t97 = -mrSges(4,1) * t131 - mrSges(4,3) * t153;
t96 = mrSges(4,2) * t131 - mrSges(4,3) * t154;
t87 = t144 * t127;
t78 = -Ifges(4,5) * t131 + (Ifges(4,1) * t130 - t158) * t127;
t77 = -Ifges(4,6) * t131 + (-Ifges(4,2) * t126 + t157) * t127;
t71 = -t126 * t163 + t91;
t69 = -mrSges(5,1) * t131 - mrSges(5,3) * t82;
t68 = mrSges(5,2) * t131 - mrSges(5,3) * t81;
t63 = Ifges(5,1) * t93 + Ifges(5,4) * t92;
t62 = Ifges(5,4) * t93 + Ifges(5,2) * t92;
t61 = -mrSges(5,1) * t92 + mrSges(5,2) * t93;
t48 = mrSges(5,1) * t81 + mrSges(5,2) * t82;
t44 = Ifges(5,1) * t82 - Ifges(5,4) * t81 - Ifges(5,5) * t131;
t43 = Ifges(5,4) * t82 - Ifges(5,2) * t81 - Ifges(5,6) * t131;
t36 = -mrSges(6,1) * t131 - mrSges(6,3) * t46;
t35 = mrSges(6,2) * t131 - mrSges(6,3) * t45;
t34 = -mrSges(7,2) * t45 - mrSges(7,3) * t131;
t30 = Ifges(6,1) * t57 - Ifges(6,4) * t56;
t29 = Ifges(7,1) * t57 + Ifges(7,5) * t56;
t28 = Ifges(6,4) * t57 - Ifges(6,2) * t56;
t27 = Ifges(7,5) * t57 + Ifges(7,3) * t56;
t26 = mrSges(6,1) * t56 + mrSges(6,2) * t57;
t25 = mrSges(7,1) * t56 - mrSges(7,3) * t57;
t17 = pkin(5) * t56 - qJ(6) * t57 + t73;
t15 = mrSges(6,1) * t45 + mrSges(6,2) * t46;
t14 = mrSges(7,1) * t45 - mrSges(7,3) * t46;
t11 = Ifges(6,1) * t46 - Ifges(6,4) * t45 - Ifges(6,5) * t131;
t10 = Ifges(7,1) * t46 - Ifges(7,4) * t131 + Ifges(7,5) * t45;
t9 = Ifges(6,4) * t46 - Ifges(6,2) * t45 - Ifges(6,6) * t131;
t8 = Ifges(7,5) * t46 - Ifges(7,6) * t131 + Ifges(7,3) * t45;
t7 = pkin(5) * t45 - qJ(6) * t46 + t59;
t1 = [0.2e1 * t7 * t14 + 0.2e1 * t59 * t15 + 0.2e1 * t2 * t34 + 0.2e1 * t3 * t37 + 0.2e1 * t32 * t69 + 0.2e1 * t33 * t68 + 0.2e1 * t6 * t35 + 0.2e1 * t5 * t36 - t81 * t43 + t82 * t44 + 0.2e1 * t98 * t48 + 0.2e1 * t71 * t97 + 0.2e1 * t72 * t96 + Ifges(2,3) + (t10 + t11) * t46 + (t8 - t9) * t45 + (t120 + t122) * mrSges(3,3) * t170 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t127 - t126 * t77 + t130 * t78 + t87 * t170) * t127 + m(3) * (pkin(1) ^ 2 + t122 * t134 + t117) + m(4) * (t71 ^ 2 + t72 ^ 2 + t117) + m(5) * (t32 ^ 2 + t33 ^ 2 + t98 ^ 2) + m(6) * (t5 ^ 2 + t59 ^ 2 + t6 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t7 ^ 2) + (0.2e1 * pkin(1) * mrSges(3,1) - t105 + (Ifges(3,2) + t146) * t131 + (Ifges(4,6) * t126 + (2 * Ifges(3,4))) * t127 + t145 + t169) * t131; (t10 / 0.2e1 + t11 / 0.2e1 + t3 * mrSges(7,2) - t5 * mrSges(6,3)) * t57 + (t8 / 0.2e1 - t9 / 0.2e1 - t6 * mrSges(6,3) - t2 * mrSges(7,2)) * t56 + (t130 * t102 / 0.2e1 - t126 * t101 / 0.2e1 + Ifges(3,5) + (t100 - mrSges(3,1)) * pkin(7)) * t127 + (t29 / 0.2e1 + t30 / 0.2e1) * t46 + (t27 / 0.2e1 - t28 / 0.2e1) * t45 + t112 * t48 + t93 * t44 / 0.2e1 + t98 * t61 - pkin(2) * t87 + t92 * t43 / 0.2e1 - t81 * t62 / 0.2e1 + t82 * t63 / 0.2e1 + t59 * t26 + t67 * t68 + t66 * t69 + t73 * t15 + t7 * t25 + t17 * t14 + (-t55 / 0.2e1 - t52 / 0.2e1 + Ifges(3,6) - t116 / 0.2e1 - t115 / 0.2e1 - t89 / 0.2e1 - t88 / 0.2e1 - t54 / 0.2e1 + t53 / 0.2e1 - pkin(7) * mrSges(3,2)) * t131 + (-t32 * t93 + t33 * t92) * mrSges(5,3) + (t77 / 0.2e1 + pkin(8) * t96 + t72 * mrSges(4,3)) * t130 + (-pkin(8) * t97 - t71 * mrSges(4,3) + t78 / 0.2e1) * t126 + (t34 + t35) * t23 + m(4) * (-pkin(2) * t118 + (-t126 * t71 + t130 * t72) * pkin(8)) + (-t36 + t37) * t21 + m(5) * (t112 * t98 + t32 * t66 + t33 * t67) + m(6) * (-t21 * t5 + t23 * t6 + t59 * t73) + m(7) * (t17 * t7 + t2 * t23 + t21 * t3); -0.2e1 * pkin(2) * t100 + t130 * t101 + t126 * t102 + 0.2e1 * t112 * t61 + 0.2e1 * t17 * t25 + 0.2e1 * t73 * t26 + t92 * t62 + t93 * t63 + Ifges(3,3) + (t29 + t30) * t57 + (t27 - t28) * t56 + m(4) * (pkin(8) ^ 2 * t152 + pkin(2) ^ 2) + m(5) * (t112 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t73 ^ 2 + t151) + m(7) * (t17 ^ 2 + t151) + 0.2e1 * (-t66 * t93 + t67 * t92) * mrSges(5,3) + 0.2e1 * t152 * pkin(8) * mrSges(4,3) + 0.2e1 * (t21 * t57 - t23 * t56) * (mrSges(7,2) + mrSges(6,3)); m(7) * (t2 * t80 + t3 * t83) + m(6) * (t5 * t84 + t6 * t85) - t146 * t131 - Ifges(4,6) * t154 + t105 + t84 * t36 + t85 * t35 + t80 * t34 + t83 * t37 + t71 * mrSges(4,1) - t72 * mrSges(4,2) + (m(5) * (t125 * t33 + t129 * t32) + t129 * t69 + t125 * t68) * pkin(3) + t138; t115 + t116 + (m(5) * (t125 * t67 + t129 * t66) + (t125 * t92 - t129 * t93) * mrSges(5,3)) * pkin(3) + t137 + m(7) * (t21 * t83 + t23 * t80) + m(6) * (-t21 * t84 + t23 * t85) + (-t56 * t80 + t57 * t83) * mrSges(7,2) + (-t56 * t85 - t57 * t84) * mrSges(6,3) - t144 * pkin(8); -0.2e1 * t162 - 0.2e1 * t161 + 0.2e1 * t74 + 0.2e1 * t79 + 0.2e1 * t142 + m(6) * (t84 ^ 2 + t85 ^ 2) + m(7) * (t80 ^ 2 + t83 ^ 2) + m(5) * (t125 ^ 2 + t129 ^ 2) * pkin(3) ^ 2 + t146; t149 * t131 + t110 * t37 + t109 * t34 + m(7) * (t109 * t2 + t110 * t3) + (t124 * t35 + t128 * t36 + m(6) * (t124 * t6 + t128 * t5)) * pkin(4) + t138; (-t109 * t56 + t110 * t57) * mrSges(7,2) + m(7) * (t109 * t23 + t110 * t21) + (m(6) * (t124 * t23 - t128 * t21) + (-t124 * t56 - t128 * t57) * mrSges(6,3)) * pkin(4) + t137; Ifges(5,3) + m(7) * (t109 * t80 + t110 * t83) + t142 + (-t83 - t110) * mrSges(7,1) + (m(6) * (t124 * t85 + t128 * t84) - t155) * pkin(4) + t143 + t168; -0.2e1 * t150 - 0.2e1 * t156 + 0.2e1 * t104 + 0.2e1 * t113 + m(7) * (t109 ^ 2 + t110 ^ 2) + m(6) * (t124 ^ 2 + t128 ^ 2) * pkin(4) ^ 2 - t149; -pkin(5) * t37 + m(7) * (-pkin(5) * t3 + qJ(6) * t2) + qJ(6) * t34 + t159 * t131 + t141; m(7) * (-pkin(5) * t21 + qJ(6) * t23) + (-pkin(5) * t57 - qJ(6) * t56) * mrSges(7,2) + t140; -t162 + m(7) * (-pkin(5) * t83 + qJ(6) * t80) + t143 + t167; -t156 + m(7) * (-pkin(5) * t110 + qJ(6) * t109) - t150 - t159 + t167 + t168; 0.2e1 * t132 + 0.2e1 * t123 + m(7) * (pkin(5) ^ 2 + qJ(6) ^ 2) - t159; m(7) * t3 + t37; m(7) * t21 + t57 * mrSges(7,2); m(7) * t83 - mrSges(7,1); m(7) * t110 - mrSges(7,1); -m(7) * pkin(5) - mrSges(7,1); m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
