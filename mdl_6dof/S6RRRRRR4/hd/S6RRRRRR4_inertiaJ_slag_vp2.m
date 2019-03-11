% Calculate joint inertia matrix for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:47
% EndTime: 2019-03-10 03:46:51
% DurationCPUTime: 1.59s
% Computational Cost: add. (4418->375), mult. (8847->543), div. (0->0), fcn. (9684->10), ass. (0->152)
t185 = 2 * pkin(7);
t137 = sin(qJ(5));
t142 = cos(qJ(5));
t138 = sin(qJ(4));
t139 = sin(qJ(3));
t143 = cos(qJ(4));
t144 = cos(qJ(3));
t106 = t138 * t144 + t139 * t143;
t140 = sin(qJ(2));
t93 = t106 * t140;
t105 = -t138 * t139 + t143 * t144;
t94 = t105 * t140;
t56 = -t137 * t94 - t142 * t93;
t57 = -t137 * t93 + t142 * t94;
t184 = -Ifges(6,5) * t57 - Ifges(6,6) * t56;
t183 = -Ifges(5,5) * t94 + Ifges(5,6) * t93;
t182 = -pkin(9) - pkin(8);
t181 = pkin(3) * t138;
t180 = pkin(4) * t137;
t179 = pkin(4) * t142;
t145 = cos(qJ(2));
t178 = pkin(7) * t145;
t131 = t140 * pkin(7);
t136 = sin(qJ(6));
t141 = cos(qJ(6));
t124 = pkin(3) * t143 + pkin(4);
t97 = t142 * t124 - t137 * t181;
t95 = pkin(5) + t97;
t99 = t124 * t137 + t142 * t181;
t62 = t136 * t95 + t141 * t99;
t177 = t62 * mrSges(7,2);
t123 = pkin(5) + t179;
t98 = t123 * t136 + t141 * t180;
t176 = t98 * mrSges(7,2);
t175 = t99 * mrSges(6,2);
t174 = -Ifges(7,3) - Ifges(6,3);
t61 = -t136 * t99 + t141 * t95;
t55 = t61 * mrSges(7,1);
t173 = Ifges(7,3) + t55;
t28 = -t136 * t57 + t141 * t56;
t29 = t136 * t56 + t141 * t57;
t172 = -Ifges(7,5) * t29 - Ifges(7,6) * t28;
t113 = -pkin(2) * t145 - pkin(8) * t140 - pkin(1);
t104 = t144 * t113;
t166 = t140 * t144;
t72 = -pkin(9) * t166 + t104 + (-pkin(7) * t139 - pkin(3)) * t145;
t167 = t139 * t140;
t85 = t139 * t113 + t144 * t178;
t78 = -pkin(9) * t167 + t85;
t46 = -t138 * t78 + t143 * t72;
t30 = -pkin(4) * t145 - pkin(10) * t94 + t46;
t47 = t138 * t72 + t143 * t78;
t36 = -pkin(10) * t93 + t47;
t14 = t137 * t30 + t142 * t36;
t117 = t182 * t139;
t118 = t182 * t144;
t80 = t143 * t117 + t118 * t138;
t64 = -pkin(10) * t106 + t80;
t81 = t138 * t117 - t143 * t118;
t65 = pkin(10) * t105 + t81;
t35 = t137 * t64 + t142 * t65;
t171 = Ifges(4,4) * t139;
t170 = Ifges(4,4) * t144;
t169 = t136 * mrSges(7,2);
t168 = t137 * mrSges(6,2);
t112 = pkin(3) * t167 + t131;
t165 = t139 ^ 2 + t144 ^ 2;
t164 = pkin(5) * t169;
t163 = pkin(4) * t168;
t162 = -Ifges(5,3) + t174;
t125 = -pkin(3) * t144 - pkin(2);
t13 = -t137 * t36 + t142 * t30;
t34 = -t137 * t65 + t142 * t64;
t161 = Ifges(4,3) - t162;
t96 = t123 * t141 - t136 * t180;
t91 = t96 * mrSges(7,1);
t160 = Ifges(7,3) + t91 - t176;
t73 = pkin(4) * t93 + t112;
t7 = -pkin(5) * t145 - pkin(11) * t57 + t13;
t8 = pkin(11) * t56 + t14;
t2 = -t136 * t8 + t141 * t7;
t3 = t136 * t7 + t141 * t8;
t159 = t2 * mrSges(7,1) - t3 * mrSges(7,2) - t172;
t70 = t105 * t142 - t106 * t137;
t71 = t105 * t137 + t106 * t142;
t40 = -t136 * t71 + t141 * t70;
t38 = Ifges(7,6) * t40;
t41 = t136 * t70 + t141 * t71;
t39 = Ifges(7,5) * t41;
t18 = -pkin(11) * t71 + t34;
t19 = pkin(11) * t70 + t35;
t5 = -t136 * t19 + t141 * t18;
t6 = t136 * t18 + t141 * t19;
t158 = t5 * mrSges(7,1) - t6 * mrSges(7,2) + t38 + t39;
t157 = mrSges(4,1) * t139 + mrSges(4,2) * t144;
t86 = -pkin(4) * t105 + t125;
t92 = t97 * mrSges(6,1);
t156 = Ifges(6,3) + t173 + t92 - t175;
t155 = (mrSges(5,1) * t143 - mrSges(5,2) * t138) * pkin(3);
t154 = t13 * mrSges(6,1) - t14 * mrSges(6,2) + t159 - t184;
t68 = Ifges(6,6) * t70;
t69 = Ifges(6,5) * t71;
t153 = t34 * mrSges(6,1) - t35 * mrSges(6,2) + t158 + t68 + t69;
t152 = t46 * mrSges(5,1) - t47 * mrSges(5,2) + t154 - t183;
t101 = Ifges(5,6) * t105;
t102 = Ifges(5,5) * t106;
t151 = t80 * mrSges(5,1) - t81 * mrSges(5,2) + t101 + t102 + t153;
t147 = pkin(7) ^ 2;
t135 = t145 ^ 2;
t133 = t140 ^ 2;
t130 = t133 * t147;
t129 = Ifges(4,5) * t139;
t128 = Ifges(4,6) * t144;
t127 = mrSges(6,1) * t179;
t126 = t141 * pkin(5) * mrSges(7,1);
t120 = Ifges(4,5) * t166;
t116 = Ifges(4,1) * t139 + t170;
t115 = Ifges(4,2) * t144 + t171;
t114 = -mrSges(4,1) * t144 + mrSges(4,2) * t139;
t111 = -mrSges(4,1) * t145 - mrSges(4,3) * t166;
t110 = mrSges(4,2) * t145 - mrSges(4,3) * t167;
t100 = t157 * t140;
t90 = -Ifges(4,5) * t145 + (Ifges(4,1) * t144 - t171) * t140;
t89 = -Ifges(4,6) * t145 + (-Ifges(4,2) * t139 + t170) * t140;
t84 = -t139 * t178 + t104;
t83 = -mrSges(5,1) * t145 - mrSges(5,3) * t94;
t82 = mrSges(5,2) * t145 - mrSges(5,3) * t93;
t77 = Ifges(5,1) * t106 + Ifges(5,4) * t105;
t76 = Ifges(5,4) * t106 + Ifges(5,2) * t105;
t75 = -mrSges(5,1) * t105 + mrSges(5,2) * t106;
t63 = mrSges(5,1) * t93 + mrSges(5,2) * t94;
t54 = Ifges(5,1) * t94 - Ifges(5,4) * t93 - Ifges(5,5) * t145;
t53 = Ifges(5,4) * t94 - Ifges(5,2) * t93 - Ifges(5,6) * t145;
t50 = -mrSges(6,1) * t145 - mrSges(6,3) * t57;
t49 = mrSges(6,2) * t145 + mrSges(6,3) * t56;
t48 = -pkin(5) * t70 + t86;
t44 = Ifges(6,1) * t71 + Ifges(6,4) * t70;
t43 = Ifges(6,4) * t71 + Ifges(6,2) * t70;
t42 = -mrSges(6,1) * t70 + mrSges(6,2) * t71;
t37 = -pkin(5) * t56 + t73;
t31 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t25 = Ifges(6,1) * t57 + Ifges(6,4) * t56 - Ifges(6,5) * t145;
t24 = Ifges(6,4) * t57 + Ifges(6,2) * t56 - Ifges(6,6) * t145;
t21 = -mrSges(7,1) * t145 - mrSges(7,3) * t29;
t20 = mrSges(7,2) * t145 + mrSges(7,3) * t28;
t17 = Ifges(7,1) * t41 + Ifges(7,4) * t40;
t16 = Ifges(7,4) * t41 + Ifges(7,2) * t40;
t15 = -mrSges(7,1) * t40 + mrSges(7,2) * t41;
t11 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t28 - Ifges(7,5) * t145;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t28 - Ifges(7,6) * t145;
t1 = [(-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t140 + t100 * t185 - t139 * t89 + t144 * t90) * t140 + (0.2e1 * pkin(1) * mrSges(3,1) - t120 + (Ifges(3,2) + t161) * t145 + (Ifges(4,6) * t139 + (2 * Ifges(3,4))) * t140 + t172 + t183 + t184) * t145 + m(3) * (pkin(1) ^ 2 + t135 * t147 + t130) + m(4) * (t84 ^ 2 + t85 ^ 2 + t130) + m(5) * (t112 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t73 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t37 ^ 2) + t28 * t9 + t29 * t10 + 0.2e1 * t73 * t31 + 0.2e1 * t47 * t82 + 0.2e1 * t46 * t83 - t93 * t53 + t94 * t54 + 0.2e1 * t85 * t110 + 0.2e1 * t84 * t111 + 0.2e1 * t112 * t63 + 0.2e1 * t37 * t11 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + Ifges(2,3) + 0.2e1 * t14 * t49 + 0.2e1 * t13 * t50 + t56 * t24 + t57 * t25 + (t133 + t135) * mrSges(3,3) * t185; m(5) * (t112 * t125 + t46 * t80 + t47 * t81) + m(6) * (t13 * t34 + t14 * t35 + t73 * t86) + m(7) * (t2 * t5 + t6 * t3 + t48 * t37) + (-t2 * t41 + t3 * t40) * mrSges(7,3) + (-t13 * t71 + t14 * t70) * mrSges(6,3) + (t47 * t105 - t46 * t106) * mrSges(5,3) + (-pkin(7) * mrSges(3,2) - t69 / 0.2e1 - t68 / 0.2e1 - t102 / 0.2e1 - t101 / 0.2e1 - t39 / 0.2e1 - t38 / 0.2e1 - t129 / 0.2e1 - t128 / 0.2e1 + Ifges(3,6)) * t145 + (t144 * t116 / 0.2e1 - t139 * t115 / 0.2e1 + Ifges(3,5) + (t114 - mrSges(3,1)) * pkin(7)) * t140 + t28 * t16 / 0.2e1 + t29 * t17 / 0.2e1 + t70 * t24 / 0.2e1 + t71 * t25 / 0.2e1 + t73 * t42 + t81 * t82 + t80 * t83 + t86 * t31 - t93 * t76 / 0.2e1 + t94 * t77 / 0.2e1 - pkin(2) * t100 + t105 * t53 / 0.2e1 + t106 * t54 / 0.2e1 + t112 * t75 + t37 * t15 + t40 * t9 / 0.2e1 + t41 * t10 / 0.2e1 + t6 * t20 + t5 * t21 + t125 * t63 + m(4) * (-pkin(2) * t131 + (-t84 * t139 + t85 * t144) * pkin(8)) + t48 * t11 + t35 * t49 + t34 * t50 + t56 * t43 / 0.2e1 + t57 * t44 / 0.2e1 + (-pkin(8) * t111 - t84 * mrSges(4,3) + t90 / 0.2e1) * t139 + (pkin(8) * t110 + t85 * mrSges(4,3) + t89 / 0.2e1) * t144; -0.2e1 * pkin(2) * t114 + t105 * t76 + t106 * t77 + t144 * t115 + t139 * t116 + 0.2e1 * t125 * t75 + 0.2e1 * t48 * t15 + t40 * t16 + t41 * t17 + 0.2e1 * t86 * t42 + t70 * t43 + t71 * t44 + Ifges(3,3) + m(4) * (t165 * pkin(8) ^ 2 + pkin(2) ^ 2) + m(5) * (t125 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(6) * (t34 ^ 2 + t35 ^ 2 + t86 ^ 2) + m(7) * (t48 ^ 2 + t5 ^ 2 + t6 ^ 2) + 0.2e1 * (t40 * t6 - t41 * t5) * mrSges(7,3) + 0.2e1 * (-t34 * t71 + t35 * t70) * mrSges(6,3) + 0.2e1 * (t105 * t81 - t106 * t80) * mrSges(5,3) + 0.2e1 * t165 * pkin(8) * mrSges(4,3); m(6) * (t13 * t97 + t14 * t99) + m(7) * (t2 * t61 + t3 * t62) - t161 * t145 - Ifges(4,6) * t167 + (m(5) * (t138 * t47 + t143 * t46) + t143 * t83 + t138 * t82) * pkin(3) + t120 + t152 + t84 * mrSges(4,1) - t85 * mrSges(4,2) + t97 * t50 + t99 * t49 + t61 * t21 + t62 * t20; m(6) * (t34 * t97 + t35 * t99) + m(7) * (t5 * t61 + t6 * t62) - t157 * pkin(8) + (t40 * t62 - t41 * t61) * mrSges(7,3) + (t70 * t99 - t71 * t97) * mrSges(6,3) + (m(5) * (t138 * t81 + t143 * t80) + (t105 * t138 - t106 * t143) * mrSges(5,3)) * pkin(3) + t128 + t129 + t151; -0.2e1 * t175 - 0.2e1 * t177 + 0.2e1 * t55 + 0.2e1 * t92 + 0.2e1 * t155 + m(7) * (t61 ^ 2 + t62 ^ 2) + m(6) * (t97 ^ 2 + t99 ^ 2) + m(5) * (t138 ^ 2 + t143 ^ 2) * pkin(3) ^ 2 + t161; (t142 * t50 + t137 * t49 + m(6) * (t13 * t142 + t137 * t14)) * pkin(4) + m(7) * (t2 * t96 + t3 * t98) + t162 * t145 + t152 + t96 * t21 + t98 * t20; (m(6) * (t137 * t35 + t142 * t34) + (t137 * t70 - t142 * t71) * mrSges(6,3)) * pkin(4) + m(7) * (t5 * t96 + t6 * t98) + (t40 * t98 - t41 * t96) * mrSges(7,3) + t151; Ifges(5,3) + m(7) * (t61 * t96 + t62 * t98) + t91 + t127 + t155 + (-t62 - t98) * mrSges(7,2) + (m(6) * (t137 * t99 + t142 * t97) - t168) * pkin(4) + t156; -0.2e1 * t163 - 0.2e1 * t176 + 0.2e1 * t127 + 0.2e1 * t91 + m(7) * (t96 ^ 2 + t98 ^ 2) + m(6) * (t137 ^ 2 + t142 ^ 2) * pkin(4) ^ 2 - t162; t174 * t145 + (t141 * t21 + m(7) * (t136 * t3 + t141 * t2) + t136 * t20) * pkin(5) + t154; (m(7) * (t136 * t6 + t141 * t5) + (t136 * t40 - t141 * t41) * mrSges(7,3)) * pkin(5) + t153; -t177 + t126 + (m(7) * (t136 * t62 + t141 * t61) - t169) * pkin(5) + t156; -t163 + Ifges(6,3) + t126 + t127 + (-t169 + m(7) * (t136 * t98 + t141 * t96)) * pkin(5) + t160; -0.2e1 * t164 + 0.2e1 * t126 + m(7) * (t136 ^ 2 + t141 ^ 2) * pkin(5) ^ 2 - t174; -Ifges(7,3) * t145 + t159; t158; t173 - t177; t160; Ifges(7,3) + t126 - t164; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
