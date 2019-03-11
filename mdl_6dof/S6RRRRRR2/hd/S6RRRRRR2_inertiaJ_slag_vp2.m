% Calculate joint inertia matrix for
% S6RRRRRR2
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:19
% EndTime: 2019-03-10 03:33:22
% DurationCPUTime: 1.29s
% Computational Cost: add. (3989->321), mult. (7441->462), div. (0->0), fcn. (8552->10), ass. (0->143)
t140 = cos(qJ(5));
t137 = sin(qJ(3));
t138 = sin(qJ(2));
t142 = cos(qJ(3));
t143 = cos(qJ(2));
t102 = -t137 * t138 + t142 * t143;
t104 = t137 * t143 + t138 * t142;
t136 = sin(qJ(4));
t141 = cos(qJ(4));
t67 = t102 * t136 + t104 * t141;
t176 = t140 * t67;
t66 = -t141 * t102 + t104 * t136;
t201 = Ifges(6,5) * t176 + Ifges(6,3) * t66;
t134 = sin(qJ(6));
t135 = sin(qJ(5));
t174 = Ifges(6,5) * t135 + Ifges(6,6) * t140;
t139 = cos(qJ(6));
t101 = -t134 * t135 + t139 * t140;
t183 = mrSges(7,3) * t101;
t200 = t134 * pkin(5) * t183 + t174;
t121 = pkin(2) * t142 + pkin(3);
t192 = pkin(2) * t137;
t90 = t136 * t121 + t141 * t192;
t88 = pkin(10) + t90;
t74 = (-pkin(11) - t88) * t135;
t129 = t140 * pkin(11);
t175 = t140 * t88;
t75 = t129 + t175;
t45 = t134 * t74 + t139 * t75;
t39 = t45 * t183;
t103 = t134 * t140 + t135 * t139;
t68 = -mrSges(7,1) * t101 + mrSges(7,2) * t103;
t188 = pkin(5) * t140;
t89 = t121 * t141 - t136 * t192;
t87 = -pkin(4) - t89;
t82 = t87 - t188;
t46 = t82 * t68;
t110 = -mrSges(6,1) * t140 + mrSges(6,2) * t135;
t71 = t87 * t110;
t130 = t135 ^ 2;
t185 = mrSges(6,3) * t130;
t80 = t88 * t185;
t132 = t140 ^ 2;
t184 = mrSges(6,3) * t132;
t81 = t88 * t184;
t84 = t89 * mrSges(5,1);
t199 = t39 + t46 + t71 + t80 + t81 + t84;
t191 = pkin(3) * t136;
t119 = pkin(10) + t191;
t107 = t119 * t185;
t108 = t119 * t184;
t190 = pkin(3) * t141;
t126 = mrSges(5,1) * t190;
t97 = (-pkin(11) - t119) * t135;
t98 = t119 * t140 + t129;
t61 = t134 * t97 + t139 * t98;
t47 = t61 * t183;
t122 = -pkin(4) - t188;
t109 = t122 - t190;
t54 = t109 * t68;
t120 = -pkin(4) - t190;
t92 = t120 * t110;
t198 = t107 + t108 + t126 + t47 + t54 + t92;
t193 = -pkin(8) - pkin(7);
t165 = t193 * t138;
t166 = t193 * t143;
t77 = t137 * t166 + t142 * t165;
t151 = -t104 * pkin(9) + t77;
t79 = t137 * t165 - t142 * t166;
t56 = pkin(9) * t102 + t79;
t31 = t136 * t56 - t141 * t151;
t197 = t31 ^ 2;
t196 = 0.2e1 * t31;
t123 = -pkin(2) * t143 - pkin(1);
t83 = -pkin(3) * t102 + t123;
t195 = 0.2e1 * t83;
t189 = pkin(4) * t110;
t187 = t90 * mrSges(5,2);
t30 = pkin(4) * t66 - pkin(10) * t67 + t83;
t33 = t136 * t151 + t141 * t56;
t13 = t135 * t30 + t140 * t33;
t186 = Ifges(7,5) * t103 + Ifges(7,6) * t101;
t182 = Ifges(6,4) * t135;
t181 = Ifges(6,4) * t140;
t180 = t103 * mrSges(7,3);
t12 = -t135 * t33 + t140 * t30;
t179 = t12 * t135;
t178 = t13 * t140;
t177 = t135 * t67;
t173 = t130 + t132;
t172 = t138 ^ 2 + t143 ^ 2;
t171 = -0.2e1 * t180;
t170 = mrSges(5,2) * t191;
t36 = t103 * t67;
t37 = t101 * t67;
t169 = Ifges(7,5) * t37 - Ifges(7,6) * t36 + Ifges(7,3) * t66;
t168 = mrSges(6,3) * t179;
t167 = t139 * t180;
t164 = t173 * t119;
t111 = Ifges(6,2) * t140 + t182;
t112 = Ifges(6,1) * t135 + t181;
t69 = Ifges(7,4) * t103 + Ifges(7,2) * t101;
t70 = Ifges(7,1) * t103 + Ifges(7,4) * t101;
t163 = t101 * t69 + t103 * t70 + t140 * t111 + t135 * t112 + Ifges(5,3);
t162 = mrSges(6,1) * t135 + mrSges(6,2) * t140;
t161 = t178 - t179;
t40 = -mrSges(6,2) * t66 - mrSges(6,3) * t177;
t41 = mrSges(6,1) * t66 - mrSges(6,3) * t176;
t160 = -t135 * t41 + t140 * t40;
t44 = -t134 * t75 + t139 * t74;
t159 = t44 * mrSges(7,1) - t45 * mrSges(7,2) + t186;
t60 = -t134 * t98 + t139 * t97;
t158 = t60 * mrSges(7,1) - t61 * mrSges(7,2) + t186;
t113 = (-pkin(11) - pkin(10)) * t135;
t114 = pkin(10) * t140 + t129;
t76 = t113 * t139 - t114 * t134;
t78 = t113 * t134 + t114 * t139;
t157 = t76 * mrSges(7,1) - t78 * mrSges(7,2) + t186;
t156 = Ifges(4,3) + t163;
t5 = pkin(5) * t66 - pkin(11) * t176 + t12;
t6 = -pkin(11) * t177 + t13;
t3 = -t134 * t6 + t139 * t5;
t4 = t134 * t5 + t139 * t6;
t155 = t3 * mrSges(7,1) - t4 * mrSges(7,2) + t169;
t154 = (mrSges(4,1) * t142 - mrSges(4,2) * t137) * pkin(2);
t153 = (mrSges(7,1) * t139 - mrSges(7,2) * t134) * pkin(5);
t124 = pkin(10) * t185;
t125 = pkin(10) * t184;
t55 = t78 * t183;
t59 = t122 * t68;
t152 = t124 + t125 + t163 + t55 + t59 - t189;
t10 = Ifges(7,1) * t37 - Ifges(7,4) * t36 + Ifges(7,5) * t66;
t16 = pkin(5) * t177 + t31;
t23 = Ifges(6,6) * t66 + (-Ifges(6,2) * t135 + t181) * t67;
t24 = Ifges(6,5) * t66 + (Ifges(6,1) * t140 - t182) * t67;
t9 = Ifges(7,4) * t37 - Ifges(7,2) * t36 + Ifges(7,6) * t66;
t150 = -t33 * mrSges(5,2) - t180 * t3 + t4 * t183 + mrSges(6,3) * t178 + t16 * t68 + t135 * t24 / 0.2e1 + t140 * t23 / 0.2e1 - t36 * t69 / 0.2e1 + t37 * t70 / 0.2e1 - t111 * t177 / 0.2e1 + t112 * t176 / 0.2e1 - Ifges(5,6) * t66 + Ifges(5,5) * t67 + t101 * t9 / 0.2e1 + t103 * t10 / 0.2e1 + (t110 - mrSges(5,1)) * t31 + (t186 + t174) * t66 / 0.2e1;
t149 = t77 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,5) * t104 + Ifges(4,6) * t102 + t150;
t38 = t162 * t67;
t18 = mrSges(7,1) * t66 - mrSges(7,3) * t37;
t17 = -mrSges(7,2) * t66 - mrSges(7,3) * t36;
t14 = mrSges(7,1) * t36 + mrSges(7,2) * t37;
t1 = [0.2e1 * t16 * t14 + 0.2e1 * t4 * t17 + 0.2e1 * t3 * t18 - t36 * t9 + t37 * t10 + t38 * t196 + 0.2e1 * t13 * t40 + 0.2e1 * t12 * t41 + Ifges(2,3) + t104 * (Ifges(4,1) * t104 + Ifges(4,4) * t102) + t102 * (Ifges(4,4) * t104 + Ifges(4,2) * t102) + 0.2e1 * t123 * (-t102 * mrSges(4,1) + t104 * mrSges(4,2)) - 0.2e1 * pkin(1) * (-t143 * mrSges(3,1) + t138 * mrSges(3,2)) + t138 * (Ifges(3,1) * t138 + Ifges(3,4) * t143) + t143 * (Ifges(3,4) * t138 + Ifges(3,2) * t143) + (mrSges(5,2) * t195 + mrSges(5,3) * t196 + Ifges(5,1) * t67 - t135 * t23 + t140 * t24) * t67 + (mrSges(5,1) * t195 - 0.2e1 * t33 * mrSges(5,3) + Ifges(5,2) * t66 + (-Ifges(6,6) * t135 - (2 * Ifges(5,4))) * t67 + t169 + t201) * t66 + m(5) * (t33 ^ 2 + t83 ^ 2 + t197) + m(6) * (t12 ^ 2 + t13 ^ 2 + t197) + m(7) * (t16 ^ 2 + t3 ^ 2 + t4 ^ 2) + m(4) * (t123 ^ 2 + t77 ^ 2 + t79 ^ 2) + m(3) * (pkin(7) ^ 2 * t172 + pkin(1) ^ 2) + 0.2e1 * (t102 * t79 - t104 * t77) * mrSges(4,3) + 0.2e1 * t172 * pkin(7) * mrSges(3,3); (-t12 * mrSges(6,3) - t88 * t41) * t135 + m(5) * (-t31 * t89 + t33 * t90) + m(7) * (t16 * t82 + t3 * t44 + t4 * t45) + m(6) * (t161 * t88 + t87 * t31) + (-mrSges(3,1) * t138 - mrSges(3,2) * t143) * pkin(7) + (-t66 * t90 - t67 * t89) * mrSges(5,3) + t149 + t40 * t175 + t44 * t18 + t45 * t17 + t82 * t14 + t87 * t38 + Ifges(3,5) * t138 + Ifges(3,6) * t143 + (m(4) * (t137 * t79 + t142 * t77) + (t102 * t137 - t104 * t142) * mrSges(4,3)) * pkin(2); m(4) * (t137 ^ 2 + t142 ^ 2) * pkin(2) ^ 2 + m(6) * (t173 * t88 ^ 2 + t87 ^ 2) + t44 * t171 + 0.2e1 * t154 + 0.2e1 * t46 + 0.2e1 * t39 + 0.2e1 * t71 + t156 + m(7) * (t44 ^ 2 + t45 ^ 2 + t82 ^ 2) + m(5) * (t89 ^ 2 + t90 ^ 2) + 0.2e1 * t84 + 0.2e1 * t81 + 0.2e1 * t80 - 0.2e1 * t187 + Ifges(3,3); (m(5) * (t136 * t33 - t141 * t31) + (-t136 * t66 - t141 * t67) * mrSges(5,3)) * pkin(3) + m(6) * (t119 * t161 + t120 * t31) + t160 * t119 + m(7) * (t109 * t16 + t3 * t60 + t4 * t61) + t149 - t168 + t60 * t18 + t61 * t17 + t109 * t14 + t120 * t38; (-t90 - t191) * mrSges(5,2) + t198 + (-t44 - t60) * t180 + m(6) * (t120 * t87 + t164 * t88) + m(7) * (t109 * t82 + t44 * t60 + t45 * t61) + t156 + t154 + m(5) * (t136 * t90 + t141 * t89) * pkin(3) + t199; -0.2e1 * t170 + t60 * t171 + 0.2e1 * t107 + 0.2e1 * t108 + 0.2e1 * t126 + 0.2e1 * t47 + 0.2e1 * t54 + 0.2e1 * t92 + m(7) * (t109 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t119 ^ 2 * t173 + t120 ^ 2) + m(5) * (t136 ^ 2 + t141 ^ 2) * pkin(3) ^ 2 + t156; m(6) * (-pkin(4) * t31 + pkin(10) * t161) + t160 * pkin(10) + m(7) * (t122 * t16 + t3 * t76 + t4 * t78) + t150 - t168 - pkin(4) * t38 + t76 * t18 + t78 * t17 + t122 * t14; m(6) * (pkin(10) * t173 * t88 - pkin(4) * t87) + m(7) * (t122 * t82 + t44 * t76 + t45 * t78) + t152 - t187 + (-t44 - t76) * t180 + t199; m(6) * (-pkin(4) * t120 + pkin(10) * t164) + m(7) * (t109 * t122 + t60 * t76 + t61 * t78) + t152 - t170 + (-t60 - t76) * t180 + t198; t76 * t171 - 0.2e1 * t189 + 0.2e1 * t124 + 0.2e1 * t125 + 0.2e1 * t55 + 0.2e1 * t59 + m(7) * (t122 ^ 2 + t76 ^ 2 + t78 ^ 2) + m(6) * (pkin(10) ^ 2 * t173 + pkin(4) ^ 2) + t163; -Ifges(6,6) * t177 + t12 * mrSges(6,1) - t13 * mrSges(6,2) + (m(7) * (t134 * t4 + t139 * t3) + t134 * t17 + t139 * t18) * pkin(5) + t155 + t201; -t162 * t88 + (m(7) * (t134 * t45 + t139 * t44) - t167) * pkin(5) + t159 + t200; -t162 * t119 + (m(7) * (t134 * t61 + t139 * t60) - t167) * pkin(5) + t158 + t200; -t162 * pkin(10) + (m(7) * (t134 * t78 + t139 * t76) - t167) * pkin(5) + t157 + t200; Ifges(6,3) + Ifges(7,3) + m(7) * (t134 ^ 2 + t139 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t153; t155; t159; t158; t157; Ifges(7,3) + t153; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
