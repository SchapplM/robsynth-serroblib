% Calculate joint inertia matrix for
% S6PRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:56
% EndTime: 2019-03-09 00:55:00
% DurationCPUTime: 1.72s
% Computational Cost: add. (2941->367), mult. (7009->543), div. (0->0), fcn. (7921->14), ass. (0->142)
t146 = sin(qJ(6));
t151 = cos(qJ(6));
t113 = -mrSges(7,1) * t151 + mrSges(7,2) * t146;
t208 = -mrSges(6,1) + t113;
t147 = sin(qJ(5));
t148 = sin(qJ(4));
t152 = cos(qJ(5));
t153 = cos(qJ(4));
t109 = t147 * t148 - t152 * t153;
t110 = t147 * t153 + t148 * t152;
t131 = -pkin(4) * t153 - pkin(3);
t67 = pkin(5) * t109 - pkin(12) * t110 + t131;
t202 = -pkin(11) - pkin(10);
t121 = t202 * t153;
t171 = t202 * t148;
t84 = -t152 * t121 + t147 * t171;
t34 = -t146 * t84 + t151 * t67;
t35 = t146 * t67 + t151 * t84;
t165 = -t146 * t34 + t151 * t35;
t144 = cos(pkin(7));
t142 = sin(pkin(7));
t149 = sin(qJ(3));
t179 = t142 * t149;
t97 = t144 * t153 - t148 * t179;
t98 = t144 * t148 + t153 * t179;
t207 = -Ifges(5,5) * t98 - Ifges(5,6) * t97;
t154 = cos(qJ(3));
t178 = t142 * t154;
t143 = sin(pkin(6));
t145 = cos(pkin(6));
t150 = sin(qJ(2));
t155 = cos(qJ(2));
t177 = t144 * t155;
t73 = t145 * t179 + (t149 * t177 + t150 * t154) * t143;
t96 = -t142 * t143 * t155 + t144 * t145;
t40 = -t148 * t73 + t153 * t96;
t41 = t148 * t96 + t153 * t73;
t21 = t147 * t41 - t152 * t40;
t206 = t21 ^ 2;
t71 = -t145 * t178 + (t149 * t150 - t154 * t177) * t143;
t70 = t71 ^ 2;
t82 = -t121 * t147 - t152 * t171;
t205 = t82 ^ 2;
t204 = 0.2e1 * t82;
t62 = t147 * t97 + t152 * t98;
t46 = -t146 * t62 - t151 * t178;
t203 = t46 / 0.2e1;
t115 = Ifges(7,5) * t146 + Ifges(7,6) * t151;
t201 = t115 / 0.2e1;
t187 = Ifges(7,4) * t151;
t118 = Ifges(7,1) * t146 + t187;
t200 = t118 / 0.2e1;
t199 = t146 / 0.2e1;
t198 = t151 / 0.2e1;
t196 = pkin(2) * t154;
t195 = pkin(12) * t146;
t194 = pkin(12) * t151;
t193 = t21 * t82;
t192 = Ifges(5,3) + Ifges(6,3);
t47 = -t146 * t178 + t151 * t62;
t24 = -mrSges(7,1) * t46 + mrSges(7,2) * t47;
t53 = -mrSges(6,1) * t178 - mrSges(6,3) * t62;
t191 = t24 - t53;
t101 = t144 * t149 * pkin(2) + pkin(9) * t178;
t90 = pkin(10) * t144 + t101;
t91 = (-pkin(3) * t154 - pkin(10) * t149 - pkin(2)) * t142;
t54 = -t148 * t90 + t153 * t91;
t32 = -pkin(4) * t178 - pkin(11) * t98 + t54;
t55 = t148 * t91 + t153 * t90;
t36 = pkin(11) * t97 + t55;
t17 = t147 * t32 + t152 * t36;
t61 = t147 * t98 - t152 * t97;
t190 = -Ifges(6,5) * t62 + Ifges(6,6) * t61;
t189 = mrSges(7,3) * t151;
t188 = Ifges(7,4) * t146;
t186 = t146 * mrSges(7,3);
t183 = t110 * t146;
t182 = t110 * t151;
t129 = pkin(4) * t147 + pkin(12);
t181 = t129 * t146;
t180 = t129 * t151;
t176 = Ifges(6,5) * t110 - Ifges(6,6) * t109;
t175 = Ifges(5,5) * t148 + Ifges(5,6) * t153;
t174 = t146 ^ 2 + t151 ^ 2;
t173 = t148 ^ 2 + t153 ^ 2;
t12 = Ifges(7,5) * t47 + Ifges(7,6) * t46 + Ifges(7,3) * t61;
t116 = Ifges(7,2) * t151 + t188;
t172 = t151 * t116 + t146 * t118 + Ifges(6,3);
t170 = Ifges(4,5) * t179 + Ifges(4,6) * t178 + Ifges(4,3) * t144;
t169 = t174 * t129;
t11 = -pkin(12) * t178 + t17;
t124 = pkin(9) * t179;
t89 = t124 + (-pkin(3) - t196) * t144;
t63 = -pkin(4) * t97 + t89;
t19 = pkin(5) * t61 - pkin(12) * t62 + t63;
t2 = -t11 * t146 + t151 * t19;
t3 = t11 * t151 + t146 * t19;
t168 = -t146 * t2 + t151 * t3;
t23 = t147 * t40 + t152 * t41;
t6 = -t146 * t23 + t151 * t71;
t7 = t146 * t71 + t151 * t23;
t167 = -t146 * t6 + t151 * t7;
t166 = mrSges(7,1) * t146 + mrSges(7,2) * t151;
t16 = -t147 * t36 + t152 * t32;
t164 = -t40 * t148 + t41 * t153;
t163 = 0.2e1 * t174 * mrSges(7,3);
t49 = Ifges(7,5) * t182 - Ifges(7,6) * t183 + Ifges(7,3) * t109;
t162 = (t152 * mrSges(6,1) - t147 * mrSges(6,2)) * pkin(4);
t161 = -t23 * mrSges(6,2) - t6 * t186 + t7 * t189 + t208 * t21;
t10 = pkin(5) * t178 - t16;
t13 = Ifges(7,4) * t47 + Ifges(7,2) * t46 + Ifges(7,6) * t61;
t14 = Ifges(7,1) * t47 + Ifges(7,4) * t46 + Ifges(7,5) * t61;
t160 = t16 * mrSges(6,1) - t17 * mrSges(6,2) + t10 * t113 + t116 * t203 + t13 * t198 + t14 * t199 - t2 * t186 + t3 * t189 + t47 * t200 + t61 * t201 - t190;
t50 = Ifges(7,6) * t109 + (-Ifges(7,2) * t146 + t187) * t110;
t51 = Ifges(7,5) * t109 + (Ifges(7,1) * t151 - t188) * t110;
t159 = -t84 * mrSges(6,2) + t50 * t198 + t51 * t199 + t176 - t116 * t183 / 0.2e1 + t182 * t200 + t109 * t201 + t208 * t82 + t165 * mrSges(7,3);
t130 = -pkin(4) * t152 - pkin(5);
t119 = Ifges(5,1) * t148 + Ifges(5,4) * t153;
t117 = Ifges(5,4) * t148 + Ifges(5,2) * t153;
t114 = -mrSges(5,1) * t153 + mrSges(5,2) * t148;
t106 = -mrSges(4,2) * t144 + mrSges(4,3) * t178;
t105 = mrSges(4,1) * t144 - mrSges(4,3) * t179;
t100 = t144 * t196 - t124;
t99 = (-mrSges(4,1) * t154 + mrSges(4,2) * t149) * t142;
t81 = -mrSges(5,1) * t178 - mrSges(5,3) * t98;
t80 = mrSges(5,2) * t178 + mrSges(5,3) * t97;
t76 = Ifges(6,1) * t110 - Ifges(6,4) * t109;
t75 = Ifges(6,4) * t110 - Ifges(6,2) * t109;
t74 = mrSges(6,1) * t109 + mrSges(6,2) * t110;
t69 = mrSges(7,1) * t109 - mrSges(7,3) * t182;
t68 = -mrSges(7,2) * t109 - mrSges(7,3) * t183;
t66 = t166 * t110;
t64 = -mrSges(5,1) * t97 + mrSges(5,2) * t98;
t57 = Ifges(5,1) * t98 + Ifges(5,4) * t97 - Ifges(5,5) * t178;
t56 = Ifges(5,4) * t98 + Ifges(5,2) * t97 - Ifges(5,6) * t178;
t52 = mrSges(6,2) * t178 - mrSges(6,3) * t61;
t29 = mrSges(6,1) * t61 + mrSges(6,2) * t62;
t28 = Ifges(6,1) * t62 - Ifges(6,4) * t61 - Ifges(6,5) * t178;
t27 = Ifges(6,4) * t62 - Ifges(6,2) * t61 - Ifges(6,6) * t178;
t26 = mrSges(7,1) * t61 - mrSges(7,3) * t47;
t25 = -mrSges(7,2) * t61 + mrSges(7,3) * t46;
t1 = [m(2) + m(7) * (t6 ^ 2 + t7 ^ 2 + t206) + m(6) * (t23 ^ 2 + t206 + t70) + m(5) * (t40 ^ 2 + t41 ^ 2 + t70) + m(4) * (t73 ^ 2 + t96 ^ 2 + t70) + m(3) * (t145 ^ 2 + (t150 ^ 2 + t155 ^ 2) * t143 ^ 2); t73 * t106 + t23 * t52 + t7 * t25 + t6 * t26 + t40 * t81 + t41 * t80 + t96 * t99 + t191 * t21 + (t155 * mrSges(3,1) - t150 * mrSges(3,2)) * t143 + (-t105 + t29 + t64) * t71 + m(7) * (t10 * t21 + t2 * t6 + t3 * t7) + m(6) * (-t16 * t21 + t17 * t23 + t63 * t71) + m(5) * (t40 * t54 + t41 * t55 + t71 * t89) + m(4) * (-pkin(2) * t142 * t96 - t100 * t71 + t101 * t73); t144 * t170 + t97 * t56 + t98 * t57 + 0.2e1 * t100 * t105 + 0.2e1 * t101 * t106 + 0.2e1 * t89 * t64 + 0.2e1 * t55 * t80 + 0.2e1 * t54 * t81 + t62 * t28 + 0.2e1 * t63 * t29 + t47 * t14 + 0.2e1 * t17 * t52 + 0.2e1 * t16 * t53 + t46 * t13 + 0.2e1 * t3 * t25 + 0.2e1 * t2 * t26 + 0.2e1 * t10 * t24 + Ifges(3,3) + (t12 - t27) * t61 + (-0.2e1 * pkin(2) * t99 + (Ifges(4,1) * t179 + Ifges(4,5) * t144) * t149 + (Ifges(4,6) * t144 + (0.2e1 * Ifges(4,4) * t149 + (Ifges(4,2) + t192) * t154) * t142 + t190 + t207) * t154) * t142 + m(4) * (pkin(2) ^ 2 * t142 ^ 2 + t100 ^ 2 + t101 ^ 2) + m(5) * (t54 ^ 2 + t55 ^ 2 + t89 ^ 2) + m(6) * (t16 ^ 2 + t17 ^ 2 + t63 ^ 2) + m(7) * (t10 ^ 2 + t2 ^ 2 + t3 ^ 2); -t73 * mrSges(4,2) + t21 * t66 + t6 * t69 + t7 * t68 + (-t23 * t109 + t21 * t110) * mrSges(6,3) + t164 * mrSges(5,3) + (-mrSges(4,1) + t114 + t74) * t71 + m(7) * (t34 * t6 + t35 * t7 + t193) + m(6) * (t131 * t71 + t23 * t84 + t193) + m(5) * (-pkin(3) * t71 + t164 * pkin(10)); t131 * t29 + t89 * t114 + t97 * t117 / 0.2e1 + t98 * t119 / 0.2e1 + t100 * mrSges(4,1) - t101 * mrSges(4,2) + t84 * t52 + m(6) * (t131 * t63 - t16 * t82 + t17 * t84) + m(7) * (t10 * t82 + t2 * t34 + t3 * t35) + (-t75 / 0.2e1 + t49 / 0.2e1) * t61 + t170 + m(5) * (-pkin(3) * t89 + (-t54 * t148 + t55 * t153) * pkin(10)) + t2 * t69 + t63 * t74 + t62 * t76 / 0.2e1 - pkin(3) * t64 + t10 * t66 + t3 * t68 + t47 * t51 / 0.2e1 + t34 * t26 + t35 * t25 + (t14 * t198 - t146 * t13 / 0.2e1 + t28 / 0.2e1 - t16 * mrSges(6,3)) * t110 + t191 * t82 + (-t27 / 0.2e1 + t12 / 0.2e1 - t17 * mrSges(6,3)) * t109 + (pkin(10) * t80 + t55 * mrSges(5,3) + t56 / 0.2e1) * t153 + (-pkin(10) * t81 - t54 * mrSges(5,3) + t57 / 0.2e1) * t148 - (t175 + t176) * t178 / 0.2e1 + t50 * t203; -0.2e1 * pkin(3) * t114 + t153 * t117 + t148 * t119 + 0.2e1 * t131 * t74 + 0.2e1 * t34 * t69 + 0.2e1 * t35 * t68 + t66 * t204 + Ifges(4,3) + 0.2e1 * t173 * pkin(10) * mrSges(5,3) + (-0.2e1 * t84 * mrSges(6,3) + t49 - t75) * t109 + m(7) * (t34 ^ 2 + t35 ^ 2 + t205) + m(6) * (t131 ^ 2 + t84 ^ 2 + t205) + m(5) * (t173 * pkin(10) ^ 2 + pkin(3) ^ 2) + (mrSges(6,3) * t204 - t146 * t50 + t151 * t51 + t76) * t110; t40 * mrSges(5,1) - t41 * mrSges(5,2) + m(7) * (t167 * t129 + t130 * t21) + m(6) * (t147 * t23 - t152 * t21) * pkin(4) + t161; t25 * t180 - t26 * t181 + t130 * t24 + m(7) * (t10 * t130 + t168 * t129) + (m(6) * (t147 * t17 + t152 * t16) + t152 * t53 + t147 * t52) * pkin(4) + t54 * mrSges(5,1) - t55 * mrSges(5,2) - t192 * t178 + t160 - t207; t68 * t180 - t69 * t181 + t130 * t66 + m(7) * (t165 * t129 + t130 * t82) + (m(6) * (t147 * t84 - t152 * t82) + (-t147 * t109 - t152 * t110) * mrSges(6,3)) * pkin(4) + t159 + (-t148 * mrSges(5,1) - t153 * mrSges(5,2)) * pkin(10) + t175; 0.2e1 * t130 * t113 + Ifges(5,3) + 0.2e1 * t162 + t129 * t163 + m(7) * (t174 * t129 ^ 2 + t130 ^ 2) + m(6) * (t147 ^ 2 + t152 ^ 2) * pkin(4) ^ 2 + t172; m(7) * (-pkin(5) * t21 + t167 * pkin(12)) + t161; -Ifges(6,3) * t178 + t25 * t194 - t26 * t195 + m(7) * (-pkin(5) * t10 + t168 * pkin(12)) - pkin(5) * t24 + t160; t68 * t194 - t69 * t195 + m(7) * (-pkin(5) * t82 + t165 * pkin(12)) - pkin(5) * t66 + t159; m(7) * (-pkin(5) * t130 + pkin(12) * t169) + (t130 - pkin(5)) * t113 + t162 + (t174 * pkin(12) + t169) * mrSges(7,3) + t172; -0.2e1 * pkin(5) * t113 + m(7) * (t174 * pkin(12) ^ 2 + pkin(5) ^ 2) + pkin(12) * t163 + t172; mrSges(7,1) * t6 - mrSges(7,2) * t7; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t34 - mrSges(7,2) * t35 + t49; -t166 * t129 + t115; -t166 * pkin(12) + t115; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
