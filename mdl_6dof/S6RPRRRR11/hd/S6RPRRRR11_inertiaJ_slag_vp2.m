% Calculate joint inertia matrix for
% S6RPRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:22
% EndTime: 2019-03-09 07:38:27
% DurationCPUTime: 2.05s
% Computational Cost: add. (6588->434), mult. (17353->641), div. (0->0), fcn. (19944->14), ass. (0->170)
t219 = 2 * pkin(10);
t156 = sin(pkin(7));
t159 = cos(pkin(7));
t164 = sin(qJ(3));
t168 = cos(qJ(3));
t155 = sin(pkin(13));
t157 = sin(pkin(6));
t158 = cos(pkin(13));
t190 = t157 * t158;
t160 = cos(pkin(6));
t205 = pkin(1) * t160;
t110 = qJ(2) * t190 + t155 * t205;
t182 = t159 * t190;
t73 = (t156 * t160 + t182) * pkin(9) + t110;
t139 = t158 * t205;
t193 = t155 * t157;
t80 = pkin(2) * t160 + t139 + (-pkin(9) * t159 - qJ(2)) * t193;
t95 = (-pkin(9) * t155 * t156 - pkin(2) * t158 - pkin(1)) * t157;
t42 = -t164 * t73 + (t156 * t95 + t159 * t80) * t168;
t218 = 2 * mrSges(3,1);
t163 = sin(qJ(4));
t167 = cos(qJ(4));
t192 = t156 * t164;
t111 = -t167 * t159 + t163 * t192;
t108 = t111 ^ 2;
t217 = 0.2e1 * t160;
t162 = sin(qJ(5));
t166 = cos(qJ(5));
t107 = -t156 * t190 + t159 * t160;
t189 = t159 * t164;
t79 = t160 * t192 + (t155 * t168 + t158 * t189) * t157;
t63 = t107 * t163 + t167 * t79;
t191 = t156 * t168;
t78 = -t160 * t191 + t164 * t193 - t168 * t182;
t46 = -t162 * t63 + t166 * t78;
t47 = t162 * t78 + t166 * t63;
t62 = -t107 * t167 + t163 * t79;
t22 = Ifges(6,1) * t47 + Ifges(6,4) * t46 + Ifges(6,5) * t62;
t216 = t22 / 0.2e1;
t161 = sin(qJ(6));
t165 = cos(qJ(6));
t28 = -t161 * t47 + t165 * t46;
t215 = t28 / 0.2e1;
t29 = t161 * t46 + t165 * t47;
t214 = t29 / 0.2e1;
t213 = -pkin(12) - pkin(11);
t199 = Ifges(6,4) * t166;
t102 = -Ifges(6,6) * t167 + (-Ifges(6,2) * t162 + t199) * t163;
t212 = t102 / 0.2e1;
t200 = Ifges(6,4) * t162;
t103 = -Ifges(6,5) * t167 + (Ifges(6,1) * t166 - t200) * t163;
t211 = t103 / 0.2e1;
t120 = t161 * t166 + t162 * t165;
t105 = t120 * t163;
t210 = -t105 / 0.2e1;
t119 = -t161 * t162 + t165 * t166;
t106 = t119 * t163;
t209 = t106 / 0.2e1;
t208 = t119 / 0.2e1;
t207 = t120 / 0.2e1;
t131 = Ifges(6,1) * t162 + t199;
t206 = t131 / 0.2e1;
t204 = pkin(10) * t163;
t203 = pkin(10) * t167;
t202 = -Ifges(7,3) - Ifges(6,3);
t57 = -t156 * t80 + t159 * t95;
t36 = pkin(3) * t78 - pkin(10) * t79 + t57;
t43 = t168 * t73 + t80 * t189 + t95 * t192;
t40 = pkin(10) * t107 + t43;
t17 = t163 * t36 + t167 * t40;
t15 = pkin(11) * t78 + t17;
t39 = -pkin(3) * t107 - t42;
t25 = pkin(4) * t62 - pkin(11) * t63 + t39;
t7 = t166 * t15 + t162 * t25;
t30 = -mrSges(6,1) * t46 + mrSges(6,2) * t47;
t49 = mrSges(5,1) * t78 - mrSges(5,3) * t63;
t201 = t30 - t49;
t198 = Ifges(7,3) * t167;
t126 = -mrSges(6,1) * t166 + mrSges(6,2) * t162;
t197 = -mrSges(5,1) + t126;
t196 = Ifges(7,5) * t106 - Ifges(7,6) * t105;
t195 = t111 * t163;
t113 = t159 * t163 + t167 * t192;
t194 = t113 * t167;
t188 = t162 * t163;
t187 = t163 * t166;
t83 = Ifges(7,5) * t120 + Ifges(7,6) * t119;
t125 = -pkin(4) * t167 - pkin(11) * t163 - pkin(3);
t97 = t162 * t125 + t166 * t203;
t128 = Ifges(6,5) * t162 + Ifges(6,6) * t166;
t186 = Ifges(5,5) * t163 + Ifges(5,6) * t167;
t185 = t162 ^ 2 + t166 ^ 2;
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t62;
t20 = Ifges(6,5) * t47 + Ifges(6,6) * t46 + Ifges(6,3) * t62;
t184 = Ifges(5,5) * t63 - Ifges(5,6) * t62 + Ifges(5,3) * t78;
t183 = Ifges(4,5) * t79 - Ifges(4,6) * t78 + Ifges(4,3) * t107;
t181 = t128 / 0.2e1 + t83 / 0.2e1;
t86 = -t113 * t162 - t166 * t191;
t87 = t113 * t166 - t162 * t191;
t55 = -t161 * t87 + t165 * t86;
t56 = t161 * t86 + t165 * t87;
t180 = t55 * mrSges(7,1) - t56 * mrSges(7,2);
t6 = -t15 * t162 + t166 * t25;
t16 = -t163 * t40 + t167 * t36;
t179 = Ifges(6,5) * t187 - Ifges(6,6) * t188;
t178 = mrSges(6,1) * t162 + mrSges(6,2) * t166;
t176 = -t162 * t86 + t166 * t87;
t118 = t166 * t125;
t77 = -pkin(12) * t187 + t118 + (-pkin(10) * t162 - pkin(5)) * t167;
t88 = -pkin(12) * t188 + t97;
t52 = -t161 * t88 + t165 * t77;
t53 = t161 * t77 + t165 * t88;
t175 = t52 * mrSges(7,1) - t53 * mrSges(7,2) + t196;
t134 = t213 * t162;
t135 = t213 * t166;
t91 = t134 * t165 + t135 * t161;
t92 = t134 * t161 - t135 * t165;
t174 = t91 * mrSges(7,1) - t92 * mrSges(7,2) + t83;
t14 = -pkin(4) * t78 - t16;
t4 = pkin(5) * t62 - pkin(12) * t47 + t6;
t5 = pkin(12) * t46 + t7;
t2 = -t161 * t5 + t165 * t4;
t3 = t161 * t4 + t165 * t5;
t173 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t172 = (mrSges(7,1) * t165 - mrSges(7,2) * t161) * pkin(5);
t170 = pkin(10) ^ 2;
t154 = t167 ^ 2;
t152 = t163 ^ 2;
t150 = t156 ^ 2;
t149 = t152 * t170;
t144 = -pkin(5) * t166 - pkin(4);
t142 = t150 * t168 ^ 2;
t137 = mrSges(3,2) * t193;
t132 = Ifges(5,1) * t163 + Ifges(5,4) * t167;
t130 = Ifges(5,4) * t163 + Ifges(5,2) * t167;
t129 = Ifges(6,2) * t166 + t200;
t127 = -mrSges(5,1) * t167 + mrSges(5,2) * t163;
t124 = (pkin(5) * t162 + pkin(10)) * t163;
t123 = -mrSges(6,1) * t167 - mrSges(6,3) * t187;
t122 = mrSges(6,2) * t167 - mrSges(6,3) * t188;
t114 = t178 * t163;
t109 = -qJ(2) * t193 + t139;
t101 = -Ifges(6,3) * t167 + t179;
t96 = -t162 * t203 + t118;
t94 = -mrSges(7,1) * t167 - mrSges(7,3) * t106;
t93 = mrSges(7,2) * t167 - mrSges(7,3) * t105;
t85 = Ifges(7,1) * t120 + Ifges(7,4) * t119;
t84 = Ifges(7,4) * t120 + Ifges(7,2) * t119;
t82 = -mrSges(7,1) * t119 + mrSges(7,2) * t120;
t69 = mrSges(7,1) * t105 + mrSges(7,2) * t106;
t68 = Ifges(7,1) * t106 - Ifges(7,4) * t105 - Ifges(7,5) * t167;
t67 = Ifges(7,4) * t106 - Ifges(7,2) * t105 - Ifges(7,6) * t167;
t66 = t196 - t198;
t65 = mrSges(4,1) * t107 - mrSges(4,3) * t79;
t64 = -mrSges(4,2) * t107 - mrSges(4,3) * t78;
t51 = mrSges(4,1) * t78 + mrSges(4,2) * t79;
t48 = -mrSges(5,2) * t78 - mrSges(5,3) * t62;
t41 = mrSges(5,1) * t62 + mrSges(5,2) * t63;
t34 = Ifges(5,1) * t63 - Ifges(5,4) * t62 + Ifges(5,5) * t78;
t33 = Ifges(5,4) * t63 - Ifges(5,2) * t62 + Ifges(5,6) * t78;
t32 = mrSges(6,1) * t62 - mrSges(6,3) * t47;
t31 = -mrSges(6,2) * t62 + mrSges(6,3) * t46;
t21 = Ifges(6,4) * t47 + Ifges(6,2) * t46 + Ifges(6,6) * t62;
t19 = mrSges(7,1) * t62 - mrSges(7,3) * t29;
t18 = -mrSges(7,2) * t62 + mrSges(7,3) * t28;
t12 = -pkin(5) * t46 + t14;
t11 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t62;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t62;
t1 = [0.2e1 * t7 * t31 + 0.2e1 * t6 * t32 + t28 * t9 + t29 * t10 + 0.2e1 * t14 * t30 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 + 0.2e1 * t12 * t11 + t63 * t34 + 0.2e1 * t43 * t64 + 0.2e1 * t42 * t65 + 0.2e1 * t57 * t51 + t47 * t22 + 0.2e1 * t17 * t48 + 0.2e1 * t16 * t49 + 0.2e1 * t39 * t41 + t46 * t21 + t79 * (Ifges(4,1) * t79 + Ifges(4,5) * t107) + t107 * t183 + (-0.2e1 * Ifges(4,4) * t79 + Ifges(4,2) * t78 - Ifges(4,6) * t107 + t184) * t78 + m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t14 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t39 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2 + t57 ^ 2) + m(3) * (pkin(1) ^ 2 * t157 ^ 2 + t109 ^ 2 + t110 ^ 2) + (t20 + t8 - t33) * t62 + Ifges(2,3) + (-0.2e1 * mrSges(3,2) * t110 + Ifges(3,3) * t160 + t109 * t218) * t160 + (-0.2e1 * pkin(1) * t137 + (-0.2e1 * mrSges(3,3) * t109 + Ifges(3,1) * t193 + Ifges(3,5) * t217) * t155 + (0.2e1 * t110 * mrSges(3,3) + Ifges(3,6) * t217 + (0.2e1 * Ifges(3,4) * t155 + Ifges(3,2) * t158 + pkin(1) * t218) * t157) * t158) * t157; t113 * t48 + t159 * t51 + t56 * t18 + t55 * t19 + t87 * t31 + t86 * t32 + t137 + (-m(3) * pkin(1) - mrSges(3,1) * t158) * t157 + (t164 * t64 + (-t41 + t65) * t168) * t156 + (t11 + t201) * t111 + m(7) * (t111 * t12 + t2 * t55 + t3 * t56) + m(6) * (t111 * t14 + t6 * t86 + t7 * t87) + m(5) * (-t111 * t16 + t113 * t17 - t191 * t39) + m(4) * (t159 * t57 + (t164 * t43 + t168 * t42) * t156); m(3) + m(7) * (t55 ^ 2 + t56 ^ 2 + t108) + m(6) * (t86 ^ 2 + t87 ^ 2 + t108) + m(5) * (t113 ^ 2 + t108 + t142) + m(4) * (t150 * t164 ^ 2 + t159 ^ 2 + t142); (t33 / 0.2e1 - t20 / 0.2e1 - t8 / 0.2e1 + t17 * mrSges(5,3) + pkin(10) * t48) * t167 + m(5) * (-pkin(3) * t39 + (-t16 * t163 + t167 * t17) * pkin(10)) + t183 + t7 * t122 + t6 * t123 + t124 * t11 + t39 * t127 + t63 * t132 / 0.2e1 + t14 * t114 + t3 * t93 + t2 * t94 + t96 * t32 + t97 * t31 + t12 * t69 + t52 * t19 + t53 * t18 - pkin(3) * t41 + t42 * mrSges(4,1) - t43 * mrSges(4,2) + m(7) * (t12 * t124 + t2 * t52 + t53 * t3) + (-t130 / 0.2e1 + t101 / 0.2e1 + t66 / 0.2e1) * t62 + t78 * t186 / 0.2e1 + m(6) * (t14 * t204 + t6 * t96 + t7 * t97) + (t34 / 0.2e1 - t16 * mrSges(5,3) - t162 * t21 / 0.2e1 + t166 * t216 + t201 * pkin(10)) * t163 + t10 * t209 + t9 * t210 + t47 * t211 + t46 * t212 + t68 * t214 + t67 * t215; mrSges(5,3) * t194 + t87 * t122 + t86 * t123 + t55 * t94 + t56 * t93 + (-t164 * mrSges(4,2) + (mrSges(4,1) - t127) * t168) * t156 + (t163 * mrSges(5,3) + t114 + t69) * t111 + m(7) * (t111 * t124 + t52 * t55 + t53 * t56) + m(6) * (pkin(10) * t195 + t86 * t96 + t87 * t97) + m(5) * (pkin(3) * t191 + (t194 + t195) * pkin(10)); -0.2e1 * pkin(3) * t127 - t105 * t67 + t106 * t68 + 0.2e1 * t97 * t122 + 0.2e1 * t96 * t123 + 0.2e1 * t124 * t69 + 0.2e1 * t52 * t94 + 0.2e1 * t53 * t93 + Ifges(4,3) + (t152 + t154) * mrSges(5,3) * t219 + (-t66 - t101 + t130) * t167 + m(7) * (t124 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(6) * (t96 ^ 2 + t97 ^ 2 + t149) + m(5) * (pkin(3) ^ 2 + t154 * t170 + t149) + (-t102 * t162 + t103 * t166 + t114 * t219 + t132) * t163; -pkin(4) * t30 + t16 * mrSges(5,1) - t17 * mrSges(5,2) + m(6) * (-pkin(4) * t14 + (-t162 * t6 + t166 * t7) * pkin(11)) + t184 + t144 * t11 + t9 * t208 + t10 * t207 + t14 * t126 + t46 * t129 / 0.2e1 + t47 * t206 + t91 * t19 + t92 * t18 + t12 * t82 + t84 * t215 + t85 * t214 + (t21 / 0.2e1 + t7 * mrSges(6,3) + pkin(11) * t31) * t166 + (-t6 * mrSges(6,3) - pkin(11) * t32 + t216) * t162 + t181 * t62 + m(7) * (t12 * t144 + t2 * t91 + t3 * t92) + (t119 * t3 - t120 * t2) * mrSges(7,3); -t113 * mrSges(5,2) + (t119 * t56 - t120 * t55) * mrSges(7,3) + t176 * mrSges(6,3) + (t82 + t197) * t111 + m(7) * (t111 * t144 + t55 * t91 + t56 * t92) + m(6) * (-pkin(4) * t111 + pkin(11) * t176); t144 * t69 + t67 * t208 + t68 * t207 + t124 * t82 - pkin(4) * t114 + t84 * t210 + t85 * t209 + t92 * t93 + t91 * t94 + t197 * t204 + (t119 * t53 - t120 * t52) * mrSges(7,3) + (t97 * mrSges(6,3) + pkin(11) * t122 + t163 * t206 + t212) * t166 + (t211 - t96 * mrSges(6,3) - pkin(11) * t123 - t163 * t129 / 0.2e1) * t162 + m(7) * (t124 * t144 + t52 * t91 + t53 * t92) + m(6) * (-pkin(4) * t204 + (-t162 * t96 + t166 * t97) * pkin(11)) + (-pkin(10) * mrSges(5,2) - t181) * t167 + t186; -0.2e1 * pkin(4) * t126 + t119 * t84 + t120 * t85 + t166 * t129 + t162 * t131 + 0.2e1 * t144 * t82 + Ifges(5,3) + m(7) * (t144 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(6) * (pkin(11) ^ 2 * t185 + pkin(4) ^ 2) + 0.2e1 * (t119 * t92 - t120 * t91) * mrSges(7,3) + 0.2e1 * t185 * pkin(11) * mrSges(6,3); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(7) * (t161 * t3 + t165 * t2) + t165 * t19 + t161 * t18) * pkin(5) + t173 + t20; t86 * mrSges(6,1) - t87 * mrSges(6,2) + m(7) * (t161 * t56 + t165 * t55) * pkin(5) + t180; t96 * mrSges(6,1) - t97 * mrSges(6,2) + t202 * t167 + (m(7) * (t161 * t53 + t165 * t52) + t165 * t94 + t161 * t93) * pkin(5) + t175 + t179; -t178 * pkin(11) + (m(7) * (t161 * t92 + t165 * t91) + (t119 * t161 - t120 * t165) * mrSges(7,3)) * pkin(5) + t174 + t128; m(7) * (t161 ^ 2 + t165 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t172 - t202; t173; t180; t175 - t198; t174; Ifges(7,3) + t172; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
