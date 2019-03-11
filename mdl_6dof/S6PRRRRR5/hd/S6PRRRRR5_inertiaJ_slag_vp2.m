% Calculate joint inertia matrix for
% S6PRRRRR5
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR5_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR5_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR5_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:04:01
% EndTime: 2019-03-09 01:04:05
% DurationCPUTime: 1.63s
% Computational Cost: add. (2918->411), mult. (7056->600), div. (0->0), fcn. (7877->14), ass. (0->157)
t202 = 2 * pkin(10);
t149 = sin(qJ(4));
t154 = cos(qJ(4));
t144 = sin(pkin(6));
t146 = cos(pkin(6));
t150 = sin(qJ(3));
t151 = sin(qJ(2));
t155 = cos(qJ(3));
t145 = cos(pkin(7));
t156 = cos(qJ(2));
t174 = t145 * t156;
t143 = sin(pkin(7));
t176 = t143 * t150;
t55 = t146 * t176 + (t150 * t174 + t151 * t155) * t144;
t92 = -t143 * t144 * t156 + t145 * t146;
t35 = t149 * t55 - t154 * t92;
t34 = t35 ^ 2;
t175 = t143 * t155;
t53 = -t146 * t175 + (t150 * t151 - t155 * t174) * t144;
t201 = t53 ^ 2;
t148 = sin(qJ(5));
t153 = cos(qJ(5));
t94 = t145 * t149 + t154 * t176;
t62 = -t148 * t94 - t153 * t175;
t63 = -t148 * t175 + t153 * t94;
t93 = -t145 * t154 + t149 * t176;
t22 = Ifges(6,1) * t63 + Ifges(6,4) * t62 + Ifges(6,5) * t93;
t200 = t22 / 0.2e1;
t147 = sin(qJ(6));
t152 = cos(qJ(6));
t28 = -t147 * t63 + t152 * t62;
t199 = t28 / 0.2e1;
t29 = t147 * t62 + t152 * t63;
t198 = t29 / 0.2e1;
t181 = Ifges(6,4) * t153;
t84 = -Ifges(6,6) * t154 + (-Ifges(6,2) * t148 + t181) * t149;
t197 = t84 / 0.2e1;
t182 = Ifges(6,4) * t148;
t85 = -Ifges(6,5) * t154 + (Ifges(6,1) * t153 - t182) * t149;
t196 = t85 / 0.2e1;
t106 = t147 * t153 + t148 * t152;
t86 = t106 * t149;
t195 = -t86 / 0.2e1;
t105 = -t147 * t148 + t152 * t153;
t87 = t105 * t149;
t194 = t87 / 0.2e1;
t193 = -pkin(12) - pkin(11);
t192 = t105 / 0.2e1;
t191 = t106 / 0.2e1;
t117 = Ifges(6,1) * t148 + t181;
t190 = t117 / 0.2e1;
t189 = pkin(2) * t155;
t188 = pkin(10) * t149;
t187 = pkin(10) * t154;
t186 = -Ifges(7,3) - Ifges(6,3);
t124 = pkin(9) * t176;
t78 = t124 + (-pkin(3) - t189) * t145;
t38 = pkin(4) * t93 - pkin(11) * t94 + t78;
t97 = pkin(2) * t145 * t150 + pkin(9) * t175;
t79 = pkin(10) * t145 + t97;
t80 = (-pkin(3) * t155 - pkin(10) * t150 - pkin(2)) * t143;
t44 = t149 * t80 + t154 * t79;
t40 = -pkin(11) * t175 + t44;
t14 = t148 * t38 + t153 * t40;
t30 = -mrSges(6,1) * t62 + mrSges(6,2) * t63;
t67 = -mrSges(5,1) * t175 - mrSges(5,3) * t94;
t185 = -t67 + t30;
t184 = Ifges(7,5) * t87 - Ifges(7,6) * t86;
t183 = -Ifges(5,5) * t94 + Ifges(5,6) * t93;
t180 = Ifges(7,3) * t154;
t179 = t35 * t149;
t37 = t149 * t92 + t154 * t55;
t178 = t37 * t154;
t112 = -mrSges(6,1) * t153 + mrSges(6,2) * t148;
t177 = -mrSges(5,1) + t112;
t59 = Ifges(7,5) * t106 + Ifges(7,6) * t105;
t173 = t148 * t149;
t172 = t149 * t153;
t111 = -pkin(4) * t154 - pkin(11) * t149 - pkin(3);
t76 = t111 * t148 + t153 * t187;
t114 = Ifges(6,5) * t148 + Ifges(6,6) * t153;
t171 = Ifges(5,5) * t149 + Ifges(5,6) * t154;
t170 = t148 ^ 2 + t153 ^ 2;
t8 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t93;
t20 = Ifges(6,5) * t63 + Ifges(6,6) * t62 + Ifges(6,3) * t93;
t15 = -t148 * t37 + t153 * t53;
t16 = t148 * t53 + t153 * t37;
t5 = -t147 * t16 + t15 * t152;
t6 = t147 * t15 + t152 * t16;
t169 = mrSges(7,1) * t5 - t6 * mrSges(7,2);
t168 = t114 / 0.2e1 + t59 / 0.2e1;
t167 = Ifges(4,5) * t176 + Ifges(4,6) * t175 + Ifges(4,3) * t145;
t13 = -t148 * t40 + t153 * t38;
t43 = -t149 * t79 + t154 * t80;
t39 = pkin(4) * t175 - t43;
t166 = Ifges(6,5) * t172 - Ifges(6,6) * t173;
t165 = mrSges(6,1) * t148 + mrSges(6,2) * t153;
t164 = -t148 * t15 + t153 * t16;
t102 = t153 * t111;
t52 = -pkin(12) * t172 + t102 + (-pkin(10) * t148 - pkin(5)) * t154;
t64 = -pkin(12) * t173 + t76;
t24 = -t147 * t64 + t152 * t52;
t25 = t147 * t52 + t152 * t64;
t163 = mrSges(7,1) * t24 - t25 * mrSges(7,2) + t184;
t120 = t193 * t148;
t121 = t193 * t153;
t68 = t120 * t152 + t121 * t147;
t69 = t120 * t147 - t121 * t152;
t162 = mrSges(7,1) * t68 - t69 * mrSges(7,2) + t59;
t11 = pkin(12) * t62 + t14;
t7 = pkin(5) * t93 - pkin(12) * t63 + t13;
t2 = -t11 * t147 + t152 * t7;
t3 = t11 * t152 + t147 * t7;
t161 = mrSges(7,1) * t2 - t3 * mrSges(7,2) + t8;
t160 = (mrSges(7,1) * t152 - mrSges(7,2) * t147) * pkin(5);
t158 = pkin(10) ^ 2;
t142 = t154 ^ 2;
t140 = t149 ^ 2;
t137 = t140 * t158;
t131 = -pkin(5) * t153 - pkin(4);
t118 = Ifges(5,1) * t149 + Ifges(5,4) * t154;
t116 = Ifges(5,4) * t149 + Ifges(5,2) * t154;
t115 = Ifges(6,2) * t153 + t182;
t113 = -mrSges(5,1) * t154 + mrSges(5,2) * t149;
t110 = (pkin(5) * t148 + pkin(10)) * t149;
t108 = -mrSges(6,1) * t154 - mrSges(6,3) * t172;
t107 = mrSges(6,2) * t154 - mrSges(6,3) * t173;
t104 = -mrSges(4,2) * t145 + mrSges(4,3) * t175;
t103 = mrSges(4,1) * t145 - mrSges(4,3) * t176;
t98 = t165 * t149;
t96 = t145 * t189 - t124;
t95 = (-mrSges(4,1) * t155 + mrSges(4,2) * t150) * t143;
t83 = -Ifges(6,3) * t154 + t166;
t75 = -t148 * t187 + t102;
t71 = -mrSges(7,1) * t154 - mrSges(7,3) * t87;
t70 = mrSges(7,2) * t154 - mrSges(7,3) * t86;
t66 = mrSges(5,2) * t175 - mrSges(5,3) * t93;
t61 = Ifges(7,1) * t106 + Ifges(7,4) * t105;
t60 = Ifges(7,4) * t106 + Ifges(7,2) * t105;
t58 = -mrSges(7,1) * t105 + mrSges(7,2) * t106;
t51 = mrSges(5,1) * t93 + mrSges(5,2) * t94;
t50 = mrSges(7,1) * t86 + mrSges(7,2) * t87;
t49 = Ifges(5,1) * t94 - Ifges(5,4) * t93 - Ifges(5,5) * t175;
t48 = Ifges(5,4) * t94 - Ifges(5,2) * t93 - Ifges(5,6) * t175;
t47 = Ifges(7,1) * t87 - Ifges(7,4) * t86 - Ifges(7,5) * t154;
t46 = Ifges(7,4) * t87 - Ifges(7,2) * t86 - Ifges(7,6) * t154;
t45 = -t180 + t184;
t42 = mrSges(6,1) * t93 - mrSges(6,3) * t63;
t41 = -mrSges(6,2) * t93 + mrSges(6,3) * t62;
t21 = Ifges(6,4) * t63 + Ifges(6,2) * t62 + Ifges(6,6) * t93;
t19 = -pkin(5) * t62 + t39;
t18 = mrSges(7,1) * t93 - mrSges(7,3) * t29;
t17 = -mrSges(7,2) * t93 + mrSges(7,3) * t28;
t12 = -mrSges(7,1) * t28 + mrSges(7,2) * t29;
t10 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t93;
t9 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t93;
t1 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t34) + m(6) * (t15 ^ 2 + t16 ^ 2 + t34) + m(5) * (t37 ^ 2 + t201 + t34) + m(4) * (t55 ^ 2 + t92 ^ 2 + t201) + m(3) * (t146 ^ 2 + (t151 ^ 2 + t156 ^ 2) * t144 ^ 2); t55 * t104 + t15 * t42 + t16 * t41 + t6 * t17 + t5 * t18 + t37 * t66 + t92 * t95 + (-t103 + t51) * t53 + (mrSges(3,1) * t156 - mrSges(3,2) * t151) * t144 + (t12 + t185) * t35 + m(7) * (t19 * t35 + t2 * t5 + t3 * t6) + m(6) * (t13 * t15 + t14 * t16 + t35 * t39) + m(5) * (-t35 * t43 + t37 * t44 + t53 * t78) + m(4) * (-pkin(2) * t143 * t92 - t53 * t96 + t55 * t97); t145 * t167 + 0.2e1 * t96 * t103 + 0.2e1 * t97 * t104 + t94 * t49 + 0.2e1 * t44 * t66 + 0.2e1 * t43 * t67 + 0.2e1 * t78 * t51 + t62 * t21 + t63 * t22 + 0.2e1 * t39 * t30 + 0.2e1 * t14 * t41 + 0.2e1 * t13 * t42 + t28 * t9 + t29 * t10 + 0.2e1 * t3 * t17 + 0.2e1 * t2 * t18 + 0.2e1 * t19 * t12 + Ifges(3,3) + (t20 + t8 - t48) * t93 + (-0.2e1 * pkin(2) * t95 + (Ifges(4,1) * t176 + Ifges(4,5) * t145) * t150 + (Ifges(4,6) * t145 + (0.2e1 * Ifges(4,4) * t150 + (Ifges(4,2) + Ifges(5,3)) * t155) * t143 + t183) * t155) * t143 + m(4) * (pkin(2) ^ 2 * t143 ^ 2 + t96 ^ 2 + t97 ^ 2) + m(5) * (t43 ^ 2 + t44 ^ 2 + t78 ^ 2) + m(6) * (t13 ^ 2 + t14 ^ 2 + t39 ^ 2) + m(7) * (t19 ^ 2 + t2 ^ 2 + t3 ^ 2); mrSges(5,3) * t178 - t55 * mrSges(4,2) + t16 * t107 + t15 * t108 + t5 * t71 + t6 * t70 + (-mrSges(4,1) + t113) * t53 + (t149 * mrSges(5,3) + t50 + t98) * t35 + m(7) * (t110 * t35 + t24 * t5 + t25 * t6) + m(6) * (pkin(10) * t179 + t15 * t75 + t16 * t76) + m(5) * (-pkin(3) * t53 + (t178 + t179) * pkin(10)); t94 * t118 / 0.2e1 + t14 * t107 + t13 * t108 + t110 * t12 + t78 * t113 + t96 * mrSges(4,1) - t97 * mrSges(4,2) + t39 * t98 + t3 * t70 + t2 * t71 + t75 * t42 + t76 * t41 + t19 * t50 - pkin(3) * t51 + t24 * t18 + t25 * t17 + m(5) * (-pkin(3) * t78 + (-t149 * t43 + t154 * t44) * pkin(10)) + (t153 * t200 - t148 * t21 / 0.2e1 - t43 * mrSges(5,3) + t49 / 0.2e1 + t185 * pkin(10)) * t149 + m(6) * (t13 * t75 + t14 * t76 + t188 * t39) - t171 * t175 / 0.2e1 + t167 + (pkin(10) * t66 + t44 * mrSges(5,3) + t48 / 0.2e1 - t20 / 0.2e1 - t8 / 0.2e1) * t154 + m(7) * (t110 * t19 + t2 * t24 + t25 * t3) + (-t116 / 0.2e1 + t83 / 0.2e1 + t45 / 0.2e1) * t93 + t10 * t194 + t9 * t195 + t63 * t196 + t62 * t197 + t47 * t198 + t46 * t199; -0.2e1 * pkin(3) * t113 + 0.2e1 * t76 * t107 + 0.2e1 * t75 * t108 + 0.2e1 * t110 * t50 + 0.2e1 * t24 * t71 + 0.2e1 * t25 * t70 - t86 * t46 + t87 * t47 + Ifges(4,3) + (t140 + t142) * mrSges(5,3) * t202 + (-t45 - t83 + t116) * t154 + m(7) * (t110 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t75 ^ 2 + t76 ^ 2 + t137) + m(5) * (pkin(3) ^ 2 + t142 * t158 + t137) + (-t148 * t84 + t153 * t85 + t202 * t98 + t118) * t149; -t37 * mrSges(5,2) + (t105 * t6 - t106 * t5) * mrSges(7,3) + t164 * mrSges(6,3) + (t58 + t177) * t35 + m(7) * (t131 * t35 + t5 * t68 + t6 * t69) + m(6) * (-pkin(4) * t35 + pkin(11) * t164); -t183 + t63 * t190 + t131 * t12 + t9 * t192 + t10 * t191 + t39 * t112 + t62 * t115 / 0.2e1 + t68 * t18 + t69 * t17 + t19 * t58 + t60 * t199 + t61 * t198 + t43 * mrSges(5,1) - t44 * mrSges(5,2) - pkin(4) * t30 + m(6) * (-pkin(4) * t39 + (-t13 * t148 + t14 * t153) * pkin(11)) - Ifges(5,3) * t175 + t168 * t93 + (pkin(11) * t41 + t14 * mrSges(6,3) + t21 / 0.2e1) * t153 + (-t13 * mrSges(6,3) - pkin(11) * t42 + t200) * t148 + m(7) * (t131 * t19 + t2 * t68 + t3 * t69) + (t105 * t3 - t106 * t2) * mrSges(7,3); t131 * t50 + t46 * t192 + t47 * t191 + t110 * t58 - pkin(4) * t98 + t60 * t195 + t61 * t194 + t69 * t70 + t68 * t71 + t177 * t188 + (t105 * t25 - t106 * t24) * mrSges(7,3) + (t76 * mrSges(6,3) + pkin(11) * t107 + t149 * t190 + t197) * t153 + (-t149 * t115 / 0.2e1 - pkin(11) * t108 - t75 * mrSges(6,3) + t196) * t148 + m(6) * (-pkin(4) * t188 + (-t148 * t75 + t153 * t76) * pkin(11)) + m(7) * (t110 * t131 + t24 * t68 + t25 * t69) + (-pkin(10) * mrSges(5,2) - t168) * t154 + t171; -0.2e1 * pkin(4) * t112 + t105 * t60 + t106 * t61 + t153 * t115 + t148 * t117 + 0.2e1 * t131 * t58 + Ifges(5,3) + m(7) * (t131 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (pkin(11) ^ 2 * t170 + pkin(4) ^ 2) + 0.2e1 * (t105 * t69 - t106 * t68) * mrSges(7,3) + 0.2e1 * t170 * pkin(11) * mrSges(6,3); t15 * mrSges(6,1) - t16 * mrSges(6,2) + m(7) * (t147 * t6 + t152 * t5) * pkin(5) + t169; t13 * mrSges(6,1) - t14 * mrSges(6,2) + (m(7) * (t147 * t3 + t152 * t2) + t147 * t17 + t152 * t18) * pkin(5) + t161 + t20; t75 * mrSges(6,1) - t76 * mrSges(6,2) + t186 * t154 + (m(7) * (t147 * t25 + t152 * t24) + t147 * t70 + t152 * t71) * pkin(5) + t163 + t166; -t165 * pkin(11) + (m(7) * (t147 * t69 + t152 * t68) + (t105 * t147 - t106 * t152) * mrSges(7,3)) * pkin(5) + t162 + t114; m(7) * (t147 ^ 2 + t152 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t160 - t186; t169; t161; t163 - t180; t162; Ifges(7,3) + t160; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
