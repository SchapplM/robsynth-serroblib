% Calculate joint inertia matrix for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRPR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:14
% EndTime: 2019-03-09 05:39:18
% DurationCPUTime: 1.92s
% Computational Cost: add. (5890->412), mult. (15573->612), div. (0->0), fcn. (17908->14), ass. (0->160)
t206 = 2 * pkin(10);
t156 = cos(pkin(13));
t161 = sin(qJ(4));
t177 = t156 * t161;
t152 = sin(pkin(13));
t182 = t152 * t161;
t110 = mrSges(6,1) * t182 + mrSges(6,2) * t177;
t160 = sin(qJ(6));
t163 = cos(qJ(6));
t117 = t152 * t163 + t156 * t160;
t101 = t117 * t161;
t116 = -t152 * t160 + t156 * t163;
t102 = t116 * t161;
t65 = mrSges(7,1) * t101 + mrSges(7,2) * t102;
t205 = t110 + t65;
t154 = sin(pkin(7));
t158 = cos(pkin(7));
t162 = sin(qJ(3));
t165 = cos(qJ(3));
t153 = sin(pkin(12));
t155 = sin(pkin(6));
t157 = cos(pkin(12));
t178 = t155 * t157;
t159 = cos(pkin(6));
t192 = pkin(1) * t159;
t106 = qJ(2) * t178 + t153 * t192;
t171 = t158 * t178;
t69 = (t154 * t159 + t171) * pkin(9) + t106;
t136 = t157 * t192;
t181 = t153 * t155;
t76 = pkin(2) * t159 + t136 + (-pkin(9) * t158 - qJ(2)) * t181;
t90 = (-pkin(9) * t153 * t154 - pkin(2) * t157 - pkin(1)) * t155;
t42 = -t162 * t69 + (t154 * t90 + t158 * t76) * t165;
t204 = 2 * mrSges(3,1);
t164 = cos(qJ(4));
t180 = t154 * t162;
t107 = -t158 * t164 + t161 * t180;
t104 = t107 ^ 2;
t203 = 0.2e1 * t159;
t103 = -t154 * t178 + t158 * t159;
t176 = t158 * t162;
t75 = t159 * t180 + (t153 * t165 + t157 * t176) * t155;
t59 = t103 * t161 + t164 * t75;
t179 = t154 * t165;
t74 = -t159 * t179 + t162 * t181 - t165 * t171;
t45 = -t152 * t59 + t156 * t74;
t46 = t152 * t74 + t156 * t59;
t58 = -t103 * t164 + t161 * t75;
t21 = Ifges(6,1) * t46 + Ifges(6,4) * t45 + Ifges(6,5) * t58;
t202 = t21 / 0.2e1;
t28 = -t160 * t46 + t163 * t45;
t201 = t28 / 0.2e1;
t29 = t160 * t45 + t163 * t46;
t200 = t29 / 0.2e1;
t186 = Ifges(6,4) * t156;
t98 = -Ifges(6,6) * t164 + (-Ifges(6,2) * t152 + t186) * t161;
t199 = t98 / 0.2e1;
t187 = Ifges(6,4) * t152;
t99 = -Ifges(6,5) * t164 + (Ifges(6,1) * t156 - t187) * t161;
t198 = t99 / 0.2e1;
t197 = -t101 / 0.2e1;
t196 = t102 / 0.2e1;
t195 = t116 / 0.2e1;
t194 = t117 / 0.2e1;
t128 = Ifges(6,1) * t152 + t186;
t193 = t128 / 0.2e1;
t191 = pkin(10) * t161;
t190 = pkin(10) * t164;
t189 = pkin(11) + qJ(5);
t54 = -t154 * t76 + t158 * t90;
t36 = pkin(3) * t74 - pkin(10) * t75 + t54;
t43 = t165 * t69 + t176 * t76 + t180 * t90;
t40 = pkin(10) * t103 + t43;
t16 = t161 * t36 + t164 * t40;
t13 = qJ(5) * t74 + t16;
t39 = -pkin(3) * t103 - t42;
t24 = pkin(4) * t58 - qJ(5) * t59 + t39;
t6 = t13 * t156 + t152 * t24;
t30 = -t45 * mrSges(6,1) + mrSges(6,2) * t46;
t48 = mrSges(5,1) * t74 - mrSges(5,3) * t59;
t188 = t30 - t48;
t124 = -mrSges(6,1) * t156 + mrSges(6,2) * t152;
t185 = -mrSges(5,1) + t124;
t184 = t107 * t161;
t109 = t158 * t161 + t164 * t180;
t183 = t109 * t164;
t78 = Ifges(7,5) * t117 + Ifges(7,6) * t116;
t122 = -pkin(4) * t164 - qJ(5) * t161 - pkin(3);
t92 = t122 * t152 + t156 * t190;
t175 = Ifges(5,5) * t161 + Ifges(5,6) * t164;
t174 = t152 ^ 2 + t156 ^ 2;
t7 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t58;
t173 = Ifges(5,5) * t59 - Ifges(5,6) * t58 + Ifges(5,3) * t74;
t172 = Ifges(4,5) * t75 - Ifges(4,6) * t74 + Ifges(4,3) * t103;
t170 = t78 / 0.2e1 + Ifges(6,5) * t152 / 0.2e1 + Ifges(6,6) * t156 / 0.2e1;
t10 = -t28 * mrSges(7,1) + mrSges(7,2) * t29;
t5 = -t13 * t152 + t156 * t24;
t15 = -t161 * t40 + t164 * t36;
t77 = -mrSges(7,1) * t116 + mrSges(7,2) * t117;
t62 = Ifges(7,5) * t102 - Ifges(7,6) * t101 - Ifges(7,3) * t164;
t81 = -t109 * t152 - t156 * t179;
t82 = t109 * t156 - t152 * t179;
t169 = -t152 * t81 + t156 * t82;
t14 = -pkin(4) * t74 - t15;
t167 = pkin(10) ^ 2;
t151 = t164 ^ 2;
t150 = t161 ^ 2;
t148 = t154 ^ 2;
t146 = t150 * t167;
t142 = -pkin(5) * t156 - pkin(4);
t140 = t148 * t165 ^ 2;
t134 = mrSges(3,2) * t181;
t131 = Ifges(5,1) * t161 + Ifges(5,4) * t164;
t130 = Ifges(5,4) * t161 + Ifges(5,2) * t164;
t129 = -mrSges(5,1) * t164 + mrSges(5,2) * t161;
t127 = Ifges(6,2) * t156 + t187;
t125 = t189 * t156;
t123 = t189 * t152;
t121 = (pkin(5) * t152 + pkin(10)) * t161;
t119 = -mrSges(6,1) * t164 - mrSges(6,3) * t177;
t118 = mrSges(6,2) * t164 - mrSges(6,3) * t182;
t115 = t156 * t122;
t105 = -qJ(2) * t181 + t136;
t97 = -Ifges(6,3) * t164 + (Ifges(6,5) * t156 - Ifges(6,6) * t152) * t161;
t91 = -t152 * t190 + t115;
t89 = -mrSges(7,1) * t164 - mrSges(7,3) * t102;
t88 = mrSges(7,2) * t164 - mrSges(7,3) * t101;
t86 = -t123 * t160 + t125 * t163;
t85 = -t123 * t163 - t125 * t160;
t84 = -pkin(11) * t182 + t92;
t80 = Ifges(7,1) * t117 + Ifges(7,4) * t116;
t79 = Ifges(7,4) * t117 + Ifges(7,2) * t116;
t70 = -pkin(11) * t177 + t115 + (-pkin(10) * t152 - pkin(5)) * t164;
t64 = Ifges(7,1) * t102 - Ifges(7,4) * t101 - Ifges(7,5) * t164;
t63 = Ifges(7,4) * t102 - Ifges(7,2) * t101 - Ifges(7,6) * t164;
t61 = mrSges(4,1) * t103 - mrSges(4,3) * t75;
t60 = -mrSges(4,2) * t103 - mrSges(4,3) * t74;
t53 = t160 * t81 + t163 * t82;
t52 = -t160 * t82 + t163 * t81;
t51 = mrSges(4,1) * t74 + mrSges(4,2) * t75;
t50 = t160 * t70 + t163 * t84;
t49 = -t160 * t84 + t163 * t70;
t47 = -mrSges(5,2) * t74 - mrSges(5,3) * t58;
t41 = mrSges(5,1) * t58 + mrSges(5,2) * t59;
t34 = Ifges(5,1) * t59 - Ifges(5,4) * t58 + Ifges(5,5) * t74;
t33 = Ifges(5,4) * t59 - Ifges(5,2) * t58 + Ifges(5,6) * t74;
t32 = mrSges(6,1) * t58 - mrSges(6,3) * t46;
t31 = -mrSges(6,2) * t58 + mrSges(6,3) * t45;
t20 = Ifges(6,4) * t46 + Ifges(6,2) * t45 + Ifges(6,6) * t58;
t19 = Ifges(6,5) * t46 + Ifges(6,6) * t45 + Ifges(6,3) * t58;
t18 = mrSges(7,1) * t58 - mrSges(7,3) * t29;
t17 = -mrSges(7,2) * t58 + mrSges(7,3) * t28;
t11 = -pkin(5) * t45 + t14;
t9 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t58;
t8 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t58;
t4 = pkin(11) * t45 + t6;
t3 = pkin(5) * t58 - pkin(11) * t46 + t5;
t2 = t160 * t3 + t163 * t4;
t1 = -t160 * t4 + t163 * t3;
t12 = [t59 * t34 + 0.2e1 * t43 * t60 + 0.2e1 * t42 * t61 + 0.2e1 * t54 * t51 + t46 * t21 + 0.2e1 * t16 * t47 + 0.2e1 * t15 * t48 + (-0.2e1 * mrSges(3,2) * t106 + Ifges(3,3) * t159 + t105 * t204) * t159 + (-0.2e1 * pkin(1) * t134 + (-0.2e1 * mrSges(3,3) * t105 + Ifges(3,1) * t181 + Ifges(3,5) * t203) * t153 + (0.2e1 * t106 * mrSges(3,3) + Ifges(3,6) * t203 + (0.2e1 * Ifges(3,4) * t153 + Ifges(3,2) * t157 + pkin(1) * t204) * t155) * t157) * t155 + m(7) * (t1 ^ 2 + t11 ^ 2 + t2 ^ 2) + m(6) * (t14 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t15 ^ 2 + t16 ^ 2 + t39 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2 + t54 ^ 2) + m(3) * (pkin(1) ^ 2 * t155 ^ 2 + t105 ^ 2 + t106 ^ 2) + (t7 + t19 - t33) * t58 + t75 * (Ifges(4,1) * t75 + Ifges(4,5) * t103) + t103 * t172 + (-0.2e1 * Ifges(4,4) * t75 + Ifges(4,2) * t74 - Ifges(4,6) * t103 + t173) * t74 + Ifges(2,3) + 0.2e1 * t11 * t10 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + t28 * t8 + t29 * t9 + 0.2e1 * t14 * t30 + 0.2e1 * t6 * t31 + 0.2e1 * t5 * t32 + 0.2e1 * t39 * t41 + t45 * t20; t109 * t47 + t158 * t51 + t53 * t17 + t52 * t18 + t82 * t31 + t81 * t32 + t134 + (-m(3) * pkin(1) - mrSges(3,1) * t157) * t155 + (t162 * t60 + (-t41 + t61) * t165) * t154 + (t10 + t188) * t107 + m(7) * (t1 * t52 + t107 * t11 + t2 * t53) + m(6) * (t107 * t14 + t5 * t81 + t6 * t82) + m(5) * (-t107 * t15 + t109 * t16 - t179 * t39) + m(4) * (t158 * t54 + (t162 * t43 + t165 * t42) * t154); m(3) + m(7) * (t52 ^ 2 + t53 ^ 2 + t104) + m(6) * (t81 ^ 2 + t82 ^ 2 + t104) + m(5) * (t109 ^ 2 + t104 + t140) + m(4) * (t148 * t162 ^ 2 + t158 ^ 2 + t140); t14 * t110 + t91 * t32 + t92 * t31 + t2 * t88 + t1 * t89 + t11 * t65 + t49 * t18 + t50 * t17 + (-t7 / 0.2e1 - t19 / 0.2e1 + t33 / 0.2e1 + t16 * mrSges(5,3) + pkin(10) * t47) * t164 + m(7) * (t1 * t49 + t11 * t121 + t2 * t50) + (t97 / 0.2e1 + t62 / 0.2e1 - t130 / 0.2e1) * t58 + t172 + m(6) * (t14 * t191 + t5 * t91 + t6 * t92) + t74 * t175 / 0.2e1 + (t34 / 0.2e1 - t15 * mrSges(5,3) - t152 * t20 / 0.2e1 + t156 * t202 + t188 * pkin(10)) * t161 + t39 * t129 + t59 * t131 / 0.2e1 + t6 * t118 + t5 * t119 + t121 * t10 + m(5) * (-pkin(3) * t39 + (-t15 * t161 + t16 * t164) * pkin(10)) + t9 * t196 + t8 * t197 + t46 * t198 + t45 * t199 + t64 * t200 + t63 * t201 - pkin(3) * t41 + t42 * mrSges(4,1) - t43 * mrSges(4,2); mrSges(5,3) * t183 + t82 * t118 + t81 * t119 + t52 * t89 + t53 * t88 + (-t162 * mrSges(4,2) + (mrSges(4,1) - t129) * t165) * t154 + (t161 * mrSges(5,3) + t205) * t107 + m(7) * (t107 * t121 + t49 * t52 + t50 * t53) + m(6) * (pkin(10) * t184 + t81 * t91 + t82 * t92) + m(5) * (pkin(3) * t179 + (t183 + t184) * pkin(10)); -0.2e1 * pkin(3) * t129 - t101 * t63 + t102 * t64 + 0.2e1 * t92 * t118 + 0.2e1 * t91 * t119 + 0.2e1 * t121 * t65 + 0.2e1 * t49 * t89 + 0.2e1 * t50 * t88 + Ifges(4,3) + (t150 + t151) * mrSges(5,3) * t206 + (-t62 + t130 - t97) * t164 + m(7) * (t121 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t91 ^ 2 + t92 ^ 2 + t146) + m(5) * (pkin(3) ^ 2 + t151 * t167 + t146) + (t110 * t206 - t152 * t98 + t156 * t99 + t131) * t161; t11 * t77 + t79 * t201 + t80 * t200 + t85 * t18 + t86 * t17 + (-t1 * t117 + t116 * t2) * mrSges(7,3) + t173 + (t20 / 0.2e1 + t6 * mrSges(6,3) + qJ(5) * t31) * t156 + (-t5 * mrSges(6,3) - qJ(5) * t32 + t202) * t152 + m(7) * (t1 * t85 + t11 * t142 + t2 * t86) + t170 * t58 + t45 * t127 / 0.2e1 + t46 * t193 + t142 * t10 + t8 * t195 + t9 * t194 + t14 * t124 + m(6) * (-pkin(4) * t14 + (-t152 * t5 + t156 * t6) * qJ(5)) + t15 * mrSges(5,1) - t16 * mrSges(5,2) - pkin(4) * t30; -t109 * mrSges(5,2) + (t116 * t53 - t117 * t52) * mrSges(7,3) + t169 * mrSges(6,3) + (t77 + t185) * t107 + m(7) * (t107 * t142 + t52 * t85 + t53 * t86) + m(6) * (-pkin(4) * t107 + qJ(5) * t169); t79 * t197 + t80 * t196 - pkin(4) * t110 + t86 * t88 + t85 * t89 + t142 * t65 + t63 * t195 + t64 * t194 + t121 * t77 + t185 * t191 + (t116 * t50 - t117 * t49) * mrSges(7,3) + (t92 * mrSges(6,3) + qJ(5) * t118 + t161 * t193 + t199) * t156 + (t198 - t91 * mrSges(6,3) - qJ(5) * t119 - t161 * t127 / 0.2e1) * t152 + m(7) * (t121 * t142 + t49 * t85 + t50 * t86) + m(6) * (-pkin(4) * t191 + (-t152 * t91 + t156 * t92) * qJ(5)) + (-mrSges(5,2) * pkin(10) - t170) * t164 + t175; -0.2e1 * pkin(4) * t124 + t116 * t79 + t117 * t80 + t156 * t127 + t152 * t128 + 0.2e1 * t142 * t77 + Ifges(5,3) + m(7) * (t142 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(6) * (qJ(5) ^ 2 * t174 + pkin(4) ^ 2) + 0.2e1 * (t116 * t86 - t117 * t85) * mrSges(7,3) + 0.2e1 * t174 * qJ(5) * mrSges(6,3); m(6) * t14 + m(7) * t11 + t10 + t30; 0.2e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * t107; m(6) * t191 + m(7) * t121 + t205; -m(6) * pkin(4) + m(7) * t142 + t124 + t77; m(6) + m(7); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t52 - mrSges(7,2) * t53; mrSges(7,1) * t49 - mrSges(7,2) * t50 + t62; mrSges(7,1) * t85 - mrSges(7,2) * t86 + t78; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t12(1) t12(2) t12(4) t12(7) t12(11) t12(16); t12(2) t12(3) t12(5) t12(8) t12(12) t12(17); t12(4) t12(5) t12(6) t12(9) t12(13) t12(18); t12(7) t12(8) t12(9) t12(10) t12(14) t12(19); t12(11) t12(12) t12(13) t12(14) t12(15) t12(20); t12(16) t12(17) t12(18) t12(19) t12(20) t12(21);];
Mq  = res;
