% Calculate joint inertia matrix for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
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
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:11:42
% EndTime: 2019-03-09 16:11:47
% DurationCPUTime: 2.20s
% Computational Cost: add. (2707->430), mult. (6129->593), div. (0->0), fcn. (6485->10), ass. (0->156)
t154 = sin(pkin(11));
t156 = cos(pkin(11));
t204 = t154 ^ 2 + t156 ^ 2;
t203 = 2 * pkin(9);
t155 = sin(pkin(6));
t201 = t155 ^ 2;
t157 = sin(qJ(6));
t160 = cos(qJ(6));
t162 = cos(qJ(2));
t182 = t155 * t162;
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t159 = sin(qJ(2));
t183 = t155 * t159;
t185 = cos(pkin(6));
t97 = t185 * t158 + t161 * t183;
t62 = t154 * t97 + t156 * t182;
t63 = -t154 * t182 + t156 * t97;
t28 = -t157 * t63 + t160 * t62;
t200 = t28 / 0.2e1;
t29 = t157 * t62 + t160 * t63;
t199 = t29 / 0.2e1;
t107 = t154 * t160 - t156 * t157;
t90 = t107 * t158;
t198 = t90 / 0.2e1;
t167 = t154 * t157 + t156 * t160;
t91 = t167 * t158;
t197 = t91 / 0.2e1;
t196 = pkin(4) + pkin(5);
t195 = -t167 / 0.2e1;
t194 = t107 / 0.2e1;
t193 = pkin(9) * t158;
t192 = pkin(9) * t161;
t191 = -pkin(10) + qJ(4);
t173 = pkin(1) * t185;
t101 = -pkin(8) * t183 + t162 * t173;
t87 = -t185 * pkin(2) - t101;
t96 = t158 * t183 - t185 * t161;
t35 = t96 * pkin(3) - t97 * qJ(4) + t87;
t165 = (-pkin(2) * t162 - pkin(9) * t159 - pkin(1)) * t155;
t102 = pkin(8) * t182 + t159 * t173;
t88 = t185 * pkin(9) + t102;
t43 = t158 * t165 + t161 * t88;
t36 = -qJ(4) * t182 + t43;
t14 = t154 * t35 + t156 * t36;
t190 = Ifges(5,4) * t154;
t189 = Ifges(5,4) * t156;
t188 = Ifges(6,5) * t154;
t187 = Ifges(6,5) * t156;
t186 = t154 * mrSges(6,3);
t184 = t154 * t158;
t181 = t156 * t158;
t59 = Ifges(7,5) * t107 - Ifges(7,6) * t167;
t110 = t161 * mrSges(6,1) + mrSges(6,2) * t181;
t100 = mrSges(5,1) * t184 + mrSges(5,2) * t181;
t113 = -pkin(3) * t161 - qJ(4) * t158 - pkin(2);
t78 = t154 * t113 + t156 * t192;
t180 = t204 * qJ(4) ^ 2;
t179 = Ifges(4,5) * t158 + Ifges(4,6) * t161;
t5 = Ifges(7,5) * t29 + Ifges(7,6) * t28 - Ifges(7,3) * t96;
t11 = t96 * qJ(5) + t14;
t44 = Ifges(7,5) * t91 + Ifges(7,6) * t90 + Ifges(7,3) * t161;
t18 = Ifges(6,5) * t63 + Ifges(6,6) * t96 + Ifges(6,3) * t62;
t21 = Ifges(5,4) * t63 - Ifges(5,2) * t62 + Ifges(5,6) * t96;
t178 = t18 / 0.2e1 - t21 / 0.2e1;
t22 = Ifges(6,1) * t63 + Ifges(6,4) * t96 + Ifges(6,5) * t62;
t23 = Ifges(5,1) * t63 - Ifges(5,4) * t62 + Ifges(5,5) * t96;
t177 = t22 / 0.2e1 + t23 / 0.2e1;
t81 = -Ifges(6,6) * t161 + (Ifges(6,3) * t154 + t187) * t158;
t84 = -Ifges(5,6) * t161 + (-Ifges(5,2) * t154 + t189) * t158;
t176 = t81 / 0.2e1 - t84 / 0.2e1;
t85 = -Ifges(6,4) * t161 + (Ifges(6,1) * t156 + t188) * t158;
t86 = -Ifges(5,5) * t161 + (Ifges(5,1) * t156 - t190) * t158;
t175 = t85 / 0.2e1 + t86 / 0.2e1;
t174 = Ifges(3,5) * t183 + Ifges(3,6) * t182 + Ifges(3,3) * t185;
t118 = -Ifges(6,3) * t156 + t188;
t121 = Ifges(5,2) * t156 + t190;
t172 = t118 / 0.2e1 - t121 / 0.2e1;
t122 = Ifges(6,1) * t154 - t187;
t123 = Ifges(5,1) * t154 + t189;
t171 = t122 / 0.2e1 + t123 / 0.2e1;
t31 = t62 * mrSges(5,1) + t63 * mrSges(5,2);
t41 = -t96 * mrSges(6,1) + t63 * mrSges(6,2);
t30 = t62 * mrSges(6,1) - t63 * mrSges(6,3);
t170 = qJ(5) * t154 + pkin(3);
t13 = -t154 * t36 + t156 * t35;
t42 = -t158 * t88 + t161 * t165;
t134 = t154 * t192;
t77 = t113 * t156 - t134;
t168 = -t59 / 0.2e1 + (Ifges(5,5) + Ifges(6,4)) * t154 / 0.2e1 + (Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1) * t156;
t10 = -t28 * mrSges(7,1) + t29 * mrSges(7,2);
t49 = -t90 * mrSges(7,1) + t91 * mrSges(7,2);
t68 = -qJ(5) * t161 + t78;
t58 = mrSges(7,1) * t167 + t107 * mrSges(7,2);
t166 = Ifges(4,5) * t97 - Ifges(4,6) * t96 - Ifges(4,3) * t182;
t37 = pkin(3) * t182 - t42;
t15 = t62 * pkin(4) - t63 * qJ(5) + t37;
t164 = pkin(9) ^ 2;
t153 = t161 ^ 2;
t152 = t158 ^ 2;
t149 = t161 * pkin(4);
t148 = t152 * t164;
t142 = t154 * mrSges(5,2);
t128 = mrSges(6,1) * t184;
t127 = qJ(5) * t181;
t126 = Ifges(4,1) * t158 + Ifges(4,4) * t161;
t125 = Ifges(4,4) * t158 + Ifges(4,2) * t161;
t124 = -mrSges(4,1) * t161 + mrSges(4,2) * t158;
t117 = t191 * t156;
t116 = -t156 * mrSges(5,1) + t142;
t115 = -t156 * mrSges(6,1) - t186;
t114 = t191 * t154;
t112 = -pkin(4) * t156 - t170;
t111 = -mrSges(6,2) * t184 - mrSges(6,3) * t161;
t109 = -mrSges(5,1) * t161 - mrSges(5,3) * t181;
t108 = mrSges(5,2) * t161 - mrSges(5,3) * t184;
t99 = -mrSges(6,3) * t181 + t128;
t98 = t196 * t156 + t170;
t89 = -t127 + (pkin(4) * t154 + pkin(9)) * t158;
t83 = -Ifges(6,2) * t161 + (Ifges(6,4) * t156 + Ifges(6,6) * t154) * t158;
t82 = -Ifges(5,3) * t161 + (Ifges(5,5) * t156 - Ifges(5,6) * t154) * t158;
t73 = mrSges(7,1) * t161 - mrSges(7,3) * t91;
t72 = -mrSges(7,2) * t161 + mrSges(7,3) * t90;
t71 = -mrSges(4,1) * t182 - mrSges(4,3) * t97;
t70 = mrSges(4,2) * t182 - mrSges(4,3) * t96;
t69 = t149 - t77;
t66 = -t127 - (-t196 * t154 - pkin(9)) * t158;
t65 = t114 * t157 + t117 * t160;
t64 = t114 * t160 - t117 * t157;
t61 = Ifges(7,1) * t107 - Ifges(7,4) * t167;
t60 = Ifges(7,4) * t107 - Ifges(7,2) * t167;
t52 = pkin(10) * t184 + t68;
t51 = mrSges(4,1) * t96 + mrSges(4,2) * t97;
t50 = pkin(5) * t161 + t134 + t149 + (-pkin(10) * t158 - t113) * t156;
t48 = Ifges(4,1) * t97 - Ifges(4,4) * t96 - Ifges(4,5) * t182;
t47 = Ifges(4,4) * t97 - Ifges(4,2) * t96 - Ifges(4,6) * t182;
t46 = Ifges(7,1) * t91 + Ifges(7,4) * t90 + Ifges(7,5) * t161;
t45 = Ifges(7,4) * t91 + Ifges(7,2) * t90 + Ifges(7,6) * t161;
t40 = mrSges(5,1) * t96 - mrSges(5,3) * t63;
t39 = -mrSges(5,2) * t96 - mrSges(5,3) * t62;
t38 = -mrSges(6,2) * t62 + mrSges(6,3) * t96;
t25 = t157 * t50 + t160 * t52;
t24 = -t157 * t52 + t160 * t50;
t20 = Ifges(6,4) * t63 + Ifges(6,2) * t96 + Ifges(6,6) * t62;
t19 = Ifges(5,5) * t63 - Ifges(5,6) * t62 + Ifges(5,3) * t96;
t17 = -mrSges(7,1) * t96 - mrSges(7,3) * t29;
t16 = mrSges(7,2) * t96 + mrSges(7,3) * t28;
t12 = -pkin(4) * t96 - t13;
t8 = pkin(5) * t62 + t15;
t7 = Ifges(7,1) * t29 + Ifges(7,4) * t28 - Ifges(7,5) * t96;
t6 = Ifges(7,4) * t29 + Ifges(7,2) * t28 - Ifges(7,6) * t96;
t4 = pkin(10) * t62 + t11;
t3 = -pkin(10) * t63 - t196 * t96 - t13;
t2 = t157 * t3 + t160 * t4;
t1 = -t157 * t4 + t160 * t3;
t9 = [t97 * t48 + 0.2e1 * t87 * t51 + 0.2e1 * t43 * t70 + 0.2e1 * t42 * t71 + 0.2e1 * t37 * t31 + 0.2e1 * t11 * t38 + 0.2e1 * t14 * t39 + 0.2e1 * t13 * t40 + 0.2e1 * t12 * t41 + t28 * t6 + t29 * t7 + 0.2e1 * t15 * t30 + 0.2e1 * t2 * t16 + 0.2e1 * t1 * t17 - 0.2e1 * t8 * t10 + m(4) * (t42 ^ 2 + t43 ^ 2 + t87 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2 + t37 ^ 2) + m(6) * (t11 ^ 2 + t12 ^ 2 + t15 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t8 ^ 2) + (t19 + t20 - t47 - t5) * t96 + (t22 + t23) * t63 + (t18 - t21) * t62 - 0.2e1 * t201 * pkin(1) * (-mrSges(3,1) * t162 + mrSges(3,2) * t159) + Ifges(2,3) - t166 * t182 + (Ifges(3,6) * t185 + (Ifges(3,4) * t159 + Ifges(3,2) * t162) * t155) * t182 + t185 * t174 + 0.2e1 * t102 * (-t185 * mrSges(3,2) + mrSges(3,3) * t182) + 0.2e1 * t101 * (t185 * mrSges(3,1) - mrSges(3,3) * t183) + (Ifges(3,5) * t185 + (Ifges(3,1) * t159 + Ifges(3,4) * t162) * t155) * t183 + m(3) * (pkin(1) ^ 2 * t201 + t101 ^ 2 + t102 ^ 2); m(6) * (t11 * t68 + t12 * t69 + t15 * t89) + m(7) * (t1 * t24 + t2 * t25 + t66 * t8) + (-t125 / 0.2e1 + t82 / 0.2e1 + t83 / 0.2e1 - t44 / 0.2e1) * t96 + t12 * t110 + t11 * t111 + t87 * t124 + t97 * t126 / 0.2e1 + t15 * t99 + t37 * t100 + t101 * mrSges(3,1) - t102 * mrSges(3,2) + t14 * t108 + t13 * t109 + t89 * t30 - t66 * t10 + t68 * t38 + t69 * t41 + t2 * t72 + t1 * t73 + t77 * t40 + t78 * t39 + t174 - t8 * t49 - pkin(2) * t51 + t24 * t17 + t25 * t16 + (pkin(9) * t70 + t43 * mrSges(4,3) + t47 / 0.2e1 - t19 / 0.2e1 - t20 / 0.2e1 + t5 / 0.2e1) * t161 + t175 * t63 + t176 * t62 + (-t42 * mrSges(4,3) + t48 / 0.2e1 + t177 * t156 + t178 * t154 + (-t71 + t31) * pkin(9)) * t158 - t179 * t182 / 0.2e1 + t7 * t197 + t6 * t198 + t46 * t199 + t45 * t200 + m(5) * (t13 * t77 + t14 * t78 + t37 * t193) + m(4) * (-pkin(2) * t87 + (-t158 * t42 + t161 * t43) * pkin(9)); -0.2e1 * pkin(2) * t124 + 0.2e1 * t78 * t108 + 0.2e1 * t77 * t109 + 0.2e1 * t69 * t110 + 0.2e1 * t68 * t111 + 0.2e1 * t24 * t73 + 0.2e1 * t25 * t72 + t90 * t45 + t91 * t46 - 0.2e1 * t66 * t49 + 0.2e1 * t89 * t99 + Ifges(3,3) + (t152 + t153) * mrSges(4,3) * t203 + (t125 - t82 - t83 + t44) * t161 + m(4) * (pkin(2) ^ 2 + t153 * t164 + t148) + m(5) * (t77 ^ 2 + t78 ^ 2 + t148) + m(6) * (t68 ^ 2 + t69 ^ 2 + t89 ^ 2) + m(7) * (t24 ^ 2 + t25 ^ 2 + t66 ^ 2) + (t100 * t203 + t126 + (t85 + t86) * t156 + (t81 - t84) * t154) * t158; (-t1 * t107 - t167 * t2) * mrSges(7,3) + t112 * t30 + t15 * t115 + t37 * t116 + t98 * t10 + t6 * t195 + t7 * t194 + t64 * t17 + t65 * t16 + t60 * t200 + t61 * t199 - t8 * t58 + t42 * mrSges(4,1) - t43 * mrSges(4,2) - pkin(3) * t31 + (-t13 * mrSges(5,3) + t12 * mrSges(6,2) + (-t40 + t41) * qJ(4) + t177) * t154 + (t14 * mrSges(5,3) + t11 * mrSges(6,2) + (t38 + t39) * qJ(4) - t178) * t156 + m(7) * (t1 * t64 + t2 * t65 - t8 * t98) + t168 * t96 + t171 * t63 + t172 * t62 + m(5) * (-pkin(3) * t37 + (-t13 * t154 + t14 * t156) * qJ(4)) + m(6) * (t112 * t15 + (t11 * t156 + t12 * t154) * qJ(4)) + t166; t112 * t99 + t89 * t115 + t98 * t49 - pkin(3) * t100 + t45 * t195 + t46 * t194 + t60 * t198 + t61 * t197 - t66 * t58 + t65 * t72 + t64 * t73 + (-mrSges(4,1) + t116) * t193 + (-t107 * t24 - t167 * t25) * mrSges(7,3) + (t78 * mrSges(5,3) + t68 * mrSges(6,2) + t171 * t158 + (t108 + t111) * qJ(4) - t176) * t156 + (t69 * mrSges(6,2) - t77 * mrSges(5,3) + t172 * t158 + (-t109 + t110) * qJ(4) + t175) * t154 + m(5) * (-pkin(3) * t193 + (-t154 * t77 + t156 * t78) * qJ(4)) + m(6) * (t112 * t89 + (t154 * t69 + t156 * t68) * qJ(4)) + m(7) * (t24 * t64 + t25 * t65 - t66 * t98) + (-pkin(9) * mrSges(4,2) - t168) * t161 + t179; -0.2e1 * pkin(3) * t116 - t167 * t60 + t107 * t61 + 0.2e1 * t112 * t115 + 0.2e1 * t98 * t58 + Ifges(4,3) + (-t118 + t121) * t156 + (t122 + t123) * t154 + 0.2e1 * (-t107 * t64 - t167 * t65) * mrSges(7,3) + m(7) * (t64 ^ 2 + t65 ^ 2 + t98 ^ 2) + m(6) * (t112 ^ 2 + t180) + m(5) * (pkin(3) ^ 2 + t180) + 0.2e1 * (mrSges(6,2) + mrSges(5,3)) * qJ(4) * t204; m(5) * t37 + m(6) * t15 + m(7) * t8 - t10 + t30 + t31; t128 + (m(5) * pkin(9) - mrSges(6,3) * t156) * t158 + m(6) * t89 + m(7) * t66 - t49 + t100; -m(5) * pkin(3) - t186 + t142 + (-mrSges(6,1) - mrSges(5,1)) * t156 + m(6) * t112 - m(7) * t98 - t58; m(5) + m(6) + m(7); t157 * t16 + t160 * t17 + m(7) * (t1 * t160 + t157 * t2) + m(6) * t12 + t41; t157 * t72 + t160 * t73 + m(7) * (t157 * t25 + t160 * t24) + m(6) * t69 + t110; m(7) * (t157 * t65 + t160 * t64) + (m(6) * qJ(4) + mrSges(6,2)) * t154 + (-t107 * t160 - t157 * t167) * mrSges(7,3); 0; m(6) + m(7) * (t157 ^ 2 + t160 ^ 2); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t5; mrSges(7,1) * t24 - mrSges(7,2) * t25 + t44; mrSges(7,1) * t64 - mrSges(7,2) * t65 + t59; 0; mrSges(7,1) * t160 - mrSges(7,2) * t157; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t9(1) t9(2) t9(4) t9(7) t9(11) t9(16); t9(2) t9(3) t9(5) t9(8) t9(12) t9(17); t9(4) t9(5) t9(6) t9(9) t9(13) t9(18); t9(7) t9(8) t9(9) t9(10) t9(14) t9(19); t9(11) t9(12) t9(13) t9(14) t9(15) t9(20); t9(16) t9(17) t9(18) t9(19) t9(20) t9(21);];
Mq  = res;
