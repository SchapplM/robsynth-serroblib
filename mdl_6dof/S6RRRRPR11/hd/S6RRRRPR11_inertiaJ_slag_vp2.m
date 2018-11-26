% Calculate joint inertia matrix for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:21
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR11_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR11_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR11_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR11_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:20:32
% EndTime: 2018-11-23 18:20:34
% DurationCPUTime: 2.10s
% Computational Cost: add. (5280->477), mult. (11945->683), div. (0->0), fcn. (13433->12), ass. (0->186)
t247 = 2 * pkin(9);
t190 = cos(pkin(6));
t193 = sin(qJ(3));
t197 = cos(qJ(3));
t188 = sin(pkin(6));
t194 = sin(qJ(2));
t217 = t188 * t194;
t136 = -t190 * t197 + t193 * t217;
t137 = t190 * t193 + t197 * t217;
t192 = sin(qJ(4));
t196 = cos(qJ(4));
t198 = cos(qJ(2));
t216 = t188 * t198;
t103 = -t137 * t192 - t196 * t216;
t104 = t137 * t196 - t192 * t216;
t187 = sin(pkin(12));
t189 = cos(pkin(12));
t58 = t103 * t189 - t104 * t187;
t59 = t103 * t187 + t104 * t189;
t20 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + Ifges(6,3) * t136;
t46 = Ifges(5,5) * t104 + Ifges(5,6) * t103 + Ifges(5,3) * t136;
t246 = t20 + t46;
t191 = sin(qJ(6));
t195 = cos(qJ(6));
t29 = -t191 * t59 + t195 * t58;
t245 = t29 / 0.2e1;
t30 = t191 * t58 + t195 * t59;
t244 = t30 / 0.2e1;
t48 = Ifges(5,1) * t104 + Ifges(5,4) * t103 + Ifges(5,5) * t136;
t243 = t48 / 0.2e1;
t242 = t58 / 0.2e1;
t241 = t59 / 0.2e1;
t145 = t187 * t196 + t189 * t192;
t127 = t145 * t193;
t144 = -t187 * t192 + t189 * t196;
t128 = t144 * t193;
t79 = -t127 * t195 - t128 * t191;
t240 = t79 / 0.2e1;
t80 = -t127 * t191 + t128 * t195;
t239 = t80 / 0.2e1;
t92 = t144 * t195 - t145 * t191;
t238 = t92 / 0.2e1;
t93 = t144 * t191 + t145 * t195;
t237 = t93 / 0.2e1;
t223 = Ifges(5,4) * t196;
t125 = -Ifges(5,6) * t197 + (-Ifges(5,2) * t192 + t223) * t193;
t236 = t125 / 0.2e1;
t224 = Ifges(5,4) * t192;
t126 = -Ifges(5,5) * t197 + (Ifges(5,1) * t196 - t224) * t193;
t235 = t126 / 0.2e1;
t234 = -t127 / 0.2e1;
t233 = t128 / 0.2e1;
t232 = t144 / 0.2e1;
t231 = t145 / 0.2e1;
t162 = Ifges(5,1) * t192 + t223;
t230 = t162 / 0.2e1;
t229 = pkin(1) * t198;
t228 = pkin(4) * t187;
t227 = pkin(9) * t197;
t182 = t193 * pkin(9);
t226 = -qJ(5) - pkin(10);
t166 = pkin(8) * t217;
t121 = t166 + (-pkin(2) - t229) * t190;
t64 = pkin(3) * t136 - pkin(10) * t137 + t121;
t139 = t190 * t194 * pkin(1) + pkin(8) * t216;
t122 = pkin(9) * t190 + t139;
t123 = (-pkin(2) * t198 - pkin(9) * t194 - pkin(1)) * t188;
t72 = t197 * t122 + t193 * t123;
t66 = -pkin(10) * t216 + t72;
t32 = -t192 * t66 + t196 * t64;
t17 = pkin(4) * t136 - qJ(5) * t104 + t32;
t33 = t192 * t64 + t196 * t66;
t25 = qJ(5) * t103 + t33;
t7 = t187 * t17 + t189 * t25;
t225 = Ifges(7,5) * t80 + Ifges(7,6) * t79;
t50 = Ifges(7,5) * t93 + Ifges(7,6) * t92;
t154 = -pkin(3) * t197 - pkin(10) * t193 - pkin(2);
t116 = t192 * t154 + t196 * t227;
t215 = t192 * t193;
t105 = -qJ(5) * t215 + t116;
t147 = t196 * t154;
t214 = t193 * t196;
t94 = -qJ(5) * t214 + t147 + (-pkin(9) * t192 - pkin(4)) * t197;
t54 = t189 * t105 + t187 * t94;
t222 = Ifges(7,3) * t197;
t174 = pkin(4) * t189 + pkin(5);
t134 = t174 * t195 - t191 * t228;
t221 = t134 * mrSges(7,1);
t135 = t174 * t191 + t195 * t228;
t220 = t135 * mrSges(7,2);
t138 = t190 * t229 - t166;
t219 = t138 * mrSges(3,1);
t218 = t139 * mrSges(3,2);
t213 = Ifges(6,5) * t128 - Ifges(6,6) * t127;
t212 = -Ifges(4,5) * t137 + Ifges(4,6) * t136;
t98 = Ifges(6,5) * t145 + Ifges(6,6) * t144;
t155 = t226 * t192;
t158 = t226 * t196;
t107 = t187 * t155 - t189 * t158;
t159 = Ifges(5,5) * t192 + Ifges(5,6) * t196;
t211 = Ifges(4,5) * t193 + Ifges(4,6) * t197;
t153 = pkin(4) * t215 + t182;
t210 = t192 ^ 2 + t196 ^ 2;
t209 = -Ifges(7,3) - Ifges(5,3) - Ifges(6,3);
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t136;
t208 = Ifges(3,5) * t217 + Ifges(3,6) * t216 + Ifges(3,3) * t190;
t175 = -pkin(4) * t196 - pkin(3);
t31 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t11 = -t29 * mrSges(7,1) + t30 * mrSges(7,2);
t37 = -t79 * mrSges(7,1) + t80 * mrSges(7,2);
t49 = -t92 * mrSges(7,1) + t93 * mrSges(7,2);
t6 = t189 * t17 - t187 * t25;
t85 = t127 * mrSges(6,1) + t128 * mrSges(6,2);
t97 = -t144 * mrSges(6,1) + t145 * mrSges(6,2);
t53 = -t105 * t187 + t189 * t94;
t71 = -t193 * t122 + t123 * t197;
t106 = t189 * t155 + t158 * t187;
t207 = t50 / 0.2e1 + t98 / 0.2e1 + t159 / 0.2e1;
t206 = Ifges(5,5) * t214 - Ifges(5,6) * t215;
t65 = pkin(3) * t216 - t71;
t205 = mrSges(5,1) * t192 + mrSges(5,2) * t196;
t39 = -pkin(5) * t197 - pkin(11) * t128 + t53;
t42 = -pkin(11) * t127 + t54;
t13 = -t191 * t42 + t195 * t39;
t14 = t191 * t39 + t195 * t42;
t204 = t13 * mrSges(7,1) - t14 * mrSges(7,2) + t225;
t83 = -pkin(11) * t145 + t106;
t84 = pkin(11) * t144 + t107;
t40 = -t191 * t84 + t195 * t83;
t41 = t191 * t83 + t195 * t84;
t203 = t40 * mrSges(7,1) - t41 * mrSges(7,2) + t50;
t4 = pkin(5) * t136 - pkin(11) * t59 + t6;
t5 = pkin(11) * t58 + t7;
t2 = -t191 * t5 + t195 * t4;
t3 = t191 * t4 + t195 * t5;
t202 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t45 = -pkin(4) * t103 + t65;
t200 = pkin(9) ^ 2;
t186 = t197 ^ 2;
t184 = t193 ^ 2;
t181 = t184 * t200;
t163 = Ifges(4,1) * t193 + Ifges(4,4) * t197;
t161 = Ifges(4,4) * t193 + Ifges(4,2) * t197;
t160 = Ifges(5,2) * t196 + t224;
t157 = -mrSges(4,1) * t197 + mrSges(4,2) * t193;
t156 = -mrSges(5,1) * t196 + mrSges(5,2) * t192;
t152 = -mrSges(5,1) * t197 - mrSges(5,3) * t214;
t151 = mrSges(5,2) * t197 - mrSges(5,3) * t215;
t143 = t205 * t193;
t124 = -Ifges(5,3) * t197 + t206;
t117 = -pkin(5) * t144 + t175;
t115 = -t192 * t227 + t147;
t111 = -mrSges(6,1) * t197 - mrSges(6,3) * t128;
t110 = mrSges(6,2) * t197 - mrSges(6,3) * t127;
t109 = -mrSges(4,1) * t216 - mrSges(4,3) * t137;
t108 = mrSges(4,2) * t216 - mrSges(4,3) * t136;
t100 = Ifges(6,1) * t145 + Ifges(6,4) * t144;
t99 = Ifges(6,4) * t145 + Ifges(6,2) * t144;
t95 = pkin(5) * t127 + t153;
t86 = mrSges(4,1) * t136 + mrSges(4,2) * t137;
t82 = Ifges(4,1) * t137 - Ifges(4,4) * t136 - Ifges(4,5) * t216;
t81 = Ifges(4,4) * t137 - Ifges(4,2) * t136 - Ifges(4,6) * t216;
t78 = Ifges(6,1) * t128 - Ifges(6,4) * t127 - Ifges(6,5) * t197;
t77 = Ifges(6,4) * t128 - Ifges(6,2) * t127 - Ifges(6,6) * t197;
t76 = -Ifges(6,3) * t197 + t213;
t70 = mrSges(5,1) * t136 - mrSges(5,3) * t104;
t69 = -mrSges(5,2) * t136 + mrSges(5,3) * t103;
t68 = -mrSges(7,1) * t197 - mrSges(7,3) * t80;
t67 = mrSges(7,2) * t197 + mrSges(7,3) * t79;
t60 = -mrSges(5,1) * t103 + mrSges(5,2) * t104;
t52 = Ifges(7,1) * t93 + Ifges(7,4) * t92;
t51 = Ifges(7,4) * t93 + Ifges(7,2) * t92;
t47 = Ifges(5,4) * t104 + Ifges(5,2) * t103 + Ifges(5,6) * t136;
t44 = mrSges(6,1) * t136 - mrSges(6,3) * t59;
t43 = -mrSges(6,2) * t136 + mrSges(6,3) * t58;
t36 = Ifges(7,1) * t80 + Ifges(7,4) * t79 - Ifges(7,5) * t197;
t35 = Ifges(7,4) * t80 + Ifges(7,2) * t79 - Ifges(7,6) * t197;
t34 = -t222 + t225;
t23 = -pkin(5) * t58 + t45;
t22 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t136;
t21 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t136;
t19 = mrSges(7,1) * t136 - mrSges(7,3) * t30;
t18 = -mrSges(7,2) * t136 + mrSges(7,3) * t29;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t136;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t136;
t1 = [m(7) * (t2 ^ 2 + t23 ^ 2 + t3 ^ 2) + m(6) * (t45 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t65 ^ 2) + m(4) * (t121 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(3) * (pkin(1) ^ 2 * t188 ^ 2 + t138 ^ 2 + t139 ^ 2) + ((-0.2e1 * t138 * mrSges(3,3) + Ifges(3,5) * t190 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t194) * t188) * t194 + (0.2e1 * t139 * mrSges(3,3) + Ifges(3,6) * t190 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t194 + (Ifges(4,3) + Ifges(3,2)) * t198) * t188 + t212) * t198) * t188 + (t208 - 0.2e1 * t218 + 0.2e1 * t219) * t190 + Ifges(2,3) + (t8 - t81 + t246) * t136 + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 + 0.2e1 * t23 * t11 + t29 * t9 + t30 * t10 + 0.2e1 * t7 * t43 + 0.2e1 * t6 * t44 + 0.2e1 * t45 * t31 + t58 * t21 + t59 * t22 + 0.2e1 * t65 * t60 + 0.2e1 * t33 * t69 + 0.2e1 * t32 * t70 + t103 * t47 + t104 * t48 + 0.2e1 * t72 * t108 + 0.2e1 * t71 * t109 + 0.2e1 * t121 * t86 + t137 * t82; -t218 + t219 + t208 + m(4) * (-pkin(2) * t121 + (-t193 * t71 + t197 * t72) * pkin(9)) + (t72 * mrSges(4,3) + pkin(9) * t108 - t8 / 0.2e1 - t20 / 0.2e1 - t46 / 0.2e1 + t81 / 0.2e1) * t197 + (t34 / 0.2e1 + t76 / 0.2e1 + t124 / 0.2e1 - t161 / 0.2e1) * t136 + m(7) * (t13 * t2 + t14 * t3 + t23 * t95) + m(6) * (t153 * t45 + t53 * t6 + t54 * t7) - t211 * t216 / 0.2e1 + m(5) * (t115 * t32 + t116 * t33 + t182 * t65) + t22 * t233 + t21 * t234 + t104 * t235 + t103 * t236 + t10 * t239 + t9 * t240 + t78 * t241 + t77 * t242 + t36 * t244 + t35 * t245 + (-t71 * mrSges(4,3) - t192 * t47 / 0.2e1 + t196 * t243 + t82 / 0.2e1 + (t60 - t109) * pkin(9)) * t193 + t14 * t18 + t13 * t19 + t23 * t37 + t53 * t44 + t54 * t43 + t3 * t67 + t2 * t68 + t45 * t85 - pkin(2) * t86 + t95 * t11 + t7 * t110 + t6 * t111 + t115 * t70 + t116 * t69 + t65 * t143 + t33 * t151 + t32 * t152 + t153 * t31 + t121 * t157 + t137 * t163 / 0.2e1; -0.2e1 * pkin(2) * t157 + 0.2e1 * t54 * t110 + 0.2e1 * t53 * t111 + 0.2e1 * t115 * t152 + 0.2e1 * t116 * t151 - t127 * t77 + t128 * t78 + 0.2e1 * t13 * t68 + 0.2e1 * t14 * t67 + 0.2e1 * t153 * t85 + t79 * t35 + t80 * t36 + 0.2e1 * t95 * t37 + Ifges(3,3) + (t184 + t186) * mrSges(4,3) * t247 + (-t34 - t76 - t124 + t161) * t197 + (-t125 * t192 + t126 * t196 + t143 * t247 + t163) * t193 + m(7) * (t13 ^ 2 + t14 ^ 2 + t95 ^ 2) + m(6) * (t153 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t115 ^ 2 + t116 ^ 2 + t181) + m(4) * (pkin(2) ^ 2 + t186 * t200 + t181); (t33 * mrSges(5,3) + pkin(10) * t69 + t47 / 0.2e1) * t196 + (-t2 * t93 + t3 * t92) * mrSges(7,3) + (t144 * t7 - t145 * t6) * mrSges(6,3) + m(7) * (t117 * t23 + t2 * t40 + t3 * t41) + m(6) * (t106 * t6 + t107 * t7 + t175 * t45) + m(5) * (-pkin(3) * t65 + (-t192 * t32 + t196 * t33) * pkin(10)) - t212 + t207 * t136 - Ifges(4,3) * t216 + t104 * t230 + t22 * t231 + t21 * t232 + t10 * t237 + t9 * t238 + t100 * t241 + t99 * t242 + t52 * t244 + t51 * t245 + (-t32 * mrSges(5,3) - pkin(10) * t70 + t243) * t192 + t40 * t19 + t41 * t18 + t23 * t49 - pkin(3) * t60 + t71 * mrSges(4,1) - t72 * mrSges(4,2) + t45 * t97 + t106 * t44 + t107 * t43 + t117 * t11 + t65 * t156 + t103 * t160 / 0.2e1 + t175 * t31; t211 + m(7) * (t117 * t95 + t13 * t40 + t14 * t41) + m(6) * (t106 * t53 + t107 * t54 + t153 * t175) + (-pkin(9) * mrSges(4,2) - t207) * t197 + m(5) * (-pkin(3) * t182 + (-t115 * t192 + t116 * t196) * pkin(10)) + (-t115 * mrSges(5,3) - pkin(10) * t152 - t193 * t160 / 0.2e1 + t235) * t192 + (t116 * mrSges(5,3) + pkin(10) * t151 + t193 * t230 + t236) * t196 + (-mrSges(4,1) + t156) * t182 + t78 * t231 + t77 * t232 + t100 * t233 + t99 * t234 + t36 * t237 + t35 * t238 + t52 * t239 + t51 * t240 + t41 * t67 + t40 * t68 + t95 * t49 + t107 * t110 + t106 * t111 + t117 * t37 - pkin(3) * t143 + t153 * t97 + t175 * t85 + (t144 * t54 - t145 * t53) * mrSges(6,3) + (-t13 * t93 + t14 * t92) * mrSges(7,3); -0.2e1 * pkin(3) * t156 + t145 * t100 + 0.2e1 * t117 * t49 + t144 * t99 + t196 * t160 + t192 * t162 + 0.2e1 * t175 * t97 + t92 * t51 + t93 * t52 + Ifges(4,3) + m(7) * (t117 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t106 ^ 2 + t107 ^ 2 + t175 ^ 2) + m(5) * (pkin(10) ^ 2 * t210 + pkin(3) ^ 2) + 0.2e1 * (-t40 * t93 + t41 * t92) * mrSges(7,3) + 0.2e1 * (-t106 * t145 + t107 * t144) * mrSges(6,3) + 0.2e1 * t210 * pkin(10) * mrSges(5,3); t202 + (m(6) * (t187 * t7 + t189 * t6) + t189 * t44 + t187 * t43) * pkin(4) + m(7) * (t134 * t2 + t135 * t3) + t6 * mrSges(6,1) - t7 * mrSges(6,2) + t32 * mrSges(5,1) - t33 * mrSges(5,2) + t134 * t19 + t135 * t18 + t246; m(7) * (t13 * t134 + t135 * t14) + t134 * t68 + t135 * t67 + t53 * mrSges(6,1) - t54 * mrSges(6,2) - t116 * mrSges(5,2) + t115 * mrSges(5,1) + t209 * t197 + (m(6) * (t187 * t54 + t189 * t53) + t187 * t110 + t189 * t111) * pkin(4) + t204 + t206 + t213; m(7) * (t134 * t40 + t135 * t41) - t107 * mrSges(6,2) + t106 * mrSges(6,1) - t205 * pkin(10) + (-t134 * t93 + t135 * t92) * mrSges(7,3) + (m(6) * (t106 * t189 + t107 * t187) + (t144 * t187 - t145 * t189) * mrSges(6,3)) * pkin(4) + t203 + t159 + t98; 0.2e1 * t221 - 0.2e1 * t220 + m(7) * (t134 ^ 2 + t135 ^ 2) - t209 + (0.2e1 * mrSges(6,1) * t189 - 0.2e1 * mrSges(6,2) * t187 + m(6) * (t187 ^ 2 + t189 ^ 2) * pkin(4)) * pkin(4); m(6) * t45 + m(7) * t23 + t11 + t31; m(6) * t153 + m(7) * t95 + t37 + t85; m(6) * t175 + m(7) * t117 + t49 + t97; 0; m(6) + m(7); t202; t204 - t222; t203; Ifges(7,3) - t220 + t221; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
