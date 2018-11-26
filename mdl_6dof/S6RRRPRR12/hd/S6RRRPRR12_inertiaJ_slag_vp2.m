% Calculate joint inertia matrix for
% S6RRRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2018-11-23 18:00
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR12_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR12_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR12_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:59:54
% EndTime: 2018-11-23 17:59:56
% DurationCPUTime: 1.97s
% Computational Cost: add. (5080->457), mult. (11524->653), div. (0->0), fcn. (13006->12), ass. (0->177)
t237 = 2 * pkin(9);
t186 = sin(qJ(6));
t190 = cos(qJ(6));
t185 = cos(pkin(6));
t188 = sin(qJ(3));
t192 = cos(qJ(3));
t183 = sin(pkin(6));
t189 = sin(qJ(2));
t209 = t183 * t189;
t133 = t185 * t188 + t192 * t209;
t182 = sin(pkin(12));
t184 = cos(pkin(12));
t193 = cos(qJ(2));
t208 = t183 * t193;
t102 = -t133 * t182 - t184 * t208;
t103 = t133 * t184 - t182 * t208;
t187 = sin(qJ(5));
t191 = cos(qJ(5));
t58 = t102 * t191 - t103 * t187;
t59 = t102 * t187 + t103 * t191;
t29 = -t186 * t59 + t190 * t58;
t236 = t29 / 0.2e1;
t30 = t186 * t58 + t190 * t59;
t235 = t30 / 0.2e1;
t132 = -t185 * t192 + t188 * t209;
t48 = Ifges(5,1) * t103 + Ifges(5,4) * t102 + Ifges(5,5) * t132;
t234 = t48 / 0.2e1;
t233 = t58 / 0.2e1;
t232 = t59 / 0.2e1;
t146 = t182 * t191 + t184 * t187;
t126 = t146 * t188;
t145 = -t182 * t187 + t184 * t191;
t127 = t145 * t188;
t79 = -t126 * t190 - t127 * t186;
t231 = t79 / 0.2e1;
t80 = -t126 * t186 + t127 * t190;
t230 = t80 / 0.2e1;
t92 = t145 * t190 - t146 * t186;
t229 = t92 / 0.2e1;
t93 = t145 * t186 + t146 * t190;
t228 = t93 / 0.2e1;
t214 = Ifges(5,4) * t184;
t121 = -Ifges(5,6) * t192 + (-Ifges(5,2) * t182 + t214) * t188;
t227 = t121 / 0.2e1;
t215 = Ifges(5,4) * t182;
t122 = -Ifges(5,5) * t192 + (Ifges(5,1) * t184 - t215) * t188;
t226 = t122 / 0.2e1;
t225 = -t126 / 0.2e1;
t224 = t127 / 0.2e1;
t223 = t145 / 0.2e1;
t222 = t146 / 0.2e1;
t156 = Ifges(5,1) * t182 + t214;
t221 = t156 / 0.2e1;
t220 = pkin(1) * t193;
t219 = pkin(9) * t192;
t177 = t188 * pkin(9);
t218 = -Ifges(7,3) - Ifges(6,3);
t217 = pkin(10) + qJ(4);
t165 = pkin(8) * t209;
t123 = t165 + (-pkin(2) - t220) * t185;
t64 = pkin(3) * t132 - qJ(4) * t133 + t123;
t136 = t185 * t189 * pkin(1) + pkin(8) * t208;
t124 = pkin(9) * t185 + t136;
t125 = (-pkin(2) * t193 - pkin(9) * t189 - pkin(1)) * t183;
t72 = t192 * t124 + t188 * t125;
t65 = -qJ(4) * t208 + t72;
t32 = -t182 * t65 + t184 * t64;
t17 = pkin(4) * t132 - pkin(10) * t103 + t32;
t33 = t182 * t64 + t184 * t65;
t25 = pkin(10) * t102 + t33;
t7 = t187 * t17 + t191 * t25;
t216 = Ifges(7,5) * t80 + Ifges(7,6) * t79;
t50 = Ifges(7,5) * t93 + Ifges(7,6) * t92;
t150 = -pkin(3) * t192 - qJ(4) * t188 - pkin(2);
t115 = t182 * t150 + t184 * t219;
t210 = t182 * t188;
t104 = -pkin(10) * t210 + t115;
t141 = t184 * t150;
t207 = t184 * t188;
t94 = -pkin(10) * t207 + t141 + (-pkin(9) * t182 - pkin(4)) * t192;
t54 = t191 * t104 + t187 * t94;
t213 = Ifges(7,3) * t192;
t135 = t185 * t220 - t165;
t212 = t135 * mrSges(3,1);
t211 = t136 * mrSges(3,2);
t206 = Ifges(6,5) * t127 - Ifges(6,6) * t126;
t205 = -Ifges(4,5) * t133 + Ifges(4,6) * t132;
t99 = Ifges(6,5) * t146 + Ifges(6,6) * t145;
t151 = t217 * t182;
t153 = t217 * t184;
t106 = -t187 * t151 + t191 * t153;
t134 = mrSges(5,1) * t210 + mrSges(5,2) * t207;
t204 = Ifges(4,5) * t188 + Ifges(4,6) * t192;
t149 = pkin(4) * t210 + t177;
t203 = t182 ^ 2 + t184 ^ 2;
t8 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t132;
t20 = Ifges(6,5) * t59 + Ifges(6,6) * t58 + Ifges(6,3) * t132;
t202 = Ifges(3,5) * t209 + Ifges(3,6) * t208 + Ifges(3,3) * t185;
t171 = -pkin(4) * t184 - pkin(3);
t31 = -t58 * mrSges(6,1) + t59 * mrSges(6,2);
t11 = -t29 * mrSges(7,1) + t30 * mrSges(7,2);
t39 = -t79 * mrSges(7,1) + t80 * mrSges(7,2);
t49 = -t92 * mrSges(7,1) + t93 * mrSges(7,2);
t60 = -t102 * mrSges(5,1) + t103 * mrSges(5,2);
t6 = t191 * t17 - t187 * t25;
t152 = -t184 * mrSges(5,1) + t182 * mrSges(5,2);
t85 = t126 * mrSges(6,1) + t127 * mrSges(6,2);
t98 = -t145 * mrSges(6,1) + t146 * mrSges(6,2);
t53 = -t104 * t187 + t191 * t94;
t71 = -t188 * t124 + t125 * t192;
t105 = -t191 * t151 - t153 * t187;
t201 = t50 / 0.2e1 + t99 / 0.2e1 + Ifges(5,5) * t182 / 0.2e1 + Ifges(5,6) * t184 / 0.2e1;
t66 = pkin(3) * t208 - t71;
t38 = -pkin(5) * t192 - pkin(11) * t127 + t53;
t42 = -pkin(11) * t126 + t54;
t13 = -t186 * t42 + t190 * t38;
t14 = t186 * t38 + t190 * t42;
t200 = t13 * mrSges(7,1) - t14 * mrSges(7,2) + t216;
t83 = -pkin(11) * t146 + t105;
t84 = pkin(11) * t145 + t106;
t40 = -t186 * t84 + t190 * t83;
t41 = t186 * t83 + t190 * t84;
t199 = t40 * mrSges(7,1) - t41 * mrSges(7,2) + t50;
t4 = pkin(5) * t132 - pkin(11) * t59 + t6;
t5 = pkin(11) * t58 + t7;
t2 = -t186 * t5 + t190 * t4;
t3 = t186 * t4 + t190 * t5;
t198 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t197 = (mrSges(7,1) * t190 - mrSges(7,2) * t186) * pkin(5);
t45 = -pkin(4) * t102 + t66;
t195 = pkin(9) ^ 2;
t181 = t192 ^ 2;
t180 = t188 ^ 2;
t176 = t180 * t195;
t159 = Ifges(4,1) * t188 + Ifges(4,4) * t192;
t158 = Ifges(4,4) * t188 + Ifges(4,2) * t192;
t157 = -mrSges(4,1) * t192 + mrSges(4,2) * t188;
t155 = Ifges(5,2) * t184 + t215;
t148 = -mrSges(5,1) * t192 - mrSges(5,3) * t207;
t147 = mrSges(5,2) * t192 - mrSges(5,3) * t210;
t120 = -Ifges(5,3) * t192 + (Ifges(5,5) * t184 - Ifges(5,6) * t182) * t188;
t116 = -pkin(5) * t145 + t171;
t114 = -t182 * t219 + t141;
t110 = -mrSges(6,1) * t192 - mrSges(6,3) * t127;
t109 = mrSges(6,2) * t192 - mrSges(6,3) * t126;
t108 = -mrSges(4,1) * t208 - mrSges(4,3) * t133;
t107 = mrSges(4,2) * t208 - mrSges(4,3) * t132;
t101 = Ifges(6,1) * t146 + Ifges(6,4) * t145;
t100 = Ifges(6,4) * t146 + Ifges(6,2) * t145;
t96 = pkin(5) * t126 + t149;
t88 = mrSges(4,1) * t132 + mrSges(4,2) * t133;
t82 = Ifges(4,1) * t133 - Ifges(4,4) * t132 - Ifges(4,5) * t208;
t81 = Ifges(4,4) * t133 - Ifges(4,2) * t132 - Ifges(4,6) * t208;
t78 = Ifges(6,1) * t127 - Ifges(6,4) * t126 - Ifges(6,5) * t192;
t77 = Ifges(6,4) * t127 - Ifges(6,2) * t126 - Ifges(6,6) * t192;
t76 = -Ifges(6,3) * t192 + t206;
t70 = mrSges(5,1) * t132 - mrSges(5,3) * t103;
t69 = -mrSges(5,2) * t132 + mrSges(5,3) * t102;
t68 = -mrSges(7,1) * t192 - mrSges(7,3) * t80;
t67 = mrSges(7,2) * t192 + mrSges(7,3) * t79;
t52 = Ifges(7,1) * t93 + Ifges(7,4) * t92;
t51 = Ifges(7,4) * t93 + Ifges(7,2) * t92;
t47 = Ifges(5,4) * t103 + Ifges(5,2) * t102 + Ifges(5,6) * t132;
t46 = Ifges(5,5) * t103 + Ifges(5,6) * t102 + Ifges(5,3) * t132;
t44 = mrSges(6,1) * t132 - mrSges(6,3) * t59;
t43 = -mrSges(6,2) * t132 + mrSges(6,3) * t58;
t36 = Ifges(7,1) * t80 + Ifges(7,4) * t79 - Ifges(7,5) * t192;
t35 = Ifges(7,4) * t80 + Ifges(7,2) * t79 - Ifges(7,6) * t192;
t34 = -t213 + t216;
t23 = -pkin(5) * t58 + t45;
t22 = Ifges(6,1) * t59 + Ifges(6,4) * t58 + Ifges(6,5) * t132;
t21 = Ifges(6,4) * t59 + Ifges(6,2) * t58 + Ifges(6,6) * t132;
t19 = mrSges(7,1) * t132 - mrSges(7,3) * t30;
t18 = -mrSges(7,2) * t132 + mrSges(7,3) * t29;
t10 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t132;
t9 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t132;
t1 = [m(6) * (t45 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(7) * (t2 ^ 2 + t23 ^ 2 + t3 ^ 2) + m(5) * (t32 ^ 2 + t33 ^ 2 + t66 ^ 2) + m(4) * (t123 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(3) * (pkin(1) ^ 2 * t183 ^ 2 + t135 ^ 2 + t136 ^ 2) + (t202 - 0.2e1 * t211 + 0.2e1 * t212) * t185 + ((-0.2e1 * t135 * mrSges(3,3) + Ifges(3,5) * t185 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t189) * t183) * t189 + (0.2e1 * t136 * mrSges(3,3) + Ifges(3,6) * t185 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t189 + (Ifges(4,3) + Ifges(3,2)) * t193) * t183 + t205) * t193) * t183 + (t8 + t20 + t46 - t81) * t132 + Ifges(2,3) + 0.2e1 * t3 * t18 + 0.2e1 * t2 * t19 + 0.2e1 * t23 * t11 + t29 * t9 + t30 * t10 + 0.2e1 * t7 * t43 + 0.2e1 * t6 * t44 + 0.2e1 * t45 * t31 + t58 * t21 + t59 * t22 + 0.2e1 * t66 * t60 + 0.2e1 * t33 * t69 + 0.2e1 * t32 * t70 + t102 * t47 + t103 * t48 + 0.2e1 * t72 * t107 + 0.2e1 * t71 * t108 + 0.2e1 * t123 * t88 + t133 * t82; (t72 * mrSges(4,3) + pkin(9) * t107 - t8 / 0.2e1 - t20 / 0.2e1 - t46 / 0.2e1 + t81 / 0.2e1) * t192 + m(7) * (t13 * t2 + t14 * t3 + t23 * t96) + m(6) * (t149 * t45 + t53 * t6 + t54 * t7) + m(4) * (-pkin(2) * t123 + (-t71 * t188 + t72 * t192) * pkin(9)) + m(5) * (t114 * t32 + t115 * t33 + t177 * t66) - t204 * t208 / 0.2e1 + t22 * t224 + t21 * t225 + t103 * t226 + t102 * t227 + (t34 / 0.2e1 + t76 / 0.2e1 + t120 / 0.2e1 - t158 / 0.2e1) * t132 - t211 + t212 + t202 + t10 * t230 + t9 * t231 + t78 * t232 + t77 * t233 + t36 * t235 + t35 * t236 + (-t71 * mrSges(4,3) - t182 * t47 / 0.2e1 + t184 * t234 + t82 / 0.2e1 + (t60 - t108) * pkin(9)) * t188 + t14 * t18 + t13 * t19 + t23 * t39 + t53 * t44 + t54 * t43 + t3 * t67 + t2 * t68 + t45 * t85 - pkin(2) * t88 + t96 * t11 + t7 * t109 + t6 * t110 + t114 * t70 + t115 * t69 + t66 * t134 + t33 * t147 + t32 * t148 + t149 * t31 + t123 * t157 + t133 * t159 / 0.2e1; -0.2e1 * pkin(2) * t157 + 0.2e1 * t54 * t109 + 0.2e1 * t53 * t110 + 0.2e1 * t114 * t148 + 0.2e1 * t115 * t147 - t126 * t77 + t127 * t78 + 0.2e1 * t13 * t68 + 0.2e1 * t14 * t67 + 0.2e1 * t149 * t85 + t79 * t35 + t80 * t36 + 0.2e1 * t96 * t39 + Ifges(3,3) + (t180 + t181) * mrSges(4,3) * t237 + (-t34 - t76 - t120 + t158) * t192 + (-t121 * t182 + t122 * t184 + t134 * t237 + t159) * t188 + m(7) * (t13 ^ 2 + t14 ^ 2 + t96 ^ 2) + m(6) * (t149 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t114 ^ 2 + t115 ^ 2 + t176) + m(4) * (pkin(2) ^ 2 + t181 * t195 + t176); (t33 * mrSges(5,3) + qJ(4) * t69 + t47 / 0.2e1) * t184 - Ifges(4,3) * t208 + t201 * t132 + t103 * t221 + t22 * t222 + t21 * t223 + t10 * t228 + t9 * t229 + (-t2 * t93 + t3 * t92) * mrSges(7,3) + (t7 * t145 - t6 * t146) * mrSges(6,3) - t205 + m(5) * (-pkin(3) * t66 + (-t32 * t182 + t33 * t184) * qJ(4)) + m(7) * (t116 * t23 + t2 * t40 + t3 * t41) + m(6) * (t105 * t6 + t106 * t7 + t171 * t45) + t101 * t232 + t100 * t233 + t52 * t235 + t51 * t236 + (-t32 * mrSges(5,3) - qJ(4) * t70 + t234) * t182 + t40 * t19 + t41 * t18 + t23 * t49 - pkin(3) * t60 + t71 * mrSges(4,1) - t72 * mrSges(4,2) + t45 * t98 + t105 * t44 + t106 * t43 + t116 * t11 + t66 * t152 + t102 * t155 / 0.2e1 + t171 * t31; (-t13 * t93 + t14 * t92) * mrSges(7,3) + (t54 * t145 - t53 * t146) * mrSges(6,3) + m(5) * (-pkin(3) * t177 + (-t114 * t182 + t115 * t184) * qJ(4)) + (-pkin(9) * mrSges(4,2) - t201) * t192 + (-mrSges(4,1) + t152) * t177 + t78 * t222 + t77 * t223 + t101 * t224 + t100 * t225 + t36 * t228 + t35 * t229 + t52 * t230 + t204 + t51 * t231 + (-t114 * mrSges(5,3) - qJ(4) * t148 - t188 * t155 / 0.2e1 + t226) * t182 + (t115 * mrSges(5,3) + qJ(4) * t147 + t188 * t221 + t227) * t184 + m(7) * (t116 * t96 + t13 * t40 + t14 * t41) + m(6) * (t105 * t53 + t106 * t54 + t149 * t171) + t41 * t67 + t40 * t68 + t96 * t49 + t106 * t109 + t105 * t110 + t116 * t39 - pkin(3) * t134 + t149 * t98 + t171 * t85; -0.2e1 * pkin(3) * t152 + t145 * t100 + t146 * t101 + 0.2e1 * t116 * t49 + t184 * t155 + t182 * t156 + 0.2e1 * t171 * t98 + t92 * t51 + t93 * t52 + Ifges(4,3) + m(7) * (t116 ^ 2 + t40 ^ 2 + t41 ^ 2) + m(6) * (t105 ^ 2 + t106 ^ 2 + t171 ^ 2) + m(5) * (qJ(4) ^ 2 * t203 + pkin(3) ^ 2) + 0.2e1 * (-t40 * t93 + t41 * t92) * mrSges(7,3) + 0.2e1 * (-t105 * t146 + t106 * t145) * mrSges(6,3) + 0.2e1 * t203 * qJ(4) * mrSges(5,3); m(5) * t66 + m(6) * t45 + m(7) * t23 + t11 + t31 + t60; m(5) * t177 + m(6) * t149 + m(7) * t96 + t134 + t39 + t85; -m(5) * pkin(3) + m(6) * t171 + m(7) * t116 + t152 + t49 + t98; m(5) + m(6) + m(7); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(7) * (t186 * t3 + t190 * t2) + t186 * t18 + t190 * t19) * pkin(5) + t198 + t20; t53 * mrSges(6,1) - t54 * mrSges(6,2) + t218 * t192 + (m(7) * (t13 * t190 + t14 * t186) + t186 * t67 + t190 * t68) * pkin(5) + t200 + t206; t105 * mrSges(6,1) - t106 * mrSges(6,2) + (m(7) * (t186 * t41 + t190 * t40) + (t186 * t92 - t190 * t93) * mrSges(7,3)) * pkin(5) + t199 + t99; 0; m(7) * (t186 ^ 2 + t190 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t197 - t218; t198; t200 - t213; t199; 0; Ifges(7,3) + t197; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
