% Calculate joint inertia matrix for
% S6RRRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2018-11-23 18:33
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRP9_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP9_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP9_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP9_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:32:51
% EndTime: 2018-11-23 18:32:53
% DurationCPUTime: 1.96s
% Computational Cost: add. (3884->458), mult. (8767->630), div. (0->0), fcn. (9569->10), ass. (0->171)
t233 = 2 * pkin(9);
t175 = cos(pkin(6));
t178 = sin(qJ(3));
t182 = cos(qJ(3));
t174 = sin(pkin(6));
t179 = sin(qJ(2));
t213 = t174 * t179;
t120 = -t175 * t182 + t178 * t213;
t176 = sin(qJ(5));
t180 = cos(qJ(5));
t121 = t175 * t178 + t182 * t213;
t177 = sin(qJ(4));
t181 = cos(qJ(4));
t183 = cos(qJ(2));
t212 = t174 * t183;
t84 = -t121 * t177 - t181 * t212;
t85 = t121 * t181 - t177 * t212;
t42 = -t176 * t85 + t180 * t84;
t43 = t176 * t84 + t180 * t85;
t10 = Ifges(7,5) * t43 + Ifges(7,6) * t42 + Ifges(7,3) * t120;
t11 = Ifges(6,5) * t43 + Ifges(6,6) * t42 + Ifges(6,3) * t120;
t232 = t10 + t11;
t231 = (t180 * mrSges(6,1) + (-mrSges(6,2) - mrSges(7,2)) * t176) * pkin(4);
t230 = 2 * mrSges(7,1);
t229 = m(7) * pkin(5);
t33 = Ifges(5,1) * t85 + Ifges(5,4) * t84 + Ifges(5,5) * t120;
t228 = t33 / 0.2e1;
t227 = -pkin(11) - pkin(10);
t218 = Ifges(5,4) * t181;
t111 = -Ifges(5,6) * t182 + (-Ifges(5,2) * t177 + t218) * t178;
t226 = t111 / 0.2e1;
t219 = Ifges(5,4) * t177;
t112 = -Ifges(5,5) * t182 + (Ifges(5,1) * t181 - t219) * t178;
t225 = t112 / 0.2e1;
t146 = Ifges(5,1) * t177 + t218;
t224 = t146 / 0.2e1;
t152 = pkin(8) * t213;
t222 = pkin(1) * t183;
t102 = t152 + (-pkin(2) - t222) * t175;
t48 = pkin(3) * t120 - pkin(10) * t121 + t102;
t123 = t175 * t179 * pkin(1) + pkin(8) * t212;
t103 = pkin(9) * t175 + t123;
t104 = (-pkin(2) * t183 - pkin(9) * t179 - pkin(1)) * t174;
t54 = t182 * t103 + t178 * t104;
t50 = -pkin(10) * t212 + t54;
t22 = t177 * t48 + t181 * t50;
t18 = pkin(11) * t84 + t22;
t21 = -t177 * t50 + t181 * t48;
t9 = pkin(4) * t120 - pkin(11) * t85 + t21;
t6 = t176 * t9 + t180 * t18;
t223 = m(7) * t176;
t221 = pkin(9) * t182;
t169 = t178 * pkin(9);
t220 = Ifges(6,3) + Ifges(7,3);
t140 = -pkin(3) * t182 - pkin(10) * t178 - pkin(2);
t131 = t181 * t140;
t210 = t178 * t181;
t71 = -pkin(11) * t210 + t131 + (-pkin(9) * t177 - pkin(4)) * t182;
t100 = t177 * t140 + t181 * t221;
t211 = t177 * t178;
t86 = -pkin(11) * t211 + t100;
t36 = t176 * t71 + t180 * t86;
t122 = t175 * t222 - t152;
t217 = t122 * mrSges(3,1);
t216 = t123 * mrSges(3,2);
t132 = -t176 * t177 + t180 * t181;
t214 = t132 * t176;
t133 = t176 * t181 + t177 * t180;
t113 = t133 * t178;
t114 = t132 * t178;
t209 = Ifges(7,5) * t114 - Ifges(7,6) * t113;
t208 = Ifges(6,5) * t114 - Ifges(6,6) * t113;
t207 = -Ifges(4,5) * t121 + Ifges(4,6) * t120;
t78 = Ifges(7,5) * t133 + Ifges(7,6) * t132;
t79 = Ifges(6,5) * t133 + Ifges(6,6) * t132;
t148 = t227 * t177;
t149 = t227 * t181;
t91 = t176 * t148 - t180 * t149;
t143 = Ifges(5,5) * t177 + Ifges(5,6) * t181;
t206 = Ifges(4,5) * t178 + Ifges(4,6) * t182;
t139 = pkin(4) * t211 + t169;
t205 = t177 ^ 2 + t181 ^ 2;
t204 = -Ifges(5,3) - t220;
t31 = Ifges(5,5) * t85 + Ifges(5,6) * t84 + Ifges(5,3) * t120;
t12 = Ifges(7,4) * t43 + Ifges(7,2) * t42 + Ifges(7,6) * t120;
t13 = Ifges(6,4) * t43 + Ifges(6,2) * t42 + Ifges(6,6) * t120;
t203 = t12 / 0.2e1 + t13 / 0.2e1;
t14 = Ifges(7,1) * t43 + Ifges(7,4) * t42 + Ifges(7,5) * t120;
t15 = Ifges(6,1) * t43 + Ifges(6,4) * t42 + Ifges(6,5) * t120;
t202 = t14 / 0.2e1 + t15 / 0.2e1;
t57 = Ifges(7,4) * t114 - Ifges(7,2) * t113 - Ifges(7,6) * t182;
t58 = Ifges(6,4) * t114 - Ifges(6,2) * t113 - Ifges(6,6) * t182;
t201 = t57 / 0.2e1 + t58 / 0.2e1;
t59 = Ifges(7,1) * t114 - Ifges(7,4) * t113 - Ifges(7,5) * t182;
t60 = Ifges(6,1) * t114 - Ifges(6,4) * t113 - Ifges(6,5) * t182;
t200 = t59 / 0.2e1 + t60 / 0.2e1;
t80 = Ifges(7,4) * t133 + Ifges(7,2) * t132;
t81 = Ifges(6,4) * t133 + Ifges(6,2) * t132;
t199 = t80 / 0.2e1 + t81 / 0.2e1;
t82 = Ifges(7,1) * t133 + Ifges(7,4) * t132;
t83 = Ifges(6,1) * t133 + Ifges(6,4) * t132;
t198 = t82 / 0.2e1 + t83 / 0.2e1;
t5 = -t176 * t18 + t180 * t9;
t2 = pkin(5) * t120 - qJ(6) * t43 + t5;
t28 = mrSges(7,1) * t120 - mrSges(7,3) * t43;
t197 = m(7) * t2 + t28;
t196 = Ifges(3,5) * t213 + Ifges(3,6) * t212 + Ifges(3,3) * t175;
t35 = -t176 * t86 + t180 * t71;
t24 = -pkin(5) * t182 - qJ(6) * t114 + t35;
t94 = -mrSges(7,1) * t182 - mrSges(7,3) * t114;
t195 = m(7) * t24 + t94;
t161 = -pkin(4) * t181 - pkin(3);
t19 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t66 = t113 * mrSges(7,1) + t114 * mrSges(7,2);
t76 = -t132 * mrSges(7,1) + t133 * mrSges(7,2);
t53 = -t178 * t103 + t104 * t182;
t90 = t180 * t148 + t149 * t176;
t193 = t78 / 0.2e1 + t79 / 0.2e1 + t143 / 0.2e1;
t192 = Ifges(5,5) * t210 - Ifges(5,6) * t211;
t49 = pkin(3) * t212 - t53;
t64 = -qJ(6) * t133 + t90;
t191 = m(7) * t64 - t133 * mrSges(7,3);
t190 = mrSges(5,1) * t177 + mrSges(5,2) * t181;
t30 = -pkin(4) * t84 + t49;
t25 = -qJ(6) * t113 + t36;
t189 = t35 * mrSges(6,1) + t24 * mrSges(7,1) - t36 * mrSges(6,2) - t25 * mrSges(7,2) + t208 + t209;
t65 = qJ(6) * t132 + t91;
t188 = t90 * mrSges(6,1) + t64 * mrSges(7,1) - t91 * mrSges(6,2) - t65 * mrSges(7,2) + t78 + t79;
t3 = qJ(6) * t42 + t6;
t187 = t5 * mrSges(6,1) + t2 * mrSges(7,1) - t6 * mrSges(6,2) - t3 * mrSges(7,2) + t232;
t186 = pkin(4) ^ 2;
t185 = pkin(9) ^ 2;
t173 = t182 ^ 2;
t171 = t178 ^ 2;
t168 = t171 * t185;
t167 = t176 ^ 2 * t186;
t160 = pkin(4) * t180 + pkin(5);
t147 = Ifges(4,1) * t178 + Ifges(4,4) * t182;
t145 = Ifges(4,4) * t178 + Ifges(4,2) * t182;
t144 = Ifges(5,2) * t181 + t219;
t142 = -mrSges(4,1) * t182 + mrSges(4,2) * t178;
t141 = -mrSges(5,1) * t181 + mrSges(5,2) * t177;
t138 = -mrSges(5,1) * t182 - mrSges(5,3) * t210;
t137 = mrSges(5,2) * t182 - mrSges(5,3) * t211;
t124 = t190 * t178;
t110 = -Ifges(5,3) * t182 + t192;
t101 = -pkin(5) * t132 + t161;
t99 = -t177 * t221 + t131;
t95 = -mrSges(6,1) * t182 - mrSges(6,3) * t114;
t93 = mrSges(6,2) * t182 - mrSges(6,3) * t113;
t92 = mrSges(7,2) * t182 - mrSges(7,3) * t113;
t89 = -mrSges(4,1) * t212 - mrSges(4,3) * t121;
t88 = mrSges(4,2) * t212 - mrSges(4,3) * t120;
t77 = -mrSges(6,1) * t132 + mrSges(6,2) * t133;
t72 = pkin(5) * t113 + t139;
t68 = mrSges(4,1) * t120 + mrSges(4,2) * t121;
t67 = mrSges(6,1) * t113 + mrSges(6,2) * t114;
t62 = Ifges(4,1) * t121 - Ifges(4,4) * t120 - Ifges(4,5) * t212;
t61 = Ifges(4,4) * t121 - Ifges(4,2) * t120 - Ifges(4,6) * t212;
t56 = -Ifges(6,3) * t182 + t208;
t55 = -Ifges(7,3) * t182 + t209;
t52 = mrSges(5,1) * t120 - mrSges(5,3) * t85;
t51 = -mrSges(5,2) * t120 + mrSges(5,3) * t84;
t44 = -mrSges(5,1) * t84 + mrSges(5,2) * t85;
t32 = Ifges(5,4) * t85 + Ifges(5,2) * t84 + Ifges(5,6) * t120;
t29 = mrSges(6,1) * t120 - mrSges(6,3) * t43;
t27 = -mrSges(6,2) * t120 + mrSges(6,3) * t42;
t26 = -mrSges(7,2) * t120 + mrSges(7,3) * t42;
t20 = -mrSges(6,1) * t42 + mrSges(6,2) * t43;
t16 = -pkin(5) * t42 + t30;
t1 = [m(7) * (t16 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t30 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t21 ^ 2 + t22 ^ 2 + t49 ^ 2) + m(4) * (t102 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(3) * (pkin(1) ^ 2 * t174 ^ 2 + t122 ^ 2 + t123 ^ 2) + 0.2e1 * t3 * t26 + 0.2e1 * t6 * t27 + 0.2e1 * t2 * t28 + 0.2e1 * t5 * t29 + 0.2e1 * t30 * t20 + 0.2e1 * t16 * t19 + (t14 + t15) * t43 + (t12 + t13) * t42 + (t196 - 0.2e1 * t216 + 0.2e1 * t217) * t175 + ((-0.2e1 * t122 * mrSges(3,3) + Ifges(3,5) * t175 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t179) * t174) * t179 + (0.2e1 * t123 * mrSges(3,3) + Ifges(3,6) * t175 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t179 + (Ifges(4,3) + Ifges(3,2)) * t183) * t174 + t207) * t183) * t174 + (t31 - t61 + t232) * t120 + 0.2e1 * t49 * t44 + 0.2e1 * t22 * t51 + 0.2e1 * t21 * t52 + t84 * t32 + t85 * t33 + 0.2e1 * t54 * t88 + 0.2e1 * t53 * t89 + Ifges(2,3) + 0.2e1 * t102 * t68 + t121 * t62; t35 * t29 + t36 * t27 + t25 * t26 + t24 * t28 + (t54 * mrSges(4,3) + pkin(9) * t88 - t31 / 0.2e1 - t10 / 0.2e1 - t11 / 0.2e1 + t61 / 0.2e1) * t182 + m(7) * (t16 * t72 + t2 * t24 + t25 * t3) + m(6) * (t139 * t30 + t35 * t5 + t36 * t6) + m(4) * (-pkin(2) * t102 + (-t53 * t178 + t54 * t182) * pkin(9)) + t196 + (-t53 * mrSges(4,3) - t177 * t32 / 0.2e1 + t181 * t228 + t62 / 0.2e1 + (t44 - t89) * pkin(9)) * t178 + m(5) * (t100 * t22 + t49 * t169 + t21 * t99) - t206 * t212 / 0.2e1 - t203 * t113 + t200 * t43 + t201 * t42 + t202 * t114 - t216 + t217 + (t55 / 0.2e1 + t56 / 0.2e1 + t110 / 0.2e1 - t145 / 0.2e1) * t120 + t16 * t66 + t30 * t67 - pkin(2) * t68 + t72 * t19 + t3 * t92 + t6 * t93 + t2 * t94 + t5 * t95 + t99 * t52 + t100 * t51 + t49 * t124 + t22 * t137 + t21 * t138 + t139 * t20 + t102 * t142 + t121 * t147 / 0.2e1 + t85 * t225 + t84 * t226; -0.2e1 * pkin(2) * t142 + 0.2e1 * t100 * t137 + 0.2e1 * t99 * t138 + 0.2e1 * t139 * t67 + 0.2e1 * t24 * t94 + 0.2e1 * t25 * t92 + 0.2e1 * t35 * t95 + 0.2e1 * t36 * t93 + 0.2e1 * t72 * t66 + Ifges(3,3) + (t59 + t60) * t114 - (t57 + t58) * t113 + (t171 + t173) * mrSges(4,3) * t233 + (-t55 - t56 - t110 + t145) * t182 + (-t111 * t177 + t112 * t181 + t124 * t233 + t147) * t178 + m(7) * (t24 ^ 2 + t25 ^ 2 + t72 ^ 2) + m(6) * (t139 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t100 ^ 2 + t99 ^ 2 + t168) + m(4) * (pkin(2) ^ 2 + t173 * t185 + t168); m(5) * (-pkin(3) * t49 + (-t21 * t177 + t22 * t181) * pkin(10)) + (t22 * mrSges(5,3) + pkin(10) * t51 + t32 / 0.2e1) * t181 + (-t21 * mrSges(5,3) - pkin(10) * t52 + t228) * t177 - Ifges(4,3) * t212 + (t6 * mrSges(6,3) + t3 * mrSges(7,3) + t203) * t132 + t198 * t43 + t199 * t42 + (-t5 * mrSges(6,3) - t2 * mrSges(7,3) + t202) * t133 + t193 * t120 + m(7) * (t101 * t16 + t2 * t64 + t3 * t65) + m(6) * (t161 * t30 + t5 * t90 + t6 * t91) - t207 - pkin(3) * t44 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + t64 * t28 + t65 * t26 + t16 * t76 + t30 * t77 + t90 * t29 + t91 * t27 + t101 * t19 + t49 * t141 + t84 * t144 / 0.2e1 + t161 * t20 + t85 * t224; m(5) * (-pkin(3) * t169 + (t100 * t181 - t177 * t99) * pkin(10)) + m(7) * (t101 * t72 + t24 * t64 + t25 * t65) + m(6) * (t139 * t161 + t35 * t90 + t36 * t91) + (-mrSges(4,1) + t141) * t169 + t198 * t114 - t199 * t113 + (-t35 * mrSges(6,3) - t24 * mrSges(7,3) + t200) * t133 + (t36 * mrSges(6,3) + t25 * mrSges(7,3) + t201) * t132 + (-pkin(9) * mrSges(4,2) - t193) * t182 + t206 + (-t99 * mrSges(5,3) - pkin(10) * t138 - t178 * t144 / 0.2e1 + t225) * t177 + t72 * t76 + t65 * t92 + (t100 * mrSges(5,3) + pkin(10) * t137 + t178 * t224 + t226) * t181 + t91 * t93 + t64 * t94 + t90 * t95 + t101 * t66 - pkin(3) * t124 + t139 * t77 + t161 * t67; -0.2e1 * pkin(3) * t141 + 0.2e1 * t101 * t76 + t181 * t144 + t177 * t146 + 0.2e1 * t161 * t77 + Ifges(4,3) + 0.2e1 * t205 * pkin(10) * mrSges(5,3) + m(7) * (t101 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t161 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(5) * (t205 * pkin(10) ^ 2 + pkin(3) ^ 2) + (-0.2e1 * mrSges(6,3) * t90 - 0.2e1 * mrSges(7,3) * t64 + t82 + t83) * t133 + (0.2e1 * mrSges(6,3) * t91 + 0.2e1 * mrSges(7,3) * t65 + t80 + t81) * t132; t21 * mrSges(5,1) - t22 * mrSges(5,2) + t197 * t160 + (t180 * t29 + (t26 + t27) * t176 + t3 * t223 + m(6) * (t176 * t6 + t180 * t5)) * pkin(4) + t187 + t31; t99 * mrSges(5,1) - t100 * mrSges(5,2) + t195 * t160 + (t180 * t95 + (t92 + t93) * t176 + t25 * t223 + m(6) * (t176 * t36 + t180 * t35)) * pkin(4) + t204 * t182 + t189 + t192; t191 * t160 - t190 * pkin(10) + (mrSges(7,3) * t214 + (-t133 * t180 + t214) * mrSges(6,3) + t65 * t223 + m(6) * (t176 * t91 + t180 * t90)) * pkin(4) + t188 + t143; t160 * t230 + m(7) * (t160 ^ 2 + t167) + m(6) * (t180 ^ 2 * t186 + t167) + 0.2e1 * t231 - t204; pkin(5) * t197 + t187; t195 * pkin(5) - t220 * t182 + t189; pkin(5) * t191 + t188; t160 * t229 + (pkin(5) + t160) * mrSges(7,1) + t231 + t220; (t230 + t229) * pkin(5) + t220; m(7) * t16 + t19; m(7) * t72 + t66; m(7) * t101 + t76; 0; 0; m(7);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
