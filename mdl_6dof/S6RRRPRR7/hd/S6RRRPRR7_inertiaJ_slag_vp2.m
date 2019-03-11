% Calculate joint inertia matrix for
% S6RRRPRR7
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:52
% EndTime: 2019-03-09 18:32:57
% DurationCPUTime: 2.08s
% Computational Cost: add. (5180->379), mult. (11437->548), div. (0->0), fcn. (13253->12), ass. (0->145)
t162 = sin(qJ(6));
t166 = cos(qJ(6));
t158 = sin(pkin(12));
t160 = cos(pkin(12));
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t127 = -t158 * t164 + t160 * t167;
t148 = -pkin(3) * t167 - pkin(2);
t108 = -pkin(4) * t127 + t148;
t128 = t158 * t167 + t160 * t164;
t163 = sin(qJ(5));
t212 = cos(qJ(5));
t96 = -t212 * t127 + t128 * t163;
t97 = t163 * t127 + t212 * t128;
t47 = pkin(5) * t96 - pkin(11) * t97 + t108;
t207 = -qJ(4) - pkin(9);
t180 = t207 * t164;
t181 = t207 * t167;
t101 = t158 * t181 + t160 * t180;
t174 = -t128 * pkin(10) + t101;
t102 = t158 * t180 - t160 * t181;
t83 = pkin(10) * t127 + t102;
t52 = t163 * t174 + t212 * t83;
t22 = -t162 * t52 + t166 * t47;
t23 = t162 * t47 + t166 * t52;
t176 = -t162 * t22 + t166 * t23;
t161 = cos(pkin(6));
t159 = sin(pkin(6));
t165 = sin(qJ(2));
t190 = t159 * t165;
t120 = t161 * t167 - t164 * t190;
t121 = t161 * t164 + t167 * t190;
t87 = t120 * t160 - t121 * t158;
t88 = t120 * t158 + t121 * t160;
t57 = t163 * t88 - t212 * t87;
t58 = t163 * t87 + t212 * t88;
t143 = pkin(8) * t190;
t168 = cos(qJ(2));
t211 = pkin(1) * t168;
t109 = t143 + (-pkin(2) - t211) * t161;
t90 = -pkin(3) * t120 + t109;
t60 = -pkin(4) * t87 + t90;
t15 = pkin(5) * t57 - pkin(11) * t58 + t60;
t189 = t159 * t168;
t123 = t161 * t165 * pkin(1) + pkin(8) * t189;
t110 = pkin(9) * t161 + t123;
t111 = (-pkin(2) * t168 - pkin(9) * t165 - pkin(1)) * t159;
t76 = -t110 * t164 + t167 * t111;
t69 = -pkin(3) * t189 - qJ(4) * t121 + t76;
t77 = t167 * t110 + t164 * t111;
t71 = qJ(4) * t120 + t77;
t29 = -t158 * t71 + t160 * t69;
t19 = -pkin(4) * t189 - pkin(10) * t88 + t29;
t30 = t158 * t69 + t160 * t71;
t27 = pkin(10) * t87 + t30;
t9 = t163 * t19 + t212 * t27;
t6 = -pkin(11) * t189 + t9;
t2 = t15 * t166 - t162 * t6;
t3 = t15 * t162 + t166 * t6;
t178 = -t162 * t2 + t166 * t3;
t225 = Ifges(4,3) + Ifges(5,3);
t224 = Ifges(4,5) * t164 + Ifges(5,5) * t128 + Ifges(4,6) * t167 + Ifges(5,6) * t127;
t195 = t162 * t97;
t64 = -mrSges(7,2) * t96 - mrSges(7,3) * t195;
t191 = t166 * t97;
t65 = mrSges(7,1) * t96 - mrSges(7,3) * t191;
t223 = -t162 * t65 + t166 * t64;
t41 = -t162 * t58 - t166 * t189;
t20 = -mrSges(7,2) * t57 + mrSges(7,3) * t41;
t42 = -t162 * t189 + t166 * t58;
t21 = mrSges(7,1) * t57 - mrSges(7,3) * t42;
t222 = -t162 * t21 + t166 * t20;
t221 = -Ifges(4,5) * t121 - Ifges(5,5) * t88 - Ifges(4,6) * t120 - Ifges(5,6) * t87;
t50 = t163 * t83 - t212 * t174;
t220 = t50 ^ 2;
t219 = 0.2e1 * t50;
t218 = t41 / 0.2e1;
t135 = Ifges(7,5) * t162 + Ifges(7,6) * t166;
t217 = t135 / 0.2e1;
t203 = Ifges(7,4) * t166;
t138 = Ifges(7,1) * t162 + t203;
t216 = t138 / 0.2e1;
t215 = t162 / 0.2e1;
t214 = t166 / 0.2e1;
t210 = pkin(3) * t158;
t206 = -Ifges(6,5) * t58 + Ifges(6,6) * t57;
t205 = Ifges(6,5) * t97 - Ifges(6,6) * t96;
t204 = Ifges(7,4) * t162;
t147 = pkin(3) * t160 + pkin(4);
t117 = t212 * t147 - t163 * t210;
t202 = t117 * mrSges(6,1);
t118 = t163 * t147 + t212 * t210;
t201 = t118 * mrSges(6,2);
t122 = t161 * t211 - t143;
t200 = t122 * mrSges(3,1);
t199 = t123 * mrSges(3,2);
t186 = t162 ^ 2 + t166 ^ 2;
t185 = t164 ^ 2 + t167 ^ 2;
t184 = Ifges(6,3) + t225;
t12 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t57;
t136 = Ifges(7,2) * t166 + t204;
t183 = t166 * t136 + t162 * t138 + Ifges(6,3);
t182 = Ifges(3,5) * t190 + Ifges(3,6) * t189 + Ifges(3,3) * t161;
t59 = -t87 * mrSges(5,1) + t88 * mrSges(5,2);
t28 = t57 * mrSges(6,1) + t58 * mrSges(6,2);
t66 = t96 * mrSges(6,1) + t97 * mrSges(6,2);
t98 = -t127 * mrSges(5,1) + t128 * mrSges(5,2);
t114 = pkin(11) + t118;
t179 = t186 * t114;
t177 = mrSges(7,1) * t162 + mrSges(7,2) * t166;
t175 = 0.2e1 * t186 * mrSges(7,3);
t33 = Ifges(7,5) * t191 - Ifges(7,6) * t195 + Ifges(7,3) * t96;
t8 = -t163 * t27 + t212 * t19;
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t57;
t133 = -mrSges(7,1) * t166 + mrSges(7,2) * t162;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t57;
t5 = pkin(5) * t189 - t8;
t173 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + mrSges(7,3) * t178 + t13 * t214 + t5 * t133 + t136 * t218 + t14 * t215 + t42 * t216 + t57 * t217 - t206;
t34 = Ifges(7,6) * t96 + (-Ifges(7,2) * t162 + t203) * t97;
t35 = Ifges(7,5) * t96 + (Ifges(7,1) * t166 - t204) * t97;
t172 = -t52 * mrSges(6,2) + t34 * t214 + t35 * t215 + t205 - t136 * t195 / 0.2e1 + t191 * t216 + t96 * t217 + (-mrSges(6,1) + t133) * t50 + t176 * mrSges(7,3);
t139 = Ifges(4,1) * t164 + Ifges(4,4) * t167;
t137 = Ifges(4,4) * t164 + Ifges(4,2) * t167;
t134 = -mrSges(4,1) * t167 + mrSges(4,2) * t164;
t113 = -pkin(5) - t117;
t104 = -mrSges(4,1) * t189 - mrSges(4,3) * t121;
t103 = mrSges(4,2) * t189 + mrSges(4,3) * t120;
t100 = Ifges(5,1) * t128 + Ifges(5,4) * t127;
t99 = Ifges(5,4) * t128 + Ifges(5,2) * t127;
t91 = -mrSges(4,1) * t120 + mrSges(4,2) * t121;
t81 = Ifges(4,1) * t121 + Ifges(4,4) * t120 - Ifges(4,5) * t189;
t80 = Ifges(4,4) * t121 + Ifges(4,2) * t120 - Ifges(4,6) * t189;
t73 = -mrSges(5,1) * t189 - mrSges(5,3) * t88;
t72 = mrSges(5,2) * t189 + mrSges(5,3) * t87;
t68 = Ifges(6,1) * t97 - Ifges(6,4) * t96;
t67 = Ifges(6,4) * t97 - Ifges(6,2) * t96;
t61 = t177 * t97;
t49 = Ifges(5,1) * t88 + Ifges(5,4) * t87 - Ifges(5,5) * t189;
t48 = Ifges(5,4) * t88 + Ifges(5,2) * t87 - Ifges(5,6) * t189;
t45 = -mrSges(6,1) * t189 - mrSges(6,3) * t58;
t44 = mrSges(6,2) * t189 - mrSges(6,3) * t57;
t26 = Ifges(6,1) * t58 - Ifges(6,4) * t57 - Ifges(6,5) * t189;
t25 = Ifges(6,4) * t58 - Ifges(6,2) * t57 - Ifges(6,6) * t189;
t16 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t1 = [m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + m(6) * (t60 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t29 ^ 2 + t30 ^ 2 + t90 ^ 2) + m(4) * (t109 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(3) * (pkin(1) ^ 2 * t159 ^ 2 + t122 ^ 2 + t123 ^ 2) + (t182 - 0.2e1 * t199 + 0.2e1 * t200) * t161 + ((-0.2e1 * t122 * mrSges(3,3) + Ifges(3,5) * t161 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t165) * t159) * t165 + (0.2e1 * t123 * mrSges(3,3) + Ifges(3,6) * t161 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t165 + (Ifges(3,2) + t184) * t168) * t159 + t206 + t221) * t168) * t159 + Ifges(2,3) + 0.2e1 * t5 * t16 + 0.2e1 * t3 * t20 + 0.2e1 * t2 * t21 + t41 * t13 + t42 * t14 + 0.2e1 * t9 * t44 + 0.2e1 * t8 * t45 + t58 * t26 + 0.2e1 * t60 * t28 + 0.2e1 * t30 * t72 + 0.2e1 * t29 * t73 + t87 * t48 + t88 * t49 + 0.2e1 * t90 * t59 + 0.2e1 * t77 * t103 + 0.2e1 * t76 * t104 + 0.2e1 * t109 * t91 + t120 * t80 + t121 * t81 + (t12 - t25) * t57; (t16 - t45) * t50 + (t77 * mrSges(4,3) + pkin(9) * t103 + t80 / 0.2e1) * t167 + (-t76 * mrSges(4,3) - pkin(9) * t104 + t81 / 0.2e1) * t164 + m(7) * (t2 * t22 + t23 * t3 + t5 * t50) + m(6) * (t108 * t60 - t50 * t8 + t52 * t9) + m(5) * (t101 * t29 + t102 * t30 + t148 * t90) + t182 + (-t9 * mrSges(6,3) + t12 / 0.2e1 - t25 / 0.2e1) * t96 + (t33 / 0.2e1 - t67 / 0.2e1) * t57 + (t30 * t127 - t29 * t128) * mrSges(5,3) - (t205 + t224) * t189 / 0.2e1 + m(4) * (-pkin(2) * t109 + (-t76 * t164 + t77 * t167) * pkin(9)) + (-t8 * mrSges(6,3) - t162 * t13 / 0.2e1 + t14 * t214 + t26 / 0.2e1) * t97 + t34 * t218 + t22 * t21 + t23 * t20 + t42 * t35 / 0.2e1 + t52 * t44 + t5 * t61 + t3 * t64 + t2 * t65 + t60 * t66 + t58 * t68 / 0.2e1 - pkin(2) * t91 + t90 * t98 + t87 * t99 / 0.2e1 + t88 * t100 / 0.2e1 + t101 * t73 + t102 * t72 + t108 * t28 + t127 * t48 / 0.2e1 + t128 * t49 / 0.2e1 + t109 * t134 + t120 * t137 / 0.2e1 + t121 * t139 / 0.2e1 + t148 * t59 - t199 + t200; -0.2e1 * pkin(2) * t134 + t128 * t100 + 0.2e1 * t108 * t66 + t127 * t99 + t167 * t137 + t164 * t139 + 0.2e1 * t148 * t98 + 0.2e1 * t22 * t65 + 0.2e1 * t23 * t64 + t61 * t219 + Ifges(3,3) + (-0.2e1 * mrSges(6,3) * t52 + t33 - t67) * t96 + (mrSges(6,3) * t219 - t162 * t34 + t166 * t35 + t68) * t97 + m(7) * (t22 ^ 2 + t23 ^ 2 + t220) + m(6) * (t108 ^ 2 + t52 ^ 2 + t220) + m(5) * (t101 ^ 2 + t102 ^ 2 + t148 ^ 2) + m(4) * (t185 * pkin(9) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (-t101 * t128 + t102 * t127) * mrSges(5,3) + 0.2e1 * t185 * pkin(9) * mrSges(4,3); -t184 * t189 + m(6) * (t117 * t8 + t118 * t9) + m(7) * (t113 * t5 + t178 * t114) + t173 + (t158 * t72 + t160 * t73 + m(5) * (t158 * t30 + t160 * t29)) * pkin(3) + t222 * t114 + t29 * mrSges(5,1) - t30 * mrSges(5,2) + t76 * mrSges(4,1) - t77 * mrSges(4,2) + t113 * t16 + t117 * t45 + t118 * t44 - t221; t223 * t114 + m(6) * (-t117 * t50 + t118 * t52) + (-mrSges(4,1) * t164 - mrSges(4,2) * t167) * pkin(9) + (-t117 * t97 - t118 * t96) * mrSges(6,3) + m(7) * (t113 * t50 + t176 * t114) + (m(5) * (t101 * t160 + t102 * t158) + (t127 * t158 - t128 * t160) * mrSges(5,3)) * pkin(3) + t172 + t101 * mrSges(5,1) - t102 * mrSges(5,2) + t113 * t61 + t224; 0.2e1 * t202 - 0.2e1 * t201 + 0.2e1 * t113 * t133 + t114 * t175 + m(7) * (t186 * t114 ^ 2 + t113 ^ 2) + m(6) * (t117 ^ 2 + t118 ^ 2) + t183 + (0.2e1 * mrSges(5,1) * t160 - 0.2e1 * mrSges(5,2) * t158 + m(5) * (t158 ^ 2 + t160 ^ 2) * pkin(3)) * pkin(3) + t225; t162 * t20 + t166 * t21 + m(7) * (t162 * t3 + t166 * t2) + m(6) * t60 + m(5) * t90 + t28 + t59; t162 * t64 + t166 * t65 + m(7) * (t162 * t23 + t166 * t22) + m(6) * t108 + m(5) * t148 + t98 + t66; 0; m(7) * t186 + m(5) + m(6); -Ifges(6,3) * t189 + t173 + (-m(7) * t5 - t16) * pkin(5) + (m(7) * t178 + t222) * pkin(11); t172 + (-m(7) * t50 - t61) * pkin(5) + (m(7) * t176 + t223) * pkin(11); m(7) * (-pkin(5) * t113 + pkin(11) * t179) - t201 + t202 + (t113 - pkin(5)) * t133 + (t186 * pkin(11) + t179) * mrSges(7,3) + t183; 0; -0.2e1 * pkin(5) * t133 + m(7) * (t186 * pkin(11) ^ 2 + pkin(5) ^ 2) + pkin(11) * t175 + t183; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t22 - mrSges(7,2) * t23 + t33; -t177 * t114 + t135; -t133; -t177 * pkin(11) + t135; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
