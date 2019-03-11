% Calculate joint inertia matrix for
% S6RRRRPR7
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

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR7_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR7_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR7_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:37
% EndTime: 2019-03-09 22:26:42
% DurationCPUTime: 2.12s
% Computational Cost: add. (5393->394), mult. (11881->578), div. (0->0), fcn. (13712->12), ass. (0->150)
t165 = sin(qJ(6));
t169 = cos(qJ(6));
t166 = sin(qJ(4));
t167 = sin(qJ(3));
t170 = cos(qJ(4));
t171 = cos(qJ(3));
t129 = -t166 * t167 + t170 * t171;
t151 = -pkin(3) * t171 - pkin(2);
t109 = -pkin(4) * t129 + t151;
t130 = t166 * t171 + t167 * t170;
t161 = sin(pkin(12));
t163 = cos(pkin(12));
t96 = -t163 * t129 + t130 * t161;
t97 = t129 * t161 + t130 * t163;
t50 = pkin(5) * t96 - pkin(11) * t97 + t109;
t221 = -pkin(10) - pkin(9);
t189 = t221 * t167;
t190 = t221 * t171;
t104 = t166 * t190 + t170 * t189;
t178 = -t130 * qJ(5) + t104;
t105 = t166 * t189 - t170 * t190;
t86 = qJ(5) * t129 + t105;
t53 = t161 * t178 + t163 * t86;
t23 = -t165 * t53 + t169 * t50;
t24 = t165 * t50 + t169 * t53;
t181 = -t165 * t23 + t169 * t24;
t164 = cos(pkin(6));
t162 = sin(pkin(6));
t168 = sin(qJ(2));
t197 = t162 * t168;
t121 = t164 * t171 - t167 * t197;
t122 = t164 * t167 + t171 * t197;
t87 = t121 * t170 - t122 * t166;
t88 = t121 * t166 + t122 * t170;
t58 = t161 * t88 - t163 * t87;
t59 = t161 * t87 + t163 * t88;
t143 = pkin(8) * t197;
t172 = cos(qJ(2));
t215 = pkin(1) * t172;
t110 = t143 + (-pkin(2) - t215) * t164;
t90 = -pkin(3) * t121 + t110;
t61 = -pkin(4) * t87 + t90;
t15 = pkin(5) * t58 - pkin(11) * t59 + t61;
t196 = t162 * t172;
t124 = t164 * t168 * pkin(1) + pkin(8) * t196;
t111 = pkin(9) * t164 + t124;
t112 = (-pkin(2) * t172 - pkin(9) * t168 - pkin(1)) * t162;
t77 = -t111 * t167 + t171 * t112;
t67 = -pkin(3) * t196 - pkin(10) * t122 + t77;
t78 = t171 * t111 + t167 * t112;
t72 = pkin(10) * t121 + t78;
t30 = -t166 * t72 + t170 * t67;
t18 = -pkin(4) * t196 - qJ(5) * t88 + t30;
t31 = t166 * t67 + t170 * t72;
t27 = qJ(5) * t87 + t31;
t9 = t161 * t18 + t163 * t27;
t6 = -pkin(11) * t196 + t9;
t2 = t15 * t169 - t165 * t6;
t3 = t15 * t165 + t169 * t6;
t184 = -t165 * t2 + t169 * t3;
t227 = -Ifges(4,5) * t122 - Ifges(4,6) * t121;
t226 = Ifges(5,5) * t130 + Ifges(6,5) * t97 + Ifges(5,6) * t129 - Ifges(6,6) * t96;
t51 = t161 * t86 - t163 * t178;
t225 = t51 ^ 2;
t224 = 0.2e1 * t51;
t133 = -mrSges(7,1) * t169 + mrSges(7,2) * t165;
t223 = 0.2e1 * t133;
t41 = -t165 * t59 - t169 * t196;
t222 = t41 / 0.2e1;
t135 = Ifges(7,5) * t165 + Ifges(7,6) * t169;
t220 = t135 / 0.2e1;
t208 = Ifges(7,4) * t169;
t138 = Ifges(7,1) * t165 + t208;
t219 = t138 / 0.2e1;
t218 = t165 / 0.2e1;
t217 = t169 / 0.2e1;
t214 = pkin(3) * t166;
t211 = Ifges(5,3) + Ifges(6,3);
t209 = Ifges(7,4) * t165;
t150 = pkin(3) * t170 + pkin(4);
t118 = t150 * t163 - t161 * t214;
t207 = t118 * mrSges(6,1);
t119 = t161 * t150 + t163 * t214;
t206 = t119 * mrSges(6,2);
t123 = t164 * t215 - t143;
t205 = t123 * mrSges(3,1);
t204 = t124 * mrSges(3,2);
t202 = t165 * t97;
t200 = t169 * t97;
t148 = pkin(4) * t161 + pkin(11);
t199 = t148 * t165;
t198 = t148 * t169;
t194 = Ifges(4,5) * t167 + Ifges(4,6) * t171;
t193 = t165 ^ 2 + t169 ^ 2;
t192 = t167 ^ 2 + t171 ^ 2;
t191 = Ifges(4,3) + t211;
t42 = -t165 * t196 + t169 * t59;
t12 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t58;
t188 = Ifges(3,5) * t197 + Ifges(3,6) * t196 + Ifges(3,3) * t164;
t28 = t58 * mrSges(6,1) + t59 * mrSges(6,2);
t68 = t96 * mrSges(6,1) + t97 * mrSges(6,2);
t187 = -Ifges(5,5) * t88 - Ifges(6,5) * t59 - Ifges(5,6) * t87 + Ifges(6,6) * t58;
t186 = t193 * t148;
t136 = Ifges(7,2) * t169 + t209;
t185 = t169 * t136 + t165 * t138 + t211;
t183 = t163 * mrSges(6,1) - t161 * mrSges(6,2);
t182 = mrSges(7,1) * t165 + mrSges(7,2) * t169;
t8 = -t161 * t27 + t163 * t18;
t180 = 0.2e1 * t193 * mrSges(7,3);
t34 = Ifges(7,5) * t200 - Ifges(7,6) * t202 + Ifges(7,3) * t96;
t179 = (mrSges(5,1) * t170 - mrSges(5,2) * t166) * pkin(3);
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t58;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t58;
t5 = pkin(5) * t196 - t8;
t177 = t30 * mrSges(5,1) + t8 * mrSges(6,1) - t31 * mrSges(5,2) - t9 * mrSges(6,2) + mrSges(7,3) * t184 + t13 * t217 + t5 * t133 + t136 * t222 + t14 * t218 + t42 * t219 + t58 * t220 - t187;
t35 = Ifges(7,6) * t96 + (-Ifges(7,2) * t165 + t208) * t97;
t36 = Ifges(7,5) * t96 + (Ifges(7,1) * t169 - t209) * t97;
t176 = -t105 * mrSges(5,2) - t53 * mrSges(6,2) + t35 * t217 + t36 * t218 + t104 * mrSges(5,1) - t136 * t202 / 0.2e1 + t200 * t219 + t96 * t220 + t226 + (-mrSges(6,1) + t133) * t51 + t181 * mrSges(7,3);
t149 = -pkin(4) * t163 - pkin(5);
t139 = Ifges(4,1) * t167 + Ifges(4,4) * t171;
t137 = Ifges(4,4) * t167 + Ifges(4,2) * t171;
t134 = -mrSges(4,1) * t171 + mrSges(4,2) * t167;
t115 = pkin(11) + t119;
t114 = -pkin(5) - t118;
t103 = -mrSges(4,1) * t196 - mrSges(4,3) * t122;
t102 = mrSges(4,2) * t196 + mrSges(4,3) * t121;
t100 = Ifges(5,1) * t130 + Ifges(5,4) * t129;
t99 = Ifges(5,4) * t130 + Ifges(5,2) * t129;
t98 = -mrSges(5,1) * t129 + mrSges(5,2) * t130;
t91 = -mrSges(4,1) * t121 + mrSges(4,2) * t122;
t82 = Ifges(4,1) * t122 + Ifges(4,4) * t121 - Ifges(4,5) * t196;
t81 = Ifges(4,4) * t122 + Ifges(4,2) * t121 - Ifges(4,6) * t196;
t74 = -mrSges(5,1) * t196 - mrSges(5,3) * t88;
t73 = mrSges(5,2) * t196 + mrSges(5,3) * t87;
t70 = Ifges(6,1) * t97 - Ifges(6,4) * t96;
t69 = Ifges(6,4) * t97 - Ifges(6,2) * t96;
t66 = mrSges(7,1) * t96 - mrSges(7,3) * t200;
t65 = -mrSges(7,2) * t96 - mrSges(7,3) * t202;
t62 = t182 * t97;
t60 = -mrSges(5,1) * t87 + mrSges(5,2) * t88;
t48 = Ifges(5,1) * t88 + Ifges(5,4) * t87 - Ifges(5,5) * t196;
t47 = Ifges(5,4) * t88 + Ifges(5,2) * t87 - Ifges(5,6) * t196;
t46 = -mrSges(6,1) * t196 - mrSges(6,3) * t59;
t45 = mrSges(6,2) * t196 - mrSges(6,3) * t58;
t26 = Ifges(6,1) * t59 - Ifges(6,4) * t58 - Ifges(6,5) * t196;
t25 = Ifges(6,4) * t59 - Ifges(6,2) * t58 - Ifges(6,6) * t196;
t20 = mrSges(7,1) * t58 - mrSges(7,3) * t42;
t19 = -mrSges(7,2) * t58 + mrSges(7,3) * t41;
t16 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t1 = [(t188 - 0.2e1 * t204 + 0.2e1 * t205) * t164 + (t12 - t25) * t58 + ((-0.2e1 * t123 * mrSges(3,3) + Ifges(3,5) * t164 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t168) * t162) * t168 + (0.2e1 * t124 * mrSges(3,3) + Ifges(3,6) * t164 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t168 + (Ifges(3,2) + t191) * t172) * t162 + t187 + t227) * t172) * t162 + m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + m(6) * (t61 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t90 ^ 2) + m(4) * (t110 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(3) * (pkin(1) ^ 2 * t162 ^ 2 + t123 ^ 2 + t124 ^ 2) + Ifges(2,3) + 0.2e1 * t5 * t16 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + t41 * t13 + t42 * t14 + 0.2e1 * t9 * t45 + 0.2e1 * t8 * t46 + t59 * t26 + 0.2e1 * t61 * t28 + 0.2e1 * t31 * t73 + 0.2e1 * t30 * t74 + t87 * t47 + t88 * t48 + 0.2e1 * t90 * t60 + 0.2e1 * t78 * t102 + 0.2e1 * t77 * t103 + 0.2e1 * t110 * t91 + t121 * t81 + t122 * t82; m(7) * (t2 * t23 + t24 * t3 + t5 * t51) + m(6) * (t109 * t61 - t51 * t8 + t53 * t9) + m(5) * (t104 * t30 + t105 * t31 + t151 * t90) + (-t9 * mrSges(6,3) + t12 / 0.2e1 - t25 / 0.2e1) * t96 + (-t77 * mrSges(4,3) - pkin(9) * t103 + t82 / 0.2e1) * t167 + (t129 * t31 - t130 * t30) * mrSges(5,3) + (t34 / 0.2e1 - t69 / 0.2e1) * t58 + m(4) * (-pkin(2) * t110 + (-t167 * t77 + t171 * t78) * pkin(9)) - (t194 + t226) * t196 / 0.2e1 + t188 - t204 + t205 + (t16 - t46) * t51 + (t78 * mrSges(4,3) + pkin(9) * t102 + t81 / 0.2e1) * t171 + (-t8 * mrSges(6,3) - t165 * t13 / 0.2e1 + t14 * t217 + t26 / 0.2e1) * t97 + t35 * t222 + t23 * t20 + t24 * t19 + t42 * t36 / 0.2e1 + t53 * t45 + t5 * t62 + t3 * t65 + t2 * t66 + t61 * t68 + t59 * t70 / 0.2e1 - pkin(2) * t91 + t90 * t98 + t87 * t99 / 0.2e1 + t88 * t100 / 0.2e1 + t104 * t74 + t105 * t73 + t109 * t28 + t129 * t47 / 0.2e1 + t130 * t48 / 0.2e1 + t110 * t134 + t121 * t137 / 0.2e1 + t122 * t139 / 0.2e1 + t151 * t60; -0.2e1 * pkin(2) * t134 + t130 * t100 + 0.2e1 * t109 * t68 + t129 * t99 + t171 * t137 + t167 * t139 + 0.2e1 * t151 * t98 + 0.2e1 * t23 * t66 + 0.2e1 * t24 * t65 + t62 * t224 + Ifges(3,3) + (-0.2e1 * mrSges(6,3) * t53 + t34 - t69) * t96 + (mrSges(6,3) * t224 - t165 * t35 + t169 * t36 + t70) * t97 + m(7) * (t23 ^ 2 + t24 ^ 2 + t225) + m(6) * (t109 ^ 2 + t53 ^ 2 + t225) + m(5) * (t104 ^ 2 + t105 ^ 2 + t151 ^ 2) + m(4) * (t192 * pkin(9) ^ 2 + pkin(2) ^ 2) + 0.2e1 * (-t104 * t130 + t105 * t129) * mrSges(5,3) + 0.2e1 * t192 * pkin(9) * mrSges(4,3); (-t165 * t20 + t169 * t19) * t115 + (t166 * t73 + t170 * t74 + m(5) * (t166 * t31 + t170 * t30)) * pkin(3) + m(7) * (t114 * t5 + t184 * t115) + t177 - t191 * t196 + t77 * mrSges(4,1) - t78 * mrSges(4,2) + t114 * t16 + t118 * t46 + t119 * t45 + m(6) * (t118 * t8 + t119 * t9) - t227; (-mrSges(4,1) * t167 - mrSges(4,2) * t171) * pkin(9) + (-t118 * t97 - t119 * t96) * mrSges(6,3) + (m(5) * (t104 * t170 + t105 * t166) + (t129 * t166 - t130 * t170) * mrSges(5,3)) * pkin(3) + m(7) * (t114 * t51 + t181 * t115) + t176 + m(6) * (-t118 * t51 + t119 * t53) + t114 * t62 + (-t165 * t66 + t169 * t65) * t115 + t194; 0.2e1 * t207 - 0.2e1 * t206 + t114 * t223 + Ifges(4,3) + 0.2e1 * t179 + t115 * t180 + m(7) * (t193 * t115 ^ 2 + t114 ^ 2) + m(6) * (t118 ^ 2 + t119 ^ 2) + m(5) * (t166 ^ 2 + t170 ^ 2) * pkin(3) ^ 2 + t185; m(7) * (t184 * t148 + t149 * t5) + (t161 * t45 + t163 * t46 + m(6) * (t161 * t9 + t163 * t8)) * pkin(4) + t177 - t211 * t196 - t20 * t199 + t19 * t198 + t149 * t16; m(7) * (t181 * t148 + t149 * t51) + (m(6) * (t161 * t53 - t163 * t51) + (-t161 * t96 - t163 * t97) * mrSges(6,3)) * pkin(4) + t176 - t66 * t199 + t65 * t198 + t149 * t62; m(7) * (t114 * t149 + t115 * t186) - t206 + t207 + (t114 + t149) * t133 + t179 + (m(6) * (t118 * t163 + t119 * t161) + t183) * pkin(4) + (t193 * t115 + t186) * mrSges(7,3) + t185; t149 * t223 + t148 * t180 + m(7) * (t193 * t148 ^ 2 + t149 ^ 2) + t185 + (0.2e1 * t183 + m(6) * (t161 ^ 2 + t163 ^ 2) * pkin(4)) * pkin(4); t165 * t19 + t169 * t20 + m(7) * (t165 * t3 + t169 * t2) + m(6) * t61 + t28; t165 * t65 + t169 * t66 + m(7) * (t165 * t24 + t169 * t23) + m(6) * t109 + t68; 0; 0; m(7) * t193 + m(6); mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t23 - mrSges(7,2) * t24 + t34; -t182 * t115 + t135; -t182 * t148 + t135; -t133; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
