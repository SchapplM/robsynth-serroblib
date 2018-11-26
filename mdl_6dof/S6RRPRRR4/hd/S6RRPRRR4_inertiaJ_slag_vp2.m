% Calculate joint inertia matrix for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2018-11-23 17:23
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:22:39
% EndTime: 2018-11-23 17:22:41
% DurationCPUTime: 1.93s
% Computational Cost: add. (4147->353), mult. (9752->517), div. (0->0), fcn. (10926->12), ass. (0->141)
t152 = sin(qJ(6));
t156 = cos(qJ(6));
t178 = t152 ^ 2 + t156 ^ 2;
t220 = mrSges(7,3) * t178;
t219 = Ifges(3,3) + Ifges(4,3);
t153 = sin(qJ(5));
t154 = sin(qJ(4));
t157 = cos(qJ(5));
t158 = cos(qJ(4));
t115 = t153 * t158 + t154 * t157;
t150 = cos(pkin(12));
t134 = -pkin(2) * t150 - pkin(3);
t117 = -pkin(4) * t158 + t134;
t113 = t153 * t154 - t157 * t158;
t203 = pkin(5) * t113;
t65 = -pkin(11) * t115 + t117 + t203;
t148 = sin(pkin(12));
t133 = pkin(2) * t148 + pkin(9);
t200 = pkin(10) + t133;
t105 = t200 * t158;
t174 = t200 * t154;
t69 = t157 * t105 - t153 * t174;
t40 = -t152 * t69 + t156 * t65;
t41 = t152 * t65 + t156 * t69;
t169 = -t152 * t40 + t156 * t41;
t151 = cos(pkin(6));
t159 = cos(qJ(2));
t204 = pkin(1) * t151;
t130 = t159 * t204;
t149 = sin(pkin(6));
t155 = sin(qJ(2));
t182 = t149 * t155;
t82 = pkin(2) * t151 + t130 + (-pkin(8) - qJ(3)) * t182;
t181 = t149 * t159;
t104 = pkin(8) * t181 + t155 * t204;
t91 = qJ(3) * t181 + t104;
t52 = -t148 * t91 + t150 * t82;
t50 = -pkin(3) * t151 - t52;
t98 = (t148 * t159 + t150 * t155) * t149;
t77 = t151 * t158 - t154 * t98;
t32 = -pkin(4) * t77 + t50;
t78 = t151 * t154 + t158 * t98;
t46 = t153 * t78 - t157 * t77;
t47 = t153 * t77 + t157 * t78;
t15 = pkin(5) * t46 - pkin(11) * t47 + t32;
t53 = t148 * t82 + t150 * t91;
t51 = pkin(9) * t151 + t53;
t116 = (-pkin(2) * t159 - pkin(1)) * t149;
t97 = t148 * t182 - t150 * t181;
t56 = pkin(3) * t97 - pkin(9) * t98 + t116;
t26 = -t154 * t51 + t158 * t56;
t18 = pkin(4) * t97 - pkin(10) * t78 + t26;
t27 = t154 * t56 + t158 * t51;
t24 = pkin(10) * t77 + t27;
t9 = t153 * t18 + t157 * t24;
t6 = pkin(11) * t97 + t9;
t2 = t15 * t156 - t152 * t6;
t3 = t15 * t152 + t156 * t6;
t171 = -t152 * t2 + t156 * t3;
t218 = Ifges(4,5) * t98 - Ifges(4,6) * t97;
t186 = t115 * t152;
t72 = -mrSges(7,2) * t113 - mrSges(7,3) * t186;
t185 = t115 * t156;
t73 = mrSges(7,1) * t113 - mrSges(7,3) * t185;
t217 = -t152 * t73 + t156 * t72;
t33 = -t152 * t47 + t156 * t97;
t19 = -mrSges(7,2) * t46 + mrSges(7,3) * t33;
t34 = t152 * t97 + t156 * t47;
t20 = mrSges(7,1) * t46 - mrSges(7,3) * t34;
t216 = -t152 * t20 + t156 * t19;
t215 = Ifges(5,5) * t78 + Ifges(5,6) * t77 + Ifges(5,3) * t97;
t67 = t105 * t153 + t157 * t174;
t214 = t67 ^ 2;
t213 = 0.2e1 * t67;
t212 = t113 ^ 2;
t211 = 0.2e1 * t116;
t210 = t33 / 0.2e1;
t120 = Ifges(7,5) * t152 + Ifges(7,6) * t156;
t208 = t120 / 0.2e1;
t196 = Ifges(7,4) * t156;
t123 = Ifges(7,1) * t152 + t196;
t207 = t123 / 0.2e1;
t206 = t152 / 0.2e1;
t205 = t156 / 0.2e1;
t16 = -mrSges(7,1) * t33 + mrSges(7,2) * t34;
t36 = mrSges(6,1) * t97 - mrSges(6,3) * t47;
t199 = t16 - t36;
t197 = Ifges(7,4) * t152;
t103 = -pkin(8) * t182 + t130;
t195 = t103 * mrSges(3,1);
t194 = t104 * mrSges(3,2);
t193 = t113 * t67;
t136 = pkin(4) * t153 + pkin(11);
t184 = t136 * t152;
t183 = t136 * t156;
t180 = Ifges(6,5) * t115 - Ifges(6,6) * t113;
t179 = Ifges(5,5) * t154 + Ifges(5,6) * t158;
t177 = t154 ^ 2 + t158 ^ 2;
t12 = Ifges(7,5) * t34 + Ifges(7,6) * t33 + Ifges(7,3) * t46;
t176 = Ifges(6,5) * t47 - Ifges(6,6) * t46 + Ifges(6,3) * t97;
t121 = Ifges(7,2) * t156 + t197;
t175 = t156 * t121 + t152 * t123 + Ifges(6,3);
t173 = t178 * pkin(11);
t79 = t113 * mrSges(6,1) + t115 * mrSges(6,2);
t172 = t178 * t136;
t119 = -t158 * mrSges(5,1) + t154 * mrSges(5,2);
t170 = mrSges(7,1) * t152 + mrSges(7,2) * t156;
t8 = -t153 * t24 + t157 * t18;
t168 = 0.2e1 * t220;
t61 = Ifges(7,5) * t185 - Ifges(7,6) * t186 + Ifges(7,3) * t113;
t167 = Ifges(3,5) * t182 + Ifges(3,6) * t181 + t151 * t219 + t218;
t166 = (mrSges(6,1) * t157 - mrSges(6,2) * t153) * pkin(4);
t118 = -mrSges(7,1) * t156 + mrSges(7,2) * t152;
t165 = t113 * t118 + t115 * t220 - t79;
t13 = Ifges(7,4) * t34 + Ifges(7,2) * t33 + Ifges(7,6) * t46;
t14 = Ifges(7,1) * t34 + Ifges(7,4) * t33 + Ifges(7,5) * t46;
t5 = -pkin(5) * t97 - t8;
t164 = t8 * mrSges(6,1) - t9 * mrSges(6,2) + mrSges(7,3) * t171 + t5 * t118 + t121 * t210 + t13 * t205 + t14 * t206 + t34 * t207 + t46 * t208 + t176;
t62 = Ifges(7,6) * t113 + (-Ifges(7,2) * t152 + t196) * t115;
t63 = Ifges(7,5) * t113 + (Ifges(7,1) * t156 - t197) * t115;
t163 = -t69 * mrSges(6,2) + t62 * t205 + t63 * t206 + t180 - t121 * t186 / 0.2e1 + t185 * t207 + t113 * t208 + (-mrSges(6,1) + t118) * t67 + t169 * mrSges(7,3);
t137 = -pkin(4) * t157 - pkin(5);
t124 = Ifges(5,1) * t154 + Ifges(5,4) * t158;
t122 = Ifges(5,4) * t154 + Ifges(5,2) * t158;
t110 = t115 ^ 2;
t92 = t98 * mrSges(4,2);
t86 = mrSges(4,1) * t151 - mrSges(4,3) * t98;
t85 = -mrSges(4,2) * t151 - mrSges(4,3) * t97;
t81 = Ifges(6,1) * t115 - Ifges(6,4) * t113;
t80 = Ifges(6,4) * t115 - Ifges(6,2) * t113;
t70 = t170 * t115;
t60 = mrSges(5,1) * t97 - mrSges(5,3) * t78;
t59 = -mrSges(5,2) * t97 + mrSges(5,3) * t77;
t49 = -mrSges(5,1) * t77 + mrSges(5,2) * t78;
t39 = Ifges(5,1) * t78 + Ifges(5,4) * t77 + Ifges(5,5) * t97;
t38 = Ifges(5,4) * t78 + Ifges(5,2) * t77 + Ifges(5,6) * t97;
t35 = -mrSges(6,2) * t97 - mrSges(6,3) * t46;
t25 = mrSges(6,1) * t46 + mrSges(6,2) * t47;
t22 = Ifges(6,1) * t47 - Ifges(6,4) * t46 + Ifges(6,5) * t97;
t21 = Ifges(6,4) * t47 - Ifges(6,2) * t46 + Ifges(6,6) * t97;
t1 = [Ifges(4,1) * t98 ^ 2 + (t167 - 0.2e1 * t194 + 0.2e1 * t195 + t218) * t151 + m(3) * (t103 ^ 2 + t104 ^ 2) + m(4) * (t116 ^ 2 + t52 ^ 2 + t53 ^ 2) + m(5) * (t26 ^ 2 + t27 ^ 2 + t50 ^ 2) + m(6) * (t32 ^ 2 + t8 ^ 2 + t9 ^ 2) + m(7) * (t2 ^ 2 + t3 ^ 2 + t5 ^ 2) + (mrSges(4,1) * t211 - 0.2e1 * Ifges(4,4) * t98 + Ifges(4,2) * t97 + t176 + t215) * t97 + t92 * t211 + (t12 - t21) * t46 + ((Ifges(3,5) * t155 + Ifges(3,6) * t159) * t151 + 0.2e1 * (-t103 * t155 + t104 * t159) * mrSges(3,3) + (-0.2e1 * pkin(1) * (-mrSges(3,1) * t159 + mrSges(3,2) * t155) + t159 * (Ifges(3,4) * t155 + Ifges(3,2) * t159) + t155 * (Ifges(3,1) * t155 + Ifges(3,4) * t159) + m(3) * pkin(1) ^ 2) * t149) * t149 + Ifges(2,3) + 0.2e1 * t5 * t16 + 0.2e1 * t3 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t32 * t25 + t33 * t13 + t34 * t14 + 0.2e1 * t9 * t35 + 0.2e1 * t8 * t36 + t47 * t22 + 0.2e1 * t50 * t49 + 0.2e1 * t27 * t59 + 0.2e1 * t26 * t60 + t77 * t38 + t78 * t39 + 0.2e1 * t53 * t85 + 0.2e1 * t52 * t86; (t27 * mrSges(5,3) + t133 * t59 + t38 / 0.2e1) * t158 + (-t26 * mrSges(5,3) - t133 * t60 + t39 / 0.2e1) * t154 + (-t9 * mrSges(6,3) + t12 / 0.2e1 - t21 / 0.2e1) * t113 + m(6) * (t117 * t32 - t67 * t8 + t69 * t9) + m(7) * (t2 * t40 + t3 * t41 + t5 * t67) + (t148 * t85 + t150 * t86 + m(4) * (t148 * t53 + t150 * t52)) * pkin(2) + (t180 + t179) * t97 / 0.2e1 + (t61 / 0.2e1 - t80 / 0.2e1) * t46 + m(5) * (t134 * t50 + (-t26 * t154 + t27 * t158) * t133) + t167 + t62 * t210 - t194 + t195 + (-t8 * mrSges(6,3) - t152 * t13 / 0.2e1 + t14 * t205 + t22 / 0.2e1) * t115 + t40 * t20 + t41 * t19 + t52 * mrSges(4,1) - t53 * mrSges(4,2) + t34 * t63 / 0.2e1 + t69 * t35 + t5 * t70 + t3 * t72 + t2 * t73 + t32 * t79 + t47 * t81 / 0.2e1 + t117 * t25 + t50 * t119 + t77 * t122 / 0.2e1 + t78 * t124 / 0.2e1 + t134 * t49 + t199 * t67; 0.2e1 * t117 * t79 + 0.2e1 * t134 * t119 + t158 * t122 + t154 * t124 + 0.2e1 * t40 * t73 + 0.2e1 * t41 * t72 + t70 * t213 + (-0.2e1 * mrSges(6,3) * t69 + t61 - t80) * t113 + (mrSges(6,3) * t213 - t152 * t62 + t156 * t63 + t81) * t115 + m(7) * (t40 ^ 2 + t41 ^ 2 + t214) + m(6) * (t117 ^ 2 + t69 ^ 2 + t214) + m(5) * (t133 ^ 2 * t177 + t134 ^ 2) + m(4) * (t148 ^ 2 + t150 ^ 2) * pkin(2) ^ 2 + 0.2e1 * (mrSges(4,1) * t150 - mrSges(4,2) * t148) * pkin(2) + 0.2e1 * t177 * t133 * mrSges(5,3) + t219; t97 * mrSges(4,1) + t154 * t59 + t158 * t60 + t92 + t199 * t113 + (t35 + t216) * t115 + m(7) * (t113 * t5 + t115 * t171) + m(6) * (-t113 * t8 + t115 * t9) + m(5) * (t154 * t27 + t158 * t26) + m(4) * t116; t113 * t70 + t217 * t115 + m(7) * (t115 * t169 + t193) + m(6) * (t115 * t69 + t193); m(4) + m(5) * t177 + m(6) * (t110 + t212) + m(7) * (t110 * t178 + t212); m(7) * (t136 * t171 + t137 * t5) + (t153 * t35 + t157 * t36 + m(6) * (t153 * t9 + t157 * t8)) * pkin(4) + t164 - t20 * t184 + t19 * t183 + t26 * mrSges(5,1) - t27 * mrSges(5,2) + t137 * t16 + t215; (-mrSges(5,1) * t154 - mrSges(5,2) * t158) * t133 + m(7) * (t136 * t169 + t137 * t67) + (m(6) * (t153 * t69 - t157 * t67) + (-t113 * t153 - t115 * t157) * mrSges(6,3)) * pkin(4) + t163 - t73 * t184 + t72 * t183 + t137 * t70 + t179; m(7) * (t137 * t113 + t115 * t172) + m(6) * (-t113 * t157 + t115 * t153) * pkin(4) + t165 - t119; 0.2e1 * t137 * t118 + Ifges(5,3) + 0.2e1 * t166 + t136 * t168 + m(7) * (t136 ^ 2 * t178 + t137 ^ 2) + m(6) * (t153 ^ 2 + t157 ^ 2) * pkin(4) ^ 2 + t175; t164 + (-m(7) * t5 - t16) * pkin(5) + (m(7) * t171 + t216) * pkin(11); t163 + (-m(7) * t67 - t70) * pkin(5) + (m(7) * t169 + t217) * pkin(11); m(7) * (t115 * t173 - t203) + t165; m(7) * (-pkin(5) * t137 + pkin(11) * t172) + (t137 - pkin(5)) * t118 + t166 + (t172 + t173) * mrSges(7,3) + t175; -0.2e1 * pkin(5) * t118 + m(7) * (pkin(11) ^ 2 * t178 + pkin(5) ^ 2) + pkin(11) * t168 + t175; mrSges(7,1) * t2 - mrSges(7,2) * t3 + t12; mrSges(7,1) * t40 - mrSges(7,2) * t41 + t61; -t70; -t136 * t170 + t120; -pkin(11) * t170 + t120; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
