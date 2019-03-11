% Calculate joint inertia matrix for
% S6RRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR4_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR4_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR4_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR4_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR4_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR4_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:06:15
% EndTime: 2019-03-09 22:06:19
% DurationCPUTime: 1.48s
% Computational Cost: add. (3746->337), mult. (7140->474), div. (0->0), fcn. (8088->10), ass. (0->132)
t217 = Ifges(5,3) + Ifges(6,3);
t156 = sin(qJ(4));
t160 = cos(qJ(4));
t157 = sin(qJ(3));
t158 = sin(qJ(2));
t161 = cos(qJ(3));
t162 = cos(qJ(2));
t129 = t157 * t158 - t161 * t162;
t130 = t157 * t162 + t158 * t161;
t145 = -pkin(2) * t162 - pkin(1);
t83 = pkin(3) * t129 - pkin(9) * t130 + t145;
t207 = -pkin(8) - pkin(7);
t138 = t207 * t162;
t178 = t207 * t158;
t99 = -t161 * t138 + t157 * t178;
t52 = -t156 * t99 + t160 * t83;
t53 = t156 * t83 + t160 * t99;
t175 = -t156 * t52 + t160 * t53;
t153 = sin(pkin(11));
t154 = cos(pkin(11));
t122 = -t153 * t156 + t154 * t160;
t123 = t153 * t160 + t154 * t156;
t155 = sin(qJ(6));
t159 = cos(qJ(6));
t84 = t122 * t159 - t123 * t155;
t85 = t122 * t155 + t123 * t159;
t197 = Ifges(7,5) * t85 + Ifges(7,6) * t84;
t216 = Ifges(5,5) * t156 + Ifges(6,5) * t123 + Ifges(5,6) * t160 + Ifges(6,6) * t122 + t197;
t146 = t160 * qJ(5);
t29 = pkin(4) * t129 - t130 * t146 + t52;
t188 = t130 * t156;
t42 = -qJ(5) * t188 + t53;
t10 = -t153 * t42 + t154 * t29;
t73 = t122 * t130;
t5 = pkin(5) * t129 - pkin(10) * t73 + t10;
t11 = t153 * t29 + t154 * t42;
t72 = t123 * t130;
t6 = -pkin(10) * t72 + t11;
t3 = -t155 * t6 + t159 * t5;
t4 = t155 * t5 + t159 * t6;
t215 = mrSges(7,1) * t3 - t4 * mrSges(7,2);
t202 = pkin(10) * t123;
t142 = pkin(2) * t157 + pkin(9);
t115 = (-qJ(5) - t142) * t156;
t116 = t142 * t160 + t146;
t76 = t115 * t154 - t116 * t153;
t59 = t76 - t202;
t114 = t122 * pkin(10);
t77 = t115 * t153 + t116 * t154;
t60 = t114 + t77;
t21 = -t155 * t60 + t159 * t59;
t22 = t155 * t59 + t159 * t60;
t214 = mrSges(7,1) * t21 - t22 * mrSges(7,2);
t133 = (-qJ(5) - pkin(9)) * t156;
t135 = pkin(9) * t160 + t146;
t93 = t133 * t154 - t135 * t153;
t67 = t93 - t202;
t94 = t133 * t153 + t135 * t154;
t68 = t114 + t94;
t33 = -t155 * t68 + t159 * t67;
t34 = t155 * t67 + t159 * t68;
t213 = mrSges(7,1) * t33 - t34 * mrSges(7,2);
t97 = -t138 * t157 - t161 * t178;
t212 = t97 ^ 2;
t46 = -t84 * mrSges(7,1) + mrSges(7,2) * t85;
t211 = 0.2e1 * t46;
t88 = -t122 * mrSges(6,1) + mrSges(6,2) * t123;
t210 = 0.2e1 * t88;
t209 = 0.2e1 * t97;
t208 = 0.2e1 * t145;
t204 = pkin(2) * t161;
t203 = pkin(4) * t153;
t199 = t84 * mrSges(7,3);
t198 = t85 * mrSges(7,3);
t196 = Ifges(5,4) * t156;
t195 = Ifges(5,4) * t160;
t140 = pkin(4) * t154 + pkin(5);
t107 = t140 * t159 - t155 * t203;
t194 = t107 * mrSges(7,1);
t108 = t140 * t155 + t159 * t203;
t193 = t108 * mrSges(7,2);
t192 = t122 * mrSges(6,3);
t191 = t123 * mrSges(6,3);
t187 = t130 * t160;
t184 = t156 ^ 2 + t160 ^ 2;
t183 = t158 ^ 2 + t162 ^ 2;
t182 = 0.2e1 * mrSges(6,3);
t181 = 0.2e1 * mrSges(7,3);
t40 = -t155 * t73 - t159 * t72;
t41 = -t155 * t72 + t159 * t73;
t180 = Ifges(7,5) * t41 + Ifges(7,6) * t40 + Ifges(7,3) * t129;
t179 = t154 * t191;
t144 = -pkin(4) * t160 - pkin(3);
t44 = t72 * mrSges(6,1) + mrSges(6,2) * t73;
t14 = -t40 * mrSges(7,1) + mrSges(7,2) * t41;
t177 = t184 * t142;
t66 = pkin(4) * t188 + t97;
t176 = mrSges(5,1) * t156 + mrSges(5,2) * t160;
t86 = -mrSges(5,2) * t129 - mrSges(5,3) * t188;
t87 = mrSges(5,1) * t129 - mrSges(5,3) * t187;
t174 = -t156 * t87 + t160 * t86;
t173 = 0.2e1 * t184 * mrSges(5,3);
t102 = -pkin(5) * t122 + t144;
t172 = (mrSges(4,1) * t161 - mrSges(4,2) * t157) * pkin(2);
t136 = Ifges(5,2) * t160 + t196;
t137 = Ifges(5,1) * t156 + t195;
t47 = Ifges(7,4) * t85 + Ifges(7,2) * t84;
t48 = Ifges(7,1) * t85 + Ifges(7,4) * t84;
t89 = Ifges(6,4) * t123 + Ifges(6,2) * t122;
t90 = Ifges(6,1) * t123 + Ifges(6,4) * t122;
t171 = t122 * t89 + t123 * t90 + t136 * t160 + t137 * t156 + t47 * t84 + t48 * t85 + Ifges(4,3);
t170 = t46 + t88;
t169 = Ifges(5,5) * t187 + Ifges(6,5) * t73 - Ifges(6,6) * t72 + t129 * t217 + t180;
t168 = -t107 * t198 + t108 * t199 + t192 * t203 + t216;
t12 = Ifges(7,4) * t41 + Ifges(7,2) * t40 + Ifges(7,6) * t129;
t13 = Ifges(7,1) * t41 + Ifges(7,4) * t40 + Ifges(7,5) * t129;
t134 = -mrSges(5,1) * t160 + mrSges(5,2) * t156;
t30 = Ifges(6,4) * t73 - Ifges(6,2) * t72 + Ifges(6,6) * t129;
t31 = Ifges(6,1) * t73 - Ifges(6,4) * t72 + Ifges(6,5) * t129;
t43 = pkin(5) * t72 + t66;
t61 = Ifges(5,6) * t129 + (-Ifges(5,2) * t156 + t195) * t130;
t62 = Ifges(5,5) * t129 + (Ifges(5,1) * t160 - t196) * t130;
t167 = t160 * t61 / 0.2e1 + t156 * t62 / 0.2e1 + Ifges(4,5) * t130 + t122 * t30 / 0.2e1 + t123 * t31 / 0.2e1 - t72 * t89 / 0.2e1 + t73 * t90 / 0.2e1 - t99 * mrSges(4,2) + t84 * t12 / 0.2e1 + t85 * t13 / 0.2e1 + t66 * t88 + t43 * t46 + t40 * t47 / 0.2e1 + t41 * t48 / 0.2e1 + t11 * t192 + t4 * t199 - t3 * t198 - t10 * t191 + t137 * t187 / 0.2e1 - t136 * t188 / 0.2e1 + (t134 - mrSges(4,1)) * t97 + t175 * mrSges(5,3) + (-Ifges(4,6) + t216 / 0.2e1) * t129;
t143 = -pkin(3) - t204;
t132 = t144 - t204;
t100 = t102 - t204;
t82 = t176 * t130;
t55 = mrSges(6,1) * t129 - mrSges(6,3) * t73;
t54 = -mrSges(6,2) * t129 - mrSges(6,3) * t72;
t28 = mrSges(7,1) * t129 - mrSges(7,3) * t41;
t27 = -mrSges(7,2) * t129 + mrSges(7,3) * t40;
t1 = [-0.2e1 * pkin(1) * (-t162 * mrSges(3,1) + t158 * mrSges(3,2)) + t162 * (Ifges(3,4) * t158 + Ifges(3,2) * t162) + t158 * (Ifges(3,1) * t158 + Ifges(3,4) * t162) + 0.2e1 * t53 * t86 + 0.2e1 * t52 * t87 + 0.2e1 * t66 * t44 - t72 * t30 + t73 * t31 + 0.2e1 * t11 * t54 + 0.2e1 * t10 * t55 + t40 * t12 + t41 * t13 + 0.2e1 * t43 * t14 + 0.2e1 * t4 * t27 + 0.2e1 * t3 * t28 + t82 * t209 + m(6) * (t10 ^ 2 + t11 ^ 2 + t66 ^ 2) + m(7) * (t3 ^ 2 + t4 ^ 2 + t43 ^ 2) + (mrSges(4,1) * t208 - 0.2e1 * t99 * mrSges(4,3) + Ifges(4,2) * t129 + (-Ifges(5,6) * t156 - (2 * Ifges(4,4))) * t130 + t169) * t129 + (mrSges(4,2) * t208 + mrSges(4,3) * t209 + Ifges(4,1) * t130 - t156 * t61 + t160 * t62) * t130 + m(4) * (t145 ^ 2 + t99 ^ 2 + t212) + m(5) * (t52 ^ 2 + t53 ^ 2 + t212) + m(3) * (pkin(7) ^ 2 * t183 + pkin(1) ^ 2) + 0.2e1 * t183 * pkin(7) * mrSges(3,3) + Ifges(2,3); Ifges(3,5) * t158 + Ifges(3,6) * t162 + t143 * t82 + t132 * t44 + t100 * t14 + t76 * t55 + t77 * t54 + t22 * t27 + t21 * t28 + t174 * t142 + m(6) * (t10 * t76 + t11 * t77 + t132 * t66) + m(7) * (t100 * t43 + t21 * t3 + t22 * t4) + (-mrSges(3,1) * t158 - mrSges(3,2) * t162) * pkin(7) + m(5) * (t142 * t175 + t143 * t97) + t167 + (m(4) * (t157 * t99 - t161 * t97) + (-t129 * t157 - t130 * t161) * mrSges(4,3)) * pkin(2); 0.2e1 * t143 * t134 + t132 * t210 + t100 * t211 + m(4) * (t157 ^ 2 + t161 ^ 2) * pkin(2) ^ 2 + 0.2e1 * t172 + (t77 * t122 - t76 * t123) * t182 + (-t21 * t85 + t22 * t84) * t181 + m(5) * (t142 ^ 2 * t184 + t143 ^ 2) + t171 + m(7) * (t100 ^ 2 + t21 ^ 2 + t22 ^ 2) + m(6) * (t132 ^ 2 + t76 ^ 2 + t77 ^ 2) + t142 * t173 + Ifges(3,3); t144 * t44 + t93 * t55 + t94 * t54 + t102 * t14 - pkin(3) * t82 + t34 * t27 + t33 * t28 + t174 * pkin(9) + m(6) * (t10 * t93 + t11 * t94 + t144 * t66) + m(7) * (t102 * t43 + t3 * t33 + t34 * t4) + m(5) * (-pkin(3) * t97 + pkin(9) * t175) + t167; (pkin(9) * t184 + t177) * mrSges(5,3) + ((-t76 - t93) * t123 + (t77 + t94) * t122) * mrSges(6,3) + ((-t21 - t33) * t85 + (t22 + t34) * t84) * mrSges(7,3) + m(5) * (-pkin(3) * t143 + pkin(9) * t177) + t171 + m(6) * (t132 * t144 + t76 * t93 + t77 * t94) + m(7) * (t100 * t102 + t21 * t33 + t22 * t34) + (t144 + t132) * t88 + (t100 + t102) * t46 + (t143 - pkin(3)) * t134 + t172; -0.2e1 * pkin(3) * t134 + t102 * t211 + t144 * t210 + (-t33 * t85 + t34 * t84) * t181 + (t94 * t122 - t93 * t123) * t182 + pkin(9) * t173 + m(7) * (t102 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t144 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(5) * (pkin(9) ^ 2 * t184 + pkin(3) ^ 2) + t171; t107 * t28 + t108 * t27 + t52 * mrSges(5,1) - t53 * mrSges(5,2) - t11 * mrSges(6,2) + t10 * mrSges(6,1) + (m(6) * (t10 * t154 + t11 * t153) + t153 * t54 + t154 * t55) * pkin(4) + t169 + m(7) * (t107 * t3 + t108 * t4) - Ifges(5,6) * t188 + t215; t76 * mrSges(6,1) - t77 * mrSges(6,2) + (m(6) * (t153 * t77 + t154 * t76) - t179) * pkin(4) + t168 - t176 * t142 + m(7) * (t107 * t21 + t108 * t22) + t214; t93 * mrSges(6,1) - t94 * mrSges(6,2) + m(7) * (t107 * t33 + t108 * t34) + (m(6) * (t153 * t94 + t154 * t93) - t179) * pkin(4) + t168 - t176 * pkin(9) + t213; 0.2e1 * t194 - 0.2e1 * t193 + Ifges(7,3) + m(7) * (t107 ^ 2 + t108 ^ 2) + (0.2e1 * mrSges(6,1) * t154 - 0.2e1 * mrSges(6,2) * t153 + m(6) * (t153 ^ 2 + t154 ^ 2) * pkin(4)) * pkin(4) + t217; m(6) * t66 + m(7) * t43 + t14 + t44; m(6) * t132 + m(7) * t100 + t170; m(6) * t144 + m(7) * t102 + t170; 0; m(6) + m(7); t180 + t215; t197 + t214; t197 + t213; Ifges(7,3) - t193 + t194; 0; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
