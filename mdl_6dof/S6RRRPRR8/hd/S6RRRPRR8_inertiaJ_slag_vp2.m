% Calculate joint inertia matrix for
% S6RRRPRR8
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR8_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_inertiaJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR8_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR8_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:11
% EndTime: 2019-03-09 18:44:16
% DurationCPUTime: 2.00s
% Computational Cost: add. (4934->424), mult. (10910->626), div. (0->0), fcn. (12496->12), ass. (0->151)
t162 = sin(pkin(6));
t172 = cos(qJ(2));
t189 = t162 * t172;
t161 = sin(pkin(12));
t163 = cos(pkin(12));
t167 = sin(qJ(3));
t171 = cos(qJ(3));
t127 = t161 * t167 - t163 * t171;
t128 = t161 * t171 + t163 * t167;
t212 = Ifges(4,5) * t167 + Ifges(5,5) * t128 + Ifges(4,6) * t171 - Ifges(5,6) * t127;
t197 = -qJ(4) - pkin(9);
t135 = t197 * t171;
t182 = t197 * t167;
t98 = -t135 * t161 - t163 * t182;
t211 = t98 ^ 2;
t210 = 0.2e1 * t98;
t166 = sin(qJ(5));
t170 = cos(qJ(5));
t164 = cos(pkin(6));
t168 = sin(qJ(2));
t190 = t162 * t168;
t113 = t164 * t171 - t167 * t190;
t114 = t164 * t167 + t171 * t190;
t76 = t113 * t161 + t114 * t163;
t56 = -t166 * t76 - t170 * t189;
t57 = -t166 * t189 + t170 * t76;
t75 = -t163 * t113 + t114 * t161;
t21 = Ifges(6,1) * t57 + Ifges(6,4) * t56 + Ifges(6,5) * t75;
t209 = t21 / 0.2e1;
t195 = Ifges(6,4) * t170;
t59 = Ifges(6,6) * t127 + (-Ifges(6,2) * t166 + t195) * t128;
t208 = t59 / 0.2e1;
t196 = Ifges(6,4) * t166;
t60 = Ifges(6,5) * t127 + (Ifges(6,1) * t170 - t196) * t128;
t207 = t60 / 0.2e1;
t165 = sin(qJ(6));
t169 = cos(qJ(6));
t130 = -t165 * t166 + t169 * t170;
t131 = t165 * t170 + t166 * t169;
t96 = Ifges(7,4) * t131 + Ifges(7,2) * t130;
t206 = t96 / 0.2e1;
t97 = Ifges(7,1) * t131 + Ifges(7,4) * t130;
t205 = t97 / 0.2e1;
t204 = t130 / 0.2e1;
t203 = t131 / 0.2e1;
t139 = Ifges(6,1) * t166 + t195;
t202 = t139 / 0.2e1;
t200 = pkin(1) * t172;
t199 = Ifges(4,3) + Ifges(5,3);
t149 = pkin(3) * t161 + pkin(10);
t198 = pkin(11) + t149;
t116 = t164 * t168 * pkin(1) + pkin(8) * t189;
t107 = pkin(9) * t164 + t116;
t108 = (-pkin(2) * t172 - pkin(9) * t168 - pkin(1)) * t162;
t63 = -t107 * t167 + t171 * t108;
t46 = -pkin(3) * t189 - qJ(4) * t114 + t63;
t64 = t171 * t107 + t167 * t108;
t51 = qJ(4) * t113 + t64;
t25 = t161 * t46 + t163 * t51;
t23 = -pkin(10) * t189 + t25;
t143 = pkin(8) * t190;
t106 = t143 + (-pkin(2) - t200) * t164;
t82 = -pkin(3) * t113 + t106;
t28 = pkin(4) * t75 - pkin(10) * t76 + t82;
t7 = t166 * t28 + t170 * t23;
t100 = -t163 * t135 + t161 * t182;
t151 = -pkin(3) * t171 - pkin(2);
t87 = pkin(4) * t127 - pkin(10) * t128 + t151;
t48 = t170 * t100 + t166 * t87;
t115 = t164 * t200 - t143;
t194 = t115 * mrSges(3,1);
t193 = t116 * mrSges(3,2);
t192 = t128 * t166;
t191 = t128 * t170;
t95 = Ifges(7,5) * t131 + Ifges(7,6) * t130;
t136 = Ifges(6,5) * t166 + Ifges(6,6) * t170;
t186 = t166 ^ 2 + t170 ^ 2;
t185 = t167 ^ 2 + t171 ^ 2;
t31 = -t165 * t57 + t169 * t56;
t32 = t165 * t56 + t169 * t57;
t8 = Ifges(7,5) * t32 + Ifges(7,6) * t31 + Ifges(7,3) * t75;
t19 = Ifges(6,5) * t57 + Ifges(6,6) * t56 + Ifges(6,3) * t75;
t77 = t131 * t128;
t78 = t130 * t128;
t35 = Ifges(7,5) * t78 - Ifges(7,6) * t77 + Ifges(7,3) * t127;
t184 = t95 / 0.2e1 + t136 / 0.2e1;
t183 = Ifges(3,5) * t190 + Ifges(3,6) * t189 + Ifges(3,3) * t164;
t150 = -pkin(3) * t163 - pkin(4);
t43 = t75 * mrSges(5,1) + t76 * mrSges(5,2);
t24 = -t161 * t51 + t163 * t46;
t6 = -t166 * t23 + t170 * t28;
t90 = t127 * mrSges(5,1) + t128 * mrSges(5,2);
t94 = -t130 * mrSges(7,1) + t131 * mrSges(7,2);
t47 = -t100 * t166 + t170 * t87;
t181 = -Ifges(4,5) * t114 - Ifges(5,5) * t76 - Ifges(4,6) * t113 + Ifges(5,6) * t75;
t22 = pkin(4) * t189 - t24;
t133 = -t170 * mrSges(6,1) + t166 * mrSges(6,2);
t180 = mrSges(6,1) * t166 + mrSges(6,2) * t170;
t122 = t198 * t166;
t123 = t198 * t170;
t84 = -t122 * t169 - t123 * t165;
t85 = -t122 * t165 + t123 * t169;
t179 = t84 * mrSges(7,1) - t85 * mrSges(7,2) + t95;
t58 = Ifges(6,5) * t191 - Ifges(6,6) * t192 + Ifges(6,3) * t127;
t4 = pkin(5) * t75 - pkin(11) * t57 + t6;
t5 = pkin(11) * t56 + t7;
t2 = -t165 * t5 + t169 * t4;
t3 = t165 * t4 + t169 * t5;
t178 = t2 * mrSges(7,1) - t3 * mrSges(7,2) + t8;
t34 = pkin(5) * t127 - pkin(11) * t191 + t47;
t42 = -pkin(11) * t192 + t48;
t14 = -t165 * t42 + t169 * t34;
t15 = t165 * t34 + t169 * t42;
t177 = t14 * mrSges(7,1) - t15 * mrSges(7,2) + t35;
t176 = (mrSges(7,1) * t169 - mrSges(7,2) * t165) * pkin(5);
t140 = Ifges(4,1) * t167 + Ifges(4,4) * t171;
t138 = Ifges(4,4) * t167 + Ifges(4,2) * t171;
t137 = Ifges(6,2) * t170 + t196;
t134 = -mrSges(4,1) * t171 + mrSges(4,2) * t167;
t132 = -pkin(5) * t170 + t150;
t102 = -mrSges(4,1) * t189 - mrSges(4,3) * t114;
t101 = mrSges(4,2) * t189 + mrSges(4,3) * t113;
t92 = Ifges(5,1) * t128 - Ifges(5,4) * t127;
t91 = Ifges(5,4) * t128 - Ifges(5,2) * t127;
t89 = mrSges(6,1) * t127 - mrSges(6,3) * t191;
t88 = -mrSges(6,2) * t127 - mrSges(6,3) * t192;
t86 = -mrSges(4,1) * t113 + mrSges(4,2) * t114;
t83 = t180 * t128;
t67 = Ifges(4,1) * t114 + Ifges(4,4) * t113 - Ifges(4,5) * t189;
t66 = Ifges(4,4) * t114 + Ifges(4,2) * t113 - Ifges(4,6) * t189;
t65 = pkin(5) * t192 + t98;
t62 = -mrSges(5,1) * t189 - mrSges(5,3) * t76;
t61 = mrSges(5,2) * t189 - mrSges(5,3) * t75;
t53 = mrSges(7,1) * t127 - mrSges(7,3) * t78;
t52 = -mrSges(7,2) * t127 - mrSges(7,3) * t77;
t44 = mrSges(7,1) * t77 + mrSges(7,2) * t78;
t41 = Ifges(5,1) * t76 - Ifges(5,4) * t75 - Ifges(5,5) * t189;
t40 = Ifges(5,4) * t76 - Ifges(5,2) * t75 - Ifges(5,6) * t189;
t39 = mrSges(6,1) * t75 - mrSges(6,3) * t57;
t38 = -mrSges(6,2) * t75 + mrSges(6,3) * t56;
t37 = Ifges(7,1) * t78 - Ifges(7,4) * t77 + Ifges(7,5) * t127;
t36 = Ifges(7,4) * t78 - Ifges(7,2) * t77 + Ifges(7,6) * t127;
t33 = -mrSges(6,1) * t56 + mrSges(6,2) * t57;
t20 = Ifges(6,4) * t57 + Ifges(6,2) * t56 + Ifges(6,6) * t75;
t17 = mrSges(7,1) * t75 - mrSges(7,3) * t32;
t16 = -mrSges(7,2) * t75 + mrSges(7,3) * t31;
t12 = -pkin(5) * t56 + t22;
t11 = -mrSges(7,1) * t31 + mrSges(7,2) * t32;
t10 = Ifges(7,1) * t32 + Ifges(7,4) * t31 + Ifges(7,5) * t75;
t9 = Ifges(7,4) * t32 + Ifges(7,2) * t31 + Ifges(7,6) * t75;
t1 = [((-0.2e1 * t115 * mrSges(3,3) + Ifges(3,5) * t164 + (-0.2e1 * pkin(1) * mrSges(3,2) + Ifges(3,1) * t168) * t162) * t168 + (0.2e1 * t116 * mrSges(3,3) + Ifges(3,6) * t164 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t168 + (Ifges(3,2) + t199) * t172) * t162 + t181) * t172) * t162 + (t183 - 0.2e1 * t193 + 0.2e1 * t194) * t164 + Ifges(2,3) + 0.2e1 * t12 * t11 + 0.2e1 * t3 * t16 + 0.2e1 * t2 * t17 + t31 * t9 + t32 * t10 + 0.2e1 * t22 * t33 + 0.2e1 * t7 * t38 + 0.2e1 * t6 * t39 + t56 * t20 + t57 * t21 + 0.2e1 * t25 * t61 + 0.2e1 * t24 * t62 + t76 * t41 + 0.2e1 * t82 * t43 + 0.2e1 * t64 * t101 + 0.2e1 * t63 * t102 + 0.2e1 * t106 * t86 + t113 * t66 + t114 * t67 + (t8 + t19 - t40) * t75 + m(7) * (t12 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t22 ^ 2 + t6 ^ 2 + t7 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2 + t82 ^ 2) + m(4) * (t106 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(3) * (pkin(1) ^ 2 * t162 ^ 2 + t115 ^ 2 + t116 ^ 2); (-t24 * mrSges(5,3) - t166 * t20 / 0.2e1 + t170 * t209 + t41 / 0.2e1) * t128 - t212 * t189 / 0.2e1 + t183 + (-t25 * mrSges(5,3) + t8 / 0.2e1 + t19 / 0.2e1 - t40 / 0.2e1) * t127 + (t33 - t62) * t98 + t57 * t207 + t56 * t208 + (t64 * mrSges(4,3) + pkin(9) * t101 + t66 / 0.2e1) * t171 + (-t63 * mrSges(4,3) - pkin(9) * t102 + t67 / 0.2e1) * t167 + m(7) * (t12 * t65 + t14 * t2 + t15 * t3) + m(6) * (t22 * t98 + t47 * t6 + t48 * t7) + m(5) * (t100 * t25 + t151 * t82 - t24 * t98) + (t35 / 0.2e1 + t58 / 0.2e1 - t91 / 0.2e1) * t75 - t193 + t194 + m(4) * (-pkin(2) * t106 + (-t63 * t167 + t64 * t171) * pkin(9)) + t15 * t16 + t14 * t17 + t31 * t36 / 0.2e1 + t32 * t37 / 0.2e1 + t12 * t44 + t47 * t39 + t48 * t38 + t3 * t52 + t2 * t53 + t65 * t11 - t77 * t9 / 0.2e1 + t78 * t10 / 0.2e1 + t22 * t83 - pkin(2) * t86 + t7 * t88 + t6 * t89 + t82 * t90 + t76 * t92 / 0.2e1 + t100 * t61 + t106 * t134 + t113 * t138 / 0.2e1 + t114 * t140 / 0.2e1 + t151 * t43; -0.2e1 * pkin(2) * t134 + t171 * t138 + 0.2e1 * t14 * t53 + t167 * t140 + 0.2e1 * t15 * t52 + 0.2e1 * t151 * t90 - t77 * t36 + t78 * t37 + 0.2e1 * t65 * t44 + 0.2e1 * t47 * t89 + 0.2e1 * t48 * t88 + t83 * t210 + Ifges(3,3) + 0.2e1 * t185 * pkin(9) * mrSges(4,3) + (mrSges(5,3) * t210 - t166 * t59 + t170 * t60 + t92) * t128 + (-0.2e1 * mrSges(5,3) * t100 + t35 + t58 - t91) * t127 + m(7) * (t14 ^ 2 + t15 ^ 2 + t65 ^ 2) + m(6) * (t47 ^ 2 + t48 ^ 2 + t211) + m(5) * (t100 ^ 2 + t151 ^ 2 + t211) + m(4) * (t185 * pkin(9) ^ 2 + pkin(2) ^ 2); (-t6 * mrSges(6,3) - t149 * t39 + t209) * t166 - t199 * t189 - t181 + (t3 * t130 - t2 * t131) * mrSges(7,3) + t9 * t204 + t32 * t205 + t31 * t206 + t57 * t202 + t10 * t203 + t184 * t75 + (t7 * mrSges(6,3) + t149 * t38 + t20 / 0.2e1) * t170 + m(7) * (t12 * t132 + t2 * t84 + t3 * t85) + (t161 * t61 + t163 * t62 + m(5) * (t161 * t25 + t163 * t24)) * pkin(3) + m(6) * (t150 * t22 + (-t6 * t166 + t7 * t170) * t149) + t24 * mrSges(5,1) - t25 * mrSges(5,2) + t63 * mrSges(4,1) - t64 * mrSges(4,2) + t84 * t17 + t85 * t16 + t12 * t94 + t132 * t11 + t22 * t133 + t56 * t137 / 0.2e1 + t150 * t33; (-mrSges(5,1) + t133) * t98 + (t48 * mrSges(6,3) + t128 * t202 + t149 * t88 + t208) * t170 + (-t47 * mrSges(6,3) - t128 * t137 / 0.2e1 - t149 * t89 + t207) * t166 + t184 * t127 + m(7) * (t132 * t65 + t14 * t84 + t15 * t85) + (-mrSges(4,1) * t167 - mrSges(4,2) * t171) * pkin(9) + (t130 * t15 - t131 * t14) * mrSges(7,3) + (m(5) * (t100 * t161 - t163 * t98) + (-t127 * t161 - t128 * t163) * mrSges(5,3)) * pkin(3) + m(6) * (t150 * t98 + (-t166 * t47 + t170 * t48) * t149) + t84 * t53 + t85 * t52 + t65 * t94 - t77 * t206 + t78 * t205 - t100 * mrSges(5,2) + t36 * t204 + t37 * t203 + t132 * t44 + t150 * t83 + t212; t130 * t96 + t131 * t97 + 0.2e1 * t132 * t94 + 0.2e1 * t150 * t133 + t170 * t137 + t166 * t139 + m(7) * (t132 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(6) * (t186 * t149 ^ 2 + t150 ^ 2) + m(5) * (t161 ^ 2 + t163 ^ 2) * pkin(3) ^ 2 + t199 + 0.2e1 * (mrSges(5,1) * t163 - mrSges(5,2) * t161) * pkin(3) + 0.2e1 * (t130 * t85 - t131 * t84) * mrSges(7,3) + 0.2e1 * t186 * t149 * mrSges(6,3); t130 * t17 + t131 * t16 + t166 * t38 + t170 * t39 + m(7) * (t130 * t2 + t131 * t3) + m(6) * (t166 * t7 + t170 * t6) + m(5) * t82 + t43; t130 * t53 + t131 * t52 + t166 * t88 + t170 * t89 + m(7) * (t130 * t14 + t131 * t15) + m(6) * (t166 * t48 + t170 * t47) + m(5) * t151 + t90; m(7) * (t130 * t84 + t131 * t85); m(5) + m(6) * t186 + m(7) * (t130 ^ 2 + t131 ^ 2); t6 * mrSges(6,1) - t7 * mrSges(6,2) + (m(7) * (t165 * t3 + t169 * t2) + t165 * t16 + t169 * t17) * pkin(5) + t178 + t19; t47 * mrSges(6,1) - t48 * mrSges(6,2) + (m(7) * (t14 * t169 + t15 * t165) + t165 * t52 + t169 * t53) * pkin(5) + t177 + t58; -t180 * t149 + (m(7) * (t165 * t85 + t169 * t84) + (t130 * t165 - t131 * t169) * mrSges(7,3)) * pkin(5) + t179 + t136; m(7) * (t130 * t169 + t131 * t165) * pkin(5) - t133 - t94; Ifges(6,3) + Ifges(7,3) + m(7) * (t165 ^ 2 + t169 ^ 2) * pkin(5) ^ 2 + 0.2e1 * t176; t178; t177; t179; -t94; Ifges(7,3) + t176; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
