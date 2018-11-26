% Calculate joint inertia matrix for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR6_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_inertiaJ_slag_vp2: pkin has to be [14x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRR6_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRR6_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:36:16
% EndTime: 2018-11-23 15:36:18
% DurationCPUTime: 2.09s
% Computational Cost: add. (4779->470), mult. (12851->701), div. (0->0), fcn. (14703->16), ass. (0->177)
t229 = 2 * pkin(12);
t163 = sin(pkin(8));
t166 = cos(pkin(8));
t171 = sin(qJ(4));
t176 = cos(qJ(4));
t172 = sin(qJ(3));
t164 = sin(pkin(7));
t177 = cos(qJ(3));
t200 = t164 * t177;
t167 = cos(pkin(7));
t214 = pkin(2) * t167;
t120 = pkin(10) * t200 + t172 * t214;
t197 = t166 * t177;
t190 = t164 * t197;
t74 = (t163 * t167 + t190) * pkin(11) + t120;
t148 = t177 * t214;
t201 = t164 * t172;
t84 = pkin(3) * t167 + t148 + (-pkin(11) * t166 - pkin(10)) * t201;
t94 = (-pkin(11) * t163 * t172 - pkin(3) * t177 - pkin(2)) * t164;
t30 = -t171 * t74 + (t163 * t94 + t166 * t84) * t176;
t169 = sin(qJ(6));
t174 = cos(qJ(6));
t130 = -mrSges(7,1) * t174 + mrSges(7,2) * t169;
t228 = -m(7) * pkin(5) - mrSges(6,1) + t130;
t170 = sin(qJ(5));
t175 = cos(qJ(5));
t165 = sin(pkin(6));
t168 = cos(pkin(6));
t178 = cos(qJ(2));
t112 = -t164 * t165 * t178 + t167 * t168;
t173 = sin(qJ(2));
t196 = t167 * t178;
t79 = t168 * t200 + (-t172 * t173 + t177 * t196) * t165;
t81 = t168 * t201 + (t172 * t196 + t173 * t177) * t165;
t34 = t176 * t81 + (t112 * t163 + t166 * t79) * t171;
t54 = t112 * t166 - t163 * t79;
t14 = t170 * t34 - t175 * t54;
t227 = t14 ^ 2;
t198 = t166 * t176;
t202 = t163 * t176;
t32 = -t112 * t202 + t171 * t81 - t79 * t198;
t226 = t32 ^ 2;
t111 = -t163 * t200 + t166 * t167;
t203 = t163 * t171;
t80 = t167 * t203 + (t171 * t197 + t172 * t176) * t164;
t55 = -t175 * t111 + t170 * t80;
t56 = t111 * t170 + t175 * t80;
t78 = -t167 * t202 + t171 * t201 - t176 * t190;
t21 = Ifges(6,1) * t56 - Ifges(6,4) * t55 + Ifges(6,5) * t78;
t225 = t21 / 0.2e1;
t113 = -t175 * t166 + t170 * t203;
t114 = t166 * t170 + t175 * t203;
t68 = Ifges(6,1) * t114 - Ifges(6,4) * t113 - Ifges(6,5) * t202;
t224 = t68 / 0.2e1;
t205 = Ifges(7,4) * t174;
t105 = -Ifges(7,6) * t175 + (-Ifges(7,2) * t169 + t205) * t170;
t223 = t105 / 0.2e1;
t206 = Ifges(7,4) * t169;
t106 = -Ifges(7,5) * t175 + (Ifges(7,1) * t174 - t206) * t170;
t222 = t106 / 0.2e1;
t132 = Ifges(7,5) * t169 + Ifges(7,6) * t174;
t221 = t132 / 0.2e1;
t134 = Ifges(7,2) * t174 + t206;
t220 = t134 / 0.2e1;
t136 = Ifges(7,1) * t169 + t205;
t219 = t136 / 0.2e1;
t137 = Ifges(6,1) * t170 + Ifges(6,4) * t175;
t218 = t137 / 0.2e1;
t217 = -t169 / 0.2e1;
t216 = t169 / 0.2e1;
t215 = t174 / 0.2e1;
t213 = pkin(3) * t163;
t212 = pkin(12) * t170;
t211 = pkin(12) * t175;
t210 = pkin(13) * t169;
t209 = pkin(13) * t174;
t37 = -t169 * t56 + t174 * t78;
t38 = t169 * t78 + t174 * t56;
t13 = -mrSges(7,1) * t37 + mrSges(7,2) * t38;
t40 = mrSges(6,1) * t78 - mrSges(6,3) * t56;
t208 = t13 - t40;
t49 = -t163 * t84 + t166 * t94;
t23 = pkin(4) * t78 - pkin(12) * t80 + t49;
t199 = t166 * t171;
t31 = t176 * t74 + t84 * t199 + t94 * t203;
t26 = pkin(12) * t111 + t31;
t8 = t170 * t23 + t175 * t26;
t85 = -t114 * t169 - t174 * t202;
t86 = t114 * t174 - t169 * t202;
t48 = -mrSges(7,1) * t85 + mrSges(7,2) * t86;
t90 = -mrSges(6,1) * t202 - mrSges(6,3) * t114;
t207 = t48 - t90;
t119 = pkin(3) * t199 + pkin(11) * t202;
t102 = pkin(12) * t166 + t119;
t103 = (-pkin(4) * t176 - pkin(12) * t171 - pkin(3)) * t163;
t65 = t175 * t102 + t170 * t103;
t204 = t14 * t170;
t195 = t169 * t170;
t194 = t170 * t174;
t133 = Ifges(6,5) * t170 + Ifges(6,6) * t175;
t193 = t169 ^ 2 + t174 ^ 2;
t9 = Ifges(7,5) * t38 + Ifges(7,6) * t37 + Ifges(7,3) * t55;
t19 = Ifges(6,5) * t56 - Ifges(6,6) * t55 + Ifges(6,3) * t78;
t20 = Ifges(6,4) * t56 - Ifges(6,2) * t55 + Ifges(6,6) * t78;
t192 = t9 / 0.2e1 - t20 / 0.2e1;
t41 = Ifges(5,5) * t80 - Ifges(5,6) * t78 + Ifges(5,3) * t111;
t44 = Ifges(7,5) * t86 + Ifges(7,6) * t85 + Ifges(7,3) * t113;
t67 = Ifges(6,4) * t114 - Ifges(6,2) * t113 - Ifges(6,6) * t202;
t191 = t44 / 0.2e1 - t67 / 0.2e1;
t97 = Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * t166;
t189 = Ifges(4,5) * t201 + Ifges(4,6) * t200 + Ifges(4,3) * t167;
t104 = Ifges(7,5) * t194 - Ifges(7,6) * t195 - Ifges(7,3) * t175;
t135 = Ifges(6,4) * t170 + Ifges(6,2) * t175;
t188 = t104 / 0.2e1 - t135 / 0.2e1;
t25 = -pkin(4) * t111 - t30;
t12 = pkin(5) * t55 - pkin(13) * t56 + t25;
t4 = pkin(13) * t78 + t8;
t1 = t12 * t174 - t169 * t4;
t2 = t12 * t169 + t174 * t4;
t187 = -t1 * t169 + t174 * t2;
t185 = mrSges(7,1) * t169 + mrSges(7,2) * t174;
t16 = t170 * t54 + t175 * t34;
t184 = t16 * t175 + t204;
t143 = pkin(11) * t203;
t101 = t143 + (-pkin(3) * t176 - pkin(4)) * t166;
t57 = pkin(5) * t113 - pkin(13) * t114 + t101;
t59 = -pkin(13) * t202 + t65;
t27 = -t169 * t59 + t174 * t57;
t28 = t169 * t57 + t174 * t59;
t182 = -t169 * t27 + t174 * t28;
t129 = -pkin(5) * t175 - pkin(13) * t170 - pkin(4);
t95 = t129 * t174 - t169 * t211;
t96 = t129 * t169 + t174 * t211;
t181 = -t169 * t95 + t174 * t96;
t7 = -t170 * t26 + t175 * t23;
t64 = -t102 * t170 + t103 * t175;
t66 = Ifges(6,5) * t114 - Ifges(6,6) * t113 - Ifges(6,3) * t202;
t180 = pkin(12) ^ 2;
t162 = t175 ^ 2;
t160 = t170 ^ 2;
t157 = t160 * t180;
t131 = -mrSges(6,1) * t175 + mrSges(6,2) * t170;
t127 = -mrSges(7,1) * t175 - mrSges(7,3) * t194;
t126 = mrSges(7,2) * t175 - mrSges(7,3) * t195;
t125 = -mrSges(4,2) * t167 + mrSges(4,3) * t200;
t124 = mrSges(4,1) * t167 - mrSges(4,3) * t201;
t123 = -mrSges(5,2) * t166 + mrSges(5,3) * t202;
t122 = mrSges(5,1) * t166 - mrSges(5,3) * t203;
t121 = t185 * t170;
t118 = -pkin(10) * t201 + t148;
t117 = pkin(3) * t198 - t143;
t116 = (-mrSges(4,1) * t177 + mrSges(4,2) * t172) * t164;
t115 = (-mrSges(5,1) * t176 + mrSges(5,2) * t171) * t163;
t99 = Ifges(5,5) * t166 + (Ifges(5,1) * t171 + Ifges(5,4) * t176) * t163;
t98 = Ifges(5,6) * t166 + (Ifges(5,4) * t171 + Ifges(5,2) * t176) * t163;
t89 = mrSges(6,2) * t202 - mrSges(6,3) * t113;
t70 = mrSges(6,1) * t113 + mrSges(6,2) * t114;
t63 = mrSges(7,1) * t113 - mrSges(7,3) * t86;
t62 = -mrSges(7,2) * t113 + mrSges(7,3) * t85;
t61 = mrSges(5,1) * t111 - mrSges(5,3) * t80;
t60 = -mrSges(5,2) * t111 - mrSges(5,3) * t78;
t58 = pkin(5) * t202 - t64;
t47 = mrSges(5,1) * t78 + mrSges(5,2) * t80;
t46 = Ifges(7,1) * t86 + Ifges(7,4) * t85 + Ifges(7,5) * t113;
t45 = Ifges(7,4) * t86 + Ifges(7,2) * t85 + Ifges(7,6) * t113;
t43 = Ifges(5,1) * t80 - Ifges(5,4) * t78 + Ifges(5,5) * t111;
t42 = Ifges(5,4) * t80 - Ifges(5,2) * t78 + Ifges(5,6) * t111;
t39 = -mrSges(6,2) * t78 - mrSges(6,3) * t55;
t29 = mrSges(6,1) * t55 + mrSges(6,2) * t56;
t18 = mrSges(7,1) * t55 - mrSges(7,3) * t38;
t17 = -mrSges(7,2) * t55 + mrSges(7,3) * t37;
t11 = Ifges(7,1) * t38 + Ifges(7,4) * t37 + Ifges(7,5) * t55;
t10 = Ifges(7,4) * t38 + Ifges(7,2) * t37 + Ifges(7,6) * t55;
t6 = t16 * t174 + t169 * t32;
t5 = -t16 * t169 + t174 * t32;
t3 = -pkin(5) * t78 - t7;
t15 = [m(2) + m(7) * (t5 ^ 2 + t6 ^ 2 + t227) + m(6) * (t16 ^ 2 + t226 + t227) + m(5) * (t34 ^ 2 + t54 ^ 2 + t226) + m(4) * (t112 ^ 2 + t79 ^ 2 + t81 ^ 2) + m(3) * (t168 ^ 2 + (t173 ^ 2 + t178 ^ 2) * t165 ^ 2); t112 * t116 + t79 * t124 + t81 * t125 + t16 * t39 + t6 * t17 + t5 * t18 + t34 * t60 + t54 * t47 + (t29 - t61) * t32 + (mrSges(3,1) * t178 - mrSges(3,2) * t173) * t165 + t208 * t14 + m(7) * (t1 * t5 + t14 * t3 + t2 * t6) + m(6) * (-t14 * t7 + t16 * t8 + t25 * t32) + m(5) * (-t30 * t32 + t31 * t34 + t49 * t54) + m(4) * (-pkin(2) * t112 * t164 + t118 * t79 + t120 * t81); m(6) * (t25 ^ 2 + t7 ^ 2 + t8 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(5) * (t30 ^ 2 + t31 ^ 2 + t49 ^ 2) + (-0.2e1 * pkin(2) * t116 + (Ifges(4,5) * t172 + Ifges(4,6) * t177) * t167 + (m(4) * pkin(2) ^ 2 + Ifges(4,2) * t177 ^ 2 + (Ifges(4,1) * t172 + 0.2e1 * Ifges(4,4) * t177) * t172) * t164) * t164 + (t19 - t42) * t78 + (t9 - t20) * t55 + t167 * t189 + m(4) * (t118 ^ 2 + t120 ^ 2) + Ifges(3,3) + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t17 + 0.2e1 * t1 * t18 + 0.2e1 * t25 * t29 + t37 * t10 + t38 * t11 + 0.2e1 * t8 * t39 + 0.2e1 * t7 * t40 + 0.2e1 * t49 * t47 + t56 * t21 + 0.2e1 * t31 * t60 + 0.2e1 * t30 * t61 + t80 * t43 + t111 * t41 + 0.2e1 * t118 * t124 + 0.2e1 * t120 * t125; t79 * mrSges(4,1) - t81 * mrSges(4,2) + t54 * t115 + t34 * t123 + t16 * t89 + t5 * t63 + t6 * t62 + (-t122 + t70) * t32 + t207 * t14 + m(7) * (t14 * t58 + t27 * t5 + t28 * t6) + m(6) * (t101 * t32 - t14 * t64 + t16 * t65) + m(5) * (-t117 * t32 + t119 * t34 - t54 * t213); t191 * t55 + m(7) * (t1 * t27 + t2 * t28 + t3 * t58) + m(6) * (t101 * t25 + t64 * t7 + t65 * t8) + (t66 / 0.2e1 - t98 / 0.2e1) * t78 + m(5) * (t117 * t30 + t119 * t31 - t49 * t213) + (-pkin(3) * t47 + t171 * t43 / 0.2e1 + (-t19 / 0.2e1 + t42 / 0.2e1) * t176) * t163 + t189 + t56 * t224 + t114 * t225 + t27 * t18 + t28 * t17 + t37 * t45 / 0.2e1 + t38 * t46 / 0.2e1 + t3 * t48 + t58 * t13 + t2 * t62 + t1 * t63 + t64 * t40 + t65 * t39 + t25 * t70 + t85 * t10 / 0.2e1 + t86 * t11 / 0.2e1 + t8 * t89 + t7 * t90 + t80 * t99 / 0.2e1 + t101 * t29 + t111 * t97 / 0.2e1 + t49 * t115 + t117 * t61 + t118 * mrSges(4,1) + t119 * t60 - t120 * mrSges(4,2) + t30 * t122 + t31 * t123 + t166 * t41 / 0.2e1 + t192 * t113; 0.2e1 * t101 * t70 + t114 * t68 + 0.2e1 * t117 * t122 + 0.2e1 * t119 * t123 + t166 * t97 + 0.2e1 * t27 * t63 + 0.2e1 * t28 * t62 + t85 * t45 + t86 * t46 + 0.2e1 * t58 * t48 + 0.2e1 * t64 * t90 + 0.2e1 * t65 * t89 + Ifges(4,3) + (t44 - t67) * t113 + (-0.2e1 * pkin(3) * t115 + t171 * t99 + (-t66 + t98) * t176) * t163 + m(7) * (t27 ^ 2 + t28 ^ 2 + t58 ^ 2) + m(6) * (t101 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(5) * (pkin(3) ^ 2 * t163 ^ 2 + t117 ^ 2 + t119 ^ 2); -t34 * mrSges(5,2) + t14 * t121 + t6 * t126 + t5 * t127 + (-mrSges(5,1) + t131) * t32 + t184 * mrSges(6,3) + m(7) * (pkin(12) * t204 + t5 * t95 + t6 * t96) + m(6) * (-pkin(4) * t32 + t184 * pkin(12)); m(6) * (-pkin(4) * t25 + (-t7 * t170 + t8 * t175) * pkin(12)) + t188 * t55 + t41 + m(7) * (t1 * t95 + t2 * t96 + t3 * t212) + (t8 * mrSges(6,3) + pkin(12) * t39 - t192) * t175 + (-t7 * mrSges(6,3) + t208 * pkin(12) + t10 * t217 + t11 * t215 + t225) * t170 - pkin(4) * t29 + t30 * mrSges(5,1) - t31 * mrSges(5,2) + t95 * t18 + t96 * t17 + t37 * t223 + t38 * t222 + t3 * t121 + t2 * t126 + t1 * t127 + t25 * t131 + t78 * t133 / 0.2e1 + t56 * t218; m(6) * (-pkin(4) * t101 + (-t64 * t170 + t65 * t175) * pkin(12)) + t188 * t113 + (-t64 * mrSges(6,3) + t207 * pkin(12) + t46 * t215 + t45 * t217 + t224) * t170 + t97 - t133 * t202 / 0.2e1 - pkin(4) * t70 + t95 * t63 + t96 * t62 + t85 * t223 + t86 * t222 + t117 * mrSges(5,1) - t119 * mrSges(5,2) + t58 * t121 + t28 * t126 + t27 * t127 + t101 * t131 + t114 * t218 + m(7) * (t58 * t212 + t27 * t95 + t28 * t96) + (t65 * mrSges(6,3) + pkin(12) * t89 - t191) * t175; -0.2e1 * pkin(4) * t131 + 0.2e1 * t96 * t126 + 0.2e1 * t95 * t127 + Ifges(5,3) + (-t104 + t135) * t175 + (t160 + t162) * mrSges(6,3) * t229 + m(7) * (t95 ^ 2 + t96 ^ 2 + t157) + m(6) * (pkin(4) ^ 2 + t162 * t180 + t157) + (-t105 * t169 + t106 * t174 + t121 * t229 + t137) * t170; -t16 * mrSges(6,2) + (m(7) * pkin(13) + mrSges(7,3)) * (-t169 * t5 + t174 * t6) + t228 * t14; m(7) * (-pkin(5) * t3 + t187 * pkin(13)) - pkin(5) * t13 + t11 * t216 - t18 * t210 + t17 * t209 + t3 * t130 + t10 * t215 + t55 * t221 + t38 * t219 + t37 * t220 - t8 * mrSges(6,2) + t7 * mrSges(6,1) + t187 * mrSges(7,3) + t19; t58 * t130 + t45 * t215 + t46 * t216 - pkin(5) * t48 + m(7) * (-pkin(5) * t58 + t182 * pkin(13)) + t62 * t209 - t63 * t210 + t85 * t220 + t113 * t221 + t86 * t219 - t65 * mrSges(6,2) + t64 * mrSges(6,1) + t182 * mrSges(7,3) + t66; t105 * t215 - pkin(5) * t121 + t106 * t216 + (m(7) * t181 + t174 * t126 - t169 * t127) * pkin(13) + (pkin(12) * t228 + t134 * t217 + t136 * t215) * t170 + (-pkin(12) * mrSges(6,2) - t132 / 0.2e1) * t175 + t181 * mrSges(7,3) + t133; Ifges(6,3) - 0.2e1 * pkin(5) * t130 + t174 * t134 + t169 * t136 + m(7) * (t193 * pkin(13) ^ 2 + pkin(5) ^ 2) + 0.2e1 * t193 * pkin(13) * mrSges(7,3); mrSges(7,1) * t5 - mrSges(7,2) * t6; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t9; mrSges(7,1) * t27 - mrSges(7,2) * t28 + t44; mrSges(7,1) * t95 - mrSges(7,2) * t96 + t104; -t185 * pkin(13) + t132; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
