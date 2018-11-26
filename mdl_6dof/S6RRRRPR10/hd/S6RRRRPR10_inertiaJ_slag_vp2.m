% Calculate joint inertia matrix for
% S6RRRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:20
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR10_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR10_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR10_inertiaJ_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR10_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR10_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR10_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:19:19
% EndTime: 2018-11-23 18:19:21
% DurationCPUTime: 2.05s
% Computational Cost: add. (3218->376), mult. (7010->520), div. (0->0), fcn. (7636->10), ass. (0->146)
t212 = Ifges(6,4) - Ifges(5,5);
t211 = Ifges(6,5) - Ifges(5,6);
t144 = cos(pkin(6));
t147 = sin(qJ(3));
t151 = cos(qJ(3));
t143 = sin(pkin(6));
t148 = sin(qJ(2));
t186 = t143 * t148;
t100 = t144 * t151 - t147 * t186;
t101 = t144 * t147 + t151 * t186;
t210 = -Ifges(4,5) * t101 - Ifges(4,6) * t100;
t145 = sin(qJ(6));
t149 = cos(qJ(6));
t179 = -t145 ^ 2 - t149 ^ 2;
t173 = t179 * mrSges(7,3);
t114 = mrSges(7,1) * t145 + mrSges(7,2) * t149;
t209 = mrSges(6,3) + t114;
t152 = cos(qJ(2));
t185 = t143 * t152;
t146 = sin(qJ(4));
t150 = cos(qJ(4));
t110 = t146 * t147 - t150 * t151;
t111 = t146 * t151 + t147 * t150;
t208 = t211 * t110 - t212 * t111;
t129 = pkin(3) * t146 + qJ(5);
t207 = t129 ^ 2;
t206 = pkin(4) + pkin(11);
t205 = -pkin(10) - pkin(9);
t136 = Ifges(7,5) * t149;
t115 = -Ifges(7,6) * t145 + t136;
t204 = t115 / 0.2e1;
t192 = Ifges(7,4) * t149;
t116 = -Ifges(7,2) * t145 + t192;
t203 = t116 / 0.2e1;
t193 = Ifges(7,4) * t145;
t118 = Ifges(7,1) * t149 - t193;
t202 = t118 / 0.2e1;
t201 = -t145 / 0.2e1;
t200 = t149 / 0.2e1;
t198 = pkin(1) * t152;
t196 = Ifges(6,1) + Ifges(5,3);
t103 = t144 * t148 * pkin(1) + pkin(8) * t185;
t93 = pkin(9) * t144 + t103;
t94 = (-pkin(2) * t152 - pkin(9) * t148 - pkin(1)) * t143;
t53 = -t147 * t93 + t151 * t94;
t31 = -pkin(3) * t185 - pkin(10) * t101 + t53;
t54 = t147 * t94 + t151 * t93;
t34 = pkin(10) * t100 + t54;
t17 = t146 * t31 + t150 * t34;
t64 = -t150 * t100 + t101 * t146;
t41 = t145 * t185 + t149 * t64;
t42 = t145 * t64 - t149 * t185;
t19 = -mrSges(7,1) * t41 + mrSges(7,2) * t42;
t48 = mrSges(6,1) * t64 + mrSges(6,3) * t185;
t195 = -t48 + t19;
t194 = mrSges(7,1) * t149;
t124 = pkin(8) * t186;
t102 = t144 * t198 - t124;
t191 = t102 * mrSges(3,1);
t190 = t103 * mrSges(3,2);
t189 = qJ(5) * t129;
t188 = t110 * t145;
t187 = t110 * t149;
t181 = Ifges(4,5) * t147 + Ifges(4,6) * t151;
t180 = t147 ^ 2 + t151 ^ 2;
t178 = Ifges(4,3) + t196;
t65 = t100 * t146 + t101 * t150;
t12 = Ifges(7,5) * t42 + Ifges(7,6) * t41 + Ifges(7,3) * t65;
t45 = Ifges(7,5) * t188 + Ifges(7,6) * t187 + Ifges(7,3) * t111;
t120 = t205 * t151;
t176 = t205 * t147;
t84 = -t120 * t146 - t150 * t176;
t86 = -t150 * t120 + t146 * t176;
t177 = t84 ^ 2 + t86 ^ 2;
t175 = Ifges(3,5) * t186 + Ifges(3,6) * t185 + Ifges(3,3) * t144;
t133 = -pkin(3) * t151 - pkin(2);
t132 = -pkin(3) * t150 - pkin(4);
t174 = m(7) * t179;
t16 = -t146 * t34 + t150 * t31;
t172 = -t211 * t64 + t212 * t65;
t171 = t179 * t206;
t170 = 0.2e1 * t209;
t11 = pkin(4) * t185 - t16;
t4 = pkin(5) * t65 + pkin(11) * t185 + t11;
t92 = t124 + (-pkin(2) - t198) * t144;
t66 = -pkin(3) * t100 + t92;
t159 = -qJ(5) * t65 + t66;
t6 = t206 * t64 + t159;
t1 = -t145 * t6 + t149 * t4;
t2 = t145 * t4 + t149 * t6;
t169 = t1 * t149 + t145 * t2;
t168 = mrSges(6,2) + t173;
t167 = -t145 * mrSges(7,2) + t194;
t161 = -qJ(5) * t111 + t133;
t52 = t110 * t206 + t161;
t62 = pkin(5) * t111 + t84;
t20 = -t145 * t52 + t149 * t62;
t21 = t145 * t62 + t149 * t52;
t166 = t145 * t21 + t149 * t20;
t22 = -mrSges(7,2) * t65 + mrSges(7,3) * t41;
t23 = mrSges(7,1) * t65 - mrSges(7,3) * t42;
t165 = t145 * t22 + t149 * t23;
t70 = mrSges(7,1) * t111 - mrSges(7,3) * t188;
t71 = -mrSges(7,2) * t111 + mrSges(7,3) * t187;
t164 = t145 * t71 + t149 * t70;
t163 = 0.2e1 * t173;
t10 = qJ(5) * t185 - t17;
t162 = -t145 * t116 + t149 * t118 + t196;
t49 = t65 * mrSges(6,1) - mrSges(6,2) * t185;
t160 = (mrSges(5,1) * t150 - mrSges(5,2) * t146) * pkin(3);
t13 = Ifges(7,4) * t42 + Ifges(7,2) * t41 + Ifges(7,6) * t65;
t14 = Ifges(7,1) * t42 + Ifges(7,4) * t41 + Ifges(7,5) * t65;
t5 = -pkin(5) * t64 - t10;
t158 = t16 * mrSges(5,1) - t17 * mrSges(5,2) + t11 * mrSges(6,2) - t10 * mrSges(6,3) - mrSges(7,3) * t169 + t5 * t114 + t13 * t201 + t14 * t200 + t42 * t202 + t41 * t203 + t65 * t204 - t172;
t46 = Ifges(7,6) * t111 + (Ifges(7,2) * t149 + t193) * t110;
t47 = Ifges(7,5) * t111 + (Ifges(7,1) * t145 + t192) * t110;
t63 = -t110 * pkin(5) + t86;
t157 = t46 * t201 + t47 * t200 + t63 * t114 + t188 * t202 + t187 * t203 + t111 * t204 + (-mrSges(5,2) + mrSges(6,3)) * t86 - t166 * mrSges(7,3) + t208 + (mrSges(6,2) - mrSges(5,1)) * t84;
t154 = qJ(5) ^ 2;
t128 = -pkin(11) + t132;
t119 = Ifges(4,1) * t147 + Ifges(4,4) * t151;
t117 = Ifges(4,4) * t147 + Ifges(4,2) * t151;
t113 = -mrSges(4,1) * t151 + mrSges(4,2) * t147;
t83 = -mrSges(4,1) * t185 - mrSges(4,3) * t101;
t82 = mrSges(4,2) * t185 + mrSges(4,3) * t100;
t77 = Ifges(5,1) * t111 - Ifges(5,4) * t110;
t76 = Ifges(5,4) * t111 - Ifges(5,2) * t110;
t75 = -Ifges(6,2) * t111 + Ifges(6,6) * t110;
t74 = -Ifges(6,6) * t111 + Ifges(6,3) * t110;
t73 = -mrSges(6,2) * t110 - mrSges(6,3) * t111;
t72 = mrSges(5,1) * t110 + mrSges(5,2) * t111;
t69 = pkin(4) * t110 + t161;
t68 = t167 * t110;
t67 = -mrSges(4,1) * t100 + mrSges(4,2) * t101;
t56 = Ifges(4,1) * t101 + Ifges(4,4) * t100 - Ifges(4,5) * t185;
t55 = Ifges(4,4) * t101 + Ifges(4,2) * t100 - Ifges(4,6) * t185;
t51 = -mrSges(5,1) * t185 - mrSges(5,3) * t65;
t50 = mrSges(5,2) * t185 - mrSges(5,3) * t64;
t29 = -mrSges(6,2) * t64 - mrSges(6,3) * t65;
t28 = mrSges(5,1) * t64 + mrSges(5,2) * t65;
t27 = Ifges(5,1) * t65 - Ifges(5,4) * t64 - Ifges(5,5) * t185;
t26 = Ifges(5,4) * t65 - Ifges(5,2) * t64 - Ifges(5,6) * t185;
t25 = -Ifges(6,4) * t185 - Ifges(6,2) * t65 + Ifges(6,6) * t64;
t24 = -Ifges(6,5) * t185 - Ifges(6,6) * t65 + Ifges(6,3) * t64;
t18 = pkin(4) * t64 + t159;
t3 = [t100 * t55 + t101 * t56 + 0.2e1 * t92 * t67 + 0.2e1 * t54 * t82 + 0.2e1 * t53 * t83 + 0.2e1 * t66 * t28 + ((-0.2e1 * t102 * mrSges(3,3) + Ifges(3,5) * t144 + (-0.2e1 * mrSges(3,2) * pkin(1) + Ifges(3,1) * t148) * t143) * t148 + (0.2e1 * t103 * mrSges(3,3) + Ifges(3,6) * t144 + (0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * Ifges(3,4) * t148 + (Ifges(3,2) + t178) * t152) * t143 + t172 + t210) * t152) * t143 + t41 * t13 + t42 * t14 + 0.2e1 * t10 * t48 + 0.2e1 * t11 * t49 + 0.2e1 * t17 * t50 + 0.2e1 * t16 * t51 + 0.2e1 * t18 * t29 + 0.2e1 * t5 * t19 + 0.2e1 * t2 * t22 + 0.2e1 * t1 * t23 + m(3) * (pkin(1) ^ 2 * t143 ^ 2 + t102 ^ 2 + t103 ^ 2) + m(4) * (t53 ^ 2 + t54 ^ 2 + t92 ^ 2) + m(5) * (t16 ^ 2 + t17 ^ 2 + t66 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2 + t18 ^ 2) + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) + (t27 + t12 - t25) * t65 + (t24 - t26) * t64 + Ifges(2,3) + (t175 - 0.2e1 * t190 + 0.2e1 * t191) * t144; (t27 / 0.2e1 + t12 / 0.2e1 - t25 / 0.2e1 - t16 * mrSges(5,3) + t11 * mrSges(6,1)) * t111 + t133 * t28 + t100 * t117 / 0.2e1 + t101 * t119 / 0.2e1 + t92 * t113 - t5 * t68 + t69 * t29 + t1 * t70 + t2 * t71 + t66 * t72 + t18 * t73 - pkin(2) * t67 + t175 + m(5) * (t133 * t66 - t16 * t84 + t17 * t86) + m(6) * (-t10 * t86 + t11 * t84 + t18 * t69) + m(7) * (t1 * t20 + t2 * t21 + t5 * t63) + (-t75 / 0.2e1 + t77 / 0.2e1 + t45 / 0.2e1) * t65 + (t74 / 0.2e1 - t76 / 0.2e1) * t64 + t63 * t19 + t41 * t46 / 0.2e1 + t42 * t47 / 0.2e1 + t21 * t22 + t20 * t23 - (t181 + t208) * t185 / 0.2e1 + m(4) * (-pkin(2) * t92 + (-t53 * t147 + t54 * t151) * pkin(9)) - t190 + t191 + (t24 / 0.2e1 - t26 / 0.2e1 + t145 * t14 / 0.2e1 + t10 * mrSges(6,1) - t17 * mrSges(5,3) + t13 * t200) * t110 + (t56 / 0.2e1 - pkin(9) * t83 - t53 * mrSges(4,3)) * t147 + (t55 / 0.2e1 + pkin(9) * t82 + t54 * mrSges(4,3)) * t151 + (t49 - t51) * t84 + (t50 - t48) * t86; -0.2e1 * pkin(2) * t113 + t151 * t117 + t147 * t119 + 0.2e1 * t133 * t72 + 0.2e1 * t20 * t70 + 0.2e1 * t21 * t71 - 0.2e1 * t63 * t68 + 0.2e1 * t69 * t73 + Ifges(3,3) + (t45 - t75 + t77) * t111 + m(5) * (t133 ^ 2 + t177) + m(6) * (t69 ^ 2 + t177) + m(7) * (t20 ^ 2 + t21 ^ 2 + t63 ^ 2) + m(4) * (pkin(9) ^ 2 * t180 + pkin(2) ^ 2) + (t145 * t47 + t149 * t46 + t74 - t76) * t110 + 0.2e1 * t180 * pkin(9) * mrSges(4,3) + 0.2e1 * (-t110 * t86 + t111 * t84) * (mrSges(6,1) + mrSges(5,3)); t132 * t49 + m(6) * (-t10 * t129 + t11 * t132) + t158 + t53 * mrSges(4,1) - t54 * mrSges(4,2) + (t146 * t50 + m(5) * (t146 * t17 + t150 * t16) + t150 * t51) * pkin(3) + m(7) * (t128 * t169 + t129 * t5) + t195 * t129 + t165 * t128 - t178 * t185 - t210; -t129 * t68 + t157 + t164 * t128 + m(6) * (t129 * t86 + t132 * t84) + (m(5) * (t146 * t86 - t150 * t84) + (-t110 * t146 - t111 * t150) * mrSges(5,3)) * pkin(3) + (-mrSges(4,1) * t147 - mrSges(4,2) * t151) * pkin(9) + m(7) * (t128 * t166 + t129 * t63) + (-t110 * t129 + t111 * t132) * mrSges(6,1) + t181; 0.2e1 * t132 * mrSges(6,2) + Ifges(4,3) + t129 * t170 + 0.2e1 * t160 + t128 * t163 + m(7) * (-t128 ^ 2 * t179 + t207) + m(6) * (t132 ^ 2 + t207) + m(5) * (t146 ^ 2 + t150 ^ 2) * pkin(3) ^ 2 + t162; t158 - pkin(4) * t49 + m(7) * (qJ(5) * t5 - t169 * t206) - t165 * t206 + t195 * qJ(5) + m(6) * (-pkin(4) * t11 - qJ(5) * t10) - t196 * t185; -qJ(5) * t68 + t157 - t164 * t206 + m(6) * (-pkin(4) * t84 + qJ(5) * t86) + m(7) * (qJ(5) * t63 - t166 * t206) + (-pkin(4) * t111 - qJ(5) * t110) * mrSges(6,1); t160 + (t132 - pkin(4)) * mrSges(6,2) + m(7) * (t128 * t171 + t189) + m(6) * (-pkin(4) * t132 + t189) + t162 + t209 * (qJ(5) + t129) + (t128 - t206) * t173; -0.2e1 * pkin(4) * mrSges(6,2) + qJ(5) * t170 - t206 * t163 + m(6) * (pkin(4) ^ 2 + t154) + m(7) * (-t179 * t206 ^ 2 + t154) + t162; m(6) * t11 + m(7) * t169 + t165 + t49; m(6) * t84 + m(7) * t166 + t111 * mrSges(6,1) + t164; m(6) * t132 - t128 * t174 + t168; -m(6) * pkin(4) + m(7) * t171 + t168; m(6) - t174; mrSges(7,1) * t1 - mrSges(7,2) * t2 + t12; mrSges(7,1) * t20 - mrSges(7,2) * t21 + t45; t128 * t167 + t115; -t206 * t194 + t136 + (mrSges(7,2) * t206 - Ifges(7,6)) * t145; t167; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t3(1) t3(2) t3(4) t3(7) t3(11) t3(16); t3(2) t3(3) t3(5) t3(8) t3(12) t3(17); t3(4) t3(5) t3(6) t3(9) t3(13) t3(18); t3(7) t3(8) t3(9) t3(10) t3(14) t3(19); t3(11) t3(12) t3(13) t3(14) t3(15) t3(20); t3(16) t3(17) t3(18) t3(19) t3(20) t3(21);];
Mq  = res;
