% Calculate joint inertia matrix for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR10V2_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:41
% EndTime: 2019-04-11 14:41:48
% DurationCPUTime: 1.73s
% Computational Cost: add. (2306->374), mult. (4922->538), div. (0->0), fcn. (5371->10), ass. (0->173)
t164 = sin(qJ(3));
t145 = pkin(2) * t164 + pkin(5);
t256 = pkin(5) + t145;
t163 = sin(qJ(4));
t157 = t163 ^ 2;
t168 = cos(qJ(4));
t160 = t168 ^ 2;
t255 = (t157 + t160) * mrSges(5,3);
t167 = cos(qJ(5));
t198 = t163 * t167;
t118 = -mrSges(6,1) * t168 - mrSges(6,3) * t198;
t161 = sin(qJ(6));
t166 = cos(qJ(6));
t103 = -t161 * t198 - t166 * t168;
t104 = -t161 * t168 + t166 * t198;
t67 = -mrSges(7,1) * t103 + mrSges(7,2) * t104;
t254 = t67 - t118;
t169 = cos(qJ(3));
t146 = -pkin(2) * t169 - pkin(3);
t162 = sin(qJ(5));
t199 = t162 * t168;
t84 = t145 * t199 - t167 * t146;
t253 = t84 ^ 2;
t201 = t162 * t163;
t79 = -mrSges(7,2) * t201 + mrSges(7,3) * t103;
t252 = 0.2e1 * t79;
t80 = mrSges(7,1) * t201 - mrSges(7,3) * t104;
t251 = 0.2e1 * t80;
t119 = pkin(3) * t167 + pkin(5) * t199;
t250 = t119 ^ 2;
t116 = mrSges(6,2) * t168 - mrSges(6,3) * t201;
t249 = 0.2e1 * t116;
t170 = cos(qJ(2));
t147 = -pkin(2) * t170 - pkin(1);
t248 = 0.2e1 * t147;
t52 = Ifges(7,1) * t104 + Ifges(7,4) * t103 + Ifges(7,5) * t201;
t247 = t52 / 0.2e1;
t165 = sin(qJ(2));
t113 = t164 * t165 - t169 * t170;
t114 = t164 * t170 + t165 * t169;
t61 = -t167 * t113 + t114 * t199;
t246 = t61 / 0.2e1;
t197 = t167 * t168;
t62 = t113 * t162 + t114 * t197;
t245 = t62 / 0.2e1;
t222 = Ifges(7,4) * t166;
t89 = -Ifges(7,6) * t167 + (-Ifges(7,2) * t161 + t222) * t162;
t244 = t89 / 0.2e1;
t223 = Ifges(7,4) * t161;
t91 = -Ifges(7,5) * t167 + (Ifges(7,1) * t166 - t223) * t162;
t243 = t91 / 0.2e1;
t242 = t103 / 0.2e1;
t241 = t104 / 0.2e1;
t130 = Ifges(7,2) * t166 + t223;
t240 = t130 / 0.2e1;
t133 = Ifges(7,1) * t161 + t222;
t239 = t133 / 0.2e1;
t238 = t161 / 0.2e1;
t237 = t162 / 0.2e1;
t236 = t166 / 0.2e1;
t235 = -t167 / 0.2e1;
t234 = -t168 / 0.2e1;
t233 = pkin(5) * t145;
t232 = pkin(6) * t167;
t231 = pkin(6) * t168;
t207 = t114 * t163;
t28 = -t161 * t62 + t166 * t207;
t29 = t161 * t207 + t166 * t62;
t230 = -mrSges(6,1) * t207 - mrSges(7,1) * t28 + mrSges(7,2) * t29 + mrSges(6,3) * t62;
t206 = t114 * t168;
t229 = -mrSges(5,1) * t113 + mrSges(6,1) * t61 + mrSges(6,2) * t62 + mrSges(5,3) * t206;
t228 = mrSges(6,2) * t167;
t227 = Ifges(5,4) * t163;
t226 = Ifges(5,4) * t168;
t225 = Ifges(6,4) * t162;
t224 = Ifges(6,4) * t167;
t221 = Ifges(5,6) * t113;
t220 = t119 * t84;
t68 = pkin(3) * t113 - pkin(5) * t114 + t147;
t65 = t68 ^ 2;
t219 = t157 * t65;
t218 = t161 * t80;
t217 = t162 * t84;
t216 = t163 * t68;
t86 = t145 * t197 + t162 * t146;
t215 = t167 * t86;
t214 = t168 * mrSges(5,2);
t213 = t168 * t68;
t106 = (mrSges(6,1) * t162 + t228) * t163;
t212 = t68 * t106;
t211 = mrSges(6,1) * t167 - mrSges(6,2) * t162 + mrSges(5,1);
t210 = Ifges(5,5) * t206 + Ifges(5,3) * t113;
t125 = -mrSges(7,1) * t166 + mrSges(7,2) * t161;
t209 = t125 - mrSges(6,1);
t205 = t119 * t162;
t121 = -pkin(3) * t162 + pkin(5) * t197;
t204 = t121 * t167;
t203 = t161 * t162;
t202 = t162 * t118;
t200 = t162 * t166;
t128 = Ifges(7,5) * t161 + Ifges(7,6) * t166;
t129 = Ifges(6,5) * t162 + Ifges(6,6) * t167;
t196 = Ifges(5,5) * t163 + Ifges(5,6) * t168;
t195 = t161 ^ 2 + t166 ^ 2;
t193 = 0.2e1 * pkin(6);
t6 = Ifges(7,5) * t29 + Ifges(7,6) * t28 + Ifges(7,3) * t61;
t19 = Ifges(6,5) * t62 - Ifges(6,6) * t61 + Ifges(6,3) * t207;
t50 = Ifges(7,5) * t104 + Ifges(7,6) * t103 + Ifges(7,3) * t201;
t192 = -t201 / 0.2e1;
t191 = t201 / 0.2e1;
t190 = t198 / 0.2e1;
t189 = Ifges(7,5) * t200 - Ifges(7,3) * t167;
t188 = 0.2e1 * t254;
t187 = m(7) * t195 * pkin(6) ^ 2;
t23 = -pkin(6) * t62 - t213;
t30 = (pkin(6) * t114 + t167 * t68) * t163;
t10 = t161 * t23 + t166 * t30;
t9 = -t161 * t30 + t166 * t23;
t186 = t10 * t166 - t161 * t9;
t185 = mrSges(7,1) * t161 + mrSges(7,2) * t166;
t100 = (t145 - t232) * t163;
t77 = t86 - t231;
t40 = t100 * t166 - t161 * t77;
t41 = t100 * t161 + t166 * t77;
t184 = -t161 * t40 + t166 * t41;
t101 = t121 - t231;
t122 = (pkin(5) - t232) * t163;
t63 = -t101 * t161 + t122 * t166;
t64 = t101 * t166 + t122 * t161;
t183 = -t161 * t63 + t166 * t64;
t182 = -t161 * t79 - t166 * t80;
t181 = (mrSges(4,1) * t169 - mrSges(4,2) * t164) * pkin(2);
t180 = -t211 * t163 - t214;
t88 = Ifges(6,5) * t198 - Ifges(6,6) * t201 - Ifges(6,3) * t168;
t179 = 0.2e1 * t163 * t106 + 0.2e1 * t255;
t132 = Ifges(5,2) * t168 + t227;
t135 = Ifges(5,1) * t163 + t226;
t51 = Ifges(7,4) * t104 + Ifges(7,2) * t103 + Ifges(7,6) * t201;
t92 = -Ifges(6,5) * t168 + (Ifges(6,1) * t167 - t225) * t163;
t178 = t103 * t51 + t104 * t52 + t163 * t135 + t92 * t198 + t50 * t201 + Ifges(4,3) + (t132 - t88) * t168;
t131 = Ifges(6,2) * t167 + t225;
t134 = Ifges(6,1) * t162 + t224;
t87 = -Ifges(7,6) * t203 + t189;
t90 = -Ifges(6,6) * t168 + (-Ifges(6,2) * t162 + t224) * t163;
t177 = t129 * t234 - t51 * t203 / 0.2e1 + t200 * t247 + t50 * t235 + t89 * t242 + t91 * t241 + t87 * t191 + t92 * t237 + t167 * t90 / 0.2e1 + t131 * t192 + t134 * t190 + t196;
t176 = t166 * pkin(6) * t79 + t103 * t240 + t104 * t239 + t128 * t191 + t51 * t236 + t52 * t238 + t88;
t175 = -t90 * t201 + t178;
t20 = Ifges(6,4) * t62 - Ifges(6,2) * t61 + Ifges(6,6) * t207;
t21 = Ifges(6,1) * t62 - Ifges(6,4) * t61 + Ifges(6,5) * t207;
t42 = t221 + (-Ifges(5,2) * t163 + t226) * t114;
t43 = Ifges(5,5) * t113 + (Ifges(5,1) * t168 - t227) * t114;
t7 = Ifges(7,4) * t29 + Ifges(7,2) * t28 + Ifges(7,6) * t61;
t8 = Ifges(7,1) * t29 + Ifges(7,4) * t28 + Ifges(7,5) * t61;
t174 = t10 * t79 + t6 * t191 + t8 * t241 + t7 * t242 + t9 * t80 + Ifges(4,5) * t114 + t28 * t51 / 0.2e1 + t29 * t247 + t20 * t192 + t21 * t190 + t19 * t234 + t50 * t246 - t61 * t90 / 0.2e1 + t92 * t245 + t163 * t43 / 0.2e1 + t168 * t42 / 0.2e1 + t135 * t206 / 0.2e1 + (t116 * t198 + t67 * t201) * t68 + (t88 / 0.2e1 - t132 / 0.2e1) * t207 + (-Ifges(4,6) + t196 / 0.2e1) * t113;
t172 = pkin(5) ^ 2;
t159 = t167 ^ 2;
t156 = t162 ^ 2;
t154 = t157 * t172;
t143 = t145 ^ 2;
t137 = t157 * t143;
t127 = -mrSges(5,1) * t168 + mrSges(5,2) * t163;
t124 = t157 * t233;
t117 = -mrSges(7,1) * t167 - mrSges(7,3) * t200;
t115 = mrSges(7,2) * t167 - mrSges(7,3) * t203;
t105 = t185 * t162;
t69 = -mrSges(5,2) * t113 - mrSges(5,3) * t207;
t66 = (mrSges(5,1) * t163 + t214) * t114;
t54 = t160 * t65;
t53 = t156 * t219;
t35 = -mrSges(6,2) * t207 - mrSges(6,3) * t61;
t17 = mrSges(7,1) * t61 - mrSges(7,3) * t29;
t16 = -mrSges(7,2) * t61 + mrSges(7,3) * t28;
t1 = [-0.2e1 * pkin(1) * (-mrSges(3,1) * t170 + mrSges(3,2) * t165) + t165 * (Ifges(3,1) * t165 + Ifges(3,4) * t170) + t170 * (Ifges(3,4) * t165 + Ifges(3,2) * t170) + m(4) * t147 ^ 2 + t62 * t21 + t28 * t7 + t29 * t8 + 0.2e1 * t10 * t16 + 0.2e1 * t9 * t17 + Ifges(2,3) + m(3) * pkin(1) ^ 2 + (mrSges(4,2) * t248 + Ifges(4,1) * t114) * t114 + (mrSges(4,1) * t248 - 0.2e1 * Ifges(4,4) * t114 + Ifges(4,2) * t113 + t210) * t113 + (t6 - t20) * t61 + (t114 * t43 - 0.2e1 * t229 * t68) * t168 + m(5) * (t54 + t219) + m(6) * (t159 * t219 + t53 + t54) + m(7) * (t10 ^ 2 + t9 ^ 2 + t53) + ((t19 - t42 - t221) * t114 + 0.2e1 * (t162 * t230 + t167 * t35 + t69) * t68) * t163; t174 + m(7) * (t10 * t41 + t40 * t9) + (t229 * t145 + (-t202 + m(6) * (-t145 * t168 + t215 + t217) + m(7) * t217) * t68) * t163 + (-t113 * t164 - t114 * t169) * pkin(2) * mrSges(4,3) + Ifges(3,6) * t170 + Ifges(3,5) * t165 + t146 * t66 + t86 * t35 + t40 * t17 + t41 * t16 + (t145 * t69 - t212) * t168 + t230 * t84; t84 * t188 + 0.2e1 * t181 + m(4) * (t164 ^ 2 + t169 ^ 2) * pkin(2) ^ 2 + t175 + t179 * t145 + 0.2e1 * t146 * t127 + t86 * t249 + t41 * t252 + t40 * t251 + Ifges(3,3) + m(5) * (t143 * t160 + t146 ^ 2 + t137) + m(6) * (t86 ^ 2 + t137 + t253) + m(7) * (t40 ^ 2 + t41 ^ 2 + t253); t174 + m(7) * (t10 * t64 + t63 * t9) + (t229 * pkin(5) + (-t202 + m(6) * (-pkin(5) * t168 + t204 + t205) + m(7) * t205) * t68) * t163 + t121 * t35 + t63 * t17 + t64 * t16 - pkin(3) * t66 + t230 * t119 + (pkin(5) * t69 - t212) * t168; m(5) * (-pkin(3) * t146 + t160 * t233 + t124) + m(7) * (t40 * t63 + t41 * t64 + t220) + m(6) * (t121 * t86 + t124 + t220) + (t63 + t40) * t80 + (t64 + t41) * t79 + (t146 - pkin(3)) * t127 + (t121 + t86) * t116 + t181 + t178 + (t256 * t106 - t162 * t90) * t163 - t254 * (-t84 - t119) + t256 * t255; t119 * t188 + t175 + m(7) * (t63 ^ 2 + t64 ^ 2 + t250) + m(6) * (t121 ^ 2 + t154 + t250) + m(5) * (pkin(3) ^ 2 + t160 * t172 + t154) + t179 * pkin(5) + t121 * t249 - 0.2e1 * pkin(3) * t127 + t64 * t252 + t63 * t251; t9 * t117 + t134 * t245 + t10 * t115 + t28 * t244 + t29 * t243 + (-t131 / 0.2e1 + t87 / 0.2e1) * t61 + t211 * t213 + (t20 / 0.2e1 - t6 / 0.2e1) * t167 + ((t129 / 0.2e1 - Ifges(5,6)) * t114 + (-mrSges(5,2) + (t156 + t159) * mrSges(6,3)) * t68) * t163 + (t21 / 0.2e1 + t8 * t236 - t161 * t7 / 0.2e1 + t105 * t216 + (m(7) * (-t10 * t161 - t166 * t9) - t166 * t17 - t161 * t16) * pkin(6)) * t162 + t210; t177 + (t84 * mrSges(6,3) + (m(7) * (-t161 * t41 - t166 * t40) + t182) * pkin(6)) * t162 + mrSges(6,3) * t215 + t40 * t117 + t41 * t115 + t84 * t105 + t180 * t145; t177 + t180 * pkin(5) + (t119 * mrSges(6,3) + (m(7) * (-t161 * t64 - t166 * t63) + t182) * pkin(6)) * t162 + mrSges(6,3) * t204 + t63 * t117 + t119 * t105 + t64 * t115; Ifges(5,3) + (-t87 + t131) * t167 + t156 * t187 + (-t161 * t89 + t166 * t91 + t134 + (-t115 * t161 - t117 * t166) * t193) * t162; t128 * t246 + t29 * t239 + t28 * t240 + t8 * t238 + t7 * t236 + t186 * mrSges(7,3) + (t209 * t162 - t228) * t216 + (m(7) * t186 + t166 * t16 - t161 * t17) * pkin(6) + t19; -t86 * mrSges(6,2) + t209 * t84 + t184 * mrSges(7,3) + (m(7) * t184 - t218) * pkin(6) + t176; -t121 * mrSges(6,2) + t209 * t119 + t183 * mrSges(7,3) + (m(7) * t183 - t218) * pkin(6) + t176; t128 * t235 + (pkin(6) * t115 + t133 * t237 + t244) * t166 + (-pkin(6) * t117 + t243 - t162 * t130 / 0.2e1) * t161 + t129; t195 * mrSges(7,3) * t193 + t166 * t130 + t161 * t133 + Ifges(6,3) + t187; mrSges(7,1) * t9 - mrSges(7,2) * t10 + t6; mrSges(7,1) * t40 - mrSges(7,2) * t41 + t50; mrSges(7,1) * t63 - mrSges(7,2) * t64 + t50; (-Ifges(7,6) * t161 + t125 * pkin(6)) * t162 + t189; -t185 * pkin(6) + t128; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
