% Calculate joint inertia matrix for
% S6RRRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRR13_inertiaJ_slag_vp2(qJ, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_inertiaJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_inertiaJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_inertiaJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR13_inertiaJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR13_inertiaJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:52:04
% EndTime: 2019-03-09 19:52:11
% DurationCPUTime: 2.76s
% Computational Cost: add. (8148->518), mult. (20928->751), div. (0->0), fcn. (24144->14), ass. (0->194)
t204 = sin(qJ(2));
t197 = sin(pkin(6));
t207 = cos(qJ(2));
t226 = t197 * t207;
t200 = cos(pkin(6));
t240 = pkin(1) * t200;
t153 = pkin(9) * t226 + t204 * t240;
t196 = sin(pkin(7));
t199 = cos(pkin(7));
t224 = t199 * t207;
t219 = t197 * t224;
t114 = (t196 * t200 + t219) * pkin(10) + t153;
t183 = t207 * t240;
t227 = t197 * t204;
t124 = pkin(2) * t200 + t183 + (-pkin(10) * t199 - pkin(9)) * t227;
t135 = (-pkin(10) * t196 * t204 - pkin(2) * t207 - pkin(1)) * t197;
t203 = sin(qJ(3));
t206 = cos(qJ(3));
t58 = -t203 * t114 + (t124 * t199 + t135 * t196) * t206;
t198 = cos(pkin(13));
t236 = pkin(11) + qJ(4);
t165 = t236 * t198;
t202 = sin(qJ(5));
t195 = sin(pkin(13));
t215 = t236 * t195;
t241 = cos(qJ(5));
t126 = t165 * t202 + t241 * t215;
t260 = t126 ^ 2;
t259 = 0.2e1 * t126;
t228 = t196 * t206;
t118 = -t200 * t228 + t203 * t227 - t206 * t219;
t229 = t196 * t203;
t119 = t200 * t229 + (t203 * t224 + t204 * t206) * t197;
t147 = -t196 * t226 + t199 * t200;
t81 = -t119 * t195 + t147 * t198;
t82 = t119 * t198 + t147 * t195;
t54 = t202 * t82 - t241 * t81;
t55 = t202 * t81 + t241 * t82;
t22 = Ifges(6,1) * t55 - Ifges(6,4) * t54 + Ifges(6,5) * t118;
t258 = t22 / 0.2e1;
t39 = Ifges(5,4) * t82 + Ifges(5,2) * t81 + Ifges(5,6) * t118;
t257 = t39 / 0.2e1;
t40 = Ifges(5,1) * t82 + Ifges(5,4) * t81 + Ifges(5,5) * t118;
t256 = t40 / 0.2e1;
t146 = -t195 * t229 + t198 * t199;
t148 = t195 * t199 + t198 * t229;
t103 = -t241 * t146 + t148 * t202;
t104 = t202 * t146 + t241 * t148;
t64 = Ifges(6,1) * t104 - Ifges(6,4) * t103 - Ifges(6,5) * t228;
t255 = t64 / 0.2e1;
t159 = t195 * t202 - t241 * t198;
t160 = t241 * t195 + t202 * t198;
t201 = sin(qJ(6));
t205 = cos(qJ(6));
t234 = Ifges(7,4) * t205;
t90 = Ifges(7,6) * t159 + (-Ifges(7,2) * t201 + t234) * t160;
t254 = t90 / 0.2e1;
t235 = Ifges(7,4) * t201;
t91 = Ifges(7,5) * t159 + (Ifges(7,1) * t205 - t235) * t160;
t253 = t91 / 0.2e1;
t97 = Ifges(5,4) * t148 + Ifges(5,2) * t146 - Ifges(5,6) * t228;
t252 = t97 / 0.2e1;
t98 = Ifges(5,1) * t148 + Ifges(5,4) * t146 - Ifges(5,5) * t228;
t251 = t98 / 0.2e1;
t123 = Ifges(6,1) * t160 - Ifges(6,4) * t159;
t250 = t123 / 0.2e1;
t167 = Ifges(5,4) * t195 + Ifges(5,2) * t198;
t249 = t167 / 0.2e1;
t168 = Ifges(5,1) * t195 + Ifges(5,4) * t198;
t248 = t168 / 0.2e1;
t170 = Ifges(7,5) * t201 + Ifges(7,6) * t205;
t247 = t170 / 0.2e1;
t171 = Ifges(7,2) * t205 + t235;
t246 = t171 / 0.2e1;
t172 = Ifges(7,1) * t201 + t234;
t245 = t172 / 0.2e1;
t244 = -t201 / 0.2e1;
t243 = t201 / 0.2e1;
t242 = t205 / 0.2e1;
t239 = pkin(2) * t206;
t238 = pkin(12) * t201;
t237 = pkin(12) * t205;
t78 = -t124 * t196 + t199 * t135;
t45 = pkin(3) * t118 - qJ(4) * t119 + t78;
t225 = t199 * t203;
t59 = t206 * t114 + t124 * t225 + t135 * t229;
t48 = qJ(4) * t147 + t59;
t23 = -t195 * t48 + t198 * t45;
t12 = pkin(4) * t118 - pkin(11) * t82 + t23;
t24 = t195 * t45 + t198 * t48;
t15 = pkin(11) * t81 + t24;
t6 = t202 * t12 + t241 * t15;
t152 = pkin(2) * t225 + pkin(10) * t228;
t139 = qJ(4) * t199 + t152;
t140 = (-pkin(3) * t206 - qJ(4) * t203 - pkin(2)) * t196;
t94 = -t139 * t195 + t198 * t140;
t72 = -pkin(4) * t228 - pkin(11) * t148 + t94;
t95 = t198 * t139 + t195 * t140;
t76 = pkin(11) * t146 + t95;
t42 = t202 * t72 + t241 * t76;
t151 = -pkin(9) * t227 + t183;
t233 = t151 * mrSges(3,1);
t232 = t153 * mrSges(3,2);
t231 = t160 * t201;
t230 = t160 * t205;
t121 = Ifges(6,5) * t160 - Ifges(6,6) * t159;
t223 = t195 ^ 2 + t198 ^ 2;
t222 = t201 ^ 2 + t205 ^ 2;
t29 = t118 * t205 - t201 * t55;
t30 = t118 * t201 + t205 * t55;
t7 = Ifges(7,5) * t30 + Ifges(7,6) * t29 + Ifges(7,3) * t54;
t21 = Ifges(6,4) * t55 - Ifges(6,2) * t54 + Ifges(6,6) * t118;
t221 = t7 / 0.2e1 - t21 / 0.2e1;
t20 = Ifges(6,5) * t55 - Ifges(6,6) * t54 + Ifges(6,3) * t118;
t85 = -t104 * t201 - t205 * t228;
t86 = t104 * t205 - t201 * t228;
t33 = Ifges(7,5) * t86 + Ifges(7,6) * t85 + Ifges(7,3) * t103;
t63 = Ifges(6,4) * t104 - Ifges(6,2) * t103 - Ifges(6,6) * t228;
t220 = t33 / 0.2e1 - t63 / 0.2e1;
t122 = Ifges(6,4) * t160 - Ifges(6,2) * t159;
t89 = Ifges(7,5) * t230 - Ifges(7,6) * t231 + Ifges(7,3) * t159;
t218 = t89 / 0.2e1 - t122 / 0.2e1;
t68 = Ifges(4,5) * t119 - Ifges(4,6) * t118 + Ifges(4,3) * t147;
t136 = Ifges(4,5) * t229 + Ifges(4,6) * t228 + Ifges(4,3) * t199;
t217 = Ifges(3,5) * t227 + Ifges(3,6) * t226 + Ifges(3,3) * t200;
t185 = -pkin(4) * t198 - pkin(3);
t216 = t121 / 0.2e1 + Ifges(5,5) * t195 / 0.2e1 + Ifges(5,6) * t198 / 0.2e1;
t56 = -t81 * mrSges(5,1) + t82 * mrSges(5,2);
t25 = t54 * mrSges(6,1) + t55 * mrSges(6,2);
t67 = t103 * mrSges(6,1) + t104 * mrSges(6,2);
t106 = -t146 * mrSges(5,1) + t148 * mrSges(5,2);
t164 = -t198 * mrSges(5,1) + t195 * mrSges(5,2);
t120 = t159 * mrSges(6,1) + t160 * mrSges(6,2);
t53 = -pkin(3) * t147 - t58;
t26 = -pkin(4) * t81 + t53;
t10 = pkin(5) * t54 - pkin(12) * t55 + t26;
t4 = pkin(12) * t118 + t6;
t1 = t10 * t205 - t201 * t4;
t2 = t10 * t201 + t205 * t4;
t214 = -t1 * t201 + t2 * t205;
t213 = mrSges(7,1) * t201 + mrSges(7,2) * t205;
t32 = -pkin(12) * t228 + t42;
t178 = pkin(10) * t229;
t142 = t178 + (-pkin(3) - t239) * t199;
t105 = -pkin(4) * t146 + t142;
t47 = pkin(5) * t103 - pkin(12) * t104 + t105;
t16 = -t201 * t32 + t205 * t47;
t17 = t201 * t47 + t205 * t32;
t212 = -t16 * t201 + t17 * t205;
t62 = Ifges(6,5) * t104 - Ifges(6,6) * t103 - Ifges(6,3) * t228;
t5 = t241 * t12 - t202 * t15;
t41 = -t202 * t76 + t241 * t72;
t169 = -mrSges(7,1) * t205 + mrSges(7,2) * t201;
t162 = -mrSges(4,2) * t199 + mrSges(4,3) * t228;
t161 = mrSges(4,1) * t199 - mrSges(4,3) * t229;
t150 = t199 * t239 - t178;
t149 = (-mrSges(4,1) * t206 + mrSges(4,2) * t203) * t196;
t138 = Ifges(4,5) * t199 + (Ifges(4,1) * t203 + Ifges(4,4) * t206) * t196;
t137 = Ifges(4,6) * t199 + (Ifges(4,4) * t203 + Ifges(4,2) * t206) * t196;
t130 = -mrSges(5,1) * t228 - mrSges(5,3) * t148;
t129 = mrSges(5,2) * t228 + mrSges(5,3) * t146;
t128 = t241 * t165 - t202 * t215;
t113 = mrSges(7,1) * t159 - mrSges(7,3) * t230;
t112 = -mrSges(7,2) * t159 - mrSges(7,3) * t231;
t111 = pkin(5) * t159 - pkin(12) * t160 + t185;
t107 = t213 * t160;
t96 = Ifges(5,5) * t148 + Ifges(5,6) * t146 - Ifges(5,3) * t228;
t93 = -mrSges(6,1) * t228 - mrSges(6,3) * t104;
t92 = mrSges(6,2) * t228 - mrSges(6,3) * t103;
t88 = mrSges(4,1) * t147 - mrSges(4,3) * t119;
t87 = -mrSges(4,2) * t147 - mrSges(4,3) * t118;
t77 = mrSges(4,1) * t118 + mrSges(4,2) * t119;
t74 = t111 * t201 + t128 * t205;
t73 = t111 * t205 - t128 * t201;
t70 = Ifges(4,1) * t119 - Ifges(4,4) * t118 + Ifges(4,5) * t147;
t69 = Ifges(4,4) * t119 - Ifges(4,2) * t118 + Ifges(4,6) * t147;
t66 = mrSges(5,1) * t118 - mrSges(5,3) * t82;
t65 = -mrSges(5,2) * t118 + mrSges(5,3) * t81;
t61 = mrSges(7,1) * t103 - mrSges(7,3) * t86;
t60 = -mrSges(7,2) * t103 + mrSges(7,3) * t85;
t57 = -mrSges(7,1) * t85 + mrSges(7,2) * t86;
t38 = Ifges(5,5) * t82 + Ifges(5,6) * t81 + Ifges(5,3) * t118;
t37 = mrSges(6,1) * t118 - mrSges(6,3) * t55;
t36 = -mrSges(6,2) * t118 - mrSges(6,3) * t54;
t35 = Ifges(7,1) * t86 + Ifges(7,4) * t85 + Ifges(7,5) * t103;
t34 = Ifges(7,4) * t86 + Ifges(7,2) * t85 + Ifges(7,6) * t103;
t31 = pkin(5) * t228 - t41;
t19 = mrSges(7,1) * t54 - mrSges(7,3) * t30;
t18 = -mrSges(7,2) * t54 + mrSges(7,3) * t29;
t13 = -mrSges(7,1) * t29 + mrSges(7,2) * t30;
t9 = Ifges(7,1) * t30 + Ifges(7,4) * t29 + Ifges(7,5) * t54;
t8 = Ifges(7,4) * t30 + Ifges(7,2) * t29 + Ifges(7,6) * t54;
t3 = -t118 * pkin(5) - t5;
t11 = [(t20 + t38 - t69) * t118 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) + m(6) * (t26 ^ 2 + t5 ^ 2 + t6 ^ 2) + m(5) * (t23 ^ 2 + t24 ^ 2 + t53 ^ 2) + m(4) * (t58 ^ 2 + t59 ^ 2 + t78 ^ 2) + m(3) * (t151 ^ 2 + t153 ^ 2) + ((Ifges(3,5) * t204 + Ifges(3,6) * t207) * t200 + 0.2e1 * (-t151 * t204 + t153 * t207) * mrSges(3,3) + (m(3) * pkin(1) ^ 2 - 0.2e1 * pkin(1) * (-mrSges(3,1) * t207 + mrSges(3,2) * t204) + t204 * (Ifges(3,1) * t204 + Ifges(3,4) * t207) + t207 * (Ifges(3,4) * t204 + Ifges(3,2) * t207)) * t197) * t197 + (t7 - t21) * t54 + (t217 - 0.2e1 * t232 + 0.2e1 * t233) * t200 + Ifges(2,3) + 0.2e1 * t3 * t13 + 0.2e1 * t2 * t18 + 0.2e1 * t1 * t19 + 0.2e1 * t26 * t25 + t29 * t8 + t30 * t9 + 0.2e1 * t6 * t36 + 0.2e1 * t5 * t37 + t55 * t22 + 0.2e1 * t53 * t56 + 0.2e1 * t24 * t65 + 0.2e1 * t23 * t66 + 0.2e1 * t78 * t77 + t81 * t39 + t82 * t40 + 0.2e1 * t59 * t87 + 0.2e1 * t58 * t88 + t119 * t70 + t147 * t68; m(7) * (t1 * t16 + t17 * t2 + t3 * t31) + m(6) * (t105 * t26 + t41 * t5 + t42 * t6) + m(5) * (t142 * t53 + t23 * t94 + t24 * t95) + m(4) * (-pkin(2) * t196 * t78 + t150 * t58 + t152 * t59) + (t62 / 0.2e1 + t96 / 0.2e1 - t137 / 0.2e1) * t118 - t232 + t233 + t217 + (-pkin(2) * t77 + t203 * t70 / 0.2e1 + (-t20 / 0.2e1 - t38 / 0.2e1 + t69 / 0.2e1) * t206) * t196 + t82 * t251 + t81 * t252 + t55 * t255 + t148 * t256 + t146 * t257 + t104 * t258 + t221 * t103 + t220 * t54 + t17 * t18 + t16 * t19 + t31 * t13 + t29 * t34 / 0.2e1 + t30 * t35 / 0.2e1 + t41 * t37 + t42 * t36 + t3 * t57 + t2 * t60 + t1 * t61 + t26 * t67 + t85 * t8 / 0.2e1 + t86 * t9 / 0.2e1 + t6 * t92 + t5 * t93 + t94 * t66 + t95 * t65 + t105 * t25 + t53 * t106 + t24 * t129 + t23 * t130 + t119 * t138 / 0.2e1 + t142 * t56 + t147 * t136 / 0.2e1 + t78 * t149 + t150 * t88 + t152 * t87 + t58 * t161 + t59 * t162 + t199 * t68 / 0.2e1; Ifges(3,3) + 0.2e1 * t31 * t57 + 0.2e1 * t17 * t60 + 0.2e1 * t16 * t61 + t85 * t34 + t86 * t35 + 0.2e1 * t42 * t92 + 0.2e1 * t41 * t93 + t104 * t64 + 0.2e1 * t105 * t67 + 0.2e1 * t95 * t129 + 0.2e1 * t94 * t130 + 0.2e1 * t142 * t106 + t146 * t97 + t148 * t98 + 0.2e1 * t150 * t161 + 0.2e1 * t152 * t162 + t199 * t136 + (t33 - t63) * t103 + (-0.2e1 * pkin(2) * t149 + t203 * t138 + (t137 - t62 - t96) * t206) * t196 + m(7) * (t16 ^ 2 + t17 ^ 2 + t31 ^ 2) + m(6) * (t105 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t142 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(4) * (pkin(2) ^ 2 * t196 ^ 2 + t150 ^ 2 + t152 ^ 2); t68 + (t13 - t37) * t126 + m(5) * (-pkin(3) * t53 + (-t195 * t23 + t198 * t24) * qJ(4)) + (-t23 * mrSges(5,3) - qJ(4) * t66 + t256) * t195 + (t24 * mrSges(5,3) + qJ(4) * t65 + t257) * t198 + (-t5 * mrSges(6,3) + t9 * t242 + t8 * t244 + t258) * t160 + t82 * t248 + t81 * t249 + t55 * t250 + t30 * t253 + t29 * t254 + (-t6 * mrSges(6,3) + t221) * t159 + t216 * t118 + t218 * t54 - pkin(3) * t56 + t58 * mrSges(4,1) - t59 * mrSges(4,2) + m(7) * (t1 * t73 + t126 * t3 + t2 * t74) + m(6) * (-t126 * t5 + t128 * t6 + t185 * t26) + t73 * t19 + t74 * t18 + t3 * t107 + t2 * t112 + t1 * t113 + t26 * t120 + t128 * t36 + t53 * t164 + t185 * t25; m(7) * (t126 * t31 + t16 * t73 + t17 * t74) + m(6) * (t105 * t185 - t126 * t41 + t128 * t42) + (t57 - t93) * t126 + t136 + m(5) * (-pkin(3) * t142 + (-t195 * t94 + t198 * t95) * qJ(4)) + (-t41 * mrSges(6,3) + t35 * t242 + t34 * t244 + t255) * t160 + (-t94 * mrSges(5,3) - qJ(4) * t130 + t251) * t195 + (t95 * mrSges(5,3) + qJ(4) * t129 + t252) * t198 + t148 * t248 + t146 * t249 + t104 * t250 + t86 * t253 + t85 * t254 - t216 * t228 + t218 * t103 + (-t42 * mrSges(6,3) + t220) * t159 + t73 * t61 + t74 * t60 - pkin(3) * t106 + t31 * t107 + t17 * t112 + t16 * t113 + t105 * t120 + t128 * t92 + t150 * mrSges(4,1) - t152 * mrSges(4,2) + t142 * t164 + t185 * t67; -0.2e1 * pkin(3) * t164 + t107 * t259 + 0.2e1 * t74 * t112 + 0.2e1 * t73 * t113 + 0.2e1 * t185 * t120 + t198 * t167 + t195 * t168 + Ifges(4,3) + 0.2e1 * t223 * qJ(4) * mrSges(5,3) + (-0.2e1 * mrSges(6,3) * t128 - t122 + t89) * t159 + m(7) * (t73 ^ 2 + t74 ^ 2 + t260) + m(6) * (t128 ^ 2 + t185 ^ 2 + t260) + m(5) * (qJ(4) ^ 2 * t223 + pkin(3) ^ 2) + (mrSges(6,3) * t259 - t201 * t90 + t205 * t91 + t123) * t160; t201 * t18 + t205 * t19 + m(7) * (t1 * t205 + t2 * t201) + m(6) * t26 + m(5) * t53 + t25 + t56; t201 * t60 + t205 * t61 + m(7) * (t16 * t205 + t17 * t201) + m(6) * t105 + m(5) * t142 + t106 + t67; -m(5) * pkin(3) + t201 * t112 + t205 * t113 + m(7) * (t201 * t74 + t205 * t73) + m(6) * t185 + t120 + t164; m(7) * t222 + m(5) + m(6); t3 * t169 + t9 * t243 + t8 * t242 + t30 * t245 + t54 * t247 + t29 * t246 + m(7) * (-pkin(5) * t3 + pkin(12) * t214) + t18 * t237 - pkin(5) * t13 - t19 * t238 - t6 * mrSges(6,2) + t5 * mrSges(6,1) + t214 * mrSges(7,3) + t20; t34 * t242 + t35 * t243 - pkin(5) * t57 + m(7) * (-pkin(5) * t31 + pkin(12) * t212) + t31 * t169 + t60 * t237 - t61 * t238 + t103 * t247 + t86 * t245 + t85 * t246 + t41 * mrSges(6,1) - t42 * mrSges(6,2) + t212 * mrSges(7,3) + t62; t159 * t247 + t112 * t237 - t113 * t238 + t90 * t242 + t91 * t243 - pkin(5) * t107 - t128 * mrSges(6,2) + (t171 * t244 + t172 * t242) * t160 + t121 + (m(7) * pkin(12) + mrSges(7,3)) * (-t201 * t73 + t205 * t74) + (-m(7) * pkin(5) - mrSges(6,1) + t169) * t126; 0; Ifges(6,3) + t205 * t171 + m(7) * (pkin(12) ^ 2 * t222 + pkin(5) ^ 2) - 0.2e1 * pkin(5) * t169 + t201 * t172 + 0.2e1 * t222 * pkin(12) * mrSges(7,3); mrSges(7,1) * t1 - mrSges(7,2) * t2 + t7; mrSges(7,1) * t16 - mrSges(7,2) * t17 + t33; mrSges(7,1) * t73 - mrSges(7,2) * t74 + t89; -t169; -pkin(12) * t213 + t170; Ifges(7,3);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
