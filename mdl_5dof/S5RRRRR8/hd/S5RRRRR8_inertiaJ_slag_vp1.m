% Calculate joint inertia matrix for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR8_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR8_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:24
% EndTime: 2019-12-31 22:24:31
% DurationCPUTime: 2.78s
% Computational Cost: add. (9350->388), mult. (9693->553), div. (0->0), fcn. (10537->10), ass. (0->205)
t200 = cos(qJ(4));
t184 = pkin(4) * t200 + pkin(3);
t203 = -pkin(9) - pkin(8);
t197 = sin(qJ(4));
t199 = sin(qJ(1));
t248 = t197 * t199;
t196 = qJ(2) + qJ(3);
t189 = cos(t196);
t202 = cos(qJ(1));
t249 = t189 * t202;
t187 = sin(t196);
t251 = t187 * t202;
t195 = qJ(4) + qJ(5);
t186 = sin(t195);
t188 = cos(t195);
t148 = -t186 * t249 + t188 * t199;
t149 = t186 * t199 + t188 * t249;
t93 = t149 * rSges(6,1) + t148 * rSges(6,2) + rSges(6,3) * t251;
t278 = pkin(4) * t248 + t184 * t249 - t203 * t251 + t93;
t239 = pkin(3) * t249 + pkin(8) * t251;
t277 = -t239 + t278;
t216 = Icges(4,5) * t189 - Icges(4,6) * t187;
t136 = -Icges(4,3) * t202 + t216 * t199;
t137 = Icges(4,3) * t199 + t216 * t202;
t250 = t189 * t199;
t146 = -t186 * t250 - t188 * t202;
t147 = -t186 * t202 + t188 * t250;
t252 = t187 * t199;
t86 = Icges(6,5) * t147 + Icges(6,6) * t146 + Icges(6,3) * t252;
t88 = Icges(6,4) * t147 + Icges(6,2) * t146 + Icges(6,6) * t252;
t90 = Icges(6,1) * t147 + Icges(6,4) * t146 + Icges(6,5) * t252;
t27 = t146 * t88 + t147 * t90 + t86 * t252;
t87 = Icges(6,5) * t149 + Icges(6,6) * t148 + Icges(6,3) * t251;
t89 = Icges(6,4) * t149 + Icges(6,2) * t148 + Icges(6,6) * t251;
t91 = Icges(6,1) * t149 + Icges(6,4) * t148 + Icges(6,5) * t251;
t28 = t146 * t89 + t147 * t91 + t87 * t252;
t15 = t199 * t28 - t202 * t27;
t194 = t202 ^ 2;
t245 = t200 * t202;
t158 = -t189 * t248 - t245;
t246 = t199 * t200;
t247 = t197 * t202;
t159 = t189 * t246 - t247;
t103 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t252;
t105 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t252;
t107 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t252;
t41 = t103 * t252 + t105 * t158 + t107 * t159;
t160 = -t189 * t247 + t246;
t161 = t189 * t245 + t248;
t104 = Icges(5,5) * t161 + Icges(5,6) * t160 + Icges(5,3) * t251;
t106 = Icges(5,4) * t161 + Icges(5,2) * t160 + Icges(5,6) * t251;
t108 = Icges(5,1) * t161 + Icges(5,4) * t160 + Icges(5,5) * t251;
t42 = t104 * t252 + t106 * t158 + t108 * t159;
t21 = t199 * t42 - t202 * t41;
t255 = Icges(4,4) * t189;
t218 = -Icges(4,2) * t187 + t255;
t139 = Icges(4,6) * t199 + t218 * t202;
t256 = Icges(4,4) * t187;
t220 = Icges(4,1) * t189 - t256;
t141 = Icges(4,5) * t199 + t220 * t202;
t214 = -t139 * t187 + t141 * t189;
t138 = -Icges(4,6) * t202 + t218 * t199;
t140 = -Icges(4,5) * t202 + t220 * t199;
t215 = t138 * t187 - t140 * t189;
t276 = -t15 - t21 - t194 * t136 - (t214 * t199 + (-t137 + t215) * t202) * t199;
t193 = t199 ^ 2;
t120 = -Icges(6,3) * t189 + (Icges(6,5) * t188 - Icges(6,6) * t186) * t187;
t121 = -Icges(6,6) * t189 + (Icges(6,4) * t188 - Icges(6,2) * t186) * t187;
t122 = -Icges(6,5) * t189 + (Icges(6,1) * t188 - Icges(6,4) * t186) * t187;
t55 = t120 * t252 + t121 * t146 + t122 * t147;
t5 = -t55 * t189 + (t199 * t27 + t202 * t28) * t187;
t29 = t148 * t88 + t149 * t90 + t86 * t251;
t30 = t148 * t89 + t149 * t91 + t87 * t251;
t56 = t120 * t251 + t121 * t148 + t122 * t149;
t6 = -t56 * t189 + (t199 * t29 + t202 * t30) * t187;
t275 = t6 * t251 + t5 * t252;
t274 = -t189 / 0.2e1;
t273 = t199 / 0.2e1;
t272 = -t202 / 0.2e1;
t198 = sin(qJ(2));
t271 = pkin(2) * t198;
t270 = pkin(3) * t189;
t269 = -pkin(3) + t184;
t268 = pkin(8) + t203;
t201 = cos(qJ(2));
t267 = rSges(3,1) * t201;
t266 = rSges(3,2) * t198;
t265 = t202 * rSges(3,3);
t36 = -t189 * t86 + (-t186 * t88 + t188 * t90) * t187;
t264 = t36 * t202;
t37 = -t189 * t87 + (-t186 * t89 + t188 * t91) * t187;
t263 = t37 * t199;
t48 = -t189 * t103 + (-t105 * t197 + t107 * t200) * t187;
t262 = t48 * t202;
t49 = -t189 * t104 + (-t106 * t197 + t108 * t200) * t187;
t261 = t49 * t199;
t114 = t187 * t188 * t122;
t254 = t121 * t186;
t62 = -t189 * t120 - t187 * t254 + t114;
t260 = t62 * t189;
t123 = -t189 * rSges(6,3) + (rSges(6,1) * t188 - rSges(6,2) * t186) * t187;
t222 = -t147 * rSges(6,1) - t146 * rSges(6,2);
t92 = rSges(6,3) * t252 - t222;
t68 = t123 * t252 + t189 * t92;
t258 = Icges(3,4) * t198;
t257 = Icges(3,4) * t201;
t131 = -Icges(5,6) * t189 + (Icges(5,4) * t200 - Icges(5,2) * t197) * t187;
t253 = t131 * t197;
t204 = -pkin(7) - pkin(6);
t244 = t202 * t204;
t119 = t269 * t187 + t268 * t189;
t243 = -t119 - t123;
t185 = pkin(2) * t201 + pkin(1);
t176 = t202 * t185;
t192 = t202 * pkin(6);
t242 = t199 * (t244 + t192 + (-pkin(1) + t185) * t199) + t202 * (-t202 * pkin(1) + t176 + (-pkin(6) - t204) * t199);
t209 = rSges(4,1) * t249 - rSges(4,2) * t251 + t199 * rSges(4,3);
t224 = rSges(4,1) * t189 - rSges(4,2) * t187;
t94 = t199 * (-t202 * rSges(4,3) + t224 * t199) + t202 * t209;
t133 = -t189 * rSges(5,3) + (rSges(5,1) * t200 - rSges(5,2) * t197) * t187;
t168 = t187 * pkin(3) - t189 * pkin(8);
t241 = -t133 - t168;
t240 = t193 * (pkin(8) * t187 + t270) + t202 * t239;
t238 = t199 * rSges(3,3) + t202 * t267;
t237 = t193 + t194;
t16 = t199 * t30 - t202 * t29;
t43 = t103 * t251 + t105 * t160 + t107 * t161;
t44 = t104 * t251 + t106 * t160 + t108 * t161;
t22 = t199 * t44 - t202 * t43;
t236 = (t193 * t137 + t16 + t22 + (t215 * t202 + (-t136 + t214) * t199) * t202) * t199;
t235 = -t168 + t243;
t110 = t161 * rSges(5,1) + t160 * rSges(5,2) + rSges(5,3) * t251;
t234 = t252 / 0.2e1;
t233 = t251 / 0.2e1;
t167 = rSges(4,1) * t187 + rSges(4,2) * t189;
t232 = -t167 - t271;
t231 = -t168 - t271;
t230 = (t36 + t55) * t234 + (t37 + t56) * t233;
t229 = -t199 * t204 + t176;
t9 = -t260 + (t199 * t36 + t202 * t37) * t187;
t228 = -t189 * t9 + t275;
t223 = -t159 * rSges(5,1) - t158 * rSges(5,2);
t109 = rSges(5,3) * t252 - t223;
t52 = t199 * t109 + t202 * t110 + t240;
t227 = t15 * t234 + t16 * t233 + t5 * t272 + t6 * t273 + (t263 - t264) * t274;
t226 = -t133 + t231;
t225 = -t266 + t267;
t221 = Icges(3,1) * t201 - t258;
t219 = -Icges(3,2) * t198 + t257;
t217 = Icges(3,5) * t201 - Icges(3,6) * t198;
t164 = Icges(4,2) * t189 + t256;
t165 = Icges(4,1) * t187 + t255;
t211 = -t164 * t187 + t165 * t189;
t210 = t231 + t243;
t101 = -pkin(4) * t247 + (-t268 * t187 + t269 * t189) * t199;
t25 = t240 + t277 * t202 + (t101 + t92) * t199;
t207 = t276 * t202 + t236;
t130 = -Icges(5,3) * t189 + (Icges(5,5) * t200 - Icges(5,6) * t197) * t187;
t132 = -Icges(5,5) * t189 + (Icges(5,1) * t200 - Icges(5,4) * t197) * t187;
t60 = t130 * t252 + t131 * t158 + t132 * t159;
t10 = -t60 * t189 + (t199 * t41 + t202 * t42) * t187;
t61 = t130 * t251 + t131 * t160 + t132 * t161;
t11 = -t61 * t189 + (t199 * t43 + t202 * t44) * t187;
t206 = t10 * t272 + t11 * t273 + t21 * t234 + t22 * t233 + t227 + (t261 - t262) * t274;
t163 = Icges(4,5) * t187 + Icges(4,6) * t189;
t205 = -t264 / 0.2e1 + t263 / 0.2e1 - t262 / 0.2e1 + t261 / 0.2e1 + (t139 * t189 + t141 * t187 + t199 * t163 + t211 * t202 + t56 + t61) * t273 + (t138 * t189 + t140 * t187 - t202 * t163 + t211 * t199 + t55 + t60) * t272;
t175 = rSges(2,1) * t202 - rSges(2,2) * t199;
t174 = -rSges(2,1) * t199 - rSges(2,2) * t202;
t173 = rSges(3,1) * t198 + rSges(3,2) * t201;
t151 = Icges(3,3) * t199 + t217 * t202;
t150 = -Icges(3,3) * t202 + t217 * t199;
t135 = t232 * t202;
t134 = t232 * t199;
t125 = t199 * pkin(6) + (pkin(1) - t266) * t202 + t238;
t124 = t265 + t192 + (-pkin(1) - t225) * t199;
t118 = t187 * t200 * t132;
t116 = t209 + t229;
t115 = (rSges(4,3) - t204) * t202 + (-t185 - t224) * t199;
t113 = t241 * t202;
t112 = t241 * t199;
t111 = t202 * (-t202 * t266 + t238) + (t225 * t199 - t265) * t199;
t98 = t226 * t202;
t97 = t226 * t199;
t80 = t92 * t251;
t77 = t229 + t110 + t239;
t76 = -t244 + (-t270 - t185 + (-rSges(5,3) - pkin(8)) * t187) * t199 + t223;
t75 = t235 * t202;
t74 = t235 * t199;
t73 = -t110 * t189 - t133 * t251;
t72 = t109 * t189 + t133 * t252;
t71 = t210 * t202;
t70 = t210 * t199;
t69 = -t123 * t251 - t189 * t93;
t67 = -t189 * t130 - t187 * t253 + t118;
t66 = t229 + t278;
t65 = (pkin(4) * t197 - t204) * t202 + (-t184 * t189 - t185 + (-rSges(6,3) + t203) * t187) * t199 + t222;
t64 = t94 + t242;
t63 = (t109 * t202 - t110 * t199) * t187;
t57 = -t93 * t252 + t80;
t39 = -t189 * t277 + t243 * t251;
t38 = t101 * t189 + t119 * t252 + t68;
t33 = t52 + t242;
t26 = t80 + (t101 * t202 - t199 * t277) * t187;
t24 = t25 + t242;
t1 = [t201 * (Icges(3,2) * t201 + t258) + t198 * (Icges(3,1) * t198 + t257) + Icges(2,3) + t114 + t118 + (-t120 - t130 + t164) * t189 + (t165 - t253 - t254) * t187 + m(6) * (t65 ^ 2 + t66 ^ 2) + m(5) * (t76 ^ 2 + t77 ^ 2) + m(4) * (t115 ^ 2 + t116 ^ 2) + m(3) * (t124 ^ 2 + t125 ^ 2) + m(2) * (t174 ^ 2 + t175 ^ 2); (t194 / 0.2e1 + t193 / 0.2e1) * (Icges(3,5) * t198 + Icges(3,6) * t201) + m(3) * (-t124 * t202 - t125 * t199) * t173 + m(6) * (t65 * t71 + t66 * t70) + m(5) * (t76 * t98 + t77 * t97) + m(4) * (t115 * t135 + t116 * t134) + t205 + (t201 * (Icges(3,6) * t199 + t219 * t202) + t198 * (Icges(3,5) * t199 + t221 * t202)) * t273 + (t201 * (-Icges(3,6) * t202 + t219 * t199) + t198 * (-Icges(3,5) * t202 + t221 * t199)) * t272; m(6) * (t24 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(5) * (t33 ^ 2 + t97 ^ 2 + t98 ^ 2) + m(4) * (t134 ^ 2 + t135 ^ 2 + t64 ^ 2) + t199 * t193 * t151 + m(3) * (t237 * t173 ^ 2 + t111 ^ 2) + t236 + (-t194 * t150 + (-t199 * t150 + t202 * t151) * t199 + t276) * t202; m(6) * (t65 * t75 + t66 * t74) + m(5) * (t112 * t77 + t113 * t76) + m(4) * (-t115 * t202 - t116 * t199) * t167 + t205; m(6) * (t25 * t24 + t70 * t74 + t71 * t75) + m(5) * (t112 * t97 + t113 * t98 + t52 * t33) + m(4) * (t94 * t64 + (-t134 * t199 - t135 * t202) * t167) + t207; m(6) * (t25 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2 + t52 ^ 2) + m(4) * (t237 * t167 ^ 2 + t94 ^ 2) + t207; (-t62 - t67) * t189 + m(6) * (t38 * t65 + t39 * t66) + m(5) * (t72 * t76 + t73 * t77) + ((t61 / 0.2e1 + t49 / 0.2e1) * t202 + (t48 / 0.2e1 + t60 / 0.2e1) * t199) * t187 + t230; m(6) * (t26 * t24 + t38 * t71 + t39 * t70) + m(5) * (t63 * t33 + t72 * t98 + t73 * t97) + t206; m(6) * (t26 * t25 + t38 * t75 + t39 * t74) + m(5) * (t112 * t73 + t113 * t72 + t63 * t52) + t206; (t67 * t189 - t9) * t189 + m(6) * (t26 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(5) * (t63 ^ 2 + t72 ^ 2 + t73 ^ 2) + (t202 * t11 + t199 * t10 - t189 * (t199 * t48 + t202 * t49)) * t187 + t275; m(6) * (t65 * t68 + t66 * t69) - t260 + t230; m(6) * (t57 * t24 + t68 * t71 + t69 * t70) + t227; m(6) * (t57 * t25 + t68 * t75 + t69 * t74) + t227; m(6) * (t57 * t26 + t38 * t68 + t39 * t69) + t228; m(6) * (t57 ^ 2 + t68 ^ 2 + t69 ^ 2) + t228;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
