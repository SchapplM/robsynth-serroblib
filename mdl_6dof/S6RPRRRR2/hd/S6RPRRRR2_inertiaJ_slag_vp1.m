% Calculate joint inertia matrix for
% S6RPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Mq [6x6]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:57:02
% EndTime: 2019-03-09 06:57:09
% DurationCPUTime: 2.96s
% Computational Cost: add. (14041->412), mult. (10139->575), div. (0->0), fcn. (10945->12), ass. (0->213)
t206 = cos(qJ(5));
t189 = pkin(5) * t206 + pkin(4);
t209 = -pkin(10) - pkin(9);
t200 = qJ(1) + pkin(11);
t194 = cos(t200);
t202 = qJ(3) + qJ(4);
t198 = cos(t202);
t258 = t194 * t198;
t196 = sin(t202);
t259 = t194 * t196;
t193 = sin(t200);
t203 = sin(qJ(5));
t260 = t193 * t203;
t201 = qJ(5) + qJ(6);
t197 = cos(t201);
t195 = sin(t201);
t254 = t195 * t198;
t150 = t193 * t197 - t194 * t254;
t253 = t197 * t198;
t151 = t193 * t195 + t194 * t253;
t94 = rSges(7,1) * t151 + rSges(7,2) * t150 + rSges(7,3) * t259;
t286 = pkin(5) * t260 + t189 * t258 - t209 * t259 + t94;
t245 = pkin(4) * t258 + pkin(9) * t259;
t285 = -t245 + t286;
t222 = Icges(5,5) * t198 - Icges(5,6) * t196;
t126 = -Icges(5,3) * t194 + t193 * t222;
t127 = Icges(5,3) * t193 + t194 * t222;
t148 = -t193 * t254 - t194 * t197;
t149 = t193 * t253 - t194 * t195;
t261 = t193 * t196;
t87 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t261;
t89 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t261;
t91 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t261;
t27 = t148 * t89 + t149 * t91 + t261 * t87;
t88 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t259;
t90 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t259;
t92 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t259;
t28 = t148 * t90 + t149 * t92 + t261 * t88;
t15 = t193 * t28 - t194 * t27;
t192 = t194 ^ 2;
t252 = t198 * t203;
t158 = -t193 * t252 - t194 * t206;
t251 = t198 * t206;
t257 = t194 * t203;
t159 = t193 * t251 - t257;
t102 = Icges(6,5) * t159 + Icges(6,6) * t158 + Icges(6,3) * t261;
t104 = Icges(6,4) * t159 + Icges(6,2) * t158 + Icges(6,6) * t261;
t106 = Icges(6,1) * t159 + Icges(6,4) * t158 + Icges(6,5) * t261;
t36 = t102 * t261 + t104 * t158 + t106 * t159;
t160 = t193 * t206 - t194 * t252;
t161 = t194 * t251 + t260;
t103 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t259;
t105 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t259;
t107 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t259;
t37 = t103 * t261 + t105 * t158 + t107 * t159;
t20 = t193 * t37 - t194 * t36;
t262 = Icges(5,4) * t198;
t224 = -Icges(5,2) * t196 + t262;
t129 = Icges(5,6) * t193 + t194 * t224;
t263 = Icges(5,4) * t196;
t226 = Icges(5,1) * t198 - t263;
t131 = Icges(5,5) * t193 + t194 * t226;
t220 = -t129 * t196 + t131 * t198;
t128 = -Icges(5,6) * t194 + t193 * t224;
t130 = -Icges(5,5) * t194 + t193 * t226;
t221 = t128 * t196 - t130 * t198;
t284 = -t15 - t20 - t192 * t126 - (t220 * t193 + (-t127 + t221) * t194) * t193;
t191 = t193 ^ 2;
t132 = -Icges(7,3) * t198 + (Icges(7,5) * t197 - Icges(7,6) * t195) * t196;
t133 = -Icges(7,6) * t198 + (Icges(7,4) * t197 - Icges(7,2) * t195) * t196;
t134 = -Icges(7,5) * t198 + (Icges(7,1) * t197 - Icges(7,4) * t195) * t196;
t56 = t132 * t261 + t133 * t148 + t134 * t149;
t5 = -t198 * t56 + (t193 * t27 + t194 * t28) * t196;
t29 = t150 * t89 + t151 * t91 + t259 * t87;
t30 = t150 * t90 + t151 * t92 + t259 * t88;
t57 = t132 * t259 + t133 * t150 + t134 * t151;
t6 = -t198 * t57 + (t193 * t29 + t194 * t30) * t196;
t283 = t259 * t6 + t261 * t5;
t282 = t193 / 0.2e1;
t281 = -t194 / 0.2e1;
t280 = -t198 / 0.2e1;
t204 = sin(qJ(3));
t279 = pkin(3) * t204;
t278 = pkin(4) * t198;
t205 = sin(qJ(1));
t277 = t205 * pkin(1);
t276 = -pkin(4) + t189;
t275 = pkin(9) + t209;
t207 = cos(qJ(3));
t274 = rSges(4,1) * t207;
t273 = rSges(4,2) * t204;
t272 = t194 * rSges(4,3);
t42 = -t198 * t87 + (-t195 * t89 + t197 * t91) * t196;
t271 = t42 * t194;
t43 = -t198 * t88 + (-t195 * t90 + t197 * t92) * t196;
t270 = t43 * t193;
t48 = -t102 * t198 + (-t104 * t203 + t106 * t206) * t196;
t269 = t48 * t194;
t49 = -t103 * t198 + (-t105 * t203 + t107 * t206) * t196;
t268 = t49 * t193;
t118 = t196 * t197 * t134;
t255 = t195 * t133;
t66 = -t198 * t132 - t196 * t255 + t118;
t267 = t66 * t198;
t135 = -rSges(7,3) * t198 + (rSges(7,1) * t197 - rSges(7,2) * t195) * t196;
t228 = -t149 * rSges(7,1) - t148 * rSges(7,2);
t93 = rSges(7,3) * t261 - t228;
t67 = t135 * t261 + t198 * t93;
t265 = Icges(4,4) * t204;
t264 = Icges(4,4) * t207;
t210 = -pkin(8) - pkin(7);
t256 = t194 * t210;
t153 = -Icges(6,6) * t198 + (Icges(6,4) * t206 - Icges(6,2) * t203) * t196;
t250 = t203 * t153;
t190 = pkin(3) * t207 + pkin(2);
t171 = t194 * t190;
t188 = t194 * pkin(7);
t249 = t193 * (t256 + t188 + (-pkin(2) + t190) * t193) + t194 * (-t194 * pkin(2) + t171 + (-pkin(7) - t210) * t193);
t215 = rSges(5,1) * t258 - rSges(5,2) * t259 + rSges(5,3) * t193;
t230 = rSges(5,1) * t198 - rSges(5,2) * t196;
t86 = t193 * (-t194 * rSges(5,3) + t193 * t230) + t194 * t215;
t125 = t196 * t276 + t198 * t275;
t248 = -t125 - t135;
t247 = t191 * (pkin(9) * t196 + t278) + t194 * t245;
t155 = -rSges(6,3) * t198 + (rSges(6,1) * t206 - rSges(6,2) * t203) * t196;
t170 = pkin(4) * t196 - pkin(9) * t198;
t246 = -t155 - t170;
t244 = rSges(4,3) * t193 + t194 * t274;
t243 = t191 + t192;
t16 = t193 * t30 - t194 * t29;
t38 = t102 * t259 + t104 * t160 + t106 * t161;
t39 = t103 * t259 + t105 * t160 + t107 * t161;
t21 = t193 * t39 - t194 * t38;
t242 = (t191 * t127 + t16 + t21 + (t221 * t194 + (-t126 + t220) * t193) * t194) * t193;
t241 = -t170 + t248;
t111 = rSges(6,1) * t161 + rSges(6,2) * t160 + rSges(6,3) * t259;
t240 = t261 / 0.2e1;
t239 = t259 / 0.2e1;
t169 = rSges(5,1) * t196 + rSges(5,2) * t198;
t238 = -t169 - t279;
t237 = -t170 - t279;
t236 = (t42 + t56) * t240 + (t43 + t57) * t239;
t229 = -t159 * rSges(6,1) - t158 * rSges(6,2);
t110 = rSges(6,3) * t261 - t229;
t50 = t110 * t193 + t111 * t194 + t247;
t11 = -t267 + (t193 * t42 + t194 * t43) * t196;
t235 = -t198 * t11 + t283;
t234 = t15 * t240 + t16 * t239 + t5 * t281 + t6 * t282 + (t270 - t271) * t280;
t233 = -t155 + t237;
t208 = cos(qJ(1));
t199 = t208 * pkin(1);
t232 = -t193 * t210 + t171 + t199;
t231 = -t273 + t274;
t227 = Icges(4,1) * t207 - t265;
t225 = -Icges(4,2) * t204 + t264;
t223 = Icges(4,5) * t207 - Icges(4,6) * t204;
t167 = Icges(5,2) * t198 + t263;
t168 = Icges(5,1) * t196 + t262;
t217 = -t167 * t196 + t168 * t198;
t216 = t237 + t248;
t100 = -pkin(5) * t257 + (-t196 * t275 + t198 * t276) * t193;
t25 = t247 + t285 * t194 + (t100 + t93) * t193;
t213 = t194 * t284 + t242;
t152 = -Icges(6,3) * t198 + (Icges(6,5) * t206 - Icges(6,6) * t203) * t196;
t154 = -Icges(6,5) * t198 + (Icges(6,1) * t206 - Icges(6,4) * t203) * t196;
t63 = t152 * t259 + t153 * t160 + t154 * t161;
t10 = -t198 * t63 + (t193 * t38 + t194 * t39) * t196;
t62 = t152 * t261 + t153 * t158 + t154 * t159;
t9 = -t198 * t62 + (t193 * t36 + t194 * t37) * t196;
t212 = t10 * t282 + t20 * t240 + t21 * t239 + t9 * t281 + (t268 - t269) * t280 + t234;
t166 = Icges(5,5) * t196 + Icges(5,6) * t198;
t211 = t270 / 0.2e1 - t271 / 0.2e1 + t268 / 0.2e1 - t269 / 0.2e1 + (t129 * t198 + t131 * t196 + t166 * t193 + t194 * t217 + t57 + t63) * t282 + (t128 * t198 + t130 * t196 - t166 * t194 + t193 * t217 + t56 + t62) * t281;
t183 = rSges(2,1) * t208 - rSges(2,2) * t205;
t182 = -rSges(2,1) * t205 - rSges(2,2) * t208;
t181 = rSges(4,1) * t204 + rSges(4,2) * t207;
t163 = rSges(3,1) * t194 - rSges(3,2) * t193 + t199;
t162 = -rSges(3,1) * t193 - rSges(3,2) * t194 - t277;
t143 = Icges(4,3) * t193 + t194 * t223;
t142 = -Icges(4,3) * t194 + t193 * t223;
t137 = t238 * t194;
t136 = t238 * t193;
t124 = t196 * t206 * t154;
t117 = t193 * pkin(7) + t199 + (pkin(2) - t273) * t194 + t244;
t116 = t272 - t277 + t188 + (-pkin(2) - t231) * t193;
t115 = t246 * t194;
t114 = t246 * t193;
t113 = t215 + t232;
t112 = -t277 + (rSges(5,3) - t210) * t194 + (-t190 - t230) * t193;
t109 = t233 * t194;
t108 = t233 * t193;
t99 = t194 * (-t194 * t273 + t244) + (t193 * t231 - t272) * t193;
t80 = t93 * t259;
t77 = t241 * t194;
t76 = t241 * t193;
t75 = -t198 * t152 - t196 * t250 + t124;
t74 = -t111 * t198 - t155 * t259;
t73 = t110 * t198 + t155 * t261;
t72 = t216 * t194;
t71 = t216 * t193;
t70 = t232 + t111 + t245;
t69 = -t277 - t256 + (-t278 - t190 + (-rSges(6,3) - pkin(9)) * t196) * t193 + t229;
t68 = -t135 * t259 - t198 * t94;
t65 = t232 + t286;
t64 = -t277 + (pkin(5) * t203 - t210) * t194 + (-t189 * t198 - t190 + (-rSges(7,3) + t209) * t196) * t193 + t228;
t61 = (t110 * t194 - t111 * t193) * t196;
t58 = t86 + t249;
t55 = -t261 * t94 + t80;
t45 = -t198 * t285 + t248 * t259;
t44 = t100 * t198 + t125 * t261 + t67;
t31 = t50 + t249;
t26 = t80 + (t100 * t194 - t193 * t285) * t196;
t24 = t25 + t249;
t1 = [t207 * (Icges(4,2) * t207 + t265) + t204 * (Icges(4,1) * t204 + t264) + Icges(2,3) + Icges(3,3) + t118 + t124 + (-t152 - t132 + t167) * t198 + (t168 - t250 - t255) * t196 + m(7) * (t64 ^ 2 + t65 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(5) * (t112 ^ 2 + t113 ^ 2) + m(4) * (t116 ^ 2 + t117 ^ 2) + m(3) * (t162 ^ 2 + t163 ^ 2) + m(2) * (t182 ^ 2 + t183 ^ 2); 0; m(3) + m(4) + m(5) + m(6) + m(7); m(4) * (-t116 * t194 - t117 * t193) * t181 + t211 + ((Icges(4,6) * t193 + t194 * t225) * t207 + (Icges(4,5) * t193 + t194 * t227) * t204) * t282 + ((-Icges(4,6) * t194 + t193 * t225) * t207 + (-Icges(4,5) * t194 + t193 * t227) * t204) * t281 + (t191 / 0.2e1 + t192 / 0.2e1) * (Icges(4,5) * t204 + Icges(4,6) * t207) + m(7) * (t64 * t72 + t65 * t71) + m(6) * (t108 * t70 + t109 * t69) + m(5) * (t112 * t137 + t113 * t136); m(4) * t99 + m(5) * t58 + m(6) * t31 + m(7) * t24; m(7) * (t24 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t108 ^ 2 + t109 ^ 2 + t31 ^ 2) + m(5) * (t136 ^ 2 + t137 ^ 2 + t58 ^ 2) + m(4) * (t181 ^ 2 * t243 + t99 ^ 2) + t193 * t191 * t143 + t242 + (-t192 * t142 + (-t193 * t142 + t194 * t143) * t193 + t284) * t194; t211 + m(5) * (-t112 * t194 - t113 * t193) * t169 + m(7) * (t64 * t77 + t65 * t76) + m(6) * (t114 * t70 + t115 * t69); m(5) * t86 + m(6) * t50 + m(7) * t25; m(7) * (t24 * t25 + t71 * t76 + t72 * t77) + m(6) * (t108 * t114 + t109 * t115 + t31 * t50) + m(5) * (t58 * t86 + (-t136 * t193 - t137 * t194) * t169) + t213; m(7) * (t25 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(6) * (t114 ^ 2 + t115 ^ 2 + t50 ^ 2) + m(5) * (t169 ^ 2 * t243 + t86 ^ 2) + t213; (-t75 - t66) * t198 + m(7) * (t44 * t64 + t45 * t65) + m(6) * (t69 * t73 + t70 * t74) + ((t63 / 0.2e1 + t49 / 0.2e1) * t194 + (t62 / 0.2e1 + t48 / 0.2e1) * t193) * t196 + t236; m(6) * t61 + m(7) * t26; m(7) * (t24 * t26 + t44 * t72 + t45 * t71) + m(6) * (t108 * t74 + t109 * t73 + t31 * t61) + t212; m(7) * (t25 * t26 + t44 * t77 + t45 * t76) + m(6) * (t114 * t74 + t115 * t73 + t50 * t61) + t212; (t75 * t198 - t11) * t198 + m(7) * (t26 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t61 ^ 2 + t73 ^ 2 + t74 ^ 2) + (t194 * t10 + t193 * t9 - t198 * (t193 * t48 + t194 * t49)) * t196 + t283; m(7) * (t64 * t67 + t65 * t68) - t267 + t236; m(7) * t55; m(7) * (t24 * t55 + t67 * t72 + t68 * t71) + t234; m(7) * (t25 * t55 + t67 * t77 + t68 * t76) + t234; m(7) * (t26 * t55 + t44 * t67 + t45 * t68) + t235; m(7) * (t55 ^ 2 + t67 ^ 2 + t68 ^ 2) + t235;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
