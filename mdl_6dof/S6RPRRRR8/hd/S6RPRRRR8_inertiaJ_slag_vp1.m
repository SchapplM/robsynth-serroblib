% Calculate joint inertia matrix for
% S6RPRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:52
% EndTime: 2019-03-09 07:19:58
% DurationCPUTime: 3.02s
% Computational Cost: add. (9797->424), mult. (10259->601), div. (0->0), fcn. (11033->10), ass. (0->213)
t208 = qJ(3) + qJ(4);
t197 = sin(t208);
t214 = cos(qJ(1));
t268 = t197 * t214;
t186 = pkin(4) * t268;
t199 = cos(t208);
t209 = sin(qJ(5));
t211 = sin(qJ(1));
t263 = t209 * t211;
t212 = cos(qJ(5));
t194 = pkin(5) * t212 + pkin(4);
t270 = t194 * t197;
t215 = -pkin(10) - pkin(9);
t284 = -pkin(9) - t215;
t207 = qJ(5) + qJ(6);
t196 = sin(t207);
t198 = cos(t207);
t267 = t198 * t211;
t148 = t196 * t268 + t267;
t266 = t198 * t214;
t149 = t196 * t211 - t197 * t266;
t235 = -rSges(7,1) * t149 - rSges(7,2) * t148;
t264 = t199 * t214;
t92 = rSges(7,3) * t264 - t235;
t294 = t92 + pkin(5) * t263 + t186 + (t284 * t199 - t270) * t214;
t293 = t211 * t214;
t292 = (rSges(5,1) * t197 + rSges(5,2) * t199) * t214;
t210 = sin(qJ(3));
t213 = cos(qJ(3));
t291 = (rSges(4,1) * t210 + rSges(4,2) * t213) * t214;
t122 = (-pkin(4) + t194) * t199 + t284 * t197;
t127 = rSges(7,3) * t197 + (rSges(7,1) * t198 - rSges(7,2) * t196) * t199;
t290 = -t122 - t127;
t205 = t211 ^ 2;
t206 = t214 ^ 2;
t269 = t197 * t211;
t146 = -t196 * t269 + t266;
t147 = t196 * t214 + t197 * t267;
t265 = t199 * t211;
t85 = Icges(7,5) * t147 + Icges(7,6) * t146 - Icges(7,3) * t265;
t87 = Icges(7,4) * t147 + Icges(7,2) * t146 - Icges(7,6) * t265;
t89 = Icges(7,1) * t147 + Icges(7,4) * t146 - Icges(7,5) * t265;
t36 = t197 * t85 + (-t196 * t87 + t198 * t89) * t199;
t86 = Icges(7,5) * t149 + Icges(7,6) * t148 + Icges(7,3) * t264;
t88 = Icges(7,4) * t149 + Icges(7,2) * t148 + Icges(7,6) * t264;
t90 = Icges(7,1) * t149 + Icges(7,4) * t148 + Icges(7,5) * t264;
t37 = t197 * t86 + (-t196 * t88 + t198 * t90) * t199;
t30 = t148 * t87 + t149 * t89 + t85 * t264;
t31 = t148 * t88 + t149 * t90 + t86 * t264;
t124 = Icges(7,3) * t197 + (Icges(7,5) * t198 - Icges(7,6) * t196) * t199;
t125 = Icges(7,6) * t197 + (Icges(7,4) * t198 - Icges(7,2) * t196) * t199;
t126 = Icges(7,5) * t197 + (Icges(7,1) * t198 - Icges(7,4) * t196) * t199;
t57 = t124 * t264 + t125 * t148 + t126 * t149;
t5 = t197 * t57 + (-t211 * t30 + t214 * t31) * t199;
t256 = t199 * t198 * t126 + t197 * t124;
t272 = t125 * t196;
t63 = (-t199 * t272 + t256) * t197;
t289 = t5 * t264 + t197 * (t63 + (-t211 * t36 + t214 * t37) * t199);
t288 = t197 / 0.2e1;
t287 = t211 / 0.2e1;
t286 = t214 / 0.2e1;
t285 = pkin(3) * t213;
t91 = t147 * rSges(7,1) + t146 * rSges(7,2) - rSges(7,3) * t265;
t185 = pkin(4) * t269;
t161 = -pkin(9) * t265 + t185;
t258 = t214 * t209;
t248 = pkin(5) * t258 + t194 * t269 + t215 * t265;
t98 = -t161 + t248;
t283 = -t91 - t98;
t281 = t36 * t214;
t280 = t37 * t211;
t259 = t212 * t214;
t163 = -t197 * t263 + t259;
t261 = t211 * t212;
t164 = t197 * t261 + t258;
t100 = Icges(6,5) * t164 + Icges(6,6) * t163 - Icges(6,3) * t265;
t102 = Icges(6,4) * t164 + Icges(6,2) * t163 - Icges(6,6) * t265;
t104 = Icges(6,1) * t164 + Icges(6,4) * t163 - Icges(6,5) * t265;
t49 = t100 * t197 + (-t102 * t209 + t104 * t212) * t199;
t279 = t49 * t214;
t165 = t197 * t258 + t261;
t166 = -t197 * t259 + t263;
t101 = Icges(6,5) * t166 + Icges(6,6) * t165 + Icges(6,3) * t264;
t103 = Icges(6,4) * t166 + Icges(6,2) * t165 + Icges(6,6) * t264;
t105 = Icges(6,1) * t166 + Icges(6,4) * t165 + Icges(6,5) * t264;
t50 = t101 * t197 + (-t103 * t209 + t105 * t212) * t199;
t278 = t50 * t211;
t121 = t211 * t127;
t68 = t199 * t121 + t197 * t91;
t236 = -rSges(6,1) * t166 - rSges(6,2) * t165;
t107 = rSges(6,3) * t264 - t236;
t145 = t214 * (pkin(9) * t264 - t186);
t277 = t214 * t107 + t145;
t276 = Icges(4,4) * t210;
t275 = Icges(4,4) * t213;
t274 = Icges(5,4) * t197;
t273 = Icges(5,4) * t199;
t131 = Icges(6,6) * t197 + (Icges(6,4) * t212 - Icges(6,2) * t209) * t199;
t271 = t131 * t209;
t262 = t210 * t211;
t116 = t211 * t122;
t133 = rSges(6,3) * t197 + (rSges(6,1) * t212 - rSges(6,2) * t209) * t199;
t128 = t211 * t133;
t260 = t211 * t213;
t254 = t164 * rSges(6,1) + t163 * rSges(6,2);
t106 = -rSges(6,3) * t265 + t254;
t257 = -t106 - t161;
t130 = Icges(6,3) * t197 + (Icges(6,5) * t212 - Icges(6,6) * t209) * t199;
t132 = Icges(6,5) * t197 + (Icges(6,1) * t212 - Icges(6,4) * t209) * t199;
t255 = t199 * t212 * t132 + t197 * t130;
t174 = pkin(4) * t199 + pkin(9) * t197;
t167 = t211 * t174;
t109 = t167 + t128;
t216 = -pkin(8) - pkin(7);
t253 = t214 * t210 * pkin(3) + t211 * t216;
t252 = t214 * pkin(1) + t211 * qJ(2);
t251 = t205 + t206;
t250 = t294 * t214 + t145;
t249 = -t161 + t283;
t77 = t167 + t116 + t121;
t142 = rSges(5,1) * t269 + rSges(5,2) * t265 + t214 * rSges(5,3);
t247 = rSges(4,1) * t262 + rSges(4,2) * t260 + t214 * rSges(4,3);
t201 = t214 * qJ(2);
t246 = t201 + t253;
t245 = t199 * (-rSges(6,3) - pkin(9));
t244 = -t265 / 0.2e1;
t243 = t264 / 0.2e1;
t242 = -t174 - t285;
t28 = t146 * t87 + t147 * t89 - t85 * t265;
t29 = t146 * t88 + t147 * t90 - t86 * t265;
t15 = t211 * t29 + t214 * t28;
t16 = t211 * t31 + t214 * t30;
t56 = -t124 * t265 + t125 * t146 + t126 * t147;
t4 = t197 * t56 + (-t211 * t28 + t214 * t29) * t199;
t241 = t15 * t244 + t16 * t243 + t4 * t286 + t5 * t287 + (t280 + t281) * t288;
t240 = t63 + (t36 + t56) * t244 + (t37 + t57) * t243;
t239 = -t4 * t265 + t289;
t228 = Icges(5,5) * t197 + Icges(5,6) * t199;
t136 = Icges(5,3) * t214 + t228 * t211;
t137 = Icges(5,3) * t211 - t228 * t214;
t42 = -t100 * t265 + t102 * t163 + t104 * t164;
t43 = -t101 * t265 + t103 * t163 + t105 * t164;
t22 = t211 * t43 + t214 * t42;
t44 = t100 * t264 + t102 * t165 + t104 * t166;
t45 = t101 * t264 + t103 * t165 + t105 * t166;
t23 = t211 * t45 + t214 * t44;
t234 = (t205 * t137 + t16 + t23) * t211 + (t206 * t136 + t15 + t22 + (t211 * t136 + t214 * t137) * t211) * t214;
t233 = Icges(4,1) * t210 + t275;
t232 = Icges(5,1) * t197 + t273;
t231 = Icges(4,2) * t213 + t276;
t230 = Icges(5,2) * t199 + t274;
t229 = Icges(4,5) * t210 + Icges(4,6) * t213;
t173 = rSges(5,1) * t199 - rSges(5,2) * t197;
t191 = pkin(3) * t260;
t134 = t173 * t211 + t191;
t135 = (-t173 - t285) * t214;
t227 = t134 * t211 - t135 * t214;
t170 = -Icges(5,2) * t197 + t273;
t171 = Icges(5,1) * t199 - t274;
t222 = t170 * t199 + t171 * t197;
t190 = pkin(3) * t262;
t221 = -t214 * t216 + t190 + t252;
t119 = t201 + t291 + (-rSges(4,3) - pkin(1) - pkin(7)) * t211;
t120 = t214 * pkin(7) + t247 + t252;
t220 = m(4) * (t119 * t211 - t120 * t214);
t111 = t292 + (-rSges(5,3) - pkin(1)) * t211 + t246;
t112 = t221 + t142;
t219 = m(5) * (t111 * t211 - t112 * t214);
t62 = t130 * t264 + t131 * t165 + t132 * t166;
t10 = t197 * t62 + (-t211 * t44 + t214 * t45) * t199;
t61 = -t130 * t265 + t131 * t163 + t132 * t164;
t9 = t197 * t61 + (-t211 * t42 + t214 * t43) * t199;
t218 = t10 * t287 + t22 * t244 + t23 * t243 + t9 * t286 + (t278 + t279) * t288 + t241;
t169 = Icges(5,5) * t199 - Icges(5,6) * t197;
t217 = t281 / 0.2e1 + t280 / 0.2e1 + t279 / 0.2e1 + t278 / 0.2e1 + (-(Icges(5,6) * t211 - t230 * t214) * t197 + (Icges(5,5) * t211 - t232 * t214) * t199 + t169 * t211 - t222 * t214 + t57 + t62) * t287 + (-(Icges(5,6) * t214 + t230 * t211) * t197 + (Icges(5,5) * t214 + t232 * t211) * t199 + t169 * t214 + t222 * t211 + t56 + t61) * t286;
t182 = rSges(2,1) * t214 - rSges(2,2) * t211;
t181 = rSges(4,1) * t213 - rSges(4,2) * t210;
t180 = -rSges(2,1) * t211 - rSges(2,2) * t214;
t162 = t190 + (-pkin(7) - t216) * t214;
t160 = -rSges(3,2) * t214 + rSges(3,3) * t211 + t252;
t159 = rSges(3,3) * t214 + t201 + (rSges(3,2) - pkin(1)) * t211;
t152 = Icges(4,3) * t211 - t229 * t214;
t151 = Icges(4,3) * t214 + t229 * t211;
t150 = t214 * (-pkin(7) * t211 - t253);
t129 = t214 * (rSges(5,3) * t211 - t292);
t115 = t127 * t264;
t110 = (-t133 - t174) * t214;
t108 = -t211 * t247 + (t211 * rSges(4,3) - t291) * t214;
t96 = (-t133 + t242) * t214;
t95 = t191 + t109;
t93 = -t142 * t211 + t129;
t78 = (-t174 + t290) * t214;
t76 = -t107 * t197 + t133 * t264;
t75 = t106 * t197 + t199 * t128;
t74 = t211 * t245 + t185 + t221 + t254;
t73 = -pkin(1) * t211 + t214 * t245 + t186 + t236 + t246;
t72 = (t242 + t290) * t214;
t71 = t191 + t77;
t70 = t129 + t150 + (-t142 - t162) * t211;
t69 = -t197 * t92 + t115;
t67 = (-t199 * t271 + t255) * t197;
t66 = t221 + t91 + t248;
t65 = (-pkin(5) * t209 - pkin(1)) * t211 + (t270 + (-rSges(7,3) + t215) * t199) * t214 + t235 + t246;
t64 = (-t106 * t214 - t107 * t211) * t199;
t58 = (-t211 * t92 - t214 * t91) * t199;
t53 = t257 * t211 + t277;
t46 = t150 + (-t162 + t257) * t211 + t277;
t39 = t122 * t264 - t197 * t294 + t115;
t38 = t199 * t116 + t197 * t98 + t68;
t27 = (-t211 * t294 + t283 * t214) * t199;
t26 = t249 * t211 + t250;
t25 = t150 + (-t162 + t249) * t211 + t250;
t1 = [-t197 * t170 - t210 * (-Icges(4,2) * t210 + t275) + t213 * (Icges(4,1) * t213 - t276) + Icges(3,1) + Icges(2,3) + (t171 - t271 - t272) * t199 + m(7) * (t65 ^ 2 + t66 ^ 2) + m(6) * (t73 ^ 2 + t74 ^ 2) + m(5) * (t111 ^ 2 + t112 ^ 2) + m(4) * (t119 ^ 2 + t120 ^ 2) + m(3) * (t159 ^ 2 + t160 ^ 2) + m(2) * (t180 ^ 2 + t182 ^ 2) + t255 + t256; m(7) * (t211 * t65 - t214 * t66) + m(6) * (t211 * t73 - t214 * t74) + t219 + t220 + m(3) * (t159 * t211 - t160 * t214); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t251; m(7) * (t65 * t71 + t66 * t72) + m(6) * (t73 * t95 + t74 * t96) + m(5) * (t111 * t134 + t112 * t135) + t181 * t220 + t217 + (t206 / 0.2e1 + t205 / 0.2e1) * (Icges(4,5) * t213 - Icges(4,6) * t210) + (-t210 * (Icges(4,6) * t211 - t231 * t214) + t213 * (Icges(4,5) * t211 - t233 * t214)) * t287 + (-t210 * (Icges(4,6) * t214 + t231 * t211) + t213 * (Icges(4,5) * t214 + t233 * t211)) * t286; m(5) * t227 + m(6) * (t211 * t95 - t214 * t96) + m(7) * (t211 * t71 - t214 * t72) + m(4) * t251 * t181; m(7) * (t25 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t46 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(5) * (t134 ^ 2 + t135 ^ 2 + t70 ^ 2) + t211 * (t151 * t293 + t205 * t152) + m(4) * (t251 * t181 ^ 2 + t108 ^ 2) + t214 * (t206 * t151 + t152 * t293) + t234; m(7) * (t65 * t77 + t66 * t78) + m(6) * (t109 * t73 + t110 * t74) + t217 + t173 * t219; m(6) * (t109 * t211 - t110 * t214) + m(7) * (t211 * t77 - t214 * t78) + m(5) * t251 * t173; m(7) * (t26 * t25 + t71 * t77 + t72 * t78) + m(6) * (t109 * t95 + t110 * t96 + t53 * t46) + m(5) * (t227 * t173 + t70 * t93) + t234; m(7) * (t26 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2 + t53 ^ 2) + m(5) * (t251 * t173 ^ 2 + t93 ^ 2) + t234; t67 + m(7) * (t38 * t66 + t39 * t65) + m(6) * (t73 * t76 + t74 * t75) + ((t62 / 0.2e1 + t50 / 0.2e1) * t214 + (-t61 / 0.2e1 - t49 / 0.2e1) * t211) * t199 + t240; m(6) * (t211 * t76 - t214 * t75) + m(7) * (t211 * t39 - t214 * t38); m(7) * (t25 * t27 + t38 * t72 + t39 * t71) + m(6) * (t64 * t46 + t75 * t96 + t76 * t95) + t218; m(7) * (t27 * t26 + t38 * t78 + t39 * t77) + m(6) * (t109 * t76 + t110 * t75 + t64 * t53) + t218; t197 * t67 + m(7) * (t27 ^ 2 + t38 ^ 2 + t39 ^ 2) + m(6) * (t64 ^ 2 + t75 ^ 2 + t76 ^ 2) + ((t197 * t50 + t10) * t214 + (-t197 * t49 - t4 - t9) * t211) * t199 + t289; m(7) * (t65 * t69 + t66 * t68) + t240; m(7) * (t211 * t69 - t214 * t68); m(7) * (t25 * t58 + t68 * t72 + t69 * t71) + t241; m(7) * (t58 * t26 + t68 * t78 + t69 * t77) + t241; m(7) * (t27 * t58 + t38 * t68 + t39 * t69) + t239; m(7) * (t58 ^ 2 + t68 ^ 2 + t69 ^ 2) + t239;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
