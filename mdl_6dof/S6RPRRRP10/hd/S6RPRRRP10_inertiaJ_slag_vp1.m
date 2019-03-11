% Calculate joint inertia matrix for
% S6RPRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
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
% Datum: 2019-03-09 06:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRRRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP10_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:30:19
% EndTime: 2019-03-09 06:30:29
% DurationCPUTime: 4.14s
% Computational Cost: add. (7578->420), mult. (11052->599), div. (0->0), fcn. (12022->8), ass. (0->206)
t214 = qJ(4) + qJ(5);
t205 = sin(t214);
t206 = cos(t214);
t220 = cos(qJ(1));
t264 = t220 * t206;
t216 = sin(qJ(3));
t217 = sin(qJ(1));
t271 = t216 * t217;
t168 = t205 * t271 - t264;
t169 = t205 * t220 + t206 * t271;
t219 = cos(qJ(3));
t268 = t217 * t219;
t92 = Icges(7,5) * t169 - Icges(7,6) * t268 + Icges(7,3) * t168;
t98 = Icges(6,4) * t169 - Icges(6,2) * t168 - Icges(6,6) * t268;
t315 = t92 - t98;
t270 = t216 * t220;
t170 = t205 * t270 + t206 * t217;
t172 = t205 * t217 - t216 * t264;
t266 = t219 * t220;
t93 = Icges(7,5) * t172 + Icges(7,6) * t266 - Icges(7,3) * t170;
t99 = Icges(6,4) * t172 + Icges(6,2) * t170 + Icges(6,6) * t266;
t314 = t93 - t99;
t94 = Icges(6,5) * t169 - Icges(6,6) * t168 - Icges(6,3) * t268;
t96 = Icges(7,4) * t169 - Icges(7,2) * t268 + Icges(7,6) * t168;
t313 = t94 + t96;
t95 = Icges(6,5) * t172 + Icges(6,6) * t170 + Icges(6,3) * t266;
t97 = Icges(7,4) * t172 + Icges(7,2) * t266 - Icges(7,6) * t170;
t312 = t95 + t97;
t100 = Icges(7,1) * t169 - Icges(7,4) * t268 + Icges(7,5) * t168;
t102 = Icges(6,1) * t169 - Icges(6,4) * t168 - Icges(6,5) * t268;
t311 = t100 + t102;
t101 = Icges(7,1) * t172 + Icges(7,4) * t266 - Icges(7,5) * t170;
t103 = Icges(6,1) * t172 + Icges(6,4) * t170 + Icges(6,5) * t266;
t310 = t101 + t103;
t309 = rSges(7,1) + pkin(5);
t308 = t315 * t168 + t311 * t169 - t313 * t268;
t307 = t314 * t168 + t310 * t169 - t312 * t268;
t306 = -t315 * t170 + t311 * t172 + t313 * t266;
t305 = -t314 * t170 + t310 * t172 + t312 * t266;
t144 = Icges(7,6) * t216 + (Icges(7,5) * t206 + Icges(7,3) * t205) * t219;
t146 = Icges(7,2) * t216 + (Icges(7,4) * t206 + Icges(7,6) * t205) * t219;
t148 = Icges(7,4) * t216 + (Icges(7,1) * t206 + Icges(7,5) * t205) * t219;
t67 = t144 * t168 - t146 * t268 + t148 * t169;
t145 = Icges(6,3) * t216 + (Icges(6,5) * t206 - Icges(6,6) * t205) * t219;
t147 = Icges(6,6) * t216 + (Icges(6,4) * t206 - Icges(6,2) * t205) * t219;
t149 = Icges(6,5) * t216 + (Icges(6,1) * t206 - Icges(6,4) * t205) * t219;
t68 = -t145 * t268 - t147 * t168 + t149 * t169;
t304 = t67 + t68;
t69 = -t144 * t170 + t146 * t266 + t148 * t172;
t70 = t145 * t266 + t147 * t170 + t149 * t172;
t303 = t70 + t69;
t302 = rSges(7,3) + qJ(6);
t301 = Icges(4,5) * t219;
t300 = Icges(4,6) * t216;
t275 = t205 * t219;
t299 = t144 * t275 + (t148 + t149) * t206 * t219 + (t145 + t146) * t216;
t298 = t301 / 0.2e1 - t300 / 0.2e1;
t261 = rSges(7,2) * t266 - t302 * t170 + t309 * t172;
t297 = (-t308 * t217 + t307 * t220) * t219 + t304 * t216;
t295 = (-t306 * t217 + t305 * t220) * t219 + t303 * t216;
t294 = t307 * t217 + t308 * t220;
t293 = t305 * t217 + t306 * t220;
t46 = t216 * t96 + (t100 * t206 + t205 * t92) * t219;
t48 = t216 * t94 + (t102 * t206 - t205 * t98) * t219;
t292 = t46 + t48;
t47 = t216 * t97 + (t101 * t206 + t205 * t93) * t219;
t49 = t216 * t95 + (t103 * t206 - t205 * t99) * t219;
t291 = t47 + t49;
t290 = (-t147 * t275 + t299) * t216;
t288 = -rSges(7,2) * t268 + t302 * t168 + t309 * t169;
t252 = rSges(7,2) * t216 + (t302 * t205 + t309 * t206) * t219;
t287 = (rSges(4,1) * t216 + rSges(4,2) * t219) * t220;
t212 = t217 ^ 2;
t213 = t220 ^ 2;
t286 = -pkin(1) - pkin(7);
t283 = t216 / 0.2e1;
t282 = t217 / 0.2e1;
t280 = t220 / 0.2e1;
t192 = rSges(4,1) * t219 - rSges(4,2) * t216;
t279 = m(4) * t192;
t105 = t169 * rSges(6,1) - t168 * rSges(6,2) - rSges(6,3) * t268;
t151 = rSges(6,3) * t216 + (rSges(6,1) * t206 - rSges(6,2) * t205) * t219;
t81 = t216 * t105 + t151 * t268;
t215 = sin(qJ(4));
t218 = cos(qJ(4));
t161 = Icges(5,6) * t216 + (Icges(5,4) * t218 - Icges(5,2) * t215) * t219;
t276 = t161 * t215;
t273 = t215 * t217;
t272 = t215 * t220;
t269 = t217 * t218;
t267 = t218 * t220;
t221 = -pkin(9) - pkin(8);
t265 = t219 * t221;
t201 = pkin(3) * t271;
t183 = -pkin(8) * t268 + t201;
t204 = pkin(4) * t218 + pkin(3);
t241 = pkin(4) * t272 + t204 * t271 + t217 * t265;
t122 = -t183 + t241;
t262 = -t105 - t122;
t107 = rSges(6,1) * t172 + rSges(6,2) * t170 + rSges(6,3) * t266;
t203 = pkin(3) * t270;
t234 = pkin(8) * t266 - t203;
t250 = t204 * t270 + t220 * t265;
t123 = pkin(4) * t273 - t234 - t250;
t260 = -t107 - t123;
t143 = (-pkin(3) + t204) * t219 + (-pkin(8) - t221) * t216;
t259 = t216 * t122 + t143 * t268;
t176 = t220 * t234;
t258 = t220 * t123 + t176;
t256 = t252 * t266;
t194 = t219 * pkin(3) + t216 * pkin(8);
t184 = t217 * t194;
t255 = t217 * t143 + t184;
t158 = Icges(5,3) * t216 + (Icges(5,5) * t218 - Icges(5,6) * t215) * t219;
t164 = Icges(5,5) * t216 + (Icges(5,1) * t218 - Icges(5,4) * t215) * t219;
t254 = t219 * t218 * t164 + t216 * t158;
t253 = -t143 - t194;
t179 = -t215 * t271 + t267;
t180 = t216 * t269 + t272;
t251 = t180 * rSges(5,1) + t179 * rSges(5,2);
t249 = t220 * pkin(1) + t217 * qJ(2);
t248 = t212 + t213;
t114 = Icges(5,5) * t180 + Icges(5,6) * t179 - Icges(5,3) * t268;
t116 = Icges(5,4) * t180 + Icges(5,2) * t179 - Icges(5,6) * t268;
t118 = Icges(5,1) * t180 + Icges(5,4) * t179 - Icges(5,5) * t268;
t58 = t114 * t216 + (-t116 * t215 + t118 * t218) * t219;
t72 = -t158 * t268 + t161 * t179 + t164 * t180;
t246 = -t58 / 0.2e1 - t72 / 0.2e1;
t181 = t215 * t270 + t269;
t182 = -t216 * t267 + t273;
t115 = Icges(5,5) * t182 + Icges(5,6) * t181 + Icges(5,3) * t266;
t117 = Icges(5,4) * t182 + Icges(5,2) * t181 + Icges(5,6) * t266;
t119 = Icges(5,1) * t182 + Icges(5,4) * t181 + Icges(5,5) * t266;
t59 = t115 * t216 + (-t117 * t215 + t119 * t218) * t219;
t73 = t158 * t266 + t161 * t181 + t164 * t182;
t245 = t59 / 0.2e1 + t73 / 0.2e1;
t244 = -t122 - t288;
t243 = -t123 - t261;
t240 = rSges(4,1) * t271 + rSges(4,2) * t268 + t220 * rSges(4,3);
t239 = t220 * pkin(7) + t249;
t238 = (-rSges(5,3) - pkin(8)) * t219;
t237 = -t268 / 0.2e1;
t236 = t266 / 0.2e1;
t235 = t295 * t266 + ((-t292 * t217 + t291 * t220) * t219 + t290) * t216;
t56 = t288 * t216 + t252 * t268;
t232 = -rSges(5,1) * t182 - rSges(5,2) * t181;
t229 = Icges(4,5) * t216 + Icges(4,6) * t219;
t226 = t239 + t241;
t225 = -t268 * t297 + t235;
t224 = (t291 * t217 + t292 * t220) * t283 + t295 * t282 + t297 * t280 + t294 * t237 + t293 * t236;
t223 = (t292 + t304) * t237 + (t291 + t303) * t236 + t290;
t208 = t220 * qJ(2);
t222 = t208 + (-pkin(4) * t215 + t286) * t217 + t250;
t193 = rSges(2,1) * t220 - rSges(2,2) * t217;
t191 = -rSges(2,1) * t217 - rSges(2,2) * t220;
t188 = -t300 + t301;
t175 = -rSges(3,2) * t220 + rSges(3,3) * t217 + t249;
t174 = rSges(3,3) * t220 + t208 + (rSges(3,2) - pkin(1)) * t217;
t167 = t216 * rSges(5,3) + (rSges(5,1) * t218 - rSges(5,2) * t215) * t219;
t160 = Icges(4,3) * t217 - t220 * t229;
t159 = Icges(4,3) * t220 + t217 * t229;
t136 = t151 * t266;
t132 = t143 * t266;
t127 = t239 + t240;
t126 = t208 + t287 + (-rSges(4,3) + t286) * t217;
t125 = (-t167 - t194) * t220;
t124 = t167 * t217 + t184;
t121 = rSges(5,3) * t266 - t232;
t120 = -rSges(5,3) * t268 + t251;
t109 = -t217 * t240 + (t217 * rSges(4,3) - t287) * t220;
t89 = (-t151 + t253) * t220;
t88 = t151 * t217 + t255;
t87 = t217 * t238 + t201 + t239 + t251;
t86 = t217 * t286 + t220 * t238 + t203 + t208 + t232;
t85 = -t121 * t216 + t167 * t266;
t84 = t120 * t216 + t167 * t268;
t83 = (-t219 * t276 + t254) * t216;
t82 = -t107 * t216 + t136;
t80 = (-t252 + t253) * t220;
t79 = t217 * t252 + t255;
t78 = t105 + t226;
t77 = -t107 + t222;
t74 = (-t120 * t220 - t121 * t217) * t219;
t71 = (-t105 * t220 - t107 * t217) * t219;
t62 = t121 * t220 + t176 + (-t120 - t183) * t217;
t61 = t226 + t288;
t60 = t222 - t261;
t57 = -t216 * t261 + t256;
t55 = t216 * t260 + t132 + t136;
t54 = t259 + t81;
t53 = t115 * t266 + t117 * t181 + t119 * t182;
t52 = t114 * t266 + t116 * t181 + t118 * t182;
t51 = -t115 * t268 + t117 * t179 + t119 * t180;
t50 = -t114 * t268 + t116 * t179 + t118 * t180;
t33 = (t217 * t260 + t220 * t262) * t219;
t32 = (-t217 * t261 - t220 * t288) * t219;
t31 = t216 * t243 + t132 + t256;
t30 = t56 + t259;
t29 = t107 * t220 + (-t183 + t262) * t217 + t258;
t28 = (t217 * t243 + t220 * t244) * t219;
t27 = t261 * t220 + (-t183 + t244) * t217 + t258;
t26 = t217 * t53 + t220 * t52;
t25 = t217 * t51 + t220 * t50;
t14 = t216 * t73 + (-t217 * t52 + t220 * t53) * t219;
t13 = t216 * t72 + (-t217 * t50 + t220 * t51) * t219;
t1 = [Icges(3,1) + Icges(2,3) + (Icges(4,1) * t219 - t147 * t205 - t276) * t219 + m(7) * (t60 ^ 2 + t61 ^ 2) + m(6) * (t77 ^ 2 + t78 ^ 2) + m(5) * (t86 ^ 2 + t87 ^ 2) + m(4) * (t126 ^ 2 + t127 ^ 2) + m(3) * (t174 ^ 2 + t175 ^ 2) + m(2) * (t191 ^ 2 + t193 ^ 2) + t254 + (-0.2e1 * Icges(4,4) * t219 + Icges(4,2) * t216) * t216 + t299; m(7) * (t217 * t60 - t220 * t61) + m(6) * (t217 * t77 - t220 * t78) + m(5) * (t217 * t86 - t220 * t87) + m(4) * (t126 * t217 - t127 * t220) + m(3) * (t174 * t217 - t175 * t220); 0.2e1 * (m(3) / 0.2e1 + m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t248; m(7) * (t60 * t79 + t61 * t80) + m(6) * (t77 * t88 + t78 * t89) + m(5) * (t124 * t86 + t125 * t87) + (t46 / 0.2e1 + t48 / 0.2e1 + t68 / 0.2e1 + t67 / 0.2e1 - t127 * t279 + t188 * t280 - t246 + t298 * t220) * t220 + (t47 / 0.2e1 + t49 / 0.2e1 + t70 / 0.2e1 + t69 / 0.2e1 + t126 * t279 + t188 * t282 + t245 + t298 * t217) * t217; m(5) * (t124 * t217 - t125 * t220) + m(6) * (t217 * t88 - t220 * t89) + m(7) * (t217 * t79 - t220 * t80) + t248 * t279; m(7) * (t27 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t29 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(5) * (t124 ^ 2 + t125 ^ 2 + t62 ^ 2) + m(4) * (t192 ^ 2 * t248 + t109 ^ 2) + (t213 * t159 + t25 + t294) * t220 + (t212 * t160 + t26 + (t217 * t159 + t220 * t160) * t220 + t293) * t217; t83 + (t217 * t246 + t220 * t245) * t219 + m(7) * (t30 * t61 + t31 * t60) + m(6) * (t54 * t78 + t55 * t77) + m(5) * (t84 * t87 + t85 * t86) + t223; m(5) * (t217 * t85 - t220 * t84) + m(6) * (t217 * t55 - t220 * t54) + m(7) * (t217 * t31 - t220 * t30); t14 * t282 + t13 * t280 + (t59 * t217 + t58 * t220) * t283 + m(7) * (t28 * t27 + t30 * t80 + t31 * t79) + m(6) * (t33 * t29 + t54 * t89 + t55 * t88) + m(5) * (t124 * t85 + t125 * t84 + t74 * t62) + (t26 * t280 - t217 * t25 / 0.2e1) * t219 + t224; t216 * t83 + m(7) * (t28 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t33 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(5) * (t74 ^ 2 + t84 ^ 2 + t85 ^ 2) + ((t216 * t59 + t14) * t220 + (-t216 * t58 - t13 - t297) * t217) * t219 + t235; m(7) * (t56 * t61 + t57 * t60) + m(6) * (t77 * t82 + t78 * t81) + t223; m(6) * (t217 * t82 - t220 * t81) + m(7) * (t217 * t57 - t220 * t56); m(7) * (t32 * t27 + t56 * t80 + t57 * t79) + m(6) * (t71 * t29 + t81 * t89 + t82 * t88) + t224; m(7) * (t28 * t32 + t30 * t56 + t31 * t57) + m(6) * (t71 * t33 + t54 * t81 + t55 * t82) + t225; m(7) * (t32 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t71 ^ 2 + t81 ^ 2 + t82 ^ 2) + t225; m(7) * (t168 * t60 - t170 * t61); m(7) * (t168 * t217 + t170 * t220); m(7) * (t168 * t79 - t170 * t80 + t27 * t275); m(7) * (t168 * t31 - t170 * t30 + t275 * t28); m(7) * (t168 * t57 - t170 * t56 + t275 * t32); m(7) * (t205 ^ 2 * t219 ^ 2 + t168 ^ 2 + t170 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
