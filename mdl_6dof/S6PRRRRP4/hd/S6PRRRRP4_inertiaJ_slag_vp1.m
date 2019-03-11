% Calculate joint inertia matrix for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:11
% EndTime: 2019-03-09 00:13:22
% DurationCPUTime: 6.89s
% Computational Cost: add. (31101->607), mult. (64382->845), div. (0->0), fcn. (82784->12), ass. (0->277)
t347 = rSges(7,1) + pkin(5);
t346 = rSges(7,3) + qJ(6);
t273 = cos(pkin(6));
t275 = sin(qJ(3));
t276 = sin(qJ(2));
t271 = sin(pkin(6));
t330 = cos(qJ(3));
t296 = t271 * t330;
t263 = t273 * t275 + t276 * t296;
t319 = qJ(4) + qJ(5);
t269 = sin(t319);
t292 = cos(t319);
t278 = cos(qJ(2));
t322 = t271 * t278;
t244 = t263 * t269 + t292 * t322;
t245 = t263 * t292 - t269 * t322;
t323 = t271 * t275;
t262 = -t273 * t330 + t276 * t323;
t186 = Icges(7,5) * t245 + Icges(7,6) * t262 + Icges(7,3) * t244;
t188 = Icges(7,4) * t245 + Icges(7,2) * t262 + Icges(7,6) * t244;
t190 = Icges(7,1) * t245 + Icges(7,4) * t262 + Icges(7,5) * t244;
t270 = sin(pkin(11));
t272 = cos(pkin(11));
t321 = t273 * t276;
t259 = t270 * t278 + t272 * t321;
t247 = t259 * t330 - t272 * t323;
t320 = t273 * t278;
t258 = t270 * t276 - t272 * t320;
t217 = t247 * t269 - t258 * t292;
t218 = t247 * t292 + t258 * t269;
t246 = t259 * t275 + t272 * t296;
t105 = t186 * t217 + t188 * t246 + t190 * t218;
t187 = Icges(6,5) * t245 - Icges(6,6) * t244 + Icges(6,3) * t262;
t189 = Icges(6,4) * t245 - Icges(6,2) * t244 + Icges(6,6) * t262;
t191 = Icges(6,1) * t245 - Icges(6,4) * t244 + Icges(6,5) * t262;
t106 = t187 * t246 - t189 * t217 + t191 * t218;
t261 = -t270 * t321 + t272 * t278;
t248 = t261 * t275 - t270 * t296;
t150 = Icges(7,5) * t218 + Icges(7,6) * t246 + Icges(7,3) * t217;
t154 = Icges(7,4) * t218 + Icges(7,2) * t246 + Icges(7,6) * t217;
t158 = Icges(7,1) * t218 + Icges(7,4) * t246 + Icges(7,5) * t217;
t81 = t150 * t217 + t154 * t246 + t158 * t218;
t249 = t261 * t330 + t270 * t323;
t260 = t270 * t320 + t272 * t276;
t219 = t249 * t269 - t260 * t292;
t220 = t249 * t292 + t260 * t269;
t151 = Icges(7,5) * t220 + Icges(7,6) * t248 + Icges(7,3) * t219;
t155 = Icges(7,4) * t220 + Icges(7,2) * t248 + Icges(7,6) * t219;
t159 = Icges(7,1) * t220 + Icges(7,4) * t248 + Icges(7,5) * t219;
t82 = t151 * t217 + t155 * t246 + t159 * t218;
t152 = Icges(6,5) * t218 - Icges(6,6) * t217 + Icges(6,3) * t246;
t156 = Icges(6,4) * t218 - Icges(6,2) * t217 + Icges(6,6) * t246;
t160 = Icges(6,1) * t218 - Icges(6,4) * t217 + Icges(6,5) * t246;
t83 = t152 * t246 - t156 * t217 + t160 * t218;
t153 = Icges(6,5) * t220 - Icges(6,6) * t219 + Icges(6,3) * t248;
t157 = Icges(6,4) * t220 - Icges(6,2) * t219 + Icges(6,6) * t248;
t161 = Icges(6,1) * t220 - Icges(6,4) * t219 + Icges(6,5) * t248;
t84 = t153 * t246 - t157 * t217 + t161 * t218;
t345 = (t105 + t106) * t262 + (t82 + t84) * t248 + (t81 + t83) * t246;
t107 = t186 * t219 + t188 * t248 + t190 * t220;
t108 = t187 * t248 - t189 * t219 + t191 * t220;
t85 = t150 * t219 + t154 * t248 + t158 * t220;
t86 = t151 * t219 + t155 * t248 + t159 * t220;
t87 = t152 * t248 - t156 * t219 + t160 * t220;
t88 = t153 * t248 - t157 * t219 + t161 * t220;
t344 = (t107 + t108) * t262 + (t86 + t88) * t248 + (t85 + t87) * t246;
t25 = -t105 * t322 + t258 * t81 + t260 * t82;
t26 = -t106 * t322 + t258 * t83 + t260 * t84;
t343 = t25 + t26;
t27 = -t107 * t322 + t258 * t85 + t260 * t86;
t28 = -t108 * t322 + t258 * t87 + t260 * t88;
t342 = t27 + t28;
t29 = t105 * t273 + (t270 * t82 - t272 * t81) * t271;
t30 = t106 * t273 + (t270 * t84 - t272 * t83) * t271;
t341 = t29 + t30;
t31 = t107 * t273 + (t270 * t86 - t272 * t85) * t271;
t32 = t108 * t273 + (t270 * t88 - t272 * t87) * t271;
t340 = t31 + t32;
t116 = t186 * t244 + t188 * t262 + t190 * t245;
t117 = t187 * t262 - t189 * t244 + t191 * t245;
t95 = t150 * t244 + t154 * t262 + t158 * t245;
t96 = t151 * t244 + t155 * t262 + t159 * t245;
t97 = t152 * t262 - t156 * t244 + t160 * t245;
t98 = t153 * t262 - t157 * t244 + t161 * t245;
t339 = (t116 + t117) * t262 + (t96 + t98) * t248 + (t95 + t97) * t246;
t49 = -t116 * t322 + t258 * t95 + t260 * t96;
t50 = -t117 * t322 + t258 * t97 + t260 * t98;
t338 = t49 + t50;
t53 = t116 * t273 + (t270 * t96 - t272 * t95) * t271;
t54 = t117 * t273 + (t270 * t98 - t272 * t97) * t271;
t337 = t53 + t54;
t314 = rSges(7,2) * t246 + t346 * t217 + t347 * t218;
t313 = rSges(7,2) * t248 + t346 * t219 + t347 * t220;
t310 = rSges(7,2) * t262 + t346 * t244 + t347 * t245;
t336 = t246 / 0.2e1;
t335 = t248 / 0.2e1;
t334 = t258 / 0.2e1;
t333 = t260 / 0.2e1;
t332 = t262 / 0.2e1;
t331 = t273 / 0.2e1;
t277 = cos(qJ(4));
t329 = pkin(4) * t277;
t274 = sin(qJ(4));
t327 = t258 * t274;
t326 = t260 * t274;
t325 = t270 * t271;
t324 = t271 * t272;
t318 = t314 * t248;
t148 = pkin(4) * t327 + pkin(10) * t246 + t329 * t247;
t215 = t247 * pkin(3) + t246 * pkin(9);
t206 = t260 * t215;
t317 = t260 * t148 + t206;
t316 = t313 * t262;
t149 = pkin(4) * t326 + pkin(10) * t248 + t329 * t249;
t165 = rSges(6,1) * t220 - rSges(6,2) * t219 + rSges(6,3) * t248;
t315 = -t149 - t165;
t223 = -t249 * t274 + t260 * t277;
t224 = t249 * t277 + t326;
t176 = rSges(5,1) * t224 + rSges(5,2) * t223 + rSges(5,3) * t248;
t216 = t249 * pkin(3) + t248 * pkin(9);
t312 = -t176 - t216;
t311 = t310 * t246;
t193 = rSges(6,1) * t245 - rSges(6,2) * t244 + rSges(6,3) * t262;
t303 = t274 * t322;
t194 = -pkin(4) * t303 + pkin(10) * t262 + t329 * t263;
t309 = -t193 - t194;
t250 = -t263 * t274 - t277 * t322;
t251 = t263 * t277 - t303;
t207 = rSges(5,1) * t251 + rSges(5,2) * t250 + rSges(5,3) * t262;
t243 = t263 * pkin(3) + t262 * pkin(9);
t308 = -t207 - t243;
t307 = t215 * t322 + t258 * t243;
t242 = pkin(2) * t261 + pkin(8) * t260;
t240 = t273 * t242;
t306 = t273 * t216 + t240;
t241 = pkin(2) * t259 + pkin(8) * t258;
t305 = -t215 - t241;
t304 = t241 * t325 + t242 * t324;
t302 = t273 * t149 + t306;
t301 = -t148 + t305;
t300 = -t149 - t313;
t299 = -t216 + t315;
t298 = -t194 - t310;
t297 = -t243 + t309;
t293 = -t322 / 0.2e1;
t237 = rSges(4,1) * t263 - rSges(4,2) * t262 - rSges(4,3) * t322;
t264 = (pkin(2) * t276 - pkin(8) * t278) * t271;
t291 = (-t237 - t264) * t271;
t290 = t148 * t322 + t258 * t194 + t307;
t289 = -t216 + t300;
t288 = -t243 + t298;
t287 = t215 * t325 + t216 * t324 + t304;
t286 = (-t264 + t308) * t271;
t285 = t246 * t345 + t344 * t248 + t339 * t262;
t284 = (-t264 + t297) * t271;
t283 = t148 * t325 + t149 * t324 + t287;
t282 = (-t264 + t288) * t271;
t281 = t339 * t293 + t338 * t332 + t344 * t333 + t334 * t345 + t342 * t335 + t343 * t336;
t280 = t341 * t336 + t340 * t335 + t337 * t332 + t339 * t331 + t344 * t325 / 0.2e1 - t345 * t324 / 0.2e1;
t257 = rSges(3,3) * t273 + (rSges(3,1) * t276 + rSges(3,2) * t278) * t271;
t256 = Icges(3,5) * t273 + (Icges(3,1) * t276 + Icges(3,4) * t278) * t271;
t255 = Icges(3,6) * t273 + (Icges(3,4) * t276 + Icges(3,2) * t278) * t271;
t254 = Icges(3,3) * t273 + (Icges(3,5) * t276 + Icges(3,6) * t278) * t271;
t236 = Icges(4,1) * t263 - Icges(4,4) * t262 - Icges(4,5) * t322;
t235 = Icges(4,4) * t263 - Icges(4,2) * t262 - Icges(4,6) * t322;
t234 = Icges(4,5) * t263 - Icges(4,6) * t262 - Icges(4,3) * t322;
t233 = rSges(3,1) * t261 - rSges(3,2) * t260 + rSges(3,3) * t325;
t232 = rSges(3,1) * t259 - rSges(3,2) * t258 - rSges(3,3) * t324;
t231 = Icges(3,1) * t261 - Icges(3,4) * t260 + Icges(3,5) * t325;
t230 = Icges(3,1) * t259 - Icges(3,4) * t258 - Icges(3,5) * t324;
t229 = Icges(3,4) * t261 - Icges(3,2) * t260 + Icges(3,6) * t325;
t228 = Icges(3,4) * t259 - Icges(3,2) * t258 - Icges(3,6) * t324;
t227 = Icges(3,5) * t261 - Icges(3,6) * t260 + Icges(3,3) * t325;
t226 = Icges(3,5) * t259 - Icges(3,6) * t258 - Icges(3,3) * t324;
t222 = t247 * t277 + t327;
t221 = -t247 * t274 + t258 * t277;
t209 = -t232 * t273 - t257 * t324;
t208 = t233 * t273 - t257 * t325;
t205 = Icges(5,1) * t251 + Icges(5,4) * t250 + Icges(5,5) * t262;
t204 = Icges(5,4) * t251 + Icges(5,2) * t250 + Icges(5,6) * t262;
t203 = Icges(5,5) * t251 + Icges(5,6) * t250 + Icges(5,3) * t262;
t202 = rSges(4,1) * t249 - rSges(4,2) * t248 + rSges(4,3) * t260;
t201 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t258;
t200 = Icges(4,1) * t249 - Icges(4,4) * t248 + Icges(4,5) * t260;
t199 = Icges(4,1) * t247 - Icges(4,4) * t246 + Icges(4,5) * t258;
t198 = Icges(4,4) * t249 - Icges(4,2) * t248 + Icges(4,6) * t260;
t197 = Icges(4,4) * t247 - Icges(4,2) * t246 + Icges(4,6) * t258;
t196 = Icges(4,5) * t249 - Icges(4,6) * t248 + Icges(4,3) * t260;
t195 = Icges(4,5) * t247 - Icges(4,6) * t246 + Icges(4,3) * t258;
t184 = (t232 * t270 + t233 * t272) * t271;
t182 = t246 * t194;
t181 = t246 * t193;
t175 = rSges(5,1) * t222 + rSges(5,2) * t221 + rSges(5,3) * t246;
t174 = Icges(5,1) * t224 + Icges(5,4) * t223 + Icges(5,5) * t248;
t173 = Icges(5,1) * t222 + Icges(5,4) * t221 + Icges(5,5) * t246;
t172 = Icges(5,4) * t224 + Icges(5,2) * t223 + Icges(5,6) * t248;
t171 = Icges(5,4) * t222 + Icges(5,2) * t221 + Icges(5,6) * t246;
t170 = Icges(5,5) * t224 + Icges(5,6) * t223 + Icges(5,3) * t248;
t169 = Icges(5,5) * t222 + Icges(5,6) * t221 + Icges(5,3) * t246;
t167 = -t202 * t322 - t237 * t260;
t166 = t201 * t322 + t237 * t258;
t163 = rSges(6,1) * t218 - rSges(6,2) * t217 + rSges(6,3) * t246;
t143 = t262 * t165;
t141 = t262 * t149;
t139 = -t234 * t322 - t235 * t262 + t236 * t263;
t138 = t248 * t163;
t136 = t248 * t148;
t135 = t201 * t260 - t202 * t258;
t134 = (-t201 - t241) * t273 + t272 * t291;
t133 = t202 * t273 + t270 * t291 + t240;
t132 = t234 * t260 - t235 * t248 + t236 * t249;
t131 = t234 * t258 - t235 * t246 + t236 * t247;
t130 = (t201 * t270 + t202 * t272) * t271 + t304;
t129 = t176 * t262 - t207 * t248;
t128 = -t175 * t262 + t207 * t246;
t127 = -t193 * t248 + t143;
t126 = -t163 * t262 + t181;
t125 = -t196 * t322 - t198 * t262 + t200 * t263;
t124 = -t195 * t322 - t197 * t262 + t199 * t263;
t123 = t203 * t262 + t204 * t250 + t205 * t251;
t122 = t196 * t260 - t198 * t248 + t200 * t249;
t121 = t195 * t260 - t197 * t248 + t199 * t249;
t120 = t196 * t258 - t198 * t246 + t200 * t247;
t119 = t195 * t258 - t197 * t246 + t199 * t247;
t118 = t175 * t248 - t176 * t246;
t115 = -t165 * t246 + t138;
t114 = t308 * t260 + t312 * t322;
t113 = t175 * t322 + t207 * t258 + t307;
t112 = t203 * t248 + t204 * t223 + t205 * t224;
t111 = t203 * t246 + t204 * t221 + t205 * t222;
t110 = (-t175 + t305) * t273 + t272 * t286;
t109 = t176 * t273 + t270 * t286 + t306;
t104 = t175 * t260 + t312 * t258 + t206;
t103 = t170 * t262 + t172 * t250 + t174 * t251;
t102 = t169 * t262 + t171 * t250 + t173 * t251;
t101 = -t310 * t248 + t316;
t100 = -t314 * t262 + t311;
t99 = (t175 * t270 + t176 * t272) * t271 + t287;
t94 = t170 * t248 + t172 * t223 + t174 * t224;
t93 = t169 * t248 + t171 * t223 + t173 * t224;
t92 = t170 * t246 + t172 * t221 + t174 * t222;
t91 = t169 * t246 + t171 * t221 + t173 * t222;
t90 = t309 * t248 + t141 + t143;
t89 = t181 + t182 + (-t148 - t163) * t262;
t80 = t297 * t260 + t299 * t322;
t79 = t163 * t322 + t193 * t258 + t290;
t78 = -t313 * t246 + t318;
t77 = (-t163 + t301) * t273 + t272 * t284;
t76 = t165 * t273 + t270 * t284 + t302;
t75 = t315 * t246 + t136 + t138;
t74 = t139 * t273 + (-t124 * t272 + t125 * t270) * t271;
t73 = t124 * t258 + t125 * t260 - t139 * t322;
t72 = t298 * t248 + t141 + t316;
t71 = t182 + (-t148 - t314) * t262 + t311;
t70 = t163 * t260 + t299 * t258 + t317;
t69 = (t163 * t270 + t165 * t272) * t271 + t283;
t68 = t288 * t260 + t289 * t322;
t67 = t310 * t258 + t314 * t322 + t290;
t66 = (t301 - t314) * t273 + t272 * t282;
t65 = t270 * t282 + t313 * t273 + t302;
t64 = t132 * t273 + (-t121 * t272 + t122 * t270) * t271;
t63 = t131 * t273 + (-t119 * t272 + t120 * t270) * t271;
t62 = t121 * t258 + t122 * t260 - t132 * t322;
t61 = t119 * t258 + t120 * t260 - t131 * t322;
t60 = t300 * t246 + t136 + t318;
t59 = t289 * t258 + t314 * t260 + t317;
t58 = (t314 * t270 + t313 * t272) * t271 + t283;
t57 = t123 * t273 + (-t102 * t272 + t103 * t270) * t271;
t56 = t102 * t258 + t103 * t260 - t123 * t322;
t55 = t102 * t246 + t103 * t248 + t123 * t262;
t38 = t112 * t273 + (t270 * t94 - t272 * t93) * t271;
t37 = t111 * t273 + (t270 * t92 - t272 * t91) * t271;
t36 = -t112 * t322 + t258 * t93 + t260 * t94;
t35 = -t111 * t322 + t258 * t91 + t260 * t92;
t34 = t112 * t262 + t246 * t93 + t248 * t94;
t33 = t111 * t262 + t246 * t91 + t248 * t92;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(3) * t184 + m(4) * t130 + m(5) * t99 + m(6) * t69 + m(7) * t58; m(7) * (t58 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t69 ^ 2 + t76 ^ 2 + t77 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t99 ^ 2) + m(4) * (t130 ^ 2 + t133 ^ 2 + t134 ^ 2) + m(3) * (t184 ^ 2 + t208 ^ 2 + t209 ^ 2) + (t38 + t64 + (t227 * t325 - t229 * t260 + t231 * t261) * t325 + t340) * t325 + (-t37 - t63 + (-t226 * t324 - t228 * t258 + t230 * t259) * t324 + (-t226 * t325 + t227 * t324 + t228 * t260 + t229 * t258 - t230 * t261 - t231 * t259) * t325 - t341) * t324 + ((t254 * t325 - t255 * t260 + t256 * t261) * t325 - (-t254 * t324 - t255 * t258 + t256 * t259) * t324 + t57 + t74 + ((t229 * t278 + t231 * t276) * t270 - (t228 * t278 + t230 * t276) * t272) * t271 ^ 2 + ((-t226 * t272 + t227 * t270 + t255 * t278 + t256 * t276) * t271 + t254 * t273) * t273 + t337) * t273; m(4) * t135 + m(5) * t104 + m(6) * t70 + m(7) * t59; (t49 / 0.2e1 + t50 / 0.2e1 + t56 / 0.2e1 + t73 / 0.2e1) * t273 + (t31 / 0.2e1 + t32 / 0.2e1 + t38 / 0.2e1 + t64 / 0.2e1) * t260 + (t29 / 0.2e1 + t30 / 0.2e1 + t37 / 0.2e1 + t63 / 0.2e1) * t258 + m(7) * (t58 * t59 + t65 * t68 + t66 * t67) + m(6) * (t69 * t70 + t76 * t80 + t77 * t79) + m(5) * (t104 * t99 + t109 * t114 + t110 * t113) + m(4) * (t130 * t135 + t133 * t167 + t134 * t166) + ((-t53 / 0.2e1 - t54 / 0.2e1 - t57 / 0.2e1 - t74 / 0.2e1) * t278 + (-t25 / 0.2e1 - t26 / 0.2e1 - t35 / 0.2e1 - t61 / 0.2e1) * t272 + (t27 / 0.2e1 + t28 / 0.2e1 + t36 / 0.2e1 + t62 / 0.2e1) * t270) * t271; (-t56 - t73 - t338) * t322 + (t36 + t62 + t342) * t260 + (t35 + t61 + t343) * t258 + m(7) * (t59 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t70 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t104 ^ 2 + t113 ^ 2 + t114 ^ 2) + m(4) * (t135 ^ 2 + t166 ^ 2 + t167 ^ 2); m(5) * t118 + m(6) * t75 + m(7) * t60; t280 + t57 * t332 + t55 * t331 + m(7) * (t58 * t60 + t65 * t72 + t66 * t71) + m(6) * (t69 * t75 + t76 * t90 + t77 * t89) + m(5) * (t109 * t129 + t110 * t128 + t118 * t99) + (t270 * t34 / 0.2e1 - t272 * t33 / 0.2e1) * t271 + t38 * t335 + t37 * t336; t281 + t36 * t335 + t35 * t336 + t34 * t333 + t56 * t332 + t33 * t334 + m(7) * (t59 * t60 + t67 * t71 + t68 * t72) + m(6) * (t70 * t75 + t79 * t89 + t80 * t90) + m(5) * (t104 * t118 + t113 * t128 + t114 * t129) + t55 * t293; t246 * t33 + t248 * t34 + t262 * t55 + m(7) * (t60 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(5) * (t118 ^ 2 + t128 ^ 2 + t129 ^ 2) + t285; m(6) * t115 + m(7) * t78; t280 + m(7) * (t100 * t66 + t101 * t65 + t58 * t78) + m(6) * (t115 * t69 + t126 * t77 + t127 * t76); t281 + m(7) * (t100 * t67 + t101 * t68 + t59 * t78) + m(6) * (t115 * t70 + t126 * t79 + t127 * t80); m(7) * (t100 * t71 + t101 * t72 + t60 * t78) + m(6) * (t115 * t75 + t126 * t89 + t127 * t90) + t285; m(7) * (t100 ^ 2 + t101 ^ 2 + t78 ^ 2) + m(6) * (t115 ^ 2 + t126 ^ 2 + t127 ^ 2) + t285; m(7) * t244; m(7) * (t217 * t65 + t219 * t66 + t244 * t58); m(7) * (t217 * t68 + t219 * t67 + t244 * t59); m(7) * (t217 * t72 + t219 * t71 + t244 * t60); m(7) * (t100 * t219 + t101 * t217 + t244 * t78); m(7) * (t217 ^ 2 + t219 ^ 2 + t244 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
