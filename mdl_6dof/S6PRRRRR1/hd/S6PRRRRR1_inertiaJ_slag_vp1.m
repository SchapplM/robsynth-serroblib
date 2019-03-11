% Calculate joint inertia matrix for
% S6PRRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR1_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:37:36
% EndTime: 2019-03-09 00:37:50
% DurationCPUTime: 6.45s
% Computational Cost: add. (42423->584), mult. (56234->837), div. (0->0), fcn. (70380->14), ass. (0->290)
t278 = sin(pkin(12));
t280 = cos(pkin(12));
t284 = sin(qJ(2));
t281 = cos(pkin(6));
t287 = cos(qJ(2));
t337 = t281 * t287;
t259 = t278 * t284 - t280 * t337;
t261 = t278 * t337 + t280 * t284;
t279 = sin(pkin(6));
t340 = t279 * t287;
t338 = t281 * t284;
t260 = t278 * t287 + t280 * t338;
t277 = qJ(3) + qJ(4);
t317 = qJ(5) + t277;
t271 = sin(t317);
t300 = cos(t317);
t344 = t279 * t280;
t238 = t260 * t300 - t271 * t344;
t282 = sin(qJ(6));
t285 = cos(qJ(6));
t206 = -t238 * t282 + t259 * t285;
t207 = t238 * t285 + t259 * t282;
t295 = t279 * t300;
t237 = t260 * t271 + t280 * t295;
t136 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t237;
t138 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t237;
t140 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t237;
t69 = t136 * t237 + t138 * t206 + t140 * t207;
t262 = -t278 * t338 + t280 * t287;
t345 = t278 * t279;
t240 = t262 * t300 + t271 * t345;
t208 = -t240 * t282 + t261 * t285;
t209 = t240 * t285 + t261 * t282;
t239 = t262 * t271 - t278 * t295;
t137 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t239;
t139 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t239;
t141 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t239;
t70 = t137 * t237 + t139 * t206 + t141 * t207;
t252 = t281 * t271 + t284 * t295;
t241 = -t252 * t282 - t285 * t340;
t242 = t252 * t285 - t282 * t340;
t342 = t279 * t284;
t251 = t271 * t342 - t281 * t300;
t159 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t251;
t160 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t251;
t161 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t251;
t82 = t159 * t237 + t160 * t206 + t161 * t207;
t11 = t259 * t69 + t261 * t70 - t340 * t82;
t210 = Icges(6,5) * t252 - Icges(6,6) * t251 - Icges(6,3) * t340;
t211 = Icges(6,4) * t252 - Icges(6,2) * t251 - Icges(6,6) * t340;
t212 = Icges(6,1) * t252 - Icges(6,4) * t251 - Icges(6,5) * t340;
t116 = t210 * t259 - t211 * t237 + t212 * t238;
t167 = Icges(6,5) * t238 - Icges(6,6) * t237 + Icges(6,3) * t259;
t169 = Icges(6,4) * t238 - Icges(6,2) * t237 + Icges(6,6) * t259;
t171 = Icges(6,1) * t238 - Icges(6,4) * t237 + Icges(6,5) * t259;
t91 = t167 * t259 - t169 * t237 + t171 * t238;
t168 = Icges(6,5) * t240 - Icges(6,6) * t239 + Icges(6,3) * t261;
t170 = Icges(6,4) * t240 - Icges(6,2) * t239 + Icges(6,6) * t261;
t172 = Icges(6,1) * t240 - Icges(6,4) * t239 + Icges(6,5) * t261;
t92 = t168 * t259 - t170 * t237 + t172 * t238;
t364 = -t116 * t340 + t259 * t91 + t261 * t92 + t11;
t117 = t210 * t261 - t211 * t239 + t212 * t240;
t71 = t136 * t239 + t138 * t208 + t140 * t209;
t72 = t137 * t239 + t139 * t208 + t141 * t209;
t83 = t159 * t239 + t160 * t208 + t161 * t209;
t12 = t259 * t71 + t261 * t72 - t340 * t83;
t93 = t167 * t261 - t169 * t239 + t171 * t240;
t94 = t168 * t261 - t170 * t239 + t172 * t240;
t363 = -t117 * t340 + t259 * t93 + t261 * t94 + t12;
t15 = t82 * t281 + (t278 * t70 - t280 * t69) * t279;
t362 = t15 + t116 * t281 + (t278 * t92 - t280 * t91) * t279;
t16 = t83 * t281 + (t278 * t72 - t280 * t71) * t279;
t361 = t16 + t117 * t281 + (t278 * t94 - t280 * t93) * t279;
t101 = -t167 * t340 - t169 * t251 + t171 * t252;
t102 = -t168 * t340 - t170 * t251 + t172 * t252;
t121 = -t210 * t340 - t211 * t251 + t212 * t252;
t73 = t136 * t251 + t138 * t241 + t140 * t242;
t74 = t137 * t251 + t139 * t241 + t141 * t242;
t88 = t159 * t251 + t160 * t241 + t161 * t242;
t21 = t259 * t73 + t261 * t74 - t340 * t88;
t360 = t101 * t259 + t102 * t261 - t121 * t340 + t21;
t23 = t88 * t281 + (t278 * t74 - t280 * t73) * t279;
t359 = t23 + t121 * t281 + (-t101 * t280 + t102 * t278) * t279;
t142 = rSges(7,1) * t207 + rSges(7,2) * t206 + rSges(7,3) * t237;
t335 = pkin(5) * t238 + pkin(11) * t237 + t142;
t162 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t251;
t358 = pkin(5) * t252 + pkin(11) * t251 + t162;
t357 = t237 / 0.2e1;
t356 = t239 / 0.2e1;
t355 = t251 / 0.2e1;
t354 = t259 / 0.2e1;
t353 = t261 / 0.2e1;
t352 = t278 / 0.2e1;
t351 = -t280 / 0.2e1;
t350 = t281 / 0.2e1;
t286 = cos(qJ(3));
t348 = t286 * pkin(3);
t283 = sin(qJ(3));
t343 = t279 * t283;
t341 = t279 * t286;
t339 = t281 * t283;
t336 = t335 * t261;
t143 = rSges(7,1) * t209 + rSges(7,2) * t208 + rSges(7,3) * t239;
t334 = pkin(5) * t240 + pkin(11) * t239 + t143;
t273 = sin(t277);
t305 = pkin(4) * t273;
t274 = cos(t277);
t322 = pkin(4) * t274;
t151 = pkin(10) * t259 + t260 * t322 - t305 * t344;
t146 = t261 * t151;
t173 = rSges(6,1) * t238 - rSges(6,2) * t237 + rSges(6,3) * t259;
t154 = t261 * t173;
t333 = t146 + t154;
t202 = t305 * t281 + (-pkin(10) * t287 + t284 * t322) * t279;
t332 = t151 * t340 + t259 * t202;
t152 = pkin(10) * t261 + t262 * t322 + t305 * t345;
t174 = rSges(6,1) * t240 - rSges(6,2) * t239 + rSges(6,3) * t261;
t331 = -t152 - t174;
t318 = t280 * t343;
t185 = -pkin(3) * t318 + pkin(9) * t259 + t260 * t348;
t231 = pkin(3) * t339 + (-pkin(9) * t287 + t284 * t348) * t279;
t329 = t185 * t340 + t259 * t231;
t319 = t278 * t343;
t186 = pkin(3) * t319 + pkin(9) * t261 + t262 * t348;
t236 = t262 * pkin(2) + t261 * pkin(8);
t234 = t281 * t236;
t328 = t281 * t186 + t234;
t245 = -t262 * t273 + t274 * t345;
t246 = t262 * t274 + t273 * t345;
t184 = rSges(5,1) * t246 + rSges(5,2) * t245 + rSges(5,3) * t261;
t327 = -t184 - t186;
t235 = t260 * pkin(2) + t259 * pkin(8);
t326 = -t185 - t235;
t213 = rSges(6,1) * t252 - rSges(6,2) * t251 - rSges(6,3) * t340;
t130 = t173 * t340 + t259 * t213;
t325 = -t202 - t213;
t243 = -t260 * t273 - t274 * t344;
t244 = t260 * t274 - t273 * t344;
t183 = rSges(5,1) * t244 + rSges(5,2) * t243 + rSges(5,3) * t259;
t257 = -t273 * t342 + t274 * t281;
t258 = t273 * t281 + t274 * t342;
t217 = rSges(5,1) * t258 + rSges(5,2) * t257 - rSges(5,3) * t340;
t134 = t183 * t340 + t259 * t217;
t324 = -t217 - t231;
t323 = t235 * t345 + t236 * t344;
t177 = Icges(5,5) * t244 + Icges(5,6) * t243 + Icges(5,3) * t259;
t179 = Icges(5,4) * t244 + Icges(5,2) * t243 + Icges(5,6) * t259;
t181 = Icges(5,1) * t244 + Icges(5,4) * t243 + Icges(5,5) * t259;
t111 = -t177 * t340 + t179 * t257 + t181 * t258;
t178 = Icges(5,5) * t246 + Icges(5,6) * t245 + Icges(5,3) * t261;
t180 = Icges(5,4) * t246 + Icges(5,2) * t245 + Icges(5,6) * t261;
t182 = Icges(5,1) * t246 + Icges(5,4) * t245 + Icges(5,5) * t261;
t112 = -t178 * t340 + t180 * t257 + t182 * t258;
t214 = Icges(5,5) * t258 + Icges(5,6) * t257 - Icges(5,3) * t340;
t215 = Icges(5,4) * t258 + Icges(5,2) * t257 - Icges(5,6) * t340;
t216 = Icges(5,1) * t258 + Icges(5,4) * t257 - Icges(5,5) * t340;
t125 = -t214 * t340 + t215 * t257 + t216 * t258;
t53 = t111 * t259 + t112 * t261 - t125 * t340;
t320 = -t53 - t360;
t316 = t146 + t336;
t315 = -t152 - t334;
t314 = t281 * t152 + t328;
t313 = -t151 + t326;
t312 = -t186 + t331;
t311 = -t202 - t358;
t310 = -t231 + t325;
t309 = t345 / 0.2e1;
t308 = -t344 / 0.2e1;
t307 = -t340 / 0.2e1;
t306 = t259 * t364 + t363 * t261;
t263 = t281 * t286 - t283 * t342;
t264 = t284 * t341 + t339;
t230 = rSges(4,1) * t264 + rSges(4,2) * t263 - rSges(4,3) * t340;
t265 = (pkin(2) * t284 - pkin(8) * t287) * t279;
t304 = (-t230 - t265) * t279;
t303 = -t186 + t315;
t302 = -t231 + t311;
t301 = t185 * t345 + t186 * t344 + t323;
t89 = t130 + t332;
t84 = t259 * t358 + t335 * t340;
t299 = (-t265 + t324) * t279;
t18 = t237 * t73 + t239 * t74 + t251 * t88;
t3 = t237 * t69 + t239 * t70 + t251 * t82;
t4 = t237 * t71 + t239 * t72 + t251 * t83;
t298 = t11 * t357 + t12 * t356 + t18 * t307 + t21 * t355 + t3 * t354 + t4 * t353;
t118 = t214 * t259 + t215 * t243 + t216 * t244;
t97 = t177 * t259 + t179 * t243 + t181 * t244;
t98 = t178 * t259 + t180 * t243 + t182 * t244;
t40 = -t118 * t340 + t259 * t97 + t261 * t98;
t100 = t178 * t261 + t180 * t245 + t182 * t246;
t119 = t214 * t261 + t215 * t245 + t216 * t246;
t99 = t177 * t261 + t179 * t245 + t181 * t246;
t41 = t100 * t261 - t119 * t340 + t259 * t99;
t297 = t259 * t40 + t261 * t41 + t306;
t296 = (-t265 + t310) * t279;
t294 = t151 * t345 + t152 * t344 + t301;
t67 = t84 + t332;
t293 = (-t265 + t302) * t279;
t292 = -t340 * t360 + t306;
t291 = t359 * t307 + t308 * t364 + t363 * t309 + t360 * t350 + t361 * t353 + t362 * t354;
t290 = t320 * t340 + t297;
t44 = t118 * t281 + (t278 * t98 - t280 * t97) * t279;
t45 = t119 * t281 + (t100 * t278 - t280 * t99) * t279;
t57 = t125 * t281 + (-t111 * t280 + t112 * t278) * t279;
t289 = t57 * t307 + t40 * t308 + t41 * t309 + t53 * t350 + t45 * t353 + t44 * t354 + t291;
t256 = t281 * rSges(3,3) + (rSges(3,1) * t284 + rSges(3,2) * t287) * t279;
t255 = Icges(3,5) * t281 + (Icges(3,1) * t284 + Icges(3,4) * t287) * t279;
t254 = Icges(3,6) * t281 + (Icges(3,4) * t284 + Icges(3,2) * t287) * t279;
t253 = Icges(3,3) * t281 + (Icges(3,5) * t284 + Icges(3,6) * t287) * t279;
t250 = t262 * t286 + t319;
t249 = -t262 * t283 + t278 * t341;
t248 = t260 * t286 - t318;
t247 = -t260 * t283 - t280 * t341;
t229 = Icges(4,1) * t264 + Icges(4,4) * t263 - Icges(4,5) * t340;
t228 = Icges(4,4) * t264 + Icges(4,2) * t263 - Icges(4,6) * t340;
t227 = Icges(4,5) * t264 + Icges(4,6) * t263 - Icges(4,3) * t340;
t226 = rSges(3,1) * t262 - rSges(3,2) * t261 + rSges(3,3) * t345;
t225 = rSges(3,1) * t260 - rSges(3,2) * t259 - rSges(3,3) * t344;
t224 = Icges(3,1) * t262 - Icges(3,4) * t261 + Icges(3,5) * t345;
t223 = Icges(3,1) * t260 - Icges(3,4) * t259 - Icges(3,5) * t344;
t222 = Icges(3,4) * t262 - Icges(3,2) * t261 + Icges(3,6) * t345;
t221 = Icges(3,4) * t260 - Icges(3,2) * t259 - Icges(3,6) * t344;
t220 = Icges(3,5) * t262 - Icges(3,6) * t261 + Icges(3,3) * t345;
t219 = Icges(3,5) * t260 - Icges(3,6) * t259 - Icges(3,3) * t344;
t198 = -t225 * t281 - t256 * t344;
t197 = t226 * t281 - t256 * t345;
t196 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t261;
t195 = rSges(4,1) * t248 + rSges(4,2) * t247 + rSges(4,3) * t259;
t193 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t261;
t192 = Icges(4,1) * t248 + Icges(4,4) * t247 + Icges(4,5) * t259;
t191 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t261;
t190 = Icges(4,4) * t248 + Icges(4,2) * t247 + Icges(4,6) * t259;
t189 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t261;
t188 = Icges(4,5) * t248 + Icges(4,6) * t247 + Icges(4,3) * t259;
t157 = (t225 * t278 + t226 * t280) * t279;
t156 = t261 * t185;
t155 = t261 * t183;
t145 = -t196 * t340 - t230 * t261;
t144 = t195 * t340 + t230 * t259;
t135 = -t184 * t340 - t261 * t217;
t132 = -t227 * t340 + t228 * t263 + t229 * t264;
t131 = -t174 * t340 - t261 * t213;
t128 = t195 * t261 - t196 * t259;
t127 = (-t195 - t235) * t281 + t280 * t304;
t126 = t281 * t196 + t278 * t304 + t234;
t124 = t227 * t261 + t228 * t249 + t229 * t250;
t123 = t227 * t259 + t228 * t247 + t229 * t248;
t122 = -t184 * t259 + t155;
t120 = -t174 * t259 + t154;
t115 = (t195 * t278 + t196 * t280) * t279 + t323;
t114 = -t189 * t340 + t191 * t263 + t193 * t264;
t113 = -t188 * t340 + t190 * t263 + t192 * t264;
t110 = t261 * t324 + t327 * t340;
t109 = t134 + t329;
t108 = t143 * t251 - t162 * t239;
t107 = -t142 * t251 + t162 * t237;
t106 = t189 * t261 + t191 * t249 + t193 * t250;
t105 = t188 * t261 + t190 * t249 + t192 * t250;
t104 = t189 * t259 + t191 * t247 + t193 * t248;
t103 = t188 * t259 + t190 * t247 + t192 * t248;
t96 = (-t183 + t326) * t281 + t280 * t299;
t95 = t281 * t184 + t278 * t299 + t328;
t90 = t261 * t325 + t331 * t340;
t87 = t142 * t239 - t143 * t237;
t86 = t259 * t327 + t155 + t156;
t85 = -t261 * t358 - t334 * t340;
t81 = (t183 * t278 + t184 * t280) * t279 + t301;
t80 = t259 * t331 + t333;
t79 = -t259 * t334 + t336;
t78 = t261 * t310 + t312 * t340;
t77 = t89 + t329;
t76 = (-t173 + t313) * t281 + t280 * t296;
t75 = t281 * t174 + t278 * t296 + t314;
t68 = t261 * t311 + t315 * t340;
t66 = t259 * t312 + t156 + t333;
t65 = (t173 * t278 + t174 * t280) * t279 + t294;
t64 = t132 * t281 + (-t113 * t280 + t114 * t278) * t279;
t63 = t113 * t259 + t114 * t261 - t132 * t340;
t62 = t261 * t302 + t303 * t340;
t61 = t67 + t329;
t60 = (t313 - t335) * t281 + t280 * t293;
t59 = t278 * t293 + t281 * t334 + t314;
t58 = t259 * t315 + t316;
t55 = t124 * t281 + (-t105 * t280 + t106 * t278) * t279;
t54 = t123 * t281 + (-t103 * t280 + t104 * t278) * t279;
t51 = t105 * t259 + t106 * t261 - t124 * t340;
t50 = t103 * t259 + t104 * t261 - t123 * t340;
t35 = t259 * t303 + t156 + t316;
t32 = (t278 * t335 + t280 * t334) * t279 + t294;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(3) * t157 + m(4) * t115 + m(5) * t81 + m(6) * t65 + m(7) * t32; m(7) * (t32 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(6) * (t65 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(5) * (t81 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(4) * (t115 ^ 2 + t126 ^ 2 + t127 ^ 2) + m(3) * (t157 ^ 2 + t197 ^ 2 + t198 ^ 2) + (t45 + t55 + (t220 * t345 - t222 * t261 + t224 * t262) * t345 + t361) * t345 + (-t44 - t54 + (-t219 * t344 - t221 * t259 + t223 * t260) * t344 + (-t219 * t345 + t220 * t344 + t221 * t261 + t222 * t259 - t223 * t262 - t224 * t260) * t345 - t362) * t344 + (-(-t253 * t344 - t259 * t254 + t260 * t255) * t344 + (t253 * t345 - t261 * t254 + t262 * t255) * t345 + t57 + t64 + ((t222 * t287 + t224 * t284) * t278 - (t221 * t287 + t223 * t284) * t280) * t279 ^ 2 + ((-t219 * t280 + t220 * t278 + t254 * t287 + t255 * t284) * t279 + t281 * t253) * t281 + t359) * t281; m(4) * t128 + m(5) * t86 + m(6) * t66 + m(7) * t35; t289 + (t51 * t352 + t50 * t351 - t287 * t64 / 0.2e1) * t279 + m(7) * (t32 * t35 + t59 * t62 + t60 * t61) + m(6) * (t65 * t66 + t75 * t78 + t76 * t77) + m(5) * (t109 * t96 + t110 * t95 + t81 * t86) + m(4) * (t115 * t128 + t126 * t145 + t127 * t144) + t54 * t354 + t55 * t353 + t63 * t350; t259 * t50 + t261 * t51 + (-t63 + t320) * t340 + m(7) * (t35 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t66 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t86 ^ 2) + m(4) * (t128 ^ 2 + t144 ^ 2 + t145 ^ 2) + t297; m(5) * t122 + m(6) * t80 + m(7) * t58; m(7) * (t32 * t58 + t59 * t68 + t60 * t67) + m(6) * (t65 * t80 + t75 * t90 + t76 * t89) + m(5) * (t122 * t81 + t134 * t96 + t135 * t95) + t289; m(7) * (t35 * t58 + t61 * t67 + t62 * t68) + m(6) * (t66 * t80 + t77 * t89 + t78 * t90) + m(5) * (t109 * t134 + t110 * t135 + t122 * t86) + t290; m(7) * (t58 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t80 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(5) * (t122 ^ 2 + t134 ^ 2 + t135 ^ 2) + t290; m(6) * t120 + m(7) * t79; m(7) * (t32 * t79 + t59 * t85 + t60 * t84) + m(6) * (t120 * t65 + t130 * t76 + t131 * t75) + t291; m(7) * (t35 * t79 + t61 * t84 + t62 * t85) + m(6) * (t120 * t66 + t130 * t77 + t131 * t78) + t292; m(7) * (t58 * t79 + t67 * t84 + t68 * t85) + m(6) * (t120 * t80 + t130 * t89 + t131 * t90) + t292; m(7) * (t79 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(6) * (t120 ^ 2 + t130 ^ 2 + t131 ^ 2) + t292; m(7) * t87; t15 * t357 + t16 * t356 + m(7) * (t107 * t60 + t108 * t59 + t32 * t87) + t18 * t350 + t23 * t355 + (t3 * t351 + t352 * t4) * t279; m(7) * (t107 * t61 + t108 * t62 + t35 * t87) + t298; m(7) * (t107 * t67 + t108 * t68 + t58 * t87) + t298; m(7) * (t107 * t84 + t108 * t85 + t79 * t87) + t298; t239 * t4 + t237 * t3 + t251 * t18 + m(7) * (t107 ^ 2 + t108 ^ 2 + t87 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
