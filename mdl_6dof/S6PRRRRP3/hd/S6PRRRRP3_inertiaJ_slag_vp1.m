% Calculate joint inertia matrix for
% S6PRRRRP3
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:25
% EndTime: 2019-03-09 00:06:41
% DurationCPUTime: 7.02s
% Computational Cost: add. (32371->612), mult. (66378->853), div. (0->0), fcn. (84930->12), ass. (0->278)
t352 = rSges(7,3) + qJ(6);
t277 = cos(pkin(6));
t279 = sin(qJ(3));
t280 = sin(qJ(2));
t275 = sin(pkin(6));
t336 = cos(qJ(3));
t300 = t275 * t336;
t261 = t277 * t279 + t280 * t300;
t273 = qJ(4) + qJ(5);
t269 = sin(t273);
t270 = cos(t273);
t282 = cos(qJ(2));
t327 = t275 * t282;
t244 = -t261 * t269 - t270 * t327;
t245 = t261 * t270 - t269 * t327;
t328 = t275 * t279;
t260 = -t277 * t336 + t280 * t328;
t187 = Icges(7,5) * t245 + Icges(7,6) * t244 + Icges(7,3) * t260;
t189 = Icges(7,4) * t245 + Icges(7,2) * t244 + Icges(7,6) * t260;
t191 = Icges(7,1) * t245 + Icges(7,4) * t244 + Icges(7,5) * t260;
t274 = sin(pkin(11));
t276 = cos(pkin(11));
t326 = t277 * t280;
t257 = t274 * t282 + t276 * t326;
t247 = t257 * t336 - t276 * t328;
t325 = t277 * t282;
t256 = t274 * t280 - t276 * t325;
t217 = -t247 * t269 + t256 * t270;
t218 = t247 * t270 + t256 * t269;
t246 = t257 * t279 + t276 * t300;
t105 = t187 * t246 + t189 * t217 + t191 * t218;
t188 = Icges(6,5) * t245 + Icges(6,6) * t244 + Icges(6,3) * t260;
t190 = Icges(6,4) * t245 + Icges(6,2) * t244 + Icges(6,6) * t260;
t192 = Icges(6,1) * t245 + Icges(6,4) * t244 + Icges(6,5) * t260;
t106 = t188 * t246 + t190 * t217 + t192 * t218;
t259 = -t274 * t326 + t276 * t282;
t248 = t259 * t279 - t274 * t300;
t154 = Icges(7,5) * t218 + Icges(7,6) * t217 + Icges(7,3) * t246;
t158 = Icges(7,4) * t218 + Icges(7,2) * t217 + Icges(7,6) * t246;
t162 = Icges(7,1) * t218 + Icges(7,4) * t217 + Icges(7,5) * t246;
t83 = t154 * t246 + t158 * t217 + t162 * t218;
t249 = t259 * t336 + t274 * t328;
t258 = t274 * t325 + t276 * t280;
t219 = -t249 * t269 + t258 * t270;
t220 = t249 * t270 + t258 * t269;
t155 = Icges(7,5) * t220 + Icges(7,6) * t219 + Icges(7,3) * t248;
t159 = Icges(7,4) * t220 + Icges(7,2) * t219 + Icges(7,6) * t248;
t163 = Icges(7,1) * t220 + Icges(7,4) * t219 + Icges(7,5) * t248;
t84 = t155 * t246 + t159 * t217 + t163 * t218;
t156 = Icges(6,5) * t218 + Icges(6,6) * t217 + Icges(6,3) * t246;
t160 = Icges(6,4) * t218 + Icges(6,2) * t217 + Icges(6,6) * t246;
t164 = Icges(6,1) * t218 + Icges(6,4) * t217 + Icges(6,5) * t246;
t85 = t156 * t246 + t160 * t217 + t164 * t218;
t157 = Icges(6,5) * t220 + Icges(6,6) * t219 + Icges(6,3) * t248;
t161 = Icges(6,4) * t220 + Icges(6,2) * t219 + Icges(6,6) * t248;
t165 = Icges(6,1) * t220 + Icges(6,4) * t219 + Icges(6,5) * t248;
t86 = t157 * t246 + t161 * t217 + t165 * t218;
t351 = (t105 + t106) * t260 + (t84 + t86) * t248 + (t83 + t85) * t246;
t107 = t187 * t248 + t189 * t219 + t191 * t220;
t108 = t188 * t248 + t190 * t219 + t192 * t220;
t87 = t154 * t248 + t158 * t219 + t162 * t220;
t88 = t155 * t248 + t159 * t219 + t163 * t220;
t89 = t156 * t248 + t160 * t219 + t164 * t220;
t90 = t157 * t248 + t161 * t219 + t165 * t220;
t350 = (t107 + t108) * t260 + (t88 + t90) * t248 + (t87 + t89) * t246;
t25 = -t105 * t327 + t256 * t83 + t258 * t84;
t26 = -t106 * t327 + t256 * t85 + t258 * t86;
t349 = t25 + t26;
t27 = -t107 * t327 + t256 * t87 + t258 * t88;
t28 = -t108 * t327 + t256 * t89 + t258 * t90;
t348 = t27 + t28;
t29 = t105 * t277 + (t274 * t84 - t276 * t83) * t275;
t30 = t106 * t277 + (t274 * t86 - t276 * t85) * t275;
t347 = t29 + t30;
t31 = t107 * t277 + (t274 * t88 - t276 * t87) * t275;
t32 = t108 * t277 + (t274 * t90 - t276 * t89) * t275;
t346 = t31 + t32;
t100 = t157 * t260 + t161 * t244 + t165 * t245;
t116 = t187 * t260 + t189 * t244 + t191 * t245;
t117 = t188 * t260 + t190 * t244 + t192 * t245;
t97 = t154 * t260 + t158 * t244 + t162 * t245;
t98 = t155 * t260 + t159 * t244 + t163 * t245;
t99 = t156 * t260 + t160 * t244 + t164 * t245;
t345 = (t116 + t117) * t260 + (t100 + t98) * t248 + (t97 + t99) * t246;
t49 = -t116 * t327 + t256 * t97 + t258 * t98;
t50 = t100 * t258 - t117 * t327 + t256 * t99;
t344 = t49 + t50;
t53 = t116 * t277 + (t274 * t98 - t276 * t97) * t275;
t54 = t117 * t277 + (t100 * t274 - t276 * t99) * t275;
t343 = t53 + t54;
t296 = pkin(5) * t269;
t309 = pkin(5) * t270;
t322 = rSges(7,1) * t218 + rSges(7,2) * t217 + t246 * t352 + t247 * t309 + t256 * t296;
t321 = rSges(7,1) * t220 + rSges(7,2) * t219 + t248 * t352 + t249 * t309 + t258 * t296;
t316 = rSges(7,1) * t245 + rSges(7,2) * t244 + t260 * t352 + t261 * t309 - t296 * t327;
t342 = t246 / 0.2e1;
t341 = t248 / 0.2e1;
t340 = t256 / 0.2e1;
t339 = t258 / 0.2e1;
t338 = t260 / 0.2e1;
t337 = t277 / 0.2e1;
t281 = cos(qJ(4));
t334 = t281 * pkin(4);
t278 = sin(qJ(4));
t332 = t256 * t278;
t331 = t258 * t278;
t330 = t274 * t275;
t329 = t275 * t276;
t324 = t322 * t248;
t323 = t321 * t260;
t152 = pkin(4) * t332 + pkin(10) * t246 + t247 * t334;
t215 = t247 * pkin(3) + t246 * pkin(9);
t207 = t258 * t215;
t320 = t258 * t152 + t207;
t153 = pkin(4) * t331 + pkin(10) * t248 + t249 * t334;
t169 = rSges(6,1) * t220 + rSges(6,2) * t219 + rSges(6,3) * t248;
t319 = -t153 - t169;
t318 = t316 * t246;
t223 = -t249 * t278 + t258 * t281;
t224 = t249 * t281 + t331;
t180 = rSges(5,1) * t224 + rSges(5,2) * t223 + rSges(5,3) * t248;
t216 = t249 * pkin(3) + t248 * pkin(9);
t317 = -t180 - t216;
t194 = rSges(6,1) * t245 + rSges(6,2) * t244 + rSges(6,3) * t260;
t307 = t278 * t327;
t195 = -pkin(4) * t307 + pkin(10) * t260 + t261 * t334;
t315 = -t194 - t195;
t250 = -t261 * t278 - t281 * t327;
t251 = t261 * t281 - t307;
t208 = rSges(5,1) * t251 + rSges(5,2) * t250 + rSges(5,3) * t260;
t243 = t261 * pkin(3) + t260 * pkin(9);
t314 = -t208 - t243;
t313 = t215 * t327 + t256 * t243;
t242 = pkin(2) * t259 + pkin(8) * t258;
t240 = t277 * t242;
t312 = t277 * t216 + t240;
t241 = pkin(2) * t257 + pkin(8) * t256;
t311 = -t215 - t241;
t310 = t241 * t330 + t242 * t329;
t306 = -t153 - t321;
t305 = t277 * t153 + t312;
t304 = -t152 + t311;
t303 = -t216 + t319;
t302 = -t195 - t316;
t301 = -t243 + t315;
t297 = -t327 / 0.2e1;
t237 = rSges(4,1) * t261 - rSges(4,2) * t260 - rSges(4,3) * t327;
t262 = (pkin(2) * t280 - pkin(8) * t282) * t275;
t295 = (-t237 - t262) * t275;
t294 = -t216 + t306;
t293 = t152 * t327 + t256 * t195 + t313;
t292 = -t243 + t302;
t291 = t215 * t330 + t216 * t329 + t310;
t290 = (-t262 + t314) * t275;
t289 = t351 * t246 + t350 * t248 + t345 * t260;
t288 = (-t262 + t301) * t275;
t287 = t152 * t330 + t153 * t329 + t291;
t286 = (-t262 + t292) * t275;
t285 = t345 * t297 + t344 * t338 + t350 * t339 + t351 * t340 + t348 * t341 + t349 * t342;
t284 = t347 * t342 + t346 * t341 + t343 * t338 + t345 * t337 + t350 * t330 / 0.2e1 - t351 * t329 / 0.2e1;
t255 = t277 * rSges(3,3) + (rSges(3,1) * t280 + rSges(3,2) * t282) * t275;
t254 = Icges(3,5) * t277 + (Icges(3,1) * t280 + Icges(3,4) * t282) * t275;
t253 = Icges(3,6) * t277 + (Icges(3,4) * t280 + Icges(3,2) * t282) * t275;
t252 = Icges(3,3) * t277 + (Icges(3,5) * t280 + Icges(3,6) * t282) * t275;
t236 = Icges(4,1) * t261 - Icges(4,4) * t260 - Icges(4,5) * t327;
t235 = Icges(4,4) * t261 - Icges(4,2) * t260 - Icges(4,6) * t327;
t234 = Icges(4,5) * t261 - Icges(4,6) * t260 - Icges(4,3) * t327;
t233 = rSges(3,1) * t259 - rSges(3,2) * t258 + rSges(3,3) * t330;
t232 = rSges(3,1) * t257 - rSges(3,2) * t256 - rSges(3,3) * t329;
t231 = Icges(3,1) * t259 - Icges(3,4) * t258 + Icges(3,5) * t330;
t230 = Icges(3,1) * t257 - Icges(3,4) * t256 - Icges(3,5) * t329;
t229 = Icges(3,4) * t259 - Icges(3,2) * t258 + Icges(3,6) * t330;
t228 = Icges(3,4) * t257 - Icges(3,2) * t256 - Icges(3,6) * t329;
t227 = Icges(3,5) * t259 - Icges(3,6) * t258 + Icges(3,3) * t330;
t226 = Icges(3,5) * t257 - Icges(3,6) * t256 - Icges(3,3) * t329;
t222 = t247 * t281 + t332;
t221 = -t247 * t278 + t256 * t281;
t210 = -t232 * t277 - t255 * t329;
t209 = t233 * t277 - t255 * t330;
t206 = Icges(5,1) * t251 + Icges(5,4) * t250 + Icges(5,5) * t260;
t205 = Icges(5,4) * t251 + Icges(5,2) * t250 + Icges(5,6) * t260;
t204 = Icges(5,5) * t251 + Icges(5,6) * t250 + Icges(5,3) * t260;
t203 = rSges(4,1) * t249 - rSges(4,2) * t248 + rSges(4,3) * t258;
t202 = rSges(4,1) * t247 - rSges(4,2) * t246 + rSges(4,3) * t256;
t201 = Icges(4,1) * t249 - Icges(4,4) * t248 + Icges(4,5) * t258;
t200 = Icges(4,1) * t247 - Icges(4,4) * t246 + Icges(4,5) * t256;
t199 = Icges(4,4) * t249 - Icges(4,2) * t248 + Icges(4,6) * t258;
t198 = Icges(4,4) * t247 - Icges(4,2) * t246 + Icges(4,6) * t256;
t197 = Icges(4,5) * t249 - Icges(4,6) * t248 + Icges(4,3) * t258;
t196 = Icges(4,5) * t247 - Icges(4,6) * t246 + Icges(4,3) * t256;
t186 = (t232 * t274 + t233 * t276) * t275;
t184 = t246 * t195;
t183 = t246 * t194;
t179 = rSges(5,1) * t222 + rSges(5,2) * t221 + rSges(5,3) * t246;
t178 = Icges(5,1) * t224 + Icges(5,4) * t223 + Icges(5,5) * t248;
t177 = Icges(5,1) * t222 + Icges(5,4) * t221 + Icges(5,5) * t246;
t176 = Icges(5,4) * t224 + Icges(5,2) * t223 + Icges(5,6) * t248;
t175 = Icges(5,4) * t222 + Icges(5,2) * t221 + Icges(5,6) * t246;
t174 = Icges(5,5) * t224 + Icges(5,6) * t223 + Icges(5,3) * t248;
t173 = Icges(5,5) * t222 + Icges(5,6) * t221 + Icges(5,3) * t246;
t171 = -t203 * t327 - t237 * t258;
t170 = t202 * t327 + t237 * t256;
t167 = rSges(6,1) * t218 + rSges(6,2) * t217 + rSges(6,3) * t246;
t147 = t260 * t169;
t145 = t260 * t153;
t143 = -t234 * t327 - t235 * t260 + t236 * t261;
t142 = t248 * t167;
t140 = t248 * t152;
t137 = t202 * t258 - t203 * t256;
t136 = (-t202 - t241) * t277 + t276 * t295;
t135 = t277 * t203 + t274 * t295 + t240;
t134 = t234 * t258 - t235 * t248 + t236 * t249;
t133 = t234 * t256 - t235 * t246 + t236 * t247;
t130 = (t202 * t274 + t203 * t276) * t275 + t310;
t129 = t180 * t260 - t208 * t248;
t128 = -t179 * t260 + t208 * t246;
t127 = -t194 * t248 + t147;
t126 = -t167 * t260 + t183;
t125 = -t197 * t327 - t199 * t260 + t201 * t261;
t124 = -t196 * t327 - t198 * t260 + t200 * t261;
t123 = t204 * t260 + t205 * t250 + t206 * t251;
t122 = t197 * t258 - t199 * t248 + t201 * t249;
t121 = t196 * t258 - t198 * t248 + t200 * t249;
t120 = t197 * t256 - t199 * t246 + t201 * t247;
t119 = t196 * t256 - t198 * t246 + t200 * t247;
t118 = t179 * t248 - t180 * t246;
t115 = -t169 * t246 + t142;
t114 = t258 * t314 + t317 * t327;
t113 = t179 * t327 + t208 * t256 + t313;
t112 = t204 * t248 + t205 * t223 + t206 * t224;
t111 = t204 * t246 + t205 * t221 + t206 * t222;
t110 = (-t179 + t311) * t277 + t276 * t290;
t109 = t277 * t180 + t274 * t290 + t312;
t104 = t258 * t179 + t256 * t317 + t207;
t103 = t174 * t260 + t176 * t250 + t178 * t251;
t102 = t173 * t260 + t175 * t250 + t177 * t251;
t101 = (t179 * t274 + t180 * t276) * t275 + t291;
t96 = t174 * t248 + t176 * t223 + t178 * t224;
t95 = t173 * t248 + t175 * t223 + t177 * t224;
t94 = t174 * t246 + t176 * t221 + t178 * t222;
t93 = t173 * t246 + t175 * t221 + t177 * t222;
t92 = t248 * t315 + t145 + t147;
t91 = t183 + t184 + (-t152 - t167) * t260;
t82 = -t248 * t316 + t323;
t81 = -t260 * t322 + t318;
t80 = t258 * t301 + t303 * t327;
t79 = t167 * t327 + t194 * t256 + t293;
t78 = (-t167 + t304) * t277 + t276 * t288;
t77 = t277 * t169 + t274 * t288 + t305;
t76 = t246 * t319 + t140 + t142;
t75 = t143 * t277 + (-t124 * t276 + t125 * t274) * t275;
t74 = t124 * t256 + t125 * t258 - t143 * t327;
t73 = -t246 * t321 + t324;
t72 = t258 * t167 + t256 * t303 + t320;
t71 = (t167 * t274 + t169 * t276) * t275 + t287;
t70 = t134 * t277 + (-t121 * t276 + t122 * t274) * t275;
t69 = t133 * t277 + (-t119 * t276 + t120 * t274) * t275;
t68 = t121 * t256 + t122 * t258 - t134 * t327;
t67 = t119 * t256 + t120 * t258 - t133 * t327;
t66 = t248 * t302 + t145 + t323;
t65 = t184 + (-t152 - t322) * t260 + t318;
t64 = t258 * t292 + t294 * t327;
t63 = t256 * t316 + t322 * t327 + t293;
t62 = (t304 - t322) * t277 + t276 * t286;
t61 = t274 * t286 + t277 * t321 + t305;
t60 = t246 * t306 + t140 + t324;
t59 = t256 * t294 + t258 * t322 + t320;
t58 = (t274 * t322 + t276 * t321) * t275 + t287;
t57 = t123 * t277 + (-t102 * t276 + t103 * t274) * t275;
t56 = t102 * t256 + t103 * t258 - t123 * t327;
t55 = t102 * t246 + t103 * t248 + t123 * t260;
t38 = t112 * t277 + (t274 * t96 - t276 * t95) * t275;
t37 = t111 * t277 + (t274 * t94 - t276 * t93) * t275;
t36 = -t112 * t327 + t256 * t95 + t258 * t96;
t35 = -t111 * t327 + t256 * t93 + t258 * t94;
t34 = t112 * t260 + t246 * t95 + t248 * t96;
t33 = t111 * t260 + t246 * t93 + t248 * t94;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(3) * t186 + m(4) * t130 + m(5) * t101 + m(6) * t71 + m(7) * t58; m(7) * (t58 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t71 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t101 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(4) * (t130 ^ 2 + t135 ^ 2 + t136 ^ 2) + m(3) * (t186 ^ 2 + t209 ^ 2 + t210 ^ 2) + (t38 + t70 + (t227 * t330 - t229 * t258 + t231 * t259) * t330 + t346) * t330 + (-t37 - t69 + (-t226 * t329 - t228 * t256 + t230 * t257) * t329 + (-t226 * t330 + t227 * t329 + t228 * t258 + t229 * t256 - t230 * t259 - t231 * t257) * t330 - t347) * t329 + ((t252 * t330 - t258 * t253 + t259 * t254) * t330 - (-t252 * t329 - t256 * t253 + t257 * t254) * t329 + t57 + t75 + ((t229 * t282 + t231 * t280) * t274 - (t228 * t282 + t230 * t280) * t276) * t275 ^ 2 + ((-t226 * t276 + t227 * t274 + t253 * t282 + t254 * t280) * t275 + t277 * t252) * t277 + t343) * t277; m(4) * t137 + m(5) * t104 + m(6) * t72 + m(7) * t59; (t74 / 0.2e1 + t49 / 0.2e1 + t50 / 0.2e1 + t56 / 0.2e1) * t277 + (t70 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 + t38 / 0.2e1) * t258 + (t69 / 0.2e1 + t29 / 0.2e1 + t30 / 0.2e1 + t37 / 0.2e1) * t256 + m(7) * (t58 * t59 + t61 * t64 + t62 * t63) + m(6) * (t71 * t72 + t77 * t80 + t78 * t79) + m(5) * (t101 * t104 + t109 * t114 + t110 * t113) + m(4) * (t130 * t137 + t135 * t171 + t136 * t170) + ((-t53 / 0.2e1 - t54 / 0.2e1 - t57 / 0.2e1 - t75 / 0.2e1) * t282 + (-t25 / 0.2e1 - t26 / 0.2e1 - t35 / 0.2e1 - t67 / 0.2e1) * t276 + (t27 / 0.2e1 + t28 / 0.2e1 + t36 / 0.2e1 + t68 / 0.2e1) * t274) * t275; (-t56 - t74 - t344) * t327 + (t36 + t68 + t348) * t258 + (t35 + t67 + t349) * t256 + m(7) * (t59 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t72 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(5) * (t104 ^ 2 + t113 ^ 2 + t114 ^ 2) + m(4) * (t137 ^ 2 + t170 ^ 2 + t171 ^ 2); m(5) * t118 + m(6) * t76 + m(7) * t60; m(7) * (t58 * t60 + t61 * t66 + t62 * t65) + m(6) * (t71 * t76 + t77 * t92 + t78 * t91) + m(5) * (t101 * t118 + t109 * t129 + t110 * t128) + (t274 * t34 / 0.2e1 - t276 * t33 / 0.2e1) * t275 + t37 * t342 + t38 * t341 + t57 * t338 + t55 * t337 + t284; m(7) * (t59 * t60 + t63 * t65 + t64 * t66) + m(6) * (t72 * t76 + t79 * t91 + t80 * t92) + m(5) * (t104 * t118 + t113 * t128 + t114 * t129) + t55 * t297 + t35 * t342 + t36 * t341 + t33 * t340 + t34 * t339 + t56 * t338 + t285; t246 * t33 + t248 * t34 + t260 * t55 + m(7) * (t60 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t76 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t118 ^ 2 + t128 ^ 2 + t129 ^ 2) + t289; m(6) * t115 + m(7) * t73; m(7) * (t58 * t73 + t61 * t82 + t62 * t81) + m(6) * (t115 * t71 + t126 * t78 + t127 * t77) + t284; m(7) * (t59 * t73 + t63 * t81 + t64 * t82) + m(6) * (t115 * t72 + t126 * t79 + t127 * t80) + t285; m(7) * (t60 * t73 + t65 * t81 + t66 * t82) + m(6) * (t115 * t76 + t126 * t91 + t127 * t92) + t289; m(7) * (t73 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(6) * (t115 ^ 2 + t126 ^ 2 + t127 ^ 2) + t289; m(7) * t260; m(7) * (t246 * t61 + t248 * t62 + t260 * t58); m(7) * (t246 * t64 + t248 * t63 + t260 * t59); m(7) * (t246 * t66 + t248 * t65 + t260 * t60); m(7) * (t246 * t82 + t248 * t81 + t260 * t73); m(7) * (t246 ^ 2 + t248 ^ 2 + t260 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
