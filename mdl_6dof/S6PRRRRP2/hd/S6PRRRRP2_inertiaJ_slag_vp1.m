% Calculate joint inertia matrix for
% S6PRRRRP2
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:18
% EndTime: 2019-03-09 00:01:31
% DurationCPUTime: 6.58s
% Computational Cost: add. (33114->572), mult. (56426->816), div. (0->0), fcn. (71934->12), ass. (0->265)
t355 = rSges(7,1) + pkin(5);
t354 = rSges(7,3) + qJ(6);
t273 = sin(pkin(11));
t275 = cos(pkin(11));
t281 = cos(qJ(2));
t276 = cos(pkin(6));
t279 = sin(qJ(2));
t325 = t276 * t279;
t262 = t273 * t281 + t275 * t325;
t323 = qJ(3) + qJ(4);
t272 = sin(t323);
t274 = sin(pkin(6));
t296 = cos(t323);
t290 = t274 * t296;
t243 = t262 * t272 + t275 * t290;
t264 = -t273 * t325 + t275 * t281;
t245 = t264 * t272 - t273 * t290;
t329 = t274 * t279;
t257 = t272 * t329 - t276 * t296;
t331 = t274 * t275;
t244 = t262 * t296 - t272 * t331;
t324 = t276 * t281;
t261 = t273 * t279 - t275 * t324;
t277 = sin(qJ(5));
t335 = cos(qJ(5));
t216 = t244 * t277 - t261 * t335;
t217 = t244 * t335 + t261 * t277;
t144 = Icges(7,5) * t217 + Icges(7,6) * t243 + Icges(7,3) * t216;
t148 = Icges(7,4) * t217 + Icges(7,2) * t243 + Icges(7,6) * t216;
t152 = Icges(7,1) * t217 + Icges(7,4) * t243 + Icges(7,5) * t216;
t79 = t144 * t216 + t148 * t243 + t152 * t217;
t332 = t273 * t274;
t246 = t264 * t296 + t272 * t332;
t263 = t273 * t324 + t275 * t279;
t218 = t246 * t277 - t263 * t335;
t219 = t246 * t335 + t263 * t277;
t145 = Icges(7,5) * t219 + Icges(7,6) * t245 + Icges(7,3) * t218;
t149 = Icges(7,4) * t219 + Icges(7,2) * t245 + Icges(7,6) * t218;
t153 = Icges(7,1) * t219 + Icges(7,4) * t245 + Icges(7,5) * t218;
t80 = t145 * t216 + t149 * t243 + t153 * t217;
t258 = t276 * t272 + t279 * t290;
t327 = t274 * t281;
t247 = t258 * t277 + t327 * t335;
t248 = t258 * t335 - t277 * t327;
t176 = Icges(7,5) * t248 + Icges(7,6) * t257 + Icges(7,3) * t247;
t178 = Icges(7,4) * t248 + Icges(7,2) * t257 + Icges(7,6) * t247;
t180 = Icges(7,1) * t248 + Icges(7,4) * t257 + Icges(7,5) * t247;
t99 = t176 * t216 + t178 * t243 + t180 * t217;
t5 = t243 * t79 + t245 * t80 + t257 * t99;
t177 = Icges(6,5) * t248 - Icges(6,6) * t247 + Icges(6,3) * t257;
t179 = Icges(6,4) * t248 - Icges(6,2) * t247 + Icges(6,6) * t257;
t181 = Icges(6,1) * t248 - Icges(6,4) * t247 + Icges(6,5) * t257;
t100 = t177 * t243 - t179 * t216 + t181 * t217;
t146 = Icges(6,5) * t217 - Icges(6,6) * t216 + Icges(6,3) * t243;
t150 = Icges(6,4) * t217 - Icges(6,2) * t216 + Icges(6,6) * t243;
t154 = Icges(6,1) * t217 - Icges(6,4) * t216 + Icges(6,5) * t243;
t81 = t146 * t243 - t150 * t216 + t154 * t217;
t147 = Icges(6,5) * t219 - Icges(6,6) * t218 + Icges(6,3) * t245;
t151 = Icges(6,4) * t219 - Icges(6,2) * t218 + Icges(6,6) * t245;
t155 = Icges(6,1) * t219 - Icges(6,4) * t218 + Icges(6,5) * t245;
t82 = t147 * t243 - t151 * t216 + t155 * t217;
t6 = t100 * t257 + t243 * t81 + t245 * t82;
t353 = t5 + t6;
t101 = t176 * t218 + t178 * t245 + t180 * t219;
t83 = t144 * t218 + t148 * t245 + t152 * t219;
t84 = t145 * t218 + t149 * t245 + t153 * t219;
t7 = t101 * t257 + t243 * t83 + t245 * t84;
t102 = t177 * t245 - t179 * t218 + t181 * t219;
t85 = t146 * t245 - t150 * t218 + t154 * t219;
t86 = t147 * t245 - t151 * t218 + t155 * t219;
t8 = t102 * t257 + t243 * t85 + t245 * t86;
t352 = t7 + t8;
t351 = (-t100 - t99) * t327 + (t80 + t82) * t263 + (t79 + t81) * t261;
t350 = (-t101 - t102) * t327 + (t84 + t86) * t263 + (t83 + t85) * t261;
t107 = t176 * t247 + t178 * t257 + t180 * t248;
t91 = t144 * t247 + t148 * t257 + t152 * t248;
t92 = t145 * t247 + t149 * t257 + t153 * t248;
t35 = t107 * t257 + t243 * t91 + t245 * t92;
t108 = t177 * t257 - t179 * t247 + t181 * t248;
t93 = t146 * t257 - t150 * t247 + t154 * t248;
t94 = t147 * t257 - t151 * t247 + t155 * t248;
t36 = t108 * t257 + t243 * t93 + t245 * t94;
t349 = t36 + t35;
t348 = (-t107 - t108) * t327 + (t92 + t94) * t263 + (t91 + t93) * t261;
t321 = rSges(7,2) * t243 + t354 * t216 + t355 * t217;
t316 = rSges(7,2) * t257 + t354 * t247 + t355 * t248;
t184 = Icges(5,5) * t244 - Icges(5,6) * t243 + Icges(5,3) * t261;
t186 = Icges(5,4) * t244 - Icges(5,2) * t243 + Icges(5,6) * t261;
t188 = Icges(5,1) * t244 - Icges(5,4) * t243 + Icges(5,5) * t261;
t111 = t184 * t261 - t186 * t243 + t188 * t244;
t185 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t263;
t187 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t263;
t189 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t263;
t112 = t185 * t261 - t187 * t243 + t189 * t244;
t220 = Icges(5,5) * t258 - Icges(5,6) * t257 - Icges(5,3) * t327;
t221 = Icges(5,4) * t258 - Icges(5,2) * t257 - Icges(5,6) * t327;
t222 = Icges(5,1) * t258 - Icges(5,4) * t257 - Icges(5,5) * t327;
t128 = t220 * t261 - t221 * t243 + t222 * t244;
t347 = t111 * t261 + t112 * t263 - t128 * t327 + t351;
t113 = t184 * t263 - t186 * t245 + t188 * t246;
t114 = t185 * t263 - t187 * t245 + t189 * t246;
t129 = t220 * t263 - t221 * t245 + t222 * t246;
t346 = t113 * t261 + t114 * t263 - t129 * t327 + t350;
t29 = t99 * t276 + (t273 * t80 - t275 * t79) * t274;
t30 = t100 * t276 + (t273 * t82 - t275 * t81) * t274;
t345 = t29 + t30 + t128 * t276 + (-t111 * t275 + t112 * t273) * t274;
t31 = t101 * t276 + (t273 * t84 - t275 * t83) * t274;
t32 = t102 * t276 + (t273 * t86 - t275 * t85) * t274;
t344 = t31 + t32 + t129 * t276 + (-t113 * t275 + t114 * t273) * t274;
t121 = -t184 * t327 - t186 * t257 + t188 * t258;
t122 = -t185 * t327 - t187 * t257 + t189 * t258;
t133 = -t220 * t327 - t221 * t257 + t222 * t258;
t343 = t121 * t261 + t122 * t263 - t133 * t327 + t348;
t45 = t107 * t276 + (t273 * t92 - t275 * t91) * t274;
t46 = t108 * t276 + (t273 * t94 - t275 * t93) * t274;
t342 = t45 + t46 + t133 * t276 + (-t121 * t275 + t122 * t273) * t274;
t338 = t261 / 0.2e1;
t337 = t263 / 0.2e1;
t336 = t276 / 0.2e1;
t280 = cos(qJ(3));
t334 = pkin(3) * t280;
t278 = sin(qJ(3));
t330 = t274 * t278;
t328 = t274 * t280;
t326 = t276 * t278;
t157 = rSges(6,1) * t217 - rSges(6,2) * t216 + rSges(6,3) * t243;
t210 = pkin(4) * t244 + pkin(10) * t243;
t194 = t263 * t210;
t322 = t263 * t157 + t194;
t320 = rSges(7,2) * t245 + t354 * t218 + t355 * t219;
t159 = rSges(6,1) * t219 - rSges(6,2) * t218 + rSges(6,3) * t245;
t211 = pkin(4) * t246 + pkin(10) * t245;
t319 = -t159 - t211;
t307 = t275 * t330;
t192 = -pkin(3) * t307 + pkin(9) * t261 + t262 * t334;
t237 = pkin(3) * t326 + (-pkin(9) * t281 + t279 * t334) * t274;
t318 = t192 * t327 + t261 * t237;
t308 = t273 * t330;
t193 = pkin(3) * t308 + pkin(9) * t263 + t264 * t334;
t242 = t264 * pkin(2) + t263 * pkin(8);
t240 = t276 * t242;
t317 = t276 * t193 + t240;
t183 = rSges(6,1) * t248 - rSges(6,2) * t247 + rSges(6,3) * t257;
t236 = pkin(4) * t258 + pkin(10) * t257;
t315 = -t183 - t236;
t191 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t263;
t314 = -t191 - t193;
t241 = t262 * pkin(2) + t261 * pkin(8);
t313 = -t192 - t241;
t312 = t210 * t327 + t261 * t236;
t190 = rSges(5,1) * t244 - rSges(5,2) * t243 + rSges(5,3) * t261;
t223 = rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t327;
t140 = t190 * t327 + t261 * t223;
t311 = -t223 - t237;
t310 = t241 * t332 + t242 * t331;
t306 = t263 * t321 + t194;
t305 = -t211 - t320;
t304 = -t193 + t319;
t303 = t276 * t211 + t317;
t302 = -t236 - t316;
t301 = -t237 + t315;
t300 = -t210 + t313;
t297 = -t327 / 0.2e1;
t265 = t276 * t280 - t278 * t329;
t266 = t279 * t328 + t326;
t235 = rSges(4,1) * t266 + rSges(4,2) * t265 - rSges(4,3) * t327;
t267 = (pkin(2) * t279 - pkin(8) * t281) * t274;
t295 = (-t235 - t267) * t274;
t294 = -t193 + t305;
t293 = t192 * t332 + t193 * t331 + t310;
t292 = -t237 + t302;
t104 = t157 * t327 + t261 * t183 + t312;
t291 = (-t267 + t311) * t274;
t289 = t261 * t347 + t263 * t346;
t288 = (-t267 + t301) * t274;
t287 = t210 * t332 + t211 * t331 + t293;
t77 = t261 * t316 + t321 * t327 + t312;
t286 = (-t267 + t292) * t274;
t285 = t351 * t243 / 0.2e1 + t350 * t245 / 0.2e1 + t348 * t257 / 0.2e1 + t353 * t338 + t352 * t337 + t349 * t297;
t284 = -t327 * t343 + t289;
t283 = t345 * t338 + t344 * t337 + t343 * t336 + t346 * t332 / 0.2e1 - t347 * t331 / 0.2e1 + t342 * t297;
t256 = t276 * rSges(3,3) + (rSges(3,1) * t279 + rSges(3,2) * t281) * t274;
t255 = Icges(3,5) * t276 + (Icges(3,1) * t279 + Icges(3,4) * t281) * t274;
t254 = Icges(3,6) * t276 + (Icges(3,4) * t279 + Icges(3,2) * t281) * t274;
t253 = Icges(3,3) * t276 + (Icges(3,5) * t279 + Icges(3,6) * t281) * t274;
t252 = t264 * t280 + t308;
t251 = -t264 * t278 + t273 * t328;
t250 = t262 * t280 - t307;
t249 = -t262 * t278 - t275 * t328;
t234 = Icges(4,1) * t266 + Icges(4,4) * t265 - Icges(4,5) * t327;
t233 = Icges(4,4) * t266 + Icges(4,2) * t265 - Icges(4,6) * t327;
t232 = Icges(4,5) * t266 + Icges(4,6) * t265 - Icges(4,3) * t327;
t231 = rSges(3,1) * t264 - rSges(3,2) * t263 + rSges(3,3) * t332;
t230 = rSges(3,1) * t262 - rSges(3,2) * t261 - rSges(3,3) * t331;
t229 = Icges(3,1) * t264 - Icges(3,4) * t263 + Icges(3,5) * t332;
t228 = Icges(3,1) * t262 - Icges(3,4) * t261 - Icges(3,5) * t331;
t227 = Icges(3,4) * t264 - Icges(3,2) * t263 + Icges(3,6) * t332;
t226 = Icges(3,4) * t262 - Icges(3,2) * t261 - Icges(3,6) * t331;
t225 = Icges(3,5) * t264 - Icges(3,6) * t263 + Icges(3,3) * t332;
t224 = Icges(3,5) * t262 - Icges(3,6) * t261 - Icges(3,3) * t331;
t208 = -t230 * t276 - t256 * t331;
t207 = t231 * t276 - t256 * t332;
t203 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t263;
t202 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t261;
t201 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t263;
t200 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t261;
t199 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t263;
t198 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t261;
t197 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t263;
t196 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t261;
t170 = (t230 * t273 + t231 * t275) * t274;
t169 = t263 * t192;
t168 = t263 * t190;
t161 = -t203 * t327 - t235 * t263;
t160 = t202 * t327 + t235 * t261;
t141 = -t191 * t327 - t263 * t223;
t137 = -t232 * t327 + t233 * t265 + t234 * t266;
t136 = t202 * t263 - t203 * t261;
t135 = (-t202 - t241) * t276 + t275 * t295;
t134 = t276 * t203 + t273 * t295 + t240;
t132 = t232 * t263 + t233 * t251 + t234 * t252;
t131 = t232 * t261 + t233 * t249 + t234 * t250;
t130 = -t191 * t261 + t168;
t127 = (t202 * t273 + t203 * t275) * t274 + t310;
t126 = -t197 * t327 + t199 * t265 + t201 * t266;
t125 = -t196 * t327 + t198 * t265 + t200 * t266;
t124 = t159 * t257 - t183 * t245;
t123 = -t157 * t257 + t183 * t243;
t120 = t263 * t311 + t314 * t327;
t119 = t140 + t318;
t118 = t197 * t263 + t199 * t251 + t201 * t252;
t117 = t196 * t263 + t198 * t251 + t200 * t252;
t116 = t197 * t261 + t199 * t249 + t201 * t250;
t115 = t196 * t261 + t198 * t249 + t200 * t250;
t110 = (-t190 + t313) * t276 + t275 * t291;
t109 = t276 * t191 + t273 * t291 + t317;
t106 = t157 * t245 - t159 * t243;
t105 = t263 * t315 + t319 * t327;
t103 = t261 * t314 + t168 + t169;
t98 = (t190 * t273 + t191 * t275) * t274 + t293;
t97 = t261 * t319 + t322;
t96 = -t245 * t316 + t257 * t320;
t95 = t243 * t316 - t257 * t321;
t90 = t263 * t301 + t304 * t327;
t89 = t104 + t318;
t88 = (-t157 + t300) * t276 + t275 * t288;
t87 = t276 * t159 + t273 * t288 + t303;
t78 = t263 * t302 + t305 * t327;
t76 = -t243 * t320 + t245 * t321;
t75 = t261 * t304 + t169 + t322;
t74 = (t157 * t273 + t159 * t275) * t274 + t287;
t73 = t137 * t276 + (-t125 * t275 + t126 * t273) * t274;
t72 = t263 * t292 + t294 * t327;
t71 = t77 + t318;
t70 = (t300 - t321) * t276 + t275 * t286;
t69 = t273 * t286 + t276 * t320 + t303;
t68 = t125 * t261 + t126 * t263 - t137 * t327;
t67 = t261 * t305 + t306;
t64 = t132 * t276 + (-t117 * t275 + t118 * t273) * t274;
t63 = t131 * t276 + (-t115 * t275 + t116 * t273) * t274;
t60 = t117 * t261 + t118 * t263 - t132 * t327;
t59 = t115 * t261 + t116 * t263 - t131 * t327;
t52 = t261 * t294 + t169 + t306;
t49 = (t273 * t321 + t275 * t320) * t274 + t287;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(3) * t170 + m(4) * t127 + m(5) * t98 + m(6) * t74 + m(7) * t49; m(7) * (t49 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t74 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t98 ^ 2) + m(4) * (t127 ^ 2 + t134 ^ 2 + t135 ^ 2) + m(3) * (t170 ^ 2 + t207 ^ 2 + t208 ^ 2) + (t64 + (t225 * t332 - t227 * t263 + t229 * t264) * t332 + t344) * t332 + (-t63 + (-t224 * t331 - t226 * t261 + t228 * t262) * t331 + (-t224 * t332 + t225 * t331 + t226 * t263 + t227 * t261 - t228 * t264 - t229 * t262) * t332 - t345) * t331 + ((t253 * t332 - t263 * t254 + t264 * t255) * t332 - (-t253 * t331 - t261 * t254 + t262 * t255) * t331 + t73 + ((t227 * t281 + t229 * t279) * t273 - (t226 * t281 + t228 * t279) * t275) * t274 ^ 2 + ((-t224 * t275 + t225 * t273 + t254 * t281 + t255 * t279) * t274 + t276 * t253) * t276 + t342) * t276; m(4) * t136 + m(5) * t103 + m(6) * t75 + m(7) * t52; m(7) * (t49 * t52 + t69 * t72 + t70 * t71) + m(6) * (t74 * t75 + t87 * t90 + t88 * t89) + m(5) * (t103 * t98 + t109 * t120 + t110 * t119) + m(4) * (t127 * t136 + t134 * t161 + t135 * t160) + (-t281 * t73 / 0.2e1 + t273 * t60 / 0.2e1 - t275 * t59 / 0.2e1) * t274 + t283 + t63 * t338 + t64 * t337 + t68 * t336; t261 * t59 + t263 * t60 + (-t68 - t343) * t327 + m(7) * (t52 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(6) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(5) * (t103 ^ 2 + t119 ^ 2 + t120 ^ 2) + m(4) * (t136 ^ 2 + t160 ^ 2 + t161 ^ 2) + t289; m(5) * t130 + m(6) * t97 + m(7) * t67; m(7) * (t49 * t67 + t69 * t78 + t70 * t77) + m(6) * (t104 * t88 + t105 * t87 + t74 * t97) + m(5) * (t109 * t141 + t110 * t140 + t130 * t98) + t283; m(7) * (t52 * t67 + t71 * t77 + t72 * t78) + m(6) * (t104 * t89 + t105 * t90 + t75 * t97) + m(5) * (t103 * t130 + t119 * t140 + t120 * t141) + t284; m(7) * (t67 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(6) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) + m(5) * (t130 ^ 2 + t140 ^ 2 + t141 ^ 2) + t284; m(6) * t106 + m(7) * t76; (t36 / 0.2e1 + t35 / 0.2e1) * t276 + (t46 / 0.2e1 + t45 / 0.2e1) * t257 + (t31 / 0.2e1 + t32 / 0.2e1) * t245 + (t29 / 0.2e1 + t30 / 0.2e1) * t243 + m(7) * (t49 * t76 + t69 * t96 + t70 * t95) + m(6) * (t106 * t74 + t123 * t88 + t124 * t87) + ((-t6 / 0.2e1 - t5 / 0.2e1) * t275 + (t8 / 0.2e1 + t7 / 0.2e1) * t273) * t274; m(7) * (t52 * t76 + t71 * t95 + t72 * t96) + m(6) * (t106 * t75 + t123 * t89 + t124 * t90) + t285; m(7) * (t67 * t76 + t77 * t95 + t78 * t96) + m(6) * (t104 * t123 + t105 * t124 + t106 * t97) + t285; t349 * t257 + t352 * t245 + t353 * t243 + m(7) * (t76 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(6) * (t106 ^ 2 + t123 ^ 2 + t124 ^ 2); m(7) * t247; m(7) * (t216 * t69 + t218 * t70 + t247 * t49); m(7) * (t216 * t72 + t218 * t71 + t247 * t52); m(7) * (t216 * t78 + t218 * t77 + t247 * t67); m(7) * (t216 * t96 + t218 * t95 + t247 * t76); m(7) * (t216 ^ 2 + t218 ^ 2 + t247 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
