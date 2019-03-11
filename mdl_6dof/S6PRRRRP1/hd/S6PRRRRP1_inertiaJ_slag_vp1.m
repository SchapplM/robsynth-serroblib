% Calculate joint inertia matrix for
% S6PRRRRP1
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:48
% EndTime: 2019-03-08 23:56:05
% DurationCPUTime: 6.92s
% Computational Cost: add. (34328->577), mult. (57478->823), div. (0->0), fcn. (73048->12), ass. (0->268)
t359 = rSges(7,3) + qJ(6);
t272 = sin(pkin(11));
t274 = cos(pkin(11));
t282 = cos(qJ(2));
t275 = cos(pkin(6));
t279 = sin(qJ(2));
t327 = t275 * t279;
t260 = t272 * t282 + t274 * t327;
t325 = qJ(3) + qJ(4);
t271 = sin(t325);
t273 = sin(pkin(6));
t297 = cos(t325);
t291 = t273 * t297;
t243 = t260 * t271 + t274 * t291;
t262 = -t272 * t327 + t274 * t282;
t245 = t262 * t271 - t272 * t291;
t331 = t273 * t279;
t257 = t271 * t331 - t275 * t297;
t333 = t273 * t274;
t244 = t260 * t297 - t271 * t333;
t326 = t275 * t282;
t259 = t272 * t279 - t274 * t326;
t277 = sin(qJ(5));
t280 = cos(qJ(5));
t216 = -t244 * t277 + t259 * t280;
t336 = t259 * t277;
t217 = t244 * t280 + t336;
t148 = Icges(7,5) * t217 + Icges(7,6) * t216 + Icges(7,3) * t243;
t152 = Icges(7,4) * t217 + Icges(7,2) * t216 + Icges(7,6) * t243;
t156 = Icges(7,1) * t217 + Icges(7,4) * t216 + Icges(7,5) * t243;
t81 = t148 * t243 + t152 * t216 + t156 * t217;
t334 = t272 * t273;
t246 = t262 * t297 + t271 * t334;
t261 = t272 * t326 + t274 * t279;
t218 = -t246 * t277 + t261 * t280;
t335 = t261 * t277;
t219 = t246 * t280 + t335;
t149 = Icges(7,5) * t219 + Icges(7,6) * t218 + Icges(7,3) * t245;
t153 = Icges(7,4) * t219 + Icges(7,2) * t218 + Icges(7,6) * t245;
t157 = Icges(7,1) * t219 + Icges(7,4) * t218 + Icges(7,5) * t245;
t82 = t149 * t243 + t153 * t216 + t157 * t217;
t258 = t275 * t271 + t279 * t291;
t329 = t273 * t282;
t247 = -t258 * t277 - t280 * t329;
t308 = t277 * t329;
t248 = t258 * t280 - t308;
t178 = Icges(7,5) * t248 + Icges(7,6) * t247 + Icges(7,3) * t257;
t180 = Icges(7,4) * t248 + Icges(7,2) * t247 + Icges(7,6) * t257;
t182 = Icges(7,1) * t248 + Icges(7,4) * t247 + Icges(7,5) * t257;
t99 = t178 * t243 + t180 * t216 + t182 * t217;
t5 = t243 * t81 + t245 * t82 + t257 * t99;
t179 = Icges(6,5) * t248 + Icges(6,6) * t247 + Icges(6,3) * t257;
t181 = Icges(6,4) * t248 + Icges(6,2) * t247 + Icges(6,6) * t257;
t183 = Icges(6,1) * t248 + Icges(6,4) * t247 + Icges(6,5) * t257;
t100 = t179 * t243 + t181 * t216 + t183 * t217;
t150 = Icges(6,5) * t217 + Icges(6,6) * t216 + Icges(6,3) * t243;
t154 = Icges(6,4) * t217 + Icges(6,2) * t216 + Icges(6,6) * t243;
t158 = Icges(6,1) * t217 + Icges(6,4) * t216 + Icges(6,5) * t243;
t83 = t150 * t243 + t154 * t216 + t158 * t217;
t151 = Icges(6,5) * t219 + Icges(6,6) * t218 + Icges(6,3) * t245;
t155 = Icges(6,4) * t219 + Icges(6,2) * t218 + Icges(6,6) * t245;
t159 = Icges(6,1) * t219 + Icges(6,4) * t218 + Icges(6,5) * t245;
t84 = t151 * t243 + t155 * t216 + t159 * t217;
t6 = t100 * t257 + t243 * t83 + t245 * t84;
t358 = t6 + t5;
t101 = t178 * t245 + t180 * t218 + t182 * t219;
t85 = t148 * t245 + t152 * t218 + t156 * t219;
t86 = t149 * t245 + t153 * t218 + t157 * t219;
t7 = t101 * t257 + t243 * t85 + t245 * t86;
t102 = t179 * t245 + t181 * t218 + t183 * t219;
t87 = t150 * t245 + t154 * t218 + t158 * t219;
t88 = t151 * t245 + t155 * t218 + t159 * t219;
t8 = t102 * t257 + t243 * t87 + t245 * t88;
t357 = t8 + t7;
t356 = (-t100 - t99) * t329 + (t82 + t84) * t261 + (t81 + t83) * t259;
t355 = (-t101 - t102) * t329 + (t86 + t88) * t261 + (t85 + t87) * t259;
t107 = t178 * t257 + t180 * t247 + t182 * t248;
t93 = t148 * t257 + t152 * t247 + t156 * t248;
t94 = t149 * t257 + t153 * t247 + t157 * t248;
t35 = t107 * t257 + t243 * t93 + t245 * t94;
t108 = t179 * t257 + t181 * t247 + t183 * t248;
t95 = t150 * t257 + t154 * t247 + t158 * t248;
t96 = t151 * t257 + t155 * t247 + t159 * t248;
t36 = t108 * t257 + t243 * t95 + t245 * t96;
t354 = t36 + t35;
t353 = (-t107 - t108) * t329 + (t94 + t96) * t261 + (t93 + t95) * t259;
t339 = pkin(5) * t280;
t323 = rSges(7,1) * t217 + rSges(7,2) * t216 + pkin(5) * t336 + t243 * t359 + t244 * t339;
t320 = rSges(7,1) * t248 + rSges(7,2) * t247 - pkin(5) * t308 + t257 * t359 + t258 * t339;
t186 = Icges(5,5) * t244 - Icges(5,6) * t243 + Icges(5,3) * t259;
t188 = Icges(5,4) * t244 - Icges(5,2) * t243 + Icges(5,6) * t259;
t190 = Icges(5,1) * t244 - Icges(5,4) * t243 + Icges(5,5) * t259;
t111 = t186 * t259 - t188 * t243 + t190 * t244;
t187 = Icges(5,5) * t246 - Icges(5,6) * t245 + Icges(5,3) * t261;
t189 = Icges(5,4) * t246 - Icges(5,2) * t245 + Icges(5,6) * t261;
t191 = Icges(5,1) * t246 - Icges(5,4) * t245 + Icges(5,5) * t261;
t112 = t187 * t259 - t189 * t243 + t191 * t244;
t220 = Icges(5,5) * t258 - Icges(5,6) * t257 - Icges(5,3) * t329;
t221 = Icges(5,4) * t258 - Icges(5,2) * t257 - Icges(5,6) * t329;
t222 = Icges(5,1) * t258 - Icges(5,4) * t257 - Icges(5,5) * t329;
t128 = t220 * t259 - t221 * t243 + t222 * t244;
t352 = t111 * t259 + t112 * t261 - t128 * t329 + t356;
t113 = t186 * t261 - t188 * t245 + t190 * t246;
t114 = t187 * t261 - t189 * t245 + t191 * t246;
t129 = t220 * t261 - t221 * t245 + t222 * t246;
t351 = t113 * t259 + t114 * t261 - t129 * t329 + t355;
t29 = t99 * t275 + (t272 * t82 - t274 * t81) * t273;
t30 = t100 * t275 + (t272 * t84 - t274 * t83) * t273;
t350 = t29 + t30 + t128 * t275 + (-t111 * t274 + t112 * t272) * t273;
t31 = t101 * t275 + (t272 * t86 - t274 * t85) * t273;
t32 = t102 * t275 + (t272 * t88 - t274 * t87) * t273;
t349 = t31 + t32 + t129 * t275 + (-t113 * t274 + t114 * t272) * t273;
t121 = -t186 * t329 - t188 * t257 + t190 * t258;
t122 = -t187 * t329 - t189 * t257 + t191 * t258;
t133 = -t220 * t329 - t221 * t257 + t222 * t258;
t348 = t121 * t259 + t122 * t261 - t133 * t329 + t353;
t45 = t107 * t275 + (t272 * t94 - t274 * t93) * t273;
t46 = t108 * t275 + (t272 * t96 - t274 * t95) * t273;
t347 = t45 + t46 + t133 * t275 + (-t121 * t274 + t122 * t272) * t273;
t343 = t259 / 0.2e1;
t342 = t261 / 0.2e1;
t341 = t275 / 0.2e1;
t281 = cos(qJ(3));
t340 = pkin(3) * t281;
t278 = sin(qJ(3));
t332 = t273 * t278;
t330 = t273 * t281;
t328 = t275 * t278;
t161 = rSges(6,1) * t217 + rSges(6,2) * t216 + rSges(6,3) * t243;
t211 = pkin(4) * t244 + pkin(10) * t243;
t196 = t261 * t211;
t324 = t261 * t161 + t196;
t322 = rSges(7,1) * t219 + rSges(7,2) * t218 + pkin(5) * t335 + t245 * t359 + t246 * t339;
t163 = rSges(6,1) * t219 + rSges(6,2) * t218 + rSges(6,3) * t245;
t212 = pkin(4) * t246 + pkin(10) * t245;
t321 = -t163 - t212;
t309 = t274 * t332;
t194 = -pkin(3) * t309 + pkin(9) * t259 + t260 * t340;
t237 = pkin(3) * t328 + (-pkin(9) * t282 + t279 * t340) * t273;
t319 = t194 * t329 + t259 * t237;
t310 = t272 * t332;
t195 = pkin(3) * t310 + pkin(9) * t261 + t262 * t340;
t242 = t262 * pkin(2) + t261 * pkin(8);
t240 = t275 * t242;
t318 = t275 * t195 + t240;
t185 = rSges(6,1) * t248 + rSges(6,2) * t247 + rSges(6,3) * t257;
t236 = pkin(4) * t258 + pkin(10) * t257;
t317 = -t185 - t236;
t193 = rSges(5,1) * t246 - rSges(5,2) * t245 + rSges(5,3) * t261;
t316 = -t193 - t195;
t241 = t260 * pkin(2) + t259 * pkin(8);
t315 = -t194 - t241;
t314 = t211 * t329 + t259 * t236;
t192 = rSges(5,1) * t244 - rSges(5,2) * t243 + rSges(5,3) * t259;
t223 = rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t329;
t144 = t192 * t329 + t259 * t223;
t313 = -t223 - t237;
t312 = t241 * t334 + t242 * t333;
t307 = t261 * t323 + t196;
t306 = -t212 - t322;
t305 = -t195 + t321;
t304 = -t236 - t320;
t303 = t275 * t212 + t318;
t302 = -t237 + t317;
t301 = -t211 + t315;
t298 = -t329 / 0.2e1;
t263 = t275 * t281 - t278 * t331;
t264 = t279 * t330 + t328;
t235 = rSges(4,1) * t264 + rSges(4,2) * t263 - rSges(4,3) * t329;
t265 = (pkin(2) * t279 - pkin(8) * t282) * t273;
t296 = (-t235 - t265) * t273;
t295 = -t195 + t306;
t294 = -t237 + t304;
t293 = t194 * t334 + t195 * t333 + t312;
t104 = t161 * t329 + t259 * t185 + t314;
t292 = (-t265 + t313) * t273;
t290 = t259 * t352 + t261 * t351;
t289 = (-t265 + t302) * t273;
t288 = t211 * t334 + t212 * t333 + t293;
t77 = t259 * t320 + t323 * t329 + t314;
t287 = (-t265 + t294) * t273;
t286 = t356 * t243 / 0.2e1 + t355 * t245 / 0.2e1 + t353 * t257 / 0.2e1 + t358 * t343 + t357 * t342 + t354 * t298;
t285 = -t329 * t348 + t290;
t284 = t350 * t343 + t349 * t342 + t348 * t341 + t351 * t334 / 0.2e1 - t352 * t333 / 0.2e1 + t347 * t298;
t256 = t275 * rSges(3,3) + (rSges(3,1) * t279 + rSges(3,2) * t282) * t273;
t255 = Icges(3,5) * t275 + (Icges(3,1) * t279 + Icges(3,4) * t282) * t273;
t254 = Icges(3,6) * t275 + (Icges(3,4) * t279 + Icges(3,2) * t282) * t273;
t253 = Icges(3,3) * t275 + (Icges(3,5) * t279 + Icges(3,6) * t282) * t273;
t252 = t262 * t281 + t310;
t251 = -t262 * t278 + t272 * t330;
t250 = t260 * t281 - t309;
t249 = -t260 * t278 - t274 * t330;
t234 = Icges(4,1) * t264 + Icges(4,4) * t263 - Icges(4,5) * t329;
t233 = Icges(4,4) * t264 + Icges(4,2) * t263 - Icges(4,6) * t329;
t232 = Icges(4,5) * t264 + Icges(4,6) * t263 - Icges(4,3) * t329;
t231 = rSges(3,1) * t262 - rSges(3,2) * t261 + rSges(3,3) * t334;
t230 = rSges(3,1) * t260 - rSges(3,2) * t259 - rSges(3,3) * t333;
t229 = Icges(3,1) * t262 - Icges(3,4) * t261 + Icges(3,5) * t334;
t228 = Icges(3,1) * t260 - Icges(3,4) * t259 - Icges(3,5) * t333;
t227 = Icges(3,4) * t262 - Icges(3,2) * t261 + Icges(3,6) * t334;
t226 = Icges(3,4) * t260 - Icges(3,2) * t259 - Icges(3,6) * t333;
t225 = Icges(3,5) * t262 - Icges(3,6) * t261 + Icges(3,3) * t334;
t224 = Icges(3,5) * t260 - Icges(3,6) * t259 - Icges(3,3) * t333;
t209 = -t230 * t275 - t256 * t333;
t208 = t231 * t275 - t256 * t334;
t204 = rSges(4,1) * t252 + rSges(4,2) * t251 + rSges(4,3) * t261;
t203 = rSges(4,1) * t250 + rSges(4,2) * t249 + rSges(4,3) * t259;
t202 = Icges(4,1) * t252 + Icges(4,4) * t251 + Icges(4,5) * t261;
t201 = Icges(4,1) * t250 + Icges(4,4) * t249 + Icges(4,5) * t259;
t200 = Icges(4,4) * t252 + Icges(4,2) * t251 + Icges(4,6) * t261;
t199 = Icges(4,4) * t250 + Icges(4,2) * t249 + Icges(4,6) * t259;
t198 = Icges(4,5) * t252 + Icges(4,6) * t251 + Icges(4,3) * t261;
t197 = Icges(4,5) * t250 + Icges(4,6) * t249 + Icges(4,3) * t259;
t171 = (t230 * t272 + t231 * t274) * t273;
t170 = t261 * t194;
t169 = t261 * t192;
t165 = -t204 * t329 - t235 * t261;
t164 = t203 * t329 + t235 * t259;
t145 = -t193 * t329 - t223 * t261;
t139 = -t232 * t329 + t233 * t263 + t234 * t264;
t136 = t203 * t261 - t204 * t259;
t135 = (-t203 - t241) * t275 + t274 * t296;
t134 = t204 * t275 + t272 * t296 + t240;
t132 = t232 * t261 + t233 * t251 + t234 * t252;
t131 = t232 * t259 + t233 * t249 + t234 * t250;
t130 = -t193 * t259 + t169;
t127 = (t203 * t272 + t204 * t274) * t273 + t312;
t126 = -t198 * t329 + t200 * t263 + t202 * t264;
t125 = -t197 * t329 + t199 * t263 + t201 * t264;
t124 = t163 * t257 - t185 * t245;
t123 = -t161 * t257 + t185 * t243;
t120 = t261 * t313 + t316 * t329;
t119 = t144 + t319;
t118 = t198 * t261 + t200 * t251 + t202 * t252;
t117 = t197 * t261 + t199 * t251 + t201 * t252;
t116 = t198 * t259 + t200 * t249 + t202 * t250;
t115 = t197 * t259 + t199 * t249 + t201 * t250;
t110 = (-t192 + t315) * t275 + t274 * t292;
t109 = t193 * t275 + t272 * t292 + t318;
t106 = t161 * t245 - t163 * t243;
t105 = t261 * t317 + t321 * t329;
t103 = t259 * t316 + t169 + t170;
t98 = (t192 * t272 + t193 * t274) * t273 + t293;
t97 = t259 * t321 + t324;
t92 = t261 * t302 + t305 * t329;
t91 = t104 + t319;
t90 = (-t161 + t301) * t275 + t274 * t289;
t89 = t163 * t275 + t272 * t289 + t303;
t80 = -t245 * t320 + t257 * t322;
t79 = t243 * t320 - t257 * t323;
t78 = t261 * t304 + t306 * t329;
t76 = t259 * t305 + t170 + t324;
t75 = (t161 * t272 + t163 * t274) * t273 + t288;
t74 = t139 * t275 + (-t125 * t274 + t126 * t272) * t273;
t73 = -t243 * t322 + t245 * t323;
t72 = t125 * t259 + t126 * t261 - t139 * t329;
t70 = t261 * t294 + t295 * t329;
t69 = t77 + t319;
t67 = (t301 - t323) * t275 + t274 * t287;
t66 = t272 * t287 + t275 * t322 + t303;
t65 = t132 * t275 + (-t117 * t274 + t118 * t272) * t273;
t64 = t131 * t275 + (-t115 * t274 + t116 * t272) * t273;
t61 = t259 * t306 + t307;
t60 = t117 * t259 + t118 * t261 - t132 * t329;
t59 = t115 * t259 + t116 * t261 - t131 * t329;
t48 = t259 * t295 + t170 + t307;
t47 = (t272 * t323 + t274 * t322) * t273 + t288;
t1 = [m(2) + m(3) + m(4) + m(5) + m(6) + m(7); m(3) * t171 + m(4) * t127 + m(5) * t98 + m(6) * t75 + m(7) * t47; m(7) * (t47 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(6) * (t75 ^ 2 + t89 ^ 2 + t90 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2 + t98 ^ 2) + m(4) * (t127 ^ 2 + t134 ^ 2 + t135 ^ 2) + m(3) * (t171 ^ 2 + t208 ^ 2 + t209 ^ 2) + (t65 + (t225 * t334 - t227 * t261 + t229 * t262) * t334 + t349) * t334 + (-t64 + (-t224 * t333 - t226 * t259 + t228 * t260) * t333 + (-t224 * t334 + t225 * t333 + t226 * t261 + t227 * t259 - t228 * t262 - t229 * t260) * t334 - t350) * t333 + (-(-t253 * t333 - t259 * t254 + t260 * t255) * t333 + (t253 * t334 - t261 * t254 + t262 * t255) * t334 + t74 + ((t227 * t282 + t229 * t279) * t272 - (t226 * t282 + t228 * t279) * t274) * t273 ^ 2 + ((-t224 * t274 + t225 * t272 + t254 * t282 + t255 * t279) * t273 + t275 * t253) * t275 + t347) * t275; m(4) * t136 + m(5) * t103 + m(6) * t76 + m(7) * t48; t284 + m(7) * (t47 * t48 + t66 * t70 + t67 * t69) + m(6) * (t75 * t76 + t89 * t92 + t90 * t91) + m(5) * (t103 * t98 + t109 * t120 + t110 * t119) + m(4) * (t127 * t136 + t134 * t165 + t135 * t164) + (-t282 * t74 / 0.2e1 - t274 * t59 / 0.2e1 + t272 * t60 / 0.2e1) * t273 + t64 * t343 + t65 * t342 + t72 * t341; t259 * t59 + t261 * t60 + (-t72 - t348) * t329 + m(7) * (t48 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t76 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t103 ^ 2 + t119 ^ 2 + t120 ^ 2) + m(4) * (t136 ^ 2 + t164 ^ 2 + t165 ^ 2) + t290; m(5) * t130 + m(6) * t97 + m(7) * t61; t284 + m(7) * (t47 * t61 + t66 * t78 + t67 * t77) + m(6) * (t104 * t90 + t105 * t89 + t75 * t97) + m(5) * (t109 * t145 + t110 * t144 + t130 * t98); m(7) * (t48 * t61 + t69 * t77 + t70 * t78) + m(6) * (t104 * t91 + t105 * t92 + t76 * t97) + m(5) * (t103 * t130 + t119 * t144 + t120 * t145) + t285; m(7) * (t61 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(6) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) + m(5) * (t130 ^ 2 + t144 ^ 2 + t145 ^ 2) + t285; m(6) * t106 + m(7) * t73; (t36 / 0.2e1 + t35 / 0.2e1) * t275 + (t46 / 0.2e1 + t45 / 0.2e1) * t257 + (t31 / 0.2e1 + t32 / 0.2e1) * t245 + (t30 / 0.2e1 + t29 / 0.2e1) * t243 + m(7) * (t47 * t73 + t66 * t80 + t67 * t79) + m(6) * (t106 * t75 + t123 * t90 + t124 * t89) + ((-t6 / 0.2e1 - t5 / 0.2e1) * t274 + (t7 / 0.2e1 + t8 / 0.2e1) * t272) * t273; m(7) * (t48 * t73 + t69 * t79 + t70 * t80) + m(6) * (t106 * t76 + t123 * t91 + t124 * t92) + t286; m(7) * (t61 * t73 + t77 * t79 + t78 * t80) + m(6) * (t104 * t123 + t105 * t124 + t106 * t97) + t286; t354 * t257 + t357 * t245 + t358 * t243 + m(7) * (t73 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t106 ^ 2 + t123 ^ 2 + t124 ^ 2); m(7) * t257; m(7) * (t243 * t66 + t245 * t67 + t257 * t47); m(7) * (t243 * t70 + t245 * t69 + t257 * t48); m(7) * (t243 * t78 + t245 * t77 + t257 * t61); m(7) * (t243 * t80 + t245 * t79 + t257 * t73); m(7) * (t243 ^ 2 + t245 ^ 2 + t257 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
