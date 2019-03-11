% Calculate joint inertia matrix for
% S6PRRRRR3
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:48:00
% EndTime: 2019-03-09 00:48:09
% DurationCPUTime: 4.99s
% Computational Cost: add. (40468->644), mult. (74508->898), div. (0->0), fcn. (95538->14), ass. (0->300)
t286 = sin(pkin(12));
t288 = cos(pkin(12));
t294 = cos(qJ(2));
t289 = cos(pkin(6));
t292 = sin(qJ(2));
t341 = t289 * t292;
t266 = t286 * t294 + t288 * t341;
t291 = sin(qJ(3));
t287 = sin(pkin(6));
t351 = cos(qJ(3));
t314 = t287 * t351;
t255 = t266 * t291 + t288 * t314;
t357 = t255 / 0.2e1;
t268 = -t286 * t341 + t288 * t294;
t257 = t268 * t291 - t286 * t314;
t356 = t257 / 0.2e1;
t340 = t289 * t294;
t265 = t286 * t292 - t288 * t340;
t355 = t265 / 0.2e1;
t267 = t286 * t340 + t288 * t292;
t354 = t267 / 0.2e1;
t343 = t287 * t291;
t269 = -t289 * t351 + t292 * t343;
t353 = t269 / 0.2e1;
t352 = t289 / 0.2e1;
t293 = cos(qJ(4));
t349 = t293 * pkin(4);
t290 = sin(qJ(4));
t347 = t265 * t290;
t346 = t267 * t290;
t345 = t286 * t287;
t344 = t287 * t288;
t342 = t287 * t294;
t285 = qJ(4) + qJ(5);
t256 = t266 * t351 - t288 * t343;
t280 = sin(t285);
t310 = pkin(5) * t280;
t281 = cos(t285);
t324 = pkin(5) * t281;
t141 = pkin(11) * t255 + t256 * t324 + t265 * t310;
t282 = qJ(6) + t285;
t277 = sin(t282);
t278 = cos(t282);
t220 = -t256 * t277 + t265 * t278;
t221 = t256 * t278 + t265 * t277;
t161 = rSges(7,1) * t221 + rSges(7,2) * t220 + rSges(7,3) * t255;
t143 = t257 * t161;
t339 = t257 * t141 + t143;
t258 = t268 * t351 + t286 * t343;
t142 = pkin(11) * t257 + t258 * t324 + t267 * t310;
t222 = -t258 * t277 + t267 * t278;
t223 = t258 * t278 + t267 * t277;
t162 = rSges(7,1) * t223 + rSges(7,2) * t222 + rSges(7,3) * t257;
t147 = t269 * t162;
t338 = t269 * t142 + t147;
t337 = t141 + t161;
t336 = t142 + t162;
t163 = pkin(4) * t347 + pkin(10) * t255 + t256 * t349;
t218 = t256 * pkin(3) + t255 * pkin(9);
t210 = t267 * t218;
t335 = t267 * t163 + t210;
t164 = pkin(4) * t346 + pkin(10) * t257 + t258 * t349;
t226 = -t258 * t280 + t267 * t281;
t227 = t258 * t281 + t267 * t280;
t172 = rSges(6,1) * t227 + rSges(6,2) * t226 + rSges(6,3) * t257;
t334 = -t164 - t172;
t270 = t289 * t291 + t292 * t314;
t251 = -t270 * t277 - t278 * t342;
t252 = t270 * t278 - t277 * t342;
t193 = rSges(7,1) * t252 + rSges(7,2) * t251 + rSges(7,3) * t269;
t184 = t255 * t193;
t185 = pkin(11) * t269 + t270 * t324 - t310 * t342;
t333 = t255 * t185 + t184;
t230 = -t258 * t290 + t267 * t293;
t231 = t258 * t293 + t346;
t183 = rSges(5,1) * t231 + rSges(5,2) * t230 + rSges(5,3) * t257;
t219 = t258 * pkin(3) + t257 * pkin(9);
t332 = -t183 - t219;
t331 = t185 + t193;
t253 = -t270 * t280 - t281 * t342;
t254 = t270 * t281 - t280 * t342;
t197 = rSges(6,1) * t254 + rSges(6,2) * t253 + rSges(6,3) * t269;
t321 = t290 * t342;
t198 = -pkin(4) * t321 + pkin(10) * t269 + t270 * t349;
t330 = -t197 - t198;
t259 = -t270 * t290 - t293 * t342;
t260 = t270 * t293 - t321;
t211 = rSges(5,1) * t260 + rSges(5,2) * t259 + rSges(5,3) * t269;
t250 = t270 * pkin(3) + t269 * pkin(9);
t329 = -t211 - t250;
t328 = t218 * t342 + t265 * t250;
t249 = pkin(2) * t268 + pkin(8) * t267;
t247 = t289 * t249;
t327 = t289 * t219 + t247;
t248 = pkin(2) * t266 + pkin(8) * t265;
t326 = -t218 - t248;
t325 = t248 * t345 + t249 * t344;
t190 = Icges(7,5) * t252 + Icges(7,6) * t251 + Icges(7,3) * t269;
t191 = Icges(7,4) * t252 + Icges(7,2) * t251 + Icges(7,6) * t269;
t192 = Icges(7,1) * t252 + Icges(7,4) * t251 + Icges(7,5) * t269;
t116 = t190 * t269 + t191 * t251 + t192 * t252;
t155 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t255;
t157 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t255;
t159 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t255;
t97 = t155 * t269 + t157 * t251 + t159 * t252;
t156 = Icges(7,5) * t223 + Icges(7,6) * t222 + Icges(7,3) * t257;
t158 = Icges(7,4) * t223 + Icges(7,2) * t222 + Icges(7,6) * t257;
t160 = Icges(7,1) * t223 + Icges(7,4) * t222 + Icges(7,5) * t257;
t98 = t156 * t269 + t158 * t251 + t160 * t252;
t38 = t116 * t269 + t255 * t97 + t257 * t98;
t105 = t190 * t255 + t191 * t220 + t192 * t221;
t83 = t155 * t255 + t157 * t220 + t159 * t221;
t84 = t156 * t255 + t158 * t220 + t160 * t221;
t7 = t105 * t269 + t255 * t83 + t257 * t84;
t106 = t190 * t257 + t191 * t222 + t192 * t223;
t85 = t155 * t257 + t157 * t222 + t159 * t223;
t86 = t156 * t257 + t158 * t222 + t160 * t223;
t8 = t106 * t269 + t255 * t85 + t257 * t86;
t322 = t255 * t7 + t257 * t8 + t269 * t38;
t320 = -t164 - t336;
t319 = t289 * t164 + t327;
t318 = -t163 + t326;
t317 = -t219 + t334;
t316 = -t198 - t331;
t315 = -t250 + t330;
t313 = t345 / 0.2e1;
t312 = -t344 / 0.2e1;
t311 = -t342 / 0.2e1;
t244 = rSges(4,1) * t270 - rSges(4,2) * t269 - rSges(4,3) * t342;
t271 = (pkin(2) * t292 - pkin(8) * t294) * t287;
t309 = (-t244 - t271) * t287;
t308 = -t219 + t320;
t307 = t163 * t342 + t265 * t198 + t328;
t306 = -t250 + t316;
t305 = t218 * t345 + t219 * t344 + t325;
t304 = (-t271 + t329) * t287;
t15 = -t105 * t342 + t265 * t83 + t267 * t84;
t16 = -t106 * t342 + t265 * t85 + t267 * t86;
t44 = -t116 * t342 + t265 * t97 + t267 * t98;
t303 = t15 * t357 + t16 * t356 + t38 * t311 + t44 * t353 + t8 * t354 + t7 * t355;
t194 = Icges(6,5) * t254 + Icges(6,6) * t253 + Icges(6,3) * t269;
t195 = Icges(6,4) * t254 + Icges(6,2) * t253 + Icges(6,6) * t269;
t196 = Icges(6,1) * t254 + Icges(6,4) * t253 + Icges(6,5) * t269;
t224 = -t256 * t280 + t265 * t281;
t225 = t256 * t281 + t265 * t280;
t107 = t194 * t255 + t195 * t224 + t196 * t225;
t165 = Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t255;
t167 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t255;
t169 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t255;
t87 = t165 * t255 + t167 * t224 + t169 * t225;
t166 = Icges(6,5) * t227 + Icges(6,6) * t226 + Icges(6,3) * t257;
t168 = Icges(6,4) * t227 + Icges(6,2) * t226 + Icges(6,6) * t257;
t170 = Icges(6,1) * t227 + Icges(6,4) * t226 + Icges(6,5) * t257;
t88 = t166 * t255 + t168 * t224 + t170 * t225;
t21 = t107 * t269 + t255 * t87 + t257 * t88;
t108 = t194 * t257 + t195 * t226 + t196 * t227;
t89 = t165 * t257 + t167 * t226 + t169 * t227;
t90 = t166 * t257 + t168 * t226 + t170 * t227;
t22 = t108 * t269 + t255 * t89 + t257 * t90;
t100 = t166 * t269 + t168 * t253 + t170 * t254;
t118 = t194 * t269 + t195 * t253 + t196 * t254;
t99 = t165 * t269 + t167 * t253 + t169 * t254;
t49 = t100 * t257 + t118 * t269 + t255 * t99;
t302 = t255 * t21 + t257 * t22 + t269 * t49 + t322;
t23 = t105 * t289 + (t286 * t84 - t288 * t83) * t287;
t24 = t106 * t289 + (t286 * t86 - t288 * t85) * t287;
t50 = t116 * t289 + (t286 * t98 - t288 * t97) * t287;
t301 = t23 * t357 + t24 * t356 + t7 * t312 + t8 * t313 + t38 * t352 + t50 * t353;
t300 = (-t271 + t315) * t287;
t299 = t163 * t345 + t164 * t344 + t305;
t298 = (-t271 + t306) * t287;
t29 = -t107 * t342 + t265 * t87 + t267 * t88;
t30 = -t108 * t342 + t265 * t89 + t267 * t90;
t52 = t100 * t267 - t118 * t342 + t265 * t99;
t297 = t21 * t355 + t22 * t354 + t29 * t357 + t30 * t356 + t49 * t311 + t52 * t353 + t303;
t31 = t107 * t289 + (t286 * t88 - t288 * t87) * t287;
t32 = t108 * t289 + (t286 * t90 - t288 * t89) * t287;
t54 = t118 * t289 + (t100 * t286 - t288 * t99) * t287;
t296 = t21 * t312 + t22 * t313 + t31 * t357 + t32 * t356 + t49 * t352 + t54 * t353 + t301;
t264 = t289 * rSges(3,3) + (rSges(3,1) * t292 + rSges(3,2) * t294) * t287;
t263 = Icges(3,5) * t289 + (Icges(3,1) * t292 + Icges(3,4) * t294) * t287;
t262 = Icges(3,6) * t289 + (Icges(3,4) * t292 + Icges(3,2) * t294) * t287;
t261 = Icges(3,3) * t289 + (Icges(3,5) * t292 + Icges(3,6) * t294) * t287;
t243 = Icges(4,1) * t270 - Icges(4,4) * t269 - Icges(4,5) * t342;
t242 = Icges(4,4) * t270 - Icges(4,2) * t269 - Icges(4,6) * t342;
t241 = Icges(4,5) * t270 - Icges(4,6) * t269 - Icges(4,3) * t342;
t240 = rSges(3,1) * t268 - rSges(3,2) * t267 + rSges(3,3) * t345;
t239 = rSges(3,1) * t266 - rSges(3,2) * t265 - rSges(3,3) * t344;
t238 = Icges(3,1) * t268 - Icges(3,4) * t267 + Icges(3,5) * t345;
t237 = Icges(3,1) * t266 - Icges(3,4) * t265 - Icges(3,5) * t344;
t236 = Icges(3,4) * t268 - Icges(3,2) * t267 + Icges(3,6) * t345;
t235 = Icges(3,4) * t266 - Icges(3,2) * t265 - Icges(3,6) * t344;
t234 = Icges(3,5) * t268 - Icges(3,6) * t267 + Icges(3,3) * t345;
t233 = Icges(3,5) * t266 - Icges(3,6) * t265 - Icges(3,3) * t344;
t229 = t256 * t293 + t347;
t228 = -t256 * t290 + t265 * t293;
t213 = -t239 * t289 - t264 * t344;
t212 = t240 * t289 - t264 * t345;
t209 = Icges(5,1) * t260 + Icges(5,4) * t259 + Icges(5,5) * t269;
t208 = Icges(5,4) * t260 + Icges(5,2) * t259 + Icges(5,6) * t269;
t207 = Icges(5,5) * t260 + Icges(5,6) * t259 + Icges(5,3) * t269;
t206 = rSges(4,1) * t258 - rSges(4,2) * t257 + rSges(4,3) * t267;
t205 = rSges(4,1) * t256 - rSges(4,2) * t255 + rSges(4,3) * t265;
t204 = Icges(4,1) * t258 - Icges(4,4) * t257 + Icges(4,5) * t267;
t203 = Icges(4,1) * t256 - Icges(4,4) * t255 + Icges(4,5) * t265;
t202 = Icges(4,4) * t258 - Icges(4,2) * t257 + Icges(4,6) * t267;
t201 = Icges(4,4) * t256 - Icges(4,2) * t255 + Icges(4,6) * t265;
t200 = Icges(4,5) * t258 - Icges(4,6) * t257 + Icges(4,3) * t267;
t199 = Icges(4,5) * t256 - Icges(4,6) * t255 + Icges(4,3) * t265;
t189 = (t239 * t286 + t240 * t288) * t287;
t187 = t255 * t198;
t186 = t255 * t197;
t182 = rSges(5,1) * t229 + rSges(5,2) * t228 + rSges(5,3) * t255;
t181 = Icges(5,1) * t231 + Icges(5,4) * t230 + Icges(5,5) * t257;
t180 = Icges(5,1) * t229 + Icges(5,4) * t228 + Icges(5,5) * t255;
t179 = Icges(5,4) * t231 + Icges(5,2) * t230 + Icges(5,6) * t257;
t178 = Icges(5,4) * t229 + Icges(5,2) * t228 + Icges(5,6) * t255;
t177 = Icges(5,5) * t231 + Icges(5,6) * t230 + Icges(5,3) * t257;
t176 = Icges(5,5) * t229 + Icges(5,6) * t228 + Icges(5,3) * t255;
t174 = -t206 * t342 - t244 * t267;
t173 = t205 * t342 + t244 * t265;
t171 = rSges(6,1) * t225 + rSges(6,2) * t224 + rSges(6,3) * t255;
t150 = t269 * t172;
t149 = t269 * t164;
t146 = -t241 * t342 - t242 * t269 + t243 * t270;
t145 = t257 * t171;
t144 = t257 * t163;
t140 = t205 * t267 - t206 * t265;
t139 = (-t205 - t248) * t289 + t288 * t309;
t138 = t206 * t289 + t286 * t309 + t247;
t137 = t241 * t267 - t242 * t257 + t243 * t258;
t136 = t241 * t265 - t242 * t255 + t243 * t256;
t133 = (t205 * t286 + t206 * t288) * t287 + t325;
t132 = t183 * t269 - t211 * t257;
t131 = -t182 * t269 + t211 * t255;
t130 = -t197 * t257 + t150;
t129 = -t171 * t269 + t186;
t128 = -t200 * t342 - t202 * t269 + t204 * t270;
t127 = -t199 * t342 - t201 * t269 + t203 * t270;
t126 = -t193 * t257 + t147;
t125 = -t161 * t269 + t184;
t124 = t207 * t269 + t208 * t259 + t209 * t260;
t123 = t200 * t267 - t202 * t257 + t204 * t258;
t122 = t199 * t267 - t201 * t257 + t203 * t258;
t121 = t200 * t265 - t202 * t255 + t204 * t256;
t120 = t199 * t265 - t201 * t255 + t203 * t256;
t119 = t182 * t257 - t183 * t255;
t117 = -t172 * t255 + t145;
t115 = t267 * t329 + t332 * t342;
t114 = t182 * t342 + t211 * t265 + t328;
t113 = -t162 * t255 + t143;
t112 = t207 * t257 + t208 * t230 + t209 * t231;
t111 = t207 * t255 + t208 * t228 + t209 * t229;
t110 = (-t182 + t326) * t289 + t288 * t304;
t109 = t183 * t289 + t286 * t304 + t327;
t104 = t182 * t267 + t265 * t332 + t210;
t103 = t177 * t269 + t179 * t259 + t181 * t260;
t102 = t176 * t269 + t178 * t259 + t180 * t260;
t101 = (t182 * t286 + t183 * t288) * t287 + t305;
t96 = t177 * t257 + t179 * t230 + t181 * t231;
t95 = t176 * t257 + t178 * t230 + t180 * t231;
t94 = t177 * t255 + t179 * t228 + t181 * t229;
t93 = t176 * t255 + t178 * t228 + t180 * t229;
t92 = t257 * t330 + t149 + t150;
t91 = t186 + t187 + (-t163 - t171) * t269;
t82 = t267 * t315 + t317 * t342;
t81 = t171 * t342 + t197 * t265 + t307;
t80 = -t257 * t331 + t338;
t79 = -t269 * t337 + t333;
t78 = (-t171 + t318) * t289 + t288 * t300;
t77 = t172 * t289 + t286 * t300 + t319;
t76 = t255 * t334 + t144 + t145;
t75 = t146 * t289 + (-t127 * t288 + t128 * t286) * t287;
t74 = t127 * t265 + t128 * t267 - t146 * t342;
t73 = -t255 * t336 + t339;
t72 = t171 * t267 + t265 * t317 + t335;
t71 = (t171 * t286 + t172 * t288) * t287 + t299;
t70 = t137 * t289 + (-t122 * t288 + t123 * t286) * t287;
t69 = t136 * t289 + (-t120 * t288 + t121 * t286) * t287;
t68 = t122 * t265 + t123 * t267 - t137 * t342;
t67 = t120 * t265 + t121 * t267 - t136 * t342;
t66 = t257 * t316 + t149 + t338;
t65 = t187 + (-t163 - t337) * t269 + t333;
t64 = t267 * t306 + t308 * t342;
t63 = t265 * t331 + t337 * t342 + t307;
t62 = (t318 - t337) * t289 + t288 * t298;
t61 = t286 * t298 + t289 * t336 + t319;
t60 = t255 * t320 + t144 + t339;
t59 = t265 * t308 + t267 * t337 + t335;
t58 = (t286 * t337 + t288 * t336) * t287 + t299;
t57 = t124 * t289 + (-t102 * t288 + t103 * t286) * t287;
t56 = t102 * t265 + t103 * t267 - t124 * t342;
t55 = t102 * t255 + t103 * t257 + t124 * t269;
t43 = t112 * t289 + (t286 * t96 - t288 * t95) * t287;
t42 = t111 * t289 + (t286 * t94 - t288 * t93) * t287;
t40 = -t112 * t342 + t265 * t95 + t267 * t96;
t39 = -t111 * t342 + t265 * t93 + t267 * t94;
t34 = t112 * t269 + t255 * t95 + t257 * t96;
t33 = t111 * t269 + t255 * t93 + t257 * t94;
t1 = [m(4) + m(5) + m(6) + m(7) + m(2) + m(3); m(3) * t189 + m(4) * t133 + m(5) * t101 + m(6) * t71 + m(7) * t58; m(7) * (t58 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t71 ^ 2 + t77 ^ 2 + t78 ^ 2) + m(5) * (t101 ^ 2 + t109 ^ 2 + t110 ^ 2) + m(4) * (t133 ^ 2 + t138 ^ 2 + t139 ^ 2) + m(3) * (t189 ^ 2 + t212 ^ 2 + t213 ^ 2) + (t24 + t32 + t43 + t70 + (t234 * t345 - t236 * t267 + t238 * t268) * t345) * t345 + (-t23 - t31 - t42 - t69 + (-t233 * t344 - t235 * t265 + t237 * t266) * t344 + (-t233 * t345 + t234 * t344 + t235 * t267 + t236 * t265 - t237 * t268 - t238 * t266) * t345) * t344 + ((t261 * t345 - t267 * t262 + t268 * t263) * t345 - (-t261 * t344 - t265 * t262 + t266 * t263) * t344 + t50 + t54 + t57 + t75 + ((t236 * t294 + t238 * t292) * t286 - (t235 * t294 + t237 * t292) * t288) * t287 ^ 2 + ((-t233 * t288 + t234 * t286 + t262 * t294 + t263 * t292) * t287 + t261 * t289) * t289) * t289; m(4) * t140 + m(5) * t104 + m(6) * t72 + m(7) * t59; (t52 / 0.2e1 + t56 / 0.2e1 + t74 / 0.2e1 + t44 / 0.2e1) * t289 + (t43 / 0.2e1 + t70 / 0.2e1 + t24 / 0.2e1 + t32 / 0.2e1) * t267 + (t69 / 0.2e1 + t23 / 0.2e1 + t31 / 0.2e1 + t42 / 0.2e1) * t265 + m(7) * (t58 * t59 + t61 * t64 + t62 * t63) + m(6) * (t71 * t72 + t77 * t82 + t78 * t81) + m(5) * (t101 * t104 + t109 * t115 + t110 * t114) + m(4) * (t133 * t140 + t138 * t174 + t139 * t173) + ((-t50 / 0.2e1 - t54 / 0.2e1 - t57 / 0.2e1 - t75 / 0.2e1) * t294 + (-t15 / 0.2e1 - t29 / 0.2e1 - t39 / 0.2e1 - t67 / 0.2e1) * t288 + (t16 / 0.2e1 + t30 / 0.2e1 + t40 / 0.2e1 + t68 / 0.2e1) * t286) * t287; (-t44 - t52 - t56 - t74) * t342 + (t16 + t30 + t40 + t68) * t267 + (t15 + t29 + t39 + t67) * t265 + m(7) * (t59 ^ 2 + t63 ^ 2 + t64 ^ 2) + m(6) * (t72 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t104 ^ 2 + t114 ^ 2 + t115 ^ 2) + m(4) * (t140 ^ 2 + t173 ^ 2 + t174 ^ 2); m(5) * t119 + m(6) * t76 + m(7) * t60; m(7) * (t58 * t60 + t61 * t66 + t62 * t65) + m(6) * (t71 * t76 + t77 * t92 + t78 * t91) + m(5) * (t101 * t119 + t109 * t132 + t110 * t131) + t296 + (t286 * t34 / 0.2e1 - t288 * t33 / 0.2e1) * t287 + t55 * t352 + t57 * t353 + t42 * t357 + t43 * t356; m(7) * (t59 * t60 + t63 * t65 + t64 * t66) + m(6) * (t72 * t76 + t81 * t91 + t82 * t92) + m(5) * (t104 * t119 + t114 * t131 + t115 * t132) + t297 + t55 * t311 + t56 * t353 + t34 * t354 + t39 * t357 + t40 * t356 + t33 * t355; t255 * t33 + t257 * t34 + t269 * t55 + m(7) * (t60 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(6) * (t76 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t119 ^ 2 + t131 ^ 2 + t132 ^ 2) + t302; m(6) * t117 + m(7) * t73; m(7) * (t58 * t73 + t61 * t80 + t62 * t79) + m(6) * (t117 * t71 + t129 * t78 + t130 * t77) + t296; t297 + m(7) * (t59 * t73 + t63 * t79 + t64 * t80) + m(6) * (t117 * t72 + t129 * t81 + t130 * t82); m(7) * (t60 * t73 + t65 * t79 + t66 * t80) + m(6) * (t117 * t76 + t129 * t91 + t130 * t92) + t302; m(7) * (t73 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t117 ^ 2 + t129 ^ 2 + t130 ^ 2) + t302; m(7) * t113; m(7) * (t113 * t58 + t125 * t62 + t126 * t61) + t301; m(7) * (t113 * t59 + t125 * t63 + t126 * t64) + t303; m(7) * (t113 * t60 + t125 * t65 + t126 * t66) + t322; m(7) * (t113 * t73 + t125 * t79 + t126 * t80) + t322; m(7) * (t113 ^ 2 + t125 ^ 2 + t126 ^ 2) + t322;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
