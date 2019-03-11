% Calculate joint inertia matrix for
% S6RRRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:33:27
% EndTime: 2019-03-10 03:33:36
% DurationCPUTime: 4.21s
% Computational Cost: add. (19659->502), mult. (15321->708), div. (0->0), fcn. (16520->12), ass. (0->261)
t241 = qJ(2) + qJ(3);
t231 = qJ(4) + t241;
t224 = cos(t231);
t240 = qJ(5) + qJ(6);
t229 = cos(t240);
t244 = sin(qJ(1));
t315 = t229 * t244;
t227 = sin(t240);
t247 = cos(qJ(1));
t316 = t227 * t247;
t179 = -t224 * t316 + t315;
t314 = t229 * t247;
t317 = t227 * t244;
t180 = t224 * t314 + t317;
t223 = sin(t231);
t319 = t223 * t247;
t109 = t180 * rSges(7,1) + t179 * rSges(7,2) + rSges(7,3) * t319;
t245 = cos(qJ(5));
t225 = pkin(5) * t245 + pkin(4);
t248 = -pkin(11) - pkin(10);
t242 = sin(qJ(5));
t313 = t242 * t244;
t318 = t224 * t247;
t351 = pkin(5) * t313 + t225 * t318 - t248 * t319 + t109;
t238 = t244 ^ 2;
t350 = t244 * pkin(7);
t303 = pkin(4) * t318 + pkin(10) * t319;
t349 = -t303 + t351;
t228 = sin(t241);
t230 = cos(t241);
t283 = rSges(4,1) * t230 - rSges(4,2) * t228;
t272 = Icges(4,5) * t230 - Icges(4,6) * t228;
t171 = -Icges(4,3) * t247 + t244 * t272;
t172 = Icges(4,3) * t244 + t247 * t272;
t239 = t247 ^ 2;
t325 = Icges(4,4) * t230;
t275 = -Icges(4,2) * t228 + t325;
t174 = Icges(4,6) * t244 + t247 * t275;
t326 = Icges(4,4) * t228;
t278 = Icges(4,1) * t230 - t326;
t176 = Icges(4,5) * t244 + t247 * t278;
t267 = -t174 * t228 + t176 * t230;
t173 = -Icges(4,6) * t247 + t244 * t275;
t175 = -Icges(4,5) * t247 + t244 * t278;
t268 = t173 * t228 - t175 * t230;
t177 = -t224 * t317 - t314;
t178 = t224 * t315 - t316;
t320 = t223 * t244;
t102 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t320;
t104 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t320;
t106 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t320;
t30 = t102 * t320 + t104 * t177 + t106 * t178;
t103 = Icges(7,5) * t180 + Icges(7,6) * t179 + Icges(7,3) * t319;
t105 = Icges(7,4) * t180 + Icges(7,2) * t179 + Icges(7,6) * t319;
t107 = Icges(7,1) * t180 + Icges(7,4) * t179 + Icges(7,5) * t319;
t31 = t103 * t320 + t105 * t177 + t107 * t178;
t15 = t244 * t31 - t247 * t30;
t271 = Icges(5,5) * t224 - Icges(5,6) * t223;
t159 = -Icges(5,3) * t247 + t244 * t271;
t160 = Icges(5,3) * t244 + t247 * t271;
t310 = t245 * t247;
t189 = -t224 * t313 - t310;
t311 = t244 * t245;
t312 = t242 * t247;
t190 = t224 * t311 - t312;
t118 = Icges(6,5) * t190 + Icges(6,6) * t189 + Icges(6,3) * t320;
t120 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t320;
t122 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t320;
t44 = t118 * t320 + t120 * t189 + t122 * t190;
t191 = -t224 * t312 + t311;
t192 = t224 * t310 + t313;
t119 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t319;
t121 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t319;
t123 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t319;
t45 = t119 * t320 + t121 * t189 + t123 * t190;
t22 = t244 * t45 - t247 * t44;
t323 = Icges(5,4) * t224;
t274 = -Icges(5,2) * t223 + t323;
t162 = Icges(5,6) * t244 + t247 * t274;
t324 = Icges(5,4) * t223;
t277 = Icges(5,1) * t224 - t324;
t164 = Icges(5,5) * t244 + t247 * t277;
t269 = -t162 * t223 + t164 * t224;
t161 = -Icges(5,6) * t247 + t244 * t274;
t163 = -Icges(5,5) * t247 + t244 * t277;
t270 = t161 * t223 - t163 * t224;
t298 = -t15 - t22 - t239 * t159 - (t269 * t244 + (-t160 + t270) * t247) * t244;
t348 = -t239 * t171 - (t267 * t244 + (-t172 + t268) * t247) * t244 + t298;
t249 = -pkin(8) - pkin(7);
t141 = -Icges(7,3) * t224 + (Icges(7,5) * t229 - Icges(7,6) * t227) * t223;
t142 = -Icges(7,6) * t224 + (Icges(7,4) * t229 - Icges(7,2) * t227) * t223;
t143 = -Icges(7,5) * t224 + (Icges(7,1) * t229 - Icges(7,4) * t227) * t223;
t60 = t141 * t320 + t142 * t177 + t143 * t178;
t5 = -t224 * t60 + (t244 * t30 + t247 * t31) * t223;
t32 = t102 * t319 + t104 * t179 + t106 * t180;
t33 = t103 * t319 + t105 * t179 + t107 * t180;
t61 = t141 * t319 + t142 * t179 + t143 * t180;
t6 = -t224 * t61 + (t244 * t32 + t247 * t33) * t223;
t347 = t6 * t319 + t5 * t320;
t346 = -t224 / 0.2e1;
t345 = t244 / 0.2e1;
t344 = -t247 / 0.2e1;
t243 = sin(qJ(2));
t343 = pkin(2) * t243;
t342 = pkin(3) * t228;
t341 = pkin(4) * t224;
t340 = -pkin(4) + t225;
t339 = pkin(10) + t248;
t246 = cos(qJ(2));
t226 = t246 * pkin(2) + pkin(1);
t338 = rSges(3,1) * t246;
t336 = rSges(3,2) * t243;
t334 = t247 * rSges(3,3);
t41 = -t102 * t224 + (-t104 * t227 + t106 * t229) * t223;
t333 = t41 * t247;
t42 = -t103 * t224 + (-t105 * t227 + t107 * t229) * t223;
t332 = t42 * t244;
t53 = -t118 * t224 + (-t120 * t242 + t122 * t245) * t223;
t331 = t53 * t247;
t54 = -t119 * t224 + (-t121 * t242 + t123 * t245) * t223;
t330 = t54 * t244;
t131 = t223 * t229 * t143;
t322 = t142 * t227;
t68 = -t224 * t141 - t223 * t322 + t131;
t329 = t68 * t224;
t280 = -rSges(7,1) * t178 - rSges(7,2) * t177;
t108 = rSges(7,3) * t320 - t280;
t144 = -rSges(7,3) * t224 + (rSges(7,1) * t229 - rSges(7,2) * t227) * t223;
t76 = t224 * t108 + t144 * t320;
t328 = Icges(3,4) * t243;
t327 = Icges(3,4) * t246;
t146 = -Icges(6,6) * t224 + (Icges(6,4) * t245 - Icges(6,2) * t242) * t223;
t321 = t146 * t242;
t206 = pkin(3) * t230 + t226;
t199 = t247 * t206;
t219 = t247 * t226;
t308 = t247 * (t199 - t219) + (t206 - t226) * t238;
t138 = t223 * t340 + t224 * t339;
t307 = -t138 - t144;
t260 = rSges(5,1) * t318 - rSges(5,2) * t319 + t244 * rSges(5,3);
t282 = rSges(5,1) * t224 - rSges(5,2) * t223;
t110 = t244 * (-rSges(5,3) * t247 + t244 * t282) + t247 * t260;
t150 = -rSges(6,3) * t224 + (rSges(6,1) * t245 - rSges(6,2) * t242) * t223;
t198 = pkin(4) * t223 - pkin(10) * t224;
t306 = -t150 - t198;
t236 = t247 * pkin(7);
t305 = t244 * (t236 + (-pkin(1) + t226) * t244) + t247 * (-t247 * pkin(1) + t219 - t350);
t261 = t244 * rSges(4,3) + t283 * t247;
t117 = t244 * (-t247 * rSges(4,3) + t244 * t283) + t247 * t261;
t304 = t238 * (pkin(10) * t223 + t341) + t247 * t303;
t302 = t244 * rSges(3,3) + t247 * t338;
t300 = t238 + t239;
t16 = t244 * t33 - t247 * t32;
t46 = t118 * t319 + t120 * t191 + t122 * t192;
t47 = t119 * t319 + t121 * t191 + t123 * t192;
t23 = t244 * t47 - t247 * t46;
t299 = (t238 * t160 + t16 + t23 + (t270 * t247 + (-t159 + t269) * t244) * t247) * t244;
t297 = -t198 + t307;
t127 = t192 * rSges(6,1) + t191 * rSges(6,2) + rSges(6,3) * t319;
t296 = t320 / 0.2e1;
t295 = t319 / 0.2e1;
t205 = rSges(4,1) * t228 + rSges(4,2) * t230;
t294 = -t205 - t343;
t197 = rSges(5,1) * t223 + rSges(5,2) * t224;
t293 = -t197 - t342;
t292 = -t198 - t342;
t291 = t244 * (t238 * t172 + (t268 * t247 + (-t171 + t267) * t244) * t247) + t299;
t290 = (t41 + t60) * t296 + (t42 + t61) * t295;
t237 = -pkin(9) + t249;
t289 = -t244 * t237 + t199;
t9 = -t329 + (t244 * t41 + t247 * t42) * t223;
t288 = -t224 * t9 + t347;
t64 = t110 + t308;
t281 = -rSges(6,1) * t190 - rSges(6,2) * t189;
t126 = rSges(6,3) * t320 - t281;
t59 = t244 * t126 + t247 * t127 + t304;
t287 = t15 * t296 + t16 * t295 + t5 * t344 + t6 * t345 + (t332 - t333) * t346;
t286 = -t150 + t292;
t285 = -t342 - t343;
t284 = -t336 + t338;
t279 = Icges(3,1) * t246 - t328;
t276 = -Icges(3,2) * t243 + t327;
t273 = Icges(3,5) * t246 - Icges(3,6) * t243;
t195 = Icges(5,2) * t224 + t324;
t196 = Icges(5,1) * t223 + t323;
t264 = -t195 * t223 + t196 * t224;
t203 = Icges(4,2) * t230 + t326;
t204 = Icges(4,1) * t228 + t325;
t263 = -t203 * t228 + t204 * t230;
t262 = t292 + t307;
t259 = -t197 + t285;
t258 = -t198 + t285;
t115 = -pkin(5) * t312 + (-t223 * t339 + t224 * t340) * t244;
t27 = t304 + t349 * t247 + (t108 + t115) * t244;
t29 = t59 + t308;
t256 = -t150 + t258;
t255 = t247 * t298 + t299;
t254 = t258 + t307;
t25 = t27 + t308;
t145 = -Icges(6,3) * t224 + (Icges(6,5) * t245 - Icges(6,6) * t242) * t223;
t147 = -Icges(6,5) * t224 + (Icges(6,1) * t245 - Icges(6,4) * t242) * t223;
t66 = t145 * t320 + t146 * t189 + t147 * t190;
t10 = -t224 * t66 + (t244 * t44 + t247 * t45) * t223;
t67 = t145 * t319 + t146 * t191 + t147 * t192;
t11 = -t224 * t67 + (t244 * t46 + t247 * t47) * t223;
t253 = t10 * t344 + t11 * t345 + t22 * t296 + t23 * t295 + (t330 - t331) * t346 + t287;
t252 = t247 * t348 + t291;
t194 = Icges(5,5) * t223 + Icges(5,6) * t224;
t251 = -t333 / 0.2e1 + t332 / 0.2e1 - t331 / 0.2e1 + t330 / 0.2e1 + (t162 * t224 + t164 * t223 + t194 * t244 + t247 * t264 + t61 + t67) * t345 + (t161 * t224 + t163 * t223 - t194 * t247 + t244 * t264 + t60 + t66) * t344;
t202 = Icges(4,5) * t228 + Icges(4,6) * t230;
t250 = t251 + (t174 * t230 + t176 * t228 + t202 * t244 + t247 * t263) * t345 + (t173 * t230 + t175 * t228 - t202 * t247 + t244 * t263) * t344;
t218 = rSges(2,1) * t247 - rSges(2,2) * t244;
t217 = -rSges(2,1) * t244 - rSges(2,2) * t247;
t216 = rSges(3,1) * t243 + rSges(3,2) * t246;
t184 = Icges(3,3) * t244 + t247 * t273;
t183 = -Icges(3,3) * t247 + t244 * t273;
t166 = t294 * t247;
t165 = t294 * t244;
t154 = t350 + (pkin(1) - t336) * t247 + t302;
t153 = t334 + t236 + (-pkin(1) - t284) * t244;
t152 = t293 * t247;
t151 = t293 * t244;
t140 = t259 * t247;
t139 = t259 * t244;
t137 = -t244 * t249 + t219 + t261;
t136 = (rSges(4,3) - t249) * t247 + (-t226 - t283) * t244;
t135 = t223 * t245 * t147;
t130 = t247 * (-t247 * t336 + t302) + (t244 * t284 - t334) * t244;
t129 = t260 + t289;
t128 = (rSges(5,3) - t237) * t247 + (-t206 - t282) * t244;
t125 = t306 * t247;
t124 = t306 * t244;
t101 = t286 * t247;
t100 = t286 * t244;
t92 = t256 * t247;
t91 = t256 * t244;
t90 = t108 * t319;
t85 = t289 + t127 + t303;
t84 = -t247 * t237 + (-t341 - t206 + (-rSges(6,3) - pkin(10)) * t223) * t244 + t281;
t83 = t297 * t247;
t82 = t297 * t244;
t81 = -t127 * t224 - t150 * t319;
t80 = t126 * t224 + t150 * t320;
t79 = t262 * t247;
t78 = t262 * t244;
t77 = -t109 * t224 - t144 * t319;
t75 = t117 + t305;
t74 = t254 * t247;
t73 = t254 * t244;
t72 = -t224 * t145 - t223 * t321 + t135;
t71 = t289 + t351;
t70 = (pkin(5) * t242 - t237) * t247 + (-t224 * t225 - t206 + (-rSges(7,3) + t248) * t223) * t244 + t280;
t69 = (t126 * t247 - t127 * t244) * t223;
t65 = -t109 * t320 + t90;
t48 = t64 + t305;
t37 = -t224 * t349 + t307 * t319;
t36 = t115 * t224 + t138 * t320 + t76;
t28 = t90 + (t115 * t247 - t244 * t349) * t223;
t26 = t29 + t305;
t18 = t25 + t305;
t1 = [t230 * t203 + t228 * t204 + t246 * (Icges(3,2) * t246 + t328) + t243 * (Icges(3,1) * t243 + t327) + Icges(2,3) + t131 + t135 + (-t141 - t145 + t195) * t224 + (t196 - t321 - t322) * t223 + m(7) * (t70 ^ 2 + t71 ^ 2) + m(6) * (t84 ^ 2 + t85 ^ 2) + m(5) * (t128 ^ 2 + t129 ^ 2) + m(4) * (t136 ^ 2 + t137 ^ 2) + m(3) * (t153 ^ 2 + t154 ^ 2) + m(2) * (t217 ^ 2 + t218 ^ 2); m(3) * (-t153 * t247 - t154 * t244) * t216 + t250 + m(5) * (t128 * t140 + t129 * t139) + m(6) * (t84 * t92 + t85 * t91) + m(7) * (t70 * t74 + t71 * t73) + m(4) * (t136 * t166 + t137 * t165) + (t239 / 0.2e1 + t238 / 0.2e1) * (Icges(3,5) * t243 + Icges(3,6) * t246) + (t246 * (-Icges(3,6) * t247 + t244 * t276) + t243 * (-Icges(3,5) * t247 + t244 * t279)) * t344 + (t246 * (Icges(3,6) * t244 + t247 * t276) + t243 * (Icges(3,5) * t244 + t247 * t279)) * t345; m(7) * (t18 ^ 2 + t73 ^ 2 + t74 ^ 2) + m(6) * (t26 ^ 2 + t91 ^ 2 + t92 ^ 2) + m(5) * (t139 ^ 2 + t140 ^ 2 + t48 ^ 2) + m(4) * (t165 ^ 2 + t166 ^ 2 + t75 ^ 2) + m(3) * (t216 ^ 2 * t300 + t130 ^ 2) + t244 * t238 * t184 + t291 + (-t239 * t183 + (-t244 * t183 + t247 * t184) * t244 + t348) * t247; t250 + m(7) * (t70 * t79 + t71 * t78) + m(6) * (t100 * t85 + t101 * t84) + m(5) * (t128 * t152 + t129 * t151) + m(4) * (-t136 * t247 - t137 * t244) * t205; m(7) * (t18 * t25 + t73 * t78 + t74 * t79) + m(6) * (t100 * t91 + t101 * t92 + t26 * t29) + m(5) * (t139 * t151 + t140 * t152 + t48 * t64) + m(4) * (t117 * t75 + (-t165 * t244 - t166 * t247) * t205) + t252; m(7) * (t25 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(6) * (t100 ^ 2 + t101 ^ 2 + t29 ^ 2) + m(5) * (t151 ^ 2 + t152 ^ 2 + t64 ^ 2) + m(4) * (t205 ^ 2 * t300 + t117 ^ 2) + t252; m(7) * (t70 * t83 + t71 * t82) + m(6) * (t124 * t85 + t125 * t84) + m(5) * (-t128 * t247 - t129 * t244) * t197 + t251; m(7) * (t18 * t27 + t73 * t82 + t74 * t83) + m(6) * (t124 * t91 + t125 * t92 + t26 * t59) + m(5) * (t110 * t48 + (-t139 * t244 - t140 * t247) * t197) + t255; m(7) * (t25 * t27 + t78 * t82 + t79 * t83) + m(6) * (t100 * t124 + t101 * t125 + t29 * t59) + m(5) * (t110 * t64 + (-t151 * t244 - t152 * t247) * t197) + t255; m(7) * (t27 ^ 2 + t82 ^ 2 + t83 ^ 2) + m(6) * (t124 ^ 2 + t125 ^ 2 + t59 ^ 2) + m(5) * (t197 ^ 2 * t300 + t110 ^ 2) + t255; (-t68 - t72) * t224 + m(7) * (t36 * t70 + t37 * t71) + m(6) * (t80 * t84 + t81 * t85) + ((t54 / 0.2e1 + t67 / 0.2e1) * t247 + (t53 / 0.2e1 + t66 / 0.2e1) * t244) * t223 + t290; m(7) * (t18 * t28 + t36 * t74 + t37 * t73) + m(6) * (t26 * t69 + t80 * t92 + t81 * t91) + t253; m(7) * (t25 * t28 + t36 * t79 + t37 * t78) + m(6) * (t100 * t81 + t101 * t80 + t29 * t69) + t253; m(7) * (t27 * t28 + t36 * t83 + t37 * t82) + m(6) * (t124 * t81 + t125 * t80 + t59 * t69) + t253; (t72 * t224 - t9) * t224 + m(7) * (t28 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t69 ^ 2 + t80 ^ 2 + t81 ^ 2) + (t247 * t11 + t244 * t10 - t224 * (t244 * t53 + t247 * t54)) * t223 + t347; -t329 + m(7) * (t70 * t76 + t71 * t77) + t290; m(7) * (t18 * t65 + t73 * t77 + t74 * t76) + t287; m(7) * (t25 * t65 + t76 * t79 + t77 * t78) + t287; m(7) * (t27 * t65 + t76 * t83 + t77 * t82) + t287; m(7) * (t28 * t65 + t36 * t76 + t37 * t77) + t288; m(7) * (t65 ^ 2 + t76 ^ 2 + t77 ^ 2) + t288;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
