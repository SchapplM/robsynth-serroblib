% Calculate joint inertia matrix for
% S6RRRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:32
% EndTime: 2019-03-09 17:29:48
% DurationCPUTime: 6.30s
% Computational Cost: add. (24132->648), mult. (49813->892), div. (0->0), fcn. (63554->12), ass. (0->295)
t283 = cos(pkin(6));
t287 = sin(qJ(1));
t288 = cos(qJ(2));
t343 = t287 * t288;
t286 = sin(qJ(2));
t289 = cos(qJ(1));
t344 = t286 * t289;
t265 = t283 * t344 + t343;
t285 = sin(qJ(3));
t281 = sin(pkin(6));
t346 = t281 * t289;
t354 = cos(qJ(3));
t244 = t265 * t354 - t285 * t346;
t342 = t288 * t289;
t345 = t286 * t287;
t264 = -t283 * t342 + t345;
t328 = pkin(11) + qJ(5);
t278 = sin(t328);
t308 = cos(t328);
t199 = t244 * t278 - t264 * t308;
t200 = t244 * t308 + t264 * t278;
t356 = rSges(7,3) + qJ(6);
t357 = rSges(7,1) + pkin(5);
t358 = -t356 * t199 - t357 * t200;
t349 = t281 * t286;
t262 = -t283 * t354 + t285 * t349;
t310 = t281 * t354;
t263 = t283 * t285 + t286 * t310;
t347 = t281 * t288;
t212 = Icges(4,5) * t263 - Icges(4,6) * t262 - Icges(4,3) * t347;
t213 = Icges(4,4) * t263 - Icges(4,2) * t262 - Icges(4,6) * t347;
t214 = Icges(4,1) * t263 - Icges(4,4) * t262 - Icges(4,5) * t347;
t106 = -t212 * t347 - t262 * t213 + t263 * t214;
t231 = t263 * t278 + t308 * t347;
t232 = t263 * t308 - t278 * t347;
t159 = Icges(7,5) * t232 + Icges(7,6) * t262 + Icges(7,3) * t231;
t161 = Icges(7,4) * t232 + Icges(7,2) * t262 + Icges(7,6) * t231;
t163 = Icges(7,1) * t232 + Icges(7,4) * t262 + Icges(7,5) * t231;
t82 = t231 * t159 + t262 * t161 + t232 * t163;
t160 = Icges(6,5) * t232 - Icges(6,6) * t231 + Icges(6,3) * t262;
t162 = Icges(6,4) * t232 - Icges(6,2) * t231 + Icges(6,6) * t262;
t164 = Icges(6,1) * t232 - Icges(6,4) * t231 + Icges(6,5) * t262;
t83 = t262 * t160 - t231 * t162 + t232 * t164;
t280 = sin(pkin(11));
t282 = cos(pkin(11));
t241 = -t263 * t280 - t282 * t347;
t317 = t280 * t347;
t242 = t263 * t282 - t317;
t169 = Icges(5,5) * t242 + Icges(5,6) * t241 + Icges(5,3) * t262;
t170 = Icges(5,4) * t242 + Icges(5,2) * t241 + Icges(5,6) * t262;
t171 = Icges(5,1) * t242 + Icges(5,4) * t241 + Icges(5,5) * t262;
t86 = t262 * t169 + t241 * t170 + t242 * t171;
t355 = -t106 - t82 - t83 - t86;
t277 = pkin(4) * t282 + pkin(3);
t353 = -pkin(3) + t277;
t216 = Icges(3,5) * t265 - Icges(3,6) * t264 - Icges(3,3) * t346;
t352 = t216 * t289;
t351 = t264 * t280;
t266 = t283 * t343 + t344;
t350 = t266 * t280;
t348 = t281 * t287;
t243 = t265 * t285 + t289 * t310;
t234 = t243 * qJ(4);
t284 = -pkin(10) - qJ(4);
t323 = pkin(4) * t351;
t112 = -t243 * t284 + t244 * t353 - t234 + t323;
t190 = pkin(3) * t244 + t234;
t181 = t266 * t190;
t341 = t266 * t112 + t181;
t267 = -t283 * t345 + t342;
t245 = t267 * t285 - t287 * t310;
t246 = t267 * t354 + t285 * t348;
t191 = t246 * pkin(3) + qJ(4) * t245;
t312 = pkin(4) * t350 - t245 * t284 + t246 * t277;
t113 = -t191 + t312;
t340 = -t113 - t191;
t339 = rSges(7,2) * t243 - t358;
t201 = t246 * t278 - t266 * t308;
t202 = t246 * t308 + t266 * t278;
t338 = t245 * rSges(7,2) + t356 * t201 + t357 * t202;
t210 = -t246 * t280 + t266 * t282;
t211 = t246 * t282 + t350;
t141 = t211 * rSges(5,1) + t210 * rSges(5,2) + t245 * rSges(5,3);
t337 = -t141 - t191;
t336 = rSges(7,2) * t262 + t356 * t231 + t357 * t232;
t167 = -pkin(4) * t317 + t353 * t263 + (-qJ(4) - t284) * t262;
t227 = pkin(3) * t263 + qJ(4) * t262;
t335 = -t167 - t227;
t172 = rSges(5,1) * t242 + rSges(5,2) * t241 + rSges(5,3) * t262;
t334 = -t172 - t227;
t333 = t190 * t347 + t264 * t227;
t230 = t267 * pkin(2) + pkin(9) * t266;
t226 = t283 * t230;
t332 = t283 * t191 + t226;
t229 = pkin(2) * t265 + t264 * pkin(9);
t331 = -t190 - t229;
t330 = t229 * t348 + t230 * t346;
t329 = t289 * pkin(1) + pkin(8) * t348;
t114 = Icges(7,5) * t200 + Icges(7,6) * t243 + Icges(7,3) * t199;
t118 = Icges(7,4) * t200 + Icges(7,2) * t243 + Icges(7,6) * t199;
t122 = Icges(7,1) * t200 + Icges(7,4) * t243 + Icges(7,5) * t199;
t44 = t114 * t199 + t118 * t243 + t122 * t200;
t115 = Icges(7,5) * t202 + Icges(7,6) * t245 + Icges(7,3) * t201;
t119 = Icges(7,4) * t202 + Icges(7,2) * t245 + Icges(7,6) * t201;
t123 = Icges(7,1) * t202 + Icges(7,4) * t245 + Icges(7,5) * t201;
t45 = t115 * t199 + t119 * t243 + t123 * t200;
t66 = t159 * t199 + t161 * t243 + t163 * t200;
t1 = t243 * t44 + t245 * t45 + t262 * t66;
t116 = Icges(6,5) * t200 - Icges(6,6) * t199 + Icges(6,3) * t243;
t120 = Icges(6,4) * t200 - Icges(6,2) * t199 + Icges(6,6) * t243;
t124 = Icges(6,1) * t200 - Icges(6,4) * t199 + Icges(6,5) * t243;
t46 = t116 * t243 - t120 * t199 + t124 * t200;
t117 = Icges(6,5) * t202 - Icges(6,6) * t201 + Icges(6,3) * t245;
t121 = Icges(6,4) * t202 - Icges(6,2) * t201 + Icges(6,6) * t245;
t125 = Icges(6,1) * t202 - Icges(6,4) * t201 + Icges(6,5) * t245;
t47 = t117 * t243 - t121 * t199 + t125 * t200;
t67 = t160 * t243 - t162 * t199 + t164 * t200;
t2 = t243 * t46 + t245 * t47 + t262 * t67;
t327 = t2 / 0.2e1 + t1 / 0.2e1;
t48 = t114 * t201 + t118 * t245 + t122 * t202;
t49 = t115 * t201 + t119 * t245 + t123 * t202;
t68 = t159 * t201 + t161 * t245 + t163 * t202;
t3 = t243 * t48 + t245 * t49 + t262 * t68;
t50 = t116 * t245 - t120 * t201 + t124 * t202;
t51 = t117 * t245 - t121 * t201 + t125 * t202;
t69 = t160 * t245 - t162 * t201 + t164 * t202;
t4 = t243 * t50 + t245 * t51 + t262 * t69;
t326 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t264 * t44 + t266 * t45 - t347 * t66;
t6 = t264 * t46 + t266 * t47 - t347 * t67;
t325 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t264 * t48 + t266 * t49 - t347 * t68;
t8 = t264 * t50 + t266 * t51 - t347 * t69;
t324 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t67 * t283 + (t287 * t47 - t289 * t46) * t281;
t9 = t66 * t283 + (t287 * t45 - t289 * t44) * t281;
t322 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t68 * t283 + (t287 * t49 - t289 * t48) * t281;
t12 = t69 * t283 + (t287 * t51 - t289 * t50) * t281;
t321 = t11 / 0.2e1 + t12 / 0.2e1;
t56 = t114 * t231 + t118 * t262 + t122 * t232;
t57 = t115 * t231 + t119 * t262 + t123 * t232;
t76 = t82 * t262;
t17 = t56 * t243 + t57 * t245 + t76;
t58 = t116 * t262 - t120 * t231 + t124 * t232;
t59 = t117 * t262 - t121 * t231 + t125 * t232;
t77 = t83 * t262;
t18 = t58 * t243 + t59 * t245 + t77;
t320 = t18 / 0.2e1 + t17 / 0.2e1;
t19 = t56 * t264 + t57 * t266 - t347 * t82;
t20 = t58 * t264 + t59 * t266 - t347 * t83;
t319 = t19 / 0.2e1 + t20 / 0.2e1;
t80 = t82 * t283;
t21 = t80 + (t57 * t287 - t56 * t289) * t281;
t81 = t83 * t283;
t22 = t81 + (t59 * t287 - t58 * t289) * t281;
t318 = t22 / 0.2e1 + t21 / 0.2e1;
t316 = t283 * t113 + t332;
t315 = -t112 + t331;
t131 = t202 * rSges(6,1) - t201 * rSges(6,2) + t245 * rSges(6,3);
t314 = -t131 + t340;
t166 = rSges(6,1) * t232 - rSges(6,2) * t231 + rSges(6,3) * t262;
t313 = -t166 + t335;
t180 = t246 * rSges(4,1) - t245 * rSges(4,2) + t266 * rSges(4,3);
t252 = Icges(3,3) * t283 + (Icges(3,5) * t286 + Icges(3,6) * t288) * t281;
t253 = Icges(3,6) * t283 + (Icges(3,4) * t286 + Icges(3,2) * t288) * t281;
t254 = Icges(3,5) * t283 + (Icges(3,1) * t286 + Icges(3,4) * t288) * t281;
t311 = t283 * t252 + t253 * t347 + t254 * t349;
t223 = t267 * rSges(3,1) - t266 * rSges(3,2) + rSges(3,3) * t348;
t309 = -t287 * pkin(1) + pkin(8) * t346;
t215 = rSges(4,1) * t263 - rSges(4,2) * t262 - rSges(4,3) * t347;
t268 = (pkin(2) * t286 - pkin(9) * t288) * t281;
t307 = t281 * (-t215 - t268);
t306 = t112 * t347 + t264 * t167 + t333;
t305 = -t338 + t340;
t304 = t335 - t336;
t303 = t190 * t348 + t191 * t346 + t330;
t302 = t281 * (-t268 + t334);
t301 = -rSges(6,1) * t200 + rSges(6,2) * t199;
t300 = t230 + t329;
t299 = t281 * (-t268 + t313);
t298 = t59 / 0.2e1 + t68 / 0.2e1 + t69 / 0.2e1 + t57 / 0.2e1;
t297 = t67 / 0.2e1 + t66 / 0.2e1 + t58 / 0.2e1 + t56 / 0.2e1;
t296 = t112 * t348 + t113 * t346 + t303;
t295 = t281 * (-t268 + t304);
t294 = -t229 + t309;
t179 = rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t264;
t208 = -t244 * t280 + t264 * t282;
t209 = t244 * t282 + t351;
t140 = rSges(5,1) * t209 + rSges(5,2) * t208 + rSges(5,3) * t243;
t222 = t265 * rSges(3,1) - t264 * rSges(3,2) - rSges(3,3) * t346;
t293 = t300 + t312;
t101 = t212 * t266 - t213 * t245 + t214 * t246;
t135 = Icges(5,5) * t211 + Icges(5,6) * t210 + Icges(5,3) * t245;
t137 = Icges(5,4) * t211 + Icges(5,2) * t210 + Icges(5,6) * t245;
t139 = Icges(5,1) * t211 + Icges(5,4) * t210 + Icges(5,5) * t245;
t64 = t135 * t262 + t137 * t241 + t139 * t242;
t75 = t169 * t245 + t170 * t210 + t171 * t211;
t174 = Icges(4,5) * t246 - Icges(4,6) * t245 + Icges(4,3) * t266;
t176 = Icges(4,4) * t246 - Icges(4,2) * t245 + Icges(4,6) * t266;
t178 = Icges(4,1) * t246 - Icges(4,4) * t245 + Icges(4,5) * t266;
t96 = -t174 * t347 - t176 * t262 + t178 * t263;
t292 = t96 / 0.2e1 + t64 / 0.2e1 + t101 / 0.2e1 + t75 / 0.2e1 + t298;
t100 = t212 * t264 - t213 * t243 + t214 * t244;
t134 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t243;
t136 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t243;
t138 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t243;
t63 = t134 * t262 + t136 * t241 + t138 * t242;
t74 = t169 * t243 + t170 * t208 + t171 * t209;
t173 = Icges(4,5) * t244 - Icges(4,6) * t243 + Icges(4,3) * t264;
t175 = Icges(4,4) * t244 - Icges(4,2) * t243 + Icges(4,6) * t264;
t177 = Icges(4,1) * t244 - Icges(4,4) * t243 + Icges(4,5) * t264;
t95 = -t173 * t347 - t175 * t262 + t177 * t263;
t291 = t74 / 0.2e1 + t95 / 0.2e1 + t63 / 0.2e1 + t100 / 0.2e1 + t297;
t290 = -t244 * t277 + t294 - t323;
t270 = rSges(2,1) * t289 - t287 * rSges(2,2);
t269 = -t287 * rSges(2,1) - rSges(2,2) * t289;
t255 = rSges(3,3) * t283 + (rSges(3,1) * t286 + rSges(3,2) * t288) * t281;
t221 = Icges(3,1) * t267 - Icges(3,4) * t266 + Icges(3,5) * t348;
t220 = Icges(3,1) * t265 - Icges(3,4) * t264 - Icges(3,5) * t346;
t219 = Icges(3,4) * t267 - Icges(3,2) * t266 + Icges(3,6) * t348;
t218 = Icges(3,4) * t265 - Icges(3,2) * t264 - Icges(3,6) * t346;
t217 = Icges(3,5) * t267 - Icges(3,6) * t266 + Icges(3,3) * t348;
t204 = t223 + t329;
t203 = -t222 + t309;
t184 = -t283 * t222 - t255 * t346;
t183 = t223 * t283 - t255 * t348;
t168 = t311 * t283;
t157 = (t222 * t287 + t223 * t289) * t281;
t156 = t252 * t348 - t253 * t266 + t254 * t267;
t155 = -t252 * t346 - t264 * t253 + t265 * t254;
t143 = t300 + t180;
t142 = -t179 + t294;
t133 = -t180 * t347 - t215 * t266;
t132 = t179 * t347 + t215 * t264;
t129 = rSges(6,3) * t243 - t301;
t127 = t217 * t283 + (t219 * t288 + t221 * t286) * t281;
t126 = t216 * t283 + (t218 * t288 + t220 * t286) * t281;
t105 = t106 * t283;
t104 = t179 * t266 - t180 * t264;
t103 = (-t179 - t229) * t283 + t289 * t307;
t102 = t180 * t283 + t287 * t307 + t226;
t99 = t300 - t337;
t98 = -t140 - t190 + t294;
t97 = (t179 * t287 + t180 * t289) * t281 + t330;
t94 = t293 + t131;
t93 = (-rSges(6,3) + t284) * t243 + t290 + t301;
t92 = t131 * t262 - t166 * t245;
t91 = -t129 * t262 + t166 * t243;
t90 = t174 * t266 - t176 * t245 + t178 * t246;
t89 = t173 * t266 - t175 * t245 + t177 * t246;
t88 = t174 * t264 - t176 * t243 + t178 * t244;
t87 = t173 * t264 - t175 * t243 + t177 * t244;
t85 = t86 * t283;
t84 = t129 * t245 - t131 * t243;
t79 = t266 * t334 + t337 * t347;
t78 = t140 * t347 + t172 * t264 + t333;
t73 = t293 + t338;
t72 = (-rSges(7,2) + t284) * t243 + t290 + t358;
t71 = (-t140 + t331) * t283 + t289 * t302;
t70 = t141 * t283 + t287 * t302 + t332;
t65 = t140 * t266 + t264 * t337 + t181;
t62 = (t140 * t287 + t141 * t289) * t281 + t303;
t61 = -t245 * t336 + t262 * t338;
t60 = t243 * t336 - t262 * t339;
t55 = t135 * t245 + t137 * t210 + t139 * t211;
t54 = t134 * t245 + t136 * t210 + t138 * t211;
t53 = t135 * t243 + t137 * t208 + t139 * t209;
t52 = t134 * t243 + t136 * t208 + t138 * t209;
t43 = t266 * t313 + t314 * t347;
t42 = t129 * t347 + t166 * t264 + t306;
t41 = -t243 * t338 + t245 * t339;
t40 = (-t129 + t315) * t283 + t289 * t299;
t39 = t131 * t283 + t287 * t299 + t316;
t38 = t105 + (t96 * t287 - t95 * t289) * t281;
t37 = -t106 * t347 + t95 * t264 + t96 * t266;
t36 = t129 * t266 + t264 * t314 + t341;
t35 = (t129 * t287 + t131 * t289) * t281 + t296;
t34 = t266 * t304 + t305 * t347;
t33 = t264 * t336 + t339 * t347 + t306;
t32 = t101 * t283 + (t287 * t90 - t289 * t89) * t281;
t31 = t100 * t283 + (t287 * t88 - t289 * t87) * t281;
t30 = (t315 - t339) * t283 + t289 * t295;
t29 = t283 * t338 + t287 * t295 + t316;
t28 = -t101 * t347 + t264 * t89 + t266 * t90;
t27 = -t100 * t347 + t264 * t87 + t266 * t88;
t26 = t264 * t305 + t266 * t339 + t341;
t25 = (t287 * t339 + t289 * t338) * t281 + t296;
t24 = t85 + (t64 * t287 - t63 * t289) * t281;
t23 = t63 * t264 + t64 * t266 - t347 * t86;
t16 = t75 * t283 + (t287 * t55 - t289 * t54) * t281;
t15 = t74 * t283 + (t287 * t53 - t289 * t52) * t281;
t14 = t264 * t54 + t266 * t55 - t347 * t75;
t13 = t264 * t52 + t266 * t53 - t347 * t74;
t107 = [m(7) * (t72 ^ 2 + t73 ^ 2) + m(6) * (t93 ^ 2 + t94 ^ 2) + m(5) * (t98 ^ 2 + t99 ^ 2) + m(4) * (t142 ^ 2 + t143 ^ 2) + m(3) * (t203 ^ 2 + t204 ^ 2) + m(2) * (t269 ^ 2 + t270 ^ 2) + Icges(2,3) + t311 - t355; t105 + t80 + t81 + t168 + t85 + m(7) * (t29 * t73 + t30 * t72) + m(6) * (t39 * t94 + t40 * t93) + m(5) * (t70 * t99 + t71 * t98) + m(4) * (t102 * t143 + t103 * t142) + m(3) * (t183 * t204 + t184 * t203) + ((-t126 / 0.2e1 - t155 / 0.2e1 - t291) * t289 + (t127 / 0.2e1 + t156 / 0.2e1 + t292) * t287) * t281; (t22 + t21 + t24 + t38 + t168) * t283 + m(7) * (t25 ^ 2 + t29 ^ 2 + t30 ^ 2) + m(6) * (t35 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(5) * (t62 ^ 2 + t70 ^ 2 + t71 ^ 2) + m(4) * (t102 ^ 2 + t103 ^ 2 + t97 ^ 2) + m(3) * (t157 ^ 2 + t183 ^ 2 + t184 ^ 2) + ((-t9 - t10 - t31 - t15 + (-t264 * t218 + t265 * t220 - t281 * t352) * t346) * t289 + (t12 + t11 + t32 + t16 + ((-t219 * t266 + t221 * t267 + (t217 * t287 - t352) * t281) * t287 + (t217 * t346 + t218 * t266 + t264 * t219 - t220 * t267 - t265 * t221) * t289) * t281) * t287 + ((-t126 - t155) * t289 + (t127 + t156) * t287) * t283) * t281; t355 * t347 + m(7) * (t33 * t72 + t34 * t73) + m(6) * (t42 * t93 + t43 * t94) + m(5) * (t78 * t98 + t79 * t99) + m(4) * (t132 * t142 + t133 * t143) + t292 * t266 + t291 * t264; (t23 / 0.2e1 + t37 / 0.2e1 + t319) * t283 + (t16 / 0.2e1 + t32 / 0.2e1 + t321) * t266 + (t15 / 0.2e1 + t31 / 0.2e1 + t322) * t264 + m(7) * (t25 * t26 + t29 * t34 + t30 * t33) + m(6) * (t35 * t36 + t39 * t43 + t40 * t42) + m(5) * (t62 * t65 + t70 * t79 + t71 * t78) + m(4) * (t102 * t133 + t103 * t132 + t104 * t97) + ((-t13 / 0.2e1 - t27 / 0.2e1 - t325) * t289 + (-t24 / 0.2e1 - t38 / 0.2e1 - t318) * t288 + (t14 / 0.2e1 + t28 / 0.2e1 + t324) * t287) * t281; (-t19 - t20 - t23 - t37) * t347 + (t7 + t8 + t28 + t14) * t266 + (t6 + t5 + t27 + t13) * t264 + m(7) * (t26 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t36 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t65 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t104 ^ 2 + t132 ^ 2 + t133 ^ 2); m(7) * (t243 * t73 + t245 * t72) + m(6) * (t243 * t94 + t245 * t93) + m(5) * (t243 * t99 + t245 * t98); m(7) * (t243 * t29 + t245 * t30 + t25 * t262) + m(6) * (t243 * t39 + t245 * t40 + t262 * t35) + m(5) * (t243 * t70 + t245 * t71 + t262 * t62); m(7) * (t243 * t34 + t245 * t33 + t26 * t262) + m(6) * (t243 * t43 + t245 * t42 + t262 * t36) + m(5) * (t243 * t79 + t245 * t78 + t262 * t65); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t243 ^ 2 + t245 ^ 2 + t262 ^ 2); t77 + t76 + m(7) * (t60 * t72 + t61 * t73) + m(6) * (t91 * t93 + t92 * t94) + t298 * t245 + t297 * t243; t320 * t283 + t318 * t262 + t321 * t245 + t322 * t243 + m(7) * (t25 * t41 + t29 * t61 + t30 * t60) + m(6) * (t35 * t84 + t39 * t92 + t40 * t91) + (t287 * t326 - t289 * t327) * t281; -t320 * t347 + t326 * t266 + t327 * t264 + t319 * t262 + t324 * t245 + t325 * t243 + m(7) * (t26 * t41 + t33 * t60 + t34 * t61) + m(6) * (t36 * t84 + t42 * t91 + t43 * t92); m(6) * (t243 * t92 + t245 * t91 + t262 * t84) + m(7) * (t243 * t61 + t245 * t60 + t262 * t41); (t17 + t18) * t262 + (t4 + t3) * t245 + (t2 + t1) * t243 + m(7) * (t41 ^ 2 + t60 ^ 2 + t61 ^ 2) + m(6) * (t84 ^ 2 + t91 ^ 2 + t92 ^ 2); m(7) * (t199 * t73 + t201 * t72); m(7) * (t199 * t29 + t201 * t30 + t231 * t25); m(7) * (t199 * t34 + t201 * t33 + t231 * t26); m(7) * (t199 * t243 + t201 * t245 + t231 * t262); m(7) * (t199 * t61 + t201 * t60 + t231 * t41); m(7) * (t199 ^ 2 + t201 ^ 2 + t231 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t107(1) t107(2) t107(4) t107(7) t107(11) t107(16); t107(2) t107(3) t107(5) t107(8) t107(12) t107(17); t107(4) t107(5) t107(6) t107(9) t107(13) t107(18); t107(7) t107(8) t107(9) t107(10) t107(14) t107(19); t107(11) t107(12) t107(13) t107(14) t107(15) t107(20); t107(16) t107(17) t107(18) t107(19) t107(20) t107(21);];
Mq  = res;
