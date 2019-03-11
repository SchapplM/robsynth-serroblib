% Calculate joint inertia matrix for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRPPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:38
% EndTime: 2019-03-09 15:38:49
% DurationCPUTime: 4.81s
% Computational Cost: add. (24504->653), mult. (38593->903), div. (0->0), fcn. (48323->14), ass. (0->307)
t327 = m(6) / 0.2e1 + m(7) / 0.2e1;
t368 = 0.2e1 * t327;
t328 = qJ(3) + pkin(11);
t281 = sin(t328);
t288 = cos(pkin(6));
t292 = sin(qJ(2));
t286 = sin(pkin(6));
t313 = cos(t328);
t304 = t286 * t313;
t251 = t288 * t281 + t292 * t304;
t284 = pkin(12) + qJ(6);
t280 = sin(t284);
t282 = cos(t284);
t295 = cos(qJ(2));
t349 = t286 * t295;
t220 = -t251 * t280 - t282 * t349;
t221 = t251 * t282 - t280 * t349;
t352 = t286 * t292;
t250 = t281 * t352 - t288 * t313;
t132 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t250;
t133 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t250;
t134 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t250;
t60 = t250 * t132 + t220 * t133 + t221 * t134;
t285 = sin(pkin(12));
t287 = cos(pkin(12));
t231 = -t251 * t285 - t287 * t349;
t324 = t285 * t349;
t232 = t251 * t287 - t324;
t140 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t250;
t141 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t250;
t142 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t250;
t63 = t250 * t140 + t231 * t141 + t232 * t142;
t367 = -t60 - t63;
t296 = cos(qJ(1));
t348 = t286 * t296;
t343 = t295 * t296;
t293 = sin(qJ(1));
t346 = t292 * t293;
t266 = -t288 * t346 + t343;
t351 = t286 * t293;
t236 = t266 * t313 + t281 * t351;
t344 = t293 * t295;
t345 = t292 * t296;
t265 = t288 * t344 + t345;
t184 = -t236 * t280 + t265 * t282;
t185 = t236 * t282 + t265 * t280;
t235 = t266 * t281 - t293 * t304;
t107 = t185 * rSges(7,1) + t184 * rSges(7,2) + t235 * rSges(7,3);
t278 = pkin(5) * t287 + pkin(4);
t290 = -pkin(10) - qJ(5);
t353 = t265 * t285;
t366 = pkin(5) * t353 - t235 * t290 + t236 * t278 + t107;
t365 = t286 ^ 2;
t264 = t288 * t345 + t344;
t233 = t264 * t281 + t296 * t304;
t234 = t264 * t313 - t281 * t348;
t263 = -t288 * t343 + t346;
t182 = -t234 * t280 + t263 * t282;
t183 = t234 * t282 + t263 * t280;
t100 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t233;
t102 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t233;
t104 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t233;
t38 = t100 * t235 + t102 * t184 + t104 * t185;
t101 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t235;
t103 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t235;
t105 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t235;
t39 = t101 * t235 + t103 * t184 + t105 * t185;
t53 = t132 * t235 + t133 * t184 + t134 * t185;
t2 = t233 * t38 + t235 * t39 + t250 * t53;
t364 = t2 / 0.2e1;
t363 = t233 / 0.2e1;
t362 = t235 / 0.2e1;
t361 = t250 / 0.2e1;
t294 = cos(qJ(3));
t279 = pkin(3) * t294 + pkin(2);
t360 = -pkin(2) + t279;
t359 = -pkin(4) + t278;
t305 = -rSges(7,1) * t183 - rSges(7,2) * t182;
t106 = rSges(7,3) * t233 - t305;
t225 = t233 * qJ(5);
t355 = t263 * t285;
t358 = pkin(5) * t355 - t233 * t290 + t359 * t234 + t106 - t225;
t178 = t236 * pkin(4) + qJ(5) * t235;
t357 = -t178 + t366;
t207 = Icges(3,5) * t264 - Icges(3,6) * t263 - Icges(3,3) * t348;
t356 = t207 * t296;
t289 = -qJ(4) - pkin(9);
t354 = t263 * t289;
t350 = t286 * t294;
t291 = sin(qJ(3));
t347 = t288 * t291;
t135 = rSges(7,1) * t221 + rSges(7,2) * t220 + rSges(7,3) * t250;
t342 = -pkin(5) * t324 + t359 * t251 + (-qJ(5) - t290) * t250 + t135;
t259 = t263 * pkin(9);
t322 = t291 * t348;
t270 = pkin(3) * t322;
t158 = t264 * t360 - t259 - t270 - t354;
t136 = t265 * t158;
t177 = pkin(4) * t234 + t225;
t341 = t265 * t177 + t136;
t205 = pkin(3) * t347 + ((pkin(9) + t289) * t295 + t360 * t292) * t286;
t340 = t158 * t349 + t263 * t205;
t224 = t266 * pkin(2) + pkin(9) * t265;
t323 = t291 * t351;
t315 = pkin(3) * t323 - t265 * t289 + t266 * t279;
t159 = -t224 + t315;
t219 = t288 * t224;
t339 = t288 * t159 + t219;
t156 = t236 * rSges(5,1) - t235 * rSges(5,2) + t265 * rSges(5,3);
t338 = -t156 - t159;
t223 = pkin(2) * t264 + t259;
t337 = -t158 - t223;
t336 = -t159 - t178;
t239 = -t264 * t291 - t294 * t348;
t240 = t264 * t294 - t322;
t167 = rSges(4,1) * t240 + rSges(4,2) * t239 + rSges(4,3) * t263;
t335 = -t167 - t223;
t198 = Icges(5,4) * t251 - Icges(5,2) * t250 - Icges(5,6) * t349;
t199 = Icges(5,1) * t251 - Icges(5,4) * t250 - Icges(5,5) * t349;
t334 = -t250 * t198 + t251 * t199;
t261 = t288 * t294 - t291 * t352;
t262 = t292 * t350 + t347;
t203 = Icges(4,4) * t262 + Icges(4,2) * t261 - Icges(4,6) * t349;
t204 = Icges(4,1) * t262 + Icges(4,4) * t261 - Icges(4,5) * t349;
t333 = t261 * t203 + t262 * t204;
t200 = rSges(5,1) * t251 - rSges(5,2) * t250 - rSges(5,3) * t349;
t332 = -t200 - t205;
t201 = pkin(4) * t251 + qJ(5) * t250;
t331 = -t201 - t205;
t330 = t223 * t351 + t224 * t348;
t329 = t296 * pkin(1) + pkin(8) * t351;
t46 = t100 * t250 + t102 * t220 + t104 * t221;
t52 = t132 * t233 + t133 * t182 + t134 * t183;
t326 = t52 / 0.2e1 + t46 / 0.2e1;
t47 = t101 * t250 + t103 * t220 + t105 * t221;
t325 = t53 / 0.2e1 + t47 / 0.2e1;
t193 = -t236 * t285 + t265 * t287;
t194 = t236 * t287 + t353;
t115 = t194 * rSges(6,1) + t193 * rSges(6,2) + t235 * rSges(6,3);
t321 = -t115 + t336;
t143 = rSges(6,1) * t232 + rSges(6,2) * t231 + rSges(6,3) * t250;
t320 = -t143 + t331;
t319 = t288 * t178 + t339;
t318 = -t177 + t337;
t241 = -t266 * t291 + t293 * t350;
t242 = t266 * t294 + t323;
t168 = t242 * rSges(4,1) + t241 * rSges(4,2) + t265 * rSges(4,3);
t247 = Icges(3,3) * t288 + (Icges(3,5) * t292 + Icges(3,6) * t295) * t286;
t248 = Icges(3,6) * t288 + (Icges(3,4) * t292 + Icges(3,2) * t295) * t286;
t249 = Icges(3,5) * t288 + (Icges(3,1) * t292 + Icges(3,4) * t295) * t286;
t316 = t288 * t247 + t248 * t349 + t249 * t352;
t214 = t266 * rSges(3,1) - t265 * rSges(3,2) + rSges(3,3) * t351;
t314 = -t293 * pkin(1) + pkin(8) * t348;
t206 = rSges(4,1) * t262 + rSges(4,2) * t261 - rSges(4,3) * t349;
t267 = (pkin(2) * t292 - pkin(9) * t295) * t286;
t312 = t286 * (-t206 - t267);
t311 = t336 - t357;
t310 = t331 - t342;
t309 = t158 * t351 + t159 * t348 + t330;
t308 = t177 * t349 + t263 * t201 + t340;
t307 = t286 * (-t267 + t332);
t306 = -rSges(5,1) * t234 + rSges(5,2) * t233;
t303 = t315 + t329;
t302 = t286 * (-t267 + t320);
t301 = t177 * t351 + t178 * t348 + t309;
t300 = t286 * (-t267 + t310);
t299 = -t264 * t279 + t270 + t314;
t191 = -t234 * t285 + t263 * t287;
t192 = t234 * t287 + t355;
t114 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t233;
t213 = t264 * rSges(3,1) - t263 * rSges(3,2) - rSges(3,3) * t348;
t108 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t233;
t110 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t233;
t112 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t233;
t50 = t108 * t250 + t110 * t231 + t112 * t232;
t55 = t140 * t233 + t141 * t191 + t142 * t192;
t149 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t263;
t151 = Icges(5,4) * t234 - Icges(5,2) * t233 + Icges(5,6) * t263;
t153 = Icges(5,1) * t234 - Icges(5,4) * t233 + Icges(5,5) * t263;
t76 = -t149 * t349 - t151 * t250 + t153 * t251;
t161 = Icges(4,5) * t240 + Icges(4,6) * t239 + Icges(4,3) * t263;
t163 = Icges(4,4) * t240 + Icges(4,2) * t239 + Icges(4,6) * t263;
t165 = Icges(4,1) * t240 + Icges(4,4) * t239 + Icges(4,5) * t263;
t84 = -t161 * t349 + t163 * t261 + t165 * t262;
t197 = Icges(5,5) * t251 - Icges(5,6) * t250 - Icges(5,3) * t349;
t87 = t197 * t263 - t198 * t233 + t199 * t234;
t202 = Icges(4,5) * t262 + Icges(4,6) * t261 - Icges(4,3) * t349;
t91 = t202 * t263 + t203 * t239 + t204 * t240;
t298 = t84 / 0.2e1 + t76 / 0.2e1 + t50 / 0.2e1 + t91 / 0.2e1 + t87 / 0.2e1 + t55 / 0.2e1 + t326;
t109 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t235;
t111 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t235;
t113 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t235;
t51 = t109 * t250 + t111 * t231 + t113 * t232;
t56 = t140 * t235 + t141 * t193 + t142 * t194;
t150 = Icges(5,5) * t236 - Icges(5,6) * t235 + Icges(5,3) * t265;
t152 = Icges(5,4) * t236 - Icges(5,2) * t235 + Icges(5,6) * t265;
t154 = Icges(5,1) * t236 - Icges(5,4) * t235 + Icges(5,5) * t265;
t77 = -t150 * t349 - t152 * t250 + t154 * t251;
t162 = Icges(4,5) * t242 + Icges(4,6) * t241 + Icges(4,3) * t265;
t164 = Icges(4,4) * t242 + Icges(4,2) * t241 + Icges(4,6) * t265;
t166 = Icges(4,1) * t242 + Icges(4,4) * t241 + Icges(4,5) * t265;
t85 = -t162 * t349 + t164 * t261 + t166 * t262;
t88 = t197 * t265 - t198 * t235 + t199 * t236;
t92 = t202 * t265 + t203 * t241 + t204 * t242;
t297 = t56 / 0.2e1 + t85 / 0.2e1 + t77 / 0.2e1 + t51 / 0.2e1 + t92 / 0.2e1 + t88 / 0.2e1 + t325;
t272 = rSges(2,1) * t296 - t293 * rSges(2,2);
t271 = -t293 * rSges(2,1) - rSges(2,2) * t296;
t252 = rSges(3,3) * t288 + (rSges(3,1) * t292 + rSges(3,2) * t295) * t286;
t212 = Icges(3,1) * t266 - Icges(3,4) * t265 + Icges(3,5) * t351;
t211 = Icges(3,1) * t264 - Icges(3,4) * t263 - Icges(3,5) * t348;
t210 = Icges(3,4) * t266 - Icges(3,2) * t265 + Icges(3,6) * t351;
t209 = Icges(3,4) * t264 - Icges(3,2) * t263 - Icges(3,6) * t348;
t208 = Icges(3,5) * t266 - Icges(3,6) * t265 + Icges(3,3) * t351;
t196 = t214 + t329;
t195 = -t213 + t314;
t176 = -t288 * t213 - t252 * t348;
t175 = t214 * t288 - t252 * t351;
t160 = t316 * t288;
t155 = rSges(5,3) * t263 - t306;
t139 = (t213 * t293 + t214 * t296) * t286;
t138 = t247 * t351 - t248 * t265 + t249 * t266;
t137 = -t247 * t348 - t263 * t248 + t264 * t249;
t128 = t224 + t168 + t329;
t127 = t314 + t335;
t121 = -t168 * t349 - t206 * t265;
t120 = t167 * t349 + t206 * t263;
t119 = t208 * t288 + (t210 * t295 + t212 * t292) * t286;
t118 = t207 * t288 + (t209 * t295 + t211 * t292) * t286;
t117 = t303 + t156;
t116 = (-rSges(5,3) + t289) * t263 + t299 + t306;
t97 = -t202 * t349 + t333;
t96 = t97 * t288;
t95 = t167 * t265 - t168 * t263;
t94 = t288 * t335 + t296 * t312;
t93 = t168 * t288 + t293 * t312 + t219;
t90 = -t197 * t349 + t334;
t89 = t90 * t288;
t86 = (t167 * t293 + t168 * t296) * t286 + t330;
t83 = t178 + t303 + t115;
t82 = -t114 - t177 + t299 + t354;
t81 = t162 * t265 + t164 * t241 + t166 * t242;
t80 = t161 * t265 + t163 * t241 + t165 * t242;
t79 = t162 * t263 + t164 * t239 + t166 * t240;
t78 = t161 * t263 + t163 * t239 + t165 * t240;
t75 = t265 * t332 + t338 * t349;
t74 = t155 * t349 + t200 * t263 + t340;
t73 = t107 * t250 - t135 * t235;
t72 = -t106 * t250 + t135 * t233;
t71 = t303 + t366;
t70 = -t234 * t278 + (-pkin(5) * t285 + t289) * t263 + (-rSges(7,3) + t290) * t233 + t299 + t305;
t69 = t150 * t265 - t152 * t235 + t154 * t236;
t68 = t149 * t265 - t151 * t235 + t153 * t236;
t67 = t150 * t263 - t152 * t233 + t154 * t234;
t66 = t149 * t263 - t151 * t233 + t153 * t234;
t65 = (-t155 + t337) * t288 + t296 * t307;
t64 = t156 * t288 + t293 * t307 + t339;
t62 = t63 * t288;
t61 = t106 * t235 - t107 * t233;
t59 = t60 * t288;
t58 = t155 * t265 + t263 * t338 + t136;
t57 = t60 * t250;
t54 = (t155 * t293 + t156 * t296) * t286 + t309;
t49 = t265 * t320 + t321 * t349;
t48 = t114 * t349 + t143 * t263 + t308;
t45 = (-t114 + t318) * t288 + t296 * t302;
t44 = t115 * t288 + t293 * t302 + t319;
t43 = t109 * t235 + t111 * t193 + t113 * t194;
t42 = t108 * t235 + t110 * t193 + t112 * t194;
t41 = t109 * t233 + t111 * t191 + t113 * t192;
t40 = t108 * t233 + t110 * t191 + t112 * t192;
t37 = t101 * t233 + t103 * t182 + t105 * t183;
t36 = t100 * t233 + t102 * t182 + t104 * t183;
t35 = t114 * t265 + t263 * t321 + t341;
t34 = (t114 * t293 + t115 * t296) * t286 + t301;
t33 = t96 + (t85 * t293 - t84 * t296) * t286;
t32 = t84 * t263 + t85 * t265 - t349 * t97;
t31 = t92 * t288 + (t293 * t81 - t296 * t80) * t286;
t30 = t91 * t288 + (t293 * t79 - t296 * t78) * t286;
t29 = t89 + (t77 * t293 - t76 * t296) * t286;
t28 = t265 * t310 + t311 * t349;
t27 = t263 * t342 + t349 * t358 + t308;
t26 = t263 * t80 + t265 * t81 - t349 * t92;
t25 = t263 * t78 + t265 * t79 - t349 * t91;
t24 = t76 * t263 + t77 * t265 - t349 * t90;
t23 = (t318 - t358) * t288 + t296 * t300;
t22 = t288 * t357 + t293 * t300 + t319;
t21 = t88 * t288 + (t293 * t69 - t296 * t68) * t286;
t20 = t87 * t288 + (t293 * t67 - t296 * t66) * t286;
t19 = t263 * t68 + t265 * t69 - t349 * t88;
t18 = t263 * t66 + t265 * t67 - t349 * t87;
t17 = t263 * t311 + t265 * t358 + t341;
t16 = (t293 * t358 + t296 * t357) * t286 + t301;
t15 = t62 + (t51 * t293 - t50 * t296) * t286;
t14 = t50 * t263 + t51 * t265 - t349 * t63;
t13 = t59 + (t47 * t293 - t46 * t296) * t286;
t12 = t46 * t263 + t47 * t265 - t349 * t60;
t11 = t46 * t233 + t47 * t235 + t57;
t10 = t56 * t288 + (t293 * t43 - t296 * t42) * t286;
t9 = t55 * t288 + (t293 * t41 - t296 * t40) * t286;
t8 = t263 * t42 + t265 * t43 - t349 * t56;
t7 = t263 * t40 + t265 * t41 - t349 * t55;
t6 = t53 * t288 + (t293 * t39 - t296 * t38) * t286;
t5 = t52 * t288 + (t293 * t37 - t296 * t36) * t286;
t4 = t263 * t38 + t265 * t39 - t349 * t53;
t3 = t263 * t36 + t265 * t37 - t349 * t52;
t1 = t233 * t36 + t235 * t37 + t250 * t52;
t98 = [(-t197 - t202) * t349 + t316 + m(7) * (t70 ^ 2 + t71 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t116 ^ 2 + t117 ^ 2) + m(4) * (t127 ^ 2 + t128 ^ 2) + m(3) * (t195 ^ 2 + t196 ^ 2) + m(2) * (t271 ^ 2 + t272 ^ 2) + Icges(2,3) + t333 + t334 - t367; t160 + t62 + t96 + t89 + t59 + m(3) * (t175 * t196 + t176 * t195) + m(7) * (t22 * t71 + t23 * t70) + m(6) * (t44 * t83 + t45 * t82) + m(5) * (t116 * t65 + t117 * t64) + m(4) * (t127 * t94 + t128 * t93) + ((-t118 / 0.2e1 - t137 / 0.2e1 - t298) * t296 + (t119 / 0.2e1 + t138 / 0.2e1 + t297) * t293) * t286; (t13 + t15 + t33 + t29 + t160) * t288 + m(7) * (t16 ^ 2 + t22 ^ 2 + t23 ^ 2) + m(5) * (t54 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(6) * (t34 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(4) * (t86 ^ 2 + t93 ^ 2 + t94 ^ 2) + m(3) * (t139 ^ 2 + t175 ^ 2 + t176 ^ 2) + ((-t5 - t9 - t30 - t20 + ((-t263 * t209 + t264 * t211) * t286 - t365 * t356) * t296) * t296 + (t6 + t10 + t31 + t21 + ((-t210 * t265 + t212 * t266 + (t208 * t293 - t356) * t286) * t293 + (t208 * t348 + t209 * t265 + t263 * t210 - t211 * t266 - t264 * t212) * t296) * t286) * t293 + ((-t118 - t137) * t296 + (t119 + t138) * t293) * t288) * t286; (-t90 - t97 + t367) * t349 + m(7) * (t27 * t70 + t28 * t71) + m(6) * (t48 * t82 + t49 * t83) + m(5) * (t116 * t74 + t117 * t75) + m(4) * (t120 * t127 + t121 * t128) + t297 * t265 + t298 * t263; (t14 / 0.2e1 + t12 / 0.2e1 + t24 / 0.2e1 + t32 / 0.2e1) * t288 + (t10 / 0.2e1 + t6 / 0.2e1 + t21 / 0.2e1 + t31 / 0.2e1) * t265 + (t9 / 0.2e1 + t5 / 0.2e1 + t20 / 0.2e1 + t30 / 0.2e1) * t263 + m(7) * (t16 * t17 + t22 * t28 + t23 * t27) + m(6) * (t34 * t35 + t44 * t49 + t45 * t48) + m(5) * (t54 * t58 + t64 * t75 + t65 * t74) + m(4) * (t120 * t94 + t121 * t93 + t86 * t95) + ((-t18 / 0.2e1 - t25 / 0.2e1 - t3 / 0.2e1 - t7 / 0.2e1) * t296 + (-t13 / 0.2e1 - t15 / 0.2e1 - t29 / 0.2e1 - t33 / 0.2e1) * t295 + (t4 / 0.2e1 + t8 / 0.2e1 + t19 / 0.2e1 + t26 / 0.2e1) * t293) * t286; (-t12 - t14 - t24 - t32) * t349 + (t4 + t8 + t19 + t26) * t265 + (t3 + t7 + t18 + t25) * t263 + m(7) * (t17 ^ 2 + t27 ^ 2 + t28 ^ 2) + m(6) * (t35 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t58 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(4) * (t120 ^ 2 + t121 ^ 2 + t95 ^ 2); m(7) * (t263 * t71 + t265 * t70) + m(6) * (t263 * t83 + t265 * t82) + m(5) * (t116 * t265 + t117 * t263); m(7) * (-t16 * t349 + t22 * t263 + t23 * t265) + m(5) * (t263 * t64 + t265 * t65 - t349 * t54) + m(6) * (t263 * t44 + t265 * t45 - t34 * t349); m(7) * (-t17 * t349 + t263 * t28 + t265 * t27) + m(6) * (t263 * t49 + t265 * t48 - t349 * t35) + m(5) * (t263 * t75 + t265 * t74 - t349 * t58); 0.2e1 * (m(5) / 0.2e1 + t327) * (t295 ^ 2 * t365 + t263 ^ 2 + t265 ^ 2); m(7) * (t233 * t71 + t235 * t70) + m(6) * (t233 * t83 + t235 * t82); m(7) * (t16 * t250 + t22 * t233 + t23 * t235) + m(6) * (t233 * t44 + t235 * t45 + t250 * t34); m(7) * (t17 * t250 + t233 * t28 + t235 * t27) + m(6) * (t233 * t49 + t235 * t48 + t250 * t35); (t233 * t263 + t235 * t265 - t250 * t349) * t368; (t233 ^ 2 + t235 ^ 2 + t250 ^ 2) * t368; m(7) * (t70 * t72 + t71 * t73) + t57 + t325 * t235 + t326 * t233; t5 * t363 + m(7) * (t16 * t61 + t22 * t73 + t23 * t72) + t13 * t361 + t288 * t11 / 0.2e1 + t6 * t362 + (-t296 * t1 / 0.2e1 + t293 * t364) * t286; -t11 * t349 / 0.2e1 + t265 * t364 + t263 * t1 / 0.2e1 + m(7) * (t17 * t61 + t27 * t72 + t28 * t73) + t3 * t363 + t4 * t362 + t12 * t361; m(7) * (t263 * t73 + t265 * t72 - t349 * t61); m(7) * (t233 * t73 + t235 * t72 + t250 * t61); t235 * t2 + t233 * t1 + t250 * t11 + m(7) * (t61 ^ 2 + t72 ^ 2 + t73 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t98(1) t98(2) t98(4) t98(7) t98(11) t98(16); t98(2) t98(3) t98(5) t98(8) t98(12) t98(17); t98(4) t98(5) t98(6) t98(9) t98(13) t98(18); t98(7) t98(8) t98(9) t98(10) t98(14) t98(19); t98(11) t98(12) t98(13) t98(14) t98(15) t98(20); t98(16) t98(17) t98(18) t98(19) t98(20) t98(21);];
Mq  = res;
