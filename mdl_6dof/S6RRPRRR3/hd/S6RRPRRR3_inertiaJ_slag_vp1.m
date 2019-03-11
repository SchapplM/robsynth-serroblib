% Calculate joint inertia matrix for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:51
% EndTime: 2019-03-09 13:22:02
% DurationCPUTime: 4.70s
% Computational Cost: add. (16585->524), mult. (13700->734), div. (0->0), fcn. (14786->12), ass. (0->260)
t349 = Icges(3,3) + Icges(4,3);
t240 = qJ(2) + pkin(11);
t231 = sin(t240);
t232 = cos(t240);
t247 = sin(qJ(2));
t250 = cos(qJ(2));
t348 = Icges(3,5) * t250 + Icges(4,5) * t232 - Icges(3,6) * t247 - Icges(4,6) * t231;
t347 = t231 / 0.2e1;
t346 = t232 / 0.2e1;
t345 = t247 / 0.2e1;
t344 = t250 / 0.2e1;
t244 = qJ(4) + qJ(5);
t235 = qJ(6) + t244;
t227 = sin(t235);
t228 = cos(t235);
t248 = sin(qJ(1));
t251 = cos(qJ(1));
t313 = t232 * t251;
t176 = -t227 * t313 + t228 * t248;
t177 = t227 * t248 + t228 * t313;
t316 = t231 * t251;
t111 = t177 * rSges(7,1) + t176 * rSges(7,2) + rSges(7,3) * t316;
t249 = cos(qJ(4));
t229 = t249 * pkin(4) + pkin(3);
t234 = cos(t244);
t206 = pkin(5) * t234 + t229;
t233 = sin(t244);
t246 = sin(qJ(4));
t207 = pkin(4) * t246 + pkin(5) * t233;
t343 = t206 * t313 + t248 * t207 + t111;
t342 = -t348 * t248 + t349 * t251;
t341 = t349 * t248 + t348 * t251;
t242 = t248 ^ 2;
t243 = t251 ^ 2;
t252 = -pkin(9) - pkin(8);
t317 = t231 * t248;
t314 = t232 * t248;
t174 = -t227 * t314 - t228 * t251;
t175 = -t227 * t251 + t228 * t314;
t104 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t317;
t106 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t317;
t108 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t317;
t38 = t104 * t317 + t106 * t174 + t108 * t175;
t105 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t316;
t107 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t316;
t109 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t316;
t39 = t105 * t317 + t107 * t174 + t109 * t175;
t148 = -Icges(7,3) * t232 + (Icges(7,5) * t228 - Icges(7,6) * t227) * t231;
t149 = -Icges(7,6) * t232 + (Icges(7,4) * t228 - Icges(7,2) * t227) * t231;
t150 = -Icges(7,5) * t232 + (Icges(7,1) * t228 - Icges(7,4) * t227) * t231;
t65 = t148 * t317 + t149 * t174 + t150 * t175;
t5 = -t232 * t65 + (t248 * t38 + t251 * t39) * t231;
t40 = t104 * t316 + t106 * t176 + t108 * t177;
t41 = t105 * t316 + t107 * t176 + t109 * t177;
t66 = t148 * t316 + t149 * t176 + t150 * t177;
t6 = -t232 * t66 + (t248 * t40 + t251 * t41) * t231;
t340 = t6 * t316 + t5 * t317;
t339 = -t232 / 0.2e1;
t338 = t248 / 0.2e1;
t337 = -t251 / 0.2e1;
t336 = pkin(2) * t247;
t335 = pkin(3) * t232;
t334 = pkin(8) * t231;
t333 = -pkin(3) + t229;
t139 = t231 * t228 * t150;
t320 = t149 * t227;
t75 = -t232 * t148 - t231 * t320 + t139;
t327 = t75 * t232;
t48 = -t232 * t104 + (-t106 * t227 + t108 * t228) * t231;
t49 = -t232 * t105 + (-t107 * t227 + t109 * t228) * t231;
t13 = -t327 + (t248 * t48 + t251 * t49) * t231;
t309 = t234 * t251;
t312 = t233 * t248;
t182 = -t232 * t312 - t309;
t310 = t234 * t248;
t311 = t233 * t251;
t183 = t232 * t310 - t311;
t113 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t317;
t115 = Icges(6,4) * t183 + Icges(6,2) * t182 + Icges(6,6) * t317;
t117 = Icges(6,1) * t183 + Icges(6,4) * t182 + Icges(6,5) * t317;
t55 = -t232 * t113 + (-t115 * t233 + t117 * t234) * t231;
t184 = -t232 * t311 + t310;
t185 = t232 * t309 + t312;
t114 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t316;
t116 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t316;
t118 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t316;
t56 = -t232 * t114 + (-t116 * t233 + t118 * t234) * t231;
t154 = -Icges(6,5) * t232 + (Icges(6,1) * t234 - Icges(6,4) * t233) * t231;
t141 = t231 * t234 * t154;
t152 = -Icges(6,3) * t232 + (Icges(6,5) * t234 - Icges(6,6) * t233) * t231;
t153 = -Icges(6,6) * t232 + (Icges(6,4) * t234 - Icges(6,2) * t233) * t231;
t319 = t153 * t233;
t80 = -t232 * t152 - t231 * t319 + t141;
t332 = -t13 + t80 * t232 - (t248 * t55 + t251 * t56) * t231;
t331 = -t75 - t80;
t330 = rSges(3,1) * t250;
t329 = rSges(3,2) * t247;
t328 = t251 * rSges(3,3);
t271 = -rSges(7,1) * t175 - rSges(7,2) * t174;
t110 = rSges(7,3) * t317 - t271;
t100 = t110 * t316;
t241 = -pkin(10) + t252;
t306 = t246 * t251;
t315 = t231 * t252;
t296 = pkin(4) * t306 + t248 * t315;
t297 = t206 - t229;
t98 = -t207 * t251 + (-t231 * t241 + t232 * t297) * t248 + t296;
t326 = t98 * t316 + t100;
t293 = t241 - t252;
t307 = t246 * t248;
t298 = -pkin(4) * t307 - t229 * t313;
t325 = -t293 * t316 + t298 + t343;
t324 = Icges(3,4) * t247;
t323 = Icges(3,4) * t250;
t322 = Icges(4,4) * t231;
t321 = Icges(4,4) * t232;
t161 = -Icges(5,6) * t232 + (Icges(5,4) * t249 - Icges(5,2) * t246) * t231;
t318 = t161 * t246;
t245 = -qJ(3) - pkin(7);
t308 = t245 * t251;
t305 = t248 * t249;
t304 = t249 * t251;
t120 = t185 * rSges(6,1) + t184 * rSges(6,2) + rSges(6,3) * t316;
t258 = -t251 * t315 - t298;
t295 = pkin(3) * t313 + pkin(8) * t316;
t128 = t258 - t295;
t303 = -t120 - t128;
t127 = (t232 * t333 - t334) * t248 - t296;
t147 = (pkin(8) + t252) * t232 + t333 * t231;
t302 = t232 * t127 + t147 * t317;
t137 = t231 * t297 + t232 * t293;
t151 = -t232 * rSges(7,3) + (rSges(7,1) * t228 - rSges(7,2) * t227) * t231;
t301 = -t137 - t151;
t83 = t232 * t110 + t151 * t317;
t272 = -rSges(6,1) * t183 - rSges(6,2) * t182;
t119 = rSges(6,3) * t317 - t272;
t155 = -t232 * rSges(6,3) + (rSges(6,1) * t234 - rSges(6,2) * t233) * t231;
t88 = t232 * t119 + t155 * t317;
t230 = pkin(2) * t250 + pkin(1);
t222 = t251 * t230;
t239 = t251 * pkin(7);
t300 = t248 * (t308 + t239 + (-pkin(1) + t230) * t248) + t251 * (-pkin(1) * t251 + t222 + (-pkin(7) - t245) * t248);
t294 = t248 * rSges(3,3) + t251 * t330;
t292 = t242 + t243;
t195 = -t232 * t307 - t304;
t196 = t232 * t305 - t306;
t129 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t317;
t131 = Icges(5,4) * t196 + Icges(5,2) * t195 + Icges(5,6) * t317;
t133 = Icges(5,1) * t196 + Icges(5,4) * t195 + Icges(5,5) * t317;
t61 = -t129 * t232 + (-t131 * t246 + t133 * t249) * t231;
t160 = -Icges(5,3) * t232 + (Icges(5,5) * t249 - Icges(5,6) * t246) * t231;
t162 = -Icges(5,5) * t232 + (Icges(5,1) * t249 - Icges(5,4) * t246) * t231;
t76 = t160 * t317 + t161 * t195 + t162 * t196;
t291 = t61 / 0.2e1 + t76 / 0.2e1;
t197 = -t232 * t306 + t305;
t198 = t232 * t304 + t307;
t130 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t316;
t132 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t316;
t134 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t316;
t62 = -t130 * t232 + (-t132 * t246 + t134 * t249) * t231;
t77 = t160 * t316 + t161 * t197 + t162 * t198;
t290 = t77 / 0.2e1 + t62 / 0.2e1;
t289 = -t128 - t325;
t44 = t113 * t317 + t115 * t182 + t117 * t183;
t45 = t114 * t317 + t116 * t182 + t118 * t183;
t71 = t152 * t317 + t153 * t182 + t154 * t183;
t11 = -t232 * t71 + (t248 * t44 + t251 * t45) * t231;
t46 = t113 * t316 + t115 * t184 + t117 * t185;
t47 = t114 * t316 + t116 * t184 + t118 * t185;
t72 = t152 * t316 + t153 * t184 + t154 * t185;
t12 = -t232 * t72 + (t248 * t46 + t251 * t47) * t231;
t288 = t11 * t317 + t12 * t316 + t340;
t136 = t198 * rSges(5,1) + t197 * rSges(5,2) + rSges(5,3) * t316;
t287 = t317 / 0.2e1;
t286 = t316 / 0.2e1;
t285 = Icges(3,5) * t345 + Icges(4,5) * t347 + Icges(3,6) * t344 + Icges(4,6) * t346;
t284 = -rSges(4,1) * t231 - rSges(4,2) * t232 - t336;
t283 = -t231 * pkin(3) + t232 * pkin(8) - t336;
t282 = (t48 + t65) * t287 + (t49 + t66) * t286;
t281 = -t248 * t245 + t222;
t36 = t137 * t317 + t232 * t98 + t83;
t280 = -t232 * t13 + t340;
t279 = t242 * (t334 + t335) + t251 * t295 + t300;
t20 = t248 * t39 - t251 * t38;
t21 = t248 * t41 - t251 * t40;
t278 = t20 * t287 + t21 * t286 + t5 * t337 + t6 * t338 + (t49 * t248 - t48 * t251) * t339;
t277 = -t147 + t283;
t163 = -t232 * rSges(5,3) + (rSges(5,1) * t249 - rSges(5,2) * t246) * t231;
t276 = -t163 + t283;
t275 = -t329 + t330;
t274 = rSges(4,1) * t232 - rSges(4,2) * t231;
t273 = -rSges(5,1) * t196 - rSges(5,2) * t195;
t270 = Icges(3,1) * t250 - t324;
t269 = Icges(4,1) * t232 - t322;
t268 = -Icges(3,2) * t247 + t323;
t267 = -Icges(4,2) * t231 + t321;
t260 = -t155 + t277;
t259 = rSges(4,1) * t313 - rSges(4,2) * t316 + t248 * rSges(4,3);
t257 = t248 * t127 + t251 * t128 + t279;
t256 = t277 + t301;
t255 = t232 * t332 + t288;
t254 = t282 + (t55 + t71) * t287 + (t56 + t72) * t286;
t24 = t248 * t45 - t251 * t44;
t25 = t248 * t47 - t251 * t46;
t253 = t11 * t337 + t12 * t338 + t24 * t287 + t25 * t286 + t278 + (t56 * t248 - t55 * t251) * t339;
t215 = rSges(2,1) * t251 - rSges(2,2) * t248;
t214 = -rSges(2,1) * t248 - rSges(2,2) * t251;
t213 = rSges(3,1) * t247 + rSges(3,2) * t250;
t167 = t284 * t251;
t166 = t284 * t248;
t159 = pkin(7) * t248 + (pkin(1) - t329) * t251 + t294;
t158 = t328 + t239 + (-pkin(1) - t275) * t248;
t146 = t231 * t249 * t162;
t144 = t259 + t281;
t143 = (rSges(4,3) - t245) * t251 + (-t230 - t274) * t248;
t138 = t251 * (-t251 * t329 + t294) + (t248 * t275 - t328) * t248;
t135 = rSges(5,3) * t317 - t273;
t124 = t276 * t251;
t123 = t276 * t248;
t112 = t127 * t316;
t102 = t119 * t316;
t97 = t281 + t136 + t295;
t96 = -t308 + (-t335 - t230 + (-rSges(5,3) - pkin(8)) * t231) * t248 + t273;
t93 = -t136 * t232 - t163 * t316;
t92 = t135 * t232 + t163 * t317;
t91 = t260 * t251;
t90 = t260 * t248;
t89 = -t120 * t232 - t155 * t316;
t87 = -t232 * t160 - t231 * t318 + t146;
t86 = t258 + t281 + t120;
t85 = -t308 + (-rSges(6,3) * t231 - t229 * t232 - t230) * t248 + t272 + t296;
t84 = -t111 * t232 - t151 * t316;
t82 = t251 * t259 + (-t251 * rSges(4,3) + t248 * t274) * t248 + t300;
t81 = (t135 * t251 - t136 * t248) * t231;
t79 = -t241 * t316 + t281 + t343;
t78 = (t207 - t245) * t251 + (-t206 * t232 - t230 + (-rSges(7,3) + t241) * t231) * t248 + t271;
t74 = -t120 * t317 + t102;
t73 = -t111 * t317 + t100;
t70 = t256 * t251;
t69 = t256 * t248;
t60 = t130 * t316 + t132 * t197 + t134 * t198;
t59 = t129 * t316 + t131 * t197 + t133 * t198;
t58 = t130 * t317 + t132 * t195 + t134 * t196;
t57 = t129 * t317 + t131 * t195 + t133 * t196;
t54 = t303 * t232 + (-t147 - t155) * t316;
t53 = t88 + t302;
t52 = t135 * t248 + t136 * t251 + t279;
t37 = -t232 * t325 + t301 * t316;
t35 = t303 * t317 + t102 + t112;
t34 = -t317 * t325 + t326;
t33 = t119 * t248 + t120 * t251 + t257;
t32 = t289 * t232 + (-t147 + t301) * t316;
t31 = t36 + t302;
t30 = t248 * t60 - t251 * t59;
t29 = t248 * t58 - t251 * t57;
t27 = t289 * t317 + t112 + t326;
t19 = t325 * t251 + (t110 + t98) * t248 + t257;
t16 = -t232 * t77 + (t248 * t59 + t251 * t60) * t231;
t15 = -t232 * t76 + (t248 * t57 + t251 * t58) * t231;
t1 = [t250 * (Icges(3,2) * t250 + t324) + t247 * (Icges(3,1) * t247 + t323) + Icges(2,3) + t139 + t141 + t146 + (Icges(4,2) * t232 - t148 - t152 - t160 + t322) * t232 + (Icges(4,1) * t231 - t318 - t319 - t320 + t321) * t231 + m(7) * (t78 ^ 2 + t79 ^ 2) + m(6) * (t85 ^ 2 + t86 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t143 ^ 2 + t144 ^ 2) + m(3) * (t158 ^ 2 + t159 ^ 2) + m(2) * (t214 ^ 2 + t215 ^ 2); (-t55 / 0.2e1 - t48 / 0.2e1 - t65 / 0.2e1 - t71 / 0.2e1 + (-Icges(4,6) * t251 + t248 * t267) * t339 - t231 * (-Icges(4,5) * t251 + t248 * t269) / 0.2e1 - t250 * (-Icges(3,6) * t251 + t248 * t268) / 0.2e1 - t247 * (-Icges(3,5) * t251 + t248 * t270) / 0.2e1 + t285 * t251 - t291) * t251 + (t56 / 0.2e1 + t49 / 0.2e1 + t66 / 0.2e1 + t72 / 0.2e1 + (Icges(4,6) * t248 + t251 * t267) * t346 + (Icges(4,5) * t248 + t251 * t269) * t347 + (Icges(3,6) * t248 + t251 * t268) * t344 + (Icges(3,5) * t248 + t251 * t270) * t345 + t285 * t248 + t290) * t248 + m(5) * (t123 * t97 + t124 * t96) + m(6) * (t85 * t91 + t86 * t90) + m(7) * (t69 * t79 + t70 * t78) + m(4) * (t143 * t167 + t144 * t166) + m(3) * (-t158 * t251 - t159 * t248) * t213; m(7) * (t19 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t33 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(5) * (t123 ^ 2 + t124 ^ 2 + t52 ^ 2) + m(4) * (t166 ^ 2 + t167 ^ 2 + t82 ^ 2) + m(3) * (t213 ^ 2 * t292 + t138 ^ 2) + (t342 * t243 - t20 - t24 - t29) * t251 + (t21 + t25 + t30 + t341 * t242 + (t342 * t248 + t341 * t251) * t251) * t248; m(7) * (t248 * t78 - t251 * t79) + m(6) * (t248 * t85 - t251 * t86) + m(5) * (t248 * t96 - t251 * t97) + m(4) * (t143 * t248 - t144 * t251); m(7) * (t248 * t70 - t251 * t69) + m(6) * (t248 * t91 - t251 * t90) + m(5) * (-t123 * t251 + t124 * t248) + m(4) * (-t166 * t251 + t167 * t248); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t292; (-t87 + t331) * t232 + m(7) * (t31 * t78 + t32 * t79) + m(6) * (t53 * t85 + t54 * t86) + m(5) * (t92 * t96 + t93 * t97) + (t248 * t291 + t251 * t290) * t231 + t254; m(7) * (t19 * t27 + t31 * t70 + t32 * t69) + m(6) * (t33 * t35 + t53 * t91 + t54 * t90) + m(5) * (t123 * t93 + t124 * t92 + t52 * t81) + (t251 * t30 / 0.2e1 + t29 * t338) * t231 + t253 + t16 * t338 + (t62 * t248 - t61 * t251) * t339 + t15 * t337; m(5) * (t248 * t92 - t251 * t93) + m(6) * (t248 * t53 - t251 * t54) + m(7) * (t248 * t31 - t251 * t32); (t87 * t232 + t332) * t232 + (t248 * t15 + t251 * t16 - t232 * (t248 * t61 + t251 * t62)) * t231 + m(7) * (t27 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t35 ^ 2 + t53 ^ 2 + t54 ^ 2) + m(5) * (t81 ^ 2 + t92 ^ 2 + t93 ^ 2) + t288; t331 * t232 + m(7) * (t36 * t78 + t37 * t79) + m(6) * (t85 * t88 + t86 * t89) + t254; m(7) * (t19 * t34 + t36 * t70 + t37 * t69) + m(6) * (t33 * t74 + t88 * t91 + t89 * t90) + t253; m(6) * (t248 * t88 - t251 * t89) + m(7) * (t248 * t36 - t251 * t37); m(7) * (t27 * t34 + t31 * t36 + t32 * t37) + m(6) * (t35 * t74 + t53 * t88 + t54 * t89) + t255; m(7) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t74 ^ 2 + t88 ^ 2 + t89 ^ 2) + t255; m(7) * (t78 * t83 + t79 * t84) - t327 + t282; m(7) * (t19 * t73 + t69 * t84 + t70 * t83) + t278; m(7) * (t248 * t83 - t251 * t84); m(7) * (t27 * t73 + t31 * t83 + t32 * t84) + t280; m(7) * (t34 * t73 + t36 * t83 + t37 * t84) + t280; m(7) * (t73 ^ 2 + t83 ^ 2 + t84 ^ 2) + t280;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
