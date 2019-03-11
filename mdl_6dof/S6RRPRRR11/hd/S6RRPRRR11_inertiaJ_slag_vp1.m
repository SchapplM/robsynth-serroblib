% Calculate joint inertia matrix for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:27
% EndTime: 2019-03-09 14:29:37
% DurationCPUTime: 4.77s
% Computational Cost: add. (11891->536), mult. (14230->741), div. (0->0), fcn. (15330->10), ass. (0->252)
t338 = Icges(4,1) + Icges(3,3);
t247 = sin(qJ(2));
t250 = cos(qJ(2));
t337 = (-Icges(4,4) + Icges(3,5)) * t250 + (Icges(4,5) - Icges(3,6)) * t247;
t248 = sin(qJ(1));
t336 = -t248 / 0.2e1;
t327 = t248 / 0.2e1;
t251 = cos(qJ(1));
t326 = -t251 / 0.2e1;
t335 = t251 / 0.2e1;
t328 = t247 / 0.2e1;
t334 = -t337 * t248 + t251 * t338;
t333 = t248 * t338 + t337 * t251;
t243 = t248 ^ 2;
t244 = t251 ^ 2;
t332 = m(4) / 0.2e1;
t331 = m(5) / 0.2e1;
t330 = m(6) / 0.2e1;
t329 = m(7) / 0.2e1;
t252 = -pkin(9) - pkin(8);
t246 = sin(qJ(4));
t325 = pkin(4) * t246;
t249 = cos(qJ(4));
t229 = t249 * pkin(4) + pkin(3);
t324 = t251 * rSges(4,1);
t323 = t251 * rSges(3,3);
t245 = qJ(4) + qJ(5);
t232 = qJ(6) + t245;
t227 = sin(t232);
t228 = cos(t232);
t313 = t247 * t248;
t163 = t227 * t251 + t228 * t313;
t164 = t227 * t313 - t228 * t251;
t274 = -t164 * rSges(7,1) - t163 * rSges(7,2);
t310 = t248 * t250;
t110 = rSges(7,3) * t310 - t274;
t308 = t250 * t251;
t100 = t110 * t308;
t241 = -pkin(10) + t252;
t230 = sin(t245);
t202 = pkin(5) * t230 + t325;
t286 = -t202 + t325;
t307 = t250 * t252;
t299 = t251 * t229 + t248 * t307;
t231 = cos(t245);
t201 = pkin(5) * t231 + t229;
t306 = t251 * t201;
t97 = -t306 + (-t241 * t250 - t247 * t286) * t248 + t299;
t322 = t97 * t308 + t100;
t312 = t247 * t251;
t161 = -t227 * t248 + t228 * t312;
t162 = t227 * t312 + t228 * t248;
t109 = t162 * rSges(7,1) + t161 * rSges(7,2) + rSges(7,3) * t308;
t101 = t247 * t109;
t295 = -t241 + t252;
t290 = t246 * t312;
t300 = -pkin(4) * t290 - t248 * t229;
t302 = t248 * t201 + t202 * t312;
t96 = t295 * t308 + t300 + t302;
t321 = t247 * t96 + t101;
t320 = t109 + t96;
t319 = t110 + t97;
t318 = Icges(3,4) * t247;
t317 = Icges(3,4) * t250;
t316 = Icges(4,6) * t247;
t315 = Icges(4,6) * t250;
t314 = qJ(3) * t247;
t311 = t248 * t249;
t309 = t249 * t251;
t140 = t247 * t295 + t250 * t286;
t152 = t247 * rSges(7,3) + (-rSges(7,1) * t227 - rSges(7,2) * t228) * t250;
t141 = t152 * t310;
t305 = t140 * t310 + t141;
t304 = -t140 - t152;
t298 = pkin(2) * t308 + qJ(3) * t312;
t303 = t243 * (pkin(2) * t250 + t314) + t251 * t298;
t210 = pkin(2) * t247 - qJ(3) * t250;
t301 = rSges(4,2) * t247 + rSges(4,3) * t250 - t210;
t297 = t251 * pkin(1) + t248 * pkin(7);
t296 = t248 * pkin(3) + pkin(8) * t308;
t294 = t243 + t244;
t103 = Icges(7,5) * t162 + Icges(7,6) * t161 + Icges(7,3) * t308;
t105 = Icges(7,4) * t162 + Icges(7,2) * t161 + Icges(7,6) * t308;
t107 = Icges(7,1) * t162 + Icges(7,4) * t161 + Icges(7,5) * t308;
t48 = t247 * t103 + (-t105 * t228 - t107 * t227) * t250;
t104 = Icges(7,5) * t164 + Icges(7,6) * t163 + Icges(7,3) * t310;
t106 = Icges(7,4) * t164 + Icges(7,2) * t163 + Icges(7,6) * t310;
t108 = Icges(7,1) * t164 + Icges(7,4) * t163 + Icges(7,5) * t310;
t49 = t247 * t104 + (-t106 * t228 - t108 * t227) * t250;
t35 = t103 * t308 + t105 * t161 + t107 * t162;
t36 = t104 * t308 + t106 * t161 + t108 * t162;
t149 = Icges(7,3) * t247 + (-Icges(7,5) * t227 - Icges(7,6) * t228) * t250;
t150 = Icges(7,6) * t247 + (-Icges(7,4) * t227 - Icges(7,2) * t228) * t250;
t151 = Icges(7,5) * t247 + (-Icges(7,1) * t227 - Icges(7,4) * t228) * t250;
t65 = t149 * t308 + t150 * t161 + t151 * t162;
t5 = t65 * t247 + (t248 * t36 + t251 * t35) * t250;
t37 = t103 * t310 + t105 * t163 + t107 * t164;
t38 = t104 * t310 + t106 * t163 + t108 * t164;
t66 = t149 * t310 + t150 * t163 + t151 * t164;
t6 = t66 * t247 + (t248 * t38 + t251 * t37) * t250;
t142 = t247 * t149;
t267 = -t150 * t228 - t151 * t227;
t77 = (t250 * t267 + t142) * t247;
t293 = t5 * t308 + t6 * t310 + t247 * (t77 + (t248 * t49 + t251 * t48) * t250);
t194 = -t246 * t248 + t247 * t309;
t195 = t290 + t311;
t124 = Icges(5,5) * t195 + Icges(5,6) * t194 + Icges(5,3) * t308;
t126 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t308;
t128 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t308;
t61 = t247 * t124 + (-t126 * t249 - t128 * t246) * t250;
t167 = Icges(5,3) * t247 + (-Icges(5,5) * t246 - Icges(5,6) * t249) * t250;
t170 = Icges(5,6) * t247 + (-Icges(5,4) * t246 - Icges(5,2) * t249) * t250;
t173 = Icges(5,5) * t247 + (-Icges(5,1) * t246 - Icges(5,4) * t249) * t250;
t80 = t167 * t308 + t170 * t194 + t173 * t195;
t292 = t61 / 0.2e1 + t80 / 0.2e1;
t196 = t246 * t251 + t247 * t311;
t197 = t246 * t313 - t309;
t125 = Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t310;
t127 = Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t310;
t129 = Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t310;
t62 = t247 * t125 + (-t127 * t249 - t129 * t246) * t250;
t81 = t167 * t310 + t170 * t196 + t173 * t197;
t291 = t81 / 0.2e1 + t62 / 0.2e1;
t183 = -t230 * t248 + t231 * t312;
t184 = t230 * t312 + t231 * t248;
t118 = t184 * rSges(6,1) + t183 * rSges(6,2) + rSges(6,3) * t308;
t133 = t195 * rSges(5,1) + t194 * rSges(5,2) + rSges(5,3) * t308;
t289 = t310 / 0.2e1;
t288 = t308 / 0.2e1;
t287 = -Icges(4,4) * t247 / 0.2e1 + Icges(3,5) * t328 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t250;
t285 = -t247 * pkin(8) - t210;
t239 = t251 * pkin(3);
t284 = t248 * (pkin(8) * t310 - t239) + t251 * t296 + t303;
t283 = t297 + t298;
t19 = t248 * t35 - t251 * t36;
t20 = t248 * t37 - t251 * t38;
t282 = t19 * t288 + t20 * t289 + t6 * t326 + t5 * t327 + (t48 * t248 - t49 * t251) * t328;
t281 = t77 + (t49 + t66) * t289 + (t48 + t65) * t288;
t182 = t247 * rSges(5,3) + (-rSges(5,1) * t246 - rSges(5,2) * t249) * t250;
t280 = -t182 + t285;
t189 = -t250 * t325 + (-pkin(8) - t252) * t247;
t279 = -t189 + t285;
t112 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t308;
t114 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t308;
t116 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t308;
t42 = t112 * t308 + t114 * t183 + t116 * t184;
t185 = t230 * t251 + t231 * t313;
t186 = t230 * t313 - t231 * t251;
t113 = Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t310;
t115 = Icges(6,4) * t186 + Icges(6,2) * t185 + Icges(6,6) * t310;
t117 = Icges(6,1) * t186 + Icges(6,4) * t185 + Icges(6,5) * t310;
t43 = t113 * t308 + t115 * t183 + t117 * t184;
t153 = Icges(6,3) * t247 + (-Icges(6,5) * t230 - Icges(6,6) * t231) * t250;
t154 = Icges(6,6) * t247 + (-Icges(6,4) * t230 - Icges(6,2) * t231) * t250;
t155 = Icges(6,5) * t247 + (-Icges(6,1) * t230 - Icges(6,4) * t231) * t250;
t72 = t153 * t308 + t154 * t183 + t155 * t184;
t11 = t72 * t247 + (t248 * t43 + t251 * t42) * t250;
t44 = t112 * t310 + t114 * t185 + t116 * t186;
t45 = t113 * t310 + t115 * t185 + t117 * t186;
t73 = t153 * t310 + t154 * t185 + t155 * t186;
t12 = t73 * t247 + (t248 * t45 + t251 * t44) * t250;
t52 = t247 * t112 + (-t114 * t231 - t116 * t230) * t250;
t53 = t247 * t113 + (-t115 * t231 - t117 * t230) * t250;
t146 = t247 * t153;
t266 = -t154 * t231 - t155 * t230;
t83 = (t250 * t266 + t146) * t247;
t278 = t11 * t308 + t12 * t310 + t247 * (t83 + (t248 * t53 + t251 * t52) * t250) + t293;
t277 = rSges(3,1) * t250 - rSges(3,2) * t247;
t276 = -t197 * rSges(5,1) - t196 * rSges(5,2);
t275 = -t186 * rSges(6,1) - t185 * rSges(6,2);
t273 = Icges(3,1) * t250 - t318;
t272 = -Icges(3,2) * t247 + t317;
t269 = -Icges(4,2) * t250 + t316;
t268 = Icges(4,3) * t247 - t315;
t265 = -t170 * t249 - t173 * t246;
t159 = t247 * rSges(6,3) + (-rSges(6,1) * t230 - rSges(6,2) * t231) * t250;
t260 = -t159 + t279;
t259 = rSges(3,1) * t308 - rSges(3,2) * t312 + t248 * rSges(3,3);
t258 = t248 * rSges(4,1) - rSges(4,2) * t308 + rSges(4,3) * t312;
t257 = -t251 * t307 - t300;
t138 = t257 - t296;
t139 = t239 + (-pkin(8) * t250 + t247 * t325) * t248 - t299;
t256 = t251 * t138 + t248 * t139 + t284;
t255 = t279 + t304;
t24 = t248 * t42 - t251 * t43;
t25 = t248 * t44 - t251 * t45;
t254 = t11 * t327 + t12 * t326 + t24 * t288 + t25 * t289 + (t52 * t248 - t53 * t251) * t328 + t282;
t253 = t281 + t83 + (t53 + t73) * t289 + (t52 + t72) * t288;
t238 = t251 * pkin(7);
t214 = rSges(2,1) * t251 - rSges(2,2) * t248;
t213 = -rSges(2,1) * t248 - rSges(2,2) * t251;
t212 = rSges(3,1) * t247 + rSges(3,2) * t250;
t160 = t247 * t167;
t156 = t189 * t310;
t148 = t301 * t251;
t147 = t301 * t248;
t145 = t259 + t297;
t144 = t323 + t238 + (-pkin(1) - t277) * t248;
t143 = t159 * t310;
t136 = t258 + t283;
t135 = t324 + t238 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t250 + (-rSges(4,3) - qJ(3)) * t247) * t248;
t134 = rSges(5,3) * t310 - t276;
t130 = t247 * t138;
t123 = t280 * t251;
t122 = t280 * t248;
t121 = t139 * t308;
t120 = t251 * t259 + (t248 * t277 - t323) * t248;
t119 = rSges(6,3) * t310 - t275;
t111 = t247 * t118;
t102 = t119 * t308;
t99 = t260 * t251;
t98 = t260 * t248;
t95 = t133 * t247 - t182 * t308;
t94 = -t134 * t247 + t182 * t310;
t91 = t283 + t133 + t296;
t90 = t238 + t239 + (-t314 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(8)) * t250) * t248 + t276;
t89 = (t250 * t265 + t160) * t247;
t88 = t251 * t258 + (-t324 + (-rSges(4,2) * t250 + rSges(4,3) * t247) * t248) * t248 + t303;
t87 = -t159 * t308 + t111;
t86 = -t119 * t247 + t143;
t85 = -t152 * t308 + t101;
t84 = -t110 * t247 + t141;
t82 = (-t133 * t248 + t134 * t251) * t250;
t79 = t257 + t283 + t118;
t78 = t238 + (-pkin(1) + (-rSges(6,3) - pkin(2)) * t250 + (-qJ(3) - t325) * t247) * t248 + t275 + t299;
t76 = t255 * t251;
t75 = t255 * t248;
t74 = -t118 * t310 + t102;
t71 = -t241 * t308 + t109 + t283 + t302;
t70 = t306 + t238 + (-pkin(1) + (-qJ(3) - t202) * t247 + (-rSges(7,3) - pkin(2) + t241) * t250) * t248 + t274;
t69 = -t109 * t310 + t100;
t60 = t111 + t130 + (-t159 - t189) * t308;
t59 = t143 + t156 + (-t119 - t139) * t247;
t58 = t133 * t251 + t134 * t248 + t284;
t57 = t125 * t310 + t127 * t196 + t129 * t197;
t56 = t124 * t310 + t126 * t196 + t128 * t197;
t55 = t125 * t308 + t127 * t194 + t129 * t195;
t54 = t124 * t308 + t126 * t194 + t128 * t195;
t41 = t102 + t121 + (-t118 - t138) * t310;
t40 = t304 * t308 + t321;
t39 = -t247 * t319 + t305;
t34 = t118 * t251 + t119 * t248 + t256;
t33 = t130 + (-t189 + t304) * t308 + t321;
t32 = t156 + (-t139 - t319) * t247 + t305;
t31 = -t310 * t320 + t322;
t30 = t248 * t56 - t251 * t57;
t29 = t248 * t54 - t251 * t55;
t27 = t121 + (-t138 - t320) * t310 + t322;
t21 = t248 * t319 + t251 * t320 + t256;
t18 = t81 * t247 + (t248 * t57 + t251 * t56) * t250;
t17 = t80 * t247 + (t248 * t55 + t251 * t54) * t250;
t1 = [Icges(2,3) + t142 + t146 + t160 + m(7) * (t70 ^ 2 + t71 ^ 2) + m(6) * (t78 ^ 2 + t79 ^ 2) + m(5) * (t90 ^ 2 + t91 ^ 2) + m(4) * (t135 ^ 2 + t136 ^ 2) + m(3) * (t144 ^ 2 + t145 ^ 2) + m(2) * (t213 ^ 2 + t214 ^ 2) + (t265 + t266 + t267 + t316 + t318 + (Icges(3,2) + Icges(4,3)) * t250) * t250 + (t315 + t317 + (Icges(3,1) + Icges(4,2)) * t247) * t247; (-t53 / 0.2e1 - t49 / 0.2e1 - t73 / 0.2e1 - t66 / 0.2e1 + t287 * t251 + (Icges(4,5) * t326 + Icges(3,6) * t335 + t268 * t327 + t272 * t336) * t250 + (Icges(4,4) * t326 + Icges(3,5) * t335 + t269 * t327 + t273 * t336) * t247 - t291) * t251 + (t52 / 0.2e1 + t48 / 0.2e1 + t72 / 0.2e1 + t65 / 0.2e1 + t287 * t248 + (Icges(4,5) * t336 + Icges(3,6) * t327 + t268 * t326 + t272 * t335) * t250 + (Icges(4,4) * t336 + Icges(3,5) * t327 + t269 * t326 + t273 * t335) * t247 + t292) * t248 + m(4) * (t135 * t148 + t136 * t147) + m(5) * (t122 * t91 + t123 * t90) + m(6) * (t78 * t99 + t79 * t98) + m(7) * (t70 * t76 + t71 * t75) + m(3) * (-t144 * t251 - t145 * t248) * t212; m(7) * (t21 ^ 2 + t75 ^ 2 + t76 ^ 2) + m(6) * (t34 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(5) * (t122 ^ 2 + t123 ^ 2 + t58 ^ 2) + m(4) * (t147 ^ 2 + t148 ^ 2 + t88 ^ 2) + m(3) * (t212 ^ 2 * t294 + t120 ^ 2) + (t334 * t244 - t20 - t25 - t30) * t251 + (t19 + t24 + t29 + t333 * t243 + (t334 * t248 + t333 * t251) * t251) * t248; 0.2e1 * ((t248 * t71 + t251 * t70) * t329 + (t248 * t79 + t251 * t78) * t330 + (t248 * t91 + t251 * t90) * t331 + (t135 * t251 + t136 * t248) * t332) * t247; m(7) * (-t250 * t21 + (t248 * t75 + t251 * t76) * t247) + m(6) * (-t250 * t34 + (t248 * t98 + t251 * t99) * t247) + m(5) * (-t250 * t58 + (t122 * t248 + t123 * t251) * t247) + m(4) * (-t250 * t88 + (t147 * t248 + t148 * t251) * t247); 0.2e1 * (t332 + t331 + t330 + t329) * (t247 ^ 2 * t294 + t250 ^ 2); (t248 * t291 + t251 * t292) * t250 + m(7) * (t32 * t70 + t33 * t71) + m(6) * (t59 * t78 + t60 * t79) + m(5) * (t90 * t94 + t91 * t95) + t253 + t89; t254 + m(7) * (t21 * t27 + t32 * t76 + t33 * t75) + m(6) * (t34 * t41 + t59 * t99 + t60 * t98) + m(5) * (t122 * t95 + t123 * t94 + t58 * t82) + (t29 * t335 + t30 * t327) * t250 + t17 * t327 + (t61 * t248 - t62 * t251) * t328 + t18 * t326; m(5) * (-t82 * t250 + (t248 * t95 + t251 * t94) * t247) + m(6) * (-t41 * t250 + (t248 * t60 + t251 * t59) * t247) + m(7) * (-t27 * t250 + (t248 * t33 + t251 * t32) * t247); t247 * t89 + (t251 * t17 + t248 * t18 + t247 * (t248 * t62 + t251 * t61)) * t250 + m(7) * (t27 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t41 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t82 ^ 2 + t94 ^ 2 + t95 ^ 2) + t278; m(7) * (t39 * t70 + t40 * t71) + m(6) * (t78 * t86 + t79 * t87) + t253; m(7) * (t21 * t31 + t39 * t76 + t40 * t75) + m(6) * (t34 * t74 + t86 * t99 + t87 * t98) + t254; m(6) * (-t74 * t250 + (t248 * t87 + t251 * t86) * t247) + m(7) * (-t31 * t250 + (t248 * t40 + t251 * t39) * t247); m(7) * (t27 * t31 + t32 * t39 + t33 * t40) + m(6) * (t41 * t74 + t59 * t86 + t60 * t87) + t278; m(7) * (t31 ^ 2 + t39 ^ 2 + t40 ^ 2) + m(6) * (t74 ^ 2 + t86 ^ 2 + t87 ^ 2) + t278; m(7) * (t70 * t84 + t71 * t85) + t281; m(7) * (t21 * t69 + t75 * t85 + t76 * t84) + t282; m(7) * (-t69 * t250 + (t248 * t85 + t251 * t84) * t247); m(7) * (t27 * t69 + t32 * t84 + t33 * t85) + t293; m(7) * (t31 * t69 + t39 * t84 + t40 * t85) + t293; m(7) * (t69 ^ 2 + t84 ^ 2 + t85 ^ 2) + t293;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
