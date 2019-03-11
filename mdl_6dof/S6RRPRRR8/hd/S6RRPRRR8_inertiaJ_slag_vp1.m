% Calculate joint inertia matrix for
% S6RRPRRR8
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:22
% EndTime: 2019-03-09 14:02:34
% DurationCPUTime: 4.76s
% Computational Cost: add. (17914->570), mult. (15502->799), div. (0->0), fcn. (16670->12), ass. (0->264)
t266 = sin(qJ(1));
t343 = t266 / 0.2e1;
t268 = cos(qJ(1));
t340 = t268 / 0.2e1;
t258 = pkin(11) + qJ(4);
t249 = qJ(5) + t258;
t246 = qJ(6) + t249;
t230 = sin(t246);
t231 = cos(t246);
t267 = cos(qJ(2));
t318 = t267 * t268;
t174 = -t230 * t318 + t266 * t231;
t175 = t266 * t230 + t231 * t318;
t265 = sin(qJ(2));
t321 = t265 * t268;
t121 = t175 * rSges(7,1) + t174 * rSges(7,2) + rSges(7,3) * t321;
t263 = cos(pkin(11));
t245 = t263 * pkin(3) + pkin(2);
t248 = cos(t258);
t217 = pkin(4) * t248 + t245;
t244 = cos(t249);
t189 = pkin(5) * t244 + t217;
t247 = sin(t258);
t262 = sin(pkin(11));
t219 = t262 * pkin(3) + pkin(4) * t247;
t243 = sin(t249);
t204 = pkin(5) * t243 + t219;
t349 = t189 * t318 + t266 * t204 + t121;
t260 = t266 ^ 2;
t261 = t268 ^ 2;
t348 = m(4) / 0.2e1;
t347 = m(5) / 0.2e1;
t346 = m(6) / 0.2e1;
t345 = m(7) / 0.2e1;
t322 = t265 * t266;
t319 = t266 * t267;
t172 = -t230 * t319 - t231 * t268;
t173 = -t230 * t268 + t231 * t319;
t113 = Icges(7,5) * t173 + Icges(7,6) * t172 + Icges(7,3) * t322;
t115 = Icges(7,4) * t173 + Icges(7,2) * t172 + Icges(7,6) * t322;
t117 = Icges(7,1) * t173 + Icges(7,4) * t172 + Icges(7,5) * t322;
t38 = t113 * t322 + t115 * t172 + t117 * t173;
t114 = Icges(7,5) * t175 + Icges(7,6) * t174 + Icges(7,3) * t321;
t116 = Icges(7,4) * t175 + Icges(7,2) * t174 + Icges(7,6) * t321;
t118 = Icges(7,1) * t175 + Icges(7,4) * t174 + Icges(7,5) * t321;
t39 = t114 * t322 + t116 * t172 + t118 * t173;
t158 = -Icges(7,3) * t267 + (Icges(7,5) * t231 - Icges(7,6) * t230) * t265;
t159 = -Icges(7,6) * t267 + (Icges(7,4) * t231 - Icges(7,2) * t230) * t265;
t160 = -Icges(7,5) * t267 + (Icges(7,1) * t231 - Icges(7,4) * t230) * t265;
t67 = t158 * t322 + t159 * t172 + t160 * t173;
t5 = -t67 * t267 + (t266 * t38 + t268 * t39) * t265;
t40 = t113 * t321 + t174 * t115 + t175 * t117;
t41 = t114 * t321 + t174 * t116 + t175 * t118;
t68 = t158 * t321 + t174 * t159 + t175 * t160;
t6 = -t68 * t267 + (t266 * t40 + t268 * t41) * t265;
t344 = t6 * t321 + t5 * t322;
t342 = -t267 / 0.2e1;
t341 = -t268 / 0.2e1;
t226 = rSges(3,1) * t265 + rSges(3,2) * t267;
t339 = m(3) * t226;
t338 = pkin(2) * t267;
t337 = -pkin(2) + t245;
t264 = -pkin(8) - qJ(3);
t150 = t265 * t231 * t160;
t327 = t159 * t230;
t83 = -t267 * t158 - t265 * t327 + t150;
t333 = t83 * t267;
t51 = -t113 * t267 + (-t115 * t230 + t117 * t231) * t265;
t52 = -t114 * t267 + (-t116 * t230 + t118 * t231) * t265;
t14 = -t333 + (t266 * t51 + t268 * t52) * t265;
t183 = -t243 * t319 - t244 * t268;
t184 = -t243 * t268 + t244 * t319;
t123 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t322;
t125 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t322;
t127 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t322;
t59 = -t123 * t267 + (-t125 * t243 + t127 * t244) * t265;
t185 = -t243 * t318 + t266 * t244;
t186 = t266 * t243 + t244 * t318;
t124 = Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t321;
t126 = Icges(6,4) * t186 + Icges(6,2) * t185 + Icges(6,6) * t321;
t128 = Icges(6,1) * t186 + Icges(6,4) * t185 + Icges(6,5) * t321;
t60 = -t124 * t267 + (-t126 * t243 + t128 * t244) * t265;
t166 = -Icges(6,5) * t267 + (Icges(6,1) * t244 - Icges(6,4) * t243) * t265;
t154 = t265 * t244 * t166;
t164 = -Icges(6,3) * t267 + (Icges(6,5) * t244 - Icges(6,6) * t243) * t265;
t165 = -Icges(6,6) * t267 + (Icges(6,4) * t244 - Icges(6,2) * t243) * t265;
t326 = t165 * t243;
t90 = -t267 * t164 - t265 * t326 + t154;
t336 = -t14 + t90 * t267 - (t266 * t59 + t268 * t60) * t265;
t335 = -t83 - t90;
t334 = t268 * rSges(3,3);
t281 = -t173 * rSges(7,1) - t172 * rSges(7,2);
t120 = rSges(7,3) * t322 - t281;
t108 = t120 * t321;
t257 = -pkin(9) + t264;
t250 = -pkin(10) + t257;
t308 = t268 * t219 + t257 * t322;
t312 = t189 - t217;
t324 = t204 * t268;
t88 = -t324 + (-t250 * t265 + t312 * t267) * t266 + t308;
t332 = t88 * t321 + t108;
t303 = t250 - t257;
t310 = t217 * t318 + t266 * t219;
t331 = -t303 * t321 - t310 + t349;
t330 = Icges(3,4) * t265;
t329 = Icges(3,4) * t267;
t328 = qJ(3) * t265;
t177 = -Icges(5,6) * t267 + (Icges(5,4) * t248 - Icges(5,2) * t247) * t265;
t325 = t177 * t247;
t323 = t262 * t268;
t320 = t266 * t262;
t305 = pkin(3) * t323 + t264 * t322;
t307 = t217 - t245;
t104 = t307 * t319 + t305 - t308;
t301 = t257 - t264;
t151 = t307 * t265 + t301 * t267;
t317 = t267 * t104 + t151 * t322;
t306 = -pkin(3) * t320 - t245 * t318;
t105 = -t301 * t321 + t306 + t310;
t130 = t186 * rSges(6,1) + t185 * rSges(6,2) + rSges(6,3) * t321;
t316 = -t105 - t130;
t131 = t312 * t265 + t303 * t267;
t161 = -rSges(7,3) * t267 + (rSges(7,1) * t231 - rSges(7,2) * t230) * t265;
t315 = -t131 - t161;
t91 = t267 * t120 + t161 * t322;
t282 = -rSges(6,1) * t184 - rSges(6,2) * t183;
t129 = rSges(6,3) * t322 - t282;
t170 = -rSges(6,3) * t267 + (rSges(6,1) * t244 - rSges(6,2) * t243) * t265;
t94 = t267 * t129 + t170 * t322;
t225 = pkin(2) * t265 - qJ(3) * t267;
t313 = -(qJ(3) + t264) * t267 - t337 * t265 - t225;
t311 = rSges(4,3) * t267 - (rSges(4,1) * t263 - rSges(4,2) * t262) * t265 - t225;
t304 = pkin(2) * t318 + qJ(3) * t321;
t309 = t260 * (t328 + t338) + t268 * t304;
t302 = t268 * pkin(1) + t266 * pkin(7);
t300 = t260 + t261;
t194 = -t247 * t319 - t248 * t268;
t195 = -t247 * t268 + t248 * t319;
t132 = Icges(5,5) * t195 + Icges(5,6) * t194 + Icges(5,3) * t322;
t134 = Icges(5,4) * t195 + Icges(5,2) * t194 + Icges(5,6) * t322;
t136 = Icges(5,1) * t195 + Icges(5,4) * t194 + Icges(5,5) * t322;
t63 = -t132 * t267 + (-t134 * t247 + t136 * t248) * t265;
t176 = -Icges(5,3) * t267 + (Icges(5,5) * t248 - Icges(5,6) * t247) * t265;
t178 = -Icges(5,5) * t267 + (Icges(5,1) * t248 - Icges(5,4) * t247) * t265;
t78 = t176 * t322 + t177 * t194 + t178 * t195;
t299 = t63 / 0.2e1 + t78 / 0.2e1;
t196 = -t247 * t318 + t266 * t248;
t197 = t266 * t247 + t248 * t318;
t133 = Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t321;
t135 = Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t321;
t137 = Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t321;
t64 = -t133 * t267 + (-t135 * t247 + t137 * t248) * t265;
t79 = t176 * t321 + t196 * t177 + t197 * t178;
t298 = t64 / 0.2e1 + t79 / 0.2e1;
t297 = -t105 - t331;
t43 = t123 * t322 + t125 * t183 + t127 * t184;
t44 = t124 * t322 + t126 * t183 + t128 * t184;
t72 = t164 * t322 + t165 * t183 + t166 * t184;
t12 = -t72 * t267 + (t266 * t43 + t268 * t44) * t265;
t45 = t123 * t321 + t185 * t125 + t186 * t127;
t46 = t124 * t321 + t185 * t126 + t186 * t128;
t73 = t164 * t321 + t185 * t165 + t186 * t166;
t13 = -t73 * t267 + (t266 * t45 + t268 * t46) * t265;
t296 = t12 * t322 + t13 * t321 + t344;
t295 = -t151 + t313;
t181 = -rSges(5,3) * t267 + (rSges(5,1) * t248 - rSges(5,2) * t247) * t265;
t294 = -t181 + t313;
t139 = t197 * rSges(5,1) + t196 * rSges(5,2) + rSges(5,3) * t321;
t213 = -t262 * t318 + t266 * t263;
t214 = t263 * t318 + t320;
t293 = t214 * rSges(4,1) + t213 * rSges(4,2) + rSges(4,3) * t321;
t292 = t322 / 0.2e1;
t291 = t321 / 0.2e1;
t290 = (t51 + t67) * t292 + (t52 + t68) * t291;
t36 = t131 * t322 + t267 * t88 + t91;
t289 = -t267 * t14 + t344;
t273 = -t264 * t321 - t306;
t288 = t266 * ((t337 * t267 - t328) * t266 - t305) + t268 * (t273 - t304) + t309;
t287 = -t170 + t295;
t21 = t39 * t266 - t268 * t38;
t22 = t41 * t266 - t268 * t40;
t286 = t21 * t292 + t22 * t291 + t5 * t341 + t6 * t343 + (t52 * t266 - t51 * t268) * t342;
t285 = rSges(3,1) * t267 - rSges(3,2) * t265;
t211 = -t262 * t319 - t263 * t268;
t212 = t263 * t319 - t323;
t284 = -rSges(4,1) * t212 - t211 * rSges(4,2);
t283 = -rSges(5,1) * t195 - rSges(5,2) * t194;
t280 = t295 + t315;
t279 = Icges(3,1) * t267 - t330;
t278 = -Icges(3,2) * t265 + t329;
t277 = Icges(3,5) * t267 - Icges(3,6) * t265;
t274 = rSges(3,1) * t318 - rSges(3,2) * t321 + t266 * rSges(3,3);
t272 = t266 * t104 + t268 * t105 + t288;
t271 = t336 * t267 + t296;
t270 = t290 + (t59 + t72) * t292 + (t60 + t73) * t291;
t25 = t44 * t266 - t268 * t43;
t26 = t46 * t266 - t268 * t45;
t269 = t12 * t341 + t13 * t343 + t25 * t292 + t26 * t291 + t286 + (t60 * t266 - t59 * t268) * t342;
t255 = t268 * pkin(7);
t228 = rSges(2,1) * t268 - t266 * rSges(2,2);
t227 = -t266 * rSges(2,1) - rSges(2,2) * t268;
t221 = Icges(3,5) * t265 + Icges(3,6) * t267;
t199 = Icges(3,3) * t266 + t277 * t268;
t198 = -Icges(3,3) * t268 + t277 * t266;
t192 = -Icges(4,5) * t267 + (Icges(4,1) * t263 - Icges(4,4) * t262) * t265;
t191 = -Icges(4,6) * t267 + (Icges(4,4) * t263 - Icges(4,2) * t262) * t265;
t163 = t274 + t302;
t162 = t334 + t255 + (-pkin(1) - t285) * t266;
t156 = t265 * t248 * t178;
t153 = t311 * t268;
t152 = t311 * t266;
t149 = Icges(4,1) * t214 + Icges(4,4) * t213 + Icges(4,5) * t321;
t148 = Icges(4,1) * t212 + Icges(4,4) * t211 + Icges(4,5) * t322;
t147 = Icges(4,4) * t214 + Icges(4,2) * t213 + Icges(4,6) * t321;
t146 = Icges(4,4) * t212 + Icges(4,2) * t211 + Icges(4,6) * t322;
t145 = Icges(4,5) * t214 + Icges(4,6) * t213 + Icges(4,3) * t321;
t144 = Icges(4,5) * t212 + Icges(4,6) * t211 + Icges(4,3) * t322;
t142 = t268 * t274 + (t285 * t266 - t334) * t266;
t138 = rSges(5,3) * t322 - t283;
t112 = t129 * t321;
t110 = t293 + t302 + t304;
t109 = t255 + (-t338 - pkin(1) + (-rSges(4,3) - qJ(3)) * t265) * t266 + t284;
t107 = t294 * t268;
t106 = t294 * t266;
t100 = t104 * t321;
t99 = -t267 * t139 - t181 * t321;
t98 = t138 * t267 + t181 * t322;
t97 = t273 + t139 + t302;
t96 = t255 + (-rSges(5,3) * t265 - t245 * t267 - pkin(1)) * t266 + t283 + t305;
t95 = -t267 * t130 - t170 * t321;
t93 = -t267 * t176 - t265 * t325 + t156;
t92 = -t267 * t121 - t161 * t321;
t86 = -t257 * t321 + t130 + t302 + t310;
t85 = t255 + (-rSges(6,3) * t265 - t217 * t267 - pkin(1)) * t266 + t282 + t308;
t82 = t287 * t268;
t81 = t287 * t266;
t80 = (t138 * t268 - t139 * t266) * t265;
t77 = -t130 * t322 + t112;
t76 = -t250 * t321 + t302 + t349;
t75 = t324 + t255 + (-t189 * t267 - pkin(1) + (-rSges(7,3) + t250) * t265) * t266 + t281;
t74 = t266 * (rSges(4,3) * t322 - t284) + t268 * t293 + t309;
t71 = -t121 * t322 + t108;
t62 = t280 * t268;
t61 = t280 * t266;
t56 = t133 * t321 + t196 * t135 + t197 * t137;
t55 = t132 * t321 + t196 * t134 + t197 * t136;
t54 = t133 * t322 + t135 * t194 + t137 * t195;
t53 = t132 * t322 + t134 * t194 + t136 * t195;
t50 = t316 * t267 + (-t151 - t170) * t321;
t49 = t94 + t317;
t42 = t266 * t138 + t139 * t268 + t288;
t37 = -t267 * t331 + t315 * t321;
t35 = t316 * t322 + t100 + t112;
t34 = -t322 * t331 + t332;
t33 = t297 * t267 + (-t151 + t315) * t321;
t32 = t36 + t317;
t30 = t266 * t129 + t130 * t268 + t272;
t29 = t56 * t266 - t268 * t55;
t28 = t54 * t266 - t268 * t53;
t17 = t297 * t322 + t100 + t332;
t16 = -t79 * t267 + (t266 * t55 + t268 * t56) * t265;
t15 = -t78 * t267 + (t266 * t53 + t268 * t54) * t265;
t7 = t331 * t268 + (t120 + t88) * t266 + t272;
t1 = [Icges(2,3) + t150 + t154 + t156 + (-t158 - t164 - t176 + t330 - (Icges(4,5) * t263 - Icges(4,6) * t262) * t265 + (Icges(3,2) + Icges(4,3)) * t267) * t267 + (Icges(3,1) * t265 - t191 * t262 + t192 * t263 - t325 - t326 - t327 + t329) * t265 + m(7) * (t75 ^ 2 + t76 ^ 2) + m(6) * (t85 ^ 2 + t86 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t109 ^ 2 + t110 ^ 2) + m(3) * (t162 ^ 2 + t163 ^ 2) + m(2) * (t227 ^ 2 + t228 ^ 2); m(7) * (t61 * t76 + t62 * t75) + m(6) * (t81 * t86 + t82 * t85) + m(5) * (t106 * t97 + t107 * t96) + m(4) * (t109 * t153 + t110 * t152) + (-t67 / 0.2e1 - t72 / 0.2e1 - t211 * t191 / 0.2e1 - t212 * t192 / 0.2e1 - t59 / 0.2e1 - t51 / 0.2e1 - t162 * t339 + t221 * t340 + (t144 / 0.2e1 + Icges(3,6) * t340 - t278 * t266 / 0.2e1) * t267 - t299) * t268 + (t213 * t191 / 0.2e1 + t214 * t192 / 0.2e1 + t60 / 0.2e1 + t52 / 0.2e1 + t68 / 0.2e1 + t73 / 0.2e1 - t163 * t339 + t221 * t343 + (-t145 / 0.2e1 + Icges(3,6) * t343 + t278 * t340) * t267 + t298) * t266 + ((Icges(3,5) * t266 - t147 * t262 + t149 * t263 + t279 * t268) * t343 + (-Icges(3,5) * t268 - t146 * t262 + t148 * t263 + t279 * t266) * t341) * t265; m(7) * (t61 ^ 2 + t62 ^ 2 + t7 ^ 2) + m(6) * (t30 ^ 2 + t81 ^ 2 + t82 ^ 2) + m(5) * (t106 ^ 2 + t107 ^ 2 + t42 ^ 2) + m(4) * (t152 ^ 2 + t153 ^ 2 + t74 ^ 2) + m(3) * (t226 ^ 2 * t300 + t142 ^ 2) + (-t261 * t198 - t21 - t25 - t28 + (t144 * t322 + t211 * t146 + t212 * t148) * t268) * t268 + (t22 + t26 + t29 + t260 * t199 + (t145 * t321 + t213 * t147 + t214 * t149) * t266 + (-t144 * t321 - t145 * t322 - t213 * t146 - t147 * t211 - t214 * t148 - t149 * t212 - t266 * t198 + t268 * t199) * t268) * t266; 0.2e1 * ((t266 * t76 + t268 * t75) * t345 + (t266 * t86 + t268 * t85) * t346 + (t266 * t97 + t268 * t96) * t347 + (t109 * t268 + t110 * t266) * t348) * t265; m(7) * (-t267 * t7 + (t266 * t61 + t268 * t62) * t265) + m(6) * (-t267 * t30 + (t266 * t81 + t268 * t82) * t265) + m(5) * (-t267 * t42 + (t106 * t266 + t107 * t268) * t265) + m(4) * (-t267 * t74 + (t152 * t266 + t153 * t268) * t265); 0.2e1 * (t348 + t347 + t346 + t345) * (t265 ^ 2 * t300 + t267 ^ 2); (-t93 + t335) * t267 + m(7) * (t32 * t75 + t33 * t76) + m(6) * (t49 * t85 + t50 * t86) + m(5) * (t96 * t98 + t97 * t99) + (t266 * t299 + t268 * t298) * t265 + t270; t269 + t15 * t341 + (t64 * t266 - t63 * t268) * t342 + m(7) * (t17 * t7 + t32 * t62 + t33 * t61) + m(6) * (t30 * t35 + t49 * t82 + t50 * t81) + m(5) * (t106 * t99 + t107 * t98 + t42 * t80) + (t28 * t343 + t29 * t340) * t265 + t16 * t343; m(5) * (-t80 * t267 + (t266 * t99 + t268 * t98) * t265) + m(6) * (-t35 * t267 + (t266 * t50 + t268 * t49) * t265) + m(7) * (-t17 * t267 + (t266 * t33 + t268 * t32) * t265); (t93 * t267 + t336) * t267 + (t268 * t16 + t266 * t15 - t267 * (t266 * t63 + t268 * t64)) * t265 + m(7) * (t17 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(5) * (t80 ^ 2 + t98 ^ 2 + t99 ^ 2) + t296; t335 * t267 + m(7) * (t36 * t75 + t37 * t76) + m(6) * (t85 * t94 + t86 * t95) + t270; m(7) * (t34 * t7 + t36 * t62 + t37 * t61) + m(6) * (t30 * t77 + t81 * t95 + t82 * t94) + t269; m(6) * (-t77 * t267 + (t266 * t95 + t268 * t94) * t265) + m(7) * (-t34 * t267 + (t266 * t37 + t268 * t36) * t265); m(7) * (t17 * t34 + t32 * t36 + t33 * t37) + m(6) * (t35 * t77 + t49 * t94 + t50 * t95) + t271; m(7) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t77 ^ 2 + t94 ^ 2 + t95 ^ 2) + t271; m(7) * (t75 * t91 + t76 * t92) - t333 + t290; m(7) * (t61 * t92 + t62 * t91 + t7 * t71) + t286; m(7) * (-t71 * t267 + (t266 * t92 + t268 * t91) * t265); m(7) * (t17 * t71 + t32 * t91 + t33 * t92) + t289; m(7) * (t34 * t71 + t36 * t91 + t37 * t92) + t289; m(7) * (t71 ^ 2 + t91 ^ 2 + t92 ^ 2) + t289;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
