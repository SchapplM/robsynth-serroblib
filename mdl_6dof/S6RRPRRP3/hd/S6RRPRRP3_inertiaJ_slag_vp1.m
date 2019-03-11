% Calculate joint inertia matrix for
% S6RRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 11:51
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP3_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:47:22
% EndTime: 2019-03-09 11:47:36
% DurationCPUTime: 5.76s
% Computational Cost: add. (12488->481), mult. (11813->674), div. (0->0), fcn. (12644->10), ass. (0->238)
t232 = qJ(2) + pkin(10);
t224 = cos(t232);
t235 = qJ(4) + qJ(5);
t226 = cos(t235);
t242 = cos(qJ(1));
t298 = t226 * t242;
t225 = sin(t235);
t239 = sin(qJ(1));
t301 = t225 * t239;
t175 = -t224 * t301 - t298;
t299 = t226 * t239;
t300 = t225 * t242;
t176 = t224 * t299 - t300;
t223 = sin(t232);
t305 = t223 * t239;
t102 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t305;
t104 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t305;
t359 = t102 + t104;
t177 = -t224 * t300 + t299;
t178 = t224 * t298 + t301;
t304 = t223 * t242;
t103 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t304;
t105 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t304;
t358 = t103 + t105;
t106 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t305;
t108 = Icges(6,4) * t176 + Icges(6,2) * t175 + Icges(6,6) * t305;
t357 = t106 + t108;
t107 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t304;
t109 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t304;
t356 = t107 + t109;
t110 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t305;
t112 = Icges(6,1) * t176 + Icges(6,4) * t175 + Icges(6,5) * t305;
t355 = t110 + t112;
t111 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t304;
t113 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t304;
t354 = t111 + t113;
t353 = t357 * t175 + t355 * t176 + t305 * t359;
t352 = t175 * t356 + t176 * t354 + t305 * t358;
t351 = t357 * t177 + t355 * t178 + t304 * t359;
t350 = t177 * t356 + t178 * t354 + t304 * t358;
t145 = -Icges(7,3) * t224 + (Icges(7,5) * t226 - Icges(7,6) * t225) * t223;
t147 = -Icges(7,6) * t224 + (Icges(7,4) * t226 - Icges(7,2) * t225) * t223;
t149 = -Icges(7,5) * t224 + (Icges(7,1) * t226 - Icges(7,4) * t225) * t223;
t69 = t145 * t305 + t147 * t175 + t149 * t176;
t146 = -Icges(6,3) * t224 + (Icges(6,5) * t226 - Icges(6,6) * t225) * t223;
t148 = -Icges(6,6) * t224 + (Icges(6,4) * t226 - Icges(6,2) * t225) * t223;
t150 = -Icges(6,5) * t224 + (Icges(6,1) * t226 - Icges(6,4) * t225) * t223;
t70 = t146 * t305 + t148 * t175 + t150 * t176;
t349 = -t70 - t69;
t71 = t145 * t304 + t147 * t177 + t149 * t178;
t72 = t146 * t304 + t148 * t177 + t150 * t178;
t348 = -t71 - t72;
t343 = (t149 + t150) * t223 * t226;
t345 = (-t147 - t148) * t225;
t346 = -t145 - t146;
t335 = t223 * t345 + t346 * t224 + t343;
t347 = t335 * t224;
t344 = Icges(3,3) + Icges(4,3);
t238 = sin(qJ(2));
t241 = cos(qJ(2));
t342 = Icges(3,5) * t241 + Icges(4,5) * t224 - Icges(3,6) * t238 - Icges(4,6) * t223;
t341 = t349 * t224 + (t353 * t239 + t352 * t242) * t223;
t340 = t348 * t224 + (t351 * t239 + t350 * t242) * t223;
t339 = t352 * t239 - t353 * t242;
t338 = t350 * t239 - t351 * t242;
t53 = -t224 * t102 + (-t106 * t225 + t110 * t226) * t223;
t55 = -t224 * t104 + (-t108 * t225 + t112 * t226) * t223;
t337 = -t53 - t55;
t54 = -t224 * t103 + (-t107 * t225 + t111 * t226) * t223;
t56 = -t224 * t105 + (-t109 * t225 + t113 * t226) * t223;
t336 = t54 + t56;
t237 = sin(qJ(4));
t200 = pkin(4) * t237 + pkin(5) * t225;
t243 = -pkin(9) - pkin(8);
t231 = -qJ(6) + t243;
t262 = -t176 * rSges(7,1) - t175 * rSges(7,2);
t295 = t237 * t242;
t303 = t223 * t243;
t284 = pkin(4) * t295 + t239 * t303;
t240 = cos(qJ(4));
t220 = t240 * pkin(4) + pkin(3);
t199 = pkin(5) * t226 + t220;
t285 = t199 - t220;
t334 = rSges(7,3) * t305 - t262 - t242 * t200 + (-t223 * t231 + t285 * t224) * t239 + t284;
t281 = t231 - t243;
t333 = (t281 - rSges(7,3)) * t224 + (rSges(7,1) * t226 - rSges(7,2) * t225 + t285) * t223;
t332 = t223 / 0.2e1;
t331 = t224 / 0.2e1;
t330 = t238 / 0.2e1;
t329 = t241 / 0.2e1;
t302 = t224 * t242;
t328 = t178 * rSges(7,1) + t177 * rSges(7,2) + rSges(7,3) * t304 + t199 * t302 + t239 * t200;
t327 = -t342 * t239 + t242 * t344;
t326 = t239 * t344 + t342 * t242;
t233 = t239 ^ 2;
t234 = t242 ^ 2;
t325 = -t224 / 0.2e1;
t324 = t239 / 0.2e1;
t323 = -t242 / 0.2e1;
t322 = pkin(2) * t238;
t321 = pkin(3) * t224;
t320 = pkin(8) * t223;
t319 = -pkin(3) + t220;
t318 = t347 + (t239 * t337 - t242 * t336) * t223;
t316 = t334 * t304;
t315 = rSges(3,1) * t241;
t314 = rSges(3,2) * t238;
t313 = t242 * rSges(3,3);
t296 = t237 * t239;
t286 = -pkin(4) * t296 - t220 * t302;
t312 = -t281 * t304 + t286 + t328;
t311 = Icges(3,4) * t238;
t310 = Icges(3,4) * t241;
t309 = Icges(4,4) * t223;
t308 = Icges(4,4) * t224;
t158 = -Icges(5,6) * t224 + (Icges(5,4) * t240 - Icges(5,2) * t237) * t223;
t297 = t237 * t158;
t294 = t239 * t240;
t293 = t240 * t242;
t236 = -qJ(3) - pkin(7);
t292 = t242 * t236;
t117 = t178 * rSges(6,1) + t177 * rSges(6,2) + rSges(6,3) * t304;
t249 = -t242 * t303 - t286;
t283 = pkin(3) * t302 + pkin(8) * t304;
t125 = t249 - t283;
t291 = -t117 - t125;
t124 = (t319 * t224 - t320) * t239 - t284;
t144 = (pkin(8) + t243) * t224 + t319 * t223;
t290 = t224 * t124 + t144 * t305;
t263 = -t176 * rSges(6,1) - t175 * rSges(6,2);
t115 = rSges(6,3) * t305 - t263;
t152 = -t224 * rSges(6,3) + (rSges(6,1) * t226 - rSges(6,2) * t225) * t223;
t85 = t224 * t115 + t152 * t305;
t221 = pkin(2) * t241 + pkin(1);
t215 = t242 * t221;
t230 = t242 * pkin(7);
t288 = t239 * (t292 + t230 + (-pkin(1) + t221) * t239) + t242 * (-t242 * pkin(1) + t215 + (-pkin(7) - t236) * t239);
t282 = t239 * rSges(3,3) + t242 * t315;
t280 = t233 + t234;
t188 = -t224 * t296 - t293;
t189 = t224 * t294 - t295;
t126 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t305;
t128 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t305;
t130 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t305;
t61 = -t224 * t126 + (-t128 * t237 + t130 * t240) * t223;
t157 = -Icges(5,3) * t224 + (Icges(5,5) * t240 - Icges(5,6) * t237) * t223;
t159 = -Icges(5,5) * t224 + (Icges(5,1) * t240 - Icges(5,4) * t237) * t223;
t74 = t157 * t305 + t158 * t188 + t159 * t189;
t279 = t74 / 0.2e1 + t61 / 0.2e1;
t190 = -t224 * t295 + t294;
t191 = t224 * t293 + t296;
t127 = Icges(5,5) * t191 + Icges(5,6) * t190 + Icges(5,3) * t304;
t129 = Icges(5,4) * t191 + Icges(5,2) * t190 + Icges(5,6) * t304;
t131 = Icges(5,1) * t191 + Icges(5,4) * t190 + Icges(5,5) * t304;
t62 = -t224 * t127 + (-t129 * t237 + t131 * t240) * t223;
t75 = t157 * t304 + t158 * t190 + t159 * t191;
t278 = t75 / 0.2e1 + t62 / 0.2e1;
t277 = -t125 - t312;
t276 = t340 * t304 + t305 * t341;
t133 = t191 * rSges(5,1) + t190 * rSges(5,2) + rSges(5,3) * t304;
t275 = t305 / 0.2e1;
t274 = t304 / 0.2e1;
t273 = Icges(3,5) * t330 + Icges(4,5) * t332 + Icges(3,6) * t329 + Icges(4,6) * t331;
t272 = -rSges(4,1) * t223 - rSges(4,2) * t224 - t322;
t271 = -t223 * pkin(3) + t224 * pkin(8) - t322;
t270 = -t239 * t236 + t215;
t36 = t224 * t334 + t305 * t333;
t269 = t233 * (t320 + t321) + t242 * t283 + t288;
t268 = -t144 + t271;
t160 = -t224 * rSges(5,3) + (rSges(5,1) * t240 - rSges(5,2) * t237) * t223;
t267 = -t160 + t271;
t266 = -t314 + t315;
t265 = rSges(4,1) * t224 - rSges(4,2) * t223;
t264 = -t189 * rSges(5,1) - t188 * rSges(5,2);
t261 = Icges(3,1) * t241 - t311;
t260 = Icges(4,1) * t224 - t309;
t259 = -Icges(3,2) * t238 + t310;
t258 = -Icges(4,2) * t223 + t308;
t251 = -t152 + t268;
t250 = rSges(4,1) * t302 - rSges(4,2) * t304 + t239 * rSges(4,3);
t248 = t239 * t124 + t242 * t125 + t269;
t247 = t268 - t333;
t246 = t318 * t224 + t276;
t245 = (-t337 - t349) * t275 + (t336 - t348) * t274;
t244 = (t239 * t336 + t242 * t337) * t325 + t340 * t324 + t341 * t323 + t339 * t275 + t338 * t274;
t208 = rSges(2,1) * t242 - rSges(2,2) * t239;
t207 = -rSges(2,1) * t239 - rSges(2,2) * t242;
t206 = rSges(3,1) * t238 + rSges(3,2) * t241;
t162 = t272 * t242;
t161 = t272 * t239;
t156 = t239 * pkin(7) + (pkin(1) - t314) * t242 + t282;
t155 = t313 + t230 + (-pkin(1) - t266) * t239;
t143 = t223 * t240 * t159;
t140 = t250 + t270;
t139 = (rSges(4,3) - t236) * t242 + (-t221 - t265) * t239;
t135 = t242 * (-t242 * t314 + t282) + (t266 * t239 - t313) * t239;
t132 = rSges(5,3) * t305 - t264;
t121 = t267 * t242;
t120 = t267 * t239;
t101 = t124 * t304;
t98 = t115 * t304;
t94 = t270 + t133 + t283;
t93 = -t292 + (-t321 - t221 + (-rSges(5,3) - pkin(8)) * t223) * t239 + t264;
t90 = -t133 * t224 - t160 * t304;
t89 = t132 * t224 + t160 * t305;
t88 = t251 * t242;
t87 = t251 * t239;
t86 = -t224 * t117 - t152 * t304;
t84 = -t224 * t157 - t223 * t297 + t143;
t83 = t249 + t270 + t117;
t82 = -t292 + (-rSges(6,3) * t223 - t220 * t224 - t221) * t239 + t263 + t284;
t81 = t242 * t250 + (-t242 * rSges(4,3) + t265 * t239) * t239 + t288;
t80 = (t132 * t242 - t133 * t239) * t223;
t77 = -t231 * t304 + t270 + t328;
t76 = (t200 - t236) * t242 + (-t199 * t224 - t221 + (-rSges(7,3) + t231) * t223) * t239 + t262;
t73 = -t117 * t305 + t98;
t68 = t247 * t242;
t67 = t247 * t239;
t60 = t127 * t304 + t129 * t190 + t131 * t191;
t59 = t126 * t304 + t128 * t190 + t130 * t191;
t58 = t127 * t305 + t129 * t188 + t131 * t189;
t57 = t126 * t305 + t128 * t188 + t130 * t189;
t52 = t291 * t224 + (-t144 - t152) * t304;
t51 = t85 + t290;
t50 = t132 * t239 + t133 * t242 + t269;
t37 = -t312 * t224 - t304 * t333;
t35 = t291 * t305 + t101 + t98;
t34 = -t312 * t305 + t316;
t33 = t115 * t239 + t117 * t242 + t248;
t32 = t277 * t224 + (-t144 - t333) * t304;
t31 = t36 + t290;
t30 = t239 * t60 - t242 * t59;
t29 = t239 * t58 - t242 * t57;
t26 = t277 * t305 + t101 + t316;
t17 = t239 * t334 + t312 * t242 + t248;
t16 = -t75 * t224 + (t239 * t59 + t242 * t60) * t223;
t15 = -t74 * t224 + (t239 * t57 + t242 * t58) * t223;
t1 = [t241 * (Icges(3,2) * t241 + t311) + t238 * (Icges(3,1) * t238 + t310) + Icges(2,3) + t143 + (Icges(4,2) * t224 - t157 + t309 + t346) * t224 + (Icges(4,1) * t223 - t297 + t308 + t345) * t223 + m(7) * (t76 ^ 2 + t77 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t93 ^ 2 + t94 ^ 2) + m(4) * (t139 ^ 2 + t140 ^ 2) + m(3) * (t155 ^ 2 + t156 ^ 2) + m(2) * (t207 ^ 2 + t208 ^ 2) + t343; (-t55 / 0.2e1 - t53 / 0.2e1 - t69 / 0.2e1 - t70 / 0.2e1 + (-Icges(4,6) * t242 + t258 * t239) * t325 - t223 * (-Icges(4,5) * t242 + t260 * t239) / 0.2e1 - t241 * (-Icges(3,6) * t242 + t259 * t239) / 0.2e1 - t238 * (-Icges(3,5) * t242 + t261 * t239) / 0.2e1 + t273 * t242 - t279) * t242 + (t56 / 0.2e1 + t54 / 0.2e1 + t71 / 0.2e1 + t72 / 0.2e1 + (Icges(4,6) * t239 + t258 * t242) * t331 + (Icges(4,5) * t239 + t260 * t242) * t332 + (Icges(3,6) * t239 + t259 * t242) * t329 + (Icges(3,5) * t239 + t261 * t242) * t330 + t273 * t239 + t278) * t239 + m(4) * (t139 * t162 + t140 * t161) + m(5) * (t120 * t94 + t121 * t93) + m(7) * (t67 * t77 + t68 * t76) + m(6) * (t82 * t88 + t83 * t87) + m(3) * (-t155 * t242 - t156 * t239) * t206; m(7) * (t17 ^ 2 + t67 ^ 2 + t68 ^ 2) + m(6) * (t33 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(5) * (t120 ^ 2 + t121 ^ 2 + t50 ^ 2) + m(4) * (t161 ^ 2 + t162 ^ 2 + t81 ^ 2) + m(3) * (t280 * t206 ^ 2 + t135 ^ 2) + (t234 * t327 - t29 - t339) * t242 + (t30 + t326 * t233 + (t327 * t239 + t242 * t326) * t242 + t338) * t239; m(7) * (t239 * t76 - t242 * t77) + m(6) * (t239 * t82 - t242 * t83) + m(5) * (t239 * t93 - t242 * t94) + m(4) * (t139 * t239 - t140 * t242); m(7) * (t239 * t68 - t242 * t67) + m(6) * (t239 * t88 - t242 * t87) + m(5) * (-t120 * t242 + t121 * t239) + m(4) * (-t161 * t242 + t162 * t239); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * t280; (-t84 - t335) * t224 + m(7) * (t31 * t76 + t32 * t77) + m(6) * (t51 * t82 + t52 * t83) + m(5) * (t89 * t93 + t90 * t94) + (t279 * t239 + t278 * t242) * t223 + t245; t244 + m(7) * (t17 * t26 + t31 * t68 + t32 * t67) + m(6) * (t33 * t35 + t51 * t88 + t52 * t87) + m(5) * (t120 * t90 + t121 * t89 + t50 * t80) + (t29 * t324 + t242 * t30 / 0.2e1) * t223 + t16 * t324 + t15 * t323 + (t62 * t239 - t61 * t242) * t325; m(5) * (t239 * t89 - t242 * t90) + m(6) * (t239 * t51 - t242 * t52) + m(7) * (t239 * t31 - t242 * t32); (t84 * t224 + t318) * t224 + (t242 * t16 + t239 * t15 - t224 * (t239 * t61 + t242 * t62)) * t223 + m(7) * (t26 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t35 ^ 2 + t51 ^ 2 + t52 ^ 2) + m(5) * (t80 ^ 2 + t89 ^ 2 + t90 ^ 2) + t276; -t347 + m(7) * (t36 * t76 + t37 * t77) + m(6) * (t82 * t85 + t83 * t86) + t245; m(7) * (t17 * t34 + t36 * t68 + t37 * t67) + m(6) * (t33 * t73 + t85 * t88 + t86 * t87) + t244; m(6) * (t239 * t85 - t242 * t86) + m(7) * (t239 * t36 - t242 * t37); m(7) * (t26 * t34 + t31 * t36 + t32 * t37) + m(6) * (t35 * t73 + t51 * t85 + t52 * t86) + t246; m(7) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t73 ^ 2 + t85 ^ 2 + t86 ^ 2) + t246; m(7) * (t239 * t77 + t242 * t76) * t223; m(7) * (-t224 * t17 + (t239 * t67 + t242 * t68) * t223); 0; m(7) * (-t224 * t26 + (t239 * t32 + t242 * t31) * t223); m(7) * (-t224 * t34 + (t239 * t37 + t242 * t36) * t223); m(7) * (t223 ^ 2 * t280 + t224 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
