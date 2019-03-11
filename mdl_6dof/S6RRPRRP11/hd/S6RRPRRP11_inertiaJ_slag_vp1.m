% Calculate joint inertia matrix for
% S6RRPRRP11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
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
% Datum: 2019-03-09 12:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP11_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP11_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP11_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP11_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP11_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP11_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:45:28
% EndTime: 2019-03-09 12:45:38
% DurationCPUTime: 5.14s
% Computational Cost: add. (8438->489), mult. (12327->673), div. (0->0), fcn. (13165->8), ass. (0->233)
t236 = qJ(4) + qJ(5);
t221 = sin(t236);
t222 = cos(t236);
t239 = sin(qJ(1));
t238 = sin(qJ(2));
t242 = cos(qJ(1));
t301 = t238 * t242;
t176 = -t221 * t239 + t222 * t301;
t177 = t221 * t301 + t222 * t239;
t241 = cos(qJ(2));
t297 = t241 * t242;
t101 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t297;
t103 = Icges(6,5) * t177 + Icges(6,6) * t176 + Icges(6,3) * t297;
t351 = t101 + t103;
t302 = t238 * t239;
t178 = t221 * t242 + t222 * t302;
t179 = t221 * t302 - t222 * t242;
t299 = t239 * t241;
t102 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t299;
t104 = Icges(6,5) * t179 + Icges(6,6) * t178 + Icges(6,3) * t299;
t350 = t102 + t104;
t105 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t297;
t107 = Icges(6,4) * t177 + Icges(6,2) * t176 + Icges(6,6) * t297;
t349 = t105 + t107;
t106 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t299;
t108 = Icges(6,4) * t179 + Icges(6,2) * t178 + Icges(6,6) * t299;
t348 = t106 + t108;
t109 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t297;
t111 = Icges(6,1) * t177 + Icges(6,4) * t176 + Icges(6,5) * t297;
t347 = t109 + t111;
t110 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t299;
t112 = Icges(6,1) * t179 + Icges(6,4) * t178 + Icges(6,5) * t299;
t346 = t110 + t112;
t345 = t349 * t176 + t347 * t177 + t351 * t297;
t344 = t348 * t176 + t346 * t177 + t350 * t297;
t343 = t349 * t178 + t347 * t179 + t351 * t299;
t342 = t348 * t178 + t346 * t179 + t350 * t299;
t146 = Icges(7,3) * t238 + (-Icges(7,5) * t221 - Icges(7,6) * t222) * t241;
t148 = Icges(7,6) * t238 + (-Icges(7,4) * t221 - Icges(7,2) * t222) * t241;
t150 = Icges(7,5) * t238 + (-Icges(7,1) * t221 - Icges(7,4) * t222) * t241;
t67 = t146 * t297 + t148 * t176 + t150 * t177;
t147 = Icges(6,3) * t238 + (-Icges(6,5) * t221 - Icges(6,6) * t222) * t241;
t149 = Icges(6,6) * t238 + (-Icges(6,4) * t221 - Icges(6,2) * t222) * t241;
t151 = Icges(6,5) * t238 + (-Icges(6,1) * t221 - Icges(6,4) * t222) * t241;
t68 = t147 * t297 + t149 * t176 + t151 * t177;
t341 = t68 + t67;
t69 = t146 * t299 + t148 * t178 + t150 * t179;
t70 = t147 * t299 + t149 * t178 + t151 * t179;
t340 = t69 + t70;
t339 = (t146 + t147) * t238;
t338 = (-t148 - t149) * t222 + (-t150 - t151) * t221;
t337 = Icges(4,1) + Icges(3,3);
t240 = cos(qJ(4));
t220 = t240 * pkin(4) + pkin(3);
t194 = pkin(5) * t222 + t220;
t237 = sin(qJ(4));
t315 = pkin(4) * t237;
t195 = pkin(5) * t221 + t315;
t336 = t177 * rSges(7,1) + t176 * rSges(7,2) + rSges(7,3) * t297 + t239 * t194 + t195 * t301;
t335 = (-Icges(4,4) + Icges(3,5)) * t241 + (Icges(4,5) - Icges(3,6)) * t238;
t334 = -rSges(7,1) * t179 - rSges(7,2) * t178 + t194 * t242;
t333 = -t239 / 0.2e1;
t317 = t239 / 0.2e1;
t316 = -t242 / 0.2e1;
t332 = t242 / 0.2e1;
t331 = (t344 * t239 + t345 * t242) * t241 + t341 * t238;
t330 = (t342 * t239 + t343 * t242) * t241 + t340 * t238;
t329 = t345 * t239 - t344 * t242;
t328 = t343 * t239 - t342 * t242;
t50 = t101 * t238 + (-t105 * t222 - t109 * t221) * t241;
t52 = t103 * t238 + (-t107 * t222 - t111 * t221) * t241;
t327 = t50 + t52;
t51 = t102 * t238 + (-t106 * t222 - t110 * t221) * t241;
t53 = t104 * t238 + (-t108 * t222 - t112 * t221) * t241;
t326 = t51 + t53;
t325 = (t338 * t241 + t339) * t238;
t243 = -pkin(9) - pkin(8);
t231 = -qJ(6) + t243;
t285 = -t231 + t243;
t281 = t237 * t301;
t290 = -pkin(4) * t281 - t239 * t220;
t310 = t285 * t297 + t290 + t336;
t277 = -t195 + t315;
t296 = t241 * t243;
t289 = t242 * t220 + t239 * t296;
t309 = rSges(7,3) * t299 + (-t231 * t241 - t238 * t277) * t239 + t289 - t334;
t324 = (-rSges(7,1) * t221 - rSges(7,2) * t222 + t277) * t241 + (t285 + rSges(7,3)) * t238;
t318 = t238 / 0.2e1;
t323 = t337 * t239 + t335 * t242;
t322 = -t335 * t239 + t337 * t242;
t233 = t239 ^ 2;
t235 = t242 ^ 2;
t321 = m(4) / 0.2e1;
t320 = m(5) / 0.2e1;
t319 = m(6) / 0.2e1;
t314 = t309 * t297;
t313 = t310 * t238;
t312 = t242 * rSges(4,1);
t311 = t242 * rSges(3,3);
t308 = Icges(3,4) * t238;
t307 = Icges(3,4) * t241;
t306 = Icges(4,6) * t238;
t305 = Icges(4,6) * t241;
t304 = qJ(3) * t238;
t300 = t239 * t240;
t298 = t240 * t242;
t295 = t324 * t299;
t288 = pkin(2) * t297 + qJ(3) * t301;
t293 = t233 * (pkin(2) * t241 + t304) + t242 * t288;
t203 = pkin(2) * t238 - qJ(3) * t241;
t291 = rSges(4,2) * t238 + rSges(4,3) * t241 - t203;
t287 = t242 * pkin(1) + t239 * pkin(7);
t286 = t239 * pkin(3) + pkin(8) * t297;
t284 = t233 + t235;
t187 = -t237 * t239 + t238 * t298;
t188 = t281 + t300;
t121 = Icges(5,5) * t188 + Icges(5,6) * t187 + Icges(5,3) * t297;
t123 = Icges(5,4) * t188 + Icges(5,2) * t187 + Icges(5,6) * t297;
t125 = Icges(5,1) * t188 + Icges(5,4) * t187 + Icges(5,5) * t297;
t61 = t121 * t238 + (-t123 * t240 - t125 * t237) * t241;
t160 = Icges(5,3) * t238 + (-Icges(5,5) * t237 - Icges(5,6) * t240) * t241;
t163 = Icges(5,6) * t238 + (-Icges(5,4) * t237 - Icges(5,2) * t240) * t241;
t166 = Icges(5,5) * t238 + (-Icges(5,1) * t237 - Icges(5,4) * t240) * t241;
t78 = t160 * t297 + t163 * t187 + t166 * t188;
t283 = t61 / 0.2e1 + t78 / 0.2e1;
t189 = t237 * t242 + t238 * t300;
t190 = t237 * t302 - t298;
t122 = Icges(5,5) * t190 + Icges(5,6) * t189 + Icges(5,3) * t299;
t124 = Icges(5,4) * t190 + Icges(5,2) * t189 + Icges(5,6) * t299;
t126 = Icges(5,1) * t190 + Icges(5,4) * t189 + Icges(5,5) * t299;
t62 = t122 * t238 + (-t124 * t240 - t126 * t237) * t241;
t79 = t160 * t299 + t163 * t189 + t166 * t190;
t282 = t79 / 0.2e1 + t62 / 0.2e1;
t114 = t177 * rSges(6,1) + t176 * rSges(6,2) + rSges(6,3) * t297;
t130 = t188 * rSges(5,1) + t187 * rSges(5,2) + rSges(5,3) * t297;
t280 = t299 / 0.2e1;
t279 = t297 / 0.2e1;
t278 = -Icges(4,4) * t238 / 0.2e1 + Icges(3,5) * t318 + (-Icges(4,5) / 0.2e1 + Icges(3,6) / 0.2e1) * t241;
t276 = -t238 * pkin(8) - t203;
t229 = t242 * pkin(3);
t275 = t239 * (pkin(8) * t299 - t229) + t242 * t286 + t293;
t274 = t287 + t288;
t175 = rSges(5,3) * t238 + (-rSges(5,1) * t237 - rSges(5,2) * t240) * t241;
t273 = -t175 + t276;
t182 = -t241 * t315 + (-pkin(8) - t243) * t238;
t272 = -t182 + t276;
t271 = t330 * t299 + t331 * t297 + ((t326 * t239 + t327 * t242) * t241 + t325) * t238;
t270 = rSges(3,1) * t241 - rSges(3,2) * t238;
t269 = -rSges(5,1) * t190 - rSges(5,2) * t189;
t268 = -rSges(6,1) * t179 - rSges(6,2) * t178;
t136 = t229 + (-pkin(8) * t241 + t238 * t315) * t239 - t289;
t152 = t182 * t299;
t32 = t152 + (-t136 - t309) * t238 + t295;
t249 = -t242 * t296 - t290;
t135 = t249 - t286;
t127 = t238 * t135;
t33 = t127 + (-t182 - t324) * t297 + t313;
t266 = t239 * t33 + t242 * t32;
t44 = -t238 * t309 + t295;
t45 = -t297 * t324 + t313;
t265 = t239 * t45 + t242 * t44;
t246 = t272 - t324;
t74 = t246 * t239;
t75 = t246 * t242;
t264 = t239 * t74 + t242 * t75;
t263 = Icges(3,1) * t241 - t308;
t262 = -Icges(3,2) * t238 + t307;
t259 = -Icges(4,2) * t241 + t306;
t258 = Icges(4,3) * t238 - t305;
t257 = -t240 * t163 - t237 * t166;
t154 = rSges(6,3) * t238 + (-rSges(6,1) * t221 - rSges(6,2) * t222) * t241;
t252 = -t154 + t272;
t251 = rSges(3,1) * t297 - rSges(3,2) * t301 + t239 * rSges(3,3);
t250 = t239 * rSges(4,1) - rSges(4,2) * t297 + rSges(4,3) * t301;
t248 = t242 * t135 + t239 * t136 + t275;
t228 = t242 * pkin(7);
t71 = t228 + (-pkin(1) + (-qJ(3) - t195) * t238 + (-rSges(7,3) - pkin(2) + t231) * t241) * t239 + t334;
t72 = -t231 * t297 + t274 + t336;
t247 = m(7) * (t239 * t72 + t242 * t71);
t245 = (t327 * t239 - t326 * t242) * t318 + t331 * t317 + t330 * t316 + t328 * t280 + t329 * t279;
t244 = (t326 + t340) * t280 + (t327 + t341) * t279 + t325;
t234 = t241 ^ 2;
t232 = t238 ^ 2;
t207 = rSges(2,1) * t242 - rSges(2,2) * t239;
t206 = -rSges(2,1) * t239 - rSges(2,2) * t242;
t205 = rSges(3,1) * t238 + rSges(3,2) * t241;
t155 = t238 * t160;
t145 = t291 * t242;
t144 = t291 * t239;
t141 = t251 + t287;
t140 = t311 + t228 + (-pkin(1) - t270) * t239;
t139 = t154 * t299;
t133 = t250 + t274;
t132 = t312 + t228 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t241 + (-rSges(4,3) - qJ(3)) * t238) * t239;
t131 = rSges(5,3) * t299 - t269;
t120 = t273 * t242;
t119 = t273 * t239;
t118 = t136 * t297;
t117 = t242 * t251 + (t239 * t270 - t311) * t239;
t116 = rSges(6,3) * t299 - t268;
t100 = t238 * t114;
t98 = t116 * t297;
t96 = t252 * t242;
t95 = t252 * t239;
t92 = t130 * t238 - t175 * t297;
t91 = -t131 * t238 + t175 * t299;
t89 = t274 + t130 + t286;
t88 = t228 + t229 + (-t304 - pkin(1) + (-rSges(5,3) - pkin(2) - pkin(8)) * t241) * t239 + t269;
t86 = (t241 * t257 + t155) * t238;
t85 = t242 * t250 + (-t312 + (-rSges(4,2) * t241 + rSges(4,3) * t238) * t239) * t239 + t293;
t84 = -t154 * t297 + t100;
t83 = -t116 * t238 + t139;
t80 = (-t130 * t239 + t131 * t242) * t241;
t77 = t249 + t274 + t114;
t76 = t228 + (-pkin(1) + (-rSges(6,3) - pkin(2)) * t241 + (-qJ(3) - t315) * t238) * t239 + t268 + t289;
t73 = -t114 * t299 + t98;
t60 = t100 + t127 + (-t154 - t182) * t297;
t59 = t139 + t152 + (-t116 - t136) * t238;
t58 = t130 * t242 + t131 * t239 + t275;
t57 = t122 * t299 + t124 * t189 + t126 * t190;
t56 = t121 * t299 + t123 * t189 + t125 * t190;
t55 = t122 * t297 + t124 * t187 + t126 * t188;
t54 = t121 * t297 + t123 * t187 + t125 * t188;
t35 = t118 + t98 + (-t114 - t135) * t299;
t34 = t114 * t242 + t116 * t239 + t248;
t31 = -t299 * t310 + t314;
t30 = t239 * t56 - t242 * t57;
t29 = t239 * t54 - t242 * t55;
t26 = t118 + (-t135 - t310) * t299 + t314;
t21 = t239 * t309 + t242 * t310 + t248;
t16 = t238 * t79 + (t239 * t57 + t242 * t56) * t241;
t15 = t238 * t78 + (t239 * t55 + t242 * t54) * t241;
t1 = [Icges(2,3) + t155 + m(7) * (t71 ^ 2 + t72 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t88 ^ 2 + t89 ^ 2) + m(4) * (t132 ^ 2 + t133 ^ 2) + m(3) * (t140 ^ 2 + t141 ^ 2) + m(2) * (t206 ^ 2 + t207 ^ 2) + (t306 + t308 + t257 + (Icges(4,3) + Icges(3,2)) * t241 + t338) * t241 + (t305 + t307 + (Icges(3,1) + Icges(4,2)) * t238) * t238 + t339; (-t69 / 0.2e1 - t70 / 0.2e1 - t53 / 0.2e1 - t51 / 0.2e1 + t278 * t242 + (Icges(4,5) * t316 + Icges(3,6) * t332 + t258 * t317 + t262 * t333) * t241 + (Icges(4,4) * t316 + Icges(3,5) * t332 + t259 * t317 + t263 * t333) * t238 - t282) * t242 + (t67 / 0.2e1 + t68 / 0.2e1 + t52 / 0.2e1 + t50 / 0.2e1 + t278 * t239 + (Icges(4,5) * t333 + Icges(3,6) * t317 + t258 * t316 + t262 * t332) * t241 + (Icges(4,4) * t333 + Icges(3,5) * t317 + t259 * t316 + t263 * t332) * t238 + t283) * t239 + m(7) * (t71 * t75 + t72 * t74) + m(6) * (t76 * t96 + t77 * t95) + m(5) * (t119 * t89 + t120 * t88) + m(4) * (t132 * t145 + t133 * t144) + m(3) * (-t140 * t242 - t141 * t239) * t205; m(7) * (t21 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(6) * (t34 ^ 2 + t95 ^ 2 + t96 ^ 2) + m(5) * (t119 ^ 2 + t120 ^ 2 + t58 ^ 2) + m(4) * (t144 ^ 2 + t145 ^ 2 + t85 ^ 2) + m(3) * (t205 ^ 2 * t284 + t117 ^ 2) + (t235 * t322 - t30 - t328) * t242 + (t29 + t323 * t233 + (t239 * t322 + t242 * t323) * t242 + t329) * t239; 0.2e1 * (t247 / 0.2e1 + (t239 * t77 + t242 * t76) * t319 + (t239 * t89 + t242 * t88) * t320 + (t132 * t242 + t133 * t239) * t321) * t238; m(7) * (-t241 * t21 + t238 * t264) + m(6) * (-t241 * t34 + (t239 * t95 + t242 * t96) * t238) + m(5) * (-t241 * t58 + (t119 * t239 + t120 * t242) * t238) + m(4) * (-t241 * t85 + (t144 * t239 + t145 * t242) * t238); 0.2e1 * (t321 + t320 + t319 + m(7) / 0.2e1) * (t232 * t284 + t234); t244 + (t239 * t282 + t242 * t283) * t241 + m(7) * (t32 * t71 + t33 * t72) + m(6) * (t59 * t76 + t60 * t77) + m(5) * (t88 * t91 + t89 * t92) + t86; t16 * t316 + (t61 * t239 - t62 * t242) * t318 + t245 + m(7) * (t26 * t21 + t32 * t75 + t33 * t74) + m(6) * (t35 * t34 + t59 * t96 + t60 * t95) + m(5) * (t119 * t92 + t120 * t91 + t58 * t80) + (t29 * t332 + t30 * t317) * t241 + t15 * t317; m(5) * (-t241 * t80 + (t239 * t92 + t242 * t91) * t238) + m(6) * (-t241 * t35 + (t239 * t60 + t242 * t59) * t238) + m(7) * (t238 * t266 - t241 * t26); t238 * t86 + (t242 * t15 + t239 * t16 + t238 * (t239 * t62 + t242 * t61)) * t241 + m(7) * (t26 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t35 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t80 ^ 2 + t91 ^ 2 + t92 ^ 2) + t271; m(7) * (t44 * t71 + t45 * t72) + m(6) * (t76 * t83 + t77 * t84) + t244; m(7) * (t31 * t21 + t44 * t75 + t45 * t74) + m(6) * (t34 * t73 + t83 * t96 + t84 * t95) + t245; m(6) * (-t241 * t73 + (t239 * t84 + t242 * t83) * t238) + m(7) * (t238 * t265 - t241 * t31); m(7) * (t26 * t31 + t32 * t44 + t33 * t45) + m(6) * (t35 * t73 + t59 * t83 + t60 * t84) + t271; m(7) * (t31 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t73 ^ 2 + t83 ^ 2 + t84 ^ 2) + t271; t241 * t247; m(7) * (t21 * t238 + t241 * t264); m(7) * (-0.1e1 + t284) * t241 * t238; m(7) * (t238 * t26 + t241 * t266); m(7) * (t238 * t31 + t241 * t265); m(7) * (t234 * t284 + t232);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
