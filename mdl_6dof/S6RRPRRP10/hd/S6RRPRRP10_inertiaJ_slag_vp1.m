% Calculate joint inertia matrix for
% S6RRPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP10_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP10_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP10_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP10_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP10_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP10_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:36:36
% EndTime: 2019-03-09 12:36:49
% DurationCPUTime: 6.27s
% Computational Cost: add. (24740->622), mult. (41613->861), div. (0->0), fcn. (52796->12), ass. (0->273)
t263 = sin(pkin(11));
t265 = cos(pkin(11));
t266 = cos(pkin(6));
t264 = sin(pkin(6));
t269 = sin(qJ(2));
t329 = t264 * t269;
t243 = -t263 * t329 + t265 * t266;
t330 = t263 * t266;
t244 = t265 * t329 + t330;
t271 = cos(qJ(2));
t327 = t264 * t271;
t187 = Icges(4,4) * t244 + Icges(4,2) * t243 - Icges(4,6) * t327;
t188 = Icges(4,1) * t244 + Icges(4,4) * t243 - Icges(4,5) * t327;
t228 = Icges(3,3) * t266 + (Icges(3,5) * t269 + Icges(3,6) * t271) * t264;
t229 = Icges(3,6) * t266 + (Icges(3,4) * t269 + Icges(3,2) * t271) * t264;
t230 = Icges(3,5) * t266 + (Icges(3,1) * t269 + Icges(3,4) * t271) * t264;
t339 = t243 * t187 + t244 * t188 + t266 * t228 + t229 * t327 + t230 * t329;
t270 = sin(qJ(1));
t323 = t270 * t271;
t272 = cos(qJ(1));
t324 = t269 * t272;
t246 = t266 * t324 + t323;
t308 = pkin(11) + qJ(4);
t261 = sin(t308);
t289 = cos(t308);
t326 = t264 * t272;
t215 = t246 * t289 - t261 * t326;
t322 = t271 * t272;
t325 = t269 * t270;
t245 = -t266 * t322 + t325;
t268 = sin(qJ(5));
t334 = cos(qJ(5));
t176 = t215 * t268 - t245 * t334;
t177 = t215 * t334 + t245 * t268;
t284 = t264 * t289;
t214 = t246 * t261 + t272 * t284;
t337 = rSges(7,3) + qJ(6);
t338 = rSges(7,1) + pkin(5);
t321 = t214 * rSges(7,2) + t337 * t176 + t177 * t338;
t232 = t266 * t261 + t269 * t284;
t212 = t232 * t268 + t327 * t334;
t213 = t232 * t334 - t268 * t327;
t231 = t261 * t329 - t266 * t289;
t126 = Icges(7,5) * t213 + Icges(7,6) * t231 + Icges(7,3) * t212;
t128 = Icges(7,4) * t213 + Icges(7,2) * t231 + Icges(7,6) * t212;
t130 = Icges(7,1) * t213 + Icges(7,4) * t231 + Icges(7,5) * t212;
t63 = t212 * t126 + t231 * t128 + t213 * t130;
t127 = Icges(6,5) * t213 - Icges(6,6) * t212 + Icges(6,3) * t231;
t129 = Icges(6,4) * t213 - Icges(6,2) * t212 + Icges(6,6) * t231;
t131 = Icges(6,1) * t213 - Icges(6,4) * t212 + Icges(6,5) * t231;
t64 = t231 * t127 - t212 * t129 + t213 * t131;
t336 = -t63 - t64;
t186 = Icges(4,5) * t244 + Icges(4,6) * t243 - Icges(4,3) * t327;
t335 = (-t186 * t327 + t339) * t266;
t260 = pkin(3) * t265 + pkin(2);
t333 = -pkin(2) + t260;
t267 = -pkin(9) - qJ(3);
t331 = t245 * t267;
t328 = t264 * t270;
t248 = -t266 * t325 + t322;
t217 = t248 * t289 + t261 * t328;
t247 = t266 * t323 + t324;
t178 = t217 * t268 - t247 * t334;
t179 = t217 * t334 + t247 * t268;
t216 = t248 * t261 - t270 * t284;
t320 = t216 * rSges(7,2) + t337 * t178 + t179 * t338;
t108 = t179 * rSges(6,1) - t178 * rSges(6,2) + t216 * rSges(6,3);
t167 = t217 * pkin(4) + pkin(10) * t216;
t319 = -t108 - t167;
t318 = rSges(7,2) * t231 + t337 * t212 + t213 * t338;
t205 = t248 * pkin(2) + qJ(3) * t247;
t298 = t263 * t328;
t291 = pkin(3) * t298 - t247 * t267 + t248 * t260;
t146 = -t205 + t291;
t203 = t266 * t205;
t317 = t266 * t146 + t203;
t237 = t245 * qJ(3);
t297 = t263 * t326;
t252 = pkin(3) * t297;
t145 = t246 * t333 - t237 - t252 - t331;
t204 = pkin(2) * t246 + t237;
t316 = -t145 - t204;
t220 = -t246 * t263 - t265 * t326;
t221 = t246 * t265 - t297;
t155 = rSges(4,1) * t221 + rSges(4,2) * t220 + rSges(4,3) * t245;
t315 = -t155 - t204;
t166 = t215 * pkin(4) + t214 * pkin(10);
t190 = pkin(4) * t232 + pkin(10) * t231;
t314 = t166 * t327 + t245 * t190;
t183 = Icges(5,4) * t232 - Icges(5,2) * t231 - Icges(5,6) * t327;
t184 = Icges(5,1) * t232 - Icges(5,4) * t231 - Icges(5,5) * t327;
t313 = -t231 * t183 + t232 * t184;
t249 = (pkin(2) * t269 - qJ(3) * t271) * t264;
t311 = -pkin(3) * t330 - ((qJ(3) + t267) * t271 + t333 * t269) * t264 - t249;
t310 = t204 * t328 + t205 * t326;
t309 = t272 * pkin(1) + pkin(8) * t328;
t101 = Icges(7,1) * t177 + Icges(7,4) * t214 + Icges(7,5) * t176;
t93 = Icges(7,5) * t177 + Icges(7,6) * t214 + Icges(7,3) * t176;
t97 = Icges(7,4) * t177 + Icges(7,2) * t214 + Icges(7,6) * t176;
t33 = t101 * t177 + t176 * t93 + t214 * t97;
t102 = Icges(7,1) * t179 + Icges(7,4) * t216 + Icges(7,5) * t178;
t94 = Icges(7,5) * t179 + Icges(7,6) * t216 + Icges(7,3) * t178;
t98 = Icges(7,4) * t179 + Icges(7,2) * t216 + Icges(7,6) * t178;
t34 = t102 * t177 + t176 * t94 + t214 * t98;
t51 = t126 * t176 + t128 * t214 + t130 * t177;
t1 = t214 * t33 + t216 * t34 + t231 * t51;
t103 = Icges(6,1) * t177 - Icges(6,4) * t176 + Icges(6,5) * t214;
t95 = Icges(6,5) * t177 - Icges(6,6) * t176 + Icges(6,3) * t214;
t99 = Icges(6,4) * t177 - Icges(6,2) * t176 + Icges(6,6) * t214;
t35 = t103 * t177 - t176 * t99 + t214 * t95;
t100 = Icges(6,4) * t179 - Icges(6,2) * t178 + Icges(6,6) * t216;
t104 = Icges(6,1) * t179 - Icges(6,4) * t178 + Icges(6,5) * t216;
t96 = Icges(6,5) * t179 - Icges(6,6) * t178 + Icges(6,3) * t216;
t36 = -t100 * t176 + t104 * t177 + t214 * t96;
t52 = t127 * t214 - t129 * t176 + t131 * t177;
t2 = t214 * t35 + t216 * t36 + t231 * t52;
t307 = t2 / 0.2e1 + t1 / 0.2e1;
t37 = t101 * t179 + t178 * t93 + t216 * t97;
t38 = t102 * t179 + t178 * t94 + t216 * t98;
t53 = t126 * t178 + t128 * t216 + t130 * t179;
t3 = t214 * t37 + t216 * t38 + t231 * t53;
t39 = t103 * t179 - t178 * t99 + t216 * t95;
t40 = -t100 * t178 + t104 * t179 + t216 * t96;
t54 = t127 * t216 - t129 * t178 + t131 * t179;
t4 = t214 * t39 + t216 * t40 + t231 * t54;
t306 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = t245 * t33 + t247 * t34 - t327 * t51;
t6 = t245 * t35 + t247 * t36 - t327 * t52;
t305 = t5 / 0.2e1 + t6 / 0.2e1;
t7 = t245 * t37 + t247 * t38 - t327 * t53;
t8 = t245 * t39 + t247 * t40 - t327 * t54;
t304 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t52 * t266 + (t270 * t36 - t272 * t35) * t264;
t9 = t51 * t266 + (t270 * t34 - t272 * t33) * t264;
t303 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t53 * t266 + (t270 * t38 - t272 * t37) * t264;
t12 = t54 * t266 + (t270 * t40 - t272 * t39) * t264;
t302 = t11 / 0.2e1 + t12 / 0.2e1;
t43 = t101 * t213 + t212 * t93 + t231 * t97;
t44 = t102 * t213 + t212 * t94 + t231 * t98;
t59 = t63 * t231;
t13 = t43 * t214 + t44 * t216 + t59;
t45 = t103 * t213 - t212 * t99 + t231 * t95;
t46 = -t100 * t212 + t104 * t213 + t231 * t96;
t60 = t64 * t231;
t14 = t45 * t214 + t46 * t216 + t60;
t301 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t43 * t245 + t44 * t247 - t327 * t63;
t16 = t45 * t245 + t46 * t247 - t327 * t64;
t300 = t15 / 0.2e1 + t16 / 0.2e1;
t61 = t63 * t266;
t17 = t61 + (t44 * t270 - t43 * t272) * t264;
t62 = t64 * t266;
t18 = t62 + (t46 * t270 - t45 * t272) * t264;
t299 = t17 / 0.2e1 + t18 / 0.2e1;
t296 = -t167 - t320;
t295 = t266 * t167 + t317;
t294 = -t166 + t316;
t293 = -t190 + t311;
t144 = t217 * rSges(5,1) - t216 * rSges(5,2) + t247 * rSges(5,3);
t222 = -t248 * t263 + t265 * t328;
t223 = t248 * t265 + t298;
t156 = t223 * rSges(4,1) + t222 * rSges(4,2) + t247 * rSges(4,3);
t199 = t248 * rSges(3,1) - t247 * rSges(3,2) + rSges(3,3) * t328;
t290 = -t270 * pkin(1) + pkin(8) * t326;
t288 = t264 * (-rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t327 - t249);
t287 = t145 * t328 + t146 * t326 + t310;
t185 = rSges(5,1) * t232 - rSges(5,2) * t231 - rSges(5,3) * t327;
t286 = t264 * (-t185 + t311);
t285 = -t215 * rSges(5,1) + t214 * rSges(5,2);
t283 = t291 + t309;
t133 = rSges(6,1) * t213 - rSges(6,2) * t212 + rSges(6,3) * t231;
t282 = t264 * (-t133 + t293);
t281 = t45 / 0.2e1 + t43 / 0.2e1 + t52 / 0.2e1 + t51 / 0.2e1;
t280 = t46 / 0.2e1 + t44 / 0.2e1 + t54 / 0.2e1 + t53 / 0.2e1;
t279 = t166 * t328 + t167 * t326 + t287;
t278 = t264 * (t293 - t318);
t277 = -t246 * t260 + t252 + t290;
t106 = rSges(6,1) * t177 - rSges(6,2) * t176 + rSges(6,3) * t214;
t198 = t246 * rSges(3,1) - t245 * rSges(3,2) - rSges(3,3) * t326;
t276 = t167 + t283;
t137 = Icges(5,5) * t215 - Icges(5,6) * t214 + Icges(5,3) * t245;
t139 = Icges(5,4) * t215 - Icges(5,2) * t214 + Icges(5,6) * t245;
t141 = Icges(5,1) * t215 - Icges(5,4) * t214 + Icges(5,5) * t245;
t72 = -t137 * t327 - t139 * t231 + t141 * t232;
t182 = Icges(5,5) * t232 - Icges(5,6) * t231 - Icges(5,3) * t327;
t81 = t182 * t245 - t183 * t214 + t184 * t215;
t275 = t81 / 0.2e1 + t72 / 0.2e1 + t281;
t138 = Icges(5,5) * t217 - Icges(5,6) * t216 + Icges(5,3) * t247;
t140 = Icges(5,4) * t217 - Icges(5,2) * t216 + Icges(5,6) * t247;
t142 = Icges(5,1) * t217 - Icges(5,4) * t216 + Icges(5,5) * t247;
t73 = -t138 * t327 - t140 * t231 + t142 * t232;
t82 = t182 * t247 - t183 * t216 + t184 * t217;
t274 = t73 / 0.2e1 + t82 / 0.2e1 + t280;
t273 = -t166 + t277 + t331;
t254 = rSges(2,1) * t272 - t270 * rSges(2,2);
t253 = -t270 * rSges(2,1) - rSges(2,2) * t272;
t233 = t266 * rSges(3,3) + (rSges(3,1) * t269 + rSges(3,2) * t271) * t264;
t197 = Icges(3,1) * t248 - Icges(3,4) * t247 + Icges(3,5) * t328;
t196 = Icges(3,1) * t246 - Icges(3,4) * t245 - Icges(3,5) * t326;
t195 = Icges(3,4) * t248 - Icges(3,2) * t247 + Icges(3,6) * t328;
t194 = Icges(3,4) * t246 - Icges(3,2) * t245 - Icges(3,6) * t326;
t193 = Icges(3,5) * t248 - Icges(3,6) * t247 + Icges(3,3) * t328;
t192 = Icges(3,5) * t246 - Icges(3,6) * t245 - Icges(3,3) * t326;
t181 = t199 + t309;
t180 = -t198 + t290;
t164 = -t266 * t198 - t233 * t326;
t163 = t199 * t266 - t233 * t328;
t154 = Icges(4,1) * t223 + Icges(4,4) * t222 + Icges(4,5) * t247;
t153 = Icges(4,1) * t221 + Icges(4,4) * t220 + Icges(4,5) * t245;
t152 = Icges(4,4) * t223 + Icges(4,2) * t222 + Icges(4,6) * t247;
t151 = Icges(4,4) * t221 + Icges(4,2) * t220 + Icges(4,6) * t245;
t150 = Icges(4,5) * t223 + Icges(4,6) * t222 + Icges(4,3) * t247;
t149 = Icges(4,5) * t221 + Icges(4,6) * t220 + Icges(4,3) * t245;
t147 = t247 * t166;
t143 = rSges(5,3) * t245 - t285;
t125 = (t198 * t270 + t199 * t272) * t264;
t124 = t228 * t328 - t229 * t247 + t230 * t248;
t123 = -t228 * t326 - t245 * t229 + t246 * t230;
t118 = t205 + t156 + t309;
t117 = t290 + t315;
t112 = t266 * t193 + (t195 * t271 + t197 * t269) * t264;
t111 = t266 * t192 + (t194 * t271 + t196 * t269) * t264;
t110 = t283 + t144;
t109 = (-rSges(5,3) + t267) * t245 + t277 + t285;
t92 = -t144 * t327 - t185 * t247;
t91 = t143 * t327 + t185 * t245;
t89 = t266 * t315 + t272 * t288;
t88 = t266 * t156 + t270 * t288 + t203;
t87 = -t182 * t327 + t313;
t86 = t143 * t247 - t144 * t245;
t85 = t87 * t266;
t84 = t186 * t247 + t187 * t222 + t188 * t223;
t83 = t186 * t245 + t187 * t220 + t188 * t221;
t80 = (t155 * t270 + t156 * t272) * t264 + t310;
t79 = -t150 * t327 + t152 * t243 + t154 * t244;
t78 = -t149 * t327 + t151 * t243 + t153 * t244;
t77 = t276 + t108;
t76 = -t106 + t273;
t75 = t108 * t231 - t133 * t216;
t74 = -t106 * t231 + t133 * t214;
t71 = t138 * t247 - t140 * t216 + t142 * t217;
t70 = t137 * t247 - t139 * t216 + t141 * t217;
t69 = t138 * t245 - t140 * t214 + t142 * t215;
t68 = t137 * t245 - t139 * t214 + t141 * t215;
t67 = (-t143 + t316) * t266 + t272 * t286;
t66 = t266 * t144 + t270 * t286 + t317;
t65 = t106 * t216 - t108 * t214;
t58 = t276 + t320;
t57 = t273 - t321;
t56 = t319 * t327 + (-t133 - t190) * t247;
t55 = t106 * t327 + t133 * t245 + t314;
t50 = (t143 * t270 + t144 * t272) * t264 + t287;
t49 = t106 * t247 + t245 * t319 + t147;
t48 = -t216 * t318 + t231 * t320;
t47 = t214 * t318 - t231 * t321;
t42 = (-t106 + t294) * t266 + t272 * t282;
t41 = t266 * t108 + t270 * t282 + t295;
t32 = t296 * t327 + (-t190 - t318) * t247;
t31 = t245 * t318 + t321 * t327 + t314;
t30 = -t214 * t320 + t216 * t321;
t29 = (t106 * t270 + t108 * t272) * t264 + t279;
t28 = t245 * t296 + t247 * t321 + t147;
t27 = (t294 - t321) * t266 + t272 * t278;
t26 = t266 * t320 + t270 * t278 + t295;
t25 = t85 + (t73 * t270 - t72 * t272) * t264;
t24 = t72 * t245 + t73 * t247 - t327 * t87;
t23 = t82 * t266 + (t270 * t71 - t272 * t70) * t264;
t22 = t81 * t266 + (t270 * t69 - t272 * t68) * t264;
t21 = t245 * t70 + t247 * t71 - t327 * t82;
t20 = t245 * t68 + t247 * t69 - t327 * t81;
t19 = (t270 * t321 + t272 * t320) * t264 + t279;
t90 = [(-t182 - t186) * t327 + m(7) * (t57 ^ 2 + t58 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t109 ^ 2 + t110 ^ 2) + m(4) * (t117 ^ 2 + t118 ^ 2) + m(3) * (t180 ^ 2 + t181 ^ 2) + m(2) * (t253 ^ 2 + t254 ^ 2) + Icges(2,3) + t313 - t336 + t339; t61 + t62 + t85 + m(6) * (t41 * t77 + t42 * t76) + m(7) * (t26 * t58 + t27 * t57) + m(5) * (t109 * t67 + t110 * t66) + m(4) * (t117 * t89 + t118 * t88) + m(3) * (t163 * t181 + t164 * t180) + ((-t111 / 0.2e1 - t78 / 0.2e1 - t83 / 0.2e1 - t123 / 0.2e1 - t275) * t272 + (t112 / 0.2e1 + t79 / 0.2e1 + t84 / 0.2e1 + t124 / 0.2e1 + t274) * t270) * t264 + t335; m(7) * (t19 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t29 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(5) * (t50 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(4) * (t80 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(3) * (t125 ^ 2 + t163 ^ 2 + t164 ^ 2) + (t12 + t11 + t23 + ((t150 * t247 + t152 * t222 + t154 * t223) * t270 - (t149 * t247 + t151 * t222 + t153 * t223) * t272) * t264 + (t193 * t328 - t195 * t247 + t248 * t197) * t328) * t328 + (-t10 - t9 - t22 - ((t150 * t245 + t152 * t220 + t154 * t221) * t270 - (t149 * t245 + t151 * t220 + t153 * t221) * t272) * t264 + (-t192 * t326 - t245 * t194 + t246 * t196) * t326 + (-t192 * t328 + t193 * t326 + t247 * t194 + t245 * t195 - t248 * t196 - t246 * t197) * t328) * t326 + (t17 + t18 + t25 + (t84 + t124) * t328 + (-t123 - t83) * t326 + ((-t111 - t78) * t272 + (t112 + t79) * t270) * t264 + t335) * t266; m(7) * (t245 * t58 + t247 * t57) + m(6) * (t245 * t77 + t247 * t76) + m(5) * (t109 * t247 + t110 * t245) + m(4) * (t117 * t247 + t118 * t245); m(7) * (-t19 * t327 + t245 * t26 + t247 * t27) + m(6) * (t245 * t41 + t247 * t42 - t29 * t327) + m(5) * (t245 * t66 + t247 * t67 - t327 * t50) + m(4) * (t245 * t88 + t247 * t89 - t327 * t80); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1 + m(4) / 0.2e1) * (t264 ^ 2 * t271 ^ 2 + t245 ^ 2 + t247 ^ 2); (-t87 + t336) * t327 + m(7) * (t31 * t57 + t32 * t58) + m(6) * (t55 * t76 + t56 * t77) + m(5) * (t109 * t91 + t110 * t92) + t274 * t247 + t275 * t245; (t24 / 0.2e1 + t300) * t266 + (t23 / 0.2e1 + t302) * t247 + (t22 / 0.2e1 + t303) * t245 + m(7) * (t19 * t28 + t26 * t32 + t27 * t31) + m(6) * (t29 * t49 + t41 * t56 + t42 * t55) + m(5) * (t50 * t86 + t66 * t92 + t67 * t91) + ((-t20 / 0.2e1 - t305) * t272 + (-t25 / 0.2e1 - t299) * t271 + (t21 / 0.2e1 + t304) * t270) * t264; m(5) * (t245 * t92 + t247 * t91 - t327 * t86) + m(6) * (t245 * t56 + t247 * t55 - t327 * t49) + m(7) * (t245 * t32 + t247 * t31 - t28 * t327); (-t15 - t16 - t24) * t327 + (t7 + t8 + t21) * t247 + (t5 + t6 + t20) * t245 + m(7) * (t28 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t49 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t86 ^ 2 + t91 ^ 2 + t92 ^ 2); t59 + t60 + m(7) * (t47 * t57 + t48 * t58) + m(6) * (t74 * t76 + t75 * t77) + t280 * t216 + t281 * t214; t301 * t266 + t299 * t231 + t302 * t216 + t303 * t214 + m(7) * (t19 * t30 + t26 * t48 + t27 * t47) + m(6) * (t29 * t65 + t41 * t75 + t42 * t74) + (t270 * t306 - t272 * t307) * t264; m(6) * (t245 * t75 + t247 * t74 - t327 * t65) + m(7) * (t245 * t48 + t247 * t47 - t30 * t327); -t301 * t327 + t306 * t247 + t307 * t245 + t300 * t231 + t304 * t216 + t305 * t214 + m(7) * (t28 * t30 + t31 * t47 + t32 * t48) + m(6) * (t49 * t65 + t55 * t74 + t56 * t75); (t13 + t14) * t231 + (t4 + t3) * t216 + (t1 + t2) * t214 + m(7) * (t30 ^ 2 + t47 ^ 2 + t48 ^ 2) + m(6) * (t65 ^ 2 + t74 ^ 2 + t75 ^ 2); m(7) * (t176 * t58 + t178 * t57); m(7) * (t176 * t26 + t178 * t27 + t19 * t212); m(7) * (t176 * t245 + t178 * t247 - t212 * t327); m(7) * (t176 * t32 + t178 * t31 + t212 * t28); m(7) * (t176 * t48 + t178 * t47 + t212 * t30); m(7) * (t176 ^ 2 + t178 ^ 2 + t212 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t90(1) t90(2) t90(4) t90(7) t90(11) t90(16); t90(2) t90(3) t90(5) t90(8) t90(12) t90(17); t90(4) t90(5) t90(6) t90(9) t90(13) t90(18); t90(7) t90(8) t90(9) t90(10) t90(14) t90(19); t90(11) t90(12) t90(13) t90(14) t90(15) t90(20); t90(16) t90(17) t90(18) t90(19) t90(20) t90(21);];
Mq  = res;
