% Calculate joint inertia matrix for
% S6RRPRRP9
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:34
% EndTime: 2019-03-09 12:27:50
% DurationCPUTime: 6.34s
% Computational Cost: add. (25575->633), mult. (42198->871), div. (0->0), fcn. (53334->12), ass. (0->278)
t343 = rSges(7,3) + qJ(6) + pkin(10);
t264 = sin(pkin(11));
t266 = cos(pkin(11));
t267 = cos(pkin(6));
t265 = sin(pkin(6));
t271 = sin(qJ(2));
t331 = t265 * t271;
t243 = -t264 * t331 + t266 * t267;
t332 = t264 * t267;
t244 = t266 * t331 + t332;
t274 = cos(qJ(2));
t329 = t265 * t274;
t186 = Icges(4,4) * t244 + Icges(4,2) * t243 - Icges(4,6) * t329;
t187 = Icges(4,1) * t244 + Icges(4,4) * t243 - Icges(4,5) * t329;
t229 = Icges(3,3) * t267 + (Icges(3,5) * t271 + Icges(3,6) * t274) * t265;
t230 = Icges(3,6) * t267 + (Icges(3,4) * t271 + Icges(3,2) * t274) * t265;
t231 = Icges(3,5) * t267 + (Icges(3,1) * t271 + Icges(3,4) * t274) * t265;
t345 = t243 * t186 + t244 * t187 + t267 * t229 + t230 * t329 + t231 * t331;
t272 = sin(qJ(1));
t325 = t272 * t274;
t275 = cos(qJ(1));
t326 = t271 * t275;
t246 = t267 * t326 + t325;
t312 = pkin(11) + qJ(4);
t262 = sin(t312);
t291 = cos(t312);
t328 = t265 * t275;
t216 = t246 * t291 - t262 * t328;
t324 = t274 * t275;
t327 = t271 * t272;
t245 = -t267 * t324 + t327;
t270 = sin(qJ(5));
t273 = cos(qJ(5));
t175 = -t216 * t270 + t245 * t273;
t334 = t245 * t270;
t176 = t216 * t273 + t334;
t285 = t265 * t291;
t215 = t246 * t262 + t275 * t285;
t344 = -rSges(7,1) * t176 - rSges(7,2) * t175 - t343 * t215;
t233 = t267 * t262 + t271 * t285;
t213 = -t233 * t270 - t273 * t329;
t299 = t270 * t329;
t214 = t233 * t273 - t299;
t232 = t262 * t331 - t267 * t291;
t127 = Icges(7,5) * t214 + Icges(7,6) * t213 + Icges(7,3) * t232;
t129 = Icges(7,4) * t214 + Icges(7,2) * t213 + Icges(7,6) * t232;
t131 = Icges(7,1) * t214 + Icges(7,4) * t213 + Icges(7,5) * t232;
t61 = t232 * t127 + t213 * t129 + t214 * t131;
t128 = Icges(6,5) * t214 + Icges(6,6) * t213 + Icges(6,3) * t232;
t130 = Icges(6,4) * t214 + Icges(6,2) * t213 + Icges(6,6) * t232;
t132 = Icges(6,1) * t214 + Icges(6,4) * t213 + Icges(6,5) * t232;
t62 = t232 * t128 + t213 * t130 + t214 * t132;
t342 = -t61 - t62;
t185 = Icges(4,5) * t244 + Icges(4,6) * t243 - Icges(4,3) * t329;
t341 = (-t185 * t329 + t345) * t267;
t248 = -t267 * t327 + t324;
t330 = t265 * t272;
t218 = t248 * t291 + t262 * t330;
t247 = t267 * t325 + t326;
t177 = -t218 * t270 + t247 * t273;
t333 = t247 * t270;
t178 = t218 * t273 + t333;
t217 = t248 * t262 - t272 * t285;
t261 = pkin(5) * t273 + pkin(4);
t340 = t178 * rSges(7,1) + t177 * rSges(7,2) + pkin(5) * t333 + t217 * t343 + t218 * t261;
t260 = pkin(3) * t266 + pkin(2);
t339 = -pkin(2) + t260;
t338 = -pkin(4) + t261;
t211 = t215 * pkin(10);
t337 = pkin(5) * t334 + t216 * t338 - t211 - t344;
t167 = t218 * pkin(4) + pkin(10) * t217;
t336 = -t167 + t340;
t269 = -pkin(9) - qJ(3);
t335 = t245 * t269;
t110 = t178 * rSges(6,1) + t177 * rSges(6,2) + t217 * rSges(6,3);
t323 = -t110 - t167;
t322 = rSges(7,1) * t214 + rSges(7,2) * t213 - pkin(5) * t299 + t338 * t233 + (-pkin(10) + t343) * t232;
t206 = t248 * pkin(2) + qJ(3) * t247;
t301 = t264 * t330;
t293 = pkin(3) * t301 - t247 * t269 + t248 * t260;
t147 = -t206 + t293;
t203 = t267 * t206;
t321 = t267 * t147 + t203;
t237 = t245 * qJ(3);
t300 = t264 * t328;
t252 = pkin(3) * t300;
t146 = t246 * t339 - t237 - t252 - t335;
t205 = pkin(2) * t246 + t237;
t320 = -t146 - t205;
t221 = -t246 * t264 - t266 * t328;
t222 = t246 * t266 - t300;
t156 = rSges(4,1) * t222 + rSges(4,2) * t221 + rSges(4,3) * t245;
t319 = -t156 - t205;
t166 = pkin(4) * t216 + t211;
t189 = t233 * pkin(4) + t232 * pkin(10);
t318 = t166 * t329 + t245 * t189;
t182 = Icges(5,4) * t233 - Icges(5,2) * t232 - Icges(5,6) * t329;
t183 = Icges(5,1) * t233 - Icges(5,4) * t232 - Icges(5,5) * t329;
t317 = -t232 * t182 + t233 * t183;
t249 = (pkin(2) * t271 - qJ(3) * t274) * t265;
t315 = -pkin(3) * t332 - ((qJ(3) + t269) * t274 + t339 * t271) * t265 - t249;
t314 = t205 * t330 + t206 * t328;
t313 = t275 * pkin(1) + pkin(8) * t330;
t103 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t215;
t95 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t215;
t99 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t215;
t35 = t103 * t176 + t175 * t99 + t215 * t95;
t100 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t217;
t104 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t217;
t96 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t217;
t36 = t100 * t175 + t104 * t176 + t215 * t96;
t51 = t127 * t215 + t129 * t175 + t131 * t176;
t1 = t215 * t35 + t217 * t36 + t232 * t51;
t101 = Icges(6,4) * t176 + Icges(6,2) * t175 + Icges(6,6) * t215;
t105 = Icges(6,1) * t176 + Icges(6,4) * t175 + Icges(6,5) * t215;
t97 = Icges(6,5) * t176 + Icges(6,6) * t175 + Icges(6,3) * t215;
t37 = t101 * t175 + t105 * t176 + t215 * t97;
t102 = Icges(6,4) * t178 + Icges(6,2) * t177 + Icges(6,6) * t217;
t106 = Icges(6,1) * t178 + Icges(6,4) * t177 + Icges(6,5) * t217;
t98 = Icges(6,5) * t178 + Icges(6,6) * t177 + Icges(6,3) * t217;
t38 = t102 * t175 + t106 * t176 + t215 * t98;
t52 = t128 * t215 + t130 * t175 + t132 * t176;
t2 = t215 * t37 + t217 * t38 + t232 * t52;
t311 = t2 / 0.2e1 + t1 / 0.2e1;
t39 = t103 * t178 + t177 * t99 + t217 * t95;
t40 = t100 * t177 + t104 * t178 + t217 * t96;
t53 = t127 * t217 + t129 * t177 + t131 * t178;
t3 = t215 * t39 + t217 * t40 + t232 * t53;
t41 = t101 * t177 + t105 * t178 + t217 * t97;
t42 = t102 * t177 + t106 * t178 + t217 * t98;
t54 = t128 * t217 + t130 * t177 + t132 * t178;
t4 = t215 * t41 + t217 * t42 + t232 * t54;
t310 = t4 / 0.2e1 + t3 / 0.2e1;
t5 = t245 * t35 + t247 * t36 - t329 * t51;
t6 = t245 * t37 + t247 * t38 - t329 * t52;
t309 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = t245 * t39 + t247 * t40 - t329 * t53;
t8 = t245 * t41 + t247 * t42 - t329 * t54;
t308 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t52 * t267 + (t272 * t38 - t275 * t37) * t265;
t9 = t51 * t267 + (t272 * t36 - t275 * t35) * t265;
t307 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t53 * t267 + (t272 * t40 - t275 * t39) * t265;
t12 = t54 * t267 + (t272 * t42 - t275 * t41) * t265;
t306 = t12 / 0.2e1 + t11 / 0.2e1;
t45 = t103 * t214 + t213 * t99 + t232 * t95;
t46 = t100 * t213 + t104 * t214 + t232 * t96;
t57 = t61 * t232;
t13 = t45 * t215 + t46 * t217 + t57;
t47 = t101 * t213 + t105 * t214 + t232 * t97;
t48 = t102 * t213 + t106 * t214 + t232 * t98;
t58 = t62 * t232;
t14 = t47 * t215 + t48 * t217 + t58;
t305 = t14 / 0.2e1 + t13 / 0.2e1;
t15 = t45 * t245 + t46 * t247 - t329 * t61;
t16 = t47 * t245 + t48 * t247 - t329 * t62;
t304 = t16 / 0.2e1 + t15 / 0.2e1;
t59 = t61 * t267;
t17 = t59 + (t46 * t272 - t45 * t275) * t265;
t60 = t62 * t267;
t18 = t60 + (t48 * t272 - t47 * t275) * t265;
t303 = t18 / 0.2e1 + t17 / 0.2e1;
t302 = -t167 - t336;
t298 = t267 * t167 + t321;
t297 = -t166 + t320;
t296 = -t189 + t315;
t145 = t218 * rSges(5,1) - t217 * rSges(5,2) + t247 * rSges(5,3);
t223 = -t248 * t264 + t266 * t330;
t224 = t248 * t266 + t301;
t157 = t224 * rSges(4,1) + t223 * rSges(4,2) + t247 * rSges(4,3);
t198 = t248 * rSges(3,1) - t247 * rSges(3,2) + rSges(3,3) * t330;
t292 = -t272 * pkin(1) + pkin(8) * t328;
t290 = t265 * (-rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t329 - t249);
t289 = t146 * t330 + t147 * t328 + t314;
t184 = rSges(5,1) * t233 - rSges(5,2) * t232 - rSges(5,3) * t329;
t288 = t265 * (-t184 + t315);
t287 = -rSges(5,1) * t216 + rSges(5,2) * t215;
t284 = t293 + t313;
t134 = rSges(6,1) * t214 + rSges(6,2) * t213 + rSges(6,3) * t232;
t283 = t265 * (-t134 + t296);
t282 = t48 / 0.2e1 + t46 / 0.2e1 + t53 / 0.2e1 + t54 / 0.2e1;
t281 = t52 / 0.2e1 + t47 / 0.2e1 + t45 / 0.2e1 + t51 / 0.2e1;
t280 = t166 * t330 + t167 * t328 + t289;
t279 = t265 * (t296 - t322);
t278 = -t246 * t260 + t252 + t292;
t108 = rSges(6,1) * t176 + rSges(6,2) * t175 + rSges(6,3) * t215;
t197 = t246 * rSges(3,1) - t245 * rSges(3,2) - rSges(3,3) * t328;
t139 = Icges(5,5) * t218 - Icges(5,6) * t217 + Icges(5,3) * t247;
t141 = Icges(5,4) * t218 - Icges(5,2) * t217 + Icges(5,6) * t247;
t143 = Icges(5,1) * t218 - Icges(5,4) * t217 + Icges(5,5) * t247;
t71 = -t139 * t329 - t141 * t232 + t143 * t233;
t181 = Icges(5,5) * t233 - Icges(5,6) * t232 - Icges(5,3) * t329;
t82 = t181 * t247 - t182 * t217 + t183 * t218;
t277 = t71 / 0.2e1 + t82 / 0.2e1 + t282;
t138 = Icges(5,5) * t216 - Icges(5,6) * t215 + Icges(5,3) * t245;
t140 = Icges(5,4) * t216 - Icges(5,2) * t215 + Icges(5,6) * t245;
t142 = Icges(5,1) * t216 - Icges(5,4) * t215 + Icges(5,5) * t245;
t70 = -t138 * t329 - t140 * t232 + t142 * t233;
t81 = t181 * t245 - t182 * t215 + t183 * t216;
t276 = t70 / 0.2e1 + t81 / 0.2e1 + t281;
t254 = rSges(2,1) * t275 - t272 * rSges(2,2);
t253 = -t272 * rSges(2,1) - rSges(2,2) * t275;
t234 = t267 * rSges(3,3) + (rSges(3,1) * t271 + rSges(3,2) * t274) * t265;
t196 = Icges(3,1) * t248 - Icges(3,4) * t247 + Icges(3,5) * t330;
t195 = Icges(3,1) * t246 - Icges(3,4) * t245 - Icges(3,5) * t328;
t194 = Icges(3,4) * t248 - Icges(3,2) * t247 + Icges(3,6) * t330;
t193 = Icges(3,4) * t246 - Icges(3,2) * t245 - Icges(3,6) * t328;
t192 = Icges(3,5) * t248 - Icges(3,6) * t247 + Icges(3,3) * t330;
t191 = Icges(3,5) * t246 - Icges(3,6) * t245 - Icges(3,3) * t328;
t180 = t198 + t313;
t179 = -t197 + t292;
t165 = -t267 * t197 - t234 * t328;
t164 = t198 * t267 - t234 * t330;
t155 = Icges(4,1) * t224 + Icges(4,4) * t223 + Icges(4,5) * t247;
t154 = Icges(4,1) * t222 + Icges(4,4) * t221 + Icges(4,5) * t245;
t153 = Icges(4,4) * t224 + Icges(4,2) * t223 + Icges(4,6) * t247;
t152 = Icges(4,4) * t222 + Icges(4,2) * t221 + Icges(4,6) * t245;
t151 = Icges(4,5) * t224 + Icges(4,6) * t223 + Icges(4,3) * t247;
t150 = Icges(4,5) * t222 + Icges(4,6) * t221 + Icges(4,3) * t245;
t148 = t247 * t166;
t144 = rSges(5,3) * t245 - t287;
t126 = (t197 * t272 + t198 * t275) * t265;
t125 = t229 * t330 - t230 * t247 + t231 * t248;
t124 = -t229 * t328 - t245 * t230 + t246 * t231;
t120 = t206 + t157 + t313;
t119 = t292 + t319;
t114 = t267 * t192 + (t194 * t274 + t196 * t271) * t265;
t113 = t267 * t191 + (t193 * t274 + t195 * t271) * t265;
t112 = t284 + t145;
t111 = (-rSges(5,3) + t269) * t245 + t278 + t287;
t94 = -t145 * t329 - t184 * t247;
t93 = t144 * t329 + t184 * t245;
t89 = t267 * t319 + t275 * t290;
t88 = t267 * t157 + t272 * t290 + t203;
t87 = -t181 * t329 + t317;
t86 = t144 * t247 - t145 * t245;
t85 = t87 * t267;
t84 = t185 * t247 + t186 * t223 + t187 * t224;
t83 = t185 * t245 + t186 * t221 + t187 * t222;
t80 = (t156 * t272 + t157 * t275) * t265 + t314;
t79 = -t151 * t329 + t153 * t243 + t155 * t244;
t78 = -t150 * t329 + t152 * t243 + t154 * t244;
t77 = t284 - t323;
t76 = -t108 - t166 + t278 + t335;
t75 = t110 * t232 - t134 * t217;
t74 = -t108 * t232 + t134 * t215;
t73 = t284 + t340;
t72 = -t216 * t261 + (-pkin(5) * t270 + t269) * t245 + t278 + t344;
t69 = t139 * t247 - t141 * t217 + t143 * t218;
t68 = t138 * t247 - t140 * t217 + t142 * t218;
t67 = t139 * t245 - t141 * t215 + t143 * t216;
t66 = t138 * t245 - t140 * t215 + t142 * t216;
t65 = (-t144 + t320) * t267 + t275 * t288;
t64 = t267 * t145 + t272 * t288 + t321;
t63 = t108 * t217 - t110 * t215;
t56 = t323 * t329 + (-t134 - t189) * t247;
t55 = t108 * t329 + t134 * t245 + t318;
t50 = (t144 * t272 + t145 * t275) * t265 + t289;
t49 = t247 * t108 + t245 * t323 + t148;
t44 = (-t108 + t297) * t267 + t275 * t283;
t43 = t110 * t267 + t272 * t283 + t298;
t34 = -t217 * t322 + t232 * t336;
t33 = t215 * t322 - t232 * t337;
t32 = t302 * t329 + (-t189 - t322) * t247;
t31 = t245 * t322 + t329 * t337 + t318;
t30 = (t108 * t272 + t110 * t275) * t265 + t280;
t29 = -t215 * t336 + t217 * t337;
t28 = t85 + (t71 * t272 - t70 * t275) * t265;
t27 = (t297 - t337) * t267 + t275 * t279;
t26 = t267 * t336 + t272 * t279 + t298;
t25 = t245 * t302 + t247 * t337 + t148;
t24 = t70 * t245 + t71 * t247 - t329 * t87;
t23 = t82 * t267 + (t272 * t69 - t275 * t68) * t265;
t22 = t81 * t267 + (t272 * t67 - t275 * t66) * t265;
t21 = t245 * t68 + t247 * t69 - t329 * t82;
t20 = t245 * t66 + t247 * t67 - t329 * t81;
t19 = (t272 * t337 + t275 * t336) * t265 + t280;
t90 = [(-t181 - t185) * t329 + m(7) * (t72 ^ 2 + t73 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t111 ^ 2 + t112 ^ 2) + m(4) * (t119 ^ 2 + t120 ^ 2) + m(3) * (t179 ^ 2 + t180 ^ 2) + m(2) * (t253 ^ 2 + t254 ^ 2) + Icges(2,3) + t317 - t342 + t345; t85 + t60 + t59 + m(7) * (t26 * t73 + t27 * t72) + m(6) * (t43 * t77 + t44 * t76) + m(5) * (t111 * t65 + t112 * t64) + m(4) * (t119 * t89 + t120 * t88) + m(3) * (t164 * t180 + t165 * t179) + ((-t113 / 0.2e1 - t78 / 0.2e1 - t83 / 0.2e1 - t124 / 0.2e1 - t276) * t275 + (t84 / 0.2e1 + t125 / 0.2e1 + t114 / 0.2e1 + t79 / 0.2e1 + t277) * t272) * t265 + t341; m(7) * (t19 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t30 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(5) * (t50 ^ 2 + t64 ^ 2 + t65 ^ 2) + m(4) * (t80 ^ 2 + t88 ^ 2 + t89 ^ 2) + m(3) * (t126 ^ 2 + t164 ^ 2 + t165 ^ 2) + (t12 + t11 + t23 + ((t151 * t247 + t223 * t153 + t224 * t155) * t272 - (t247 * t150 + t223 * t152 + t224 * t154) * t275) * t265 + (t192 * t330 - t194 * t247 + t248 * t196) * t330) * t330 + (-t10 - t9 - t22 - ((t245 * t151 + t221 * t153 + t222 * t155) * t272 - (t150 * t245 + t152 * t221 + t154 * t222) * t275) * t265 + (-t191 * t328 - t245 * t193 + t246 * t195) * t328 + (-t191 * t330 + t192 * t328 + t247 * t193 + t245 * t194 - t248 * t195 - t246 * t196) * t330) * t328 + (t18 + t17 + t28 + (t84 + t125) * t330 + (-t83 - t124) * t328 + ((-t113 - t78) * t275 + (t114 + t79) * t272) * t265 + t341) * t267; m(7) * (t245 * t73 + t247 * t72) + m(6) * (t245 * t77 + t247 * t76) + m(5) * (t111 * t247 + t112 * t245) + m(4) * (t119 * t247 + t120 * t245); m(7) * (-t19 * t329 + t245 * t26 + t247 * t27) + m(6) * (t245 * t43 + t247 * t44 - t30 * t329) + m(5) * (t245 * t64 + t247 * t65 - t329 * t50) + m(4) * (t245 * t88 + t247 * t89 - t329 * t80); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1 + m(4) / 0.2e1) * (t265 ^ 2 * t274 ^ 2 + t245 ^ 2 + t247 ^ 2); (-t87 + t342) * t329 + m(7) * (t31 * t72 + t32 * t73) + m(6) * (t55 * t76 + t56 * t77) + m(5) * (t111 * t93 + t112 * t94) + t277 * t247 + t276 * t245; (t24 / 0.2e1 + t304) * t267 + (t23 / 0.2e1 + t306) * t247 + (t22 / 0.2e1 + t307) * t245 + m(7) * (t19 * t25 + t26 * t32 + t27 * t31) + m(6) * (t30 * t49 + t43 * t56 + t44 * t55) + m(5) * (t50 * t86 + t64 * t94 + t65 * t93) + ((-t20 / 0.2e1 - t309) * t275 + (-t28 / 0.2e1 - t303) * t274 + (t21 / 0.2e1 + t308) * t272) * t265; m(5) * (t245 * t94 + t247 * t93 - t329 * t86) + m(6) * (t245 * t56 + t247 * t55 - t329 * t49) + m(7) * (t245 * t32 + t247 * t31 - t25 * t329); (-t15 - t16 - t24) * t329 + (t8 + t7 + t21) * t247 + (t6 + t5 + t20) * t245 + m(7) * (t25 ^ 2 + t31 ^ 2 + t32 ^ 2) + m(6) * (t49 ^ 2 + t55 ^ 2 + t56 ^ 2) + m(5) * (t86 ^ 2 + t93 ^ 2 + t94 ^ 2); t58 + t57 + m(7) * (t33 * t72 + t34 * t73) + m(6) * (t74 * t76 + t75 * t77) + t282 * t217 + t281 * t215; t305 * t267 + t303 * t232 + t306 * t217 + t307 * t215 + m(7) * (t19 * t29 + t26 * t34 + t27 * t33) + m(6) * (t30 * t63 + t43 * t75 + t44 * t74) + (t272 * t310 - t275 * t311) * t265; m(6) * (t245 * t75 + t247 * t74 - t329 * t63) + m(7) * (t245 * t34 + t247 * t33 - t29 * t329); -t305 * t329 + t310 * t247 + t311 * t245 + t304 * t232 + t308 * t217 + t309 * t215 + m(7) * (t25 * t29 + t31 * t33 + t32 * t34) + m(6) * (t49 * t63 + t55 * t74 + t56 * t75); (t14 + t13) * t232 + (t4 + t3) * t217 + (t1 + t2) * t215 + m(7) * (t29 ^ 2 + t33 ^ 2 + t34 ^ 2) + m(6) * (t63 ^ 2 + t74 ^ 2 + t75 ^ 2); m(7) * (t215 * t73 + t217 * t72); m(7) * (t19 * t232 + t215 * t26 + t217 * t27); m(7) * (t215 * t245 + t217 * t247 - t232 * t329); m(7) * (t215 * t32 + t217 * t31 + t232 * t25); m(7) * (t215 * t34 + t217 * t33 + t232 * t29); m(7) * (t215 ^ 2 + t217 ^ 2 + t232 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t90(1) t90(2) t90(4) t90(7) t90(11) t90(16); t90(2) t90(3) t90(5) t90(8) t90(12) t90(17); t90(4) t90(5) t90(6) t90(9) t90(13) t90(18); t90(7) t90(8) t90(9) t90(10) t90(14) t90(19); t90(11) t90(12) t90(13) t90(14) t90(15) t90(20); t90(16) t90(17) t90(18) t90(19) t90(20) t90(21);];
Mq  = res;
