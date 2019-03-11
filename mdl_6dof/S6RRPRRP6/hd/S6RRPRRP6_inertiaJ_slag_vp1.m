% Calculate joint inertia matrix for
% S6RRPRRP6
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:50
% EndTime: 2019-03-09 12:08:07
% DurationCPUTime: 7.02s
% Computational Cost: add. (27807->619), mult. (71372->868), div. (0->0), fcn. (94006->12), ass. (0->277)
t270 = cos(pkin(6));
t329 = sin(pkin(11));
t330 = cos(pkin(11));
t333 = sin(qJ(2));
t336 = cos(qJ(2));
t277 = t329 * t336 + t330 * t333;
t246 = t277 * t270;
t256 = -t333 * t329 + t336 * t330;
t273 = sin(qJ(1));
t274 = cos(qJ(1));
t230 = t246 * t274 + t273 * t256;
t272 = sin(qJ(4));
t269 = sin(pkin(6));
t326 = t269 * t274;
t335 = cos(qJ(4));
t200 = t230 * t335 - t272 * t326;
t275 = t270 * t256;
t229 = -t273 * t277 + t274 * t275;
t271 = sin(qJ(5));
t334 = cos(qJ(5));
t161 = t200 * t271 + t229 * t334;
t162 = t200 * t334 - t229 * t271;
t300 = t269 * t335;
t199 = t230 * t272 + t274 * t300;
t343 = rSges(7,3) + qJ(6);
t344 = rSges(7,1) + pkin(5);
t324 = t199 * rSges(7,2) + t343 * t161 + t162 * t344;
t244 = t256 * t269;
t245 = t277 * t269;
t206 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t270;
t207 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t270;
t208 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t270;
t240 = Icges(3,3) * t270 + (Icges(3,5) * t333 + Icges(3,6) * t336) * t269;
t241 = Icges(3,6) * t270 + (Icges(3,4) * t333 + Icges(3,2) * t336) * t269;
t242 = Icges(3,5) * t270 + (Icges(3,1) * t333 + Icges(3,4) * t336) * t269;
t299 = t269 * t333;
t342 = t269 * t336 * t241 + t244 * t207 + t245 * t208 + t242 * t299 + (t206 + t240) * t270;
t341 = t269 ^ 2;
t340 = m(4) / 0.2e1;
t339 = m(5) / 0.2e1;
t338 = m(6) / 0.2e1;
t337 = m(7) / 0.2e1;
t332 = pkin(1) * t274;
t247 = t270 * t333 * pkin(2) + (-pkin(8) - qJ(3)) * t269;
t328 = t247 * t274;
t327 = t269 * t273;
t267 = pkin(2) * t336 + pkin(1);
t325 = t273 * t267;
t232 = -t273 * t246 + t256 * t274;
t202 = t232 * t335 + t272 * t327;
t231 = -t273 * t275 - t274 * t277;
t163 = t202 * t271 + t231 * t334;
t164 = t202 * t334 - t231 * t271;
t201 = t232 * t272 - t273 * t300;
t323 = t201 * rSges(7,2) + t343 * t163 + t164 * t344;
t322 = t342 * t270;
t234 = t245 * t335 + t270 * t272;
t189 = t234 * t271 + t244 * t334;
t190 = t234 * t334 - t244 * t271;
t233 = t245 * t272 - t270 * t335;
t321 = rSges(7,2) * t233 + t343 * t189 + t190 * t344;
t182 = t232 * pkin(3) - pkin(9) * t231;
t261 = t274 * t267;
t223 = -t332 + t261 + (-t269 * pkin(8) - t247) * t273;
t210 = t270 * t223;
t320 = t270 * t182 + t210;
t181 = t230 * pkin(3) - t229 * pkin(9);
t265 = pkin(8) * t326;
t222 = t328 + t265 + (-pkin(1) + t267) * t273;
t319 = -t181 - t222;
t318 = t222 * t327 + t223 * t326;
t166 = Icges(4,5) * t230 + Icges(4,6) * t229 - Icges(4,3) * t326;
t296 = t274 * t336;
t297 = t273 * t333;
t251 = t270 * t296 - t297;
t295 = t274 * t333;
t298 = t273 * t336;
t252 = t270 * t295 + t298;
t211 = Icges(3,5) * t252 + Icges(3,6) * t251 - Icges(3,3) * t326;
t317 = -t211 - t166;
t167 = Icges(4,5) * t232 + Icges(4,6) * t231 + Icges(4,3) * t327;
t253 = -t270 * t298 - t295;
t254 = -t270 * t297 + t296;
t212 = Icges(3,5) * t254 + Icges(3,6) * t253 + Icges(3,3) * t327;
t316 = t212 + t167;
t257 = pkin(2) * t299 + t270 * qJ(3);
t315 = -pkin(3) * t245 + pkin(9) * t244 - t257;
t101 = Icges(7,1) * t162 + Icges(7,4) * t199 + Icges(7,5) * t161;
t93 = Icges(7,5) * t162 + Icges(7,6) * t199 + Icges(7,3) * t161;
t97 = Icges(7,4) * t162 + Icges(7,2) * t199 + Icges(7,6) * t161;
t33 = t101 * t162 + t161 * t93 + t199 * t97;
t102 = Icges(7,1) * t164 + Icges(7,4) * t201 + Icges(7,5) * t163;
t94 = Icges(7,5) * t164 + Icges(7,6) * t201 + Icges(7,3) * t163;
t98 = Icges(7,4) * t164 + Icges(7,2) * t201 + Icges(7,6) * t163;
t34 = t102 * t162 + t161 * t94 + t199 * t98;
t131 = Icges(7,5) * t190 + Icges(7,6) * t233 + Icges(7,3) * t189;
t133 = Icges(7,4) * t190 + Icges(7,2) * t233 + Icges(7,6) * t189;
t135 = Icges(7,1) * t190 + Icges(7,4) * t233 + Icges(7,5) * t189;
t52 = t131 * t161 + t133 * t199 + t135 * t162;
t1 = t199 * t33 + t201 * t34 + t233 * t52;
t103 = Icges(6,1) * t162 - Icges(6,4) * t161 + Icges(6,5) * t199;
t95 = Icges(6,5) * t162 - Icges(6,6) * t161 + Icges(6,3) * t199;
t99 = Icges(6,4) * t162 - Icges(6,2) * t161 + Icges(6,6) * t199;
t35 = t103 * t162 - t161 * t99 + t199 * t95;
t100 = Icges(6,4) * t164 - Icges(6,2) * t163 + Icges(6,6) * t201;
t104 = Icges(6,1) * t164 - Icges(6,4) * t163 + Icges(6,5) * t201;
t96 = Icges(6,5) * t164 - Icges(6,6) * t163 + Icges(6,3) * t201;
t36 = -t100 * t161 + t104 * t162 + t199 * t96;
t132 = Icges(6,5) * t190 - Icges(6,6) * t189 + Icges(6,3) * t233;
t134 = Icges(6,4) * t190 - Icges(6,2) * t189 + Icges(6,6) * t233;
t136 = Icges(6,1) * t190 - Icges(6,4) * t189 + Icges(6,5) * t233;
t53 = t132 * t199 - t134 * t161 + t136 * t162;
t2 = t199 * t35 + t201 * t36 + t233 * t53;
t314 = -t1 / 0.2e1 - t2 / 0.2e1;
t37 = t101 * t164 + t163 * t93 + t201 * t97;
t38 = t102 * t164 + t163 * t94 + t201 * t98;
t54 = t131 * t163 + t133 * t201 + t135 * t164;
t3 = t199 * t37 + t201 * t38 + t233 * t54;
t39 = t103 * t164 - t163 * t99 + t201 * t95;
t40 = -t100 * t163 + t104 * t164 + t201 * t96;
t55 = t132 * t201 - t134 * t163 + t136 * t164;
t4 = t199 * t39 + t201 * t40 + t233 * t55;
t313 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = -t229 * t33 - t231 * t34 - t244 * t52;
t6 = -t229 * t35 - t231 * t36 - t244 * t53;
t312 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = -t229 * t37 - t231 * t38 - t244 * t54;
t8 = -t229 * t39 - t231 * t40 - t244 * t55;
t311 = t8 / 0.2e1 + t7 / 0.2e1;
t10 = t53 * t270 + (t273 * t36 - t274 * t35) * t269;
t9 = t52 * t270 + (t273 * t34 - t274 * t33) * t269;
t310 = t9 / 0.2e1 + t10 / 0.2e1;
t11 = t54 * t270 + (t273 * t38 - t274 * t37) * t269;
t12 = t55 * t270 + (t273 * t40 - t274 * t39) * t269;
t309 = t11 / 0.2e1 + t12 / 0.2e1;
t41 = t101 * t190 + t189 * t93 + t233 * t97;
t42 = t102 * t190 + t189 * t94 + t233 * t98;
t70 = t189 * t131 + t233 * t133 + t190 * t135;
t59 = t70 * t233;
t13 = t41 * t199 + t42 * t201 + t59;
t43 = t103 * t190 - t189 * t99 + t233 * t95;
t44 = -t100 * t189 + t104 * t190 + t233 * t96;
t71 = t233 * t132 - t189 * t134 + t190 * t136;
t60 = t71 * t233;
t14 = t43 * t199 + t44 * t201 + t60;
t308 = t13 / 0.2e1 + t14 / 0.2e1;
t61 = t70 * t244;
t15 = -t41 * t229 - t42 * t231 - t61;
t62 = t71 * t244;
t16 = -t43 * t229 - t44 * t231 - t62;
t307 = t16 / 0.2e1 + t15 / 0.2e1;
t67 = t70 * t270;
t17 = t67 + (t42 * t273 - t41 * t274) * t269;
t68 = t71 * t270;
t18 = t68 + (t44 * t273 - t43 * t274) * t269;
t306 = t17 / 0.2e1 + t18 / 0.2e1;
t177 = Icges(5,5) * t234 - Icges(5,6) * t233 - Icges(5,3) * t244;
t178 = Icges(5,4) * t234 - Icges(5,2) * t233 - Icges(5,6) * t244;
t179 = Icges(5,1) * t234 - Icges(5,4) * t233 - Icges(5,5) * t244;
t85 = -t244 * t177 - t233 * t178 + t234 * t179;
t151 = t202 * pkin(4) + pkin(10) * t201;
t305 = t270 * t151 + t320;
t150 = t200 * pkin(4) + t199 * pkin(10);
t304 = -t150 + t319;
t108 = t164 * rSges(6,1) - t163 * rSges(6,2) + t201 * rSges(6,3);
t188 = pkin(4) * t234 + pkin(10) * t233;
t302 = -t188 + t315;
t129 = t202 * rSges(5,1) - t201 * rSges(5,2) - t231 * rSges(5,3);
t173 = t232 * rSges(4,1) + t231 * rSges(4,2) + rSges(4,3) * t327;
t220 = t254 * rSges(3,1) + t253 * rSges(3,2) + rSges(3,3) * t327;
t294 = t269 * (-rSges(4,1) * t245 - rSges(4,2) * t244 - rSges(4,3) * t270 - t257);
t293 = -t273 * t247 + t261;
t292 = t181 * t327 + t182 * t326 + t318;
t180 = rSges(5,1) * t234 - rSges(5,2) * t233 - rSges(5,3) * t244;
t291 = t269 * (-t180 + t315);
t288 = -t230 * rSges(4,1) - t229 * rSges(4,2);
t138 = rSges(6,1) * t190 - rSges(6,2) * t189 + rSges(6,3) * t233;
t287 = t269 * (-t138 + t302);
t286 = t43 / 0.2e1 + t41 / 0.2e1 + t53 / 0.2e1 + t52 / 0.2e1;
t285 = t44 / 0.2e1 + t42 / 0.2e1 + t55 / 0.2e1 + t54 / 0.2e1;
t284 = t150 * t327 + t151 * t326 + t292;
t283 = t269 * (t302 - t321);
t282 = t182 + t293;
t128 = t200 * rSges(5,1) - t199 * rSges(5,2) - t229 * rSges(5,3);
t106 = t162 * rSges(6,1) - t161 * rSges(6,2) + t199 * rSges(6,3);
t219 = t252 * rSges(3,1) + t251 * rSges(3,2) - rSges(3,3) * t326;
t281 = -t181 - t325 - t328;
t122 = Icges(5,5) * t200 - Icges(5,6) * t199 - Icges(5,3) * t229;
t124 = Icges(5,4) * t200 - Icges(5,2) * t199 - Icges(5,6) * t229;
t126 = Icges(5,1) * t200 - Icges(5,4) * t199 - Icges(5,5) * t229;
t72 = -t122 * t244 - t124 * t233 + t126 * t234;
t80 = -t177 * t229 - t178 * t199 + t179 * t200;
t280 = -t72 / 0.2e1 - t80 / 0.2e1 - t286;
t123 = Icges(5,5) * t202 - Icges(5,6) * t201 - Icges(5,3) * t231;
t125 = Icges(5,4) * t202 - Icges(5,2) * t201 - Icges(5,6) * t231;
t127 = Icges(5,1) * t202 - Icges(5,4) * t201 - Icges(5,5) * t231;
t73 = -t123 * t244 - t125 * t233 + t127 * t234;
t81 = -t177 * t231 - t178 * t201 + t179 * t202;
t279 = -t81 / 0.2e1 - t73 / 0.2e1 - t285;
t278 = t151 + t282;
t276 = -t150 + t281;
t259 = rSges(2,1) * t274 - t273 * rSges(2,2);
t258 = -t273 * rSges(2,1) - rSges(2,2) * t274;
t243 = t270 * rSges(3,3) + (rSges(3,1) * t333 + rSges(3,2) * t336) * t269;
t216 = Icges(3,1) * t254 + Icges(3,4) * t253 + Icges(3,5) * t327;
t215 = Icges(3,1) * t252 + Icges(3,4) * t251 - Icges(3,5) * t326;
t214 = Icges(3,4) * t254 + Icges(3,2) * t253 + Icges(3,6) * t327;
t213 = Icges(3,4) * t252 + Icges(3,2) * t251 - Icges(3,6) * t326;
t198 = pkin(8) * t327 + t220 + t332;
t197 = -t273 * pkin(1) - t219 + t265;
t185 = -t270 * t219 - t243 * t326;
t184 = t220 * t270 - t243 * t327;
t172 = -rSges(4,3) * t326 - t288;
t171 = Icges(4,1) * t232 + Icges(4,4) * t231 + Icges(4,5) * t327;
t170 = Icges(4,1) * t230 + Icges(4,4) * t229 - Icges(4,5) * t326;
t169 = Icges(4,4) * t232 + Icges(4,2) * t231 + Icges(4,6) * t327;
t168 = Icges(4,4) * t230 + Icges(4,2) * t229 - Icges(4,6) * t326;
t165 = (t219 * t273 + t220 * t274) * t269;
t160 = t240 * t327 + t241 * t253 + t242 * t254;
t159 = -t240 * t326 + t251 * t241 + t252 * t242;
t153 = t229 * t188;
t143 = t293 + t173;
t142 = -t325 + (rSges(4,3) * t269 - t247) * t274 + t288;
t141 = t244 * t151;
t140 = t270 * t212 + (t214 * t336 + t216 * t333) * t269;
t139 = t270 * t211 + (t213 * t336 + t215 * t333) * t269;
t130 = t231 * t150;
t112 = (-t172 - t222) * t270 + t274 * t294;
t111 = t173 * t270 + t273 * t294 + t210;
t110 = t206 * t327 + t207 * t231 + t208 * t232;
t109 = -t206 * t326 + t229 * t207 + t230 * t208;
t92 = t282 + t129;
t91 = -t128 + t281;
t90 = (t172 * t273 + t173 * t274) * t269 + t318;
t89 = -t129 * t244 + t180 * t231;
t88 = t128 * t244 - t180 * t229;
t87 = t167 * t270 + t169 * t244 + t171 * t245;
t86 = t166 * t270 + t168 * t244 + t170 * t245;
t84 = t85 * t270;
t83 = t85 * t244;
t82 = -t128 * t231 + t129 * t229;
t79 = (-t128 + t319) * t270 + t274 * t291;
t78 = t129 * t270 + t273 * t291 + t320;
t77 = t278 + t108;
t76 = -t106 + t276;
t75 = t108 * t233 - t138 * t201;
t74 = -t106 * t233 + t138 * t199;
t69 = (t128 * t273 + t129 * t274) * t269 + t292;
t66 = -t123 * t231 - t125 * t201 + t127 * t202;
t65 = -t122 * t231 - t124 * t201 + t126 * t202;
t64 = -t123 * t229 - t125 * t199 + t127 * t200;
t63 = -t122 * t229 - t124 * t199 + t126 * t200;
t58 = t106 * t201 - t108 * t199;
t57 = t278 + t323;
t56 = t276 - t324;
t51 = -t108 * t244 - t141 + (t138 + t188) * t231;
t50 = -t138 * t229 - t153 - (-t106 - t150) * t244;
t49 = (-t106 + t304) * t270 + t274 * t287;
t48 = t108 * t270 + t273 * t287 + t305;
t47 = -t201 * t321 + t233 * t323;
t46 = t199 * t321 - t233 * t324;
t45 = -t106 * t231 - t130 + (t108 + t151) * t229;
t32 = (t106 * t273 + t108 * t274) * t269 + t284;
t31 = -t141 - t323 * t244 + (t188 + t321) * t231;
t30 = -t153 - t321 * t229 - (-t150 - t324) * t244;
t29 = (t304 - t324) * t270 + t274 * t283;
t28 = t270 * t323 + t273 * t283 + t305;
t27 = -t199 * t323 + t201 * t324;
t26 = -t130 - t324 * t231 + (t151 + t323) * t229;
t25 = t84 + (t73 * t273 - t72 * t274) * t269;
t24 = (t273 * t324 + t274 * t323) * t269 + t284;
t23 = -t72 * t229 - t73 * t231 - t83;
t22 = t81 * t270 + (t273 * t66 - t274 * t65) * t269;
t21 = t80 * t270 + (t273 * t64 - t274 * t63) * t269;
t20 = -t229 * t65 - t231 * t66 - t244 * t81;
t19 = -t229 * t63 - t231 * t64 - t244 * t80;
t105 = [m(7) * (t56 ^ 2 + t57 ^ 2) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(5) * (t91 ^ 2 + t92 ^ 2) + m(4) * (t142 ^ 2 + t143 ^ 2) + m(3) * (t197 ^ 2 + t198 ^ 2) + m(2) * (t258 ^ 2 + t259 ^ 2) + Icges(2,3) + t85 + t71 + t70 + t342; t84 + t68 + t67 + m(7) * (t28 * t57 + t29 * t56) + m(6) * (t48 * t77 + t49 * t76) + m(5) * (t78 * t92 + t79 * t91) + m(4) * (t111 * t143 + t112 * t142) + m(3) * (t184 * t198 + t185 * t197) + ((-t86 / 0.2e1 - t139 / 0.2e1 - t109 / 0.2e1 - t159 / 0.2e1 + t280) * t274 + (t87 / 0.2e1 + t140 / 0.2e1 + t110 / 0.2e1 + t160 / 0.2e1 - t279) * t273) * t269 + t322; (t17 + t18 + t25 + t322) * t270 + m(7) * (t24 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t69 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t111 ^ 2 + t112 ^ 2 + t90 ^ 2) + m(3) * (t165 ^ 2 + t184 ^ 2 + t185 ^ 2) + ((-t10 - t21 - t9 + ((t229 * t168 + t230 * t170 + t251 * t213 + t252 * t215) * t269 + t317 * t341 * t274) * t274 + (-t109 - t139 - t159 - t86) * t270) * t274 + (t11 + t12 + t22 + ((t169 * t231 + t171 * t232 + t214 * t253 + t216 * t254) * t269 + t316 * t341 * t273) * t273 + (t160 + t110 + t87 + t140) * t270 + (-t168 * t231 - t229 * t169 - t170 * t232 - t230 * t171 - t213 * t253 - t251 * t214 - t215 * t254 - t252 * t216 + (t273 * t317 + t274 * t316) * t269) * t326) * t273) * t269; 0.2e1 * ((t273 * t56 - t274 * t57) * t337 + (t273 * t76 - t274 * t77) * t338 + (t273 * t91 - t274 * t92) * t339 + (t142 * t273 - t143 * t274) * t340) * t269; m(7) * (t270 * t24 + (t273 * t29 - t274 * t28) * t269) + m(6) * (t270 * t32 + (t273 * t49 - t274 * t48) * t269) + m(5) * (t270 * t69 + (t273 * t79 - t274 * t78) * t269) + m(4) * (t270 * t90 + (-t111 * t274 + t112 * t273) * t269); 0.2e1 * (t340 + t339 + t338 + t337) * (t270 ^ 2 + (t273 ^ 2 + t274 ^ 2) * t341); -t62 - t61 - t83 + m(7) * (t30 * t56 + t31 * t57) + m(6) * (t50 * t76 + t51 * t77) + m(5) * (t88 * t91 + t89 * t92) + t279 * t231 + t280 * t229; (t23 / 0.2e1 + t307) * t270 - (t25 / 0.2e1 + t306) * t244 + (-t22 / 0.2e1 - t309) * t231 + (-t21 / 0.2e1 - t310) * t229 + m(7) * (t24 * t26 + t28 * t31 + t29 * t30) + m(6) * (t32 * t45 + t48 * t51 + t49 * t50) + m(5) * (t69 * t82 + t78 * t89 + t79 * t88) + ((-t19 / 0.2e1 - t312) * t274 + (t20 / 0.2e1 + t311) * t273) * t269; m(5) * (t82 * t270 + (t273 * t88 - t274 * t89) * t269) + m(6) * (t45 * t270 + (t273 * t50 - t274 * t51) * t269) + m(7) * (t26 * t270 + (t273 * t30 - t274 * t31) * t269); -(t16 + t15 + t23) * t244 + (-t7 - t8 - t20) * t231 + (-t5 - t6 - t19) * t229 + m(7) * (t26 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t45 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t82 ^ 2 + t88 ^ 2 + t89 ^ 2); t60 + t59 + m(7) * (t46 * t56 + t47 * t57) + m(6) * (t74 * t76 + t75 * t77) + t285 * t201 + t286 * t199; t308 * t270 + t306 * t233 + t309 * t201 + t310 * t199 + m(7) * (t24 * t27 + t28 * t47 + t29 * t46) + m(6) * (t32 * t58 + t48 * t75 + t49 * t74) + (t273 * t313 + t274 * t314) * t269; m(6) * (t58 * t270 + (t273 * t74 - t274 * t75) * t269) + m(7) * (t27 * t270 + (t273 * t46 - t274 * t47) * t269); -t308 * t244 + t307 * t233 - t313 * t231 + t314 * t229 + t311 * t201 + t312 * t199 + m(7) * (t26 * t27 + t30 * t46 + t31 * t47) + m(6) * (t45 * t58 + t50 * t74 + t51 * t75); (t13 + t14) * t233 + (t4 + t3) * t201 + (t1 + t2) * t199 + m(7) * (t27 ^ 2 + t46 ^ 2 + t47 ^ 2) + m(6) * (t58 ^ 2 + t74 ^ 2 + t75 ^ 2); m(7) * (t161 * t57 + t163 * t56); m(7) * (t161 * t28 + t163 * t29 + t189 * t24); m(7) * (t189 * t270 + (-t161 * t274 + t163 * t273) * t269); m(7) * (t161 * t31 + t163 * t30 + t189 * t26); m(7) * (t161 * t47 + t163 * t46 + t189 * t27); m(7) * (t161 ^ 2 + t163 ^ 2 + t189 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t105(1) t105(2) t105(4) t105(7) t105(11) t105(16); t105(2) t105(3) t105(5) t105(8) t105(12) t105(17); t105(4) t105(5) t105(6) t105(9) t105(13) t105(18); t105(7) t105(8) t105(9) t105(10) t105(14) t105(19); t105(11) t105(12) t105(13) t105(14) t105(15) t105(20); t105(16) t105(17) t105(18) t105(19) t105(20) t105(21);];
Mq  = res;
