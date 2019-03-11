% Calculate joint inertia matrix for
% S6RRPRRP5
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:52
% EndTime: 2019-03-09 11:58:10
% DurationCPUTime: 7.11s
% Computational Cost: add. (28349->628), mult. (72140->876), div. (0->0), fcn. (94794->12), ass. (0->284)
t352 = rSges(7,3) + qJ(6) + pkin(10);
t270 = cos(pkin(6));
t335 = sin(pkin(11));
t336 = cos(pkin(11));
t341 = sin(qJ(2));
t343 = cos(qJ(2));
t278 = t335 * t343 + t336 * t341;
t245 = t278 * t270;
t255 = -t341 * t335 + t343 * t336;
t274 = sin(qJ(1));
t276 = cos(qJ(1));
t230 = t245 * t276 + t274 * t255;
t273 = sin(qJ(4));
t269 = sin(pkin(6));
t329 = t269 * t276;
t342 = cos(qJ(4));
t201 = t230 * t342 - t273 * t329;
t277 = t270 * t255;
t229 = -t274 * t278 + t276 * t277;
t272 = sin(qJ(5));
t275 = cos(qJ(5));
t160 = -t201 * t272 - t229 * t275;
t334 = t229 * t272;
t161 = t201 * t275 - t334;
t301 = t269 * t342;
t200 = t230 * t273 + t276 * t301;
t353 = t161 * rSges(7,1) + t160 * rSges(7,2) - pkin(5) * t334 + t200 * t352;
t243 = t255 * t269;
t244 = t278 * t269;
t207 = Icges(4,5) * t244 + Icges(4,6) * t243 + Icges(4,3) * t270;
t208 = Icges(4,4) * t244 + Icges(4,2) * t243 + Icges(4,6) * t270;
t209 = Icges(4,1) * t244 + Icges(4,4) * t243 + Icges(4,5) * t270;
t239 = Icges(3,3) * t270 + (Icges(3,5) * t341 + Icges(3,6) * t343) * t269;
t240 = Icges(3,6) * t270 + (Icges(3,4) * t341 + Icges(3,2) * t343) * t269;
t241 = Icges(3,5) * t270 + (Icges(3,1) * t341 + Icges(3,4) * t343) * t269;
t300 = t269 * t341;
t351 = t269 * t343 * t240 + t243 * t208 + t244 * t209 + t241 * t300 + (t207 + t239) * t270;
t232 = -t274 * t245 + t255 * t276;
t330 = t269 * t274;
t203 = t232 * t342 + t273 * t330;
t231 = -t274 * t277 - t276 * t278;
t162 = -t203 * t272 - t231 * t275;
t333 = t231 * t272;
t163 = t203 * t275 - t333;
t202 = t232 * t273 - t274 * t301;
t266 = pkin(5) * t275 + pkin(4);
t350 = t163 * rSges(7,1) + t162 * rSges(7,2) - pkin(5) * t333 + t352 * t202 + t203 * t266;
t348 = t269 ^ 2;
t347 = m(4) / 0.2e1;
t346 = m(5) / 0.2e1;
t345 = m(6) / 0.2e1;
t344 = m(7) / 0.2e1;
t340 = pkin(1) * t276;
t339 = -pkin(4) + t266;
t196 = t200 * pkin(10);
t338 = t201 * t339 - t196 + t353;
t151 = t203 * pkin(4) + pkin(10) * t202;
t337 = -t151 + t350;
t332 = t243 * t272;
t246 = t270 * t341 * pkin(2) + (-pkin(8) - qJ(3)) * t269;
t331 = t246 * t276;
t267 = pkin(2) * t343 + pkin(1);
t328 = t274 * t267;
t108 = t161 * rSges(6,1) + t160 * rSges(6,2) + t200 * rSges(6,3);
t150 = t201 * pkin(4) + t196;
t327 = -t108 - t150;
t110 = t163 * rSges(6,1) + t162 * rSges(6,2) + t202 * rSges(6,3);
t326 = t110 + t151;
t325 = t351 * t270;
t234 = t244 * t342 + t270 * t273;
t190 = -t234 * t272 - t243 * t275;
t191 = t234 * t275 - t332;
t233 = t244 * t273 - t270 * t342;
t324 = rSges(7,1) * t191 + rSges(7,2) * t190 - pkin(5) * t332 + t339 * t234 + (-pkin(10) + t352) * t233;
t181 = t232 * pkin(3) - pkin(9) * t231;
t260 = t276 * t267;
t223 = -t340 + t260 + (-t269 * pkin(8) - t246) * t274;
t211 = t270 * t223;
t323 = t270 * t181 + t211;
t180 = t230 * pkin(3) - t229 * pkin(9);
t264 = pkin(8) * t329;
t222 = t331 + t264 + (-pkin(1) + t267) * t274;
t322 = -t180 - t222;
t321 = t222 * t330 + t223 * t329;
t165 = Icges(4,5) * t230 + Icges(4,6) * t229 - Icges(4,3) * t329;
t297 = t276 * t343;
t298 = t274 * t341;
t250 = t270 * t297 - t298;
t296 = t276 * t341;
t299 = t274 * t343;
t251 = t270 * t296 + t299;
t213 = Icges(3,5) * t251 + Icges(3,6) * t250 - Icges(3,3) * t329;
t320 = -t213 - t165;
t166 = Icges(4,5) * t232 + Icges(4,6) * t231 + Icges(4,3) * t330;
t252 = -t270 * t299 - t296;
t253 = -t270 * t298 + t297;
t214 = Icges(3,5) * t253 + Icges(3,6) * t252 + Icges(3,3) * t330;
t319 = t214 + t166;
t256 = pkin(2) * t300 + t270 * qJ(3);
t318 = -pkin(3) * t244 + pkin(9) * t243 - t256;
t103 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t200;
t95 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t200;
t99 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t200;
t33 = t103 * t161 + t160 * t99 + t200 * t95;
t100 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t202;
t104 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t202;
t96 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t202;
t34 = t100 * t160 + t104 * t161 + t200 * t96;
t132 = Icges(7,5) * t191 + Icges(7,6) * t190 + Icges(7,3) * t233;
t134 = Icges(7,4) * t191 + Icges(7,2) * t190 + Icges(7,6) * t233;
t136 = Icges(7,1) * t191 + Icges(7,4) * t190 + Icges(7,5) * t233;
t52 = t132 * t200 + t134 * t160 + t136 * t161;
t1 = t200 * t33 + t202 * t34 + t233 * t52;
t101 = Icges(6,4) * t161 + Icges(6,2) * t160 + Icges(6,6) * t200;
t105 = Icges(6,1) * t161 + Icges(6,4) * t160 + Icges(6,5) * t200;
t97 = Icges(6,5) * t161 + Icges(6,6) * t160 + Icges(6,3) * t200;
t35 = t101 * t160 + t105 * t161 + t200 * t97;
t102 = Icges(6,4) * t163 + Icges(6,2) * t162 + Icges(6,6) * t202;
t106 = Icges(6,1) * t163 + Icges(6,4) * t162 + Icges(6,5) * t202;
t98 = Icges(6,5) * t163 + Icges(6,6) * t162 + Icges(6,3) * t202;
t36 = t102 * t160 + t106 * t161 + t200 * t98;
t133 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t233;
t135 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t233;
t137 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t233;
t53 = t133 * t200 + t135 * t160 + t137 * t161;
t2 = t200 * t35 + t202 * t36 + t233 * t53;
t317 = -t1 / 0.2e1 - t2 / 0.2e1;
t37 = t103 * t163 + t162 * t99 + t202 * t95;
t38 = t100 * t162 + t104 * t163 + t202 * t96;
t54 = t132 * t202 + t134 * t162 + t136 * t163;
t3 = t200 * t37 + t202 * t38 + t233 * t54;
t39 = t101 * t162 + t105 * t163 + t202 * t97;
t40 = t102 * t162 + t106 * t163 + t202 * t98;
t55 = t133 * t202 + t135 * t162 + t137 * t163;
t4 = t200 * t39 + t202 * t40 + t233 * t55;
t316 = t3 / 0.2e1 + t4 / 0.2e1;
t5 = -t229 * t33 - t231 * t34 - t243 * t52;
t6 = -t229 * t35 - t231 * t36 - t243 * t53;
t315 = t6 / 0.2e1 + t5 / 0.2e1;
t7 = -t229 * t37 - t231 * t38 - t243 * t54;
t8 = -t229 * t39 - t231 * t40 - t243 * t55;
t314 = t7 / 0.2e1 + t8 / 0.2e1;
t10 = t53 * t270 + (t274 * t36 - t276 * t35) * t269;
t9 = t52 * t270 + (t274 * t34 - t276 * t33) * t269;
t312 = t10 / 0.2e1 + t9 / 0.2e1;
t11 = t54 * t270 + (t274 * t38 - t276 * t37) * t269;
t12 = t55 * t270 + (t274 * t40 - t276 * t39) * t269;
t311 = t11 / 0.2e1 + t12 / 0.2e1;
t43 = t103 * t191 + t190 * t99 + t233 * t95;
t44 = t100 * t190 + t104 * t191 + t233 * t96;
t68 = t233 * t132 + t190 * t134 + t191 * t136;
t57 = t68 * t233;
t13 = t43 * t200 + t44 * t202 + t57;
t45 = t101 * t190 + t105 * t191 + t233 * t97;
t46 = t102 * t190 + t106 * t191 + t233 * t98;
t69 = t233 * t133 + t190 * t135 + t191 * t137;
t58 = t69 * t233;
t14 = t45 * t200 + t46 * t202 + t58;
t310 = t14 / 0.2e1 + t13 / 0.2e1;
t59 = t68 * t243;
t15 = -t43 * t229 - t44 * t231 - t59;
t60 = t69 * t243;
t16 = -t45 * t229 - t46 * t231 - t60;
t309 = t16 / 0.2e1 + t15 / 0.2e1;
t65 = t68 * t270;
t17 = t65 + (t44 * t274 - t43 * t276) * t269;
t66 = t69 * t270;
t18 = t66 + (t46 * t274 - t45 * t276) * t269;
t308 = t18 / 0.2e1 + t17 / 0.2e1;
t176 = Icges(5,5) * t234 - Icges(5,6) * t233 - Icges(5,3) * t243;
t177 = Icges(5,4) * t234 - Icges(5,2) * t233 - Icges(5,6) * t243;
t178 = Icges(5,1) * t234 - Icges(5,4) * t233 - Icges(5,5) * t243;
t85 = -t243 * t176 - t233 * t177 + t234 * t178;
t307 = t270 * t151 + t323;
t306 = -t150 + t322;
t187 = pkin(4) * t234 + pkin(10) * t233;
t304 = -t187 + t318;
t130 = t203 * rSges(5,1) - t202 * rSges(5,2) - t231 * rSges(5,3);
t172 = t232 * rSges(4,1) + t231 * rSges(4,2) + rSges(4,3) * t330;
t220 = t253 * rSges(3,1) + t252 * rSges(3,2) + rSges(3,3) * t330;
t295 = t269 * (-rSges(4,1) * t244 - rSges(4,2) * t243 - rSges(4,3) * t270 - t256);
t294 = -t274 * t246 + t260;
t293 = t180 * t330 + t181 * t329 + t321;
t179 = rSges(5,1) * t234 - rSges(5,2) * t233 - rSges(5,3) * t243;
t292 = t269 * (-t179 + t318);
t289 = -t230 * rSges(4,1) - t229 * rSges(4,2);
t139 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t233;
t287 = t269 * (-t139 + t304);
t286 = t44 / 0.2e1 + t54 / 0.2e1 + t55 / 0.2e1 + t46 / 0.2e1;
t285 = t52 / 0.2e1 + t53 / 0.2e1 + t45 / 0.2e1 + t43 / 0.2e1;
t284 = t150 * t330 + t151 * t329 + t293;
t283 = t269 * (t304 - t324);
t282 = t181 + t294;
t129 = t201 * rSges(5,1) - t200 * rSges(5,2) - t229 * rSges(5,3);
t219 = t251 * rSges(3,1) + t250 * rSges(3,2) - rSges(3,3) * t329;
t281 = -t180 - t328 - t331;
t124 = Icges(5,5) * t203 - Icges(5,6) * t202 - Icges(5,3) * t231;
t126 = Icges(5,4) * t203 - Icges(5,2) * t202 - Icges(5,6) * t231;
t128 = Icges(5,1) * t203 - Icges(5,4) * t202 - Icges(5,5) * t231;
t71 = -t124 * t243 - t126 * t233 + t128 * t234;
t81 = -t176 * t231 - t177 * t202 + t178 * t203;
t280 = -t81 / 0.2e1 - t71 / 0.2e1 - t286;
t123 = Icges(5,5) * t201 - Icges(5,6) * t200 - Icges(5,3) * t229;
t125 = Icges(5,4) * t201 - Icges(5,2) * t200 - Icges(5,6) * t229;
t127 = Icges(5,1) * t201 - Icges(5,4) * t200 - Icges(5,5) * t229;
t70 = -t123 * t243 - t125 * t233 + t127 * t234;
t80 = -t176 * t229 - t177 * t200 + t178 * t201;
t279 = -t80 / 0.2e1 - t70 / 0.2e1 - t285;
t258 = rSges(2,1) * t276 - t274 * rSges(2,2);
t257 = -t274 * rSges(2,1) - rSges(2,2) * t276;
t242 = t270 * rSges(3,3) + (rSges(3,1) * t341 + rSges(3,2) * t343) * t269;
t218 = Icges(3,1) * t253 + Icges(3,4) * t252 + Icges(3,5) * t330;
t217 = Icges(3,1) * t251 + Icges(3,4) * t250 - Icges(3,5) * t329;
t216 = Icges(3,4) * t253 + Icges(3,2) * t252 + Icges(3,6) * t330;
t215 = Icges(3,4) * t251 + Icges(3,2) * t250 - Icges(3,6) * t329;
t199 = pkin(8) * t330 + t220 + t340;
t198 = -t274 * pkin(1) - t219 + t264;
t184 = -t270 * t219 - t242 * t329;
t183 = t220 * t270 - t242 * t330;
t171 = -rSges(4,3) * t329 - t289;
t170 = Icges(4,1) * t232 + Icges(4,4) * t231 + Icges(4,5) * t330;
t169 = Icges(4,1) * t230 + Icges(4,4) * t229 - Icges(4,5) * t329;
t168 = Icges(4,4) * t232 + Icges(4,2) * t231 + Icges(4,6) * t330;
t167 = Icges(4,4) * t230 + Icges(4,2) * t229 - Icges(4,6) * t329;
t164 = (t219 * t274 + t220 * t276) * t269;
t159 = t239 * t330 + t240 * t252 + t241 * t253;
t158 = -t239 * t329 + t250 * t240 + t251 * t241;
t153 = t229 * t187;
t144 = t294 + t172;
t143 = -t328 + (rSges(4,3) * t269 - t246) * t276 + t289;
t142 = t243 * t151;
t141 = t270 * t214 + (t216 * t343 + t218 * t341) * t269;
t140 = t270 * t213 + (t215 * t343 + t217 * t341) * t269;
t131 = t231 * t150;
t114 = (-t171 - t222) * t270 + t276 * t295;
t113 = t172 * t270 + t274 * t295 + t211;
t112 = t207 * t330 + t208 * t231 + t209 * t232;
t111 = -t207 * t329 + t229 * t208 + t230 * t209;
t94 = t282 + t130;
t93 = -t129 + t281;
t90 = (t171 * t274 + t172 * t276) * t269 + t321;
t89 = -t130 * t243 + t179 * t231;
t88 = t129 * t243 - t179 * t229;
t87 = t166 * t270 + t168 * t243 + t170 * t244;
t86 = t165 * t270 + t167 * t243 + t169 * t244;
t84 = t85 * t270;
t83 = t85 * t243;
t82 = -t129 * t231 + t130 * t229;
t79 = (-t129 + t322) * t270 + t276 * t292;
t78 = t130 * t270 + t274 * t292 + t323;
t77 = t282 + t326;
t76 = t281 + t327;
t75 = t110 * t233 - t139 * t202;
t74 = -t108 * t233 + t139 * t200;
t73 = t282 + t350;
t72 = -t201 * t266 + t281 - t353;
t67 = (t129 * t274 + t130 * t276) * t269 + t293;
t64 = -t124 * t231 - t126 * t202 + t128 * t203;
t63 = -t123 * t231 - t125 * t202 + t127 * t203;
t62 = -t124 * t229 - t126 * t200 + t128 * t201;
t61 = -t123 * t229 - t125 * t200 + t127 * t201;
t56 = t108 * t202 - t110 * t200;
t51 = -t110 * t243 - t142 + (t139 + t187) * t231;
t50 = -t139 * t229 - t243 * t327 - t153;
t49 = (-t108 + t306) * t270 + t276 * t287;
t48 = t110 * t270 + t274 * t287 + t307;
t47 = -t108 * t231 + t229 * t326 - t131;
t42 = -t202 * t324 + t233 * t337;
t41 = t200 * t324 - t233 * t338;
t32 = (t108 * t274 + t110 * t276) * t269 + t284;
t31 = -t142 - t337 * t243 + (t187 + t324) * t231;
t30 = -t153 - t324 * t229 - (-t150 - t338) * t243;
t29 = (t306 - t338) * t270 + t276 * t283;
t28 = t270 * t337 + t274 * t283 + t307;
t27 = -t200 * t337 + t202 * t338;
t26 = t84 + (t71 * t274 - t70 * t276) * t269;
t25 = -t70 * t229 - t71 * t231 - t83;
t24 = -t131 - t338 * t231 + (t151 + t337) * t229;
t23 = (t274 * t338 + t276 * t337) * t269 + t284;
t22 = t81 * t270 + (t274 * t64 - t276 * t63) * t269;
t21 = t80 * t270 + (t274 * t62 - t276 * t61) * t269;
t20 = -t229 * t63 - t231 * t64 - t243 * t81;
t19 = -t229 * t61 - t231 * t62 - t243 * t80;
t91 = [Icges(2,3) + m(6) * (t76 ^ 2 + t77 ^ 2) + m(7) * (t72 ^ 2 + t73 ^ 2) + m(5) * (t93 ^ 2 + t94 ^ 2) + m(4) * (t143 ^ 2 + t144 ^ 2) + m(3) * (t198 ^ 2 + t199 ^ 2) + m(2) * (t257 ^ 2 + t258 ^ 2) + t85 + t69 + t68 + t351; t65 + t84 + t66 + m(7) * (t28 * t73 + t29 * t72) + m(6) * (t48 * t77 + t49 * t76) + m(5) * (t78 * t94 + t79 * t93) + m(4) * (t113 * t144 + t114 * t143) + m(3) * (t183 * t199 + t184 * t198) + ((-t140 / 0.2e1 - t86 / 0.2e1 - t111 / 0.2e1 - t158 / 0.2e1 + t279) * t276 + (t141 / 0.2e1 + t87 / 0.2e1 + t112 / 0.2e1 + t159 / 0.2e1 - t280) * t274) * t269 + t325; (t17 + t18 + t26 + t325) * t270 + m(7) * (t23 ^ 2 + t28 ^ 2 + t29 ^ 2) + m(6) * (t32 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t67 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(4) * (t113 ^ 2 + t114 ^ 2 + t90 ^ 2) + m(3) * (t164 ^ 2 + t183 ^ 2 + t184 ^ 2) + ((-t10 - t21 - t9 + ((t229 * t167 + t230 * t169 + t250 * t215 + t251 * t217) * t269 + t320 * t348 * t276) * t276 + (-t111 - t140 - t158 - t86) * t270) * t276 + (t11 + t12 + t22 + ((t168 * t231 + t170 * t232 + t216 * t252 + t218 * t253) * t269 + t319 * t348 * t274) * t274 + (t159 + t112 + t141 + t87) * t270 + (-t167 * t231 - t229 * t168 - t169 * t232 - t230 * t170 - t215 * t252 - t250 * t216 - t217 * t253 - t251 * t218 + (t274 * t320 + t319 * t276) * t269) * t329) * t274) * t269; 0.2e1 * ((t274 * t76 - t276 * t77) * t345 + (t274 * t72 - t276 * t73) * t344 + (t274 * t93 - t276 * t94) * t346 + (t143 * t274 - t144 * t276) * t347) * t269; m(7) * (t270 * t23 + (t274 * t29 - t276 * t28) * t269) + m(6) * (t270 * t32 + (t274 * t49 - t276 * t48) * t269) + m(5) * (t270 * t67 + (t274 * t79 - t276 * t78) * t269) + m(4) * (t270 * t90 + (-t113 * t276 + t114 * t274) * t269); 0.2e1 * (t347 + t346 + t345 + t344) * (t270 ^ 2 + (t274 ^ 2 + t276 ^ 2) * t348); -t59 - t60 - t83 + m(6) * (t50 * t76 + t51 * t77) + m(7) * (t30 * t72 + t31 * t73) + m(5) * (t88 * t93 + t89 * t94) + t280 * t231 + t279 * t229; (t25 / 0.2e1 + t309) * t270 - (t26 / 0.2e1 + t308) * t243 + (-t22 / 0.2e1 - t311) * t231 + (-t21 / 0.2e1 - t312) * t229 + m(7) * (t23 * t24 + t28 * t31 + t29 * t30) + m(6) * (t32 * t47 + t48 * t51 + t49 * t50) + m(5) * (t67 * t82 + t78 * t89 + t79 * t88) + ((-t19 / 0.2e1 - t315) * t276 + (t20 / 0.2e1 + t314) * t274) * t269; m(5) * (t82 * t270 + (t274 * t88 - t276 * t89) * t269) + m(6) * (t47 * t270 + (t274 * t50 - t276 * t51) * t269) + m(7) * (t24 * t270 + (t274 * t30 - t276 * t31) * t269); -(t16 + t15 + t25) * t243 + (-t7 - t8 - t20) * t231 + (-t5 - t6 - t19) * t229 + m(7) * (t24 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t47 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t82 ^ 2 + t88 ^ 2 + t89 ^ 2); t57 + t58 + m(6) * (t74 * t76 + t75 * t77) + m(7) * (t41 * t72 + t42 * t73) + t286 * t202 + t285 * t200; t310 * t270 + t308 * t233 + t311 * t202 + t312 * t200 + m(7) * (t23 * t27 + t28 * t42 + t29 * t41) + m(6) * (t32 * t56 + t48 * t75 + t49 * t74) + (t274 * t316 + t276 * t317) * t269; m(6) * (t56 * t270 + (t274 * t74 - t276 * t75) * t269) + m(7) * (t27 * t270 + (t274 * t41 - t276 * t42) * t269); -t310 * t243 + t309 * t233 - t316 * t231 + t317 * t229 + t314 * t202 + t315 * t200 + m(7) * (t24 * t27 + t30 * t41 + t31 * t42) + m(6) * (t47 * t56 + t50 * t74 + t51 * t75); (t13 + t14) * t233 + (t3 + t4) * t202 + (t1 + t2) * t200 + m(7) * (t27 ^ 2 + t41 ^ 2 + t42 ^ 2) + m(6) * (t56 ^ 2 + t74 ^ 2 + t75 ^ 2); m(7) * (t200 * t73 + t202 * t72); m(7) * (t200 * t28 + t202 * t29 + t23 * t233); m(7) * (t233 * t270 + (-t200 * t276 + t202 * t274) * t269); m(7) * (t200 * t31 + t202 * t30 + t233 * t24); m(7) * (t200 * t42 + t202 * t41 + t233 * t27); m(7) * (t200 ^ 2 + t202 ^ 2 + t233 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t91(1) t91(2) t91(4) t91(7) t91(11) t91(16); t91(2) t91(3) t91(5) t91(8) t91(12) t91(17); t91(4) t91(5) t91(6) t91(9) t91(13) t91(18); t91(7) t91(8) t91(9) t91(10) t91(14) t91(19); t91(11) t91(12) t91(13) t91(14) t91(15) t91(20); t91(16) t91(17) t91(18) t91(19) t91(20) t91(21);];
Mq  = res;
