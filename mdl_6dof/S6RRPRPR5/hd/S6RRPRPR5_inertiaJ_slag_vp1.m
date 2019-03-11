% Calculate joint inertia matrix for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRPR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:29:03
% EndTime: 2019-03-09 10:29:17
% DurationCPUTime: 6.94s
% Computational Cost: add. (25011->619), mult. (58485->870), div. (0->0), fcn. (76581->14), ass. (0->288)
t331 = sin(pkin(11));
t332 = cos(pkin(11));
t337 = sin(qJ(2));
t339 = cos(qJ(2));
t256 = -t337 * t331 + t339 * t332;
t274 = sin(pkin(6));
t244 = t256 * t274;
t282 = t331 * t339 + t332 * t337;
t245 = t282 * t274;
t276 = cos(pkin(6));
t208 = Icges(4,5) * t245 + Icges(4,6) * t244 + Icges(4,3) * t276;
t209 = Icges(4,4) * t245 + Icges(4,2) * t244 + Icges(4,6) * t276;
t210 = Icges(4,1) * t245 + Icges(4,4) * t244 + Icges(4,5) * t276;
t240 = Icges(3,3) * t276 + (Icges(3,5) * t337 + Icges(3,6) * t339) * t274;
t241 = Icges(3,6) * t276 + (Icges(3,4) * t337 + Icges(3,2) * t339) * t274;
t242 = Icges(3,5) * t276 + (Icges(3,1) * t337 + Icges(3,4) * t339) * t274;
t302 = t274 * t337;
t352 = t274 * t339 * t241 + t244 * t209 + t245 * t210 + t242 * t302 + (t208 + t240) * t276;
t345 = m(7) / 0.2e1;
t346 = m(6) / 0.2e1;
t313 = t346 + t345;
t351 = 0.2e1 * t313;
t246 = t282 * t276;
t279 = sin(qJ(1));
t280 = cos(qJ(1));
t233 = -t279 * t246 + t256 * t280;
t278 = sin(qJ(4));
t338 = cos(qJ(4));
t303 = t274 * t338;
t203 = t233 * t278 - t279 * t303;
t327 = t274 * t279;
t204 = t233 * t338 + t278 * t327;
t275 = cos(pkin(12));
t267 = pkin(5) * t275 + pkin(4);
t277 = -pkin(10) - qJ(5);
t281 = t276 * t256;
t232 = -t279 * t281 - t280 * t282;
t273 = sin(pkin(12));
t329 = t232 * t273;
t272 = pkin(12) + qJ(6);
t269 = sin(t272);
t270 = cos(t272);
t152 = -t204 * t269 - t232 * t270;
t153 = t204 * t270 - t232 * t269;
t93 = t153 * rSges(7,1) + t152 * rSges(7,2) + t203 * rSges(7,3);
t350 = -pkin(5) * t329 - t203 * t277 + t204 * t267 + t93;
t349 = t274 ^ 2;
t348 = m(4) / 0.2e1;
t347 = m(5) / 0.2e1;
t231 = t246 * t280 + t279 * t256;
t201 = t231 * t278 + t280 * t303;
t234 = t245 * t278 - t276 * t338;
t326 = t274 * t280;
t202 = t231 * t338 - t278 * t326;
t230 = -t279 * t282 + t280 * t281;
t150 = -t202 * t269 - t230 * t270;
t151 = t202 * t270 - t230 * t269;
t86 = Icges(7,5) * t151 + Icges(7,6) * t150 + Icges(7,3) * t201;
t88 = Icges(7,4) * t151 + Icges(7,2) * t150 + Icges(7,6) * t201;
t90 = Icges(7,1) * t151 + Icges(7,4) * t150 + Icges(7,5) * t201;
t28 = t150 * t88 + t151 * t90 + t201 * t86;
t87 = Icges(7,5) * t153 + Icges(7,6) * t152 + Icges(7,3) * t203;
t89 = Icges(7,4) * t153 + Icges(7,2) * t152 + Icges(7,6) * t203;
t91 = Icges(7,1) * t153 + Icges(7,4) * t152 + Icges(7,5) * t203;
t29 = t150 * t89 + t151 * t91 + t201 * t87;
t235 = t245 * t338 + t276 * t278;
t187 = -t235 * t269 - t244 * t270;
t188 = t235 * t270 - t244 * t269;
t115 = Icges(7,5) * t188 + Icges(7,6) * t187 + Icges(7,3) * t234;
t116 = Icges(7,4) * t188 + Icges(7,2) * t187 + Icges(7,6) * t234;
t117 = Icges(7,1) * t188 + Icges(7,4) * t187 + Icges(7,5) * t234;
t44 = t115 * t201 + t116 * t150 + t117 * t151;
t1 = t201 * t28 + t203 * t29 + t234 * t44;
t344 = -t1 / 0.2e1;
t37 = t187 * t88 + t188 * t90 + t234 * t86;
t38 = t187 * t89 + t188 * t91 + t234 * t87;
t54 = t234 * t115 + t187 * t116 + t188 * t117;
t50 = t54 * t234;
t11 = t37 * t201 + t38 * t203 + t50;
t343 = t11 / 0.2e1;
t342 = t201 / 0.2e1;
t341 = t203 / 0.2e1;
t340 = t234 / 0.2e1;
t336 = t280 * pkin(1);
t335 = -pkin(4) + t267;
t193 = t201 * qJ(5);
t330 = t230 * t273;
t312 = pkin(5) * t330;
t290 = -t151 * rSges(7,1) - t150 * rSges(7,2);
t92 = t201 * rSges(7,3) - t290;
t334 = -t201 * t277 + t202 * t335 - t193 - t312 + t92;
t146 = t204 * pkin(4) + qJ(5) * t203;
t333 = -t146 + t350;
t328 = t244 * t273;
t268 = pkin(2) * t339 + pkin(1);
t325 = t279 * t268;
t247 = t276 * t337 * pkin(2) + (-pkin(8) - qJ(3)) * t274;
t324 = t280 * t247;
t159 = -t202 * t273 - t230 * t275;
t160 = t202 * t275 - t330;
t102 = t160 * rSges(6,1) + t159 * rSges(6,2) + t201 * rSges(6,3);
t145 = t202 * pkin(4) + t193;
t323 = -t102 - t145;
t161 = -t204 * t273 - t232 * t275;
t162 = t204 * t275 - t329;
t103 = t162 * rSges(6,1) + t161 * rSges(6,2) + t203 * rSges(6,3);
t322 = t103 + t146;
t321 = t352 * t276;
t118 = rSges(7,1) * t188 + rSges(7,2) * t187 + rSges(7,3) * t234;
t320 = -t118 + pkin(5) * t328 - t335 * t235 - (-qJ(5) - t277) * t234;
t180 = t233 * pkin(3) - pkin(9) * t232;
t261 = t280 * t268;
t224 = -t336 + t261 + (-t274 * pkin(8) - t247) * t279;
t212 = t276 * t224;
t319 = t276 * t180 + t212;
t179 = t231 * pkin(3) - t230 * pkin(9);
t265 = pkin(8) * t326;
t223 = t324 + t265 + (-pkin(1) + t268) * t279;
t318 = -t179 - t223;
t317 = t223 * t327 + t224 * t326;
t164 = Icges(4,5) * t231 + Icges(4,6) * t230 - Icges(4,3) * t326;
t299 = t280 * t339;
t300 = t279 * t337;
t251 = t276 * t299 - t300;
t298 = t280 * t337;
t301 = t279 * t339;
t252 = t276 * t298 + t301;
t214 = Icges(3,5) * t252 + Icges(3,6) * t251 - Icges(3,3) * t326;
t316 = -t214 - t164;
t165 = Icges(4,5) * t233 + Icges(4,6) * t232 + Icges(4,3) * t327;
t253 = -t276 * t301 - t298;
t254 = -t276 * t300 + t299;
t215 = Icges(3,5) * t254 + Icges(3,6) * t253 + Icges(3,3) * t327;
t315 = t215 + t165;
t257 = pkin(2) * t302 + t276 * qJ(3);
t314 = -pkin(3) * t245 + pkin(9) * t244 - t257;
t311 = t37 / 0.2e1 + t44 / 0.2e1;
t45 = t115 * t203 + t116 * t152 + t117 * t153;
t310 = t38 / 0.2e1 + t45 / 0.2e1;
t190 = -t235 * t273 - t244 * t275;
t191 = t235 * t275 - t328;
t130 = Icges(6,5) * t191 + Icges(6,6) * t190 + Icges(6,3) * t234;
t131 = Icges(6,4) * t191 + Icges(6,2) * t190 + Icges(6,6) * t234;
t132 = Icges(6,1) * t191 + Icges(6,4) * t190 + Icges(6,5) * t234;
t61 = t234 * t130 + t190 * t131 + t191 * t132;
t175 = Icges(5,5) * t235 - Icges(5,6) * t234 - Icges(5,3) * t244;
t176 = Icges(5,4) * t235 - Icges(5,2) * t234 - Icges(5,6) * t244;
t177 = Icges(5,1) * t235 - Icges(5,4) * t234 - Icges(5,5) * t244;
t78 = -t244 * t175 - t234 * t176 + t235 * t177;
t309 = t276 * t146 + t319;
t308 = -t145 + t318;
t186 = t235 * pkin(4) + t234 * qJ(5);
t306 = -t186 + t314;
t128 = t204 * rSges(5,1) - t203 * rSges(5,2) - t232 * rSges(5,3);
t171 = t233 * rSges(4,1) + t232 * rSges(4,2) + rSges(4,3) * t327;
t221 = t254 * rSges(3,1) + t253 * rSges(3,2) + rSges(3,3) * t327;
t297 = t274 * (-rSges(4,1) * t245 - rSges(4,2) * t244 - rSges(4,3) * t276 - t257);
t296 = -t279 * t247 + t261;
t295 = t179 * t327 + t180 * t326 + t317;
t178 = rSges(5,1) * t235 - rSges(5,2) * t234 - rSges(5,3) * t244;
t294 = t274 * (-t178 + t314);
t291 = -t231 * rSges(4,1) - t230 * rSges(4,2);
t133 = rSges(6,1) * t191 + rSges(6,2) * t190 + rSges(6,3) * t234;
t289 = t274 * (-t133 + t306);
t288 = t145 * t327 + t146 * t326 + t295;
t287 = t274 * (t306 + t320);
t286 = t180 + t296;
t127 = t202 * rSges(5,1) - t201 * rSges(5,2) - t230 * rSges(5,3);
t220 = t252 * rSges(3,1) + t251 * rSges(3,2) - rSges(3,3) * t326;
t285 = -t179 - t324 - t325;
t100 = Icges(6,1) * t160 + Icges(6,4) * t159 + Icges(6,5) * t201;
t96 = Icges(6,5) * t160 + Icges(6,6) * t159 + Icges(6,3) * t201;
t98 = Icges(6,4) * t160 + Icges(6,2) * t159 + Icges(6,6) * t201;
t39 = t100 * t191 + t190 * t98 + t234 * t96;
t46 = t130 * t201 + t131 * t159 + t132 * t160;
t121 = Icges(5,5) * t202 - Icges(5,6) * t201 - Icges(5,3) * t230;
t123 = Icges(5,4) * t202 - Icges(5,2) * t201 - Icges(5,6) * t230;
t125 = Icges(5,1) * t202 - Icges(5,4) * t201 - Icges(5,5) * t230;
t65 = -t121 * t244 - t123 * t234 + t125 * t235;
t73 = -t175 * t230 - t176 * t201 + t177 * t202;
t284 = -t73 / 0.2e1 - t46 / 0.2e1 - t65 / 0.2e1 - t39 / 0.2e1 - t311;
t101 = Icges(6,1) * t162 + Icges(6,4) * t161 + Icges(6,5) * t203;
t97 = Icges(6,5) * t162 + Icges(6,6) * t161 + Icges(6,3) * t203;
t99 = Icges(6,4) * t162 + Icges(6,2) * t161 + Icges(6,6) * t203;
t40 = t101 * t191 + t190 * t99 + t234 * t97;
t47 = t130 * t203 + t131 * t161 + t132 * t162;
t122 = Icges(5,5) * t204 - Icges(5,6) * t203 - Icges(5,3) * t232;
t124 = Icges(5,4) * t204 - Icges(5,2) * t203 - Icges(5,6) * t232;
t126 = Icges(5,1) * t204 - Icges(5,4) * t203 - Icges(5,5) * t232;
t66 = -t122 * t244 - t124 * t234 + t126 * t235;
t74 = -t175 * t232 - t176 * t203 + t177 * t204;
t283 = -t66 / 0.2e1 - t40 / 0.2e1 - t74 / 0.2e1 - t47 / 0.2e1 - t310;
t259 = rSges(2,1) * t280 - t279 * rSges(2,2);
t258 = -t279 * rSges(2,1) - rSges(2,2) * t280;
t243 = t276 * rSges(3,3) + (rSges(3,1) * t337 + rSges(3,2) * t339) * t274;
t219 = Icges(3,1) * t254 + Icges(3,4) * t253 + Icges(3,5) * t327;
t218 = Icges(3,1) * t252 + Icges(3,4) * t251 - Icges(3,5) * t326;
t217 = Icges(3,4) * t254 + Icges(3,2) * t253 + Icges(3,6) * t327;
t216 = Icges(3,4) * t252 + Icges(3,2) * t251 - Icges(3,6) * t326;
t200 = pkin(8) * t327 + t221 + t336;
t199 = -t279 * pkin(1) - t220 + t265;
t183 = -t276 * t220 - t243 * t326;
t182 = t221 * t276 - t243 * t327;
t170 = -rSges(4,3) * t326 - t291;
t169 = Icges(4,1) * t233 + Icges(4,4) * t232 + Icges(4,5) * t327;
t168 = Icges(4,1) * t231 + Icges(4,4) * t230 - Icges(4,5) * t326;
t167 = Icges(4,4) * t233 + Icges(4,2) * t232 + Icges(4,6) * t327;
t166 = Icges(4,4) * t231 + Icges(4,2) * t230 - Icges(4,6) * t326;
t163 = (t220 * t279 + t221 * t280) * t274;
t158 = t240 * t327 + t241 * t253 + t242 * t254;
t157 = -t240 * t326 + t251 * t241 + t252 * t242;
t154 = t230 * t186;
t139 = t296 + t171;
t138 = -t325 + (rSges(4,3) * t274 - t247) * t280 + t291;
t137 = t244 * t146;
t135 = t276 * t215 + (t217 * t339 + t219 * t337) * t274;
t134 = t276 * t214 + (t216 * t339 + t218 * t337) * t274;
t129 = t232 * t145;
t109 = (-t170 - t223) * t276 + t280 * t297;
t108 = t276 * t171 + t279 * t297 + t212;
t105 = t208 * t327 + t209 * t232 + t210 * t233;
t104 = -t208 * t326 + t230 * t209 + t231 * t210;
t95 = t286 + t128;
t94 = -t127 + t285;
t83 = (t170 * t279 + t171 * t280) * t274 + t317;
t82 = -t128 * t244 + t178 * t232;
t81 = t127 * t244 - t178 * t230;
t80 = t165 * t276 + t167 * t244 + t169 * t245;
t79 = t164 * t276 + t166 * t244 + t168 * t245;
t77 = t78 * t276;
t76 = t78 * t244;
t75 = -t127 * t232 + t128 * t230;
t72 = (-t127 + t318) * t276 + t280 * t294;
t71 = t276 * t128 + t279 * t294 + t319;
t70 = t286 + t322;
t69 = t285 + t323;
t68 = -t118 * t203 + t234 * t93;
t67 = t118 * t201 - t234 * t92;
t64 = t286 + t350;
t63 = t312 - t202 * t267 + (-rSges(7,3) + t277) * t201 + t285 + t290;
t62 = (t127 * t279 + t128 * t280) * t274 + t295;
t60 = t61 * t276;
t59 = -t122 * t232 - t124 * t203 + t126 * t204;
t58 = -t121 * t232 - t123 * t203 + t125 * t204;
t57 = -t122 * t230 - t124 * t201 + t126 * t202;
t56 = -t121 * t230 - t123 * t201 + t125 * t202;
t55 = t61 * t244;
t53 = t54 * t276;
t52 = -t201 * t93 + t203 * t92;
t51 = t54 * t244;
t49 = -t244 * t103 - t137 + (t133 + t186) * t232;
t48 = -t230 * t133 - t244 * t323 - t154;
t43 = (-t102 + t308) * t276 + t280 * t289;
t42 = t276 * t103 + t279 * t289 + t309;
t41 = -t232 * t102 + t230 * t322 - t129;
t36 = t101 * t162 + t161 * t99 + t203 * t97;
t35 = t100 * t162 + t161 * t98 + t203 * t96;
t34 = t101 * t160 + t159 * t99 + t201 * t97;
t33 = t100 * t160 + t159 * t98 + t201 * t96;
t32 = (t102 * t279 + t103 * t280) * t274 + t288;
t31 = t152 * t89 + t153 * t91 + t203 * t87;
t30 = t152 * t88 + t153 * t90 + t203 * t86;
t27 = -t137 - t333 * t244 + (t186 - t320) * t232;
t26 = -t154 + t320 * t230 - (-t145 - t334) * t244;
t25 = (t308 - t334) * t276 + t280 * t287;
t24 = t276 * t333 + t279 * t287 + t309;
t23 = t77 + (t66 * t279 - t65 * t280) * t274;
t22 = -t65 * t230 - t66 * t232 - t76;
t21 = -t129 - t334 * t232 + (t146 + t333) * t230;
t20 = (t279 * t334 + t280 * t333) * t274 + t288;
t19 = t74 * t276 + (t279 * t59 - t280 * t58) * t274;
t18 = t73 * t276 + (t279 * t57 - t280 * t56) * t274;
t17 = -t230 * t58 - t232 * t59 - t244 * t74;
t16 = -t230 * t56 - t232 * t57 - t244 * t73;
t15 = t60 + (t40 * t279 - t39 * t280) * t274;
t14 = -t39 * t230 - t40 * t232 - t55;
t13 = t53 + (t38 * t279 - t37 * t280) * t274;
t12 = -t37 * t230 - t38 * t232 - t51;
t10 = t47 * t276 + (t279 * t36 - t280 * t35) * t274;
t9 = t46 * t276 + (t279 * t34 - t280 * t33) * t274;
t8 = -t230 * t35 - t232 * t36 - t244 * t47;
t7 = -t230 * t33 - t232 * t34 - t244 * t46;
t6 = t45 * t276 + (t279 * t31 - t280 * t30) * t274;
t5 = t44 * t276 + (t279 * t29 - t28 * t280) * t274;
t4 = -t230 * t30 - t232 * t31 - t244 * t45;
t3 = -t230 * t28 - t232 * t29 - t244 * t44;
t2 = t201 * t30 + t203 * t31 + t234 * t45;
t84 = [m(7) * (t63 ^ 2 + t64 ^ 2) + m(5) * (t94 ^ 2 + t95 ^ 2) + m(6) * (t69 ^ 2 + t70 ^ 2) + m(4) * (t138 ^ 2 + t139 ^ 2) + m(3) * (t199 ^ 2 + t200 ^ 2) + m(2) * (t258 ^ 2 + t259 ^ 2) + Icges(2,3) + t78 + t61 + t54 + t352; t60 + t53 + t77 + m(7) * (t24 * t64 + t25 * t63) + m(6) * (t42 * t70 + t43 * t69) + m(5) * (t71 * t95 + t72 * t94) + m(4) * (t108 * t139 + t109 * t138) + m(3) * (t182 * t200 + t183 * t199) + ((-t134 / 0.2e1 - t79 / 0.2e1 - t104 / 0.2e1 - t157 / 0.2e1 + t284) * t280 + (t135 / 0.2e1 + t80 / 0.2e1 + t105 / 0.2e1 + t158 / 0.2e1 - t283) * t279) * t274 + t321; (t13 + t15 + t23 + t321) * t276 + m(7) * (t20 ^ 2 + t24 ^ 2 + t25 ^ 2) + m(6) * (t32 ^ 2 + t42 ^ 2 + t43 ^ 2) + m(5) * (t62 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(4) * (t108 ^ 2 + t109 ^ 2 + t83 ^ 2) + m(3) * (t163 ^ 2 + t182 ^ 2 + t183 ^ 2) + ((-t18 - t5 - t9 + ((t230 * t166 + t231 * t168 + t251 * t216 + t252 * t218) * t274 + t316 * t349 * t280) * t280 + (-t104 - t134 - t157 - t79) * t276) * t280 + (t6 + t19 + t10 + ((t232 * t167 + t233 * t169 + t253 * t217 + t254 * t219) * t274 + t315 * t349 * t279) * t279 + (t158 + t105 + t135 + t80) * t276 + (-t232 * t166 - t230 * t167 - t233 * t168 - t231 * t169 - t253 * t216 - t251 * t217 - t254 * t218 - t252 * t219 + (t279 * t316 + t280 * t315) * t274) * t326) * t279) * t274; 0.2e1 * ((t279 * t63 - t280 * t64) * t345 + (t279 * t94 - t280 * t95) * t347 + (t279 * t69 - t280 * t70) * t346 + (t138 * t279 - t139 * t280) * t348) * t274; m(7) * (t276 * t20 + (-t24 * t280 + t25 * t279) * t274) + m(6) * (t276 * t32 + (t279 * t43 - t280 * t42) * t274) + m(5) * (t276 * t62 + (t279 * t72 - t280 * t71) * t274) + m(4) * (t276 * t83 + (-t108 * t280 + t109 * t279) * t274); 0.2e1 * (t348 + t347 + t313) * (t276 ^ 2 + (t279 ^ 2 + t280 ^ 2) * t349); -t51 - t55 - t76 + m(7) * (t26 * t63 + t27 * t64) + m(5) * (t81 * t94 + t82 * t95) + m(6) * (t48 * t69 + t49 * t70) + t283 * t232 + t284 * t230; (t12 / 0.2e1 + t14 / 0.2e1 + t22 / 0.2e1) * t276 - (t13 / 0.2e1 + t23 / 0.2e1 + t15 / 0.2e1) * t244 + (-t6 / 0.2e1 - t10 / 0.2e1 - t19 / 0.2e1) * t232 + (-t5 / 0.2e1 - t18 / 0.2e1 - t9 / 0.2e1) * t230 + m(7) * (t20 * t21 + t24 * t27 + t25 * t26) + m(6) * (t32 * t41 + t42 * t49 + t43 * t48) + m(5) * (t62 * t75 + t71 * t82 + t72 * t81) + ((-t3 / 0.2e1 - t16 / 0.2e1 - t7 / 0.2e1) * t280 + (t4 / 0.2e1 + t17 / 0.2e1 + t8 / 0.2e1) * t279) * t274; m(5) * (t75 * t276 + (t279 * t81 - t280 * t82) * t274) + m(6) * (t41 * t276 + (t279 * t48 - t280 * t49) * t274) + m(7) * (t21 * t276 + (t26 * t279 - t27 * t280) * t274); -(t12 + t14 + t22) * t244 + (-t4 - t17 - t8) * t232 + (-t3 - t7 - t16) * t230 + m(7) * (t21 ^ 2 + t26 ^ 2 + t27 ^ 2) + m(6) * (t41 ^ 2 + t48 ^ 2 + t49 ^ 2) + m(5) * (t75 ^ 2 + t81 ^ 2 + t82 ^ 2); m(7) * (t201 * t64 + t203 * t63) + m(6) * (t201 * t70 + t203 * t69); m(7) * (t20 * t234 + t201 * t24 + t203 * t25) + m(6) * (t201 * t42 + t203 * t43 + t234 * t32); (t234 * t276 + (-t201 * t280 + t203 * t279) * t274) * t351; m(7) * (t201 * t27 + t203 * t26 + t21 * t234) + m(6) * (t201 * t49 + t203 * t48 + t234 * t41); (t201 ^ 2 + t203 ^ 2 + t234 ^ 2) * t351; m(7) * (t63 * t67 + t64 * t68) + t50 + t310 * t203 + t311 * t201; t276 * t343 + t5 * t342 + t6 * t341 + t13 * t340 + m(7) * (t20 * t52 + t24 * t68 + t25 * t67) + (t279 * t2 / 0.2e1 + t280 * t344) * t274; m(7) * (t52 * t276 + (t279 * t67 - t280 * t68) * t274); t230 * t344 + t12 * t340 + m(7) * (t21 * t52 + t26 * t67 + t27 * t68) + t4 * t341 - t232 * t2 / 0.2e1 - t244 * t343 + t3 * t342; m(7) * (t201 * t68 + t203 * t67 + t234 * t52); t203 * t2 + t201 * t1 + t234 * t11 + m(7) * (t52 ^ 2 + t67 ^ 2 + t68 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t84(1) t84(2) t84(4) t84(7) t84(11) t84(16); t84(2) t84(3) t84(5) t84(8) t84(12) t84(17); t84(4) t84(5) t84(6) t84(9) t84(13) t84(18); t84(7) t84(8) t84(9) t84(10) t84(14) t84(19); t84(11) t84(12) t84(13) t84(14) t84(15) t84(20); t84(16) t84(17) t84(18) t84(19) t84(20) t84(21);];
Mq  = res;
