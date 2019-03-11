% Calculate joint inertia matrix for
% S6RRPRRP8
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRP8_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP8_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP8_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:22:00
% EndTime: 2019-03-09 12:22:14
% DurationCPUTime: 6.22s
% Computational Cost: add. (12645->512), mult. (12871->731), div. (0->0), fcn. (13924->10), ass. (0->242)
t243 = pkin(10) + qJ(4);
t236 = qJ(5) + t243;
t231 = sin(t236);
t232 = cos(t236);
t253 = cos(qJ(1));
t251 = sin(qJ(1));
t252 = cos(qJ(2));
t301 = t251 * t252;
t175 = t231 * t301 + t232 * t253;
t176 = -t231 * t253 + t232 * t301;
t250 = sin(qJ(2));
t304 = t250 * t251;
t107 = Icges(7,5) * t176 + Icges(7,6) * t304 + Icges(7,3) * t175;
t113 = Icges(6,4) * t176 - Icges(6,2) * t175 + Icges(6,6) * t304;
t354 = t107 - t113;
t300 = t252 * t253;
t177 = t231 * t300 - t251 * t232;
t178 = t251 * t231 + t232 * t300;
t303 = t250 * t253;
t108 = Icges(7,5) * t178 + Icges(7,6) * t303 + Icges(7,3) * t177;
t114 = Icges(6,4) * t178 - Icges(6,2) * t177 + Icges(6,6) * t303;
t353 = t108 - t114;
t109 = Icges(6,5) * t176 - Icges(6,6) * t175 + Icges(6,3) * t304;
t111 = Icges(7,4) * t176 + Icges(7,2) * t304 + Icges(7,6) * t175;
t352 = t109 + t111;
t110 = Icges(6,5) * t178 - Icges(6,6) * t177 + Icges(6,3) * t303;
t112 = Icges(7,4) * t178 + Icges(7,2) * t303 + Icges(7,6) * t177;
t351 = t110 + t112;
t115 = Icges(7,1) * t176 + Icges(7,4) * t304 + Icges(7,5) * t175;
t117 = Icges(6,1) * t176 - Icges(6,4) * t175 + Icges(6,5) * t304;
t350 = t115 + t117;
t116 = Icges(7,1) * t178 + Icges(7,4) * t303 + Icges(7,5) * t177;
t118 = Icges(6,1) * t178 - Icges(6,4) * t177 + Icges(6,5) * t303;
t349 = t116 + t118;
t330 = rSges(7,3) + qJ(6);
t338 = rSges(7,1) + pkin(5);
t348 = -t330 * t175 - t338 * t176;
t347 = t175 * t354 + t350 * t176 + t352 * t304;
t346 = t175 * t353 + t176 * t349 + t304 * t351;
t345 = t177 * t354 + t350 * t178 + t352 * t303;
t344 = t177 * t353 + t178 * t349 + t303 * t351;
t157 = -Icges(7,6) * t252 + (Icges(7,5) * t232 + Icges(7,3) * t231) * t250;
t159 = -Icges(7,2) * t252 + (Icges(7,4) * t232 + Icges(7,6) * t231) * t250;
t161 = -Icges(7,4) * t252 + (Icges(7,1) * t232 + Icges(7,5) * t231) * t250;
t71 = t157 * t175 + t159 * t304 + t161 * t176;
t158 = -Icges(6,3) * t252 + (Icges(6,5) * t232 - Icges(6,6) * t231) * t250;
t160 = -Icges(6,6) * t252 + (Icges(6,4) * t232 - Icges(6,2) * t231) * t250;
t162 = -Icges(6,5) * t252 + (Icges(6,1) * t232 - Icges(6,4) * t231) * t250;
t72 = t158 * t304 - t160 * t175 + t162 * t176;
t343 = -t72 - t71;
t73 = t177 * t157 + t159 * t303 + t178 * t161;
t74 = t158 * t303 - t177 * t160 + t178 * t162;
t342 = -t73 - t74;
t307 = t231 * t250;
t339 = t157 * t307 + (t161 + t162) * t232 * t250;
t340 = -t158 - t159;
t331 = -t160 * t307 + t340 * t252 + t339;
t341 = t331 * t252;
t323 = t251 / 0.2e1;
t320 = t253 / 0.2e1;
t337 = t343 * t252 + (t347 * t251 + t346 * t253) * t250;
t336 = t342 * t252 + (t345 * t251 + t344 * t253) * t250;
t335 = t346 * t251 - t347 * t253;
t334 = t344 * t251 - t345 * t253;
t55 = -t252 * t111 + (t107 * t231 + t115 * t232) * t250;
t57 = -t252 * t109 + (-t113 * t231 + t117 * t232) * t250;
t333 = -t55 - t57;
t56 = -t252 * t112 + (t108 * t231 + t116 * t232) * t250;
t58 = -t252 * t110 + (-t114 * t231 + t118 * t232) * t250;
t332 = t56 + t58;
t329 = rSges(7,2) * t304 - t348;
t328 = -t252 * rSges(7,2) + (t330 * t231 + t338 * t232) * t250;
t245 = t251 ^ 2;
t246 = t253 ^ 2;
t327 = m(4) / 0.2e1;
t326 = m(5) / 0.2e1;
t325 = m(6) / 0.2e1;
t324 = m(7) / 0.2e1;
t322 = -t252 / 0.2e1;
t321 = -t253 / 0.2e1;
t217 = rSges(3,1) * t250 + rSges(3,2) * t252;
t319 = m(3) * t217;
t318 = pkin(2) * t252;
t248 = cos(pkin(10));
t233 = t248 * pkin(3) + pkin(2);
t317 = -pkin(2) + t233;
t249 = -pkin(8) - qJ(3);
t316 = t341 + (t251 * t333 - t253 * t332) * t250;
t314 = t253 * rSges(3,3);
t122 = t178 * rSges(6,1) - t177 * rSges(6,2) + rSges(6,3) * t303;
t242 = -pkin(9) + t249;
t285 = t242 - t249;
t247 = sin(pkin(10));
t302 = t251 * t247;
t289 = -pkin(3) * t302 - t233 * t300;
t235 = cos(t243);
t208 = pkin(4) * t235 + t233;
t234 = sin(t243);
t210 = pkin(3) * t247 + pkin(4) * t234;
t293 = t208 * t300 + t251 * t210;
t98 = -t285 * t303 + t289 + t293;
t313 = -t122 - t98;
t290 = t208 - t233;
t145 = t250 * t290 + t252 * t285;
t305 = t247 * t253;
t288 = pkin(3) * t305 + t249 * t304;
t291 = -t253 * t210 - t242 * t304;
t97 = t290 * t301 + t288 + t291;
t312 = t145 * t304 + t252 * t97;
t311 = Icges(3,4) * t250;
t310 = Icges(3,4) * t252;
t309 = qJ(3) * t250;
t167 = -Icges(5,6) * t252 + (Icges(5,4) * t235 - Icges(5,2) * t234) * t250;
t308 = t167 * t234;
t299 = t329 * t303;
t298 = rSges(7,2) * t303 + t177 * t330 + t178 * t338;
t267 = -t176 * rSges(6,1) + t175 * rSges(6,2);
t120 = rSges(6,3) * t304 - t267;
t164 = -t252 * rSges(6,3) + (rSges(6,1) * t232 - rSges(6,2) * t231) * t250;
t87 = t252 * t120 + t164 * t304;
t216 = t250 * pkin(2) - t252 * qJ(3);
t295 = -(qJ(3) + t249) * t252 - t317 * t250 - t216;
t294 = t252 * rSges(4,3) - (rSges(4,1) * t248 - rSges(4,2) * t247) * t250 - t216;
t287 = pkin(2) * t300 + qJ(3) * t303;
t292 = t245 * (t309 + t318) + t253 * t287;
t286 = t253 * pkin(1) + t251 * pkin(7);
t284 = t245 + t246;
t186 = -t234 * t301 - t235 * t253;
t187 = -t234 * t253 + t235 * t301;
t125 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t304;
t127 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t304;
t129 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t304;
t59 = -t252 * t125 + (-t127 * t234 + t129 * t235) * t250;
t166 = -Icges(5,3) * t252 + (Icges(5,5) * t235 - Icges(5,6) * t234) * t250;
t168 = -Icges(5,5) * t252 + (Icges(5,1) * t235 - Icges(5,4) * t234) * t250;
t77 = t166 * t304 + t167 * t186 + t168 * t187;
t283 = t77 / 0.2e1 + t59 / 0.2e1;
t188 = -t234 * t300 + t251 * t235;
t189 = t251 * t234 + t235 * t300;
t126 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t303;
t128 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t303;
t130 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t303;
t60 = -t252 * t126 + (-t128 * t234 + t130 * t235) * t250;
t78 = t166 * t303 + t188 * t167 + t189 * t168;
t282 = t78 / 0.2e1 + t60 / 0.2e1;
t281 = -t98 - t298;
t280 = t303 * t336 + t304 * t337;
t279 = -t145 + t295;
t174 = -t252 * rSges(5,3) + (rSges(5,1) * t235 - rSges(5,2) * t234) * t250;
t278 = -t174 + t295;
t132 = t189 * rSges(5,1) + t188 * rSges(5,2) + rSges(5,3) * t303;
t204 = -t247 * t300 + t251 * t248;
t205 = t248 * t300 + t302;
t277 = t205 * rSges(4,1) + t204 * rSges(4,2) + rSges(4,3) * t303;
t240 = t253 * pkin(7);
t276 = t240 - t291;
t275 = t304 / 0.2e1;
t274 = t303 / 0.2e1;
t273 = -t208 * t252 - pkin(1);
t259 = -t249 * t303 - t289;
t272 = t251 * ((t252 * t317 - t309) * t251 - t288) + t253 * (t259 - t287) + t292;
t271 = -t164 + t279;
t61 = t252 * t329 + t304 * t328;
t270 = rSges(3,1) * t252 - rSges(3,2) * t250;
t202 = -t247 * t301 - t248 * t253;
t203 = t248 * t301 - t305;
t269 = -t203 * rSges(4,1) - t202 * rSges(4,2);
t268 = -t187 * rSges(5,1) - t186 * rSges(5,2);
t266 = t279 - t328;
t265 = Icges(3,1) * t252 - t311;
t264 = -Icges(3,2) * t250 + t310;
t263 = Icges(3,5) * t252 - Icges(3,6) * t250;
t260 = rSges(3,1) * t300 - rSges(3,2) * t303 + t251 * rSges(3,3);
t258 = t251 * t97 + t253 * t98 + t272;
t257 = t252 * t316 + t280;
t256 = (-t333 - t343) * t275 + (t332 - t342) * t274;
t255 = -t242 * t303 + t286 + t293;
t254 = t336 * t323 + (t251 * t332 + t253 * t333) * t322 + t337 * t321 + t335 * t275 + t334 * t274;
t244 = t250 ^ 2;
t219 = rSges(2,1) * t253 - t251 * rSges(2,2);
t218 = -t251 * rSges(2,1) - rSges(2,2) * t253;
t212 = Icges(3,5) * t250 + Icges(3,6) * t252;
t191 = Icges(3,3) * t251 + t253 * t263;
t190 = -Icges(3,3) * t253 + t251 * t263;
t184 = -Icges(4,5) * t252 + (Icges(4,1) * t248 - Icges(4,4) * t247) * t250;
t183 = -Icges(4,6) * t252 + (Icges(4,4) * t248 - Icges(4,2) * t247) * t250;
t156 = t260 + t286;
t155 = t314 + t240 + (-pkin(1) - t270) * t251;
t151 = t250 * t235 * t168;
t147 = t294 * t253;
t146 = t294 * t251;
t144 = Icges(4,1) * t205 + Icges(4,4) * t204 + Icges(4,5) * t303;
t143 = Icges(4,1) * t203 + Icges(4,4) * t202 + Icges(4,5) * t304;
t142 = Icges(4,4) * t205 + Icges(4,2) * t204 + Icges(4,6) * t303;
t141 = Icges(4,4) * t203 + Icges(4,2) * t202 + Icges(4,6) * t304;
t140 = Icges(4,5) * t205 + Icges(4,6) * t204 + Icges(4,3) * t303;
t139 = Icges(4,5) * t203 + Icges(4,6) * t202 + Icges(4,3) * t304;
t137 = t253 * t260 + (t251 * t270 - t314) * t251;
t131 = rSges(5,3) * t304 - t268;
t104 = t120 * t303;
t102 = t277 + t286 + t287;
t101 = t240 + (-t318 - pkin(1) + (-rSges(4,3) - qJ(3)) * t250) * t251 + t269;
t100 = t278 * t253;
t99 = t278 * t251;
t93 = t97 * t303;
t92 = -t252 * t132 - t174 * t303;
t91 = t131 * t252 + t174 * t304;
t90 = t259 + t132 + t286;
t89 = t240 + (-rSges(5,3) * t250 - t233 * t252 - pkin(1)) * t251 + t268 + t288;
t88 = -t252 * t122 - t164 * t303;
t86 = -t252 * t166 - t250 * t308 + t151;
t83 = t255 + t122;
t82 = (-rSges(6,3) * t250 + t273) * t251 + t267 + t276;
t81 = t271 * t253;
t80 = t271 * t251;
t79 = (t131 * t253 - t132 * t251) * t250;
t76 = -t122 * t304 + t104;
t75 = t251 * (rSges(4,3) * t304 - t269) + t253 * t277 + t292;
t70 = t266 * t253;
t69 = t266 * t251;
t64 = t255 + t298;
t63 = (-rSges(7,2) * t250 + t273) * t251 + t276 + t348;
t62 = -t252 * t298 - t303 * t328;
t50 = t126 * t303 + t188 * t128 + t189 * t130;
t49 = t125 * t303 + t188 * t127 + t189 * t129;
t48 = t126 * t304 + t128 * t186 + t130 * t187;
t47 = t125 * t304 + t127 * t186 + t129 * t187;
t46 = t313 * t252 + (-t145 - t164) * t303;
t45 = t87 + t312;
t36 = -t298 * t304 + t299;
t35 = t251 * t131 + t132 * t253 + t272;
t34 = t304 * t313 + t104 + t93;
t33 = t281 * t252 + (-t145 - t328) * t303;
t32 = t61 + t312;
t31 = t281 * t304 + t299 + t93;
t28 = t251 * t120 + t122 * t253 + t258;
t27 = t50 * t251 - t253 * t49;
t26 = t48 * t251 - t253 * t47;
t17 = t251 * t329 + t298 * t253 + t258;
t14 = -t78 * t252 + (t251 * t49 + t253 * t50) * t250;
t13 = -t77 * t252 + (t251 * t47 + t253 * t48) * t250;
t1 = [Icges(2,3) + t151 + (-t166 - (Icges(4,5) * t248 - Icges(4,6) * t247) * t250 + t311 + (Icges(4,3) + Icges(3,2)) * t252 + t340) * t252 + (Icges(3,1) * t250 - t160 * t231 - t183 * t247 + t184 * t248 - t308 + t310) * t250 + m(7) * (t63 ^ 2 + t64 ^ 2) + m(6) * (t82 ^ 2 + t83 ^ 2) + m(5) * (t89 ^ 2 + t90 ^ 2) + m(4) * (t101 ^ 2 + t102 ^ 2) + m(3) * (t155 ^ 2 + t156 ^ 2) + m(2) * (t218 ^ 2 + t219 ^ 2) + t339; m(4) * (t101 * t147 + t102 * t146) + m(7) * (t63 * t70 + t64 * t69) + m(6) * (t80 * t83 + t81 * t82) + m(5) * (t100 * t89 + t90 * t99) + (-t202 * t183 / 0.2e1 - t203 * t184 / 0.2e1 - t155 * t319 - t57 / 0.2e1 - t55 / 0.2e1 - t71 / 0.2e1 - t72 / 0.2e1 + t212 * t320 + (t139 / 0.2e1 + Icges(3,6) * t320 - t251 * t264 / 0.2e1) * t252 - t283) * t253 + (t204 * t183 / 0.2e1 + t205 * t184 / 0.2e1 - t156 * t319 + t58 / 0.2e1 + t56 / 0.2e1 + t73 / 0.2e1 + t74 / 0.2e1 + t212 * t323 + (-t140 / 0.2e1 + Icges(3,6) * t323 + t264 * t320) * t252 + t282) * t251 + ((Icges(3,5) * t251 - t142 * t247 + t144 * t248 + t253 * t265) * t323 + (-Icges(3,5) * t253 - t141 * t247 + t143 * t248 + t251 * t265) * t321) * t250; m(7) * (t17 ^ 2 + t69 ^ 2 + t70 ^ 2) + m(6) * (t28 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(5) * (t100 ^ 2 + t35 ^ 2 + t99 ^ 2) + m(4) * (t146 ^ 2 + t147 ^ 2 + t75 ^ 2) + m(3) * (t217 ^ 2 * t284 + t137 ^ 2) + (-t246 * t190 - t26 + (t139 * t304 + t202 * t141 + t203 * t143) * t253 - t335) * t253 + (t27 + t245 * t191 + (t140 * t303 + t204 * t142 + t205 * t144) * t251 + (-t139 * t303 - t140 * t304 - t204 * t141 - t142 * t202 - t205 * t143 - t144 * t203 - t251 * t190 + t253 * t191) * t253 + t334) * t251; 0.2e1 * ((t251 * t64 + t253 * t63) * t324 + (t251 * t83 + t253 * t82) * t325 + (t251 * t90 + t253 * t89) * t326 + (t101 * t253 + t102 * t251) * t327) * t250; m(7) * (-t252 * t17 + (t251 * t69 + t253 * t70) * t250) + m(6) * (-t252 * t28 + (t251 * t80 + t253 * t81) * t250) + m(5) * (-t252 * t35 + (t100 * t253 + t251 * t99) * t250) + m(4) * (-t252 * t75 + (t146 * t251 + t147 * t253) * t250); 0.2e1 * (t327 + t326 + t325 + t324) * (t244 * t284 + t252 ^ 2); (-t86 - t331) * t252 + m(7) * (t32 * t63 + t33 * t64) + m(6) * (t45 * t82 + t46 * t83) + m(5) * (t89 * t91 + t90 * t92) + (t251 * t283 + t253 * t282) * t250 + t256; m(7) * (t17 * t31 + t32 * t70 + t33 * t69) + m(6) * (t28 * t34 + t45 * t81 + t46 * t80) + m(5) * (t100 * t91 + t35 * t79 + t92 * t99) + (t26 * t323 + t27 * t320) * t250 + t254 + t14 * t323 + t13 * t321 + (t60 * t251 - t59 * t253) * t322; m(5) * (-t79 * t252 + (t251 * t92 + t253 * t91) * t250) + m(6) * (-t34 * t252 + (t251 * t46 + t253 * t45) * t250) + m(7) * (-t31 * t252 + (t251 * t33 + t253 * t32) * t250); (t86 * t252 + t316) * t252 + (t253 * t14 + t251 * t13 - t252 * (t251 * t59 + t253 * t60)) * t250 + m(7) * (t31 ^ 2 + t32 ^ 2 + t33 ^ 2) + m(6) * (t34 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(5) * (t79 ^ 2 + t91 ^ 2 + t92 ^ 2) + t280; -t341 + m(7) * (t61 * t63 + t62 * t64) + m(6) * (t82 * t87 + t83 * t88) + t256; m(7) * (t17 * t36 + t61 * t70 + t62 * t69) + m(6) * (t28 * t76 + t80 * t88 + t81 * t87) + t254; m(6) * (-t76 * t252 + (t251 * t88 + t253 * t87) * t250) + m(7) * (-t36 * t252 + (t251 * t62 + t253 * t61) * t250); m(7) * (t31 * t36 + t32 * t61 + t33 * t62) + m(6) * (t34 * t76 + t45 * t87 + t46 * t88) + t257; m(7) * (t36 ^ 2 + t61 ^ 2 + t62 ^ 2) + m(6) * (t76 ^ 2 + t87 ^ 2 + t88 ^ 2) + t257; m(7) * (t175 * t64 + t177 * t63); m(7) * (t17 * t307 + t175 * t69 + t177 * t70); m(7) * (t175 * t251 + t177 * t253 - t231 * t252) * t250; m(7) * (t175 * t33 + t177 * t32 + t307 * t31); m(7) * (t175 * t62 + t177 * t61 + t307 * t36); m(7) * (t231 ^ 2 * t244 + t175 ^ 2 + t177 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
