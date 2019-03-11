% Calculate joint inertia matrix for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP3_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP3_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP3_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:50
% EndTime: 2019-03-09 20:54:02
% DurationCPUTime: 5.06s
% Computational Cost: add. (9377->467), mult. (12432->664), div. (0->0), fcn. (13592->8), ass. (0->226)
t241 = qJ(2) + qJ(3);
t235 = cos(t241);
t245 = cos(qJ(4));
t247 = cos(qJ(1));
t300 = t247 * t245;
t242 = sin(qJ(4));
t244 = sin(qJ(1));
t304 = t244 * t242;
t203 = t235 * t304 + t300;
t301 = t247 * t242;
t303 = t244 * t245;
t204 = t235 * t303 - t301;
t234 = sin(t241);
t309 = t234 * t244;
t115 = Icges(7,5) * t309 + Icges(7,6) * t203 + Icges(7,3) * t204;
t121 = Icges(6,4) * t309 - Icges(6,2) * t204 + Icges(6,6) * t203;
t131 = Icges(5,1) * t204 - Icges(5,4) * t203 + Icges(5,5) * t309;
t355 = t115 - t121 + t131;
t205 = t235 * t301 - t303;
t206 = t235 * t300 + t304;
t307 = t234 * t247;
t116 = Icges(7,5) * t307 + Icges(7,6) * t205 + Icges(7,3) * t206;
t122 = Icges(6,4) * t307 - Icges(6,2) * t206 + Icges(6,6) * t205;
t132 = Icges(5,1) * t206 - Icges(5,4) * t205 + Icges(5,5) * t307;
t354 = t116 - t122 + t132;
t117 = Icges(6,5) * t309 - Icges(6,6) * t204 + Icges(6,3) * t203;
t119 = Icges(7,4) * t309 + Icges(7,2) * t203 + Icges(7,6) * t204;
t129 = Icges(5,4) * t204 - Icges(5,2) * t203 + Icges(5,6) * t309;
t353 = t117 + t119 - t129;
t118 = Icges(6,5) * t307 - Icges(6,6) * t206 + Icges(6,3) * t205;
t120 = Icges(7,4) * t307 + Icges(7,2) * t205 + Icges(7,6) * t206;
t130 = Icges(5,4) * t206 - Icges(5,2) * t205 + Icges(5,6) * t307;
t352 = t118 + t120 - t130;
t123 = Icges(7,1) * t309 + Icges(7,4) * t203 + Icges(7,5) * t204;
t125 = Icges(6,1) * t309 - Icges(6,4) * t204 + Icges(6,5) * t203;
t127 = Icges(5,5) * t204 - Icges(5,6) * t203 + Icges(5,3) * t309;
t351 = t123 + t125 + t127;
t124 = Icges(7,1) * t307 + Icges(7,4) * t205 + Icges(7,5) * t206;
t126 = Icges(6,1) * t307 - Icges(6,4) * t206 + Icges(6,5) * t205;
t128 = Icges(5,5) * t206 - Icges(5,6) * t205 + Icges(5,3) * t307;
t350 = t124 + t126 + t128;
t349 = rSges(7,1) + pkin(5);
t348 = rSges(7,3) + qJ(6);
t347 = t353 * t203 + t355 * t204 + t351 * t309;
t346 = t352 * t203 + t354 * t204 + t350 * t309;
t345 = t353 * t205 + t355 * t206 + t351 * t307;
t344 = t352 * t205 + t354 * t206 + t350 * t307;
t169 = -Icges(7,5) * t235 + (Icges(7,6) * t242 + Icges(7,3) * t245) * t234;
t171 = -Icges(7,4) * t235 + (Icges(7,2) * t242 + Icges(7,6) * t245) * t234;
t173 = -Icges(7,1) * t235 + (Icges(7,4) * t242 + Icges(7,5) * t245) * t234;
t79 = t206 * t169 + t205 * t171 + t173 * t307;
t170 = -Icges(6,5) * t235 + (-Icges(6,6) * t245 + Icges(6,3) * t242) * t234;
t172 = -Icges(6,4) * t235 + (-Icges(6,2) * t245 + Icges(6,6) * t242) * t234;
t174 = -Icges(6,1) * t235 + (-Icges(6,4) * t245 + Icges(6,5) * t242) * t234;
t80 = t205 * t170 - t206 * t172 + t174 * t307;
t166 = -Icges(5,3) * t235 + (Icges(5,5) * t245 - Icges(5,6) * t242) * t234;
t167 = -Icges(5,6) * t235 + (Icges(5,4) * t245 - Icges(5,2) * t242) * t234;
t168 = -Icges(5,5) * t235 + (Icges(5,1) * t245 - Icges(5,4) * t242) * t234;
t82 = t166 * t307 - t205 * t167 + t206 * t168;
t343 = -t80 - t82 - t79;
t77 = t204 * t169 + t203 * t171 + t173 * t309;
t78 = t203 * t170 - t204 * t172 + t174 * t309;
t81 = t166 * t309 - t203 * t167 + t204 * t168;
t342 = -t81 - t77 - t78;
t323 = t203 * rSges(7,2);
t298 = t348 * t204 + t349 * t309 + t323;
t340 = t205 * rSges(7,2) + t348 * t206 + t349 * t307;
t339 = -t242 * t167 - t245 * t172;
t338 = t342 * t235 + (t244 * t347 + t346 * t247) * t234;
t337 = t343 * t235 + (t244 * t345 + t247 * t344) * t234;
t336 = t346 * t244 - t247 * t347;
t335 = t244 * t344 - t247 * t345;
t334 = -t166 - t173 - t174;
t308 = t234 * t245;
t310 = t234 * t242;
t332 = (t170 + t171) * t310 + (t168 + t169) * t308;
t262 = Icges(4,5) * t235 - Icges(4,6) * t234;
t180 = -Icges(4,3) * t247 + t244 * t262;
t181 = Icges(4,3) * t244 + t247 * t262;
t240 = t247 ^ 2;
t311 = Icges(4,4) * t235;
t264 = -Icges(4,2) * t234 + t311;
t183 = Icges(4,6) * t244 + t247 * t264;
t312 = Icges(4,4) * t234;
t266 = Icges(4,1) * t235 - t312;
t185 = Icges(4,5) * t244 + t247 * t266;
t260 = -t183 * t234 + t185 * t235;
t182 = -Icges(4,6) * t247 + t244 * t264;
t184 = -Icges(4,5) * t247 + t244 * t266;
t261 = t182 * t234 - t184 * t235;
t331 = -t240 * t180 - (t260 * t244 + (-t181 + t261) * t247) * t244 - t336;
t239 = t244 ^ 2;
t329 = t244 / 0.2e1;
t328 = -t247 / 0.2e1;
t243 = sin(qJ(2));
t327 = pkin(2) * t243;
t326 = pkin(3) * t235;
t246 = cos(qJ(2));
t325 = rSges(3,1) * t246;
t324 = rSges(3,2) * t243;
t322 = t203 * rSges(6,3);
t321 = t247 * rSges(3,3);
t60 = -t235 * t127 + (-t129 * t242 + t131 * t245) * t234;
t320 = t60 * t247;
t61 = -t235 * t128 + (-t130 * t242 + t132 * t245) * t234;
t319 = t61 * t244;
t62 = -t235 * t123 + (t115 * t245 + t119 * t242) * t234;
t318 = t62 * t247;
t63 = -t235 * t124 + (t116 * t245 + t120 * t242) * t234;
t317 = t63 * t244;
t64 = -t235 * t125 + (t117 * t242 - t121 * t245) * t234;
t316 = t64 * t247;
t65 = -t235 * t126 + (t118 * t242 - t122 * t245) * t234;
t315 = t65 * t244;
t314 = Icges(3,4) * t243;
t313 = Icges(3,4) * t246;
t306 = t235 * t247;
t248 = -pkin(8) - pkin(7);
t299 = t247 * t248;
t136 = rSges(6,1) * t307 - t206 * rSges(6,2) + t205 * rSges(6,3);
t148 = t206 * pkin(4) + t205 * qJ(5);
t297 = -t136 - t148;
t194 = t203 * qJ(5);
t147 = t204 * pkin(4) + t194;
t202 = (pkin(4) * t245 + qJ(5) * t242) * t234;
t296 = t235 * t147 + t202 * t309;
t232 = t246 * pkin(2) + pkin(1);
t221 = t247 * t232;
t238 = t247 * pkin(7);
t294 = t244 * (t299 + t238 + (-pkin(1) + t232) * t244) + t247 * (-t247 * pkin(1) + t221 + (-pkin(7) - t248) * t244);
t255 = rSges(4,1) * t306 - rSges(4,2) * t307 + t244 * rSges(4,3);
t269 = rSges(4,1) * t235 - rSges(4,2) * t234;
t106 = t244 * (-t247 * rSges(4,3) + t244 * t269) + t247 * t255;
t175 = -t235 * rSges(5,3) + (rSges(5,1) * t245 - rSges(5,2) * t242) * t234;
t213 = t234 * pkin(3) - t235 * pkin(9);
t293 = -t175 - t213;
t292 = (rSges(7,2) * t242 + rSges(7,3) * t245) * t234 + qJ(6) * t308 - t349 * t235;
t177 = -t235 * rSges(6,1) + (-rSges(6,2) * t245 + rSges(6,3) * t242) * t234;
t291 = -t177 - t202;
t289 = pkin(3) * t306 + pkin(9) * t307;
t290 = t239 * (pkin(9) * t234 + t326) + t247 * t289;
t288 = t244 * rSges(3,3) + t247 * t325;
t287 = t239 + t240;
t286 = t339 * t234 + t334 * t235 + t332;
t285 = -t148 - t340;
t284 = -t202 - t292;
t283 = -t213 + t291;
t138 = t206 * rSges(5,1) - t205 * rSges(5,2) + rSges(5,3) * t307;
t212 = t234 * rSges(4,1) + rSges(4,2) * t235;
t280 = -t212 - t327;
t279 = -t213 - t327;
t278 = -t232 - t326;
t277 = (t239 * t181 + (t261 * t247 + (-t180 + t260) * t244) * t247 + t335) * t244;
t276 = -t194 - t299;
t275 = -t244 * t248 + t221;
t274 = t244 * t147 + t247 * t148 + t290;
t273 = -t213 + t284;
t268 = -t204 * rSges(5,1) + t203 * rSges(5,2);
t137 = rSges(5,3) * t309 - t268;
t70 = t244 * t137 + t247 * t138 + t290;
t272 = -t175 + t279;
t271 = -t202 + t279;
t270 = -t324 + t325;
t267 = Icges(3,1) * t246 - t314;
t265 = -Icges(3,2) * t243 + t313;
t263 = Icges(3,5) * t246 - Icges(3,6) * t243;
t210 = Icges(4,2) * t235 + t312;
t211 = Icges(4,1) * t234 + t311;
t257 = -t210 * t234 + t211 * t235;
t256 = -t177 + t271;
t254 = t275 + t289;
t134 = rSges(6,1) * t309 - t204 * rSges(6,2) + t322;
t35 = t244 * t134 + t247 * t136 + t274;
t253 = t271 - t292;
t32 = t298 * t244 + t340 * t247 + t274;
t252 = t331 * t247 + t277;
t251 = t148 + t254;
t250 = -(t319 - t320 + t317 - t318 + t315 - t316) * t235 / 0.2e1 + t337 * t329 + t338 * t328 + t336 * t309 / 0.2e1 + t335 * t307 / 0.2e1;
t209 = Icges(4,5) * t234 + Icges(4,6) * t235;
t249 = -t320 / 0.2e1 + t319 / 0.2e1 - t318 / 0.2e1 + t317 / 0.2e1 - t316 / 0.2e1 + t315 / 0.2e1 + (t235 * t183 + t234 * t185 + t244 * t209 + t247 * t257 - t343) * t329 + (t235 * t182 + t234 * t184 - t247 * t209 + t244 * t257 - t342) * t328;
t233 = t234 ^ 2;
t220 = t247 * rSges(2,1) - t244 * rSges(2,2);
t219 = -t244 * rSges(2,1) - t247 * rSges(2,2);
t218 = t243 * rSges(3,1) + t246 * rSges(3,2);
t189 = Icges(3,3) * t244 + t247 * t263;
t188 = -Icges(3,3) * t247 + t244 * t263;
t179 = t280 * t247;
t178 = t280 * t244;
t160 = t244 * pkin(7) + (pkin(1) - t324) * t247 + t288;
t159 = t321 + t238 + (-pkin(1) - t270) * t244;
t150 = t255 + t275;
t149 = (rSges(4,3) - t248) * t247 + (-t232 - t269) * t244;
t142 = t293 * t247;
t141 = t293 * t244;
t140 = t247 * (-t247 * t324 + t288) + (t244 * t270 - t321) * t244;
t139 = t147 * t307;
t108 = t272 * t247;
t107 = t272 * t244;
t101 = t283 * t247;
t100 = t283 * t244;
t99 = t256 * t247;
t98 = t256 * t244;
t97 = t254 + t138;
t96 = -t299 + ((-rSges(5,3) - pkin(9)) * t234 + t278) * t244 + t268;
t95 = -t235 * t138 - t175 * t307;
t94 = t235 * t137 + t175 * t309;
t93 = t273 * t247;
t92 = t273 * t244;
t91 = t253 * t247;
t90 = t253 * t244;
t86 = t106 + t294;
t85 = (t137 * t247 - t138 * t244) * t234;
t84 = t251 + t136;
t83 = -t322 + (rSges(6,2) - pkin(4)) * t204 + ((-rSges(6,1) - pkin(9)) * t234 + t278) * t244 + t276;
t69 = t251 + t340;
t68 = -t323 + (-pkin(4) - t348) * t204 + ((-pkin(9) - t349) * t234 + t278) * t244 + t276;
t67 = t235 * t297 + t291 * t307;
t66 = t235 * t134 + t177 * t309 + t296;
t39 = t70 + t294;
t38 = t139 + (t134 * t247 + t244 * t297) * t234;
t37 = t235 * t285 + t284 * t307;
t36 = t235 * t298 + t292 * t309 + t296;
t34 = t139 + (t244 * t285 + t247 * t298) * t234;
t33 = t35 + t294;
t28 = t32 + t294;
t1 = [t246 * (Icges(3,2) * t246 + t314) + t243 * (Icges(3,1) * t243 + t313) + Icges(2,3) + (t211 + t339) * t234 + (t210 + t334) * t235 + m(6) * (t83 ^ 2 + t84 ^ 2) + m(7) * (t68 ^ 2 + t69 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t149 ^ 2 + t150 ^ 2) + m(3) * (t159 ^ 2 + t160 ^ 2) + m(2) * (t219 ^ 2 + t220 ^ 2) + t332; (t246 * (Icges(3,6) * t244 + t247 * t265) + t243 * (Icges(3,5) * t244 + t247 * t267)) * t329 + m(3) * (-t159 * t247 - t160 * t244) * t218 + m(7) * (t68 * t91 + t69 * t90) + m(6) * (t83 * t99 + t84 * t98) + m(5) * (t107 * t97 + t108 * t96) + m(4) * (t149 * t179 + t150 * t178) + (t240 / 0.2e1 + t239 / 0.2e1) * (Icges(3,5) * t243 + Icges(3,6) * t246) + (t246 * (-Icges(3,6) * t247 + t244 * t265) + t243 * (-Icges(3,5) * t247 + t244 * t267)) * t328 + t249; m(7) * (t28 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(6) * (t33 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(5) * (t107 ^ 2 + t108 ^ 2 + t39 ^ 2) + m(4) * (t178 ^ 2 + t179 ^ 2 + t86 ^ 2) + t244 * t239 * t189 + m(3) * (t218 ^ 2 * t287 + t140 ^ 2) + t277 + (-t240 * t188 + (-t244 * t188 + t247 * t189) * t244 + t331) * t247; m(6) * (t100 * t84 + t101 * t83) + m(7) * (t68 * t93 + t69 * t92) + m(5) * (t141 * t97 + t142 * t96) + m(4) * (-t149 * t247 - t150 * t244) * t212 + t249; m(7) * (t32 * t28 + t90 * t92 + t91 * t93) + m(6) * (t100 * t98 + t101 * t99 + t35 * t33) + m(5) * (t107 * t141 + t108 * t142 + t70 * t39) + m(4) * (t106 * t86 + (-t178 * t244 - t179 * t247) * t212) + t252; m(7) * (t32 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(6) * (t100 ^ 2 + t101 ^ 2 + t35 ^ 2) + m(5) * (t141 ^ 2 + t142 ^ 2 + t70 ^ 2) + m(4) * (t212 ^ 2 * t287 + t106 ^ 2) + t252; -t286 * t235 + m(6) * (t66 * t83 + t67 * t84) + m(7) * (t36 * t68 + t37 * t69) + m(5) * (t94 * t96 + t95 * t97) + ((t65 / 0.2e1 + t63 / 0.2e1 + t82 / 0.2e1 + t80 / 0.2e1 + t61 / 0.2e1 + t79 / 0.2e1) * t247 + (t81 / 0.2e1 + t78 / 0.2e1 + t77 / 0.2e1 + t64 / 0.2e1 + t62 / 0.2e1 + t60 / 0.2e1) * t244) * t234; m(7) * (t34 * t28 + t36 * t91 + t37 * t90) + m(6) * (t38 * t33 + t66 * t99 + t67 * t98) + m(5) * (t107 * t95 + t108 * t94 + t39 * t85) + t250; m(7) * (t34 * t32 + t36 * t93 + t37 * t92) + m(6) * (t100 * t67 + t101 * t66 + t38 * t35) + m(5) * (t141 * t95 + t142 * t94 + t70 * t85) + t250; m(7) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t38 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t85 ^ 2 + t94 ^ 2 + t95 ^ 2) + t286 * t235 ^ 2 + (t337 * t247 + t338 * t244 + ((-t61 - t63 - t65) * t247 + (-t60 - t62 - t64) * t244) * t235) * t234; m(6) * (t203 * t84 + t205 * t83) + m(7) * (t203 * t69 + t205 * t68); m(7) * (t203 * t90 + t205 * t91 + t28 * t310) + m(6) * (t203 * t98 + t205 * t99 + t310 * t33); m(7) * (t203 * t92 + t205 * t93 + t310 * t32) + m(6) * (t203 * t100 + t205 * t101 + t310 * t35); m(7) * (t203 * t37 + t205 * t36 + t310 * t34) + m(6) * (t203 * t67 + t205 * t66 + t310 * t38); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t233 * t242 ^ 2 + t203 ^ 2 + t205 ^ 2); m(7) * (t204 * t69 + t206 * t68); m(7) * (t204 * t90 + t206 * t91 + t28 * t308); m(7) * (t204 * t92 + t206 * t93 + t308 * t32); m(7) * (t204 * t37 + t206 * t36 + t308 * t34); m(7) * (t233 * t245 * t242 + t204 * t203 + t206 * t205); m(7) * (t233 * t245 ^ 2 + t204 ^ 2 + t206 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
