% Calculate joint inertia matrix for
% S6RRRRPP2
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:49:03
% EndTime: 2019-03-09 20:49:15
% DurationCPUTime: 5.19s
% Computational Cost: add. (9366->467), mult. (12404->664), div. (0->0), fcn. (13552->8), ass. (0->228)
t241 = qJ(2) + qJ(3);
t235 = cos(t241);
t245 = cos(qJ(4));
t247 = cos(qJ(1));
t301 = t245 * t247;
t242 = sin(qJ(4));
t244 = sin(qJ(1));
t304 = t242 * t244;
t206 = t235 * t304 + t301;
t302 = t244 * t245;
t303 = t242 * t247;
t207 = t235 * t302 - t303;
t234 = sin(t241);
t309 = t234 * t244;
t115 = Icges(7,5) * t207 + Icges(7,6) * t206 - Icges(7,3) * t309;
t119 = Icges(5,5) * t207 - Icges(5,6) * t206 + Icges(5,3) * t309;
t123 = Icges(6,4) * t207 + Icges(6,2) * t309 + Icges(6,6) * t206;
t356 = -t115 + t119 + t123;
t208 = t235 * t303 - t302;
t209 = t235 * t301 + t304;
t307 = t234 * t247;
t116 = Icges(7,5) * t209 + Icges(7,6) * t208 - Icges(7,3) * t307;
t120 = Icges(5,5) * t209 - Icges(5,6) * t208 + Icges(5,3) * t307;
t124 = Icges(6,4) * t209 + Icges(6,2) * t307 + Icges(6,6) * t208;
t355 = -t116 + t120 + t124;
t117 = Icges(6,5) * t207 + Icges(6,6) * t309 + Icges(6,3) * t206;
t121 = Icges(7,4) * t207 + Icges(7,2) * t206 - Icges(7,6) * t309;
t125 = Icges(5,4) * t207 - Icges(5,2) * t206 + Icges(5,6) * t309;
t354 = t117 + t121 - t125;
t118 = Icges(6,5) * t209 + Icges(6,6) * t307 + Icges(6,3) * t208;
t122 = Icges(7,4) * t209 + Icges(7,2) * t208 - Icges(7,6) * t307;
t126 = Icges(5,4) * t209 - Icges(5,2) * t208 + Icges(5,6) * t307;
t353 = t118 + t122 - t126;
t127 = Icges(7,1) * t207 + Icges(7,4) * t206 - Icges(7,5) * t309;
t129 = Icges(6,1) * t207 + Icges(6,4) * t309 + Icges(6,5) * t206;
t131 = Icges(5,1) * t207 - Icges(5,4) * t206 + Icges(5,5) * t309;
t352 = t127 + t129 + t131;
t128 = Icges(7,1) * t209 + Icges(7,4) * t208 - Icges(7,5) * t307;
t130 = Icges(6,1) * t209 + Icges(6,4) * t307 + Icges(6,5) * t208;
t132 = Icges(5,1) * t209 - Icges(5,4) * t208 + Icges(5,5) * t307;
t351 = t128 + t130 + t132;
t350 = rSges(7,1) + pkin(5);
t342 = rSges(7,3) + qJ(6);
t341 = t208 * rSges(7,2) + t209 * t350 - t307 * t342;
t349 = t354 * t206 + t352 * t207 + t309 * t356;
t348 = t206 * t353 + t207 * t351 + t309 * t355;
t347 = t354 * t208 + t352 * t209 + t307 * t356;
t346 = t208 * t353 + t209 * t351 + t307 * t355;
t168 = Icges(7,3) * t235 + (Icges(7,5) * t245 + Icges(7,6) * t242) * t234;
t171 = Icges(7,6) * t235 + (Icges(7,4) * t245 + Icges(7,2) * t242) * t234;
t174 = Icges(7,5) * t235 + (Icges(7,1) * t245 + Icges(7,4) * t242) * t234;
t77 = -t168 * t309 + t171 * t206 + t174 * t207;
t169 = -Icges(6,6) * t235 + (Icges(6,5) * t245 + Icges(6,3) * t242) * t234;
t172 = -Icges(6,2) * t235 + (Icges(6,4) * t245 + Icges(6,6) * t242) * t234;
t175 = -Icges(6,4) * t235 + (Icges(6,1) * t245 + Icges(6,5) * t242) * t234;
t78 = t169 * t206 + t172 * t309 + t175 * t207;
t170 = -Icges(5,3) * t235 + (Icges(5,5) * t245 - Icges(5,6) * t242) * t234;
t173 = -Icges(5,6) * t235 + (Icges(5,4) * t245 - Icges(5,2) * t242) * t234;
t176 = -Icges(5,5) * t235 + (Icges(5,1) * t245 - Icges(5,4) * t242) * t234;
t79 = t170 * t309 - t173 * t206 + t176 * t207;
t345 = -t77 - t78 - t79;
t80 = -t168 * t307 + t171 * t208 + t174 * t209;
t81 = t169 * t208 + t172 * t307 + t175 * t209;
t82 = t170 * t307 - t173 * t208 + t176 * t209;
t344 = -t80 - t81 - t82;
t324 = t206 * rSges(7,2);
t299 = t350 * t207 - t342 * t309 + t324;
t340 = t170 + t172;
t339 = t345 * t235 + (t244 * t349 + t348 * t247) * t234;
t338 = t344 * t235 + (t244 * t347 + t247 * t346) * t234;
t337 = t348 * t244 - t247 * t349;
t336 = t244 * t346 - t247 * t347;
t308 = t234 * t245;
t310 = t234 * t242;
t335 = t235 * t168 + (t169 + t171) * t310 + (t174 + t175 + t176) * t308;
t262 = Icges(4,5) * t235 - Icges(4,6) * t234;
t182 = -Icges(4,3) * t247 + t244 * t262;
t183 = Icges(4,3) * t244 + t247 * t262;
t240 = t247 ^ 2;
t311 = Icges(4,4) * t235;
t264 = -Icges(4,2) * t234 + t311;
t185 = Icges(4,6) * t244 + t247 * t264;
t312 = Icges(4,4) * t234;
t266 = Icges(4,1) * t235 - t312;
t187 = Icges(4,5) * t244 + t247 * t266;
t260 = -t185 * t234 + t187 * t235;
t184 = -Icges(4,6) * t247 + t244 * t264;
t186 = -Icges(4,5) * t247 + t244 * t266;
t261 = t184 * t234 - t186 * t235;
t334 = -t240 * t182 - (t260 * t244 + (-t183 + t261) * t247) * t244 - t337;
t333 = t235 ^ 2;
t239 = t244 ^ 2;
t331 = t244 / 0.2e1;
t330 = -t247 / 0.2e1;
t329 = m(7) * t234;
t243 = sin(qJ(2));
t328 = pkin(2) * t243;
t327 = pkin(3) * t235;
t246 = cos(qJ(2));
t326 = rSges(3,1) * t246;
t325 = rSges(3,2) * t243;
t323 = t206 * rSges(6,3);
t322 = t247 * rSges(3,3);
t60 = t115 * t235 + (t121 * t242 + t127 * t245) * t234;
t321 = t60 * t247;
t61 = t116 * t235 + (t122 * t242 + t128 * t245) * t234;
t320 = t61 * t244;
t62 = -t123 * t235 + (t117 * t242 + t129 * t245) * t234;
t319 = t62 * t247;
t63 = -t124 * t235 + (t118 * t242 + t130 * t245) * t234;
t318 = t63 * t244;
t64 = -t119 * t235 + (-t125 * t242 + t131 * t245) * t234;
t317 = t64 * t247;
t65 = -t120 * t235 + (-t126 * t242 + t132 * t245) * t234;
t316 = t65 * t244;
t314 = Icges(3,4) * t243;
t313 = Icges(3,4) * t246;
t306 = t235 * t247;
t305 = t242 * t173;
t248 = -pkin(8) - pkin(7);
t300 = t247 * t248;
t137 = t209 * rSges(6,1) + rSges(6,2) * t307 + t208 * rSges(6,3);
t148 = t209 * pkin(4) + t208 * qJ(5);
t298 = -t137 - t148;
t196 = t206 * qJ(5);
t147 = t207 * pkin(4) + t196;
t205 = (pkin(4) * t245 + qJ(5) * t242) * t234;
t297 = t235 * t147 + t205 * t309;
t232 = pkin(2) * t246 + pkin(1);
t224 = t247 * t232;
t238 = t247 * pkin(7);
t295 = t244 * (t300 + t238 + (-pkin(1) + t232) * t244) + t247 * (-t247 * pkin(1) + t224 + (-pkin(7) - t248) * t244);
t255 = rSges(4,1) * t306 - rSges(4,2) * t307 + t244 * rSges(4,3);
t269 = rSges(4,1) * t235 - rSges(4,2) * t234;
t106 = t244 * (-t247 * rSges(4,3) + t244 * t269) + t247 * t255;
t294 = (rSges(7,1) * t245 + rSges(7,2) * t242) * t234 + pkin(5) * t308 + t342 * t235;
t178 = -rSges(6,2) * t235 + (rSges(6,1) * t245 + rSges(6,3) * t242) * t234;
t293 = -t178 - t205;
t179 = -rSges(5,3) * t235 + (rSges(5,1) * t245 - rSges(5,2) * t242) * t234;
t216 = pkin(3) * t234 - pkin(9) * t235;
t292 = -t179 - t216;
t289 = pkin(3) * t306 + pkin(9) * t307;
t291 = t239 * (pkin(9) * t234 + t327) + t247 * t289;
t288 = t244 * rSges(3,3) + t247 * t326;
t287 = t239 + t240;
t286 = t234 * t305 + t340 * t235 - t335;
t285 = -t148 - t341;
t284 = -t205 - t294;
t283 = -t216 + t293;
t138 = t209 * rSges(5,1) - t208 * rSges(5,2) + rSges(5,3) * t307;
t215 = rSges(4,1) * t234 + rSges(4,2) * t235;
t280 = -t215 - t328;
t279 = -t216 - t328;
t278 = -t232 - t327;
t277 = (t239 * t183 + (t261 * t247 + (-t182 + t260) * t244) * t247 + t336) * t244;
t276 = -t196 - t300;
t275 = -t244 * t248 + t224;
t274 = t244 * t147 + t247 * t148 + t291;
t273 = -t216 + t284;
t268 = -t207 * rSges(5,1) + t206 * rSges(5,2);
t135 = rSges(5,3) * t309 - t268;
t70 = t244 * t135 + t247 * t138 + t291;
t272 = -t179 + t279;
t271 = -t205 + t279;
t270 = -t325 + t326;
t267 = Icges(3,1) * t246 - t314;
t265 = -Icges(3,2) * t243 + t313;
t263 = Icges(3,5) * t246 - Icges(3,6) * t243;
t213 = Icges(4,2) * t235 + t312;
t214 = Icges(4,1) * t234 + t311;
t257 = -t213 * t234 + t214 * t235;
t256 = -t178 + t271;
t254 = t275 + t289;
t134 = t207 * rSges(6,1) + rSges(6,2) * t309 + t323;
t35 = t244 * t134 + t247 * t137 + t274;
t253 = t271 - t294;
t32 = t244 * t299 + t247 * t341 + t274;
t252 = t334 * t247 + t277;
t251 = t148 + t254;
t250 = -(t320 - t321 + t318 - t319 + t316 - t317) * t235 / 0.2e1 + t338 * t331 + t339 * t330 + t337 * t309 / 0.2e1 + t336 * t307 / 0.2e1;
t212 = Icges(4,5) * t234 + Icges(4,6) * t235;
t249 = -t321 / 0.2e1 + t320 / 0.2e1 - t319 / 0.2e1 + t318 / 0.2e1 - t317 / 0.2e1 + t316 / 0.2e1 + (t185 * t235 + t187 * t234 + t212 * t244 + t247 * t257 - t344) * t331 + (t184 * t235 + t186 * t234 - t212 * t247 + t244 * t257 - t345) * t330;
t233 = t234 ^ 2;
t223 = rSges(2,1) * t247 - rSges(2,2) * t244;
t222 = -rSges(2,1) * t244 - rSges(2,2) * t247;
t221 = rSges(3,1) * t243 + rSges(3,2) * t246;
t191 = Icges(3,3) * t244 + t247 * t263;
t190 = -Icges(3,3) * t247 + t244 * t263;
t181 = t280 * t247;
t180 = t280 * t244;
t162 = t244 * pkin(7) + (pkin(1) - t325) * t247 + t288;
t161 = t322 + t238 + (-pkin(1) - t270) * t244;
t150 = t255 + t275;
t149 = (rSges(4,3) - t248) * t247 + (-t232 - t269) * t244;
t142 = t292 * t247;
t141 = t292 * t244;
t140 = t247 * (-t247 * t325 + t288) + (t244 * t270 - t322) * t244;
t139 = t147 * t307;
t108 = t272 * t247;
t107 = t272 * t244;
t101 = t283 * t247;
t100 = t283 * t244;
t99 = t256 * t247;
t98 = t256 * t244;
t97 = t254 + t138;
t96 = -t300 + ((-rSges(5,3) - pkin(9)) * t234 + t278) * t244 + t268;
t95 = -t138 * t235 - t179 * t307;
t94 = t135 * t235 + t179 * t309;
t93 = t273 * t247;
t92 = t273 * t244;
t91 = t253 * t247;
t90 = t253 * t244;
t86 = t106 + t295;
t85 = (t135 * t247 - t138 * t244) * t234;
t84 = t251 + t137;
t83 = -t323 + (-rSges(6,1) - pkin(4)) * t207 + ((-rSges(6,2) - pkin(9)) * t234 + t278) * t244 + t276;
t69 = t251 + t341;
t68 = -t324 + (-pkin(4) - t350) * t207 + ((-pkin(9) + t342) * t234 + t278) * t244 + t276;
t67 = t235 * t298 + t293 * t307;
t66 = t134 * t235 + t178 * t309 + t297;
t39 = t70 + t295;
t38 = t139 + (t134 * t247 + t244 * t298) * t234;
t37 = t235 * t285 + t284 * t307;
t36 = t235 * t299 + t294 * t309 + t297;
t34 = t139 + (t244 * t285 + t247 * t299) * t234;
t33 = t35 + t295;
t28 = t32 + t295;
t1 = [t246 * (Icges(3,2) * t246 + t314) + t243 * (Icges(3,1) * t243 + t313) + Icges(2,3) + (t214 - t305) * t234 + (t213 - t340) * t235 + m(7) * (t68 ^ 2 + t69 ^ 2) + m(6) * (t83 ^ 2 + t84 ^ 2) + m(5) * (t96 ^ 2 + t97 ^ 2) + m(4) * (t149 ^ 2 + t150 ^ 2) + m(3) * (t161 ^ 2 + t162 ^ 2) + m(2) * (t222 ^ 2 + t223 ^ 2) + t335; t249 + m(7) * (t68 * t91 + t69 * t90) + m(6) * (t83 * t99 + t84 * t98) + m(5) * (t107 * t97 + t108 * t96) + m(4) * (t149 * t181 + t150 * t180) + (t240 / 0.2e1 + t239 / 0.2e1) * (Icges(3,5) * t243 + Icges(3,6) * t246) + m(3) * (-t161 * t247 - t162 * t244) * t221 + ((-Icges(3,6) * t247 + t244 * t265) * t246 + (-Icges(3,5) * t247 + t244 * t267) * t243) * t330 + ((Icges(3,6) * t244 + t247 * t265) * t246 + (Icges(3,5) * t244 + t247 * t267) * t243) * t331; m(6) * (t33 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(7) * (t28 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(5) * (t107 ^ 2 + t108 ^ 2 + t39 ^ 2) + m(4) * (t180 ^ 2 + t181 ^ 2 + t86 ^ 2) + m(3) * (t221 ^ 2 * t287 + t140 ^ 2) + t244 * t191 * t239 + t277 + (-t190 * t240 + (-t190 * t244 + t191 * t247) * t244 + t334) * t247; t249 + m(7) * (t68 * t93 + t69 * t92) + m(6) * (t100 * t84 + t101 * t83) + m(5) * (t141 * t97 + t142 * t96) + m(4) * (-t149 * t247 - t150 * t244) * t215; m(6) * (t100 * t98 + t101 * t99 + t35 * t33) + m(7) * (t32 * t28 + t90 * t92 + t91 * t93) + m(5) * (t107 * t141 + t108 * t142 + t70 * t39) + m(4) * (t106 * t86 + (-t180 * t244 - t181 * t247) * t215) + t252; m(7) * (t32 ^ 2 + t92 ^ 2 + t93 ^ 2) + m(6) * (t100 ^ 2 + t101 ^ 2 + t35 ^ 2) + m(5) * (t141 ^ 2 + t142 ^ 2 + t70 ^ 2) + m(4) * (t215 ^ 2 * t287 + t106 ^ 2) + t252; t286 * t235 + m(7) * (t36 * t68 + t37 * t69) + m(6) * (t66 * t83 + t67 * t84) + m(5) * (t94 * t96 + t95 * t97) + ((t82 / 0.2e1 + t81 / 0.2e1 + t80 / 0.2e1 + t63 / 0.2e1 + t61 / 0.2e1 + t65 / 0.2e1) * t247 + (t78 / 0.2e1 + t77 / 0.2e1 + t62 / 0.2e1 + t60 / 0.2e1 + t64 / 0.2e1 + t79 / 0.2e1) * t244) * t234; t250 + m(6) * (t38 * t33 + t66 * t99 + t67 * t98) + m(7) * (t34 * t28 + t36 * t91 + t37 * t90) + m(5) * (t107 * t95 + t108 * t94 + t39 * t85); t250 + m(7) * (t34 * t32 + t36 * t93 + t37 * t92) + m(6) * (t100 * t67 + t101 * t66 + t38 * t35) + m(5) * (t141 * t95 + t142 * t94 + t70 * t85); m(7) * (t34 ^ 2 + t36 ^ 2 + t37 ^ 2) + m(6) * (t38 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t85 ^ 2 + t94 ^ 2 + t95 ^ 2) - t286 * t333 + (t338 * t247 + t339 * t244 + ((-t61 - t63 - t65) * t247 + (-t60 - t62 - t64) * t244) * t235) * t234; m(7) * (t206 * t69 + t208 * t68) + m(6) * (t206 * t84 + t208 * t83); m(6) * (t206 * t98 + t208 * t99 + t310 * t33) + m(7) * (t206 * t90 + t208 * t91 + t28 * t310); m(7) * (t206 * t92 + t208 * t93 + t310 * t32) + m(6) * (t100 * t206 + t101 * t208 + t310 * t35); m(7) * (t206 * t37 + t208 * t36 + t310 * t34) + m(6) * (t206 * t67 + t208 * t66 + t310 * t38); 0.2e1 * (m(6) / 0.2e1 + m(7) / 0.2e1) * (t233 * t242 ^ 2 + t206 ^ 2 + t208 ^ 2); (-t244 * t69 - t247 * t68) * t329; m(7) * (t235 * t28 + (-t244 * t90 - t247 * t91) * t234); m(7) * (t235 * t32 + (-t244 * t92 - t247 * t93) * t234); m(7) * (t235 * t34 + (-t244 * t37 - t247 * t36) * t234); (-t206 * t244 - t208 * t247 + t235 * t242) * t329; m(7) * (t233 * t287 + t333);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
