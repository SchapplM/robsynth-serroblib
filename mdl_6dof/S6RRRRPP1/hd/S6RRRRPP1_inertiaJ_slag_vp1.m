% Calculate joint inertia matrix for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPP1_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_inertiaJ_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP1_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP1_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:44:14
% EndTime: 2019-03-09 20:44:26
% DurationCPUTime: 5.55s
% Computational Cost: add. (12609->481), mult. (12709->688), div. (0->0), fcn. (13704->10), ass. (0->241)
t247 = qJ(4) + pkin(10);
t241 = cos(t247);
t250 = qJ(2) + qJ(3);
t243 = cos(t250);
t257 = cos(qJ(1));
t240 = sin(t247);
t254 = sin(qJ(1));
t321 = t240 * t254;
t193 = t241 * t257 + t243 * t321;
t314 = t254 * t241;
t194 = -t240 * t257 + t243 * t314;
t242 = sin(t250);
t319 = t242 * t254;
t111 = Icges(7,5) * t194 + Icges(7,6) * t319 + Icges(7,3) * t193;
t117 = Icges(6,4) * t194 - Icges(6,2) * t193 + Icges(6,6) * t319;
t370 = t111 - t117;
t317 = t243 * t257;
t195 = t240 * t317 - t314;
t196 = t241 * t317 + t321;
t318 = t242 * t257;
t112 = Icges(7,5) * t196 + Icges(7,6) * t318 + Icges(7,3) * t195;
t118 = Icges(6,4) * t196 - Icges(6,2) * t195 + Icges(6,6) * t318;
t369 = t112 - t118;
t119 = Icges(7,1) * t194 + Icges(7,4) * t319 + Icges(7,5) * t193;
t121 = Icges(6,1) * t194 - Icges(6,4) * t193 + Icges(6,5) * t319;
t368 = t119 + t121;
t120 = Icges(7,1) * t196 + Icges(7,4) * t318 + Icges(7,5) * t195;
t122 = Icges(6,1) * t196 - Icges(6,4) * t195 + Icges(6,5) * t318;
t367 = t120 + t122;
t113 = Icges(6,5) * t194 - Icges(6,6) * t193 + Icges(6,3) * t319;
t115 = Icges(7,4) * t194 + Icges(7,2) * t319 + Icges(7,6) * t193;
t255 = cos(qJ(4));
t312 = t255 * t257;
t252 = sin(qJ(4));
t316 = t252 * t254;
t207 = -t243 * t316 - t312;
t313 = t254 * t255;
t315 = t252 * t257;
t208 = t243 * t313 - t315;
t139 = Icges(5,5) * t208 + Icges(5,6) * t207 + Icges(5,3) * t319;
t366 = t113 + t115 + t139;
t114 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t318;
t116 = Icges(7,4) * t196 + Icges(7,2) * t318 + Icges(7,6) * t195;
t209 = -t243 * t315 + t313;
t210 = t243 * t312 + t316;
t140 = Icges(5,5) * t210 + Icges(5,6) * t209 + Icges(5,3) * t318;
t365 = t114 + t116 + t140;
t356 = rSges(7,3) + qJ(6);
t357 = rSges(7,1) + pkin(5);
t364 = -t193 * t356 - t194 * t357;
t141 = Icges(5,4) * t208 + Icges(5,2) * t207 + Icges(5,6) * t319;
t143 = Icges(5,1) * t208 + Icges(5,4) * t207 + Icges(5,5) * t319;
t363 = t141 * t207 + t143 * t208 + t193 * t370 + t194 * t368 + t319 * t366;
t142 = Icges(5,4) * t210 + Icges(5,2) * t209 + Icges(5,6) * t318;
t144 = Icges(5,1) * t210 + Icges(5,4) * t209 + Icges(5,5) * t318;
t362 = t142 * t207 + t144 * t208 + t193 * t369 + t194 * t367 + t319 * t365;
t361 = t141 * t209 + t143 * t210 + t195 * t370 + t196 * t368 + t318 * t366;
t360 = t142 * t209 + t144 * t210 + t195 * t369 + t196 * t367 + t318 * t365;
t160 = -Icges(7,6) * t243 + (Icges(7,5) * t241 + Icges(7,3) * t240) * t242;
t162 = -Icges(7,2) * t243 + (Icges(7,4) * t241 + Icges(7,6) * t240) * t242;
t164 = -Icges(7,4) * t243 + (Icges(7,1) * t241 + Icges(7,5) * t240) * t242;
t75 = t160 * t193 + t162 * t319 + t164 * t194;
t161 = -Icges(6,3) * t243 + (Icges(6,5) * t241 - Icges(6,6) * t240) * t242;
t163 = -Icges(6,6) * t243 + (Icges(6,4) * t241 - Icges(6,2) * t240) * t242;
t165 = -Icges(6,5) * t243 + (Icges(6,1) * t241 - Icges(6,4) * t240) * t242;
t76 = t161 * t319 - t163 * t193 + t165 * t194;
t174 = -Icges(5,3) * t243 + (Icges(5,5) * t255 - Icges(5,6) * t252) * t242;
t175 = -Icges(5,6) * t243 + (Icges(5,4) * t255 - Icges(5,2) * t252) * t242;
t176 = -Icges(5,5) * t243 + (Icges(5,1) * t255 - Icges(5,4) * t252) * t242;
t81 = t174 * t319 + t175 * t207 + t176 * t208;
t359 = -t81 - t75 - t76;
t77 = t160 * t195 + t162 * t318 + t164 * t196;
t78 = t161 * t318 - t163 * t195 + t165 * t196;
t82 = t174 * t318 + t175 * t209 + t176 * t210;
t358 = -t82 - t77 - t78;
t310 = rSges(7,2) * t319 - t364;
t355 = rSges(7,2) * t318 + t195 * t356 + t196 * t357;
t354 = t359 * t243 + (t254 * t363 + t257 * t362) * t242;
t353 = t358 * t243 + (t254 * t361 + t257 * t360) * t242;
t352 = t254 * t362 - t257 * t363;
t351 = t254 * t360 - t257 * t361;
t350 = t161 + t162 + t174;
t322 = t240 * t242;
t348 = t160 * t322 + (t255 * t176 + (t164 + t165) * t241) * t242;
t272 = Icges(4,5) * t243 - Icges(4,6) * t242;
t180 = -Icges(4,3) * t257 + t254 * t272;
t181 = Icges(4,3) * t254 + t257 * t272;
t249 = t257 ^ 2;
t324 = Icges(4,4) * t243;
t274 = -Icges(4,2) * t242 + t324;
t183 = Icges(4,6) * t254 + t257 * t274;
t325 = Icges(4,4) * t242;
t276 = Icges(4,1) * t243 - t325;
t185 = Icges(4,5) * t254 + t257 * t276;
t270 = -t183 * t242 + t185 * t243;
t182 = -Icges(4,6) * t257 + t254 * t274;
t184 = -Icges(4,5) * t257 + t254 * t276;
t271 = t182 * t242 - t184 * t243;
t347 = -t249 * t180 - (t270 * t254 + (-t181 + t271) * t257) * t254 - t352;
t346 = t243 ^ 2;
t248 = t254 ^ 2;
t345 = m(6) / 0.2e1;
t344 = m(7) / 0.2e1;
t342 = t254 / 0.2e1;
t341 = -t257 / 0.2e1;
t253 = sin(qJ(2));
t340 = pkin(2) * t253;
t339 = pkin(3) * t243;
t338 = pkin(9) * t242;
t237 = pkin(4) * t255 + pkin(3);
t337 = -pkin(3) + t237;
t256 = cos(qJ(2));
t336 = rSges(3,1) * t256;
t335 = rSges(3,2) * t253;
t334 = t257 * rSges(3,3);
t52 = -t243 * t115 + (t111 * t240 + t119 * t241) * t242;
t333 = t52 * t257;
t53 = -t243 * t116 + (t112 * t240 + t120 * t241) * t242;
t332 = t53 * t254;
t54 = -t243 * t113 + (-t117 * t240 + t121 * t241) * t242;
t331 = t54 * t257;
t55 = -t243 * t114 + (-t118 * t240 + t122 * t241) * t242;
t330 = t55 * t254;
t66 = -t243 * t139 + (-t141 * t252 + t143 * t255) * t242;
t329 = t66 * t257;
t67 = -t243 * t140 + (-t142 * t252 + t144 * t255) * t242;
t328 = t67 * t254;
t327 = Icges(3,4) * t253;
t326 = Icges(3,4) * t256;
t323 = t175 * t252;
t258 = -pkin(8) - pkin(7);
t311 = t257 * t258;
t126 = rSges(6,1) * t196 - rSges(6,2) * t195 + rSges(6,3) * t318;
t251 = -qJ(5) - pkin(9);
t264 = pkin(4) * t316 + t237 * t317 - t251 * t318;
t300 = pkin(3) * t317 + pkin(9) * t318;
t138 = t264 - t300;
t309 = -t126 - t138;
t301 = -pkin(4) * t315 - t251 * t319;
t137 = (t243 * t337 - t338) * t254 + t301;
t159 = (pkin(9) + t251) * t243 + t337 * t242;
t308 = t137 * t243 + t159 * t319;
t167 = -t243 * rSges(6,3) + (rSges(6,1) * t241 - rSges(6,2) * t240) * t242;
t306 = -t159 - t167;
t305 = -t243 * rSges(7,2) + (t240 * t356 + t241 * t357) * t242;
t238 = pkin(2) * t256 + pkin(1);
t226 = t257 * t238;
t246 = t257 * pkin(7);
t304 = t254 * (t311 + t246 + (-pkin(1) + t238) * t254) + t257 * (-t257 * pkin(1) + t226 + (-pkin(7) - t258) * t254);
t265 = rSges(4,1) * t317 - rSges(4,2) * t318 + rSges(4,3) * t254;
t280 = rSges(4,1) * t243 - rSges(4,2) * t242;
t128 = t254 * (-t257 * rSges(4,3) + t254 * t280) + t257 * t265;
t177 = -t243 * rSges(5,3) + (rSges(5,1) * t255 - rSges(5,2) * t252) * t242;
t217 = t242 * pkin(3) - t243 * pkin(9);
t303 = -t177 - t217;
t302 = t248 * (t338 + t339) + t257 * t300;
t299 = rSges(3,3) * t254 + t257 * t336;
t298 = t248 + t249;
t297 = t163 * t322 + t242 * t323 + t243 * t350 - t348;
t296 = -t138 - t355;
t295 = -t159 - t305;
t294 = -t217 + t306;
t148 = rSges(5,1) * t210 + rSges(5,2) * t209 + rSges(5,3) * t318;
t216 = rSges(4,1) * t242 + rSges(4,2) * t243;
t291 = -t216 - t340;
t290 = -t217 - t340;
t289 = (t248 * t181 + (t271 * t257 + (-t180 + t270) * t254) * t257 + t351) * t254;
t288 = -t254 * t258 + t226;
t287 = -t237 * t243 - t238;
t286 = t137 * t254 + t138 * t257 + t302;
t285 = -t217 + t295;
t279 = -t208 * rSges(5,1) - t207 * rSges(5,2);
t147 = rSges(5,3) * t319 - t279;
t74 = t147 * t254 + t148 * t257 + t302;
t284 = -t159 + t290;
t283 = -t177 + t290;
t282 = -t301 - t311;
t281 = -t335 + t336;
t278 = -t194 * rSges(6,1) + t193 * rSges(6,2);
t277 = Icges(3,1) * t256 - t327;
t275 = -Icges(3,2) * t253 + t326;
t273 = Icges(3,5) * t256 - Icges(3,6) * t253;
t213 = Icges(4,2) * t243 + t325;
t214 = Icges(4,1) * t242 + t324;
t267 = -t213 * t242 + t214 * t243;
t266 = -t167 + t284;
t124 = rSges(6,3) * t319 - t278;
t37 = t124 * t254 + t126 * t257 + t286;
t263 = t284 - t305;
t31 = t254 * t310 + t257 * t355 + t286;
t262 = t257 * t347 + t289;
t261 = t264 + t288;
t260 = -(t332 - t333 + t330 - t331 + t328 - t329) * t243 / 0.2e1 + t353 * t342 + t354 * t341 + t352 * t319 / 0.2e1 + t351 * t318 / 0.2e1;
t212 = Icges(4,5) * t242 + Icges(4,6) * t243;
t259 = -t333 / 0.2e1 + t332 / 0.2e1 - t331 / 0.2e1 + t330 / 0.2e1 - t329 / 0.2e1 + t328 / 0.2e1 + (t183 * t243 + t185 * t242 + t254 * t212 + t257 * t267 - t358) * t342 + (t182 * t243 + t184 * t242 - t257 * t212 + t254 * t267 - t359) * t341;
t239 = t242 ^ 2;
t225 = rSges(2,1) * t257 - rSges(2,2) * t254;
t224 = -rSges(2,1) * t254 - rSges(2,2) * t257;
t223 = rSges(3,1) * t253 + rSges(3,2) * t256;
t200 = Icges(3,3) * t254 + t257 * t273;
t199 = -Icges(3,3) * t257 + t254 * t273;
t179 = t291 * t257;
t178 = t291 * t254;
t169 = t254 * pkin(7) + (pkin(1) - t335) * t257 + t299;
t168 = t334 + t246 + (-pkin(1) - t281) * t254;
t157 = t265 + t288;
t156 = (rSges(4,3) - t258) * t257 + (-t238 - t280) * t254;
t151 = t303 * t257;
t150 = t303 * t254;
t149 = t257 * (-t257 * t335 + t299) + (t254 * t281 - t334) * t254;
t132 = t283 * t257;
t131 = t283 * t254;
t110 = t137 * t318;
t101 = t288 + t148 + t300;
t100 = -t311 + (-t339 - t238 + (-rSges(5,3) - pkin(9)) * t242) * t254 + t279;
t99 = t294 * t257;
t98 = t294 * t254;
t97 = -t148 * t243 - t177 * t318;
t96 = t147 * t243 + t177 * t319;
t95 = t266 * t257;
t94 = t266 * t254;
t92 = t261 + t126;
t91 = (-rSges(6,3) * t242 + t287) * t254 + t278 + t282;
t90 = t128 + t304;
t89 = (t147 * t257 - t148 * t254) * t242;
t88 = t285 * t257;
t87 = t285 * t254;
t84 = t263 * t257;
t83 = t263 * t254;
t69 = t261 + t355;
t68 = (-rSges(7,2) * t242 + t287) * t254 + t282 + t364;
t57 = t243 * t309 + t306 * t318;
t56 = t124 * t243 + t167 * t319 + t308;
t51 = t74 + t304;
t38 = t110 + (t124 * t257 + t254 * t309) * t242;
t36 = t243 * t296 + t295 * t318;
t35 = t243 * t310 + t305 * t319 + t308;
t34 = t37 + t304;
t32 = t110 + (t254 * t296 + t257 * t310) * t242;
t23 = t31 + t304;
t1 = [t256 * (Icges(3,2) * t256 + t327) + t253 * (Icges(3,1) * t253 + t326) + Icges(2,3) + (-t163 * t240 + t214 - t323) * t242 + (t213 - t350) * t243 + m(7) * (t68 ^ 2 + t69 ^ 2) + m(6) * (t91 ^ 2 + t92 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2) + m(4) * (t156 ^ 2 + t157 ^ 2) + m(3) * (t168 ^ 2 + t169 ^ 2) + m(2) * (t224 ^ 2 + t225 ^ 2) + t348; m(3) * (-t168 * t257 - t169 * t254) * t223 + (t249 / 0.2e1 + t248 / 0.2e1) * (Icges(3,5) * t253 + Icges(3,6) * t256) + m(5) * (t100 * t132 + t101 * t131) + m(6) * (t91 * t95 + t92 * t94) + m(7) * (t68 * t84 + t69 * t83) + m(4) * (t156 * t179 + t157 * t178) + (t256 * (Icges(3,6) * t254 + t257 * t275) + t253 * (Icges(3,5) * t254 + t257 * t277)) * t342 + (t256 * (-Icges(3,6) * t257 + t254 * t275) + t253 * (-Icges(3,5) * t257 + t254 * t277)) * t341 + t259; m(7) * (t23 ^ 2 + t83 ^ 2 + t84 ^ 2) + m(6) * (t34 ^ 2 + t94 ^ 2 + t95 ^ 2) + m(5) * (t131 ^ 2 + t132 ^ 2 + t51 ^ 2) + m(4) * (t178 ^ 2 + t179 ^ 2 + t90 ^ 2) + t254 * t248 * t200 + m(3) * (t223 ^ 2 * t298 + t149 ^ 2) + t289 + (-t249 * t199 + (-t254 * t199 + t257 * t200) * t254 + t347) * t257; m(4) * (-t156 * t257 - t157 * t254) * t216 + m(7) * (t68 * t88 + t69 * t87) + m(6) * (t91 * t99 + t92 * t98) + m(5) * (t100 * t151 + t101 * t150) + t259; m(7) * (t23 * t31 + t83 * t87 + t84 * t88) + m(6) * (t34 * t37 + t94 * t98 + t95 * t99) + m(5) * (t131 * t150 + t132 * t151 + t51 * t74) + m(4) * (t128 * t90 + (-t178 * t254 - t179 * t257) * t216) + t262; m(7) * (t31 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(6) * (t37 ^ 2 + t98 ^ 2 + t99 ^ 2) + m(5) * (t150 ^ 2 + t151 ^ 2 + t74 ^ 2) + m(4) * (t216 ^ 2 * t298 + t128 ^ 2) + t262; t297 * t243 + m(7) * (t35 * t68 + t36 * t69) + m(6) * (t56 * t91 + t57 * t92) + m(5) * (t100 * t96 + t101 * t97) + ((t53 / 0.2e1 + t67 / 0.2e1 + t82 / 0.2e1 + t78 / 0.2e1 + t77 / 0.2e1 + t55 / 0.2e1) * t257 + (t81 / 0.2e1 + t76 / 0.2e1 + t75 / 0.2e1 + t54 / 0.2e1 + t52 / 0.2e1 + t66 / 0.2e1) * t254) * t242; m(7) * (t23 * t32 + t35 * t84 + t36 * t83) + m(6) * (t34 * t38 + t56 * t95 + t57 * t94) + m(5) * (t131 * t97 + t132 * t96 + t51 * t89) + t260; m(7) * (t31 * t32 + t35 * t88 + t36 * t87) + m(6) * (t37 * t38 + t56 * t99 + t57 * t98) + m(5) * (t150 * t97 + t151 * t96 + t74 * t89) + t260; m(6) * (t38 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(7) * (t32 ^ 2 + t35 ^ 2 + t36 ^ 2) + m(5) * (t89 ^ 2 + t96 ^ 2 + t97 ^ 2) - t297 * t346 + (t353 * t257 + t354 * t254 + ((-t53 - t55 - t67) * t257 + (-t52 - t54 - t66) * t254) * t243) * t242; 0.2e1 * ((t254 * t69 + t257 * t68) * t344 + (t254 * t92 + t257 * t91) * t345) * t242; m(7) * (-t243 * t23 + (t254 * t83 + t257 * t84) * t242) + m(6) * (-t243 * t34 + (t254 * t94 + t257 * t95) * t242); m(7) * (-t243 * t31 + (t254 * t87 + t257 * t88) * t242) + m(6) * (-t243 * t37 + (t254 * t98 + t257 * t99) * t242); m(6) * (-t243 * t38 + (t254 * t57 + t257 * t56) * t242) + m(7) * (-t243 * t32 + (t254 * t36 + t257 * t35) * t242); 0.2e1 * (t345 + t344) * (t239 * t298 + t346); m(7) * (t193 * t69 + t195 * t68); m(7) * (t193 * t83 + t195 * t84 + t23 * t322); m(7) * (t193 * t87 + t195 * t88 + t31 * t322); m(7) * (t193 * t36 + t195 * t35 + t32 * t322); m(7) * (t193 * t254 + t195 * t257 - t240 * t243) * t242; m(7) * (t239 * t240 ^ 2 + t193 ^ 2 + t195 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
