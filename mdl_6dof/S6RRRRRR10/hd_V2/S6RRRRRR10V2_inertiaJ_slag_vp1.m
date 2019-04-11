% Calculate joint inertia matrix for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
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
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR10V2_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10V2_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:41:42
% EndTime: 2019-04-11 14:42:01
% DurationCPUTime: 5.42s
% Computational Cost: add. (28419->578), mult. (40855->838), div. (0->0), fcn. (51308->12), ass. (0->278)
t280 = qJ(2) + qJ(3);
t275 = cos(t280);
t287 = cos(qJ(4));
t289 = cos(qJ(1));
t337 = t287 * t289;
t283 = sin(qJ(4));
t285 = sin(qJ(1));
t339 = t285 * t283;
t247 = t275 * t339 + t337;
t338 = t285 * t287;
t340 = t283 * t289;
t249 = t275 * t340 - t338;
t274 = sin(t280);
t345 = t274 * t283;
t248 = t275 * t338 - t340;
t282 = sin(qJ(5));
t363 = cos(qJ(5));
t325 = t274 * t363;
t210 = t248 * t282 - t285 * t325;
t344 = t274 * t285;
t211 = t248 * t363 + t282 * t344;
t153 = Icges(6,5) * t211 - Icges(6,6) * t210 + Icges(6,3) * t247;
t155 = Icges(6,4) * t211 - Icges(6,2) * t210 + Icges(6,6) * t247;
t157 = Icges(6,1) * t211 - Icges(6,4) * t210 + Icges(6,5) * t247;
t81 = t153 * t247 - t155 * t210 + t157 * t211;
t250 = t275 * t337 + t339;
t212 = t250 * t282 - t289 * t325;
t342 = t274 * t289;
t213 = t250 * t363 + t282 * t342;
t154 = Icges(6,5) * t213 - Icges(6,6) * t212 + Icges(6,3) * t249;
t156 = Icges(6,4) * t213 - Icges(6,2) * t212 + Icges(6,6) * t249;
t158 = Icges(6,1) * t213 - Icges(6,4) * t212 + Icges(6,5) * t249;
t82 = t154 * t247 - t156 * t210 + t158 * t211;
t343 = t274 * t287;
t240 = t275 * t363 + t282 * t343;
t241 = -t275 * t282 + t287 * t325;
t185 = Icges(6,5) * t241 - Icges(6,6) * t240 + Icges(6,3) * t345;
t186 = Icges(6,4) * t241 - Icges(6,2) * t240 + Icges(6,6) * t345;
t187 = Icges(6,1) * t241 - Icges(6,4) * t240 + Icges(6,5) * t345;
t96 = t185 * t247 - t186 * t210 + t187 * t211;
t31 = t247 * t81 + t249 * t82 + t345 * t96;
t281 = sin(qJ(6));
t286 = cos(qJ(6));
t178 = -t211 * t281 + t247 * t286;
t179 = t211 * t286 + t247 * t281;
t122 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t210;
t124 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t210;
t126 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t210;
t49 = t122 * t210 + t124 * t178 + t126 * t179;
t180 = -t213 * t281 + t249 * t286;
t181 = t213 * t286 + t249 * t281;
t123 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t212;
t125 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t212;
t127 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t212;
t50 = t123 * t210 + t125 * t178 + t127 * t179;
t202 = -t241 * t281 + t286 * t345;
t203 = t241 * t286 + t281 * t345;
t146 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t240;
t147 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t240;
t148 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t240;
t65 = t146 * t210 + t147 * t178 + t148 * t179;
t7 = t247 * t49 + t249 * t50 + t345 * t65;
t384 = t7 + t31;
t83 = t153 * t249 - t155 * t212 + t157 * t213;
t84 = t154 * t249 - t156 * t212 + t158 * t213;
t97 = t185 * t249 - t186 * t212 + t187 * t213;
t32 = t247 * t83 + t249 * t84 + t345 * t97;
t51 = t122 * t212 + t124 * t180 + t126 * t181;
t52 = t123 * t212 + t125 * t180 + t127 * t181;
t66 = t146 * t212 + t147 * t180 + t148 * t181;
t8 = t247 * t51 + t249 * t52 + t345 * t66;
t383 = t8 + t32;
t23 = t50 * t285 - t289 * t49;
t382 = t82 * t285 - t289 * t81 + t23;
t24 = t52 * t285 - t289 * t51;
t381 = t84 * t285 - t289 * t83 + t24;
t56 = t123 * t240 + t125 * t202 + t127 * t203;
t354 = t56 * t285;
t55 = t122 * t240 + t124 * t202 + t126 * t203;
t355 = t55 * t289;
t28 = t354 - t355;
t90 = t154 * t345 - t156 * t240 + t158 * t241;
t352 = t90 * t285;
t89 = t153 * t345 - t155 * t240 + t157 * t241;
t353 = t89 * t289;
t380 = t28 + t352 - t353;
t191 = Icges(5,5) * t248 - Icges(5,6) * t247 + Icges(5,3) * t344;
t193 = Icges(5,4) * t248 - Icges(5,2) * t247 + Icges(5,6) * t344;
t195 = Icges(5,1) * t248 - Icges(5,4) * t247 + Icges(5,5) * t344;
t106 = t191 * t344 - t193 * t247 + t195 * t248;
t192 = Icges(5,5) * t250 - Icges(5,6) * t249 + Icges(5,3) * t342;
t194 = Icges(5,4) * t250 - Icges(5,2) * t249 + Icges(5,6) * t342;
t196 = Icges(5,1) * t250 - Icges(5,4) * t249 + Icges(5,5) * t342;
t107 = t192 * t344 - t194 * t247 + t196 * t248;
t11 = -t65 * t275 + (t285 * t49 + t289 * t50) * t274;
t216 = -Icges(5,3) * t275 + (Icges(5,5) * t287 - Icges(5,6) * t283) * t274;
t217 = -Icges(5,6) * t275 + (Icges(5,4) * t287 - Icges(5,2) * t283) * t274;
t218 = -Icges(5,5) * t275 + (Icges(5,1) * t287 - Icges(5,4) * t283) * t274;
t139 = t216 * t344 - t217 * t247 + t218 * t248;
t35 = -t96 * t275 + (t285 * t81 + t289 * t82) * t274;
t379 = t11 + t35 - t139 * t275 + (t106 * t285 + t107 * t289) * t274;
t108 = t191 * t342 - t249 * t193 + t250 * t195;
t109 = t192 * t342 - t249 * t194 + t250 * t196;
t12 = -t66 * t275 + (t285 * t51 + t289 * t52) * t274;
t140 = t216 * t342 - t249 * t217 + t250 * t218;
t36 = -t97 * t275 + (t285 * t83 + t289 * t84) * t274;
t378 = t12 + t36 - t140 * t275 + (t108 * t285 + t109 * t289) * t274;
t377 = -t106 * t289 + t107 * t285 + t382;
t376 = -t108 * t289 + t109 * t285 + t381;
t103 = t185 * t345 - t240 * t186 + t241 * t187;
t76 = t240 * t146 + t202 * t147 + t203 * t148;
t375 = -t103 - t76;
t159 = t211 * rSges(6,1) - t210 * rSges(6,2) + t247 * rSges(6,3);
t160 = t213 * rSges(6,1) - t212 * rSges(6,2) + t249 * rSges(6,3);
t374 = t285 * t159 + t289 * t160;
t310 = -t248 * rSges(5,1) + t247 * rSges(5,2);
t197 = rSges(5,3) * t344 - t310;
t198 = t250 * rSges(5,1) - t249 * rSges(5,2) + rSges(5,3) * t342;
t373 = t285 * t197 + t289 * t198;
t303 = Icges(4,5) * t275 - Icges(4,6) * t274;
t222 = -Icges(4,3) * t289 + t285 * t303;
t223 = Icges(4,3) * t285 + t289 * t303;
t279 = t289 ^ 2;
t348 = Icges(4,4) * t275;
t305 = -Icges(4,2) * t274 + t348;
t225 = Icges(4,6) * t285 + t289 * t305;
t349 = Icges(4,4) * t274;
t307 = Icges(4,1) * t275 - t349;
t227 = Icges(4,5) * t285 + t289 * t307;
t301 = -t225 * t274 + t227 * t275;
t224 = -Icges(4,6) * t289 + t285 * t305;
t226 = -Icges(4,5) * t289 + t285 * t307;
t302 = t224 * t274 - t226 * t275;
t372 = -t279 * t222 - (t301 * t285 + (-t223 + t302) * t289) * t285 - t377;
t278 = t285 ^ 2;
t371 = t210 / 0.2e1;
t370 = t212 / 0.2e1;
t369 = t240 / 0.2e1;
t368 = t247 / 0.2e1;
t367 = t249 / 0.2e1;
t366 = -t275 / 0.2e1;
t365 = t285 / 0.2e1;
t364 = -t289 / 0.2e1;
t284 = sin(qJ(2));
t362 = pkin(2) * t284;
t361 = pkin(3) * t275;
t360 = pkin(6) * t240;
t288 = cos(qJ(2));
t359 = rSges(3,1) * t288;
t358 = rSges(3,2) * t284;
t357 = t289 * rSges(3,3);
t356 = t289 * rSges(4,3);
t351 = Icges(3,4) * t284;
t350 = Icges(3,4) * t288;
t114 = -t275 * t191 + (-t193 * t283 + t195 * t287) * t274;
t347 = t114 * t289;
t115 = -t275 * t192 + (-t194 * t283 + t196 * t287) * t274;
t346 = t115 * t285;
t341 = t275 * t289;
t149 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t240;
t256 = pkin(3) * t274 - pkin(5) * t275;
t336 = -t149 - t256;
t188 = rSges(6,1) * t241 - rSges(6,2) * t240 + rSges(6,3) * t345;
t335 = -t188 - t256;
t297 = rSges(4,1) * t341 - rSges(4,2) * t342 + t285 * rSges(4,3);
t311 = rSges(4,1) * t275 - rSges(4,2) * t274;
t182 = t285 * (t285 * t311 - t356) + t289 * t297;
t219 = -t275 * rSges(5,3) + (rSges(5,1) * t287 - rSges(5,2) * t283) * t274;
t334 = -t219 - t256;
t313 = pkin(5) * t274 + t361;
t331 = pkin(3) * t341 + pkin(5) * t342;
t333 = t278 * t313 + t289 * t331;
t273 = pkin(2) * t288 + pkin(1);
t266 = t289 * t273;
t332 = t278 * (-pkin(1) + t273) + t289 * (-pkin(1) * t289 + t266);
t330 = t285 * rSges(3,3) + t289 * t359;
t329 = t278 + t279;
t328 = t56 / 0.2e1 + t66 / 0.2e1;
t327 = t65 / 0.2e1 + t55 / 0.2e1;
t129 = t181 * rSges(7,1) + t180 * rSges(7,2) + t212 * rSges(7,3);
t326 = t266 + t331;
t324 = t345 / 0.2e1;
t255 = rSges(4,1) * t274 + rSges(4,2) * t275;
t321 = -t255 - t362;
t320 = -t256 - t362;
t319 = (t278 * t223 + (t302 * t289 + (-t222 + t301) * t285) * t289 + t376) * t285;
t318 = t332 + t333;
t3 = t210 * t49 + t212 * t50 + t240 * t65;
t4 = t210 * t51 + t212 * t52 + t240 * t66;
t317 = t23 * t371 + t24 * t370 + t28 * t369 + t3 * t364 + t4 * t365;
t316 = -t149 + t320;
t315 = -t188 + t320;
t314 = -t219 + t320;
t312 = -t358 + t359;
t309 = -t179 * rSges(7,1) - t178 * rSges(7,2);
t308 = Icges(3,1) * t288 - t351;
t306 = -Icges(3,2) * t284 + t350;
t304 = Icges(3,5) * t288 - Icges(3,6) * t284;
t253 = Icges(4,2) * t275 + t349;
t254 = Icges(4,1) * t274 + t348;
t298 = -t253 * t274 + t254 * t275;
t296 = t97 / 0.2e1 + t90 / 0.2e1 + t328;
t295 = t96 / 0.2e1 + t89 / 0.2e1 + t327;
t294 = (-t273 - t313) * t285;
t293 = t380 * t324 + t384 * t364 + t383 * t365 + t381 * t367 + t382 * t368;
t292 = t289 * t372 + t319;
t128 = t210 * rSges(7,3) - t309;
t67 = t285 * t128 + t289 * t129 + (t210 * t285 + t212 * t289) * pkin(6) + t333;
t291 = (t346 - t347 + t380) * t366 + t378 * t365 + t379 * t364 + t377 * t344 / 0.2e1 + t376 * t342 / 0.2e1;
t252 = Icges(4,5) * t274 + Icges(4,6) * t275;
t290 = -t347 / 0.2e1 + t346 / 0.2e1 - t355 / 0.2e1 + t354 / 0.2e1 - t353 / 0.2e1 + t352 / 0.2e1 + (t225 * t275 + t227 * t274 + t285 * t252 + t289 * t298 + t140 + t66 + t97) * t365 + (t224 * t275 + t226 * t274 - t289 * t252 + t285 * t298 + t139 + t65 + t96) * t364;
t265 = rSges(2,1) * t289 - t285 * rSges(2,2);
t264 = -t285 * rSges(2,1) - rSges(2,2) * t289;
t263 = rSges(3,1) * t284 + rSges(3,2) * t288;
t235 = Icges(3,3) * t285 + t289 * t304;
t234 = -Icges(3,3) * t289 + t285 * t304;
t231 = t289 * t360;
t230 = t285 * t360;
t229 = (pkin(1) - t358) * t289 + t330;
t228 = t357 + (-pkin(1) - t312) * t285;
t221 = t321 * t289;
t220 = t321 * t285;
t209 = t266 + t297;
t208 = t356 + (-t273 - t311) * t285;
t207 = t218 * t343;
t201 = t334 * t289;
t200 = t334 * t285;
t199 = t289 * (-t289 * t358 + t330) + (t285 * t312 - t357) * t285;
t184 = t314 * t289;
t183 = t314 * t285;
t168 = t335 * t289;
t167 = t335 * t285;
t166 = t326 + t198;
t165 = (-t361 - t273 + (-rSges(5,3) - pkin(5)) * t274) * t285 + t310;
t164 = t315 * t289;
t163 = t315 * t285;
t162 = -t275 * t198 - t219 * t342;
t161 = t197 * t275 + t219 * t344;
t152 = t332 + t182;
t145 = -t275 * t216 - t217 * t345 + t207;
t144 = (t197 * t289 - t198 * t285) * t274;
t142 = t326 + t160;
t141 = t294 - t159;
t134 = t289 * t336 - t231;
t133 = t285 * t336 - t230;
t132 = t289 * t316 - t231;
t131 = t285 * t316 - t230;
t130 = t333 + t373;
t121 = -t275 * t160 - t188 * t342;
t120 = t159 * t275 + t188 * t344;
t117 = t160 * t345 - t188 * t249;
t116 = -t159 * t345 + t188 * t247;
t111 = t318 + t373;
t105 = (t159 * t289 - t160 * t285) * t274;
t102 = t159 * t249 - t160 * t247;
t101 = pkin(6) * t212 + t129 + t326;
t100 = (-rSges(7,3) - pkin(6)) * t210 + t294 + t309;
t99 = t103 * t345;
t98 = t333 + t374;
t93 = t318 + t374;
t92 = t129 * t240 - t149 * t212;
t91 = -t128 * t240 + t149 * t210;
t86 = -t149 * t342 - t275 * t129 + (-t212 * t275 - t240 * t342) * pkin(6);
t85 = t149 * t344 + t275 * t128 + (t210 * t275 + t240 * t344) * pkin(6);
t79 = t129 * t345 - t249 * t149 + (t212 * t345 - t240 * t249) * pkin(6);
t78 = -t128 * t345 + t247 * t149 + (-t210 * t345 + t240 * t247) * pkin(6);
t77 = t128 * t212 - t129 * t210;
t75 = t76 * t345;
t71 = t76 * t240;
t68 = (t128 * t289 - t129 * t285 + (t210 * t289 - t212 * t285) * pkin(6)) * t274;
t62 = t249 * t128 - t247 * t129 + (t210 * t249 - t212 * t247) * pkin(6);
t61 = t67 + t332;
t38 = -t103 * t275 + (t89 * t285 + t90 * t289) * t274;
t37 = t89 * t247 + t90 * t249 + t99;
t15 = -t76 * t275 + (t55 * t285 + t56 * t289) * t274;
t14 = t55 * t247 + t56 * t249 + t75;
t13 = t55 * t210 + t56 * t212 + t71;
t1 = [t288 * (Icges(3,2) * t288 + t351) + t284 * (Icges(3,1) * t284 + t350) + Icges(2,3) + t207 + (-t216 + t253) * t275 + (-t217 * t283 + t254) * t274 + m(7) * (t100 ^ 2 + t101 ^ 2) + m(6) * (t141 ^ 2 + t142 ^ 2) + m(5) * (t165 ^ 2 + t166 ^ 2) + m(4) * (t208 ^ 2 + t209 ^ 2) + m(3) * (t228 ^ 2 + t229 ^ 2) + m(2) * (t264 ^ 2 + t265 ^ 2) - t375; t290 + m(7) * (t100 * t132 + t101 * t131) + m(5) * (t165 * t184 + t166 * t183) + m(6) * (t141 * t164 + t142 * t163) + m(4) * (t208 * t221 + t209 * t220) + m(3) * (-t228 * t289 - t229 * t285) * t263 + (t279 / 0.2e1 + t278 / 0.2e1) * (Icges(3,5) * t284 + Icges(3,6) * t288) + (t288 * (Icges(3,6) * t285 + t289 * t306) + t284 * (Icges(3,5) * t285 + t289 * t308)) * t365 + (t288 * (-Icges(3,6) * t289 + t285 * t306) + t284 * (-Icges(3,5) * t289 + t285 * t308)) * t364; m(7) * (t131 ^ 2 + t132 ^ 2 + t61 ^ 2) + m(6) * (t163 ^ 2 + t164 ^ 2 + t93 ^ 2) + m(5) * (t111 ^ 2 + t183 ^ 2 + t184 ^ 2) + m(4) * (t152 ^ 2 + t220 ^ 2 + t221 ^ 2) + t285 * t278 * t235 + m(3) * (t263 ^ 2 * t329 + t199 ^ 2) + t319 + (-t279 * t234 + (-t285 * t234 + t289 * t235) * t285 + t372) * t289; m(4) * (-t208 * t289 - t209 * t285) * t255 + t290 + m(7) * (t100 * t134 + t101 * t133) + m(6) * (t141 * t168 + t142 * t167) + m(5) * (t165 * t201 + t166 * t200); m(7) * (t131 * t133 + t132 * t134 + t61 * t67) + m(6) * (t163 * t167 + t164 * t168 + t93 * t98) + m(5) * (t111 * t130 + t183 * t200 + t184 * t201) + m(4) * (t182 * t152 + (-t220 * t285 - t221 * t289) * t255) + t292; m(7) * (t133 ^ 2 + t134 ^ 2 + t67 ^ 2) + m(6) * (t167 ^ 2 + t168 ^ 2 + t98 ^ 2) + m(5) * (t130 ^ 2 + t200 ^ 2 + t201 ^ 2) + m(4) * (t255 ^ 2 * t329 + t182 ^ 2) + t292; (-t145 + t375) * t275 + m(7) * (t100 * t85 + t101 * t86) + m(6) * (t120 * t141 + t121 * t142) + m(5) * (t161 * t165 + t162 * t166) + ((t140 / 0.2e1 + t115 / 0.2e1 + t296) * t289 + (t114 / 0.2e1 + t139 / 0.2e1 + t295) * t285) * t274; m(7) * (t131 * t86 + t132 * t85 + t61 * t68) + m(6) * (t105 * t93 + t120 * t164 + t121 * t163) + m(5) * (t111 * t144 + t161 * t184 + t162 * t183) + t291; m(7) * (t133 * t86 + t134 * t85 + t67 * t68) + m(6) * (t105 * t98 + t120 * t168 + t121 * t167) + m(5) * (t130 * t144 + t161 * t201 + t162 * t200) + t291; (t145 * t275 - t15 - t38) * t275 + m(7) * (t68 ^ 2 + t85 ^ 2 + t86 ^ 2) + m(6) * (t105 ^ 2 + t120 ^ 2 + t121 ^ 2) + m(5) * (t144 ^ 2 + t161 ^ 2 + t162 ^ 2) + ((-t115 * t275 + t378) * t289 + (-t114 * t275 + t379) * t285) * t274; t75 + t99 + m(7) * (t100 * t78 + t101 * t79) + m(6) * (t116 * t141 + t117 * t142) + t296 * t249 + t295 * t247; m(7) * (t131 * t79 + t132 * t78 + t61 * t62) + m(6) * (t102 * t93 + t116 * t164 + t117 * t163) + t293; m(7) * (t133 * t79 + t134 * t78 + t62 * t67) + m(6) * (t102 * t98 + t116 * t168 + t117 * t167) + t293; (-t14 / 0.2e1 - t37 / 0.2e1) * t275 + (t12 / 0.2e1 + t36 / 0.2e1) * t249 + (t11 / 0.2e1 + t35 / 0.2e1) * t247 + m(7) * (t62 * t68 + t78 * t85 + t79 * t86) + m(6) * (t102 * t105 + t116 * t120 + t117 * t121) + ((t8 / 0.2e1 + t32 / 0.2e1) * t289 + (t7 / 0.2e1 + t31 / 0.2e1) * t285 + (t15 / 0.2e1 + t38 / 0.2e1) * t283) * t274; (t14 + t37) * t345 + t383 * t249 + t384 * t247 + m(7) * (t62 ^ 2 + t78 ^ 2 + t79 ^ 2) + m(6) * (t102 ^ 2 + t116 ^ 2 + t117 ^ 2); m(7) * (t100 * t91 + t101 * t92) + t71 + t328 * t212 + t327 * t210; m(7) * (t131 * t92 + t132 * t91 + t61 * t77) + t317; m(7) * (t133 * t92 + t134 * t91 + t67 * t77) + t317; t11 * t371 + t12 * t370 + m(7) * (t68 * t77 + t85 * t91 + t86 * t92) + t15 * t369 + t13 * t366 + (t289 * t4 / 0.2e1 + t3 * t365) * t274; t13 * t324 + m(7) * (t62 * t77 + t78 * t91 + t79 * t92) + t4 * t367 + t3 * t368 + t7 * t371 + t8 * t370 + t14 * t369; t212 * t4 + t210 * t3 + t240 * t13 + m(7) * (t77 ^ 2 + t91 ^ 2 + t92 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
