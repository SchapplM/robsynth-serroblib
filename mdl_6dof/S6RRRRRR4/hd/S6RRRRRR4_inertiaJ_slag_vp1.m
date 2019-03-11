% Calculate joint inertia matrix for
% S6RRRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 03:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR4_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR4_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR4_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 03:46:46
% EndTime: 2019-03-10 03:46:59
% DurationCPUTime: 6.13s
% Computational Cost: add. (27172->626), mult. (24508->843), div. (0->0), fcn. (26590->12), ass. (0->309)
t307 = sin(qJ(2));
t407 = Icges(3,5) * t307;
t406 = t407 / 0.2e1;
t305 = qJ(3) + qJ(4);
t296 = qJ(5) + t305;
t293 = qJ(6) + t296;
t277 = sin(t293);
t278 = cos(t293);
t308 = sin(qJ(1));
t310 = cos(qJ(2));
t311 = cos(qJ(1));
t378 = t310 * t311;
t222 = -t277 * t378 + t278 * t308;
t223 = t277 * t308 + t278 * t378;
t381 = t307 * t311;
t157 = t223 * rSges(7,1) + t222 * rSges(7,2) + rSges(7,3) * t381;
t309 = cos(qJ(3));
t292 = t309 * pkin(3) + pkin(2);
t295 = cos(t305);
t264 = pkin(4) * t295 + t292;
t291 = cos(t296);
t246 = pkin(5) * t291 + t264;
t294 = sin(t305);
t306 = sin(qJ(3));
t266 = t306 * pkin(3) + pkin(4) * t294;
t290 = sin(t296);
t252 = pkin(5) * t290 + t266;
t405 = t246 * t378 + t308 * t252 + t157;
t404 = t308 ^ 2;
t403 = t311 ^ 2;
t312 = -pkin(9) - pkin(8);
t382 = t307 * t308;
t379 = t308 * t310;
t220 = -t277 * t379 - t278 * t311;
t221 = -t277 * t311 + t278 * t379;
t150 = Icges(7,5) * t221 + Icges(7,6) * t220 + Icges(7,3) * t382;
t152 = Icges(7,4) * t221 + Icges(7,2) * t220 + Icges(7,6) * t382;
t154 = Icges(7,1) * t221 + Icges(7,4) * t220 + Icges(7,5) * t382;
t56 = t150 * t382 + t152 * t220 + t154 * t221;
t151 = Icges(7,5) * t223 + Icges(7,6) * t222 + Icges(7,3) * t381;
t153 = Icges(7,4) * t223 + Icges(7,2) * t222 + Icges(7,6) * t381;
t155 = Icges(7,1) * t223 + Icges(7,4) * t222 + Icges(7,5) * t381;
t57 = t151 * t382 + t153 * t220 + t155 * t221;
t205 = -Icges(7,3) * t310 + (Icges(7,5) * t278 - Icges(7,6) * t277) * t307;
t206 = -Icges(7,6) * t310 + (Icges(7,4) * t278 - Icges(7,2) * t277) * t307;
t207 = -Icges(7,5) * t310 + (Icges(7,1) * t278 - Icges(7,4) * t277) * t307;
t96 = t205 * t382 + t206 * t220 + t207 * t221;
t5 = -t310 * t96 + (t308 * t56 + t311 * t57) * t307;
t58 = t150 * t381 + t152 * t222 + t154 * t223;
t59 = t151 * t381 + t153 * t222 + t155 * t223;
t97 = t205 * t381 + t206 * t222 + t207 * t223;
t6 = -t310 * t97 + (t308 * t58 + t311 * t59) * t307;
t402 = t6 * t381 + t5 * t382;
t401 = t308 / 0.2e1;
t400 = -t310 / 0.2e1;
t399 = -t311 / 0.2e1;
t398 = t311 / 0.2e1;
t397 = pkin(2) * t310;
t396 = pkin(8) * t307;
t395 = -pkin(2) + t292;
t213 = -Icges(6,5) * t310 + (Icges(6,1) * t291 - Icges(6,4) * t290) * t307;
t199 = t307 * t291 * t213;
t211 = -Icges(6,3) * t310 + (Icges(6,5) * t291 - Icges(6,6) * t290) * t307;
t212 = -Icges(6,6) * t310 + (Icges(6,4) * t291 - Icges(6,2) * t290) * t307;
t388 = t212 * t290;
t119 = -t310 * t211 - t307 * t388 + t199;
t195 = t307 * t278 * t207;
t389 = t206 * t277;
t114 = -t310 * t205 - t307 * t389 + t195;
t390 = t114 * t310;
t70 = -t150 * t310 + (-t152 * t277 + t154 * t278) * t307;
t71 = -t151 * t310 + (-t153 * t277 + t155 * t278) * t307;
t15 = -t390 + (t308 * t70 + t311 * t71) * t307;
t231 = -t290 * t379 - t291 * t311;
t232 = -t290 * t311 + t291 * t379;
t161 = Icges(6,5) * t232 + Icges(6,6) * t231 + Icges(6,3) * t382;
t163 = Icges(6,4) * t232 + Icges(6,2) * t231 + Icges(6,6) * t382;
t165 = Icges(6,1) * t232 + Icges(6,4) * t231 + Icges(6,5) * t382;
t78 = -t161 * t310 + (-t163 * t290 + t165 * t291) * t307;
t233 = -t290 * t378 + t291 * t308;
t234 = t290 * t308 + t291 * t378;
t162 = Icges(6,5) * t234 + Icges(6,6) * t233 + Icges(6,3) * t381;
t164 = Icges(6,4) * t234 + Icges(6,2) * t233 + Icges(6,6) * t381;
t166 = Icges(6,1) * t234 + Icges(6,4) * t233 + Icges(6,5) * t381;
t79 = -t162 * t310 + (-t164 * t290 + t166 * t291) * t307;
t394 = -t15 + t119 * t310 - (t308 * t78 + t311 * t79) * t307;
t393 = t311 * rSges(3,3);
t391 = Icges(3,4) * t310;
t225 = -Icges(5,6) * t310 + (Icges(5,4) * t295 - Icges(5,2) * t294) * t307;
t387 = t225 * t294;
t240 = -Icges(4,6) * t310 + (Icges(4,4) * t309 - Icges(4,2) * t306) * t307;
t386 = t240 * t306;
t385 = t252 * t311;
t384 = t306 * t308;
t383 = t306 * t311;
t380 = t307 * t312;
t377 = -t114 - t119;
t304 = -pkin(10) + t312;
t298 = -pkin(11) + t304;
t361 = t311 * t266 + t304 * t382;
t365 = t246 - t264;
t120 = -t385 + (-t298 * t307 + t310 * t365) * t308 + t361;
t328 = -rSges(7,1) * t221 - rSges(7,2) * t220;
t156 = rSges(7,3) * t382 - t328;
t145 = t156 * t381;
t376 = t120 * t381 + t145;
t356 = t298 - t304;
t363 = t264 * t378 + t308 * t266;
t375 = -t356 * t381 - t363 + t405;
t358 = pkin(3) * t383 + t308 * t380;
t360 = t264 - t292;
t141 = t360 * t379 + t358 - t361;
t134 = t141 * t381;
t329 = -rSges(6,1) * t232 - rSges(6,2) * t231;
t167 = rSges(6,3) * t382 - t329;
t149 = t167 * t381;
t374 = t134 + t149;
t354 = t304 - t312;
t196 = t307 * t360 + t310 * t354;
t373 = t310 * t141 + t196 * t382;
t359 = -pkin(3) * t384 - t292 * t378;
t142 = -t354 * t381 + t359 + t363;
t168 = t234 * rSges(6,1) + t233 * rSges(6,2) + rSges(6,3) * t381;
t372 = -t142 - t168;
t169 = t307 * t365 + t310 * t356;
t208 = -rSges(7,3) * t310 + (rSges(7,1) * t278 - rSges(7,2) * t277) * t307;
t371 = -t169 - t208;
t250 = -t294 * t378 + t295 * t308;
t251 = t294 * t308 + t295 * t378;
t178 = t251 * rSges(5,1) + t250 * rSges(5,2) + rSges(5,3) * t381;
t320 = -t311 * t380 - t359;
t357 = pkin(2) * t378 + pkin(8) * t381;
t194 = t320 - t357;
t370 = -t178 - t194;
t193 = (t310 * t395 - t396) * t308 - t358;
t219 = (pkin(8) + t312) * t310 + t395 * t307;
t369 = t310 * t193 + t219 * t382;
t217 = -rSges(6,3) * t310 + (rSges(6,1) * t291 - rSges(6,2) * t290) * t307;
t368 = -t196 - t217;
t124 = t310 * t156 + t208 * t382;
t128 = t310 * t167 + t217 * t382;
t248 = -t294 * t379 - t295 * t311;
t249 = -t294 * t311 + t295 * t379;
t330 = -rSges(5,1) * t249 - rSges(5,2) * t248;
t177 = rSges(5,3) * t382 - t330;
t229 = -rSges(5,3) * t310 + (rSges(5,1) * t295 - rSges(5,2) * t294) * t307;
t132 = t310 * t177 + t229 * t382;
t366 = -t219 - t229;
t247 = -rSges(4,3) * t310 + (rSges(4,1) * t309 - rSges(4,2) * t306) * t307;
t275 = t307 * pkin(2) - t310 * pkin(8);
t364 = -t247 - t275;
t362 = t404 * (t396 + t397) + t311 * t357;
t355 = t311 * pkin(1) + t308 * pkin(7);
t226 = -Icges(5,5) * t310 + (Icges(5,1) * t295 - Icges(5,4) * t294) * t307;
t201 = t307 * t295 * t226;
t224 = -Icges(5,3) * t310 + (Icges(5,5) * t295 - Icges(5,6) * t294) * t307;
t127 = -t310 * t224 - t307 * t387 + t201;
t171 = Icges(5,5) * t249 + Icges(5,6) * t248 + Icges(5,3) * t382;
t173 = Icges(5,4) * t249 + Icges(5,2) * t248 + Icges(5,6) * t382;
t175 = Icges(5,1) * t249 + Icges(5,4) * t248 + Icges(5,5) * t382;
t82 = -t171 * t310 + (-t173 * t294 + t175 * t295) * t307;
t172 = Icges(5,5) * t251 + Icges(5,6) * t250 + Icges(5,3) * t381;
t174 = Icges(5,4) * t251 + Icges(5,2) * t250 + Icges(5,6) * t381;
t176 = Icges(5,1) * t251 + Icges(5,4) * t250 + Icges(5,5) * t381;
t83 = -t172 * t310 + (-t174 * t294 + t176 * t295) * t307;
t353 = t127 * t310 - (t308 * t82 + t311 * t83) * t307 + t394;
t237 = -Icges(4,3) * t310 + (Icges(4,5) * t309 - Icges(4,6) * t306) * t307;
t243 = -Icges(4,5) * t310 + (Icges(4,1) * t309 - Icges(4,4) * t306) * t307;
t258 = -t306 * t379 - t309 * t311;
t259 = t309 * t379 - t383;
t122 = t237 * t382 + t240 * t258 + t243 * t259;
t185 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t382;
t187 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t382;
t189 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t382;
t94 = -t185 * t310 + (-t187 * t306 + t189 * t309) * t307;
t352 = t122 / 0.2e1 + t94 / 0.2e1;
t260 = -t306 * t378 + t308 * t309;
t261 = t309 * t378 + t384;
t123 = t237 * t381 + t240 * t260 + t243 * t261;
t186 = Icges(4,5) * t261 + Icges(4,6) * t260 + Icges(4,3) * t381;
t188 = Icges(4,4) * t261 + Icges(4,2) * t260 + Icges(4,6) * t381;
t190 = Icges(4,1) * t261 + Icges(4,4) * t260 + Icges(4,5) * t381;
t95 = -t186 * t310 + (-t188 * t306 + t190 * t309) * t307;
t351 = t95 / 0.2e1 + t123 / 0.2e1;
t350 = -t127 + t377;
t349 = t134 + t376;
t348 = -t142 - t375;
t347 = -t194 + t372;
t346 = -t196 + t371;
t345 = -t219 + t368;
t344 = -t275 + t366;
t192 = t261 * rSges(4,1) + t260 * rSges(4,2) + rSges(4,3) * t381;
t101 = t211 * t382 + t212 * t231 + t213 * t232;
t61 = t161 * t382 + t163 * t231 + t165 * t232;
t62 = t162 * t382 + t164 * t231 + t166 * t232;
t13 = -t101 * t310 + (t308 * t61 + t311 * t62) * t307;
t102 = t211 * t381 + t212 * t233 + t213 * t234;
t63 = t161 * t381 + t163 * t233 + t165 * t234;
t64 = t162 * t381 + t164 * t233 + t166 * t234;
t14 = -t102 * t310 + (t308 * t63 + t311 * t64) * t307;
t343 = t13 * t382 + t14 * t381 + t402;
t342 = t382 / 0.2e1;
t341 = t381 / 0.2e1;
t340 = (t70 + t96) * t342 + (t71 + t97) * t341;
t339 = -t310 * t15 + t402;
t338 = -t194 + t348;
t54 = t310 * t120 + t169 * t382 + t124;
t337 = -t219 + t346;
t336 = t308 * t193 + t311 * t194 + t362;
t68 = t128 + t373;
t335 = -t275 + t345;
t26 = t308 * t57 - t311 * t56;
t27 = t308 * t59 - t311 * t58;
t334 = t26 * t342 + t27 * t341 + t5 * t399 + t6 * t401 + (t71 * t308 - t70 * t311) * t400;
t109 = t224 * t382 + t225 * t248 + t226 * t249;
t72 = t171 * t382 + t173 * t248 + t175 * t249;
t73 = t172 * t382 + t174 * t248 + t176 * t249;
t20 = -t109 * t310 + (t308 * t72 + t311 * t73) * t307;
t110 = t224 * t381 + t225 * t250 + t226 * t251;
t74 = t171 * t381 + t173 * t250 + t175 * t251;
t75 = t172 * t381 + t174 * t250 + t176 * t251;
t21 = -t110 * t310 + (t308 * t74 + t311 * t75) * t307;
t333 = t20 * t382 + t21 * t381 + t343;
t332 = rSges(3,1) * t310 - rSges(3,2) * t307;
t331 = -rSges(4,1) * t259 - rSges(4,2) * t258;
t327 = -t275 + t337;
t325 = -Icges(3,2) * t307 + t391;
t324 = Icges(3,5) * t310 - Icges(3,6) * t307;
t321 = rSges(3,1) * t378 - rSges(3,2) * t381 + t308 * rSges(3,3);
t319 = t308 * t141 + t311 * t142 + t336;
t44 = t54 + t373;
t318 = t310 * t394 + t343;
t317 = t340 + (t101 + t78) * t342 + (t102 + t79) * t341;
t35 = t308 * t62 - t311 * t61;
t36 = t308 * t64 - t311 * t63;
t316 = t13 * t399 + t14 * t401 + t36 * t341 + t35 * t342 + t334 + (t79 * t308 - t78 * t311) * t400;
t315 = t310 * t353 + t333;
t314 = t317 + (t109 + t82) * t342 + (t110 + t83) * t341;
t40 = t308 * t73 - t311 * t72;
t41 = t308 * t75 - t311 * t74;
t313 = t20 * t399 + t21 * t401 + t41 * t341 + t40 * t342 + t316 + (t83 * t308 - t82 * t311) * t400;
t302 = t311 * pkin(7);
t273 = rSges(2,1) * t311 - rSges(2,2) * t308;
t272 = -rSges(2,1) * t308 - rSges(2,2) * t311;
t271 = rSges(3,1) * t307 + rSges(3,2) * t310;
t268 = Icges(3,6) * t310 + t407;
t239 = Icges(3,3) * t308 + t311 * t324;
t238 = -Icges(3,3) * t311 + t308 * t324;
t216 = t307 * t309 * t243;
t210 = t321 + t355;
t209 = t393 + t302 + (-pkin(1) - t332) * t308;
t198 = t364 * t311;
t197 = t364 * t308;
t191 = rSges(4,3) * t382 - t331;
t180 = t311 * t321 + (t308 * t332 - t393) * t308;
t179 = t193 * t381;
t160 = t177 * t381;
t147 = t192 + t355 + t357;
t146 = t302 + (-t397 - pkin(1) + (-rSges(4,3) - pkin(8)) * t307) * t308 + t331;
t144 = t344 * t311;
t143 = t344 * t308;
t140 = -t192 * t310 - t247 * t381;
t139 = t191 * t310 + t247 * t382;
t135 = -t310 * t237 - t307 * t386 + t216;
t133 = -t178 * t310 - t229 * t381;
t131 = t320 + t178 + t355;
t130 = t302 + (-rSges(5,3) * t307 - t292 * t310 - pkin(1)) * t308 + t330 + t358;
t129 = -t168 * t310 - t217 * t381;
t126 = (t191 * t311 - t192 * t308) * t307;
t125 = -t157 * t310 - t208 * t381;
t117 = -t304 * t381 + t168 + t355 + t363;
t116 = t302 + (-rSges(6,3) * t307 - t264 * t310 - pkin(1)) * t308 + t329 + t361;
t113 = t335 * t311;
t112 = t335 * t308;
t111 = -t178 * t382 + t160;
t108 = -t168 * t382 + t149;
t105 = -t298 * t381 + t355 + t405;
t104 = t385 + t302 + (-t246 * t310 - pkin(1) + (-rSges(7,3) + t298) * t307) * t308 + t328;
t103 = t191 * t308 + t192 * t311 + t362;
t100 = -t157 * t382 + t145;
t91 = t310 * t370 + t366 * t381;
t90 = t132 + t369;
t89 = t186 * t381 + t188 * t260 + t190 * t261;
t88 = t185 * t381 + t187 * t260 + t189 * t261;
t87 = t186 * t382 + t188 * t258 + t190 * t259;
t86 = t185 * t382 + t187 * t258 + t189 * t259;
t85 = t327 * t311;
t84 = t327 * t308;
t69 = t310 * t372 + t368 * t381;
t67 = t370 * t382 + t160 + t179;
t60 = t177 * t308 + t178 * t311 + t336;
t55 = -t310 * t375 + t371 * t381;
t53 = t372 * t382 + t374;
t52 = t310 * t347 + t345 * t381;
t51 = t68 + t369;
t50 = t308 * t89 - t311 * t88;
t49 = t308 * t87 - t311 * t86;
t48 = -t375 * t382 + t376;
t46 = t347 * t382 + t179 + t374;
t45 = t310 * t348 + t346 * t381;
t42 = t167 * t308 + t168 * t311 + t319;
t34 = -t123 * t310 + (t308 * t88 + t311 * t89) * t307;
t33 = -t122 * t310 + (t308 * t86 + t311 * t87) * t307;
t29 = t310 * t338 + t337 * t381;
t28 = t44 + t369;
t22 = t348 * t382 + t349;
t8 = t338 * t382 + t179 + t349;
t7 = t375 * t311 + (t120 + t156) * t308 + t319;
t1 = [Icges(2,3) + t195 + t199 + t201 + t216 + (Icges(3,4) * t307 + Icges(3,2) * t310 - t205 - t211 - t224 - t237) * t310 + (Icges(3,1) * t307 - t386 - t387 - t388 - t389 + t391) * t307 + m(7) * (t104 ^ 2 + t105 ^ 2) + m(6) * (t116 ^ 2 + t117 ^ 2) + m(5) * (t130 ^ 2 + t131 ^ 2) + m(4) * (t146 ^ 2 + t147 ^ 2) + m(3) * (t209 ^ 2 + t210 ^ 2) + m(2) * (t272 ^ 2 + t273 ^ 2); (-t70 / 0.2e1 + (-Icges(3,6) * t311 + t308 * t325) * t400 + t311 * t406 - t109 / 0.2e1 - t96 / 0.2e1 - t101 / 0.2e1 - t82 / 0.2e1 - t78 / 0.2e1 + t268 * t398 - t352) * t311 + ((Icges(3,6) * t308 + t311 * t325) * t310 / 0.2e1 + t308 * t406 + t83 / 0.2e1 + t79 / 0.2e1 + t71 / 0.2e1 + t102 / 0.2e1 + t110 / 0.2e1 + t97 / 0.2e1 + t268 * t401 + t351) * t308 + m(4) * (t146 * t198 + t147 * t197) + m(7) * (t104 * t85 + t105 * t84) + m(6) * (t112 * t117 + t113 * t116) + m(5) * (t130 * t144 + t131 * t143) + m(3) * (-t209 * t311 - t210 * t308) * t271; m(7) * (t7 ^ 2 + t84 ^ 2 + t85 ^ 2) + m(6) * (t112 ^ 2 + t113 ^ 2 + t42 ^ 2) + m(5) * (t143 ^ 2 + t144 ^ 2 + t60 ^ 2) + m(4) * (t103 ^ 2 + t197 ^ 2 + t198 ^ 2) + m(3) * (t180 ^ 2 + (t403 + t404) * t271 ^ 2) + (-t403 * t238 - t26 - t35 - t40 - t49) * t311 + (t404 * t239 + t27 + t36 + t41 + t50 + (-t308 * t238 + t311 * t239) * t311) * t308; t314 + (t308 * t352 + t311 * t351) * t307 + (-t135 + t350) * t310 + m(7) * (t104 * t28 + t105 * t29) + m(6) * (t116 * t51 + t117 * t52) + m(5) * (t130 * t90 + t131 * t91) + m(4) * (t139 * t146 + t140 * t147); m(7) * (t28 * t85 + t29 * t84 + t7 * t8) + m(6) * (t112 * t52 + t113 * t51 + t42 * t46) + m(5) * (t143 * t91 + t144 * t90 + t60 * t67) + m(4) * (t103 * t126 + t139 * t198 + t140 * t197) + (t398 * t50 + t401 * t49) * t307 + t33 * t399 + (t95 * t308 - t94 * t311) * t400 + t34 * t401 + t313; (t308 * t33 + t311 * t34) * t307 + (t135 * t310 + (-t308 * t94 - t311 * t95) * t307 + t353) * t310 + m(5) * (t67 ^ 2 + t90 ^ 2 + t91 ^ 2) + m(4) * (t126 ^ 2 + t139 ^ 2 + t140 ^ 2) + m(7) * (t28 ^ 2 + t29 ^ 2 + t8 ^ 2) + m(6) * (t46 ^ 2 + t51 ^ 2 + t52 ^ 2) + t333; t350 * t310 + m(7) * (t104 * t44 + t105 * t45) + m(6) * (t116 * t68 + t117 * t69) + m(5) * (t130 * t132 + t131 * t133) + t314; m(7) * (t22 * t7 + t44 * t85 + t45 * t84) + m(6) * (t112 * t69 + t113 * t68 + t42 * t53) + m(5) * (t111 * t60 + t132 * t144 + t133 * t143) + t313; m(7) * (t22 * t8 + t28 * t44 + t29 * t45) + m(6) * (t46 * t53 + t51 * t68 + t52 * t69) + m(5) * (t111 * t67 + t132 * t90 + t133 * t91) + t315; m(7) * (t22 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t53 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(5) * (t111 ^ 2 + t132 ^ 2 + t133 ^ 2) + t315; t377 * t310 + m(7) * (t104 * t54 + t105 * t55) + m(6) * (t116 * t128 + t117 * t129) + t317; m(7) * (t48 * t7 + t54 * t85 + t55 * t84) + m(6) * (t108 * t42 + t112 * t129 + t113 * t128) + t316; m(7) * (t28 * t54 + t29 * t55 + t48 * t8) + m(6) * (t108 * t46 + t128 * t51 + t129 * t52) + t318; m(7) * (t22 * t48 + t44 * t54 + t45 * t55) + m(6) * (t108 * t53 + t128 * t68 + t129 * t69) + t318; m(7) * (t48 ^ 2 + t54 ^ 2 + t55 ^ 2) + m(6) * (t108 ^ 2 + t128 ^ 2 + t129 ^ 2) + t318; -t390 + m(7) * (t104 * t124 + t105 * t125) + t340; m(7) * (t100 * t7 + t124 * t85 + t125 * t84) + t334; m(7) * (t100 * t8 + t124 * t28 + t125 * t29) + t339; m(7) * (t100 * t22 + t124 * t44 + t125 * t45) + t339; m(7) * (t100 * t48 + t124 * t54 + t125 * t55) + t339; m(7) * (t100 ^ 2 + t124 ^ 2 + t125 ^ 2) + t339;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
