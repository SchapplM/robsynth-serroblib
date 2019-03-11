% Calculate joint inertia matrix for
% S6RRPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR9_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR9_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR9_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR9_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR9_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR9_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:08:23
% EndTime: 2019-03-09 14:08:40
% DurationCPUTime: 7.12s
% Computational Cost: add. (33082->649), mult. (42775->885), div. (0->0), fcn. (53302->14), ass. (0->302)
t313 = sin(pkin(12));
t315 = cos(pkin(12));
t316 = cos(pkin(6));
t314 = sin(pkin(6));
t319 = sin(qJ(2));
t384 = t314 * t319;
t286 = -t313 * t384 + t315 * t316;
t385 = t313 * t316;
t287 = t315 * t384 + t385;
t322 = cos(qJ(2));
t382 = t314 * t322;
t222 = Icges(4,4) * t287 + Icges(4,2) * t286 - Icges(4,6) * t382;
t223 = Icges(4,1) * t287 + Icges(4,4) * t286 - Icges(4,5) * t382;
t269 = Icges(3,3) * t316 + (Icges(3,5) * t319 + Icges(3,6) * t322) * t314;
t270 = Icges(3,6) * t316 + (Icges(3,4) * t319 + Icges(3,2) * t322) * t314;
t271 = Icges(3,5) * t316 + (Icges(3,1) * t319 + Icges(3,4) * t322) * t314;
t411 = t286 * t222 + t287 * t223 + t316 * t269 + t270 * t382 + t271 * t384;
t312 = pkin(12) + qJ(4);
t342 = qJ(5) + t312;
t305 = sin(t342);
t338 = cos(t342);
t266 = t305 * t384 - t316 * t338;
t330 = t314 * t338;
t267 = t316 * t305 + t319 * t330;
t210 = Icges(6,5) * t267 - Icges(6,6) * t266 - Icges(6,3) * t382;
t211 = Icges(6,4) * t267 - Icges(6,2) * t266 - Icges(6,6) * t382;
t212 = Icges(6,1) * t267 - Icges(6,4) * t266 - Icges(6,5) * t382;
t320 = sin(qJ(1));
t378 = t320 * t322;
t323 = cos(qJ(1));
t379 = t319 * t323;
t289 = t316 * t379 + t378;
t244 = t289 * t305 + t323 * t330;
t381 = t314 * t323;
t245 = t289 * t338 - t305 * t381;
t377 = t322 * t323;
t380 = t319 * t320;
t288 = -t316 * t377 + t380;
t100 = t210 * t288 - t211 * t244 + t212 * t245;
t290 = t316 * t378 + t379;
t318 = sin(qJ(6));
t321 = cos(qJ(6));
t206 = -t245 * t318 + t288 * t321;
t207 = t245 * t321 + t288 * t318;
t124 = Icges(7,5) * t207 + Icges(7,6) * t206 + Icges(7,3) * t244;
t126 = Icges(7,4) * t207 + Icges(7,2) * t206 + Icges(7,6) * t244;
t128 = Icges(7,1) * t207 + Icges(7,4) * t206 + Icges(7,5) * t244;
t51 = t124 * t244 + t126 * t206 + t128 * t207;
t291 = -t316 * t380 + t377;
t383 = t314 * t320;
t247 = t291 * t338 + t305 * t383;
t208 = -t247 * t318 + t290 * t321;
t209 = t247 * t321 + t290 * t318;
t246 = t291 * t305 - t320 * t330;
t125 = Icges(7,5) * t209 + Icges(7,6) * t208 + Icges(7,3) * t246;
t127 = Icges(7,4) * t209 + Icges(7,2) * t208 + Icges(7,6) * t246;
t129 = Icges(7,1) * t209 + Icges(7,4) * t208 + Icges(7,5) * t246;
t52 = t125 * t244 + t127 * t206 + t129 * t207;
t242 = -t267 * t318 - t321 * t382;
t243 = t267 * t321 - t318 * t382;
t152 = Icges(7,5) * t243 + Icges(7,6) * t242 + Icges(7,3) * t266;
t153 = Icges(7,4) * t243 + Icges(7,2) * t242 + Icges(7,6) * t266;
t154 = Icges(7,1) * t243 + Icges(7,4) * t242 + Icges(7,5) * t266;
t65 = t152 * t244 + t153 * t206 + t154 * t207;
t11 = t288 * t51 + t290 * t52 - t382 * t65;
t160 = Icges(6,5) * t245 - Icges(6,6) * t244 + Icges(6,3) * t288;
t162 = Icges(6,4) * t245 - Icges(6,2) * t244 + Icges(6,6) * t288;
t164 = Icges(6,1) * t245 - Icges(6,4) * t244 + Icges(6,5) * t288;
t76 = t160 * t288 - t162 * t244 + t164 * t245;
t161 = Icges(6,5) * t247 - Icges(6,6) * t246 + Icges(6,3) * t290;
t163 = Icges(6,4) * t247 - Icges(6,2) * t246 + Icges(6,6) * t290;
t165 = Icges(6,1) * t247 - Icges(6,4) * t246 + Icges(6,5) * t290;
t77 = t161 * t288 - t163 * t244 + t165 * t245;
t410 = -t100 * t382 + t288 * t76 + t290 * t77 + t11;
t101 = t210 * t290 - t211 * t246 + t212 * t247;
t53 = t124 * t246 + t126 * t208 + t128 * t209;
t54 = t125 * t246 + t127 * t208 + t129 * t209;
t66 = t152 * t246 + t153 * t208 + t154 * t209;
t12 = t288 * t53 + t290 * t54 - t382 * t66;
t78 = t160 * t290 - t162 * t246 + t164 * t247;
t79 = t161 * t290 - t163 * t246 + t165 * t247;
t409 = -t101 * t382 + t288 * t78 + t290 * t79 + t12;
t15 = t65 * t316 + (t320 * t52 - t323 * t51) * t314;
t408 = t15 + t100 * t316 + (t320 * t77 - t323 * t76) * t314;
t16 = t66 * t316 + (t320 * t54 - t323 * t53) * t314;
t407 = t16 + t101 * t316 + (t320 * t79 - t323 * t78) * t314;
t368 = -t266 * t211 + t267 * t212;
t106 = -t210 * t382 + t368;
t58 = t125 * t266 + t127 * t242 + t129 * t243;
t389 = t58 * t290;
t57 = t124 * t266 + t126 * t242 + t128 * t243;
t390 = t57 * t288;
t72 = t266 * t152 + t242 * t153 + t243 * t154;
t21 = -t382 * t72 + t389 + t390;
t89 = -t161 * t382 - t163 * t266 + t165 * t267;
t387 = t89 * t290;
t88 = -t160 * t382 - t162 * t266 + t164 * t267;
t388 = t88 * t288;
t406 = -t106 * t382 + t21 + t387 + t388;
t105 = t106 * t316;
t71 = t72 * t316;
t23 = t71 + (t58 * t320 - t57 * t323) * t314;
t405 = t23 + t105 + (t89 * t320 - t88 * t323) * t314;
t333 = -rSges(7,1) * t207 - rSges(7,2) * t206;
t130 = rSges(7,3) * t244 - t333;
t393 = pkin(5) * t245;
t375 = pkin(11) * t244 + t130 + t393;
t155 = rSges(7,1) * t243 + rSges(7,2) * t242 + rSges(7,3) * t266;
t404 = pkin(5) * t267 + pkin(11) * t266 + t155;
t221 = Icges(4,5) * t287 + Icges(4,6) * t286 - Icges(4,3) * t382;
t403 = (-t221 * t382 + t411) * t316;
t402 = t244 / 0.2e1;
t401 = t246 / 0.2e1;
t400 = t266 / 0.2e1;
t399 = t288 / 0.2e1;
t398 = t290 / 0.2e1;
t397 = t316 / 0.2e1;
t396 = t320 / 0.2e1;
t395 = -t323 / 0.2e1;
t394 = pkin(3) * t313;
t306 = t315 * pkin(3) + pkin(2);
t392 = -pkin(2) + t306;
t317 = -pkin(9) - qJ(3);
t386 = -t106 - t72;
t376 = t375 * t290;
t131 = t209 * rSges(7,1) + t208 * rSges(7,2) + t246 * rSges(7,3);
t374 = t247 * pkin(5) + pkin(11) * t246 + t131;
t307 = sin(t312);
t295 = pkin(4) * t307 + t394;
t278 = t295 * t381;
t356 = t313 * t381;
t297 = pkin(3) * t356;
t311 = -pkin(10) + t317;
t360 = t311 - t317;
t308 = cos(t312);
t293 = pkin(4) * t308 + t306;
t362 = t293 - t306;
t146 = -t288 * t360 + t289 * t362 - t278 + t297;
t197 = (t295 - t394) * t316 + (t319 * t362 + t322 * t360) * t314;
t373 = t146 * t382 + t288 * t197;
t357 = t313 * t383;
t347 = -pkin(3) * t357 + t290 * t317 - t291 * t306;
t349 = -t290 * t311 + t291 * t293 + t295 * t383;
t147 = t347 + t349;
t167 = t247 * rSges(6,1) - t246 * rSges(6,2) + t290 * rSges(6,3);
t372 = -t147 - t167;
t249 = t291 * pkin(2) + qJ(3) * t290;
t180 = -t249 - t347;
t237 = t316 * t249;
t370 = t316 * t180 + t237;
t279 = t288 * qJ(3);
t179 = -t288 * t317 + t289 * t392 - t279 - t297;
t248 = pkin(2) * t289 + t279;
t369 = -t179 - t248;
t258 = -t289 * t313 - t315 * t381;
t259 = t289 * t315 - t356;
t191 = rSges(4,1) * t259 + rSges(4,2) * t258 + rSges(4,3) * t288;
t367 = -t191 - t248;
t272 = -t307 * t384 + t308 * t316;
t273 = t307 * t316 + t308 * t384;
t217 = Icges(5,4) * t273 + Icges(5,2) * t272 - Icges(5,6) * t382;
t218 = Icges(5,1) * t273 + Icges(5,4) * t272 - Icges(5,5) * t382;
t366 = t272 * t217 + t273 * t218;
t334 = -rSges(6,1) * t245 + rSges(6,2) * t244;
t166 = rSges(6,3) * t288 - t334;
t213 = rSges(6,1) * t267 - rSges(6,2) * t266 - rSges(6,3) * t382;
t117 = t166 * t382 + t288 * t213;
t292 = (pkin(2) * t319 - qJ(3) * t322) * t314;
t364 = -pkin(3) * t385 - ((qJ(3) + t317) * t322 + t392 * t319) * t314 - t292;
t363 = t248 * t383 + t249 * t381;
t361 = t323 * pkin(1) + pkin(8) * t383;
t359 = t57 / 0.2e1 + t65 / 0.2e1;
t358 = t58 / 0.2e1 + t66 / 0.2e1;
t216 = Icges(5,5) * t273 + Icges(5,6) * t272 - Icges(5,3) * t382;
t252 = -t289 * t307 - t308 * t381;
t253 = t289 * t308 - t307 * t381;
t103 = t216 * t288 + t217 * t252 + t218 * t253;
t171 = Icges(5,5) * t253 + Icges(5,6) * t252 + Icges(5,3) * t288;
t173 = Icges(5,4) * t253 + Icges(5,2) * t252 + Icges(5,6) * t288;
t175 = Icges(5,1) * t253 + Icges(5,4) * t252 + Icges(5,5) * t288;
t92 = -t171 * t382 + t173 * t272 + t175 * t273;
t355 = t92 / 0.2e1 + t103 / 0.2e1;
t254 = -t291 * t307 + t308 * t383;
t255 = t291 * t308 + t307 * t383;
t104 = t216 * t290 + t217 * t254 + t218 * t255;
t172 = Icges(5,5) * t255 + Icges(5,6) * t254 + Icges(5,3) * t290;
t174 = Icges(5,4) * t255 + Icges(5,2) * t254 + Icges(5,6) * t290;
t176 = Icges(5,1) * t255 + Icges(5,4) * t254 + Icges(5,5) * t290;
t93 = -t172 * t382 + t174 * t272 + t176 * t273;
t354 = t93 / 0.2e1 + t104 / 0.2e1;
t353 = -t147 - t374;
t352 = t316 * t147 + t370;
t351 = -t146 + t369;
t350 = -t197 + t364;
t178 = t255 * rSges(5,1) + t254 * rSges(5,2) + t290 * rSges(5,3);
t260 = -t291 * t313 + t315 * t383;
t261 = t291 * t315 + t357;
t192 = t261 * rSges(4,1) + t260 * rSges(4,2) + t290 * rSges(4,3);
t233 = t291 * rSges(3,1) - t290 * rSges(3,2) + rSges(3,3) * t383;
t345 = -t382 / 0.2e1;
t343 = t288 * t410 + t409 * t290;
t341 = -t320 * pkin(1) + pkin(8) * t381;
t340 = t314 * (-rSges(4,1) * t287 - rSges(4,2) * t286 + rSges(4,3) * t382 - t292);
t339 = t179 * t383 + t180 * t381 + t363;
t68 = t288 * t404 + t375 * t382;
t219 = rSges(5,1) * t273 + rSges(5,2) * t272 - rSges(5,3) * t382;
t337 = t314 * (-t219 + t364);
t70 = t72 * t266;
t18 = t57 * t244 + t58 * t246 + t70;
t3 = t244 * t51 + t246 * t52 + t266 * t65;
t4 = t244 * t53 + t246 * t54 + t266 * t66;
t336 = t11 * t402 + t12 * t401 + t18 * t345 + t21 * t400 + t3 * t399 + t4 * t398;
t335 = -rSges(5,1) * t253 - rSges(5,2) * t252;
t332 = t349 + t361;
t331 = t314 * (-t213 + t350);
t329 = t146 * t383 + t147 * t381 + t339;
t328 = t314 * (t350 - t404);
t327 = t390 / 0.2e1 + t389 / 0.2e1 + t388 / 0.2e1 + t387 / 0.2e1 + (t65 + t100) * t399 + (t66 + t101) * t398;
t326 = -t289 * t293 + t278 + t341;
t325 = -t382 * t406 + t343;
t232 = t289 * rSges(3,1) - t288 * rSges(3,2) - rSges(3,3) * t381;
t324 = t408 * t399 + t407 * t398 + t406 * t397 + t409 * t383 / 0.2e1 + t405 * t345 - t410 * t381 / 0.2e1;
t300 = rSges(2,1) * t323 - t320 * rSges(2,2);
t299 = -t320 * rSges(2,1) - rSges(2,2) * t323;
t274 = rSges(3,3) * t316 + (rSges(3,1) * t319 + rSges(3,2) * t322) * t314;
t231 = Icges(3,1) * t291 - Icges(3,4) * t290 + Icges(3,5) * t383;
t230 = Icges(3,1) * t289 - Icges(3,4) * t288 - Icges(3,5) * t381;
t229 = Icges(3,4) * t291 - Icges(3,2) * t290 + Icges(3,6) * t383;
t228 = Icges(3,4) * t289 - Icges(3,2) * t288 - Icges(3,6) * t381;
t227 = Icges(3,5) * t291 - Icges(3,6) * t290 + Icges(3,3) * t383;
t226 = Icges(3,5) * t289 - Icges(3,6) * t288 - Icges(3,3) * t381;
t215 = t233 + t361;
t214 = -t232 + t341;
t200 = -t316 * t232 - t274 * t381;
t199 = t233 * t316 - t274 * t383;
t190 = Icges(4,1) * t261 + Icges(4,4) * t260 + Icges(4,5) * t290;
t189 = Icges(4,1) * t259 + Icges(4,4) * t258 + Icges(4,5) * t288;
t188 = Icges(4,4) * t261 + Icges(4,2) * t260 + Icges(4,6) * t290;
t187 = Icges(4,4) * t259 + Icges(4,2) * t258 + Icges(4,6) * t288;
t186 = Icges(4,5) * t261 + Icges(4,6) * t260 + Icges(4,3) * t290;
t185 = Icges(4,5) * t259 + Icges(4,6) * t258 + Icges(4,3) * t288;
t177 = rSges(5,3) * t288 - t335;
t157 = (t232 * t320 + t233 * t323) * t314;
t151 = t269 * t383 - t270 * t290 + t271 * t291;
t150 = -t269 * t381 - t288 * t270 + t289 * t271;
t149 = t290 * t166;
t140 = t249 + t192 + t361;
t139 = t341 + t367;
t138 = t290 * t146;
t135 = t227 * t316 + (t229 * t322 + t231 * t319) * t314;
t134 = t226 * t316 + (t228 * t322 + t230 * t319) * t314;
t133 = -t347 + t178 + t361;
t132 = -t289 * t306 + t297 + (-rSges(5,3) + t317) * t288 + t335 + t341;
t123 = -t178 * t382 - t219 * t290;
t122 = t177 * t382 + t219 * t288;
t120 = t332 + t167;
t119 = (-rSges(6,3) + t311) * t288 + t326 + t334;
t118 = -t167 * t382 - t213 * t290;
t114 = t316 * t367 + t323 * t340;
t113 = t192 * t316 + t320 * t340 + t237;
t112 = -t216 * t382 + t366;
t111 = t177 * t290 - t178 * t288;
t110 = t112 * t316;
t109 = t221 * t290 + t222 * t260 + t223 * t261;
t108 = t221 * t288 + t222 * t258 + t223 * t259;
t107 = -t167 * t288 + t149;
t102 = (t191 * t320 + t192 * t323) * t314 + t363;
t97 = -t186 * t382 + t188 * t286 + t190 * t287;
t96 = -t185 * t382 + t187 * t286 + t189 * t287;
t95 = t332 + t374;
t94 = -t393 + t288 * t311 + (-rSges(7,3) - pkin(11)) * t244 + t326 + t333;
t91 = t131 * t266 - t155 * t246;
t90 = -t130 * t266 + t155 * t244;
t87 = t172 * t290 + t174 * t254 + t176 * t255;
t86 = t171 * t290 + t173 * t254 + t175 * t255;
t85 = t172 * t288 + t174 * t252 + t176 * t253;
t84 = t171 * t288 + t173 * t252 + t175 * t253;
t81 = (-t177 + t369) * t316 + t323 * t337;
t80 = t178 * t316 + t320 * t337 + t370;
t75 = t372 * t382 + (-t197 - t213) * t290;
t74 = t117 + t373;
t73 = t130 * t246 - t131 * t244;
t69 = -t290 * t404 - t374 * t382;
t67 = (t177 * t320 + t178 * t323) * t314 + t339;
t62 = t288 * t372 + t138 + t149;
t61 = -t288 * t374 + t376;
t60 = (-t166 + t351) * t316 + t323 * t331;
t59 = t167 * t316 + t320 * t331 + t352;
t50 = t353 * t382 + (-t197 - t404) * t290;
t49 = t68 + t373;
t48 = (t166 * t320 + t167 * t323) * t314 + t329;
t47 = t288 * t353 + t138 + t376;
t46 = (t351 - t375) * t316 + t323 * t328;
t45 = t316 * t374 + t320 * t328 + t352;
t44 = t110 + (t93 * t320 - t92 * t323) * t314;
t43 = -t112 * t382 + t92 * t288 + t93 * t290;
t38 = t104 * t316 + (t320 * t87 - t323 * t86) * t314;
t37 = t103 * t316 + (t320 * t85 - t323 * t84) * t314;
t36 = -t104 * t382 + t288 * t86 + t290 * t87;
t35 = -t103 * t382 + t288 * t84 + t290 * t85;
t34 = (t320 * t375 + t323 * t374) * t314 + t329;
t1 = [Icges(2,3) + (-t210 - t216 - t221) * t382 + m(7) * (t94 ^ 2 + t95 ^ 2) + m(6) * (t119 ^ 2 + t120 ^ 2) + m(5) * (t132 ^ 2 + t133 ^ 2) + m(4) * (t139 ^ 2 + t140 ^ 2) + m(3) * (t214 ^ 2 + t215 ^ 2) + m(2) * (t299 ^ 2 + t300 ^ 2) + t72 + t366 + t368 + t411; t71 + t110 + t105 + m(7) * (t45 * t95 + t46 * t94) + m(6) * (t119 * t60 + t120 * t59) + m(5) * (t132 * t81 + t133 * t80) + m(4) * (t113 * t140 + t114 * t139) + m(3) * (t199 * t215 + t200 * t214) + ((-t88 / 0.2e1 - t96 / 0.2e1 - t134 / 0.2e1 - t150 / 0.2e1 - t100 / 0.2e1 - t108 / 0.2e1 - t355 - t359) * t323 + (t97 / 0.2e1 + t135 / 0.2e1 + t151 / 0.2e1 + t101 / 0.2e1 + t109 / 0.2e1 + t89 / 0.2e1 + t354 + t358) * t320) * t314 + t403; m(7) * (t34 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t48 ^ 2 + t59 ^ 2 + t60 ^ 2) + m(5) * (t67 ^ 2 + t80 ^ 2 + t81 ^ 2) + m(4) * (t102 ^ 2 + t113 ^ 2 + t114 ^ 2) + m(3) * (t157 ^ 2 + t199 ^ 2 + t200 ^ 2) + (t38 + ((t186 * t290 + t188 * t260 + t190 * t261) * t320 - (t185 * t290 + t187 * t260 + t189 * t261) * t323) * t314 + (t227 * t383 - t229 * t290 + t231 * t291) * t383 + t407) * t383 + (-t37 - ((t186 * t288 + t188 * t258 + t190 * t259) * t320 - (t185 * t288 + t187 * t258 + t189 * t259) * t323) * t314 + (-t226 * t381 - t288 * t228 + t289 * t230) * t381 + (-t226 * t383 + t227 * t381 + t228 * t290 + t288 * t229 - t230 * t291 - t289 * t231) * t383 - t408) * t381 + (t44 + (t109 + t151) * t383 + (-t150 - t108) * t381 + ((-t134 - t96) * t323 + (t135 + t97) * t320) * t314 + t403 + t405) * t316; m(7) * (t288 * t95 + t290 * t94) + m(6) * (t119 * t290 + t120 * t288) + m(5) * (t132 * t290 + t133 * t288) + m(4) * (t139 * t290 + t140 * t288); m(7) * (t288 * t45 + t290 * t46 - t34 * t382) + m(6) * (t288 * t59 + t290 * t60 - t382 * t48) + m(5) * (t288 * t80 + t290 * t81 - t382 * t67) + m(4) * (-t102 * t382 + t113 * t288 + t114 * t290); 0.2e1 * (m(4) / 0.2e1 + m(5) / 0.2e1 + m(6) / 0.2e1 + m(7) / 0.2e1) * (t314 ^ 2 * t322 ^ 2 + t288 ^ 2 + t290 ^ 2); t354 * t290 + t355 * t288 + (-t112 + t386) * t382 + m(7) * (t49 * t94 + t50 * t95) + m(6) * (t119 * t74 + t120 * t75) + m(5) * (t122 * t132 + t123 * t133) + t327; t324 + t43 * t397 + m(7) * (t34 * t47 + t45 * t50 + t46 * t49) + m(6) * (t48 * t62 + t59 * t75 + t60 * t74) + m(5) * (t111 * t67 + t122 * t81 + t123 * t80) + t37 * t399 + (t36 * t396 + t35 * t395 - t322 * t44 / 0.2e1) * t314 + t38 * t398; m(5) * (-t111 * t382 + t122 * t290 + t123 * t288) + m(6) * (t288 * t75 + t290 * t74 - t382 * t62) + m(7) * (t288 * t50 + t290 * t49 - t382 * t47); t288 * t35 + t290 * t36 + (-t43 - t406) * t382 + m(7) * (t47 ^ 2 + t49 ^ 2 + t50 ^ 2) + m(6) * (t62 ^ 2 + t74 ^ 2 + t75 ^ 2) + m(5) * (t111 ^ 2 + t122 ^ 2 + t123 ^ 2) + t343; t386 * t382 + m(7) * (t68 * t94 + t69 * t95) + m(6) * (t117 * t119 + t118 * t120) + t327; t324 + m(7) * (t34 * t61 + t45 * t69 + t46 * t68) + m(6) * (t107 * t48 + t117 * t60 + t118 * t59); m(6) * (-t107 * t382 + t117 * t290 + t118 * t288) + m(7) * (t288 * t69 + t290 * t68 - t382 * t61); m(7) * (t47 * t61 + t49 * t68 + t50 * t69) + m(6) * (t107 * t62 + t117 * t74 + t118 * t75) + t325; m(7) * (t61 ^ 2 + t68 ^ 2 + t69 ^ 2) + m(6) * (t107 ^ 2 + t117 ^ 2 + t118 ^ 2) + t325; t70 + m(7) * (t90 * t94 + t91 * t95) + t358 * t246 + t359 * t244; t15 * t402 + m(7) * (t34 * t73 + t45 * t91 + t46 * t90) + t16 * t401 + t23 * t400 + t18 * t397 + (t3 * t395 + t396 * t4) * t314; m(7) * (t288 * t91 + t290 * t90 - t382 * t73); m(7) * (t47 * t73 + t49 * t90 + t50 * t91) + t336; m(7) * (t61 * t73 + t68 * t90 + t69 * t91) + t336; t246 * t4 + t244 * t3 + t266 * t18 + m(7) * (t73 ^ 2 + t90 ^ 2 + t91 ^ 2);];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
