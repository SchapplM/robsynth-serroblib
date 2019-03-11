% Calculate joint inertia matrix for
% S6RRRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR6_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR6_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR6_inertiaJ_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR6_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR6_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR6_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:18:43
% EndTime: 2019-03-09 22:19:01
% DurationCPUTime: 7.89s
% Computational Cost: add. (20886->590), mult. (19836->811), div. (0->0), fcn. (21355->12), ass. (0->293)
t303 = qJ(3) + qJ(4);
t289 = pkin(11) + t303;
t283 = sin(t289);
t284 = cos(t289);
t309 = cos(qJ(1));
t306 = sin(qJ(1));
t308 = cos(qJ(2));
t372 = t306 * t308;
t224 = -t283 * t372 - t284 * t309;
t225 = -t283 * t309 + t284 * t372;
t305 = sin(qJ(2));
t375 = t305 * t306;
t154 = Icges(6,5) * t225 + Icges(6,6) * t224 + Icges(6,3) * t375;
t290 = sin(t303);
t291 = cos(t303);
t242 = -t290 * t372 - t291 * t309;
t243 = -t290 * t309 + t291 * t372;
t165 = Icges(5,5) * t243 + Icges(5,6) * t242 + Icges(5,3) * t375;
t416 = t154 + t165;
t371 = t308 * t309;
t226 = -t283 * t371 + t284 * t306;
t227 = t283 * t306 + t284 * t371;
t374 = t305 * t309;
t155 = Icges(6,5) * t227 + Icges(6,6) * t226 + Icges(6,3) * t374;
t244 = -t290 * t371 + t291 * t306;
t245 = t290 * t306 + t291 * t371;
t166 = Icges(5,5) * t245 + Icges(5,6) * t244 + Icges(5,3) * t374;
t415 = t155 + t166;
t156 = Icges(6,4) * t225 + Icges(6,2) * t224 + Icges(6,6) * t375;
t158 = Icges(6,1) * t225 + Icges(6,4) * t224 + Icges(6,5) * t375;
t167 = Icges(5,4) * t243 + Icges(5,2) * t242 + Icges(5,6) * t375;
t169 = Icges(5,1) * t243 + Icges(5,4) * t242 + Icges(5,5) * t375;
t414 = t156 * t224 + t158 * t225 + t167 * t242 + t169 * t243 + t375 * t416;
t157 = Icges(6,4) * t227 + Icges(6,2) * t226 + Icges(6,6) * t374;
t159 = Icges(6,1) * t227 + Icges(6,4) * t226 + Icges(6,5) * t374;
t168 = Icges(5,4) * t245 + Icges(5,2) * t244 + Icges(5,6) * t374;
t170 = Icges(5,1) * t245 + Icges(5,4) * t244 + Icges(5,5) * t374;
t413 = t157 * t224 + t159 * t225 + t168 * t242 + t170 * t243 + t375 * t415;
t412 = t156 * t226 + t158 * t227 + t167 * t244 + t169 * t245 + t374 * t416;
t411 = t157 * t226 + t159 * t227 + t168 * t244 + t170 * t245 + t374 * t415;
t220 = -Icges(5,3) * t308 + (Icges(5,5) * t291 - Icges(5,6) * t290) * t305;
t221 = -Icges(5,6) * t308 + (Icges(5,4) * t291 - Icges(5,2) * t290) * t305;
t222 = -Icges(5,5) * t308 + (Icges(5,1) * t291 - Icges(5,4) * t290) * t305;
t105 = t220 * t375 + t221 * t242 + t222 * t243;
t205 = -Icges(6,3) * t308 + (Icges(6,5) * t284 - Icges(6,6) * t283) * t305;
t206 = -Icges(6,6) * t308 + (Icges(6,4) * t284 - Icges(6,2) * t283) * t305;
t207 = -Icges(6,5) * t308 + (Icges(6,1) * t284 - Icges(6,4) * t283) * t305;
t98 = t205 * t375 + t206 * t224 + t207 * t225;
t410 = -t98 - t105;
t106 = t220 * t374 + t221 * t244 + t222 * t245;
t99 = t205 * t374 + t206 * t226 + t207 * t227;
t409 = -t99 - t106;
t408 = -t205 - t220;
t407 = -t206 * t283 - t221 * t290;
t406 = Icges(3,5) * t305;
t405 = (t207 * t284 + t222 * t291) * t305;
t404 = t406 / 0.2e1;
t403 = t410 * t308 + (t414 * t306 + t413 * t309) * t305;
t402 = t409 * t308 + (t412 * t306 + t411 * t309) * t305;
t401 = t413 * t306 - t414 * t309;
t400 = t411 * t306 - t412 * t309;
t75 = -t308 * t154 + (-t156 * t283 + t158 * t284) * t305;
t81 = -t308 * t165 + (-t167 * t290 + t169 * t291) * t305;
t399 = -t75 - t81;
t76 = -t308 * t155 + (-t157 * t283 + t159 * t284) * t305;
t82 = -t308 * t166 + (-t168 * t290 + t170 * t291) * t305;
t398 = t76 + t82;
t397 = t407 * t305 + t408 * t308 + t405;
t287 = qJ(6) + t289;
t272 = sin(t287);
t273 = cos(t287);
t215 = -t272 * t371 + t273 * t306;
t216 = t272 * t306 + t273 * t371;
t151 = t216 * rSges(7,1) + t215 * rSges(7,2) + rSges(7,3) * t374;
t307 = cos(qJ(3));
t288 = t307 * pkin(3) + pkin(2);
t259 = pkin(4) * t291 + t288;
t231 = pkin(5) * t284 + t259;
t304 = sin(qJ(3));
t261 = t304 * pkin(3) + pkin(4) * t290;
t246 = pkin(5) * t283 + t261;
t396 = t231 * t371 + t306 * t246 + t151;
t301 = t306 ^ 2;
t302 = t309 ^ 2;
t395 = m(6) / 0.2e1;
t394 = m(7) / 0.2e1;
t310 = -pkin(9) - pkin(8);
t213 = -t272 * t372 - t273 * t309;
t214 = -t272 * t309 + t273 * t372;
t144 = Icges(7,5) * t214 + Icges(7,6) * t213 + Icges(7,3) * t375;
t146 = Icges(7,4) * t214 + Icges(7,2) * t213 + Icges(7,6) * t375;
t148 = Icges(7,1) * t214 + Icges(7,4) * t213 + Icges(7,5) * t375;
t53 = t144 * t375 + t146 * t213 + t148 * t214;
t145 = Icges(7,5) * t216 + Icges(7,6) * t215 + Icges(7,3) * t374;
t147 = Icges(7,4) * t216 + Icges(7,2) * t215 + Icges(7,6) * t374;
t149 = Icges(7,1) * t216 + Icges(7,4) * t215 + Icges(7,5) * t374;
t54 = t145 * t375 + t147 * t213 + t149 * t214;
t198 = -Icges(7,3) * t308 + (Icges(7,5) * t273 - Icges(7,6) * t272) * t305;
t199 = -Icges(7,6) * t308 + (Icges(7,4) * t273 - Icges(7,2) * t272) * t305;
t200 = -Icges(7,5) * t308 + (Icges(7,1) * t273 - Icges(7,4) * t272) * t305;
t91 = t198 * t375 + t199 * t213 + t200 * t214;
t5 = -t91 * t308 + (t306 * t53 + t309 * t54) * t305;
t55 = t144 * t374 + t146 * t215 + t148 * t216;
t56 = t145 * t374 + t147 * t215 + t149 * t216;
t92 = t198 * t374 + t199 * t215 + t200 * t216;
t6 = -t92 * t308 + (t306 * t55 + t309 * t56) * t305;
t393 = t6 * t374 + t5 * t375;
t392 = t306 / 0.2e1;
t391 = -t308 / 0.2e1;
t390 = -t309 / 0.2e1;
t389 = t309 / 0.2e1;
t388 = pkin(2) * t308;
t387 = pkin(8) * t305;
t386 = -pkin(2) + t288;
t385 = t309 * rSges(3,3);
t383 = Icges(3,4) * t308;
t189 = t305 * t273 * t200;
t381 = t199 * t272;
t107 = -t308 * t198 - t305 * t381 + t189;
t382 = t107 * t308;
t235 = -Icges(4,6) * t308 + (Icges(4,4) * t307 - Icges(4,2) * t304) * t305;
t378 = t235 * t304;
t377 = t304 * t306;
t376 = t304 * t309;
t373 = t305 * t310;
t370 = t309 * t246;
t299 = -qJ(5) + t310;
t292 = -pkin(10) + t299;
t351 = t292 - t299;
t358 = t259 * t371 + t306 * t261;
t369 = -t351 * t374 - t358 + t396;
t353 = pkin(3) * t376 + t306 * t373;
t355 = t259 - t288;
t356 = t309 * t261 + t299 * t375;
t135 = t355 * t372 + t353 - t356;
t128 = t135 * t374;
t324 = -t225 * rSges(6,1) - t224 * rSges(6,2);
t160 = rSges(6,3) * t375 - t324;
t368 = t160 * t374 + t128;
t349 = t299 - t310;
t190 = t305 * t355 + t308 * t349;
t367 = t308 * t135 + t190 * t375;
t354 = -pkin(3) * t377 - t288 * t371;
t136 = -t349 * t374 + t354 + t358;
t161 = t227 * rSges(6,1) + t226 * rSges(6,2) + rSges(6,3) * t374;
t366 = -t136 - t161;
t172 = t245 * rSges(5,1) + t244 * rSges(5,2) + rSges(5,3) * t374;
t315 = -t309 * t373 - t354;
t352 = pkin(2) * t371 + pkin(8) * t374;
t188 = t315 - t352;
t365 = -t172 - t188;
t187 = (t308 * t386 - t387) * t306 - t353;
t217 = (pkin(8) + t310) * t308 + t386 * t305;
t364 = t308 * t187 + t217 * t375;
t210 = -t308 * rSges(6,3) + (rSges(6,1) * t284 - rSges(6,2) * t283) * t305;
t363 = -t190 - t210;
t323 = -t214 * rSges(7,1) - t213 * rSges(7,2);
t150 = rSges(7,3) * t375 - t323;
t202 = -t308 * rSges(7,3) + (rSges(7,1) * t273 - rSges(7,2) * t272) * t305;
t118 = t308 * t150 + t202 * t375;
t325 = -t243 * rSges(5,1) - t242 * rSges(5,2);
t171 = rSges(5,3) * t375 - t325;
t223 = -t308 * rSges(5,3) + (rSges(5,1) * t291 - rSges(5,2) * t290) * t305;
t126 = t308 * t171 + t223 * t375;
t361 = -t217 - t223;
t360 = t231 - t259;
t241 = -t308 * rSges(4,3) + (rSges(4,1) * t307 - rSges(4,2) * t304) * t305;
t270 = t305 * pkin(2) - t308 * pkin(8);
t359 = -t241 - t270;
t357 = t301 * (t387 + t388) + t309 * t352;
t350 = t309 * pkin(1) + t306 * pkin(7);
t348 = t301 + t302;
t64 = -t308 * t144 + (-t146 * t272 + t148 * t273) * t305;
t65 = -t308 * t145 + (-t147 * t272 + t149 * t273) * t305;
t15 = -t382 + (t306 * t64 + t309 * t65) * t305;
t347 = -t15 + t397 * t308 + (t399 * t306 - t398 * t309) * t305;
t232 = -Icges(4,3) * t308 + (Icges(4,5) * t307 - Icges(4,6) * t304) * t305;
t238 = -Icges(4,5) * t308 + (Icges(4,1) * t307 - Icges(4,4) * t304) * t305;
t253 = -t304 * t372 - t307 * t309;
t254 = t307 * t372 - t376;
t120 = t232 * t375 + t235 * t253 + t238 * t254;
t179 = Icges(4,5) * t254 + Icges(4,6) * t253 + Icges(4,3) * t375;
t181 = Icges(4,4) * t254 + Icges(4,2) * t253 + Icges(4,6) * t375;
t183 = Icges(4,1) * t254 + Icges(4,4) * t253 + Icges(4,5) * t375;
t93 = -t308 * t179 + (-t181 * t304 + t183 * t307) * t305;
t346 = t93 / 0.2e1 + t120 / 0.2e1;
t255 = -t304 * t371 + t306 * t307;
t256 = t307 * t371 + t377;
t121 = t232 * t374 + t235 * t255 + t238 * t256;
t180 = Icges(4,5) * t256 + Icges(4,6) * t255 + Icges(4,3) * t374;
t182 = Icges(4,4) * t256 + Icges(4,2) * t255 + Icges(4,6) * t374;
t184 = Icges(4,1) * t256 + Icges(4,4) * t255 + Icges(4,5) * t374;
t94 = -t308 * t180 + (-t182 * t304 + t184 * t307) * t305;
t345 = t94 / 0.2e1 + t121 / 0.2e1;
t344 = -t107 - t397;
t116 = -t370 + (-t292 * t305 + t308 * t360) * t306 + t356;
t139 = t150 * t374;
t343 = t116 * t374 + t128 + t139;
t342 = -t136 - t369;
t341 = -t188 + t366;
t163 = t305 * t360 + t308 * t351;
t340 = -t163 - t190 - t202;
t339 = -t217 + t363;
t338 = -t270 + t361;
t186 = t256 * rSges(4,1) + t255 * rSges(4,2) + rSges(4,3) * t374;
t337 = t375 / 0.2e1;
t336 = t374 / 0.2e1;
t335 = (t64 + t91) * t337 + (t65 + t92) * t336;
t334 = -t308 * t15 + t393;
t333 = -t188 + t342;
t332 = -t217 + t340;
t331 = t306 * t187 + t309 * t188 + t357;
t66 = t308 * t160 + t210 * t375 + t367;
t330 = -t270 + t339;
t26 = t306 * t54 - t309 * t53;
t27 = t306 * t56 - t309 * t55;
t329 = t26 * t337 + t27 * t336 + t5 * t390 + t6 * t392 + (t65 * t306 - t64 * t309) * t391;
t328 = t402 * t374 + t403 * t375 + t393;
t327 = rSges(3,1) * t308 - rSges(3,2) * t305;
t326 = -t254 * rSges(4,1) - t253 * rSges(4,2);
t322 = -t270 + t332;
t320 = -Icges(3,2) * t305 + t383;
t319 = Icges(3,5) * t308 - Icges(3,6) * t305;
t316 = rSges(3,1) * t371 - rSges(3,2) * t374 + t306 * rSges(3,3);
t314 = t306 * t135 + t309 * t136 + t331;
t44 = t308 * t116 + t163 * t375 + t118 + t367;
t313 = t308 * t347 + t328;
t312 = t335 + (-t399 - t410) * t337 + (t398 - t409) * t336;
t311 = t329 + t402 * t392 + (t398 * t306 + t399 * t309) * t391 + t403 * t390 + t401 * t337 + t400 * t336;
t297 = t309 * pkin(7);
t269 = rSges(2,1) * t309 - rSges(2,2) * t306;
t268 = -rSges(2,1) * t306 - rSges(2,2) * t309;
t267 = rSges(3,1) * t305 + rSges(3,2) * t308;
t263 = Icges(3,6) * t308 + t406;
t234 = Icges(3,3) * t306 + t309 * t319;
t233 = -Icges(3,3) * t309 + t306 * t319;
t211 = t305 * t307 * t238;
t204 = t316 + t350;
t203 = t385 + t297 + (-pkin(1) - t327) * t306;
t193 = t359 * t309;
t192 = t359 * t306;
t185 = rSges(4,3) * t375 - t326;
t174 = t309 * t316 + (t306 * t327 - t385) * t306;
t173 = t187 * t374;
t162 = t171 * t374;
t142 = t186 + t350 + t352;
t141 = t297 + (-t388 - pkin(1) + (-rSges(4,3) - pkin(8)) * t305) * t306 + t326;
t138 = t338 * t309;
t137 = t338 * t306;
t134 = -t186 * t308 - t241 * t374;
t133 = t185 * t308 + t241 * t375;
t129 = -t308 * t232 - t305 * t378 + t211;
t127 = -t308 * t172 - t223 * t374;
t125 = t315 + t172 + t350;
t124 = t297 + (-rSges(5,3) * t305 - t288 * t308 - pkin(1)) * t306 + t325 + t353;
t122 = (t185 * t309 - t186 * t306) * t305;
t119 = -t308 * t151 - t202 * t374;
t113 = -t299 * t374 + t161 + t350 + t358;
t112 = t297 + (-rSges(6,3) * t305 - t259 * t308 - pkin(1)) * t306 + t324 + t356;
t110 = t330 * t309;
t109 = t330 * t306;
t108 = -t172 * t375 + t162;
t102 = t185 * t306 + t186 * t309 + t357;
t101 = -t292 * t374 + t350 + t396;
t100 = t370 + t297 + (-t231 * t308 - pkin(1) + (-rSges(7,3) + t292) * t305) * t306 + t323;
t97 = -t151 * t375 + t139;
t88 = t308 * t365 + t361 * t374;
t87 = t126 + t364;
t86 = t180 * t374 + t182 * t255 + t184 * t256;
t85 = t179 * t374 + t181 * t255 + t183 * t256;
t84 = t180 * t375 + t182 * t253 + t184 * t254;
t83 = t179 * t375 + t181 * t253 + t183 * t254;
t80 = t322 * t309;
t79 = t322 * t306;
t68 = t365 * t375 + t162 + t173;
t67 = t308 * t366 + t363 * t374;
t57 = t171 * t306 + t172 * t309 + t331;
t52 = t366 * t375 + t368;
t51 = t308 * t341 + t339 * t374;
t50 = t66 + t364;
t49 = t306 * t86 - t309 * t85;
t48 = t306 * t84 - t309 * t83;
t46 = t341 * t375 + t173 + t368;
t45 = t308 * t342 + t340 * t374;
t43 = t160 * t306 + t161 * t309 + t314;
t36 = -t121 * t308 + (t306 * t85 + t309 * t86) * t305;
t35 = -t120 * t308 + (t306 * t83 + t309 * t84) * t305;
t31 = t308 * t333 + t332 * t374;
t30 = t44 + t364;
t23 = t342 * t375 + t343;
t12 = t333 * t375 + t173 + t343;
t7 = t369 * t309 + (t116 + t150) * t306 + t314;
t1 = [Icges(2,3) + t189 + t211 + (Icges(3,4) * t305 + Icges(3,2) * t308 - t198 - t232 + t408) * t308 + (Icges(3,1) * t305 - t378 - t381 + t383 + t407) * t305 + m(7) * (t100 ^ 2 + t101 ^ 2) + m(6) * (t112 ^ 2 + t113 ^ 2) + m(5) * (t124 ^ 2 + t125 ^ 2) + m(4) * (t141 ^ 2 + t142 ^ 2) + m(3) * (t203 ^ 2 + t204 ^ 2) + m(2) * (t268 ^ 2 + t269 ^ 2) + t405; (-t81 / 0.2e1 - t75 / 0.2e1 - t64 / 0.2e1 - t105 / 0.2e1 - t98 / 0.2e1 - t91 / 0.2e1 + (-Icges(3,6) * t309 + t306 * t320) * t391 + t309 * t404 + t263 * t389 - t346) * t309 + (t82 / 0.2e1 + t76 / 0.2e1 + t65 / 0.2e1 + t106 / 0.2e1 + t99 / 0.2e1 + t92 / 0.2e1 + t308 * (Icges(3,6) * t306 + t309 * t320) / 0.2e1 + t306 * t404 + t263 * t392 + t345) * t306 + m(7) * (t100 * t80 + t101 * t79) + m(6) * (t109 * t113 + t110 * t112) + m(5) * (t124 * t138 + t125 * t137) + m(4) * (t141 * t193 + t142 * t192) + m(3) * (-t203 * t309 - t204 * t306) * t267; m(7) * (t7 ^ 2 + t79 ^ 2 + t80 ^ 2) + m(6) * (t109 ^ 2 + t110 ^ 2 + t43 ^ 2) + m(5) * (t137 ^ 2 + t138 ^ 2 + t57 ^ 2) + m(4) * (t102 ^ 2 + t192 ^ 2 + t193 ^ 2) + m(3) * (t267 ^ 2 * t348 + t174 ^ 2) + (-t302 * t233 - t26 - t401 - t48) * t309 + (t301 * t234 + t27 + t49 + (-t306 * t233 + t309 * t234) * t309 + t400) * t306; t312 + (t306 * t346 + t309 * t345) * t305 + m(4) * (t133 * t141 + t134 * t142) + m(7) * (t100 * t30 + t101 * t31) + m(6) * (t112 * t50 + t113 * t51) + m(5) * (t124 * t87 + t125 * t88) + (-t129 + t344) * t308; t311 + m(7) * (t12 * t7 + t30 * t80 + t31 * t79) + m(6) * (t109 * t51 + t110 * t50 + t43 * t46) + m(5) * (t137 * t88 + t138 * t87 + t57 * t68) + m(4) * (t102 * t122 + t133 * t193 + t134 * t192) + (t389 * t49 + t392 * t48) * t305 + t36 * t392 + t35 * t390 + (t94 * t306 - t93 * t309) * t391; (t306 * t35 + t309 * t36) * t305 + (t129 * t308 + (-t306 * t93 - t309 * t94) * t305 + t347) * t308 + m(7) * (t12 ^ 2 + t30 ^ 2 + t31 ^ 2) + m(6) * (t46 ^ 2 + t50 ^ 2 + t51 ^ 2) + m(5) * (t68 ^ 2 + t87 ^ 2 + t88 ^ 2) + m(4) * (t122 ^ 2 + t133 ^ 2 + t134 ^ 2) + t328; t312 + t344 * t308 + m(7) * (t100 * t44 + t101 * t45) + m(6) * (t112 * t66 + t113 * t67) + m(5) * (t124 * t126 + t125 * t127); t311 + m(7) * (t23 * t7 + t44 * t80 + t45 * t79) + m(6) * (t109 * t67 + t110 * t66 + t43 * t52) + m(5) * (t108 * t57 + t126 * t138 + t127 * t137); m(7) * (t12 * t23 + t30 * t44 + t31 * t45) + m(6) * (t46 * t52 + t50 * t66 + t51 * t67) + m(5) * (t108 * t68 + t126 * t87 + t127 * t88) + t313; m(7) * (t23 ^ 2 + t44 ^ 2 + t45 ^ 2) + m(6) * (t52 ^ 2 + t66 ^ 2 + t67 ^ 2) + m(5) * (t108 ^ 2 + t126 ^ 2 + t127 ^ 2) + t313; 0.2e1 * ((t100 * t309 + t101 * t306) * t394 + (t112 * t309 + t113 * t306) * t395) * t305; m(7) * (-t308 * t7 + (t306 * t79 + t309 * t80) * t305) + m(6) * (-t308 * t43 + (t109 * t306 + t110 * t309) * t305); m(7) * (-t308 * t12 + (t30 * t309 + t306 * t31) * t305) + m(6) * (-t308 * t46 + (t306 * t51 + t309 * t50) * t305); m(7) * (-t308 * t23 + (t306 * t45 + t309 * t44) * t305) + m(6) * (-t308 * t52 + (t306 * t67 + t309 * t66) * t305); 0.2e1 * (t395 + t394) * (t305 ^ 2 * t348 + t308 ^ 2); m(7) * (t100 * t118 + t101 * t119) - t382 + t335; m(7) * (t118 * t80 + t119 * t79 + t7 * t97) + t329; m(7) * (t118 * t30 + t119 * t31 + t12 * t97) + t334; m(7) * (t118 * t44 + t119 * t45 + t23 * t97) + t334; m(7) * (-t97 * t308 + (t118 * t309 + t119 * t306) * t305); m(7) * (t118 ^ 2 + t119 ^ 2 + t97 ^ 2) + t334;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
