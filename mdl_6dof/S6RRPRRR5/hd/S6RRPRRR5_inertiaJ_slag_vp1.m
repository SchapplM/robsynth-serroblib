% Calculate joint inertia matrix for
% S6RRPRRR5
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
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRPRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaJ_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR5_inertiaJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR5_inertiaJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR5_inertiaJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:38:06
% EndTime: 2019-03-09 13:38:25
% DurationCPUTime: 8.12s
% Computational Cost: add. (36753->664), mult. (85917->922), div. (0->0), fcn. (113238->14), ass. (0->310)
t366 = sin(pkin(12));
t367 = cos(pkin(12));
t372 = sin(qJ(2));
t374 = cos(qJ(2));
t286 = -t372 * t366 + t374 * t367;
t303 = sin(pkin(6));
t274 = t286 * t303;
t312 = t366 * t374 + t367 * t372;
t275 = t312 * t303;
t304 = cos(pkin(6));
t238 = Icges(4,5) * t275 + Icges(4,6) * t274 + Icges(4,3) * t304;
t239 = Icges(4,4) * t275 + Icges(4,2) * t274 + Icges(4,6) * t304;
t240 = Icges(4,1) * t275 + Icges(4,4) * t274 + Icges(4,5) * t304;
t270 = Icges(3,3) * t304 + (Icges(3,5) * t372 + Icges(3,6) * t374) * t303;
t271 = Icges(3,6) * t304 + (Icges(3,4) * t372 + Icges(3,2) * t374) * t303;
t272 = Icges(3,5) * t304 + (Icges(3,1) * t372 + Icges(3,4) * t374) * t303;
t335 = t303 * t372;
t387 = t303 * t374 * t271 + t274 * t239 + t275 * t240 + t272 * t335 + (t238 + t270) * t304;
t386 = t303 ^ 2;
t385 = m(4) / 0.2e1;
t384 = m(5) / 0.2e1;
t383 = m(6) / 0.2e1;
t382 = m(7) / 0.2e1;
t276 = t312 * t304;
t307 = sin(qJ(1));
t309 = cos(qJ(1));
t261 = t276 * t309 + t286 * t307;
t306 = sin(qJ(4));
t373 = cos(qJ(4));
t336 = t303 * t373;
t231 = t261 * t306 + t309 * t336;
t381 = t231 / 0.2e1;
t263 = -t276 * t307 + t286 * t309;
t233 = t263 * t306 - t307 * t336;
t380 = t233 / 0.2e1;
t311 = t304 * t286;
t260 = -t307 * t312 + t309 * t311;
t379 = -t260 / 0.2e1;
t262 = -t307 * t311 - t309 * t312;
t378 = -t262 / 0.2e1;
t264 = t275 * t306 - t304 * t373;
t377 = t264 / 0.2e1;
t376 = -t274 / 0.2e1;
t375 = t304 / 0.2e1;
t371 = pkin(1) * t309;
t308 = cos(qJ(5));
t297 = pkin(5) * t308 + pkin(4);
t370 = -pkin(4) + t297;
t361 = t303 * t309;
t232 = t261 * t373 - t306 * t361;
t302 = qJ(5) + qJ(6);
t299 = sin(t302);
t300 = cos(t302);
t180 = -t232 * t299 - t260 * t300;
t181 = t232 * t300 - t260 * t299;
t117 = Icges(7,5) * t181 + Icges(7,6) * t180 + Icges(7,3) * t231;
t119 = Icges(7,4) * t181 + Icges(7,2) * t180 + Icges(7,6) * t231;
t121 = Icges(7,1) * t181 + Icges(7,4) * t180 + Icges(7,5) * t231;
t265 = t275 * t373 + t304 * t306;
t217 = -t265 * t299 - t274 * t300;
t218 = t265 * t300 - t274 * t299;
t60 = t117 * t264 + t119 * t217 + t121 * t218;
t369 = t60 * t231;
t362 = t303 * t307;
t234 = t263 * t373 + t306 * t362;
t182 = -t234 * t299 - t262 * t300;
t183 = t234 * t300 - t262 * t299;
t118 = Icges(7,5) * t183 + Icges(7,6) * t182 + Icges(7,3) * t233;
t120 = Icges(7,4) * t183 + Icges(7,2) * t182 + Icges(7,6) * t233;
t122 = Icges(7,1) * t183 + Icges(7,4) * t182 + Icges(7,5) * t233;
t61 = t118 * t264 + t120 * t217 + t122 * t218;
t368 = t61 * t233;
t305 = sin(qJ(5));
t365 = t260 * t305;
t364 = t262 * t305;
t363 = t274 * t305;
t298 = pkin(2) * t374 + pkin(1);
t360 = t307 * t298;
t277 = t304 * t372 * pkin(2) + (-pkin(8) - qJ(3)) * t303;
t359 = t309 * t277;
t227 = t231 * pkin(10);
t310 = -pkin(11) - pkin(10);
t345 = pkin(5) * t365;
t115 = -t231 * t310 + t232 * t370 - t227 - t345;
t320 = -t181 * rSges(7,1) - t180 * rSges(7,2);
t123 = t231 * rSges(7,3) - t320;
t358 = t115 + t123;
t176 = t234 * pkin(4) + pkin(10) * t233;
t338 = -pkin(5) * t364 - t233 * t310 + t234 * t297;
t116 = -t176 + t338;
t124 = t183 * rSges(7,1) + t182 * rSges(7,2) + t233 * rSges(7,3);
t357 = t116 + t124;
t189 = -t232 * t305 - t260 * t308;
t190 = t232 * t308 - t365;
t133 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t231;
t175 = t232 * pkin(4) + t227;
t356 = -t133 - t175;
t191 = -t234 * t305 - t262 * t308;
t192 = t234 * t308 - t364;
t134 = t192 * rSges(6,1) + t191 * rSges(6,2) + t233 * rSges(6,3);
t355 = t134 + t176;
t354 = t387 * t304;
t150 = rSges(7,1) * t218 + rSges(7,2) * t217 + rSges(7,3) * t264;
t151 = -pkin(5) * t363 + t370 * t265 + (-pkin(10) - t310) * t264;
t353 = -t150 - t151;
t210 = t263 * pkin(3) - pkin(9) * t262;
t291 = t309 * t298;
t254 = -t371 + t291 + (-t303 * pkin(8) - t277) * t307;
t242 = t304 * t254;
t352 = t304 * t210 + t242;
t209 = t261 * pkin(3) - t260 * pkin(9);
t295 = pkin(8) * t361;
t253 = t359 + t295 + (-pkin(1) + t298) * t307;
t351 = -t209 - t253;
t350 = t253 * t362 + t254 * t361;
t194 = Icges(4,5) * t261 + Icges(4,6) * t260 - Icges(4,3) * t361;
t332 = t309 * t374;
t333 = t307 * t372;
t281 = t304 * t332 - t333;
t331 = t309 * t372;
t334 = t307 * t374;
t282 = t304 * t331 + t334;
t244 = Icges(3,5) * t282 + Icges(3,6) * t281 - Icges(3,3) * t361;
t349 = -t244 - t194;
t195 = Icges(4,5) * t263 + Icges(4,6) * t262 + Icges(4,3) * t362;
t283 = -t304 * t334 - t331;
t284 = -t304 * t333 + t332;
t245 = Icges(3,5) * t284 + Icges(3,6) * t283 + Icges(3,3) * t362;
t348 = t245 + t195;
t287 = pkin(2) * t335 + t304 * qJ(3);
t347 = -pkin(3) * t275 + pkin(9) * t274 - t287;
t49 = t117 * t233 + t119 * t182 + t121 * t183;
t50 = t118 * t233 + t120 * t182 + t122 * t183;
t147 = Icges(7,5) * t218 + Icges(7,6) * t217 + Icges(7,3) * t264;
t148 = Icges(7,4) * t218 + Icges(7,2) * t217 + Icges(7,6) * t264;
t149 = Icges(7,1) * t218 + Icges(7,4) * t217 + Icges(7,5) * t264;
t70 = t147 * t233 + t148 * t182 + t149 * t183;
t10 = t231 * t49 + t233 * t50 + t264 * t70;
t80 = t264 * t147 + t217 * t148 + t218 * t149;
t75 = t80 * t264;
t26 = t368 + t75 + t369;
t47 = t117 * t231 + t119 * t180 + t121 * t181;
t48 = t118 * t231 + t120 * t180 + t122 * t181;
t69 = t147 * t231 + t148 * t180 + t149 * t181;
t9 = t231 * t47 + t233 * t48 + t264 * t69;
t346 = t233 * t10 + t231 * t9 + t264 * t26;
t127 = Icges(6,5) * t190 + Icges(6,6) * t189 + Icges(6,3) * t231;
t129 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t231;
t131 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t231;
t220 = -t265 * t305 - t274 * t308;
t221 = t265 * t308 - t363;
t62 = t127 * t264 + t129 * t220 + t131 * t221;
t161 = Icges(6,5) * t221 + Icges(6,6) * t220 + Icges(6,3) * t264;
t162 = Icges(6,4) * t221 + Icges(6,2) * t220 + Icges(6,6) * t264;
t163 = Icges(6,1) * t221 + Icges(6,4) * t220 + Icges(6,5) * t264;
t73 = t161 * t231 + t162 * t189 + t163 * t190;
t344 = t62 / 0.2e1 + t73 / 0.2e1;
t128 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t233;
t130 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t233;
t132 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t233;
t63 = t128 * t264 + t130 * t220 + t132 * t221;
t74 = t161 * t233 + t162 * t191 + t163 * t192;
t343 = t63 / 0.2e1 + t74 / 0.2e1;
t89 = t264 * t161 + t220 * t162 + t221 * t163;
t205 = Icges(5,5) * t265 - Icges(5,6) * t264 - Icges(5,3) * t274;
t206 = Icges(5,4) * t265 - Icges(5,2) * t264 - Icges(5,6) * t274;
t207 = Icges(5,1) * t265 - Icges(5,4) * t264 - Icges(5,5) * t274;
t107 = -t274 * t205 - t264 * t206 + t265 * t207;
t342 = t304 * t176 + t352;
t341 = -t175 + t351;
t216 = t265 * pkin(4) + t264 * pkin(10);
t339 = -t216 + t347;
t159 = t234 * rSges(5,1) - t233 * rSges(5,2) - t262 * rSges(5,3);
t201 = t263 * rSges(4,1) + t262 * rSges(4,2) + rSges(4,3) * t362;
t251 = t284 * rSges(3,1) + t283 * rSges(3,2) + rSges(3,3) * t362;
t330 = t303 * (-rSges(4,1) * t275 - rSges(4,2) * t274 - rSges(4,3) * t304 - t287);
t329 = -t277 * t307 + t291;
t328 = t209 * t362 + t210 * t361 + t350;
t327 = t369 / 0.2e1 + t368 / 0.2e1 + t69 * t381 + t70 * t380 + t75;
t208 = rSges(5,1) * t265 - rSges(5,2) * t264 - rSges(5,3) * t274;
t326 = t303 * (-t208 + t347);
t11 = -t260 * t47 - t262 * t48 - t274 * t69;
t12 = -t260 * t49 - t262 * t50 - t274 * t70;
t77 = t80 * t274;
t28 = -t60 * t260 - t61 * t262 - t77;
t323 = t10 * t378 + t11 * t381 + t12 * t380 + t26 * t376 + t28 * t377 + t9 * t379;
t17 = t304 * t69 + (t307 * t48 - t309 * t47) * t303;
t18 = t304 * t70 + (t307 * t50 - t309 * t49) * t303;
t79 = t80 * t304;
t31 = t79 + (t61 * t307 - t60 * t309) * t303;
t322 = t17 * t381 + t18 * t380 + t26 * t375 + t31 * t377 + t10 * t362 / 0.2e1 - t9 * t361 / 0.2e1;
t321 = -rSges(4,1) * t261 - rSges(4,2) * t260;
t164 = rSges(6,1) * t221 + rSges(6,2) * t220 + rSges(6,3) * t264;
t319 = t303 * (-t164 + t339);
t318 = t175 * t362 + t176 * t361 + t328;
t317 = t303 * (t339 + t353);
t316 = t210 + t329;
t158 = rSges(5,1) * t232 - rSges(5,2) * t231 - rSges(5,3) * t260;
t250 = rSges(3,1) * t282 + rSges(3,2) * t281 - rSges(3,3) * t361;
t315 = -t209 - t359 - t360;
t103 = -t205 * t262 - t206 * t233 + t207 * t234;
t153 = Icges(5,5) * t234 - Icges(5,6) * t233 - Icges(5,3) * t262;
t155 = Icges(5,4) * t234 - Icges(5,2) * t233 - Icges(5,6) * t262;
t157 = Icges(5,1) * t234 - Icges(5,4) * t233 - Icges(5,5) * t262;
t93 = -t153 * t274 - t155 * t264 + t157 * t265;
t314 = -t61 / 0.2e1 - t70 / 0.2e1 - t103 / 0.2e1 - t93 / 0.2e1 - t343;
t102 = -t205 * t260 - t206 * t231 + t207 * t232;
t152 = Icges(5,5) * t232 - Icges(5,6) * t231 - Icges(5,3) * t260;
t154 = Icges(5,4) * t232 - Icges(5,2) * t231 - Icges(5,6) * t260;
t156 = Icges(5,1) * t232 - Icges(5,4) * t231 - Icges(5,5) * t260;
t92 = -t152 * t274 - t154 * t264 + t156 * t265;
t313 = -t69 / 0.2e1 - t60 / 0.2e1 - t92 / 0.2e1 - t102 / 0.2e1 - t344;
t289 = rSges(2,1) * t309 - rSges(2,2) * t307;
t288 = -rSges(2,1) * t307 - rSges(2,2) * t309;
t273 = t304 * rSges(3,3) + (rSges(3,1) * t372 + rSges(3,2) * t374) * t303;
t249 = Icges(3,1) * t284 + Icges(3,4) * t283 + Icges(3,5) * t362;
t248 = Icges(3,1) * t282 + Icges(3,4) * t281 - Icges(3,5) * t361;
t247 = Icges(3,4) * t284 + Icges(3,2) * t283 + Icges(3,6) * t362;
t246 = Icges(3,4) * t282 + Icges(3,2) * t281 - Icges(3,6) * t361;
t230 = pkin(8) * t362 + t251 + t371;
t229 = -pkin(1) * t307 - t250 + t295;
t213 = -t250 * t304 - t273 * t361;
t212 = t251 * t304 - t273 * t362;
t200 = -rSges(4,3) * t361 - t321;
t199 = Icges(4,1) * t263 + Icges(4,4) * t262 + Icges(4,5) * t362;
t198 = Icges(4,1) * t261 + Icges(4,4) * t260 - Icges(4,5) * t361;
t197 = Icges(4,4) * t263 + Icges(4,2) * t262 + Icges(4,6) * t362;
t196 = Icges(4,4) * t261 + Icges(4,2) * t260 - Icges(4,6) * t361;
t193 = (t250 * t307 + t251 * t309) * t303;
t188 = t270 * t362 + t271 * t283 + t272 * t284;
t187 = -t270 * t361 + t271 * t281 + t272 * t282;
t184 = t260 * t216;
t169 = t329 + t201;
t168 = -t360 + (rSges(4,3) * t303 - t277) * t309 + t321;
t167 = t274 * t176;
t166 = t304 * t245 + (t247 * t374 + t249 * t372) * t303;
t165 = t304 * t244 + (t246 * t374 + t248 * t372) * t303;
t160 = t262 * t175;
t141 = t231 * t150;
t140 = (-t200 - t253) * t304 + t309 * t330;
t139 = t201 * t304 + t307 * t330 + t242;
t136 = t238 * t362 + t239 * t262 + t240 * t263;
t135 = -t238 * t361 + t239 * t260 + t240 * t261;
t126 = t316 + t159;
t125 = -t158 + t315;
t114 = (t200 * t307 + t201 * t309) * t303 + t350;
t113 = t264 * t124;
t112 = -t159 * t274 + t208 * t262;
t111 = t158 * t274 - t208 * t260;
t110 = t195 * t304 + t197 * t274 + t199 * t275;
t109 = t194 * t304 + t196 * t274 + t198 * t275;
t108 = t233 * t123;
t106 = t107 * t304;
t105 = t107 * t274;
t104 = -t158 * t262 + t159 * t260;
t101 = (-t158 + t351) * t304 + t309 * t326;
t100 = t159 * t304 + t307 * t326 + t352;
t99 = t316 + t355;
t98 = t315 + t356;
t97 = t134 * t264 - t164 * t233;
t96 = -t133 * t264 + t164 * t231;
t95 = -t150 * t233 + t113;
t94 = -t123 * t264 + t141;
t91 = t316 + t338 + t124;
t90 = t345 - t232 * t297 + (-rSges(7,3) + t310) * t231 + t315 + t320;
t88 = (t158 * t307 + t159 * t309) * t303 + t328;
t87 = t89 * t304;
t86 = -t153 * t262 - t155 * t233 + t157 * t234;
t85 = -t152 * t262 - t154 * t233 + t156 * t234;
t84 = -t153 * t260 - t155 * t231 + t157 * t232;
t83 = -t152 * t260 - t154 * t231 + t156 * t232;
t82 = t89 * t274;
t81 = t89 * t264;
t78 = t133 * t233 - t134 * t231;
t76 = -t124 * t231 + t108;
t72 = -t134 * t274 - t167 + (t164 + t216) * t262;
t71 = -t164 * t260 - t274 * t356 - t184;
t66 = (-t133 + t341) * t304 + t309 * t319;
t65 = t134 * t304 + t307 * t319 + t342;
t64 = -t133 * t262 + t260 * t355 - t160;
t57 = t116 * t264 + t233 * t353 + t113;
t56 = t151 * t231 - t264 * t358 + t141;
t55 = t128 * t233 + t130 * t191 + t132 * t192;
t54 = t127 * t233 + t129 * t191 + t131 * t192;
t53 = t128 * t231 + t130 * t189 + t132 * t190;
t52 = t127 * t231 + t129 * t189 + t131 * t190;
t51 = (t133 * t307 + t134 * t309) * t303 + t318;
t46 = -t167 - t357 * t274 + (t216 - t353) * t262;
t45 = -t184 + t353 * t260 - (-t175 - t358) * t274;
t44 = (t341 - t358) * t304 + t309 * t317;
t43 = t304 * t357 + t307 * t317 + t342;
t42 = t115 * t233 - t231 * t357 + t108;
t41 = t106 + (t93 * t307 - t92 * t309) * t303;
t40 = -t92 * t260 - t93 * t262 - t105;
t39 = -t160 - t358 * t262 + (t176 + t357) * t260;
t38 = (t307 * t358 + t309 * t357) * t303 + t318;
t37 = t103 * t304 + (t307 * t86 - t309 * t85) * t303;
t36 = t102 * t304 + (t307 * t84 - t309 * t83) * t303;
t35 = -t103 * t274 - t260 * t85 - t262 * t86;
t34 = -t102 * t274 - t260 * t83 - t262 * t84;
t33 = t87 + (t63 * t307 - t62 * t309) * t303;
t32 = -t62 * t260 - t63 * t262 - t82;
t30 = t62 * t231 + t63 * t233 + t81;
t22 = t304 * t74 + (t307 * t55 - t309 * t54) * t303;
t21 = t304 * t73 + (t307 * t53 - t309 * t52) * t303;
t20 = -t260 * t54 - t262 * t55 - t274 * t74;
t19 = -t260 * t52 - t262 * t53 - t274 * t73;
t16 = t231 * t54 + t233 * t55 + t264 * t74;
t15 = t231 * t52 + t233 * t53 + t264 * t73;
t1 = [m(7) * (t90 ^ 2 + t91 ^ 2) + m(6) * (t98 ^ 2 + t99 ^ 2) + m(5) * (t125 ^ 2 + t126 ^ 2) + m(4) * (t168 ^ 2 + t169 ^ 2) + m(3) * (t229 ^ 2 + t230 ^ 2) + m(2) * (t288 ^ 2 + t289 ^ 2) + Icges(2,3) + t107 + t89 + t80 + t387; t106 + t79 + t87 + m(7) * (t43 * t91 + t44 * t90) + m(6) * (t65 * t99 + t66 * t98) + m(5) * (t100 * t126 + t101 * t125) + m(4) * (t139 * t169 + t140 * t168) + m(3) * (t212 * t230 + t213 * t229) + ((-t109 / 0.2e1 - t165 / 0.2e1 - t135 / 0.2e1 - t187 / 0.2e1 + t313) * t309 + (t110 / 0.2e1 + t166 / 0.2e1 + t136 / 0.2e1 + t188 / 0.2e1 - t314) * t307) * t303 + t354; (t31 + t33 + t41 + t354) * t304 + m(7) * (t38 ^ 2 + t43 ^ 2 + t44 ^ 2) + m(6) * (t51 ^ 2 + t65 ^ 2 + t66 ^ 2) + m(5) * (t100 ^ 2 + t101 ^ 2 + t88 ^ 2) + m(4) * (t114 ^ 2 + t139 ^ 2 + t140 ^ 2) + m(3) * (t193 ^ 2 + t212 ^ 2 + t213 ^ 2) + ((-t17 - t21 - t36 + ((t196 * t260 + t198 * t261 + t246 * t281 + t248 * t282) * t303 + t349 * t386 * t309) * t309 + (-t109 - t135 - t165 - t187) * t304) * t309 + (t18 + t22 + t37 + ((t197 * t262 + t199 * t263 + t247 * t283 + t249 * t284) * t303 + t348 * t386 * t307) * t307 + (t188 + t136 + t110 + t166) * t304 + (-t196 * t262 - t197 * t260 - t198 * t263 - t199 * t261 - t246 * t283 - t247 * t281 - t248 * t284 - t249 * t282 + (t307 * t349 + t309 * t348) * t303) * t361) * t307) * t303; 0.2e1 * ((t307 * t90 - t309 * t91) * t382 + (t307 * t98 - t309 * t99) * t383 + (t125 * t307 - t126 * t309) * t384 + (t168 * t307 - t169 * t309) * t385) * t303; m(7) * (t304 * t38 + (t307 * t44 - t309 * t43) * t303) + m(6) * (t304 * t51 + (t307 * t66 - t309 * t65) * t303) + m(5) * (t304 * t88 + (-t100 * t309 + t101 * t307) * t303) + m(4) * (t114 * t304 + (-t139 * t309 + t140 * t307) * t303); 0.2e1 * (t385 + t384 + t383 + t382) * (t304 ^ 2 + (t307 ^ 2 + t309 ^ 2) * t386); -t77 - t82 - t105 + m(7) * (t45 * t90 + t46 * t91) + m(6) * (t71 * t98 + t72 * t99) + m(5) * (t111 * t125 + t112 * t126) + t314 * t262 + t313 * t260; (t28 / 0.2e1 + t32 / 0.2e1 + t40 / 0.2e1) * t304 - (t31 / 0.2e1 + t33 / 0.2e1 + t41 / 0.2e1) * t274 + (-t18 / 0.2e1 - t22 / 0.2e1 - t37 / 0.2e1) * t262 + (-t17 / 0.2e1 - t21 / 0.2e1 - t36 / 0.2e1) * t260 + m(7) * (t38 * t39 + t43 * t46 + t44 * t45) + m(6) * (t51 * t64 + t65 * t72 + t66 * t71) + m(5) * (t100 * t112 + t101 * t111 + t104 * t88) + ((-t11 / 0.2e1 - t19 / 0.2e1 - t34 / 0.2e1) * t309 + (t12 / 0.2e1 + t20 / 0.2e1 + t35 / 0.2e1) * t307) * t303; m(5) * (t104 * t304 + (t111 * t307 - t112 * t309) * t303) + m(6) * (t64 * t304 + (t307 * t71 - t309 * t72) * t303) + m(7) * (t39 * t304 + (t307 * t45 - t309 * t46) * t303); -(t28 + t32 + t40) * t274 + (-t12 - t20 - t35) * t262 + (-t11 - t19 - t34) * t260 + m(7) * (t39 ^ 2 + t45 ^ 2 + t46 ^ 2) + m(6) * (t64 ^ 2 + t71 ^ 2 + t72 ^ 2) + m(5) * (t104 ^ 2 + t111 ^ 2 + t112 ^ 2); t81 + t343 * t233 + t344 * t231 + m(7) * (t56 * t90 + t57 * t91) + m(6) * (t96 * t98 + t97 * t99) + t327; t33 * t377 + t30 * t375 + t22 * t380 + t21 * t381 + (t307 * t16 / 0.2e1 - t309 * t15 / 0.2e1) * t303 + m(7) * (t38 * t42 + t43 * t57 + t44 * t56) + m(6) * (t51 * t78 + t65 * t97 + t66 * t96) + t322; m(6) * (t304 * t78 + (t307 * t96 - t309 * t97) * t303) + m(7) * (t304 * t42 + (t307 * t56 - t309 * t57) * t303); t32 * t377 + t15 * t379 + t30 * t376 + t20 * t380 + t19 * t381 + t16 * t378 + m(7) * (t39 * t42 + t45 * t56 + t46 * t57) + m(6) * (t64 * t78 + t71 * t96 + t72 * t97) + t323; t231 * t15 + t233 * t16 + t264 * t30 + m(7) * (t42 ^ 2 + t56 ^ 2 + t57 ^ 2) + m(6) * (t78 ^ 2 + t96 ^ 2 + t97 ^ 2) + t346; m(7) * (t90 * t94 + t91 * t95) + t327; m(7) * (t38 * t76 + t43 * t95 + t44 * t94) + t322; m(7) * (t304 * t76 + (t307 * t94 - t309 * t95) * t303); m(7) * (t39 * t76 + t45 * t94 + t46 * t95) + t323; m(7) * (t42 * t76 + t56 * t94 + t57 * t95) + t346; m(7) * (t76 ^ 2 + t94 ^ 2 + t95 ^ 2) + t346;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
