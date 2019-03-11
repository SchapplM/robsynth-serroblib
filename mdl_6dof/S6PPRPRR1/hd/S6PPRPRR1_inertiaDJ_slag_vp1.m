% Calculate time derivative of joint inertia matrix for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PPRPRR1_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRPRR1_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:07
% EndTime: 2019-03-08 18:42:32
% DurationCPUTime: 17.70s
% Computational Cost: add. (111270->819), mult. (326689->1145), div. (0->0), fcn. (406375->16), ass. (0->363)
t372 = sin(pkin(13));
t373 = sin(pkin(7));
t325 = t373 * t372;
t375 = cos(pkin(13));
t327 = t375 * t373;
t385 = sin(qJ(3));
t387 = cos(qJ(3));
t279 = t325 * t387 + t327 * t385;
t299 = cos(pkin(7));
t313 = t372 * t387 + t375 * t385;
t281 = t313 * t299;
t296 = sin(pkin(12));
t297 = sin(pkin(11));
t300 = cos(pkin(6));
t298 = cos(pkin(12));
t376 = cos(pkin(11));
t342 = t376 * t298;
t285 = -t297 * t296 + t300 * t342;
t343 = t376 * t296;
t286 = t297 * t298 + t300 * t343;
t312 = t372 * t385 - t375 * t387;
t374 = sin(pkin(6));
t328 = t376 * t374;
t235 = -t279 * t328 + t285 * t281 - t286 * t312;
t302 = sin(qJ(5));
t311 = t285 * t373 + t299 * t328;
t386 = cos(qJ(5));
t215 = t235 * t386 - t302 * t311;
t310 = -t325 * t385 + t327 * t387;
t275 = t310 * qJD(3);
t280 = t312 * t299;
t277 = qJD(3) * t280;
t290 = t313 * qJD(3);
t231 = -t275 * t328 - t285 * t277 - t286 * t290;
t158 = qJD(5) * t215 + t231 * t302;
t320 = -t235 * t302 - t311 * t386;
t159 = qJD(5) * t320 + t231 * t386;
t276 = t299 * t290;
t289 = t312 * qJD(3);
t309 = qJD(3) * t279;
t305 = t374 * t309;
t230 = -t276 * t285 + t286 * t289 + t305 * t376;
t120 = rSges(6,1) * t159 - rSges(6,2) * t158 - rSges(6,3) * t230;
t371 = t297 * t300;
t287 = -t298 * t371 - t343;
t288 = -t296 * t371 + t342;
t347 = t297 * t374;
t237 = t279 * t347 + t287 * t281 - t288 * t312;
t315 = t287 * t373 - t299 * t347;
t217 = t237 * t386 - t302 * t315;
t233 = t275 * t347 - t287 * t277 - t288 * t290;
t160 = qJD(5) * t217 + t233 * t302;
t319 = -t237 * t302 - t315 * t386;
t161 = qJD(5) * t319 + t233 * t386;
t232 = -t287 * t276 + t288 * t289 - t297 * t305;
t121 = rSges(6,1) * t161 - rSges(6,2) * t160 - rSges(6,3) * t232;
t307 = t310 * t374;
t234 = -t280 * t285 - t286 * t313 - t307 * t376;
t144 = rSges(6,1) * t215 + rSges(6,2) * t320 - rSges(6,3) * t234;
t236 = -t287 * t280 - t288 * t313 + t297 * t307;
t145 = rSges(6,1) * t217 + rSges(6,2) * t319 - rSges(6,3) * t236;
t57 = -t120 * t236 + t121 * t234 - t144 * t232 + t145 * t230;
t412 = m(6) * t57;
t411 = 2 * m(4);
t410 = 2 * m(5);
t346 = t298 * t374;
t348 = t296 * t374;
t246 = -t276 * t346 + t289 * t348 - t300 * t309;
t254 = -t280 * t346 + t300 * t310 - t313 * t348;
t301 = sin(qJ(6));
t303 = cos(qJ(6));
t163 = -t215 * t301 - t234 * t303;
t164 = t215 * t303 - t234 * t301;
t104 = Icges(7,5) * t164 + Icges(7,6) * t163 - Icges(7,3) * t320;
t106 = Icges(7,4) * t164 + Icges(7,2) * t163 - Icges(7,6) * t320;
t108 = Icges(7,1) * t164 + Icges(7,4) * t163 - Icges(7,5) * t320;
t122 = -qJD(6) * t164 - t159 * t301 - t230 * t303;
t123 = qJD(6) * t163 + t159 * t303 - t230 * t301;
t79 = Icges(7,5) * t123 + Icges(7,6) * t122 + Icges(7,3) * t158;
t81 = Icges(7,4) * t123 + Icges(7,2) * t122 + Icges(7,6) * t158;
t83 = Icges(7,1) * t123 + Icges(7,4) * t122 + Icges(7,5) * t158;
t16 = t104 * t158 + t106 * t122 + t108 * t123 + t163 * t81 + t164 * t83 - t320 * t79;
t165 = -t217 * t301 - t236 * t303;
t166 = t217 * t303 - t236 * t301;
t105 = Icges(7,5) * t166 + Icges(7,6) * t165 - Icges(7,3) * t319;
t107 = Icges(7,4) * t166 + Icges(7,2) * t165 - Icges(7,6) * t319;
t109 = Icges(7,1) * t166 + Icges(7,4) * t165 - Icges(7,5) * t319;
t124 = -qJD(6) * t166 - t161 * t301 - t232 * t303;
t125 = qJD(6) * t165 + t161 * t303 - t232 * t301;
t80 = Icges(7,5) * t125 + Icges(7,6) * t124 + Icges(7,3) * t160;
t82 = Icges(7,4) * t125 + Icges(7,2) * t124 + Icges(7,6) * t160;
t84 = Icges(7,1) * t125 + Icges(7,4) * t124 + Icges(7,5) * t160;
t17 = t105 * t158 + t107 * t122 + t109 * t123 + t163 * t82 + t164 * t84 - t320 * t80;
t255 = t300 * t279 + t281 * t346 - t312 * t348;
t326 = t374 * t373;
t314 = t298 * t326 - t300 * t299;
t240 = t255 * t386 - t302 * t314;
t193 = -t240 * t301 - t254 * t303;
t194 = t240 * t303 - t254 * t301;
t318 = -t255 * t302 - t314 * t386;
t128 = Icges(7,5) * t194 + Icges(7,6) * t193 - Icges(7,3) * t318;
t129 = Icges(7,4) * t194 + Icges(7,2) * t193 - Icges(7,6) * t318;
t130 = Icges(7,1) * t194 + Icges(7,4) * t193 - Icges(7,5) * t318;
t247 = t300 * t275 - t277 * t346 - t290 * t348;
t191 = qJD(5) * t318 + t247 * t386;
t136 = -qJD(6) * t194 - t191 * t301 - t246 * t303;
t137 = qJD(6) * t193 + t191 * t303 - t246 * t301;
t190 = qJD(5) * t240 + t247 * t302;
t95 = Icges(7,5) * t137 + Icges(7,6) * t136 + Icges(7,3) * t190;
t96 = Icges(7,4) * t137 + Icges(7,2) * t136 + Icges(7,6) * t190;
t97 = Icges(7,1) * t137 + Icges(7,4) * t136 + Icges(7,5) * t190;
t30 = t122 * t129 + t123 * t130 + t128 * t158 + t163 * t96 + t164 * t97 - t320 * t95;
t53 = -t104 * t320 + t106 * t163 + t108 * t164;
t54 = -t105 * t320 + t107 * t163 + t109 * t164;
t66 = -t128 * t320 + t129 * t163 + t130 * t164;
t3 = -t16 * t234 - t17 * t236 - t230 * t53 - t232 * t54 - t246 * t66 - t254 * t30;
t114 = Icges(6,5) * t159 - Icges(6,6) * t158 - Icges(6,3) * t230;
t116 = Icges(6,4) * t159 - Icges(6,2) * t158 - Icges(6,6) * t230;
t118 = Icges(6,1) * t159 - Icges(6,4) * t158 - Icges(6,5) * t230;
t138 = Icges(6,5) * t215 + Icges(6,6) * t320 - Icges(6,3) * t234;
t140 = Icges(6,4) * t215 + Icges(6,2) * t320 - Icges(6,6) * t234;
t142 = Icges(6,1) * t215 + Icges(6,4) * t320 - Icges(6,5) * t234;
t43 = -t114 * t234 + t116 * t320 + t118 * t215 - t138 * t230 - t140 * t158 + t142 * t159;
t115 = Icges(6,5) * t161 - Icges(6,6) * t160 - Icges(6,3) * t232;
t117 = Icges(6,4) * t161 - Icges(6,2) * t160 - Icges(6,6) * t232;
t119 = Icges(6,1) * t161 - Icges(6,4) * t160 - Icges(6,5) * t232;
t139 = Icges(6,5) * t217 + Icges(6,6) * t319 - Icges(6,3) * t236;
t141 = Icges(6,4) * t217 + Icges(6,2) * t319 - Icges(6,6) * t236;
t143 = Icges(6,1) * t217 + Icges(6,4) * t319 - Icges(6,5) * t236;
t44 = -t115 * t234 + t117 * t320 + t119 * t215 - t139 * t230 - t141 * t158 + t143 * t159;
t132 = Icges(6,5) * t191 - Icges(6,6) * t190 - Icges(6,3) * t246;
t133 = Icges(6,4) * t191 - Icges(6,2) * t190 - Icges(6,6) * t246;
t134 = Icges(6,1) * t191 - Icges(6,4) * t190 - Icges(6,5) * t246;
t150 = Icges(6,5) * t240 + Icges(6,6) * t318 - Icges(6,3) * t254;
t151 = Icges(6,4) * t240 + Icges(6,2) * t318 - Icges(6,6) * t254;
t152 = Icges(6,1) * t240 + Icges(6,4) * t318 - Icges(6,5) * t254;
t49 = -t132 * t234 + t133 * t320 + t134 * t215 - t150 * t230 - t151 * t158 + t152 * t159;
t75 = -t138 * t234 + t140 * t320 + t142 * t215;
t76 = -t139 * t234 + t141 * t320 + t143 * t215;
t93 = -t150 * t234 + t151 * t320 + t152 * t215;
t409 = -t230 * t75 - t232 * t76 - t234 * t43 - t236 * t44 - t246 * t93 - t254 * t49 + t3;
t18 = t104 * t160 + t106 * t124 + t108 * t125 + t165 * t81 + t166 * t83 - t319 * t79;
t19 = t105 * t160 + t107 * t124 + t109 * t125 + t165 * t82 + t166 * t84 - t319 * t80;
t31 = t124 * t129 + t125 * t130 + t128 * t160 + t165 * t96 + t166 * t97 - t319 * t95;
t55 = -t104 * t319 + t106 * t165 + t108 * t166;
t56 = -t105 * t319 + t107 * t165 + t109 * t166;
t67 = -t128 * t319 + t129 * t165 + t130 * t166;
t4 = -t18 * t234 - t19 * t236 - t230 * t55 - t232 * t56 - t246 * t67 - t254 * t31;
t45 = -t114 * t236 + t116 * t319 + t118 * t217 - t138 * t232 - t140 * t160 + t142 * t161;
t46 = -t115 * t236 + t117 * t319 + t119 * t217 - t139 * t232 - t141 * t160 + t143 * t161;
t50 = -t132 * t236 + t133 * t319 + t134 * t217 - t150 * t232 - t151 * t160 + t152 * t161;
t77 = -t138 * t236 + t140 * t319 + t142 * t217;
t78 = -t139 * t236 + t141 * t319 + t143 * t217;
t94 = -t150 * t236 + t151 * t319 + t152 * t217;
t408 = -t230 * t77 - t232 * t78 - t234 * t45 - t236 * t46 - t246 * t94 - t254 * t50 + t4;
t99 = -t150 * t254 + t151 * t318 + t152 * t240;
t380 = t246 * t99;
t47 = -t114 * t254 + t116 * t318 + t118 * t240 - t138 * t246 - t140 * t190 + t142 * t191;
t48 = -t115 * t254 + t117 * t318 + t119 * t240 - t139 * t246 - t141 * t190 + t143 * t191;
t51 = -t132 * t254 + t133 * t318 + t134 * t240 - t150 * t246 - t151 * t190 + t152 * t191;
t20 = t104 * t190 + t106 * t136 + t108 * t137 + t193 * t81 + t194 * t83 - t318 * t79;
t21 = t105 * t190 + t107 * t136 + t109 * t137 + t193 * t82 + t194 * t84 - t318 * t80;
t34 = t128 * t190 + t129 * t136 + t130 * t137 + t193 * t96 + t194 * t97 - t318 * t95;
t58 = -t104 * t318 + t106 * t193 + t108 * t194;
t59 = -t105 * t318 + t107 * t193 + t109 * t194;
t70 = -t128 * t318 + t129 * t193 + t130 * t194;
t6 = -t20 * t234 - t21 * t236 - t230 * t58 - t232 * t59 - t246 * t70 - t254 * t34;
t87 = -t138 * t254 + t140 * t318 + t142 * t240;
t88 = -t139 * t254 + t141 * t318 + t143 * t240;
t407 = -t230 * t87 - t232 * t88 - t234 * t47 - t236 * t48 - t254 * t51 - t380 + t6;
t382 = pkin(3) * t387;
t406 = pkin(8) + qJ(4);
t405 = t298 * t299;
t330 = t374 * t385;
t322 = t296 * t330;
t331 = t387 * t373;
t332 = t387 * t374;
t269 = t300 * t331 + t332 * t405 - t322;
t265 = t269 * qJD(3);
t329 = t373 * t385;
t304 = -t296 * t332 - t300 * t329 - t330 * t405;
t266 = t304 * qJD(3);
t404 = Icges(4,5) * t265 + Icges(5,5) * t247 + Icges(4,6) * t266 + Icges(5,6) * t246;
t316 = t376 * t326;
t350 = t299 * t387;
t260 = t285 * t350 - t286 * t385 - t316 * t387;
t256 = t260 * qJD(3);
t349 = t299 * t385;
t306 = -t285 * t349 - t286 * t387 + t316 * t385;
t257 = t306 * qJD(3);
t403 = Icges(4,5) * t256 + Icges(5,5) * t231 + Icges(4,6) * t257 + Icges(5,6) * t230;
t321 = t297 * t326;
t262 = t287 * t350 - t288 * t385 + t321 * t387;
t258 = t262 * qJD(3);
t308 = -t287 * t349 - t288 * t387 - t321 * t385;
t259 = t308 * qJD(3);
t402 = Icges(4,5) * t258 + Icges(5,5) * t233 + Icges(4,6) * t259 + Icges(5,6) * t232;
t111 = rSges(7,1) * t166 + rSges(7,2) * t165 - rSges(7,3) * t319;
t369 = pkin(5) * t217 - pkin(10) * t319 + t111;
t110 = rSges(7,1) * t164 + rSges(7,2) * t163 - rSges(7,3) * t320;
t370 = pkin(5) * t215 - pkin(10) * t320 + t110;
t86 = rSges(7,1) * t125 + rSges(7,2) * t124 + rSges(7,3) * t160;
t378 = pkin(5) * t161 + pkin(10) * t160 + t86;
t85 = rSges(7,1) * t123 + rSges(7,2) * t122 + rSges(7,3) * t158;
t379 = pkin(5) * t159 + pkin(10) * t158 + t85;
t32 = t369 * t230 - t370 * t232 + t378 * t234 - t379 * t236;
t401 = 0.2e1 * t311;
t400 = m(7) / 0.2e1;
t399 = t158 / 0.2e1;
t398 = t160 / 0.2e1;
t397 = t190 / 0.2e1;
t396 = -t320 / 0.2e1;
t395 = -t319 / 0.2e1;
t394 = -t230 / 0.2e1;
t393 = -t232 / 0.2e1;
t392 = -t318 / 0.2e1;
t391 = -t246 / 0.2e1;
t390 = -t311 / 0.2e1;
t389 = -t315 / 0.2e1;
t388 = -t314 / 0.2e1;
t182 = rSges(5,1) * t231 + rSges(5,2) * t230;
t202 = rSges(5,1) * t247 + rSges(5,2) * t246;
t381 = pkin(3) * qJD(3);
t291 = t299 * qJD(4) + t331 * t381;
t292 = -qJD(4) * t373 + t350 * t381;
t264 = t300 * t291 + t292 * t346 - t322 * t381;
t245 = t311 * t264;
t335 = t385 * t381;
t248 = t285 * t292 - t286 * t335 - t291 * t328;
t112 = -t202 * t311 - t245 - (-t182 - t248) * t314;
t384 = m(5) * t112;
t183 = rSges(5,1) * t233 + rSges(5,2) * t232;
t249 = t287 * t292 - t288 * t335 + t291 * t347;
t238 = t314 * t249;
t113 = -t183 * t314 - t238 - (-t202 - t264) * t315;
t383 = m(5) * t113;
t98 = rSges(7,1) * t137 + rSges(7,2) * t136 + rSges(7,3) * t190;
t377 = pkin(5) * t191 + pkin(10) * t190 + t98;
t131 = rSges(7,1) * t194 + rSges(7,2) * t193 - rSges(7,3) * t318;
t368 = pkin(5) * t240 - pkin(10) * t318 + t131;
t184 = pkin(4) * t231 - pkin(9) * t230;
t162 = t315 * t184;
t229 = t315 * t248;
t367 = -t162 - t229;
t185 = pkin(4) * t233 - pkin(9) * t232;
t366 = -t185 * t314 - t238;
t186 = pkin(4) * t235 - pkin(9) * t234;
t282 = pkin(3) * t329 + t299 * t406;
t283 = pkin(3) * t349 - t373 * t406;
t212 = pkin(8) * t311 - t282 * t328 + t285 * t283 + t286 * t382;
t195 = t315 * t212;
t365 = -t186 * t315 - t195;
t187 = pkin(4) * t237 - pkin(9) * t236;
t213 = pkin(8) * t315 + t282 * t347 + t287 * t283 + t288 * t382;
t204 = t314 * t213;
t364 = -t187 * t314 - t204;
t363 = -t184 - t248;
t362 = -t185 - t249;
t361 = -t186 - t212;
t360 = -t187 - t213;
t203 = pkin(4) * t247 - pkin(9) * t246;
t359 = -t203 * t311 - t245;
t211 = pkin(4) * t255 - pkin(9) * t254;
t244 = t314 * pkin(8) + t300 * t282 + t283 * t346 + t348 * t382;
t226 = t311 * t244;
t358 = -t211 * t311 - t226;
t357 = -t203 - t264;
t356 = -t211 - t244;
t355 = t249 * t401 - 0.2e1 * t229;
t352 = m(6) * t374;
t351 = m(7) * t374;
t340 = 0.2e1 * m(6);
t338 = 0.2e1 * m(7);
t333 = m(7) * t32 + t412;
t324 = m(6) * t328;
t323 = m(7) * t328;
t103 = -t182 * t315 - t229 - (-t183 - t249) * t311;
t38 = -t379 * t315 - (t362 - t378) * t311 + t367;
t65 = -t120 * t315 - (-t121 + t362) * t311 + t367;
t317 = m(5) * t103 + m(6) * t65 + m(7) * t38;
t37 = t110 * t160 - t111 * t158 - t319 * t85 + t320 * t86;
t253 = rSges(4,1) * t265 + rSges(4,2) * t266;
t252 = Icges(4,1) * t265 + Icges(4,4) * t266;
t251 = Icges(4,4) * t265 + Icges(4,2) * t266;
t243 = -rSges(4,1) * t304 + rSges(4,2) * t269 - rSges(4,3) * t314;
t242 = -Icges(4,1) * t304 + Icges(4,4) * t269 - Icges(4,5) * t314;
t241 = -Icges(4,4) * t304 + Icges(4,2) * t269 - Icges(4,6) * t314;
t225 = rSges(4,1) * t258 + rSges(4,2) * t259;
t224 = rSges(4,1) * t256 + rSges(4,2) * t257;
t223 = Icges(4,1) * t258 + Icges(4,4) * t259;
t222 = Icges(4,1) * t256 + Icges(4,4) * t257;
t221 = Icges(4,4) * t258 + Icges(4,2) * t259;
t220 = Icges(4,4) * t256 + Icges(4,2) * t257;
t210 = -rSges(4,1) * t308 + rSges(4,2) * t262 - rSges(4,3) * t315;
t209 = -rSges(4,1) * t306 + rSges(4,2) * t260 - rSges(4,3) * t311;
t208 = -Icges(4,1) * t308 + Icges(4,4) * t262 - Icges(4,5) * t315;
t207 = -Icges(4,1) * t306 + Icges(4,4) * t260 - Icges(4,5) * t311;
t206 = -Icges(4,4) * t308 + Icges(4,2) * t262 - Icges(4,6) * t315;
t205 = -Icges(4,4) * t306 + Icges(4,2) * t260 - Icges(4,6) * t311;
t201 = Icges(5,1) * t247 + Icges(5,4) * t246;
t200 = Icges(5,4) * t247 + Icges(5,2) * t246;
t198 = rSges(5,1) * t255 + rSges(5,2) * t254 - rSges(5,3) * t314;
t197 = Icges(5,1) * t255 + Icges(5,4) * t254 - Icges(5,5) * t314;
t196 = Icges(5,4) * t255 + Icges(5,2) * t254 - Icges(5,6) * t314;
t181 = Icges(5,1) * t233 + Icges(5,4) * t232;
t180 = Icges(5,1) * t231 + Icges(5,4) * t230;
t179 = Icges(5,4) * t233 + Icges(5,2) * t232;
t178 = Icges(5,4) * t231 + Icges(5,2) * t230;
t173 = rSges(5,1) * t237 + rSges(5,2) * t236 - rSges(5,3) * t315;
t172 = rSges(5,1) * t235 + rSges(5,2) * t234 - rSges(5,3) * t311;
t171 = Icges(5,1) * t237 + Icges(5,4) * t236 - Icges(5,5) * t315;
t170 = Icges(5,1) * t235 + Icges(5,4) * t234 - Icges(5,5) * t311;
t169 = Icges(5,4) * t237 + Icges(5,2) * t236 - Icges(5,6) * t315;
t168 = Icges(5,4) * t235 + Icges(5,2) * t234 - Icges(5,6) * t311;
t155 = -t225 * t314 + t253 * t315;
t154 = t224 * t314 - t253 * t311;
t153 = rSges(6,1) * t240 + rSges(6,2) * t318 - rSges(6,3) * t254;
t147 = -t224 * t315 + t225 * t311;
t135 = rSges(6,1) * t191 - rSges(6,2) * t190 - rSges(6,3) * t246;
t102 = -t145 * t254 + t153 * t236;
t101 = t144 * t254 - t153 * t234;
t100 = -t144 * t236 + t145 * t234;
t92 = -t111 * t318 + t131 * t319;
t91 = t110 * t318 - t131 * t320;
t90 = -t145 * t314 - (-t153 + t356) * t315 + t364;
t89 = -t153 * t311 - (-t144 + t361) * t314 + t358;
t74 = -t110 * t319 + t111 * t320;
t73 = -t144 * t315 - (-t145 + t360) * t311 + t365;
t72 = -t121 * t314 - (-t135 + t357) * t315 + t366;
t71 = -t135 * t311 - (-t120 + t363) * t314 + t359;
t69 = t236 * t368 - t254 * t369;
t68 = -t234 * t368 + t254 * t370;
t64 = -t121 * t254 + t135 * t236 - t145 * t246 + t153 * t232;
t63 = t120 * t254 - t135 * t234 + t144 * t246 - t153 * t230;
t62 = t234 * t369 - t236 * t370;
t61 = -t369 * t314 - (t356 - t368) * t315 + t364;
t60 = -t368 * t311 - (t361 - t370) * t314 + t358;
t52 = -t370 * t315 - (t360 - t369) * t311 + t365;
t42 = -t378 * t314 - (t357 - t377) * t315 + t366;
t41 = -t377 * t311 - (t363 - t379) * t314 + t359;
t40 = t111 * t190 - t131 * t160 - t318 * t86 + t319 * t98;
t39 = -t110 * t190 + t131 * t158 + t318 * t85 - t320 * t98;
t36 = t232 * t368 + t236 * t377 - t246 * t369 - t254 * t378;
t35 = -t230 * t368 - t234 * t377 + t246 * t370 + t254 * t379;
t33 = -t311 * t58 - t314 * t70 - t315 * t59;
t29 = -t234 * t58 - t236 * t59 - t254 * t70;
t28 = -t318 * t70 - t319 * t59 - t320 * t58;
t27 = -t311 * t55 - t314 * t67 - t315 * t56;
t26 = -t311 * t53 - t314 * t66 - t315 * t54;
t25 = -t234 * t55 - t236 * t56 - t254 * t67;
t24 = -t234 * t53 - t236 * t54 - t254 * t66;
t23 = -t318 * t67 - t319 * t56 - t320 * t55;
t22 = -t318 * t66 - t319 * t54 - t320 * t53;
t15 = -t311 * t47 - t314 * t51 - t315 * t48;
t14 = -t311 * t45 - t314 * t50 - t315 * t46;
t13 = -t311 * t43 - t314 * t49 - t315 * t44;
t9 = -t20 * t311 - t21 * t315 - t314 * t34;
t8 = -t18 * t311 - t19 * t315 - t31 * t314;
t7 = -t16 * t311 - t17 * t315 - t30 * t314;
t5 = t158 * t58 + t160 * t59 + t190 * t70 - t20 * t320 - t21 * t319 - t318 * t34;
t2 = t158 * t55 + t160 * t56 - t18 * t320 - t19 * t319 + t190 * t67 - t31 * t318;
t1 = t158 * t53 - t16 * t320 + t160 * t54 - t17 * t319 + t190 * t66 - t30 * t318;
t10 = [0; 0; 0; m(5) * t355 / 0.2e1 - (m(4) * t224 + m(5) * t182 + m(6) * t120 + m(7) * t379) * t315 - (-m(4) * t225 - m(5) * t183 - m(6) * t121 - m(7) * t378) * t311 + (m(6) / 0.2e1 + t400) * (t185 * t401 - 0.2e1 * t162 + t355); -t72 * t324 - t42 * t323 + (m(4) * t147 + t317) * t300 + (t351 * t41 + t352 * t71 + (m(4) * t154 + t384) * t374) * t297 + (-m(4) * t155 - t383) * t328; (t38 * t52 + t41 * t60 + t42 * t61) * t338 + (t65 * t73 + t71 * t89 + t72 * t90) * t340 + (-t195 * t103 - t226 * t112 - t113 * t204) * t410 - (t9 + t15 - (t196 * t246 + t197 * t247 + t200 * t254 + t201 * t255 + t241 * t266 + t242 * t265 + t251 * t269 - t252 * t304 - t314 * t404) * t314 + ((-t172 - t212) * t112 + t173 * t113) * t410 + (-t154 * t209 + t155 * t210) * t411) * t314 - (t8 + t14 - (t169 * t232 + t171 * t233 + t179 * t236 + t181 * t237 + t206 * t259 + t208 * t258 + t221 * t262 - t223 * t308 - t315 * t402) * t315 + ((-t198 - t244) * t113 + t172 * t103) * t410 + (t147 * t209 - t155 * t243) * t411 - (t169 * t246 + t171 * t247 + t179 * t254 + t181 * t255 + t206 * t266 + t208 * t265 + t221 * t269 - t223 * t304 + t241 * t259 + t242 * t258 + t251 * t262 - t252 * t308 + t196 * t232 + t197 * t233 + t200 * t236 + t201 * t237 - t404 * t315 - t402 * t314) * t314) * t315 - (t7 + t13 - (t168 * t230 + t170 * t231 + t178 * t234 + t180 * t235 + t205 * t257 + t207 * t256 + t260 * t220 - t222 * t306 - t311 * t403) * t311 + (t198 * t112 + (-t173 - t213) * t103) * t410 + (-t147 * t210 + t154 * t243) * t411 - (t241 * t257 + t242 * t256 + t251 * t260 - t252 * t306 + t168 * t246 + t170 * t247 + t178 * t254 + t180 * t255 + t205 * t266 + t207 * t265 + t220 * t269 - t222 * t304 + t196 * t230 + t197 * t231 + t200 * t234 + t201 * t235 - t404 * t311 - t403 * t314) * t314 - (t206 * t257 + t208 * t256 + t221 * t260 - t223 * t306 + t169 * t230 + t171 * t231 + t179 * t234 + t181 * t235 + t205 * t259 + t207 * t258 + t220 * t262 - t222 * t308 + t168 * t232 + t170 * t233 + t178 * t236 + t180 * t237 - t402 * t311 - t403 * t315) * t315) * t311; 0; 0; -t317 * t314 - (m(6) * t71 + m(7) * t41 + t384) * t315 - (m(6) * t72 + m(7) * t42 + t383) * t311; 0; 0.2e1 * t32 * t400 + t412; -t64 * t324 - t36 * t323 + t333 * t300 + (t35 * t351 + t352 * t63) * t297; (t32 * t52 + t35 * t60 + t36 * t61 + t38 * t62 + t41 * t68 + t42 * t69) * m(7) + (t100 * t65 + t101 * t71 + t102 * t72 + t57 * t73 + t63 * t89 + t64 * t90) * m(6) + (-t9 / 0.2e1 - t15 / 0.2e1) * t254 + (-t8 / 0.2e1 - t14 / 0.2e1) * t236 + (-t7 / 0.2e1 - t13 / 0.2e1) * t234 + (-t311 * t75 - t314 * t93 - t315 * t76 + t26) * t394 + (-t311 * t77 - t314 * t94 - t315 * t78 + t27) * t393 + (-t311 * t87 - t314 * t99 - t315 * t88 + t33) * t391 + t409 * t390 + t408 * t389 + t407 * t388; -t333 * t314 - (m(6) * t63 + m(7) * t35) * t315 - (m(6) * t64 + m(7) * t36) * t311; (t32 * t62 + t35 * t68 + t36 * t69) * t338 - t246 * t29 + (t100 * t57 + t101 * t63 + t102 * t64) * t340 + (t380 - t407) * t254 + (t246 * t88 - t408) * t236 + (t246 * t87 - t409) * t234 + (t234 * t77 + t236 * t78 + t254 * t94 - t25) * t232 + (t234 * t75 + t236 * t76 + t254 * t93 - t24) * t230; t37 * m(7); (t37 * t300 - t328 * t40 + t347 * t39) * m(7); t1 * t390 + t2 * t389 + (t37 * t52 + t38 * t74 + t39 * t60 + t40 * t61 + t41 * t91 + t42 * t92) * m(7) + t33 * t397 + t9 * t392 + t5 * t388 + t26 * t399 + t7 * t396 + t27 * t398 + t8 * t395; (-t311 * t40 - t314 * t37 - t315 * t39) * m(7); (t32 * t74 + t35 * t91 + t36 * t92 + t37 * t62 + t39 * t68 + t40 * t69) * m(7) + t29 * t397 + t6 * t392 + t22 * t394 - t234 * t1 / 0.2e1 + t25 * t398 + t4 * t395 + t24 * t399 + t3 * t396 + t28 * t391 - t254 * t5 / 0.2e1 + t23 * t393 - t236 * t2 / 0.2e1; t160 * t23 - t319 * t2 + t158 * t22 - t320 * t1 + t190 * t28 - t318 * t5 + (t37 * t74 + t39 * t91 + t40 * t92) * t338;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t10(1) t10(2) t10(4) t10(7) t10(11) t10(16); t10(2) t10(3) t10(5) t10(8) t10(12) t10(17); t10(4) t10(5) t10(6) t10(9) t10(13) t10(18); t10(7) t10(8) t10(9) t10(10) t10(14) t10(19); t10(11) t10(12) t10(13) t10(14) t10(15) t10(20); t10(16) t10(17) t10(18) t10(19) t10(20) t10(21);];
Mq  = res;
