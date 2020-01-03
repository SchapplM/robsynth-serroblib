% Calculate time derivative of joint inertia matrix for
% S4RRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRR4_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR4_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR4_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR4_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:25:42
% EndTime: 2019-12-31 17:25:53
% DurationCPUTime: 6.46s
% Computational Cost: add. (15989->575), mult. (22534->827), div. (0->0), fcn. (21363->8), ass. (0->317)
t251 = -pkin(6) - pkin(5);
t246 = sin(qJ(2));
t336 = qJD(2) * t246;
t332 = pkin(2) * t336;
t421 = qJD(1) * t251 + t332;
t244 = qJ(2) + qJ(3);
t236 = sin(t244);
t241 = qJD(2) + qJD(3);
t248 = cos(qJ(4));
t245 = sin(qJ(4));
t381 = Icges(5,4) * t248;
t288 = -Icges(5,2) * t245 + t381;
t237 = cos(t244);
t362 = t237 * t241;
t382 = Icges(5,4) * t245;
t108 = t288 * t362 + (Icges(5,6) * t241 + (-Icges(5,2) * t248 - t382) * qJD(4)) * t236;
t167 = -Icges(5,6) * t237 + t236 * t288;
t292 = Icges(5,1) * t248 - t382;
t168 = -Icges(5,5) * t237 + t236 * t292;
t420 = -t245 * t108 + (-t167 * t248 - t168 * t245) * qJD(4);
t247 = sin(qJ(1));
t249 = cos(qJ(2));
t222 = rSges(3,1) * t246 + rSges(3,2) * t249;
t268 = qJD(2) * t222;
t419 = t247 * t268;
t250 = cos(qJ(1));
t385 = Icges(3,4) * t249;
t291 = -Icges(3,2) * t246 + t385;
t193 = Icges(3,6) * t247 + t250 * t291;
t386 = Icges(3,4) * t246;
t295 = Icges(3,1) * t249 - t386;
t195 = Icges(3,5) * t247 + t250 * t295;
t276 = t193 * t246 - t195 * t249;
t418 = t247 * t276;
t383 = Icges(4,4) * t237;
t289 = -Icges(4,2) * t236 + t383;
t177 = Icges(4,6) * t247 + t250 * t289;
t384 = Icges(4,4) * t236;
t293 = Icges(4,1) * t237 - t384;
t179 = Icges(4,5) * t247 + t250 * t293;
t278 = t177 * t236 - t179 * t237;
t417 = t247 * t278;
t233 = pkin(2) * t249 + pkin(1);
t394 = pkin(1) - t233;
t416 = t247 * t394;
t192 = -Icges(3,6) * t250 + t247 * t291;
t194 = -Icges(3,5) * t250 + t247 * t295;
t277 = t192 * t246 - t194 * t249;
t415 = t250 * t277;
t176 = -Icges(4,6) * t250 + t247 * t289;
t178 = -Icges(4,5) * t250 + t247 * t293;
t279 = t176 * t236 - t178 * t237;
t414 = t250 * t279;
t285 = Icges(5,5) * t248 - Icges(5,6) * t245;
t107 = t285 * t362 + (Icges(5,3) * t241 + (-Icges(5,5) * t245 - Icges(5,6) * t248) * qJD(4)) * t236;
t374 = t167 * t245;
t413 = -t241 * t374 - t107;
t352 = t248 * t250;
t355 = t245 * t247;
t203 = -t237 * t355 - t352;
t353 = t247 * t248;
t354 = t245 * t250;
t204 = t237 * t353 - t354;
t302 = -t204 * rSges(5,1) - t203 * rSges(5,2);
t364 = t236 * t247;
t141 = rSges(5,3) * t364 - t302;
t205 = -t237 * t354 + t353;
t206 = t237 * t352 + t355;
t363 = t236 * t250;
t142 = t206 * rSges(5,1) + t205 * rSges(5,2) + rSges(5,3) * t363;
t412 = -t247 * t141 - t250 * t142;
t208 = Icges(4,2) * t237 + t384;
t360 = t241 * t208;
t411 = t293 * t241 - t360;
t286 = Icges(4,5) * t237 - Icges(4,6) * t236;
t174 = -Icges(4,3) * t250 + t247 * t286;
t410 = qJD(1) * t174;
t287 = Icges(3,5) * t249 - Icges(3,6) * t246;
t190 = -Icges(3,3) * t250 + t247 * t287;
t308 = qJD(1) * t237 - qJD(4);
t357 = t241 * t250;
t326 = t236 * t357;
t409 = t247 * t308 + t326;
t209 = Icges(4,1) * t236 + t383;
t275 = t208 * t236 - t209 * t237;
t407 = qJD(1) * t275 + t286 * t241;
t406 = 2 * m(3);
t405 = 2 * m(4);
t404 = 2 * m(5);
t242 = t247 ^ 2;
t243 = t250 ^ 2;
t403 = t247 / 0.2e1;
t402 = -t250 / 0.2e1;
t401 = -rSges(5,3) - pkin(7);
t400 = m(3) * t222;
t210 = rSges(4,1) * t236 + rSges(4,2) * t237;
t399 = m(4) * t210;
t398 = pkin(2) * t246;
t397 = pkin(3) * t236;
t396 = pkin(3) * t237;
t395 = t247 * pkin(5);
t240 = t250 * pkin(5);
t393 = -pkin(5) - t251;
t392 = rSges(3,1) * t249;
t391 = rSges(4,1) * t237;
t390 = rSges(3,2) * t246;
t389 = rSges(3,3) * t250;
t135 = Icges(5,5) * t204 + Icges(5,6) * t203 + Icges(5,3) * t364;
t137 = Icges(5,4) * t204 + Icges(5,2) * t203 + Icges(5,6) * t364;
t139 = Icges(5,1) * t204 + Icges(5,4) * t203 + Icges(5,5) * t364;
t284 = -t137 * t245 + t139 * t248;
t309 = -qJD(4) * t237 + qJD(1);
t359 = t241 * t247;
t327 = t236 * t359;
t117 = t309 * t353 + (-t250 * t308 + t327) * t245;
t358 = t241 * t248;
t118 = t308 * t352 + (-t236 * t358 + t245 * t309) * t247;
t338 = qJD(1) * t250;
t259 = t236 * t338 + t237 * t359;
t69 = Icges(5,5) * t118 + Icges(5,6) * t117 + Icges(5,3) * t259;
t71 = Icges(5,4) * t118 + Icges(5,2) * t117 + Icges(5,6) * t259;
t73 = Icges(5,1) * t118 + Icges(5,4) * t117 + Icges(5,5) * t259;
t19 = (t241 * t284 - t69) * t237 + (t135 * t241 - t245 * t71 + t248 * t73 + (-t137 * t248 - t139 * t245) * qJD(4)) * t236;
t388 = t19 * t250;
t136 = Icges(5,5) * t206 + Icges(5,6) * t205 + Icges(5,3) * t363;
t138 = Icges(5,4) * t206 + Icges(5,2) * t205 + Icges(5,6) * t363;
t140 = Icges(5,1) * t206 + Icges(5,4) * t205 + Icges(5,5) * t363;
t283 = -t138 * t245 + t140 * t248;
t274 = t309 * t250;
t115 = t409 * t245 + t248 * t274;
t116 = t245 * t274 - t409 * t248;
t339 = qJD(1) * t247;
t319 = t236 * t339;
t325 = t237 * t357;
t258 = -t319 + t325;
t68 = Icges(5,5) * t116 + Icges(5,6) * t115 + Icges(5,3) * t258;
t70 = Icges(5,4) * t116 + Icges(5,2) * t115 + Icges(5,6) * t258;
t72 = Icges(5,1) * t116 + Icges(5,4) * t115 + Icges(5,5) * t258;
t20 = (t241 * t283 - t68) * t237 + (t136 * t241 - t245 * t70 + t248 * t72 + (-t138 * t248 - t140 * t245) * qJD(4)) * t236;
t387 = t20 * t247;
t239 = t247 * rSges(3,3);
t238 = t247 * rSges(4,3);
t304 = -rSges(4,2) * t236 + t391;
t188 = t304 * t241;
t371 = t188 * t247;
t370 = t192 * t249;
t369 = t193 * t249;
t368 = t194 * t246;
t367 = t195 * t246;
t366 = t209 * t241;
t365 = t236 * t241;
t361 = t237 * t250;
t351 = t250 * t251;
t301 = rSges(5,1) * t248 - rSges(5,2) * t245;
t110 = t301 * t362 + (rSges(5,3) * t241 + (-rSges(5,1) * t245 - rSges(5,2) * t248) * qJD(4)) * t236;
t306 = pkin(7) * t236 + t396;
t350 = -t306 * t241 - t110;
t201 = pkin(3) * t361 + pkin(7) * t363;
t349 = -t142 - t201;
t169 = -rSges(5,3) * t237 + t236 * t301;
t157 = t169 * t339;
t211 = -pkin(7) * t237 + t397;
t348 = t211 * t339 + t157;
t172 = t240 + t351 - t416;
t223 = t250 * t233;
t173 = -t250 * pkin(1) + t247 * t393 + t223;
t347 = t247 * t172 + t250 * t173;
t180 = -t250 * rSges(4,3) + t247 * t304;
t181 = rSges(4,1) * t361 - rSges(4,2) * t363 + t238;
t119 = t247 * t180 + t250 * t181;
t346 = -t169 - t211;
t345 = rSges(4,2) * t319 + rSges(4,3) * t338;
t344 = t421 * t247;
t343 = t250 * t392 + t239;
t342 = t242 + t243;
t175 = Icges(4,3) * t247 + t250 * t286;
t341 = qJD(1) * t175;
t191 = Icges(3,3) * t247 + t250 * t287;
t340 = qJD(1) * t191;
t335 = qJD(2) * t249;
t333 = t250 * t390;
t331 = pkin(2) * t335;
t61 = -t135 * t237 + t236 * t284;
t166 = -Icges(5,3) * t237 + t236 * t285;
t80 = t166 * t364 + t167 * t203 + t168 * t204;
t330 = t61 / 0.2e1 + t80 / 0.2e1;
t62 = -t136 * t237 + t236 * t283;
t81 = t166 * t363 + t167 * t205 + t168 * t206;
t329 = t62 / 0.2e1 + t81 / 0.2e1;
t109 = t292 * t362 + (Icges(5,5) * t241 + (-Icges(5,1) * t245 - t381) * qJD(4)) * t236;
t324 = t236 * t248 * t109 + t237 * t168 * t358 + t166 * t365;
t260 = -t237 * t339 - t326;
t269 = t210 * t241;
t323 = t247 * (-t247 * t269 + (t250 * t304 + t238) * qJD(1)) + t250 * (rSges(4,1) * t260 - rSges(4,2) * t325 + t345) + t180 * t338;
t322 = t116 * rSges(5,1) + t115 * rSges(5,2) + rSges(5,3) * t325;
t321 = t247 * ((-t250 * t394 - t395) * qJD(1) - t344) + t250 * (-t250 * t332 + (t250 * t393 + t416) * qJD(1)) + t172 * t338;
t320 = t246 * t339;
t317 = t339 / 0.2e1;
t316 = t338 / 0.2e1;
t315 = -t210 - t398;
t144 = t346 * t250;
t122 = -t176 * qJD(1) - t250 * t360;
t314 = t179 * t241 + t122;
t123 = qJD(1) * t177 - t247 * t360;
t313 = t178 * t241 + t123;
t124 = -t178 * qJD(1) - t250 * t366;
t312 = -t177 * t241 + t124;
t125 = qJD(1) * t179 - t247 * t366;
t311 = t176 * t241 - t125;
t310 = -t247 * t251 + t223;
t200 = t306 * t247;
t65 = t247 * t200 + t250 * t201 - t412;
t307 = t346 - t398;
t305 = -t390 + t392;
t303 = rSges(5,1) * t118 + rSges(5,2) * t117;
t52 = t135 * t364 + t137 * t203 + t139 * t204;
t53 = t136 * t364 + t138 * t203 + t140 * t204;
t38 = t247 * t53 - t250 * t52;
t300 = t247 * t52 + t250 * t53;
t54 = t135 * t363 + t137 * t205 + t139 * t206;
t55 = t136 * t363 + t138 * t205 + t140 * t206;
t39 = t247 * t55 - t250 * t54;
t299 = t247 * t54 + t250 * t55;
t298 = t62 * t247 - t61 * t250;
t297 = t61 * t247 + t62 * t250;
t207 = Icges(4,5) * t236 + Icges(4,6) * t237;
t265 = t241 * t207;
t120 = -t250 * t265 - t410;
t121 = -t247 * t265 + t341;
t13 = t115 * t137 + t116 * t139 + t135 * t258 + t205 * t71 + t206 * t73 + t363 * t69;
t14 = t115 * t138 + t116 * t140 + t136 * t258 + t205 * t70 + t206 * t72 + t363 * t68;
t8 = qJD(1) * t299 - t13 * t250 + t14 * t247;
t86 = -t174 * t250 - t247 * t279;
t87 = -t175 * t250 - t417;
t88 = t174 * t247 - t414;
t89 = t175 * t247 - t250 * t278;
t296 = t39 * t338 + t38 * t339 + (-t88 * t338 - t86 * t339) * t250 + (t8 + (t89 * qJD(1) + (t123 * t236 - t125 * t237 + t176 * t362 + t178 * t365 - t410) * t250) * t250 + t87 * t339 + t89 * t338 + ((t88 + t417) * qJD(1) + (-t121 + t312 * t237 - t314 * t236 + (t175 - t279) * qJD(1)) * t250 + t247 * t120) * t247) * t247;
t294 = Icges(3,1) * t246 + t385;
t290 = Icges(3,2) * t249 + t386;
t282 = t141 * t250 - t142 * t247;
t273 = -t331 + t350;
t272 = -pkin(1) - t305;
t217 = pkin(7) * t325;
t74 = -rSges(5,3) * t319 + t322;
t75 = rSges(5,3) * t259 + t303;
t271 = (t141 + t200) * t338 + (pkin(3) * t260 - pkin(7) * t319 + t217 + t74) * t250 + (t75 + t259 * pkin(7) + (t237 * t338 - t327) * pkin(3)) * t247;
t128 = t307 * t250;
t270 = -t233 - t304;
t264 = qJD(2) * t294;
t263 = qJD(2) * t290;
t262 = qJD(2) * (-Icges(3,5) * t246 - Icges(3,6) * t249);
t261 = t236 * t401 - t233 - t396;
t12 = (t250 * t121 + (t87 + t414) * qJD(1)) * t250 + (t86 * qJD(1) + (-t122 * t236 + t124 * t237 - t177 * t362 - t179 * t365 + t341) * t247 + (-t120 + t311 * t237 + t313 * t236 + (-t174 - t278) * qJD(1)) * t250) * t247;
t15 = t117 * t137 + t118 * t139 + t135 * t259 + t203 * t71 + t204 * t73 + t364 * t69;
t16 = t117 * t138 + t118 * t140 + t136 * t259 + t203 * t70 + t204 * t72 + t364 * t68;
t9 = qJD(1) * t300 - t15 * t250 + t16 * t247;
t257 = (-t12 - t9) * t250 + t296;
t25 = t236 * t300 - t237 * t80;
t26 = t236 * t299 - t237 * t81;
t30 = t107 * t363 + t108 * t205 + t109 * t206 + t115 * t167 + t116 * t168 + t166 * t258;
t3 = (t241 * t299 - t30) * t237 + (-qJD(1) * t39 + t13 * t247 + t14 * t250 + t241 * t81) * t236;
t31 = t107 * t364 + t108 * t203 + t109 * t204 + t117 * t167 + t118 * t168 + t166 * t259;
t4 = (t241 * t300 - t31) * t237 + (-qJD(1) * t38 + t15 * t247 + t16 * t250 + t241 * t80) * t236;
t256 = t3 * t403 + t9 * t364 / 0.2e1 + t4 * t402 - t237 * (qJD(1) * t297 + t387 - t388) / 0.2e1 + t25 * t317 - t39 * t319 / 0.2e1 + t298 * t365 / 0.2e1 + t8 * t363 / 0.2e1 + (t247 * t38 + t250 * t39) * t362 / 0.2e1 + (t236 * t38 + t26) * t316;
t255 = rSges(3,2) * t320 + rSges(3,3) * t338 - t250 * t268;
t254 = t247 * t261 - t351;
t186 = t289 * t241;
t253 = qJD(1) * t207 + t411 * t237 + (-t186 - t366) * t236;
t252 = -t388 / 0.2e1 + t387 / 0.2e1 + (t236 * t312 + t237 * t314 + t407 * t247 + t253 * t250 + t30) * t403 + (-t236 * t311 + t237 * t313 + t253 * t247 - t407 * t250 + t31) * t402 + (t176 * t237 + t178 * t236 - t207 * t250 - t247 * t275 + t61 + t80) * t317 + (t177 * t237 + t179 * t236 + t207 * t247 - t250 * t275 + t62 + t81) * t316;
t230 = pkin(2) * t320;
t216 = t305 * qJD(2);
t197 = -t333 + t343;
t196 = t247 * t305 - t389;
t171 = t315 * t250;
t170 = t315 * t247;
t161 = t395 + (pkin(1) - t390) * t250 + t343;
t160 = t247 * t272 + t240 + t389;
t154 = t181 + t310;
t153 = (rSges(4,3) - t251) * t250 + t270 * t247;
t148 = t247 * t262 + t340;
t147 = -qJD(1) * t190 + t250 * t262;
t143 = t346 * t247;
t130 = t419 + ((-rSges(3,3) - pkin(5)) * t247 + t272 * t250) * qJD(1);
t129 = (t240 + (-pkin(1) - t392) * t247) * qJD(1) + t255;
t127 = t307 * t247;
t106 = -t210 * t338 - t371 + (-t246 * t338 - t247 * t335) * pkin(2);
t105 = t210 * t339 + t230 + (-t188 - t331) * t250;
t99 = t191 * t247 - t250 * t276;
t98 = t190 * t247 - t415;
t97 = -t191 * t250 - t418;
t96 = -t190 * t250 - t247 * t277;
t95 = t210 * t359 + (t250 * t270 - t238) * qJD(1) + t344;
t94 = (-t233 - t391) * t339 + (-t269 - t421) * t250 + t345;
t93 = t310 - t349;
t92 = t254 + t302;
t91 = -t142 * t237 - t169 * t363;
t90 = t141 * t237 + t169 * t364;
t85 = -t166 * t237 + (t168 * t248 - t374) * t236;
t84 = t85 * t365;
t83 = t119 + t347;
t82 = t282 * t236;
t77 = qJD(1) * t144 + t247 * t350;
t76 = t250 * t350 + t348;
t64 = qJD(1) * t128 + t247 * t273;
t63 = t250 * t273 + t230 + t348;
t60 = -t181 * t339 + t323;
t49 = t65 + t347;
t46 = (t237 * t401 + t397) * t359 + t261 * t338 - t303 + t344;
t45 = t217 + (-pkin(3) * t365 - t332) * t250 + t254 * qJD(1) + t322;
t44 = (-t173 - t181) * t339 + t321 + t323;
t43 = (t169 * t359 + t75) * t237 + (t110 * t247 - t141 * t241 + t169 * t338) * t236;
t42 = (-t169 * t357 - t74) * t237 + (-t110 * t250 + t142 * t241 + t157) * t236;
t40 = t420 * t236 + t413 * t237 + t324;
t27 = t282 * t362 + (t412 * qJD(1) - t247 * t74 + t250 * t75) * t236;
t24 = t339 * t349 + t271;
t21 = (-t173 + t349) * t339 + t271 + t321;
t1 = [(t45 * t93 + t46 * t92) * t404 + (t153 * t95 + t154 * t94) * t405 + t209 * t362 + (t129 * t161 + t130 * t160) * t406 + t324 + (t295 - t290) * t336 + (t294 + t291) * t335 + (t186 + t413) * t237 + (t411 + t420) * t236; (t242 / 0.2e1 + t243 / 0.2e1) * t287 * qJD(2) + m(3) * ((-t129 * t247 - t130 * t250) * t222 + (-t160 * t250 - t161 * t247) * t216) + t252 + ((-t161 * t400 + t369 / 0.2e1 + t367 / 0.2e1) * t250 + (t370 / 0.2e1 + t368 / 0.2e1 + t160 * t400) * t247) * qJD(1) + m(5) * (t127 * t45 + t128 * t46 + t63 * t92 + t64 * t93) + m(4) * (t105 * t153 + t106 * t154 + t170 * t94 + t171 * t95) + (-qJD(2) * t276 + (-t192 * qJD(1) - t250 * t263) * t249 + (-t194 * qJD(1) - t250 * t264) * t246) * t403 + (-qJD(2) * t277 + (qJD(1) * t193 - t247 * t263) * t249 + (qJD(1) * t195 - t247 * t264) * t246) * t402; -t250 * t9 + (t127 * t64 + t128 * t63 + t49 * t21) * t404 - t250 * t12 + (t105 * t171 + t106 * t170 + t44 * t83) * t405 + ((t196 * t247 + t197 * t250) * ((qJD(1) * t196 + t255) * t250 + (-t419 + (-t197 - t333 + t239) * qJD(1)) * t247) + t342 * t222 * t216) * t406 + (t247 * t97 - t250 * t96) * t339 - t250 * ((t250 * t148 + (t97 + t415) * qJD(1)) * t250 + (t96 * qJD(1) + (-t193 * t335 - t195 * t336 + t340) * t247 + (-t147 + (t368 + t370) * qJD(2) - t276 * qJD(1)) * t250) * t247) + (t247 * t99 - t250 * t98) * t338 + t247 * ((t247 * t147 + (t98 + t418) * qJD(1)) * t247 + (t99 * qJD(1) + (t192 * t335 + t194 * t336) * t250 + (-t148 + (-t367 - t369) * qJD(2) + (t191 - t277) * qJD(1)) * t247) * t250) + t296; (-t247 * t94 - t250 * t95 + (t153 * t247 - t154 * t250) * qJD(1)) * t399 + m(4) * (-t153 * t250 - t154 * t247) * t188 + t252 + m(5) * (t143 * t45 + t144 * t46 + t76 * t92 + t77 * t93); m(5) * (t127 * t77 + t128 * t76 + t143 * t64 + t144 * t63 + t21 * t65 + t24 * t49) + m(4) * (-t171 * t188 * t250 + t119 * t44 - t170 * t371 + t60 * t83) + (-t105 * t250 - t106 * t247 + (-t170 * t250 + t171 * t247) * qJD(1)) * t399 + t257; (t188 * t210 * t342 + t119 * t60) * t405 + (t143 * t77 + t144 * t76 + t24 * t65) * t404 + t257; m(5) * (t42 * t93 + t43 * t92 + t45 * t91 + t46 * t90) + t84 + (-t40 + (t247 * t330 + t250 * t329) * t241) * t237 + ((t20 / 0.2e1 + t30 / 0.2e1) * t250 + (t19 / 0.2e1 + t31 / 0.2e1) * t247 + (-t247 * t329 + t250 * t330) * qJD(1)) * t236; m(5) * (t127 * t42 + t128 * t43 + t21 * t82 + t27 * t49 + t63 * t90 + t64 * t91) + t256; m(5) * (t143 * t42 + t144 * t43 + t24 * t82 + t27 * t65 + t76 * t90 + t77 * t91) + t256; (t27 * t82 + t42 * t91 + t43 * t90) * t404 + (t40 * t237 - t84 + (-t237 * t297 + t247 * t25 + t250 * t26) * t241) * t237 + (t250 * t3 + t247 * t4 + t297 * t365 + (-t19 * t247 - t20 * t250 - t241 * t85) * t237 + (t237 * t298 - t247 * t26 + t250 * t25) * qJD(1)) * t236;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
