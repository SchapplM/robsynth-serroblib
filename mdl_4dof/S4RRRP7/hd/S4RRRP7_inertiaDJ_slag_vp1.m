% Calculate time derivative of joint inertia matrix for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP7_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP7_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:10
% EndTime: 2019-12-31 17:20:24
% DurationCPUTime: 8.36s
% Computational Cost: add. (9193->634), mult. (25166->908), div. (0->0), fcn. (24492->6), ass. (0->314)
t231 = sin(qJ(2));
t234 = cos(qJ(2));
t230 = sin(qJ(3));
t233 = cos(qJ(3));
t268 = Icges(4,5) * t233 - Icges(4,6) * t230;
t167 = -Icges(4,3) * t234 + t231 * t268;
t270 = Icges(5,4) * t233 + Icges(5,6) * t230;
t170 = -Icges(5,2) * t234 + t231 * t270;
t404 = t167 + t170;
t235 = cos(qJ(1));
t337 = t233 * t235;
t232 = sin(qJ(1));
t338 = t232 * t234;
t190 = t230 * t338 + t337;
t191 = -t230 * t235 + t233 * t338;
t398 = rSges(5,3) + qJ(4);
t400 = rSges(5,1) + pkin(3);
t403 = t398 * t190 + t400 * t191;
t321 = qJD(2) * t232;
t301 = t231 * t321;
t316 = qJD(3) * t233;
t318 = qJD(3) * t230;
t323 = qJD(1) * t235;
t325 = qJD(1) * t232;
t119 = -t230 * t301 - t235 * t318 - t233 * t325 + (t230 * t323 + t232 * t316) * t234;
t256 = (-qJD(3) * t234 + qJD(1)) * t230;
t324 = qJD(1) * t234;
t290 = -qJD(3) + t324;
t322 = qJD(2) * t231;
t120 = t290 * t337 + (-t233 * t322 + t256) * t232;
t402 = -t190 * qJD(4) - t398 * t119 - t400 * t120;
t401 = -qJD(1) * t231 / 0.2e1;
t354 = Icges(5,5) * t233;
t267 = Icges(5,3) * t230 + t354;
t166 = -Icges(5,6) * t234 + t231 * t267;
t357 = Icges(4,4) * t233;
t271 = -Icges(4,2) * t230 + t357;
t171 = -Icges(4,6) * t234 + t231 * t271;
t355 = Icges(5,5) * t230;
t274 = Icges(5,1) * t233 + t355;
t174 = -Icges(5,4) * t234 + t231 * t274;
t358 = Icges(4,4) * t230;
t275 = Icges(4,1) * t233 - t358;
t175 = -Icges(4,5) * t234 + t231 * t275;
t399 = t404 * t234 + ((-t174 - t175) * t233 + (-t166 + t171) * t230) * t231;
t297 = t231 * t316;
t320 = qJD(2) * t234;
t302 = t230 * t320;
t397 = t297 + t302;
t336 = t234 * t235;
t192 = t230 * t336 - t232 * t233;
t311 = t233 * t336;
t193 = t232 * t230 + t311;
t339 = t231 * t235;
t126 = Icges(5,4) * t193 + Icges(5,2) * t339 + Icges(5,6) * t192;
t122 = Icges(5,5) * t193 + Icges(5,6) * t339 + Icges(5,3) * t192;
t130 = Icges(5,1) * t193 + Icges(5,4) * t339 + Icges(5,5) * t192;
t265 = t122 * t230 + t130 * t233;
t60 = -t126 * t234 + t231 * t265;
t124 = Icges(4,5) * t193 - Icges(4,6) * t192 + Icges(4,3) * t339;
t128 = Icges(4,4) * t193 - Icges(4,2) * t192 + Icges(4,6) * t339;
t132 = Icges(4,1) * t193 - Icges(4,4) * t192 + Icges(4,5) * t339;
t263 = -t128 * t230 + t132 * t233;
t62 = -t124 * t234 + t231 * t263;
t366 = t60 + t62;
t341 = t231 * t232;
t125 = Icges(5,4) * t191 + Icges(5,2) * t341 + Icges(5,6) * t190;
t121 = Icges(5,5) * t191 + Icges(5,6) * t341 + Icges(5,3) * t190;
t129 = Icges(5,1) * t191 + Icges(5,4) * t341 + Icges(5,5) * t190;
t266 = t121 * t230 + t129 * t233;
t59 = -t125 * t234 + t231 * t266;
t123 = Icges(4,5) * t191 - Icges(4,6) * t190 + Icges(4,3) * t341;
t127 = Icges(4,4) * t191 - Icges(4,2) * t190 + Icges(4,6) * t341;
t131 = Icges(4,1) * t191 - Icges(4,4) * t190 + Icges(4,5) * t341;
t264 = -t127 * t230 + t131 * t233;
t61 = -t123 * t234 + t231 * t264;
t367 = t59 + t61;
t396 = t232 * t367 + t235 * t366;
t296 = t232 * t320;
t241 = t231 * t323 + t296;
t67 = Icges(5,5) * t120 + Icges(5,6) * t241 + Icges(5,3) * t119;
t71 = Icges(5,4) * t120 + Icges(5,2) * t241 + Icges(5,6) * t119;
t75 = Icges(5,1) * t120 + Icges(5,4) * t241 + Icges(5,5) * t119;
t19 = (qJD(2) * t266 - t71) * t234 + (qJD(2) * t125 + t230 * t67 + t233 * t75 + (t121 * t233 - t129 * t230) * qJD(3)) * t231;
t69 = Icges(4,5) * t120 - Icges(4,6) * t119 + Icges(4,3) * t241;
t73 = Icges(4,4) * t120 - Icges(4,2) * t119 + Icges(4,6) * t241;
t77 = Icges(4,1) * t120 - Icges(4,4) * t119 + Icges(4,5) * t241;
t21 = (qJD(2) * t264 - t69) * t234 + (qJD(2) * t123 - t230 * t73 + t233 * t77 + (-t127 * t233 - t131 * t230) * qJD(3)) * t231;
t395 = -t19 - t21;
t319 = qJD(2) * t235;
t299 = t231 * t319;
t117 = qJD(1) * t190 - qJD(3) * t311 + t230 * t299 - t232 * t318;
t118 = t235 * t256 + (-t232 * t290 - t299) * t233;
t295 = t234 * t319;
t303 = t231 * t325;
t240 = t295 - t303;
t66 = Icges(5,5) * t118 + Icges(5,6) * t240 - Icges(5,3) * t117;
t70 = Icges(5,4) * t118 + Icges(5,2) * t240 - Icges(5,6) * t117;
t74 = Icges(5,1) * t118 + Icges(5,4) * t240 - Icges(5,5) * t117;
t20 = (qJD(2) * t265 - t70) * t234 + (qJD(2) * t126 + t230 * t66 + t233 * t74 + (t122 * t233 - t130 * t230) * qJD(3)) * t231;
t68 = Icges(4,5) * t118 + Icges(4,6) * t117 + Icges(4,3) * t240;
t72 = Icges(4,4) * t118 + Icges(4,2) * t117 + Icges(4,6) * t240;
t76 = Icges(4,1) * t118 + Icges(4,4) * t117 + Icges(4,5) * t240;
t22 = (qJD(2) * t263 - t68) * t234 + (qJD(2) * t124 - t230 * t72 + t233 * t76 + (-t128 * t233 - t132 * t230) * qJD(3)) * t231;
t394 = t20 + t22;
t49 = t123 * t341 - t127 * t190 + t131 * t191;
t50 = t124 * t341 - t128 * t190 + t132 * t191;
t280 = t232 * t49 + t235 * t50;
t47 = t121 * t190 + t125 * t341 + t129 * t191;
t48 = t122 * t190 + t126 * t341 + t130 * t191;
t281 = t232 * t47 + t235 * t48;
t86 = t166 * t190 + t170 * t341 + t174 * t191;
t87 = t167 * t341 - t171 * t190 + t175 * t191;
t393 = (-t86 - t87) * t234 + (t280 + t281) * t231;
t53 = t123 * t339 - t192 * t127 + t193 * t131;
t54 = t124 * t339 - t192 * t128 + t193 * t132;
t278 = t232 * t53 + t235 * t54;
t51 = t192 * t121 + t125 * t339 + t193 * t129;
t52 = t192 * t122 + t126 * t339 + t193 * t130;
t279 = t232 * t51 + t235 * t52;
t88 = t192 * t166 + t170 * t339 + t193 * t174;
t89 = t167 * t339 - t192 * t171 + t193 * t175;
t392 = (-t88 - t89) * t234 + (t278 + t279) * t231;
t317 = qJD(3) * t231;
t138 = (Icges(5,3) * t233 - t355) * t317 + (Icges(5,6) * t231 + t234 * t267) * qJD(2);
t139 = (-Icges(4,5) * t230 - Icges(4,6) * t233) * t317 + (Icges(4,3) * t231 + t234 * t268) * qJD(2);
t142 = (-Icges(5,4) * t230 + Icges(5,6) * t233) * t317 + (Icges(5,2) * t231 + t234 * t270) * qJD(2);
t146 = (-Icges(5,1) * t230 + t354) * t317 + (Icges(5,4) * t231 + t234 * t274) * qJD(2);
t147 = (-Icges(4,1) * t230 - t357) * t317 + (Icges(4,5) * t231 + t234 * t275) * qJD(2);
t298 = t230 * t317;
t300 = t233 * t320;
t342 = t230 * t231;
t391 = -t171 * t297 + t175 * t300 + t138 * t342 + (-t298 + t300) * t174 + t397 * t166 + (t147 + t146) * t231 * t233 + t404 * t322 + (-t139 - t142) * t234;
t390 = rSges(5,2) * t295 + t192 * qJD(4) - t398 * t117 + t400 * t118;
t360 = Icges(3,4) * t231;
t277 = Icges(3,1) * t234 - t360;
t177 = Icges(3,5) * t232 + t235 * t277;
t343 = t177 * t234;
t359 = Icges(3,4) * t234;
t273 = -Icges(3,2) * t231 + t359;
t173 = Icges(3,6) * t232 + t235 * t273;
t348 = t173 * t231;
t257 = -t343 + t348;
t389 = t232 * t257;
t176 = -Icges(3,5) * t235 + t232 * t277;
t345 = t176 * t234;
t172 = -Icges(3,6) * t235 + t232 * t273;
t350 = t172 * t231;
t258 = -t345 + t350;
t388 = t235 * t258;
t387 = -rSges(3,2) * t339 + t232 * rSges(3,3);
t269 = Icges(3,5) * t234 - Icges(3,6) * t231;
t168 = -Icges(3,3) * t235 + t232 * t269;
t386 = t234 * t366 - t392;
t333 = rSges(5,2) * t339 + t398 * t192 + t400 * t193;
t334 = rSges(5,2) * t341 + t403;
t385 = -t232 * t334 - t235 * t333;
t384 = 2 * m(3);
t383 = 2 * m(4);
t382 = 2 * m(5);
t381 = t232 ^ 2;
t380 = t235 ^ 2;
t378 = -t234 / 0.2e1;
t375 = -rSges(5,2) - pkin(6);
t374 = -rSges(4,3) - pkin(6);
t208 = rSges(3,1) * t231 + rSges(3,2) * t234;
t373 = m(3) * t208;
t372 = pkin(2) * t234;
t227 = t232 * pkin(5);
t143 = (-Icges(4,2) * t233 - t358) * t317 + (Icges(4,6) * t231 + t234 * t271) * qJD(2);
t368 = (-t171 * t320 + (-qJD(3) * t175 - t143) * t231) * t230 + t391;
t365 = -rSges(5,2) * t303 + t390;
t364 = rSges(5,2) * t241 - t402;
t363 = t399 * t322;
t362 = rSges(3,3) * t235;
t284 = -rSges(4,1) * t191 + rSges(4,2) * t190;
t135 = rSges(4,3) * t341 - t284;
t351 = t135 * t235;
t349 = t172 * t234;
t347 = t173 * t234;
t346 = t176 * t231;
t344 = t177 * t231;
t282 = rSges(5,1) * t233 + rSges(5,3) * t230;
t335 = (pkin(3) * t320 + qJ(4) * t317) * t233 + (qJ(4) * t320 + (-pkin(3) * qJD(3) + qJD(4)) * t231) * t230 + (-rSges(5,1) * t230 + rSges(5,3) * t233) * t317 + (rSges(5,2) * t231 + t234 * t282) * qJD(2);
t283 = rSges(4,1) * t233 - rSges(4,2) * t230;
t151 = (-rSges(4,1) * t230 - rSges(4,2) * t233) * t317 + (rSges(4,3) * t231 + t234 * t283) * qJD(2);
t288 = pkin(6) * t231 + t372;
t202 = t288 * qJD(2);
t332 = -t151 - t202;
t331 = -rSges(5,2) * t234 + (pkin(3) * t233 + qJ(4) * t230 + t282) * t231;
t179 = -rSges(4,3) * t234 + t231 * t283;
t209 = pkin(2) * t231 - pkin(6) * t234;
t330 = -t179 - t209;
t195 = t288 * t232;
t222 = pkin(2) * t336;
t196 = pkin(6) * t339 + t222;
t329 = t232 * t195 + t235 * t196;
t328 = rSges(3,2) * t303 + rSges(3,3) * t323;
t327 = t235 * pkin(1) + t227;
t169 = Icges(3,3) * t232 + t235 * t269;
t326 = qJD(1) * t169;
t35 = t48 * t232 - t235 * t47;
t36 = t50 * t232 - t235 * t49;
t313 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t52 * t232 - t235 * t51;
t38 = t54 * t232 - t235 * t53;
t312 = t37 / 0.2e1 + t38 / 0.2e1;
t310 = t118 * rSges(4,1) + t117 * rSges(4,2) + rSges(4,3) * t295;
t308 = -t202 - t335;
t213 = pkin(2) * t301;
t214 = pkin(6) * t295;
t242 = -t232 * t324 - t299;
t307 = t232 * (pkin(6) * t241 + qJD(1) * t222 - t213) + t235 * (pkin(2) * t242 - pkin(6) * t303 + t214) + t195 * t323;
t306 = -t209 - t331;
t137 = t193 * rSges(4,1) - t192 * rSges(4,2) + rSges(4,3) * t339;
t305 = -pkin(1) - t372;
t304 = t179 * t325;
t294 = t235 * t331;
t293 = t334 * t235;
t292 = t331 * t232;
t155 = t330 * t235;
t291 = qJD(1) * t331;
t289 = t327 + t196;
t106 = t306 * t235;
t287 = t232 * t291;
t286 = rSges(3,1) * t234 - rSges(3,2) * t231;
t285 = t120 * rSges(4,1) - t119 * rSges(4,2);
t272 = Icges(3,2) * t234 + t360;
t262 = -t137 * t232 + t351;
t261 = -t232 * t135 - t137 * t235;
t181 = rSges(3,1) * t336 + t387;
t255 = -t234 * t367 + t393;
t254 = -pkin(1) - t286;
t33 = t119 * t166 + t120 * t174 + t190 * t138 + t142 * t341 + t191 * t146 + t170 * t241;
t34 = -t119 * t171 + t120 * t175 + t139 * t341 - t190 * t143 + t191 * t147 + t167 * t241;
t253 = t21 / 0.2e1 + t19 / 0.2e1 + t33 / 0.2e1 + t34 / 0.2e1;
t31 = -t117 * t166 + t118 * t174 + t192 * t138 + t142 * t339 + t193 * t146 + t170 * t240;
t32 = t117 * t171 + t118 * t175 + t139 * t339 - t192 * t143 + t193 * t147 + t167 * t240;
t252 = t31 / 0.2e1 + t32 / 0.2e1 + t22 / 0.2e1 + t20 / 0.2e1;
t251 = t61 / 0.2e1 + t59 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t250 = t62 / 0.2e1 + t60 / 0.2e1 + t88 / 0.2e1 + t89 / 0.2e1;
t225 = pkin(5) * t323;
t249 = -pkin(2) * t299 + t214 + t225;
t248 = qJD(2) * t208;
t247 = t231 * t375 + t305;
t246 = t231 * t374 + t305;
t244 = qJD(2) * t272;
t243 = qJD(2) * (-Icges(3,5) * t231 - Icges(3,6) * t234);
t239 = t247 * t232;
t238 = t246 * t232;
t237 = -t232 * t333 + t293;
t228 = t235 * pkin(5);
t201 = t286 * qJD(2);
t197 = t209 * t325;
t180 = t232 * t286 - t362;
t159 = t181 + t327;
t158 = t232 * t254 + t228 + t362;
t154 = t330 * t232;
t141 = t232 * t243 + t326;
t140 = -qJD(1) * t168 + t235 * t243;
t108 = t208 * t321 + ((-rSges(3,3) - pkin(5)) * t232 + t254 * t235) * qJD(1);
t107 = rSges(3,1) * t242 - rSges(3,2) * t295 - pkin(1) * t325 + t225 + t328;
t105 = t306 * t232;
t104 = t289 + t137;
t103 = t228 + t238 + t284;
t102 = -t234 * t137 - t179 * t339;
t101 = t135 * t234 + t179 * t341;
t100 = t232 * t169 - t257 * t235;
t99 = t232 * t168 - t388;
t98 = -t169 * t235 - t389;
t97 = -t168 * t235 - t232 * t258;
t92 = t262 * t231;
t91 = qJD(1) * t155 + t232 * t332;
t90 = t235 * t332 + t197 + t304;
t85 = t289 + t333;
t84 = t228 + t239 - t403;
t83 = rSges(4,3) * t241 + t285;
t81 = -rSges(4,3) * t303 + t310;
t65 = -t261 + t329;
t64 = -t231 * t294 - t234 * t333;
t63 = t231 * t292 + t234 * t334;
t58 = t213 + t374 * t296 + (t235 * t246 - t227) * qJD(1) - t285;
t57 = qJD(1) * t238 + t249 + t310;
t56 = qJD(1) * t106 + t232 * t308;
t55 = t235 * t308 + t197 + t287;
t46 = t237 * t231;
t45 = t329 - t385;
t44 = (t179 * t321 + t83) * t234 + (-qJD(2) * t135 + t232 * t151 + t179 * t323) * t231;
t43 = (-t179 * t319 - t81) * t234 + (qJD(2) * t137 - t151 * t235 + t304) * t231;
t42 = t213 + t375 * t296 + (t235 * t247 - t227) * qJD(1) + t402;
t41 = qJD(1) * t239 + t249 + t390;
t30 = t262 * t320 + (qJD(1) * t261 - t232 * t81 + t235 * t83) * t231;
t29 = t232 * t83 + t235 * t81 + (t351 + (-t137 - t196) * t232) * qJD(1) + t307;
t24 = (qJD(2) * t292 + t364) * t234 + (-qJD(2) * t334 + t232 * t335 + t235 * t291) * t231;
t23 = (-qJD(2) * t294 - t365) * t234 + (qJD(2) * t333 - t235 * t335 + t287) * t231;
t18 = -t119 * t128 + t120 * t132 + t124 * t241 - t190 * t72 + t191 * t76 + t341 * t68;
t17 = -t119 * t127 + t120 * t131 + t123 * t241 - t190 * t73 + t191 * t77 + t341 * t69;
t16 = t119 * t122 + t120 * t130 + t126 * t241 + t190 * t66 + t191 * t74 + t341 * t70;
t15 = t119 * t121 + t120 * t129 + t125 * t241 + t190 * t67 + t191 * t75 + t341 * t71;
t14 = t117 * t128 + t118 * t132 + t124 * t240 - t192 * t72 + t193 * t76 + t339 * t68;
t13 = t117 * t127 + t118 * t131 + t123 * t240 - t192 * t73 + t193 * t77 + t339 * t69;
t12 = -t117 * t122 + t118 * t130 + t126 * t240 + t192 * t66 + t193 * t74 + t339 * t70;
t11 = -t117 * t121 + t118 * t129 + t125 * t240 + t192 * t67 + t193 * t75 + t339 * t71;
t10 = t365 * t235 + t364 * t232 + (t293 + (-t196 - t333) * t232) * qJD(1) + t307;
t9 = t237 * t320 + (qJD(1) * t385 - t365 * t232 + t364 * t235) * t231;
t8 = qJD(1) * t280 - t17 * t235 + t18 * t232;
t7 = qJD(1) * t281 - t15 * t235 + t16 * t232;
t6 = qJD(1) * t278 - t13 * t235 + t14 * t232;
t5 = qJD(1) * t279 - t11 * t235 + t12 * t232;
t4 = (qJD(2) * t280 - t34) * t234 + (-qJD(1) * t36 + qJD(2) * t87 + t17 * t232 + t18 * t235) * t231;
t3 = (qJD(2) * t281 - t33) * t234 + (-qJD(1) * t35 + qJD(2) * t86 + t15 * t232 + t16 * t235) * t231;
t2 = (qJD(2) * t278 - t32) * t234 + (-qJD(1) * t38 + qJD(2) * t89 + t13 * t232 + t14 * t235) * t231;
t1 = (qJD(2) * t279 - t31) * t234 + (-qJD(1) * t37 + qJD(2) * t88 + t11 * t232 + t12 * t235) * t231;
t25 = [-t143 * t342 - t171 * t302 - t175 * t298 + (t41 * t85 + t42 * t84) * t382 + (t103 * t58 + t104 * t57) * t383 + (t107 * t159 + t108 * t158) * t384 + (-t272 + t277) * t322 + (Icges(3,1) * t231 + t273 + t359) * t320 + t391; m(4) * (t103 * t90 + t104 * t91 + t154 * t57 + t155 * t58) + m(5) * (t105 * t41 + t106 * t42 + t55 * t84 + t56 * t85) + (t381 / 0.2e1 + t380 / 0.2e1) * t269 * qJD(2) + ((qJD(1) * t173 - t232 * t244) * t378 + t177 * t401 + (t350 / 0.2e1 - t345 / 0.2e1) * qJD(2) - t253 + m(3) * (-t108 * t208 - t158 * t201) + (-t159 * t373 + t347 / 0.2e1 + t344 / 0.2e1 + t250) * qJD(1)) * t235 + ((-qJD(1) * t172 - t235 * t244) * t234 / 0.2e1 + t176 * t401 + (-t348 / 0.2e1 + t343 / 0.2e1) * qJD(2) + t252 + m(3) * (-t107 * t208 - t159 * t201) + (t158 * t373 + t349 / 0.2e1 + t346 / 0.2e1 + t251) * qJD(1)) * t232; -t235 * t8 + t232 * t5 + (t10 * t45 + t105 * t56 + t106 * t55) * t382 + t232 * t6 - t235 * t7 + (t154 * t91 + t155 * t90 + t29 * t65) * t383 + ((t232 * t180 + t181 * t235) * ((qJD(1) * t180 - t235 * t248 + t328) * t235 + (-t232 * t248 + (-t181 + t387) * qJD(1)) * t232) + (t380 + t381) * t208 * t201) * t384 - t235 * ((t235 * t141 + (t98 + t388) * qJD(1)) * t235 + (t97 * qJD(1) + (-t173 * t320 - t177 * t322 + t326) * t232 + (-t140 + (t346 + t349) * qJD(2) - t257 * qJD(1)) * t235) * t232) + t232 * ((t232 * t140 + (t99 + t389) * qJD(1)) * t232 + (t100 * qJD(1) + (t172 * t320 + t176 * t322) * t235 + (-t141 + (-t344 - t347) * qJD(2) + (t169 - t258) * qJD(1)) * t232) * t235) + (t98 * t232 - t97 * t235 + t35 + t36) * t325 + (t100 * t232 - t235 * t99 + t37 + t38) * t323; m(4) * (t101 * t58 + t102 * t57 + t103 * t44 + t104 * t43) + m(5) * (t23 * t85 + t24 * t84 + t41 * t64 + t42 * t63) + ((t232 * t251 + t235 * t250) * qJD(2) - t368) * t234 + (t252 * t235 + t253 * t232 + (-t232 * t250 + t235 * t251) * qJD(1)) * t231 - t363; m(4) * (t101 * t90 + t102 * t91 + t154 * t43 + t155 * t44 + t29 * t92 + t30 * t65) + m(5) * (t10 * t46 + t105 * t23 + t106 * t24 + t45 * t9 + t55 * t63 + t56 * t64) + (-t3 / 0.2e1 - t4 / 0.2e1 + t312 * t320) * t235 + (t1 / 0.2e1 + t2 / 0.2e1 + t313 * t320) * t232 + ((-t232 * t312 + t235 * t313) * qJD(1) + (t7 + t8) * t232 / 0.2e1 + (t5 + t6) * t235 / 0.2e1 + (t366 * t232 - t367 * t235) * qJD(2) / 0.2e1) * t231 + (t396 * qJD(1) + t394 * t232 + t395 * t235) * t378 + (t393 * t232 + t392 * t235) * qJD(1) / 0.2e1; (t23 * t64 + t24 * t63 + t46 * t9) * t382 + (t101 * t44 + t102 * t43 + t30 * t92) * t383 + (t368 * t234 + (t255 * t232 - t235 * t386) * qJD(2) + t363) * t234 + ((-t394 * t234 + t1 + t2) * t235 + (t395 * t234 + t3 + t4) * t232 + (t396 * t231 + t399 * t234) * qJD(2) + (t232 * t386 + t255 * t235) * qJD(1)) * t231; m(5) * (-t117 * t84 + t119 * t85 + t190 * t41 + t192 * t42); m(5) * (t45 * t297 + t105 * t119 - t106 * t117 + t190 * t56 + t192 * t55 + (t10 * t231 + t320 * t45) * t230); m(5) * (t46 * t297 - t117 * t63 + t119 * t64 + t190 * t23 + t192 * t24 + (t231 * t9 + t320 * t46) * t230); (-t117 * t192 + t119 * t190 + t397 * t342) * t382;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t25(1), t25(2), t25(4), t25(7); t25(2), t25(3), t25(5), t25(8); t25(4), t25(5), t25(6), t25(9); t25(7), t25(8), t25(9), t25(10);];
Mq = res;
