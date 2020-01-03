% Calculate time derivative of joint inertia matrix for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP6_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_inertiaDJ_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP6_inertiaDJ_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP6_inertiaDJ_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:10
% EndTime: 2019-12-31 17:18:26
% DurationCPUTime: 9.86s
% Computational Cost: add. (9586->658), mult. (25213->935), div. (0->0), fcn. (24176->6), ass. (0->329)
t232 = cos(qJ(3));
t230 = sin(qJ(2));
t233 = cos(qJ(2));
t229 = sin(qJ(3));
t367 = Icges(5,4) * t229;
t275 = Icges(5,1) * t232 - t367;
t171 = -Icges(5,5) * t233 + t230 * t275;
t369 = Icges(4,4) * t229;
t276 = Icges(4,1) * t232 - t369;
t172 = -Icges(4,5) * t233 + t230 * t276;
t417 = t171 + t172;
t421 = t417 * t232;
t268 = Icges(5,5) * t232 - Icges(5,6) * t229;
t163 = -Icges(5,3) * t233 + t230 * t268;
t269 = Icges(4,5) * t232 - Icges(4,6) * t229;
t164 = -Icges(4,3) * t233 + t230 * t269;
t420 = t163 + t164;
t366 = Icges(5,4) * t232;
t271 = -Icges(5,2) * t229 + t366;
t167 = -Icges(5,6) * t233 + t230 * t271;
t368 = Icges(4,4) * t232;
t272 = -Icges(4,2) * t229 + t368;
t168 = -Icges(4,6) * t233 + t230 * t272;
t418 = -t168 - t167;
t228 = -qJ(4) - pkin(6);
t419 = rSges(5,3) - t228;
t416 = -qJD(1) * t230 / 0.2e1;
t415 = t420 * t233 + (-t418 * t229 - t421) * t230;
t323 = qJD(3) * t230;
t142 = (-Icges(5,2) * t232 - t367) * t323 + (Icges(5,6) * t230 + t233 * t271) * qJD(2);
t143 = (-Icges(4,2) * t232 - t369) * t323 + (Icges(4,6) * t230 + t233 * t272) * qJD(2);
t414 = -t143 - t142;
t325 = qJD(2) * t233;
t412 = t418 * t325;
t231 = sin(qJ(1));
t234 = cos(qJ(1));
t343 = t233 * t234;
t346 = t231 * t232;
t188 = -t229 * t343 + t346;
t347 = t231 * t229;
t189 = t232 * t343 + t347;
t348 = t230 * t234;
t121 = Icges(5,5) * t189 + Icges(5,6) * t188 + Icges(5,3) * t348;
t125 = Icges(5,4) * t189 + Icges(5,2) * t188 + Icges(5,6) * t348;
t129 = Icges(5,1) * t189 + Icges(5,4) * t188 + Icges(5,5) * t348;
t266 = -t125 * t229 + t129 * t232;
t62 = -t121 * t233 + t230 * t266;
t123 = Icges(4,5) * t189 + Icges(4,6) * t188 + Icges(4,3) * t348;
t127 = Icges(4,4) * t189 + Icges(4,2) * t188 + Icges(4,6) * t348;
t131 = Icges(4,1) * t189 + Icges(4,4) * t188 + Icges(4,5) * t348;
t264 = -t127 * t229 + t131 * t232;
t64 = -t123 * t233 + t230 * t264;
t378 = t62 + t64;
t344 = t232 * t234;
t345 = t231 * t233;
t186 = -t229 * t345 - t344;
t351 = t229 * t234;
t187 = t232 * t345 - t351;
t350 = t230 * t231;
t120 = Icges(5,5) * t187 + Icges(5,6) * t186 + Icges(5,3) * t350;
t124 = Icges(5,4) * t187 + Icges(5,2) * t186 + Icges(5,6) * t350;
t128 = Icges(5,1) * t187 + Icges(5,4) * t186 + Icges(5,5) * t350;
t267 = -t124 * t229 + t128 * t232;
t61 = -t120 * t233 + t230 * t267;
t122 = Icges(4,5) * t187 + Icges(4,6) * t186 + Icges(4,3) * t350;
t126 = Icges(4,4) * t187 + Icges(4,2) * t186 + Icges(4,6) * t350;
t130 = Icges(4,1) * t187 + Icges(4,4) * t186 + Icges(4,5) * t350;
t265 = -t126 * t229 + t130 * t232;
t63 = -t122 * t233 + t230 * t265;
t379 = t61 + t63;
t411 = t231 * t379 + t234 * t378;
t329 = qJD(1) * t233;
t293 = -qJD(3) + t329;
t294 = qJD(3) * t233 - qJD(1);
t326 = qJD(2) * t231;
t306 = t230 * t326;
t117 = -t294 * t346 + (-t234 * t293 + t306) * t229;
t257 = t294 * t229;
t327 = qJD(2) * t230;
t118 = t293 * t344 + (-t232 * t327 - t257) * t231;
t301 = t231 * t325;
t328 = qJD(1) * t234;
t241 = t230 * t328 + t301;
t69 = Icges(5,5) * t118 + Icges(5,6) * t117 + Icges(5,3) * t241;
t73 = Icges(5,4) * t118 + Icges(5,2) * t117 + Icges(5,6) * t241;
t77 = Icges(5,1) * t118 + Icges(5,4) * t117 + Icges(5,5) * t241;
t19 = (qJD(2) * t267 - t69) * t233 + (qJD(2) * t120 - t229 * t73 + t232 * t77 + (-t124 * t232 - t128 * t229) * qJD(3)) * t230;
t71 = Icges(4,5) * t118 + Icges(4,6) * t117 + Icges(4,3) * t241;
t75 = Icges(4,4) * t118 + Icges(4,2) * t117 + Icges(4,6) * t241;
t79 = Icges(4,1) * t118 + Icges(4,4) * t117 + Icges(4,5) * t241;
t21 = (qJD(2) * t265 - t71) * t233 + (qJD(2) * t122 - t229 * t75 + t232 * t79 + (-t126 * t232 - t130 * t229) * qJD(3)) * t230;
t410 = -t19 - t21;
t256 = t294 * t234;
t324 = qJD(2) * t234;
t304 = t230 * t324;
t398 = t231 * t293 + t304;
t115 = t398 * t229 - t232 * t256;
t116 = -t229 * t256 - t398 * t232;
t300 = t233 * t324;
t330 = qJD(1) * t231;
t308 = t230 * t330;
t240 = t300 - t308;
t68 = Icges(5,5) * t116 + Icges(5,6) * t115 + Icges(5,3) * t240;
t72 = Icges(5,4) * t116 + Icges(5,2) * t115 + Icges(5,6) * t240;
t76 = Icges(5,1) * t116 + Icges(5,4) * t115 + Icges(5,5) * t240;
t20 = (qJD(2) * t266 - t68) * t233 + (qJD(2) * t121 - t229 * t72 + t232 * t76 + (-t125 * t232 - t129 * t229) * qJD(3)) * t230;
t70 = Icges(4,5) * t116 + Icges(4,6) * t115 + Icges(4,3) * t240;
t74 = Icges(4,4) * t116 + Icges(4,2) * t115 + Icges(4,6) * t240;
t78 = Icges(4,1) * t116 + Icges(4,4) * t115 + Icges(4,5) * t240;
t22 = (qJD(2) * t264 - t70) * t233 + (qJD(2) * t123 - t229 * t74 + t232 * t78 + (-t127 * t232 - t131 * t229) * qJD(3)) * t230;
t409 = t20 + t22;
t51 = t122 * t350 + t126 * t186 + t130 * t187;
t52 = t123 * t350 + t127 * t186 + t131 * t187;
t281 = t231 * t51 + t234 * t52;
t49 = t120 * t350 + t124 * t186 + t128 * t187;
t50 = t121 * t350 + t125 * t186 + t129 * t187;
t282 = t231 * t49 + t234 * t50;
t84 = t163 * t350 + t167 * t186 + t171 * t187;
t85 = t164 * t350 + t168 * t186 + t172 * t187;
t408 = (-t84 - t85) * t233 + (t281 + t282) * t230;
t55 = t122 * t348 + t188 * t126 + t189 * t130;
t56 = t123 * t348 + t188 * t127 + t189 * t131;
t279 = t231 * t55 + t234 * t56;
t53 = t120 * t348 + t188 * t124 + t189 * t128;
t54 = t121 * t348 + t188 * t125 + t189 * t129;
t280 = t231 * t53 + t234 * t54;
t86 = t163 * t348 + t188 * t167 + t189 * t171;
t87 = t164 * t348 + t188 * t168 + t189 * t172;
t407 = (-t86 - t87) * t233 + (t279 + t280) * t230;
t371 = Icges(3,4) * t230;
t278 = Icges(3,1) * t233 - t371;
t174 = Icges(3,5) * t231 + t234 * t278;
t354 = t174 * t233;
t370 = Icges(3,4) * t233;
t274 = -Icges(3,2) * t230 + t370;
t170 = Icges(3,6) * t231 + t234 * t274;
t359 = t170 * t230;
t258 = -t354 + t359;
t406 = t231 * t258;
t382 = pkin(6) + t228;
t405 = t233 * t382;
t173 = -Icges(3,5) * t234 + t231 * t278;
t356 = t173 * t233;
t169 = -Icges(3,6) * t234 + t231 * t274;
t361 = t169 * t230;
t259 = -t356 + t361;
t404 = t234 * t259;
t219 = pkin(3) * t232 + pkin(2);
t403 = t189 * rSges(5,1) + t188 * rSges(5,2) + pkin(3) * t347 + t219 * t343 + t419 * t348;
t138 = (-Icges(5,5) * t229 - Icges(5,6) * t232) * t323 + (Icges(5,3) * t230 + t233 * t268) * qJD(2);
t139 = (-Icges(4,5) * t229 - Icges(4,6) * t232) * t323 + (Icges(4,3) * t230 + t233 * t269) * qJD(2);
t146 = (-Icges(5,1) * t229 - t366) * t323 + (Icges(5,5) * t230 + t233 * t275) * qJD(2);
t147 = (-Icges(4,1) * t229 - t368) * t323 + (Icges(4,5) * t230 + t233 * t276) * qJD(2);
t322 = qJD(3) * t232;
t402 = t420 * t327 + t325 * t421 + (-t138 - t139) * t233 + ((t146 + t147) * t232 + t418 * t322) * t230;
t320 = pkin(3) * t351;
t401 = -t187 * rSges(5,1) - t186 * rSges(5,2) + t320;
t317 = pkin(3) * t322;
t321 = qJD(4) * t230;
t400 = t116 * rSges(5,1) + t115 * rSges(5,2) + rSges(5,3) * t300 + qJD(1) * t320 + t228 * t308 + t231 * t317 + t234 * t321;
t399 = -rSges(3,2) * t348 + t231 * rSges(3,3);
t270 = Icges(3,5) * t233 - Icges(3,6) * t230;
t165 = -Icges(3,3) * t234 + t231 * t270;
t397 = t233 * t378 - t407;
t218 = pkin(2) * t343;
t191 = pkin(6) * t348 + t218;
t340 = -t191 + t403;
t383 = pkin(2) - t219;
t239 = -t230 * t382 - t233 * t383;
t341 = rSges(5,3) * t350 + t231 * t239 - t401;
t396 = -t231 * t341 - t234 * t340;
t395 = 2 * m(3);
t394 = 2 * m(4);
t393 = 2 * m(5);
t226 = t231 ^ 2;
t227 = t234 ^ 2;
t391 = -t233 / 0.2e1;
t389 = -rSges(4,3) - pkin(6);
t203 = rSges(3,1) * t230 + rSges(3,2) * t233;
t388 = m(3) * t203;
t387 = pkin(2) * t233;
t386 = pkin(3) * t229;
t223 = t231 * pkin(5);
t380 = t402 + (t412 + (-qJD(3) * t417 + t414) * t230) * t229;
t210 = pkin(6) * t300;
t318 = qJD(3) * t386;
t377 = -t210 + (pkin(6) * t330 + t324 * t383) * t230 + ((-qJD(2) * t228 - t318) * t234 + t383 * t330) * t233 - rSges(5,3) * t308 + t400;
t209 = pkin(2) * t306;
t285 = t118 * rSges(5,1) + t117 * rSges(5,2);
t353 = t219 * t230;
t376 = t209 + (qJD(1) * t239 - t317) * t234 + (t321 - pkin(3) * t257 + (-t353 - t405) * qJD(2)) * t231 + rSges(5,3) * t241 + t285;
t375 = t415 * t327;
t374 = rSges(3,3) * t234;
t373 = rSges(5,3) * t230;
t287 = -rSges(4,1) * t187 - rSges(4,2) * t186;
t135 = rSges(4,3) * t350 - t287;
t362 = t135 * t234;
t360 = t169 * t233;
t358 = t170 * t233;
t357 = t173 * t230;
t355 = t174 * t230;
t283 = rSges(5,1) * t232 - rSges(5,2) * t229;
t303 = t229 * t323;
t342 = -pkin(3) * t303 - qJD(4) * t233 + (-rSges(5,1) * t229 - rSges(5,2) * t232) * t323 + (t233 * t283 + t239 + t373) * qJD(2);
t137 = t189 * rSges(4,1) + t188 * rSges(4,2) + rSges(4,3) * t348;
t339 = -t137 - t191;
t286 = rSges(4,1) * t232 - rSges(4,2) * t229;
t151 = (-rSges(4,1) * t229 - rSges(4,2) * t232) * t323 + (rSges(4,3) * t230 + t233 * t286) * qJD(2);
t291 = pkin(6) * t230 + t387;
t197 = t291 * qJD(2);
t338 = -t151 - t197;
t337 = -rSges(5,3) * t233 + t405 + (t283 - t383) * t230;
t176 = -rSges(4,3) * t233 + t230 * t286;
t204 = pkin(2) * t230 - pkin(6) * t233;
t336 = -t176 - t204;
t190 = t291 * t231;
t335 = t231 * t190 + t234 * t191;
t334 = rSges(3,2) * t308 + rSges(3,3) * t328;
t333 = t234 * pkin(1) + t223;
t332 = t226 + t227;
t166 = Icges(3,3) * t231 + t234 * t270;
t331 = qJD(1) * t166;
t35 = t50 * t231 - t234 * t49;
t36 = t52 * t231 - t234 * t51;
t316 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t54 * t231 - t234 * t53;
t38 = t56 * t231 - t234 * t55;
t315 = t37 / 0.2e1 + t38 / 0.2e1;
t313 = t116 * rSges(4,1) + t115 * rSges(4,2) + rSges(4,3) * t300;
t312 = -t197 - t342;
t242 = -t231 * t329 - t304;
t311 = t231 * (pkin(6) * t241 + qJD(1) * t218 - t209) + t234 * (pkin(2) * t242 - pkin(6) * t308 + t210) + t190 * t328;
t310 = -t204 - t337;
t309 = t176 * t330;
t299 = -t219 * t233 - pkin(1);
t298 = t231 * t337;
t297 = t234 * t337;
t296 = t341 * t234;
t155 = t336 * t234;
t295 = qJD(1) * t337;
t106 = t310 * t234;
t290 = t231 * t295;
t289 = rSges(3,1) * t233 - rSges(3,2) * t230;
t288 = t118 * rSges(4,1) + t117 * rSges(4,2);
t273 = Icges(3,2) * t233 + t371;
t263 = -t137 * t231 + t362;
t262 = -t231 * t135 - t137 * t234;
t178 = rSges(3,1) * t343 + t399;
t255 = -t233 * t379 + t408;
t254 = -pkin(1) - t289;
t31 = t115 * t167 + t116 * t171 + t138 * t348 + t188 * t142 + t189 * t146 + t163 * t240;
t32 = t115 * t168 + t116 * t172 + t139 * t348 + t188 * t143 + t189 * t147 + t164 * t240;
t252 = t22 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1 + t20 / 0.2e1;
t33 = t117 * t167 + t118 * t171 + t138 * t350 + t186 * t142 + t187 * t146 + t163 * t241;
t34 = t117 * t168 + t118 * t172 + t139 * t350 + t186 * t143 + t187 * t147 + t164 * t241;
t251 = t33 / 0.2e1 + t34 / 0.2e1 + t19 / 0.2e1 + t21 / 0.2e1;
t250 = t61 / 0.2e1 + t63 / 0.2e1 + t84 / 0.2e1 + t85 / 0.2e1;
t249 = t62 / 0.2e1 + t64 / 0.2e1 + t86 / 0.2e1 + t87 / 0.2e1;
t248 = qJD(2) * t203;
t247 = t230 * t389 - pkin(1) - t387;
t245 = qJD(2) * t273;
t244 = qJD(2) * (-Icges(3,5) * t230 - Icges(3,6) * t233);
t243 = -t230 * t419 + t299;
t238 = t247 * t231;
t237 = -t231 * t340 + t296;
t224 = t234 * pkin(5);
t221 = pkin(5) * t328;
t196 = t289 * qJD(2);
t192 = t204 * t330;
t177 = t231 * t289 - t374;
t157 = t178 + t333;
t156 = t231 * t254 + t224 + t374;
t154 = t336 * t231;
t141 = t231 * t244 + t331;
t140 = -qJD(1) * t165 + t234 * t244;
t108 = t203 * t326 + ((-rSges(3,3) - pkin(5)) * t231 + t254 * t234) * qJD(1);
t107 = rSges(3,1) * t242 - rSges(3,2) * t300 - pkin(1) * t330 + t221 + t334;
t105 = t310 * t231;
t104 = t333 - t339;
t103 = t224 + t238 + t287;
t102 = -t233 * t137 - t176 * t348;
t101 = t135 * t233 + t176 * t350;
t100 = t231 * t166 - t258 * t234;
t99 = t231 * t165 - t404;
t98 = -t234 * t166 - t406;
t97 = -t234 * t165 - t231 * t259;
t94 = t333 + t403;
t93 = t231 * t243 + t224 + t401;
t90 = t263 * t230;
t89 = qJD(1) * t155 + t231 * t338;
t88 = t234 * t338 + t192 + t309;
t83 = rSges(4,3) * t241 + t288;
t81 = -rSges(4,3) * t308 + t313;
t65 = -t262 + t335;
t60 = -t230 * t297 - t233 * t340;
t59 = t230 * t298 + t233 * t341;
t58 = t209 + t389 * t301 + (t234 * t247 - t223) * qJD(1) - t288;
t57 = -pkin(2) * t304 + qJD(1) * t238 + t210 + t221 + t313;
t48 = qJD(1) * t106 + t231 * t312;
t47 = t234 * t312 + t192 + t290;
t46 = t237 * t230;
t45 = (qJD(1) * t243 + t317) * t234 + (-pkin(5) * qJD(1) - t321 + t294 * t386 + (-t233 * t419 + t353) * qJD(2)) * t231 - t285;
t44 = t221 + (-t233 * t318 + (-t228 * t233 - t353) * qJD(2)) * t234 + (t299 - t373) * t330 + t400;
t43 = t335 - t396;
t42 = (t176 * t326 + t83) * t233 + (-qJD(2) * t135 + t231 * t151 + t176 * t328) * t230;
t41 = (-t176 * t324 - t81) * t233 + (qJD(2) * t137 - t151 * t234 + t309) * t230;
t30 = t263 * t325 + (qJD(1) * t262 - t231 * t81 + t234 * t83) * t230;
t29 = t231 * t83 + t234 * t81 + (t231 * t339 + t362) * qJD(1) + t311;
t24 = (qJD(2) * t298 + t376) * t233 + (-qJD(2) * t341 + t231 * t342 + t234 * t295) * t230;
t23 = (-qJD(2) * t297 - t377) * t233 + (qJD(2) * t340 - t234 * t342 + t290) * t230;
t18 = t117 * t127 + t118 * t131 + t123 * t241 + t186 * t74 + t187 * t78 + t350 * t70;
t17 = t117 * t126 + t118 * t130 + t122 * t241 + t186 * t75 + t187 * t79 + t350 * t71;
t16 = t117 * t125 + t118 * t129 + t121 * t241 + t186 * t72 + t187 * t76 + t350 * t68;
t15 = t117 * t124 + t118 * t128 + t120 * t241 + t186 * t73 + t187 * t77 + t350 * t69;
t14 = t115 * t127 + t116 * t131 + t123 * t240 + t188 * t74 + t189 * t78 + t348 * t70;
t13 = t115 * t126 + t116 * t130 + t122 * t240 + t188 * t75 + t189 * t79 + t348 * t71;
t12 = t115 * t125 + t116 * t129 + t121 * t240 + t188 * t72 + t189 * t76 + t348 * t68;
t11 = t115 * t124 + t116 * t128 + t120 * t240 + t188 * t73 + t189 * t77 + t348 * t69;
t10 = t377 * t234 + t376 * t231 + (t296 + (-t191 - t340) * t231) * qJD(1) + t311;
t9 = t237 * t325 + (t396 * qJD(1) - t377 * t231 + t376 * t234) * t230;
t8 = qJD(1) * t281 - t17 * t234 + t18 * t231;
t7 = qJD(1) * t282 - t15 * t234 + t16 * t231;
t6 = qJD(1) * t279 - t13 * t234 + t14 * t231;
t5 = qJD(1) * t280 - t11 * t234 + t12 * t231;
t4 = (qJD(2) * t281 - t34) * t233 + (-qJD(1) * t36 + qJD(2) * t85 + t17 * t231 + t18 * t234) * t230;
t3 = (qJD(2) * t282 - t33) * t233 + (-qJD(1) * t35 + qJD(2) * t84 + t15 * t231 + t16 * t234) * t230;
t2 = (qJD(2) * t279 - t32) * t233 + (-qJD(1) * t38 + qJD(2) * t87 + t13 * t231 + t14 * t234) * t230;
t1 = (qJD(2) * t280 - t31) * t233 + (-qJD(1) * t37 + qJD(2) * t86 + t11 * t231 + t12 * t234) * t230;
t25 = [(t44 * t94 + t45 * t93) * t393 + (t103 * t58 + t104 * t57) * t394 + (t107 * t157 + t108 * t156) * t395 + (t278 - t273) * t327 + (Icges(3,1) * t230 + t274 + t370) * t325 - t417 * t303 + t402 + (t230 * t414 + t412) * t229; m(4) * (t103 * t88 + t104 * t89 + t154 * t57 + t155 * t58) + m(5) * (t105 * t44 + t106 * t45 + t47 * t93 + t48 * t94) + (t226 / 0.2e1 + t227 / 0.2e1) * t270 * qJD(2) + ((qJD(1) * t170 - t231 * t245) * t391 + t174 * t416 + (t361 / 0.2e1 - t356 / 0.2e1) * qJD(2) - t251 + m(3) * (-t108 * t203 - t156 * t196) + (t358 / 0.2e1 + t355 / 0.2e1 - t157 * t388 + t249) * qJD(1)) * t234 + ((-t169 * qJD(1) - t234 * t245) * t233 / 0.2e1 + t173 * t416 + (-t359 / 0.2e1 + t354 / 0.2e1) * qJD(2) + t252 + m(3) * (-t107 * t203 - t157 * t196) + (t156 * t388 + t360 / 0.2e1 + t357 / 0.2e1 + t250) * qJD(1)) * t231; (t10 * t43 + t105 * t48 + t106 * t47) * t393 + t231 * t5 - t234 * t7 + (t154 * t89 + t155 * t88 + t29 * t65) * t394 + t231 * t6 - t234 * t8 + ((t231 * t177 + t178 * t234) * ((qJD(1) * t177 - t234 * t248 + t334) * t234 + (-t231 * t248 + (-t178 + t399) * qJD(1)) * t231) + t332 * t203 * t196) * t395 + t231 * ((t231 * t140 + (t99 + t406) * qJD(1)) * t231 + (t100 * qJD(1) + (t169 * t325 + t173 * t327) * t234 + (-t141 + (-t355 - t358) * qJD(2) + (t166 - t259) * qJD(1)) * t231) * t234) - t234 * ((t234 * t141 + (t98 + t404) * qJD(1)) * t234 + (t97 * qJD(1) + (-t170 * t325 - t174 * t327 + t331) * t231 + (-t140 + (t357 + t360) * qJD(2) - t258 * qJD(1)) * t234) * t231) + (t98 * t231 - t97 * t234 + t35 + t36) * t330 + (t100 * t231 - t234 * t99 + t37 + t38) * t328; m(4) * (t101 * t58 + t102 * t57 + t103 * t42 + t104 * t41) + m(5) * (t23 * t94 + t24 * t93 + t44 * t60 + t45 * t59) + ((t231 * t250 + t234 * t249) * qJD(2) - t380) * t233 + (t252 * t234 + t251 * t231 + (-t231 * t249 + t234 * t250) * qJD(1)) * t230 - t375; m(4) * (t101 * t88 + t102 * t89 + t154 * t41 + t155 * t42 + t29 * t90 + t30 * t65) + m(5) * (t10 * t46 + t105 * t23 + t106 * t24 + t43 * t9 + t47 * t59 + t48 * t60) + (-t3 / 0.2e1 - t4 / 0.2e1 + t315 * t325) * t234 + (t1 / 0.2e1 + t2 / 0.2e1 + t316 * t325) * t231 + ((-t231 * t315 + t234 * t316) * qJD(1) + (t7 + t8) * t231 / 0.2e1 + (t5 + t6) * t234 / 0.2e1 + (t378 * t231 - t379 * t234) * qJD(2) / 0.2e1) * t230 + (qJD(1) * t411 + t409 * t231 + t410 * t234) * t391 + (t408 * t231 + t407 * t234) * qJD(1) / 0.2e1; (t23 * t60 + t24 * t59 + t46 * t9) * t393 + (t101 * t42 + t102 * t41 + t30 * t90) * t394 + (t380 * t233 + (t255 * t231 - t397 * t234) * qJD(2) + t375) * t233 + ((-t409 * t233 + t1 + t2) * t234 + (t410 * t233 + t3 + t4) * t231 + (t411 * t230 + t233 * t415) * qJD(2) + (t397 * t231 + t255 * t234) * qJD(1)) * t230; m(5) * ((t231 * t94 + t234 * t93) * t325 + (t231 * t44 + t234 * t45 + (-t231 * t93 + t234 * t94) * qJD(1)) * t230); m(5) * ((-t10 + (t105 * t231 + t106 * t234) * qJD(2)) * t233 + (qJD(2) * t43 + t231 * t48 + t234 * t47 + (t105 * t234 - t106 * t231) * qJD(1)) * t230); m(5) * ((-t9 + (t231 * t60 + t234 * t59) * qJD(2)) * t233 + (qJD(2) * t46 + t23 * t231 + t234 * t24 + (-t231 * t59 + t234 * t60) * qJD(1)) * t230); (-0.1e1 + t332) * t230 * t325 * t393;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t25(1), t25(2), t25(4), t25(7); t25(2), t25(3), t25(5), t25(8); t25(4), t25(5), t25(6), t25(9); t25(7), t25(8), t25(9), t25(10);];
Mq = res;
