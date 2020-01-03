% Calculate vector of inverse dynamics joint torques for
% S4RPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR3_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:05
% EndTime: 2019-12-31 16:49:17
% DurationCPUTime: 10.36s
% Computational Cost: add. (9969->546), mult. (9178->729), div. (0->0), fcn. (7207->8), ass. (0->303)
t237 = cos(qJ(1));
t230 = t237 * pkin(1);
t232 = qJ(1) + pkin(7);
t224 = sin(t232);
t225 = cos(t232);
t156 = rSges(3,1) * t224 + rSges(3,2) * t225;
t235 = sin(qJ(1));
t398 = t235 * pkin(1);
t147 = -t156 - t398;
t240 = qJD(1) ^ 2;
t319 = t240 * t230;
t233 = qJ(3) + qJ(4);
t227 = cos(t233);
t221 = Icges(5,4) * t227;
t226 = sin(t233);
t164 = Icges(5,1) * t226 + t221;
t274 = -Icges(5,2) * t226 + t221;
t422 = t164 + t274;
t360 = t224 * t226;
t177 = Icges(5,4) * t360;
t359 = t224 * t227;
t376 = Icges(5,5) * t225;
t110 = Icges(5,1) * t359 - t177 - t376;
t372 = Icges(5,6) * t225;
t108 = Icges(5,4) * t359 - Icges(5,2) * t360 - t372;
t369 = t108 * t226;
t273 = -t110 * t227 + t369;
t259 = t273 * t224;
t161 = Icges(5,5) * t227 - Icges(5,6) * t226;
t107 = Icges(5,3) * t224 + t161 * t225;
t378 = Icges(5,4) * t226;
t165 = Icges(5,1) * t227 - t378;
t111 = Icges(5,5) * t224 + t165 * t225;
t354 = t225 * t227;
t393 = t224 * t107 + t111 * t354;
t421 = -t259 - t393;
t180 = rSges(5,2) * t360;
t112 = rSges(5,1) * t359 - t225 * rSges(5,3) - t180;
t213 = t224 * rSges(5,3);
t355 = t225 * t226;
t113 = rSges(5,1) * t354 - rSges(5,2) * t355 + t213;
t389 = rSges(5,2) * t227;
t166 = rSges(5,1) * t226 + t389;
t130 = t166 * t224;
t131 = t166 * t225;
t231 = qJD(3) + qJD(4);
t154 = t224 * t231;
t155 = t225 * t231;
t219 = t225 * pkin(5);
t158 = pkin(2) * t224 - t219;
t238 = -pkin(6) - pkin(5);
t204 = t225 * t238;
t236 = cos(qJ(3));
t229 = t236 * pkin(3);
t395 = t229 + pkin(2);
t337 = -t224 * t395 - t204;
t100 = t158 + t337;
t159 = t225 * pkin(2) + t224 * pkin(5);
t294 = -t224 * t238 + t225 * t395;
t101 = t294 - t159;
t37 = t112 * t154 + t113 * t155 + qJD(2) + (-t100 * t224 + t101 * t225) * qJD(3);
t234 = sin(qJ(3));
t324 = qJD(3) * t234;
t313 = t225 * t324;
t292 = pkin(3) * t313;
t255 = -t155 * t166 - t292;
t419 = -t112 - t398;
t290 = t100 + t419;
t265 = -t158 + t290;
t41 = qJD(1) * t265 + t255;
t222 = t227 * rSges(5,1);
t418 = -rSges(5,2) * t226 + t222;
t317 = pkin(3) * t324;
t183 = t224 * t317;
t304 = t159 + t230;
t346 = t101 + t113;
t42 = -t154 * t166 - t183 + (t304 + t346) * qJD(1);
t322 = qJD(1) * qJD(3);
t150 = qJDD(3) * t224 + t225 * t322;
t321 = qJD(1) * qJD(4);
t104 = qJDD(4) * t224 + t225 * t321 + t150;
t200 = t224 * t322;
t105 = t224 * t321 + t200 + (-qJDD(3) - qJDD(4)) * t225;
t151 = -qJDD(3) * t225 + t200;
t318 = t231 * t389;
t329 = qJD(1) * t224;
t328 = qJD(1) * t225;
t338 = rSges(5,3) * t328 + qJD(1) * t180;
t351 = t226 * t231;
t67 = -t225 * t318 + (-t225 * t351 - t227 * t329) * rSges(5,1) + t338;
t263 = t166 * t231;
t68 = -t224 * t263 + (t225 * t418 + t213) * qJD(1);
t205 = pkin(5) * t328;
t84 = -t292 - t205 + (-t224 * t229 - t204) * qJD(1);
t85 = -t183 + (t225 * t229 + (-pkin(5) - t238) * t224) * qJD(1);
t8 = -t100 * t150 - t101 * t151 + t104 * t112 - t105 * t113 + t154 * t68 + t155 * t67 + qJDD(2) + (t224 * t85 + t225 * t84) * qJD(3);
t420 = -(qJD(1) * t130 - t155 * t418) * t41 - t37 * (-t154 * t130 - t131 * t155) - t42 * (-qJD(1) * t131 - t154 * t418) + t8 * (t224 * t112 + t225 * t113);
t214 = t224 * rSges(4,3);
t352 = t225 * t236;
t353 = t225 * t234;
t122 = rSges(4,1) * t352 - rSges(4,2) * t353 + t214;
t87 = t122 + t304;
t157 = t225 * rSges(3,1) - rSges(3,2) * t224;
t148 = t157 + t230;
t228 = Icges(4,4) * t236;
t275 = -Icges(4,2) * t234 + t228;
t189 = Icges(4,1) * t234 + t228;
t186 = Icges(4,5) * t236 - Icges(4,6) * t234;
t185 = Icges(4,5) * t234 + Icges(4,6) * t236;
t256 = qJD(3) * t185;
t379 = Icges(4,4) * t234;
t190 = Icges(4,1) * t236 - t379;
t120 = Icges(4,5) * t224 + t190 * t225;
t118 = Icges(4,6) * t224 + t225 * t275;
t366 = t118 * t234;
t270 = -t120 * t236 + t366;
t371 = Icges(4,3) * t225;
t416 = -t225 * t256 + (-t186 * t224 + t270 + t371) * qJD(1);
t358 = t224 * t234;
t197 = Icges(4,4) * t358;
t357 = t224 * t236;
t377 = Icges(4,5) * t225;
t119 = Icges(4,1) * t357 - t197 - t377;
t373 = Icges(4,6) * t225;
t117 = Icges(4,4) * t357 - Icges(4,2) * t358 - t373;
t367 = t117 * t234;
t271 = -t119 * t236 + t367;
t116 = Icges(4,3) * t224 + t186 * t225;
t331 = qJD(1) * t116;
t415 = qJD(1) * t271 - t224 * t256 + t331;
t160 = Icges(5,5) * t226 + Icges(5,6) * t227;
t260 = t160 * t231;
t109 = Icges(5,6) * t224 + t225 * t274;
t368 = t109 * t226;
t370 = Icges(5,3) * t225;
t414 = -t225 * t260 + (-t111 * t227 - t161 * t224 + t368 + t370) * qJD(1);
t332 = qJD(1) * t107;
t413 = qJD(1) * t273 - t224 * t260 + t332;
t115 = Icges(4,5) * t357 - Icges(4,6) * t358 - t371;
t47 = -t225 * t115 - t224 * t271;
t162 = Icges(5,2) * t227 + t378;
t268 = t162 * t226 - t164 * t227;
t412 = qJD(1) * t268 + t161 * t231;
t187 = Icges(4,2) * t236 + t379;
t266 = t187 * t234 - t189 * t236;
t411 = qJD(1) * t266 + t186 * qJD(3);
t341 = -Icges(4,2) * t357 + t119 - t197;
t343 = t189 * t224 + t117;
t410 = -t234 * t341 - t236 * t343;
t409 = qJD(1) * t422 + t154 * (-t162 * t225 + t111) - t155 * (-Icges(5,2) * t359 + t110 - t177);
t408 = t104 / 0.2e1;
t407 = t105 / 0.2e1;
t406 = t150 / 0.2e1;
t405 = t151 / 0.2e1;
t404 = -t154 / 0.2e1;
t403 = t154 / 0.2e1;
t402 = -t155 / 0.2e1;
t401 = t155 / 0.2e1;
t400 = t224 / 0.2e1;
t399 = -t225 / 0.2e1;
t397 = -qJD(1) / 0.2e1;
t396 = qJD(1) / 0.2e1;
t106 = Icges(5,5) * t359 - Icges(5,6) * t360 - t370;
t394 = -t224 * t106 - t110 * t354;
t392 = rSges(4,1) * t236;
t191 = rSges(4,1) * t234 + rSges(4,2) * t236;
t145 = t191 * t225;
t326 = qJD(3) * t224;
t60 = qJD(1) * t87 - t191 * t326;
t388 = t145 * t60;
t333 = rSges(4,2) * t358 + t225 * rSges(4,3);
t121 = rSges(4,1) * t357 - t333;
t305 = -t121 - t398;
t289 = -t158 + t305;
t325 = qJD(3) * t225;
t314 = t191 * t325;
t59 = qJD(1) * t289 - t314;
t387 = t224 * t59;
t386 = t225 * t42;
t385 = t225 * t59;
t384 = qJDD(1) / 0.2e1;
t383 = -t224 * t115 - t119 * t352;
t382 = t224 * t116 + t120 * t352;
t365 = t160 * t224;
t364 = t160 * t225;
t363 = t162 * t231;
t362 = t185 * t224;
t361 = t185 * t225;
t350 = t236 * qJD(3) ^ 2;
t57 = -t224 * t268 - t364;
t348 = t57 * qJD(1);
t80 = -t224 * t266 - t361;
t347 = t80 * qJD(1);
t342 = -t189 * t225 - t118;
t340 = -t187 * t225 + t120;
t327 = qJD(1) * t234;
t336 = rSges(4,2) * t224 * t327 + rSges(4,3) * t328;
t335 = -t187 + t190;
t334 = t189 + t275;
t330 = qJD(1) * t186;
t323 = qJD(3) * t236;
t320 = t112 * t328 + t224 * t68 + t225 * t67;
t312 = -pkin(2) - t392;
t311 = t329 / 0.2e1;
t310 = t328 / 0.2e1;
t309 = -t326 / 0.2e1;
t308 = t326 / 0.2e1;
t307 = -t325 / 0.2e1;
t306 = t325 / 0.2e1;
t254 = -pkin(3) * t234 - t166;
t262 = t164 * t231;
t303 = qJD(1) * t111 - t108 * t231 - t224 * t262;
t302 = -t109 * t231 - t225 * t262 + (-t165 * t224 + t376) * qJD(1);
t301 = qJD(1) * t109 + t110 * t231 - t224 * t363;
t300 = t111 * t231 - t225 * t363 + (-t224 * t274 + t372) * qJD(1);
t97 = t120 * t357;
t299 = t225 * t116 - t97;
t298 = -t106 + t368;
t297 = -t115 + t366;
t296 = t422 * t231;
t295 = t165 * t231 - t363;
t291 = qJDD(1) * t230 - t240 * t398;
t137 = t418 * t231;
t288 = -pkin(3) * t323 - t137;
t88 = t111 * t359;
t287 = t109 * t360 - t88;
t194 = rSges(2,1) * t237 - rSges(2,2) * t235;
t192 = rSges(2,1) * t235 + rSges(2,2) * t237;
t193 = -rSges(4,2) * t234 + t392;
t70 = t118 * t236 + t120 * t234;
t257 = qJD(3) * t187;
t76 = -t225 * t257 + (-t224 * t275 + t373) * qJD(1);
t258 = qJD(3) * t189;
t78 = -t225 * t258 + (-t190 * t224 + t377) * qJD(1);
t243 = -qJD(3) * t70 - t234 * t76 + t236 * t78 + t331;
t69 = t117 * t236 + t119 * t234;
t77 = qJD(1) * t118 - t224 * t257;
t79 = qJD(1) * t120 - t224 * t258;
t244 = qJD(1) * t115 - qJD(3) * t69 - t234 * t77 + t236 * t79;
t284 = -(t224 * t415 + t244 * t225) * t225 + (t224 * t416 + t243 * t225) * t224;
t283 = -(t244 * t224 - t225 * t415) * t225 + (t243 * t224 - t225 * t416) * t224;
t282 = -t224 * t42 - t225 * t41;
t48 = -t118 * t358 - t299;
t281 = t224 * t48 - t225 * t47;
t49 = -t117 * t353 - t383;
t50 = -t118 * t353 + t382;
t280 = t224 * t50 - t225 * t49;
t279 = -t224 * t60 - t385;
t82 = -rSges(4,2) * t225 * t323 + (-t236 * t329 - t313) * rSges(4,1) + t336;
t144 = t191 * t224;
t83 = -qJD(3) * t144 + (t193 * t225 + t214) * qJD(1);
t278 = t224 * t83 + t225 * t82;
t51 = t108 * t227 + t110 * t226;
t269 = t121 * t224 + t122 * t225;
t267 = t187 * t236 + t189 * t234;
t264 = qJD(1) * (-pkin(2) * t329 + t205) + qJDD(1) * t159 + t291;
t253 = qJD(1) * t161 - t154 * t364 + t155 * t365;
t252 = -t234 * t340 + t236 * t342;
t251 = (-t234 * t334 + t236 * t335) * qJD(1);
t246 = -t226 * t300 + t227 * t302 + t332;
t10 = t224 * t414 + t246 * t225;
t247 = qJD(1) * t106 - t226 * t301 + t227 * t303;
t11 = t247 * t224 - t225 * t413;
t12 = t246 * t224 - t225 * t414;
t43 = -t106 * t225 - t259;
t44 = -t107 * t225 - t287;
t22 = t154 * t44 - t155 * t43 + t348;
t45 = -t108 * t355 - t394;
t46 = -t109 * t355 + t393;
t58 = -t225 * t268 + t365;
t53 = t58 * qJD(1);
t23 = t154 * t46 - t155 * t45 + t53;
t248 = (-t164 * t225 - t109) * t154 - (-t164 * t224 - t108) * t155 + (-t162 + t165) * qJD(1);
t241 = -t226 * t409 + t248 * t227;
t28 = t226 * t303 + t227 * t301;
t29 = t226 * t302 + t227 * t300;
t245 = qJD(1) * t160 - t226 * t296 + t227 * t295;
t30 = t224 * t412 + t245 * t225;
t31 = t245 * t224 - t225 * t412;
t52 = t109 * t227 + t111 * t226;
t9 = t224 * t413 + t247 * t225;
t250 = (qJD(1) * t30 + qJDD(1) * t58 + t10 * t154 + t104 * t46 + t105 * t45 - t155 * t9) * t400 + (t248 * t226 + t227 * t409) * t397 + t22 * t311 + t23 * t310 + (qJD(1) * t31 + qJDD(1) * t57 + t104 * t44 + t105 * t43 - t11 * t155 + t12 * t154) * t399 + (t224 * t46 - t225 * t45) * t408 + (t224 * t44 - t225 * t43) * t407 + (t10 * t224 - t225 * t9 + (t224 * t45 + t225 * t46) * qJD(1)) * t403 + (t224 * t52 - t225 * t51) * t384 + (-t11 * t225 + t12 * t224 + (t224 * t43 + t225 * t44) * qJD(1)) * t402 + (t224 * t29 - t225 * t28 + (t224 * t51 + t225 * t52) * qJD(1)) * t396 + (t224 * t253 + t225 * t241) * t404 + (t224 * t241 - t225 * t253) * t401;
t169 = t275 * qJD(3);
t170 = t190 * qJD(3);
t242 = qJD(1) * t185 - qJD(3) * t267 - t169 * t234 + t170 * t236;
t172 = t193 * qJD(3);
t153 = qJD(1) * t158;
t149 = t159 * qJD(1);
t81 = -t225 * t266 + t362;
t71 = t81 * qJD(1);
t56 = qJD(3) * t269 + qJD(2);
t40 = t242 * t224 - t225 * t411;
t39 = t224 * t411 + t242 * t225;
t36 = qJD(1) * t82 + qJDD(1) * t122 - t150 * t191 - t172 * t326 + t264;
t35 = -t319 - t172 * t325 + t151 * t191 + (-t149 - t83) * qJD(1) + t289 * qJDD(1);
t34 = -qJD(3) * t270 + t234 * t78 + t236 * t76;
t33 = -t271 * qJD(3) + t234 * t79 + t236 * t77;
t32 = qJD(3) * t278 + t121 * t150 - t122 * t151 + qJDD(2);
t27 = qJD(3) * t280 + t71;
t26 = qJD(3) * t281 + t347;
t21 = -t104 * t166 - t137 * t154 + t346 * qJDD(1) + (t67 + t84) * qJD(1) + (-t150 * t234 - t224 * t350) * pkin(3) + t264;
t20 = -t319 + t105 * t166 - t137 * t155 + (t151 * t234 - t225 * t350) * pkin(3) + (-t149 - t68 - t85) * qJD(1) + t265 * qJDD(1);
t1 = [(t53 + (t44 + (t107 + t369) * t225 + t287 + t394) * t155 + (-t225 * t298 - t421 + t43) * t154) * t401 + (t71 + ((t48 - t97 + (t116 + t367) * t225 + t383) * t225 + t382 * t224) * qJD(3)) * t306 - m(2) * (-g(1) * t192 + g(2) * t194) + (t52 + t58) * t408 + (t51 + t57) * t407 + (t70 + t81) * t406 + (t69 + t80) * t405 + (-t348 + (t46 + t421) * t155 + (t224 * t298 + t45 - t88) * t154 + ((t107 + t273) * t154 + t298 * t155) * t225 + t22) * t404 + (t30 + t29) * t403 + (-t347 + ((t225 * t297 - t382 + t50) * t225 + (t224 * t297 + t299 + t49) * t224) * qJD(3) + t26) * t309 + (t34 + t39) * t308 + (-qJD(3) * t266 + t169 * t236 + t170 * t234 + t226 * t295 + t227 * t296) * qJD(1) + ((-t156 * t240 - g(2) + t291) * t148 + (-t319 + (-0.2e1 * t157 - t230 + t148) * t240 - g(1)) * t147) * m(3) + (t28 + t31 + t23) * t402 + (t33 + t40 + t27) * t307 + (t41 * (t183 + (rSges(5,1) * t351 + t318) * t224) + t42 * t338 + (-t263 - t317) * t386 + ((-t235 * t42 - t237 * t41) * pkin(1) + (t41 * (-t418 - t395) - t42 * t238) * t225 + (t41 * (t238 - rSges(5,3)) + t42 * (-t222 - t395)) * t224) * qJD(1) - (qJD(1) * t290 - t153 + t255 - t41) * t42 + (t21 - g(2)) * (t113 + t230 + t294) + (t20 - g(1)) * (t337 + t419)) * m(5) + (t60 * (t205 + t336) + (t191 * t387 - t388) * qJD(3) + ((-t235 * t60 - t237 * t59) * pkin(1) + (-pkin(2) - t193) * t385 + (t59 * (-rSges(4,3) - pkin(5)) + t60 * t312) * t224) * qJD(1) - (qJD(1) * t305 - t153 - t314 - t59) * t60 + (t36 - g(2)) * t87 + (t35 - g(1)) * (t312 * t224 + t219 + t333 - t398)) * m(4) + (m(3) * (t147 ^ 2 + t157 * t148) + t267 + m(2) * (t192 ^ 2 + t194 ^ 2) + t162 * t227 + t164 * t226 + Icges(2,3) + Icges(3,3)) * qJDD(1); m(3) * qJDD(2) + m(4) * t32 + m(5) * t8 + (-m(3) - m(4) - m(5)) * g(3); ((t49 * t224 + t50 * t225) * qJD(1) + t284) * t308 + ((t47 * t224 + t48 * t225) * qJD(1) + t283) * t307 + (t224 * t34 - t225 * t33 + (t69 * t224 + t225 * t70) * qJD(1)) * t396 + ((t234 * t335 + t236 * t334) * qJD(1) + ((t224 * t340 - t225 * t341) * t236 + (t224 * t342 + t225 * t343) * t234) * qJD(3)) * t397 + ((-t325 * t362 - t330) * t225 + (t251 + (t252 * t224 + (t361 - t410) * t225) * qJD(3)) * t224) * t306 + t250 + ((-t326 * t361 + t330) * t224 + (t251 + (-t410 * t225 + (t362 + t252) * t224) * qJD(3)) * t225) * t309 + (qJD(1) * t40 + qJD(3) * t283 + qJDD(1) * t80 + t150 * t48 + t151 * t47) * t399 + (qJD(1) * t39 + qJD(3) * t284 + qJDD(1) * t81 + t150 * t50 + t151 * t49) * t400 + t26 * t311 + t27 * t310 + (t224 * t70 - t225 * t69) * t384 + t280 * t406 + t281 * t405 + (t37 * t320 + (t20 * t254 + t41 * t288 + t8 * t101 + t37 * t84 + (-t37 * t100 + t254 * t42) * qJD(1)) * t225 + (t21 * t254 + t42 * t288 - t8 * t100 + t37 * t85 + (t41 * t166 - t346 * t37) * qJD(1)) * t224 - g(3) * (t418 + t229) - (g(1) * t225 + g(2) * t224) * t254 - (-t327 * t386 + (t282 * t236 + t37 * (-t224 ^ 2 - t225 ^ 2) * t234) * qJD(3)) * pkin(3) + t420) * m(5) + (t32 * t269 + t56 * ((t121 * t225 - t122 * t224) * qJD(1) + t278) + t279 * t172 + (-t36 * t224 - t35 * t225 + (-t225 * t60 + t387) * qJD(1)) * t191 - (t144 * t59 - t388) * qJD(1) - (t56 * (-t144 * t224 - t145 * t225) + t279 * t193) * qJD(3) + g(1) * t145 + g(2) * t144 - g(3) * t193) * m(4); t250 + (t37 * (-t113 * t329 + t320) + t282 * t137 + (-t20 * t225 - t21 * t224 + (t224 * t41 - t386) * qJD(1)) * t166 + g(1) * t131 + g(2) * t130 - g(3) * t418 + t420) * m(5);];
tau = t1;
