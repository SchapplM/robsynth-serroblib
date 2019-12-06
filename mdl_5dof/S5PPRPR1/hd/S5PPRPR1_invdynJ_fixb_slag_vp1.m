% Calculate vector of inverse dynamics joint torques for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR1_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR1_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:53
% EndTime: 2019-12-05 15:01:18
% DurationCPUTime: 20.88s
% Computational Cost: add. (12513->579), mult. (13753->878), div. (0->0), fcn. (13018->8), ass. (0->301)
t250 = sin(pkin(7));
t248 = pkin(8) + qJ(3);
t242 = sin(t248);
t244 = cos(t248);
t311 = -Icges(4,5) * t242 - Icges(4,6) * t244;
t187 = t311 * t250;
t252 = cos(pkin(7));
t188 = t311 * t252;
t245 = t250 ^ 2;
t246 = t252 ^ 2;
t374 = t245 + t246;
t459 = qJD(3) * t374;
t312 = Icges(4,5) * t244 - Icges(4,6) * t242;
t458 = (-Icges(4,3) * t252 + t250 * t312) * t252;
t455 = 0.2e1 * qJD(3);
t454 = 2 * qJDD(3);
t247 = pkin(9) + qJ(5);
t241 = sin(t247);
t414 = rSges(6,2) * t241;
t243 = cos(t247);
t416 = rSges(6,1) * t243;
t333 = -t414 + t416;
t143 = t242 * rSges(6,3) + t244 * t333;
t249 = sin(pkin(9));
t415 = rSges(5,2) * t249;
t251 = cos(pkin(9));
t417 = rSges(5,1) * t251;
t334 = -t415 + t417;
t453 = t242 * rSges(5,3) + t244 * t334;
t390 = t251 * t252;
t393 = t249 * t250;
t197 = -t244 * t393 - t390;
t391 = t250 * t251;
t392 = t249 * t252;
t198 = t244 * t391 - t392;
t398 = t242 * t250;
t80 = Icges(5,5) * t198 + Icges(5,6) * t197 + Icges(5,3) * t398;
t199 = -t244 * t392 + t391;
t200 = t244 * t390 + t393;
t397 = t242 * t252;
t81 = Icges(5,5) * t200 + Icges(5,6) * t199 + Icges(5,3) * t397;
t82 = Icges(5,4) * t198 + Icges(5,2) * t197 + Icges(5,6) * t398;
t83 = Icges(5,4) * t200 + Icges(5,2) * t199 + Icges(5,6) * t397;
t84 = Icges(5,1) * t198 + Icges(5,4) * t197 + Icges(5,5) * t398;
t85 = Icges(5,1) * t200 + Icges(5,4) * t199 + Icges(5,5) * t397;
t452 = (-(t249 * t82 - t251 * t84) * t242 - t80 * t244 + t188) * t252 + ((t249 * t83 - t251 * t85) * t242 + t81 * t244 + t187) * t250;
t438 = rSges(5,3) * t244 - t242 * t334;
t451 = t438 * t459;
t363 = qJD(3) * qJD(4);
t448 = qJDD(4) * t242 + t244 * t363;
t372 = qJD(3) * t242;
t371 = qJD(3) * t244;
t440 = t244 * pkin(3) + t242 * qJ(4);
t439 = g(1) * t252 + g(2) * t250;
t310 = Icges(6,5) * t243 - Icges(6,6) * t241;
t135 = -Icges(6,3) * t244 + t242 * t310;
t403 = Icges(6,4) * t243;
t313 = -Icges(6,2) * t241 + t403;
t137 = -Icges(6,6) * t244 + t242 * t313;
t404 = Icges(6,4) * t241;
t316 = Icges(6,1) * t243 - t404;
t139 = -Icges(6,5) * t244 + t242 * t316;
t365 = qJD(5) * t242;
t370 = qJD(3) * t250;
t203 = t252 * t365 + t370;
t369 = qJD(3) * t252;
t204 = t250 * t365 - t369;
t389 = t252 * t241;
t184 = t243 * t250 - t244 * t389;
t395 = t244 * t252;
t185 = t241 * t250 + t243 * t395;
t405 = Icges(6,4) * t185;
t73 = Icges(6,2) * t184 + Icges(6,6) * t397 + t405;
t166 = Icges(6,4) * t184;
t75 = Icges(6,1) * t185 + Icges(6,5) * t397 + t166;
t327 = -t241 * t73 + t243 * t75;
t396 = t244 * t250;
t182 = -t241 * t396 - t243 * t252;
t183 = t243 * t396 - t389;
t406 = Icges(6,4) * t183;
t72 = Icges(6,2) * t182 + Icges(6,6) * t398 + t406;
t165 = Icges(6,4) * t182;
t74 = Icges(6,1) * t183 + Icges(6,5) * t398 + t165;
t328 = -t241 * t72 + t243 * t74;
t437 = -(-t135 * t252 - t327) * t203 - (-t135 * t250 - t328) * t204;
t253 = -pkin(6) - qJ(4);
t388 = qJ(4) + t253;
t238 = pkin(4) * t251 + pkin(3);
t420 = pkin(3) - t238;
t436 = t242 * t420 - t244 * t388;
t168 = (-Icges(6,2) * t243 - t404) * t242;
t364 = qJD(5) * t244;
t259 = t203 * (-Icges(6,2) * t185 + t166 + t75) + t204 * (-Icges(6,2) * t183 + t165 + t74) - t364 * (t139 + t168);
t435 = -m(5) - m(6);
t362 = qJD(3) * qJD(5);
t284 = qJDD(5) * t242 + t244 * t362;
t360 = qJDD(3) * t250;
t144 = t252 * t284 + t360;
t434 = t144 / 0.2e1;
t359 = qJDD(3) * t252;
t145 = t250 * t284 - t359;
t433 = t145 / 0.2e1;
t432 = -t374 * t372 / 0.2e1;
t202 = -qJDD(5) * t244 + t242 * t362;
t431 = t202 / 0.2e1;
t430 = -t203 / 0.2e1;
t429 = t203 / 0.2e1;
t428 = -t204 / 0.2e1;
t427 = t204 / 0.2e1;
t426 = -t244 / 0.2e1;
t421 = t242 * pkin(3);
t71 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t397;
t25 = t182 * t73 + t183 * t75 + t398 * t71;
t412 = t25 * t252;
t70 = Icges(6,5) * t183 + Icges(6,6) * t182 + Icges(6,3) * t398;
t26 = t184 * t72 + t185 * t74 + t397 * t70;
t411 = t26 * t250;
t43 = t135 * t398 + t137 * t182 + t139 * t183;
t410 = t43 * t242;
t44 = t135 * t397 + t137 * t184 + t139 * t185;
t409 = t44 * t242;
t399 = t135 * t244;
t394 = t244 * t253;
t368 = qJD(4) * t244;
t170 = qJD(3) * t440 - t368;
t279 = -t242 * t388 - t244 * t420;
t387 = -t279 * qJD(3) - t170;
t367 = qJD(4) * t250;
t219 = t242 * t367;
t207 = -t244 * qJ(4) + t421;
t293 = qJD(3) * t207;
t133 = -t250 * t293 + t219;
t366 = qJD(4) * t252;
t221 = t242 * t366;
t134 = -t252 * t293 + t221;
t386 = t250 * t133 + t252 * t134;
t385 = -t207 + t436;
t213 = t244 * t238;
t384 = t242 * t253 - t213;
t382 = -qJD(3) * t453 - t170;
t381 = -t207 + t438;
t194 = t440 * t250;
t196 = t440 * t252;
t380 = t250 * t194 + t252 * t196;
t354 = t242 * t414;
t379 = rSges(6,3) * t396 + t250 * t354;
t378 = rSges(6,3) * t395 + t252 * t354;
t355 = t242 * t415;
t377 = rSges(5,3) * t396 + t250 * t355;
t376 = rSges(5,3) * t395 + t252 * t355;
t240 = qJD(2) * t250;
t375 = t221 + t240;
t373 = qJD(2) * t252;
t361 = qJDD(2) * t252;
t357 = t242 * t417;
t356 = t242 * t416;
t171 = (-rSges(6,1) * t241 - rSges(6,2) * t243) * t242;
t69 = qJD(3) * t143 + qJD(5) * t171;
t353 = -t69 + t387;
t352 = -m(3) - m(4) + t435;
t142 = -rSges(6,3) * t244 + t242 * t333;
t351 = -t142 + t385;
t224 = qJ(4) * t396;
t225 = qJ(4) * t395;
t350 = (-pkin(3) * t398 + t224) * t370 + (-pkin(3) * t397 + t225) * t369 + qJD(4) * t242;
t239 = qJDD(2) * t250;
t349 = t252 * t448 + t239;
t348 = t242 * t370;
t347 = t242 * t369;
t346 = t244 * t370;
t345 = t244 * t369;
t341 = t369 / 0.2e1;
t340 = -t364 / 0.2e1;
t339 = t364 / 0.2e1;
t337 = qJD(3) * t385;
t336 = qJD(3) * t381;
t335 = t219 - t373;
t210 = rSges(4,1) * t244 - rSges(4,2) * t242;
t208 = rSges(4,1) * t242 + rSges(4,2) * t244;
t277 = qJD(3) * t387 + qJDD(3) * t385;
t100 = -qJD(5) * t183 + t241 * t348;
t101 = qJD(5) * t182 - t243 * t348;
t59 = rSges(6,1) * t101 + rSges(6,2) * t100 + rSges(6,3) * t346;
t76 = rSges(6,1) * t183 + rSges(6,2) * t182 + rSges(6,3) * t398;
t13 = t142 * t145 - t202 * t76 + t204 * t69 + t252 * t277 + t364 * t59 + t349;
t319 = t250 * t448 - t361;
t102 = -qJD(5) * t185 + t241 * t347;
t103 = qJD(5) * t184 - t243 * t347;
t60 = rSges(6,1) * t103 + rSges(6,2) * t102 + rSges(6,3) * t345;
t77 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t397;
t14 = -t142 * t144 + t202 * t77 - t203 * t69 + t250 * t277 - t364 * t60 + t319;
t331 = t13 * t250 - t14 * t252;
t330 = t71 * t203 + t70 * t204;
t24 = t182 * t72 + t183 * t74 + t398 * t70;
t329 = t24 * t250 + t412;
t27 = t184 * t73 + t185 * t75 + t397 * t71;
t324 = t27 * t252 + t411;
t28 = t242 * t328 - t244 * t70;
t29 = t242 * t327 - t244 * t71;
t323 = t28 * t250 + t29 * t252;
t322 = -t250 * t77 + t252 * t76;
t78 = -pkin(4) * t392 + t250 * t279;
t79 = pkin(4) * t393 + t252 * t279;
t321 = t250 * t78 + t252 * t79;
t88 = rSges(5,1) * t198 + rSges(5,2) * t197 + rSges(5,3) * t398;
t89 = rSges(5,1) * t200 + rSges(5,2) * t199 + rSges(5,3) * t397;
t320 = t250 * t88 + t252 * t89;
t309 = -t137 * t241 + t139 * t243;
t306 = -(-t208 * t369 + t240) * t252 - (-t208 * t370 - t373) * t250;
t305 = t374 * t210;
t304 = t208 * t459;
t201 = t210 * qJD(3);
t303 = -qJD(3) * t201 - qJDD(3) * t208;
t300 = t194 * t370 + t196 * t369 + qJD(1) - t368;
t53 = Icges(6,5) * t101 + Icges(6,6) * t100 + Icges(6,3) * t346;
t299 = t242 * t53 + t371 * t70;
t54 = Icges(6,5) * t103 + Icges(6,6) * t102 + Icges(6,3) * t345;
t298 = t242 * t54 + t371 * t71;
t136 = Icges(6,3) * t242 + t244 * t310;
t167 = (-Icges(6,5) * t241 - Icges(6,6) * t243) * t242;
t66 = qJD(3) * t136 + qJD(5) * t167;
t297 = t135 * t371 + t242 * t66;
t296 = t136 - t309;
t23 = qJD(3) * t321 + t203 * t76 - t204 * t77 + t300;
t295 = t23 * t322;
t169 = (-Icges(6,1) * t241 - t403) * t242;
t281 = qJD(3) * t311;
t280 = -t242 * t238 - t394 + t421;
t278 = t167 * t364 - t203 * (Icges(6,5) * t184 - Icges(6,6) * t185) - t204 * (Icges(6,5) * t182 - Icges(6,6) * t183);
t276 = qJD(3) * t382 + qJDD(3) * t381;
t275 = -qJDD(4) * t244 + t133 * t370 + t134 * t369 + t194 * t360 + t196 * t359 + t242 * t363 + qJDD(1);
t274 = qJD(3) * t436;
t273 = Icges(5,5) * t244 + (-Icges(5,1) * t251 + Icges(5,4) * t249) * t242;
t140 = Icges(6,5) * t242 + t244 * t316;
t271 = Icges(5,6) * t244 + (-Icges(5,4) * t251 + Icges(5,2) * t249) * t242;
t138 = Icges(6,6) * t242 + t244 * t313;
t267 = t242 * t278;
t265 = qJD(3) * t273;
t264 = qJD(3) * t271;
t260 = (Icges(6,1) * t184 - t405 - t73) * t203 + (Icges(6,1) * t182 - t406 - t72) * t204 - (-t137 + t169) * t364;
t254 = (-t296 * t364 - t437) * t242;
t222 = t244 * t366;
t220 = t244 * t367;
t195 = t208 * t252;
t193 = t208 * t250;
t175 = t252 * t281;
t174 = t250 * t281;
t150 = Icges(4,3) * t250 + t252 * t312;
t130 = t273 * t252;
t129 = t273 * t250;
t128 = t271 * t252;
t127 = t271 * t250;
t120 = -t252 * t356 + t378;
t119 = -t250 * t356 + t379;
t115 = t139 * t252;
t114 = t139 * t250;
t113 = t137 * t252;
t112 = t137 * t250;
t109 = t252 * t265;
t108 = t250 * t265;
t107 = t252 * t264;
t106 = t250 * t264;
t99 = t250 * t303 - t361;
t98 = t252 * t303 + t239;
t97 = rSges(6,1) * t184 - rSges(6,2) * t185;
t96 = rSges(6,1) * t182 - rSges(6,2) * t183;
t87 = t252 * t274;
t86 = t250 * t274;
t68 = qJD(3) * t140 + qJD(5) * t169;
t67 = qJD(3) * t138 + qJD(5) * t168;
t65 = qJD(3) * t305 + qJD(1);
t64 = t250 * t336 + t335;
t63 = t252 * t336 + t375;
t58 = Icges(6,1) * t103 + Icges(6,4) * t102 + Icges(6,5) * t345;
t57 = Icges(6,1) * t101 + Icges(6,4) * t100 + Icges(6,5) * t346;
t56 = Icges(6,4) * t103 + Icges(6,2) * t102 + Icges(6,6) * t345;
t55 = Icges(6,4) * t101 + Icges(6,2) * t100 + Icges(6,6) * t346;
t48 = -qJD(3) * t304 + qJDD(3) * t305 + qJDD(1);
t47 = t242 * t309 - t399;
t46 = t250 * t276 + t319;
t45 = t252 * t276 + t349;
t38 = qJD(3) * t320 + t300;
t37 = -t142 * t203 + t250 * t337 - t364 * t77 + t335;
t36 = t142 * t204 + t252 * t337 + t364 * t76 + t375;
t22 = qJD(3) * t451 + t320 * qJDD(3) + t275;
t17 = (qJD(3) * t309 - t66) * t244 + (qJD(3) * t135 - t241 * t67 + t243 * t68 + (-t137 * t243 - t139 * t241) * qJD(5)) * t242;
t16 = t102 * t137 + t103 * t139 + t184 * t67 + t185 * t68 + t252 * t297;
t15 = t100 * t137 + t101 * t139 + t182 * t67 + t183 * t68 + t250 * t297;
t12 = t203 * t29 + t204 * t28 - t364 * t47;
t11 = t102 * t73 + t103 * t75 + t184 * t56 + t185 * t58 + t252 * t298;
t10 = t102 * t72 + t103 * t74 + t184 * t55 + t185 * t57 + t252 * t299;
t9 = t100 * t73 + t101 * t75 + t182 * t56 + t183 * t58 + t250 * t298;
t8 = t100 * t72 + t101 * t74 + t182 * t55 + t183 * t57 + t250 * t299;
t7 = t203 * t27 + t204 * t26 - t364 * t44;
t6 = t203 * t25 + t204 * t24 - t364 * t43;
t5 = (qJD(3) * t327 - t54) * t244 + (qJD(3) * t71 - t241 * t56 + t243 * t58 + (-t241 * t75 - t243 * t73) * qJD(5)) * t242;
t4 = (qJD(3) * t328 - t53) * t244 + (qJD(3) * t70 - t241 * t55 + t243 * t57 + (-t241 * t74 - t243 * t72) * qJD(5)) * t242;
t3 = t144 * t76 - t145 * t77 + t203 * t59 - t204 * t60 + t321 * qJDD(3) + (t250 * t86 + t252 * t87) * qJD(3) + t275;
t2 = t10 * t204 + t203 * t11 + t144 * t27 + t145 * t26 - t16 * t364 + t202 * t44;
t1 = t144 * t25 + t145 * t24 - t15 * t364 + t202 * t43 + t203 * t9 + t204 * t8;
t18 = [(m(2) + m(3)) * qJDD(1) + m(4) * t48 + m(5) * t22 + m(6) * t3 + (-m(2) + t352) * g(3); t352 * (g(1) * t250 - g(2) * t252) + m(4) * (t250 * t98 - t252 * t99) + m(5) * (t250 * t45 - t252 * t46) + m(6) * t331 + m(3) * t374 * qJDD(2); (t250 * t9 - t252 * t8) * t427 + ((-t113 * t182 - t115 * t183) * t203 + (-t112 * t182 - t114 * t183) * t204 + (t410 + (-t138 * t182 - t140 * t183 + t412) * t244) * qJD(5) + (((t24 - t399) * qJD(5) + t330) * t244 + t254) * t250) * t428 + (-t10 * t252 + t11 * t250) * t429 + ((-t113 * t184 - t115 * t185) * t203 + (-t112 * t184 - t114 * t185) * t204 + (t409 + (-t138 * t184 - t140 * t185 + t411) * t244) * qJD(5) + (((t27 - t399) * qJD(5) + t330) * t244 + t254) * t252) * t430 + (t250 * t29 - t252 * t28) * t431 + (-t24 * t252 + t25 * t250) * t433 + (t250 * t27 - t252 * t26) * t434 - t12 * t365 / 0.2e1 + (((t113 * t241 - t115 * t243 + t71) * t203 + (t112 * t241 - t114 * t243 + t70) * t204 + t47 * qJD(5)) * t242 + ((t296 * t244 + (t138 * t241 - t140 * t243 - t135) * t242 + t323) * qJD(5) + t437) * t244) * t339 - (t188 * qJD(3) * t245 + (t128 * t199 + t130 * t200) * t370 + (-t199 * t127 - t200 * t129 - t250 * t187 + t452) * t369) * t370 / 0.2e1 + (t187 * qJD(3) * t246 - (t127 * t197 + t129 * t198) * t369 + (t197 * t128 + t198 * t130 - t252 * t188 + t452) * t370) * t341 + ((-t4 + t7) * t252 + (t5 + t6) * t250) * t340 + (-g(1) * (-t252 * t394 + t378) - g(2) * (-t250 * t394 + t379) - g(3) * (t143 + t213) - (-g(3) * t253 + t439 * (-t238 - t416)) * t242 - t36 * (t143 * t204 + t222) - t37 * (-t143 * t203 + t220) - t23 * (t203 * t119 - t204 * t120 + t350) - ((t36 * t384 + (t252 * t280 - t225) * t23) * t252 + (t37 * t384 + (t250 * t280 - t224) * t23) * t250) * qJD(3) - ((-t36 * t76 + t37 * t77) * t242 + (t36 * (t142 * t250 + t119) + t37 * (-t142 * t252 - t120) + t295) * t244) * qJD(5) + t3 * t380 + t23 * t386 + (t13 * t351 + t36 * t353 + t3 * (t77 + t79) + t23 * (t60 + t87)) * t252 + (t14 * t351 + t37 * t353 + t3 * (t76 + t78) + t23 * (t59 + t86)) * t250) * m(6) + (-g(1) * (t225 + t376) - g(2) * (t224 + t377) - t439 * t242 * (-pkin(3) - t417) - t63 * t222 - t64 * t220 + t22 * t380 + (t22 * t89 + t381 * t45 + t382 * t63) * t252 + (t22 * t88 + t381 * t46 + t382 * t64) * t250 + (-(-t64 * t250 - t63 * t252) * qJD(3) - g(3)) * (t440 + t453) + (-t350 - (t250 * (-t250 * t357 + t377) + t252 * (-t252 * t357 + t376)) * qJD(3) + t386 + t451) * t38) * m(5) + (-(t65 * (-t193 * t250 - t195 * t252) + t306 * t210) * qJD(3) + t48 * t305 - t65 * t304 + (-t250 * t99 - t252 * t98) * t208 + t306 * t201 + g(1) * t195 + g(2) * t193 - g(3) * t210) * m(4) + (t2 + ((t107 * t199 + t109 * t200 + t175 * t250) * t250 + (-t106 * t199 - t108 * t200 - t174 * t250) * t252) * t455 + ((-t199 * t82 - t200 * t84 - t397 * t80) * t252 + (t150 * t250 + t199 * t83 + t200 * t85 + t397 * t81 - t458) * t250) * t454) * t250 / 0.2e1 - (t1 + ((-t106 * t197 - t108 * t198 + t174 * t252) * t252 + (t107 * t197 + t109 * t198 - t175 * t252) * t250) * t455 + ((-t197 * t82 - t198 * t84 - t398 * t80 + t458) * t252 + (-t150 * t252 + t197 * t83 + t198 * t85 + t398 * t81) * t250) * t454) * t252 / 0.2e1; -t435 * g(3) * t244 + 0.2e1 * (t23 * t432 + t3 * t426) * m(6) + 0.2e1 * (t22 * t426 + t38 * t432) * m(5) + (t435 * t439 + m(5) * (qJD(3) * t38 + t250 * t46 + t252 * t45) + m(6) * (qJD(3) * t23 + t13 * t252 + t14 * t250)) * t242; t244 * t7 * t341 + t2 * t397 / 0.2e1 + (t242 * t324 - t244 * t44) * t434 + (-t16 * t244 + (t10 * t250 + t11 * t252) * t242 + (t244 * t324 + t409) * qJD(3)) * t429 + t6 * t346 / 0.2e1 + t1 * t398 / 0.2e1 + (t242 * t329 - t244 * t43) * t433 + (-t15 * t244 + (t250 * t8 + t252 * t9) * t242 + (t244 * t329 + t410) * qJD(3)) * t427 + t12 * t372 / 0.2e1 + (t29 * t144 + t28 * t145 - t17 * t364 + t47 * t202 + t5 * t203 + t4 * t204) * t426 + (t242 * t323 - t244 * t47) * t431 + (-t17 * t244 + (t250 * t4 + t252 * t5) * t242 + (t47 * t242 + t244 * t323) * qJD(3)) * t340 + (t184 * t259 + t185 * t260 - t252 * t267) * t430 + (t182 * t259 + t183 * t260 - t250 * t267) * t428 + (t278 * t244 + (-t241 * t259 + t260 * t243) * t242) * t339 + ((t13 * t76 - t14 * t77 + t36 * t59 - t37 * t60 + (t295 + (t250 * t36 - t252 * t37) * t142) * qJD(3)) * t244 + (t36 * (-qJD(3) * t76 + t250 * t69) + t37 * (qJD(3) * t77 - t252 * t69) + t3 * t322 + t23 * (-t250 * t60 + t252 * t59) + t331 * t142) * t242 - t36 * (t171 * t204 + t364 * t96) - t37 * (-t171 * t203 - t364 * t97) - t23 * (t203 * t96 - t204 * t97) - g(1) * t97 - g(2) * t96 - g(3) * t171) * m(6);];
tau = t18;
