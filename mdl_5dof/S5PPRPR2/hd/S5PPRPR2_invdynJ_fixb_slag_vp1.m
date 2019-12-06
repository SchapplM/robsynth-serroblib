% Calculate vector of inverse dynamics joint torques for
% S5PPRPR2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPRPR2_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:02:59
% EndTime: 2019-12-05 15:03:21
% DurationCPUTime: 17.43s
% Computational Cost: add. (8533->549), mult. (12390->826), div. (0->0), fcn. (11559->6), ass. (0->278)
t455 = Icges(5,4) - Icges(4,5);
t454 = Icges(5,5) - Icges(4,6);
t236 = pkin(8) + qJ(3);
t232 = sin(t236);
t233 = cos(t236);
t449 = t455 * t232 + t454 * t233;
t453 = Icges(5,1) + Icges(4,3);
t452 = t449 * qJD(3);
t451 = t454 * t232 - t455 * t233;
t239 = sin(qJ(5));
t240 = cos(qJ(5));
t323 = rSges(6,1) * t239 + rSges(6,2) * t240;
t432 = t233 * t323;
t237 = sin(pkin(7));
t238 = cos(pkin(7));
t450 = t237 * t238;
t235 = t238 ^ 2;
t448 = t451 * t237 - t453 * t238;
t447 = t453 * t237 + t451 * t238;
t446 = t452 * t237;
t445 = t452 * t238;
t115 = t233 * rSges(6,3) + t232 * t323;
t234 = t237 ^ 2;
t365 = t234 + t235;
t388 = Icges(5,6) * t233;
t296 = Icges(5,3) * t232 - t388;
t124 = -Icges(5,5) * t238 + t237 * t296;
t384 = t232 * t237;
t216 = Icges(5,6) * t384;
t382 = t233 * t237;
t126 = -Icges(5,4) * t238 - Icges(5,2) * t382 + t216;
t362 = qJD(3) * t233;
t363 = qJD(3) * t232;
t394 = Icges(4,4) * t233;
t395 = Icges(4,4) * t232;
t305 = -Icges(4,2) * t232 + t394;
t119 = -Icges(4,6) * t238 + t237 * t305;
t308 = Icges(4,1) * t233 - t395;
t121 = -Icges(4,5) * t238 + t237 * t308;
t429 = t119 * t233 + t121 * t232;
t443 = ((-Icges(5,6) * t232 - t395 + (-Icges(4,2) - Icges(5,3)) * t233) * t363 + (t388 + t394 + (Icges(4,1) + Icges(5,2)) * t232) * t362) * t237 + (-t124 * t233 - t126 * t232 + t429) * qJD(3);
t441 = (-t121 + t126) * t233 + (t119 - t124) * t232;
t440 = 0.2e1 * qJD(3);
t439 = 2 * qJDD(3);
t437 = t365 * t363;
t436 = t449 * t237;
t435 = t449 * t238;
t381 = t233 * t238;
t383 = t232 * t238;
t434 = ((Icges(5,5) * t237 + t238 * t296) * t237 - t124 * t238) * t233 + ((Icges(5,4) * t237 + 0.2e1 * Icges(5,6) * t383 + (-Icges(5,2) + Icges(5,3)) * t381) * t237 - (Icges(5,3) * t382 + t126 + t216) * t238) * t232 + t237 * (-(Icges(4,6) * t237 + t238 * t305) * t233 - (Icges(4,5) * t237 + t238 * t308) * t232) + t238 * t429;
t352 = qJD(3) * qJD(4);
t433 = qJDD(4) * t232 + t233 * t352;
t322 = -rSges(5,2) * t233 + t232 * rSges(5,3);
t423 = t233 * pkin(3) + t232 * qJ(4);
t422 = g(1) * t238 + g(2) * t237;
t298 = Icges(6,5) * t239 + Icges(6,6) * t240;
t254 = -Icges(6,3) * t232 + t233 * t298;
t391 = Icges(6,4) * t239;
t301 = Icges(6,2) * t240 + t391;
t255 = -Icges(6,6) * t232 + t233 * t301;
t390 = Icges(6,4) * t240;
t306 = Icges(6,1) * t239 + t390;
t256 = -Icges(6,5) * t232 + t233 * t306;
t377 = t238 * t240;
t380 = t237 * t239;
t183 = t232 * t377 - t380;
t173 = Icges(6,4) * t183;
t378 = t238 * t239;
t379 = t237 * t240;
t185 = t232 * t379 + t378;
t174 = Icges(6,4) * t185;
t180 = (Icges(6,2) * t239 - t390) * t233;
t184 = t232 * t378 + t379;
t186 = -t232 * t380 + t377;
t353 = qJD(5) * t233;
t361 = qJD(3) * t237;
t190 = t238 * t353 + t361;
t359 = qJD(3) * t238;
t191 = t237 * t353 - t359;
t354 = qJD(5) * t232;
t74 = Icges(6,1) * t184 + Icges(6,5) * t381 + t173;
t75 = -Icges(6,1) * t186 + Icges(6,5) * t382 + t174;
t246 = t190 * (-Icges(6,2) * t184 + t173 + t74) + t191 * (Icges(6,2) * t186 + t174 + t75) + t354 * (-t256 + t180);
t181 = (-Icges(6,1) * t240 + t391) * t233;
t392 = Icges(6,4) * t186;
t393 = Icges(6,4) * t184;
t72 = Icges(6,2) * t183 + Icges(6,6) * t381 + t393;
t73 = Icges(6,2) * t185 + Icges(6,6) * t382 - t392;
t247 = t190 * (-Icges(6,1) * t183 + t393 + t72) + t191 * (-Icges(6,1) * t185 - t392 + t73) + t354 * (-t255 - t181);
t421 = -m(5) - m(6);
t351 = qJD(3) * qJD(5);
t267 = qJDD(5) * t233 - t232 * t351;
t349 = qJDD(3) * t237;
t106 = t238 * t267 + t349;
t420 = t106 / 0.2e1;
t348 = qJDD(3) * t238;
t107 = t237 * t267 - t348;
t419 = t107 / 0.2e1;
t418 = -t437 / 0.2e1;
t189 = qJDD(5) * t232 + t233 * t351;
t417 = t189 / 0.2e1;
t416 = -t190 / 0.2e1;
t415 = t190 / 0.2e1;
t414 = -t191 / 0.2e1;
t413 = t191 / 0.2e1;
t412 = -t233 / 0.2e1;
t356 = qJD(4) * t237;
t209 = t232 * t356;
t194 = pkin(3) * t232 - qJ(4) * t233;
t276 = qJD(3) * t194;
t104 = -t237 * t276 + t209;
t355 = qJD(4) * t238;
t211 = t232 * t355;
t105 = -t238 * t276 + t211;
t403 = t237 * t104 + t238 * t105;
t71 = -Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t382;
t21 = t183 * t73 + t184 * t75 + t381 * t71;
t399 = t21 * t237;
t70 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t381;
t22 = t185 * t72 - t186 * t74 + t382 * t70;
t398 = t22 * t238;
t38 = -t183 * t255 - t184 * t256 - t254 * t381;
t397 = t38 * t233;
t39 = -t185 * t255 + t186 * t256 - t254 * t382;
t396 = t39 * t233;
t385 = t254 * t232;
t357 = qJD(4) * t233;
t141 = qJD(3) * t423 - t357;
t374 = -t322 * qJD(3) - t141;
t176 = t423 * t237;
t178 = t423 * t238;
t373 = t237 * t176 + t238 * t178;
t321 = rSges(5,2) * t232 + rSges(5,3) * t233;
t372 = -t194 + t321;
t371 = -t423 - t322;
t370 = t432 * t237;
t369 = t432 * t238;
t231 = qJD(2) * t237;
t368 = t211 + t231;
t367 = rSges(5,2) * t384 + rSges(5,3) * t382;
t366 = rSges(5,2) * t383 + rSges(5,3) * t381;
t364 = qJD(2) * t238;
t350 = qJDD(2) * t238;
t342 = -m(3) - m(4) + t421;
t215 = qJ(4) * t381;
t325 = -pkin(3) * t383 + t215;
t214 = qJ(4) * t382;
t326 = -pkin(3) * t384 + t214;
t341 = qJD(4) * t232 + t325 * t359 + t326 * t361;
t230 = qJDD(2) * t237;
t340 = t238 * t433 + t230;
t339 = t232 * t361;
t338 = t232 * t359;
t337 = t239 * t362;
t336 = t240 * t362;
t334 = -t361 / 0.2e1;
t332 = -t354 / 0.2e1;
t331 = t354 / 0.2e1;
t330 = -pkin(6) * t232 - t194;
t329 = t365 * t232;
t328 = qJD(3) * t372;
t327 = t209 - t364;
t116 = rSges(6,3) * t232 - t432;
t324 = -t116 + t330;
t199 = rSges(4,1) * t233 - rSges(4,2) * t232;
t196 = rSges(4,1) * t232 + rSges(4,2) * t233;
t241 = qJD(3) ^ 2;
t248 = -qJD(3) * t141 - qJDD(3) * t194 + (-qJDD(3) * t232 - t233 * t241) * pkin(6);
t102 = qJD(5) * t185 + t237 * t337;
t103 = qJD(5) * t186 + t237 * t336;
t60 = rSges(6,1) * t102 + rSges(6,2) * t103 - rSges(6,3) * t339;
t182 = (-rSges(6,1) * t240 + rSges(6,2) * t239) * t233;
t69 = qJD(3) * t115 + qJD(5) * t182;
t77 = -rSges(6,1) * t186 + rSges(6,2) * t185 + rSges(6,3) * t382;
t13 = t107 * t116 - t189 * t77 + t191 * t69 + t238 * t248 - t354 * t60 + t340;
t310 = t237 * t433 - t350;
t100 = qJD(5) * t183 + t238 * t337;
t101 = -qJD(5) * t184 + t238 * t336;
t59 = rSges(6,1) * t100 + rSges(6,2) * t101 - rSges(6,3) * t338;
t76 = rSges(6,1) * t184 + rSges(6,2) * t183 + rSges(6,3) * t381;
t14 = -t106 * t116 + t189 * t76 - t190 * t69 + t237 * t248 + t354 * t59 + t310;
t319 = t13 * t237 - t14 * t238;
t318 = -t70 * t190 - t71 * t191;
t20 = t183 * t72 + t184 * t74 + t381 * t70;
t317 = t20 * t238 + t399;
t23 = t185 * t73 - t186 * t75 + t382 * t71;
t316 = t23 * t237 + t398;
t312 = t239 * t74 + t240 * t72;
t24 = t232 * t70 - t233 * t312;
t311 = t239 * t75 + t240 * t73;
t25 = t232 * t71 - t233 * t311;
t315 = t25 * t237 + t24 * t238;
t309 = qJD(3) * t330;
t36 = t116 * t191 + t238 * t309 - t354 * t77 + t368;
t37 = -t116 * t190 + t237 * t309 + t354 * t76 + t327;
t314 = -t237 * t37 - t238 * t36;
t313 = t237 * t76 - t238 * t77;
t295 = -t239 * t256 - t240 * t255;
t290 = -(-t196 * t359 + t231) * t238 - (-t196 * t361 - t364) * t237;
t289 = t365 * t199;
t139 = rSges(5,1) * t237 + t238 * t322;
t140 = -rSges(5,1) * t238 + t237 * t322;
t288 = t139 * t238 + t140 * t237;
t287 = t365 * qJD(3) * t196;
t192 = pkin(4) * t237 + pkin(6) * t381;
t193 = -pkin(4) * t238 + pkin(6) * t382;
t286 = t192 * t238 + t193 * t237;
t285 = -pkin(6) * t362 - t141 - t69;
t188 = t199 * qJD(3);
t284 = -qJD(3) * t188 - qJDD(3) * t196;
t283 = t176 * t361 + t178 * t359 + qJD(1) - t357;
t53 = Icges(6,5) * t100 + Icges(6,6) * t101 - Icges(6,3) * t338;
t282 = t233 * t53 - t363 * t70;
t54 = Icges(6,5) * t102 + Icges(6,6) * t103 - Icges(6,3) * t339;
t281 = t233 * t54 - t363 * t71;
t108 = Icges(6,3) * t233 + t232 * t298;
t179 = (-Icges(6,5) * t240 + Icges(6,6) * t239) * t233;
t66 = qJD(3) * t108 + qJD(5) * t179;
t280 = t233 * t66 + t254 * t363;
t19 = qJD(3) * t286 + t190 * t77 - t191 * t76 + t283;
t279 = t19 * t313;
t277 = qJD(3) * t321;
t260 = (t108 + t295) * t232;
t259 = t179 * t354 + t190 * (Icges(6,5) * t183 - Icges(6,6) * t184) + t191 * (Icges(6,5) * t185 + Icges(6,6) * t186);
t258 = qJD(3) * t374 + qJDD(3) * t372;
t257 = -qJDD(4) * t233 + t104 * t361 + t105 * t359 + t176 * t349 + t178 * t348 + t232 * t352 + qJDD(1);
t112 = Icges(6,5) * t233 + t232 * t306;
t110 = Icges(6,6) * t233 + t232 * t301;
t253 = t233 * t259;
t245 = (t254 * t238 + t312) * t190 + (t254 * t237 + t311) * t191;
t242 = (qJD(5) * t260 + t245) * t233;
t212 = t233 * t355;
t210 = t233 * t356;
t177 = t196 * t238;
t175 = t196 * t237;
t158 = t238 * t277;
t156 = t237 * t277;
t99 = -rSges(6,3) * t383 + t369;
t98 = -rSges(6,3) * t384 + t370;
t97 = t256 * t238;
t96 = t256 * t237;
t95 = t255 * t238;
t94 = t255 * t237;
t87 = rSges(6,1) * t185 + rSges(6,2) * t186;
t86 = rSges(6,1) * t183 - rSges(6,2) * t184;
t79 = t237 * t284 - t350;
t78 = t238 * t284 + t230;
t68 = qJD(3) * t112 + qJD(5) * t181;
t67 = qJD(3) * t110 + qJD(5) * t180;
t65 = t237 * t328 + t327;
t64 = t238 * t328 + t368;
t63 = qJD(3) * t289 + qJD(1);
t58 = Icges(6,1) * t102 + Icges(6,4) * t103 - Icges(6,5) * t339;
t57 = Icges(6,1) * t100 + Icges(6,4) * t101 - Icges(6,5) * t338;
t56 = Icges(6,4) * t102 + Icges(6,2) * t103 - Icges(6,6) * t339;
t55 = Icges(6,4) * t100 + Icges(6,2) * t101 - Icges(6,6) * t338;
t44 = -t233 * t295 - t385;
t43 = t237 * t258 + t310;
t42 = t238 * t258 + t340;
t41 = -qJD(3) * t287 + qJDD(3) * t289 + qJDD(1);
t40 = qJD(3) * t288 + t283;
t18 = t288 * qJDD(3) + (t156 * t237 + t158 * t238) * qJD(3) + t257;
t17 = (qJD(3) * t295 + t66) * t232 + (-qJD(3) * t254 - t239 * t68 - t240 * t67 + (-t239 * t255 + t240 * t256) * qJD(5)) * t233;
t16 = -t102 * t256 - t103 * t255 + t185 * t67 - t186 * t68 + t237 * t280;
t15 = -t100 * t256 - t101 * t255 + t183 * t67 + t184 * t68 + t238 * t280;
t12 = t190 * t24 + t191 * t25 + t354 * t44;
t11 = t102 * t75 + t103 * t73 + t185 * t56 - t186 * t58 + t237 * t281;
t10 = t102 * t74 + t103 * t72 + t185 * t55 - t186 * t57 + t237 * t282;
t9 = t100 * t75 + t101 * t73 + t183 * t56 + t184 * t58 + t238 * t281;
t8 = t100 * t74 + t101 * t72 + t183 * t55 + t184 * t57 + t238 * t282;
t7 = -pkin(6) * t241 * t329 + qJDD(3) * t286 + t106 * t77 - t107 * t76 + t190 * t60 - t191 * t59 + t257;
t6 = t190 * t22 + t191 * t23 + t354 * t39;
t5 = t190 * t20 + t191 * t21 + t354 * t38;
t4 = (qJD(3) * t311 + t54) * t232 + (qJD(3) * t71 - t239 * t58 - t240 * t56 + (t239 * t73 - t240 * t75) * qJD(5)) * t233;
t3 = (qJD(3) * t312 + t53) * t232 + (qJD(3) * t70 - t239 * t57 - t240 * t55 + (t239 * t72 - t240 * t74) * qJD(5)) * t233;
t2 = t10 * t190 + t106 * t22 + t107 * t23 + t11 * t191 + t16 * t354 + t189 * t39;
t1 = t106 * t20 + t107 * t21 + t15 * t354 + t189 * t38 + t190 * t8 + t191 * t9;
t26 = [(m(2) + m(3)) * qJDD(1) + m(4) * t41 + m(5) * t18 + m(6) * t7 + (-m(2) + t342) * g(3); t342 * (g(1) * t237 - g(2) * t238) + m(4) * (t237 * t78 - t238 * t79) + m(5) * (t237 * t42 - t238 * t43) + m(6) * t319 + m(3) * t365 * qJDD(2); -t12 * t353 / 0.2e1 + (t10 * t237 - t11 * t238) * t413 + ((t185 * t95 - t186 * t97) * t190 + (t185 * t94 - t186 * t96) * t191 + (t396 + (t110 * t185 - t112 * t186 - t398) * t232) * qJD(5) + (((-t23 + t385) * qJD(5) + t318) * t232 + t242) * t237) * t414 + (t237 * t8 - t238 * t9) * t415 + ((t183 * t95 + t184 * t97) * t190 + (t183 * t94 + t184 * t96) * t191 + (t397 + (t110 * t183 + t112 * t184 - t399) * t232) * qJD(5) + (((-t20 + t385) * qJD(5) + t318) * t232 + t242) * t238) * t416 + (t237 * t24 - t238 * t25) * t417 + (t22 * t237 - t23 * t238) * t419 + (t20 * t237 - t21 * t238) * t420 + (((-t239 * t97 - t240 * t95 + t70) * t190 + (-t239 * t96 - t240 * t94 + t71) * t191 + t44 * qJD(5)) * t233 + ((t260 + (-t110 * t240 - t112 * t239 - t254) * t233 - t315) * qJD(5) + t245) * t232) * t332 + (t435 * qJD(3) * t234 + (-t237 * t436 + t434) * t359) * t334 + ((-t238 * t435 + t434) * t361 + t436 * qJD(3) * t235) * t359 / 0.2e1 + ((-t4 + t5) * t238 + (t3 + t6) * t237) * t331 + (-g(1) * (t215 + t369) - g(2) * (t214 + t370) - g(3) * (pkin(6) * t233 + t115 + t423) - t422 * t232 * (-rSges(6,3) - pkin(3) - pkin(6)) + t7 * t373 + t19 * (-pkin(6) * t437 + t403) + (t13 * t324 + t36 * t285 + t7 * (t192 + t76) + t19 * t59) * t238 + (t14 * t324 + t37 * t285 + t7 * (t193 + t77) + t19 * t60) * t237 - t36 * (t115 * t191 + t212) - t37 * (-t115 * t190 + t210) - t19 * (t190 * t98 - t191 * t99 + t341) - (t314 * t423 + (-t19 * t329 + t233 * t314) * pkin(6)) * qJD(3) - ((-t36 * t77 + t37 * t76) * t233 + (t36 * (-t116 * t237 - t98) + t37 * (t116 * t238 + t99) + t279) * t232) * qJD(5)) * m(6) + (-g(1) * (t325 + t366) - g(2) * (t326 + t367) + g(3) * t371 - t64 * t212 - t65 * t210 - t40 * t341 - ((t366 * t40 + t371 * t64) * t238 + (t367 * t40 + t371 * t65) * t237) * qJD(3) + t18 * t373 + t40 * t403 + (t18 * t139 + t40 * t158 + t372 * t42 + t374 * t64) * t238 + (t18 * t140 + t40 * t156 + t372 * t43 + t374 * t65) * t237) * m(5) + (-(t63 * (-t175 * t237 - t177 * t238) + t290 * t199) * qJD(3) + t41 * t289 - t63 * t287 + (-t237 * t79 - t238 * t78) * t196 + t290 * t188 + g(1) * t177 + g(2) * t175 - g(3) * t199) * m(4) + (t1 + (t443 * t235 + (t445 * t237 - t238 * t446) * t237) * t440 + (t441 * t235 + (t447 * t237 - t238 * t448) * t237) * t439) * t237 / 0.2e1 - (t2 + (t446 * t235 + (t443 - t445) * t450) * t440 + (t448 * t235 + (t441 - t447) * t450) * t439) * t238 / 0.2e1; -t421 * g(3) * t233 + 0.2e1 * (t19 * t418 + t412 * t7) * m(6) + 0.2e1 * (t18 * t412 + t40 * t418) * m(5) + (t421 * t422 + m(5) * (qJD(3) * t40 + t237 * t43 + t238 * t42) + m(6) * (qJD(3) * t19 + t13 * t238 + t14 * t237)) * t232; -t5 * t338 / 0.2e1 + t1 * t381 / 0.2e1 + (t232 * t38 + t233 * t317) * t420 + (t15 * t232 + (t237 * t9 + t238 * t8) * t233 + (-t232 * t317 + t397) * qJD(3)) * t415 + t232 * t6 * t334 + t2 * t382 / 0.2e1 + (t232 * t39 + t233 * t316) * t419 + (t16 * t232 + (t10 * t238 + t11 * t237) * t233 + (-t232 * t316 + t396) * qJD(3)) * t413 + t12 * t362 / 0.2e1 + t232 * (t106 * t24 + t107 * t25 + t17 * t354 + t189 * t44 + t190 * t3 + t191 * t4) / 0.2e1 + (t232 * t44 + t233 * t315) * t417 + (t17 * t232 + (t237 * t4 + t238 * t3) * t233 + (-t232 * t315 + t44 * t233) * qJD(3)) * t331 + (t246 * t183 - t184 * t247 + t238 * t253) * t416 + (t185 * t246 + t186 * t247 + t237 * t253) * t414 + (t259 * t232 + (t247 * t239 - t240 * t246) * t233) * t332 + ((-t13 * t77 + t14 * t76 - t36 * t60 + t37 * t59 + (t279 + (-t237 * t36 + t238 * t37) * t116) * qJD(3)) * t232 + (t36 * (-qJD(3) * t77 + t237 * t69) + t37 * (qJD(3) * t76 - t238 * t69) - t7 * t313 + t19 * (-t237 * t59 + t238 * t60) + t319 * t116) * t233 - t36 * (t182 * t191 - t354 * t87) - t37 * (-t182 * t190 + t354 * t86) - t19 * (t190 * t87 - t191 * t86) - g(1) * t86 - g(2) * t87 - g(3) * t182) * m(6);];
tau = t26;
