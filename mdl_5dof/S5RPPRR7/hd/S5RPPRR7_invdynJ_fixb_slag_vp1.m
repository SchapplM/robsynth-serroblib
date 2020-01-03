% Calculate vector of inverse dynamics joint torques for
% S5RPPRR7
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR7_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:52
% DurationCPUTime: 17.60s
% Computational Cost: add. (14447->781), mult. (18887->1039), div. (0->0), fcn. (17498->8), ass. (0->367)
t285 = cos(qJ(1));
t278 = t285 * pkin(1);
t279 = qJ(1) + pkin(8);
t274 = sin(t279);
t275 = cos(t279);
t205 = rSges(3,1) * t274 + rSges(3,2) * t275;
t282 = sin(qJ(1));
t469 = pkin(1) * t282;
t188 = -t205 - t469;
t284 = cos(qJ(4));
t280 = sin(qJ(5));
t457 = rSges(6,2) * t280;
t283 = cos(qJ(5));
t460 = rSges(6,1) * t283;
t348 = -t457 + t460;
t281 = sin(qJ(4));
t455 = rSges(6,3) * t281;
t183 = t284 * t348 + t455;
t394 = qJD(5) * t284;
t370 = t275 * t394;
t399 = qJD(4) * t274;
t190 = t370 + t399;
t465 = t284 * pkin(4);
t242 = pkin(7) * t281 + t465;
t395 = qJD(5) * t281;
t261 = qJD(1) + t395;
t426 = t280 * t281;
t166 = t274 * t283 + t275 * t426;
t425 = t281 * t283;
t167 = -t274 * t280 + t275 * t425;
t417 = t167 * rSges(6,1) - t166 * rSges(6,2);
t427 = t275 * t284;
t95 = rSges(6,3) * t427 - t417;
t512 = t183 * t190 + t242 * t399 - t261 * t95;
t286 = qJD(1) ^ 2;
t511 = t286 * t278;
t159 = Icges(6,4) * t167;
t89 = Icges(6,2) * t166 + Icges(6,6) * t427 - t159;
t158 = Icges(6,4) * t166;
t91 = Icges(6,1) * t167 - Icges(6,5) * t427 - t158;
t345 = -t166 * t89 - t167 * t91;
t164 = -t274 * t426 + t275 * t283;
t165 = t274 * t425 + t275 * t280;
t429 = t274 * t284;
t438 = Icges(6,4) * t165;
t87 = Icges(6,2) * t164 - Icges(6,6) * t429 + t438;
t157 = Icges(6,4) * t164;
t90 = Icges(6,1) * t165 - Icges(6,5) * t429 + t157;
t461 = t164 * t87 + t165 * t90;
t84 = Icges(6,5) * t165 + Icges(6,6) * t164 - Icges(6,3) * t429;
t86 = -Icges(6,5) * t167 + Icges(6,6) * t166 + Icges(6,3) * t427;
t510 = t345 + t461 + (-t274 * t84 - t275 * t86) * t284;
t398 = qJD(4) * t275;
t191 = -t274 * t394 + t398;
t29 = t166 * t87 - t167 * t90 + t84 * t427;
t325 = Icges(6,5) * t283 - Icges(6,6) * t280;
t169 = Icges(6,3) * t281 + t284 * t325;
t436 = Icges(6,4) * t283;
t327 = -Icges(6,2) * t280 + t436;
t171 = Icges(6,6) * t281 + t284 * t327;
t437 = Icges(6,4) * t280;
t329 = Icges(6,1) * t283 - t437;
t173 = Icges(6,5) * t281 + t284 * t329;
t53 = t166 * t171 - t167 * t173 + t169 * t427;
t509 = -t191 * t29 - t261 * t53;
t207 = -rSges(4,2) * t275 + t274 * rSges(4,3);
t497 = t275 * pkin(2) + t274 * qJ(3);
t381 = t278 + t497;
t144 = t207 + t381;
t246 = pkin(7) * t427;
t507 = t417 - t246;
t331 = t280 * t89 + t283 * t91;
t36 = t281 * t86 - t284 * t331;
t28 = t164 * t89 - t165 * t91 - t429 * t86;
t396 = qJD(4) * t284;
t401 = qJD(1) * t281;
t504 = t274 * t396 + t275 * t401;
t439 = Icges(5,4) * t284;
t230 = -Icges(5,2) * t281 + t439;
t330 = Icges(5,1) * t281 + t439;
t412 = t230 + t330;
t440 = Icges(5,4) * t281;
t232 = Icges(5,1) * t284 - t440;
t328 = Icges(5,2) * t284 + t440;
t413 = -t328 + t232;
t503 = (t281 * t412 - t284 * t413) * qJD(1);
t52 = t164 * t171 + t165 * t173 - t169 * t429;
t501 = t190 * t28 + t52 * t261;
t470 = rSges(6,3) + pkin(7);
t500 = t281 * t470;
t168 = Icges(6,3) * t284 - t281 * t325;
t319 = t171 * t280 - t173 * t283;
t332 = t280 * t87 - t283 * t90;
t287 = t190 * (-t169 * t275 + t331) + t191 * (t169 * t274 + t332) + t261 * (t168 + t319);
t499 = t287 * t284;
t148 = -Icges(5,6) * t274 + t275 * t328;
t150 = -Icges(5,5) * t274 + t275 * t330;
t321 = t148 * t284 + t150 * t281;
t498 = t321 * t275;
t263 = t275 * qJ(3);
t468 = pkin(2) * t274;
t203 = -t263 + t468;
t354 = -pkin(6) * t274 - t469;
t314 = -t203 + t354;
t208 = t275 * rSges(3,1) - rSges(3,2) * t274;
t189 = t208 + t278;
t277 = t284 * pkin(7);
t467 = pkin(4) * t281;
t241 = t277 - t467;
t391 = qJD(1) * qJD(3);
t403 = qJD(1) * t274;
t259 = qJD(3) * t274;
t402 = qJD(1) * t275;
t410 = qJ(3) * t402 + t259;
t495 = qJD(1) * (-pkin(2) * t403 + t410) + qJDD(1) * t497 + t274 * t391;
t147 = Icges(5,6) * t275 + t274 * t328;
t177 = t230 * t275;
t100 = qJD(1) * t147 - qJD(4) * t177;
t179 = t232 * t275;
t435 = Icges(5,5) * t275;
t102 = -qJD(4) * t179 + (t274 * t330 + t435) * qJD(1);
t326 = Icges(5,5) * t281 + Icges(5,6) * t284;
t146 = -Icges(5,3) * t274 + t275 * t326;
t405 = qJD(1) * t146;
t76 = t148 * t281 - t150 * t284;
t494 = qJD(4) * t76 + t100 * t284 + t102 * t281 + t405;
t492 = t261 * t280 - t283 * t396;
t491 = t261 * t283 + t280 * t396;
t212 = t328 * qJD(4);
t213 = t330 * qJD(4);
t228 = Icges(5,5) * t284 - Icges(5,6) * t281;
t317 = t230 * t281 - t232 * t284;
t490 = qJD(1) * t228 + qJD(4) * t317 + t212 * t284 + t213 * t281;
t101 = qJD(1) * t148 + t230 * t399;
t178 = t232 * t274;
t103 = qJD(1) * t150 + qJD(4) * t178;
t237 = Icges(5,4) * t429;
t430 = t274 * t281;
t149 = Icges(5,1) * t430 + t237 + t435;
t322 = t147 * t281 - t149 * t284;
t145 = Icges(5,3) * t275 + t274 * t326;
t406 = qJD(1) * t145;
t489 = qJD(4) * t322 - t101 * t284 - t103 * t281 + t406;
t419 = t150 + t177;
t421 = t148 - t179;
t487 = t281 * t421 - t284 * t419;
t420 = -Icges(5,2) * t430 + t149 + t237;
t422 = t147 - t178;
t486 = t281 * t422 - t284 * t420;
t200 = (-Icges(6,2) * t283 - t437) * t284;
t290 = t190 * (Icges(6,2) * t167 + t158 - t91) + t191 * (-Icges(6,2) * t165 + t157 + t90) + t261 * (t173 + t200);
t201 = (-Icges(6,1) * t280 - t436) * t284;
t289 = t190 * (Icges(6,1) * t166 + t159 - t89) + t191 * (Icges(6,1) * t164 - t438 - t87) - t261 * (t171 - t201);
t485 = -pkin(2) - pkin(6);
t390 = qJD(1) * qJD(4);
t193 = qJDD(4) * t274 + t275 * t390;
t397 = qJD(4) * t281;
t374 = t275 * t397;
t400 = qJD(1) * t284;
t379 = t274 * t400;
t301 = -t374 - t379;
t389 = qJDD(5) * t284;
t110 = qJD(5) * t301 + t275 * t389 + t193;
t484 = t110 / 0.2e1;
t256 = qJDD(4) * t275;
t111 = -qJD(1) * t370 + t256 + (-t389 + (-qJD(1) + t395) * qJD(4)) * t274;
t483 = t111 / 0.2e1;
t482 = -t190 / 0.2e1;
t481 = t190 / 0.2e1;
t480 = -t191 / 0.2e1;
t479 = t191 / 0.2e1;
t478 = t193 / 0.2e1;
t194 = -t274 * t390 + t256;
t477 = t194 / 0.2e1;
t210 = qJD(4) * t394 + qJDD(5) * t281 + qJDD(1);
t476 = t210 / 0.2e1;
t475 = -t261 / 0.2e1;
t474 = t261 / 0.2e1;
t473 = t274 / 0.2e1;
t472 = -t275 / 0.2e1;
t466 = g(2) * t275;
t271 = t275 * pkin(6);
t376 = t274 * t397;
t377 = t275 * t400;
t302 = t376 - t377;
t359 = qJD(5) + t401;
t315 = t359 * t280;
t81 = -t274 * t491 - t275 * t315;
t316 = t283 * t359;
t82 = -t274 * t492 + t275 * t316;
t44 = Icges(6,5) * t82 + Icges(6,6) * t81 + Icges(6,3) * t302;
t46 = Icges(6,4) * t82 + Icges(6,2) * t81 + Icges(6,6) * t302;
t48 = Icges(6,1) * t82 + Icges(6,4) * t81 + Icges(6,5) * t302;
t7 = (qJD(4) * t332 + t44) * t281 + (qJD(4) * t84 - t280 * t46 + t283 * t48 + (-t280 * t90 - t283 * t87) * qJD(5)) * t284;
t464 = t7 * t191;
t79 = -t274 * t315 + t275 * t491;
t80 = t274 * t316 + t275 * t492;
t43 = Icges(6,5) * t80 + Icges(6,6) * t79 + Icges(6,3) * t301;
t45 = Icges(6,4) * t80 + Icges(6,2) * t79 + Icges(6,6) * t301;
t47 = Icges(6,1) * t80 + Icges(6,4) * t79 + Icges(6,5) * t301;
t8 = (qJD(4) * t331 + t43) * t281 + (qJD(4) * t86 - t280 * t45 + t283 * t47 + (t280 * t91 - t283 * t89) * qJD(5)) * t284;
t463 = t8 * t190;
t199 = (-Icges(6,5) * t280 - Icges(6,6) * t283) * t284;
t124 = qJD(4) * t168 + qJD(5) * t199;
t170 = Icges(6,6) * t284 - t281 * t327;
t125 = qJD(4) * t170 + qJD(5) * t200;
t172 = Icges(6,5) * t284 - t281 * t329;
t126 = qJD(4) * t172 + qJD(5) * t201;
t24 = (qJD(4) * t319 + t124) * t281 + (qJD(4) * t169 - t125 * t280 + t126 * t283 + (-t171 * t283 - t173 * t280) * qJD(5)) * t284;
t66 = t169 * t281 - t284 * t319;
t462 = t66 * t210 + t24 * t261;
t456 = rSges(4,3) * t275;
t235 = rSges(5,1) * t284 - rSges(5,2) * t281;
t181 = t235 * t275;
t269 = t275 * rSges(5,3);
t350 = rSges(5,1) * t281 + rSges(5,2) * t284;
t106 = -qJD(4) * t181 + (t274 * t350 + t269) * qJD(1);
t260 = qJD(3) * t275;
t154 = qJD(1) * t497 - t260;
t214 = t350 * qJD(4);
t353 = -t278 - t271;
t411 = qJDD(3) * t274 + t275 * t391;
t292 = t286 * t353 + t411;
t152 = -t274 * rSges(5,3) + t275 * t350;
t305 = t152 + t314;
t25 = -t214 * t399 + t193 * t235 + (-t106 - t154) * qJD(1) + t305 * qJDD(1) + t292;
t453 = t25 * t274;
t383 = rSges(5,1) * t504 + rSges(5,2) * t377;
t107 = (-rSges(5,2) * t397 - rSges(5,3) * qJD(1)) * t274 + t383;
t151 = rSges(5,1) * t430 + rSges(5,2) * t429 + t269;
t273 = qJDD(1) * t278;
t291 = qJDD(1) * t271 + t286 * t354 + t273 + t495;
t26 = qJD(1) * t107 + qJDD(1) * t151 - t194 * t235 + (qJD(4) * t214 - qJDD(3)) * t275 + t291;
t452 = t26 * t275;
t192 = t235 * t399;
t64 = qJD(1) * t305 + t192 + t259;
t449 = t275 * t64;
t276 = t284 * rSges(6,3);
t35 = t281 * t84 - t284 * t332;
t448 = t35 * t111;
t447 = t36 * t110;
t244 = pkin(4) * t430;
t184 = -pkin(7) * t429 + t244;
t418 = t165 * rSges(6,1) + t164 * rSges(6,2);
t93 = -rSges(6,3) * t429 + t418;
t442 = -t184 - t93;
t428 = t275 * t281;
t186 = pkin(4) * t428 - t246;
t441 = t186 - t95;
t431 = t228 * t275;
t174 = t274 * t228;
t202 = (-rSges(6,1) * t280 - rSges(6,2) * t283) * t284;
t127 = qJD(5) * t202 + (-t281 * t348 + t276) * qJD(4);
t217 = t241 * qJD(4);
t423 = t127 + t217;
t414 = t183 + t242;
t185 = pkin(4) * t429 + pkin(7) * t430;
t409 = rSges(4,2) * t403 + rSges(4,3) * t402;
t408 = rSges(6,2) * t426 + t276;
t407 = -qJD(1) * t203 + t259;
t404 = qJD(1) * t326;
t318 = t284 * t230 + t232 * t281;
t105 = t275 * t318 - t174;
t393 = t105 * qJD(1);
t392 = -m(4) - m(5) - m(6);
t388 = -rSges(5,3) + t485;
t387 = t82 * rSges(6,1) + t81 * rSges(6,2) + rSges(6,3) * t376;
t386 = t284 * t460;
t385 = t284 * t457;
t56 = t275 * t145 + t147 * t429 + t149 * t430;
t57 = -t275 * t146 - t148 * t429 - t150 * t430;
t382 = pkin(4) * t504 + pkin(7) * t376;
t380 = t470 * t284;
t373 = t275 * t396;
t369 = -pkin(4) - t460;
t368 = -t403 / 0.2e1;
t366 = -t400 / 0.2e1;
t365 = -t399 / 0.2e1;
t364 = t399 / 0.2e1;
t363 = -t398 / 0.2e1;
t362 = t398 / 0.2e1;
t361 = rSges(4,2) * t274 + t456 - t469;
t360 = t263 - t469;
t358 = t271 + t381;
t357 = -t286 * t469 + t273;
t355 = rSges(6,1) * t80 + rSges(6,2) * t79;
t352 = -t276 + t467;
t236 = rSges(2,1) * t285 - rSges(2,2) * t282;
t234 = rSges(2,1) * t282 + rSges(2,2) * t285;
t323 = t147 * t284 + t149 * t281;
t294 = qJD(1) * t323 + qJD(4) * t174 + t405;
t295 = -qJD(1) * t321 - qJD(4) * t431 + t406;
t344 = (t294 * t274 + t275 * t489) * t275 + (t295 * t274 - t275 * t494) * t274;
t343 = (-t274 * t489 + t294 * t275) * t275 + (t274 * t494 + t295 * t275) * t274;
t27 = -t429 * t84 + t461;
t342 = t27 * t275 + t274 * t28;
t341 = t27 * t274 - t275 * t28;
t30 = t427 * t86 - t345;
t340 = t274 * t30 + t275 * t29;
t339 = t274 * t29 - t275 * t30;
t338 = t274 * t36 + t275 * t35;
t337 = t274 * t35 - t275 * t36;
t336 = t274 * t57 + t275 * t56;
t140 = t274 * t145;
t58 = -t323 * t275 + t140;
t59 = -t146 * t274 + t498;
t335 = t274 * t59 + t275 * t58;
t313 = t497 - t353;
t65 = -t235 * t398 - t260 + (t151 + t313) * qJD(1);
t334 = t274 * t64 - t275 * t65;
t333 = t274 * t95 + t275 * t93;
t324 = t106 * t275 - t107 * t274;
t320 = -t151 * t274 - t152 * t275;
t310 = -t284 * t43 + t397 * t86;
t309 = -t284 * t44 + t397 * t84;
t138 = rSges(6,3) * t430 + (-t385 + t386) * t274;
t304 = t186 + t314;
t303 = t485 * t274 + t360;
t300 = t169 * t261 + t190 * t86 + t191 * t84;
t299 = (Icges(6,5) * t164 - Icges(6,6) * t165) * t191 + (Icges(6,5) * t166 + Icges(6,6) * t167) * t190 + t199 * t261;
t293 = t318 * qJD(1) - t326 * qJD(4);
t32 = -t190 * t93 + t191 * t95 + qJD(2) + (-t184 * t274 - t186 * t275) * qJD(4);
t39 = qJD(1) * t304 + t259 + t512;
t40 = -t242 * t398 - t183 * t191 + t261 * t93 - t260 + (t184 + t313) * qJD(1);
t288 = t32 * t333 + (-t274 * t40 - t275 * t39) * t183;
t219 = t275 * t385;
t187 = t242 * t275;
t182 = -rSges(6,1) * t425 + t408;
t180 = t235 * t274;
t139 = t219 + (-t386 - t455) * t275;
t137 = t173 * t275;
t136 = t173 * t274;
t135 = t171 * t275;
t134 = t171 * t274;
t123 = -pkin(7) * t377 + t382;
t122 = t301 * pkin(7) + (t274 * t401 - t373) * pkin(4);
t121 = rSges(6,1) * t166 + rSges(6,2) * t167;
t120 = rSges(6,1) * t164 - rSges(6,2) * t165;
t104 = t274 * t318 + t431;
t83 = t104 * qJD(1);
t73 = qJD(4) * t320 + qJD(2);
t55 = qJD(1) * t409 + qJDD(1) * t207 - qJDD(3) * t275 + t357 + t495;
t54 = -t511 + (-t203 + t361) * qJDD(1) + (-qJD(1) * t207 - t154) * qJD(1) + t411;
t50 = -rSges(6,3) * t377 + t387;
t49 = rSges(6,3) * t301 + t355;
t42 = -t274 * t490 + t293 * t275;
t41 = t293 * t274 + t275 * t490;
t38 = t321 * qJD(4) - t100 * t281 + t102 * t284;
t37 = -qJD(4) * t323 - t101 * t281 + t103 * t284;
t31 = qJD(4) * t324 - t151 * t193 - t152 * t194 + qJDD(2);
t23 = qJD(4) * t335 - t393;
t22 = qJD(4) * t336 + t83;
t16 = -t124 * t429 + t125 * t164 + t126 * t165 + t169 * t302 + t171 * t81 + t173 * t82;
t15 = t124 * t427 + t125 * t166 - t126 * t167 + t169 * t301 + t171 * t79 + t173 * t80;
t14 = t190 * t36 + t191 * t35 + t261 * t66;
t13 = qJD(1) * t123 + qJDD(1) * t184 - t111 * t183 - t127 * t191 - t194 * t242 + t210 * t93 + t261 * t50 + (-qJD(4) * t217 - qJDD(3)) * t275 + t291;
t12 = t217 * t399 + t110 * t183 + t127 * t190 + t193 * t242 - t210 * t95 - t261 * t49 + (-t122 - t154) * qJD(1) + t304 * qJDD(1) + t292;
t11 = t190 * t30 - t509;
t10 = t191 * t27 + t501;
t9 = -t110 * t93 + t111 * t95 - t184 * t193 - t186 * t194 - t190 * t50 + t191 * t49 + qJDD(2) + (t122 * t275 - t123 * t274) * qJD(4);
t6 = t164 * t45 + t165 * t47 + t274 * t310 - t377 * t86 + t81 * t89 - t82 * t91;
t5 = t164 * t46 + t165 * t48 + t274 * t309 - t377 * t84 + t81 * t87 + t82 * t90;
t4 = t166 * t45 - t167 * t47 - t275 * t310 - t379 * t86 + t79 * t89 - t80 * t91;
t3 = t166 * t46 - t167 * t48 - t275 * t309 - t379 * t84 + t79 * t87 + t80 * t90;
t2 = t110 * t28 + t111 * t27 + t16 * t261 + t190 * t6 + t191 * t5 + t210 * t52;
t1 = t110 * t30 + t111 * t29 + t15 * t261 + t190 * t4 + t191 * t3 + t210 * t53;
t17 = [(t38 + t41 + t22) * t364 + t448 / 0.2e1 + (t15 + t10) * t481 + (t37 + t42) * t362 + t447 / 0.2e1 + ((-t205 * t286 - g(2) + t357) * t189 + (-t511 + (-0.2e1 * t208 - t278 + t189) * t286 - g(1)) * t188) * m(3) + t462 + t52 * t483 + t53 * t484 + t76 * t478 + t16 * t479 + ((t30 + t510) * t191 + t501) * t482 + (t11 + (-t27 + t510) * t190 + t509) * t480 + (-(-t39 + t407 + (t186 + t354) * qJD(1) + t512) * t40 - g(1) * (t275 * t352 + t303 + t507) + t12 * (t314 + t507) + t39 * (t260 - t355) + t40 * (t382 + t387 + t410) + (t12 * t352 + t39 * (t465 + t500) * qJD(4)) * t275 + ((-t282 * t40 - t285 * t39) * pkin(1) + (-t380 * t40 + t39 * t485) * t275 + (t39 * (-qJ(3) - t352 + t277) + t40 * t485) * t274) * qJD(1) + (-g(2) + t13) * (-t274 * t380 + t244 + t358 + t418)) * m(6) + (t23 + t393 + (t146 * t274 ^ 2 + (-t140 + t57 + (t146 + t323) * t275) * t275) * qJD(4)) * t363 + t463 / 0.2e1 + t464 / 0.2e1 + ((t54 - g(1)) * (t456 + (rSges(4,2) - pkin(2)) * t274 + t360) + (t55 - g(2)) * t144 + (-t407 + t409 + t410 + (-t361 - t468 - t469) * qJD(1)) * (qJD(1) * t144 - t260)) * m(4) + (-t322 + t104) * t477 + (-t317 + m(2) * (t234 ^ 2 + t236 ^ 2) + m(3) * (t188 ^ 2 + t208 * t189) + Icges(2,3) + Icges(3,3) + Icges(4,1)) * qJDD(1) + (-qJD(4) * t318 + t212 * t281 - t213 * t284) * qJD(1) + (t83 + ((-t58 + t140 + t57) * t274 + (t59 - t498 + (t146 - t323) * t274 + t56) * t275) * qJD(4)) * t365 + (t64 * (rSges(5,1) * t373 - rSges(5,2) * t374 + t260) + t65 * (-rSges(5,2) * t376 + t383 + t410) + ((-t282 * t65 - t285 * t64) * pkin(1) + t388 * t449 + (t64 * (-qJ(3) - t350) + t65 * t388) * t274) * qJD(1) - (t192 - t64 + (t152 + t354) * qJD(1) + t407) * t65 + (t26 - g(2)) * (t358 + t151) + (t25 - g(1)) * (t152 + t303)) * m(5) - t193 * t105 / 0.2e1 - m(2) * (-g(1) * t234 + g(2) * t236); (m(3) + m(4)) * qJDD(2) + m(5) * t31 + m(6) * t9 + (-m(3) + t392) * g(3); t392 * (g(1) * t274 - t466) + 0.2e1 * (t12 * t473 + t13 * t472) * m(6) + 0.2e1 * (t453 / 0.2e1 - t452 / 0.2e1) * m(5) + 0.2e1 * (t472 * t55 + t473 * t54) * m(4); -qJD(1) * ((-t281 * t413 - t284 * t412) * qJD(1) + ((t274 * t421 - t275 * t422) * t284 + (t274 * t419 - t275 * t420) * t281) * qJD(4)) / 0.2e1 + qJD(1) * (t274 * t38 + t275 * t37 + (t274 * t322 + t76 * t275) * qJD(1)) / 0.2e1 + qJDD(1) * (t274 * t76 - t275 * t322) / 0.2e1 + (qJD(1) * t41 + qJD(4) * t344 - qJDD(1) * t105 + t193 * t59 + t194 * t58 + t1) * t473 + (t31 * t320 + t73 * ((-t151 * t275 + t152 * t274) * qJD(1) + t324) - t334 * t214 + (t453 - t452 + (t274 * t65 + t449) * qJD(1)) * t235 - (t180 * t65 + t181 * t64) * qJD(1) - (t73 * (-t180 * t274 - t181 * t275) - t334 * t350) * qJD(4) - g(1) * t180 + g(2) * t181 + g(3) * t350) * m(5) + (((-t134 * t280 + t136 * t283 + t84) * t191 + (t135 * t280 - t137 * t283 + t86) * t190 + (-t170 * t280 + t172 * t283 + t169) * t261 + t66 * qJD(5)) * t284 + (qJD(5) * t337 + t287) * t281) * t475 + ((t134 * t164 + t136 * t165) * t191 + (-t135 * t164 - t137 * t165) * t190 + (t164 * t170 + t165 * t172) * t261 + (-t28 * t428 + t284 * t52) * qJD(5) + ((qJD(5) * t27 + t300) * t281 - t499) * t274) * t480 + ((t134 * t166 - t136 * t167) * t191 + (-t135 * t166 + t137 * t167) * t190 + (t166 * t170 - t167 * t172) * t261 + (t284 * t53 + t29 * t430) * qJD(5) + ((-qJD(5) * t30 - t300) * t281 + t499) * t275) * t482 + (t23 + t11) * t402 / 0.2e1 + (qJD(1) * t42 + qJD(4) * t343 + qJDD(1) * t104 + t11 * t395 + t193 * t57 + t194 * t56 + t2) * t275 / 0.2e1 + (t368 - t274 * t395 / 0.2e1) * t10 + t22 * t368 + ((-t58 * t274 + t59 * t275) * qJD(1) + t344) * t364 + ((-t56 * t274 + t57 * t275) * qJD(1) + t343) * t362 + (-qJD(1) * t339 + t274 * t4 + t275 * t3) * t481 + t342 * t483 + t340 * t484 - t14 * t394 / 0.2e1 + (-qJD(1) * t337 + t274 * t8 + t275 * t7) * t474 + t338 * t476 + t336 * t477 + t335 * t478 + (-qJD(1) * t341 + t274 * t6 + t275 * t5) * t479 + ((t174 * t398 - t404) * t275 + (-t503 + (t487 * t274 + (-t431 - t486) * t275) * qJD(4)) * t274) * t363 + ((-t399 * t431 - t404) * t274 + (t503 + (t486 * t275 + (t174 - t487) * t274) * qJD(4)) * t275) * t365 + (-t39 * (qJD(1) * t187 - t139 * t261 + t182 * t190 + t241 * t399) - t40 * (qJD(1) * t185 + t138 * t261 - t182 * t191 - t241 * t398) - t32 * (-t138 * t190 + t139 * t191 - t185 * t399 - t187 * t398) - ((-t39 * t95 + t40 * t93) * t284 + t288 * t281) * qJD(5) + (-t13 * t414 - t40 * t423 - t9 * t441 + t32 * (t122 + t49) + (t32 * t442 + t39 * t414) * qJD(1)) * t275 + (t12 * t414 + t39 * t423 + t9 * t442 + t32 * (-t123 - t50) + (t32 * t441 + t40 * t414) * qJD(1)) * t274 - g(1) * (t138 + t185) - g(2) * t219 - g(3) * (t281 * t369 + t277 + t408) - (t284 * t369 - t500) * t466) * m(6); -t2 * t429 / 0.2e1 + (t281 * t52 - t284 * t341) * t483 + ((qJD(4) * t341 + t16) * t281 + (-qJD(1) * t342 + qJD(4) * t52 - t274 * t5 + t275 * t6) * t284) * t479 + t1 * t427 / 0.2e1 + (t281 * t53 - t284 * t339) * t484 + ((qJD(4) * t339 + t15) * t281 + (-qJD(1) * t340 + qJD(4) * t53 - t274 * t3 + t275 * t4) * t284) * t481 + t14 * t396 / 0.2e1 + t281 * (t447 + t448 + t462 + t463 + t464) / 0.2e1 + (t281 * t66 - t284 * t337) * t476 + ((qJD(4) * t337 + t24) * t281 + (-qJD(1) * t338 + qJD(4) * t66 - t274 * t7 + t275 * t8) * t284) * t474 + (t164 * t290 + t165 * t289 - t299 * t429) * t480 + (t290 * t166 - t167 * t289 + t299 * t427) * t482 + (t299 * t281 + (-t280 * t290 + t289 * t283) * t284) * t475 + (t274 * t366 + t281 * t363) * t11 + (t275 * t366 + t281 * t364) * t10 + ((qJD(4) * t288 - t12 * t95 + t13 * t93 - t39 * t49 + t40 * t50) * t281 + (t39 * (-qJD(4) * t95 + t127 * t275) + t40 * (qJD(4) * t93 + t127 * t274) - t9 * t333 + t32 * (-t274 * t49 - t275 * t50 - t402 * t95 + t403 * t93) + (t12 * t275 + t13 * t274 + (-t274 * t39 + t275 * t40) * qJD(1)) * t183) * t284 - t39 * (-t121 * t261 + t190 * t202) - t40 * (t120 * t261 - t191 * t202) - t32 * (-t120 * t190 + t121 * t191) - g(1) * t120 - g(2) * t121 - g(3) * t202) * m(6);];
tau = t17;
