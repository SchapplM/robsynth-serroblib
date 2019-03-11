% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPPPRR2_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR2_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:31:39
% EndTime: 2019-03-09 01:31:59
% DurationCPUTime: 17.64s
% Computational Cost: add. (20871->755), mult. (19783->1013), div. (0->0), fcn. (17987->10), ass. (0->377)
t277 = qJ(1) + pkin(9);
t274 = cos(t277);
t259 = qJD(3) * t274;
t272 = sin(t277);
t369 = qJD(4) * t272 - t259;
t276 = pkin(10) + qJ(5);
t273 = cos(t276);
t271 = sin(t276);
t479 = rSges(6,2) * t271;
t221 = rSges(6,1) * t273 - t479;
t411 = qJD(5) * t274;
t390 = t221 * t411;
t539 = t369 - t390;
t492 = pkin(2) * t272;
t283 = cos(qJ(6));
t448 = t274 * t283;
t281 = sin(qJ(6));
t453 = t272 * t281;
t192 = -t271 * t453 + t448;
t449 = t274 * t281;
t452 = t272 * t283;
t193 = t271 * t452 + t449;
t455 = t272 * t273;
t100 = Icges(7,5) * t193 + Icges(7,6) * t192 - Icges(7,3) * t455;
t194 = t271 * t449 + t452;
t195 = t271 * t448 - t453;
t451 = t273 * t274;
t102 = -Icges(7,5) * t195 + Icges(7,6) * t194 + Icges(7,3) * t451;
t180 = Icges(7,4) * t195;
t105 = Icges(7,2) * t194 + Icges(7,6) * t451 - t180;
t179 = Icges(7,4) * t194;
t107 = Icges(7,1) * t195 - Icges(7,5) * t451 - t179;
t334 = -t105 * t194 - t195 * t107;
t467 = Icges(7,4) * t193;
t103 = Icges(7,2) * t192 - Icges(7,6) * t455 + t467;
t178 = Icges(7,4) * t192;
t106 = Icges(7,1) * t193 - Icges(7,5) * t455 + t178;
t483 = t103 * t192 + t106 * t193;
t538 = t334 + t483 + (-t100 * t272 - t102 * t274) * t273;
t280 = -pkin(7) - qJ(4);
t446 = qJ(4) + t280;
t278 = sin(pkin(10));
t450 = t274 * t278;
t170 = pkin(4) * t450 + t272 * t446;
t257 = qJD(4) * t274;
t261 = t274 * qJ(3);
t218 = -t261 + t492;
t258 = qJD(3) * t272;
t421 = -qJD(1) * t218 + t258;
t395 = t257 + t421;
t415 = qJD(1) * t274;
t253 = qJ(3) * t415;
t422 = t257 + t258;
t396 = t253 + t422;
t391 = t278 * t415;
t414 = qJD(1) * t280;
t425 = pkin(4) * t391 + t272 * t414;
t537 = -qJD(1) * t170 - t395 + t396 + t425;
t282 = sin(qJ(1));
t493 = pkin(1) * t282;
t353 = -qJ(4) * t272 - t493;
t536 = -t353 - t493;
t409 = qJD(6) * t273;
t203 = -t272 * t409 + t411;
t410 = qJD(6) * t271;
t248 = qJD(1) + t410;
t30 = t100 * t451 + t103 * t194 - t106 * t195;
t336 = Icges(7,5) * t283 - Icges(7,6) * t281;
t162 = Icges(7,3) * t271 + t273 * t336;
t465 = Icges(7,4) * t283;
t338 = -Icges(7,2) * t281 + t465;
t164 = Icges(7,6) * t271 + t273 * t338;
t466 = Icges(7,4) * t281;
t340 = Icges(7,1) * t283 - t466;
t166 = Icges(7,5) * t271 + t273 * t340;
t54 = t162 * t451 + t164 * t194 - t166 * t195;
t535 = -t203 * t30 - t248 * t54;
t370 = -rSges(4,2) * t274 + rSges(4,3) * t272;
t284 = cos(qJ(1));
t275 = t284 * pkin(1);
t270 = t274 * pkin(2);
t518 = qJ(3) * t272 + t270;
t394 = t275 + t518;
t530 = t370 + t394;
t333 = t105 * t281 + t107 * t283;
t41 = t102 * t271 - t273 * t333;
t29 = -t102 * t455 + t105 * t192 - t107 * t193;
t285 = qJD(1) ^ 2;
t412 = qJD(5) * t273;
t529 = t271 * t415 + t272 * t412;
t233 = Icges(6,4) * t455;
t457 = t271 * t272;
t464 = Icges(6,5) * t274;
t154 = Icges(6,1) * t457 + t233 + t464;
t468 = Icges(6,4) * t273;
t341 = Icges(6,1) * t271 + t468;
t155 = -Icges(6,5) * t272 + t274 * t341;
t214 = -Icges(6,2) * t271 + t468;
t175 = t214 * t274;
t297 = t272 * (t155 + t175) - t274 * (-Icges(6,2) * t457 + t154 + t233);
t469 = Icges(6,4) * t271;
t339 = Icges(6,2) * t273 + t469;
t152 = Icges(6,6) * t274 + t272 * t339;
t153 = -Icges(6,6) * t272 + t274 * t339;
t216 = Icges(6,1) * t273 - t469;
t176 = t216 * t272;
t177 = t216 * t274;
t298 = t272 * (t153 - t177) - t274 * (t152 - t176);
t528 = -t298 * t271 + t297 * t273;
t427 = t214 + t341;
t428 = -t339 + t216;
t527 = (t271 * t427 - t273 * t428) * qJD(1);
t356 = rSges(7,1) * t195 - rSges(7,2) * t194;
t399 = rSges(7,3) * t451;
t111 = -t356 + t399;
t355 = rSges(7,1) * t283 - rSges(7,2) * t281;
t168 = rSges(7,3) * t271 + t273 * t355;
t408 = qJD(6) * t274;
t413 = qJD(5) * t272;
t202 = t273 * t408 + t413;
t488 = t273 * pkin(5);
t224 = pkin(8) * t271 + t488;
t526 = -t111 * t248 + t168 * t202 + t224 * t413;
t524 = 0.2e1 * qJD(5);
t53 = -t162 * t455 + t164 * t192 + t166 * t193;
t523 = t202 * t29 + t248 * t53;
t161 = Icges(7,3) * t273 - t271 * t336;
t327 = t164 * t281 - t166 * t283;
t335 = t103 * t281 - t106 * t283;
t286 = t202 * (-t162 * t274 + t333) + t203 * (t162 * t272 + t335) + t248 * (t161 + t327);
t522 = t286 * t273;
t329 = t153 * t273 + t155 * t271;
t519 = t329 * t274;
t371 = rSges(3,1) * t274 - rSges(3,2) * t272;
t517 = t275 + t371;
t454 = t272 * t278;
t480 = rSges(5,2) * cos(pkin(10));
t515 = rSges(5,1) * t454 + rSges(5,3) * t274 + t272 * t480;
t337 = Icges(6,5) * t271 + Icges(6,6) * t273;
t151 = -Icges(6,3) * t272 + t274 * t337;
t419 = qJD(1) * t151;
t79 = t153 * t271 - t155 * t273;
t94 = qJD(1) * t152 - qJD(5) * t175;
t96 = -qJD(5) * t177 + (t272 * t341 + t464) * qJD(1);
t514 = qJD(5) * t79 + t271 * t96 + t273 * t94 + t419;
t512 = t248 * t281 - t283 * t412;
t511 = t248 * t283 + t281 * t412;
t206 = t339 * qJD(5);
t207 = t341 * qJD(5);
t212 = Icges(6,5) * t273 - Icges(6,6) * t271;
t510 = qJD(1) * t212 + qJD(5) * (t214 * t271 - t216 * t273) + t206 * t273 + t207 * t271;
t330 = t152 * t271 - t154 * t273;
t150 = Icges(6,3) * t274 + t272 * t337;
t420 = qJD(1) * t150;
t95 = qJD(1) * t153 + t214 * t413;
t97 = qJD(1) * t155 + qJD(5) * t176;
t509 = qJD(5) * t330 - t271 * t97 - t273 * t95 + t420;
t200 = (-Icges(7,2) * t283 - t466) * t273;
t290 = t202 * (Icges(7,2) * t195 - t107 + t179) + t203 * (-Icges(7,2) * t193 + t106 + t178) + t248 * (t166 + t200);
t201 = (-Icges(7,1) * t281 - t465) * t273;
t507 = t202 * (-Icges(7,1) * t194 + t105 - t180) + t203 * (-Icges(7,1) * t192 + t103 + t467) + t248 * (t164 - t201);
t405 = qJD(5) * qJD(6);
t384 = t271 * t405;
t147 = -qJD(1) * t202 + t272 * t384;
t506 = t147 / 0.2e1;
t148 = qJD(1) * t203 - t274 * t384;
t505 = t148 / 0.2e1;
t504 = -t202 / 0.2e1;
t503 = t202 / 0.2e1;
t502 = -t203 / 0.2e1;
t501 = t203 / 0.2e1;
t500 = -t248 / 0.2e1;
t499 = t248 / 0.2e1;
t498 = t271 / 0.2e1;
t497 = t272 / 0.2e1;
t496 = -t274 / 0.2e1;
t494 = rSges(7,3) + pkin(8);
t491 = pkin(4) * t278;
t490 = pkin(5) * t271;
t489 = pkin(8) * t273;
t388 = t271 * t413;
t392 = t273 * t415;
t307 = t388 - t392;
t368 = qJD(1) * t271 + qJD(6);
t324 = t281 * t368;
t90 = -t272 * t511 - t274 * t324;
t323 = t283 * t368;
t91 = -t272 * t512 + t274 * t323;
t45 = Icges(7,5) * t91 + Icges(7,6) * t90 + Icges(7,3) * t307;
t47 = Icges(7,4) * t91 + Icges(7,2) * t90 + Icges(7,6) * t307;
t49 = Icges(7,1) * t91 + Icges(7,4) * t90 + Icges(7,5) * t307;
t7 = (qJD(5) * t335 + t45) * t271 + (qJD(5) * t100 - t281 * t47 + t283 * t49 + (-t103 * t283 - t106 * t281) * qJD(6)) * t273;
t487 = t7 * t203;
t416 = qJD(1) * t273;
t306 = -t271 * t411 - t272 * t416;
t88 = -t272 * t324 + t274 * t511;
t89 = t272 * t323 + t274 * t512;
t44 = Icges(7,5) * t89 + Icges(7,6) * t88 + Icges(7,3) * t306;
t46 = Icges(7,4) * t89 + Icges(7,2) * t88 + Icges(7,6) * t306;
t48 = Icges(7,1) * t89 + Icges(7,4) * t88 + Icges(7,5) * t306;
t8 = (qJD(5) * t333 + t44) * t271 + (qJD(5) * t102 - t281 * t46 + t283 * t48 + (-t105 * t283 + t107 * t281) * qJD(6)) * t273;
t486 = t8 * t202;
t485 = -qJD(1) / 0.2e1;
t199 = (-Icges(7,5) * t281 - Icges(7,6) * t283) * t273;
t116 = qJD(5) * t161 + qJD(6) * t199;
t163 = Icges(7,6) * t273 - t271 * t338;
t117 = qJD(5) * t163 + qJD(6) * t200;
t165 = Icges(7,5) * t273 - t271 * t340;
t118 = qJD(5) * t165 + qJD(6) * t201;
t22 = (qJD(5) * t327 + t116) * t271 + (qJD(5) * t162 - t117 * t281 + t118 * t283 + (-t164 * t283 - t166 * t281) * qJD(6)) * t273;
t383 = t273 * t405;
t65 = t162 * t271 - t273 * t327;
t484 = t22 * t248 + t65 * t383;
t478 = rSges(4,3) * t274;
t476 = rSges(7,3) * t273;
t357 = rSges(6,1) * t271 + rSges(6,2) * t273;
t208 = t357 * qJD(5);
t407 = qJD(1) * qJD(3);
t252 = t274 * t407;
t460 = qJ(4) * t274;
t352 = -t460 - t275;
t406 = qJD(1) * qJD(4);
t295 = -0.2e1 * t272 * t406 + t285 * t352 + t252;
t181 = qJD(1) * t518 - t259;
t242 = t274 * t414;
t247 = pkin(4) * t454;
t434 = t242 - (t247 - t460) * qJD(1) - t181;
t183 = t221 * t274;
t267 = t274 * rSges(6,3);
t98 = -qJD(5) * t183 + (t272 * t357 + t267) * qJD(1);
t43 = -t208 * t413 + (-t98 + t390 + t434) * qJD(1) + t295;
t473 = t272 * t43;
t196 = t221 * t413;
t417 = qJD(1) * t272;
t424 = t253 + t258;
t430 = qJD(1) * (-pkin(2) * t417 + t424) + t272 * t407;
t294 = 0.2e1 * t274 * t406 + t285 * t353 + t430;
t293 = qJD(1) * (qJ(4) * t417 + t425) + t294;
t398 = rSges(6,1) * t529 + rSges(6,2) * t392;
t99 = (-rSges(6,3) * qJD(1) - qJD(5) * t479) * t272 + t398;
t42 = t208 * t411 + (t99 + t196) * qJD(1) + t293;
t472 = t274 * t42;
t40 = t100 * t271 - t273 * t335;
t471 = t40 * t147;
t470 = t41 * t148;
t458 = t212 * t274;
t456 = t271 * t274;
t172 = t272 * t212;
t326 = t214 * t273 + t216 * t271;
t85 = t274 * t326 - t172;
t447 = t85 * qJD(1);
t429 = rSges(7,1) * t193 + rSges(7,2) * t192;
t109 = -rSges(7,3) * t455 + t429;
t236 = pkin(5) * t457;
t186 = -pkin(8) * t455 + t236;
t441 = -t109 - t186;
t237 = pkin(8) * t451;
t188 = pkin(5) * t456 - t237;
t440 = t111 - t188;
t167 = -t271 * t355 + t476;
t204 = (-rSges(7,1) * t281 - rSges(7,2) * t283) * t273;
t121 = qJD(5) * t167 + qJD(6) * t204;
t223 = t489 - t490;
t209 = qJD(5) * t223;
t439 = t121 + t209;
t431 = t168 + t224;
t400 = t274 * t480;
t426 = rSges(5,1) * t391 + qJD(1) * t400;
t423 = rSges(4,2) * t417 + rSges(4,3) * t415;
t418 = qJD(1) * t337;
t404 = -rSges(5,3) - pkin(2) - qJ(4);
t403 = t285 * t493;
t402 = t285 * t275;
t401 = rSges(7,1) * t91 + rSges(7,2) * t90 + rSges(7,3) * t388;
t60 = t150 * t274 + t152 * t455 + t154 * t457;
t61 = -t151 * t274 - t153 * t455 - t155 * t457;
t397 = pkin(5) * t529 + pkin(8) * t388;
t156 = rSges(6,1) * t457 + rSges(6,2) * t455 + t267;
t389 = t224 * t411;
t381 = -t416 / 0.2e1;
t379 = -t413 / 0.2e1;
t378 = t413 / 0.2e1;
t377 = t412 / 0.2e1;
t376 = -t411 / 0.2e1;
t373 = t261 - t493;
t366 = t247 + t394;
t364 = qJD(6) * t377;
t363 = rSges(7,1) * t89 + rSges(7,2) * t88;
t362 = -t492 - t493;
t361 = -t275 - t270;
t220 = rSges(3,1) * t272 + rSges(3,2) * t274;
t358 = rSges(5,1) * t278 + t480;
t120 = -pkin(8) * t392 + t397;
t51 = -rSges(7,3) * t392 + t401;
t13 = qJD(1) * t120 - t121 * t203 - t147 * t168 + t248 * t51 + (t109 * t409 - t209 * t274 + t224 * t417) * qJD(5) + t293;
t119 = t306 * pkin(8) + (t271 * t417 - t273 * t411) * pkin(5);
t50 = rSges(7,3) * t306 + t363;
t14 = t121 * t202 + t148 * t168 - t248 * t50 + (-t111 * t409 + t209 * t272) * qJD(5) + (-t119 + t389 + t434) * qJD(1) + t295;
t350 = t13 * t272 + t14 * t274;
t28 = -t100 * t455 + t483;
t349 = t272 * t29 + t274 * t28;
t348 = t272 * t28 - t274 * t29;
t31 = t102 * t451 - t334;
t347 = t272 * t31 + t274 * t30;
t346 = t272 * t30 - t274 * t31;
t345 = t272 * t41 + t274 * t40;
t344 = t272 * t40 - t274 * t41;
t265 = t272 * rSges(6,3);
t157 = t274 * t357 - t265;
t322 = -t218 + t353;
t312 = t170 + t322;
t58 = t196 + (t157 + t312) * qJD(1) + t422;
t321 = t518 - t352;
t311 = -t274 * t446 + t247 + t321;
t59 = (t156 + t311) * qJD(1) + t539;
t343 = t272 * t58 - t274 * t59;
t332 = t109 * t274 + t111 * t272;
t331 = t152 * t273 + t154 * t271;
t328 = -t156 * t272 - t157 * t274;
t314 = (t272 * t61 + t274 * t60) * qJD(5);
t144 = t272 * t150;
t62 = -t274 * t331 + t144;
t63 = -t151 * t272 + t519;
t313 = (t272 * t63 + t274 * t62) * qJD(5);
t310 = -t476 + t490 + t491;
t308 = t357 + t491;
t303 = t100 * t203 + t102 * t202 + t162 * t248;
t302 = (Icges(7,5) * t192 - Icges(7,6) * t193) * t203 + (Icges(7,5) * t194 + Icges(7,6) * t195) * t202 + t199 * t248;
t300 = -qJD(1) * t329 - qJD(5) * t458 + t420;
t299 = qJD(1) * t331 + qJD(5) * t172 + t419;
t296 = t326 * qJD(1) - qJD(5) * t337;
t292 = -t272 * t99 + t274 * t98 + (-t156 * t274 + t157 * t272) * qJD(1);
t32 = (t188 + t312) * qJD(1) + t422 + t526;
t33 = -t389 + t109 * t248 - t168 * t203 + (t186 + t311) * qJD(1) + t369;
t34 = -t109 * t202 + t111 * t203 + qJD(2) + (-t186 * t272 - t188 * t274) * qJD(5);
t287 = t34 * t332 + (-t272 * t33 - t274 * t32) * t168;
t189 = t224 * t274;
t187 = t224 * t272;
t182 = t221 * t272;
t160 = -rSges(5,3) * t272 + t274 * t358;
t143 = t168 * t274;
t142 = t168 * t272;
t141 = t166 * t274;
t140 = t166 * t272;
t139 = t164 * t274;
t138 = t164 * t272;
t129 = rSges(7,1) * t194 + rSges(7,2) * t195;
t128 = rSges(7,1) * t192 - rSges(7,2) * t193;
t113 = -t402 + t252 + (-qJD(1) * t370 - t181) * qJD(1);
t112 = qJD(1) * t423 - t403 + t430;
t84 = t272 * t326 + t458;
t82 = t84 * qJD(1);
t81 = (t321 + t515) * qJD(1) + t369;
t80 = (t160 + t322) * qJD(1) + t422;
t77 = qJD(5) * t328 + qJD(2);
t67 = (-qJD(1) * t515 - t181) * qJD(1) + t295;
t66 = qJD(1) * (-rSges(5,3) * t417 + t426) + t294;
t39 = -t272 * t510 + t274 * t296;
t38 = t272 * t296 + t274 * t510;
t37 = qJD(5) * t329 - t271 * t94 + t273 * t96;
t36 = -qJD(5) * t331 - t271 * t95 + t273 * t97;
t35 = t292 * qJD(5);
t24 = t313 - t447;
t23 = t82 + t314;
t16 = -t116 * t455 + t117 * t192 + t118 * t193 + t162 * t307 + t164 * t90 + t166 * t91;
t15 = t116 * t451 + t117 * t194 - t118 * t195 + t162 * t306 + t164 * t88 + t166 * t89;
t12 = t202 * t41 + t203 * t40 + t248 * t65;
t11 = t202 * t31 - t535;
t10 = t203 * t28 + t523;
t9 = -t109 * t148 + t111 * t147 - t202 * t51 + t203 * t50 + (t119 * t274 - t120 * t272 + (-t186 * t274 + t188 * t272) * qJD(1)) * qJD(5);
t6 = t102 * t307 + t105 * t90 - t107 * t91 + t192 * t46 + t193 * t48 - t44 * t455;
t5 = t100 * t307 + t103 * t90 + t106 * t91 + t192 * t47 + t193 * t49 - t45 * t455;
t4 = t102 * t306 + t105 * t88 - t107 * t89 + t194 * t46 - t195 * t48 + t44 * t451;
t3 = t100 * t306 + t103 * t88 + t106 * t89 + t194 * t47 - t195 * t49 + t45 * t451;
t2 = t147 * t28 + t148 * t29 + t16 * t248 + t202 * t6 + t203 * t5 + t383 * t53;
t1 = t147 * t30 + t148 * t31 + t15 * t248 + t202 * t4 + t203 * t3 + t383 * t54;
t17 = [((t31 + t538) * t203 + t523) * t504 + (-qJD(5) * t326 + t206 * t271 - t207 * t273) * qJD(1) + (t82 + ((-t62 + t144 + t61) * t272 + (t63 - t519 + (t151 - t331) * t272 + t60) * t274) * qJD(5)) * t379 + m(3) * ((-t220 * t285 - t403) * t517 + (-t402 + (-0.2e1 * t371 - t275 + t517) * t285) * (-t220 - t493)) + t54 * t505 + t53 * t506 + t16 * t501 + t470 / 0.2e1 + t471 / 0.2e1 + t484 + t486 / 0.2e1 + t487 / 0.2e1 + (t15 + t10) * t503 + (t11 + (-t28 + t538) * t202 + t535) * t502 + (t447 + (t272 ^ 2 * t151 + (-t144 + t61 + (t151 + t331) * t274) * t274) * qJD(5) + t24) * t376 + (t14 * (-t237 + t356 + t373) + t13 * (t236 + t366 + t429) + (t14 * (-pkin(2) + t280) - t13 * t494 * t273) * t272 + (-t13 * t280 + t14 * t310) * t274 + (t242 - t363 + (t271 * t494 + t488) * t411 + (t361 + (-qJ(3) - t310 + t489) * t272) * qJD(1) - t369) * t32 + (t397 + t401 + t32 - t526 + (-t188 - t237 - t399 + t536 - t492) * qJD(1) + t537) * t33) * m(7) + (t43 * (t272 * t280 + t261 - t265 + t362) + t42 * (t366 + t156) + (-t42 * t280 + t308 * t43) * t274 + (t242 + (t361 - t267 + (-qJ(3) - t308) * t272) * qJD(1) - t539) * t58 + (-rSges(6,2) * t388 - t196 + t398 + t58 + ((-rSges(6,3) - pkin(2)) * t272 - t157 + t536) * qJD(1) + t537) * t59) * m(6) + (-(-t80 + (t160 + t353) * qJD(1) + t395) * t81 + t67 * (rSges(5,1) * t450 + t373 + t400) + t80 * t259 + t66 * (t394 + t460 + t515) + t81 * (t396 + t426) + (-t80 * qJD(4) + t404 * t67) * t272 + ((-t282 * t81 - t284 * t80) * pkin(1) + t80 * t404 * t274 + (t80 * (-qJ(3) - t358) + t81 * t404) * t272) * qJD(1)) * m(5) + (t113 * (t478 + (rSges(4,2) - pkin(2)) * t272 + t373) + t112 * t530 + (-t421 + t423 + t424 + (-rSges(4,2) * t272 + t362 - t478 + t493) * qJD(1)) * (qJD(1) * t530 - t259)) * m(4) + (t37 + t38 + t23) * t378 + (qJD(1) * t79 + t36 + t39) * t411 / 0.2e1 + (t274 * t85 + (-t330 + t84) * t272) * qJD(5) * t485; m(6) * t35 + m(7) * t9; 0.2e1 * (t13 * t496 + t14 * t497) * m(7) + 0.2e1 * (t473 / 0.2e1 - t472 / 0.2e1) * m(6) + 0.2e1 * (t496 * t66 + t497 * t67) * m(5) + 0.2e1 * (t112 * t496 + t113 * t497) * m(4); m(5) * (t272 * t66 + t274 * t67) + m(6) * (t272 * t42 + t274 * t43) + m(7) * t350; t345 * t364 - t272 * t10 * t410 / 0.2e1 + qJD(1) * (t272 * t37 + t274 * t36 + (t272 * t330 + t274 * t79) * qJD(1)) / 0.2e1 + ((t138 * t194 - t140 * t195) * t203 + (-t139 * t194 + t141 * t195) * t202 + (t163 * t194 - t165 * t195) * t248 + (t273 * t54 + t30 * t457) * qJD(6) + ((-qJD(6) * t31 - t303) * t271 + t522) * t274) * t504 + t347 * t505 + t349 * t506 + t11 * t408 * t498 + (-qJD(1) * t344 + t272 * t8 + t274 * t7) * t499 + (((-t138 * t281 + t140 * t283 + t100) * t203 + (t139 * t281 - t141 * t283 + t102) * t202 + (-t163 * t281 + t165 * t283 + t162) * t248 + t65 * qJD(6)) * t273 + (qJD(6) * t344 + t286) * t271) * t500 + (-qJD(1) * t348 + t272 * t6 + t274 * t5) * t501 + ((t138 * t192 + t140 * t193) * t203 + (-t139 * t192 - t141 * t193) * t202 + (t163 * t192 + t165 * t193) * t248 + (t273 * t53 - t29 * t456) * qJD(6) + ((qJD(6) * t28 + t303) * t271 - t522) * t272) * t502 + (-qJD(1) * t346 + t272 * t4 + t274 * t3) * t503 + ((-t271 * t428 - t273 * t427) * qJD(1) + (t271 * t297 + t273 * t298) * qJD(5)) * t485 - t12 * t409 / 0.2e1 + ((t172 * t411 - t418) * t274 + (-t527 + (-t274 * t458 - t528) * qJD(5)) * t272) * t376 + ((-t413 * t458 - t418) * t272 + (t527 + (t272 * t172 + t528) * qJD(5)) * t274) * t379 + (-t32 * (qJD(1) * t189 + t143 * t248 + t167 * t202 + t223 * t413) - t33 * (qJD(1) * t187 + t142 * t248 - t167 * t203 - t223 * t411) - t34 * (-t142 * t202 - t143 * t203 - t187 * t413 - t189 * t411) - ((t109 * t33 - t111 * t32) * t273 + t287 * t271) * qJD(6) + (-t13 * t431 - t33 * t439 + t9 * t440 + t34 * (t119 + t50) + (t32 * t431 + t34 * t441) * qJD(1)) * t274 + (t14 * t431 + t32 * t439 + t9 * t441 + t34 * (-t120 - t51) + (t33 * t431 - t34 * t440) * qJD(1)) * t272) * m(7) + (-(t182 * t59 + t183 * t58) * qJD(1) - (t77 * (-t182 * t272 - t183 * t274) - t343 * t357) * qJD(5) + t35 * t328 + t77 * t292 - t343 * t208 + (t473 - t472 + (t272 * t59 + t274 * t58) * qJD(1)) * t221) * m(6) + (qJD(1) * t38 + t1 + ((t272 * t299 + t274 * t509) * t274 + (t272 * t300 - t274 * t514) * t272 + (-t272 * t62 + t274 * t63) * qJD(1)) * t524) * t497 + (qJD(1) * t39 + t2 + ((-t272 * t509 + t274 * t299) * t274 + (t272 * t514 + t274 * t300) * t272 + (-t272 * t60 + t274 * t61) * qJD(1)) * t524) * t274 / 0.2e1 - (t314 + t23 + t10) * t417 / 0.2e1 + (t313 + t24 + t11) * t415 / 0.2e1; -t2 * t455 / 0.2e1 + (t271 * t53 - t273 * t348) * t506 + ((qJD(5) * t348 + t16) * t271 + (-qJD(1) * t349 + qJD(5) * t53 - t272 * t5 + t274 * t6) * t273) * t501 + t1 * t451 / 0.2e1 + (t271 * t54 - t273 * t346) * t505 + ((qJD(5) * t346 + t15) * t271 + (-qJD(1) * t347 + qJD(5) * t54 - t272 * t3 + t274 * t4) * t273) * t503 + t12 * t377 + (t470 + t471 + t484 + t486 + t487) * t498 + (t271 * t65 - t273 * t344) * t364 + ((qJD(5) * t344 + t22) * t271 + (-qJD(1) * t345 + qJD(5) * t65 - t272 * t7 + t274 * t8) * t273) * t499 + (t192 * t290 - t193 * t507 - t302 * t455) * t502 + (t290 * t194 + t195 * t507 + t302 * t451) * t504 + (t302 * t271 + (-t281 * t290 - t283 * t507) * t273) * t500 + (t271 * t376 + t272 * t381) * t11 + (t271 * t378 + t274 * t381) * t10 + ((qJD(5) * t287 + t109 * t13 - t111 * t14 - t32 * t50 + t33 * t51) * t271 + (t32 * (-qJD(5) * t111 + t121 * t274) + t33 * (qJD(5) * t109 + t121 * t272) - t9 * t332 + t34 * (t109 * t417 - t111 * t415 - t272 * t50 - t274 * t51) + ((-t272 * t32 + t274 * t33) * qJD(1) + t350) * t168) * t273 - t32 * (-t129 * t248 + t202 * t204) - t33 * (t128 * t248 - t203 * t204) - t34 * (-t128 * t202 + t129 * t203)) * m(7);];
tauc  = t17(:);
