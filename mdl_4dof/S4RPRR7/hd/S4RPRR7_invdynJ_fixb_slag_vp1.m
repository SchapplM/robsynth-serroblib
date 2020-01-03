% Calculate vector of inverse dynamics joint torques for
% S4RPRR7
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR7_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR7_invdynJ_fixb_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:42
% EndTime: 2019-12-31 16:54:03
% DurationCPUTime: 17.78s
% Computational Cost: add. (13093->747), mult. (18717->1016), div. (0->0), fcn. (17508->8), ass. (0->356)
t289 = cos(qJ(1));
t287 = sin(qJ(1));
t282 = pkin(7) + qJ(3);
t267 = cos(t282);
t255 = t267 * pkin(3);
t266 = sin(t282);
t490 = t266 * pkin(6) + t255;
t189 = t490 * t287;
t273 = t289 * qJ(2);
t227 = pkin(1) * t287 - t273;
t284 = cos(pkin(7));
t256 = pkin(2) * t284 + pkin(1);
t285 = -pkin(5) - qJ(2);
t261 = t289 * t285;
t397 = -t287 * t256 - t261;
t156 = t227 + t397;
t408 = t156 - t227;
t372 = -t189 + t408;
t288 = cos(qJ(4));
t421 = t288 * t289;
t286 = sin(qJ(4));
t423 = t286 * t287;
t193 = t267 * t423 + t421;
t419 = t289 * t286;
t422 = t287 * t288;
t194 = t267 * t422 - t419;
t349 = t194 * rSges(5,1) - t193 * rSges(5,2);
t427 = t266 * t287;
t109 = -rSges(5,3) * t427 - t349;
t455 = rSges(5,2) * t286;
t457 = rSges(5,1) * t288;
t348 = -t455 + t457;
t154 = -rSges(5,3) * t267 + t266 * t348;
t385 = qJD(4) * t266;
t386 = qJD(3) * t289;
t204 = -t287 * t385 + t386;
t464 = t266 * pkin(3);
t213 = -pkin(6) * t267 + t464;
t384 = qJD(4) * t267;
t243 = qJD(1) - t384;
t270 = qJD(2) * t287;
t501 = -t109 * t243 + t154 * t204 + t213 * t386 - t270;
t35 = qJD(1) * t372 - t501;
t502 = t289 * t35;
t178 = Icges(5,4) * t194;
t102 = -Icges(5,2) * t193 + Icges(5,6) * t427 + t178;
t177 = Icges(5,4) * t193;
t106 = -Icges(5,1) * t194 - Icges(5,5) * t427 + t177;
t497 = t102 * t286 + t106 * t288;
t99 = Icges(5,5) * t194 - Icges(5,6) * t193 + Icges(5,3) * t427;
t40 = -t266 * t497 - t267 * t99;
t253 = t266 * rSges(5,3);
t155 = t267 * t348 + t253;
t387 = qJD(3) * t287;
t203 = t289 * t385 + t387;
t27 = -t102 * t193 - t106 * t194 + t427 * t99;
t195 = -t267 * t419 + t422;
t196 = t267 * t421 + t423;
t426 = t266 * t289;
t101 = Icges(5,5) * t196 + Icges(5,6) * t195 + Icges(5,3) * t426;
t441 = Icges(5,4) * t196;
t104 = Icges(5,2) * t195 + Icges(5,6) * t426 + t441;
t179 = Icges(5,4) * t195;
t107 = Icges(5,1) * t196 + Icges(5,5) * t426 + t179;
t28 = t101 * t427 - t193 * t104 + t194 * t107;
t331 = Icges(5,5) * t288 - Icges(5,6) * t286;
t148 = -Icges(5,3) * t267 + t266 * t331;
t439 = Icges(5,4) * t288;
t332 = -Icges(5,2) * t286 + t439;
t150 = -Icges(5,6) * t267 + t266 * t332;
t440 = Icges(5,4) * t286;
t334 = Icges(5,1) * t288 - t440;
t152 = -Icges(5,5) * t267 + t266 * t334;
t51 = t148 * t427 - t150 * t193 + t152 * t194;
t12 = t203 * t28 - t204 * t27 + t243 * t51;
t29 = t195 * t102 - t106 * t196 + t99 * t426;
t30 = t101 * t426 + t195 * t104 + t196 * t107;
t52 = t148 * t426 + t150 * t195 + t152 * t196;
t13 = t203 * t30 - t204 * t29 + t52 * t243;
t222 = qJD(1) * t227;
t491 = qJD(1) * t156 - t222;
t252 = Icges(4,4) * t267;
t333 = -Icges(4,2) * t266 + t252;
t209 = Icges(4,1) * t266 + t252;
t272 = t287 * qJ(2);
t229 = t289 * pkin(1) + t272;
t381 = qJD(1) * qJD(3);
t220 = -qJDD(3) * t289 + t287 * t381;
t368 = t267 * t387;
t389 = qJD(1) * t289;
t309 = t266 * t389 + t368;
t380 = qJDD(4) * t266;
t112 = qJD(4) * t309 + t287 * t380 + t220;
t370 = t266 * t387;
t116 = t309 * pkin(6) + (t267 * t389 - t370) * pkin(3);
t192 = qJD(3) * t385 - qJDD(4) * t267 + qJDD(1);
t202 = t490 * qJD(3);
t382 = qJD(1) * qJD(2);
t393 = qJDD(2) * t287 + t289 * t382;
t271 = qJD(2) * t289;
t200 = qJD(1) * t229 - t271;
t390 = qJD(1) * t287;
t249 = t285 * t390;
t461 = pkin(1) - t256;
t411 = t249 - (-t289 * t461 - t272) * qJD(1) - t200;
t354 = qJD(1) * t267 - qJD(4);
t88 = t243 * t422 + (-t289 * t354 + t370) * t286;
t388 = qJD(3) * t266;
t89 = t354 * t421 + (t243 * t286 - t288 * t388) * t287;
t353 = rSges(5,1) * t89 + rSges(5,2) * t88;
t50 = rSges(5,3) * t309 + t353;
t180 = (-rSges(5,1) * t286 - rSges(5,2) * t288) * t266;
t90 = qJD(3) * t155 + qJD(4) * t180;
t10 = -t202 * t386 + t109 * t192 + t112 * t154 - t204 * t90 + t213 * t220 - t243 * t50 + t372 * qJDD(1) + (-t116 + t411) * qJD(1) + t393;
t110 = t196 * rSges(5,1) + t195 * rSges(5,2) + rSges(5,3) * t426;
t219 = qJDD(3) * t287 + t289 * t381;
t367 = t267 * t386;
t371 = t266 * t390;
t308 = t367 - t371;
t111 = qJD(4) * t308 + t289 * t380 + t219;
t225 = pkin(6) * t367;
t369 = t266 * t386;
t310 = -t267 * t390 - t369;
t115 = pkin(3) * t310 - pkin(6) * t371 + t225;
t424 = t267 * t289;
t191 = pkin(3) * t424 + pkin(6) * t426;
t355 = t289 * t256 - t285 * t287;
t157 = t355 - t229;
t262 = qJ(2) * t389;
t392 = t262 + t270;
t320 = -qJDD(2) * t289 + qJD(1) * (-pkin(1) * t390 + t392) + qJDD(1) * t229 + t287 * t382;
t305 = qJD(1) * (-t262 + (t287 * t461 - t261) * qJD(1)) + qJDD(1) * t157 + t320;
t321 = t243 * t289;
t484 = t287 * t354 + t369;
t86 = t286 * t484 + t288 * t321;
t87 = t286 * t321 - t288 * t484;
t379 = t87 * rSges(5,1) + t86 * rSges(5,2) + rSges(5,3) * t367;
t49 = -rSges(5,3) * t371 + t379;
t11 = qJD(1) * t115 + qJDD(1) * t191 + t110 * t192 - t111 * t154 - t202 * t387 - t203 * t90 - t213 * t219 + t243 * t49 + t305;
t488 = t10 * t287 - t11 * t289;
t283 = sin(pkin(7));
t456 = rSges(3,2) * t283;
t459 = rSges(3,1) * t284;
t167 = t287 * rSges(3,3) + (-t456 + t459) * t289;
t465 = g(1) * t287;
t487 = -g(2) * t289 + t465;
t206 = Icges(4,5) * t267 - Icges(4,6) * t266;
t205 = Icges(4,5) * t266 + Icges(4,6) * t267;
t311 = qJD(3) * t205;
t442 = Icges(4,4) * t266;
t210 = Icges(4,1) * t267 - t442;
t163 = Icges(4,5) * t287 + t210 * t289;
t161 = Icges(4,6) * t287 + t289 * t333;
t430 = t161 * t266;
t325 = -t163 * t267 + t430;
t433 = Icges(4,3) * t289;
t486 = -t289 * t311 + (-t206 * t287 + t325 + t433) * qJD(1);
t233 = Icges(4,4) * t427;
t425 = t267 * t287;
t438 = Icges(4,5) * t289;
t162 = Icges(4,1) * t425 - t233 - t438;
t435 = Icges(4,6) * t289;
t160 = Icges(4,4) * t425 - Icges(4,2) * t427 - t435;
t431 = t160 * t266;
t326 = -t162 * t267 + t431;
t159 = Icges(4,3) * t287 + t206 * t289;
t391 = qJD(1) * t159;
t485 = qJD(1) * t326 - t287 * t311 + t391;
t158 = Icges(4,5) * t425 - Icges(4,6) * t427 - t433;
t57 = -t289 * t158 - t287 * t326;
t207 = Icges(4,2) * t267 + t442;
t322 = t207 * t266 - t209 * t267;
t483 = t322 * qJD(1) + t206 * qJD(3);
t482 = t287 * (-t207 * t289 + t163) - t289 * (-Icges(4,2) * t425 + t162 - t233);
t149 = Icges(5,3) * t266 + t267 * t331;
t327 = -t150 * t286 + t152 * t288;
t329 = -t104 * t286 + t107 * t288;
t481 = t203 * (-t148 * t289 - t329) - t204 * (-t148 * t287 + t497) + t243 * (t149 - t327);
t171 = (-Icges(5,2) * t288 - t440) * t266;
t480 = t203 * (-Icges(5,2) * t196 + t107 + t179) - t204 * (-Icges(5,2) * t194 - t106 - t177) + t243 * (t152 + t171);
t479 = t111 / 0.2e1;
t478 = t112 / 0.2e1;
t477 = t192 / 0.2e1;
t476 = -t203 / 0.2e1;
t475 = t203 / 0.2e1;
t474 = -t204 / 0.2e1;
t473 = t204 / 0.2e1;
t472 = t219 / 0.2e1;
t471 = t220 / 0.2e1;
t470 = -t243 / 0.2e1;
t469 = t243 / 0.2e1;
t468 = t287 / 0.2e1;
t467 = -t289 / 0.2e1;
t466 = -rSges(5,3) - pkin(6);
t44 = Icges(5,5) * t89 + Icges(5,6) * t88 + Icges(5,3) * t309;
t46 = Icges(5,4) * t89 + Icges(5,2) * t88 + Icges(5,6) * t309;
t48 = Icges(5,1) * t89 + Icges(5,4) * t88 + Icges(5,5) * t309;
t7 = (-qJD(3) * t497 - t44) * t267 + (qJD(3) * t99 - t286 * t46 + t288 * t48 + (-t102 * t288 + t106 * t286) * qJD(4)) * t266;
t463 = t7 * t204;
t43 = Icges(5,5) * t87 + Icges(5,6) * t86 + Icges(5,3) * t308;
t45 = Icges(5,4) * t87 + Icges(5,2) * t86 + Icges(5,6) * t308;
t47 = Icges(5,1) * t87 + Icges(5,4) * t86 + Icges(5,5) * t308;
t8 = (qJD(3) * t329 - t43) * t267 + (qJD(3) * t101 - t286 * t45 + t288 * t47 + (-t104 * t288 - t107 * t286) * qJD(4)) * t266;
t462 = t8 * t203;
t168 = (-Icges(5,5) * t286 - Icges(5,6) * t288) * t266;
t83 = qJD(3) * t149 + qJD(4) * t168;
t151 = Icges(5,6) * t266 + t267 * t332;
t84 = qJD(3) * t151 + qJD(4) * t171;
t153 = Icges(5,5) * t266 + t267 * t334;
t174 = (-Icges(5,1) * t286 - t439) * t266;
t85 = qJD(3) * t153 + qJD(4) * t174;
t18 = (qJD(3) * t327 - t83) * t267 + (qJD(3) * t148 - t286 * t84 + t288 * t85 + (-t150 * t288 - t152 * t286) * qJD(4)) * t266;
t54 = -t148 * t267 + t266 * t327;
t460 = t18 * t243 + t54 * t192;
t458 = rSges(4,1) * t267;
t211 = rSges(4,1) * t266 + rSges(4,2) * t267;
t182 = t211 * t289;
t277 = t287 * rSges(4,3);
t165 = rSges(4,1) * t424 - rSges(4,2) * t426 + t277;
t407 = t157 + t229;
t64 = -t211 * t387 - t271 + (t165 + t407) * qJD(1);
t452 = t182 * t64;
t36 = -t213 * t387 + t110 * t243 - t154 * t203 - t271 + (t191 + t407) * qJD(1);
t448 = t287 * t36;
t351 = -t211 * t386 + t270;
t164 = rSges(4,1) * t425 - rSges(4,2) * t427 - t289 * rSges(4,3);
t373 = -t164 + t408;
t63 = qJD(1) * t373 + t351;
t447 = t287 * t63;
t446 = t40 * t112;
t41 = -t101 * t267 + t266 * t329;
t445 = t41 * t111;
t444 = -t202 - t90;
t429 = t205 * t287;
t428 = t205 * t289;
t73 = -t287 * t322 - t428;
t418 = t73 * qJD(1);
t415 = -t109 + t189;
t414 = t110 + t191;
t413 = -t287 * t158 - t162 * t424;
t412 = t287 * t159 + t163 * t424;
t409 = t154 + t213;
t143 = t167 + t229;
t403 = -t207 + t210;
t402 = t209 + t333;
t375 = t266 * t455;
t401 = rSges(5,3) * t425 + t287 * t375;
t400 = rSges(5,3) * t424 + t289 * t375;
t399 = rSges(4,2) * t371 + rSges(4,3) * t389;
t378 = t287 * t459;
t250 = t287 * t456;
t394 = t289 * rSges(3,3) + t250;
t166 = t378 - t394;
t398 = -t227 - t166;
t396 = rSges(3,3) * t389 + qJD(1) * t250;
t395 = t249 + t271;
t383 = t206 * qJD(1);
t377 = t266 * t457;
t366 = -pkin(1) - t459;
t364 = t390 / 0.2e1;
t363 = t389 / 0.2e1;
t362 = -t387 / 0.2e1;
t361 = t387 / 0.2e1;
t360 = -t386 / 0.2e1;
t359 = t386 / 0.2e1;
t358 = -t384 / 0.2e1;
t138 = t163 * t425;
t357 = t289 * t159 - t138;
t356 = -t158 + t430;
t230 = rSges(2,1) * t289 - rSges(2,2) * t287;
t228 = rSges(2,1) * t287 + rSges(2,2) * t289;
t212 = -rSges(4,2) * t266 + t458;
t76 = t161 * t267 + t163 * t266;
t312 = qJD(3) * t207;
t93 = -t289 * t312 + (-t287 * t333 + t435) * qJD(1);
t313 = qJD(3) * t209;
t95 = -t289 * t313 + (-t210 * t287 + t438) * qJD(1);
t294 = -qJD(3) * t76 - t266 * t93 + t267 * t95 + t391;
t75 = t160 * t267 + t162 * t266;
t94 = qJD(1) * t161 - t287 * t312;
t96 = qJD(1) * t163 - t287 * t313;
t295 = qJD(1) * t158 - qJD(3) * t75 - t266 * t94 + t267 * t96;
t347 = -(t287 * t485 + t295 * t289) * t289 + (t287 * t486 + t294 * t289) * t287;
t346 = -(t295 * t287 - t289 * t485) * t289 + (t294 * t287 - t289 * t486) * t287;
t345 = t27 * t289 - t28 * t287;
t344 = t27 * t287 + t28 * t289;
t343 = t287 * t30 - t289 * t29;
t342 = t287 * t29 + t289 * t30;
t341 = t287 * t41 - t289 * t40;
t340 = t287 * t40 + t289 * t41;
t58 = -t161 * t427 - t357;
t339 = t287 * t58 - t289 * t57;
t59 = -t160 * t426 - t413;
t60 = -t161 * t426 + t412;
t338 = t287 * t60 - t289 * t59;
t337 = -t287 * t64 - t289 * t63;
t97 = rSges(4,1) * t310 - rSges(4,2) * t367 + t399;
t181 = t211 * t287;
t98 = -qJD(3) * t181 + (t212 * t289 + t277) * qJD(1);
t336 = t287 * t98 + t289 * t97;
t328 = -t109 * t289 - t110 * t287;
t324 = t164 * t287 + t165 * t289;
t323 = t207 * t267 + t209 * t266;
t319 = t266 * t466 - t255;
t307 = t101 * t203 + t148 * t243 - t204 * t99;
t306 = (-Icges(5,5) * t193 - Icges(5,6) * t194) * t204 - (Icges(5,5) * t195 - Icges(5,6) * t196) * t203 - t168 * t243;
t304 = t160 * t289 - t161 * t287;
t303 = t266 * t306;
t299 = (-t266 * t402 + t267 * t403) * qJD(1);
t298 = (Icges(5,1) * t195 - t104 - t441) * t203 - (-Icges(5,1) * t193 - t102 - t178) * t204 + (-t150 + t174) * t243;
t37 = -t109 * t203 + t110 * t204 + (t189 * t287 + t191 * t289) * qJD(3);
t296 = t37 * t328 + (t287 * t35 - t289 * t36) * t154;
t198 = t333 * qJD(3);
t199 = t210 * qJD(3);
t293 = qJD(1) * t205 - qJD(3) * t323 - t198 * t266 + t199 * t267;
t292 = -t266 * t482 + t304 * t267;
t291 = t481 * t266;
t246 = pkin(6) * t424;
t244 = pkin(6) * t425;
t201 = t212 * qJD(3);
t190 = -pkin(3) * t426 + t246;
t188 = -pkin(3) * t427 + t244;
t136 = -t289 * t377 + t400;
t135 = -t287 * t377 + t401;
t134 = t152 * t289;
t133 = t152 * t287;
t132 = t150 * t289;
t131 = t150 * t287;
t128 = qJD(1) * t143 - t271;
t127 = qJD(1) * t398 + t270;
t126 = rSges(5,1) * t195 - rSges(5,2) * t196;
t125 = -rSges(5,1) * t193 - rSges(5,2) * t194;
t77 = t324 * qJD(3);
t74 = -t289 * t322 + t429;
t72 = t74 * qJD(1);
t62 = qJDD(1) * t167 + qJD(1) * (-qJD(1) * t378 + t396) + t320;
t61 = t398 * qJDD(1) + (-qJD(1) * t167 - t200) * qJD(1) + t393;
t39 = -qJD(3) * t325 + t266 * t95 + t267 * t93;
t38 = -t326 * qJD(3) + t266 * t96 + t267 * t94;
t34 = t293 * t287 - t289 * t483;
t33 = t287 * t483 + t293 * t289;
t26 = qJD(1) * t97 + qJDD(1) * t165 - t201 * t387 - t211 * t219 + t305;
t25 = -t201 * t386 + t211 * t220 + t373 * qJDD(1) + (-t98 + t411) * qJD(1) + t393;
t24 = qJD(3) * t338 + t72;
t23 = qJD(3) * t339 + t418;
t16 = t148 * t309 + t150 * t88 + t152 * t89 - t193 * t84 + t194 * t85 + t427 * t83;
t15 = t148 * t308 + t150 * t86 + t152 * t87 + t195 * t84 + t196 * t85 + t426 * t83;
t14 = t203 * t41 - t204 * t40 + t243 * t54;
t9 = -t109 * t111 - t110 * t112 + t189 * t219 - t191 * t220 + t203 * t50 + t204 * t49 + (t115 * t289 + t116 * t287) * qJD(3);
t6 = t101 * t309 + t104 * t88 + t107 * t89 - t193 * t45 + t194 * t47 + t427 * t43;
t5 = t99 * t368 + t102 * t88 - t106 * t89 - t193 * t46 + t194 * t48 + (t287 * t44 + t389 * t99) * t266;
t4 = t101 * t308 + t104 * t86 + t107 * t87 + t195 * t45 + t196 * t47 + t426 * t43;
t3 = t99 * t367 + t102 * t86 - t106 * t87 + t195 * t46 + t196 * t48 + (t289 * t44 - t390 * t99) * t266;
t2 = t111 * t28 + t112 * t27 + t16 * t243 + t192 * t51 + t203 * t6 - t204 * t5;
t1 = t111 * t30 + t112 * t29 + t15 * t243 + t192 * t52 + t203 * t4 - t204 * t3;
t17 = [t445 / 0.2e1 + t446 / 0.2e1 - t463 / 0.2e1 + t462 / 0.2e1 + (-qJD(3) * t322 + t198 * t267 + t199 * t266) * qJD(1) + (t72 + ((t58 - t138 + (t159 + t431) * t289 + t413) * t289 + t412 * t287) * qJD(3)) * t359 + t460 + t13 * t473 + t15 * t475 + t51 * t478 + t52 * t479 - m(2) * (-g(1) * t228 + g(2) * t230) + (t16 + t13) * t474 + (t76 + t74) * t472 + (t75 + t73) * t471 + (-t418 + ((t289 * t356 - t412 + t60) * t289 + (t287 * t356 + t357 + t59) * t287) * qJD(3) + t23) * t362 + (t39 + t33) * t361 + (t38 + t34 + t24) * t360 + (t35 * (-t353 + t395) + (t10 * t319 + t35 * (t267 * t466 + t464) * qJD(3)) * t287 + ((-t256 + t319) * t448 + (-t256 - t490 - t253) * t502) * qJD(1) - t319 * t465 + (t11 - g(2)) * (t355 + t414) + (t10 - g(1)) * (-t349 + t397) + (-pkin(3) * t369 + t225 + t270 + t35 + t379 - t491 + (t189 - t261) * qJD(1) + t501) * t36) * m(5) + (-(-qJD(1) * t164 + t351 + t491 - t63) * t64 + t63 * t395 + t64 * (t270 + t399) + (t211 * t447 - t452) * qJD(3) + ((-t63 * rSges(4,3) + t64 * (-t256 - t458)) * t287 + (t63 * (-t212 - t256) - t64 * t285) * t289) * qJD(1) + (t26 - g(2)) * (t165 + t355) + (t25 - g(1)) * (-t164 + t397)) * m(4) + (t127 * t271 + t128 * (t392 + t396) + (t127 * (t366 + t456) * t289 + (t127 * (-rSges(3,3) - qJ(2)) + t128 * t366) * t287) * qJD(1) - (-qJD(1) * t166 - t127 - t222 + t270) * t128 + (t62 - g(2)) * t143 + (t61 - g(1)) * (t366 * t287 + t273 + t394)) * m(3) + (Icges(3,2) * t284 ^ 2 + (Icges(3,1) * t283 + 0.2e1 * Icges(3,4) * t284) * t283 + m(2) * (t228 ^ 2 + t230 ^ 2) + t323 + Icges(2,3)) * qJDD(1); (-m(3) - m(4)) * t487 + 0.2e1 * (t25 * t468 + t26 * t467) * m(4) + 0.2e1 * (t467 * t62 + t468 * t61) * m(3) + (-t487 + t488) * m(5); (qJD(1) * t33 + qJD(3) * t347 + qJDD(1) * t74 + t219 * t60 + t220 * t59 + t1) * t468 + ((-t10 * t409 + t35 * t444 + t9 * t414 + t37 * (t115 + t49) + (-t36 * t409 + t37 * t415) * qJD(1)) * t289 + (-t11 * t409 + t36 * t444 + t9 * t415 + t37 * (t116 + t50) + (t35 * t409 - t37 * t414) * qJD(1)) * t287 - g(1) * (t246 + t400) - g(2) * (t244 + t401) - g(3) * (t155 + t490) - (g(1) * t289 + g(2) * t287) * t266 * (-pkin(3) - t457) - t35 * (-qJD(1) * t188 - t135 * t243 - t155 * t204 - t386 * t490) - t36 * (qJD(1) * t190 + t136 * t243 - t155 * t203 - t387 * t490) - t37 * (t135 * t203 + t136 * t204 + t188 * t387 + t190 * t386) - ((t109 * t35 + t110 * t36) * t266 + t296 * t267) * qJD(4)) * m(5) - t14 * t385 / 0.2e1 + (((t132 * t286 - t134 * t288 + t101) * t203 - (t131 * t286 - t133 * t288 + t99) * t204 + (-t151 * t286 + t153 * t288 + t148) * t243 + t54 * qJD(4)) * t266 + (qJD(4) * t340 - t481) * t267) * t470 + ((t132 * t193 - t134 * t194) * t203 - (t131 * t193 - t133 * t194) * t204 + (-t151 * t193 + t153 * t194) * t243 + (t266 * t51 + t28 * t424) * qJD(4) + ((qJD(4) * t27 + t307) * t267 + t291) * t287) * t473 + ((-t132 * t195 - t134 * t196) * t203 - (-t131 * t195 - t133 * t196) * t204 + (t151 * t195 + t153 * t196) * t243 + (t266 * t52 + t29 * t425) * qJD(4) + ((qJD(4) * t30 + t307) * t267 + t291) * t289) * t476 + ((t57 * t287 + t58 * t289) * qJD(1) + t346) * t360 + ((t59 * t287 + t60 * t289) * qJD(1) + t347) * t361 - qJD(1) * ((t266 * t403 + t267 * t402) * qJD(1) + (t304 * t266 + t267 * t482) * qJD(3)) / 0.2e1 + qJDD(1) * (t287 * t76 - t289 * t75) / 0.2e1 + (qJD(1) * t34 + qJD(3) * t346 + qJDD(1) * t73 + t219 * t58 + t220 * t57 + t2) * t467 - t112 * t345 / 0.2e1 + ((qJD(3) * t336 + t164 * t219 - t165 * t220) * t324 + t77 * ((t164 * t289 - t165 * t287) * qJD(1) + t336) + t337 * t201 + (-t25 * t289 - t26 * t287 + (-t289 * t64 + t447) * qJD(1)) * t211 + g(1) * t182 + g(2) * t181 - g(3) * t212 - (t181 * t63 - t452) * qJD(1) - (t77 * (-t181 * t287 - t182 * t289) + t337 * t212) * qJD(3)) * m(4) + qJD(1) * (t287 * t39 - t289 * t38 + (t75 * t287 + t289 * t76) * qJD(1)) / 0.2e1 + ((-t386 * t429 - t383) * t289 + (t299 + (t289 * t428 + t292) * qJD(3)) * t287) * t359 + ((-t387 * t428 + t383) * t287 + (t299 + (t287 * t429 + t292) * qJD(3)) * t289) * t362 + t24 * t363 + t23 * t364 + (t287 * t358 + t364) * t12 + (t289 * t358 + t363) * t13 + (qJD(1) * t340 + t287 * t8 - t289 * t7) * t469 + t339 * t471 + t338 * t472 + (qJD(1) * t344 + t287 * t6 - t289 * t5) * t474 + (qJD(1) * t342 + t287 * t4 - t289 * t3) * t475 + t341 * t477 + t343 * t479; t1 * t426 / 0.2e1 + (t266 * t342 - t267 * t52) * t479 + ((qJD(3) * t342 - t15) * t267 + (-qJD(1) * t343 + qJD(3) * t52 + t287 * t3 + t289 * t4) * t266) * t475 + t2 * t427 / 0.2e1 + (t266 * t344 - t267 * t51) * t478 + ((qJD(3) * t344 - t16) * t267 + (qJD(1) * t345 + qJD(3) * t51 + t287 * t5 + t289 * t6) * t266) * t474 + t14 * t388 / 0.2e1 - t267 * (t445 + t446 + t460 + t462 - t463) / 0.2e1 + (t266 * t340 - t267 * t54) * t477 + ((qJD(3) * t340 - t18) * t267 + (-qJD(1) * t341 + qJD(3) * t54 + t287 * t7 + t289 * t8) * t266) * t469 + (t195 * t480 + t298 * t196 - t289 * t303) * t476 + (-t193 * t480 + t194 * t298 - t287 * t303) * t473 + (t306 * t267 + (-t286 * t480 + t288 * t298) * t266) * t470 + (-t371 / 0.2e1 + t267 * t359) * t13 + (t266 * t363 + t267 * t361) * t12 + ((qJD(3) * t296 - t10 * t109 - t11 * t110 + t35 * t50 - t36 * t49) * t267 + (t35 * (qJD(3) * t109 + t287 * t90) + t36 * (qJD(3) * t110 - t289 * t90) + t9 * t328 + t37 * (t109 * t390 - t110 * t389 - t287 * t49 + t289 * t50) + ((t448 + t502) * qJD(1) + t488) * t154) * t266 - t35 * (-t125 * t243 - t180 * t204) - t36 * (t126 * t243 - t180 * t203) - t37 * (t125 * t203 + t126 * t204) - g(1) * t126 - g(2) * t125 - g(3) * t180) * m(5);];
tau = t17;
