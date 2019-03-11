% Calculate time derivative of joint inertia matrix for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
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
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP7_inertiaDJ_slag_vp11(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRP7_inertiaDJ_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:12:42
% EndTime: 2019-03-09 02:12:58
% DurationCPUTime: 12.15s
% Computational Cost: add. (19653->758), mult. (28271->1057), div. (0->0), fcn. (26618->8), ass. (0->365)
t260 = sin(qJ(5));
t253 = pkin(9) + qJ(4);
t244 = sin(t253);
t245 = cos(t253);
t262 = cos(qJ(5));
t418 = Icges(7,4) * t262;
t301 = -Icges(7,2) * t260 + t418;
t173 = Icges(7,6) * t244 + t245 * t301;
t420 = Icges(6,4) * t262;
t302 = -Icges(6,2) * t260 + t420;
t174 = Icges(6,6) * t244 + t245 * t302;
t480 = -t174 - t173;
t485 = t260 * t480;
t297 = Icges(7,5) * t262 - Icges(7,6) * t260;
t171 = Icges(7,3) * t244 + t245 * t297;
t298 = Icges(6,5) * t262 - Icges(6,6) * t260;
t172 = Icges(6,3) * t244 + t245 * t298;
t484 = t172 + t171;
t419 = Icges(7,4) * t260;
t305 = Icges(7,1) * t262 - t419;
t175 = Icges(7,5) * t244 + t245 * t305;
t421 = Icges(6,4) * t260;
t306 = Icges(6,1) * t262 - t421;
t176 = Icges(6,5) * t244 + t245 * t306;
t479 = -t176 - t175;
t258 = -qJ(6) - pkin(8);
t424 = rSges(7,3) - t258;
t261 = sin(qJ(1));
t263 = cos(qJ(1));
t390 = t262 * t263;
t392 = t261 * t260;
t199 = -t244 * t392 + t390;
t391 = t261 * t262;
t394 = t260 * t263;
t200 = t244 * t391 + t394;
t399 = t245 * t261;
t142 = Icges(6,5) * t200 + Icges(6,6) * t199 - Icges(6,3) * t399;
t146 = Icges(6,4) * t200 + Icges(6,2) * t199 - Icges(6,6) * t399;
t150 = Icges(6,1) * t200 + Icges(6,4) * t199 - Icges(6,5) * t399;
t53 = -t142 * t399 + t146 * t199 + t150 * t200;
t201 = t244 * t394 + t391;
t202 = -t244 * t390 + t392;
t397 = t245 * t263;
t143 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t397;
t147 = Icges(6,4) * t202 + Icges(6,2) * t201 + Icges(6,6) * t397;
t151 = Icges(6,1) * t202 + Icges(6,4) * t201 + Icges(6,5) * t397;
t54 = -t143 * t399 + t147 * t199 + t151 * t200;
t316 = t261 * t53 - t263 * t54;
t140 = Icges(7,5) * t200 + Icges(7,6) * t199 - Icges(7,3) * t399;
t144 = Icges(7,4) * t200 + Icges(7,2) * t199 - Icges(7,6) * t399;
t148 = Icges(7,1) * t200 + Icges(7,4) * t199 - Icges(7,5) * t399;
t51 = -t140 * t399 + t144 * t199 + t148 * t200;
t141 = Icges(7,5) * t202 + Icges(7,6) * t201 + Icges(7,3) * t397;
t145 = Icges(7,4) * t202 + Icges(7,2) * t201 + Icges(7,6) * t397;
t149 = Icges(7,1) * t202 + Icges(7,4) * t201 + Icges(7,5) * t397;
t52 = -t141 * t399 + t145 * t199 + t149 * t200;
t317 = t261 * t51 - t263 * t52;
t84 = -t171 * t399 + t173 * t199 + t175 * t200;
t85 = -t172 * t399 + t174 * t199 + t176 * t200;
t483 = (-t316 - t317) * t245 + (t84 + t85) * t244;
t57 = t142 * t397 + t201 * t146 + t202 * t150;
t58 = t143 * t397 + t201 * t147 + t202 * t151;
t314 = t261 * t57 - t263 * t58;
t55 = t140 * t397 + t201 * t144 + t202 * t148;
t56 = t141 * t397 + t201 * t145 + t202 * t149;
t315 = t261 * t55 - t263 * t56;
t86 = t171 * t397 + t201 * t173 + t202 * t175;
t87 = t172 * t397 + t201 * t174 + t202 * t176;
t482 = (-t314 - t315) * t245 + (t86 + t87) * t244;
t481 = (-t262 * t479 + t485) * t245 + t484 * t244;
t371 = qJD(5) * t245;
t118 = (-Icges(7,5) * t260 - Icges(7,6) * t262) * t371 + (Icges(7,3) * t245 - t244 * t297) * qJD(4);
t119 = (-Icges(6,5) * t260 - Icges(6,6) * t262) * t371 + (Icges(6,3) * t245 - t244 * t298) * qJD(4);
t122 = (-Icges(7,1) * t260 - t418) * t371 + (Icges(7,5) * t245 - t244 * t305) * qJD(4);
t123 = (-Icges(6,1) * t260 - t420) * t371 + (Icges(6,5) * t245 - t244 * t306) * qJD(4);
t373 = qJD(4) * t262;
t375 = qJD(4) * t245;
t376 = qJD(4) * t244;
t478 = (t122 + t123) * t245 * t262 + t484 * t375 - t376 * t485 + (t373 * t479 + t118 + t119) * t244;
t374 = qJD(4) * t261;
t344 = t245 * t374;
t377 = qJD(1) * t263;
t477 = t244 * t377 + t344;
t372 = qJD(4) * t263;
t345 = t244 * t372;
t378 = qJD(1) * t261;
t476 = t245 * t378 + t345;
t347 = t244 * t374;
t350 = t245 * t377;
t271 = t347 - t350;
t120 = (-Icges(7,2) * t262 - t419) * t371 + (Icges(7,6) * t245 - t244 * t301) * qJD(4);
t121 = (-Icges(6,2) * t262 - t421) * t371 + (Icges(6,6) * t245 - t244 * t302) * qJD(4);
t475 = (-t120 - t121) * t260;
t368 = qJD(6) * t245;
t332 = qJD(5) * t244 + qJD(1);
t465 = t332 * t260;
t474 = pkin(5) * t465 + t368;
t291 = t145 * t260 - t149 * t262;
t62 = t141 * t244 - t245 * t291;
t289 = t147 * t260 - t151 * t262;
t64 = t143 * t244 - t245 * t289;
t431 = t62 + t64;
t292 = t144 * t260 - t148 * t262;
t61 = t140 * t244 - t245 * t292;
t290 = t146 * t260 - t150 * t262;
t63 = t142 * t244 - t245 * t290;
t432 = t61 + t63;
t473 = -t261 * t432 + t263 * t431;
t472 = t261 * t431 + t263 * t432;
t331 = qJD(1) * t244 + qJD(5);
t126 = -t332 * t391 + (-t263 * t331 - t344) * t260;
t127 = t331 * t390 + (t245 * t373 - t465) * t261;
t69 = Icges(7,5) * t127 + Icges(7,6) * t126 + Icges(7,3) * t271;
t73 = Icges(7,4) * t127 + Icges(7,2) * t126 + Icges(7,6) * t271;
t77 = Icges(7,1) * t127 + Icges(7,4) * t126 + Icges(7,5) * t271;
t19 = (qJD(4) * t292 + t69) * t244 + (qJD(4) * t140 - t260 * t73 + t262 * t77 + (-t144 * t262 - t148 * t260) * qJD(5)) * t245;
t71 = Icges(6,5) * t127 + Icges(6,6) * t126 + Icges(6,3) * t271;
t75 = Icges(6,4) * t127 + Icges(6,2) * t126 + Icges(6,6) * t271;
t79 = Icges(6,1) * t127 + Icges(6,4) * t126 + Icges(6,5) * t271;
t21 = (qJD(4) * t290 + t71) * t244 + (qJD(4) * t142 - t260 * t75 + t262 * t79 + (-t146 * t262 - t150 * t260) * qJD(5)) * t245;
t471 = t19 + t21;
t343 = t245 * t372;
t460 = t261 * t331 - t343;
t124 = -t460 * t260 + t332 * t390;
t125 = t460 * t262 + t263 * t465;
t68 = Icges(7,5) * t125 + Icges(7,6) * t124 - Icges(7,3) * t476;
t72 = Icges(7,4) * t125 + Icges(7,2) * t124 - Icges(7,6) * t476;
t76 = Icges(7,1) * t125 + Icges(7,4) * t124 - Icges(7,5) * t476;
t20 = (qJD(4) * t291 + t68) * t244 + (qJD(4) * t141 - t260 * t72 + t262 * t76 + (-t145 * t262 - t149 * t260) * qJD(5)) * t245;
t70 = Icges(6,5) * t125 + Icges(6,6) * t124 - Icges(6,3) * t476;
t74 = Icges(6,4) * t125 + Icges(6,2) * t124 - Icges(6,6) * t476;
t78 = Icges(6,1) * t125 + Icges(6,4) * t124 - Icges(6,5) * t476;
t22 = (qJD(4) * t289 + t70) * t244 + (qJD(4) * t143 - t260 * t74 + t262 * t78 + (-t147 * t262 - t151 * t260) * qJD(5)) * t245;
t470 = t20 + t22;
t243 = pkin(5) * t262 + pkin(4);
t369 = qJD(5) * t262;
t365 = pkin(5) * t369;
t469 = t127 * rSges(7,1) + t126 * rSges(7,2) + t271 * rSges(7,3) + t477 * t243 + t258 * t350 + t263 * t365;
t436 = pkin(4) - t243;
t468 = t244 * t436;
t422 = Icges(5,4) * t245;
t307 = Icges(5,1) * t244 + t422;
t183 = Icges(5,5) * t263 + t261 * t307;
t406 = t183 * t244;
t423 = Icges(5,4) * t244;
t303 = Icges(5,2) * t245 + t423;
t181 = Icges(5,6) * t263 + t261 * t303;
t409 = t181 * t245;
t287 = t406 + t409;
t467 = t263 * t287;
t327 = rSges(5,1) * t244 + rSges(5,2) * t245;
t276 = t263 * t327;
t440 = pkin(4) * t244;
t233 = t263 * t440;
t322 = -t202 * rSges(7,1) - t201 * rSges(7,2);
t435 = -pkin(8) - t258;
t340 = t435 * t245;
t402 = t243 * t244;
t387 = pkin(5) * t392 + t233 + (t340 - t402) * t263 + rSges(7,3) * t397 - t322;
t466 = t263 * t387;
t400 = t244 * t261;
t464 = t200 * rSges(7,1) + t199 * rSges(7,2) + pkin(5) * t394 + t243 * t400 - t399 * t424;
t256 = sin(pkin(9));
t441 = pkin(3) * t256;
t239 = t263 * t441;
t250 = t263 * qJ(2);
t259 = -pkin(7) - qJ(3);
t356 = t261 * t259 + t239 + t250;
t443 = -rSges(5,3) - pkin(1);
t159 = t261 * t443 + t276 + t356;
t251 = t263 * rSges(5,3);
t185 = rSges(5,1) * t400 + rSges(5,2) * t399 + t251;
t252 = t263 * pkin(1);
t381 = t261 * qJ(2) + t252;
t283 = -t259 * t263 + t261 * t441 + t381;
t160 = t283 + t185;
t463 = -t159 * t261 + t160 * t263;
t382 = qJ(2) * t377 + qJD(2) * t261;
t355 = qJD(3) * t263 + t382;
t309 = qJD(1) * t239 + t259 * t378 + t355;
t359 = t477 * rSges(5,1) + rSges(5,2) * t350;
t105 = (-rSges(5,2) * t376 + qJD(1) * t443) * t261 + t309 + t359;
t428 = rSges(5,2) * t244;
t214 = rSges(5,1) * t245 - t428;
t237 = t259 * t377;
t248 = qJD(2) * t263;
t333 = -qJD(3) * t261 + t248;
t329 = t237 + t333;
t339 = -qJ(2) - t441;
t106 = t214 * t372 + (t443 * t263 + (-t327 + t339) * t261) * qJD(1) + t329;
t462 = -t105 * t263 + t261 * t106;
t461 = -t245 * t424 + t402;
t299 = Icges(5,5) * t244 + Icges(5,6) * t245;
t459 = -Icges(5,3) * t261 + t263 * t299;
t458 = -Icges(5,6) * t261 + t263 * t303;
t457 = -Icges(5,5) * t261 + t263 * t307;
t216 = pkin(4) * t245 + pkin(8) * t244;
t206 = t261 * t216;
t321 = rSges(7,1) * t262 - rSges(7,2) * t260;
t385 = (t321 - t436) * t245 + (rSges(7,3) + t435) * t244;
t336 = t385 * t261;
t107 = t206 + t336;
t108 = (-t216 - t385) * t263;
t195 = t216 * t378;
t208 = (pkin(8) * t245 - t440) * qJD(4);
t335 = qJD(1) * t385;
t370 = qJD(5) * t260;
t342 = t245 * t370;
t389 = -pkin(5) * t342 + qJD(6) * t244 + (-rSges(7,1) * t260 - rSges(7,2) * t262) * t371 + (rSges(7,3) * t245 - t244 * t321 + t340 + t468) * qJD(4);
t47 = t195 + t261 * t335 + (-t208 - t389) * t263;
t266 = t261 * t389 + t263 * t335;
t384 = t261 * t208 + t216 * t377;
t48 = t266 + t384;
t456 = t48 * t261 - t263 * t47 + (t107 * t263 + t108 * t261) * qJD(1);
t264 = -t258 * t376 - t474;
t43 = (-pkin(1) * qJD(1) + t264) * t261 + t309 + t469;
t323 = t125 * rSges(7,1) + t124 * rSges(7,2);
t354 = -pkin(5) * t260 - pkin(1);
t401 = t243 * t245;
t44 = t237 + t248 + (-qJD(3) - t365) * t261 + (-pkin(5) * t244 * t370 - t368 + (t244 * t424 + t401) * qJD(4)) * t263 + (t354 * t263 + (t339 - t461) * t261) * qJD(1) - t323;
t93 = t354 * t261 + t461 * t263 + t322 + t356;
t94 = t283 + t464;
t455 = t261 * t44 - t263 * t43 + (t261 * t94 + t263 * t93) * qJD(1);
t334 = qJD(4) * t385;
t357 = pkin(4) * t343 + t476 * pkin(8);
t430 = (t365 + (t245 * t258 - t468) * qJD(1)) * t261 + ((t244 * t258 - t401) * qJD(4) + t474) * t263 + t357 - rSges(7,3) * t476 + t323;
t23 = (-t263 * t334 - t430) * t244 + (-qJD(4) * t387 + t263 * t389 - t378 * t385) * t245;
t232 = pkin(4) * t400;
t193 = -pkin(8) * t399 + t232;
t388 = -t193 + t464;
t358 = t477 * pkin(4) + pkin(8) * t347;
t269 = pkin(8) * t350 - t358;
t429 = -t261 * t264 - t269 - t469;
t24 = (-t261 * t334 - t429) * t244 + (qJD(4) * t388 + t266) * t245;
t59 = t244 * t388 + t245 * t336;
t60 = -t244 * t387 + t385 * t397;
t454 = t23 * t261 - t24 * t263 + (t261 * t59 + t263 * t60) * qJD(1);
t453 = t261 * t387 + t263 * t388;
t452 = 2 * m(5);
t451 = 2 * m(6);
t450 = 2 * m(7);
t254 = t261 ^ 2;
t255 = t263 ^ 2;
t449 = -t244 / 0.2e1;
t447 = t245 / 0.2e1;
t444 = rSges(3,2) - pkin(1);
t442 = m(5) * t214;
t427 = rSges(6,3) * t245;
t426 = t261 * rSges(5,3);
t425 = rSges(4,3) + qJ(3);
t325 = -t202 * rSges(6,1) - t201 * rSges(6,2);
t155 = rSges(6,3) * t397 - t325;
t413 = t155 * t263;
t410 = t181 * t244;
t408 = t458 * t244;
t407 = t458 * t245;
t405 = t183 * t245;
t404 = t457 * t244;
t403 = t457 * t245;
t383 = t200 * rSges(6,1) + t199 * rSges(6,2);
t153 = -rSges(6,3) * t399 + t383;
t386 = -t153 - t193;
t380 = t254 + t255;
t179 = Icges(5,3) * t263 + t261 * t299;
t379 = qJD(1) * t179;
t367 = -pkin(1) - t425;
t366 = m(7) * t376;
t35 = t52 * t261 + t263 * t51;
t36 = t54 * t261 + t263 * t53;
t364 = t35 / 0.2e1 + t36 / 0.2e1;
t37 = t56 * t261 + t263 * t55;
t38 = t58 * t261 + t263 * t57;
t363 = -t37 / 0.2e1 - t38 / 0.2e1;
t362 = t127 * rSges(6,1) + t126 * rSges(6,2) + rSges(6,3) * t347;
t361 = -t193 - t388;
t353 = t245 * (-rSges(6,3) - pkin(8));
t324 = rSges(6,1) * t262 - rSges(6,2) * t260;
t178 = rSges(6,3) * t244 + t245 * t324;
t349 = t178 * t378;
t338 = t481 * t375 + (((t479 * t260 + t480 * t262) * qJD(5) + t475) * t245 + t478) * t244;
t207 = t327 * qJD(4);
t337 = t207 * t380;
t328 = rSges(4,1) * t256 + rSges(4,2) * cos(pkin(9));
t326 = t125 * rSges(6,1) + t124 * rSges(6,2);
t313 = t261 * t60 - t263 * t59;
t310 = t261 * t93 - t263 * t94;
t308 = Icges(5,1) * t245 - t423;
t304 = -Icges(5,2) * t244 + t422;
t300 = Icges(5,5) * t245 - Icges(5,6) * t244;
t295 = t107 * t261 - t108 * t263;
t288 = t153 * t263 + t155 * t261;
t286 = -t404 - t407;
t33 = -t118 * t399 + t199 * t120 + t200 * t122 + t126 * t173 + t127 * t175 + t171 * t271;
t34 = -t119 * t399 + t199 * t121 + t200 * t123 + t126 * t174 + t127 * t176 + t172 * t271;
t282 = -t21 / 0.2e1 - t19 / 0.2e1 - t33 / 0.2e1 - t34 / 0.2e1;
t31 = t118 * t397 + t201 * t120 + t202 * t122 + t124 * t173 + t125 * t175 - t171 * t476;
t32 = t119 * t397 + t201 * t121 + t202 * t123 + t124 * t174 + t125 * t176 - t172 * t476;
t281 = t22 / 0.2e1 + t20 / 0.2e1 + t31 / 0.2e1 + t32 / 0.2e1;
t280 = -t64 / 0.2e1 - t62 / 0.2e1 - t86 / 0.2e1 - t87 / 0.2e1;
t279 = t84 / 0.2e1 + t85 / 0.2e1 + t63 / 0.2e1 + t61 / 0.2e1;
t278 = rSges(3,3) * t263 + t261 * t444;
t129 = (-rSges(6,1) * t260 - rSges(6,2) * t262) * t371 + (-t244 * t324 + t427) * qJD(4);
t277 = t261 * t129 + t178 * t377;
t275 = t286 * t261;
t274 = qJD(4) * t308;
t273 = qJD(4) * t304;
t272 = -t261 * pkin(1) + t263 * t353;
t265 = t261 * t367 + t263 * t328;
t198 = -rSges(3,2) * t263 + t261 * rSges(3,3) + t381;
t197 = t250 + t278;
t194 = pkin(8) * t397 - t233;
t187 = t263 * t194;
t186 = t426 - t276;
t170 = t248 + (t444 * t263 + (-rSges(3,3) - qJ(2)) * t261) * qJD(1);
t169 = qJD(1) * t278 + t382;
t168 = t261 * t328 + t263 * t425 + t381;
t167 = t250 + t265;
t158 = t263 * (qJD(1) * t232 - t357);
t157 = (-t178 - t216) * t263;
t156 = t178 * t261 + t206;
t135 = t459 * qJD(1) + t300 * t374;
t134 = -t300 * t372 + t379;
t131 = (t367 * t263 + (-qJ(2) - t328) * t261) * qJD(1) + t333;
t130 = qJD(1) * t265 + t355;
t104 = t261 * t353 + t232 + t283 + t383;
t103 = t233 + t272 + t325 + t356;
t102 = -t244 * t155 + t178 * t397;
t101 = t153 * t244 + t178 * t399;
t100 = -t261 * t459 - t263 * t286;
t99 = t261 * t179 - t467;
t98 = -t263 * t459 + t275;
t97 = t179 * t263 + t261 * t287;
t90 = t288 * t245;
t89 = t277 + t384;
t88 = t349 + t195 + (-t129 - t208) * t263;
t83 = -rSges(6,3) * t350 + t362;
t81 = -rSges(6,3) * t476 + t326;
t65 = t261 * t386 + t187 + t413;
t50 = rSges(6,3) * t345 + (-t252 + (t339 + t427 - t440) * t261) * qJD(1) - t326 + t329 + t357;
t49 = qJD(1) * t272 + t309 + t358 + t362;
t46 = t453 * t245;
t45 = t261 * t361 + t187 + t466;
t42 = (-t178 * t374 + t83) * t244 + (qJD(4) * t153 + t277) * t245;
t41 = (-t178 * t372 - t81) * t244 + (-qJD(4) * t155 + t129 * t263 - t349) * t245;
t30 = t288 * t376 + (-t261 * t81 - t263 * t83 + (t261 * t153 - t413) * qJD(1)) * t245;
t29 = t263 * t81 + t158 + (t269 - t83) * t261 + (t386 * t263 + (-t155 - t194) * t261) * qJD(1);
t18 = t126 * t147 + t127 * t151 + t143 * t271 + t199 * t74 + t200 * t78 - t399 * t70;
t17 = t126 * t146 + t127 * t150 + t142 * t271 + t199 * t75 + t200 * t79 - t399 * t71;
t16 = t126 * t145 + t127 * t149 + t141 * t271 + t199 * t72 + t200 * t76 - t399 * t68;
t15 = t126 * t144 + t127 * t148 + t140 * t271 + t199 * t73 + t200 * t77 - t399 * t69;
t14 = t124 * t147 + t125 * t151 - t143 * t476 + t201 * t74 + t202 * t78 + t397 * t70;
t13 = t124 * t146 + t125 * t150 - t142 * t476 + t201 * t75 + t202 * t79 + t397 * t71;
t12 = t124 * t145 + t125 * t149 - t141 * t476 + t201 * t72 + t202 * t76 + t397 * t68;
t11 = t124 * t144 + t125 * t148 - t140 * t476 + t201 * t73 + t202 * t77 + t397 * t69;
t10 = t158 + t430 * t263 + (t269 + t429) * t261 + (t361 * t263 + (-t194 - t387) * t261) * qJD(1);
t9 = t453 * t376 + (t429 * t263 - t430 * t261 + (t261 * t388 - t466) * qJD(1)) * t245;
t8 = -qJD(1) * t316 + t17 * t263 + t18 * t261;
t7 = -qJD(1) * t317 + t15 * t263 + t16 * t261;
t6 = -qJD(1) * t314 + t13 * t263 + t14 * t261;
t5 = -qJD(1) * t315 + t11 * t263 + t12 * t261;
t4 = (qJD(4) * t316 + t34) * t244 + (-qJD(1) * t36 + qJD(4) * t85 - t17 * t261 + t18 * t263) * t245;
t3 = (qJD(4) * t317 + t33) * t244 + (-qJD(1) * t35 + qJD(4) * t84 - t15 * t261 + t16 * t263) * t245;
t2 = (qJD(4) * t314 + t32) * t244 + (-qJD(1) * t38 + qJD(4) * t87 - t13 * t261 + t14 * t263) * t245;
t1 = (qJD(4) * t315 + t31) * t244 + (-qJD(1) * t37 + qJD(4) * t86 - t11 * t261 + t12 * t263) * t245;
t25 = [-t244 * t274 - t307 * t375 + t303 * t376 + 0.2e1 * m(3) * (t169 * t198 + t170 * t197) + (t105 * t160 + t106 * t159) * t452 + 0.2e1 * m(4) * (t130 * t168 + t131 * t167) + (t103 * t50 + t104 * t49) * t451 + (t43 * t94 + t44 * t93) * t450 + t479 * t342 + (t480 * t369 - t273 + t475) * t245 + t478; m(7) * t455 + m(6) * (t261 * t50 - t263 * t49 + (t103 * t263 + t104 * t261) * qJD(1)) + m(5) * ((t159 * t263 + t160 * t261) * qJD(1) + t462) + m(3) * (-t169 * t263 + t261 * t170 + (t197 * t263 + t198 * t261) * qJD(1)) + m(4) * (-t130 * t263 + t261 * t131 + (t167 * t263 + t168 * t261) * qJD(1)); 0; m(7) * (-qJD(1) * t310 + t261 * t43 + t263 * t44) + m(6) * (t261 * t49 + t263 * t50 + (-t103 * t261 + t104 * t263) * qJD(1)) + m(5) * (t463 * qJD(1) + t261 * t105 + t106 * t263) + m(4) * (t261 * t130 + t131 * t263 + (-t167 * t261 + t168 * t263) * qJD(1)); 0; 0; ((t458 * qJD(1) + t261 * t273) * t449 + (t457 * qJD(1) + t261 * t274) * t447 + (-t409 / 0.2e1 - t406 / 0.2e1) * qJD(4) - t282) * t263 + ((qJD(1) * t181 - t304 * t372) * t449 + (qJD(1) * t183 - t308 * t372) * t447 + (t407 / 0.2e1 + t404 / 0.2e1) * qJD(4) + t281) * t261 + m(5) * (t463 * t207 + t462 * t214) + m(6) * (t103 * t89 + t104 * t88 + t156 * t50 + t157 * t49) + m(7) * (t107 * t44 + t108 * t43 + t47 * t94 + t48 * t93) - (t255 / 0.2e1 + t254 / 0.2e1) * t299 * qJD(4) + ((t160 * t442 + t410 / 0.2e1 - t405 / 0.2e1 - t279) * t261 + (t159 * t442 + t408 / 0.2e1 - t403 / 0.2e1 - t280) * t263) * qJD(1); m(6) * (t89 * t261 - t263 * t88 + (t156 * t263 + t157 * t261) * qJD(1)) + m(7) * t456 - m(5) * t337; m(6) * (t88 * t261 + t263 * t89 + (-t156 * t261 + t157 * t263) * qJD(1)) + m(7) * (-qJD(1) * t295 + t47 * t261 + t263 * t48); (t10 * t45 + t107 * t48 + t108 * t47) * t450 + t261 * t6 + t261 * t5 + t263 * t8 + (t156 * t89 + t157 * t88 + t29 * t65) * t451 + t263 * t7 + t263 * ((t263 * t135 + (t98 + t467) * qJD(1)) * t263 + (-t97 * qJD(1) + (-t375 * t457 + t376 * t458) * t261 + (t134 + (t405 - t410) * qJD(4) + (-t179 + t286) * qJD(1)) * t263) * t261) + t261 * ((t261 * t134 + (-t99 + t275) * qJD(1)) * t261 + (t100 * qJD(1) + (t181 * t376 - t183 * t375 + t379) * t263 + (t135 + (t403 - t408) * qJD(4) + t287 * qJD(1)) * t261) * t263) + ((-t261 * t185 + t186 * t263) * (-t261 * t359 + (-t214 * t255 + t254 * t428) * qJD(4) + ((-t185 + t251) * t263 + (-t186 + t276 + t426) * t261) * qJD(1)) - t214 * t337) * t452 + (-t98 * t261 - t263 * t97 - t35 - t36) * t378 + (t100 * t261 + t263 * t99 + t37 + t38) * t377; m(6) * (t101 * t49 + t102 * t50 + t103 * t41 + t104 * t42) + m(7) * (t23 * t93 + t24 * t94 + t43 * t59 + t44 * t60) + (t261 * t279 + t263 * t280) * t376 + (t281 * t263 + t282 * t261 + (t261 * t280 - t263 * t279) * qJD(1)) * t245 + t338; m(6) * (t41 * t261 - t263 * t42 + (t101 * t261 + t102 * t263) * qJD(1)) + m(7) * t454; m(6) * (t42 * t261 + t263 * t41 + (t101 * t263 - t102 * t261) * qJD(1)) + m(7) * (-qJD(1) * t313 + t23 * t263 + t24 * t261); m(6) * (t101 * t88 + t102 * t89 + t156 * t41 + t157 * t42 - t29 * t90 + t30 * t65) + m(7) * (-t10 * t46 + t107 * t23 + t108 * t24 + t45 * t9 + t47 * t59 + t48 * t60) + (t3 / 0.2e1 + t4 / 0.2e1 + t363 * t376 + t482 * qJD(1) / 0.2e1) * t263 + (t1 / 0.2e1 + t2 / 0.2e1 + t364 * t376 - t483 * qJD(1) / 0.2e1) * t261 + ((t261 * t363 - t263 * t364) * qJD(1) - (t7 + t8) * t261 / 0.2e1 + (t5 + t6) * t263 / 0.2e1 + t472 * qJD(4) / 0.2e1) * t245 + (t473 * qJD(1) + t470 * t261 + t471 * t263) * t244 / 0.2e1; (t23 * t60 + t24 * t59 - t46 * t9) * t450 + (t101 * t42 + t102 * t41 - t30 * t90) * t451 + (((-t244 * t431 - t482) * t263 + (t244 * t432 + t483) * t261) * qJD(4) + t338) * t244 + ((t1 + t2) * t263 + (-t3 - t4) * t261 + (-t471 * t261 + t470 * t263) * t244 + (t481 * t244 + t473 * t245) * qJD(4) + (-t472 * t244 - t261 * t482 - t263 * t483) * qJD(1)) * t245; m(7) * (-t455 * t245 + t310 * t376); t380 * t366; 0; m(7) * ((qJD(4) * t295 + t10) * t244 + (qJD(4) * t45 - t456) * t245); m(7) * ((qJD(4) * t313 + t9) * t244 + (-qJD(4) * t46 - t454) * t245); 0.2e1 * (0.1e1 - t380) * t245 * t366;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t25(1) t25(2) t25(4) t25(7) t25(11) t25(16); t25(2) t25(3) t25(5) t25(8) t25(12) t25(17); t25(4) t25(5) t25(6) t25(9) t25(13) t25(18); t25(7) t25(8) t25(9) t25(10) t25(14) t25(19); t25(11) t25(12) t25(13) t25(14) t25(15) t25(20); t25(16) t25(17) t25(18) t25(19) t25(20) t25(21);];
Mq  = res;
