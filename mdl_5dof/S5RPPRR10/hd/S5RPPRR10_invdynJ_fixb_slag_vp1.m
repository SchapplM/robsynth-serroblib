% Calculate vector of inverse dynamics joint torques for
% S5RPPRR10
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
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR10_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR10_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:03:59
% EndTime: 2019-12-31 18:04:16
% DurationCPUTime: 13.54s
% Computational Cost: add. (12290->726), mult. (20951->914), div. (0->0), fcn. (20611->8), ass. (0->326)
t332 = cos(qJ(1));
t509 = rSges(6,3) * t332;
t446 = qJ(4) + qJ(5);
t318 = sin(t446);
t328 = cos(pkin(8));
t327 = sin(pkin(8));
t391 = cos(t446);
t375 = t327 * t391;
t253 = -t328 * t318 + t375;
t347 = t327 * t318 + t328 * t391;
t188 = rSges(6,1) * t253 - rSges(6,2) * t347;
t329 = sin(qJ(4));
t449 = t328 * t329;
t405 = pkin(4) * t449;
t331 = cos(qJ(4));
t309 = pkin(4) * t331 + pkin(3);
t469 = pkin(3) - t309;
t226 = -t327 * t469 - t405;
t326 = qJD(4) + qJD(5);
t280 = t326 * t332;
t414 = qJD(4) * t332;
t508 = -t188 * t280 - t226 * t414;
t330 = sin(qJ(1));
t448 = t328 * t330;
t222 = t318 * t448 - t330 * t375;
t223 = t347 * t330;
t130 = -t223 * rSges(6,1) + t222 * rSges(6,2) - t509;
t333 = -pkin(7) - pkin(6);
t408 = pkin(3) * t448;
t453 = t327 * t329;
t407 = pkin(4) * t453;
t427 = -t309 * t448 - t330 * t407;
t177 = t408 + (pkin(6) + t333) * t332 + t427;
t440 = t130 + t177;
t355 = t328 * t331 + t453;
t246 = t355 * t330;
t451 = t327 * t331;
t499 = -t329 * t448 + t330 * t451;
t147 = Icges(5,5) * t246 + Icges(5,6) * t499 + Icges(5,3) * t332;
t228 = Icges(5,4) * t246;
t150 = Icges(5,2) * t499 + Icges(5,6) * t332 + t228;
t227 = Icges(5,4) * t499;
t154 = -Icges(5,1) * t246 - Icges(5,5) * t332 - t227;
t266 = -t449 + t451;
t247 = t266 * t332;
t248 = t355 * t332;
t445 = t247 * t150 - t154 * t248;
t52 = -t147 * t330 + t445;
t454 = t326 * t330;
t120 = Icges(6,5) * t223 - Icges(6,6) * t222 + Icges(6,3) * t332;
t208 = Icges(6,4) * t223;
t123 = -Icges(6,2) * t222 + Icges(6,6) * t332 + t208;
t207 = Icges(6,4) * t222;
t127 = -Icges(6,1) * t223 - Icges(6,5) * t332 + t207;
t224 = t253 * t332;
t225 = t347 * t332;
t48 = -t120 * t330 + t224 * t123 - t127 * t225;
t122 = Icges(6,5) * t225 + Icges(6,6) * t224 - Icges(6,3) * t330;
t459 = Icges(6,4) * t225;
t125 = Icges(6,2) * t224 - Icges(6,6) * t330 + t459;
t209 = Icges(6,4) * t224;
t128 = Icges(6,1) * t225 - Icges(6,5) * t330 + t209;
t49 = -t122 * t330 + t224 * t125 + t225 * t128;
t182 = Icges(6,5) * t253 - Icges(6,6) * t347;
t458 = Icges(6,4) * t253;
t184 = -Icges(6,2) * t347 + t458;
t242 = Icges(6,4) * t347;
t186 = Icges(6,1) * t253 - t242;
t56 = -t182 * t330 + t184 * t224 + t186 * t225;
t17 = t56 * qJD(1) + t280 * t48 - t454 * t49;
t62 = -t123 * t347 - t127 * t253;
t505 = t150 * t499 - t154 * t246;
t69 = -t150 * t355 - t154 * t266;
t46 = t120 * t332 - t123 * t222 - t127 * t223;
t149 = Icges(5,5) * t248 + Icges(5,6) * t247 - Icges(5,3) * t330;
t461 = Icges(5,4) * t248;
t152 = Icges(5,2) * t247 - Icges(5,6) * t330 + t461;
t229 = Icges(5,4) * t247;
t155 = Icges(5,1) * t248 - Icges(5,5) * t330 + t229;
t444 = t247 * t152 + t248 * t155;
t376 = t149 * t330 - t444;
t50 = t147 * t332 + t505;
t504 = t376 - t50;
t406 = pkin(4) * t451;
t500 = (-t405 + t406) * qJD(4);
t497 = t152 * t499 + t246 * t155;
t131 = t225 * rSges(6,1) + t224 * rSges(6,2) - rSges(6,3) * t330;
t145 = -t222 * rSges(6,1) - t223 * rSges(6,2);
t146 = t224 * rSges(6,1) - t225 * rSges(6,2);
t419 = qJD(1) * t330;
t496 = t131 * t419 + t145 * t454 + t146 * t280;
t287 = t332 * pkin(1) + t330 * qJ(2);
t447 = t328 * t332;
t450 = t327 * t332;
t238 = rSges(4,1) * t447 + t330 * rSges(4,2) + rSges(4,3) * t450;
t259 = pkin(2) * t447 + qJ(3) * t450;
t377 = t259 + t287;
t493 = t377 + t238;
t491 = qJD(4) * t266;
t47 = t332 * t122 - t222 * t125 + t223 * t128;
t157 = -t246 * rSges(5,1) - rSges(5,2) * t499 - rSges(5,3) * t332;
t239 = rSges(3,1) * t447 - rSges(3,2) * t450 + t330 * rSges(3,3);
t489 = t328 * t469 - t407;
t488 = t330 * (-Icges(5,2) * t248 + t155 + t229) - t332 * (-Icges(5,2) * t246 - t154 + t227);
t335 = qJD(1) * (-Icges(6,2) * t253 + t186 - t242) - t454 * (-Icges(6,2) * t225 + t128 + t209) + t280 * (-Icges(6,2) * t223 - t127 - t207);
t487 = t328 ^ 2;
t380 = qJD(1) * t326;
t220 = (-qJDD(4) - qJDD(5)) * t330 - t332 * t380;
t486 = t220 / 0.2e1;
t314 = qJDD(4) * t332;
t221 = qJDD(5) * t332 - t330 * t380 + t314;
t485 = t221 / 0.2e1;
t411 = qJD(1) * qJD(4);
t273 = -qJDD(4) * t330 - t332 * t411;
t484 = t273 / 0.2e1;
t274 = -t330 * t411 + t314;
t483 = t274 / 0.2e1;
t481 = -t454 / 0.2e1;
t480 = -t280 / 0.2e1;
t479 = t280 / 0.2e1;
t478 = -t330 / 0.2e1;
t477 = t330 / 0.2e1;
t476 = -t332 / 0.2e1;
t475 = t332 / 0.2e1;
t474 = -rSges(5,3) - pkin(6);
t473 = pkin(6) * t332;
t472 = g(1) * t330;
t471 = -qJD(1) / 0.2e1;
t470 = qJD(1) / 0.2e1;
t418 = qJD(1) * t332;
t210 = t347 * t326;
t116 = -t210 * t332 - t253 * t419;
t117 = t224 * t326 - t347 * t419;
t443 = t117 * rSges(6,1) + t116 * rSges(6,2);
t77 = -rSges(6,3) * t418 + t443;
t349 = t332 * t500 + t333 * t418;
t99 = (t330 * t489 + t473) * qJD(1) + t349;
t468 = t77 + t99;
t463 = qJDD(1) / 0.2e1;
t312 = t330 * t333;
t313 = pkin(6) * t419;
t348 = t330 * t491;
t100 = t313 + pkin(4) * t348 + (-t332 * t489 + t312) * qJD(1);
t118 = qJD(1) * t224 - t210 * t330;
t119 = qJD(1) * t225 + t253 * t454;
t371 = -rSges(6,1) * t119 - rSges(6,2) * t118;
t78 = -rSges(6,3) * t419 - t371;
t462 = -t100 - t78;
t460 = Icges(5,4) * t266;
t457 = qJ(3) * t327;
t195 = Icges(5,5) * t266 - Icges(5,6) * t355;
t197 = -Icges(5,2) * t355 + t460;
t260 = Icges(5,4) * t355;
t199 = Icges(5,1) * t266 - t260;
t67 = t195 * t332 + t197 * t499 + t199 * t246;
t455 = qJD(1) * t67;
t452 = t327 * t330;
t302 = pkin(3) * t447;
t270 = -pkin(6) * t330 + t302;
t404 = pkin(4) * t329 * t332;
t397 = t309 * t447 + t327 * t404 + t312;
t178 = -t270 + t397;
t439 = t131 + t178;
t211 = t253 * t326;
t135 = -rSges(6,1) * t210 - rSges(6,2) * t211;
t256 = t355 * qJD(4);
t249 = pkin(4) * t256;
t438 = t135 - t249;
t230 = t499 * pkin(4);
t437 = -t145 - t230;
t165 = -t256 * t332 - t266 * t419;
t166 = t332 * t491 - t355 * t419;
t434 = t166 * rSges(5,1) + t165 * rSges(5,2);
t431 = -Icges(5,2) * t266 + t199 - t260;
t430 = -Icges(5,1) * t355 - t197 - t460;
t429 = t248 * rSges(5,1) + t247 * rSges(5,2);
t203 = t239 + t287;
t187 = -rSges(6,1) * t347 - t253 * rSges(6,2);
t370 = pkin(2) * t328 + t457;
t258 = t370 * t330;
t320 = t332 * qJ(2);
t285 = pkin(1) * t330 - t320;
t428 = -t258 - t285;
t403 = rSges(3,1) * t448;
t298 = rSges(3,2) * t452;
t423 = t332 * rSges(3,3) + t298;
t237 = t403 - t423;
t426 = -t285 - t237;
t425 = rSges(3,3) * t418 + qJD(1) * t298;
t417 = qJD(3) * t327;
t296 = t332 * t417;
t316 = qJD(2) * t330;
t424 = t296 + t316;
t412 = qJD(1) * qJD(2);
t422 = qJDD(2) * t330 + t332 * t412;
t421 = qJ(2) * t418 + t316;
t420 = -qJD(1) * t285 + t316;
t416 = qJD(3) * t328;
t415 = qJD(4) * t330;
t413 = -m(4) - m(5) - m(6);
t410 = qJDD(3) * t327;
t409 = qJDD(3) * t328;
t51 = t332 * t149 + t497;
t269 = t408 + t473;
t399 = -t269 + t428;
t398 = t270 + t377;
t324 = t332 * rSges(4,2);
t373 = rSges(4,1) * t328 + rSges(4,3) * t327;
t236 = t330 * t373 - t324;
t396 = -t236 + t428;
t395 = t332 * t410 + t422;
t394 = t296 + t421;
t393 = t330 * t417;
t392 = -rSges(3,1) * t328 - pkin(1);
t390 = -t419 / 0.2e1;
t389 = -t418 / 0.2e1;
t388 = -t415 / 0.2e1;
t387 = t415 / 0.2e1;
t386 = -t414 / 0.2e1;
t385 = t414 / 0.2e1;
t383 = qJD(1) * t146 + t187 * t454;
t379 = t157 + t399;
t378 = -qJD(1) * t258 + t296 + t420;
t317 = qJD(2) * t332;
t374 = t317 - t393;
t288 = rSges(2,1) * t332 - rSges(2,2) * t330;
t286 = rSges(2,1) * t330 + rSges(2,2) * t332;
t167 = qJD(1) * t247 - t256 * t330;
t168 = qJD(1) * t248 + t348;
t372 = rSges(5,1) * t168 + rSges(5,2) * t167;
t81 = Icges(5,5) * t166 + Icges(5,6) * t165 - Icges(5,3) * t418;
t82 = Icges(5,5) * t168 + Icges(5,6) * t167 - Icges(5,3) * t419;
t83 = Icges(5,4) * t166 + Icges(5,2) * t165 - Icges(5,6) * t418;
t84 = Icges(5,4) * t168 + Icges(5,2) * t167 - Icges(5,6) * t419;
t85 = Icges(5,1) * t166 + Icges(5,4) * t165 - Icges(5,5) * t418;
t86 = Icges(5,1) * t168 + Icges(5,4) * t167 - Icges(5,5) * t419;
t369 = (-t147 * t418 + t150 * t165 - t154 * t166 + t247 * t84 + t248 * t86 - t330 * t82) * t332 - (-t149 * t418 + t152 * t165 + t155 * t166 + t247 * t83 + t248 * t85 - t330 * t81) * t330;
t368 = (-t147 * t419 + t150 * t167 - t154 * t168 + t246 * t86 + t332 * t82 + t499 * t84) * t332 - (-t149 * t419 + t152 * t167 + t155 * t168 + t246 * t85 + t332 * t81 + t499 * t83) * t330;
t367 = -t330 * t51 + t332 * t50;
t366 = t330 * t376 + t332 * t52;
t201 = rSges(5,1) * t266 - rSges(5,2) * t355;
t193 = t201 * t414;
t64 = qJD(1) * t379 + t193 + t424;
t158 = -rSges(5,3) * t330 + t429;
t65 = -t317 + (qJD(4) * t201 + t417) * t330 + (t158 + t398) * qJD(1);
t365 = t330 * t65 + t332 * t64;
t87 = -rSges(5,3) * t418 + t434;
t88 = -rSges(5,3) * t419 + t372;
t364 = -t330 * t88 - t332 * t87;
t363 = t399 + t440;
t362 = -qJD(1) * t269 + t378;
t357 = t157 * t330 - t158 * t332;
t356 = (Icges(5,5) * t499 - Icges(5,6) * t246) * t332 - (Icges(5,5) * t247 - Icges(5,6) * t248) * t330;
t231 = -t328 * t404 + t332 * t406;
t354 = -pkin(1) - t370;
t254 = qJD(1) * t287 - t317;
t353 = -t370 * t418 - t254 - 0.2e1 * t393;
t352 = -qJDD(2) * t332 + qJD(1) * (-pkin(1) * t419 + t421) + qJDD(1) * t287 + t330 * t412;
t261 = t355 * pkin(4);
t351 = -qJD(1) * t302 + t313 + t353;
t55 = t182 * t332 - t184 * t222 + t186 * t223;
t346 = qJD(1) * (-Icges(6,5) * t347 - Icges(6,6) * t253) + (-Icges(6,5) * t222 - Icges(6,6) * t223) * t280 - (Icges(6,5) * t224 - Icges(6,6) * t225) * t454;
t345 = -pkin(3) * t328 + t354;
t343 = -(Icges(5,1) * t247 - t152 - t461) * t330 + (Icges(5,1) * t499 - t150 - t228) * t332;
t341 = -pkin(1) + (-rSges(4,1) - pkin(2)) * t328 + (-rSges(4,3) - qJ(3)) * t327;
t340 = qJDD(1) * t259 + t330 * t410 + t352 + (-t370 * t419 + 0.2e1 * t296) * qJD(1);
t71 = Icges(6,5) * t117 + Icges(6,6) * t116 - Icges(6,3) * t418;
t73 = Icges(6,4) * t117 + Icges(6,2) * t116 - Icges(6,6) * t418;
t75 = Icges(6,1) * t117 + Icges(6,4) * t116 - Icges(6,5) * t418;
t10 = t116 * t125 + t117 * t128 - t122 * t418 + t224 * t73 + t225 * t75 - t330 * t71;
t72 = Icges(6,5) * t119 + Icges(6,6) * t118 - Icges(6,3) * t419;
t74 = Icges(6,4) * t119 + Icges(6,2) * t118 - Icges(6,6) * t419;
t76 = Icges(6,1) * t119 + Icges(6,4) * t118 - Icges(6,5) * t419;
t11 = t118 * t123 - t119 * t127 - t120 * t419 - t222 * t74 + t223 * t76 + t332 * t72;
t12 = t118 * t125 + t119 * t128 - t122 * t419 - t222 * t73 + t223 * t75 + t332 * t71;
t132 = -Icges(6,5) * t210 - Icges(6,6) * t211;
t133 = -Icges(6,4) * t210 - Icges(6,2) * t211;
t134 = -Icges(6,1) * t210 - Icges(6,4) * t211;
t28 = t116 * t184 + t117 * t186 - t132 * t330 + t133 * t224 + t134 * t225 - t182 * t418;
t29 = t118 * t184 + t119 * t186 + t132 * t332 - t133 * t222 + t134 * t223 - t182 * t419;
t30 = -t123 * t211 + t127 * t210 + t253 * t76 - t347 * t74;
t31 = -t125 * t211 - t128 * t210 + t253 * t75 - t347 * t73;
t336 = -(Icges(6,1) * t224 - t125 - t459) * t454 + (-Icges(6,1) * t222 - t123 - t208) * t280 + (-Icges(6,1) * t347 - t184 - t458) * qJD(1);
t63 = -t125 * t347 + t128 * t253;
t9 = t116 * t123 - t117 * t127 - t120 * t418 + t224 * t74 + t225 * t76 - t330 * t72;
t339 = (qJD(1) * t28 + qJDD(1) * t56 - t10 * t454 + t220 * t49 + t221 * t48 + t280 * t9) * t478 + (qJD(1) * t55 + t280 * t46 - t454 * t47) * t390 + t17 * t389 + (qJD(1) * t29 + qJDD(1) * t55 + t11 * t280 - t12 * t454 + t220 * t47 + t221 * t46) * t475 + (-t330 * t49 + t332 * t48) * t486 + (-t330 * t47 + t332 * t46) * t485 + (-t10 * t330 + t332 * t9 + (-t330 * t48 - t332 * t49) * qJD(1)) * t481 + (-t330 * t63 + t332 * t62) * t463 + (t11 * t332 - t12 * t330 + (-t330 * t46 - t332 * t47) * qJD(1)) * t479 + (t30 * t332 - t31 * t330 + (-t330 * t62 - t332 * t63) * qJD(1)) * t470 + (t224 * t335 + t225 * t336 - t330 * t346) * t454 / 0.2e1 + (-t222 * t335 + t336 * t223 + t346 * t332) * t480 + (t253 * t336 - t335 * t347) * t471;
t338 = -t309 * t328 + t354 - t407;
t337 = -qJD(1) ^ 2 * t269 + qJDD(1) * t270 + t340;
t311 = rSges(4,2) * t418;
t200 = -rSges(5,1) * t355 - rSges(5,2) * t266;
t194 = -Icges(5,5) * t355 - Icges(5,6) * t266;
t192 = -rSges(5,1) * t256 - rSges(5,2) * t491;
t191 = -Icges(5,1) * t256 - Icges(5,4) * t491;
t190 = -Icges(5,4) * t256 - Icges(5,2) * t491;
t189 = -Icges(5,5) * t256 - Icges(5,6) * t491;
t180 = qJD(1) * t203 - t317;
t179 = qJD(1) * t426 + t316;
t176 = rSges(5,1) * t247 - rSges(5,2) * t248;
t175 = rSges(5,1) * t499 - rSges(5,2) * t246;
t144 = t280 * t187;
t110 = qJD(1) * t493 - t374;
t109 = qJD(1) * t396 + t424;
t90 = qJDD(1) * t239 + qJD(1) * (-qJD(1) * t403 + t425) + t352;
t89 = t426 * qJDD(1) + (-qJD(1) * t239 - t254) * qJD(1) + t422;
t79 = qJD(4) * t357 - t416;
t70 = -t152 * t355 + t155 * t266;
t68 = -t195 * t330 + t197 * t247 + t199 * t248;
t66 = t68 * qJD(1);
t58 = qJDD(1) * t238 + (-t373 * t419 + t311) * qJD(1) + t340;
t57 = t396 * qJDD(1) + (-qJD(1) * t238 + t353) * qJD(1) + t395;
t45 = -t416 + t130 * t454 - t131 * t280 + (t177 * t330 - t178 * t332) * qJD(4);
t44 = t188 * t454 - t317 + (qJD(4) * t226 + t417) * t330 + (t398 + t439) * qJD(1);
t43 = qJD(1) * t363 + t424 - t508;
t38 = qJD(4) * t364 - t157 * t273 - t158 * t274 - t409;
t37 = t167 * t197 + t168 * t199 + t189 * t332 + t190 * t499 + t191 * t246 - t195 * t419;
t36 = t165 * t197 + t166 * t199 - t189 * t330 + t190 * t247 + t191 * t248 - t195 * t418;
t35 = -t152 * t491 - t155 * t256 + t266 * t85 - t355 * t83;
t34 = -t150 * t491 + t154 * t256 + t266 * t86 - t355 * t84;
t33 = qJD(1) * t87 + qJDD(1) * t158 + t192 * t415 - t201 * t273 + t337;
t32 = t192 * t414 + t201 * t274 + t379 * qJDD(1) + (t351 - t88) * qJD(1) + t395;
t25 = qJD(4) * t366 + t66;
t24 = qJD(4) * t367 + t455;
t19 = qJD(1) * t468 + qJDD(1) * t439 + t135 * t454 - t188 * t220 - t226 * t273 - t249 * t415 + t337;
t18 = -t249 * t414 + t135 * t280 + t188 * t221 + t226 * t274 + t363 * qJDD(1) + (t351 + t462) * qJD(1) + t395;
t13 = -t409 - t130 * t220 - t131 * t221 - t177 * t273 - t178 * t274 - t454 * t78 - t280 * t77 + (-t100 * t330 - t332 * t99) * qJD(4);
t1 = [-m(2) * (-g(1) * t286 + g(2) * t288) + t17 * t480 + (t66 + (t445 * t332 + (t504 + t505) * t330) * qJD(4)) * t386 + (t63 + t56) * t486 + (t62 + t55) * t485 + (t70 + t68) * t484 + (t69 + t67) * t483 + (t31 + t28) * t481 + (t35 + t36) * t388 + (-t455 + ((t444 + t504) * t332 + t497 * t330) * qJD(4) + t24) * t387 + (t30 + t29 + t17) * t479 + (t34 + t37 + t25) * t385 + (-t133 * t347 + t134 * t253 - t184 * t211 - t186 * t210 - t190 * t355 + t191 * t266 - t197 * t491 - t199 * t256) * qJD(1) + ((t19 - g(2)) * (t131 + t377 + t397) + (t18 - g(1)) * (t332 * t333 + t130 + t320 + t427) + (t18 * t330 - t472) * t354 + (t317 + t371 + (-t417 - t500) * t330 + (t338 * t332 + (rSges(6,3) - qJ(2) - t333) * t330) * qJD(1)) * t43 + (-t362 + t43 + t349 + t394 + t443 + (t338 * t330 - t440 - t509) * qJD(1) + t508) * t44) * m(6) + (-(-t457 - pkin(1) + (-pkin(2) - pkin(3)) * t328) * t472 + t64 * (t313 + t317 - t372) + t65 * (t394 + t434) + (t32 * t345 - t64 * t417) * t330 + ((t345 * t64 + t474 * t65) * t332 + (t64 * (rSges(5,3) - qJ(2)) + t65 * t345) * t330) * qJD(1) - (qJD(1) * t157 + t193 + t362 - t64) * t65 + (-g(2) + t33) * (t330 * t474 + t302 + t377 + t429) + (-g(1) + t32) * (t157 + t320 - t473)) * m(5) + (-(-qJD(1) * t236 - t109 + t378) * t110 + t109 * t374 + t110 * (t311 + t394) + (t109 * t341 * t332 + (t109 * (-rSges(4,2) - qJ(2)) + t110 * (t354 - t373)) * t330) * qJD(1) + (t58 - g(2)) * t493 + (t57 - g(1)) * (t330 * t341 + t320 + t324)) * m(4) + (-(-qJD(1) * t237 - t179 + t420) * t180 + t179 * t317 + t180 * (t421 + t425) + (t179 * (rSges(3,2) * t327 + t392) * t332 + (t179 * (-rSges(3,3) - qJ(2)) + t180 * t392) * t330) * qJD(1) + (t90 - g(2)) * t203 + (t89 - g(1)) * (t392 * t330 + t320 + t423)) * m(3) + (-t197 * t355 + t199 * t266 - t184 * t347 + t186 * t253 + m(2) * (t286 ^ 2 + t288 ^ 2) + Icges(2,3) + (Icges(4,3) + Icges(3,2)) * t487 + ((Icges(3,1) + Icges(4,1)) * t327 + 0.2e1 * (Icges(3,4) - Icges(4,5)) * t328) * t327) * qJDD(1); (-m(3) + t413) * (-g(2) * t332 + t472) + 0.2e1 * (t18 * t477 + t19 * t476) * m(6) + 0.2e1 * (t32 * t477 + t33 * t476) * m(5) + 0.2e1 * (t476 * t58 + t477 * t57) * m(4) + 0.2e1 * (t476 * t90 + t477 * t89) * m(3); t413 * (-g(3) * t328 + (g(1) * t332 + g(2) * t330) * t327) + m(4) * (qJDD(3) * t487 + t450 * t57 + t452 * t58) + m(5) * (t32 * t450 - t328 * t38 + t33 * t452) + m(6) * (-t13 * t328 + t18 * t450 + t19 * t452); t367 * t483 + t366 * t484 + (-t330 * t70 + t332 * t69) * t463 + ((-t194 * t330 + t247 * t431 + t248 * t430) * qJD(1) + (-t247 * t488 + t248 * t343 - t330 * t356) * qJD(4)) * t387 + (qJD(1) * t37 + qJD(4) * t368 + qJDD(1) * t67 + t273 * t51 + t274 * t50) * t475 + (qJD(1) * t36 + qJD(4) * t369 + qJDD(1) * t68 - t273 * t376 + t274 * t52) * t478 + ((t194 * t332 + t246 * t430 + t431 * t499) * qJD(1) + (t343 * t246 + t356 * t332 - t488 * t499) * qJD(4)) * t386 + ((-t52 * t330 + t332 * t376) * qJD(1) + t369) * t388 + ((-t50 * t330 - t51 * t332) * qJD(1) + t368) * t385 + t339 + (-t330 * t35 + t332 * t34 + (-t69 * t330 - t332 * t70) * qJD(1)) * t470 + t24 * t390 + t25 * t389 + ((t266 * t430 - t355 * t431) * qJD(1) + (t266 * t343 + t355 * t488) * qJD(4)) * t471 + (-g(1) * (t231 + t146) + g(2) * t437 - g(3) * (-t261 + t187) - t43 * (qJD(1) * t437 - t261 * t414 + t144) - t44 * (qJD(1) * t231 - t261 * t415 + t383) + (t13 * t440 + t44 * t438) * t330 + (-t13 * t439 + t43 * t438) * t332 + ((-qJD(1) * t43 + t19) * t330 + (qJD(1) * t44 + t18) * t332) * (t188 + t226) + (-(-t230 * t330 - t231 * t332) * qJD(4) + (qJD(1) * t178 + t462) * t330 + (qJD(1) * t440 - t468) * t332 + t496) * t45) * m(6) + (-g(1) * t176 - g(2) * t175 - g(3) * t200 - (-t175 * t64 + t176 * t65) * qJD(1) - (t79 * (-t175 * t330 - t176 * t332) + t365 * t200) * qJD(4) + t38 * t357 + t79 * ((t157 * t332 + t158 * t330) * qJD(1) + t364) + t365 * t192 + (t32 * t332 + t33 * t330 + (-t330 * t64 + t332 * t65) * qJD(1)) * t201) * m(5); t339 + (-g(1) * t146 - g(2) * t145 - g(3) * t187 + t13 * (t130 * t330 - t131 * t332) + (t330 * t44 + t332 * t43) * t135 + (t18 * t332 + t19 * t330 + (-t330 * t43 + t332 * t44) * qJD(1)) * t188 - t43 * (-qJD(1) * t145 + t144) - t44 * t383 + (-t330 * t78 + (qJD(1) * t130 - t77) * t332 + t496) * t45) * m(6);];
tau = t1;
