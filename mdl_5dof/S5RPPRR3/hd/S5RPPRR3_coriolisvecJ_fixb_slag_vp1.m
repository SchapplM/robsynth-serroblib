% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:08
% EndTime: 2022-01-23 09:14:23
% DurationCPUTime: 11.88s
% Computational Cost: add. (13713->546), mult. (9886->710), div. (0->0), fcn. (7620->10), ass. (0->327)
t255 = qJ(1) + pkin(8);
t246 = sin(t255);
t253 = pkin(9) + qJ(4);
t247 = cos(t253);
t391 = t246 * t247;
t353 = rSges(5,1) * t391;
t248 = cos(t255);
t258 = -pkin(6) - qJ(3);
t217 = t248 * t258;
t257 = cos(pkin(9));
t244 = t257 * pkin(3) + pkin(2);
t371 = t246 * t244 + t217;
t259 = sin(qJ(1));
t439 = pkin(1) * t259;
t475 = -t353 - t439 - t371;
t226 = t248 * qJ(3);
t185 = pkin(2) * t246 - t226;
t114 = t185 - t371;
t174 = qJD(1) * t185;
t223 = qJD(3) * t246;
t474 = -qJD(1) * t114 + t174 + t223;
t234 = t246 * rSges(6,3);
t249 = qJ(5) + t253;
t243 = cos(t249);
t394 = t243 * t248;
t373 = rSges(6,1) * t394 + t234;
t242 = sin(t249);
t396 = t242 * t248;
t123 = -rSges(6,2) * t396 + t373;
t368 = -t248 * t244 + t246 * t258;
t438 = pkin(4) * t247;
t197 = t244 + t438;
t252 = pkin(7) - t258;
t213 = t252 * t246;
t468 = t248 * t197 + t213;
t459 = t368 + t468;
t473 = -t123 - t459;
t262 = qJD(1) ^ 2;
t397 = t242 * t246;
t195 = rSges(6,2) * t397;
t374 = t248 * rSges(6,3) + t195;
t395 = t243 * t246;
t122 = rSges(6,1) * t395 - t374;
t94 = t197 * t246 - t248 * t252 - t371;
t472 = -t122 - t94;
t409 = Icges(6,6) * t248;
t118 = Icges(6,4) * t395 - Icges(6,2) * t397 - t409;
t222 = Icges(6,4) * t243;
t170 = Icges(6,1) * t242 + t222;
t471 = -t170 * t246 - t118;
t300 = -Icges(6,2) * t242 + t222;
t119 = Icges(6,6) * t246 + t248 * t300;
t470 = -t170 * t248 - t119;
t415 = Icges(6,4) * t242;
t171 = Icges(6,1) * t243 - t415;
t121 = Icges(6,5) * t246 + t171 * t248;
t168 = Icges(6,2) * t243 + t415;
t469 = -t168 * t248 + t121;
t467 = -t168 + t171;
t466 = t170 + t300;
t430 = rSges(4,2) * sin(pkin(9));
t434 = rSges(4,1) * t257;
t290 = t246 * rSges(4,3) + (-t430 + t434) * t248;
t225 = t246 * qJ(3);
t188 = t248 * pkin(2) + t225;
t260 = cos(qJ(1));
t251 = t260 * pkin(1);
t335 = t188 + t251;
t465 = t290 + t335;
t314 = -t188 - t368 + t335;
t464 = -t314 + t473;
t361 = qJD(4) * t248;
t463 = 0.2e1 * qJD(4);
t245 = sin(t253);
t462 = rSges(5,2) * t245;
t393 = t245 * t246;
t460 = -rSges(5,2) * t393 - t248 * rSges(5,3);
t233 = Icges(5,4) * t247;
t301 = -Icges(5,2) * t245 + t233;
t181 = Icges(5,1) * t245 + t233;
t331 = t248 * rSges(3,1) - rSges(3,2) * t246;
t458 = t251 + t331;
t178 = Icges(5,5) * t247 - Icges(5,6) * t245;
t177 = Icges(5,5) * t245 + Icges(5,6) * t247;
t281 = qJD(4) * t177;
t416 = Icges(5,4) * t245;
t182 = Icges(5,1) * t247 - t416;
t130 = Icges(5,5) * t246 + t182 * t248;
t128 = Icges(5,6) * t246 + t248 * t301;
t402 = t128 * t245;
t296 = -t130 * t247 + t402;
t408 = Icges(5,3) * t248;
t456 = -t248 * t281 + (-t178 * t246 + t296 + t408) * qJD(1);
t201 = Icges(5,4) * t393;
t414 = Icges(5,5) * t248;
t129 = Icges(5,1) * t391 - t201 - t414;
t410 = Icges(5,6) * t248;
t127 = Icges(5,4) * t391 - Icges(5,2) * t393 - t410;
t403 = t127 * t245;
t297 = -t129 * t247 + t403;
t126 = Icges(5,3) * t246 + t178 * t248;
t365 = qJD(1) * t126;
t455 = qJD(1) * t297 - t246 * t281 + t365;
t167 = Icges(6,5) * t243 - Icges(6,6) * t242;
t254 = qJD(4) + qJD(5);
t166 = Icges(6,5) * t242 + Icges(6,6) * t243;
t400 = t166 * t248;
t405 = t119 * t242;
t407 = Icges(6,3) * t248;
t454 = -t254 * t400 + (-t121 * t243 - t167 * t246 + t405 + t407) * qJD(1);
t192 = Icges(6,4) * t397;
t413 = Icges(6,5) * t248;
t120 = Icges(6,1) * t395 - t192 - t413;
t299 = t118 * t242 - t120 * t243;
t117 = Icges(6,3) * t246 + t167 * t248;
t366 = qJD(1) * t117;
t401 = t166 * t246;
t453 = qJD(1) * t299 - t254 * t401 + t366;
t125 = Icges(5,5) * t391 - Icges(5,6) * t393 - t408;
t50 = -t125 * t248 - t246 * t297;
t294 = t168 * t242 - t170 * t243;
t452 = qJD(1) * t294 + t167 * t254;
t179 = Icges(5,2) * t247 + t416;
t293 = t245 * t179 - t247 * t181;
t451 = t293 * qJD(1) + t178 * qJD(4);
t427 = rSges(6,2) * t243;
t433 = rSges(6,1) * t242;
t172 = t427 + t433;
t141 = t172 * t246;
t142 = t172 * t248;
t428 = rSges(6,2) * t242;
t432 = rSges(6,1) * t243;
t173 = -t428 + t432;
t175 = t246 * t254;
t176 = t248 * t254;
t35 = t122 * t175 + t123 * t176 + qJD(2) + (t246 * t94 + t248 * t459) * qJD(4);
t336 = -t185 - t439;
t315 = t114 + t336;
t346 = t245 * t361;
t318 = pkin(4) * t346;
t320 = -t172 * t176 + t223;
t40 = -t318 + (t315 + t472) * qJD(1) + t320;
t426 = pkin(4) * qJD(4);
t350 = t245 * t426;
t193 = t246 * t350;
t224 = qJD(3) * t248;
t375 = t193 + t224;
t310 = t172 * t175 + t375;
t41 = -t464 * qJD(1) - t310;
t317 = qJD(1) * t254;
t160 = t246 * t317;
t161 = t248 * t317;
t351 = t254 * t427;
t364 = qJD(1) * t246;
t363 = qJD(1) * t248;
t377 = rSges(6,3) * t363 + qJD(1) * t195;
t77 = -t248 * t351 + (-t176 * t242 - t243 * t364) * rSges(6,1) + t377;
t78 = -t254 * t141 + (t173 * t248 + t234) * qJD(1);
t207 = t252 * t363;
t372 = t197 - t244;
t79 = -t318 + t207 + (-t246 * t372 + t217) * qJD(1);
t208 = t258 * t364;
t80 = -t193 + t208 + (t248 * t372 + t213) * qJD(1);
t8 = t122 * t161 - t123 * t160 + t175 * t78 + t176 * t77 + (t246 * t80 + t248 * t79 + (-t246 * t459 + t248 * t94) * qJD(1)) * qJD(4);
t450 = -(qJD(1) * t141 - t176 * t173) * t40 - t35 * (-t175 * t141 - t142 * t176) - t41 * (-qJD(1) * t142 - t173 * t175) + t8 * (t246 * t122 + t248 * t123);
t449 = t246 * (-t179 * t248 + t130) - t248 * (-Icges(5,2) * t391 + t129 - t201);
t448 = qJD(1) * t466 + t175 * t469 - t176 * (-Icges(6,2) * t395 + t120 - t192);
t447 = t160 / 0.2e1;
t446 = t161 / 0.2e1;
t445 = -t175 / 0.2e1;
t444 = t175 / 0.2e1;
t443 = -t176 / 0.2e1;
t442 = t176 / 0.2e1;
t441 = t246 / 0.2e1;
t440 = -t248 / 0.2e1;
t437 = -qJD(1) / 0.2e1;
t436 = qJD(1) / 0.2e1;
t435 = pkin(2) - t244;
t429 = rSges(5,2) * t247;
t184 = rSges(5,1) * t245 + t429;
t152 = t184 * t248;
t235 = t246 * rSges(5,3);
t390 = t247 * t248;
t392 = t245 * t248;
t132 = rSges(5,1) * t390 - rSges(5,2) * t392 + t235;
t362 = qJD(4) * t246;
t347 = t184 * t362;
t55 = -t347 - t224 + (t132 + t314) * qJD(1);
t425 = t152 * t55;
t424 = t245 * t41;
t131 = t353 + t460;
t345 = t184 * t361;
t309 = t223 - t345;
t54 = (-t131 + t315) * qJD(1) + t309;
t423 = t246 * t54;
t422 = t40 * t172;
t421 = t40 * t248;
t116 = Icges(6,5) * t395 - Icges(6,6) * t397 - t407;
t420 = -t246 * t116 - t120 * t394;
t419 = t246 * t117 + t121 * t394;
t406 = t116 * t248;
t399 = t177 * t246;
t398 = t177 * t248;
t59 = -t246 * t294 - t400;
t389 = t59 * qJD(1);
t75 = -t246 * t293 - t398;
t388 = t75 * qJD(1);
t387 = -t246 * t125 - t129 * t390;
t386 = t246 * t126 + t130 * t390;
t150 = qJD(1) * t188 - t224;
t385 = t208 - (-t248 * t435 - t225) * qJD(1) - t150;
t379 = -t179 + t182;
t378 = t181 + t301;
t376 = rSges(5,3) * t363 + t364 * t462;
t209 = t246 * t430;
t370 = rSges(4,3) * t363 + qJD(1) * t209;
t369 = t248 * rSges(4,3) + t209;
t218 = qJ(3) * t363;
t367 = t218 + t223;
t360 = t178 * qJD(1);
t359 = qJD(1) * qJD(3);
t358 = qJD(4) ^ 2 * t438;
t357 = t262 * t439;
t356 = t262 * t251;
t355 = t122 * t363 + t246 * t78 + t248 * t77;
t354 = t246 * t434;
t344 = -pkin(2) - t434;
t342 = t364 / 0.2e1;
t341 = t363 / 0.2e1;
t340 = -t362 / 0.2e1;
t337 = t361 / 0.2e1;
t334 = -pkin(4) * t245 - t172;
t332 = -t197 - t432;
t97 = t121 * t395;
t330 = t117 * t248 - t97;
t329 = qJD(1) * t121 + t254 * t471;
t328 = (-t171 * t246 + t413) * qJD(1) + t470 * t254;
t327 = qJD(1) * t119 + t120 * t254 - t168 * t175;
t326 = (-t246 * t300 + t409) * qJD(1) + t469 * t254;
t102 = t130 * t391;
t325 = t126 * t248 - t102;
t324 = -t116 + t405;
t323 = -t125 + t402;
t322 = t466 * t254;
t321 = t467 * t254;
t316 = t248 * t359 - t356;
t156 = t173 * t254;
t311 = -t247 * t426 - t156;
t186 = rSges(3,1) * t246 + rSges(3,2) * t248;
t306 = rSges(5,1) * t247 - t462;
t305 = -t41 * t246 - t421;
t304 = -t55 * t246 - t54 * t248;
t57 = t118 * t243 + t120 * t242;
t64 = t127 * t247 + t129 * t245;
t65 = t128 * t247 + t130 * t245;
t295 = t131 * t246 + t132 * t248;
t292 = qJD(1) * (-pkin(2) * t364 + t367) + t246 * t359 - t357;
t288 = qJD(1) * (-t218 + (t246 * t435 - t217) * qJD(1)) + t292;
t151 = t184 * t246;
t51 = -t128 * t393 - t325;
t286 = (t246 * t51 - t248 * t50) * qJD(4);
t52 = -t127 * t392 - t387;
t53 = -t128 * t392 + t386;
t285 = (t246 * t53 - t248 * t52) * qJD(4);
t284 = t299 * t246;
t283 = qJD(4) * t181;
t282 = qJD(4) * t179;
t279 = qJD(1) * t167 - t175 * t400 + t176 * t401;
t278 = t127 * t248 - t128 * t246;
t269 = -t242 * t326 + t243 * t328 + t366;
t10 = t246 * t454 + t269 * t248;
t270 = qJD(1) * t116 - t242 * t327 + t243 * t329;
t11 = t270 * t246 - t248 * t453;
t12 = t269 * t246 - t248 * t454;
t44 = -t284 - t406;
t45 = -t119 * t397 - t330;
t20 = t175 * t45 - t176 * t44 + t389;
t46 = -t118 * t396 - t420;
t47 = -t119 * t396 + t419;
t60 = -t248 * t294 + t401;
t56 = t60 * qJD(1);
t21 = t175 * t47 - t176 * t46 + t56;
t273 = qJD(1) * t467 + t175 * t470 - t176 * t471;
t263 = -t242 * t448 + t273 * t243;
t268 = qJD(1) * t166 - t242 * t322 + t243 * t321;
t28 = t246 * t452 + t268 * t248;
t29 = t268 * t246 - t248 * t452;
t30 = t242 * t329 + t243 * t327;
t31 = t242 * t328 + t243 * t326;
t58 = t119 * t243 + t121 * t242;
t9 = t246 * t453 + t270 * t248;
t277 = (qJD(1) * t28 + t10 * t175 + t160 * t46 + t161 * t47 - t176 * t9) * t441 + (t273 * t242 + t243 * t448) * t437 + t20 * t342 + t21 * t341 + (qJD(1) * t29 - t11 * t176 + t12 * t175 + t160 * t44 + t161 * t45) * t440 + (t246 * t45 - t248 * t44) * t447 + (t246 * t47 - t248 * t46) * t446 + (t10 * t246 - t248 * t9 + (t246 * t46 + t248 * t47) * qJD(1)) * t444 + (-t11 * t248 + t12 * t246 + (t246 * t44 + t248 * t45) * qJD(1)) * t443 + (t246 * t31 - t248 * t30 + (t246 * t57 + t248 * t58) * qJD(1)) * t436 + (t246 * t279 + t248 * t263) * t445 + (t246 * t263 - t248 * t279) * t442;
t276 = (-t245 * t378 + t247 * t379) * qJD(1);
t89 = -t361 * t429 + (-t247 * t364 - t346) * rSges(5,1) + t376;
t90 = -qJD(4) * t151 + (t248 * t306 + t235) * qJD(1);
t274 = t246 * t90 + t248 * t89 + (t131 * t248 - t132 * t246) * qJD(1);
t86 = qJD(1) * t128 - t246 * t282;
t88 = qJD(1) * t130 - t246 * t283;
t267 = qJD(1) * t125 - qJD(4) * t64 - t245 * t86 + t247 * t88;
t85 = -t248 * t282 + (-t246 * t301 + t410) * qJD(1);
t87 = -t248 * t283 + (-t182 * t246 + t414) * qJD(1);
t266 = -qJD(4) * t65 - t245 * t85 + t247 * t87 + t365;
t163 = t301 * qJD(4);
t164 = t182 * qJD(4);
t265 = qJD(1) * t177 - t163 * t245 + t164 * t247 + (-t179 * t247 - t181 * t245) * qJD(4);
t264 = -t245 * t449 + t278 * t247;
t165 = t306 * qJD(4);
t134 = t354 - t369;
t92 = qJD(1) * t465 - t224;
t91 = t223 + (-t134 + t336) * qJD(1);
t76 = -t248 * t293 + t399;
t68 = t76 * qJD(1);
t67 = (-qJD(1) * t290 - t150) * qJD(1) + t316;
t66 = qJD(1) * (-qJD(1) * t354 + t370) + t292;
t61 = qJD(4) * t295 + qJD(2);
t43 = -t165 * t361 + (-t90 + t347 + t385) * qJD(1) + t316;
t42 = -t165 * t362 + (t89 - t345) * qJD(1) + t288;
t37 = t265 * t246 - t248 * t451;
t36 = t246 * t451 + t265 * t248;
t34 = -qJD(4) * t296 + t245 * t87 + t247 * t85;
t33 = -t297 * qJD(4) + t245 * t88 + t247 * t86;
t32 = t274 * qJD(4);
t25 = -t248 * t358 - t156 * t176 + t160 * t172 + (-t78 - t80 + t193 + t385) * qJD(1) + t316;
t24 = -t246 * t358 - t156 * t175 - t161 * t172 + (t77 + t79 - t318) * qJD(1) + t288;
t23 = t68 + t285;
t22 = t286 + t388;
t1 = [m(3) * ((-t186 * t262 - t357) * t458 + (-t356 + (-0.2e1 * t331 - t251 + t458) * t262) * (-t186 - t439)) + (t68 + ((t51 - t102 + (t126 + t403) * t248 + t387) * t248 + t386 * t246) * qJD(4)) * t337 + (t56 + (t45 + (t118 * t248 + t119 * t246) * t242 + t330 + t420) * t176 + (-t120 * t395 + t406 + t44 + (t118 * t246 - t119 * t248) * t242 + t419) * t175) * t442 + (t57 + t59) * t447 + (t58 + t60) * t446 + (-t389 + (t47 - t284 - t419) * t176 + (t324 * t246 + t46 - t97) * t175 + ((t117 + t299) * t175 + t324 * t176) * t248 + t20) * t445 + (t31 + t28) * t444 + (-t388 + ((t248 * t323 - t386 + t53) * t248 + (t246 * t323 + t325 + t52) * t246) * qJD(4) + t22) * t340 + (t34 + t36) * t362 / 0.2e1 + (-t293 * qJD(4) + t163 * t247 + t164 * t245 + t321 * t242 + t322 * t243) * qJD(1) + (pkin(4) * t424 * t361 + t25 * (t374 - t439) + t24 * (t251 + t373 + t468) + (-t24 * t428 + t25 * t252) * t248 + (t25 * t332 + t254 * t422) * t246 + (-t320 + t207 + t377 + (-t254 * t433 - t350 - t351) * t248 + t474) * t41 + (t375 - t310) * t40 + (-t40 * t464 - t41 * (-t439 + t472) + (-t259 * t41 - t40 * t260) * pkin(1) + (-t173 - t197) * t421 + (t40 * (-rSges(6,3) - t252) + t41 * t332) * t246) * qJD(1)) * m(6) + (t43 * (-t460 + t475) + t42 * (t132 + t251 - t368) + (t184 * t423 - t425) * qJD(4) + (t208 + t224 + (-t235 - t251 + (-t244 - t306) * t248) * qJD(1)) * t54 + (t54 - t309 + t376 + (t131 + t439 + t475) * qJD(1) + t474) * t55) * m(5) + (t67 * (t246 * t344 + t226 + t369 - t439) + t91 * t224 + t66 * t465 + t92 * (t367 + t370) + ((-t259 * t92 - t91 * t260) * pkin(1) + t91 * (t344 + t430) * t248 + (t91 * (-rSges(4,3) - qJ(3)) + t92 * t344) * t246) * qJD(1) - (-t174 + t223 - t91 + (-t134 - t439) * qJD(1)) * t92) * m(4) + (t30 + t29 + t21) * t443 - (t33 + t37 + t23) * t361 / 0.2e1 + ((t64 + t75) * t246 + (t65 + t76) * t248) * qJD(4) * t436; m(5) * t32 + m(6) * t8; 0.2e1 * (t24 * t440 + t25 * t441) * m(6) + 0.2e1 * (t42 * t440 + t43 * t441) * m(5) + 0.2e1 * (t440 * t66 + t441 * t67) * m(4); ((t245 * t379 + t247 * t378) * qJD(1) + (t278 * t245 + t247 * t449) * qJD(4)) * t437 + (t246 * t34 - t248 * t33 + (t64 * t246 + t248 * t65) * qJD(1)) * t436 + ((-t361 * t399 - t360) * t248 + (t276 + (t248 * t398 + t264) * qJD(4)) * t246) * t337 + ((-t362 * t398 + t360) * t246 + (t276 + (t246 * t399 + t264) * qJD(4)) * t248) * t340 + t277 + (qJD(1) * t36 + (-(t246 * t455 + t267 * t248) * t248 + (t246 * t456 + t266 * t248) * t246 + (t52 * t246 + t53 * t248) * qJD(1)) * t463) * t441 + (qJD(1) * t37 + (-(t267 * t246 - t248 * t455) * t248 + (t266 * t246 - t248 * t456) * t246 + (t50 * t246 + t51 * t248) * qJD(1)) * t463) * t440 + (t286 + t22) * t342 + (t23 + t285) * t341 + (t35 * t355 + (t25 * t334 + t40 * t311 + t8 * t459 + t35 * t79 + (t334 * t41 + t35 * t94) * qJD(1)) * t248 + (t24 * t334 + t41 * t311 + t8 * t94 + t35 * t80 + (t35 * t473 + t422) * qJD(1)) * t246 - (-t363 * t424 + (t305 * t247 + t35 * (-t246 ^ 2 - t248 ^ 2) * t245) * qJD(4)) * pkin(4) + t450) * m(6) + (t32 * t295 + t61 * t274 + t304 * t165 + (-t42 * t246 - t43 * t248 + (-t248 * t55 + t423) * qJD(1)) * t184 - (t151 * t54 - t425) * qJD(1) - (t61 * (-t151 * t246 - t152 * t248) + t304 * t306) * qJD(4)) * m(5); t277 + (t35 * (-t123 * t364 + t355) + t305 * t156 + (-t24 * t246 - t25 * t248 + (t246 * t40 - t248 * t41) * qJD(1)) * t172 + t450) * m(6);];
tauc = t1(:);
