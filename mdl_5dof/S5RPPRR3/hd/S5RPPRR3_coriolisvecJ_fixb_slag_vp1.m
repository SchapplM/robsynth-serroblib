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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:42
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:41:37
% EndTime: 2019-12-05 17:41:55
% DurationCPUTime: 10.34s
% Computational Cost: add. (13681->548), mult. (9868->702), div. (0->0), fcn. (7604->10), ass. (0->322)
t233 = pkin(9) + qJ(4);
t230 = qJ(5) + t233;
t223 = sin(t230);
t224 = cos(t230);
t168 = rSges(6,1) * t223 + rSges(6,2) * t224;
t235 = qJ(1) + pkin(8);
t229 = cos(t235);
t385 = t223 * t229;
t200 = rSges(6,2) * t385;
t227 = sin(t235);
t383 = t224 * t229;
t277 = -rSges(6,1) * t383 - rSges(6,3) * t227;
t117 = -t200 - t277;
t421 = rSges(6,2) * t223;
t424 = rSges(6,1) * t224;
t169 = -t421 + t424;
t453 = -t229 * rSges(6,3) + t169 * t227;
t481 = t229 * t117 + t227 * t453;
t386 = t223 * t227;
t194 = Icges(6,4) * t386;
t384 = t224 * t227;
t401 = Icges(6,5) * t229;
t115 = -Icges(6,1) * t384 + t194 + t401;
t195 = Icges(6,4) * t385;
t402 = Icges(6,5) * t227;
t116 = Icges(6,1) * t383 - t195 + t402;
t234 = qJD(4) + qJD(5);
t175 = t227 * t234;
t176 = t229 * t234;
t217 = Icges(6,4) * t224;
t290 = -Icges(6,2) * t223 + t217;
t463 = Icges(6,1) * t223 + t217;
t471 = t463 + t290;
t247 = qJD(1) * t471 + t175 * (-Icges(6,2) * t383 + t116 - t195) + t176 * (Icges(6,2) * t384 + t115 + t194);
t405 = Icges(6,4) * t223;
t164 = Icges(6,2) * t224 + t405;
t167 = Icges(6,1) * t224 - t405;
t475 = t164 - t167;
t397 = Icges(6,6) * t227;
t114 = Icges(6,4) * t383 - Icges(6,2) * t385 + t397;
t476 = t229 * t463 + t114;
t113 = Icges(6,6) * t229 - t227 * t290;
t477 = -t227 * t463 + t113;
t445 = qJD(1) * t475 + t175 * t476 + t176 * t477;
t480 = t247 * t223 + t224 * t445;
t226 = sin(t233);
t228 = cos(t233);
t382 = t226 * t227;
t204 = Icges(5,4) * t382;
t379 = t227 * t228;
t403 = Icges(5,5) * t229;
t124 = -Icges(5,1) * t379 + t204 + t403;
t381 = t226 * t229;
t205 = Icges(5,4) * t381;
t378 = t228 * t229;
t404 = Icges(5,5) * t227;
t125 = Icges(5,1) * t378 - t205 + t404;
t259 = t227 * (-Icges(5,2) * t378 + t125 - t205) + t229 * (Icges(5,2) * t379 + t124 + t204);
t221 = Icges(5,4) * t228;
t291 = -Icges(5,2) * t226 + t221;
t122 = Icges(5,6) * t229 - t227 * t291;
t398 = Icges(5,6) * t227;
t123 = Icges(5,4) * t378 - Icges(5,2) * t381 + t398;
t461 = Icges(5,1) * t226 + t221;
t143 = t461 * t227;
t144 = t461 * t229;
t449 = t227 * (t123 + t144) + t229 * (t122 - t143);
t479 = -t259 * t226 - t228 * t449;
t478 = 2 * qJD(4);
t184 = rSges(5,1) * t226 + rSges(5,2) * t228;
t349 = qJD(4) * t229;
t332 = t184 * t349;
t136 = t168 * t227;
t353 = qJD(1) * t227;
t474 = qJD(1) * t136 - t168 * t353 - t169 * t176;
t389 = t168 * t229;
t249 = qJD(1) * t453;
t74 = -t234 * t389 - t249;
t472 = t168 * t234;
t338 = qJD(1) * t200 + t227 * t472;
t75 = qJD(1) * t277 + t338;
t473 = t136 * t175 + t176 * t389 - t227 * t75 + t229 * t74;
t362 = t461 + t291;
t406 = Icges(5,4) * t226;
t179 = Icges(5,2) * t228 + t406;
t182 = Icges(5,1) * t228 - t406;
t363 = t179 - t182;
t470 = (t226 * t362 + t228 * t363) * qJD(1);
t151 = t169 * t234;
t352 = qJD(1) * t229;
t351 = qJD(4) * t227;
t331 = t226 * t351;
t196 = pkin(4) * t331;
t220 = qJD(3) * t229;
t360 = t196 + t220;
t301 = t168 * t175 + t360;
t429 = pkin(6) + qJ(3);
t215 = t227 * t429;
t237 = cos(pkin(9));
t377 = t229 * t237;
t380 = t227 * qJ(3);
t276 = -pkin(3) * t377 + t380;
t110 = -t215 + t276;
t186 = pkin(2) * t229 + t380;
t239 = cos(qJ(1));
t432 = t239 * pkin(1);
t322 = -t186 - t432;
t306 = t110 + t322;
t435 = pkin(3) * t237;
t336 = pkin(2) + t435;
t282 = pkin(4) * t228 + t336;
t155 = t229 * t282;
t232 = -pkin(7) - t429;
t303 = t229 * t336;
t90 = t227 * t232 - t155 + t215 + t303;
t468 = -t117 + t90;
t40 = (t306 + t468) * qJD(1) + t301;
t469 = (-qJD(1) * t389 + t151 * t227 + t168 * t352 - t169 * t175) * t40;
t395 = Icges(6,3) * t227;
t112 = Icges(6,5) * t383 - Icges(6,6) * t385 + t395;
t457 = -t112 * t229 + t116 * t384;
t44 = t114 * t386 - t457;
t286 = t164 * t223 - t224 * t463;
t162 = Icges(6,5) * t223 + Icges(6,6) * t224;
t390 = t162 * t229;
t57 = t227 * t286 + t390;
t467 = qJD(1) * t57 + t175 * t44;
t287 = t123 * t226 - t125 * t228;
t466 = t229 * t287;
t219 = qJD(3) * t227;
t318 = qJD(1) * t429;
t361 = -t229 * t318 + t336 * t353;
t462 = -t219 + t361;
t300 = -rSges(3,1) * t229 - t432;
t460 = rSges(3,2) * t227 + t300;
t348 = pkin(4) * t379;
t265 = (-t232 - t429) * t229 - t348;
t253 = qJD(1) * t265;
t459 = -t229 * t253 + t353 * t90;
t422 = rSges(5,2) * t226;
t426 = rSges(5,1) * t228;
t294 = -t422 + t426;
t417 = t229 * rSges(5,3);
t452 = t227 * t294 - t417;
t250 = qJD(1) * t452;
t238 = sin(qJ(1));
t433 = t238 * pkin(1);
t341 = qJD(1) * t433;
t280 = -qJD(1) * (-pkin(2) * t227 + t229 * qJ(3)) - t219 + t341;
t274 = -qJD(1) * (pkin(6) * t229 - t227 * t435) + t280;
t52 = t250 + t274 + t332;
t206 = rSges(5,2) * t381;
t278 = -rSges(5,1) * t378 - t227 * rSges(5,3);
t126 = -t206 - t278;
t154 = t184 * t351;
t53 = t154 + t220 + (-t126 + t306) * qJD(1);
t456 = t227 * t53 + t229 * t52;
t396 = Icges(5,3) * t227;
t121 = Icges(5,5) * t378 - Icges(5,6) * t381 + t396;
t62 = t123 * t228 + t125 * t226;
t81 = qJD(1) * t122 - t179 * t349;
t83 = -qJD(4) * t144 + (-t182 * t227 + t403) * qJD(1);
t455 = -qJD(1) * t121 + qJD(4) * t62 + t226 * t81 - t228 * t83;
t178 = Icges(5,5) * t228 - Icges(5,6) * t226;
t120 = Icges(5,3) * t229 - t178 * t227;
t356 = qJD(1) * t120;
t61 = t122 * t228 + t124 * t226;
t82 = t179 * t351 + (-t229 * t291 - t398) * qJD(1);
t84 = qJD(4) * t143 + (-t182 * t229 - t404) * qJD(1);
t454 = qJD(4) * t61 + t226 * t82 - t228 * t84 - t356;
t159 = t291 * qJD(4);
t160 = t182 * qJD(4);
t177 = Icges(5,5) * t226 + Icges(5,6) * t228;
t451 = -qJD(1) * t177 + qJD(4) * (t179 * t228 + t226 * t461) + t159 * t226 - t160 * t228;
t309 = t475 * t234;
t310 = t471 * t234;
t448 = -qJD(1) * t162 + t223 * t310 + t224 * t309;
t314 = qJD(1) * t113 + t116 * t234 - t164 * t176;
t316 = -(-t167 * t227 + t401) * qJD(1) + t476 * t234;
t447 = -qJD(1) * t112 + t223 * t314 + t224 * t316;
t315 = t115 * t234 + t164 * t175 + (-t229 * t290 - t397) * qJD(1);
t317 = -(-t167 * t229 - t402) * qJD(1) + t477 * t234;
t163 = Icges(6,5) * t224 - Icges(6,6) * t223;
t111 = Icges(6,3) * t229 - t163 * t227;
t358 = qJD(1) * t111;
t446 = t223 * t315 + t224 * t317 - t358;
t444 = t229 ^ 2;
t241 = qJD(1) ^ 2;
t156 = t234 * t353;
t443 = -t156 / 0.2e1;
t157 = qJD(1) * t176;
t442 = t157 / 0.2e1;
t441 = -t175 / 0.2e1;
t440 = t175 / 0.2e1;
t439 = -t176 / 0.2e1;
t438 = t176 / 0.2e1;
t437 = t227 / 0.2e1;
t436 = t229 / 0.2e1;
t434 = pkin(4) * t226;
t431 = -qJD(1) / 0.2e1;
t430 = qJD(1) / 0.2e1;
t427 = rSges(4,1) * t237;
t423 = rSges(4,2) * sin(pkin(9));
t414 = rSges(4,3) + qJ(3);
t413 = -t111 * t227 - t115 * t383;
t412 = t111 * t229 + t113 * t386;
t411 = t120 * t227 + t124 * t378;
t410 = t120 * t229 + t122 * t382;
t409 = -t232 + rSges(6,3);
t394 = t112 * t227;
t393 = t122 * t226;
t392 = t124 * t228;
t130 = t162 * t227;
t139 = t177 * t227;
t388 = t177 * t229;
t146 = t184 * t227;
t387 = t184 * t229;
t58 = -t229 * t286 + t130;
t376 = t58 * qJD(1);
t284 = t179 * t226 - t228 * t461;
t73 = -t229 * t284 + t139;
t375 = t73 * qJD(1);
t174 = qJD(1) * t186;
t364 = t174 - t220;
t354 = qJD(1) * t178;
t350 = qJD(4) * t228;
t347 = pkin(4) * t378;
t346 = t241 * t432;
t343 = qJD(4) * t434;
t342 = pkin(4) * t349;
t211 = t229 * t423;
t337 = rSges(5,1) * t331 + (t226 * t352 + t227 * t350) * rSges(5,2);
t329 = -pkin(2) - t427;
t328 = -t353 / 0.2e1;
t327 = t352 / 0.2e1;
t326 = -t351 / 0.2e1;
t324 = -t349 / 0.2e1;
t320 = t168 + t434;
t312 = -t121 - t392;
t307 = t226 * t342;
t218 = pkin(2) * t353;
t299 = qJ(3) * t352 - t218;
t145 = t219 + t299;
t302 = -t145 + t299 + t462;
t296 = rSges(3,1) * t227 + rSges(3,2) * t229;
t295 = -t423 + t427;
t56 = t114 * t224 + t116 * t223;
t289 = t114 * t223 - t116 * t224;
t288 = -t392 + t393;
t283 = -t346 + (t220 - t364) * qJD(1);
t281 = -t336 - t426;
t45 = -t113 * t385 - t413;
t49 = t121 * t229 + t123 * t382 - t125 * t379;
t279 = -rSges(4,1) * t377 - rSges(4,3) * t227;
t209 = t227 * t318;
t275 = t283 + qJD(1) * (qJD(1) * t276 - t209);
t48 = -t124 * t379 + t410;
t273 = (t227 * t49 + t229 * t48) * qJD(4);
t50 = -t122 * t381 + t411;
t51 = t121 * t227 - t466;
t272 = (t227 * t51 + t229 * t50) * qJD(4);
t268 = -t282 - t424;
t267 = qJD(1) * t163 + t130 * t176 - t175 * t390;
t264 = qJD(1) * t289 - t234 * t390 + t358;
t263 = t234 * t130 + (t113 * t223 - t115 * t224 - t163 * t229 - t395) * qJD(1);
t262 = qJD(1) * t287 - qJD(4) * t388 + t356;
t261 = qJD(4) * t139 + (-t178 * t229 + t288 - t396) * qJD(1);
t258 = qJD(1) * t286 + t163 * t234;
t257 = t284 * qJD(1) + qJD(4) * t178;
t10 = t227 * t264 - t229 * t447;
t11 = t227 * t446 + t229 * t263;
t12 = t227 * t447 + t229 * t264;
t43 = -t115 * t384 + t412;
t20 = t176 * t43 + t467;
t46 = -t229 * t289 + t394;
t21 = t175 * t46 + t176 * t45 + t376;
t28 = t227 * t258 - t229 * t448;
t29 = t227 * t448 + t229 * t258;
t30 = -t223 * t317 + t224 * t315;
t31 = -t223 * t316 + t224 * t314;
t55 = t113 * t224 + t115 * t223;
t9 = t227 * t263 - t229 * t446;
t256 = (qJD(1) * t28 + t10 * t175 - t156 * t45 + t157 * t46 + t176 * t9) * t437 + (-t223 * t445 + t224 * t247) * t431 + t20 * t328 + t21 * t327 + (qJD(1) * t29 + t11 * t176 + t12 * t175 - t156 * t43 + t157 * t44) * t436 + (t227 * t44 + t229 * t43) * t443 + (t227 * t46 + t229 * t45) * t442 + (t10 * t227 + t229 * t9 + (-t227 * t45 + t229 * t46) * qJD(1)) * t440 + (t11 * t229 + t12 * t227 + (-t227 * t43 + t229 * t44) * qJD(1)) * t438 + (t227 * t31 + t229 * t30 + (-t227 * t55 + t229 * t56) * qJD(1)) * t430 + (t267 * t227 - t229 * t480) * t441 + (t227 * t480 + t267 * t229) * t439;
t254 = t227 * t265;
t252 = -t268 - t421;
t244 = t229 * t126 + t227 * t452;
t243 = t114 * t385 - t394 + (-t115 * t227 - t116 * t229) * t224 + t412;
t85 = -qJD(4) * t387 - t250;
t86 = qJD(1) * t278 + t337;
t242 = -t227 * t86 + t229 * t85 + (-t444 * rSges(5,3) + (t229 * t294 - t126) * t227) * qJD(1);
t240 = qJD(4) ^ 2;
t225 = t241 * t433;
t216 = rSges(3,2) * t353;
t208 = t232 * t353;
t207 = qJD(1) * t211;
t161 = t294 * qJD(4);
t129 = -t211 - t279;
t108 = qJD(1) * t110;
t88 = t220 + (-t129 + t322) * qJD(1);
t87 = -qJD(1) * (rSges(4,3) * t229 - t227 * t295) + t280;
t77 = -qJD(1) * t347 + t196 + t208 + t209;
t76 = -t307 + (-t227 * t282 - t229 * t232) * qJD(1) + t361;
t72 = t227 * t284 + t388;
t65 = t72 * qJD(1);
t64 = (qJD(1) * t279 + t207) * qJD(1) + t283;
t63 = t225 + (-rSges(4,3) * t352 - t145 + (qJD(1) * t295 - qJD(3)) * t227) * qJD(1);
t59 = qJD(4) * t244 + qJD(2);
t42 = -t161 * t349 + (t86 + t154) * qJD(1) + t275;
t41 = t161 * t351 + t225 + (t302 - t85 + t332) * qJD(1);
t39 = t168 * t176 + t249 - t253 + t274 + t307;
t37 = t227 * t451 + t229 * t257;
t36 = t227 * t257 - t229 * t451;
t35 = -qJD(4) * t254 + t117 * t176 + t175 * t453 - t349 * t90 + qJD(2);
t34 = -qJD(4) * t287 + t226 * t83 + t228 * t81;
t33 = -qJD(4) * t288 + t226 * t84 + t228 * t82;
t32 = t242 * qJD(4);
t25 = -t240 * t347 - t151 * t176 + t156 * t168 + (t75 + t77 + t196) * qJD(1) + t275;
t24 = t240 * t348 + t151 * t175 + t157 * t168 + t225 + (t302 - t74 - t76 + t307) * qJD(1);
t23 = t272 + t375;
t22 = t65 + t273;
t8 = qJD(4) * t459 - t156 * t117 + t157 * t453 - t175 * t75 + t176 * t74 + t349 * t76 - t351 * t77;
t1 = [(-t284 * qJD(4) + t159 * t228 + t160 * t226) * qJD(1) + ((t46 + t243) * t176 + t467) * t441 + (-t309 * t223 + t310 * t224) * qJD(1) + m(3) * ((t241 * t296 + t225) * t460 - (qJD(1) * t296 + t341) * (qJD(1) * t300 + t216) + (t346 + (rSges(3,1) * t352 + qJD(1) * t460 - t216) * qJD(1)) * (t296 + t433)) + (t55 + t57) * t443 + (t56 + t58) * t442 + (-t376 + (t44 + (t113 * t229 - t114 * t227) * t223 + t413 + t457) * t176 + (-t43 + t243) * t175 + t21) * t439 + (t30 + t29) * t438 + (t23 - t375 + ((t49 + (-t121 + t393) * t229 - t411) * t229 + (t227 * t312 + t410 - t48) * t227) * qJD(4)) * t324 + (-(-t108 + t174 + t40 + (t432 - t468) * qJD(1) - t301) * t39 + t24 * (-t155 + t200 - t432) - t40 * t219 - t25 * t433 - t39 * (t208 + t338 + t360) + (-t24 * t424 + t40 * (t343 + t472) + t25 * t409) * t229 + (-t24 * t409 - t25 * t252) * t227 + ((t238 * t40 + t239 * t39) * pkin(1) + (-t268 * t39 - t40 * t409) * t229 + (rSges(6,3) * t39 + t252 * t40) * t227) * qJD(1)) * m(6) + (t41 * (t206 - t215 - t432) - t42 * t433 + (-t41 * rSges(5,3) + t42 * (t281 + t422)) * t227 + (t41 * t281 + t42 * (rSges(5,3) + t429)) * t229 + (t462 + t332 + (rSges(5,1) * t379 - rSges(5,2) * t382 - t417 + t433) * qJD(1)) * t53 + (t108 + t154 - t53 - t364 + t209 - t220 - t337 + (-t126 - t278 + t303) * qJD(1)) * t52) * m(5) + (-(t88 + (t129 + t432) * qJD(1) + t364) * t87 + t63 * (t211 - t432) + t88 * (t218 - t219) - t64 * t433 - t87 * (t207 + t220) + (t63 * t329 + t64 * t414) * t229 + (-t63 * t414 + t64 * (-pkin(2) - t295)) * t227 + ((t238 * t88 + t239 * t87) * pkin(1) + (-t329 * t87 - t414 * t88) * t229 + (t295 * t88 + t414 * t87) * t227) * qJD(1)) * m(4) + (t31 + t28 + t20) * t440 + (t34 + t36 + t22) * t351 / 0.2e1 + (t65 + ((t410 + t51 + t466) * t229 + (-t50 + (t312 - t393) * t229 + t49 + t411) * t227) * qJD(4) + (t61 + t72) * qJD(1)) * t326 + (t33 + t37 + (t62 + t73) * qJD(1)) * t349 / 0.2e1; m(5) * t32 + m(6) * t8; m(4) * (t227 * t64 + t229 * t63) + m(5) * (t227 * t42 + t229 * t41) + m(6) * (t227 * t25 + t229 * t24); ((-t226 * t363 + t228 * t362) * qJD(1) + (-t226 * t449 + t228 * t259) * qJD(4)) * t431 + ((t139 * t349 + t354) * t229 + (t470 + (-t229 * t388 - t479) * qJD(4)) * t227) * t324 + t256 + (t227 * t34 + t229 * t33 + (-t227 * t61 + t229 * t62) * qJD(1)) * t430 + ((-t351 * t388 + t354) * t227 + (-t470 + (t227 * t139 + t479) * qJD(4)) * t229) * t326 + (qJD(1) * t36 + ((t227 * t261 - t229 * t454) * t229 + (t227 * t262 - t229 * t455) * t227 + (-t227 * t50 + t229 * t51) * qJD(1)) * t478) * t437 + (qJD(1) * t37 + ((t227 * t454 + t229 * t261) * t229 + (t227 * t455 + t229 * t262) * t227 + (-t227 * t48 + t229 * t49) * qJD(1)) * t478) * t436 + (t273 + t22) * t328 + (t272 + t23) * t327 + (t8 * (-t254 + t481) + t24 * t320 * t227 + (-t25 * t320 - t8 * t90) * t229 + (-(-t227 ^ 2 - t444) * t343 - t117 * t353 - t227 * t77 + t459 + (t249 + t76) * t229 + t473) * t35 + t469 + (-t228 * t342 + (pkin(4) * t350 + t151) * t229 + t474) * t39) * m(6) + (-(-t146 * t52 + t387 * t53) * qJD(1) - (t59 * (-t146 * t227 - t229 * t387) + t456 * t294) * qJD(4) + t41 * t146 + t59 * t242 + t32 * t244 - t42 * t387 + (t352 * t53 - t353 * t52) * t184 + t456 * t161) * m(5); t256 + (t8 * t481 + t24 * t136 - t25 * t389 + ((-t444 * rSges(6,3) + (t169 * t229 - t117) * t227) * qJD(1) + t473) * t35 + t469 + (t151 * t229 + t474) * t39) * m(6);];
tauc = t1(:);
