% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 18:26
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:26:35
% EndTime: 2018-11-23 18:26:46
% DurationCPUTime: 11.65s
% Computational Cost: add. (18885->659), mult. (48237->853), div. (0->0), fcn. (35573->8), ass. (0->316)
t507 = Ifges(6,4) + Ifges(7,4);
t508 = Ifges(6,1) + Ifges(7,1);
t499 = Ifges(6,5) + Ifges(7,5);
t312 = qJD(2) + qJD(3);
t309 = qJD(4) + t312;
t313 = sin(qJ(5));
t317 = cos(qJ(5));
t315 = sin(qJ(3));
t316 = sin(qJ(2));
t319 = cos(qJ(3));
t320 = cos(qJ(2));
t395 = t319 * t320;
t278 = -t315 * t316 + t395;
t266 = t278 * qJD(1);
t279 = t315 * t320 + t319 * t316;
t267 = t279 * qJD(1);
t314 = sin(qJ(4));
t318 = cos(qJ(4));
t338 = t266 * t314 + t318 * t267;
t196 = t309 * t317 - t313 * t338;
t510 = t507 * t196;
t509 = t338 * Ifges(5,1) / 0.2e1;
t506 = Ifges(6,2) + Ifges(7,2);
t498 = Ifges(7,6) + Ifges(6,6);
t197 = t309 * t313 + t317 * t338;
t364 = t318 * t266 - t267 * t314;
t213 = qJD(5) - t364;
t495 = t197 * t508 + t499 * t213 + t510;
t473 = t364 * t313;
t204 = pkin(5) * t473;
t383 = qJD(5) * t313;
t307 = pkin(5) * t383;
t505 = t307 - t204;
t504 = qJ(6) * t473 + t317 * qJD(6);
t503 = t507 * t197;
t462 = -pkin(8) - pkin(7);
t295 = t462 * t320;
t284 = qJD(1) * t295;
t268 = t315 * t284;
t294 = t462 * t316;
t283 = qJD(1) * t294;
t274 = qJD(2) * pkin(2) + t283;
t228 = t319 * t274 + t268;
t260 = t267 * pkin(9);
t194 = t228 - t260;
t183 = pkin(3) * t312 + t194;
t271 = t319 * t284;
t229 = t274 * t315 - t271;
t438 = pkin(9) * t266;
t195 = t229 + t438;
t184 = t314 * t195;
t119 = t183 * t318 - t184;
t304 = -pkin(2) * t320 - pkin(1);
t293 = qJD(1) * t304;
t239 = -t266 * pkin(3) + t293;
t444 = -t317 / 0.2e1;
t108 = -pkin(4) * t309 - t119;
t353 = mrSges(7,1) * t313 + mrSges(7,2) * t317;
t355 = mrSges(6,1) * t313 + mrSges(6,2) * t317;
t85 = -pkin(5) * t196 + qJD(6) + t108;
t489 = t108 * t355 + t85 * t353;
t502 = t239 * mrSges(5,2) - t119 * mrSges(5,3) + Ifges(5,4) * t364 + t309 * Ifges(5,5) - t495 * t444 + t489 + t509;
t501 = Ifges(5,2) / 0.2e1;
t237 = t312 * t278;
t224 = t237 * qJD(1);
t238 = t312 * t279;
t225 = t238 * qJD(1);
t117 = qJD(4) * t338 + t224 * t314 + t318 * t225;
t116 = qJD(4) * t364 + t224 * t318 - t225 * t314;
t83 = qJD(5) * t196 + t116 * t317;
t84 = -qJD(5) * t197 - t116 * t313;
t497 = t498 * t117 + t506 * t84 + t507 * t83;
t496 = t499 * t117 + t507 * t84 + t508 * t83;
t478 = t506 * t196 + t498 * t213 + t503;
t311 = t317 * qJ(6);
t357 = pkin(5) * t338 - t311 * t364;
t300 = pkin(3) * t314 + pkin(10);
t393 = -qJ(6) - t300;
t362 = qJD(5) * t393;
t384 = qJD(4) * t318;
t380 = pkin(3) * t384;
t126 = t194 * t318 - t184;
t169 = pkin(4) * t338 - pkin(10) * t364;
t441 = pkin(3) * t267;
t141 = t169 + t441;
t55 = -t126 * t313 + t317 * t141;
t494 = -t357 - t55 + (-qJD(6) - t380) * t313 + t317 * t362;
t302 = pkin(2) * t319 + pkin(3);
t385 = qJD(4) * t314;
t399 = t314 * t315;
t222 = t302 * t384 + (-t315 * t385 + (t318 * t319 - t399) * qJD(3)) * pkin(2);
t233 = -t283 * t315 + t271;
t198 = t233 - t438;
t234 = t319 * t283 + t268;
t199 = -t260 + t234;
t135 = t198 * t314 + t199 * t318;
t390 = qJD(1) * t316;
t306 = pkin(2) * t390;
t136 = t141 + t306;
t58 = t317 * t135 + t313 * t136;
t493 = t222 * t317 - t58;
t56 = t317 * t126 + t313 * t141;
t492 = t313 * t362 + t317 * t380 + t504 - t56;
t23 = -mrSges(6,1) * t84 + mrSges(6,2) * t83;
t185 = t318 * t195;
t120 = t183 * t314 + t185;
t374 = qJD(2) * t462;
t386 = qJD(3) * t319;
t387 = qJD(3) * t315;
t171 = t267 * t374 + t274 * t386 + t284 * t387;
t124 = -pkin(9) * t225 + t171;
t280 = t315 * t294;
t326 = (t395 * t462 - t280) * qJD(2) * qJD(1);
t172 = -qJD(3) * t229 + t326;
t436 = t224 * pkin(9);
t28 = qJD(4) * t120 + t124 * t314 - t318 * (t172 - t436);
t491 = m(6) * t28 + t23;
t125 = t194 * t314 + t185;
t490 = pkin(3) * t385 - t125 + t505;
t398 = t315 * t318;
t471 = -t318 * t198 + t199 * t314 - t302 * t385 - (t315 * t384 + (t314 * t319 + t398) * qJD(3)) * pkin(2);
t488 = t499 * t313 + t498 * t317;
t342 = Ifges(7,5) * t317 - Ifges(7,6) * t313;
t344 = Ifges(6,5) * t317 - Ifges(6,6) * t313;
t487 = t344 + t342;
t425 = Ifges(7,4) * t313;
t427 = Ifges(6,4) * t313;
t486 = t506 * t317 + t425 + t427;
t424 = Ifges(7,4) * t317;
t346 = -Ifges(7,2) * t313 + t424;
t426 = Ifges(6,4) * t317;
t348 = -Ifges(6,2) * t313 + t426;
t485 = t348 + t346;
t484 = t508 * t313 + t424 + t426;
t350 = Ifges(7,1) * t317 - t425;
t352 = Ifges(6,1) * t317 - t427;
t483 = t352 + t350;
t482 = t309 * Ifges(5,6) / 0.2e1;
t109 = pkin(10) * t309 + t120;
t131 = -pkin(4) * t364 - pkin(10) * t338 + t239;
t52 = -t109 * t313 + t317 * t131;
t38 = -qJ(6) * t197 + t52;
t29 = pkin(5) * t213 + t38;
t53 = t109 * t317 + t131 * t313;
t39 = qJ(6) * t196 + t53;
t419 = t364 * Ifges(5,2);
t428 = Ifges(5,4) * t338;
t467 = t482 + t428 / 0.2e1;
t480 = -t239 * mrSges(5,1) - t52 * mrSges(6,1) - t29 * mrSges(7,1) + t53 * mrSges(6,2) + t39 * mrSges(7,2) + t419 / 0.2e1 + t467;
t451 = -t338 / 0.2e1;
t262 = pkin(2) * t398 + t314 * t302;
t256 = pkin(10) + t262;
t394 = -qJ(6) - t256;
t363 = qJD(5) * t394;
t57 = -t135 * t313 + t317 * t136;
t477 = -t357 - t57 + (-qJD(6) - t222) * t313 + t317 * t363;
t433 = -qJ(6) - pkin(10);
t365 = qJD(5) * t433;
t59 = -t119 * t313 + t317 * t169;
t476 = -qJD(6) * t313 + t317 * t365 - t357 - t59;
t475 = t313 * t363 + t493 + t504;
t60 = t317 * t119 + t313 * t169;
t474 = t313 * t365 + t504 - t60;
t472 = -t471 + t505;
t392 = -mrSges(5,1) * t309 - mrSges(6,1) * t196 + mrSges(6,2) * t197 + mrSges(5,3) * t338;
t470 = -t135 + t222;
t242 = t319 * t294 + t295 * t315;
t210 = -pkin(9) * t279 + t242;
t243 = -t319 * t295 + t280;
t211 = pkin(9) * t278 + t243;
t158 = t210 * t314 + t211 * t318;
t152 = t317 * t158;
t232 = t278 * t314 + t279 * t318;
t248 = -t278 * pkin(3) + t304;
t337 = t318 * t278 - t279 * t314;
t153 = -pkin(4) * t337 - t232 * pkin(10) + t248;
t73 = t313 * t153 + t152;
t469 = t318 * t210 - t211 * t314;
t415 = t313 * t52;
t468 = t317 * t53 - t415;
t328 = t213 * t342;
t329 = t213 * t344;
t330 = t197 * t350;
t331 = t197 * t352;
t332 = t196 * t346;
t333 = t196 * t348;
t339 = t313 * t53 + t317 * t52;
t445 = t313 / 0.2e1;
t466 = -(t29 * t317 + t313 * t39) * mrSges(7,3) - t339 * mrSges(6,3) + t333 / 0.2e1 + t331 / 0.2e1 + t332 / 0.2e1 + t330 / 0.2e1 + t328 / 0.2e1 + t329 / 0.2e1 - t478 * t445 + t502;
t377 = -Ifges(7,3) / 0.2e1 - Ifges(6,3) / 0.2e1;
t378 = Ifges(7,6) / 0.2e1 + Ifges(6,6) / 0.2e1;
t379 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t95 = t197 * Ifges(7,5) + t196 * Ifges(7,6) + t213 * Ifges(7,3);
t96 = t197 * Ifges(6,5) + t196 * Ifges(6,6) + t213 * Ifges(6,3);
t465 = t378 * t196 + t379 * t197 - t377 * t213 - t120 * mrSges(5,3) + t95 / 0.2e1 + t96 / 0.2e1 - t467 - t480;
t464 = t83 / 0.2e1;
t463 = t84 / 0.2e1;
t460 = pkin(1) * mrSges(3,1);
t459 = pkin(1) * mrSges(3,2);
t458 = t117 / 0.2e1;
t455 = -t196 / 0.2e1;
t454 = -t197 / 0.2e1;
t453 = t197 / 0.2e1;
t452 = -t213 / 0.2e1;
t449 = t266 / 0.2e1;
t448 = -t267 / 0.2e1;
t447 = t267 / 0.2e1;
t443 = t317 / 0.2e1;
t442 = m(4) * t293;
t440 = pkin(3) * t318;
t439 = pkin(5) * t317;
t437 = pkin(10) * t317;
t27 = t318 * t124 + t119 * qJD(4) + (-t274 * t387 + t284 * t386 + t326 - t436) * t314;
t382 = qJD(5) * t317;
t202 = pkin(3) * t225 + qJD(2) * t306;
t40 = pkin(4) * t117 - pkin(10) * t116 + t202;
t7 = -t109 * t383 + t131 * t382 + t317 * t27 + t313 * t40;
t435 = t317 * t7;
t432 = mrSges(4,3) * t266;
t431 = mrSges(4,3) * t267;
t430 = Ifges(3,4) * t316;
t429 = Ifges(4,4) * t267;
t422 = t469 * t28;
t413 = Ifges(3,5) * qJD(2);
t412 = Ifges(3,6) * qJD(2);
t411 = qJD(2) * mrSges(3,1);
t410 = qJD(2) * mrSges(3,2);
t409 = t120 * t338;
t404 = t364 * t317;
t402 = t232 * t313;
t400 = t300 * t317;
t389 = qJD(1) * t320;
t388 = qJD(2) * t316;
t286 = t316 * t374;
t287 = t320 * t374;
t179 = t319 * t286 + t315 * t287 + t294 * t386 + t295 * t387;
t150 = -pkin(9) * t238 + t179;
t180 = -qJD(3) * t243 - t286 * t315 + t319 * t287;
t151 = -pkin(9) * t237 + t180;
t43 = qJD(4) * t469 + t150 * t318 + t151 * t314;
t137 = qJD(4) * t337 + t237 * t318 - t238 * t314;
t138 = qJD(4) * t232 + t237 * t314 + t318 * t238;
t219 = pkin(2) * t388 + pkin(3) * t238;
t50 = pkin(4) * t138 - pkin(10) * t137 + t219;
t381 = t153 * t382 + t313 * t50 + t317 * t43;
t303 = -pkin(4) - t439;
t373 = t232 * t382;
t372 = t413 / 0.2e1;
t371 = -t412 / 0.2e1;
t22 = -t84 * mrSges(7,1) + t83 * mrSges(7,2);
t366 = -t313 * t43 + t317 * t50;
t72 = t317 * t153 - t158 * t313;
t261 = -pkin(2) * t399 + t302 * t318;
t8 = -qJD(5) * t53 - t27 * t313 + t317 * t40;
t1 = pkin(5) * t117 - qJ(6) * t83 - qJD(6) * t197 + t8;
t360 = -t8 * mrSges(6,3) - t1 * mrSges(7,3);
t255 = -pkin(4) - t261;
t359 = -mrSges(6,3) * t52 - mrSges(7,3) * t29;
t358 = -mrSges(6,3) * t53 - t39 * mrSges(7,3);
t356 = mrSges(6,1) * t317 - mrSges(6,2) * t313;
t354 = mrSges(7,1) * t317 - mrSges(7,2) * t313;
t336 = -qJ(6) * t137 - qJD(6) * t232;
t3 = qJ(6) * t84 + qJD(6) * t196 + t7;
t327 = t8 * mrSges(6,1) + t1 * mrSges(7,1) - t7 * mrSges(6,2) - t3 * mrSges(7,2);
t44 = qJD(4) * t158 + t150 * t314 - t318 * t151;
t325 = m(6) * (-qJD(5) * t339 - t313 * t8 + t435);
t12 = -pkin(5) * t84 + t28;
t323 = t3 * t317 * mrSges(7,3) - t27 * mrSges(5,2) + mrSges(6,3) * t435 + Ifges(5,5) * t116 - Ifges(5,6) * t117 - t12 * t354 + t484 * t464 + t486 * t463 + t488 * t458 + t496 * t445 + t497 * t443 - t478 * t383 / 0.2e1 + t495 * t382 / 0.2e1 + (-t356 - mrSges(5,1)) * t28 + t489 * qJD(5) + (t333 + t332 + t331 + t330 + t329 + t328) * qJD(5) / 0.2e1;
t207 = t266 * Ifges(4,2) + t312 * Ifges(4,6) + t429;
t257 = Ifges(4,4) * t266;
t208 = t267 * Ifges(4,1) + t312 * Ifges(4,5) + t257;
t322 = -t293 * (mrSges(4,1) * t267 + mrSges(4,2) * t266) + t323 + t207 * t447 + (Ifges(4,1) * t266 - t429) * t448 + t228 * t432 - t171 * mrSges(4,2) + t172 * mrSges(4,1) + Ifges(4,5) * t224 - Ifges(4,6) * t225 - t312 * (Ifges(4,5) * t266 - Ifges(4,6) * t267) / 0.2e1 - (-Ifges(4,2) * t267 + t208 + t257) * t266 / 0.2e1 + t478 * t473 / 0.2e1 + (t29 * t404 + t39 * t473) * mrSges(7,3) + (t404 * t52 + t473 * t53) * mrSges(6,3) + (-t428 + t95 + t96) * t451 + (t482 + t498 * t455 + t499 * t454 + (Ifges(7,3) + Ifges(6,3)) * t452 + t480) * t338 + (Ifges(5,1) * t451 + t338 * t501 + t487 * t452 + t483 * t454 + t485 * t455 - t502) * t364;
t305 = Ifges(3,4) * t389;
t292 = t311 + t437;
t291 = t433 * t313;
t290 = mrSges(3,3) * t389 - t410;
t289 = -mrSges(3,3) * t390 + t411;
t288 = t303 - t440;
t276 = t311 + t400;
t275 = t393 * t313;
t265 = Ifges(3,1) * t390 + t305 + t413;
t264 = t412 + (t320 * Ifges(3,2) + t430) * qJD(1);
t247 = t255 - t439;
t246 = mrSges(4,1) * t312 - t431;
t245 = -mrSges(4,2) * t312 + t432;
t244 = t306 + t441;
t241 = t256 * t317 + t311;
t240 = t394 * t313;
t227 = -mrSges(4,1) * t266 + mrSges(4,2) * t267;
t200 = -mrSges(5,2) * t309 + mrSges(5,3) * t364;
t168 = -mrSges(5,1) * t364 + mrSges(5,2) * t338;
t145 = mrSges(6,1) * t213 - mrSges(6,3) * t197;
t144 = mrSges(7,1) * t213 - mrSges(7,3) * t197;
t143 = -mrSges(6,2) * t213 + mrSges(6,3) * t196;
t142 = -mrSges(7,2) * t213 + mrSges(7,3) * t196;
t132 = -mrSges(7,1) * t196 + mrSges(7,2) * t197;
t113 = Ifges(6,3) * t117;
t112 = Ifges(7,3) * t117;
t107 = pkin(5) * t402 - t469;
t88 = t120 + t204;
t82 = Ifges(6,5) * t83;
t81 = Ifges(7,5) * t83;
t80 = Ifges(6,6) * t84;
t79 = Ifges(7,6) * t84;
t54 = -qJ(6) * t402 + t73;
t47 = -pkin(5) * t337 - t232 * t311 + t72;
t33 = -mrSges(6,2) * t117 + mrSges(6,3) * t84;
t32 = -mrSges(7,2) * t117 + mrSges(7,3) * t84;
t31 = mrSges(6,1) * t117 - mrSges(6,3) * t83;
t30 = mrSges(7,1) * t117 - mrSges(7,3) * t83;
t21 = (t137 * t313 + t373) * pkin(5) + t44;
t10 = -qJD(5) * t73 + t366;
t9 = -t158 * t383 + t381;
t5 = -qJ(6) * t373 + (-qJD(5) * t158 + t336) * t313 + t381;
t4 = pkin(5) * t138 + t336 * t317 + (-t152 + (qJ(6) * t232 - t153) * t313) * qJD(5) + t366;
t2 = [(t509 + t466) * t137 + (t265 / 0.2e1 - pkin(7) * t289 + t372 + (-0.2e1 * t459 + 0.3e1 / 0.2e1 * Ifges(3,4) * t320) * qJD(1)) * t320 * qJD(2) + t293 * (mrSges(4,1) * t238 + mrSges(4,2) * t237) + t312 * (Ifges(4,5) * t237 - Ifges(4,6) * t238) / 0.2e1 + m(6) * (t10 * t52 + t108 * t44 + t53 * t9 + t7 * t73 + t72 * t8 - t422) + m(5) * (-t119 * t44 + t120 * t43 + t158 * t27 + t202 * t248 + t219 * t239 - t422) - (-t27 * mrSges(5,3) + t81 / 0.2e1 + t79 / 0.2e1 + t112 / 0.2e1 + t82 / 0.2e1 + t80 / 0.2e1 + t113 / 0.2e1 - Ifges(5,4) * t116 + t202 * mrSges(5,1) + t378 * t84 + t379 * t83 + (Ifges(5,2) - t377) * t117 + t327) * t337 + m(4) * (t171 * t243 + t172 * t242 + t179 * t229 + t180 * t228) + (-t116 * t469 - t117 * t158) * mrSges(5,3) - t469 * t23 + (t171 * t278 - t172 * t279 - t224 * t242 - t225 * t243 - t228 * t237 - t229 * t238) * mrSges(4,3) + (-t278 * t225 - t238 * t449) * Ifges(4,2) + (t278 * t224 - t279 * t225 + t237 * t449 - t238 * t447) * Ifges(4,4) + t304 * (mrSges(4,1) * t225 + mrSges(4,2) * t224) + (-t419 / 0.2e1 + t465) * t138 + m(7) * (t1 * t47 + t107 * t12 + t21 * t85 + t29 * t4 + t3 * t54 + t39 * t5) + t47 * t30 + t54 * t32 + (t202 * mrSges(5,2) + t12 * t353 + Ifges(5,1) * t116 - Ifges(5,4) * t117 + (mrSges(5,3) + t355) * t28 + (-t1 * t317 - t3 * t313) * mrSges(7,3) + (-t313 * t7 - t317 * t8) * mrSges(6,3) + (t85 * t354 + t108 * t356 + (t29 * t313 - t317 * t39) * mrSges(7,3) - t468 * mrSges(6,3) + t486 * t455 + t484 * t454 + t488 * t452 + t478 * t444) * qJD(5) + t483 * t464 + t485 * t463 + t487 * t458 - (t495 * qJD(5) + t497) * t313 / 0.2e1 + t496 * t443) * t232 + t72 * t31 + t73 * t33 + t107 * t22 + t21 * t132 + t5 * t142 + t9 * t143 + t4 * t144 + t10 * t145 + t392 * t44 + t43 * t200 + t219 * t168 + t237 * t208 / 0.2e1 - t238 * t207 / 0.2e1 + (t279 * t224 + t237 * t447) * Ifges(4,1) + t179 * t245 + t180 * t246 + t248 * (mrSges(5,1) * t117 + mrSges(5,2) * t116) + (-t264 / 0.2e1 - pkin(7) * t290 + t371 + (-0.2e1 * t460 - 0.3e1 / 0.2e1 * t430 + (-0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1)) * t320) * qJD(1) + (0.2e1 * t442 + t227 + qJD(1) * (-mrSges(4,1) * t278 + mrSges(4,2) * t279)) * pkin(2)) * t388; t322 + (-t222 * t145 - t256 * t31 + (-t143 * t256 + t358) * qJD(5) + t360) * t313 + (t222 * t143 + t256 * t33 + (-t145 * t256 + t359) * qJD(5)) * t317 + ((t245 * t319 - t246 * t315) * qJD(3) + (-t224 * t319 - t225 * t315) * mrSges(4,3)) * pkin(2) + (-t228 * t233 - t229 * t234 + (t171 * t315 + t172 * t319 + (-t228 * t315 + t229 * t319) * qJD(3)) * pkin(2)) * m(4) + t229 * t431 + t470 * t200 + (t119 * t471 + t120 * t470 - t239 * t244 - t261 * t28 + t262 * t27) * m(5) + (-t471 * t108 - t222 * t415 + t255 * t28 + t493 * t53 - t52 * t57) * m(6) - t392 * t471 + t472 * t132 + t256 * t325 + t475 * t142 + t477 * t144 + (t1 * t240 + t12 * t247 + t241 * t3 + t29 * t477 + t39 * t475 + t472 * t85) * m(7) - t58 * t143 - t57 * t145 + (-t116 * t261 - t117 * t262 + t409) * mrSges(5,3) + t240 * t30 + t241 * t32 - t244 * t168 - t234 * t245 - t233 * t246 + t247 * t22 + t255 * t23 + ((-t305 / 0.2e1 - t265 / 0.2e1 + t372 + qJD(1) * t459 + (t289 - t411) * pkin(7)) * t320 + (t264 / 0.2e1 + t371 + (t460 + t430 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t320) * qJD(1) + (t290 + t410) * pkin(7) + (-t227 - t442) * pkin(2)) * t316) * qJD(1); t322 + (-t300 * t31 + t360) * t313 + ((-t145 * t300 + t359) * t317 + (-t143 * t300 + t358) * t313) * qJD(5) + t490 * t132 + t492 * t142 - m(5) * (-t119 * t125 + t120 * t126) - m(6) * (t108 * t125 + t52 * t55 + t53 * t56) + t33 * t400 + mrSges(5,3) * t409 + (-t267 * t168 + (-t116 * t318 - t117 * t314) * mrSges(5,3) + (t392 * t314 + (t143 * t317 - t145 * t313 + t200) * t318 + m(6) * (t108 * t314 + t318 * t468)) * qJD(4) + (t27 * t314 - t28 * t318 + 0.2e1 * t239 * t448 + (-t119 * t314 + t120 * t318) * qJD(4)) * m(5)) * pkin(3) + t300 * t325 - t56 * t143 - t55 * t145 - t392 * t125 + t494 * t144 - t126 * t200 + (t246 + t431) * t229 - t228 * t245 + t275 * t30 + t276 * t32 + t288 * t22 + t491 * (-pkin(4) - t440) + (t1 * t275 + t12 * t288 + t276 * t3 + t29 * t494 + t39 * t492 + t490 * t85) * m(7); ((t501 - Ifges(5,1) / 0.2e1) * t338 - t466) * t364 + t323 + t476 * t144 + (-pkin(10) * t31 + t360) * t313 + t474 * t142 - m(6) * (t108 * t120 + t52 * t59 + t53 * t60) - t392 * t120 + pkin(10) * t325 + t33 * t437 + ((-pkin(10) * t145 + t359) * t317 + (pkin(5) * t132 - pkin(10) * t143 + t358) * t313) * qJD(5) - t88 * t132 - t60 * t143 - t59 * t145 - t119 * t200 - t465 * t338 + t291 * t30 + t292 * t32 + t303 * t22 - t491 * pkin(4) + (t1 * t291 + t12 * t303 + t292 * t3 + (t307 - t88) * t85 + t474 * t39 + t476 * t29) * m(7); t327 + t79 + t113 + t112 + t82 + t81 + t80 + (-t132 * t197 + t30) * pkin(5) + (-(-t29 + t38) * t39 + (-t197 * t85 + t1) * pkin(5)) * m(7) - t38 * t142 - t52 * t143 + t39 * t144 + t53 * t145 - t85 * (mrSges(7,1) * t197 + mrSges(7,2) * t196) - t108 * (mrSges(6,1) * t197 + mrSges(6,2) * t196) + (t196 * t52 + t197 * t53) * mrSges(6,3) + (t196 * t29 + t197 * t39) * mrSges(7,3) + (t508 * t196 - t503) * t454 + t478 * t453 + (t196 * t499 - t197 * t498) * t452 + (-t506 * t197 + t495 + t510) * t455; -t196 * t142 + t197 * t144 + 0.2e1 * (t12 / 0.2e1 + t39 * t455 + t29 * t453) * m(7) + t22;];
tauc  = t2(:);
