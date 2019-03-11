% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:37
% EndTime: 2019-03-09 10:56:49
% DurationCPUTime: 39.06s
% Computational Cost: add. (22115->879), mult. (59025->1249), div. (0->0), fcn. (48473->12), ass. (0->384)
t342 = sin(pkin(11));
t345 = cos(pkin(11));
t348 = sin(qJ(4));
t351 = cos(qJ(4));
t315 = t342 * t351 + t345 * t348;
t343 = sin(pkin(6));
t352 = cos(qJ(2));
t417 = t343 * t352;
t360 = t315 * t417;
t268 = qJD(1) * t360;
t307 = t315 * qJD(4);
t407 = t268 - t307;
t313 = t342 * t348 - t351 * t345;
t359 = t313 * t417;
t269 = qJD(1) * t359;
t306 = t313 * qJD(4);
t406 = -t269 + t306;
t349 = sin(qJ(2));
t372 = pkin(2) * t349 - qJ(3) * t352;
t405 = qJD(1) * t343;
t297 = t372 * t405;
t346 = cos(pkin(6));
t404 = qJD(1) * t346;
t397 = pkin(1) * t404;
t385 = t352 * t397;
t394 = t349 * t405;
t298 = -pkin(8) * t394 + t385;
t233 = t345 * t297 - t342 * t298;
t414 = t345 * t352;
t362 = (pkin(3) * t349 - pkin(9) * t414) * t343;
t201 = qJD(1) * t362 + t233;
t234 = t342 * t297 + t345 * t298;
t393 = t352 * t405;
t383 = t342 * t393;
t211 = -pkin(9) * t383 + t234;
t457 = pkin(9) + qJ(3);
t322 = t457 * t342;
t324 = t457 * t345;
t499 = -t351 * t322 - t324 * t348;
t508 = -t313 * qJD(3) + qJD(4) * t499 - t348 * t201 - t351 * t211;
t529 = qJ(5) * t394 - t508;
t299 = pkin(8) * t393 + t349 * t397;
t259 = pkin(3) * t383 + t299;
t528 = -t407 * pkin(4) + t406 * qJ(5) - qJD(5) * t315 - t259;
t341 = sin(pkin(12));
t344 = cos(pkin(12));
t328 = qJD(4) - t393;
t384 = qJD(2) + t404;
t272 = qJ(3) * t384 + t299;
t292 = (-pkin(2) * t352 - qJ(3) * t349 - pkin(1)) * t343;
t278 = qJD(1) * t292;
t202 = -t342 * t272 + t345 * t278;
t285 = t342 * t384 + t345 * t394;
t152 = -pkin(3) * t393 - t285 * pkin(9) + t202;
t275 = (qJD(2) * t372 - qJD(3) * t349) * t343;
t255 = qJD(1) * t275;
t403 = qJD(2) * t343;
t389 = qJD(1) * t403;
t381 = t349 * t389;
t289 = -pkin(8) * t381 + qJD(2) * t385;
t256 = qJD(3) * t384 + t289;
t194 = t345 * t255 - t342 * t256;
t358 = qJD(2) * t362;
t155 = qJD(1) * t358 + t194;
t203 = t345 * t272 + t342 * t278;
t354 = t342 * t394 - t345 * t384;
t172 = -pkin(9) * t354 + t203;
t195 = t342 * t255 + t345 * t256;
t380 = t352 * t389;
t365 = t342 * t380;
t173 = -pkin(9) * t365 + t195;
t400 = qJD(4) * t351;
t401 = qJD(4) * t348;
t38 = t152 * t400 + t348 * t155 - t172 * t401 + t351 * t173;
t35 = qJ(5) * t381 + qJD(5) * t328 + t38;
t221 = t348 * t285 + t351 * t354;
t356 = qJD(2) * t359;
t170 = -qJD(1) * t356 - qJD(4) * t221;
t353 = t351 * t285 - t348 * t354;
t357 = qJD(2) * t360;
t171 = qJD(1) * t357 + qJD(4) * t353;
t309 = t346 * t349 * pkin(1) + pkin(8) * t417;
t301 = t309 * qJD(2);
t290 = qJD(1) * t301;
t249 = pkin(3) * t365 + t290;
t72 = t171 * pkin(4) - t170 * qJ(5) - qJD(5) * t353 + t249;
t18 = -t341 * t35 + t344 * t72;
t19 = t341 * t72 + t344 * t35;
t143 = -t170 * t341 + t344 * t381;
t13 = pkin(10) * t143 + t19;
t347 = sin(qJ(6));
t350 = cos(qJ(6));
t190 = t328 * t341 + t344 * t353;
t261 = -pkin(2) * t384 + qJD(3) - t298;
t218 = pkin(3) * t354 + t261;
t102 = t221 * pkin(4) - qJ(5) * t353 + t218;
t97 = t152 * t348 + t172 * t351;
t90 = qJ(5) * t328 + t97;
t40 = t344 * t102 - t341 * t90;
t25 = pkin(5) * t221 - pkin(10) * t190 + t40;
t388 = t344 * t328 - t341 * t353;
t41 = t341 * t102 + t344 * t90;
t31 = pkin(10) * t388 + t41;
t5 = t25 * t350 - t31 * t347;
t144 = t170 * t344 + t341 * t381;
t7 = pkin(5) * t171 - pkin(10) * t144 + t18;
t1 = qJD(6) * t5 + t13 * t350 + t347 * t7;
t6 = t25 * t347 + t31 * t350;
t2 = -qJD(6) * t6 - t13 * t347 + t350 * t7;
t379 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t430 = t170 * Ifges(5,4);
t485 = t171 / 0.2e1;
t486 = t144 / 0.2e1;
t487 = t143 / 0.2e1;
t112 = t190 * t350 + t347 * t388;
t46 = -qJD(6) * t112 + t143 * t350 - t144 * t347;
t43 = Ifges(7,6) * t46;
t506 = -t190 * t347 + t350 * t388;
t45 = qJD(6) * t506 + t143 * t347 + t144 * t350;
t44 = Ifges(7,5) * t45;
t8 = Ifges(7,3) * t171 + t43 + t44;
t527 = t379 + t18 * mrSges(6,1) + t249 * mrSges(5,1) + Ifges(6,3) * t485 + t8 / 0.2e1 - t19 * mrSges(6,2) + 0.2e1 * Ifges(6,5) * t486 + 0.2e1 * Ifges(6,6) * t487 - t430 / 0.2e1;
t235 = t269 * t341 + t344 * t394;
t420 = t341 * t306;
t387 = t235 - t420;
t236 = -t269 * t344 + t341 * t394;
t416 = t344 * t306;
t526 = t236 + t416;
t513 = t341 * t529 + t528 * t344;
t512 = t528 * t341 - t344 * t529;
t525 = -pkin(5) * t407 + pkin(10) * t526 + t513;
t524 = -pkin(10) * t387 + t512;
t267 = -t322 * t348 + t324 * t351;
t509 = -qJD(3) * t315 - qJD(4) * t267 - t201 * t351 + t348 * t211;
t521 = t171 * Ifges(5,2);
t520 = t249 * mrSges(5,2);
t510 = pkin(4) * t394 - t509;
t216 = Ifges(5,4) * t221;
t445 = Ifges(5,5) * t328;
t128 = Ifges(5,1) * t353 - t216 + t445;
t96 = t152 * t351 - t348 * t172;
t519 = t218 * mrSges(5,2) - t96 * mrSges(5,3) + t128 / 0.2e1;
t438 = Ifges(5,6) * t328;
t452 = Ifges(5,4) * t353;
t127 = -Ifges(5,2) * t221 + t438 + t452;
t517 = t41 * mrSges(6,2) + t6 * mrSges(7,2) + t97 * mrSges(5,3) + t127 / 0.2e1 - t218 * mrSges(5,1) - t40 * mrSges(6,1) - t5 * mrSges(7,1);
t340 = -pkin(3) * t345 - pkin(2);
t246 = pkin(4) * t313 - qJ(5) * t315 + t340;
t180 = t344 * t246 - t267 * t341;
t462 = pkin(10) * t344;
t141 = pkin(5) * t313 - t315 * t462 + t180;
t181 = t341 * t246 + t344 * t267;
t422 = t315 * t341;
t156 = -pkin(10) * t422 + t181;
t83 = t141 * t350 - t156 * t347;
t516 = qJD(6) * t83 + t347 * t525 + t350 * t524;
t84 = t141 * t347 + t156 * t350;
t515 = -qJD(6) * t84 - t347 * t524 + t350 * t525;
t217 = qJD(6) + t221;
t434 = Ifges(7,3) * t217;
t435 = Ifges(7,6) * t506;
t442 = Ifges(7,5) * t112;
t48 = t434 + t435 + t442;
t437 = Ifges(6,6) * t388;
t444 = Ifges(6,5) * t190;
t86 = Ifges(6,3) * t221 + t437 + t444;
t514 = t86 + t48;
t511 = pkin(5) * t387 + t510;
t314 = t341 * t350 + t344 * t347;
t129 = t314 * t221;
t305 = t314 * qJD(6);
t413 = t129 + t305;
t366 = t341 * t347 - t344 * t350;
t131 = t366 * t221;
t304 = t366 * qJD(6);
t412 = t131 + t304;
t426 = t328 * Ifges(5,3);
t427 = t353 * Ifges(5,5);
t507 = t426 + t427;
t495 = t45 / 0.2e1;
t494 = t46 / 0.2e1;
t456 = pkin(10) + qJ(5);
t321 = t456 * t341;
t323 = t456 * t344;
t264 = -t321 * t350 - t323 * t347;
t137 = pkin(4) * t353 + qJ(5) * t221;
t58 = t344 * t137 - t341 * t96;
t32 = pkin(5) * t353 + t221 * t462 + t58;
t423 = t221 * t341;
t59 = t341 * t137 + t344 * t96;
t47 = pkin(10) * t423 + t59;
t505 = -qJD(5) * t366 + qJD(6) * t264 - t32 * t347 - t350 * t47;
t266 = -t321 * t347 + t323 * t350;
t504 = -qJD(5) * t314 - qJD(6) * t266 - t32 * t350 + t347 * t47;
t498 = -0.2e1 * pkin(1);
t497 = Ifges(7,4) * t495 + Ifges(7,2) * t494 + Ifges(7,6) * t485;
t496 = Ifges(7,1) * t495 + Ifges(7,4) * t494 + Ifges(7,5) * t485;
t65 = Ifges(6,1) * t144 + Ifges(6,4) * t143 + Ifges(6,5) * t171;
t493 = t65 / 0.2e1;
t492 = -t506 / 0.2e1;
t491 = t506 / 0.2e1;
t490 = -t112 / 0.2e1;
t489 = t112 / 0.2e1;
t483 = t388 / 0.2e1;
t482 = t190 / 0.2e1;
t481 = -t217 / 0.2e1;
t480 = t217 / 0.2e1;
t479 = -t221 / 0.2e1;
t478 = t221 / 0.2e1;
t477 = t353 / 0.2e1;
t415 = t345 * t346;
t418 = t343 * t349;
t302 = -t342 * t418 + t415;
t303 = t342 * t346 + t345 * t418;
t367 = t351 * t302 - t303 * t348;
t476 = t367 / 0.2e1;
t238 = t302 * t348 + t303 * t351;
t474 = t238 / 0.2e1;
t473 = t302 / 0.2e1;
t472 = t303 / 0.2e1;
t471 = t328 / 0.2e1;
t470 = -t341 / 0.2e1;
t469 = -t342 / 0.2e1;
t468 = t344 / 0.2e1;
t467 = t345 / 0.2e1;
t465 = t349 / 0.2e1;
t463 = pkin(1) * t352;
t461 = t38 * mrSges(5,2);
t39 = -t152 * t401 + t155 * t351 - t172 * t400 - t348 * t173;
t460 = t39 * mrSges(5,1);
t459 = t96 * mrSges(5,1);
t458 = t97 * mrSges(5,2);
t402 = qJD(2) * t349;
t392 = t343 * t402;
t399 = t346 * t463;
t300 = -pkin(8) * t392 + qJD(2) * t399;
t284 = qJD(3) * t346 + t300;
t209 = t345 * t275 - t342 * t284;
t177 = t209 + t358;
t291 = qJ(3) * t346 + t309;
t226 = -t342 * t291 + t345 * t292;
t179 = -pkin(3) * t417 - t303 * pkin(9) + t226;
t210 = t342 * t275 + t345 * t284;
t419 = t342 * t352;
t382 = t403 * t419;
t193 = -pkin(9) * t382 + t210;
t227 = t345 * t291 + t342 * t292;
t198 = pkin(9) * t302 + t227;
t56 = t348 * t177 + t179 * t400 + t351 * t193 - t198 * t401;
t53 = (qJ(5) * t402 - qJD(5) * t352) * t343 + t56;
t185 = qJD(4) * t367 - t356;
t186 = qJD(4) * t238 + t357;
t260 = pkin(3) * t382 + t301;
t78 = t186 * pkin(4) - t185 * qJ(5) - t238 * qJD(5) + t260;
t24 = t341 * t78 + t344 * t53;
t455 = mrSges(4,2) * t345;
t454 = Ifges(4,1) * t345;
t453 = Ifges(3,4) * t349;
t451 = Ifges(6,4) * t341;
t450 = Ifges(6,4) * t344;
t449 = Ifges(7,4) * t112;
t448 = Ifges(3,5) * t352;
t447 = Ifges(4,5) * t285;
t446 = Ifges(4,5) * t349;
t443 = Ifges(6,5) * t344;
t441 = Ifges(3,6) * t346;
t440 = Ifges(3,6) * t349;
t439 = Ifges(4,6) * t349;
t436 = Ifges(6,6) * t341;
t431 = t170 * Ifges(5,1);
t429 = t171 * Ifges(5,4);
t428 = t221 * Ifges(5,6);
t425 = t342 * Ifges(4,6);
t424 = t352 * Ifges(3,2);
t108 = t348 * t179 + t351 * t198;
t104 = -qJ(5) * t417 + t108;
t334 = pkin(8) * t418;
t294 = t334 + (-pkin(2) - t463) * t346;
t243 = -t302 * pkin(3) + t294;
t125 = -pkin(4) * t367 - t238 * qJ(5) + t243;
t67 = t344 * t104 + t341 * t125;
t150 = t235 * t350 - t236 * t347;
t160 = t304 * t315 + t306 * t314;
t411 = t150 - t160;
t151 = t235 * t347 + t236 * t350;
t159 = -t305 * t315 + t306 * t366;
t410 = t151 - t159;
t114 = -mrSges(6,1) * t388 + mrSges(6,2) * t190;
t200 = mrSges(5,1) * t328 - mrSges(5,3) * t353;
t409 = t200 - t114;
t408 = -mrSges(3,1) * t384 + mrSges(4,1) * t354 + t285 * mrSges(4,2) + mrSges(3,3) * t394;
t273 = mrSges(4,1) * t365 + t380 * t455;
t398 = Ifges(4,4) * t419;
t396 = -Ifges(5,2) / 0.2e1 - Ifges(6,3) / 0.2e1;
t395 = Ifges(5,5) * t170 - Ifges(5,6) * t171 + Ifges(5,3) * t381;
t391 = t419 / 0.2e1;
t16 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t23 = -t341 * t53 + t344 * t78;
t99 = t171 * mrSges(5,1) + t170 * mrSges(5,2);
t80 = -t143 * mrSges(6,1) + t144 * mrSges(6,2);
t66 = -t104 * t341 + t344 * t125;
t107 = t179 * t351 - t348 * t198;
t105 = pkin(4) * t417 - t107;
t378 = mrSges(6,1) * t341 + mrSges(6,2) * t344;
t377 = -Ifges(4,4) * t342 + t454;
t376 = Ifges(6,1) * t344 - t451;
t375 = Ifges(4,4) * t345 - Ifges(4,2) * t342;
t374 = -Ifges(6,2) * t341 + t450;
t373 = -t436 + t443;
t371 = t18 * t344 + t19 * t341;
t370 = -t18 * t341 + t19 * t344;
t369 = -t341 * t40 + t344 * t41;
t215 = t344 * t238 - t341 * t417;
t36 = -pkin(5) * t367 - pkin(10) * t215 + t66;
t214 = -t341 * t238 - t344 * t417;
t51 = pkin(10) * t214 + t67;
t14 = -t347 * t51 + t350 * t36;
t15 = t347 * t36 + t350 * t51;
t123 = -mrSges(6,2) * t221 + mrSges(6,3) * t388;
t124 = mrSges(6,1) * t221 - mrSges(6,3) * t190;
t368 = t123 * t344 - t124 * t341;
t133 = t214 * t350 - t215 * t347;
t134 = t214 * t347 + t215 * t350;
t57 = t177 * t351 - t179 * t401 - t348 * t193 - t198 * t400;
t364 = mrSges(4,1) * t349 - mrSges(4,3) * t414;
t363 = -mrSges(4,2) * t349 - mrSges(4,3) * t419;
t89 = -pkin(4) * t328 + qJD(5) - t96;
t361 = t375 * t467;
t355 = Ifges(4,4) * t414 - Ifges(4,2) * t419 + t439;
t54 = -pkin(4) * t392 - t57;
t37 = -pkin(4) * t381 - t39;
t339 = -pkin(5) * t344 - pkin(4);
t329 = Ifges(3,4) * t393;
t326 = Ifges(3,5) * t380;
t308 = -t334 + t399;
t296 = -mrSges(3,2) * t384 + mrSges(3,3) * t393;
t280 = t364 * t389;
t279 = t363 * t389;
t271 = Ifges(3,1) * t394 + Ifges(3,5) * t384 + t329;
t270 = Ifges(3,6) * qJD(2) + (t441 + (t424 + t453) * t343) * qJD(1);
t252 = -mrSges(4,1) * t393 - t285 * mrSges(4,3);
t251 = mrSges(4,2) * t393 - mrSges(4,3) * t354;
t245 = (t352 * t377 + t446) * t389;
t244 = (t352 * t375 + t439) * t389;
t240 = t366 * t315;
t239 = t314 * t315;
t232 = pkin(5) * t422 - t499;
t208 = Ifges(4,1) * t285 - Ifges(4,4) * t354 - Ifges(4,5) * t393;
t207 = Ifges(4,4) * t285 - Ifges(4,2) * t354 - Ifges(4,6) * t393;
t206 = -Ifges(4,6) * t354 - Ifges(4,3) * t393 + t447;
t199 = -mrSges(5,2) * t328 - mrSges(5,3) * t221;
t167 = t185 * t344 + t341 * t392;
t166 = -t185 * t341 + t344 * t392;
t146 = -mrSges(5,2) * t381 - mrSges(5,3) * t171;
t145 = mrSges(5,1) * t381 - mrSges(5,3) * t170;
t138 = mrSges(5,1) * t221 + mrSges(5,2) * t353;
t126 = -t428 + t507;
t106 = Ifges(7,4) * t506;
t95 = Ifges(5,5) * t381 - t429 + t431;
t94 = Ifges(5,6) * t381 + t430 - t521;
t92 = mrSges(6,1) * t171 - mrSges(6,3) * t144;
t91 = -mrSges(6,2) * t171 + mrSges(6,3) * t143;
t88 = t190 * Ifges(6,1) + Ifges(6,4) * t388 + Ifges(6,5) * t221;
t87 = t190 * Ifges(6,4) + Ifges(6,2) * t388 + Ifges(6,6) * t221;
t82 = mrSges(7,1) * t217 - mrSges(7,3) * t112;
t81 = -mrSges(7,2) * t217 + mrSges(7,3) * t506;
t79 = -pkin(5) * t214 + t105;
t73 = -pkin(5) * t423 + t97;
t69 = -pkin(5) * t388 + t89;
t64 = t144 * Ifges(6,4) + t143 * Ifges(6,2) + t171 * Ifges(6,6);
t62 = -qJD(6) * t134 + t166 * t350 - t167 * t347;
t61 = qJD(6) * t133 + t166 * t347 + t167 * t350;
t60 = -mrSges(7,1) * t506 + mrSges(7,2) * t112;
t50 = Ifges(7,1) * t112 + Ifges(7,5) * t217 + t106;
t49 = Ifges(7,2) * t506 + Ifges(7,6) * t217 + t449;
t33 = -pkin(5) * t166 + t54;
t30 = -mrSges(7,2) * t171 + mrSges(7,3) * t46;
t29 = mrSges(7,1) * t171 - mrSges(7,3) * t45;
t28 = -pkin(5) * t143 + t37;
t20 = pkin(10) * t166 + t24;
t17 = pkin(5) * t186 - pkin(10) * t167 + t23;
t4 = -qJD(6) * t15 + t17 * t350 - t20 * t347;
t3 = qJD(6) * t14 + t17 * t347 + t20 * t350;
t9 = [(Ifges(6,4) * t215 + Ifges(6,2) * t214) * t487 + (Ifges(6,4) * t167 + Ifges(6,2) * t166) * t483 + (t1 * t133 - t134 * t2 - t5 * t61 + t6 * t62) * mrSges(7,3) + (Ifges(7,5) * t61 + Ifges(7,6) * t62) * t480 + (Ifges(6,3) * t478 - t517 + t514 / 0.2e1 - Ifges(5,6) * t471 - Ifges(5,4) * t477 - Ifges(5,2) * t479 + Ifges(6,5) * t482 + Ifges(6,6) * t483 + Ifges(7,3) * t480 + Ifges(7,5) * t489 + Ifges(7,6) * t491) * t186 + t238 * t520 + (Ifges(5,1) * t477 + Ifges(5,4) * t479 + Ifges(5,5) * t471 + t519) * t185 + (Ifges(7,4) * t134 + Ifges(7,2) * t133) * t494 + (Ifges(7,4) * t61 + Ifges(7,2) * t62) * t491 + (Ifges(7,1) * t134 + Ifges(7,4) * t133) * t495 + (Ifges(7,1) * t61 + Ifges(7,4) * t62) * t489 + (t166 * t41 - t167 * t40 - t18 * t215 + t19 * t214) * mrSges(6,3) - t171 * (Ifges(5,4) * t238 - Ifges(5,6) * t417) / 0.2e1 + t170 * (Ifges(5,1) * t238 - Ifges(5,5) * t417) / 0.2e1 + (Ifges(6,1) * t215 + Ifges(6,4) * t214) * t486 + (Ifges(6,1) * t167 + Ifges(6,4) * t166) * t482 + (Ifges(6,5) * t167 + Ifges(6,6) * t166) * t478 + m(7) * (t1 * t15 + t14 * t2 + t28 * t79 + t3 * t6 + t33 * t69 + t4 * t5) + m(6) * (t105 * t37 + t18 * t66 + t19 * t67 + t23 * t40 + t24 * t41 + t54 * t89) + m(5) * (t107 * t39 + t108 * t38 + t218 * t260 + t243 * t249 + t56 * t97 + t57 * t96) + m(4) * (t194 * t226 + t195 * t227 + t202 * t209 + t203 * t210 + t261 * t301 + t290 * t294) + m(3) * (t289 * t309 - t290 * t308 - t298 * t301 + t299 * t300) + t346 * t326 / 0.2e1 + ((-t309 * mrSges(3,3) + Ifges(4,5) * t472 + Ifges(5,5) * t474 + Ifges(5,6) * t476 - 0.3e1 / 0.2e1 * t441 + (t473 + t415 / 0.2e1) * Ifges(4,6) + (mrSges(3,1) * t498 + (-0.3e1 / 0.2e1 * Ifges(3,4) - t425 / 0.2e1) * t349) * t343) * t349 + ((Ifges(4,1) * t303 + Ifges(4,4) * t302) * t467 - t308 * mrSges(3,3) + (Ifges(4,4) * t303 + Ifges(4,2) * t302) * t469 + (Ifges(3,5) + t361) * t346 + (mrSges(3,2) * t498 + (-0.3e1 / 0.2e1 * t345 * Ifges(4,5) + 0.3e1 / 0.2e1 * t425 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t352) * t343 + (-0.3e1 / 0.2e1 * Ifges(4,3) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(3,1) - Ifges(5,3) / 0.2e1 + t375 * t469) * t418) * t352) * t389 + t214 * t64 / 0.2e1 + t37 * (-mrSges(6,1) * t214 + mrSges(6,2) * t215) + t57 * t200 + ((t271 / 0.2e1 - t298 * mrSges(3,3) + t207 * t469 + t208 * t467 + t261 * (mrSges(4,1) * t342 + t455) + t285 * t377 / 0.2e1 + (Ifges(3,5) / 0.2e1 + t361) * qJD(2) + (-t202 * t345 - t203 * t342) * mrSges(4,3)) * t352 + (-t270 / 0.2e1 + t206 / 0.2e1 + t126 / 0.2e1 - t299 * mrSges(3,3) - t428 / 0.2e1 + t427 / 0.2e1 + t459 - t458 + t426 / 0.2e1 - t203 * mrSges(4,2) + t202 * mrSges(4,1) + t447 / 0.2e1 + (-Ifges(3,6) / 0.2e1 + Ifges(4,6) * t467) * qJD(2)) * t349) * t403 + t56 * t199 + (Ifges(6,5) * t215 + Ifges(7,5) * t134 + Ifges(6,6) * t214 + Ifges(7,6) * t133) * t485 + t417 * t461 - t238 * t39 * mrSges(5,3) - t417 * t460 + t166 * t87 / 0.2e1 + t89 * (-mrSges(6,1) * t166 + mrSges(6,2) * t167) + t167 * t88 / 0.2e1 + t108 * t146 + t107 * t145 + t28 * (-mrSges(7,1) * t133 + mrSges(7,2) * t134) + t24 * t123 + t23 * t124 + t54 * t114 + t105 * t80 + t66 * t92 - t395 * t417 / 0.2e1 + t194 * (-mrSges(4,1) * t417 - t303 * mrSges(4,3)) + t289 * (-t346 * mrSges(3,2) + mrSges(3,3) * t417) + t195 * (mrSges(4,2) * t417 + t302 * mrSges(4,3)) + (-mrSges(3,1) * t346 - mrSges(4,1) * t302 + mrSges(4,2) * t303 + mrSges(3,3) * t418) * t290 + t67 * t91 + t4 * t82 + t79 * t16 + t3 * t81 + t69 * (-mrSges(7,1) * t62 + mrSges(7,2) * t61) + t408 * t301 + t33 * t60 + t61 * t50 / 0.2e1 + t62 * t49 / 0.2e1 + t14 * t29 + t15 * t30 + t245 * t472 + t244 * t473 + t95 * t474 + t94 * t476 + t243 * t99 + t210 * t251 + t209 * t252 + t260 * t138 + t227 * t279 + t226 * t280 + ((-Ifges(6,3) - Ifges(7,3)) * t485 + t38 * mrSges(5,3) - t521 / 0.2e1 - Ifges(7,6) * t494 - Ifges(7,5) * t495 - t527) * t367 + t294 * t273 + t300 * t296 + t215 * t493 + t134 * t496 + t133 * t497; (-Ifges(5,4) * t306 + Ifges(6,5) * t236 - Ifges(5,2) * t307 + Ifges(6,6) * t235 + Ifges(6,3) * t268) * t479 + (t65 * t468 + t95 / 0.2e1 - t429 / 0.2e1 + t431 / 0.2e1 + t520 + t37 * t378 + t373 * t485 + t374 * t487 + t376 * t486 + t64 * t470 - t371 * mrSges(6,3)) * t315 + (t239 * t28 - t407 * t5 + t411 * t69) * mrSges(7,1) + (-t1 * t239 + t2 * t240 + t410 * t5 - t411 * t6) * mrSges(7,3) + (-Ifges(7,5) * t240 - Ifges(7,6) * t239) * t485 + (-Ifges(7,4) * t240 - Ifges(7,2) * t239) * t494 + (-Ifges(7,1) * t240 - Ifges(7,4) * t239) * t495 + (-t240 * t28 + t407 * t6 - t410 * t69) * mrSges(7,2) + (-pkin(2) * t290 + (-t202 * t342 + t203 * t345) * qJD(3) + (-t194 * t342 + t195 * t345) * qJ(3) - t202 * t233 - t203 * t234 - t261 * t299) * m(4) - t388 * (Ifges(6,4) * t236 + Ifges(6,2) * t235 + Ifges(6,6) * t268) / 0.2e1 + t326 + (-Ifges(5,5) * t306 - Ifges(5,6) * t307) * t471 + (-Ifges(5,1) * t306 - Ifges(5,4) * t307) * t477 + (Ifges(6,5) * t307 - t306 * t376) * t482 + (Ifges(6,6) * t307 - t306 * t374) * t483 + (t160 / 0.2e1 - t150 / 0.2e1) * t49 + (-t416 / 0.2e1 - t236 / 0.2e1) * t88 + (-t235 / 0.2e1 + t420 / 0.2e1) * t87 + (t195 * mrSges(4,3) + qJD(3) * t251 + qJ(3) * t279 + t244 / 0.2e1 - t290 * mrSges(4,1)) * t345 + (t245 / 0.2e1 + t290 * mrSges(4,2) - t194 * mrSges(4,3) - qJD(3) * t252 - qJ(3) * t280) * t342 + (t159 / 0.2e1 - t151 / 0.2e1) * t50 - (t80 - t145) * t499 + (-mrSges(6,1) * t407 + mrSges(6,3) * t526) * t40 + (mrSges(6,1) * t387 - mrSges(6,2) * t526) * t89 - t353 * (-Ifges(5,1) * t269 - Ifges(5,4) * t268) / 0.2e1 + (t269 / 0.2e1 - t306 / 0.2e1) * t128 - t328 * (-Ifges(5,5) * t269 - Ifges(5,6) * t268) / 0.2e1 + t180 * t92 + t181 * t91 + t83 * t29 + t84 * t30 + (-Ifges(5,4) * t269 - Ifges(5,2) * t268 + Ifges(6,3) * t307 - t306 * t373) * t478 + (-t313 * t38 - t315 * t39 + t406 * t96 + t407 * t97) * mrSges(5,3) + (mrSges(6,2) * t407 - mrSges(6,3) * t387) * t41 + (-mrSges(5,1) * t407 - mrSges(5,2) * t406) * t218 - t408 * t299 + (Ifges(7,5) * t159 + Ifges(7,6) * t160 + Ifges(7,3) * t307) * t480 + (Ifges(7,5) * t151 + Ifges(7,6) * t150 + Ifges(7,3) * t268) * t481 + t232 * t16 - t234 * t251 - t233 * t252 - t259 * t138 + (t270 * t465 - t261 * (mrSges(4,1) * t419 + mrSges(4,2) * t414) - t203 * t363 - t202 * t364 - t285 * (Ifges(4,1) * t414 - t398 + t446) / 0.2e1 + (-t346 * (-t440 + t448) / 0.2e1 + (t352 * (Ifges(4,5) * t414 - Ifges(4,6) * t419 + Ifges(4,3) * t349) / 0.2e1 + pkin(1) * (mrSges(3,1) * t349 + mrSges(3,2) * t352) + (t342 * t355 + t424) * t465) * t343 - t355 * t415 / 0.2e1) * qJD(1) + t207 * t391 - t208 * t414 / 0.2e1 + t428 * t465 - t349 * t459 + t349 * t458 + (t298 * t352 + t299 * t349) * mrSges(3,3) + (-t448 / 0.2e1 - t440 / 0.2e1 + (Ifges(5,5) * t315 - Ifges(5,6) * t313) * t465 + (t446 / 0.2e1 - t398 / 0.2e1) * t342 + t391 * t454) * qJD(2) - (t271 + t329) * t352 / 0.2e1 - (t206 + t126 + (Ifges(3,1) * t352 - t453) * t405 + t507) * t349 / 0.2e1) * t405 + t267 * t146 - t190 * (Ifges(6,1) * t236 + Ifges(6,4) * t235 + Ifges(6,5) * t268) / 0.2e1 - pkin(2) * t273 + (t44 / 0.2e1 + t43 / 0.2e1 - t94 / 0.2e1 + (Ifges(7,3) / 0.2e1 - t396) * t171 + t527) * t313 - t289 * mrSges(3,2) - t290 * mrSges(3,1) - t298 * t296 + t508 * t199 + t509 * t200 + (-t218 * t259 + t249 * t340 + t267 * t38 + t499 * t39 + t508 * t97 + t509 * t96) * m(5) + t510 * t114 + (Ifges(7,1) * t159 + Ifges(7,4) * t160 + Ifges(7,5) * t307) * t489 + (Ifges(7,1) * t151 + Ifges(7,4) * t150 + Ifges(7,5) * t268) * t490 + (Ifges(7,4) * t159 + Ifges(7,2) * t160 + Ifges(7,6) * t307) * t491 + (Ifges(7,4) * t151 + Ifges(7,2) * t150 + Ifges(7,6) * t268) * t492 - t240 * t496 - t239 * t497 + t340 * t99 + t511 * t60 + t512 * t123 + t513 * t124 + (t18 * t180 + t181 * t19 - t499 * t37 + t40 * t513 + t41 * t512 + t510 * t89) * m(6) + (t127 - t514) * (t268 / 0.2e1 - t307 / 0.2e1) + t515 * t82 + t516 * t81 + (t1 * t84 + t2 * t83 + t232 * t28 + t515 * t5 + t511 * t69 + t516 * t6) * m(7); -t413 * t82 - t412 * t81 + t354 * t251 + (-t60 + t409) * t353 + t99 + t285 * t252 - t366 * t29 + t314 * t30 + t341 * t91 + t344 * t92 - (-t199 - t368) * t221 + t273 + (t1 * t314 - t2 * t366 - t353 * t69 - t412 * t6 - t413 * t5) * m(7) + (t221 * t369 - t353 * t89 + t371) * m(6) + (t221 * t97 + t353 * t96 + t249) * m(5) + (t202 * t285 + t203 * t354 + t290) * m(4); (t438 / 0.2e1 - t444 / 0.2e1 - t437 / 0.2e1 - t434 / 0.2e1 - t435 / 0.2e1 - t442 / 0.2e1 - t86 / 0.2e1 - t48 / 0.2e1 + t452 / 0.2e1 + t517) * t353 + t504 * t82 + (t1 * t266 + t2 * t264 + t28 * t339 + t5 * t504 + t505 * t6 - t69 * t73) * m(7) + t505 * t81 + (-t216 / 0.2e1 + t445 / 0.2e1 + t89 * t378 + t376 * t482 + t374 * t483 + t87 * t470 + t88 * t468 + (t443 / 0.2e1 - t436 / 0.2e1) * t221 + (-t341 * t41 - t344 * t40) * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + t396) * t353 + t519) * t221 + (-t131 / 0.2e1 - t304 / 0.2e1) * t50 + t395 + t460 - t461 + (-t341 * t92 + t344 * t91) * qJ(5) + (-Ifges(7,5) * t304 - Ifges(7,6) * t305) * t480 + (-Ifges(7,1) * t304 - Ifges(7,4) * t305) * t489 + (-Ifges(7,4) * t304 - Ifges(7,2) * t305) * t491 + (-t129 / 0.2e1 - t305 / 0.2e1) * t49 - t96 * t199 + t64 * t468 + (-t1 * t366 - t2 * t314 + t412 * t5 - t413 * t6) * mrSges(7,3) + (Ifges(6,5) * t341 + Ifges(7,5) * t314 + Ifges(6,6) * t344 - Ifges(7,6) * t366) * t485 + (Ifges(7,4) * t314 - Ifges(7,2) * t366) * t494 + (Ifges(7,1) * t314 - Ifges(7,4) * t366) * t495 + t28 * (mrSges(7,1) * t366 + mrSges(7,2) * t314) - t366 * t497 - t58 * t124 - t59 * t123 - pkin(4) * t80 + (mrSges(7,1) * t413 - mrSges(7,2) * t412) * t69 - t73 * t60 + t409 * t97 + (Ifges(7,5) * t131 + Ifges(7,6) * t129) * t481 + (-pkin(4) * t37 + qJ(5) * t370 + qJD(5) * t369 - t40 * t58 - t41 * t59 - t89 * t97) * m(6) + t264 * t29 + t266 * t30 + (Ifges(6,1) * t341 + t450) * t486 + (Ifges(6,2) * t344 + t451) * t487 + (Ifges(7,1) * t131 + Ifges(7,4) * t129) * t490 + (Ifges(7,4) * t131 + Ifges(7,2) * t129) * t492 + t341 * t493 + t314 * t496 + t339 * t16 + t37 * (-mrSges(6,1) * t344 + mrSges(6,2) * t341) + t368 * qJD(5) + t370 * mrSges(6,3); t112 * t82 - t506 * t81 - t388 * t123 + t190 * t124 + t16 + t80 + (t112 * t5 - t506 * t6 + t28) * m(7) + (t190 * t40 - t388 * t41 + t37) * m(6); -t69 * (mrSges(7,1) * t112 + mrSges(7,2) * t506) + (Ifges(7,1) * t506 - t449) * t490 + t49 * t489 + (Ifges(7,5) * t506 - Ifges(7,6) * t112) * t481 - t5 * t81 + t6 * t82 + (t112 * t6 + t5 * t506) * mrSges(7,3) + t379 + t8 + (-Ifges(7,2) * t112 + t106 + t50) * t492;];
tauc  = t9(:);
