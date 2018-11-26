% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4]';
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
% Datum: 2018-11-23 18:10
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:09:51
% EndTime: 2018-11-23 18:10:14
% DurationCPUTime: 24.08s
% Computational Cost: add. (12118->801), mult. (32042->995), div. (0->0), fcn. (23896->8), ass. (0->343)
t491 = Ifges(6,5) + Ifges(7,4);
t492 = Ifges(5,5) + Ifges(6,4);
t470 = -Ifges(7,5) + t492;
t489 = Ifges(7,2) + Ifges(6,3);
t487 = Ifges(6,6) - Ifges(7,6);
t469 = -Ifges(5,4) + t491;
t471 = Ifges(7,1) + Ifges(6,1) + Ifges(5,1);
t525 = Ifges(4,5) / 0.2e1;
t296 = cos(qJ(2));
t289 = sin(pkin(6));
t383 = qJD(1) * t289;
t363 = t296 * t383;
t325 = -qJD(3) + t363;
t524 = -t325 / 0.2e1;
t490 = Ifges(6,2) + Ifges(5,3);
t515 = Ifges(5,6) - Ifges(6,6);
t295 = cos(qJ(3));
t292 = sin(qJ(3));
t359 = t292 * t363;
t379 = qJD(3) * t292;
t293 = sin(qJ(2));
t367 = t293 * t383;
t290 = cos(pkin(6));
t432 = pkin(1) * t296;
t371 = t290 * t432;
t243 = -pkin(8) * t367 + qJD(1) * t371;
t319 = (pkin(2) * t293 - pkin(9) * t296) * t289;
t244 = qJD(1) * t319;
t165 = t295 * t243 + t292 * t244;
t152 = pkin(10) * t367 + t165;
t283 = t290 * t293 * pkin(1);
t354 = pkin(3) * t292 - pkin(10) * t295;
t394 = t289 * t296;
t169 = (t283 + (pkin(8) + t354) * t394) * qJD(1);
t291 = sin(qJ(4));
t294 = cos(qJ(4));
t88 = t294 * t152 + t291 * t169;
t523 = -qJD(5) * t295 - t88 + (-t359 + t379) * qJ(5);
t279 = qJD(1) * t290 + qJD(2);
t222 = t279 * t292 + t295 * t367;
t171 = t294 * t222 - t291 * t325;
t522 = t491 * t171;
t170 = t291 * t222 + t294 * t325;
t521 = t491 * t170;
t384 = pkin(8) * t394 + t283;
t248 = t384 * qJD(2);
t235 = qJD(1) * t248;
t380 = qJD(2) * t296;
t364 = t295 * t380;
t378 = qJD(3) * t295;
t179 = t279 * t378 + (-t293 * t379 + t364) * t383;
t404 = t179 * Ifges(4,4);
t365 = t292 * t380;
t180 = t279 * t379 + (t293 * t378 + t365) * t383;
t460 = pkin(4) + pkin(5);
t196 = -pkin(2) * t279 - t243;
t221 = t295 * t279 - t292 * t367;
t119 = -pkin(3) * t221 - pkin(10) * t222 + t196;
t246 = t384 * qJD(1);
t197 = pkin(9) * t279 + t246;
t238 = (-pkin(2) * t296 - pkin(9) * t293 - pkin(1)) * t289;
t212 = qJD(1) * t238;
t142 = t295 * t197 + t292 * t212;
t123 = -pkin(10) * t325 + t142;
t376 = qJD(4) * t294;
t377 = qJD(4) * t291;
t382 = qJD(2) * t289;
t361 = qJD(1) * t382;
t358 = t293 * t361;
t245 = qJD(2) * t319;
t233 = qJD(1) * t245;
t395 = t289 * t293;
t280 = pkin(8) * t395;
t256 = -t280 + t371;
t247 = t256 * qJD(2);
t234 = qJD(1) * t247;
t71 = -t197 * t379 + t212 * t378 + t292 * t233 + t295 * t234;
t61 = pkin(10) * t358 + t71;
t90 = pkin(3) * t180 - pkin(10) * t179 + t235;
t9 = -t119 * t377 - t123 * t376 - t291 * t61 + t294 * t90;
t97 = -qJD(4) * t170 + t294 * t179 + t291 * t358;
t1 = -qJ(6) * t97 - qJD(6) * t171 - t180 * t460 - t9;
t215 = qJD(4) - t221;
t8 = t119 * t376 - t123 * t377 + t291 * t90 + t294 * t61;
t6 = t180 * qJ(5) + t215 * qJD(5) + t8;
t98 = qJD(4) * t171 + t291 * t179 - t294 * t358;
t2 = qJ(6) * t98 + qJD(6) * t170 + t6;
t7 = -pkin(4) * t180 - t9;
t495 = t9 * mrSges(5,1) - t7 * mrSges(6,1) - t1 * mrSges(7,1) - t8 * mrSges(5,2) + t2 * mrSges(7,2) + t6 * mrSges(6,3);
t520 = t495 - t404 / 0.2e1 + t235 * mrSges(4,1) - t71 * mrSges(4,3);
t141 = -t292 * t197 + t295 * t212;
t122 = pkin(3) * t325 - t141;
t41 = t294 * t119 - t291 * t123;
t29 = qJ(6) * t171 + t41;
t479 = qJD(5) - t29;
t17 = -t215 * t460 + t479;
t307 = -t171 * qJ(5) + t122;
t28 = -t170 * t460 + qJD(6) - t307;
t478 = qJD(5) - t41;
t36 = -pkin(4) * t215 + t478;
t40 = t170 * pkin(4) + t307;
t168 = Ifges(5,4) * t170;
t467 = t471 * t171 + t215 * t470 - t168 + t521;
t519 = t17 * mrSges(7,3) + t40 * mrSges(6,3) + t41 * mrSges(5,3) - t467 / 0.2e1 - t122 * mrSges(5,2) - t28 * mrSges(7,2) - t36 * mrSges(6,2);
t202 = t215 * qJ(5);
t42 = t291 * t119 + t294 * t123;
t30 = qJ(6) * t170 + t42;
t27 = t202 + t30;
t37 = t202 + t42;
t484 = t170 * t489 + t215 * t487 + t522;
t417 = Ifges(5,4) * t171;
t78 = -Ifges(5,2) * t170 + Ifges(5,6) * t215 + t417;
t518 = t28 * mrSges(7,1) + t37 * mrSges(6,2) + t42 * mrSges(5,3) - t484 / 0.2e1 + t78 / 0.2e1 - t122 * mrSges(5,1) - t27 * mrSges(7,3) - t40 * mrSges(6,1);
t455 = -t170 / 0.2e1;
t454 = t170 / 0.2e1;
t453 = -t171 / 0.2e1;
t452 = t171 / 0.2e1;
t442 = -t215 / 0.2e1;
t517 = Ifges(4,6) * t524;
t516 = t325 * t525;
t392 = t294 * t296;
t207 = (t291 * t293 + t295 * t392) * t383;
t263 = t354 * qJD(3);
t431 = pkin(9) * t291;
t368 = -pkin(4) - t431;
t372 = qJD(6) * t294;
t273 = -pkin(3) * t295 - pkin(10) * t292 - pkin(2);
t285 = t294 * t295 * pkin(9);
t385 = qJD(4) * t285 + t273 * t377;
t87 = -t291 * t152 + t169 * t294;
t514 = qJ(6) * t207 + t359 * t460 + (-qJ(6) * t378 - t263) * t294 + (qJ(6) * t377 - t372 + (-pkin(5) + t368) * qJD(3)) * t292 + t385 + t87;
t393 = t291 * t295;
t206 = -t294 * t367 + t363 * t393;
t386 = t291 * t263 + t273 * t376;
t513 = -qJ(6) * t206 + (-pkin(9) * qJD(3) + qJ(6) * qJD(4)) * t294 * t292 + (qJD(6) * t292 + (-pkin(9) * qJD(4) + qJ(6) * qJD(3)) * t295) * t291 + t386 + t523;
t164 = -t292 * t243 + t295 * t244;
t314 = pkin(3) * t367 + t164;
t309 = qJ(5) * t207 + t314;
t398 = qJ(5) * t294;
t318 = -t291 * t460 + t398;
t312 = -pkin(9) + t318;
t374 = qJD(5) * t294;
t399 = qJ(5) * t291;
t465 = -t294 * t460 - t399;
t512 = t206 * t460 + (qJD(4) * t465 + t374) * t292 + t312 * t378 - t309;
t158 = (-t294 * t379 - t295 * t377) * pkin(9) + t386;
t511 = t158 + t523;
t353 = -t263 * t294 + t385;
t510 = pkin(4) * t359 + t368 * t379 + t353 + t87;
t326 = pkin(4) * t291 - t398;
t321 = pkin(9) + t326;
t327 = pkin(4) * t294 + t399;
t509 = -pkin(4) * t206 + (qJD(4) * t327 - t374) * t292 + t321 * t378 + t309;
t508 = t491 * t294;
t507 = t491 * t291;
t214 = Ifges(4,4) * t221;
t420 = Ifges(4,1) * t222;
t505 = t214 / 0.2e1 + t420 / 0.2e1;
t409 = Ifges(4,2) * t221;
t418 = Ifges(4,4) * t222;
t504 = -t409 / 0.2e1 - t418 / 0.2e1;
t72 = -t197 * t378 - t212 * t379 + t295 * t233 - t292 * t234;
t503 = t72 * mrSges(4,3);
t486 = t180 * t487 + t489 * t98 + t491 * t97;
t502 = qJD(5) * t291 + t142;
t501 = t291 * t489 + t508;
t500 = -t291 * t515 + t294 * t492;
t355 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1;
t356 = Ifges(7,6) / 0.2e1 + Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t357 = Ifges(7,5) / 0.2e1 - Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1;
t494 = -t196 * mrSges(4,1) - t41 * mrSges(5,1) + t36 * mrSges(6,1) + t17 * mrSges(7,1) + t42 * mrSges(5,2) - t27 * mrSges(7,2) + t142 * mrSges(4,3) - t37 * mrSges(6,3) + Ifges(7,5) * t452 + Ifges(6,6) * t455 - t504 + t517 + (Ifges(7,6) + Ifges(5,6)) * t454 + t492 * t453 + (Ifges(7,3) + t490) * t442;
t300 = t494 + t418 / 0.2e1 - t355 * t215 + t356 * t170 + t357 * t171;
t499 = t300 + t409 / 0.2e1;
t468 = t180 * t470 + t469 * t98 + t471 * t97;
t416 = Ifges(5,4) * t291;
t498 = t294 * t471 - t416 + t507;
t496 = -t196 * mrSges(4,2) + t141 * mrSges(4,3) - t505 + t516;
t493 = -t279 * Ifges(3,6) / 0.2e1;
t20 = Ifges(5,5) * t97 - Ifges(5,6) * t98 + Ifges(5,3) * t180;
t22 = Ifges(6,4) * t97 + Ifges(6,2) * t180 + Ifges(6,6) * t98;
t485 = t22 + t20;
t483 = t215 * t318 + t502;
t397 = qJ(6) * t221;
t423 = pkin(10) - qJ(6);
t160 = pkin(3) * t222 - pkin(10) * t221;
t67 = t294 * t141 + t291 * t160;
t46 = t222 * qJ(5) + t67;
t482 = -t291 * t397 - t377 * t423 - t372 - t46;
t481 = t215 * t326 - t502;
t131 = t291 * t141;
t275 = t423 * t294;
t480 = qJD(4) * t275 - qJD(6) * t291 - t131 - (-t160 - t397) * t294 + t460 * t222;
t476 = -t214 / 0.2e1 + t496;
t475 = -t294 * t489 + t507;
t474 = t291 * t492 + t294 * t515;
t473 = -t291 * t9 + t294 * t8;
t472 = t291 * t7 + t294 * t6;
t415 = Ifges(5,4) * t294;
t466 = t291 * t471 + t415 - t508;
t322 = t291 * t42 + t294 * t41;
t323 = t291 * t37 - t294 * t36;
t329 = Ifges(7,5) * t294 + Ifges(7,6) * t291;
t340 = -Ifges(5,2) * t291 + t415;
t347 = -mrSges(7,1) * t291 + mrSges(7,2) * t294;
t349 = mrSges(6,1) * t291 - mrSges(6,3) * t294;
t351 = mrSges(5,1) * t291 + mrSges(5,2) * t294;
t433 = t294 / 0.2e1;
t435 = t291 / 0.2e1;
t436 = -t291 / 0.2e1;
t441 = t215 / 0.2e1;
t464 = t323 * mrSges(6,2) + t322 * mrSges(5,3) + (t17 * t294 - t27 * t291) * mrSges(7,3) - t122 * t351 - t28 * t347 - t329 * t442 - t340 * t455 - t349 * t40 - t436 * t78 - t501 * t454 - t500 * t441 - t484 * t435 - t498 * t452 - t467 * t433;
t463 = t97 / 0.2e1;
t462 = -t98 / 0.2e1;
t461 = t98 / 0.2e1;
t459 = pkin(1) * mrSges(3,1);
t458 = pkin(1) * mrSges(3,2);
t451 = -t180 / 0.2e1;
t450 = t180 / 0.2e1;
t253 = -t290 * t295 + t292 * t395;
t440 = -t253 / 0.2e1;
t254 = t290 * t292 + t295 * t395;
t438 = t254 / 0.2e1;
t437 = t290 / 0.2e1;
t434 = -t294 / 0.2e1;
t426 = t71 * mrSges(4,2);
t425 = t72 * mrSges(4,1);
t424 = qJD(3) / 0.2e1;
t422 = mrSges(5,3) * t170;
t421 = mrSges(5,3) * t171;
t419 = Ifges(3,4) * t293;
t405 = t179 * Ifges(4,1);
t403 = t180 * Ifges(4,4);
t400 = qJ(5) * t170;
t396 = qJ(6) * t292;
t124 = -mrSges(6,2) * t170 + mrSges(6,3) * t215;
t126 = -mrSges(5,2) * t215 - t422;
t391 = -t124 - t126;
t125 = mrSges(7,2) * t215 + mrSges(7,3) * t170;
t390 = t124 + t125;
t128 = mrSges(5,1) * t215 - t421;
t129 = -mrSges(6,1) * t215 + mrSges(6,2) * t171;
t389 = -t128 + t129;
t236 = t280 + (-pkin(2) - t432) * t290;
t147 = pkin(3) * t253 - pkin(10) * t254 + t236;
t237 = pkin(9) * t290 + t384;
t162 = t295 * t237 + t292 * t238;
t149 = -pkin(10) * t394 + t162;
t64 = t291 * t147 + t294 * t149;
t109 = mrSges(5,1) * t170 + mrSges(5,2) * t171;
t182 = -mrSges(4,1) * t325 - t222 * mrSges(4,3);
t388 = t182 - t109;
t161 = -t292 * t237 + t295 * t238;
t387 = mrSges(3,1) * t279 + mrSges(4,1) * t221 - mrSges(4,2) * t222 - mrSges(3,3) * t367;
t232 = t291 * t273 + t285;
t381 = qJD(2) * t293;
t370 = t291 * t394;
t49 = t253 * qJ(5) + t64;
t369 = Ifges(4,5) * t179 - Ifges(4,6) * t180 + Ifges(4,3) * t358;
t148 = pkin(3) * t394 - t161;
t366 = t289 * t381;
t362 = -t394 / 0.2e1;
t32 = -t98 * mrSges(7,1) + t97 * mrSges(7,2);
t360 = t381 / 0.2e1;
t53 = -t180 * mrSges(6,1) + t97 * mrSges(6,2);
t51 = -t180 * mrSges(7,1) - t97 * mrSges(7,3);
t66 = t160 * t294 - t131;
t63 = t147 * t294 - t291 * t149;
t284 = pkin(9) * t393;
t231 = t273 * t294 - t284;
t100 = -t237 * t378 - t238 * t379 + t295 * t245 - t292 * t247;
t199 = -qJ(5) * t295 + t232;
t352 = mrSges(5,1) * t294 - mrSges(5,2) * t291;
t350 = mrSges(6,1) * t294 + mrSges(6,3) * t291;
t348 = mrSges(7,1) * t294 + mrSges(7,2) * t291;
t339 = Ifges(5,2) * t294 + t416;
t187 = -qJD(3) * t253 + t289 * t364;
t188 = qJD(3) * t254 + t289 * t365;
t334 = Ifges(4,5) * t187 - Ifges(4,6) * t188;
t328 = Ifges(7,5) * t291 - Ifges(7,6) * t294;
t110 = pkin(3) * t188 - pkin(10) * t187 + t248;
t99 = -t237 * t379 + t238 * t378 + t292 * t245 + t295 * t247;
t85 = pkin(10) * t366 + t99;
t15 = t110 * t294 - t147 * t377 - t149 * t376 - t291 * t85;
t190 = t254 * t294 - t370;
t320 = qJ(5) * t190 - t148;
t189 = t254 * t291 + t289 * t392;
t14 = t291 * t110 + t147 * t376 - t149 * t377 + t294 * t85;
t276 = Ifges(3,4) * t363;
t311 = t279 * Ifges(3,5) - t243 * mrSges(3,3) + Ifges(3,1) * t367 / 0.2e1 + t276 / 0.2e1;
t18 = Ifges(7,5) * t97 + Ifges(7,6) * t98 - Ifges(7,3) * t180;
t310 = pkin(3) * t366 + t100;
t11 = t188 * qJ(5) + t253 * qJD(5) + t14;
t308 = pkin(3) * t358 + t72;
t306 = -t420 / 0.2e1 + t476;
t116 = -qJD(4) * t189 + t187 * t294 + t291 * t366;
t304 = qJ(5) * t116 + qJD(5) * t190 + t310;
t303 = qJ(5) * t97 + qJD(5) * t171 + t308;
t302 = t141 * mrSges(4,1) + Ifges(4,6) * t221 + Ifges(4,5) * t222 + t493 - (Ifges(3,2) * t296 + t419) * t383 / 0.2e1 - t142 * mrSges(4,2) - t246 * mrSges(3,3) + (t524 + t424) * Ifges(4,3);
t288 = t295 * pkin(4);
t274 = t423 * t291;
t272 = Ifges(3,5) * t296 * t361;
t267 = -pkin(3) - t327;
t260 = pkin(3) - t465;
t241 = -mrSges(3,2) * t279 + mrSges(3,3) * t363;
t239 = t321 * t292;
t200 = -t231 + t288;
t198 = t312 * t292;
t183 = t291 * t396 + t199;
t181 = mrSges(4,2) * t325 + t221 * mrSges(4,3);
t178 = pkin(5) * t295 + t284 + t288 + (-t273 - t396) * t294;
t159 = t379 * t431 - t353;
t155 = -mrSges(4,2) * t358 - mrSges(4,3) * t180;
t154 = mrSges(4,1) * t358 - mrSges(4,3) * t179;
t127 = -mrSges(7,1) * t215 - mrSges(7,3) * t171;
t117 = -qJD(4) * t370 + t187 * t291 + t254 * t376 - t294 * t366;
t113 = mrSges(4,1) * t180 + mrSges(4,2) * t179;
t108 = -mrSges(7,1) * t170 + mrSges(7,2) * t171;
t107 = mrSges(6,1) * t170 - mrSges(6,3) * t171;
t106 = pkin(4) * t171 + t400;
t102 = Ifges(4,5) * t358 - t403 + t405;
t101 = -t180 * Ifges(4,2) + Ifges(4,6) * t358 + t404;
t65 = pkin(4) * t189 - t320;
t58 = -t171 * t460 - t400;
t56 = -mrSges(6,2) * t98 + mrSges(6,3) * t180;
t55 = -mrSges(5,2) * t180 - mrSges(5,3) * t98;
t54 = mrSges(7,2) * t180 + mrSges(7,3) * t98;
t52 = mrSges(5,1) * t180 - mrSges(5,3) * t97;
t50 = -pkin(4) * t253 - t63;
t48 = -pkin(4) * t222 - t66;
t44 = -t189 * t460 + t320;
t38 = qJ(6) * t189 + t49;
t35 = -qJ(6) * t190 - t253 * t460 - t63;
t33 = mrSges(5,1) * t98 + mrSges(5,2) * t97;
t31 = mrSges(6,1) * t98 - mrSges(6,3) * t97;
t23 = t97 * Ifges(5,4) - t98 * Ifges(5,2) + t180 * Ifges(5,6);
t16 = pkin(4) * t117 - t304;
t13 = -pkin(4) * t188 - t15;
t12 = pkin(4) * t98 - t303;
t10 = -t117 * t460 + t304;
t5 = -t460 * t98 + t303;
t4 = qJ(6) * t117 + qJD(6) * t189 + t11;
t3 = -qJ(6) * t116 - qJD(6) * t190 - t188 * t460 - t15;
t19 = [(-t496 + t505) * t187 + (-mrSges(3,1) * t290 + mrSges(4,2) * t254 + mrSges(3,3) * t395) * t235 + (-t296 * t334 / 0.2e1 + ((Ifges(3,5) * t437 - t256 * mrSges(3,3) + (-0.2e1 * t458 + 0.3e1 / 0.2e1 * Ifges(3,4) * t296) * t289) * t296 + (-Ifges(3,6) * t290 + Ifges(4,5) * t438 + Ifges(4,6) * t440 - t384 * mrSges(3,3) + (-0.2e1 * t459 - 0.3e1 / 0.2e1 * t419 + (-Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t296) * t289) * t293) * qJD(2)) * t383 + m(3) * (t234 * t384 - t235 * t256 - t243 * t248 + t246 * t247) + t179 * (Ifges(4,1) * t254 - Ifges(4,5) * t394) / 0.2e1 + (Ifges(4,4) * t254 - Ifges(4,6) * t394) * t451 + (Ifges(5,4) * t455 + Ifges(7,5) * t442 + t492 * t441 + t471 * t452 + t491 * t454 - t519) * t116 + t272 * t437 + t102 * t438 + t101 * t440 + t18 * t440 + t334 * t424 + t394 * t426 - t394 * t425 - t387 * t248 + m(7) * (t1 * t35 + t10 * t28 + t17 * t3 + t2 * t38 + t27 * t4 + t44 * t5) + m(6) * (t11 * t37 + t12 * t65 + t13 * t36 + t16 * t40 + t49 * t6 + t50 * t7) + m(4) * (t100 * t141 + t142 * t99 + t161 * t72 + t162 * t71 + t196 * t248 + t235 * t236) + t234 * (-mrSges(3,2) * t290 + mrSges(3,3) * t394) + (t470 * t463 + (-Ifges(4,2) - Ifges(7,3)) * t451 + Ifges(5,6) * t462 + t485 / 0.2e1 + t487 * t461 + t490 * t450 + t520) * t253 + t369 * t362 + (Ifges(5,6) * t455 - Ifges(7,3) * t442 + t490 * t441 + t470 * t452 + t487 * t454 - t494 + t504) * t188 + (t468 / 0.2e1 + t471 * t463 + Ifges(7,5) * t451 + Ifges(5,4) * t462 - t9 * mrSges(5,3) + t7 * mrSges(6,2) - t1 * mrSges(7,3) - t308 * mrSges(5,2) - t12 * mrSges(6,3) + t5 * mrSges(7,2) + t491 * t461 + t492 * t450) * t190 + m(5) * (-t122 * t310 + t14 * t42 - t148 * t308 + t15 * t41 + t63 * t9 + t64 * t8) - t310 * t109 + (t469 * t463 + Ifges(7,6) * t451 - Ifges(5,2) * t462 - t6 * mrSges(6,2) + t2 * mrSges(7,3) - t8 * mrSges(5,3) - t23 / 0.2e1 - t308 * mrSges(5,1) + t12 * mrSges(6,1) - t5 * mrSges(7,1) + t486 / 0.2e1 + t489 * t461 - t515 * t450) * t189 + (-Ifges(5,2) * t455 + Ifges(7,6) * t442 - t441 * t515 + t469 * t452 + t489 * t454 - t518) * t117 + t247 * t241 + t236 * t113 + t99 * t181 + t100 * t182 + t161 * t154 + t162 * t155 + t148 * t33 + t11 * t124 + t4 * t125 + t14 * t126 + t3 * t127 + t15 * t128 + t13 * t129 + t16 * t107 + t10 * t108 + t63 * t52 + t64 * t55 + t65 * t31 + t35 * t51 + t50 * t53 + t38 * t54 + t49 * t56 + t44 * t32 - t254 * t503 + (t311 * t296 + (t493 + t302) * t293) * t382; (((-m(4) * t142 - t181) * pkin(9) - qJD(3) * Ifges(4,6) / 0.2e1 - t499) * qJD(3) + (-m(4) * t72 - m(5) * t308 - t154 + t33) * pkin(9) + t23 * t436 + (-t291 * t6 + t294 * t7) * mrSges(6,2) - t308 * t351 + t5 * t347 + t12 * t349 - t503 + t235 * mrSges(4,2) + t329 * t451 + t405 / 0.2e1 - t403 / 0.2e1 + t340 * t462 + t102 / 0.2e1 + (-t1 * t294 + t2 * t291) * mrSges(7,3) + (Ifges(4,5) * t360 + ((qJD(1) * t362 + qJD(3)) * Ifges(4,6) + t499) * t296) * t383 + (-t291 * t8 - t294 * t9) * mrSges(5,3) + t501 * t461 + t500 * t450 + t486 * t435 + t498 * t463 + t468 * t433 + (t78 * t434 + t40 * t350 - t28 * t348 + t328 * t441 + t339 * t454 + t122 * t352 + (t17 * t291 + t27 * t294) * mrSges(7,3) + (t291 * t41 - t294 * t42) * mrSges(5,3) + (-t291 * t36 - t294 * t37) * mrSges(6,2) + t475 * t455 + t474 * t442 + t484 * t433 + t466 * t453 + t467 * t436) * qJD(4)) * t292 + (-m(4) * t235 - t113) * pkin(2) + (t159 - t87) * t128 + ((-t276 / 0.2e1 + t383 * t458 - t311) * t296 + ((t459 + t419 / 0.2e1 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t296) * t383 + (t279 / 0.2e1 - qJD(2)) * Ifges(3,6) - t302) * t293) * t383 + ((-Ifges(4,2) / 0.2e1 - t355) * t180 + (-t464 + (-m(4) * t141 + m(5) * t122 - t388) * pkin(9) + (t424 - t363) * Ifges(4,5) - t306) * qJD(3) - t22 / 0.2e1 + (m(4) * t71 + t155) * pkin(9) + t356 * t98 + t357 * t97 - t20 / 0.2e1 + t101 / 0.2e1 + t18 / 0.2e1 + (Ifges(4,6) * t360 + (t363 * t525 + t306) * t296) * t383 - t520) * t295 - m(4) * (t141 * t164 + t142 * t165 + t196 * t246) + m(5) * (t158 * t42 + t159 * t41 + t231 * t9 + t232 * t8) + (Ifges(5,4) * t454 + Ifges(7,5) * t441 + t442 * t492 + t471 * t453 + t455 * t491 + t519) * t207 + (t158 - t88) * t126 + t387 * t246 + t509 * t107 + t510 * t129 + (t12 * t239 + t199 * t6 + t200 * t7 + t36 * t510 + t37 * t511 + t40 * t509) * m(6) + t511 * t124 + t512 * t108 + t513 * t125 + (t1 * t178 + t17 * t514 + t183 * t2 + t198 * t5 + t27 * t513 + t28 * t512) * m(7) + t514 * t127 + (-Ifges(5,2) * t454 + Ifges(7,6) * t441 - t442 * t515 + t453 * t469 + t455 * t489 + t518) * t206 + t272 - m(5) * (-t122 * t314 + t41 * t87 + t42 * t88) + t314 * t109 + t239 * t31 - t243 * t241 + t231 * t52 + t232 * t55 - t234 * mrSges(3,2) - t235 * mrSges(3,1) + t199 * t56 + t200 * t53 + t198 * t32 - t165 * t181 - t164 * t182 + t183 * t54 + t178 * t51; t472 * mrSges(6,2) + ((-m(5) * t322 - m(6) * t323 + t291 * t391 + t294 * t389) * qJD(4) + m(5) * t473 + m(6) * t472 + (t55 + t56) * t294 + (-t52 + t53) * t291) * pkin(10) + t473 * mrSges(5,3) + t474 * t450 + t475 * t461 + t468 * t435 + t466 * t463 - t464 * qJD(4) + t23 * t433 + t369 + t388 * t142 + t5 * t348 - t12 * t350 + t328 * t451 + (pkin(3) * t308 - t122 * t142 - t41 * t66 - t42 * t67) * m(5) + t308 * t352 + t425 - t426 + t339 * t462 + (-t1 * t291 - t2 * t294) * mrSges(7,3) + t267 * t31 + t274 * t51 + t275 * t54 + t260 * t32 - t141 * t181 - t46 * t124 - t67 * t126 - t66 * t128 - t48 * t129 - pkin(3) * t33 + t480 * t127 + (t12 * t267 - t36 * t48 - t37 * t46 + t40 * t481) * m(6) + t481 * t107 + t482 * t125 + t483 * t108 + (t1 * t274 + t17 * t480 + t2 * t275 + t260 * t5 + t27 * t482 + t28 * t483) * m(7) + (t516 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t222 + t464 + t476) * t221 + (t300 + t517) * t222 + t486 * t434; (-pkin(4) * t7 + qJ(5) * t6 - t106 * t40 - t36 * t42 + t37 * t478) * m(6) + (-Ifges(5,2) * t171 - t168 + t467) * t454 + t495 + (-Ifges(7,5) * t170 + Ifges(7,6) * t171) * t441 + t390 * qJD(5) + (-t17 * t170 - t171 * t27) * mrSges(7,3) + (-t389 + t421) * t42 + (t391 - t422) * t41 + t78 * t452 + (t54 + t56) * qJ(5) + (t170 * t36 + t171 * t37) * mrSges(6,2) + (-t170 * t471 - t417 + t484 + t522) * t453 + t485 + (t2 * qJ(5) - t1 * t460 - t17 * t30 + t27 * t479 - t28 * t58) * m(7) - t460 * t51 + (-t170 * t492 - t171 * t515) * t442 - t122 * (mrSges(5,1) * t171 - mrSges(5,2) * t170) - t40 * (mrSges(6,1) * t171 + mrSges(6,3) * t170) - t28 * (-mrSges(7,1) * t171 - mrSges(7,2) * t170) - t29 * t125 - t30 * t127 - t106 * t107 - t58 * t108 - pkin(4) * t53 + (t171 * t489 - t521) * t455 - t18; -t390 * t215 + (t107 - t108) * t171 + t51 + t53 + (-t171 * t28 - t215 * t27 + t1) * m(7) + (t171 * t40 - t215 * t37 + t7) * m(6); -t170 * t125 + t171 * t127 + 0.2e1 * (t5 / 0.2e1 + t17 * t452 + t27 * t455) * m(7) + t32;];
tauc  = t19(:);
