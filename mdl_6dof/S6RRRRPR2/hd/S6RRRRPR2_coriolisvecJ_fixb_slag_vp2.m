% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:25
% EndTime: 2019-03-09 21:57:50
% DurationCPUTime: 12.58s
% Computational Cost: add. (23648->697), mult. (60760->954), div. (0->0), fcn. (45911->10), ass. (0->332)
t347 = sin(pkin(11));
t348 = cos(pkin(11));
t349 = sin(qJ(6));
t353 = cos(qJ(6));
t369 = t347 * t349 - t348 * t353;
t293 = t369 * qJD(6);
t351 = sin(qJ(3));
t352 = sin(qJ(2));
t354 = cos(qJ(3));
t355 = cos(qJ(2));
t397 = t354 * t355;
t368 = t351 * t352 - t397;
t366 = t368 * qJD(1);
t447 = cos(qJ(4));
t282 = t447 * t366;
t312 = t351 * t355 + t354 * t352;
t298 = t312 * qJD(1);
t350 = sin(qJ(4));
t245 = t298 * t350 + t282;
t496 = t369 * t245;
t504 = t496 + t293;
t498 = t245 * t347;
t234 = pkin(5) * t498;
t389 = pkin(10) * t498;
t310 = t347 * t353 + t348 * t349;
t497 = t310 * t245;
t503 = Ifges(7,1) * t496 + Ifges(7,4) * t497;
t502 = Ifges(7,4) * t496 + Ifges(7,2) * t497;
t501 = Ifges(7,5) * t496 + Ifges(7,6) * t497;
t466 = -pkin(8) - pkin(7);
t328 = t466 * t355;
t317 = qJD(1) * t328;
t299 = t351 * t317;
t327 = t466 * t352;
t316 = qJD(1) * t327;
t424 = qJD(2) * pkin(2);
t307 = t316 + t424;
t260 = t354 * t307 + t299;
t290 = t298 * pkin(9);
t220 = t260 - t290;
t346 = qJD(2) + qJD(3);
t213 = pkin(3) * t346 + t220;
t302 = t354 * t317;
t261 = t351 * t307 - t302;
t365 = pkin(9) * t366;
t221 = t261 - t365;
t214 = t350 * t221;
t143 = t213 * t447 - t214;
t390 = -qJD(3) - qJD(4);
t344 = qJD(2) - t390;
t133 = -t344 * pkin(4) + qJD(5) - t143;
t363 = t350 * t366;
t361 = t298 * t447 - t363;
t380 = t348 * t344 - t347 * t361;
t106 = -pkin(5) * t380 + t133;
t267 = t346 * t368;
t256 = t267 * qJD(1);
t268 = t346 * t312;
t257 = t268 * qJD(1);
t140 = -t447 * t256 + (-qJD(4) * t298 - t257) * t350 - qJD(4) * t282;
t381 = qJD(4) * t447;
t141 = -qJD(4) * t363 - t350 * t256 + t257 * t447 + t298 * t381;
t385 = qJD(2) * t466;
t392 = qJD(3) * t354;
t393 = qJD(3) * t351;
t202 = t298 * t385 + t307 * t392 + t317 * t393;
t150 = -pkin(9) * t257 + t202;
t313 = t351 * t327;
t203 = -t261 * qJD(3) + (t397 * t466 - t313) * qJD(2) * qJD(1);
t359 = t256 * pkin(9) + t203;
t391 = qJD(4) * t350;
t53 = t447 * t150 + t213 * t381 - t221 * t391 + t350 * t359;
t51 = qJD(5) * t344 + t53;
t395 = qJD(1) * t352;
t342 = pkin(2) * t395;
t232 = pkin(3) * t257 + qJD(2) * t342;
t60 = pkin(4) * t141 - qJ(5) * t140 - qJD(5) * t361 + t232;
t16 = t347 * t60 + t348 * t51;
t412 = t140 * t347;
t11 = -pkin(10) * t412 + t16;
t225 = t344 * t347 + t348 * t361;
t215 = t447 * t221;
t144 = t350 * t213 + t215;
t134 = t344 * qJ(5) + t144;
t340 = -t355 * pkin(2) - pkin(1);
t280 = pkin(3) * t368 + t340;
t273 = t280 * qJD(1);
t163 = t245 * pkin(4) - qJ(5) * t361 + t273;
t87 = -t134 * t347 + t348 * t163;
t59 = pkin(5) * t245 - pkin(10) * t225 + t87;
t88 = t348 * t134 + t347 * t163;
t66 = pkin(10) * t380 + t88;
t19 = -t349 * t66 + t353 * t59;
t15 = -t347 * t51 + t348 * t60;
t411 = t140 * t348;
t6 = pkin(5) * t141 - pkin(10) * t411 + t15;
t2 = qJD(6) * t19 + t11 * t353 + t349 * t6;
t20 = t349 * t59 + t353 * t66;
t294 = t310 * qJD(6);
t54 = qJD(4) * t144 + t350 * t150 - t447 * t359;
t34 = pkin(5) * t412 + t54;
t430 = Ifges(6,2) * t347;
t436 = Ifges(6,4) * t348;
t373 = -t430 + t436;
t40 = t141 * Ifges(6,6) + t140 * t373;
t437 = Ifges(6,4) * t347;
t374 = Ifges(6,1) * t348 - t437;
t41 = t141 * Ifges(6,5) + t140 * t374;
t421 = t16 * t348;
t448 = t348 / 0.2e1;
t242 = qJD(6) + t245;
t455 = t242 / 0.2e1;
t157 = t225 * t353 + t349 * t380;
t459 = t157 / 0.2e1;
t495 = -t225 * t349 + t353 * t380;
t461 = t495 / 0.2e1;
t463 = t141 / 0.2e1;
t149 = Ifges(7,4) * t495;
t79 = Ifges(7,1) * t157 + Ifges(7,5) * t242 + t149;
t467 = t79 / 0.2e1;
t435 = Ifges(7,4) * t157;
t78 = Ifges(7,2) * t495 + Ifges(7,6) * t242 + t435;
t469 = t78 / 0.2e1;
t46 = -qJD(6) * t157 - t140 * t310;
t471 = t46 / 0.2e1;
t45 = qJD(6) * t495 - t140 * t369;
t472 = t45 / 0.2e1;
t473 = Ifges(7,1) * t472 + Ifges(7,4) * t471 + Ifges(7,5) * t463;
t474 = Ifges(7,4) * t472 + Ifges(7,2) * t471 + Ifges(7,6) * t463;
t500 = -t369 * t474 + (Ifges(6,5) * t347 + Ifges(7,5) * t310 + Ifges(6,6) * t348 - Ifges(7,6) * t369) * t463 + (-t2 * t369 - t20 * t294) * mrSges(7,3) + t34 * (mrSges(7,1) * t369 + mrSges(7,2) * t310) + (Ifges(7,4) * t310 - Ifges(7,2) * t369) * t471 + (Ifges(7,1) * t310 - Ifges(7,4) * t369) * t472 + (-Ifges(7,5) * t293 - Ifges(7,6) * t294) * t455 + (-Ifges(7,1) * t293 - Ifges(7,4) * t294) * t459 + (-Ifges(7,4) * t293 - Ifges(7,2) * t294) * t461 - t496 * t79 / 0.2e1 - t497 * t78 / 0.2e1 + (-mrSges(6,1) * t348 + mrSges(6,2) * t347 - mrSges(5,1)) * t54 + t347 * t41 / 0.2e1 - (Ifges(6,2) * t348 + t437) * t412 / 0.2e1 + (Ifges(6,1) * t347 + t436) * t411 / 0.2e1 + mrSges(6,3) * t421 - t53 * mrSges(5,2) + t40 * t448 - t293 * t467 - t294 * t469 + t310 * t473 + Ifges(5,5) * t140 - Ifges(5,6) * t141 + (-t504 * mrSges(7,2) + (t294 + t497) * mrSges(7,1)) * t106;
t452 = t245 / 0.2e1;
t499 = -t347 / 0.2e1;
t238 = Ifges(5,4) * t245;
t195 = pkin(4) * t361 + qJ(5) * t245;
t494 = -t361 / 0.2e1;
t493 = (t225 * Ifges(6,4) + Ifges(6,2) * t380 + Ifges(6,6) * t245) * t499;
t491 = pkin(5) * t361;
t438 = Ifges(5,4) * t361;
t445 = pkin(3) * t350;
t336 = qJ(5) + t445;
t305 = (-pkin(10) - t336) * t347;
t345 = t348 * pkin(10);
t403 = t336 * t348;
t306 = t345 + t403;
t249 = t305 * t353 - t306 * t349;
t379 = pkin(3) * t381;
t332 = t379 + qJD(5);
t406 = t245 * t348;
t376 = pkin(10) * t406 + t491;
t158 = t220 * t447 - t214;
t446 = pkin(3) * t298;
t171 = t195 + t446;
t91 = -t158 * t347 + t348 * t171;
t64 = t376 + t91;
t92 = t348 * t158 + t347 * t171;
t80 = t389 + t92;
t490 = qJD(6) * t249 - t332 * t369 - t349 * t64 - t353 * t80;
t250 = t305 * t349 + t306 * t353;
t489 = -qJD(6) * t250 - t310 * t332 + t349 * t80 - t353 * t64;
t322 = (-pkin(10) - qJ(5)) * t347;
t415 = qJ(5) * t348;
t323 = t345 + t415;
t269 = t322 * t353 - t323 * t349;
t96 = -t143 * t347 + t348 * t195;
t67 = t245 * t345 + t491 + t96;
t97 = t348 * t143 + t347 * t195;
t84 = t97 + t389;
t488 = -qJD(5) * t369 + qJD(6) * t269 - t349 * t67 - t353 * t84;
t270 = t322 * t349 + t323 * t353;
t487 = -qJD(5) * t310 - qJD(6) * t270 + t349 * t84 - t353 * t67;
t339 = pkin(2) * t354 + pkin(3);
t384 = t447 * t351;
t292 = pkin(2) * t384 + t350 * t339;
t287 = qJ(5) + t292;
t271 = (-pkin(10) - t287) * t347;
t404 = t287 * t348;
t272 = t345 + t404;
t207 = t271 * t349 + t272 * t353;
t334 = t350 * t351 * pkin(2);
t254 = t447 * pkin(2) * t392 + t334 * t390 + t339 * t381;
t253 = qJD(5) + t254;
t265 = -t351 * t316 + t302;
t226 = t265 + t365;
t266 = t354 * t316 + t299;
t227 = -t290 + t266;
t165 = t350 * t226 + t227 * t447;
t166 = t171 + t342;
t93 = -t165 * t347 + t348 * t166;
t65 = t376 + t93;
t94 = t348 * t165 + t347 * t166;
t81 = t389 + t94;
t486 = -qJD(6) * t207 - t253 * t310 + t349 * t81 - t353 * t65;
t206 = t271 * t353 - t272 * t349;
t485 = qJD(6) * t206 - t253 * t369 - t349 * t65 - t353 * t81;
t396 = mrSges(5,1) * t344 + mrSges(6,1) * t380 - mrSges(6,2) * t225 - mrSges(5,3) * t361;
t89 = -mrSges(7,1) * t495 + mrSges(7,2) * t157;
t484 = t89 - t396;
t375 = mrSges(6,1) * t347 + mrSges(6,2) * t348;
t483 = t133 * t375;
t274 = t354 * t327 + t328 * t351;
t239 = -pkin(9) * t312 + t274;
t275 = -t354 * t328 + t313;
t240 = -pkin(9) * t368 + t275;
t480 = t447 * t239 - t350 * t240;
t479 = t254 - t165;
t164 = -t447 * t226 + t227 * t350;
t255 = t339 * t391 + (t351 * t381 + (t350 * t354 + t384) * qJD(3)) * pkin(2);
t478 = t255 - t164;
t3 = -qJD(6) * t20 - t11 * t349 + t353 * t6;
t477 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t45 + Ifges(7,6) * t46;
t428 = Ifges(6,6) * t380;
t433 = Ifges(6,5) * t225;
t119 = Ifges(6,3) * t245 + t428 + t433;
t429 = Ifges(5,6) * t344;
t186 = -Ifges(5,2) * t245 + t429 + t438;
t434 = Ifges(5,5) * t344;
t441 = Ifges(5,1) * t361;
t187 = -t238 + t434 + t441;
t439 = Ifges(4,4) * t298;
t236 = -Ifges(4,2) * t366 + Ifges(4,6) * t346 + t439;
t289 = Ifges(4,4) * t366;
t237 = t298 * Ifges(4,1) + t346 * Ifges(4,5) - t289;
t326 = t340 * qJD(1);
t364 = mrSges(4,3) * t366;
t427 = Ifges(6,6) * t347;
t432 = Ifges(6,5) * t348;
t372 = -t427 + t432;
t401 = t348 * (t225 * Ifges(6,1) + Ifges(6,4) * t380 + t245 * Ifges(6,5));
t418 = t298 * mrSges(4,3);
t423 = t143 * mrSges(5,3);
t443 = t3 * t310;
t449 = t298 / 0.2e1;
t453 = -t245 / 0.2e1;
t456 = -t242 / 0.2e1;
t460 = -t157 / 0.2e1;
t462 = -t495 / 0.2e1;
t425 = Ifges(7,3) * t242;
t426 = Ifges(7,6) * t495;
t431 = Ifges(7,5) * t157;
t77 = t425 + t426 + t431;
t476 = (-Ifges(4,2) * t298 + t237 - t289) * t366 / 0.2e1 + (t19 * t293 - t443) * mrSges(7,3) - t273 * (mrSges(5,1) * t361 - mrSges(5,2) * t245) - t380 * (Ifges(6,6) * t361 - t245 * t373) / 0.2e1 - t344 * (-Ifges(5,5) * t245 - Ifges(5,6) * t361) / 0.2e1 - t225 * (Ifges(6,5) * t361 - t245 * t374) / 0.2e1 + (Ifges(6,3) * t361 - t245 * t372) * t453 + t245 * t483 - t245 * t423 + t245 * t493 - t20 * (-mrSges(7,2) * t361 + mrSges(7,3) * t497) - t19 * (mrSges(7,1) * t361 - mrSges(7,3) * t496) - t88 * (-mrSges(6,2) * t361 + mrSges(6,3) * t498) + (-Ifges(5,1) * t245 + t119 - t438 + t77) * t494 + (-Ifges(5,2) * t361 + t187 - t238 + t401) * t452 + t500 - t326 * (t298 * mrSges(4,1) - mrSges(4,2) * t366) - t298 * (-Ifges(4,1) * t366 - t439) / 0.2e1 - t346 * (-Ifges(4,5) * t366 - Ifges(4,6) * t298) / 0.2e1 - t87 * (mrSges(6,1) * t361 + mrSges(6,3) * t406) + t361 * t186 / 0.2e1 - Ifges(4,5) * t256 - Ifges(4,6) * t257 - t202 * mrSges(4,2) + t203 * mrSges(4,1) - t260 * t364 + t261 * t418 + t236 * t449 + (Ifges(7,3) * t361 + t501) * t456 + (Ifges(7,6) * t361 + t502) * t462 + (Ifges(7,5) * t361 + t503) * t460;
t475 = t380 * t373 / 0.2e1 + t225 * t374 / 0.2e1 + t273 * mrSges(5,2) + (-t347 * t88 - t348 * t87) * mrSges(6,3) + t187 / 0.2e1 + t434 / 0.2e1 + t493 + t401 / 0.2e1 + t483;
t465 = pkin(1) * mrSges(3,1);
t464 = pkin(1) * mrSges(3,2);
t444 = t15 * mrSges(6,3);
t442 = t348 * pkin(5);
t367 = -t350 * t312 - t368 * t447;
t167 = qJD(4) * t367 - t267 * t447 - t350 * t268;
t264 = t312 * t447 - t350 * t368;
t168 = qJD(4) * t264 - t350 * t267 + t268 * t447;
t343 = t352 * t424;
t248 = pkin(3) * t268 + t343;
t73 = pkin(4) * t168 - qJ(5) * t167 - qJD(5) * t264 + t248;
t318 = t352 * t385;
t319 = t355 * t385;
t209 = t354 * t318 + t351 * t319 + t327 * t392 + t328 * t393;
t177 = -pkin(9) * t268 + t209;
t210 = -qJD(3) * t275 - t318 * t351 + t354 * t319;
t178 = pkin(9) * t267 + t210;
t74 = qJD(4) * t480 + t447 * t177 + t350 * t178;
t22 = t347 * t73 + t348 * t74;
t440 = Ifges(3,4) * t352;
t422 = t15 * t347;
t420 = t480 * t54;
t417 = Ifges(3,5) * qJD(2);
t416 = Ifges(3,6) * qJD(2);
t414 = qJD(2) * mrSges(3,1);
t413 = qJD(2) * mrSges(3,2);
t410 = t144 * t361;
t409 = t167 * t347;
t405 = t264 * t347;
t83 = mrSges(6,1) * t412 + mrSges(6,2) * t411;
t184 = -pkin(4) * t367 - t264 * qJ(5) + t280;
t189 = t350 * t239 + t240 * t447;
t102 = t347 * t184 + t348 * t189;
t394 = qJD(1) * t355;
t388 = t447 * pkin(3);
t386 = Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1;
t383 = t417 / 0.2e1;
t382 = -t416 / 0.2e1;
t13 = -t46 * mrSges(7,1) + t45 * mrSges(7,2);
t21 = -t347 * t74 + t348 * t73;
t101 = t348 * t184 - t189 * t347;
t155 = t220 * t350 + t215;
t338 = -t388 - pkin(4);
t291 = t339 * t447 - t334;
t371 = t421 - t422;
t370 = -t347 * t87 + t348 * t88;
t82 = -pkin(5) * t367 - t264 * t345 + t101;
t90 = -pkin(10) * t405 + t102;
t32 = -t349 * t90 + t353 * t82;
t33 = t349 * t82 + t353 * t90;
t288 = -pkin(4) - t291;
t75 = qJD(4) * t189 + t350 * t177 - t447 * t178;
t357 = t19 * mrSges(7,1) + t273 * mrSges(5,1) + t87 * mrSges(6,1) + t119 / 0.2e1 - t186 / 0.2e1 + t77 / 0.2e1 + t433 / 0.2e1 + t431 / 0.2e1 - t429 / 0.2e1 + t428 / 0.2e1 + t426 / 0.2e1 + t425 / 0.2e1 - t20 * mrSges(7,2) - t88 * mrSges(6,2);
t341 = Ifges(3,4) * t394;
t337 = -pkin(4) - t442;
t325 = mrSges(3,3) * t394 - t413;
t324 = -mrSges(3,3) * t395 + t414;
t321 = t338 - t442;
t297 = Ifges(3,1) * t395 + t341 + t417;
t296 = t416 + (Ifges(3,2) * t355 + t440) * qJD(1);
t279 = t288 - t442;
t278 = mrSges(4,1) * t346 - t418;
t277 = -t346 * mrSges(4,2) - t364;
t276 = t342 + t446;
t259 = mrSges(4,1) * t366 + t298 * mrSges(4,2);
t228 = -mrSges(5,2) * t344 - mrSges(5,3) * t245;
t199 = t369 * t264;
t198 = t310 * t264;
t197 = mrSges(5,1) * t245 + mrSges(5,2) * t361;
t174 = mrSges(6,1) * t245 - mrSges(6,3) * t225;
t173 = -mrSges(6,2) * t245 + mrSges(6,3) * t380;
t137 = Ifges(7,3) * t141;
t130 = pkin(5) * t405 - t480;
t118 = t164 - t234;
t112 = t155 - t234;
t111 = mrSges(7,1) * t242 - mrSges(7,3) * t157;
t110 = -mrSges(7,2) * t242 + mrSges(7,3) * t495;
t108 = t144 - t234;
t86 = mrSges(6,1) * t141 - mrSges(6,3) * t411;
t85 = -mrSges(6,2) * t141 - mrSges(6,3) * t412;
t63 = -t167 * t310 + t264 * t293;
t62 = -t167 * t369 - t264 * t294;
t37 = pkin(5) * t409 + t75;
t31 = -mrSges(7,2) * t141 + mrSges(7,3) * t46;
t30 = mrSges(7,1) * t141 - mrSges(7,3) * t45;
t17 = -pkin(10) * t409 + t22;
t12 = pkin(5) * t168 - t167 * t345 + t21;
t5 = -qJD(6) * t33 + t12 * t353 - t17 * t349;
t4 = qJD(6) * t32 + t12 * t349 + t17 * t353;
t1 = [t326 * (mrSges(4,1) * t268 - mrSges(4,2) * t267) + (-t368 * (-Ifges(4,4) * t267 - Ifges(4,2) * t268) / 0.2e1 + ((-0.2e1 * t464 + 0.3e1 / 0.2e1 * Ifges(3,4) * t355) * t355 + (-0.2e1 * t465 - 0.3e1 / 0.2e1 * t440 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t355 + (m(4) * t340 + mrSges(4,1) * t368 + mrSges(4,2) * t312) * pkin(2)) * t352) * qJD(2)) * qJD(1) + t346 * (-Ifges(4,5) * t267 - Ifges(4,6) * t268) / 0.2e1 + m(6) * (t101 * t15 + t102 * t16 + t133 * t75 + t21 * t87 + t22 * t88 - t420) + m(5) * (-t143 * t75 + t144 * t74 + t189 * t53 + t232 * t280 + t248 * t273 - t420) + (t140 * t367 - t141 * t264 + t167 * t453 + t168 * t494) * Ifges(5,4) - (t15 * mrSges(6,1) - t16 * mrSges(6,2) + t232 * mrSges(5,1) + t137 / 0.2e1 + t372 * t140 + (Ifges(5,2) + Ifges(6,3) + Ifges(7,3) / 0.2e1) * t141 + t477) * t367 + (t372 * t463 + t54 * t375 + t232 * mrSges(5,2) + t40 * t499 + t41 * t448 + (-t15 * t348 - t16 * t347) * mrSges(6,3) + (Ifges(5,1) + Ifges(6,1) * t348 ^ 2 / 0.2e1 + (-t436 + t430 / 0.2e1) * t347) * t140) * t264 + (-t256 * t312 - t267 * t449) * Ifges(4,1) + (t256 * t368 - t257 * t312 - t268 * t449) * Ifges(4,4) + (-t202 * t368 - t203 * t312 + t256 * t274 - t257 * t275 + t260 * t267 - t261 * t268) * mrSges(4,3) + t340 * (mrSges(4,1) * t257 - mrSges(4,2) * t256) + (-t140 * t480 - t141 * t189 - t143 * t167 - t144 * t168 + t264 * t54 + t367 * t53) * mrSges(5,3) - t480 * t83 + t34 * (mrSges(7,1) * t198 - mrSges(7,2) * t199) + (-t19 * t62 - t198 * t2 + t199 * t3 + t20 * t63) * mrSges(7,3) + (-Ifges(7,5) * t199 - Ifges(7,6) * t198) * t463 + (-Ifges(7,4) * t199 - Ifges(7,2) * t198) * t471 + (-Ifges(7,1) * t199 - Ifges(7,4) * t198) * t472 + (t441 / 0.2e1 + t372 * t452 + t475) * t167 + t368 * Ifges(4,2) * t257 + m(7) * (t106 * t37 + t130 * t34 + t19 * t5 + t2 * t33 + t20 * t4 + t3 * t32) + t209 * t277 + t210 * t278 + t280 * (mrSges(5,1) * t141 + mrSges(5,2) * t140) - t267 * t237 / 0.2e1 - t268 * t236 / 0.2e1 + t248 * t197 + t74 * t228 + t22 * t173 + t21 * t174 + t32 * t30 + t33 * t31 + (Ifges(7,5) * t62 + Ifges(7,6) * t63) * t455 + (Ifges(7,1) * t62 + Ifges(7,4) * t63) * t459 + (Ifges(7,4) * t62 + Ifges(7,2) * t63) * t461 + t62 * t467 + t63 * t469 - t199 * t473 - t198 * t474 + ((t297 / 0.2e1 - pkin(7) * t324 + t383) * t355 + (-t296 / 0.2e1 + pkin(2) * t259 - pkin(7) * t325 + t382) * t352) * qJD(2) + t37 * t89 + (t245 * t386 + t357) * t168 + m(4) * (t202 * t275 + t203 * t274 + t209 * t261 + t210 * t260 + t326 * t343) + t101 * t86 + t102 * t85 + t106 * (-mrSges(7,1) * t63 + mrSges(7,2) * t62) + t4 * t110 + t5 * t111 - t396 * t75 + t130 * t13; ((t277 * t354 - t278 * t351) * qJD(3) + (t256 * t354 - t257 * t351) * mrSges(4,3)) * pkin(2) - m(4) * (t260 * t265 + t261 * t266) + t476 + t288 * t83 - t276 * t197 - t266 * t277 - t265 * t278 + t279 * t13 + (t253 * t348 - t94) * t173 + ((t383 - t297 / 0.2e1 - t341 / 0.2e1 + qJD(1) * t464 + (t324 - t414) * pkin(7)) * t355 + (t382 + t296 / 0.2e1 + (t465 + t440 / 0.2e1 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t355) * qJD(1) + (t325 + t413) * pkin(7) + (-m(4) * t326 - t259) * pkin(2)) * t352) * qJD(1) + (-t253 * t174 - t287 * t86 - t444) * t347 + t206 * t30 + t207 * t31 - t93 * t174 + (t133 * t478 + t253 * t370 + t287 * t371 + t288 * t54 - t87 * t93 - t88 * t94) * m(6) + t479 * t228 + (-t143 * t478 + t144 * t479 - t273 * t276 - t291 * t54 + t292 * t53) * m(5) + t85 * t404 + t484 * t255 + t485 * t110 + t486 * t111 + (t2 * t207 + t206 * t3 + t279 * t34 + t485 * t20 + t486 * t19 + (t255 - t118) * t106) * m(7) + m(4) * (t202 * t351 + t203 * t354 + (-t260 * t351 + t261 * t354) * qJD(3)) * pkin(2) - t118 * t89 + t396 * t164 + (-t140 * t291 - t141 * t292 + t410) * mrSges(5,3); (m(6) * t133 + m(7) * t106 + t484) * pkin(3) * t391 + t490 * t110 + t396 * t155 + t476 + (t143 * t155 - t144 * t158 - t273 * t446 + (-t447 * t54 + t350 * t53 + (-t143 * t350 + t144 * t447) * qJD(4)) * pkin(3)) * m(5) + (-t158 + t379) * t228 + (-t140 * t388 - t141 * t445 + t410) * mrSges(5,3) + (-t106 * t112 + t19 * t489 + t2 * t250 + t20 * t490 + t249 * t3 + t321 * t34) * m(7) + (t332 * t348 - t92) * t173 + (-t133 * t155 + t332 * t370 + t336 * t371 + t338 * t54 - t87 * t91 - t88 * t92) * m(6) + t489 * t111 + (-t332 * t347 - t91) * t174 + t338 * t83 + t321 * t13 - t260 * t277 + t261 * t278 + t249 * t30 + t250 * t31 - t197 * t446 - mrSges(6,3) * t422 - t347 * t336 * t86 + t85 * t403 - t112 * t89; (t19 * t504 - t20 * t497 - t443) * mrSges(7,3) + t501 * t456 + t503 * t460 + t502 * t462 + (-t238 / 0.2e1 - t423 + (t432 / 0.2e1 - t427 / 0.2e1) * t245 + (Ifges(5,1) / 0.2e1 - t386) * t361 + t475) * t245 + (qJD(5) * t348 - t97) * t173 + t500 + (t144 * mrSges(5,3) - t357 + t438 / 0.2e1) * t361 + (-pkin(4) * t54 + qJ(5) * t371 + qJD(5) * t370 - t133 * t144 - t87 * t96 - t88 * t97) * m(6) + t337 * t13 + t269 * t30 + t270 * t31 - t143 * t228 + (-qJ(5) * t86 - qJD(5) * t174 - t444) * t347 - t96 * t174 + t85 * t415 + t487 * t111 + t488 * t110 + (-t106 * t108 + t19 * t487 + t2 * t270 + t20 * t488 + t269 * t3 + t337 * t34) * m(7) - pkin(4) * t83 - t108 * t89 + t396 * t144; -t495 * t110 + t157 * t111 - t380 * t173 + t225 * t174 + t13 + t83 + (t157 * t19 - t20 * t495 + t34) * m(7) + (t225 * t87 - t380 * t88 + t54) * m(6); t137 - t106 * (mrSges(7,1) * t157 + mrSges(7,2) * t495) + (Ifges(7,1) * t495 - t435) * t460 + t78 * t459 + (Ifges(7,5) * t495 - Ifges(7,6) * t157) * t456 - t19 * t110 + t20 * t111 + (t157 * t20 + t19 * t495) * mrSges(7,3) + (-Ifges(7,2) * t157 + t149 + t79) * t462 + t477;];
tauc  = t1(:);
