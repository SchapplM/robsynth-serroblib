% Calculate time derivative of joint inertia matrix for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_inertiaDJ_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRR6_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRR6_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:51
% EndTime: 2019-03-10 04:12:14
% DurationCPUTime: 10.82s
% Computational Cost: add. (24258->818), mult. (61201->1184), div. (0->0), fcn. (61894->12), ass. (0->338)
t352 = sin(qJ(5));
t357 = cos(qJ(5));
t358 = cos(qJ(4));
t485 = (t352 ^ 2 + t357 ^ 2) * t358;
t353 = sin(qJ(4));
t354 = sin(qJ(3));
t359 = cos(qJ(3));
t306 = t353 * t354 - t358 * t359;
t484 = qJD(3) + qJD(4);
t253 = t484 * t306;
t308 = t353 * t359 + t354 * t358;
t401 = qJD(5) * t357;
t367 = -t352 * t253 + t308 * t401;
t402 = qJD(5) * t352;
t411 = t357 * t253;
t366 = t308 * t402 + t411;
t254 = t484 * t308;
t406 = qJD(3) * t354;
t394 = pkin(3) * t406;
t174 = pkin(4) * t254 + pkin(11) * t253 + t394;
t469 = -pkin(10) - pkin(9);
t388 = qJD(3) * t469;
t320 = t354 * t388;
t378 = t359 * t388;
t329 = t469 * t354;
t331 = t469 * t359;
t486 = t358 * t329 + t331 * t353;
t200 = t486 * qJD(4) + t358 * t320 + t353 * t378;
t343 = -pkin(3) * t359 - pkin(2);
t239 = pkin(4) * t306 - pkin(11) * t308 + t343;
t275 = t329 * t353 - t331 * t358;
t75 = t352 * t174 + t357 * t200 + t239 * t401 - t275 * t402;
t266 = t357 * t275;
t180 = t352 * t239 + t266;
t379 = t357 * t174 - t200 * t352;
t76 = -t180 * qJD(5) + t379;
t492 = -t76 * t352 + t357 * t75;
t350 = cos(pkin(6));
t355 = sin(qJ(2));
t349 = sin(pkin(6));
t360 = cos(qJ(2));
t415 = t349 * t360;
t300 = t350 * t355 * pkin(1) + pkin(8) * t415;
t286 = pkin(9) * t350 + t300;
t287 = (-pkin(2) * t360 - pkin(9) * t355 - pkin(1)) * t349;
t215 = -t354 * t286 + t359 * t287;
t416 = t349 * t355;
t295 = t350 * t354 + t359 * t416;
t171 = -pkin(3) * t415 - t295 * pkin(10) + t215;
t216 = t359 * t286 + t354 * t287;
t294 = t350 * t359 - t354 * t416;
t190 = pkin(10) * t294 + t216;
t110 = t353 * t171 + t358 * t190;
t103 = -pkin(11) * t415 + t110;
t226 = t294 * t353 + t295 * t358;
t335 = pkin(8) * t416;
t445 = pkin(1) * t360;
t285 = t335 + (-pkin(2) - t445) * t350;
t233 = -t294 * pkin(3) + t285;
t373 = t358 * t294 - t295 * t353;
t126 = -pkin(4) * t373 - t226 * pkin(11) + t233;
t407 = qJD(2) * t355;
t386 = t349 * t407;
t408 = qJD(2) * t349;
t291 = (pkin(2) * t355 - pkin(9) * t360) * t408;
t299 = t350 * t445 - t335;
t292 = t299 * qJD(2);
t155 = -t216 * qJD(3) + t359 * t291 - t292 * t354;
t385 = t360 * t408;
t265 = qJD(3) * t294 + t359 * t385;
t117 = pkin(3) * t386 - pkin(10) * t265 + t155;
t405 = qJD(3) * t359;
t154 = -t286 * t406 + t287 * t405 + t354 * t291 + t359 * t292;
t264 = -qJD(3) * t295 - t354 * t385;
t123 = pkin(10) * t264 + t154;
t403 = qJD(4) * t358;
t404 = qJD(4) * t353;
t43 = t353 * t117 + t358 * t123 + t171 * t403 - t190 * t404;
t40 = pkin(11) * t386 + t43;
t139 = qJD(4) * t373 + t264 * t353 + t265 * t358;
t140 = qJD(4) * t226 - t358 * t264 + t265 * t353;
t293 = t300 * qJD(2);
t214 = -t264 * pkin(3) + t293;
t67 = t140 * pkin(4) - t139 * pkin(11) + t214;
t13 = -t103 * t402 + t126 * t401 + t352 * t67 + t357 * t40;
t60 = t357 * t103 + t352 * t126;
t14 = -qJD(5) * t60 - t352 * t40 + t357 * t67;
t491 = t13 * t357 - t14 * t352;
t447 = t357 / 0.2e1;
t448 = t352 / 0.2e1;
t351 = sin(qJ(6));
t356 = cos(qJ(6));
t307 = t351 * t357 + t352 * t356;
t455 = t307 / 0.2e1;
t372 = t351 * t352 - t356 * t357;
t456 = -t372 / 0.2e1;
t490 = Ifges(6,5) * t448 + Ifges(7,5) * t455 + Ifges(6,6) * t447 + Ifges(7,6) * t456;
t483 = qJD(5) + qJD(6);
t251 = t483 * t372;
t252 = t483 * t307;
t185 = -Ifges(7,5) * t251 - Ifges(7,6) * t252;
t344 = Ifges(6,5) * t401;
t382 = -t402 / 0.2e1;
t489 = t185 / 0.2e1 + Ifges(6,6) * t382 + t344 / 0.2e1;
t230 = t372 * t308;
t340 = pkin(3) * t353 + pkin(11);
t443 = -pkin(12) - t340;
t303 = t443 * t352;
t346 = t357 * pkin(12);
t304 = t340 * t357 + t346;
t234 = t303 * t356 - t304 * t351;
t380 = qJD(5) * t443;
t393 = pkin(3) * t403;
t276 = t352 * t380 + t357 * t393;
t277 = -t352 * t393 + t357 * t380;
t149 = qJD(6) * t234 + t276 * t356 + t277 * t351;
t235 = t303 * t351 + t304 * t356;
t150 = -qJD(6) * t235 - t276 * t351 + t277 * t356;
t488 = t150 * mrSges(7,1) - t149 * mrSges(7,2);
t468 = -pkin(12) - pkin(11);
t328 = t468 * t352;
t330 = pkin(11) * t357 + t346;
t272 = t328 * t356 - t330 * t351;
t387 = qJD(5) * t468;
t317 = t352 * t387;
t318 = t357 * t387;
t198 = qJD(6) * t272 + t317 * t356 + t318 * t351;
t274 = t328 * t351 + t330 * t356;
t199 = -qJD(6) * t274 - t317 * t351 + t318 * t356;
t487 = t199 * mrSges(7,1) - t198 * mrSges(7,2);
t482 = Ifges(4,5) * t265 + Ifges(4,6) * t264 + Ifges(4,3) * t386;
t481 = 2 * m(5);
t480 = 2 * m(6);
t479 = 2 * m(7);
t478 = -2 * mrSges(3,3);
t477 = -2 * mrSges(5,3);
t183 = mrSges(7,1) * t252 - mrSges(7,2) * t251;
t476 = 0.2e1 * t183;
t201 = qJD(4) * t275 + t320 * t353 - t358 * t378;
t475 = 0.2e1 * t201;
t257 = mrSges(7,1) * t372 + mrSges(7,2) * t307;
t474 = 0.2e1 * t257;
t473 = -0.2e1 * t486;
t472 = 0.2e1 * t293;
t471 = m(6) * pkin(4);
t368 = -t357 * t226 + t352 * t415;
t91 = qJD(5) * t368 - t352 * t139 + t357 * t386;
t470 = t91 / 0.2e1;
t207 = -t352 * t226 - t357 * t415;
t130 = t207 * t356 + t351 * t368;
t467 = t130 / 0.2e1;
t131 = t207 * t351 - t356 * t368;
t466 = t131 / 0.2e1;
t464 = t207 / 0.2e1;
t229 = t307 * t308;
t463 = -t229 / 0.2e1;
t462 = -t230 / 0.2e1;
t461 = -t251 / 0.2e1;
t460 = -t252 / 0.2e1;
t260 = Ifges(7,4) * t307 - Ifges(7,2) * t372;
t458 = t260 / 0.2e1;
t262 = Ifges(7,1) * t307 - Ifges(7,4) * t372;
t457 = t262 / 0.2e1;
t438 = Ifges(6,4) * t352;
t375 = Ifges(6,1) * t357 - t438;
t315 = t375 * qJD(5);
t453 = t315 / 0.2e1;
t437 = Ifges(6,4) * t357;
t326 = Ifges(6,1) * t352 + t437;
t450 = t326 / 0.2e1;
t449 = -t352 / 0.2e1;
t444 = pkin(3) * t358;
t341 = -pkin(4) - t444;
t446 = m(6) * t341;
t442 = mrSges(7,3) * t252;
t441 = mrSges(7,3) * t372;
t440 = Ifges(4,4) * t354;
t439 = Ifges(4,4) * t359;
t436 = Ifges(6,6) * t352;
t435 = pkin(3) * qJD(4);
t434 = pkin(5) * qJD(6);
t429 = t292 * mrSges(3,2);
t428 = t293 * mrSges(3,1);
t427 = t353 * mrSges(5,1);
t425 = t358 * mrSges(5,2);
t423 = t201 * t486;
t422 = t486 * t353;
t421 = t308 * t352;
t420 = t308 * t357;
t413 = t352 * t358;
t322 = -mrSges(6,1) * t357 + mrSges(6,2) * t352;
t412 = t353 * t322;
t410 = t357 * t358;
t133 = -mrSges(6,1) * t207 - mrSges(6,2) * t368;
t213 = -mrSges(5,1) * t415 - t226 * mrSges(5,3);
t409 = t133 - t213;
t400 = qJD(6) * t351;
t399 = qJD(6) * t356;
t398 = 0.2e1 * mrSges(7,3);
t397 = 0.2e1 * t349;
t396 = mrSges(7,3) * t434;
t90 = qJD(5) * t207 + t357 * t139 + t352 * t386;
t31 = qJD(6) * t130 + t351 * t91 + t356 * t90;
t32 = -qJD(6) * t131 - t351 * t90 + t356 * t91;
t9 = Ifges(7,5) * t31 + Ifges(7,6) * t32 + Ifges(7,3) * t140;
t35 = Ifges(6,5) * t90 + Ifges(6,6) * t91 + Ifges(6,3) * t140;
t100 = -t252 * t308 + t253 * t372;
t101 = t483 * t230 + t307 * t253;
t49 = Ifges(7,5) * t100 + Ifges(7,6) * t101 + Ifges(7,3) * t254;
t392 = pkin(5) * t402;
t391 = Ifges(4,6) * t415;
t390 = t356 * t251 * mrSges(7,3);
t389 = Ifges(5,5) * t139 - Ifges(5,6) * t140 + Ifges(5,3) * t386;
t342 = -pkin(5) * t357 - pkin(4);
t381 = t401 / 0.2e1;
t59 = -t103 * t352 + t357 * t126;
t109 = t171 * t358 - t353 * t190;
t179 = t357 * t239 - t275 * t352;
t377 = mrSges(6,3) * t485;
t102 = pkin(4) * t415 - t109;
t376 = mrSges(6,1) * t352 + mrSges(6,2) * t357;
t374 = -Ifges(6,2) * t352 + t437;
t45 = -pkin(5) * t373 + pkin(12) * t368 + t59;
t52 = pkin(12) * t207 + t60;
t21 = -t351 * t52 + t356 * t45;
t22 = t351 * t45 + t356 * t52;
t145 = pkin(5) * t306 - pkin(12) * t420 + t179;
t160 = -pkin(12) * t421 + t180;
t78 = t145 * t356 - t160 * t351;
t79 = t145 * t351 + t160 * t356;
t44 = t117 * t358 - t353 * t123 - t171 * t404 - t190 * t403;
t5 = pkin(5) * t140 - pkin(12) * t90 + t14;
t6 = pkin(12) * t91 + t13;
t3 = qJD(6) * t21 + t351 * t5 + t356 * t6;
t4 = -qJD(6) * t22 - t351 * t6 + t356 * t5;
t371 = t4 * mrSges(7,1) - t3 * mrSges(7,2) + t9;
t56 = pkin(12) * t411 + pkin(5) * t254 + (-t266 + (pkin(12) * t308 - t239) * t352) * qJD(5) + t379;
t68 = -pkin(12) * t367 + t75;
t18 = qJD(6) * t78 + t351 * t56 + t356 * t68;
t19 = -qJD(6) * t79 - t351 * t68 + t356 * t56;
t370 = t19 * mrSges(7,1) - t18 * mrSges(7,2) + t49;
t369 = -t356 * t372 * t396 + t185 + t344 + (-pkin(5) * t442 + t307 * t396) * t351;
t41 = -pkin(4) * t386 - t44;
t186 = -Ifges(7,4) * t251 - Ifges(7,2) * t252;
t188 = -Ifges(7,1) * t251 - Ifges(7,4) * t252;
t313 = t374 * qJD(5);
t324 = Ifges(6,2) * t357 + t438;
t365 = -t186 * t372 + t307 * t188 - t251 * t262 - t252 * t260 + t357 * t313 + t352 * t315 - t324 * t402 + t326 * t401;
t112 = -Ifges(6,5) * t366 - Ifges(6,6) * t367 + Ifges(6,3) * t254;
t147 = mrSges(6,2) * t373 + mrSges(6,3) * t207;
t148 = -mrSges(6,1) * t373 + mrSges(6,3) * t368;
t57 = mrSges(6,1) * t140 - mrSges(6,3) * t90;
t58 = -mrSges(6,2) * t140 + mrSges(6,3) * t91;
t364 = m(6) * (-t401 * t59 - t402 * t60 + t491) + t357 * t58 - t352 * t57 - t147 * t402 - t148 * t401;
t158 = mrSges(6,1) * t254 + mrSges(6,3) * t366;
t159 = -mrSges(6,2) * t254 - mrSges(6,3) * t367;
t241 = -mrSges(6,2) * t306 - mrSges(6,3) * t421;
t242 = mrSges(6,1) * t306 - mrSges(6,3) * t420;
t363 = m(6) * (-t179 * t401 - t180 * t402 + t492) + t357 * t159 - t352 * t158 - t241 * t402 - t242 * t401;
t10 = Ifges(7,4) * t31 + Ifges(7,2) * t32 + Ifges(7,6) * t140;
t105 = -Ifges(6,4) * t368 + Ifges(6,2) * t207 - Ifges(6,6) * t373;
t106 = -Ifges(6,1) * t368 + Ifges(6,4) * t207 - Ifges(6,5) * t373;
t11 = Ifges(7,1) * t31 + Ifges(7,4) * t32 + Ifges(7,5) * t140;
t24 = -pkin(5) * t91 + t41;
t310 = t376 * qJD(5);
t36 = Ifges(6,4) * t90 + Ifges(6,2) * t91 + Ifges(6,6) * t140;
t37 = Ifges(6,1) * t90 + Ifges(6,4) * t91 + Ifges(6,5) * t140;
t63 = Ifges(7,4) * t131 + Ifges(7,2) * t130 - Ifges(7,6) * t373;
t64 = Ifges(7,1) * t131 + Ifges(7,4) * t130 - Ifges(7,5) * t373;
t77 = -pkin(5) * t207 + t102;
t362 = -t489 * t373 + t490 * t140 + ((-t352 * t60 - t357 * t59) * qJD(5) + t491) * mrSges(6,3) + t106 * t381 + t105 * t382 + (t21 * t251 - t307 * t4) * mrSges(7,3) - t368 * t453 + t313 * t464 + t188 * t466 + t186 * t467 + t324 * t470 + t389 + t41 * t322 + t102 * t310 + t24 * t257 + t77 * t183 + t44 * mrSges(5,1) - t43 * mrSges(5,2) - t3 * t441 - t22 * t442 + t36 * t447 + t37 * t448 + t90 * t450 + t11 * t455 + t10 * t456 + t31 * t457 + t32 * t458 + t63 * t460 + t64 * t461;
t113 = -Ifges(6,4) * t366 - Ifges(6,2) * t367 + Ifges(6,6) * t254;
t114 = -Ifges(6,1) * t366 - Ifges(6,4) * t367 + Ifges(6,5) * t254;
t132 = pkin(5) * t367 + t201;
t152 = -Ifges(7,4) * t230 - Ifges(7,2) * t229 + Ifges(7,6) * t306;
t153 = -Ifges(7,1) * t230 - Ifges(7,4) * t229 + Ifges(7,5) * t306;
t210 = Ifges(6,6) * t306 + t308 * t374;
t211 = Ifges(6,5) * t306 + t308 * t375;
t221 = pkin(5) * t421 - t486;
t246 = Ifges(5,6) * t254;
t248 = Ifges(5,5) * t253;
t50 = Ifges(7,4) * t100 + Ifges(7,2) * t101 + Ifges(7,6) * t254;
t51 = Ifges(7,1) * t100 + Ifges(7,4) * t101 + Ifges(7,5) * t254;
t361 = t489 * t306 + t490 * t254 + ((-t179 * t357 - t180 * t352) * qJD(5) + t492) * mrSges(6,3) - t367 * t324 / 0.2e1 + t211 * t381 + (-t19 * t307 + t251 * t78) * mrSges(7,3) - t313 * t421 / 0.2e1 - t248 + t188 * t462 + t186 * t463 - t486 * t310 + t132 * t257 + t221 * t183 - t200 * mrSges(5,2) + (t308 * t326 + t210) * t382 + (t322 - mrSges(5,1)) * t201 - t18 * t441 - t79 * t442 + t113 * t447 + t114 * t448 - t411 * t450 + t420 * t453 + t51 * t455 + t50 * t456 + t100 * t457 + t101 * t458 + t152 * t460 + t153 * t461 - t246;
t345 = Ifges(4,5) * t405;
t334 = Ifges(3,5) * t385;
t327 = Ifges(4,1) * t354 + t439;
t325 = Ifges(4,2) * t359 + t440;
t321 = t342 - t444;
t319 = pkin(3) * t404 + t392;
t316 = (Ifges(4,1) * t359 - t440) * qJD(3);
t314 = (-Ifges(4,2) * t354 + t439) * qJD(3);
t311 = (mrSges(4,1) * t354 + mrSges(4,2) * t359) * qJD(3);
t302 = (-mrSges(7,1) * t351 - mrSges(7,2) * t356) * t434;
t271 = -mrSges(4,1) * t415 - t295 * mrSges(4,3);
t270 = mrSges(4,2) * t415 + t294 * mrSges(4,3);
t263 = Ifges(5,1) * t308 - Ifges(5,4) * t306;
t261 = Ifges(5,4) * t308 - Ifges(5,2) * t306;
t258 = mrSges(5,1) * t306 + mrSges(5,2) * t308;
t237 = t376 * t308;
t224 = mrSges(4,1) * t386 - mrSges(4,3) * t265;
t223 = -mrSges(4,2) * t386 + mrSges(4,3) * t264;
t220 = Ifges(4,1) * t295 + Ifges(4,4) * t294 - Ifges(4,5) * t415;
t219 = Ifges(4,4) * t295 + Ifges(4,2) * t294 - t391;
t212 = mrSges(5,2) * t415 + mrSges(5,3) * t373;
t209 = Ifges(6,3) * t306 + (Ifges(6,5) * t357 - t436) * t308;
t205 = mrSges(7,1) * t306 + mrSges(7,3) * t230;
t204 = -mrSges(7,2) * t306 - mrSges(7,3) * t229;
t191 = -mrSges(4,1) * t264 + mrSges(4,2) * t265;
t189 = -Ifges(5,1) * t253 - Ifges(5,4) * t254;
t187 = -Ifges(5,4) * t253 - Ifges(5,2) * t254;
t184 = mrSges(5,1) * t254 - mrSges(5,2) * t253;
t173 = Ifges(4,1) * t265 + Ifges(4,4) * t264 + Ifges(4,5) * t386;
t172 = Ifges(4,4) * t265 + Ifges(4,2) * t264 + Ifges(4,6) * t386;
t162 = mrSges(7,1) * t229 - mrSges(7,2) * t230;
t161 = -mrSges(5,1) * t373 + mrSges(5,2) * t226;
t157 = Ifges(5,1) * t226 + Ifges(5,4) * t373 - Ifges(5,5) * t415;
t156 = Ifges(5,4) * t226 + Ifges(5,2) * t373 - Ifges(5,6) * t415;
t151 = -Ifges(7,5) * t230 - Ifges(7,6) * t229 + Ifges(7,3) * t306;
t141 = mrSges(6,1) * t367 - mrSges(6,2) * t366;
t129 = -mrSges(5,2) * t386 - mrSges(5,3) * t140;
t128 = mrSges(5,1) * t386 - mrSges(5,3) * t139;
t104 = -Ifges(6,5) * t368 + Ifges(6,6) * t207 - Ifges(6,3) * t373;
t96 = -mrSges(7,1) * t373 - mrSges(7,3) * t131;
t95 = mrSges(7,2) * t373 + mrSges(7,3) * t130;
t81 = -mrSges(7,2) * t254 + mrSges(7,3) * t101;
t80 = mrSges(7,1) * t254 - mrSges(7,3) * t100;
t73 = mrSges(5,1) * t140 + mrSges(5,2) * t139;
t71 = Ifges(5,1) * t139 - Ifges(5,4) * t140 + Ifges(5,5) * t386;
t70 = Ifges(5,4) * t139 - Ifges(5,2) * t140 + Ifges(5,6) * t386;
t69 = -mrSges(7,1) * t130 + mrSges(7,2) * t131;
t62 = Ifges(7,5) * t131 + Ifges(7,6) * t130 - Ifges(7,3) * t373;
t53 = -mrSges(7,1) * t101 + mrSges(7,2) * t100;
t48 = -mrSges(6,1) * t91 + mrSges(6,2) * t90;
t26 = -mrSges(7,2) * t140 + mrSges(7,3) * t32;
t25 = mrSges(7,1) * t140 - mrSges(7,3) * t31;
t15 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t1 = [(mrSges(3,3) * t355 * t472 + (0.2e1 * mrSges(3,3) * t292 - t389 - t482) * t360 + ((t299 * t478 + Ifges(3,5) * t350 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t360) * t397) * t360 + (t300 * t478 + Ifges(4,5) * t295 + Ifges(5,5) * t226 - 0.2e1 * Ifges(3,6) * t350 + Ifges(4,6) * t294 + Ifges(5,6) * t373 + (-mrSges(3,1) * pkin(1) - Ifges(3,4) * t355) * t397 + ((2 * Ifges(3,1)) - (2 * Ifges(3,2)) - Ifges(4,3) - Ifges(5,3)) * t415) * t355) * qJD(2)) * t349 - (t35 + t9 - t70) * t373 - t368 * t37 + (t334 - 0.2e1 * t428 - 0.2e1 * t429) * t350 + 0.2e1 * m(4) * (t154 * t216 + t155 * t215 + t285 * t293) + 0.2e1 * m(3) * (t292 * t300 - t293 * t299) + (t21 * t4 + t22 * t3 + t24 * t77) * t479 + (t102 * t41 + t13 * t60 + t14 * t59) * t480 + (t109 * t44 + t110 * t43 + t214 * t233) * t481 + (-mrSges(4,1) * t294 + mrSges(4,2) * t295) * t472 + (-t156 + t104 + t62) * t140 + t294 * t172 + t295 * t173 + 0.2e1 * t285 * t191 + 0.2e1 * t154 * t270 + 0.2e1 * t155 * t271 + t264 * t219 + t265 * t220 + 0.2e1 * t233 * t73 + t226 * t71 + 0.2e1 * t216 * t223 + 0.2e1 * t215 * t224 + t207 * t36 + 0.2e1 * t43 * t212 + 0.2e1 * t44 * t213 + 0.2e1 * t214 * t161 + t139 * t157 + 0.2e1 * t13 * t147 + 0.2e1 * t14 * t148 + t130 * t10 + t131 * t11 + 0.2e1 * t41 * t133 + 0.2e1 * t109 * t128 + 0.2e1 * t110 * t129 + t91 * t105 + t90 * t106 + 0.2e1 * t102 * t48 + 0.2e1 * t3 * t95 + 0.2e1 * t4 * t96 + 0.2e1 * t77 * t15 + 0.2e1 * t24 * t69 + 0.2e1 * t59 * t57 + 0.2e1 * t60 * t58 + t32 * t63 + t31 * t64 + 0.2e1 * t21 * t25 + 0.2e1 * t22 * t26; -(t112 / 0.2e1 + t49 / 0.2e1 - t187 / 0.2e1) * t373 - t368 * t114 / 0.2e1 - (t157 / 0.2e1 + t106 * t447 + t105 * t449 - t109 * mrSges(5,3)) * t253 + (t71 / 0.2e1 + t37 * t447 + t36 * t449 - t44 * mrSges(5,3) + (t106 * t449 - t357 * t105 / 0.2e1) * qJD(5)) * t308 + t409 * t201 + ((t220 / 0.2e1 - pkin(9) * t271 - t215 * mrSges(4,3)) * t359 + (t391 / 0.2e1 - t219 / 0.2e1 - pkin(9) * t270 - t216 * mrSges(4,3) + pkin(3) * t161) * t354) * qJD(3) + t334 + m(7) * (t132 * t77 + t18 * t22 + t19 * t21 + t221 * t24 + t3 * t79 + t4 * t78) + ((t154 * t359 - t155 * t354 + (-t215 * t359 - t216 * t354) * qJD(3)) * pkin(9) - pkin(2) * t293) * m(4) + t11 * t462 + t10 * t463 + t113 * t464 + t51 * t466 + t50 * t467 + t210 * t470 + m(5) * (-t109 * t201 + t110 * t200 + t214 * t343 + t233 * t394 + t275 * t43 + t44 * t486) + m(6) * (t102 * t201 + t13 * t180 + t14 * t179 - t41 * t486 + t59 * t76 + t60 * t75) - (t48 - t128) * t486 + ((-t345 / 0.2e1 + t248 / 0.2e1 + t246 / 0.2e1) * t360 + (Ifges(5,5) * t308 / 0.2e1 - Ifges(5,6) * t306 / 0.2e1 - Ifges(3,6) + Ifges(4,5) * t354 / 0.2e1 + Ifges(4,6) * t359 / 0.2e1) * t407) * t349 + t343 * t73 + t264 * t325 / 0.2e1 + t265 * t327 / 0.2e1 + t285 * t311 + t294 * t314 / 0.2e1 + t295 * t316 / 0.2e1 + t275 * t129 + t214 * t258 + t139 * t263 / 0.2e1 + t233 * t184 + t41 * t237 + t13 * t241 + t14 * t242 + t226 * t189 / 0.2e1 + t221 * t15 + t90 * t211 / 0.2e1 + t200 * t212 + t3 * t204 + t4 * t205 - pkin(2) * t191 + t179 * t57 + t180 * t58 + t24 * t162 + t32 * t152 / 0.2e1 + t31 * t153 / 0.2e1 + t59 * t158 + t60 * t159 + t75 * t147 + t76 * t148 + t102 * t141 + t132 * t69 + t100 * t64 / 0.2e1 + t101 * t63 / 0.2e1 + t18 * t95 + t19 * t96 + t79 * t26 + t21 * t80 + t22 * t81 + t77 * t53 + t78 * t25 + (t104 / 0.2e1 + t62 / 0.2e1 - t156 / 0.2e1 - t110 * mrSges(5,3)) * t254 + (-t261 / 0.2e1 + t209 / 0.2e1 + t151 / 0.2e1) * t140 + (t35 / 0.2e1 + t9 / 0.2e1 - t70 / 0.2e1 - t43 * mrSges(5,3)) * t306 + (t293 * mrSges(4,2) + t173 / 0.2e1 - t155 * mrSges(4,3) - pkin(9) * t224) * t354 + (-t293 * mrSges(4,1) + t172 / 0.2e1 + t154 * mrSges(4,3) + pkin(9) * t223) * t359 - t428 - t429; (t179 * t76 + t180 * t75 - t423) * t480 + (t200 * t275 + t343 * t394 - t423) * t481 + (mrSges(5,3) * t475 - t352 * t113 + t357 * t114 + t189 + (-t210 * t357 - t211 * t352) * qJD(5)) * t308 + (t275 * t477 + t151 + t209 - t261) * t254 + (t200 * t477 + t112 - t187 + t49) * t306 - (mrSges(5,3) * t473 - t210 * t352 + t211 * t357 + t263) * t253 + (t132 * t221 + t18 * t79 + t19 * t78) * t479 + t141 * t473 + t237 * t475 + t359 * t314 + t354 * t316 + 0.2e1 * t343 * t184 - 0.2e1 * pkin(2) * t311 - t230 * t51 + 0.2e1 * t75 * t241 + 0.2e1 * t76 * t242 - t229 * t50 + 0.2e1 * t221 * t53 + 0.2e1 * t18 * t204 + 0.2e1 * t19 * t205 + 0.2e1 * t179 * t158 + 0.2e1 * t180 * t159 + 0.2e1 * t132 * t162 + t101 * t152 + t100 * t153 + 0.2e1 * t78 * t80 + 0.2e1 * t79 * t81 + (t359 * t327 + (0.2e1 * pkin(3) * t258 - t325) * t354) * qJD(3); (m(5) * (t353 * t43 + t358 * t44) + t358 * t128 + t353 * t129 + (t409 * t353 + (t147 * t357 - t148 * t352 + t212) * t358 + m(6) * (t102 * t353 + t410 * t60 - t413 * t59) + m(5) * (-t109 * t353 + t110 * t358)) * qJD(4)) * pkin(3) + t364 * t340 + t362 + t341 * t48 + t319 * t69 + t321 * t15 + t234 * t25 + t235 * t26 - t154 * mrSges(4,2) + t155 * mrSges(4,1) + t149 * t95 + t150 * t96 + t41 * t446 + m(7) * (t149 * t22 + t150 * t21 + t234 * t4 + t235 * t3 + t24 * t321 + t319 * t77) + t482; t363 * t340 + (m(5) * (t200 * t353 - t201 * t358) + (t253 * t358 - t254 * t353) * mrSges(5,3) + ((t308 * mrSges(5,3) + t237) * t353 + (-t306 * mrSges(5,3) + t241 * t357 - t242 * t352) * t358 + m(6) * (-t179 * t413 + t180 * t410 - t422) + m(5) * (t275 * t358 - t422)) * qJD(4)) * pkin(3) + m(7) * (t132 * t321 + t149 * t79 + t150 * t78 + t18 * t235 + t19 * t234 + t221 * t319) + t361 + t345 + t341 * t141 + t319 * t162 + t321 * t53 + t234 * t80 + t235 * t81 + t149 * t204 + t150 * t205 + t201 * t446 + (-Ifges(4,6) * t354 + (-mrSges(4,1) * t359 + mrSges(4,2) * t354) * pkin(9)) * qJD(3); t319 * t474 + t321 * t476 + (t149 * t235 + t150 * t234 + t319 * t321) * t479 + 0.2e1 * t341 * t310 + (-t149 * t372 - t150 * t307 + t234 * t251 - t235 * t252) * t398 + (-0.2e1 * t425 - 0.2e1 * t427 + 0.2e1 * t412 + (t485 * t340 + t341 * t353) * t480 + 0.2e1 * t377) * t435 + t365; t364 * pkin(11) + t362 + t342 * t15 + t272 * t25 + t274 * t26 + t198 * t95 + t199 * t96 - pkin(4) * t48 + m(7) * (t198 * t22 + t199 * t21 + t24 * t342 + t272 * t4 + t274 * t3 + t392 * t77) - t41 * t471 + t69 * t392; t363 * pkin(11) + t361 + m(7) * (t132 * t342 + t18 * t274 + t19 * t272 + t198 * t79 + t199 * t78 + t221 * t392) + t342 * t53 + t272 * t80 + t274 * t81 + t198 * t204 + t199 * t205 - pkin(4) * t141 + t162 * t392 - t201 * t471; m(7) * (t149 * t274 + t150 * t272 + t198 * t235 + t199 * t234 + t319 * t342 + t321 * t392) + (t341 - pkin(4)) * t310 + (t319 + t392) * t257 + (t342 + t321) * t183 + (m(6) * (-pkin(4) * t353 + t485 * pkin(11)) - t427 + t412 - t425 + t377) * t435 + ((-t150 - t199) * t307 - (t149 + t198) * t372 - (t235 + t274) * t252 - (-t234 - t272) * t251) * mrSges(7,3) + t365; t392 * t474 + t342 * t476 - 0.2e1 * pkin(4) * t310 + (t198 * t274 + t199 * t272 + t342 * t392) * t479 + (-t198 * t372 - t199 * t307 + t251 * t272 - t252 * t274) * t398 + t365; t14 * mrSges(6,1) - t13 * mrSges(6,2) + (t95 * t399 + t351 * t26 + m(7) * (-t21 * t400 + t22 * t399 + t3 * t351 + t356 * t4) - t96 * t400 + t356 * t25) * pkin(5) + t371 + t35; t76 * mrSges(6,1) - t75 * mrSges(6,2) + (-t205 * t400 + t356 * t80 + m(7) * (t18 * t351 + t19 * t356 + t399 * t79 - t400 * t78) + t204 * t399 + t351 * t81) * pkin(5) + t112 + t370; -t376 * t393 + (t322 * t340 - t436) * qJD(5) + (t390 + m(7) * (t149 * t351 + t150 * t356 - t234 * t400 + t235 * t399)) * pkin(5) + t369 + t488; (pkin(11) * t322 - t436) * qJD(5) + (t390 + m(7) * (t198 * t351 + t199 * t356 - t272 * t400 + t274 * t399)) * pkin(5) + t369 + t487; 0.2e1 * t302; t371; t370; t185 + t488; t185 + t487; t302; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t1(1) t1(2) t1(4) t1(7) t1(11) t1(16); t1(2) t1(3) t1(5) t1(8) t1(12) t1(17); t1(4) t1(5) t1(6) t1(9) t1(13) t1(18); t1(7) t1(8) t1(9) t1(10) t1(14) t1(19); t1(11) t1(12) t1(13) t1(14) t1(15) t1(20); t1(16) t1(17) t1(18) t1(19) t1(20) t1(21);];
Mq  = res;
