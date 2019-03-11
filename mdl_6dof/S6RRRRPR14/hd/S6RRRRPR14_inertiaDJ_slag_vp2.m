% Calculate time derivative of joint inertia matrix for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RRRRPR14_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_inertiaDJ_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR14_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR14_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:49
% EndTime: 2019-03-10 00:09:15
% DurationCPUTime: 12.05s
% Computational Cost: add. (27364->999), mult. (82393->1458), div. (0->0), fcn. (83866->14), ass. (0->398)
t371 = sin(pkin(6));
t498 = 0.2e1 * t371;
t382 = cos(qJ(2));
t374 = cos(pkin(6));
t443 = pkin(1) * t374;
t362 = t382 * t443;
t354 = qJD(2) * t362;
t378 = sin(qJ(2));
t373 = cos(pkin(7));
t394 = t371 * (-pkin(10) * t373 - pkin(9));
t384 = t378 * t394;
t247 = pkin(2) * t374 + t362 + t384;
t428 = t247 * t373;
t497 = qJD(2) * t384 + qJD(3) * t428 + t354;
t369 = sin(pkin(13));
t372 = cos(pkin(13));
t496 = mrSges(6,1) * t369 + mrSges(6,2) * t372;
t495 = 2 * pkin(11);
t375 = sin(qJ(6));
t379 = cos(qJ(6));
t385 = t369 * t375 - t372 * t379;
t455 = -t385 / 0.2e1;
t323 = t369 * t379 + t372 * t375;
t454 = t323 / 0.2e1;
t447 = t372 / 0.2e1;
t381 = cos(qJ(3));
t416 = t378 * t381;
t377 = sin(qJ(3));
t417 = t377 * t382;
t383 = t373 * t417 + t416;
t370 = sin(pkin(7));
t411 = qJD(3) * t377;
t401 = t370 * t411;
t178 = t374 * t401 + (t383 * qJD(3) + (t373 * t416 + t417) * qJD(2)) * t371;
t442 = pkin(10) * t370;
t274 = (-pkin(2) * t382 - t378 * t442 - pkin(1)) * t371;
t174 = -t247 * t370 + t373 * t274;
t424 = t370 * t381;
t415 = t381 * t382;
t418 = t377 * t378;
t492 = t373 * t415 - t418;
t234 = -t371 * t492 - t374 * t424;
t425 = t370 * t377;
t235 = t371 * t383 + t374 * t425;
t121 = pkin(3) * t234 - pkin(11) * t235 + t174;
t361 = t378 * t443;
t423 = t371 * t382;
t318 = pkin(9) * t423 + t361;
t231 = (t370 * t374 + t373 * t423) * pkin(10) + t318;
t420 = t373 * t377;
t137 = t381 * t231 + t247 * t420 + t274 * t425;
t307 = -t370 * t423 + t374 * t373;
t130 = pkin(11) * t307 + t137;
t376 = sin(qJ(4));
t380 = cos(qJ(4));
t408 = qJD(4) * t380;
t409 = qJD(4) * t376;
t413 = qJD(2) * t371;
t403 = t378 * t413;
t397 = t370 * t403;
t251 = (t382 * t394 - t361) * qJD(2);
t280 = (pkin(2) * t378 - t382 * t442) * t413;
t410 = qJD(3) * t381;
t400 = t370 * t410;
t77 = -t231 * t411 + t251 * t420 + t274 * t400 + t280 * t425 + t381 * t497;
t75 = pkin(11) * t397 + t77;
t179 = t374 * t400 + (t492 * qJD(3) + (-t373 * t418 + t415) * qJD(2)) * t371;
t182 = -t251 * t370 + t373 * t280;
t90 = pkin(3) * t178 - pkin(11) * t179 + t182;
t23 = -t121 * t409 - t130 * t408 - t376 * t75 + t380 * t90;
t20 = -pkin(4) * t178 - t23;
t185 = t235 * t376 - t307 * t380;
t117 = -qJD(4) * t185 + t179 * t380 + t376 * t397;
t80 = -t117 * t369 + t178 * t372;
t81 = t117 * t372 + t178 * t369;
t41 = -t80 * mrSges(6,1) + t81 * mrSges(6,2);
t494 = -m(6) * t20 - t41;
t317 = pkin(2) * t420 + pkin(10) * t424;
t290 = pkin(11) * t373 + t317;
t291 = (-pkin(3) * t381 - pkin(11) * t377 - pkin(2)) * t370;
t412 = qJD(3) * t370;
t299 = (pkin(3) * t377 - pkin(11) * t381) * t412;
t355 = pkin(10) * t425;
t419 = t373 * t381;
t315 = pkin(2) * t419 - t355;
t300 = t315 * qJD(3);
t149 = -t290 * t408 - t291 * t409 + t299 * t380 - t376 * t300;
t140 = -pkin(4) * t401 - t149;
t310 = -t380 * t373 + t376 * t425;
t252 = -qJD(4) * t310 + t380 * t400;
t204 = -t252 * t369 + t372 * t401;
t205 = t252 * t372 + t369 * t401;
t147 = -t204 * mrSges(6,1) + t205 * mrSges(6,2);
t493 = -m(6) * t140 - t147;
t308 = t385 * qJD(6);
t136 = -t377 * t231 + t381 * (t274 * t370 + t428);
t491 = 2 * m(4);
t490 = 2 * m(5);
t489 = 0.2e1 * m(6);
t488 = 2 * m(7);
t487 = -2 * mrSges(3,3);
t486 = -2 * mrSges(4,3);
t186 = t235 * t380 + t307 * t376;
t116 = qJD(4) * t186 + t179 * t376 - t380 * t397;
t34 = Ifges(6,4) * t81 + Ifges(6,2) * t80 + Ifges(6,6) * t116;
t485 = t34 / 0.2e1;
t35 = Ifges(6,1) * t81 + Ifges(6,4) * t80 + Ifges(6,5) * t116;
t484 = t35 / 0.2e1;
t52 = Ifges(5,1) * t117 - Ifges(5,4) * t116 + Ifges(5,5) * t178;
t483 = t52 / 0.2e1;
t145 = -t186 * t369 + t234 * t372;
t146 = t186 * t372 + t234 * t369;
t86 = t145 * t379 - t146 * t375;
t482 = t86 / 0.2e1;
t87 = t145 * t375 + t146 * t379;
t481 = t87 / 0.2e1;
t107 = Ifges(5,1) * t186 - Ifges(5,4) * t185 + Ifges(5,5) * t234;
t479 = t107 / 0.2e1;
t311 = t373 * t376 + t380 * t425;
t253 = qJD(4) * t311 + t376 * t400;
t124 = Ifges(6,4) * t205 + Ifges(6,2) * t204 + Ifges(6,6) * t253;
t478 = t124 / 0.2e1;
t125 = Ifges(6,1) * t205 + Ifges(6,4) * t204 + Ifges(6,5) * t253;
t477 = t125 / 0.2e1;
t167 = Ifges(5,1) * t252 - Ifges(5,4) * t253 + Ifges(5,5) * t401;
t476 = t167 / 0.2e1;
t248 = -t311 * t369 - t372 * t424;
t249 = t311 * t372 - t369 * t424;
t170 = t248 * t379 - t249 * t375;
t475 = t170 / 0.2e1;
t171 = t248 * t375 + t249 * t379;
t474 = t171 / 0.2e1;
t292 = t323 * t376;
t293 = t385 * t376;
t201 = -Ifges(7,4) * t293 - Ifges(7,2) * t292 - Ifges(7,6) * t380;
t473 = t201 / 0.2e1;
t202 = -Ifges(7,1) * t293 - Ifges(7,4) * t292 - Ifges(7,5) * t380;
t472 = t202 / 0.2e1;
t211 = Ifges(5,1) * t311 - Ifges(5,4) * t310 - Ifges(5,5) * t424;
t471 = t211 / 0.2e1;
t309 = t323 * qJD(6);
t212 = -t309 * t376 - t385 * t408;
t470 = t212 / 0.2e1;
t213 = t308 * t376 - t323 * t408;
t469 = t213 / 0.2e1;
t222 = -Ifges(7,5) * t308 - Ifges(7,6) * t309;
t468 = t222 / 0.2e1;
t223 = -Ifges(7,4) * t308 - Ifges(7,2) * t309;
t467 = t223 / 0.2e1;
t224 = -Ifges(7,1) * t308 - Ifges(7,4) * t309;
t466 = t224 / 0.2e1;
t238 = Ifges(7,4) * t323 - Ifges(7,2) * t385;
t465 = t238 / 0.2e1;
t239 = Ifges(7,1) * t323 - Ifges(7,4) * t385;
t464 = t239 / 0.2e1;
t431 = Ifges(6,4) * t372;
t388 = -Ifges(6,2) * t369 + t431;
t269 = (Ifges(6,6) * t376 + t380 * t388) * qJD(4);
t463 = t269 / 0.2e1;
t432 = Ifges(6,4) * t369;
t389 = Ifges(6,1) * t372 - t432;
t270 = (Ifges(6,5) * t376 + t380 * t389) * qJD(4);
t462 = t270 / 0.2e1;
t287 = -Ifges(6,6) * t380 + t376 * t388;
t461 = t287 / 0.2e1;
t288 = -Ifges(6,5) * t380 + t376 * t389;
t460 = t288 / 0.2e1;
t459 = -t292 / 0.2e1;
t458 = -t293 / 0.2e1;
t457 = -t308 / 0.2e1;
t456 = -t309 / 0.2e1;
t434 = Ifges(5,4) * t376;
t332 = (Ifges(5,1) * t380 - t434) * qJD(4);
t453 = t332 / 0.2e1;
t341 = Ifges(6,2) * t372 + t432;
t452 = t341 / 0.2e1;
t342 = Ifges(6,1) * t369 + t431;
t451 = t342 / 0.2e1;
t450 = Ifges(5,5) * t376 / 0.2e1 + Ifges(5,6) * t380 / 0.2e1;
t433 = Ifges(5,4) * t380;
t346 = Ifges(5,1) * t376 + t433;
t449 = t346 / 0.2e1;
t448 = -t369 / 0.2e1;
t446 = t373 / 0.2e1;
t444 = m(6) * t380;
t441 = pkin(11) * t376;
t440 = pkin(11) * t380;
t439 = pkin(12) + qJ(5);
t22 = t121 * t408 - t130 * t409 + t376 * t90 + t380 * t75;
t17 = qJ(5) * t178 + qJD(5) * t234 + t22;
t398 = -t231 * t410 - t274 * t401 - t377 * t497;
t76 = -t251 * t419 + (-pkin(3) * t403 - t280 * t381) * t370 - t398;
t32 = pkin(4) * t116 - qJ(5) * t117 - qJD(5) * t186 + t76;
t11 = t372 * t17 + t369 * t32;
t64 = t376 * t121 + t380 * t130;
t57 = qJ(5) * t234 + t64;
t129 = -pkin(3) * t307 - t136;
t70 = pkin(4) * t185 - qJ(5) * t186 + t129;
t37 = t369 * t70 + t372 * t57;
t436 = Ifges(4,4) * t377;
t435 = Ifges(4,4) * t381;
t302 = -pkin(9) * t403 + t354;
t430 = t302 * mrSges(3,2);
t303 = t318 * qJD(2);
t429 = t303 * mrSges(3,1);
t427 = t369 * t376;
t426 = t369 * t380;
t422 = t372 * t376;
t421 = t372 * t380;
t148 = -t290 * t409 + t291 * t408 + t376 * t299 + t380 * t300;
t135 = (qJ(5) * t411 - qJD(5) * t381) * t370 + t148;
t301 = t317 * qJD(3);
t141 = pkin(4) * t253 - qJ(5) * t252 - qJD(5) * t311 + t301;
t74 = t372 * t135 + t369 * t141;
t289 = t355 + (-pkin(2) * t381 - pkin(3)) * t373;
t184 = pkin(4) * t310 - qJ(5) * t311 + t289;
t196 = t380 * t290 + t376 * t291;
t187 = -qJ(5) * t424 + t196;
t132 = t369 * t184 + t372 * t187;
t295 = t496 * t408;
t306 = -qJD(5) * t376 + (pkin(4) * t376 - qJ(5) * t380) * qJD(4);
t406 = pkin(11) * t409;
t245 = t372 * t306 + t369 * t406;
t336 = -pkin(4) * t380 - qJ(5) * t376 - pkin(3);
t276 = pkin(11) * t421 + t369 * t336;
t28 = qJD(6) * t86 + t375 * t80 + t379 * t81;
t29 = -qJD(6) * t87 - t375 * t81 + t379 * t80;
t5 = Ifges(7,5) * t28 + Ifges(7,6) * t29 + Ifges(7,3) * t116;
t95 = qJD(6) * t170 + t204 * t375 + t205 * t379;
t96 = -qJD(6) * t171 + t204 * t379 - t205 * t375;
t44 = Ifges(7,5) * t95 + Ifges(7,6) * t96 + Ifges(7,3) * t253;
t50 = Ifges(5,5) * t117 - Ifges(5,6) * t116 + Ifges(5,3) * t178;
t113 = Ifges(4,5) * t179 - Ifges(4,6) * t178 + Ifges(4,3) * t397;
t142 = Ifges(7,5) * t212 + Ifges(7,6) * t213 + Ifges(7,3) * t409;
t165 = Ifges(5,5) * t252 - Ifges(5,6) * t253 + Ifges(5,3) * t401;
t404 = pkin(5) * t369 + pkin(11);
t399 = Ifges(6,5) * t369 / 0.2e1 + Ifges(6,6) * t447 + Ifges(7,5) * t454 + Ifges(7,6) * t455;
t12 = -t29 * mrSges(7,1) + t28 * mrSges(7,2);
t49 = -t96 * mrSges(7,1) + t95 * mrSges(7,2);
t10 = -t17 * t369 + t372 * t32;
t36 = -t369 * t57 + t372 * t70;
t152 = -t213 * mrSges(7,1) + t212 * mrSges(7,2);
t63 = t121 * t380 - t376 * t130;
t73 = -t135 * t369 + t372 * t141;
t131 = t372 * t184 - t187 * t369;
t195 = -t376 * t290 + t291 * t380;
t33 = Ifges(6,5) * t81 + Ifges(6,6) * t80 + Ifges(6,3) * t116;
t51 = Ifges(5,4) * t117 - Ifges(5,2) * t116 + Ifges(5,6) * t178;
t396 = t5 / 0.2e1 - t51 / 0.2e1 + t33 / 0.2e1;
t106 = Ifges(5,4) * t186 - Ifges(5,2) * t185 + Ifges(5,6) * t234;
t38 = Ifges(7,5) * t87 + Ifges(7,6) * t86 + Ifges(7,3) * t185;
t65 = Ifges(6,5) * t146 + Ifges(6,6) * t145 + Ifges(6,3) * t185;
t395 = t65 / 0.2e1 + t38 / 0.2e1 - t106 / 0.2e1;
t123 = Ifges(6,5) * t205 + Ifges(6,6) * t204 + Ifges(6,3) * t253;
t166 = Ifges(5,4) * t252 - Ifges(5,2) * t253 + Ifges(5,6) * t401;
t393 = t123 / 0.2e1 + t44 / 0.2e1 - t166 / 0.2e1;
t158 = Ifges(6,5) * t249 + Ifges(6,6) * t248 + Ifges(6,3) * t310;
t210 = Ifges(5,4) * t311 - Ifges(5,2) * t310 - Ifges(5,6) * t424;
t99 = Ifges(7,5) * t171 + Ifges(7,6) * t170 + Ifges(7,3) * t310;
t392 = -t210 / 0.2e1 + t158 / 0.2e1 + t99 / 0.2e1;
t188 = pkin(4) * t424 - t195;
t387 = Ifges(6,5) * t372 - Ifges(6,6) * t369;
t268 = (Ifges(6,3) * t376 + t380 * t387) * qJD(4);
t331 = (-Ifges(5,2) * t376 + t433) * qJD(4);
t391 = -t331 / 0.2e1 + t268 / 0.2e1 + t142 / 0.2e1;
t200 = -Ifges(7,5) * t293 - Ifges(7,6) * t292 - Ifges(7,3) * t380;
t286 = -Ifges(6,3) * t380 + t376 * t387;
t345 = Ifges(5,2) * t380 + t434;
t390 = -t345 / 0.2e1 + t286 / 0.2e1 + t200 / 0.2e1;
t21 = pkin(5) * t185 - pkin(12) * t146 + t36;
t24 = pkin(12) * t145 + t37;
t8 = t21 * t379 - t24 * t375;
t9 = t21 * t375 + t24 * t379;
t104 = pkin(12) * t248 + t132;
t94 = pkin(5) * t310 - pkin(12) * t249 + t131;
t54 = t104 * t379 + t375 * t94;
t53 = -t104 * t375 + t379 * t94;
t321 = t372 * t336;
t230 = -pkin(12) * t422 + t321 + (-pkin(11) * t369 - pkin(5)) * t380;
t254 = -pkin(12) * t427 + t276;
t168 = t230 * t379 - t254 * t375;
t169 = t230 * t375 + t254 * t379;
t337 = t439 * t369;
t339 = t439 * t372;
t256 = -t337 * t379 - t339 * t375;
t257 = -t337 * t375 + t339 * t379;
t58 = -pkin(4) * t234 - t63;
t366 = Ifges(5,5) * t408;
t364 = -pkin(5) * t372 - pkin(4);
t352 = Ifges(3,5) * t382 * t413;
t351 = Ifges(4,5) * t400;
t343 = -mrSges(5,1) * t380 + mrSges(5,2) * t376;
t338 = -mrSges(6,1) * t372 + mrSges(6,2) * t369;
t333 = t404 * t376;
t330 = -Ifges(5,6) * t409 + t366;
t329 = (mrSges(5,1) * t376 + mrSges(5,2) * t380) * qJD(4);
t328 = -mrSges(6,1) * t380 - mrSges(6,3) * t422;
t327 = mrSges(6,2) * t380 - mrSges(6,3) * t427;
t326 = -mrSges(4,2) * t373 + mrSges(4,3) * t424;
t325 = mrSges(4,1) * t373 - mrSges(4,3) * t425;
t319 = t404 * t408;
t316 = -pkin(9) * t371 * t378 + t362;
t314 = (mrSges(6,1) * t376 - mrSges(6,3) * t421) * qJD(4);
t313 = (-mrSges(6,2) * t376 - mrSges(6,3) * t426) * qJD(4);
t312 = t496 * t376;
t298 = (Ifges(4,1) * t381 - t436) * t412;
t297 = (-Ifges(4,2) * t377 + t435) * t412;
t296 = -Ifges(4,6) * t401 + t351;
t294 = (mrSges(4,1) * t377 + mrSges(4,2) * t381) * t412;
t284 = Ifges(4,5) * t373 + (Ifges(4,1) * t377 + t435) * t370;
t283 = Ifges(4,6) * t373 + (Ifges(4,2) * t381 + t436) * t370;
t281 = t369 * t306;
t275 = -pkin(11) * t426 + t321;
t263 = -mrSges(7,1) * t380 + mrSges(7,3) * t293;
t262 = mrSges(7,2) * t380 - mrSges(7,3) * t292;
t261 = -mrSges(5,1) * t424 - mrSges(5,3) * t311;
t260 = mrSges(5,2) * t424 - mrSges(5,3) * t310;
t246 = -t372 * t406 + t281;
t236 = mrSges(7,1) * t385 + mrSges(7,2) * t323;
t225 = mrSges(5,1) * t310 + mrSges(5,2) * t311;
t221 = mrSges(7,1) * t309 - mrSges(7,2) * t308;
t217 = t281 + (-pkin(11) * t422 - pkin(12) * t426) * qJD(4);
t216 = -mrSges(5,2) * t401 - mrSges(5,3) * t253;
t215 = mrSges(5,1) * t401 - mrSges(5,3) * t252;
t214 = mrSges(7,1) * t292 - mrSges(7,2) * t293;
t209 = Ifges(5,5) * t311 - Ifges(5,6) * t310 - Ifges(5,3) * t424;
t203 = (pkin(5) * t376 - pkin(12) * t421) * qJD(4) + t245;
t198 = -qJD(5) * t323 - qJD(6) * t257;
t197 = -qJD(5) * t385 + qJD(6) * t256;
t194 = mrSges(6,1) * t310 - mrSges(6,3) * t249;
t193 = -mrSges(6,2) * t310 + mrSges(6,3) * t248;
t192 = -mrSges(7,2) * t409 + mrSges(7,3) * t213;
t191 = mrSges(7,1) * t409 - mrSges(7,3) * t212;
t190 = mrSges(4,1) * t307 - mrSges(4,3) * t235;
t189 = -mrSges(4,2) * t307 - mrSges(4,3) * t234;
t173 = mrSges(5,1) * t253 + mrSges(5,2) * t252;
t172 = -mrSges(6,1) * t248 + mrSges(6,2) * t249;
t164 = mrSges(4,1) * t397 - mrSges(4,3) * t179;
t163 = -mrSges(4,2) * t397 - mrSges(4,3) * t178;
t162 = mrSges(6,1) * t253 - mrSges(6,3) * t205;
t161 = -mrSges(6,2) * t253 + mrSges(6,3) * t204;
t160 = Ifges(6,1) * t249 + Ifges(6,4) * t248 + Ifges(6,5) * t310;
t159 = Ifges(6,4) * t249 + Ifges(6,2) * t248 + Ifges(6,6) * t310;
t157 = Ifges(4,1) * t235 - Ifges(4,4) * t234 + Ifges(4,5) * t307;
t156 = Ifges(4,4) * t235 - Ifges(4,2) * t234 + Ifges(4,6) * t307;
t155 = -pkin(5) * t248 + t188;
t154 = mrSges(7,1) * t310 - mrSges(7,3) * t171;
t153 = -mrSges(7,2) * t310 + mrSges(7,3) * t170;
t151 = mrSges(5,1) * t234 - mrSges(5,3) * t186;
t150 = -mrSges(5,2) * t234 - mrSges(5,3) * t185;
t144 = Ifges(7,1) * t212 + Ifges(7,4) * t213 + Ifges(7,5) * t409;
t143 = Ifges(7,4) * t212 + Ifges(7,2) * t213 + Ifges(7,6) * t409;
t133 = mrSges(5,1) * t185 + mrSges(5,2) * t186;
t128 = mrSges(4,1) * t178 + mrSges(4,2) * t179;
t115 = Ifges(4,1) * t179 - Ifges(4,4) * t178 + Ifges(4,5) * t397;
t114 = Ifges(4,4) * t179 - Ifges(4,2) * t178 + Ifges(4,6) * t397;
t109 = -mrSges(7,1) * t170 + mrSges(7,2) * t171;
t108 = -pkin(5) * t204 + t140;
t105 = Ifges(5,5) * t186 - Ifges(5,6) * t185 + Ifges(5,3) * t234;
t103 = mrSges(6,1) * t185 - mrSges(6,3) * t146;
t102 = -mrSges(6,2) * t185 + mrSges(6,3) * t145;
t101 = Ifges(7,1) * t171 + Ifges(7,4) * t170 + Ifges(7,5) * t310;
t100 = Ifges(7,4) * t171 + Ifges(7,2) * t170 + Ifges(7,6) * t310;
t98 = -qJD(6) * t169 + t203 * t379 - t217 * t375;
t97 = qJD(6) * t168 + t203 * t375 + t217 * t379;
t88 = -mrSges(6,1) * t145 + mrSges(6,2) * t146;
t85 = mrSges(5,1) * t178 - mrSges(5,3) * t117;
t84 = -mrSges(5,2) * t178 - mrSges(5,3) * t116;
t83 = -mrSges(7,2) * t253 + mrSges(7,3) * t96;
t82 = mrSges(7,1) * t253 - mrSges(7,3) * t95;
t78 = (t251 * t373 + t280 * t370) * t381 + t398;
t67 = Ifges(6,1) * t146 + Ifges(6,4) * t145 + Ifges(6,5) * t185;
t66 = Ifges(6,4) * t146 + Ifges(6,2) * t145 + Ifges(6,6) * t185;
t62 = mrSges(7,1) * t185 - mrSges(7,3) * t87;
t61 = -mrSges(7,2) * t185 + mrSges(7,3) * t86;
t60 = pkin(12) * t204 + t74;
t59 = mrSges(5,1) * t116 + mrSges(5,2) * t117;
t55 = pkin(5) * t253 - pkin(12) * t205 + t73;
t48 = mrSges(6,1) * t116 - mrSges(6,3) * t81;
t47 = -mrSges(6,2) * t116 + mrSges(6,3) * t80;
t46 = Ifges(7,1) * t95 + Ifges(7,4) * t96 + Ifges(7,5) * t253;
t45 = Ifges(7,4) * t95 + Ifges(7,2) * t96 + Ifges(7,6) * t253;
t43 = -pkin(5) * t145 + t58;
t42 = -mrSges(7,1) * t86 + mrSges(7,2) * t87;
t40 = Ifges(7,1) * t87 + Ifges(7,4) * t86 + Ifges(7,5) * t185;
t39 = Ifges(7,4) * t87 + Ifges(7,2) * t86 + Ifges(7,6) * t185;
t19 = -mrSges(7,2) * t116 + mrSges(7,3) * t29;
t18 = mrSges(7,1) * t116 - mrSges(7,3) * t28;
t15 = -qJD(6) * t54 - t375 * t60 + t379 * t55;
t14 = qJD(6) * t53 + t375 * t55 + t379 * t60;
t13 = -pkin(5) * t80 + t20;
t7 = Ifges(7,1) * t28 + Ifges(7,4) * t29 + Ifges(7,5) * t116;
t6 = Ifges(7,4) * t28 + Ifges(7,2) * t29 + Ifges(7,6) * t116;
t4 = pkin(12) * t80 + t11;
t3 = pkin(5) * t116 - pkin(12) * t81 + t10;
t2 = -qJD(6) * t9 + t3 * t379 - t375 * t4;
t1 = qJD(6) * t8 + t3 * t375 + t379 * t4;
t16 = [(t352 - 0.2e1 * t429 - 0.2e1 * t430) * t374 + 0.2e1 * m(3) * (t302 * t318 - t303 * t316) + (t1 * t9 + t13 * t43 + t2 * t8) * t488 + (t10 * t36 + t11 * t37 + t20 * t58) * t489 + (t129 * t76 + t22 * t64 + t23 * t63) * t490 + (t136 * t78 + t137 * t77 + t174 * t182) * t491 + (t33 + t5 - t51) * t185 + (t105 - t156) * t178 + (t65 + t38 - t106) * t116 + t307 * t113 - t234 * t114 + 0.2e1 * t182 * (mrSges(4,1) * t234 + mrSges(4,2) * t235) + t235 * t115 + t234 * t50 + 0.2e1 * t77 * t189 + 0.2e1 * t78 * t190 + t186 * t52 + t179 * t157 + 0.2e1 * t174 * t128 + 0.2e1 * t137 * t163 + 0.2e1 * t136 * t164 + t146 * t35 + 0.2e1 * t22 * t150 + 0.2e1 * t23 * t151 + t145 * t34 + 0.2e1 * t76 * t133 + 0.2e1 * t129 * t59 + t117 * t107 + 0.2e1 * t11 * t102 + 0.2e1 * t10 * t103 + t86 * t6 + t87 * t7 + 0.2e1 * t20 * t88 + t81 * t67 + 0.2e1 * t64 * t84 + 0.2e1 * t63 * t85 + t80 * t66 + 0.2e1 * t1 * t61 + 0.2e1 * t2 * t62 + 0.2e1 * t58 * t41 + 0.2e1 * t37 * t47 + 0.2e1 * t36 * t48 + 0.2e1 * t13 * t42 + 0.2e1 * t43 * t12 + t29 * t39 + t28 * t40 + 0.2e1 * t9 * t19 + 0.2e1 * t8 * t18 + (0.2e1 * (t302 * t382 + t303 * t378) * mrSges(3,3) + ((t316 * t487 + Ifges(3,5) * t374 + (-mrSges(3,2) * pkin(1) + Ifges(3,4) * t382) * t498) * t382 + (t370 * (Ifges(4,5) * t235 - Ifges(4,6) * t234 + Ifges(4,3) * t307) + t318 * t487 - 0.2e1 * Ifges(3,6) * t374 + (-pkin(1) * mrSges(3,1) - Ifges(3,4) * t378 + (Ifges(3,1) - Ifges(3,2)) * t382) * t498) * t378) * qJD(2)) * t371; (t133 - t190) * t301 + m(4) * (-t136 * t301 + t137 * t300 + t315 * t78 + t317 * t77) + t6 * t475 + t186 * t476 + t146 * t477 + t145 * t478 + t252 * t479 + t46 * t481 + t45 * t482 + t311 * t483 + t249 * t484 + t248 * t485 + t117 * t471 + t7 * t474 + t113 * t446 + m(7) * (t1 * t54 + t108 * t43 + t13 * t155 + t14 * t9 + t15 * t8 + t2 * t53) + m(6) * (t10 * t131 + t11 * t132 + t140 * t58 + t188 * t20 + t36 * t73 + t37 * t74) + m(5) * (t129 * t301 + t148 * t64 + t149 * t63 + t195 * t23 + t196 * t22 + t289 * t76) + t352 + (-t297 / 0.2e1 + t165 / 0.2e1) * t234 + t78 * t325 + t77 * t326 + t315 * t164 + t317 * t163 + t307 * t296 / 0.2e1 + t235 * t298 / 0.2e1 + t300 * t189 + t179 * t284 / 0.2e1 + t289 * t59 + t174 * t294 + t22 * t260 + t23 * t261 + t63 * t215 + t64 * t216 + t76 * t225 + t204 * t66 / 0.2e1 + t205 * t67 / 0.2e1 + t196 * t84 + t188 * t41 + t11 * t193 + t10 * t194 + t195 * t85 + t20 * t172 + t129 * t173 + t1 * t153 + t2 * t154 + t155 * t12 + t80 * t159 / 0.2e1 + t81 * t160 / 0.2e1 + t37 * t161 + t36 * t162 + t58 * t147 + t148 * t150 + t149 * t151 + t140 * t88 + t132 * t47 + t131 * t48 + t108 * t42 + t13 * t109 + t29 * t100 / 0.2e1 + t28 * t101 / 0.2e1 + t74 * t102 + t73 * t103 + t95 * t40 / 0.2e1 + t96 * t39 / 0.2e1 + t8 * t82 + t9 * t83 + t14 * t61 + t15 * t62 + t53 * t18 + (-t283 / 0.2e1 + t209 / 0.2e1) * t178 + t54 * t19 + t43 * t49 + ((t182 * mrSges(4,2) + t115 / 0.2e1) * t377 + (-t182 * mrSges(4,1) - t50 / 0.2e1 + t114 / 0.2e1) * t381 + (Ifges(4,3) * t446 + (Ifges(4,5) * t377 + Ifges(4,6) * t381) * t370 / 0.2e1) * t403 + (-m(4) * t182 - t128) * pkin(2) + ((-t136 * mrSges(4,3) + t157 / 0.2e1) * t381 + (-t156 / 0.2e1 + t105 / 0.2e1 - t137 * mrSges(4,3)) * t377) * qJD(3)) * t370 - t429 - t430 + t392 * t116 + t393 * t185 + t395 * t253 + t396 * t310 - Ifges(3,6) * t403; (t108 * t155 + t14 * t54 + t15 * t53) * t488 + (t131 * t73 + t132 * t74 + t140 * t188) * t489 + (t148 * t196 + t149 * t195 + t289 * t301) * t490 + (t300 * t317 - t301 * t315) * t491 + t373 * t296 + 0.2e1 * t300 * t326 + t311 * t167 + 0.2e1 * t289 * t173 + 0.2e1 * t148 * t260 + 0.2e1 * t149 * t261 + t248 * t124 + t249 * t125 + t252 * t211 + 0.2e1 * t195 * t215 + 0.2e1 * t196 * t216 + t204 * t159 + t205 * t160 + 0.2e1 * t188 * t147 + 0.2e1 * t74 * t193 + 0.2e1 * t73 * t194 + t170 * t45 + t171 * t46 + 0.2e1 * t140 * t172 + 0.2e1 * t14 * t153 + 0.2e1 * t15 * t154 + 0.2e1 * t155 * t49 + 0.2e1 * t132 * t161 + 0.2e1 * t131 * t162 + 0.2e1 * t108 * t109 + t96 * t100 + t95 * t101 + 0.2e1 * t53 * t82 + 0.2e1 * t54 * t83 + 0.2e1 * (-t325 + t225) * t301 + (-0.2e1 * pkin(2) * t294 + t377 * t298 + (-t165 + t297) * t381 + ((t315 * t486 + t284) * t381 + (t317 * t486 + t209 - t283) * t377) * qJD(3)) * t370 + (t158 + t99 - t210) * t253 + (-t166 + t123 + t44) * t310; t113 + m(7) * (t1 * t169 + t13 * t333 + t168 * t2 + t319 * t43 + t8 * t98 + t9 * t97) + t144 * t481 + t143 * t482 + t6 * t459 + t81 * t460 + t80 * t461 + t146 * t462 + t145 * t463 + t39 * t469 + t40 * t470 + t28 * t472 + t29 * t473 + t117 * t449 + t178 * t450 + t186 * t453 + t7 * t458 + t333 * t12 + t76 * t343 + t319 * t42 + t11 * t327 + t10 * t328 + t129 * t329 + t234 * t330 / 0.2e1 + t20 * t312 + t37 * t313 + t36 * t314 + t58 * t295 + t1 * t262 + t2 * t263 + t275 * t48 + t276 * t47 + t245 * t103 + t246 * t102 + t13 * t214 + t8 * t191 + t9 * t192 + t168 * t18 + t169 * t19 + t43 * t152 + t97 * t61 + t98 * t62 - t77 * mrSges(4,2) + t78 * mrSges(4,1) - pkin(3) * t59 + (t483 - t23 * mrSges(5,3) + t34 * t448 + t35 * t447 + (t41 - t85) * pkin(11)) * t376 + ((-t63 * mrSges(5,3) + t67 * t447 + t66 * t448 + t479) * t380 + (-t64 * mrSges(5,3) + t395) * t376 + (-t376 * t150 + (-t151 + t88) * t380 + m(5) * (-t376 * t64 - t380 * t63) + t58 * t444) * pkin(11)) * qJD(4) + m(5) * (-pkin(3) * t76 + t22 * t440 - t23 * t441) + m(6) * (t10 * t275 + t11 * t276 + t20 * t441 + t245 * t36 + t246 * t37) + t390 * t116 + t391 * t185 + (t22 * mrSges(5,3) + pkin(11) * t84 - t396) * t380; t45 * t459 + t205 * t460 + t204 * t461 + t249 * t462 + t248 * t463 + t100 * t469 + t101 * t470 + t95 * t472 + t96 * t473 + t144 * t474 + t143 * t475 + t252 * t449 + t311 * t453 + t46 * t458 + t351 + (t343 - mrSges(4,1)) * t301 + t333 * t49 + t319 * t109 + t74 * t327 + t73 * t328 + t289 * t329 + t140 * t312 + t132 * t313 + t131 * t314 + t188 * t295 - t300 * mrSges(4,2) + t14 * t262 + t15 * t263 + t275 * t162 + t276 * t161 + t245 * t194 + t246 * t193 + t108 * t214 + t53 * t191 + t54 * t192 - pkin(3) * t173 + t168 * t82 + t169 * t83 + t97 * t153 + t98 * t154 + t155 * t152 + ((-t195 * mrSges(5,3) + t159 * t448 + t160 * t447 + t471) * t380 + (-t196 * mrSges(5,3) + t392) * t376 + (-t376 * t260 + (t172 - t261) * t380 + t188 * t444 + m(5) * (-t195 * t380 - t196 * t376)) * pkin(11)) * qJD(4) + m(7) * (t108 * t333 + t14 * t169 + t15 * t168 + t155 * t319 + t53 * t98 + t54 * t97) + (t476 - t149 * mrSges(5,3) + t124 * t448 + t125 * t447 + (t147 - t215) * pkin(11)) * t376 + (-t381 * t330 / 0.2e1 + (-Ifges(4,6) + t450) * t411) * t370 + m(5) * (-pkin(3) * t301 + t148 * t440 - t149 * t441) + m(6) * (t131 * t245 + t132 * t246 + t140 * t441 + t275 * t73 + t276 * t74) + t390 * t253 + t391 * t310 + (t148 * mrSges(5,3) + pkin(11) * t216 - t393) * t380; -0.2e1 * pkin(3) * t329 - t292 * t143 - t293 * t144 + 0.2e1 * t333 * t152 + 0.2e1 * t168 * t191 + 0.2e1 * t169 * t192 + t213 * t201 + t212 * t202 + 0.2e1 * t319 * t214 + 0.2e1 * t245 * t328 + 0.2e1 * t246 * t327 + 0.2e1 * t97 * t262 + 0.2e1 * t98 * t263 + 0.2e1 * t275 * t314 + 0.2e1 * t276 * t313 + (t168 * t98 + t169 * t97 + t319 * t333) * t488 + (t245 * t275 + t246 * t276) * t489 + (t331 - t268 - t142) * t380 + (-t269 * t369 + t270 * t372 + t295 * t495 + t332) * t376 + ((t200 + t286 - t345) * t376 + (-t369 * t287 + t372 * t288 + t346 + (m(6) * t441 + t312) * t495) * t380) * qJD(4); (-t1 * t385 - t2 * t323 + t308 * t8 - t309 * t9) * mrSges(7,3) + t50 + t28 * t464 + t29 * t465 + t87 * t466 + t86 * t467 + t185 * t468 + t81 * t451 + t80 * t452 + t7 * t454 + t6 * t455 + t39 * t456 + t40 * t457 + t364 * t12 + t20 * t338 + t256 * t18 + t257 * t19 + t13 * t236 + t43 * t221 + t197 * t61 + t198 * t62 - t22 * mrSges(5,2) + t23 * mrSges(5,1) + (t484 - t10 * mrSges(6,3) - qJ(5) * t48 - qJD(5) * t103 + m(6) * (-qJ(5) * t10 - qJD(5) * t36)) * t369 + (t485 + t11 * mrSges(6,3) + qJ(5) * t47 + qJD(5) * t102 + m(6) * (qJ(5) * t11 + qJD(5) * t37)) * t372 + m(7) * (t1 * t257 + t13 * t364 + t197 * t9 + t198 * t8 + t2 * t256) + t494 * pkin(4) + t399 * t116; (-t14 * t385 - t15 * t323 + t308 * t53 - t309 * t54) * mrSges(7,3) + m(7) * (t108 * t364 + t14 * t257 + t15 * t256 + t197 * t54 + t198 * t53) + t165 + t95 * t464 + t96 * t465 + t171 * t466 + t170 * t467 + t310 * t468 + t205 * t451 + t204 * t452 + t46 * t454 + t45 * t455 + t100 * t456 + t101 * t457 + t364 * t49 + t140 * t338 + t256 * t82 + t257 * t83 + t108 * t236 + t155 * t221 + t197 * t153 + t198 * t154 - t148 * mrSges(5,2) + t149 * mrSges(5,1) + (t477 - t73 * mrSges(6,3) - qJ(5) * t162 - qJD(5) * t194 + m(6) * (-qJ(5) * t73 - qJD(5) * t131)) * t369 + (t478 + t74 * mrSges(6,3) + qJ(5) * t161 + qJD(5) * t193 + m(6) * (qJ(5) * t74 + qJD(5) * t132)) * t372 + t493 * pkin(4) + t399 * t253; t366 - t380 * t222 / 0.2e1 + t364 * t152 + t333 * t221 + t143 * t455 + t144 * t454 + t319 * t236 + t202 * t457 + t201 * t456 - pkin(4) * t295 + t223 * t459 + t224 * t458 + t197 * t262 + t198 * t263 + t256 * t191 + t257 * t192 + t213 * t465 + t212 * t464 + m(7) * (t168 * t198 + t169 * t197 + t256 * t98 + t257 * t97 + t319 * t364) + (t463 + m(6) * (qJ(5) * t246 + qJD(5) * t276) + t246 * mrSges(6,3) + qJ(5) * t313 + qJD(5) * t327) * t372 + (t462 + m(6) * (-qJ(5) * t245 - qJD(5) * t275) - t245 * mrSges(6,3) - qJ(5) * t314 - qJD(5) * t328) * t369 + (t168 * t308 - t169 * t309 - t323 * t98 - t385 * t97) * mrSges(7,3) + ((pkin(11) * mrSges(5,2) - Ifges(5,6) + t399) * t376 + (t342 * t447 + t341 * t448 + (-m(6) * pkin(4) - mrSges(5,1) + t338) * pkin(11)) * t380) * qJD(4); (t197 * t257 + t198 * t256) * t488 - t308 * t239 + t323 * t224 - t309 * t238 - t385 * t223 + 0.2e1 * t364 * t221 + 0.2e1 * (-t197 * t385 - t198 * t323 + t256 * t308 - t257 * t309) * mrSges(7,3) + (qJ(5) * t489 + 0.2e1 * mrSges(6,3)) * qJD(5) * (t369 ^ 2 + t372 ^ 2); m(7) * t13 + t12 - t494; m(7) * t108 + t49 - t493; m(6) * pkin(11) * t408 + m(7) * t319 + t152 + t295; t221; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 + t5; mrSges(7,1) * t15 - mrSges(7,2) * t14 + t44; mrSges(7,1) * t98 - mrSges(7,2) * t97 + t142; mrSges(7,1) * t198 - mrSges(7,2) * t197 + t222; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t16(1) t16(2) t16(4) t16(7) t16(11) t16(16); t16(2) t16(3) t16(5) t16(8) t16(12) t16(17); t16(4) t16(5) t16(6) t16(9) t16(13) t16(18); t16(7) t16(8) t16(9) t16(10) t16(14) t16(19); t16(11) t16(12) t16(13) t16(14) t16(15) t16(20); t16(16) t16(17) t16(18) t16(19) t16(20) t16(21);];
Mq  = res;
