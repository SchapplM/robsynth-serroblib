% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRPRRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:24
% EndTime: 2019-03-08 20:09:35
% DurationCPUTime: 6.69s
% Computational Cost: add. (11167->498), mult. (25741->682), div. (0->0), fcn. (28270->10), ass. (0->251)
t295 = sin(qJ(5));
t290 = t295 ^ 2;
t296 = cos(qJ(5));
t291 = t296 ^ 2;
t509 = t290 + t291;
t476 = -m(7) / 0.2e1;
t292 = sin(pkin(11));
t294 = cos(pkin(11));
t447 = sin(qJ(4));
t449 = cos(qJ(4));
t270 = t292 * t447 - t294 * t449;
t272 = -t292 * t449 - t294 * t447;
t401 = t272 * t295;
t198 = -t270 * mrSges(6,2) + mrSges(6,3) * t401;
t433 = t270 * mrSges(7,3);
t203 = mrSges(7,2) * t401 + t433;
t492 = t198 + t203;
t499 = mrSges(7,2) + mrSges(6,3);
t393 = t292 ^ 2 + t294 ^ 2;
t508 = t393 * mrSges(4,3);
t400 = t272 * t296;
t201 = t270 * mrSges(6,1) + mrSges(6,3) * t400;
t202 = -t270 * mrSges(7,1) - mrSges(7,2) * t400;
t507 = t201 - t202;
t506 = t509 * pkin(9);
t293 = sin(pkin(6));
t448 = sin(qJ(2));
t378 = t293 * t448;
t420 = cos(pkin(6));
t255 = t292 * t420 + t294 * t378;
t311 = t292 * t378 - t294 * t420;
t162 = t255 * t449 - t311 * t447;
t297 = cos(qJ(2));
t399 = t293 * t297;
t125 = t295 * t162 + t296 * t399;
t126 = t162 * t296 - t295 * t399;
t161 = t255 * t447 + t311 * t449;
t337 = -pkin(5) * t296 - qJ(6) * t295;
t188 = t337 * t272;
t363 = t202 / 0.2e1 - t201 / 0.2e1;
t439 = pkin(8) + qJ(3);
t275 = t439 * t294;
t359 = t439 * t292;
t209 = t275 * t449 - t359 * t447;
t396 = t296 * t209;
t285 = -pkin(3) * t294 - pkin(2);
t444 = pkin(9) * t272;
t195 = pkin(4) * t270 + t285 + t444;
t397 = t295 * t195;
t103 = t396 + t397;
t406 = t270 * qJ(6);
t73 = t103 + t406;
t421 = t103 - t73;
t102 = t195 * t296 - t209 * t295;
t74 = -pkin(5) * t270 - t102;
t422 = t102 + t74;
t505 = (t125 * t421 + t126 * t422 + t161 * t188) * t476 - t363 * t126;
t259 = t270 * mrSges(5,2);
t403 = t270 * t296;
t199 = -mrSges(6,1) * t272 + mrSges(6,3) * t403;
t386 = mrSges(7,2) * t403;
t431 = t272 * mrSges(7,1);
t200 = -t386 + t431;
t364 = t200 / 0.2e1 - t199 / 0.2e1;
t504 = (-t272 * mrSges(5,1) - t259) * t399 / 0.2e1 - t364 * t125;
t226 = t270 * t399;
t182 = -t226 * t295 - t296 * t378;
t183 = -t296 * t226 + t295 * t378;
t478 = -m(6) / 0.2e1;
t503 = -(t476 + t478) * (-t182 * t296 + t183 * t295) + (m(4) + m(5)) * t378 / 0.2e1;
t501 = m(6) + m(7);
t336 = pkin(5) * t295 - qJ(6) * t296;
t321 = m(7) * t336;
t500 = mrSges(6,1) + mrSges(7,1);
t498 = Ifges(7,4) + Ifges(6,5);
t497 = Ifges(7,2) + Ifges(6,3);
t496 = Ifges(6,6) - Ifges(7,6);
t404 = t270 * t295;
t197 = mrSges(6,2) * t272 + mrSges(6,3) * t404;
t204 = mrSges(7,2) * t404 - mrSges(7,3) * t272;
t493 = t197 + t204;
t276 = -t296 * mrSges(7,1) - t295 * mrSges(7,3);
t277 = -t296 * mrSges(6,1) + t295 * mrSges(6,2);
t491 = t276 + t277;
t490 = t294 * t255 + t292 * t311;
t424 = t296 * mrSges(7,3);
t428 = t295 * mrSges(7,1);
t346 = -t424 + t428;
t425 = t296 * mrSges(6,2);
t429 = t295 * mrSges(6,1);
t347 = t425 + t429;
t489 = t346 + t347;
t487 = t509 * t444 / 0.2e1;
t445 = pkin(9) * t270;
t446 = pkin(4) * t272;
t206 = t445 - t446;
t208 = t275 * t447 + t449 * t359;
t108 = t206 * t296 + t208 * t295;
t109 = t295 * t206 - t296 * t208;
t486 = -t108 * t295 + t109 * t296;
t80 = -qJ(6) * t272 + t109;
t81 = pkin(5) * t272 - t108;
t485 = t295 * t81 + t296 * t80;
t361 = -t276 / 0.2e1 - t277 / 0.2e1;
t191 = t346 * t272;
t192 = t347 * t272;
t484 = t272 * mrSges(5,3) + t191 + t192;
t286 = m(7) * qJ(6) + mrSges(7,3);
t436 = Ifges(7,5) * t295;
t343 = t296 * Ifges(7,1) + t436;
t134 = t270 * Ifges(7,4) - t272 * t343;
t438 = Ifges(6,4) * t295;
t345 = t296 * Ifges(6,1) - t438;
t136 = t270 * Ifges(6,5) - t272 * t345;
t383 = Ifges(7,4) / 0.2e1 + Ifges(6,5) / 0.2e1;
t483 = -t383 * t270 - t134 / 0.2e1 - t136 / 0.2e1;
t452 = -t296 / 0.2e1;
t454 = t295 / 0.2e1;
t455 = -t295 / 0.2e1;
t482 = t201 * t454 + t202 * t455 + t492 * t452;
t472 = mrSges(7,3) / 0.2e1;
t473 = -mrSges(6,2) / 0.2e1;
t475 = m(7) / 0.2e1;
t481 = (-pkin(5) * t182 + qJ(6) * t183) * t475 + (t473 + t472) * t183;
t480 = 0.2e1 * m(7);
t479 = 2 * qJD(4);
t477 = m(6) / 0.2e1;
t474 = -mrSges(7,1) / 0.2e1;
t471 = t103 / 0.2e1;
t469 = t126 / 0.2e1;
t468 = t161 / 0.2e1;
t466 = -t191 / 0.2e1;
t461 = t204 / 0.2e1;
t459 = -t272 / 0.2e1;
t457 = t276 / 0.2e1;
t451 = t296 / 0.2e1;
t440 = Ifges(7,5) - Ifges(6,4);
t437 = Ifges(6,4) * t296;
t435 = Ifges(7,5) * t296;
t434 = t270 * mrSges(5,3);
t432 = t270 * Ifges(6,6);
t303 = (t295 * t422 - t296 * t421) * t476 + t482;
t315 = t321 / 0.2e1;
t305 = t315 - t424 / 0.2e1 + t425 / 0.2e1 + t429 / 0.2e1 + t428 / 0.2e1;
t304 = t305 * t270;
t12 = t304 + t303;
t417 = t12 * qJD(2);
t416 = t126 * t296;
t225 = t272 * t399;
t107 = t161 * t225;
t415 = t161 * t295;
t394 = t506 * t270;
t273 = -pkin(4) + t337;
t402 = t272 * t273;
t301 = (-t394 + t446) * t477 + (-t394 - t402) * t475 + t499 * t270 * (-t291 / 0.2e1 - t290 / 0.2e1);
t308 = (t108 * t296 + t109 * t295) * t478 + (t295 * t80 - t296 * t81) * t476;
t362 = -t204 / 0.2e1 - t197 / 0.2e1;
t17 = t259 + t364 * t296 + t362 * t295 + (mrSges(5,1) + t361) * t272 + t301 + t308;
t414 = t17 * qJD(2);
t18 = t501 * (-t125 * t415 + (t162 - t416) * t161);
t413 = t18 * qJD(1);
t412 = t182 * t295;
t411 = t183 * t296;
t379 = t293 ^ 2 * t448;
t19 = m(5) * (-t162 * t226 - t297 * t379 - t107) + m(4) * (t293 * t490 - t379) * t297 + t501 * (t125 * t182 + t126 * t183 - t107);
t410 = t19 * qJD(1);
t409 = t208 * t225;
t408 = t208 * t272;
t405 = t270 * t126;
t124 = t272 * t161;
t395 = t506 * t161;
t391 = qJD(5) * t295;
t185 = m(7) * t404;
t390 = t185 * qJD(2);
t389 = t81 * t475;
t382 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t381 = t432 / 0.2e1;
t376 = t401 / 0.2e1;
t374 = -t400 / 0.2e1;
t373 = -t400 / 0.4e1;
t372 = t400 / 0.2e1;
t252 = Ifges(7,5) * t400;
t130 = Ifges(7,6) * t270 - Ifges(7,3) * t401 - t252;
t341 = -Ifges(6,2) * t295 + t437;
t132 = -t272 * t341 + t432;
t367 = -t130 / 0.2e1 + t132 / 0.2e1;
t365 = t198 / 0.2e1 + t203 / 0.2e1;
t358 = qJ(6) * mrSges(7,2) + Ifges(6,6);
t356 = t393 * qJ(3);
t353 = -mrSges(7,2) * pkin(5) + t498;
t349 = (t411 + t412) * pkin(9);
t344 = Ifges(6,1) * t295 + t437;
t342 = Ifges(7,1) * t295 - t435;
t340 = Ifges(6,2) * t296 + t438;
t339 = Ifges(7,3) * t295 + t435;
t338 = -Ifges(7,3) * t296 + t436;
t189 = t346 * t270;
t190 = t347 * t270;
t299 = -(-mrSges(5,1) / 0.2e1 - t361) * t225 + (pkin(4) * t225 + t349) * t477 + (-t225 * t273 + t349) * t475 + t226 * mrSges(5,2) / 0.2e1 + t499 * (t411 / 0.2e1 + t412 / 0.2e1);
t312 = -t108 * t125 + t109 * t126 + t162 * t208;
t111 = -t272 * t336 + t208;
t314 = t111 * t162 + t125 * t81 + t126 * t80;
t329 = t102 * t295 - t103 * t296;
t320 = t209 + t329;
t112 = -t270 * t336 + t209;
t333 = -t295 * t74 - t296 * t73;
t323 = t112 + t333;
t2 = (t191 / 0.2e1 + t192 / 0.2e1) * t162 + t362 * t126 + t312 * t478 + t314 * t476 + t299 + (t189 / 0.2e1 + t190 / 0.2e1 + t365 * t296 + t363 * t295 + t320 * t478 + t323 * t476) * t161 + t504;
t129 = -t272 * Ifges(7,6) - t270 * t339;
t131 = -t272 * Ifges(6,6) - t270 * t341;
t133 = -t272 * Ifges(7,4) - t270 * t343;
t135 = -t272 * Ifges(6,5) - t270 * t345;
t3 = t102 * t199 + t103 * t197 + t108 * t201 + t109 * t198 - t111 * t189 - t112 * t191 - t208 * t190 - t209 * t192 + t74 * t200 + t81 * t202 + t80 * t203 + t73 * t204 - t285 * t259 + m(6) * (t102 * t108 + t103 * t109 + t208 * t209) + m(7) * (t111 * t112 + t73 * t80 + t74 * t81) + (Ifges(5,4) * t270 + t483 * t296 + (-t270 * t382 + t367) * t295) * t270 + (-t285 * mrSges(5,1) - Ifges(5,4) * t272 + (-t135 / 0.2e1 - t133 / 0.2e1 + t383 * t272) * t296 + (t131 / 0.2e1 - t129 / 0.2e1 + t382 * t272) * t295 + (-Ifges(5,2) + Ifges(5,1) - t497) * t270) * t272;
t335 = -t2 * qJD(1) + t3 * qJD(2);
t251 = Ifges(7,6) * t400;
t6 = -t270 * t251 / 0.2e1 + t103 * t202 + m(7) * (t102 * t73 + t103 * t74 + t111 * t188) + t102 * t203 - t188 * t191 + t102 * t198 - t103 * t201 + ((t381 - t208 * mrSges(6,1) + t73 * mrSges(7,2) - t111 * mrSges(7,1) + t103 * mrSges(6,3) + t252 / 0.2e1 + Ifges(6,4) * t374 + t367) * t296 + (t208 * mrSges(6,2) + t74 * mrSges(7,2) - t111 * mrSges(7,3) - t102 * mrSges(6,3) + (Ifges(6,4) / 0.2e1 - Ifges(7,5) / 0.2e1) * t401 + (Ifges(6,2) / 0.2e1 - Ifges(7,1) / 0.2e1 + Ifges(7,3) / 0.2e1 - Ifges(6,1) / 0.2e1) * t400 - t483) * t295) * t272;
t7 = (-mrSges(6,1) / 0.2e1 + t474) * t182 + t365 * t125 + (t499 * (t125 * t455 - t416 / 0.2e1) + t361 * t161) * t272 + t481 + t505;
t334 = -t7 * qJD(1) + t6 * qJD(2);
t300 = m(5) * (-t162 * t270 - t124) / 0.2e1 + m(4) * t490 / 0.2e1 + (t477 + t475) * (-t125 * t404 - t126 * t403 - t124);
t15 = t300 - t503;
t9 = t508 + t484 * t272 + (t507 * t295 - t492 * t296 + t434) * t270 + m(7) * (-t111 * t272 + t270 * t333) + m(6) * (t270 * t329 - t408) + m(5) * (-t209 * t270 - t408) + m(4) * t356;
t332 = qJD(1) * t15 + qJD(2) * t9;
t37 = -t336 * t276 + (pkin(4) * mrSges(6,2) + t440 * t296) * t296 + (pkin(4) * mrSges(6,1) - t440 * t295 + (-Ifges(6,1) - Ifges(7,1) + Ifges(6,2) + Ifges(7,3)) * t296) * t295 + (-t321 - t346) * t273;
t313 = t272 * t457;
t317 = t272 * t277;
t298 = t273 * t313 + t208 * t347 / 0.2e1 + t111 * t346 / 0.2e1 - pkin(4) * t317 / 0.2e1 + t336 * t466 + (t336 * t111 + t273 * t188) * t475 + t188 * t457 + t344 * t376 + t340 * t372 + t338 * t374 - t295 * t132 / 0.4e1 - t339 * t401 / 0.4e1 + (-t295 * t496 + t296 * t498) * t270 / 0.4e1 + (Ifges(7,1) * t401 + t130 - t252) * t295 / 0.4e1 + (t136 + t134) * t296 / 0.4e1 + (t342 + t341) * t401 / 0.4e1 + (t345 + t343) * t373 + t487 * mrSges(6,3) + ((t295 * t421 + t296 * t422) * t475 + t363 * t296 + t492 * t455) * pkin(9) + ((t471 - t73 / 0.2e1) * t295 + t422 * t451 + t487) * mrSges(7,2);
t302 = (-pkin(5) * t81 + qJ(6) * t80) * t475 - pkin(5) * t200 / 0.2e1 + qJ(6) * t461 + t108 * mrSges(6,1) / 0.2e1 + t109 * t473 + t80 * t472 + t81 * t474;
t4 = -t298 + t295 * t381 - Ifges(7,6) * t404 / 0.2e1 + t302 + t497 * t459 - t498 * t403 / 0.2e1;
t331 = t4 * qJD(2) + t37 * qJD(4);
t29 = m(7) * (t111 * t400 + t270 * t73) + t270 * t203 - t191 * t400;
t39 = (t182 / 0.4e1 + t161 * t373 - t405 / 0.4e1) * t480;
t330 = -qJD(1) * t39 + qJD(2) * t29;
t207 = (m(7) * t273 + t276) * t295;
t306 = (-t111 * t295 + (t402 + t445) * t296) * t475 - t191 * t455;
t27 = -t386 + (mrSges(7,1) / 0.2e1 + t276 * t452) * t272 + t389 - t306;
t327 = qJD(2) * t27 + qJD(4) * t207;
t36 = t433 + (t406 / 0.2e1 + t397 / 0.4e1 + t396 / 0.4e1 - t103 / 0.4e1) * t480;
t326 = qJD(2) * t36 + qJD(5) * t286;
t274 = (m(7) * pkin(9) + mrSges(7,2)) * t296;
t58 = m(7) * t415;
t40 = (t161 * t400 + t182 + t405) * t475;
t35 = (t103 + 0.2e1 * t406) * t475 + m(7) * t471 + t203;
t28 = t276 * t372 + t431 / 0.2e1 + t389 + t306;
t20 = t199 * t451 + t200 * t452 + t361 * t272 + t454 * t493 + t301 - t308;
t14 = t300 + t503;
t13 = t304 - t303;
t11 = t489 * t468 + (t305 + t315) * t161;
t8 = t317 * t468 + t161 * t313 - t492 * t125 / 0.2e1 - t500 * t182 / 0.2e1 + t481 + t499 * (t125 * t376 + t126 * t372) - t505;
t5 = t298 + (-t295 * t382 - t296 * t383) * t270 + (-Ifges(7,2) / 0.2e1 - Ifges(6,3) / 0.2e1) * t272 + t302;
t1 = t299 + t126 * t461 + t197 * t469 + t314 * t475 + t312 * t477 + (-t189 - t190) * t468 + (t466 - t192 / 0.2e1) * t162 + (t320 * t477 + t323 * t475 + t482) * t161 - t504;
t10 = [qJD(2) * t19 + qJD(4) * t18, t14 * qJD(3) + t1 * qJD(4) + t8 * qJD(5) + t40 * qJD(6) + t410 + (-m(6) * t409 + m(5) * (-t209 * t226 - t409) + t226 * t434 + m(4) * (-pkin(2) * t448 + t297 * t356) * t293 + (m(5) * t285 - mrSges(4,1) * t294 + mrSges(5,1) * t270 + mrSges(4,2) * t292 - mrSges(5,2) * t272 - mrSges(3,1)) * t378 - (m(7) * t111 - t484) * t225 + (m(6) * t103 + m(7) * t73 + t492) * t183 + (-m(6) * t102 + m(7) * t74 - t507) * t182 + (-mrSges(3,2) + t508) * t399) * qJD(2), qJD(2) * t14, t413 + t1 * qJD(2) + t11 * qJD(5) - t58 * qJD(6) + ((-pkin(4) * t162 - t395) * t477 + (t162 * t273 - t395) * t475) * t479 + ((-mrSges(5,1) + t491) * t162 + (-t499 * t509 + mrSges(5,2)) * t161) * qJD(4), t8 * qJD(2) + t11 * qJD(4) + (-t500 * t126 + (mrSges(6,2) - mrSges(7,3)) * t125) * qJD(5) + ((-pkin(5) * t126 - qJ(6) * t125) * qJD(5) / 0.2e1 + qJD(6) * t469) * t480, m(7) * t126 * qJD(5) + t40 * qJD(2) - t58 * qJD(4); qJD(3) * t15 - qJD(4) * t2 - qJD(5) * t7 - qJD(6) * t39 - t410, qJD(3) * t9 + qJD(4) * t3 + qJD(5) * t6 + qJD(6) * t29, qJD(4) * t20 + qJD(5) * t13 + t332, t20 * qJD(3) + t5 * qJD(5) + t28 * qJD(6) + ((-pkin(4) * t209 + pkin(9) * t486) * t477 + (pkin(9) * t485 + t112 * t273) * t475) * t479 + t335 + (t131 * t451 - t273 * t189 + t112 * t276 + Ifges(5,6) * t272 + t208 * mrSges(5,2) + pkin(4) * t190 + (t493 * t296 + (-t199 + t200) * t295) * pkin(9) + (t338 * t455 - Ifges(5,5)) * t270 + (t295 * t498 + t296 * t496) * t459 + (t270 * t340 + t133 + t135) * t454 + (t129 + (t342 + t344) * t270) * t452 + (t277 - mrSges(5,1)) * t209 + t486 * mrSges(6,3) + t485 * mrSges(7,2)) * qJD(4), t13 * qJD(3) + t5 * qJD(4) + t35 * qJD(6) + t334 + (-t251 + (-m(7) * pkin(5) - t500) * t103 + (-mrSges(6,2) + t286) * t102 + (t295 * t353 + t296 * t358) * t272) * qJD(5), qJD(4) * t28 + qJD(5) * t35 + t330; -qJD(2) * t15, -qJD(4) * t17 - qJD(5) * t12 + qJD(6) * t185 - t332, 0, -t414, -t417 + (-t489 - t321) * qJD(5) + m(7) * t295 * qJD(6), m(7) * t391 + t390; qJD(2) * t2 - t413, qJD(3) * t17 - qJD(5) * t4 - qJD(6) * t27 - t335, t414, -qJD(5) * t37 - qJD(6) * t207, t274 * qJD(6) + (Ifges(7,6) - t358) * t391 - t331 + (t353 * t296 + (m(7) * t337 + t491) * pkin(9)) * qJD(5), qJD(5) * t274 - t327; qJD(2) * t7, qJD(3) * t12 + qJD(4) * t4 + qJD(6) * t36 - t334, t417, t331, t286 * qJD(6), t326; qJD(2) * t39, -qJD(3) * t185 + qJD(4) * t27 - qJD(5) * t36 - t330, -t390, t327, -t326, 0;];
Cq  = t10;
