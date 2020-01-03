% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR16_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_coriolismatJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR16_coriolismatJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_coriolismatJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR16_coriolismatJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR16_coriolismatJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:48
% EndTime: 2019-12-31 20:45:03
% DurationCPUTime: 6.98s
% Computational Cost: add. (12743->652), mult. (30706->915), div. (0->0), fcn. (30217->8), ass. (0->344)
t571 = -Ifges(6,3) / 0.2e1;
t359 = sin(pkin(5));
t361 = sin(qJ(5));
t364 = cos(qJ(5));
t366 = cos(qJ(2));
t362 = sin(qJ(4));
t363 = sin(qJ(2));
t464 = t362 * t363;
t245 = (-t361 * t464 + t364 * t366) * t359;
t543 = t245 / 0.2e1;
t462 = t363 * t364;
t246 = (t361 * t366 + t362 * t462) * t359;
t542 = t246 / 0.2e1;
t485 = t364 * mrSges(6,1);
t487 = t361 * mrSges(6,2);
t416 = t485 - t487;
t570 = t362 * t416;
t555 = m(6) * pkin(4);
t569 = -t555 - mrSges(5,1) - t416;
t365 = cos(qJ(4));
t461 = t363 * t365;
t446 = t359 * t461;
t568 = t446 * t571;
t360 = cos(pkin(5));
t471 = t359 * t366;
t291 = t360 * t365 - t362 * t471;
t491 = t291 * mrSges(5,3);
t200 = -t291 * t361 + t359 * t462;
t472 = t359 * t363;
t201 = t291 * t364 + t361 * t472;
t104 = -mrSges(6,1) * t200 + mrSges(6,2) * t201;
t215 = mrSges(5,1) * t472 - t491;
t566 = -t104 + t215;
t355 = Ifges(6,5) * t364;
t502 = Ifges(6,6) * t361;
t565 = t355 - t502;
t356 = Ifges(6,4) * t364;
t564 = -Ifges(6,2) * t361 + t356;
t331 = Ifges(6,1) * t361 + t356;
t182 = mrSges(6,2) * t446 + mrSges(6,3) * t245;
t183 = -mrSges(6,1) * t446 - mrSges(6,3) * t246;
t524 = t364 / 0.2e1;
t530 = -t361 / 0.2e1;
t563 = t182 * t524 + t183 * t530;
t562 = t502 / 0.2e1 - t355 / 0.2e1;
t337 = pkin(4) * t365 + pkin(9) * t362;
t367 = -pkin(2) - pkin(8);
t454 = t365 * t367;
t243 = t364 * t337 - t361 * t454;
t244 = t361 * t337 + t364 * t454;
t405 = -t243 * t361 + t244 * t364;
t561 = Ifges(6,5) * t542 + Ifges(6,6) * t543;
t519 = pkin(1) * t360;
t349 = t366 * t519;
t549 = pkin(3) + pkin(7);
t235 = -t549 * t472 + t349;
t560 = m(4) / 0.2e1;
t559 = m(5) / 0.2e1;
t558 = -m(6) / 0.2e1;
t557 = m(6) / 0.2e1;
t556 = -pkin(4) / 0.2e1;
t554 = -mrSges(6,1) / 0.2e1;
t553 = mrSges(6,2) / 0.2e1;
t552 = -mrSges(6,3) / 0.2e1;
t551 = mrSges(6,3) / 0.2e1;
t192 = t367 * t360 - t235;
t425 = -qJ(3) * t363 - pkin(1);
t212 = (t367 * t366 + t425) * t359;
t94 = t192 * t365 - t212 * t362;
t76 = -pkin(4) * t472 - t94;
t550 = t76 / 0.2e1;
t290 = t360 * t362 + t365 * t471;
t497 = t201 * mrSges(6,3);
t141 = mrSges(6,1) * t290 - t497;
t548 = -t141 / 0.2e1;
t547 = t200 / 0.2e1;
t546 = t201 / 0.2e1;
t518 = pkin(4) * t362;
t321 = -pkin(9) * t365 + qJ(3) + t518;
t463 = t362 * t367;
t230 = t364 * t321 - t361 * t463;
t545 = t230 / 0.2e1;
t544 = t244 / 0.2e1;
t541 = t290 / 0.4e1;
t540 = t291 / 0.2e1;
t484 = t364 * mrSges(6,2);
t488 = t361 * mrSges(6,1);
t324 = t484 + t488;
t298 = t365 * t324;
t539 = t298 / 0.2e1;
t465 = t361 * t365;
t311 = -mrSges(6,2) * t362 - mrSges(6,3) * t465;
t538 = -t311 / 0.2e1;
t537 = t311 / 0.2e1;
t456 = t364 * t365;
t313 = mrSges(6,1) * t362 - mrSges(6,3) * t456;
t536 = -t313 / 0.2e1;
t535 = t416 / 0.2e1;
t534 = t324 / 0.2e1;
t325 = Ifges(6,5) * t361 + Ifges(6,6) * t364;
t533 = t325 / 0.2e1;
t509 = Ifges(6,4) * t361;
t327 = Ifges(6,2) * t364 + t509;
t532 = -t327 / 0.4e1;
t531 = -t360 / 0.2e1;
t529 = t361 / 0.2e1;
t528 = -t362 / 0.2e1;
t527 = t362 / 0.2e1;
t526 = t362 / 0.4e1;
t525 = -t364 / 0.2e1;
t523 = t364 / 0.4e1;
t522 = -t365 / 0.2e1;
t521 = t365 / 0.2e1;
t520 = -t367 / 0.2e1;
t295 = pkin(7) * t471 + t363 * t519;
t353 = t360 * qJ(3);
t248 = -t295 - t353;
t348 = pkin(3) * t471;
t211 = -t248 + t348;
t102 = pkin(4) * t290 - pkin(9) * t291 + t211;
t95 = t192 * t362 + t212 * t365;
t77 = pkin(9) * t472 + t95;
t43 = t102 * t364 - t361 * t77;
t517 = t43 * mrSges(6,3);
t44 = t102 * t361 + t364 * t77;
t516 = t44 * mrSges(6,3);
t515 = m(6) * qJD(3);
t513 = Ifges(5,4) * t291;
t512 = Ifges(5,4) * t362;
t511 = Ifges(5,4) * t365;
t510 = Ifges(6,4) * t201;
t507 = Ifges(6,5) * t290;
t506 = Ifges(6,5) * t362;
t503 = Ifges(6,6) * t290;
t501 = Ifges(6,6) * t362;
t500 = Ifges(6,3) * t291;
t499 = Ifges(6,3) * t365;
t498 = t200 * mrSges(6,3);
t496 = t230 * mrSges(6,3);
t231 = t361 * t321 + t364 * t463;
t495 = t231 * mrSges(6,3);
t494 = t290 * mrSges(5,1);
t493 = t290 * mrSges(5,3);
t492 = t291 * mrSges(5,2);
t294 = pkin(7) * t472 - t349;
t490 = t294 * mrSges(3,2);
t347 = pkin(2) * t472;
t234 = t347 + (pkin(8) * t363 - qJ(3) * t366) * t359;
t236 = t295 + t348;
t132 = -t234 * t362 + t236 * t365;
t112 = -pkin(4) * t471 - t132;
t118 = Ifges(6,4) * t246 + Ifges(6,2) * t245 - Ifges(6,6) * t446;
t119 = Ifges(6,1) * t246 + Ifges(6,4) * t245 - Ifges(6,5) * t446;
t133 = t365 * t234 + t362 * t236;
t140 = -mrSges(6,2) * t290 + t498;
t145 = -mrSges(6,1) * t245 + mrSges(6,2) * t246;
t283 = Ifges(5,4) * t290;
t449 = Ifges(5,5) * t472;
t152 = Ifges(5,1) * t291 - t283 + t449;
t176 = t492 + t494;
t415 = Ifges(5,1) * t362 + t511;
t205 = (t366 * Ifges(5,5) + t415 * t363) * t359;
t214 = -mrSges(5,2) * t472 - t493;
t323 = mrSges(5,1) * t365 - mrSges(5,2) * t362;
t247 = t323 * t472;
t249 = (-pkin(2) * t366 + t425) * t359;
t253 = -pkin(2) * t360 + t294;
t276 = (mrSges(5,1) * t366 - mrSges(5,3) * t464) * t359;
t277 = (-mrSges(5,2) * t366 + mrSges(5,3) * t461) * t359;
t292 = (mrSges(4,2) * t366 - mrSges(4,3) * t363) * t359;
t293 = -qJ(3) * t471 + t347;
t354 = t360 * mrSges(4,3);
t448 = mrSges(4,1) * t471;
t309 = -t354 - t448;
t344 = Ifges(4,5) * t472;
t345 = Ifges(3,5) * t471;
t414 = Ifges(5,2) * t365 + t512;
t433 = t568 - (t366 * Ifges(5,6) + t414 * t363) * t359 / 0.2e1 + t561;
t151 = -Ifges(5,2) * t290 + Ifges(5,6) * t472 + t513;
t71 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t290;
t440 = t151 / 0.2e1 - t71 / 0.2e1;
t113 = pkin(9) * t471 + t133;
t153 = t349 + (-t337 - t549) * t472;
t57 = -t113 * t361 + t153 * t364;
t58 = t113 * t364 + t153 * t361;
t72 = Ifges(6,2) * t200 + t503 + t510;
t195 = Ifges(6,4) * t200;
t73 = Ifges(6,1) * t201 + t195 + t507;
t3 = ((t253 * mrSges(4,1) - t249 * mrSges(4,3) + Ifges(5,5) * t540 - Ifges(5,6) * t290 / 0.2e1 + Ifges(3,4) * t471 + (-Ifges(4,4) + Ifges(3,5) / 0.2e1) * t360 + (-mrSges(3,2) * pkin(1) + Ifges(4,6) * t366) * t359) * t366 + (t152 * t527 - t249 * mrSges(4,2) + t440 * t365 + (-Ifges(3,6) + Ifges(4,5) / 0.2e1) * t360 + (-pkin(1) * mrSges(3,1) + (Ifges(5,5) * t527 + Ifges(5,6) * t521 - Ifges(3,4) - Ifges(4,6)) * t363) * t359 + (t248 + t295) * mrSges(4,1) + (Ifges(3,1) - Ifges(3,2) + Ifges(4,2) - Ifges(4,3) + Ifges(5,3)) * t471) * t363) * t359 + (t490 + t344 / 0.2e1 + t345 / 0.2e1 + (mrSges(4,2) - mrSges(3,1)) * t295) * t360 + t433 * t290 + t294 * t309 + t293 * t292 + t94 * t276 + t95 * t277 - t211 * t247 + t235 * t176 + t133 * t214 + t132 * t215 + t44 * t182 + t43 * t183 + t76 * t145 + t58 * t140 + t57 * t141 + t112 * t104 + t205 * t540 + t73 * t542 + t72 * t543 + t119 * t546 + t118 * t547 + m(6) * (t112 * t76 + t43 * t57 + t44 * t58) + m(5) * (t132 * t94 + t133 * t95 + t211 * t235) + m(4) * (t248 * t294 + t249 * t293 + t253 * t295);
t489 = t3 * qJD(1);
t486 = t361 * t72;
t483 = t364 * mrSges(6,3);
t482 = t364 * t73;
t114 = -t290 * t565 + t500;
t115 = Ifges(6,6) * t291 - t290 * t564;
t332 = Ifges(6,1) * t364 - t509;
t116 = Ifges(6,5) * t291 - t332 * t290;
t154 = t324 * t290;
t168 = mrSges(6,3) * t290 * t361 - mrSges(6,2) * t291;
t169 = mrSges(6,1) * t291 + t290 * t483;
t175 = mrSges(5,1) * t291 - mrSges(5,2) * t290;
t177 = -Ifges(5,2) * t291 - t283;
t178 = -Ifges(5,1) * t290 - t513;
t437 = t472 / 0.2e1;
t453 = -Ifges(5,5) * t290 - Ifges(5,6) * t291;
t179 = pkin(4) * t291 + pkin(9) * t290;
t55 = t179 * t364 - t361 * t94;
t56 = t179 * t361 + t364 * t94;
t4 = t453 * t437 + t56 * t140 + t44 * t168 + t55 * t141 + t43 * t169 - t76 * t154 + m(6) * (t43 * t55 + t44 * t56) + t116 * t546 + t115 * t547 + t94 * t214 + t211 * t175 + (t178 / 0.2e1 - t440) * t291 + (-t482 / 0.2e1 + t486 / 0.2e1 + t94 * mrSges(5,3) + t114 / 0.2e1 - t152 / 0.2e1 - t177 / 0.2e1) * t290 + (m(6) * t76 - t491 - t566) * t95;
t481 = t4 * qJD(1);
t480 = t57 * t361;
t479 = t58 * t364;
t103 = mrSges(6,1) * t201 + mrSges(6,2) * t200;
t105 = Ifges(6,5) * t200 - Ifges(6,6) * t201;
t106 = -Ifges(6,2) * t201 + t195;
t107 = Ifges(6,1) * t200 - t510;
t7 = t290 * t105 / 0.2e1 + t43 * t140 - t44 * t141 + t76 * t103 + (-t516 + t107 / 0.2e1 - t72 / 0.2e1) * t201 + (-t517 + t73 / 0.2e1 + t106 / 0.2e1) * t200;
t478 = t7 * qJD(1);
t274 = t360 * t364 + t361 * t446;
t275 = t360 * t361 - t364 * t446;
t447 = t359 * t464;
t455 = t365 * t214;
t14 = t275 * t140 + t274 * t141 + (t176 - t309) * t360 + (t362 * t566 - t292 - t455) * t472 + m(6) * (t274 * t43 + t275 * t44 - t76 * t447) + m(5) * (t211 * t360 + (t362 * t94 - t365 * t95) * t472) + m(4) * (-t248 * t360 - t249 * t472);
t477 = t14 * qJD(1);
t474 = t274 * t361;
t473 = t275 * t364;
t470 = t361 * t169;
t268 = t365 * t564 + t501;
t469 = t361 * t268;
t468 = t361 * t311;
t312 = mrSges(6,1) * t365 + t362 * t483;
t467 = t361 * t312;
t466 = t361 * t362;
t460 = t364 * t168;
t270 = t332 * t365 + t506;
t459 = t364 * t270;
t310 = -mrSges(6,2) * t365 + mrSges(6,3) * t466;
t458 = t364 * t310;
t457 = t364 * t313;
t357 = t361 ^ 2;
t358 = t364 ^ 2;
t452 = t357 + t358;
t451 = m(6) * t550;
t450 = t95 * t558;
t443 = mrSges(5,3) * t528;
t442 = mrSges(5,3) * t521;
t441 = t107 / 0.4e1 - t72 / 0.4e1;
t439 = t73 / 0.4e1 + t106 / 0.4e1;
t438 = -t472 / 0.4e1;
t434 = t104 / 0.2e1 - t215 / 0.2e1;
t432 = -t145 / 0.2e1 + t276 / 0.2e1;
t431 = t154 / 0.2e1 + t214 / 0.2e1;
t301 = t331 * t365;
t430 = -t268 / 0.4e1 - t301 / 0.4e1;
t300 = t327 * t365;
t429 = t270 / 0.4e1 - t300 / 0.4e1;
t266 = Ifges(6,3) * t362 + t365 * t565;
t330 = -Ifges(5,2) * t362 + t511;
t428 = t330 / 0.2e1 - t266 / 0.2e1;
t427 = t331 / 0.4e1 + t564 / 0.4e1;
t426 = t332 / 0.4e1 + t532;
t424 = t452 * t365;
t422 = -t447 / 0.2e1;
t421 = t362 * t437;
t420 = t365 * t437;
t419 = mrSges(6,3) * (-t358 / 0.2e1 - t357 / 0.2e1);
t418 = -t330 / 0.4e1 - t415 / 0.4e1 + t266 / 0.4e1;
t417 = t362 * mrSges(5,1) + t365 * mrSges(5,2);
t411 = -t361 * t55 + t364 * t56;
t410 = t479 - t480;
t296 = t416 * t365;
t299 = t365 * t325;
t24 = t230 * t311 - t231 * t313 - t299 * t527 + (-t367 * t296 + (-t495 - t268 / 0.2e1 - t301 / 0.2e1) * t364 + (t496 + t300 / 0.2e1 - t270 / 0.2e1) * t361) * t365;
t370 = (-t496 / 0.2e1 + t429) * t200 + (-t495 / 0.2e1 + t430) * t201 + t140 * t545 + t231 * t548 - t299 * t541 + t105 * t526 + t43 * t537 + t44 * t536 + t296 * t550;
t376 = t103 * t520 + (-t516 / 0.2e1 + t441) * t364 + (t517 / 0.2e1 - t439) * t361;
t379 = t58 * t553 + t57 * t554 - t561;
t5 = (Ifges(6,3) * t437 + t376) * t365 + t370 + t379;
t409 = t5 * qJD(1) + t24 * qJD(2);
t388 = t295 + 0.2e1 * t353;
t392 = t140 * t530 + t141 * t525;
t368 = -t354 - m(4) * t388 / 0.2e1 - m(5) * (t348 + t388) / 0.2e1 + (t230 * t274 + t231 * t275 + t361 * t44 + t364 * t43 + t446 * t463) * t558 + t274 * t536 + t275 * t538 - t494 / 0.2e1 - t492 / 0.2e1 + t392;
t406 = t132 * t365 + t133 * t362;
t372 = (-t112 * t365 + t410 * t362) * t557 + t406 * t559 + t295 * t560;
t386 = t277 / 0.2e1 + t563;
t13 = (mrSges(5,2) * t531 + t432) * t365 + (mrSges(5,1) * t531 + t298 * t437 + t386) * t362 + t368 + t372;
t59 = m(6) * (t230 * t364 + t231 * t361) + t468 + t457 + mrSges(4,3) + (m(5) + m(4)) * qJ(3) + t417;
t408 = -qJD(1) * t13 + qJD(2) * t59;
t378 = (t200 * t529 + t201 * t525) * mrSges(6,3) + t392;
t373 = t103 * t522 + t378 * t362;
t399 = t274 * mrSges(6,1) / 0.2e1 - t275 * mrSges(6,2) / 0.2e1;
t18 = t373 - t399;
t377 = (-t468 / 0.2e1 - t457 / 0.2e1 + t365 * t419) * t362 + t296 * t522;
t398 = t487 / 0.2e1 - t485 / 0.2e1;
t38 = t377 + t398;
t407 = t18 * qJD(1) + t38 * qJD(2);
t404 = t473 - t474;
t403 = -t361 * t43 + t364 * t44 - t95;
t402 = t411 + t76;
t401 = t296 * t556 + t526 * t565;
t400 = mrSges(6,2) * t544 + t243 * t554;
t397 = -t484 / 0.2e1 - t488 / 0.2e1;
t396 = pkin(9) * t538 + t430;
t395 = pkin(9) * t536 + t429;
t394 = -t491 / 0.2e1 + t434;
t393 = t493 / 0.2e1 + t431;
t391 = t327 * t529 + t331 * t525;
t390 = t404 * pkin(9);
t267 = Ifges(6,6) * t365 - t362 * t564;
t269 = Ifges(6,5) * t365 - t332 * t362;
t297 = t324 * t362;
t369 = (t230 * t55 + t231 * t56 + t243 * t43 + t244 * t44) * t557 + qJ(3) * t175 / 0.2e1 + t200 * t267 / 0.4e1 + t201 * t269 / 0.4e1 + t211 * t323 / 0.2e1 + t169 * t545 + t231 * t168 / 0.2e1 + t243 * t141 / 0.2e1 + t140 * t544 + t43 * t312 / 0.2e1 + t44 * t310 / 0.2e1 + t55 * t313 / 0.2e1 + t56 * t537 - t297 * t550 + t95 * t539;
t371 = (-pkin(4) * t112 + t410 * pkin(9)) * t558 + pkin(4) * t145 / 0.2e1 + t112 * t535 - t132 * mrSges(5,1) / 0.2e1 + t133 * mrSges(5,2) / 0.2e1 + t245 * t532 - t246 * t331 / 0.4e1 - Ifges(5,3) * t471 / 0.2e1;
t265 = -t362 * t565 + t499;
t334 = Ifges(5,1) * t365 - t512;
t380 = t414 / 0.4e1 - t334 / 0.4e1 + t265 / 0.4e1 + t469 / 0.4e1 - t459 / 0.4e1;
t381 = t178 / 0.4e1 - t151 / 0.4e1 + t71 / 0.4e1 + t116 * t523 - t361 * t115 / 0.4e1;
t382 = t114 / 0.4e1 - t152 / 0.4e1 - t177 / 0.4e1 + t486 / 0.4e1 - t482 / 0.4e1;
t2 = (-t119 / 0.4e1 + pkin(9) * t183 / 0.2e1 + t57 * t551) * t361 + t418 * t291 + (-t118 / 0.4e1 - pkin(9) * t182 / 0.2e1 + t58 * t552) * t364 + t380 * t290 + ((t325 / 0.4e1 - 0.3e1 / 0.4e1 * Ifges(5,6)) * t472 + (t450 + t393) * t367 + t381) * t365 + (-0.3e1 / 0.4e1 * t449 + (t451 + t394) * t367 + t382) * t362 + t369 + t371;
t20 = m(6) * (t230 * t243 + t231 * t244) + t244 * t311 + t231 * t310 + t243 * t313 + t230 * t312 + qJ(3) * t323 + (-t415 / 0.2e1 + t267 * t530 + t269 * t524 + t367 * t297 - t428) * t365 + (t265 / 0.2e1 - t334 / 0.2e1 + t414 / 0.2e1 + t469 / 0.2e1 - t459 / 0.2e1 + (-m(6) * t454 + t298) * t367) * t362;
t30 = ((-t230 * t361 + t231 * t364) * t557 + t311 * t524 + t313 * t530 + t297 / 0.2e1) * t365 + ((t405 - 0.2e1 * t454) * t557 + t458 / 0.2e1 - t467 / 0.2e1 + t539) * t362;
t389 = t2 * qJD(1) + t20 * qJD(2) + t30 * qJD(3);
t387 = (t473 / 0.2e1 - t474 / 0.2e1) * mrSges(6,3);
t10 = t390 * t557 + t387 + (t402 * t558 - t460 / 0.2e1 + t470 / 0.2e1 + (t535 + t555 / 0.2e1 + mrSges(5,1) / 0.2e1) * t472 - t394) * t362 + (mrSges(5,2) * t437 + t140 * t525 + t141 * t529 + t403 * t558 - t393) * t365;
t227 = (-0.1e1 + t452) * t365 * t362;
t385 = -t10 * qJD(1) + t30 * qJD(2) + t227 * t515;
t384 = -t500 / 0.2e1 + t55 * t554 + t56 * t553;
t198 = (t534 + t397) * t365;
t375 = pkin(9) * t419 + t324 * t520 - t427 * t361 + t426 * t364;
t22 = (t506 / 0.2e1 + t395) * t364 + (-t501 / 0.2e1 + t396) * t361 + (t571 + t375) * t365 + t400 + t401;
t374 = t103 * t556 + t427 * t200 + t426 * t201 + t76 * t534 + t541 * t565;
t9 = (t507 / 0.2e1 + (-t497 / 0.2e1 + t548) * pkin(9) + t439) * t364 + (-t503 / 0.2e1 + (-t140 / 0.2e1 + t498 / 0.2e1) * pkin(9) + t441) * t361 + t374 + t384;
t96 = pkin(4) * t324 + (-t331 / 0.2e1 - t564 / 0.2e1) * t364 + (-t332 / 0.2e1 + t327 / 0.2e1) * t361;
t383 = t9 * qJD(1) + t22 * qJD(2) - t198 * qJD(3) - t96 * qJD(4);
t199 = t324 * t522 + t397 * t365;
t37 = t377 - t398;
t29 = t30 * qJD(4);
t21 = t499 / 0.2e1 + t395 * t364 + t396 * t361 + t375 * t365 - t400 + t401 + t562 * t362;
t19 = t373 + t399;
t12 = t386 * t362 + t432 * t365 + t298 * t422 + t360 * t417 / 0.2e1 - t368 + t448 + t372;
t11 = t140 * t456 / 0.2e1 + t465 * t548 - t169 * t466 / 0.2e1 - t154 * t522 + t455 / 0.2e1 + t291 * t443 + t215 * t528 + t290 * t442 - t416 * t422 + mrSges(5,2) * t420 + mrSges(5,1) * t421 + t387 + (pkin(4) * t447 + t402 * t362 + t403 * t365 + t390) * t557 + (t460 + t104) * t527;
t8 = t378 * pkin(9) + t290 * t562 + t441 * t361 + t439 * t364 + t374 - t384;
t6 = t376 * t365 + t370 - t379 + t568;
t1 = -t371 + (t367 * t443 + t418) * t291 + (t367 * t442 + t380) * t290 + Ifges(5,6) * t420 + Ifges(5,5) * t421 + t118 * t523 + t361 * t119 / 0.4e1 + (Ifges(5,5) * t438 + (t451 + t434) * t367 + t382) * t362 + t369 + t480 * t552 + t479 * t551 + ((t450 + t431) * t367 + t381 + (t325 + Ifges(5,6)) * t438) * t365 + t563 * pkin(9);
t15 = [qJD(2) * t3 + qJD(3) * t14 + qJD(4) * t4 + qJD(5) * t7, t12 * qJD(3) + t1 * qJD(4) + t6 * qJD(5) + t489 + (-t295 * mrSges(3,1) + t295 * mrSges(4,2) - t294 * mrSges(4,3) - qJ(3) * t247 + t112 * t298 + t231 * t182 + t230 * t183 + t268 * t543 + t270 * t542 + t58 * t311 + t57 * t313 + t344 + t345 + t490 + (t235 * mrSges(5,1) - t133 * mrSges(5,3) + t367 * t277 + t433) * t362 + (t235 * mrSges(5,2) + t205 / 0.2e1 - t132 * mrSges(5,3) + t118 * t530 + t119 * t524 + (-t145 + t276) * t367) * t365 + 0.2e1 * (-pkin(2) * t295 - qJ(3) * t294) * t560 + 0.2e1 * (-t112 * t454 + t230 * t57 + t231 * t58) * t557 + 0.2e1 * (qJ(3) * t235 + t406 * t367) * t559 + ((-pkin(2) * mrSges(4,1) + Ifges(5,5) * t521 + Ifges(5,6) * t528 - Ifges(4,4)) * t366 + (-qJ(3) * mrSges(4,1) + t334 * t527 + t428 * t365 - Ifges(3,6)) * t363) * t359) * qJD(2), t477 + t12 * qJD(2) + t11 * qJD(4) + t19 * qJD(5) + (t404 + t446) * t362 * t515, t481 + t1 * qJD(2) + t11 * qJD(3) + (-t94 * mrSges(5,2) + t411 * mrSges(6,3) + pkin(4) * t154 + t115 * t524 + t116 * t529 + t391 * t290 + t291 * t533 + t453 + t569 * t95 + (m(6) * t411 + t460 - t470) * pkin(9)) * qJD(4) + t8 * qJD(5), t478 + t6 * qJD(2) + t19 * qJD(3) + t8 * qJD(4) + (-mrSges(6,1) * t44 - mrSges(6,2) * t43 + t105) * qJD(5); -qJD(3) * t13 + qJD(4) * t2 + qJD(5) * t5 - t489, qJD(3) * t59 + qJD(4) * t20 + qJD(5) * t24, qJD(5) * t37 + t29 + t408, t21 * qJD(5) + t389 + (pkin(4) * t297 + t267 * t524 + t269 * t529 + (m(6) * t405 + t458 - t467) * pkin(9) + (t569 * t367 - Ifges(5,5) + t391) * t362 + (-t367 * mrSges(5,2) - Ifges(5,6) + t533) * t365 + t405 * mrSges(6,3)) * qJD(4), t37 * qJD(3) + t21 * qJD(4) + (-mrSges(6,1) * t231 - mrSges(6,2) * t230 - t299) * qJD(5) + t409; qJD(2) * t13 - qJD(4) * t10 + qJD(5) * t18 - t477, qJD(5) * t38 + t29 - t408, m(6) * t227 * qJD(4), (-t570 + m(6) * (pkin(9) * t424 - t518) + mrSges(6,3) * t424 - t417) * qJD(4) + t199 * qJD(5) + t385, t199 * qJD(4) - qJD(5) * t570 + t407; -qJD(2) * t2 + qJD(3) * t10 + qJD(5) * t9 - t481, qJD(5) * t22 - t389, -t198 * qJD(5) - t385, -t96 * qJD(5), (-pkin(9) * t416 + t565) * qJD(5) + t383; -qJD(2) * t5 - qJD(3) * t18 - qJD(4) * t9 - t478, -qJD(3) * t38 - qJD(4) * t22 - t409, qJD(4) * t198 - t407, -t383, 0;];
Cq = t15;
