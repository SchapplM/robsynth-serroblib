% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR12_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_coriolismatJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:16
% EndTime: 2019-03-09 04:18:31
% DurationCPUTime: 7.94s
% Computational Cost: add. (13787->560), mult. (25371->725), div. (0->0), fcn. (24160->6), ass. (0->321)
t625 = mrSges(6,1) / 0.2e1;
t350 = sin(qJ(3));
t354 = -pkin(1) - pkin(7);
t334 = t350 * t354;
t299 = -pkin(4) * t350 + t334;
t351 = cos(qJ(5));
t282 = t351 * t299;
t352 = cos(qJ(3));
t447 = pkin(3) * t352 + qJ(4) * t350;
t284 = pkin(8) * t352 + t447;
t349 = sin(qJ(5));
t186 = -t284 * t349 + t282;
t187 = t351 * t284 + t349 * t299;
t348 = sin(qJ(6));
t456 = t351 * t352;
t412 = t456 / 0.2e1;
t465 = t349 * t352;
t415 = t465 / 0.2e1;
t539 = cos(qJ(6));
t438 = pkin(5) * t539;
t580 = m(7) * pkin(5);
t442 = t580 / 0.2e1;
t535 = pkin(5) * t348;
t420 = t539 * t349;
t386 = t348 * t351 + t420;
t251 = t386 * t352;
t211 = -mrSges(7,1) * t350 - mrSges(7,3) * t251;
t565 = t211 / 0.2e1;
t419 = t539 * t351;
t249 = t348 * t465 - t352 * t419;
t505 = t249 * mrSges(7,3);
t209 = mrSges(7,2) * t350 - t505;
t567 = t209 / 0.2e1;
t579 = -mrSges(6,2) / 0.2e1;
t559 = t251 / 0.2e1;
t561 = -t249 / 0.2e1;
t389 = Ifges(7,5) * t559 + Ifges(7,6) * t561;
t544 = t350 / 0.2e1;
t149 = -pkin(5) * t350 + t282 + (-pkin(9) * t352 - t284) * t349;
t166 = pkin(9) * t456 + t187;
t84 = t149 * t539 - t348 * t166;
t85 = t348 * t149 + t166 * t539;
t586 = Ifges(7,3) * t544 + t85 * mrSges(7,2) / 0.2e1 - t84 * mrSges(7,1) / 0.2e1 - t389;
t624 = -Ifges(6,3) * t350 / 0.2e1 + t186 * t625 + t187 * t579 + (t348 * t85 + t539 * t84) * t442 + Ifges(6,5) * t415 + Ifges(6,6) * t412 + t535 * t567 + t438 * t565 - t586;
t385 = -t348 * t349 + t419;
t599 = t385 * t535 - t386 * t438;
t390 = m(7) * t599;
t500 = t351 * mrSges(6,2);
t501 = t349 * mrSges(6,1);
t274 = t386 * mrSges(7,1);
t503 = t385 * mrSges(7,2);
t594 = t503 / 0.2e1 + t274 / 0.2e1;
t623 = -t390 / 0.2e1 + t501 / 0.2e1 + t500 / 0.2e1 + t594;
t461 = t350 * t351;
t250 = -t348 * t461 - t350 * t420;
t504 = t250 * mrSges(7,2);
t248 = t385 * t350;
t508 = t248 * mrSges(7,1);
t450 = t508 / 0.2e1 + t504 / 0.2e1;
t622 = t461 * t625 + t450;
t391 = m(7) * (t249 * t535 + t251 * t438);
t224 = t251 * mrSges(7,1);
t506 = t249 * mrSges(7,2);
t592 = t506 / 0.2e1 - t224 / 0.2e1;
t621 = t592 - t391 / 0.2e1;
t473 = t385 * t251;
t620 = m(7) * (-t249 * t386 - t473);
t353 = -pkin(3) - pkin(8);
t464 = t349 * t353;
t298 = -t349 * pkin(9) + t464;
t300 = (-pkin(9) + t353) * t351;
t204 = t298 * t539 + t348 * t300;
t403 = -t298 * t348 + t539 * t300;
t448 = -Ifges(7,5) * t386 - Ifges(7,6) * t385;
t37 = -t204 * mrSges(7,1) - t403 * mrSges(7,2) + t448;
t619 = t37 * qJD(6);
t302 = mrSges(6,1) * t351 - mrSges(6,2) * t349;
t462 = t350 * t302;
t618 = -t462 / 0.2e1;
t337 = t352 * qJ(4);
t446 = -t350 * pkin(3) + t337;
t617 = m(5) * t446;
t534 = pkin(5) * t349;
t325 = qJ(4) + t534;
t536 = m(7) * t325;
t344 = t352 * pkin(4);
t301 = t352 * t354 - t344;
t283 = t351 * t301;
t297 = qJ(2) - t446;
t268 = pkin(8) * t350 + t297;
t406 = pkin(9) * t350 + t268;
t155 = -t349 * t406 - t283;
t148 = pkin(5) * t352 + t155;
t470 = t301 * t349;
t156 = t351 * t406 - t470;
t422 = t539 * t156;
t79 = t348 * t148 + t422;
t86 = -t348 * t155 - t422;
t616 = t86 + t79;
t346 = t349 ^ 2;
t347 = t351 ^ 2;
t407 = t346 / 0.2e1 + t347 / 0.2e1;
t401 = mrSges(6,3) * t407;
t614 = t350 * t401;
t303 = t500 + t501;
t404 = -t274 - t503;
t376 = -t303 + t404;
t304 = -t350 * mrSges(5,2) - t352 * mrSges(5,3);
t613 = m(5) * t297 + t304;
t517 = Ifges(7,4) * t385;
t199 = -Ifges(7,2) * t386 + t517;
t200 = -Ifges(7,1) * t386 - t517;
t409 = t199 / 0.4e1 - t200 / 0.4e1;
t548 = t299 / 0.2e1;
t271 = t303 * t350;
t558 = t271 / 0.2e1;
t196 = mrSges(7,1) * t385 - mrSges(7,2) * t386;
t533 = pkin(5) * t351;
t252 = t334 + (-pkin(4) - t533) * t350;
t540 = t352 / 0.4e1;
t210 = mrSges(7,1) * t352 + t250 * mrSges(7,3);
t566 = -t210 / 0.2e1;
t208 = -mrSges(7,2) * t352 + t248 * mrSges(7,3);
t568 = t208 / 0.2e1;
t150 = -mrSges(7,1) * t250 + mrSges(7,2) * t248;
t574 = t150 / 0.2e1;
t609 = t204 * t566 + t403 * t568 + t252 * t196 / 0.2e1 + t325 * t574 + t448 * t540;
t612 = qJ(4) * t558 + t409 * t250 + t302 * t548 + t609;
t180 = -t268 * t349 - t283;
t181 = t268 * t351 - t470;
t498 = t352 * mrSges(6,2);
t290 = mrSges(6,3) * t461 - t498;
t458 = t351 * t290;
t466 = t349 * t350;
t499 = t352 * mrSges(6,1);
t288 = -mrSges(6,3) * t466 + t499;
t467 = t349 * t288;
t610 = -m(6) * (-t180 * t349 + t181 * t351) - t458 + t467 - t613;
t608 = t404 * qJD(6);
t405 = t224 - t506;
t607 = qJD(6) * t405;
t514 = Ifges(6,6) * t351;
t516 = Ifges(6,5) * t349;
t605 = Ifges(5,6) - t516 / 0.2e1 - t514 / 0.2e1 + Ifges(4,4);
t445 = t346 + t347;
t602 = t445 * t350 * t353;
t227 = Ifges(7,4) * t248;
t142 = -Ifges(7,1) * t250 + t352 * Ifges(7,5) + t227;
t153 = Ifges(7,2) * t250 + t227;
t518 = Ifges(7,4) * t250;
t140 = Ifges(7,2) * t248 + t352 * Ifges(7,6) - t518;
t154 = Ifges(7,1) * t248 + t518;
t399 = (-t154 / 0.4e1 + t140 / 0.4e1) * t385;
t596 = -t399 - (t142 / 0.4e1 + t153 / 0.4e1) * t386;
t485 = t204 * t250;
t486 = t403 * t248;
t595 = t485 / 0.2e1 - t486 / 0.2e1;
t421 = t539 * t248;
t478 = t250 * t348;
t591 = t478 - t421;
t479 = t250 * t386;
t483 = t248 * t385;
t590 = t479 - t483;
t589 = t385 * t84 + t386 * t85;
t324 = Ifges(6,4) * t461;
t272 = -Ifges(6,2) * t466 + t324;
t550 = -t288 / 0.2e1;
t588 = t353 * t550 - t272 / 0.4e1;
t584 = 2 * qJD(3);
t583 = m(6) / 0.2e1;
t582 = -m(7) / 0.2e1;
t581 = m(7) / 0.2e1;
t578 = -mrSges(7,3) / 0.2e1;
t151 = -t504 - t508;
t573 = t151 / 0.2e1;
t569 = -t403 / 0.2e1;
t563 = t248 / 0.2e1;
t560 = -t250 / 0.2e1;
t520 = Ifges(6,4) * t349;
t309 = Ifges(6,1) * t351 - t520;
t273 = t350 * t309;
t556 = t273 / 0.4e1;
t554 = -t385 / 0.2e1;
t553 = t386 / 0.2e1;
t552 = -t386 / 0.2e1;
t549 = t290 / 0.2e1;
t519 = Ifges(6,4) * t351;
t307 = -Ifges(6,2) * t349 + t519;
t547 = -t307 / 0.4e1;
t397 = Ifges(6,1) * t349 + t519;
t546 = -t397 / 0.4e1;
t545 = t349 / 0.2e1;
t543 = -t351 / 0.2e1;
t542 = t351 / 0.2e1;
t541 = t352 / 0.2e1;
t537 = m(7) * t252;
t469 = t348 * t156;
t78 = t148 * t539 - t469;
t532 = t78 * mrSges(7,2);
t531 = t79 * mrSges(7,1);
t528 = t86 * mrSges(7,1);
t87 = t155 * t539 - t469;
t527 = t87 * mrSges(7,2);
t522 = t79 * t249 + t78 * t251;
t68 = t78 * t386;
t521 = t385 * t79 - t68;
t515 = Ifges(6,5) * t352;
t507 = t248 * t78;
t141 = Ifges(7,4) * t251 - Ifges(7,2) * t249 - Ifges(7,6) * t350;
t143 = Ifges(7,1) * t251 - Ifges(7,4) * t249 - Ifges(7,5) * t350;
t152 = mrSges(7,1) * t249 + mrSges(7,2) * t251;
t396 = Ifges(6,2) * t351 + t520;
t244 = -t350 * Ifges(6,6) + t352 * t396;
t246 = -t350 * Ifges(6,5) + t352 * t397;
t253 = -t344 + (t354 - t533) * t352;
t270 = t302 * t352;
t289 = -mrSges(6,1) * t350 - mrSges(6,3) * t465;
t291 = mrSges(6,2) * t350 + mrSges(6,3) * t456;
t497 = t352 * Ifges(6,6);
t243 = t350 * t396 + t497;
t459 = t351 * t243;
t245 = Ifges(6,1) * t466 + t324 + t515;
t468 = t349 * t245;
t3 = -t299 * t270 - t301 * t462 + t186 * t288 + t180 * t289 + t187 * t290 + t181 * t291 + t252 * t152 + t253 * t151 + t143 * t560 + t142 * t559 + t141 * t563 + t140 * t561 + t85 * t208 + t79 * t209 + t84 * t210 + t78 * t211 + t613 * t447 + m(6) * (t180 * t186 + t181 * t187 + t299 * t301) + m(7) * (t252 * t253 + t78 * t84 + t79 * t85) + (t244 * t542 + t246 * t545 + Ifges(7,5) * t250 / 0.2e1 - Ifges(7,6) * t248 / 0.2e1 - qJ(2) * mrSges(4,2) + t297 * mrSges(5,3) + t605 * t350) * t350 + (t459 / 0.2e1 + t468 / 0.2e1 + qJ(2) * mrSges(4,1) - t297 * mrSges(5,2) - t605 * t352 + (-Ifges(6,3) + Ifges(5,3) - Ifges(4,1) + Ifges(4,2) - Ifges(5,2) - Ifges(7,3)) * t350 + t389) * t352;
t502 = t3 * qJD(1);
t323 = Ifges(6,5) * t461;
t449 = Ifges(7,5) * t248 + Ifges(7,6) * t250;
t365 = t252 * t150 + (t142 / 0.2e1 + t153 / 0.2e1) * t248 + (t79 * mrSges(7,3) - t154 / 0.2e1 + t140 / 0.2e1) * t250 - mrSges(7,3) * t507 + t449 * t541;
t6 = m(7) * (t78 * t86 + t79 * t87) + t87 * t208 + t86 * t210 + t180 * t290 - t181 * t288 + t299 * t271 + t323 * t541 + ((-t180 * mrSges(6,3) + t245 / 0.2e1 + t272 / 0.2e1) * t351 + (-t181 * mrSges(6,3) - t497 / 0.2e1 + t273 / 0.2e1 - t243 / 0.2e1 + (t151 + t537) * pkin(5)) * t349) * t350 + t365;
t496 = t6 * qJD(1);
t7 = t78 * t208 - t79 * t210 + t365;
t495 = t7 * qJD(1);
t476 = t251 * t250;
t477 = t251 * t210;
t481 = t249 * t248;
t482 = t249 * t208;
t368 = (-t476 / 0.2e1 - t481 / 0.2e1) * mrSges(7,3) + t482 / 0.2e1 + t477 / 0.2e1;
t388 = t467 / 0.2e1 - t458 / 0.2e1;
t357 = (t388 + t614) * t352 + (t350 ^ 2 * t534 + t249 * t86 - t251 * t87 + t522) * t581 + t368;
t11 = (t574 + t558) * t350 + t357 + t623;
t489 = t11 * qJD(1);
t414 = t150 * t544;
t360 = t414 + t368;
t15 = t360 + t594;
t488 = t15 * qJD(1);
t472 = t386 * t210;
t474 = t385 * t208;
t369 = (t479 / 0.2e1 - t483 / 0.2e1) * mrSges(7,3) + t474 / 0.2e1 - t472 / 0.2e1;
t21 = t369 + t592;
t484 = t21 * qJD(1);
t26 = m(7) * t522 + t352 * t610 + t477 + t482;
t475 = t26 * qJD(1);
t398 = -t350 * mrSges(4,1) - t352 * mrSges(4,2);
t29 = t474 - t472 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * t521 - t398 - t610;
t471 = t29 * qJD(1);
t463 = t350 * t196;
t327 = t350 * t352;
t457 = t351 * t307;
t444 = qJD(3) * t350;
t443 = qJD(3) * t352;
t434 = -t78 / 0.2e1 + t87 / 0.2e1;
t433 = t86 / 0.2e1 + t79 / 0.2e1;
t425 = mrSges(7,3) * t554;
t424 = mrSges(7,3) * t552;
t177 = t251 * t425;
t423 = t463 / 0.2e1 + t177 - t473 * t578;
t416 = t466 / 0.2e1;
t410 = t353 * t549;
t277 = Ifges(7,4) * t386;
t198 = -Ifges(7,2) * t385 - t277;
t201 = Ifges(7,1) * t385 - t277;
t408 = t201 / 0.4e1 + t198 / 0.4e1;
t395 = -t514 - t516;
t61 = m(6) * (-t445 + 0.1e1) * t327 + m(7) * (t327 + t476 + t481);
t387 = -t349 * t291 / 0.2e1 + t289 * t543;
t392 = t186 * t351 + t187 * t349;
t355 = (t573 + t618 + t387) * t352 + (t152 / 0.2e1 + t290 * t545 + t288 * t542 - t270 / 0.2e1) * t350 + ((t299 - t392) * t352 + (t180 * t351 + t181 * t349 + t301) * t350) * t583 + (t249 * t84 - t250 * t79 - t251 * t85 + t252 * t352 + t253 * t350 + t507) * t581 + t210 * t563 + t249 * t565 + t208 * t560 - t251 * t209 / 0.2e1;
t384 = m(7) * (t204 * t385 - t386 * t403);
t9 = -t384 / 0.2e1 + t355;
t394 = t9 * qJD(1) + t61 * qJD(2);
t358 = -t614 + (t385 * t86 + t386 * t87 + t521) * t581 + t369;
t13 = (t549 - t498 / 0.2e1) * t351 + (t550 - t499 / 0.2e1) * t349 + t358 + t621;
t393 = t13 * qJD(1);
t380 = t590 * t581;
t33 = t325 * t196 - (t198 / 0.2e1 + t201 / 0.2e1) * t386 - (t199 / 0.2e1 - t200 / 0.2e1) * t385;
t35 = t423 - t450;
t356 = (mrSges(7,3) * t569 + t408) * t248 + (t204 * mrSges(7,3) / 0.2e1 + t409) * t250 + t596 + t609;
t5 = t356 + t586;
t379 = t5 * qJD(1) + t35 * qJD(2) + t33 * qJD(3);
t112 = t536 + mrSges(5,3) + (m(6) + m(5)) * qJ(4) - t376;
t361 = -m(6) * t392 / 0.2e1 + t589 * t582 + t211 * t554 + t387;
t364 = t177 + mrSges(6,2) * t416 + m(6) * t548 + (t204 * t249 + t403 * t251 + t252) * t581 - t622;
t19 = -(t505 / 0.2e1 + t567) * t386 + t361 + t364;
t73 = -t380 + (t582 + (-0.1e1 / 0.2e1 + t407) * m(6)) * t350;
t378 = -qJD(1) * t19 + qJD(2) * t73 - qJD(3) * t112;
t366 = (t539 * t568 + t348 * t566 + (t478 / 0.2e1 - t421 / 0.2e1) * mrSges(7,3)) * pkin(5);
t17 = -mrSges(7,1) * t433 + mrSges(7,2) * t434 + t366;
t292 = (mrSges(7,1) * t348 + mrSges(7,2) * t539) * pkin(5);
t38 = (t569 + t403 / 0.2e1) * mrSges(7,2);
t377 = -t17 * qJD(1) - t38 * qJD(3) + t292 * qJD(5);
t375 = t616 * t403 + (-t78 + t87) * t204;
t2 = (t201 + t198) * t248 / 0.4e1 - t459 / 0.4e1 - (t153 + t142) * t386 / 0.4e1 + (t309 - t396) * t461 / 0.4e1 - t468 / 0.4e1 - mrSges(6,3) * t602 / 0.2e1 + t595 * mrSges(7,3) + t588 * t349 + t616 * t425 + t533 * t573 - t68 * t578 + ((t252 * t351 + t325 * t466) * pkin(5) + t375) * t581 + t351 * t556 + t612 + t395 * t540 + t466 * t546 + t466 * t547 - t624 + t351 * t410 - pkin(5) * t404 * t416 + t87 * t424 - t399;
t20 = qJ(4) * t302 - t397 * t542 - t457 / 0.2e1 + (-t309 / 0.2e1 + t396 / 0.2e1) * t349 + t33 + (t536 - t404) * t533;
t367 = pkin(5) * t461 * t582 + t618;
t370 = -t591 * t442 + t466 * t579 + t622;
t31 = -t463 / 0.2e1 + t367 + t370;
t374 = t2 * qJD(1) - t31 * qJD(2) + t20 * qJD(3);
t285 = t292 * qJD(6);
t72 = -t380 + (t445 * t583 + m(5)) * t350 + (m(6) + m(7)) * t544;
t36 = t423 + t450;
t32 = -t367 + t370 + t423;
t22 = t369 - t592;
t18 = t249 * t424 + t209 * t553 + (m(5) * t354 - mrSges(5,1)) * t350 - t361 + t364;
t16 = t360 - t594;
t14 = -t532 / 0.2e1 - t531 / 0.2e1 - t527 / 0.2e1 + t528 / 0.2e1 + t366 + t449;
t12 = mrSges(6,1) * t415 + mrSges(6,2) * t412 + t358 - t388 - t621;
t10 = t271 * t544 + t357 + t414 - t623;
t8 = t384 / 0.2e1 + t355;
t4 = t356 - t586;
t1 = (-t385 * t433 - t386 * t434 + t595) * mrSges(7,3) + (-t515 / 0.4e1 - t245 / 0.4e1 + t588) * t349 + t408 * t248 + (t410 - t497 / 0.4e1 + t556 - t243 / 0.4e1 + (t573 + t537 / 0.2e1) * pkin(5)) * t351 + t375 * t581 + ((t309 / 0.4e1 - t396 / 0.4e1) * t351 - t353 * t401 + (t546 + t547 + (t536 / 0.2e1 - t404 / 0.2e1) * pkin(5)) * t349) * t350 + t596 + t612 + t624;
t23 = [qJD(2) * t29 + qJD(3) * t3 + qJD(4) * t26 + qJD(5) * t6 + qJD(6) * t7, qJD(2) * t620 + t8 * qJD(3) + t10 * qJD(5) + t16 * qJD(6) + t471, t502 + t8 * qJD(2) + t18 * qJD(4) + t1 * qJD(5) + t4 * qJD(6) + ((qJ(4) * t301 + t353 * t392) * t583 + (t204 * t85 + t253 * t325 + t403 * t84) * t581) * t584 + (pkin(3) * mrSges(5,1) + Ifges(6,5) * t543 + Ifges(7,5) * t554 + Ifges(6,6) * t545 + Ifges(7,6) * t553 + Ifges(5,4) - Ifges(4,5)) * t444 + (-qJ(4) * mrSges(5,1) + t457 / 0.2e1 + t309 * t545 - Ifges(4,6) + Ifges(5,5)) * t443 + (t325 * t152 + t301 * t303 + t141 * t552 + t385 * t143 / 0.2e1 - qJ(4) * t270 - t253 * t404 + t201 * t559 + t199 * t561 + t204 * t209 + t403 * t211 + (t353 * t289 - t186 * mrSges(6,3) + t246 / 0.2e1) * t351 + (t353 * t291 - t187 * mrSges(6,3) - t244 / 0.2e1) * t349 + (-t304 + t398 + t617) * t354 - t589 * mrSges(7,3)) * qJD(3), t18 * qJD(3) - qJD(4) * t620 + t12 * qJD(5) + t22 * qJD(6) + t475, t496 + t10 * qJD(2) + t1 * qJD(3) + t12 * qJD(4) + ((t348 * t87 + t539 * t86) * t580 - t527 + t528 + t323 - Ifges(6,6) * t466 - t180 * mrSges(6,2) - t181 * mrSges(6,1) + t449 + t591 * mrSges(7,3) * pkin(5)) * qJD(5) + t14 * qJD(6), t495 + t16 * qJD(2) + t4 * qJD(3) + t22 * qJD(4) + t14 * qJD(5) + (t449 - t531 - t532) * qJD(6); qJD(3) * t9 + qJD(5) * t11 + qJD(6) * t15 - t471, t61 * qJD(3) (mrSges(7,3) * t590 - t304) * qJD(3) + t72 * qJD(4) + t32 * qJD(5) + t36 * qJD(6) + (-mrSges(4,2) - t376) * t443 + (-mrSges(6,3) * t445 - mrSges(4,1)) * t444 + ((t325 * t352 - t485 + t486) * t581 + (t337 + t602) * t583 + t617 / 0.2e1) * t584 + t394, qJD(3) * t72, t489 + t32 * qJD(3) + (t352 * t303 + t391 + t405) * qJD(5) + t607, t36 * qJD(3) + qJD(5) * t405 + t488 + t607; -qJD(2) * t9 + qJD(4) * t19 + qJD(5) * t2 + qJD(6) * t5 - t502, -qJD(4) * t73 - qJD(5) * t31 + qJD(6) * t35 - t394, qJD(4) * t112 + qJD(5) * t20 + qJD(6) * t33, -t378 ((-t204 * t539 + t348 * t403) * t580 - mrSges(6,1) * t464 - t353 * t500 + t395 - t599 * mrSges(7,3) + t37) * qJD(5) + t619 + t374, t37 * qJD(5) + t379 + t619; -qJD(3) * t19 + qJD(5) * t13 + qJD(6) * t21 - t475, qJD(3) * t73, t378, 0 (t390 + t376) * qJD(5) + t393 + t608, qJD(5) * t404 + t484 + t608; -qJD(2) * t11 - qJD(3) * t2 - qJD(4) * t13 + qJD(6) * t17 - t496, t31 * qJD(3) - t489, qJD(6) * t38 - t374, -t393, -t285, -t285 - t377; -qJD(2) * t15 - qJD(3) * t5 - qJD(4) * t21 - qJD(5) * t17 - t495, -t35 * qJD(3) - t488, -qJD(5) * t38 - t379, -t484, t377, 0;];
Cq  = t23;
