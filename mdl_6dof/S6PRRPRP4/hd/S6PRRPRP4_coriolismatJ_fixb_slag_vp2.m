% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRPRP4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:30
% EndTime: 2019-03-08 21:40:40
% DurationCPUTime: 6.94s
% Computational Cost: add. (6492->504), mult. (15026->688), div. (0->0), fcn. (13431->8), ass. (0->260)
t329 = cos(qJ(5));
t481 = Ifges(7,4) + Ifges(6,4);
t553 = t481 * t329;
t326 = sin(qJ(5));
t321 = t326 ^ 2;
t323 = t329 ^ 2;
t432 = t321 + t323;
t499 = -t326 / 0.2e1;
t479 = Ifges(6,2) + Ifges(7,2);
t483 = Ifges(6,1) + Ifges(7,1);
t552 = t326 * t483;
t325 = sin(pkin(6));
t327 = sin(qJ(3));
t328 = sin(qJ(2));
t331 = cos(qJ(2));
t439 = t326 * t331;
t178 = (t327 * t439 + t328 * t329) * t325;
t449 = t178 * t326;
t434 = t329 * t331;
t177 = (-t326 * t328 + t327 * t434) * t325;
t450 = t177 * t329;
t538 = mrSges(6,3) + mrSges(7,3);
t551 = t538 * (t450 / 0.2e1 + t449 / 0.2e1);
t330 = cos(qJ(3));
t550 = t432 * t330 / 0.2e1;
t514 = pkin(4) + pkin(8);
t446 = t325 * t328;
t480 = Ifges(6,5) + Ifges(7,5);
t546 = t326 * t480;
t478 = Ifges(6,6) + Ifges(7,6);
t545 = t329 * t478;
t324 = t330 ^ 2;
t544 = t327 ^ 2 + t324;
t392 = -pkin(3) * t327 + qJ(4) * t330;
t379 = m(5) * t392;
t475 = mrSges(5,2) * t327;
t543 = t475 + t379;
t427 = m(7) / 0.4e1 + m(6) / 0.4e1;
t541 = 0.4e1 * t427;
t518 = m(7) * pkin(5);
t419 = mrSges(7,1) + t518;
t332 = -pkin(3) - pkin(9);
t535 = -qJ(6) + t332;
t274 = t535 * t326;
t397 = -qJ(4) * t327 - pkin(2);
t231 = t330 * t332 + t397;
t285 = t514 * t327;
t124 = -t231 * t326 + t285 * t329;
t440 = t326 * t330;
t94 = qJ(6) * t440 + t124;
t91 = pkin(5) * t327 + t94;
t537 = -t274 * t330 - t91;
t459 = cos(pkin(6));
t224 = t327 * t446 - t330 * t459;
t129 = t224 * t329 + t325 * t439;
t435 = t329 * t330;
t462 = t327 * mrSges(7,2);
t266 = -mrSges(7,3) * t435 - t462;
t267 = -t327 * mrSges(6,2) - mrSges(6,3) * t435;
t402 = t267 / 0.2e1 + t266 / 0.2e1;
t536 = t402 * t129;
t398 = t480 * t327;
t463 = t327 * mrSges(7,1);
t262 = mrSges(7,3) * t440 + t463;
t263 = mrSges(6,1) * t327 + mrSges(6,3) * t440;
t534 = t263 + t262;
t533 = t266 + t267;
t460 = t329 * mrSges(7,2);
t282 = t326 * mrSges(7,1) + t460;
t283 = t326 * mrSges(6,1) + t329 * mrSges(6,2);
t532 = t282 + t283;
t316 = t326 * mrSges(7,2);
t531 = t329 * mrSges(7,1) - t316;
t395 = -t479 + t483;
t530 = -t326 * t481 + t329 * t395;
t464 = t327 * mrSges(4,1);
t284 = mrSges(4,2) * t330 + t464;
t506 = t263 / 0.2e1;
t507 = t262 / 0.2e1;
t404 = t506 + t507;
t527 = Ifges(4,4) + Ifges(5,6) - t545;
t494 = t329 / 0.2e1;
t526 = mrSges(6,3) * t550 + t263 * t499 + t267 * t494;
t445 = t325 * t330;
t225 = t327 * t459 + t328 * t445;
t525 = 0.2e1 * t225;
t524 = 2 * qJD(3);
t523 = m(5) / 0.2e1;
t522 = -m(6) / 0.2e1;
t521 = m(6) / 0.2e1;
t520 = -m(7) / 0.2e1;
t519 = m(7) / 0.2e1;
t517 = mrSges(6,1) / 0.2e1;
t516 = mrSges(7,2) / 0.2e1;
t254 = pkin(9) * t327 - t392;
t286 = t514 * t330;
t258 = t329 * t286;
t92 = pkin(5) * t330 + t258 + (-qJ(6) * t327 - t254) * t326;
t515 = t92 / 0.2e1;
t513 = qJ(4) / 0.2e1;
t223 = pkin(5) * t435 + t286;
t512 = t223 / 0.2e1;
t229 = t531 * t330;
t510 = -t229 / 0.2e1;
t441 = t326 * t327;
t260 = mrSges(7,1) * t330 - mrSges(7,3) * t441;
t509 = t260 / 0.2e1;
t508 = -t262 / 0.2e1;
t503 = -t274 / 0.2e1;
t502 = -t284 / 0.2e1;
t501 = t286 / 0.2e1;
t312 = pkin(5) * t326 + qJ(4);
t500 = t312 / 0.2e1;
t498 = -t326 / 0.4e1;
t497 = t326 / 0.2e1;
t495 = -t329 / 0.4e1;
t493 = -t330 / 0.2e1;
t487 = t329 * pkin(5);
t222 = (-t487 - t514) * t327;
t491 = m(7) * t222;
t490 = m(7) * t223;
t489 = m(7) * t312;
t488 = m(7) * t330;
t486 = mrSges(4,2) - mrSges(5,3);
t485 = mrSges(5,2) - mrSges(4,1);
t484 = mrSges(6,2) + mrSges(7,2);
t477 = t91 - t94;
t474 = mrSges(5,3) * t330;
t472 = Ifges(6,4) * t326;
t470 = Ifges(7,4) * t326;
t468 = t223 * mrSges(7,1);
t275 = t535 * t329;
t467 = t275 * mrSges(7,3);
t466 = t286 * mrSges(6,1);
t465 = t286 * mrSges(6,2);
t458 = qJ(4) * t224;
t437 = t329 * t129;
t130 = -t224 * t326 + t325 * t434;
t443 = t326 * t130;
t367 = t437 - t443;
t11 = (-t224 + t367) * t225 * t541;
t456 = t11 * qJD(1);
t444 = t325 * t331;
t417 = t330 * t444;
t150 = t225 * t417;
t357 = (t224 * t327 - t446) * t325;
t12 = (t129 * t177 - t130 * t178 + t150) * t541 + m(5) * (t225 * t445 + t357) * t331 + (t331 * t357 + t150) * m(4);
t455 = t12 * qJD(1);
t126 = -t254 * t326 + t258;
t454 = t126 * t329;
t453 = t129 * t326;
t418 = m(7) * t477;
t360 = t418 / 0.2e1 + t507;
t423 = mrSges(7,1) / 0.2e1 + t517;
t362 = t518 / 0.2e1 + t423;
t401 = -t323 / 0.2e1 - t321 / 0.2e1;
t382 = t401 * mrSges(6,3);
t422 = t516 + mrSges(6,2) / 0.2e1;
t13 = (mrSges(7,3) * t401 + t382) * t330 + (t327 * t422 - t402) * t329 + (t327 * t362 + t360 + t506) * t326;
t452 = t13 * qJD(2);
t451 = t130 * t329;
t447 = t275 * t266;
t442 = t326 * t262;
t438 = t327 * t329;
t436 = t329 * t266;
t127 = t254 * t329 + t286 * t326;
t433 = t544 * pkin(8) * t444;
t430 = qJD(5) * t326;
t429 = qJD(5) * t329;
t253 = (-0.1e1 / 0.2e1 + t401) * m(7);
t428 = t253 * qJD(3);
t424 = t487 / 0.2e1;
t421 = Ifges(6,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t420 = Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1;
t410 = t444 / 0.2e1;
t408 = -t440 / 0.2e1;
t406 = t435 / 0.2e1;
t261 = mrSges(6,1) * t330 - mrSges(6,3) * t441;
t405 = t509 + t261 / 0.2e1;
t264 = -mrSges(7,2) * t330 + mrSges(7,3) * t438;
t265 = -mrSges(6,2) * t330 + mrSges(6,3) * t438;
t403 = -t264 / 0.2e1 - t265 / 0.2e1;
t396 = qJD(5) * t484;
t393 = t432 * t225;
t391 = mrSges(6,1) + t419;
t390 = mrSges(7,3) * pkin(5) - t480;
t385 = t330 * t410;
t384 = t422 * t326;
t383 = t489 / 0.2e1 + t282 / 0.2e1;
t366 = t451 + t453;
t381 = 0.2e1 * t427 * t366 * t327;
t378 = mrSges(6,1) * t329 - mrSges(6,2) * t326;
t230 = t378 * t330;
t96 = qJ(6) * t438 + t127;
t334 = (t510 - t230 / 0.2e1) * t224 + t405 * t129 + t403 * t130 + (t126 * t129 - t127 * t130 - t224 * t286) * t521 + (t129 * t92 - t130 * t96 - t223 * t224) * t519;
t227 = t531 * t327;
t228 = t378 * t327;
t339 = -t227 / 0.2e1 - t228 / 0.2e1 + t404 * t329 + t402 * t326;
t364 = t449 + t450;
t347 = t364 * t522;
t125 = t231 * t329 + t285 * t326;
t349 = m(6) * (t124 * t329 + t125 * t326 - t285);
t95 = -qJ(6) * t435 + t125;
t350 = m(7) * (t326 * t95 + t329 * t91 + t222);
t365 = t177 * t275 + t178 * t274;
t1 = (t502 + t464 / 0.2e1) * t444 + t339 * t225 + t365 * t520 + t332 * t347 + (t349 / 0.4e1 + t350 / 0.4e1) * t525 + (mrSges(4,2) / 0.2e1 - t283 / 0.2e1 + qJ(4) * t522 - t383) * t417 + t334 + t551;
t278 = -pkin(3) * t330 + t397;
t281 = mrSges(5,2) * t330 - mrSges(5,3) * t327;
t5 = pkin(2) * t284 - t222 * t229 + t223 * t227 + t285 * t230 + t286 * t228 + t392 * t281 - t91 * t260 - t92 * t262 - t95 * t264 - t96 * t266 - t124 * t261 - t125 * t265 - t126 * t263 - t127 * t267 - m(6) * (t124 * t126 + t125 * t127 - t285 * t286) - m(7) * (t222 * t223 + t91 * t92 + t95 * t96) + (-t527 + t546) * t324 + (t527 * t327 + (t479 * t323 - Ifges(4,1) + Ifges(4,2) - Ifges(5,2) + Ifges(5,3) - Ifges(6,3) - Ifges(7,3)) * t330 + (-t398 + (0.2e1 * t553 + t552) * t330) * t326) * t327 + (t474 + t543) * t278;
t372 = t1 * qJD(1) - t5 * qJD(2);
t341 = t326 * t362 + t329 * t422;
t342 = t177 * t362 - t178 * t422;
t6 = -t536 + (t477 * t520 - t404) * t130 + (t225 * t341 + t538 * (-t437 / 0.2e1 + t443 / 0.2e1)) * t330 + t342;
t308 = Ifges(7,6) * t440;
t309 = Ifges(6,6) * t440;
t8 = t124 * t267 - t125 * t263 + t94 * t266 + (t309 / 0.2e1 + t308 / 0.2e1) * t327 + (-t418 - t262) * t95 + ((-mrSges(7,2) * t223 + mrSges(6,3) * t124 + t91 * mrSges(7,3) + t330 * t553 - t398 - t465) * t329 + (t95 * mrSges(7,3) + t125 * mrSges(6,3) - t468 - t466 + t421 * t327 + (-t229 - t490) * pkin(5) + t530 * t330) * t326) * t330;
t371 = -qJD(1) * t6 + qJD(2) * t8;
t358 = m(7) * (t326 * t91 - t329 * t95);
t20 = (-m(5) * t278 - t281 - t533 * t329 + t534 * t326 + t358 + m(6) * (t124 * t326 - t125 * t329)) * t327;
t337 = t364 * t520 + t347;
t22 = t381 + t337;
t369 = -qJD(1) * t22 - qJD(2) * t20;
t36 = (t358 - t436 + t442) * t330;
t46 = (t410 - t453 / 0.2e1 - t451 / 0.2e1) * t488;
t368 = -qJD(1) * t46 + qJD(2) * t36;
t363 = t274 * t326 + t275 * t329;
t166 = (t326 * t419 + t460) * t330;
t232 = t329 * t419 - t316;
t361 = qJD(2) * t166 - qJD(3) * t232;
t359 = Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 - Ifges(6,2) / 0.2e1 - Ifges(7,2) / 0.2e1;
t26 = -t312 * t531 - t282 * t487 - qJ(4) * t378 - t481 * t321 + (-pkin(5) * t489 + t326 * t395 + t553) * t329;
t335 = (m(7) * t515 + t509) * pkin(5) + t126 * t517 - t127 * mrSges(6,2) / 0.2e1 + mrSges(7,1) * t515 - t96 * mrSges(7,2) / 0.2e1;
t3 = -t447 / 0.2e1 + t316 * t512 + t360 * t274 + (t332 * t382 + t420) * t330 + (-t466 / 0.2e1 - t332 * t267 / 0.2e1 - t468 / 0.2e1 + t478 * t327 + (-t490 / 0.2e1 + t510) * pkin(5) + (mrSges(6,2) * t513 + mrSges(7,2) * t500 - t467 / 0.2e1 + t359 * t329) * t330) * t329 + (t332 * t506 + t465 / 0.2e1 + t398 + (t94 / 0.2e1 - t91 / 0.2e1) * mrSges(7,3) + (mrSges(6,1) * t513 + mrSges(7,1) * t500 + mrSges(7,3) * t503 + pkin(5) * t383 - t326 * t359 - 0.2e1 * t553) * t330) * t326 + t335;
t348 = -t3 * qJD(2) - t26 * qJD(3);
t346 = t283 * t493;
t345 = t282 * t493;
t340 = (t537 * t329 + (t275 * t330 - t95) * t326) * t519;
t29 = (t508 + t463 / 0.2e1) * t329 + (-t266 / 0.2e1 - t462 / 0.2e1) * t326 + t340 - t491 / 0.2e1;
t43 = 0.2e1 * (-t224 / 0.4e1 + t437 / 0.4e1 - t443 / 0.4e1) * m(7);
t80 = -m(7) * t363 + mrSges(7,3) * t432;
t344 = -qJD(1) * t43 + qJD(2) * t29 + qJD(3) * t80;
t104 = t489 + mrSges(5,3) + t484 * t329 + (mrSges(6,1) + mrSges(7,1)) * t326 + (m(6) + m(5)) * qJ(4);
t336 = m(6) * t501 + ((-t274 * t329 + t275 * t326) * t327 + t223) * t519;
t338 = (t127 * t326 + t454) * t522 + (t326 * t96 + t329 * t92) * t520;
t17 = (t330 * t423 - t405) * t329 + (-t330 * t422 + t403) * t326 + t336 + t338;
t38 = (-t432 + 0.1e1) * t427 * t525;
t343 = qJD(1) * t38 + qJD(2) * t17 + qJD(3) * t104;
t252 = t432 * t520 + t519;
t47 = t366 * t488 / 0.2e1 + m(7) * t385;
t44 = (t224 + t367) * t520;
t33 = (t523 + t427) * t525 + (m(7) + m(6)) * t393 / 0.2e1;
t28 = t340 + t329 * t508 + t266 * t499 + t491 / 0.2e1 + t441 * t516 - mrSges(7,1) * t438 / 0.2e1;
t21 = m(5) * t327 * t444 - t337 + t381;
t15 = (m(5) * pkin(8) + t329 * t423 + mrSges(5,1) - t384) * t330 + t336 - t338 + (t265 + t264) * t497 + (t261 + t260) * t494;
t14 = -t418 * t497 - t442 / 0.2e1 + t436 / 0.2e1 + t341 * t327 + t550 * mrSges(7,3) + t526;
t10 = (t378 / 0.2e1 + t531 / 0.2e1 + m(7) * t424 + t329 * t362 - t384) * t225;
t7 = t342 + (-pkin(5) * t440 * t519 + t345 + t346) * t225 + (t477 * t519 + t404) * t130 + t536 + t538 * (t129 * t406 + t130 * t408);
t4 = t312 * t345 + qJ(4) * t346 + (-t477 * t274 + (t223 * t329 - t312 * t440) * pkin(5)) * t519 + t531 * t512 + t378 * t501 + t262 * t503 + t406 * t467 + pkin(5) * t282 * t408 + t335 + t420 * t330 + t447 / 0.2e1 + t229 * t424 + (-t545 - t546) * t327 / 0.4e1 + (t472 / 0.2e1 + t470 / 0.2e1 - t483 * t329 / 0.2e1) * t435 + (t421 * t329 + (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t326 + t480 * t498 + t478 * t495) * t327 + t526 * t332 + 0.2e1 * ((-t553 - t552) * t498 + (-t329 * t479 - t470 - t472) * t495) * t330 + (t94 + t537) * mrSges(7,3) * t499 + (-t479 * t326 + t553) * t440 / 0.2e1;
t2 = (qJ(4) * t417 + t332 * t364) * t521 + (t349 / 0.2e1 + t350 / 0.2e1 + t339) * t225 + t334 + (t312 * t417 + t365) * t519 + (t502 + t475 / 0.2e1 + t474 / 0.2e1) * t444 + t543 * t410 - (-t379 + t284) * t444 / 0.2e1 + (mrSges(5,3) + t532) * t385 - t551;
t9 = [qJD(2) * t12 + qJD(3) * t11, t2 * qJD(3) + t21 * qJD(4) + t7 * qJD(5) + t47 * qJD(6) + t455 + (t534 * t177 + t533 * t178 + ((-mrSges(4,1) * t330 + mrSges(4,2) * t327 - mrSges(3,1) + t281) * t328 + (-mrSges(3,2) + (t229 + t230) * t330 + (mrSges(4,3) + mrSges(5,1)) * t544) * t331) * t325 + 0.2e1 * (t177 * t91 + t178 * t95 + t223 * t417) * t519 + 0.2e1 * (t124 * t177 + t125 * t178 + t286 * t417) * t521 + 0.2e1 * (t278 * t446 + t433) * t523 + m(4) * (-pkin(2) * t446 + t433)) * qJD(2), t456 + t2 * qJD(2) + t33 * qJD(4) + t10 * qJD(5) + t44 * qJD(6) + ((-t224 * t312 + t225 * t363) * t519 + (t332 * t393 - t458) * t521 + (-pkin(3) * t225 - t458) * t523) * t524 + ((t486 - t532) * t224 + (-t432 * t538 + t485) * t225) * qJD(3), qJD(2) * t21 + qJD(3) * t33, qJD(5) * t130 * t391 + t7 * qJD(2) + t10 * qJD(3) - t129 * t396, qJD(2) * t47 + qJD(3) * t44; qJD(3) * t1 + qJD(4) * t22 - qJD(5) * t6 - qJD(6) * t46 - t455, -qJD(3) * t5 + qJD(4) * t20 + qJD(5) * t8 + qJD(6) * t36, t15 * qJD(4) + t4 * qJD(5) + t28 * qJD(6) + ((-qJ(4) * t285 + t332 * t454) * t521 + (t222 * t312 + t274 * t96 + t275 * t92) * t519) * t524 + t372 + (-qJ(4) * t228 + t222 * t282 - t312 * t227 + t275 * t260 + t274 * t264 - t285 * t283 + (-t126 * mrSges(6,3) - t92 * mrSges(7,3) + t332 * t261) * t329 + (-pkin(3) * mrSges(5,1) - Ifges(5,4) + Ifges(4,5) + t480 * t329 + (-m(5) * pkin(3) + t485) * pkin(8)) * t330 + (-qJ(4) * mrSges(5,1) + Ifges(5,5) - Ifges(4,6) + t481 * t323 + (-m(5) * qJ(4) + t486) * pkin(8)) * t327 + (-t96 * mrSges(7,3) + t332 * t265 - t478 * t330 + (m(6) * t332 - mrSges(6,3)) * t127 + t530 * t327) * t326) * qJD(3), qJD(3) * t15 + qJD(5) * t14 - t369, t4 * qJD(3) + t14 * qJD(4) + (-mrSges(6,1) * t125 - mrSges(6,2) * t124 - mrSges(7,2) * t94 - t419 * t95 + t308 + t309) * qJD(5) + t390 * t330 * t429 + t371, qJD(3) * t28 + t368; -qJD(2) * t1 + qJD(4) * t38 - qJD(6) * t43 - t456, qJD(4) * t17 - qJD(5) * t3 + qJD(6) * t29 - t372, qJD(4) * t104 - qJD(5) * t26 + qJD(6) * t80, qJD(6) * t252 + t343 (-mrSges(7,2) * t275 - t274 * t419) * qJD(5) + (-mrSges(6,2) * t332 - t478) * t429 + (-mrSges(6,1) * t332 + t390) * t430 + t348, qJD(4) * t252 + t344; -qJD(2) * t22 - qJD(3) * t38, -qJD(3) * t17 - qJD(5) * t13 + t369, qJD(6) * t253 - t343, 0, -t329 * t396 - t391 * t430 - t452, t428; qJD(2) * t6, qJD(3) * t3 + qJD(4) * t13 + qJD(6) * t166 - t371, -qJD(6) * t232 - t348, t452, 0, t361; qJD(2) * t46 + qJD(3) * t43, -qJD(3) * t29 - qJD(5) * t166 - t368, -qJD(4) * t253 + qJD(5) * t232 - t344, -t428, -t361, 0;];
Cq  = t9;
