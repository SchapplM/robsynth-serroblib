% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-03-08 22:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPP3_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP3_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP3_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP3_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP3_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP3_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP3_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:54:42
% EndTime: 2019-03-08 22:54:54
% DurationCPUTime: 7.51s
% Computational Cost: add. (7572->667), mult. (18833->871), div. (0->0), fcn. (16773->8), ass. (0->329)
t602 = m(6) + m(7);
t606 = m(5) + t602;
t605 = Ifges(7,4) + Ifges(6,5);
t578 = m(5) / 0.2e1;
t576 = m(6) / 0.2e1;
t574 = m(7) / 0.2e1;
t595 = Ifges(5,5) + Ifges(7,5);
t377 = sin(qJ(4));
t371 = t377 ^ 2;
t380 = cos(qJ(4));
t373 = t380 ^ 2;
t604 = t373 + t371;
t378 = sin(qJ(3));
t381 = cos(qJ(3));
t539 = pkin(9) * t381;
t332 = pkin(3) * t378 - t539;
t288 = t377 * t332;
t541 = pkin(8) * t380;
t451 = -qJ(5) + t541;
t167 = t451 * t378 - t288;
t496 = t377 * t381;
t128 = -pkin(5) * t496 - t167;
t501 = t332 * t380;
t542 = pkin(8) * t377;
t174 = -t501 + (-pkin(4) - t542) * t378;
t497 = t377 * t378;
t207 = pkin(8) * t497 + t501;
t495 = t378 * t380;
t208 = -pkin(8) * t495 + t288;
t376 = pkin(4) + qJ(6);
t294 = mrSges(6,1) * t496 - mrSges(6,3) * t378;
t523 = t378 * mrSges(7,2);
t295 = -mrSges(7,1) * t496 + t523;
t455 = t294 / 0.2e1 - t295 / 0.2e1;
t494 = t380 * t381;
t296 = mrSges(6,1) * t494 + t378 * mrSges(6,2);
t557 = t296 / 0.2e1;
t293 = mrSges(7,1) * t494 - t378 * mrSges(7,3);
t558 = t293 / 0.2e1;
t568 = -t128 / 0.2e1;
t450 = qJ(6) + t542;
t96 = (pkin(5) * t381 - t332) * t380 + (-pkin(4) - t450) * t378;
t570 = t96 / 0.2e1;
t572 = mrSges(5,2) / 0.2e1;
t575 = -m(7) / 0.2e1;
t577 = -m(6) / 0.2e1;
t603 = (-pkin(4) * t174 - qJ(5) * t167) * t577 + (qJ(5) * t128 - t376 * t96) * t575 + pkin(4) * t557 + mrSges(7,2) * t568 + t167 * mrSges(6,3) / 0.2e1 - t174 * mrSges(6,2) / 0.2e1 - t207 * mrSges(5,1) / 0.2e1 + t208 * t572 + t376 * t558 + mrSges(7,3) * t570 + t455 * qJ(5);
t601 = -t378 / 0.2e1;
t549 = t378 / 0.2e1;
t600 = -t381 / 0.2e1;
t599 = t381 / 0.2e1;
t598 = mrSges(5,1) + mrSges(7,3);
t597 = mrSges(6,1) + mrSges(5,3);
t596 = mrSges(6,3) + mrSges(7,2);
t305 = -t380 * mrSges(7,2) + t377 * mrSges(7,3);
t268 = t305 * t378;
t310 = -t377 * mrSges(6,2) - t380 * mrSges(6,3);
t273 = t310 * t378;
t592 = t268 + t273;
t306 = t380 * mrSges(6,2) - t377 * mrSges(6,3);
t430 = mrSges(7,2) * t377 + mrSges(7,3) * t380;
t591 = t306 - t430;
t590 = Ifges(6,4) * t497 + t605 * t495;
t360 = Ifges(7,6) * t377;
t589 = Ifges(7,3) * t380 + t360;
t318 = -Ifges(7,2) * t380 + t360;
t365 = Ifges(5,4) * t380;
t588 = -Ifges(5,2) * t377 + t365;
t323 = Ifges(5,1) * t377 + t365;
t446 = t572 - mrSges(6,3) / 0.2e1 - mrSges(7,2) / 0.2e1;
t571 = mrSges(6,2) / 0.2e1;
t473 = -mrSges(7,3) / 0.2e1 + t571;
t448 = mrSges(5,1) / 0.2e1 - t473;
t587 = t448 * t377 + t446 * t380;
t522 = t381 * mrSges(4,2);
t312 = t378 * mrSges(4,1) + t522;
t585 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t527 = Ifges(6,6) * t380;
t316 = Ifges(6,3) * t377 - t527;
t526 = Ifges(7,6) * t380;
t317 = Ifges(7,2) * t377 + t526;
t472 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t444 = Ifges(5,6) / 0.2e1 + t472;
t584 = -t378 * t444 + Ifges(5,6) * t601 + t588 * t600 + (t316 + t317) * t599 + t605 * t549;
t538 = t381 * pkin(3);
t569 = pkin(5) + pkin(9);
t402 = (-t569 * t378 - pkin(2) - t538) * t377;
t126 = t451 * t381 + t402;
t375 = sin(pkin(6));
t379 = sin(qJ(2));
t500 = t375 * t379;
t514 = cos(pkin(6));
t261 = t514 * t378 + t381 * t500;
t382 = cos(qJ(2));
t493 = t380 * t382;
t135 = t261 * t377 + t375 * t493;
t499 = t375 * t382;
t470 = t377 * t499;
t136 = t261 * t380 - t470;
t510 = qJ(5) * t380;
t423 = qJ(6) * t377 - t510;
t368 = t378 * pkin(8);
t483 = pkin(4) * t497 + t368;
t158 = t378 * t423 + t483;
t369 = t381 * pkin(8);
t482 = pkin(4) * t496 + t369;
t159 = t381 * t423 + t482;
t441 = -pkin(3) * t377 + t541;
t467 = -pkin(9) * t378 - pkin(2);
t443 = t377 * t467;
t164 = -t443 + (qJ(5) - t441) * t381;
t416 = t467 - t538;
t287 = t380 * t416;
t202 = pkin(8) * t496 - t287;
t370 = t381 * pkin(4);
t165 = t202 + t370;
t476 = pkin(8) * t494;
t203 = t377 * t416 + t476;
t233 = -qJ(5) * t495 + t483;
t234 = -qJ(5) * t494 + t482;
t260 = t378 * t500 - t514 * t381;
t271 = t305 * t381;
t311 = t377 * mrSges(5,1) + t380 * mrSges(5,2);
t274 = t311 * t378;
t275 = t310 * t381;
t276 = t311 * t381;
t290 = -mrSges(5,2) * t378 - mrSges(5,3) * t496;
t292 = mrSges(5,1) * t378 - mrSges(5,3) * t494;
t359 = t381 * mrSges(7,3);
t474 = mrSges(7,1) * t495;
t297 = t359 + t474;
t291 = -mrSges(5,1) * t381 - mrSges(5,3) * t495;
t300 = mrSges(6,1) * t495 - mrSges(6,2) * t381;
t456 = t291 / 0.2e1 - t300 / 0.2e1;
t434 = -t297 / 0.2e1 + t456;
t521 = t381 * mrSges(7,2);
t299 = -mrSges(7,1) * t497 - t521;
t289 = mrSges(5,2) * t381 - mrSges(5,3) * t497;
t520 = t381 * mrSges(6,3);
t298 = mrSges(6,1) * t497 + t520;
t457 = -t289 / 0.2e1 + t298 / 0.2e1;
t435 = -t299 / 0.2e1 + t457;
t559 = t273 / 0.2e1;
t562 = t268 / 0.2e1;
t458 = t559 + t562;
t477 = pkin(5) * t495;
t415 = t287 - t370 - t477;
t94 = t450 * t381 - t415;
t583 = (t274 / 0.2e1 + t458) * t261 + (-t292 / 0.2e1 + t558 + t557) * t135 + (t290 / 0.2e1 - t455) * t136 + (-t135 * t207 + t136 * t208 + t261 * t368) * t578 + (t135 * t174 - t136 * t167 + t233 * t261) * t576 + (t128 * t136 + t135 * t96 + t158 * t261) * t574 + ((-t202 * t377 - t203 * t380 + t369) * t578 + (t164 * t380 - t165 * t377 + t234) * t576 + (-t126 * t380 - t377 * t94 + t159) * t574 + t275 / 0.2e1 + t276 / 0.2e1 + t271 / 0.2e1 + t435 * t380 + t434 * t377) * t260;
t582 = 0.2e1 * m(7);
t581 = -0.2e1 * qJ(5);
t580 = 2 * qJD(3);
t579 = 2 * qJD(4);
t573 = -mrSges(7,1) / 0.2e1;
t139 = t402 + t476;
t567 = t139 / 0.2e1;
t201 = (t377 * t379 + t381 * t493) * t375;
t565 = t201 / 0.2e1;
t564 = t203 / 0.2e1;
t269 = t378 * t306;
t561 = -t269 / 0.2e1;
t272 = t430 * t378;
t560 = t272 / 0.2e1;
t556 = t305 / 0.2e1;
t555 = t306 / 0.2e1;
t554 = -t430 / 0.2e1;
t553 = t310 / 0.2e1;
t552 = t311 / 0.2e1;
t330 = t569 * t377;
t551 = -t330 / 0.2e1;
t511 = qJ(5) * t377;
t426 = -pkin(4) * t380 - t511;
t303 = -pkin(3) + t426;
t547 = m(6) * t303;
t367 = t377 * pkin(4);
t308 = t367 - t510;
t546 = m(6) * t308;
t138 = -t202 - t477;
t545 = m(7) * t138;
t266 = -t376 * t380 - pkin(3) - t511;
t544 = m(7) * t266;
t286 = t367 + t423;
t543 = m(7) * t286;
t537 = mrSges(6,1) + mrSges(7,1);
t535 = m(7) * qJD(5);
t534 = m(7) * qJD(6);
t533 = mrSges(4,2) * t378;
t531 = Ifges(5,4) * t377;
t528 = Ifges(6,6) * t377;
t519 = t381 * Ifges(6,4);
t518 = t381 * Ifges(5,6);
t517 = t138 + t94;
t431 = t380 * mrSges(5,1) - t377 * mrSges(5,2);
t515 = -t431 - mrSges(4,1);
t513 = qJ(5) * t135;
t512 = qJ(5) * t201;
t502 = t260 * t380;
t503 = t260 * t377;
t11 = t606 * (-t135 * t503 - t136 * t502 + t260 * t261);
t509 = t11 * qJD(1);
t200 = -t380 * t500 + t381 * t470;
t469 = t378 * t499;
t12 = m(4) * (t260 * t378 + t261 * t381 - t500) * t499 + t606 * (t135 * t200 + t136 * t201 + t260 * t469);
t508 = t12 * qJD(1);
t505 = t200 * t377;
t504 = t201 * t380;
t498 = t376 * t377;
t492 = t381 * t135;
t491 = -t126 + t139;
t490 = t164 + t203;
t489 = t165 - t202;
t488 = t604 * pkin(9) * t260;
t486 = t298 - t299;
t485 = t300 - t291;
t267 = pkin(4) * t495 + qJ(5) * t497;
t481 = qJD(4) * t135;
t480 = t576 + t574;
t479 = m(6) / 0.4e1 + m(7) / 0.4e1;
t471 = Ifges(5,5) / 0.2e1 + Ifges(7,5) / 0.2e1;
t468 = t260 * t497;
t461 = t430 * t549;
t454 = t555 + t554;
t452 = -t430 + t544;
t447 = -mrSges(6,1) / 0.2e1 + t573 - mrSges(5,3) / 0.2e1;
t445 = Ifges(6,4) / 0.2e1 - t471;
t442 = -mrSges(6,1) * pkin(4) - mrSges(7,1) * t376;
t440 = 0.2e1 * t479 * t200;
t439 = -t537 * qJ(5) - Ifges(5,6);
t221 = t378 * t588 - t518;
t228 = -Ifges(6,5) * t381 + t378 * t316;
t344 = Ifges(7,6) * t495;
t230 = -Ifges(7,4) * t381 + Ifges(7,2) * t497 + t344;
t438 = -t221 / 0.2e1 + t228 / 0.2e1 + t230 / 0.2e1;
t324 = Ifges(5,1) * t380 - t531;
t223 = -t381 * Ifges(5,5) + t378 * t324;
t226 = -t381 * Ifges(7,5) + t378 * t589;
t345 = Ifges(6,6) * t497;
t232 = -Ifges(6,2) * t495 + t345 - t519;
t436 = t223 / 0.2e1 + t226 / 0.2e1 - t232 / 0.2e1;
t314 = Ifges(7,3) * t377 - t526;
t428 = Ifges(6,2) * t377 + t527;
t433 = t314 / 0.2e1 + t428 / 0.2e1 + t323 / 0.2e1;
t315 = -Ifges(6,3) * t380 - t528;
t321 = Ifges(5,2) * t380 + t531;
t432 = t315 / 0.2e1 + t318 / 0.2e1 - t321 / 0.2e1;
t320 = -Ifges(6,2) * t380 + t528;
t331 = t569 * t380;
t406 = (t504 + t505) * pkin(9);
t386 = -m(5) * (-pkin(3) * t469 + t406) / 0.2e1 + (t303 * t469 + t406) * t577 + (t330 * t200 + t331 * t201 + t266 * t469) * t575;
t411 = t380 * t447;
t412 = t377 * t447;
t4 = t201 * t411 + t200 * t412 + (-t312 / 0.2e1 + t522 / 0.2e1 + (mrSges(4,1) / 0.2e1 + t431 / 0.2e1 - t454) * t378) * t499 + t386 + t583;
t400 = Ifges(6,4) * t601 + t320 * t600 - t445 * t378 + (t324 + t589) * t599 + t595 * t549;
t5 = -pkin(2) * t312 - t202 * t292 + t94 * t293 + t164 * t294 + t126 * t295 + t165 * t296 + t96 * t297 + t167 * t298 + t128 * t299 + t174 * t300 + t208 * t289 + t203 * t290 + t207 * t291 + t234 * t273 + t233 * t275 + t159 * t268 + t158 * t271 + m(5) * (-t202 * t207 + t203 * t208) + m(6) * (t164 * t167 + t165 * t174 + t233 * t234) + m(7) * (t126 * t128 + t158 * t159 + t94 * t96) + (pkin(8) * t274 + Ifges(4,4) * t381 + (t381 * t445 + t436) * t380 + (t381 * t444 + t438) * t377) * t381 + (-Ifges(4,4) * t378 + pkin(8) * t276 + t400 * t380 + t584 * t377 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - t585) * t381) * t378;
t425 = t4 * qJD(1) + t5 * qJD(2);
t185 = qJ(6) * t495 + t267;
t270 = t431 * t378;
t384 = (t135 * t412 + t136 * t411) * t378 + (t561 + t270 / 0.2e1 + t560) * t260 - t434 * t136 + t435 * t135 + (t490 * t135 + t489 * t136 + t260 * t267) * t576 + (t491 * t135 + t517 * t136 + t185 * t260) * t574;
t396 = (-pkin(4) * t200 + t512) * t577 + (-t200 * t376 + t512) * t575;
t6 = t448 * t200 + t446 * t201 + t384 + t396;
t277 = -Ifges(7,3) * t497 + t344;
t278 = Ifges(6,3) * t495 + t345;
t279 = t378 * t318;
t280 = t428 * t378;
t281 = t378 * t321;
t282 = t378 * t323;
t8 = t139 * t297 + t138 * t299 + t267 * t273 + t185 * t268 - t233 * t269 + t158 * t272 + t485 * t203 + (t298 - t289) * t202 + m(6) * (t164 * t202 + t165 * t203 + t233 * t267) + m(7) * (t126 * t138 + t139 * t94 + t158 * t185) + (pkin(8) * t270 + (t277 / 0.2e1 - t280 / 0.2e1 - t282 / 0.2e1 + t518 / 0.2e1 - t126 * mrSges(7,1) + t164 * mrSges(6,1) - t203 * mrSges(5,3) + t438) * t380 + (t278 / 0.2e1 - t279 / 0.2e1 + t281 / 0.2e1 - t94 * mrSges(7,1) - t165 * mrSges(6,1) - t202 * mrSges(5,3) + t471 * t381 - t436) * t377) * t378 + t590 * t600;
t424 = t6 * qJD(1) + t8 * qJD(2);
t24 = t486 * t381 - t592 * t495 + m(6) * (t164 * t381 - t233 * t495) + m(7) * (-t126 * t381 - t158 * t495);
t393 = (t577 + t575) * (-t136 * t381 - t260 * t495);
t26 = t440 + t393;
t422 = -qJD(1) * t26 + qJD(2) * t24;
t33 = t268 * t497 + m(7) * (t158 * t497 + t381 * t94) + t381 * t297;
t37 = (t201 / 0.4e1 - t468 / 0.4e1 - t492 / 0.4e1) * t582;
t421 = -qJD(1) * t37 + qJD(2) * t33;
t391 = (-t233 * t377 + (-t303 * t378 - t539) * t380) * t577 + (-t158 * t377 - t266 * t495 - t331 * t381) * t575;
t414 = m(7) * t570 + t174 * t576;
t17 = t537 * t494 + t458 * t377 + (t454 * t380 + t473) * t378 + t391 + t414;
t49 = (t306 + t452 + t547) * t377;
t420 = qJD(2) * t17 + qJD(3) * t49;
t122 = t452 * t380;
t397 = (-t158 * t380 + t266 * t497 + t330 * t381) * t574 - t380 * t268 / 0.2e1;
t408 = m(7) * t568 - t523 / 0.2e1;
t29 = (t381 * mrSges(7,1) - t461) * t377 + t397 + t408;
t419 = -qJD(2) * t29 + qJD(3) * t122;
t388 = (t443 + (t581 + t441) * t381) * t577 + ((t581 + t541) * t381 + t402) * t575;
t413 = m(6) * t564 + m(7) * t567;
t32 = t388 + t413 + t520 + t521;
t336 = qJ(5) * t602 + t596;
t418 = -qJD(2) * t32 + qJD(4) * t336;
t335 = m(7) * t376 + mrSges(7,3);
t394 = -t359 + ((-t376 - t450) * t381 + t415) * t574;
t41 = -t545 / 0.2e1 + t394;
t417 = qJD(2) * t41 + qJD(4) * t335;
t361 = Ifges(7,5) * t380;
t362 = Ifges(6,5) * t377;
t363 = Ifges(5,5) * t380;
t364 = Ifges(7,4) * t377;
t383 = (-t364 / 0.4e1 - t361 / 0.4e1 - t363 / 0.4e1 - t362 / 0.4e1) * t381 + (t233 * t308 + t267 * t303) * t576 + (t158 * t286 + t185 * t266 + t491 * t330 + t517 * t331) * t574 - pkin(3) * t270 / 0.2e1 + t158 * t556 + t185 * t554 + t233 * t553 + t266 * t560 + t267 * t555 + t286 * t562 + t303 * t561 + t308 * t559 + t299 * t551 + t331 * t297 / 0.2e1;
t389 = (t567 - t126 / 0.2e1) * mrSges(7,1) + (t564 + t164 / 0.2e1) * mrSges(6,1) + (t490 * t576 + t457) * pkin(9) - t221 / 0.4e1 + t228 / 0.4e1 + t230 / 0.4e1 + t277 / 0.4e1 - t280 / 0.4e1 - t282 / 0.4e1;
t390 = (t138 / 0.2e1 + t94 / 0.2e1) * mrSges(7,1) + (t165 / 0.2e1 - t202 / 0.2e1) * mrSges(6,1) + (t489 * t576 - t456) * pkin(9) + t223 / 0.4e1 + t226 / 0.4e1 - t232 / 0.4e1 - t278 / 0.4e1 + t279 / 0.4e1 - t281 / 0.4e1;
t395 = pkin(8) * t552 + t597 * pkin(9) * (-t373 / 0.2e1 - t371 / 0.2e1);
t398 = t324 / 0.4e1 - t321 / 0.4e1 - t320 / 0.4e1 + t318 / 0.4e1 + t315 / 0.4e1 + t589 / 0.4e1 + t331 * t573;
t399 = -t323 / 0.4e1 - t588 / 0.4e1 - t428 / 0.4e1 + t317 / 0.4e1 + t316 / 0.4e1 - t314 / 0.4e1 + mrSges(7,1) * t551;
t1 = (-Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1 - Ifges(5,3) / 0.2e1 + t395) * t378 + t383 + ((0.3e1 / 0.4e1 * Ifges(6,4) - t471) * t381 + t398 * t378 + t390) * t380 + ((0.3e1 / 0.4e1 * Ifges(5,6) + t472) * t381 + t399 * t378 + t389) * t377 + t603;
t13 = -pkin(3) * t311 - t286 * t430 + t308 * t306 + (t310 + t546) * t303 + (t305 + t543) * t266 + (-t316 / 0.2e1 + t588 / 0.2e1 - t317 / 0.2e1 + t433) * t380 + (-t320 / 0.2e1 + t324 / 0.2e1 + t589 / 0.2e1 + t432) * t377;
t204 = qJ(5) * t502;
t9 = -0.2e1 * t479 * t204 + (-t305 / 0.2e1 - t311 / 0.2e1 - t310 / 0.2e1 + (t498 / 0.4e1 - t286 / 0.4e1) * t582 + 0.2e1 * (t367 / 0.4e1 - t308 / 0.4e1) * m(6) + t587) * t260;
t403 = t9 * qJD(1) - t1 * qJD(2) - t13 * qJD(3);
t374 = t381 ^ 2;
t372 = t378 ^ 2;
t325 = (qJD(2) * t381 - qJD(4)) * m(7);
t304 = t372 * pkin(8) * t499;
t235 = -m(7) * t330 - t377 * mrSges(7,1);
t137 = m(7) * t331 + (m(6) * pkin(9) + t537) * t380;
t80 = m(7) * t502;
t47 = 0.2e1 * t480 * t136;
t39 = (-t480 - t602 / 0.2e1) * t503;
t38 = (t468 + t492) * t574 + m(7) * t565;
t36 = -t474 + t545 / 0.2e1 + t394;
t31 = -t377 * t461 + t397 - t408;
t30 = -t381 * t596 - t537 * t497 - t388 + t413;
t27 = t440 - t393;
t19 = t473 * t378 - t391 + t414 - t592 * t377 / 0.2e1 - t591 * t495 / 0.2e1;
t10 = (pkin(4) * t503 - t204) * t576 - t204 * t574 + (t546 + t543) * t260 / 0.2e1 + (t498 * t574 + t552 + t553 + t556 + t587) * t260;
t7 = -t201 * mrSges(5,2) / 0.2e1 + t384 - t396 + t596 * t565 + (t571 - t598 / 0.2e1) * t200;
t3 = -t386 - t312 * t499 + (-t431 + t591) * t469 / 0.2e1 + t583 + (t505 / 0.2e1 + t504 / 0.2e1) * (mrSges(7,1) + t597);
t2 = t383 + (t377 * t399 + t380 * t398 + t395) * t378 + (t518 / 0.4e1 + t389) * t377 + (t519 / 0.4e1 + t390) * t380 + t585 * t549 + (-Ifges(6,4) / 0.2e1 + t595 / 0.2e1) * t494 - t603 + (-Ifges(5,6) / 0.2e1 + t605 / 0.2e1) * t496;
t14 = [qJD(2) * t12 + qJD(3) * t11, t3 * qJD(3) + t7 * qJD(4) + t27 * qJD(5) + t38 * qJD(6) + t508 + ((t289 - t486) * t201 + (t297 + t485) * t200 + ((-mrSges(4,1) * t381 - mrSges(3,1) + t533) * t379 + (-mrSges(3,2) + (t372 + t374) * mrSges(4,3) + (t274 + t592) * t378) * t382) * t375 + 0.2e1 * (-t164 * t201 + t165 * t200 + t233 * t469) * t576 + 0.2e1 * (t126 * t201 + t158 * t469 + t94 * t200) * t574 + 0.2e1 * (t200 * t202 + t201 * t203 + t304) * t578 + m(4) * (t304 + (pkin(8) * t374 * t382 - pkin(2) * t379) * t375)) * qJD(2), t509 + t3 * qJD(2) + t10 * qJD(4) + t39 * qJD(5) - t80 * qJD(6) + ((t261 * t303 - t488) * t576 + t261 * t544 / 0.2e1 + (-pkin(3) * t261 - t488) * t578) * t580 + ((t515 + t591) * t261 + (m(7) * (-t330 * t377 - t331 * t380) + mrSges(4,2) + t604 * (-mrSges(5,3) - t537)) * t260) * qJD(3), t7 * qJD(2) + t10 * qJD(3) + t47 * qJD(5) + (mrSges(6,2) - t598) * qJD(4) * t136 + (mrSges(5,2) - t596) * t481 + ((-pkin(4) * t136 - t513) * t576 + (-t136 * t376 - t513) * t574) * t579 - t135 * t534, qJD(2) * t27 + qJD(3) * t39 + qJD(4) * t47, -m(7) * t481 + t38 * qJD(2) - t80 * qJD(3); qJD(3) * t4 + qJD(4) * t6 - qJD(5) * t26 - qJD(6) * t37 - t508, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t24 + qJD(6) * t33, t2 * qJD(4) + t19 * qJD(5) + t31 * qJD(6) + ((t128 * t331 + t159 * t266 + t330 * t96) * t574 + t234 * t547 / 0.2e1) * t580 + t425 + (-Ifges(4,6) * t378 - pkin(3) * t276 + pkin(8) * t533 - t159 * t430 + t234 * t306 + t266 * t271 + t303 * t275 + t330 * t293 + t331 * t295 + (-t167 * mrSges(6,1) + t128 * mrSges(7,1) + t208 * mrSges(5,3) - t584) * t380 + (t174 * mrSges(6,1) + t96 * mrSges(7,1) - t207 * mrSges(5,3) + t400) * t377 + ((t290 - t294) * t380 + (-t292 + t296) * t377 + m(5) * (-t207 * t377 + t208 * t380) + m(6) * (-t167 * t380 + t174 * t377)) * pkin(9) + (Ifges(4,5) + t433 * t380 + t432 * t377 + (-m(5) * pkin(3) + t515) * pkin(8)) * t381) * qJD(3), t2 * qJD(3) + t30 * qJD(5) + t36 * qJD(6) + ((-pkin(4) * t203 - qJ(5) * t202) * t576 + (qJ(5) * t138 - t139 * t376) * t574) * t579 + t424 + (t138 * mrSges(7,2) - t139 * mrSges(7,3) + (t439 * t380 + (-t442 - t595) * t377) * t378 + (-mrSges(5,1) + mrSges(6,2)) * t203 + (mrSges(5,2) - mrSges(6,3)) * t202 + t590) * qJD(4), qJD(3) * t19 + qJD(4) * t30 + t422, qJD(3) * t31 + qJD(4) * t36 + t421; -qJD(2) * t4 - qJD(4) * t9 - t509, qJD(4) * t1 - qJD(5) * t17 + qJD(6) * t29 - t425, qJD(4) * t13 - qJD(5) * t49 - qJD(6) * t122, t137 * qJD(5) + t235 * qJD(6) - t403 + (m(7) * (-qJ(5) * t330 - t331 * t376) + t363 + t364 + t361 + t362 - t330 * mrSges(7,2) - t331 * mrSges(7,3) + (-Ifges(6,4) + t442) * t380 + t439 * t377 + (m(6) * t426 + t306 - t431) * pkin(9)) * qJD(4), qJD(4) * t137 - t420, qJD(4) * t235 - t419; -qJD(2) * t6 + qJD(3) * t9, -qJD(3) * t1 - qJD(5) * t32 + qJD(6) * t41 - t424, t403, qJD(5) * t336 + qJD(6) * t335, t418, t417; qJD(2) * t26, qJD(3) * t17 + qJD(4) * t32 + t381 * t534 - t422, t420, -t418 - t534, 0, t325; qJD(2) * t37, -qJD(3) * t29 - qJD(4) * t41 - t381 * t535 - t421, t419, -t417 + t535, -t325, 0;];
Cq  = t14;
