% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP8_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP8_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP8_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:57
% EndTime: 2019-03-09 04:54:09
% DurationCPUTime: 5.70s
% Computational Cost: add. (6515->622), mult. (13297->770), div. (0->0), fcn. (10148->4), ass. (0->290)
t544 = m(7) / 0.2e1;
t378 = sin(qJ(4));
t365 = t378 * qJ(5);
t377 = pkin(4) + qJ(6);
t380 = cos(qJ(4));
t554 = t377 * t380 + t365;
t566 = t554 * t544;
t373 = t378 ^ 2;
t375 = t380 ^ 2;
t553 = t373 + t375;
t381 = cos(qJ(3));
t379 = sin(qJ(3));
t513 = pkin(8) * t379;
t319 = pkin(3) * t381 + t513;
t382 = -pkin(1) - pkin(7);
t476 = t381 * t382;
t169 = t378 * t319 + t380 * t476;
t138 = -qJ(5) * t381 - t169;
t335 = t378 * t476;
t168 = t319 * t380 - t335;
t139 = -pkin(4) * t381 - t168;
t482 = t378 * t379;
t282 = -mrSges(6,1) * t482 - mrSges(6,3) * t381;
t341 = mrSges(7,1) * t482;
t494 = t381 * mrSges(7,2);
t283 = t341 + t494;
t441 = t282 / 0.2e1 - t283 / 0.2e1;
t479 = t379 * t380;
t495 = t381 * mrSges(6,2);
t284 = -mrSges(6,1) * t479 + t495;
t531 = t284 / 0.2e1;
t493 = t381 * mrSges(7,3);
t281 = -mrSges(7,1) * t479 - t493;
t532 = t281 / 0.2e1;
t535 = -t139 / 0.2e1;
t540 = mrSges(7,3) / 0.2e1;
t541 = mrSges(6,3) / 0.2e1;
t542 = -mrSges(7,2) / 0.2e1;
t545 = -m(7) / 0.2e1;
t547 = -m(6) / 0.2e1;
t71 = t335 - t377 * t381 + (-pkin(5) * t379 - t319) * t380;
t92 = pkin(5) * t482 - t138;
t565 = (-pkin(4) * t139 - qJ(5) * t138) * t547 + (qJ(5) * t92 - t377 * t71) * t545 + pkin(4) * t531 + t138 * t541 + mrSges(6,2) * t535 - t168 * mrSges(5,1) / 0.2e1 + t169 * mrSges(5,2) / 0.2e1 + t377 * t532 + t71 * t540 + t92 * t542 + t441 * qJ(5);
t497 = t380 * mrSges(6,2);
t502 = t378 * mrSges(6,3);
t295 = t497 - t502;
t496 = t380 * mrSges(7,3);
t503 = t378 * mrSges(7,2);
t416 = t496 + t503;
t564 = t416 / 0.2e1 - t295 / 0.2e1;
t438 = pkin(4) * t380 + t365;
t293 = -pkin(3) - t438;
t300 = t378 * mrSges(5,1) + t380 * mrSges(5,2);
t366 = Ifges(7,6) * t378;
t301 = Ifges(7,3) * t380 + t366;
t505 = Ifges(7,6) * t380;
t302 = Ifges(7,3) * t378 - t505;
t501 = t378 * Ifges(6,6);
t303 = -t380 * Ifges(6,3) - t501;
t506 = Ifges(6,6) * t380;
t304 = Ifges(6,3) * t378 - t506;
t305 = Ifges(7,2) * t378 + t505;
t306 = -t380 * Ifges(7,2) + t366;
t307 = -Ifges(6,2) * t378 - t506;
t308 = -Ifges(6,2) * t380 + t501;
t507 = Ifges(5,4) * t378;
t309 = t380 * Ifges(5,2) + t507;
t371 = Ifges(5,4) * t380;
t310 = -Ifges(5,2) * t378 + t371;
t311 = Ifges(5,1) * t378 + t371;
t312 = Ifges(5,1) * t380 - t507;
t538 = pkin(5) + pkin(8);
t320 = t538 * t380;
t396 = -Ifges(5,2) / 0.4e1 - Ifges(7,2) / 0.4e1 - Ifges(6,3) / 0.4e1 + Ifges(6,2) / 0.4e1 + Ifges(7,3) / 0.4e1 + Ifges(5,1) / 0.4e1;
t451 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t318 = t538 * t378;
t527 = -t318 / 0.2e1;
t247 = -pkin(3) - t554;
t534 = t247 / 0.2e1;
t543 = -mrSges(6,2) / 0.2e1;
t561 = mrSges(6,1) + mrSges(5,3);
t563 = (mrSges(7,1) * t527 + t293 * t541 + mrSges(7,2) * t534 - t311 / 0.4e1 - t310 / 0.4e1 + t307 / 0.4e1 + t305 / 0.4e1 + t304 / 0.4e1 - t302 / 0.4e1 - t396 * t378 + (-Ifges(5,4) - t451) * t380) * t378 + (-t320 * mrSges(7,1) / 0.2e1 + t293 * t543 + mrSges(7,3) * t534 + t312 / 0.4e1 - t309 / 0.4e1 - t308 / 0.4e1 + t306 / 0.4e1 + t303 / 0.4e1 + t301 / 0.4e1 + t396 * t380) * t380 - t382 * t300 / 0.2e1 + t561 * pkin(8) * (-t375 / 0.2e1 - t373 / 0.2e1);
t539 = m(6) + m(7);
t562 = -mrSges(5,1) + mrSges(6,2);
t560 = Ifges(7,4) + Ifges(6,5);
t559 = Ifges(5,5) + Ifges(7,5);
t477 = t380 * t381;
t454 = mrSges(7,1) * t477;
t498 = t379 * mrSges(7,3);
t277 = t454 - t498;
t456 = mrSges(6,1) * t477;
t280 = mrSges(6,2) * t379 + t456;
t557 = t280 + t277;
t417 = t380 * mrSges(7,2) - t378 * mrSges(7,3);
t297 = t378 * pkin(4) - qJ(5) * t380;
t487 = qJ(6) * t378;
t268 = t297 + t487;
t518 = m(7) * t268;
t556 = -t417 + t518;
t418 = t378 * mrSges(6,2) + t380 * mrSges(6,3);
t520 = m(6) * t297;
t555 = -t418 + t520;
t552 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t376 = t381 ^ 2;
t550 = 2 * qJD(3);
t549 = 2 * qJD(4);
t548 = m(5) / 0.2e1;
t546 = m(6) / 0.2e1;
t537 = m(7) * t92;
t512 = pkin(8) * t381;
t515 = pkin(3) * t379;
t292 = qJ(2) - t512 + t515;
t270 = t380 * t292;
t360 = t379 * t382;
t157 = t360 * t378 - t270;
t457 = pkin(5) * t477;
t103 = -t157 - t457;
t536 = t103 / 0.2e1;
t481 = t378 * t381;
t249 = pkin(4) * t477 + qJ(5) * t481;
t533 = t249 / 0.2e1;
t530 = -t417 / 0.2e1;
t529 = -t416 / 0.2e1;
t528 = -t418 / 0.2e1;
t342 = Ifges(7,6) * t477;
t526 = t342 / 0.2e1;
t343 = Ifges(6,6) * t481;
t525 = -t343 / 0.2e1;
t524 = -t380 / 0.2e1;
t523 = -t381 / 0.2e1;
t521 = m(6) * t293;
t519 = m(7) * t247;
t517 = m(7) * t377;
t516 = m(7) * t378;
t511 = mrSges(6,1) + mrSges(7,1);
t510 = mrSges(6,3) + mrSges(7,2);
t509 = m(7) * qJD(5);
t508 = m(7) * qJD(6);
t471 = qJ(5) * t479 + t360;
t112 = -t377 * t482 + t471;
t339 = qJ(5) * t477;
t469 = -pkin(4) * t481 + t339;
t113 = (-t382 + t487) * t381 - t469;
t158 = t378 * t292 + t380 * t360;
t480 = t379 * qJ(5);
t115 = -t158 - t480;
t437 = t378 * t382 - pkin(4);
t117 = t379 * t437 - t270;
t171 = -pkin(4) * t482 + t471;
t172 = -t469 - t476;
t248 = t417 * t379;
t250 = t417 * t381;
t252 = t418 * t379;
t253 = t300 * t379;
t254 = t418 * t381;
t273 = -mrSges(5,2) * t381 + mrSges(5,3) * t482;
t455 = mrSges(5,3) * t481;
t499 = t379 * mrSges(5,2);
t274 = -t455 - t499;
t275 = mrSges(5,1) * t381 + mrSges(5,3) * t479;
t500 = t379 * mrSges(5,1);
t276 = -mrSges(5,3) * t477 + t500;
t278 = mrSges(6,1) * t481 - t379 * mrSges(6,3);
t279 = -mrSges(7,1) * t481 + t379 * mrSges(7,2);
t461 = Ifges(6,4) - t559;
t436 = t461 * t379;
t458 = Ifges(6,2) + Ifges(7,3) + Ifges(5,1);
t459 = Ifges(5,2) + Ifges(7,2) + Ifges(6,3);
t460 = -Ifges(5,6) + t560;
t427 = -qJ(6) + t437;
t428 = -t270 + t457;
t68 = t379 * t427 + t428;
t104 = -pkin(5) * t481 + t158;
t89 = t104 + t480;
t3 = -t89 * t283 - t71 * t277 - t68 * t281 - t168 * t276 - t169 * t274 + t157 * t275 - t158 * t273 - t172 * t252 + t171 * t254 - t139 * t280 - t138 * t278 - t117 * t284 - t115 * t282 - t113 * t248 + t112 * t250 - t92 * t279 - m(5) * (-t157 * t168 + t158 * t169) - m(6) * (t115 * t138 + t117 * t139 + t171 * t172) - m(7) * (t112 * t113 + t68 * t71 + t89 * t92) + (-qJ(2) * mrSges(4,1) + Ifges(4,4) * t381 - t382 * t253 - t460 * t481 + t461 * t477) * t381 + (-t300 * t476 + qJ(2) * mrSges(4,2) + (t525 - t436) * t380 + (m(5) * t382 ^ 2 + t458 * t375 + Ifges(4,1) - Ifges(4,2) - t552) * t381 + (t526 + (t459 * t378 + (-0.3e1 / 0.2e1 * Ifges(6,6) + 0.3e1 / 0.2e1 * Ifges(7,6) - 0.2e1 * Ifges(5,4)) * t380) * t381) * t378 + (t378 * t460 - Ifges(4,4)) * t379) * t379;
t504 = t3 * qJD(1);
t167 = qJ(6) * t477 + t249;
t419 = t380 * mrSges(5,1) - t378 * mrSges(5,2);
t251 = t419 * t381;
t4 = -t115 * t456 + t167 * t250 - t104 * t277 - t103 * t279 - m(7) * (t103 * t89 + t104 * t68 + t113 * t167) + (-m(6) * t172 + t254) * t249 + (-m(6) * t117 + t276 - t280) * t158 + (t172 * t295 - t113 * t416 + t382 * t251 + (t89 * mrSges(7,1) + t158 * mrSges(5,3) - t342 + (Ifges(6,6) + Ifges(5,4)) * t477 - t460 * t379) * t380 + (t117 * mrSges(6,1) + t68 * mrSges(7,1) - t343 + (Ifges(7,6) - Ifges(5,4)) * t481 - t436 + (t458 - t459) * t477) * t378) * t381 + (-m(6) * t115 + t274 - t278 + t455) * t157;
t492 = t4 * qJD(1);
t430 = t541 + mrSges(7,2) / 0.2e1 - mrSges(5,2) / 0.2e1;
t452 = t543 + t540;
t432 = mrSges(5,1) / 0.2e1 + t452;
t389 = t430 * t378 + t432 * t380 + t438 * t546 + t566;
t409 = -t157 * t380 + t158 * t378;
t411 = -t115 * t378 - t117 * t380;
t415 = t89 * t378 - t68 * t380;
t442 = t276 / 0.2e1 - t280 / 0.2e1;
t422 = -t277 / 0.2e1 + t442;
t443 = t274 / 0.2e1 - t278 / 0.2e1;
t423 = t279 / 0.2e1 + t443;
t7 = (t251 / 0.2e1 + (t503 / 0.2e1 + t496 / 0.2e1 - t497 / 0.2e1 + t502 / 0.2e1) * t381 + m(6) * t533 + t167 * t544) * t381 + (t422 * t380 + t423 * t378 + (t409 - t411) * t547 + (t103 * t380 + t104 * t378 - t415) * t545 + t553 * t381 * (mrSges(6,1) / 0.2e1 + mrSges(7,1) / 0.2e1 + mrSges(5,3) / 0.2e1)) * t379 + t389;
t491 = t7 * qJD(1);
t490 = t103 + t68;
t489 = t104 - t89;
t488 = -t419 - mrSges(4,1);
t472 = -t278 + t279;
t10 = t379 * mrSges(4,1) + t381 * mrSges(4,2) + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + (t276 - t557) * t380 + (t274 + t472) * t378 + m(7) * t415 + m(6) * t411 + m(5) * t409;
t486 = qJD(1) * t10;
t15 = (t250 + t254) * t477 + t472 * t379 + m(7) * (-t113 * t477 + t379 * t89) + m(6) * (-t115 * t379 - t172 * t477);
t485 = qJD(1) * t15;
t31 = m(7) * (t113 * t481 - t379 * t68) - t379 * t277 - t250 * t481;
t484 = qJD(1) * t31;
t478 = t379 * t381;
t462 = m(7) / 0.4e1 + m(6) / 0.4e1;
t374 = t379 ^ 2;
t468 = t374 + t376;
t425 = 0.2e1 * t462 * t468 * t380;
t51 = (t544 + t546) * t380 + t425;
t475 = t51 * qJD(1);
t474 = t115 + t158;
t473 = t117 - t157;
t470 = t553 * t512;
t467 = qJD(3) * t379;
t466 = qJD(3) * t381;
t465 = qJD(4) * t378;
t464 = qJD(4) * t380;
t160 = (-t374 / 0.2e1 - t376 / 0.2e1 - 0.1e1 / 0.2e1) * t516;
t463 = t160 * qJD(1);
t450 = t378 * t529;
t449 = -t482 / 0.2e1;
t448 = t482 / 0.2e1;
t447 = -t479 / 0.2e1;
t446 = t479 / 0.2e1;
t444 = t254 / 0.2e1 + t250 / 0.2e1;
t439 = -t416 + t519;
t429 = -Ifges(5,4) / 0.2e1 - t451;
t426 = t564 * t381;
t424 = -qJ(5) * t511 - Ifges(5,6);
t421 = -t302 / 0.2e1 + t307 / 0.2e1 - t311 / 0.2e1;
t420 = -t303 / 0.2e1 - t306 / 0.2e1 + t309 / 0.2e1;
t36 = 0.4e1 * (m(5) / 0.4e1 + t462) * (-0.1e1 + t553) * t478;
t401 = m(5) * (-t168 * t378 + t169 * t380);
t410 = -t138 * t380 + t139 * t378;
t383 = ((t273 / 0.2e1 - t441) * t380 + (t531 + t532 - t275 / 0.2e1) * t378 + (t378 * t71 + t380 * t92 + t113) * t544 + (t172 + t410) * t546 + t401 / 0.2e1 - t444) * t379 + (-t252 / 0.2e1 - t248 / 0.2e1 + t253 / 0.2e1 + (t499 / 0.2e1 + t423) * t380 + (t500 / 0.2e1 - t422) * t378 + (t378 * t68 + t380 * t89 - t112) * t544 + (-t115 * t380 + t117 * t378 - t171) * t546 + (t157 * t378 + t158 * t380 - 0.2e1 * t360) * t548) * t381;
t400 = m(7) * (-t318 * t380 + t320 * t378);
t6 = -t400 / 0.2e1 + t383;
t414 = t6 * qJD(1) + t36 * qJD(2);
t388 = t444 * t378 + (-t172 * t378 + (-t293 * t381 + t513) * t380) * t546 + (-t113 * t378 - t247 * t477 + t320 * t379) * t544;
t405 = m(6) * t535 + t545 * t71;
t11 = t452 * t381 + (t379 * t511 + t426) * t380 + t388 + t405;
t42 = (t295 + t439 + t521) * t378;
t413 = -qJD(1) * t11 + qJD(3) * t42;
t395 = (-t113 * t380 + t247 * t481 - t318 * t379) * t544 - t250 * t524;
t25 = -t341 + (t450 + t542) * t381 - t537 / 0.2e1 + t395;
t94 = t439 * t380;
t412 = -qJD(1) * t25 + qJD(3) * t94;
t364 = 0.2e1 * t480;
t392 = t510 * t379 + (t364 + t158) * t546 + (t104 + t364) * t544;
t404 = t104 * t545 + t158 * t547;
t29 = t392 + t404;
t322 = qJ(5) * t539 + t510;
t408 = qJD(1) * t29 + qJD(4) * t322;
t321 = mrSges(7,3) + t517;
t397 = (t377 - t427) * t544;
t35 = t498 + 0.2e1 * (-t457 / 0.4e1 + t270 / 0.4e1 - t103 / 0.4e1) * m(7) + t379 * t397;
t407 = qJD(1) * t35 + qJD(4) * t321;
t406 = -mrSges(6,1) * pkin(4) - mrSges(7,1) * t377 - Ifges(6,4);
t403 = m(6) * t438;
t367 = Ifges(7,5) * t380;
t368 = Ifges(6,5) * t378;
t369 = Ifges(5,5) * t380;
t370 = Ifges(7,4) * t378;
t402 = t369 / 0.4e1 + t368 / 0.4e1 + t370 / 0.4e1 + t367 / 0.4e1;
t399 = t339 * t544 + t469 * t546;
t384 = (t172 * t297 + t249 * t293) * t546 + (t113 * t268 + t167 * t247 + t318 * t489 + t320 * t490) * t544 - pkin(3) * t251 / 0.2e1 + t113 * t530 + t167 * t529 + t172 * t528 + t295 * t533 - t268 * t250 / 0.2e1 - t297 * t254 / 0.2e1 + t279 * t527 + t320 * t277 / 0.2e1;
t390 = t525 + (t536 + t68 / 0.2e1) * mrSges(7,1) + (t117 / 0.2e1 - t157 / 0.2e1) * mrSges(6,1) + (t473 * t546 - t442) * pkin(8);
t391 = t526 + (-t89 / 0.2e1 + t104 / 0.2e1) * mrSges(7,1) + (t158 / 0.2e1 + t115 / 0.2e1) * mrSges(6,1) + (t474 * t546 - t443) * pkin(8);
t1 = (-Ifges(6,1) / 0.2e1 - Ifges(7,1) / 0.2e1 - Ifges(5,3) / 0.2e1 + t563) * t381 + t384 + ((-Ifges(6,4) + 0.3e1 / 0.4e1 * Ifges(5,5) + 0.3e1 / 0.4e1 * Ifges(7,5)) * t380 + (0.3e1 / 0.4e1 * Ifges(7,4) + 0.3e1 / 0.4e1 * Ifges(6,5) - Ifges(5,6)) * t378 + t402) * t379 + t390 * t380 + t391 * t378 + t565;
t393 = t430 * t380 + (-t517 / 0.2e1 - t432) * t378;
t17 = (t520 / 0.2e1 + t518 / 0.2e1 + t300 / 0.2e1 + t528 + t530 + t393) * t381 + t399;
t9 = -pkin(3) * t300 - t268 * t416 + t297 * t295 + t555 * t293 + t556 * t247 + (-t305 / 0.2e1 - t304 / 0.2e1 + t310 / 0.2e1 - t421) * t380 + (t301 / 0.2e1 - t308 / 0.2e1 + t312 / 0.2e1 - t420) * t378;
t398 = t1 * qJD(1) - t17 * qJD(2) + t9 * qJD(3);
t313 = (qJD(1) * t379 + qJD(4)) * m(7);
t256 = t539 * t479;
t255 = t539 * t481;
t214 = -m(7) * t318 - t378 * mrSges(7,1);
t159 = (-t468 + 0.1e1) * t516 / 0.2e1;
t111 = m(7) * t320 + (m(6) * pkin(8) + t511) * t380;
t52 = t524 * t539 + t425;
t34 = t428 * t545 - t454 + m(7) * t536 + (t397 + mrSges(7,3)) * t379;
t27 = t381 * t450 + t537 / 0.2e1 + t494 / 0.2e1 + t395;
t24 = -t481 * t511 + t392 - t404;
t18 = t393 * t381 + t399 + (t300 + t555 + t556) * t523;
t12 = t495 / 0.2e1 - t493 / 0.2e1 + t380 * t426 + t388 - t405;
t8 = t389 + t276 * t447 + t278 * t448 + (-t167 * t381 + (t378 * t489 + t380 * t490) * t379) * t544 + (-t249 * t381 + (t378 * t474 + t380 * t473) * t379) * t546 + t251 * t523 + (t279 + t274) * t449 + t557 * t446 - t564 * t376 - t553 * t478 * (mrSges(7,1) + t561) / 0.2e1;
t5 = t400 / 0.2e1 + t383;
t2 = t402 * t379 + t384 + Ifges(6,4) * t446 + Ifges(5,6) * t448 + ((Ifges(7,5) / 0.4e1 + Ifges(5,5) / 0.4e1 - Ifges(6,4) / 0.2e1) * t379 + t390) * t380 + ((Ifges(7,4) / 0.4e1 + Ifges(6,5) / 0.4e1 - Ifges(5,6) / 0.2e1) * t379 + t391) * t378 + t560 * t449 + t559 * t447 + (t552 / 0.2e1 + t563) * t381 - t565;
t13 = [qJD(2) * t10 - qJD(3) * t3 - qJD(4) * t4 + qJD(5) * t15 + qJD(6) * t31, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t52 + qJD(6) * t159 + t486, -t504 + t5 * qJD(2) + t2 * qJD(4) + t12 * qJD(5) + t27 * qJD(6) + (t171 * t521 / 0.2e1 + (t112 * t247 + t318 * t71 + t320 * t92) * t544) * t550 + (-Ifges(4,5) + (t380 * t429 + t421) * t380 + (-t429 * t378 + (Ifges(6,3) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1 - Ifges(7,3) / 0.2e1 - Ifges(6,2) / 0.2e1) * t380 + t420) * t378 + (-m(5) * pkin(3) + t488) * t382) * t467 + (-t382 * mrSges(4,2) - t378 * t461 - t380 * t460 - Ifges(4,6)) * t466 + (pkin(3) * t253 - t112 * t416 + t171 * t295 + t247 * t248 + t293 * t252 + t318 * t281 + t320 * t283 + (-mrSges(6,1) * t138 + mrSges(7,1) * t92 + mrSges(5,3) * t169) * t380 + (t139 * mrSges(6,1) + t71 * mrSges(7,1) - t168 * mrSges(5,3)) * t378 + ((t273 - t282) * t380 + (-t275 + t284) * t378 + t401 + m(6) * t410) * pkin(8)) * qJD(3), -t492 + t8 * qJD(2) + t2 * qJD(3) + t24 * qJD(5) + t34 * qJD(6) + ((qJ(5) * t103 - t104 * t377) * t544 + (-pkin(4) * t158 - qJ(5) * t157) * t546) * t549 + (t103 * mrSges(7,2) - t104 * mrSges(7,3) + ((t424 + t560) * t380 + (-t406 - t559) * t378) * t381 + t562 * t158 + (mrSges(5,2) - mrSges(6,3)) * t157) * qJD(4), qJD(2) * t52 + qJD(3) * t12 + qJD(4) * t24 + t485, qJD(2) * t159 + qJD(3) * t27 + qJD(4) * t34 + t484; qJD(3) * t6 - qJD(4) * t7 + qJD(5) * t51 + qJD(6) * t160 - t486, t36 * qJD(3), t18 * qJD(4) + t255 * qJD(5) + (t295 - t416 + t488) * t467 + (t379 * t519 / 0.2e1 + (t470 - t515) * t548 + (t293 * t379 + t470) * t546) * t550 + (t380 * t508 + (m(7) * (t318 * t378 + t320 * t380) - mrSges(4,2) + t553 * (mrSges(5,3) + t511)) * qJD(3)) * t381 + t414, -t491 + t18 * qJD(3) + t256 * qJD(5) + ((-mrSges(7,3) + t562) * t464 + ((mrSges(5,2) - t510) * qJD(4) - t508) * t378 + (-t403 / 0.2e1 - t566) * t549) * t379, qJD(3) * t255 + qJD(4) * t256 + t475, t463 + (-t379 * t465 + t380 * t466) * m(7); -qJD(2) * t6 + qJD(4) * t1 + qJD(5) * t11 + qJD(6) * t25 + t504, -qJD(4) * t17 - t414, qJD(4) * t9 - qJD(5) * t42 - qJD(6) * t94, t111 * qJD(5) + t214 * qJD(6) + t406 * t464 + t424 * t465 + t398 + (m(7) * (-qJ(5) * t318 - t320 * t377) + t370 + t367 + t368 + t369 - t320 * mrSges(7,3) - t318 * mrSges(7,2) + (t295 - t403 - t419) * pkin(8)) * qJD(4), qJD(4) * t111 - t413, qJD(4) * t214 - t412; qJD(2) * t7 - qJD(3) * t1 + qJD(5) * t29 + qJD(6) * t35 + t492, qJD(3) * t17 + t491, -t398, qJD(5) * t322 + qJD(6) * t321, t408, t407; -qJD(2) * t51 - qJD(3) * t11 - qJD(4) * t29 - t379 * t508 - t485, -t475, t413, -t408 - t508, 0, -t313; -qJD(2) * t160 - qJD(3) * t25 - qJD(4) * t35 + t379 * t509 - t484, -t463, t412, -t407 + t509, t313, 0;];
Cq  = t13;
