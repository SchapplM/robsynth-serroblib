% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PRRRPP2
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
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PRRRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPP2_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:48:53
% EndTime: 2019-03-08 22:49:12
% DurationCPUTime: 10.58s
% Computational Cost: add. (7595->667), mult. (18936->881), div. (0->0), fcn. (16856->8), ass. (0->323)
t395 = sin(qJ(4));
t580 = pkin(4) + pkin(5);
t506 = t395 * t580;
t398 = cos(qJ(4));
t510 = t398 * qJ(5);
t632 = t506 - t510;
t634 = m(7) * t632;
t421 = t634 / 0.2e1;
t334 = pkin(4) * t395 - t510;
t430 = m(6) * t334;
t636 = t421 + t430 / 0.2e1;
t635 = -t634 / 0.2e1;
t396 = sin(qJ(3));
t379 = t398 * mrSges(7,2);
t540 = t395 * mrSges(7,1);
t610 = t379 - t540;
t280 = t610 * t396;
t633 = -t280 / 0.2e1;
t615 = Ifges(6,4) + Ifges(5,5);
t378 = t396 * qJ(5);
t363 = t395 * t378;
t511 = t396 * t398;
t483 = t511 * t580;
t200 = -t363 - t483;
t594 = m(5) / 0.2e1;
t592 = m(6) / 0.2e1;
t590 = m(7) / 0.2e1;
t631 = m(6) + m(7);
t630 = mrSges(6,2) + mrSges(5,3);
t614 = Ifges(5,6) + Ifges(7,6);
t399 = cos(qJ(3));
t394 = sin(pkin(6));
t397 = sin(qJ(2));
t517 = t394 * t397;
t533 = cos(pkin(6));
t269 = t396 * t533 + t399 * t517;
t400 = cos(qJ(2));
t507 = t398 * t400;
t140 = t269 * t395 + t394 * t507;
t516 = t394 * t400;
t484 = t395 * t516;
t141 = t269 * t398 - t484;
t204 = -t398 * t517 + t399 * t484;
t205 = (t395 * t397 + t399 * t507) * t394;
t629 = t140 * t204 + t141 * t205;
t389 = t399 * pkin(4);
t328 = -pkin(3) * t399 - pkin(9) * t396 - pkin(2);
t513 = t395 * t399;
t499 = pkin(8) * t513 - t398 * t328;
t170 = t389 + t499;
t628 = t170 - t499;
t390 = t395 ^ 2;
t392 = t398 ^ 2;
t608 = t390 + t392;
t544 = Ifges(6,5) * t398;
t545 = Ifges(7,4) * t398;
t627 = t544 + t545 + (Ifges(7,2) + Ifges(6,3)) * t395;
t453 = t398 * mrSges(6,1) + t395 * mrSges(6,3);
t454 = t398 * mrSges(5,1) - t395 * mrSges(5,2);
t626 = -t454 - t453;
t509 = t398 * t141;
t515 = t395 * t140;
t625 = -t509 - t515;
t385 = Ifges(7,4) * t395;
t449 = Ifges(7,1) * t398 + t385;
t383 = Ifges(6,5) * t395;
t450 = Ifges(6,1) * t398 + t383;
t546 = Ifges(5,4) * t395;
t451 = Ifges(5,1) * t398 - t546;
t624 = t449 + t450 + t451;
t508 = t398 * t399;
t365 = mrSges(7,3) * t508;
t319 = -t396 * mrSges(7,1) - t365;
t623 = t319 / 0.2e1;
t622 = -t396 / 0.2e1;
t621 = t396 / 0.2e1;
t620 = -t399 / 0.2e1;
t619 = t399 / 0.2e1;
t617 = -mrSges(5,1) - mrSges(6,1);
t616 = mrSges(7,2) + mrSges(6,3);
t613 = pkin(9) - qJ(6);
t332 = mrSges(7,1) * t398 + mrSges(7,2) * t395;
t611 = -t453 - t332;
t609 = Ifges(6,6) * t395 + t398 * t615;
t387 = Ifges(5,4) * t398;
t448 = -Ifges(5,2) * t395 + t387;
t345 = Ifges(5,1) * t395 + t387;
t524 = t205 * t398;
t526 = t204 * t395;
t605 = t524 + t526;
t604 = Ifges(6,3) * t398 - t383;
t603 = Ifges(7,2) * t398 - t385;
t602 = mrSges(7,1) - t617;
t601 = Ifges(6,2) + Ifges(5,3) + Ifges(7,3);
t486 = Ifges(5,6) / 0.2e1 + Ifges(7,6) / 0.2e1;
t582 = -Ifges(6,6) / 0.2e1;
t461 = t582 + t486;
t600 = Ifges(6,6) * t621 - t396 * t461 + t448 * t620 + t614 * t622 + t619 * t627;
t537 = t399 * Ifges(7,5);
t233 = t396 * t449 + t537;
t235 = -t399 * Ifges(6,4) + t396 * t450;
t237 = -t399 * Ifges(5,5) + t396 * t451;
t487 = -Ifges(6,4) / 0.2e1 - Ifges(5,5) / 0.2e1;
t583 = -Ifges(7,5) / 0.2e1;
t462 = t583 - t487;
t599 = -t462 * t399 + t233 / 0.2e1 + t235 / 0.2e1 + t237 / 0.2e1;
t519 = t328 * t395;
t169 = t519 + (pkin(8) * t398 - qJ(5)) * t399;
t514 = t395 * t396;
t362 = qJ(6) * t514;
t127 = t362 + t169;
t554 = pkin(9) * t399;
t350 = pkin(3) * t396 - t554;
t214 = -pkin(8) * t511 + t395 * t350;
t173 = t214 + t378;
t130 = qJ(6) * t513 + t173;
t425 = -pkin(8) - t632;
t166 = t425 * t396;
t167 = t425 * t399;
t482 = -pkin(8) * t395 - pkin(4);
t518 = t350 * t398;
t180 = t396 * t482 - t518;
t208 = pkin(8) * t508 + t519;
t213 = pkin(8) * t514 + t518;
t434 = pkin(8) + t334;
t239 = t434 * t396;
t240 = t434 * t399;
t268 = t396 * t517 - t399 * t533;
t338 = t395 * mrSges(5,1) + t398 * mrSges(5,2);
t281 = t338 * t396;
t336 = t395 * mrSges(6,1) - t398 * mrSges(6,3);
t282 = t336 * t399;
t283 = t610 * t399;
t284 = t338 * t399;
t315 = -mrSges(5,2) * t396 - mrSges(5,3) * t513;
t320 = mrSges(5,1) * t396 - mrSges(5,3) * t508;
t313 = mrSges(5,2) * t399 - mrSges(5,3) * t514;
t380 = t399 * mrSges(6,3);
t492 = mrSges(6,2) * t514;
t323 = -t380 - t492;
t473 = -t313 / 0.2e1 - t323 / 0.2e1;
t364 = mrSges(7,3) * t514;
t381 = t399 * mrSges(7,2);
t312 = t364 - t381;
t572 = -t312 / 0.2e1;
t455 = t572 + t473;
t539 = t399 * mrSges(7,1);
t316 = -mrSges(7,3) * t511 + t539;
t317 = -mrSges(5,1) * t399 - mrSges(5,3) * t511;
t318 = mrSges(6,1) * t399 + mrSges(6,2) * t511;
t472 = t317 / 0.2e1 - t318 / 0.2e1;
t456 = -t316 / 0.2e1 + t472;
t314 = mrSges(7,2) * t396 + mrSges(7,3) * t513;
t322 = -mrSges(6,2) * t513 + mrSges(6,3) * t396;
t471 = -t322 / 0.2e1 - t314 / 0.2e1;
t279 = t336 * t396;
t573 = t279 / 0.2e1;
t474 = t573 + t633;
t491 = mrSges(6,2) * t508;
t321 = -t396 * mrSges(6,1) + t491;
t571 = t321 / 0.2e1;
t143 = qJ(6) * t511 - t499;
t98 = pkin(5) * t399 - t143 + t389;
t99 = (-qJ(6) * t399 - t350) * t398 + (-pkin(5) + t482) * t396;
t598 = (t281 / 0.2e1 + t474) * t269 + (t623 - t320 / 0.2e1 + t571) * t140 + (t315 / 0.2e1 - t471) * t141 + (pkin(8) * t269 * t396 - t140 * t213 + t141 * t214) * t594 + (t140 * t180 + t141 * t173 + t239 * t269) * t592 + (t130 * t141 + t140 * t99 - t166 * t269) * t590 + ((pkin(8) * t399 - t208 * t398 - t395 * t499) * t594 + (-t169 * t398 - t170 * t395 + t240) * t592 + (-t127 * t398 - t395 * t98 - t167) * t590 + t282 / 0.2e1 - t283 / 0.2e1 + t284 / 0.2e1 + t455 * t398 + t456 * t395) * t268;
t597 = 0.2e1 * m(7);
t596 = 2 * qJD(3);
t595 = 2 * qJD(4);
t593 = -m(6) / 0.2e1;
t591 = -m(7) / 0.2e1;
t589 = m(7) / 0.4e1;
t588 = mrSges(5,2) / 0.2e1;
t587 = -mrSges(7,2) / 0.2e1;
t586 = -mrSges(6,3) / 0.2e1;
t584 = mrSges(7,3) / 0.2e1;
t581 = t99 / 0.2e1;
t144 = t208 + t362;
t578 = -t144 / 0.2e1;
t577 = t180 / 0.2e1;
t276 = t453 * t396;
t574 = t276 / 0.2e1;
t570 = -t453 / 0.2e1;
t569 = t332 / 0.2e1;
t335 = t613 * t398;
t568 = t335 / 0.2e1;
t567 = t336 / 0.2e1;
t566 = t338 / 0.2e1;
t565 = -t395 / 0.2e1;
t564 = t395 / 0.2e1;
t562 = -t398 / 0.2e1;
t532 = qJ(5) * t395;
t445 = -pkin(4) * t398 - t532;
t327 = -pkin(3) + t445;
t559 = m(6) * t327;
t557 = m(7) * t395;
t556 = m(7) * t396;
t553 = -mrSges(6,2) + mrSges(7,3);
t550 = mrSges(4,2) * t396;
t538 = t399 * mrSges(4,2);
t536 = t143 + t98;
t534 = -t454 - mrSges(4,1);
t11 = (m(5) + t631) * (t269 + t625) * t268;
t531 = t11 * qJD(1);
t485 = t396 * t516;
t521 = t268 * t396;
t12 = m(7) * t629 + 0.4e1 * (t521 * t589 + m(4) * (t269 * t399 + t521) / 0.4e1 - m(4) * t517 / 0.4e1) * t516 + (m(5) + m(6)) * (t268 * t485 + t629);
t530 = t12 * qJD(1);
t528 = t140 * t398;
t527 = t141 * t395;
t525 = t205 * qJ(5);
t522 = t268 * t395;
t512 = t396 * t332;
t505 = -t127 + t144;
t504 = -t169 + t208;
t503 = t608 * pkin(9) * t268;
t502 = -t279 + t280;
t501 = t318 - t317;
t500 = t323 + t313;
t497 = qJD(3) * t395;
t494 = m(7) * t511;
t490 = t167 * t590;
t489 = t269 * t591;
t488 = -mrSges(6,1) / 0.2e1 - mrSges(7,1) / 0.2e1;
t477 = -t516 / 0.2e1;
t470 = -t332 / 0.2e1 + t570;
t467 = t396 * t477;
t465 = mrSges(5,1) / 0.2e1 - t488;
t464 = t588 + t587 + t586;
t463 = -mrSges(6,2) / 0.2e1 - mrSges(5,3) / 0.2e1 + t584;
t460 = 0.2e1 * (t589 + m(6) / 0.4e1) * t204;
t367 = Ifges(6,5) * t511;
t227 = -Ifges(6,6) * t399 + Ifges(6,3) * t514 + t367;
t368 = Ifges(7,4) * t511;
t229 = Ifges(7,2) * t514 + t399 * Ifges(7,6) + t368;
t231 = -t399 * Ifges(5,6) + t396 * t448;
t459 = t227 / 0.2e1 + t229 / 0.2e1 - t231 / 0.2e1;
t342 = Ifges(5,2) * t398 + t546;
t339 = t396 * mrSges(4,1) + t538;
t291 = t398 * t580 + pkin(3) + t532;
t330 = t613 * t395;
t426 = t605 * pkin(9);
t405 = -m(5) * (-pkin(3) * t485 + t426) / 0.2e1 + (t327 * t485 + t426) * t593 + (t204 * t330 + t205 * t335 - t291 * t485) * t591;
t431 = t398 * t463;
t432 = t395 * t463;
t3 = t205 * t431 + t204 * t432 + (-t339 / 0.2e1 + t538 / 0.2e1 + (mrSges(4,1) / 0.2e1 + t454 / 0.2e1 - t470) * t396) * t516 + t405 + t598;
t416 = Ifges(7,5) * t622 + t462 * t396 + t615 * t621 + t619 * t624;
t5 = -pkin(2) * t339 + t180 * t318 + t98 * t319 - t499 * t320 + t170 * t321 + t169 * t322 + t173 * t323 + t130 * t312 + t214 * t313 + t127 * t314 + t208 * t315 + t99 * t316 + t213 * t317 + t240 * t279 + t167 * t280 + t239 * t282 + t166 * t283 + m(5) * (t208 * t214 - t213 * t499) + m(6) * (t169 * t173 + t170 * t180 + t239 * t240) + m(7) * (t127 * t130 + t166 * t167 + t98 * t99) + (pkin(8) * t281 + Ifges(4,4) * t399 + t599 * t398 + (t399 * t461 + t459) * t395) * t399 + (-Ifges(4,4) * t396 + pkin(8) * t284 + t416 * t398 + t600 * t395 + (m(5) * pkin(8) ^ 2 + Ifges(4,1) - Ifges(4,2) - t601) * t399) * t396;
t444 = t3 * qJD(1) + t5 * qJD(2);
t275 = pkin(4) * t511 + t363;
t278 = t454 * t396;
t407 = (t278 / 0.2e1 + t574 + t512 / 0.2e1) * t268 + (t140 * t504 + t141 * t628 + t268 * t275) * t592 + (t140 * t505 + t141 * t536 - t200 * t268) * t590;
t414 = (-pkin(4) * t204 + t525) * t593 + (-t204 * t580 + t525) * t591;
t7 = t464 * t205 + t465 * t204 + (t396 * t431 - t456) * t141 + (t396 * t432 + t455) * t140 + t407 + t414;
t285 = t604 * t396;
t286 = t603 * t396;
t287 = t396 * t342;
t288 = -Ifges(7,1) * t514 + t368;
t289 = -Ifges(6,1) * t514 + t367;
t290 = t345 * t396;
t366 = Ifges(6,6) * t511;
t8 = t366 * t620 + t143 * t312 + t144 * t316 - t166 * t512 + t275 * t279 + t200 * t280 + t239 * t276 + t501 * t208 - t500 * t499 + m(6) * (-t169 * t499 + t170 * t208 + t239 * t275) + m(7) * (t127 * t143 + t144 * t98 + t166 * t200) + (pkin(8) * t278 + (t288 / 0.2e1 + t289 / 0.2e1 - t290 / 0.2e1 - t208 * mrSges(5,3) + t127 * mrSges(7,3) - t169 * mrSges(6,2) + t486 * t399 + t459) * t398 + (t285 / 0.2e1 + t286 / 0.2e1 + t287 / 0.2e1 + t98 * mrSges(7,3) - t499 * mrSges(5,3) - t170 * mrSges(6,2) - t599) * t395) * t396;
t443 = t7 * qJD(1) + t8 * qJD(2);
t24 = (-t312 - t323) * t399 + t502 * t511 + m(6) * (-t169 * t399 - t239 * t511) + m(7) * (-t127 * t399 + t166 * t511);
t413 = (t591 + t593) * (-t141 * t399 - t268 * t511);
t28 = t460 + t413;
t442 = -qJD(1) * t28 + qJD(2) * t24;
t33 = (t395 * t312 - t398 * t316 + m(7) * (t127 * t395 - t98 * t398)) * t396;
t40 = (t477 + t528 / 0.2e1 - t527 / 0.2e1) * t556;
t441 = -qJD(1) * t40 + qJD(2) * t33;
t410 = (-t239 * t395 + (-t327 * t396 - t554) * t398) * t593 + (t166 * t395 + t291 * t511 - t335 * t399) * t591;
t436 = m(6) * t577 + m(7) * t581;
t17 = t491 - t365 + t474 * t395 + (t398 * t470 + t488) * t396 + t410 + t436;
t55 = -t291 * t557 + (t559 + t611) * t395;
t440 = -qJD(2) * t17 - qJD(3) * t55;
t52 = t512 + (t483 / 0.4e1 + t363 / 0.4e1 - t200 / 0.4e1) * t597;
t83 = (-t510 / 0.4e1 + t506 / 0.4e1 + t632 / 0.4e1) * t597 - t610;
t439 = qJD(2) * t52 + qJD(3) * t83;
t424 = -0.2e1 * t399 * qJ(5) + t208;
t411 = -t380 - t381 + t424 * t592 + (t362 + t424) * t590;
t433 = m(7) * t578 + t208 * t593;
t32 = t411 + t433;
t353 = qJ(5) * t631 + t616;
t438 = qJD(2) * t32 + qJD(4) * t353;
t437 = -mrSges(6,2) * pkin(4) + mrSges(7,3) * t580 - Ifges(7,5);
t435 = qJ(5) * t553 - t614;
t427 = m(7) * (-t330 * t395 - t335 * t398);
t402 = (t239 * t334 + t275 * t327) * t592 + (-t166 * t632 + t200 * t291 + t330 * t505 + t335 * t536) * t590 - pkin(3) * t278 / 0.2e1 + t166 * t610 / 0.2e1 + t200 * t569 + t239 * t567 + t275 * t570 - t291 * t512 / 0.2e1 + t632 * t633 + t327 * t574 + t330 * t572 + t334 * t573 + t316 * t568 - t609 * t399 / 0.4e1;
t343 = Ifges(7,1) * t395 - t545;
t344 = Ifges(6,1) * t395 - t544;
t403 = pkin(8) * t566 + (t383 / 0.4e1 + t385 / 0.4e1 - t342 / 0.4e1 - t603 / 0.4e1 - t604 / 0.4e1 + mrSges(7,3) * t568 + (Ifges(5,1) / 0.4e1 + Ifges(6,1) / 0.4e1 + Ifges(7,1) / 0.4e1) * t398) * t398 + (-t345 / 0.4e1 - t344 / 0.4e1 - t343 / 0.4e1 - t387 / 0.4e1 + t330 * t584 + (Ifges(7,2) / 0.4e1 + Ifges(6,3) / 0.4e1 + Ifges(5,2) / 0.4e1) * t395 + (Ifges(7,4) / 0.4e1 + Ifges(6,5) / 0.4e1 - Ifges(5,4) / 0.4e1) * t398) * t395 + t630 * pkin(9) * (-t392 / 0.2e1 - t390 / 0.2e1);
t406 = (-pkin(4) * t180 + qJ(5) * t173) * t593 + (qJ(5) * t130 - t580 * t99) * t591 + pkin(4) * t571 + t130 * t587 + t173 * t586 + mrSges(6,1) * t577 - t213 * mrSges(5,1) / 0.2e1 + t214 * t588 + t580 * t623 + mrSges(7,1) * t581;
t408 = (t578 + t127 / 0.2e1) * mrSges(7,3) + (t208 / 0.2e1 - t169 / 0.2e1) * mrSges(6,2) + (t504 * t592 + t473) * pkin(9) + t227 / 0.4e1 + t229 / 0.4e1 - t231 / 0.4e1 + t288 / 0.4e1 + t289 / 0.4e1 - t290 / 0.4e1;
t409 = (-t143 / 0.2e1 - t98 / 0.2e1) * mrSges(7,3) + (-t499 / 0.2e1 + t170 / 0.2e1) * mrSges(6,2) + (t592 * t628 - t472) * pkin(9) + t233 / 0.4e1 + t235 / 0.4e1 + t237 / 0.4e1 - t285 / 0.4e1 - t286 / 0.4e1 - t287 / 0.4e1;
t1 = (-Ifges(6,2) / 0.2e1 - Ifges(5,3) / 0.2e1 - Ifges(7,3) / 0.2e1 + t403) * t396 + t406 + t402 + ((0.3e1 / 0.4e1 * Ifges(7,5) + t487) * t399 + t409) * t398 + ((0.3e1 / 0.4e1 * Ifges(5,6) + t582 + 0.3e1 / 0.4e1 * Ifges(7,6)) * t399 + t408) * t395 + t471 * qJ(5);
t417 = t566 + t567 - t610 / 0.2e1 + t636;
t10 = (-t464 * t398 - t465 * t395 + t635 - t430 / 0.2e1 + t417) * t268;
t13 = pkin(3) * t338 + t632 * t332 + t334 * t453 + t342 * t564 + t627 * t398 / 0.2e1 + (-t336 - t430) * t327 + (-t610 + t634) * t291 + (t345 + t343 + t344 + t448) * t562 + (-t603 - t604 + t624) * t565;
t420 = t10 * qJD(1) + t1 * qJD(2) - t13 * qJD(3);
t415 = m(7) * ((-t330 * t396 - t127) * t398 + (t335 * t396 - t98) * t395);
t26 = (t381 / 0.2e1 + t312 / 0.2e1) * t398 + (-t539 / 0.2e1 + t316 / 0.2e1) * t395 + t490 - t415 / 0.2e1;
t38 = (-t269 / 0.4e1 + t515 / 0.4e1 + t509 / 0.4e1) * t597;
t87 = mrSges(7,3) * t608 + t427;
t419 = -qJD(1) * t38 - qJD(2) * t26 + qJD(3) * t87;
t393 = t399 ^ 2;
t391 = t396 ^ 2;
t329 = t391 * pkin(8) * t516;
t292 = (qJD(2) * t511 + t497) * m(7);
t142 = m(7) * t335 + (m(6) * pkin(9) - t553) * t398;
t134 = t140 * qJ(5);
t107 = t421 + t635;
t51 = (0.2e1 * t593 - m(7)) * t522;
t49 = 0.2e1 * (t590 + t592) * t141;
t41 = (t527 - t528) * t556 / 0.2e1 + m(7) * t467;
t39 = t590 * t625 + t489;
t31 = t364 + t411 - t433 - t492;
t29 = t460 - t413;
t27 = t316 * t565 + t312 * t562 + t415 / 0.2e1 + t490 + (t379 / 0.2e1 - t540 / 0.2e1) * t399;
t19 = t279 * t565 + t280 * t564 + t488 * t396 - t410 + t436 + (t569 + t453 / 0.2e1) * t511;
t9 = t602 * t522 / 0.2e1 + (t417 + (t588 - t616 / 0.2e1) * t398 + t636) * t268;
t6 = -t456 * t141 + t455 * t140 + (t140 * t432 + t141 * t431) * t396 + t407 - t414 - t602 * t204 / 0.2e1 + (-mrSges(5,2) / 0.2e1 + t616 / 0.2e1) * t205;
t4 = -t405 - t605 * mrSges(7,3) / 0.2e1 + (t339 + t538) * t477 + (mrSges(4,1) + t332) * t467 + t626 * t485 / 0.2e1 + t598 + t630 * (t524 / 0.2e1 + t526 / 0.2e1);
t2 = t403 * t396 - t406 + t402 + (t537 / 0.4e1 + t409) * t398 + ((Ifges(5,6) / 0.4e1 + Ifges(7,6) / 0.4e1) * t399 + t408) * t395 + (t322 + t314) * qJ(5) / 0.2e1 + t601 * t621 + (Ifges(6,6) / 0.2e1 - t614 / 0.2e1) * t513 + (t583 + t615 / 0.2e1) * t508;
t14 = [qJD(2) * t12 + qJD(3) * t11, t4 * qJD(3) + t6 * qJD(4) + t29 * qJD(5) + t41 * qJD(6) + t530 + ((t312 + t500) * t205 + (t316 + t501) * t204 + ((-mrSges(4,1) * t399 - mrSges(3,1) + t550) * t397 + (-mrSges(3,2) + (t391 + t393) * mrSges(4,3) + (t281 - t502) * t396) * t400) * t394 + 0.2e1 * (t127 * t205 - t166 * t485 + t204 * t98) * t590 + 0.2e1 * (t204 * t499 + t205 * t208 + t329) * t594 + 0.2e1 * (t169 * t205 + t170 * t204 + t239 * t485) * t592 + m(4) * (t329 + (pkin(8) * t393 * t400 - pkin(2) * t397) * t394)) * qJD(2), t531 + t4 * qJD(2) + t9 * qJD(4) + t51 * qJD(5) + t39 * qJD(6) + (t291 * t489 + (-pkin(3) * t269 - t503) * t594 + (t269 * t327 - t503) * t592) * t596 + ((t534 + t611) * t269 + (mrSges(4,2) + t427 + t608 * (-mrSges(5,3) + t553)) * t268) * qJD(3), t6 * qJD(2) + t9 * qJD(3) + t49 * qJD(5) + ((-t141 * t580 - t134) * t590 + (-pkin(4) * t141 - t134) * t592) * t595 + (-t602 * t141 + (mrSges(5,2) - t616) * t140) * qJD(4), qJD(2) * t29 + qJD(3) * t51 + qJD(4) * t49, qJD(2) * t41 + qJD(3) * t39; qJD(3) * t3 + qJD(4) * t7 - qJD(5) * t28 - qJD(6) * t40 - t530, qJD(3) * t5 + qJD(4) * t8 + qJD(5) * t24 + qJD(6) * t33, t2 * qJD(4) + t19 * qJD(5) + t27 * qJD(6) + ((t130 * t335 + t167 * t291 + t330 * t99) * t590 + t240 * t559 / 0.2e1) * t596 + (t180 * mrSges(6,2) - t213 * mrSges(5,3) - t99 * mrSges(7,3) + t416) * t497 + t444 + (-Ifges(4,6) * t396 - pkin(3) * t284 + pkin(8) * t550 + t167 * t332 - t240 * t453 + t327 * t282 + t291 * t283 + t335 * t314 + t330 * t319 + (t173 * mrSges(6,2) + t214 * mrSges(5,3) - t130 * mrSges(7,3) - t600) * t398 + ((t315 + t322) * t398 + (-t320 + t321) * t395 + m(5) * (-t213 * t395 + t214 * t398) + m(6) * (t173 * t398 + t180 * t395)) * pkin(9) + (Ifges(4,5) + (t343 / 0.2e1 + t344 / 0.2e1 + t345 / 0.2e1) * t398 + (-t604 / 0.2e1 - t603 / 0.2e1 - t342 / 0.2e1) * t395 + (-m(5) * pkin(3) + t534) * pkin(8)) * t399) * qJD(3), t2 * qJD(3) + t31 * qJD(5) + ((-pkin(4) * t208 - qJ(5) * t499) * t592 + (qJ(5) * t143 - t144 * t580) * t590) * t595 + t443 + (-t144 * mrSges(7,1) + t143 * mrSges(7,2) + t366 + (t435 * t398 + (-t437 - t615) * t395) * t396 + t617 * t208 - (-mrSges(5,2) + mrSges(6,3)) * t499) * qJD(4), qJD(3) * t19 + qJD(4) * t31 + t442, qJD(3) * t27 + t441; -qJD(2) * t3 + qJD(4) * t10 - qJD(6) * t38 - t531, qJD(4) * t1 - qJD(5) * t17 - qJD(6) * t26 - t444, -qJD(4) * t13 - qJD(5) * t55 + qJD(6) * t87, t142 * qJD(5) + t107 * qJD(6) + t420 + (m(7) * (-qJ(5) * t330 - t335 * t580) - t335 * mrSges(7,1) - t330 * mrSges(7,2) + t437 * t398 + t435 * t395 + (m(6) * t445 + t626) * pkin(9) + t609) * qJD(4), qJD(4) * t142 + t440, qJD(4) * t107 + t419; -qJD(2) * t7 - qJD(3) * t10, -qJD(3) * t1 + qJD(5) * t32 + qJD(6) * t52 - t443, qJD(6) * t83 - t420, t353 * qJD(5), t438, t439; qJD(2) * t28, qJD(3) * t17 - qJD(4) * t32 - qJD(6) * t494 - t442, -qJD(6) * t557 - t440, -t438, 0, -t292; qJD(2) * t40 + qJD(3) * t38, qJD(3) * t26 - qJD(4) * t52 + qJD(5) * t494 - t441, -qJD(4) * t83 + qJD(5) * t557 - t419, -t439, t292, 0;];
Cq  = t14;
