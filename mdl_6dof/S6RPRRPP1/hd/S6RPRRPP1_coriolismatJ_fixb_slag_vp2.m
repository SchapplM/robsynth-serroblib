% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:28:00
% EndTime: 2019-03-09 04:28:17
% DurationCPUTime: 9.19s
% Computational Cost: add. (12876->658), mult. (27570->873), div. (0->0), fcn. (26679->8), ass. (0->329)
t379 = cos(qJ(4));
t377 = sin(qJ(4));
t501 = sin(pkin(10));
t437 = t501 * t377;
t502 = cos(pkin(10));
t319 = -t502 * t379 + t437;
t361 = -pkin(4) * t379 - pkin(3);
t439 = t502 * t377;
t397 = t379 * t501 + t439;
t186 = pkin(5) * t319 - qJ(6) * t397 + t361;
t210 = mrSges(7,1) * t319 - mrSges(7,3) * t397;
t622 = m(7) * t186 + t210;
t612 = Ifges(7,4) + Ifges(6,5);
t620 = Ifges(6,6) - Ifges(7,6);
t378 = sin(qJ(3));
t380 = cos(qJ(3));
t458 = -cos(pkin(9)) * pkin(1) - pkin(2);
t316 = -pkin(3) * t380 - t378 * pkin(8) + t458;
t359 = sin(pkin(9)) * pkin(1) + pkin(7);
t340 = t380 * t359;
t463 = t379 * t340;
t394 = t377 * (qJ(5) * t378 - t316) - t463;
t619 = t394 * t501;
t618 = t394 * t502;
t613 = -mrSges(6,3) - mrSges(7,2);
t524 = Ifges(7,5) * t319;
t217 = Ifges(7,1) * t397 + t524;
t312 = Ifges(6,4) * t319;
t219 = Ifges(6,1) * t397 - t312;
t617 = t219 + t217;
t339 = t378 * t359;
t535 = t378 * pkin(3);
t536 = pkin(8) * t380;
t343 = t535 - t536;
t236 = t377 * t339 + t379 * t343;
t488 = t378 * t379;
t237 = t377 * t343 - t359 * t488;
t411 = -t236 * t377 + t237 * t379;
t287 = t319 * t378;
t288 = t397 * t378;
t470 = pkin(4) * t488;
t121 = -pkin(5) * t287 + qJ(6) * t288 + t470;
t160 = mrSges(7,1) * t288 + mrSges(7,3) * t287;
t537 = pkin(4) * t377;
t193 = pkin(5) * t397 + qJ(6) * t319 + t537;
t208 = mrSges(7,1) * t397 + t319 * mrSges(7,3);
t209 = mrSges(6,1) * t397 - t319 * mrSges(6,2);
t368 = t379 * qJ(5);
t474 = t379 * pkin(8) + t368;
t533 = -qJ(5) - pkin(8);
t225 = -t533 * t439 + t474 * t501;
t490 = t377 * t378;
t297 = pkin(4) * t490 + t339;
t596 = t533 * t437 + t474 * t502;
t293 = t379 * t316;
t220 = -t340 * t377 + t293;
t462 = t378 * t368;
t179 = t220 - t462;
t66 = t179 * t501 - t618;
t67 = t179 * t502 + t619;
t440 = t225 * t66 + t596 * t67;
t232 = -mrSges(6,1) * t380 + t287 * mrSges(6,3);
t507 = t380 * mrSges(7,1);
t521 = t287 * mrSges(7,2);
t233 = t507 - t521;
t442 = t233 / 0.2e1 - t232 / 0.2e1;
t369 = t380 * mrSges(7,3);
t520 = t288 * mrSges(7,2);
t229 = -t369 - t520;
t506 = t380 * mrSges(6,2);
t230 = -t288 * mrSges(6,3) + t506;
t443 = t230 / 0.2e1 + t229 / 0.2e1;
t539 = -t380 / 0.4e1;
t155 = -t462 + t293 + (-t359 * t377 - pkin(4)) * t380;
t62 = t501 * t155 - t618;
t57 = -qJ(6) * t380 + t62;
t571 = t193 / 0.2e1;
t575 = t121 / 0.2e1;
t61 = t155 * t502 + t619;
t58 = t380 * pkin(5) - t61;
t586 = m(7) / 0.2e1;
t602 = -t319 * t612 - t397 * t620;
t89 = pkin(5) * t288 + qJ(6) * t287 + t297;
t616 = (t121 * t186 + t193 * t89 - t225 * t57 + t58 * t596 + t440) * t586 + t210 * t575 + t160 * t571 + t297 * t209 / 0.2e1 + t89 * t208 / 0.2e1 + t602 * t539 + t442 * t596 - t443 * t225;
t615 = mrSges(6,3) / 0.2e1;
t614 = m(6) + m(7);
t370 = Ifges(5,5) * t379;
t522 = Ifges(5,6) * t377;
t611 = Ifges(4,4) - t370 / 0.2e1 + t522 / 0.2e1;
t466 = t615 + mrSges(7,2) / 0.2e1;
t608 = t466 * t287;
t607 = t466 * t288;
t606 = t208 + t209;
t480 = t229 + t230;
t371 = Ifges(5,4) * t379;
t605 = -Ifges(5,2) * t377 + t371;
t337 = Ifges(5,1) * t377 + t371;
t374 = t377 ^ 2;
t375 = t379 ^ 2;
t473 = t374 + t375;
t528 = Ifges(5,4) * t377;
t335 = Ifges(5,2) * t379 + t528;
t314 = t378 * t335;
t330 = -mrSges(5,1) * t380 - mrSges(5,3) * t488;
t585 = -pkin(8) / 0.2e1;
t604 = t330 * t585 - t314 / 0.4e1;
t315 = t378 * t337;
t328 = t380 * mrSges(5,2) - mrSges(5,3) * t490;
t603 = t328 * t585 - t315 / 0.4e1;
t601 = t287 * t620 - t288 * t612;
t599 = -mrSges(5,1) * t379 + mrSges(5,2) * t377;
t124 = m(7) * t287;
t598 = qJD(6) * t124;
t597 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t456 = t501 * pkin(4);
t354 = t456 + qJ(6);
t325 = m(7) * t354 + mrSges(7,3);
t290 = t380 * t319;
t515 = t290 * mrSges(7,3);
t517 = t290 * mrSges(6,2);
t289 = t397 * t380;
t518 = t289 * mrSges(7,1);
t519 = t289 * mrSges(6,1);
t594 = t515 / 0.2e1 - t517 / 0.2e1 + t518 / 0.2e1 + t519 / 0.2e1;
t457 = t502 * pkin(4);
t360 = -t457 - pkin(5);
t584 = m(6) * pkin(4);
t593 = m(7) * t360 - t502 * t584 - mrSges(6,1) - mrSges(7,1);
t592 = t501 * t584 - mrSges(6,2) + t325;
t591 = 2 * qJD(3);
t590 = m(5) / 0.2e1;
t589 = -m(6) / 0.2e1;
t588 = m(6) / 0.2e1;
t587 = -m(7) / 0.2e1;
t583 = -mrSges(5,2) / 0.2e1;
t582 = t57 / 0.2e1;
t581 = t61 / 0.2e1;
t580 = t62 / 0.2e1;
t579 = -t66 / 0.2e1;
t578 = -t67 / 0.2e1;
t484 = t379 * t380;
t182 = t378 * pkin(4) - qJ(5) * t484 + t236;
t489 = t377 * t380;
t203 = -qJ(5) * t489 + t237;
t78 = t182 * t502 - t203 * t501;
t69 = -t378 * pkin(5) - t78;
t577 = -t69 / 0.2e1;
t576 = m(7) * t66;
t527 = Ifges(6,4) * t287;
t147 = -Ifges(6,2) * t288 - t380 * Ifges(6,6) - t527;
t574 = t147 / 0.4e1;
t161 = mrSges(6,1) * t288 - mrSges(6,2) * t287;
t573 = t161 / 0.2e1;
t525 = Ifges(7,5) * t288;
t164 = -Ifges(7,3) * t287 - t525;
t572 = t164 / 0.4e1;
t212 = Ifges(7,3) * t397 - t524;
t570 = t212 / 0.4e1;
t526 = Ifges(6,4) * t397;
t215 = -Ifges(6,2) * t319 + t526;
t569 = t215 / 0.4e1;
t509 = t378 * mrSges(7,1);
t516 = t290 * mrSges(7,2);
t235 = -t509 - t516;
t564 = t235 / 0.2e1;
t263 = t288 * mrSges(7,3);
t563 = t263 / 0.2e1;
t561 = t287 / 0.2e1;
t559 = -t289 / 0.2e1;
t558 = t289 / 0.2e1;
t557 = -t290 / 0.2e1;
t554 = -t319 / 0.2e1;
t552 = t319 / 0.2e1;
t551 = t397 / 0.2e1;
t549 = -t397 / 0.2e1;
t548 = -t335 / 0.4e1;
t338 = Ifges(5,1) * t379 - t528;
t547 = t338 / 0.4e1;
t546 = -t361 / 0.2e1;
t545 = t361 / 0.2e1;
t544 = -t377 / 0.2e1;
t543 = t377 / 0.2e1;
t542 = t378 / 0.2e1;
t541 = t379 / 0.2e1;
t540 = -t380 / 0.2e1;
t538 = t380 / 0.2e1;
t530 = mrSges(5,3) * t378;
t514 = t319 * mrSges(7,2);
t513 = t319 * mrSges(6,3);
t512 = t397 * mrSges(7,2);
t511 = t397 * mrSges(6,3);
t510 = t377 * mrSges(5,1);
t508 = t379 * mrSges(5,2);
t505 = t380 * Ifges(5,5);
t504 = t380 * Ifges(5,6);
t503 = -mrSges(4,1) + t599;
t24 = t380 * t229 - m(7) * (t89 * t287 - t380 * t57) - t287 * t160;
t497 = t24 * qJD(1);
t496 = t287 * t397;
t493 = t297 * t377;
t281 = t378 * t605 - t504;
t492 = t377 * t281;
t331 = t378 * mrSges(5,1) - mrSges(5,3) * t484;
t491 = t377 * t331;
t487 = t378 * t380;
t283 = t338 * t378 - t505;
t486 = t379 * t283;
t329 = -t378 * mrSges(5,2) - mrSges(5,3) * t489;
t485 = t379 * t329;
t334 = t508 + t510;
t313 = t380 * t334;
t79 = t501 * t182 + t502 * t203;
t452 = t288 * t552;
t453 = -t496 / 0.2e1;
t483 = mrSges(7,2) * t452 + mrSges(6,3) * t453;
t479 = t233 - t232;
t353 = pkin(4) * t489;
t298 = t340 + t353;
t472 = qJD(3) * t380;
t471 = t584 / 0.2e1;
t469 = t537 / 0.2e1;
t465 = -Ifges(7,4) / 0.2e1 - Ifges(6,5) / 0.2e1;
t464 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t461 = -t514 / 0.2e1;
t451 = t397 * t540;
t450 = t359 * t334 / 0.2e1;
t448 = -t489 / 0.2e1;
t447 = t484 / 0.2e1;
t446 = t319 * t538;
t266 = Ifges(7,5) * t287;
t145 = -t380 * Ifges(7,6) + Ifges(7,3) * t288 - t266;
t445 = t145 / 0.2e1 - t147 / 0.2e1;
t149 = -Ifges(7,1) * t287 - t380 * Ifges(7,4) + t525;
t269 = Ifges(6,4) * t288;
t151 = -Ifges(6,1) * t287 - t380 * Ifges(6,5) - t269;
t444 = t149 / 0.2e1 + t151 / 0.2e1;
t441 = t375 / 0.2e1 + t374 / 0.2e1;
t436 = -t287 * mrSges(7,1) + t263;
t262 = t287 * mrSges(6,1);
t435 = -t288 * mrSges(6,2) - t262;
t434 = t370 - t522;
t431 = t225 * t289 - t290 * t596;
t426 = t470 / 0.2e1;
t425 = mrSges(6,3) * t457;
t424 = mrSges(6,3) * t456;
t146 = -Ifges(7,5) * t290 + t378 * Ifges(7,6) + Ifges(7,3) * t289;
t148 = -Ifges(6,4) * t290 - Ifges(6,2) * t289 + t378 * Ifges(6,6);
t150 = -Ifges(7,1) * t290 + t378 * Ifges(7,4) + Ifges(7,5) * t289;
t152 = -Ifges(6,1) * t290 - Ifges(6,4) * t289 + t378 * Ifges(6,5);
t162 = t515 + t518;
t163 = -t517 + t519;
t221 = t377 * t316 + t463;
t228 = -mrSges(7,2) * t289 + mrSges(7,3) * t378;
t231 = -mrSges(6,2) * t378 - mrSges(6,3) * t289;
t234 = mrSges(6,1) * t378 + mrSges(6,3) * t290;
t282 = Ifges(5,6) * t378 + t380 * t605;
t284 = Ifges(5,5) * t378 + t338 * t380;
t68 = qJ(6) * t378 + t79;
t90 = pkin(5) * t289 + qJ(6) * t290 + t298;
t4 = (t458 * mrSges(4,2) - t492 / 0.2e1 + t486 / 0.2e1 - t465 * t290 - t464 * t289 + t611 * t380) * t380 + (t458 * mrSges(4,1) + t359 * t313 + t282 * t544 + t284 * t541 - t611 * t378 + t464 * t288 + t465 * t287 + (-Ifges(4,2) + Ifges(4,1) + (m(5) * t359 + t334) * t359 - t597) * t380) * t378 + m(7) * (t57 * t68 + t58 * t69 + t89 * t90) + m(6) * (t297 * t298 + t61 * t78 + t62 * t79) + t445 * t289 - t444 * t290 + (-t148 / 0.2e1 + t146 / 0.2e1) * t288 + (-t150 / 0.2e1 - t152 / 0.2e1) * t287 + t237 * t328 + t221 * t329 + t236 * t330 + t220 * t331 + t297 * t163 + t298 * t161 + t69 * t233 + t61 * t234 + t58 * t235 + t57 * t228 + t68 * t229 + t79 * t230 + t62 * t231 + t78 * t232 + t90 * t160 + t89 * t162 + m(5) * (t220 * t236 + t221 * t237);
t6 = -t443 * t290 + t442 * t289 + (t564 - t234 / 0.2e1) * t288 + (-t231 / 0.2e1 - t228 / 0.2e1) * t287 + (t328 * t541 + t330 * t544 - t313 / 0.2e1 - t162 / 0.2e1 - t163 / 0.2e1) * t380 + (t485 / 0.2e1 - t491 / 0.2e1 + t160 / 0.2e1 + t573 + (t510 / 0.2e1 + t508 / 0.2e1) * t378) * t378 + (-t79 * t287 - t78 * t288 - t61 * t289 - t62 * t290 + t297 * t378 - t298 * t380) * t588 + (-t68 * t287 + t69 * t288 + t58 * t289 - t57 * t290 + t89 * t378 - t380 * t90) * t586 + ((-t220 * t377 + t221 * t379 - t340) * t380 + (t411 + t339) * t378) * t590;
t419 = t4 * qJD(1) + t6 * qJD(2);
t165 = Ifges(6,2) * t287 - t269;
t166 = -Ifges(7,1) * t288 - t266;
t167 = -Ifges(6,1) * t288 + t527;
t5 = t220 * t328 - t221 * t330 - t297 * t262 + t89 * t263 + t121 * t160 + t480 * t67 + t479 * t66 + m(6) * (-t61 * t66 + t62 * t67) + m(7) * (t121 * t89 + t57 * t67 + t58 * t66) + (-t297 * mrSges(6,2) - t165 / 0.2e1 + t164 / 0.2e1 - t58 * mrSges(7,2) + t61 * mrSges(6,3) - t444) * t288 + (-t89 * mrSges(7,1) - t166 / 0.2e1 - t167 / 0.2e1 + t62 * mrSges(6,3) + t57 * mrSges(7,2) - t445) * t287 + ((t505 / 0.2e1 + t220 * mrSges(5,3) - t283 / 0.2e1 + t314 / 0.2e1 - mrSges(5,2) * t339) * t377 + (t504 / 0.2e1 - t221 * mrSges(5,3) - t281 / 0.2e1 - t315 / 0.2e1 + mrSges(5,1) * t339 + (m(6) * t297 + t161) * pkin(4)) * t379) * t378 + t601 * t540;
t417 = t61 * t287 - t62 * t288;
t7 = (t563 - t262 / 0.2e1) * t380 + (-t506 / 0.2e1 + t607 + t443) * t288 + (-t507 / 0.2e1 + t608 + t442) * t287 + (-t67 * t287 + t66 * t288 + t417) * t589 + (-t121 * t380 + (-t57 + t66) * t288 + (-t58 - t67) * t287) * t587 + (t328 * t543 + t330 * t541 + t441 * t530 + t447 * t584 - t538 * t599) * t378;
t418 = t5 * qJD(1) - t7 * qJD(2);
t30 = m(5) * (-0.1e1 + t473) * t487 + t614 * (t287 * t290 + t288 * t289 - t487);
t416 = t6 * qJD(1) + t30 * qJD(2);
t415 = t7 * qJD(1);
t14 = -t480 * t288 - t479 * t287 + m(7) * (-t58 * t287 - t57 * t288) + m(6) * t417;
t414 = qJD(1) * t14;
t390 = (t287 * t502 - t288 * t501) * t584;
t405 = m(7) * (-t360 * t287 - t354 * t288);
t387 = t405 / 0.2e1 + t390 / 0.2e1;
t402 = m(6) * t426 + m(7) * t575;
t403 = t435 + t436;
t26 = -t387 + t402 + t403;
t386 = (-t319 * t354 + t360 * t397) * t586 + (-t319 * t501 - t397 * t502) * t471;
t404 = m(6) * t469 + m(7) * t571;
t32 = -t386 + t404 + t606;
t413 = qJD(1) * t26 + qJD(3) * t32;
t412 = -t225 * t287 - t288 * t596;
t188 = m(7) * t397;
t410 = qJD(1) * t124 - qJD(3) * t188;
t408 = m(7) * t577 + t509 / 0.2e1;
t383 = -t353 * t589 - t193 * t380 * t587 + t313 / 0.2e1 + t606 * t538;
t384 = (t289 * t360 - t290 * t354) * t586 + (-t289 * t502 - t290 * t501) * t471 + mrSges(5,1) * t448 + t484 * t583 - t594;
t11 = t383 + t384;
t211 = mrSges(6,1) * t319 + mrSges(6,2) * t397;
t309 = Ifges(7,5) * t397;
t213 = Ifges(7,3) * t319 + t309;
t214 = -Ifges(6,2) * t397 - t312;
t216 = -Ifges(7,1) * t319 + t309;
t218 = -Ifges(6,1) * t319 - t526;
t381 = (t354 * t68 + t360 * t69) * t586 + Ifges(6,6) * t559 + Ifges(7,6) * t558 + t236 * mrSges(5,1) / 0.2e1 + t237 * t583 + t354 * t228 / 0.2e1 + t360 * t564 + t68 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t577 + t78 * mrSges(6,1) / 0.2e1 - t79 * mrSges(6,2) / 0.2e1 + (t501 * t79 + t502 * t78) * t471 + Ifges(5,5) * t447 + Ifges(5,6) * t448 + t234 * t457 / 0.2e1 + t231 * t456 / 0.2e1 + t612 * t557 + t597 * t542;
t396 = -t225 * t62 - t596 * t61 + t440;
t2 = -(t165 + t151 + t149) * t319 / 0.4e1 - (t218 + t216 + t213) * t287 / 0.4e1 + t603 * t377 + t604 * t379 - t381 - t613 * (t596 * t561 - t225 * t288 / 0.2e1 + t67 * t554 + t66 * t551) + (t167 + t166 + t145) * t397 / 0.4e1 - t397 * t574 - (t337 + t605) * t490 / 0.4e1 + t58 * t461 + t161 * t469 + t378 * t450 + t616 + t186 * t436 / 0.2e1 + t599 * t535 / 0.2e1 + t211 * t426 + t486 / 0.4e1 + t473 * t530 * t585 - (t214 + t617) * t288 / 0.4e1 - t492 / 0.4e1 + t434 * t539 + t435 * t545 + t488 * t547 + t488 * t548 + t287 * t569 + t288 * t570 + t319 * t572 - t511 * t580 + t513 * t581 - t512 * t582 + ((t361 * t488 + t493) * pkin(4) + t396) * t588;
t8 = -pkin(3) * t334 + t186 * t208 + (t337 / 0.2e1 + t605 / 0.2e1) * t379 + (t338 / 0.2e1 - t335 / 0.2e1 + pkin(4) * t211) * t377 - (t215 / 0.2e1 - t213 / 0.2e1 - t216 / 0.2e1 - t218 / 0.2e1) * t397 + (-t214 / 0.2e1 - t217 / 0.2e1 - t219 / 0.2e1 + t212 / 0.2e1) * t319 + (m(6) * t537 + t209) * t361 + t622 * t193;
t401 = t2 * qJD(1) - t11 * qJD(2) + t8 * qJD(3);
t25 = -(t319 ^ 2 + t397 ^ 2) * t613 + t614 * (t225 * t397 - t319 * t596);
t389 = (t588 + t586) * (t319 * t287 + t288 * t397);
t42 = (t589 + t587) * t378 + t389;
t385 = (-t319 * t62 - t397 * t61 + t412) * t589 + (-t319 * t57 + t397 * t58 + t412) * t587;
t391 = t298 * t588 + t90 * t586 + t594;
t9 = -(t442 - t608) * t397 + (t443 - t607) * t319 + t385 + t391;
t400 = -qJD(1) * t9 + qJD(2) * t42 + qJD(3) * t25;
t157 = (-t451 + t559) * m(7);
t388 = (t186 * t287 - t380 * t596 - t397 * t89) * t586 + t210 * t561 + t160 * t549;
t21 = (t446 + t290 / 0.2e1) * mrSges(7,2) + t388 + t408;
t54 = t622 * t397;
t399 = qJD(1) * t21 + qJD(2) * t157 - qJD(3) * t54;
t392 = -t369 + ((-qJ(6) - t354) * t380 + t62) * t586;
t28 = -t576 / 0.2e1 + t392;
t398 = qJD(1) * t28 + qJD(4) * t325;
t156 = (-t451 + t558) * m(7);
t56 = 0.2e1 * t586 * t596 - t514;
t43 = t386 + t404;
t41 = t542 * t614 + t389;
t36 = t387 + t402;
t27 = -t520 + t576 / 0.2e1 + t392;
t20 = mrSges(7,2) * t446 - t516 / 0.2e1 + t388 - t408;
t12 = t288 * t461 + t496 * t615 - t383 + t384 + t483;
t10 = mrSges(7,2) * t453 + mrSges(6,3) * t452 + t232 * t549 + t233 * t551 + t480 * t554 - t385 + t391 + t483;
t3 = t6 * qJD(3) - t7 * qJD(4);
t1 = t381 + (t504 / 0.4e1 - t281 / 0.4e1 + pkin(4) * t573 + t603) * t377 + (-t165 / 0.4e1 - t151 / 0.4e1 - t149 / 0.4e1 + t572 + (t578 + t581) * mrSges(6,3) + (t578 - t58 / 0.2e1) * mrSges(7,2)) * t319 - (t574 - t167 / 0.4e1 - t166 / 0.4e1 - t145 / 0.4e1 + (t580 + t579) * mrSges(6,3) + (t579 + t582) * mrSges(7,2)) * t397 + (pkin(4) * t493 + t396) * t588 + (t283 / 0.4e1 + t604) * t379 + (t450 - t441 * pkin(8) * mrSges(5,3) + (pkin(3) * mrSges(5,2) / 0.2e1 - t337 / 0.4e1 - t605 / 0.4e1) * t377 + (-pkin(3) * mrSges(5,1) / 0.2e1 + t547 + t548 + (t211 / 0.2e1 + m(6) * t545) * pkin(4)) * t379) * t378 + (-t218 / 0.4e1 - t216 / 0.4e1 - t213 / 0.4e1 - t186 * mrSges(7,1) / 0.2e1 + t569 + t466 * t596) * t287 + (mrSges(6,2) * t546 - t219 / 0.4e1 - t217 / 0.4e1 - t214 / 0.4e1 + t570 - t466 * t225) * t288 + t370 * t539 + t186 * t563 + t262 * t546 + t616;
t13 = [qJD(3) * t4 + qJD(4) * t5 + qJD(5) * t14 - qJD(6) * t24, t3, t1 * qJD(4) + t10 * qJD(5) + t20 * qJD(6) + (t335 * t544 + t337 * t541 + t359 * t503 + Ifges(4,5)) * t472 + ((t186 * t90 + t225 * t69 + t596 * t68) * t586 + (-t225 * t78 + t298 * t361 + t596 * t79) * t588 + (-pkin(3) * t340 + pkin(8) * t411) * t590) * t591 + t419 + (mrSges(4,2) * t339 - t78 * t511 - t79 * t513 - t68 * t514 + t69 * t512 - Ifges(4,6) * t378 + t361 * t163 - pkin(3) * t313 + t298 * t211 + t90 * t210 + t186 * t162 + t282 * t541 + t284 * t543 + t146 * t552 + t148 * t554 + t213 * t558 + t215 * t559 + (Ifges(5,5) * t377 + Ifges(5,6) * t379 - t319 * t620 + t612 * t397) * t542 + t617 * t557 + (t152 + t150) * t551 + (t228 + t231) * t596 + (-t234 + t235) * t225 + (t485 - t491) * pkin(8) + t411 * mrSges(5,3)) * qJD(3), t1 * qJD(3) + (-t221 * mrSges(5,1) - t220 * mrSges(5,2) - Ifges(5,5) * t490 - Ifges(5,6) * t488 + t287 * t424 + t288 * t425 + t354 * t521 - t360 * t520 + t592 * t67 + t593 * t66 + t601) * qJD(4) + t36 * qJD(5) + t27 * qJD(6) + t418, qJD(3) * t10 + qJD(4) * t36 + t414, t20 * qJD(3) + t27 * qJD(4) - t497; t3, qJD(3) * t30, t12 * qJD(4) + t41 * qJD(5) + t156 * qJD(6) + (mrSges(5,3) * t473 - mrSges(4,2)) * t472 + ((t361 * t378 + t431) * t588 + (t186 * t378 + t431) * t586 + (t473 * t536 - t535) * t590) * t591 + t416 + ((t210 + t211 + t503) * t378 - t613 * (t289 * t397 + t290 * t319)) * qJD(3), t12 * qJD(3) + (-mrSges(5,1) * t488 + mrSges(5,2) * t490 + t390 - t403 + t405) * qJD(4) - t598 - t415, qJD(3) * t41, t156 * qJD(3) - qJD(4) * t124; qJD(4) * t2 - qJD(5) * t9 + qJD(6) * t21 - t419, -qJD(4) * t11 + qJD(5) * t42 + qJD(6) * t157 - t416, qJD(4) * t8 + qJD(5) * t25 - qJD(6) * t54 (pkin(8) * t599 - t225 * t592 + t319 * t425 - t354 * t512 - t360 * t514 - t397 * t424 + t593 * t596 + t434 + t602) * qJD(4) + t43 * qJD(5) + t56 * qJD(6) + t401, qJD(4) * t43 + t400, qJD(4) * t56 + t399; -qJD(3) * t2 - qJD(5) * t26 + qJD(6) * t28 - t418, qJD(3) * t11 + t415, -qJD(5) * t32 - t401, t325 * qJD(6), -t413, t398; qJD(3) * t9 + qJD(4) * t26 - t414 + t598, -qJD(3) * t42, qJD(4) * t32 - qJD(6) * t188 - t400, t413, 0, t410; -t21 * qJD(3) - t28 * qJD(4) - t124 * qJD(5) + t497, -t157 * qJD(3), qJD(5) * t188 - t399, -t398, -t410, 0;];
Cq  = t13;
