% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRRR4_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:25:32
% EndTime: 2019-03-09 02:25:47
% DurationCPUTime: 8.77s
% Computational Cost: add. (17233->605), mult. (32076->856), div. (0->0), fcn. (31787->8), ass. (0->324)
t638 = -mrSges(6,1) / 0.2e1;
t637 = mrSges(6,2) / 0.2e1;
t394 = cos(pkin(10));
t396 = sin(qJ(5));
t400 = cos(qJ(4));
t500 = t396 * t400;
t393 = sin(pkin(10));
t399 = cos(qJ(5));
t509 = t393 * t399;
t318 = -t394 * t500 + t509;
t493 = t399 * t400;
t320 = t393 * t396 + t394 * t493;
t395 = sin(qJ(6));
t398 = cos(qJ(6));
t198 = t318 * t398 - t320 * t395;
t201 = t318 * t395 + t320 * t398;
t606 = m(7) * pkin(5);
t481 = t606 / 0.2e1;
t621 = t198 * mrSges(7,1) / 0.2e1 - t201 * mrSges(7,2) / 0.2e1;
t636 = t318 * t638 + t320 * t637 - (t198 * t398 + t201 * t395) * t481 - t621;
t397 = sin(qJ(4));
t307 = -t395 * t500 + t398 * t493;
t590 = t307 / 0.2e1;
t434 = t395 * t399 + t398 * t396;
t305 = t434 * t400;
t593 = -t305 / 0.2e1;
t432 = Ifges(7,5) * t590 + Ifges(7,6) * t593;
t634 = mrSges(7,2) / 0.2e1;
t635 = -mrSges(7,1) / 0.2e1;
t401 = -pkin(1) - pkin(2);
t486 = t394 * qJ(2) + t393 * t401;
t337 = -pkin(7) + t486;
t571 = pkin(8) * t400;
t574 = pkin(4) * t397;
t361 = t571 - t574;
t501 = t396 * t397;
t213 = t337 * t501 + t399 * t361;
t572 = pkin(5) * t397;
t184 = pkin(9) * t493 + t213 - t572;
t499 = t397 * t399;
t214 = -t337 * t499 + t396 * t361;
t203 = pkin(9) * t500 + t214;
t94 = t184 * t398 - t203 * t395;
t95 = t184 * t395 + t203 * t398;
t617 = t95 * t634 + t94 * t635 + t432;
t612 = -Ifges(7,3) * t397 / 0.2e1 - t617;
t451 = -t393 * qJ(2) + t394 * t401;
t449 = pkin(3) - t451;
t569 = t400 * pkin(4);
t570 = t397 * pkin(8);
t274 = t449 + t569 + t570;
t254 = t399 * t274;
t377 = pkin(9) * t499;
t133 = t254 + t377 + (-t337 * t396 + pkin(5)) * t400;
t471 = t337 * t493;
t144 = t471 + (pkin(9) * t397 + t274) * t396;
t525 = t144 * t395;
t81 = t133 * t398 - t525;
t182 = -t337 * t500 + t254;
t143 = t182 + t377;
t86 = t143 * t398 - t525;
t633 = -t81 + t86;
t317 = -t393 * t500 - t394 * t399;
t319 = t393 * t493 - t396 * t394;
t200 = t317 * t395 + t319 * t398;
t453 = t398 * t317 - t319 * t395;
t47 = -t200 * mrSges(7,1) - t453 * mrSges(7,2);
t632 = t47 * qJD(6);
t604 = -pkin(9) - pkin(8);
t359 = t604 * t396;
t360 = t604 * t399;
t245 = t359 * t395 - t360 * t398;
t452 = t398 * t359 + t360 * t395;
t327 = Ifges(7,6) * t434;
t495 = t398 * t399;
t433 = t395 * t396 - t495;
t328 = Ifges(7,5) * t433;
t488 = -t328 - t327;
t49 = -t245 * mrSges(7,1) - t452 * mrSges(7,2) + t488;
t631 = t49 * qJD(6);
t304 = t434 * t397;
t306 = -t395 * t501 + t397 * t495;
t620 = -t306 * mrSges(7,1) + t304 * mrSges(7,2);
t630 = qJD(6) * t620;
t382 = -pkin(5) * t399 - pkin(4);
t575 = m(7) * t382;
t524 = t144 * t398;
t82 = t133 * t395 + t524;
t85 = -t143 * t395 - t524;
t629 = t82 + t85;
t514 = t320 * t399;
t516 = t318 * t396;
t435 = t514 - t516;
t508 = t394 * t397;
t624 = (t516 / 0.2e1 - t514 / 0.2e1) * mrSges(6,3) - m(6) * (-pkin(4) * t508 + pkin(8) * t435) / 0.2e1 - m(7) * (t198 * t452 + t201 * t245 + t382 * t508) / 0.2e1;
t534 = t434 * mrSges(7,3);
t387 = Ifges(6,5) * t399;
t556 = Ifges(6,6) * t396;
t454 = t387 - t556;
t388 = Ifges(6,4) * t399;
t354 = Ifges(6,1) * t396 + t388;
t389 = t396 ^ 2;
t391 = t399 ^ 2;
t485 = t389 + t391;
t618 = -t213 * t396 + t214 * t399;
t529 = t399 * mrSges(6,2);
t531 = t396 * mrSges(6,1);
t351 = t529 + t531;
t557 = Ifges(6,2) * t396;
t444 = t557 - t388;
t350 = -mrSges(6,1) * t399 + mrSges(6,2) * t396;
t536 = t307 * mrSges(7,2);
t538 = t305 * mrSges(7,1);
t490 = -t538 / 0.2e1 - t536 / 0.2e1;
t255 = t393 * t304;
t510 = t393 * t397;
t256 = t433 * t510;
t430 = t255 * t635 + t256 * t634;
t615 = mrSges(6,3) * t485 - mrSges(5,2);
t537 = t306 * Ifges(7,4);
t176 = t304 * Ifges(7,2) + t400 * Ifges(7,6) - t537;
t283 = Ifges(7,4) * t304;
t177 = -t306 * Ifges(7,1) + Ifges(7,5) * t400 + t283;
t194 = Ifges(7,2) * t306 + t283;
t195 = Ifges(7,1) * t304 + t537;
t573 = pkin(5) * t396;
t456 = t337 - t573;
t261 = t456 * t397;
t539 = t304 * mrSges(7,3);
t246 = -mrSges(7,2) * t400 + t539;
t597 = t246 / 0.2e1;
t225 = mrSges(7,1) * t434 - mrSges(7,2) * t433;
t601 = t225 / 0.2e1;
t603 = t620 / 0.2e1;
t614 = t452 * t597 - (t177 / 0.4e1 + t194 / 0.4e1) * t433 - (-t195 / 0.4e1 + t176 / 0.4e1) * t434 + t261 * t601 + t382 * t603;
t591 = t306 / 0.2e1;
t594 = -t304 / 0.2e1;
t613 = t453 * t597 + (t200 * t591 + t453 * t594) * mrSges(7,3);
t611 = m(6) / 0.2e1;
t610 = m(6) / 0.4e1;
t609 = -m(7) / 0.4e1;
t608 = m(7) / 0.2e1;
t605 = mrSges(7,3) / 0.2e1;
t602 = t200 / 0.2e1;
t226 = mrSges(7,1) * t433 + mrSges(7,2) * t434;
t600 = -t226 / 0.2e1;
t533 = t434 * Ifges(7,4);
t228 = -Ifges(7,2) * t433 + t533;
t599 = t228 / 0.2e1;
t598 = -t452 / 0.2e1;
t248 = mrSges(7,1) * t400 + t306 * mrSges(7,3);
t596 = -t248 / 0.2e1;
t595 = t248 / 0.2e1;
t592 = -t306 / 0.2e1;
t589 = t317 / 0.2e1;
t325 = t397 * t354;
t588 = t325 / 0.4e1;
t587 = t434 / 0.2e1;
t345 = mrSges(6,1) * t400 + mrSges(6,3) * t499;
t586 = -t345 / 0.2e1;
t585 = t351 / 0.2e1;
t584 = -t396 / 0.2e1;
t583 = t396 / 0.2e1;
t582 = -t399 / 0.2e1;
t581 = t399 / 0.2e1;
t580 = -t400 / 0.2e1;
t579 = t400 / 0.2e1;
t578 = t400 / 0.4e1;
t576 = m(7) * t261;
t568 = t81 * mrSges(7,2);
t567 = t82 * mrSges(7,1);
t566 = t85 * mrSges(7,1);
t565 = t86 * mrSges(7,2);
t559 = Ifges(6,4) * t396;
t553 = pkin(5) * qJD(5);
t535 = t433 * mrSges(7,3);
t532 = t395 * t95;
t530 = t398 * t94;
t528 = t400 * mrSges(5,2);
t527 = t350 - mrSges(5,1);
t507 = t394 * t400;
t62 = (-t198 * t304 + t201 * t306) * t608 + 0.2e1 * (t507 * t609 + (t435 - t507) * t610) * t397;
t526 = qJD(1) * t62;
t465 = t510 / 0.2e1;
t410 = t200 * t596 + t465 * t620 + t613;
t15 = t410 - t621;
t523 = t15 * qJD(1);
t522 = t201 * t433;
t22 = t246 * t594 + t248 * t592 + t620 * t580 + (t306 ^ 2 / 0.2e1 + t304 ^ 2 / 0.2e1) * mrSges(7,3);
t519 = t22 * qJD(1);
t518 = t261 * t396;
t517 = t317 * t396;
t515 = t319 * t399;
t513 = t337 * t394;
t512 = t337 * t400;
t511 = t393 * t394;
t247 = mrSges(7,2) * t397 + mrSges(7,3) * t305;
t506 = t395 * t247;
t505 = t395 * t306;
t298 = Ifges(6,6) * t400 + t397 * t444;
t504 = t396 * t298;
t346 = -mrSges(6,1) * t397 + mrSges(6,3) * t493;
t503 = t396 * t346;
t352 = Ifges(6,2) * t399 + t559;
t502 = t396 * t352;
t498 = t397 * t400;
t249 = -mrSges(7,1) * t397 + mrSges(7,3) * t307;
t497 = t398 * t249;
t496 = t398 * t304;
t344 = mrSges(6,2) * t397 + mrSges(6,3) * t500;
t494 = t399 * t344;
t489 = Ifges(7,5) * t304 + Ifges(7,6) * t306;
t487 = Ifges(6,5) * t501 + Ifges(6,6) * t499;
t390 = t397 ^ 2;
t392 = t400 ^ 2;
t484 = t390 - t392;
t483 = qJD(4) * t397;
t482 = qJD(4) * t400;
t436 = -t515 + t517;
t41 = (t200 * t307 - t255 * t304 + t256 * t306 - t305 * t453) * t608 + m(6) * t436 * t580 + 0.2e1 * (m(7) * t484 / 0.4e1 + (-t390 * t485 + t484) * t610) * t393;
t480 = t41 * qJD(4);
t479 = mrSges(5,1) * t508;
t474 = -t81 / 0.2e1 + t86 / 0.2e1;
t473 = t82 / 0.2e1 + t85 / 0.2e1;
t472 = t226 + t527;
t470 = t387 / 0.2e1;
t467 = -t534 / 0.2e1;
t466 = t225 * t580 + t306 * t467 - t534 * t592;
t192 = -mrSges(7,1) * t304 - mrSges(7,2) * t306;
t463 = t192 * t583;
t343 = -mrSges(6,2) * t400 + mrSges(6,3) * t501;
t462 = t343 * t584;
t460 = t602 - t200 / 0.2e1;
t229 = -Ifges(7,1) * t433 - t533;
t459 = -t229 / 0.4e1 + t228 / 0.4e1;
t329 = Ifges(7,4) * t433;
t227 = -Ifges(7,2) * t434 - t329;
t230 = Ifges(7,1) * t434 - t329;
t458 = t230 / 0.4e1 + t227 / 0.4e1;
t322 = t397 * t351;
t457 = -t322 / 0.2e1 + t192 / 0.2e1;
t450 = t225 * t465;
t321 = -mrSges(6,1) * t499 + mrSges(6,2) * t501;
t448 = mrSges(6,3) * (t391 / 0.2e1 + t389 / 0.2e1);
t447 = t488 * t578;
t446 = -t397 * mrSges(5,1) - t528;
t355 = Ifges(6,1) * t399 - t559;
t445 = -Ifges(7,1) * t307 + Ifges(7,4) * t305;
t443 = -Ifges(7,4) * t307 + Ifges(7,2) * t305;
t183 = t274 * t396 + t471;
t262 = t456 * t400;
t193 = -t536 - t538;
t323 = t400 * t351;
t421 = -t323 / 0.2e1 + t193 / 0.2e1 + t345 * t583 + t343 * t582;
t13 = t246 * t590 + t247 * t591 + t248 * t593 + t249 * t594 - t421 * t400 + (t494 / 0.2e1 - t503 / 0.2e1 + t457) * t397 + (t261 * t397 - t262 * t400 - t304 * t94 - t305 * t81 + t306 * t95 + t307 * t82) * t608 + ((t183 * t400 + t214 * t397) * t399 + (-t182 * t400 - t213 * t397) * t396 + t484 * t337) * t611;
t376 = Ifges(6,4) * t501;
t299 = -Ifges(6,1) * t499 + Ifges(6,5) * t400 + t376;
t4 = -t95 * t246 - t94 * t248 + t322 * t512 - t82 * t247 - t81 * t249 + t443 * t594 - t183 * t344 - t182 * t346 + t445 * t591 - t262 * t192 - t261 * t193 - t214 * t343 - t449 * t446 - t213 * t345 + t177 * t590 + t176 * t593 - m(7) * (t261 * t262 + t81 * t94 + t82 * t95) - m(6) * (t182 * t213 + t183 * t214) + (t299 * t581 - t504 / 0.2e1 + (-Ifges(5,4) + t470 - t556 / 0.2e1) * t400 + t432) * t400 + (t337 * t323 - Ifges(7,5) * t306 + Ifges(7,6) * t304 + (-Ifges(5,1) + Ifges(5,2) + Ifges(6,3) + Ifges(7,3) - t391 * Ifges(6,1) / 0.2e1 - m(6) * t337 ^ 2 + (t388 - t557 / 0.2e1) * t396) * t400 + (Ifges(5,4) - t454) * t397) * t397;
t442 = -t4 * qJD(1) + t13 * qJD(3);
t16 = (-t629 * t304 + t633 * t306 + t493 * t572) * t608 + t321 * t580 + t397 * t462 + t499 * t586 + t390 * t448 + t22;
t324 = Ifges(6,2) * t499 + t376;
t408 = t261 * t620 + (t177 / 0.2e1 + t194 / 0.2e1) * t304 + (t82 * mrSges(7,3) - t195 / 0.2e1 + t176 / 0.2e1) * t306 - t81 * t539 + t489 * t579;
t5 = m(7) * (t81 * t85 + t82 * t86) + t86 * t246 + t85 * t248 + t487 * t579 + t182 * t343 - t183 * t345 + (t337 * t321 + (t299 / 0.2e1 + t324 / 0.2e1 - t182 * mrSges(6,3)) * t396 + (-t325 / 0.2e1 + t298 / 0.2e1 + t183 * mrSges(6,3) + (-t192 - t576) * pkin(5)) * t399) * t397 + t408;
t441 = t5 * qJD(1) + t16 * qJD(3);
t440 = t16 * qJD(1);
t12 = t81 * t246 - t82 * t248 + t408;
t439 = t12 * qJD(1) + t22 * qJD(3);
t302 = t390 * t513;
t17 = m(3) * qJ(2) + t198 * t248 + t201 * t246 + t318 * t345 + t320 * t343 + mrSges(3,3) + (t400 * mrSges(5,1) - t397 * mrSges(5,2) + mrSges(4,1)) * t393 + (mrSges(4,2) + (t192 - t322) * t397 + (-t390 - t392) * mrSges(5,3)) * t394 + m(7) * (t198 * t81 + t201 * t82 + t261 * t508) + m(6) * (t182 * t318 + t183 * t320 + t302) + m(5) * (t392 * t513 + t393 * t449 + t302) + m(4) * (-t393 * t451 + t394 * t486);
t438 = t17 * qJD(1) + t62 * qJD(3);
t437 = t213 * t317 + t214 * t319;
t431 = t213 * t638 + t214 * t637;
t429 = -t529 / 0.2e1 - t531 / 0.2e1;
t33 = t450 - t430;
t428 = (t255 * t398 + t256 * t395) * t606;
t362 = t393 ^ 2 * t498;
t59 = m(7) * (t200 * t256 + t255 * t453 + t362) + m(6) * (t436 * t510 + t362);
t409 = t453 * t249 / 0.2e1 + t247 * t602 + t255 * t595 + t256 * t597 + t346 * t589 + t319 * t344 / 0.2e1;
t420 = t200 * t95 + t255 * t81 + t256 * t82 + t453 * t94;
t423 = t182 * t396 - t183 * t399 + 0.2e1 * t512;
t9 = (t522 / 0.2e1 + t198 * t587) * mrSges(7,3) + (t528 + (-t350 / 0.2e1 + t600) * t397) * t394 + (t397 * t421 + t400 * t457) * t393 + (t423 * t510 + t437) * t611 + ((t261 * t400 + t262 * t397) * t393 + t420) * t608 + t479 + t409 + t624;
t427 = t9 * qJD(1) + t59 * qJD(2) + t41 * qJD(3);
t403 = ((t603 + t321 / 0.2e1) * t393 + (t515 / 0.2e1 - t517 / 0.2e1) * mrSges(6,3)) * t397 + (-pkin(5) * t390 * t509 + t629 * t453) * t608 + t343 * t589 + t319 * t586 + t613 + (t633 * t608 - t595) * t200;
t11 = t403 + t636;
t425 = t11 * qJD(1);
t80 = m(7) * (t304 * t305 + t306 * t307) + 0.4e1 * ((-0.1e1 + t485) * t610 + t609) * t498;
t424 = t13 * qJD(1) + t41 * qJD(2) + t80 * qJD(3);
t422 = t585 + t601 + t429;
t419 = t633 * t245 + t629 * t452;
t402 = t459 * t306 + t458 * t304 + (t324 / 0.4e1 + t299 / 0.4e1 + pkin(8) * t586) * t399 + (t245 * t591 - t433 * t474 - t434 * t473 + t452 * t594) * mrSges(7,3) - pkin(4) * t321 / 0.2e1 - t245 * t595 + t614;
t406 = t337 * t585 + (t354 / 0.4e1 - t444 / 0.4e1) * t396 + pkin(8) * t448 + (-t355 / 0.4e1 + t352 / 0.4e1 + (-t575 / 0.2e1 + t600) * pkin(5)) * t399;
t2 = t419 * t608 + (-t506 / 0.2e1 + t463 - t497 / 0.2e1 + 0.2e1 * (t518 / 0.4e1 - t532 / 0.4e1 - t530 / 0.4e1) * m(7)) * pkin(5) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t406) * t397 + (t588 - t298 / 0.4e1 - pkin(8) * t343 / 0.2e1) * t396 + t402 + (-0.3e1 / 0.4e1 * t556 + t387 / 0.4e1 - t328 / 0.4e1 - t327 / 0.4e1 + t470) * t400 + t431 + t617;
t405 = pkin(5) * t393 * t501 * t608 - t460 * t534;
t21 = t422 * t510 - t428 / 0.2e1 + t405 + t430;
t37 = t382 * t225 - (-t229 / 0.2e1 + t599) * t434 - (t230 / 0.2e1 + t227 / 0.2e1) * t433;
t27 = -t502 / 0.2e1 + t355 * t583 - pkin(4) * t351 + (-t444 / 0.2e1 + t354 / 0.2e1) * t399 + t37 + (t575 + t226) * t573;
t415 = t500 * t606;
t416 = (-t305 * t398 + t307 * t395) * t481 + t490;
t35 = t422 * t400 + t415 / 0.2e1 + t416;
t418 = t2 * qJD(1) + t21 * qJD(2) - t35 * qJD(3) + t27 * qJD(4);
t32 = t450 + t430;
t38 = t466 - t490;
t404 = (mrSges(7,3) * t598 + t458) * t304 + (t245 * t605 + t459) * t306 + t245 * t596 + t447 + t614;
t7 = t404 - t612;
t417 = t7 * qJD(1) + t32 * qJD(2) + t38 * qJD(3) + t37 * qJD(4);
t412 = (t398 * t597 + t395 * t596 + (t505 / 0.2e1 - t496 / 0.2e1) * mrSges(7,3)) * pkin(5);
t19 = -mrSges(7,1) * t473 + mrSges(7,2) * t474 + t412;
t349 = (mrSges(7,1) * t395 + mrSges(7,2) * t398) * pkin(5);
t46 = mrSges(7,1) * t460;
t57 = (t598 + t452 / 0.2e1) * mrSges(7,2);
t413 = -t19 * qJD(1) + t46 * qJD(2) - t57 * qJD(4) + t349 * qJD(5);
t358 = t390 * t511;
t336 = t349 * qJD(6);
t39 = t466 + t490;
t36 = t351 * t580 - t415 / 0.2e1 + t429 * t400 + t416 + t466;
t20 = t428 / 0.2e1 + t405 + t33 + 0.2e1 * t351 * t465;
t18 = -t567 / 0.2e1 - t568 / 0.2e1 + t566 / 0.2e1 - t565 / 0.2e1 + t412 + t489;
t14 = t410 + t621;
t10 = t403 - t636;
t8 = -t522 * t605 + t198 * t467 + t409 + ((t576 / 0.2e1 + t457) * t400 + (t262 * t608 + t423 * t611 + t421) * t397) * t393 - t394 * t446 / 0.2e1 - t479 / 0.2e1 - mrSges(5,2) * t507 / 0.2e1 + t437 * t611 + t420 * t608 + (t350 + t226) * t508 / 0.2e1 - t624;
t6 = t404 + t612;
t3 = qJD(2) * t62 + qJD(4) * t13 + qJD(5) * t16 + qJD(6) * t22;
t1 = (t506 + t497) * pkin(5) / 0.2e1 + t612 + (t406 - Ifges(6,3) / 0.2e1) * t397 + pkin(8) * t462 + pkin(5) * t463 + t402 + (pkin(5) * t518 + t419) * t608 + (t530 + t532) * t481 + Ifges(6,6) * t500 / 0.2e1 - Ifges(6,5) * t493 / 0.2e1 + t447 - t431 - t504 / 0.4e1 + t454 * t578 + t396 * t588;
t23 = [qJD(2) * t17 - qJD(4) * t4 + qJD(5) * t5 + qJD(6) * t12, 0.2e1 * ((t198 * t453 + t200 * t201 + t358) * t608 + (t317 * t318 + t319 * t320 + t358) * t611 + m(5) * (t358 + (t392 - 0.1e1) * t511) / 0.2e1) * qJD(2) + t8 * qJD(4) + t10 * qJD(5) + t14 * qJD(6) + t438, t3, t8 * qJD(2) + t1 * qJD(5) + t6 * qJD(6) + (-Ifges(5,5) + t355 * t584 + t444 * t581 + t502 / 0.2e1 + t354 * t582 + (-m(6) * pkin(4) + t527) * t337) * t482 + (mrSges(5,2) * t337 - Ifges(6,5) * t396 - Ifges(7,5) * t434 - Ifges(6,6) * t399 + Ifges(7,6) * t433 + Ifges(5,6)) * t483 + t442 + (-t94 * t534 - t95 * t535 - t433 * t443 / 0.2e1 + t445 * t587 + t382 * t193 + pkin(4) * t323 - t307 * t230 / 0.2e1 + t305 * t599 + t262 * t226 + t245 * t247 + t452 * t249 + m(7) * (t245 * t95 + t262 * t382 + t452 * t94) + (m(6) * t618 + t494 - t503) * pkin(8) + t618 * mrSges(6,3)) * qJD(4), t10 * qJD(2) + t1 * qJD(4) + (-t183 * mrSges(6,1) - t182 * mrSges(6,2) + t487 + t489 - t565 + t566) * qJD(5) + t18 * qJD(6) + (m(7) * (t395 * t86 + t398 * t85) + (-t496 + t505) * mrSges(7,3)) * t553 + t441, t14 * qJD(2) + t6 * qJD(4) + t18 * qJD(5) + (t489 - t567 - t568) * qJD(6) + t439; qJD(4) * t9 + qJD(5) * t11 + qJD(6) * t15 - t438, t59 * qJD(4), t480 - t526, t20 * qJD(5) + t33 * qJD(6) + t427 + (m(7) * (t245 * t256 + t255 * t452) - t256 * t535 - t255 * t534 + (-t615 * t397 + m(6) * (-t485 * t570 - t569) + (t472 + t575) * t400) * t393) * qJD(4), t20 * qJD(4) + (-t319 * mrSges(6,1) - t317 * mrSges(6,2) + (-t200 * t398 + t395 * t453) * t606 + t47) * qJD(5) + t632 + t425, t33 * qJD(4) + t47 * qJD(5) + t523 + t632; t3, t480 + t526, t80 * qJD(4), t36 * qJD(5) + t39 * qJD(6) + (t305 * t434 - t307 * t433) * qJD(4) * mrSges(7,3) + t615 * t482 + t472 * t483 + 0.2e1 * ((t485 * t571 - t574) * t611 + (t245 * t307 - t305 * t452 + t382 * t397) * t608) * qJD(4) + t424, t36 * qJD(4) + (t321 + t620 + (-t304 * t395 - t306 * t398) * t606) * qJD(5) + t630 + t440, t39 * qJD(4) + qJD(5) * t620 + t519 + t630; -qJD(2) * t9 + qJD(5) * t2 + qJD(6) * t7 - t442, qJD(5) * t21 + qJD(6) * t32 - t427, -qJD(5) * t35 + qJD(6) * t38 - t424, qJD(5) * t27 + qJD(6) * t37 (pkin(8) * t350 + t454 + t49) * qJD(5) + t631 + (m(7) * (-t245 * t398 + t395 * t452) + (-t395 * t434 + t398 * t433) * mrSges(7,3)) * t553 + t418, t49 * qJD(5) + t417 + t631; -qJD(2) * t11 - qJD(4) * t2 + qJD(6) * t19 - t441, -qJD(4) * t21 - qJD(6) * t46 - t425, t35 * qJD(4) - t440, qJD(6) * t57 - t418, -t336, -t336 - t413; -qJD(2) * t15 - qJD(4) * t7 - qJD(5) * t19 - t439, -t32 * qJD(4) + t46 * qJD(5) - t523, -t38 * qJD(4) - t519, -qJD(5) * t57 - t417, t413, 0;];
Cq  = t23;
