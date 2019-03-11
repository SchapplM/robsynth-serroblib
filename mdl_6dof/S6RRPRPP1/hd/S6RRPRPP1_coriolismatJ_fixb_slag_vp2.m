% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RRPRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:15
% EndTime: 2019-03-09 09:45:36
% DurationCPUTime: 11.25s
% Computational Cost: add. (22210->653), mult. (44518->867), div. (0->0), fcn. (50117->8), ass. (0->322)
t378 = sin(qJ(4));
t379 = cos(qJ(4));
t510 = sin(pkin(10));
t512 = cos(pkin(10));
t342 = t510 * t378 - t512 * t379;
t344 = -t378 * t512 - t379 * t510;
t511 = sin(pkin(9));
t513 = cos(pkin(9));
t545 = sin(qJ(2));
t546 = cos(qJ(2));
t347 = -t511 * t546 - t513 * t545;
t250 = t342 * t347;
t345 = t511 * t545 - t513 * t546;
t181 = mrSges(6,1) * t345 - t250 * mrSges(6,3);
t528 = t250 * mrSges(7,2);
t182 = -mrSges(7,1) * t345 + t528;
t252 = t344 * t347;
t178 = -mrSges(6,2) * t345 - t252 * mrSges(6,3);
t333 = t345 * mrSges(7,3);
t525 = t252 * mrSges(7,2);
t183 = t333 - t525;
t486 = t178 + t183;
t557 = -t344 / 0.2e1;
t560 = -t342 / 0.2e1;
t625 = t344 / 0.2e1;
t424 = t181 * t625 + t182 * t557 + t486 * t560;
t373 = -pkin(2) * t546 - pkin(1);
t274 = t345 * pkin(3) + t347 * pkin(8) + t373;
t471 = t545 * pkin(7);
t353 = -qJ(3) * t545 - t471;
t473 = t546 * pkin(7);
t355 = qJ(3) * t546 + t473;
t615 = t511 * t353 + t513 * t355;
t500 = t615 * t379;
t394 = t378 * (qJ(5) * t347 + t274) + t500;
t628 = t512 * t394;
t167 = t379 * t274 - t378 * t615;
t374 = t379 * qJ(5);
t125 = t347 * t374 + t167;
t93 = t345 * pkin(4) + t125;
t60 = t510 * t93 + t628;
t47 = qJ(6) * t345 + t60;
t629 = t394 * t510;
t59 = t512 * t93 - t629;
t48 = -t345 * pkin(5) - t59;
t493 = t347 * t378;
t278 = -mrSges(5,2) * t345 + mrSges(5,3) * t493;
t487 = t379 * t278;
t492 = t347 * t379;
t280 = t345 * mrSges(5,1) + mrSges(5,3) * t492;
t489 = t378 * t280;
t592 = m(7) / 0.2e1;
t594 = m(6) / 0.2e1;
t558 = t342 / 0.2e1;
t618 = mrSges(6,3) + mrSges(7,2);
t636 = t618 * (t250 * t557 + t252 * t558);
t66 = t125 * t510 + t628;
t67 = t125 * t512 - t629;
t638 = ((t59 - t67) * t344 + (-t60 + t66) * t342) * t594 + ((-t48 - t67) * t344 + (-t47 + t66) * t342) * t592 - t489 / 0.2e1 + t487 / 0.2e1 + t424 - t636;
t460 = t511 * pkin(2);
t370 = t460 + pkin(8);
t434 = (-qJ(5) - t370) * t378;
t360 = t379 * t370;
t478 = t360 + t374;
t264 = -t512 * t434 + t478 * t510;
t605 = t510 * t434 + t478 * t512;
t412 = t250 * t264 - t252 * t605;
t593 = -m(7) / 0.2e1;
t595 = -m(6) / 0.2e1;
t637 = t636 - (-t342 * t60 + t344 * t59 + t412) * t595 - (-t342 * t47 - t344 * t48 + t412) * t593;
t462 = t513 * pkin(2);
t372 = -t462 - pkin(3);
t352 = -t379 * pkin(4) + t372;
t249 = t342 * pkin(5) + t344 * qJ(6) + t352;
t283 = mrSges(7,1) * t342 + mrSges(7,3) * t344;
t635 = m(7) * t249 + t283;
t634 = -t264 / 0.2e1;
t472 = t545 * pkin(2);
t633 = m(4) * t472;
t632 = m(6) * t352;
t617 = Ifges(7,4) + Ifges(6,5);
t630 = Ifges(6,6) - Ifges(7,6);
t626 = -t347 * mrSges(4,1) - t345 * mrSges(4,2);
t624 = t605 / 0.2e1;
t251 = t344 * t345;
t177 = mrSges(6,2) * t347 - t251 * mrSges(6,3);
t184 = -t251 * mrSges(7,2) - mrSges(7,3) * t347;
t623 = t177 + t184;
t375 = Ifges(5,5) * t379;
t530 = Ifges(5,6) * t378;
t616 = Ifges(4,4) - t375 / 0.2e1 + t530 / 0.2e1;
t376 = Ifges(5,4) * t379;
t358 = Ifges(5,1) * t378 + t376;
t613 = -t342 * t617 + t344 * t630;
t612 = -t250 * t630 - t252 * t617;
t416 = Ifges(5,2) * t378 - t376;
t271 = t347 * t358;
t549 = -t370 / 0.2e1;
t611 = t278 * t549 + t271 / 0.4e1;
t538 = Ifges(5,4) * t378;
t356 = Ifges(5,2) * t379 + t538;
t270 = t347 * t356;
t610 = t280 * t549 + t270 / 0.4e1;
t608 = t378 ^ 2 + t379 ^ 2;
t273 = m(7) * t344;
t607 = qJD(6) * t273;
t606 = Ifges(7,2) + Ifges(5,3) + Ifges(6,3);
t459 = t510 * pkin(4);
t361 = t459 + qJ(6);
t351 = m(7) * t361 + mrSges(7,3);
t533 = Ifges(7,5) * t342;
t285 = -Ifges(7,3) * t344 - t533;
t340 = Ifges(6,4) * t342;
t287 = Ifges(6,2) * t344 - t340;
t290 = -Ifges(7,1) * t344 + t533;
t292 = -Ifges(6,1) * t344 - t340;
t603 = -t292 / 0.4e1 - t290 / 0.4e1 - t287 / 0.4e1 + t285 / 0.4e1;
t337 = Ifges(7,5) * t344;
t286 = Ifges(7,3) * t342 - t337;
t536 = Ifges(6,4) * t344;
t288 = -Ifges(6,2) * t342 - t536;
t289 = -Ifges(7,1) * t342 - t337;
t291 = -Ifges(6,1) * t342 + t536;
t602 = -t291 / 0.4e1 - t289 / 0.4e1 + t288 / 0.4e1 - t286 / 0.4e1;
t253 = t345 * t342;
t522 = t253 * mrSges(7,3);
t524 = t253 * mrSges(6,2);
t526 = t251 * mrSges(7,1);
t527 = t251 * mrSges(6,1);
t601 = -t522 / 0.2e1 + t524 / 0.2e1 + t526 / 0.2e1 + t527 / 0.2e1;
t461 = t512 * pkin(4);
t371 = -t461 - pkin(5);
t590 = m(6) * pkin(4);
t600 = m(7) * t371 - t512 * t590 - mrSges(6,1) - mrSges(7,1);
t599 = t510 * t590 - mrSges(6,2) + t351;
t596 = m(5) / 0.2e1;
t591 = m(4) * pkin(2);
t589 = t47 / 0.2e1;
t588 = -t48 / 0.2e1;
t276 = -t347 * pkin(3) + t345 * pkin(8) + t472;
t302 = -t513 * t353 + t355 * t511;
t170 = t378 * t276 - t302 * t379;
t495 = t345 * t378;
t126 = qJ(5) * t495 + t170;
t169 = t379 * t276 + t302 * t378;
t94 = -t347 * pkin(4) + t345 * t374 + t169;
t61 = -t126 * t510 + t512 * t94;
t50 = t347 * pkin(5) - t61;
t587 = -t50 / 0.2e1;
t586 = t59 / 0.2e1;
t585 = t60 / 0.2e1;
t584 = -t66 / 0.2e1;
t583 = -t67 / 0.2e1;
t537 = Ifges(6,4) * t250;
t107 = -Ifges(6,2) * t252 + Ifges(6,6) * t345 + t537;
t582 = t107 / 0.4e1;
t534 = Ifges(7,5) * t252;
t133 = Ifges(7,3) * t250 - t534;
t581 = t133 / 0.4e1;
t516 = t347 * mrSges(7,1);
t523 = t253 * mrSges(7,2);
t180 = t516 + t523;
t580 = t180 / 0.2e1;
t579 = -t181 / 0.2e1;
t578 = -t251 / 0.2e1;
t577 = t251 / 0.2e1;
t576 = t253 / 0.2e1;
t544 = pkin(4) * t378;
t275 = -pkin(5) * t344 + qJ(6) * t342 + t544;
t570 = t275 / 0.2e1;
t569 = t283 / 0.2e1;
t555 = -t345 / 0.2e1;
t553 = t345 / 0.4e1;
t552 = -t347 / 0.2e1;
t551 = t356 / 0.4e1;
t359 = Ifges(5,1) * t379 - t538;
t550 = -t359 / 0.4e1;
t548 = t378 / 0.2e1;
t547 = t379 / 0.2e1;
t541 = mrSges(4,3) * t345;
t540 = mrSges(4,3) * t347;
t535 = Ifges(5,5) * t345;
t531 = Ifges(5,6) * t345;
t179 = -mrSges(6,1) * t347 - mrSges(6,3) * t253;
t417 = -mrSges(5,1) * t379 + t378 * mrSges(5,2);
t267 = t347 * t417;
t284 = mrSges(6,1) * t342 - mrSges(6,2) * t344;
t410 = t251 * t264 + t253 * t605;
t381 = (-t345 * t370 * t608 - t372 * t347) * t596 + (-t347 * t352 + t410) * t594 + (-t249 * t347 + t410) * t592 - t267 / 0.2e1 + (-t345 * t511 + t347 * t513) * t591 / 0.2e1 + (t283 + t284) * t552 + t608 * mrSges(5,3) * t555 + t618 * (t251 * t557 + t253 * t560);
t277 = mrSges(5,2) * t347 + mrSges(5,3) * t495;
t494 = t345 * t379;
t279 = -t347 * mrSges(5,1) + mrSges(5,3) * t494;
t62 = t512 * t126 + t510 * t94;
t49 = -qJ(6) * t347 + t62;
t383 = (t169 * t379 + t378 * t170) * t596 + (-t342 * t61 - t344 * t62) * t594 + (t342 * t50 - t344 * t49) * t592 + t277 * t548 + t279 * t547 + t633 / 0.2e1;
t7 = t179 * t560 + t180 * t558 + t557 * t623 - t381 + t383 + t626;
t529 = qJD(1) * t7;
t104 = Ifges(7,5) * t253 - Ifges(7,6) * t347 + Ifges(7,3) * t251;
t106 = Ifges(6,4) * t253 - Ifges(6,2) * t251 - Ifges(6,6) * t347;
t108 = Ifges(7,1) * t253 - Ifges(7,4) * t347 + Ifges(7,5) * t251;
t110 = Ifges(6,1) * t253 - Ifges(6,4) * t251 - Ifges(6,5) * t347;
t129 = -t522 + t526;
t130 = t524 + t527;
t131 = mrSges(7,1) * t252 - mrSges(7,3) * t250;
t132 = mrSges(6,1) * t252 + mrSges(6,2) * t250;
t168 = t378 * t274 + t500;
t201 = -Ifges(5,6) * t347 + t345 * t416;
t203 = -Ifges(5,5) * t347 - t345 * t359;
t215 = -pkin(4) * t495 + t615;
t216 = -pkin(4) * t493 + t302;
t354 = t378 * mrSges(5,1) + t379 * mrSges(5,2);
t268 = t354 * t345;
t269 = t354 * t347;
t109 = Ifges(7,1) * t250 + Ifges(7,4) * t345 + t534;
t244 = Ifges(6,4) * t252;
t111 = Ifges(6,1) * t250 + Ifges(6,5) * t345 - t244;
t442 = t109 / 0.2e1 + t111 / 0.2e1;
t241 = Ifges(7,5) * t250;
t105 = Ifges(7,6) * t345 + Ifges(7,3) * t252 + t241;
t443 = t105 / 0.2e1 - t107 / 0.2e1;
t467 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1;
t468 = Ifges(6,5) / 0.2e1 + Ifges(7,4) / 0.2e1;
t204 = -t347 * t359 + t535;
t488 = t379 * t204;
t202 = t347 * t416 + t531;
t490 = t378 * t202;
t75 = t251 * pkin(5) - t253 * qJ(6) + t215;
t76 = t252 * pkin(5) - t250 * qJ(6) + t216;
t3 = (-t379 * t203 / 0.2e1 + t201 * t548 - mrSges(4,2) * t472 - t616 * t347 + t467 * t252 - t468 * t250 + (Ifges(4,1) - Ifges(4,2) - t606) * t345) * t347 - t615 * t269 + (Ifges(3,1) - Ifges(3,2)) * t546 * t545 - (-t108 / 0.2e1 - t110 / 0.2e1) * t250 + m(5) * (t167 * t169 + t168 * t170 + t302 * t615) + (-t545 ^ 2 + t546 ^ 2) * Ifges(3,4) + (-t106 / 0.2e1 + t104 / 0.2e1) * t252 - pkin(1) * (mrSges(3,1) * t545 + mrSges(3,2) * t546) + t442 * t253 + t443 * t251 + m(6) * (t215 * t216 + t59 * t61 + t60 * t62) + m(7) * (t47 * t49 + t48 * t50 + t75 * t76) + (-t488 / 0.2e1 + t490 / 0.2e1 + mrSges(4,1) * t472 + t468 * t253 - t467 * t251 + t616 * t345) * t345 + (t626 + t633) * t373 + t76 * t129 + t75 * t131 + t60 * t177 + t62 * t178 + t59 * t179 + t48 * t180 + t61 * t181 + t50 * t182 + t49 * t183 + t47 * t184 + t215 * t132 + t216 * t130 + t168 * t277 + t170 * t278 + t167 * t279 + t169 * t280 - t302 * t268;
t521 = t3 * qJD(1);
t520 = t342 * mrSges(7,2);
t519 = t342 * mrSges(6,3);
t518 = t344 * mrSges(7,2);
t517 = t344 * mrSges(6,3);
t127 = t250 * mrSges(7,1) + t252 * mrSges(7,3);
t128 = t250 * mrSges(6,1) - t252 * mrSges(6,2);
t134 = -Ifges(6,2) * t250 - t244;
t135 = -Ifges(7,1) * t252 + t241;
t136 = -Ifges(6,1) * t252 - t537;
t485 = t181 - t182;
t475 = pkin(4) * t492;
t95 = t250 * pkin(5) + t252 * qJ(6) - t475;
t4 = t76 * t127 + t95 * t131 + t216 * t128 + t167 * t278 - t168 * t280 + t302 * t267 + t486 * t67 - t485 * t66 + m(7) * (t47 * t67 + t48 * t66 + t76 * t95) + m(6) * (-t59 * t66 + t60 * t67) + (t59 * mrSges(6,3) - t48 * mrSges(7,2) + t133 / 0.2e1 - t134 / 0.2e1 - t442) * t252 - (t60 * mrSges(6,3) + t47 * mrSges(7,2) - t135 / 0.2e1 - t136 / 0.2e1 - t443) * t250 + ((t270 / 0.2e1 + t204 / 0.2e1 + t535 / 0.2e1 - t167 * mrSges(5,3)) * t378 + (-t271 / 0.2e1 + t202 / 0.2e1 + t531 / 0.2e1 + t168 * mrSges(5,3) + (-m(6) * t216 - t132) * pkin(4)) * t379) * t347 + t612 * t345 / 0.2e1;
t515 = t4 * qJD(1);
t501 = t302 * t347;
t9 = t486 * t253 - t485 * t251 + (-t487 + t489 + t541) * t345 + (-t131 - t132 + t269 + t540) * t347 + m(6) * (-t216 * t347 - t251 * t59 + t253 * t60) + m(7) * (t251 * t48 + t253 * t47 - t347 * t76) + m(5) * (-t501 + (t167 * t378 - t168 * t379) * t345) + m(4) * (-t345 * t615 - t501);
t514 = t9 * qJD(1);
t13 = -t486 * t252 - t485 * t250 + m(6) * (-t250 * t59 - t252 * t60) + m(7) * (t250 * t48 - t252 * t47);
t509 = qJD(1) * t13;
t21 = t345 * t183 + m(7) * (-t250 * t76 + t345 * t47) - t250 * t131;
t508 = qJD(1) * t21;
t505 = t216 * t378;
t491 = t370 * t378;
t450 = t344 * t555;
t146 = (t450 + t578) * m(7);
t477 = t146 * qJD(1);
t476 = t590 / 0.2e1;
t474 = t66 * t592;
t470 = t544 / 0.2e1;
t469 = mrSges(6,3) / 0.2e1 + mrSges(7,2) / 0.2e1;
t465 = mrSges(5,3) * t370 / 0.2e1;
t449 = t342 * t555;
t448 = t495 / 0.2e1;
t447 = -t494 / 0.2e1;
t441 = t178 / 0.2e1 + t183 / 0.2e1;
t440 = t579 + t182 / 0.2e1;
t439 = t264 * t66 + t605 * t67;
t281 = -t344 * mrSges(7,1) + t342 * mrSges(7,3);
t282 = -t344 * mrSges(6,1) - t342 * mrSges(6,2);
t435 = t375 - t530;
t428 = -t475 / 0.2e1;
t427 = mrSges(6,3) * t461;
t426 = mrSges(6,3) * t459;
t420 = 0.2e1 * (-m(6) / 0.4e1 - m(7) / 0.4e1) * t347;
t389 = (-t249 * t250 + t344 * t76 + t345 * t605) * t592 - t250 * t569 + t131 * t625;
t408 = m(7) * t587 - t516 / 0.2e1;
t19 = (t449 - t253 / 0.2e1) * mrSges(7,2) + t389 + t408;
t96 = t635 * t344;
t415 = qJD(1) * t19 + qJD(2) * t96;
t388 = (t250 * t371 - t252 * t361) * t592 + (-t250 * t512 - t252 * t510) * t476;
t402 = m(6) * t428 + t95 * t592;
t25 = t127 + t128 - t388 + t402;
t391 = (-t342 * t510 + t344 * t512) * t590;
t405 = m(7) * (-t361 * t342 - t371 * t344);
t387 = t405 / 0.2e1 + t391 / 0.2e1;
t403 = t281 + t282;
t404 = m(6) * t470 + m(7) * t570;
t63 = -t387 + t403 + t404;
t414 = qJD(1) * t25 + qJD(2) * t63;
t390 = (t593 + t595) * (t250 * t342 + t252 * t344);
t35 = t420 + t390;
t413 = qJD(1) * t35;
t124 = m(7) * t250;
t409 = -qJD(1) * t124 + qJD(2) * t273;
t380 = (t361 * t49 + t371 * t50) * t592 + Ifges(6,6) * t578 + Ifges(7,6) * t577 + t169 * mrSges(5,1) / 0.2e1 - t170 * mrSges(5,2) / 0.2e1 + t361 * t184 / 0.2e1 + t371 * t580 + t49 * mrSges(7,3) / 0.2e1 + mrSges(7,1) * t587 + t61 * mrSges(6,1) / 0.2e1 - t62 * mrSges(6,2) / 0.2e1 + (t510 * t62 + t512 * t61) * t476 + Ifges(5,5) * t447 + Ifges(5,6) * t448 + t179 * t461 / 0.2e1 + t177 * t459 / 0.2e1 + t617 * t576 + t606 * t552;
t382 = (t249 * t95 - t264 * t47 + t275 * t76 + t48 * t605 + t439) * t592 + t216 * t282 / 0.2e1 + t249 * t127 / 0.2e1 + t131 * t570 + t302 * t354 / 0.2e1 + t352 * t128 / 0.2e1 + t372 * t267 / 0.2e1 + t76 * t281 / 0.2e1 + t95 * t569 + t613 * t553;
t396 = -t264 * t60 - t59 * t605 + t439;
t1 = t610 * t379 + t611 * t378 + t605 * t579 - (t136 + t135 + t105) * t344 / 0.4e1 - (t134 + t111 + t109) * t342 / 0.4e1 - t380 - t490 / 0.4e1 + t488 / 0.4e1 - t602 * t250 + t603 * t252 + t182 * t624 + (t358 - t416) * t493 / 0.4e1 + t284 * t428 + t486 * t634 + ((-t352 * t492 + t505) * pkin(4) + t396) * t594 + t342 * t581 + t344 * t582 + t517 * t585 + t519 * t586 + t520 * t588 + t518 * t589 + t435 * t553 + t492 * t550 + t492 * t551 + t132 * t470 + t382 + t608 * t347 * t465 + t618 * (-t250 * t624 + t252 * t634 + t66 * t557 + t67 * t560);
t12 = t249 * t281 + t352 * t282 + t372 * t354 + (-t416 / 0.2e1 + t358 / 0.2e1) * t379 + (pkin(4) * t284 - t356 / 0.2e1 + t359 / 0.2e1) * t378 + t544 * t632 + (-t286 / 0.2e1 - t289 / 0.2e1 - t291 / 0.2e1 + t288 / 0.2e1) * t344 + (-t287 / 0.2e1 - t290 / 0.2e1 - t292 / 0.2e1 + t285 / 0.2e1) * t342 + t635 * t275;
t401 = t1 * qJD(1) + t12 * qJD(2);
t384 = (t251 * t371 + t253 * t361) * t592 + (-t251 * t512 + t253 * t510) * t476 + mrSges(5,1) * t448 + mrSges(5,2) * t494 / 0.2e1 - t601;
t5 = t384 - t638;
t400 = t5 * qJD(1);
t392 = t215 * t594 + t592 * t75 + t601;
t10 = t342 * t441 + t344 * t440 + t392 - t637;
t33 = (t342 ^ 2 + t344 ^ 2) * t618 + (m(6) + m(7)) * (-t264 * t344 - t605 * t342);
t398 = -qJD(1) * t10 + qJD(2) * t33;
t393 = -t333 + ((qJ(6) + t361) * t345 + t60) * t593;
t23 = t474 + t393;
t397 = qJD(1) * t23 - qJD(4) * t351;
t145 = (t450 + t577) * m(7);
t82 = 0.2e1 * t624 * m(7) - t520;
t74 = t387 + t404;
t36 = t420 - t390;
t32 = t388 + t402;
t24 = -t393 + t474 - t525;
t18 = mrSges(7,2) * t449 + t523 / 0.2e1 + t389 - t408;
t11 = t392 + t424 + t637;
t8 = t381 + (-t179 / 0.2e1 + t580) * t342 + (-t177 / 0.2e1 - t184 / 0.2e1) * t344 + t383;
t6 = t384 + t638;
t2 = t380 + (t204 / 0.4e1 + t610) * t379 + ((t358 / 0.4e1 - t416 / 0.4e1 + t378 * t465) * t378 + (t550 + t551 + t379 * t465 + (t352 * t595 - t284 / 0.2e1) * pkin(4)) * t379) * t347 + t375 * t553 + t382 + (pkin(4) * t505 + t396) * t594 + t440 * t605 - t441 * t264 + (pkin(4) * t132 / 0.2e1 - t531 / 0.4e1 - t202 / 0.4e1 + t611) * t378 + (-t264 * t469 + t603) * t252 - (t469 * t605 + t602) * t250 + (-t136 / 0.4e1 - t135 / 0.4e1 - t105 / 0.4e1 + t582 + (t584 + t585) * mrSges(6,3) + (t584 + t589) * mrSges(7,2)) * t344 + (-t111 / 0.4e1 - t109 / 0.4e1 + t581 - t134 / 0.4e1 + (t583 + t586) * mrSges(6,3) + (t588 + t583) * mrSges(7,2)) * t342;
t14 = [qJD(2) * t3 + qJD(3) * t9 + qJD(4) * t4 + qJD(5) * t13 + qJD(6) * t21, t8 * qJD(3) + t2 * qJD(4) + t11 * qJD(5) + t18 * qJD(6) + t521 + ((m(5) * t372 - t513 * t591 - mrSges(4,1) + t417) * t615 + (-t511 * t591 + mrSges(4,2)) * t302 + (m(5) * t370 + mrSges(5,3)) * (-t169 * t378 + t170 * t379) - Ifges(3,6) * t545 + Ifges(3,5) * t546 - t62 * t519 - t49 * t520 - t50 * t518 - t279 * t491 - mrSges(3,1) * t473 + t106 * t560 + t286 * t577 + t288 * t578 + t104 * t558 + t201 * t547 + t203 * t548 + t460 * t540 + t462 * t541 + (t292 + t290) * t576 + (t110 + t108) * t557 + (-m(6) * t61 + m(7) * t50 - t179 + t180) * t264 + (m(7) * t75 + t129) * t249 - t372 * t268 + t352 * t130 - Ifges(4,5) * t345 + Ifges(4,6) * t347 + (m(6) * t62 + m(7) * t49 + t623) * t605 + (Ifges(5,5) * t378 + Ifges(5,6) * t379 - t342 * t630 - t617 * t344) * t552 + (t284 + t632) * t215 + t358 * t447 + t356 * t448 + mrSges(3,2) * t471 + t277 * t360 + t61 * t517 + t75 * t283) * qJD(2), t8 * qJD(2) + t6 * qJD(4) + t36 * qJD(5) + t145 * qJD(6) + t514 + 0.2e1 * (t592 + t594) * qJD(3) * (t251 * t342 - t344 * t253) t515 + t2 * qJD(2) + t6 * qJD(3) + (-t168 * mrSges(5,1) - t167 * mrSges(5,2) + Ifges(5,5) * t493 + Ifges(5,6) * t492 - t250 * t426 + t252 * t427 - t361 * t528 - t371 * t525 + t599 * t67 + t600 * t66 + t612) * qJD(4) + t32 * qJD(5) + t24 * qJD(6), qJD(2) * t11 + qJD(3) * t36 + qJD(4) * t32 + t509, qJD(2) * t18 + qJD(3) * t145 + qJD(4) * t24 + t508; -qJD(3) * t7 + qJD(4) * t1 - qJD(5) * t10 + qJD(6) * t19 - t521, qJD(4) * t12 + qJD(5) * t33 + qJD(6) * t96, -t529 (-mrSges(5,1) * t360 + mrSges(5,2) * t491 - t264 * t599 + t342 * t427 + t344 * t426 + t361 * t518 - t371 * t520 + t600 * t605 + t435 + t613) * qJD(4) + t74 * qJD(5) + t82 * qJD(6) + t401, qJD(4) * t74 + t398, qJD(4) * t82 + t415; qJD(2) * t7 - qJD(4) * t5 - qJD(5) * t35 + qJD(6) * t146 - t514, t529, 0 (-t354 + t391 - t403 + t405) * qJD(4) - t607 - t400, -t413, -qJD(4) * t273 + t477; -qJD(2) * t1 + qJD(3) * t5 - qJD(5) * t25 - qJD(6) * t23 - t515, -qJD(5) * t63 - t401, t400, t351 * qJD(6), -t414, -t397; qJD(2) * t10 + qJD(3) * t35 + qJD(4) * t25 - qJD(6) * t124 - t509, qJD(4) * t63 - t398 + t607, t413, t414, 0, t409; -qJD(2) * t19 - qJD(3) * t146 + qJD(4) * t23 + qJD(5) * t124 - t508, -qJD(5) * t273 - t415, -t477, t397, -t409, 0;];
Cq  = t14;
