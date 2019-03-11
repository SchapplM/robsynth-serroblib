% Calculate vector of inverse dynamics joint torques for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRRRPR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR7_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:39:03
% EndTime: 2019-03-08 23:40:17
% DurationCPUTime: 45.52s
% Computational Cost: add. (17699->975), mult. (44399->1403), div. (0->0), fcn. (37670->18), ass. (0->423)
t367 = sin(pkin(6));
t373 = sin(qJ(3));
t376 = cos(qJ(3));
t508 = cos(pkin(7));
t528 = sin(qJ(2));
t433 = t508 * t528;
t529 = cos(qJ(2));
t389 = -t373 * t433 + t376 * t529;
t274 = t389 * t367;
t253 = qJD(1) * t274;
t452 = t376 * t508;
t366 = sin(pkin(7));
t498 = t366 * t373;
t316 = pkin(2) * t452 - pkin(9) * t498;
t300 = t316 * qJD(3);
t587 = -t253 + t300;
t453 = t373 * t508;
t496 = t366 * t376;
t317 = pkin(2) * t453 + pkin(9) * t496;
t285 = pkin(10) * t508 + t317;
t427 = -pkin(3) * t376 - pkin(10) * t373;
t286 = (-pkin(2) + t427) * t366;
t406 = (pkin(3) * t373 - pkin(10) * t376) * t366;
t299 = qJD(3) * t406;
t372 = sin(qJ(4));
t375 = cos(qJ(4));
t473 = t367 * t528;
t441 = qJD(1) * t473;
t415 = t366 * t441;
t483 = qJD(4) * t375;
t484 = qJD(4) * t372;
t594 = -t285 * t484 + t286 * t483 + t587 * t375 + (t299 - t415) * t372;
t388 = t373 * t529 + t376 * t433;
t273 = t388 * t367;
t639 = qJD(1) * t273 - t317 * qJD(3);
t487 = qJD(3) * t373;
t638 = -(qJ(5) * t487 - qJD(5) * t376) * t366 - t594;
t314 = t372 * t498 - t375 * t508;
t486 = qJD(3) * t376;
t467 = t366 * t486;
t240 = -qJD(4) * t314 + t375 * t467;
t497 = t366 * t375;
t315 = t372 * t508 + t373 * t497;
t241 = qJD(4) * t315 + t372 * t467;
t637 = t241 * pkin(4) - t240 * qJ(5) - t315 * qJD(5) - t639;
t417 = pkin(4) * t372 - qJ(5) * t375;
t474 = t367 * t529;
t442 = qJD(1) * t474;
t334 = qJD(2) * pkin(2) + t442;
t450 = t508 * t334;
t488 = qJD(2) * t376;
t369 = cos(pkin(6));
t490 = qJD(1) * t369;
t489 = qJD(2) * t366;
t319 = pkin(9) * t489 + t441;
t492 = t376 * t319;
t636 = -t373 * t450 - t492 - (t373 * t490 + t417 * t488) * t366 + qJD(4) * t417 - qJD(5) * t372;
t365 = sin(pkin(13));
t368 = cos(pkin(13));
t604 = t638 * t365 + t637 * t368;
t603 = t637 * t365 - t638 * t368;
t394 = t366 * t490 + t450;
t195 = -t373 * t319 + t376 * t394;
t298 = qJD(2) * t406;
t144 = t375 * t195 + t372 * t298;
t470 = t373 * t489;
t133 = qJ(5) * t470 + t144;
t635 = -t368 * t133 + t636 * t365;
t481 = pkin(10) * t484;
t599 = t636 * t368 + (t133 + t481) * t365;
t616 = -m(6) - m(5);
t468 = t366 * t487;
t207 = t240 * t368 + t365 * t468;
t634 = pkin(5) * t241 - pkin(11) * t207 + t604;
t206 = -t240 * t365 + t368 * t468;
t633 = -pkin(11) * t206 - t603;
t493 = t375 * t376;
t251 = (t365 * t373 + t368 * t493) * t489;
t469 = t366 * t488;
t443 = t372 * t469;
t494 = t368 * t375;
t632 = -pkin(5) * t443 + pkin(11) * t251 + (pkin(5) * t372 - pkin(11) * t494) * qJD(4) + t599;
t250 = (-t365 * t493 + t368 * t373) * t489;
t495 = t368 * t372;
t500 = t365 * t375;
t631 = -pkin(11) * t250 + (-pkin(10) * t495 - pkin(11) * t500) * qJD(4) + t635;
t361 = pkin(5) * t368 + pkin(4);
t364 = pkin(13) + qJ(6);
t362 = sin(t364);
t363 = cos(t364);
t425 = -mrSges(6,1) * t368 + mrSges(6,2) * t365;
t630 = m(6) * pkin(4) + m(7) * t361 + mrSges(7,1) * t363 - mrSges(7,2) * t362 - t425;
t522 = pkin(11) + qJ(5);
t629 = m(6) * qJ(5) + m(7) * t522 + mrSges(6,3) + mrSges(7,3);
t412 = qJD(2) * t508 + qJD(3);
t272 = t372 * t412 + t375 * t470;
t342 = qJD(4) - t469;
t220 = t272 * t368 + t342 * t365;
t271 = t372 * t470 - t375 * t412;
t178 = -pkin(3) * t412 - t195;
t109 = t271 * pkin(4) - t272 * qJ(5) + t178;
t196 = t373 * t394 + t492;
t179 = pkin(10) * t412 + t196;
t454 = t369 * t508;
t350 = qJD(1) * t454;
t216 = t350 + (qJD(2) * t427 - t334) * t366;
t108 = t179 * t375 + t216 * t372;
t95 = qJ(5) * t342 + t108;
t46 = t368 * t109 - t365 * t95;
t32 = pkin(5) * t271 - pkin(11) * t220 + t46;
t449 = -t272 * t365 + t368 * t342;
t47 = t365 * t109 + t368 * t95;
t35 = pkin(11) * t449 + t47;
t371 = sin(qJ(6));
t374 = cos(qJ(6));
t14 = t32 * t374 - t35 * t371;
t15 = t32 * t371 + t35 * t374;
t513 = t108 * mrSges(5,3);
t606 = t342 * Ifges(5,6);
t628 = t513 - t178 * mrSges(5,1) - t46 * mrSges(6,1) - t14 * mrSges(7,1) + t47 * mrSges(6,2) + t15 * mrSges(7,2) + t606 / 0.2e1;
t627 = Ifges(6,6) * t449;
t626 = t220 * Ifges(6,5);
t125 = t220 * t374 + t371 * t449;
t261 = qJD(6) + t271;
t621 = -t220 * t371 + t374 * t449;
t605 = t125 * Ifges(7,5) + Ifges(7,6) * t621 + t271 * Ifges(6,3) + t261 * Ifges(7,3) + t626 + t627;
t408 = qJD(2) * t441;
t306 = qJDD(1) * t474 - t408;
t280 = qJDD(2) * pkin(2) + t306;
t625 = qJDD(1) * t366 * t369 + t508 * t280;
t426 = -mrSges(5,1) * t375 + mrSges(5,2) * t372;
t577 = m(7) - t616;
t624 = pkin(3) * t577 + t372 * t629 + t375 * t630 + mrSges(4,1) - t426;
t346 = qJD(2) * t442;
t307 = qJDD(1) * t473 + t346;
t623 = pkin(9) * qJDD(2) * t366 + qJD(3) * t394 + t307;
t424 = t365 * mrSges(6,1) + t368 * mrSges(6,2);
t475 = pkin(5) * t365 + pkin(10);
t568 = -m(7) * t475 - t362 * mrSges(7,1) - t363 * mrSges(7,2) + mrSges(4,2) - mrSges(5,3) - t424;
t448 = mrSges(4,3) * t470;
t572 = -m(5) * t178 + mrSges(4,1) * t412 - mrSges(5,1) * t271 - mrSges(5,2) * t272 - t448;
t622 = -m(4) * t195 - t572;
t482 = qJD(2) * qJD(3);
t305 = (qJDD(2) * t373 + t376 * t482) * t366;
t353 = qJDD(2) * t508 + qJDD(3);
t485 = qJD(4) * t271;
t176 = t375 * t305 + t372 * t353 - t485;
t304 = (-qJDD(2) * t376 + t373 * t482) * t366;
t289 = qJDD(4) + t304;
t131 = t176 * t368 + t289 * t365;
t177 = qJD(4) * t272 + t372 * t305 - t375 * t353;
t229 = qJDD(1) * t454 - t280 * t366;
t142 = pkin(3) * t304 - pkin(10) * t305 + t229;
t93 = -t319 * t487 + t373 * t625 + t376 * t623;
t86 = pkin(10) * t353 + t93;
t29 = t372 * t142 - t179 * t484 + t216 * t483 + t375 * t86;
t23 = qJ(5) * t289 + qJD(5) * t342 + t29;
t94 = -t319 * t486 - t373 * t623 + t376 * t625;
t87 = -t353 * pkin(3) - t94;
t36 = t177 * pkin(4) - t176 * qJ(5) - t272 * qJD(5) + t87;
t7 = -t23 * t365 + t368 * t36;
t5 = pkin(5) * t177 - pkin(11) * t131 + t7;
t130 = -t176 * t365 + t289 * t368;
t8 = t368 * t23 + t365 * t36;
t6 = pkin(11) * t130 + t8;
t1 = qJD(6) * t14 + t371 * t5 + t374 * t6;
t2 = -qJD(6) * t15 - t371 * t6 + t374 * t5;
t620 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t542 = t177 / 0.2e1;
t547 = t131 / 0.2e1;
t619 = Ifges(6,1) * t547 + Ifges(6,5) * t542;
t511 = t272 * Ifges(5,4);
t159 = -t271 * Ifges(5,2) + t511 + t606;
t537 = t261 / 0.2e1;
t549 = t125 / 0.2e1;
t551 = t621 / 0.2e1;
t618 = Ifges(7,5) * t549 + Ifges(7,6) * t551 + Ifges(7,3) * t537 - t159 / 0.2e1;
t40 = qJD(6) * t621 + t130 * t371 + t131 * t374;
t563 = t40 / 0.2e1;
t41 = -qJD(6) * t125 + t130 * t374 - t131 * t371;
t562 = t41 / 0.2e1;
t548 = t130 / 0.2e1;
t168 = qJDD(6) + t177;
t545 = t168 / 0.2e1;
t544 = t176 / 0.2e1;
t543 = -t177 / 0.2e1;
t532 = t289 / 0.2e1;
t237 = t315 * t368 - t365 * t496;
t284 = -pkin(3) * t508 - t316;
t184 = t314 * pkin(4) - t315 * qJ(5) + t284;
t198 = t375 * t285 + t372 * t286;
t187 = -qJ(5) * t496 + t198;
t89 = t368 * t184 - t187 * t365;
t63 = pkin(5) * t314 - pkin(11) * t237 + t89;
t236 = -t315 * t365 - t368 * t496;
t90 = t365 * t184 + t368 * t187;
t75 = pkin(11) * t236 + t90;
t28 = t371 * t63 + t374 * t75;
t615 = -qJD(6) * t28 + t371 * t633 + t374 * t634;
t27 = -t371 * t75 + t374 * t63;
t614 = qJD(6) * t27 + t371 * t634 - t374 * t633;
t613 = t93 * mrSges(4,2);
t612 = t94 * mrSges(4,1);
t11 = Ifges(7,5) * t40 + Ifges(7,6) * t41 + Ifges(7,3) * t168;
t610 = t131 * Ifges(6,5) + t130 * Ifges(6,6) + t177 * Ifges(6,3) + t11;
t338 = -pkin(4) * t375 - qJ(5) * t372 - pkin(3);
t321 = t368 * t338;
t230 = -pkin(11) * t495 + t321 + (-pkin(10) * t365 - pkin(5)) * t375;
t270 = pkin(10) * t494 + t365 * t338;
t501 = t365 * t372;
t242 = -pkin(11) * t501 + t270;
t145 = t230 * t374 - t242 * t371;
t609 = qJD(6) * t145 + t371 * t632 + t374 * t631;
t146 = t230 * t371 + t242 * t374;
t608 = -qJD(6) * t146 - t371 * t631 + t374 * t632;
t607 = t342 * Ifges(5,5);
t339 = t522 * t365;
t340 = t522 * t368;
t244 = -t339 * t371 + t340 * t374;
t324 = t365 * t374 + t368 * t371;
t504 = t271 * t368;
t107 = -t372 * t179 + t216 * t375;
t193 = pkin(4) * t272 + qJ(5) * t271;
t68 = -t107 * t365 + t368 * t193;
t50 = pkin(5) * t272 + pkin(11) * t504 + t68;
t505 = t271 * t365;
t69 = t368 * t107 + t365 * t193;
t59 = pkin(11) * t505 + t69;
t602 = -qJD(5) * t324 - qJD(6) * t244 + t371 * t59 - t374 * t50;
t243 = -t339 * t374 - t340 * t371;
t416 = t365 * t371 - t368 * t374;
t601 = -qJD(5) * t416 + qJD(6) * t243 - t371 * t50 - t374 * t59;
t135 = mrSges(5,1) * t289 - mrSges(5,3) * t176;
t62 = -t130 * mrSges(6,1) + t131 * mrSges(6,2);
t600 = t62 - t135;
t598 = -t368 * t481 + t635;
t143 = -t372 * t195 + t298 * t375;
t134 = -pkin(4) * t470 - t143;
t595 = pkin(5) * t250 + t475 * t483 - t134;
t132 = -mrSges(6,1) * t449 + mrSges(6,2) * t220;
t226 = mrSges(5,1) * t342 - mrSges(5,3) * t272;
t593 = -t132 + t226;
t170 = t250 * t374 - t251 * t371;
t312 = t416 * qJD(6);
t213 = t312 * t372 - t324 * t483;
t592 = t170 - t213;
t171 = t250 * t371 + t251 * t374;
t313 = t324 * qJD(6);
t212 = -t313 * t372 - t416 * t483;
t591 = t171 - t212;
t174 = t324 * t271;
t590 = t174 + t313;
t175 = t416 * t271;
t589 = t175 + t312;
t506 = sin(pkin(12));
t429 = t506 * t528;
t507 = cos(pkin(12));
t432 = t507 * t529;
t385 = -t369 * t432 + t429;
t456 = t367 * t507;
t586 = t366 * t456 + t385 * t508;
t430 = t506 * t529;
t431 = t507 * t528;
t386 = t369 * t430 + t431;
t455 = t367 * t506;
t585 = -t366 * t455 + t386 * t508;
t434 = t508 * t529;
t584 = -t528 * t373 + t376 * t434;
t583 = t443 - t484;
t103 = t220 * Ifges(6,4) + Ifges(6,2) * t449 + Ifges(6,6) * t271;
t554 = -t103 / 0.2e1;
t581 = -t107 * mrSges(5,3) + t365 * t554;
t30 = t142 * t375 - t179 * t483 - t216 * t484 - t372 * t86;
t24 = -pkin(4) * t289 + qJDD(5) - t30;
t91 = -pkin(4) * t342 + qJD(5) - t107;
t580 = t24 * t372 + t91 * t483;
t579 = t29 * t375 - t30 * t372;
t578 = -t365 * t7 + t368 * t8;
t258 = -t334 * t366 + t350;
t576 = (-t258 * (mrSges(4,1) * t373 + mrSges(4,2) * t376) - t412 * (Ifges(4,5) * t376 - Ifges(4,6) * t373) / 0.2e1) * t366;
t575 = mrSges(5,1) + t630;
t570 = mrSges(5,2) - t629;
t571 = pkin(10) * t616 + t568;
t61 = -mrSges(7,1) * t621 + mrSges(7,2) * t125;
t72 = -pkin(5) * t449 + t91;
t569 = m(5) * t107 - m(6) * t91 - m(7) * t72 + t593 - t61;
t566 = t366 ^ 2;
t377 = qJD(2) ^ 2;
t565 = Ifges(7,4) * t563 + Ifges(7,2) * t562 + Ifges(7,6) * t545;
t564 = Ifges(7,1) * t563 + Ifges(7,4) * t562 + Ifges(7,5) * t545;
t44 = t131 * Ifges(6,4) + t130 * Ifges(6,2) + t177 * Ifges(6,6);
t561 = t44 / 0.2e1;
t560 = Ifges(6,4) * t548 + t619;
t515 = Ifges(7,4) * t125;
t52 = Ifges(7,2) * t621 + t261 * Ifges(7,6) + t515;
t559 = -t52 / 0.2e1;
t558 = t52 / 0.2e1;
t122 = Ifges(7,4) * t621;
t53 = t125 * Ifges(7,1) + t261 * Ifges(7,5) + t122;
t557 = -t53 / 0.2e1;
t556 = t53 / 0.2e1;
t555 = Ifges(5,1) * t544 + Ifges(5,4) * t543 + Ifges(5,5) * t532;
t104 = t220 * Ifges(6,1) + Ifges(6,4) * t449 + Ifges(6,5) * t271;
t553 = t104 / 0.2e1;
t552 = -t621 / 0.2e1;
t550 = -t125 / 0.2e1;
t541 = -t449 / 0.2e1;
t540 = -t220 / 0.2e1;
t538 = -t261 / 0.2e1;
t536 = -t271 / 0.2e1;
t535 = t271 / 0.2e1;
t533 = t272 / 0.2e1;
t527 = pkin(2) * t366;
t521 = Ifges(4,4) * t373;
t520 = Ifges(4,4) * t376;
t519 = Ifges(5,4) * t372;
t518 = Ifges(5,4) * t375;
t517 = Ifges(6,4) * t365;
t516 = Ifges(6,4) * t368;
t310 = t369 * t431 + t430;
t503 = t310 * t366;
t311 = -t369 * t429 + t432;
t502 = t311 * t366;
t499 = t366 * t372;
t445 = t366 * t473;
t491 = pkin(2) * t474 + pkin(9) * t445;
t480 = pkin(10) * t483;
t478 = Ifges(5,5) * t176 - Ifges(5,6) * t177 + Ifges(5,3) * t289;
t477 = t274 * pkin(3) + t491;
t476 = Ifges(4,5) * t305 - Ifges(4,6) * t304 + Ifges(4,3) * t353;
t465 = t498 / 0.2e1;
t16 = -t41 * mrSges(7,1) + t40 * mrSges(7,2);
t461 = -t489 / 0.2e1;
t457 = t483 / 0.2e1;
t197 = -t372 * t285 + t286 * t375;
t447 = mrSges(4,3) * t469;
t439 = t376 * t461;
t438 = -t385 * pkin(2) + pkin(9) * t503;
t437 = -t386 * pkin(2) + pkin(9) * t502;
t188 = pkin(4) * t496 - t197;
t423 = Ifges(5,1) * t375 - t519;
t422 = Ifges(6,1) * t368 - t517;
t421 = -Ifges(5,2) * t372 + t518;
t420 = -Ifges(6,2) * t365 + t516;
t419 = Ifges(5,5) * t375 - Ifges(5,6) * t372;
t418 = Ifges(6,5) * t368 - Ifges(6,6) * t365;
t390 = t373 * t434 + t376 * t528;
t233 = t367 * t390 + t369 * t498;
t392 = -t366 * t474 + t454;
t186 = t233 * t375 + t372 * t392;
t232 = -t367 * t584 - t369 * t496;
t114 = -t186 * t365 + t232 * t368;
t115 = t186 * t368 + t232 * t365;
t57 = t114 * t374 - t115 * t371;
t58 = t114 * t371 + t115 * t374;
t152 = t236 * t374 - t237 * t371;
t153 = t236 * t371 + t237 * t374;
t414 = qJD(2) * t445;
t209 = -t310 * t453 - t376 * t385;
t411 = t209 * pkin(3) + t438;
t211 = -t311 * t453 - t376 * t386;
t410 = t211 * pkin(3) + t437;
t121 = -t285 * t483 - t286 * t484 + t299 * t375 - t372 * t300;
t405 = t178 * (mrSges(5,1) * t372 + mrSges(5,2) * t375);
t404 = t366 * (-mrSges(4,1) * t376 + mrSges(4,2) * t373);
t181 = t310 * t376 - t373 * t586;
t378 = t366 * t385 - t456 * t508;
t116 = t181 * t372 - t375 * t378;
t183 = t311 * t376 - t373 * t585;
t379 = t366 * t386 + t455 * t508;
t118 = t183 * t372 - t375 * t379;
t185 = t233 * t372 - t375 * t392;
t400 = -g(1) * t118 - g(2) * t116 - g(3) * t185;
t397 = t373 * t566 * (Ifges(4,1) * t376 - t521);
t112 = -pkin(4) * t468 - t121;
t387 = Ifges(4,6) * t508 + (t376 * Ifges(4,2) + t521) * t366;
t347 = Ifges(4,4) * t469;
t329 = t475 * t372;
t297 = qJD(2) * t404;
t296 = -mrSges(4,2) * t412 + t447;
t294 = t416 * t372;
t293 = t324 * t372;
t269 = -pkin(10) * t500 + t321;
t259 = Ifges(5,4) * t271;
t246 = Ifges(4,1) * t470 + Ifges(4,5) * t412 + t347;
t245 = Ifges(4,6) * qJD(3) + qJD(2) * t387;
t239 = mrSges(4,1) * t353 - mrSges(4,3) * t305;
t238 = -mrSges(4,2) * t353 - mrSges(4,3) * t304;
t225 = -mrSges(5,2) * t342 - mrSges(5,3) * t271;
t217 = mrSges(4,1) * t304 + mrSges(4,2) * t305;
t210 = t311 * t452 - t373 * t386;
t208 = t310 * t452 - t373 * t385;
t182 = t311 * t373 + t376 * t585;
t180 = t310 * t373 + t376 * t586;
t163 = t369 * t467 + (t389 * qJD(2) + qJD(3) * t584) * t367;
t162 = t369 * t468 + (qJD(2) * t388 + qJD(3) * t390) * t367;
t160 = t272 * Ifges(5,1) - t259 + t607;
t158 = t272 * Ifges(5,5) - t271 * Ifges(5,6) + t342 * Ifges(5,3);
t155 = mrSges(6,1) * t271 - mrSges(6,3) * t220;
t154 = -mrSges(6,2) * t271 + mrSges(6,3) * t449;
t136 = -mrSges(5,2) * t289 - mrSges(5,3) * t177;
t128 = -pkin(5) * t236 + t188;
t119 = t183 * t375 + t372 * t379;
t117 = t181 * t375 + t372 * t378;
t97 = mrSges(7,1) * t261 - mrSges(7,3) * t125;
t96 = -mrSges(7,2) * t261 + mrSges(7,3) * t621;
t92 = mrSges(5,1) * t177 + mrSges(5,2) * t176;
t83 = -pkin(5) * t505 + t108;
t82 = -qJD(4) * t185 + t163 * t375 + t372 * t414;
t79 = t176 * Ifges(5,4) - t177 * Ifges(5,2) + t289 * Ifges(5,6);
t78 = mrSges(6,1) * t177 - mrSges(6,3) * t131;
t77 = -mrSges(6,2) * t177 + mrSges(6,3) * t130;
t76 = -pkin(5) * t206 + t112;
t65 = -qJD(6) * t153 + t206 * t374 - t207 * t371;
t64 = qJD(6) * t152 + t206 * t371 + t207 * t374;
t55 = t162 * t365 + t368 * t82;
t54 = t162 * t368 - t365 * t82;
t26 = -mrSges(7,2) * t168 + mrSges(7,3) * t41;
t25 = mrSges(7,1) * t168 - mrSges(7,3) * t40;
t17 = -pkin(5) * t130 + t24;
t10 = -qJD(6) * t58 - t371 * t55 + t374 * t54;
t9 = qJD(6) * t57 + t371 * t54 + t374 * t55;
t3 = [m(4) * (t196 * t163 + t229 * t392 + t93 * t233 + t258 * t414) + m(5) * (t108 * t82 + t186 * t29) + (-m(2) - m(3) - m(4) - t577) * g(3) + m(3) * (t369 ^ 2 * qJDD(1) + (t306 * t529 + t307 * t528) * t367) - t569 * (qJD(4) * t186 + t163 * t372 - t375 * t414) + (-m(4) * t94 + m(5) * t87 - t239 + t92) * t232 + (qJDD(2) * t474 - t377 * t473) * mrSges(3,1) + t55 * t154 + t54 * t155 + t392 * t217 + t115 * t77 + t114 * t78 + t9 * t96 + t10 * t97 + t57 * t25 + t58 * t26 + m(7) * (t1 * t58 + t10 * t14 + t15 * t9 + t2 * t57) + m(6) * (t114 * t7 + t115 * t8 + t46 * t54 + t47 * t55) + (-qJDD(2) * t473 - t377 * t474) * mrSges(3,2) + (-m(5) * t30 + m(6) * t24 + m(7) * t17 + t16 + t600) * t185 + t186 * t136 + t82 * t225 + t233 * t238 + t622 * t162 + t163 * t296 + t297 * t414 + m(2) * qJDD(1); (-m(4) * t437 - m(7) * t410 + mrSges(3,1) * t386 - t211 * mrSges(4,1) + t311 * mrSges(3,2) - mrSges(4,3) * t502 + t616 * (pkin(10) * t210 + t410) - t575 * (t211 * t375 + t311 * t499) + t568 * t210 + t570 * (t211 * t372 - t311 * t497)) * g(1) + (-m(4) * t438 - m(7) * t411 + mrSges(3,1) * t385 - t209 * mrSges(4,1) + t310 * mrSges(3,2) - mrSges(4,3) * t503 + t616 * (pkin(10) * t208 + t411) - t575 * (t209 * t375 + t310 * t499) + t568 * t208 + t570 * (t209 * t372 - t310 * t497)) * g(2) + (-m(4) * t491 - t274 * mrSges(4,1) - mrSges(4,3) * t445 - (mrSges(3,1) * t529 - mrSges(3,2) * t528) * t367 - m(7) * t477 + t616 * (pkin(10) * t273 + t477) - t575 * (t274 * t375 + t372 * t445) + t568 * t273 + t570 * (t274 * t372 - t375 * t445)) * g(3) + (Ifges(6,1) * t237 + Ifges(6,4) * t236) * t547 + t220 * (Ifges(6,1) * t207 + Ifges(6,4) * t206) / 0.2e1 + (Ifges(5,1) * t240 + Ifges(5,5) * t468) * t533 + (Ifges(5,1) * t315 - Ifges(5,5) * t496) * t544 + (Ifges(7,1) * t64 + Ifges(7,4) * t65) * t549 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t563 + (t1 * t152 - t14 * t64 + t15 * t65 - t153 * t2) * mrSges(7,3) + (t206 * t47 - t207 * t46 + t236 * t8 - t237 * t7) * mrSges(6,3) + (-t108 * t468 + t178 * t240 + t29 * t496 + t315 * t87) * mrSges(5,2) + (Ifges(5,4) * t240 + Ifges(5,6) * t468) * t536 + (Ifges(5,4) * t315 - Ifges(5,6) * t496) * t543 + (t627 / 0.2e1 - Ifges(5,4) * t533 + Ifges(6,3) * t535 - Ifges(5,2) * t536 + t605 / 0.2e1 + t626 / 0.2e1 + t618 - t628) * t241 + (t408 + t306) * mrSges(3,1) + (t107 * t121 + t108 * t594 + t197 * t30 + t198 * t29 + t284 * t87) * m(5) - t622 * t639 + (t158 * t465 - t576) * qJD(3) + t397 * t482 / 0.2e1 + t587 * t296 - t217 * t527 - t478 * t496 / 0.2e1 + t30 * (-mrSges(5,1) * t496 - mrSges(5,3) * t315) + (t196 * t587 - t229 * t527 - t258 * t415 + t316 * t94 + t317 * t93) * m(4) + t569 * (t253 * t372 - t375 * t415) + t508 * t612 + t17 * (-mrSges(7,1) * t152 + mrSges(7,2) * t153) + (t566 * qJD(2) * (-Ifges(4,2) * t373 + t520) + t366 * t246) * t486 / 0.2e1 + t229 * t404 + (Ifges(4,1) * t305 - Ifges(4,4) * t304 + Ifges(4,5) * t353) * t465 + (Ifges(4,4) * t305 - Ifges(4,2) * t304 + Ifges(4,6) * t353) * t496 / 0.2e1 + t128 * t16 + t112 * t132 + t305 * (Ifges(4,5) * t508 + (t373 * Ifges(4,1) + t520) * t366) / 0.2e1 - t297 * t415 + t508 * t476 / 0.2e1 + t353 * (Ifges(4,3) * t508 + (Ifges(4,5) * t373 + Ifges(4,6) * t376) * t366) / 0.2e1 - t304 * t387 / 0.2e1 + (Ifges(7,5) * t153 + Ifges(7,6) * t152) * t545 + (Ifges(7,5) * t64 + Ifges(7,6) * t65) * t537 - t245 * t468 / 0.2e1 + t107 * (mrSges(5,1) * t468 - mrSges(5,3) * t240) + t89 * t78 + t90 * t77 + t76 * t61 + t72 * (-mrSges(7,1) * t65 + mrSges(7,2) * t64) + (-t195 * t467 - t196 * t468 + t496 * t93 - t498 * t94) * mrSges(4,3) + Ifges(3,3) * qJDD(2) + t27 * t25 + t28 * t26 + (Ifges(6,5) * t207 + Ifges(6,6) * t206) * t535 + (Ifges(6,5) * t237 + Ifges(6,6) * t236) * t542 + (Ifges(7,4) * t64 + Ifges(7,2) * t65) * t551 + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t562 + t449 * (Ifges(6,4) * t207 + Ifges(6,2) * t206) / 0.2e1 + (Ifges(6,4) * t237 + Ifges(6,2) * t236) * t548 + (t346 - t307) * mrSges(3,2) + t342 * (Ifges(5,5) * t240 + Ifges(5,3) * t468) / 0.2e1 + (Ifges(5,5) * t315 - Ifges(5,3) * t496) * t532 + t603 * t154 + t604 * t155 + (t112 * t91 + t188 * t24 + t46 * t604 + t47 * t603 + t7 * t89 + t8 * t90) * m(6) + t188 * t62 + t594 * t225 + t197 * t135 + t198 * t136 + t206 * t103 / 0.2e1 + t91 * (-mrSges(6,1) * t206 + mrSges(6,2) * t207) - t508 * t613 + t121 * t226 + t614 * t96 + t615 * t97 + (t1 * t28 + t128 * t17 + t14 * t615 + t15 * t614 + t2 * t27 + t72 * t76) * m(7) + t207 * t553 + t315 * t555 + t64 * t556 + t65 * t558 + t237 * t560 + t236 * t561 + t153 * t564 + t152 * t565 + (-t29 * mrSges(5,3) - Ifges(5,6) * t532 + t610 / 0.2e1 + Ifges(6,3) * t542 - Ifges(5,2) * t543 - Ifges(5,4) * t544 + Ifges(7,3) * t545 + Ifges(6,5) * t547 + Ifges(6,6) * t548 + Ifges(7,6) * t562 + Ifges(7,5) * t563 - t8 * mrSges(6,2) + t7 * mrSges(6,1) - t79 / 0.2e1 + t87 * mrSges(5,1) + t620) * t314 + t24 * (-mrSges(6,1) * t236 + mrSges(6,2) * t237) + t240 * t160 / 0.2e1 + t284 * t92 + t316 * t239 + t317 * t238; (t342 * (Ifges(5,3) * t373 + t376 * t419) + t272 * (Ifges(5,5) * t373 + t376 * t423) + t373 * t158) * t461 + (-Ifges(4,2) * t470 + t375 * t160 + t246 + t347) * t439 + (-t134 + t480) * t132 + (Ifges(7,4) * t212 + Ifges(7,2) * t213) * t551 + (-pkin(3) * t87 - t107 * t143 - t108 * t144) * m(5) + (-t134 * t91 + t269 * t7 + t270 * t8 + t46 * t599 + t47 * t598) * m(6) + (((-t107 * t375 - t108 * t372) * qJD(4) + t579) * m(5) + t600 * t372 + t375 * t136 + t580 * m(6)) * pkin(10) + (-t480 - t143) * t226 + t579 * mrSges(5,3) + t580 * t424 + t581 * t483 + (Ifges(7,1) * t212 + Ifges(7,4) * t213) * t549 + t576 * qJD(2) + (-t108 * (-mrSges(5,3) * t372 * t376 - mrSges(5,2) * t373) - t107 * (mrSges(5,1) * t373 - mrSges(5,3) * t493)) * t489 + (-t481 - t144) * t225 - t44 * t501 / 0.2e1 + t8 * (mrSges(6,2) * t375 - mrSges(6,3) * t501) + (Ifges(6,3) * t372 + t375 * t418) * t485 / 0.2e1 - t421 * t485 / 0.2e1 + t612 - t613 + (t368 * t457 - t251 / 0.2e1) * t104 + (-t296 + t447) * t195 + t160 * t457 + t145 * t25 + t146 * t26 + (t271 * (Ifges(5,6) * t373 + t376 * t421) + t373 * t245) * t489 / 0.2e1 + (t220 * (Ifges(6,5) * t372 + t375 * t422) + t449 * (Ifges(6,6) * t372 + t375 * t420) + t342 * t419 + t272 * t423) * qJD(4) / 0.2e1 + t1 * (mrSges(7,2) * t375 - mrSges(7,3) * t293) + (-Ifges(7,5) * t294 - Ifges(7,6) * t293 - Ifges(7,3) * t375) * t545 + (-Ifges(7,4) * t294 - Ifges(7,2) * t293 - Ifges(7,6) * t375) * t562 + (-Ifges(7,1) * t294 - Ifges(7,4) * t293 - Ifges(7,5) * t375) * t563 + t17 * (mrSges(7,1) * t293 - mrSges(7,2) * t294) + t2 * (-mrSges(7,1) * t375 + mrSges(7,3) * t294) + (t448 + t572) * t196 + (Ifges(7,5) * t212 + Ifges(7,6) * t213) * t537 + t476 + t7 * (-mrSges(6,1) * t375 - mrSges(6,3) * t495) - t377 * t397 / 0.2e1 - t405 * t469 - pkin(3) * t92 - t47 * (-mrSges(6,2) * t443 + mrSges(6,3) * t250) - t46 * (mrSges(6,1) * t443 - mrSges(6,3) * t251) + t159 * t443 / 0.2e1 + (t47 * (-mrSges(6,2) * t372 - mrSges(6,3) * t500) + t46 * (mrSges(6,1) * t372 - mrSges(6,3) * t494) + t405) * qJD(4) + t375 * t79 / 0.2e1 + t605 * (t372 * t439 + t484 / 0.2e1) + (Ifges(5,5) * t372 + Ifges(5,6) * t375) * t532 + (Ifges(6,5) * t251 + Ifges(6,6) * t250 + Ifges(6,3) * t443) * t536 + (Ifges(7,5) * t171 + Ifges(7,6) * t170 + Ifges(7,3) * t443) * t538 + t598 * t154 + t599 * t155 + (-t513 + t618) * t484 + (-mrSges(7,1) * t583 + mrSges(7,3) * t591) * t14 + (mrSges(7,1) * t592 - mrSges(7,2) * t591) * t72 + (mrSges(7,2) * t583 - mrSges(7,3) * t592) * t15 + t595 * t61 + t608 * t97 + t609 * t96 + (t1 * t146 + t14 * t608 + t145 * t2 + t15 * t609 + t17 * t329 + t595 * t72) * m(7) - t610 * t375 / 0.2e1 + (Ifges(6,1) * t251 + Ifges(6,4) * t250 + Ifges(6,5) * t443) * t540 + (Ifges(6,4) * t251 + Ifges(6,2) * t250 + Ifges(6,6) * t443) * t541 + (-Ifges(6,3) * t375 + t372 * t418) * t542 + (Ifges(5,2) * t375 + t519) * t543 + (Ifges(5,1) * t372 + t518) * t544 + (-Ifges(6,5) * t375 + t372 * t422) * t547 + (-Ifges(6,6) * t375 + t372 * t420) * t548 + (Ifges(7,1) * t171 + Ifges(7,4) * t170 + Ifges(7,5) * t443) * t550 + (Ifges(7,4) * t171 + Ifges(7,2) * t170 + Ifges(7,6) * t443) * t552 + t250 * t554 + t372 * t555 + t212 * t556 + t171 * t557 + t213 * t558 + t170 * t559 + t495 * t560 - t294 * t564 - t293 * t565 - t91 * (-mrSges(6,1) * t250 + mrSges(6,2) * t251) + t87 * t426 + t269 * t78 + t270 * t77 + (t182 * t624 + t183 * t571) * g(1) + (t232 * t624 + t233 * t571) * g(3) + (t180 * t624 + t181 * t571) * g(2) + t329 * t16; (Ifges(7,5) * t324 - Ifges(7,6) * t416) * t545 + (Ifges(7,4) * t324 - Ifges(7,2) * t416) * t562 + (Ifges(7,1) * t324 - Ifges(7,4) * t416) * t563 + t17 * (mrSges(7,1) * t416 + mrSges(7,2) * t324) + (-t1 * t416 + t14 * t589 - t15 * t590 - t2 * t324) * mrSges(7,3) - t416 * t565 + (Ifges(6,5) * t540 + Ifges(7,5) * t550 - Ifges(5,2) * t535 + Ifges(6,6) * t541 + Ifges(7,6) * t552 + Ifges(6,3) * t536 + Ifges(7,3) * t538 + t628) * t272 + (Ifges(7,4) * t175 + Ifges(7,2) * t174) * t552 + (Ifges(7,1) * t175 + Ifges(7,4) * t174) * t550 + (Ifges(6,2) * t548 + Ifges(6,6) * t542 + qJ(5) * t77 + qJD(5) * t154 + t561) * t368 + (-t46 * t504 - t47 * t505 + t578) * mrSges(6,3) + (-pkin(4) * t24 + (-t365 * t46 + t368 * t47) * qJD(5) + t578 * qJ(5) - t108 * t91 - t46 * t68 - t47 * t69) * m(6) + (t118 * t575 + t119 * t570) * g(1) + (t116 * t575 + t117 * t570) * g(2) + (t185 * t575 + t186 * t570) * g(3) + (-Ifges(7,1) * t312 - Ifges(7,4) * t313) * t549 + (-Ifges(7,4) * t312 - Ifges(7,2) * t313) * t551 + (-Ifges(7,5) * t312 - Ifges(7,6) * t313) * t537 + (-t259 + t160) * t535 + (mrSges(7,1) * t590 - mrSges(7,2) * t589) * t72 + t478 - t68 * t155 - t69 * t154 + t516 * t547 + t517 * t548 - t83 * t61 - pkin(4) * t62 - t29 * mrSges(5,2) + t30 * mrSges(5,1) - t361 * t16 + t159 * t533 + t601 * t96 + t602 * t97 + (t1 * t244 + t14 * t602 + t15 * t601 - t17 * t361 + t2 * t243 - t72 * t83) * m(7) - (-Ifges(5,1) * t271 - t511 + t605) * t272 / 0.2e1 + (t91 * t424 + t607 / 0.2e1 - t418 * t536 - t422 * t540 - t420 * t541 + t178 * mrSges(5,2) + t581) * t271 + (Ifges(7,5) * t175 + Ifges(7,6) * t174) * t538 + t593 * t108 - t107 * t225 + (-qJ(5) * t78 - qJD(5) * t155 + t560 + t619) * t365 + t504 * t553 - t312 * t556 + t175 * t557 - t313 * t558 + t174 * t559 + t324 * t564 + t243 * t25 + t244 * t26 + t24 * t425; t125 * t97 - t621 * t96 - t449 * t154 + t220 * t155 + t16 + t62 + (t125 * t14 - t15 * t621 + t17 + t400) * m(7) + (t220 * t46 - t449 * t47 + t24 + t400) * m(6); -t72 * (mrSges(7,1) * t125 + mrSges(7,2) * t621) + (Ifges(7,1) * t621 - t515) * t550 + t52 * t549 + (Ifges(7,5) * t621 - Ifges(7,6) * t125) * t538 - t14 * t96 + t15 * t97 - g(1) * ((-t119 * t362 + t182 * t363) * mrSges(7,1) + (-t119 * t363 - t182 * t362) * mrSges(7,2)) - g(2) * ((-t117 * t362 + t180 * t363) * mrSges(7,1) + (-t117 * t363 - t180 * t362) * mrSges(7,2)) - g(3) * ((-t186 * t362 + t232 * t363) * mrSges(7,1) + (-t186 * t363 - t232 * t362) * mrSges(7,2)) + (t125 * t15 + t14 * t621) * mrSges(7,3) + t11 + (-Ifges(7,2) * t125 + t122 + t53) * t552 + t620;];
tau  = t3;
