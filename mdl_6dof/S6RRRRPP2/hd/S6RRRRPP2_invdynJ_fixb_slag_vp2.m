% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP2_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:48:58
% EndTime: 2019-03-09 20:49:47
% DurationCPUTime: 28.40s
% Computational Cost: add. (11631->807), mult. (25231->976), div. (0->0), fcn. (17390->10), ass. (0->361)
t677 = Ifges(7,4) + Ifges(6,5);
t681 = -Ifges(5,4) + t677;
t634 = Ifges(7,2) + Ifges(6,3);
t632 = Ifges(6,6) - Ifges(7,6);
t664 = Ifges(6,4) + Ifges(5,5);
t680 = -Ifges(7,5) + t664;
t361 = sin(qJ(2));
t365 = cos(qJ(2));
t477 = qJD(1) * qJD(2);
t299 = qJDD(1) * t365 - t361 * t477;
t300 = qJDD(1) * t361 + t365 * t477;
t360 = sin(qJ(3));
t364 = cos(qJ(3));
t292 = t360 * t361 - t364 * t365;
t386 = t292 * qJD(3);
t167 = -qJD(1) * t386 + t299 * t360 + t300 * t364;
t357 = qJDD(2) + qJDD(3);
t359 = sin(qJ(4));
t363 = cos(qJ(4));
t293 = t360 * t365 + t361 * t364;
t278 = t293 * qJD(1);
t475 = qJD(2) + qJD(3);
t242 = t359 * t278 - t363 * t475;
t482 = qJD(4) * t242;
t100 = t363 * t167 + t359 * t357 - t482;
t584 = t100 / 0.2e1;
t645 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t679 = t645 * t584;
t243 = t363 * t278 + t359 * t475;
t101 = qJD(4) * t243 + t359 * t167 - t363 * t357;
t582 = t101 / 0.2e1;
t387 = t293 * qJD(3);
t168 = -qJD(1) * t387 + t364 * t299 - t300 * t360;
t166 = qJDD(4) - t168;
t581 = -t166 / 0.2e1;
t580 = t166 / 0.2e1;
t678 = mrSges(6,2) + mrSges(5,3);
t635 = Ifges(6,2) + Ifges(5,3);
t663 = Ifges(5,6) - Ifges(6,6);
t676 = t677 * t243;
t480 = qJD(4) * t359;
t486 = qJD(1) * t365;
t487 = qJD(1) * t361;
t277 = -t360 * t487 + t364 * t486;
t514 = t277 * t359;
t649 = t480 - t514;
t639 = -m(6) - m(7);
t273 = qJD(4) - t277;
t659 = t242 * t634 + t273 * t632 + t676;
t546 = mrSges(6,2) * t242;
t171 = mrSges(6,3) * t273 - t546;
t544 = mrSges(5,3) * t242;
t173 = -mrSges(5,2) * t273 - t544;
t673 = -t173 - t171;
t543 = mrSges(5,3) * t243;
t175 = mrSges(5,1) * t273 - t543;
t545 = mrSges(6,2) * t243;
t176 = -mrSges(6,1) * t273 + t545;
t672 = t176 - t175;
t367 = -pkin(8) - pkin(7);
t311 = t367 * t365;
t296 = qJD(1) * t311;
t280 = t364 * t296;
t310 = t367 * t361;
t295 = qJD(1) * t310;
t222 = t295 * t360 - t280;
t484 = qJD(3) * t360;
t671 = pkin(2) * t484 - t222;
t670 = t677 * t242;
t669 = t677 * t363;
t668 = t677 * t359;
t667 = -mrSges(5,1) - mrSges(6,1) - mrSges(7,1);
t602 = -mrSges(7,2) - mrSges(6,3) + mrSges(5,2);
t368 = -pkin(4) - pkin(5);
t356 = t365 * pkin(2);
t345 = t356 + pkin(1);
t309 = t345 * qJD(1);
t181 = -pkin(3) * t277 - pkin(9) * t278 - t309;
t283 = qJD(2) * pkin(2) + t295;
t218 = t360 * t283 - t280;
t198 = pkin(9) * t475 + t218;
t479 = qJD(4) * t363;
t517 = qJDD(1) * pkin(1);
t268 = -pkin(2) * t299 - t517;
t53 = -pkin(3) * t168 - pkin(9) * t167 + t268;
t291 = t300 * pkin(7);
t234 = qJDD(2) * pkin(2) - pkin(8) * t300 - t291;
t290 = t299 * pkin(7);
t241 = pkin(8) * t299 + t290;
t483 = qJD(3) * t364;
t67 = t360 * t234 + t364 * t241 + t283 * t483 + t296 * t484;
t64 = pkin(9) * t357 + t67;
t14 = -t181 * t480 - t198 * t479 - t359 * t64 + t363 * t53;
t388 = qJDD(5) - t14;
t1 = -qJ(6) * t100 - qJD(6) * t243 + t166 * t368 + t388;
t666 = -mrSges(7,3) * t1 + t680 * t580 + t582 * t681 + t679;
t13 = t181 * t479 - t198 * t480 + t359 * t53 + t363 * t64;
t6 = t166 * qJ(5) + t273 * qJD(5) + t13;
t2 = qJ(6) * t101 + qJD(6) * t242 + t6;
t665 = t2 * mrSges(7,3) - Ifges(5,4) * t100 / 0.2e1 + Ifges(5,6) * t581 + t632 * t580 + t677 * t584 + (Ifges(5,2) + t634) * t582;
t661 = t475 * Ifges(4,5);
t660 = t475 * Ifges(4,6);
t658 = -t242 * t663 + t243 * t664 + t273 * t635;
t538 = Ifges(4,4) * t278;
t657 = t243 * Ifges(7,5) + Ifges(4,2) * t277 + t242 * Ifges(7,6) - t273 * Ifges(7,3) + t538 + t660;
t525 = t278 * mrSges(4,3);
t614 = mrSges(4,1) * t475 - mrSges(5,1) * t242 - mrSges(5,2) * t243 - t525;
t271 = pkin(4) * t480 - qJ(5) * t479 - t359 * qJD(5);
t656 = -pkin(5) * t649 - t271;
t358 = qJ(2) + qJ(3);
t351 = sin(t358);
t655 = t678 * t351;
t338 = qJ(6) * t480;
t654 = -qJ(6) * t514 + t338;
t653 = t359 * t634 + t669;
t652 = -t359 * t663 + t363 * t664;
t513 = t277 * t363;
t491 = pkin(4) * t514 - qJ(5) * t513;
t650 = -t491 + t671;
t648 = -t513 + t479;
t647 = t13 * t363 - t14 * t359;
t8 = -pkin(4) * t166 + t388;
t646 = t359 * t8 + t363 * t6;
t240 = Ifges(5,4) * t242;
t597 = t243 * t645 + t273 * t680 - t240 + t670;
t536 = Ifges(5,4) * t359;
t643 = t363 * t645 - t536 + t668;
t353 = t359 * qJ(5);
t415 = -t363 * mrSges(6,1) - t359 * mrSges(6,3);
t424 = m(7) * t368 - mrSges(7,1);
t504 = t351 * t363;
t523 = t359 * mrSges(7,2);
t401 = t363 * pkin(4) + t353;
t631 = -pkin(3) - t401;
t642 = mrSges(5,1) * t504 + (-m(7) * (-pkin(3) - t353) + t523 - t424 * t363 - m(6) * t631 - t415) * t351;
t39 = -mrSges(6,2) * t101 + mrSges(6,3) * t166;
t41 = mrSges(5,1) * t166 - mrSges(5,3) * t100;
t42 = -t166 * mrSges(6,1) + t100 * mrSges(6,2);
t44 = -mrSges(5,2) * t166 - mrSges(5,3) * t101;
t92 = t363 * t181 - t359 * t198;
t59 = -pkin(4) * t273 + qJD(5) - t92;
t257 = t273 * qJ(5);
t93 = t359 * t181 + t363 * t198;
t62 = t257 + t93;
t641 = t363 * (t44 + t39) + t359 * (-t41 + t42) + m(5) * ((-t359 * t93 - t363 * t92) * qJD(4) + t647) + m(6) * ((-t359 * t62 + t363 * t59) * qJD(4) + t646) + t673 * t480 + t672 * t479;
t279 = t360 * t296;
t217 = t364 * t283 + t279;
t197 = -pkin(3) * t475 - t217;
t522 = t363 * mrSges(7,2);
t413 = -t359 * mrSges(7,1) + t522;
t521 = t363 * mrSges(6,3);
t414 = mrSges(6,1) * t359 - t521;
t416 = mrSges(5,1) * t359 + mrSges(5,2) * t363;
t380 = -t243 * qJ(5) + t197;
t49 = t242 * t368 + qJD(6) - t380;
t84 = t242 * pkin(4) + t380;
t640 = t197 * t416 + t413 * t49 + t414 * t84;
t638 = t299 / 0.2e1;
t563 = t365 / 0.2e1;
t627 = t365 * Ifges(3,2);
t560 = pkin(2) * t360;
t343 = pkin(9) + t560;
t464 = pkin(2) * t483;
t420 = -qJD(6) + t464;
t213 = pkin(3) * t278 - pkin(9) * t277;
t467 = pkin(2) * t487;
t191 = t213 + t467;
t223 = t295 * t364 + t279;
t118 = t359 * t191 + t363 * t223;
t263 = t278 * qJ(5);
t76 = t263 + t118;
t624 = -t343 * t480 + t363 * t420 + t654 - t76;
t214 = t359 * t223;
t492 = -qJ(6) + t343;
t289 = t492 * t363;
t461 = t368 * t278;
t518 = qJ(6) * t277;
t623 = qJD(4) * t289 + t359 * t420 - t214 - (-t191 - t518) * t363 - t461;
t622 = -t650 + t656;
t134 = t218 + t491;
t621 = t134 + t656;
t124 = t359 * t213 + t363 * t217;
t85 = t263 + t124;
t620 = -pkin(9) * t480 - qJD(6) * t363 + t654 - t85;
t209 = t359 * t217;
t547 = pkin(9) - qJ(6);
t308 = t547 * t363;
t619 = qJD(4) * t308 - qJD(6) * t359 - t209 - (-t213 - t518) * t363 - t461;
t51 = qJ(6) * t243 + t92;
t617 = qJD(5) - t51;
t542 = mrSges(7,3) * t242;
t172 = mrSges(7,2) * t273 + t542;
t616 = t171 + t172;
t231 = qJD(2) * t293 + t387;
t615 = t231 * qJ(5) + t292 * qJD(5);
t613 = t271 + t650;
t612 = t271 - t134;
t611 = t364 * t310 + t311 * t360;
t352 = cos(t358);
t610 = t352 * mrSges(4,1) - t351 * mrSges(4,2);
t607 = t100 * t664 - t101 * t663 + t166 * t635;
t554 = pkin(7) * t365;
t555 = pkin(7) * t361;
t606 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t487) * t554 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t486) * t555;
t605 = t290 * t365 + t291 * t361;
t362 = sin(qJ(1));
t366 = cos(qJ(1));
t604 = g(1) * t366 + g(2) * t362;
t599 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t307 = -mrSges(3,1) * t365 + mrSges(3,2) * t361;
t598 = -m(3) * pkin(1) - mrSges(2,1) + t307 - t610;
t557 = pkin(3) * t351;
t559 = pkin(2) * t361;
t596 = -m(7) * (-qJ(6) * t352 - t559) + t352 * mrSges(7,3) + m(6) * t559 - m(5) * (-t557 - t559) + t642;
t593 = -m(7) * qJ(6) - mrSges(7,3);
t595 = -t352 * t593 + t642;
t499 = t352 * t366;
t328 = pkin(9) * t499;
t505 = t351 * t359;
t468 = mrSges(5,2) * t505;
t594 = t328 * t639 - t366 * t468 - t499 * t678;
t245 = t310 * t360 - t311 * t364;
t457 = qJD(2) * t367;
t297 = t361 * t457;
t298 = t365 * t457;
t147 = t245 * qJD(3) + t297 * t360 - t364 * t298;
t500 = t352 * t363;
t590 = t352 * t359 * t602 + t351 * mrSges(7,3) + t500 * t667 - t610 - t655;
t589 = m(7) * pkin(5) - t667;
t588 = (-t468 + (-t678 + (-m(5) + t639) * pkin(9)) * t352) * t362;
t586 = t14 * mrSges(5,1) - t8 * mrSges(6,1) - t1 * mrSges(7,1) - t13 * mrSges(5,2) + t2 * mrSges(7,2) + t6 * mrSges(6,3);
t583 = -t101 / 0.2e1;
t577 = -t242 / 0.2e1;
t576 = t242 / 0.2e1;
t575 = -t243 / 0.2e1;
t574 = t243 / 0.2e1;
t573 = -t273 / 0.2e1;
t572 = t273 / 0.2e1;
t571 = -t277 / 0.2e1;
t570 = t277 / 0.2e1;
t568 = t278 / 0.2e1;
t558 = pkin(2) * t364;
t556 = pkin(4) * t278;
t336 = t351 * pkin(9);
t337 = t352 * pkin(3);
t541 = mrSges(7,3) * t243;
t540 = Ifges(3,4) * t361;
t539 = Ifges(3,4) * t365;
t537 = Ifges(5,4) * t243;
t535 = Ifges(5,4) * t363;
t526 = t277 * mrSges(4,3);
t520 = qJ(5) * t242;
t519 = qJ(5) * t363;
t230 = -qJD(2) * t292 - t386;
t516 = t230 * t359;
t515 = t230 * t363;
t510 = t293 * t359;
t503 = t351 * t366;
t130 = -Ifges(5,2) * t242 + Ifges(5,6) * t273 + t537;
t498 = t359 * t130;
t497 = t359 * t362;
t496 = t362 * t363;
t495 = t363 * t366;
t494 = t366 * t359;
t493 = t366 * t367;
t216 = pkin(3) * t292 - pkin(9) * t293 - t345;
t139 = t359 * t216 + t363 * t245;
t488 = t337 + t336;
t485 = qJD(2) * t361;
t478 = qJD(5) * t363;
t466 = pkin(2) * t485;
t137 = pkin(3) * t231 - pkin(9) * t230 + t466;
t146 = qJD(3) * t611 + t297 * t364 + t298 * t360;
t460 = t359 * t137 + t363 * t146 + t216 * t479;
t459 = t359 * t146 + t216 * t480 + t245 * t479;
t109 = t292 * qJ(5) + t139;
t344 = -pkin(3) - t558;
t456 = t293 * t479;
t440 = -t480 / 0.2e1;
t438 = t479 / 0.2e1;
t35 = -t101 * mrSges(7,1) + t100 * mrSges(7,2);
t435 = -t345 - t337;
t434 = t477 / 0.2e1;
t40 = -t166 * mrSges(7,1) - t100 * mrSges(7,3);
t117 = t191 * t363 - t214;
t123 = t213 * t363 - t209;
t235 = t359 * t245;
t138 = t216 * t363 - t235;
t428 = t366 * t345 - t362 * t367;
t427 = t359 * t464;
t426 = t363 * t464;
t68 = t364 * t234 - t360 * t241 - t283 * t484 + t296 * t483;
t425 = pkin(4) * t500 + t352 * t353 + t488;
t52 = qJ(6) * t242 + t93;
t419 = mrSges(3,1) * t361 + mrSges(3,2) * t365;
t417 = mrSges(4,1) * t351 + mrSges(4,2) * t352;
t409 = t540 + t627;
t408 = -Ifges(5,2) * t359 + t535;
t405 = Ifges(3,5) * t365 - Ifges(3,6) * t361;
t402 = Ifges(7,5) * t363 + Ifges(7,6) * t359;
t400 = pkin(4) * t359 - t519;
t285 = t344 - t401;
t399 = -qJ(6) * t230 - qJD(6) * t293;
t32 = t137 * t363 - t459;
t396 = pkin(3) * t499 + pkin(9) * t503 + t428;
t395 = pkin(1) * t419;
t394 = t359 * t368 + t519;
t392 = t456 + t516;
t391 = t293 * t480 - t515;
t390 = pkin(3) * t357 + t68;
t389 = t361 * (Ifges(3,1) * t365 - t540);
t31 = -t245 * t480 + t460;
t384 = Ifges(7,5) * t100 + Ifges(7,6) * t101 - Ifges(7,3) * t166;
t383 = pkin(5) * t500 - qJ(6) * t351 + t425;
t264 = t352 * t497 + t495;
t266 = t352 * t494 - t496;
t381 = -g(1) * t266 - g(2) * t264 - g(3) * t505;
t377 = qJ(5) * t100 + qJD(5) * t243 + t390;
t15 = pkin(4) * t101 - t377;
t269 = Ifges(4,4) * t277;
t193 = Ifges(4,1) * t278 + t269 + t661;
t37 = t273 * t368 + t617;
t4 = t101 * t368 + qJDD(6) + t377;
t46 = t257 + t52;
t369 = (t643 * t277 + t278 * t680) * t575 + (t535 - t669) * t584 + (t309 * mrSges(4,2) + t402 * t572 + t408 * t576 - t661 / 0.2e1 + t653 * t577 + t652 * t573 - t640) * t277 + (t390 * mrSges(5,1) + t4 * mrSges(7,1) + Ifges(5,2) * t583 - Ifges(7,6) * t581 + t663 * t580 - t634 * t582 - t665) * t363 + t536 * t583 + (-t46 * mrSges(7,2) + t93 * mrSges(5,2) + t37 * mrSges(7,1) - t92 * mrSges(5,1) + t309 * mrSges(4,1) - Ifges(4,2) * t571 - Ifges(7,3) * t572 + Ifges(5,6) * t576 + t660 / 0.2e1 + t59 * mrSges(6,1) - t62 * mrSges(6,3) + t632 * t577 + t635 * t573) * t278 + (-t390 * mrSges(5,2) + Ifges(7,5) * t581 + t664 * t580 + t666 + t679) * t359 + (-t37 * t648 + t46 * t649) * mrSges(7,3) + (t269 + t193) * t571 + ((t652 / 0.2e1 - t402 / 0.2e1) * t273 + t574 * t643 + t640) * qJD(4) + t668 * t582 + t657 * t568 - (Ifges(4,1) * t277 - t538 + t658) * t278 / 0.2e1 + (t59 * t648 - t62 * t649 + t646) * mrSges(6,2) + (-t648 * t92 - t649 * t93 + t647) * mrSges(5,3) + t498 * t570 + t130 * t440 + t659 * (-t514 / 0.2e1 + t480 / 0.2e1) + (-t513 / 0.2e1 + t438) * t597 + Ifges(4,3) * t357 + t218 * t525 + t217 * t526 + Ifges(4,5) * t167 + Ifges(4,6) * t168 + (t653 / 0.2e1 - t408 / 0.2e1) * t482 - t67 * mrSges(4,2) + t68 * mrSges(4,1) + t15 * t415 + t4 * t523;
t354 = t363 * pkin(5);
t347 = Ifges(3,4) * t486;
t306 = t547 * t359;
t288 = t492 * t359;
t286 = t354 - t631;
t276 = Ifges(3,1) * t487 + Ifges(3,5) * qJD(2) + t347;
t275 = Ifges(3,6) * qJD(2) + qJD(1) * t409;
t267 = t352 * t495 + t497;
t265 = t352 * t496 - t494;
t256 = -t285 + t354;
t246 = -mrSges(4,2) * t475 + t526;
t212 = -mrSges(4,1) * t277 + mrSges(4,2) * t278;
t174 = -mrSges(7,1) * t273 - t541;
t155 = -mrSges(4,2) * t357 + mrSges(4,3) * t168;
t154 = mrSges(4,1) * t357 - mrSges(4,3) * t167;
t152 = -mrSges(7,1) * t242 + mrSges(7,2) * t243;
t151 = mrSges(6,1) * t242 - mrSges(6,3) * t243;
t150 = pkin(4) * t243 + t520;
t145 = t293 * t400 - t611;
t119 = t293 * t394 + t611;
t116 = -pkin(4) * t292 - t138;
t108 = t243 * t368 - t520;
t87 = -t123 - t556;
t77 = -t117 - t556;
t75 = qJ(6) * t510 + t109;
t58 = t235 + (-qJ(6) * t293 - t216) * t363 + t368 * t292;
t43 = mrSges(7,2) * t166 + mrSges(7,3) * t101;
t36 = mrSges(5,1) * t101 + mrSges(5,2) * t100;
t34 = mrSges(6,1) * t101 - mrSges(6,3) * t100;
t33 = t400 * t230 + (qJD(4) * t401 - t478) * t293 + t147;
t24 = t394 * t230 + (t478 + (t363 * t368 - t353) * qJD(4)) * t293 - t147;
t17 = -pkin(4) * t231 - t32;
t16 = t31 + t615;
t12 = qJ(6) * t456 + (-qJD(4) * t245 - t399) * t359 + t460 + t615;
t9 = t293 * t338 + t368 * t231 + (-t137 + t399) * t363 + t459;
t3 = [(t607 / 0.2e1 - t384 / 0.2e1 - Ifges(4,6) * t357 - Ifges(4,2) * t168 - Ifges(4,4) * t167 + t268 * mrSges(4,1) - t67 * mrSges(4,3) - Ifges(7,3) * t581 + Ifges(5,6) * t583 + t632 * t582 + t635 * t580 + t680 * t584 + t586) * t292 + (t231 * t632 - t391 * t677 + t392 * t634) * t576 + (-m(4) * t428 - m(5) * t396 + t639 * (t267 * pkin(4) + qJ(5) * t266 + t396) + (-t593 - t678) * t503 + t598 * t366 + t599 * t362 - t589 * t267 + t602 * t266) * g(2) + (t659 * t438 + t268 * mrSges(4,2) - t68 * mrSges(4,3) + Ifges(4,1) * t167 + Ifges(4,4) * t168 + Ifges(4,5) * t357 + t15 * t414 - t416 * t390 + t4 * t413 + t402 * t581 + t408 * t583 + t597 * t440 + t580 * t652 + t582 * t653 + t584 * t643 + (mrSges(6,2) * t8 - mrSges(5,3) * t14 + t666) * t363) * t293 + t659 * t516 / 0.2e1 + (-t6 * mrSges(6,2) - t13 * mrSges(5,3) + t665) * t510 + m(4) * (t146 * t218 + t245 * t67 - t268 * t345 - t309 * t466) + t409 * t638 + m(5) * (t13 * t139 + t138 * t14 + t31 * t93 + t32 * t92) + (t365 * t539 + t389) * t434 + (t231 * t635 - t391 * t664 - t392 * t663) * t572 + (m(5) * t493 + t639 * (-t265 * pkin(4) - qJ(5) * t264 - t493) + (m(4) * t367 + t599) * t366 + t589 * t265 - t602 * t264 + (-m(7) * t435 - (-m(7) * t547 + mrSges(7,3)) * t351 + m(4) * t345 + (-m(5) - m(6)) * (t435 - t336) - t598 + t655) * t362) * g(1) - t657 * t231 / 0.2e1 + t658 * t231 / 0.2e1 - t130 * t456 / 0.2e1 + (Ifges(3,4) * t300 + Ifges(3,2) * t299) * t563 + (Ifges(4,1) * t230 - Ifges(4,4) * t231) * t568 + (Ifges(4,4) * t230 - Ifges(4,2) * t231) * t570 + (-Ifges(7,5) * t391 + Ifges(7,6) * t392 - Ifges(7,3) * t231) * t573 + (-Ifges(5,4) * t391 - Ifges(5,2) * t392 + Ifges(5,6) * t231) * t577 + t475 * (Ifges(4,5) * t230 - Ifges(4,6) * t231) / 0.2e1 - t395 * t477 - t275 * t485 / 0.2e1 - t307 * t517 + t212 * t466 + m(7) * (t1 * t58 + t119 * t4 + t12 * t46 + t2 * t75 + t24 * t49 + t37 * t9) + m(6) * (t109 * t6 + t116 * t8 + t145 * t15 + t16 * t62 + t17 * t59 + t33 * t84) - (-m(4) * t68 - m(5) * t390 - t154 + t36) * t611 + Ifges(2,3) * qJDD(1) - t345 * (-mrSges(4,1) * t168 + mrSges(4,2) * t167) - t309 * (mrSges(4,1) * t231 + mrSges(4,2) * t230) - pkin(1) * (-mrSges(3,1) * t299 + mrSges(3,2) * t300) + t245 * t155 + t146 * t246 + t230 * t193 / 0.2e1 + t31 * t173 + t9 * t174 + t32 * t175 + t17 * t176 + t16 * t171 + t12 * t172 + t33 * t151 + t24 * t152 + t139 * t44 + t145 * t34 + t138 * t41 + t116 * t42 + t119 * t35 + t109 * t39 + t75 * t43 + t58 * t40 + t597 * t515 / 0.2e1 + (-t217 * t230 - t218 * t231) * mrSges(4,3) - t230 * t498 / 0.2e1 + (t680 * t231 - t645 * t391 + t392 * t681) * t574 + (t299 * t554 + t300 * t555 + t605) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t605) + (t276 * t563 + t405 * qJD(2) / 0.2e1 - t606) * qJD(2) + (-m(4) * t217 + m(5) * t197 - t614) * t147 + (-mrSges(3,1) * t555 - mrSges(3,2) * t554 + 0.2e1 * Ifges(3,6) * t563) * qJDD(2) + (Ifges(3,1) * t300 + Ifges(3,4) * t638 + Ifges(3,5) * qJDD(2) - t434 * t627) * t361 + t37 * (-mrSges(7,1) * t231 + mrSges(7,3) * t391) + t92 * (mrSges(5,1) * t231 + mrSges(5,3) * t391) + t59 * (-mrSges(6,1) * t231 - mrSges(6,2) * t391) + t300 * t539 / 0.2e1 + t46 * (mrSges(7,2) * t231 + mrSges(7,3) * t392) + t84 * (mrSges(6,1) * t392 + mrSges(6,3) * t391) + t197 * (mrSges(5,1) * t392 - mrSges(5,2) * t391) + t62 * (-mrSges(6,2) * t392 + mrSges(6,3) * t231) + t93 * (-mrSges(5,2) * t231 - mrSges(5,3) * t392) + t49 * (-mrSges(7,1) * t392 - mrSges(7,2) * t391); ((t360 * t67 + t364 * t68 + (-t217 * t360 + t218 * t364) * qJD(3)) * pkin(2) + t217 * t222 - t218 * t223 + t309 * t467) * m(4) + (-t344 * t390 + (t197 * t360 + (-t359 * t92 + t363 * t93) * t364) * qJD(3) * pkin(2) - t117 * t92 - t118 * t93 - t197 * t222) * m(5) - t614 * t671 + (-m(4) * t356 - m(6) * (t356 + t425) + t307 - m(5) * (t356 + t488) - m(7) * (t356 + t383) + t590) * g(3) - (-Ifges(3,2) * t487 + t276 + t347) * t486 / 0.2e1 + t641 * t343 + (-t77 + t427) * t176 + (-t427 - t117) * t175 + t154 * t558 + t155 * t560 + (t606 + (t395 - t389 / 0.2e1) * qJD(1)) * qJD(1) + (m(4) * t559 + t417 + t419) * t604 + (-t118 + t426) * t173 + (t464 - t223) * t246 + t369 - t405 * t477 / 0.2e1 + t275 * t487 / 0.2e1 + (-t76 + t426) * t171 - t212 * t467 + t344 * t36 + Ifges(3,6) * t299 + Ifges(3,5) * t300 - t290 * mrSges(3,2) - t291 * mrSges(3,1) + t285 * t34 + t288 * t40 + t289 * t43 + t256 * t35 + Ifges(3,3) * qJDD(2) + (t362 * t596 + t588) * g(2) + (-m(5) * t328 + t366 * t596 + t594) * g(1) + (t15 * t285 + (t359 * t59 + t363 * t62) * t464 - t59 * t77 - t62 * t76 + t613 * t84) * m(6) + t613 * t151 + t622 * t152 + t623 * t174 + t624 * t172 + (t1 * t288 + t2 * t289 + t256 * t4 + t37 * t623 + t46 * t624 + t49 * t622) * m(7); (-m(5) * t488 - m(6) * t425 - m(7) * t383 + t590) * g(3) + (pkin(3) * t390 - t123 * t92 - t124 * t93 - t197 * t218) * m(5) + t369 + t631 * t34 + t306 * t40 + t308 * t43 + t286 * t35 - t217 * t246 - t124 * t173 - t123 * t175 - t87 * t176 - t85 * t171 - pkin(3) * t36 + ((m(5) * t557 + t595) * t362 + t588) * g(2) + (-m(5) * (-pkin(3) * t503 + t328) + t595 * t366 + t594) * g(1) + t604 * t417 + (t15 * t631 - t59 * t87 + t612 * t84 - t62 * t85) * m(6) + t612 * t151 + t614 * t218 + t619 * t174 + t620 * t172 + t621 * t152 + (t1 * t306 + t2 * t308 + t286 * t4 + t37 * t619 + t46 * t620 + t49 * t621) * m(7) + t641 * pkin(9); (t243 * t634 - t670) * t577 + (-t242 * t645 - t537 + t659 + t676) * t575 + (-t242 * t664 - t243 * t663) * t573 + t586 + (-Ifges(7,5) * t242 + Ifges(7,6) * t243) * t572 + t130 * t574 + t62 * t545 + t59 * t546 - t46 * t541 - t37 * t542 + (-pkin(4) * t8 + qJ(5) * t6 + qJD(5) * t62 - t150 * t84) * m(6) + (-m(6) * t59 + t543 - t672) * t93 + (-m(6) * t62 - t544 + t673) * t92 - t384 + t607 + t368 * t40 - t197 * (mrSges(5,1) * t243 - mrSges(5,2) * t242) - t84 * (mrSges(6,1) * t243 + mrSges(6,3) * t242) - t49 * (-mrSges(7,1) * t243 - mrSges(7,2) * t242) + (t39 + t43) * qJ(5) - t52 * t174 - t51 * t172 - t150 * t151 - t108 * t152 - pkin(4) * t42 + (-Ifges(5,2) * t243 - t240 + t597) * t576 + t616 * qJD(5) + (t2 * qJ(5) + t1 * t368 - t108 * t49 - t37 * t52 + t46 * t617) * m(7) + (t639 * qJ(5) * t504 + (-t521 - t522 + t416 + (m(6) * pkin(4) + mrSges(6,1) - t424) * t359) * t351) * g(3) + (t639 * (-t266 * pkin(4) + qJ(5) * t267) + t602 * t267 + t589 * t266) * g(1) + (t639 * (-t264 * pkin(4) + qJ(5) * t265) + t602 * t265 + t589 * t264) * g(2); -t616 * t273 + (t151 - t152) * t243 + t40 + t42 + (-t243 * t49 - t273 * t46 + t1 + t381) * m(7) + (t243 * t84 - t273 * t62 + t381 + t8) * m(6); -t242 * t172 + t243 * t174 + (-g(3) * t352 - t242 * t46 + t37 * t243 + t351 * t604 + t4) * m(7) + t35;];
tau  = t3;
