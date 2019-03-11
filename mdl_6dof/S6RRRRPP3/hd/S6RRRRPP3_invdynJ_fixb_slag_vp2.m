% Calculate vector of inverse dynamics joint torques for
% S6RRRRPP3
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRPP3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:50
% EndTime: 2019-03-09 20:54:37
% DurationCPUTime: 27.18s
% Computational Cost: add. (11649->807), mult. (25266->978), div. (0->0), fcn. (17402->10), ass. (0->370)
t672 = Ifges(5,1) + Ifges(7,3);
t636 = Ifges(7,4) + Ifges(6,5);
t635 = Ifges(5,5) + Ifges(7,5);
t671 = -Ifges(6,3) - Ifges(7,2);
t604 = Ifges(6,4) - t635;
t682 = Ifges(5,6) - t636;
t366 = sin(qJ(3));
t367 = sin(qJ(2));
t370 = cos(qJ(3));
t371 = cos(qJ(2));
t289 = t366 * t371 + t367 * t370;
t276 = t289 * qJD(1);
t365 = sin(qJ(4));
t369 = cos(qJ(4));
t486 = qJD(1) * qJD(2);
t297 = qJDD(1) * t371 - t367 * t486;
t298 = qJDD(1) * t367 + t371 * t486;
t288 = t366 * t367 - t370 * t371;
t392 = t288 * qJD(3);
t380 = -qJD(1) * t392 + t297 * t366 + t298 * t370;
t484 = qJD(2) + qJD(3);
t439 = t369 * t484;
t483 = qJDD(2) + qJDD(3);
t488 = qJD(4) * t365;
t103 = -qJD(4) * t439 + t276 * t488 - t365 * t483 - t369 * t380;
t589 = -t103 / 0.2e1;
t588 = t103 / 0.2e1;
t241 = t369 * t276 + t365 * t484;
t489 = qJD(4) * t241;
t104 = t365 * t380 - t369 * t483 + t489;
t587 = -t104 / 0.2e1;
t586 = t104 / 0.2e1;
t393 = t289 * qJD(3);
t172 = -qJD(1) * t393 + t297 * t370 - t298 * t366;
t171 = qJDD(4) - t172;
t585 = t171 / 0.2e1;
t688 = mrSges(6,1) + mrSges(5,3);
t687 = -Ifges(5,4) + Ifges(7,6);
t686 = Ifges(6,6) - Ifges(7,6);
t487 = qJD(4) * t369;
t494 = qJD(1) * t371;
t495 = qJD(1) * t367;
t275 = -t366 * t495 + t370 * t494;
t521 = t275 * t369;
t653 = t487 - t521;
t522 = t275 * t365;
t652 = t488 - t522;
t605 = Ifges(6,1) + Ifges(7,1) + Ifges(5,3);
t641 = -m(6) - m(7);
t683 = -t171 / 0.2e1;
t240 = t276 * t365 - t439;
t237 = Ifges(5,4) * t240;
t270 = qJD(4) - t275;
t539 = Ifges(7,6) * t240;
t664 = t241 * t672 + t635 * t270 - t237 + t539;
t235 = Ifges(7,6) * t241;
t533 = t241 * Ifges(6,6);
t663 = -t240 * t671 + t270 * t636 + t235 - t533;
t553 = mrSges(6,1) * t240;
t176 = -mrSges(6,3) * t270 + t553;
t548 = mrSges(5,3) * t240;
t179 = -mrSges(5,2) * t270 - t548;
t681 = -t179 + t176;
t552 = mrSges(6,1) * t241;
t178 = mrSges(6,2) * t270 + t552;
t547 = mrSges(5,3) * t241;
t180 = mrSges(5,1) * t270 - t547;
t680 = -t180 + t178;
t373 = -pkin(8) - pkin(7);
t309 = t373 * t367;
t291 = qJD(1) * t309;
t311 = t373 * t371;
t292 = qJD(1) * t311;
t501 = t370 * t292;
t222 = t291 * t366 - t501;
t492 = qJD(3) * t366;
t679 = pkin(2) * t492 - t222;
t678 = -m(5) + t641;
t282 = qJD(2) * pkin(2) + t291;
t219 = t366 * t282 - t501;
t198 = pkin(9) * t484 + t219;
t525 = qJDD(1) * pkin(1);
t267 = -t297 * pkin(2) - t525;
t54 = -t172 * pkin(3) - pkin(9) * t380 + t267;
t362 = t371 * pkin(2);
t351 = t362 + pkin(1);
t307 = t351 * qJD(1);
t386 = -pkin(3) * t275 - pkin(9) * t276 - t307;
t287 = t298 * pkin(7);
t232 = qJDD(2) * pkin(2) - pkin(8) * t298 - t287;
t286 = t297 * pkin(7);
t239 = pkin(8) * t297 + t286;
t491 = qJD(3) * t370;
t69 = t366 * t232 + t370 * t239 + t282 * t491 + t292 * t492;
t666 = pkin(9) * t483 + qJD(4) * t386 + t69;
t13 = -t198 * t488 + t365 * t54 + t369 * t666;
t7 = -qJ(5) * t171 - qJD(5) * t270 - t13;
t4 = -pkin(5) * t104 + qJDD(6) - t7;
t677 = -t4 * mrSges(7,1) + Ifges(5,2) * t586 + Ifges(5,6) * t683 + t585 * t636 + t587 * t671 + (Ifges(5,4) + t686) * t588;
t14 = -t198 * t487 - t365 * t666 + t369 * t54;
t394 = qJDD(5) - t14;
t554 = pkin(4) + qJ(6);
t2 = -pkin(5) * t103 - qJD(6) * t270 - t171 * t554 + t394;
t676 = mrSges(7,1) * t2 + Ifges(6,4) * t683 + Ifges(6,6) * t587 + t635 * t585 + t687 * t586 + (Ifges(6,2) + t672) * t589;
t675 = -mrSges(5,1) + mrSges(6,2) - mrSges(7,3);
t606 = mrSges(7,2) - mrSges(5,2) + mrSges(6,3);
t363 = qJ(2) + qJ(3);
t358 = sin(t363);
t447 = -qJ(5) * t365 - pkin(3);
t303 = -pkin(4) * t369 + t447;
t430 = -m(7) * t554 - mrSges(7,3);
t529 = t365 * mrSges(6,3);
t530 = t365 * mrSges(7,2);
t613 = (-m(6) * t303 - m(7) * t447 - t369 * t430 + t529 + t530) * t358;
t674 = m(5) * pkin(3) * t358 + t613;
t673 = -m(7) * pkin(5) - mrSges(7,1) - t688;
t668 = t484 * Ifges(4,5);
t667 = t484 * Ifges(4,6);
t213 = pkin(3) * t276 - pkin(9) * t275;
t479 = pkin(2) * t495;
t191 = t213 + t479;
t277 = t366 * t292;
t223 = t291 * t370 + t277;
t121 = t365 * t191 + t369 * t223;
t477 = pkin(2) * t491;
t437 = t369 * t477;
t665 = -t121 + t437;
t440 = pkin(4) * t488 - qJD(5) * t365;
t662 = (-qJ(5) * qJD(4) - qJD(6)) * t369 + t440 + t652 * qJ(6);
t531 = t276 * mrSges(4,3);
t619 = mrSges(4,1) * t484 - mrSges(5,1) * t240 - mrSges(5,2) * t241 - t531;
t661 = t688 * t358;
t660 = pkin(5) * t653 + t276 * t554;
t429 = pkin(4) * t522 - qJ(5) * t521;
t659 = -t429 + t679;
t538 = Ifges(7,6) * t365;
t543 = Ifges(5,4) * t365;
t658 = t369 * t672 + t538 - t543;
t537 = Ifges(7,6) * t369;
t540 = Ifges(6,6) * t369;
t657 = -t365 * t671 + t537 - t540;
t120 = t191 * t369 - t365 * t223;
t438 = t365 * t477;
t656 = -t438 - t120;
t228 = -qJD(2) * t288 - t392;
t524 = t228 * t365;
t397 = t289 * t487 + t524;
t236 = Ifges(6,6) * t240;
t132 = Ifges(6,4) * t270 - Ifges(6,2) * t241 + t236;
t534 = t241 * Ifges(5,4);
t133 = -t240 * Ifges(5,2) + t270 * Ifges(5,6) + t534;
t651 = t369 * t132 + t365 * t133;
t650 = t13 * t369 - t14 * t365;
t9 = -pkin(4) * t171 + t394;
t649 = t365 * t9 - t369 * t7;
t94 = t198 * t365 - t369 * t386;
t401 = pkin(5) * t241 + t94;
t647 = qJD(5) + t401;
t646 = t240 * pkin(5) - qJD(6);
t644 = -t240 * t682 - t241 * t604 + t270 * t605;
t643 = -t365 * t682 - t369 * t604;
t218 = t370 * t282 + t277;
t197 = -pkin(3) * t484 - t218;
t528 = t369 * mrSges(7,2);
t419 = mrSges(7,3) * t365 - t528;
t527 = t369 * mrSges(6,3);
t420 = -mrSges(6,2) * t365 - t527;
t421 = mrSges(5,1) * t365 + mrSges(5,2) * t369;
t385 = -t241 * qJ(5) + t197;
t50 = t240 * t554 + t385;
t86 = t240 * pkin(4) + t385;
t642 = t197 * t421 + t419 * t50 + t420 * t86;
t638 = t297 / 0.2e1;
t569 = t371 / 0.2e1;
t41 = mrSges(6,1) * t104 - mrSges(6,3) * t171;
t45 = -mrSges(5,2) * t171 - mrSges(5,3) * t104;
t633 = -t41 + t45;
t43 = -t103 * mrSges(6,1) + t171 * mrSges(6,2);
t44 = mrSges(5,1) * t171 + mrSges(5,3) * t103;
t632 = t43 - t44;
t631 = t371 * Ifges(3,2);
t628 = t659 + t662;
t137 = t429 + t219;
t627 = -t137 + t662;
t568 = pkin(2) * t366;
t349 = pkin(9) + t568;
t262 = t276 * qJ(5);
t433 = -pkin(5) * t522 + t262;
t626 = (-pkin(5) - t349) * t488 - t433 + t665;
t625 = t349 * t487 - t656 + t660;
t127 = t365 * t213 + t369 * t218;
t590 = -pkin(5) - pkin(9);
t624 = t590 * t488 - t127 - t433;
t126 = t213 * t369 - t365 * t218;
t623 = pkin(9) * t487 + t126 + t660;
t95 = t369 * t198 + t365 * t386;
t621 = -t95 + t646;
t551 = mrSges(7,1) * t240;
t177 = mrSges(7,2) * t270 - t551;
t620 = t177 - t176;
t269 = -qJ(5) * t487 + t440;
t618 = t269 + t659;
t617 = t269 - t137;
t616 = t370 * t309 + t311 * t366;
t359 = cos(t363);
t615 = t359 * mrSges(4,1) - mrSges(4,2) * t358;
t562 = pkin(7) * t371;
t563 = pkin(7) * t367;
t610 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t495) * t562 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t494) * t563;
t609 = t286 * t371 + t287 * t367;
t368 = sin(qJ(1));
t372 = cos(qJ(1));
t608 = g(1) * t372 + g(2) * t368;
t602 = -m(3) * pkin(7) + mrSges(2,2) - mrSges(3,3) - mrSges(4,3);
t306 = -mrSges(3,1) * t371 + mrSges(3,2) * t367;
t601 = -m(3) * pkin(1) - mrSges(2,1) + t306 - t615;
t512 = t358 * t369;
t567 = pkin(2) * t367;
t600 = mrSges(5,1) * t512 - t567 * t678 + t674;
t507 = t359 * t372;
t330 = pkin(9) * t507;
t502 = t369 * t372;
t471 = t358 * t502;
t513 = t358 * t365;
t480 = mrSges(5,2) * t513;
t599 = -mrSges(6,2) * t471 + t330 * t641 - t372 * t480 + t673 * t507;
t598 = t103 * t604 - t104 * t682 + t605 * t171;
t508 = t359 * t369;
t510 = t359 * t365;
t597 = -t358 * mrSges(7,1) + t508 * t675 - t510 * t606 - t615 - t661;
t596 = m(7) * qJ(6) - t675;
t504 = t368 * t369;
t472 = t358 * t504;
t595 = -mrSges(6,2) * t472 + (-t480 + (pkin(9) * t678 + t673) * t359) * t368;
t61 = -pkin(4) * t270 + qJD(5) + t94;
t64 = -t270 * qJ(5) - t95;
t594 = m(5) * ((-t365 * t95 + t369 * t94) * qJD(4) + t650) + m(6) * ((t365 * t64 + t369 * t61) * qJD(4) + t649) + t681 * t488 + t680 * t487;
t593 = t14 * mrSges(5,1) - t13 * mrSges(5,2) + t9 * mrSges(6,2) + t4 * mrSges(7,2) - t7 * mrSges(6,3) - t2 * mrSges(7,3);
t582 = -t240 / 0.2e1;
t581 = t240 / 0.2e1;
t580 = -t241 / 0.2e1;
t579 = t241 / 0.2e1;
t578 = -t270 / 0.2e1;
t577 = t270 / 0.2e1;
t576 = -t275 / 0.2e1;
t575 = t275 / 0.2e1;
t573 = t276 / 0.2e1;
t566 = pkin(2) * t370;
t564 = pkin(4) * t276;
t561 = pkin(9) * t365;
t560 = pkin(9) * t369;
t342 = t358 * pkin(9);
t344 = t359 * pkin(3);
t550 = mrSges(7,1) * t241;
t546 = Ifges(3,4) * t367;
t545 = Ifges(3,4) * t371;
t544 = Ifges(4,4) * t276;
t542 = Ifges(5,4) * t369;
t541 = Ifges(6,6) * t365;
t532 = t275 * mrSges(4,3);
t526 = qJ(5) * t240;
t523 = t228 * t369;
t518 = t289 * t365;
t517 = t289 * t369;
t515 = t349 * t365;
t514 = t349 * t369;
t511 = t358 * t372;
t505 = t365 * t368;
t500 = t372 * t365;
t499 = t372 * t373;
t217 = pkin(3) * t288 - pkin(9) * t289 - t351;
t243 = t309 * t366 - t311 * t370;
t142 = t365 * t217 + t369 * t243;
t496 = t344 + t342;
t493 = qJD(2) * t367;
t478 = pkin(2) * t493;
t469 = qJD(2) * t373;
t468 = t289 * t488;
t451 = -t488 / 0.2e1;
t450 = t488 / 0.2e1;
t449 = -t487 / 0.2e1;
t448 = t487 / 0.2e1;
t40 = -t103 * mrSges(7,1) - t171 * mrSges(7,3);
t446 = -t351 - t344;
t445 = t486 / 0.2e1;
t42 = -t104 * mrSges(7,1) + t171 * mrSges(7,2);
t233 = t365 * t243;
t141 = t217 * t369 - t233;
t441 = t372 * t351 - t368 * t373;
t436 = pkin(4) * t508 + qJ(5) * t510 + t496;
t112 = -qJ(5) * t288 - t142;
t428 = pkin(4) * t518 - t616;
t424 = mrSges(3,1) * t367 + mrSges(3,2) * t371;
t422 = mrSges(4,1) * t358 + mrSges(4,2) * t359;
t417 = t546 + t631;
t416 = -Ifges(5,2) * t365 + t542;
t413 = Ifges(3,5) * t371 - Ifges(3,6) * t367;
t411 = -Ifges(6,2) * t369 + t541;
t407 = -qJ(5) * t369 + qJ(6) * t365;
t229 = qJD(2) * t289 + t393;
t140 = pkin(3) * t229 - pkin(9) * t228 + t478;
t295 = t367 * t469;
t296 = t371 * t469;
t149 = qJD(3) * t616 + t295 * t370 + t296 * t366;
t33 = t140 * t369 - t365 * t149 - t217 * t488 - t243 * t487;
t70 = t370 * t232 - t366 * t239 - t282 * t492 + t292 * t491;
t402 = pkin(3) * t507 + pkin(9) * t511 + t441;
t400 = t358 * pkin(5) + qJ(6) * t508 + t436;
t399 = pkin(1) * t424;
t396 = t468 - t523;
t395 = t367 * (Ifges(3,1) * t371 - t546);
t32 = t365 * t140 + t369 * t149 + t217 * t487 - t243 * t488;
t281 = -t369 * t554 + t447;
t263 = t359 * t505 + t502;
t265 = t359 * t500 - t504;
t389 = -g(1) * t265 - g(2) * t263 - g(3) * t513;
t67 = -pkin(3) * t483 - t70;
t150 = qJD(3) * t243 + t295 * t366 - t370 * t296;
t18 = -qJ(5) * t229 - qJD(5) * t288 - t32;
t379 = pkin(4) * t397 + qJ(5) * t468 + t150;
t377 = t103 * qJ(5) - t241 * qJD(5) + t67;
t15 = t104 * pkin(4) + t377;
t192 = Ifges(4,2) * t275 + t544 + t667;
t268 = Ifges(4,4) * t275;
t193 = Ifges(4,1) * t276 + t268 + t668;
t39 = -t270 * t554 + t647;
t46 = -t64 - t646;
t6 = t240 * qJD(6) + t104 * t554 + t377;
t374 = ((-t416 / 0.2e1 + t657 / 0.2e1) * t240 + t577 * t643 + t642) * qJD(4) - t15 * t529 - t540 * t588 - t6 * t530 + (-t411 / 0.2e1 + t658 / 0.2e1) * t489 + t543 * t587 + (t268 + t193) * t576 + Ifges(4,5) * t380 + (-t67 * mrSges(5,1) + t15 * mrSges(6,2) - t6 * mrSges(7,3) + Ifges(5,2) * t587 + t585 * t682 + t671 * t586 - t677) * t369 + (t67 * mrSges(5,2) - Ifges(6,2) * t588 - t585 * t604 + t672 * t589 + t676) * t365 + (-t541 + t538) * t586 + t663 * (t450 - t522 / 0.2e1) + t664 * (t448 - t521 / 0.2e1) + Ifges(4,3) * t483 + (-t537 + t542) * t589 + t219 * t531 + t218 * t532 + Ifges(4,6) * t172 - t69 * mrSges(4,2) + t70 * mrSges(4,1) + t192 * t573 + t132 * t449 + t133 * t451 - (Ifges(4,1) * t275 - t544 + t644) * t276 / 0.2e1 + (t275 * t643 + t276 * t605) * t578 + t651 * t575 + (t61 * t653 + t64 * t652 + t649) * mrSges(6,1) + (t39 * t653 - t46 * t652) * mrSges(7,1) + (-t652 * t95 + t653 * t94 + t650) * mrSges(5,3) + (-t668 / 0.2e1 + t307 * mrSges(4,2) + t411 * t579 + t416 * t581 + t657 * t582 + t658 * t580 - t642) * t275 + (t667 / 0.2e1 + t95 * mrSges(5,2) + t94 * mrSges(5,1) + t39 * mrSges(7,3) - t46 * mrSges(7,2) + t307 * mrSges(4,1) + t64 * mrSges(6,3) - t61 * mrSges(6,2) - Ifges(4,2) * t576 + Ifges(6,4) * t579 + Ifges(5,6) * t581 + t636 * t582 + t635 * t580) * t276;
t361 = t369 * pkin(5);
t360 = t365 * pkin(5);
t353 = Ifges(3,4) * t494;
t350 = -pkin(3) - t566;
t310 = t361 + t560;
t308 = t360 + t561;
t285 = t361 + t514;
t284 = t360 + t515;
t283 = t303 - t566;
t274 = Ifges(3,1) * t495 + Ifges(3,5) * qJD(2) + t353;
t273 = Ifges(3,6) * qJD(2) + qJD(1) * t417;
t266 = t359 * t502 + t505;
t264 = t359 * t504 - t500;
t254 = t281 - t566;
t246 = -mrSges(4,2) * t484 + t532;
t212 = -mrSges(4,1) * t275 + mrSges(4,2) * t276;
t175 = -mrSges(7,3) * t270 + t550;
t158 = -mrSges(4,2) * t483 + t172 * mrSges(4,3);
t157 = mrSges(4,1) * t483 - mrSges(4,3) * t380;
t156 = -mrSges(6,2) * t240 - mrSges(6,3) * t241;
t154 = pkin(4) * t241 + t526;
t153 = -mrSges(7,2) * t241 + mrSges(7,3) * t240;
t148 = -qJ(5) * t517 + t428;
t122 = t289 * t407 + t428;
t119 = -pkin(4) * t288 - t141;
t111 = t241 * t554 + t526;
t89 = -t126 - t564;
t88 = -t262 - t127;
t79 = -t120 - t564;
t78 = -t262 - t121;
t77 = -pkin(5) * t518 - t112;
t60 = t233 + (pkin(5) * t289 - t217) * t369 - t554 * t288;
t37 = mrSges(7,2) * t103 + mrSges(7,3) * t104;
t36 = mrSges(5,1) * t104 - mrSges(5,2) * t103;
t35 = -mrSges(6,2) * t104 + mrSges(6,3) * t103;
t34 = (-qJ(5) * t228 - qJD(5) * t289) * t369 + t379;
t19 = -pkin(4) * t229 - t33;
t17 = t407 * t228 + (qJD(6) * t365 + (qJ(6) * qJD(4) - qJD(5)) * t369) * t289 + t379;
t16 = -pkin(5) * t397 - t18;
t10 = -pkin(5) * t396 - qJD(6) * t288 - t229 * t554 - t33;
t1 = [(t636 * t229 + t396 * t686 - t671 * t397) * t581 + (t635 * t229 - t672 * t396 + t397 * t687) * t579 + (-m(4) * t441 - m(5) * t402 + t641 * (t266 * pkin(4) + t265 * qJ(5) + t402) + t673 * t511 + t601 * t372 + t602 * t368 - t596 * t266 - t606 * t265) * g(2) + (t598 / 0.2e1 + t267 * mrSges(4,1) - t69 * mrSges(4,3) - Ifges(4,4) * t380 + Ifges(6,4) * t588 - Ifges(4,2) * t172 - Ifges(4,6) * t483 + Ifges(5,6) * t587 + t585 * t605 + t586 * t636 + t589 * t635 + t593) * t288 + m(4) * (t149 * t219 + t243 * t69 - t267 * t351 - t307 * t478) + t298 * t545 / 0.2e1 + (-t218 * t228 - t219 * t229) * mrSges(4,3) + (-m(4) * t218 + m(5) * t197 - t619) * t150 + (t297 * t562 + t298 * t563 + t609) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t609) + (t274 * t569 + t413 * qJD(2) / 0.2e1 - t610) * qJD(2) + (t371 * t545 + t395) * t445 - t351 * (-t172 * mrSges(4,1) + mrSges(4,2) * t380) - (-m(4) * t70 + m(5) * t67 - t157 + t36) * t616 + (mrSges(6,1) * t9 - mrSges(5,3) * t14 + t676) * t517 + (t7 * mrSges(6,1) - t13 * mrSges(5,3) + t677) * t518 + t212 * t478 + m(7) * (t10 * t39 + t122 * t6 + t16 * t46 + t17 * t50 + t2 * t60 + t4 * t77) + m(6) * (t112 * t7 + t119 * t9 + t148 * t15 + t18 * t64 + t19 * t61 + t34 * t86) - t273 * t493 / 0.2e1 + t484 * (Ifges(4,5) * t228 - Ifges(4,6) * t229) / 0.2e1 + m(5) * (t13 * t142 + t14 * t141 + t32 * t95 - t33 * t94) + (Ifges(3,4) * t298 + Ifges(3,2) * t297) * t569 + (t229 * t605 + t396 * t604 - t397 * t682) * t577 - t306 * t525 + t417 * t638 + Ifges(2,3) * qJDD(1) - t307 * (mrSges(4,1) * t229 + mrSges(4,2) * t228) - pkin(1) * (-mrSges(3,1) * t297 + mrSges(3,2) * t298) + t243 * t158 + t149 * t246 - t229 * t192 / 0.2e1 + t228 * t193 / 0.2e1 + t16 * t177 + t19 * t178 + t32 * t179 + t33 * t180 + t10 * t175 + t18 * t176 + t34 * t156 + t148 * t35 + t17 * t153 + t141 * t44 + t142 * t45 - t399 * t486 + t112 * t41 + t119 * t43 + t122 * t37 + t77 * t42 + t60 * t40 + (Ifges(4,1) * t228 - Ifges(4,4) * t229) * t573 + (Ifges(4,4) * t228 - Ifges(4,2) * t229) * t575 + t644 * t229 / 0.2e1 - t651 * t228 / 0.2e1 + (-mrSges(3,1) * t563 - mrSges(3,2) * t562 + 0.2e1 * Ifges(3,6) * t569) * qJDD(2) + (Ifges(3,1) * t298 + Ifges(3,4) * t638 + Ifges(3,5) * qJDD(2) - t445 * t631) * t367 + (t267 * mrSges(4,2) - t70 * mrSges(4,3) + Ifges(4,1) * t380 + Ifges(4,4) * t172 + Ifges(4,5) * t483 + t132 * t450 + t133 * t449 + t15 * t420 + t411 * t588 + t416 * t587 + t6 * t419 + t67 * t421 + t643 * t585 + t657 * t586 + t658 * t589) * t289 + (m(5) * t499 + t641 * (-t264 * pkin(4) - t263 * qJ(5) - t499) + (m(4) * t373 + t602) * t372 + t596 * t264 + t606 * t263 + (-m(7) * t446 - (m(7) * t590 - mrSges(7,1)) * t358 + m(4) * t351 + (-m(5) - m(6)) * (t446 - t342) - t601 + t661) * t368) * g(1) + t663 * (t289 * t448 + t524 / 0.2e1) + t664 * (t289 * t451 + t523 / 0.2e1) - t94 * (mrSges(5,1) * t229 + mrSges(5,3) * t396) + t61 * (-mrSges(6,1) * t396 + mrSges(6,2) * t229) + t39 * (-mrSges(7,1) * t396 - mrSges(7,3) * t229) + t64 * (mrSges(6,1) * t397 - mrSges(6,3) * t229) + t50 * (mrSges(7,2) * t396 + mrSges(7,3) * t397) + t197 * (mrSges(5,1) * t397 - mrSges(5,2) * t396) + t46 * (-mrSges(7,1) * t397 + mrSges(7,2) * t229) + t95 * (-mrSges(5,2) * t229 - mrSges(5,3) * t397) + t86 * (-mrSges(6,2) * t397 + mrSges(6,3) * t396) + (Ifges(6,4) * t229 + Ifges(6,2) * t396 + Ifges(6,6) * t397) * t580 + (-Ifges(5,4) * t396 - Ifges(5,2) * t397 + Ifges(5,6) * t229) * t582; (t218 * t222 - t219 * t223 + t307 * t479 + (t366 * t69 + t370 * t70 + (-t218 * t366 + t219 * t370) * qJD(3)) * pkin(2)) * m(4) + t594 * t349 + (t120 * t94 - t121 * t95 - t197 * t222 + t350 * t67 + (t197 * t366 + (t365 * t94 + t369 * t95) * t370) * qJD(3) * pkin(2)) * m(5) + t625 * t175 + t626 * t177 + (t2 * t284 + t254 * t6 + t285 * t4 + t39 * t625 + t46 * t626 + t50 * t628) * m(7) + t628 * t153 + t618 * t156 + (-t61 * t79 - t64 * t78 + t15 * t283 + (t365 * t61 - t369 * t64) * t477 + t618 * t86) * m(6) + (-m(5) * t330 + t372 * t600 + t599) * g(1) + (t368 * t600 + t595) * g(2) + (-m(4) * t362 - m(5) * (t362 + t496) - m(7) * (t362 + t400) - m(6) * (t362 + t436) + t306 + t597) * g(3) + (-t79 + t438) * t178 + (t477 - t223) * t246 + t157 * t566 + t158 * t568 - t619 * t679 - t212 * t479 + t273 * t495 / 0.2e1 + t374 + (-t437 - t78) * t176 - (-Ifges(3,2) * t495 + t274 + t353) * t494 / 0.2e1 + (m(4) * t567 + t422 + t424) * t608 + (t610 + (-t395 / 0.2e1 + t399) * qJD(1)) * qJD(1) + t350 * t36 + Ifges(3,5) * t298 + Ifges(3,6) * t297 - t287 * mrSges(3,1) + t283 * t35 + t284 * t40 + t285 * t42 - t286 * mrSges(3,2) + t254 * t37 + Ifges(3,3) * qJDD(2) - t413 * t486 / 0.2e1 + t632 * t515 + t633 * t514 + t656 * t180 + t665 * t179; t594 * pkin(9) + (mrSges(5,1) * t472 + t368 * t674 + t595) * g(2) + (-pkin(3) * t67 + t126 * t94 - t127 * t95 - t197 * t219) * m(5) + t624 * t177 + (t2 * t308 + t281 * t6 + t310 * t4 + t39 * t623 + t46 * t624 + t50 * t627) * m(7) + t627 * t153 + t623 * t175 + t619 * t219 + (-m(5) * (-pkin(3) * t511 + t330) + mrSges(5,1) * t471 + t613 * t372 + t599) * g(1) + (t15 * t303 - t61 * t89 + t617 * t86 - t64 * t88) * m(6) + t617 * t156 + t608 * t422 + (-m(5) * t496 - m(6) * t436 - m(7) * t400 + t597) * g(3) + t374 + t308 * t40 + t310 * t42 + t303 * t35 + t281 * t37 - t218 * t246 - t89 * t178 - t127 * t179 - t126 * t180 - t88 * t176 - pkin(3) * t36 + t632 * t561 + t633 * t560; t401 * t177 - t554 * t40 + (Ifges(6,2) * t240 + t133 + t533) * t579 + (-t41 + t42) * qJ(5) + (-pkin(4) * t9 - qJ(5) * t7 - qJD(5) * t64 - t154 * t86) * m(6) + t620 * qJD(5) + t621 * t175 + t46 * t550 + t39 * t551 + t61 * t553 + (-m(6) * t61 + t547 - t680) * t95 + (-m(6) * t64 + t548 - t681) * t94 + t593 - t64 * t552 + (t240 * t604 - t241 * t682) * t578 + t598 + (-t241 * t671 + t132 + t236 - t539) * t582 + (-t240 * t672 + t235 - t534 + t663) * t580 - t197 * (mrSges(5,1) * t241 - mrSges(5,2) * t240) - t50 * (mrSges(7,2) * t240 + mrSges(7,3) * t241) - t86 * (-mrSges(6,2) * t241 + mrSges(6,3) * t240) - t154 * t156 - t111 * t153 - pkin(4) * t43 + (qJ(5) * t4 - t111 * t50 - t2 * t554 + t621 * t39 + t46 * t647) * m(7) + (t641 * (-t265 * pkin(4) + qJ(5) * t266) - t606 * t266 + t596 * t265) * g(1) + (t641 * (-t263 * pkin(4) + qJ(5) * t264) - t606 * t264 + t596 * t263) * g(2) + (t641 * qJ(5) * t512 + (t421 - t527 - t528 + (m(6) * pkin(4) - mrSges(6,2) - t430) * t365) * t358) * g(3) + (-Ifges(5,2) * t241 - t237 + t664) * t581; -t620 * t270 + (t153 + t156) * t241 + t43 + t40 + (t241 * t50 - t270 * t46 + t2 + t389) * m(7) + (t241 * t86 + t270 * t64 + t389 + t9) * m(6); -t240 * t153 + t270 * t175 + (-g(1) * t266 - g(2) * t264 - g(3) * t512 - t240 * t50 + t270 * t39 + t4) * m(7) + t42;];
tau  = t1;
