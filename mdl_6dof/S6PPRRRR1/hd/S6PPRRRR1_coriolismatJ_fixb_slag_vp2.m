% Calculate matrix of centrifugal and coriolis load on the joints for
% S6PPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6PPRRRR1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:59
% EndTime: 2019-03-08 19:00:13
% DurationCPUTime: 10.99s
% Computational Cost: add. (21409->542), mult. (56832->764), div. (0->0), fcn. (67577->14), ass. (0->344)
t345 = sin(qJ(6));
t347 = cos(qJ(6));
t346 = sin(qJ(4));
t348 = cos(qJ(4));
t344 = sin(pkin(7));
t586 = sin(qJ(3));
t474 = t586 * t344;
t539 = cos(pkin(7));
t306 = t346 * t539 + t348 * t474;
t385 = t346 * t474 - t348 * t539;
t585 = sin(qJ(5));
t587 = cos(qJ(5));
t371 = t306 * t587 - t385 * t585;
t588 = cos(qJ(3));
t478 = t344 * t588;
t201 = -t345 * t478 + t347 * t371;
t241 = t306 * t585 + t385 * t587;
t547 = t347 * mrSges(7,2);
t551 = t345 * mrSges(7,1);
t409 = t551 / 0.2e1 + t547 / 0.2e1;
t391 = t241 * t409;
t441 = t547 + t551;
t402 = t241 * t441;
t550 = t345 * mrSges(7,3);
t482 = t550 / 0.2e1;
t483 = -t550 / 0.2e1;
t643 = t482 + t483;
t667 = t391 + t402 / 0.2e1 + t643 * t201;
t538 = sin(pkin(6));
t442 = t538 * sin(pkin(13));
t443 = cos(pkin(13)) * t538;
t540 = cos(pkin(6));
t644 = t344 * t540 + t443 * t539;
t364 = t442 * t588 + t586 * t644;
t376 = -t344 * t443 + t539 * t540;
t190 = t346 * t376 + t348 * t364;
t357 = t346 * t364 - t348 * t376;
t123 = t190 * t585 + t357 * t587;
t527 = t201 * t347;
t200 = -t345 * t371 - t347 * t478;
t528 = t200 * t345;
t407 = t371 - t527 + t528;
t352 = t190 * t587 - t357 * t585;
t268 = t442 * t586 - t588 * t644;
t85 = t268 * t345 + t347 * t352;
t541 = t85 * t347;
t84 = t268 * t347 - t345 * t352;
t543 = t84 * t345;
t415 = t352 - t541 + t543;
t607 = m(7) / 0.2e1;
t688 = (t123 * t407 + t241 * t415) * t607;
t691 = t407 * m(7) * t241;
t692 = qJD(1) * t688 + qJD(2) * t691;
t698 = qJD(6) * t667 + t692;
t392 = t123 * t409;
t403 = t123 * t441;
t666 = t392 + t403 / 0.2e1 + t643 * t85;
t690 = t415 * m(7) * t123;
t693 = qJD(1) * t690 + qJD(2) * t688;
t697 = qJD(6) * t666 + t693;
t689 = qJD(4) + qJD(5);
t450 = t588 * t585;
t451 = t588 * t587;
t292 = (t346 * t451 + t348 * t450) * t344;
t293 = (-t346 * t450 + t348 * t451) * t344;
t323 = -mrSges(7,1) * t347 + mrSges(7,2) * t345;
t661 = t323 / 0.2e1 - mrSges(6,1) / 0.2e1;
t696 = t293 * mrSges(6,2) / 0.2e1 - t661 * t292;
t473 = t585 * t346;
t315 = -t348 * t587 + t473;
t475 = t587 * t346;
t316 = -t348 * t585 - t475;
t335 = -pkin(4) * t348 - pkin(3);
t262 = pkin(5) * t315 + pkin(11) * t316 + t335;
t657 = pkin(10) + pkin(9);
t503 = t657 * t348;
t618 = -t473 * t657 + t503 * t587;
t165 = t262 * t347 - t345 * t618;
t166 = t262 * t345 + t347 * t618;
t423 = -t165 * t345 + t166 * t347;
t408 = t618 - t423;
t695 = t408 * t607;
t257 = t293 * t347 + t345 * t474;
t521 = t257 * t347;
t256 = -t293 * t345 + t347 * t474;
t522 = t256 * t345;
t694 = (t522 / 0.2e1 - t521 / 0.2e1) * mrSges(7,3) + t696;
t152 = t316 * t268;
t153 = t315 * t268;
t687 = t661 * t152 - t153 * mrSges(6,2) / 0.2e1;
t589 = t347 / 0.2e1;
t623 = t346 ^ 2 + t348 ^ 2;
t684 = mrSges(5,3) * t623 - mrSges(4,2);
t573 = Ifges(7,2) * t345;
t576 = Ifges(7,4) * t347;
t438 = -t573 + t576;
t205 = -Ifges(7,6) * t316 - t315 * t438;
t577 = Ifges(7,4) * t345;
t440 = Ifges(7,1) * t347 - t577;
t206 = -Ifges(7,5) * t316 - t315 * t440;
t401 = t316 * (Ifges(7,5) * t345 + Ifges(7,6) * t347);
t591 = t345 / 0.2e1;
t439 = Ifges(7,1) * t345 + t576;
t590 = -t347 / 0.2e1;
t675 = Ifges(7,2) * t589 + t577 / 0.2e1;
t640 = t345 * t675 + t439 * t590;
t397 = t205 * t589 + t206 * t591 - t401 / 0.2e1 + Ifges(6,6) * t316 + (-Ifges(6,5) + t640) * t315;
t645 = t618 * t323;
t648 = t618 * mrSges(6,1);
t285 = t475 * t657 + t503 * t585;
t673 = t285 * mrSges(6,2);
t681 = t397 + t645 - t648 + t673;
t678 = -t645 / 0.2e1 + t648 / 0.2e1 - t673 / 0.2e1;
t674 = m(6) * t335;
t671 = t285 * t345;
t670 = t285 * t347;
t629 = t285 * t352;
t628 = t285 * t371;
t669 = t285 * t585;
t518 = t618 * t285;
t261 = t441 * t316;
t554 = t316 * mrSges(6,3);
t668 = -t261 - t554;
t340 = t345 ^ 2;
t342 = t347 ^ 2;
t502 = t340 + t342;
t626 = t352 * t323;
t631 = t352 * mrSges(6,1);
t650 = t123 * mrSges(6,2);
t665 = t626 - t631 + t650;
t625 = t371 * t323;
t630 = t371 * mrSges(6,1);
t649 = t241 * mrSges(6,2);
t664 = t625 - t630 + t649;
t663 = -t626 / 0.2e1 + t631 / 0.2e1 - t650 / 0.2e1;
t662 = -t625 / 0.2e1 + t630 / 0.2e1 - t649 / 0.2e1;
t260 = t441 * t315;
t264 = mrSges(7,2) * t316 + t315 * t550;
t546 = t347 * mrSges(7,3);
t266 = -mrSges(7,1) * t316 + t315 * t546;
t269 = -mrSges(6,1) * t316 - mrSges(6,2) * t315;
t571 = Ifges(7,6) * t345;
t574 = Ifges(7,5) * t347;
t414 = -t574 / 0.2e1 + t571 / 0.2e1;
t436 = -t571 + t574;
t627 = t315 * t436;
t658 = (Ifges(6,4) * t315 - t627 + (-Ifges(7,3) + Ifges(6,1) - Ifges(6,2) + t342 * Ifges(7,1) / 0.2e1 + (-t576 + t573 / 0.2e1) * t345) * t316) * t315 - t618 * t261 + (t205 * t591 + t206 * t590 + (-Ifges(6,4) - t414) * t316) * t316 + t165 * t266 + t166 * t264 - t285 * t260 + t335 * t269;
t598 = -t260 / 0.2e1;
t656 = t264 / 0.2e1;
t655 = t266 / 0.2e1;
t654 = t268 / 0.2e1;
t545 = t348 * mrSges(5,2);
t325 = mrSges(5,1) * t346 + t545;
t653 = -t325 / 0.2e1;
t652 = pkin(5) * t618;
t651 = pkin(9) * t623;
t496 = t587 * pkin(4);
t334 = -t496 - pkin(5);
t646 = t334 * t618;
t270 = mrSges(6,1) * t315 - mrSges(6,2) * t316;
t324 = -mrSges(5,1) * t348 + mrSges(5,2) * t346;
t636 = -mrSges(4,1) + t270 + t324 + t674;
t634 = -m(6) / 0.2e1;
t633 = t441 / 0.2e1;
t608 = -m(7) / 0.2e1;
t632 = pkin(5) * t608;
t495 = t585 * pkin(4);
t555 = t315 * mrSges(6,3);
t333 = t495 + pkin(11);
t624 = t502 * t333;
t512 = t316 * t345;
t265 = -mrSges(7,2) * t315 + mrSges(7,3) * t512;
t504 = t347 * t265;
t511 = t316 * t347;
t267 = mrSges(7,1) * t315 + mrSges(7,3) * t511;
t507 = t345 * t267;
t404 = -t504 / 0.2e1 + t507 / 0.2e1;
t552 = t342 * mrSges(7,3);
t553 = t340 * mrSges(7,3);
t622 = -t552 / 0.2e1 - t553 / 0.2e1;
t578 = mrSges(7,3) * t316;
t621 = t267 * t590 + t502 * t578 / 0.2e1;
t620 = (-Ifges(7,1) + Ifges(7,2)) * t347;
t597 = -t261 / 0.2e1;
t619 = t200 * t655 + t201 * t656 + t371 * t597;
t617 = pkin(11) * t502 * t608 + t622;
t606 = m(6) * pkin(4);
t616 = m(7) * t334 - t587 * t606;
t129 = -t153 * t345 + t347 * t364;
t130 = t153 * t347 + t345 * t364;
t481 = t546 / 0.2e1;
t614 = t129 * t483 + t130 * t481 + t687;
t613 = t256 * t483 + t257 * t481 - t696;
t612 = t123 * t598 + t269 * t654 + t655 * t84 + t656 * t85;
t611 = m(7) * t624 + t585 * t606 + t552 + t553;
t610 = (-t624 + t495) * t607 + t622;
t609 = m(6) / 0.2e1;
t605 = m(7) * pkin(4);
t604 = mrSges(7,1) / 0.2e1;
t603 = -mrSges(7,2) / 0.2e1;
t602 = -Ifges(7,3) / 0.2e1;
t600 = t241 / 0.2e1;
t259 = t323 * t316;
t599 = t259 / 0.2e1;
t596 = t265 / 0.2e1;
t595 = -t267 / 0.2e1;
t594 = -t285 / 0.2e1;
t593 = t333 / 0.2e1;
t592 = -t345 / 0.2e1;
t584 = pkin(4) * t346;
t583 = pkin(5) * t259;
t582 = pkin(5) * t260;
t581 = pkin(11) * t265;
t575 = Ifges(7,5) * t315;
t572 = Ifges(7,6) * t315;
t535 = t123 * t152;
t532 = t129 * t345;
t531 = t130 * t347;
t271 = -pkin(5) * t316 + pkin(11) * t315;
t263 = t271 + t584;
t171 = t263 * t347 + t671;
t530 = t171 * t345;
t172 = t263 * t345 - t670;
t529 = t172 * t347;
t524 = t241 * t292;
t519 = t268 * t346;
t516 = t285 * t152;
t515 = t285 * t292;
t510 = t334 * t259;
t509 = t334 * t260;
t508 = t345 * t266;
t505 = t347 * t264;
t501 = t606 / 0.2e1;
t500 = pkin(11) * t508;
t499 = pkin(11) * t505;
t498 = -t588 / 0.2e1;
t494 = mrSges(7,3) * t529;
t493 = Ifges(7,2) / 0.2e1 - Ifges(7,1) / 0.2e1;
t492 = t333 * t505;
t487 = -t554 / 0.2e1;
t486 = t554 / 0.2e1;
t477 = t345 * t587;
t476 = t347 * t587;
t469 = t333 * t592;
t466 = t504 / 0.2e1;
t462 = t502 * t123;
t461 = t502 * t241;
t456 = t352 * t487;
t453 = t438 * t589 + t440 * t591 - t640;
t452 = t348 * t478;
t448 = -t477 / 0.2e1;
t447 = t269 * t498;
t446 = t325 * t498;
t445 = t352 * t486 + t612;
t444 = mrSges(7,3) * (-t340 / 0.2e1 - t342 / 0.2e1);
t356 = t357 * t346;
t361 = t364 * t478;
t384 = t385 * t346;
t425 = t123 * t292 + t152 * t241;
t13 = (t129 * t200 + t130 * t201 + t256 * t84 + t257 * t85 + t425) * t607 + (t153 * t371 + t268 * t474 + t293 * t352 - t361 + t425) * t609 + m(5) * (t190 * t452 + t356 * t478 - t361 + (-t306 * t348 - t384 + t474) * t268) / 0.2e1;
t179 = t268 * t364;
t20 = m(7) * (t129 * t84 + t130 * t85 + t535) + m(6) * (t153 * t352 + t179 + t535) + m(5) * (t179 + (-t190 * t348 - t356) * t268);
t431 = qJD(1) * t20 + qJD(2) * t13;
t417 = t344 ^ 2 * t588 * t586;
t55 = m(7) * (t200 * t256 + t201 * t257 + t524) + m(6) * (t293 * t371 - t417 + t524) + m(5) * (t306 * t452 + t384 * t478 - t417);
t430 = t13 * qJD(1) + t55 * qJD(2);
t426 = t123 * t618 + t629;
t424 = t531 - t532;
t422 = t529 - t530;
t181 = t271 * t347 + t671;
t182 = t271 * t345 - t670;
t421 = -t181 * t345 + t182 * t347;
t418 = t521 - t522;
t416 = t344 * t447;
t413 = t129 * t604 + t130 * t603;
t412 = t171 * t604 + t172 * t603;
t411 = t181 * t604 + t182 * t603;
t410 = -t256 * mrSges(7,1) / 0.2e1 + t257 * mrSges(7,2) / 0.2e1;
t406 = t577 + t620;
t400 = t334 * t441;
t395 = t502 * t587;
t14 = -pkin(3) * t325 + t171 * t267 + t172 * t265 + (-Ifges(5,4) * t346 + pkin(4) * t270) * t346 + m(7) * (t165 * t171 + t166 * t172 + t518) + t584 * t674 + (Ifges(5,4) * t348 + (Ifges(5,1) - Ifges(5,2)) * t346) * t348 + t658;
t350 = (pkin(4) * t519 + t426 - t629) * t634 + (t171 * t84 + t172 * t85 + t426) * t608 + t268 * t653 + (t261 / 0.2e1 + t486) * t352 + (-t423 * t608 - t618 * t634 - t404) * t123;
t359 = (t152 * t334 + t333 * t424) * t607 + (-t152 * t587 + t153 * t585) * t501 + mrSges(5,1) * t519 / 0.2e1 + t545 * t654 + t614;
t3 = t359 + t350 + t456 - t612;
t351 = -t478 * t584 * t609 + (t171 * t200 + t172 * t201 + t628) * t607 + t619 + (t598 + t695) * t241;
t360 = (t292 * t334 + t333 * t418) * t607 + (-t292 * t587 + t293 * t585) * t501 + t478 * t653;
t9 = t351 + t344 * t446 + t416 + t507 * t600 - t360 + (t486 + t487) * t371 - t466 * t241 + t694;
t394 = -t3 * qJD(1) + t9 * qJD(2) + t14 * qJD(3);
t393 = t181 * t84 + t182 * t85 + t629;
t374 = t404 + t695;
t353 = (t598 + t374) * t241 + (t181 * t200 + t182 * t201 + t628) * t607 + t416 + t619;
t379 = m(7) * (-pkin(5) * t292 + pkin(11) * t418);
t16 = -t379 / 0.2e1 + t353 + t694;
t19 = m(7) * (t165 * t181 + t166 * t182 + t518) + t182 * t265 + t181 * t267 + t658;
t380 = m(7) * (-pkin(5) * t152 + pkin(11) * t424);
t5 = (t597 + t487) * t352 + t404 * t123 + (-t531 / 0.2e1 + t532 / 0.2e1) * mrSges(7,3) + (t123 * t408 + t393) * t607 - t380 / 0.2e1 + t445 - t687;
t390 = qJD(1) * t5 + qJD(2) * t16 + qJD(3) * t19;
t369 = (t541 / 0.2e1 - t543 / 0.2e1) * t578 + t123 * t599 + t84 * t596 + t85 * t595;
t10 = t369 - t413;
t35 = t165 * t265 - t166 * t267 + t285 * t259 + ((mrSges(7,3) * t166 - Ifges(7,4) * t511 + t572) * t347 + (-mrSges(7,3) * t165 + t316 * t406 + t575) * t345) * t316;
t368 = (-t528 / 0.2e1 + t527 / 0.2e1) * t578 + t200 * t596 + t201 * t595 + t241 * t599;
t41 = t368 + t410;
t389 = qJD(1) * t10 + qJD(2) * t41 + qJD(3) * t35;
t375 = pkin(5) * t633 - t400 / 0.2e1;
t377 = (mrSges(7,1) * t448 + t476 * t603) * pkin(4);
t132 = t377 + (-t342 + t340) * Ifges(7,4) + t375 + t345 * t620;
t173 = (mrSges(7,2) * pkin(5) - t576) * t347 + (mrSges(7,1) * pkin(5) + t406) * t345;
t373 = t572 + (t345 * t493 - 0.2e1 * t576) * t316 + mrSges(7,1) * t594;
t378 = mrSges(7,2) * t594 - t493 * t511 - t575;
t38 = t583 / 0.2e1 + (pkin(11) * t444 + t602) * t316 + (pkin(11) * t267 / 0.2e1 + t378) * t347 + (t581 / 0.2e1 + t373) * t345 + t411;
t388 = qJD(3) * t38 + qJD(4) * t132 + qJD(5) * t173;
t365 = (t323 - mrSges(6,1)) * t495 + (mrSges(7,3) * t502 - mrSges(6,2)) * t496;
t141 = (t333 * t395 + t334 * t585) * t605 + t365;
t349 = (t646 + t421 * t333 + (-t165 * t477 + t166 * t476 + t669) * pkin(4)) * t607 - t509 / 0.2e1 + t181 * t483 + t182 * t481 + t266 * t469 + t492 / 0.2e1 + t495 * t597 + pkin(4) * t267 * t448 + t466 * t496 - t678;
t358 = (pkin(11) * t422 - t652) * t608 - t582 / 0.2e1 + t171 * t482 - t494 / 0.2e1 + t500 / 0.2e1 - t499 / 0.2e1 + t678;
t22 = t349 + t358;
t354 = (t334 * t371 + (-t200 * t477 + t201 * t476) * pkin(4)) * t607 + t610 * t241 - t662;
t362 = -t241 * t617 - t371 * t632 + t662;
t37 = t354 + t362;
t355 = (t334 * t352 + (t476 * t85 - t477 * t84) * pkin(4)) * t607 + t610 * t123 - t663;
t363 = -t123 * t617 - t352 * t632 + t663;
t7 = t355 + t363;
t383 = qJD(1) * t7 + qJD(2) * t37 + qJD(3) * t22 + qJD(4) * t141;
t158 = t400 + t453;
t28 = -t403 / 0.2e1 + t392;
t32 = -t510 / 0.2e1 + (t333 * t444 + t602) * t316 + (t267 * t593 + t378) * t347 + (t265 * t593 + t373) * t345 + t412;
t57 = -t402 / 0.2e1 + t391;
t381 = -qJD(1) * t28 - qJD(2) * t57 - qJD(3) * t32 + qJD(4) * t158;
t367 = t414 * t315 + t316 * t602 - t345 * (-t316 * t438 + t572) / 0.4e1 + t347 * (-t316 * t440 + t575) / 0.4e1 + t285 * t633 + t627 / 0.4e1 + (t675 - t440 / 0.4e1) * t511 + (0.2e1 * t439 + t438) * t512 / 0.4e1 + t643 * t166;
t133 = t377 - t375 + t453;
t40 = t368 - t410;
t39 = t367 - t583 / 0.2e1 + t581 * t592 + t411 + t621 * pkin(11);
t36 = t354 - t362;
t33 = t367 + t510 / 0.2e1 + t265 * t469 + t412 + t621 * t333;
t18 = -t358 + t349 + t397;
t15 = t379 / 0.2e1 + t353 + t613;
t11 = t369 + t413;
t8 = t351 + (-t241 / 0.2e1 + t600) * t555 + t404 * t241 + (t447 + t446) * t344 + t360 + t613;
t6 = t355 - t363;
t4 = t393 * t607 + t352 * t597 + t456 + t380 / 0.2e1 + t374 * t123 + t445 + t614;
t2 = t359 - t350 + t445;
t1 = t13 * qJD(3) + t688 * t689;
t12 = [t20 * qJD(3) + t689 * t690, t1, t2 * qJD(4) + t4 * qJD(5) + t11 * qJD(6) + t431 + (t130 * t265 + t129 * t267 + m(7) * (t129 * t165 + t130 * t166 + t516) - t153 * t555 + m(6) * (t153 * t618 + t516) + (-m(5) * pkin(3) + t636) * t364 + t668 * t152 + (-m(5) * t651 - t684) * t268) * qJD(3), t2 * qJD(3) + (-t190 * mrSges(5,1) + t357 * mrSges(5,2) - t123 * t611 + t352 * t616 + t665) * qJD(4) + t6 * qJD(5) + t697, t4 * qJD(3) + t6 * qJD(4) + (m(7) * (-pkin(5) * t352 - pkin(11) * t462) - mrSges(7,3) * t462 + t665) * qJD(5) + t697, t11 * qJD(3) + (-mrSges(7,1) * t85 - mrSges(7,2) * t84) * qJD(6) + t689 * t666; t1, t55 * qJD(3) + t689 * t691, t8 * qJD(4) + t15 * qJD(5) + t40 * qJD(6) + t430 + (m(7) * (t165 * t256 + t166 * t257 + t515) + t256 * t267 + t257 * t265 - t293 * t555 + m(6) * (t293 * t618 + t515) + m(5) * (-t586 * pkin(3) + t588 * t651) * t344 + t636 * t474 + t668 * t292 + t684 * t478) * qJD(3), t8 * qJD(3) + (-t306 * mrSges(5,1) + t385 * mrSges(5,2) - t241 * t611 + t371 * t616 + t664) * qJD(4) + t36 * qJD(5) + t698, t15 * qJD(3) + t36 * qJD(4) + (m(7) * (-pkin(5) * t371 - pkin(11) * t461) - mrSges(7,3) * t461 + t664) * qJD(5) + t698, t40 * qJD(3) + (-mrSges(7,1) * t201 - mrSges(7,2) * t200) * qJD(6) + t689 * t667; -qJD(4) * t3 + qJD(5) * t5 + qJD(6) * t10 - t431, qJD(4) * t9 + qJD(5) * t16 + qJD(6) * t41 - t430, qJD(4) * t14 + qJD(5) * t19 + qJD(6) * t35 ((-t587 * t618 - t669) * t606 + Ifges(5,5) * t348 - Ifges(5,6) * t346 - t509 + t496 * t555 + t495 * t554 + m(7) * (t333 * t422 + t646) - t333 * t508 - mrSges(7,3) * t530 + t492 + t494 + t324 * pkin(9) + t681) * qJD(4) + t18 * qJD(5) + t33 * qJD(6) + t394, t18 * qJD(4) + (m(7) * (pkin(11) * t421 - t652) + t499 - t500 + t582 + t421 * mrSges(7,3) + t681) * qJD(5) + t39 * qJD(6) + t390, t33 * qJD(4) + t39 * qJD(5) + (-mrSges(7,1) * t166 - mrSges(7,2) * t165 + t401) * qJD(6) + t389; qJD(3) * t3 + qJD(5) * t7 - qJD(6) * t28 - t693, -qJD(3) * t9 + qJD(5) * t37 - qJD(6) * t57 - t692, qJD(5) * t22 - qJD(6) * t32 - t394, qJD(5) * t141 + qJD(6) * t158 ((-pkin(5) * t585 + pkin(11) * t395) * t605 + t365) * qJD(5) + t133 * qJD(6) + t383, t133 * qJD(5) + (t323 * t333 + t436) * qJD(6) + t381; -qJD(3) * t5 - qJD(4) * t7 - t693, -qJD(3) * t16 - qJD(4) * t37 - t692, -qJD(4) * t22 - qJD(6) * t38 - t390, -qJD(6) * t132 - t383, -t173 * qJD(6) (pkin(11) * t323 + t436) * qJD(6) - t388; -qJD(3) * t10 + qJD(4) * t28, -qJD(3) * t41 + qJD(4) * t57, qJD(4) * t32 + qJD(5) * t38 - t389, qJD(5) * t132 - t381, t388, 0;];
Cq  = t12;
