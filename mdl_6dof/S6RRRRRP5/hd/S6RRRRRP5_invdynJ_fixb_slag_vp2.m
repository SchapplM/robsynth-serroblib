% Calculate vector of inverse dynamics joint torques for
% S6RRRRRP5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRRRP5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:11
% EndTime: 2019-03-10 01:19:39
% DurationCPUTime: 51.93s
% Computational Cost: add. (23632->947), mult. (52528->1226), div. (0->0), fcn. (37995->14), ass. (0->426)
t387 = sin(qJ(2));
t392 = cos(qJ(2));
t435 = pkin(2) * t387 - pkin(8) * t392;
t329 = t435 * qJD(1);
t391 = cos(qJ(3));
t386 = sin(qJ(3));
t478 = qJD(1) * t387;
t454 = t386 * t478;
t244 = pkin(7) * t454 + t391 * t329;
t488 = t391 * t392;
t415 = pkin(3) * t387 - pkin(9) * t488;
t394 = -pkin(9) - pkin(8);
t455 = qJD(3) * t394;
t684 = -qJD(1) * t415 + t391 * t455 - t244;
t300 = t386 * t329;
t491 = t387 * t391;
t493 = t386 * t392;
t683 = t300 + (-pkin(7) * t491 - pkin(9) * t493) * qJD(1) - t386 * t455;
t385 = sin(qJ(4));
t390 = cos(qJ(4));
t324 = t385 * t391 + t386 * t390;
t609 = qJD(3) + qJD(4);
t239 = t609 * t324;
t411 = t324 * t392;
t268 = qJD(1) * t411;
t682 = t239 - t268;
t645 = -mrSges(7,1) - mrSges(6,1);
t644 = -mrSges(7,2) - mrSges(6,2);
t345 = t394 * t386;
t346 = t394 * t391;
t469 = qJD(4) * t390;
t470 = qJD(4) * t385;
t620 = t345 * t469 + t346 * t470 + t385 * t684 - t683 * t390;
t243 = t385 * t345 - t390 * t346;
t619 = -qJD(4) * t243 + t683 * t385 + t390 * t684;
t475 = qJD(2) * t391;
t321 = -t454 + t475;
t452 = t391 * t478;
t322 = qJD(2) * t386 + t452;
t228 = t321 * t385 + t322 * t390;
t384 = sin(qJ(5));
t389 = cos(qJ(5));
t437 = t390 * t321 - t322 * t385;
t162 = t228 * t389 + t384 * t437;
t681 = pkin(5) * t162;
t642 = Ifges(6,4) + Ifges(7,4);
t654 = qJ(6) * t162;
t416 = t385 * t386 - t390 * t391;
t238 = t609 * t416;
t410 = t416 * t392;
t269 = qJD(1) * t410;
t680 = -pkin(4) * t478 + t619 + (t238 - t269) * pkin(10);
t679 = t682 * pkin(10) - t620;
t477 = qJD(1) * t392;
t356 = qJD(3) - t477;
t348 = qJD(4) + t356;
t339 = qJD(5) + t348;
t548 = t339 / 0.2e1;
t567 = t162 / 0.2e1;
t648 = -t228 * t384 + t389 * t437;
t571 = t648 / 0.2e1;
t369 = pkin(7) * t478;
t343 = -qJD(2) * pkin(2) + t369;
t248 = -pkin(3) * t321 + t343;
t176 = -pkin(4) * t437 + t248;
t436 = pkin(2) * t392 + pkin(8) * t387;
t338 = -pkin(1) - t436;
t308 = t338 * qJD(1);
t370 = pkin(7) * t477;
t344 = qJD(2) * pkin(8) + t370;
t233 = t391 * t308 - t344 * t386;
t191 = -pkin(9) * t322 + t233;
t178 = pkin(3) * t356 + t191;
t234 = t308 * t386 + t344 * t391;
t192 = pkin(9) * t321 + t234;
t184 = t385 * t192;
t116 = t390 * t178 - t184;
t655 = pkin(10) * t228;
t92 = t116 - t655;
t84 = pkin(4) * t348 + t92;
t186 = t390 * t192;
t117 = t178 * t385 + t186;
t646 = pkin(10) * t437;
t93 = t117 + t646;
t87 = t384 * t93;
t46 = t389 * t84 - t87;
t29 = t46 - t654;
t25 = pkin(5) * t339 + t29;
t99 = -pkin(5) * t648 + qJD(6) + t176;
t596 = -t176 * mrSges(6,2) - t99 * mrSges(7,2) + mrSges(6,3) * t46 + mrSges(7,3) * t25;
t641 = Ifges(6,5) + Ifges(7,5);
t643 = Ifges(6,1) + Ifges(7,1);
t653 = t642 * t648;
t629 = t162 * t643 + t339 * t641 + t653;
t678 = -t596 + t548 * t641 + t567 * t643 + t571 * t642 + t629 / 0.2e1;
t89 = t389 * t93;
t47 = t384 * t84 + t89;
t622 = qJ(6) * t648;
t30 = t47 + t622;
t597 = -t176 * mrSges(6,1) - t99 * mrSges(7,1) + mrSges(6,3) * t47 + mrSges(7,3) * t30;
t639 = Ifges(6,6) + Ifges(7,6);
t640 = Ifges(6,2) + Ifges(7,2);
t652 = t642 * t162;
t630 = t339 * t639 + t640 * t648 + t652;
t677 = t597 + t548 * t639 + t567 * t642 + t571 * t640 + t630 / 0.2e1;
t676 = mrSges(4,3) + mrSges(5,3) + mrSges(6,3) + mrSges(7,3);
t466 = qJD(1) * qJD(2);
t333 = qJDD(1) * t392 - t387 * t466;
t318 = qJDD(3) - t333;
t307 = qJDD(4) + t318;
t291 = qJDD(5) + t307;
t334 = qJDD(1) * t387 + t392 * t466;
t218 = qJD(3) * t321 + qJDD(2) * t386 + t334 * t391;
t219 = -qJD(3) * t322 + qJDD(2) * t391 - t334 * t386;
t111 = qJD(4) * t437 + t218 * t390 + t219 * t385;
t112 = -qJD(4) * t228 - t218 * t385 + t219 * t390;
t44 = qJD(5) * t648 + t111 * t389 + t112 * t384;
t45 = -qJD(5) * t162 - t111 * t384 + t112 * t389;
t638 = Ifges(6,3) + Ifges(7,3);
t615 = t291 * t638 + t44 * t641 + t45 * t639;
t504 = qJDD(1) * pkin(1);
t237 = -pkin(2) * t333 - pkin(8) * t334 - t504;
t316 = t333 * pkin(7);
t286 = qJDD(2) * pkin(8) + t316;
t471 = qJD(3) * t391;
t473 = qJD(3) * t386;
t137 = t386 * t237 + t391 * t286 + t308 * t471 - t344 * t473;
t105 = pkin(9) * t219 + t137;
t138 = -qJD(3) * t234 + t391 * t237 - t286 * t386;
t95 = pkin(3) * t318 - pkin(9) * t218 + t138;
t28 = -qJD(4) * t117 - t105 * t385 + t390 * t95;
t19 = pkin(4) * t307 - pkin(10) * t111 + t28;
t27 = t390 * t105 + t178 * t469 - t192 * t470 + t385 * t95;
t21 = pkin(10) * t112 + t27;
t6 = -qJD(5) * t47 + t389 * t19 - t21 * t384;
t2 = pkin(5) * t291 - qJ(6) * t44 - qJD(6) * t162 + t6;
t467 = qJD(5) * t389;
t468 = qJD(5) * t384;
t5 = t384 * t19 + t389 * t21 + t84 * t467 - t468 * t93;
t3 = qJ(6) * t45 + qJD(6) * t648 + t5;
t665 = -t6 * mrSges(6,1) - t2 * mrSges(7,1) + t5 * mrSges(6,2) + t3 * mrSges(7,2);
t400 = t615 - t665;
t457 = Ifges(5,5) * t111 + Ifges(5,6) * t112 + Ifges(5,3) * t307;
t549 = -t339 / 0.2e1;
t568 = -t162 / 0.2e1;
t572 = -t648 / 0.2e1;
t659 = -t630 / 0.2e1;
t668 = -t28 * mrSges(5,1) + t27 * mrSges(5,2);
t675 = -(t549 * t639 + t568 * t642 + t572 * t640 - t597 + t659) * t162 + t400 + t457 - t668;
t554 = t291 / 0.2e1;
t588 = t45 / 0.2e1;
t589 = t44 / 0.2e1;
t656 = t554 * t639 + t588 * t640 + t589 * t642;
t657 = t554 * t641 + t588 * t642 + t589 * t643;
t674 = m(4) * pkin(7);
t673 = t334 / 0.2e1;
t230 = t324 * t389 - t384 * t416;
t121 = -qJD(5) * t230 + t238 * t384 - t239 * t389;
t187 = -t268 * t389 + t269 * t384;
t672 = t121 - t187;
t453 = t386 * t477;
t616 = -t370 + (-t453 + t473) * pkin(3);
t433 = mrSges(3,1) * t387 + mrSges(3,2) * t392;
t518 = Ifges(3,4) * t387;
t671 = t387 * (Ifges(3,1) * t392 - t518) / 0.2e1 - pkin(1) * t433;
t383 = qJ(3) + qJ(4);
t374 = cos(t383);
t360 = pkin(4) * t374;
t379 = t391 * pkin(3);
t337 = t360 + t379;
t376 = qJ(5) + t383;
t362 = cos(t376);
t354 = pkin(5) * t362;
t279 = t354 + t337;
t273 = pkin(2) + t279;
t328 = pkin(2) + t337;
t361 = sin(t376);
t365 = t379 + pkin(2);
t373 = sin(t383);
t432 = -mrSges(4,1) * t391 + mrSges(4,2) * t386;
t670 = m(4) * pkin(2) + m(5) * t365 + m(6) * t328 + m(7) * t273 + t374 * mrSges(5,1) - t373 * mrSges(5,2) + t644 * t361 - t362 * t645 - t432;
t382 = -pkin(10) + t394;
t375 = -qJ(6) + t382;
t669 = -m(4) * pkin(8) + m(5) * t394 + m(6) * t382 + m(7) * t375 - t676;
t667 = -t138 * mrSges(4,1) + t137 * mrSges(4,2);
t666 = -t233 * mrSges(4,1) + t234 * mrSges(4,2);
t425 = t392 * Ifges(3,2) + t518;
t662 = t117 * mrSges(5,2) + t30 * mrSges(7,2) + t47 * mrSges(6,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t425 / 0.2e1 - t116 * mrSges(5,1) - t25 * mrSges(7,1) - t46 * mrSges(6,1);
t661 = -t629 / 0.2e1;
t537 = pkin(4) * t228;
t242 = t390 * t345 + t346 * t385;
t203 = -pkin(10) * t324 + t242;
t204 = -pkin(10) * t416 + t243;
t632 = t203 * t467 - t204 * t468 + t384 * t680 - t389 * t679;
t145 = t384 * t203 + t389 * t204;
t631 = -qJD(5) * t145 + t384 * t679 + t389 * t680;
t123 = -t191 * t385 - t186;
t100 = t123 - t646;
t124 = t390 * t191 - t184;
t101 = t124 - t655;
t538 = pkin(3) * t390;
t364 = pkin(4) + t538;
t497 = t384 * t385;
t625 = -t384 * t100 - t389 * t101 + t364 * t467 + (-t385 * t468 + (t389 * t390 - t497) * qJD(4)) * pkin(3);
t496 = t385 * t389;
t623 = -t389 * t100 + t101 * t384 - t364 * t468 + (-t385 * t467 + (-t384 * t390 - t496) * qJD(4)) * pkin(3);
t618 = t682 * pkin(4) + t616;
t474 = qJD(2) * t392;
t404 = t386 * t474 + t387 * t471;
t651 = mrSges(6,1) * t361 - t362 * t644;
t388 = sin(qJ(1));
t393 = cos(qJ(1));
t650 = g(1) * t393 + g(2) * t388;
t592 = m(5) * pkin(3);
t583 = t111 / 0.2e1;
t582 = t112 / 0.2e1;
t562 = t218 / 0.2e1;
t561 = t219 / 0.2e1;
t553 = t307 / 0.2e1;
t552 = t318 / 0.2e1;
t229 = -t324 * t384 - t389 * t416;
t635 = qJ(6) * t672 + qJD(6) * t229 + t632;
t120 = qJD(5) * t229 - t238 * t389 - t239 * t384;
t188 = -t268 * t384 - t269 * t389;
t634 = -pkin(5) * t478 - qJD(6) * t230 + t631 + (-t120 + t188) * qJ(6);
t37 = -mrSges(7,2) * t291 + mrSges(7,3) * t45;
t38 = -mrSges(6,2) * t291 + mrSges(6,3) * t45;
t633 = t38 + t37;
t628 = mrSges(4,1) + t592;
t627 = -pkin(5) * t672 + t618;
t626 = t625 + t654;
t624 = t623 + t622;
t284 = t416 * t387;
t320 = t391 * t338;
t232 = -pkin(9) * t491 + t320 + (-pkin(7) * t386 - pkin(3)) * t392;
t358 = pkin(7) * t488;
t263 = t386 * t338 + t358;
t495 = t386 * t387;
t241 = -pkin(9) * t495 + t263;
t167 = t390 * t232 - t241 * t385;
t133 = -pkin(4) * t392 + pkin(10) * t284 + t167;
t168 = t385 * t232 + t390 * t241;
t283 = t324 * t387;
t146 = -pkin(10) * t283 + t168;
t72 = t384 * t133 + t389 * t146;
t617 = -qJD(2) * mrSges(3,1) - mrSges(4,1) * t321 + mrSges(4,2) * t322 + mrSges(3,3) * t478;
t487 = t392 * t393;
t260 = -t361 * t487 + t362 * t388;
t261 = t361 * t388 + t362 * t487;
t614 = t260 * t645 - t261 * t644;
t490 = t388 * t392;
t258 = t361 * t490 + t362 * t393;
t259 = t361 * t393 - t362 * t490;
t613 = -t258 * t645 + t259 * t644;
t312 = Ifges(4,4) * t321;
t207 = t322 * Ifges(4,1) + t356 * Ifges(4,5) + t312;
t368 = Ifges(3,4) * t477;
t612 = Ifges(3,1) * t478 + Ifges(3,5) * qJD(2) + t391 * t207 + t368;
t317 = t334 * pkin(7);
t611 = t316 * t392 + t317 * t387;
t610 = t137 * t391 - t138 * t386;
t608 = mrSges(5,1) * t373 + mrSges(7,1) * t361 + mrSges(5,2) * t374 + t651;
t276 = -t373 * t487 + t374 * t388;
t277 = t373 * t388 + t374 * t487;
t607 = -t276 * mrSges(5,1) + t277 * mrSges(5,2) + t614;
t274 = t373 * t490 + t374 * t393;
t275 = t373 * t393 - t374 * t490;
t606 = t274 * mrSges(5,1) - t275 * mrSges(5,2) + t613;
t605 = -mrSges(5,1) * t248 + mrSges(5,3) * t117;
t604 = mrSges(5,2) * t248 - mrSges(5,3) * t116;
t603 = -m(3) - m(7) - m(6) - m(5) - m(4);
t602 = t322 * Ifges(4,5) + t228 * Ifges(5,5) + t321 * Ifges(4,6) + Ifges(5,6) * t437 + t356 * Ifges(4,3) + t348 * Ifges(5,3) + t162 * t641 + t339 * t638 + t639 * t648;
t536 = pkin(4) * t373;
t313 = -pkin(5) * t361 - t536;
t539 = pkin(3) * t386;
t278 = -t313 + t539;
t336 = t536 + t539;
t601 = -m(6) * t336 - m(7) * t278;
t434 = t392 * mrSges(3,1) - mrSges(3,2) * t387;
t599 = t387 * t676 + mrSges(2,1) + t434;
t598 = mrSges(2,2) - mrSges(3,3) + t601;
t594 = t549 * t641 + t568 * t643 + t572 * t642 + t596;
t591 = m(6) * pkin(4);
t590 = m(7) * pkin(5);
t587 = Ifges(5,4) * t583 + Ifges(5,2) * t582 + Ifges(5,6) * t553;
t586 = Ifges(5,1) * t583 + Ifges(5,4) * t582 + Ifges(5,5) * t553;
t579 = Ifges(4,1) * t562 + Ifges(4,4) * t561 + Ifges(4,5) * t552;
t514 = Ifges(5,4) * t228;
t149 = Ifges(5,2) * t437 + t348 * Ifges(5,6) + t514;
t578 = -t149 / 0.2e1;
t577 = t149 / 0.2e1;
t222 = Ifges(5,4) * t437;
t150 = t228 * Ifges(5,1) + t348 * Ifges(5,5) + t222;
t576 = -t150 / 0.2e1;
t575 = t150 / 0.2e1;
t560 = -t437 / 0.2e1;
t559 = t437 / 0.2e1;
t558 = -t228 / 0.2e1;
t557 = t228 / 0.2e1;
t550 = t322 / 0.2e1;
t547 = -t348 / 0.2e1;
t546 = t348 / 0.2e1;
t540 = pkin(3) * t322;
t534 = pkin(4) * t389;
t530 = g(3) * t387;
t377 = t387 * pkin(7);
t49 = t389 * t92 - t87;
t524 = mrSges(5,3) * t437;
t523 = mrSges(5,3) * t228;
t522 = mrSges(6,3) * t648;
t521 = mrSges(6,3) * t162;
t520 = mrSges(7,3) * t648;
t519 = mrSges(7,3) * t162;
t517 = Ifges(3,4) * t392;
t516 = Ifges(4,4) * t386;
t515 = Ifges(4,4) * t391;
t511 = t233 * mrSges(4,3);
t510 = t234 * mrSges(4,3);
t509 = t322 * Ifges(4,4);
t494 = t386 * t388;
t492 = t386 * t393;
t332 = t435 * qJD(2);
t476 = qJD(2) * t387;
t459 = pkin(7) * t476;
t480 = t391 * t332 + t386 * t459;
t335 = pkin(3) * t495 + t377;
t472 = qJD(3) * t387;
t461 = pkin(4) * t468;
t460 = pkin(4) * t467;
t372 = pkin(7) * t474;
t456 = Ifges(4,5) * t218 + Ifges(4,6) * t219 + Ifges(4,3) * t318;
t257 = pkin(3) * t404 + t372;
t206 = t321 * Ifges(4,2) + t356 * Ifges(4,6) + t509;
t449 = -t386 * t206 / 0.2e1;
t15 = -t45 * mrSges(7,1) + t44 * mrSges(7,2);
t48 = -t384 * t92 - t89;
t71 = t389 * t133 - t146 * t384;
t144 = t389 * t203 - t204 * t384;
t270 = pkin(4) * t416 - t365;
t288 = -pkin(3) * t497 + t389 * t364;
t235 = pkin(4) * t283 + t335;
t183 = t540 + t537;
t287 = -qJDD(2) * pkin(2) + t317;
t431 = mrSges(4,1) * t386 + mrSges(4,2) * t391;
t427 = Ifges(4,1) * t391 - t516;
t426 = Ifges(4,1) * t386 + t515;
t424 = -Ifges(4,2) * t386 + t515;
t423 = Ifges(4,2) * t391 + t516;
t422 = Ifges(3,5) * t392 - Ifges(3,6) * t387;
t421 = Ifges(4,5) * t391 - Ifges(4,6) * t386;
t420 = Ifges(4,5) * t386 + Ifges(4,6) * t391;
t419 = t273 * t392 - t375 * t387;
t196 = -t283 * t389 + t284 * t384;
t197 = -t283 * t384 - t284 * t389;
t418 = t328 * t392 - t382 * t387;
t417 = t392 * t365 - t387 * t394;
t172 = -qJD(2) * t411 + t284 * t609;
t147 = -pkin(4) * t172 + t257;
t297 = -t386 * t487 + t388 * t391;
t295 = t386 * t490 + t391 * t393;
t413 = t343 * t431;
t171 = -qJD(2) * t410 - t239 * t387;
t166 = t415 * qJD(2) + (-t358 + (pkin(9) * t387 - t338) * t386) * qJD(3) + t480;
t189 = t386 * t332 + t338 * t471 + (-t387 * t475 - t392 * t473) * pkin(7);
t170 = -pkin(9) * t404 + t189;
t70 = -qJD(4) * t168 + t390 * t166 - t170 * t385;
t54 = pkin(4) * t476 - pkin(10) * t171 + t70;
t69 = t385 * t166 + t390 * t170 + t232 * t469 - t241 * t470;
t58 = pkin(10) * t172 + t69;
t9 = t133 * t467 - t146 * t468 + t384 * t54 + t389 * t58;
t180 = -pkin(3) * t219 + t287;
t405 = -t386 * t472 + t391 * t474;
t75 = -pkin(4) * t112 + t180;
t403 = Ifges(4,5) * t387 + t392 * t427;
t402 = Ifges(4,6) * t387 + t392 * t424;
t401 = Ifges(4,3) * t387 + t392 * t421;
t10 = -qJD(5) * t72 - t384 * t58 + t389 * t54;
t363 = pkin(5) + t534;
t341 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t477;
t314 = t354 + t360;
t309 = t431 * t387;
t298 = t391 * t487 + t494;
t296 = -t388 * t488 + t492;
t289 = pkin(3) * t496 + t364 * t384;
t285 = pkin(5) + t288;
t262 = -pkin(7) * t493 + t320;
t247 = mrSges(4,1) * t356 - mrSges(4,3) * t322;
t246 = -mrSges(4,2) * t356 + mrSges(4,3) * t321;
t245 = -pkin(7) * t452 + t300;
t194 = mrSges(5,1) * t348 - t523;
t193 = -mrSges(5,2) * t348 + t524;
t190 = -qJD(3) * t263 + t480;
t182 = -mrSges(4,2) * t318 + mrSges(4,3) * t219;
t181 = mrSges(4,1) * t318 - mrSges(4,3) * t218;
t179 = -pkin(5) * t229 + t270;
t165 = -mrSges(5,1) * t437 + mrSges(5,2) * t228;
t154 = -mrSges(4,1) * t219 + mrSges(4,2) * t218;
t153 = -pkin(5) * t196 + t235;
t143 = mrSges(6,1) * t339 - t521;
t142 = mrSges(7,1) * t339 - t519;
t141 = -mrSges(6,2) * t339 + t522;
t140 = -mrSges(7,2) * t339 + t520;
t134 = t218 * Ifges(4,4) + t219 * Ifges(4,2) + t318 * Ifges(4,6);
t113 = t537 + t681;
t108 = qJ(6) * t229 + t145;
t107 = -qJ(6) * t230 + t144;
t106 = t183 + t681;
t103 = -mrSges(5,2) * t307 + mrSges(5,3) * t112;
t102 = mrSges(5,1) * t307 - mrSges(5,3) * t111;
t86 = -mrSges(6,1) * t648 + mrSges(6,2) * t162;
t85 = -mrSges(7,1) * t648 + mrSges(7,2) * t162;
t74 = -qJD(5) * t197 - t171 * t384 + t172 * t389;
t73 = qJD(5) * t196 + t171 * t389 + t172 * t384;
t66 = qJ(6) * t196 + t72;
t65 = -mrSges(5,1) * t112 + mrSges(5,2) * t111;
t64 = -pkin(5) * t392 - qJ(6) * t197 + t71;
t61 = -pkin(5) * t74 + t147;
t36 = mrSges(6,1) * t291 - mrSges(6,3) * t44;
t35 = mrSges(7,1) * t291 - mrSges(7,3) * t44;
t32 = t49 - t654;
t31 = t48 - t622;
t22 = -pkin(5) * t45 + qJDD(6) + t75;
t16 = -mrSges(6,1) * t45 + mrSges(6,2) * t44;
t8 = qJ(6) * t74 + qJD(6) * t196 + t9;
t7 = pkin(5) * t476 - qJ(6) * t73 - qJD(6) * t197 + t10;
t1 = [(t602 / 0.2e1 + Ifges(5,3) * t546 + Ifges(5,5) * t557 + Ifges(5,6) * t559 + t638 * t548 + t639 * t571 + t641 * t567 - t662 - t666) * t476 + m(4) * (t137 * t263 + t138 * t262 + t189 * t234 + t190 * t233) + t517 * t673 + t677 * t74 + t678 * t73 + (-t116 * t171 + t117 * t172 - t27 * t283 + t28 * t284) * mrSges(5,3) + (t518 + t425) * t333 / 0.2e1 + (-Ifges(5,5) * t284 - Ifges(5,6) * t283) * t553 + (Ifges(5,5) * t171 + Ifges(5,6) * t172) * t546 - (t206 * t391 + t207 * t386) * t472 / 0.2e1 + (mrSges(6,2) * t75 + mrSges(7,2) * t22 - mrSges(6,3) * t6 - mrSges(7,3) * t2 + 0.2e1 * t657) * t197 + (-mrSges(6,1) * t75 - mrSges(7,1) * t22 + mrSges(6,3) * t5 + mrSges(7,3) * t3 + 0.2e1 * t656) * t196 + (-Ifges(5,1) * t284 - Ifges(5,4) * t283) * t583 + (Ifges(5,1) * t171 + Ifges(5,4) * t172) * t557 + t180 * (mrSges(5,1) * t283 - mrSges(5,2) * t284) - t134 * t495 / 0.2e1 + t671 * t466 + t356 * (qJD(2) * t401 - t420 * t472) / 0.2e1 + m(5) * (t116 * t70 + t117 * t69 + t167 * t28 + t168 * t27 + t180 * t335 + t248 * t257) + m(6) * (t10 * t46 + t147 * t176 + t235 * t75 + t47 * t9 + t5 * t72 + t6 * t71) + m(7) * (t153 * t22 + t2 * t64 + t25 * t7 + t3 * t66 + t30 * t8 + t61 * t99) + t434 * t504 + (-qJDD(2) * mrSges(3,1) + mrSges(3,3) * t334 + t154) * t377 + t321 * (qJD(2) * t402 - t423 * t472) / 0.2e1 + (-Ifges(5,4) * t284 - Ifges(5,2) * t283) * t582 + (Ifges(5,4) * t171 + Ifges(5,2) * t172) * t559 + (-t137 * t495 - t138 * t491 - t233 * t405 - t234 * t404) * mrSges(4,3) + t611 * mrSges(3,3) - t341 * t459 + Ifges(2,3) * qJDD(1) + qJD(2) ^ 2 * t422 / 0.2e1 - pkin(1) * (-mrSges(3,1) * t333 + mrSges(3,2) * t334) + t335 * t65 + t287 * t309 + t262 * t181 + t263 * t182 + t257 * t165 + t189 * t246 + t190 * t247 + t248 * (-mrSges(5,1) * t172 + mrSges(5,2) * t171) + t235 * t16 + (-t492 * t592 - t296 * mrSges(4,1) - t275 * mrSges(5,1) - t295 * mrSges(4,2) - t274 * mrSges(5,2) + t645 * t259 + t644 * t258 + (-m(7) * (-pkin(1) - t419) - m(6) * (-pkin(1) - t418) - m(5) * (-pkin(1) - t417) - m(4) * t338 + m(3) * pkin(1) + t599) * t388 + (pkin(7) * t603 + t598) * t393) * g(1) + t69 * t193 + t70 * t194 + t167 * t102 + t168 * t103 + t343 * (mrSges(4,1) * t404 + mrSges(4,2) * t405) + t147 * t86 + t153 * t15 + t8 * t140 + t9 * t141 + t7 * t142 + t10 * t143 + (-t494 * t592 - t298 * mrSges(4,1) - t277 * mrSges(5,1) - t297 * mrSges(4,2) - t276 * mrSges(5,2) + t603 * (t393 * pkin(1) + t388 * pkin(7)) + t598 * t388 + t645 * t261 + t644 * t260 + (-m(4) * t436 - m(5) * t417 - m(6) * t418 - m(7) * t419 - t599) * t393) * g(2) + t61 * t85 + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t611) + t617 * t372 + ((-Ifges(3,2) * t387 + t517) * t466 / 0.2e1 - Ifges(5,3) * t553 - Ifges(5,6) * t582 - Ifges(5,5) * t583 + (-mrSges(3,2) * pkin(7) + Ifges(3,6)) * qJDD(2) + (pkin(7) * mrSges(3,3) + Ifges(3,2) / 0.2e1) * t333 + t665 - t457 / 0.2e1 - t615 / 0.2e1 - t456 / 0.2e1 - Ifges(4,3) * t552 - Ifges(4,6) * t561 - Ifges(4,5) * t562 - t638 * t554 - t639 * t588 - t641 * t589 + Ifges(3,4) * t673 + t667 + t668) * t392 + (t343 * t674 + t449 + t612 / 0.2e1) * t474 + (t334 * Ifges(3,1) + Ifges(3,5) * qJDD(2) + t287 * t674 + t421 * t552 + t424 * t561 + t427 * t562) * t387 + (qJD(2) * t403 - t426 * t472) * t550 + t171 * t575 + t172 * t577 + t491 * t579 - t284 * t586 - t283 * t587 + t64 * t35 + t66 * t37 + t71 * t36 + t72 * t38; t650 * t433 - t434 * g(3) + (t449 + t413) * qJD(3) + (-g(3) * t670 + t650 * t669) * t392 + (Ifges(5,5) * t558 + Ifges(5,6) * t560 + Ifges(5,3) * t547 + t638 * t549 + t641 * t568 + t639 * t572 + t662) * t478 + t180 * (mrSges(5,1) * t416 + mrSges(5,2) * t324) + (Ifges(5,5) * t324 - Ifges(5,6) * t416) * t553 + t677 * t121 + t678 * t120 + (-t187 * t30 + t188 * t25 - t2 * t230 + t229 * t3) * mrSges(7,3) + (t321 * t424 + t322 * t427 + t356 * t421) * qJD(3) / 0.2e1 - (t321 * t402 + t322 * t403 + t356 * t401) * qJD(1) / 0.2e1 - t602 * t478 / 0.2e1 + (-pkin(2) * t287 - t233 * t244 - t234 * t245 - t343 * t370) * m(4) + ((t233 * t488 + t234 * t493) * qJD(1) + t610) * mrSges(4,3) + t206 * t453 / 0.2e1 - t473 * t510 + (-t116 * t269 + t117 * t268 - t27 * t416 - t28 * t324) * mrSges(5,3) + (-Ifges(5,4) * t269 - Ifges(5,2) * t268) * t560 - t248 * (mrSges(5,1) * t268 - mrSges(5,2) * t269) + (-Ifges(5,5) * t269 - Ifges(5,6) * t268) * t547 - t422 * t466 / 0.2e1 - t671 * qJD(1) ^ 2 + (Ifges(5,4) * t324 - Ifges(5,2) * t416) * t582 + (Ifges(5,1) * t324 - Ifges(5,4) * t416) * t583 - t416 * t587 + (-Ifges(5,1) * t269 - Ifges(5,4) * t268) * t558 + (g(3) * t669 + qJD(1) * t666 + t650 * t670) * t387 + (-t187 * t47 + t188 * t46 + t229 * t5 - t230 * t6) * mrSges(6,3) + t229 * t656 + t230 * t657 + t187 * t659 + t188 * t661 - t413 * t477 + t341 * t369 - (Ifges(5,1) * t557 + Ifges(5,4) * t559 + Ifges(5,5) * t546 + t575 + t604) * t238 - (Ifges(5,4) * t557 + Ifges(5,2) * t559 + Ifges(5,6) * t546 + t577 + t605) * t239 + t287 * t432 + t391 * t134 / 0.2e1 - t365 * t65 + Ifges(3,6) * t333 + Ifges(3,5) * t334 + Ifges(3,3) * qJDD(2) - t316 * mrSges(3,2) - t317 * mrSges(3,1) + t270 * t16 + t242 * t102 + t243 * t103 - t245 * t246 - t244 * t247 + t22 * (-mrSges(7,1) * t229 + mrSges(7,2) * t230) + t75 * (-mrSges(6,1) * t229 + mrSges(6,2) * t230) + (-t511 + t207 / 0.2e1) * t471 - t99 * (-mrSges(7,1) * t187 + mrSges(7,2) * t188) - t176 * (-mrSges(6,1) * t187 + mrSges(6,2) * t188) + t179 * t15 + t145 * t38 - pkin(2) * t154 + t144 * t36 + t107 * t35 + t108 * t37 + (-t386 * t181 + m(4) * ((-t233 * t391 - t234 * t386) * qJD(3) + t610) - t247 * t471 - t246 * t473 + t391 * t182) * pkin(8) - (-Ifges(3,2) * t478 + t368 + t612) * t477 / 0.2e1 + t616 * t165 - t617 * t370 + t618 * t86 + t619 * t194 + t620 * t193 + (t116 * t619 + t117 * t620 - t180 * t365 + t242 * t28 + t243 * t27 + t248 * t616) * m(5) + t627 * t85 + t420 * t552 + t423 * t561 + t426 * t562 - t269 * t576 - t268 * t578 + t386 * t579 + t324 * t586 + t631 * t143 + t632 * t141 + (t144 * t6 + t145 * t5 + t176 * t618 + t270 * t75 + t46 * t631 + t47 * t632) * m(6) + t634 * t142 + t635 * t140 + (t107 * t2 + t108 * t3 + t179 * t22 + t25 * t634 + t30 * t635 + t627 * t99) * m(7) + (t229 * t639 + t230 * t641) * t554 + (t187 * t639 + t188 * t641) * t549 + (t187 * t640 + t188 * t642) * t572 + (t229 * t640 + t230 * t642) * t588 + (t229 * t642 + t230 * t643) * t589 + (t187 * t642 + t188 * t643) * t568; -(Ifges(5,4) * t558 + Ifges(5,2) * t560 + Ifges(5,6) * t547 + t578 - t605) * t228 + (t594 + t661) * t648 - (-Ifges(4,2) * t322 + t207 + t312) * t321 / 0.2e1 + t456 - t322 * (Ifges(4,1) * t321 - t509) / 0.2e1 + t675 - t165 * t540 - m(5) * (t116 * t123 + t117 * t124 + t248 * t540) + (Ifges(5,1) * t558 + Ifges(5,4) * t560 + Ifges(5,5) * t547 + t576 - t604) * t437 + (t103 * t385 + t193 * t469 - t194 * t470) * pkin(3) - t356 * (Ifges(4,5) * t321 - Ifges(4,6) * t322) / 0.2e1 - t343 * (mrSges(4,1) * t322 + mrSges(4,2) * t321) + g(3) * t309 + t288 * t36 + t285 * t35 - t233 * t246 + t234 * t247 - t124 * t193 - t123 * t194 - t183 * t86 - t106 * t85 + t322 * t510 + t321 * t511 + (m(5) * t539 - t601 + t608) * t530 + t623 * t143 + t624 * t142 + t625 * t141 + (-t176 * t183 + t288 * t6 + t289 * t5 + t46 * t623 + t47 * t625) * m(6) + t626 * t140 + (-t106 * t99 + t2 * t285 + t25 * t624 + t289 * t3 + t30 * t626) * m(7) + (-m(6) * (-t336 * t490 - t337 * t393) - m(7) * (-t278 * t490 - t279 * t393) - mrSges(4,2) * t296 + t628 * t295 + t606) * g(2) + (mrSges(4,2) * t298 - m(6) * (-t336 * t487 + t337 * t388) - m(7) * (-t278 * t487 + t279 * t388) - t628 * t297 + t607) * g(1) + t102 * t538 + t206 * t550 + (t27 * t385 + t28 * t390 + (-t116 * t385 + t117 * t390) * qJD(4)) * t592 + t633 * t289 - t667; t594 * t648 - t248 * (mrSges(5,1) * t228 + mrSges(5,2) * t437) + (Ifges(5,5) * t437 - Ifges(5,6) * t228) * t547 + (Ifges(5,1) * t437 - t514) * t558 + t633 * pkin(4) * t384 + (-t193 + t524) * t116 + (-Ifges(5,2) * t228 + t150 + t222) * t560 + (t2 * t363 + (t3 * t384 + (-t25 * t384 + t30 * t389) * qJD(5)) * pkin(4) - t113 * t99 - t25 * t31 - t30 * t32) * m(7) - t86 * t537 - m(6) * (t176 * t537 + t46 * t48 + t47 * t49) + (t194 + t523) * t117 + (-t461 - t48) * t143 + (-t49 + t460) * t141 + (-t32 + t460) * t140 + t675 + (-t461 - t31) * t142 + t363 * t35 - t113 * t85 + t629 * t572 + (-m(7) * (t313 * t490 - t314 * t393) + t274 * t591 + t606) * g(2) + (-t276 * t591 - m(7) * (t313 * t487 + t314 * t388) + t607) * g(1) + (m(6) * t536 - m(7) * t313 + t608) * t530 + t36 * t534 + t149 * t557 + (t384 * t5 + t389 * t6 + (-t384 * t46 + t389 * t47) * qJD(5)) * t591; -t176 * (mrSges(6,1) * t162 + mrSges(6,2) * t648) - t99 * (mrSges(7,1) * t162 + mrSges(7,2) * t648) - t29 * t140 + t25 * t520 + t2 * t590 + t400 + (t643 * t648 - t652) * t568 + t630 * t567 + (-t162 * t639 + t641 * t648) * t549 + (-(-mrSges(7,1) - t590) * t361 + t651) * t530 + (t143 + t521) * t47 + (-t141 + t522) * t46 + (-m(7) * (-t25 + t29) + t142 + t519) * t30 + (t258 * t590 + t613) * g(2) + (-t260 * t590 + t614) * g(1) + (-t162 * t640 + t629 + t653) * t572 + (t35 + (-m(7) * t99 - t85) * t162) * pkin(5); -t648 * t140 + t162 * t142 + (g(3) * t392 + t162 * t25 - t30 * t648 - t387 * t650 + t22) * m(7) + t15;];
tau  = t1;
