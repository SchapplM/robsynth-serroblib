% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:22:04
% EndTime: 2019-03-09 17:22:57
% DurationCPUTime: 37.24s
% Computational Cost: add. (9810->839), mult. (20785->1020), div. (0->0), fcn. (13726->8), ass. (0->383)
t668 = Ifges(4,1) + Ifges(5,1);
t667 = Ifges(5,4) + Ifges(4,5);
t666 = Ifges(4,6) - Ifges(5,6);
t331 = cos(qJ(2));
t464 = qJD(1) * t331;
t288 = qJD(3) - t464;
t279 = qJD(5) - t288;
t326 = sin(qJ(5));
t527 = cos(qJ(5));
t333 = -pkin(3) - pkin(4);
t328 = sin(qJ(2));
t468 = t331 * pkin(2) + t328 * pkin(8);
t625 = -pkin(1) - t468;
t225 = t625 * qJD(1);
t308 = pkin(7) * t464;
t268 = qJD(2) * pkin(8) + t308;
t327 = sin(qJ(3));
t330 = cos(qJ(3));
t147 = t330 * t225 - t327 * t268;
t465 = qJD(1) * t328;
t437 = t330 * t465;
t243 = qJD(2) * t327 + t437;
t96 = pkin(9) * t243 + t147;
t663 = qJD(4) - t96;
t77 = t288 * t333 + t663;
t272 = t288 * qJ(4);
t148 = t327 * t225 + t330 * t268;
t439 = t327 * t465;
t456 = t330 * qJD(2);
t242 = t439 - t456;
t399 = pkin(9) * t242 + t148;
t84 = t272 + t399;
t28 = -t326 * t84 + t527 * t77;
t604 = -t28 + qJD(6);
t22 = -t279 * pkin(5) + t604;
t145 = t326 * t242 + t243 * t527;
t361 = t242 * t527 - t326 * t243;
t307 = pkin(7) * t465;
t267 = -qJD(2) * pkin(2) + t307;
t111 = t242 * pkin(3) - t243 * qJ(4) + t267;
t90 = -pkin(4) * t242 - t111;
t30 = -pkin(5) * t361 - qJ(6) * t145 + t90;
t140 = Ifges(6,4) * t361;
t504 = Ifges(7,5) * t361;
t628 = Ifges(7,4) + Ifges(6,5);
t631 = Ifges(6,1) + Ifges(7,1);
t618 = t145 * t631 + t279 * t628 + t140 - t504;
t677 = t90 * mrSges(6,2) + mrSges(7,2) * t22 - mrSges(6,3) * t28 - t30 * mrSges(7,3) + t618 / 0.2e1;
t676 = t288 / 0.2e1;
t533 = t279 / 0.2e1;
t545 = t145 / 0.2e1;
t548 = -t361 / 0.2e1;
t549 = t361 / 0.2e1;
t675 = Ifges(6,4) * t549 + Ifges(7,5) * t548 + t628 * t533 + t631 * t545 + t677;
t534 = -t279 / 0.2e1;
t546 = -t145 / 0.2e1;
t674 = -Ifges(6,4) * t548 - Ifges(7,5) * t549 - t628 * t534 - t631 * t546 + t677;
t626 = Ifges(6,3) + Ifges(7,2);
t673 = t626 * t533 + t628 * t545;
t455 = qJD(1) * qJD(2);
t255 = qJDD(1) * t328 + t331 * t455;
t461 = qJD(3) * t242;
t136 = qJDD(2) * t327 + t255 * t330 - t461;
t137 = qJD(3) * t243 - t330 * qJDD(2) + t255 * t327;
t37 = qJD(5) * t361 + t136 * t527 + t326 * t137;
t672 = -t37 / 0.2e1;
t38 = qJD(5) * t145 + t326 * t136 - t137 * t527;
t560 = -t38 / 0.2e1;
t671 = m(7) + m(6);
t552 = t136 / 0.2e1;
t550 = t137 / 0.2e1;
t254 = t331 * qJDD(1) - t328 * t455;
t240 = qJDD(3) - t254;
t224 = qJDD(5) - t240;
t670 = -t224 / 0.2e1;
t540 = t240 / 0.2e1;
t538 = t242 / 0.2e1;
t669 = -t243 / 0.2e1;
t629 = -Ifges(6,4) + Ifges(7,5);
t665 = Ifges(6,6) - Ifges(7,6);
t664 = -Ifges(4,3) - Ifges(5,2);
t600 = -t327 * t666 + t330 * t667;
t506 = Ifges(5,5) * t327;
t509 = Ifges(4,4) * t327;
t597 = t330 * t668 + t506 - t509;
t388 = t330 * mrSges(5,1) + t327 * mrSges(5,3);
t390 = mrSges(4,1) * t330 - mrSges(4,2) * t327;
t662 = -t390 - t388;
t436 = t330 * t464;
t458 = qJD(3) * t330;
t661 = -t327 * qJD(4) - t308 + (t436 - t458) * qJ(4);
t583 = qJD(4) - t147;
t105 = -pkin(3) * t288 + t583;
t585 = t105 * mrSges(5,2) - t147 * mrSges(4,3);
t660 = t267 * mrSges(4,2) - t111 * mrSges(5,3) + t585;
t107 = t272 + t148;
t586 = -t107 * mrSges(5,2) - t148 * mrSges(4,3);
t659 = t267 * mrSges(4,1) + t111 * mrSges(5,1) + t586;
t29 = t326 * t77 + t527 * t84;
t23 = t279 * qJ(6) + t29;
t139 = Ifges(7,5) * t145;
t56 = t279 * Ifges(7,6) - Ifges(7,3) * t361 + t139;
t507 = Ifges(6,4) * t145;
t59 = Ifges(6,2) * t361 + t279 * Ifges(6,6) + t507;
t658 = t56 / 0.2e1 + mrSges(6,1) * t90 + mrSges(7,1) * t30 - t59 / 0.2e1 - mrSges(7,2) * t23 - mrSges(6,3) * t29;
t657 = mrSges(6,3) + mrSges(7,2) - mrSges(4,3) - mrSges(5,2);
t551 = -t137 / 0.2e1;
t635 = t254 / 0.2e1;
t656 = t255 / 0.2e1;
t654 = t667 * t540 + (-Ifges(4,4) + Ifges(5,5)) * t550 + t668 * t552;
t443 = -pkin(7) * t327 - pkin(3);
t475 = t330 * t331;
t346 = -pkin(9) * t475 + (-pkin(4) + t443) * t328;
t397 = pkin(2) * t328 - pkin(8) * t331;
t251 = t397 * qJD(1);
t483 = t251 * t330;
t102 = qJD(1) * t346 - t483;
t218 = t327 * t251;
t298 = qJ(4) * t465;
t478 = t328 * t330;
t479 = t327 * t331;
t110 = t218 + t298 + (-pkin(7) * t478 + pkin(9) * t479) * qJD(1);
t553 = pkin(8) - pkin(9);
t269 = t553 * t327;
t270 = t553 * t330;
t160 = t326 * t269 + t270 * t527;
t460 = qJD(3) * t327;
t252 = t553 * t460;
t430 = qJD(3) * t527;
t404 = t330 * t430;
t619 = -qJD(5) * t160 - t102 * t527 + t404 * t553 + (t110 + t252) * t326;
t100 = -mrSges(7,1) * t279 + mrSges(7,2) * t145;
t512 = mrSges(6,3) * t145;
t99 = mrSges(6,1) * t279 - t512;
t653 = t99 - t100;
t234 = Ifges(4,4) * t242;
t499 = t242 * Ifges(5,5);
t652 = t243 * t668 + t288 * t667 - t234 + t499;
t438 = t327 * t464;
t446 = t333 * t327;
t608 = qJD(3) * t446 - t333 * t438 - t661;
t329 = sin(qJ(1));
t332 = cos(qJ(1));
t473 = t332 * t327;
t216 = -t329 * t330 + t331 * t473;
t474 = t331 * t332;
t217 = t327 * t329 + t330 * t474;
t127 = t216 * t326 + t217 * t527;
t363 = -t216 * t527 + t217 * t326;
t406 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t409 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t650 = -t127 * t406 + t363 * t409;
t649 = -Ifges(6,2) * t548 + Ifges(7,3) * t549 - t534 * t665 + t546 * t629 - t658;
t476 = t329 * t331;
t214 = t327 * t476 + t330 * t332;
t215 = t329 * t475 - t473;
t364 = t214 * t326 + t215 * t527;
t607 = -t214 * t527 + t215 * t326;
t647 = -t364 * t406 + t409 * t607;
t72 = pkin(5) * t145 - qJ(6) * t361;
t487 = qJDD(1) * pkin(1);
t150 = -pkin(2) * t254 - pkin(8) * t255 - t487;
t238 = t254 * pkin(7);
t207 = qJDD(2) * pkin(8) + t238;
t54 = t150 * t330 - t327 * t207 - t225 * t460 - t268 * t458;
t354 = qJDD(4) - t54;
t18 = -pkin(9) * t136 + t240 * t333 + t354;
t53 = t327 * t150 + t330 * t207 + t225 * t458 - t268 * t460;
t39 = t240 * qJ(4) + t288 * qJD(4) + t53;
t20 = pkin(9) * t137 + t39;
t4 = -qJD(5) * t29 + t18 * t527 - t326 * t20;
t2 = -t224 * pkin(5) + qJDD(6) - t4;
t541 = t224 / 0.2e1;
t559 = t38 / 0.2e1;
t561 = t37 / 0.2e1;
t643 = -mrSges(7,2) * t2 + mrSges(6,3) * t4 - Ifges(7,5) * t559 + (-t561 + t672) * t631 + (-t541 + t670) * t628 + (-Ifges(6,4) + t629) * t560;
t642 = -Ifges(6,2) * t549 + Ifges(7,3) * t548 - t533 * t665 + t545 * t629 + t658;
t511 = Ifges(3,4) * t328;
t616 = t331 * Ifges(3,2);
t382 = t511 + t616;
t638 = t23 * mrSges(7,3) + t28 * mrSges(6,1) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t382 / 0.2e1 + t665 * t549 + t664 * t676 + t667 * t669 + t666 * t538 - t22 * mrSges(7,1) - t29 * mrSges(6,2) + t673;
t462 = qJD(2) * t331;
t449 = pkin(7) * t462;
t634 = mrSges(5,1) + mrSges(4,1);
t633 = -mrSges(3,3) + mrSges(2,2);
t632 = -mrSges(5,3) + mrSges(4,2);
t24 = mrSges(6,1) * t224 - mrSges(6,3) * t37;
t25 = -t224 * mrSges(7,1) + t37 * mrSges(7,2);
t624 = t25 - t24;
t26 = -mrSges(6,2) * t224 - mrSges(6,3) * t38;
t27 = -mrSges(7,2) * t38 + mrSges(7,3) * t224;
t623 = t26 + t27;
t440 = t527 * t330;
t359 = -t326 * t327 - t440;
t582 = -qJD(3) + qJD(5);
t151 = t582 * t359;
t429 = qJD(5) * t527;
t433 = t326 * t458;
t457 = qJD(5) * t326;
t152 = -t330 * t457 + t433 + (t429 - t430) * t327;
t441 = t527 * t327;
t405 = t331 * t441;
t177 = -qJD(1) * t405 + t326 * t436;
t197 = t359 * t331;
t178 = qJD(1) * t197;
t245 = -t326 * t330 + t441;
t622 = -qJD(6) * t245 + t608 + (-t151 - t178) * qJ(6) + (t152 - t177) * pkin(5);
t621 = -pkin(5) * t465 - t619;
t97 = mrSges(7,2) * t361 + mrSges(7,3) * t279;
t513 = mrSges(6,3) * t361;
t98 = -mrSges(6,2) * t279 + t513;
t617 = -t98 - t97;
t392 = mrSges(3,1) * t331 - mrSges(3,2) * t328;
t615 = -t392 - mrSges(2,1);
t349 = t328 * t245;
t289 = pkin(7) * t479;
t320 = t331 * pkin(3);
t518 = pkin(9) * t328;
t135 = pkin(4) * t331 + t289 + t320 + (-t625 - t518) * t330;
t291 = pkin(7) * t475;
t176 = t327 * t625 + t291;
t161 = -qJ(4) * t331 + t176;
t480 = t327 * t328;
t146 = pkin(9) * t480 + t161;
t609 = t326 * t135 + t527 * t146;
t606 = (-t438 + t460) * pkin(3) + t661;
t605 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t242 + mrSges(4,2) * t243 + mrSges(3,3) * t465;
t257 = t527 * qJ(4) + t326 * t333;
t603 = -t328 * t664 + t331 * t600;
t602 = t328 * t667 + t331 * t597;
t491 = t330 * mrSges(5,3);
t387 = t327 * mrSges(5,1) - t491;
t389 = mrSges(4,1) * t327 + mrSges(4,2) * t330;
t601 = -t111 * t387 - t267 * t389;
t599 = t327 * t667 + t330 * t666;
t505 = Ifges(5,5) * t330;
t508 = Ifges(4,4) * t330;
t598 = t327 * t668 - t505 + t508;
t459 = qJD(3) * t328;
t434 = t327 * t459;
t435 = t331 * t456;
t352 = -t434 + t435;
t594 = t136 * t667 - t137 * t666 - t240 * t664;
t593 = t224 * t626 + t37 * t628 - t38 * t665;
t233 = Ifges(5,5) * t243;
t112 = t288 * Ifges(5,6) + t242 * Ifges(5,3) + t233;
t305 = Ifges(3,4) * t464;
t592 = Ifges(3,1) * t465 + Ifges(3,5) * qJD(2) + t327 * t112 + t305;
t239 = t255 * pkin(7);
t590 = t238 * t331 + t239 * t328;
t194 = -t326 * t480 - t328 * t440;
t195 = t326 * t478 - t328 * t441;
t589 = pkin(5) * t195 + qJ(6) * t194;
t588 = -t327 * t54 + t330 * t53;
t41 = -pkin(3) * t240 + t354;
t587 = t327 * t41 + t330 * t39;
t314 = t327 * qJ(4);
t258 = -t330 * pkin(3) - pkin(2) - t314;
t580 = -t136 * Ifges(5,5) / 0.2e1 - t240 * Ifges(5,6) / 0.2e1 + Ifges(4,4) * t552 + Ifges(4,6) * t540 + (Ifges(5,3) + Ifges(4,2)) * t551;
t579 = -m(4) - m(5) - t671;
t578 = -g(1) * t216 - g(2) * t214 - g(3) * t480;
t253 = t397 * qJD(2);
t393 = qJD(3) * t291 - t253 * t330 + t460 * t625;
t70 = pkin(9) * t434 + qJD(2) * t346 + t393;
t463 = qJD(2) * t328;
t300 = qJ(4) * t463;
t472 = t327 * t253 + t458 * t625;
t71 = t300 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t478 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t327) * t331 + t472;
t13 = -qJD(5) * t609 - t326 * t71 + t527 * t70;
t577 = t657 * t328;
t3 = t326 * t18 + t527 * t20 + t77 * t429 - t457 * t84;
t1 = qJ(6) * t224 + qJD(6) * t279 + t3;
t570 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t569 = -t54 * mrSges(4,1) + t41 * mrSges(5,1) + t53 * mrSges(4,2) - t39 * mrSges(5,3);
t566 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 + 0.2e1 * Ifges(7,3) * t559 + Ifges(6,4) * t672 + Ifges(6,6) * t670 + (t629 + Ifges(7,5)) * t561 + (-t665 + Ifges(7,6)) * t541 + (-t560 + t559) * Ifges(6,2);
t539 = -t242 / 0.2e1;
t536 = t243 / 0.2e1;
t532 = -t288 / 0.2e1;
t521 = pkin(7) * t328;
t515 = -qJD(1) / 0.2e1;
t514 = qJD(3) / 0.2e1;
t510 = Ifges(3,4) * t331;
t498 = t243 * Ifges(4,4);
t52 = t326 * t102 + t527 * t110;
t477 = t328 * t332;
t155 = t243 * pkin(3) + t242 * qJ(4);
t470 = qJ(4) * t435 + qJD(4) * t478;
t467 = t332 * pkin(1) + t329 * pkin(7);
t450 = pkin(7) * t463;
t448 = pkin(8) * t460;
t447 = pkin(8) * t458;
t208 = -qJDD(2) * pkin(2) + t239;
t115 = -t242 * Ifges(4,2) + t288 * Ifges(4,6) + t498;
t432 = -t327 * t115 / 0.2e1;
t422 = -t464 / 0.2e1;
t419 = t462 / 0.2e1;
t416 = -t459 / 0.2e1;
t415 = t458 / 0.2e1;
t414 = t455 / 0.2e1;
t86 = -t240 * mrSges(5,1) + t136 * mrSges(5,2);
t413 = t195 * mrSges(6,1) - t194 * mrSges(6,2);
t412 = t195 * mrSges(7,1) + t194 * mrSges(7,3);
t411 = -t214 * pkin(3) + qJ(4) * t215;
t410 = -t216 * pkin(3) + qJ(4) * t217;
t175 = t330 * t625 - t289;
t232 = t330 * pkin(4) - t258;
t408 = pkin(3) * t475 + qJ(4) * t479 + t468;
t407 = pkin(2) * t474 + pkin(8) * t477 + t467;
t400 = -pkin(7) + t446;
t101 = -pkin(4) * t243 - t155;
t398 = t443 * t328;
t322 = t332 * pkin(7);
t395 = -t215 * pkin(3) - qJ(4) * t214 + t322;
t391 = mrSges(3,1) * t328 + mrSges(3,2) * t331;
t381 = -Ifges(4,2) * t327 + t508;
t380 = Ifges(4,2) * t330 + t509;
t377 = Ifges(3,5) * t331 - Ifges(3,6) * t328;
t374 = Ifges(5,3) * t327 + t505;
t373 = -Ifges(5,3) * t330 + t506;
t164 = -pkin(7) * t437 + t218;
t284 = qJ(4) * t478;
t371 = t328 * t446 + t284;
t187 = -qJ(4) * t457 + t527 * qJD(4) + t333 * t429;
t40 = t137 * pkin(3) - t136 * qJ(4) - t243 * qJD(4) + t208;
t366 = pkin(1) * t391;
t365 = t330 * t333 - t314;
t256 = -t326 * qJ(4) + t333 * t527;
t66 = t135 * t527 - t326 * t146;
t360 = t269 * t527 - t326 * t270;
t356 = t328 * (Ifges(3,1) * t331 - t511);
t12 = t135 * t429 - t146 * t457 + t326 * t70 + t527 * t71;
t158 = t328 * t400 + t284;
t351 = t327 * t462 + t328 * t458;
t350 = t217 * pkin(3) + qJ(4) * t216 + t407;
t21 = -pkin(4) * t137 - t40;
t342 = Ifges(4,6) * t328 + t331 * t381;
t341 = Ifges(5,6) * t328 + t331 * t374;
t93 = (-t328 * t456 - t331 * t460) * pkin(7) + t472;
t337 = t570 + t593;
t78 = t365 * t459 + t400 * t462 + t470;
t265 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t464;
t250 = pkin(5) - t256;
t249 = -qJ(6) + t257;
t227 = t389 * t328;
t186 = -t284 + (pkin(3) * t327 + pkin(7)) * t328;
t179 = -qJD(6) + t187;
t168 = -mrSges(5,2) * t242 + mrSges(5,3) * t288;
t167 = -mrSges(5,1) * t288 + mrSges(5,2) * t243;
t166 = mrSges(4,1) * t288 - mrSges(4,3) * t243;
t165 = -mrSges(4,2) * t288 - mrSges(4,3) * t242;
t163 = pkin(7) * t439 + t483;
t162 = -t175 + t320;
t156 = mrSges(5,1) * t242 - mrSges(5,3) * t243;
t154 = qJD(1) * t398 - t483;
t153 = t164 + t298;
t103 = -pkin(5) * t359 - qJ(6) * t245 + t232;
t94 = t327 * t450 - t393;
t92 = pkin(3) * t351 + qJ(4) * t434 + t449 - t470;
t89 = qJD(2) * t398 + t393;
t88 = -mrSges(5,2) * t137 + mrSges(5,3) * t240;
t87 = -mrSges(4,2) * t240 - mrSges(4,3) * t137;
t85 = mrSges(4,1) * t240 - mrSges(4,3) * t136;
t83 = -qJD(4) * t331 + t300 + t93;
t81 = qJD(5) * t360 - t252 * t527 + t433 * t553;
t80 = -qJD(2) * t197 + t349 * t582;
t79 = -qJD(2) * t405 + (-qJD(5) * t359 - t404) * t328 + t352 * t326;
t76 = t158 + t589;
t74 = -mrSges(6,1) * t361 + mrSges(6,2) * t145;
t73 = -mrSges(7,1) * t361 - mrSges(7,3) * t145;
t69 = mrSges(4,1) * t137 + mrSges(4,2) * t136;
t68 = mrSges(5,1) * t137 - mrSges(5,3) * t136;
t62 = -t331 * pkin(5) - t66;
t55 = qJ(6) * t331 + t609;
t45 = -qJ(6) * t465 + t52;
t44 = t326 * t399 + t527 * t96;
t36 = t101 - t72;
t16 = pkin(5) * t79 - qJ(6) * t80 + qJD(6) * t194 + t78;
t15 = mrSges(6,1) * t38 + mrSges(6,2) * t37;
t14 = mrSges(7,1) * t38 - mrSges(7,3) * t37;
t11 = pkin(5) * t463 - t13;
t10 = -qJ(6) * t463 + qJD(6) * t331 + t12;
t5 = pkin(5) * t38 - qJ(6) * t37 - qJD(6) * t145 + t21;
t6 = [t382 * t635 + (Ifges(3,4) * t656 + Ifges(3,2) * t635 + Ifges(3,6) * qJDD(2) + t593 / 0.2e1 - t594 / 0.2e1 - Ifges(5,6) * t550 - Ifges(4,6) * t551 + t510 * t414 + pkin(7) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t254) + Ifges(7,6) * t559 + Ifges(6,6) * t560 + t628 * t561 - t667 * t552 + t626 * t541 + t664 * t540 + t569 + t570) * t331 + (-m(5) * t395 - t671 * (-t215 * pkin(4) + t329 * t518 + t395) + t633 * t332 + (-m(4) - m(3)) * t322 + t634 * t215 - t632 * t214 + t409 * t364 + t406 * t607 + (m(3) * pkin(1) + t579 * t625 - t577 - t615) * t329) * g(1) + (-m(3) * t467 - m(4) * t407 - m(5) * t350 - t671 * (t217 * pkin(4) - pkin(9) * t477 + t350) + t615 * t332 + t633 * t329 - t634 * t217 + t632 * t216 - t409 * t127 - t406 * t363 + t657 * t477) * g(2) + qJD(2) ^ 2 * t377 / 0.2e1 + (qJD(2) * t603 - t459 * t599) * t676 + t592 * t419 + (qJD(2) * t602 - t459 * t598) * t536 + t605 * t449 + m(6) * (t12 * t29 + t13 * t28 + t158 * t21 + t3 * t609 + t4 * t66 + t78 * t90) + t609 * t26 + t392 * t487 + (-qJDD(2) * mrSges(3,1) + t69) * t521 + m(7) * (t1 * t55 + t10 * t23 + t11 * t22 + t16 * t30 + t2 * t62 + t5 * t76) + m(5) * (t105 * t89 + t107 * t83 + t111 * t92 + t161 * t39 + t162 * t41 + t186 * t40) + t356 * t414 + (t147 * mrSges(4,1) - t105 * mrSges(5,1) - t148 * mrSges(4,2) + t107 * mrSges(5,3) - Ifges(6,6) * t549 - Ifges(7,6) * t548 - t638 - t673) * t463 + (qJD(2) * t341 - t373 * t459) * t538 + (qJD(2) * t342 - t380 * t459) * t539 + t675 * t80 + t432 * t462 + t642 * t79 + t643 * t194 + (-t39 * mrSges(5,2) - t53 * mrSges(4,3) - t580) * t480 + (t255 * t521 + t590) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t590) - t366 * t455 + t5 * t412 + t21 * t413 - t265 * t450 + t566 * t195 + Ifges(2,3) * qJDD(1) + t510 * t656 - pkin(1) * (-mrSges(3,1) * t254 + mrSges(3,2) * t255) + t208 * t227 + t55 * t27 + t62 * t25 + t66 * t24 + m(4) * (t147 * t94 + t148 * t93 + t175 * t54 + t176 * t53 + t267 * t449) + t16 * t73 + t330 * t115 * t416 + t76 * t14 + (m(4) * t208 * pkin(7) + Ifges(3,1) * t255 + Ifges(3,4) * t635 + Ifges(3,5) * qJDD(2) + t112 * t415 + t374 * t550 + t381 * t551 + t40 * t387 - t414 * t616 + t540 * t600 + t552 * t597) * t328 + t78 * t74 + t652 * (t327 * t416 + t330 * t419) + t10 * t97 + t12 * t98 + t13 * t99 + t11 * t100 + t659 * t351 + t660 * t352 + t92 * t156 + t158 * t15 + (mrSges(5,2) * t41 - mrSges(4,3) * t54 + t654) * t478 + t161 * t88 + t162 * t86 + t93 * t165 + t94 * t166 + t89 * t167 + t83 * t168 + t175 * t85 + t176 * t87 + t186 * t68; ((g(1) * t474 + g(2) * t476) * t579 + ((t105 * t330 - t107 * t327) * qJD(3) + t587) * m(5) + (t86 - t85) * t327 + (t88 + t87) * t330 + ((-t147 * t330 - t148 * t327) * qJD(3) + t588) * m(4)) * pkin(8) + (-t448 - t153) * t168 + (-m(4) * t468 - m(5) * t408 - t392 - t671 * (pkin(4) * t475 + t408 - t518) + t662 * t331 + t409 * t197 - t406 * (t326 * t475 - t405) + t577) * g(3) + t598 * t552 + t599 * t540 + (t432 - t601) * qJD(3) + t601 * t464 - t605 * t308 + t606 * t156 + t608 * t74 + (t360 * t4 + t160 * t3 + t21 * t232 + t608 * t90 + (-t52 + t81) * t29 + t619 * t28) * m(6) + (t1 * t160 + t103 * t5 - t360 * t2 + t622 * t30 + (-t45 + t81) * t23 + t621 * t22) * m(7) - t624 * t360 - (mrSges(6,1) * t21 + mrSges(7,1) * t5 + t566) * t359 + (t447 - t154) * t167 + (-t381 / 0.2e1 + t374 / 0.2e1) * t461 + (-t447 - t163) * t166 + t674 * t178 + t675 * t151 + (t305 + t592) * t422 - t617 * t81 + t619 * t99 + (-t448 - t164) * t165 + (-Ifges(3,2) * t422 - Ifges(6,6) * t548 - Ifges(7,6) * t549 - t626 * t534 - t628 * t546 + t638) * t465 + ((t342 / 0.2e1 - t341 / 0.2e1) * t242 - t147 * (mrSges(4,1) * t328 - mrSges(4,3) * t475) - t105 * (-mrSges(5,1) * t328 + mrSges(5,2) * t475) - t148 * (-mrSges(4,2) * t328 - mrSges(4,3) * t479) - t107 * (-mrSges(5,2) * t479 + mrSges(5,3) * t328) + (t366 - t356 / 0.2e1) * qJD(1)) * qJD(1) + t373 * t550 + t642 * t152 + (mrSges(6,2) * t21 - mrSges(7,3) * t5 - t643) * t245 + t580 * t330 + (-t105 * t154 - t107 * t153 + t606 * t111 + t258 * t40) * m(5) + (-pkin(2) * t208 - t147 * t163 - t148 * t164 - t267 * t308) * m(4) + t585 * t458 + (t112 / 0.2e1 + t586) * t460 + t587 * mrSges(5,2) + t588 * mrSges(4,3) + t115 * t438 / 0.2e1 + (t514 * t600 + t515 * t603) * t288 - t377 * t455 / 0.2e1 - t40 * t388 - t208 * t390 + (t514 * t597 + t515 * t602) * t243 + t327 * t654 + Ifges(3,3) * qJDD(2) + Ifges(3,6) * t254 + Ifges(3,5) * t255 + (g(1) * t332 + g(2) * t329) * (-t349 * t406 + t391 + (pkin(9) * t671 + t657) * t331 + (-t359 * t409 - t671 * (-pkin(2) + t365) + m(4) * pkin(2) - m(5) * t258 - t662) * t328) + t258 * t68 - t238 * mrSges(3,2) - t239 * mrSges(3,1) + t232 * t15 + t380 * t551 + t621 * t100 + t622 * t73 + t623 * t160 - pkin(2) * t69 + t265 * t307 + t649 * t177 + t652 * (t330 * t422 + t415) - t45 * t97 - t52 * t98 + t103 * t14; (-t532 * t667 + t660) * t242 + (-m(5) * t411 - t671 * (-t214 * pkin(4) + t411) + t632 * t215 + t634 * t214 - t647) * g(2) + (-m(5) * t410 - t671 * (-t216 * pkin(4) + t410) + t632 * t217 + t634 * t216 - t650) * g(1) + (-t242 * t668 + t112 + t233 - t498) * t669 + (-t165 - t168) * t147 - t499 * t539 + t115 * t536 + t674 * t361 + (t1 * t249 + t2 * t250 - t30 * t36 + (-t44 + t179) * t23) * m(7) + (-t101 * t90 + t256 * t4 + t257 * t3 + (-t44 + t187) * t29) * m(6) + t617 * t44 + (t166 - t167) * t148 - t569 + (-m(7) * (t371 + t589) - t412 - m(6) * t371 - t413 - m(5) * t284 - (t491 + (-m(5) * pkin(3) - mrSges(5,1)) * t327) * t328 + t227) * g(3) + t594 + (-pkin(3) * t41 + qJ(4) * t39 - t105 * t148 + t107 * t583 - t111 * t155) * m(5) + t256 * t24 + t257 * t26 + t249 * t27 + t250 * t25 - t36 * t73 - t649 * t145 - pkin(3) * t86 + qJ(4) * t88 + (-Ifges(4,2) * t243 - t234 + t652) * t538 - t101 * t74 + (-m(6) * t28 + m(7) * t22 - t653) * (qJD(5) * t257 + t326 * t663 + t399 * t527) - t155 * t156 + qJD(4) * t168 + (Ifges(5,3) * t539 - t532 * t666 - t659) * t243 - t337 + t179 * t97 + t187 * t98; t86 - t624 * t527 - t617 * t429 + t623 * t326 + (t156 - t73 - t74) * t243 + (t1 * t326 - t2 * t527 - t30 * t243 + t279 * (t22 * t326 + t23 * t527) + t578) * m(7) + (-t90 * t243 + t3 * t326 + t4 * t527 + t279 * (-t28 * t326 + t29 * t527) + t578) * m(6) + (t111 * t243 + t41 + t578) * m(5) + t653 * (t288 * t326 - t457) + (-m(5) * t107 + t527 * t617 - t168) * t288; t59 * t545 + (Ifges(7,3) * t145 + t504) * t549 + (t512 + t653) * t29 + (t513 + t617) * t28 - pkin(5) * t25 + qJ(6) * t27 + (t194 * t406 + t195 * t409) * g(3) + t650 * g(1) + t647 * g(2) + (t145 * t23 - t22 * t361) * mrSges(7,2) - t72 * t73 + qJD(6) * t97 - t30 * (mrSges(7,1) * t145 - mrSges(7,3) * t361) - t90 * (mrSges(6,1) * t145 + mrSges(6,2) * t361) + t337 + (-t145 * t665 + t361 * t628) * t534 + (-pkin(5) * t2 + qJ(6) * t1 - t22 * t29 + t23 * t604 - t30 * t72) * m(7) + (-Ifges(6,2) * t145 + t140 + t618) * t548 + (t361 * t631 + t139 - t507 + t56) * t546; t145 * t73 - t279 * t97 + (-g(1) * t363 - g(2) * t607 - g(3) * t195 + t30 * t145 - t23 * t279 + t2) * m(7) + t25;];
tau  = t6;
