% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP11_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP11_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:41:35
% EndTime: 2019-03-09 17:42:54
% DurationCPUTime: 49.16s
% Computational Cost: add. (13818->969), mult. (33110->1228), div. (0->0), fcn. (25529->10), ass. (0->436)
t668 = Ifges(6,4) + Ifges(7,4);
t669 = Ifges(6,1) + Ifges(7,1);
t711 = Ifges(4,4) + Ifges(5,6);
t666 = Ifges(6,5) + Ifges(7,5);
t665 = Ifges(6,2) + Ifges(7,2);
t663 = Ifges(6,6) + Ifges(7,6);
t346 = cos(qJ(2));
t514 = cos(pkin(6));
t438 = t514 * qJD(1);
t424 = pkin(1) * t438;
t338 = sin(pkin(6));
t342 = sin(qJ(2));
t488 = qJD(1) * t342;
t464 = t338 * t488;
t257 = -pkin(8) * t464 + t346 * t424;
t371 = (pkin(2) * t342 - pkin(9) * t346) * t338;
t258 = qJD(1) * t371;
t341 = sin(qJ(3));
t345 = cos(qJ(3));
t167 = -t341 * t257 + t258 * t345;
t483 = qJD(3) * t345;
t489 = qJD(1) * t338;
t495 = t345 * t346;
t587 = pkin(4) + pkin(9);
t588 = pkin(3) + pkin(10);
t710 = -(pkin(4) * t495 - t342 * t588) * t489 + t167 + t587 * t483;
t340 = sin(qJ(5));
t344 = cos(qJ(5));
t496 = t344 * t346;
t229 = (-t340 * t342 + t341 * t496) * t338;
t217 = qJD(1) * t229;
t480 = qJD(5) * t345;
t459 = t340 * t480;
t709 = t217 - t459;
t477 = qJD(1) * qJD(2);
t367 = qJDD(1) * t342 + t346 * t477;
t356 = t338 * t367;
t379 = t438 + qJD(2);
t363 = qJD(3) * t379;
t435 = t514 * qJDD(1);
t374 = t435 + qJDD(2);
t505 = t338 * t342;
t471 = t341 * t505;
t425 = qJD(3) * t471;
t141 = qJD(1) * t425 - t341 * t374 + (-t356 - t363) * t345;
t581 = -t141 / 0.2e1;
t462 = t345 * t488;
t142 = t341 * t363 - t345 * t374 + (qJD(3) * t462 + t341 * t367) * t338;
t579 = -t142 / 0.2e1;
t264 = (-qJDD(1) * t346 + t342 * t477) * t338;
t251 = qJDD(3) + t264;
t562 = t251 / 0.2e1;
t667 = Ifges(4,5) - Ifges(5,4);
t664 = Ifges(4,6) - Ifges(5,5);
t662 = Ifges(4,3) + Ifges(5,1);
t661 = Ifges(7,3) + Ifges(6,3);
t227 = t341 * t464 - t345 * t379;
t463 = t346 * t489;
t304 = -qJD(3) + t463;
t177 = t227 * t344 + t304 * t340;
t708 = t668 * t177;
t484 = qJD(3) * t341;
t644 = -t344 * t484 + t709;
t426 = t345 * t463;
t707 = t426 - t483;
t178 = t227 * t340 - t304 * t344;
t706 = t668 * t178;
t138 = qJDD(5) - t141;
t696 = -pkin(8) * t338 * t477 + pkin(1) * t435;
t476 = qJDD(1) * t338;
t697 = pkin(8) * t476 + qJD(2) * t424;
t179 = t342 * t696 + t346 * t697;
t159 = pkin(9) * t374 + t179;
t166 = t264 * pkin(2) + (-qJDD(1) * pkin(1) - pkin(9) * t367) * t338;
t260 = pkin(8) * t463 + t342 * t424;
t212 = pkin(9) * t379 + t260;
t219 = (-pkin(2) * t346 - pkin(9) * t342 - pkin(1)) * t489;
t43 = -t341 * t159 + t166 * t345 - t212 * t483 - t219 * t484;
t364 = qJDD(4) - t43;
t19 = -pkin(4) * t141 - t251 * t588 + t364;
t180 = -t342 * t697 + t346 * t696;
t160 = -pkin(2) * t374 - t180;
t228 = t338 * t462 + t341 * t379;
t348 = t141 * qJ(4) - t228 * qJD(4) + t160;
t26 = t142 * t588 + t348;
t139 = t212 * t341 - t345 * t219;
t372 = pkin(4) * t228 + t139;
t691 = qJD(4) + t372;
t82 = t304 * t588 + t691;
t211 = -pkin(2) * t379 - t257;
t349 = -t228 * qJ(4) + t211;
t87 = t227 * t588 + t349;
t29 = t340 * t82 + t344 * t87;
t4 = -qJD(5) * t29 + t344 * t19 - t26 * t340;
t67 = qJD(5) * t177 + t142 * t340 + t251 * t344;
t1 = pkin(5) * t138 - qJ(6) * t67 - qJD(6) * t178 + t4;
t580 = t141 / 0.2e1;
t582 = t138 / 0.2e1;
t68 = -qJD(5) * t178 + t142 * t344 - t251 * t340;
t589 = t68 / 0.2e1;
t590 = t67 / 0.2e1;
t481 = qJD(5) * t344;
t482 = qJD(5) * t340;
t3 = t340 * t19 + t344 * t26 + t82 * t481 - t482 * t87;
t2 = qJ(6) * t68 + qJD(6) * t177 + t3;
t603 = t4 * mrSges(6,1) - t3 * mrSges(6,2) - t2 * mrSges(7,2);
t660 = t138 * t661 + t663 * t68 + t666 * t67;
t702 = -t251 / 0.2e1;
t705 = t603 + 0.2e1 * Ifges(4,1) * t581 + Ifges(5,4) * t702 + t582 * t661 + t589 * t663 + t590 * t666 + (-t580 + t581) * Ifges(5,2) + t1 * mrSges(7,1) + t660 / 0.2e1 + t711 * t579 + (t667 + Ifges(4,5)) * t562;
t114 = t227 * pkin(3) + t349;
t140 = t345 * t212 + t341 * t219;
t291 = t304 * qJ(4);
t120 = t291 - t140;
t623 = t120 * mrSges(5,1) - t140 * mrSges(4,3);
t704 = t211 * mrSges(4,1) - t114 * mrSges(5,2) + t623;
t343 = sin(qJ(1));
t555 = cos(qJ(1));
t412 = t514 * t555;
t275 = t342 * t412 + t343 * t346;
t465 = t338 * t555;
t198 = t275 * t341 + t345 * t465;
t440 = t342 * t514;
t277 = -t343 * t440 + t346 * t555;
t503 = t338 * t345;
t202 = t277 * t341 - t343 * t503;
t272 = -t345 * t514 + t471;
t602 = g(1) * t202 + g(2) * t198 + g(3) * t272;
t703 = t1 * t344 + t2 * t340 - t602;
t578 = t142 / 0.2e1;
t658 = t138 * t666 + t668 * t68 + t669 * t67;
t701 = -t658 / 0.2e1;
t659 = t138 * t663 + t665 * t68 + t668 * t67;
t700 = -t659 / 0.2e1;
t222 = qJD(5) + t228;
t655 = t177 * t665 + t222 * t663 + t706;
t654 = t178 * t669 + t666 * t222 + t708;
t383 = pkin(10) * t341 - qJ(4) * t345;
t427 = t341 * t463;
t466 = pkin(3) * t427 + t260;
t148 = t383 * t463 + t466;
t699 = t148 * t340 + t344 * t710;
t434 = pkin(3) * t484 - qJD(4) * t341;
t226 = qJD(3) * t383 + t434;
t512 = qJ(4) * t341;
t443 = -pkin(2) - t512;
t283 = -t345 * t588 + t443;
t307 = t587 * t341;
t646 = -t283 * t482 + t307 * t481 + (-t148 + t226) * t344 + t710 * t340;
t339 = -qJ(6) - pkin(10);
t402 = mrSges(5,2) * t345 - mrSges(5,3) * t341;
t410 = mrSges(4,1) * t345 - mrSges(4,2) * t341;
t501 = t340 * t341;
t698 = -m(7) * (pkin(5) * t501 - t339 * t345) - t345 * mrSges(7,3) - t410 + t402;
t620 = -qJD(4) - t139;
t118 = pkin(3) * t304 - t620;
t695 = t118 * mrSges(5,1) + t139 * mrSges(4,3);
t694 = t668 * t344;
t693 = t668 * t340;
t684 = m(6) * pkin(10);
t605 = -m(7) * t339 + mrSges(4,1) - mrSges(5,2) + t684;
t683 = -m(7) - m(5);
t619 = m(6) - t683;
t692 = pkin(3) * t619 + t605;
t274 = t342 * t343 - t346 * t412;
t690 = -t198 * t340 - t274 * t344;
t689 = t198 * t344 - t274 * t340;
t688 = Ifges(4,6) * t702 + 0.2e1 * Ifges(5,3) * t578 + t711 * t580 + (t578 - t579) * Ifges(4,2) + (-t664 + Ifges(5,5)) * t562;
t28 = -t340 * t87 + t344 * t82;
t23 = -qJ(6) * t178 + t28;
t16 = pkin(5) * t222 + t23;
t24 = qJ(6) * t177 + t29;
t687 = t28 * mrSges(6,1) + t16 * mrSges(7,1) - t29 * mrSges(6,2) - t24 * mrSges(7,2);
t220 = Ifges(5,6) * t227;
t130 = -t304 * Ifges(5,4) - t228 * Ifges(5,2) + t220;
t221 = Ifges(4,4) * t227;
t615 = t228 * Ifges(4,1) - t304 * Ifges(4,5) + t177 * t663 + t178 * t666 + t222 * t661 - t221;
t686 = -t130 / 0.2e1 + t615 / 0.2e1 + t695;
t682 = -t264 / 0.2e1;
t681 = t356 / 0.2e1;
t680 = t374 / 0.2e1;
t332 = pkin(5) * t340 + qJ(4);
t679 = m(7) * t332;
t333 = pkin(5) * t344 + pkin(4);
t678 = m(7) * t333;
t673 = -mrSges(5,1) - mrSges(4,3);
t672 = -mrSges(7,1) - mrSges(6,1);
t671 = mrSges(4,2) - mrSges(5,3);
t670 = mrSges(6,2) + mrSges(7,2);
t498 = t341 * t346;
t230 = (t340 * t498 + t342 * t344) * t338;
t218 = qJD(1) * t230;
t436 = qJ(6) * t345 - t283;
t657 = qJ(6) * t218 + t436 * t481 + (-qJ(6) * t484 - qJD(5) * t307 + qJD(6) * t345 - t226) * t340 + t699 - t707 * pkin(5);
t478 = t344 * qJD(6);
t656 = -qJ(6) * t644 - t345 * t478 + t646;
t653 = m(7) * pkin(5) + mrSges(7,1);
t184 = mrSges(5,1) * t227 + mrSges(5,3) * t304;
t98 = -mrSges(6,1) * t177 + mrSges(6,2) * t178;
t652 = -t184 + t98;
t168 = t345 * t257 + t341 * t258;
t147 = -qJ(4) * t464 - t168;
t125 = -pkin(4) * t427 - t147;
t542 = pkin(9) + t333;
t651 = pkin(5) * t709 - t484 * t542 - t125;
t513 = qJ(4) * t227;
t119 = t228 * t588 + t513;
t494 = qJ(6) + t588;
t100 = -pkin(4) * t227 + t140;
t93 = t344 * t100;
t650 = t482 * t494 - t478 + pkin(5) * t227 - t93 - (-qJ(6) * t228 - t119) * t340;
t294 = t494 * t344;
t45 = t340 * t100 + t344 * t119;
t510 = t228 * t344;
t649 = -qJ(6) * t510 - qJD(5) * t294 - t340 * qJD(6) - t45;
t648 = pkin(5) * t481 + t228 * t333 - t620;
t188 = t344 * t283 + t340 * t307;
t647 = -qJD(5) * t188 - t226 * t340 + t699;
t645 = -t227 * t664 + t228 * t667 - t304 * t662;
t643 = -t340 * t484 + t344 * t480 + t218;
t507 = t274 * t345;
t642 = -pkin(3) * t507 - t274 * t512;
t439 = t346 * t514;
t276 = t342 * t555 + t343 * t439;
t506 = t276 * t345;
t641 = -pkin(3) * t506 - t276 * t512;
t640 = qJ(4) * t707 + t434 - t466;
t639 = -t587 * t484 - t125;
t281 = pkin(1) * t439 - pkin(8) * t505;
t404 = t340 * mrSges(7,1) + t344 * mrSges(7,2);
t407 = mrSges(6,1) * t340 + mrSges(6,2) * t344;
t638 = -m(6) * qJ(4) - t404 - t407 - t679;
t637 = -t114 * (-mrSges(5,2) * t341 - mrSges(5,3) * t345) - t211 * (mrSges(4,1) * t341 + mrSges(4,2) * t345);
t636 = t340 * t666 + t344 * t663;
t635 = -t340 * t663 + t344 * t666;
t634 = t344 * t665 + t693;
t633 = -t340 * t665 + t694;
t632 = -t341 * t664 + t345 * t667;
t631 = t340 * t669 + t694;
t630 = t344 * t669 - t693;
t628 = -t141 * t667 - t142 * t664 + t251 * t662;
t627 = -t481 - t510;
t511 = t228 * t340;
t626 = t482 + t511;
t42 = t345 * t159 + t341 * t166 - t212 * t484 + t219 * t483;
t625 = -t341 * t43 + t345 * t42;
t32 = -t251 * qJ(4) + t304 * qJD(4) - t42;
t40 = -pkin(3) * t251 + t364;
t624 = -t32 * t345 + t341 * t40;
t622 = -t3 * t340 - t344 * t4;
t537 = Ifges(3,4) * t342;
t594 = t338 ^ 2;
t618 = (t342 * (t346 * Ifges(3,1) - t537) / 0.2e1 - pkin(1) * (mrSges(3,1) * t342 + mrSges(3,2) * t346)) * t594;
t616 = -mrSges(6,1) - t653;
t611 = mrSges(3,1) - t698;
t610 = t671 - t679;
t432 = mrSges(3,3) * t464;
t609 = -m(4) * t211 + mrSges(3,1) * t379 - mrSges(4,1) * t227 - mrSges(4,2) * t228 - t432;
t405 = mrSges(7,1) * t344 - mrSges(7,2) * t340;
t408 = mrSges(6,1) * t344 - mrSges(6,2) * t340;
t88 = t100 - t291;
t57 = -pkin(5) * t177 + qJD(6) + t88;
t607 = t57 * t405 + t88 * t408;
t606 = -m(6) * t587 + mrSges(3,2) + t673;
t604 = -m(5) * qJ(4) + t638 + t671;
t601 = t606 - t678;
t600 = m(7) * t542 - t606;
t598 = -mrSges(6,3) - mrSges(7,3) - t605;
t597 = -t43 * mrSges(4,1) + t42 * mrSges(4,2) - t40 * mrSges(5,2) + t32 * mrSges(5,3);
t596 = t211 * mrSges(4,2) - t114 * mrSges(5,3) + t687;
t521 = t228 * Ifges(4,4);
t131 = -t227 * Ifges(4,2) - t304 * Ifges(4,6) + t521;
t583 = -t131 / 0.2e1;
t577 = -t177 / 0.2e1;
t576 = t177 / 0.2e1;
t575 = -t178 / 0.2e1;
t574 = t178 / 0.2e1;
t568 = -t222 / 0.2e1;
t567 = t222 / 0.2e1;
t566 = -t227 / 0.2e1;
t565 = t227 / 0.2e1;
t564 = -t228 / 0.2e1;
t563 = t228 / 0.2e1;
t560 = -t304 / 0.2e1;
t559 = t304 / 0.2e1;
t554 = pkin(1) * t338;
t502 = t338 * t346;
t196 = t272 * t344 + t340 * t502;
t552 = pkin(5) * t196;
t551 = pkin(9) * t274;
t550 = pkin(9) * t276;
t546 = t272 * pkin(10);
t335 = t345 * pkin(9);
t541 = mrSges(6,3) * t177;
t540 = mrSges(6,3) * t178;
t539 = mrSges(7,3) * t177;
t538 = mrSges(7,3) * t178;
t536 = Ifges(3,4) * t346;
t535 = Ifges(4,4) * t341;
t534 = Ifges(4,4) * t345;
t529 = Ifges(5,6) * t341;
t528 = Ifges(5,6) * t345;
t520 = t228 * Ifges(5,6);
t282 = pkin(1) * t440 + pkin(8) * t502;
t247 = pkin(9) * t514 + t282;
t491 = pkin(2) * t502 + pkin(9) * t505;
t248 = -t491 - t554;
t164 = -t341 * t247 + t248 * t345;
t326 = pkin(3) * t502;
t146 = -t164 + t326;
t273 = t341 * t514 + t342 * t503;
t105 = pkin(4) * t273 + pkin(10) * t502 + t146;
t246 = -pkin(2) * t514 - t281;
t265 = t272 * pkin(3);
t437 = t273 * qJ(4) - t265;
t144 = t246 - t437;
t115 = t144 + t546;
t47 = t340 * t105 + t344 * t115;
t504 = t338 * t343;
t500 = t340 * t345;
t499 = t341 * t344;
t497 = t344 * t345;
t165 = t345 * t247 + t341 * t248;
t308 = t345 * pkin(4) + t335;
t490 = t555 * pkin(1) + pkin(8) * t504;
t487 = qJD(2) * t346;
t479 = qJD(5) * t588;
t475 = pkin(9) * t484;
t474 = pkin(9) * t483;
t472 = qJ(4) * t502;
t468 = Ifges(3,5) * t356 - Ifges(3,6) * t264 + Ifges(3,3) * t374;
t467 = t277 * pkin(2) + t490;
t461 = qJD(2) * t505;
t460 = t338 * t487;
t456 = t505 / 0.2e1;
t20 = -t68 * mrSges(7,1) + t67 * mrSges(7,2);
t452 = -t489 / 0.2e1;
t451 = t489 / 0.2e1;
t448 = t484 / 0.2e1;
t445 = -t481 / 0.2e1;
t444 = -pkin(1) * t343 + pkin(8) * t465;
t266 = t274 * pkin(2);
t442 = pkin(9) * t275 - t266;
t268 = t276 * pkin(2);
t441 = pkin(9) * t277 - t268;
t104 = -t141 * mrSges(5,1) + t251 * mrSges(5,2);
t46 = t344 * t105 - t115 * t340;
t199 = t275 * t345 - t341 * t465;
t431 = mrSges(3,3) * t463;
t203 = t277 * t345 + t341 * t504;
t430 = t203 * pkin(3) + t467;
t421 = t346 * t452;
t420 = t346 * t451;
t414 = -t275 * pkin(2) + t444;
t411 = mrSges(4,1) * t272 + mrSges(4,2) * t273;
t197 = -t272 * t340 + t338 * t496;
t409 = mrSges(6,1) * t196 + mrSges(6,2) * t197;
t406 = -t196 * mrSges(7,1) - t197 * mrSges(7,2);
t403 = -t272 * mrSges(5,2) - t273 * mrSges(5,3);
t401 = Ifges(4,1) * t345 - t535;
t396 = -Ifges(4,2) * t341 + t534;
t385 = -Ifges(5,2) * t345 + t529;
t384 = Ifges(5,3) * t341 - t528;
t381 = t28 * t340 - t29 * t344;
t154 = t202 * t344 - t276 * t340;
t378 = -pkin(3) * t199 + t414;
t145 = t472 - t165;
t259 = qJD(2) * t371;
t261 = t281 * qJD(2);
t90 = -t247 * t483 - t248 * t484 + t259 * t345 - t341 * t261;
t195 = -t425 + (qJD(3) * t514 + t460) * t345;
t51 = pkin(4) * t195 - t461 * t588 - t90;
t194 = qJD(3) * t273 + t341 * t460;
t262 = t282 * qJD(2);
t350 = -t195 * qJ(4) - t273 * qJD(4) + t262;
t61 = t194 * t588 + t350;
t7 = t105 * t481 - t115 * t482 + t340 * t51 + t344 * t61;
t368 = qJ(4) * t202 + t430;
t89 = -t247 * t484 + t248 * t483 + t341 * t259 + t345 * t261;
t116 = -pkin(4) * t272 - t145;
t357 = -qJ(4) * t198 + t378;
t8 = -qJD(5) * t47 - t340 * t61 + t344 * t51;
t353 = Ifges(3,6) * t514 + (Ifges(3,2) * t346 + t537) * t338;
t22 = -pkin(4) * t142 - t32;
t72 = -qJ(4) * t461 + qJD(4) * t502 - t89;
t352 = -qJD(5) * t381 - t622;
t351 = t338 * t379 * (Ifges(3,5) * t346 - Ifges(3,6) * t342);
t52 = -pkin(4) * t194 - t72;
t318 = Ifges(3,4) * t463;
t298 = -pkin(3) * t345 + t443;
t293 = t494 * t340;
t286 = t344 * t307;
t278 = (-mrSges(3,1) * t346 + mrSges(3,2) * t342) * t338;
t271 = pkin(5) * t497 + t308;
t256 = -mrSges(3,2) * t379 + t431;
t208 = Ifges(3,1) * t464 + Ifges(3,5) * t379 + t318;
t207 = Ifges(3,6) * qJD(2) + qJD(1) * t353;
t187 = -t283 * t340 + t286;
t185 = mrSges(5,1) * t228 - mrSges(5,2) * t304;
t183 = -mrSges(4,1) * t304 - mrSges(4,3) * t228;
t182 = mrSges(4,2) * t304 - mrSges(4,3) * t227;
t181 = -qJ(6) * t497 + t188;
t169 = pkin(5) * t341 + t340 * t436 + t286;
t163 = -mrSges(5,2) * t227 - mrSges(5,3) * t228;
t161 = pkin(3) * t228 + t513;
t155 = t202 * t340 + t276 * t344;
t149 = -pkin(3) * t464 - t167;
t128 = -t304 * Ifges(5,5) + t227 * Ifges(5,3) - t520;
t124 = mrSges(6,1) * t222 - t540;
t123 = mrSges(7,1) * t222 - t538;
t122 = -mrSges(6,2) * t222 + t541;
t121 = -mrSges(7,2) * t222 + t539;
t112 = qJD(5) * t196 + t194 * t340 + t344 * t461;
t111 = qJD(5) * t197 + t194 * t344 - t340 * t461;
t103 = mrSges(5,1) * t142 - mrSges(5,3) * t251;
t102 = -mrSges(4,2) * t251 - mrSges(4,3) * t142;
t101 = mrSges(4,1) * t251 + mrSges(4,3) * t141;
t97 = -mrSges(7,1) * t177 + mrSges(7,2) * t178;
t86 = t194 * pkin(3) + t350;
t84 = -pkin(3) * t461 - t90;
t83 = t116 - t552;
t70 = mrSges(4,1) * t142 - mrSges(4,2) * t141;
t69 = -mrSges(5,2) * t142 + mrSges(5,3) * t141;
t44 = -t119 * t340 + t93;
t38 = -mrSges(6,2) * t138 + mrSges(6,3) * t68;
t37 = -mrSges(7,2) * t138 + mrSges(7,3) * t68;
t36 = mrSges(6,1) * t138 - mrSges(6,3) * t67;
t35 = mrSges(7,1) * t138 - mrSges(7,3) * t67;
t34 = t142 * pkin(3) + t348;
t33 = qJ(6) * t196 + t47;
t31 = pkin(5) * t273 + qJ(6) * t197 + t46;
t27 = -pkin(5) * t111 + t52;
t21 = -mrSges(6,1) * t68 + mrSges(6,2) * t67;
t9 = -pkin(5) * t68 + qJDD(6) + t22;
t6 = qJ(6) * t111 + qJD(6) * t196 + t7;
t5 = pkin(5) * t195 - qJ(6) * t112 + qJD(6) * t197 + t8;
t10 = [t514 * t468 / 0.2e1 + (Ifges(4,1) * t563 + Ifges(4,4) * t566 - Ifges(5,2) * t564 - Ifges(5,6) * t565 + t667 * t560 + t661 * t567 + t666 * t574 + t663 * t576 + t596 + t686) * t195 + (t594 * qJD(1) * (-Ifges(3,2) * t342 + t536) + t338 * t208) * t487 / 0.2e1 + m(3) * (pkin(1) ^ 2 * qJDD(1) * t594 + t179 * t282 + t180 * t281 + t260 * t261) - t22 * t409 + t160 * t411 + (t111 * t29 - t112 * t28 + t196 * t3 + t197 * t4) * mrSges(6,3) + (t128 / 0.2e1 - Ifges(4,4) * t563 + Ifges(5,6) * t564 + Ifges(5,3) * t565 - Ifges(4,2) * t566 + t583 - t664 * t560 + t704) * t194 + (t40 * mrSges(5,1) - t43 * mrSges(4,3) + Ifges(4,4) * t579 - Ifges(5,6) * t578 + t705) * t273 + (-t180 * t505 - t257 * t460 - t260 * t461 - t264 * t282 - t281 * t356) * mrSges(3,3) + (t1 * t197 + t111 * t24 - t112 * t16 + t196 * t2) * mrSges(7,3) + (t645 * t456 + t351 / 0.2e1) * qJD(2) + (t32 * mrSges(5,1) - t42 * mrSges(4,3) - Ifges(4,4) * t581 + Ifges(5,6) * t580 + t688) * t272 + (-m(7) * t378 - m(6) * t357 - m(3) * t444 + t275 * mrSges(3,1) - mrSges(3,3) * t465 - m(4) * (t414 - t551) - m(5) * (t357 - t551) + t343 * mrSges(2,1) + mrSges(2,2) * t555 - t610 * t198 + t672 * t690 + t670 * t689 + t600 * t274 - t598 * t199) * g(1) + t34 * t403 + t9 * t406 + (-t179 * t514 - t282 * t374 - t356 * t554) * mrSges(3,2) + (Ifges(3,4) * t681 - Ifges(5,4) * t580 - Ifges(4,5) * t581 - Ifges(5,5) * t578 + Ifges(3,2) * t682 + Ifges(3,6) * t680 - Ifges(4,6) * t579 - t662 * t562 + t597 + t179 * mrSges(3,3) - t628 / 0.2e1) * t502 + t33 * t37 + t31 * t35 - pkin(1) * t278 * t476 + m(4) * (-t139 * t90 + t140 * t89 + t160 * t246 + t164 * t43 + t165 * t42) + m(5) * (t114 * t86 + t118 * t84 + t120 * t72 + t144 * t34 + t145 * t32 + t146 * t40) + m(6) * (t116 * t22 + t28 * t8 + t29 * t7 + t3 * t47 + t4 * t46 + t52 * t88) + m(7) * (t1 * t31 + t16 * t5 + t2 * t33 + t24 * t6 + t27 * t57 + t83 * t9) + Ifges(2,3) * qJDD(1) + t261 * t256 + t246 * t70 + (-m(7) * t430 - m(3) * t490 - t277 * mrSges(3,1) - mrSges(3,3) * t504 - m(6) * t368 - m(5) * (t368 + t550) - m(4) * (t467 + t550) - mrSges(2,1) * t555 + t343 * mrSges(2,2) + t610 * t202 + t672 * t155 - t670 * t154 - t600 * t276 + t598 * t203) * g(2) + (t111 * t663 + t112 * t666) * t567 + (t196 * t663 - t197 * t666) * t582 + (t111 * t665 + t112 * t668) * t576 + (t196 * t665 - t197 * t668) * t589 + (t111 * t668 + t112 * t669) * t574 + (t196 * t668 - t197 * t669) * t590 + t659 * t196 / 0.2e1 + (-t120 * mrSges(5,3) - t140 * mrSges(4,2) + t118 * mrSges(5,2) - t139 * mrSges(4,1) - t207 / 0.2e1 + Ifges(4,5) * t563 + Ifges(5,4) * t564 + Ifges(5,5) * t565 + Ifges(4,6) * t566 + t662 * t560) * t461 + t654 * t112 / 0.2e1 + t655 * t111 / 0.2e1 + t90 * t183 + t72 * t184 + t84 * t185 + t89 * t182 + (t180 * t514 - t264 * t554 + t281 * t374) * mrSges(3,1) + (Ifges(3,1) * t356 - Ifges(3,4) * t264 + Ifges(3,5) * t374) * t456 + t46 * t36 + t47 * t38 + (-m(3) * t257 - t609) * t262 + t83 * t20 + t618 * t477 + t27 * t97 + t52 * t98 + (Ifges(3,3) * t514 + (Ifges(3,5) * t342 + Ifges(3,6) * t346) * t338) * t680 + (Ifges(3,5) * t514 + (t342 * Ifges(3,1) + t536) * t338) * t681 + t353 * t682 + t57 * (-mrSges(7,1) * t111 + mrSges(7,2) * t112) + t88 * (-mrSges(6,1) * t111 + mrSges(6,2) * t112) + t116 * t21 + t6 * t121 + t7 * t122 + t5 * t123 + t8 * t124 + t144 * t69 + t145 * t103 + t146 * t104 + t86 * t163 + t164 * t101 + t165 * t102 + t197 * t701; -t529 * t578 + t535 * t579 + (-t475 - t168) * t182 + (t686 + t687) * t483 + (t130 * t420 + t22 * t408 + t9 * t405 - t582 * t636 - t589 * t634 - t590 * t631 - t688) * t345 + (-t474 - t167) * t183 + t654 * (t340 * t448 + t345 * t445 - t218 / 0.2e1) + t655 * (t459 / 0.2e1 + t344 * t448 - t217 / 0.2e1) + (-m(4) * t491 + t278 + t672 * t230 - t670 * t229 - t619 * (t345 * t326 + t491) + ((-mrSges(6,3) - t684) * t495 + t698 * t346 + (-m(6) * pkin(4) + t673 - t678) * t342) * t338) * g(3) - t528 * t580 + t534 * t581 + (-pkin(2) * t160 + pkin(9) * t625 + t139 * t167 - t140 * t168) * m(4) + (pkin(9) * t624 + t640 * t114 - t118 * t149 - t120 * t147 + t298 * t34) * m(5) - t160 * t410 + (t24 * t426 - t57 * t643) * mrSges(7,2) + (t131 * t420 + t128 * t421 + (t104 - t101) * pkin(9) - t619 * t472 * g(3) + t705) * t341 + (-t28 * t426 + t644 * t88) * mrSges(6,1) + (-t351 / 0.2e1 - t618 * qJD(1)) * qJD(1) + (t228 * (Ifges(4,5) * t342 + t346 * t401) + t227 * (Ifges(5,5) * t342 + t346 * t384) + t645 * t342) * t452 + (-t140 * (-mrSges(4,2) * t342 - mrSges(4,3) * t498) - t120 * (mrSges(5,1) * t498 - mrSges(5,3) * t342) + t139 * (mrSges(4,1) * t342 - mrSges(4,3) * t495) - t118 * (mrSges(5,1) * t495 + mrSges(5,2) * t342)) * t489 + (-Ifges(3,2) * t464 + t345 * t615 + t208 + t318) * t421 + (-t16 * t426 + t57 * t644) * mrSges(7,1) + (t342 * t207 + t228 * (Ifges(5,4) * t342 + t346 * t385) + t227 * (Ifges(4,6) * t342 + t346 * t396) + (t342 * t662 + t346 * t632) * t304) * t451 + (-t256 + t431) * t257 + (-t103 + t102) * t335 + (t474 - t149) * t185 + t34 * t402 + (t475 - t147) * t184 + (t29 * t426 - t643 * t88) * mrSges(6,2) + (t632 * t560 + (t401 / 0.2e1 - t385 / 0.2e1) * t228 + (-t396 / 0.2e1 + t384 / 0.2e1) * t227 + (t341 * t634 + t345 * t663) * t576 + (t341 * t631 + t345 * t666) * t574 + (t341 * t636 + t345 * t661) * t567 - t637 + ((t139 * t345 - t140 * t341) * m(4) + (t118 * t345 + t120 * t341) * m(5)) * pkin(9)) * qJD(3) + t468 + t308 * t21 + t298 * t69 + t271 * t20 + (-m(6) * (-pkin(10) * t506 - t268 + t641) + mrSges(6,3) * t506 - m(4) * t441 + t683 * (t441 + t641) + t672 * (-t276 * t501 + t277 * t344) - t670 * (-t276 * t499 - t277 * t340) + t601 * t277 + t611 * t276) * g(1) + (-m(6) * (-pkin(10) * t507 - t266 + t642) + mrSges(6,3) * t507 - m(4) * t442 + t683 * (t442 + t642) + t672 * (-t274 * t501 + t275 * t344) - t670 * (-t274 * t499 - t275 * t340) + t601 * t275 + t611 * t274) * g(2) + (t217 * t665 + t218 * t668 + t426 * t663) * t577 + (t217 * t668 + t218 * t669 + t426 * t666) * t575 + (t217 * t663 + t218 * t666 + t426 * t661) * t568 + t656 * t121 + t657 * t123 + (t1 * t169 + t16 * t657 + t181 * t2 + t24 * t656 + t271 * t9 + t57 * t651) * m(7) + t187 * t36 + t188 * t38 - t179 * mrSges(3,2) + t180 * mrSges(3,1) + t181 * t37 + (t432 + t609) * t260 - pkin(2) * t70 + (t583 + t623) * t484 + t624 * mrSges(5,1) + t128 * t448 + t625 * mrSges(4,3) + (-t567 * t635 - t574 * t630 - t576 * t633) * t480 + t637 * t463 + t639 * t98 + t640 * t163 + (t28 * t643 - t29 * t644 - t3 * t497 + t4 * t500) * mrSges(6,3) + (t1 * t500 + t16 * t643 - t2 * t497 - t24 * t644) * mrSges(7,3) + t169 * t35 + t646 * t122 + t647 * t124 + (t187 * t4 + t188 * t3 + t22 * t308 + t28 * t647 + t29 * t646 + t639 * t88) * m(6) + t651 * t97 + t497 * t700 + t500 * t701; (-t482 / 0.2e1 - t511 / 0.2e1) * t654 + (t183 - t185) * t140 + (-t184 + t182) * t139 + (t198 * t692 + t199 * t604) * g(2) + (t202 * t692 + t203 * t604) * g(1) + (-Ifges(4,1) * t564 + Ifges(5,2) * t563 - t559 * t667 - t568 * t661 - t575 * t666 - t577 * t663 + t596 + t695) * t227 + (t22 * qJ(4) - t28 * t44 - t29 * t45 - t588 * t352 + t691 * t88) * m(6) + (t16 * t626 + t24 * t627 - t703) * mrSges(7,3) + (-Ifges(4,2) * t565 + Ifges(5,3) * t566 - t559 * t664 + t568 * t636 + t575 * t631 + t577 * t634 + t607 - t704) * t228 + (-t103 + t21) * qJ(4) + (t520 + t131) * t563 + t372 * t98 - (t340 * t38 + t344 * t36) * t588 - t597 + t9 * t404 + t22 * t407 + (t445 - t510 / 0.2e1) * t655 + t607 * qJD(5) + (-t521 + t128) * t564 + (t220 + t130) * t566 + t628 + t332 * t20 - t293 * t37 - t294 * t35 + t658 * t344 / 0.2e1 + (t340 * t479 - t44) * t124 + (-t344 * t479 - t45) * t122 + (-t221 + t615) * t565 + (-pkin(3) * t40 - qJ(4) * t32 - t114 * t161 - t118 * t140 + t120 * t620) * m(5) - pkin(3) * t104 + (t28 * t626 + t29 * t627 + t602 + t622) * mrSges(6,3) + t630 * t590 + t633 * t589 + t635 * t582 - (t177 * t634 + t178 * t631 + t222 * t636) * qJD(5) / 0.2e1 + (-m(7) * (t272 * t339 - t265) + t411 - m(6) * (-t265 - t546) - m(5) * t437 + t403 + t638 * t273) * g(3) - t161 * t163 + t648 * t97 + t649 * t121 + t650 * t123 + (-t1 * t294 + t16 * t650 - t2 * t293 + t24 * t649 + t332 * t9 + t57 * t648) * m(7) + t652 * qJD(4) + t340 * t700; t228 * t163 + (t97 + t652) * t304 + (t35 + t36 + t222 * (t121 + t122)) * t344 + (t37 + t38 + t222 * (-t123 - t124)) * t340 + t104 + (t304 * t57 - t222 * (t16 * t340 - t24 * t344) + t703) * m(7) + (-t228 * t381 + t304 * t88 + t352 - t602) * m(6) + (t114 * t228 - t120 * t304 + t40 - t602) * m(5); (t177 * t669 - t706) * t575 + (t541 - t122) * t28 + t603 + (-t178 * t665 + t654 + t708) * t577 + (t154 * t616 + t155 * t670) * g(1) + t655 * t574 + (t540 + t124) * t29 + (t616 * t689 - t670 * t690) * g(2) + (t177 * t666 - t178 * t663) * t568 + (-m(7) * t552 + t406 - t409) * g(3) + t653 * t1 + (t538 + t123 - m(7) * (-t16 + t23)) * t24 - t57 * (mrSges(7,1) * t178 + mrSges(7,2) * t177) - t88 * (mrSges(6,1) * t178 + mrSges(6,2) * t177) + t16 * t539 - t23 * t121 + t660 + ((-m(7) * t57 - t97) * t178 + t35) * pkin(5); -t177 * t121 + t178 * t123 + (-g(1) * t203 - g(2) * t199 - g(3) * t273 + t16 * t178 - t24 * t177 + t9) * m(7) + t20;];
tau  = t10;
