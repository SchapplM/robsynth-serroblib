% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:29:35
% EndTime: 2019-03-09 17:31:10
% DurationCPUTime: 64.98s
% Computational Cost: add. (22486->1026), mult. (54449->1372), div. (0->0), fcn. (43836->14), ass. (0->425)
t419 = sin(qJ(2));
t547 = cos(pkin(6));
t499 = pkin(1) * t547;
t403 = t419 * t499;
t418 = sin(qJ(3));
t421 = cos(qJ(3));
t454 = pkin(3) * t418 - qJ(4) * t421;
t414 = sin(pkin(6));
t422 = cos(qJ(2));
t531 = t414 * t422;
t719 = (t403 + (pkin(8) + t454) * t531) * qJD(1) - qJD(3) * t454 + qJD(4) * t418;
t480 = t547 * qJD(1);
t471 = pkin(1) * t480;
t517 = qJD(1) * t414;
t495 = t419 * t517;
t336 = -pkin(8) * t495 + t422 * t471;
t441 = (pkin(2) * t419 - pkin(9) * t422) * t414;
t337 = qJD(1) * t441;
t229 = t421 * t336 + t418 * t337;
t200 = qJ(4) * t495 + t229;
t413 = sin(pkin(11));
t415 = cos(pkin(11));
t718 = t415 * t200 + t413 * t719;
t513 = qJD(3) * t418;
t506 = pkin(9) * t513;
t646 = -t719 * t415 + (t200 + t506) * t413;
t412 = pkin(11) + qJ(5);
t408 = sin(t412);
t409 = cos(t412);
t478 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t708 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t717 = t408 * t708 + t409 * t478;
t526 = t421 * t422;
t292 = (t413 * t419 + t415 * t526) * t517;
t494 = t422 * t517;
t472 = t418 * t494;
t529 = t415 * t421;
t715 = -pkin(4) * t472 + t292 * pkin(10) + (pkin(4) * t418 - pkin(10) * t529) * qJD(3) + t646;
t291 = (-t413 * t526 + t415 * t419) * t517;
t530 = t415 * t418;
t533 = t413 * t421;
t714 = pkin(10) * t291 - (-pkin(9) * t530 - pkin(10) * t533) * qJD(3) + t718;
t444 = t480 + qJD(2);
t305 = t418 * t495 - t421 * t444;
t306 = t418 * t444 + t421 * t495;
t551 = t306 * Ifges(4,4);
t380 = qJD(3) - t494;
t661 = t380 * Ifges(4,6);
t173 = -t305 * Ifges(4,2) + t551 + t661;
t298 = qJD(5) + t305;
t582 = t298 / 0.2e1;
t237 = t306 * t415 + t380 * t413;
t417 = sin(qJ(5));
t420 = cos(qJ(5));
t447 = t306 * t413 - t380 * t415;
t693 = t420 * t237 - t417 * t447;
t597 = t693 / 0.2e1;
t142 = t417 * t237 + t420 * t447;
t599 = t142 / 0.2e1;
t600 = -t142 / 0.2e1;
t665 = Ifges(5,6) * t447;
t673 = Ifges(6,3) + Ifges(7,2);
t674 = -Ifges(6,6) + Ifges(7,6);
t675 = Ifges(7,4) + Ifges(6,5);
t703 = Ifges(5,5) * t237;
t635 = Ifges(5,3) * t305 + t142 * t674 + t298 * t673 + t675 * t693 - t665 + t703;
t713 = Ifges(6,6) * t600 + Ifges(7,6) * t599 + t673 * t582 + t675 * t597 - t173 / 0.2e1 + t635 / 0.2e1;
t677 = Ifges(6,1) + Ifges(7,1);
t676 = -Ifges(6,4) + Ifges(7,5);
t183 = -t420 * t291 + t292 * t417;
t366 = t413 * t420 + t415 * t417;
t509 = qJD(5) * t420;
t510 = qJD(5) * t418;
t512 = qJD(3) * t421;
t535 = t413 * t417;
t233 = t366 * t512 + t509 * t530 - t510 * t535;
t712 = -t183 + t233;
t184 = t291 * t417 + t292 * t420;
t365 = -t420 * t415 + t535;
t232 = -t365 * t512 - t366 * t510;
t650 = t184 - t232;
t185 = t366 * t305;
t350 = t366 * qJD(5);
t649 = t185 + t350;
t186 = t365 * t305;
t349 = t365 * qJD(5);
t648 = t186 + t349;
t228 = -t418 * t336 + t337 * t421;
t202 = -pkin(3) * t495 - t228;
t407 = pkin(9) * t512;
t711 = t407 - t202;
t643 = t472 - t513;
t623 = Ifges(6,4) * t600 + Ifges(7,5) * t599 + t675 * t582 + t677 * t597;
t140 = Ifges(6,4) * t142;
t556 = Ifges(7,5) * t142;
t670 = t675 * t298 + t677 * t693 - t140 + t556;
t710 = t623 + t670 / 0.2e1;
t284 = -pkin(2) * t444 - t336;
t159 = t305 * pkin(3) - t306 * qJ(4) + t284;
t519 = pkin(8) * t531 + t403;
t326 = t547 * pkin(9) + t519;
t285 = qJD(2) * pkin(9) + qJD(1) * t326;
t296 = (-pkin(2) * t422 - pkin(9) * t419 - pkin(1)) * t517;
t181 = t421 * t285 + t418 * t296;
t164 = qJ(4) * t380 + t181;
t101 = t415 * t159 - t164 * t413;
t102 = t413 * t159 + t415 * t164;
t65 = pkin(4) * t305 - pkin(10) * t237 + t101;
t78 = -pkin(10) * t447 + t102;
t25 = -t417 * t78 + t420 * t65;
t651 = qJD(6) - t25;
t23 = -pkin(5) * t298 + t651;
t26 = t417 * t65 + t420 * t78;
t24 = qJ(6) * t298 + t26;
t554 = t181 * mrSges(4,3);
t709 = -t554 - t661 / 0.2e1 + t284 * mrSges(4,1) + t101 * mrSges(5,1) + t25 * mrSges(6,1) - t23 * mrSges(7,1) - t102 * mrSges(5,2) - t26 * mrSges(6,2) + t24 * mrSges(7,3);
t706 = -t670 / 0.2e1;
t180 = -t418 * t285 + t421 * t296;
t161 = -t380 * pkin(3) + qJD(4) - t180;
t126 = pkin(4) * t447 + t161;
t43 = t142 * pkin(5) - qJ(6) * t693 + t126;
t681 = t43 * mrSges(7,3);
t686 = -m(6) - m(7);
t704 = -m(4) + t686;
t372 = -pkin(3) * t421 - qJ(4) * t418 - pkin(2);
t363 = t415 * t372;
t253 = -pkin(10) * t530 + t363 + (-pkin(9) * t413 - pkin(4)) * t421;
t304 = pkin(9) * t529 + t413 * t372;
t534 = t413 * t418;
t269 = -pkin(10) * t534 + t304;
t511 = qJD(5) * t417;
t667 = t253 * t509 - t269 * t511 + t417 * t715 - t420 * t714;
t647 = t417 * t253 + t420 * t269;
t666 = -qJD(5) * t647 + t417 * t714 + t420 * t715;
t212 = pkin(3) * t306 + qJ(4) * t305;
t124 = t415 * t180 + t413 * t212;
t546 = t305 * t413;
t109 = pkin(10) * t546 + t124;
t566 = pkin(10) + qJ(4);
t374 = t566 * t413;
t375 = t566 * t415;
t446 = -t420 * t374 - t375 * t417;
t123 = -t180 * t413 + t415 * t212;
t545 = t305 * t415;
t95 = pkin(4) * t306 + pkin(10) * t545 + t123;
t654 = -qJD(4) * t365 + qJD(5) * t446 - t420 * t109 - t417 * t95;
t280 = -t374 * t417 + t375 * t420;
t653 = -qJD(4) * t366 - qJD(5) * t280 + t109 * t417 - t420 * t95;
t568 = pkin(4) * t413;
t644 = pkin(4) * t291 + t512 * t568 + t711;
t463 = -mrSges(5,1) * t415 + mrSges(5,2) * t413;
t438 = m(5) * pkin(3) - t463;
t464 = t421 * mrSges(4,1) - t418 * mrSges(4,2);
t488 = m(5) * qJ(4) + mrSges(5,3);
t702 = -t418 * t488 - t421 * t438 - t464;
t507 = qJDD(1) * t414;
t701 = pkin(8) * t507 + qJD(2) * t471;
t479 = t547 * qJDD(1);
t508 = qJD(1) * qJD(2);
t700 = -pkin(8) * t414 * t508 + pkin(1) * t479;
t699 = mrSges(7,2) * t24 + mrSges(6,3) * t26;
t698 = t681 + t706;
t139 = Ifges(7,5) * t693;
t68 = t298 * Ifges(7,6) + t142 * Ifges(7,3) + t139;
t557 = Ifges(6,4) * t693;
t71 = -t142 * Ifges(6,2) + t298 * Ifges(6,6) + t557;
t695 = t43 * mrSges(7,1) + t68 / 0.2e1 - t71 / 0.2e1;
t621 = -Ifges(6,2) * t600 + Ifges(7,3) * t599 + t582 * t674 + t597 * t676 + t695;
t697 = t621 - t699;
t571 = cos(qJ(1));
t467 = t547 * t571;
t570 = sin(qJ(1));
t354 = t419 * t467 + t422 * t570;
t498 = t414 * t571;
t271 = t354 * t421 - t418 * t498;
t353 = t419 * t570 - t422 * t467;
t189 = t271 * t408 - t353 * t409;
t696 = t271 * t409 + t353 * t408;
t343 = (qJDD(1) * t419 + t422 * t508) * t414;
t393 = t479 + qJDD(2);
t188 = qJD(3) * t306 + t418 * t343 - t421 * t393;
t587 = t188 / 0.2e1;
t514 = qJD(3) * t305;
t187 = t421 * t343 + t418 * t393 - t514;
t342 = (-qJDD(1) * t422 + t419 * t508) * t414;
t331 = qJDD(3) + t342;
t150 = t187 * t415 + t331 * t413;
t595 = t150 / 0.2e1;
t694 = Ifges(5,1) * t595 + Ifges(5,5) * t587;
t179 = qJDD(5) + t188;
t240 = t419 * t700 + t422 * t701;
t210 = pkin(9) * t393 + t240;
t505 = pkin(1) * t507;
t222 = pkin(2) * t342 - pkin(9) * t343 - t505;
t98 = t421 * t210 + t418 * t222 - t285 * t513 + t296 * t512;
t77 = qJ(4) * t331 + qJD(4) * t380 + t98;
t241 = -t419 * t701 + t422 * t700;
t211 = -t393 * pkin(2) - t241;
t81 = t188 * pkin(3) - t187 * qJ(4) - t306 * qJD(4) + t211;
t27 = -t413 * t77 + t415 * t81;
t18 = pkin(4) * t188 - pkin(10) * t150 + t27;
t149 = -t187 * t413 + t331 * t415;
t28 = t413 * t81 + t415 * t77;
t22 = pkin(10) * t149 + t28;
t4 = -qJD(5) * t26 + t18 * t420 - t22 * t417;
t2 = -pkin(5) * t179 + qJDD(6) - t4;
t99 = -t418 * t210 + t222 * t421 - t285 * t512 - t296 * t513;
t82 = -pkin(3) * t331 + qJDD(4) - t99;
t48 = -pkin(4) * t149 + t82;
t592 = t179 / 0.2e1;
t58 = qJD(5) * t693 - t420 * t149 + t417 * t150;
t611 = t58 / 0.2e1;
t612 = -t58 / 0.2e1;
t57 = -qJD(5) * t142 + t417 * t149 + t420 * t150;
t613 = t57 / 0.2e1;
t671 = t675 * t179 + t57 * t677 + t676 * t58;
t7 = pkin(5) * t58 - qJ(6) * t57 - qJD(6) * t693 + t48;
t692 = t48 * mrSges(6,2) - t7 * mrSges(7,3) + Ifges(6,4) * t612 + Ifges(7,5) * t611 + t675 * t592 + t677 * t613 - mrSges(6,3) * t4 + mrSges(7,2) * t2 + t671 / 0.2e1;
t660 = t415 * mrSges(5,2);
t462 = t413 * mrSges(5,1) + t660;
t679 = mrSges(3,2) - mrSges(4,3);
t691 = -m(5) * pkin(9) - t408 * t478 + t409 * t708 + t568 * t686 - t462 + t679;
t406 = pkin(4) * t415 + pkin(3);
t678 = mrSges(6,3) + mrSges(7,2);
t689 = mrSges(3,1) - t702 + (-t566 * t686 + t678) * t418 + (-t406 * t686 + t717) * t421;
t687 = -m(5) - m(4);
t596 = t149 / 0.2e1;
t589 = t187 / 0.2e1;
t588 = -t188 / 0.2e1;
t577 = t331 / 0.2e1;
t685 = -t447 / 0.2e1;
t680 = t98 * mrSges(4,2);
t672 = t179 * t673 + t57 * t675 + t58 * t674;
t669 = -qJ(6) * t643 - qJD(6) * t421 + t667;
t668 = pkin(5) * t643 - t666;
t662 = t380 * Ifges(4,5);
t659 = mrSges(4,2) - t488;
t532 = t414 * t419;
t358 = -pkin(8) * t532 + t422 * t499;
t325 = -pkin(2) * t547 - t358;
t351 = t418 * t532 - t421 * t547;
t352 = t418 * t547 + t421 * t532;
t196 = t351 * pkin(3) - t352 * qJ(4) + t325;
t520 = pkin(2) * t531 + pkin(9) * t532;
t569 = pkin(1) * t414;
t327 = -t520 - t569;
t215 = t421 * t326 + t418 * t327;
t197 = -qJ(4) * t531 + t215;
t117 = t413 * t196 + t415 * t197;
t265 = -t352 * t413 - t415 * t531;
t103 = pkin(10) * t265 + t117;
t116 = t415 * t196 - t197 * t413;
t266 = t352 * t415 - t413 * t531;
t91 = pkin(4) * t351 - pkin(10) * t266 + t116;
t658 = t420 * t103 + t417 * t91;
t333 = t365 * t418;
t657 = pkin(5) * t712 + qJ(6) * t650 + qJD(6) * t333 + t644;
t148 = -pkin(4) * t546 + t181;
t656 = pkin(5) * t649 + qJ(6) * t648 - qJD(6) * t366 - t148;
t655 = -qJ(6) * t306 + t654;
t652 = pkin(5) * t306 - t653;
t645 = -t415 * t506 - t718;
t642 = t161 * t512 + t418 * t82;
t641 = -t418 * t99 + t421 * t98;
t640 = -t27 * t413 + t28 * t415;
t129 = Ifges(5,4) * t237 - t447 * Ifges(5,2) + Ifges(5,6) * t305;
t602 = -t129 / 0.2e1;
t639 = -t180 * mrSges(4,3) + t413 * t602;
t563 = Ifges(3,4) * t419;
t617 = t414 ^ 2;
t638 = (t419 * (Ifges(3,1) * t422 - t563) / 0.2e1 - pkin(1) * (mrSges(3,1) * t419 + mrSges(3,2) * t422)) * t617;
t637 = t150 * Ifges(5,5) + t149 * Ifges(5,6) + t188 * Ifges(5,3) + t672;
t633 = t438 + t717;
t515 = qJD(2) * t422;
t492 = t414 * t515;
t268 = -qJD(3) * t351 + t421 * t492;
t516 = qJD(2) * t419;
t493 = t414 * t516;
t227 = t268 * t415 + t413 * t493;
t267 = qJD(3) * t352 + t418 * t492;
t338 = qJD(2) * t441;
t340 = t358 * qJD(2);
t137 = -t326 * t513 + t327 * t512 + t418 * t338 + t421 * t340;
t127 = (qJ(4) * t516 - qJD(4) * t422) * t414 + t137;
t341 = t519 * qJD(2);
t136 = t267 * pkin(3) - t268 * qJ(4) - t352 * qJD(4) + t341;
t66 = -t127 * t413 + t415 * t136;
t42 = pkin(4) * t267 - pkin(10) * t227 + t66;
t226 = -t268 * t413 + t415 * t493;
t67 = t415 * t127 + t413 * t136;
t59 = pkin(10) * t226 + t67;
t9 = -qJD(5) * t658 - t417 * t59 + t42 * t420;
t631 = mrSges(4,1) + t633;
t630 = t660 - t679;
t476 = mrSges(3,3) * t495;
t629 = -m(4) * t284 + mrSges(3,1) * t444 - mrSges(4,1) * t305 - mrSges(4,2) * t306 - t476;
t628 = t659 - t678;
t466 = t547 * t570;
t356 = -t419 * t466 + t422 * t571;
t497 = t414 * t570;
t275 = t356 * t421 + t418 * t497;
t627 = -g(1) * t275 - g(2) * t271 - g(3) * t352;
t626 = -mrSges(4,1) - t438;
t3 = t417 * t18 + t420 * t22 + t65 * t509 - t511 * t78;
t1 = qJ(6) * t179 + qJD(6) * t298 + t3;
t625 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t583 = -t298 / 0.2e1;
t598 = -t693 / 0.2e1;
t622 = -Ifges(6,2) * t599 + Ifges(7,3) * t600 + t583 * t674 + t598 * t676 - t695;
t619 = t48 * mrSges(6,1) + t7 * mrSges(7,1) + 0.2e1 * Ifges(7,3) * t611 - t57 * Ifges(6,4) / 0.2e1 - t179 * Ifges(6,6) / 0.2e1 + (t676 + Ifges(7,5)) * t613 + (t674 + Ifges(7,6)) * t592 + (-t612 + t611) * Ifges(6,2);
t618 = -mrSges(7,2) * t1 - mrSges(6,3) * t3 + t619;
t61 = t150 * Ifges(5,4) + t149 * Ifges(5,2) + t188 * Ifges(5,6);
t610 = t61 / 0.2e1;
t609 = Ifges(5,4) * t596 + t694;
t603 = Ifges(4,1) * t589 + Ifges(4,4) * t588 + Ifges(4,5) * t577;
t130 = Ifges(5,1) * t237 - t447 * Ifges(5,4) + Ifges(5,5) * t305;
t601 = t130 / 0.2e1;
t585 = -t237 / 0.2e1;
t581 = -t305 / 0.2e1;
t580 = t305 / 0.2e1;
t578 = t306 / 0.2e1;
t410 = t418 * pkin(9);
t565 = mrSges(6,3) * t142;
t564 = mrSges(6,3) * t693;
t562 = Ifges(3,4) * t422;
t561 = Ifges(4,4) * t418;
t560 = Ifges(4,4) * t421;
t559 = Ifges(5,4) * t413;
t558 = Ifges(5,4) * t415;
t542 = t353 * t413;
t355 = t419 * t571 + t422 * t466;
t540 = t355 * t413;
t527 = t418 * t422;
t119 = -mrSges(7,2) * t142 + mrSges(7,3) * t298;
t120 = -mrSges(6,2) * t298 - t565;
t525 = -t119 - t120;
t121 = mrSges(6,1) * t298 - t564;
t122 = -mrSges(7,1) * t298 + mrSges(7,2) * t693;
t524 = t121 - t122;
t367 = pkin(4) * t534 + t410;
t518 = t571 * pkin(1) + pkin(8) * t497;
t503 = t414 * t527;
t502 = t414 * t526;
t501 = Ifges(4,5) * t187 - Ifges(4,6) * t188 + Ifges(4,3) * t331;
t500 = Ifges(3,5) * t343 - Ifges(3,6) * t342 + Ifges(3,3) * t393;
t489 = t532 / 0.2e1;
t20 = t58 * mrSges(6,1) + t57 * mrSges(6,2);
t19 = t58 * mrSges(7,1) - t57 * mrSges(7,3);
t483 = t512 / 0.2e1;
t35 = -t179 * mrSges(7,1) + t57 * mrSges(7,2);
t87 = -t149 * mrSges(5,1) + t150 * mrSges(5,2);
t214 = -t418 * t326 + t327 * t421;
t475 = mrSges(3,3) * t494;
t468 = -pkin(1) * t570 + pkin(8) * t498;
t198 = pkin(3) * t531 - t214;
t465 = mrSges(4,1) * t351 + mrSges(4,2) * t352;
t460 = Ifges(4,1) * t421 - t561;
t459 = Ifges(5,1) * t415 - t559;
t458 = -Ifges(4,2) * t418 + t560;
t457 = -Ifges(5,2) * t413 + t558;
t456 = Ifges(4,5) * t421 - Ifges(4,6) * t418;
t455 = Ifges(5,5) * t415 - Ifges(5,6) * t413;
t37 = -t103 * t417 + t420 * t91;
t162 = t253 * t420 - t269 * t417;
t448 = t420 * t265 - t266 * t417;
t166 = t265 * t417 + t266 * t420;
t443 = t356 * pkin(2) + pkin(9) * t355 + t518;
t138 = -t326 * t512 - t327 * t513 + t338 * t421 - t418 * t340;
t8 = -t103 * t511 + t417 * t42 + t420 * t59 + t91 * t509;
t440 = t284 * (mrSges(4,1) * t418 + mrSges(4,2) * t421);
t270 = t354 * t418 + t421 * t498;
t274 = t356 * t418 - t421 * t497;
t436 = -g(1) * t274 - g(2) * t270 - g(3) * t351;
t147 = -pkin(4) * t265 + t198;
t433 = -t354 * pkin(2) - t353 * pkin(9) + t468;
t135 = -pkin(3) * t493 - t138;
t426 = Ifges(3,6) * t547 + (t422 * Ifges(3,2) + t563) * t414;
t424 = t414 * t444 * (Ifges(3,5) * t422 - Ifges(3,6) * t419);
t104 = -pkin(4) * t226 + t135;
t389 = Ifges(3,4) * t494;
t357 = (-mrSges(3,1) * t422 + mrSges(3,2) * t419) * t414;
t346 = t355 * pkin(2);
t344 = t353 * pkin(2);
t339 = t519 * qJD(1);
t335 = -mrSges(3,2) * t444 + t475;
t332 = t366 * t418;
t303 = -pkin(9) * t533 + t363;
t297 = Ifges(4,4) * t305;
t282 = Ifges(3,1) * t495 + Ifges(3,5) * t444 + t389;
t281 = Ifges(3,6) * qJD(2) + qJD(1) * t426;
t251 = t352 * t408 + t409 * t531;
t246 = mrSges(4,1) * t380 - mrSges(4,3) * t306;
t245 = -mrSges(4,2) * t380 - mrSges(4,3) * t305;
t244 = pkin(5) * t365 - qJ(6) * t366 - t406;
t199 = pkin(5) * t332 + qJ(6) * t333 + t367;
t194 = t275 * t409 + t355 * t408;
t193 = t275 * t408 - t355 * t409;
t174 = t306 * Ifges(4,1) - t297 + t662;
t172 = t306 * Ifges(4,5) - t305 * Ifges(4,6) + t380 * Ifges(4,3);
t168 = mrSges(5,1) * t305 - mrSges(5,3) * t237;
t167 = -t305 * mrSges(5,2) - mrSges(5,3) * t447;
t158 = pkin(5) * t421 - t162;
t157 = -qJ(6) * t421 + t647;
t153 = -mrSges(4,2) * t331 - mrSges(4,3) * t188;
t152 = mrSges(4,1) * t331 - mrSges(4,3) * t187;
t151 = mrSges(5,1) * t447 + t237 * mrSges(5,2);
t118 = mrSges(4,1) * t188 + mrSges(4,2) * t187;
t111 = t187 * Ifges(4,4) - t188 * Ifges(4,2) + t331 * Ifges(4,6);
t106 = mrSges(5,1) * t188 - mrSges(5,3) * t150;
t105 = -mrSges(5,2) * t188 + mrSges(5,3) * t149;
t94 = qJD(5) * t166 - t420 * t226 + t227 * t417;
t93 = qJD(5) * t448 + t226 * t417 + t227 * t420;
t85 = mrSges(6,1) * t142 + mrSges(6,2) * t693;
t84 = mrSges(7,1) * t142 - mrSges(7,3) * t693;
t83 = pkin(5) * t693 + qJ(6) * t142;
t63 = -pkin(5) * t448 - qJ(6) * t166 + t147;
t36 = -mrSges(6,2) * t179 - mrSges(6,3) * t58;
t34 = mrSges(6,1) * t179 - mrSges(6,3) * t57;
t33 = -mrSges(7,2) * t58 + mrSges(7,3) * t179;
t30 = -pkin(5) * t351 - t37;
t29 = qJ(6) * t351 + t658;
t16 = pkin(5) * t94 - qJ(6) * t93 - qJD(6) * t166 + t104;
t6 = -pkin(5) * t267 - t9;
t5 = qJ(6) * t267 + qJD(6) * t351 + t8;
t10 = [m(6) * (t104 * t126 + t147 * t48 + t25 * t9 + t26 * t8 + t3 * t658 + t37 * t4) + t658 * t36 - t618 * t448 + (t424 / 0.2e1 + t172 * t489) * qJD(2) + (Ifges(5,6) * t685 - Ifges(4,2) * t581 - Ifges(4,4) * t578 + Ifges(5,3) * t580 + t703 / 0.2e1 + t709 + t713) * t267 - t342 * t426 / 0.2e1 + t547 * t500 / 0.2e1 + t393 * (Ifges(3,3) * t547 + (Ifges(3,5) * t419 + Ifges(3,6) * t422) * t414) / 0.2e1 + (Ifges(5,5) * t266 + Ifges(5,6) * t265) * t587 + (Ifges(5,5) * t227 + Ifges(5,6) * t226) * t580 - t501 * t531 / 0.2e1 + t99 * (-mrSges(4,1) * t531 - t352 * mrSges(4,3)) + (t126 * mrSges(6,2) + t23 * mrSges(7,2) - t25 * mrSges(6,3) - t681 + t710) * t93 + t380 * (Ifges(4,5) * t268 + Ifges(4,3) * t493) / 0.2e1 + t638 * t508 - t357 * t505 + t692 * t166 + t531 * t680 + (Ifges(4,5) * t352 - Ifges(4,3) * t531) * t577 + (Ifges(4,4) * t268 + Ifges(4,6) * t493) * t581 + m(4) * (t137 * t181 + t138 * t180 + t211 * t325 + t214 * t99 + t215 * t98) + (-m(3) * t336 - t629) * t341 + m(7) * (t1 * t29 + t16 * t43 + t2 * t30 + t23 * t6 + t24 * t5 + t63 * t7) + m(5) * (t101 * t66 + t102 * t67 + t116 * t27 + t117 * t28 + t135 * t161 + t198 * t82) + t237 * (Ifges(5,1) * t227 + Ifges(5,4) * t226) / 0.2e1 + (-t181 * t493 + t268 * t284) * mrSges(4,2) + (-m(3) * t518 - mrSges(2,1) * t571 - t356 * mrSges(3,1) - t540 * mrSges(5,1) + mrSges(2,2) * t570 - mrSges(3,3) * t497 + t687 * t443 + t686 * (pkin(4) * t540 + t274 * t566 + t275 * t406 + t443) + t626 * t275 - t630 * t355 - t478 * t194 - t708 * t193 + t628 * t274) * g(2) + (-m(3) * t468 + mrSges(2,1) * t570 + t354 * mrSges(3,1) + t542 * mrSges(5,1) + mrSges(2,2) * t571 - mrSges(3,3) * t498 + t687 * t433 + t686 * (-pkin(4) * t542 - t270 * t566 - t271 * t406 + t433) - t626 * t271 + t630 * t353 + t478 * t696 + t708 * t189 - t628 * t270) * g(1) + (Ifges(5,4) * t266 + Ifges(5,2) * t265) * t596 + (Ifges(4,4) * t352 - Ifges(4,6) * t531) * t588 + (Ifges(4,1) * t352 - Ifges(4,5) * t531) * t589 + (Ifges(5,4) * t227 + Ifges(5,2) * t226) * t685 + t343 * (Ifges(3,5) * t547 + (t419 * Ifges(3,1) + t562) * t414) / 0.2e1 - t281 * t493 / 0.2e1 + t180 * (mrSges(4,1) * t493 - mrSges(4,3) * t268) + (-t101 * t227 + t102 * t226 + t265 * t28 - t266 * t27) * mrSges(5,3) + (-t98 * mrSges(4,3) + Ifges(6,6) * t612 + Ifges(7,6) * t611 + Ifges(5,5) * t595 + Ifges(5,6) * t596 + Ifges(5,3) * t587 - Ifges(4,2) * t588 - Ifges(4,4) * t589 - t111 / 0.2e1 - t28 * mrSges(5,2) + t27 * mrSges(5,1) - Ifges(4,6) * t577 + t675 * t613 + t673 * t592 + t625 + t637 / 0.2e1) * t351 + (Ifges(4,1) * t268 + Ifges(4,5) * t493) * t578 + (Ifges(5,1) * t266 + Ifges(5,4) * t265) * t595 + t227 * t601 + t352 * t603 + t266 * t609 + t265 * t610 + Ifges(2,3) * qJDD(1) + (t241 * t547 - t342 * t569 + t358 * t393) * mrSges(3,1) + (Ifges(3,4) * t343 - Ifges(3,2) * t342 + Ifges(3,6) * t393) * t531 / 0.2e1 + (Ifges(3,1) * t343 - Ifges(3,4) * t342 + Ifges(3,5) * t393) * t489 + (t240 * t531 - t241 * t532 - t336 * t492 - t339 * t493 - t342 * t519 - t343 * t358) * mrSges(3,3) + m(3) * (pkin(1) ^ 2 * qJDD(1) * t617 + t240 * t519 + t241 * t358 + t339 * t340) + (-t240 * t547 - t343 * t569 - t393 * t519) * mrSges(3,2) + t340 * t335 + (t617 * qJD(1) * (-Ifges(3,2) * t419 + t562) + t414 * t282) * t515 / 0.2e1 + t211 * t465 + t325 * t118 + t268 * t174 / 0.2e1 + t82 * (-mrSges(5,1) * t265 + mrSges(5,2) * t266) + t137 * t245 + t138 * t246 + t161 * (-mrSges(5,1) * t226 + mrSges(5,2) * t227) + t226 * t129 / 0.2e1 + t214 * t152 + t215 * t153 + t198 * t87 + t67 * t167 + t66 * t168 + t147 * t20 + t135 * t151 + t5 * t119 + t8 * t120 + t9 * t121 + t6 * t122 + t117 * t105 + t116 * t106 + t104 * t85 + t16 * t84 + (mrSges(6,1) * t126 + t697) * t94 + t29 * t33 + t30 * t35 + t37 * t34 + t63 * t19; t27 * (-mrSges(5,1) * t421 - mrSges(5,3) * t530) + (t126 * t644 + t162 * t4 + t25 * t666 + t26 * t667 + t3 * t647 + t367 * t48) * m(6) + t647 * t36 + t657 * t84 + (-t407 - t228) * t246 + (-t180 * (mrSges(4,1) * t419 - mrSges(4,3) * t526) - t181 * (-mrSges(4,2) * t419 - mrSges(4,3) * t527)) * t517 + (-t554 + t713) * t513 + t666 * t121 + t667 * t120 + t668 * t122 + (t1 * t157 + t158 * t2 + t199 * t7 + t23 * t668 + t24 * t669 + t43 * t657) * m(7) + t669 * t119 - t440 * t494 + t618 * t332 + (-mrSges(6,1) * t643 + mrSges(6,3) * t650) * t25 + (-t1 * t421 - t24 * t643 + t333 * t7 + t43 * t650) * mrSges(7,3) + (mrSges(7,1) * t643 - mrSges(7,2) * t650) * t23 - t671 * t333 / 0.2e1 + (t184 * t675 + t472 * t673) * t583 + (-t333 * t675 - t421 * t673) * t592 + (-t152 + t87) * t410 - t637 * t421 / 0.2e1 + (Ifges(7,5) * t184 + Ifges(7,6) * t472) * t600 + (t26 * t643 + t3 * t421 - t333 * t48) * mrSges(6,2) + t644 * t85 + t645 * t167 + t646 * t168 + (pkin(9) * t642 + t101 * t646 + t102 * t645 - t161 * t202 + t27 * t303 + t28 * t304) * m(5) + (t476 + t629) * t339 + t174 * t483 + (-t424 / 0.2e1 - t638 * qJD(1)) * qJD(1) - t61 * t534 / 0.2e1 + t28 * (mrSges(5,2) * t421 - mrSges(5,3) * t534) + (-t335 + t475) * t336 + (mrSges(6,1) * t712 - t650 * mrSges(6,2)) * t126 + (t357 + t687 * t520 - t678 * t503 + t686 * (t406 * t502 + t503 * t566 + t532 * t568 + t520) - t708 * (t408 * t502 - t409 * t532) + (t702 * t422 + (-t462 - mrSges(4,3)) * t419 - t478 * (t408 * t419 + t409 * t526)) * t414) * g(3) + (-t506 - t229) * t245 + (t184 * t677 + t472 * t675) * t598 + (-t333 * t677 - t421 * t675) * t613 + (m(4) * (-t180 * t421 - t181 * t418) * pkin(9) + (Ifges(5,6) * t418 + t421 * t457) * t685 + t101 * (mrSges(5,1) * t418 - mrSges(5,3) * t529) + t102 * (-mrSges(5,2) * t418 - mrSges(5,3) * t533) + t440) * qJD(3) - ((-Ifges(3,2) * t495 + t421 * t174 + t418 * t635 + t282 + t389) * t422 + t380 * (Ifges(4,3) * t419 + t422 * t456) + t306 * (Ifges(4,5) * t419 + t422 * t460) + t419 * t172) * t517 / 0.2e1 + t184 * t706 + (Ifges(5,3) * t418 + t421 * t455) * t514 / 0.2e1 - t458 * t514 / 0.2e1 + t710 * t232 + t711 * t151 + (-t292 / 0.2e1 + t415 * t483) * t130 + t639 * t512 + (-pkin(2) * t211 + pkin(9) * t641 - t180 * t228 - t181 * t229) * m(4) + t641 * mrSges(4,3) + t642 * t462 - t102 * (-mrSges(5,2) * t472 + t291 * mrSges(5,3)) + t447 * (Ifges(5,4) * t292 + Ifges(5,2) * t291 + Ifges(5,6) * t472) / 0.2e1 - t101 * (mrSges(5,1) * t472 - t292 * mrSges(5,3)) + t173 * t472 / 0.2e1 + (Ifges(6,4) * t184 + Ifges(6,6) * t472) * t599 + t421 * pkin(9) * t153 + t291 * t602 + t418 * t603 + t530 * t609 + (-Ifges(5,5) * t421 + t418 * t459) * t595 + (-Ifges(5,6) * t421 + t418 * t457) * t596 + (Ifges(5,5) * t292 + Ifges(5,6) * t291 + Ifges(5,3) * t472) * t581 + (Ifges(5,1) * t292 + Ifges(5,4) * t291 + Ifges(5,5) * t472) * t585 + (-Ifges(5,3) * t421 + t418 * t455) * t587 + (Ifges(4,2) * t421 + t561) * t588 + (Ifges(4,1) * t418 + t560) * t589 + t421 * t111 / 0.2e1 + t367 * t20 + (-Ifges(7,5) * t333 - Ifges(7,6) * t421) * t611 + (-Ifges(6,4) * t333 - Ifges(6,6) * t421) * t612 + t4 * (-mrSges(6,1) * t421 + mrSges(6,3) * t333) + t2 * (mrSges(7,1) * t421 - mrSges(7,2) * t333) + (t305 * (Ifges(4,6) * t419 + t422 * t458) + t419 * t281) * t517 / 0.2e1 + (t237 * (Ifges(5,5) * t418 + t421 * t459) + t380 * t456 + t306 * t460) * qJD(3) / 0.2e1 - t211 * t464 + t303 * t106 + t304 * t105 - t161 * (-mrSges(5,1) * t291 + mrSges(5,2) * t292) + t500 - t240 * mrSges(3,2) + t241 * mrSges(3,1) + t199 * t19 + t162 * t34 + t157 * t33 + t158 * t35 - pkin(2) * t118 + t697 * t233 + (t622 + t699) * t183 + (m(5) * t346 + t704 * (pkin(9) * t356 - t346) + t691 * t356 + t689 * t355) * g(1) + (m(5) * t344 + t704 * (pkin(9) * t354 - t344) + t691 * t354 + t689 * t353) * g(2) + (Ifges(4,5) * t418 + Ifges(4,6) * t421) * t577; (-t126 * t148 + t25 * t653 + t26 * t654 + t280 * t3 + t4 * t446 - t406 * t48) * m(6) + (t1 * t280 - t2 * t446 + t23 * t652 + t24 * t655 + t244 * t7 + t43 * t656) * m(7) - (t35 - t34) * t446 + (t686 * (-t274 * t406 + t275 * t566) + t659 * t275 + t631 * t274) * g(1) + (t686 * (-t270 * t406 + t271 * t566) + t659 * t271 + t631 * t270) * g(2) + (-t352 * t488 + t465 + t686 * (-t351 * t406 + t352 * t566) + t633 * t351) * g(3) + t656 * t84 + t619 * t365 + t621 * t350 - t622 * t185 + (mrSges(6,1) * t649 - mrSges(6,2) * t648) * t126 + (t33 + t36) * t280 - (-Ifges(4,1) * t305 - t551 + t635) * t306 / 0.2e1 + t652 * t122 + t653 * t121 + t654 * t120 + t655 * t119 + t692 * t366 + (t25 * t648 - t26 * t649 - t3 * t365 + t627) * mrSges(6,3) + (Ifges(5,2) * t596 + Ifges(5,6) * t587 + qJ(4) * t105 + qJD(4) * t167 + t610) * t415 + (t246 - t151) * t181 + (-t297 + t174) * t580 + (t673 * t583 + t675 * t598 + Ifges(7,6) * t600 + Ifges(6,6) * t599 + Ifges(5,3) * t581 + Ifges(5,5) * t585 - Ifges(4,2) * t580 + t665 / 0.2e1 - t709) * t306 + t559 * t596 + (-qJ(4) * t106 - qJD(4) * t168 + t609 + t694) * t413 + t558 * t595 + (t161 * t462 - t455 * t581 - t459 * t585 + t662 / 0.2e1 + t457 * t685 + t284 * mrSges(4,2) + t639) * t305 + (-t101 * t545 - t102 * t546 + t640) * mrSges(5,3) + (-pkin(3) * t82 + (-t101 * t413 + t102 * t415) * qJD(4) + t640 * qJ(4) - t101 * t123 - t102 * t124 - t161 * t181) * m(5) + t545 * t601 - t680 - t406 * t20 + t82 * t463 + t501 + t244 * t19 - t180 * t245 - t124 * t167 - t123 * t168 - t148 * t85 + t99 * mrSges(4,1) - pkin(3) * t87 + (-t1 * t365 - t23 * t648 - t24 * t649 + t627) * mrSges(7,2) + (Ifges(6,4) * t599 + Ifges(7,5) * t600 + t675 * t583 + t677 * t598 + t698) * t186 + (-t623 + t698) * t349 + t173 * t578; t524 * t693 - t525 * t142 + t447 * t167 + t237 * t168 + t19 + t20 + t87 + (t142 * t24 - t23 * t693 + t436 + t7) * m(7) + (t142 * t26 + t25 * t693 + t436 + t48) * m(6) + (t101 * t237 + t102 * t447 + t436 + t82) * m(5); (-t142 * t677 + t139 - t557 + t68) * t598 + (-t708 * (t352 * t409 - t408 * t531) + t478 * t251) * g(3) + t625 + (t193 * t478 - t194 * t708) * g(1) + (t478 * t189 - t696 * t708) * g(2) + (-pkin(5) * t2 + qJ(6) * t1 - t23 * t26 + t24 * t651 - t43 * t83) * m(7) + (-t142 * t675 + t674 * t693) * t583 + (t524 + t564) * t26 + (t525 - t565) * t25 + (-Ifges(6,2) * t693 - t140 + t670) * t599 + (Ifges(7,3) * t693 - t556) * t600 + t71 * t597 + (t142 * t23 + t24 * t693) * mrSges(7,2) - t43 * (mrSges(7,1) * t693 + mrSges(7,3) * t142) - t126 * (mrSges(6,1) * t693 - mrSges(6,2) * t142) + qJD(6) * t119 - t83 * t84 + qJ(6) * t33 - pkin(5) * t35 + t672; -t298 * t119 + t693 * t84 + (-g(1) * t193 - g(2) * t189 - g(3) * t251 - t24 * t298 + t43 * t693 + t2) * m(7) + t35;];
tau  = t10;
