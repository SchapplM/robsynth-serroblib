% Calculate vector of inverse dynamics joint torques for
% S6RRRPRR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRR10_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:15:40
% EndTime: 2019-03-09 19:16:58
% DurationCPUTime: 45.66s
% Computational Cost: add. (16397->1005), mult. (34446->1294), div. (0->0), fcn. (23700->12), ass. (0->455)
t375 = sin(qJ(2));
t361 = t375 * pkin(8);
t380 = cos(qJ(2));
t366 = t380 * pkin(2);
t484 = -pkin(1) - t366;
t410 = t484 - t361;
t274 = t410 * qJD(1);
t357 = t380 * qJD(1);
t352 = pkin(7) * t357;
t319 = qJD(2) * pkin(8) + t352;
t374 = sin(qJ(3));
t379 = cos(qJ(3));
t196 = t379 * t274 - t374 * t319;
t643 = qJD(4) - t196;
t604 = pkin(8) - pkin(9);
t321 = t604 * t379;
t483 = -pkin(7) * t374 - pkin(3);
t518 = t379 * t380;
t395 = -pkin(9) * t518 + (-pkin(4) + t483) * t375;
t446 = pkin(2) * t375 - pkin(8) * t380;
t304 = t446 * qJD(1);
t529 = t304 * t379;
t727 = -qJD(1) * t395 + qJD(3) * t321 + t529;
t266 = t374 * t304;
t505 = qJD(1) * t375;
t341 = qJ(4) * t505;
t501 = qJD(3) * t374;
t522 = t375 * t379;
t523 = t374 * t380;
t726 = t266 + t341 + (-pkin(7) * t522 + pkin(9) * t523) * qJD(1) + t604 * t501;
t480 = t379 * t505;
t293 = qJD(2) * t374 + t480;
t725 = -pkin(9) * t293 + t643;
t718 = Ifges(4,1) + Ifges(5,1);
t717 = Ifges(5,4) + Ifges(4,5);
t716 = Ifges(4,6) - Ifges(5,6);
t378 = cos(qJ(5));
t373 = sin(qJ(5));
t528 = t373 * t374;
t411 = t378 * t379 + t528;
t641 = qJD(3) - qJD(5);
t204 = t641 * t411;
t402 = t380 * t411;
t234 = qJD(1) * t402;
t724 = t204 - t234;
t524 = t374 * t378;
t527 = t373 * t379;
t412 = -t524 + t527;
t205 = t641 * t412;
t401 = t412 * t380;
t233 = qJD(1) * t401;
t710 = t205 - t233;
t333 = -t357 + qJD(3);
t327 = qJD(5) - t333;
t314 = qJD(6) + t327;
t723 = t333 / 0.2e1;
t686 = mrSges(6,3) + mrSges(7,3);
t197 = t374 * t274 + t379 * t319;
t482 = t374 * t505;
t496 = t379 * qJD(2);
t292 = t482 - t496;
t143 = pkin(9) * t292 + t197;
t383 = -pkin(3) - pkin(4);
t310 = -qJ(4) * t373 + t378 * t383;
t675 = qJD(5) * t310 - t373 * t143 + t378 * t725;
t311 = t378 * qJ(4) + t373 * t383;
t674 = -qJD(5) * t311 - t378 * t143 - t373 * t725;
t320 = t604 * t374;
t216 = t373 * t320 + t378 * t321;
t673 = -qJD(5) * t216 + t726 * t373 + t378 * t727;
t497 = qJD(5) * t378;
t498 = qJD(5) * t373;
t672 = t320 * t497 - t321 * t498 + t373 * t727 - t726 * t378;
t376 = sin(qJ(1));
t381 = cos(qJ(1));
t516 = t381 * t374;
t262 = t376 * t518 - t516;
t519 = t376 * t380;
t261 = t374 * t519 + t379 * t381;
t534 = t261 * t373;
t416 = t262 * t378 + t534;
t369 = qJ(5) + qJ(6);
t358 = sin(t369);
t359 = cos(t369);
t418 = t261 * t358 + t262 * t359;
t650 = t261 * t359 - t262 * t358;
t515 = mrSges(7,1) * t650 - t418 * mrSges(7,2);
t722 = t416 * mrSges(6,2) - t515;
t582 = t327 / 0.2e1;
t584 = t314 / 0.2e1;
t415 = t292 * t373 + t293 * t378;
t593 = t415 / 0.2e1;
t190 = t292 * t378 - t293 * t373;
t595 = t190 / 0.2e1;
t372 = sin(qJ(6));
t377 = cos(qJ(6));
t112 = t190 * t372 + t377 * t415;
t600 = t112 / 0.2e1;
t456 = t377 * t190 - t415 * t372;
t602 = t456 / 0.2e1;
t721 = Ifges(6,5) * t593 + Ifges(7,5) * t600 + Ifges(6,6) * t595 + Ifges(7,6) * t602 + Ifges(6,3) * t582 + Ifges(7,3) * t584;
t495 = qJD(1) * qJD(2);
t309 = qJDD(1) * t375 + t380 * t495;
t502 = qJD(3) * t292;
t184 = qJDD(2) * t374 + t309 * t379 - t502;
t599 = t184 / 0.2e1;
t185 = qJD(3) * t293 - t379 * qJDD(2) + t309 * t374;
t597 = t185 / 0.2e1;
t308 = t380 * qJDD(1) - t375 * t495;
t290 = qJDD(3) - t308;
t590 = t290 / 0.2e1;
t588 = t292 / 0.2e1;
t720 = -t293 / 0.2e1;
t719 = mrSges(5,2) + mrSges(4,3);
t715 = -Ifges(4,3) - Ifges(5,2);
t714 = pkin(5) * t505 - pkin(10) * t724 + t673;
t713 = pkin(10) * t710 - t672;
t689 = pkin(10) * t415;
t712 = -t689 + t675;
t702 = pkin(10) * t190;
t711 = -t702 + t674;
t660 = -t374 * t716 + t379 * t717;
t556 = Ifges(5,5) * t374;
t560 = Ifges(4,4) * t374;
t657 = t379 * t718 + t556 - t560;
t499 = qJD(3) * t379;
t708 = -t374 * qJD(4) - t352 + (t357 * t379 - t499) * qJ(4);
t64 = qJD(5) * t190 + t184 * t378 + t185 * t373;
t65 = -qJD(5) * t415 - t184 * t373 + t185 * t378;
t21 = qJD(6) * t456 + t372 * t65 + t377 * t64;
t22 = -qJD(6) * t112 - t372 * t64 + t377 * t65;
t273 = qJDD(5) - t290;
t258 = qJDD(6) + t273;
t493 = Ifges(7,5) * t21 + Ifges(7,6) * t22 + Ifges(7,3) * t258;
t116 = t333 * t383 + t725;
t323 = t333 * qJ(4);
t123 = t143 + t323;
t536 = qJDD(1) * pkin(1);
t199 = -pkin(2) * t308 - pkin(8) * t309 - t536;
t288 = t308 * pkin(7);
t255 = qJDD(2) * pkin(8) + t288;
t85 = t199 * t379 - t374 * t255 - t274 * t501 - t319 * t499;
t400 = qJDD(4) - t85;
t49 = -pkin(9) * t184 + t290 * t383 + t400;
t84 = t374 * t199 + t379 * t255 + t274 * t499 - t319 * t501;
t66 = t290 * qJ(4) + t333 * qJD(4) + t84;
t52 = pkin(9) * t185 + t66;
t11 = t116 * t497 - t123 * t498 + t373 * t49 + t378 * t52;
t10 = pkin(10) * t65 + t11;
t58 = t378 * t116 - t123 * t373;
t45 = t58 - t689;
t41 = pkin(5) * t327 + t45;
t59 = t116 * t373 + t123 * t378;
t46 = t59 + t702;
t546 = t372 * t46;
t13 = t377 * t41 - t546;
t12 = -qJD(5) * t59 - t373 * t52 + t378 * t49;
t9 = pkin(5) * t273 - pkin(10) * t64 + t12;
t2 = qJD(6) * t13 + t10 * t377 + t372 * t9;
t541 = t377 * t46;
t14 = t372 * t41 + t541;
t3 = -qJD(6) * t14 - t10 * t372 + t377 * t9;
t631 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t408 = t493 + t631;
t492 = Ifges(6,5) * t64 + Ifges(6,6) * t65 + Ifges(6,3) * t273;
t585 = -t314 / 0.2e1;
t601 = -t112 / 0.2e1;
t603 = -t456 / 0.2e1;
t104 = Ifges(7,4) * t456;
t44 = Ifges(7,1) * t112 + Ifges(7,5) * t314 + t104;
t612 = -t44 / 0.2e1;
t557 = Ifges(7,4) * t112;
t43 = Ifges(7,2) * t456 + Ifges(7,6) * t314 + t557;
t614 = -t43 / 0.2e1;
t630 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t351 = pkin(7) * t505;
t318 = -qJD(2) * pkin(2) + t351;
t166 = t292 * pkin(3) - t293 * qJ(4) + t318;
t131 = -pkin(4) * t292 - t166;
t86 = -pkin(5) * t190 + t131;
t634 = mrSges(7,2) * t86 - mrSges(7,3) * t13;
t635 = -mrSges(7,1) * t86 + mrSges(7,3) * t14;
t707 = t408 + t492 + t630 + (Ifges(7,1) * t601 + Ifges(7,4) * t603 + Ifges(7,5) * t585 + t612 - t634) * t456 - (Ifges(7,4) * t601 + Ifges(7,2) * t603 + Ifges(7,6) * t585 + t614 - t635) * t112;
t706 = m(6) * pkin(9) + t686 - t719;
t413 = t358 * t374 + t359 * t379;
t414 = t358 * t379 - t359 * t374;
t439 = t379 * mrSges(5,1) + t374 * mrSges(5,3);
t441 = mrSges(4,1) * t379 - mrSges(4,2) * t374;
t694 = t411 * mrSges(6,1) - t412 * mrSges(6,2);
t705 = t413 * mrSges(7,1) - t414 * mrSges(7,2) + t439 + t441 + t694;
t598 = -t185 / 0.2e1;
t690 = t308 / 0.2e1;
t704 = t309 / 0.2e1;
t703 = t717 * t590 + (-Ifges(4,4) + Ifges(5,5)) * t597 + t718 * t599;
t283 = Ifges(4,4) * t292;
t548 = t292 * Ifges(5,5);
t699 = t293 * t718 + t333 * t717 - t283 + t548;
t481 = t374 * t357;
t487 = t383 * t374;
t666 = qJD(3) * t487 - t383 * t481 - t708;
t263 = -t376 * t379 + t380 * t516;
t517 = t380 * t381;
t264 = t374 * t376 + t379 * t517;
t173 = t263 * t378 - t264 * t373;
t174 = t263 * t373 + t264 * t378;
t155 = t263 * t359 - t264 * t358;
t156 = t263 * t358 + t264 * t359;
t514 = t155 * mrSges(7,1) - t156 * mrSges(7,2);
t698 = mrSges(6,1) * t173 - mrSges(6,2) * t174 + t514;
t248 = t373 * t522 - t375 * t524;
t249 = t411 * t375;
t459 = t248 * mrSges(6,1) + t249 * mrSges(6,2);
t513 = (-mrSges(7,1) * t414 - mrSges(7,2) * t413) * t375;
t697 = -t513 + t459;
t363 = t379 * pkin(4);
t360 = t374 * qJ(4);
t639 = -t379 * pkin(3) - pkin(2) - t360;
t696 = t639 - t363;
t562 = Ifges(3,4) * t375;
t679 = Ifges(3,2) * t380;
t432 = t562 + t679;
t695 = t58 * mrSges(6,1) + t13 * mrSges(7,1) - t59 * mrSges(6,2) - t14 * mrSges(7,2) + Ifges(3,6) * qJD(2) / 0.2e1 + qJD(1) * t432 / 0.2e1 + t715 * t723 + t717 * t720 + t716 * t588 + t721;
t618 = m(7) * pkin(5);
t617 = t21 / 0.2e1;
t616 = t22 / 0.2e1;
t693 = -t64 * Ifges(6,4) / 0.2e1 - t65 * Ifges(6,2) / 0.2e1 - t273 * Ifges(6,6) / 0.2e1;
t610 = t64 / 0.2e1;
t609 = t65 / 0.2e1;
t692 = -m(6) - m(5);
t691 = -m(6) - m(7);
t592 = t258 / 0.2e1;
t591 = t273 / 0.2e1;
t571 = pkin(5) * t415;
t503 = qJD(2) * t380;
t490 = pkin(7) * t503;
t688 = -mrSges(4,2) + mrSges(5,3);
t687 = -mrSges(3,3) + mrSges(2,2);
t303 = -pkin(5) + t310;
t201 = t303 * t372 + t311 * t377;
t684 = -qJD(6) * t201 - t372 * t712 + t377 * t711;
t200 = t303 * t377 - t311 * t372;
t683 = qJD(6) * t200 + t372 * t711 + t377 * t712;
t215 = t378 * t320 - t321 * t373;
t161 = pkin(10) * t412 + t215;
t162 = -pkin(10) * t411 + t216;
t90 = t161 * t372 + t162 * t377;
t682 = -qJD(6) * t90 + t372 * t713 + t377 * t714;
t89 = t161 * t377 - t162 * t372;
t681 = qJD(6) * t89 + t372 * t714 - t377 * t713;
t443 = mrSges(3,1) * t380 - mrSges(3,2) * t375;
t676 = mrSges(2,1) + t443;
t670 = pkin(5) * t710 + t666;
t508 = t366 + t361;
t313 = -pkin(1) - t508;
t334 = pkin(7) * t523;
t365 = t380 * pkin(3);
t567 = pkin(9) * t375;
t183 = pkin(4) * t380 + t334 + t365 + (-t313 - t567) * t379;
t335 = pkin(7) * t518;
t232 = t374 * t313 + t335;
t217 = -qJ(4) * t380 + t232;
t525 = t374 * t375;
t195 = pkin(9) * t525 + t217;
t101 = t373 * t183 + t378 * t195;
t296 = t372 * t378 + t373 * t377;
t668 = t314 * t296;
t294 = -t372 * t373 + t377 * t378;
t667 = t314 * t294;
t665 = (-t481 + t501) * pkin(3) + t708;
t664 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t292 + mrSges(4,2) * t293 + mrSges(3,3) * t505;
t663 = -t375 * t715 + t380 * t660;
t662 = t375 * t717 + t380 * t657;
t540 = t379 * mrSges(5,3);
t438 = t374 * mrSges(5,1) - t540;
t440 = mrSges(4,1) * t374 + mrSges(4,2) * t379;
t661 = -t166 * t438 - t318 * t440;
t659 = t374 * t717 + t379 * t716;
t555 = Ifges(5,5) * t379;
t559 = Ifges(4,4) * t379;
t658 = t374 * t718 - t555 + t559;
t654 = t184 * t717 - t185 * t716 - t290 * t715;
t282 = Ifges(5,5) * t293;
t167 = t333 * Ifges(5,6) + t292 * Ifges(5,3) + t282;
t349 = Ifges(3,4) * t357;
t653 = Ifges(3,1) * t505 + Ifges(3,5) * qJD(2) + t374 * t167 + t349;
t289 = t309 * pkin(7);
t652 = t288 * t380 + t289 * t375;
t651 = t261 * t378 - t262 * t373;
t649 = t719 * t375;
t648 = -t374 * t85 + t379 * t84;
t69 = -pkin(3) * t290 + t400;
t647 = t374 * t69 + t379 * t66;
t160 = t323 + t197;
t646 = -t160 * mrSges(5,2) - t197 * mrSges(4,3);
t158 = -pkin(3) * t333 + t643;
t645 = t158 * mrSges(5,2) - t196 * mrSges(4,3);
t642 = -m(5) + t691;
t637 = -Ifges(5,5) * t184 / 0.2e1 - Ifges(5,6) * t290 / 0.2e1 + Ifges(4,4) * t599 + Ifges(4,6) * t590 + (Ifges(5,3) + Ifges(4,2)) * t598;
t633 = -m(7) * (pkin(5) * t373 + qJ(4)) - t688;
t632 = pkin(8) * (-m(4) + t642);
t382 = -pkin(10) - pkin(9);
t346 = pkin(5) * t378 + pkin(4);
t407 = pkin(5) * t528 + t346 * t379;
t442 = mrSges(3,1) * t375 + mrSges(3,2) * t380;
t627 = t442 + (-m(7) * t382 + t706) * t380 + (-m(7) * (-t407 + t639) - m(6) * t696 - m(5) * t639 + m(4) * pkin(2) + t705) * t375;
t626 = m(6) * pkin(4) + m(7) * t346 + mrSges(4,1) + mrSges(5,1);
t625 = -t85 * mrSges(4,1) + t69 * mrSges(5,1) + t84 * mrSges(4,2) - t66 * mrSges(5,3);
t620 = Ifges(7,4) * t617 + Ifges(7,2) * t616 + Ifges(7,6) * t592;
t619 = Ifges(7,1) * t617 + Ifges(7,4) * t616 + Ifges(7,5) * t592;
t615 = Ifges(6,1) * t610 + Ifges(6,4) * t609 + Ifges(6,5) * t591;
t613 = t43 / 0.2e1;
t611 = t44 / 0.2e1;
t558 = Ifges(6,4) * t415;
t92 = Ifges(6,2) * t190 + Ifges(6,6) * t327 + t558;
t608 = -t92 / 0.2e1;
t607 = t92 / 0.2e1;
t187 = Ifges(6,4) * t190;
t93 = Ifges(6,1) * t415 + Ifges(6,5) * t327 + t187;
t606 = -t93 / 0.2e1;
t605 = t93 / 0.2e1;
t596 = -t190 / 0.2e1;
t594 = -t415 / 0.2e1;
t589 = -t292 / 0.2e1;
t586 = t293 / 0.2e1;
t583 = -t327 / 0.2e1;
t581 = -t333 / 0.2e1;
t576 = m(7) * t375;
t575 = mrSges(6,3) * t58;
t574 = mrSges(6,3) * t59;
t570 = pkin(7) * t375;
t566 = -qJD(1) / 0.2e1;
t565 = qJD(3) / 0.2e1;
t564 = mrSges(6,3) * t190;
t563 = mrSges(6,3) * t415;
t561 = Ifges(3,4) * t380;
t547 = t293 * Ifges(4,4);
t521 = t375 * t381;
t520 = t375 * t382;
t307 = t446 * qJD(2);
t512 = t374 * t307 + t313 * t499;
t208 = t293 * pkin(3) + t292 * qJ(4);
t479 = t380 * t496;
t510 = qJ(4) * t479 + qJD(4) * t522;
t507 = t381 * pkin(1) + t376 * pkin(7);
t504 = qJD(2) * t375;
t500 = qJD(3) * t375;
t491 = pkin(7) * t504;
t489 = pkin(8) * t501;
t488 = pkin(8) * t499;
t256 = -qJDD(2) * pkin(2) + t289;
t478 = t374 * t500;
t170 = -t292 * Ifges(4,2) + t333 * Ifges(4,6) + t547;
t477 = -t374 * t170 / 0.2e1;
t467 = t503 / 0.2e1;
t464 = -t500 / 0.2e1;
t463 = t499 / 0.2e1;
t462 = -t357 / 0.2e1;
t460 = t495 / 0.2e1;
t127 = -t290 * mrSges(5,1) + t184 * mrSges(5,2);
t100 = t378 * t183 - t195 * t373;
t231 = t313 * t379 - t334;
t453 = pkin(3) * t518 + qJ(4) * t523 + t508;
t452 = pkin(2) * t517 + pkin(8) * t521 + t507;
t448 = -pkin(7) + t487;
t147 = -pkin(4) * t293 - t208;
t447 = t375 * t483;
t444 = qJD(3) * t335 - t307 * t379 + t313 * t501;
t431 = -Ifges(4,2) * t374 + t559;
t430 = Ifges(4,2) * t379 + t560;
t427 = Ifges(3,5) * t380 - Ifges(3,6) * t375;
t424 = Ifges(5,3) * t374 + t555;
t423 = -Ifges(5,3) * t379 + t556;
t72 = pkin(5) * t380 - pkin(10) * t249 + t100;
t75 = -pkin(10) * t248 + t101;
t37 = -t372 * t75 + t377 * t72;
t38 = t372 * t72 + t377 * t75;
t421 = t264 * pkin(3) + t452;
t144 = -mrSges(6,2) * t327 + t564;
t145 = mrSges(6,1) * t327 - t563;
t420 = t144 * t378 - t145 * t373;
t153 = -t248 * t377 - t249 * t372;
t419 = t248 * t372 - t249 * t377;
t193 = t372 * t412 - t377 * t411;
t194 = -t372 * t411 - t377 * t412;
t220 = -pkin(7) * t480 + t266;
t68 = t185 * pkin(3) - t184 * qJ(4) - t293 * qJD(4) + t256;
t406 = pkin(1) * t442;
t403 = t375 * (Ifges(3,1) * t380 - t562);
t105 = pkin(9) * t478 + qJD(2) * t395 + t444;
t343 = qJ(4) * t504;
t106 = t343 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t522 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t374) * t380 + t512;
t34 = t373 * t105 + t378 * t106 + t183 * t497 - t195 * t498;
t330 = qJ(4) * t522;
t214 = t375 * t448 + t330;
t399 = -t478 + t479;
t398 = t374 * t503 + t375 * t499;
t396 = -g(1) * t263 - g(2) * t261 - g(3) * t525;
t53 = -pkin(4) * t185 - t68;
t391 = Ifges(4,6) * t375 + t380 * t431;
t390 = Ifges(5,6) * t375 + t380 * t424;
t35 = -qJD(5) * t101 + t378 * t105 - t106 * t373;
t138 = (-t375 * t496 - t380 * t501) * pkin(7) + t512;
t117 = (t379 * t383 - t360) * t500 + t448 * t503 + t510;
t367 = t381 * pkin(7);
t316 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t357;
t276 = t440 * t375;
t253 = t263 * pkin(3);
t251 = t261 * pkin(3);
t241 = -t330 + (pkin(3) * t374 + pkin(7)) * t375;
t224 = -mrSges(5,2) * t292 + mrSges(5,3) * t333;
t223 = -mrSges(5,1) * t333 + mrSges(5,2) * t293;
t222 = mrSges(4,1) * t333 - mrSges(4,3) * t293;
t221 = -mrSges(4,2) * t333 - mrSges(4,3) * t292;
t219 = pkin(7) * t482 + t529;
t218 = -t231 + t365;
t211 = pkin(5) * t411 - t696;
t209 = mrSges(5,1) * t292 - mrSges(5,3) * t293;
t207 = qJD(1) * t447 - t529;
t206 = t220 + t341;
t146 = pkin(5) * t248 + t214;
t139 = t374 * t491 - t444;
t137 = pkin(3) * t398 + qJ(4) * t478 + t490 - t510;
t136 = -t233 * t372 + t234 * t377;
t135 = -t233 * t377 - t234 * t372;
t130 = qJD(2) * t447 + t444;
t129 = -mrSges(5,2) * t185 + mrSges(5,3) * t290;
t128 = -mrSges(4,2) * t290 - mrSges(4,3) * t185;
t126 = mrSges(4,1) * t290 - mrSges(4,3) * t184;
t122 = -qJD(4) * t380 + t138 + t343;
t119 = qJD(2) * t402 + t205 * t375;
t118 = -qJD(2) * t401 + t204 * t375;
t113 = -mrSges(6,1) * t190 + mrSges(6,2) * t415;
t103 = mrSges(4,1) * t185 + mrSges(4,2) * t184;
t102 = mrSges(5,1) * t185 - mrSges(5,3) * t184;
t94 = t147 - t571;
t88 = mrSges(7,1) * t314 - mrSges(7,3) * t112;
t87 = -mrSges(7,2) * t314 + mrSges(7,3) * t456;
t60 = -pkin(5) * t118 + t117;
t57 = -mrSges(6,2) * t273 + mrSges(6,3) * t65;
t56 = mrSges(6,1) * t273 - mrSges(6,3) * t64;
t51 = -mrSges(7,1) * t456 + mrSges(7,2) * t112;
t40 = qJD(6) * t419 + t118 * t377 - t119 * t372;
t39 = qJD(6) * t153 + t118 * t372 + t119 * t377;
t36 = -mrSges(6,1) * t65 + mrSges(6,2) * t64;
t29 = -pkin(5) * t65 + t53;
t24 = pkin(10) * t118 + t34;
t23 = -pkin(5) * t504 - pkin(10) * t119 + t35;
t18 = t377 * t45 - t546;
t17 = -t372 * t45 - t541;
t16 = -mrSges(7,2) * t258 + mrSges(7,3) * t22;
t15 = mrSges(7,1) * t258 - mrSges(7,3) * t21;
t8 = -mrSges(7,1) * t22 + mrSges(7,2) * t21;
t5 = -qJD(6) * t38 + t23 * t377 - t24 * t372;
t4 = qJD(6) * t37 + t23 * t372 + t24 * t377;
t1 = [t522 * t703 + t561 * t704 + t432 * t690 + m(7) * (t13 * t5 + t14 * t4 + t146 * t29 + t2 * t38 + t3 * t37 + t60 * t86) + m(6) * (t100 * t12 + t101 * t11 + t117 * t131 + t214 * t53 + t34 * t59 + t35 * t58) + m(5) * (t122 * t160 + t130 * t158 + t137 * t166 + t217 * t66 + t218 * t69 + t241 * t68) + (Ifges(7,4) * t39 + Ifges(7,2) * t40) * t602 + (t534 * t618 + t416 * mrSges(6,1) + t418 * mrSges(7,1) + t651 * mrSges(6,2) + t650 * mrSges(7,2) + t642 * (-t262 * pkin(3) - qJ(4) * t261 + t367) + t687 * t381 + (-m(3) - m(4)) * t367 + t626 * t262 + t688 * t261 + (m(3) * pkin(1) + t691 * t484 + (-m(4) - m(5)) * t410 + (m(6) * t604 - m(7) * (-pkin(8) - t382) - t686) * t375 + t649 + t676) * t376) * g(1) + (t158 * t399 - t160 * t398 + t522 * t69) * mrSges(5,2) + (Ifges(7,1) * t39 + Ifges(7,4) * t40) * t600 + (-m(3) * t507 - m(4) * t452 - m(7) * t421 - t174 * mrSges(6,1) - t156 * mrSges(7,1) - t173 * mrSges(6,2) - t155 * mrSges(7,2) + t692 * (qJ(4) * t263 + t421) + (-m(7) * t520 - t676) * t381 + t687 * t376 - t626 * t264 + t633 * t263 + t706 * t521) * g(2) + qJD(2) ^ 2 * t427 / 0.2e1 + (-t196 * t399 - t197 * t398 - t522 * t85) * mrSges(4,3) + t443 * t536 + (-t11 * t248 + t118 * t59 - t119 * t58 - t12 * t249) * mrSges(6,3) + (Ifges(6,5) * t119 + Ifges(6,6) * t118) * t582 + (Ifges(6,5) * t249 - Ifges(6,6) * t248) * t591 + (t196 * mrSges(4,1) - t158 * mrSges(5,1) - t197 * mrSges(4,2) + t160 * mrSges(5,3) - t695 - t721) * t504 + t318 * (mrSges(4,1) * t398 + mrSges(4,2) * t399) + t166 * (mrSges(5,1) * t398 - mrSges(5,3) * t399) + (Ifges(6,1) * t249 - Ifges(6,4) * t248) * t610 + t37 * t15 + t38 * t16 + (Ifges(7,5) * t39 + Ifges(7,6) * t40) * t584 + t248 * t693 + (Ifges(6,1) * t119 + Ifges(6,4) * t118) * t593 + (Ifges(6,4) * t119 + Ifges(6,2) * t118) * t595 + (Ifges(6,4) * t249 - Ifges(6,2) * t248) * t609 + (-t13 * t39 + t14 * t40 + t153 * t2 + t3 * t419) * mrSges(7,3) + (-Ifges(7,4) * t419 + Ifges(7,2) * t153) * t616 + (-Ifges(7,1) * t419 + Ifges(7,4) * t153) * t617 + (-Ifges(7,5) * t419 + Ifges(7,6) * t153) * t592 + t29 * (-mrSges(7,1) * t153 - mrSges(7,2) * t419) - t419 * t619 + t477 * t503 + Ifges(2,3) * qJDD(1) + t699 * (t374 * t464 + t379 * t467) + (qJD(2) * t663 - t500 * t659) * t723 - t406 * t495 - t316 * t491 - pkin(1) * (-mrSges(3,1) * t308 + mrSges(3,2) * t309) + t256 * t276 + t241 * t102 + t231 * t126 + t232 * t128 + t214 * t36 + t217 * t129 + t218 * t127 + (m(4) * t256 * pkin(7) + Ifges(3,1) * t309 + Ifges(3,4) * t690 + Ifges(3,5) * qJDD(2) + t167 * t463 + t424 * t597 + t431 * t598 + t68 * t438 - t460 * t679 + t590 * t660 + t599 * t657) * t375 + t138 * t221 + t139 * t222 + t130 * t223 + t122 * t224 + m(4) * (t138 * t197 + t139 * t196 + t231 * t85 + t232 * t84 + t318 * t490) + t137 * t209 + t60 * t51 + (-t66 * mrSges(5,2) - t84 * mrSges(4,3) - t637) * t525 + t379 * t170 * t464 + (Ifges(3,4) * t704 + Ifges(3,2) * t690 + t492 / 0.2e1 + t493 / 0.2e1 + Ifges(3,6) * qJDD(2) - t654 / 0.2e1 - t717 * t599 + pkin(7) * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t308) + Ifges(6,3) * t591 + Ifges(7,3) * t592 - Ifges(5,6) * t597 - Ifges(4,6) * t598 + Ifges(6,6) * t609 + Ifges(6,5) * t610 + Ifges(7,6) * t616 + Ifges(7,5) * t617 + t715 * t590 + t625 + t561 * t460 + t630 + t631) * t380 + t86 * (-mrSges(7,1) * t40 + mrSges(7,2) * t39) + t4 * t87 + t5 * t88 + (t309 * t570 + t652) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(7) * t652) + t653 * t467 + t100 * t56 + t101 * t57 + (qJD(2) * t662 - t500 * t658) * t586 + t664 * t490 + t403 * t460 + (qJD(2) * t390 - t423 * t500) * t588 + (qJD(2) * t391 - t430 * t500) * t589 + t119 * t605 + t118 * t607 + t39 * t611 + t40 * t613 + t249 * t615 + t153 * t620 + t117 * t113 + t53 * t459 + (-qJDD(2) * mrSges(3,1) + t103) * t570 + t131 * (-mrSges(6,1) * t118 + mrSges(6,2) * t119) + t34 * t144 + t35 * t145 + t146 * t8; (t349 + t653) * t462 + t374 * t703 + (-t489 - t206) * t224 + (Ifges(7,1) * t136 + Ifges(7,4) * t135) * t601 + (-Ifges(6,5) * t412 - Ifges(6,6) * t411) * t591 + (-Ifges(6,4) * t412 - Ifges(6,2) * t411) * t609 + (-Ifges(6,1) * t412 - Ifges(6,4) * t411) * t610 - t412 * t615 + (-t158 * t207 - t160 * t206 + t665 * t166 + t639 * t68) * m(5) + t639 * t102 + (Ifges(7,1) * t600 + Ifges(7,4) * t602 + Ifges(7,5) * t584 + t611 + t634) * (qJD(6) * t193 + t204 * t377 - t205 * t372) + (Ifges(7,4) * t600 + Ifges(7,2) * t602 + Ifges(7,6) * t584 + t613 + t635) * (-qJD(6) * t194 - t204 * t372 - t205 * t377) + (Ifges(6,4) * t234 - Ifges(6,2) * t233) * t596 + (-t11 * t411 + t12 * t412 + t233 * t59 + t234 * t58) * mrSges(6,3) + (Ifges(6,1) * t234 - Ifges(6,4) * t233) * t594 + (Ifges(6,5) * t234 - Ifges(6,6) * t233) * t583 + (t376 * t627 + t519 * t632) * g(2) + (t381 * t627 + t517 * t632) * g(1) - t68 * t439 - t256 * t441 + (-m(6) * (t453 - t567) - m(5) * t453 - m(4) * t508 - m(7) * (t453 + t520) - t443 + t686 * t375 + (-m(6) * t363 - m(7) * t407 - t705) * t380 - t649) * g(3) + (Ifges(7,4) * t136 + Ifges(7,2) * t135) * t603 + (-t488 - t219) * t222 - (Ifges(6,4) * t593 + Ifges(6,2) * t595 + Ifges(6,6) * t582 + t574 + t607) * t205 + (t565 * t660 + t566 * t663) * t333 + (t710 * mrSges(6,1) + mrSges(6,2) * t724) * t131 + (Ifges(6,1) * t593 + Ifges(6,4) * t595 + Ifges(6,5) * t582 - t575 + t605) * t204 + (-t489 - t220) * t221 + (t11 * t216 + t12 * t215 + t131 * t666 - t53 * t696 + t58 * t673 + t59 * t672) * m(6) - t696 * t36 + (t565 * t657 + t566 * t662) * t293 + t411 * t693 + (t13 * t136 - t135 * t14 + t193 * t2 - t194 * t3) * mrSges(7,3) + (-t207 + t488) * t223 + (-t197 * (-mrSges(4,2) * t375 - mrSges(4,3) * t523) - t160 * (-mrSges(5,2) * t523 + mrSges(5,3) * t375) - t196 * (mrSges(4,1) * t375 - mrSges(4,3) * t518) - t158 * (-mrSges(5,1) * t375 + mrSges(5,2) * t518) + (t391 / 0.2e1 - t390 / 0.2e1) * t292 + (t406 - t403 / 0.2e1) * qJD(1)) * qJD(1) + t699 * (t379 * t462 + t463) + (-Ifges(6,5) * t594 - Ifges(7,5) * t601 - Ifges(3,2) * t462 - Ifges(6,6) * t596 - Ifges(7,6) * t603 - Ifges(6,3) * t583 - Ifges(7,3) * t585 + t695) * t505 + (-pkin(2) * t256 - t196 * t219 - t197 * t220 - t318 * t352) * m(4) + ((t129 + t128) * t379 + (-t126 + t127) * t374 + ((-t196 * t379 - t197 * t374) * qJD(3) + t648) * m(4) + ((t158 * t379 - t160 * t374) * qJD(3) + t647) * m(5)) * pkin(8) + t53 * t694 + (Ifges(7,5) * t136 + Ifges(7,6) * t135) * t585 - t427 * t495 / 0.2e1 + Ifges(3,3) * qJDD(2) + (-t431 / 0.2e1 + t424 / 0.2e1) * t502 + Ifges(3,6) * t308 + Ifges(3,5) * t309 - t288 * mrSges(3,2) - t289 * mrSges(3,1) + t170 * t481 / 0.2e1 + t215 * t56 + t216 * t57 + t211 * t8 + t29 * (-mrSges(7,1) * t193 + mrSges(7,2) * t194) + t637 * t379 + t316 * t351 + t645 * t499 + (t167 / 0.2e1 + t646) * t501 + t647 * mrSges(5,2) + t648 * mrSges(4,3) + t89 * t15 + t90 * t16 + t658 * t599 + t659 * t590 + (t477 - t661) * qJD(3) + t661 * t357 - pkin(2) * t103 - t664 * t352 + t665 * t209 + t666 * t113 + t670 * t51 + (Ifges(7,5) * t194 + Ifges(7,6) * t193) * t592 + t423 * t597 + t430 * t598 + t234 * t606 - t233 * t608 + t136 * t612 + t135 * t614 + (Ifges(7,4) * t194 + Ifges(7,2) * t193) * t616 + (Ifges(7,1) * t194 + Ifges(7,4) * t193) * t617 + t194 * t619 + t193 * t620 - t86 * (-mrSges(7,1) * t135 + mrSges(7,2) * t136) + t672 * t144 + t673 * t145 + t681 * t87 + t682 * t88 + (t13 * t682 + t14 * t681 + t2 * t90 + t211 * t29 + t3 * t89 + t670 * t86) * m(7); (mrSges(6,1) * t131 + Ifges(6,4) * t594 + Ifges(6,2) * t596 + Ifges(6,6) * t583 - t574 + t608) * t415 + (-t318 * mrSges(4,1) - t166 * mrSges(5,1) + Ifges(5,3) * t589 - t581 * t716 - t646) * t293 - t625 - t707 + t654 + (-t221 - t224) * t196 - (-mrSges(6,2) * t131 + Ifges(6,1) * t594 + Ifges(6,4) * t596 + Ifges(6,5) * t583 + t575 + t606) * t190 + (t222 - t223) * t197 + (m(7) * t253 + t692 * (qJ(4) * t264 - t253) + t633 * t264 + t626 * t263 + t698) * g(1) + (-Ifges(4,2) * t293 - t283 + t699) * t588 + (-(pkin(5) * t527 + (-pkin(3) - t346) * t374) * t576 + t276 + (-t540 - (-m(5) * pkin(3) - mrSges(5,1)) * t374 - m(6) * t487) * t375 + t642 * t330 - t697) * g(3) + t310 * t56 + t311 * t57 + qJD(4) * t224 + t200 * t15 + t201 * t16 - t208 * t209 - t548 * t589 + (-t292 * t718 + t167 + t282 - t547) * t720 + (-pkin(3) * t69 + qJ(4) * t66 - t158 * t197 + t160 * t643 - t166 * t208) * m(5) + (t318 * mrSges(4,2) - t166 * mrSges(5,3) - t581 * t717 + t645) * t292 + (m(7) * t251 + t651 * mrSges(6,1) + t692 * (qJ(4) * t262 - t251) + t633 * t262 + t626 * t261 - t722) * g(2) - t94 * t51 + t170 * t586 - pkin(3) * t127 + qJ(4) * t129 + t674 * t145 + t675 * t144 + (t11 * t311 + t12 * t310 - t131 * t147 + t58 * t674 + t59 * t675) * m(6) - t147 * t113 + t683 * t87 + t684 * t88 + (t13 * t684 + t14 * t683 + t2 * t201 + t200 * t3 - t86 * t94) * m(7); t294 * t15 + t296 * t16 + t373 * t57 + t378 * t56 - t668 * t88 + t667 * t87 + t420 * qJD(5) + (-t224 - t420) * t333 + (t209 - t51 - t113) * t293 + t127 + (-t13 * t668 + t14 * t667 + t2 * t296 - t293 * t86 + t294 * t3 + t396) * m(7) + (t11 * t373 + t12 * t378 - t131 * t293 + t396 + t327 * (-t373 * t58 + t378 * t59)) * m(6) + (-t160 * t333 + t166 * t293 + t396 + t69) * m(5); -m(7) * (t13 * t17 + t14 * t18 + t571 * t86) - t51 * t571 - t131 * (mrSges(6,1) * t415 + mrSges(6,2) * t190) - t18 * t87 - t17 * t88 + (Ifges(6,5) * t190 - Ifges(6,6) * t415) * t583 + t92 * t593 + (Ifges(6,1) * t190 - t558) * t594 + (t2 * t372 + t3 * t377 + (-t13 * t372 + t14 * t377) * qJD(6)) * t618 + (-Ifges(6,2) * t415 + t187 + t93) * t596 + (t563 + t145) * t59 + (t564 - t144) * t58 + t697 * g(3) + ((-mrSges(6,1) - t618) * t651 + t722) * g(2) + (-t173 * t618 - t698) * g(1) + ((-t372 * t88 + t377 * t87) * qJD(6) + g(3) * t412 * t576 + t15 * t377 + t16 * t372) * pkin(5) + t707; -t86 * (mrSges(7,1) * t112 + mrSges(7,2) * t456) + (Ifges(7,1) * t456 - t557) * t601 + t43 * t600 + (Ifges(7,5) * t456 - Ifges(7,6) * t112) * t585 - t13 * t87 + t14 * t88 - g(1) * t514 - g(2) * t515 - g(3) * t513 + (t112 * t14 + t13 * t456) * mrSges(7,3) + t408 + (-Ifges(7,2) * t112 + t104 + t44) * t603;];
tau  = t1;
