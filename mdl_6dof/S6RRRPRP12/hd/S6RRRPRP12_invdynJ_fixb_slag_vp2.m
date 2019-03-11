% Calculate vector of inverse dynamics joint torques for
% S6RRRPRP12
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
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RRRPRP12_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP12_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:51:43
% EndTime: 2019-03-09 17:53:03
% DurationCPUTime: 49.92s
% Computational Cost: add. (13884->987), mult. (33208->1250), div. (0->0), fcn. (25594->10), ass. (0->422)
t351 = cos(qJ(2));
t344 = sin(pkin(6));
t510 = qJD(1) * t344;
t481 = t351 * t510;
t311 = -qJD(3) + t481;
t347 = sin(qJ(2));
t535 = cos(pkin(6));
t452 = t535 * qJD(1);
t438 = pkin(1) * t452;
t264 = pkin(8) * t481 + t347 * t438;
t346 = sin(qJ(3));
t441 = t346 * t481;
t505 = qJD(3) * t346;
t702 = -qJD(4) * t346 - t264 + (-t441 + t505) * pkin(3);
t670 = Ifges(6,1) + Ifges(7,1);
t701 = Ifges(4,4) + Ifges(5,6);
t669 = Ifges(7,4) + Ifges(6,5);
t509 = qJD(1) * t347;
t482 = t344 * t509;
t261 = -pkin(8) * t482 + t351 * t438;
t381 = t344 * (pkin(2) * t347 - pkin(9) * t351);
t262 = qJD(1) * t381;
t350 = cos(qJ(3));
t163 = -t346 * t261 + t262 * t350;
t504 = qJD(3) * t350;
t515 = t350 * t351;
t598 = pkin(4) + pkin(9);
t599 = pkin(3) + pkin(10);
t700 = -(pkin(4) * t515 - t347 * t599) * t510 + t163 + t598 * t504;
t699 = -t702 + t311 * (pkin(10) * t346 - qJ(4) * t350);
t390 = t452 + qJD(2);
t373 = qJD(3) * t390;
t498 = qJD(1) * qJD(2);
t377 = qJDD(1) * t347 + t351 * t498;
t450 = t535 * qJDD(1);
t385 = t450 + qJDD(2);
t480 = t350 * t509;
t137 = t346 * t373 - t350 * t385 + (qJD(3) * t480 + t346 * t377) * t344;
t268 = (-qJDD(1) * t351 + t347 * t498) * t344;
t256 = qJDD(3) + t268;
t345 = sin(qJ(5));
t349 = cos(qJ(5));
t228 = t346 * t482 - t350 * t390;
t393 = t349 * t228 + t311 * t345;
t66 = qJD(5) * t393 + t137 * t345 + t256 * t349;
t603 = t66 / 0.2e1;
t173 = t228 * t345 - t311 * t349;
t67 = qJD(5) * t173 - t349 * t137 + t256 * t345;
t601 = t67 / 0.2e1;
t365 = t344 * t377;
t526 = t344 * t347;
t492 = t346 * t526;
t439 = qJD(3) * t492;
t136 = qJD(1) * t439 - t346 * t385 + (-t365 - t373) * t350;
t133 = qJDD(5) - t136;
t595 = t133 / 0.2e1;
t594 = -t136 / 0.2e1;
t592 = -t137 / 0.2e1;
t577 = t256 / 0.2e1;
t668 = Ifges(4,5) - Ifges(5,4);
t667 = Ifges(7,5) - Ifges(6,4);
t666 = Ifges(4,6) - Ifges(5,5);
t665 = Ifges(6,6) - Ifges(7,6);
t664 = Ifges(4,3) + Ifges(5,1);
t663 = Ifges(6,3) + Ifges(7,2);
t440 = t350 * t481;
t698 = t440 - t504;
t593 = t136 / 0.2e1;
t602 = -t67 / 0.2e1;
t229 = t344 * t480 + t346 * t390;
t223 = qJD(5) + t229;
t691 = -pkin(8) * t344 * t498 + pkin(1) * t450;
t497 = qJDD(1) * t344;
t692 = pkin(8) * t497 + qJD(2) * t438;
t174 = t347 * t691 + t351 * t692;
t155 = pkin(9) * t385 + t174;
t162 = t268 * pkin(2) + (-qJDD(1) * pkin(1) - pkin(9) * t377) * t344;
t211 = pkin(9) * t390 + t264;
t219 = (-pkin(2) * t351 - pkin(9) * t347 - pkin(1)) * t510;
t41 = -t346 * t155 + t162 * t350 - t211 * t504 - t219 * t505;
t374 = qJDD(4) - t41;
t18 = -pkin(4) * t136 - t256 * t599 + t374;
t175 = -t347 * t692 + t351 * t691;
t156 = -pkin(2) * t385 - t175;
t353 = t136 * qJ(4) - t229 * qJD(4) + t156;
t23 = t137 * t599 + t353;
t501 = qJD(5) * t349;
t502 = qJD(5) * t345;
t134 = t211 * t346 - t350 * t219;
t383 = pkin(4) * t229 + t134;
t689 = qJD(4) + t383;
t79 = t311 * t599 + t689;
t210 = -pkin(2) * t390 - t261;
t354 = -t229 * qJ(4) + t210;
t82 = t228 * t599 + t354;
t3 = t345 * t18 + t349 * t23 + t79 * t501 - t502 * t82;
t1 = qJ(6) * t133 + qJD(6) * t223 + t3;
t27 = t345 * t79 + t349 * t82;
t4 = -qJD(5) * t27 + t18 * t349 - t23 * t345;
t2 = -pkin(5) * t133 + qJDD(6) - t4;
t612 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t662 = t133 * t663 + t66 * t669 - t665 * t67;
t694 = -t256 / 0.2e1;
t697 = t612 + 0.2e1 * Ifges(4,1) * t594 + Ifges(5,4) * t694 + Ifges(6,6) * t602 + Ifges(7,6) * t601 + t595 * t663 + t603 * t669 + (-t593 + t594) * Ifges(5,2) + t662 / 0.2e1 + t701 * t592 + (t668 + Ifges(4,5)) * t577;
t108 = t228 * pkin(3) + t354;
t135 = t350 * t211 + t346 * t219;
t297 = t311 * qJ(4);
t114 = t297 - t135;
t630 = t114 * mrSges(5,1) - t135 * mrSges(4,3);
t696 = t210 * mrSges(4,1) - t108 * mrSges(5,2) + t630;
t94 = -pkin(4) * t228 + t135;
t83 = -t297 + t94;
t39 = -pkin(5) * t393 - qJ(6) * t173 + t83;
t424 = mrSges(7,1) * t349 + mrSges(7,3) * t345;
t426 = mrSges(6,1) * t349 - mrSges(6,2) * t345;
t556 = Ifges(6,4) * t173;
t76 = Ifges(6,2) * t393 + t223 * Ifges(6,6) + t556;
t538 = t349 * t76;
t170 = Ifges(7,5) * t173;
t73 = t223 * Ifges(7,6) - Ifges(7,3) * t393 + t170;
t539 = t349 * t73;
t583 = -t223 / 0.2e1;
t641 = t345 * t669 + t349 * t665;
t695 = t583 * t641 + t39 * t424 + t83 * t426 - t538 / 0.2e1 + t539 / 0.2e1;
t591 = t137 / 0.2e1;
t587 = t173 / 0.2e1;
t661 = t133 * t669 + t66 * t670 + t667 * t67;
t693 = -t661 / 0.2e1;
t171 = Ifges(6,4) * t393;
t553 = Ifges(7,5) * t393;
t656 = t173 * t670 + t669 * t223 + t171 - t553;
t533 = qJ(4) * t346;
t460 = -pkin(2) - t533;
t289 = -t350 * t599 + t460;
t315 = t598 * t346;
t646 = t349 * t289 + t345 * t315;
t652 = -qJD(5) * t646 + t345 * t699 + t349 * t700;
t651 = -t289 * t502 + t315 * t501 + t345 * t700 - t349 * t699;
t626 = -qJD(4) - t134;
t112 = pkin(3) * t311 - t626;
t690 = t112 * mrSges(5,1) + t134 * mrSges(4,3);
t627 = -t66 * Ifges(6,4) / 0.2e1 - t133 * Ifges(6,6) / 0.2e1 + Ifges(7,5) * t603 + Ifges(7,6) * t595 + (Ifges(6,2) + Ifges(7,3)) * t601;
t348 = sin(qJ(1));
t571 = cos(qJ(1));
t429 = t535 * t571;
t282 = t347 * t429 + t348 * t351;
t483 = t344 * t571;
t198 = t282 * t346 + t350 * t483;
t281 = t347 * t348 - t351 * t429;
t687 = t198 * t345 + t281 * t349;
t686 = t198 * t349 - t281 * t345;
t685 = Ifges(4,6) * t694 + 0.2e1 * Ifges(5,3) * t591 + t701 * t593 + (t591 - t592) * Ifges(4,2) + (-t666 + Ifges(5,5)) * t577;
t26 = -t345 * t82 + t349 * t79;
t650 = qJD(6) - t26;
t24 = -pkin(5) * t223 + t650;
t25 = qJ(6) * t223 + t27;
t684 = t26 * mrSges(6,1) - t24 * mrSges(7,1) - t27 * mrSges(6,2) + t25 * mrSges(7,3);
t220 = Ifges(5,6) * t228;
t125 = -t311 * Ifges(5,4) - t229 * Ifges(5,2) + t220;
t221 = Ifges(4,4) * t228;
t622 = t229 * Ifges(4,1) - t311 * Ifges(4,5) + t173 * t669 + t223 * t663 + t393 * t665 - t221;
t683 = -t125 / 0.2e1 + t622 / 0.2e1 + t690;
t682 = -m(7) - m(6);
t681 = -t268 / 0.2e1;
t680 = t365 / 0.2e1;
t679 = t385 / 0.2e1;
t674 = -mrSges(5,1) - mrSges(4,3);
t673 = mrSges(4,2) - mrSges(5,3);
t672 = -mrSges(5,2) + mrSges(4,1);
t671 = mrSges(6,3) + mrSges(7,2);
t30 = -mrSges(7,2) * t67 + mrSges(7,3) * t133;
t33 = -mrSges(6,2) * t133 - mrSges(6,3) * t67;
t660 = t30 + t33;
t31 = mrSges(6,1) * t133 - mrSges(6,3) * t66;
t32 = -t133 * mrSges(7,1) + t66 * mrSges(7,2);
t659 = t31 - t32;
t658 = -qJ(6) * t698 + qJD(6) * t346 + t651;
t657 = pkin(5) * t698 - t652;
t164 = t350 * t261 + t346 * t262;
t143 = -qJ(4) * t482 - t164;
t119 = -pkin(4) * t441 - t143;
t217 = t345 * t482 - t349 * t441;
t518 = t346 * t351;
t233 = (t345 * t518 + t347 * t349) * t344;
t218 = qJD(1) * t233;
t402 = pkin(5) * t349 + qJ(6) * t345;
t384 = -pkin(4) - t402;
t401 = -pkin(5) * t345 + qJ(6) * t349;
t655 = -pkin(5) * t217 + qJ(6) * t218 - t119 + (qJD(5) * t401 + qJD(6) * t345) * t350 + (-pkin(9) + t384) * t505;
t179 = mrSges(5,1) * t228 + mrSges(5,3) * t311;
t92 = -mrSges(6,1) * t393 + mrSges(6,2) * t173;
t654 = -t179 + t92;
t653 = qJD(5) * t402 - qJD(6) * t349 - t229 * t384 - t626;
t649 = -t228 * t666 + t229 * t668 - t311 * t664;
t500 = qJD(5) * t350;
t477 = t345 * t500;
t648 = t349 * t505 + t217 + t477;
t647 = -t345 * t505 + t349 * t500 + t218;
t645 = qJ(4) * t698 + t702;
t644 = -t598 * t505 - t119;
t453 = t351 * t535;
t287 = pkin(1) * t453 - pkin(8) * t526;
t302 = qJ(4) - t401;
t423 = t345 * mrSges(7,1) - t349 * mrSges(7,3);
t425 = mrSges(6,1) * t345 + mrSges(6,2) * t349;
t643 = -m(6) * qJ(4) - m(7) * t302 - t423 - t425;
t642 = -t108 * (-mrSges(5,2) * t346 - mrSges(5,3) * t350) - t210 * (mrSges(4,1) * t346 + mrSges(4,2) * t350);
t640 = -t346 * t666 + t350 * t668;
t639 = -t345 * t665 + t349 * t669;
t551 = Ifges(7,5) * t349;
t554 = Ifges(6,4) * t349;
t638 = -t345 * t670 + t551 - t554;
t552 = Ifges(7,5) * t345;
t555 = Ifges(6,4) * t345;
t637 = t349 * t670 + t552 - t555;
t421 = mrSges(5,2) * t350 - mrSges(5,3) * t346;
t427 = mrSges(4,1) * t350 - mrSges(4,2) * t346;
t636 = t421 - t427;
t635 = -t136 * t668 - t137 * t666 + t256 * t664;
t634 = -t229 * t349 - t501;
t532 = t229 * t345;
t633 = t502 + t532;
t40 = t350 * t155 + t346 * t162 - t211 * t505 + t219 * t504;
t632 = -t346 * t41 + t350 * t40;
t28 = -t256 * qJ(4) + t311 * qJD(4) - t40;
t34 = -pkin(3) * t256 + t374;
t631 = -t28 * t350 + t34 * t346;
t629 = -t3 * t345 - t349 * t4;
t628 = -t1 * t345 + t2 * t349;
t560 = Ifges(3,4) * t347;
t608 = t344 ^ 2;
t625 = (t347 * (Ifges(3,1) * t351 - t560) / 0.2e1 - pkin(1) * (mrSges(3,1) * t347 + mrSges(3,2) * t351)) * t608;
t624 = mrSges(3,2) + t674;
t621 = mrSges(3,1) - t636;
t619 = -t671 - t672;
t508 = qJD(2) * t351;
t478 = t344 * t508;
t195 = -t439 + (qJD(3) * t535 + t478) * t350;
t479 = qJD(2) * t526;
t454 = t347 * t535;
t523 = t344 * t351;
t288 = pkin(1) * t454 + pkin(8) * t523;
t252 = pkin(9) * t535 + t288;
t512 = pkin(2) * t523 + pkin(9) * t526;
t570 = pkin(1) * t344;
t253 = -t512 - t570;
t263 = qJD(2) * t381;
t265 = t287 * qJD(2);
t85 = -t252 * t504 - t253 * t505 + t263 * t350 - t346 * t265;
t49 = pkin(4) * t195 - t479 * t599 - t85;
t160 = -t346 * t252 + t253 * t350;
t141 = pkin(3) * t523 - t160;
t524 = t344 * t350;
t280 = t346 * t535 + t347 * t524;
t100 = pkin(4) * t280 + pkin(10) * t523 + t141;
t251 = -pkin(2) * t535 - t287;
t279 = -t350 * t535 + t492;
t269 = t279 * pkin(3);
t451 = t280 * qJ(4) - t269;
t139 = t251 - t451;
t566 = t279 * pkin(10);
t109 = t139 + t566;
t536 = t345 * t100 + t349 * t109;
t194 = qJD(3) * t280 + t346 * t478;
t266 = t288 * qJD(2);
t355 = -t195 * qJ(4) - t280 * qJD(4) + t266;
t60 = t194 * t599 + t355;
t9 = -qJD(5) * t536 - t345 * t60 + t349 * t49;
t448 = m(7) * pkin(5) + mrSges(6,1) + mrSges(7,1);
t443 = m(7) * qJ(6) - mrSges(6,2) + mrSges(7,3);
t446 = mrSges(3,3) * t482;
t618 = -m(4) * t210 + mrSges(3,1) * t390 - mrSges(4,1) * t228 - mrSges(4,2) * t229 - t446;
t616 = -m(5) * qJ(4) + t643 + t673;
t284 = -t348 * t454 + t351 * t571;
t202 = t284 * t346 - t348 * t524;
t615 = g(1) * t202 + g(2) * t198 + g(3) * t279;
t611 = -t41 * mrSges(4,1) + t40 * mrSges(4,2) - t34 * mrSges(5,2) + t28 * mrSges(5,3);
t610 = t210 * mrSges(4,2) - t108 * mrSges(5,3) + t684;
t600 = t76 / 0.2e1;
t590 = t393 / 0.2e1;
t589 = -t393 / 0.2e1;
t588 = -t173 / 0.2e1;
t582 = t223 / 0.2e1;
t581 = -t228 / 0.2e1;
t580 = t228 / 0.2e1;
t579 = -t229 / 0.2e1;
t578 = t229 / 0.2e1;
t575 = -t311 / 0.2e1;
t574 = t311 / 0.2e1;
t341 = t350 * pkin(9);
t562 = mrSges(6,3) * t393;
t561 = mrSges(6,3) * t173;
t559 = Ifges(3,4) * t351;
t558 = Ifges(4,4) * t346;
t557 = Ifges(4,4) * t350;
t550 = Ifges(5,6) * t346;
t549 = Ifges(5,6) * t350;
t544 = t229 * Ifges(4,4);
t543 = t229 * Ifges(5,6);
t534 = qJ(4) * t228;
t113 = t229 * t599 + t534;
t43 = t349 * t113 + t345 * t94;
t528 = t281 * t350;
t283 = t347 * t571 + t348 * t453;
t527 = t283 * t350;
t525 = t344 * t348;
t522 = t345 * t346;
t521 = t345 * t350;
t519 = t346 * t349;
t517 = t349 * t350;
t115 = mrSges(7,2) * t393 + mrSges(7,3) * t223;
t116 = -mrSges(6,2) * t223 + t562;
t514 = t115 + t116;
t117 = mrSges(6,1) * t223 - t561;
t118 = -mrSges(7,1) * t223 + mrSges(7,2) * t173;
t513 = -t117 + t118;
t161 = t350 * t252 + t346 * t253;
t316 = t350 * pkin(4) + t341;
t511 = t571 * pkin(1) + pkin(8) * t525;
t499 = qJD(5) * t599;
t496 = pkin(9) * t505;
t495 = pkin(9) * t504;
t493 = qJ(4) * t523;
t491 = t349 * t523;
t490 = t344 * t515;
t485 = Ifges(3,5) * t365 - Ifges(3,6) * t268 + Ifges(3,3) * t385;
t476 = t345 * t499;
t475 = t349 * t499;
t473 = t526 / 0.2e1;
t470 = -t510 / 0.2e1;
t469 = t510 / 0.2e1;
t465 = t505 / 0.2e1;
t462 = -t500 / 0.2e1;
t461 = -pkin(1) * t348 + pkin(8) * t483;
t459 = -t281 * pkin(2) + pkin(9) * t282;
t458 = -t283 * pkin(2) + pkin(9) * t284;
t99 = -t136 * mrSges(5,1) + t256 * mrSges(5,2);
t199 = t282 * t350 - t346 * t483;
t445 = mrSges(3,3) * t481;
t444 = pkin(3) * t490 + t346 * t493 + t512;
t435 = t351 * t470;
t434 = t351 * t469;
t428 = mrSges(4,1) * t279 + mrSges(4,2) * t280;
t422 = -t279 * mrSges(5,2) - t280 * mrSges(5,3);
t420 = Ifges(4,1) * t350 - t558;
t415 = -Ifges(4,2) * t346 + t557;
t413 = -Ifges(6,2) * t345 + t554;
t412 = Ifges(6,2) * t349 + t555;
t406 = Ifges(7,3) * t345 + t551;
t405 = -Ifges(7,3) * t349 + t552;
t404 = -Ifges(5,2) * t350 + t550;
t403 = Ifges(5,3) * t346 - t549;
t398 = t24 * t345 + t25 * t349;
t397 = t26 * t345 - t27 * t349;
t42 = -t113 * t345 + t349 * t94;
t44 = t100 * t349 - t109 * t345;
t187 = -t289 * t345 + t315 * t349;
t389 = -pkin(3) * t528 - t281 * t533 + t459;
t388 = -pkin(3) * t527 - t283 * t533 + t458;
t387 = t284 * pkin(2) + pkin(9) * t283 + t511;
t140 = t493 - t161;
t196 = t279 * t349 + t345 * t523;
t8 = t100 * t501 - t109 * t502 + t345 * t49 + t349 * t60;
t378 = -t282 * pkin(2) - pkin(9) * t281 + t461;
t84 = -t252 * t505 + t253 * t504 + t346 * t263 + t350 * t265;
t110 = -pkin(4) * t279 - t140;
t203 = t284 * t350 + t346 * t525;
t364 = t203 * pkin(3) + qJ(4) * t202 + t387;
t363 = -pkin(3) * t199 - qJ(4) * t198 + t378;
t362 = Ifges(3,6) * t535 + (Ifges(3,2) * t351 + t560) * t344;
t21 = -pkin(4) * t137 - t28;
t72 = -qJ(4) * t479 + qJD(4) * t523 - t84;
t361 = qJD(5) * t398 - t628;
t360 = -qJD(5) * t397 - t629;
t358 = t344 * t390 * (Ifges(3,5) * t351 - Ifges(3,6) * t347);
t50 = -pkin(4) * t194 - t72;
t326 = Ifges(3,4) * t481;
t304 = -pkin(3) * t350 + t460;
t285 = (-mrSges(3,1) * t351 + mrSges(3,2) * t347) * t344;
t260 = -mrSges(3,2) * t390 + t445;
t222 = t350 * t402 + t316;
t207 = Ifges(3,1) * t482 + Ifges(3,5) * t390 + t326;
t206 = Ifges(3,6) * qJD(2) + qJD(1) * t362;
t197 = -t279 * t345 + t491;
t191 = t202 * pkin(3);
t189 = t198 * pkin(3);
t181 = -pkin(5) * t346 - t187;
t180 = mrSges(5,1) * t229 - mrSges(5,2) * t311;
t178 = -mrSges(4,1) * t311 - mrSges(4,3) * t229;
t177 = mrSges(4,2) * t311 - mrSges(4,3) * t228;
t176 = qJ(6) * t346 + t646;
t159 = -mrSges(5,2) * t228 - mrSges(5,3) * t229;
t157 = pkin(3) * t229 + t534;
t151 = t202 * t345 + t283 * t349;
t150 = -t202 * t349 + t283 * t345;
t145 = -pkin(3) * t482 - t163;
t126 = -t228 * Ifges(4,2) - t311 * Ifges(4,6) + t544;
t123 = -t311 * Ifges(5,5) + t228 * Ifges(5,3) - t543;
t106 = qJD(5) * t196 + t194 * t345 + t349 * t479;
t105 = -t194 * t349 - qJD(5) * t491 + (qJD(5) * t279 + t479) * t345;
t98 = mrSges(5,1) * t137 - mrSges(5,3) * t256;
t97 = -mrSges(4,2) * t256 - mrSges(4,3) * t137;
t96 = mrSges(4,1) * t256 + mrSges(4,3) * t136;
t91 = -mrSges(7,1) * t393 - mrSges(7,3) * t173;
t90 = pkin(5) * t173 - qJ(6) * t393;
t81 = t194 * pkin(3) + t355;
t80 = -pkin(3) * t479 - t85;
t69 = mrSges(4,1) * t137 - mrSges(4,2) * t136;
t68 = -mrSges(5,2) * t137 + mrSges(5,3) * t136;
t55 = -pkin(5) * t196 + qJ(6) * t197 + t110;
t38 = -pkin(5) * t280 - t44;
t37 = qJ(6) * t280 + t536;
t36 = pkin(5) * t228 - t42;
t35 = -qJ(6) * t228 + t43;
t29 = t137 * pkin(3) + t353;
t20 = mrSges(6,1) * t67 + mrSges(6,2) * t66;
t19 = mrSges(7,1) * t67 - mrSges(7,3) * t66;
t10 = pkin(5) * t105 - qJ(6) * t106 + qJD(6) * t197 + t50;
t7 = -pkin(5) * t195 - t9;
t6 = qJ(6) * t195 + qJD(6) * t280 + t8;
t5 = pkin(5) * t67 - qJ(6) * t66 - qJD(6) * t173 + t21;
t11 = [m(3) * (pkin(1) ^ 2 * qJDD(1) * t608 + t174 * t288 + t175 * t287 + t264 * t265) + (t649 * t473 + t358 / 0.2e1) * qJD(2) + t535 * t485 / 0.2e1 + (mrSges(6,2) * t83 + mrSges(7,2) * t24 - mrSges(6,3) * t26 - mrSges(7,3) * t39 + Ifges(6,4) * t590 + Ifges(7,5) * t589 + t669 * t582 + t670 * t587 + t656 / 0.2e1) * t106 + m(6) * (t110 * t21 + t26 * t9 + t27 * t8 + t3 * t536 + t4 * t44 + t50 * t83) + t536 * t33 + (t608 * qJD(1) * (-Ifges(3,2) * t347 + t559) + t344 * t207) * t508 / 0.2e1 + (t28 * mrSges(5,1) - t40 * mrSges(4,3) - Ifges(4,4) * t594 + Ifges(5,6) * t593 + t685) * t279 + (-m(3) * t461 - m(4) * t378 - m(5) * t363 + t348 * mrSges(2,1) + t282 * mrSges(3,1) + mrSges(2,2) * t571 - mrSges(3,3) * t483 + t682 * (-t281 * pkin(4) - pkin(10) * t199 + t363) - t673 * t198 + t448 * t687 - t443 * t686 - t624 * t281 - t619 * t199) * g(1) + (Ifges(4,1) * t578 + Ifges(4,4) * t581 - Ifges(5,2) * t579 - Ifges(5,6) * t580 + Ifges(6,6) * t590 + Ifges(7,6) * t589 + t668 * t575 + t663 * t582 + t669 * t587 + t610 + t683) * t195 + (-t175 * t526 - t261 * t478 - t264 * t479 - t268 * t288 - t287 * t365) * mrSges(3,3) + (-t21 * mrSges(6,1) - t5 * mrSges(7,1) + t1 * mrSges(7,2) + t3 * mrSges(6,3) + Ifges(6,2) * t602 - Ifges(7,3) * t601 + t665 * t595 - t667 * t603 - t627) * t196 + (t123 / 0.2e1 - t126 / 0.2e1 - Ifges(4,4) * t578 + Ifges(5,6) * t579 + Ifges(5,3) * t580 - Ifges(4,2) * t581 - t666 * t575 + t696) * t194 + (t34 * mrSges(5,1) - t41 * mrSges(4,3) + Ifges(4,4) * t592 - Ifges(5,6) * t591 + t697) * t280 + (t174 * mrSges(3,3) + Ifges(3,4) * t680 - Ifges(5,4) * t593 - Ifges(4,5) * t594 - Ifges(5,5) * t591 + Ifges(3,2) * t681 + Ifges(3,6) * t679 - Ifges(4,6) * t592 - t577 * t664 + t611 - t635 / 0.2e1) * t523 + (t175 * t535 - t268 * t570 + t287 * t385) * mrSges(3,1) + (Ifges(3,1) * t365 - Ifges(3,4) * t268 + Ifges(3,5) * t385) * t473 + (-mrSges(6,2) * t21 - mrSges(7,2) * t2 + mrSges(6,3) * t4 + mrSges(7,3) * t5 - Ifges(6,4) * t602 - Ifges(7,5) * t601 - t669 * t595 - t670 * t603 + t693) * t197 + (Ifges(3,3) * t535 + (Ifges(3,5) * t347 + Ifges(3,6) * t351) * t344) * t679 + (Ifges(3,5) * t535 + (t347 * Ifges(3,1) + t559) * t344) * t680 + t362 * t681 + (-m(3) * t511 - m(4) * t387 - m(5) * t364 - mrSges(2,1) * t571 - t284 * mrSges(3,1) + t348 * mrSges(2,2) - mrSges(3,3) * t525 + t682 * (t283 * pkin(4) + pkin(10) * t203 + t364) + t673 * t202 - t448 * t151 - t443 * t150 + t624 * t283 + t619 * t203) * g(2) + t29 * t422 + (-t25 * mrSges(7,2) - t27 * mrSges(6,3) + t39 * mrSges(7,1) + t83 * mrSges(6,1) + t73 / 0.2e1 - t76 / 0.2e1 + Ifges(7,3) * t589 - Ifges(6,2) * t590 + t667 * t587 - t665 * t582) * t105 + (-t206 / 0.2e1 - t114 * mrSges(5,3) - t135 * mrSges(4,2) + t112 * mrSges(5,2) - t134 * mrSges(4,1) + Ifges(4,5) * t578 + Ifges(5,4) * t579 + Ifges(5,5) * t580 + Ifges(4,6) * t581 + t664 * t575) * t479 + t625 * t498 + t156 * t428 + (-m(3) * t261 - t618) * t266 + m(5) * (t108 * t81 + t112 * t80 + t114 * t72 + t139 * t29 + t140 * t28 + t141 * t34) + m(7) * (t1 * t37 + t10 * t39 + t2 * t38 + t24 * t7 + t25 * t6 + t5 * t55) + m(4) * (-t134 * t85 + t135 * t84 + t156 * t251 + t160 * t41 + t161 * t40) + Ifges(2,3) * qJDD(1) + t265 * t260 + t251 * t69 + t85 * t178 + t72 * t179 + t80 * t180 + t84 * t177 + t81 * t159 + t160 * t96 + t161 * t97 + t139 * t68 + t140 * t98 + t141 * t99 + t6 * t115 + t8 * t116 + t9 * t117 + t7 * t118 + t110 * t20 - pkin(1) * t285 * t497 + t37 * t30 + t38 * t32 + t44 * t31 + t55 * t19 + t10 * t91 + t50 * t92 + (-t174 * t535 - t288 * t385 - t365 * t570) * mrSges(3,2); (-t98 + t97) * t341 + (-t260 + t445) * t261 + (-t143 + t496) * t179 + (-t495 - t163) * t178 + (t27 * t440 - t647 * t83) * mrSges(6,2) + (-t496 - t164) * t177 + (t229 * (Ifges(4,5) * t347 + t351 * t420) + t228 * (Ifges(5,5) * t347 + t351 * t403) + t649 * t347) * t470 + ((-t415 / 0.2e1 + t403 / 0.2e1) * t228 + (t420 / 0.2e1 - t404 / 0.2e1) * t229 + t640 * t575 + (Ifges(7,6) * t350 + t346 * t405) * t589 + (Ifges(6,6) * t350 + t346 * t412) * t590 + (-t346 * t638 + t350 * t669) * t587 + (t346 * t641 + t350 * t663) * t582 - t642 + ((t112 * t350 + t114 * t346) * m(5) + (t134 * t350 - t135 * t346) * m(4)) * pkin(9)) * qJD(3) + (t229 * (Ifges(5,4) * t347 + t351 * t404) + t228 * (Ifges(4,6) * t347 + t351 * t415) + t347 * t206 + (t347 * t664 + t351 * t640) * t311) * t469 - (t126 + t539) * t505 / 0.2e1 + t656 * (t345 * t465 + t349 * t462 - t218 / 0.2e1) + (t187 * t4 + t21 * t316 + t26 * t652 + t27 * t651 + t3 * t646 + t644 * t83) * m(6) + t646 * t33 - t549 * t593 + t557 * t594 + (t683 + t684) * t504 + (t125 * t434 + t21 * t426 - t405 * t601 - t412 * t602 + t5 * t424 + t435 * t622 - t595 * t641 + t603 * t638 - t685) * t350 + (-Ifges(3,2) * t482 + t207 + t326) * t435 + (pkin(9) * t631 + t645 * t108 - t112 * t145 - t114 * t143 + t29 * t304) * m(5) + (-pkin(2) * t156 + pkin(9) * t632 + t134 * t163 - t135 * t164) * m(4) + (t134 * (mrSges(4,1) * t347 - mrSges(4,3) * t515) - t112 * (mrSges(5,1) * t515 + mrSges(5,2) * t347) - t135 * (-mrSges(4,2) * t347 - mrSges(4,3) * t518) - t114 * (mrSges(5,1) * t518 - mrSges(5,3) * t347)) * t510 + (t123 * t435 + (-t96 + t99) * pkin(9) + t443 * t491 * g(3) + t126 * t434 + t697) * t346 + (-t25 * t440 + t39 * t647) * mrSges(7,3) + (-t145 + t495) * t180 + (t24 * t440 - t39 * t648) * mrSges(7,1) + (t123 + t538) * t465 + (-m(4) * t512 - m(5) * t444 + t285 - t671 * t490 + t682 * (pkin(4) * t526 + pkin(10) * t490 + t444) + (t347 * t674 + t351 * t636) * t344 - t448 * t233 - t443 * t345 * t526) * g(3) + (-t26 * t440 - t648 * t83) * mrSges(6,1) + (-m(4) * t458 - m(5) * t388 + t671 * t527 + t682 * (t284 * pkin(4) - pkin(10) * t527 + t388) - t448 * (-t283 * t522 + t284 * t349) - t443 * (t283 * t519 + t284 * t345) + t624 * t284 + t621 * t283) * g(1) + (-m(4) * t459 - m(5) * t389 + t671 * t528 + t682 * (t282 * pkin(4) - pkin(10) * t528 + t389) - t448 * (-t281 * t522 + t282 * t349) - t443 * (t281 * t519 + t282 * t345) + t624 * t282 + t621 * t281) * g(2) + t651 * t116 + t652 * t117 + t655 * t91 + t345 * t73 * t462 + t29 * t421 + (-t1 * t517 - t2 * t521 - t24 * t647 + t25 * t648) * mrSges(7,2) + (t26 * t647 + t27 * t648 - t3 * t517 + t4 * t521) * mrSges(6,3) + (-t406 * t589 - t413 * t590 - t582 * t639 - t587 * t637) * t500 + t642 * t481 + t644 * t92 + t645 * t159 + t627 * t517 + t630 * t505 + t631 * mrSges(5,1) + t632 * mrSges(4,3) + (-t73 / 0.2e1 - Ifges(6,2) * t589 + Ifges(7,3) * t590 + t600 + t667 * t588 - t665 * t583) * t217 + (t218 * t669 + t440 * t663) * t583 + (t218 * t670 + t440 * t669) * t588 + t657 * t118 + t658 * t115 + (t1 * t176 + t181 * t2 + t222 * t5 + t24 * t657 + t25 * t658 + t39 * t655) * m(7) - t156 * t427 - t550 * t591 + t558 * t592 + (t446 + t618) * t264 + t485 + (Ifges(6,4) * t218 + Ifges(6,6) * t440) * t589 + (Ifges(7,5) * t218 + Ifges(7,6) * t440) * t590 + t521 * t693 + t316 * t20 + t304 * t68 + t222 * t19 + t181 * t32 + t187 * t31 - t174 * mrSges(3,2) + t175 * mrSges(3,1) + t176 * t30 - pkin(2) * t69 + (-t358 / 0.2e1 - t625 * qJD(1)) * qJD(1) + t477 * t600; (t20 - t98) * qJ(4) + (-Ifges(4,1) * t579 + Ifges(5,2) * t578 - Ifges(6,6) * t589 - Ifges(7,6) * t590 - t574 * t668 - t583 * t663 - t588 * t669 + t610 + t690) * t228 + (t543 + t126) * t578 + t383 * t92 + (-(t412 / 0.2e1 - t405 / 0.2e1) * t393 + t638 * t587 + t695) * qJD(5) + (t178 - t180) * t135 + (-t659 * t599 + t661 / 0.2e1) * t349 + (-t24 * t36 - t25 * t35 + t5 * t302 - t361 * t599 + t39 * t653) * m(7) + (-t599 * t660 + t627) * t345 + (-t475 - t43) * t116 + (-t475 - t35) * t115 + (-t476 - t36) * t118 + (t476 - t42) * t117 - t611 + (-t544 + t123) * t579 + (t21 * qJ(4) - t26 * t42 - t27 * t43 - t599 * t360 + t689 * t83) * m(6) + t635 + (-Ifges(4,2) * t580 + Ifges(5,3) * t581 + t405 * t590 + t412 * t589 - t574 * t666 - t588 * t638 + t695 - t696) * t229 + (-t502 / 0.2e1 - t532 / 0.2e1) * t656 + (-m(5) * t451 + t422 + t428 + t682 * (-t269 - t566) + t643 * t280) * g(3) + (m(5) * t191 + t682 * (-pkin(10) * t202 - t191) + t672 * t202 + t616 * t203) * g(1) + (m(5) * t189 + t682 * (-pkin(10) * t198 - t189) + t672 * t198 + t616 * t199) * g(2) + t653 * t91 + t654 * qJD(4) + t5 * t423 + t21 * t425 + (t220 + t125) * t581 + t639 * t595 + (t26 * t633 + t27 * t634 + t615 + t629) * mrSges(6,3) + (-t24 * t633 + t25 * t634 + t615 + t628) * mrSges(7,2) + t637 * t603 + (-pkin(3) * t34 - qJ(4) * t28 - t108 * t157 - t112 * t135 + t114 * t626) * m(5) + (-t221 + t622) * t580 + (-t179 + t177) * t134 + t302 * t19 - t157 * t159 - pkin(3) * t99 + t406 * t601 + t413 * t602; t229 * t159 + (t91 + t654) * t311 + (t223 * t514 + t659) * t349 + (t223 * t513 + t660) * t345 + t99 + (t229 * t398 + t311 * t39 + t361 - t615) * m(7) + (-t229 * t397 + t311 * t83 + t360 - t615) * m(6) + (t108 * t229 - t114 * t311 + t34 - t615) * m(5); (-t513 + t561) * t27 + (-t443 * t687 - t448 * t686) * g(2) + (t173 * t25 - t24 * t393) * mrSges(7,2) - t39 * (mrSges(7,1) * t173 - mrSges(7,3) * t393) - t83 * (mrSges(6,1) * t173 + mrSges(6,2) * t393) + (-t173 * t665 + t393 * t669) * t583 + (t393 * t670 + t170 - t556 + t73) * t588 + (t150 * t448 - t151 * t443) * g(1) + (-t196 * t448 + t197 * t443) * g(3) + t612 + (-pkin(5) * t2 + qJ(6) * t1 - t24 * t27 + t25 * t650 - t39 * t90) * m(7) + (Ifges(7,3) * t173 + t553) * t590 + (-t514 + t562) * t26 + qJD(6) * t115 + (-Ifges(6,2) * t173 + t171 + t656) * t589 + qJ(6) * t30 - pkin(5) * t32 - t90 * t91 + t76 * t587 + t662; -t223 * t115 + t173 * t91 + (-g(1) * t150 + g(2) * t686 + g(3) * t196 + t39 * t173 - t25 * t223 + t2) * m(7) + t32;];
tau  = t11;
