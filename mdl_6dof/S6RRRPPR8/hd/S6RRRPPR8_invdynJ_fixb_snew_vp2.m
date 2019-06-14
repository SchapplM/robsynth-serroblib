% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-07 06:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRRPPR8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPPR8_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:09:53
% EndTime: 2019-05-07 06:10:01
% DurationCPUTime: 3.66s
% Computational Cost: add. (31174->323), mult. (67243->389), div. (0->0), fcn. (50251->10), ass. (0->140)
t599 = cos(pkin(6));
t595 = qJD(1) * t599 + qJD(2);
t601 = sin(qJ(3));
t602 = sin(qJ(2));
t598 = sin(pkin(6));
t635 = qJD(1) * t598;
t625 = t602 * t635;
t646 = cos(qJ(3));
t567 = t601 * t595 + t625 * t646;
t605 = cos(qJ(2));
t634 = qJD(1) * t605;
t624 = t598 * t634;
t588 = -qJD(3) + t624;
t549 = pkin(4) * t588 - qJ(5) * t567;
t656 = (2 * qJD(4)) + t549;
t655 = Ifges(4,1) + Ifges(5,1) + Ifges(6,2);
t632 = Ifges(4,5) + Ifges(5,4) + Ifges(6,6);
t631 = Ifges(5,5) - Ifges(6,4) - Ifges(4,4);
t630 = Ifges(5,6) - Ifges(6,5) - Ifges(4,6);
t654 = Ifges(5,3) + Ifges(6,1) + Ifges(4,2);
t653 = Ifges(6,3) + Ifges(4,3) + Ifges(5,2);
t566 = -t595 * t646 + t601 * t625;
t581 = (qJD(2) * t634 + qJDD(1) * t602) * t598;
t594 = qJDD(1) * t599 + qJDD(2);
t533 = -t566 * qJD(3) + t581 * t646 + t601 * t594;
t536 = mrSges(6,1) * t567 + mrSges(6,2) * t566;
t538 = mrSges(5,1) * t566 - mrSges(5,3) * t567;
t633 = qJDD(1) * t598;
t582 = -qJD(2) * t625 + t605 * t633;
t574 = -qJDD(3) + t582;
t580 = (-pkin(2) * t605 - pkin(9) * t602) * t635;
t593 = t595 ^ 2;
t607 = qJD(1) ^ 2;
t603 = sin(qJ(1));
t606 = cos(qJ(1));
t621 = -g(1) * t606 - g(2) * t603;
t578 = -pkin(1) * t607 + pkin(8) * t633 + t621;
t623 = t603 * g(1) - g(2) * t606;
t644 = pkin(8) * t598;
t577 = qJDD(1) * pkin(1) + t607 * t644 + t623;
t640 = t577 * t599;
t636 = t605 * t578 + t602 * t640;
t501 = -pkin(2) * t593 + pkin(9) * t594 + (-g(3) * t602 + t580 * t634) * t598 + t636;
t643 = g(3) * t599;
t502 = -pkin(2) * t582 - pkin(9) * t581 - t643 + (-t577 + (pkin(2) * t602 - pkin(9) * t605) * t595 * qJD(1)) * t598;
t485 = -t601 * t501 + t502 * t646;
t537 = pkin(3) * t566 - qJ(4) * t567;
t651 = t588 ^ 2;
t617 = -qJ(4) * t651 + t567 * t537 + qJDD(4) - t485;
t483 = t574 * pkin(3) + t617;
t553 = -mrSges(5,2) * t566 - mrSges(5,3) * t588;
t532 = qJD(3) * t567 + t581 * t601 - t594 * t646;
t565 = t566 ^ 2;
t638 = t598 * t605;
t534 = -g(3) * t638 - t602 * t578 + t605 * t640;
t500 = -t594 * pkin(2) - t593 * pkin(9) + t580 * t625 - t534;
t641 = t566 * t588;
t618 = t532 * pkin(3) + t500 + (-t533 - t641) * qJ(4);
t612 = -qJ(5) * t565 + t656 * t567 + qJDD(5) - t618;
t647 = -pkin(4) - pkin(10);
t473 = t612 + (pkin(5) * t566 + (pkin(3) + pkin(10)) * t567) * t588 + t647 * t532 + pkin(5) * t533;
t540 = pkin(5) * t567 - pkin(10) * t566;
t648 = -2 * qJD(5);
t610 = (-t533 + t641) * qJ(5) + t617 + (pkin(4) * t566 + t648) * t567;
t474 = -t651 * pkin(5) - t567 * t540 + (pkin(3) - t647) * t574 + t610;
t600 = sin(qJ(6));
t604 = cos(qJ(6));
t470 = t473 * t604 - t474 * t600;
t545 = -t566 * t600 + t588 * t604;
t492 = qJD(6) * t545 + t532 * t604 + t574 * t600;
t546 = t566 * t604 + t588 * t600;
t503 = -mrSges(7,1) * t545 + mrSges(7,2) * t546;
t564 = qJD(6) + t567;
t507 = -mrSges(7,2) * t564 + mrSges(7,3) * t545;
t525 = qJDD(6) + t533;
t467 = m(7) * t470 + mrSges(7,1) * t525 - mrSges(7,3) * t492 - t503 * t546 + t507 * t564;
t471 = t473 * t600 + t474 * t604;
t491 = -qJD(6) * t546 - t532 * t600 + t574 * t604;
t508 = mrSges(7,1) * t564 - mrSges(7,3) * t546;
t468 = m(7) * t471 - mrSges(7,2) * t525 + mrSges(7,3) * t491 + t503 * t545 - t508 * t564;
t457 = -t600 * t467 + t468 * t604;
t477 = (pkin(3) + pkin(4)) * t574 + t610;
t547 = mrSges(6,1) * t588 - mrSges(6,3) * t566;
t619 = m(6) * t477 - t574 * mrSges(6,2) - t588 * t547 + t457;
t614 = m(5) * t483 + t574 * mrSges(5,1) + t588 * t553 + t619;
t453 = (-t536 + t538) * t567 + (mrSges(5,2) - mrSges(6,3)) * t533 + t614;
t454 = -mrSges(6,3) * t533 - t536 * t567 + t619;
t650 = -2 * qJD(4);
t583 = t588 * t650;
t486 = t646 * t501 + t601 * t502;
t620 = pkin(3) * t651 + t574 * qJ(4) + t566 * t537 - t486;
t615 = pkin(4) * t565 - qJ(5) * t532 + t620;
t476 = -pkin(5) * t574 - pkin(10) * t651 - t549 * t588 + t583 + ((2 * qJD(5)) + t540) * t566 - t615;
t493 = Ifges(7,5) * t546 + Ifges(7,6) * t545 + Ifges(7,3) * t564;
t495 = Ifges(7,1) * t546 + Ifges(7,4) * t545 + Ifges(7,5) * t564;
t461 = -mrSges(7,1) * t476 + mrSges(7,3) * t471 + Ifges(7,4) * t492 + Ifges(7,2) * t491 + Ifges(7,6) * t525 - t493 * t546 + t495 * t564;
t494 = Ifges(7,4) * t546 + Ifges(7,2) * t545 + Ifges(7,6) * t564;
t462 = mrSges(7,2) * t476 - mrSges(7,3) * t470 + Ifges(7,1) * t492 + Ifges(7,4) * t491 + Ifges(7,5) * t525 + t493 * t545 - t494 * t564;
t472 = -m(7) * t476 + mrSges(7,1) * t491 - t492 * mrSges(7,2) + t507 * t545 - t546 * t508;
t478 = t566 * t648 + t656 * t588 + t615;
t482 = t583 - t620;
t552 = -mrSges(6,2) * t588 - mrSges(6,3) * t567;
t551 = mrSges(5,1) * t588 + mrSges(5,2) * t567;
t613 = -m(6) * t478 + t532 * mrSges(6,3) + t566 * t536 - t472;
t611 = m(5) * t482 - t574 * mrSges(5,3) - t588 * t551 + t613;
t626 = t631 * t566 + t655 * t567 - t632 * t588;
t628 = t654 * t566 + t631 * t567 - t630 * t588;
t652 = t532 * t630 + t533 * t632 + t566 * t626 - t567 * t628 - t653 * t574 + mrSges(4,1) * t485 - mrSges(5,1) * t483 - mrSges(6,1) * t478 - mrSges(4,2) * t486 + mrSges(6,2) * t477 + mrSges(5,3) * t482 - pkin(3) * t453 - pkin(4) * t454 - pkin(5) * t472 - pkin(10) * t457 + qJ(4) * (-mrSges(6,1) * t574 - mrSges(5,2) * t532 - t538 * t566 - t552 * t588 + t611) - t604 * t461 - t600 * t462;
t645 = pkin(3) * t588;
t642 = -mrSges(4,3) - mrSges(5,2);
t639 = t598 * t602;
t548 = mrSges(4,2) * t588 - mrSges(4,3) * t566;
t637 = -mrSges(4,1) * t566 - mrSges(4,2) * t567 - t538;
t451 = m(4) * t485 - mrSges(4,1) * t574 - t548 * t588 + (t536 + t637) * t567 + (mrSges(6,3) + t642) * t533 - t614;
t550 = -mrSges(4,1) * t588 - mrSges(4,3) * t567;
t460 = t611 + (t550 - t552) * t588 + t637 * t566 + (mrSges(4,2) - mrSges(6,1)) * t574 + t642 * t532 + m(4) * t486;
t448 = t646 * t451 + t601 * t460;
t456 = t604 * t467 + t600 * t468;
t627 = -t630 * t566 - t632 * t567 + t653 * t588;
t622 = -t451 * t601 + t646 * t460;
t480 = -pkin(4) * t532 + t567 * t645 + t612;
t455 = m(6) * t480 + t533 * mrSges(6,1) + t532 * mrSges(6,2) + t566 * t547 + t567 * t552 + t456;
t616 = mrSges(7,1) * t470 - mrSges(7,2) * t471 + Ifges(7,5) * t492 + Ifges(7,6) * t491 + Ifges(7,3) * t525 + t546 * t494 - t545 * t495;
t484 = (t650 - t645) * t567 + t618;
t452 = m(5) * t484 + t532 * mrSges(5,1) - t533 * mrSges(5,3) - t567 * t551 + t566 * t553 - t455;
t609 = -m(4) * t500 - t532 * mrSges(4,1) - t533 * mrSges(4,2) - t566 * t548 - t567 * t550 - t452;
t579 = (-mrSges(3,1) * t605 + mrSges(3,2) * t602) * t635;
t576 = -mrSges(3,2) * t595 + mrSges(3,3) * t624;
t575 = mrSges(3,1) * t595 - mrSges(3,3) * t625;
t558 = -t577 * t598 - t643;
t557 = Ifges(3,5) * t595 + (Ifges(3,1) * t602 + Ifges(3,4) * t605) * t635;
t556 = Ifges(3,6) * t595 + (Ifges(3,4) * t602 + Ifges(3,2) * t605) * t635;
t555 = Ifges(3,3) * t595 + (Ifges(3,5) * t602 + Ifges(3,6) * t605) * t635;
t535 = -g(3) * t639 + t636;
t449 = m(3) * t534 + t594 * mrSges(3,1) - t581 * mrSges(3,3) + t595 * t576 - t579 * t625 + t609;
t447 = m(3) * t535 - mrSges(3,2) * t594 + mrSges(3,3) * t582 - t575 * t595 + t579 * t624 + t622;
t446 = mrSges(6,1) * t480 + mrSges(4,2) * t500 + mrSges(5,2) * t483 - mrSges(4,3) * t485 - mrSges(5,3) * t484 - mrSges(6,3) * t477 + pkin(5) * t456 - qJ(4) * t452 - qJ(5) * t454 + t631 * t532 + t655 * t533 + t627 * t566 - t632 * t574 - t628 * t588 + t616;
t445 = -t604 * t462 + t600 * t461 - qJ(5) * t613 - mrSges(4,1) * t500 - mrSges(5,1) * t484 + mrSges(4,3) * t486 - mrSges(6,2) * t480 + mrSges(5,2) * t482 + mrSges(6,3) * t478 + pkin(10) * t456 + pkin(4) * t455 - pkin(3) * t452 + (qJ(5) * t552 - t626) * t588 + (mrSges(6,1) * qJ(5) + t630) * t574 + t627 * t567 - t631 * t533 - t654 * t532;
t444 = Ifges(3,5) * t581 + Ifges(3,6) * t582 + Ifges(3,3) * t594 + mrSges(3,1) * t534 - mrSges(3,2) * t535 + t601 * t446 + t646 * t445 + pkin(2) * t609 + pkin(9) * t622 + (t556 * t602 - t557 * t605) * t635;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t623 - mrSges(2,2) * t621 + (mrSges(3,2) * t558 - mrSges(3,3) * t534 + Ifges(3,1) * t581 + Ifges(3,4) * t582 + Ifges(3,5) * t594 - pkin(9) * t448 - t601 * t445 + t446 * t646 + t555 * t624 - t595 * t556) * t639 + (-mrSges(3,1) * t558 + mrSges(3,3) * t535 + Ifges(3,4) * t581 + Ifges(3,2) * t582 + Ifges(3,6) * t594 - pkin(2) * t448 - t555 * t625 + t595 * t557 - t652) * t638 + t599 * t444 + pkin(1) * ((t447 * t602 + t449 * t605) * t599 + (-m(3) * t558 + t582 * mrSges(3,1) - t581 * mrSges(3,2) + (-t575 * t602 + t576 * t605) * t635 - t448) * t598) + (t447 * t605 - t449 * t602) * t644; t444; t652; t453; t455; t616;];
tauJ  = t1;
