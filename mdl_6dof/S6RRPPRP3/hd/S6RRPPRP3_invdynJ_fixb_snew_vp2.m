% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
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
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:14:37
% EndTime: 2019-05-06 09:14:41
% DurationCPUTime: 2.54s
% Computational Cost: add. (6668->300), mult. (13784->337), div. (0->0), fcn. (6817->6), ass. (0->125)
t639 = Ifges(6,4) + Ifges(7,4);
t661 = Ifges(6,2) + Ifges(7,2);
t658 = Ifges(6,6) + Ifges(7,6);
t660 = Ifges(6,1) + Ifges(7,1);
t659 = Ifges(6,5) + Ifges(7,5);
t657 = Ifges(6,3) + Ifges(7,3);
t590 = sin(qJ(5));
t593 = cos(qJ(5));
t594 = cos(qJ(2));
t626 = qJD(1) * t594;
t551 = -qJD(2) * t593 + t590 * t626;
t552 = qJD(2) * t590 + t593 * t626;
t591 = sin(qJ(2));
t627 = qJD(1) * t591;
t573 = qJD(5) + t627;
t656 = -t661 * t551 + t639 * t552 - t658 * t573;
t655 = Ifges(3,1) + Ifges(4,1) + Ifges(5,2);
t622 = Ifges(3,4) - Ifges(4,5) + Ifges(5,4);
t621 = Ifges(3,5) + Ifges(4,4) + Ifges(5,6);
t620 = Ifges(3,6) - Ifges(4,6) + Ifges(5,5);
t654 = Ifges(3,3) + Ifges(4,2) + Ifges(5,3);
t653 = Ifges(4,3) + Ifges(5,1) + Ifges(3,2);
t652 = -pkin(4) - qJ(3);
t625 = qJD(1) * qJD(2);
t611 = t591 * t625;
t559 = qJDD(1) * t594 - t611;
t564 = -qJD(2) * pkin(3) - qJ(4) * t627;
t624 = qJD(1) * qJD(4);
t586 = t594 ^ 2;
t597 = qJD(1) ^ 2;
t635 = t586 * t597;
t592 = sin(qJ(1));
t595 = cos(qJ(1));
t609 = -g(1) * t595 - g(2) * t592;
t538 = -pkin(1) * t597 + qJDD(1) * pkin(7) + t609;
t521 = -g(3) * t591 + t594 * t538;
t553 = (-pkin(2) * t594 - qJ(3) * t591) * qJD(1);
t649 = qJDD(2) * qJ(3) + t553 * t626 + t521;
t651 = pkin(3) * t635 + t559 * qJ(4) - qJD(2) * t564 + 0.2e1 * t594 * t624 - t649;
t510 = qJD(5) * t552 - qJDD(2) * t593 + t559 * t590;
t515 = -mrSges(7,2) * t573 + mrSges(7,3) * t551;
t650 = -t510 * mrSges(7,1) - t551 * t515;
t613 = t594 * t625;
t558 = qJDD(1) * t591 + t613;
t648 = -0.2e1 * t591 * t624 + (-t558 + t613) * qJ(4);
t623 = qJD(3) * qJD(2);
t578 = 0.2e1 * t623;
t596 = qJD(2) ^ 2;
t643 = pkin(2) * t596;
t492 = t578 - t643 + t649;
t554 = (-mrSges(4,1) * t594 - mrSges(4,3) * t591) * qJD(1);
t567 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t627;
t647 = m(4) * t492 + qJDD(2) * mrSges(4,3) + qJD(2) * t567 + t554 * t626;
t646 = 2 * qJD(6);
t645 = -pkin(2) - pkin(8);
t644 = pkin(3) + pkin(8);
t642 = -mrSges(4,1) + mrSges(5,2);
t640 = mrSges(4,2) - mrSges(5,3);
t634 = t594 * t597;
t628 = t592 * g(1) - t595 * g(2);
t537 = -qJDD(1) * pkin(1) - t597 * pkin(7) - t628;
t606 = -t559 * pkin(2) + t537 + (-t558 - t613) * qJ(3);
t601 = -qJ(4) * t635 + qJDD(4) - t606 + (0.2e1 * qJD(3) + t564) * t627;
t480 = t644 * t559 + pkin(4) * t558 + t601 + (pkin(4) * t594 + t591 * t645) * t625;
t557 = (pkin(4) * t591 + pkin(8) * t594) * qJD(1);
t520 = -t594 * g(3) - t591 * t538;
t610 = t553 * t627 + qJDD(3) - t520;
t484 = t652 * t596 + (-pkin(3) * t634 - qJD(1) * t557) * t591 + (-pkin(2) - t644) * qJDD(2) + t610 + t648;
t474 = t593 * t480 - t484 * t590;
t511 = qJD(5) * t551 - qJDD(2) * t590 - t559 * t593;
t513 = -mrSges(7,1) * t551 - mrSges(7,2) * t552;
t514 = -mrSges(6,1) * t551 - mrSges(6,2) * t552;
t516 = -mrSges(6,2) * t573 + mrSges(6,3) * t551;
t549 = qJDD(5) + t558;
t469 = t552 * t646 + (t551 * t573 - t511) * qJ(6) + (-t551 * t552 + t549) * pkin(5) + t474;
t619 = m(7) * t469 + t549 * mrSges(7,1) + t573 * t515;
t461 = m(6) * t474 + mrSges(6,1) * t549 + t516 * t573 + (t513 + t514) * t552 + (-mrSges(6,3) - mrSges(7,3)) * t511 + t619;
t475 = t590 * t480 + t593 * t484;
t518 = mrSges(7,1) * t573 + mrSges(7,3) * t552;
t519 = mrSges(6,1) * t573 + mrSges(6,3) * t552;
t517 = pkin(5) * t573 + qJ(6) * t552;
t548 = t551 ^ 2;
t471 = -pkin(5) * t548 + qJ(6) * t510 - t517 * t573 + t551 * t646 + t475;
t618 = m(7) * t471 + t510 * mrSges(7,3) + t551 * t513;
t464 = m(6) * t475 + mrSges(6,3) * t510 + t514 * t551 + (-t518 - t519) * t573 + (-mrSges(6,2) - mrSges(7,2)) * t549 + t618;
t458 = t593 * t461 + t590 * t464;
t633 = -t590 * t461 + t593 * t464;
t632 = t551 * t658 - t552 * t659 + t573 * t657;
t631 = -t639 * t551 + t552 * t660 - t659 * t573;
t565 = qJD(2) * mrSges(5,2) - mrSges(5,3) * t627;
t629 = -t565 - t567;
t483 = qJDD(2) * pkin(4) - t557 * t626 + t596 * t645 + t578 - t651;
t477 = -pkin(5) * t510 - qJ(6) * t548 - t517 * t552 + qJDD(6) + t483;
t617 = -m(7) * t477 - t511 * mrSges(7,2) + t552 * t518;
t616 = -t654 * qJD(2) + (-t591 * t621 - t594 * t620) * qJD(1);
t615 = -t620 * qJD(2) + (-t591 * t622 - t653 * t594) * qJD(1);
t614 = t621 * qJD(2) + (t591 * t655 + t594 * t622) * qJD(1);
t608 = -m(6) * t483 - t511 * mrSges(6,2) + t552 * t519 + t617;
t485 = -pkin(2) * t611 + pkin(3) * t559 + t601;
t568 = -qJD(2) * mrSges(5,1) + mrSges(5,3) * t626;
t607 = m(5) * t485 + t558 * mrSges(5,1) - t568 * t626 + t458;
t488 = (pkin(2) * qJD(2) - 0.2e1 * qJD(3)) * t627 + t606;
t605 = m(4) * t488 - t607;
t493 = -qJDD(2) * pkin(2) - qJ(3) * t596 + t610;
t487 = (-t591 * t634 - qJDD(2)) * pkin(3) + t493 + t648;
t556 = (mrSges(5,1) * t591 - mrSges(5,2) * t594) * qJD(1);
t603 = m(5) * t487 + qJDD(2) * mrSges(5,2) + qJD(2) * t568 - t556 * t627 + t633;
t486 = -0.2e1 * t623 + t643 + t651;
t602 = -m(5) * t486 + qJDD(2) * mrSges(5,1) - t559 * mrSges(5,3) + qJD(2) * t565 - t608;
t570 = mrSges(4,2) * t626 + qJD(2) * mrSges(4,3);
t600 = m(4) * t493 - qJDD(2) * mrSges(4,1) - qJD(2) * t570 + t603;
t466 = -mrSges(7,3) * t511 + t513 * t552 + t619;
t599 = mrSges(6,1) * t474 + mrSges(7,1) * t469 - mrSges(6,2) * t475 - mrSges(7,2) * t471 + pkin(5) * t466 + t510 * t658 + t511 * t659 + t549 * t657 + t551 * t631 + t552 * t656;
t598 = -t510 * mrSges(6,1) - t551 * t516 + t602 + t650;
t569 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t626;
t566 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t627;
t555 = (-mrSges(3,1) * t594 + mrSges(3,2) * t591) * qJD(1);
t472 = -t617 + t650;
t457 = mrSges(6,2) * t483 + mrSges(7,2) * t477 - mrSges(6,3) * t474 - mrSges(7,3) * t469 - qJ(6) * t466 + t639 * t510 + t511 * t660 + t659 * t549 + t632 * t551 + t656 * t573;
t456 = -mrSges(5,3) * t558 + t603;
t455 = -mrSges(5,2) * t559 + t565 * t627 + t607;
t454 = t554 * t627 + t558 * t640 + t600;
t453 = -mrSges(4,3) * t558 + t642 * t559 + (-t570 * t594 + t629 * t591) * qJD(1) + t605;
t452 = -mrSges(6,1) * t483 + mrSges(6,3) * t475 - mrSges(7,1) * t477 + mrSges(7,3) * t471 - pkin(5) * t472 + qJ(6) * t618 + (-qJ(6) * t518 - t631) * t573 + t632 * t552 + (-mrSges(7,2) * qJ(6) + t658) * t549 + t639 * t511 + t661 * t510;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t628 - mrSges(2,2) * t609 + t591 * (mrSges(5,1) * t485 + mrSges(3,2) * t537 + mrSges(4,2) * t493 - mrSges(3,3) * t520 - mrSges(4,3) * t488 - mrSges(5,3) * t487 + pkin(4) * t458 - qJ(3) * t453 - qJ(4) * t456 + t615 * qJD(2) + t621 * qJDD(2) + t655 * t558 + t622 * t559 - t616 * t626 + t599) + t594 * (-t593 * t457 - qJ(4) * t598 + t590 * t452 - mrSges(3,1) * t537 + mrSges(3,3) * t521 - mrSges(5,2) * t485 + mrSges(5,3) * t486 - mrSges(4,1) * t488 + mrSges(4,2) * t492 + pkin(8) * t458 + pkin(3) * t455 - pkin(2) * t453 + t653 * t559 + t622 * t558 + t620 * qJDD(2) + t614 * qJD(2) + (qJ(4) * t556 * t594 + t591 * t616) * qJD(1)) + pkin(1) * (-m(3) * t537 + (-mrSges(3,2) + mrSges(4,3)) * t558 + (mrSges(3,1) - t642) * t559 + ((t569 + t570) * t594 + (-t566 - t629) * t591) * qJD(1) - t605) + pkin(7) * (t594 * (m(3) * t521 - qJDD(2) * mrSges(3,2) - qJD(2) * t566 + t598 + (mrSges(3,3) + mrSges(4,2)) * t559 + t647) + (t555 - t556) * t586 * qJD(1) + (-m(3) * t520 - qJDD(2) * mrSges(3,1) - qJD(2) * t569 + t600 - (-mrSges(3,3) - t640) * t558 + (t554 + t555) * t627) * t591); -t593 * t452 + qJ(3) * (t602 + t647) - t590 * t457 - pkin(4) * t608 + mrSges(3,1) * t520 - mrSges(3,2) * t521 - mrSges(4,1) * t493 - mrSges(5,1) * t486 + mrSges(5,2) * t487 + mrSges(4,3) * t492 - pkin(8) * t633 - pkin(3) * t456 - pkin(2) * t454 + (qJ(3) * mrSges(4,2) + t620) * t559 + t621 * t558 + t654 * qJDD(2) + (-t615 * t591 + (-qJ(3) * t556 - t614) * t594) * qJD(1) + ((mrSges(6,1) + mrSges(7,1)) * t510 + (t515 + t516) * t551) * t652; t454; t455; t599; t472;];
tauJ  = t1;
