% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:32:11
% EndTime: 2019-12-31 18:32:15
% DurationCPUTime: 2.94s
% Computational Cost: add. (24594->272), mult. (58978->327), div. (0->0), fcn. (39736->8), ass. (0->117)
t612 = Ifges(4,1) + Ifges(5,2);
t606 = Ifges(4,4) + Ifges(5,6);
t605 = Ifges(4,5) - Ifges(5,4);
t611 = Ifges(4,2) + Ifges(5,3);
t604 = Ifges(4,6) - Ifges(5,5);
t610 = -Ifges(4,3) - Ifges(5,1);
t572 = qJD(1) ^ 2;
t609 = -2 * qJD(4);
t608 = cos(qJ(3));
t565 = cos(pkin(8));
t607 = pkin(2) * t565;
t564 = sin(pkin(8));
t603 = mrSges(3,2) * t564;
t562 = t565 ^ 2;
t602 = t562 * t572;
t568 = sin(qJ(1));
t570 = cos(qJ(1));
t552 = -t570 * g(1) - t568 * g(2);
t548 = -t572 * pkin(1) + qJDD(1) * qJ(2) + t552;
t593 = qJD(1) * qJD(2);
t589 = -t565 * g(3) - 0.2e1 * t564 * t593;
t508 = (-pkin(6) * qJDD(1) + t572 * t607 - t548) * t564 + t589;
t531 = -t564 * g(3) + (t548 + 0.2e1 * t593) * t565;
t591 = qJDD(1) * t565;
t509 = -pkin(2) * t602 + pkin(6) * t591 + t531;
t567 = sin(qJ(3));
t491 = t608 * t508 - t567 * t509;
t590 = t565 * t608;
t594 = t564 * qJD(1);
t546 = -qJD(1) * t590 + t567 * t594;
t580 = t608 * t564 + t565 * t567;
t547 = t580 * qJD(1);
t520 = t546 * mrSges(4,1) + t547 * mrSges(4,2);
t596 = t546 * qJD(3);
t529 = t580 * qJDD(1) - t596;
t535 = -qJD(3) * mrSges(4,2) - t546 * mrSges(4,3);
t537 = t546 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t592 = qJDD(1) * t564;
t595 = t547 * qJD(3);
t528 = -qJDD(1) * t590 + t567 * t592 + t595;
t539 = t547 * pkin(4) - qJD(3) * pkin(7);
t545 = t546 ^ 2;
t561 = t564 ^ 2;
t551 = t568 * g(1) - t570 * g(2);
t585 = qJDD(2) - t551;
t527 = (-pkin(1) - t607) * qJDD(1) + (-qJ(2) + (-t561 - t562) * pkin(6)) * t572 + t585;
t573 = pkin(3) * t595 + t547 * t609 + (-t529 + t596) * qJ(4) + t527;
t483 = -t545 * pkin(4) - t547 * t539 + (pkin(3) + pkin(7)) * t528 + t573;
t519 = t546 * pkin(3) - t547 * qJ(4);
t571 = qJD(3) ^ 2;
t490 = -qJDD(3) * pkin(3) - t571 * qJ(4) + t547 * t519 + qJDD(4) - t491;
t484 = (t546 * t547 - qJDD(3)) * pkin(7) + (t529 + t596) * pkin(4) + t490;
t566 = sin(qJ(5));
t569 = cos(qJ(5));
t481 = -t566 * t483 + t569 * t484;
t532 = -t566 * qJD(3) + t569 * t546;
t499 = t532 * qJD(5) + t569 * qJDD(3) + t566 * t528;
t533 = t569 * qJD(3) + t566 * t546;
t500 = -t532 * mrSges(6,1) + t533 * mrSges(6,2);
t543 = qJD(5) + t547;
t504 = -t543 * mrSges(6,2) + t532 * mrSges(6,3);
t526 = qJDD(5) + t529;
t479 = m(6) * t481 + t526 * mrSges(6,1) - t499 * mrSges(6,3) - t533 * t500 + t543 * t504;
t482 = t569 * t483 + t566 * t484;
t498 = -t533 * qJD(5) - t566 * qJDD(3) + t569 * t528;
t505 = t543 * mrSges(6,1) - t533 * mrSges(6,3);
t480 = m(6) * t482 - t526 * mrSges(6,2) + t498 * mrSges(6,3) + t532 * t500 - t543 * t505;
t471 = t569 * t479 + t566 * t480;
t521 = -t546 * mrSges(5,2) - t547 * mrSges(5,3);
t578 = -m(5) * t490 - t529 * mrSges(5,1) - t547 * t521 - t471;
t469 = m(4) * t491 - t529 * mrSges(4,3) - t547 * t520 + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t535 - t537) * qJD(3) + t578;
t492 = t567 * t508 + t608 * t509;
t536 = qJD(3) * mrSges(4,1) - t547 * mrSges(4,3);
t577 = -t571 * pkin(3) + qJDD(3) * qJ(4) - t546 * t519 + t492;
t489 = qJD(3) * t609 - t577;
t538 = t547 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t486 = -t528 * pkin(4) - t545 * pkin(7) + ((2 * qJD(4)) + t539) * qJD(3) + t577;
t579 = -m(6) * t486 + t498 * mrSges(6,1) - t499 * mrSges(6,2) + t532 * t504 - t533 * t505;
t576 = -m(5) * t489 + qJDD(3) * mrSges(5,3) + qJD(3) * t538 - t579;
t476 = m(4) * t492 - qJDD(3) * mrSges(4,2) - qJD(3) * t536 + (-t520 - t521) * t546 + (-mrSges(4,3) - mrSges(5,1)) * t528 + t576;
t465 = t608 * t469 + t567 * t476;
t530 = -t564 * t548 + t589;
t581 = mrSges(3,3) * qJDD(1) + t572 * (-mrSges(3,1) * t565 + t603);
t463 = m(3) * t530 - t581 * t564 + t465;
t586 = -t567 * t469 + t608 * t476;
t464 = m(3) * t531 + t581 * t565 + t586;
t587 = -t564 * t463 + t565 * t464;
t457 = m(2) * t552 - t572 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t587;
t544 = -qJDD(1) * pkin(1) - t572 * qJ(2) + t585;
t488 = t528 * pkin(3) + t573;
t600 = -t566 * t479 + t569 * t480;
t470 = m(5) * t488 - t528 * mrSges(5,2) - t529 * mrSges(5,3) - t546 * t537 - t547 * t538 + t600;
t575 = m(4) * t527 + t528 * mrSges(4,1) + t529 * mrSges(4,2) + t546 * t535 + t547 * t536 + t470;
t574 = -m(3) * t544 + mrSges(3,1) * t591 - t575 + (t561 * t572 + t602) * mrSges(3,3);
t467 = -t572 * mrSges(2,2) + m(2) * t551 + (mrSges(2,1) - t603) * qJDD(1) + t574;
t601 = t568 * t457 + t570 * t467;
t458 = t565 * t463 + t564 * t464;
t599 = t610 * qJD(3) + t604 * t546 - t605 * t547;
t598 = -t604 * qJD(3) + t611 * t546 - t606 * t547;
t597 = t605 * qJD(3) - t606 * t546 + t612 * t547;
t588 = t570 * t457 - t568 * t467;
t584 = Ifges(3,1) * t564 + Ifges(3,4) * t565;
t583 = Ifges(3,4) * t564 + Ifges(3,2) * t565;
t582 = Ifges(3,5) * t564 + Ifges(3,6) * t565;
t550 = t582 * qJD(1);
t495 = Ifges(6,1) * t533 + Ifges(6,4) * t532 + Ifges(6,5) * t543;
t494 = Ifges(6,4) * t533 + Ifges(6,2) * t532 + Ifges(6,6) * t543;
t493 = Ifges(6,5) * t533 + Ifges(6,6) * t532 + Ifges(6,3) * t543;
t473 = mrSges(6,2) * t486 - mrSges(6,3) * t481 + Ifges(6,1) * t499 + Ifges(6,4) * t498 + Ifges(6,5) * t526 + t532 * t493 - t543 * t494;
t472 = -mrSges(6,1) * t486 + mrSges(6,3) * t482 + Ifges(6,4) * t499 + Ifges(6,2) * t498 + Ifges(6,6) * t526 - t533 * t493 + t543 * t495;
t459 = mrSges(5,1) * t490 + mrSges(6,1) * t481 + mrSges(4,2) * t527 - mrSges(6,2) * t482 - mrSges(4,3) * t491 - mrSges(5,3) * t488 + Ifges(6,5) * t499 + Ifges(6,6) * t498 + Ifges(6,3) * t526 + pkin(4) * t471 - qJ(4) * t470 + t533 * t494 - t532 * t495 + t599 * t546 + t612 * t529 - t606 * t528 + t605 * qJDD(3) + t598 * qJD(3);
t454 = -mrSges(4,1) * t527 - mrSges(5,1) * t489 + mrSges(5,2) * t488 + mrSges(4,3) * t492 - pkin(3) * t470 - pkin(4) * t579 - pkin(7) * t600 + t597 * qJD(3) + t604 * qJDD(3) - t569 * t472 - t566 * t473 - t611 * t528 + t606 * t529 + t599 * t547;
t453 = t565 * qJD(1) * t550 + mrSges(3,2) * t544 - mrSges(3,3) * t530 - pkin(6) * t465 + t584 * qJDD(1) - t567 * t454 + t608 * t459;
t452 = t566 * t472 - pkin(3) * (-qJD(3) * t537 + t578) - t569 * t473 + mrSges(2,3) * t552 - qJ(4) * t576 - mrSges(3,1) * t530 + mrSges(3,2) * t531 + pkin(7) * t471 - pkin(2) * t465 + mrSges(2,1) * g(3) - mrSges(5,2) * t490 - mrSges(4,1) * t491 + mrSges(4,2) * t492 + mrSges(5,3) * t489 - pkin(1) * t458 + t598 * t547 + (qJ(4) * t521 - t597) * t546 - t605 * t529 + (qJ(4) * mrSges(5,1) + t604) * t528 + (pkin(3) * mrSges(5,2) + t610) * qJDD(3) + (Ifges(2,6) - t582) * qJDD(1) + (-t564 * t583 + t565 * t584 + Ifges(2,5)) * t572;
t451 = -mrSges(3,1) * t544 + mrSges(3,3) * t531 - pkin(2) * t575 + pkin(6) * t586 + t583 * qJDD(1) + t608 * t454 + t567 * t459 - t550 * t594;
t450 = -mrSges(2,2) * g(3) - mrSges(2,3) * t551 + Ifges(2,5) * qJDD(1) - t572 * Ifges(2,6) - qJ(2) * t458 - t564 * t451 + t565 * t453;
t1 = [-m(1) * g(1) + t588; -m(1) * g(2) + t601; (-m(1) - m(2)) * g(3) + t458; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t601 + t570 * t450 - t568 * t452; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t588 + t568 * t450 + t570 * t452; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t551 - mrSges(2,2) * t552 + t564 * t453 + t565 * t451 + pkin(1) * (-mrSges(3,2) * t592 + t574) + qJ(2) * t587;];
tauB = t1;
