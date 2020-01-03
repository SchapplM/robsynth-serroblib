% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPP4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:55:12
% EndTime: 2019-12-31 20:55:17
% DurationCPUTime: 4.66s
% Computational Cost: add. (49233->290), mult. (110242->358), div. (0->0), fcn. (74646->8), ass. (0->113)
t615 = Ifges(5,1) + Ifges(6,1);
t609 = Ifges(5,4) - Ifges(6,5);
t608 = Ifges(6,4) + Ifges(5,5);
t614 = Ifges(5,2) + Ifges(6,3);
t613 = -Ifges(6,2) - Ifges(5,3);
t607 = Ifges(5,6) - Ifges(6,6);
t612 = -2 * qJD(4);
t585 = qJD(1) ^ 2;
t611 = pkin(2) * t585;
t610 = -mrSges(5,3) - mrSges(6,2);
t606 = cos(pkin(8));
t581 = sin(qJ(1));
t584 = cos(qJ(1));
t570 = -t584 * g(1) - t581 * g(2);
t559 = -t585 * pkin(1) + qJDD(1) * pkin(6) + t570;
t580 = sin(qJ(2));
t605 = t580 * t559;
t583 = cos(qJ(2));
t597 = qJD(1) * qJD(2);
t564 = t580 * qJDD(1) + t583 * t597;
t523 = qJDD(2) * pkin(2) - t564 * pkin(7) - t605 + (pkin(7) * t597 + t580 * t611 - g(3)) * t583;
t546 = -t580 * g(3) + t583 * t559;
t565 = t583 * qJDD(1) - t580 * t597;
t599 = qJD(1) * t580;
t568 = qJD(2) * pkin(2) - pkin(7) * t599;
t577 = t583 ^ 2;
t524 = t565 * pkin(7) - qJD(2) * t568 - t577 * t611 + t546;
t579 = sin(qJ(3));
t582 = cos(qJ(3));
t499 = t582 * t523 - t579 * t524;
t556 = (-t579 * t580 + t582 * t583) * qJD(1);
t530 = t556 * qJD(3) + t582 * t564 + t579 * t565;
t557 = (t579 * t583 + t580 * t582) * qJD(1);
t575 = qJDD(2) + qJDD(3);
t576 = qJD(2) + qJD(3);
t494 = (t556 * t576 - t530) * qJ(4) + (t556 * t557 + t575) * pkin(3) + t499;
t500 = t579 * t523 + t582 * t524;
t529 = -t557 * qJD(3) - t579 * t564 + t582 * t565;
t548 = t576 * pkin(3) - t557 * qJ(4);
t552 = t556 ^ 2;
t496 = -t552 * pkin(3) + t529 * qJ(4) - t576 * t548 + t500;
t578 = sin(pkin(8));
t542 = -t606 * t556 + t578 * t557;
t492 = t578 * t494 + t606 * t496 + t542 * t612;
t505 = -t606 * t529 + t578 * t530;
t543 = t578 * t556 + t606 * t557;
t533 = t576 * mrSges(5,1) - t543 * mrSges(5,3);
t517 = t542 * pkin(4) - t543 * qJ(5);
t574 = t576 ^ 2;
t487 = -t574 * pkin(4) + t575 * qJ(5) + 0.2e1 * qJD(5) * t576 - t542 * t517 + t492;
t534 = -t576 * mrSges(6,1) + t543 * mrSges(6,2);
t596 = m(6) * t487 + t575 * mrSges(6,3) + t576 * t534;
t518 = t542 * mrSges(6,1) - t543 * mrSges(6,3);
t600 = -t542 * mrSges(5,1) - t543 * mrSges(5,2) - t518;
t482 = m(5) * t492 - t575 * mrSges(5,2) + t610 * t505 - t576 * t533 + t600 * t542 + t596;
t589 = t606 * t494 - t578 * t496;
t491 = t543 * t612 + t589;
t506 = t578 * t529 + t606 * t530;
t532 = -t576 * mrSges(5,2) - t542 * mrSges(5,3);
t488 = -t575 * pkin(4) - t574 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t517) * t543 - t589;
t535 = -t542 * mrSges(6,2) + t576 * mrSges(6,3);
t591 = -m(6) * t488 + t575 * mrSges(6,1) + t576 * t535;
t484 = m(5) * t491 + t575 * mrSges(5,1) + t610 * t506 + t576 * t532 + t600 * t543 + t591;
t477 = t578 * t482 + t606 * t484;
t544 = -t556 * mrSges(4,1) + t557 * mrSges(4,2);
t547 = -t576 * mrSges(4,2) + t556 * mrSges(4,3);
t473 = m(4) * t499 + t575 * mrSges(4,1) - t530 * mrSges(4,3) - t557 * t544 + t576 * t547 + t477;
t549 = t576 * mrSges(4,1) - t557 * mrSges(4,3);
t592 = t606 * t482 - t578 * t484;
t474 = m(4) * t500 - t575 * mrSges(4,2) + t529 * mrSges(4,3) + t556 * t544 - t576 * t549 + t592;
t469 = t582 * t473 + t579 * t474;
t545 = -t583 * g(3) - t605;
t563 = (-mrSges(3,1) * t583 + mrSges(3,2) * t580) * qJD(1);
t598 = qJD(1) * t583;
t567 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t598;
t467 = m(3) * t545 + qJDD(2) * mrSges(3,1) - t564 * mrSges(3,3) + qJD(2) * t567 - t563 * t599 + t469;
t566 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t599;
t593 = -t579 * t473 + t582 * t474;
t468 = m(3) * t546 - qJDD(2) * mrSges(3,2) + t565 * mrSges(3,3) - qJD(2) * t566 + t563 * t598 + t593;
t594 = -t580 * t467 + t583 * t468;
t460 = m(2) * t570 - t585 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t594;
t569 = t581 * g(1) - t584 * g(2);
t590 = -qJDD(1) * pkin(1) - t569;
t558 = -t585 * pkin(6) + t590;
t531 = -t565 * pkin(2) + t568 * t599 + (-pkin(7) * t577 - pkin(6)) * t585 + t590;
t498 = -t529 * pkin(3) - t552 * qJ(4) + t557 * t548 + qJDD(4) + t531;
t490 = -0.2e1 * qJD(5) * t543 + (t542 * t576 - t506) * qJ(5) + (t543 * t576 + t505) * pkin(4) + t498;
t485 = m(6) * t490 + t505 * mrSges(6,1) - t506 * mrSges(6,3) - t543 * t534 + t542 * t535;
t588 = m(5) * t498 + t505 * mrSges(5,1) + t506 * mrSges(5,2) + t542 * t532 + t543 * t533 + t485;
t587 = m(4) * t531 - t529 * mrSges(4,1) + t530 * mrSges(4,2) - t556 * t547 + t557 * t549 + t588;
t586 = -m(3) * t558 + t565 * mrSges(3,1) - t564 * mrSges(3,2) - t566 * t599 + t567 * t598 - t587;
t479 = m(2) * t569 + qJDD(1) * mrSges(2,1) - t585 * mrSges(2,2) + t586;
t604 = t581 * t460 + t584 * t479;
t461 = t583 * t467 + t580 * t468;
t603 = t614 * t542 - t609 * t543 - t607 * t576;
t602 = t607 * t542 - t608 * t543 + t613 * t576;
t601 = -t609 * t542 + t615 * t543 + t608 * t576;
t595 = t584 * t460 - t581 * t479;
t555 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t580 + Ifges(3,4) * t583) * qJD(1);
t554 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t580 + Ifges(3,2) * t583) * qJD(1);
t553 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t580 + Ifges(3,6) * t583) * qJD(1);
t539 = Ifges(4,1) * t557 + Ifges(4,4) * t556 + Ifges(4,5) * t576;
t538 = Ifges(4,4) * t557 + Ifges(4,2) * t556 + Ifges(4,6) * t576;
t537 = Ifges(4,5) * t557 + Ifges(4,6) * t556 + Ifges(4,3) * t576;
t476 = mrSges(5,2) * t498 + mrSges(6,2) * t488 - mrSges(5,3) * t491 - mrSges(6,3) * t490 - qJ(5) * t485 - t609 * t505 + t615 * t506 + t602 * t542 + t608 * t575 + t603 * t576;
t475 = -mrSges(5,1) * t498 - mrSges(6,1) * t490 + mrSges(6,2) * t487 + mrSges(5,3) * t492 - pkin(4) * t485 - t614 * t505 + t609 * t506 + t602 * t543 + t607 * t575 + t601 * t576;
t463 = mrSges(4,2) * t531 - mrSges(4,3) * t499 + Ifges(4,1) * t530 + Ifges(4,4) * t529 + Ifges(4,5) * t575 - qJ(4) * t477 - t578 * t475 + t606 * t476 + t556 * t537 - t576 * t538;
t462 = -mrSges(4,1) * t531 + mrSges(4,3) * t500 + Ifges(4,4) * t530 + Ifges(4,2) * t529 + Ifges(4,6) * t575 - pkin(3) * t588 + qJ(4) * t592 + t606 * t475 + t578 * t476 - t557 * t537 + t576 * t539;
t457 = mrSges(3,2) * t558 - mrSges(3,3) * t545 + Ifges(3,1) * t564 + Ifges(3,4) * t565 + Ifges(3,5) * qJDD(2) - pkin(7) * t469 - qJD(2) * t554 - t579 * t462 + t582 * t463 + t553 * t598;
t456 = (pkin(4) * mrSges(6,2) - t608) * t506 + (qJ(5) * mrSges(6,2) + t607) * t505 + (qJ(5) * t518 - t601) * t542 + (pkin(4) * t518 + t603) * t543 - pkin(4) * t591 + Ifges(2,6) * qJDD(1) + (-Ifges(4,3) + t613) * t575 - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + (-t580 * t554 + t583 * t555) * qJD(1) + t585 * Ifges(2,5) - Ifges(3,6) * t565 + mrSges(2,3) * t570 - t557 * t538 - Ifges(3,5) * t564 - mrSges(3,1) * t545 + mrSges(3,2) * t546 + t556 * t539 - Ifges(4,6) * t529 - Ifges(4,5) * t530 - mrSges(4,1) * t499 + mrSges(4,2) * t500 + mrSges(5,2) * t492 + mrSges(6,1) * t488 - mrSges(5,1) * t491 - mrSges(6,3) * t487 - pkin(3) * t477 - pkin(2) * t469 - pkin(1) * t461 - qJ(5) * t596;
t455 = -mrSges(3,1) * t558 + mrSges(3,3) * t546 + Ifges(3,4) * t564 + Ifges(3,2) * t565 + Ifges(3,6) * qJDD(2) - pkin(2) * t587 + pkin(7) * t593 + qJD(2) * t555 + t582 * t462 + t579 * t463 - t553 * t599;
t454 = -mrSges(2,2) * g(3) - mrSges(2,3) * t569 + Ifges(2,5) * qJDD(1) - t585 * Ifges(2,6) - pkin(6) * t461 - t580 * t455 + t583 * t457;
t1 = [-m(1) * g(1) + t595; -m(1) * g(2) + t604; (-m(1) - m(2)) * g(3) + t461; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t604 + t584 * t454 - t581 * t456; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t595 + t581 * t454 + t584 * t456; -mrSges(1,1) * g(2) + mrSges(2,1) * t569 + mrSges(1,2) * g(1) - mrSges(2,2) * t570 + Ifges(2,3) * qJDD(1) + pkin(1) * t586 + pkin(6) * t594 + t583 * t455 + t580 * t457;];
tauB = t1;
