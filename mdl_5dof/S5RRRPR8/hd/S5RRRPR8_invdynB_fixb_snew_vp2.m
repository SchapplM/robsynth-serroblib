% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR8_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR8_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:26
% EndTime: 2019-12-31 21:19:32
% DurationCPUTime: 3.47s
% Computational Cost: add. (32439->294), mult. (67352->355), div. (0->0), fcn. (43646->8), ass. (0->116)
t617 = Ifges(4,1) + Ifges(5,2);
t616 = -Ifges(5,1) - Ifges(4,3);
t611 = Ifges(4,4) + Ifges(5,6);
t610 = Ifges(4,5) - Ifges(5,4);
t615 = Ifges(4,2) + Ifges(5,3);
t609 = Ifges(4,6) - Ifges(5,5);
t614 = -2 * qJD(4);
t613 = cos(qJ(3));
t587 = qJD(1) ^ 2;
t612 = pkin(2) * t587;
t581 = sin(qJ(3));
t585 = cos(qJ(2));
t600 = qJD(1) * t585;
t582 = sin(qJ(2));
t601 = qJD(1) * t582;
t557 = t581 * t601 - t613 * t600;
t578 = qJD(2) + qJD(3);
t608 = t557 * t578;
t583 = sin(qJ(1));
t586 = cos(qJ(1));
t572 = -t586 * g(1) - t583 * g(2);
t560 = -t587 * pkin(1) + qJDD(1) * pkin(6) + t572;
t607 = t582 * t560;
t599 = qJD(1) * qJD(2);
t566 = t582 * qJDD(1) + t585 * t599;
t513 = qJDD(2) * pkin(2) - t566 * pkin(7) - t607 + (pkin(7) * t599 + t582 * t612 - g(3)) * t585;
t543 = -t582 * g(3) + t585 * t560;
t567 = t585 * qJDD(1) - t582 * t599;
t570 = qJD(2) * pkin(2) - pkin(7) * t601;
t579 = t585 ^ 2;
t514 = t567 * pkin(7) - qJD(2) * t570 - t579 * t612 + t543;
t500 = t613 * t513 - t581 * t514;
t523 = -t557 * qJD(3) + t613 * t566 + t581 * t567;
t558 = (t581 * t585 + t613 * t582) * qJD(1);
t537 = t557 * mrSges(4,1) + t558 * mrSges(4,2);
t544 = -t578 * mrSges(4,2) - t557 * mrSges(4,3);
t546 = t557 * mrSges(5,1) - t578 * mrSges(5,3);
t577 = qJDD(2) + qJDD(3);
t522 = t558 * qJD(3) + t581 * t566 - t613 * t567;
t548 = t558 * pkin(4) - t578 * pkin(8);
t553 = t557 ^ 2;
t571 = t583 * g(1) - t586 * g(2);
t595 = -qJDD(1) * pkin(1) - t571;
t524 = -t567 * pkin(2) + t570 * t601 + (-pkin(7) * t579 - pkin(6)) * t587 + t595;
t589 = (-t523 + t608) * qJ(4) + t524 + (t578 * pkin(3) + t614) * t558;
t492 = -t553 * pkin(4) - t558 * t548 + (pkin(3) + pkin(8)) * t522 + t589;
t536 = t557 * pkin(3) - t558 * qJ(4);
t576 = t578 ^ 2;
t499 = -t577 * pkin(3) - t576 * qJ(4) + t558 * t536 + qJDD(4) - t500;
t493 = (t557 * t558 - t577) * pkin(8) + (t523 + t608) * pkin(4) + t499;
t580 = sin(qJ(5));
t584 = cos(qJ(5));
t490 = -t580 * t492 + t584 * t493;
t540 = t584 * t557 - t580 * t578;
t504 = t540 * qJD(5) + t580 * t522 + t584 * t577;
t541 = t580 * t557 + t584 * t578;
t512 = -t540 * mrSges(6,1) + t541 * mrSges(6,2);
t521 = qJDD(5) + t523;
t552 = qJD(5) + t558;
t525 = -t552 * mrSges(6,2) + t540 * mrSges(6,3);
t488 = m(6) * t490 + t521 * mrSges(6,1) - t504 * mrSges(6,3) - t541 * t512 + t552 * t525;
t491 = t584 * t492 + t580 * t493;
t503 = -t541 * qJD(5) + t584 * t522 - t580 * t577;
t526 = t552 * mrSges(6,1) - t541 * mrSges(6,3);
t489 = m(6) * t491 - t521 * mrSges(6,2) + t503 * mrSges(6,3) + t540 * t512 - t552 * t526;
t480 = t584 * t488 + t580 * t489;
t538 = -t557 * mrSges(5,2) - t558 * mrSges(5,3);
t593 = -m(5) * t499 - t523 * mrSges(5,1) - t558 * t538 - t480;
t478 = m(4) * t500 - t523 * mrSges(4,3) - t558 * t537 + (t544 - t546) * t578 + (mrSges(4,1) - mrSges(5,2)) * t577 + t593;
t501 = t581 * t513 + t613 * t514;
t545 = t578 * mrSges(4,1) - t558 * mrSges(4,3);
t592 = -t576 * pkin(3) + t577 * qJ(4) - t557 * t536 + t501;
t498 = t578 * t614 - t592;
t547 = t558 * mrSges(5,1) + t578 * mrSges(5,2);
t495 = -t522 * pkin(4) - t553 * pkin(8) + ((2 * qJD(4)) + t548) * t578 + t592;
t594 = -m(6) * t495 + t503 * mrSges(6,1) - t504 * mrSges(6,2) + t540 * t525 - t541 * t526;
t591 = -m(5) * t498 + t577 * mrSges(5,3) + t578 * t547 - t594;
t485 = m(4) * t501 - t577 * mrSges(4,2) - t578 * t545 + (-t537 - t538) * t557 + (-mrSges(4,3) - mrSges(5,1)) * t522 + t591;
t474 = t613 * t478 + t581 * t485;
t542 = -t585 * g(3) - t607;
t565 = (-mrSges(3,1) * t585 + mrSges(3,2) * t582) * qJD(1);
t569 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t600;
t472 = m(3) * t542 + qJDD(2) * mrSges(3,1) - t566 * mrSges(3,3) + qJD(2) * t569 - t565 * t601 + t474;
t568 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t601;
t596 = -t581 * t478 + t613 * t485;
t473 = m(3) * t543 - qJDD(2) * mrSges(3,2) + t567 * mrSges(3,3) - qJD(2) * t568 + t565 * t600 + t596;
t597 = -t582 * t472 + t585 * t473;
t466 = m(2) * t572 - t587 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t597;
t559 = -t587 * pkin(6) + t595;
t497 = t522 * pkin(3) + t589;
t605 = -t580 * t488 + t584 * t489;
t479 = m(5) * t497 - t522 * mrSges(5,2) - t523 * mrSges(5,3) - t557 * t546 - t558 * t547 + t605;
t590 = m(4) * t524 + t522 * mrSges(4,1) + t523 * mrSges(4,2) + t557 * t544 + t558 * t545 + t479;
t588 = -m(3) * t559 + t567 * mrSges(3,1) - t566 * mrSges(3,2) - t568 * t601 + t569 * t600 - t590;
t476 = m(2) * t571 + qJDD(1) * mrSges(2,1) - t587 * mrSges(2,2) + t588;
t606 = t583 * t466 + t586 * t476;
t467 = t585 * t472 + t582 * t473;
t604 = t609 * t557 - t610 * t558 + t616 * t578;
t603 = t615 * t557 - t611 * t558 - t609 * t578;
t602 = -t611 * t557 + t617 * t558 + t610 * t578;
t598 = t586 * t466 - t583 * t476;
t556 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t582 + Ifges(3,4) * t585) * qJD(1);
t555 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t582 + Ifges(3,2) * t585) * qJD(1);
t554 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t582 + Ifges(3,6) * t585) * qJD(1);
t507 = Ifges(6,1) * t541 + Ifges(6,4) * t540 + Ifges(6,5) * t552;
t506 = Ifges(6,4) * t541 + Ifges(6,2) * t540 + Ifges(6,6) * t552;
t505 = Ifges(6,5) * t541 + Ifges(6,6) * t540 + Ifges(6,3) * t552;
t482 = mrSges(6,2) * t495 - mrSges(6,3) * t490 + Ifges(6,1) * t504 + Ifges(6,4) * t503 + Ifges(6,5) * t521 + t540 * t505 - t552 * t506;
t481 = -mrSges(6,1) * t495 + mrSges(6,3) * t491 + Ifges(6,4) * t504 + Ifges(6,2) * t503 + Ifges(6,6) * t521 - t541 * t505 + t552 * t507;
t468 = mrSges(5,1) * t499 + mrSges(6,1) * t490 + mrSges(4,2) * t524 - mrSges(6,2) * t491 - mrSges(4,3) * t500 - mrSges(5,3) * t497 + Ifges(6,5) * t504 + Ifges(6,6) * t503 + Ifges(6,3) * t521 + pkin(4) * t480 - qJ(4) * t479 + t541 * t506 - t540 * t507 + t603 * t578 + t610 * t577 + t604 * t557 + t617 * t523 - t611 * t522;
t463 = -mrSges(4,1) * t524 - mrSges(5,1) * t498 + mrSges(5,2) * t497 + mrSges(4,3) * t501 - pkin(3) * t479 - pkin(4) * t594 - pkin(8) * t605 - t584 * t481 - t580 * t482 - t615 * t522 + t611 * t523 + t604 * t558 + t609 * t577 + t602 * t578;
t462 = mrSges(3,2) * t559 - mrSges(3,3) * t542 + Ifges(3,1) * t566 + Ifges(3,4) * t567 + Ifges(3,5) * qJDD(2) - pkin(7) * t474 - qJD(2) * t555 - t581 * t463 + t613 * t468 + t554 * t600;
t461 = (pkin(3) * mrSges(5,2) + t616) * t577 - mrSges(3,1) * t542 + mrSges(3,2) * t543 - pkin(2) * t474 - Ifges(3,6) * t567 + mrSges(2,3) * t572 + Ifges(2,6) * qJDD(1) - pkin(1) * t467 - t584 * t482 + t587 * Ifges(2,5) - Ifges(3,3) * qJDD(2) + mrSges(2,1) * g(3) + t580 * t481 + pkin(8) * t480 + (qJ(4) * mrSges(5,1) + t609) * t522 - t610 * t523 + (qJ(4) * t538 - t602) * t557 + t603 * t558 - qJ(4) * t591 - pkin(3) * (-t578 * t546 + t593) - Ifges(3,5) * t566 + mrSges(5,3) * t498 - mrSges(5,2) * t499 - mrSges(4,1) * t500 + mrSges(4,2) * t501 + (-t582 * t555 + t585 * t556) * qJD(1);
t460 = -mrSges(3,1) * t559 + mrSges(3,3) * t543 + Ifges(3,4) * t566 + Ifges(3,2) * t567 + Ifges(3,6) * qJDD(2) - pkin(2) * t590 + pkin(7) * t596 + qJD(2) * t556 + t613 * t463 + t581 * t468 - t554 * t601;
t459 = -mrSges(2,2) * g(3) - mrSges(2,3) * t571 + Ifges(2,5) * qJDD(1) - t587 * Ifges(2,6) - pkin(6) * t467 - t582 * t460 + t585 * t462;
t1 = [-m(1) * g(1) + t598; -m(1) * g(2) + t606; (-m(1) - m(2)) * g(3) + t467; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t606 + t586 * t459 - t583 * t461; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t598 + t583 * t459 + t586 * t461; -mrSges(1,1) * g(2) + mrSges(2,1) * t571 + mrSges(1,2) * g(1) - mrSges(2,2) * t572 + Ifges(2,3) * qJDD(1) + pkin(1) * t588 + pkin(6) * t597 + t585 * t460 + t582 * t462;];
tauB = t1;
