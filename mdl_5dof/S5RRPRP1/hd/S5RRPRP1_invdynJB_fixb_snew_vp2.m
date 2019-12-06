% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5RRPRP1_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:12
% EndTime: 2019-12-05 18:22:14
% DurationCPUTime: 1.80s
% Computational Cost: add. (23252->208), mult. (29959->247), div. (0->0), fcn. (14555->8), ass. (0->93)
t627 = Ifges(5,1) + Ifges(6,1);
t619 = Ifges(5,4) + Ifges(6,4);
t618 = Ifges(5,5) + Ifges(6,5);
t626 = Ifges(5,2) + Ifges(6,2);
t617 = Ifges(5,6) + Ifges(6,6);
t625 = Ifges(5,3) + Ifges(6,3);
t580 = qJD(1) + qJD(2);
t587 = sin(qJ(4));
t590 = cos(qJ(4));
t552 = (-mrSges(6,1) * t590 + mrSges(6,2) * t587) * t580;
t579 = qJDD(1) + qJDD(2);
t609 = qJD(4) * t580;
t605 = t590 * t609;
t554 = t587 * t579 + t605;
t589 = sin(qJ(1));
t592 = cos(qJ(1));
t569 = t592 * g(2) + t589 * g(3);
t560 = qJDD(1) * pkin(1) + t569;
t568 = t589 * g(2) - t592 * g(3);
t593 = qJD(1) ^ 2;
t561 = -t593 * pkin(1) + t568;
t588 = sin(qJ(2));
t591 = cos(qJ(2));
t538 = t591 * t560 - t588 * t561;
t535 = t579 * pkin(2) + t538;
t539 = t588 * t560 + t591 * t561;
t578 = t580 ^ 2;
t536 = -t578 * pkin(2) + t539;
t585 = sin(pkin(8));
t586 = cos(pkin(8));
t531 = t585 * t535 + t586 * t536;
t528 = -t578 * pkin(3) + t579 * pkin(7) + t531;
t584 = -g(1) + qJDD(3);
t571 = t590 * t584;
t608 = qJD(5) * t580;
t621 = pkin(4) * t578;
t521 = qJDD(4) * pkin(4) + t571 + (-t554 + t605) * qJ(5) + (t590 * t621 - t528 - 0.2e1 * t608) * t587;
t614 = t580 * t590;
t565 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t614;
t607 = m(6) * t521 + qJDD(4) * mrSges(6,1) + qJD(4) * t565;
t615 = t580 * t587;
t518 = -t554 * mrSges(6,3) - t552 * t615 + t607;
t525 = t590 * t528 + t587 * t584;
t555 = t590 * t579 - t587 * t609;
t562 = qJD(4) * pkin(4) - qJ(5) * t615;
t583 = t590 ^ 2;
t522 = t555 * qJ(5) - qJD(4) * t562 - t583 * t621 + 0.2e1 * t590 * t608 + t525;
t524 = -t587 * t528 + t571;
t611 = (t627 * t587 + t619 * t590) * t580 + t618 * qJD(4);
t612 = (t619 * t587 + t626 * t590) * t580 + t617 * qJD(4);
t624 = mrSges(5,1) * t524 + mrSges(6,1) * t521 - mrSges(5,2) * t525 - mrSges(6,2) * t522 + pkin(4) * t518 + t625 * qJDD(4) + t618 * t554 + t617 * t555 + (t612 * t587 - t611 * t590) * t580;
t620 = -mrSges(5,2) - mrSges(6,2);
t553 = (-mrSges(5,1) * t590 + mrSges(5,2) * t587) * t580;
t566 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t614;
t515 = m(5) * t524 + qJDD(4) * mrSges(5,1) + qJD(4) * t566 + (-t552 - t553) * t615 + (-mrSges(5,3) - mrSges(6,3)) * t554 + t607;
t606 = m(6) * t522 + t555 * mrSges(6,3) + t552 * t614;
t563 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t615;
t610 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t615 - t563;
t516 = m(5) * t525 + t555 * mrSges(5,3) + t610 * qJD(4) + t620 * qJDD(4) + t553 * t614 + t606;
t601 = -t587 * t515 + t590 * t516;
t505 = m(4) * t531 - t578 * mrSges(4,1) - t579 * mrSges(4,2) + t601;
t530 = t586 * t535 - t585 * t536;
t598 = -t579 * pkin(3) - t530;
t527 = -t578 * pkin(7) + t598;
t523 = t562 * t615 - t555 * pkin(4) + qJDD(5) + (-qJ(5) * t583 - pkin(7)) * t578 + t598;
t600 = -m(6) * t523 + t555 * mrSges(6,1) + t565 * t614;
t595 = -m(5) * t527 + t555 * mrSges(5,1) + t620 * t554 + t566 * t614 + t610 * t615 + t600;
t510 = m(4) * t530 + t579 * mrSges(4,1) - t578 * mrSges(4,2) + t595;
t498 = t585 * t505 + t586 * t510;
t495 = m(3) * t538 + t579 * mrSges(3,1) - t578 * mrSges(3,2) + t498;
t602 = t586 * t505 - t585 * t510;
t496 = m(3) * t539 - t578 * mrSges(3,1) - t579 * mrSges(3,2) + t602;
t490 = t591 * t495 + t588 * t496;
t508 = t590 * t515 + t587 * t516;
t613 = (t618 * t587 + t617 * t590) * t580 + t625 * qJD(4);
t506 = m(4) * t584 + t508;
t603 = -t588 * t495 + t591 * t496;
t487 = m(2) * t568 - t593 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t603;
t488 = m(2) * t569 + qJDD(1) * mrSges(2,1) - t593 * mrSges(2,2) + t490;
t604 = t592 * t487 - t589 * t488;
t599 = -t589 * t487 - t592 * t488;
t517 = t554 * mrSges(6,2) + t563 * t615 - t600;
t500 = -mrSges(5,1) * t527 + mrSges(5,3) * t525 - mrSges(6,1) * t523 + mrSges(6,3) * t522 - pkin(4) * t517 + qJ(5) * t606 - t613 * t615 + t626 * t555 + t619 * t554 + (-qJ(5) * mrSges(6,2) + t617) * qJDD(4) + (-qJ(5) * t563 + t611) * qJD(4);
t502 = mrSges(5,2) * t527 + mrSges(6,2) * t523 - mrSges(5,3) * t524 - mrSges(6,3) * t521 - qJ(5) * t518 - t612 * qJD(4) + t618 * qJDD(4) + t627 * t554 + t619 * t555 + t613 * t614;
t596 = mrSges(3,1) * t538 + mrSges(4,1) * t530 - mrSges(3,2) * t539 - mrSges(4,2) * t531 + pkin(2) * t498 + pkin(3) * t595 + pkin(7) * t601 + t590 * t500 + t587 * t502 + (Ifges(4,3) + Ifges(3,3)) * t579;
t594 = mrSges(2,1) * t569 - mrSges(2,2) * t568 + Ifges(2,3) * qJDD(1) + pkin(1) * t490 + t596;
t491 = -mrSges(4,1) * t584 + mrSges(4,3) * t531 + t578 * Ifges(4,5) + Ifges(4,6) * t579 - pkin(3) * t508 - t624;
t485 = mrSges(4,2) * t584 - mrSges(4,3) * t530 + Ifges(4,5) * t579 - t578 * Ifges(4,6) - pkin(7) * t508 - t587 * t500 + t590 * t502;
t484 = -mrSges(3,2) * g(1) - mrSges(3,3) * t538 + Ifges(3,5) * t579 - t578 * Ifges(3,6) - qJ(3) * t498 + t586 * t485 - t585 * t491;
t483 = mrSges(3,1) * g(1) + mrSges(3,3) * t539 + t578 * Ifges(3,5) + Ifges(3,6) * t579 - pkin(2) * t506 + qJ(3) * t602 + t585 * t485 + t586 * t491;
t482 = -mrSges(2,2) * g(1) - mrSges(2,3) * t569 + Ifges(2,5) * qJDD(1) - t593 * Ifges(2,6) - pkin(6) * t490 - t588 * t483 + t591 * t484;
t481 = Ifges(2,6) * qJDD(1) + t593 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t568 + t588 * t484 + t591 * t483 - pkin(1) * (-m(3) * g(1) + t506) + pkin(6) * t603;
t1 = [(-m(1) - m(2) - m(3)) * g(1) + t506; -m(1) * g(2) + t599; -m(1) * g(3) + t604; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t604 - t592 * t481 - t589 * t482; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t599 - t589 * t481 + t592 * t482; t594; t596; t506; t624; t517;];
tauJB = t1;
