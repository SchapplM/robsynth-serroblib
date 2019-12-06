% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:05:55
% EndTime: 2019-12-05 18:06:01
% DurationCPUTime: 3.40s
% Computational Cost: add. (28740->267), mult. (69192->337), div. (0->0), fcn. (45196->8), ass. (0->113)
t578 = sin(qJ(1));
t581 = cos(qJ(1));
t560 = t578 * g(2) - t581 * g(3);
t582 = qJD(1) ^ 2;
t623 = -t582 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t560;
t622 = Ifges(5,1) + Ifges(6,1);
t618 = Ifges(5,4) + Ifges(6,4);
t617 = Ifges(5,5) + Ifges(6,5);
t621 = Ifges(5,2) + Ifges(6,2);
t616 = -Ifges(5,6) - Ifges(6,6);
t620 = -Ifges(5,3) - Ifges(6,3);
t574 = sin(pkin(8));
t575 = cos(pkin(8));
t533 = -t575 * g(1) - t623 * t574;
t605 = t575 * qJD(1);
t564 = qJD(3) - t605;
t577 = sin(qJ(3));
t606 = t574 * qJD(1);
t598 = t577 * t606;
t546 = -t564 * mrSges(4,2) - mrSges(4,3) * t598;
t580 = cos(qJ(3));
t597 = t580 * t606;
t547 = t564 * mrSges(4,1) - mrSges(4,3) * t597;
t619 = -t546 * t577 - t547 * t580;
t615 = mrSges(3,2) * t574;
t612 = t574 ^ 2 * t582;
t534 = -t574 * g(1) + t623 * t575;
t555 = (-mrSges(3,1) * t575 + t615) * qJD(1);
t591 = -pkin(2) * t575 - pkin(6) * t574;
t557 = t591 * qJD(1);
t522 = t557 * t605 + t534;
t561 = t581 * g(2) + t578 * g(3);
t586 = -t582 * qJ(2) + qJDD(2) - t561;
t535 = (-pkin(1) + t591) * qJDD(1) + t586;
t532 = t580 * t535;
t603 = qJD(1) * qJD(3);
t551 = (qJDD(1) * t580 - t577 * t603) * t574;
t602 = t575 * qJDD(1);
t563 = qJDD(3) - t602;
t499 = t563 * pkin(3) - t551 * pkin(7) + t532 + (-pkin(3) * t580 * t612 - pkin(7) * t564 * t606 - t522) * t577;
t503 = t580 * t522 + t577 * t535;
t549 = t564 * pkin(3) - pkin(7) * t597;
t550 = (-qJDD(1) * t577 - t580 * t603) * t574;
t601 = t577 ^ 2 * t612;
t500 = -pkin(3) * t601 + t550 * pkin(7) - t564 * t549 + t503;
t576 = sin(qJ(4));
t579 = cos(qJ(4));
t492 = t579 * t499 - t576 * t500;
t541 = (-t576 * t580 - t577 * t579) * t606;
t510 = t541 * qJD(4) + t576 * t550 + t579 * t551;
t542 = (-t576 * t577 + t579 * t580) * t606;
t523 = -t541 * mrSges(6,1) + t542 * mrSges(6,2);
t524 = -t541 * mrSges(5,1) + t542 * mrSges(5,2);
t562 = qJD(4) + t564;
t527 = -t562 * mrSges(5,2) + t541 * mrSges(5,3);
t559 = qJDD(4) + t563;
t489 = -0.2e1 * qJD(5) * t542 + (t541 * t562 - t510) * qJ(5) + (t541 * t542 + t559) * pkin(4) + t492;
t526 = -t562 * mrSges(6,2) + t541 * mrSges(6,3);
t600 = m(6) * t489 + t559 * mrSges(6,1) + t562 * t526;
t481 = m(5) * t492 + t559 * mrSges(5,1) + t562 * t527 + (-t523 - t524) * t542 + (-mrSges(5,3) - mrSges(6,3)) * t510 + t600;
t493 = t576 * t499 + t579 * t500;
t509 = -t542 * qJD(4) + t579 * t550 - t576 * t551;
t529 = t562 * mrSges(6,1) - t542 * mrSges(6,3);
t530 = t562 * mrSges(5,1) - t542 * mrSges(5,3);
t528 = t562 * pkin(4) - t542 * qJ(5);
t540 = t541 ^ 2;
t491 = -t540 * pkin(4) + t509 * qJ(5) + 0.2e1 * qJD(5) * t541 - t562 * t528 + t493;
t599 = m(6) * t491 + t509 * mrSges(6,3) + t541 * t523;
t484 = m(5) * t493 + t509 * mrSges(5,3) + t541 * t524 + (-t529 - t530) * t562 + (-mrSges(5,2) - mrSges(6,2)) * t559 + t599;
t479 = t579 * t481 + t576 * t484;
t502 = -t577 * t522 + t532;
t548 = (mrSges(4,1) * t577 + mrSges(4,2) * t580) * t606;
t476 = m(4) * t502 + t563 * mrSges(4,1) - t551 * mrSges(4,3) + t564 * t546 - t548 * t597 + t479;
t592 = -t576 * t481 + t579 * t484;
t477 = m(4) * t503 - t563 * mrSges(4,2) + t550 * mrSges(4,3) - t564 * t547 - t548 * t598 + t592;
t593 = -t577 * t476 + t580 * t477;
t607 = qJDD(1) * mrSges(3,3);
t472 = m(3) * t534 + (qJD(1) * t555 + t607) * t575 + t593;
t521 = t557 * t606 - t533;
t501 = -t550 * pkin(3) - pkin(7) * t601 + t549 * t597 + t521;
t495 = -t509 * pkin(4) - t540 * qJ(5) + t542 * t528 + qJDD(5) + t501;
t588 = m(6) * t495 - t509 * mrSges(6,1) + t510 * mrSges(6,2) - t541 * t526 + t542 * t529;
t584 = m(5) * t501 - t509 * mrSges(5,1) + t510 * mrSges(5,2) - t541 * t527 + t542 * t530 + t588;
t583 = -m(4) * t521 + t550 * mrSges(4,1) - t551 * mrSges(4,2) - t584;
t486 = (-t607 + (-t555 + t619) * qJD(1)) * t574 + m(3) * t533 + t583;
t468 = t574 * t472 + t575 * t486;
t611 = t616 * t541 - t617 * t542 + t620 * t562;
t610 = -t621 * t541 - t618 * t542 + t616 * t562;
t609 = t618 * t541 + t622 * t542 + t617 * t562;
t594 = t575 * t472 - t574 * t486;
t467 = m(2) * t560 - t582 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t594;
t473 = t580 * t476 + t577 * t477;
t553 = -qJDD(1) * pkin(1) + t586;
t585 = -m(3) * t553 + mrSges(3,1) * t602 - t473 + (t575 ^ 2 * t582 + t612) * mrSges(3,3);
t469 = m(2) * t561 - t582 * mrSges(2,2) + (mrSges(2,1) - t615) * qJDD(1) + t585;
t595 = t581 * t467 - t578 * t469;
t590 = Ifges(3,1) * t574 + Ifges(3,4) * t575;
t589 = Ifges(3,5) * t574 + Ifges(3,6) * t575;
t587 = -t578 * t467 - t581 * t469;
t556 = t589 * qJD(1);
t538 = Ifges(4,5) * t564 + (Ifges(4,1) * t580 - Ifges(4,4) * t577) * t606;
t537 = Ifges(4,6) * t564 + (Ifges(4,4) * t580 - Ifges(4,2) * t577) * t606;
t536 = Ifges(4,3) * t564 + (Ifges(4,5) * t580 - Ifges(4,6) * t577) * t606;
t487 = -t510 * mrSges(6,3) - t542 * t523 + t600;
t478 = mrSges(5,2) * t501 + mrSges(6,2) * t495 - mrSges(5,3) * t492 - mrSges(6,3) * t489 - qJ(5) * t487 + t618 * t509 + t622 * t510 - t611 * t541 + t617 * t559 + t610 * t562;
t474 = -mrSges(5,1) * t501 + mrSges(5,3) * t493 - mrSges(6,1) * t495 + mrSges(6,3) * t491 - pkin(4) * t588 + qJ(5) * t599 + (-qJ(5) * t529 + t609) * t562 + (-qJ(5) * mrSges(6,2) - t616) * t559 + t611 * t542 + t618 * t510 + t621 * t509;
t465 = mrSges(4,2) * t521 - mrSges(4,3) * t502 + Ifges(4,1) * t551 + Ifges(4,4) * t550 + Ifges(4,5) * t563 - pkin(7) * t479 - t576 * t474 + t579 * t478 - t536 * t598 - t564 * t537;
t464 = -mrSges(4,1) * t521 + mrSges(4,3) * t503 + Ifges(4,4) * t551 + Ifges(4,2) * t550 + Ifges(4,6) * t563 - pkin(3) * t584 + pkin(7) * t592 + t579 * t474 + t576 * t478 - t536 * t597 + t564 * t538;
t463 = Ifges(3,2) * t602 - mrSges(3,1) * t553 - mrSges(4,1) * t502 - mrSges(5,1) * t492 - mrSges(6,1) * t489 + mrSges(4,2) * t503 + mrSges(5,2) * t493 + mrSges(6,2) * t491 + mrSges(3,3) * t534 - Ifges(4,5) * t551 - Ifges(4,6) * t550 - Ifges(4,3) * t563 - pkin(2) * t473 - pkin(3) * t479 - pkin(4) * t487 + t620 * t559 + t610 * t542 + t609 * t541 - t617 * t510 + t616 * t509 + (Ifges(3,4) * qJDD(1) + (-t537 * t580 - t538 * t577 - t556) * qJD(1)) * t574;
t462 = mrSges(3,2) * t553 - mrSges(3,3) * t533 - pkin(6) * t473 + t590 * qJDD(1) - t577 * t464 + t580 * t465 + t556 * t605;
t461 = t582 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t560 - mrSges(3,1) * t533 + mrSges(3,2) * t534 - t577 * t465 - t580 * t464 - pkin(2) * t583 - pkin(6) * t593 - pkin(1) * t468 + (Ifges(2,6) - t589) * qJDD(1) + (-pkin(2) * t619 * t574 + (-t574 * (Ifges(3,4) * t574 + Ifges(3,2) * t575) + t575 * t590) * qJD(1)) * qJD(1);
t460 = -mrSges(2,2) * g(1) - mrSges(2,3) * t561 + Ifges(2,5) * qJDD(1) - t582 * Ifges(2,6) - qJ(2) * t468 + t575 * t462 - t574 * t463;
t1 = [(-m(1) - m(2)) * g(1) + t468; -m(1) * g(2) + t587; -m(1) * g(3) + t595; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t561 - mrSges(2,2) * t560 + t574 * t462 + t575 * t463 + pkin(1) * (-qJDD(1) * t615 + t585) + qJ(2) * t594; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t595 - t578 * t460 - t581 * t461; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t587 + t581 * t460 - t578 * t461;];
tauB = t1;
