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
% m [6x1]
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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:32:06
% EndTime: 2022-01-23 09:32:11
% DurationCPUTime: 3.52s
% Computational Cost: add. (28740->267), mult. (69192->337), div. (0->0), fcn. (45196->8), ass. (0->113)
t580 = sin(qJ(1));
t583 = cos(qJ(1));
t563 = -g(1) * t583 - g(2) * t580;
t584 = qJD(1) ^ 2;
t625 = -pkin(1) * t584 + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t563;
t624 = Ifges(5,1) + Ifges(6,1);
t620 = Ifges(5,4) + Ifges(6,4);
t619 = Ifges(5,5) + Ifges(6,5);
t623 = Ifges(5,2) + Ifges(6,2);
t618 = -Ifges(5,6) - Ifges(6,6);
t622 = -Ifges(5,3) - Ifges(6,3);
t576 = sin(pkin(8));
t577 = cos(pkin(8));
t535 = -t577 * g(3) - t576 * t625;
t607 = qJD(1) * t577;
t566 = qJD(3) - t607;
t579 = sin(qJ(3));
t608 = qJD(1) * t576;
t599 = t579 * t608;
t548 = -mrSges(4,2) * t566 - mrSges(4,3) * t599;
t582 = cos(qJ(3));
t598 = t582 * t608;
t549 = mrSges(4,1) * t566 - mrSges(4,3) * t598;
t621 = -t548 * t579 - t549 * t582;
t617 = mrSges(3,2) * t576;
t614 = t576 ^ 2 * t584;
t536 = -g(3) * t576 + t577 * t625;
t557 = (-mrSges(3,1) * t577 + t617) * qJD(1);
t592 = -pkin(2) * t577 - pkin(6) * t576;
t559 = t592 * qJD(1);
t524 = t559 * t607 + t536;
t562 = g(1) * t580 - t583 * g(2);
t588 = -qJ(2) * t584 + qJDD(2) - t562;
t537 = (-pkin(1) + t592) * qJDD(1) + t588;
t534 = t582 * t537;
t604 = qJD(1) * qJD(3);
t553 = (qJDD(1) * t582 - t579 * t604) * t576;
t603 = qJDD(1) * t577;
t565 = qJDD(3) - t603;
t501 = pkin(3) * t565 - pkin(7) * t553 + t534 + (-pkin(3) * t582 * t614 - pkin(7) * t566 * t608 - t524) * t579;
t505 = t524 * t582 + t579 * t537;
t551 = pkin(3) * t566 - pkin(7) * t598;
t552 = (-qJDD(1) * t579 - t582 * t604) * t576;
t602 = t579 ^ 2 * t614;
t502 = -pkin(3) * t602 + pkin(7) * t552 - t551 * t566 + t505;
t578 = sin(qJ(4));
t581 = cos(qJ(4));
t494 = t501 * t581 - t502 * t578;
t543 = (-t578 * t582 - t579 * t581) * t608;
t512 = qJD(4) * t543 + t552 * t578 + t553 * t581;
t544 = (-t578 * t579 + t581 * t582) * t608;
t525 = -mrSges(6,1) * t543 + mrSges(6,2) * t544;
t526 = -mrSges(5,1) * t543 + mrSges(5,2) * t544;
t564 = qJD(4) + t566;
t529 = -mrSges(5,2) * t564 + mrSges(5,3) * t543;
t561 = qJDD(4) + t565;
t491 = -0.2e1 * qJD(5) * t544 + (t543 * t564 - t512) * qJ(5) + (t543 * t544 + t561) * pkin(4) + t494;
t528 = -mrSges(6,2) * t564 + mrSges(6,3) * t543;
t601 = m(6) * t491 + t561 * mrSges(6,1) + t564 * t528;
t483 = m(5) * t494 + mrSges(5,1) * t561 + t529 * t564 + (-t525 - t526) * t544 + (-mrSges(5,3) - mrSges(6,3)) * t512 + t601;
t495 = t501 * t578 + t502 * t581;
t511 = -qJD(4) * t544 + t552 * t581 - t553 * t578;
t531 = mrSges(6,1) * t564 - mrSges(6,3) * t544;
t532 = mrSges(5,1) * t564 - mrSges(5,3) * t544;
t530 = pkin(4) * t564 - qJ(5) * t544;
t542 = t543 ^ 2;
t493 = -pkin(4) * t542 + qJ(5) * t511 + 0.2e1 * qJD(5) * t543 - t530 * t564 + t495;
t600 = m(6) * t493 + mrSges(6,3) * t511 + t525 * t543;
t486 = m(5) * t495 + mrSges(5,3) * t511 + t526 * t543 + (-t531 - t532) * t564 + (-mrSges(5,2) - mrSges(6,2)) * t561 + t600;
t481 = t483 * t581 + t486 * t578;
t504 = -t524 * t579 + t534;
t550 = (mrSges(4,1) * t579 + mrSges(4,2) * t582) * t608;
t478 = m(4) * t504 + mrSges(4,1) * t565 - mrSges(4,3) * t553 + t548 * t566 - t550 * t598 + t481;
t593 = -t483 * t578 + t486 * t581;
t479 = m(4) * t505 - mrSges(4,2) * t565 + mrSges(4,3) * t552 - t549 * t566 - t550 * t599 + t593;
t594 = -t478 * t579 + t479 * t582;
t606 = qJDD(1) * mrSges(3,3);
t474 = m(3) * t536 + (qJD(1) * t557 + t606) * t577 + t594;
t523 = t559 * t608 - t535;
t503 = -pkin(3) * t552 - pkin(7) * t602 + t551 * t598 + t523;
t497 = -pkin(4) * t511 - qJ(5) * t542 + t530 * t544 + qJDD(5) + t503;
t589 = m(6) * t497 - mrSges(6,1) * t511 + mrSges(6,2) * t512 - t528 * t543 + t531 * t544;
t586 = m(5) * t503 - mrSges(5,1) * t511 + mrSges(5,2) * t512 - t529 * t543 + t532 * t544 + t589;
t585 = -m(4) * t523 + t552 * mrSges(4,1) - mrSges(4,2) * t553 - t586;
t488 = m(3) * t535 + (-t606 + (-t557 + t621) * qJD(1)) * t576 + t585;
t595 = t474 * t577 - t488 * t576;
t468 = m(2) * t563 - mrSges(2,1) * t584 - qJDD(1) * mrSges(2,2) + t595;
t475 = t478 * t582 + t479 * t579;
t555 = -qJDD(1) * pkin(1) + t588;
t587 = -m(3) * t555 + mrSges(3,1) * t603 - t475 + (t577 ^ 2 * t584 + t614) * mrSges(3,3);
t471 = m(2) * t562 - mrSges(2,2) * t584 + (mrSges(2,1) - t617) * qJDD(1) + t587;
t613 = t468 * t580 + t471 * t583;
t469 = t474 * t576 + t488 * t577;
t612 = t543 * t618 - t544 * t619 + t564 * t622;
t611 = -t543 * t623 - t544 * t620 + t564 * t618;
t610 = t543 * t620 + t544 * t624 + t564 * t619;
t596 = t468 * t583 - t471 * t580;
t591 = Ifges(3,1) * t576 + Ifges(3,4) * t577;
t590 = Ifges(3,5) * t576 + Ifges(3,6) * t577;
t558 = t590 * qJD(1);
t540 = Ifges(4,5) * t566 + (Ifges(4,1) * t582 - Ifges(4,4) * t579) * t608;
t539 = Ifges(4,6) * t566 + (Ifges(4,4) * t582 - Ifges(4,2) * t579) * t608;
t538 = Ifges(4,3) * t566 + (Ifges(4,5) * t582 - Ifges(4,6) * t579) * t608;
t489 = -mrSges(6,3) * t512 - t525 * t544 + t601;
t480 = mrSges(5,2) * t503 + mrSges(6,2) * t497 - mrSges(5,3) * t494 - mrSges(6,3) * t491 - qJ(5) * t489 + t620 * t511 + t512 * t624 - t612 * t543 + t619 * t561 + t611 * t564;
t476 = -mrSges(5,1) * t503 + mrSges(5,3) * t495 - mrSges(6,1) * t497 + mrSges(6,3) * t493 - pkin(4) * t589 + qJ(5) * t600 + (-qJ(5) * t531 + t610) * t564 + (-mrSges(6,2) * qJ(5) - t618) * t561 + t612 * t544 + t620 * t512 + t623 * t511;
t465 = mrSges(4,2) * t523 - mrSges(4,3) * t504 + Ifges(4,1) * t553 + Ifges(4,4) * t552 + Ifges(4,5) * t565 - pkin(7) * t481 - t476 * t578 + t480 * t581 - t538 * t599 - t539 * t566;
t464 = -mrSges(4,1) * t523 + mrSges(4,3) * t505 + Ifges(4,4) * t553 + Ifges(4,2) * t552 + Ifges(4,6) * t565 - pkin(3) * t586 + pkin(7) * t593 + t581 * t476 + t578 * t480 - t538 * t598 + t566 * t540;
t463 = Ifges(3,2) * t603 - mrSges(3,1) * t555 - mrSges(4,1) * t504 - mrSges(5,1) * t494 - mrSges(6,1) * t491 + mrSges(4,2) * t505 + mrSges(5,2) * t495 + mrSges(6,2) * t493 + mrSges(3,3) * t536 - Ifges(4,5) * t553 - Ifges(4,6) * t552 - Ifges(4,3) * t565 - pkin(2) * t475 - pkin(3) * t481 - pkin(4) * t489 + t622 * t561 + t611 * t544 + t610 * t543 - t619 * t512 + t618 * t511 + (Ifges(3,4) * qJDD(1) + (-t539 * t582 - t540 * t579 - t558) * qJD(1)) * t576;
t462 = mrSges(3,2) * t555 - mrSges(3,3) * t535 - pkin(6) * t475 + qJDD(1) * t591 - t464 * t579 + t465 * t582 + t558 * t607;
t461 = t584 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t563 - mrSges(3,1) * t535 + mrSges(3,2) * t536 - t579 * t465 - t582 * t464 - pkin(2) * t585 - pkin(6) * t594 - pkin(1) * t469 + (Ifges(2,6) - t590) * qJDD(1) + (-pkin(2) * t621 * t576 + (-t576 * (Ifges(3,4) * t576 + Ifges(3,2) * t577) + t577 * t591) * qJD(1)) * qJD(1);
t460 = -mrSges(2,2) * g(3) - mrSges(2,3) * t562 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t584 - qJ(2) * t469 + t462 * t577 - t463 * t576;
t1 = [-m(1) * g(1) + t596; -m(1) * g(2) + t613; (-m(1) - m(2)) * g(3) + t469; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t613 + t460 * t583 - t461 * t580; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t596 + t580 * t460 + t583 * t461; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + mrSges(2,1) * t562 - mrSges(2,2) * t563 + t576 * t462 + t577 * t463 + pkin(1) * (-qJDD(1) * t617 + t587) + qJ(2) * t595;];
tauB = t1;
