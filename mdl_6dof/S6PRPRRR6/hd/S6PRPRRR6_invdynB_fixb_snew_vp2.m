% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PRPRRR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PRPRRR6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:31:07
% EndTime: 2019-05-05 01:31:17
% DurationCPUTime: 7.56s
% Computational Cost: add. (102822->292), mult. (192779->366), div. (0->0), fcn. (128458->12), ass. (0->126)
t595 = sin(pkin(11));
t597 = cos(pkin(11));
t579 = g(1) * t595 - g(2) * t597;
t580 = -g(1) * t597 - g(2) * t595;
t592 = -g(3) + qJDD(1);
t606 = cos(qJ(2));
t598 = cos(pkin(6));
t602 = sin(qJ(2));
t627 = t598 * t602;
t596 = sin(pkin(6));
t628 = t596 * t602;
t542 = t579 * t627 + t606 * t580 + t592 * t628;
t634 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t542;
t541 = -t602 * t580 + (t579 * t598 + t592 * t596) * t606;
t633 = -pkin(2) - pkin(8);
t632 = mrSges(3,1) - mrSges(4,2);
t631 = (-Ifges(4,4) + Ifges(3,5));
t630 = Ifges(4,5) - Ifges(3,6);
t608 = qJD(2) ^ 2;
t611 = -t608 * qJ(3) + qJDD(3) - t541;
t536 = t633 * qJDD(2) + t611;
t559 = -t579 * t596 + t592 * t598;
t601 = sin(qJ(4));
t605 = cos(qJ(4));
t530 = t601 * t536 + t605 * t559;
t575 = (mrSges(5,1) * t601 + mrSges(5,2) * t605) * qJD(2);
t624 = qJD(2) * qJD(4);
t586 = t605 * t624;
t577 = -t601 * qJDD(2) - t586;
t625 = qJD(2) * t605;
t582 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t625;
t588 = t601 * qJD(2);
t576 = (pkin(4) * t601 - pkin(9) * t605) * qJD(2);
t607 = qJD(4) ^ 2;
t522 = -pkin(4) * t607 + qJDD(4) * pkin(9) - t576 * t588 + t530;
t535 = t633 * t608 - t634;
t622 = t601 * t624;
t578 = qJDD(2) * t605 - t622;
t525 = (-t578 + t622) * pkin(9) + (-t577 + t586) * pkin(4) + t535;
t600 = sin(qJ(5));
t604 = cos(qJ(5));
t514 = -t600 * t522 + t604 * t525;
t573 = qJD(4) * t604 - t600 * t625;
t549 = qJD(5) * t573 + qJDD(4) * t600 + t578 * t604;
t570 = qJDD(5) - t577;
t574 = qJD(4) * t600 + t604 * t625;
t585 = t588 + qJD(5);
t512 = (t573 * t585 - t549) * pkin(10) + (t573 * t574 + t570) * pkin(5) + t514;
t515 = t604 * t522 + t600 * t525;
t548 = -qJD(5) * t574 + qJDD(4) * t604 - t578 * t600;
t557 = pkin(5) * t585 - pkin(10) * t574;
t569 = t573 ^ 2;
t513 = -pkin(5) * t569 + pkin(10) * t548 - t557 * t585 + t515;
t599 = sin(qJ(6));
t603 = cos(qJ(6));
t510 = t512 * t603 - t513 * t599;
t550 = t573 * t603 - t574 * t599;
t519 = qJD(6) * t550 + t548 * t599 + t549 * t603;
t551 = t573 * t599 + t574 * t603;
t532 = -mrSges(7,1) * t550 + mrSges(7,2) * t551;
t584 = qJD(6) + t585;
t539 = -mrSges(7,2) * t584 + mrSges(7,3) * t550;
t565 = qJDD(6) + t570;
t508 = m(7) * t510 + mrSges(7,1) * t565 - mrSges(7,3) * t519 - t532 * t551 + t539 * t584;
t511 = t512 * t599 + t513 * t603;
t518 = -qJD(6) * t551 + t548 * t603 - t549 * t599;
t540 = mrSges(7,1) * t584 - mrSges(7,3) * t551;
t509 = m(7) * t511 - mrSges(7,2) * t565 + mrSges(7,3) * t518 + t532 * t550 - t540 * t584;
t501 = t603 * t508 + t599 * t509;
t552 = -mrSges(6,1) * t573 + mrSges(6,2) * t574;
t555 = -mrSges(6,2) * t585 + mrSges(6,3) * t573;
t499 = m(6) * t514 + mrSges(6,1) * t570 - mrSges(6,3) * t549 - t552 * t574 + t555 * t585 + t501;
t556 = mrSges(6,1) * t585 - mrSges(6,3) * t574;
t618 = -t508 * t599 + t603 * t509;
t500 = m(6) * t515 - mrSges(6,2) * t570 + mrSges(6,3) * t548 + t552 * t573 - t556 * t585 + t618;
t619 = -t499 * t600 + t604 * t500;
t494 = m(5) * t530 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t577 - qJD(4) * t582 - t575 * t588 + t619;
t529 = t536 * t605 - t601 * t559;
t581 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t588;
t521 = -qJDD(4) * pkin(4) - pkin(9) * t607 + t576 * t625 - t529;
t516 = -pkin(5) * t548 - pkin(10) * t569 + t557 * t574 + t521;
t612 = m(7) * t516 - t518 * mrSges(7,1) + mrSges(7,2) * t519 - t550 * t539 + t540 * t551;
t609 = -m(6) * t521 + t548 * mrSges(6,1) - mrSges(6,2) * t549 + t573 * t555 - t556 * t574 - t612;
t504 = m(5) * t529 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t578 + qJD(4) * t581 - t575 * t625 + t609;
t487 = t601 * t494 + t605 * t504;
t538 = -qJDD(2) * pkin(2) + t611;
t614 = -m(4) * t538 + (t608 * mrSges(4,3)) - t487;
t482 = m(3) * t541 - (t608 * mrSges(3,2)) + t632 * qJDD(2) + t614;
t629 = t482 * t606;
t620 = t605 * t494 - t504 * t601;
t485 = m(4) * t559 + t620;
t484 = m(3) * t559 + t485;
t537 = t608 * pkin(2) + t634;
t495 = t604 * t499 + t600 * t500;
t613 = -m(5) * t535 + mrSges(5,1) * t577 - t578 * mrSges(5,2) - t581 * t588 - t582 * t625 - t495;
t610 = -m(4) * t537 + (t608 * mrSges(4,2)) + qJDD(2) * mrSges(4,3) - t613;
t492 = m(3) * t542 - (mrSges(3,1) * t608) - qJDD(2) * mrSges(3,2) + t610;
t473 = -t484 * t596 + t492 * t627 + t598 * t629;
t471 = m(2) * t579 + t473;
t478 = -t482 * t602 + t606 * t492;
t477 = m(2) * t580 + t478;
t626 = t597 * t471 + t595 * t477;
t472 = t598 * t484 + t492 * t628 + t596 * t629;
t621 = -t471 * t595 + t597 * t477;
t526 = Ifges(7,5) * t551 + Ifges(7,6) * t550 + Ifges(7,3) * t584;
t528 = Ifges(7,1) * t551 + Ifges(7,4) * t550 + Ifges(7,5) * t584;
t502 = -mrSges(7,1) * t516 + mrSges(7,3) * t511 + Ifges(7,4) * t519 + Ifges(7,2) * t518 + Ifges(7,6) * t565 - t526 * t551 + t528 * t584;
t527 = Ifges(7,4) * t551 + Ifges(7,2) * t550 + Ifges(7,6) * t584;
t503 = mrSges(7,2) * t516 - mrSges(7,3) * t510 + Ifges(7,1) * t519 + Ifges(7,4) * t518 + Ifges(7,5) * t565 + t526 * t550 - t527 * t584;
t543 = Ifges(6,5) * t574 + Ifges(6,6) * t573 + Ifges(6,3) * t585;
t545 = Ifges(6,1) * t574 + Ifges(6,4) * t573 + Ifges(6,5) * t585;
t486 = -mrSges(6,1) * t521 + mrSges(6,3) * t515 + Ifges(6,4) * t549 + Ifges(6,2) * t548 + Ifges(6,6) * t570 - pkin(5) * t612 + pkin(10) * t618 + t603 * t502 + t599 * t503 - t574 * t543 + t585 * t545;
t544 = Ifges(6,4) * t574 + Ifges(6,2) * t573 + Ifges(6,6) * t585;
t488 = mrSges(6,2) * t521 - mrSges(6,3) * t514 + Ifges(6,1) * t549 + Ifges(6,4) * t548 + Ifges(6,5) * t570 - pkin(10) * t501 - t502 * t599 + t503 * t603 + t543 * t573 - t544 * t585;
t562 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t605 - Ifges(5,6) * t601) * qJD(2);
t563 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t605 - Ifges(5,2) * t601) * qJD(2);
t474 = mrSges(5,2) * t535 - mrSges(5,3) * t529 + Ifges(5,1) * t578 + Ifges(5,4) * t577 + Ifges(5,5) * qJDD(4) - pkin(9) * t495 - qJD(4) * t563 - t486 * t600 + t488 * t604 - t562 * t588;
t564 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t605 - Ifges(5,4) * t601) * qJD(2);
t479 = Ifges(5,4) * t578 + Ifges(5,2) * t577 + Ifges(5,6) * qJDD(4) - t562 * t625 + qJD(4) * t564 - mrSges(5,1) * t535 + mrSges(5,3) * t530 - Ifges(6,5) * t549 - Ifges(6,6) * t548 - Ifges(6,3) * t570 - t574 * t544 + t573 * t545 - mrSges(6,1) * t514 + mrSges(6,2) * t515 - Ifges(7,5) * t519 - Ifges(7,6) * t518 - Ifges(7,3) * t565 - t551 * t527 + t550 * t528 - mrSges(7,1) * t510 + mrSges(7,2) * t511 - pkin(5) * t501 - pkin(4) * t495;
t468 = -mrSges(4,1) * t537 + mrSges(3,3) * t542 - pkin(2) * t485 - pkin(3) * t613 - pkin(8) * t620 - t630 * qJDD(2) - t601 * t474 - t605 * t479 - t632 * t559 + (t631 * t608);
t469 = -qJ(3) * t485 - mrSges(3,3) * t541 + pkin(3) * t487 + mrSges(4,1) * t538 + pkin(9) * t619 + t600 * t488 + t604 * t486 + pkin(4) * t609 + mrSges(5,1) * t529 - mrSges(5,2) * t530 + Ifges(5,5) * t578 + Ifges(5,6) * t577 + Ifges(5,3) * qJDD(4) + t630 * t608 + (mrSges(3,2) - mrSges(4,3)) * t559 + t631 * qJDD(2) + (t563 * t605 + t564 * t601) * qJD(2);
t615 = pkin(7) * t478 + t468 * t606 + t469 * t602;
t467 = mrSges(3,1) * t541 - mrSges(3,2) * t542 + mrSges(4,2) * t538 - mrSges(4,3) * t537 + t605 * t474 - t601 * t479 - pkin(8) * t487 + pkin(2) * t614 + qJ(3) * t610 + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * qJDD(2);
t466 = mrSges(2,2) * t592 - mrSges(2,3) * t579 - t602 * t468 + t606 * t469 + (-t472 * t596 - t473 * t598) * pkin(7);
t465 = -mrSges(2,1) * t592 + mrSges(2,3) * t580 - pkin(1) * t472 - t596 * t467 + t615 * t598;
t1 = [-m(1) * g(1) + t621; -m(1) * g(2) + t626; -m(1) * g(3) + m(2) * t592 + t472; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t626 - t595 * t465 + t597 * t466; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t621 + t597 * t465 + t595 * t466; -mrSges(1,1) * g(2) + mrSges(2,1) * t579 + mrSges(1,2) * g(1) - mrSges(2,2) * t580 + pkin(1) * t473 + t598 * t467 + t615 * t596;];
tauB  = t1;
