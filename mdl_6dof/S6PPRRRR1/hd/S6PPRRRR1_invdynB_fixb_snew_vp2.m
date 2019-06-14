% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6PPRRRR1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:39:43
% EndTime: 2019-05-04 20:40:04
% DurationCPUTime: 20.30s
% Computational Cost: add. (377849->295), mult. (696817->391), div. (0->0), fcn. (547544->16), ass. (0->134)
t597 = sin(pkin(12));
t601 = cos(pkin(12));
t587 = -g(1) * t601 - g(2) * t597;
t596 = sin(pkin(13));
t600 = cos(pkin(13));
t586 = g(1) * t597 - g(2) * t601;
t595 = -g(3) + qJDD(1);
t599 = sin(pkin(6));
t603 = cos(pkin(6));
t620 = t586 * t603 + t595 * t599;
t560 = -t596 * t587 + t600 * t620;
t561 = t600 * t587 + t596 * t620;
t573 = -t586 * t599 + t595 * t603 + qJDD(2);
t611 = cos(qJ(3));
t602 = cos(pkin(7));
t607 = sin(qJ(3));
t631 = t602 * t607;
t598 = sin(pkin(7));
t632 = t598 * t607;
t539 = t560 * t631 + t611 * t561 + t573 * t632;
t612 = qJD(3) ^ 2;
t537 = -pkin(3) * t612 + qJDD(3) * pkin(9) + t539;
t549 = -t560 * t598 + t573 * t602;
t606 = sin(qJ(4));
t610 = cos(qJ(4));
t532 = -t606 * t537 + t610 * t549;
t627 = qJD(3) * qJD(4);
t626 = t610 * t627;
t584 = qJDD(3) * t606 + t626;
t530 = (-t584 + t626) * pkin(10) + (t606 * t610 * t612 + qJDD(4)) * pkin(4) + t532;
t533 = t610 * t537 + t606 * t549;
t585 = qJDD(3) * t610 - t606 * t627;
t629 = qJD(3) * t606;
t590 = qJD(4) * pkin(4) - pkin(10) * t629;
t594 = t610 ^ 2;
t531 = -pkin(4) * t594 * t612 + pkin(10) * t585 - qJD(4) * t590 + t533;
t605 = sin(qJ(5));
t609 = cos(qJ(5));
t526 = t605 * t530 + t609 * t531;
t579 = (t605 * t610 + t606 * t609) * qJD(3);
t553 = -qJD(5) * t579 - t584 * t605 + t585 * t609;
t578 = (t605 * t606 - t609 * t610) * qJD(3);
t566 = mrSges(6,1) * t578 + mrSges(6,2) * t579;
t593 = qJD(4) + qJD(5);
t572 = mrSges(6,1) * t593 - mrSges(6,3) * t579;
t592 = qJDD(4) + qJDD(5);
t567 = pkin(5) * t578 - pkin(11) * t579;
t591 = t593 ^ 2;
t524 = -pkin(5) * t591 + pkin(11) * t592 - t567 * t578 + t526;
t538 = -t607 * t561 + (t560 * t602 + t573 * t598) * t611;
t615 = -qJDD(3) * pkin(3) - t538;
t534 = -t585 * pkin(4) + t590 * t629 + (-pkin(10) * t594 - pkin(9)) * t612 + t615;
t554 = -qJD(5) * t578 + t584 * t609 + t585 * t605;
t527 = (t578 * t593 - t554) * pkin(11) + (t579 * t593 - t553) * pkin(5) + t534;
t604 = sin(qJ(6));
t608 = cos(qJ(6));
t521 = -t524 * t604 + t527 * t608;
t569 = -t579 * t604 + t593 * t608;
t542 = qJD(6) * t569 + t554 * t608 + t592 * t604;
t570 = t579 * t608 + t593 * t604;
t550 = -mrSges(7,1) * t569 + mrSges(7,2) * t570;
t552 = qJDD(6) - t553;
t574 = qJD(6) + t578;
t555 = -mrSges(7,2) * t574 + mrSges(7,3) * t569;
t519 = m(7) * t521 + mrSges(7,1) * t552 - mrSges(7,3) * t542 - t550 * t570 + t555 * t574;
t522 = t524 * t608 + t527 * t604;
t541 = -qJD(6) * t570 - t554 * t604 + t592 * t608;
t556 = mrSges(7,1) * t574 - mrSges(7,3) * t570;
t520 = m(7) * t522 - mrSges(7,2) * t552 + mrSges(7,3) * t541 + t550 * t569 - t556 * t574;
t622 = -t519 * t604 + t608 * t520;
t510 = m(6) * t526 - mrSges(6,2) * t592 + mrSges(6,3) * t553 - t566 * t578 - t572 * t593 + t622;
t525 = t530 * t609 - t531 * t605;
t571 = -mrSges(6,2) * t593 - mrSges(6,3) * t578;
t523 = -pkin(5) * t592 - pkin(11) * t591 + t567 * t579 - t525;
t616 = -m(7) * t523 + t541 * mrSges(7,1) - mrSges(7,2) * t542 + t569 * t555 - t556 * t570;
t515 = m(6) * t525 + mrSges(6,1) * t592 - mrSges(6,3) * t554 - t566 * t579 + t571 * t593 + t616;
t504 = t605 * t510 + t609 * t515;
t583 = (-mrSges(5,1) * t610 + mrSges(5,2) * t606) * qJD(3);
t628 = qJD(3) * t610;
t589 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t628;
t502 = m(5) * t532 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t584 + qJD(4) * t589 - t583 * t629 + t504;
t588 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t629;
t623 = t609 * t510 - t515 * t605;
t503 = m(5) * t533 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t585 - qJD(4) * t588 + t583 * t628 + t623;
t624 = -t502 * t606 + t610 * t503;
t493 = m(4) * t539 - mrSges(4,1) * t612 - qJDD(3) * mrSges(4,2) + t624;
t496 = t610 * t502 + t606 * t503;
t495 = m(4) * t549 + t496;
t536 = -t612 * pkin(9) + t615;
t511 = t608 * t519 + t604 * t520;
t614 = m(6) * t534 - t553 * mrSges(6,1) + mrSges(6,2) * t554 + t578 * t571 + t572 * t579 + t511;
t613 = -m(5) * t536 + t585 * mrSges(5,1) - mrSges(5,2) * t584 - t588 * t629 + t589 * t628 - t614;
t507 = m(4) * t538 + qJDD(3) * mrSges(4,1) - mrSges(4,2) * t612 + t613;
t633 = t507 * t611;
t482 = t493 * t631 - t495 * t598 + t602 * t633;
t478 = m(3) * t560 + t482;
t489 = t611 * t493 - t507 * t607;
t488 = m(3) * t561 + t489;
t636 = t478 * t600 + t488 * t596;
t481 = t493 * t632 + t602 * t495 + t598 * t633;
t480 = m(3) * t573 + t481;
t468 = -t480 * t599 + t603 * t636;
t466 = m(2) * t586 + t468;
t474 = -t478 * t596 + t600 * t488;
t473 = m(2) * t587 + t474;
t630 = t601 * t466 + t597 * t473;
t467 = t603 * t480 + t599 * t636;
t625 = -t466 * t597 + t601 * t473;
t543 = Ifges(7,5) * t570 + Ifges(7,6) * t569 + Ifges(7,3) * t574;
t545 = Ifges(7,1) * t570 + Ifges(7,4) * t569 + Ifges(7,5) * t574;
t512 = -mrSges(7,1) * t523 + mrSges(7,3) * t522 + Ifges(7,4) * t542 + Ifges(7,2) * t541 + Ifges(7,6) * t552 - t543 * t570 + t545 * t574;
t544 = Ifges(7,4) * t570 + Ifges(7,2) * t569 + Ifges(7,6) * t574;
t513 = mrSges(7,2) * t523 - mrSges(7,3) * t521 + Ifges(7,1) * t542 + Ifges(7,4) * t541 + Ifges(7,5) * t552 + t543 * t569 - t544 * t574;
t562 = Ifges(6,5) * t579 - Ifges(6,6) * t578 + Ifges(6,3) * t593;
t563 = Ifges(6,4) * t579 - Ifges(6,2) * t578 + Ifges(6,6) * t593;
t497 = mrSges(6,2) * t534 - mrSges(6,3) * t525 + Ifges(6,1) * t554 + Ifges(6,4) * t553 + Ifges(6,5) * t592 - pkin(11) * t511 - t512 * t604 + t513 * t608 - t562 * t578 - t563 * t593;
t564 = Ifges(6,1) * t579 - Ifges(6,4) * t578 + Ifges(6,5) * t593;
t498 = -mrSges(6,1) * t534 - mrSges(7,1) * t521 + mrSges(7,2) * t522 + mrSges(6,3) * t526 + Ifges(6,4) * t554 - Ifges(7,5) * t542 + Ifges(6,2) * t553 + Ifges(6,6) * t592 - Ifges(7,6) * t541 - Ifges(7,3) * t552 - pkin(5) * t511 - t544 * t570 + t545 * t569 - t562 * t579 + t564 * t593;
t575 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t606 + Ifges(5,6) * t610) * qJD(3);
t577 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t606 + Ifges(5,4) * t610) * qJD(3);
t483 = -mrSges(5,1) * t536 + mrSges(5,3) * t533 + Ifges(5,4) * t584 + Ifges(5,2) * t585 + Ifges(5,6) * qJDD(4) - pkin(4) * t614 + pkin(10) * t623 + qJD(4) * t577 + t605 * t497 + t609 * t498 - t575 * t629;
t576 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t606 + Ifges(5,2) * t610) * qJD(3);
t484 = mrSges(5,2) * t536 - mrSges(5,3) * t532 + Ifges(5,1) * t584 + Ifges(5,4) * t585 + Ifges(5,5) * qJDD(4) - pkin(10) * t504 - qJD(4) * t576 + t497 * t609 - t498 * t605 + t575 * t628;
t470 = mrSges(4,2) * t549 - mrSges(4,3) * t538 + Ifges(4,5) * qJDD(3) - Ifges(4,6) * t612 - pkin(9) * t496 - t483 * t606 + t484 * t610;
t475 = t612 * Ifges(4,5) - t579 * t563 - t578 * t564 - pkin(3) * t496 + mrSges(4,3) * t539 - mrSges(4,1) * t549 + Ifges(4,6) * qJDD(3) - Ifges(5,5) * t584 - Ifges(5,6) * t585 - Ifges(5,3) * qJDD(4) - mrSges(5,1) * t532 + mrSges(5,2) * t533 - pkin(4) * t504 - mrSges(6,1) * t525 + mrSges(6,2) * t526 - t604 * t513 - t608 * t512 - pkin(5) * t616 - pkin(11) * t622 - Ifges(6,5) * t554 - Ifges(6,6) * t553 - Ifges(6,3) * t592 + (-t576 * t606 + t577 * t610) * qJD(3);
t618 = pkin(8) * t489 + t470 * t607 + t475 * t611;
t469 = mrSges(4,1) * t538 - mrSges(4,2) * t539 + Ifges(4,3) * qJDD(3) + pkin(3) * t613 + pkin(9) * t624 + t610 * t483 + t606 * t484;
t463 = -mrSges(3,1) * t573 + mrSges(3,3) * t561 - pkin(2) * t481 - t598 * t469 + t602 * t618;
t464 = mrSges(3,2) * t573 - mrSges(3,3) * t560 + t611 * t470 - t607 * t475 + (-t481 * t598 - t482 * t602) * pkin(8);
t617 = qJ(2) * t474 + t463 * t600 + t464 * t596;
t462 = mrSges(3,1) * t560 - mrSges(3,2) * t561 + pkin(2) * t482 + t602 * t469 + t598 * t618;
t461 = mrSges(2,2) * t595 - mrSges(2,3) * t586 - t596 * t463 + t600 * t464 + (-t467 * t599 - t468 * t603) * qJ(2);
t460 = -mrSges(2,1) * t595 + mrSges(2,3) * t587 - pkin(1) * t467 - t599 * t462 + t603 * t617;
t1 = [-m(1) * g(1) + t625; -m(1) * g(2) + t630; -m(1) * g(3) + m(2) * t595 + t467; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t630 - t597 * t460 + t601 * t461; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t625 + t601 * t460 + t597 * t461; -mrSges(1,1) * g(2) + mrSges(2,1) * t586 + mrSges(1,2) * g(1) - mrSges(2,2) * t587 + pkin(1) * t468 + t603 * t462 + t599 * t617;];
tauB  = t1;
