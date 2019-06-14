% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-05-05 15:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:31:08
% EndTime: 2019-05-05 15:31:16
% DurationCPUTime: 4.97s
% Computational Cost: add. (63504->294), mult. (118761->357), div. (0->0), fcn. (71380->10), ass. (0->119)
t572 = sin(qJ(1));
t576 = cos(qJ(1));
t549 = t572 * g(1) - g(2) * t576;
t541 = qJDD(1) * pkin(1) + t549;
t550 = -g(1) * t576 - g(2) * t572;
t578 = qJD(1) ^ 2;
t543 = -pkin(1) * t578 + t550;
t567 = sin(pkin(10));
t568 = cos(pkin(10));
t521 = t567 * t541 + t568 * t543;
t601 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t521;
t600 = -pkin(2) - pkin(7);
t599 = mrSges(3,1) - mrSges(4,2);
t598 = -Ifges(4,4) + Ifges(3,5);
t597 = Ifges(4,5) - Ifges(3,6);
t520 = t568 * t541 - t567 * t543;
t584 = -t578 * qJ(3) + qJDD(3) - t520;
t508 = qJDD(1) * t600 + t584;
t564 = -g(3) + qJDD(2);
t571 = sin(qJ(4));
t575 = cos(qJ(4));
t502 = t571 * t508 + t575 * t564;
t542 = (mrSges(5,1) * t571 + mrSges(5,2) * t575) * qJD(1);
t594 = qJD(1) * qJD(4);
t555 = t575 * t594;
t545 = -t571 * qJDD(1) - t555;
t595 = qJD(1) * t575;
t548 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t595;
t557 = t571 * qJD(1);
t505 = t578 * t600 - t601;
t591 = t571 * t594;
t546 = qJDD(1) * t575 - t591;
t492 = (-t546 + t591) * pkin(8) + (-t545 + t555) * pkin(4) + t505;
t544 = (pkin(4) * t571 - pkin(8) * t575) * qJD(1);
t577 = qJD(4) ^ 2;
t499 = -pkin(4) * t577 + qJDD(4) * pkin(8) - t544 * t557 + t502;
t570 = sin(qJ(5));
t574 = cos(qJ(5));
t484 = t574 * t492 - t570 * t499;
t539 = qJD(4) * t574 - t570 * t595;
t517 = qJD(5) * t539 + qJDD(4) * t570 + t546 * t574;
t538 = qJDD(5) - t545;
t540 = qJD(4) * t570 + t574 * t595;
t552 = t557 + qJD(5);
t482 = (t539 * t552 - t517) * pkin(9) + (t539 * t540 + t538) * pkin(5) + t484;
t485 = t570 * t492 + t574 * t499;
t516 = -qJD(5) * t540 + qJDD(4) * t574 - t546 * t570;
t525 = pkin(5) * t552 - pkin(9) * t540;
t537 = t539 ^ 2;
t483 = -pkin(5) * t537 + pkin(9) * t516 - t525 * t552 + t485;
t569 = sin(qJ(6));
t573 = cos(qJ(6));
t480 = t482 * t573 - t483 * t569;
t518 = t539 * t573 - t540 * t569;
t489 = qJD(6) * t518 + t516 * t569 + t517 * t573;
t519 = t539 * t569 + t540 * t573;
t500 = -mrSges(7,1) * t518 + mrSges(7,2) * t519;
t551 = qJD(6) + t552;
t506 = -mrSges(7,2) * t551 + mrSges(7,3) * t518;
t533 = qJDD(6) + t538;
t478 = m(7) * t480 + mrSges(7,1) * t533 - mrSges(7,3) * t489 - t500 * t519 + t506 * t551;
t481 = t482 * t569 + t483 * t573;
t488 = -qJD(6) * t519 + t516 * t573 - t517 * t569;
t507 = mrSges(7,1) * t551 - mrSges(7,3) * t519;
t479 = m(7) * t481 - mrSges(7,2) * t533 + mrSges(7,3) * t488 + t500 * t518 - t507 * t551;
t471 = t573 * t478 + t569 * t479;
t522 = -mrSges(6,1) * t539 + mrSges(6,2) * t540;
t523 = -mrSges(6,2) * t552 + mrSges(6,3) * t539;
t469 = m(6) * t484 + mrSges(6,1) * t538 - mrSges(6,3) * t517 - t522 * t540 + t523 * t552 + t471;
t524 = mrSges(6,1) * t552 - mrSges(6,3) * t540;
t586 = -t478 * t569 + t573 * t479;
t470 = m(6) * t485 - mrSges(6,2) * t538 + mrSges(6,3) * t516 + t522 * t539 - t524 * t552 + t586;
t587 = -t469 * t570 + t574 * t470;
t464 = m(5) * t502 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t545 - qJD(4) * t548 - t542 * t557 + t587;
t501 = t508 * t575 - t571 * t564;
t547 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t557;
t498 = -qJDD(4) * pkin(4) - pkin(8) * t577 + t544 * t595 - t501;
t486 = -pkin(5) * t516 - pkin(9) * t537 + t525 * t540 + t498;
t581 = m(7) * t486 - t488 * mrSges(7,1) + mrSges(7,2) * t489 - t518 * t506 + t507 * t519;
t579 = -m(6) * t498 + t516 * mrSges(6,1) - mrSges(6,2) * t517 + t539 * t523 - t524 * t540 - t581;
t474 = m(5) * t501 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t546 + qJD(4) * t547 - t542 * t595 + t579;
t458 = t571 * t464 + t575 * t474;
t510 = -qJDD(1) * pkin(2) + t584;
t583 = -m(4) * t510 + t578 * mrSges(4,3) - t458;
t455 = m(3) * t520 - t578 * mrSges(3,2) + qJDD(1) * t599 + t583;
t509 = t578 * pkin(2) + t601;
t465 = t574 * t469 + t570 * t470;
t582 = -m(5) * t505 + mrSges(5,1) * t545 - t546 * mrSges(5,2) - t547 * t557 - t548 * t595 - t465;
t580 = -m(4) * t509 + t578 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t582;
t462 = m(3) * t521 - mrSges(3,1) * t578 - qJDD(1) * mrSges(3,2) + t580;
t452 = t568 * t455 + t567 * t462;
t450 = m(2) * t549 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t578 + t452;
t589 = -t455 * t567 + t568 * t462;
t451 = m(2) * t550 - mrSges(2,1) * t578 - qJDD(1) * mrSges(2,2) + t589;
t596 = t576 * t450 + t572 * t451;
t590 = -t450 * t572 + t576 * t451;
t588 = t575 * t464 - t571 * t474;
t457 = m(4) * t564 + t588;
t585 = m(3) * t564 + t457;
t532 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t575 - Ifges(5,4) * t571) * qJD(1);
t531 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t575 - Ifges(5,2) * t571) * qJD(1);
t530 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t575 - Ifges(5,6) * t571) * qJD(1);
t513 = Ifges(6,1) * t540 + Ifges(6,4) * t539 + Ifges(6,5) * t552;
t512 = Ifges(6,4) * t540 + Ifges(6,2) * t539 + Ifges(6,6) * t552;
t511 = Ifges(6,5) * t540 + Ifges(6,6) * t539 + Ifges(6,3) * t552;
t496 = Ifges(7,1) * t519 + Ifges(7,4) * t518 + Ifges(7,5) * t551;
t495 = Ifges(7,4) * t519 + Ifges(7,2) * t518 + Ifges(7,6) * t551;
t494 = Ifges(7,5) * t519 + Ifges(7,6) * t518 + Ifges(7,3) * t551;
t473 = mrSges(7,2) * t486 - mrSges(7,3) * t480 + Ifges(7,1) * t489 + Ifges(7,4) * t488 + Ifges(7,5) * t533 + t494 * t518 - t495 * t551;
t472 = -mrSges(7,1) * t486 + mrSges(7,3) * t481 + Ifges(7,4) * t489 + Ifges(7,2) * t488 + Ifges(7,6) * t533 - t494 * t519 + t496 * t551;
t459 = mrSges(6,2) * t498 - mrSges(6,3) * t484 + Ifges(6,1) * t517 + Ifges(6,4) * t516 + Ifges(6,5) * t538 - pkin(9) * t471 - t472 * t569 + t473 * t573 + t511 * t539 - t512 * t552;
t456 = -mrSges(6,1) * t498 + mrSges(6,3) * t485 + Ifges(6,4) * t517 + Ifges(6,2) * t516 + Ifges(6,6) * t538 - pkin(5) * t581 + pkin(9) * t586 + t573 * t472 + t569 * t473 - t540 * t511 + t552 * t513;
t453 = Ifges(5,4) * t546 + Ifges(5,2) * t545 + Ifges(5,6) * qJDD(4) - t530 * t595 + qJD(4) * t532 - mrSges(5,1) * t505 + mrSges(5,3) * t502 - Ifges(6,5) * t517 - Ifges(6,6) * t516 - Ifges(6,3) * t538 - t540 * t512 + t539 * t513 - mrSges(6,1) * t484 + mrSges(6,2) * t485 - Ifges(7,5) * t489 - Ifges(7,6) * t488 - Ifges(7,3) * t533 - t519 * t495 + t518 * t496 - mrSges(7,1) * t480 + mrSges(7,2) * t481 - pkin(5) * t471 - pkin(4) * t465;
t446 = mrSges(5,2) * t505 - mrSges(5,3) * t501 + Ifges(5,1) * t546 + Ifges(5,4) * t545 + Ifges(5,5) * qJDD(4) - pkin(8) * t465 - qJD(4) * t531 - t456 * t570 + t459 * t574 - t530 * t557;
t445 = -qJ(3) * t457 - mrSges(3,3) * t520 + pkin(3) * t458 + mrSges(4,1) * t510 + t574 * t456 + pkin(4) * t579 + pkin(8) * t587 + t570 * t459 + Ifges(5,5) * t546 + Ifges(5,6) * t545 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t501 - mrSges(5,2) * t502 + t597 * t578 + (mrSges(3,2) - mrSges(4,3)) * t564 + t598 * qJDD(1) + (t531 * t575 + t532 * t571) * qJD(1);
t444 = -mrSges(4,1) * t509 + mrSges(3,3) * t521 - pkin(2) * t457 - pkin(3) * t582 - pkin(7) * t588 - qJDD(1) * t597 - t571 * t446 - t575 * t453 - t564 * t599 + t578 * t598;
t443 = -mrSges(2,2) * g(3) - mrSges(2,3) * t549 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t578 - qJ(2) * t452 - t444 * t567 + t445 * t568;
t442 = mrSges(2,1) * g(3) + mrSges(2,3) * t550 + t578 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t585 + qJ(2) * t589 + t568 * t444 + t567 * t445;
t1 = [-m(1) * g(1) + t590; -m(1) * g(2) + t596; (-m(1) - m(2)) * g(3) + t585; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t596 - t572 * t442 + t576 * t443; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t590 + t576 * t442 + t572 * t443; pkin(1) * t452 + mrSges(2,1) * t549 - mrSges(2,2) * t550 + pkin(2) * t583 + qJ(3) * t580 + mrSges(3,1) * t520 - mrSges(3,2) * t521 + t575 * t446 - t571 * t453 - pkin(7) * t458 + mrSges(4,2) * t510 - mrSges(4,3) * t509 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(2,3) + Ifges(3,3)) * qJDD(1);];
tauB  = t1;
