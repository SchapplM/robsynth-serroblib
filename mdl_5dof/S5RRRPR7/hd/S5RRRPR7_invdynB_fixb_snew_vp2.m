% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR7_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:36
% EndTime: 2019-12-31 21:16:44
% DurationCPUTime: 7.53s
% Computational Cost: add. (92559->313), mult. (192826->399), div. (0->0), fcn. (133906->10), ass. (0->121)
t594 = cos(qJ(3));
t577 = qJD(1) ^ 2;
t593 = pkin(2) * t577;
t573 = sin(qJ(1));
t576 = cos(qJ(1));
t561 = -t576 * g(1) - t573 * g(2);
t549 = -t577 * pkin(1) + qJDD(1) * pkin(6) + t561;
t572 = sin(qJ(2));
t592 = t572 * t549;
t575 = cos(qJ(2));
t588 = qJD(1) * qJD(2);
t555 = t572 * qJDD(1) + t575 * t588;
t515 = qJDD(2) * pkin(2) - t555 * pkin(7) - t592 + (pkin(7) * t588 + t572 * t593 - g(3)) * t575;
t538 = -t572 * g(3) + t575 * t549;
t556 = t575 * qJDD(1) - t572 * t588;
t590 = qJD(1) * t572;
t559 = qJD(2) * pkin(2) - pkin(7) * t590;
t567 = t575 ^ 2;
t516 = t556 * pkin(7) - qJD(2) * t559 - t567 * t593 + t538;
t571 = sin(qJ(3));
t499 = t571 * t515 + t594 * t516;
t547 = (t571 * t575 + t594 * t572) * qJD(1);
t520 = t547 * qJD(3) + t571 * t555 - t594 * t556;
t589 = qJD(1) * t575;
t546 = t571 * t590 - t594 * t589;
t531 = t546 * mrSges(4,1) + t547 * mrSges(4,2);
t566 = qJD(2) + qJD(3);
t540 = t566 * mrSges(4,1) - t547 * mrSges(4,3);
t565 = qJDD(2) + qJDD(3);
t521 = -t546 * qJD(3) + t594 * t555 + t571 * t556;
t560 = t573 * g(1) - t576 * g(2);
t582 = -qJDD(1) * pkin(1) - t560;
t522 = -t556 * pkin(2) + t559 * t590 + (-pkin(7) * t567 - pkin(6)) * t577 + t582;
t489 = (t546 * t566 - t521) * qJ(4) + (t547 * t566 + t520) * pkin(3) + t522;
t530 = t546 * pkin(3) - t547 * qJ(4);
t564 = t566 ^ 2;
t492 = -t564 * pkin(3) + t565 * qJ(4) - t546 * t530 + t499;
t568 = sin(pkin(9));
t569 = cos(pkin(9));
t536 = t569 * t547 + t568 * t566;
t481 = -0.2e1 * qJD(4) * t536 + t569 * t489 - t568 * t492;
t510 = t569 * t521 + t568 * t565;
t535 = -t568 * t547 + t569 * t566;
t479 = (t535 * t546 - t510) * pkin(8) + (t535 * t536 + t520) * pkin(4) + t481;
t482 = 0.2e1 * qJD(4) * t535 + t568 * t489 + t569 * t492;
t509 = -t568 * t521 + t569 * t565;
t525 = t546 * pkin(4) - t536 * pkin(8);
t534 = t535 ^ 2;
t480 = -t534 * pkin(4) + t509 * pkin(8) - t546 * t525 + t482;
t570 = sin(qJ(5));
t574 = cos(qJ(5));
t477 = t574 * t479 - t570 * t480;
t507 = t574 * t535 - t570 * t536;
t488 = t507 * qJD(5) + t570 * t509 + t574 * t510;
t508 = t570 * t535 + t574 * t536;
t497 = -t507 * mrSges(6,1) + t508 * mrSges(6,2);
t542 = qJD(5) + t546;
t500 = -t542 * mrSges(6,2) + t507 * mrSges(6,3);
t519 = qJDD(5) + t520;
t475 = m(6) * t477 + t519 * mrSges(6,1) - t488 * mrSges(6,3) - t508 * t497 + t542 * t500;
t478 = t570 * t479 + t574 * t480;
t487 = -t508 * qJD(5) + t574 * t509 - t570 * t510;
t501 = t542 * mrSges(6,1) - t508 * mrSges(6,3);
t476 = m(6) * t478 - t519 * mrSges(6,2) + t487 * mrSges(6,3) + t507 * t497 - t542 * t501;
t467 = t574 * t475 + t570 * t476;
t512 = -t535 * mrSges(5,1) + t536 * mrSges(5,2);
t523 = -t546 * mrSges(5,2) + t535 * mrSges(5,3);
t465 = m(5) * t481 + t520 * mrSges(5,1) - t510 * mrSges(5,3) - t536 * t512 + t546 * t523 + t467;
t524 = t546 * mrSges(5,1) - t536 * mrSges(5,3);
t583 = -t570 * t475 + t574 * t476;
t466 = m(5) * t482 - t520 * mrSges(5,2) + t509 * mrSges(5,3) + t535 * t512 - t546 * t524 + t583;
t584 = -t568 * t465 + t569 * t466;
t460 = m(4) * t499 - t565 * mrSges(4,2) - t520 * mrSges(4,3) - t546 * t531 - t566 * t540 + t584;
t498 = t594 * t515 - t571 * t516;
t539 = -t566 * mrSges(4,2) - t546 * mrSges(4,3);
t491 = -t565 * pkin(3) - t564 * qJ(4) + t547 * t530 + qJDD(4) - t498;
t483 = -t509 * pkin(4) - t534 * pkin(8) + t536 * t525 + t491;
t581 = m(6) * t483 - t487 * mrSges(6,1) + t488 * mrSges(6,2) - t507 * t500 + t508 * t501;
t579 = -m(5) * t491 + t509 * mrSges(5,1) - t510 * mrSges(5,2) + t535 * t523 - t536 * t524 - t581;
t471 = m(4) * t498 + t565 * mrSges(4,1) - t521 * mrSges(4,3) - t547 * t531 + t566 * t539 + t579;
t453 = t571 * t460 + t594 * t471;
t537 = -t575 * g(3) - t592;
t554 = (-mrSges(3,1) * t575 + mrSges(3,2) * t572) * qJD(1);
t558 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t589;
t451 = m(3) * t537 + qJDD(2) * mrSges(3,1) - t555 * mrSges(3,3) + qJD(2) * t558 - t554 * t590 + t453;
t557 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t590;
t585 = t594 * t460 - t571 * t471;
t452 = m(3) * t538 - qJDD(2) * mrSges(3,2) + t556 * mrSges(3,3) - qJD(2) * t557 + t554 * t589 + t585;
t586 = -t572 * t451 + t575 * t452;
t445 = m(2) * t561 - t577 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t586;
t548 = -t577 * pkin(6) + t582;
t461 = t569 * t465 + t568 * t466;
t580 = m(4) * t522 + t520 * mrSges(4,1) + t521 * mrSges(4,2) + t546 * t539 + t547 * t540 + t461;
t578 = -m(3) * t548 + t556 * mrSges(3,1) - t555 * mrSges(3,2) - t557 * t590 + t558 * t589 - t580;
t457 = m(2) * t560 + qJDD(1) * mrSges(2,1) - t577 * mrSges(2,2) + t578;
t591 = t573 * t445 + t576 * t457;
t446 = t575 * t451 + t572 * t452;
t587 = t576 * t445 - t573 * t457;
t545 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t572 + Ifges(3,4) * t575) * qJD(1);
t544 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t572 + Ifges(3,2) * t575) * qJD(1);
t543 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t572 + Ifges(3,6) * t575) * qJD(1);
t528 = Ifges(4,1) * t547 - Ifges(4,4) * t546 + Ifges(4,5) * t566;
t527 = Ifges(4,4) * t547 - Ifges(4,2) * t546 + Ifges(4,6) * t566;
t526 = Ifges(4,5) * t547 - Ifges(4,6) * t546 + Ifges(4,3) * t566;
t504 = Ifges(5,1) * t536 + Ifges(5,4) * t535 + Ifges(5,5) * t546;
t503 = Ifges(5,4) * t536 + Ifges(5,2) * t535 + Ifges(5,6) * t546;
t502 = Ifges(5,5) * t536 + Ifges(5,6) * t535 + Ifges(5,3) * t546;
t495 = Ifges(6,1) * t508 + Ifges(6,4) * t507 + Ifges(6,5) * t542;
t494 = Ifges(6,4) * t508 + Ifges(6,2) * t507 + Ifges(6,6) * t542;
t493 = Ifges(6,5) * t508 + Ifges(6,6) * t507 + Ifges(6,3) * t542;
t469 = mrSges(6,2) * t483 - mrSges(6,3) * t477 + Ifges(6,1) * t488 + Ifges(6,4) * t487 + Ifges(6,5) * t519 + t507 * t493 - t542 * t494;
t468 = -mrSges(6,1) * t483 + mrSges(6,3) * t478 + Ifges(6,4) * t488 + Ifges(6,2) * t487 + Ifges(6,6) * t519 - t508 * t493 + t542 * t495;
t455 = mrSges(5,2) * t491 - mrSges(5,3) * t481 + Ifges(5,1) * t510 + Ifges(5,4) * t509 + Ifges(5,5) * t520 - pkin(8) * t467 - t570 * t468 + t574 * t469 + t535 * t502 - t546 * t503;
t454 = -mrSges(5,1) * t491 + mrSges(5,3) * t482 + Ifges(5,4) * t510 + Ifges(5,2) * t509 + Ifges(5,6) * t520 - pkin(4) * t581 + pkin(8) * t583 + t574 * t468 + t570 * t469 - t536 * t502 + t546 * t504;
t447 = Ifges(4,4) * t521 + Ifges(4,6) * t565 - t547 * t526 + t566 * t528 - mrSges(4,1) * t522 + mrSges(4,3) * t499 - Ifges(5,5) * t510 - Ifges(5,6) * t509 - t536 * t503 + t535 * t504 - mrSges(5,1) * t481 + mrSges(5,2) * t482 - Ifges(6,5) * t488 - Ifges(6,6) * t487 - Ifges(6,3) * t519 - t508 * t494 + t507 * t495 - mrSges(6,1) * t477 + mrSges(6,2) * t478 - pkin(4) * t467 - pkin(3) * t461 + (-Ifges(4,2) - Ifges(5,3)) * t520;
t442 = mrSges(4,2) * t522 - mrSges(4,3) * t498 + Ifges(4,1) * t521 - Ifges(4,4) * t520 + Ifges(4,5) * t565 - qJ(4) * t461 - t568 * t454 + t569 * t455 - t546 * t526 - t566 * t527;
t441 = mrSges(3,2) * t548 - mrSges(3,3) * t537 + Ifges(3,1) * t555 + Ifges(3,4) * t556 + Ifges(3,5) * qJDD(2) - pkin(7) * t453 - qJD(2) * t544 + t594 * t442 - t571 * t447 + t543 * t589;
t440 = -pkin(1) * t446 + mrSges(2,3) * t561 - pkin(2) * t453 - Ifges(3,5) * t555 - Ifges(3,6) * t556 - Ifges(3,3) * qJDD(2) - mrSges(3,1) * t537 + mrSges(3,2) * t538 - t568 * t455 - t569 * t454 - pkin(3) * t579 - qJ(4) * t584 - Ifges(4,5) * t521 + Ifges(4,6) * t520 - Ifges(4,3) * t565 - mrSges(4,1) * t498 + mrSges(4,2) * t499 + mrSges(2,1) * g(3) + t577 * Ifges(2,5) - t547 * t527 - t546 * t528 + Ifges(2,6) * qJDD(1) + (-t572 * t544 + t575 * t545) * qJD(1);
t439 = -mrSges(3,1) * t548 + mrSges(3,3) * t538 + Ifges(3,4) * t555 + Ifges(3,2) * t556 + Ifges(3,6) * qJDD(2) - pkin(2) * t580 + pkin(7) * t585 + qJD(2) * t545 + t571 * t442 + t594 * t447 - t543 * t590;
t438 = -mrSges(2,2) * g(3) - mrSges(2,3) * t560 + Ifges(2,5) * qJDD(1) - t577 * Ifges(2,6) - pkin(6) * t446 - t572 * t439 + t575 * t441;
t1 = [-m(1) * g(1) + t587; -m(1) * g(2) + t591; (-m(1) - m(2)) * g(3) + t446; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t591 + t576 * t438 - t573 * t440; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t587 + t573 * t438 + t576 * t440; -mrSges(1,1) * g(2) + mrSges(2,1) * t560 + mrSges(1,2) * g(1) - mrSges(2,2) * t561 + Ifges(2,3) * qJDD(1) + pkin(1) * t578 + pkin(6) * t586 + t575 * t439 + t572 * t441;];
tauB = t1;
