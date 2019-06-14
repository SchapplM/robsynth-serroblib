% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 14:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S6RPPRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_invdynB_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:57:33
% EndTime: 2019-05-05 14:57:35
% DurationCPUTime: 1.83s
% Computational Cost: add. (14008->267), mult. (25896->300), div. (0->0), fcn. (13191->6), ass. (0->99)
t602 = Ifges(6,1) + Ifges(7,1);
t594 = Ifges(6,4) + Ifges(7,4);
t593 = Ifges(6,5) + Ifges(7,5);
t601 = Ifges(6,2) + Ifges(7,2);
t600 = Ifges(6,6) + Ifges(7,6);
t599 = Ifges(6,3) + Ifges(7,3);
t559 = sin(qJ(1));
t562 = cos(qJ(1));
t540 = t559 * g(1) - t562 * g(2);
t564 = qJD(1) ^ 2;
t520 = -qJDD(1) * pkin(1) - t564 * qJ(2) + qJDD(2) - t540;
t512 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t520;
t541 = -t562 * g(1) - t559 * g(2);
t598 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t541;
t597 = -m(3) - m(4);
t596 = mrSges(2,1) - mrSges(3,2);
t595 = -mrSges(6,2) - mrSges(7,2);
t519 = t564 * pkin(1) - t598;
t513 = qJDD(3) + (-pkin(1) - qJ(3)) * t564 + t598;
t510 = -qJDD(1) * pkin(7) + t513;
t558 = sin(qJ(4));
t561 = cos(qJ(4));
t503 = -t561 * g(3) + t558 * t510;
t534 = (mrSges(5,1) * t558 + mrSges(5,2) * t561) * qJD(1);
t584 = qJD(1) * qJD(4);
t576 = t561 * t584;
t536 = -t558 * qJDD(1) - t576;
t586 = qJD(1) * t561;
t539 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t586;
t509 = -t564 * pkin(7) - t512;
t577 = t558 * t584;
t537 = t561 * qJDD(1) - t577;
t484 = (-t537 + t577) * pkin(8) + (-t536 + t576) * pkin(4) + t509;
t535 = (pkin(4) * t558 - pkin(8) * t561) * qJD(1);
t563 = qJD(4) ^ 2;
t585 = t558 * qJD(1);
t487 = -t563 * pkin(4) + qJDD(4) * pkin(8) - t535 * t585 + t503;
t557 = sin(qJ(5));
t560 = cos(qJ(5));
t480 = t560 * t484 - t557 * t487;
t532 = t560 * qJD(4) - t557 * t586;
t501 = t532 * qJD(5) + t557 * qJDD(4) + t560 * t537;
t533 = t557 * qJD(4) + t560 * t586;
t506 = -t532 * mrSges(7,1) + t533 * mrSges(7,2);
t507 = -t532 * mrSges(6,1) + t533 * mrSges(6,2);
t542 = qJD(5) + t585;
t515 = -t542 * mrSges(6,2) + t532 * mrSges(6,3);
t531 = qJDD(5) - t536;
t476 = -0.2e1 * qJD(6) * t533 + (t532 * t542 - t501) * qJ(6) + (t532 * t533 + t531) * pkin(5) + t480;
t514 = -t542 * mrSges(7,2) + t532 * mrSges(7,3);
t579 = m(7) * t476 + t531 * mrSges(7,1) + t542 * t514;
t468 = m(6) * t480 + t531 * mrSges(6,1) + t542 * t515 + (-t506 - t507) * t533 + (-mrSges(6,3) - mrSges(7,3)) * t501 + t579;
t481 = t557 * t484 + t560 * t487;
t500 = -t533 * qJD(5) + t560 * qJDD(4) - t557 * t537;
t516 = t542 * pkin(5) - t533 * qJ(6);
t530 = t532 ^ 2;
t478 = -t530 * pkin(5) + t500 * qJ(6) + 0.2e1 * qJD(6) * t532 - t542 * t516 + t481;
t578 = m(7) * t478 + t500 * mrSges(7,3) + t532 * t506;
t517 = t542 * mrSges(7,1) - t533 * mrSges(7,3);
t587 = -t542 * mrSges(6,1) + t533 * mrSges(6,3) - t517;
t471 = m(6) * t481 + t500 * mrSges(6,3) + t532 * t507 + t595 * t531 + t587 * t542 + t578;
t573 = -t557 * t468 + t560 * t471;
t464 = m(5) * t503 - qJDD(4) * mrSges(5,2) + t536 * mrSges(5,3) - qJD(4) * t539 - t534 * t585 + t573;
t502 = t558 * g(3) + t561 * t510;
t538 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t585;
t486 = -qJDD(4) * pkin(4) - t563 * pkin(8) + t535 * t586 - t502;
t479 = -t500 * pkin(5) - t530 * qJ(6) + t533 * t516 + qJDD(6) + t486;
t571 = m(7) * t479 - t500 * mrSges(7,1) - t532 * t514;
t565 = -m(6) * t486 + t500 * mrSges(6,1) + t595 * t501 + t532 * t515 + t587 * t533 - t571;
t473 = m(5) * t502 + qJDD(4) * mrSges(5,1) - t537 * mrSges(5,3) + qJD(4) * t538 - t534 * t586 + t565;
t457 = t558 * t464 + t561 * t473;
t572 = -m(4) * t513 - qJDD(1) * mrSges(4,2) - t457;
t567 = -m(3) * t519 + t564 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t572;
t455 = m(2) * t541 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t564 + t567;
t466 = t560 * t468 + t557 * t471;
t568 = -m(5) * t509 + t536 * mrSges(5,1) - t537 * mrSges(5,2) - t538 * t585 - t539 * t586 - t466;
t461 = m(4) * t512 - t564 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t568;
t566 = -m(3) * t520 + t564 * mrSges(3,3) - t461;
t460 = m(2) * t540 - t564 * mrSges(2,2) + t596 * qJDD(1) + t566;
t591 = t559 * t455 + t562 * t460;
t590 = t600 * t532 + t593 * t533 + t599 * t542;
t589 = -t601 * t532 - t594 * t533 - t600 * t542;
t588 = t594 * t532 + t602 * t533 + t593 * t542;
t581 = Ifges(2,5) - Ifges(3,4) + Ifges(4,5);
t580 = Ifges(2,6) - Ifges(3,5) - Ifges(4,4);
t575 = t562 * t455 - t559 * t460;
t574 = t561 * t464 - t558 * t473;
t524 = (Ifges(5,5) * qJD(4)) + (Ifges(5,1) * t561 - Ifges(5,4) * t558) * qJD(1);
t523 = (Ifges(5,6) * qJD(4)) + (Ifges(5,4) * t561 - Ifges(5,2) * t558) * qJD(1);
t522 = (Ifges(5,3) * qJD(4)) + (Ifges(5,5) * t561 - Ifges(5,6) * t558) * qJD(1);
t474 = -t501 * mrSges(7,3) - t533 * t506 + t579;
t465 = mrSges(6,2) * t486 + mrSges(7,2) * t479 - mrSges(6,3) * t480 - mrSges(7,3) * t476 - qJ(6) * t474 + t594 * t500 + t602 * t501 + t593 * t531 + t590 * t532 + t589 * t542;
t458 = -mrSges(6,1) * t486 + mrSges(6,3) * t481 - mrSges(7,1) * t479 + mrSges(7,3) * t478 - pkin(5) * t571 + qJ(6) * t578 + (-qJ(6) * t517 + t588) * t542 + (-pkin(5) * t517 - t590) * t533 + (-qJ(6) * mrSges(7,2) + t600) * t531 + (-pkin(5) * mrSges(7,2) + t594) * t501 + t601 * t500;
t456 = t597 * g(3) + t574;
t452 = -t522 * t586 - mrSges(5,1) * t509 - mrSges(6,1) * t480 - mrSges(7,1) * t476 + mrSges(6,2) * t481 + mrSges(7,2) * t478 + mrSges(5,3) * t503 + Ifges(5,4) * t537 + Ifges(5,2) * t536 + Ifges(5,6) * qJDD(4) - pkin(4) * t466 - pkin(5) * t474 + qJD(4) * t524 + t589 * t533 + t588 * t532 - t599 * t531 - t593 * t501 - t600 * t500;
t451 = mrSges(5,2) * t509 - mrSges(5,3) * t502 + Ifges(5,1) * t537 + Ifges(5,4) * t536 + Ifges(5,5) * qJDD(4) - pkin(8) * t466 - qJD(4) * t523 - t557 * t458 + t560 * t465 - t522 * t585;
t450 = -pkin(1) * t456 + pkin(3) * t457 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t502 - mrSges(5,2) * t503 + mrSges(4,1) * t513 - mrSges(3,1) * t519 + pkin(4) * t565 + Ifges(5,6) * t536 + Ifges(5,5) * t537 + mrSges(2,3) * t541 + t557 * t465 + pkin(8) * t573 - qJ(3) * t574 + t560 * t458 - pkin(2) * t572 + (t561 * t523 + t558 * t524) * qJD(1) + (-pkin(2) * mrSges(4,3) + t581) * t564 + t580 * qJDD(1) + (qJ(3) * m(4) + mrSges(4,3) + t596) * g(3);
t449 = -qJ(2) * t456 - mrSges(2,3) * t540 + pkin(2) * t461 + mrSges(3,1) * t520 + t558 * t451 + t561 * t452 + pkin(3) * t568 + pkin(7) * t574 + mrSges(4,1) * t512 - t580 * t564 + t581 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3) + mrSges(4,2)) * g(3);
t1 = [-m(1) * g(1) + t575; -m(1) * g(2) + t591; (-m(1) - m(2) + t597) * g(3) + t574; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t591 + t562 * t449 - t559 * t450; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t575 + t559 * t449 + t562 * t450; pkin(1) * t566 + mrSges(2,1) * t540 - mrSges(2,2) * t541 + qJ(2) * (-t564 * mrSges(4,3) + t567) - qJ(3) * t461 - mrSges(3,3) * t519 + mrSges(3,2) * t520 - t558 * t452 - pkin(7) * t457 + t561 * t451 + mrSges(4,2) * t513 - mrSges(4,3) * t512 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1);];
tauB  = t1;
