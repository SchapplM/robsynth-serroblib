% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPR9_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:34
% EndTime: 2019-12-31 17:09:37
% DurationCPUTime: 2.31s
% Computational Cost: add. (22928->237), mult. (49261->304), div. (0->0), fcn. (30809->8), ass. (0->95)
t530 = sin(qJ(1));
t533 = cos(qJ(1));
t519 = t530 * g(1) - t533 * g(2);
t535 = qJD(1) ^ 2;
t504 = -qJDD(1) * pkin(1) - t535 * pkin(5) - t519;
t529 = sin(qJ(2));
t532 = cos(qJ(2));
t547 = qJD(1) * qJD(2);
t546 = t532 * t547;
t514 = t529 * qJDD(1) + t546;
t522 = t529 * t547;
t515 = t532 * qJDD(1) - t522;
t476 = (-t514 - t546) * qJ(3) + (-t515 + t522) * pkin(2) + t504;
t520 = -t533 * g(1) - t530 * g(2);
t505 = -t535 * pkin(1) + qJDD(1) * pkin(5) + t520;
t490 = -t529 * g(3) + t532 * t505;
t512 = (-pkin(2) * t532 - qJ(3) * t529) * qJD(1);
t534 = qJD(2) ^ 2;
t548 = t532 * qJD(1);
t479 = -t534 * pkin(2) + qJDD(2) * qJ(3) + t512 * t548 + t490;
t526 = sin(pkin(7));
t527 = cos(pkin(7));
t549 = qJD(1) * t529;
t509 = t526 * qJD(2) + t527 * t549;
t463 = -0.2e1 * qJD(3) * t509 + t527 * t476 - t526 * t479;
t494 = t526 * qJDD(2) + t527 * t514;
t508 = t527 * qJD(2) - t526 * t549;
t461 = (-t508 * t548 - t494) * pkin(6) + (t508 * t509 - t515) * pkin(3) + t463;
t464 = 0.2e1 * qJD(3) * t508 + t526 * t476 + t527 * t479;
t493 = t527 * qJDD(2) - t526 * t514;
t495 = -pkin(3) * t548 - t509 * pkin(6);
t507 = t508 ^ 2;
t462 = -t507 * pkin(3) + t493 * pkin(6) + t495 * t548 + t464;
t528 = sin(qJ(4));
t531 = cos(qJ(4));
t460 = t528 * t461 + t531 * t462;
t489 = -t532 * g(3) - t529 * t505;
t478 = -qJDD(2) * pkin(2) - t534 * qJ(3) + t512 * t549 + qJDD(3) - t489;
t465 = -t493 * pkin(3) - t507 * pkin(6) + t509 * t495 + t478;
t487 = t528 * t508 + t531 * t509;
t467 = -t487 * qJD(4) + t531 * t493 - t528 * t494;
t486 = t531 * t508 - t528 * t509;
t468 = t486 * qJD(4) + t528 * t493 + t531 * t494;
t521 = qJD(4) - t548;
t469 = Ifges(5,5) * t487 + Ifges(5,6) * t486 + Ifges(5,3) * t521;
t471 = Ifges(5,1) * t487 + Ifges(5,4) * t486 + Ifges(5,5) * t521;
t511 = qJDD(4) - t515;
t449 = -mrSges(5,1) * t465 + mrSges(5,3) * t460 + Ifges(5,4) * t468 + Ifges(5,2) * t467 + Ifges(5,6) * t511 - t487 * t469 + t521 * t471;
t459 = t531 * t461 - t528 * t462;
t470 = Ifges(5,4) * t487 + Ifges(5,2) * t486 + Ifges(5,6) * t521;
t450 = mrSges(5,2) * t465 - mrSges(5,3) * t459 + Ifges(5,1) * t468 + Ifges(5,4) * t467 + Ifges(5,5) * t511 + t486 * t469 - t521 * t470;
t482 = Ifges(4,5) * t509 + Ifges(4,6) * t508 - Ifges(4,3) * t548;
t484 = Ifges(4,1) * t509 + Ifges(4,4) * t508 - Ifges(4,5) * t548;
t480 = -t521 * mrSges(5,2) + t486 * mrSges(5,3);
t481 = t521 * mrSges(5,1) - t487 * mrSges(5,3);
t539 = m(5) * t465 - t467 * mrSges(5,1) + t468 * mrSges(5,2) - t486 * t480 + t487 * t481;
t473 = -t486 * mrSges(5,1) + t487 * mrSges(5,2);
t456 = m(5) * t459 + t511 * mrSges(5,1) - t468 * mrSges(5,3) - t487 * t473 + t521 * t480;
t457 = m(5) * t460 - t511 * mrSges(5,2) + t467 * mrSges(5,3) + t486 * t473 - t521 * t481;
t543 = -t528 * t456 + t531 * t457;
t430 = -mrSges(4,1) * t478 + mrSges(4,3) * t464 + Ifges(4,4) * t494 + Ifges(4,2) * t493 - Ifges(4,6) * t515 - pkin(3) * t539 + pkin(6) * t543 + t531 * t449 + t528 * t450 - t509 * t482 - t484 * t548;
t448 = t531 * t456 + t528 * t457;
t483 = Ifges(4,4) * t509 + Ifges(4,2) * t508 - Ifges(4,6) * t548;
t436 = mrSges(4,2) * t478 - mrSges(4,3) * t463 + Ifges(4,1) * t494 + Ifges(4,4) * t493 - Ifges(4,5) * t515 - pkin(6) * t448 - t528 * t449 + t531 * t450 + t508 * t482 + t483 * t548;
t488 = -t508 * mrSges(4,1) + t509 * mrSges(4,2);
t541 = mrSges(4,2) * t548 + t508 * mrSges(4,3);
t446 = m(4) * t463 - t515 * mrSges(4,1) - t494 * mrSges(4,3) - t509 * t488 - t541 * t548 + t448;
t492 = -mrSges(4,1) * t548 - t509 * mrSges(4,3);
t447 = m(4) * t464 + t515 * mrSges(4,2) + t493 * mrSges(4,3) + t508 * t488 + t492 * t548 + t543;
t444 = -t526 * t446 + t527 * t447;
t458 = m(4) * t478 - t493 * mrSges(4,1) + t494 * mrSges(4,2) + t509 * t492 - t508 * t541 + t539;
t502 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t529 + Ifges(3,2) * t532) * qJD(1);
t503 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t529 + Ifges(3,4) * t532) * qJD(1);
t551 = mrSges(3,1) * t489 - mrSges(3,2) * t490 + Ifges(3,5) * t514 + Ifges(3,6) * t515 + Ifges(3,3) * qJDD(2) - pkin(2) * t458 + qJ(3) * t444 + t527 * t430 + t526 * t436 + (t529 * t502 - t532 * t503) * qJD(1);
t513 = (-mrSges(3,1) * t532 + mrSges(3,2) * t529) * qJD(1);
t517 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t549;
t442 = m(3) * t490 - qJDD(2) * mrSges(3,2) + t515 * mrSges(3,3) - qJD(2) * t517 + t513 * t548 + t444;
t518 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t548;
t452 = m(3) * t489 + qJDD(2) * mrSges(3,1) - t514 * mrSges(3,3) + qJD(2) * t518 - t513 * t549 - t458;
t544 = t532 * t442 - t529 * t452;
t433 = m(2) * t520 - t535 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t544;
t443 = t527 * t446 + t526 * t447;
t538 = -m(3) * t504 + t515 * mrSges(3,1) - t514 * mrSges(3,2) - t517 * t549 + t518 * t548 - t443;
t438 = m(2) * t519 + qJDD(1) * mrSges(2,1) - t535 * mrSges(2,2) + t538;
t550 = t530 * t433 + t533 * t438;
t435 = t529 * t442 + t532 * t452;
t545 = t533 * t433 - t530 * t438;
t501 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t529 + Ifges(3,6) * t532) * qJD(1);
t427 = mrSges(3,2) * t504 - mrSges(3,3) * t489 + Ifges(3,1) * t514 + Ifges(3,4) * t515 + Ifges(3,5) * qJDD(2) - qJ(3) * t443 - qJD(2) * t502 - t526 * t430 + t527 * t436 + t501 * t548;
t537 = mrSges(5,1) * t459 - mrSges(5,2) * t460 + Ifges(5,5) * t468 + Ifges(5,6) * t467 + Ifges(5,3) * t511 + t487 * t470 - t486 * t471;
t429 = -t501 * t549 + (Ifges(3,2) + Ifges(4,3)) * t515 + Ifges(3,6) * qJDD(2) + Ifges(3,4) * t514 + qJD(2) * t503 - mrSges(3,1) * t504 + t508 * t484 - t509 * t483 + mrSges(3,3) * t490 - Ifges(4,6) * t493 - Ifges(4,5) * t494 - mrSges(4,1) * t463 + mrSges(4,2) * t464 - pkin(3) * t448 - t537 - pkin(2) * t443;
t540 = mrSges(2,1) * t519 - mrSges(2,2) * t520 + Ifges(2,3) * qJDD(1) + pkin(1) * t538 + pkin(5) * t544 + t529 * t427 + t532 * t429;
t425 = mrSges(2,1) * g(3) + mrSges(2,3) * t520 + t535 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t435 - t551;
t424 = -mrSges(2,2) * g(3) - mrSges(2,3) * t519 + Ifges(2,5) * qJDD(1) - t535 * Ifges(2,6) - pkin(5) * t435 + t532 * t427 - t529 * t429;
t1 = [-m(1) * g(1) + t545; -m(1) * g(2) + t550; (-m(1) - m(2)) * g(3) + t435; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t550 + t533 * t424 - t530 * t425; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t545 + t530 * t424 + t533 * t425; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t540; t540; t551; t458; t537;];
tauJB = t1;
