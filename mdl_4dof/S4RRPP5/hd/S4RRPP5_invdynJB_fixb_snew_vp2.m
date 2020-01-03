% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRPP5
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRPP5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:16
% EndTime: 2019-12-31 17:00:18
% DurationCPUTime: 0.82s
% Computational Cost: add. (3603->196), mult. (7350->228), div. (0->0), fcn. (2931->4), ass. (0->80)
t582 = Ifges(3,1) + Ifges(4,2) + Ifges(5,3);
t566 = Ifges(3,4) + Ifges(4,6) - Ifges(5,6);
t565 = Ifges(3,5) - Ifges(4,4) + Ifges(5,5);
t581 = Ifges(3,2) + Ifges(4,3) + Ifges(5,2);
t564 = Ifges(3,6) - Ifges(4,5) - Ifges(5,4);
t580 = Ifges(3,3) + Ifges(4,1) + Ifges(5,1);
t541 = cos(qJ(2));
t568 = qJD(1) * t541;
t524 = -mrSges(4,1) * t568 - qJD(2) * mrSges(4,3);
t540 = sin(qJ(1));
t542 = cos(qJ(1));
t528 = -t542 * g(1) - t540 * g(2);
t544 = qJD(1) ^ 2;
t500 = -t544 * pkin(1) + qJDD(1) * pkin(5) + t528;
t539 = sin(qJ(2));
t483 = -t541 * g(3) - t539 * t500;
t510 = (-pkin(2) * t541 - qJ(3) * t539) * qJD(1);
t543 = qJD(2) ^ 2;
t569 = qJD(1) * t539;
t482 = -qJDD(2) * pkin(2) - t543 * qJ(3) + t510 * t569 + qJDD(3) - t483;
t567 = qJD(1) * qJD(2);
t558 = t541 * t567;
t514 = t539 * qJDD(1) + t558;
t574 = -2 * qJD(4);
t478 = qJD(2) * t574 + (-t539 * t541 * t544 - qJDD(2)) * qJ(4) + (t514 - t558) * pkin(3) + t482;
t525 = mrSges(5,1) * t568 + qJD(2) * mrSges(5,2);
t554 = -m(5) * t478 + qJDD(2) * mrSges(5,3) + qJD(2) * t525;
t551 = m(4) * t482 + t514 * mrSges(4,1) - t554;
t511 = (mrSges(4,2) * t541 - mrSges(4,3) * t539) * qJD(1);
t513 = (-mrSges(5,2) * t539 - mrSges(5,3) * t541) * qJD(1);
t570 = t511 + t513;
t572 = t514 * mrSges(5,1);
t472 = qJDD(2) * mrSges(4,2) + qJD(2) * t524 + t570 * t569 + t551 + t572;
t473 = t513 * t569 - t554 + t572;
t559 = t539 * t567;
t515 = t541 * qJDD(1) - t559;
t522 = pkin(3) * t569 - qJD(2) * qJ(4);
t538 = t541 ^ 2;
t484 = -t539 * g(3) + t541 * t500;
t547 = -t543 * pkin(2) + qJDD(2) * qJ(3) + t510 * t568 + t484;
t479 = -t538 * t544 * qJ(4) + t515 * pkin(3) + qJDD(4) + ((2 * qJD(3)) + t522) * qJD(2) + t547;
t575 = -2 * qJD(3);
t481 = qJD(2) * t575 - t547;
t526 = mrSges(4,1) * t569 + qJD(2) * mrSges(4,2);
t523 = mrSges(5,1) * t569 - qJD(2) * mrSges(5,3);
t555 = m(5) * t479 + qJDD(2) * mrSges(5,2) + qJD(2) * t523 + t513 * t568;
t548 = -m(4) * t481 + qJDD(2) * mrSges(4,3) + qJD(2) * t526 + t511 * t568 + t555;
t560 = t565 * qJD(2) + (t582 * t539 + t566 * t541) * qJD(1);
t561 = t564 * qJD(2) + (t566 * t539 + t581 * t541) * qJD(1);
t579 = (t561 * t539 - t560 * t541) * qJD(1) + t580 * qJDD(2) + t565 * t514 + t564 * t515 + mrSges(3,1) * t483 - mrSges(3,2) * t484 + mrSges(4,2) * t482 + mrSges(5,2) * t479 - mrSges(4,3) * t481 - mrSges(5,3) * t478 - pkin(2) * t472 + qJ(3) * ((mrSges(4,1) + mrSges(5,1)) * t515 + t548) - qJ(4) * t473;
t578 = pkin(2) * t559 + t569 * t575;
t573 = -mrSges(5,1) - mrSges(3,3);
t512 = (-mrSges(3,1) * t541 + mrSges(3,2) * t539) * qJD(1);
t520 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t569;
t468 = t512 * t568 + m(3) * t484 - qJDD(2) * mrSges(3,2) - qJD(2) * t520 + (mrSges(4,1) - t573) * t515 + t548;
t521 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t568;
t469 = m(3) * t483 + t573 * t514 + (mrSges(3,1) - mrSges(4,2)) * qJDD(2) + (t521 - t524) * qJD(2) + (-t512 - t570) * t569 - t551;
t556 = t541 * t468 - t539 * t469;
t459 = m(2) * t528 - t544 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t556;
t527 = t540 * g(1) - t542 * g(2);
t552 = -qJDD(1) * pkin(1) - t527;
t499 = -t544 * pkin(5) + t552;
t480 = -t515 * pkin(2) + (-t514 - t558) * qJ(3) + t499 + t578;
t476 = -t514 * qJ(3) + (-pkin(3) * t538 - pkin(5)) * t544 + (-pkin(2) - qJ(4)) * t515 + (-t522 * t539 + (-qJ(3) * qJD(2) + t574) * t541) * qJD(1) + t552 + t578;
t553 = m(5) * t476 - t514 * mrSges(5,2) - t515 * mrSges(5,3) - t523 * t569 - t525 * t568;
t549 = -m(4) * t480 - t515 * mrSges(4,2) + t526 * t569 - t553;
t545 = -m(3) * t499 + t521 * t568 + t515 * mrSges(3,1) + (-mrSges(3,2) + mrSges(4,3)) * t514 + (-t520 * t539 - t524 * t541) * qJD(1) + t549;
t463 = m(2) * t527 + qJDD(1) * mrSges(2,1) - t544 * mrSges(2,2) + t545;
t571 = t540 * t459 + t542 * t463;
t461 = t539 * t468 + t541 * t469;
t562 = t580 * qJD(2) + (t565 * t539 + t564 * t541) * qJD(1);
t557 = t542 * t459 - t540 * t463;
t470 = -t514 * mrSges(4,3) + t524 * t568 - t549;
t474 = t515 * mrSges(5,1) + t555;
t454 = -mrSges(3,1) * t499 - mrSges(4,1) * t481 + mrSges(5,1) * t479 + mrSges(4,2) * t480 + mrSges(3,3) * t484 - mrSges(5,3) * t476 - pkin(2) * t470 + pkin(3) * t474 - qJ(4) * t553 + t560 * qJD(2) + t564 * qJDD(2) + t566 * t514 + t581 * t515 - t562 * t569;
t456 = mrSges(4,1) * t482 + mrSges(5,1) * t478 + mrSges(3,2) * t499 - mrSges(5,2) * t476 - mrSges(3,3) * t483 - mrSges(4,3) * t480 + pkin(3) * t473 - qJ(3) * t470 - t561 * qJD(2) + t565 * qJDD(2) + t582 * t514 + t566 * t515 + t562 * t568;
t550 = mrSges(2,1) * t527 - mrSges(2,2) * t528 + Ifges(2,3) * qJDD(1) + pkin(1) * t545 + pkin(5) * t556 + t541 * t454 + t539 * t456;
t452 = mrSges(2,1) * g(3) + mrSges(2,3) * t528 + t544 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t461 - t579;
t451 = -mrSges(2,2) * g(3) - mrSges(2,3) * t527 + Ifges(2,5) * qJDD(1) - t544 * Ifges(2,6) - pkin(5) * t461 - t539 * t454 + t541 * t456;
t1 = [-m(1) * g(1) + t557; -m(1) * g(2) + t571; (-m(1) - m(2)) * g(3) + t461; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t571 + t542 * t451 - t540 * t452; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t557 + t540 * t451 + t542 * t452; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t550; t550; t579; t472; t474;];
tauJB = t1;
