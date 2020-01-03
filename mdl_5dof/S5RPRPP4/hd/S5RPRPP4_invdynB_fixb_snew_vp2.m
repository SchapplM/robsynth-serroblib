% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRPP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:24
% EndTime: 2019-12-31 18:14:26
% DurationCPUTime: 1.50s
% Computational Cost: add. (10799->241), mult. (22831->285), div. (0->0), fcn. (12616->6), ass. (0->94)
t518 = Ifges(5,1) + Ifges(6,1);
t511 = Ifges(5,4) - Ifges(6,5);
t509 = Ifges(5,5) + Ifges(6,4);
t517 = Ifges(5,2) + Ifges(6,3);
t516 = (-Ifges(6,2) - Ifges(5,3));
t507 = (Ifges(5,6) - Ifges(6,6));
t479 = sin(qJ(1));
t481 = cos(qJ(1));
t464 = -t481 * g(1) - t479 * g(2);
t489 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t464;
t515 = -2 * qJD(4);
t514 = -pkin(1) - pkin(6);
t513 = mrSges(2,1) - mrSges(3,2);
t512 = -mrSges(5,3) - mrSges(6,2);
t510 = Ifges(2,5) - Ifges(3,4);
t508 = (-Ifges(2,6) + Ifges(3,5));
t463 = t479 * g(1) - t481 * g(2);
t483 = qJD(1) ^ 2;
t488 = -t483 * qJ(2) + qJDD(2) - t463;
t443 = t514 * qJDD(1) + t488;
t478 = sin(qJ(3));
t480 = cos(qJ(3));
t431 = t478 * g(3) + t480 * t443;
t499 = qJD(1) * qJD(3);
t496 = t478 * t499;
t459 = t480 * qJDD(1) - t496;
t412 = (-t459 - t496) * qJ(4) + (-t478 * t480 * t483 + qJDD(3)) * pkin(3) + t431;
t432 = -t480 * g(3) + t478 * t443;
t458 = -t478 * qJDD(1) - t480 * t499;
t500 = qJD(1) * t480;
t461 = qJD(3) * pkin(3) - qJ(4) * t500;
t473 = t478 ^ 2;
t413 = -t473 * t483 * pkin(3) + t458 * qJ(4) - qJD(3) * t461 + t432;
t476 = sin(pkin(7));
t477 = cos(pkin(7));
t447 = (t476 * t480 + t477 * t478) * qJD(1);
t409 = t476 * t412 + t477 * t413 + t447 * t515;
t429 = -t477 * t458 + t476 * t459;
t501 = qJD(1) * t478;
t448 = -t476 * t501 + t477 * t500;
t440 = (qJD(3) * mrSges(5,1)) - t448 * mrSges(5,3);
t424 = t447 * pkin(4) - t448 * qJ(5);
t482 = qJD(3) ^ 2;
t404 = -t482 * pkin(4) + qJDD(3) * qJ(5) + (2 * qJD(5) * qJD(3)) - t447 * t424 + t409;
t441 = -(qJD(3) * mrSges(6,1)) + t448 * mrSges(6,2);
t497 = m(6) * t404 + qJDD(3) * mrSges(6,3) + qJD(3) * t441;
t425 = t447 * mrSges(6,1) - t448 * mrSges(6,3);
t502 = -t447 * mrSges(5,1) - t448 * mrSges(5,2) - t425;
t400 = m(5) * t409 - qJDD(3) * mrSges(5,2) - qJD(3) * t440 + t512 * t429 + t502 * t447 + t497;
t491 = -t477 * t412 + t476 * t413;
t408 = t448 * t515 - t491;
t430 = t476 * t458 + t477 * t459;
t439 = -qJD(3) * mrSges(5,2) - t447 * mrSges(5,3);
t405 = -qJDD(3) * pkin(4) - t482 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t424) * t448 + t491;
t442 = -t447 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t492 = -m(6) * t405 + qJDD(3) * mrSges(6,1) + qJD(3) * t442;
t401 = m(5) * t408 + qJDD(3) * mrSges(5,1) + qJD(3) * t439 + t512 * t430 + t502 * t448 + t492;
t393 = t476 * t400 + t477 * t401;
t457 = (mrSges(4,1) * t478 + mrSges(4,2) * t480) * qJD(1);
t460 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t501;
t391 = m(4) * t431 + qJDD(3) * mrSges(4,1) - t459 * mrSges(4,3) + qJD(3) * t460 - t457 * t500 + t393;
t462 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t500;
t493 = t477 * t400 - t476 * t401;
t392 = m(4) * t432 - qJDD(3) * mrSges(4,2) + t458 * mrSges(4,3) - qJD(3) * t462 - t457 * t501 + t493;
t387 = t480 * t391 + t478 * t392;
t446 = -qJDD(1) * pkin(1) + t488;
t487 = -m(3) * t446 + (t483 * mrSges(3,3)) - t387;
t385 = m(2) * t463 - (t483 * mrSges(2,2)) + t513 * qJDD(1) + t487;
t444 = t483 * pkin(1) - t489;
t438 = t514 * t483 + t489;
t415 = -t458 * pkin(3) + qJDD(4) + t461 * t500 + (-qJ(4) * t473 + t514) * t483 + t489;
t407 = -0.2e1 * qJD(5) * t448 + (qJD(3) * t447 - t430) * qJ(5) + (qJD(3) * t448 + t429) * pkin(4) + t415;
t402 = m(6) * t407 + t429 * mrSges(6,1) - t430 * mrSges(6,3) - t448 * t441 + t447 * t442;
t486 = m(5) * t415 + t429 * mrSges(5,1) + t430 * mrSges(5,2) + t447 * t439 + t448 * t440 + t402;
t485 = -m(4) * t438 + t458 * mrSges(4,1) - t459 * mrSges(4,2) - t460 * t501 - t462 * t500 - t486;
t484 = -m(3) * t444 + (t483 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t485;
t396 = m(2) * t464 - (t483 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t484;
t506 = t481 * t385 + t479 * t396;
t505 = -(t507 * qJD(3)) + t517 * t447 - t511 * t448;
t504 = (t516 * qJD(3)) + t507 * t447 - t509 * t448;
t503 = t509 * qJD(3) - t511 * t447 + t518 * t448;
t495 = -t479 * t385 + t481 * t396;
t494 = -t478 * t391 + t480 * t392;
t451 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t480 - Ifges(4,4) * t478) * qJD(1);
t450 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t480 - Ifges(4,2) * t478) * qJD(1);
t449 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t480 - Ifges(4,6) * t478) * qJD(1);
t389 = mrSges(5,2) * t415 + mrSges(6,2) * t405 - mrSges(5,3) * t408 - mrSges(6,3) * t407 - qJ(5) * t402 + t505 * qJD(3) + t509 * qJDD(3) - t511 * t429 + t518 * t430 + t504 * t447;
t388 = -mrSges(5,1) * t415 - mrSges(6,1) * t407 + mrSges(6,2) * t404 + mrSges(5,3) * t409 - pkin(4) * t402 + t503 * qJD(3) + t507 * qJDD(3) - t517 * t429 + t511 * t430 + t504 * t448;
t386 = -m(3) * g(3) + t494;
t383 = mrSges(4,2) * t438 - mrSges(4,3) * t431 + Ifges(4,1) * t459 + Ifges(4,4) * t458 + Ifges(4,5) * qJDD(3) - qJ(4) * t393 - qJD(3) * t450 - t476 * t388 + t477 * t389 - t449 * t501;
t382 = -mrSges(4,1) * t438 + mrSges(4,3) * t432 + Ifges(4,4) * t459 + Ifges(4,2) * t458 + Ifges(4,6) * qJDD(3) - pkin(3) * t486 + qJ(4) * t493 + qJD(3) * t451 + t477 * t388 + t476 * t389 - t449 * t500;
t381 = qJ(5) * t497 + pkin(4) * t492 + Ifges(4,6) * t458 + Ifges(4,5) * t459 - mrSges(2,3) * t463 + mrSges(3,1) * t446 + mrSges(4,1) * t431 - mrSges(4,2) * t432 + mrSges(5,1) * t408 - mrSges(5,2) * t409 + mrSges(6,3) * t404 - mrSges(6,1) * t405 + pkin(3) * t393 + pkin(2) * t387 - qJ(2) * t386 + (t508 * t483) + (-pkin(4) * t425 - t505) * t448 + (-qJ(5) * t425 + t503) * t447 + (-pkin(4) * mrSges(6,2) + t509) * t430 + (-qJ(5) * mrSges(6,2) - t507) * t429 + t510 * qJDD(1) + (t480 * t450 + t478 * t451) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (Ifges(4,3) - t516) * qJDD(3);
t380 = -mrSges(3,1) * t444 + mrSges(2,3) * t464 - pkin(1) * t386 - pkin(2) * t485 - pkin(6) * t494 + t513 * g(3) - t508 * qJDD(1) - t480 * t382 - t478 * t383 + t510 * t483;
t1 = [-m(1) * g(1) + t495; -m(1) * g(2) + t506; (-m(1) - m(2) - m(3)) * g(3) + t494; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t506 - t479 * t380 + t481 * t381; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t495 + t481 * t380 + t479 * t381; pkin(1) * t487 + qJ(2) * t484 + t480 * t383 - t478 * t382 - pkin(6) * t387 + mrSges(2,1) * t463 - mrSges(2,2) * t464 + mrSges(3,2) * t446 - mrSges(3,3) * t444 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
