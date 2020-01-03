% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRRP3
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-12-31 16:27
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRRP3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:26:48
% EndTime: 2019-12-31 16:26:49
% DurationCPUTime: 0.73s
% Computational Cost: add. (4310->164), mult. (8369->201), div. (0->0), fcn. (4132->6), ass. (0->76)
t486 = Ifges(4,1) + Ifges(5,1);
t478 = Ifges(4,4) + Ifges(5,4);
t477 = Ifges(4,5) + Ifges(5,5);
t485 = Ifges(4,2) + Ifges(5,2);
t476 = Ifges(4,6) + Ifges(5,6);
t484 = Ifges(4,3) + Ifges(5,3);
t448 = sin(qJ(3));
t450 = cos(qJ(3));
t426 = (-mrSges(5,1) * t450 + mrSges(5,2) * t448) * qJD(2);
t467 = qJD(2) * qJD(3);
t462 = t450 * t467;
t428 = t448 * qJDD(2) + t462;
t446 = sin(pkin(6));
t447 = cos(pkin(6));
t431 = t446 * g(1) - t447 * g(2);
t432 = -t447 * g(1) - t446 * g(2);
t449 = sin(qJ(2));
t451 = cos(qJ(2));
t410 = t449 * t431 + t451 * t432;
t452 = qJD(2) ^ 2;
t407 = -t452 * pkin(2) + qJDD(2) * pkin(5) + t410;
t445 = -g(3) + qJDD(1);
t439 = t450 * t445;
t466 = qJD(2) * qJD(4);
t480 = pkin(3) * t452;
t400 = qJDD(3) * pkin(3) + t439 + (-t428 + t462) * qJ(4) + (t450 * t480 - t407 - 0.2e1 * t466) * t448;
t468 = qJD(2) * t450;
t436 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t468;
t465 = m(5) * t400 + qJDD(3) * mrSges(5,1) + qJD(3) * t436;
t469 = qJD(2) * t448;
t396 = -t428 * mrSges(5,3) - t426 * t469 + t465;
t404 = t450 * t407 + t448 * t445;
t429 = t450 * qJDD(2) - t448 * t467;
t433 = qJD(3) * pkin(3) - qJ(4) * t469;
t444 = t450 ^ 2;
t401 = t429 * qJ(4) - qJD(3) * t433 - t444 * t480 + 0.2e1 * t450 * t466 + t404;
t403 = -t448 * t407 + t439;
t471 = t477 * qJD(3) + (t486 * t448 + t478 * t450) * qJD(2);
t472 = t476 * qJD(3) + (t478 * t448 + t485 * t450) * qJD(2);
t483 = mrSges(4,1) * t403 + mrSges(5,1) * t400 - mrSges(4,2) * t404 - mrSges(5,2) * t401 + pkin(3) * t396 + (t472 * t448 - t471 * t450) * qJD(2) + t484 * qJDD(3) + t477 * t428 + t476 * t429;
t479 = -mrSges(4,2) - mrSges(5,2);
t427 = (-mrSges(4,1) * t450 + mrSges(4,2) * t448) * qJD(2);
t437 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t468;
t393 = m(4) * t403 + qJDD(3) * mrSges(4,1) + qJD(3) * t437 + (-mrSges(4,3) - mrSges(5,3)) * t428 + (-t426 - t427) * t469 + t465;
t464 = m(5) * t401 + t429 * mrSges(5,3) + t426 * t468;
t434 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t469;
t470 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t469 - t434;
t394 = m(4) * t404 + t429 * mrSges(4,3) + t470 * qJD(3) + t479 * qJDD(3) + t427 * t468 + t464;
t459 = -t448 * t393 + t450 * t394;
t385 = m(3) * t410 - t452 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t459;
t409 = t451 * t431 - t449 * t432;
t456 = -qJDD(2) * pkin(2) - t409;
t406 = -t452 * pkin(5) + t456;
t402 = t433 * t469 - t429 * pkin(3) + qJDD(4) + (-qJ(4) * t444 - pkin(5)) * t452 + t456;
t457 = -m(5) * t402 + t429 * mrSges(5,1) + t436 * t468;
t453 = -m(4) * t406 + t429 * mrSges(4,1) + t479 * t428 + t437 * t468 + t470 * t469 + t457;
t389 = m(3) * t409 + qJDD(2) * mrSges(3,1) - t452 * mrSges(3,2) + t453;
t378 = t449 * t385 + t451 * t389;
t376 = m(2) * t431 + t378;
t460 = t451 * t385 - t449 * t389;
t377 = m(2) * t432 + t460;
t474 = t447 * t376 + t446 * t377;
t387 = t450 * t393 + t448 * t394;
t473 = t484 * qJD(3) + (t477 * t448 + t476 * t450) * qJD(2);
t463 = m(3) * t445 + t387;
t461 = -t446 * t376 + t447 * t377;
t458 = m(2) * t445 + t463;
t397 = t428 * mrSges(5,2) + t434 * t469 - t457;
t380 = -mrSges(4,1) * t406 + mrSges(4,3) * t404 - mrSges(5,1) * t402 + mrSges(5,3) * t401 - pkin(3) * t397 + qJ(4) * t464 + t485 * t429 + t478 * t428 + (-qJ(4) * mrSges(5,2) + t476) * qJDD(3) + (-qJ(4) * t434 + t471) * qJD(3) - t473 * t469;
t382 = mrSges(4,2) * t406 + mrSges(5,2) * t402 - mrSges(4,3) * t403 - mrSges(5,3) * t400 - qJ(4) * t396 - t472 * qJD(3) + t477 * qJDD(3) + t486 * t428 + t478 * t429 + t473 * t468;
t455 = mrSges(3,1) * t409 - mrSges(3,2) * t410 + Ifges(3,3) * qJDD(2) + pkin(2) * t453 + pkin(5) * t459 + t450 * t380 + t448 * t382;
t372 = -mrSges(3,1) * t445 + mrSges(3,3) * t410 + t452 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t387 - t483;
t371 = mrSges(3,2) * t445 - mrSges(3,3) * t409 + Ifges(3,5) * qJDD(2) - t452 * Ifges(3,6) - pkin(5) * t387 - t448 * t380 + t450 * t382;
t370 = mrSges(2,2) * t445 - mrSges(2,3) * t431 - pkin(4) * t378 + t451 * t371 - t449 * t372;
t369 = -mrSges(2,1) * t445 + mrSges(2,3) * t432 - pkin(1) * t463 + pkin(4) * t460 + t449 * t371 + t451 * t372;
t1 = [-m(1) * g(1) + t461; -m(1) * g(2) + t474; -m(1) * g(3) + t458; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t474 - t446 * t369 + t447 * t370; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t461 + t447 * t369 + t446 * t370; -mrSges(1,1) * g(2) + mrSges(2,1) * t431 + mrSges(1,2) * g(1) - mrSges(2,2) * t432 + pkin(1) * t378 + t455; t458; t455; t483; t397;];
tauJB = t1;
