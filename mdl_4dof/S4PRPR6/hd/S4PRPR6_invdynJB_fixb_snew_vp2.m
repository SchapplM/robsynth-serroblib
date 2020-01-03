% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PRPR6
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PRPR6_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:25
% EndTime: 2019-12-31 16:24:26
% DurationCPUTime: 0.96s
% Computational Cost: add. (7937->164), mult. (16964->213), div. (0->0), fcn. (10687->8), ass. (0->80)
t439 = qJD(2) ^ 2;
t432 = sin(pkin(7));
t464 = t432 ^ 2;
t434 = cos(pkin(7));
t463 = pkin(3) * t434;
t462 = mrSges(4,2) * t432;
t461 = cos(pkin(6));
t430 = t434 ^ 2;
t460 = t430 * t439;
t433 = sin(pkin(6));
t421 = -t461 * g(1) - t433 * g(2);
t431 = -g(3) + qJDD(1);
t436 = sin(qJ(2));
t438 = cos(qJ(2));
t411 = t438 * t421 + t436 * t431;
t407 = -t439 * pkin(2) + qJDD(2) * qJ(3) + t411;
t420 = t433 * g(1) - t461 * g(2);
t456 = qJD(2) * qJD(3);
t458 = -t434 * t420 - 0.2e1 * t432 * t456;
t392 = (-pkin(5) * qJDD(2) + t439 * t463 - t407) * t432 + t458;
t395 = -t432 * t420 + (t407 + 0.2e1 * t456) * t434;
t455 = qJDD(2) * t434;
t393 = -pkin(3) * t460 + pkin(5) * t455 + t395;
t435 = sin(qJ(4));
t437 = cos(qJ(4));
t390 = t437 * t392 - t435 * t393;
t445 = -t432 * t435 + t434 * t437;
t412 = t445 * qJD(2);
t446 = t432 * t437 + t434 * t435;
t413 = t446 * qJD(2);
t401 = -t412 * mrSges(5,1) + t413 * mrSges(5,2);
t404 = t412 * qJD(4) + t446 * qJDD(2);
t408 = -qJD(4) * mrSges(5,2) + t412 * mrSges(5,3);
t387 = m(5) * t390 + qJDD(4) * mrSges(5,1) - t404 * mrSges(5,3) + qJD(4) * t408 - t413 * t401;
t391 = t435 * t392 + t437 * t393;
t403 = -t413 * qJD(4) + t445 * qJDD(2);
t409 = qJD(4) * mrSges(5,1) - t413 * mrSges(5,3);
t388 = m(5) * t391 - qJDD(4) * mrSges(5,2) + t403 * mrSges(5,3) - qJD(4) * t409 + t412 * t401;
t379 = t437 * t387 + t435 * t388;
t394 = -t432 * t407 + t458;
t444 = mrSges(4,3) * qJDD(2) + t439 * (-mrSges(4,1) * t434 + t462);
t377 = m(4) * t394 - t444 * t432 + t379;
t451 = -t435 * t387 + t437 * t388;
t378 = m(4) * t395 + t444 * t434 + t451;
t375 = -t432 * t377 + t434 * t378;
t371 = m(3) * t411 - t439 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t375;
t410 = -t436 * t421 + t438 * t431;
t447 = qJDD(3) - t410;
t406 = -qJDD(2) * pkin(2) - t439 * qJ(3) + t447;
t396 = (-pkin(2) - t463) * qJDD(2) + (-qJ(3) + (-t430 - t464) * pkin(5)) * t439 + t447;
t443 = m(5) * t396 - t403 * mrSges(5,1) + t404 * mrSges(5,2) - t412 * t408 + t413 * t409;
t442 = -m(4) * t406 + mrSges(4,1) * t455 - t443 + (t439 * t464 + t460) * mrSges(4,3);
t383 = m(3) * t410 - t439 * mrSges(3,2) + (mrSges(3,1) - t462) * qJDD(2) + t442;
t452 = t438 * t371 - t436 * t383;
t366 = m(2) * t421 + t452;
t374 = t434 * t377 + t432 * t378;
t373 = (m(2) + m(3)) * t420 - t374;
t459 = t433 * t366 + t461 * t373;
t367 = t436 * t371 + t438 * t383;
t448 = Ifges(4,5) * t432 + Ifges(4,6) * t434;
t457 = t439 * t448;
t454 = m(2) * t431 + t367;
t453 = t461 * t366 - t433 * t373;
t450 = Ifges(4,1) * t432 + Ifges(4,4) * t434;
t449 = Ifges(4,4) * t432 + Ifges(4,2) * t434;
t397 = Ifges(5,5) * t413 + Ifges(5,6) * t412 + Ifges(5,3) * qJD(4);
t399 = Ifges(5,1) * t413 + Ifges(5,4) * t412 + Ifges(5,5) * qJD(4);
t380 = -mrSges(5,1) * t396 + mrSges(5,3) * t391 + Ifges(5,4) * t404 + Ifges(5,2) * t403 + Ifges(5,6) * qJDD(4) + qJD(4) * t399 - t413 * t397;
t398 = Ifges(5,4) * t413 + Ifges(5,2) * t412 + Ifges(5,6) * qJD(4);
t381 = mrSges(5,2) * t396 - mrSges(5,3) * t390 + Ifges(5,1) * t404 + Ifges(5,4) * t403 + Ifges(5,5) * qJDD(4) - qJD(4) * t398 + t412 * t397;
t363 = -mrSges(4,1) * t406 + mrSges(4,3) * t395 - pkin(3) * t443 + pkin(5) * t451 + t449 * qJDD(2) + t437 * t380 + t435 * t381 - t432 * t457;
t368 = mrSges(4,2) * t406 - mrSges(4,3) * t394 - pkin(5) * t379 + t450 * qJDD(2) - t435 * t380 + t437 * t381 + t434 * t457;
t389 = qJDD(2) * t462 - t442;
t441 = mrSges(3,1) * t410 - mrSges(3,2) * t411 + Ifges(3,3) * qJDD(2) - pkin(2) * t389 + qJ(3) * t375 + t434 * t363 + t432 * t368;
t440 = mrSges(5,1) * t390 - mrSges(5,2) * t391 + Ifges(5,5) * t404 + Ifges(5,6) * t403 + Ifges(5,3) * qJDD(4) + t413 * t398 - t412 * t399;
t362 = mrSges(3,1) * t420 - mrSges(4,1) * t394 + mrSges(4,2) * t395 + mrSges(3,3) * t411 - pkin(2) * t374 - pkin(3) * t379 + (Ifges(3,6) - t448) * qJDD(2) - t440 + (-t432 * t449 + t434 * t450 + Ifges(3,5)) * t439;
t361 = -mrSges(3,2) * t420 - mrSges(3,3) * t410 + Ifges(3,5) * qJDD(2) - t439 * Ifges(3,6) - qJ(3) * t374 - t432 * t363 + t434 * t368;
t360 = -mrSges(2,1) * t431 + mrSges(2,3) * t421 - pkin(1) * t367 - t441;
t359 = mrSges(2,2) * t431 - mrSges(2,3) * t420 - pkin(4) * t367 + t438 * t361 - t436 * t362;
t1 = [-m(1) * g(1) + t453; -m(1) * g(2) + t459; -m(1) * g(3) + t454; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t459 + t461 * t359 - t433 * t360; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t453 + t433 * t359 + t461 * t360; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t420 - mrSges(2,2) * t421 + t436 * t361 + t438 * t362 + pkin(1) * (m(3) * t420 - t374) + pkin(4) * t452; t454; t441; t389; t440;];
tauJB = t1;
