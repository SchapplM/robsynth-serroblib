% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RRRP2
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RRRP2_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:58
% EndTime: 2019-12-31 17:12:59
% DurationCPUTime: 0.85s
% Computational Cost: add. (7807->177), mult. (9964->212), div. (0->0), fcn. (4470->6), ass. (0->79)
t514 = Ifges(4,1) + Ifges(5,1);
t505 = Ifges(4,4) + Ifges(5,4);
t504 = Ifges(4,5) + Ifges(5,5);
t513 = Ifges(4,2) + Ifges(5,2);
t503 = Ifges(4,6) + Ifges(5,6);
t512 = Ifges(4,3) + Ifges(5,3);
t470 = qJD(1) + qJD(2);
t474 = sin(qJ(3));
t477 = cos(qJ(3));
t448 = (-mrSges(5,1) * t477 + mrSges(5,2) * t474) * t470;
t469 = qJDD(1) + qJDD(2);
t494 = qJD(3) * t470;
t490 = t477 * t494;
t450 = t474 * t469 + t490;
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t464 = t476 * g(1) - t479 * g(2);
t456 = qJDD(1) * pkin(1) + t464;
t465 = -t479 * g(1) - t476 * g(2);
t480 = qJD(1) ^ 2;
t457 = -t480 * pkin(1) + t465;
t475 = sin(qJ(2));
t478 = cos(qJ(2));
t435 = t475 * t456 + t478 * t457;
t468 = t470 ^ 2;
t432 = -t468 * pkin(2) + t469 * pkin(6) + t435;
t493 = qJD(4) * t470;
t507 = t477 * g(3);
t508 = pkin(3) * t468;
t425 = qJDD(3) * pkin(3) - t507 + (-t450 + t490) * qJ(4) + (t477 * t508 - t432 - 0.2e1 * t493) * t474;
t500 = t470 * t477;
t461 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t500;
t492 = m(5) * t425 + qJDD(3) * mrSges(5,1) + qJD(3) * t461;
t501 = t470 * t474;
t422 = -t450 * mrSges(5,3) - t448 * t501 + t492;
t429 = -t474 * g(3) + t477 * t432;
t451 = t477 * t469 - t474 * t494;
t458 = qJD(3) * pkin(3) - qJ(4) * t501;
t473 = t477 ^ 2;
t426 = t451 * qJ(4) - qJD(3) * t458 - t473 * t508 + 0.2e1 * t477 * t493 + t429;
t428 = -t474 * t432 - t507;
t496 = (t514 * t474 + t505 * t477) * t470 + t504 * qJD(3);
t497 = (t505 * t474 + t513 * t477) * t470 + t503 * qJD(3);
t511 = mrSges(4,1) * t428 + mrSges(5,1) * t425 - mrSges(4,2) * t429 - mrSges(5,2) * t426 + pkin(3) * t422 + t512 * qJDD(3) + t504 * t450 + t503 * t451 + (t497 * t474 - t496 * t477) * t470;
t506 = -mrSges(4,2) - mrSges(5,2);
t449 = (-mrSges(4,1) * t477 + mrSges(4,2) * t474) * t470;
t462 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t500;
t419 = m(4) * t428 + qJDD(3) * mrSges(4,1) + qJD(3) * t462 + (-t448 - t449) * t501 + (-mrSges(4,3) - mrSges(5,3)) * t450 + t492;
t491 = m(5) * t426 + t451 * mrSges(5,3) + t448 * t500;
t459 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t501;
t495 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t501 - t459;
t420 = m(4) * t429 + t451 * mrSges(4,3) + t495 * qJD(3) + t506 * qJDD(3) + t449 * t500 + t491;
t487 = -t474 * t419 + t477 * t420;
t410 = m(3) * t435 - t468 * mrSges(3,1) - t469 * mrSges(3,2) + t487;
t434 = t478 * t456 - t475 * t457;
t485 = -t469 * pkin(2) - t434;
t431 = -t468 * pkin(6) + t485;
t427 = t458 * t501 - t451 * pkin(3) + qJDD(4) + (-qJ(4) * t473 - pkin(6)) * t468 + t485;
t486 = -m(5) * t427 + t451 * mrSges(5,1) + t461 * t500;
t481 = -m(4) * t431 + t451 * mrSges(4,1) + t506 * t450 + t462 * t500 + t495 * t501 + t486;
t414 = m(3) * t434 + t469 * mrSges(3,1) - t468 * mrSges(3,2) + t481;
t403 = t475 * t410 + t478 * t414;
t400 = m(2) * t464 + qJDD(1) * mrSges(2,1) - t480 * mrSges(2,2) + t403;
t488 = t478 * t410 - t475 * t414;
t401 = m(2) * t465 - t480 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t488;
t499 = t479 * t400 + t476 * t401;
t412 = t477 * t419 + t474 * t420;
t498 = (t504 * t474 + t503 * t477) * t470 + t512 * qJD(3);
t489 = -t476 * t400 + t479 * t401;
t421 = t450 * mrSges(5,2) + t459 * t501 - t486;
t405 = -mrSges(4,1) * t431 + mrSges(4,3) * t429 - mrSges(5,1) * t427 + mrSges(5,3) * t426 - pkin(3) * t421 + qJ(4) * t491 - t498 * t501 + t513 * t451 + t505 * t450 + (-qJ(4) * mrSges(5,2) + t503) * qJDD(3) + (-qJ(4) * t459 + t496) * qJD(3);
t407 = mrSges(4,2) * t431 + mrSges(5,2) * t427 - mrSges(4,3) * t428 - mrSges(5,3) * t425 - qJ(4) * t422 - t497 * qJD(3) + t504 * qJDD(3) + t514 * t450 + t505 * t451 + t498 * t500;
t484 = mrSges(3,1) * t434 - mrSges(3,2) * t435 + Ifges(3,3) * t469 + pkin(2) * t481 + pkin(6) * t487 + t477 * t405 + t474 * t407;
t482 = mrSges(2,1) * t464 - mrSges(2,2) * t465 + Ifges(2,3) * qJDD(1) + pkin(1) * t403 + t484;
t396 = mrSges(3,1) * g(3) + mrSges(3,3) * t435 + t468 * Ifges(3,5) + Ifges(3,6) * t469 - pkin(2) * t412 - t511;
t395 = -mrSges(3,2) * g(3) - mrSges(3,3) * t434 + Ifges(3,5) * t469 - t468 * Ifges(3,6) - pkin(6) * t412 - t474 * t405 + t477 * t407;
t394 = -mrSges(2,2) * g(3) - mrSges(2,3) * t464 + Ifges(2,5) * qJDD(1) - t480 * Ifges(2,6) - pkin(5) * t403 + t478 * t395 - t475 * t396;
t393 = Ifges(2,6) * qJDD(1) + t480 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t465 + t475 * t395 + t478 * t396 - pkin(1) * (-m(3) * g(3) + t412) + pkin(5) * t488;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t499; (-m(1) - m(2) - m(3)) * g(3) + t412; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t499 - t476 * t393 + t479 * t394; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t489 + t479 * t393 + t476 * t394; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t482; t482; t484; t511; t421;];
tauJB = t1;
