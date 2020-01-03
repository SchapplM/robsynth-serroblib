% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRP4
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRP4_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:31
% EndTime: 2019-12-31 17:34:32
% DurationCPUTime: 0.78s
% Computational Cost: add. (5330->178), mult. (9693->210), div. (0->0), fcn. (4790->6), ass. (0->80)
t496 = Ifges(5,1) + Ifges(6,1);
t487 = Ifges(5,4) + Ifges(6,4);
t486 = Ifges(5,5) + Ifges(6,5);
t495 = Ifges(5,2) + Ifges(6,2);
t485 = Ifges(5,6) + Ifges(6,6);
t494 = Ifges(5,3) + Ifges(6,3);
t459 = sin(qJ(4));
t461 = cos(qJ(4));
t436 = (-mrSges(6,1) * t461 + mrSges(6,2) * t459) * qJD(3);
t476 = qJD(3) * qJD(4);
t472 = t461 * t476;
t438 = t459 * qJDD(3) + t472;
t457 = sin(pkin(7));
t458 = cos(pkin(7));
t443 = t457 * g(1) - t458 * g(2);
t441 = qJDD(2) - t443;
t444 = -t458 * g(1) - t457 * g(2);
t460 = sin(qJ(3));
t462 = cos(qJ(3));
t420 = t460 * t441 + t462 * t444;
t463 = qJD(3) ^ 2;
t418 = -t463 * pkin(3) + qJDD(3) * pkin(6) + t420;
t456 = g(3) - qJDD(1);
t451 = t461 * t456;
t475 = qJD(3) * qJD(5);
t490 = pkin(4) * t463;
t411 = qJDD(4) * pkin(4) + t451 + (-t438 + t472) * qJ(5) + (t461 * t490 - t418 - 0.2e1 * t475) * t459;
t477 = qJD(3) * t461;
t448 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t477;
t474 = m(6) * t411 + qJDD(4) * mrSges(6,1) + qJD(4) * t448;
t478 = qJD(3) * t459;
t407 = -t438 * mrSges(6,3) - t436 * t478 + t474;
t415 = t461 * t418 + t459 * t456;
t439 = t461 * qJDD(3) - t459 * t476;
t445 = qJD(4) * pkin(4) - qJ(5) * t478;
t455 = t461 ^ 2;
t412 = t439 * qJ(5) - qJD(4) * t445 - t455 * t490 + 0.2e1 * t461 * t475 + t415;
t414 = -t459 * t418 + t451;
t480 = t486 * qJD(4) + (t496 * t459 + t487 * t461) * qJD(3);
t481 = t485 * qJD(4) + (t487 * t459 + t495 * t461) * qJD(3);
t493 = mrSges(5,1) * t414 + mrSges(6,1) * t411 - mrSges(5,2) * t415 - mrSges(6,2) * t412 + pkin(4) * t407 + (t481 * t459 - t480 * t461) * qJD(3) + t494 * qJDD(4) + t486 * t438 + t485 * t439;
t489 = -mrSges(2,2) + mrSges(3,3);
t488 = -mrSges(5,2) - mrSges(6,2);
t437 = (-mrSges(5,1) * t461 + mrSges(5,2) * t459) * qJD(3);
t449 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t477;
t404 = m(5) * t414 + qJDD(4) * mrSges(5,1) + qJD(4) * t449 + (-mrSges(5,3) - mrSges(6,3)) * t438 + (-t436 - t437) * t478 + t474;
t473 = m(6) * t412 + t439 * mrSges(6,3) + t436 * t477;
t446 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t478;
t479 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t478 - t446;
t405 = m(5) * t415 + t439 * mrSges(5,3) + t479 * qJD(4) + t488 * qJDD(4) + t437 * t477 + t473;
t401 = -t459 * t404 + t461 * t405;
t398 = m(4) * t420 - t463 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t401;
t419 = t462 * t441 - t460 * t444;
t467 = -qJDD(3) * pkin(3) - t419;
t417 = -t463 * pkin(6) + t467;
t413 = t445 * t478 - t439 * pkin(4) + qJDD(5) + (-qJ(5) * t455 - pkin(6)) * t463 + t467;
t469 = -m(6) * t413 + t439 * mrSges(6,1) + t448 * t477;
t406 = -m(5) * t417 + t439 * mrSges(5,1) + t488 * t438 + t449 * t477 + t479 * t478 + t469;
t402 = m(4) * t419 + qJDD(3) * mrSges(4,1) - t463 * mrSges(4,2) + t406;
t394 = t460 * t398 + t462 * t402;
t393 = m(3) * t441 + t394;
t391 = m(2) * t443 - t393;
t470 = t462 * t398 - t460 * t402;
t468 = m(3) * t444 + t470;
t392 = m(2) * t444 + t468;
t483 = t458 * t391 + t457 * t392;
t482 = t494 * qJD(4) + (t486 * t459 + t485 * t461) * qJD(3);
t471 = -t457 * t391 + t458 * t392;
t400 = t461 * t404 + t459 * t405;
t399 = -t400 + (-m(3) - m(4)) * t456;
t466 = -m(2) * t456 + t399;
t408 = t438 * mrSges(6,2) + t446 * t478 - t469;
t395 = -mrSges(5,1) * t417 + mrSges(5,3) * t415 - mrSges(6,1) * t413 + mrSges(6,3) * t412 - pkin(4) * t408 + qJ(5) * t473 + t495 * t439 + t487 * t438 + (-qJ(5) * mrSges(6,2) + t485) * qJDD(4) + (-qJ(5) * t446 + t480) * qJD(4) - t482 * t478;
t396 = mrSges(5,2) * t417 + mrSges(6,2) * t413 - mrSges(5,3) * t414 - mrSges(6,3) * t411 - qJ(5) * t407 - t481 * qJD(4) + t486 * qJDD(4) + t496 * t438 + t487 * t439 + t482 * t477;
t464 = mrSges(4,1) * t419 - mrSges(4,2) * t420 + Ifges(4,3) * qJDD(3) + pkin(3) * t406 + pkin(6) * t401 + t461 * t395 + t459 * t396;
t387 = -mrSges(4,1) * t456 + mrSges(4,3) * t420 + t463 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t400 - t493;
t386 = mrSges(4,2) * t456 - mrSges(4,3) * t419 + Ifges(4,5) * qJDD(3) - t463 * Ifges(4,6) - pkin(6) * t400 - t459 * t395 + t461 * t396;
t385 = mrSges(3,2) * t441 - mrSges(2,3) * t443 - pkin(5) * t394 - qJ(2) * t399 + t462 * t386 - t460 * t387 + t489 * t456;
t384 = -t460 * t386 - t462 * t387 + pkin(2) * t400 - pkin(5) * t470 - pkin(1) * t399 + (pkin(2) * m(4) + mrSges(2,1) + mrSges(3,1)) * t456 + (mrSges(3,2) + mrSges(2,3)) * t444;
t1 = [-m(1) * g(1) + t471; -m(1) * g(2) + t483; -m(1) * g(3) + t466; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t483 - t457 * t384 + t458 * t385; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t471 + t458 * t384 + t457 * t385; -mrSges(1,1) * g(2) + mrSges(2,1) * t443 - mrSges(3,1) * t441 + mrSges(1,2) * g(1) - pkin(1) * t393 - pkin(2) * t394 + qJ(2) * t468 + t489 * t444 - t464; t466; t393; t464; t493; t408;];
tauJB = t1;
