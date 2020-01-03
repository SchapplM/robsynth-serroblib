% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR3
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR3_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR3_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:49:11
% EndTime: 2019-12-31 16:49:12
% DurationCPUTime: 1.24s
% Computational Cost: add. (12121->197), mult. (23436->253), div. (0->0), fcn. (13204->8), ass. (0->83)
t474 = sin(qJ(1));
t477 = cos(qJ(1));
t457 = t474 * g(1) - t477 * g(2);
t448 = qJDD(1) * pkin(1) + t457;
t458 = -t477 * g(1) - t474 * g(2);
t478 = qJD(1) ^ 2;
t450 = -t478 * pkin(1) + t458;
t470 = sin(pkin(7));
t471 = cos(pkin(7));
t435 = t470 * t448 + t471 * t450;
t431 = -t478 * pkin(2) + qJDD(1) * pkin(5) + t435;
t469 = -g(3) + qJDD(2);
t473 = sin(qJ(3));
t476 = cos(qJ(3));
t421 = -t473 * t431 + t476 * t469;
t491 = qJD(1) * qJD(3);
t490 = t476 * t491;
t451 = t473 * qJDD(1) + t490;
t414 = (-t451 + t490) * pkin(6) + (t473 * t476 * t478 + qJDD(3)) * pkin(3) + t421;
t422 = t476 * t431 + t473 * t469;
t452 = t476 * qJDD(1) - t473 * t491;
t493 = qJD(1) * t473;
t456 = qJD(3) * pkin(3) - pkin(6) * t493;
t468 = t476 ^ 2;
t415 = -t468 * t478 * pkin(3) + t452 * pkin(6) - qJD(3) * t456 + t422;
t472 = sin(qJ(4));
t475 = cos(qJ(4));
t412 = t475 * t414 - t472 * t415;
t444 = (-t472 * t473 + t475 * t476) * qJD(1);
t424 = t444 * qJD(4) + t475 * t451 + t472 * t452;
t445 = (t472 * t476 + t473 * t475) * qJD(1);
t432 = -t444 * mrSges(5,1) + t445 * mrSges(5,2);
t465 = qJD(3) + qJD(4);
t436 = -t465 * mrSges(5,2) + t444 * mrSges(5,3);
t464 = qJDD(3) + qJDD(4);
t409 = m(5) * t412 + t464 * mrSges(5,1) - t424 * mrSges(5,3) - t445 * t432 + t465 * t436;
t413 = t472 * t414 + t475 * t415;
t423 = -t445 * qJD(4) - t472 * t451 + t475 * t452;
t437 = t465 * mrSges(5,1) - t445 * mrSges(5,3);
t410 = m(5) * t413 - t464 * mrSges(5,2) + t423 * mrSges(5,3) + t444 * t432 - t465 * t437;
t400 = t475 * t409 + t472 * t410;
t442 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t473 + Ifges(4,2) * t476) * qJD(1);
t443 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t473 + Ifges(4,4) * t476) * qJD(1);
t427 = Ifges(5,4) * t445 + Ifges(5,2) * t444 + Ifges(5,6) * t465;
t428 = Ifges(5,1) * t445 + Ifges(5,4) * t444 + Ifges(5,5) * t465;
t482 = -mrSges(5,1) * t412 + mrSges(5,2) * t413 - Ifges(5,5) * t424 - Ifges(5,6) * t423 - Ifges(5,3) * t464 - t445 * t427 + t444 * t428;
t495 = mrSges(4,1) * t421 - mrSges(4,2) * t422 + Ifges(4,5) * t451 + Ifges(4,6) * t452 + Ifges(4,3) * qJDD(3) + pkin(3) * t400 + (t473 * t442 - t476 * t443) * qJD(1) - t482;
t449 = (-mrSges(4,1) * t476 + mrSges(4,2) * t473) * qJD(1);
t492 = qJD(1) * t476;
t455 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t492;
t398 = m(4) * t421 + qJDD(3) * mrSges(4,1) - t451 * mrSges(4,3) + qJD(3) * t455 - t449 * t493 + t400;
t454 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t493;
t486 = -t472 * t409 + t475 * t410;
t399 = m(4) * t422 - qJDD(3) * mrSges(4,2) + t452 * mrSges(4,3) - qJD(3) * t454 + t449 * t492 + t486;
t487 = -t473 * t398 + t476 * t399;
t391 = m(3) * t435 - t478 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t487;
t434 = t471 * t448 - t470 * t450;
t484 = -qJDD(1) * pkin(2) - t434;
t430 = -t478 * pkin(5) + t484;
t416 = t456 * t493 - t452 * pkin(3) + (-pkin(6) * t468 - pkin(5)) * t478 + t484;
t483 = m(5) * t416 - t423 * mrSges(5,1) + t424 * mrSges(5,2) - t444 * t436 + t445 * t437;
t480 = -m(4) * t430 + t452 * mrSges(4,1) - t451 * mrSges(4,2) - t454 * t493 + t455 * t492 - t483;
t404 = m(3) * t434 + qJDD(1) * mrSges(3,1) - t478 * mrSges(3,2) + t480;
t386 = t470 * t391 + t471 * t404;
t383 = m(2) * t457 + qJDD(1) * mrSges(2,1) - t478 * mrSges(2,2) + t386;
t488 = t471 * t391 - t470 * t404;
t384 = m(2) * t458 - t478 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t488;
t494 = t477 * t383 + t474 * t384;
t394 = t476 * t398 + t473 * t399;
t392 = m(3) * t469 + t394;
t489 = -t474 * t383 + t477 * t384;
t426 = Ifges(5,5) * t445 + Ifges(5,6) * t444 + Ifges(5,3) * t465;
t401 = -mrSges(5,1) * t416 + mrSges(5,3) * t413 + Ifges(5,4) * t424 + Ifges(5,2) * t423 + Ifges(5,6) * t464 - t445 * t426 + t465 * t428;
t402 = mrSges(5,2) * t416 - mrSges(5,3) * t412 + Ifges(5,1) * t424 + Ifges(5,4) * t423 + Ifges(5,5) * t464 + t444 * t426 - t465 * t427;
t441 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t473 + Ifges(4,6) * t476) * qJD(1);
t379 = -mrSges(4,1) * t430 + mrSges(4,3) * t422 + Ifges(4,4) * t451 + Ifges(4,2) * t452 + Ifges(4,6) * qJDD(3) - pkin(3) * t483 + pkin(6) * t486 + qJD(3) * t443 + t475 * t401 + t472 * t402 - t441 * t493;
t388 = mrSges(4,2) * t430 - mrSges(4,3) * t421 + Ifges(4,1) * t451 + Ifges(4,4) * t452 + Ifges(4,5) * qJDD(3) - pkin(6) * t400 - qJD(3) * t442 - t472 * t401 + t475 * t402 + t441 * t492;
t481 = mrSges(2,1) * t457 + mrSges(3,1) * t434 - mrSges(2,2) * t458 - mrSges(3,2) * t435 + pkin(1) * t386 + pkin(2) * t480 + pkin(5) * t487 + t476 * t379 + t473 * t388 + (Ifges(2,3) + Ifges(3,3)) * qJDD(1);
t377 = -mrSges(3,1) * t469 + mrSges(3,3) * t435 + t478 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t394 - t495;
t376 = mrSges(3,2) * t469 - mrSges(3,3) * t434 + Ifges(3,5) * qJDD(1) - t478 * Ifges(3,6) - pkin(5) * t394 - t473 * t379 + t476 * t388;
t375 = -mrSges(2,2) * g(3) - mrSges(2,3) * t457 + Ifges(2,5) * qJDD(1) - t478 * Ifges(2,6) - qJ(2) * t386 + t471 * t376 - t470 * t377;
t374 = mrSges(2,1) * g(3) + mrSges(2,3) * t458 + t478 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t392 + qJ(2) * t488 + t470 * t376 + t471 * t377;
t1 = [-m(1) * g(1) + t489; -m(1) * g(2) + t494; (-m(1) - m(2)) * g(3) + t392; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t494 - t474 * t374 + t477 * t375; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t489 + t477 * t374 + t474 * t375; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t481; t481; t392; t495; -t482;];
tauJB = t1;
