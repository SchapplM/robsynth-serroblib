% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4RPRR8_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:10
% EndTime: 2019-12-31 16:55:11
% DurationCPUTime: 0.86s
% Computational Cost: add. (6751->191), mult. (12998->236), div. (0->0), fcn. (6984->6), ass. (0->80)
t467 = sin(qJ(1));
t470 = cos(qJ(1));
t450 = -t470 * g(1) - t467 * g(2);
t480 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t450;
t494 = -pkin(1) - pkin(5);
t493 = mrSges(2,1) - mrSges(3,2);
t492 = Ifges(2,5) - Ifges(3,4);
t491 = (-Ifges(2,6) + Ifges(3,5));
t449 = t467 * g(1) - t470 * g(2);
t471 = qJD(1) ^ 2;
t479 = -t471 * qJ(2) + qJDD(2) - t449;
t426 = t494 * qJDD(1) + t479;
t466 = sin(qJ(3));
t469 = cos(qJ(3));
t418 = t466 * g(3) + t469 * t426;
t487 = qJD(1) * qJD(3);
t485 = t466 * t487;
t444 = t469 * qJDD(1) - t485;
t402 = (-t444 - t485) * pkin(6) + (-t466 * t469 * t471 + qJDD(3)) * pkin(3) + t418;
t419 = -t469 * g(3) + t466 * t426;
t443 = -t466 * qJDD(1) - t469 * t487;
t488 = qJD(1) * t469;
t448 = qJD(3) * pkin(3) - pkin(6) * t488;
t462 = t466 ^ 2;
t403 = -t462 * t471 * pkin(3) + t443 * pkin(6) - qJD(3) * t448 + t419;
t465 = sin(qJ(4));
t468 = cos(qJ(4));
t400 = t468 * t402 - t465 * t403;
t435 = (-t465 * t469 - t466 * t468) * qJD(1);
t411 = t435 * qJD(4) + t465 * t443 + t468 * t444;
t436 = (-t465 * t466 + t468 * t469) * qJD(1);
t416 = -t435 * mrSges(5,1) + t436 * mrSges(5,2);
t456 = qJD(3) + qJD(4);
t423 = -t456 * mrSges(5,2) + t435 * mrSges(5,3);
t455 = qJDD(3) + qJDD(4);
t397 = m(5) * t400 + t455 * mrSges(5,1) - t411 * mrSges(5,3) - t436 * t416 + t456 * t423;
t401 = t465 * t402 + t468 * t403;
t410 = -t436 * qJD(4) + t468 * t443 - t465 * t444;
t424 = t456 * mrSges(5,1) - t436 * mrSges(5,3);
t398 = m(5) * t401 - t455 * mrSges(5,2) + t410 * mrSges(5,3) + t435 * t416 - t456 * t424;
t387 = t468 * t397 + t465 * t398;
t442 = (mrSges(4,1) * t466 + mrSges(4,2) * t469) * qJD(1);
t489 = qJD(1) * t466;
t446 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t489;
t384 = m(4) * t418 + qJDD(3) * mrSges(4,1) - t444 * mrSges(4,3) + qJD(3) * t446 - t442 * t488 + t387;
t447 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t488;
t482 = -t465 * t397 + t468 * t398;
t385 = m(4) * t419 - qJDD(3) * mrSges(4,2) + t443 * mrSges(4,3) - qJD(3) * t447 - t442 * t489 + t482;
t382 = t469 * t384 + t466 * t385;
t431 = -qJDD(1) * pkin(1) + t479;
t477 = -m(3) * t431 + (t471 * mrSges(3,3)) - t382;
t378 = m(2) * t449 - (t471 * mrSges(2,2)) + t493 * qJDD(1) + t477;
t429 = t471 * pkin(1) - t480;
t425 = t494 * t471 + t480;
t406 = t448 * t488 - t443 * pkin(3) + (-pkin(6) * t462 + t494) * t471 + t480;
t478 = m(5) * t406 - t410 * mrSges(5,1) + t411 * mrSges(5,2) - t435 * t423 + t436 * t424;
t475 = -m(4) * t425 + t443 * mrSges(4,1) - t444 * mrSges(4,2) - t446 * t489 - t447 * t488 - t478;
t473 = -m(3) * t429 + (t471 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t475;
t392 = m(2) * t450 - (t471 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t473;
t490 = t470 * t378 + t467 * t392;
t484 = -t467 * t378 + t470 * t392;
t483 = -t466 * t384 + t469 * t385;
t413 = Ifges(5,4) * t436 + Ifges(5,2) * t435 + Ifges(5,6) * t456;
t414 = Ifges(5,1) * t436 + Ifges(5,4) * t435 + Ifges(5,5) * t456;
t476 = mrSges(5,1) * t400 - mrSges(5,2) * t401 + Ifges(5,5) * t411 + Ifges(5,6) * t410 + Ifges(5,3) * t455 + t436 * t413 - t435 * t414;
t412 = Ifges(5,5) * t436 + Ifges(5,6) * t435 + Ifges(5,3) * t456;
t388 = -mrSges(5,1) * t406 + mrSges(5,3) * t401 + Ifges(5,4) * t411 + Ifges(5,2) * t410 + Ifges(5,6) * t455 - t436 * t412 + t456 * t414;
t389 = mrSges(5,2) * t406 - mrSges(5,3) * t400 + Ifges(5,1) * t411 + Ifges(5,4) * t410 + Ifges(5,5) * t455 + t435 * t412 - t456 * t413;
t432 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t469 - Ifges(4,6) * t466) * qJD(1);
t434 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t469 - Ifges(4,4) * t466) * qJD(1);
t374 = -mrSges(4,1) * t425 + mrSges(4,3) * t419 + Ifges(4,4) * t444 + Ifges(4,2) * t443 + Ifges(4,6) * qJDD(3) - pkin(3) * t478 + pkin(6) * t482 + qJD(3) * t434 + t468 * t388 + t465 * t389 - t432 * t488;
t433 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t469 - Ifges(4,2) * t466) * qJD(1);
t376 = mrSges(4,2) * t425 - mrSges(4,3) * t418 + Ifges(4,1) * t444 + Ifges(4,4) * t443 + Ifges(4,5) * qJDD(3) - pkin(6) * t387 - qJD(3) * t433 - t465 * t388 + t468 * t389 - t432 * t489;
t380 = qJDD(1) * mrSges(3,2) - t477;
t474 = mrSges(2,1) * t449 - mrSges(2,2) * t450 + mrSges(3,2) * t431 - mrSges(3,3) * t429 - pkin(1) * t380 - pkin(5) * t382 + qJ(2) * t473 - t466 * t374 + t469 * t376 + (Ifges(3,1) + Ifges(2,3)) * qJDD(1);
t472 = mrSges(4,1) * t418 - mrSges(4,2) * t419 + Ifges(4,5) * t444 + Ifges(4,6) * t443 + Ifges(4,3) * qJDD(3) + pkin(3) * t387 + t433 * t488 + t434 * t489 + t476;
t381 = -m(3) * g(3) + t483;
t373 = pkin(2) * t382 + mrSges(3,1) * t431 - mrSges(2,3) * t449 - qJ(2) * t381 + t472 + t492 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t491 * t471);
t372 = -mrSges(3,1) * t429 + mrSges(2,3) * t450 - pkin(1) * t381 - pkin(2) * t475 - pkin(5) * t483 + t493 * g(3) - t491 * qJDD(1) - t469 * t374 - t466 * t376 + t492 * t471;
t1 = [-m(1) * g(1) + t484; -m(1) * g(2) + t490; (-m(1) - m(2) - m(3)) * g(3) + t483; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t490 - t467 * t372 + t470 * t373; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t484 + t470 * t372 + t467 * t373; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t474; t474; t380; t472; t476;];
tauJB = t1;
