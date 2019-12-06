% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:51:00
% EndTime: 2019-12-05 16:51:03
% DurationCPUTime: 2.05s
% Computational Cost: add. (19651->236), mult. (38880->291), div. (0->0), fcn. (23986->8), ass. (0->95)
t515 = Ifges(5,1) + Ifges(6,1);
t510 = Ifges(5,4) - Ifges(6,5);
t509 = Ifges(6,4) + Ifges(5,5);
t514 = Ifges(5,2) + Ifges(6,3);
t513 = -Ifges(6,2) - Ifges(5,3);
t508 = Ifges(5,6) - Ifges(6,6);
t512 = cos(qJ(4));
t511 = -mrSges(5,3) - mrSges(6,2);
t507 = cos(pkin(8));
t482 = sin(pkin(8));
t468 = -t507 * g(1) - t482 * g(2);
t481 = -g(3) + qJDD(1);
t485 = sin(qJ(2));
t487 = cos(qJ(2));
t450 = t487 * t468 + t485 * t481;
t488 = qJD(2) ^ 2;
t444 = -t488 * pkin(2) + qJDD(2) * pkin(6) + t450;
t467 = t482 * g(1) - t507 * g(2);
t484 = sin(qJ(3));
t486 = cos(qJ(3));
t426 = -t484 * t444 - t486 * t467;
t499 = qJD(2) * qJD(3);
t497 = t486 * t499;
t465 = t484 * qJDD(2) + t497;
t418 = (-t465 + t497) * pkin(7) + (t484 * t486 * t488 + qJDD(3)) * pkin(3) + t426;
t427 = t486 * t444 - t484 * t467;
t466 = t486 * qJDD(2) - t484 * t499;
t501 = qJD(2) * t484;
t471 = qJD(3) * pkin(3) - pkin(7) * t501;
t480 = t486 ^ 2;
t419 = -t480 * t488 * pkin(3) + t466 * pkin(7) - qJD(3) * t471 + t427;
t483 = sin(qJ(4));
t415 = t483 * t418 + t512 * t419;
t455 = (t483 * t486 + t512 * t484) * qJD(2);
t424 = t455 * qJD(4) + t483 * t465 - t512 * t466;
t479 = qJD(3) + qJD(4);
t446 = t479 * mrSges(5,1) - t455 * mrSges(5,3);
t500 = qJD(2) * t486;
t454 = t483 * t501 - t512 * t500;
t478 = qJDD(3) + qJDD(4);
t437 = t454 * pkin(4) - t455 * qJ(5);
t477 = t479 ^ 2;
t410 = -t477 * pkin(4) + t478 * qJ(5) + 0.2e1 * qJD(5) * t479 - t454 * t437 + t415;
t447 = -t479 * mrSges(6,1) + t455 * mrSges(6,2);
t498 = m(6) * t410 + t478 * mrSges(6,3) + t479 * t447;
t438 = t454 * mrSges(6,1) - t455 * mrSges(6,3);
t502 = -t454 * mrSges(5,1) - t455 * mrSges(5,2) - t438;
t405 = m(5) * t415 - t478 * mrSges(5,2) + t511 * t424 - t479 * t446 + t502 * t454 + t498;
t414 = t512 * t418 - t483 * t419;
t425 = -t454 * qJD(4) + t512 * t465 + t483 * t466;
t445 = -t479 * mrSges(5,2) - t454 * mrSges(5,3);
t411 = -t478 * pkin(4) - t477 * qJ(5) + t455 * t437 + qJDD(5) - t414;
t448 = -t454 * mrSges(6,2) + t479 * mrSges(6,3);
t492 = -m(6) * t411 + t478 * mrSges(6,1) + t479 * t448;
t407 = m(5) * t414 + t478 * mrSges(5,1) + t511 * t425 + t479 * t445 + t502 * t455 + t492;
t400 = t483 * t405 + t512 * t407;
t464 = (-mrSges(4,1) * t486 + mrSges(4,2) * t484) * qJD(2);
t470 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t500;
t398 = m(4) * t426 + qJDD(3) * mrSges(4,1) - t465 * mrSges(4,3) + qJD(3) * t470 - t464 * t501 + t400;
t469 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t501;
t493 = t512 * t405 - t483 * t407;
t399 = m(4) * t427 - qJDD(3) * mrSges(4,2) + t466 * mrSges(4,3) - qJD(3) * t469 + t464 * t500 + t493;
t494 = -t484 * t398 + t486 * t399;
t391 = m(3) * t450 - t488 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t494;
t449 = -t485 * t468 + t487 * t481;
t491 = -qJDD(2) * pkin(2) - t449;
t443 = -t488 * pkin(6) + t491;
t420 = -t466 * pkin(3) + t471 * t501 + (-pkin(7) * t480 - pkin(6)) * t488 + t491;
t413 = -0.2e1 * qJD(5) * t455 + (t454 * t479 - t425) * qJ(5) + (t455 * t479 + t424) * pkin(4) + t420;
t408 = m(6) * t413 + t424 * mrSges(6,1) - t425 * mrSges(6,3) - t455 * t447 + t454 * t448;
t490 = m(5) * t420 + t424 * mrSges(5,1) + t425 * mrSges(5,2) + t454 * t445 + t455 * t446 + t408;
t489 = -m(4) * t443 + t466 * mrSges(4,1) - t465 * mrSges(4,2) - t469 * t501 + t470 * t500 - t490;
t402 = m(3) * t449 + qJDD(2) * mrSges(3,1) - t488 * mrSges(3,2) + t489;
t495 = t487 * t391 - t485 * t402;
t387 = m(2) * t468 + t495;
t394 = t486 * t398 + t484 * t399;
t393 = (m(2) + m(3)) * t467 - t394;
t506 = t482 * t387 + t507 * t393;
t388 = t485 * t391 + t487 * t402;
t505 = t514 * t454 - t510 * t455 - t508 * t479;
t504 = t508 * t454 - t509 * t455 + t513 * t479;
t503 = -t510 * t454 + t515 * t455 + t509 * t479;
t496 = t507 * t387 - t482 * t393;
t453 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t484 + Ifges(4,4) * t486) * qJD(2);
t452 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t484 + Ifges(4,2) * t486) * qJD(2);
t451 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t484 + Ifges(4,6) * t486) * qJD(2);
t396 = mrSges(5,2) * t420 + mrSges(6,2) * t411 - mrSges(5,3) * t414 - mrSges(6,3) * t413 - qJ(5) * t408 - t510 * t424 + t515 * t425 + t504 * t454 + t509 * t478 + t505 * t479;
t395 = -mrSges(5,1) * t420 - mrSges(6,1) * t413 + mrSges(6,2) * t410 + mrSges(5,3) * t415 - pkin(4) * t408 - t514 * t424 + t510 * t425 + t504 * t455 + t508 * t478 + t503 * t479;
t384 = mrSges(4,2) * t443 - mrSges(4,3) * t426 + Ifges(4,1) * t465 + Ifges(4,4) * t466 + Ifges(4,5) * qJDD(3) - pkin(7) * t400 - qJD(3) * t452 - t483 * t395 + t512 * t396 + t451 * t500;
t383 = -mrSges(4,1) * t443 + mrSges(4,3) * t427 + Ifges(4,4) * t465 + Ifges(4,2) * t466 + Ifges(4,6) * qJDD(3) - pkin(3) * t490 + pkin(7) * t493 + qJD(3) * t453 + t512 * t395 + t483 * t396 - t451 * t501;
t382 = -Ifges(4,3) * qJDD(3) + Ifges(3,6) * qJDD(2) - qJ(5) * t498 + t488 * Ifges(3,5) - pkin(4) * t492 - Ifges(4,5) * t465 - Ifges(4,6) * t466 + mrSges(3,1) * t467 + mrSges(3,3) * t450 - mrSges(4,1) * t426 + mrSges(4,2) * t427 - mrSges(5,1) * t414 + mrSges(5,2) * t415 - mrSges(6,3) * t410 + mrSges(6,1) * t411 - pkin(3) * t400 - pkin(2) * t394 + t513 * t478 + (pkin(4) * t438 + t505) * t455 + (qJ(5) * t438 - t503) * t454 + (pkin(4) * mrSges(6,2) - t509) * t425 + (qJ(5) * mrSges(6,2) + t508) * t424 + (-t484 * t452 + t486 * t453) * qJD(2);
t381 = -mrSges(3,2) * t467 - mrSges(3,3) * t449 + Ifges(3,5) * qJDD(2) - t488 * Ifges(3,6) - pkin(6) * t394 - t484 * t383 + t486 * t384;
t380 = -mrSges(2,1) * t481 - mrSges(3,1) * t449 + mrSges(3,2) * t450 + mrSges(2,3) * t468 - Ifges(3,3) * qJDD(2) - pkin(1) * t388 - pkin(2) * t489 - pkin(6) * t494 - t486 * t383 - t484 * t384;
t379 = mrSges(2,2) * t481 - mrSges(2,3) * t467 - pkin(5) * t388 + t487 * t381 - t485 * t382;
t1 = [-m(1) * g(1) + t496; -m(1) * g(2) + t506; -m(1) * g(3) + m(2) * t481 + t388; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t506 + t507 * t379 - t482 * t380; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t496 + t482 * t379 + t507 * t380; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t467 - mrSges(2,2) * t468 + t485 * t381 + t487 * t382 + pkin(1) * (m(3) * t467 - t394) + pkin(5) * t495;];
tauB = t1;
