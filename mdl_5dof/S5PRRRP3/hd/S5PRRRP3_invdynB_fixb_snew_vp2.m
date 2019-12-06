% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP3
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:46
% EndTime: 2019-12-05 16:43:50
% DurationCPUTime: 2.36s
% Computational Cost: add. (22202->237), mult. (44961->292), div. (0->0), fcn. (28134->8), ass. (0->96)
t525 = Ifges(5,1) + Ifges(6,1);
t522 = Ifges(5,4) + Ifges(6,4);
t521 = Ifges(5,5) + Ifges(6,5);
t524 = Ifges(5,2) + Ifges(6,2);
t520 = -Ifges(5,6) - Ifges(6,6);
t523 = -Ifges(5,3) - Ifges(6,3);
t492 = sin(pkin(8));
t493 = cos(pkin(8));
t479 = t492 * g(1) - t493 * g(2);
t480 = -t493 * g(1) - t492 * g(2);
t496 = sin(qJ(2));
t499 = cos(qJ(2));
t458 = t496 * t479 + t499 * t480;
t500 = qJD(2) ^ 2;
t455 = -t500 * pkin(2) + qJDD(2) * pkin(6) + t458;
t491 = -g(3) + qJDD(1);
t495 = sin(qJ(3));
t498 = cos(qJ(3));
t440 = -t495 * t455 + t498 * t491;
t513 = qJD(2) * qJD(3);
t509 = t498 * t513;
t477 = t495 * qJDD(2) + t509;
t430 = (-t477 + t509) * pkin(7) + (t495 * t498 * t500 + qJDD(3)) * pkin(3) + t440;
t441 = t498 * t455 + t495 * t491;
t478 = t498 * qJDD(2) - t495 * t513;
t515 = qJD(2) * t495;
t483 = qJD(3) * pkin(3) - pkin(7) * t515;
t490 = t498 ^ 2;
t431 = -t490 * t500 * pkin(3) + t478 * pkin(7) - qJD(3) * t483 + t441;
t494 = sin(qJ(4));
t497 = cos(qJ(4));
t425 = t497 * t430 - t494 * t431;
t468 = (-t494 * t495 + t497 * t498) * qJD(2);
t439 = t468 * qJD(4) + t497 * t477 + t494 * t478;
t469 = (t494 * t498 + t495 * t497) * qJD(2);
t452 = -t468 * mrSges(6,1) + t469 * mrSges(6,2);
t453 = -t468 * mrSges(5,1) + t469 * mrSges(5,2);
t489 = qJD(3) + qJD(4);
t460 = -t489 * mrSges(5,2) + t468 * mrSges(5,3);
t488 = qJDD(3) + qJDD(4);
t420 = -0.2e1 * qJD(5) * t469 + (t468 * t489 - t439) * qJ(5) + (t468 * t469 + t488) * pkin(4) + t425;
t459 = -t489 * mrSges(6,2) + t468 * mrSges(6,3);
t512 = m(6) * t420 + t488 * mrSges(6,1) + t489 * t459;
t414 = m(5) * t425 + t488 * mrSges(5,1) + t489 * t460 + (-t452 - t453) * t469 + (-mrSges(5,3) - mrSges(6,3)) * t439 + t512;
t426 = t494 * t430 + t497 * t431;
t438 = -t469 * qJD(4) - t494 * t477 + t497 * t478;
t462 = t489 * mrSges(6,1) - t469 * mrSges(6,3);
t463 = t489 * mrSges(5,1) - t469 * mrSges(5,3);
t461 = t489 * pkin(4) - t469 * qJ(5);
t464 = t468 ^ 2;
t422 = -t464 * pkin(4) + t438 * qJ(5) + 0.2e1 * qJD(5) * t468 - t489 * t461 + t426;
t511 = m(6) * t422 + t438 * mrSges(6,3) + t468 * t452;
t417 = m(5) * t426 + t438 * mrSges(5,3) + t468 * t453 + (-t462 - t463) * t489 + (-mrSges(5,2) - mrSges(6,2)) * t488 + t511;
t410 = t497 * t414 + t494 * t417;
t476 = (-mrSges(4,1) * t498 + mrSges(4,2) * t495) * qJD(2);
t514 = qJD(2) * t498;
t482 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t514;
t408 = m(4) * t440 + qJDD(3) * mrSges(4,1) - t477 * mrSges(4,3) + qJD(3) * t482 - t476 * t515 + t410;
t481 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t515;
t505 = -t494 * t414 + t497 * t417;
t409 = m(4) * t441 - qJDD(3) * mrSges(4,2) + t478 * mrSges(4,3) - qJD(3) * t481 + t476 * t514 + t505;
t506 = -t495 * t408 + t498 * t409;
t401 = m(3) * t458 - t500 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t506;
t457 = t499 * t479 - t496 * t480;
t503 = -qJDD(2) * pkin(2) - t457;
t454 = -t500 * pkin(6) + t503;
t432 = -t478 * pkin(3) + t483 * t515 + (-pkin(7) * t490 - pkin(6)) * t500 + t503;
t424 = -t438 * pkin(4) - t464 * qJ(5) + t469 * t461 + qJDD(5) + t432;
t504 = m(6) * t424 - t438 * mrSges(6,1) + t439 * mrSges(6,2) - t468 * t459 + t469 * t462;
t502 = m(5) * t432 - t438 * mrSges(5,1) + t439 * mrSges(5,2) - t468 * t460 + t469 * t463 + t504;
t501 = -m(4) * t454 + t478 * mrSges(4,1) - t477 * mrSges(4,2) - t481 * t515 + t482 * t514 - t502;
t412 = m(3) * t457 + qJDD(2) * mrSges(3,1) - t500 * mrSges(3,2) + t501;
t398 = t496 * t401 + t499 * t412;
t396 = m(2) * t479 + t398;
t507 = t499 * t401 - t496 * t412;
t397 = m(2) * t480 + t507;
t519 = t493 * t396 + t492 * t397;
t402 = t498 * t408 + t495 * t409;
t518 = t520 * t468 - t521 * t469 + t523 * t489;
t517 = -t524 * t468 - t522 * t469 + t520 * t489;
t516 = t522 * t468 + t525 * t469 + t521 * t489;
t510 = m(3) * t491 + t402;
t508 = -t492 * t396 + t493 * t397;
t467 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t495 + Ifges(4,4) * t498) * qJD(2);
t466 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t495 + Ifges(4,2) * t498) * qJD(2);
t465 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t495 + Ifges(4,6) * t498) * qJD(2);
t418 = -t439 * mrSges(6,3) - t469 * t452 + t512;
t404 = mrSges(5,2) * t432 + mrSges(6,2) * t424 - mrSges(5,3) * t425 - mrSges(6,3) * t420 - qJ(5) * t418 + t522 * t438 + t525 * t439 - t518 * t468 + t521 * t488 + t517 * t489;
t403 = -mrSges(5,1) * t432 + mrSges(5,3) * t426 - mrSges(6,1) * t424 + mrSges(6,3) * t422 - pkin(4) * t504 + qJ(5) * t511 + (-qJ(5) * t462 + t516) * t489 + (-qJ(5) * mrSges(6,2) - t520) * t488 + t518 * t469 + t522 * t439 + t524 * t438;
t392 = mrSges(4,2) * t454 - mrSges(4,3) * t440 + Ifges(4,1) * t477 + Ifges(4,4) * t478 + Ifges(4,5) * qJDD(3) - pkin(7) * t410 - qJD(3) * t466 - t494 * t403 + t497 * t404 + t465 * t514;
t391 = -mrSges(4,1) * t454 + mrSges(4,3) * t441 + Ifges(4,4) * t477 + Ifges(4,2) * t478 + Ifges(4,6) * qJDD(3) - pkin(3) * t502 + pkin(7) * t505 + qJD(3) * t467 + t497 * t403 + t494 * t404 - t465 * t515;
t390 = Ifges(3,6) * qJDD(2) - Ifges(4,3) * qJDD(3) + t500 * Ifges(3,5) - mrSges(3,1) * t491 - Ifges(4,5) * t477 - Ifges(4,6) * t478 + mrSges(3,3) * t458 - mrSges(4,1) * t440 + mrSges(4,2) * t441 - mrSges(5,1) * t425 + mrSges(5,2) * t426 - mrSges(6,1) * t420 + mrSges(6,2) * t422 - pkin(4) * t418 - pkin(3) * t410 - pkin(2) * t402 + t523 * t488 + t517 * t469 + t516 * t468 - t521 * t439 + t520 * t438 + (-t495 * t466 + t498 * t467) * qJD(2);
t389 = mrSges(3,2) * t491 - mrSges(3,3) * t457 + Ifges(3,5) * qJDD(2) - t500 * Ifges(3,6) - pkin(6) * t402 - t495 * t391 + t498 * t392;
t388 = mrSges(2,2) * t491 - mrSges(2,3) * t479 - pkin(5) * t398 + t499 * t389 - t496 * t390;
t387 = -mrSges(2,1) * t491 + mrSges(2,3) * t480 - pkin(1) * t510 + pkin(5) * t507 + t496 * t389 + t499 * t390;
t1 = [-m(1) * g(1) + t508; -m(1) * g(2) + t519; -m(1) * g(3) + m(2) * t491 + t510; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t519 - t492 * t387 + t493 * t388; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t508 + t493 * t387 + t492 * t388; -mrSges(1,1) * g(2) + mrSges(2,1) * t479 + mrSges(3,1) * t457 + mrSges(1,2) * g(1) - mrSges(2,2) * t480 - mrSges(3,2) * t458 + Ifges(3,3) * qJDD(2) + pkin(1) * t398 + pkin(2) * t501 + pkin(6) * t506 + t498 * t391 + t495 * t392;];
tauB = t1;
