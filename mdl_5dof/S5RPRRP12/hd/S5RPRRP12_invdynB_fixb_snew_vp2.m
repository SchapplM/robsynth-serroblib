% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP12
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:56:25
% EndTime: 2019-12-31 18:56:28
% DurationCPUTime: 1.49s
% Computational Cost: add. (10965->242), mult. (20735->280), div. (0->0), fcn. (11269->6), ass. (0->96)
t517 = Ifges(5,1) + Ifges(6,1);
t508 = Ifges(5,4) + Ifges(6,4);
t506 = Ifges(5,5) + Ifges(6,5);
t516 = Ifges(5,2) + Ifges(6,2);
t515 = Ifges(5,6) + Ifges(6,6);
t514 = Ifges(5,3) + Ifges(6,3);
t476 = sin(qJ(1));
t479 = cos(qJ(1));
t466 = -t479 * g(1) - t476 * g(2);
t513 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t466;
t512 = -pkin(1) - pkin(6);
t511 = mrSges(2,1) - mrSges(3,2);
t510 = -mrSges(5,2) - mrSges(6,2);
t509 = -Ifges(3,4) + Ifges(2,5);
t507 = (Ifges(3,5) - Ifges(2,6));
t465 = t476 * g(1) - t479 * g(2);
t481 = qJD(1) ^ 2;
t486 = -t481 * qJ(2) + qJDD(2) - t465;
t444 = t512 * qJDD(1) + t486;
t475 = sin(qJ(3));
t478 = cos(qJ(3));
t435 = -t478 * g(3) + t475 * t444;
t459 = (mrSges(4,1) * t475 + mrSges(4,2) * t478) * qJD(1);
t497 = qJD(1) * qJD(3);
t492 = t478 * t497;
t461 = -t475 * qJDD(1) - t492;
t499 = qJD(1) * t478;
t464 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t499;
t443 = t512 * t481 - t513;
t493 = t475 * t497;
t462 = t478 * qJDD(1) - t493;
t413 = (-t462 + t493) * pkin(7) + (-t461 + t492) * pkin(3) + t443;
t460 = (pkin(3) * t475 - pkin(7) * t478) * qJD(1);
t480 = qJD(3) ^ 2;
t498 = t475 * qJD(1);
t416 = -t480 * pkin(3) + qJDD(3) * pkin(7) - t460 * t498 + t435;
t474 = sin(qJ(4));
t477 = cos(qJ(4));
t409 = t477 * t413 - t474 * t416;
t457 = t477 * qJD(3) - t474 * t499;
t430 = t457 * qJD(4) + t474 * qJDD(3) + t477 * t462;
t458 = t474 * qJD(3) + t477 * t499;
t432 = -t457 * mrSges(6,1) + t458 * mrSges(6,2);
t433 = -t457 * mrSges(5,1) + t458 * mrSges(5,2);
t467 = qJD(4) + t498;
t437 = -t467 * mrSges(5,2) + t457 * mrSges(5,3);
t456 = qJDD(4) - t461;
t405 = -0.2e1 * qJD(5) * t458 + (t457 * t467 - t430) * qJ(5) + (t457 * t458 + t456) * pkin(4) + t409;
t436 = -t467 * mrSges(6,2) + t457 * mrSges(6,3);
t495 = m(6) * t405 + t456 * mrSges(6,1) + t467 * t436;
t398 = m(5) * t409 + t456 * mrSges(5,1) + t467 * t437 + (-t432 - t433) * t458 + (-mrSges(5,3) - mrSges(6,3)) * t430 + t495;
t410 = t474 * t413 + t477 * t416;
t429 = -t458 * qJD(4) + t477 * qJDD(3) - t474 * t462;
t438 = t467 * pkin(4) - t458 * qJ(5);
t455 = t457 ^ 2;
t407 = -t455 * pkin(4) + t429 * qJ(5) + 0.2e1 * qJD(5) * t457 - t467 * t438 + t410;
t494 = m(6) * t407 + t429 * mrSges(6,3) + t457 * t432;
t439 = t467 * mrSges(6,1) - t458 * mrSges(6,3);
t500 = -t467 * mrSges(5,1) + t458 * mrSges(5,3) - t439;
t401 = m(5) * t410 + t429 * mrSges(5,3) + t457 * t433 + t510 * t456 + t500 * t467 + t494;
t489 = -t474 * t398 + t477 * t401;
t395 = m(4) * t435 - qJDD(3) * mrSges(4,2) + t461 * mrSges(4,3) - qJD(3) * t464 - t459 * t498 + t489;
t434 = t475 * g(3) + t478 * t444;
t463 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t498;
t415 = -qJDD(3) * pkin(3) - t480 * pkin(7) + t460 * t499 - t434;
t408 = -t429 * pkin(4) - t455 * qJ(5) + t458 * t438 + qJDD(5) + t415;
t488 = m(6) * t408 - t429 * mrSges(6,1) - t457 * t436;
t482 = -m(5) * t415 + t429 * mrSges(5,1) + t510 * t430 + t457 * t437 + t500 * t458 - t488;
t402 = m(4) * t434 + qJDD(3) * mrSges(4,1) - t462 * mrSges(4,3) + qJD(3) * t463 - t459 * t499 + t482;
t388 = t475 * t395 + t478 * t402;
t446 = -qJDD(1) * pkin(1) + t486;
t485 = -m(3) * t446 + (t481 * mrSges(3,3)) - t388;
t386 = m(2) * t465 - (t481 * mrSges(2,2)) + t511 * qJDD(1) + t485;
t445 = t481 * pkin(1) + t513;
t396 = t477 * t398 + t474 * t401;
t484 = -m(4) * t443 + t461 * mrSges(4,1) - t462 * mrSges(4,2) - t463 * t498 - t464 * t499 - t396;
t483 = -m(3) * t445 + (t481 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t484;
t392 = m(2) * t466 - (t481 * mrSges(2,1)) - qJDD(1) * mrSges(2,2) + t483;
t504 = t479 * t386 + t476 * t392;
t503 = t515 * t457 + t506 * t458 + t514 * t467;
t502 = -t516 * t457 - t508 * t458 - t515 * t467;
t501 = t508 * t457 + t517 * t458 + t506 * t467;
t491 = -t476 * t386 + t479 * t392;
t490 = t478 * t395 - t475 * t402;
t450 = (Ifges(4,5) * qJD(3)) + (Ifges(4,1) * t478 - Ifges(4,4) * t475) * qJD(1);
t449 = (Ifges(4,6) * qJD(3)) + (Ifges(4,4) * t478 - Ifges(4,2) * t475) * qJD(1);
t448 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t478 - Ifges(4,6) * t475) * qJD(1);
t403 = -t430 * mrSges(6,3) - t458 * t432 + t495;
t393 = mrSges(5,2) * t415 + mrSges(6,2) * t408 - mrSges(5,3) * t409 - mrSges(6,3) * t405 - qJ(5) * t403 + t508 * t429 + t517 * t430 + t506 * t456 + t503 * t457 + t502 * t467;
t389 = -mrSges(5,1) * t415 + mrSges(5,3) * t410 - mrSges(6,1) * t408 + mrSges(6,3) * t407 - pkin(4) * t488 + qJ(5) * t494 + (-qJ(5) * t439 + t501) * t467 + (-pkin(4) * t439 - t503) * t458 + (-qJ(5) * mrSges(6,2) + t515) * t456 + (-pkin(4) * mrSges(6,2) + t508) * t430 + t516 * t429;
t387 = -m(3) * g(3) + t490;
t384 = -t448 * t499 - mrSges(4,1) * t443 - mrSges(5,1) * t409 - mrSges(6,1) * t405 + mrSges(5,2) * t410 + mrSges(6,2) * t407 + mrSges(4,3) * t435 + Ifges(4,4) * t462 + Ifges(4,2) * t461 + Ifges(4,6) * qJDD(3) - pkin(3) * t396 - pkin(4) * t403 + qJD(3) * t450 + t502 * t458 + t501 * t457 - t514 * t456 - t506 * t430 - t515 * t429;
t383 = mrSges(4,2) * t443 - mrSges(4,3) * t434 + Ifges(4,1) * t462 + Ifges(4,4) * t461 + Ifges(4,5) * qJDD(3) - pkin(7) * t396 - qJD(3) * t449 - t474 * t389 + t477 * t393 - t448 * t498;
t382 = -qJ(2) * t387 - mrSges(2,3) * t465 + pkin(2) * t388 + mrSges(3,1) * t446 + t477 * t389 + pkin(3) * t482 + pkin(7) * t489 + t474 * t393 + Ifges(4,5) * t462 + Ifges(4,6) * t461 + Ifges(4,3) * qJDD(3) + mrSges(4,1) * t434 - mrSges(4,2) * t435 + (t507 * t481) + t509 * qJDD(1) + (t478 * t449 + t475 * t450) * qJD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3);
t381 = -mrSges(3,1) * t445 + mrSges(2,3) * t466 - pkin(1) * t387 - pkin(2) * t484 - pkin(6) * t490 + t511 * g(3) - t507 * qJDD(1) - t475 * t383 - t478 * t384 + t509 * t481;
t1 = [-m(1) * g(1) + t491; -m(1) * g(2) + t504; (-m(1) - m(2) - m(3)) * g(3) + t490; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t504 - t476 * t381 + t479 * t382; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t491 + t479 * t381 + t476 * t382; qJ(2) * t483 + pkin(1) * t485 + t478 * t383 - t475 * t384 - pkin(6) * t388 + mrSges(2,1) * t465 - mrSges(2,2) * t466 + mrSges(3,2) * t446 - mrSges(3,3) * t445 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-pkin(1) * mrSges(3,2) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
