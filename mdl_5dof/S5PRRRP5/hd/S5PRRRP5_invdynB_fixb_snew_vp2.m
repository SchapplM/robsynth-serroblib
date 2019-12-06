% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRP5
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRP5_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:48:00
% EndTime: 2019-12-05 16:48:04
% DurationCPUTime: 2.18s
% Computational Cost: add. (20602->238), mult. (41475->292), div. (0->0), fcn. (25861->8), ass. (0->95)
t520 = Ifges(5,1) + Ifges(6,1);
t517 = Ifges(5,4) + Ifges(6,4);
t516 = Ifges(5,5) + Ifges(6,5);
t519 = Ifges(5,2) + Ifges(6,2);
t515 = -Ifges(5,6) - Ifges(6,6);
t518 = -Ifges(5,3) - Ifges(6,3);
t514 = cos(pkin(8));
t488 = sin(pkin(8));
t477 = -t514 * g(1) - t488 * g(2);
t487 = -g(3) + qJDD(1);
t491 = sin(qJ(2));
t494 = cos(qJ(2));
t459 = t494 * t477 + t491 * t487;
t495 = qJD(2) ^ 2;
t452 = -t495 * pkin(2) + qJDD(2) * pkin(6) + t459;
t476 = t488 * g(1) - t514 * g(2);
t490 = sin(qJ(3));
t493 = cos(qJ(3));
t436 = -t490 * t452 - t493 * t476;
t507 = qJD(2) * qJD(3);
t504 = t493 * t507;
t474 = t490 * qJDD(2) + t504;
t426 = (-t474 + t504) * pkin(7) + (t490 * t493 * t495 + qJDD(3)) * pkin(3) + t436;
t437 = t493 * t452 - t490 * t476;
t475 = t493 * qJDD(2) - t490 * t507;
t509 = qJD(2) * t490;
t480 = qJD(3) * pkin(3) - pkin(7) * t509;
t486 = t493 ^ 2;
t427 = -t486 * t495 * pkin(3) + t475 * pkin(7) - qJD(3) * t480 + t437;
t489 = sin(qJ(4));
t492 = cos(qJ(4));
t421 = t492 * t426 - t489 * t427;
t464 = (-t489 * t490 + t492 * t493) * qJD(2);
t435 = t464 * qJD(4) + t492 * t474 + t489 * t475;
t465 = (t489 * t493 + t490 * t492) * qJD(2);
t447 = -t464 * mrSges(6,1) + t465 * mrSges(6,2);
t448 = -t464 * mrSges(5,1) + t465 * mrSges(5,2);
t485 = qJD(3) + qJD(4);
t454 = -t485 * mrSges(5,2) + t464 * mrSges(5,3);
t484 = qJDD(3) + qJDD(4);
t416 = -0.2e1 * qJD(5) * t465 + (t464 * t485 - t435) * qJ(5) + (t464 * t465 + t484) * pkin(4) + t421;
t453 = -t485 * mrSges(6,2) + t464 * mrSges(6,3);
t506 = m(6) * t416 + t484 * mrSges(6,1) + t485 * t453;
t408 = m(5) * t421 + t484 * mrSges(5,1) + t485 * t454 + (-t447 - t448) * t465 + (-mrSges(5,3) - mrSges(6,3)) * t435 + t506;
t422 = t489 * t426 + t492 * t427;
t434 = -t465 * qJD(4) - t489 * t474 + t492 * t475;
t456 = t485 * mrSges(6,1) - t465 * mrSges(6,3);
t457 = t485 * mrSges(5,1) - t465 * mrSges(5,3);
t455 = t485 * pkin(4) - t465 * qJ(5);
t460 = t464 ^ 2;
t418 = -t460 * pkin(4) + t434 * qJ(5) + 0.2e1 * qJD(5) * t464 - t485 * t455 + t422;
t505 = m(6) * t418 + t434 * mrSges(6,3) + t464 * t447;
t413 = m(5) * t422 + t434 * mrSges(5,3) + t464 * t448 + (-t456 - t457) * t485 + (-mrSges(5,2) - mrSges(6,2)) * t484 + t505;
t406 = t492 * t408 + t489 * t413;
t473 = (-mrSges(4,1) * t493 + mrSges(4,2) * t490) * qJD(2);
t508 = qJD(2) * t493;
t479 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t508;
t404 = m(4) * t436 + qJDD(3) * mrSges(4,1) - t474 * mrSges(4,3) + qJD(3) * t479 - t473 * t509 + t406;
t478 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t509;
t500 = -t489 * t408 + t492 * t413;
t405 = m(4) * t437 - qJDD(3) * mrSges(4,2) + t475 * mrSges(4,3) - qJD(3) * t478 + t473 * t508 + t500;
t501 = -t490 * t404 + t493 * t405;
t397 = m(3) * t459 - t495 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t501;
t458 = -t491 * t477 + t494 * t487;
t498 = -qJDD(2) * pkin(2) - t458;
t451 = -t495 * pkin(6) + t498;
t428 = -t475 * pkin(3) + t480 * t509 + (-pkin(7) * t486 - pkin(6)) * t495 + t498;
t420 = -t434 * pkin(4) - t460 * qJ(5) + t465 * t455 + qJDD(5) + t428;
t499 = m(6) * t420 - t434 * mrSges(6,1) + t435 * mrSges(6,2) - t464 * t453 + t465 * t456;
t497 = m(5) * t428 - t434 * mrSges(5,1) + t435 * mrSges(5,2) - t464 * t454 + t465 * t457 + t499;
t496 = -m(4) * t451 + t475 * mrSges(4,1) - t474 * mrSges(4,2) - t478 * t509 + t479 * t508 - t497;
t412 = m(3) * t458 + qJDD(2) * mrSges(3,1) - t495 * mrSges(3,2) + t496;
t502 = t494 * t397 - t491 * t412;
t393 = m(2) * t477 + t502;
t400 = t493 * t404 + t490 * t405;
t399 = (m(2) + m(3)) * t476 - t400;
t513 = t488 * t393 + t514 * t399;
t394 = t491 * t397 + t494 * t412;
t512 = t515 * t464 - t516 * t465 + t518 * t485;
t511 = -t519 * t464 - t517 * t465 + t515 * t485;
t510 = t517 * t464 + t520 * t465 + t516 * t485;
t503 = t514 * t393 - t488 * t399;
t463 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t490 + Ifges(4,4) * t493) * qJD(2);
t462 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t490 + Ifges(4,2) * t493) * qJD(2);
t461 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t490 + Ifges(4,6) * t493) * qJD(2);
t414 = -t435 * mrSges(6,3) - t465 * t447 + t506;
t403 = mrSges(5,2) * t428 + mrSges(6,2) * t420 - mrSges(5,3) * t421 - mrSges(6,3) * t416 - qJ(5) * t414 + t517 * t434 + t520 * t435 - t512 * t464 + t516 * t484 + t511 * t485;
t401 = -mrSges(5,1) * t428 + mrSges(5,3) * t422 - mrSges(6,1) * t420 + mrSges(6,3) * t418 - pkin(4) * t499 + qJ(5) * t505 + (-qJ(5) * t456 + t510) * t485 + (-qJ(5) * mrSges(6,2) - t515) * t484 + t512 * t465 + t517 * t435 + t519 * t434;
t390 = mrSges(4,2) * t451 - mrSges(4,3) * t436 + Ifges(4,1) * t474 + Ifges(4,4) * t475 + Ifges(4,5) * qJDD(3) - pkin(7) * t406 - qJD(3) * t462 - t489 * t401 + t492 * t403 + t461 * t508;
t389 = -mrSges(4,1) * t451 + mrSges(4,3) * t437 + Ifges(4,4) * t474 + Ifges(4,2) * t475 + Ifges(4,6) * qJDD(3) - pkin(3) * t497 + pkin(7) * t500 + qJD(3) * t463 + t492 * t401 + t489 * t403 - t461 * t509;
t388 = -Ifges(4,3) * qJDD(3) + Ifges(3,6) * qJDD(2) + t495 * Ifges(3,5) + mrSges(3,1) * t476 - Ifges(4,5) * t474 - Ifges(4,6) * t475 + mrSges(3,3) * t459 - mrSges(4,1) * t436 + mrSges(4,2) * t437 + mrSges(6,2) * t418 - mrSges(5,1) * t421 + mrSges(5,2) * t422 - pkin(4) * t414 - mrSges(6,1) * t416 - pkin(3) * t406 - pkin(2) * t400 + t518 * t484 + t511 * t465 + t510 * t464 - t516 * t435 + t515 * t434 + (-t490 * t462 + t493 * t463) * qJD(2);
t387 = -mrSges(3,2) * t476 - mrSges(3,3) * t458 + Ifges(3,5) * qJDD(2) - t495 * Ifges(3,6) - pkin(6) * t400 - t490 * t389 + t493 * t390;
t386 = -mrSges(2,1) * t487 - mrSges(3,1) * t458 + mrSges(3,2) * t459 + mrSges(2,3) * t477 - Ifges(3,3) * qJDD(2) - pkin(1) * t394 - pkin(2) * t496 - pkin(6) * t501 - t493 * t389 - t490 * t390;
t385 = mrSges(2,2) * t487 - mrSges(2,3) * t476 - pkin(5) * t394 + t494 * t387 - t491 * t388;
t1 = [-m(1) * g(1) + t503; -m(1) * g(2) + t513; -m(1) * g(3) + m(2) * t487 + t394; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t513 + t514 * t385 - t488 * t386; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t503 + t488 * t385 + t514 * t386; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + mrSges(2,1) * t476 - mrSges(2,2) * t477 + t491 * t387 + t494 * t388 + pkin(1) * (m(3) * t476 - t400) + pkin(5) * t502;];
tauB = t1;
