% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRPRP1
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
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:29
% EndTime: 2019-12-05 15:28:31
% DurationCPUTime: 2.06s
% Computational Cost: add. (16812->214), mult. (37344->263), div. (0->0), fcn. (24046->8), ass. (0->99)
t525 = Ifges(5,1) + Ifges(6,1);
t519 = Ifges(5,4) - Ifges(6,5);
t518 = Ifges(5,5) + Ifges(6,4);
t524 = Ifges(5,2) + Ifges(6,3);
t517 = Ifges(5,6) - Ifges(6,6);
t523 = -Ifges(5,3) - Ifges(6,2);
t486 = qJD(2) ^ 2;
t522 = cos(qJ(4));
t480 = cos(pkin(8));
t521 = pkin(3) * t480;
t520 = -mrSges(5,3) - mrSges(6,2);
t478 = sin(pkin(8));
t516 = mrSges(4,2) * t478;
t474 = t480 ^ 2;
t515 = t474 * t486;
t479 = sin(pkin(7));
t481 = cos(pkin(7));
t461 = t479 * g(1) - t481 * g(2);
t462 = -t481 * g(1) - t479 * g(2);
t483 = sin(qJ(2));
t484 = cos(qJ(2));
t447 = t483 * t461 + t484 * t462;
t445 = -t486 * pkin(2) + qJDD(2) * qJ(3) + t447;
t477 = -g(3) + qJDD(1);
t505 = qJD(2) * qJD(3);
t509 = t480 * t477 - 0.2e1 * t478 * t505;
t421 = (-pkin(6) * qJDD(2) + t486 * t521 - t445) * t478 + t509;
t425 = t478 * t477 + (t445 + 0.2e1 * t505) * t480;
t503 = qJDD(2) * t480;
t422 = -pkin(3) * t515 + pkin(6) * t503 + t425;
t482 = sin(qJ(4));
t418 = t482 * t421 + t522 * t422;
t500 = t480 * t522;
t504 = qJDD(2) * t478;
t489 = t522 * t478 + t480 * t482;
t455 = t489 * qJD(2);
t507 = t455 * qJD(4);
t443 = -qJDD(2) * t500 + t482 * t504 + t507;
t451 = qJD(4) * mrSges(5,1) - t455 * mrSges(5,3);
t506 = t478 * qJD(2);
t454 = -qJD(2) * t500 + t482 * t506;
t436 = t454 * pkin(4) - t455 * qJ(5);
t485 = qJD(4) ^ 2;
t413 = -t485 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t454 * t436 + t418;
t452 = -qJD(4) * mrSges(6,1) + t455 * mrSges(6,2);
t502 = m(6) * t413 + qJDD(4) * mrSges(6,3) + qJD(4) * t452;
t437 = t454 * mrSges(6,1) - t455 * mrSges(6,3);
t510 = -t454 * mrSges(5,1) - t455 * mrSges(5,2) - t437;
t409 = m(5) * t418 - qJDD(4) * mrSges(5,2) - qJD(4) * t451 + t520 * t443 + t510 * t454 + t502;
t417 = t522 * t421 - t482 * t422;
t508 = t454 * qJD(4);
t444 = t489 * qJDD(2) - t508;
t450 = -qJD(4) * mrSges(5,2) - t454 * mrSges(5,3);
t414 = -qJDD(4) * pkin(4) - t485 * qJ(5) + t455 * t436 + qJDD(5) - t417;
t453 = -t454 * mrSges(6,2) + qJD(4) * mrSges(6,3);
t495 = -m(6) * t414 + qJDD(4) * mrSges(6,1) + qJD(4) * t453;
t410 = m(5) * t417 + qJDD(4) * mrSges(5,1) + qJD(4) * t450 + t520 * t444 + t510 * t455 + t495;
t403 = t482 * t409 + t522 * t410;
t424 = -t478 * t445 + t509;
t490 = mrSges(4,3) * qJDD(2) + t486 * (-mrSges(4,1) * t480 + t516);
t401 = m(4) * t424 - t490 * t478 + t403;
t496 = t522 * t409 - t482 * t410;
t402 = m(4) * t425 + t490 * t480 + t496;
t497 = -t478 * t401 + t480 * t402;
t394 = m(3) * t447 - t486 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t497;
t446 = t484 * t461 - t483 * t462;
t491 = qJDD(3) - t446;
t442 = -qJDD(2) * pkin(2) - t486 * qJ(3) + t491;
t473 = t478 ^ 2;
t423 = (-pkin(2) - t521) * qJDD(2) + (-qJ(3) + (-t473 - t474) * pkin(6)) * t486 + t491;
t416 = -0.2e1 * qJD(5) * t455 + (-t444 + t508) * qJ(5) + (t443 + t507) * pkin(4) + t423;
t411 = m(6) * t416 + t443 * mrSges(6,1) - t444 * mrSges(6,3) - t455 * t452 + t454 * t453;
t488 = m(5) * t423 + t443 * mrSges(5,1) + t444 * mrSges(5,2) + t454 * t450 + t455 * t451 + t411;
t487 = -m(4) * t442 + mrSges(4,1) * t503 - t488 + (t473 * t486 + t515) * mrSges(4,3);
t405 = -t486 * mrSges(3,2) + m(3) * t446 + t487 + (mrSges(3,1) - t516) * qJDD(2);
t391 = t483 * t394 + t484 * t405;
t389 = m(2) * t461 + t391;
t498 = t484 * t394 - t483 * t405;
t390 = m(2) * t462 + t498;
t514 = t481 * t389 + t479 * t390;
t395 = t480 * t401 + t478 * t402;
t513 = -t517 * qJD(4) + t524 * t454 - t519 * t455;
t512 = t523 * qJD(4) + t517 * t454 - t518 * t455;
t511 = t518 * qJD(4) - t519 * t454 + t525 * t455;
t501 = m(3) * t477 + t395;
t499 = -t479 * t389 + t481 * t390;
t494 = Ifges(4,1) * t478 + Ifges(4,4) * t480;
t493 = Ifges(4,4) * t478 + Ifges(4,2) * t480;
t492 = Ifges(4,5) * t478 + Ifges(4,6) * t480;
t460 = t492 * qJD(2);
t397 = mrSges(5,2) * t423 + mrSges(6,2) * t414 - mrSges(5,3) * t417 - mrSges(6,3) * t416 - qJ(5) * t411 + t513 * qJD(4) + t518 * qJDD(4) - t519 * t443 + t525 * t444 + t512 * t454;
t396 = -mrSges(5,1) * t423 - mrSges(6,1) * t416 + mrSges(6,2) * t413 + mrSges(5,3) * t418 - pkin(4) * t411 + t511 * qJD(4) + t517 * qJDD(4) - t524 * t443 + t519 * t444 + t512 * t455;
t385 = t480 * qJD(2) * t460 + mrSges(4,2) * t442 - mrSges(4,3) * t424 - pkin(6) * t403 + t494 * qJDD(2) - t482 * t396 + t522 * t397;
t384 = -mrSges(4,1) * t442 + mrSges(4,3) * t425 - pkin(3) * t488 + pkin(6) * t496 + t493 * qJDD(2) + t522 * t396 + t482 * t397 - t460 * t506;
t383 = -qJ(5) * t502 - pkin(4) * t495 - mrSges(3,1) * t477 + mrSges(3,3) * t447 - mrSges(4,1) * t424 + mrSges(4,2) * t425 - mrSges(5,1) * t417 + mrSges(5,2) * t418 - mrSges(6,3) * t413 + mrSges(6,1) * t414 - pkin(3) * t403 - pkin(2) * t395 + (pkin(4) * t437 + t513) * t455 + (qJ(5) * t437 - t511) * t454 + (pkin(4) * mrSges(6,2) - t518) * t444 + (qJ(5) * mrSges(6,2) + t517) * t443 + t523 * qJDD(4) + (Ifges(3,6) - t492) * qJDD(2) + (-t478 * t493 + t480 * t494 + Ifges(3,5)) * t486;
t382 = mrSges(3,2) * t477 - mrSges(3,3) * t446 + Ifges(3,5) * qJDD(2) - t486 * Ifges(3,6) - qJ(3) * t395 - t478 * t384 + t480 * t385;
t381 = mrSges(2,2) * t477 - mrSges(2,3) * t461 - pkin(5) * t391 + t484 * t382 - t483 * t383;
t380 = -mrSges(2,1) * t477 + mrSges(2,3) * t462 - pkin(1) * t501 + pkin(5) * t498 + t483 * t382 + t484 * t383;
t1 = [-m(1) * g(1) + t499; -m(1) * g(2) + t514; -m(1) * g(3) + m(2) * t477 + t501; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t514 - t479 * t380 + t481 * t381; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t499 + t481 * t380 + t479 * t381; pkin(1) * t391 + mrSges(2,1) * t461 - mrSges(2,2) * t462 + t478 * t385 + t480 * t384 + pkin(2) * (-mrSges(4,2) * t504 + t487) + qJ(3) * t497 + mrSges(3,1) * t446 - mrSges(3,2) * t447 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * qJDD(2);];
tauB = t1;
