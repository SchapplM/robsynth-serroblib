% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 02:49
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 02:44:30
% EndTime: 2019-05-05 02:44:32
% DurationCPUTime: 1.87s
% Computational Cost: add. (12310->273), mult. (26772->333), div. (0->0), fcn. (18268->12), ass. (0->117)
t534 = -2 * qJD(4);
t533 = Ifges(5,1) + Ifges(6,2);
t532 = -Ifges(6,1) - Ifges(5,3);
t528 = Ifges(5,4) + Ifges(6,6);
t527 = Ifges(5,5) - Ifges(6,4);
t531 = Ifges(5,2) + Ifges(6,3);
t526 = Ifges(5,6) - Ifges(6,5);
t488 = sin(pkin(10));
t490 = cos(pkin(10));
t477 = g(1) * t488 - g(2) * t490;
t486 = -g(3) + qJDD(1);
t489 = sin(pkin(6));
t491 = cos(pkin(6));
t530 = t477 * t491 + t486 * t489;
t478 = -g(1) * t490 - g(2) * t488;
t494 = sin(qJ(2));
t497 = cos(qJ(2));
t430 = t497 * t478 + t494 * t530;
t499 = qJD(2) ^ 2;
t422 = -pkin(2) * t499 + qJDD(2) * pkin(8) + t430;
t455 = -t477 * t489 + t486 * t491;
t493 = sin(qJ(3));
t496 = cos(qJ(3));
t404 = -t493 * t422 + t496 * t455;
t513 = qJD(2) * qJD(3);
t512 = t496 * t513;
t475 = qJDD(2) * t493 + t512;
t400 = (-t475 + t512) * qJ(4) + (t493 * t496 * t499 + qJDD(3)) * pkin(3) + t404;
t405 = t496 * t422 + t493 * t455;
t476 = qJDD(2) * t496 - t493 * t513;
t517 = qJD(2) * t493;
t479 = qJD(3) * pkin(3) - qJ(4) * t517;
t485 = t496 ^ 2;
t401 = -pkin(3) * t485 * t499 + qJ(4) * t476 - qJD(3) * t479 + t405;
t487 = sin(pkin(11));
t525 = cos(pkin(11));
t463 = (t487 * t496 + t493 * t525) * qJD(2);
t393 = t400 * t525 - t487 * t401 + t463 * t534;
t429 = -t494 * t478 + t497 * t530;
t529 = -2 * qJD(5);
t516 = qJD(2) * t496;
t462 = t487 * t517 - t516 * t525;
t435 = mrSges(5,1) * t462 + mrSges(5,2) * t463;
t444 = t475 * t525 + t487 * t476;
t450 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t462;
t452 = mrSges(6,1) * t462 - qJD(3) * mrSges(6,3);
t434 = pkin(4) * t462 - qJ(5) * t463;
t498 = qJD(3) ^ 2;
t391 = -qJDD(3) * pkin(4) - t498 * qJ(5) + t463 * t434 + qJDD(5) - t393;
t515 = qJD(3) * t462;
t387 = (t462 * t463 - qJDD(3)) * pkin(9) + (t444 + t515) * pkin(5) + t391;
t443 = t487 * t475 - t476 * t525;
t454 = pkin(5) * t463 - qJD(3) * pkin(9);
t461 = t462 ^ 2;
t503 = -qJDD(2) * pkin(2) - t429;
t403 = -t476 * pkin(3) + qJDD(4) + t479 * t517 + (-qJ(4) * t485 - pkin(8)) * t499 + t503;
t500 = (-t444 + t515) * qJ(5) + t403 + (pkin(4) * qJD(3) + t529) * t463;
t392 = -t461 * pkin(5) - t463 * t454 + (pkin(4) + pkin(9)) * t443 + t500;
t492 = sin(qJ(6));
t495 = cos(qJ(6));
t385 = t387 * t495 - t392 * t492;
t445 = -qJD(3) * t492 + t462 * t495;
t414 = qJD(6) * t445 + qJDD(3) * t495 + t443 * t492;
t446 = qJD(3) * t495 + t462 * t492;
t415 = -mrSges(7,1) * t445 + mrSges(7,2) * t446;
t460 = qJD(6) + t463;
t419 = -mrSges(7,2) * t460 + mrSges(7,3) * t445;
t442 = qJDD(6) + t444;
t382 = m(7) * t385 + mrSges(7,1) * t442 - mrSges(7,3) * t414 - t415 * t446 + t419 * t460;
t386 = t387 * t492 + t392 * t495;
t413 = -qJD(6) * t446 - qJDD(3) * t492 + t443 * t495;
t420 = mrSges(7,1) * t460 - mrSges(7,3) * t446;
t383 = m(7) * t386 - mrSges(7,2) * t442 + mrSges(7,3) * t413 + t415 * t445 - t420 * t460;
t374 = t495 * t382 + t492 * t383;
t436 = -mrSges(6,2) * t462 - mrSges(6,3) * t463;
t505 = -m(6) * t391 - t444 * mrSges(6,1) - t463 * t436 - t374;
t370 = m(5) * t393 - t444 * mrSges(5,3) - t463 * t435 + (mrSges(5,1) - mrSges(6,2)) * qJDD(3) + (t450 - t452) * qJD(3) + t505;
t458 = t462 * t534;
t521 = t487 * t400 + t525 * t401;
t394 = t458 + t521;
t451 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t463;
t507 = t498 * pkin(4) - qJDD(3) * qJ(5) - t521;
t390 = qJD(3) * t529 + ((2 * qJD(4)) + t434) * t462 + t507;
t453 = mrSges(6,1) * t463 + qJD(3) * mrSges(6,2);
t389 = -t443 * pkin(5) - t461 * pkin(9) - t462 * t434 + t458 + ((2 * qJD(5)) + t454) * qJD(3) - t507;
t506 = -m(7) * t389 + t413 * mrSges(7,1) - t414 * mrSges(7,2) + t445 * t419 - t446 * t420;
t502 = -m(6) * t390 + qJDD(3) * mrSges(6,3) + qJD(3) * t453 - t506;
t379 = m(5) * t394 - qJDD(3) * mrSges(5,2) - qJD(3) * t451 + (-t435 - t436) * t462 + (-mrSges(5,3) - mrSges(6,1)) * t443 + t502;
t368 = t525 * t370 + t487 * t379;
t522 = -t492 * t382 + t495 * t383;
t520 = qJD(3) * t532 + t462 * t526 - t463 * t527;
t519 = -qJD(3) * t526 + t462 * t531 - t463 * t528;
t518 = qJD(3) * t527 - t462 * t528 + t463 * t533;
t474 = (-mrSges(4,1) * t496 + mrSges(4,2) * t493) * qJD(2);
t481 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t516;
t366 = m(4) * t404 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t475 + qJD(3) * t481 - t474 * t517 + t368;
t480 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t517;
t510 = -t370 * t487 + t525 * t379;
t367 = m(4) * t405 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t476 - qJD(3) * t480 + t474 * t516 + t510;
t511 = -t366 * t493 + t496 * t367;
t396 = t443 * pkin(4) + t500;
t373 = m(6) * t396 - t443 * mrSges(6,2) - t444 * mrSges(6,3) - t462 * t452 - t463 * t453 + t522;
t407 = Ifges(7,4) * t446 + Ifges(7,2) * t445 + Ifges(7,6) * t460;
t408 = Ifges(7,1) * t446 + Ifges(7,4) * t445 + Ifges(7,5) * t460;
t504 = mrSges(7,1) * t385 - mrSges(7,2) * t386 + Ifges(7,5) * t414 + Ifges(7,6) * t413 + Ifges(7,3) * t442 + t446 * t407 - t445 * t408;
t371 = m(5) * t403 + t443 * mrSges(5,1) + mrSges(5,2) * t444 + t462 * t450 + t451 * t463 + t373;
t421 = -t499 * pkin(8) + t503;
t501 = -m(4) * t421 + t476 * mrSges(4,1) - mrSges(4,2) * t475 - t480 * t517 + t481 * t516 - t371;
t467 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t493 + Ifges(4,4) * t496) * qJD(2);
t466 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t493 + Ifges(4,2) * t496) * qJD(2);
t406 = Ifges(7,5) * t446 + Ifges(7,6) * t445 + Ifges(7,3) * t460;
t378 = mrSges(7,2) * t389 - mrSges(7,3) * t385 + Ifges(7,1) * t414 + Ifges(7,4) * t413 + Ifges(7,5) * t442 + t406 * t445 - t407 * t460;
t377 = -mrSges(7,1) * t389 + mrSges(7,3) * t386 + Ifges(7,4) * t414 + Ifges(7,2) * t413 + Ifges(7,6) * t442 - t406 * t446 + t408 * t460;
t372 = qJDD(3) * mrSges(6,2) + qJD(3) * t452 - t505;
t364 = mrSges(6,1) * t391 + mrSges(5,2) * t403 - mrSges(5,3) * t393 - mrSges(6,3) * t396 + pkin(5) * t374 - qJ(5) * t373 + t519 * qJD(3) + t527 * qJDD(3) - t528 * t443 + t444 * t533 + t520 * t462 + t504;
t363 = -mrSges(5,1) * t403 - mrSges(6,1) * t390 + mrSges(6,2) * t396 + mrSges(5,3) * t394 - pkin(4) * t373 - pkin(5) * t506 - pkin(9) * t522 + t518 * qJD(3) + t526 * qJDD(3) - t495 * t377 - t492 * t378 - t443 * t531 + t528 * t444 + t520 * t463;
t1 = [m(2) * t486 + t491 * (m(3) * t455 + t366 * t496 + t367 * t493) + (t494 * (m(3) * t430 - mrSges(3,1) * t499 - qJDD(2) * mrSges(3,2) + t511) + t497 * (m(3) * t429 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t499 + t501)) * t489; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t429 - mrSges(3,2) * t430 + t493 * (mrSges(4,2) * t421 - mrSges(4,3) * t404 + Ifges(4,1) * t475 + Ifges(4,4) * t476 + Ifges(4,5) * qJDD(3) - qJ(4) * t368 - qJD(3) * t466 - t487 * t363 + t364 * t525) + t496 * (-mrSges(4,1) * t421 + mrSges(4,3) * t405 + Ifges(4,4) * t475 + Ifges(4,2) * t476 + Ifges(4,6) * qJDD(3) - pkin(3) * t371 + qJ(4) * t510 + qJD(3) * t467 + t363 * t525 + t487 * t364) + pkin(2) * t501 + pkin(8) * t511; -t492 * t377 + t495 * t378 + Ifges(4,6) * t476 + qJ(5) * t502 + Ifges(4,5) * t475 + mrSges(4,1) * t404 - mrSges(4,2) * t405 + mrSges(5,1) * t393 - mrSges(5,2) * t394 - mrSges(6,3) * t390 + mrSges(6,2) * t391 - pkin(9) * t374 - pkin(4) * t372 + pkin(3) * t368 - t519 * t463 + (-qJ(5) * t436 + t518) * t462 + t527 * t444 + (-mrSges(6,1) * qJ(5) - t526) * t443 + (t466 * t493 - t467 * t496) * qJD(2) + (Ifges(4,3) - t532) * qJDD(3); t371; t372; t504;];
tauJ  = t1;
