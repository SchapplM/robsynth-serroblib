% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPRRP3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:47:10
% EndTime: 2020-01-03 11:47:14
% DurationCPUTime: 2.69s
% Computational Cost: add. (23998->248), mult. (47662->301), div. (0->0), fcn. (28134->8), ass. (0->98)
t541 = Ifges(5,1) + Ifges(6,1);
t538 = Ifges(5,4) + Ifges(6,4);
t537 = Ifges(5,5) + Ifges(6,5);
t540 = Ifges(5,2) + Ifges(6,2);
t536 = -Ifges(5,6) - Ifges(6,6);
t539 = -Ifges(5,3) - Ifges(6,3);
t512 = sin(qJ(1));
t515 = cos(qJ(1));
t497 = -t512 * g(2) + t515 * g(3);
t516 = qJD(1) ^ 2;
t498 = -t515 * g(2) - t512 * g(3);
t489 = qJDD(1) * pkin(1) + t498;
t491 = -t516 * pkin(1) + t497;
t508 = sin(pkin(8));
t509 = cos(pkin(8));
t470 = t508 * t489 + t509 * t491;
t466 = -t516 * pkin(2) + qJDD(1) * pkin(6) + t470;
t507 = -g(1) + qJDD(2);
t511 = sin(qJ(3));
t514 = cos(qJ(3));
t451 = -t511 * t466 + t514 * t507;
t529 = qJD(1) * qJD(3);
t525 = t514 * t529;
t492 = t511 * qJDD(1) + t525;
t443 = (-t492 + t525) * pkin(7) + (t511 * t514 * t516 + qJDD(3)) * pkin(3) + t451;
t452 = t514 * t466 + t511 * t507;
t493 = t514 * qJDD(1) - t511 * t529;
t531 = qJD(1) * t511;
t496 = qJD(3) * pkin(3) - pkin(7) * t531;
t506 = t514 ^ 2;
t444 = -t506 * t516 * pkin(3) + t493 * pkin(7) - qJD(3) * t496 + t452;
t510 = sin(qJ(4));
t513 = cos(qJ(4));
t438 = t513 * t443 - t510 * t444;
t484 = (-t510 * t511 + t513 * t514) * qJD(1);
t454 = t484 * qJD(4) + t513 * t492 + t510 * t493;
t485 = (t510 * t514 + t511 * t513) * qJD(1);
t467 = -t484 * mrSges(6,1) + t485 * mrSges(6,2);
t468 = -t484 * mrSges(5,1) + t485 * mrSges(5,2);
t505 = qJD(3) + qJD(4);
t473 = -t505 * mrSges(5,2) + t484 * mrSges(5,3);
t504 = qJDD(3) + qJDD(4);
t433 = -0.2e1 * qJD(5) * t485 + (t484 * t505 - t454) * qJ(5) + (t484 * t485 + t504) * pkin(4) + t438;
t472 = -t505 * mrSges(6,2) + t484 * mrSges(6,3);
t528 = m(6) * t433 + t504 * mrSges(6,1) + t505 * t472;
t427 = m(5) * t438 + t504 * mrSges(5,1) + t505 * t473 + (-t467 - t468) * t485 + (-mrSges(5,3) - mrSges(6,3)) * t454 + t528;
t439 = t510 * t443 + t513 * t444;
t453 = -t485 * qJD(4) - t510 * t492 + t513 * t493;
t475 = t505 * mrSges(6,1) - t485 * mrSges(6,3);
t476 = t505 * mrSges(5,1) - t485 * mrSges(5,3);
t474 = t505 * pkin(4) - t485 * qJ(5);
t477 = t484 ^ 2;
t435 = -t477 * pkin(4) + t453 * qJ(5) + 0.2e1 * qJD(5) * t484 - t505 * t474 + t439;
t527 = m(6) * t435 + t453 * mrSges(6,3) + t484 * t467;
t430 = m(5) * t439 + t453 * mrSges(5,3) + t484 * t468 + (-t475 - t476) * t505 + (-mrSges(5,2) - mrSges(6,2)) * t504 + t527;
t423 = t513 * t427 + t510 * t430;
t490 = (-mrSges(4,1) * t514 + mrSges(4,2) * t511) * qJD(1);
t530 = qJD(1) * t514;
t495 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t530;
t421 = m(4) * t451 + qJDD(3) * mrSges(4,1) - t492 * mrSges(4,3) + qJD(3) * t495 - t490 * t531 + t423;
t494 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t531;
t521 = -t510 * t427 + t513 * t430;
t422 = m(4) * t452 - qJDD(3) * mrSges(4,2) + t493 * mrSges(4,3) - qJD(3) * t494 + t490 * t530 + t521;
t522 = -t511 * t421 + t514 * t422;
t414 = m(3) * t470 - t516 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t522;
t469 = t509 * t489 - t508 * t491;
t519 = -qJDD(1) * pkin(2) - t469;
t465 = -t516 * pkin(6) + t519;
t445 = -t493 * pkin(3) + t496 * t531 + (-pkin(7) * t506 - pkin(6)) * t516 + t519;
t437 = -t453 * pkin(4) - t477 * qJ(5) + t485 * t474 + qJDD(5) + t445;
t520 = m(6) * t437 - t453 * mrSges(6,1) + t454 * mrSges(6,2) - t484 * t472 + t485 * t475;
t518 = m(5) * t445 - t453 * mrSges(5,1) + t454 * mrSges(5,2) - t484 * t473 + t485 * t476 + t520;
t517 = -m(4) * t465 + t493 * mrSges(4,1) - t492 * mrSges(4,2) - t494 * t531 + t495 * t530 - t518;
t425 = m(3) * t469 + qJDD(1) * mrSges(3,1) - t516 * mrSges(3,2) + t517;
t523 = t509 * t414 - t508 * t425;
t409 = m(2) * t497 - t516 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t523;
t411 = t508 * t414 + t509 * t425;
t410 = m(2) * t498 + qJDD(1) * mrSges(2,1) - t516 * mrSges(2,2) + t411;
t535 = t512 * t409 + t515 * t410;
t415 = t514 * t421 + t511 * t422;
t534 = t536 * t484 - t537 * t485 + t539 * t505;
t533 = -t540 * t484 - t538 * t485 + t536 * t505;
t532 = t538 * t484 + t541 * t485 + t537 * t505;
t526 = m(3) * t507 + t415;
t524 = -t515 * t409 + t512 * t410;
t483 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t511 + Ifges(4,4) * t514) * qJD(1);
t482 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t511 + Ifges(4,2) * t514) * qJD(1);
t481 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t511 + Ifges(4,6) * t514) * qJD(1);
t431 = -t454 * mrSges(6,3) - t485 * t467 + t528;
t417 = mrSges(5,2) * t445 + mrSges(6,2) * t437 - mrSges(5,3) * t438 - mrSges(6,3) * t433 - qJ(5) * t431 + t538 * t453 + t541 * t454 - t534 * t484 + t537 * t504 + t533 * t505;
t416 = -mrSges(5,1) * t445 + mrSges(5,3) * t439 - mrSges(6,1) * t437 + mrSges(6,3) * t435 - pkin(4) * t520 + qJ(5) * t527 + (-qJ(5) * t475 + t532) * t505 + (-qJ(5) * mrSges(6,2) - t536) * t504 + t534 * t485 + t538 * t454 + t540 * t453;
t405 = mrSges(4,2) * t465 - mrSges(4,3) * t451 + Ifges(4,1) * t492 + Ifges(4,4) * t493 + Ifges(4,5) * qJDD(3) - pkin(7) * t423 - qJD(3) * t482 - t510 * t416 + t513 * t417 + t481 * t530;
t404 = -mrSges(4,1) * t465 + mrSges(4,3) * t452 + Ifges(4,4) * t492 + Ifges(4,2) * t493 + Ifges(4,6) * qJDD(3) - pkin(3) * t518 + pkin(7) * t521 + qJD(3) * t483 + t513 * t416 + t510 * t417 - t481 * t531;
t403 = -Ifges(4,3) * qJDD(3) + Ifges(3,6) * qJDD(1) + t516 * Ifges(3,5) - mrSges(3,1) * t507 - Ifges(4,5) * t492 - Ifges(4,6) * t493 + mrSges(3,3) * t470 - mrSges(4,1) * t451 + mrSges(4,2) * t452 - mrSges(5,1) * t438 + mrSges(5,2) * t439 - pkin(4) * t431 - mrSges(6,1) * t433 + mrSges(6,2) * t435 - pkin(3) * t423 - pkin(2) * t415 + t539 * t504 + t533 * t485 + t532 * t484 - t537 * t454 + t536 * t453 + (-t511 * t482 + t514 * t483) * qJD(1);
t402 = mrSges(3,2) * t507 - mrSges(3,3) * t469 + Ifges(3,5) * qJDD(1) - t516 * Ifges(3,6) - pkin(6) * t415 - t511 * t404 + t514 * t405;
t401 = -mrSges(2,2) * g(1) - mrSges(2,3) * t498 + Ifges(2,5) * qJDD(1) - t516 * Ifges(2,6) - qJ(2) * t411 + t509 * t402 - t508 * t403;
t400 = mrSges(2,1) * g(1) + mrSges(2,3) * t497 + t516 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t526 + qJ(2) * t523 + t508 * t402 + t509 * t403;
t1 = [(-m(1) - m(2)) * g(1) + t526; -m(1) * g(2) + t535; -m(1) * g(3) + t524; pkin(1) * t411 + t511 * t405 + t514 * t404 + pkin(2) * t517 + pkin(6) * t522 + mrSges(3,1) * t469 - mrSges(3,2) * t470 + mrSges(2,1) * t498 - mrSges(2,2) * t497 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t524 + t515 * t400 + t512 * t401; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t535 + t512 * t400 - t515 * t401;];
tauB = t1;
