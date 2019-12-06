% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP1_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:48
% EndTime: 2019-12-05 17:35:51
% DurationCPUTime: 2.07s
% Computational Cost: add. (15273->226), mult. (31695->283), div. (0->0), fcn. (17540->8), ass. (0->101)
t497 = sin(qJ(1));
t499 = cos(qJ(1));
t476 = t499 * g(2) + t497 * g(3);
t472 = qJDD(1) * pkin(1) + t476;
t475 = t497 * g(2) - t499 * g(3);
t500 = qJD(1) ^ 2;
t473 = -t500 * pkin(1) + t475;
t493 = sin(pkin(7));
t495 = cos(pkin(7));
t445 = t493 * t472 + t495 * t473;
t545 = -t500 * pkin(2) + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t445;
t544 = Ifges(5,1) + Ifges(6,1);
t537 = Ifges(5,4) + Ifges(6,4);
t536 = Ifges(5,5) + Ifges(6,5);
t543 = Ifges(5,2) + Ifges(6,2);
t542 = Ifges(5,6) + Ifges(6,6);
t541 = -Ifges(5,3) - Ifges(6,3);
t491 = -g(1) + qJDD(2);
t492 = sin(pkin(8));
t494 = cos(pkin(8));
t435 = t494 * t491 - t545 * t492;
t508 = -pkin(3) * t494 - pkin(6) * t492;
t471 = t508 * qJD(1);
t526 = t492 * qJD(1);
t433 = t471 * t526 - t435;
t496 = sin(qJ(4));
t498 = cos(qJ(4));
t523 = qJD(1) * qJD(4);
t463 = (-qJDD(1) * t496 - t498 * t523) * t492;
t464 = (qJDD(1) * t498 - t496 * t523) * t492;
t525 = t494 * qJD(1);
t478 = qJD(4) - t525;
t517 = t498 * t526;
t458 = t478 * pkin(4) - qJ(5) * t517;
t533 = t492 ^ 2 * t500;
t521 = t496 ^ 2 * t533;
t431 = -t463 * pkin(4) - qJ(5) * t521 + t458 * t517 + qJDD(5) + t433;
t515 = m(6) * t431 - t463 * mrSges(6,1);
t538 = -mrSges(5,2) - mrSges(6,2);
t540 = -m(5) * t433 + t463 * mrSges(5,1) + t538 * t464 - t515;
t539 = t494 ^ 2;
t534 = mrSges(4,2) * t492;
t436 = t492 * t491 + t545 * t494;
t469 = (-mrSges(4,1) * t494 + t534) * qJD(1);
t434 = t471 * t525 + t436;
t444 = t495 * t472 - t493 * t473;
t503 = -t500 * qJ(3) + qJDD(3) - t444;
t439 = (-pkin(2) + t508) * qJDD(1) + t503;
t438 = t498 * t439;
t429 = -t496 * t434 + t438;
t518 = t496 * t526;
t457 = -t478 * mrSges(5,2) - mrSges(5,3) * t518;
t522 = t494 * qJDD(1);
t477 = qJDD(4) - t522;
t461 = (mrSges(6,1) * t496 + mrSges(6,2) * t498) * t526;
t507 = (-t461 - (mrSges(5,1) * t496 + mrSges(5,2) * t498) * t526) * t526;
t509 = -0.2e1 * qJD(5) * t526;
t426 = t498 * t509 + t477 * pkin(4) - t464 * qJ(5) + t438 + (-pkin(4) * t498 * t533 - qJ(5) * t478 * t526 - t434) * t496;
t456 = -t478 * mrSges(6,2) - mrSges(6,3) * t518;
t520 = m(6) * t426 + t477 * mrSges(6,1) + t478 * t456;
t420 = m(5) * t429 + t477 * mrSges(5,1) + t478 * t457 + (-mrSges(5,3) - mrSges(6,3)) * t464 + t498 * t507 + t520;
t430 = t498 * t434 + t496 * t439;
t459 = t478 * mrSges(6,1) - mrSges(6,3) * t517;
t528 = -t478 * mrSges(5,1) + mrSges(5,3) * t517 - t459;
t428 = -pkin(4) * t521 + t463 * qJ(5) - t478 * t458 + t496 * t509 + t430;
t532 = m(6) * t428 + t463 * mrSges(6,3);
t421 = m(5) * t430 + t463 * mrSges(5,3) + t538 * t477 + t528 * t478 + t496 * t507 + t532;
t511 = -t496 * t420 + t498 * t421;
t527 = qJDD(1) * mrSges(4,3);
t416 = m(4) * t436 + (qJD(1) * t469 + t527) * t494 + t511;
t502 = t528 * t498 + (-t456 - t457) * t496;
t423 = m(4) * t435 + (-t527 + (-t469 + t502) * qJD(1)) * t492 + t540;
t512 = t494 * t416 - t492 * t423;
t409 = m(3) * t445 - t500 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t512;
t418 = t498 * t420 + t496 * t421;
t442 = -qJDD(1) * pkin(2) + t503;
t501 = -m(4) * t442 + mrSges(4,1) * t522 - t418 + (t500 * t539 + t533) * mrSges(4,3);
t413 = m(3) * t444 - t500 * mrSges(3,2) + (mrSges(3,1) - t534) * qJDD(1) + t501;
t405 = t493 * t409 + t495 * t413;
t410 = t492 * t416 + t494 * t423;
t531 = (t542 * t496 - t536 * t498) * t526 + t541 * t478;
t530 = (t543 * t496 - t537 * t498) * t526 - t542 * t478;
t529 = (t537 * t496 - t544 * t498) * t526 - t536 * t478;
t519 = m(3) * t491 + t410;
t513 = t495 * t409 - t493 * t413;
t403 = m(2) * t475 - t500 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t513;
t404 = m(2) * t476 + qJDD(1) * mrSges(2,1) - t500 * mrSges(2,2) + t405;
t514 = t499 * t403 - t497 * t404;
t506 = Ifges(4,5) * t492 + Ifges(4,6) * t494;
t505 = -t497 * t403 - t499 * t404;
t470 = t506 * qJD(1);
t424 = -t464 * mrSges(6,3) - t461 * t517 + t520;
t417 = mrSges(5,2) * t433 + mrSges(6,2) * t431 - mrSges(5,3) * t429 - mrSges(6,3) * t426 - qJ(5) * t424 + t537 * t463 + t544 * t464 + t536 * t477 + t530 * t478 + t531 * t518;
t411 = -mrSges(5,1) * t433 + mrSges(5,3) * t430 - mrSges(6,1) * t431 + mrSges(6,3) * t428 - pkin(4) * t515 + qJ(5) * t532 + (-qJ(5) * t459 - t529) * t478 + (-qJ(5) * mrSges(6,2) + t542) * t477 + (-pkin(4) * mrSges(6,2) + t537) * t464 + t543 * t463 + ((-pkin(4) * t456 - qJ(5) * t461) * t496 + (-pkin(4) * t459 + t531) * t498) * t526;
t406 = Ifges(4,2) * t522 - mrSges(4,1) * t442 - mrSges(5,1) * t429 - mrSges(6,1) * t426 + mrSges(5,2) * t430 + mrSges(6,2) * t428 + mrSges(4,3) * t436 - pkin(3) * t418 - pkin(4) * t424 + t541 * t477 - t536 * t464 - t542 * t463 + (Ifges(4,4) * qJDD(1) + (t529 * t496 + t530 * t498 - t470) * qJD(1)) * t492;
t401 = t470 * t525 + mrSges(4,2) * t442 - mrSges(4,3) * t435 - pkin(6) * t418 - t496 * t411 + t498 * t417 + (Ifges(4,1) * t492 + Ifges(4,4) * t494) * qJDD(1);
t400 = t500 * Ifges(3,5) - mrSges(3,1) * t491 + mrSges(3,3) * t445 - mrSges(4,1) * t435 + mrSges(4,2) * t436 - t496 * t417 - t498 * t411 - pkin(3) * t540 - pkin(6) * t511 - pkin(2) * t410 + (Ifges(3,6) - t506) * qJDD(1) + (Ifges(4,4) * t539 * qJD(1) + (-pkin(3) * t502 + (-Ifges(4,4) * t492 + (Ifges(4,1) - Ifges(4,2)) * t494) * qJD(1)) * t492) * qJD(1);
t399 = mrSges(3,2) * t491 - mrSges(3,3) * t444 + Ifges(3,5) * qJDD(1) - t500 * Ifges(3,6) - qJ(3) * t410 + t494 * t401 - t492 * t406;
t398 = -mrSges(2,2) * g(1) - mrSges(2,3) * t476 + Ifges(2,5) * qJDD(1) - t500 * Ifges(2,6) - qJ(2) * t405 + t495 * t399 - t493 * t400;
t397 = mrSges(2,1) * g(1) + mrSges(2,3) * t475 + t500 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t519 + qJ(2) * t513 + t493 * t399 + t495 * t400;
t1 = [(-m(1) - m(2)) * g(1) + t519; -m(1) * g(2) + t505; -m(1) * g(3) + t514; pkin(1) * t405 + qJ(3) * t512 + t492 * t401 + t494 * t406 + pkin(2) * t501 + mrSges(3,1) * t444 - mrSges(3,2) * t445 + mrSges(2,1) * t476 - mrSges(2,2) * t475 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (-pkin(2) * t534 + Ifges(2,3) + Ifges(3,3)) * qJDD(1); mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t514 - t499 * t397 - t497 * t398; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t505 - t497 * t397 + t499 * t398;];
tauB = t1;
