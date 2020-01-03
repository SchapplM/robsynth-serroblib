% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RPPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RPPRP6_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:05
% EndTime: 2019-12-31 17:55:06
% DurationCPUTime: 1.33s
% Computational Cost: add. (9412->216), mult. (20533->250), div. (0->0), fcn. (11988->6), ass. (0->94)
t535 = Ifges(5,1) + Ifges(6,1);
t526 = Ifges(5,4) - Ifges(6,5);
t525 = Ifges(5,5) + Ifges(6,4);
t534 = Ifges(5,2) + Ifges(6,3);
t523 = Ifges(5,6) - Ifges(6,6);
t533 = Ifges(5,3) + Ifges(6,2);
t485 = sin(qJ(1));
t487 = cos(qJ(1));
t464 = g(1) * t485 - t487 * g(2);
t489 = qJD(1) ^ 2;
t495 = -qJ(2) * t489 + qJDD(2) - t464;
t522 = -pkin(1) - qJ(3);
t532 = -(2 * qJD(1) * qJD(3)) + t522 * qJDD(1) + t495;
t482 = sin(pkin(7));
t474 = t482 ^ 2;
t483 = cos(pkin(7));
t515 = t483 ^ 2 + t474;
t507 = t515 * mrSges(4,3);
t465 = -g(1) * t487 - g(2) * t485;
t531 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t465;
t530 = pkin(3) * t489;
t529 = mrSges(2,1) - mrSges(3,2);
t528 = -mrSges(5,3) - mrSges(6,2);
t527 = -Ifges(3,4) + Ifges(2,5);
t524 = -Ifges(2,6) + Ifges(3,5);
t521 = mrSges(4,2) * t483;
t443 = t482 * g(3) + t532 * t483;
t426 = (-pkin(6) * qJDD(1) - t482 * t530) * t483 + t443;
t444 = -g(3) * t483 + t532 * t482;
t510 = qJDD(1) * t482;
t427 = -pkin(6) * t510 - t474 * t530 + t444;
t484 = sin(qJ(4));
t486 = cos(qJ(4));
t423 = t484 * t426 + t486 * t427;
t499 = t482 * t486 + t483 * t484;
t498 = -t482 * t484 + t483 * t486;
t461 = t498 * qJD(1);
t512 = qJD(4) * t461;
t445 = t499 * qJDD(1) + t512;
t454 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t461;
t460 = t499 * qJD(1);
t438 = pkin(4) * t460 - qJ(5) * t461;
t488 = qJD(4) ^ 2;
t420 = -pkin(4) * t488 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t438 * t460 + t423;
t455 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t461;
t508 = m(6) * t420 + qJDD(4) * mrSges(6,3) + qJD(4) * t455;
t439 = mrSges(6,1) * t460 - mrSges(6,3) * t461;
t516 = -mrSges(5,1) * t460 - mrSges(5,2) * t461 - t439;
t414 = m(5) * t423 - qJDD(4) * mrSges(5,2) - qJD(4) * t454 + t528 * t445 + t516 * t460 + t508;
t422 = t426 * t486 - t427 * t484;
t513 = qJD(4) * t460;
t446 = t498 * qJDD(1) - t513;
t453 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t460;
t421 = -qJDD(4) * pkin(4) - qJ(5) * t488 + t438 * t461 + qJDD(5) - t422;
t456 = -mrSges(6,2) * t460 + qJD(4) * mrSges(6,3);
t503 = -m(6) * t421 + qJDD(4) * mrSges(6,1) + qJD(4) * t456;
t415 = m(5) * t422 + qJDD(4) * mrSges(5,1) + qJD(4) * t453 + t528 * t446 + t516 * t461 + t503;
t407 = t484 * t414 + t486 * t415;
t497 = -qJDD(1) * mrSges(4,3) - t489 * (mrSges(4,1) * t482 + t521);
t405 = m(4) * t443 + t497 * t483 + t407;
t504 = t486 * t414 - t415 * t484;
t406 = m(4) * t444 + t497 * t482 + t504;
t401 = t405 * t483 + t406 * t482;
t459 = -qJDD(1) * pkin(1) + t495;
t493 = -m(3) * t459 + t489 * mrSges(3,3) - t401;
t399 = m(2) * t464 - mrSges(2,2) * t489 + t529 * qJDD(1) + t493;
t458 = pkin(1) * t489 - t531;
t494 = qJDD(3) + t531;
t450 = t522 * t489 + t494;
t429 = pkin(3) * t510 + (-t515 * pkin(6) + t522) * t489 + t494;
t418 = -0.2e1 * qJD(5) * t461 + (-t446 + t513) * qJ(5) + (t445 + t512) * pkin(4) + t429;
t416 = m(6) * t418 + t445 * mrSges(6,1) - mrSges(6,3) * t446 - t455 * t461 + t460 * t456;
t492 = m(5) * t429 + mrSges(5,1) * t445 + t446 * mrSges(5,2) + t453 * t460 + t461 * t454 + t416;
t491 = m(4) * t450 + mrSges(4,1) * t510 + qJDD(1) * t521 + t492;
t490 = -m(3) * t458 + t489 * mrSges(3,2) + qJDD(1) * mrSges(3,3) + t491;
t410 = -qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t507) * t489 + m(2) * t465 + t490;
t520 = t487 * t399 + t485 * t410;
t519 = -t523 * qJD(4) + t534 * t460 - t526 * t461;
t518 = -t533 * qJD(4) + t523 * t460 - t525 * t461;
t517 = t525 * qJD(4) - t526 * t460 + t535 * t461;
t500 = Ifges(4,5) * t483 - Ifges(4,6) * t482;
t514 = t489 * t500;
t506 = -t399 * t485 + t487 * t410;
t505 = -t482 * t405 + t483 * t406;
t502 = Ifges(4,1) * t483 - Ifges(4,4) * t482;
t501 = Ifges(4,4) * t483 - Ifges(4,2) * t482;
t403 = mrSges(5,2) * t429 + mrSges(6,2) * t421 - mrSges(5,3) * t422 - mrSges(6,3) * t418 - qJ(5) * t416 + t519 * qJD(4) + t525 * qJDD(4) - t526 * t445 + t535 * t446 + t518 * t460;
t402 = -mrSges(5,1) * t429 - mrSges(6,1) * t418 + mrSges(6,2) * t420 + mrSges(5,3) * t423 - pkin(4) * t416 + t517 * qJD(4) + t523 * qJDD(4) - t534 * t445 + t526 * t446 + t518 * t461;
t400 = -m(3) * g(3) + t505;
t397 = mrSges(4,2) * t450 - mrSges(4,3) * t443 - pkin(6) * t407 + t502 * qJDD(1) - t402 * t484 + t403 * t486 - t482 * t514;
t396 = -mrSges(4,1) * t450 + mrSges(4,3) * t444 - pkin(3) * t492 + pkin(6) * t504 + t501 * qJDD(1) + t486 * t402 + t484 * t403 - t483 * t514;
t395 = -mrSges(2,3) * t464 + qJ(5) * t508 + pkin(4) * t503 + mrSges(3,1) * t459 + mrSges(4,1) * t443 - mrSges(4,2) * t444 + mrSges(6,3) * t420 - mrSges(6,1) * t421 + mrSges(5,1) * t422 - mrSges(5,2) * t423 + pkin(3) * t407 + pkin(2) * t401 - qJ(2) * t400 + (-pkin(4) * t439 - t519) * t461 + (-qJ(5) * t439 + t517) * t460 + (-pkin(4) * mrSges(6,2) + t525) * t446 + (-qJ(5) * mrSges(6,2) - t523) * t445 + t533 * qJDD(4) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (t500 + t527) * qJDD(1) + (t482 * t502 + t483 * t501 + t524) * t489;
t394 = mrSges(2,3) * t465 - mrSges(3,1) * t458 - t482 * t397 - t483 * t396 + pkin(2) * t491 - qJ(3) * t505 - pkin(1) * t400 - t524 * qJDD(1) + t529 * g(3) + (-pkin(2) * t507 + t527) * t489;
t1 = [-m(1) * g(1) + t506; -m(1) * g(2) + t520; (-m(1) - m(2) - m(3)) * g(3) + t505; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t520 - t485 * t394 + t487 * t395; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t506 + t487 * t394 + t485 * t395; pkin(1) * t493 + qJ(2) * (-t489 * t507 + t490) + t483 * t397 - t482 * t396 - qJ(3) * t401 + mrSges(2,1) * t464 - mrSges(2,2) * t465 + mrSges(3,2) * t459 - mrSges(3,3) * t458 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
