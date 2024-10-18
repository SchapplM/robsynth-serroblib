% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% m [6x1]
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR13_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:50
% EndTime: 2024-09-27 17:30:59
% DurationCPUTime: 3.05s
% Computational Cost: add. (153428->243), mult. (179152->322), div. (0->0), fcn. (113846->12), ass. (0->112)
t539 = -m(3) - m(4);
t506 = sin(pkin(5));
t538 = pkin(9) * t506;
t507 = cos(pkin(5));
t537 = t507 * g(3);
t504 = qJD(1) + qJD(2);
t499 = qJD(3) + t504;
t497 = t499 ^ 2;
t536 = t497 * t506 ^ 2;
t535 = t499 * t506;
t509 = sin(qJ(4));
t534 = t506 * t509;
t514 = cos(qJ(4));
t533 = t506 * t514;
t532 = t507 * t509;
t531 = t507 * t514;
t502 = qJDD(1) + qJDD(2);
t498 = qJDD(3) + t502;
t529 = qJD(4) * t499;
t482 = (t498 * t509 + t514 * t529) * t506;
t491 = t507 * t498 + qJDD(4);
t492 = t507 * t499 + qJD(4);
t512 = sin(qJ(1));
t517 = cos(qJ(1));
t493 = t512 * g(1) - t517 * g(2);
t488 = qJDD(1) * pkin(1) + t493;
t494 = -t517 * g(1) - t512 * g(2);
t518 = qJD(1) ^ 2;
t489 = -t518 * pkin(1) + t494;
t511 = sin(qJ(2));
t516 = cos(qJ(2));
t474 = t516 * t488 - t511 * t489;
t472 = t502 * pkin(2) + t474;
t475 = t511 * t488 + t516 * t489;
t501 = t504 ^ 2;
t473 = -t501 * pkin(2) + t475;
t510 = sin(qJ(3));
t515 = cos(qJ(3));
t457 = t515 * t472 - t510 * t473;
t452 = t498 * pkin(3) + t497 * t538 + t457;
t458 = t510 * t472 + t515 * t473;
t453 = -t497 * pkin(3) + t498 * t538 + t458;
t521 = t452 * t531 - t509 * t453;
t442 = t491 * pkin(4) - t482 * pkin(10) + (pkin(4) * t509 * t536 + (pkin(10) * t492 * t499 - g(3)) * t506) * t514 + t521;
t445 = -g(3) * t534 + t452 * t532 + t514 * t453;
t527 = t499 * t534;
t481 = t492 * pkin(4) - pkin(10) * t527;
t483 = (t498 * t514 - t509 * t529) * t506;
t528 = t514 ^ 2 * t536;
t443 = -pkin(4) * t528 + t483 * pkin(10) - t492 * t481 + t445;
t508 = sin(qJ(5));
t513 = cos(qJ(5));
t440 = t513 * t442 - t508 * t443;
t476 = (-t508 * t509 + t513 * t514) * t535;
t456 = t476 * qJD(5) + t513 * t482 + t508 * t483;
t477 = (t508 * t514 + t509 * t513) * t535;
t463 = -t476 * mrSges(6,1) + t477 * mrSges(6,2);
t490 = qJD(5) + t492;
t464 = -t490 * mrSges(6,2) + t476 * mrSges(6,3);
t487 = qJDD(5) + t491;
t436 = m(6) * t440 + t487 * mrSges(6,1) - t456 * mrSges(6,3) - t477 * t463 + t490 * t464;
t441 = t508 * t442 + t513 * t443;
t455 = -t477 * qJD(5) - t508 * t482 + t513 * t483;
t465 = t490 * mrSges(6,1) - t477 * mrSges(6,3);
t437 = m(6) * t441 - t487 * mrSges(6,2) + t455 * mrSges(6,3) + t476 * t463 - t490 * t465;
t430 = t513 * t436 + t508 * t437;
t444 = -g(3) * t533 + t521;
t526 = t499 * t533;
t479 = -t492 * mrSges(5,2) + mrSges(5,3) * t526;
t480 = (-mrSges(5,1) * t514 + mrSges(5,2) * t509) * t535;
t428 = m(5) * t444 + t491 * mrSges(5,1) - t482 * mrSges(5,3) + t492 * t479 - t480 * t527 + t430;
t478 = t492 * mrSges(5,1) - mrSges(5,3) * t527;
t522 = -t508 * t436 + t513 * t437;
t429 = m(5) * t445 - t491 * mrSges(5,2) + t483 * mrSges(5,3) - t492 * t478 + t480 * t526 + t522;
t448 = -t506 * t452 - t537;
t447 = -pkin(10) * t528 - t483 * pkin(4) - t537 + (t481 * t499 * t509 - t452) * t506;
t519 = m(6) * t447 - t455 * mrSges(6,1) + t456 * mrSges(6,2) - t476 * t464 + t477 * t465;
t439 = m(5) * t448 - t483 * mrSges(5,1) + t482 * mrSges(5,2) + (t478 * t509 - t479 * t514) * t535 + t519;
t416 = t428 * t531 + t429 * t532 - t506 * t439;
t414 = m(4) * t457 + t498 * mrSges(4,1) - t497 * mrSges(4,2) + t416;
t421 = -t509 * t428 + t514 * t429;
t420 = m(4) * t458 - t497 * mrSges(4,1) - t498 * mrSges(4,2) + t421;
t411 = t515 * t414 + t510 * t420;
t409 = m(3) * t474 + t502 * mrSges(3,1) - t501 * mrSges(3,2) + t411;
t523 = -t510 * t414 + t515 * t420;
t410 = m(3) * t475 - t501 * mrSges(3,1) - t502 * mrSges(3,2) + t523;
t405 = t516 * t409 + t511 * t410;
t403 = m(2) * t493 + qJDD(1) * mrSges(2,1) - t518 * mrSges(2,2) + t405;
t524 = -t511 * t409 + t516 * t410;
t404 = m(2) * t494 - t518 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t524;
t530 = t517 * t403 + t512 * t404;
t415 = t428 * t533 + t429 * t534 + t507 * t439;
t525 = -t512 * t403 + t517 * t404;
t459 = Ifges(6,5) * t477 + Ifges(6,6) * t476 + Ifges(6,3) * t490;
t461 = Ifges(6,1) * t477 + Ifges(6,4) * t476 + Ifges(6,5) * t490;
t431 = -mrSges(6,1) * t447 + mrSges(6,3) * t441 + Ifges(6,4) * t456 + Ifges(6,2) * t455 + Ifges(6,6) * t487 - t477 * t459 + t490 * t461;
t460 = Ifges(6,4) * t477 + Ifges(6,2) * t476 + Ifges(6,6) * t490;
t432 = mrSges(6,2) * t447 - mrSges(6,3) * t440 + Ifges(6,1) * t456 + Ifges(6,4) * t455 + Ifges(6,5) * t487 + t476 * t459 - t490 * t460;
t469 = Ifges(5,3) * t492 + (Ifges(5,5) * t509 + Ifges(5,6) * t514) * t535;
t471 = Ifges(5,5) * t492 + (Ifges(5,1) * t509 + Ifges(5,4) * t514) * t535;
t412 = -mrSges(5,1) * t448 + mrSges(5,3) * t445 + Ifges(5,4) * t482 + Ifges(5,2) * t483 + Ifges(5,6) * t491 - pkin(4) * t519 + pkin(10) * t522 + t513 * t431 + t508 * t432 - t469 * t527 + t492 * t471;
t470 = Ifges(5,6) * t492 + (Ifges(5,4) * t509 + Ifges(5,2) * t514) * t535;
t417 = mrSges(5,2) * t448 - mrSges(5,3) * t444 + Ifges(5,1) * t482 + Ifges(5,4) * t483 + Ifges(5,5) * t491 - pkin(10) * t430 - t508 * t431 + t513 * t432 + t469 * t526 - t492 * t470;
t520 = pkin(9) * t421 + t514 * t412 + t509 * t417;
t422 = mrSges(5,1) * t444 + mrSges(6,1) * t440 - mrSges(5,2) * t445 - mrSges(6,2) * t441 + Ifges(5,5) * t482 + Ifges(6,5) * t456 + Ifges(5,6) * t483 + Ifges(6,6) * t455 + Ifges(5,3) * t491 + Ifges(6,3) * t487 + pkin(4) * t430 + t477 * t460 - t476 * t461 + (t470 * t509 - t471 * t514) * t535;
t399 = -mrSges(4,2) * g(3) - mrSges(4,3) * t457 + Ifges(4,5) * t498 - t497 * Ifges(4,6) - t509 * t412 + t514 * t417 + (-t415 * t506 - t416 * t507) * pkin(9);
t398 = mrSges(4,1) * g(3) + mrSges(4,3) * t458 + t497 * Ifges(4,5) + Ifges(4,6) * t498 - pkin(3) * t415 - t506 * t422 + t520 * t507;
t397 = -mrSges(3,2) * g(3) - mrSges(3,3) * t474 + Ifges(3,5) * t502 - t501 * Ifges(3,6) - pkin(8) * t411 - t510 * t398 + t515 * t399;
t396 = Ifges(3,6) * t502 + t501 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t475 + t510 * t399 + t515 * t398 - pkin(2) * (-m(4) * g(3) + t415) + pkin(8) * t523;
t395 = -mrSges(2,2) * g(3) - mrSges(2,3) * t493 + Ifges(2,5) * qJDD(1) - t518 * Ifges(2,6) - pkin(7) * t405 - t511 * t396 + t516 * t397;
t394 = Ifges(2,6) * qJDD(1) + t518 * Ifges(2,5) + mrSges(2,3) * t494 + t511 * t397 + t516 * t396 - pkin(1) * t415 + pkin(7) * t524 + (-pkin(1) * t539 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t525; -m(1) * g(2) + t530; (-m(1) - m(2) + t539) * g(3) + t415; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(6) * t530 - t512 * t394 + t517 * t395; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(6) * t525 + t517 * t394 + t512 * t395; -mrSges(1,1) * g(2) + mrSges(2,1) * t493 + mrSges(3,1) * t474 + mrSges(4,1) * t457 + mrSges(1,2) * g(1) - mrSges(2,2) * t494 - mrSges(3,2) * t475 - mrSges(4,2) * t458 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t502 + Ifges(4,3) * t498 + pkin(1) * t405 + pkin(2) * t411 + pkin(3) * t416 + t507 * t422 + t520 * t506;];
tauB = t1;
