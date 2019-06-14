% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6PRRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-05-05 05:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6PRRPRR4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR4_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:06:51
% EndTime: 2019-05-05 05:06:55
% DurationCPUTime: 2.12s
% Computational Cost: add. (12686->271), mult. (25126->335), div. (0->0), fcn. (16174->12), ass. (0->119)
t509 = sin(qJ(3));
t513 = cos(qJ(3));
t543 = Ifges(4,4) - Ifges(5,5);
t551 = t509 * (Ifges(4,1) + Ifges(5,1)) + t513 * t543;
t550 = -t509 * t543 - t513 * (Ifges(4,2) + Ifges(5,3));
t542 = Ifges(4,5) + Ifges(5,4);
t541 = Ifges(4,6) - Ifges(5,6);
t502 = sin(pkin(11));
t504 = cos(pkin(11));
t480 = g(1) * t502 - g(2) * t504;
t500 = -g(3) + qJDD(1);
t503 = sin(pkin(6));
t505 = cos(pkin(6));
t548 = t480 * t505 + t500 * t503;
t508 = sin(qJ(5));
t512 = cos(qJ(5));
t461 = (t508 * t509 + t512 * t513) * qJD(2);
t546 = t509 * (qJD(2) * t550 - t541 * qJD(3)) + t513 * (qJD(2) * t551 + t542 * qJD(3));
t545 = 2 * qJD(4);
t544 = mrSges(4,3) + mrSges(5,2);
t531 = qJD(2) * qJD(3);
t529 = t509 * t531;
t477 = qJDD(2) * t513 - t529;
t540 = mrSges(5,1) * t477;
t448 = -t480 * t503 + t500 * t505;
t539 = t448 * t513;
t516 = qJD(2) ^ 2;
t537 = t513 ^ 2 * t516;
t481 = -g(1) * t504 - g(2) * t502;
t510 = sin(qJ(2));
t514 = cos(qJ(2));
t438 = t514 * t481 + t548 * t510;
t433 = -pkin(2) * t516 + qJDD(2) * pkin(8) + t438;
t421 = t513 * t433 + t509 * t448;
t473 = (-t513 * pkin(3) - t509 * qJ(4)) * qJD(2);
t515 = qJD(3) ^ 2;
t532 = qJD(2) * t513;
t410 = -pkin(3) * t515 + qJDD(3) * qJ(4) + qJD(3) * t545 + t473 * t532 + t421;
t533 = qJD(2) * t509;
t488 = -qJD(3) * pkin(4) - pkin(9) * t533;
t405 = -pkin(4) * t537 - pkin(9) * t477 + qJD(3) * t488 + t410;
t530 = t513 * t531;
t476 = qJDD(2) * t509 + t530;
t426 = t509 * t433;
t525 = -qJ(4) * t515 + t473 * t533 + qJDD(4) + t426;
t406 = -pkin(9) * t476 + (-pkin(3) - pkin(4)) * qJDD(3) + (-pkin(4) * t509 * t516 + pkin(9) * t531 - t448) * t513 + t525;
t401 = t512 * t405 + t508 * t406;
t462 = (-t513 * t508 + t509 * t512) * qJD(2);
t442 = pkin(5) * t461 - pkin(10) * t462;
t495 = -qJD(3) + qJD(5);
t493 = t495 ^ 2;
t494 = -qJDD(3) + qJDD(5);
t399 = -pkin(5) * t493 + pkin(10) * t494 - t442 * t461 + t401;
t437 = -t510 * t481 + t548 * t514;
t432 = -qJDD(2) * pkin(2) - t516 * pkin(8) - t437;
t522 = -t477 * pkin(3) + t432 + (-t476 - t530) * qJ(4);
t408 = -pkin(3) * t529 + pkin(4) * t477 - pkin(9) * t537 - t522 + (t488 + t545) * t533;
t428 = -qJD(5) * t462 - t476 * t508 - t477 * t512;
t429 = -qJD(5) * t461 + t476 * t512 - t477 * t508;
t402 = (t461 * t495 - t429) * pkin(10) + (t462 * t495 - t428) * pkin(5) + t408;
t507 = sin(qJ(6));
t511 = cos(qJ(6));
t396 = -t399 * t507 + t402 * t511;
t443 = -t462 * t507 + t495 * t511;
t415 = qJD(6) * t443 + t429 * t511 + t494 * t507;
t444 = t462 * t511 + t495 * t507;
t422 = -mrSges(7,1) * t443 + mrSges(7,2) * t444;
t425 = qJDD(6) - t428;
t452 = qJD(6) + t461;
t430 = -mrSges(7,2) * t452 + mrSges(7,3) * t443;
t393 = m(7) * t396 + mrSges(7,1) * t425 - mrSges(7,3) * t415 - t422 * t444 + t430 * t452;
t397 = t399 * t511 + t402 * t507;
t414 = -qJD(6) * t444 - t429 * t507 + t494 * t511;
t431 = mrSges(7,1) * t452 - mrSges(7,3) * t444;
t394 = m(7) * t397 - mrSges(7,2) * t425 + mrSges(7,3) * t414 + t422 * t443 - t431 * t452;
t385 = t511 * t393 + t507 * t394;
t475 = (-t513 * mrSges(4,1) + t509 * mrSges(4,2)) * qJD(2);
t482 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t533;
t474 = (-t513 * mrSges(5,1) - t509 * mrSges(5,3)) * qJD(2);
t483 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t533;
t386 = -t393 * t507 + t511 * t394;
t441 = mrSges(6,1) * t461 + mrSges(6,2) * t462;
t447 = mrSges(6,1) * t495 - mrSges(6,3) * t462;
t384 = m(6) * t401 - mrSges(6,2) * t494 + mrSges(6,3) * t428 - t441 * t461 - t447 * t495 + t386;
t400 = -t405 * t508 + t406 * t512;
t398 = -pkin(5) * t494 - pkin(10) * t493 + t442 * t462 - t400;
t395 = -m(7) * t398 + t414 * mrSges(7,1) - mrSges(7,2) * t415 + t443 * t430 - t431 * t444;
t446 = -mrSges(6,2) * t495 - mrSges(6,3) * t461;
t389 = m(6) * t400 + mrSges(6,1) * t494 - mrSges(6,3) * t429 - t441 * t462 + t446 * t495 + t395;
t527 = t512 * t384 - t389 * t508;
t523 = m(5) * t410 + qJDD(3) * mrSges(5,3) + qJD(3) * t483 + t474 * t532 + t527;
t378 = m(4) * t421 - qJDD(3) * mrSges(4,2) - qJD(3) * t482 + t475 * t532 + t544 * t477 + t523;
t420 = -t426 + t539;
t484 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t532;
t381 = t384 * t508 + t389 * t512;
t411 = -qJDD(3) * pkin(3) + t525 - t539;
t485 = mrSges(5,2) * t532 + qJD(3) * mrSges(5,3);
t521 = -m(5) * t411 + qJDD(3) * mrSges(5,1) + qJD(3) * t485 - t381;
t379 = m(4) * t420 + qJDD(3) * mrSges(4,1) + qJD(3) * t484 - t544 * t476 + (-t474 - t475) * t533 + t521;
t528 = t513 * t378 - t379 * t509;
t524 = -m(6) * t408 + t428 * mrSges(6,1) - t429 * mrSges(6,2) - t461 * t446 - t462 * t447 - t385;
t412 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t533 + t522;
t520 = m(5) * t412 - t476 * mrSges(5,3) - t483 * t533 - t485 * t532 + t524;
t417 = Ifges(7,4) * t444 + Ifges(7,2) * t443 + Ifges(7,6) * t452;
t418 = Ifges(7,1) * t444 + Ifges(7,4) * t443 + Ifges(7,5) * t452;
t519 = mrSges(7,1) * t396 - mrSges(7,2) * t397 + Ifges(7,5) * t415 + Ifges(7,6) * t414 + Ifges(7,3) * t425 + t417 * t444 - t443 * t418;
t518 = -m(4) * t432 + t477 * mrSges(4,1) - t482 * t533 + t484 * t532 - t520;
t416 = Ifges(7,5) * t444 + Ifges(7,6) * t443 + Ifges(7,3) * t452;
t387 = -mrSges(7,1) * t398 + mrSges(7,3) * t397 + Ifges(7,4) * t415 + Ifges(7,2) * t414 + Ifges(7,6) * t425 - t416 * t444 + t418 * t452;
t388 = mrSges(7,2) * t398 - mrSges(7,3) * t396 + Ifges(7,1) * t415 + Ifges(7,4) * t414 + Ifges(7,5) * t425 + t416 * t443 - t417 * t452;
t435 = Ifges(6,4) * t462 - Ifges(6,2) * t461 + Ifges(6,6) * t495;
t436 = Ifges(6,1) * t462 - Ifges(6,4) * t461 + Ifges(6,5) * t495;
t517 = mrSges(6,1) * t400 - mrSges(6,2) * t401 + Ifges(6,5) * t429 + Ifges(6,6) * t428 + Ifges(6,3) * t494 + pkin(5) * t395 + pkin(10) * t386 + t511 * t387 + t507 * t388 + t462 * t435 + t461 * t436;
t434 = Ifges(6,5) * t462 - Ifges(6,6) * t461 + Ifges(6,3) * t495;
t382 = t520 - t540;
t380 = mrSges(5,2) * t476 + t474 * t533 - t521;
t376 = -mrSges(6,1) * t408 + mrSges(6,3) * t401 + Ifges(6,4) * t429 + Ifges(6,2) * t428 + Ifges(6,6) * t494 - pkin(5) * t385 - t434 * t462 + t436 * t495 - t519;
t375 = mrSges(6,2) * t408 - mrSges(6,3) * t400 + Ifges(6,1) * t429 + Ifges(6,4) * t428 + Ifges(6,5) * t494 - pkin(10) * t385 - t387 * t507 + t388 * t511 - t434 * t461 - t435 * t495;
t1 = [m(2) * t500 + t505 * (m(3) * t448 + t378 * t509 + t379 * t513) + (t510 * (m(3) * t438 - mrSges(3,1) * t516 - qJDD(2) * mrSges(3,2) + t528) + t514 * (m(3) * t437 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t516 - mrSges(4,2) * t476 + t518 + t540)) * t503; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t437 - mrSges(3,2) * t438 + t509 * (mrSges(4,2) * t432 + mrSges(5,2) * t411 - mrSges(4,3) * t420 - mrSges(5,3) * t412 - pkin(9) * t381 - qJ(4) * t382 + t512 * t375 - t508 * t376) + t513 * (-mrSges(4,1) * t432 - mrSges(5,1) * t412 + mrSges(5,2) * t410 + mrSges(4,3) * t421 - pkin(3) * t382 - pkin(4) * t524 - pkin(9) * t527 - t508 * t375 - t512 * t376) + pkin(2) * t518 + pkin(8) * t528 + (pkin(2) * mrSges(5,1) - t550) * t477 + (-pkin(2) * mrSges(4,2) + t551) * t476 + (t509 * t542 + t513 * t541) * qJDD(3) + t546 * qJD(3); qJ(4) * t523 + (qJ(4) * mrSges(5,2) + t541) * t477 - pkin(4) * t381 - pkin(3) * t380 + mrSges(4,1) * t420 - mrSges(4,2) * t421 + mrSges(5,3) * t410 - mrSges(5,1) * t411 + (Ifges(4,3) + Ifges(5,2)) * qJDD(3) + t542 * t476 - t546 * qJD(2) - t517; t380; t517; t519;];
tauJ  = t1;
