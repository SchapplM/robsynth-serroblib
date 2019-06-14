% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RRPRPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-05-06 13:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RRPRPR2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR2_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 13:02:09
% EndTime: 2019-05-06 13:02:18
% DurationCPUTime: 5.34s
% Computational Cost: add. (43931->329), mult. (103355->403), div. (0->0), fcn. (75477->10), ass. (0->133)
t567 = Ifges(5,4) + Ifges(6,6);
t579 = -Ifges(5,2) - Ifges(6,3);
t575 = Ifges(5,6) - Ifges(6,5);
t535 = sin(qJ(2));
t538 = cos(qJ(2));
t559 = qJD(1) * qJD(2);
t519 = qJDD(1) * t538 - t535 * t559;
t561 = qJD(1) * t535;
t520 = qJD(2) * pkin(2) - qJ(3) * t561;
t530 = t538 ^ 2;
t540 = qJD(1) ^ 2;
t536 = sin(qJ(1));
t539 = cos(qJ(1));
t557 = t536 * g(1) - g(2) * t539;
t551 = -qJDD(1) * pkin(1) - t557;
t477 = -t519 * pkin(2) + qJDD(3) + t520 * t561 + (-qJ(3) * t530 - pkin(7)) * t540 + t551;
t518 = qJDD(1) * t535 + t538 * t559;
t531 = sin(pkin(10));
t532 = cos(pkin(10));
t497 = -t518 * t531 + t519 * t532;
t510 = (t531 * t538 + t532 * t535) * qJD(1);
t503 = qJD(2) * pkin(3) - pkin(8) * t510;
t509 = (-t531 * t535 + t532 * t538) * qJD(1);
t508 = t509 ^ 2;
t439 = -t497 * pkin(3) - t508 * pkin(8) + t503 * t510 + t477;
t534 = sin(qJ(4));
t570 = cos(qJ(4));
t493 = t509 * t534 + t510 * t570;
t498 = t518 * t532 + t519 * t531;
t455 = qJD(4) * t493 - t497 * t570 + t498 * t534;
t492 = -t509 * t570 + t510 * t534;
t456 = -qJD(4) * t492 + t497 * t534 + t498 * t570;
t529 = qJD(2) + qJD(4);
t481 = mrSges(5,1) * t529 - mrSges(5,3) * t493;
t566 = t492 * t529;
t571 = -2 * qJD(5);
t541 = (-t456 + t566) * qJ(5) + t439 + (pkin(4) * t529 + t571) * t493;
t419 = t455 * pkin(4) + t541;
t483 = mrSges(6,1) * t493 + mrSges(6,2) * t529;
t552 = -g(1) * t539 - g(2) * t536;
t515 = -pkin(1) * t540 + qJDD(1) * pkin(7) + t552;
t565 = t515 * t535;
t569 = pkin(2) * t540;
t474 = qJDD(2) * pkin(2) - qJ(3) * t518 - t565 + (qJ(3) * t559 + t535 * t569 - g(3)) * t538;
t500 = -g(3) * t535 + t515 * t538;
t475 = qJ(3) * t519 - qJD(2) * t520 - t530 * t569 + t500;
t440 = -0.2e1 * qJD(3) * t510 + t474 * t532 - t475 * t531;
t424 = (qJD(2) * t509 - t498) * pkin(8) + (t509 * t510 + qJDD(2)) * pkin(3) + t440;
t441 = 0.2e1 * qJD(3) * t509 + t474 * t531 + t475 * t532;
t427 = -pkin(3) * t508 + pkin(8) * t497 - qJD(2) * t503 + t441;
t421 = t424 * t570 - t427 * t534;
t468 = pkin(4) * t492 - qJ(5) * t493;
t527 = t529 ^ 2;
t528 = qJDD(2) + qJDD(4);
t417 = -t528 * pkin(4) - t527 * qJ(5) + t468 * t493 + qJDD(5) - t421;
t411 = (t492 * t493 - t528) * pkin(9) + (t456 + t566) * pkin(5) + t417;
t484 = pkin(5) * t493 - pkin(9) * t529;
t488 = t492 ^ 2;
t414 = (pkin(4) + pkin(9)) * t455 + t541 - t493 * t484 - t488 * pkin(5);
t533 = sin(qJ(6));
t537 = cos(qJ(6));
t409 = t411 * t537 - t414 * t533;
t478 = t492 * t537 - t529 * t533;
t433 = qJD(6) * t478 + t455 * t533 + t528 * t537;
t454 = qJDD(6) + t456;
t479 = t492 * t533 + t529 * t537;
t457 = -mrSges(7,1) * t478 + mrSges(7,2) * t479;
t487 = qJD(6) + t493;
t458 = -mrSges(7,2) * t487 + mrSges(7,3) * t478;
t406 = m(7) * t409 + mrSges(7,1) * t454 - mrSges(7,3) * t433 - t457 * t479 + t458 * t487;
t410 = t411 * t533 + t414 * t537;
t432 = -qJD(6) * t479 + t455 * t537 - t528 * t533;
t459 = mrSges(7,1) * t487 - mrSges(7,3) * t479;
t407 = m(7) * t410 - mrSges(7,2) * t454 + mrSges(7,3) * t432 + t457 * t478 - t459 * t487;
t553 = -t533 * t406 + t407 * t537;
t550 = -m(6) * t419 + mrSges(6,3) * t456 + t483 * t493 - t553;
t482 = mrSges(6,1) * t492 - mrSges(6,3) * t529;
t562 = -mrSges(5,2) * t529 - mrSges(5,3) * t492 - t482;
t568 = mrSges(6,2) - mrSges(5,1);
t578 = m(5) * t439 + mrSges(5,2) * t456 - t568 * t455 + t481 * t493 + t562 * t492 - t550;
t577 = Ifges(5,1) + Ifges(6,2);
t576 = Ifges(5,5) - Ifges(6,4);
t574 = Ifges(5,3) + Ifges(6,1);
t573 = t492 * t579 + t493 * t567 + t529 * t575;
t501 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t509;
t502 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t510;
t391 = m(4) * t477 - t497 * mrSges(4,1) + t498 * mrSges(4,2) - t509 * t501 + t510 * t502 + t578;
t469 = mrSges(5,1) * t492 + mrSges(5,2) * t493;
t397 = t406 * t537 + t407 * t533;
t470 = -mrSges(6,2) * t492 - mrSges(6,3) * t493;
t548 = -m(6) * t417 - mrSges(6,1) * t456 - t470 * t493 - t397;
t393 = m(5) * t421 - mrSges(5,3) * t456 - t469 * t493 - t528 * t568 + t529 * t562 + t548;
t422 = t424 * t534 + t427 * t570;
t547 = -pkin(4) * t527 + qJ(5) * t528 - t468 * t492 + t422;
t415 = t529 * t571 - t547;
t413 = -pkin(5) * t455 - pkin(9) * t488 + ((2 * qJD(5)) + t484) * t529 + t547;
t549 = -m(7) * t413 + mrSges(7,1) * t432 - mrSges(7,2) * t433 + t458 * t478 - t459 * t479;
t544 = -m(6) * t415 + mrSges(6,3) * t528 + t483 * t529 - t549;
t403 = m(5) * t422 - mrSges(5,2) * t528 - t481 * t529 + (-t469 - t470) * t492 + (-mrSges(5,3) - mrSges(6,1)) * t455 + t544;
t390 = t393 * t570 + t403 * t534;
t495 = -mrSges(4,1) * t509 + mrSges(4,2) * t510;
t388 = m(4) * t440 + qJDD(2) * mrSges(4,1) - mrSges(4,3) * t498 + qJD(2) * t501 - t495 * t510 + t390;
t554 = -t393 * t534 + t403 * t570;
t389 = m(4) * t441 - qJDD(2) * mrSges(4,2) + mrSges(4,3) * t497 - qJD(2) * t502 + t495 * t509 + t554;
t383 = t388 * t532 + t389 * t531;
t564 = t492 * t575 - t493 * t576 - t529 * t574;
t563 = -t492 * t567 + t493 * t577 + t529 * t576;
t560 = qJD(1) * t538;
t555 = -t388 * t531 + t389 * t532;
t435 = Ifges(7,4) * t479 + Ifges(7,2) * t478 + Ifges(7,6) * t487;
t436 = Ifges(7,1) * t479 + Ifges(7,4) * t478 + Ifges(7,5) * t487;
t545 = mrSges(7,1) * t409 - mrSges(7,2) * t410 + Ifges(7,5) * t433 + Ifges(7,6) * t432 + Ifges(7,3) * t454 + t435 * t479 - t478 * t436;
t396 = mrSges(6,2) * t528 + t482 * t529 - t548;
t434 = Ifges(7,5) * t479 + Ifges(7,6) * t478 + Ifges(7,3) * t487;
t399 = -mrSges(7,1) * t413 + mrSges(7,3) * t410 + Ifges(7,4) * t433 + Ifges(7,2) * t432 + Ifges(7,6) * t454 - t434 * t479 + t436 * t487;
t400 = mrSges(7,2) * t413 - mrSges(7,3) * t409 + Ifges(7,1) * t433 + Ifges(7,4) * t432 + Ifges(7,5) * t454 + t434 * t478 - t435 * t487;
t542 = -mrSges(5,2) * t422 - mrSges(6,3) * t415 - pkin(4) * t396 - pkin(9) * t397 - t533 * t399 + t537 * t400 + t563 * t492 + qJ(5) * (-t470 * t492 + t544) + mrSges(6,2) * t417 + mrSges(5,1) * t421 + t574 * t528 + t573 * t493 + t576 * t456 + (-mrSges(6,1) * qJ(5) - t575) * t455;
t522 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t560;
t521 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t561;
t517 = (-mrSges(3,1) * t538 + mrSges(3,2) * t535) * qJD(1);
t514 = -pkin(7) * t540 + t551;
t513 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t535 + Ifges(3,4) * t538) * qJD(1);
t512 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t535 + Ifges(3,2) * t538) * qJD(1);
t499 = -g(3) * t538 - t565;
t491 = Ifges(4,1) * t510 + Ifges(4,4) * t509 + Ifges(4,5) * qJD(2);
t490 = Ifges(4,4) * t510 + Ifges(4,2) * t509 + Ifges(4,6) * qJD(2);
t489 = Ifges(4,5) * t510 + Ifges(4,6) * t509 + Ifges(4,3) * qJD(2);
t394 = -t455 * mrSges(6,2) - t492 * t482 - t550;
t384 = mrSges(6,1) * t417 + mrSges(5,2) * t439 - mrSges(5,3) * t421 - mrSges(6,3) * t419 + pkin(5) * t397 - qJ(5) * t394 - t455 * t567 + t456 * t577 + t492 * t564 + t528 * t576 - t529 * t573 + t545;
t382 = -mrSges(5,1) * t439 - mrSges(6,1) * t415 + mrSges(6,2) * t419 + mrSges(5,3) * t422 - pkin(4) * t394 - pkin(5) * t549 - pkin(9) * t553 - t537 * t399 - t533 * t400 + t455 * t579 + t456 * t567 + t493 * t564 + t528 * t575 + t529 * t563;
t381 = mrSges(4,2) * t477 - mrSges(4,3) * t440 + Ifges(4,1) * t498 + Ifges(4,4) * t497 + Ifges(4,5) * qJDD(2) - pkin(8) * t390 - qJD(2) * t490 - t382 * t534 + t384 * t570 + t489 * t509;
t380 = -mrSges(4,1) * t477 + mrSges(4,3) * t441 + Ifges(4,4) * t498 + Ifges(4,2) * t497 + Ifges(4,6) * qJDD(2) - pkin(3) * t578 + pkin(8) * t554 + qJD(2) * t491 + t382 * t570 + t534 * t384 - t510 * t489;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t557 - mrSges(2,2) * t552 + t535 * (mrSges(3,2) * t514 - mrSges(3,3) * t499 + Ifges(3,1) * t518 + Ifges(3,4) * t519 + Ifges(3,5) * qJDD(2) - qJ(3) * t383 - qJD(2) * t512 - t531 * t380 + t532 * t381) + t538 * (-mrSges(3,1) * t514 + mrSges(3,3) * t500 + Ifges(3,4) * t518 + Ifges(3,2) * t519 + Ifges(3,6) * qJDD(2) - pkin(2) * t391 + qJ(3) * t555 + qJD(2) * t513 + t532 * t380 + t531 * t381) + pkin(1) * ((-t521 * t535 + t522 * t538) * qJD(1) - m(3) * t514 - t518 * mrSges(3,2) + t519 * mrSges(3,1) - t391) + pkin(7) * (t538 * (m(3) * t500 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t519 - qJD(2) * t521 + t517 * t560 + t555) - t535 * (m(3) * t499 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t518 + qJD(2) * t522 - t517 * t561 + t383)); pkin(2) * t383 + t542 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t535 * t512 - t538 * t513) * qJD(1) + t510 * t490 + Ifges(3,5) * t518 + Ifges(3,6) * t519 - t509 * t491 + Ifges(4,6) * t497 + Ifges(4,5) * t498 + mrSges(3,1) * t499 - mrSges(3,2) * t500 - mrSges(4,2) * t441 + mrSges(4,1) * t440 + pkin(3) * t390; t391; t542; t396; t545;];
tauJ  = t1;
