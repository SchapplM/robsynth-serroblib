% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPRPRR11
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 20:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPRPRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR11_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:14:26
% EndTime: 2019-05-05 20:14:51
% DurationCPUTime: 24.66s
% Computational Cost: add. (353516->355), mult. (1121769->490), div. (0->0), fcn. (963914->16), ass. (0->165)
t589 = cos(pkin(6));
t592 = sin(qJ(3));
t588 = cos(pkin(7));
t633 = cos(qJ(3));
t618 = t588 * t633;
t584 = sin(pkin(7));
t619 = t584 * t633;
t585 = sin(pkin(6));
t587 = cos(pkin(12));
t625 = t585 * t587;
t583 = sin(pkin(12));
t627 = t583 * t585;
t635 = -t589 * t619 + t592 * t627 - t618 * t625;
t634 = 0.2e1 * t587;
t632 = pkin(9) * t583;
t631 = pkin(9) * t584;
t630 = pkin(9) * qJDD(1);
t629 = qJ(2) * t585;
t597 = qJD(1) ^ 2;
t593 = sin(qJ(1));
t596 = cos(qJ(1));
t616 = t593 * g(1) - g(2) * t596;
t570 = qJDD(1) * pkin(1) + t597 * t629 + t616;
t628 = t570 * t589;
t626 = t584 * t592;
t624 = t585 * t588;
t623 = t588 * t592;
t567 = (t584 * t589 + t587 * t624) * qJD(1) * pkin(9);
t611 = -g(1) * t596 - g(2) * t593;
t571 = -pkin(1) * t597 + qJDD(1) * t629 + t611;
t607 = -pkin(2) * t587 - t583 * t631;
t622 = qJD(1) * t585;
t603 = qJD(1) * t607 * t622 + t588 * t630;
t617 = qJD(2) * t622;
t608 = -g(3) * t625 - 0.2e1 * t583 * t617 + t587 * t628;
t520 = (pkin(2) * qJDD(1) + qJD(1) * t567) * t589 + (-t603 * t585 - t571) * t583 + t608;
t572 = (pkin(2) * t589 - t624 * t632) * qJD(1);
t620 = t587 * t571 + t583 * t628 + t617 * t634;
t521 = (-qJD(1) * t572 + t584 * t630) * t589 + (-g(3) * t583 + t603 * t587) * t585 + t620;
t615 = -g(3) * t589 + qJDD(2);
t529 = (-t570 + t607 * qJDD(1) + (-t567 * t587 + t572 * t583) * qJD(1)) * t585 + t615;
t493 = t520 * t623 + t633 * t521 + t529 * t626;
t555 = t635 * qJD(1);
t600 = t589 * t626 + (t633 * t583 + t587 * t623) * t585;
t556 = t600 * qJD(1);
t540 = pkin(3) * t555 - qJ(4) * t556;
t604 = -t584 * t625 + t588 * t589;
t568 = t604 * qJD(1) + qJD(3);
t564 = t568 ^ 2;
t565 = t604 * qJDD(1) + qJDD(3);
t484 = -pkin(3) * t564 + qJ(4) * t565 - t540 * t555 + t493;
t507 = -t520 * t584 + t588 * t529;
t542 = qJD(3) * t556 + t635 * qJDD(1);
t543 = -t555 * qJD(3) + t600 * qJDD(1);
t487 = (t555 * t568 - t543) * qJ(4) + (t556 * t568 + t542) * pkin(3) + t507;
t582 = sin(pkin(13));
t586 = cos(pkin(13));
t550 = t556 * t586 + t568 * t582;
t476 = -0.2e1 * qJD(4) * t550 - t484 * t582 + t586 * t487;
t535 = t543 * t586 + t565 * t582;
t549 = -t556 * t582 + t568 * t586;
t473 = (t549 * t555 - t535) * pkin(10) + (t549 * t550 + t542) * pkin(4) + t476;
t477 = 0.2e1 * qJD(4) * t549 + t586 * t484 + t582 * t487;
t533 = pkin(4) * t555 - pkin(10) * t550;
t534 = -t543 * t582 + t565 * t586;
t546 = t549 ^ 2;
t475 = -pkin(4) * t546 + pkin(10) * t534 - t533 * t555 + t477;
t591 = sin(qJ(5));
t595 = cos(qJ(5));
t470 = t591 * t473 + t595 * t475;
t523 = t549 * t595 - t550 * t591;
t524 = t549 * t591 + t550 * t595;
t506 = -pkin(5) * t523 - pkin(11) * t524;
t539 = qJDD(5) + t542;
t554 = qJD(5) + t555;
t553 = t554 ^ 2;
t468 = -pkin(5) * t553 + pkin(11) * t539 + t506 * t523 + t470;
t492 = t520 * t618 - t592 * t521 + t529 * t619;
t483 = -t565 * pkin(3) - t564 * qJ(4) + t556 * t540 + qJDD(4) - t492;
t478 = -t534 * pkin(4) - t546 * pkin(10) + t550 * t533 + t483;
t497 = -qJD(5) * t524 + t534 * t595 - t535 * t591;
t498 = qJD(5) * t523 + t534 * t591 + t535 * t595;
t471 = (-t523 * t554 - t498) * pkin(11) + (t524 * t554 - t497) * pkin(5) + t478;
t590 = sin(qJ(6));
t594 = cos(qJ(6));
t465 = -t468 * t590 + t471 * t594;
t508 = -t524 * t590 + t554 * t594;
t481 = qJD(6) * t508 + t498 * t594 + t539 * t590;
t509 = t524 * t594 + t554 * t590;
t494 = -mrSges(7,1) * t508 + mrSges(7,2) * t509;
t496 = qJDD(6) - t497;
t522 = qJD(6) - t523;
t499 = -mrSges(7,2) * t522 + mrSges(7,3) * t508;
t462 = m(7) * t465 + mrSges(7,1) * t496 - mrSges(7,3) * t481 - t494 * t509 + t499 * t522;
t466 = t468 * t594 + t471 * t590;
t480 = -qJD(6) * t509 - t498 * t590 + t539 * t594;
t500 = mrSges(7,1) * t522 - mrSges(7,3) * t509;
t463 = m(7) * t466 - mrSges(7,2) * t496 + mrSges(7,3) * t480 + t494 * t508 - t500 * t522;
t454 = -t462 * t590 + t594 * t463;
t505 = -mrSges(6,1) * t523 + mrSges(6,2) * t524;
t511 = mrSges(6,1) * t554 - mrSges(6,3) * t524;
t448 = m(6) * t470 - mrSges(6,2) * t539 + mrSges(6,3) * t497 + t505 * t523 - t511 * t554 + t454;
t469 = t473 * t595 - t475 * t591;
t467 = -pkin(5) * t539 - pkin(11) * t553 + t506 * t524 - t469;
t464 = -m(7) * t467 + t480 * mrSges(7,1) - mrSges(7,2) * t481 + t508 * t499 - t500 * t509;
t510 = -mrSges(6,2) * t554 + mrSges(6,3) * t523;
t458 = m(6) * t469 + mrSges(6,1) * t539 - mrSges(6,3) * t498 - t505 * t524 + t510 * t554 + t464;
t445 = t591 * t448 + t595 * t458;
t525 = -mrSges(5,1) * t549 + mrSges(5,2) * t550;
t531 = -mrSges(5,2) * t555 + mrSges(5,3) * t549;
t443 = m(5) * t476 + mrSges(5,1) * t542 - mrSges(5,3) * t535 - t525 * t550 + t531 * t555 + t445;
t532 = mrSges(5,1) * t555 - mrSges(5,3) * t550;
t613 = t595 * t448 - t458 * t591;
t444 = m(5) * t477 - mrSges(5,2) * t542 + mrSges(5,3) * t534 + t525 * t549 - t532 * t555 + t613;
t437 = t586 * t443 + t582 * t444;
t453 = t594 * t462 + t590 * t463;
t541 = mrSges(4,1) * t555 + mrSges(4,2) * t556;
t552 = mrSges(4,1) * t568 - mrSges(4,3) * t556;
t614 = -t443 * t582 + t586 * t444;
t434 = m(4) * t493 - mrSges(4,2) * t565 - mrSges(4,3) * t542 - t541 * t555 - t552 * t568 + t614;
t551 = -mrSges(4,2) * t568 - mrSges(4,3) * t555;
t436 = m(4) * t507 + mrSges(4,1) * t542 + mrSges(4,2) * t543 + t551 * t555 + t552 * t556 + t437;
t601 = m(6) * t478 - t497 * mrSges(6,1) + mrSges(6,2) * t498 - t523 * t510 + t511 * t524 + t453;
t452 = m(5) * t483 - t534 * mrSges(5,1) + mrSges(5,2) * t535 - t549 * t531 + t532 * t550 + t601;
t451 = m(4) * t492 + mrSges(4,1) * t565 - mrSges(4,3) * t543 - t541 * t556 + t551 * t568 - t452;
t425 = t434 * t626 + t588 * t436 + t451 * t619;
t430 = t633 * t434 - t451 * t592;
t426 = t434 * t623 - t436 * t584 + t451 * t618;
t610 = -mrSges(3,1) * t587 + mrSges(3,2) * t583;
t606 = mrSges(3,1) * t589 - mrSges(3,3) * t627;
t605 = -mrSges(3,2) * t589 + mrSges(3,3) * t625;
t488 = Ifges(7,5) * t509 + Ifges(7,6) * t508 + Ifges(7,3) * t522;
t490 = Ifges(7,1) * t509 + Ifges(7,4) * t508 + Ifges(7,5) * t522;
t455 = -mrSges(7,1) * t467 + mrSges(7,3) * t466 + Ifges(7,4) * t481 + Ifges(7,2) * t480 + Ifges(7,6) * t496 - t488 * t509 + t490 * t522;
t489 = Ifges(7,4) * t509 + Ifges(7,2) * t508 + Ifges(7,6) * t522;
t456 = mrSges(7,2) * t467 - mrSges(7,3) * t465 + Ifges(7,1) * t481 + Ifges(7,4) * t480 + Ifges(7,5) * t496 + t488 * t508 - t489 * t522;
t501 = Ifges(6,5) * t524 + Ifges(6,6) * t523 + Ifges(6,3) * t554;
t502 = Ifges(6,4) * t524 + Ifges(6,2) * t523 + Ifges(6,6) * t554;
t438 = mrSges(6,2) * t478 - mrSges(6,3) * t469 + Ifges(6,1) * t498 + Ifges(6,4) * t497 + Ifges(6,5) * t539 - pkin(11) * t453 - t455 * t590 + t456 * t594 + t501 * t523 - t502 * t554;
t503 = Ifges(6,1) * t524 + Ifges(6,4) * t523 + Ifges(6,5) * t554;
t599 = mrSges(7,1) * t465 - mrSges(7,2) * t466 + Ifges(7,5) * t481 + Ifges(7,6) * t480 + Ifges(7,3) * t496 + t489 * t509 - t490 * t508;
t439 = -mrSges(6,1) * t478 + mrSges(6,3) * t470 + Ifges(6,4) * t498 + Ifges(6,2) * t497 + Ifges(6,6) * t539 - pkin(5) * t453 - t501 * t524 + t503 * t554 - t599;
t512 = Ifges(5,5) * t550 + Ifges(5,6) * t549 + Ifges(5,3) * t555;
t514 = Ifges(5,1) * t550 + Ifges(5,4) * t549 + Ifges(5,5) * t555;
t427 = -mrSges(5,1) * t483 + mrSges(5,3) * t477 + Ifges(5,4) * t535 + Ifges(5,2) * t534 + Ifges(5,6) * t542 - pkin(4) * t601 + pkin(10) * t613 + t591 * t438 + t595 * t439 - t550 * t512 + t555 * t514;
t513 = Ifges(5,4) * t550 + Ifges(5,2) * t549 + Ifges(5,6) * t555;
t428 = mrSges(5,2) * t483 - mrSges(5,3) * t476 + Ifges(5,1) * t535 + Ifges(5,4) * t534 + Ifges(5,5) * t542 - pkin(10) * t445 + t438 * t595 - t439 * t591 + t512 * t549 - t513 * t555;
t536 = Ifges(4,5) * t556 - Ifges(4,6) * t555 + Ifges(4,3) * t568;
t537 = Ifges(4,4) * t556 - Ifges(4,2) * t555 + Ifges(4,6) * t568;
t421 = mrSges(4,2) * t507 - mrSges(4,3) * t492 + Ifges(4,1) * t543 - Ifges(4,4) * t542 + Ifges(4,5) * t565 - qJ(4) * t437 - t427 * t582 + t428 * t586 - t536 * t555 - t537 * t568;
t538 = Ifges(4,1) * t556 - Ifges(4,4) * t555 + Ifges(4,5) * t568;
t598 = mrSges(6,1) * t469 - mrSges(6,2) * t470 + Ifges(6,5) * t498 + Ifges(6,6) * t497 + Ifges(6,3) * t539 + pkin(5) * t464 + pkin(11) * t454 + t594 * t455 + t590 * t456 + t524 * t502 - t523 * t503;
t422 = -t598 + (-Ifges(5,3) - Ifges(4,2)) * t542 - pkin(4) * t445 - pkin(3) * t437 + Ifges(4,6) * t565 + t568 * t538 - t556 * t536 + t549 * t514 - t550 * t513 + Ifges(4,4) * t543 - Ifges(5,6) * t534 - Ifges(5,5) * t535 - mrSges(4,1) * t507 + mrSges(4,3) * t493 + mrSges(5,2) * t477 - mrSges(5,1) * t476;
t602 = pkin(9) * t430 + t592 * t421 + t633 * t422;
t574 = t605 * qJD(1);
t573 = t606 * qJD(1);
t569 = t610 * t622;
t557 = -t570 * t585 + t615;
t548 = -g(3) * t627 + t620;
t547 = -t571 * t583 + t608;
t429 = m(3) * t548 + t605 * qJDD(1) + (t569 * t625 - t573 * t589) * qJD(1) + t430;
t424 = m(3) * t557 + (t610 * qJDD(1) + (t573 * t583 - t574 * t587) * qJD(1)) * t585 + t425;
t423 = m(3) * t547 + t606 * qJDD(1) + (-t569 * t627 + t574 * t589) * qJD(1) + t426;
t420 = mrSges(4,1) * t492 - mrSges(4,2) * t493 + Ifges(4,5) * t543 - Ifges(4,6) * t542 + Ifges(4,3) * t565 - pkin(3) * t452 + qJ(4) * t614 + t586 * t427 + t582 * t428 + t556 * t537 + t555 * t538;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t616 - mrSges(2,2) * t611 + (mrSges(3,1) * t547 - mrSges(3,2) * t548 + pkin(2) * t426 + t588 * t420 + pkin(1) * (t423 * t587 + t429 * t583) + Ifges(3,3) * t589 * qJDD(1) + t602 * t584) * t589 + (t583 * (mrSges(3,2) * t557 - mrSges(3,3) * t547 + t633 * t421 - t592 * t422 - t425 * t631) + t587 * (-mrSges(3,1) * t557 + mrSges(3,3) * t548 - pkin(2) * t425 - t584 * t420) - pkin(1) * t424 + qJ(2) * (-t423 * t583 + t429 * t587) + (-t426 * t632 + t587 * t602) * t588 + ((Ifges(3,2) * t587 ^ 2 + (Ifges(3,1) * t583 + Ifges(3,4) * t634) * t583) * t585 + 0.2e1 * t589 * (Ifges(3,5) * t583 + Ifges(3,6) * t587)) * qJDD(1)) * t585; t424; t420; t452; t598; t599;];
tauJ  = t1;
