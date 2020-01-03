% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR12
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR12_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:20
% EndTime: 2019-12-31 20:29:24
% DurationCPUTime: 2.93s
% Computational Cost: add. (26274->290), mult. (54017->356), div. (0->0), fcn. (31749->8), ass. (0->114)
t633 = Ifges(3,1) + Ifges(4,1);
t628 = Ifges(3,4) - Ifges(4,5);
t627 = Ifges(3,5) + Ifges(4,4);
t632 = Ifges(3,2) + Ifges(4,3);
t626 = Ifges(3,6) - Ifges(4,6);
t631 = Ifges(3,3) + Ifges(4,2);
t596 = sin(qJ(4));
t597 = sin(qJ(2));
t600 = cos(qJ(4));
t601 = cos(qJ(2));
t553 = (t596 * t597 + t600 * t601) * qJD(1);
t630 = 2 * qJD(3);
t629 = mrSges(3,3) + mrSges(4,2);
t604 = qJD(1) ^ 2;
t625 = t601 ^ 2 * t604;
t598 = sin(qJ(1));
t602 = cos(qJ(1));
t578 = -t602 * g(1) - t598 * g(2);
t556 = -t604 * pkin(1) + qJDD(1) * pkin(6) + t578;
t538 = -t597 * g(3) + t601 * t556;
t567 = (-mrSges(3,1) * t601 + mrSges(3,2) * t597) * qJD(1);
t618 = qJD(1) * qJD(2);
t617 = t597 * t618;
t569 = t601 * qJDD(1) - t617;
t620 = qJD(1) * t597;
t572 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t620;
t565 = (-pkin(2) * t601 - qJ(3) * t597) * qJD(1);
t603 = qJD(2) ^ 2;
t619 = qJD(1) * t601;
t519 = -t603 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t630 + t565 * t619 + t538;
t566 = (-mrSges(4,1) * t601 - mrSges(4,3) * t597) * qJD(1);
t573 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t620;
t576 = -qJD(2) * pkin(3) - pkin(7) * t620;
t510 = -pkin(3) * t625 - t569 * pkin(7) + qJD(2) * t576 + t519;
t537 = -t601 * g(3) - t597 * t556;
t523 = -qJDD(2) * pkin(2) - t603 * qJ(3) + t565 * t620 + qJDD(3) - t537;
t616 = t601 * t618;
t568 = t597 * qJDD(1) + t616;
t511 = (-t568 + t616) * pkin(7) + (-t597 * t601 * t604 - qJDD(2)) * pkin(3) + t523;
t502 = t600 * t510 + t596 * t511;
t554 = (-t596 * t601 + t597 * t600) * qJD(1);
t524 = -t554 * qJD(4) - t596 * t568 - t600 * t569;
t533 = t553 * mrSges(5,1) + t554 * mrSges(5,2);
t588 = -qJD(2) + qJD(4);
t540 = t588 * mrSges(5,1) - t554 * mrSges(5,3);
t587 = -qJDD(2) + qJDD(4);
t577 = t598 * g(1) - t602 * g(2);
t555 = -qJDD(1) * pkin(1) - t604 * pkin(6) - t577;
t609 = -t569 * pkin(2) + t555 + (-t568 - t616) * qJ(3);
t504 = -pkin(2) * t617 + t569 * pkin(3) - pkin(7) * t625 - t609 + (t576 + t630) * t620;
t525 = -t553 * qJD(4) + t600 * t568 - t596 * t569;
t498 = t504 + (t553 * t588 - t525) * pkin(8) + (t554 * t588 - t524) * pkin(4);
t534 = t553 * pkin(4) - t554 * pkin(8);
t586 = t588 ^ 2;
t500 = -t586 * pkin(4) + t587 * pkin(8) - t553 * t534 + t502;
t595 = sin(qJ(5));
t599 = cos(qJ(5));
t496 = t599 * t498 - t595 * t500;
t535 = -t595 * t554 + t599 * t588;
t507 = t535 * qJD(5) + t599 * t525 + t595 * t587;
t536 = t599 * t554 + t595 * t588;
t517 = -t535 * mrSges(6,1) + t536 * mrSges(6,2);
t522 = qJDD(5) - t524;
t546 = qJD(5) + t553;
t526 = -t546 * mrSges(6,2) + t535 * mrSges(6,3);
t494 = m(6) * t496 + t522 * mrSges(6,1) - t507 * mrSges(6,3) - t536 * t517 + t546 * t526;
t497 = t595 * t498 + t599 * t500;
t506 = -t536 * qJD(5) - t595 * t525 + t599 * t587;
t527 = t546 * mrSges(6,1) - t536 * mrSges(6,3);
t495 = m(6) * t497 - t522 * mrSges(6,2) + t506 * mrSges(6,3) + t535 * t517 - t546 * t527;
t612 = -t595 * t494 + t599 * t495;
t486 = m(5) * t502 - t587 * mrSges(5,2) + t524 * mrSges(5,3) - t553 * t533 - t588 * t540 + t612;
t501 = -t596 * t510 + t600 * t511;
t539 = -t588 * mrSges(5,2) - t553 * mrSges(5,3);
t499 = -t587 * pkin(4) - t586 * pkin(8) + t554 * t534 - t501;
t607 = -m(6) * t499 + t506 * mrSges(6,1) - t507 * mrSges(6,2) + t535 * t526 - t536 * t527;
t490 = m(5) * t501 + t587 * mrSges(5,1) - t525 * mrSges(5,3) - t554 * t533 + t588 * t539 + t607;
t613 = t600 * t486 - t596 * t490;
t608 = m(4) * t519 + qJDD(2) * mrSges(4,3) + qJD(2) * t573 + t566 * t619 + t613;
t479 = m(3) * t538 - qJDD(2) * mrSges(3,2) - qJD(2) * t572 + t567 * t619 + t629 * t569 + t608;
t574 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t619;
t481 = t596 * t486 + t600 * t490;
t575 = mrSges(4,2) * t619 + qJD(2) * mrSges(4,3);
t606 = -m(4) * t523 + qJDD(2) * mrSges(4,1) + qJD(2) * t575 - t481;
t480 = m(3) * t537 + qJDD(2) * mrSges(3,1) + qJD(2) * t574 - t629 * t568 + (-t566 - t567) * t620 + t606;
t614 = t601 * t479 - t597 * t480;
t472 = m(2) * t578 - t604 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t614;
t516 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t620 + t609;
t487 = t599 * t494 + t595 * t495;
t610 = -m(5) * t504 + t524 * mrSges(5,1) - t525 * mrSges(5,2) - t553 * t539 - t554 * t540 - t487;
t484 = m(4) * t516 - t569 * mrSges(4,1) - t568 * mrSges(4,3) - t573 * t620 - t575 * t619 + t610;
t605 = -m(3) * t555 + t569 * mrSges(3,1) - t568 * mrSges(3,2) - t572 * t620 + t574 * t619 - t484;
t483 = m(2) * t577 + qJDD(1) * mrSges(2,1) - t604 * mrSges(2,2) + t605;
t624 = t598 * t472 + t602 * t483;
t473 = t597 * t479 + t601 * t480;
t623 = t631 * qJD(2) + (t627 * t597 + t626 * t601) * qJD(1);
t622 = -t626 * qJD(2) + (-t628 * t597 - t632 * t601) * qJD(1);
t621 = t627 * qJD(2) + (t633 * t597 + t628 * t601) * qJD(1);
t615 = t602 * t472 - t598 * t483;
t530 = Ifges(5,1) * t554 - Ifges(5,4) * t553 + Ifges(5,5) * t588;
t529 = Ifges(5,4) * t554 - Ifges(5,2) * t553 + Ifges(5,6) * t588;
t528 = Ifges(5,5) * t554 - Ifges(5,6) * t553 + Ifges(5,3) * t588;
t514 = Ifges(6,1) * t536 + Ifges(6,4) * t535 + Ifges(6,5) * t546;
t513 = Ifges(6,4) * t536 + Ifges(6,2) * t535 + Ifges(6,6) * t546;
t512 = Ifges(6,5) * t536 + Ifges(6,6) * t535 + Ifges(6,3) * t546;
t489 = mrSges(6,2) * t499 - mrSges(6,3) * t496 + Ifges(6,1) * t507 + Ifges(6,4) * t506 + Ifges(6,5) * t522 + t535 * t512 - t546 * t513;
t488 = -mrSges(6,1) * t499 + mrSges(6,3) * t497 + Ifges(6,4) * t507 + Ifges(6,2) * t506 + Ifges(6,6) * t522 - t536 * t512 + t546 * t514;
t475 = -mrSges(5,1) * t504 - mrSges(6,1) * t496 + mrSges(6,2) * t497 + mrSges(5,3) * t502 + Ifges(5,4) * t525 - Ifges(6,5) * t507 + Ifges(5,2) * t524 + Ifges(5,6) * t587 - Ifges(6,6) * t506 - Ifges(6,3) * t522 - pkin(4) * t487 - t536 * t513 + t535 * t514 - t554 * t528 + t588 * t530;
t474 = mrSges(5,2) * t504 - mrSges(5,3) * t501 + Ifges(5,1) * t525 + Ifges(5,4) * t524 + Ifges(5,5) * t587 - pkin(8) * t487 - t595 * t488 + t599 * t489 - t553 * t528 - t588 * t529;
t469 = mrSges(3,2) * t555 + mrSges(4,2) * t523 - mrSges(3,3) * t537 - mrSges(4,3) * t516 - pkin(7) * t481 - qJ(3) * t484 + t622 * qJD(2) + t627 * qJDD(2) + t600 * t474 - t596 * t475 + t633 * t568 + t628 * t569 + t623 * t619;
t468 = -mrSges(3,1) * t555 - mrSges(4,1) * t516 + mrSges(4,2) * t519 + mrSges(3,3) * t538 - pkin(2) * t484 - pkin(3) * t610 - pkin(7) * t613 + t621 * qJD(2) + t626 * qJDD(2) - t596 * t474 - t600 * t475 + t628 * t568 + t632 * t569 - t623 * t620;
t467 = t554 * t529 - mrSges(3,1) * t537 + mrSges(3,2) * t538 + t553 * t530 + t595 * t489 + t599 * t488 + Ifges(5,3) * t587 + mrSges(2,3) * t578 - t631 * qJDD(2) + (-qJ(3) * mrSges(4,2) - t626) * t569 + (pkin(2) * mrSges(4,2) - t627) * t568 + pkin(3) * t481 - qJ(3) * t608 - pkin(2) * t606 + pkin(4) * t607 - mrSges(4,3) * t519 + mrSges(4,1) * t523 + Ifges(5,6) * t524 + Ifges(5,5) * t525 + mrSges(5,1) * t501 - mrSges(5,2) * t502 - pkin(1) * t473 + Ifges(2,6) * qJDD(1) + mrSges(2,1) * g(3) + pkin(8) * t612 + (t621 * t601 + (pkin(2) * t566 + t622) * t597) * qJD(1) + t604 * Ifges(2,5);
t466 = -mrSges(2,2) * g(3) - mrSges(2,3) * t577 + Ifges(2,5) * qJDD(1) - t604 * Ifges(2,6) - pkin(6) * t473 - t597 * t468 + t601 * t469;
t1 = [-m(1) * g(1) + t615; -m(1) * g(2) + t624; (-m(1) - m(2)) * g(3) + t473; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t624 + t602 * t466 - t598 * t467; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t615 + t598 * t466 + t602 * t467; -mrSges(1,1) * g(2) + mrSges(2,1) * t577 + mrSges(1,2) * g(1) - mrSges(2,2) * t578 + Ifges(2,3) * qJDD(1) + pkin(1) * t605 + pkin(6) * t614 + t601 * t468 + t597 * t469;];
tauB = t1;
