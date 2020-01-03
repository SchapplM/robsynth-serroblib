% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:29
% EndTime: 2019-12-31 20:03:31
% DurationCPUTime: 1.32s
% Computational Cost: add. (4140->238), mult. (8731->277), div. (0->0), fcn. (4758->6), ass. (0->94)
t474 = Ifges(3,1) + Ifges(4,1);
t473 = Ifges(5,1) + Ifges(6,1);
t465 = Ifges(3,4) - Ifges(4,5);
t464 = Ifges(5,4) + Ifges(6,4);
t463 = Ifges(3,5) + Ifges(4,4);
t462 = Ifges(5,5) + Ifges(6,5);
t472 = Ifges(3,2) + Ifges(4,3);
t471 = Ifges(5,2) + Ifges(6,2);
t461 = Ifges(3,6) - Ifges(4,6);
t460 = Ifges(5,6) + Ifges(6,6);
t470 = Ifges(3,3) + Ifges(4,2);
t469 = Ifges(5,3) + Ifges(6,3);
t428 = sin(qJ(4));
t429 = sin(qJ(2));
t431 = cos(qJ(4));
t432 = cos(qJ(2));
t393 = (-t429 * t428 - t432 * t431) * qJD(1);
t448 = qJD(1) * qJD(2);
t445 = t432 * t448;
t403 = qJDD(1) * t429 + t445;
t444 = t429 * t448;
t404 = qJDD(1) * t432 - t444;
t361 = qJD(4) * t393 + t403 * t431 - t404 * t428;
t394 = (-t432 * t428 + t429 * t431) * qJD(1);
t371 = -mrSges(6,1) * t393 + mrSges(6,2) * t394;
t435 = qJD(1) ^ 2;
t430 = sin(qJ(1));
t433 = cos(qJ(1));
t442 = -g(1) * t433 - g(2) * t430;
t396 = -pkin(1) * t435 + qJDD(1) * pkin(6) + t442;
t375 = -g(3) * t429 + t432 * t396;
t400 = (-t432 * pkin(2) - t429 * qJ(3)) * qJD(1);
t434 = qJD(2) ^ 2;
t449 = t432 * qJD(1);
t467 = 2 * qJD(3);
t354 = -pkin(2) * t434 + qJDD(2) * qJ(3) + qJD(2) * t467 + t400 * t449 + t375;
t450 = t429 * qJD(1);
t411 = -qJD(2) * pkin(3) - pkin(7) * t450;
t458 = t432 ^ 2 * t435;
t349 = -pkin(3) * t458 - pkin(7) * t404 + qJD(2) * t411 + t354;
t374 = -t432 * g(3) - t429 * t396;
t359 = -qJDD(2) * pkin(2) - qJ(3) * t434 + t400 * t450 + qJDD(3) - t374;
t350 = (-t403 + t445) * pkin(7) + (-t429 * t432 * t435 - qJDD(2)) * pkin(3) + t359;
t343 = -t349 * t428 + t431 * t350;
t420 = -qJDD(2) + qJDD(4);
t421 = -qJD(2) + qJD(4);
t338 = -0.2e1 * qJD(5) * t394 + (t393 * t421 - t361) * qJ(5) + (t393 * t394 + t420) * pkin(4) + t343;
t376 = -mrSges(6,2) * t421 + mrSges(6,3) * t393;
t447 = m(6) * t338 + t420 * mrSges(6,1) + t421 * t376;
t335 = -mrSges(6,3) * t361 - t371 * t394 + t447;
t344 = t431 * t349 + t428 * t350;
t360 = -qJD(4) * t394 - t403 * t428 - t404 * t431;
t378 = pkin(4) * t421 - qJ(5) * t394;
t386 = t393 ^ 2;
t340 = -pkin(4) * t386 + qJ(5) * t360 + 0.2e1 * qJD(5) * t393 - t378 * t421 + t344;
t455 = t393 * t464 + t394 * t473 + t421 * t462;
t456 = -t393 * t471 - t394 * t464 - t421 * t460;
t468 = mrSges(5,1) * t343 + mrSges(6,1) * t338 - mrSges(5,2) * t344 - mrSges(6,2) * t340 + pkin(4) * t335 + t460 * t360 + t462 * t361 - t455 * t393 - t456 * t394 + t420 * t469;
t466 = mrSges(3,3) + mrSges(4,2);
t457 = -t393 * t460 - t394 * t462 - t421 * t469;
t454 = t470 * qJD(2) + (t429 * t463 + t432 * t461) * qJD(1);
t453 = -t461 * qJD(2) + (-t429 * t465 - t472 * t432) * qJD(1);
t452 = t463 * qJD(2) + (t429 * t474 + t432 * t465) * qJD(1);
t451 = t430 * g(1) - t433 * g(2);
t446 = m(6) * t340 + t360 * mrSges(6,3) + t393 * t371;
t372 = -mrSges(5,1) * t393 + mrSges(5,2) * t394;
t377 = -mrSges(5,2) * t421 + mrSges(5,3) * t393;
t332 = m(5) * t343 + mrSges(5,1) * t420 + t377 * t421 + (-t371 - t372) * t394 + (-mrSges(5,3) - mrSges(6,3)) * t361 + t447;
t379 = mrSges(6,1) * t421 - mrSges(6,3) * t394;
t380 = mrSges(5,1) * t421 - mrSges(5,3) * t394;
t334 = m(5) * t344 + mrSges(5,3) * t360 + t372 * t393 + (-t379 - t380) * t421 + (-mrSges(5,2) - mrSges(6,2)) * t420 + t446;
t443 = -t428 * t332 + t431 * t334;
t395 = -qJDD(1) * pkin(1) - t435 * pkin(6) - t451;
t441 = -t404 * pkin(2) + t395 + (-t403 - t445) * qJ(3);
t345 = -pkin(2) * t444 + pkin(3) * t404 - pkin(7) * t458 - t441 + (t411 + t467) * t450;
t342 = -pkin(4) * t360 - qJ(5) * t386 + t378 * t394 + qJDD(5) + t345;
t336 = m(6) * t342 - t360 * mrSges(6,1) + t361 * mrSges(6,2) - t393 * t376 + t394 * t379;
t330 = t431 * t332 + t428 * t334;
t401 = (-t432 * mrSges(4,1) - t429 * mrSges(4,3)) * qJD(1);
t408 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t450;
t440 = m(4) * t354 + qJDD(2) * mrSges(4,3) + qJD(2) * t408 + t401 * t449 + t443;
t410 = mrSges(4,2) * t449 + qJD(2) * mrSges(4,3);
t439 = m(4) * t359 - qJDD(2) * mrSges(4,1) - qJD(2) * t410 + t330;
t437 = m(5) * t345 - t360 * mrSges(5,1) + mrSges(5,2) * t361 - t393 * t377 + t380 * t394 + t336;
t351 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t450 + t441;
t436 = m(4) * t351 - t437;
t409 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t449;
t407 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t450;
t402 = (-t432 * mrSges(3,1) + t429 * mrSges(3,2)) * qJD(1);
t331 = (-t408 * t429 - t410 * t432) * qJD(1) - mrSges(4,1) * t404 - mrSges(4,3) * t403 + t436;
t329 = t403 * mrSges(4,2) + t401 * t450 + t439;
t328 = mrSges(5,2) * t345 + mrSges(6,2) * t342 - mrSges(5,3) * t343 - mrSges(6,3) * t338 - qJ(5) * t335 + t464 * t360 + t361 * t473 - t457 * t393 + t462 * t420 + t456 * t421;
t327 = -mrSges(5,1) * t345 + mrSges(5,3) * t344 - mrSges(6,1) * t342 + mrSges(6,3) * t340 - pkin(4) * t336 + qJ(5) * t446 + (-qJ(5) * t379 + t455) * t421 + (-qJ(5) * mrSges(6,2) + t460) * t420 + t457 * t394 + t464 * t361 + t471 * t360;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t451 - mrSges(2,2) * t442 + t429 * (mrSges(3,2) * t395 + mrSges(4,2) * t359 - mrSges(3,3) * t374 - mrSges(4,3) * t351 - pkin(7) * t330 - qJ(3) * t331 + t453 * qJD(2) + t463 * qJDD(2) - t428 * t327 + t431 * t328 + t474 * t403 + t465 * t404 + t454 * t449) + t432 * (-mrSges(3,1) * t395 - mrSges(4,1) * t351 + mrSges(4,2) * t354 + mrSges(3,3) * t375 - pkin(2) * t331 + pkin(3) * t437 - pkin(7) * t443 + t452 * qJD(2) + t461 * qJDD(2) - t431 * t327 - t428 * t328 + t465 * t403 + t472 * t404 - t454 * t450) + pkin(1) * ((mrSges(3,1) + mrSges(4,1)) * t404 + (-mrSges(3,2) + mrSges(4,3)) * t403 + ((t409 + t410) * t432 + (-t407 + t408) * t429) * qJD(1) - m(3) * t395 - t436) + pkin(6) * (t432 * (m(3) * t375 - qJDD(2) * mrSges(3,2) - qJD(2) * t407 + t402 * t449 + t466 * t404 + t440) + (-m(3) * t374 - qJDD(2) * mrSges(3,1) - qJD(2) * t409 + t466 * t403 + (t401 + t402) * t450 + t439) * t429); mrSges(3,1) * t374 - mrSges(3,2) * t375 + mrSges(4,3) * t354 - mrSges(4,1) * t359 + t463 * t403 + qJ(3) * t440 + (qJ(3) * mrSges(4,2) + t461) * t404 + (-t453 * t429 - t452 * t432) * qJD(1) + t470 * qJDD(2) - pkin(2) * t329 - pkin(3) * t330 - t468; t329; t468; t336;];
tauJ = t1;
