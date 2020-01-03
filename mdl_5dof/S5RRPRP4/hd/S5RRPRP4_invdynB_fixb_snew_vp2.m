% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRP4
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRP4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:42
% EndTime: 2019-12-31 19:52:44
% DurationCPUTime: 1.03s
% Computational Cost: add. (8932->200), mult. (10601->235), div. (0->0), fcn. (4472->6), ass. (0->82)
t443 = Ifges(5,1) + Ifges(6,1);
t434 = Ifges(5,4) - Ifges(6,5);
t432 = (Ifges(5,5) + Ifges(6,4));
t442 = Ifges(5,2) + Ifges(6,3);
t430 = (Ifges(5,6) - Ifges(6,6));
t441 = (Ifges(5,3) + Ifges(6,2));
t405 = sin(qJ(1));
t408 = cos(qJ(1));
t390 = t405 * g(1) - g(2) * t408;
t384 = qJDD(1) * pkin(1) + t390;
t391 = -g(1) * t408 - g(2) * t405;
t410 = qJD(1) ^ 2;
t385 = -pkin(1) * t410 + t391;
t404 = sin(qJ(2));
t407 = cos(qJ(2));
t362 = t404 * t384 + t407 * t385;
t399 = qJDD(1) + qJDD(2);
t400 = (qJD(1) + qJD(2));
t440 = -t399 * qJ(3) - (2 * qJD(3) * t400) - t362;
t439 = -m(3) - m(4);
t438 = -pkin(2) - pkin(7);
t403 = sin(qJ(4));
t437 = t403 * g(3);
t436 = mrSges(3,1) - mrSges(4,2);
t435 = -mrSges(5,3) - mrSges(6,2);
t433 = Ifges(3,5) - Ifges(4,4);
t431 = (-Ifges(3,6) + Ifges(4,5));
t429 = t400 * t403;
t406 = cos(qJ(4));
t428 = t400 * t406;
t361 = t407 * t384 - t404 * t385;
t398 = t400 ^ 2;
t414 = -t398 * qJ(3) + qJDD(3) - t361;
t358 = t399 * t438 + t414;
t354 = -g(3) * t406 + t403 * t358;
t423 = qJD(4) * t400;
t377 = t399 * t403 + t406 * t423;
t387 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t428;
t375 = (mrSges(6,1) * t403 - mrSges(6,3) * t406) * t400;
t419 = t400 * (-t375 - (mrSges(5,1) * t403 + mrSges(5,2) * t406) * t400);
t374 = (pkin(4) * t403 - qJ(5) * t406) * t400;
t409 = qJD(4) ^ 2;
t351 = -pkin(4) * t409 + qJDD(4) * qJ(5) + (2 * qJD(5) * qJD(4)) - t374 * t429 + t354;
t388 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t428;
t422 = m(6) * t351 + qJDD(4) * mrSges(6,3) + qJD(4) * t388;
t345 = m(5) * t354 - qJDD(4) * mrSges(5,2) - qJD(4) * t387 + t377 * t435 + t403 * t419 + t422;
t353 = t358 * t406 + t437;
t378 = t399 * t406 - t403 * t423;
t386 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t429;
t352 = -qJDD(4) * pkin(4) - t437 - t409 * qJ(5) + qJDD(5) + (t374 * t400 - t358) * t406;
t389 = -mrSges(6,2) * t429 + qJD(4) * mrSges(6,3);
t415 = -m(6) * t352 + qJDD(4) * mrSges(6,1) + qJD(4) * t389;
t346 = m(5) * t353 + qJDD(4) * mrSges(5,1) + qJD(4) * t386 + t378 * t435 + t406 * t419 + t415;
t340 = t403 * t345 + t406 * t346;
t360 = -pkin(2) * t399 + t414;
t413 = -m(4) * t360 + (t398 * mrSges(4,3)) - t340;
t338 = m(3) * t361 - (t398 * mrSges(3,2)) + t399 * t436 + t413;
t359 = t398 * pkin(2) + t440;
t356 = t398 * t438 - t440;
t349 = t377 * pkin(4) - t378 * qJ(5) + (-0.2e1 * qJD(5) * t406 + (pkin(4) * t406 + qJ(5) * t403) * qJD(4)) * t400 + t356;
t347 = m(6) * t349 + t377 * mrSges(6,1) - mrSges(6,3) * t378 - t388 * t428 + t389 * t429;
t412 = -m(5) * t356 - mrSges(5,1) * t377 - t378 * mrSges(5,2) - t386 * t429 - t387 * t428 - t347;
t411 = -m(4) * t359 + (t398 * mrSges(4,2)) + t399 * mrSges(4,3) - t412;
t343 = m(3) * t362 - (mrSges(3,1) * t398) - mrSges(3,2) * t399 + t411;
t334 = t338 * t407 + t404 * t343;
t332 = m(2) * t390 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t410 + t334;
t417 = -t338 * t404 + t407 * t343;
t333 = m(2) * t391 - mrSges(2,1) * t410 - qJDD(1) * mrSges(2,2) + t417;
t427 = t332 * t408 + t333 * t405;
t426 = (t403 * t442 - t406 * t434) * t400 - (t430 * qJD(4));
t425 = (-t403 * t434 + t406 * t443) * t400 + (t432 * qJD(4));
t420 = t400 * ((t403 * t430 - t406 * t432) * t400 - (t441 * qJD(4)));
t418 = -t332 * t405 + t333 * t408;
t416 = -t406 * t345 + t403 * t346;
t339 = -m(4) * g(3) - t416;
t336 = mrSges(5,2) * t356 + mrSges(6,2) * t352 - mrSges(5,3) * t353 - mrSges(6,3) * t349 - qJ(5) * t347 + t426 * qJD(4) + t432 * qJDD(4) - t434 * t377 + t378 * t443 + t403 * t420;
t335 = -mrSges(5,1) * t356 - mrSges(6,1) * t349 + mrSges(6,2) * t351 + mrSges(5,3) * t354 - pkin(4) * t347 + t425 * qJD(4) + t430 * qJDD(4) - t377 * t442 + t434 * t378 + t406 * t420;
t328 = qJ(5) * t422 + pkin(4) * t415 + mrSges(6,3) * t351 - mrSges(6,1) * t352 + mrSges(5,1) * t353 - mrSges(5,2) * t354 + mrSges(4,1) * t360 - mrSges(3,3) * t361 - qJ(3) * t339 + pkin(3) * t340 + t433 * t399 + (t431 * t398) + (-mrSges(6,2) * pkin(4) + t432) * t378 + (-mrSges(6,2) * qJ(5) - t430) * t377 + t441 * qJDD(4) + (-mrSges(3,2) + mrSges(4,3)) * g(3) + ((-pkin(4) * t375 - t426) * t406 + (-qJ(5) * t375 + t425) * t403) * t400;
t327 = -mrSges(4,1) * t359 + mrSges(3,3) * t362 - pkin(2) * t339 - pkin(3) * t412 + pkin(7) * t416 + g(3) * t436 - t406 * t335 - t403 * t336 + t398 * t433 - t399 * t431;
t326 = -mrSges(2,2) * g(3) - mrSges(2,3) * t390 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t410 - pkin(6) * t334 - t327 * t404 + t328 * t407;
t325 = Ifges(2,6) * qJDD(1) + t410 * Ifges(2,5) + mrSges(2,3) * t391 + t404 * t328 + t407 * t327 + pkin(1) * t416 + pkin(6) * t417 + (-pkin(1) * t439 + mrSges(2,1)) * g(3);
t1 = [-m(1) * g(1) + t418; -m(1) * g(2) + t427; (-m(1) - m(2) + t439) * g(3) - t416; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t427 - t325 * t405 + t326 * t408; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t418 + t408 * t325 + t405 * t326; pkin(1) * t334 + mrSges(2,1) * t390 - mrSges(2,2) * t391 + pkin(2) * t413 + qJ(3) * t411 - mrSges(3,2) * t362 - t403 * t335 - pkin(7) * t340 + mrSges(3,1) * t361 + mrSges(4,2) * t360 - mrSges(4,3) * t359 + t406 * t336 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(2,3) * qJDD(1) + (-mrSges(4,2) * pkin(2) + Ifges(4,1) + Ifges(3,3)) * t399;];
tauB = t1;
