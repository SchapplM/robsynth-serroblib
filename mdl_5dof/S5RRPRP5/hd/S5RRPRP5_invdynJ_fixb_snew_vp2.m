% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:19
% EndTime: 2019-12-31 19:54:23
% DurationCPUTime: 2.05s
% Computational Cost: add. (12502->251), mult. (29177->314), div. (0->0), fcn. (19785->8), ass. (0->99)
t445 = Ifges(5,1) + Ifges(6,1);
t435 = Ifges(5,4) - Ifges(6,5);
t443 = Ifges(6,4) + Ifges(5,5);
t444 = Ifges(5,2) + Ifges(6,3);
t442 = Ifges(5,6) - Ifges(6,6);
t441 = Ifges(5,3) + Ifges(6,2);
t411 = sin(pkin(8));
t412 = cos(pkin(8));
t414 = sin(qJ(2));
t416 = cos(qJ(2));
t388 = (-t414 * t411 + t416 * t412) * qJD(1);
t389 = (t416 * t411 + t414 * t412) * qJD(1);
t413 = sin(qJ(4));
t438 = cos(qJ(4));
t371 = -t388 * t438 + t413 * t389;
t372 = t413 * t388 + t389 * t438;
t409 = qJD(2) + qJD(4);
t440 = t444 * t371 - t435 * t372 - t442 * t409;
t439 = -t435 * t371 + t445 * t372 + t443 * t409;
t418 = qJD(1) ^ 2;
t437 = pkin(2) * t418;
t436 = -mrSges(5,3) - mrSges(6,2);
t415 = sin(qJ(1));
t417 = cos(qJ(1));
t423 = -t417 * g(1) - t415 * g(2);
t394 = -t418 * pkin(1) + qJDD(1) * pkin(6) + t423;
t434 = t414 * t394;
t429 = qJD(1) * qJD(2);
t397 = t414 * qJDD(1) + t416 * t429;
t359 = qJDD(2) * pkin(2) - t397 * qJ(3) - t434 + (qJ(3) * t429 + t414 * t437 - g(3)) * t416;
t379 = -t414 * g(3) + t416 * t394;
t398 = t416 * qJDD(1) - t414 * t429;
t431 = qJD(1) * t414;
t399 = qJD(2) * pkin(2) - qJ(3) * t431;
t410 = t416 ^ 2;
t360 = t398 * qJ(3) - qJD(2) * t399 - t410 * t437 + t379;
t333 = -0.2e1 * qJD(3) * t389 + t412 * t359 - t411 * t360;
t377 = t412 * t397 + t411 * t398;
t328 = (qJD(2) * t388 - t377) * pkin(7) + (t388 * t389 + qJDD(2)) * pkin(3) + t333;
t334 = 0.2e1 * qJD(3) * t388 + t411 * t359 + t412 * t360;
t376 = -t411 * t397 + t412 * t398;
t382 = qJD(2) * pkin(3) - t389 * pkin(7);
t387 = t388 ^ 2;
t330 = -t387 * pkin(3) + t376 * pkin(7) - qJD(2) * t382 + t334;
t326 = t413 * t328 + t438 * t330;
t343 = t372 * qJD(4) - t376 * t438 + t413 * t377;
t365 = t409 * mrSges(5,1) - t372 * mrSges(5,3);
t408 = qJDD(2) + qJDD(4);
t353 = t371 * pkin(4) - t372 * qJ(5);
t407 = t409 ^ 2;
t320 = -t407 * pkin(4) + t408 * qJ(5) + 0.2e1 * qJD(5) * t409 - t371 * t353 + t326;
t366 = -t409 * mrSges(6,1) + t372 * mrSges(6,2);
t428 = m(6) * t320 + t408 * mrSges(6,3) + t409 * t366;
t354 = t371 * mrSges(6,1) - t372 * mrSges(6,3);
t432 = -t371 * mrSges(5,1) - t372 * mrSges(5,2) - t354;
t311 = m(5) * t326 - t408 * mrSges(5,2) + t343 * t436 - t409 * t365 + t371 * t432 + t428;
t325 = t328 * t438 - t413 * t330;
t344 = -t371 * qJD(4) + t413 * t376 + t377 * t438;
t364 = -t409 * mrSges(5,2) - t371 * mrSges(5,3);
t321 = -t408 * pkin(4) - t407 * qJ(5) + t372 * t353 + qJDD(5) - t325;
t367 = -t371 * mrSges(6,2) + t409 * mrSges(6,3);
t424 = -m(6) * t321 + t408 * mrSges(6,1) + t409 * t367;
t313 = m(5) * t325 + t408 * mrSges(5,1) + t344 * t436 + t409 * t364 + t372 * t432 + t424;
t307 = t413 * t311 + t438 * t313;
t374 = -t388 * mrSges(4,1) + t389 * mrSges(4,2);
t380 = -qJD(2) * mrSges(4,2) + t388 * mrSges(4,3);
t303 = m(4) * t333 + qJDD(2) * mrSges(4,1) - t377 * mrSges(4,3) + qJD(2) * t380 - t389 * t374 + t307;
t381 = qJD(2) * mrSges(4,1) - t389 * mrSges(4,3);
t425 = t438 * t311 - t413 * t313;
t304 = m(4) * t334 - qJDD(2) * mrSges(4,2) + t376 * mrSges(4,3) - qJD(2) * t381 + t388 * t374 + t425;
t299 = t412 * t303 + t411 * t304;
t433 = t442 * t371 - t443 * t372 - t441 * t409;
t430 = qJD(1) * t416;
t427 = t415 * g(1) - t417 * g(2);
t426 = -t411 * t303 + t412 * t304;
t422 = -qJDD(1) * pkin(1) - t427;
t363 = -t398 * pkin(2) + qJDD(3) + t399 * t431 + (-qJ(3) * t410 - pkin(6)) * t418 + t422;
t332 = -t376 * pkin(3) - t387 * pkin(7) + t389 * t382 + t363;
t323 = -0.2e1 * qJD(5) * t372 + (t371 * t409 - t344) * qJ(5) + (t372 * t409 + t343) * pkin(4) + t332;
t421 = -m(6) * t323 - t343 * mrSges(6,1) + t344 * mrSges(6,3) + t372 * t366 - t371 * t367;
t420 = m(5) * t332 + t343 * mrSges(5,1) + t344 * mrSges(5,2) + t371 * t364 + t372 * t365 - t421;
t317 = t344 * mrSges(6,2) + t372 * t354 - t424;
t419 = mrSges(5,1) * t325 - mrSges(6,1) * t321 - mrSges(5,2) * t326 + mrSges(6,3) * t320 - pkin(4) * t317 + qJ(5) * t428 + t441 * t408 - t440 * t372 + (-qJ(5) * t354 + t439) * t371 + t443 * t344 + (-qJ(5) * mrSges(6,2) - t442) * t343;
t308 = m(4) * t363 - t376 * mrSges(4,1) + t377 * mrSges(4,2) - t388 * t380 + t389 * t381 + t420;
t401 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t430;
t400 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t431;
t396 = (-t416 * mrSges(3,1) + t414 * mrSges(3,2)) * qJD(1);
t393 = -t418 * pkin(6) + t422;
t392 = Ifges(3,5) * qJD(2) + (t414 * Ifges(3,1) + t416 * Ifges(3,4)) * qJD(1);
t391 = Ifges(3,6) * qJD(2) + (t414 * Ifges(3,4) + t416 * Ifges(3,2)) * qJD(1);
t378 = -t416 * g(3) - t434;
t370 = Ifges(4,1) * t389 + Ifges(4,4) * t388 + Ifges(4,5) * qJD(2);
t369 = Ifges(4,4) * t389 + Ifges(4,2) * t388 + Ifges(4,6) * qJD(2);
t368 = Ifges(4,5) * t389 + Ifges(4,6) * t388 + Ifges(4,3) * qJD(2);
t306 = mrSges(5,2) * t332 + mrSges(6,2) * t321 - mrSges(5,3) * t325 - mrSges(6,3) * t323 + qJ(5) * t421 - t435 * t343 + t445 * t344 + t433 * t371 + t443 * t408 + t440 * t409;
t305 = -mrSges(5,1) * t332 - mrSges(6,1) * t323 + mrSges(6,2) * t320 + mrSges(5,3) * t326 + pkin(4) * t421 - t444 * t343 + t435 * t344 + t433 * t372 + t442 * t408 + t439 * t409;
t298 = mrSges(4,2) * t363 - mrSges(4,3) * t333 + Ifges(4,1) * t377 + Ifges(4,4) * t376 + Ifges(4,5) * qJDD(2) - pkin(7) * t307 - qJD(2) * t369 - t413 * t305 + t306 * t438 + t388 * t368;
t297 = -mrSges(4,1) * t363 + mrSges(4,3) * t334 + Ifges(4,4) * t377 + Ifges(4,2) * t376 + Ifges(4,6) * qJDD(2) - pkin(3) * t420 + pkin(7) * t425 + qJD(2) * t370 + t305 * t438 + t413 * t306 - t389 * t368;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t427 - mrSges(2,2) * t423 + t414 * (mrSges(3,2) * t393 - mrSges(3,3) * t378 + Ifges(3,1) * t397 + Ifges(3,4) * t398 + Ifges(3,5) * qJDD(2) - qJ(3) * t299 - qJD(2) * t391 - t411 * t297 + t412 * t298) + t416 * (-mrSges(3,1) * t393 + mrSges(3,3) * t379 + Ifges(3,4) * t397 + Ifges(3,2) * t398 + Ifges(3,6) * qJDD(2) - pkin(2) * t308 + qJ(3) * t426 + qJD(2) * t392 + t412 * t297 + t411 * t298) + pkin(1) * (-t308 - t397 * mrSges(3,2) + t398 * mrSges(3,1) - m(3) * t393 + (-t414 * t400 + t416 * t401) * qJD(1)) + pkin(6) * (t416 * (m(3) * t379 - qJDD(2) * mrSges(3,2) + t398 * mrSges(3,3) - qJD(2) * t400 + t396 * t430 + t426) - t414 * (m(3) * t378 + qJDD(2) * mrSges(3,1) - t397 * mrSges(3,3) + qJD(2) * t401 - t396 * t431 + t299)); t419 + Ifges(3,5) * t397 + Ifges(3,6) * t398 - t388 * t370 + t389 * t369 + Ifges(4,5) * t377 + mrSges(3,1) * t378 - mrSges(3,2) * t379 + Ifges(4,6) * t376 + mrSges(4,1) * t333 - mrSges(4,2) * t334 + pkin(3) * t307 + pkin(2) * t299 + (Ifges(3,3) + Ifges(4,3)) * qJDD(2) + (t414 * t391 - t416 * t392) * qJD(1); t308; t419; t317;];
tauJ = t1;
