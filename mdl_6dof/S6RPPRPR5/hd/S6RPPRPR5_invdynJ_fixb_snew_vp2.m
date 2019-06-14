% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPRPR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-05-05 14:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPRPR5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:21:12
% EndTime: 2019-05-05 14:21:14
% DurationCPUTime: 1.22s
% Computational Cost: add. (8461->230), mult. (16983->286), div. (0->0), fcn. (9578->8), ass. (0->88)
t416 = qJD(1) ^ 2;
t411 = sin(qJ(1));
t414 = cos(qJ(1));
t434 = t411 * g(1) - t414 * g(2);
t380 = -qJDD(1) * pkin(1) - t416 * qJ(2) + qJDD(2) - t434;
t371 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t380;
t368 = -t416 * pkin(7) - t371;
t410 = sin(qJ(4));
t413 = cos(qJ(4));
t431 = qJD(1) * qJD(4);
t427 = t413 * t431;
t393 = t410 * qJDD(1) + t427;
t428 = t410 * t431;
t394 = t413 * qJDD(1) - t428;
t352 = (-t394 + t428) * qJ(5) + (t393 + t427) * pkin(4) + t368;
t424 = -t414 * g(1) - t411 * g(2);
t436 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t424;
t372 = qJDD(3) + (-pkin(1) - qJ(3)) * t416 + t436;
t369 = -qJDD(1) * pkin(7) + t372;
t363 = -t413 * g(3) + t410 * t369;
t391 = (t410 * pkin(4) - t413 * qJ(5)) * qJD(1);
t415 = qJD(4) ^ 2;
t433 = t410 * qJD(1);
t355 = -t415 * pkin(4) + qJDD(4) * qJ(5) - t391 * t433 + t363;
t407 = sin(pkin(9));
t408 = cos(pkin(9));
t432 = t413 * qJD(1);
t389 = t407 * qJD(4) + t408 * t432;
t339 = -0.2e1 * qJD(5) * t389 + t408 * t352 - t407 * t355;
t377 = t407 * qJDD(4) + t408 * t394;
t388 = t408 * qJD(4) - t407 * t432;
t337 = (t388 * t433 - t377) * pkin(8) + (t388 * t389 + t393) * pkin(5) + t339;
t340 = 0.2e1 * qJD(5) * t388 + t407 * t352 + t408 * t355;
t376 = t408 * qJDD(4) - t407 * t394;
t378 = pkin(5) * t433 - t389 * pkin(8);
t387 = t388 ^ 2;
t338 = -t387 * pkin(5) + t376 * pkin(8) - t378 * t433 + t340;
t409 = sin(qJ(6));
t412 = cos(qJ(6));
t335 = t412 * t337 - t409 * t338;
t364 = t412 * t388 - t409 * t389;
t344 = t364 * qJD(6) + t409 * t376 + t412 * t377;
t365 = t409 * t388 + t412 * t389;
t349 = -t364 * mrSges(7,1) + t365 * mrSges(7,2);
t397 = qJD(6) + t433;
t356 = -t397 * mrSges(7,2) + t364 * mrSges(7,3);
t390 = qJDD(6) + t393;
t332 = m(7) * t335 + t390 * mrSges(7,1) - t344 * mrSges(7,3) - t365 * t349 + t397 * t356;
t336 = t409 * t337 + t412 * t338;
t343 = -t365 * qJD(6) + t412 * t376 - t409 * t377;
t357 = t397 * mrSges(7,1) - t365 * mrSges(7,3);
t333 = m(7) * t336 - t390 * mrSges(7,2) + t343 * mrSges(7,3) + t364 * t349 - t397 * t357;
t325 = t412 * t332 + t409 * t333;
t366 = -t388 * mrSges(6,1) + t389 * mrSges(6,2);
t374 = -mrSges(6,2) * t433 + t388 * mrSges(6,3);
t323 = m(6) * t339 + t393 * mrSges(6,1) - t377 * mrSges(6,3) - t389 * t366 + t374 * t433 + t325;
t375 = mrSges(6,1) * t433 - t389 * mrSges(6,3);
t425 = -t409 * t332 + t412 * t333;
t324 = m(6) * t340 - t393 * mrSges(6,2) + t376 * mrSges(6,3) + t388 * t366 - t375 * t433 + t425;
t319 = t408 * t323 + t407 * t324;
t395 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t433;
t396 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t432;
t437 = m(4) * t371 - m(5) * t368 - t393 * mrSges(5,1) - t394 * mrSges(5,2) - t319 + (-t395 * t410 - t396 * t413) * qJD(1);
t362 = t410 * g(3) + t413 * t369;
t354 = -qJDD(4) * pkin(4) - t415 * qJ(5) + t391 * t432 + qJDD(5) - t362;
t341 = -t376 * pkin(5) - t387 * pkin(8) + t389 * t378 + t354;
t419 = m(7) * t341 - t343 * mrSges(7,1) + t344 * mrSges(7,2) - t364 * t356 + t365 * t357;
t334 = m(6) * t354 - t376 * mrSges(6,1) + t377 * mrSges(6,2) - t388 * t374 + t389 * t375 + t419;
t392 = (t410 * mrSges(5,1) + t413 * mrSges(5,2)) * qJD(1);
t426 = -t407 * t323 + t408 * t324;
t435 = t410 * (m(5) * t363 - qJDD(4) * mrSges(5,2) - t393 * mrSges(5,3) - qJD(4) * t396 - t392 * t433 + t426) + t413 * (m(5) * t362 + qJDD(4) * mrSges(5,1) - t394 * mrSges(5,3) + qJD(4) * t395 - t392 * t432 - t334);
t421 = m(4) * t372 + qJDD(1) * mrSges(4,2) - t416 * mrSges(4,3) + t435;
t346 = Ifges(7,4) * t365 + Ifges(7,2) * t364 + Ifges(7,6) * t397;
t347 = Ifges(7,1) * t365 + Ifges(7,4) * t364 + Ifges(7,5) * t397;
t417 = mrSges(7,1) * t335 - mrSges(7,2) * t336 + Ifges(7,5) * t344 + Ifges(7,6) * t343 + Ifges(7,3) * t390 + t365 * t346 - t364 * t347;
t386 = (Ifges(5,5) * qJD(4)) + (t413 * Ifges(5,1) - t410 * Ifges(5,4)) * qJD(1);
t385 = (Ifges(5,6) * qJD(4)) + (t413 * Ifges(5,4) - Ifges(5,2) * t410) * qJD(1);
t379 = t416 * pkin(1) - t436;
t360 = Ifges(6,1) * t389 + Ifges(6,4) * t388 + Ifges(6,5) * t433;
t359 = Ifges(6,4) * t389 + Ifges(6,2) * t388 + Ifges(6,6) * t433;
t358 = Ifges(6,5) * t389 + Ifges(6,6) * t388 + Ifges(6,3) * t433;
t345 = Ifges(7,5) * t365 + Ifges(7,6) * t364 + Ifges(7,3) * t397;
t327 = mrSges(7,2) * t341 - mrSges(7,3) * t335 + Ifges(7,1) * t344 + Ifges(7,4) * t343 + Ifges(7,5) * t390 + t364 * t345 - t397 * t346;
t326 = -mrSges(7,1) * t341 + mrSges(7,3) * t336 + Ifges(7,4) * t344 + Ifges(7,2) * t343 + Ifges(7,6) * t390 - t365 * t345 + t397 * t347;
t317 = m(3) * t380 + (-mrSges(4,2) - mrSges(3,3)) * t416 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t437;
t316 = mrSges(6,2) * t354 - mrSges(6,3) * t339 + Ifges(6,1) * t377 + Ifges(6,4) * t376 + Ifges(6,5) * t393 - pkin(8) * t325 - t409 * t326 + t412 * t327 + t388 * t358 - t359 * t433;
t315 = -mrSges(6,1) * t354 + mrSges(6,3) * t340 + Ifges(6,4) * t377 + Ifges(6,2) * t376 + Ifges(6,6) * t393 - pkin(5) * t419 + pkin(8) * t425 + t412 * t326 + t409 * t327 - t389 * t358 + t360 * t433;
t1 = [qJ(2) * (-m(3) * t379 + t416 * mrSges(3,2) + t421) - pkin(1) * t317 + mrSges(2,1) * t434 - mrSges(2,2) * t424 + mrSges(3,2) * t380 - mrSges(3,3) * t379 + t413 * (mrSges(5,2) * t368 - mrSges(5,3) * t362 + Ifges(5,1) * t394 - Ifges(5,4) * t393 + Ifges(5,5) * qJDD(4) - qJ(5) * t319 - qJD(4) * t385 - t407 * t315 + t408 * t316) - t410 * (-mrSges(5,1) * t368 - mrSges(6,1) * t339 + mrSges(6,2) * t340 + mrSges(5,3) * t363 + Ifges(5,4) * t394 - Ifges(6,5) * t377 + Ifges(5,6) * qJDD(4) - Ifges(6,6) * t376 - pkin(4) * t319 - pkin(5) * t325 + qJD(4) * t386 - t389 * t359 + t388 * t360 - t417 + (-Ifges(5,2) - Ifges(6,3)) * t393) - pkin(7) * t435 + mrSges(4,2) * t372 - mrSges(4,3) * t371 + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t416 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t437) * qJ(3); t317; t421; Ifges(5,5) * t394 - Ifges(5,6) * t393 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t362 - mrSges(5,2) * t363 + t407 * t316 + t408 * t315 - pkin(4) * t334 + qJ(5) * t426 + (t413 * t385 + t410 * t386) * qJD(1); t334; t417;];
tauJ  = t1;
