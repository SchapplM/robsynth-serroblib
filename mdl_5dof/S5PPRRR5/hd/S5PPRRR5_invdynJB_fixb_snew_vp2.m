% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S5PPRRR5
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
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1]';
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
% tauJB [(6+5)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S5PPRRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR5_invdynJB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:35:43
% EndTime: 2019-12-31 17:35:44
% DurationCPUTime: 0.82s
% Computational Cost: add. (9836->158), mult. (12550->196), div. (0->0), fcn. (7568->8), ass. (0->72)
t388 = sin(pkin(8));
t389 = cos(pkin(8));
t379 = t388 * g(1) - t389 * g(2);
t377 = qJDD(2) - t379;
t380 = -t389 * g(1) - t388 * g(2);
t392 = sin(qJ(3));
t395 = cos(qJ(3));
t362 = t395 * t377 - t392 * t380;
t360 = qJDD(3) * pkin(3) + t362;
t363 = t392 * t377 + t395 * t380;
t396 = qJD(3) ^ 2;
t361 = -t396 * pkin(3) + t363;
t391 = sin(qJ(4));
t394 = cos(qJ(4));
t357 = t391 * t360 + t394 * t361;
t386 = qJD(3) + qJD(4);
t384 = t386 ^ 2;
t385 = qJDD(3) + qJDD(4);
t354 = -t384 * pkin(4) + t385 * pkin(7) + t357;
t387 = g(3) - qJDD(1);
t390 = sin(qJ(5));
t393 = cos(qJ(5));
t351 = -t390 * t354 + t393 * t387;
t352 = t393 * t354 + t390 * t387;
t365 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t390 + Ifges(6,2) * t393) * t386;
t366 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t390 + Ifges(6,4) * t393) * t386;
t408 = qJD(5) * t386;
t370 = t390 * t385 + t393 * t408;
t371 = t393 * t385 - t390 * t408;
t414 = mrSges(6,1) * t351 - mrSges(6,2) * t352 + Ifges(6,5) * t370 + Ifges(6,6) * t371 + Ifges(6,3) * qJDD(5) + (t365 * t390 - t366 * t393) * t386;
t413 = -m(4) - m(5);
t412 = -mrSges(2,2) + mrSges(3,3);
t411 = t386 * t390;
t410 = t386 * t393;
t369 = (-mrSges(6,1) * t393 + mrSges(6,2) * t390) * t386;
t376 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t410;
t349 = m(6) * t351 + qJDD(5) * mrSges(6,1) - t370 * mrSges(6,3) + qJD(5) * t376 - t369 * t411;
t375 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t411;
t350 = m(6) * t352 - qJDD(5) * mrSges(6,2) + t371 * mrSges(6,3) - qJD(5) * t375 + t369 * t410;
t404 = -t390 * t349 + t393 * t350;
t333 = m(5) * t357 - t384 * mrSges(5,1) - t385 * mrSges(5,2) + t404;
t356 = t394 * t360 - t391 * t361;
t353 = -t385 * pkin(4) - t384 * pkin(7) - t356;
t399 = -m(6) * t353 + t371 * mrSges(6,1) - t370 * mrSges(6,2) - t375 * t411 + t376 * t410;
t343 = m(5) * t356 + t385 * mrSges(5,1) - t384 * mrSges(5,2) + t399;
t330 = t391 * t333 + t394 * t343;
t328 = m(4) * t362 + qJDD(3) * mrSges(4,1) - t396 * mrSges(4,2) + t330;
t405 = t394 * t333 - t391 * t343;
t329 = m(4) * t363 - t396 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t405;
t324 = t395 * t328 + t392 * t329;
t323 = m(3) * t377 + t324;
t321 = m(2) * t379 - t323;
t406 = -t392 * t328 + t395 * t329;
t403 = m(3) * t380 + t406;
t322 = m(2) * t380 + t403;
t409 = t389 * t321 + t388 * t322;
t336 = t393 * t349 + t390 * t350;
t407 = -t388 * t321 + t389 * t322;
t334 = -t336 + (-m(3) + t413) * t387;
t401 = -m(2) * t387 + t334;
t364 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t390 + Ifges(6,6) * t393) * t386;
t340 = -mrSges(6,1) * t353 + mrSges(6,3) * t352 + Ifges(6,4) * t370 + Ifges(6,2) * t371 + Ifges(6,6) * qJDD(5) + qJD(5) * t366 - t364 * t411;
t341 = mrSges(6,2) * t353 - mrSges(6,3) * t351 + Ifges(6,1) * t370 + Ifges(6,4) * t371 + Ifges(6,5) * qJDD(5) - qJD(5) * t365 + t364 * t410;
t400 = -mrSges(5,1) * t356 + mrSges(5,2) * t357 - Ifges(5,3) * t385 - pkin(4) * t399 - pkin(7) * t404 - t393 * t340 - t390 * t341;
t397 = mrSges(4,1) * t362 - mrSges(4,2) * t363 + Ifges(4,3) * qJDD(3) + pkin(3) * t330 - t400;
t327 = -mrSges(5,1) * t387 + mrSges(5,3) * t357 + t384 * Ifges(5,5) + Ifges(5,6) * t385 - pkin(4) * t336 - t414;
t325 = mrSges(5,2) * t387 - mrSges(5,3) * t356 + Ifges(5,5) * t385 - t384 * Ifges(5,6) - pkin(7) * t336 - t390 * t340 + t393 * t341;
t317 = mrSges(4,2) * t387 - mrSges(4,3) * t362 + Ifges(4,5) * qJDD(3) - t396 * Ifges(4,6) - pkin(6) * t330 + t394 * t325 - t391 * t327;
t316 = Ifges(4,6) * qJDD(3) + t396 * Ifges(4,5) - mrSges(4,1) * t387 + mrSges(4,3) * t363 + t391 * t325 + t394 * t327 - pkin(3) * (m(5) * t387 + t336) + pkin(6) * t405;
t315 = mrSges(3,2) * t377 - mrSges(2,3) * t379 - pkin(5) * t324 - qJ(2) * t334 - t392 * t316 + t395 * t317 + t412 * t387;
t314 = -t392 * t317 - t395 * t316 + pkin(2) * t336 - pkin(5) * t406 - pkin(1) * t334 + (mrSges(3,2) + mrSges(2,3)) * t380 + (-pkin(2) * t413 + mrSges(2,1) + mrSges(3,1)) * t387;
t1 = [-m(1) * g(1) + t407; -m(1) * g(2) + t409; -m(1) * g(3) + t401; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t409 - t388 * t314 + t389 * t315; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t407 + t389 * t314 + t388 * t315; -mrSges(1,1) * g(2) + mrSges(2,1) * t379 - mrSges(3,1) * t377 + mrSges(1,2) * g(1) - pkin(1) * t323 - pkin(2) * t324 + qJ(2) * t403 + t412 * t380 - t397; t401; t323; t397; -t400; t414;];
tauJB = t1;
