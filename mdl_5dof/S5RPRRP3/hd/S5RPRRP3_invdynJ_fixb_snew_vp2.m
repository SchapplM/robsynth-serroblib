% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP3
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% tauJ [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:30
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:29:45
% EndTime: 2022-01-23 09:29:46
% DurationCPUTime: 1.10s
% Computational Cost: add. (4142->197), mult. (8267->242), div. (0->0), fcn. (4901->8), ass. (0->82)
t393 = Ifges(5,4) + Ifges(6,4);
t399 = Ifges(5,2) + Ifges(6,2);
t396 = Ifges(5,6) + Ifges(6,6);
t398 = Ifges(5,1) + Ifges(6,1);
t397 = Ifges(5,5) + Ifges(6,5);
t395 = Ifges(5,3) + Ifges(6,3);
t370 = sin(qJ(4));
t371 = sin(qJ(3));
t373 = cos(qJ(4));
t374 = cos(qJ(3));
t345 = (-t371 * t370 + t374 * t373) * qJD(1);
t346 = (t374 * t370 + t371 * t373) * qJD(1);
t365 = qJD(3) + qJD(4);
t394 = t399 * t345 + t393 * t346 + t396 * t365;
t372 = sin(qJ(1));
t375 = cos(qJ(1));
t384 = t372 * g(1) - t375 * g(2);
t350 = qJDD(1) * pkin(1) + t384;
t376 = qJD(1) ^ 2;
t381 = -t375 * g(1) - t372 * g(2);
t352 = -t376 * pkin(1) + t381;
t368 = sin(pkin(8));
t369 = cos(pkin(8));
t331 = t368 * t350 + t369 * t352;
t327 = -t376 * pkin(2) + qJDD(1) * pkin(6) + t331;
t367 = -g(3) + qJDD(2);
t312 = -t371 * t327 + t374 * t367;
t388 = qJD(1) * qJD(3);
t385 = t374 * t388;
t353 = t371 * qJDD(1) + t385;
t298 = (-t353 + t385) * pkin(7) + (t371 * t374 * t376 + qJDD(3)) * pkin(3) + t312;
t313 = t374 * t327 + t371 * t367;
t354 = t374 * qJDD(1) - t371 * t388;
t389 = t371 * qJD(1);
t357 = qJD(3) * pkin(3) - pkin(7) * t389;
t366 = t374 ^ 2;
t299 = -t366 * t376 * pkin(3) + t354 * pkin(7) - qJD(3) * t357 + t313;
t293 = t373 * t298 - t370 * t299;
t315 = t345 * qJD(4) + t373 * t353 + t370 * t354;
t328 = -t345 * mrSges(6,1) + t346 * mrSges(6,2);
t329 = -t345 * mrSges(5,1) + t346 * mrSges(5,2);
t334 = -t365 * mrSges(5,2) + t345 * mrSges(5,3);
t364 = qJDD(3) + qJDD(4);
t287 = -0.2e1 * qJD(5) * t346 + (t345 * t365 - t315) * qJ(5) + (t345 * t346 + t364) * pkin(4) + t293;
t333 = -t365 * mrSges(6,2) + t345 * mrSges(6,3);
t387 = m(6) * t287 + t364 * mrSges(6,1) + t365 * t333;
t278 = m(5) * t293 + t364 * mrSges(5,1) + t365 * t334 + (-t328 - t329) * t346 + (-mrSges(5,3) - mrSges(6,3)) * t315 + t387;
t294 = t370 * t298 + t373 * t299;
t314 = -t346 * qJD(4) - t370 * t353 + t373 * t354;
t336 = t365 * mrSges(6,1) - t346 * mrSges(6,3);
t337 = t365 * mrSges(5,1) - t346 * mrSges(5,3);
t335 = t365 * pkin(4) - t346 * qJ(5);
t338 = t345 ^ 2;
t289 = -t338 * pkin(4) + t314 * qJ(5) + 0.2e1 * qJD(5) * t345 - t365 * t335 + t294;
t386 = m(6) * t289 + t314 * mrSges(6,3) + t345 * t328;
t281 = m(5) * t294 + t314 * mrSges(5,3) + t345 * t329 + (-t336 - t337) * t365 + (-mrSges(5,2) - mrSges(6,2)) * t364 + t386;
t276 = t373 * t278 + t370 * t281;
t392 = -t396 * t345 - t397 * t346 - t395 * t365;
t391 = -t393 * t345 - t398 * t346 - t397 * t365;
t390 = qJD(1) * t374;
t351 = (-t374 * mrSges(4,1) + t371 * mrSges(4,2)) * qJD(1);
t356 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t390;
t274 = m(4) * t312 + qJDD(3) * mrSges(4,1) - t353 * mrSges(4,3) + qJD(3) * t356 - t351 * t389 + t276;
t355 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t389;
t382 = -t370 * t278 + t373 * t281;
t275 = m(4) * t313 - qJDD(3) * mrSges(4,2) + t354 * mrSges(4,3) - qJD(3) * t355 + t351 * t390 + t382;
t383 = -t371 * t274 + t374 * t275;
t330 = t369 * t350 - t368 * t352;
t380 = -qJDD(1) * pkin(2) - t330;
t300 = -t354 * pkin(3) + t357 * t389 + (-pkin(7) * t366 - pkin(6)) * t376 + t380;
t291 = -t314 * pkin(4) - t338 * qJ(5) + t346 * t335 + qJDD(5) + t300;
t284 = m(6) * t291 - t314 * mrSges(6,1) + t315 * mrSges(6,2) - t345 * t333 + t346 * t336;
t379 = m(5) * t300 - t314 * mrSges(5,1) + t315 * mrSges(5,2) - t345 * t334 + t346 * t337 + t284;
t283 = -t315 * mrSges(6,3) - t346 * t328 + t387;
t378 = mrSges(5,1) * t293 + mrSges(6,1) * t287 - mrSges(5,2) * t294 - mrSges(6,2) * t289 + pkin(4) * t283 + t396 * t314 + t397 * t315 + t391 * t345 + t394 * t346 + t395 * t364;
t326 = -t376 * pkin(6) + t380;
t377 = -m(4) * t326 + t354 * mrSges(4,1) - t353 * mrSges(4,2) - t355 * t389 + t356 * t390 - t379;
t344 = Ifges(4,5) * qJD(3) + (t371 * Ifges(4,1) + t374 * Ifges(4,4)) * qJD(1);
t343 = Ifges(4,6) * qJD(3) + (t371 * Ifges(4,4) + t374 * Ifges(4,2)) * qJD(1);
t272 = mrSges(5,2) * t300 + mrSges(6,2) * t291 - mrSges(5,3) * t293 - mrSges(6,3) * t287 - qJ(5) * t283 + t393 * t314 + t398 * t315 - t392 * t345 + t397 * t364 - t394 * t365;
t271 = -mrSges(5,1) * t300 + mrSges(5,3) * t294 - mrSges(6,1) * t291 + mrSges(6,3) * t289 - pkin(4) * t284 + qJ(5) * t386 + (-qJ(5) * t336 - t391) * t365 + (-qJ(5) * mrSges(6,2) + t396) * t364 + t392 * t346 + t393 * t315 + t399 * t314;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t384 - mrSges(2,2) * t381 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t330 - mrSges(3,2) * t331 + t371 * (mrSges(4,2) * t326 - mrSges(4,3) * t312 + Ifges(4,1) * t353 + Ifges(4,4) * t354 + Ifges(4,5) * qJDD(3) - pkin(7) * t276 - qJD(3) * t343 - t370 * t271 + t373 * t272) + t374 * (-mrSges(4,1) * t326 + mrSges(4,3) * t313 + Ifges(4,4) * t353 + Ifges(4,2) * t354 + Ifges(4,6) * qJDD(3) - pkin(3) * t379 + pkin(7) * t382 + qJD(3) * t344 + t373 * t271 + t370 * t272) + pkin(2) * t377 + pkin(6) * t383 + pkin(1) * (t368 * (m(3) * t331 - t376 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t383) + t369 * (m(3) * t330 + qJDD(1) * mrSges(3,1) - t376 * mrSges(3,2) + t377)); m(3) * t367 + t374 * t274 + t371 * t275; (t371 * t343 - t374 * t344) * qJD(1) + mrSges(4,1) * t312 - mrSges(4,2) * t313 + t378 + Ifges(4,3) * qJDD(3) + pkin(3) * t276 + Ifges(4,5) * t353 + Ifges(4,6) * t354; t378; t284;];
tauJ = t1;
