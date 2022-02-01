% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2022-01-20 11:44
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPR3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:42:46
% EndTime: 2022-01-20 11:42:48
% DurationCPUTime: 1.54s
% Computational Cost: add. (20826->222), mult. (27803->290), div. (0->0), fcn. (17752->10), ass. (0->92)
t367 = qJD(1) + qJD(2);
t363 = t367 ^ 2;
t393 = pkin(3) * t363;
t372 = sin(qJ(3));
t392 = t367 * t372;
t376 = cos(qJ(3));
t391 = t367 * t376;
t374 = sin(qJ(1));
t378 = cos(qJ(1));
t388 = t374 * g(1) - t378 * g(2);
t355 = qJDD(1) * pkin(1) + t388;
t384 = -t378 * g(1) - t374 * g(2);
t356 = -qJD(1) ^ 2 * pkin(1) + t384;
t373 = sin(qJ(2));
t377 = cos(qJ(2));
t334 = t373 * t355 + t377 * t356;
t365 = qJDD(1) + qJDD(2);
t328 = -t363 * pkin(2) + t365 * pkin(7) + t334;
t390 = t372 * t328;
t389 = qJD(3) * t367;
t350 = t372 * t365 + t376 * t389;
t312 = qJDD(3) * pkin(3) - t350 * qJ(4) - t390 + (qJ(4) * t389 + t372 * t393 - g(3)) * t376;
t318 = -t372 * g(3) + t376 * t328;
t351 = t376 * t365 - t372 * t389;
t357 = qJD(3) * pkin(3) - qJ(4) * t392;
t368 = t376 ^ 2;
t313 = t351 * qJ(4) - qJD(3) * t357 - t368 * t393 + t318;
t369 = sin(pkin(9));
t370 = cos(pkin(9));
t342 = (t369 * t376 + t370 * t372) * t367;
t292 = -0.2e1 * qJD(4) * t342 + t370 * t312 - t369 * t313;
t331 = t370 * t350 + t369 * t351;
t341 = (-t369 * t372 + t370 * t376) * t367;
t290 = (qJD(3) * t341 - t331) * pkin(8) + (t341 * t342 + qJDD(3)) * pkin(4) + t292;
t293 = 0.2e1 * qJD(4) * t341 + t369 * t312 + t370 * t313;
t330 = -t369 * t350 + t370 * t351;
t337 = qJD(3) * pkin(4) - t342 * pkin(8);
t340 = t341 ^ 2;
t291 = -t340 * pkin(4) + t330 * pkin(8) - qJD(3) * t337 + t293;
t371 = sin(qJ(5));
t375 = cos(qJ(5));
t288 = t375 * t290 - t371 * t291;
t322 = t375 * t341 - t371 * t342;
t302 = t322 * qJD(5) + t371 * t330 + t375 * t331;
t323 = t371 * t341 + t375 * t342;
t308 = -t322 * mrSges(6,1) + t323 * mrSges(6,2);
t366 = qJD(3) + qJD(5);
t315 = -t366 * mrSges(6,2) + t322 * mrSges(6,3);
t364 = qJDD(3) + qJDD(5);
t284 = m(6) * t288 + t364 * mrSges(6,1) - t302 * mrSges(6,3) - t323 * t308 + t366 * t315;
t289 = t371 * t290 + t375 * t291;
t301 = -t323 * qJD(5) + t375 * t330 - t371 * t331;
t316 = t366 * mrSges(6,1) - t323 * mrSges(6,3);
t285 = m(6) * t289 - t364 * mrSges(6,2) + t301 * mrSges(6,3) + t322 * t308 - t366 * t316;
t277 = t375 * t284 + t371 * t285;
t326 = -t341 * mrSges(5,1) + t342 * mrSges(5,2);
t335 = -qJD(3) * mrSges(5,2) + t341 * mrSges(5,3);
t275 = m(5) * t292 + qJDD(3) * mrSges(5,1) - t331 * mrSges(5,3) + qJD(3) * t335 - t342 * t326 + t277;
t336 = qJD(3) * mrSges(5,1) - t342 * mrSges(5,3);
t385 = -t371 * t284 + t375 * t285;
t276 = m(5) * t293 - qJDD(3) * mrSges(5,2) + t330 * mrSges(5,3) - qJD(3) * t336 + t341 * t326 + t385;
t271 = t370 * t275 + t369 * t276;
t317 = -t376 * g(3) - t390;
t349 = (-mrSges(4,1) * t376 + mrSges(4,2) * t372) * t367;
t358 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t392;
t359 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t391;
t386 = -t369 * t275 + t370 * t276;
t387 = -t372 * (m(4) * t317 + qJDD(3) * mrSges(4,1) - t350 * mrSges(4,3) + qJD(3) * t359 - t349 * t392 + t271) + t376 * (m(4) * t318 - qJDD(3) * mrSges(4,2) + t351 * mrSges(4,3) - qJD(3) * t358 + t349 * t391 + t386);
t333 = t377 * t355 - t373 * t356;
t382 = -t365 * pkin(2) - t333;
t314 = -t351 * pkin(3) + qJDD(4) + t357 * t392 + (-qJ(4) * t368 - pkin(7)) * t363 + t382;
t295 = -t330 * pkin(4) - t340 * pkin(8) + t342 * t337 + t314;
t383 = m(6) * t295 - t301 * mrSges(6,1) + t302 * mrSges(6,2) - t322 * t315 + t323 * t316;
t303 = Ifges(6,5) * t323 + Ifges(6,6) * t322 + Ifges(6,3) * t366;
t305 = Ifges(6,1) * t323 + Ifges(6,4) * t322 + Ifges(6,5) * t366;
t278 = -mrSges(6,1) * t295 + mrSges(6,3) * t289 + Ifges(6,4) * t302 + Ifges(6,2) * t301 + Ifges(6,6) * t364 - t323 * t303 + t366 * t305;
t304 = Ifges(6,4) * t323 + Ifges(6,2) * t322 + Ifges(6,6) * t366;
t279 = mrSges(6,2) * t295 - mrSges(6,3) * t288 + Ifges(6,1) * t302 + Ifges(6,4) * t301 + Ifges(6,5) * t364 + t322 * t303 - t366 * t304;
t319 = Ifges(5,5) * t342 + Ifges(5,6) * t341 + Ifges(5,3) * qJD(3);
t321 = Ifges(5,1) * t342 + Ifges(5,4) * t341 + Ifges(5,5) * qJD(3);
t267 = -mrSges(5,1) * t314 + mrSges(5,3) * t293 + Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * qJDD(3) - pkin(4) * t383 + pkin(8) * t385 + qJD(3) * t321 + t375 * t278 + t371 * t279 - t342 * t319;
t320 = Ifges(5,4) * t342 + Ifges(5,2) * t341 + Ifges(5,6) * qJD(3);
t268 = mrSges(5,2) * t314 - mrSges(5,3) * t292 + Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * qJDD(3) - pkin(8) * t277 - qJD(3) * t320 - t371 * t278 + t375 * t279 + t341 * t319;
t286 = m(5) * t314 - t330 * mrSges(5,1) + t331 * mrSges(5,2) - t341 * t335 + t342 * t336 + t383;
t327 = -t363 * pkin(7) + t382;
t343 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t372 + Ifges(4,6) * t376) * t367;
t344 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t372 + Ifges(4,2) * t376) * t367;
t345 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t372 + Ifges(4,4) * t376) * t367;
t379 = -m(4) * t327 + t351 * mrSges(4,1) - t350 * mrSges(4,2) - t358 * t392 + t359 * t391 - t286;
t381 = -mrSges(3,2) * t334 + t376 * (-mrSges(4,1) * t327 + mrSges(4,3) * t318 + Ifges(4,4) * t350 + Ifges(4,2) * t351 + Ifges(4,6) * qJDD(3) - pkin(3) * t286 + qJ(4) * t386 + qJD(3) * t345 + t370 * t267 + t369 * t268 - t343 * t392) + t372 * (mrSges(4,2) * t327 - mrSges(4,3) * t317 + Ifges(4,1) * t350 + Ifges(4,4) * t351 + Ifges(4,5) * qJDD(3) - qJ(4) * t271 - qJD(3) * t344 - t369 * t267 + t370 * t268 + t343 * t391) + pkin(7) * t387 + pkin(2) * t379 + mrSges(3,1) * t333 + Ifges(3,3) * t365;
t380 = mrSges(6,1) * t288 - mrSges(6,2) * t289 + Ifges(6,5) * t302 + Ifges(6,6) * t301 + Ifges(6,3) * t364 + t323 * t304 - t322 * t305;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t388 - mrSges(2,2) * t384 + pkin(1) * (t373 * (m(3) * t334 - t363 * mrSges(3,1) - t365 * mrSges(3,2) + t387) + t377 * (m(3) * t333 + t365 * mrSges(3,1) - t363 * mrSges(3,2) + t379)) + t381; t381; t380 + (t372 * t344 - t376 * t345) * t367 + (Ifges(4,3) + Ifges(5,3)) * qJDD(3) + Ifges(4,5) * t350 + Ifges(4,6) * t351 - t341 * t321 + t342 * t320 + Ifges(5,6) * t330 + Ifges(5,5) * t331 + mrSges(4,1) * t317 - mrSges(4,2) * t318 + mrSges(5,1) * t292 - mrSges(5,2) * t293 + pkin(4) * t277 + pkin(3) * t271; t286; t380;];
tauJ = t1;
