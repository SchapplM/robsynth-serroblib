% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Datum: 2019-12-31 20:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:49:34
% EndTime: 2019-12-31 20:49:36
% DurationCPUTime: 0.96s
% Computational Cost: add. (7792->199), mult. (10302->247), div. (0->0), fcn. (5925->8), ass. (0->84)
t389 = Ifges(5,1) + Ifges(6,1);
t383 = Ifges(5,4) - Ifges(6,5);
t382 = Ifges(5,5) + Ifges(6,4);
t388 = -Ifges(5,2) - Ifges(6,3);
t387 = -Ifges(6,2) - Ifges(5,3);
t381 = Ifges(5,6) - Ifges(6,6);
t386 = -2 * qJD(4);
t350 = qJD(1) + qJD(2);
t348 = t350 ^ 2;
t385 = pkin(3) * t348;
t384 = -mrSges(5,3) - mrSges(6,2);
t380 = cos(pkin(8));
t355 = sin(qJ(3));
t379 = t350 * t355;
t358 = cos(qJ(3));
t378 = t350 * t358;
t357 = sin(qJ(1));
t360 = cos(qJ(1));
t370 = t357 * g(1) - t360 * g(2);
t340 = qJDD(1) * pkin(1) + t370;
t366 = -t360 * g(1) - t357 * g(2);
t341 = -qJD(1) ^ 2 * pkin(1) + t366;
t356 = sin(qJ(2));
t359 = cos(qJ(2));
t317 = t356 * t340 + t359 * t341;
t349 = qJDD(1) + qJDD(2);
t309 = -t348 * pkin(2) + t349 * pkin(7) + t317;
t377 = t355 * t309;
t372 = qJD(3) * t350;
t335 = t355 * t349 + t358 * t372;
t290 = qJDD(3) * pkin(3) - t335 * qJ(4) - t377 + (qJ(4) * t372 + t355 * t385 - g(3)) * t358;
t294 = -t355 * g(3) + t358 * t309;
t336 = t358 * t349 - t355 * t372;
t342 = qJD(3) * pkin(3) - qJ(4) * t379;
t353 = t358 ^ 2;
t291 = t336 * qJ(4) - qJD(3) * t342 - t353 * t385 + t294;
t354 = sin(pkin(8));
t325 = t354 * t379 - t380 * t378;
t287 = t354 * t290 + t380 * t291 + t325 * t386;
t313 = t354 * t335 - t380 * t336;
t326 = (t354 * t358 + t380 * t355) * t350;
t321 = qJD(3) * mrSges(5,1) - t326 * mrSges(5,3);
t305 = t325 * pkin(4) - t326 * qJ(5);
t361 = qJD(3) ^ 2;
t282 = -t361 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t325 * t305 + t287;
t322 = -qJD(3) * mrSges(6,1) + t326 * mrSges(6,2);
t371 = m(6) * t282 + qJDD(3) * mrSges(6,3) + qJD(3) * t322;
t306 = t325 * mrSges(6,1) - t326 * mrSges(6,3);
t373 = -t325 * mrSges(5,1) - t326 * mrSges(5,2) - t306;
t276 = m(5) * t287 - qJDD(3) * mrSges(5,2) - qJD(3) * t321 + t384 * t313 + t373 * t325 + t371;
t364 = t380 * t290 - t354 * t291;
t286 = t326 * t386 + t364;
t314 = t380 * t335 + t354 * t336;
t320 = -qJD(3) * mrSges(5,2) - t325 * mrSges(5,3);
t283 = -qJDD(3) * pkin(4) - t361 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t305) * t326 - t364;
t323 = -t325 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t367 = -m(6) * t283 + qJDD(3) * mrSges(6,1) + qJD(3) * t323;
t277 = m(5) * t286 + qJDD(3) * mrSges(5,1) + qJD(3) * t320 + t384 * t314 + t373 * t326 + t367;
t271 = t354 * t276 + t380 * t277;
t376 = t387 * qJD(3) + t381 * t325 - t382 * t326;
t375 = t381 * qJD(3) + t388 * t325 + t383 * t326;
t374 = t382 * qJD(3) - t383 * t325 + t389 * t326;
t293 = -t358 * g(3) - t377;
t334 = (-mrSges(4,1) * t358 + mrSges(4,2) * t355) * t350;
t343 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t379;
t344 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t378;
t368 = t380 * t276 - t354 * t277;
t369 = -t355 * (m(4) * t293 + qJDD(3) * mrSges(4,1) - t335 * mrSges(4,3) + qJD(3) * t344 - t334 * t379 + t271) + t358 * (m(4) * t294 - qJDD(3) * mrSges(4,2) + t336 * mrSges(4,3) - qJD(3) * t343 + t334 * t378 + t368);
t316 = t359 * t340 - t356 * t341;
t365 = -t349 * pkin(2) - t316;
t292 = -t336 * pkin(3) + qJDD(4) + t342 * t379 + (-qJ(4) * t353 - pkin(7)) * t348 + t365;
t285 = -0.2e1 * qJD(5) * t326 + (qJD(3) * t325 - t314) * qJ(5) + (qJD(3) * t326 + t313) * pkin(4) + t292;
t280 = m(6) * t285 + t313 * mrSges(6,1) - t314 * mrSges(6,3) - t326 * t322 + t325 * t323;
t269 = -mrSges(5,1) * t292 - mrSges(6,1) * t285 + mrSges(6,2) * t282 + mrSges(5,3) * t287 - pkin(4) * t280 + t374 * qJD(3) + t381 * qJDD(3) + t388 * t313 + t383 * t314 + t376 * t326;
t270 = mrSges(5,2) * t292 + mrSges(6,2) * t283 - mrSges(5,3) * t286 - mrSges(6,3) * t285 - qJ(5) * t280 - t375 * qJD(3) + t382 * qJDD(3) - t383 * t313 + t389 * t314 + t376 * t325;
t278 = m(5) * t292 + t313 * mrSges(5,1) + t314 * mrSges(5,2) + t325 * t320 + t326 * t321 + t280;
t308 = -t348 * pkin(7) + t365;
t328 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t355 + Ifges(4,6) * t358) * t350;
t329 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t355 + Ifges(4,2) * t358) * t350;
t330 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t355 + Ifges(4,4) * t358) * t350;
t362 = -m(4) * t308 + t336 * mrSges(4,1) - t335 * mrSges(4,2) - t343 * t379 + t344 * t378 - t278;
t363 = -mrSges(3,2) * t317 + t358 * (-mrSges(4,1) * t308 + mrSges(4,3) * t294 + Ifges(4,4) * t335 + Ifges(4,2) * t336 + Ifges(4,6) * qJDD(3) - pkin(3) * t278 + qJ(4) * t368 + qJD(3) * t330 + t380 * t269 + t354 * t270 - t328 * t379) + t355 * (mrSges(4,2) * t308 - mrSges(4,3) * t293 + Ifges(4,1) * t335 + Ifges(4,4) * t336 + Ifges(4,5) * qJDD(3) - qJ(4) * t271 - qJD(3) * t329 - t354 * t269 + t380 * t270 + t328 * t378) + pkin(7) * t369 + pkin(2) * t362 + mrSges(3,1) * t316 + Ifges(3,3) * t349;
t279 = t314 * mrSges(6,2) + t326 * t306 - t367;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t370 - mrSges(2,2) * t366 + pkin(1) * (t356 * (m(3) * t317 - t348 * mrSges(3,1) - t349 * mrSges(3,2) + t369) + t359 * (m(3) * t316 + t349 * mrSges(3,1) - t348 * mrSges(3,2) + t362)) + t363; t363; Ifges(4,5) * t335 + Ifges(4,6) * t336 + mrSges(4,1) * t293 - mrSges(4,2) * t294 + mrSges(5,1) * t286 - mrSges(5,2) * t287 - mrSges(6,1) * t283 + mrSges(6,3) * t282 - pkin(4) * t279 + qJ(5) * t371 + pkin(3) * t271 + (t355 * t329 - t358 * t330) * t350 + t375 * t326 + (-qJ(5) * t306 + t374) * t325 + t382 * t314 + (-qJ(5) * mrSges(6,2) - t381) * t313 + (Ifges(4,3) - t387) * qJDD(3); t278; t279;];
tauJ = t1;
