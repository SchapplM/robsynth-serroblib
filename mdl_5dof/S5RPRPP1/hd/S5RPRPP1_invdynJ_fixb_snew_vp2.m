% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRPP1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:48
% EndTime: 2019-12-31 18:08:49
% DurationCPUTime: 0.82s
% Computational Cost: add. (3451->195), mult. (7276->242), div. (0->0), fcn. (4179->8), ass. (0->82)
t382 = Ifges(5,1) + Ifges(6,1);
t377 = Ifges(5,4) - Ifges(6,5);
t376 = Ifges(5,5) + Ifges(6,4);
t381 = -Ifges(5,2) - Ifges(6,3);
t380 = -Ifges(6,2) - Ifges(5,3);
t375 = Ifges(5,6) - Ifges(6,6);
t379 = -2 * qJD(4);
t378 = -mrSges(5,3) - mrSges(6,2);
t374 = cos(pkin(8));
t352 = sin(qJ(1));
t354 = cos(qJ(1));
t364 = t352 * g(1) - t354 * g(2);
t332 = qJDD(1) * pkin(1) + t364;
t356 = qJD(1) ^ 2;
t360 = -t354 * g(1) - t352 * g(2);
t334 = -t356 * pkin(1) + t360;
t349 = sin(pkin(7));
t350 = cos(pkin(7));
t310 = t349 * t332 + t350 * t334;
t299 = -t356 * pkin(2) + qJDD(1) * pkin(6) + t310;
t347 = -g(3) + qJDD(2);
t351 = sin(qJ(3));
t353 = cos(qJ(3));
t289 = -t351 * t299 + t353 * t347;
t367 = qJD(1) * qJD(3);
t365 = t353 * t367;
t335 = t351 * qJDD(1) + t365;
t286 = (-t335 + t365) * qJ(4) + (t351 * t353 * t356 + qJDD(3)) * pkin(3) + t289;
t290 = t353 * t299 + t351 * t347;
t336 = t353 * qJDD(1) - t351 * t367;
t368 = t351 * qJD(1);
t337 = qJD(3) * pkin(3) - qJ(4) * t368;
t346 = t353 ^ 2;
t287 = -t346 * t356 * pkin(3) + t336 * qJ(4) - qJD(3) * t337 + t290;
t348 = sin(pkin(8));
t369 = qJD(1) * t353;
t320 = t348 * t368 - t374 * t369;
t283 = t348 * t286 + t374 * t287 + t320 * t379;
t311 = t348 * t335 - t374 * t336;
t321 = (t348 * t353 + t374 * t351) * qJD(1);
t316 = qJD(3) * mrSges(5,1) - t321 * mrSges(5,3);
t303 = t320 * pkin(4) - t321 * qJ(5);
t355 = qJD(3) ^ 2;
t278 = -t355 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t320 * t303 + t283;
t317 = -qJD(3) * mrSges(6,1) + t321 * mrSges(6,2);
t366 = m(6) * t278 + qJDD(3) * mrSges(6,3) + qJD(3) * t317;
t304 = t320 * mrSges(6,1) - t321 * mrSges(6,3);
t370 = -t320 * mrSges(5,1) - t321 * mrSges(5,2) - t304;
t272 = m(5) * t283 - qJDD(3) * mrSges(5,2) - qJD(3) * t316 + t378 * t311 + t370 * t320 + t366;
t358 = t374 * t286 - t348 * t287;
t282 = t321 * t379 + t358;
t312 = t374 * t335 + t348 * t336;
t315 = -qJD(3) * mrSges(5,2) - t320 * mrSges(5,3);
t279 = -qJDD(3) * pkin(4) - t355 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t303) * t321 - t358;
t318 = -t320 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t361 = -m(6) * t279 + qJDD(3) * mrSges(6,1) + qJD(3) * t318;
t273 = m(5) * t282 + qJDD(3) * mrSges(5,1) + qJD(3) * t315 + t378 * t312 + t370 * t321 + t361;
t268 = t348 * t272 + t374 * t273;
t373 = t380 * qJD(3) + t375 * t320 - t376 * t321;
t372 = t375 * qJD(3) + t381 * t320 + t377 * t321;
t371 = t376 * qJD(3) - t377 * t320 + t382 * t321;
t333 = (-t353 * mrSges(4,1) + t351 * mrSges(4,2)) * qJD(1);
t339 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t369;
t264 = m(4) * t289 + qJDD(3) * mrSges(4,1) - t335 * mrSges(4,3) + qJD(3) * t339 - t333 * t368 + t268;
t338 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t368;
t362 = t374 * t272 - t348 * t273;
t265 = m(4) * t290 - qJDD(3) * mrSges(4,2) + t336 * mrSges(4,3) - qJD(3) * t338 + t333 * t369 + t362;
t363 = -t351 * t264 + t353 * t265;
t309 = t350 * t332 - t349 * t334;
t359 = -qJDD(1) * pkin(2) - t309;
t288 = -t336 * pkin(3) + qJDD(4) + t337 * t368 + (-qJ(4) * t346 - pkin(6)) * t356 + t359;
t281 = -0.2e1 * qJD(5) * t321 + (qJD(3) * t320 - t312) * qJ(5) + (qJD(3) * t321 + t311) * pkin(4) + t288;
t276 = m(6) * t281 + t311 * mrSges(6,1) - t312 * mrSges(6,3) - t321 * t317 + t320 * t318;
t274 = m(5) * t288 + t311 * mrSges(5,1) + t312 * mrSges(5,2) + t320 * t315 + t321 * t316 + t276;
t298 = -t356 * pkin(6) + t359;
t357 = -m(4) * t298 + t336 * mrSges(4,1) - t335 * mrSges(4,2) - t338 * t368 + t339 * t369 - t274;
t327 = Ifges(4,5) * qJD(3) + (t351 * Ifges(4,1) + t353 * Ifges(4,4)) * qJD(1);
t326 = Ifges(4,6) * qJD(3) + (t351 * Ifges(4,4) + t353 * Ifges(4,2)) * qJD(1);
t275 = t312 * mrSges(6,2) + t321 * t304 - t361;
t267 = mrSges(5,2) * t288 + mrSges(6,2) * t279 - mrSges(5,3) * t282 - mrSges(6,3) * t281 - qJ(5) * t276 - t372 * qJD(3) + t376 * qJDD(3) - t377 * t311 + t382 * t312 + t373 * t320;
t266 = -mrSges(5,1) * t288 - mrSges(6,1) * t281 + mrSges(6,2) * t278 + mrSges(5,3) * t283 - pkin(4) * t276 + t371 * qJD(3) + t375 * qJDD(3) + t381 * t311 + t377 * t312 + t373 * t321;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t364 - mrSges(2,2) * t360 + Ifges(3,3) * qJDD(1) + mrSges(3,1) * t309 - mrSges(3,2) * t310 + t351 * (mrSges(4,2) * t298 - mrSges(4,3) * t289 + Ifges(4,1) * t335 + Ifges(4,4) * t336 + Ifges(4,5) * qJDD(3) - qJ(4) * t268 - qJD(3) * t326 - t348 * t266 + t374 * t267) + t353 * (-mrSges(4,1) * t298 + mrSges(4,3) * t290 + Ifges(4,4) * t335 + Ifges(4,2) * t336 + Ifges(4,6) * qJDD(3) - pkin(3) * t274 + qJ(4) * t362 + qJD(3) * t327 + t374 * t266 + t348 * t267) + pkin(2) * t357 + pkin(6) * t363 + pkin(1) * (t349 * (m(3) * t310 - t356 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t363) + t350 * (m(3) * t309 + qJDD(1) * mrSges(3,1) - t356 * mrSges(3,2) + t357)); m(3) * t347 + t353 * t264 + t351 * t265; Ifges(4,5) * t335 + Ifges(4,6) * t336 + mrSges(4,1) * t289 - mrSges(4,2) * t290 + mrSges(5,1) * t282 - mrSges(5,2) * t283 - mrSges(6,1) * t279 + mrSges(6,3) * t278 - pkin(4) * t275 + qJ(5) * t366 + pkin(3) * t268 + t372 * t321 + (-qJ(5) * t304 + t371) * t320 + t376 * t312 + (-qJ(5) * mrSges(6,2) - t375) * t311 + (t351 * t326 - t353 * t327) * qJD(1) + (Ifges(4,3) - t380) * qJDD(3); t274; t275;];
tauJ = t1;
