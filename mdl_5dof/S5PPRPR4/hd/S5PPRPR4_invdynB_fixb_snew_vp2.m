% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PPRPR4
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-31 17:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PPRPR4_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR4_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:32:21
% EndTime: 2019-12-31 17:32:22
% DurationCPUTime: 1.03s
% Computational Cost: add. (9295->178), mult. (18823->223), div. (0->0), fcn. (11933->8), ass. (0->82)
t353 = qJD(3) ^ 2;
t347 = cos(pkin(8));
t378 = pkin(4) * t347;
t377 = -mrSges(2,2) + mrSges(3,3);
t345 = sin(pkin(8));
t376 = mrSges(5,2) * t345;
t343 = t347 ^ 2;
t375 = t343 * t353;
t346 = sin(pkin(7));
t348 = cos(pkin(7));
t332 = t346 * g(1) - t348 * g(2);
t330 = qJDD(2) - t332;
t333 = -t348 * g(1) - t346 * g(2);
t350 = sin(qJ(3));
t352 = cos(qJ(3));
t320 = t350 * t330 + t352 * t333;
t316 = -t353 * pkin(3) + qJDD(3) * qJ(4) + t320;
t344 = g(3) - qJDD(1);
t371 = qJD(3) * qJD(4);
t373 = t347 * t344 - 0.2e1 * t345 * t371;
t303 = (-pkin(6) * qJDD(3) + t353 * t378 - t316) * t345 + t373;
t307 = t345 * t344 + (t316 + 0.2e1 * t371) * t347;
t370 = qJDD(3) * t347;
t304 = -pkin(4) * t375 + pkin(6) * t370 + t307;
t349 = sin(qJ(5));
t351 = cos(qJ(5));
t301 = t351 * t303 - t349 * t304;
t359 = -t345 * t349 + t347 * t351;
t323 = t359 * qJD(3);
t360 = t345 * t351 + t347 * t349;
t324 = t360 * qJD(3);
t314 = -t323 * mrSges(6,1) + t324 * mrSges(6,2);
t318 = t323 * qJD(5) + t360 * qJDD(3);
t321 = -qJD(5) * mrSges(6,2) + t323 * mrSges(6,3);
t299 = m(6) * t301 + qJDD(5) * mrSges(6,1) - t318 * mrSges(6,3) + qJD(5) * t321 - t324 * t314;
t302 = t349 * t303 + t351 * t304;
t317 = -t324 * qJD(5) + t359 * qJDD(3);
t322 = qJD(5) * mrSges(6,1) - t324 * mrSges(6,3);
t300 = m(6) * t302 - qJDD(5) * mrSges(6,2) + t317 * mrSges(6,3) - qJD(5) * t322 + t323 * t314;
t292 = t351 * t299 + t349 * t300;
t306 = -t345 * t316 + t373;
t358 = mrSges(5,3) * qJDD(3) + t353 * (-mrSges(5,1) * t347 + t376);
t290 = m(5) * t306 - t358 * t345 + t292;
t366 = -t349 * t299 + t351 * t300;
t291 = m(5) * t307 + t358 * t347 + t366;
t367 = -t345 * t290 + t347 * t291;
t286 = m(4) * t320 - t353 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t367;
t319 = t352 * t330 - t350 * t333;
t361 = qJDD(4) - t319;
t313 = -qJDD(3) * pkin(3) - t353 * qJ(4) + t361;
t342 = t345 ^ 2;
t305 = (-pkin(3) - t378) * qJDD(3) + (-qJ(4) + (-t342 - t343) * pkin(6)) * t353 + t361;
t355 = m(6) * t305 - t317 * mrSges(6,1) + t318 * mrSges(6,2) - t323 * t321 + t324 * t322;
t354 = -m(5) * t313 + mrSges(5,1) * t370 - t355 + (t342 * t353 + t375) * mrSges(5,3);
t295 = m(4) * t319 - t353 * mrSges(4,2) + (mrSges(4,1) - t376) * qJDD(3) + t354;
t283 = t350 * t286 + t352 * t295;
t356 = -m(3) * t330 - t283;
t281 = m(2) * t332 + t356;
t368 = t352 * t286 - t350 * t295;
t365 = m(3) * t333 + t368;
t282 = m(2) * t333 + t365;
t374 = t348 * t281 + t346 * t282;
t362 = Ifges(5,5) * t345 + Ifges(5,6) * t347;
t372 = t353 * t362;
t369 = -t346 * t281 + t348 * t282;
t364 = Ifges(5,1) * t345 + Ifges(5,4) * t347;
t363 = Ifges(5,4) * t345 + Ifges(5,2) * t347;
t288 = t347 * t290 + t345 * t291;
t357 = -m(3) * t344 - t288;
t310 = Ifges(6,1) * t324 + Ifges(6,4) * t323 + Ifges(6,5) * qJD(5);
t309 = Ifges(6,4) * t324 + Ifges(6,2) * t323 + Ifges(6,6) * qJD(5);
t308 = Ifges(6,5) * t324 + Ifges(6,6) * t323 + Ifges(6,3) * qJD(5);
t294 = mrSges(6,2) * t305 - mrSges(6,3) * t301 + Ifges(6,1) * t318 + Ifges(6,4) * t317 + Ifges(6,5) * qJDD(5) - qJD(5) * t309 + t323 * t308;
t293 = -mrSges(6,1) * t305 + mrSges(6,3) * t302 + Ifges(6,4) * t318 + Ifges(6,2) * t317 + Ifges(6,6) * qJDD(5) + qJD(5) * t310 - t324 * t308;
t287 = -m(4) * t344 + t357;
t284 = mrSges(5,2) * t313 - mrSges(5,3) * t306 - pkin(6) * t292 + t364 * qJDD(3) - t349 * t293 + t351 * t294 + t347 * t372;
t277 = -mrSges(5,1) * t313 + mrSges(5,3) * t307 - pkin(4) * t355 + pkin(6) * t366 + t363 * qJDD(3) + t351 * t293 + t349 * t294 - t345 * t372;
t276 = -mrSges(4,1) * t344 - mrSges(5,1) * t306 - mrSges(6,1) * t301 + mrSges(5,2) * t307 + mrSges(6,2) * t302 + mrSges(4,3) * t320 - Ifges(6,5) * t318 - Ifges(6,6) * t317 - Ifges(6,3) * qJDD(5) - pkin(3) * t288 - pkin(4) * t292 - t324 * t309 + t323 * t310 + (Ifges(4,6) - t362) * qJDD(3) + (-t345 * t363 + t347 * t364 + Ifges(4,5)) * t353;
t275 = mrSges(4,2) * t344 - mrSges(4,3) * t319 + Ifges(4,5) * qJDD(3) - t353 * Ifges(4,6) - qJ(4) * t288 - t345 * t277 + t347 * t284;
t274 = mrSges(3,2) * t330 - mrSges(2,3) * t332 - pkin(5) * t283 - qJ(2) * t287 + t352 * t275 - t350 * t276 + t377 * t344;
t273 = -t350 * t275 - t352 * t276 + pkin(2) * t288 - pkin(5) * t368 - pkin(1) * t287 + (pkin(2) * m(4) + mrSges(2,1) + mrSges(3,1)) * t344 + (mrSges(3,2) + mrSges(2,3)) * t333;
t1 = [-m(1) * g(1) + t369; -m(1) * g(2) + t374; -m(1) * g(3) + (-m(2) - m(4)) * t344 + t357; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t374 - t346 * t273 + t348 * t274; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t369 + t348 * t273 + t346 * t274; pkin(1) * t356 + qJ(2) * t365 + mrSges(2,1) * t332 - pkin(2) * t283 - mrSges(3,1) * t330 - t345 * t284 - t347 * t277 - pkin(3) * (-qJDD(3) * t376 + t354) - qJ(4) * t367 - mrSges(4,1) * t319 + mrSges(4,2) * t320 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - Ifges(4,3) * qJDD(3) + t377 * t333;];
tauB = t1;
