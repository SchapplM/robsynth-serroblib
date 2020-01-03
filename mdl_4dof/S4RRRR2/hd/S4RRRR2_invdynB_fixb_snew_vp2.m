% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR2_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:13
% EndTime: 2019-12-31 17:23:14
% DurationCPUTime: 1.13s
% Computational Cost: add. (15694->199), mult. (20169->256), div. (0->0), fcn. (11357->8), ass. (0->82)
t356 = qJD(1) + qJD(2);
t352 = t356 ^ 2;
t379 = pkin(3) * t352;
t359 = sin(qJ(3));
t378 = t356 * t359;
t363 = cos(qJ(3));
t377 = t356 * t363;
t361 = sin(qJ(1));
t365 = cos(qJ(1));
t349 = t361 * g(1) - t365 * g(2);
t344 = qJDD(1) * pkin(1) + t349;
t350 = -t365 * g(1) - t361 * g(2);
t366 = qJD(1) ^ 2;
t345 = -t366 * pkin(1) + t350;
t360 = sin(qJ(2));
t364 = cos(qJ(2));
t328 = t360 * t344 + t364 * t345;
t354 = qJDD(1) + qJDD(2);
t326 = -t352 * pkin(2) + t354 * pkin(6) + t328;
t376 = t359 * t326;
t374 = qJD(3) * t356;
t339 = t359 * t354 + t363 * t374;
t311 = qJDD(3) * pkin(3) - t339 * pkin(7) - t376 + (pkin(7) * t374 + t359 * t379 - g(3)) * t363;
t318 = -t359 * g(3) + t363 * t326;
t340 = t363 * t354 - t359 * t374;
t348 = qJD(3) * pkin(3) - pkin(7) * t378;
t357 = t363 ^ 2;
t312 = t340 * pkin(7) - qJD(3) * t348 - t357 * t379 + t318;
t358 = sin(qJ(4));
t362 = cos(qJ(4));
t309 = t362 * t311 - t358 * t312;
t334 = (-t358 * t359 + t362 * t363) * t356;
t316 = t334 * qJD(4) + t362 * t339 + t358 * t340;
t335 = (t358 * t363 + t359 * t362) * t356;
t324 = -t334 * mrSges(5,1) + t335 * mrSges(5,2);
t355 = qJD(3) + qJD(4);
t329 = -t355 * mrSges(5,2) + t334 * mrSges(5,3);
t353 = qJDD(3) + qJDD(4);
t307 = m(5) * t309 + t353 * mrSges(5,1) - t316 * mrSges(5,3) - t335 * t324 + t355 * t329;
t310 = t358 * t311 + t362 * t312;
t315 = -t335 * qJD(4) - t358 * t339 + t362 * t340;
t330 = t355 * mrSges(5,1) - t335 * mrSges(5,3);
t308 = m(5) * t310 - t353 * mrSges(5,2) + t315 * mrSges(5,3) + t334 * t324 - t355 * t330;
t299 = t362 * t307 + t358 * t308;
t317 = -t363 * g(3) - t376;
t338 = (-mrSges(4,1) * t363 + mrSges(4,2) * t359) * t356;
t347 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t377;
t297 = m(4) * t317 + qJDD(3) * mrSges(4,1) - t339 * mrSges(4,3) + qJD(3) * t347 - t338 * t378 + t299;
t346 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t378;
t370 = -t358 * t307 + t362 * t308;
t298 = m(4) * t318 - qJDD(3) * mrSges(4,2) + t340 * mrSges(4,3) - qJD(3) * t346 + t338 * t377 + t370;
t371 = -t359 * t297 + t363 * t298;
t292 = m(3) * t328 - t352 * mrSges(3,1) - t354 * mrSges(3,2) + t371;
t327 = t364 * t344 - t360 * t345;
t369 = -t354 * pkin(2) - t327;
t325 = -t352 * pkin(6) + t369;
t313 = t348 * t378 - t340 * pkin(3) + (-pkin(7) * t357 - pkin(6)) * t352 + t369;
t368 = m(5) * t313 - t315 * mrSges(5,1) + t316 * mrSges(5,2) - t334 * t329 + t335 * t330;
t367 = -m(4) * t325 + t340 * mrSges(4,1) - t339 * mrSges(4,2) - t346 * t378 + t347 * t377 - t368;
t303 = m(3) * t327 + t354 * mrSges(3,1) - t352 * mrSges(3,2) + t367;
t288 = t360 * t292 + t364 * t303;
t286 = m(2) * t349 + qJDD(1) * mrSges(2,1) - t366 * mrSges(2,2) + t288;
t372 = t364 * t292 - t360 * t303;
t287 = m(2) * t350 - t366 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t372;
t375 = t365 * t286 + t361 * t287;
t293 = t363 * t297 + t359 * t298;
t373 = -t361 * t286 + t365 * t287;
t333 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t359 + Ifges(4,4) * t363) * t356;
t332 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t359 + Ifges(4,2) * t363) * t356;
t331 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t359 + Ifges(4,6) * t363) * t356;
t321 = Ifges(5,1) * t335 + Ifges(5,4) * t334 + Ifges(5,5) * t355;
t320 = Ifges(5,4) * t335 + Ifges(5,2) * t334 + Ifges(5,6) * t355;
t319 = Ifges(5,5) * t335 + Ifges(5,6) * t334 + Ifges(5,3) * t355;
t301 = mrSges(5,2) * t313 - mrSges(5,3) * t309 + Ifges(5,1) * t316 + Ifges(5,4) * t315 + Ifges(5,5) * t353 + t334 * t319 - t355 * t320;
t300 = -mrSges(5,1) * t313 + mrSges(5,3) * t310 + Ifges(5,4) * t316 + Ifges(5,2) * t315 + Ifges(5,6) * t353 - t335 * t319 + t355 * t321;
t289 = mrSges(4,2) * t325 - mrSges(4,3) * t317 + Ifges(4,1) * t339 + Ifges(4,4) * t340 + Ifges(4,5) * qJDD(3) - pkin(7) * t299 - qJD(3) * t332 - t358 * t300 + t362 * t301 + t331 * t377;
t282 = -mrSges(4,1) * t325 + mrSges(4,3) * t318 + Ifges(4,4) * t339 + Ifges(4,2) * t340 + Ifges(4,6) * qJDD(3) - pkin(3) * t368 + pkin(7) * t370 + qJD(3) * t333 + t362 * t300 + t358 * t301 - t331 * t378;
t281 = Ifges(3,6) * t354 + t352 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t328 - Ifges(4,5) * t339 - Ifges(4,6) * t340 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t317 + mrSges(4,2) * t318 - Ifges(5,5) * t316 - Ifges(5,6) * t315 - Ifges(5,3) * t353 - t335 * t320 + t334 * t321 - mrSges(5,1) * t309 + mrSges(5,2) * t310 - pkin(3) * t299 - pkin(2) * t293 + (-t359 * t332 + t363 * t333) * t356;
t280 = -mrSges(3,2) * g(3) - mrSges(3,3) * t327 + Ifges(3,5) * t354 - t352 * Ifges(3,6) - pkin(6) * t293 - t359 * t282 + t363 * t289;
t279 = -mrSges(2,2) * g(3) - mrSges(2,3) * t349 + Ifges(2,5) * qJDD(1) - t366 * Ifges(2,6) - pkin(5) * t288 + t364 * t280 - t360 * t281;
t278 = Ifges(2,6) * qJDD(1) + t366 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t350 + t360 * t280 + t364 * t281 - pkin(1) * (-m(3) * g(3) + t293) + pkin(5) * t372;
t1 = [-m(1) * g(1) + t373; -m(1) * g(2) + t375; (-m(1) - m(2) - m(3)) * g(3) + t293; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t375 - t361 * t278 + t365 * t279; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t373 + t365 * t278 + t361 * t279; -mrSges(1,1) * g(2) + mrSges(2,1) * t349 + mrSges(3,1) * t327 + mrSges(1,2) * g(1) - mrSges(2,2) * t350 - mrSges(3,2) * t328 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t354 + pkin(1) * t288 + pkin(2) * t367 + pkin(6) * t371 + t363 * t282 + t359 * t289;];
tauB = t1;
