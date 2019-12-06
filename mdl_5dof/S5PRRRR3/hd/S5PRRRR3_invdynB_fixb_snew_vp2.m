% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5PRRRR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5PRRRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:13
% EndTime: 2019-12-05 17:06:15
% DurationCPUTime: 1.94s
% Computational Cost: add. (31009->174), mult. (39113->224), div. (0->0), fcn. (25560->10), ass. (0->80)
t349 = qJD(2) + qJD(3);
t343 = qJD(4) + t349;
t353 = sin(qJ(5));
t374 = t343 * t353;
t357 = cos(qJ(5));
t373 = t343 * t357;
t351 = sin(pkin(9));
t352 = cos(pkin(9));
t339 = t351 * g(1) - t352 * g(2);
t340 = -t352 * g(1) - t351 * g(2);
t356 = sin(qJ(2));
t360 = cos(qJ(2));
t324 = t360 * t339 - t356 * t340;
t322 = qJDD(2) * pkin(2) + t324;
t325 = t356 * t339 + t360 * t340;
t361 = qJD(2) ^ 2;
t323 = -t361 * pkin(2) + t325;
t355 = sin(qJ(3));
t359 = cos(qJ(3));
t317 = t359 * t322 - t355 * t323;
t348 = qJDD(2) + qJDD(3);
t315 = t348 * pkin(3) + t317;
t318 = t355 * t322 + t359 * t323;
t347 = t349 ^ 2;
t316 = -t347 * pkin(3) + t318;
t354 = sin(qJ(4));
t358 = cos(qJ(4));
t312 = t354 * t315 + t358 * t316;
t341 = t343 ^ 2;
t342 = qJDD(4) + t348;
t310 = -t341 * pkin(4) + t342 * pkin(8) + t312;
t350 = -g(3) + qJDD(1);
t307 = -t353 * t310 + t357 * t350;
t331 = (-mrSges(6,1) * t357 + mrSges(6,2) * t353) * t343;
t371 = qJD(5) * t343;
t332 = t353 * t342 + t357 * t371;
t335 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t373;
t305 = m(6) * t307 + qJDD(5) * mrSges(6,1) - t332 * mrSges(6,3) + qJD(5) * t335 - t331 * t374;
t308 = t357 * t310 + t353 * t350;
t333 = t357 * t342 - t353 * t371;
t334 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t374;
t306 = m(6) * t308 - qJDD(5) * mrSges(6,2) + t333 * mrSges(6,3) - qJD(5) * t334 + t331 * t373;
t365 = -t353 * t305 + t357 * t306;
t296 = m(5) * t312 - t341 * mrSges(5,1) - t342 * mrSges(5,2) + t365;
t311 = t358 * t315 - t354 * t316;
t309 = -t342 * pkin(4) - t341 * pkin(8) - t311;
t362 = -m(6) * t309 + t333 * mrSges(6,1) - t332 * mrSges(6,2) - t334 * t374 + t335 * t373;
t301 = m(5) * t311 + t342 * mrSges(5,1) - t341 * mrSges(5,2) + t362;
t293 = t354 * t296 + t358 * t301;
t290 = m(4) * t317 + t348 * mrSges(4,1) - t347 * mrSges(4,2) + t293;
t366 = t358 * t296 - t354 * t301;
t291 = m(4) * t318 - t347 * mrSges(4,1) - t348 * mrSges(4,2) + t366;
t285 = t359 * t290 + t355 * t291;
t283 = m(3) * t324 + qJDD(2) * mrSges(3,1) - t361 * mrSges(3,2) + t285;
t367 = -t355 * t290 + t359 * t291;
t284 = m(3) * t325 - t361 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t367;
t277 = t360 * t283 + t356 * t284;
t275 = m(2) * t339 + t277;
t368 = -t356 * t283 + t360 * t284;
t276 = m(2) * t340 + t368;
t372 = t352 * t275 + t351 * t276;
t297 = t357 * t305 + t353 * t306;
t370 = m(5) * t350 + t297;
t369 = -t351 * t275 + t352 * t276;
t364 = m(4) * t350 + t370;
t363 = m(3) * t350 + t364;
t328 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t353 + Ifges(6,4) * t357) * t343;
t327 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t353 + Ifges(6,2) * t357) * t343;
t326 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t353 + Ifges(6,6) * t357) * t343;
t299 = mrSges(6,2) * t309 - mrSges(6,3) * t307 + Ifges(6,1) * t332 + Ifges(6,4) * t333 + Ifges(6,5) * qJDD(5) - qJD(5) * t327 + t326 * t373;
t298 = -mrSges(6,1) * t309 + mrSges(6,3) * t308 + Ifges(6,4) * t332 + Ifges(6,2) * t333 + Ifges(6,6) * qJDD(5) + qJD(5) * t328 - t326 * t374;
t292 = -mrSges(5,1) * t350 - mrSges(6,1) * t307 + mrSges(6,2) * t308 + mrSges(5,3) * t312 + t341 * Ifges(5,5) - Ifges(6,5) * t332 + Ifges(5,6) * t342 - Ifges(6,6) * t333 - Ifges(6,3) * qJDD(5) - pkin(4) * t297 + (-t327 * t353 + t328 * t357) * t343;
t286 = mrSges(5,2) * t350 - mrSges(5,3) * t311 + Ifges(5,5) * t342 - t341 * Ifges(5,6) - pkin(8) * t297 - t353 * t298 + t357 * t299;
t279 = mrSges(4,2) * t350 - mrSges(4,3) * t317 + Ifges(4,5) * t348 - t347 * Ifges(4,6) - pkin(7) * t293 + t358 * t286 - t354 * t292;
t278 = -mrSges(4,1) * t350 + mrSges(4,3) * t318 + t347 * Ifges(4,5) + Ifges(4,6) * t348 - pkin(3) * t370 + pkin(7) * t366 + t354 * t286 + t358 * t292;
t271 = mrSges(3,2) * t350 - mrSges(3,3) * t324 + Ifges(3,5) * qJDD(2) - t361 * Ifges(3,6) - pkin(6) * t285 - t355 * t278 + t359 * t279;
t270 = -mrSges(3,1) * t350 + mrSges(3,3) * t325 + t361 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t364 + pkin(6) * t367 + t359 * t278 + t355 * t279;
t269 = mrSges(2,2) * t350 - mrSges(2,3) * t339 - pkin(5) * t277 - t356 * t270 + t360 * t271;
t268 = -mrSges(2,1) * t350 + mrSges(2,3) * t340 - pkin(1) * t363 + pkin(5) * t368 + t360 * t270 + t356 * t271;
t1 = [-m(1) * g(1) + t369; -m(1) * g(2) + t372; -m(1) * g(3) + m(2) * t350 + t363; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t372 - t351 * t268 + t352 * t269; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t369 + t352 * t268 + t351 * t269; pkin(1) * t277 + mrSges(2,1) * t339 - mrSges(2,2) * t340 + pkin(2) * t285 - mrSges(3,2) * t325 + mrSges(3,1) * t324 + pkin(3) * t293 + mrSges(4,1) * t317 - mrSges(4,2) * t318 + pkin(4) * t362 + pkin(8) * t365 + mrSges(5,1) * t311 - mrSges(5,2) * t312 + t353 * t299 + t357 * t298 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(5,3) * t342 + Ifges(4,3) * t348 + Ifges(3,3) * qJDD(2);];
tauB = t1;
