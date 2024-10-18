% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRRRR13
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
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
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRRRR13_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR13_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:30:49
% EndTime: 2024-09-27 17:30:50
% DurationCPUTime: 0.70s
% Computational Cost: add. (23657->175), mult. (27699->243), div. (0->0), fcn. (17587->12), ass. (0->93)
t359 = sin(pkin(5));
t393 = pkin(9) * t359;
t360 = cos(pkin(5));
t392 = t360 * g(3);
t357 = qJD(1) + qJD(2);
t352 = qJD(3) + t357;
t350 = t352 ^ 2;
t391 = t350 * t359 ^ 2;
t390 = t352 * t359;
t362 = sin(qJ(4));
t389 = t359 * t362;
t367 = cos(qJ(4));
t388 = t359 * t367;
t387 = t360 * t362;
t386 = t360 * t367;
t365 = sin(qJ(1));
t370 = cos(qJ(1));
t380 = t365 * g(1) - t370 * g(2);
t341 = qJDD(1) * pkin(1) + t380;
t376 = -t370 * g(1) - t365 * g(2);
t342 = -qJD(1) ^ 2 * pkin(1) + t376;
t364 = sin(qJ(2));
t369 = cos(qJ(2));
t326 = t369 * t341 - t364 * t342;
t355 = qJDD(1) + qJDD(2);
t323 = t355 * pkin(2) + t326;
t327 = t364 * t341 + t369 * t342;
t354 = t357 ^ 2;
t324 = -t354 * pkin(2) + t327;
t363 = sin(qJ(3));
t368 = cos(qJ(3));
t308 = t368 * t323 - t363 * t324;
t351 = qJDD(3) + t355;
t384 = qJD(4) * t352;
t334 = (t351 * t362 + t367 * t384) * t359;
t344 = t360 * t351 + qJDD(4);
t345 = t360 * t352 + qJD(4);
t300 = t351 * pkin(3) + t350 * t393 + t308;
t309 = t363 * t323 + t368 * t324;
t301 = -t350 * pkin(3) + t351 * t393 + t309;
t377 = t300 * t386 - t362 * t301;
t289 = t344 * pkin(4) - t334 * pkin(10) + (pkin(4) * t362 * t391 + (pkin(10) * t345 * t352 - g(3)) * t359) * t367 + t377;
t292 = -g(3) * t389 + t300 * t387 + t367 * t301;
t382 = t352 * t389;
t333 = t345 * pkin(4) - pkin(10) * t382;
t335 = (t351 * t367 - t362 * t384) * t359;
t383 = t367 ^ 2 * t391;
t290 = -pkin(4) * t383 + t335 * pkin(10) - t345 * t333 + t292;
t361 = sin(qJ(5));
t366 = cos(qJ(5));
t287 = t366 * t289 - t361 * t290;
t328 = (-t361 * t362 + t366 * t367) * t390;
t306 = t328 * qJD(5) + t366 * t334 + t361 * t335;
t329 = (t361 * t367 + t362 * t366) * t390;
t314 = -t328 * mrSges(6,1) + t329 * mrSges(6,2);
t343 = qJD(5) + t345;
t315 = -t343 * mrSges(6,2) + t328 * mrSges(6,3);
t340 = qJDD(5) + t344;
t283 = m(6) * t287 + t340 * mrSges(6,1) - t306 * mrSges(6,3) - t329 * t314 + t343 * t315;
t288 = t361 * t289 + t366 * t290;
t305 = -t329 * qJD(5) - t361 * t334 + t366 * t335;
t316 = t343 * mrSges(6,1) - t329 * mrSges(6,3);
t284 = m(6) * t288 - t340 * mrSges(6,2) + t305 * mrSges(6,3) + t328 * t314 - t343 * t316;
t277 = t366 * t283 + t361 * t284;
t291 = -g(3) * t388 + t377;
t381 = t352 * t388;
t331 = -t345 * mrSges(5,2) + mrSges(5,3) * t381;
t332 = (-mrSges(5,1) * t367 + mrSges(5,2) * t362) * t390;
t275 = m(5) * t291 + t344 * mrSges(5,1) - t334 * mrSges(5,3) + t345 * t331 - t332 * t382 + t277;
t330 = t345 * mrSges(5,1) - mrSges(5,3) * t382;
t378 = -t361 * t283 + t366 * t284;
t276 = m(5) * t292 - t344 * mrSges(5,2) + t335 * mrSges(5,3) - t345 * t330 + t332 * t381 + t378;
t295 = -t359 * t300 - t392;
t294 = -pkin(10) * t383 - t335 * pkin(4) - t392 + (t333 * t352 * t362 - t300) * t359;
t374 = m(6) * t294 - t305 * mrSges(6,1) + t306 * mrSges(6,2) - t328 * t315 + t329 * t316;
t375 = t275 * t386 + t276 * t387 - t359 * (m(5) * t295 - t335 * mrSges(5,1) + t334 * mrSges(5,2) + (t330 * t362 - t331 * t367) * t390 + t374);
t264 = m(4) * t308 + t351 * mrSges(4,1) - t350 * mrSges(4,2) + t375;
t379 = -t362 * t275 + t367 * t276;
t268 = m(4) * t309 - t350 * mrSges(4,1) - t351 * mrSges(4,2) + t379;
t385 = t368 * t264 + t363 * t268;
t321 = Ifges(5,6) * t345 + (Ifges(5,4) * t362 + Ifges(5,2) * t367) * t390;
t322 = Ifges(5,5) * t345 + (Ifges(5,1) * t362 + Ifges(5,4) * t367) * t390;
t311 = Ifges(6,4) * t329 + Ifges(6,2) * t328 + Ifges(6,6) * t343;
t312 = Ifges(6,1) * t329 + Ifges(6,4) * t328 + Ifges(6,5) * t343;
t372 = mrSges(6,1) * t287 - mrSges(6,2) * t288 + Ifges(6,5) * t306 + Ifges(6,6) * t305 + Ifges(6,3) * t340 + t329 * t311 - t328 * t312;
t271 = mrSges(5,1) * t291 - mrSges(5,2) * t292 + Ifges(5,5) * t334 + Ifges(5,6) * t335 + Ifges(5,3) * t344 + pkin(4) * t277 + (t321 * t362 - t322 * t367) * t390 + t372;
t310 = Ifges(6,5) * t329 + Ifges(6,6) * t328 + Ifges(6,3) * t343;
t278 = -mrSges(6,1) * t294 + mrSges(6,3) * t288 + Ifges(6,4) * t306 + Ifges(6,2) * t305 + Ifges(6,6) * t340 - t329 * t310 + t343 * t312;
t279 = mrSges(6,2) * t294 - mrSges(6,3) * t287 + Ifges(6,1) * t306 + Ifges(6,4) * t305 + Ifges(6,5) * t340 + t328 * t310 - t343 * t311;
t320 = Ifges(5,3) * t345 + (Ifges(5,5) * t362 + Ifges(5,6) * t367) * t390;
t373 = -mrSges(4,2) * t309 + (-mrSges(5,1) * t295 + mrSges(5,3) * t292 + Ifges(5,4) * t334 + Ifges(5,2) * t335 + Ifges(5,6) * t344 - pkin(4) * t374 + pkin(10) * t378 + t366 * t278 + t361 * t279 - t320 * t382 + t345 * t322) * t388 + pkin(3) * t375 + (mrSges(5,2) * t295 - mrSges(5,3) * t291 + Ifges(5,1) * t334 + Ifges(5,4) * t335 + Ifges(5,5) * t344 - pkin(10) * t277 - t361 * t278 + t366 * t279 + t320 * t381 - t345 * t321) * t389 + t379 * t393 + t360 * t271 + mrSges(4,1) * t308 + Ifges(4,3) * t351;
t371 = mrSges(3,1) * t326 - mrSges(3,2) * t327 + Ifges(3,3) * t355 + pkin(2) * t385 + t373;
t1 = [Ifges(2,3) * qJDD(1) + pkin(1) * (t364 * (m(3) * t327 - t354 * mrSges(3,1) - t355 * mrSges(3,2) - t363 * t264 + t368 * t268) + t369 * (m(3) * t326 + t355 * mrSges(3,1) - t354 * mrSges(3,2) + t385)) + mrSges(2,1) * t380 - mrSges(2,2) * t376 + t371; t371; t373; t271; t372;];
tauJ = t1;
