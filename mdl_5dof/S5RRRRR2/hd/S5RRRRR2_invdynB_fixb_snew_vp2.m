% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRRR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:08
% EndTime: 2019-12-05 18:53:11
% DurationCPUTime: 1.66s
% Computational Cost: add. (19192->237), mult. (26009->306), div. (0->0), fcn. (17830->10), ass. (0->90)
t365 = qJD(1) + qJD(2);
t367 = sin(qJ(4));
t368 = sin(qJ(3));
t372 = cos(qJ(4));
t373 = cos(qJ(3));
t349 = (t367 * t368 - t372 * t373) * t365;
t387 = t365 * t368;
t386 = t365 * t373;
t370 = sin(qJ(1));
t375 = cos(qJ(1));
t359 = -t375 * g(1) - t370 * g(2);
t376 = qJD(1) ^ 2;
t355 = -t376 * pkin(1) + t359;
t369 = sin(qJ(2));
t374 = cos(qJ(2));
t358 = t370 * g(1) - t375 * g(2);
t379 = qJDD(1) * pkin(1) + t358;
t340 = t374 * t355 + t369 * t379;
t336 = -t373 * g(3) - t368 * t340;
t351 = (-mrSges(4,1) * t373 + mrSges(4,2) * t368) * t365;
t363 = qJDD(1) + qJDD(2);
t382 = qJD(3) * t365;
t352 = t368 * t363 + t373 * t382;
t357 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t386;
t337 = -t368 * g(3) + t373 * t340;
t361 = t365 ^ 2;
t331 = (-t361 * t373 ^ 2 - qJD(3) ^ 2) * pkin(2) + t337;
t378 = (t361 * t368 * t373 + qJDD(3)) * pkin(2) + t336;
t317 = t372 * t331 + t367 * t378;
t339 = t369 * t355 - t374 * t379;
t381 = t368 * t382;
t353 = t373 * t363 - t381;
t327 = (-t353 + t381) * pkin(2) + t339;
t366 = sin(qJ(5));
t371 = cos(qJ(5));
t314 = -t366 * t317 + t371 * t327;
t326 = -t349 * qJD(4) + t372 * t352 + t367 * t353;
t350 = (t367 * t373 + t368 * t372) * t365;
t364 = qJD(3) + qJD(4);
t341 = -t366 * t350 + t371 * t364;
t362 = qJDD(3) + qJDD(4);
t319 = t341 * qJD(5) + t371 * t326 + t366 * t362;
t342 = t371 * t350 + t366 * t364;
t323 = -t341 * mrSges(6,1) + t342 * mrSges(6,2);
t325 = -t350 * qJD(4) - t367 * t352 + t372 * t353;
t324 = qJDD(5) - t325;
t345 = qJD(5) + t349;
t329 = -t345 * mrSges(6,2) + t341 * mrSges(6,3);
t312 = m(6) * t314 + t324 * mrSges(6,1) - t319 * mrSges(6,3) - t342 * t323 + t345 * t329;
t315 = t371 * t317 + t366 * t327;
t318 = -t342 * qJD(5) - t366 * t326 + t371 * t362;
t330 = t345 * mrSges(6,1) - t342 * mrSges(6,3);
t313 = m(6) * t315 - t324 * mrSges(6,2) + t318 * mrSges(6,3) + t341 * t323 - t345 * t330;
t335 = t349 * mrSges(5,1) + t350 * mrSges(5,2);
t344 = t364 * mrSges(5,1) - t350 * mrSges(5,3);
t306 = m(5) * t317 - t362 * mrSges(5,2) + t325 * mrSges(5,3) - t366 * t312 + t371 * t313 - t349 * t335 - t364 * t344;
t316 = t367 * t331 - t372 * t378;
t343 = -t364 * mrSges(5,2) - t349 * mrSges(5,3);
t311 = t362 * mrSges(5,1) + t318 * mrSges(6,1) - t319 * mrSges(6,2) - t326 * mrSges(5,3) + t341 * t329 - t342 * t330 - t350 * t335 + t364 * t343 + (-m(5) - m(6)) * t316;
t383 = t367 * t306 + t372 * t311;
t300 = m(4) * t336 + qJDD(3) * mrSges(4,1) - t352 * mrSges(4,3) + qJD(3) * t357 - t351 * t387 + t383;
t356 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t387;
t301 = m(4) * t337 - qJDD(3) * mrSges(4,2) + t353 * mrSges(4,3) - qJD(3) * t356 + t372 * t306 - t367 * t311 + t351 * t386;
t295 = m(3) * t340 - t361 * mrSges(3,1) - t363 * mrSges(3,2) - t368 * t300 + t373 * t301;
t377 = m(5) * t327 - t325 * mrSges(5,1) + t326 * mrSges(5,2) + t371 * t312 + t366 * t313 + t349 * t343 + t350 * t344;
t304 = t363 * mrSges(3,1) + t353 * mrSges(4,1) - t361 * mrSges(3,2) - t352 * mrSges(4,2) + (-t356 * t368 + t357 * t373) * t365 + (-m(3) - m(4)) * t339 - t377;
t385 = t369 * t295 + t374 * t304;
t384 = t373 * t300 + t368 * t301;
t348 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t368 + Ifges(4,4) * t373) * t365;
t347 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t368 + Ifges(4,2) * t373) * t365;
t346 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t368 + Ifges(4,6) * t373) * t365;
t334 = Ifges(5,1) * t350 - Ifges(5,4) * t349 + Ifges(5,5) * t364;
t333 = Ifges(5,4) * t350 - Ifges(5,2) * t349 + Ifges(5,6) * t364;
t332 = Ifges(5,5) * t350 - Ifges(5,6) * t349 + Ifges(5,3) * t364;
t322 = Ifges(6,1) * t342 + Ifges(6,4) * t341 + Ifges(6,5) * t345;
t321 = Ifges(6,4) * t342 + Ifges(6,2) * t341 + Ifges(6,6) * t345;
t320 = Ifges(6,5) * t342 + Ifges(6,6) * t341 + Ifges(6,3) * t345;
t309 = mrSges(6,2) * t316 - mrSges(6,3) * t314 + Ifges(6,1) * t319 + Ifges(6,4) * t318 + Ifges(6,5) * t324 + t341 * t320 - t345 * t321;
t308 = -mrSges(6,1) * t316 + mrSges(6,3) * t315 + Ifges(6,4) * t319 + Ifges(6,2) * t318 + Ifges(6,6) * t324 - t342 * t320 + t345 * t322;
t307 = -mrSges(5,1) * t327 - mrSges(6,1) * t314 + mrSges(6,2) * t315 + mrSges(5,3) * t317 + Ifges(5,4) * t326 - Ifges(6,5) * t319 + Ifges(5,2) * t325 + Ifges(5,6) * t362 - Ifges(6,6) * t318 - Ifges(6,3) * t324 - t342 * t321 + t341 * t322 - t350 * t332 + t364 * t334;
t302 = mrSges(5,2) * t327 + mrSges(5,3) * t316 + Ifges(5,1) * t326 + Ifges(5,4) * t325 + Ifges(5,5) * t362 - t366 * t308 + t371 * t309 - t349 * t332 - t364 * t333;
t297 = mrSges(4,2) * t339 - mrSges(4,3) * t336 + Ifges(4,1) * t352 + Ifges(4,4) * t353 + Ifges(4,5) * qJDD(3) - qJD(3) * t347 + t372 * t302 - t367 * t307 + t346 * t386;
t296 = Ifges(3,6) * t363 + t361 * Ifges(3,5) + mrSges(3,1) * g(3) + mrSges(3,3) * t340 - Ifges(4,5) * t352 - Ifges(4,6) * t353 - Ifges(4,3) * qJDD(3) - mrSges(4,1) * t336 + mrSges(4,2) * t337 - Ifges(5,5) * t326 - Ifges(5,6) * t325 - Ifges(5,3) * t362 - t350 * t333 - t349 * t334 + mrSges(5,1) * t316 + mrSges(5,2) * t317 - t366 * t309 - t371 * t308 - pkin(2) * t383 + (-t368 * t347 + t373 * t348) * t365;
t293 = -mrSges(4,1) * t339 + mrSges(4,3) * t337 + Ifges(4,4) * t352 + Ifges(4,2) * t353 + Ifges(4,6) * qJDD(3) - pkin(2) * t377 + qJD(3) * t348 + t367 * t302 + t372 * t307 - t346 * t387;
t292 = m(2) * t359 - t376 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t374 * t295 - t369 * t304;
t291 = m(2) * t358 + qJDD(1) * mrSges(2,1) - t376 * mrSges(2,2) + t385;
t290 = -mrSges(3,2) * g(3) + mrSges(3,3) * t339 + Ifges(3,5) * t363 - t361 * Ifges(3,6) - t368 * t293 + t373 * t297;
t289 = -mrSges(2,2) * g(3) - mrSges(2,3) * t358 + Ifges(2,5) * qJDD(1) - t376 * Ifges(2,6) + t374 * t290 - t369 * t296;
t288 = Ifges(2,6) * qJDD(1) + t376 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t359 + t369 * t290 + t374 * t296 - pkin(1) * (-m(3) * g(3) + t384);
t1 = [-m(1) * g(1) - t370 * t291 + t375 * t292; -m(1) * g(2) + t375 * t291 + t370 * t292; (-m(1) - m(2) - m(3)) * g(3) + t384; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - t370 * t288 + t375 * t289; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t375 * t288 + t370 * t289; -mrSges(1,1) * g(2) + mrSges(2,1) * t358 - mrSges(3,1) * t339 + mrSges(1,2) * g(1) - mrSges(2,2) * t359 - mrSges(3,2) * t340 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t363 + pkin(1) * t385 + t373 * t293 + t368 * t297;];
tauB = t1;
