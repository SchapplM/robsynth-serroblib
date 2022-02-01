% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RRPRR6_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:06
% EndTime: 2022-01-20 11:17:08
% DurationCPUTime: 1.40s
% Computational Cost: add. (12557->186), mult. (17171->250), div. (0->0), fcn. (10524->10), ass. (0->89)
t363 = sin(qJ(1));
t367 = cos(qJ(1));
t379 = t363 * g(1) - t367 * g(2);
t339 = qJDD(1) * pkin(1) + t379;
t375 = -t367 * g(1) - t363 * g(2);
t340 = -qJD(1) ^ 2 * pkin(1) + t375;
t362 = sin(qJ(2));
t366 = cos(qJ(2));
t322 = t362 * t339 + t366 * t340;
t356 = qJD(1) + qJD(2);
t353 = t356 ^ 2;
t354 = qJDD(1) + qJDD(2);
t393 = -t353 * pkin(2) + t354 * qJ(3) + 0.2e1 * qJD(3) * t356 + t322;
t358 = sin(pkin(9));
t359 = cos(pkin(9));
t308 = -t359 * g(3) - t393 * t358;
t392 = mrSges(4,2) * t358;
t391 = mrSges(4,3) * t354;
t390 = t353 * t358 ^ 2;
t389 = t356 * t358;
t388 = t359 * t354;
t387 = t359 * t356;
t309 = -t358 * g(3) + t393 * t359;
t376 = -pkin(3) * t359 - pkin(7) * t358;
t334 = t376 * t356;
t297 = t334 * t387 + t309;
t321 = t366 * t339 - t362 * t340;
t373 = -t353 * qJ(3) + qJDD(3) - t321;
t307 = (-pkin(2) + t376) * t354 + t373;
t365 = cos(qJ(4));
t306 = t365 * t307;
t361 = sin(qJ(4));
t384 = qJD(4) * t356;
t331 = (t354 * t365 - t361 * t384) * t358;
t343 = qJDD(4) - t388;
t344 = qJD(4) - t387;
t289 = t343 * pkin(4) - t331 * pkin(8) + t306 + (-pkin(4) * t365 * t390 - pkin(8) * t344 * t389 - t297) * t361;
t292 = t365 * t297 + t361 * t307;
t381 = t365 * t389;
t329 = t344 * pkin(4) - pkin(8) * t381;
t330 = (-t354 * t361 - t365 * t384) * t358;
t383 = t361 ^ 2 * t390;
t290 = -pkin(4) * t383 + t330 * pkin(8) - t344 * t329 + t292;
t360 = sin(qJ(5));
t364 = cos(qJ(5));
t287 = t364 * t289 - t360 * t290;
t323 = (-t365 * t360 - t361 * t364) * t389;
t301 = t323 * qJD(5) + t360 * t330 + t364 * t331;
t324 = (-t361 * t360 + t365 * t364) * t389;
t310 = -t323 * mrSges(6,1) + t324 * mrSges(6,2);
t342 = qJD(5) + t344;
t315 = -t342 * mrSges(6,2) + t323 * mrSges(6,3);
t341 = qJDD(5) + t343;
t284 = m(6) * t287 + t341 * mrSges(6,1) - t301 * mrSges(6,3) - t324 * t310 + t342 * t315;
t288 = t360 * t289 + t364 * t290;
t300 = -t324 * qJD(5) + t364 * t330 - t360 * t331;
t316 = t342 * mrSges(6,1) - t324 * mrSges(6,3);
t285 = m(6) * t288 - t341 * mrSges(6,2) + t300 * mrSges(6,3) + t323 * t310 - t342 * t316;
t277 = t364 * t284 + t360 * t285;
t382 = t361 * t389;
t291 = -t361 * t297 + t306;
t326 = -t344 * mrSges(5,2) - mrSges(5,3) * t382;
t328 = (t361 * mrSges(5,1) + t365 * mrSges(5,2)) * t389;
t275 = m(5) * t291 + t343 * mrSges(5,1) - t331 * mrSges(5,3) + t344 * t326 - t328 * t381 + t277;
t327 = t344 * mrSges(5,1) - mrSges(5,3) * t381;
t377 = -t360 * t284 + t364 * t285;
t276 = m(5) * t292 - t343 * mrSges(5,2) + t330 * mrSges(5,3) - t344 * t327 - t328 * t382 + t377;
t296 = t334 * t389 - t308;
t332 = (-mrSges(4,1) * t359 + t392) * t356;
t293 = -t330 * pkin(4) - pkin(8) * t383 + t329 * t381 + t296;
t369 = m(6) * t293 - t300 * mrSges(6,1) + t301 * mrSges(6,2) - t323 * t315 + t324 * t316;
t378 = -t358 * (m(4) * t308 - m(5) * t296 + t330 * mrSges(5,1) - t331 * mrSges(5,2) + (-t391 + (-t326 * t361 - t327 * t365 - t332) * t356) * t358 - t369) + t359 * (m(4) * t309 - t361 * t275 + t365 * t276 + (t332 * t356 + t391) * t359);
t274 = t365 * t275 + t361 * t276;
t318 = Ifges(5,6) * t344 + (t365 * Ifges(5,4) - t361 * Ifges(5,2)) * t389;
t319 = Ifges(5,5) * t344 + (t365 * Ifges(5,1) - t361 * Ifges(5,4)) * t389;
t374 = t365 * t318 + t361 * t319;
t313 = -t354 * pkin(2) + t373;
t370 = -m(4) * t313 + mrSges(4,1) * t388 - t274 + (t353 * t359 ^ 2 + t390) * mrSges(4,3);
t273 = t354 * t392 - t370;
t302 = Ifges(6,5) * t324 + Ifges(6,6) * t323 + Ifges(6,3) * t342;
t304 = Ifges(6,1) * t324 + Ifges(6,4) * t323 + Ifges(6,5) * t342;
t278 = -mrSges(6,1) * t293 + mrSges(6,3) * t288 + Ifges(6,4) * t301 + Ifges(6,2) * t300 + Ifges(6,6) * t341 - t324 * t302 + t342 * t304;
t303 = Ifges(6,4) * t324 + Ifges(6,2) * t323 + Ifges(6,6) * t342;
t279 = mrSges(6,2) * t293 - mrSges(6,3) * t287 + Ifges(6,1) * t301 + Ifges(6,4) * t300 + Ifges(6,5) * t341 + t323 * t302 - t342 * t303;
t333 = (Ifges(4,5) * t358 + Ifges(4,6) * t359) * t356;
t371 = -mrSges(6,1) * t287 + mrSges(6,2) * t288 - Ifges(6,5) * t301 - Ifges(6,6) * t300 - Ifges(6,3) * t341 - t324 * t303 + t323 * t304;
t368 = mrSges(5,1) * t291 - mrSges(5,2) * t292 + Ifges(5,5) * t331 + Ifges(5,6) * t330 + Ifges(5,3) * t343 + pkin(4) * t277 - t371;
t372 = -mrSges(3,2) * t322 + t358 * (t333 * t387 + mrSges(4,2) * t313 - mrSges(4,3) * t308 + t365 * (mrSges(5,2) * t296 - mrSges(5,3) * t291 + Ifges(5,1) * t331 + Ifges(5,4) * t330 + Ifges(5,5) * t343 - pkin(8) * t277 - t360 * t278 + t364 * t279 - t344 * t318) - t361 * (-mrSges(5,1) * t296 + mrSges(5,3) * t292 + Ifges(5,4) * t331 + Ifges(5,2) * t330 + Ifges(5,6) * t343 - pkin(4) * t369 + pkin(8) * t377 + t364 * t278 + t360 * t279 + t344 * t319) - pkin(7) * t274 + (Ifges(4,1) * t358 + Ifges(4,4) * t359) * t354) + t359 * (-t368 + Ifges(4,2) * t388 + (Ifges(4,4) * t354 + (-t333 - t374) * t356) * t358 + mrSges(4,3) * t309 - mrSges(4,1) * t313 - pkin(3) * t274) + qJ(3) * t378 - pkin(2) * t273 + mrSges(3,1) * t321 + Ifges(3,3) * t354;
t1 = [Ifges(2,3) * qJDD(1) + mrSges(2,1) * t379 - mrSges(2,2) * t375 + pkin(1) * (t362 * (m(3) * t322 - t353 * mrSges(3,1) - t354 * mrSges(3,2) + t378) + t366 * (m(3) * t321 - t353 * mrSges(3,2) + (mrSges(3,1) - t392) * t354 + t370)) + t372; t372; t273; t374 * t389 + t368; -t371;];
tauJ = t1;
