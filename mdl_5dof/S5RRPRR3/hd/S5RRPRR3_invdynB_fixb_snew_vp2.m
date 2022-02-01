% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRPRR3
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
% tauB [6x1]
%   base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:34
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRPRR3_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR3_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:34:10
% EndTime: 2022-01-20 10:34:12
% DurationCPUTime: 2.45s
% Computational Cost: add. (39679->186), mult. (47778->233), div. (0->0), fcn. (25560->10), ass. (0->81)
t366 = qJD(1) + qJD(2);
t360 = qJD(4) + t366;
t370 = sin(qJ(5));
t390 = t360 * t370;
t374 = cos(qJ(5));
t389 = t360 * t374;
t373 = sin(qJ(1));
t377 = cos(qJ(1));
t356 = t373 * g(1) - t377 * g(2);
t354 = qJDD(1) * pkin(1) + t356;
t357 = -t377 * g(1) - t373 * g(2);
t378 = qJD(1) ^ 2;
t355 = -t378 * pkin(1) + t357;
t372 = sin(qJ(2));
t376 = cos(qJ(2));
t339 = t376 * t354 - t372 * t355;
t365 = qJDD(1) + qJDD(2);
t337 = t365 * pkin(2) + t339;
t340 = t372 * t354 + t376 * t355;
t364 = t366 ^ 2;
t338 = -t364 * pkin(2) + t340;
t368 = sin(pkin(9));
t369 = cos(pkin(9));
t332 = t369 * t337 - t368 * t338;
t330 = t365 * pkin(3) + t332;
t333 = t368 * t337 + t369 * t338;
t331 = -t364 * pkin(3) + t333;
t371 = sin(qJ(4));
t375 = cos(qJ(4));
t327 = t371 * t330 + t375 * t331;
t358 = t360 ^ 2;
t359 = qJDD(4) + t365;
t325 = -t358 * pkin(4) + t359 * pkin(8) + t327;
t367 = -g(3) + qJDD(3);
t322 = -t370 * t325 + t374 * t367;
t346 = (-mrSges(6,1) * t374 + mrSges(6,2) * t370) * t360;
t387 = qJD(5) * t360;
t347 = t370 * t359 + t374 * t387;
t353 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t389;
t320 = m(6) * t322 + qJDD(5) * mrSges(6,1) - t347 * mrSges(6,3) + qJD(5) * t353 - t346 * t390;
t323 = t374 * t325 + t370 * t367;
t348 = t374 * t359 - t370 * t387;
t352 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t390;
t321 = m(6) * t323 - qJDD(5) * mrSges(6,2) + t348 * mrSges(6,3) - qJD(5) * t352 + t346 * t389;
t381 = -t370 * t320 + t374 * t321;
t311 = m(5) * t327 - t358 * mrSges(5,1) - t359 * mrSges(5,2) + t381;
t326 = t375 * t330 - t371 * t331;
t324 = -t359 * pkin(4) - t358 * pkin(8) - t326;
t379 = -m(6) * t324 + t348 * mrSges(6,1) - t347 * mrSges(6,2) - t352 * t390 + t353 * t389;
t316 = m(5) * t326 + t359 * mrSges(5,1) - t358 * mrSges(5,2) + t379;
t308 = t371 * t311 + t375 * t316;
t305 = m(4) * t332 + t365 * mrSges(4,1) - t364 * mrSges(4,2) + t308;
t382 = t375 * t311 - t371 * t316;
t306 = m(4) * t333 - t364 * mrSges(4,1) - t365 * mrSges(4,2) + t382;
t300 = t369 * t305 + t368 * t306;
t298 = m(3) * t339 + t365 * mrSges(3,1) - t364 * mrSges(3,2) + t300;
t383 = -t368 * t305 + t369 * t306;
t299 = m(3) * t340 - t364 * mrSges(3,1) - t365 * mrSges(3,2) + t383;
t292 = t376 * t298 + t372 * t299;
t290 = m(2) * t356 + qJDD(1) * mrSges(2,1) - t378 * mrSges(2,2) + t292;
t384 = -t372 * t298 + t376 * t299;
t291 = m(2) * t357 - t378 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t384;
t388 = t377 * t290 + t373 * t291;
t312 = t374 * t320 + t370 * t321;
t386 = m(5) * t367 + t312;
t385 = -t373 * t290 + t377 * t291;
t380 = m(4) * t367 + t386;
t343 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t370 + Ifges(6,4) * t374) * t360;
t342 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t370 + Ifges(6,2) * t374) * t360;
t341 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t370 + Ifges(6,6) * t374) * t360;
t314 = mrSges(6,2) * t324 - mrSges(6,3) * t322 + Ifges(6,1) * t347 + Ifges(6,4) * t348 + Ifges(6,5) * qJDD(5) - qJD(5) * t342 + t341 * t389;
t313 = -mrSges(6,1) * t324 + mrSges(6,3) * t323 + Ifges(6,4) * t347 + Ifges(6,2) * t348 + Ifges(6,6) * qJDD(5) + qJD(5) * t343 - t341 * t390;
t307 = -mrSges(5,1) * t367 - mrSges(6,1) * t322 + mrSges(6,2) * t323 + mrSges(5,3) * t327 + t358 * Ifges(5,5) - Ifges(6,5) * t347 + Ifges(5,6) * t359 - Ifges(6,6) * t348 - Ifges(6,3) * qJDD(5) - pkin(4) * t312 + (-t342 * t370 + t343 * t374) * t360;
t301 = mrSges(5,2) * t367 - mrSges(5,3) * t326 + Ifges(5,5) * t359 - t358 * Ifges(5,6) - pkin(8) * t312 - t370 * t313 + t374 * t314;
t294 = mrSges(4,2) * t367 - mrSges(4,3) * t332 + Ifges(4,5) * t365 - t364 * Ifges(4,6) - pkin(7) * t308 + t375 * t301 - t371 * t307;
t293 = -mrSges(4,1) * t367 + mrSges(4,3) * t333 + t364 * Ifges(4,5) + Ifges(4,6) * t365 - pkin(3) * t386 + pkin(7) * t382 + t371 * t301 + t375 * t307;
t286 = -mrSges(3,2) * g(3) - mrSges(3,3) * t339 + Ifges(3,5) * t365 - t364 * Ifges(3,6) - qJ(3) * t300 - t368 * t293 + t369 * t294;
t285 = mrSges(3,1) * g(3) + mrSges(3,3) * t340 + t364 * Ifges(3,5) + Ifges(3,6) * t365 - pkin(2) * t380 + qJ(3) * t383 + t369 * t293 + t368 * t294;
t284 = -mrSges(2,2) * g(3) - mrSges(2,3) * t356 + Ifges(2,5) * qJDD(1) - t378 * Ifges(2,6) - pkin(6) * t292 - t372 * t285 + t376 * t286;
t283 = Ifges(2,6) * qJDD(1) + t378 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t357 + t372 * t286 + t376 * t285 - pkin(1) * (-m(3) * g(3) + t380) + pkin(6) * t384;
t1 = [-m(1) * g(1) + t385; -m(1) * g(2) + t388; (-m(1) - m(2) - m(3)) * g(3) + t380; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(5) * t388 - t373 * t283 + t377 * t284; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(5) * t385 + t377 * t283 + t373 * t284; pkin(1) * t292 + mrSges(2,1) * t356 - mrSges(2,2) * t357 + pkin(2) * t300 + mrSges(3,1) * t339 - mrSges(3,2) * t340 + pkin(3) * t308 + mrSges(4,1) * t332 - mrSges(4,2) * t333 + pkin(8) * t381 + mrSges(5,1) * t326 - mrSges(5,2) * t327 + t370 * t314 + t374 * t313 + pkin(4) * t379 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(5,3) * t359 + Ifges(2,3) * qJDD(1) + (Ifges(3,3) + Ifges(4,3)) * t365;];
tauB = t1;
