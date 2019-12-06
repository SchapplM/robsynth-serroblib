% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S5RRRPR2
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
% Datum: 2019-12-05 18:41
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S5RRRPR2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR2_invdynB_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:40:53
% EndTime: 2019-12-05 18:40:55
% DurationCPUTime: 2.18s
% Computational Cost: add. (41161->187), mult. (47778->233), div. (0->0), fcn. (25560->10), ass. (0->81)
t392 = -m(3) - m(4);
t368 = qJD(1) + qJD(2);
t361 = qJD(3) + t368;
t372 = sin(qJ(5));
t391 = t361 * t372;
t376 = cos(qJ(5));
t390 = t361 * t376;
t375 = sin(qJ(1));
t379 = cos(qJ(1));
t358 = t379 * g(2) + t375 * g(3);
t355 = qJDD(1) * pkin(1) + t358;
t357 = t375 * g(2) - t379 * g(3);
t380 = qJD(1) ^ 2;
t356 = -t380 * pkin(1) + t357;
t374 = sin(qJ(2));
t378 = cos(qJ(2));
t340 = t378 * t355 - t374 * t356;
t367 = qJDD(1) + qJDD(2);
t338 = t367 * pkin(2) + t340;
t341 = t374 * t355 + t378 * t356;
t366 = t368 ^ 2;
t339 = -t366 * pkin(2) + t341;
t373 = sin(qJ(3));
t377 = cos(qJ(3));
t333 = t377 * t338 - t373 * t339;
t360 = qJDD(3) + t367;
t331 = t360 * pkin(3) + t333;
t334 = t373 * t338 + t377 * t339;
t359 = t361 ^ 2;
t332 = -t359 * pkin(3) + t334;
t370 = sin(pkin(9));
t371 = cos(pkin(9));
t328 = t370 * t331 + t371 * t332;
t326 = -t359 * pkin(4) + t360 * pkin(8) + t328;
t369 = -g(1) + qJDD(4);
t323 = -t372 * t326 + t376 * t369;
t347 = (-mrSges(6,1) * t376 + mrSges(6,2) * t372) * t361;
t389 = qJD(5) * t361;
t348 = t372 * t360 + t376 * t389;
t354 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t390;
t321 = m(6) * t323 + qJDD(5) * mrSges(6,1) - t348 * mrSges(6,3) + qJD(5) * t354 - t347 * t391;
t324 = t376 * t326 + t372 * t369;
t349 = t376 * t360 - t372 * t389;
t353 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t391;
t322 = m(6) * t324 - qJDD(5) * mrSges(6,2) + t349 * mrSges(6,3) - qJD(5) * t353 + t347 * t390;
t383 = -t372 * t321 + t376 * t322;
t312 = m(5) * t328 - t359 * mrSges(5,1) - t360 * mrSges(5,2) + t383;
t327 = t371 * t331 - t370 * t332;
t325 = -t360 * pkin(4) - t359 * pkin(8) - t327;
t381 = -m(6) * t325 + t349 * mrSges(6,1) - t348 * mrSges(6,2) - t353 * t391 + t354 * t390;
t317 = m(5) * t327 + t360 * mrSges(5,1) - t359 * mrSges(5,2) + t381;
t309 = t370 * t312 + t371 * t317;
t306 = m(4) * t333 + t360 * mrSges(4,1) - t359 * mrSges(4,2) + t309;
t384 = t371 * t312 - t370 * t317;
t307 = m(4) * t334 - t359 * mrSges(4,1) - t360 * mrSges(4,2) + t384;
t301 = t377 * t306 + t373 * t307;
t299 = m(3) * t340 + t367 * mrSges(3,1) - t366 * mrSges(3,2) + t301;
t385 = -t373 * t306 + t377 * t307;
t300 = m(3) * t341 - t366 * mrSges(3,1) - t367 * mrSges(3,2) + t385;
t293 = t378 * t299 + t374 * t300;
t313 = t376 * t321 + t372 * t322;
t388 = m(5) * t369 + t313;
t386 = -t374 * t299 + t378 * t300;
t291 = m(2) * t357 - t380 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t386;
t292 = m(2) * t358 + qJDD(1) * mrSges(2,1) - t380 * mrSges(2,2) + t293;
t387 = t379 * t291 - t375 * t292;
t382 = -t375 * t291 - t379 * t292;
t344 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t372 + Ifges(6,4) * t376) * t361;
t343 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t372 + Ifges(6,2) * t376) * t361;
t342 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t372 + Ifges(6,6) * t376) * t361;
t315 = mrSges(6,2) * t325 - mrSges(6,3) * t323 + Ifges(6,1) * t348 + Ifges(6,4) * t349 + Ifges(6,5) * qJDD(5) - qJD(5) * t343 + t342 * t390;
t314 = -mrSges(6,1) * t325 + mrSges(6,3) * t324 + Ifges(6,4) * t348 + Ifges(6,2) * t349 + Ifges(6,6) * qJDD(5) + qJD(5) * t344 - t342 * t391;
t308 = -mrSges(5,1) * t369 - mrSges(6,1) * t323 + mrSges(6,2) * t324 + mrSges(5,3) * t328 + t359 * Ifges(5,5) - Ifges(6,5) * t348 + Ifges(5,6) * t360 - Ifges(6,6) * t349 - Ifges(6,3) * qJDD(5) - pkin(4) * t313 + (-t343 * t372 + t344 * t376) * t361;
t302 = mrSges(5,2) * t369 - mrSges(5,3) * t327 + Ifges(5,5) * t360 - t359 * Ifges(5,6) - pkin(8) * t313 - t372 * t314 + t376 * t315;
t295 = -mrSges(4,2) * g(1) - mrSges(4,3) * t333 + Ifges(4,5) * t360 - t359 * Ifges(4,6) - qJ(4) * t309 + t371 * t302 - t370 * t308;
t294 = mrSges(4,1) * g(1) + mrSges(4,3) * t334 + t359 * Ifges(4,5) + Ifges(4,6) * t360 - pkin(3) * t388 + qJ(4) * t384 + t370 * t302 + t371 * t308;
t289 = -mrSges(3,2) * g(1) - mrSges(3,3) * t340 + Ifges(3,5) * t367 - t366 * Ifges(3,6) - pkin(7) * t301 - t373 * t294 + t377 * t295;
t288 = Ifges(3,6) * t367 + t366 * Ifges(3,5) + mrSges(3,1) * g(1) + mrSges(3,3) * t341 + t373 * t295 + t377 * t294 - pkin(2) * (-m(4) * g(1) + t388) + pkin(7) * t385;
t287 = -mrSges(2,2) * g(1) - mrSges(2,3) * t358 + Ifges(2,5) * qJDD(1) - t380 * Ifges(2,6) - pkin(6) * t293 - t374 * t288 + t378 * t289;
t286 = Ifges(2,6) * qJDD(1) + t380 * Ifges(2,5) + mrSges(2,3) * t357 + t374 * t289 + t378 * t288 - pkin(1) * t388 + pkin(6) * t386 + (-pkin(1) * t392 + mrSges(2,1)) * g(1);
t1 = [(-m(1) - m(2) + t392) * g(1) + t388; -m(1) * g(2) + t382; -m(1) * g(3) + t387; pkin(1) * t293 + pkin(2) * t301 - mrSges(3,2) * t341 + mrSges(3,1) * t340 + pkin(3) * t309 + mrSges(4,1) * t333 - mrSges(4,2) * t334 + t372 * t315 + t376 * t314 + pkin(4) * t381 + pkin(8) * t383 + mrSges(5,1) * t327 - mrSges(5,2) * t328 - mrSges(2,2) * t357 + mrSges(2,1) * t358 + Ifges(2,3) * qJDD(1) + Ifges(3,3) * t367 - mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + (Ifges(5,3) + Ifges(4,3)) * t360; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - pkin(5) * t387 - t379 * t286 - t375 * t287; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + pkin(5) * t382 - t375 * t286 + t379 * t287;];
tauB = t1;
