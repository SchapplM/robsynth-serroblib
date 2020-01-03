% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RRRP2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RRRP2_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:57
% EndTime: 2019-12-31 17:12:58
% DurationCPUTime: 0.89s
% Computational Cost: add. (6405->179), mult. (8185->215), div. (0->0), fcn. (3674->6), ass. (0->76)
t392 = Ifges(4,1) + Ifges(5,1);
t384 = Ifges(4,4) + Ifges(5,4);
t383 = Ifges(4,5) + Ifges(5,5);
t391 = Ifges(4,2) + Ifges(5,2);
t390 = Ifges(4,6) + Ifges(5,6);
t389 = Ifges(4,3) + Ifges(5,3);
t354 = qJD(1) + qJD(2);
t352 = t354 ^ 2;
t359 = sin(qJ(1));
t362 = cos(qJ(1));
t349 = t359 * g(1) - t362 * g(2);
t342 = qJDD(1) * pkin(1) + t349;
t350 = -t362 * g(1) - t359 * g(2);
t363 = qJD(1) ^ 2;
t343 = -t363 * pkin(1) + t350;
t358 = sin(qJ(2));
t361 = cos(qJ(2));
t320 = t361 * t342 - t358 * t343;
t353 = qJDD(1) + qJDD(2);
t365 = -t353 * pkin(2) - t320;
t318 = -t352 * pkin(6) + t365;
t357 = sin(qJ(3));
t360 = cos(qJ(3));
t374 = qJD(3) * t354;
t370 = t360 * t374;
t336 = t357 * t353 + t370;
t337 = t360 * t353 - t357 * t374;
t380 = t354 * t360;
t348 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t380;
t381 = t354 * t357;
t344 = qJD(3) * pkin(3) - qJ(4) * t381;
t356 = t360 ^ 2;
t314 = t344 * t381 - t337 * pkin(3) + qJDD(4) + (-qJ(4) * t356 - pkin(6)) * t352 + t365;
t347 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t380;
t366 = m(5) * t314 - t337 * mrSges(5,1) - t347 * t380;
t385 = -mrSges(4,2) - mrSges(5,2);
t388 = -m(4) * t318 + t337 * mrSges(4,1) + t385 * t336 + t348 * t380 - t366;
t387 = pkin(3) * t352;
t386 = t360 * g(3);
t321 = t358 * t342 + t361 * t343;
t319 = -t352 * pkin(2) + t353 * pkin(6) + t321;
t316 = -t357 * g(3) + t360 * t319;
t335 = (-mrSges(4,1) * t360 + mrSges(4,2) * t357) * t354;
t373 = qJD(4) * t354;
t313 = t337 * qJ(4) - qJD(3) * t344 - t356 * t387 + 0.2e1 * t360 * t373 + t316;
t334 = (-mrSges(5,1) * t360 + mrSges(5,2) * t357) * t354;
t371 = m(5) * t313 + t337 * mrSges(5,3) + t334 * t380;
t345 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t381;
t375 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t381 - t345;
t308 = m(4) * t316 + t337 * mrSges(4,3) + t375 * qJD(3) + t385 * qJDD(3) + t335 * t380 + t371;
t306 = t360 * t308;
t315 = -t357 * t319 - t386;
t312 = qJDD(3) * pkin(3) - t386 + (-t336 + t370) * qJ(4) + (t360 * t387 - t319 - 0.2e1 * t373) * t357;
t372 = m(5) * t312 + qJDD(3) * mrSges(5,1) + qJD(3) * t347;
t307 = m(4) * t315 + qJDD(3) * mrSges(4,1) + qJD(3) * t348 + (-t334 - t335) * t381 + (-mrSges(4,3) - mrSges(5,3)) * t336 + t372;
t300 = m(3) * t321 - t352 * mrSges(3,1) - t353 * mrSges(3,2) - t357 * t307 + t306;
t369 = t354 * t375;
t303 = m(3) * t320 + t353 * mrSges(3,1) - t352 * mrSges(3,2) + t357 * t369 + t388;
t295 = t358 * t300 + t361 * t303;
t293 = m(2) * t349 + qJDD(1) * mrSges(2,1) - t363 * mrSges(2,2) + t295;
t367 = t361 * t300 - t358 * t303;
t294 = m(2) * t350 - t363 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t367;
t379 = t362 * t293 + t359 * t294;
t301 = t360 * t307 + t357 * t308;
t378 = (t383 * t357 + t390 * t360) * t354 + t389 * qJD(3);
t377 = (-t384 * t357 - t391 * t360) * t354 - t390 * qJD(3);
t376 = (t392 * t357 + t384 * t360) * t354 + t383 * qJD(3);
t368 = -t359 * t293 + t362 * t294;
t309 = -t336 * mrSges(5,3) - t334 * t381 + t372;
t297 = mrSges(4,2) * t318 + mrSges(5,2) * t314 - mrSges(4,3) * t315 - mrSges(5,3) * t312 - qJ(4) * t309 + t377 * qJD(3) + t383 * qJDD(3) + t392 * t336 + t384 * t337 + t378 * t380;
t296 = -mrSges(4,1) * t318 + mrSges(4,3) * t316 - mrSges(5,1) * t314 + mrSges(5,3) * t313 - pkin(3) * t366 + qJ(4) * t371 + (-pkin(3) * t345 - t378) * t381 + t391 * t337 + (-pkin(3) * mrSges(5,2) + t384) * t336 + (-qJ(4) * mrSges(5,2) + t390) * qJDD(3) + (-qJ(4) * t345 + t376) * qJD(3);
t289 = mrSges(3,1) * g(3) - mrSges(4,1) * t315 - mrSges(5,1) * t312 + mrSges(4,2) * t316 + mrSges(5,2) * t313 + mrSges(3,3) * t321 + t352 * Ifges(3,5) + Ifges(3,6) * t353 - pkin(2) * t301 - pkin(3) * t309 - t390 * t337 - t383 * t336 - t389 * qJDD(3) + (t377 * t357 + t376 * t360) * t354;
t288 = -mrSges(3,2) * g(3) - mrSges(3,3) * t320 + Ifges(3,5) * t353 - t352 * Ifges(3,6) - pkin(6) * t301 - t357 * t296 + t360 * t297;
t287 = -mrSges(2,2) * g(3) - mrSges(2,3) * t349 + Ifges(2,5) * qJDD(1) - t363 * Ifges(2,6) - pkin(5) * t295 + t361 * t288 - t358 * t289;
t286 = Ifges(2,6) * qJDD(1) + t363 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t350 + t358 * t288 + t361 * t289 - pkin(1) * (-m(3) * g(3) + t301) + pkin(5) * t367;
t1 = [-m(1) * g(1) + t368; -m(1) * g(2) + t379; (-m(1) - m(2) - m(3)) * g(3) + t301; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t379 - t359 * t286 + t362 * t287; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t368 + t362 * t286 + t359 * t287; pkin(1) * t295 + mrSges(2,1) * t349 - mrSges(2,2) * t350 + t360 * t296 + pkin(2) * t388 + pkin(6) * t306 + mrSges(3,1) * t320 - mrSges(3,2) * t321 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + Ifges(3,3) * t353 + Ifges(2,3) * qJDD(1) + (pkin(2) * t369 - pkin(6) * t307 + t297) * t357;];
tauB = t1;
