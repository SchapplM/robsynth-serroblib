% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRRPP4_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:51
% EndTime: 2019-12-31 17:40:52
% DurationCPUTime: 0.74s
% Computational Cost: add. (1015->168), mult. (2009->198), div. (0->0), fcn. (926->6), ass. (0->74)
t355 = sin(qJ(3));
t357 = cos(qJ(3));
t374 = Ifges(4,4) - Ifges(5,5) - Ifges(6,4);
t396 = t355 * (Ifges(4,1) + Ifges(5,1) + Ifges(6,1)) + t357 * t374;
t395 = -t357 * (Ifges(4,2) + Ifges(5,3) + Ifges(6,2)) - t355 * t374;
t373 = Ifges(4,5) + Ifges(5,4) - Ifges(6,5);
t372 = Ifges(4,6) - Ifges(5,6) + Ifges(6,6);
t390 = t355 * (t395 * qJD(2) - t372 * qJD(3));
t389 = 2 * qJD(4);
t388 = pkin(3) + pkin(4);
t360 = qJD(2) ^ 2;
t387 = pkin(4) * t360;
t386 = pkin(6) * t360;
t385 = -mrSges(6,2) - mrSges(5,3);
t384 = mrSges(4,3) + mrSges(5,2);
t383 = qJ(4) * t357;
t352 = -g(3) + qJDD(1);
t382 = t352 * t357;
t353 = sin(pkin(7));
t354 = cos(pkin(7));
t334 = g(1) * t353 - g(2) * t354;
t335 = -g(1) * t354 - g(2) * t353;
t356 = sin(qJ(2));
t358 = cos(qJ(2));
t380 = t356 * t334 + t358 * t335;
t303 = -pkin(2) * t360 + qJDD(2) * pkin(6) + t380;
t299 = t357 * t303 + t355 * t352;
t381 = t358 * t334 - t356 * t335;
t377 = qJD(2) * t355;
t337 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t377;
t339 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t377;
t379 = -t337 - t339;
t376 = qJD(2) * t357;
t340 = qJD(3) * mrSges(6,2) - mrSges(6,3) * t376;
t342 = mrSges(5,2) * t376 + qJD(3) * mrSges(5,3);
t378 = -t340 - t342;
t375 = qJD(2) * qJD(3);
t371 = -0.2e1 * qJD(2) * qJD(5);
t369 = t396 * qJD(2) + t373 * qJD(3);
t329 = qJDD(2) * t357 - t355 * t375;
t336 = -qJD(3) * pkin(4) - qJ(5) * t377;
t351 = t357 ^ 2;
t328 = qJDD(2) * t355 + t357 * t375;
t366 = qJDD(2) * pkin(2) + t381;
t363 = -qJ(4) * t328 - t366;
t290 = qJDD(5) + (-qJ(5) * t351 + pkin(6)) * t360 + t388 * t329 + (qJD(3) * t383 + (-pkin(3) * qJD(3) + t336 + t389) * t355) * qJD(2) - t363;
t368 = m(6) * t290 + t329 * mrSges(6,1);
t324 = (-t357 * pkin(3) - t355 * qJ(4)) * qJD(2);
t359 = qJD(3) ^ 2;
t296 = -pkin(3) * t359 + qJDD(3) * qJ(4) + qJD(3) * t389 + t324 * t376 + t299;
t292 = -qJ(5) * t329 + qJD(3) * t336 - t351 * t387 + t357 * t371 + t296;
t367 = m(6) * t292 + qJDD(3) * mrSges(6,2) - t329 * mrSges(6,3) + qJD(3) * t337;
t294 = -pkin(3) * t329 - t386 + (-0.2e1 * qJD(4) * t355 + (pkin(3) * t355 - t383) * qJD(3)) * qJD(2) + t363;
t365 = m(5) * t294 - t368;
t300 = t355 * t303;
t364 = -qJ(4) * t359 + t324 * t377 + qJDD(4) + t300;
t293 = t355 * t371 - qJ(5) * t328 - t388 * qJDD(3) + (qJ(5) * t375 - t355 * t387 - t352) * t357 + t364;
t326 = (t357 * mrSges(6,1) + t355 * mrSges(6,2)) * qJD(2);
t289 = m(6) * t293 - qJDD(3) * mrSges(6,1) - t328 * mrSges(6,3) - qJD(3) * t340 - t326 * t377;
t325 = (-t357 * mrSges(5,1) - t355 * mrSges(5,3)) * qJD(2);
t362 = m(5) * t296 + qJDD(3) * mrSges(5,3) + qJD(3) * t339 + t325 * t376 + t367;
t297 = -qJDD(3) * pkin(3) + t364 - t382;
t361 = -m(5) * t297 + qJDD(3) * mrSges(5,1) + qJD(3) * t342 - t289;
t341 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t376;
t338 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t377;
t327 = (-t357 * mrSges(4,1) + t355 * mrSges(4,2)) * qJD(2);
t302 = -t366 - t386;
t298 = -t300 + t382;
t288 = mrSges(6,2) * t328 + (t337 * t355 + t340 * t357) * qJD(2) + t368;
t287 = mrSges(5,2) * t328 + t325 * t377 - t361;
t286 = -mrSges(5,1) * t329 + t385 * t328 + (t379 * t355 + t378 * t357) * qJD(2) + t365;
t285 = m(4) * t298 + qJDD(3) * mrSges(4,1) + qJD(3) * t341 - t384 * t328 + (-t325 - t327) * t377 + t361;
t284 = m(4) * t299 - qJDD(3) * mrSges(4,2) - qJD(3) * t338 + t384 * t329 + (-t326 + t327) * t376 + t362;
t1 = [t284 * t355 + t285 * t357 + (m(2) + m(3)) * t352; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t381 - mrSges(3,2) * t380 + t355 * (mrSges(4,2) * t302 + mrSges(5,2) * t297 + mrSges(6,2) * t290 - mrSges(4,3) * t298 - mrSges(5,3) * t294 - mrSges(6,3) * t293 - qJ(4) * t286 - qJ(5) * t289) + t357 * (-mrSges(4,1) * t302 - mrSges(5,1) * t294 + mrSges(6,1) * t290 + mrSges(5,2) * t296 + mrSges(4,3) * t299 - mrSges(6,3) * t292 - pkin(3) * t286 + pkin(4) * t288 - qJ(5) * t367) + pkin(2) * (-m(4) * t302 - t365) + pkin(6) * (t284 * t357 - t285 * t355) + (t355 * t373 + t357 * t372) * qJDD(3) + (t357 * t369 + t390) * qJD(3) + (pkin(2) * (mrSges(4,1) + mrSges(5,1)) - t395) * t329 + (pkin(2) * (-mrSges(4,2) - t385) + t396) * t328 + (pkin(2) * (-t338 - t379) * t355 + (pkin(2) * (t341 - t378) + qJ(5) * t326 * t357) * t357) * qJD(2); mrSges(4,1) * t298 - mrSges(4,2) * t299 - mrSges(5,1) * t297 + mrSges(5,3) * t296 - mrSges(6,1) * t293 + mrSges(6,2) * t292 - pkin(4) * t289 - pkin(3) * t287 + qJ(4) * t362 + (qJ(4) * mrSges(5,2) + t372) * t329 + t373 * t328 + (Ifges(4,3) + Ifges(5,2) + Ifges(6,3)) * qJDD(3) + (-t390 + (-qJ(4) * t326 - t369) * t357) * qJD(2); t287; t288;];
tauJ = t1;
