% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauJ [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 13:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S6RPPPRR1_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_invdynJ_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:32:02
% EndTime: 2019-05-05 13:32:03
% DurationCPUTime: 0.63s
% Computational Cost: add. (2671->177), mult. (4702->216), div. (0->0), fcn. (2272->8), ass. (0->74)
t370 = sin(qJ(1));
t373 = cos(qJ(1));
t385 = t370 * g(1) - g(2) * t373;
t342 = qJDD(1) * pkin(1) + t385;
t375 = qJD(1) ^ 2;
t383 = -g(1) * t373 - g(2) * t370;
t344 = -pkin(1) * t375 + t383;
t366 = sin(pkin(9));
t367 = cos(pkin(9));
t324 = t367 * t342 - t366 * t344;
t316 = -qJDD(1) * pkin(2) - t375 * qJ(3) + qJDD(3) - t324;
t313 = -qJDD(1) * qJ(4) - (2 * qJD(4) * qJD(1)) + t316;
t325 = t366 * t342 + t367 * t344;
t396 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t325;
t361 = -g(3) + qJDD(2);
t369 = sin(qJ(5));
t395 = t361 * t369;
t314 = qJDD(4) + (-pkin(2) - qJ(4)) * t375 + t396;
t311 = -qJDD(1) * pkin(7) + t314;
t372 = cos(qJ(5));
t307 = t369 * t311 + t372 * t361;
t343 = (t369 * mrSges(6,1) + t372 * mrSges(6,2)) * qJD(1);
t391 = qJD(1) * qJD(5);
t386 = t372 * t391;
t346 = -qJDD(1) * t369 - t386;
t392 = qJD(1) * t372;
t349 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t392;
t310 = -pkin(7) * t375 - t313;
t387 = t369 * t391;
t347 = qJDD(1) * t372 - t387;
t303 = (-t347 + t387) * pkin(8) + (-t346 + t386) * pkin(5) + t310;
t345 = (t369 * pkin(5) - t372 * pkin(8)) * qJD(1);
t374 = qJD(5) ^ 2;
t393 = qJD(1) * t369;
t305 = -pkin(5) * t374 + qJDD(5) * pkin(8) - t345 * t393 + t307;
t368 = sin(qJ(6));
t371 = cos(qJ(6));
t301 = t303 * t371 - t305 * t368;
t340 = qJD(5) * t371 - t368 * t392;
t323 = qJD(6) * t340 + qJDD(5) * t368 + t347 * t371;
t341 = qJD(5) * t368 + t371 * t392;
t326 = -mrSges(7,1) * t340 + mrSges(7,2) * t341;
t350 = qJD(6) + t393;
t327 = -mrSges(7,2) * t350 + mrSges(7,3) * t340;
t339 = qJDD(6) - t346;
t299 = m(7) * t301 + mrSges(7,1) * t339 - mrSges(7,3) * t323 - t326 * t341 + t327 * t350;
t302 = t303 * t368 + t305 * t371;
t322 = -qJD(6) * t341 + qJDD(5) * t371 - t347 * t368;
t328 = mrSges(7,1) * t350 - mrSges(7,3) * t341;
t300 = m(7) * t302 - mrSges(7,2) * t339 + mrSges(7,3) * t322 + t326 * t340 - t328 * t350;
t384 = -t299 * t368 + t371 * t300;
t290 = m(6) * t307 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t346 - qJD(5) * t349 - t343 * t393 + t384;
t306 = t311 * t372 - t395;
t348 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t393;
t304 = -qJDD(5) * pkin(5) - pkin(8) * t374 + t395 + (qJD(1) * t345 - t311) * t372;
t380 = -m(7) * t304 + t322 * mrSges(7,1) - mrSges(7,2) * t323 + t340 * t327 - t328 * t341;
t295 = m(6) * t306 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t347 + qJD(5) * t348 - t343 * t392 + t380;
t394 = t369 * t290 + t372 * t295;
t291 = t371 * t299 + t368 * t300;
t381 = m(5) * t314 + qJDD(1) * mrSges(5,2) - mrSges(5,3) * t375 + t394;
t315 = pkin(2) * t375 - t396;
t379 = -m(4) * t315 + t375 * mrSges(4,2) + qJDD(1) * mrSges(4,3) + t381;
t378 = m(5) * t313 - m(6) * t310 + mrSges(6,1) * t346 - t375 * mrSges(5,2) - t347 * mrSges(6,2) - qJDD(1) * mrSges(5,3) - t348 * t393 - t349 * t392 - t291;
t377 = m(4) * t316 - t375 * mrSges(4,3) + t378;
t318 = Ifges(7,4) * t341 + Ifges(7,2) * t340 + Ifges(7,6) * t350;
t319 = Ifges(7,1) * t341 + Ifges(7,4) * t340 + Ifges(7,5) * t350;
t376 = mrSges(7,1) * t301 - mrSges(7,2) * t302 + Ifges(7,5) * t323 + Ifges(7,6) * t322 + Ifges(7,3) * t339 + t318 * t341 - t340 * t319;
t335 = (Ifges(6,5) * qJD(5)) + (t372 * Ifges(6,1) - t369 * Ifges(6,4)) * qJD(1);
t334 = (Ifges(6,6) * qJD(5)) + (t372 * Ifges(6,4) - t369 * Ifges(6,2)) * qJD(1);
t317 = Ifges(7,5) * t341 + Ifges(7,6) * t340 + Ifges(7,3) * t350;
t293 = mrSges(7,2) * t304 - mrSges(7,3) * t301 + Ifges(7,1) * t323 + Ifges(7,4) * t322 + Ifges(7,5) * t339 + t317 * t340 - t318 * t350;
t292 = -mrSges(7,1) * t304 + mrSges(7,3) * t302 + Ifges(7,4) * t323 + Ifges(7,2) * t322 + Ifges(7,6) * t339 - t317 * t341 + t319 * t350;
t288 = qJDD(1) * mrSges(4,2) + t377;
t1 = [pkin(1) * (t366 * (m(3) * t325 - mrSges(3,1) * t375 + t379) + t367 * (m(3) * t324 - mrSges(3,2) * t375 - t377)) - mrSges(2,2) * t383 + mrSges(2,1) * t385 - pkin(2) * t288 + mrSges(3,1) * t324 - mrSges(3,2) * t325 + qJ(3) * t379 - qJ(4) * t378 + mrSges(4,2) * t316 - mrSges(4,3) * t315 + t372 * (mrSges(6,2) * t310 - mrSges(6,3) * t306 + Ifges(6,1) * t347 + Ifges(6,4) * t346 + Ifges(6,5) * qJDD(5) - pkin(8) * t291 - qJD(5) * t334 - t292 * t368 + t293 * t371) - t369 * (-mrSges(6,1) * t310 + mrSges(6,3) * t307 + Ifges(6,4) * t347 + Ifges(6,2) * t346 + Ifges(6,6) * qJDD(5) - pkin(5) * t291 + qJD(5) * t335 - t376) - pkin(7) * t394 + mrSges(5,2) * t314 - mrSges(5,3) * t313 + (pkin(1) * (-t366 * mrSges(3,2) + t367 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(3,3) + Ifges(4,1) + Ifges(5,1)) * qJDD(1); t290 * t372 - t295 * t369 + (m(3) + m(4) + m(5)) * t361; t288; t381; Ifges(6,5) * t347 + Ifges(6,6) * t346 + Ifges(6,3) * qJDD(5) + mrSges(6,1) * t306 - mrSges(6,2) * t307 + t368 * t293 + t371 * t292 + pkin(5) * t380 + pkin(8) * t384 + (t372 * t334 + t369 * t335) * qJD(1); t376;];
tauJ  = t1;
