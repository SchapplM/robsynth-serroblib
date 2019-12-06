% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRP5_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP5_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:37:45
% EndTime: 2019-12-05 15:37:47
% DurationCPUTime: 0.62s
% Computational Cost: add. (2278->161), mult. (5065->198), div. (0->0), fcn. (3255->8), ass. (0->79)
t381 = Ifges(5,1) + Ifges(6,1);
t374 = Ifges(5,4) - Ifges(6,5);
t373 = Ifges(5,5) + Ifges(6,4);
t380 = -Ifges(5,2) - Ifges(6,3);
t372 = Ifges(5,6) - Ifges(6,6);
t379 = Ifges(5,3) + Ifges(6,2);
t348 = qJD(2) ^ 2;
t342 = cos(pkin(8));
t336 = t342 ^ 2;
t378 = 0.2e1 * t342;
t377 = cos(qJ(4));
t376 = pkin(3) * t342;
t375 = -mrSges(5,3) - mrSges(6,2);
t340 = sin(pkin(8));
t371 = t340 * mrSges(4,2);
t370 = t336 * t348;
t341 = sin(pkin(7));
t343 = cos(pkin(7));
t326 = -g(1) * t343 - g(2) * t341;
t339 = -g(3) + qJDD(1);
t345 = sin(qJ(2));
t346 = cos(qJ(2));
t316 = t346 * t326 + t345 * t339;
t308 = -pkin(2) * t348 + qJDD(2) * qJ(3) + t316;
t325 = -g(1) * t341 + g(2) * t343;
t361 = qJD(2) * qJD(3);
t364 = t342 * t325 - 0.2e1 * t340 * t361;
t284 = (-pkin(6) * qJDD(2) + t348 * t376 - t308) * t340 + t364;
t287 = t342 * t308 + t340 * t325 + t361 * t378;
t359 = t342 * qJDD(2);
t285 = -pkin(3) * t370 + pkin(6) * t359 + t287;
t344 = sin(qJ(4));
t281 = t344 * t284 + t377 * t285;
t357 = t342 * t377;
t360 = t340 * qJDD(2);
t351 = t377 * t340 + t342 * t344;
t318 = t351 * qJD(2);
t363 = qJD(4) * t318;
t304 = -qJDD(2) * t357 + t344 * t360 + t363;
t312 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t318;
t317 = (t340 * t344 - t357) * qJD(2);
t298 = pkin(4) * t317 - qJ(5) * t318;
t347 = qJD(4) ^ 2;
t276 = -pkin(4) * t347 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t298 * t317 + t281;
t313 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t318;
t358 = m(6) * t276 + qJDD(4) * mrSges(6,3) + qJD(4) * t313;
t299 = mrSges(6,1) * t317 - mrSges(6,3) * t318;
t365 = -mrSges(5,1) * t317 - mrSges(5,2) * t318 - t299;
t271 = m(5) * t281 - qJDD(4) * mrSges(5,2) - qJD(4) * t312 + t375 * t304 + t365 * t317 + t358;
t280 = t377 * t284 - t344 * t285;
t362 = t317 * qJD(4);
t305 = t351 * qJDD(2) - t362;
t311 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t317;
t277 = -qJDD(4) * pkin(4) - t347 * qJ(5) + t318 * t298 + qJDD(5) - t280;
t314 = -mrSges(6,2) * t317 + qJD(4) * mrSges(6,3);
t354 = -m(6) * t277 + qJDD(4) * mrSges(6,1) + qJD(4) * t314;
t272 = m(5) * t280 + qJDD(4) * mrSges(5,1) + qJD(4) * t311 + t375 * t305 + t365 * t318 + t354;
t369 = t344 * t271 + t377 * t272;
t368 = -t379 * qJD(4) + t372 * t317 - t373 * t318;
t367 = t372 * qJD(4) + t380 * t317 + t374 * t318;
t366 = t373 * qJD(4) - t374 * t317 + t381 * t318;
t286 = -t308 * t340 + t364;
t352 = mrSges(4,3) * qJDD(2) + t348 * (-t342 * mrSges(4,1) + t371);
t355 = t377 * t271 - t272 * t344;
t356 = -(m(4) * t286 - t352 * t340 + t369) * t340 + t342 * (m(4) * t287 + t352 * t342 + t355);
t315 = -t345 * t326 + t339 * t346;
t353 = qJDD(3) - t315;
t335 = t340 ^ 2;
t288 = (-pkin(2) - t376) * qJDD(2) + (-qJ(3) + (-t335 - t336) * pkin(6)) * t348 + t353;
t279 = -0.2e1 * qJD(5) * t318 + (-t305 + t362) * qJ(5) + (t304 + t363) * pkin(4) + t288;
t273 = m(6) * t279 + t304 * mrSges(6,1) - t305 * mrSges(6,3) - t318 * t313 + t317 * t314;
t350 = m(5) * t288 + t304 * mrSges(5,1) + mrSges(5,2) * t305 + t317 * t311 + t312 * t318 + t273;
t307 = -qJDD(2) * pkin(2) - qJ(3) * t348 + t353;
t349 = -m(4) * t307 + mrSges(4,1) * t359 - t350 + (t335 * t348 + t370) * mrSges(4,3);
t274 = mrSges(6,2) * t305 + t299 * t318 - t354;
t267 = mrSges(4,2) * t360 - t349;
t264 = mrSges(5,2) * t288 + mrSges(6,2) * t277 - mrSges(5,3) * t280 - mrSges(6,3) * t279 - qJ(5) * t273 - t367 * qJD(4) + t373 * qJDD(4) - t374 * t304 + t381 * t305 + t368 * t317;
t263 = -mrSges(5,1) * t288 - mrSges(6,1) * t279 + mrSges(6,2) * t276 + mrSges(5,3) * t281 - pkin(4) * t273 + t366 * qJD(4) + t372 * qJDD(4) + t380 * t304 + t374 * t305 + t368 * t318;
t1 = [m(2) * t339 + t345 * (m(3) * t316 - mrSges(3,1) * t348 - qJDD(2) * mrSges(3,2) + t356) + t346 * (t349 + (mrSges(3,1) - t371) * qJDD(2) + m(3) * t315 - mrSges(3,2) * t348); mrSges(3,1) * t315 - mrSges(3,2) * t316 + t340 * (mrSges(4,2) * t307 - mrSges(4,3) * t286 - pkin(6) * t369 - t344 * t263 + t377 * t264) + t342 * (-mrSges(4,1) * t307 + mrSges(4,3) * t287 - pkin(3) * t350 + pkin(6) * t355 + t377 * t263 + t344 * t264) - pkin(2) * t267 + qJ(3) * t356 + (Ifges(4,2) * t336 + Ifges(3,3) + (Ifges(4,1) * t340 + Ifges(4,4) * t378) * t340) * qJDD(2); t267; mrSges(5,1) * t280 - mrSges(5,2) * t281 - mrSges(6,1) * t277 + mrSges(6,3) * t276 - pkin(4) * t274 + qJ(5) * t358 + t367 * t318 + (-qJ(5) * t299 + t366) * t317 + t373 * t305 + (-qJ(5) * mrSges(6,2) - t372) * t304 + t379 * qJDD(4); t274;];
tauJ = t1;
