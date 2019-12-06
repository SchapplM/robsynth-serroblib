% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5PRPRP2
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5PRPRP2_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:42
% EndTime: 2019-12-05 15:30:44
% DurationCPUTime: 1.03s
% Computational Cost: add. (1702->148), mult. (3652->195), div. (0->0), fcn. (2194->8), ass. (0->77)
t333 = sin(qJ(4));
t335 = cos(qJ(4));
t367 = Ifges(5,4) + Ifges(6,4);
t375 = t335 * (Ifges(5,1) + Ifges(6,1)) - t333 * t367;
t374 = t335 * t367 - t333 * (Ifges(5,2) + Ifges(6,2));
t366 = Ifges(5,5) + Ifges(6,5);
t365 = Ifges(5,6) + Ifges(6,6);
t331 = cos(pkin(8));
t356 = qJD(2) * t331;
t322 = qJD(4) - t356;
t329 = sin(pkin(8));
t357 = qJD(2) * t329;
t346 = (t365 * t322 + t374 * t357) * t335;
t337 = qJD(2) ^ 2;
t330 = sin(pkin(7));
t332 = cos(pkin(7));
t318 = g(1) * t330 - g(2) * t332;
t319 = -g(1) * t332 - g(2) * t330;
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t358 = t334 * t318 + t336 * t319;
t293 = -pkin(2) * t337 + qJDD(2) * qJ(3) + t358;
t328 = -g(3) + qJDD(1);
t354 = qJD(2) * qJD(3);
t285 = t328 * t331 + (-t293 - 0.2e1 * t354) * t329;
t307 = (t333 * mrSges(6,1) + t335 * mrSges(6,2)) * t357;
t353 = qJD(2) * qJD(4);
t310 = (qJDD(2) * t335 - t333 * t353) * t329;
t348 = t335 * t357;
t370 = 0.2e1 * t331;
t286 = t331 * t293 + t329 * t328 + t354 * t370;
t342 = -pkin(3) * t331 - pkin(6) * t329;
t317 = t342 * qJD(2);
t284 = t317 * t356 + t286;
t344 = t318 * t336 - t334 * t319;
t339 = -qJ(3) * t337 + qJDD(3) - t344;
t289 = (-pkin(2) + t342) * qJDD(2) + t339;
t288 = t335 * t289;
t352 = t331 * qJDD(2);
t321 = qJDD(4) - t352;
t343 = -0.2e1 * qJD(5) * t357;
t326 = t329 ^ 2;
t363 = t326 * t337;
t276 = t335 * t343 + pkin(4) * t321 - qJ(5) * t310 + t288 + (-pkin(4) * t335 * t363 - qJ(5) * t322 * t357 - t284) * t333;
t349 = t333 * t357;
t302 = -mrSges(6,2) * t322 - mrSges(6,3) * t349;
t350 = m(6) * t276 + t321 * mrSges(6,1) + t322 * t302;
t273 = -mrSges(6,3) * t310 - t307 * t348 + t350;
t281 = t335 * t284 + t333 * t289;
t304 = pkin(4) * t322 - qJ(5) * t348;
t309 = (-qJDD(2) * t333 - t335 * t353) * t329;
t351 = t333 ^ 2 * t363;
t278 = -pkin(4) * t351 + qJ(5) * t309 - t304 * t322 + t333 * t343 + t281;
t280 = -t284 * t333 + t288;
t371 = mrSges(5,1) * t280 + mrSges(6,1) * t276 - mrSges(5,2) * t281 - mrSges(6,2) * t278 + pkin(4) * t273 + t365 * t309 + t366 * t310 + (Ifges(5,3) + Ifges(6,3)) * t321;
t283 = t317 * t357 - t285;
t279 = -pkin(4) * t309 - qJ(5) * t351 + t304 * t348 + qJDD(5) + t283;
t369 = m(6) * t279;
t368 = -mrSges(5,2) - mrSges(6,2);
t362 = m(6) * t278 + t309 * mrSges(6,3);
t360 = t366 * t322 + t375 * t357;
t305 = mrSges(6,1) * t322 - mrSges(6,3) * t348;
t359 = -mrSges(5,1) * t322 + mrSges(5,3) * t348 - t305;
t355 = qJDD(2) * mrSges(4,3);
t341 = (-t307 - (t333 * mrSges(5,1) + t335 * mrSges(5,2)) * t357) * t357;
t340 = -mrSges(4,1) * t331 + mrSges(4,2) * t329;
t303 = -mrSges(5,2) * t322 - mrSges(5,3) * t349;
t270 = m(5) * t280 + mrSges(5,1) * t321 + t303 * t322 + (-mrSges(5,3) - mrSges(6,3)) * t310 + t335 * t341 + t350;
t271 = m(5) * t281 + mrSges(5,3) * t309 + t368 * t321 + t359 * t322 + t333 * t341 + t362;
t269 = t270 * t335 + t271 * t333;
t315 = t340 * qJD(2);
t292 = -qJDD(2) * pkin(2) + t339;
t274 = t369 - mrSges(6,1) * t309 + mrSges(6,2) * t310 + (t302 * t333 + t305 * t335) * t357;
t272 = m(4) * t285 - m(5) * t283 - t369 + t368 * t310 + (mrSges(5,1) + mrSges(6,1)) * t309 + (-t355 + (-t315 + t359 * t335 + (-t302 - t303) * t333) * qJD(2)) * t329;
t268 = m(4) * t292 + t340 * qJDD(2) + (-t331 ^ 2 - t326) * t337 * mrSges(4,3) + t269;
t267 = m(4) * t286 - t270 * t333 + t271 * t335 + (qJD(2) * t315 + t355) * t331;
t1 = [t267 * t329 + t272 * t331 + (m(2) + m(3)) * t328; Ifges(3,3) * qJDD(2) + mrSges(3,1) * t344 - mrSges(3,2) * t358 - pkin(2) * t268 + (-mrSges(4,1) * t292 + mrSges(4,3) * t286 + Ifges(4,2) * t352 - pkin(3) * t269 + qJ(3) * t267 - t371) * t331 + (mrSges(4,2) * t292 - mrSges(4,3) * t285 + t335 * (mrSges(5,2) * t283 + mrSges(6,2) * t279 - mrSges(5,3) * t280 - mrSges(6,3) * t276 - qJ(5) * t273) - t333 * (-mrSges(5,1) * t283 - mrSges(6,1) * t279 + mrSges(5,3) * t281 + mrSges(6,3) * t278 - pkin(4) * t274 + qJ(5) * t362) - pkin(6) * t269 - qJ(3) * t272 + (-t346 - t333 * (-qJ(5) * t305 + t360)) * t322 + (t335 * t366 - t333 * (-qJ(5) * mrSges(6,2) + t365)) * t321 + t375 * t310 + t374 * t309 + (Ifges(4,1) * t329 + Ifges(4,4) * t370) * qJDD(2) + (-t331 * t346 + (qJ(5) * t307 * t329 * t333 - t331 * t360) * t333) * qJD(2)) * t329; t268; (t360 * t333 + t346) * t357 + t371; t274;];
tauJ = t1;
