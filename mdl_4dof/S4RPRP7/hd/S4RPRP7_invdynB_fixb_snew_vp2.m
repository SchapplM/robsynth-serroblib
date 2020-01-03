% Calculate vector of inverse dynamics base forces with Newton-Euler for
% S4RPRP7
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB = S4RPRP7_invdynB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP7_invdynB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:47:06
% EndTime: 2019-12-31 16:47:06
% DurationCPUTime: 0.56s
% Computational Cost: add. (2096->167), mult. (3845->198), div. (0->0), fcn. (1465->4), ass. (0->66)
t355 = Ifges(4,1) + Ifges(5,1);
t347 = Ifges(4,4) - Ifges(5,5);
t346 = (Ifges(5,4) + Ifges(4,5));
t354 = Ifges(4,2) + Ifges(5,3);
t343 = (Ifges(4,6) - Ifges(5,6));
t353 = (Ifges(4,3) + Ifges(5,2));
t320 = sin(qJ(1));
t322 = cos(qJ(1));
t310 = -g(1) * t322 - g(2) * t320;
t352 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t310;
t351 = -pkin(1) - pkin(5);
t319 = sin(qJ(3));
t350 = g(3) * t319;
t349 = mrSges(2,1) - mrSges(3,2);
t348 = -mrSges(4,3) - mrSges(5,2);
t345 = Ifges(2,5) - Ifges(3,4);
t344 = (Ifges(2,6) - Ifges(3,5));
t309 = g(1) * t320 - t322 * g(2);
t324 = qJD(1) ^ 2;
t328 = -qJ(2) * t324 + qJDD(2) - t309;
t284 = t351 * qJDD(1) + t328;
t321 = cos(qJ(3));
t280 = -g(3) * t321 + t319 * t284;
t337 = qJD(1) * qJD(3);
t301 = qJDD(1) * t319 + t321 * t337;
t338 = qJD(1) * t321;
t306 = (qJD(3) * mrSges(4,1)) - mrSges(4,3) * t338;
t299 = (mrSges(5,1) * t319 - mrSges(5,3) * t321) * qJD(1);
t331 = qJD(1) * (-t299 - (mrSges(4,1) * t319 + mrSges(4,2) * t321) * qJD(1));
t298 = (pkin(3) * t319 - qJ(4) * t321) * qJD(1);
t323 = qJD(3) ^ 2;
t339 = qJD(1) * t319;
t277 = -pkin(3) * t323 + qJDD(3) * qJ(4) + (2 * qJD(4) * qJD(3)) - t298 * t339 + t280;
t307 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t338;
t335 = m(5) * t277 + qJDD(3) * mrSges(5,3) + qJD(3) * t307;
t271 = m(4) * t280 - qJDD(3) * mrSges(4,2) - qJD(3) * t306 + t348 * t301 + t319 * t331 + t335;
t279 = t284 * t321 + t350;
t302 = qJDD(1) * t321 - t319 * t337;
t305 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t339;
t278 = -qJDD(3) * pkin(3) - t350 - qJ(4) * t323 + qJDD(4) + (qJD(1) * t298 - t284) * t321;
t308 = -mrSges(5,2) * t339 + qJD(3) * mrSges(5,3);
t330 = -m(5) * t278 + qJDD(3) * mrSges(5,1) + qJD(3) * t308;
t272 = m(4) * t279 + qJDD(3) * mrSges(4,1) + qJD(3) * t305 + t348 * t302 + t321 * t331 + t330;
t266 = t271 * t319 + t272 * t321;
t286 = -qJDD(1) * pkin(1) + t328;
t327 = -m(3) * t286 + (t324 * mrSges(3,3)) - t266;
t264 = m(2) * t309 - (mrSges(2,2) * t324) + t349 * qJDD(1) + t327;
t285 = pkin(1) * t324 + t352;
t283 = t351 * t324 - t352;
t275 = pkin(3) * t301 - qJ(4) * t302 + (-0.2e1 * qJD(4) * t321 + (pkin(3) * t321 + qJ(4) * t319) * qJD(3)) * qJD(1) + t283;
t273 = m(5) * t275 + t301 * mrSges(5,1) - mrSges(5,3) * t302 - t307 * t338 + t308 * t339;
t326 = -m(4) * t283 - t301 * mrSges(4,1) - t302 * mrSges(4,2) - t305 * t339 - t306 * t338 - t273;
t325 = -m(3) * t285 + (t324 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) - t326;
t269 = m(2) * t310 - (mrSges(2,1) * t324) - qJDD(1) * mrSges(2,2) + t325;
t342 = t322 * t264 + t320 * t269;
t341 = -(t343 * qJD(3)) + (t354 * t319 - t347 * t321) * qJD(1);
t340 = (t346 * qJD(3)) + (-t347 * t319 + t355 * t321) * qJD(1);
t334 = -t264 * t320 + t322 * t269;
t333 = t321 * t271 - t319 * t272;
t332 = qJD(1) * (-(t353 * qJD(3)) + (t343 * t319 - t346 * t321) * qJD(1));
t265 = -m(3) * g(3) + t333;
t262 = mrSges(4,2) * t283 + mrSges(5,2) * t278 - mrSges(4,3) * t279 - mrSges(5,3) * t275 - qJ(4) * t273 + t341 * qJD(3) + t346 * qJDD(3) - t347 * t301 + t355 * t302 + t319 * t332;
t261 = -mrSges(4,1) * t283 - mrSges(5,1) * t275 + mrSges(5,2) * t277 + mrSges(4,3) * t280 - pkin(3) * t273 + t340 * qJD(3) + t343 * qJDD(3) - t354 * t301 + t347 * t302 + t321 * t332;
t260 = qJ(4) * t335 + pkin(3) * t330 - mrSges(2,3) * t309 + mrSges(5,3) * t277 - mrSges(5,1) * t278 + mrSges(4,1) * t279 - mrSges(4,2) * t280 + mrSges(3,1) * t286 - qJ(2) * t265 + pkin(2) * t266 - (t344 * t324) + (-pkin(3) * mrSges(5,2) + t346) * t302 + (-mrSges(5,2) * qJ(4) - t343) * t301 + t353 * qJDD(3) + t345 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) + ((-pkin(3) * t299 - t341) * t321 + (-qJ(4) * t299 + t340) * t319) * qJD(1);
t259 = -mrSges(3,1) * t285 + mrSges(2,3) * t310 - pkin(1) * t265 - pkin(2) * t326 - pkin(5) * t333 + t349 * g(3) + t344 * qJDD(1) - t321 * t261 - t319 * t262 + t345 * t324;
t1 = [-m(1) * g(1) + t334; -m(1) * g(2) + t342; (-m(1) - m(2) - m(3)) * g(3) + t333; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - pkin(4) * t342 - t320 * t259 + t322 * t260; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + pkin(4) * t334 + t322 * t259 + t320 * t260; pkin(1) * t327 + qJ(2) * t325 + t321 * t262 - t319 * t261 - pkin(5) * t266 + mrSges(2,1) * t309 - mrSges(2,2) * t310 + mrSges(3,2) * t286 - mrSges(3,3) * t285 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-mrSges(3,2) * pkin(1) + Ifges(3,1) + Ifges(2,3)) * qJDD(1);];
tauB = t1;
