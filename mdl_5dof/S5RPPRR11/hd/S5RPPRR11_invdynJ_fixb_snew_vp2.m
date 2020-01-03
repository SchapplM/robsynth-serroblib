% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRR11
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRR11_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:36
% EndTime: 2019-12-31 18:05:37
% DurationCPUTime: 0.47s
% Computational Cost: add. (1693->158), mult. (3067->196), div. (0->0), fcn. (1428->6), ass. (0->63)
t310 = qJD(1) ^ 2;
t305 = sin(qJ(1));
t308 = cos(qJ(1));
t327 = t305 * g(1) - t308 * g(2);
t280 = -qJDD(1) * pkin(1) - t310 * qJ(2) + qJDD(2) - t327;
t275 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t280;
t272 = -t310 * pkin(6) - t275;
t304 = sin(qJ(4));
t307 = cos(qJ(4));
t324 = qJD(1) * qJD(4);
t320 = t307 * t324;
t289 = -t304 * qJDD(1) - t320;
t321 = t304 * t324;
t290 = t307 * qJDD(1) - t321;
t258 = (-t290 + t321) * pkin(7) + (-t289 + t320) * pkin(4) + t272;
t318 = -t308 * g(1) - t305 * g(2);
t330 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t318;
t276 = qJDD(3) + (-pkin(1) - qJ(3)) * t310 + t330;
t273 = -qJDD(1) * pkin(6) + t276;
t269 = -t307 * g(3) + t304 * t273;
t288 = (t304 * pkin(4) - t307 * pkin(7)) * qJD(1);
t309 = qJD(4) ^ 2;
t325 = t304 * qJD(1);
t260 = -t309 * pkin(4) + qJDD(4) * pkin(7) - t288 * t325 + t269;
t303 = sin(qJ(5));
t306 = cos(qJ(5));
t256 = t306 * t258 - t303 * t260;
t326 = qJD(1) * t307;
t285 = t306 * qJD(4) - t303 * t326;
t267 = t285 * qJD(5) + t303 * qJDD(4) + t306 * t290;
t286 = t303 * qJD(4) + t306 * t326;
t271 = -t285 * mrSges(6,1) + t286 * mrSges(6,2);
t293 = qJD(5) + t325;
t277 = -t293 * mrSges(6,2) + t285 * mrSges(6,3);
t284 = qJDD(5) - t289;
t254 = m(6) * t256 + t284 * mrSges(6,1) - t267 * mrSges(6,3) - t286 * t271 + t293 * t277;
t257 = t303 * t258 + t306 * t260;
t266 = -t286 * qJD(5) + t306 * qJDD(4) - t303 * t290;
t278 = t293 * mrSges(6,1) - t286 * mrSges(6,3);
t255 = m(6) * t257 - t284 * mrSges(6,2) + t266 * mrSges(6,3) + t285 * t271 - t293 * t278;
t247 = t306 * t254 + t303 * t255;
t291 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t325;
t292 = (qJD(4) * mrSges(5,1)) - mrSges(5,3) * t326;
t331 = m(4) * t275 - m(5) * t272 + t289 * mrSges(5,1) - t290 * mrSges(5,2) - t247 + (-t291 * t304 - t292 * t307) * qJD(1);
t329 = t304 * g(3);
t268 = t307 * t273 + t329;
t287 = (t304 * mrSges(5,1) + t307 * mrSges(5,2)) * qJD(1);
t259 = -qJDD(4) * pkin(4) - t309 * pkin(7) - t329 + (qJD(1) * t288 - t273) * t307;
t313 = -m(6) * t259 + t266 * mrSges(6,1) - t267 * mrSges(6,2) + t285 * t277 - t286 * t278;
t319 = -t303 * t254 + t306 * t255;
t328 = t304 * (m(5) * t269 - qJDD(4) * mrSges(5,2) + t289 * mrSges(5,3) - qJD(4) * t292 - t287 * t325 + t319) + t307 * (m(5) * t268 + qJDD(4) * mrSges(5,1) - t290 * mrSges(5,3) + qJD(4) * t291 - t287 * t326 + t313);
t315 = m(4) * t276 + qJDD(1) * mrSges(4,2) - t310 * mrSges(4,3) + t328;
t262 = Ifges(6,4) * t286 + Ifges(6,2) * t285 + Ifges(6,6) * t293;
t263 = Ifges(6,1) * t286 + Ifges(6,4) * t285 + Ifges(6,5) * t293;
t311 = mrSges(6,1) * t256 - mrSges(6,2) * t257 + Ifges(6,5) * t267 + Ifges(6,6) * t266 + Ifges(6,3) * t284 + t286 * t262 - t285 * t263;
t283 = (Ifges(5,5) * qJD(4)) + (t307 * Ifges(5,1) - t304 * Ifges(5,4)) * qJD(1);
t282 = (Ifges(5,6) * qJD(4)) + (t307 * Ifges(5,4) - t304 * Ifges(5,2)) * qJD(1);
t279 = t310 * pkin(1) - t330;
t261 = Ifges(6,5) * t286 + Ifges(6,6) * t285 + Ifges(6,3) * t293;
t249 = mrSges(6,2) * t259 - mrSges(6,3) * t256 + Ifges(6,1) * t267 + Ifges(6,4) * t266 + Ifges(6,5) * t284 + t285 * t261 - t293 * t262;
t248 = -mrSges(6,1) * t259 + mrSges(6,3) * t257 + Ifges(6,4) * t267 + Ifges(6,2) * t266 + Ifges(6,6) * t284 - t286 * t261 + t293 * t263;
t245 = m(3) * t280 + (-mrSges(4,2) - mrSges(3,3)) * t310 + (mrSges(3,2) - mrSges(4,3)) * qJDD(1) + t331;
t1 = [qJ(2) * (-m(3) * t279 + t310 * mrSges(3,2) + t315) - pkin(1) * t245 + mrSges(2,1) * t327 - mrSges(2,2) * t318 + mrSges(3,2) * t280 - mrSges(3,3) * t279 + mrSges(4,2) * t276 - mrSges(4,3) * t275 + t307 * (mrSges(5,2) * t272 - mrSges(5,3) * t268 + Ifges(5,1) * t290 + Ifges(5,4) * t289 + Ifges(5,5) * qJDD(4) - pkin(7) * t247 - qJD(4) * t282 - t303 * t248 + t306 * t249) - t304 * (-mrSges(5,1) * t272 + mrSges(5,3) * t269 + Ifges(5,4) * t290 + Ifges(5,2) * t289 + Ifges(5,6) * qJDD(4) - pkin(4) * t247 + qJD(4) * t283 - t311) - pkin(6) * t328 + (qJ(2) * mrSges(3,3) + Ifges(3,1) + Ifges(4,1) + Ifges(2,3)) * qJDD(1) + (t310 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t331) * qJ(3); t245; t315; Ifges(5,5) * t290 + Ifges(5,6) * t289 + Ifges(5,3) * qJDD(4) + mrSges(5,1) * t268 - mrSges(5,2) * t269 + t303 * t249 + t306 * t248 + pkin(4) * t313 + pkin(7) * t319 + (t307 * t282 + t304 * t283) * qJD(1); t311;];
tauJ = t1;
