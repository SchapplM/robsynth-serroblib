% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPPRP3
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
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
% Datum: 2019-12-31 17:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPPRP3_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP3_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:50:48
% EndTime: 2019-12-31 17:50:49
% DurationCPUTime: 0.80s
% Computational Cost: add. (992->141), mult. (1777->165), div. (0->0), fcn. (758->6), ass. (0->61)
t300 = sin(qJ(4));
t302 = cos(qJ(4));
t326 = Ifges(5,4) + Ifges(6,4);
t334 = t302 * (Ifges(5,1) + Ifges(6,1)) - t300 * t326;
t333 = t302 * t326 - t300 * (Ifges(5,2) + Ifges(6,2));
t325 = Ifges(5,5) + Ifges(6,5);
t324 = Ifges(5,6) + Ifges(6,6);
t330 = (t333 * qJD(1) + t324 * qJD(4)) * t302;
t301 = sin(qJ(1));
t303 = cos(qJ(1));
t313 = t301 * g(1) - g(2) * t303;
t276 = qJDD(1) * pkin(1) + t313;
t304 = qJD(1) ^ 2;
t309 = -g(1) * t303 - g(2) * t301;
t279 = -pkin(1) * t304 + t309;
t298 = sin(pkin(7));
t299 = cos(pkin(7));
t259 = t298 * t276 + t299 * t279;
t310 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t259;
t329 = -pkin(2) - pkin(6);
t328 = pkin(4) * t304;
t258 = t276 * t299 - t298 * t279;
t307 = -qJ(3) * t304 + qJDD(3) - t258;
t255 = t329 * qJDD(1) + t307;
t295 = -g(3) + qJDD(2);
t250 = t300 * t255 + t302 * t295;
t318 = qJD(1) * qJD(4);
t280 = -qJDD(1) * t300 - t302 * t318;
t319 = t302 * qJD(1);
t285 = qJD(4) * pkin(4) - qJ(5) * t319;
t294 = t300 ^ 2;
t316 = -0.2e1 * qJD(1) * qJD(5);
t246 = qJ(5) * t280 - qJD(4) * t285 - t294 * t328 + t300 * t316 + t250;
t323 = m(6) * t246 + t280 * mrSges(6,3);
t321 = t334 * qJD(1) + t325 * qJD(4);
t320 = qJD(1) * t300;
t252 = t302 * t255;
t281 = qJDD(1) * t302 - t300 * t318;
t245 = t302 * t316 + qJDD(4) * pkin(4) - qJ(5) * t281 + t252 + (-qJ(5) * t318 - t302 * t328 - t295) * t300;
t283 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t320;
t315 = m(6) * t245 + qJDD(4) * mrSges(6,1) + qJD(4) * t283;
t277 = (t300 * mrSges(6,1) + t302 * mrSges(6,2)) * qJD(1);
t312 = qJD(1) * (-t277 - (t300 * mrSges(5,1) + t302 * mrSges(5,2)) * qJD(1));
t248 = t285 * t319 - pkin(4) * t280 + qJDD(5) + (-qJ(5) * t294 + t329) * t304 + t310;
t286 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t319;
t311 = m(6) * t248 + t281 * mrSges(6,2) + t283 * t320 + t286 * t319;
t249 = -t295 * t300 + t252;
t284 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t320;
t239 = m(5) * t249 + qJDD(4) * mrSges(5,1) + qJD(4) * t284 + (-mrSges(5,3) - mrSges(6,3)) * t281 + t302 * t312 + t315;
t287 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t319;
t240 = m(5) * t250 + mrSges(5,3) * t280 + (-mrSges(5,2) - mrSges(6,2)) * qJDD(4) + (-t286 - t287) * qJD(4) + t300 * t312 + t323;
t308 = t239 * t302 + t240 * t300;
t257 = -qJDD(1) * pkin(2) + t307;
t306 = m(4) * t257 - t304 * mrSges(4,3) + t308;
t254 = t329 * t304 + t310;
t256 = pkin(2) * t304 - t310;
t305 = -m(4) * t256 + m(5) * t254 + t304 * mrSges(4,2) + t281 * mrSges(5,2) + qJDD(1) * mrSges(4,3) + t284 * t320 + t287 * t319 + t311;
t242 = -mrSges(6,1) * t280 + t311;
t241 = -mrSges(6,3) * t281 - t277 * t319 + t315;
t238 = qJDD(1) * mrSges(4,2) + t306;
t1 = [pkin(1) * (t298 * (m(3) * t259 - mrSges(3,1) * t304 + t305) + t299 * (m(3) * t258 - mrSges(3,2) * t304 - t306)) + mrSges(2,1) * t313 - mrSges(2,2) * t309 - pkin(2) * t238 + qJ(3) * t305 + mrSges(4,2) * t257 - mrSges(4,3) * t256 + t302 * (mrSges(5,2) * t254 + mrSges(6,2) * t248 - mrSges(5,3) * t249 - mrSges(6,3) * t245 - qJ(5) * t241) - t300 * (-mrSges(5,1) * t254 + mrSges(5,3) * t250 - mrSges(6,1) * t248 + mrSges(6,3) * t246 - pkin(4) * t242 + qJ(5) * (-t277 * t320 + t323)) - pkin(6) * t308 + mrSges(3,1) * t258 - mrSges(3,2) * t259 + t334 * t281 + (t302 * t325 - t300 * (-qJ(5) * mrSges(6,2) + t324)) * qJDD(4) + (-t330 - t300 * (-qJ(5) * t286 + t321)) * qJD(4) + ((pkin(1) * t298 + qJ(3)) * (-mrSges(5,1) - mrSges(6,1)) + t333) * t280 + (pkin(1) * (-t298 * mrSges(3,2) + t299 * (mrSges(3,1) - mrSges(4,2))) + Ifges(2,3) + Ifges(4,1) + Ifges(3,3)) * qJDD(1); -t239 * t300 + t240 * t302 + (m(3) + m(4)) * t295; t238; mrSges(5,1) * t249 + mrSges(6,1) * t245 - mrSges(5,2) * t250 - mrSges(6,2) * t246 + pkin(4) * t241 + t325 * t281 + t324 * t280 + (Ifges(5,3) + Ifges(6,3)) * qJDD(4) + (t321 * t300 + t330) * qJD(1); t242;];
tauJ = t1;
