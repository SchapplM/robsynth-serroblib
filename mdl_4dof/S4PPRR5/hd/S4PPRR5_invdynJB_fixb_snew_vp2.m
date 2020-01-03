% Calculate vector of inverse dynamics joint torques and base forces with Newton-Euler
% S4PPRR5
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% tauJB [(6+4)x1]
%   joint torques and base forces of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJB = S4PPRR5_invdynJB_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR5_invdynJB_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJB_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:19:44
% EndTime: 2019-12-31 16:19:44
% DurationCPUTime: 0.40s
% Computational Cost: add. (2237->124), mult. (3661->159), div. (0->0), fcn. (1922->6), ass. (0->55)
t279 = sin(pkin(6));
t280 = cos(pkin(6));
t270 = t279 * g(1) - t280 * g(2);
t268 = qJDD(2) - t270;
t278 = -g(3) + qJDD(1);
t282 = sin(qJ(3));
t284 = cos(qJ(3));
t258 = t282 * t268 + t284 * t278;
t285 = qJD(3) ^ 2;
t255 = -t285 * pkin(3) + qJDD(3) * pkin(5) + t258;
t271 = t280 * g(1) + t279 * g(2);
t281 = sin(qJ(4));
t283 = cos(qJ(4));
t252 = -t281 * t255 - t283 * t271;
t253 = t283 * t255 - t281 * t271;
t260 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t281 + Ifges(5,2) * t283) * qJD(3);
t261 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t281 + Ifges(5,4) * t283) * qJD(3);
t295 = qJD(3) * qJD(4);
t266 = t281 * qJDD(3) + t283 * t295;
t267 = t283 * qJDD(3) - t281 * t295;
t300 = mrSges(5,1) * t252 - mrSges(5,2) * t253 + Ifges(5,5) * t266 + Ifges(5,6) * t267 + Ifges(5,3) * qJDD(4) + (t260 * t281 - t261 * t283) * qJD(3);
t299 = mrSges(2,2) - mrSges(3,3);
t265 = (-mrSges(5,1) * t283 + mrSges(5,2) * t281) * qJD(3);
t296 = qJD(3) * t283;
t273 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t296;
t297 = qJD(3) * t281;
t249 = m(5) * t252 + qJDD(4) * mrSges(5,1) - t266 * mrSges(5,3) + qJD(4) * t273 - t265 * t297;
t272 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t297;
t250 = m(5) * t253 - qJDD(4) * mrSges(5,2) + t267 * mrSges(5,3) - qJD(4) * t272 + t265 * t296;
t291 = -t281 * t249 + t283 * t250;
t238 = m(4) * t258 - t285 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t291;
t257 = t284 * t268 - t282 * t278;
t254 = -qJDD(3) * pkin(3) - t285 * pkin(5) - t257;
t287 = -m(5) * t254 + t267 * mrSges(5,1) - t266 * mrSges(5,2) - t272 * t297 + t273 * t296;
t245 = m(4) * t257 + qJDD(3) * mrSges(4,1) - t285 * mrSges(4,2) + t287;
t233 = t282 * t238 + t284 * t245;
t231 = m(3) * t268 + t233;
t230 = m(2) * t270 - t231;
t240 = t283 * t249 + t281 * t250;
t294 = m(4) * t271 - t240;
t236 = (-m(2) - m(3)) * t271 - t294;
t298 = t280 * t230 + t279 * t236;
t293 = -t279 * t230 + t280 * t236;
t292 = t284 * t238 - t282 * t245;
t232 = m(3) * t278 + t292;
t289 = m(2) * t278 + t232;
t259 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t281 + Ifges(5,6) * t283) * qJD(3);
t243 = -mrSges(5,1) * t254 + mrSges(5,3) * t253 + Ifges(5,4) * t266 + Ifges(5,2) * t267 + Ifges(5,6) * qJDD(4) + qJD(4) * t261 - t259 * t297;
t244 = mrSges(5,2) * t254 - mrSges(5,3) * t252 + Ifges(5,1) * t266 + Ifges(5,4) * t267 + Ifges(5,5) * qJDD(4) - qJD(4) * t260 + t259 * t296;
t288 = mrSges(4,1) * t257 - mrSges(4,2) * t258 + Ifges(4,3) * qJDD(3) + pkin(3) * t287 + pkin(5) * t291 + t283 * t243 + t281 * t244;
t228 = mrSges(4,1) * t271 + mrSges(4,3) * t258 + t285 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t240 - t300;
t227 = -mrSges(4,2) * t271 - mrSges(4,3) * t257 + Ifges(4,5) * qJDD(3) - t285 * Ifges(4,6) - pkin(5) * t240 - t281 * t243 + t283 * t244;
t226 = mrSges(3,1) * t268 - mrSges(2,3) * t270 + pkin(2) * t233 - qJ(2) * t232 + t299 * t278 + t288;
t225 = -t282 * t227 - t284 * t228 - pkin(2) * t294 - pkin(4) * t292 - pkin(1) * t232 + (-mrSges(2,1) + mrSges(3,2)) * t278 + (-mrSges(3,1) - mrSges(2,3)) * t271;
t1 = [-m(1) * g(1) + t293; -m(1) * g(2) + t298; -m(1) * g(3) + t289; -mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * t298 - t279 * t225 + t280 * t226; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + qJ(1) * t293 + t280 * t225 + t279 * t226; -pkin(1) * t231 - qJ(2) * t294 + t284 * t227 - t282 * t228 - pkin(4) * t233 + mrSges(2,1) * t270 + mrSges(3,2) * t268 - mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + (-qJ(2) * m(3) + t299) * t271; t289; t231; t288; t300;];
tauJB = t1;
