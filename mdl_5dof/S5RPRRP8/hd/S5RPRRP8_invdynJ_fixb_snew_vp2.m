% Calculate vector of inverse dynamics joint torques for with Newton-Euler
% S5RPRRP8
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
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
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ = S5RPRRP8_invdynJ_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP8_invdynJ_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_tauJ_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:47:13
% EndTime: 2019-12-31 18:47:14
% DurationCPUTime: 0.52s
% Computational Cost: add. (2432->143), mult. (3093->175), div. (0->0), fcn. (1109->6), ass. (0->64)
t310 = Ifges(5,1) + Ifges(6,1);
t304 = Ifges(5,4) - Ifges(6,5);
t303 = Ifges(5,5) + Ifges(6,4);
t309 = Ifges(5,2) + Ifges(6,3);
t302 = Ifges(5,6) - Ifges(6,6);
t308 = Ifges(5,3) + Ifges(6,2);
t307 = -pkin(1) - pkin(2);
t283 = cos(qJ(4));
t306 = t283 * g(3);
t305 = mrSges(5,3) + mrSges(6,2);
t276 = -qJD(1) + qJD(3);
t280 = sin(qJ(4));
t301 = t276 * t280;
t300 = t276 * t283;
t287 = qJD(1) ^ 2;
t282 = sin(qJ(1));
t285 = cos(qJ(1));
t292 = -t285 * g(1) - t282 * g(2);
t290 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t292;
t244 = t307 * t287 + t290;
t295 = t282 * g(1) - t285 * g(2);
t289 = -t287 * qJ(2) + qJDD(2) - t295;
t245 = t307 * qJDD(1) + t289;
t281 = sin(qJ(3));
t284 = cos(qJ(3));
t239 = t284 * t244 + t281 * t245;
t274 = t276 ^ 2;
t275 = -qJDD(1) + qJDD(3);
t237 = -t274 * pkin(3) + t275 * pkin(7) + t239;
t234 = t280 * g(3) + t283 * t237;
t299 = (-t280 * t304 - t309 * t283) * t276 - t302 * qJD(4);
t298 = (t280 * t303 + t283 * t302) * t276 + t308 * qJD(4);
t297 = (t310 * t280 + t283 * t304) * t276 + t303 * qJD(4);
t296 = qJD(4) * t276;
t233 = -t280 * t237 + t306;
t260 = (-mrSges(6,1) * t283 - mrSges(6,3) * t280) * t276;
t261 = (-mrSges(5,1) * t283 + mrSges(5,2) * t280) * t276;
t262 = t280 * t275 + t283 * t296;
t263 = t283 * t275 - t280 * t296;
t266 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t301;
t268 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t300;
t259 = (-pkin(4) * t283 - qJ(5) * t280) * t276;
t286 = qJD(4) ^ 2;
t232 = -qJDD(4) * pkin(4) - t306 - t286 * qJ(5) + qJDD(5) + (t259 * t276 + t237) * t280;
t269 = mrSges(6,2) * t300 + qJD(4) * mrSges(6,3);
t293 = -m(6) * t232 + qJDD(4) * mrSges(6,1) + qJD(4) * t269;
t231 = -t286 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t259 * t300 + t234;
t267 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t301;
t294 = m(6) * t231 + qJDD(4) * mrSges(6,3) + qJD(4) * t267 + t260 * t300;
t222 = -t280 * (m(5) * t233 + qJDD(4) * mrSges(5,1) + qJD(4) * t268 + (-t260 - t261) * t301 - t305 * t262 + t293) + t283 * (m(5) * t234 - qJDD(4) * mrSges(5,2) - qJD(4) * t266 + t261 * t300 + t305 * t263 + t294);
t238 = -t281 * t244 + t284 * t245;
t221 = m(4) * t239 - t274 * mrSges(4,1) - t275 * mrSges(4,2) + t222;
t236 = -t275 * pkin(3) - t274 * pkin(7) - t238;
t229 = -t263 * pkin(4) - t262 * qJ(5) + (-0.2e1 * qJD(5) * t280 + (pkin(4) * t280 - qJ(5) * t283) * qJD(4)) * t276 + t236;
t227 = m(6) * t229 - t263 * mrSges(6,1) - t262 * mrSges(6,3) - t267 * t301 - t269 * t300;
t224 = -m(5) * t236 + t263 * mrSges(5,1) - t262 * mrSges(5,2) - t266 * t301 + t268 * t300 - t227;
t223 = m(4) * t238 + t275 * mrSges(4,1) - t274 * mrSges(4,2) + t224;
t291 = t281 * t221 + t284 * t223;
t288 = mrSges(4,1) * t238 - mrSges(4,2) * t239 + Ifges(4,3) * t275 + pkin(3) * t224 + pkin(7) * t222 + t283 * (-mrSges(5,1) * t236 - mrSges(6,1) * t229 + mrSges(6,2) * t231 + mrSges(5,3) * t234 - pkin(4) * t227 + t297 * qJD(4) + t302 * qJDD(4) + t304 * t262 + t309 * t263 - t298 * t301) + t280 * (mrSges(5,2) * t236 + mrSges(6,2) * t232 - mrSges(5,3) * t233 - mrSges(6,3) * t229 - qJ(5) * t227 + t299 * qJD(4) + t303 * qJDD(4) + t310 * t262 + t304 * t263 + t298 * t300);
t258 = -qJDD(1) * pkin(1) + t289;
t252 = -t287 * pkin(1) + t290;
t228 = t262 * mrSges(6,2) + t260 * t301 - t293;
t218 = m(3) * t258 - qJDD(1) * mrSges(3,1) - t287 * mrSges(3,3) + t291;
t1 = [-pkin(1) * t218 + qJ(2) * (m(3) * t252 - t287 * mrSges(3,1) + t284 * t221 - t281 * t223) + mrSges(2,1) * t295 - mrSges(2,2) * t292 - pkin(2) * t291 - mrSges(3,1) * t258 + mrSges(3,3) * t252 + (qJ(2) * mrSges(3,3) + Ifges(3,2) + Ifges(2,3)) * qJDD(1) - t288; t218; t288; mrSges(5,1) * t233 - mrSges(5,2) * t234 - mrSges(6,1) * t232 + mrSges(6,3) * t231 - pkin(4) * t228 + qJ(5) * t294 + (qJ(5) * mrSges(6,2) + t302) * t263 + t303 * t262 + t308 * qJDD(4) + (-t299 * t280 - t297 * t283) * t276; t228;];
tauJ = t1;
