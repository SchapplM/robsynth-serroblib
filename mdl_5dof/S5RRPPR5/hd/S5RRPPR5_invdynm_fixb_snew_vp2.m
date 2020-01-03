% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% m [3x6]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:28:59
% EndTime: 2019-12-31 19:29:09
% DurationCPUTime: 4.94s
% Computational Cost: add. (53113->310), mult. (123710->378), div. (0->0), fcn. (79328->8), ass. (0->117)
t319 = -2 * qJD(3);
t278 = sin(qJ(2));
t281 = cos(qJ(2));
t304 = qJD(1) * qJD(2);
t254 = qJDD(1) * t281 - t278 * t304;
t308 = qJD(1) * t278;
t256 = qJD(2) * pkin(2) - qJ(3) * t308;
t275 = t281 ^ 2;
t284 = qJD(1) ^ 2;
t279 = sin(qJ(1));
t282 = cos(qJ(1));
t259 = t279 * g(1) - t282 * g(2);
t299 = qJDD(1) * pkin(1) + t259;
t199 = -pkin(2) * t254 + t256 * t308 - (qJ(3) * t275 + pkin(6)) * t284 + qJDD(3) - t299;
t253 = qJDD(1) * t278 + t281 * t304;
t276 = sin(pkin(8));
t312 = cos(pkin(8));
t227 = t312 * t253 + t276 * t254;
t307 = qJD(1) * t281;
t241 = t276 * t308 - t312 * t307;
t306 = qJD(2) * t241;
t318 = t199 + (-t227 + t306) * qJ(4);
t260 = -g(1) * t282 - g(2) * t279;
t248 = -pkin(1) * t284 + qJDD(1) * pkin(6) + t260;
t311 = t248 * t278;
t314 = pkin(2) * t284;
t195 = qJDD(2) * pkin(2) - qJ(3) * t253 - t311 + (qJ(3) * t304 + t278 * t314 - g(3)) * t281;
t231 = -g(3) * t278 + t281 * t248;
t196 = qJ(3) * t254 - qJD(2) * t256 - t275 * t314 + t231;
t242 = (t276 * t281 + t312 * t278) * qJD(1);
t174 = t312 * t195 - t276 * t196 + t242 * t319;
t175 = t276 * t195 + t312 * t196 + t241 * t319;
t226 = t253 * t276 - t312 * t254;
t233 = qJD(2) * mrSges(4,1) - mrSges(4,3) * t242;
t214 = pkin(3) * t241 - qJ(4) * t242;
t283 = qJD(2) ^ 2;
t170 = -qJDD(2) * pkin(3) - t283 * qJ(4) + t242 * t214 + qJDD(4) - t174;
t162 = (-t227 - t306) * pkin(7) + (t241 * t242 - qJDD(2)) * pkin(4) + t170;
t315 = 2 * qJD(4);
t168 = -pkin(3) * t283 + qJDD(2) * qJ(4) + qJD(2) * t315 - t241 * t214 + t175;
t236 = -qJD(2) * pkin(4) - pkin(7) * t242;
t240 = t241 ^ 2;
t163 = -pkin(4) * t240 + pkin(7) * t226 + qJD(2) * t236 + t168;
t277 = sin(qJ(5));
t280 = cos(qJ(5));
t160 = t162 * t280 - t163 * t277;
t209 = t241 * t280 - t242 * t277;
t183 = qJD(5) * t209 + t226 * t277 + t227 * t280;
t210 = t241 * t277 + t242 * t280;
t189 = -mrSges(6,1) * t209 + mrSges(6,2) * t210;
t267 = -qJD(2) + qJD(5);
t200 = -mrSges(6,2) * t267 + mrSges(6,3) * t209;
t266 = -qJDD(2) + qJDD(5);
t155 = m(6) * t160 + mrSges(6,1) * t266 - mrSges(6,3) * t183 - t189 * t210 + t200 * t267;
t161 = t162 * t277 + t163 * t280;
t182 = -qJD(5) * t210 + t226 * t280 - t227 * t277;
t201 = mrSges(6,1) * t267 - mrSges(6,3) * t210;
t156 = m(6) * t161 - mrSges(6,2) * t266 + mrSges(6,3) * t182 + t189 * t209 - t201 * t267;
t147 = -t155 * t277 + t280 * t156;
t234 = -qJD(2) * mrSges(5,1) + mrSges(5,2) * t242;
t296 = m(5) * t168 + qJDD(2) * mrSges(5,3) + qJD(2) * t234 + t147;
t215 = mrSges(5,1) * t241 - mrSges(5,3) * t242;
t309 = -mrSges(4,1) * t241 - mrSges(4,2) * t242 - t215;
t313 = -mrSges(4,3) - mrSges(5,2);
t140 = m(4) * t175 - qJDD(2) * mrSges(4,2) - qJD(2) * t233 + t313 * t226 + t309 * t241 + t296;
t232 = -qJD(2) * mrSges(4,2) - mrSges(4,3) * t241;
t146 = t155 * t280 + t156 * t277;
t235 = -mrSges(5,2) * t241 + qJD(2) * mrSges(5,3);
t292 = -m(5) * t170 + qJDD(2) * mrSges(5,1) + qJD(2) * t235 - t146;
t141 = m(4) * t174 + qJDD(2) * mrSges(4,1) + qJD(2) * t232 + t313 * t227 + t309 * t242 + t292;
t136 = t276 * t140 + t312 * t141;
t230 = -g(3) * t281 - t311;
t244 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t278 + Ifges(3,2) * t281) * qJD(1);
t245 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t278 + Ifges(3,4) * t281) * qJD(1);
t206 = Ifges(4,4) * t242 - Ifges(4,2) * t241 + Ifges(4,6) * qJD(2);
t208 = Ifges(4,1) * t242 - Ifges(4,4) * t241 + Ifges(4,5) * qJD(2);
t203 = Ifges(5,5) * t242 + Ifges(5,6) * qJD(2) + Ifges(5,3) * t241;
t207 = Ifges(5,1) * t242 + Ifges(5,4) * qJD(2) + Ifges(5,5) * t241;
t185 = Ifges(6,4) * t210 + Ifges(6,2) * t209 + Ifges(6,6) * t267;
t186 = Ifges(6,1) * t210 + Ifges(6,4) * t209 + Ifges(6,5) * t267;
t295 = mrSges(6,1) * t160 - mrSges(6,2) * t161 + Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t266 + t210 * t185 - t209 * t186;
t288 = mrSges(5,1) * t170 - mrSges(5,3) * t168 - Ifges(5,4) * t227 - Ifges(5,2) * qJDD(2) - Ifges(5,6) * t226 + pkin(4) * t146 + t242 * t203 - t241 * t207 + t295;
t286 = mrSges(4,2) * t175 - t241 * t208 - qJ(4) * (-mrSges(5,2) * t226 - t215 * t241 + t296) - pkin(3) * (-mrSges(5,2) * t227 - t215 * t242 + t292) - mrSges(4,1) * t174 - t242 * t206 + Ifges(4,6) * t226 - Ifges(4,5) * t227 - Ifges(4,3) * qJDD(2) + t288;
t316 = mrSges(3,1) * t230 - mrSges(3,2) * t231 + Ifges(3,5) * t253 + Ifges(3,6) * t254 + Ifges(3,3) * qJDD(2) + pkin(2) * t136 + (t278 * t244 - t281 * t245) * qJD(1) - t286;
t205 = Ifges(5,4) * t242 + Ifges(5,2) * qJD(2) + Ifges(5,6) * t241;
t310 = -Ifges(4,5) * t242 + Ifges(4,6) * t241 - Ifges(4,3) * qJD(2) - t205;
t252 = (-mrSges(3,1) * t281 + mrSges(3,2) * t278) * qJD(1);
t258 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t307;
t134 = m(3) * t230 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t253 + qJD(2) * t258 - t252 * t308 + t136;
t257 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t308;
t301 = t312 * t140 - t276 * t141;
t135 = m(3) * t231 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t254 - qJD(2) * t257 + t252 * t307 + t301;
t302 = -t134 * t278 + t281 * t135;
t165 = -pkin(7) * t240 + (-pkin(3) - pkin(4)) * t226 + (-pkin(3) * qJD(2) + t236 + t315) * t242 - t318;
t157 = -m(6) * t165 + t182 * mrSges(6,1) - t183 * mrSges(6,2) + t209 * t200 - t210 * t201;
t172 = -0.2e1 * qJD(4) * t242 + (qJD(2) * t242 + t226) * pkin(3) + t318;
t153 = m(5) * t172 + t226 * mrSges(5,1) - t227 * mrSges(5,3) - t242 * t234 + t241 * t235 + t157;
t184 = Ifges(6,5) * t210 + Ifges(6,6) * t209 + Ifges(6,3) * t267;
t149 = -mrSges(6,1) * t165 + mrSges(6,3) * t161 + Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t266 - t184 * t210 + t186 * t267;
t150 = mrSges(6,2) * t165 - mrSges(6,3) * t160 + Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t266 + t184 * t209 - t185 * t267;
t290 = -mrSges(5,1) * t172 + mrSges(5,2) * t168 - pkin(4) * t157 - pkin(7) * t147 - t280 * t149 - t277 * t150;
t131 = -mrSges(4,1) * t199 + mrSges(4,3) * t175 - pkin(3) * t153 + t310 * t242 + (Ifges(4,4) - Ifges(5,5)) * t227 + (-Ifges(4,2) - Ifges(5,3)) * t226 + (Ifges(4,6) - Ifges(5,6)) * qJDD(2) + (t208 + t207) * qJD(2) + t290;
t291 = mrSges(5,2) * t170 - mrSges(5,3) * t172 + Ifges(5,1) * t227 + Ifges(5,4) * qJDD(2) + Ifges(5,5) * t226 - pkin(7) * t146 + qJD(2) * t203 - t149 * t277 + t280 * t150;
t132 = mrSges(4,2) * t199 - mrSges(4,3) * t174 + Ifges(4,1) * t227 - Ifges(4,4) * t226 + Ifges(4,5) * qJDD(2) - qJ(4) * t153 - qJD(2) * t206 + t310 * t241 + t291;
t243 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t278 + Ifges(3,6) * t281) * qJD(1);
t247 = -pkin(6) * t284 - t299;
t289 = m(4) * t199 + t226 * mrSges(4,1) + t227 * mrSges(4,2) + t241 * t232 + t242 * t233 + t153;
t125 = -mrSges(3,1) * t247 + mrSges(3,3) * t231 + Ifges(3,4) * t253 + Ifges(3,2) * t254 + Ifges(3,6) * qJDD(2) - pkin(2) * t289 + qJ(3) * t301 + qJD(2) * t245 + t312 * t131 + t276 * t132 - t243 * t308;
t127 = mrSges(3,2) * t247 - mrSges(3,3) * t230 + Ifges(3,1) * t253 + Ifges(3,4) * t254 + Ifges(3,5) * qJDD(2) - qJ(3) * t136 - qJD(2) * t244 - t276 * t131 + t312 * t132 + t243 * t307;
t287 = -m(3) * t247 + t254 * mrSges(3,1) - t253 * mrSges(3,2) - t257 * t308 + t258 * t307 - t289;
t294 = mrSges(2,1) * t259 - mrSges(2,2) * t260 + Ifges(2,3) * qJDD(1) + pkin(1) * t287 + pkin(6) * t302 + t281 * t125 + t278 * t127;
t151 = m(2) * t259 + qJDD(1) * mrSges(2,1) - t284 * mrSges(2,2) + t287;
t130 = t134 * t281 + t135 * t278;
t128 = m(2) * t260 - mrSges(2,1) * t284 - qJDD(1) * mrSges(2,2) + t302;
t123 = mrSges(2,1) * g(3) + mrSges(2,3) * t260 + t284 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t130 - t316;
t122 = -mrSges(2,2) * g(3) - mrSges(2,3) * t259 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t284 - pkin(6) * t130 - t125 * t278 + t127 * t281;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t282 * t122 - t279 * t123 - pkin(5) * (t128 * t279 + t151 * t282), t122, t127, t132, -t205 * t241 + t291, t150; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t279 * t122 + t282 * t123 + pkin(5) * (t128 * t282 - t151 * t279), t123, t125, t131, -t288, t149; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t294, t294, t316, -t286, Ifges(5,5) * t227 + Ifges(5,6) * qJDD(2) + Ifges(5,3) * t226 - qJD(2) * t207 + t242 * t205 - t290, t295;];
m_new = t1;
