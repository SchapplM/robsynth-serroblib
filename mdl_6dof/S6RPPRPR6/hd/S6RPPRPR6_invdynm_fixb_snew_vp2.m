% Calculate vector of cutting torques with Newton-Euler for
% S6RPPRPR6
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
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
% m [3x7]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:27
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S6RPPRPR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR6_invdynm_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:26:31
% EndTime: 2019-05-05 14:26:36
% DurationCPUTime: 2.41s
% Computational Cost: add. (24275->305), mult. (46307->342), div. (0->0), fcn. (20299->6), ass. (0->115)
t332 = 2 * qJD(1);
t280 = sin(qJ(1));
t283 = cos(qJ(1));
t248 = -g(1) * t283 - g(2) * t280;
t331 = qJDD(1) * qJ(2) + (qJD(2) * t332) + t248;
t247 = t280 * g(1) - t283 * g(2);
t286 = qJD(1) ^ 2;
t210 = -qJDD(1) * pkin(1) - t286 * qJ(2) + qJDD(2) - t247;
t307 = qJDD(1) * qJ(3) + (qJD(3) * t332) - t210;
t330 = -2 * qJD(5);
t329 = pkin(4) + pkin(8);
t328 = pkin(8) * t286;
t279 = sin(qJ(4));
t327 = g(3) * t279;
t326 = t286 * pkin(7);
t325 = mrSges(3,2) - mrSges(4,3);
t324 = (Ifges(3,4) - Ifges(4,5));
t323 = Ifges(5,4) + Ifges(6,6);
t322 = Ifges(2,6) - Ifges(3,5);
t321 = t286 * mrSges(4,3);
t201 = qJDD(3) + (-pkin(1) - qJ(3)) * t286 + t331;
t196 = -qJDD(1) * pkin(7) + t201;
t282 = cos(qJ(4));
t320 = t196 * t282;
t190 = t320 + t327;
t317 = qJD(1) * qJD(4);
t253 = t279 * t317;
t239 = qJDD(1) * t282 - t253;
t319 = qJD(1) * t279;
t242 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t319;
t244 = mrSges(6,1) * t319 - (qJD(4) * mrSges(6,3));
t314 = t282 * t317;
t238 = qJDD(1) * t279 + t314;
t318 = qJD(1) * t282;
t246 = pkin(5) * t318 - (qJD(4) * pkin(8));
t272 = t279 ^ 2;
t295 = pkin(4) * t314 + t318 * t330 + t307 + (-t239 + t253) * qJ(5);
t168 = t295 + t329 * t238 + (-pkin(5) * t272 - pkin(7)) * t286 - t246 * t318;
t235 = (pkin(4) * t279 - qJ(5) * t282) * qJD(1);
t285 = qJD(4) ^ 2;
t306 = -t285 * qJ(5) + t235 * t318 + qJDD(5) - t320;
t171 = pkin(5) * t239 - t329 * qJDD(4) + (pkin(5) * t317 + t282 * t328 - g(3)) * t279 + t306;
t278 = sin(qJ(6));
t281 = cos(qJ(6));
t166 = -t168 * t278 + t171 * t281;
t233 = -qJD(4) * t278 + t281 * t319;
t189 = qJD(6) * t233 + qJDD(4) * t281 + t238 * t278;
t234 = qJD(4) * t281 + t278 * t319;
t193 = -mrSges(7,1) * t233 + mrSges(7,2) * t234;
t251 = qJD(6) + t318;
t202 = -mrSges(7,2) * t251 + mrSges(7,3) * t233;
t232 = qJDD(6) + t239;
t162 = m(7) * t166 + mrSges(7,1) * t232 - t189 * mrSges(7,3) - t193 * t234 + t202 * t251;
t167 = t168 * t281 + t171 * t278;
t188 = -qJD(6) * t234 - qJDD(4) * t278 + t238 * t281;
t203 = mrSges(7,1) * t251 - mrSges(7,3) * t234;
t163 = m(7) * t167 - mrSges(7,2) * t232 + t188 * mrSges(7,3) + t193 * t233 - t203 * t251;
t149 = t162 * t281 + t163 * t278;
t177 = -qJDD(4) * pkin(4) + t306 - t327;
t303 = -m(6) * t177 - t239 * mrSges(6,1) - t149;
t236 = (-mrSges(6,2) * t279 - mrSges(6,3) * t282) * qJD(1);
t311 = qJD(1) * (-t236 - (mrSges(5,1) * t279 + mrSges(5,2) * t282) * qJD(1));
t144 = m(5) * t190 - mrSges(5,3) * t239 + (mrSges(5,1) - mrSges(6,2)) * qJDD(4) + (t242 - t244) * qJD(4) + t282 * t311 + t303;
t191 = -g(3) * t282 + t279 * t196;
t243 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t318;
t298 = -t285 * pkin(4) + qJDD(4) * qJ(5) - t235 * t319 + t191;
t170 = -t272 * t328 - pkin(5) * t238 + ((2 * qJD(5)) + t246) * qJD(4) + t298;
t164 = -m(7) * t170 + t188 * mrSges(7,1) - t189 * mrSges(7,2) + t202 * t233 - t234 * t203;
t175 = (qJD(4) * t330) - t298;
t245 = mrSges(6,1) * t318 + qJD(4) * mrSges(6,2);
t297 = -m(6) * t175 + qJDD(4) * mrSges(6,3) + qJD(4) * t245 - t164;
t157 = m(5) * t191 - qJDD(4) * mrSges(5,2) - qJD(4) * t243 + (-mrSges(5,3) - mrSges(6,1)) * t238 + t279 * t311 + t297;
t136 = t282 * t144 + t279 * t157;
t150 = -t278 * t162 + t281 * t163;
t313 = -t144 * t279 + t282 * t157;
t217 = (Ifges(6,1) * qJD(4)) + (-Ifges(6,4) * t282 + Ifges(6,5) * t279) * qJD(1);
t312 = qJD(1) * (-(Ifges(5,3) * qJD(4)) - (Ifges(5,5) * t282 - Ifges(5,6) * t279) * qJD(1) - t217);
t310 = m(4) * t201 + qJDD(1) * mrSges(4,2) + t136;
t173 = pkin(4) * t238 + t295 - t326;
t145 = m(6) * t173 - t238 * mrSges(6,2) - t239 * mrSges(6,3) - t244 * t319 - t245 * t318 + t150;
t195 = t307 - t326;
t214 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t282 - Ifges(5,4) * t279) * qJD(1);
t216 = Ifges(6,4) * qJD(4) + (-Ifges(6,2) * t282 + Ifges(6,6) * t279) * qJD(1);
t179 = Ifges(7,5) * t234 + Ifges(7,6) * t233 + Ifges(7,3) * t251;
t181 = Ifges(7,1) * t234 + Ifges(7,4) * t233 + Ifges(7,5) * t251;
t153 = -mrSges(7,1) * t170 + mrSges(7,3) * t167 + Ifges(7,4) * t189 + Ifges(7,2) * t188 + Ifges(7,6) * t232 - t179 * t234 + t181 * t251;
t180 = Ifges(7,4) * t234 + Ifges(7,2) * t233 + Ifges(7,6) * t251;
t154 = mrSges(7,2) * t170 - mrSges(7,3) * t166 + Ifges(7,1) * t189 + Ifges(7,4) * t188 + Ifges(7,5) * t232 + t179 * t233 - t180 * t251;
t294 = -mrSges(6,1) * t175 + mrSges(6,2) * t173 - pkin(5) * t164 - pkin(8) * t150 - t281 * t153 - t278 * t154;
t126 = -mrSges(5,1) * t195 + mrSges(5,3) * t191 - pkin(4) * t145 + t323 * t239 + (-Ifges(5,2) - Ifges(6,3)) * t238 + (Ifges(5,6) - Ifges(6,5)) * qJDD(4) + (t214 - t216) * qJD(4) + t282 * t312 + t294;
t213 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t282 - Ifges(5,2) * t279) * qJD(1);
t215 = Ifges(6,5) * qJD(4) + (-Ifges(6,6) * t282 + Ifges(6,3) * t279) * qJD(1);
t301 = mrSges(7,1) * t166 - mrSges(7,2) * t167 + Ifges(7,5) * t189 + Ifges(7,6) * t188 + Ifges(7,3) * t232 + t234 * t180 - t233 * t181;
t293 = mrSges(6,1) * t177 - mrSges(6,3) * t173 + pkin(5) * t149 + t301;
t129 = t293 + (Ifges(5,1) + Ifges(6,2)) * t239 - t323 * t238 + (Ifges(5,5) - Ifges(6,4)) * qJDD(4) + (-t213 + t215) * qJD(4) - mrSges(5,3) * t190 + mrSges(5,2) * t195 - qJ(5) * t145 + t279 * t312;
t299 = -m(5) * t195 - t238 * mrSges(5,1) - t239 * mrSges(5,2) - t242 * t319 - t243 * t318 - t145;
t305 = -mrSges(4,1) * t307 + mrSges(4,2) * g(3) + (t286 * Ifges(4,4)) + Ifges(4,5) * qJDD(1) + pkin(3) * t299 + pkin(7) * t313 + t282 * t126 + t279 * t129;
t208 = t286 * pkin(1) - t331;
t304 = -m(3) * t208 + (t286 * mrSges(3,2)) + qJDD(1) * mrSges(3,3) + t310;
t302 = mrSges(6,2) * t177 - mrSges(6,3) * t175 + Ifges(6,1) * qJDD(4) - Ifges(6,4) * t239 + Ifges(6,5) * t238 - pkin(8) * t149 - t278 * t153 + t281 * t154 - t215 * t318 - t216 * t319;
t300 = mrSges(4,2) * t201 + mrSges(4,3) * t307 + Ifges(4,1) * qJDD(1) - pkin(7) * t136 - t126 * t279 + t282 * t129;
t140 = -m(4) * t307 - t286 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t299;
t296 = mrSges(3,1) * t210 + pkin(2) * t140 + t305;
t292 = mrSges(3,2) * t210 - mrSges(3,3) * t208 + Ifges(3,1) * qJDD(1) - qJ(3) * t140 + t300;
t291 = t302 + t214 * t319 + t213 * t318 - mrSges(5,2) * t191 + mrSges(5,1) * t190 + pkin(4) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t244 - t236 * t318 + t303) + Ifges(5,3) * qJDD(4) + qJ(5) * (-mrSges(6,1) * t238 - t236 * t319 + t297) + Ifges(5,5) * t239 - Ifges(5,6) * t238;
t290 = -m(3) * t210 + t286 * mrSges(3,3) - t140;
t289 = -mrSges(2,2) * t248 + qJ(2) * (t304 - t321) + pkin(1) * (-qJDD(1) * mrSges(3,2) + t290) + mrSges(2,1) * t247 + Ifges(2,3) * qJDD(1) + t292;
t288 = mrSges(4,1) * t201 - Ifges(4,4) * qJDD(1) + pkin(3) * t136 + t291;
t287 = mrSges(3,1) * t208 + pkin(2) * (-t310 + t321) + qJ(3) * (-m(4) * g(3) + t313) - t288;
t138 = t290 + (mrSges(2,1) - mrSges(3,2)) * qJDD(1) - t286 * mrSges(2,2) + m(2) * t247;
t133 = (-m(3) - m(4)) * g(3) + t313;
t130 = m(2) * t248 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - mrSges(4,3)) * t286 + t304;
t124 = -t287 + ((Ifges(2,5) - t324) * t286) + t322 * qJDD(1) + (mrSges(2,1) - t325) * g(3) + mrSges(2,3) * t248 - pkin(1) * t133;
t123 = t296 - t322 * t286 + (-Ifges(3,4) + Ifges(2,5)) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t247 - qJ(2) * t133;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t283 * t123 - t280 * t124 - pkin(6) * (t130 * t280 + t138 * t283), t123, t292, t300, t129, t302, t154; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t280 * t123 + t283 * t124 + pkin(6) * (t130 * t283 - t138 * t280), t124, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t286 * Ifges(3,5)) - t296, -mrSges(4,3) * g(3) - t286 * Ifges(4,5) - t288, t126, Ifges(6,4) * qJDD(4) - Ifges(6,2) * t239 + Ifges(6,6) * t238 - qJD(4) * t215 + t217 * t319 - t293, t153; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t289, t289, Ifges(3,5) * qJDD(1) + g(3) * t325 + (t286 * t324) + t287, t305, t291, Ifges(6,5) * qJDD(4) - Ifges(6,6) * t239 + Ifges(6,3) * t238 + qJD(4) * t216 + t217 * t318 - t294, t301;];
m_new  = t1;
