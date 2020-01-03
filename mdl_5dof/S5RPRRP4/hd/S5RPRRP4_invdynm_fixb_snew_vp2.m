% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:22
% EndTime: 2020-01-03 11:49:33
% DurationCPUTime: 5.06s
% Computational Cost: add. (49974->283), mult. (120788->358), div. (0->0), fcn. (79109->8), ass. (0->116)
t270 = sin(qJ(1));
t273 = cos(qJ(1));
t249 = -t270 * g(2) + t273 * g(3);
t274 = qJD(1) ^ 2;
t315 = -t274 * pkin(1) + qJDD(1) * qJ(2) + (2 * qJD(1) * qJD(2)) + t249;
t268 = sin(qJ(4));
t269 = sin(qJ(3));
t271 = cos(qJ(4));
t272 = cos(qJ(3));
t266 = sin(pkin(8));
t307 = qJD(1) * t266;
t223 = (-t268 * t272 - t269 * t271) * t307;
t303 = qJD(1) * qJD(3);
t232 = (-qJDD(1) * t269 - t272 * t303) * t266;
t233 = (qJDD(1) * t272 - t269 * t303) * t266;
t192 = t223 * qJD(4) + t268 * t232 + t271 * t233;
t224 = (-t268 * t269 + t271 * t272) * t307;
t205 = -t223 * mrSges(6,1) + t224 * mrSges(6,2);
t267 = cos(pkin(8));
t216 = -t266 * g(1) + t315 * t267;
t292 = -pkin(2) * t267 - pkin(6) * t266;
t241 = t292 * qJD(1);
t305 = t267 * qJD(1);
t204 = t241 * t305 + t216;
t250 = -t273 * g(2) - t270 * g(3);
t285 = -t274 * qJ(2) + qJDD(2) - t250;
t217 = (-pkin(1) + t292) * qJDD(1) + t285;
t214 = t272 * t217;
t302 = t267 * qJDD(1);
t252 = qJDD(3) - t302;
t253 = qJD(3) - t305;
t310 = t266 ^ 2 * t274;
t169 = t252 * pkin(3) - t233 * pkin(7) + t214 + (-pkin(3) * t272 * t310 - pkin(7) * t253 * t307 - t204) * t269;
t173 = t272 * t204 + t269 * t217;
t297 = t272 * t307;
t231 = t253 * pkin(3) - pkin(7) * t297;
t301 = t269 ^ 2 * t310;
t170 = -pkin(3) * t301 + t232 * pkin(7) - t253 * t231 + t173;
t161 = t271 * t169 - t268 * t170;
t248 = qJDD(4) + t252;
t251 = qJD(4) + t253;
t156 = -0.2e1 * qJD(5) * t224 + (t223 * t251 - t192) * qJ(5) + (t223 * t224 + t248) * pkin(4) + t161;
t208 = -t251 * mrSges(6,2) + t223 * mrSges(6,3);
t300 = m(6) * t156 + t248 * mrSges(6,1) + t251 * t208;
t153 = -t192 * mrSges(6,3) - t224 * t205 + t300;
t162 = t268 * t169 + t271 * t170;
t191 = -t224 * qJD(4) + t271 * t232 - t268 * t233;
t196 = Ifges(5,4) * t224 + Ifges(5,2) * t223 + Ifges(5,6) * t251;
t197 = Ifges(6,1) * t224 + Ifges(6,4) * t223 + Ifges(6,5) * t251;
t198 = Ifges(5,1) * t224 + Ifges(5,4) * t223 + Ifges(5,5) * t251;
t210 = t251 * pkin(4) - t224 * qJ(5);
t222 = t223 ^ 2;
t159 = -t222 * pkin(4) + t191 * qJ(5) + 0.2e1 * qJD(5) * t223 - t251 * t210 + t162;
t195 = Ifges(6,4) * t224 + Ifges(6,2) * t223 + Ifges(6,6) * t251;
t283 = -mrSges(6,1) * t156 + mrSges(6,2) * t159 - Ifges(6,5) * t192 - Ifges(6,6) * t191 - Ifges(6,3) * t248 - t224 * t195;
t314 = mrSges(5,1) * t161 - mrSges(5,2) * t162 + Ifges(5,5) * t192 + Ifges(5,6) * t191 + Ifges(5,3) * t248 + pkin(4) * t153 + t224 * t196 - t283 + (-t198 - t197) * t223;
t206 = -t223 * mrSges(5,1) + t224 * mrSges(5,2);
t209 = -t251 * mrSges(5,2) + t223 * mrSges(5,3);
t146 = m(5) * t161 + t248 * mrSges(5,1) + t251 * t209 + (-t205 - t206) * t224 + (-mrSges(5,3) - mrSges(6,3)) * t192 + t300;
t211 = t251 * mrSges(6,1) - t224 * mrSges(6,3);
t212 = t251 * mrSges(5,1) - t224 * mrSges(5,3);
t299 = m(6) * t159 + t191 * mrSges(6,3) + t223 * t205;
t149 = m(5) * t162 + t191 * mrSges(5,3) + t223 * t206 + (-t211 - t212) * t251 + (-mrSges(5,2) - mrSges(6,2)) * t248 + t299;
t144 = t271 * t146 + t268 * t149;
t172 = -t269 * t204 + t214;
t313 = -mrSges(4,1) * t172 + mrSges(4,2) * t173 - Ifges(4,5) * t233 - Ifges(4,6) * t232 - Ifges(4,3) * t252 - pkin(3) * t144 - t314;
t215 = -t267 * g(1) - t315 * t266;
t203 = t241 * t307 - t215;
t171 = -t232 * pkin(3) - pkin(7) * t301 + t231 * t297 + t203;
t193 = Ifges(6,5) * t224 + Ifges(6,6) * t223 + Ifges(6,3) * t251;
t194 = Ifges(5,5) * t224 + Ifges(5,6) * t223 + Ifges(5,3) * t251;
t165 = -t191 * pkin(4) - t222 * qJ(5) + t224 * t210 + qJDD(5) + t171;
t284 = -mrSges(6,1) * t165 + mrSges(6,3) * t159 + Ifges(6,4) * t192 + Ifges(6,2) * t191 + Ifges(6,6) * t248 + t251 * t197;
t289 = m(6) * t165 - t191 * mrSges(6,1) + t192 * mrSges(6,2) - t223 * t208 + t224 * t211;
t139 = Ifges(5,4) * t192 + Ifges(5,2) * t191 + Ifges(5,6) * t248 + t251 * t198 - mrSges(5,1) * t171 + mrSges(5,3) * t162 - pkin(4) * t289 + qJ(5) * (-t248 * mrSges(6,2) - t251 * t211 + t299) + (-t194 - t193) * t224 + t284;
t282 = mrSges(6,2) * t165 - mrSges(6,3) * t156 + Ifges(6,1) * t192 + Ifges(6,4) * t191 + Ifges(6,5) * t248 + t223 * t193;
t143 = mrSges(5,2) * t171 - mrSges(5,3) * t161 + Ifges(5,1) * t192 + Ifges(5,4) * t191 + Ifges(5,5) * t248 - qJ(5) * t153 + t223 * t194 + (-t195 - t196) * t251 + t282;
t218 = Ifges(4,3) * t253 + (Ifges(4,5) * t272 - Ifges(4,6) * t269) * t307;
t220 = Ifges(4,5) * t253 + (Ifges(4,1) * t272 - Ifges(4,4) * t269) * t307;
t279 = m(5) * t171 - t191 * mrSges(5,1) + t192 * mrSges(5,2) - t223 * t209 + t224 * t212 + t289;
t293 = -t268 * t146 + t271 * t149;
t128 = -mrSges(4,1) * t203 + mrSges(4,3) * t173 + Ifges(4,4) * t233 + Ifges(4,2) * t232 + Ifges(4,6) * t252 - pkin(3) * t279 + pkin(7) * t293 + t271 * t139 + t268 * t143 - t218 * t297 + t253 * t220;
t219 = Ifges(4,6) * t253 + (Ifges(4,4) * t272 - Ifges(4,2) * t269) * t307;
t298 = t269 * t307;
t129 = mrSges(4,2) * t203 - mrSges(4,3) * t172 + Ifges(4,1) * t233 + Ifges(4,4) * t232 + Ifges(4,5) * t252 - pkin(7) * t144 - t268 * t139 + t271 * t143 - t218 * t298 - t253 * t219;
t228 = -t253 * mrSges(4,2) - mrSges(4,3) * t298;
t230 = (mrSges(4,1) * t269 + mrSges(4,2) * t272) * t307;
t141 = m(4) * t172 + t252 * mrSges(4,1) - t233 * mrSges(4,3) + t253 * t228 - t230 * t297 + t144;
t229 = t253 * mrSges(4,1) - mrSges(4,3) * t297;
t142 = m(4) * t173 - t252 * mrSges(4,2) + t232 * mrSges(4,3) - t253 * t229 - t230 * t298 + t293;
t138 = -t269 * t141 + t272 * t142;
t276 = -m(4) * t203 + t232 * mrSges(4,1) - t233 * mrSges(4,2) - t279;
t287 = -t228 * t269 - t229 * t272;
t291 = Ifges(3,1) * t266 + Ifges(3,4) * t267;
t312 = -((Ifges(3,4) * t266 + Ifges(3,2) * t267) * t307 - t291 * t305) * qJD(1) - mrSges(3,1) * t215 + mrSges(3,2) * t216 - pkin(2) * (t287 * t307 + t276) - pkin(6) * t138 - t272 * t128 - t269 * t129;
t311 = mrSges(3,2) * t266;
t306 = qJDD(1) * mrSges(3,3);
t237 = (-mrSges(3,1) * t267 + t311) * qJD(1);
t135 = m(3) * t216 + (qJD(1) * t237 + t306) * t267 + t138;
t150 = m(3) * t215 + t276 + (-t306 + (-t237 + t287) * qJD(1)) * t266;
t294 = t267 * t135 - t266 * t150;
t290 = Ifges(3,5) * t266 + Ifges(3,6) * t267;
t137 = t272 * t141 + t269 * t142;
t288 = t219 * t272 + t220 * t269;
t235 = -qJDD(1) * pkin(1) + t285;
t238 = t290 * qJD(1);
t125 = mrSges(3,2) * t235 - mrSges(3,3) * t215 - pkin(6) * t137 + t291 * qJDD(1) - t269 * t128 + t272 * t129 + t238 * t305;
t127 = Ifges(3,2) * t302 + (Ifges(3,4) * qJDD(1) + (-t238 - t288) * qJD(1)) * t266 - mrSges(3,1) * t235 + mrSges(3,3) * t216 - pkin(2) * t137 + t313;
t280 = -m(3) * t235 + mrSges(3,1) * t302 - t137 + (t267 ^ 2 * t274 + t310) * mrSges(3,3);
t281 = -mrSges(2,2) * t249 + qJ(2) * t294 + t266 * t125 + t267 * t127 + pkin(1) * (-qJDD(1) * t311 + t280) + mrSges(2,1) * t250 + Ifges(2,3) * qJDD(1);
t133 = m(2) * t250 - t274 * mrSges(2,2) + (mrSges(2,1) - t311) * qJDD(1) + t280;
t132 = t266 * t135 + t267 * t150;
t130 = m(2) * t249 - t274 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t294;
t123 = mrSges(2,1) * g(1) + mrSges(2,3) * t249 + t274 * Ifges(2,5) - pkin(1) * t132 + (Ifges(2,6) - t290) * qJDD(1) + t312;
t122 = -mrSges(2,2) * g(1) - mrSges(2,3) * t250 + Ifges(2,5) * qJDD(1) - t274 * Ifges(2,6) - qJ(2) * t132 + t267 * t125 - t266 * t127;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t281, t122, t125, t129, t143, -t251 * t195 + t282; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t270 * t122 + t273 * t123 - pkin(5) * (-t273 * t130 + t270 * t133), t123, t127, t128, t139, -t224 * t193 + t284; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t273 * t122 + t270 * t123 + pkin(5) * (t270 * t130 + t273 * t133), t281, t290 * qJDD(1) - t312, t288 * t307 - t313, t314, -t223 * t197 - t283;];
m_new = t1;
