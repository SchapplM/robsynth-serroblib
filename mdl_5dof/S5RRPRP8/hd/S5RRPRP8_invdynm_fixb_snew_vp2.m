% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP8_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:03:24
% EndTime: 2019-12-31 20:03:29
% DurationCPUTime: 2.64s
% Computational Cost: add. (24552->302), mult. (51723->363), div. (0->0), fcn. (28156->6), ass. (0->105)
t292 = sin(qJ(1));
t295 = cos(qJ(1));
t266 = -g(1) * t295 - g(2) * t292;
t297 = qJD(1) ^ 2;
t239 = -pkin(1) * t297 + qJDD(1) * pkin(6) + t266;
t291 = sin(qJ(2));
t294 = cos(qJ(2));
t215 = -t294 * g(3) - t291 * t239;
t216 = -g(3) * t291 + t294 * t239;
t230 = Ifges(4,6) * qJD(2) + (Ifges(4,5) * t291 - Ifges(4,3) * t294) * qJD(1);
t233 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t291 + Ifges(3,2) * t294) * qJD(1);
t253 = (-mrSges(4,1) * t294 - mrSges(4,3) * t291) * qJD(1);
t317 = qJD(1) * qJD(2);
t313 = t294 * t317;
t255 = qJDD(1) * t291 + t313;
t314 = t291 * t317;
t256 = qJDD(1) * t294 - t314;
t252 = (-pkin(2) * t294 - qJ(3) * t291) * qJD(1);
t296 = qJD(2) ^ 2;
t318 = qJD(1) * t294;
t324 = 2 * qJD(3);
t183 = -pkin(2) * t296 + qJDD(2) * qJ(3) + qJD(2) * t324 + t252 * t318 + t216;
t319 = qJD(1) * t291;
t264 = -qJD(2) * pkin(3) - pkin(7) * t319;
t322 = t294 ^ 2 * t297;
t172 = -pkin(3) * t322 - pkin(7) * t256 + qJD(2) * t264 + t183;
t198 = -qJDD(2) * pkin(2) - t296 * qJ(3) + t252 * t319 + qJDD(3) - t215;
t173 = (-t255 + t313) * pkin(7) + (-t291 * t294 * t297 - qJDD(2)) * pkin(3) + t198;
t290 = sin(qJ(4));
t293 = cos(qJ(4));
t165 = -t290 * t172 + t173 * t293;
t236 = (-t290 * t291 - t293 * t294) * qJD(1);
t200 = qJD(4) * t236 + t255 * t293 - t256 * t290;
t237 = (-t290 * t294 + t291 * t293) * qJD(1);
t212 = -mrSges(6,1) * t236 + mrSges(6,2) * t237;
t213 = -mrSges(5,1) * t236 + mrSges(5,2) * t237;
t280 = -qJD(2) + qJD(4);
t218 = -mrSges(5,2) * t280 + mrSges(5,3) * t236;
t279 = -qJDD(2) + qJDD(4);
t155 = -0.2e1 * qJD(5) * t237 + (t236 * t280 - t200) * qJ(5) + (t236 * t237 + t279) * pkin(4) + t165;
t217 = -mrSges(6,2) * t280 + mrSges(6,3) * t236;
t316 = m(6) * t155 + t279 * mrSges(6,1) + t280 * t217;
t147 = m(5) * t165 + t279 * mrSges(5,1) + t280 * t218 + (-t212 - t213) * t237 + (-mrSges(5,3) - mrSges(6,3)) * t200 + t316;
t166 = t172 * t293 + t173 * t290;
t199 = -qJD(4) * t237 - t255 * t290 - t256 * t293;
t220 = mrSges(6,1) * t280 - mrSges(6,3) * t237;
t221 = mrSges(5,1) * t280 - mrSges(5,3) * t237;
t219 = pkin(4) * t280 - qJ(5) * t237;
t229 = t236 ^ 2;
t159 = -pkin(4) * t229 + qJ(5) * t199 + 0.2e1 * qJD(5) * t236 - t219 * t280 + t166;
t315 = m(6) * t159 + mrSges(6,3) * t199 + t212 * t236;
t149 = m(5) * t166 + t199 * mrSges(5,3) + t236 * t213 + (-t220 - t221) * t280 + (-mrSges(5,2) - mrSges(6,2)) * t279 + t315;
t142 = t293 * t147 + t290 * t149;
t152 = -t200 * mrSges(6,3) - t237 * t212 + t316;
t205 = Ifges(5,4) * t237 + Ifges(5,2) * t236 + Ifges(5,6) * t280;
t207 = Ifges(5,1) * t237 + Ifges(5,4) * t236 + Ifges(5,5) * t280;
t204 = Ifges(6,4) * t237 + Ifges(6,2) * t236 + Ifges(6,6) * t280;
t206 = Ifges(6,1) * t237 + Ifges(6,4) * t236 + Ifges(6,5) * t280;
t310 = mrSges(6,1) * t155 - mrSges(6,2) * t159 + Ifges(6,5) * t200 + Ifges(6,6) * t199 + Ifges(6,3) * t279 + t204 * t237 - t206 * t236;
t301 = mrSges(5,1) * t165 - mrSges(5,2) * t166 + Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * t279 + pkin(4) * t152 + t205 * t237 - t207 * t236 + t310;
t300 = -mrSges(4,1) * t198 + mrSges(4,3) * t183 + Ifges(4,4) * t255 + Ifges(4,2) * qJDD(2) - Ifges(4,6) * t256 - pkin(3) * t142 - t301;
t263 = mrSges(4,2) * t318 + qJD(2) * mrSges(4,3);
t304 = -m(4) * t198 + qJDD(2) * mrSges(4,1) + qJD(2) * t263 - t142;
t143 = -t290 * t147 + t149 * t293;
t261 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t319;
t306 = m(4) * t183 + qJDD(2) * mrSges(4,3) + qJD(2) * t261 + t253 * t318 + t143;
t234 = Ifges(4,4) * qJD(2) + (Ifges(4,1) * t291 - Ifges(4,5) * t294) * qJD(1);
t320 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t291 + Ifges(3,4) * t294) * qJD(1) + t234;
t326 = -((t230 - t233) * t291 + t294 * t320) * qJD(1) + mrSges(3,1) * t215 - mrSges(3,2) * t216 + Ifges(3,5) * t255 + Ifges(3,6) * t256 + Ifges(3,3) * qJDD(2) + pkin(2) * (-t255 * mrSges(4,2) - t253 * t319 + t304) + qJ(3) * (mrSges(4,2) * t256 + t306) + t300;
t323 = mrSges(3,3) + mrSges(4,2);
t265 = t292 * g(1) - t295 * g(2);
t254 = (-mrSges(3,1) * t294 + mrSges(3,2) * t291) * qJD(1);
t260 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t319;
t136 = m(3) * t216 - qJDD(2) * mrSges(3,2) - qJD(2) * t260 + t254 * t318 + t256 * t323 + t306;
t262 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t318;
t137 = m(3) * t215 + qJDD(2) * mrSges(3,1) + qJD(2) * t262 - t323 * t255 + (-t253 - t254) * t319 + t304;
t312 = t136 * t294 - t137 * t291;
t238 = -qJDD(1) * pkin(1) - t297 * pkin(6) - t265;
t309 = -t256 * pkin(2) + t238 + (-t255 - t313) * qJ(3);
t168 = -pkin(2) * t314 + t256 * pkin(3) - pkin(7) * t322 - t309 + (t264 + t324) * t319;
t162 = -pkin(4) * t199 - qJ(5) * t229 + t219 * t237 + qJDD(5) + t168;
t311 = m(6) * t162 - mrSges(6,1) * t199 + mrSges(6,2) * t200 - t236 * t217 + t237 * t220;
t308 = -mrSges(6,1) * t162 + mrSges(6,3) * t159 + Ifges(6,4) * t200 + Ifges(6,2) * t199 + Ifges(6,6) * t279 + t280 * t206;
t202 = Ifges(6,5) * t237 + Ifges(6,6) * t236 + Ifges(6,3) * t280;
t307 = mrSges(6,2) * t162 - mrSges(6,3) * t155 + Ifges(6,1) * t200 + Ifges(6,4) * t199 + Ifges(6,5) * t279 + t202 * t236;
t150 = -m(5) * t168 + mrSges(5,1) * t199 - mrSges(5,2) * t200 + t236 * t218 - t237 * t221 - t311;
t174 = (pkin(2) * qJD(2) - (2 * qJD(3))) * t319 + t309;
t146 = m(4) * t174 - mrSges(4,1) * t256 - t255 * mrSges(4,3) - t261 * t319 - t263 * t318 + t150;
t231 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t291 + Ifges(3,6) * t294) * qJD(1);
t232 = Ifges(4,2) * qJD(2) + (Ifges(4,4) * t291 - Ifges(4,6) * t294) * qJD(1);
t203 = Ifges(5,5) * t237 + Ifges(5,6) * t236 + Ifges(5,3) * t280;
t134 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * t279 + t280 * t207 - mrSges(5,1) * t168 + mrSges(5,3) * t166 - pkin(4) * t311 + qJ(5) * (-t279 * mrSges(6,2) - t280 * t220 + t315) + (-t203 - t202) * t237 + t308;
t139 = mrSges(5,2) * t168 - mrSges(5,3) * t165 + Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * t279 - qJ(5) * t152 + t236 * t203 + (-t204 - t205) * t280 + t307;
t302 = -mrSges(4,1) * t174 + mrSges(4,2) * t183 - pkin(3) * t150 - pkin(7) * t143 - t293 * t134 - t290 * t139;
t128 = -mrSges(3,1) * t238 + mrSges(3,3) * t216 - pkin(2) * t146 + (Ifges(3,2) + Ifges(4,3)) * t256 + (Ifges(3,4) - Ifges(4,5)) * t255 + (Ifges(3,6) - Ifges(4,6)) * qJDD(2) + t320 * qJD(2) + (-t231 - t232) * t319 + t302;
t303 = mrSges(4,2) * t198 - mrSges(4,3) * t174 + Ifges(4,1) * t255 + Ifges(4,4) * qJDD(2) - Ifges(4,5) * t256 - pkin(7) * t142 + qJD(2) * t230 - t134 * t290 + t139 * t293 + t232 * t318;
t130 = mrSges(3,2) * t238 - mrSges(3,3) * t215 + Ifges(3,1) * t255 + Ifges(3,4) * t256 + Ifges(3,5) * qJDD(2) - qJ(3) * t146 - qJD(2) * t233 + t231 * t318 + t303;
t299 = -m(3) * t238 + t256 * mrSges(3,1) - mrSges(3,2) * t255 - t260 * t319 + t262 * t318 - t146;
t305 = mrSges(2,1) * t265 - mrSges(2,2) * t266 + Ifges(2,3) * qJDD(1) + pkin(1) * t299 + pkin(6) * t312 + t128 * t294 + t130 * t291;
t144 = m(2) * t265 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t297 + t299;
t133 = t136 * t291 + t137 * t294;
t131 = m(2) * t266 - mrSges(2,1) * t297 - qJDD(1) * mrSges(2,2) + t312;
t126 = mrSges(2,1) * g(3) + mrSges(2,3) * t266 + t297 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t133 - t326;
t125 = -mrSges(2,2) * g(3) - mrSges(2,3) * t265 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t297 - pkin(6) * t133 - t128 * t291 + t130 * t294;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t295 * t125 - t292 * t126 - pkin(5) * (t131 * t292 + t144 * t295), t125, t130, t303, t139, -t204 * t280 + t307; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t292 * t125 + t295 * t126 + pkin(5) * (t131 * t295 - t144 * t292), t126, t128, t300 + (-t291 * t230 - t294 * t234) * qJD(1), t134, -t237 * t202 + t308; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t305, t305, t326, Ifges(4,5) * t255 + Ifges(4,6) * qJDD(2) - Ifges(4,3) * t256 - qJD(2) * t234 + t232 * t319 - t302, t301, t310;];
m_new = t1;
