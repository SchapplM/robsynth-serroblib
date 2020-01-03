% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPP5
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPP5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:57:39
% EndTime: 2019-12-31 20:57:45
% DurationCPUTime: 3.00s
% Computational Cost: add. (29345->310), mult. (61636->358), div. (0->0), fcn. (37556->6), ass. (0->110)
t275 = sin(qJ(3));
t278 = cos(qJ(2));
t302 = qJD(1) * t278;
t276 = sin(qJ(2));
t303 = qJD(1) * t276;
t312 = cos(qJ(3));
t240 = t275 * t303 - t312 * t302;
t301 = qJD(1) * qJD(2);
t249 = t276 * qJDD(1) + t278 * t301;
t250 = t278 * qJDD(1) - t276 * t301;
t202 = -t240 * qJD(3) + t312 * t249 + t275 * t250;
t254 = qJD(2) * pkin(2) - pkin(7) * t303;
t274 = t278 ^ 2;
t280 = qJD(1) ^ 2;
t277 = sin(qJ(1));
t279 = cos(qJ(1));
t255 = t277 * g(1) - t279 * g(2);
t294 = qJDD(1) * pkin(1) + t255;
t203 = -t250 * pkin(2) + t254 * t303 - (pkin(7) * t274 + pkin(6)) * t280 - t294;
t272 = qJD(2) + qJD(3);
t308 = t240 * t272;
t319 = t203 + (-t202 + t308) * qJ(4);
t256 = -t279 * g(1) - t277 * g(2);
t243 = -t280 * pkin(1) + qJDD(1) * pkin(6) + t256;
t306 = t276 * t243;
t310 = pkin(2) * t280;
t176 = qJDD(2) * pkin(2) - t249 * pkin(7) - t306 + (pkin(7) * t301 + t276 * t310 - g(3)) * t278;
t227 = -t276 * g(3) + t278 * t243;
t177 = t250 * pkin(7) - qJD(2) * t254 - t274 * t310 + t227;
t167 = t275 * t176 + t312 * t177;
t241 = (t275 * t278 + t312 * t276) * qJD(1);
t220 = t240 * pkin(3) - t241 * qJ(4);
t270 = t272 ^ 2;
t271 = qJDD(2) + qJDD(3);
t313 = 2 * qJD(4);
t162 = -t270 * pkin(3) + t271 * qJ(4) - t240 * t220 + t272 * t313 + t167;
t233 = -t272 * mrSges(5,1) + t241 * mrSges(5,2);
t317 = m(5) * t162 + t271 * mrSges(5,3) + t272 * t233;
t166 = t312 * t176 - t275 * t177;
t164 = -t271 * pkin(3) - t270 * qJ(4) + t241 * t220 + qJDD(4) - t166;
t234 = -t240 * mrSges(5,2) + t272 * mrSges(5,3);
t316 = -m(5) * t164 + t271 * mrSges(5,1) + t272 * t234;
t228 = t272 * mrSges(6,2) + t240 * mrSges(6,3);
t229 = -t272 * mrSges(4,2) - t240 * mrSges(4,3);
t314 = -0.2e1 * t241;
t153 = qJD(5) * t314 + (-t202 - t308) * qJ(5) + (t240 * t241 - t271) * pkin(4) + t164;
t222 = -t240 * mrSges(6,1) + t241 * mrSges(6,2);
t295 = -m(6) * t153 + t202 * mrSges(6,3) + t241 * t222;
t221 = t240 * mrSges(5,1) - t241 * mrSges(5,3);
t304 = -t240 * mrSges(4,1) - t241 * mrSges(4,2) - t221;
t309 = -mrSges(4,3) - mrSges(5,2);
t137 = m(4) * t166 + (t228 + t229) * t272 + (mrSges(4,1) + mrSges(6,1)) * t271 + t304 * t241 + t309 * t202 + t295 + t316;
t201 = t241 * qJD(3) + t275 * t249 - t312 * t250;
t231 = -t272 * mrSges(6,1) - t241 * mrSges(6,3);
t232 = t272 * mrSges(4,1) - t241 * mrSges(4,3);
t230 = -t272 * pkin(4) - t241 * qJ(5);
t236 = t240 ^ 2;
t157 = -t236 * pkin(4) + t201 * qJ(5) + 0.2e1 * qJD(5) * t240 + t272 * t230 + t162;
t300 = m(6) * t157 + t201 * mrSges(6,3) + t240 * t222;
t140 = m(4) * t167 + (t231 - t232) * t272 + (-mrSges(4,2) + mrSges(6,2)) * t271 + t304 * t240 + t309 * t201 + t300 + t317;
t133 = t312 * t137 + t275 * t140;
t226 = -t278 * g(3) - t306;
t238 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t276 + Ifges(3,2) * t278) * qJD(1);
t239 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t276 + Ifges(3,4) * t278) * qJD(1);
t147 = -t271 * mrSges(6,1) - t272 * t228 - t295;
t148 = t271 * mrSges(6,2) + t272 * t231 + t300;
t211 = Ifges(4,4) * t241 - Ifges(4,2) * t240 + Ifges(4,6) * t272;
t214 = Ifges(4,1) * t241 - Ifges(4,4) * t240 + Ifges(4,5) * t272;
t207 = Ifges(5,5) * t241 + Ifges(5,6) * t272 + Ifges(5,3) * t240;
t213 = Ifges(5,1) * t241 + Ifges(5,4) * t272 + Ifges(5,5) * t240;
t209 = Ifges(6,4) * t241 + Ifges(6,2) * t240 - Ifges(6,6) * t272;
t212 = Ifges(6,1) * t241 + Ifges(6,4) * t240 - Ifges(6,5) * t272;
t292 = mrSges(6,1) * t153 - mrSges(6,2) * t157 + Ifges(6,5) * t202 + Ifges(6,6) * t201 - Ifges(6,3) * t271 + t241 * t209 - t240 * t212;
t284 = mrSges(5,1) * t164 - mrSges(5,3) * t162 - Ifges(5,4) * t202 - Ifges(5,2) * t271 - Ifges(5,6) * t201 + pkin(4) * t147 + t241 * t207 - t240 * t213 + t292;
t282 = mrSges(4,2) * t167 - t240 * t214 - pkin(3) * (-t202 * mrSges(5,2) - t241 * t221 - t147 + t316) - qJ(4) * (-t201 * mrSges(5,2) - t240 * t221 + t148 + t317) - mrSges(4,1) * t166 - t241 * t211 + Ifges(4,6) * t201 - Ifges(4,5) * t202 - Ifges(4,3) * t271 + t284;
t315 = mrSges(3,1) * t226 - mrSges(3,2) * t227 + Ifges(3,5) * t249 + Ifges(3,6) * t250 + Ifges(3,3) * qJDD(2) + pkin(2) * t133 + (t276 * t238 - t278 * t239) * qJD(1) - t282;
t307 = t272 * t209;
t210 = Ifges(5,4) * t241 + Ifges(5,2) * t272 + Ifges(5,6) * t240;
t305 = -Ifges(4,5) * t241 + Ifges(4,6) * t240 - Ifges(4,3) * t272 - t210;
t248 = (-mrSges(3,1) * t278 + mrSges(3,2) * t276) * qJD(1);
t253 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t302;
t131 = m(3) * t226 + qJDD(2) * mrSges(3,1) - t249 * mrSges(3,3) + qJD(2) * t253 - t248 * t303 + t133;
t252 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t303;
t297 = -t275 * t137 + t312 * t140;
t132 = m(3) * t227 - qJDD(2) * mrSges(3,2) + t250 * mrSges(3,3) - qJD(2) * t252 + t248 * t302 + t297;
t298 = -t276 * t131 + t278 * t132;
t151 = -t236 * qJ(5) + qJDD(5) + (-pkin(3) - pkin(4)) * t201 + (-pkin(3) * t272 + t230 + t313) * t241 - t319;
t146 = -m(6) * t151 + t201 * mrSges(6,1) - t202 * mrSges(6,2) + t240 * t228 - t241 * t231;
t291 = mrSges(6,1) * t151 - mrSges(6,3) * t157 - Ifges(6,4) * t202 - Ifges(6,2) * t201 + Ifges(6,6) * t271 + t272 * t212;
t206 = Ifges(6,5) * t241 + Ifges(6,6) * t240 - Ifges(6,3) * t272;
t290 = mrSges(6,2) * t151 - mrSges(6,3) * t153 + Ifges(6,1) * t202 + Ifges(6,4) * t201 - Ifges(6,5) * t271 + t240 * t206;
t159 = qJD(4) * t314 + (t241 * t272 + t201) * pkin(3) + t319;
t141 = m(5) * t159 + t201 * mrSges(5,1) - t202 * mrSges(5,3) - t241 * t233 + t240 * t234 + t146;
t287 = mrSges(5,1) * t159 - mrSges(5,2) * t162 + pkin(4) * t146 + qJ(5) * t148 - t291;
t128 = (Ifges(4,6) - Ifges(5,6)) * t271 - mrSges(4,1) * t203 + mrSges(4,3) * t167 - pkin(3) * t141 + (Ifges(4,4) - Ifges(5,5)) * t202 + (-Ifges(4,2) - Ifges(5,3)) * t201 + (t206 + t305) * t241 - t287 + (t214 + t213) * t272;
t286 = mrSges(5,2) * t164 - mrSges(5,3) * t159 + Ifges(5,1) * t202 + Ifges(5,4) * t271 + Ifges(5,5) * t201 - qJ(5) * t147 + t272 * t207 + t290;
t129 = (-t211 + t209) * t272 + t305 * t240 + Ifges(4,5) * t271 + mrSges(4,2) * t203 - Ifges(4,4) * t201 + Ifges(4,1) * t202 - mrSges(4,3) * t166 - qJ(4) * t141 + t286;
t237 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t276 + Ifges(3,6) * t278) * qJD(1);
t242 = -t280 * pkin(6) - t294;
t285 = m(4) * t203 + t201 * mrSges(4,1) + t202 * mrSges(4,2) + t240 * t229 + t241 * t232 + t141;
t122 = -mrSges(3,1) * t242 + mrSges(3,3) * t227 + Ifges(3,4) * t249 + Ifges(3,2) * t250 + Ifges(3,6) * qJDD(2) - pkin(2) * t285 + pkin(7) * t297 + qJD(2) * t239 + t312 * t128 + t275 * t129 - t237 * t303;
t124 = mrSges(3,2) * t242 - mrSges(3,3) * t226 + Ifges(3,1) * t249 + Ifges(3,4) * t250 + Ifges(3,5) * qJDD(2) - pkin(7) * t133 - qJD(2) * t238 - t275 * t128 + t312 * t129 + t237 * t302;
t283 = -m(3) * t242 + t250 * mrSges(3,1) - t249 * mrSges(3,2) - t252 * t303 + t253 * t302 - t285;
t288 = mrSges(2,1) * t255 - mrSges(2,2) * t256 + Ifges(2,3) * qJDD(1) + pkin(1) * t283 + pkin(6) * t298 + t278 * t122 + t276 * t124;
t134 = m(2) * t255 + qJDD(1) * mrSges(2,1) - t280 * mrSges(2,2) + t283;
t127 = t278 * t131 + t276 * t132;
t125 = m(2) * t256 - t280 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t298;
t120 = mrSges(2,1) * g(3) + mrSges(2,3) * t256 + t280 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t127 - t315;
t119 = -mrSges(2,2) * g(3) - mrSges(2,3) * t255 + Ifges(2,5) * qJDD(1) - t280 * Ifges(2,6) - pkin(6) * t127 - t276 * t122 + t278 * t124;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t279 * t119 - t277 * t120 - pkin(5) * (t277 * t125 + t279 * t134), t119, t124, t129, -t240 * t210 + t286 + t307, t290 + t307; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t277 * t119 + t279 * t120 + pkin(5) * (t279 * t125 - t277 * t134), t120, t122, t128, -t284, -t241 * t206 - t291; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t288, t288, t315, -t282, -t272 * t213 + Ifges(5,6) * t271 + Ifges(5,3) * t201 + Ifges(5,5) * t202 + (-t206 + t210) * t241 + t287, t292;];
m_new = t1;
