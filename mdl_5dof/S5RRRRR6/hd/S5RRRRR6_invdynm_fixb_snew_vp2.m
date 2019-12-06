% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 19:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 19:00:12
% EndTime: 2019-12-05 19:00:21
% DurationCPUTime: 6.43s
% Computational Cost: add. (160198->270), mult. (205381->343), div. (0->0), fcn. (133069->10), ass. (0->113)
t232 = qJDD(1) + qJDD(2);
t239 = sin(qJ(3));
t244 = cos(qJ(3));
t234 = qJD(1) + qJD(2);
t262 = qJD(3) * t234;
t209 = t232 * t239 + t244 * t262;
t241 = sin(qJ(1));
t246 = cos(qJ(1));
t221 = t246 * g(2) + t241 * g(3);
t214 = qJDD(1) * pkin(1) + t221;
t220 = t241 * g(2) - g(3) * t246;
t247 = qJD(1) ^ 2;
t215 = -pkin(1) * t247 + t220;
t240 = sin(qJ(2));
t245 = cos(qJ(2));
t195 = t240 * t214 + t245 * t215;
t230 = t234 ^ 2;
t192 = -pkin(2) * t230 + pkin(7) * t232 + t195;
t263 = t239 * t192;
t266 = pkin(3) * t230;
t170 = qJDD(3) * pkin(3) - t209 * pkin(8) - t263 + (pkin(8) * t262 + t239 * t266 - g(1)) * t244;
t182 = -g(1) * t239 + t244 * t192;
t210 = t232 * t244 - t239 * t262;
t265 = t234 * t239;
t218 = qJD(3) * pkin(3) - pkin(8) * t265;
t236 = t244 ^ 2;
t171 = pkin(8) * t210 - qJD(3) * t218 - t236 * t266 + t182;
t238 = sin(qJ(4));
t243 = cos(qJ(4));
t152 = t243 * t170 - t238 * t171;
t203 = (-t238 * t239 + t243 * t244) * t234;
t178 = qJD(4) * t203 + t209 * t243 + t210 * t238;
t204 = (t238 * t244 + t239 * t243) * t234;
t231 = qJDD(3) + qJDD(4);
t233 = qJD(3) + qJD(4);
t147 = (t203 * t233 - t178) * pkin(9) + (t203 * t204 + t231) * pkin(4) + t152;
t153 = t238 * t170 + t243 * t171;
t177 = -qJD(4) * t204 - t209 * t238 + t210 * t243;
t198 = pkin(4) * t233 - pkin(9) * t204;
t199 = t203 ^ 2;
t148 = -pkin(4) * t199 + pkin(9) * t177 - t198 * t233 + t153;
t237 = sin(qJ(5));
t242 = cos(qJ(5));
t145 = t147 * t242 - t148 * t237;
t187 = t203 * t242 - t204 * t237;
t159 = qJD(5) * t187 + t177 * t237 + t178 * t242;
t188 = t203 * t237 + t204 * t242;
t166 = -mrSges(6,1) * t187 + mrSges(6,2) * t188;
t226 = qJD(5) + t233;
t179 = -mrSges(6,2) * t226 + mrSges(6,3) * t187;
t225 = qJDD(5) + t231;
t142 = m(6) * t145 + mrSges(6,1) * t225 - mrSges(6,3) * t159 - t166 * t188 + t179 * t226;
t146 = t147 * t237 + t148 * t242;
t158 = -qJD(5) * t188 + t177 * t242 - t178 * t237;
t180 = mrSges(6,1) * t226 - mrSges(6,3) * t188;
t143 = m(6) * t146 - mrSges(6,2) * t225 + mrSges(6,3) * t158 + t166 * t187 - t180 * t226;
t133 = t142 * t242 + t143 * t237;
t190 = -mrSges(5,1) * t203 + mrSges(5,2) * t204;
t196 = -mrSges(5,2) * t233 + mrSges(5,3) * t203;
t130 = m(5) * t152 + mrSges(5,1) * t231 - mrSges(5,3) * t178 - t190 * t204 + t196 * t233 + t133;
t197 = mrSges(5,1) * t233 - mrSges(5,3) * t204;
t258 = -t142 * t237 + t143 * t242;
t131 = m(5) * t153 - mrSges(5,2) * t231 + mrSges(5,3) * t177 + t190 * t203 - t197 * t233 + t258;
t126 = t130 * t243 + t131 * t238;
t181 = -t244 * g(1) - t263;
t201 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t239 + Ifges(4,2) * t244) * t234;
t202 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t239 + Ifges(4,4) * t244) * t234;
t184 = Ifges(5,4) * t204 + Ifges(5,2) * t203 + Ifges(5,6) * t233;
t185 = Ifges(5,1) * t204 + Ifges(5,4) * t203 + Ifges(5,5) * t233;
t162 = Ifges(6,4) * t188 + Ifges(6,2) * t187 + Ifges(6,6) * t226;
t163 = Ifges(6,1) * t188 + Ifges(6,4) * t187 + Ifges(6,5) * t226;
t253 = -mrSges(6,1) * t145 + mrSges(6,2) * t146 - Ifges(6,5) * t159 - Ifges(6,6) * t158 - Ifges(6,3) * t225 - t162 * t188 + t187 * t163;
t250 = -mrSges(5,1) * t152 + mrSges(5,2) * t153 - Ifges(5,5) * t178 - Ifges(5,6) * t177 - Ifges(5,3) * t231 - pkin(4) * t133 - t204 * t184 + t203 * t185 + t253;
t267 = mrSges(4,1) * t181 - mrSges(4,2) * t182 + Ifges(4,5) * t209 + Ifges(4,6) * t210 + Ifges(4,3) * qJDD(3) + pkin(3) * t126 + (t201 * t239 - t202 * t244) * t234 - t250;
t264 = t234 * t244;
t208 = (-mrSges(4,1) * t244 + mrSges(4,2) * t239) * t234;
t217 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t264;
t124 = m(4) * t181 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t209 + qJD(3) * t217 - t208 * t265 + t126;
t216 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t259 = -t130 * t238 + t131 * t243;
t125 = m(4) * t182 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t210 - qJD(3) * t216 + t208 * t264 + t259;
t260 = -t124 * t239 + t125 * t244;
t116 = m(3) * t195 - mrSges(3,1) * t230 - mrSges(3,2) * t232 + t260;
t194 = t245 * t214 - t240 * t215;
t255 = -t232 * pkin(2) - t194;
t191 = -pkin(7) * t230 + t255;
t172 = -t210 * pkin(3) + t218 * t265 + (-pkin(8) * t236 - pkin(7)) * t230 + t255;
t150 = -t177 * pkin(4) - t199 * pkin(9) + t204 * t198 + t172;
t257 = m(6) * t150 - mrSges(6,1) * t158 + mrSges(6,2) * t159 - t187 * t179 + t188 * t180;
t251 = m(5) * t172 - t177 * mrSges(5,1) + mrSges(5,2) * t178 - t203 * t196 + t197 * t204 + t257;
t249 = -m(4) * t191 + t210 * mrSges(4,1) - mrSges(4,2) * t209 - t216 * t265 + t217 * t264 - t251;
t137 = m(3) * t194 + mrSges(3,1) * t232 - mrSges(3,2) * t230 + t249;
t113 = t116 * t240 + t137 * t245;
t118 = t124 * t244 + t125 * t239;
t261 = t116 * t245 - t137 * t240;
t161 = Ifges(6,5) * t188 + Ifges(6,6) * t187 + Ifges(6,3) * t226;
t134 = -mrSges(6,1) * t150 + mrSges(6,3) * t146 + Ifges(6,4) * t159 + Ifges(6,2) * t158 + Ifges(6,6) * t225 - t161 * t188 + t163 * t226;
t135 = mrSges(6,2) * t150 - mrSges(6,3) * t145 + Ifges(6,1) * t159 + Ifges(6,4) * t158 + Ifges(6,5) * t225 + t161 * t187 - t162 * t226;
t183 = Ifges(5,5) * t204 + Ifges(5,6) * t203 + Ifges(5,3) * t233;
t119 = -mrSges(5,1) * t172 + mrSges(5,3) * t153 + Ifges(5,4) * t178 + Ifges(5,2) * t177 + Ifges(5,6) * t231 - pkin(4) * t257 + pkin(9) * t258 + t242 * t134 + t237 * t135 - t204 * t183 + t233 * t185;
t120 = mrSges(5,2) * t172 - mrSges(5,3) * t152 + Ifges(5,1) * t178 + Ifges(5,4) * t177 + Ifges(5,5) * t231 - pkin(9) * t133 - t134 * t237 + t135 * t242 + t183 * t203 - t184 * t233;
t200 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t239 + Ifges(4,6) * t244) * t234;
t107 = -mrSges(4,1) * t191 + mrSges(4,3) * t182 + Ifges(4,4) * t209 + Ifges(4,2) * t210 + Ifges(4,6) * qJDD(3) - pkin(3) * t251 + pkin(8) * t259 + qJD(3) * t202 + t243 * t119 + t238 * t120 - t200 * t265;
t109 = mrSges(4,2) * t191 - mrSges(4,3) * t181 + Ifges(4,1) * t209 + Ifges(4,4) * t210 + Ifges(4,5) * qJDD(3) - pkin(8) * t126 - qJD(3) * t201 - t119 * t238 + t120 * t243 + t200 * t264;
t254 = mrSges(3,1) * t194 - mrSges(3,2) * t195 + Ifges(3,3) * t232 + pkin(2) * t249 + pkin(7) * t260 + t107 * t244 + t109 * t239;
t252 = mrSges(2,1) * t221 - mrSges(2,2) * t220 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t254;
t111 = m(2) * t221 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t247 + t113;
t110 = m(2) * t220 - mrSges(2,1) * t247 - qJDD(1) * mrSges(2,2) + t261;
t105 = mrSges(3,1) * g(1) + mrSges(3,3) * t195 + t230 * Ifges(3,5) + Ifges(3,6) * t232 - pkin(2) * t118 - t267;
t104 = -mrSges(3,2) * g(1) - mrSges(3,3) * t194 + Ifges(3,5) * t232 - Ifges(3,6) * t230 - pkin(7) * t118 - t107 * t239 + t109 * t244;
t103 = -mrSges(2,2) * g(1) - mrSges(2,3) * t221 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t247 - pkin(6) * t113 + t104 * t245 - t105 * t240;
t102 = Ifges(2,6) * qJDD(1) + t247 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t220 + t240 * t104 + t245 * t105 - pkin(1) * (-m(3) * g(1) + t118) + pkin(6) * t261;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t252, t103, t104, t109, t120, t135; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t241 * t103 - t246 * t102 - pkin(5) * (t110 * t246 - t111 * t241), t102, t105, t107, t119, t134; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t246 * t103 - t241 * t102 + pkin(5) * (-t110 * t241 - t111 * t246), t252, t254, t267, -t250, -t253;];
m_new = t1;
