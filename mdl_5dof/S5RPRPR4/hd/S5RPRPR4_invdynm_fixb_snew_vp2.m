% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:19
% EndTime: 2019-12-05 17:53:30
% DurationCPUTime: 7.16s
% Computational Cost: add. (91738->267), mult. (197708->343), div. (0->0), fcn. (125769->10), ass. (0->108)
t243 = sin(qJ(1));
t246 = cos(qJ(1));
t223 = t246 * g(2) + t243 * g(3);
t213 = qJDD(1) * pkin(1) + t223;
t222 = t243 * g(2) - g(3) * t246;
t247 = qJD(1) ^ 2;
t215 = -pkin(1) * t247 + t222;
t238 = sin(pkin(8));
t240 = cos(pkin(8));
t193 = t238 * t213 + t240 * t215;
t185 = -pkin(2) * t247 + qJDD(1) * pkin(6) + t193;
t236 = -g(1) + qJDD(2);
t242 = sin(qJ(3));
t245 = cos(qJ(3));
t174 = -t242 * t185 + t245 * t236;
t263 = qJD(1) * qJD(3);
t262 = t245 * t263;
t216 = qJDD(1) * t242 + t262;
t170 = (-t216 + t262) * qJ(4) + (t242 * t245 * t247 + qJDD(3)) * pkin(3) + t174;
t175 = t245 * t185 + t242 * t236;
t217 = qJDD(1) * t245 - t242 * t263;
t265 = qJD(1) * t242;
t219 = qJD(3) * pkin(3) - qJ(4) * t265;
t235 = t245 ^ 2;
t171 = -pkin(3) * t235 * t247 + qJ(4) * t217 - qJD(3) * t219 + t175;
t237 = sin(pkin(9));
t239 = cos(pkin(9));
t203 = (t237 * t245 + t239 * t242) * qJD(1);
t150 = -0.2e1 * qJD(4) * t203 + t239 * t170 - t237 * t171;
t195 = t216 * t239 + t217 * t237;
t202 = (-t237 * t242 + t239 * t245) * qJD(1);
t147 = (qJD(3) * t202 - t195) * pkin(7) + (t202 * t203 + qJDD(3)) * pkin(4) + t150;
t151 = 0.2e1 * qJD(4) * t202 + t237 * t170 + t239 * t171;
t194 = -t216 * t237 + t217 * t239;
t198 = qJD(3) * pkin(4) - pkin(7) * t203;
t201 = t202 ^ 2;
t148 = -pkin(4) * t201 + pkin(7) * t194 - qJD(3) * t198 + t151;
t241 = sin(qJ(5));
t244 = cos(qJ(5));
t145 = t147 * t244 - t148 * t241;
t182 = t202 * t244 - t203 * t241;
t160 = qJD(5) * t182 + t194 * t241 + t195 * t244;
t183 = t202 * t241 + t203 * t244;
t169 = -mrSges(6,1) * t182 + mrSges(6,2) * t183;
t231 = qJD(3) + qJD(5);
t176 = -mrSges(6,2) * t231 + mrSges(6,3) * t182;
t230 = qJDD(3) + qJDD(5);
t142 = m(6) * t145 + mrSges(6,1) * t230 - mrSges(6,3) * t160 - t169 * t183 + t176 * t231;
t146 = t147 * t241 + t148 * t244;
t159 = -qJD(5) * t183 + t194 * t244 - t195 * t241;
t177 = mrSges(6,1) * t231 - mrSges(6,3) * t183;
t143 = m(6) * t146 - mrSges(6,2) * t230 + mrSges(6,3) * t159 + t169 * t182 - t177 * t231;
t133 = t244 * t142 + t241 * t143;
t187 = -mrSges(5,1) * t202 + mrSges(5,2) * t203;
t196 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t202;
t130 = m(5) * t150 + qJDD(3) * mrSges(5,1) - mrSges(5,3) * t195 + qJD(3) * t196 - t187 * t203 + t133;
t197 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t203;
t258 = -t142 * t241 + t244 * t143;
t131 = m(5) * t151 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t194 - qJD(3) * t197 + t187 * t202 + t258;
t126 = t239 * t130 + t237 * t131;
t208 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t242 + Ifges(4,2) * t245) * qJD(1);
t209 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t242 + Ifges(4,4) * t245) * qJD(1);
t180 = Ifges(5,4) * t203 + Ifges(5,2) * t202 + Ifges(5,6) * qJD(3);
t181 = Ifges(5,1) * t203 + Ifges(5,4) * t202 + Ifges(5,5) * qJD(3);
t162 = Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t231;
t163 = Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t231;
t253 = -mrSges(6,1) * t145 + mrSges(6,2) * t146 - Ifges(6,5) * t160 - Ifges(6,6) * t159 - Ifges(6,3) * t230 - t183 * t162 + t182 * t163;
t250 = -mrSges(5,1) * t150 + mrSges(5,2) * t151 - Ifges(5,5) * t195 - Ifges(5,6) * t194 - Ifges(5,3) * qJDD(3) - pkin(4) * t133 - t203 * t180 + t202 * t181 + t253;
t266 = mrSges(4,1) * t174 - mrSges(4,2) * t175 + Ifges(4,5) * t216 + Ifges(4,6) * t217 + Ifges(4,3) * qJDD(3) + pkin(3) * t126 + (t208 * t242 - t209 * t245) * qJD(1) - t250;
t214 = (-mrSges(4,1) * t245 + mrSges(4,2) * t242) * qJD(1);
t264 = qJD(1) * t245;
t221 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t264;
t124 = m(4) * t174 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t216 + qJD(3) * t221 - t214 * t265 + t126;
t220 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t259 = -t130 * t237 + t239 * t131;
t125 = m(4) * t175 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t217 - qJD(3) * t220 + t214 * t264 + t259;
t260 = -t124 * t242 + t245 * t125;
t116 = m(3) * t193 - mrSges(3,1) * t247 - qJDD(1) * mrSges(3,2) + t260;
t192 = t240 * t213 - t238 * t215;
t255 = -qJDD(1) * pkin(2) - t192;
t184 = -pkin(6) * t247 + t255;
t172 = -t217 * pkin(3) + qJDD(4) + t219 * t265 + (-qJ(4) * t235 - pkin(6)) * t247 + t255;
t153 = -pkin(4) * t194 - pkin(7) * t201 + t198 * t203 + t172;
t257 = m(6) * t153 - t159 * mrSges(6,1) + t160 * mrSges(6,2) - t182 * t176 + t183 * t177;
t251 = m(5) * t172 - t194 * mrSges(5,1) + mrSges(5,2) * t195 - t202 * t196 + t197 * t203 + t257;
t249 = -m(4) * t184 + t217 * mrSges(4,1) - mrSges(4,2) * t216 - t220 * t265 + t221 * t264 - t251;
t137 = m(3) * t192 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t247 + t249;
t113 = t238 * t116 + t240 * t137;
t118 = t245 * t124 + t242 * t125;
t261 = t240 * t116 - t238 * t137;
t161 = Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t231;
t134 = -mrSges(6,1) * t153 + mrSges(6,3) * t146 + Ifges(6,4) * t160 + Ifges(6,2) * t159 + Ifges(6,6) * t230 - t161 * t183 + t163 * t231;
t135 = mrSges(6,2) * t153 - mrSges(6,3) * t145 + Ifges(6,1) * t160 + Ifges(6,4) * t159 + Ifges(6,5) * t230 + t161 * t182 - t162 * t231;
t179 = Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * qJD(3);
t119 = -mrSges(5,1) * t172 + mrSges(5,3) * t151 + Ifges(5,4) * t195 + Ifges(5,2) * t194 + Ifges(5,6) * qJDD(3) - pkin(4) * t257 + pkin(7) * t258 + qJD(3) * t181 + t244 * t134 + t241 * t135 - t203 * t179;
t120 = mrSges(5,2) * t172 - mrSges(5,3) * t150 + Ifges(5,1) * t195 + Ifges(5,4) * t194 + Ifges(5,5) * qJDD(3) - pkin(7) * t133 - qJD(3) * t180 - t134 * t241 + t135 * t244 + t179 * t202;
t207 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t242 + Ifges(4,6) * t245) * qJD(1);
t107 = -mrSges(4,1) * t184 + mrSges(4,3) * t175 + Ifges(4,4) * t216 + Ifges(4,2) * t217 + Ifges(4,6) * qJDD(3) - pkin(3) * t251 + qJ(4) * t259 + qJD(3) * t209 + t239 * t119 + t237 * t120 - t207 * t265;
t109 = mrSges(4,2) * t184 - mrSges(4,3) * t174 + Ifges(4,1) * t216 + Ifges(4,4) * t217 + Ifges(4,5) * qJDD(3) - qJ(4) * t126 - qJD(3) * t208 - t119 * t237 + t120 * t239 + t207 * t264;
t254 = mrSges(3,1) * t192 - mrSges(3,2) * t193 + Ifges(3,3) * qJDD(1) + pkin(2) * t249 + pkin(6) * t260 + t245 * t107 + t242 * t109;
t252 = mrSges(2,1) * t223 - mrSges(2,2) * t222 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t254;
t111 = m(2) * t223 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t247 + t113;
t110 = m(2) * t222 - mrSges(2,1) * t247 - qJDD(1) * mrSges(2,2) + t261;
t105 = -mrSges(3,1) * t236 + mrSges(3,3) * t193 + t247 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t118 - t266;
t104 = mrSges(3,2) * t236 - mrSges(3,3) * t192 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t247 - pkin(6) * t118 - t107 * t242 + t109 * t245;
t103 = -mrSges(2,2) * g(1) - mrSges(2,3) * t223 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t247 - qJ(2) * t113 + t104 * t240 - t105 * t238;
t102 = Ifges(2,6) * qJDD(1) + t247 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t222 + t238 * t104 + t240 * t105 - pkin(1) * (m(3) * t236 + t118) + qJ(2) * t261;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t252, t103, t104, t109, t120, t135; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t243 * t103 - t246 * t102 - pkin(5) * (t110 * t246 - t111 * t243), t102, t105, t107, t119, t134; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t246 * t103 - t243 * t102 + pkin(5) * (-t110 * t243 - t111 * t246), t252, t254, t266, -t250, -t253;];
m_new = t1;
