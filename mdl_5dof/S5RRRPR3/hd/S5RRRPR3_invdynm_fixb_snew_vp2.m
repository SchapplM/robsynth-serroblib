% Calculate vector of cutting torques with Newton-Euler for
% S5RRRPR3
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:27
% EndTime: 2019-12-05 18:42:36
% DurationCPUTime: 6.20s
% Computational Cost: add. (147781->269), mult. (197708->345), div. (0->0), fcn. (125769->10), ass. (0->111)
t231 = qJDD(1) + qJDD(2);
t240 = sin(qJ(3));
t244 = cos(qJ(3));
t233 = qJD(1) + qJD(2);
t262 = qJD(3) * t233;
t211 = t240 * t231 + t244 * t262;
t242 = sin(qJ(1));
t246 = cos(qJ(1));
t223 = t246 * g(2) + t242 * g(3);
t216 = qJDD(1) * pkin(1) + t223;
t222 = t242 * g(2) - t246 * g(3);
t247 = qJD(1) ^ 2;
t217 = -t247 * pkin(1) + t222;
t241 = sin(qJ(2));
t245 = cos(qJ(2));
t195 = t241 * t216 + t245 * t217;
t229 = t233 ^ 2;
t187 = -t229 * pkin(2) + t231 * pkin(7) + t195;
t263 = t240 * t187;
t266 = pkin(3) * t229;
t170 = qJDD(3) * pkin(3) - t211 * qJ(4) - t263 + (qJ(4) * t262 + t240 * t266 - g(1)) * t244;
t177 = -t240 * g(1) + t244 * t187;
t212 = t244 * t231 - t240 * t262;
t265 = t233 * t240;
t218 = qJD(3) * pkin(3) - qJ(4) * t265;
t236 = t244 ^ 2;
t171 = t212 * qJ(4) - qJD(3) * t218 - t236 * t266 + t177;
t237 = sin(pkin(9));
t238 = cos(pkin(9));
t203 = (t237 * t244 + t238 * t240) * t233;
t150 = -0.2e1 * qJD(4) * t203 + t238 * t170 - t237 * t171;
t192 = t238 * t211 + t237 * t212;
t202 = (-t237 * t240 + t238 * t244) * t233;
t147 = (qJD(3) * t202 - t192) * pkin(8) + (t202 * t203 + qJDD(3)) * pkin(4) + t150;
t151 = 0.2e1 * qJD(4) * t202 + t237 * t170 + t238 * t171;
t191 = -t237 * t211 + t238 * t212;
t198 = qJD(3) * pkin(4) - t203 * pkin(8);
t201 = t202 ^ 2;
t148 = -t201 * pkin(4) + t191 * pkin(8) - qJD(3) * t198 + t151;
t239 = sin(qJ(5));
t243 = cos(qJ(5));
t145 = t243 * t147 - t239 * t148;
t181 = t243 * t202 - t239 * t203;
t160 = t181 * qJD(5) + t239 * t191 + t243 * t192;
t182 = t239 * t202 + t243 * t203;
t166 = -t181 * mrSges(6,1) + t182 * mrSges(6,2);
t232 = qJD(3) + qJD(5);
t174 = -t232 * mrSges(6,2) + t181 * mrSges(6,3);
t230 = qJDD(3) + qJDD(5);
t142 = m(6) * t145 + t230 * mrSges(6,1) - t160 * mrSges(6,3) - t182 * t166 + t232 * t174;
t146 = t239 * t147 + t243 * t148;
t159 = -t182 * qJD(5) + t243 * t191 - t239 * t192;
t175 = t232 * mrSges(6,1) - t182 * mrSges(6,3);
t143 = m(6) * t146 - t230 * mrSges(6,2) + t159 * mrSges(6,3) + t181 * t166 - t232 * t175;
t133 = t243 * t142 + t239 * t143;
t185 = -t202 * mrSges(5,1) + t203 * mrSges(5,2);
t196 = -qJD(3) * mrSges(5,2) + t202 * mrSges(5,3);
t130 = m(5) * t150 + qJDD(3) * mrSges(5,1) - t192 * mrSges(5,3) + qJD(3) * t196 - t203 * t185 + t133;
t197 = qJD(3) * mrSges(5,1) - t203 * mrSges(5,3);
t258 = -t239 * t142 + t243 * t143;
t131 = m(5) * t151 - qJDD(3) * mrSges(5,2) + t191 * mrSges(5,3) - qJD(3) * t197 + t202 * t185 + t258;
t126 = t238 * t130 + t237 * t131;
t176 = -t244 * g(1) - t263;
t205 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t240 + Ifges(4,2) * t244) * t233;
t206 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t240 + Ifges(4,4) * t244) * t233;
t179 = Ifges(5,4) * t203 + Ifges(5,2) * t202 + Ifges(5,6) * qJD(3);
t180 = Ifges(5,1) * t203 + Ifges(5,4) * t202 + Ifges(5,5) * qJD(3);
t162 = Ifges(6,4) * t182 + Ifges(6,2) * t181 + Ifges(6,6) * t232;
t163 = Ifges(6,1) * t182 + Ifges(6,4) * t181 + Ifges(6,5) * t232;
t253 = -mrSges(6,1) * t145 + mrSges(6,2) * t146 - Ifges(6,5) * t160 - Ifges(6,6) * t159 - Ifges(6,3) * t230 - t182 * t162 + t181 * t163;
t250 = -mrSges(5,1) * t150 + mrSges(5,2) * t151 - Ifges(5,5) * t192 - Ifges(5,6) * t191 - Ifges(5,3) * qJDD(3) - pkin(4) * t133 - t203 * t179 + t202 * t180 + t253;
t267 = mrSges(4,1) * t176 - mrSges(4,2) * t177 + Ifges(4,5) * t211 + Ifges(4,6) * t212 + Ifges(4,3) * qJDD(3) + pkin(3) * t126 + (t240 * t205 - t244 * t206) * t233 - t250;
t264 = t233 * t244;
t210 = (-mrSges(4,1) * t244 + mrSges(4,2) * t240) * t233;
t220 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t264;
t124 = m(4) * t176 + qJDD(3) * mrSges(4,1) - t211 * mrSges(4,3) + qJD(3) * t220 - t210 * t265 + t126;
t219 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t259 = -t237 * t130 + t238 * t131;
t125 = m(4) * t177 - qJDD(3) * mrSges(4,2) + t212 * mrSges(4,3) - qJD(3) * t219 + t210 * t264 + t259;
t260 = -t240 * t124 + t244 * t125;
t116 = m(3) * t195 - t229 * mrSges(3,1) - t231 * mrSges(3,2) + t260;
t194 = t245 * t216 - t241 * t217;
t255 = -t231 * pkin(2) - t194;
t186 = -t229 * pkin(7) + t255;
t172 = -t212 * pkin(3) + qJDD(4) + t218 * t265 + (-qJ(4) * t236 - pkin(7)) * t229 + t255;
t153 = -t191 * pkin(4) - t201 * pkin(8) + t203 * t198 + t172;
t257 = m(6) * t153 - t159 * mrSges(6,1) + t160 * mrSges(6,2) - t181 * t174 + t182 * t175;
t251 = m(5) * t172 - t191 * mrSges(5,1) + t192 * mrSges(5,2) - t202 * t196 + t203 * t197 + t257;
t249 = -m(4) * t186 + t212 * mrSges(4,1) - t211 * mrSges(4,2) - t219 * t265 + t220 * t264 - t251;
t137 = m(3) * t194 + t231 * mrSges(3,1) - t229 * mrSges(3,2) + t249;
t113 = t241 * t116 + t245 * t137;
t118 = t244 * t124 + t240 * t125;
t261 = t245 * t116 - t241 * t137;
t161 = Ifges(6,5) * t182 + Ifges(6,6) * t181 + Ifges(6,3) * t232;
t134 = -mrSges(6,1) * t153 + mrSges(6,3) * t146 + Ifges(6,4) * t160 + Ifges(6,2) * t159 + Ifges(6,6) * t230 - t182 * t161 + t232 * t163;
t135 = mrSges(6,2) * t153 - mrSges(6,3) * t145 + Ifges(6,1) * t160 + Ifges(6,4) * t159 + Ifges(6,5) * t230 + t181 * t161 - t232 * t162;
t178 = Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * qJD(3);
t119 = -mrSges(5,1) * t172 + mrSges(5,3) * t151 + Ifges(5,4) * t192 + Ifges(5,2) * t191 + Ifges(5,6) * qJDD(3) - pkin(4) * t257 + pkin(8) * t258 + qJD(3) * t180 + t243 * t134 + t239 * t135 - t203 * t178;
t120 = mrSges(5,2) * t172 - mrSges(5,3) * t150 + Ifges(5,1) * t192 + Ifges(5,4) * t191 + Ifges(5,5) * qJDD(3) - pkin(8) * t133 - qJD(3) * t179 - t239 * t134 + t243 * t135 + t202 * t178;
t204 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t240 + Ifges(4,6) * t244) * t233;
t107 = -mrSges(4,1) * t186 + mrSges(4,3) * t177 + Ifges(4,4) * t211 + Ifges(4,2) * t212 + Ifges(4,6) * qJDD(3) - pkin(3) * t251 + qJ(4) * t259 + qJD(3) * t206 + t238 * t119 + t237 * t120 - t204 * t265;
t109 = mrSges(4,2) * t186 - mrSges(4,3) * t176 + Ifges(4,1) * t211 + Ifges(4,4) * t212 + Ifges(4,5) * qJDD(3) - qJ(4) * t126 - qJD(3) * t205 - t237 * t119 + t238 * t120 + t204 * t264;
t254 = mrSges(3,1) * t194 - mrSges(3,2) * t195 + Ifges(3,3) * t231 + pkin(2) * t249 + pkin(7) * t260 + t244 * t107 + t240 * t109;
t252 = mrSges(2,1) * t223 - mrSges(2,2) * t222 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t254;
t111 = m(2) * t223 + qJDD(1) * mrSges(2,1) - t247 * mrSges(2,2) + t113;
t110 = m(2) * t222 - t247 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t261;
t105 = mrSges(3,1) * g(1) + mrSges(3,3) * t195 + t229 * Ifges(3,5) + Ifges(3,6) * t231 - pkin(2) * t118 - t267;
t104 = -mrSges(3,2) * g(1) - mrSges(3,3) * t194 + Ifges(3,5) * t231 - t229 * Ifges(3,6) - pkin(7) * t118 - t240 * t107 + t244 * t109;
t103 = -mrSges(2,2) * g(1) - mrSges(2,3) * t223 + Ifges(2,5) * qJDD(1) - t247 * Ifges(2,6) - pkin(6) * t113 + t245 * t104 - t241 * t105;
t102 = Ifges(2,6) * qJDD(1) + t247 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t222 + t241 * t104 + t245 * t105 - pkin(1) * (-m(3) * g(1) + t118) + pkin(6) * t261;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t252, t103, t104, t109, t120, t135; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t242 * t103 - t246 * t102 - pkin(5) * (t246 * t110 - t242 * t111), t102, t105, t107, t119, t134; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t246 * t103 - t242 * t102 + pkin(5) * (-t242 * t110 - t246 * t111), t252, t254, t267, -t250, -t253;];
m_new = t1;
