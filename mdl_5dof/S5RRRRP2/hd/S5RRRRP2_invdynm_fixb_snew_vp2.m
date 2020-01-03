% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:11:10
% EndTime: 2020-01-03 12:11:17
% DurationCPUTime: 3.82s
% Computational Cost: add. (63584->264), mult. (80564->325), div. (0->0), fcn. (47650->8), ass. (0->102)
t233 = qJD(1) + qJD(2);
t236 = sin(qJ(4));
t237 = sin(qJ(3));
t240 = cos(qJ(4));
t241 = cos(qJ(3));
t203 = (-t236 * t237 + t240 * t241) * t233;
t231 = qJDD(1) + qJDD(2);
t263 = qJD(3) * t233;
t209 = t237 * t231 + t241 * t263;
t210 = t241 * t231 - t237 * t263;
t172 = t203 * qJD(4) + t240 * t209 + t236 * t210;
t204 = (t236 * t241 + t237 * t240) * t233;
t186 = -t203 * mrSges(6,1) + t204 * mrSges(6,2);
t239 = sin(qJ(1));
t243 = cos(qJ(1));
t221 = -t243 * g(2) - t239 * g(3);
t214 = qJDD(1) * pkin(1) + t221;
t220 = -t239 * g(2) + t243 * g(3);
t244 = qJD(1) ^ 2;
t215 = -t244 * pkin(1) + t220;
t238 = sin(qJ(2));
t242 = cos(qJ(2));
t192 = t238 * t214 + t242 * t215;
t229 = t233 ^ 2;
t189 = -t229 * pkin(2) + t231 * pkin(7) + t192;
t265 = t237 * t189;
t268 = pkin(3) * t229;
t152 = qJDD(3) * pkin(3) - t209 * pkin(8) - t265 + (pkin(8) * t263 + t237 * t268 - g(1)) * t241;
t175 = -t237 * g(1) + t241 * t189;
t267 = t233 * t237;
t218 = qJD(3) * pkin(3) - pkin(8) * t267;
t235 = t241 ^ 2;
t153 = t210 * pkin(8) - qJD(3) * t218 - t235 * t268 + t175;
t147 = t240 * t152 - t236 * t153;
t230 = qJDD(3) + qJDD(4);
t232 = qJD(3) + qJD(4);
t139 = -0.2e1 * qJD(5) * t204 + (t203 * t232 - t172) * qJ(5) + (t203 * t204 + t230) * pkin(4) + t147;
t194 = -t232 * mrSges(6,2) + t203 * mrSges(6,3);
t262 = m(6) * t139 + t230 * mrSges(6,1) + t232 * t194;
t136 = -t172 * mrSges(6,3) - t204 * t186 + t262;
t148 = t236 * t152 + t240 * t153;
t171 = -t204 * qJD(4) - t236 * t209 + t240 * t210;
t179 = Ifges(5,4) * t204 + Ifges(5,2) * t203 + Ifges(5,6) * t232;
t180 = Ifges(6,1) * t204 + Ifges(6,4) * t203 + Ifges(6,5) * t232;
t181 = Ifges(5,1) * t204 + Ifges(5,4) * t203 + Ifges(5,5) * t232;
t196 = t232 * pkin(4) - t204 * qJ(5);
t199 = t203 ^ 2;
t142 = -t199 * pkin(4) + t171 * qJ(5) + 0.2e1 * qJD(5) * t203 - t232 * t196 + t148;
t178 = Ifges(6,4) * t204 + Ifges(6,2) * t203 + Ifges(6,6) * t232;
t252 = -mrSges(6,1) * t139 + mrSges(6,2) * t142 - Ifges(6,5) * t172 - Ifges(6,6) * t171 - Ifges(6,3) * t230 - t204 * t178;
t271 = mrSges(5,1) * t147 - mrSges(5,2) * t148 + Ifges(5,5) * t172 + Ifges(5,6) * t171 + Ifges(5,3) * t230 + pkin(4) * t136 + t204 * t179 - t252 + (-t181 - t180) * t203;
t187 = -t203 * mrSges(5,1) + t204 * mrSges(5,2);
t195 = -t232 * mrSges(5,2) + t203 * mrSges(5,3);
t131 = m(5) * t147 + t230 * mrSges(5,1) + t232 * t195 + (-t186 - t187) * t204 + (-mrSges(5,3) - mrSges(6,3)) * t172 + t262;
t197 = t232 * mrSges(6,1) - t204 * mrSges(6,3);
t198 = t232 * mrSges(5,1) - t204 * mrSges(5,3);
t261 = m(6) * t142 + t171 * mrSges(6,3) + t203 * t186;
t134 = m(5) * t148 + t171 * mrSges(5,3) + t203 * t187 + (-t197 - t198) * t232 + (-mrSges(5,2) - mrSges(6,2)) * t230 + t261;
t126 = t240 * t131 + t236 * t134;
t174 = -t241 * g(1) - t265;
t201 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t237 + Ifges(4,2) * t241) * t233;
t202 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t237 + Ifges(4,4) * t241) * t233;
t270 = mrSges(4,1) * t174 - mrSges(4,2) * t175 + Ifges(4,5) * t209 + Ifges(4,6) * t210 + Ifges(4,3) * qJDD(3) + pkin(3) * t126 + (t237 * t201 - t241 * t202) * t233 + t271;
t266 = t233 * t241;
t208 = (-mrSges(4,1) * t241 + mrSges(4,2) * t237) * t233;
t217 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t266;
t124 = m(4) * t174 + qJDD(3) * mrSges(4,1) - t209 * mrSges(4,3) + qJD(3) * t217 - t208 * t267 + t126;
t216 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t267;
t257 = -t236 * t131 + t240 * t134;
t125 = m(4) * t175 - qJDD(3) * mrSges(4,2) + t210 * mrSges(4,3) - qJD(3) * t216 + t208 * t266 + t257;
t258 = -t237 * t124 + t241 * t125;
t116 = m(3) * t192 - t229 * mrSges(3,1) - t231 * mrSges(3,2) + t258;
t191 = t242 * t214 - t238 * t215;
t254 = -t231 * pkin(2) - t191;
t188 = -t229 * pkin(7) + t254;
t154 = -t210 * pkin(3) + t218 * t267 + (-pkin(8) * t235 - pkin(7)) * t229 + t254;
t145 = -t171 * pkin(4) - t199 * qJ(5) + t204 * t196 + qJDD(5) + t154;
t256 = m(6) * t145 - t171 * mrSges(6,1) + t172 * mrSges(6,2) - t203 * t194 + t204 * t197;
t248 = m(5) * t154 - t171 * mrSges(5,1) + t172 * mrSges(5,2) - t203 * t195 + t204 * t198 + t256;
t246 = -m(4) * t188 + t210 * mrSges(4,1) - t209 * mrSges(4,2) - t216 * t267 + t217 * t266 - t248;
t128 = m(3) * t191 + t231 * mrSges(3,1) - t229 * mrSges(3,2) + t246;
t113 = t238 * t116 + t242 * t128;
t118 = t241 * t124 + t237 * t125;
t259 = t242 * t116 - t238 * t128;
t253 = -mrSges(6,1) * t145 + mrSges(6,3) * t142 + Ifges(6,4) * t172 + Ifges(6,2) * t171 + Ifges(6,6) * t230 + t232 * t180;
t176 = Ifges(6,5) * t204 + Ifges(6,6) * t203 + Ifges(6,3) * t232;
t251 = mrSges(6,2) * t145 - mrSges(6,3) * t139 + Ifges(6,1) * t172 + Ifges(6,4) * t171 + Ifges(6,5) * t230 + t203 * t176;
t177 = Ifges(5,5) * t204 + Ifges(5,6) * t203 + Ifges(5,3) * t232;
t119 = Ifges(5,4) * t172 + Ifges(5,2) * t171 + Ifges(5,6) * t230 + t232 * t181 - mrSges(5,1) * t154 + mrSges(5,3) * t148 - pkin(4) * t256 + qJ(5) * (-t230 * mrSges(6,2) - t232 * t197 + t261) + (-t177 - t176) * t204 + t253;
t120 = mrSges(5,2) * t154 - mrSges(5,3) * t147 + Ifges(5,1) * t172 + Ifges(5,4) * t171 + Ifges(5,5) * t230 - qJ(5) * t136 + t203 * t177 + (-t178 - t179) * t232 + t251;
t200 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t237 + Ifges(4,6) * t241) * t233;
t107 = -mrSges(4,1) * t188 + mrSges(4,3) * t175 + Ifges(4,4) * t209 + Ifges(4,2) * t210 + Ifges(4,6) * qJDD(3) - pkin(3) * t248 + pkin(8) * t257 + qJD(3) * t202 + t240 * t119 + t236 * t120 - t200 * t267;
t109 = mrSges(4,2) * t188 - mrSges(4,3) * t174 + Ifges(4,1) * t209 + Ifges(4,4) * t210 + Ifges(4,5) * qJDD(3) - pkin(8) * t126 - qJD(3) * t201 - t236 * t119 + t240 * t120 + t200 * t266;
t250 = mrSges(3,1) * t191 - mrSges(3,2) * t192 + Ifges(3,3) * t231 + pkin(2) * t246 + pkin(7) * t258 + t241 * t107 + t237 * t109;
t249 = mrSges(2,1) * t221 - mrSges(2,2) * t220 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t250;
t111 = m(2) * t221 + qJDD(1) * mrSges(2,1) - t244 * mrSges(2,2) + t113;
t110 = m(2) * t220 - t244 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t259;
t105 = mrSges(3,1) * g(1) + mrSges(3,3) * t192 + t229 * Ifges(3,5) + Ifges(3,6) * t231 - pkin(2) * t118 - t270;
t104 = -mrSges(3,2) * g(1) - mrSges(3,3) * t191 + Ifges(3,5) * t231 - t229 * Ifges(3,6) - pkin(7) * t118 - t237 * t107 + t241 * t109;
t103 = -mrSges(2,2) * g(1) - mrSges(2,3) * t221 + Ifges(2,5) * qJDD(1) - t244 * Ifges(2,6) - pkin(6) * t113 + t242 * t104 - t238 * t105;
t102 = Ifges(2,6) * qJDD(1) + t244 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t220 + t238 * t104 + t242 * t105 - pkin(1) * (-m(3) * g(1) + t118) + pkin(6) * t259;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t249, t103, t104, t109, t120, -t232 * t178 + t251; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t239 * t103 + t243 * t102 - pkin(5) * (-t243 * t110 + t239 * t111), t102, t105, t107, t119, -t204 * t176 + t253; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t243 * t103 + t239 * t102 + pkin(5) * (t239 * t110 + t243 * t111), t249, t250, t270, t271, -t203 * t180 - t252;];
m_new = t1;
