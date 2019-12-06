% Calculate vector of cutting torques with Newton-Euler for
% S5PRPPR1
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRPPR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:54
% EndTime: 2019-12-05 15:22:00
% DurationCPUTime: 4.59s
% Computational Cost: add. (47276->224), mult. (114767->317), div. (0->0), fcn. (74901->10), ass. (0->115)
t220 = sin(pkin(7));
t223 = cos(pkin(7));
t203 = t220 * g(1) - t223 * g(2);
t204 = -t223 * g(1) - t220 * g(2);
t225 = sin(qJ(2));
t227 = cos(qJ(2));
t181 = t225 * t203 + t227 * t204;
t228 = qJD(2) ^ 2;
t272 = -t228 * pkin(2) + qJDD(2) * qJ(3) + (2 * qJD(2) * qJD(3)) + t181;
t180 = t227 * t203 - t225 * t204;
t239 = -t228 * qJ(3) + qJDD(3) - t180;
t219 = sin(pkin(8));
t222 = cos(pkin(8));
t248 = -pkin(3) * t222 - qJ(4) * t219;
t264 = qJD(2) * t219;
t271 = (-pkin(2) + t248) * qJDD(2) + t239 - 0.2e1 * qJD(4) * t264;
t217 = -g(3) + qJDD(1);
t164 = t222 * t217 - t272 * t219;
t263 = t222 * qJD(2);
t165 = t219 * t217 + t272 * t222;
t195 = t248 * qJD(2);
t157 = t195 * t263 + t165;
t215 = t219 ^ 2;
t218 = sin(pkin(9));
t221 = cos(pkin(9));
t266 = t219 * t221;
t244 = -pkin(4) * t222 - pkin(6) * t266;
t265 = t271 * t221;
t149 = t244 * qJDD(2) + (-t157 + (-pkin(4) * t215 * t221 + pkin(6) * t219 * t222) * t228) * t218 + t265;
t152 = t221 * t157 + t271 * t218;
t191 = t244 * qJD(2);
t268 = t215 * t228;
t259 = t218 ^ 2 * t268;
t261 = qJDD(2) * t219;
t150 = -t218 * pkin(6) * t261 - pkin(4) * t259 + t191 * t263 + t152;
t224 = sin(qJ(5));
t226 = cos(qJ(5));
t148 = t224 * t149 + t226 * t150;
t156 = t195 * t264 + qJDD(4) - t164;
t153 = -pkin(6) * t259 + (pkin(4) * qJDD(2) * t218 + qJD(2) * t191 * t221) * t219 + t156;
t241 = (-t218 * t226 - t221 * t224) * t219;
t185 = qJD(2) * t241;
t240 = (-t218 * t224 + t221 * t226) * t219;
t186 = qJD(2) * t240;
t207 = qJD(5) - t263;
t158 = Ifges(6,5) * t186 + Ifges(6,6) * t185 + Ifges(6,3) * t207;
t160 = Ifges(6,1) * t186 + Ifges(6,4) * t185 + Ifges(6,5) * t207;
t171 = -t186 * qJD(5) + qJDD(2) * t241;
t172 = t185 * qJD(5) + qJDD(2) * t240;
t260 = t222 * qJDD(2);
t206 = qJDD(5) - t260;
t136 = -mrSges(6,1) * t153 + mrSges(6,3) * t148 + Ifges(6,4) * t172 + Ifges(6,2) * t171 + Ifges(6,6) * t206 - t186 * t158 + t207 * t160;
t147 = t226 * t149 - t224 * t150;
t159 = Ifges(6,4) * t186 + Ifges(6,2) * t185 + Ifges(6,6) * t207;
t137 = mrSges(6,2) * t153 - mrSges(6,3) * t147 + Ifges(6,1) * t172 + Ifges(6,4) * t171 + Ifges(6,5) * t206 + t185 * t158 - t207 * t159;
t249 = Ifges(5,5) * t221 - Ifges(5,6) * t218;
t182 = (-Ifges(5,3) * t222 + t249 * t219) * qJD(2);
t236 = -Ifges(5,5) * t222 + (Ifges(5,1) * t221 - Ifges(5,4) * t218) * t219;
t184 = t236 * qJD(2);
t235 = -Ifges(5,6) * t222 + (Ifges(5,4) * t221 - Ifges(5,2) * t218) * t219;
t177 = -t207 * mrSges(6,2) + t185 * mrSges(6,3);
t178 = t207 * mrSges(6,1) - t186 * mrSges(6,3);
t237 = m(6) * t153 - t171 * mrSges(6,1) + t172 * mrSges(6,2) - t185 * t177 + t186 * t178;
t166 = -t185 * mrSges(6,1) + t186 * mrSges(6,2);
t143 = m(6) * t147 + t206 * mrSges(6,1) - t172 * mrSges(6,3) - t186 * t166 + t207 * t177;
t144 = m(6) * t148 - t206 * mrSges(6,2) + t171 * mrSges(6,3) + t185 * t166 - t207 * t178;
t255 = -t224 * t143 + t226 * t144;
t122 = -mrSges(5,1) * t156 + mrSges(5,3) * t152 + t224 * t137 + t226 * t136 - pkin(4) * t237 + pkin(6) * t255 + (-t182 * t266 - t222 * t184) * qJD(2) + t235 * qJDD(2);
t135 = t226 * t143 + t224 * t144;
t151 = -t218 * t157 + t265;
t183 = t235 * qJD(2);
t267 = t218 * t219;
t123 = mrSges(5,2) * t156 - mrSges(5,3) * t151 - pkin(6) * t135 - t224 * t136 + t226 * t137 + (-t182 * t267 + t183 * t222) * qJD(2) + t236 * qJDD(2);
t252 = mrSges(5,1) * t218 + mrSges(5,2) * t221;
t187 = t252 * t264;
t242 = mrSges(5,2) * t222 - mrSges(5,3) * t267;
t189 = t242 * qJD(2);
t243 = -mrSges(5,1) * t222 - mrSges(5,3) * t266;
t133 = m(5) * t151 + t243 * qJDD(2) + (-t187 * t266 - t189 * t222) * qJD(2) + t135;
t190 = t243 * qJD(2);
t134 = m(5) * t152 + t242 * qJDD(2) + (-t187 * t267 + t190 * t222) * qJD(2) + t255;
t131 = -t218 * t133 + t221 * t134;
t232 = -m(5) * t156 - t237;
t246 = -t189 * t218 - t190 * t221;
t251 = Ifges(4,1) * t219 + Ifges(4,4) * t222;
t270 = -((Ifges(4,4) * t219 + Ifges(4,2) * t222) * t264 - t251 * t263) * qJD(2) - mrSges(4,1) * t164 + mrSges(4,2) * t165 - pkin(3) * ((t246 * qJD(2) - t252 * qJDD(2)) * t219 + t232) - qJ(4) * t131 - t221 * t122 - t218 * t123;
t269 = mrSges(4,2) * t219;
t196 = (-mrSges(4,1) * t222 + t269) * qJD(2);
t128 = m(4) * t165 + (qJDD(2) * mrSges(4,3) + qJD(2) * t196) * t222 + t131;
t139 = m(4) * t164 + ((-mrSges(4,3) - t252) * qJDD(2) + (-t196 + t246) * qJD(2)) * t219 + t232;
t256 = t222 * t128 - t219 * t139;
t119 = m(3) * t181 - t228 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t256;
t130 = t221 * t133 + t218 * t134;
t175 = -qJDD(2) * pkin(2) + t239;
t233 = -m(4) * t175 + mrSges(4,1) * t260 - t130 + (t222 ^ 2 * t228 + t268) * mrSges(4,3);
t125 = m(3) * t180 - t228 * mrSges(3,2) + (mrSges(3,1) - t269) * qJDD(2) + t233;
t114 = t225 * t119 + t227 * t125;
t121 = t219 * t128 + t222 * t139;
t257 = t227 * t119 - t225 * t125;
t250 = Ifges(4,5) * t219 + Ifges(4,6) * t222;
t247 = t183 * t221 + t184 * t218;
t197 = t250 * qJD(2);
t110 = mrSges(4,2) * t175 - mrSges(4,3) * t164 - qJ(4) * t130 + t251 * qJDD(2) - t218 * t122 + t221 * t123 + t197 * t263;
t234 = -mrSges(6,1) * t147 + mrSges(6,2) * t148 - Ifges(6,5) * t172 - Ifges(6,6) * t171 - Ifges(6,3) * t206 - t186 * t159 + t185 * t160;
t229 = mrSges(5,1) * t151 - mrSges(5,2) * t152 + pkin(4) * t135 - t234;
t116 = ((Ifges(4,4) - t249) * qJDD(2) + (-t197 - t247) * qJD(2)) * t219 + (Ifges(4,2) + Ifges(5,3)) * t260 + mrSges(4,3) * t165 - mrSges(4,1) * t175 - pkin(3) * t130 - t229;
t238 = -mrSges(3,2) * t181 + qJ(3) * t256 + t219 * t110 + t222 * t116 + pkin(2) * (-mrSges(4,2) * t261 + t233) + mrSges(3,1) * t180 + Ifges(3,3) * qJDD(2);
t231 = mrSges(2,1) * t203 - mrSges(2,2) * t204 + pkin(1) * t114 + t238;
t112 = m(2) * t204 + t257;
t111 = m(2) * t203 + t114;
t108 = -mrSges(3,1) * t217 + mrSges(3,3) * t181 + t228 * Ifges(3,5) - pkin(2) * t121 + (Ifges(3,6) - t250) * qJDD(2) + t270;
t107 = mrSges(3,2) * t217 - mrSges(3,3) * t180 + Ifges(3,5) * qJDD(2) - t228 * Ifges(3,6) - qJ(3) * t121 + t222 * t110 - t219 * t116;
t106 = mrSges(2,2) * t217 - mrSges(2,3) * t203 - pkin(5) * t114 + t227 * t107 - t225 * t108;
t105 = -mrSges(2,1) * t217 + mrSges(2,3) * t204 + t225 * t107 + t227 * t108 - pkin(1) * (m(3) * t217 + t121) + pkin(5) * t257;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t223 * t106 - t220 * t105 - qJ(1) * (t223 * t111 + t220 * t112), t106, t107, t110, t123, t137; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t220 * t106 + t223 * t105 + qJ(1) * (-t220 * t111 + t223 * t112), t105, t108, t116, t122, t136; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t231, t231, t238, t250 * qJDD(2) - t270, -Ifges(5,3) * t260 + (t247 * qJD(2) + t249 * qJDD(2)) * t219 + t229, -t234;];
m_new = t1;
