% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR2
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
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
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:39
% EndTime: 2019-12-05 17:39:47
% DurationCPUTime: 4.44s
% Computational Cost: add. (50285->246), mult. (113186->300), div. (0->0), fcn. (76079->8), ass. (0->109)
t236 = qJD(1) ^ 2;
t228 = sin(pkin(8));
t218 = t228 ^ 2;
t229 = cos(pkin(8));
t267 = t229 ^ 2 + t218;
t259 = t267 * mrSges(4,3);
t275 = t236 * t259;
t232 = sin(qJ(1));
t235 = cos(qJ(1));
t205 = t232 * g(1) - t235 * g(2);
t250 = -t236 * qJ(2) + qJDD(2) - t205;
t269 = -pkin(1) - qJ(3);
t274 = -(2 * qJD(1) * qJD(3)) + t269 * qJDD(1) + t250;
t206 = -t235 * g(1) - t232 * g(2);
t273 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t206;
t272 = pkin(3) * t236;
t271 = mrSges(2,1) - mrSges(3,2);
t270 = -Ifges(2,6) + Ifges(3,5);
t268 = Ifges(4,6) * t228;
t180 = t228 * g(3) + t274 * t229;
t163 = (-pkin(6) * qJDD(1) - t228 * t272) * t229 + t180;
t181 = -g(3) * t229 + t274 * t228;
t263 = qJDD(1) * t228;
t164 = -pkin(6) * t263 - t218 * t272 + t181;
t231 = sin(qJ(4));
t234 = cos(qJ(4));
t152 = t234 * t163 - t231 * t164;
t253 = -t228 * t231 + t229 * t234;
t254 = -t228 * t234 - t229 * t231;
t198 = t254 * qJD(1);
t265 = t198 * qJD(4);
t183 = t253 * qJDD(1) + t265;
t199 = t253 * qJD(1);
t141 = (-t183 + t265) * pkin(7) + (t198 * t199 + qJDD(4)) * pkin(4) + t152;
t153 = t231 * t163 + t234 * t164;
t182 = -t199 * qJD(4) + t254 * qJDD(1);
t190 = qJD(4) * pkin(4) - pkin(7) * t199;
t197 = t198 ^ 2;
t142 = -pkin(4) * t197 + pkin(7) * t182 - qJD(4) * t190 + t153;
t230 = sin(qJ(5));
t233 = cos(qJ(5));
t139 = t141 * t233 - t142 * t230;
t172 = t198 * t233 - t199 * t230;
t150 = qJD(5) * t172 + t182 * t230 + t183 * t233;
t173 = t198 * t230 + t199 * t233;
t158 = -mrSges(6,1) * t172 + mrSges(6,2) * t173;
t220 = qJD(4) + qJD(5);
t165 = -mrSges(6,2) * t220 + mrSges(6,3) * t172;
t217 = qJDD(4) + qJDD(5);
t136 = m(6) * t139 + mrSges(6,1) * t217 - mrSges(6,3) * t150 - t158 * t173 + t165 * t220;
t140 = t141 * t230 + t142 * t233;
t149 = -qJD(5) * t173 + t182 * t233 - t183 * t230;
t166 = mrSges(6,1) * t220 - mrSges(6,3) * t173;
t137 = m(6) * t140 - mrSges(6,2) * t217 + mrSges(6,3) * t149 + t158 * t172 - t166 * t220;
t127 = t233 * t136 + t230 * t137;
t176 = -mrSges(5,1) * t198 + mrSges(5,2) * t199;
t188 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t198;
t124 = m(5) * t152 + qJDD(4) * mrSges(5,1) - mrSges(5,3) * t183 + qJD(4) * t188 - t176 * t199 + t127;
t189 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t199;
t257 = -t136 * t230 + t233 * t137;
t125 = m(5) * t153 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t182 - qJD(4) * t189 + t176 * t198 + t257;
t120 = t234 * t124 + t231 * t125;
t266 = t236 * (Ifges(4,5) * t229 - t268);
t262 = qJDD(1) * t229;
t260 = Ifges(3,4) + t268;
t252 = -qJDD(1) * mrSges(4,3) - t236 * (mrSges(4,1) * t228 + mrSges(4,2) * t229);
t117 = m(4) * t180 + t252 * t229 + t120;
t258 = -t231 * t124 + t234 * t125;
t118 = m(4) * t181 + t252 * t228 + t258;
t113 = -t117 * t228 + t229 * t118;
t256 = Ifges(4,1) * t229 - Ifges(4,4) * t228;
t255 = Ifges(4,4) * t229 - Ifges(4,2) * t228;
t112 = t229 * t117 + t228 * t118;
t249 = qJDD(3) + t273;
t168 = pkin(3) * t263 + (-t267 * pkin(6) + t269) * t236 + t249;
t144 = -t182 * pkin(4) - t197 * pkin(7) + t199 * t190 + t168;
t248 = -m(6) * t144 + t149 * mrSges(6,1) - t150 * mrSges(6,2) + t172 * t165 - t173 * t166;
t196 = -qJDD(1) * pkin(1) + t250;
t247 = -m(3) * t196 + t236 * mrSges(3,3) - t112;
t155 = Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t220;
t156 = Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t220;
t246 = -mrSges(6,1) * t139 + mrSges(6,2) * t140 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * t217 - t173 * t155 + t172 * t156;
t154 = Ifges(6,5) * t173 + Ifges(6,6) * t172 + Ifges(6,3) * t220;
t128 = -mrSges(6,1) * t144 + mrSges(6,3) * t140 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t217 - t154 * t173 + t156 * t220;
t129 = mrSges(6,2) * t144 - mrSges(6,3) * t139 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t217 + t154 * t172 - t155 * t220;
t169 = Ifges(5,5) * t199 + Ifges(5,6) * t198 + Ifges(5,3) * qJD(4);
t171 = Ifges(5,1) * t199 + Ifges(5,4) * t198 + Ifges(5,5) * qJD(4);
t114 = -mrSges(5,1) * t168 + mrSges(5,3) * t153 + Ifges(5,4) * t183 + Ifges(5,2) * t182 + Ifges(5,6) * qJDD(4) + pkin(4) * t248 + pkin(7) * t257 + qJD(4) * t171 + t233 * t128 + t230 * t129 - t199 * t169;
t170 = Ifges(5,4) * t199 + Ifges(5,2) * t198 + Ifges(5,6) * qJD(4);
t115 = mrSges(5,2) * t168 - mrSges(5,3) * t152 + Ifges(5,1) * t183 + Ifges(5,4) * t182 + Ifges(5,5) * qJDD(4) - pkin(7) * t127 - qJD(4) * t170 - t128 * t230 + t129 * t233 + t169 * t198;
t187 = t269 * t236 + t249;
t243 = -m(5) * t168 + t182 * mrSges(5,1) - t183 * mrSges(5,2) + t198 * t188 - t199 * t189 + t248;
t106 = -mrSges(4,1) * t187 + mrSges(4,3) * t181 + pkin(3) * t243 + pkin(6) * t258 + t255 * qJDD(1) + t234 * t114 + t231 * t115 - t229 * t266;
t108 = mrSges(4,2) * t187 - mrSges(4,3) * t180 - pkin(6) * t120 + t256 * qJDD(1) - t231 * t114 + t234 * t115 - t228 * t266;
t192 = t236 * pkin(1) - t273;
t245 = mrSges(3,2) * t196 - mrSges(3,3) * t192 + Ifges(3,1) * qJDD(1) - qJ(3) * t112 - t106 * t228 + t229 * t108;
t241 = -m(4) * t187 - mrSges(4,1) * t263 - mrSges(4,2) * t262 + t243;
t244 = -mrSges(3,1) * t192 - pkin(2) * (t241 + t275) - qJ(3) * t113 - t229 * t106 - t228 * t108;
t239 = -m(3) * t192 + t236 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t241;
t242 = -mrSges(2,2) * t206 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t247) + qJ(2) * (t239 - t275) + mrSges(2,1) * t205 + Ifges(2,3) * qJDD(1) + t245;
t240 = -mrSges(5,1) * t152 + mrSges(5,2) * t153 - Ifges(5,5) * t183 - Ifges(5,6) * t182 - Ifges(5,3) * qJDD(4) - pkin(4) * t127 - t199 * t170 + t198 * t171 + t246;
t238 = -mrSges(4,1) * t180 + mrSges(4,2) * t181 - Ifges(4,5) * t262 - pkin(3) * t120 + t240 + (-t228 * t256 - t229 * t255) * t236;
t237 = -mrSges(3,1) * t196 - pkin(2) * t112 + t238;
t130 = t239 + m(2) * t206 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t259) * t236;
t111 = -m(3) * g(3) + t113;
t109 = m(2) * t205 - t236 * mrSges(2,2) + t271 * qJDD(1) + t247;
t105 = -t237 - mrSges(2,3) * t205 - qJ(2) * t111 + (-mrSges(2,2) + mrSges(3,3)) * g(3) + (Ifges(2,5) - t260) * qJDD(1) + t270 * t236;
t104 = mrSges(2,3) * t206 - pkin(1) * t111 + (-Ifges(3,4) + Ifges(2,5)) * t236 - t270 * qJDD(1) + t271 * g(3) + t244;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t235 * t105 - t232 * t104 - pkin(5) * (t109 * t235 + t130 * t232), t105, t245, t108, t115, t129; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t232 * t105 + t235 * t104 + pkin(5) * (-t109 * t232 + t130 * t235), t104, -mrSges(3,3) * g(3) - t236 * Ifges(3,5) + t260 * qJDD(1) + t237, t106, t114, t128; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t242, t242, mrSges(3,2) * g(3) + t236 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t244, -Ifges(4,6) * t263 - t238, -t240, -t246;];
m_new = t1;
