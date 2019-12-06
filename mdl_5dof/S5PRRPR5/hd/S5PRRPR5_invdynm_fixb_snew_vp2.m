% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:25:37
% EndTime: 2019-12-05 16:25:56
% DurationCPUTime: 7.20s
% Computational Cost: add. (113338->266), mult. (239592->351), div. (0->0), fcn. (165026->12), ass. (0->117)
t234 = sin(pkin(9));
t237 = cos(pkin(9));
t224 = t234 * g(1) - t237 * g(2);
t225 = -t237 * g(1) - t234 * g(2);
t232 = -g(3) + qJDD(1);
t244 = cos(qJ(2));
t238 = cos(pkin(5));
t241 = sin(qJ(2));
t266 = t238 * t241;
t235 = sin(pkin(5));
t267 = t235 * t241;
t190 = t224 * t266 + t244 * t225 + t232 * t267;
t246 = qJD(2) ^ 2;
t185 = -t246 * pkin(2) + qJDD(2) * pkin(7) + t190;
t206 = -t235 * t224 + t238 * t232;
t240 = sin(qJ(3));
t243 = cos(qJ(3));
t170 = -t240 * t185 + t243 * t206;
t263 = qJD(2) * qJD(3);
t262 = t243 * t263;
t221 = t240 * qJDD(2) + t262;
t167 = (-t221 + t262) * qJ(4) + (t240 * t243 * t246 + qJDD(3)) * pkin(3) + t170;
t171 = t243 * t185 + t240 * t206;
t222 = t243 * qJDD(2) - t240 * t263;
t265 = qJD(2) * t240;
t226 = qJD(3) * pkin(3) - qJ(4) * t265;
t231 = t243 ^ 2;
t168 = -t231 * t246 * pkin(3) + t222 * qJ(4) - qJD(3) * t226 + t171;
t233 = sin(pkin(10));
t236 = cos(pkin(10));
t209 = (t233 * t240 - t236 * t243) * qJD(2);
t271 = 2 * qJD(4);
t163 = t233 * t167 + t236 * t168 - t209 * t271;
t210 = (t233 * t243 + t236 * t240) * qJD(2);
t192 = t209 * mrSges(5,1) + t210 * mrSges(5,2);
t198 = -t233 * t221 + t236 * t222;
t205 = qJD(3) * mrSges(5,1) - t210 * mrSges(5,3);
t193 = t209 * pkin(4) - t210 * pkin(8);
t245 = qJD(3) ^ 2;
t160 = -t245 * pkin(4) + qJDD(3) * pkin(8) - t209 * t193 + t163;
t189 = -t241 * t225 + (t224 * t238 + t232 * t235) * t244;
t252 = -qJDD(2) * pkin(2) - t189;
t169 = -t222 * pkin(3) + qJDD(4) + t226 * t265 + (-qJ(4) * t231 - pkin(7)) * t246 + t252;
t199 = t236 * t221 + t233 * t222;
t164 = (qJD(3) * t209 - t199) * pkin(8) + (qJD(3) * t210 - t198) * pkin(4) + t169;
t239 = sin(qJ(5));
t242 = cos(qJ(5));
t157 = -t239 * t160 + t242 * t164;
t200 = t242 * qJD(3) - t239 * t210;
t178 = t200 * qJD(5) + t239 * qJDD(3) + t242 * t199;
t201 = t239 * qJD(3) + t242 * t210;
t180 = -t200 * mrSges(6,1) + t201 * mrSges(6,2);
t208 = qJD(5) + t209;
t182 = -t208 * mrSges(6,2) + t200 * mrSges(6,3);
t197 = qJDD(5) - t198;
t153 = m(6) * t157 + t197 * mrSges(6,1) - t178 * mrSges(6,3) - t201 * t180 + t208 * t182;
t158 = t242 * t160 + t239 * t164;
t177 = -t201 * qJD(5) + t242 * qJDD(3) - t239 * t199;
t183 = t208 * mrSges(6,1) - t201 * mrSges(6,3);
t154 = m(6) * t158 - t197 * mrSges(6,2) + t177 * mrSges(6,3) + t200 * t180 - t208 * t183;
t259 = -t239 * t153 + t242 * t154;
t140 = m(5) * t163 - qJDD(3) * mrSges(5,2) + t198 * mrSges(5,3) - qJD(3) * t205 - t209 * t192 + t259;
t258 = -t236 * t167 + t233 * t168;
t162 = -0.2e1 * qJD(4) * t210 - t258;
t204 = -qJD(3) * mrSges(5,2) - t209 * mrSges(5,3);
t159 = -qJDD(3) * pkin(4) - t245 * pkin(8) + (t271 + t193) * t210 + t258;
t253 = -m(6) * t159 + t177 * mrSges(6,1) - t178 * mrSges(6,2) + t200 * t182 - t201 * t183;
t149 = m(5) * t162 + qJDD(3) * mrSges(5,1) - t199 * mrSges(5,3) + qJD(3) * t204 - t210 * t192 + t253;
t135 = t233 * t140 + t236 * t149;
t220 = (-mrSges(4,1) * t243 + mrSges(4,2) * t240) * qJD(2);
t264 = qJD(2) * t243;
t228 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t264;
t133 = m(4) * t170 + qJDD(3) * mrSges(4,1) - t221 * mrSges(4,3) + qJD(3) * t228 - t220 * t265 + t135;
t227 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t265;
t260 = t236 * t140 - t233 * t149;
t134 = m(4) * t171 - qJDD(3) * mrSges(4,2) + t222 * mrSges(4,3) - qJD(3) * t227 + t220 * t264 + t260;
t127 = t243 * t133 + t240 * t134;
t213 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t240 + Ifges(4,2) * t243) * qJD(2);
t214 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t240 + Ifges(4,4) * t243) * qJD(2);
t172 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t208;
t174 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t208;
t146 = -mrSges(6,1) * t159 + mrSges(6,3) * t158 + Ifges(6,4) * t178 + Ifges(6,2) * t177 + Ifges(6,6) * t197 - t201 * t172 + t208 * t174;
t173 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t208;
t147 = mrSges(6,2) * t159 - mrSges(6,3) * t157 + Ifges(6,1) * t178 + Ifges(6,4) * t177 + Ifges(6,5) * t197 + t200 * t172 - t208 * t173;
t187 = Ifges(5,4) * t210 - Ifges(5,2) * t209 + Ifges(5,6) * qJD(3);
t188 = Ifges(5,1) * t210 - Ifges(5,4) * t209 + Ifges(5,5) * qJD(3);
t250 = -mrSges(5,1) * t162 + mrSges(5,2) * t163 - Ifges(5,5) * t199 - Ifges(5,6) * t198 - Ifges(5,3) * qJDD(3) - pkin(4) * t253 - pkin(8) * t259 - t242 * t146 - t239 * t147 - t210 * t187 - t209 * t188;
t272 = mrSges(4,1) * t170 - mrSges(4,2) * t171 + Ifges(4,5) * t221 + Ifges(4,6) * t222 + Ifges(4,3) * qJDD(3) + pkin(3) * t135 + (t240 * t213 - t243 * t214) * qJD(2) - t250;
t113 = -mrSges(3,1) * t206 + mrSges(3,3) * t190 + t246 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t127 - t272;
t261 = -t240 * t133 + t243 * t134;
t125 = m(3) * t190 - t246 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t261;
t184 = -t246 * pkin(7) + t252;
t142 = t242 * t153 + t239 * t154;
t251 = m(5) * t169 - t198 * mrSges(5,1) + t199 * mrSges(5,2) + t209 * t204 + t210 * t205 + t142;
t248 = -m(4) * t184 + t222 * mrSges(4,1) - t221 * mrSges(4,2) - t227 * t265 + t228 * t264 - t251;
t137 = m(3) * t189 + qJDD(2) * mrSges(3,1) - t246 * mrSges(3,2) + t248;
t122 = t244 * t125 - t241 * t137;
t273 = pkin(6) * t122 + t113 * t244;
t268 = t137 * t244;
t126 = m(3) * t206 + t127;
t117 = t125 * t266 - t235 * t126 + t238 * t268;
t186 = Ifges(5,5) * t210 - Ifges(5,6) * t209 + Ifges(5,3) * qJD(3);
t128 = mrSges(5,2) * t169 - mrSges(5,3) * t162 + Ifges(5,1) * t199 + Ifges(5,4) * t198 + Ifges(5,5) * qJDD(3) - pkin(8) * t142 - qJD(3) * t187 - t239 * t146 + t242 * t147 - t209 * t186;
t249 = mrSges(6,1) * t157 - mrSges(6,2) * t158 + Ifges(6,5) * t178 + Ifges(6,6) * t177 + Ifges(6,3) * t197 + t201 * t173 - t200 * t174;
t129 = -mrSges(5,1) * t169 + mrSges(5,3) * t163 + Ifges(5,4) * t199 + Ifges(5,2) * t198 + Ifges(5,6) * qJDD(3) - pkin(4) * t142 + qJD(3) * t188 - t210 * t186 - t249;
t212 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t240 + Ifges(4,6) * t243) * qJD(2);
t118 = -mrSges(4,1) * t184 + mrSges(4,3) * t171 + Ifges(4,4) * t221 + Ifges(4,2) * t222 + Ifges(4,6) * qJDD(3) - pkin(3) * t251 + qJ(4) * t260 + qJD(3) * t214 + t233 * t128 + t236 * t129 - t212 * t265;
t119 = mrSges(4,2) * t184 - mrSges(4,3) * t170 + Ifges(4,1) * t221 + Ifges(4,4) * t222 + Ifges(4,5) * qJDD(3) - qJ(4) * t135 - qJD(3) * t213 + t236 * t128 - t233 * t129 + t212 * t264;
t109 = mrSges(3,1) * t189 - mrSges(3,2) * t190 + Ifges(3,3) * qJDD(2) + pkin(2) * t248 + pkin(7) * t261 + t243 * t118 + t240 * t119;
t111 = mrSges(3,2) * t206 - mrSges(3,3) * t189 + Ifges(3,5) * qJDD(2) - t246 * Ifges(3,6) - pkin(7) * t127 - t240 * t118 + t243 * t119;
t254 = mrSges(2,1) * t224 - mrSges(2,2) * t225 + pkin(1) * t117 + t238 * t109 + t111 * t267 + t273 * t235;
t120 = m(2) * t225 + t122;
t116 = t238 * t126 + (t125 * t241 + t268) * t235;
t114 = m(2) * t224 + t117;
t107 = mrSges(2,2) * t232 - mrSges(2,3) * t224 + t244 * t111 - t241 * t113 + (-t116 * t235 - t117 * t238) * pkin(6);
t106 = -mrSges(2,1) * t232 + mrSges(2,3) * t225 - pkin(1) * t116 - t235 * t109 + (t111 * t241 + t273) * t238;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t237 * t107 - t234 * t106 - qJ(1) * (t237 * t114 + t234 * t120), t107, t111, t119, t128, t147; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t234 * t107 + t237 * t106 + qJ(1) * (-t234 * t114 + t237 * t120), t106, t113, t118, t129, t146; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t254, t254, t109, t272, -t250, t249;];
m_new = t1;
