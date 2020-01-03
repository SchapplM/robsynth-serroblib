% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta3]';
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
% Datum: 2019-12-31 17:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:55:00
% EndTime: 2019-12-31 17:55:03
% DurationCPUTime: 1.91s
% Computational Cost: add. (18596->242), mult. (40840->281), div. (0->0), fcn. (23954->6), ass. (0->99)
t227 = qJD(1) ^ 2;
t220 = sin(pkin(7));
t207 = t220 ^ 2;
t221 = cos(pkin(7));
t262 = t221 ^ 2 + t207;
t252 = t262 * mrSges(4,3);
t273 = t227 * t252;
t223 = sin(qJ(1));
t225 = cos(qJ(1));
t196 = g(1) * t223 - t225 * g(2);
t242 = -qJ(2) * t227 + qJDD(2) - t196;
t266 = -pkin(1) - qJ(3);
t272 = -(2 * qJD(1) * qJD(3)) + t266 * qJDD(1) + t242;
t197 = -g(1) * t225 - g(2) * t223;
t271 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t197;
t270 = pkin(3) * t227;
t269 = mrSges(2,1) - mrSges(3,2);
t268 = -mrSges(5,3) - mrSges(6,2);
t267 = -Ifges(2,6) + Ifges(3,5);
t265 = Ifges(4,6) * t220;
t169 = t220 * g(3) + t272 * t221;
t144 = (-pkin(6) * qJDD(1) - t220 * t270) * t221 + t169;
t170 = -g(3) * t221 + t272 * t220;
t257 = qJDD(1) * t220;
t145 = -pkin(6) * t257 - t207 * t270 + t170;
t222 = sin(qJ(4));
t224 = cos(qJ(4));
t139 = t222 * t144 + t224 * t145;
t246 = t220 * t224 + t221 * t222;
t245 = -t220 * t222 + t221 * t224;
t190 = t245 * qJD(1);
t259 = qJD(4) * t190;
t171 = t246 * qJDD(1) + t259;
t180 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t190;
t189 = t246 * qJD(1);
t158 = pkin(4) * t189 - qJ(5) * t190;
t226 = qJD(4) ^ 2;
t134 = -pkin(4) * t226 + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) - t158 * t189 + t139;
t181 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t190;
t254 = m(6) * t134 + qJDD(4) * mrSges(6,3) + qJD(4) * t181;
t159 = mrSges(6,1) * t189 - mrSges(6,3) * t190;
t263 = -mrSges(5,1) * t189 - mrSges(5,2) * t190 - t159;
t123 = m(5) * t139 - qJDD(4) * mrSges(5,2) - qJD(4) * t180 + t268 * t171 + t263 * t189 + t254;
t138 = t144 * t224 - t145 * t222;
t260 = qJD(4) * t189;
t172 = t245 * qJDD(1) - t260;
t179 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t189;
t136 = -qJDD(4) * pkin(4) - qJ(5) * t226 + t158 * t190 + qJDD(5) - t138;
t182 = -mrSges(6,2) * t189 + qJD(4) * mrSges(6,3);
t250 = -m(6) * t136 + qJDD(4) * mrSges(6,1) + qJD(4) * t182;
t124 = m(5) * t138 + qJDD(4) * mrSges(5,1) + qJD(4) * t179 + t268 * t172 + t263 * t190 + t250;
t116 = t222 * t123 + t224 * t124;
t151 = Ifges(6,4) * t190 + Ifges(6,2) * qJD(4) + Ifges(6,6) * t189;
t264 = -Ifges(5,5) * t190 + Ifges(5,6) * t189 - Ifges(5,3) * qJD(4) - t151;
t261 = t227 * (Ifges(4,5) * t221 - t265);
t256 = qJDD(1) * t221;
t253 = Ifges(3,4) + t265;
t244 = -qJDD(1) * mrSges(4,3) - t227 * (mrSges(4,1) * t220 + mrSges(4,2) * t221);
t113 = m(4) * t169 + t244 * t221 + t116;
t251 = t224 * t123 - t222 * t124;
t114 = m(4) * t170 + t244 * t220 + t251;
t109 = -t113 * t220 + t221 * t114;
t240 = qJDD(3) + t271;
t147 = pkin(3) * t257 + (-t262 * pkin(6) + t266) * t227 + t240;
t131 = -0.2e1 * qJD(5) * t190 + (-t172 + t260) * qJ(5) + (t171 + t259) * pkin(4) + t147;
t249 = -mrSges(6,1) * t131 + mrSges(6,2) * t134;
t248 = Ifges(4,1) * t221 - Ifges(4,4) * t220;
t247 = Ifges(4,4) * t221 - Ifges(4,2) * t220;
t108 = t113 * t221 + t114 * t220;
t149 = Ifges(6,5) * t190 + Ifges(6,6) * qJD(4) + Ifges(6,3) * t189;
t241 = mrSges(6,2) * t136 - mrSges(6,3) * t131 + Ifges(6,1) * t172 + Ifges(6,4) * qJDD(4) + Ifges(6,5) * t171 + qJD(4) * t149;
t239 = -m(6) * t131 - t171 * mrSges(6,1) + mrSges(6,3) * t172 + t181 * t190 - t189 * t182;
t188 = -qJDD(1) * pkin(1) + t242;
t238 = -m(3) * t188 + t227 * mrSges(3,3) - t108;
t153 = Ifges(6,1) * t190 + Ifges(6,4) * qJD(4) + Ifges(6,5) * t189;
t237 = mrSges(6,1) * t136 - mrSges(6,3) * t134 - Ifges(6,4) * t172 - Ifges(6,2) * qJDD(4) - Ifges(6,6) * t171 + t190 * t149 - t189 * t153;
t154 = Ifges(5,1) * t190 - Ifges(5,4) * t189 + Ifges(5,5) * qJD(4);
t110 = -mrSges(5,1) * t147 + mrSges(5,3) * t139 + pkin(4) * t239 + t264 * t190 + (Ifges(5,4) - Ifges(6,5)) * t172 + (-Ifges(5,2) - Ifges(6,3)) * t171 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + (t153 + t154) * qJD(4) + t249;
t152 = Ifges(5,4) * t190 - Ifges(5,2) * t189 + Ifges(5,6) * qJD(4);
t111 = mrSges(5,2) * t147 - mrSges(5,3) * t138 + Ifges(5,1) * t172 - Ifges(5,4) * t171 + Ifges(5,5) * qJDD(4) + qJ(5) * t239 - qJD(4) * t152 + t264 * t189 + t241;
t176 = t266 * t227 + t240;
t234 = m(5) * t147 + mrSges(5,1) * t171 + t172 * mrSges(5,2) + t179 * t189 + t190 * t180 - t239;
t102 = -mrSges(4,1) * t176 + mrSges(4,3) * t170 - pkin(3) * t234 + pkin(6) * t251 + t247 * qJDD(1) + t224 * t110 + t222 * t111 - t221 * t261;
t104 = mrSges(4,2) * t176 - mrSges(4,3) * t169 - pkin(6) * t116 + t248 * qJDD(1) - t110 * t222 + t111 * t224 - t220 * t261;
t184 = pkin(1) * t227 - t271;
t236 = mrSges(3,2) * t188 - mrSges(3,3) * t184 + Ifges(3,1) * qJDD(1) - qJ(3) * t108 - t102 * t220 + t221 * t104;
t232 = -m(4) * t176 - mrSges(4,1) * t257 - mrSges(4,2) * t256 - t234;
t235 = -mrSges(3,1) * t184 - pkin(2) * (t232 + t273) - qJ(3) * t109 - t102 * t221 - t104 * t220;
t230 = -m(3) * t184 + t227 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t232;
t233 = -mrSges(2,2) * t197 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t238) + qJ(2) * (t230 - t273) + mrSges(2,1) * t196 + Ifges(2,3) * qJDD(1) + t236;
t231 = mrSges(5,2) * t139 - t189 * t154 - qJ(5) * (-mrSges(6,2) * t171 - t159 * t189 + t254) - pkin(4) * (-mrSges(6,2) * t172 - t159 * t190 + t250) - mrSges(5,1) * t138 - t190 * t152 + Ifges(5,6) * t171 - Ifges(5,5) * t172 - Ifges(5,3) * qJDD(4) + t237;
t229 = -mrSges(4,1) * t169 + mrSges(4,2) * t170 - Ifges(4,5) * t256 - pkin(3) * t116 + t231 + (-t220 * t248 - t221 * t247) * t227;
t228 = -mrSges(3,1) * t188 - pkin(2) * t108 + t229;
t117 = t230 - qJDD(1) * mrSges(2,2) + (-mrSges(2,1) - t252) * t227 + m(2) * t197;
t107 = -m(3) * g(3) + t109;
t105 = m(2) * t196 - mrSges(2,2) * t227 + t269 * qJDD(1) + t238;
t101 = t267 * t227 + (Ifges(2,5) - t253) * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t196 - qJ(2) * t107 - t228;
t100 = mrSges(2,3) * t197 - pkin(1) * t107 + (-Ifges(3,4) + Ifges(2,5)) * t227 - t267 * qJDD(1) + t269 * g(3) + t235;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t225 * t101 - t223 * t100 - pkin(5) * (t105 * t225 + t117 * t223), t101, t236, t104, t111, -t151 * t189 + t241; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t223 * t101 + t225 * t100 + pkin(5) * (-t105 * t223 + t117 * t225), t100, -mrSges(3,3) * g(3) - t227 * Ifges(3,5) + t253 * qJDD(1) + t228, t102, t110, -t237; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t233, t233, mrSges(3,2) * g(3) + Ifges(3,4) * t227 + Ifges(3,5) * qJDD(1) - t235, -Ifges(4,6) * t257 - t229, -t231, Ifges(6,5) * t172 + Ifges(6,6) * qJDD(4) + Ifges(6,3) * t171 - qJD(4) * t153 + t151 * t190 - t249;];
m_new = t1;
