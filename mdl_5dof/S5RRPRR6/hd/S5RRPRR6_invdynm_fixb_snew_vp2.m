% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:48
% EndTime: 2019-12-05 18:35:55
% DurationCPUTime: 5.53s
% Computational Cost: add. (92924->246), mult. (127422->323), div. (0->0), fcn. (78313->10), ass. (0->114)
t233 = sin(qJ(1));
t237 = cos(qJ(1));
t213 = t237 * g(2) + t233 * g(3);
t203 = qJDD(1) * pkin(1) + t213;
t212 = t233 * g(2) - t237 * g(3);
t238 = qJD(1) ^ 2;
t204 = -t238 * pkin(1) + t212;
t232 = sin(qJ(2));
t236 = cos(qJ(2));
t183 = t232 * t203 + t236 * t204;
t225 = (qJD(1) + qJD(2));
t222 = t225 ^ 2;
t223 = qJDD(1) + qJDD(2);
t271 = -t222 * pkin(2) + t223 * qJ(3) + (2 * qJD(3) * t225) + t183;
t228 = sin(pkin(9));
t266 = t225 * t228;
t229 = cos(pkin(9));
t264 = t229 * t225;
t169 = -t229 * g(1) - t271 * t228;
t170 = -t228 * g(1) + t271 * t229;
t253 = -pkin(3) * t229 - pkin(7) * t228;
t198 = t253 * t225;
t156 = t198 * t264 + t170;
t182 = t236 * t203 - t232 * t204;
t247 = -t222 * qJ(3) + qJDD(3) - t182;
t168 = (-pkin(2) + t253) * t223 + t247;
t235 = cos(qJ(4));
t167 = t235 * t168;
t231 = sin(qJ(4));
t261 = qJD(4) * t225;
t193 = (t223 * t235 - t231 * t261) * t228;
t265 = t229 * t223;
t208 = qJDD(4) - t265;
t209 = qJD(4) - t264;
t267 = t222 * t228 ^ 2;
t148 = t208 * pkin(4) - t193 * pkin(8) + t167 + (-pkin(4) * t235 * t267 - pkin(8) * t209 * t266 - t156) * t231;
t151 = t235 * t156 + t231 * t168;
t258 = t235 * t266;
t191 = t209 * pkin(4) - pkin(8) * t258;
t192 = (-t223 * t231 - t235 * t261) * t228;
t260 = t231 ^ 2 * t267;
t149 = -pkin(4) * t260 + t192 * pkin(8) - t209 * t191 + t151;
t230 = sin(qJ(5));
t234 = cos(qJ(5));
t147 = t230 * t148 + t234 * t149;
t155 = t198 * t266 - t169;
t152 = -t192 * pkin(4) - pkin(8) * t260 + t191 * t258 + t155;
t185 = (-t230 * t231 + t234 * t235) * t266;
t160 = -t185 * qJD(5) + t234 * t192 - t230 * t193;
t184 = (-t230 * t235 - t231 * t234) * t266;
t161 = t184 * qJD(5) + t230 * t192 + t234 * t193;
t207 = qJD(5) + t209;
t162 = Ifges(6,5) * t185 + Ifges(6,6) * t184 + Ifges(6,3) * t207;
t164 = Ifges(6,1) * t185 + Ifges(6,4) * t184 + Ifges(6,5) * t207;
t205 = qJDD(5) + t208;
t135 = -mrSges(6,1) * t152 + mrSges(6,3) * t147 + Ifges(6,4) * t161 + Ifges(6,2) * t160 + Ifges(6,6) * t205 - t185 * t162 + t207 * t164;
t146 = t234 * t148 - t230 * t149;
t163 = Ifges(6,4) * t185 + Ifges(6,2) * t184 + Ifges(6,6) * t207;
t136 = mrSges(6,2) * t152 - mrSges(6,3) * t146 + Ifges(6,1) * t161 + Ifges(6,4) * t160 + Ifges(6,5) * t205 + t184 * t162 - t207 * t163;
t178 = Ifges(5,3) * t209 + (Ifges(5,5) * t235 - Ifges(5,6) * t231) * t266;
t180 = Ifges(5,5) * t209 + (Ifges(5,1) * t235 - Ifges(5,4) * t231) * t266;
t176 = -t207 * mrSges(6,2) + t184 * mrSges(6,3);
t177 = t207 * mrSges(6,1) - t185 * mrSges(6,3);
t245 = m(6) * t152 - t160 * mrSges(6,1) + t161 * mrSges(6,2) - t184 * t176 + t185 * t177;
t171 = -t184 * mrSges(6,1) + t185 * mrSges(6,2);
t142 = m(6) * t146 + t205 * mrSges(6,1) - t161 * mrSges(6,3) - t185 * t171 + t207 * t176;
t143 = m(6) * t147 - t205 * mrSges(6,2) + t160 * mrSges(6,3) + t184 * t171 - t207 * t177;
t254 = -t230 * t142 + t234 * t143;
t119 = -mrSges(5,1) * t155 + mrSges(5,3) * t151 + Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * t208 - pkin(4) * t245 + pkin(8) * t254 + t234 * t135 + t230 * t136 - t178 * t258 + t209 * t180;
t134 = t234 * t142 + t230 * t143;
t150 = -t231 * t156 + t167;
t179 = Ifges(5,6) * t209 + (Ifges(5,4) * t235 - Ifges(5,2) * t231) * t266;
t259 = t231 * t266;
t122 = mrSges(5,2) * t155 - mrSges(5,3) * t150 + Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * t208 - pkin(8) * t134 - t230 * t135 + t234 * t136 - t178 * t259 - t209 * t179;
t188 = -t209 * mrSges(5,2) - mrSges(5,3) * t259;
t190 = (mrSges(5,1) * t231 + mrSges(5,2) * t235) * t266;
t132 = m(5) * t150 + t208 * mrSges(5,1) - t193 * mrSges(5,3) + t209 * t188 - t190 * t258 + t134;
t189 = t209 * mrSges(5,1) - mrSges(5,3) * t258;
t133 = m(5) * t151 - t208 * mrSges(5,2) + t192 * mrSges(5,3) - t209 * t189 - t190 * t259 + t254;
t130 = -t231 * t132 + t235 * t133;
t240 = -m(5) * t155 + t192 * mrSges(5,1) - t193 * mrSges(5,2) - t245;
t249 = -t188 * t231 - t189 * t235;
t252 = Ifges(4,1) * t228 + Ifges(4,4) * t229;
t270 = -((Ifges(4,4) * t228 + Ifges(4,2) * t229) * t266 - t252 * t264) * t225 - mrSges(4,1) * t169 + mrSges(4,2) * t170 - pkin(3) * (t249 * t266 + t240) - pkin(7) * t130 - t235 * t119 - t231 * t122;
t269 = mrSges(4,2) * t228;
t268 = mrSges(4,3) * t223;
t194 = (-mrSges(4,1) * t229 + t269) * t225;
t127 = m(4) * t170 + (t194 * t225 + t268) * t229 + t130;
t138 = m(4) * t169 + (-t268 + (-t194 + t249) * t225) * t228 + t240;
t255 = t229 * t127 - t228 * t138;
t118 = m(3) * t183 - t222 * mrSges(3,1) - t223 * mrSges(3,2) + t255;
t129 = t235 * t132 + t231 * t133;
t174 = -t223 * pkin(2) + t247;
t243 = -m(4) * t174 + mrSges(4,1) * t265 - t129 + (t222 * t229 ^ 2 + t267) * mrSges(4,3);
t124 = m(3) * t182 - t222 * mrSges(3,2) + (mrSges(3,1) - t269) * t223 + t243;
t113 = t232 * t118 + t236 * t124;
t121 = t228 * t127 + t229 * t138;
t256 = t236 * t118 - t232 * t124;
t251 = Ifges(4,5) * t228 + Ifges(4,6) * t229;
t250 = t179 * t235 + t180 * t231;
t195 = t251 * t225;
t109 = mrSges(4,2) * t174 - mrSges(4,3) * t169 - pkin(7) * t129 - t231 * t119 + t235 * t122 + t195 * t264 + t252 * t223;
t244 = -mrSges(6,1) * t146 + mrSges(6,2) * t147 - Ifges(6,5) * t161 - Ifges(6,6) * t160 - Ifges(6,3) * t205 - t185 * t163 + t184 * t164;
t239 = mrSges(5,1) * t150 - mrSges(5,2) * t151 + Ifges(5,5) * t193 + Ifges(5,6) * t192 + Ifges(5,3) * t208 + pkin(4) * t134 - t244;
t115 = -t239 + Ifges(4,2) * t265 + (Ifges(4,4) * t223 + (-t195 - t250) * t225) * t228 + mrSges(4,3) * t170 - mrSges(4,1) * t174 - pkin(3) * t129;
t246 = -mrSges(3,2) * t183 + qJ(3) * t255 + t228 * t109 + t229 * t115 + pkin(2) * (-t223 * t269 + t243) + mrSges(3,1) * t182 + Ifges(3,3) * t223;
t242 = mrSges(2,1) * t213 - mrSges(2,2) * t212 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t246;
t111 = m(2) * t213 + qJDD(1) * mrSges(2,1) - t238 * mrSges(2,2) + t113;
t110 = m(2) * t212 - t238 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t256;
t107 = mrSges(3,1) * g(1) + mrSges(3,3) * t183 + t222 * Ifges(3,5) - pkin(2) * t121 + (Ifges(3,6) - t251) * t223 + t270;
t106 = -mrSges(3,2) * g(1) - mrSges(3,3) * t182 + Ifges(3,5) * t223 - t222 * Ifges(3,6) - qJ(3) * t121 + t229 * t109 - t228 * t115;
t105 = -mrSges(2,2) * g(1) - mrSges(2,3) * t213 + Ifges(2,5) * qJDD(1) - t238 * Ifges(2,6) - pkin(6) * t113 + t236 * t106 - t232 * t107;
t104 = Ifges(2,6) * qJDD(1) + t238 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t212 + t232 * t106 + t236 * t107 - pkin(1) * (-m(3) * g(1) + t121) + pkin(6) * t256;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t242, t105, t106, t109, t122, t136; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t233 * t105 - t237 * t104 - pkin(5) * (t237 * t110 - t233 * t111), t104, t107, t115, t119, t135; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t237 * t105 - t233 * t104 + pkin(5) * (-t233 * t110 - t237 * t111), t242, t246, t251 * t223 - t270, t250 * t266 + t239, -t244;];
m_new = t1;
