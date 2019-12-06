% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:35:40
% EndTime: 2019-12-05 17:35:45
% DurationCPUTime: 2.68s
% Computational Cost: add. (25524->241), mult. (53064->308), div. (0->0), fcn. (29382->8), ass. (0->103)
t230 = sin(qJ(1));
t232 = cos(qJ(1));
t208 = t232 * g(2) + t230 * g(3);
t200 = qJDD(1) * pkin(1) + t208;
t207 = t230 * g(2) - g(3) * t232;
t233 = qJD(1) ^ 2;
t201 = -pkin(1) * t233 + t207;
t226 = sin(pkin(7));
t228 = cos(pkin(7));
t165 = t226 * t200 + t228 * t201;
t274 = -pkin(2) * t233 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t165;
t224 = -g(1) + qJDD(2);
t225 = sin(pkin(8));
t227 = cos(pkin(8));
t151 = t224 * t227 - t274 * t225;
t266 = qJD(1) * t225;
t265 = qJD(1) * t227;
t152 = t225 * t224 + t274 * t227;
t250 = -pkin(3) * t227 - pkin(6) * t225;
t199 = t250 * qJD(1);
t150 = t199 * t265 + t152;
t164 = t228 * t200 - t226 * t201;
t243 = -t233 * qJ(3) + qJDD(3) - t164;
t155 = (-pkin(2) + t250) * qJDD(1) + t243;
t229 = sin(qJ(4));
t231 = cos(qJ(4));
t145 = t231 * t150 + t229 * t155;
t149 = t199 * t266 - t151;
t210 = qJD(4) - t265;
t171 = Ifges(5,5) * t210 + (Ifges(5,1) * t231 - Ifges(5,4) * t229) * t266;
t258 = t229 * t266;
t182 = -mrSges(6,2) * t210 - mrSges(6,3) * t258;
t257 = t231 * t266;
t185 = mrSges(6,1) * t210 - mrSges(6,3) * t257;
t187 = (mrSges(6,1) * t229 + mrSges(6,2) * t231) * t266;
t262 = qJD(1) * qJD(4);
t189 = (-qJDD(1) * t229 - t231 * t262) * t225;
t190 = (qJDD(1) * t231 - t229 * t262) * t225;
t261 = qJDD(1) * t227;
t209 = qJDD(4) - t261;
t184 = pkin(4) * t210 - qJ(5) * t257;
t251 = -0.2e1 * qJD(5) * t266;
t270 = t225 ^ 2 * t233;
t260 = t229 ^ 2 * t270;
t143 = -pkin(4) * t260 + qJ(5) * t189 - t184 * t210 + t229 * t251 + t145;
t147 = -pkin(4) * t189 - qJ(5) * t260 + t184 * t257 + qJDD(5) + t149;
t170 = Ifges(6,5) * t210 + (Ifges(6,1) * t231 - Ifges(6,4) * t229) * t266;
t242 = -mrSges(6,1) * t147 + mrSges(6,3) * t143 + Ifges(6,4) * t190 + Ifges(6,2) * t189 + Ifges(6,6) * t209 + t210 * t170;
t255 = -m(6) * t147 + t189 * mrSges(6,1);
t166 = Ifges(6,3) * t210 + (Ifges(6,5) * t231 - Ifges(6,6) * t229) * t266;
t268 = -t166 - Ifges(5,3) * t210 - (Ifges(5,5) * t231 - Ifges(5,6) * t229) * t266;
t269 = m(6) * t143 + t189 * mrSges(6,3);
t121 = Ifges(5,4) * t190 + Ifges(5,2) * t189 + Ifges(5,6) * t209 + t210 * t171 - mrSges(5,1) * t149 + mrSges(5,3) * t145 - pkin(4) * (t190 * mrSges(6,2) - t255) + qJ(5) * (-mrSges(6,2) * t209 - t185 * t210 + t269) + ((-pkin(4) * t182 - qJ(5) * t187) * t229 + (-pkin(4) * t185 + t268) * t231) * t266 + t242;
t154 = t231 * t155;
t141 = t231 * t251 + t209 * pkin(4) - t190 * qJ(5) + t154 + (-pkin(4) * t231 * t270 - qJ(5) * t210 * t266 - t150) * t229;
t259 = m(6) * t141 + t209 * mrSges(6,1) + t210 * t182;
t137 = -mrSges(6,3) * t190 - t187 * t257 + t259;
t144 = -t229 * t150 + t154;
t168 = Ifges(6,6) * t210 + (Ifges(6,4) * t231 - Ifges(6,2) * t229) * t266;
t169 = Ifges(5,6) * t210 + (Ifges(5,4) * t231 - Ifges(5,2) * t229) * t266;
t244 = mrSges(6,2) * t147 - mrSges(6,3) * t141 + Ifges(6,1) * t190 + Ifges(6,4) * t189 + Ifges(6,5) * t209;
t128 = mrSges(5,2) * t149 - mrSges(5,3) * t144 + Ifges(5,1) * t190 + Ifges(5,4) * t189 + Ifges(5,5) * t209 - qJ(5) * t137 + (-t168 - t169) * t210 + t268 * t258 + t244;
t183 = -mrSges(5,2) * t210 - mrSges(5,3) * t258;
t249 = (-t187 - (mrSges(5,1) * t229 + mrSges(5,2) * t231) * t266) * t266;
t132 = m(5) * t144 + mrSges(5,1) * t209 + t183 * t210 + (-mrSges(5,3) - mrSges(6,3)) * t190 + t231 * t249 + t259;
t267 = -mrSges(5,1) * t210 + mrSges(5,3) * t257 - t185;
t272 = -mrSges(5,2) - mrSges(6,2);
t133 = m(5) * t145 + mrSges(5,3) * t189 + t272 * t209 + t267 * t210 + t229 * t249 + t269;
t130 = -t229 * t132 + t231 * t133;
t238 = -m(5) * t149 + t189 * mrSges(5,1) + t272 * t190 + t255;
t241 = t267 * t231 + (-t182 - t183) * t229;
t248 = Ifges(4,1) * t225 + Ifges(4,4) * t227;
t273 = -((Ifges(4,4) * t225 + Ifges(4,2) * t227) * t266 - t248 * t265) * qJD(1) - mrSges(4,1) * t151 + mrSges(4,2) * t152 - pkin(3) * (t241 * t266 + t238) - pkin(6) * t130 - t121 * t231 - t128 * t229;
t271 = mrSges(4,2) * t225;
t195 = (-mrSges(4,1) * t227 + t271) * qJD(1);
t264 = qJDD(1) * mrSges(4,3);
t126 = m(4) * t152 + (qJD(1) * t195 + t264) * t227 + t130;
t135 = m(4) * t151 + (-t264 + (-t195 + t241) * qJD(1)) * t225 + t238;
t253 = t227 * t126 - t135 * t225;
t118 = m(3) * t165 - mrSges(3,1) * t233 - qJDD(1) * mrSges(3,2) + t253;
t129 = t231 * t132 + t229 * t133;
t158 = -qJDD(1) * pkin(2) + t243;
t237 = -m(4) * t158 + mrSges(4,1) * t261 - t129 + (t227 ^ 2 * t233 + t270) * mrSges(4,3);
t123 = m(3) * t164 - t233 * mrSges(3,2) + (mrSges(3,1) - t271) * qJDD(1) + t237;
t113 = t226 * t118 + t228 * t123;
t120 = t225 * t126 + t227 * t135;
t254 = t228 * t118 - t226 * t123;
t247 = Ifges(4,5) * t225 + Ifges(4,6) * t227;
t246 = t169 * t231 + t171 * t229;
t196 = t247 * qJD(1);
t109 = mrSges(4,2) * t158 - mrSges(4,3) * t151 - pkin(6) * t129 + t248 * qJDD(1) - t229 * t121 + t231 * t128 + t196 * t265;
t239 = -mrSges(6,1) * t141 + mrSges(6,2) * t143 - Ifges(6,5) * t190 - Ifges(6,6) * t189 - Ifges(6,3) * t209 - t168 * t257 - t170 * t258;
t234 = mrSges(5,1) * t144 - mrSges(5,2) * t145 + Ifges(5,5) * t190 + Ifges(5,6) * t189 + Ifges(5,3) * t209 + pkin(4) * t137 - t239;
t115 = -t234 + mrSges(4,3) * t152 - mrSges(4,1) * t158 - pkin(3) * t129 + Ifges(4,2) * t261 + (Ifges(4,4) * qJDD(1) + (-t196 - t246) * qJD(1)) * t225;
t240 = -mrSges(3,2) * t165 + qJ(3) * t253 + t225 * t109 + t227 * t115 + pkin(2) * (-qJDD(1) * t271 + t237) + mrSges(3,1) * t164 + Ifges(3,3) * qJDD(1);
t236 = mrSges(2,1) * t208 - mrSges(2,2) * t207 + Ifges(2,3) * qJDD(1) + pkin(1) * t113 + t240;
t111 = m(2) * t208 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t233 + t113;
t110 = m(2) * t207 - mrSges(2,1) * t233 - qJDD(1) * mrSges(2,2) + t254;
t107 = -mrSges(3,1) * t224 + mrSges(3,3) * t165 + Ifges(3,5) * t233 - pkin(2) * t120 + (Ifges(3,6) - t247) * qJDD(1) + t273;
t106 = mrSges(3,2) * t224 - mrSges(3,3) * t164 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t233 - qJ(3) * t120 + t109 * t227 - t115 * t225;
t105 = -mrSges(2,2) * g(1) - mrSges(2,3) * t208 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t233 - qJ(2) * t113 + t106 * t228 - t107 * t226;
t104 = Ifges(2,6) * qJDD(1) + t233 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t207 + t226 * t106 + t228 * t107 - pkin(1) * (m(3) * t224 + t120) + qJ(2) * t254;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t236, t105, t106, t109, t128, -t166 * t258 - t168 * t210 + t244; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t230 * t105 - t232 * t104 - pkin(5) * (t110 * t232 - t111 * t230), t104, t107, t115, t121, -t166 * t257 + t242; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232 * t105 - t230 * t104 + pkin(5) * (-t110 * t230 - t111 * t232), t236, t240, t247 * qJDD(1) - t273, t246 * t266 + t234, -t239;];
m_new = t1;
