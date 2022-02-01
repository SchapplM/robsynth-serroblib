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
% m [6x1]
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
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:12:16
% EndTime: 2022-01-23 09:12:21
% DurationCPUTime: 2.70s
% Computational Cost: add. (25524->241), mult. (53064->309), div. (0->0), fcn. (29382->8), ass. (0->104)
t226 = sin(qJ(1));
t228 = cos(qJ(1));
t205 = t226 * g(1) - g(2) * t228;
t198 = qJDD(1) * pkin(1) + t205;
t206 = -g(1) * t228 - g(2) * t226;
t229 = qJD(1) ^ 2;
t199 = -pkin(1) * t229 + t206;
t222 = sin(pkin(7));
t224 = cos(pkin(7));
t163 = t222 * t198 + t224 * t199;
t271 = -pkin(2) * t229 + qJDD(1) * qJ(3) + (2 * qJD(1) * qJD(3)) + t163;
t220 = -g(3) + qJDD(2);
t221 = sin(pkin(8));
t223 = cos(pkin(8));
t149 = t220 * t223 - t271 * t221;
t263 = qJD(1) * t221;
t262 = qJD(1) * t223;
t150 = t221 * t220 + t271 * t223;
t246 = -pkin(3) * t223 - pkin(6) * t221;
t197 = t246 * qJD(1);
t148 = t197 * t262 + t150;
t162 = t198 * t224 - t222 * t199;
t239 = -qJ(3) * t229 + qJDD(3) - t162;
t153 = (-pkin(2) + t246) * qJDD(1) + t239;
t225 = sin(qJ(4));
t227 = cos(qJ(4));
t143 = t227 * t148 + t225 * t153;
t147 = t197 * t263 - t149;
t208 = qJD(4) - t262;
t169 = Ifges(5,5) * t208 + (Ifges(5,1) * t227 - Ifges(5,4) * t225) * t263;
t255 = t225 * t263;
t180 = -mrSges(6,2) * t208 - mrSges(6,3) * t255;
t254 = t227 * t263;
t183 = mrSges(6,1) * t208 - mrSges(6,3) * t254;
t185 = (mrSges(6,1) * t225 + mrSges(6,2) * t227) * t263;
t259 = qJD(1) * qJD(4);
t187 = (-qJDD(1) * t225 - t227 * t259) * t221;
t188 = (qJDD(1) * t227 - t225 * t259) * t221;
t258 = qJDD(1) * t223;
t207 = qJDD(4) - t258;
t182 = pkin(4) * t208 - qJ(5) * t254;
t247 = -0.2e1 * qJD(5) * t263;
t267 = t221 ^ 2 * t229;
t257 = t225 ^ 2 * t267;
t141 = -pkin(4) * t257 + qJ(5) * t187 - t182 * t208 + t225 * t247 + t143;
t145 = -pkin(4) * t187 - qJ(5) * t257 + t182 * t254 + qJDD(5) + t147;
t168 = Ifges(6,5) * t208 + (Ifges(6,1) * t227 - Ifges(6,4) * t225) * t263;
t238 = -mrSges(6,1) * t145 + mrSges(6,3) * t141 + Ifges(6,4) * t188 + Ifges(6,2) * t187 + Ifges(6,6) * t207 + t208 * t168;
t251 = -m(6) * t145 + t187 * mrSges(6,1);
t164 = Ifges(6,3) * t208 + (Ifges(6,5) * t227 - Ifges(6,6) * t225) * t263;
t265 = -t164 - Ifges(5,3) * t208 - (Ifges(5,5) * t227 - Ifges(5,6) * t225) * t263;
t266 = m(6) * t141 + t187 * mrSges(6,3);
t119 = Ifges(5,4) * t188 + Ifges(5,2) * t187 + Ifges(5,6) * t207 + t208 * t169 - mrSges(5,1) * t147 + mrSges(5,3) * t143 - pkin(4) * (t188 * mrSges(6,2) - t251) + qJ(5) * (-mrSges(6,2) * t207 - t183 * t208 + t266) + ((-pkin(4) * t180 - qJ(5) * t185) * t225 + (-pkin(4) * t183 + t265) * t227) * t263 + t238;
t152 = t227 * t153;
t139 = t227 * t247 + pkin(4) * t207 - qJ(5) * t188 + t152 + (-pkin(4) * t227 * t267 - qJ(5) * t208 * t263 - t148) * t225;
t256 = m(6) * t139 + t207 * mrSges(6,1) + t208 * t180;
t135 = -mrSges(6,3) * t188 - t185 * t254 + t256;
t142 = -t148 * t225 + t152;
t166 = Ifges(6,6) * t208 + (Ifges(6,4) * t227 - Ifges(6,2) * t225) * t263;
t167 = Ifges(5,6) * t208 + (Ifges(5,4) * t227 - Ifges(5,2) * t225) * t263;
t240 = mrSges(6,2) * t145 - mrSges(6,3) * t139 + Ifges(6,1) * t188 + Ifges(6,4) * t187 + Ifges(6,5) * t207;
t126 = mrSges(5,2) * t147 - mrSges(5,3) * t142 + Ifges(5,1) * t188 + Ifges(5,4) * t187 + Ifges(5,5) * t207 - qJ(5) * t135 + (-t166 - t167) * t208 + t265 * t255 + t240;
t181 = -mrSges(5,2) * t208 - mrSges(5,3) * t255;
t245 = (-t185 - (mrSges(5,1) * t225 + mrSges(5,2) * t227) * t263) * t263;
t130 = m(5) * t142 + mrSges(5,1) * t207 + t181 * t208 + (-mrSges(5,3) - mrSges(6,3)) * t188 + t227 * t245 + t256;
t264 = -mrSges(5,1) * t208 + mrSges(5,3) * t254 - t183;
t269 = -mrSges(5,2) - mrSges(6,2);
t131 = m(5) * t143 + mrSges(5,3) * t187 + t269 * t207 + t264 * t208 + t225 * t245 + t266;
t128 = -t130 * t225 + t227 * t131;
t234 = -m(5) * t147 + t187 * mrSges(5,1) + t269 * t188 + t251;
t237 = t264 * t227 + (-t180 - t181) * t225;
t244 = Ifges(4,1) * t221 + Ifges(4,4) * t223;
t270 = -((Ifges(4,4) * t221 + Ifges(4,2) * t223) * t263 - t244 * t262) * qJD(1) - mrSges(4,1) * t149 + mrSges(4,2) * t150 - pkin(3) * (t237 * t263 + t234) - pkin(6) * t128 - t119 * t227 - t126 * t225;
t268 = mrSges(4,2) * t221;
t193 = (-mrSges(4,1) * t223 + t268) * qJD(1);
t261 = qJDD(1) * mrSges(4,3);
t124 = m(4) * t150 + (qJD(1) * t193 + t261) * t223 + t128;
t133 = m(4) * t149 + (-t261 + (-t193 + t237) * qJD(1)) * t221 + t234;
t249 = t223 * t124 - t133 * t221;
t116 = m(3) * t163 - mrSges(3,1) * t229 - qJDD(1) * mrSges(3,2) + t249;
t127 = t130 * t227 + t131 * t225;
t156 = -qJDD(1) * pkin(2) + t239;
t233 = -m(4) * t156 + mrSges(4,1) * t258 - t127 + (t223 ^ 2 * t229 + t267) * mrSges(4,3);
t121 = m(3) * t162 - mrSges(3,2) * t229 + (mrSges(3,1) - t268) * qJDD(1) + t233;
t111 = t222 * t116 + t224 * t121;
t118 = t221 * t124 + t223 * t133;
t253 = t164 * t263;
t250 = t224 * t116 - t121 * t222;
t243 = Ifges(4,5) * t221 + Ifges(4,6) * t223;
t242 = t167 * t227 + t169 * t225;
t194 = t243 * qJD(1);
t107 = mrSges(4,2) * t156 - mrSges(4,3) * t149 - pkin(6) * t127 + t244 * qJDD(1) - t119 * t225 + t126 * t227 + t194 * t262;
t235 = -mrSges(6,1) * t139 + mrSges(6,2) * t141 - Ifges(6,5) * t188 - Ifges(6,6) * t187 - Ifges(6,3) * t207 - t166 * t254 - t168 * t255;
t230 = mrSges(5,1) * t142 - mrSges(5,2) * t143 + Ifges(5,5) * t188 + Ifges(5,6) * t187 + Ifges(5,3) * t207 + pkin(4) * t135 - t235;
t113 = -t230 + mrSges(4,3) * t150 - mrSges(4,1) * t156 - pkin(3) * t127 + Ifges(4,2) * t258 + (Ifges(4,4) * qJDD(1) + (-t194 - t242) * qJD(1)) * t221;
t236 = -mrSges(3,2) * t163 + qJ(3) * t249 + t221 * t107 + t223 * t113 + pkin(2) * (-qJDD(1) * t268 + t233) + mrSges(3,1) * t162 + Ifges(3,3) * qJDD(1);
t232 = mrSges(2,1) * t205 - mrSges(2,2) * t206 + Ifges(2,3) * qJDD(1) + pkin(1) * t111 + t236;
t109 = m(2) * t206 - mrSges(2,1) * t229 - qJDD(1) * mrSges(2,2) + t250;
t108 = m(2) * t205 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t229 + t111;
t105 = -mrSges(3,1) * t220 + mrSges(3,3) * t163 + Ifges(3,5) * t229 - pkin(2) * t118 + (Ifges(3,6) - t243) * qJDD(1) + t270;
t104 = mrSges(3,2) * t220 - mrSges(3,3) * t162 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t229 - qJ(3) * t118 + t107 * t223 - t113 * t221;
t103 = -mrSges(2,2) * g(3) - mrSges(2,3) * t205 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t229 - qJ(2) * t111 + t104 * t224 - t105 * t222;
t102 = Ifges(2,6) * qJDD(1) + t229 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t206 + t222 * t104 + t224 * t105 - pkin(1) * (m(3) * t220 + t118) + qJ(2) * t250;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t228 * t103 - t226 * t102 - pkin(5) * (t108 * t228 + t109 * t226), t103, t104, t107, t126, -t166 * t208 - t225 * t253 + t240; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t226 * t103 + t228 * t102 + pkin(5) * (-t108 * t226 + t109 * t228), t102, t105, t113, t119, -t227 * t253 + t238; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232, t232, t236, t243 * qJDD(1) - t270, t242 * t263 + t230, -t235;];
m_new = t1;
