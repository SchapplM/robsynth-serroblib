% Calculate vector of cutting torques with Newton-Euler for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPRR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:41
% EndTime: 2020-01-03 11:27:57
% DurationCPUTime: 7.09s
% Computational Cost: add. (82225->247), mult. (182364->311), div. (0->0), fcn. (122991->10), ass. (0->110)
t238 = qJD(1) ^ 2;
t234 = sin(qJ(1));
t237 = cos(qJ(1));
t210 = -t237 * g(2) - t234 * g(3);
t206 = qJDD(1) * pkin(1) + t210;
t209 = -t234 * g(2) + t237 * g(3);
t207 = -t238 * pkin(1) + t209;
t229 = sin(pkin(8));
t231 = cos(pkin(8));
t192 = t229 * t206 + t231 * t207;
t182 = -t238 * pkin(2) + qJDD(1) * qJ(3) + t192;
t228 = sin(pkin(9));
t227 = -g(1) + qJDD(2);
t230 = cos(pkin(9));
t260 = qJD(1) * qJD(3);
t263 = t230 * t227 - 0.2e1 * t228 * t260;
t266 = pkin(3) * t230;
t167 = (-pkin(6) * qJDD(1) + t238 * t266 - t182) * t228 + t263;
t171 = t228 * t227 + (t182 + 0.2e1 * t260) * t230;
t259 = qJDD(1) * t230;
t222 = t230 ^ 2;
t264 = t222 * t238;
t168 = -pkin(3) * t264 + pkin(6) * t259 + t171;
t233 = sin(qJ(4));
t236 = cos(qJ(4));
t149 = t236 * t167 - t233 * t168;
t248 = t228 * t236 + t230 * t233;
t247 = -t228 * t233 + t230 * t236;
t197 = t247 * qJD(1);
t261 = t197 * qJD(4);
t189 = t248 * qJDD(1) + t261;
t198 = t248 * qJD(1);
t144 = (-t189 + t261) * pkin(7) + (t197 * t198 + qJDD(4)) * pkin(4) + t149;
t150 = t233 * t167 + t236 * t168;
t188 = -t198 * qJD(4) + t247 * qJDD(1);
t195 = qJD(4) * pkin(4) - t198 * pkin(7);
t196 = t197 ^ 2;
t145 = -t196 * pkin(4) + t188 * pkin(7) - qJD(4) * t195 + t150;
t232 = sin(qJ(5));
t235 = cos(qJ(5));
t142 = t235 * t144 - t232 * t145;
t180 = t235 * t197 - t232 * t198;
t157 = t180 * qJD(5) + t232 * t188 + t235 * t189;
t181 = t232 * t197 + t235 * t198;
t163 = -t180 * mrSges(6,1) + t181 * mrSges(6,2);
t223 = qJD(4) + qJD(5);
t173 = -t223 * mrSges(6,2) + t180 * mrSges(6,3);
t220 = qJDD(4) + qJDD(5);
t139 = m(6) * t142 + t220 * mrSges(6,1) - t157 * mrSges(6,3) - t181 * t163 + t223 * t173;
t143 = t232 * t144 + t235 * t145;
t156 = -t181 * qJD(5) + t235 * t188 - t232 * t189;
t174 = t223 * mrSges(6,1) - t181 * mrSges(6,3);
t140 = m(6) * t143 - t220 * mrSges(6,2) + t156 * mrSges(6,3) + t180 * t163 - t223 * t174;
t130 = t235 * t139 + t232 * t140;
t184 = -t197 * mrSges(5,1) + t198 * mrSges(5,2);
t193 = -qJD(4) * mrSges(5,2) + t197 * mrSges(5,3);
t127 = m(5) * t149 + qJDD(4) * mrSges(5,1) - t189 * mrSges(5,3) + qJD(4) * t193 - t198 * t184 + t130;
t194 = qJD(4) * mrSges(5,1) - t198 * mrSges(5,3);
t255 = -t232 * t139 + t235 * t140;
t128 = m(5) * t150 - qJDD(4) * mrSges(5,2) + t188 * mrSges(5,3) - qJD(4) * t194 + t197 * t184 + t255;
t123 = t236 * t127 + t233 * t128;
t170 = -t228 * t182 + t263;
t178 = Ifges(5,4) * t198 + Ifges(5,2) * t197 + Ifges(5,6) * qJD(4);
t179 = Ifges(5,1) * t198 + Ifges(5,4) * t197 + Ifges(5,5) * qJD(4);
t159 = Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t223;
t160 = Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t223;
t244 = -mrSges(6,1) * t142 + mrSges(6,2) * t143 - Ifges(6,5) * t157 - Ifges(6,6) * t156 - Ifges(6,3) * t220 - t181 * t159 + t180 * t160;
t240 = -mrSges(5,1) * t149 + mrSges(5,2) * t150 - Ifges(5,5) * t189 - Ifges(5,6) * t188 - Ifges(5,3) * qJDD(4) - pkin(4) * t130 - t198 * t178 + t197 * t179 + t244;
t253 = Ifges(4,4) * t228 + Ifges(4,2) * t230;
t254 = Ifges(4,1) * t228 + Ifges(4,4) * t230;
t267 = -mrSges(4,1) * t170 + mrSges(4,2) * t171 - pkin(3) * t123 - (t228 * t253 - t230 * t254) * t238 + t240;
t265 = mrSges(4,2) * t228;
t246 = mrSges(4,3) * qJDD(1) + t238 * (-mrSges(4,1) * t230 + t265);
t121 = m(4) * t170 - t246 * t228 + t123;
t256 = -t233 * t127 + t236 * t128;
t122 = m(4) * t171 + t246 * t230 + t256;
t257 = -t228 * t121 + t230 * t122;
t113 = m(3) * t192 - t238 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t257;
t191 = t231 * t206 - t229 * t207;
t251 = qJDD(3) - t191;
t176 = -qJDD(1) * pkin(2) - t238 * qJ(3) + t251;
t221 = t228 ^ 2;
t169 = (-pkin(2) - t266) * qJDD(1) + (-qJ(3) + (-t221 - t222) * pkin(6)) * t238 + t251;
t147 = -t188 * pkin(4) - t196 * pkin(7) + t198 * t195 + t169;
t250 = m(6) * t147 - t156 * mrSges(6,1) + t157 * mrSges(6,2) - t180 * t173 + t181 * t174;
t242 = m(5) * t169 - t188 * mrSges(5,1) + t189 * mrSges(5,2) - t197 * t193 + t198 * t194 + t250;
t241 = -m(4) * t176 + mrSges(4,1) * t259 - t242 + (t221 * t238 + t264) * mrSges(4,3);
t134 = (mrSges(3,1) - t265) * qJDD(1) - t238 * mrSges(3,2) + m(3) * t191 + t241;
t110 = t229 * t113 + t231 * t134;
t115 = t230 * t121 + t228 * t122;
t252 = Ifges(4,5) * t228 + Ifges(4,6) * t230;
t262 = t238 * t252;
t258 = t231 * t113 - t229 * t134;
t158 = Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t223;
t131 = -mrSges(6,1) * t147 + mrSges(6,3) * t143 + Ifges(6,4) * t157 + Ifges(6,2) * t156 + Ifges(6,6) * t220 - t181 * t158 + t223 * t160;
t132 = mrSges(6,2) * t147 - mrSges(6,3) * t142 + Ifges(6,1) * t157 + Ifges(6,4) * t156 + Ifges(6,5) * t220 + t180 * t158 - t223 * t159;
t177 = Ifges(5,5) * t198 + Ifges(5,6) * t197 + Ifges(5,3) * qJD(4);
t116 = -mrSges(5,1) * t169 + mrSges(5,3) * t150 + Ifges(5,4) * t189 + Ifges(5,2) * t188 + Ifges(5,6) * qJDD(4) - pkin(4) * t250 + pkin(7) * t255 + qJD(4) * t179 + t235 * t131 + t232 * t132 - t198 * t177;
t117 = mrSges(5,2) * t169 - mrSges(5,3) * t149 + Ifges(5,1) * t189 + Ifges(5,4) * t188 + Ifges(5,5) * qJDD(4) - pkin(7) * t130 - qJD(4) * t178 - t232 * t131 + t235 * t132 + t197 * t177;
t104 = -mrSges(4,1) * t176 + mrSges(4,3) * t171 - pkin(3) * t242 + pkin(6) * t256 + t253 * qJDD(1) + t236 * t116 + t233 * t117 - t228 * t262;
t106 = mrSges(4,2) * t176 - mrSges(4,3) * t170 - pkin(6) * t123 + t254 * qJDD(1) - t233 * t116 + t236 * t117 + t230 * t262;
t245 = -mrSges(3,2) * t192 + qJ(3) * t257 + t230 * t104 + t228 * t106 + pkin(2) * (-qJDD(1) * t265 + t241) + mrSges(3,1) * t191 + Ifges(3,3) * qJDD(1);
t243 = mrSges(2,1) * t210 - mrSges(2,2) * t209 + Ifges(2,3) * qJDD(1) + pkin(1) * t110 + t245;
t108 = m(2) * t210 + qJDD(1) * mrSges(2,1) - t238 * mrSges(2,2) + t110;
t107 = m(2) * t209 - t238 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t258;
t102 = (Ifges(3,6) - t252) * qJDD(1) + t238 * Ifges(3,5) - mrSges(3,1) * t227 + mrSges(3,3) * t192 - pkin(2) * t115 + t267;
t101 = mrSges(3,2) * t227 - mrSges(3,3) * t191 + Ifges(3,5) * qJDD(1) - t238 * Ifges(3,6) - qJ(3) * t115 - t228 * t104 + t230 * t106;
t100 = -mrSges(2,2) * g(1) - mrSges(2,3) * t210 + Ifges(2,5) * qJDD(1) - t238 * Ifges(2,6) - qJ(2) * t110 + t231 * t101 - t229 * t102;
t99 = Ifges(2,6) * qJDD(1) + t238 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t209 + t229 * t101 + t231 * t102 - pkin(1) * (m(3) * t227 + t115) + qJ(2) * t258;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t243, t100, t101, t106, t117, t132; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t234 * t100 + t237 * t99 - pkin(5) * (-t237 * t107 + t234 * t108), t99, t102, t104, t116, t131; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t237 * t100 + t234 * t99 + pkin(5) * (t234 * t107 + t237 * t108), t243, t245, t252 * qJDD(1) - t267, -t240, -t244;];
m_new = t1;
