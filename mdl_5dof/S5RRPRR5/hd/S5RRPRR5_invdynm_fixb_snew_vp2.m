% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRR5
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
% Datum: 2022-01-20 11:03
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR5_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:02:24
% EndTime: 2022-01-20 11:02:34
% DurationCPUTime: 6.16s
% Computational Cost: add. (131840->248), mult. (182364->312), div. (0->0), fcn. (122991->10), ass. (0->112)
t225 = qJD(1) + qJD(2);
t219 = t225 ^ 2;
t233 = sin(qJ(1));
t237 = cos(qJ(1));
t211 = t233 * g(1) - t237 * g(2);
t206 = qJDD(1) * pkin(1) + t211;
t212 = -t237 * g(1) - t233 * g(2);
t238 = qJD(1) ^ 2;
t207 = -t238 * pkin(1) + t212;
t232 = sin(qJ(2));
t236 = cos(qJ(2));
t192 = t232 * t206 + t236 * t207;
t221 = qJDD(1) + qJDD(2);
t189 = -t219 * pkin(2) + t221 * qJ(3) + t192;
t228 = sin(pkin(9));
t229 = cos(pkin(9));
t261 = qJD(3) * t225;
t259 = -t229 * g(3) - 0.2e1 * t228 * t261;
t266 = pkin(3) * t229;
t167 = (-pkin(7) * t221 + t219 * t266 - t189) * t228 + t259;
t172 = -t228 * g(3) + (t189 + 0.2e1 * t261) * t229;
t263 = t221 * t229;
t223 = t229 ^ 2;
t264 = t219 * t223;
t168 = -pkin(3) * t264 + pkin(7) * t263 + t172;
t231 = sin(qJ(4));
t235 = cos(qJ(4));
t149 = t235 * t167 - t231 * t168;
t247 = t228 * t235 + t229 * t231;
t246 = -t228 * t231 + t229 * t235;
t197 = t246 * t225;
t260 = t197 * qJD(4);
t188 = t247 * t221 + t260;
t198 = t247 * t225;
t144 = (-t188 + t260) * pkin(8) + (t197 * t198 + qJDD(4)) * pkin(4) + t149;
t150 = t231 * t167 + t235 * t168;
t187 = -t198 * qJD(4) + t246 * t221;
t195 = qJD(4) * pkin(4) - t198 * pkin(8);
t196 = t197 ^ 2;
t145 = -t196 * pkin(4) + t187 * pkin(8) - qJD(4) * t195 + t150;
t230 = sin(qJ(5));
t234 = cos(qJ(5));
t142 = t234 * t144 - t230 * t145;
t178 = t234 * t197 - t230 * t198;
t157 = t178 * qJD(5) + t230 * t187 + t234 * t188;
t179 = t230 * t197 + t234 * t198;
t163 = -t178 * mrSges(6,1) + t179 * mrSges(6,2);
t224 = qJD(4) + qJD(5);
t173 = -t224 * mrSges(6,2) + t178 * mrSges(6,3);
t220 = qJDD(4) + qJDD(5);
t139 = m(6) * t142 + t220 * mrSges(6,1) - t157 * mrSges(6,3) - t179 * t163 + t224 * t173;
t143 = t230 * t144 + t234 * t145;
t156 = -t179 * qJD(5) + t234 * t187 - t230 * t188;
t174 = t224 * mrSges(6,1) - t179 * mrSges(6,3);
t140 = m(6) * t143 - t220 * mrSges(6,2) + t156 * mrSges(6,3) + t178 * t163 - t224 * t174;
t130 = t234 * t139 + t230 * t140;
t185 = -t197 * mrSges(5,1) + t198 * mrSges(5,2);
t193 = -qJD(4) * mrSges(5,2) + t197 * mrSges(5,3);
t127 = m(5) * t149 + qJDD(4) * mrSges(5,1) - t188 * mrSges(5,3) + qJD(4) * t193 - t198 * t185 + t130;
t194 = qJD(4) * mrSges(5,1) - t198 * mrSges(5,3);
t255 = -t230 * t139 + t234 * t140;
t128 = m(5) * t150 - qJDD(4) * mrSges(5,2) + t187 * mrSges(5,3) - qJD(4) * t194 + t197 * t185 + t255;
t123 = t235 * t127 + t231 * t128;
t171 = -t228 * t189 + t259;
t176 = Ifges(5,4) * t198 + Ifges(5,2) * t197 + Ifges(5,6) * qJD(4);
t177 = Ifges(5,1) * t198 + Ifges(5,4) * t197 + Ifges(5,5) * qJD(4);
t159 = Ifges(6,4) * t179 + Ifges(6,2) * t178 + Ifges(6,6) * t224;
t160 = Ifges(6,1) * t179 + Ifges(6,4) * t178 + Ifges(6,5) * t224;
t244 = -mrSges(6,1) * t142 + mrSges(6,2) * t143 - Ifges(6,5) * t157 - Ifges(6,6) * t156 - Ifges(6,3) * t220 - t179 * t159 + t178 * t160;
t240 = -mrSges(5,1) * t149 + mrSges(5,2) * t150 - Ifges(5,5) * t188 - Ifges(5,6) * t187 - Ifges(5,3) * qJDD(4) - pkin(4) * t130 - t198 * t176 + t197 * t177 + t244;
t253 = Ifges(4,4) * t228 + Ifges(4,2) * t229;
t254 = Ifges(4,1) * t228 + Ifges(4,4) * t229;
t267 = -mrSges(4,1) * t171 + mrSges(4,2) * t172 - pkin(3) * t123 - (t228 * t253 - t229 * t254) * t219 + t240;
t265 = mrSges(4,2) * t228;
t252 = Ifges(4,5) * t228 + Ifges(4,6) * t229;
t262 = t219 * t252;
t249 = mrSges(4,3) * t221 + (-mrSges(4,1) * t229 + t265) * t219;
t121 = m(4) * t171 - t249 * t228 + t123;
t256 = -t231 * t127 + t235 * t128;
t122 = m(4) * t172 + t249 * t229 + t256;
t257 = -t228 * t121 + t229 * t122;
t113 = m(3) * t192 - t219 * mrSges(3,1) - t221 * mrSges(3,2) + t257;
t191 = t236 * t206 - t232 * t207;
t251 = qJDD(3) - t191;
t186 = -t221 * pkin(2) - t219 * qJ(3) + t251;
t222 = t228 ^ 2;
t170 = (-pkin(2) - t266) * t221 + (-qJ(3) + (-t222 - t223) * pkin(7)) * t219 + t251;
t147 = -t187 * pkin(4) - t196 * pkin(8) + t198 * t195 + t170;
t250 = m(6) * t147 - t156 * mrSges(6,1) + t157 * mrSges(6,2) - t178 * t173 + t179 * t174;
t242 = m(5) * t170 - t187 * mrSges(5,1) + t188 * mrSges(5,2) - t197 * t193 + t198 * t194 + t250;
t241 = -m(4) * t186 + mrSges(4,1) * t263 - t242 + (t219 * t222 + t264) * mrSges(4,3);
t134 = (mrSges(3,1) - t265) * t221 + t241 - t219 * mrSges(3,2) + m(3) * t191;
t110 = t232 * t113 + t236 * t134;
t115 = t229 * t121 + t228 * t122;
t258 = t236 * t113 - t232 * t134;
t158 = Ifges(6,5) * t179 + Ifges(6,6) * t178 + Ifges(6,3) * t224;
t131 = -mrSges(6,1) * t147 + mrSges(6,3) * t143 + Ifges(6,4) * t157 + Ifges(6,2) * t156 + Ifges(6,6) * t220 - t179 * t158 + t224 * t160;
t132 = mrSges(6,2) * t147 - mrSges(6,3) * t142 + Ifges(6,1) * t157 + Ifges(6,4) * t156 + Ifges(6,5) * t220 + t178 * t158 - t224 * t159;
t175 = Ifges(5,5) * t198 + Ifges(5,6) * t197 + Ifges(5,3) * qJD(4);
t116 = -mrSges(5,1) * t170 + mrSges(5,3) * t150 + Ifges(5,4) * t188 + Ifges(5,2) * t187 + Ifges(5,6) * qJDD(4) - pkin(4) * t250 + pkin(8) * t255 + qJD(4) * t177 + t234 * t131 + t230 * t132 - t198 * t175;
t117 = mrSges(5,2) * t170 - mrSges(5,3) * t149 + Ifges(5,1) * t188 + Ifges(5,4) * t187 + Ifges(5,5) * qJDD(4) - pkin(8) * t130 - qJD(4) * t176 - t230 * t131 + t234 * t132 + t197 * t175;
t104 = -mrSges(4,1) * t186 + mrSges(4,3) * t172 - pkin(3) * t242 + pkin(7) * t256 + t235 * t116 + t231 * t117 + t253 * t221 - t228 * t262;
t106 = mrSges(4,2) * t186 - mrSges(4,3) * t171 - pkin(7) * t123 - t231 * t116 + t235 * t117 + t254 * t221 + t229 * t262;
t245 = -mrSges(3,2) * t192 + qJ(3) * t257 + t229 * t104 + t228 * t106 + pkin(2) * (-t221 * t265 + t241) + mrSges(3,1) * t191 + Ifges(3,3) * t221;
t243 = mrSges(2,1) * t211 - mrSges(2,2) * t212 + Ifges(2,3) * qJDD(1) + pkin(1) * t110 + t245;
t108 = m(2) * t212 - t238 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t258;
t107 = m(2) * t211 + qJDD(1) * mrSges(2,1) - t238 * mrSges(2,2) + t110;
t102 = mrSges(3,1) * g(3) + (Ifges(3,6) - t252) * t221 + t219 * Ifges(3,5) + mrSges(3,3) * t192 - pkin(2) * t115 + t267;
t101 = -mrSges(3,2) * g(3) - mrSges(3,3) * t191 + Ifges(3,5) * t221 - t219 * Ifges(3,6) - qJ(3) * t115 - t228 * t104 + t229 * t106;
t100 = -mrSges(2,2) * g(3) - mrSges(2,3) * t211 + Ifges(2,5) * qJDD(1) - t238 * Ifges(2,6) - pkin(6) * t110 + t236 * t101 - t232 * t102;
t99 = Ifges(2,6) * qJDD(1) + t238 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t212 + t232 * t101 + t236 * t102 - pkin(1) * (-m(3) * g(3) + t115) + pkin(6) * t258;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t237 * t100 - t233 * t99 - pkin(5) * (t237 * t107 + t233 * t108), t100, t101, t106, t117, t132; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t233 * t100 + t237 * t99 + pkin(5) * (-t233 * t107 + t237 * t108), t99, t102, t104, t116, t131; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t243, t243, t245, t252 * t221 - t267, -t240, -t244;];
m_new = t1;
