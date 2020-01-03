% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR16
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
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
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR16_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR16_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:38:47
% EndTime: 2019-12-31 18:38:51
% DurationCPUTime: 1.82s
% Computational Cost: add. (16399->272), mult. (32060->319), div. (0->0), fcn. (14991->6), ass. (0->110)
t224 = sin(qJ(1));
t227 = cos(qJ(1));
t203 = t224 * g(1) - t227 * g(2);
t229 = qJD(1) ^ 2;
t247 = -t229 * qJ(2) + qJDD(2) - t203;
t273 = -pkin(1) - pkin(6);
t165 = t273 * qJDD(1) + t247;
t226 = cos(qJ(3));
t262 = t226 * t165;
t223 = sin(qJ(3));
t269 = t223 * g(3);
t158 = t262 + t269;
t259 = qJD(1) * qJD(3);
t209 = t223 * t259;
t195 = qJDD(1) * t226 - t209;
t261 = qJD(1) * t223;
t198 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t261;
t200 = mrSges(5,1) * t261 - qJD(3) * mrSges(5,3);
t256 = t226 * t259;
t194 = qJDD(1) * t223 + t256;
t260 = qJD(1) * t226;
t202 = pkin(4) * t260 - qJD(3) * pkin(7);
t219 = t223 ^ 2;
t204 = -t227 * g(1) - t224 * g(2);
t249 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t204;
t274 = -2 * qJD(4);
t239 = pkin(3) * t256 + t260 * t274 + t249 + (-t195 + t209) * qJ(4);
t272 = pkin(3) + pkin(7);
t137 = -t202 * t260 + t272 * t194 + (-pkin(4) * t219 + t273) * t229 + t239;
t191 = (pkin(3) * t223 - qJ(4) * t226) * qJD(1);
t228 = qJD(3) ^ 2;
t246 = -t228 * qJ(4) + t191 * t260 + qJDD(4) - t262;
t270 = pkin(7) * t229;
t140 = t195 * pkin(4) - t272 * qJDD(3) + (pkin(4) * t259 + t226 * t270 - g(3)) * t223 + t246;
t222 = sin(qJ(5));
t225 = cos(qJ(5));
t135 = -t137 * t222 + t140 * t225;
t189 = -qJD(3) * t222 + t225 * t261;
t155 = t189 * qJD(5) + qJDD(3) * t225 + t194 * t222;
t190 = qJD(3) * t225 + t222 * t261;
t157 = -mrSges(6,1) * t189 + mrSges(6,2) * t190;
t207 = qJD(5) + t260;
t160 = -mrSges(6,2) * t207 + t189 * mrSges(6,3);
t188 = qJDD(5) + t195;
t131 = m(6) * t135 + t188 * mrSges(6,1) - t155 * mrSges(6,3) - t190 * t157 + t160 * t207;
t136 = t137 * t225 + t140 * t222;
t154 = -t190 * qJD(5) - qJDD(3) * t222 + t194 * t225;
t161 = mrSges(6,1) * t207 - t190 * mrSges(6,3);
t132 = m(6) * t136 - t188 * mrSges(6,2) + t154 * mrSges(6,3) + t189 * t157 - t161 * t207;
t121 = t225 * t131 + t222 * t132;
t145 = -qJDD(3) * pkin(3) + t246 - t269;
t244 = -m(5) * t145 - t195 * mrSges(5,1) - t121;
t192 = (-mrSges(5,2) * t223 - mrSges(5,3) * t226) * qJD(1);
t254 = qJD(1) * (-t192 - (mrSges(4,1) * t223 + mrSges(4,2) * t226) * qJD(1));
t267 = mrSges(4,1) - mrSges(5,2);
t117 = m(4) * t158 - t195 * mrSges(4,3) + t267 * qJDD(3) + (t198 - t200) * qJD(3) + t226 * t254 + t244;
t159 = -t226 * g(3) + t223 * t165;
t199 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t260;
t240 = -t228 * pkin(3) + qJDD(3) * qJ(4) - t191 * t261 + t159;
t139 = -t219 * t270 - t194 * pkin(4) + ((2 * qJD(4)) + t202) * qJD(3) + t240;
t133 = -m(6) * t139 + t154 * mrSges(6,1) - t155 * mrSges(6,2) + t189 * t160 - t190 * t161;
t143 = qJD(3) * t274 - t240;
t201 = mrSges(5,1) * t260 + qJD(3) * mrSges(5,2);
t238 = -m(5) * t143 + qJDD(3) * mrSges(5,3) + qJD(3) * t201 - t133;
t127 = m(4) * t159 - qJDD(3) * mrSges(4,2) - qJD(3) * t199 + (-mrSges(4,3) - mrSges(5,1)) * t194 + t223 * t254 + t238;
t112 = t226 * t117 + t223 * t127;
t170 = -qJDD(1) * pkin(1) + t247;
t173 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t226 - Ifges(4,2) * t223) * qJD(1);
t174 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t226 - Ifges(4,4) * t223) * qJD(1);
t147 = Ifges(6,5) * t190 + Ifges(6,6) * t189 + Ifges(6,3) * t207;
t149 = Ifges(6,1) * t190 + Ifges(6,4) * t189 + Ifges(6,5) * t207;
t124 = -mrSges(6,1) * t139 + mrSges(6,3) * t136 + Ifges(6,4) * t155 + Ifges(6,2) * t154 + Ifges(6,6) * t188 - t190 * t147 + t149 * t207;
t148 = Ifges(6,4) * t190 + Ifges(6,2) * t189 + Ifges(6,6) * t207;
t125 = mrSges(6,2) * t139 - mrSges(6,3) * t135 + Ifges(6,1) * t155 + Ifges(6,4) * t154 + Ifges(6,5) * t188 + t189 * t147 - t148 * t207;
t175 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t226 + Ifges(5,3) * t223) * qJD(1);
t176 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t226 + Ifges(5,6) * t223) * qJD(1);
t277 = mrSges(5,2) * t145 - mrSges(5,3) * t143 + Ifges(5,1) * qJDD(3) - Ifges(5,4) * t195 + Ifges(5,5) * t194 - pkin(7) * t121 - t222 * t124 + t225 * t125 - (t175 * t226 + t176 * t223) * qJD(1);
t230 = -mrSges(4,2) * t159 + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t200 - t192 * t260 + t244) + qJ(4) * (-t194 * mrSges(5,1) - t192 * t261 + t238) + mrSges(4,1) * t158 + t174 * t261 + t173 * t260 - Ifges(4,6) * t194 + Ifges(4,5) * t195 + Ifges(4,3) * qJDD(3) + t277;
t279 = mrSges(3,1) * t170 + pkin(2) * t112 + t230;
t257 = t273 * t229;
t164 = t257 + t249;
t122 = -t222 * t131 + t225 * t132;
t142 = t194 * pkin(3) + t239 + t257;
t236 = m(5) * t142 - t195 * mrSges(5,3) - (t200 * t223 + t201 * t226) * qJD(1) + t122;
t278 = -m(4) * t164 - t195 * mrSges(4,2) - t267 * t194 - t198 * t261 - t199 * t260 - t236;
t268 = mrSges(2,1) - mrSges(3,2);
t266 = Ifges(4,4) + Ifges(5,6);
t265 = Ifges(2,5) - Ifges(3,4);
t264 = -Ifges(2,6) + Ifges(3,5);
t113 = -t117 * t223 + t226 * t127;
t177 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t226 + Ifges(5,5) * t223) * qJD(1);
t255 = qJD(1) * (-Ifges(4,3) * qJD(3) - (Ifges(4,5) * t226 - Ifges(4,6) * t223) * qJD(1) - t177);
t245 = -m(3) * t170 + t229 * mrSges(3,3) - t112;
t243 = mrSges(6,1) * t135 - mrSges(6,2) * t136 + Ifges(6,5) * t155 + Ifges(6,6) * t154 + Ifges(6,3) * t188 + t190 * t148 - t189 * t149;
t118 = -t194 * mrSges(5,2) + t236;
t235 = -mrSges(5,1) * t143 + mrSges(5,2) * t142 - pkin(4) * t133 - pkin(7) * t122 - t225 * t124 - t222 * t125;
t106 = -mrSges(4,1) * t164 + mrSges(4,3) * t159 - pkin(3) * t118 + t266 * t195 + (-Ifges(4,2) - Ifges(5,3)) * t194 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t174 - t176) * qJD(3) + t226 * t255 + t235;
t233 = mrSges(5,1) * t145 - mrSges(5,3) * t142 + pkin(4) * t121 + t243;
t108 = (Ifges(4,1) + Ifges(5,2)) * t195 - t266 * t194 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) + (-t173 + t175) * qJD(3) - mrSges(4,3) * t158 + mrSges(4,2) * t164 - qJ(4) * t118 + t223 * t255 + t233;
t168 = pkin(1) * t229 - t249;
t242 = mrSges(3,2) * t170 - mrSges(3,3) * t168 + Ifges(3,1) * qJDD(1) - pkin(6) * t112 - t106 * t223 + t226 * t108;
t241 = -mrSges(3,1) * t168 - pkin(2) * t278 - pkin(6) * t113 - t226 * t106 - t223 * t108;
t231 = -m(3) * t168 + t229 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t278;
t234 = -mrSges(2,2) * t204 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t245) + qJ(2) * t231 + mrSges(2,1) * t203 + Ifges(2,3) * qJDD(1) + t242;
t114 = m(2) * t204 - mrSges(2,1) * t229 - qJDD(1) * mrSges(2,2) + t231;
t111 = -m(3) * g(3) + t113;
t109 = m(2) * t203 - t229 * mrSges(2,2) + t268 * qJDD(1) + t245;
t105 = t264 * t229 + t265 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t203 - qJ(2) * t111 + t279;
t104 = mrSges(2,3) * t204 - pkin(1) * t111 + t268 * g(3) - t264 * qJDD(1) + t265 * t229 + t241;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t227 * t105 - t224 * t104 - pkin(5) * (t109 * t227 + t114 * t224), t105, t242, t108, t277, t125; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t224 * t105 + t227 * t104 + pkin(5) * (-t109 * t224 + t114 * t227), t104, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - t229 * Ifges(3,5) - t279, t106, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t195 + Ifges(5,6) * t194 - qJD(3) * t175 + t177 * t261 - t233, t124; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t234, t234, mrSges(3,2) * g(3) + t229 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t241, t230, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t195 + Ifges(5,3) * t194 + qJD(3) * t176 + t177 * t260 - t235, t243;];
m_new = t1;
