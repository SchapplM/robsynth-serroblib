% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR14_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:35
% EndTime: 2019-12-31 18:34:41
% DurationCPUTime: 3.60s
% Computational Cost: add. (43347->268), mult. (92843->331), div. (0->0), fcn. (56509->8), ass. (0->109)
t233 = sin(qJ(1));
t236 = cos(qJ(1));
t215 = -t236 * g(1) - t233 * g(2);
t251 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t215;
t229 = sin(pkin(8));
t230 = cos(pkin(8));
t232 = sin(qJ(3));
t235 = cos(qJ(3));
t196 = (t229 * t235 + t230 * t232) * qJD(1);
t266 = 2 * qJD(4);
t265 = -pkin(1) - pkin(6);
t264 = mrSges(2,1) - mrSges(3,2);
t263 = Ifges(2,5) - Ifges(3,4);
t262 = (-Ifges(2,6) + Ifges(3,5));
t214 = t233 * g(1) - t236 * g(2);
t238 = qJD(1) ^ 2;
t250 = -t238 * qJ(2) + qJDD(2) - t214;
t188 = t265 * qJDD(1) + t250;
t178 = t232 * g(3) + t235 * t188;
t259 = qJD(1) * qJD(3);
t257 = t232 * t259;
t209 = t235 * qJDD(1) - t257;
t157 = (-t209 - t257) * qJ(4) + (-t232 * t235 * t238 + qJDD(3)) * pkin(3) + t178;
t179 = -t235 * g(3) + t232 * t188;
t208 = -t232 * qJDD(1) - t235 * t259;
t260 = qJD(1) * t235;
t212 = qJD(3) * pkin(3) - qJ(4) * t260;
t226 = t232 ^ 2;
t158 = -t226 * t238 * pkin(3) + t208 * qJ(4) - qJD(3) * t212 + t179;
t147 = t229 * t157 + t230 * t158 - t196 * t266;
t261 = qJD(1) * t232;
t197 = -t229 * t261 + t230 * t260;
t169 = t196 * mrSges(5,1) + t197 * mrSges(5,2);
t175 = t230 * t208 - t229 * t209;
t187 = qJD(3) * mrSges(5,1) - t197 * mrSges(5,3);
t170 = t196 * pkin(4) - t197 * pkin(7);
t237 = qJD(3) ^ 2;
t143 = -t237 * pkin(4) + qJDD(3) * pkin(7) - t196 * t170 + t147;
t160 = -t208 * pkin(3) + qJDD(4) + t212 * t260 + (-qJ(4) * t226 + t265) * t238 + t251;
t176 = t229 * t208 + t230 * t209;
t144 = (qJD(3) * t196 - t176) * pkin(7) + (qJD(3) * t197 - t175) * pkin(4) + t160;
t231 = sin(qJ(5));
t234 = cos(qJ(5));
t140 = -t231 * t143 + t234 * t144;
t180 = t234 * qJD(3) - t231 * t197;
t154 = t180 * qJD(5) + t231 * qJDD(3) + t234 * t176;
t181 = t231 * qJD(3) + t234 * t197;
t162 = -t180 * mrSges(6,1) + t181 * mrSges(6,2);
t194 = qJD(5) + t196;
t163 = -t194 * mrSges(6,2) + t180 * mrSges(6,3);
t174 = qJDD(5) - t175;
t136 = m(6) * t140 + t174 * mrSges(6,1) - t154 * mrSges(6,3) - t181 * t162 + t194 * t163;
t141 = t234 * t143 + t231 * t144;
t153 = -t181 * qJD(5) + t234 * qJDD(3) - t231 * t176;
t164 = t194 * mrSges(6,1) - t181 * mrSges(6,3);
t137 = m(6) * t141 - t174 * mrSges(6,2) + t153 * mrSges(6,3) + t180 * t162 - t194 * t164;
t255 = -t231 * t136 + t234 * t137;
t123 = m(5) * t147 - qJDD(3) * mrSges(5,2) + t175 * mrSges(5,3) - qJD(3) * t187 - t196 * t169 + t255;
t254 = -t230 * t157 + t229 * t158;
t146 = -0.2e1 * qJD(4) * t197 - t254;
t186 = -qJD(3) * mrSges(5,2) - t196 * mrSges(5,3);
t142 = -qJDD(3) * pkin(4) - t237 * pkin(7) + (t266 + t170) * t197 + t254;
t248 = -m(6) * t142 + t153 * mrSges(6,1) - t154 * mrSges(6,2) + t180 * t163 - t181 * t164;
t132 = m(5) * t146 + qJDD(3) * mrSges(5,1) - t176 * mrSges(5,3) + qJD(3) * t186 - t197 * t169 + t248;
t117 = t229 * t123 + t230 * t132;
t125 = t234 * t136 + t231 * t137;
t207 = (mrSges(4,1) * t232 + mrSges(4,2) * t235) * qJD(1);
t211 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t261;
t114 = m(4) * t178 + qJDD(3) * mrSges(4,1) - t209 * mrSges(4,3) + qJD(3) * t211 - t207 * t260 + t117;
t213 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t260;
t256 = t230 * t123 - t229 * t132;
t115 = m(4) * t179 - qJDD(3) * mrSges(4,2) + t208 * mrSges(4,3) - qJD(3) * t213 - t207 * t261 + t256;
t110 = -t232 * t114 + t235 * t115;
t109 = t235 * t114 + t232 * t115;
t195 = -qJDD(1) * pkin(1) + t250;
t249 = -m(3) * t195 + (t238 * mrSges(3,3)) - t109;
t247 = m(5) * t160 - t175 * mrSges(5,1) + t176 * mrSges(5,2) + t196 * t186 + t197 * t187 + t125;
t148 = Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t194;
t150 = Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t194;
t129 = -mrSges(6,1) * t142 + mrSges(6,3) * t141 + Ifges(6,4) * t154 + Ifges(6,2) * t153 + Ifges(6,6) * t174 - t181 * t148 + t194 * t150;
t149 = Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t194;
t130 = mrSges(6,2) * t142 - mrSges(6,3) * t140 + Ifges(6,1) * t154 + Ifges(6,4) * t153 + Ifges(6,5) * t174 + t180 * t148 - t194 * t149;
t165 = Ifges(5,5) * t197 - Ifges(5,6) * t196 + (Ifges(5,3) * qJD(3));
t166 = Ifges(5,4) * t197 - Ifges(5,2) * t196 + Ifges(5,6) * qJD(3);
t111 = mrSges(5,2) * t160 - mrSges(5,3) * t146 + Ifges(5,1) * t176 + Ifges(5,4) * t175 + Ifges(5,5) * qJDD(3) - pkin(7) * t125 - qJD(3) * t166 - t231 * t129 + t234 * t130 - t196 * t165;
t167 = Ifges(5,1) * t197 - Ifges(5,4) * t196 + Ifges(5,5) * qJD(3);
t242 = mrSges(6,1) * t140 - mrSges(6,2) * t141 + Ifges(6,5) * t154 + Ifges(6,6) * t153 + Ifges(6,3) * t174 + t181 * t149 - t180 * t150;
t112 = -mrSges(5,1) * t160 + mrSges(5,3) * t147 + Ifges(5,4) * t176 + Ifges(5,2) * t175 + Ifges(5,6) * qJDD(3) - pkin(4) * t125 + qJD(3) * t167 - t197 * t165 - t242;
t185 = t265 * t238 + t251;
t198 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t235 - Ifges(4,6) * t232) * qJD(1);
t200 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t235 - Ifges(4,4) * t232) * qJD(1);
t103 = -mrSges(4,1) * t185 + mrSges(4,3) * t179 + Ifges(4,4) * t209 + Ifges(4,2) * t208 + Ifges(4,6) * qJDD(3) - pkin(3) * t247 + qJ(4) * t256 + qJD(3) * t200 + t229 * t111 + t230 * t112 - t198 * t260;
t199 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t235 - Ifges(4,2) * t232) * qJD(1);
t105 = mrSges(4,2) * t185 - mrSges(4,3) * t178 + Ifges(4,1) * t209 + Ifges(4,4) * t208 + Ifges(4,5) * qJDD(3) - qJ(4) * t117 - qJD(3) * t199 + t230 * t111 - t229 * t112 - t198 * t261;
t191 = t238 * pkin(1) - t251;
t246 = mrSges(3,2) * t195 - mrSges(3,3) * t191 + Ifges(3,1) * qJDD(1) - pkin(6) * t109 - t232 * t103 + t235 * t105;
t120 = -m(4) * t185 + t208 * mrSges(4,1) - t209 * mrSges(4,2) - t211 * t261 - t213 * t260 - t247;
t245 = -mrSges(3,1) * t191 - pkin(2) * t120 - pkin(6) * t110 - t235 * t103 - t232 * t105;
t244 = -mrSges(5,1) * t146 + mrSges(5,2) * t147 - Ifges(5,5) * t176 - Ifges(5,6) * t175 - Ifges(5,3) * qJDD(3) - pkin(4) * t248 - pkin(7) * t255 - t234 * t129 - t231 * t130 - t197 * t166 - t196 * t167;
t241 = -m(3) * t191 + t238 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t120;
t243 = -mrSges(2,2) * t215 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t249) + qJ(2) * t241 + mrSges(2,1) * t214 + Ifges(2,3) * qJDD(1) + t246;
t240 = -mrSges(4,1) * t178 + mrSges(4,2) * t179 - Ifges(4,5) * t209 - Ifges(4,6) * t208 - Ifges(4,3) * qJDD(3) - pkin(3) * t117 - t199 * t260 - t200 * t261 + t244;
t239 = -mrSges(3,1) * t195 - pkin(2) * t109 + t240;
t118 = m(2) * t215 - t238 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t241;
t108 = -m(3) * g(3) + t110;
t106 = m(2) * t214 - t238 * mrSges(2,2) + t264 * qJDD(1) + t249;
t102 = -t239 + (t262 * t238) + t263 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t214 - qJ(2) * t108;
t101 = mrSges(2,3) * t215 - pkin(1) * t108 + t264 * g(3) - t262 * qJDD(1) + t263 * t238 + t245;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t236 * t102 - t233 * t101 - pkin(5) * (t236 * t106 + t233 * t118), t102, t246, t105, t111, t130; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t233 * t102 + t236 * t101 + pkin(5) * (-t233 * t106 + t236 * t118), t101, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t238 * Ifges(3,5)) + t239, t103, t112, t129; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t243, t243, mrSges(3,2) * g(3) + t238 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t245, -t240, -t244, t242;];
m_new = t1;
