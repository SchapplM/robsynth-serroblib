% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:46
% EndTime: 2019-12-31 18:23:50
% DurationCPUTime: 2.54s
% Computational Cost: add. (26658->274), mult. (51840->328), div. (0->0), fcn. (26079->8), ass. (0->110)
t230 = sin(qJ(1));
t233 = cos(qJ(1));
t212 = t230 * g(1) - t233 * g(2);
t198 = qJDD(1) * pkin(1) + t212;
t213 = -t233 * g(1) - t230 * g(2);
t235 = qJD(1) ^ 2;
t202 = -t235 * pkin(1) + t213;
t226 = sin(pkin(8));
t227 = cos(pkin(8));
t171 = t226 * t198 + t227 * t202;
t159 = -t235 * pkin(2) + qJDD(1) * pkin(6) + t171;
t229 = sin(qJ(3));
t156 = t229 * t159;
t225 = -g(3) + qJDD(2);
t232 = cos(qJ(3));
t260 = t232 * t225;
t154 = -t156 + t260;
t155 = t232 * t159 + t229 * t225;
t183 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t229 + Ifges(4,4) * t232) * qJD(1);
t200 = (mrSges(5,2) * t232 - mrSges(5,3) * t229) * qJD(1);
t255 = qJD(1) * qJD(3);
t253 = t232 * t255;
t203 = t229 * qJDD(1) + t253;
t254 = t229 * t255;
t204 = t232 * qJDD(1) - t254;
t257 = qJD(1) * t232;
t209 = -mrSges(5,1) * t257 - qJD(3) * mrSges(5,3);
t256 = t229 * qJD(1);
t211 = pkin(4) * t256 - qJD(3) * pkin(7);
t224 = t232 ^ 2;
t170 = t227 * t198 - t226 * t202;
t249 = -qJDD(1) * pkin(2) - t170;
t266 = -2 * qJD(4);
t239 = pkin(3) * t254 + t256 * t266 + (-t203 - t253) * qJ(4) + t249;
t265 = -pkin(3) - pkin(7);
t145 = -t211 * t256 + (-pkin(4) * t224 - pkin(6)) * t235 + t265 * t204 + t239;
t199 = (-pkin(3) * t232 - qJ(4) * t229) * qJD(1);
t234 = qJD(3) ^ 2;
t250 = -t234 * qJ(4) + t199 * t256 + qJDD(4) + t156;
t264 = pkin(7) * t235;
t148 = t203 * pkin(4) + t265 * qJDD(3) + (-pkin(4) * t255 - t229 * t264 - t225) * t232 + t250;
t228 = sin(qJ(5));
t231 = cos(qJ(5));
t143 = -t228 * t145 + t231 * t148;
t196 = -t228 * qJD(3) - t231 * t257;
t168 = t196 * qJD(5) + t231 * qJDD(3) - t228 * t204;
t197 = t231 * qJD(3) - t228 * t257;
t172 = -t196 * mrSges(6,1) + t197 * mrSges(6,2);
t215 = qJD(5) + t256;
t173 = -t215 * mrSges(6,2) + t196 * mrSges(6,3);
t195 = qJDD(5) + t203;
t139 = m(6) * t143 + t195 * mrSges(6,1) - t168 * mrSges(6,3) - t197 * t172 + t215 * t173;
t144 = t231 * t145 + t228 * t148;
t167 = -t197 * qJD(5) - t228 * qJDD(3) - t231 * t204;
t174 = t215 * mrSges(6,1) - t197 * mrSges(6,3);
t140 = m(6) * t144 - t195 * mrSges(6,2) + t167 * mrSges(6,3) + t196 * t172 - t215 * t174;
t128 = t231 * t139 + t228 * t140;
t244 = -t234 * pkin(3) + qJDD(3) * qJ(4) + t199 * t257 + t155;
t147 = -t224 * t264 + t204 * pkin(4) + ((2 * qJD(4)) + t211) * qJD(3) + t244;
t160 = Ifges(6,5) * t197 + Ifges(6,6) * t196 + Ifges(6,3) * t215;
t162 = Ifges(6,1) * t197 + Ifges(6,4) * t196 + Ifges(6,5) * t215;
t131 = -mrSges(6,1) * t147 + mrSges(6,3) * t144 + Ifges(6,4) * t168 + Ifges(6,2) * t167 + Ifges(6,6) * t195 - t197 * t160 + t215 * t162;
t161 = Ifges(6,4) * t197 + Ifges(6,2) * t196 + Ifges(6,6) * t215;
t132 = mrSges(6,2) * t147 - mrSges(6,3) * t143 + Ifges(6,1) * t168 + Ifges(6,4) * t167 + Ifges(6,5) * t195 + t196 * t160 - t215 * t161;
t150 = qJD(3) * t266 - t244;
t152 = -qJDD(3) * pkin(3) + t250 - t260;
t185 = Ifges(5,4) * qJD(3) + (-Ifges(5,2) * t229 - Ifges(5,6) * t232) * qJD(1);
t242 = -mrSges(5,2) * t152 + mrSges(5,3) * t150 - Ifges(5,1) * qJDD(3) + Ifges(5,4) * t203 + Ifges(5,5) * t204 + pkin(7) * t128 + t228 * t131 - t231 * t132 - t185 * t257;
t141 = -m(6) * t147 + t167 * mrSges(6,1) - t168 * mrSges(6,2) + t196 * t173 - t197 * t174;
t210 = mrSges(5,1) * t256 + qJD(3) * mrSges(5,2);
t243 = -m(5) * t150 + qJDD(3) * mrSges(5,3) + qJD(3) * t210 + t200 * t257 - t141;
t247 = -m(5) * t152 - t203 * mrSges(5,1) - t128;
t184 = Ifges(5,5) * qJD(3) + (-Ifges(5,6) * t229 - Ifges(5,3) * t232) * qJD(1);
t258 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t229 + Ifges(4,2) * t232) * qJD(1) - t184;
t268 = (-t232 * t183 + t258 * t229) * qJD(1) + mrSges(4,1) * t154 - mrSges(4,2) * t155 + Ifges(4,5) * t203 + Ifges(4,6) * t204 + Ifges(4,3) * qJDD(3) + pkin(3) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t209 - t200 * t256 + t247) + qJ(4) * (t204 * mrSges(5,1) + t243) - t242;
t263 = t235 * pkin(6);
t262 = Ifges(4,4) + Ifges(5,6);
t201 = (-mrSges(4,1) * t232 + mrSges(4,2) * t229) * qJD(1);
t208 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t257;
t125 = m(4) * t154 - t203 * mrSges(4,3) + (mrSges(4,1) - mrSges(5,2)) * qJDD(3) + (t208 - t209) * qJD(3) + (-t200 - t201) * t256 + t247;
t207 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t256;
t135 = t201 * t257 + m(4) * t155 - qJDD(3) * mrSges(4,2) - qJD(3) * t207 + (mrSges(4,3) + mrSges(5,1)) * t204 + t243;
t251 = -t229 * t125 + t232 * t135;
t118 = m(3) * t171 - t235 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t251;
t158 = t249 - t263;
t129 = -t228 * t139 + t231 * t140;
t149 = -t204 * pkin(3) + t239 - t263;
t248 = -m(5) * t149 - t204 * mrSges(5,2) + t210 * t256 - t129;
t237 = -m(4) * t158 + t208 * t257 + t204 * mrSges(4,1) + (-mrSges(4,2) + mrSges(5,3)) * t203 + (-t207 * t229 - t209 * t232) * qJD(1) + t248;
t122 = m(3) * t170 + qJDD(1) * mrSges(3,1) - t235 * mrSges(3,2) + t237;
t115 = t226 * t118 + t227 * t122;
t120 = t232 * t125 + t229 * t135;
t186 = Ifges(5,1) * qJD(3) + (-Ifges(5,4) * t229 - Ifges(5,5) * t232) * qJD(1);
t259 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t229 + Ifges(4,6) * t232) * qJD(1) + t186;
t252 = t227 * t118 - t226 * t122;
t126 = -t203 * mrSges(5,3) + t209 * t257 - t248;
t240 = -mrSges(5,1) * t150 + mrSges(5,2) * t149 - pkin(4) * t141 - pkin(7) * t129 - t231 * t131 - t228 * t132;
t109 = -mrSges(4,1) * t158 + mrSges(4,3) * t155 - pkin(3) * t126 + (Ifges(4,2) + Ifges(5,3)) * t204 + t262 * t203 + (Ifges(4,6) - Ifges(5,5)) * qJDD(3) + (t183 - t185) * qJD(3) - t259 * t256 + t240;
t245 = mrSges(6,1) * t143 - mrSges(6,2) * t144 + Ifges(6,5) * t168 + Ifges(6,6) * t167 + Ifges(6,3) * t195 + t197 * t161 - t196 * t162;
t238 = mrSges(5,1) * t152 - mrSges(5,3) * t149 + pkin(4) * t128 + t245;
t111 = t238 - mrSges(4,3) * t154 + mrSges(4,2) * t158 + t259 * t257 - qJ(4) * t126 + t262 * t204 + (Ifges(4,1) + Ifges(5,2)) * t203 + (Ifges(4,5) - Ifges(5,4)) * qJDD(3) - t258 * qJD(3);
t246 = mrSges(3,1) * t170 - mrSges(3,2) * t171 + Ifges(3,3) * qJDD(1) + pkin(2) * t237 + pkin(6) * t251 + t232 * t109 + t229 * t111;
t241 = mrSges(2,1) * t212 - mrSges(2,2) * t213 + Ifges(2,3) * qJDD(1) + pkin(1) * t115 + t246;
t113 = m(2) * t213 - t235 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t252;
t112 = m(2) * t212 + qJDD(1) * mrSges(2,1) - t235 * mrSges(2,2) + t115;
t107 = -mrSges(3,1) * t225 + mrSges(3,3) * t171 + t235 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t120 - t268;
t106 = mrSges(3,2) * t225 - mrSges(3,3) * t170 + Ifges(3,5) * qJDD(1) - t235 * Ifges(3,6) - pkin(6) * t120 - t229 * t109 + t232 * t111;
t105 = -mrSges(2,2) * g(3) - mrSges(2,3) * t212 + Ifges(2,5) * qJDD(1) - t235 * Ifges(2,6) - qJ(2) * t115 + t227 * t106 - t226 * t107;
t104 = Ifges(2,6) * qJDD(1) + t235 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t213 + t226 * t106 + t227 * t107 - pkin(1) * (m(3) * t225 + t120) + qJ(2) * t252;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t233 * t105 - t230 * t104 - pkin(5) * (t233 * t112 + t230 * t113), t105, t106, t111, -t184 * t256 - t242, t132; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t230 * t105 + t233 * t104 + pkin(5) * (-t230 * t112 + t233 * t113), t104, t107, t109, Ifges(5,4) * qJDD(3) - Ifges(5,2) * t203 - Ifges(5,6) * t204 - qJD(3) * t184 - t186 * t257 - t238, t131; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t241, t241, t246, t268, Ifges(5,5) * qJDD(3) - Ifges(5,6) * t203 - Ifges(5,3) * t204 + qJD(3) * t185 + t186 * t256 - t240, t245;];
m_new = t1;
