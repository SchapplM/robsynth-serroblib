% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:44
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:43:37
% EndTime: 2019-12-05 16:43:44
% DurationCPUTime: 3.44s
% Computational Cost: add. (37505->251), mult. (76098->312), div. (0->0), fcn. (47650->8), ass. (0->97)
t233 = sin(qJ(4));
t234 = sin(qJ(3));
t236 = cos(qJ(4));
t237 = cos(qJ(3));
t202 = (-t233 * t234 + t236 * t237) * qJD(2);
t259 = qJD(2) * qJD(3);
t256 = t237 * t259;
t211 = t234 * qJDD(2) + t256;
t212 = t237 * qJDD(2) - t234 * t259;
t171 = t202 * qJD(4) + t236 * t211 + t233 * t212;
t203 = (t233 * t237 + t234 * t236) * qJD(2);
t185 = -t202 * mrSges(6,1) + t203 * mrSges(6,2);
t231 = sin(pkin(8));
t232 = cos(pkin(8));
t214 = t231 * g(1) - t232 * g(2);
t215 = -t232 * g(1) - t231 * g(2);
t235 = sin(qJ(2));
t238 = cos(qJ(2));
t192 = t235 * t214 + t238 * t215;
t239 = qJD(2) ^ 2;
t188 = -t239 * pkin(2) + qJDD(2) * pkin(6) + t192;
t230 = -g(3) + qJDD(1);
t173 = -t234 * t188 + t237 * t230;
t151 = (-t211 + t256) * pkin(7) + (t234 * t237 * t239 + qJDD(3)) * pkin(3) + t173;
t174 = t237 * t188 + t234 * t230;
t261 = qJD(2) * t234;
t218 = qJD(3) * pkin(3) - pkin(7) * t261;
t229 = t237 ^ 2;
t152 = -t229 * t239 * pkin(3) + t212 * pkin(7) - qJD(3) * t218 + t174;
t146 = t236 * t151 - t233 * t152;
t226 = qJDD(3) + qJDD(4);
t227 = qJD(3) + qJD(4);
t138 = -0.2e1 * qJD(5) * t203 + (t202 * t227 - t171) * qJ(5) + (t202 * t203 + t226) * pkin(4) + t146;
t193 = -t227 * mrSges(6,2) + t202 * mrSges(6,3);
t258 = m(6) * t138 + t226 * mrSges(6,1) + t227 * t193;
t135 = -t171 * mrSges(6,3) - t203 * t185 + t258;
t147 = t233 * t151 + t236 * t152;
t170 = -t203 * qJD(4) - t233 * t211 + t236 * t212;
t178 = Ifges(5,4) * t203 + Ifges(5,2) * t202 + Ifges(5,6) * t227;
t179 = Ifges(6,1) * t203 + Ifges(6,4) * t202 + Ifges(6,5) * t227;
t180 = Ifges(5,1) * t203 + Ifges(5,4) * t202 + Ifges(5,5) * t227;
t195 = t227 * pkin(4) - t203 * qJ(5);
t198 = t202 ^ 2;
t141 = -t198 * pkin(4) + t170 * qJ(5) + 0.2e1 * qJD(5) * t202 - t227 * t195 + t147;
t177 = Ifges(6,4) * t203 + Ifges(6,2) * t202 + Ifges(6,6) * t227;
t247 = -mrSges(6,1) * t138 + mrSges(6,2) * t141 - Ifges(6,5) * t171 - Ifges(6,6) * t170 - Ifges(6,3) * t226 - t203 * t177;
t265 = mrSges(5,1) * t146 - mrSges(5,2) * t147 + Ifges(5,5) * t171 + Ifges(5,6) * t170 + Ifges(5,3) * t226 + pkin(4) * t135 + t203 * t178 - t247 + (-t180 - t179) * t202;
t186 = -t202 * mrSges(5,1) + t203 * mrSges(5,2);
t194 = -t227 * mrSges(5,2) + t202 * mrSges(5,3);
t129 = m(5) * t146 + t226 * mrSges(5,1) + t227 * t194 + (-t185 - t186) * t203 + (-mrSges(5,3) - mrSges(6,3)) * t171 + t258;
t196 = t227 * mrSges(6,1) - t203 * mrSges(6,3);
t197 = t227 * mrSges(5,1) - t203 * mrSges(5,3);
t257 = m(6) * t141 + t170 * mrSges(6,3) + t202 * t185;
t132 = m(5) * t147 + t170 * mrSges(5,3) + t202 * t186 + (-t196 - t197) * t227 + (-mrSges(5,2) - mrSges(6,2)) * t226 + t257;
t125 = t236 * t129 + t233 * t132;
t200 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t234 + Ifges(4,2) * t237) * qJD(2);
t201 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t234 + Ifges(4,4) * t237) * qJD(2);
t264 = mrSges(4,1) * t173 - mrSges(4,2) * t174 + Ifges(4,5) * t211 + Ifges(4,6) * t212 + Ifges(4,3) * qJDD(3) + pkin(3) * t125 + (t234 * t200 - t237 * t201) * qJD(2) + t265;
t210 = (-mrSges(4,1) * t237 + mrSges(4,2) * t234) * qJD(2);
t260 = qJD(2) * t237;
t217 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t260;
t123 = m(4) * t173 + qJDD(3) * mrSges(4,1) - t211 * mrSges(4,3) + qJD(3) * t217 - t210 * t261 + t125;
t216 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t261;
t252 = -t233 * t129 + t236 * t132;
t124 = m(4) * t174 - qJDD(3) * mrSges(4,2) + t212 * mrSges(4,3) - qJD(3) * t216 + t210 * t260 + t252;
t253 = -t234 * t123 + t237 * t124;
t115 = m(3) * t192 - t239 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t253;
t191 = t238 * t214 - t235 * t215;
t249 = -qJDD(2) * pkin(2) - t191;
t187 = -t239 * pkin(6) + t249;
t153 = -t212 * pkin(3) + t218 * t261 + (-pkin(7) * t229 - pkin(6)) * t239 + t249;
t144 = -t170 * pkin(4) - t198 * qJ(5) + t203 * t195 + qJDD(5) + t153;
t251 = m(6) * t144 - t170 * mrSges(6,1) + t171 * mrSges(6,2) - t202 * t193 + t203 * t196;
t243 = m(5) * t153 - t170 * mrSges(5,1) + t171 * mrSges(5,2) - t202 * t194 + t203 * t197 + t251;
t241 = -m(4) * t187 + t212 * mrSges(4,1) - t211 * mrSges(4,2) - t216 * t261 + t217 * t260 - t243;
t127 = m(3) * t191 + qJDD(2) * mrSges(3,1) - t239 * mrSges(3,2) + t241;
t112 = t235 * t115 + t238 * t127;
t117 = t237 * t123 + t234 * t124;
t254 = t238 * t115 - t235 * t127;
t248 = -mrSges(6,1) * t144 + mrSges(6,3) * t141 + Ifges(6,4) * t171 + Ifges(6,2) * t170 + Ifges(6,6) * t226 + t227 * t179;
t175 = Ifges(6,5) * t203 + Ifges(6,6) * t202 + Ifges(6,3) * t227;
t246 = mrSges(6,2) * t144 - mrSges(6,3) * t138 + Ifges(6,1) * t171 + Ifges(6,4) * t170 + Ifges(6,5) * t226 + t202 * t175;
t176 = Ifges(5,5) * t203 + Ifges(5,6) * t202 + Ifges(5,3) * t227;
t118 = Ifges(5,4) * t171 + Ifges(5,2) * t170 + Ifges(5,6) * t226 + t227 * t180 - mrSges(5,1) * t153 + mrSges(5,3) * t147 - pkin(4) * t251 + qJ(5) * (-t226 * mrSges(6,2) - t227 * t196 + t257) + (-t176 - t175) * t203 + t248;
t119 = mrSges(5,2) * t153 - mrSges(5,3) * t146 + Ifges(5,1) * t171 + Ifges(5,4) * t170 + Ifges(5,5) * t226 - qJ(5) * t135 + t202 * t176 + (-t177 - t178) * t227 + t246;
t199 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t234 + Ifges(4,6) * t237) * qJD(2);
t106 = -mrSges(4,1) * t187 + mrSges(4,3) * t174 + Ifges(4,4) * t211 + Ifges(4,2) * t212 + Ifges(4,6) * qJDD(3) - pkin(3) * t243 + pkin(7) * t252 + qJD(3) * t201 + t236 * t118 + t233 * t119 - t199 * t261;
t108 = mrSges(4,2) * t187 - mrSges(4,3) * t173 + Ifges(4,1) * t211 + Ifges(4,4) * t212 + Ifges(4,5) * qJDD(3) - pkin(7) * t125 - qJD(3) * t200 - t233 * t118 + t236 * t119 + t199 * t260;
t245 = mrSges(3,1) * t191 - mrSges(3,2) * t192 + Ifges(3,3) * qJDD(2) + pkin(2) * t241 + pkin(6) * t253 + t237 * t106 + t234 * t108;
t244 = mrSges(2,1) * t214 - mrSges(2,2) * t215 + pkin(1) * t112 + t245;
t110 = m(2) * t215 + t254;
t109 = m(2) * t214 + t112;
t104 = -mrSges(3,1) * t230 + mrSges(3,3) * t192 + t239 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t117 - t264;
t103 = mrSges(3,2) * t230 - mrSges(3,3) * t191 + Ifges(3,5) * qJDD(2) - t239 * Ifges(3,6) - pkin(6) * t117 - t234 * t106 + t237 * t108;
t102 = mrSges(2,2) * t230 - mrSges(2,3) * t214 - pkin(5) * t112 + t238 * t103 - t235 * t104;
t101 = -mrSges(2,1) * t230 + mrSges(2,3) * t215 + t235 * t103 + t238 * t104 - pkin(1) * (m(3) * t230 + t117) + pkin(5) * t254;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t232 * t102 - t231 * t101 - qJ(1) * (t232 * t109 + t231 * t110), t102, t103, t108, t119, -t227 * t177 + t246; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t231 * t102 + t232 * t101 + qJ(1) * (-t231 * t109 + t232 * t110), t101, t104, t106, t118, -t203 * t175 + t248; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t244, t244, t245, t264, t265, -t202 * t179 - t247;];
m_new = t1;
