% Calculate vector of cutting torques with Newton-Euler for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:28
% EndTime: 2020-01-03 11:38:43
% DurationCPUTime: 7.31s
% Computational Cost: add. (91738->267), mult. (197708->343), div. (0->0), fcn. (125769->10), ass. (0->108)
t239 = sin(qJ(1));
t242 = cos(qJ(1));
t221 = -t242 * g(2) - t239 * g(3);
t211 = qJDD(1) * pkin(1) + t221;
t220 = -t239 * g(2) + t242 * g(3);
t243 = qJD(1) ^ 2;
t213 = -t243 * pkin(1) + t220;
t234 = sin(pkin(8));
t236 = cos(pkin(8));
t191 = t234 * t211 + t236 * t213;
t183 = -t243 * pkin(2) + qJDD(1) * pkin(6) + t191;
t232 = -g(1) + qJDD(2);
t238 = sin(qJ(3));
t241 = cos(qJ(3));
t172 = -t238 * t183 + t241 * t232;
t259 = qJD(1) * qJD(3);
t258 = t241 * t259;
t214 = t238 * qJDD(1) + t258;
t168 = (-t214 + t258) * qJ(4) + (t238 * t241 * t243 + qJDD(3)) * pkin(3) + t172;
t173 = t241 * t183 + t238 * t232;
t215 = t241 * qJDD(1) - t238 * t259;
t261 = qJD(1) * t238;
t217 = qJD(3) * pkin(3) - qJ(4) * t261;
t231 = t241 ^ 2;
t169 = -t231 * t243 * pkin(3) + t215 * qJ(4) - qJD(3) * t217 + t173;
t233 = sin(pkin(9));
t235 = cos(pkin(9));
t201 = (t233 * t241 + t235 * t238) * qJD(1);
t148 = -0.2e1 * qJD(4) * t201 + t235 * t168 - t233 * t169;
t193 = t235 * t214 + t233 * t215;
t200 = (-t233 * t238 + t235 * t241) * qJD(1);
t145 = (qJD(3) * t200 - t193) * pkin(7) + (t200 * t201 + qJDD(3)) * pkin(4) + t148;
t149 = 0.2e1 * qJD(4) * t200 + t233 * t168 + t235 * t169;
t192 = -t233 * t214 + t235 * t215;
t196 = qJD(3) * pkin(4) - t201 * pkin(7);
t199 = t200 ^ 2;
t146 = -t199 * pkin(4) + t192 * pkin(7) - qJD(3) * t196 + t149;
t237 = sin(qJ(5));
t240 = cos(qJ(5));
t143 = t240 * t145 - t237 * t146;
t180 = t240 * t200 - t237 * t201;
t158 = t180 * qJD(5) + t237 * t192 + t240 * t193;
t181 = t237 * t200 + t240 * t201;
t167 = -t180 * mrSges(6,1) + t181 * mrSges(6,2);
t227 = qJD(3) + qJD(5);
t174 = -t227 * mrSges(6,2) + t180 * mrSges(6,3);
t226 = qJDD(3) + qJDD(5);
t140 = m(6) * t143 + t226 * mrSges(6,1) - t158 * mrSges(6,3) - t181 * t167 + t227 * t174;
t144 = t237 * t145 + t240 * t146;
t157 = -t181 * qJD(5) + t240 * t192 - t237 * t193;
t175 = t227 * mrSges(6,1) - t181 * mrSges(6,3);
t141 = m(6) * t144 - t226 * mrSges(6,2) + t157 * mrSges(6,3) + t180 * t167 - t227 * t175;
t131 = t240 * t140 + t237 * t141;
t185 = -t200 * mrSges(5,1) + t201 * mrSges(5,2);
t194 = -qJD(3) * mrSges(5,2) + t200 * mrSges(5,3);
t128 = m(5) * t148 + qJDD(3) * mrSges(5,1) - t193 * mrSges(5,3) + qJD(3) * t194 - t201 * t185 + t131;
t195 = qJD(3) * mrSges(5,1) - t201 * mrSges(5,3);
t254 = -t237 * t140 + t240 * t141;
t129 = m(5) * t149 - qJDD(3) * mrSges(5,2) + t192 * mrSges(5,3) - qJD(3) * t195 + t200 * t185 + t254;
t124 = t235 * t128 + t233 * t129;
t206 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t238 + Ifges(4,2) * t241) * qJD(1);
t207 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t238 + Ifges(4,4) * t241) * qJD(1);
t178 = Ifges(5,4) * t201 + Ifges(5,2) * t200 + Ifges(5,6) * qJD(3);
t179 = Ifges(5,1) * t201 + Ifges(5,4) * t200 + Ifges(5,5) * qJD(3);
t160 = Ifges(6,4) * t181 + Ifges(6,2) * t180 + Ifges(6,6) * t227;
t161 = Ifges(6,1) * t181 + Ifges(6,4) * t180 + Ifges(6,5) * t227;
t249 = -mrSges(6,1) * t143 + mrSges(6,2) * t144 - Ifges(6,5) * t158 - Ifges(6,6) * t157 - Ifges(6,3) * t226 - t181 * t160 + t180 * t161;
t246 = -mrSges(5,1) * t148 + mrSges(5,2) * t149 - Ifges(5,5) * t193 - Ifges(5,6) * t192 - Ifges(5,3) * qJDD(3) - pkin(4) * t131 - t201 * t178 + t200 * t179 + t249;
t262 = mrSges(4,1) * t172 - mrSges(4,2) * t173 + Ifges(4,5) * t214 + Ifges(4,6) * t215 + Ifges(4,3) * qJDD(3) + pkin(3) * t124 + (t238 * t206 - t241 * t207) * qJD(1) - t246;
t212 = (-mrSges(4,1) * t241 + mrSges(4,2) * t238) * qJD(1);
t260 = qJD(1) * t241;
t219 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t260;
t122 = m(4) * t172 + qJDD(3) * mrSges(4,1) - t214 * mrSges(4,3) + qJD(3) * t219 - t212 * t261 + t124;
t218 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t261;
t255 = -t233 * t128 + t235 * t129;
t123 = m(4) * t173 - qJDD(3) * mrSges(4,2) + t215 * mrSges(4,3) - qJD(3) * t218 + t212 * t260 + t255;
t256 = -t238 * t122 + t241 * t123;
t114 = m(3) * t191 - t243 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t256;
t190 = t236 * t211 - t234 * t213;
t251 = -qJDD(1) * pkin(2) - t190;
t182 = -t243 * pkin(6) + t251;
t170 = -t215 * pkin(3) + qJDD(4) + t217 * t261 + (-qJ(4) * t231 - pkin(6)) * t243 + t251;
t151 = -t192 * pkin(4) - t199 * pkin(7) + t201 * t196 + t170;
t253 = m(6) * t151 - t157 * mrSges(6,1) + t158 * mrSges(6,2) - t180 * t174 + t181 * t175;
t247 = m(5) * t170 - t192 * mrSges(5,1) + t193 * mrSges(5,2) - t200 * t194 + t201 * t195 + t253;
t245 = -m(4) * t182 + t215 * mrSges(4,1) - t214 * mrSges(4,2) - t218 * t261 + t219 * t260 - t247;
t135 = m(3) * t190 + qJDD(1) * mrSges(3,1) - t243 * mrSges(3,2) + t245;
t111 = t234 * t114 + t236 * t135;
t116 = t241 * t122 + t238 * t123;
t257 = t236 * t114 - t234 * t135;
t159 = Ifges(6,5) * t181 + Ifges(6,6) * t180 + Ifges(6,3) * t227;
t132 = -mrSges(6,1) * t151 + mrSges(6,3) * t144 + Ifges(6,4) * t158 + Ifges(6,2) * t157 + Ifges(6,6) * t226 - t181 * t159 + t227 * t161;
t133 = mrSges(6,2) * t151 - mrSges(6,3) * t143 + Ifges(6,1) * t158 + Ifges(6,4) * t157 + Ifges(6,5) * t226 + t180 * t159 - t227 * t160;
t177 = Ifges(5,5) * t201 + Ifges(5,6) * t200 + Ifges(5,3) * qJD(3);
t117 = -mrSges(5,1) * t170 + mrSges(5,3) * t149 + Ifges(5,4) * t193 + Ifges(5,2) * t192 + Ifges(5,6) * qJDD(3) - pkin(4) * t253 + pkin(7) * t254 + qJD(3) * t179 + t240 * t132 + t237 * t133 - t201 * t177;
t118 = mrSges(5,2) * t170 - mrSges(5,3) * t148 + Ifges(5,1) * t193 + Ifges(5,4) * t192 + Ifges(5,5) * qJDD(3) - pkin(7) * t131 - qJD(3) * t178 - t237 * t132 + t240 * t133 + t200 * t177;
t205 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t238 + Ifges(4,6) * t241) * qJD(1);
t105 = -mrSges(4,1) * t182 + mrSges(4,3) * t173 + Ifges(4,4) * t214 + Ifges(4,2) * t215 + Ifges(4,6) * qJDD(3) - pkin(3) * t247 + qJ(4) * t255 + qJD(3) * t207 + t235 * t117 + t233 * t118 - t205 * t261;
t107 = mrSges(4,2) * t182 - mrSges(4,3) * t172 + Ifges(4,1) * t214 + Ifges(4,4) * t215 + Ifges(4,5) * qJDD(3) - qJ(4) * t124 - qJD(3) * t206 - t233 * t117 + t235 * t118 + t205 * t260;
t250 = mrSges(3,1) * t190 - mrSges(3,2) * t191 + Ifges(3,3) * qJDD(1) + pkin(2) * t245 + pkin(6) * t256 + t241 * t105 + t238 * t107;
t248 = mrSges(2,1) * t221 - mrSges(2,2) * t220 + Ifges(2,3) * qJDD(1) + pkin(1) * t111 + t250;
t109 = m(2) * t221 + qJDD(1) * mrSges(2,1) - t243 * mrSges(2,2) + t111;
t108 = m(2) * t220 - t243 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t257;
t103 = -mrSges(3,1) * t232 + mrSges(3,3) * t191 + t243 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t116 - t262;
t102 = mrSges(3,2) * t232 - mrSges(3,3) * t190 + Ifges(3,5) * qJDD(1) - t243 * Ifges(3,6) - pkin(6) * t116 - t238 * t105 + t241 * t107;
t101 = -mrSges(2,2) * g(1) - mrSges(2,3) * t221 + Ifges(2,5) * qJDD(1) - t243 * Ifges(2,6) - qJ(2) * t111 + t236 * t102 - t234 * t103;
t100 = Ifges(2,6) * qJDD(1) + t243 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t220 + t234 * t102 + t236 * t103 - pkin(1) * (m(3) * t232 + t116) + qJ(2) * t257;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t248, t101, t102, t107, t118, t133; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t239 * t101 + t242 * t100 - pkin(5) * (-t242 * t108 + t239 * t109), t100, t103, t105, t117, t132; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t242 * t101 + t239 * t100 + pkin(5) * (t239 * t108 + t242 * t109), t248, t250, t262, -t246, -t249;];
m_new = t1;
