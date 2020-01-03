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
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:05:36
% EndTime: 2020-01-03 12:05:43
% DurationCPUTime: 5.50s
% Computational Cost: add. (92924->246), mult. (127422->323), div. (0->0), fcn. (78313->10), ass. (0->114)
t229 = sin(qJ(1));
t233 = cos(qJ(1));
t211 = -t233 * g(2) - t229 * g(3);
t201 = qJDD(1) * pkin(1) + t211;
t210 = -t229 * g(2) + t233 * g(3);
t234 = qJD(1) ^ 2;
t202 = -t234 * pkin(1) + t210;
t228 = sin(qJ(2));
t232 = cos(qJ(2));
t181 = t228 * t201 + t232 * t202;
t221 = (qJD(1) + qJD(2));
t218 = t221 ^ 2;
t219 = qJDD(1) + qJDD(2);
t267 = -t218 * pkin(2) + t219 * qJ(3) + (2 * qJD(3) * t221) + t181;
t224 = sin(pkin(9));
t262 = t221 * t224;
t225 = cos(pkin(9));
t260 = t225 * t221;
t167 = -t225 * g(1) - t267 * t224;
t168 = -t224 * g(1) + t267 * t225;
t249 = -pkin(3) * t225 - pkin(7) * t224;
t196 = t249 * t221;
t154 = t196 * t260 + t168;
t180 = t232 * t201 - t228 * t202;
t243 = -t218 * qJ(3) + qJDD(3) - t180;
t166 = (-pkin(2) + t249) * t219 + t243;
t231 = cos(qJ(4));
t165 = t231 * t166;
t227 = sin(qJ(4));
t257 = qJD(4) * t221;
t191 = (t219 * t231 - t227 * t257) * t224;
t261 = t225 * t219;
t206 = qJDD(4) - t261;
t207 = qJD(4) - t260;
t263 = t218 * t224 ^ 2;
t146 = t206 * pkin(4) - t191 * pkin(8) + t165 + (-pkin(4) * t231 * t263 - pkin(8) * t207 * t262 - t154) * t227;
t149 = t231 * t154 + t227 * t166;
t254 = t231 * t262;
t189 = t207 * pkin(4) - pkin(8) * t254;
t190 = (-t219 * t227 - t231 * t257) * t224;
t256 = t227 ^ 2 * t263;
t147 = -pkin(4) * t256 + t190 * pkin(8) - t207 * t189 + t149;
t226 = sin(qJ(5));
t230 = cos(qJ(5));
t145 = t226 * t146 + t230 * t147;
t153 = t196 * t262 - t167;
t150 = -t190 * pkin(4) - pkin(8) * t256 + t189 * t254 + t153;
t183 = (-t226 * t227 + t230 * t231) * t262;
t158 = -t183 * qJD(5) + t230 * t190 - t226 * t191;
t182 = (-t226 * t231 - t227 * t230) * t262;
t159 = t182 * qJD(5) + t226 * t190 + t230 * t191;
t205 = qJD(5) + t207;
t160 = Ifges(6,5) * t183 + Ifges(6,6) * t182 + Ifges(6,3) * t205;
t162 = Ifges(6,1) * t183 + Ifges(6,4) * t182 + Ifges(6,5) * t205;
t203 = qJDD(5) + t206;
t133 = -mrSges(6,1) * t150 + mrSges(6,3) * t145 + Ifges(6,4) * t159 + Ifges(6,2) * t158 + Ifges(6,6) * t203 - t183 * t160 + t205 * t162;
t144 = t230 * t146 - t226 * t147;
t161 = Ifges(6,4) * t183 + Ifges(6,2) * t182 + Ifges(6,6) * t205;
t134 = mrSges(6,2) * t150 - mrSges(6,3) * t144 + Ifges(6,1) * t159 + Ifges(6,4) * t158 + Ifges(6,5) * t203 + t182 * t160 - t205 * t161;
t176 = Ifges(5,3) * t207 + (Ifges(5,5) * t231 - Ifges(5,6) * t227) * t262;
t178 = Ifges(5,5) * t207 + (Ifges(5,1) * t231 - Ifges(5,4) * t227) * t262;
t174 = -t205 * mrSges(6,2) + t182 * mrSges(6,3);
t175 = t205 * mrSges(6,1) - t183 * mrSges(6,3);
t241 = m(6) * t150 - t158 * mrSges(6,1) + t159 * mrSges(6,2) - t182 * t174 + t183 * t175;
t169 = -t182 * mrSges(6,1) + t183 * mrSges(6,2);
t140 = m(6) * t144 + t203 * mrSges(6,1) - t159 * mrSges(6,3) - t183 * t169 + t205 * t174;
t141 = m(6) * t145 - t203 * mrSges(6,2) + t158 * mrSges(6,3) + t182 * t169 - t205 * t175;
t250 = -t226 * t140 + t230 * t141;
t117 = -mrSges(5,1) * t153 + mrSges(5,3) * t149 + Ifges(5,4) * t191 + Ifges(5,2) * t190 + Ifges(5,6) * t206 - pkin(4) * t241 + pkin(8) * t250 + t230 * t133 + t226 * t134 - t176 * t254 + t207 * t178;
t132 = t230 * t140 + t226 * t141;
t148 = -t227 * t154 + t165;
t177 = Ifges(5,6) * t207 + (Ifges(5,4) * t231 - Ifges(5,2) * t227) * t262;
t255 = t227 * t262;
t120 = mrSges(5,2) * t153 - mrSges(5,3) * t148 + Ifges(5,1) * t191 + Ifges(5,4) * t190 + Ifges(5,5) * t206 - pkin(8) * t132 - t226 * t133 + t230 * t134 - t176 * t255 - t207 * t177;
t186 = -t207 * mrSges(5,2) - mrSges(5,3) * t255;
t188 = (mrSges(5,1) * t227 + mrSges(5,2) * t231) * t262;
t130 = m(5) * t148 + t206 * mrSges(5,1) - t191 * mrSges(5,3) + t207 * t186 - t188 * t254 + t132;
t187 = t207 * mrSges(5,1) - mrSges(5,3) * t254;
t131 = m(5) * t149 - t206 * mrSges(5,2) + t190 * mrSges(5,3) - t207 * t187 - t188 * t255 + t250;
t128 = -t227 * t130 + t231 * t131;
t236 = -m(5) * t153 + t190 * mrSges(5,1) - t191 * mrSges(5,2) - t241;
t245 = -t186 * t227 - t187 * t231;
t248 = Ifges(4,1) * t224 + Ifges(4,4) * t225;
t266 = -((Ifges(4,4) * t224 + Ifges(4,2) * t225) * t262 - t248 * t260) * t221 - mrSges(4,1) * t167 + mrSges(4,2) * t168 - pkin(3) * (t245 * t262 + t236) - pkin(7) * t128 - t231 * t117 - t227 * t120;
t265 = mrSges(4,2) * t224;
t264 = mrSges(4,3) * t219;
t192 = (-mrSges(4,1) * t225 + t265) * t221;
t125 = m(4) * t168 + (t192 * t221 + t264) * t225 + t128;
t136 = m(4) * t167 + (-t264 + (-t192 + t245) * t221) * t224 + t236;
t251 = t225 * t125 - t224 * t136;
t116 = m(3) * t181 - t218 * mrSges(3,1) - t219 * mrSges(3,2) + t251;
t127 = t231 * t130 + t227 * t131;
t172 = -t219 * pkin(2) + t243;
t239 = -m(4) * t172 + mrSges(4,1) * t261 - t127 + (t218 * t225 ^ 2 + t263) * mrSges(4,3);
t122 = m(3) * t180 - t218 * mrSges(3,2) + (mrSges(3,1) - t265) * t219 + t239;
t111 = t228 * t116 + t232 * t122;
t119 = t224 * t125 + t225 * t136;
t252 = t232 * t116 - t228 * t122;
t247 = Ifges(4,5) * t224 + Ifges(4,6) * t225;
t246 = t177 * t231 + t178 * t227;
t193 = t247 * t221;
t107 = mrSges(4,2) * t172 - mrSges(4,3) * t167 - pkin(7) * t127 - t227 * t117 + t231 * t120 + t193 * t260 + t248 * t219;
t240 = -mrSges(6,1) * t144 + mrSges(6,2) * t145 - Ifges(6,5) * t159 - Ifges(6,6) * t158 - Ifges(6,3) * t203 - t183 * t161 + t182 * t162;
t235 = mrSges(5,1) * t148 - mrSges(5,2) * t149 + Ifges(5,5) * t191 + Ifges(5,6) * t190 + Ifges(5,3) * t206 + pkin(4) * t132 - t240;
t113 = Ifges(4,2) * t261 - t235 + (Ifges(4,4) * t219 + (-t193 - t246) * t221) * t224 + mrSges(4,3) * t168 - mrSges(4,1) * t172 - pkin(3) * t127;
t242 = -mrSges(3,2) * t181 + qJ(3) * t251 + t224 * t107 + t225 * t113 + pkin(2) * (-t219 * t265 + t239) + mrSges(3,1) * t180 + Ifges(3,3) * t219;
t238 = mrSges(2,1) * t211 - mrSges(2,2) * t210 + Ifges(2,3) * qJDD(1) + pkin(1) * t111 + t242;
t109 = m(2) * t211 + qJDD(1) * mrSges(2,1) - t234 * mrSges(2,2) + t111;
t108 = m(2) * t210 - t234 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t252;
t105 = mrSges(3,1) * g(1) + mrSges(3,3) * t181 + t218 * Ifges(3,5) - pkin(2) * t119 + (Ifges(3,6) - t247) * t219 + t266;
t104 = -mrSges(3,2) * g(1) - mrSges(3,3) * t180 + Ifges(3,5) * t219 - t218 * Ifges(3,6) - qJ(3) * t119 + t225 * t107 - t224 * t113;
t103 = -mrSges(2,2) * g(1) - mrSges(2,3) * t211 + Ifges(2,5) * qJDD(1) - t234 * Ifges(2,6) - pkin(6) * t111 + t232 * t104 - t228 * t105;
t102 = Ifges(2,6) * qJDD(1) + t234 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t210 + t228 * t104 + t232 * t105 - pkin(1) * (-m(3) * g(1) + t119) + pkin(6) * t252;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t238, t103, t104, t107, t120, t134; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t229 * t103 + t233 * t102 - pkin(5) * (-t233 * t108 + t229 * t109), t102, t105, t113, t117, t133; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t233 * t103 + t229 * t102 + pkin(5) * (t229 * t108 + t233 * t109), t238, t242, t247 * t219 - t266, t246 * t262 + t235, -t240;];
m_new = t1;
