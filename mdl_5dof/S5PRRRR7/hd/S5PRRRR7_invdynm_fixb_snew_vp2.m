% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR7
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:56
% EndTime: 2019-12-05 17:12:09
% DurationCPUTime: 7.00s
% Computational Cost: add. (90309->258), mult. (187520->330), div. (0->0), fcn. (127537->10), ass. (0->108)
t224 = sin(pkin(9));
t250 = cos(pkin(9));
t210 = -t250 * g(1) - t224 * g(2);
t223 = -g(3) + qJDD(1);
t228 = sin(qJ(2));
t232 = cos(qJ(2));
t191 = t232 * t210 + t228 * t223;
t233 = qJD(2) ^ 2;
t186 = -t233 * pkin(2) + qJDD(2) * pkin(6) + t191;
t209 = t224 * g(1) - t250 * g(2);
t227 = sin(qJ(3));
t231 = cos(qJ(3));
t175 = -t227 * t186 - t231 * t209;
t247 = qJD(2) * qJD(3);
t246 = t231 * t247;
t206 = t227 * qJDD(2) + t246;
t164 = (-t206 + t246) * pkin(7) + (t227 * t231 * t233 + qJDD(3)) * pkin(3) + t175;
t176 = t231 * t186 - t227 * t209;
t207 = t231 * qJDD(2) - t227 * t247;
t249 = qJD(2) * t227;
t213 = qJD(3) * pkin(3) - pkin(7) * t249;
t222 = t231 ^ 2;
t165 = -t222 * t233 * pkin(3) + t207 * pkin(7) - qJD(3) * t213 + t176;
t226 = sin(qJ(4));
t230 = cos(qJ(4));
t146 = t230 * t164 - t226 * t165;
t196 = (-t226 * t227 + t230 * t231) * qJD(2);
t172 = t196 * qJD(4) + t230 * t206 + t226 * t207;
t197 = (t226 * t231 + t227 * t230) * qJD(2);
t220 = qJDD(3) + qJDD(4);
t221 = qJD(3) + qJD(4);
t141 = (t196 * t221 - t172) * pkin(8) + (t196 * t197 + t220) * pkin(4) + t146;
t147 = t226 * t164 + t230 * t165;
t171 = -t197 * qJD(4) - t226 * t206 + t230 * t207;
t189 = t221 * pkin(4) - t197 * pkin(8);
t192 = t196 ^ 2;
t142 = -t192 * pkin(4) + t171 * pkin(8) - t221 * t189 + t147;
t225 = sin(qJ(5));
t229 = cos(qJ(5));
t139 = t229 * t141 - t225 * t142;
t181 = t229 * t196 - t225 * t197;
t153 = t181 * qJD(5) + t225 * t171 + t229 * t172;
t182 = t225 * t196 + t229 * t197;
t160 = -t181 * mrSges(6,1) + t182 * mrSges(6,2);
t218 = qJD(5) + t221;
t173 = -t218 * mrSges(6,2) + t181 * mrSges(6,3);
t217 = qJDD(5) + t220;
t136 = m(6) * t139 + t217 * mrSges(6,1) - t153 * mrSges(6,3) - t182 * t160 + t218 * t173;
t140 = t225 * t141 + t229 * t142;
t152 = -t182 * qJD(5) + t229 * t171 - t225 * t172;
t174 = t218 * mrSges(6,1) - t182 * mrSges(6,3);
t137 = m(6) * t140 - t217 * mrSges(6,2) + t152 * mrSges(6,3) + t181 * t160 - t218 * t174;
t128 = t229 * t136 + t225 * t137;
t183 = -t196 * mrSges(5,1) + t197 * mrSges(5,2);
t187 = -t221 * mrSges(5,2) + t196 * mrSges(5,3);
t125 = m(5) * t146 + t220 * mrSges(5,1) - t172 * mrSges(5,3) - t197 * t183 + t221 * t187 + t128;
t188 = t221 * mrSges(5,1) - t197 * mrSges(5,3);
t243 = -t225 * t136 + t229 * t137;
t126 = m(5) * t147 - t220 * mrSges(5,2) + t171 * mrSges(5,3) + t196 * t183 - t221 * t188 + t243;
t121 = t230 * t125 + t226 * t126;
t194 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t227 + Ifges(4,2) * t231) * qJD(2);
t195 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t227 + Ifges(4,4) * t231) * qJD(2);
t178 = Ifges(5,4) * t197 + Ifges(5,2) * t196 + Ifges(5,6) * t221;
t179 = Ifges(5,1) * t197 + Ifges(5,4) * t196 + Ifges(5,5) * t221;
t156 = Ifges(6,4) * t182 + Ifges(6,2) * t181 + Ifges(6,6) * t218;
t157 = Ifges(6,1) * t182 + Ifges(6,4) * t181 + Ifges(6,5) * t218;
t238 = -mrSges(6,1) * t139 + mrSges(6,2) * t140 - Ifges(6,5) * t153 - Ifges(6,6) * t152 - Ifges(6,3) * t217 - t182 * t156 + t181 * t157;
t235 = -mrSges(5,1) * t146 + mrSges(5,2) * t147 - Ifges(5,5) * t172 - Ifges(5,6) * t171 - Ifges(5,3) * t220 - pkin(4) * t128 - t197 * t178 + t196 * t179 + t238;
t251 = mrSges(4,1) * t175 - mrSges(4,2) * t176 + Ifges(4,5) * t206 + Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t121 + (t227 * t194 - t231 * t195) * qJD(2) - t235;
t248 = qJD(2) * t231;
t205 = (-mrSges(4,1) * t231 + mrSges(4,2) * t227) * qJD(2);
t212 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t248;
t119 = m(4) * t175 + qJDD(3) * mrSges(4,1) - t206 * mrSges(4,3) + qJD(3) * t212 - t205 * t249 + t121;
t211 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t249;
t244 = -t226 * t125 + t230 * t126;
t120 = m(4) * t176 - qJDD(3) * mrSges(4,2) + t207 * mrSges(4,3) - qJD(3) * t211 + t205 * t248 + t244;
t115 = -t227 * t119 + t231 * t120;
t111 = m(3) * t191 - t233 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t115;
t190 = -t228 * t210 + t232 * t223;
t240 = -qJDD(2) * pkin(2) - t190;
t185 = -t233 * pkin(6) + t240;
t166 = -t207 * pkin(3) + t213 * t249 + (-pkin(7) * t222 - pkin(6)) * t233 + t240;
t144 = -t171 * pkin(4) - t192 * pkin(8) + t197 * t189 + t166;
t242 = m(6) * t144 - t152 * mrSges(6,1) + t153 * mrSges(6,2) - t181 * t173 + t182 * t174;
t237 = m(5) * t166 - t171 * mrSges(5,1) + t172 * mrSges(5,2) - t196 * t187 + t197 * t188 + t242;
t132 = -m(4) * t185 + t207 * mrSges(4,1) - t206 * mrSges(4,2) - t211 * t249 + t212 * t248 - t237;
t131 = m(3) * t190 + qJDD(2) * mrSges(3,1) - t233 * mrSges(3,2) + t132;
t245 = t232 * t111 - t228 * t131;
t114 = t231 * t119 + t227 * t120;
t155 = Ifges(6,5) * t182 + Ifges(6,6) * t181 + Ifges(6,3) * t218;
t129 = -mrSges(6,1) * t144 + mrSges(6,3) * t140 + Ifges(6,4) * t153 + Ifges(6,2) * t152 + Ifges(6,6) * t217 - t182 * t155 + t218 * t157;
t130 = mrSges(6,2) * t144 - mrSges(6,3) * t139 + Ifges(6,1) * t153 + Ifges(6,4) * t152 + Ifges(6,5) * t217 + t181 * t155 - t218 * t156;
t177 = Ifges(5,5) * t197 + Ifges(5,6) * t196 + Ifges(5,3) * t221;
t116 = -mrSges(5,1) * t166 + mrSges(5,3) * t147 + Ifges(5,4) * t172 + Ifges(5,2) * t171 + Ifges(5,6) * t220 - pkin(4) * t242 + pkin(8) * t243 + t229 * t129 + t225 * t130 - t197 * t177 + t221 * t179;
t117 = mrSges(5,2) * t166 - mrSges(5,3) * t146 + Ifges(5,1) * t172 + Ifges(5,4) * t171 + Ifges(5,5) * t220 - pkin(8) * t128 - t225 * t129 + t229 * t130 + t196 * t177 - t221 * t178;
t193 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t227 + Ifges(4,6) * t231) * qJD(2);
t105 = -mrSges(4,1) * t185 + mrSges(4,3) * t176 + Ifges(4,4) * t206 + Ifges(4,2) * t207 + Ifges(4,6) * qJDD(3) - pkin(3) * t237 + pkin(7) * t244 + qJD(3) * t195 + t230 * t116 + t226 * t117 - t193 * t249;
t106 = mrSges(4,2) * t185 - mrSges(4,3) * t175 + Ifges(4,1) * t206 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - pkin(7) * t121 - qJD(3) * t194 - t226 * t116 + t230 * t117 + t193 * t248;
t102 = -mrSges(3,2) * t209 - mrSges(3,3) * t190 + Ifges(3,5) * qJDD(2) - t233 * Ifges(3,6) - pkin(6) * t114 - t227 * t105 + t231 * t106;
t104 = mrSges(3,1) * t209 + mrSges(3,3) * t191 + t233 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t114 - t251;
t239 = -mrSges(2,2) * t210 + pkin(5) * t245 + t228 * t102 + t232 * t104 + pkin(1) * (m(3) * t209 - t114) + mrSges(2,1) * t209;
t236 = mrSges(3,1) * t190 - mrSges(3,2) * t191 + Ifges(3,3) * qJDD(2) + pkin(2) * t132 + pkin(6) * t115 + t231 * t105 + t227 * t106;
t112 = (m(2) + m(3)) * t209 - t114;
t109 = t228 * t111 + t232 * t131;
t107 = m(2) * t210 + t245;
t100 = -mrSges(2,1) * t223 + mrSges(2,3) * t210 - pkin(1) * t109 - t236;
t99 = mrSges(2,2) * t223 - mrSges(2,3) * t209 - pkin(5) * t109 + t232 * t102 - t228 * t104;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t250 * t99 - t224 * t100 - qJ(1) * (t224 * t107 + t250 * t112), t99, t102, t106, t117, t130; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t224 * t99 + t250 * t100 + qJ(1) * (t250 * t107 - t224 * t112), t100, t104, t105, t116, t129; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t239, t239, t236, t251, -t235, -t238;];
m_new = t1;
