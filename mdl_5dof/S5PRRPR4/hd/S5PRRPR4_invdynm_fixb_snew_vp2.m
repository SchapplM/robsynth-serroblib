% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPR4
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:39
% EndTime: 2019-12-05 16:21:54
% DurationCPUTime: 6.83s
% Computational Cost: add. (81853->257), mult. (180062->332), div. (0->0), fcn. (120457->10), ass. (0->106)
t225 = sin(pkin(8));
t250 = cos(pkin(8));
t212 = -t250 * g(1) - t225 * g(2);
t223 = -g(3) + qJDD(1);
t229 = sin(qJ(2));
t232 = cos(qJ(2));
t191 = t232 * t212 + t229 * t223;
t233 = qJD(2) ^ 2;
t186 = -t233 * pkin(2) + qJDD(2) * pkin(6) + t191;
t211 = t225 * g(1) - t250 * g(2);
t228 = sin(qJ(3));
t231 = cos(qJ(3));
t170 = -t228 * t186 - t231 * t211;
t247 = qJD(2) * qJD(3);
t246 = t231 * t247;
t208 = t228 * qJDD(2) + t246;
t164 = (-t208 + t246) * qJ(4) + (t228 * t231 * t233 + qJDD(3)) * pkin(3) + t170;
t171 = t231 * t186 - t228 * t211;
t209 = t231 * qJDD(2) - t228 * t247;
t249 = qJD(2) * t228;
t213 = qJD(3) * pkin(3) - qJ(4) * t249;
t222 = t231 ^ 2;
t165 = -t222 * t233 * pkin(3) + t209 * qJ(4) - qJD(3) * t213 + t171;
t224 = sin(pkin(9));
t226 = cos(pkin(9));
t196 = (t224 * t231 + t226 * t228) * qJD(2);
t144 = -0.2e1 * qJD(4) * t196 + t226 * t164 - t224 * t165;
t183 = t226 * t208 + t224 * t209;
t195 = (-t224 * t228 + t226 * t231) * qJD(2);
t141 = (qJD(3) * t195 - t183) * pkin(7) + (t195 * t196 + qJDD(3)) * pkin(4) + t144;
t145 = 0.2e1 * qJD(4) * t195 + t224 * t164 + t226 * t165;
t182 = -t224 * t208 + t226 * t209;
t189 = qJD(3) * pkin(4) - t196 * pkin(7);
t194 = t195 ^ 2;
t142 = -t194 * pkin(4) + t182 * pkin(7) - qJD(3) * t189 + t145;
t227 = sin(qJ(5));
t230 = cos(qJ(5));
t139 = t230 * t141 - t227 * t142;
t175 = t230 * t195 - t227 * t196;
t154 = t175 * qJD(5) + t227 * t182 + t230 * t183;
t176 = t227 * t195 + t230 * t196;
t160 = -t175 * mrSges(6,1) + t176 * mrSges(6,2);
t220 = qJD(3) + qJD(5);
t168 = -t220 * mrSges(6,2) + t175 * mrSges(6,3);
t219 = qJDD(3) + qJDD(5);
t136 = m(6) * t139 + t219 * mrSges(6,1) - t154 * mrSges(6,3) - t176 * t160 + t220 * t168;
t140 = t227 * t141 + t230 * t142;
t153 = -t176 * qJD(5) + t230 * t182 - t227 * t183;
t169 = t220 * mrSges(6,1) - t176 * mrSges(6,3);
t137 = m(6) * t140 - t219 * mrSges(6,2) + t153 * mrSges(6,3) + t175 * t160 - t220 * t169;
t128 = t230 * t136 + t227 * t137;
t178 = -t195 * mrSges(5,1) + t196 * mrSges(5,2);
t187 = -qJD(3) * mrSges(5,2) + t195 * mrSges(5,3);
t125 = m(5) * t144 + qJDD(3) * mrSges(5,1) - t183 * mrSges(5,3) + qJD(3) * t187 - t196 * t178 + t128;
t188 = qJD(3) * mrSges(5,1) - t196 * mrSges(5,3);
t243 = -t227 * t136 + t230 * t137;
t126 = m(5) * t145 - qJDD(3) * mrSges(5,2) + t182 * mrSges(5,3) - qJD(3) * t188 + t195 * t178 + t243;
t121 = t226 * t125 + t224 * t126;
t198 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t228 + Ifges(4,2) * t231) * qJD(2);
t199 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t228 + Ifges(4,4) * t231) * qJD(2);
t173 = Ifges(5,4) * t196 + Ifges(5,2) * t195 + Ifges(5,6) * qJD(3);
t174 = Ifges(5,1) * t196 + Ifges(5,4) * t195 + Ifges(5,5) * qJD(3);
t156 = Ifges(6,4) * t176 + Ifges(6,2) * t175 + Ifges(6,6) * t220;
t157 = Ifges(6,1) * t176 + Ifges(6,4) * t175 + Ifges(6,5) * t220;
t238 = -mrSges(6,1) * t139 + mrSges(6,2) * t140 - Ifges(6,5) * t154 - Ifges(6,6) * t153 - Ifges(6,3) * t219 - t176 * t156 + t175 * t157;
t235 = -mrSges(5,1) * t144 + mrSges(5,2) * t145 - Ifges(5,5) * t183 - Ifges(5,6) * t182 - Ifges(5,3) * qJDD(3) - pkin(4) * t128 - t196 * t173 + t195 * t174 + t238;
t251 = mrSges(4,1) * t170 - mrSges(4,2) * t171 + Ifges(4,5) * t208 + Ifges(4,6) * t209 + Ifges(4,3) * qJDD(3) + pkin(3) * t121 + (t228 * t198 - t231 * t199) * qJD(2) - t235;
t248 = qJD(2) * t231;
t207 = (-mrSges(4,1) * t231 + mrSges(4,2) * t228) * qJD(2);
t215 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t248;
t119 = m(4) * t170 + qJDD(3) * mrSges(4,1) - t208 * mrSges(4,3) + qJD(3) * t215 - t207 * t249 + t121;
t214 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t249;
t244 = -t224 * t125 + t226 * t126;
t120 = m(4) * t171 - qJDD(3) * mrSges(4,2) + t209 * mrSges(4,3) - qJD(3) * t214 + t207 * t248 + t244;
t115 = -t228 * t119 + t231 * t120;
t111 = m(3) * t191 - t233 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t115;
t190 = -t229 * t212 + t232 * t223;
t240 = -qJDD(2) * pkin(2) - t190;
t185 = -t233 * pkin(6) + t240;
t166 = -t209 * pkin(3) + qJDD(4) + t213 * t249 + (-qJ(4) * t222 - pkin(6)) * t233 + t240;
t147 = -t182 * pkin(4) - t194 * pkin(7) + t196 * t189 + t166;
t242 = m(6) * t147 - t153 * mrSges(6,1) + t154 * mrSges(6,2) - t175 * t168 + t176 * t169;
t237 = m(5) * t166 - t182 * mrSges(5,1) + t183 * mrSges(5,2) - t195 * t187 + t196 * t188 + t242;
t132 = -m(4) * t185 + t209 * mrSges(4,1) - t208 * mrSges(4,2) - t214 * t249 + t215 * t248 - t237;
t131 = m(3) * t190 + qJDD(2) * mrSges(3,1) - t233 * mrSges(3,2) + t132;
t245 = t232 * t111 - t229 * t131;
t114 = t231 * t119 + t228 * t120;
t155 = Ifges(6,5) * t176 + Ifges(6,6) * t175 + Ifges(6,3) * t220;
t129 = -mrSges(6,1) * t147 + mrSges(6,3) * t140 + Ifges(6,4) * t154 + Ifges(6,2) * t153 + Ifges(6,6) * t219 - t176 * t155 + t220 * t157;
t130 = mrSges(6,2) * t147 - mrSges(6,3) * t139 + Ifges(6,1) * t154 + Ifges(6,4) * t153 + Ifges(6,5) * t219 + t175 * t155 - t220 * t156;
t172 = Ifges(5,5) * t196 + Ifges(5,6) * t195 + Ifges(5,3) * qJD(3);
t116 = -mrSges(5,1) * t166 + mrSges(5,3) * t145 + Ifges(5,4) * t183 + Ifges(5,2) * t182 + Ifges(5,6) * qJDD(3) - pkin(4) * t242 + pkin(7) * t243 + qJD(3) * t174 + t230 * t129 + t227 * t130 - t196 * t172;
t117 = mrSges(5,2) * t166 - mrSges(5,3) * t144 + Ifges(5,1) * t183 + Ifges(5,4) * t182 + Ifges(5,5) * qJDD(3) - pkin(7) * t128 - qJD(3) * t173 - t227 * t129 + t230 * t130 + t195 * t172;
t197 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t228 + Ifges(4,6) * t231) * qJD(2);
t105 = -mrSges(4,1) * t185 + mrSges(4,3) * t171 + Ifges(4,4) * t208 + Ifges(4,2) * t209 + Ifges(4,6) * qJDD(3) - pkin(3) * t237 + qJ(4) * t244 + qJD(3) * t199 + t226 * t116 + t224 * t117 - t197 * t249;
t106 = mrSges(4,2) * t185 - mrSges(4,3) * t170 + Ifges(4,1) * t208 + Ifges(4,4) * t209 + Ifges(4,5) * qJDD(3) - qJ(4) * t121 - qJD(3) * t198 - t224 * t116 + t226 * t117 + t197 * t248;
t102 = -mrSges(3,2) * t211 - mrSges(3,3) * t190 + Ifges(3,5) * qJDD(2) - t233 * Ifges(3,6) - pkin(6) * t114 - t228 * t105 + t231 * t106;
t104 = mrSges(3,1) * t211 + mrSges(3,3) * t191 + t233 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t114 - t251;
t239 = -mrSges(2,2) * t212 + pkin(5) * t245 + t229 * t102 + t232 * t104 + pkin(1) * (m(3) * t211 - t114) + mrSges(2,1) * t211;
t236 = mrSges(3,1) * t190 - mrSges(3,2) * t191 + Ifges(3,3) * qJDD(2) + pkin(2) * t132 + pkin(6) * t115 + t231 * t105 + t228 * t106;
t112 = (m(2) + m(3)) * t211 - t114;
t109 = t229 * t111 + t232 * t131;
t107 = m(2) * t212 + t245;
t100 = -mrSges(2,1) * t223 + mrSges(2,3) * t212 - pkin(1) * t109 - t236;
t99 = mrSges(2,2) * t223 - mrSges(2,3) * t211 - pkin(5) * t109 + t232 * t102 - t229 * t104;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t250 * t99 - t225 * t100 - qJ(1) * (t225 * t107 + t250 * t112), t99, t102, t106, t117, t130; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t225 * t99 + t250 * t100 + qJ(1) * (t250 * t107 - t225 * t112), t100, t104, t105, t116, t129; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t239, t239, t236, t251, -t235, -t238;];
m_new = t1;
