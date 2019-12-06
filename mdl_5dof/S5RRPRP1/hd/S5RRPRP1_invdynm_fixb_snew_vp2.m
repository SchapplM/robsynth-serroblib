% Calculate vector of cutting torques with Newton-Euler for
% S5RRPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPRP1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:04
% EndTime: 2019-12-05 18:22:08
% DurationCPUTime: 2.15s
% Computational Cost: add. (32487->221), mult. (41797->272), div. (0->0), fcn. (20254->8), ass. (0->91)
t197 = qJD(1) + qJD(2);
t207 = sin(qJ(4));
t210 = cos(qJ(4));
t169 = (-mrSges(6,1) * t210 + mrSges(6,2) * t207) * t197;
t196 = qJDD(1) + qJDD(2);
t232 = qJD(4) * t197;
t227 = t210 * t232;
t171 = t207 * t196 + t227;
t209 = sin(qJ(1));
t212 = cos(qJ(1));
t186 = t212 * g(2) + t209 * g(3);
t177 = qJDD(1) * pkin(1) + t186;
t185 = t209 * g(2) - t212 * g(3);
t213 = qJD(1) ^ 2;
t178 = -t213 * pkin(1) + t185;
t208 = sin(qJ(2));
t211 = cos(qJ(2));
t146 = t211 * t177 - t208 * t178;
t143 = t196 * pkin(2) + t146;
t147 = t208 * t177 + t211 * t178;
t195 = t197 ^ 2;
t144 = -t195 * pkin(2) + t147;
t205 = sin(pkin(8));
t206 = cos(pkin(8));
t139 = t205 * t143 + t206 * t144;
t136 = -t195 * pkin(3) + t196 * pkin(7) + t139;
t204 = -g(1) + qJDD(3);
t188 = t210 * t204;
t231 = qJD(5) * t197;
t239 = pkin(4) * t195;
t128 = qJDD(4) * pkin(4) + t188 + (-t171 + t227) * qJ(5) + (t210 * t239 - t136 - 0.2e1 * t231) * t207;
t236 = t197 * t210;
t182 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t236;
t230 = m(6) * t128 + qJDD(4) * mrSges(6,1) + qJD(4) * t182;
t237 = t197 * t207;
t123 = -t171 * mrSges(6,3) - t169 * t237 + t230;
t132 = -t207 * t136 + t188;
t133 = t210 * t136 + t207 * t204;
t155 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t207 + Ifges(5,2) * t210) * t197;
t156 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t207 + Ifges(6,4) * t210) * t197;
t157 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t207 + Ifges(5,4) * t210) * t197;
t172 = t210 * t196 - t207 * t232;
t179 = qJD(4) * pkin(4) - qJ(5) * t237;
t203 = t210 ^ 2;
t129 = t172 * qJ(5) - qJD(4) * t179 - t203 * t239 + 0.2e1 * t210 * t231 + t133;
t154 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t207 + Ifges(6,2) * t210) * t197;
t220 = -mrSges(6,1) * t128 + mrSges(6,2) * t129 - Ifges(6,5) * t171 - Ifges(6,6) * t172 - Ifges(6,3) * qJDD(4) - t154 * t237;
t241 = mrSges(5,1) * t132 - mrSges(5,2) * t133 + Ifges(5,5) * t171 + Ifges(5,6) * t172 + Ifges(5,3) * qJDD(4) + pkin(4) * t123 - (-t207 * t155 + (t156 + t157) * t210) * t197 - t220;
t238 = -mrSges(5,2) - mrSges(6,2);
t170 = (-mrSges(5,1) * t210 + mrSges(5,2) * t207) * t197;
t183 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t236;
t121 = m(5) * t132 + qJDD(4) * mrSges(5,1) + qJD(4) * t183 + (-t169 - t170) * t237 + (-mrSges(5,3) - mrSges(6,3)) * t171 + t230;
t229 = m(6) * t129 + t172 * mrSges(6,3) + t169 * t236;
t180 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t237;
t233 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t237 - t180;
t122 = m(5) * t133 + t172 * mrSges(5,3) + t233 * qJD(4) + t238 * qJDD(4) + t170 * t236 + t229;
t224 = -t207 * t121 + t210 * t122;
t112 = m(4) * t139 - t195 * mrSges(4,1) - t196 * mrSges(4,2) + t224;
t138 = t206 * t143 - t205 * t144;
t222 = -t196 * pkin(3) - t138;
t135 = -t195 * pkin(7) + t222;
t131 = t179 * t237 - t172 * pkin(4) + qJDD(5) + (-qJ(5) * t203 - pkin(7)) * t195 + t222;
t223 = -m(6) * t131 + t172 * mrSges(6,1) + t182 * t236;
t216 = -m(5) * t135 + t172 * mrSges(5,1) + t238 * t171 + t183 * t236 + t233 * t237 + t223;
t116 = m(4) * t138 + t196 * mrSges(4,1) - t195 * mrSges(4,2) + t216;
t105 = t205 * t112 + t206 * t116;
t102 = m(3) * t146 + t196 * mrSges(3,1) - t195 * mrSges(3,2) + t105;
t225 = t206 * t112 - t205 * t116;
t103 = m(3) * t147 - t195 * mrSges(3,1) - t196 * mrSges(3,2) + t225;
t97 = t211 * t102 + t208 * t103;
t114 = t210 * t121 + t207 * t122;
t228 = m(4) * t204 + t114;
t226 = -t208 * t102 + t211 * t103;
t221 = -mrSges(6,1) * t131 + mrSges(6,3) * t129 + Ifges(6,4) * t171 + Ifges(6,2) * t172 + Ifges(6,6) * qJDD(4) + qJD(4) * t156;
t152 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t207 + Ifges(6,6) * t210) * t197;
t219 = mrSges(6,2) * t131 - mrSges(6,3) * t128 + Ifges(6,1) * t171 + Ifges(6,4) * t172 + Ifges(6,5) * qJDD(4) + t152 * t236;
t153 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t207 + Ifges(5,6) * t210) * t197;
t107 = Ifges(5,4) * t171 + Ifges(5,2) * t172 + Ifges(5,6) * qJDD(4) + qJD(4) * t157 - mrSges(5,1) * t135 + mrSges(5,3) * t133 - pkin(4) * (t171 * mrSges(6,2) - t223) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t180 + t229) + (-pkin(4) * t180 - t152 - t153) * t237 + t221;
t109 = t153 * t236 + mrSges(5,2) * t135 - mrSges(5,3) * t132 + Ifges(5,1) * t171 + Ifges(5,4) * t172 + Ifges(5,5) * qJDD(4) - qJ(5) * t123 + (-t154 - t155) * qJD(4) + t219;
t218 = mrSges(4,1) * t138 - mrSges(4,2) * t139 + Ifges(4,3) * t196 + pkin(3) * t216 + pkin(7) * t224 + t210 * t107 + t207 * t109;
t217 = mrSges(3,1) * t146 - mrSges(3,2) * t147 + Ifges(3,3) * t196 + pkin(2) * t105 + t218;
t215 = mrSges(2,1) * t186 - mrSges(2,2) * t185 + Ifges(2,3) * qJDD(1) + pkin(1) * t97 + t217;
t98 = -mrSges(4,1) * t204 + mrSges(4,3) * t139 + t195 * Ifges(4,5) + Ifges(4,6) * t196 - pkin(3) * t114 - t241;
t95 = m(2) * t186 + qJDD(1) * mrSges(2,1) - t213 * mrSges(2,2) + t97;
t94 = m(2) * t185 - t213 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t226;
t93 = mrSges(4,2) * t204 - mrSges(4,3) * t138 + Ifges(4,5) * t196 - t195 * Ifges(4,6) - pkin(7) * t114 - t207 * t107 + t210 * t109;
t92 = -mrSges(3,2) * g(1) - mrSges(3,3) * t146 + Ifges(3,5) * t196 - t195 * Ifges(3,6) - qJ(3) * t105 - t205 * t98 + t206 * t93;
t91 = mrSges(3,1) * g(1) + mrSges(3,3) * t147 + t195 * Ifges(3,5) + Ifges(3,6) * t196 - pkin(2) * t228 + qJ(3) * t225 + t205 * t93 + t206 * t98;
t90 = -mrSges(2,2) * g(1) - mrSges(2,3) * t186 + Ifges(2,5) * qJDD(1) - t213 * Ifges(2,6) - pkin(6) * t97 - t208 * t91 + t211 * t92;
t89 = Ifges(2,6) * qJDD(1) + t213 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t185 + t208 * t92 + t211 * t91 - pkin(1) * (-m(3) * g(1) + t228) + pkin(6) * t226;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t215, t90, t92, t93, t109, -qJD(4) * t154 + t219; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t209 * t90 - t212 * t89 - pkin(5) * (-t209 * t95 + t212 * t94), t89, t91, t98, t107, -t152 * t237 + t221; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t212 * t90 - t209 * t89 + pkin(5) * (-t209 * t94 - t212 * t95), t215, t217, t218, t241, -t156 * t236 - t220;];
m_new = t1;
