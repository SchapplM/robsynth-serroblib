% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:01:29
% EndTime: 2019-12-05 18:01:32
% DurationCPUTime: 2.07s
% Computational Cost: add. (30713->221), mult. (41797->272), div. (0->0), fcn. (20254->8), ass. (0->91)
t197 = qJD(1) + qJD(3);
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t170 = (-mrSges(6,1) * t211 + mrSges(6,2) * t208) * t197;
t196 = qJDD(1) + qJDD(3);
t233 = qJD(4) * t197;
t228 = t211 * t233;
t172 = t196 * t208 + t228;
t210 = sin(qJ(1));
t213 = cos(qJ(1));
t187 = t213 * g(2) + t210 * g(3);
t178 = qJDD(1) * pkin(1) + t187;
t186 = t210 * g(2) - g(3) * t213;
t214 = qJD(1) ^ 2;
t179 = -pkin(1) * t214 + t186;
t206 = sin(pkin(8));
t207 = cos(pkin(8));
t147 = t207 * t178 - t179 * t206;
t144 = qJDD(1) * pkin(2) + t147;
t148 = t206 * t178 + t207 * t179;
t145 = -pkin(2) * t214 + t148;
t209 = sin(qJ(3));
t212 = cos(qJ(3));
t140 = t209 * t144 + t212 * t145;
t195 = t197 ^ 2;
t137 = -pkin(3) * t195 + pkin(7) * t196 + t140;
t205 = -g(1) + qJDD(2);
t189 = t211 * t205;
t232 = qJD(5) * t197;
t240 = pkin(4) * t195;
t129 = qJDD(4) * pkin(4) + t189 + (-t172 + t228) * qJ(5) + (t211 * t240 - t137 - 0.2e1 * t232) * t208;
t237 = t197 * t211;
t183 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t237;
t231 = m(6) * t129 + qJDD(4) * mrSges(6,1) + qJD(4) * t183;
t238 = t197 * t208;
t124 = -t172 * mrSges(6,3) - t170 * t238 + t231;
t133 = -t208 * t137 + t189;
t134 = t211 * t137 + t208 * t205;
t156 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t208 + Ifges(5,2) * t211) * t197;
t157 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t208 + Ifges(6,4) * t211) * t197;
t158 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t208 + Ifges(5,4) * t211) * t197;
t173 = t196 * t211 - t208 * t233;
t180 = qJD(4) * pkin(4) - qJ(5) * t238;
t204 = t211 ^ 2;
t130 = qJ(5) * t173 - qJD(4) * t180 - t204 * t240 + 0.2e1 * t211 * t232 + t134;
t155 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t208 + Ifges(6,2) * t211) * t197;
t221 = -mrSges(6,1) * t129 + mrSges(6,2) * t130 - Ifges(6,5) * t172 - Ifges(6,6) * t173 - Ifges(6,3) * qJDD(4) - t155 * t238;
t242 = mrSges(5,1) * t133 - mrSges(5,2) * t134 + Ifges(5,5) * t172 + Ifges(5,6) * t173 + Ifges(5,3) * qJDD(4) + pkin(4) * t124 - (-t208 * t156 + (t157 + t158) * t211) * t197 - t221;
t239 = -mrSges(5,2) - mrSges(6,2);
t171 = (-mrSges(5,1) * t211 + mrSges(5,2) * t208) * t197;
t184 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t237;
t122 = m(5) * t133 + qJDD(4) * mrSges(5,1) + qJD(4) * t184 + (-t170 - t171) * t238 + (-mrSges(5,3) - mrSges(6,3)) * t172 + t231;
t230 = m(6) * t130 + t173 * mrSges(6,3) + t170 * t237;
t181 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t238;
t234 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t238 - t181;
t123 = m(5) * t134 + mrSges(5,3) * t173 + t234 * qJD(4) + t239 * qJDD(4) + t171 * t237 + t230;
t225 = -t122 * t208 + t211 * t123;
t113 = m(4) * t140 - mrSges(4,1) * t195 - mrSges(4,2) * t196 + t225;
t139 = t212 * t144 - t209 * t145;
t223 = -t196 * pkin(3) - t139;
t136 = -pkin(7) * t195 + t223;
t132 = t180 * t238 - t173 * pkin(4) + qJDD(5) + (-qJ(5) * t204 - pkin(7)) * t195 + t223;
t224 = -m(6) * t132 + t173 * mrSges(6,1) + t183 * t237;
t217 = -m(5) * t136 + t173 * mrSges(5,1) + t239 * t172 + t184 * t237 + t234 * t238 + t224;
t117 = m(4) * t139 + mrSges(4,1) * t196 - mrSges(4,2) * t195 + t217;
t106 = t209 * t113 + t212 * t117;
t103 = m(3) * t147 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t214 + t106;
t226 = t212 * t113 - t117 * t209;
t104 = m(3) * t148 - mrSges(3,1) * t214 - qJDD(1) * mrSges(3,2) + t226;
t98 = t207 * t103 + t206 * t104;
t115 = t211 * t122 + t208 * t123;
t229 = m(4) * t205 + t115;
t227 = -t206 * t103 + t207 * t104;
t222 = -mrSges(6,1) * t132 + mrSges(6,3) * t130 + Ifges(6,4) * t172 + Ifges(6,2) * t173 + Ifges(6,6) * qJDD(4) + qJD(4) * t157;
t153 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t208 + Ifges(6,6) * t211) * t197;
t220 = mrSges(6,2) * t132 - mrSges(6,3) * t129 + Ifges(6,1) * t172 + Ifges(6,4) * t173 + Ifges(6,5) * qJDD(4) + t153 * t237;
t154 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t208 + Ifges(5,6) * t211) * t197;
t108 = Ifges(5,4) * t172 + Ifges(5,2) * t173 + Ifges(5,6) * qJDD(4) + qJD(4) * t158 - mrSges(5,1) * t136 + mrSges(5,3) * t134 - pkin(4) * (t172 * mrSges(6,2) - t224) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t181 + t230) + (-pkin(4) * t181 - t153 - t154) * t238 + t222;
t110 = t154 * t237 + mrSges(5,2) * t136 - mrSges(5,3) * t133 + Ifges(5,1) * t172 + Ifges(5,4) * t173 + Ifges(5,5) * qJDD(4) - qJ(5) * t124 + (-t155 - t156) * qJD(4) + t220;
t219 = mrSges(4,1) * t139 - mrSges(4,2) * t140 + Ifges(4,3) * t196 + pkin(3) * t217 + pkin(7) * t225 + t211 * t108 + t208 * t110;
t218 = mrSges(3,1) * t147 - mrSges(3,2) * t148 + Ifges(3,3) * qJDD(1) + pkin(2) * t106 + t219;
t216 = mrSges(2,1) * t187 - mrSges(2,2) * t186 + Ifges(2,3) * qJDD(1) + pkin(1) * t98 + t218;
t99 = -mrSges(4,1) * t205 + mrSges(4,3) * t140 + t195 * Ifges(4,5) + Ifges(4,6) * t196 - pkin(3) * t115 - t242;
t96 = m(2) * t187 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t214 + t98;
t95 = m(2) * t186 - mrSges(2,1) * t214 - qJDD(1) * mrSges(2,2) + t227;
t94 = mrSges(4,2) * t205 - mrSges(4,3) * t139 + Ifges(4,5) * t196 - Ifges(4,6) * t195 - pkin(7) * t115 - t108 * t208 + t110 * t211;
t93 = mrSges(3,2) * t205 - mrSges(3,3) * t147 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t214 - pkin(6) * t106 - t209 * t99 + t212 * t94;
t92 = -mrSges(3,1) * t205 + mrSges(3,3) * t148 + t214 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t229 + pkin(6) * t226 + t209 * t94 + t212 * t99;
t91 = -mrSges(2,2) * g(1) - mrSges(2,3) * t187 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t214 - qJ(2) * t98 - t206 * t92 + t207 * t93;
t90 = Ifges(2,6) * qJDD(1) + t214 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t186 + t206 * t93 + t207 * t92 - pkin(1) * (m(3) * t205 + t229) + qJ(2) * t227;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t216, t91, t93, t94, t110, -qJD(4) * t155 + t220; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) - t210 * t91 - t213 * t90 - pkin(5) * (-t210 * t96 + t213 * t95), t90, t92, t99, t108, -t153 * t238 + t222; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t213 * t91 - t210 * t90 + pkin(5) * (-t210 * t95 - t213 * t96), t216, t218, t219, t242, -t157 * t237 - t221;];
m_new = t1;
