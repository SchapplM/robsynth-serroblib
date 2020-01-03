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
% Datum: 2020-01-03 11:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:45:11
% EndTime: 2020-01-03 11:45:15
% DurationCPUTime: 2.20s
% Computational Cost: add. (30713->221), mult. (41797->272), div. (0->0), fcn. (20254->8), ass. (0->91)
t193 = qJD(1) + qJD(3);
t204 = sin(qJ(4));
t207 = cos(qJ(4));
t168 = (-mrSges(6,1) * t207 + mrSges(6,2) * t204) * t193;
t192 = qJDD(1) + qJDD(3);
t229 = qJD(4) * t193;
t224 = t207 * t229;
t170 = t192 * t204 + t224;
t206 = sin(qJ(1));
t209 = cos(qJ(1));
t185 = -g(2) * t209 - g(3) * t206;
t176 = qJDD(1) * pkin(1) + t185;
t184 = -g(2) * t206 + t209 * g(3);
t210 = qJD(1) ^ 2;
t177 = -pkin(1) * t210 + t184;
t202 = sin(pkin(8));
t203 = cos(pkin(8));
t145 = t203 * t176 - t177 * t202;
t142 = qJDD(1) * pkin(2) + t145;
t146 = t202 * t176 + t203 * t177;
t143 = -pkin(2) * t210 + t146;
t205 = sin(qJ(3));
t208 = cos(qJ(3));
t138 = t205 * t142 + t208 * t143;
t191 = t193 ^ 2;
t135 = -pkin(3) * t191 + pkin(7) * t192 + t138;
t201 = -g(1) + qJDD(2);
t187 = t207 * t201;
t228 = qJD(5) * t193;
t236 = pkin(4) * t191;
t127 = qJDD(4) * pkin(4) + t187 + (-t170 + t224) * qJ(5) + (t207 * t236 - t135 - 0.2e1 * t228) * t204;
t233 = t193 * t207;
t181 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t233;
t227 = m(6) * t127 + qJDD(4) * mrSges(6,1) + qJD(4) * t181;
t234 = t193 * t204;
t122 = -t170 * mrSges(6,3) - t168 * t234 + t227;
t131 = -t204 * t135 + t187;
t132 = t135 * t207 + t204 * t201;
t154 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t204 + Ifges(5,2) * t207) * t193;
t155 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t204 + Ifges(6,4) * t207) * t193;
t156 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t204 + Ifges(5,4) * t207) * t193;
t171 = t192 * t207 - t204 * t229;
t178 = qJD(4) * pkin(4) - qJ(5) * t234;
t200 = t207 ^ 2;
t128 = qJ(5) * t171 - qJD(4) * t178 - t200 * t236 + 0.2e1 * t207 * t228 + t132;
t153 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t204 + Ifges(6,2) * t207) * t193;
t217 = -mrSges(6,1) * t127 + mrSges(6,2) * t128 - Ifges(6,5) * t170 - Ifges(6,6) * t171 - Ifges(6,3) * qJDD(4) - t153 * t234;
t238 = mrSges(5,1) * t131 - mrSges(5,2) * t132 + Ifges(5,5) * t170 + Ifges(5,6) * t171 + Ifges(5,3) * qJDD(4) + pkin(4) * t122 - (-t204 * t154 + (t155 + t156) * t207) * t193 - t217;
t235 = -mrSges(5,2) - mrSges(6,2);
t169 = (-mrSges(5,1) * t207 + mrSges(5,2) * t204) * t193;
t182 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t233;
t120 = m(5) * t131 + qJDD(4) * mrSges(5,1) + qJD(4) * t182 + (-t168 - t169) * t234 + (-mrSges(5,3) - mrSges(6,3)) * t170 + t227;
t226 = m(6) * t128 + t171 * mrSges(6,3) + t168 * t233;
t179 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t234;
t230 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t234 - t179;
t121 = m(5) * t132 + t171 * mrSges(5,3) + qJD(4) * t230 + qJDD(4) * t235 + t169 * t233 + t226;
t221 = -t120 * t204 + t121 * t207;
t111 = m(4) * t138 - mrSges(4,1) * t191 - mrSges(4,2) * t192 + t221;
t137 = t208 * t142 - t205 * t143;
t219 = -t192 * pkin(3) - t137;
t134 = -pkin(7) * t191 + t219;
t130 = t178 * t234 - t171 * pkin(4) + qJDD(5) + (-qJ(5) * t200 - pkin(7)) * t191 + t219;
t220 = -m(6) * t130 + t171 * mrSges(6,1) + t181 * t233;
t213 = -m(5) * t134 + t171 * mrSges(5,1) + t170 * t235 + t182 * t233 + t230 * t234 + t220;
t115 = m(4) * t137 + t192 * mrSges(4,1) - t191 * mrSges(4,2) + t213;
t104 = t111 * t205 + t115 * t208;
t101 = m(3) * t145 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t210 + t104;
t222 = t111 * t208 - t115 * t205;
t102 = m(3) * t146 - mrSges(3,1) * t210 - qJDD(1) * mrSges(3,2) + t222;
t96 = t101 * t203 + t102 * t202;
t113 = t120 * t207 + t121 * t204;
t225 = m(4) * t201 + t113;
t223 = -t101 * t202 + t102 * t203;
t218 = -mrSges(6,1) * t130 + mrSges(6,3) * t128 + Ifges(6,4) * t170 + Ifges(6,2) * t171 + Ifges(6,6) * qJDD(4) + qJD(4) * t155;
t151 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t204 + Ifges(6,6) * t207) * t193;
t216 = mrSges(6,2) * t130 - mrSges(6,3) * t127 + Ifges(6,1) * t170 + Ifges(6,4) * t171 + Ifges(6,5) * qJDD(4) + t151 * t233;
t152 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t204 + Ifges(5,6) * t207) * t193;
t106 = Ifges(5,4) * t170 + Ifges(5,2) * t171 + Ifges(5,6) * qJDD(4) + qJD(4) * t156 - mrSges(5,1) * t134 + mrSges(5,3) * t132 - pkin(4) * (t170 * mrSges(6,2) - t220) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t179 + t226) + (-pkin(4) * t179 - t151 - t152) * t234 + t218;
t108 = t152 * t233 + mrSges(5,2) * t134 - mrSges(5,3) * t131 + Ifges(5,1) * t170 + Ifges(5,4) * t171 + Ifges(5,5) * qJDD(4) - qJ(5) * t122 + (-t153 - t154) * qJD(4) + t216;
t215 = mrSges(4,1) * t137 - mrSges(4,2) * t138 + Ifges(4,3) * t192 + pkin(3) * t213 + pkin(7) * t221 + t106 * t207 + t108 * t204;
t214 = mrSges(3,1) * t145 - mrSges(3,2) * t146 + Ifges(3,3) * qJDD(1) + pkin(2) * t104 + t215;
t212 = mrSges(2,1) * t185 - mrSges(2,2) * t184 + Ifges(2,3) * qJDD(1) + pkin(1) * t96 + t214;
t97 = -mrSges(4,1) * t201 + mrSges(4,3) * t138 + t191 * Ifges(4,5) + Ifges(4,6) * t192 - pkin(3) * t113 - t238;
t94 = m(2) * t185 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t210 + t96;
t93 = m(2) * t184 - mrSges(2,1) * t210 - qJDD(1) * mrSges(2,2) + t223;
t92 = mrSges(4,2) * t201 - mrSges(4,3) * t137 + Ifges(4,5) * t192 - Ifges(4,6) * t191 - pkin(7) * t113 - t106 * t204 + t108 * t207;
t91 = mrSges(3,2) * t201 - mrSges(3,3) * t145 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t210 - pkin(6) * t104 - t205 * t97 + t208 * t92;
t90 = -mrSges(3,1) * t201 + mrSges(3,3) * t146 + t210 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t225 + pkin(6) * t222 + t205 * t92 + t208 * t97;
t89 = -mrSges(2,2) * g(1) - mrSges(2,3) * t185 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t210 - qJ(2) * t96 - t202 * t90 + t203 * t91;
t88 = Ifges(2,6) * qJDD(1) + t210 * Ifges(2,5) + mrSges(2,1) * g(1) + mrSges(2,3) * t184 + t202 * t91 + t203 * t90 - pkin(1) * (m(3) * t201 + t225) + qJ(2) * t223;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t212, t89, t91, t92, t108, -qJD(4) * t153 + t216; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t206 * t89 + t209 * t88 - pkin(5) * (t206 * t94 - t209 * t93), t88, t90, t97, t106, -t151 * t234 + t218; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) - t209 * t89 + t206 * t88 + pkin(5) * (t206 * t93 + t209 * t94), t212, t214, t215, t238, -t155 * t233 - t217;];
m_new = t1;
