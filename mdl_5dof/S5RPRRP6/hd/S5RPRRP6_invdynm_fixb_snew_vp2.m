% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP6
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
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:42:14
% EndTime: 2019-12-31 18:42:19
% DurationCPUTime: 2.82s
% Computational Cost: add. (35776->261), mult. (67561->318), div. (0->0), fcn. (38740->8), ass. (0->100)
t221 = sin(qJ(4));
t224 = cos(qJ(4));
t222 = sin(qJ(3));
t247 = qJD(1) * t222;
t200 = t224 * qJD(3) - t221 * t247;
t225 = cos(qJ(3));
t245 = qJD(1) * qJD(3);
t241 = t225 * t245;
t206 = t222 * qJDD(1) + t241;
t173 = t200 * qJD(4) + t221 * qJDD(3) + t224 * t206;
t201 = t221 * qJD(3) + t224 * t247;
t178 = -t200 * mrSges(6,1) + t201 * mrSges(6,2);
t223 = sin(qJ(1));
t226 = cos(qJ(1));
t211 = t223 * g(1) - t226 * g(2);
t202 = qJDD(1) * pkin(1) + t211;
t212 = -t226 * g(1) - t223 * g(2);
t228 = qJD(1) ^ 2;
t204 = -t228 * pkin(1) + t212;
t219 = sin(pkin(8));
t220 = cos(pkin(8));
t175 = t220 * t202 - t219 * t204;
t153 = -qJDD(1) * pkin(2) - t228 * pkin(6) - t175;
t242 = t222 * t245;
t207 = t225 * qJDD(1) - t242;
t141 = (-t206 - t241) * pkin(7) + (-t207 + t242) * pkin(3) + t153;
t176 = t219 * t202 + t220 * t204;
t154 = -t228 * pkin(2) + qJDD(1) * pkin(6) + t176;
t218 = -g(3) + qJDD(2);
t148 = t225 * t154 + t222 * t218;
t205 = (-pkin(3) * t225 - pkin(7) * t222) * qJD(1);
t227 = qJD(3) ^ 2;
t246 = t225 * qJD(1);
t144 = -t227 * pkin(3) + qJDD(3) * pkin(7) + t205 * t246 + t148;
t137 = t224 * t141 - t221 * t144;
t199 = qJDD(4) - t207;
t213 = qJD(4) - t246;
t131 = -0.2e1 * qJD(5) * t201 + (t200 * t213 - t173) * qJ(5) + (t200 * t201 + t199) * pkin(4) + t137;
t180 = -t213 * mrSges(6,2) + t200 * mrSges(6,3);
t244 = m(6) * t131 + t199 * mrSges(6,1) + t213 * t180;
t128 = -t173 * mrSges(6,3) - t201 * t178 + t244;
t138 = t221 * t141 + t224 * t144;
t158 = Ifges(5,4) * t201 + Ifges(5,2) * t200 + Ifges(5,6) * t213;
t159 = Ifges(6,1) * t201 + Ifges(6,4) * t200 + Ifges(6,5) * t213;
t160 = Ifges(5,1) * t201 + Ifges(5,4) * t200 + Ifges(5,5) * t213;
t172 = -t201 * qJD(4) + t224 * qJDD(3) - t221 * t206;
t182 = t213 * pkin(4) - t201 * qJ(5);
t198 = t200 ^ 2;
t134 = -t198 * pkin(4) + t172 * qJ(5) + 0.2e1 * qJD(5) * t200 - t213 * t182 + t138;
t157 = Ifges(6,4) * t201 + Ifges(6,2) * t200 + Ifges(6,6) * t213;
t235 = -mrSges(6,1) * t131 + mrSges(6,2) * t134 - Ifges(6,5) * t173 - Ifges(6,6) * t172 - Ifges(6,3) * t199 - t201 * t157;
t252 = mrSges(5,1) * t137 - mrSges(5,2) * t138 + Ifges(5,5) * t173 + Ifges(5,6) * t172 + Ifges(5,3) * t199 + pkin(4) * t128 + t201 * t158 - (t160 + t159) * t200 - t235;
t147 = -t222 * t154 + t225 * t218;
t143 = -qJDD(3) * pkin(3) - t227 * pkin(7) + t205 * t247 - t147;
t155 = Ifges(6,5) * t201 + Ifges(6,6) * t200 + Ifges(6,3) * t213;
t156 = Ifges(5,5) * t201 + Ifges(5,6) * t200 + Ifges(5,3) * t213;
t183 = t213 * mrSges(6,1) - t201 * mrSges(6,3);
t136 = -t172 * pkin(4) - t198 * qJ(5) + t201 * t182 + qJDD(5) + t143;
t236 = -mrSges(6,1) * t136 + mrSges(6,3) * t134 + Ifges(6,4) * t173 + Ifges(6,2) * t172 + Ifges(6,6) * t199 + t213 * t159;
t238 = -m(6) * t136 + t172 * mrSges(6,1) + t200 * t180;
t243 = m(6) * t134 + t172 * mrSges(6,3) + t200 * t178;
t112 = Ifges(5,4) * t173 + Ifges(5,2) * t172 + Ifges(5,6) * t199 + t213 * t160 - mrSges(5,1) * t143 + mrSges(5,3) * t138 - pkin(4) * (t173 * mrSges(6,2) - t238) + qJ(5) * (-t199 * mrSges(6,2) - t213 * t183 + t243) + (-pkin(4) * t183 - t155 - t156) * t201 + t236;
t234 = mrSges(6,2) * t136 - mrSges(6,3) * t131 + Ifges(6,1) * t173 + Ifges(6,4) * t172 + Ifges(6,5) * t199 + t200 * t155;
t119 = mrSges(5,2) * t143 - mrSges(5,3) * t137 + Ifges(5,1) * t173 + Ifges(5,4) * t172 + Ifges(5,5) * t199 - qJ(5) * t128 + t200 * t156 + (-t157 - t158) * t213 + t234;
t179 = -t200 * mrSges(5,1) + t201 * mrSges(5,2);
t181 = -t213 * mrSges(5,2) + t200 * mrSges(5,3);
t122 = m(5) * t137 + t199 * mrSges(5,1) + t213 * t181 + (-t178 - t179) * t201 + (-mrSges(5,3) - mrSges(6,3)) * t173 + t244;
t248 = -t213 * mrSges(5,1) + t201 * mrSges(5,3) - t183;
t250 = -mrSges(5,2) - mrSges(6,2);
t124 = m(5) * t138 + t172 * mrSges(5,3) + t200 * t179 + t250 * t199 + t248 * t213 + t243;
t121 = -t221 * t122 + t224 * t124;
t127 = -m(5) * t143 + t172 * mrSges(5,1) + t250 * t173 + t200 * t181 + t248 * t201 + t238;
t190 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t222 + Ifges(4,2) * t225) * qJD(1);
t191 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t222 + Ifges(4,4) * t225) * qJD(1);
t251 = mrSges(4,1) * t147 - mrSges(4,2) * t148 + Ifges(4,5) * t206 + Ifges(4,6) * t207 + Ifges(4,3) * qJDD(3) + pkin(3) * t127 + pkin(7) * t121 + t224 * t112 + t221 * t119 + (t222 * t190 - t225 * t191) * qJD(1);
t203 = (-mrSges(4,1) * t225 + mrSges(4,2) * t222) * qJD(1);
t209 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t247;
t118 = m(4) * t148 - qJDD(3) * mrSges(4,2) + t207 * mrSges(4,3) - qJD(3) * t209 + t203 * t246 + t121;
t210 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t246;
t126 = m(4) * t147 + qJDD(3) * mrSges(4,1) - t206 * mrSges(4,3) + qJD(3) * t210 - t203 * t247 + t127;
t239 = t225 * t118 - t222 * t126;
t109 = m(3) * t176 - t228 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t239;
t120 = t224 * t122 + t221 * t124;
t231 = -m(4) * t153 + t207 * mrSges(4,1) - t206 * mrSges(4,2) - t209 * t247 + t210 * t246 - t120;
t114 = m(3) * t175 + qJDD(1) * mrSges(3,1) - t228 * mrSges(3,2) + t231;
t104 = t219 * t109 + t220 * t114;
t111 = t222 * t118 + t225 * t126;
t240 = t220 * t109 - t219 * t114;
t189 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t222 + Ifges(4,6) * t225) * qJD(1);
t100 = mrSges(4,2) * t153 - mrSges(4,3) * t147 + Ifges(4,1) * t206 + Ifges(4,4) * t207 + Ifges(4,5) * qJDD(3) - pkin(7) * t120 - qJD(3) * t190 - t221 * t112 + t224 * t119 + t189 * t246;
t106 = -mrSges(4,1) * t153 + mrSges(4,3) * t148 + Ifges(4,4) * t206 + Ifges(4,2) * t207 + Ifges(4,6) * qJDD(3) - pkin(3) * t120 + qJD(3) * t191 - t189 * t247 - t252;
t233 = mrSges(3,1) * t175 - mrSges(3,2) * t176 + Ifges(3,3) * qJDD(1) + pkin(2) * t231 + pkin(6) * t239 + t222 * t100 + t225 * t106;
t232 = mrSges(2,1) * t211 - mrSges(2,2) * t212 + Ifges(2,3) * qJDD(1) + pkin(1) * t104 + t233;
t102 = m(2) * t212 - t228 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t240;
t101 = m(2) * t211 + qJDD(1) * mrSges(2,1) - t228 * mrSges(2,2) + t104;
t98 = -mrSges(3,1) * t218 + mrSges(3,3) * t176 + t228 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t111 - t251;
t97 = mrSges(3,2) * t218 - mrSges(3,3) * t175 + Ifges(3,5) * qJDD(1) - t228 * Ifges(3,6) - pkin(6) * t111 + t225 * t100 - t222 * t106;
t96 = -mrSges(2,2) * g(3) - mrSges(2,3) * t211 + Ifges(2,5) * qJDD(1) - t228 * Ifges(2,6) - qJ(2) * t104 - t219 * t98 + t220 * t97;
t95 = Ifges(2,6) * qJDD(1) + t228 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t212 + t219 * t97 + t220 * t98 - pkin(1) * (m(3) * t218 + t111) + qJ(2) * t240;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t226 * t96 - t223 * t95 - pkin(5) * (t226 * t101 + t223 * t102), t96, t97, t100, t119, -t213 * t157 + t234; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t223 * t96 + t226 * t95 + pkin(5) * (-t223 * t101 + t226 * t102), t95, t98, t106, t112, -t201 * t155 + t236; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232, t232, t233, t251, t252, -t200 * t159 - t235;];
m_new = t1;
