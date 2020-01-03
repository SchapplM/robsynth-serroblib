% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRP7
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
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:44:38
% EndTime: 2019-12-31 18:44:43
% DurationCPUTime: 2.78s
% Computational Cost: add. (34745->262), mult. (65045->318), div. (0->0), fcn. (37063->8), ass. (0->99)
t221 = sin(qJ(1));
t223 = cos(qJ(1));
t206 = t221 * g(1) - g(2) * t223;
t197 = qJDD(1) * pkin(1) + t206;
t207 = -g(1) * t223 - g(2) * t221;
t225 = qJD(1) ^ 2;
t199 = -pkin(1) * t225 + t207;
t217 = sin(pkin(8));
t218 = cos(pkin(8));
t170 = t218 * t197 - t217 * t199;
t151 = -qJDD(1) * pkin(2) - t225 * pkin(6) - t170;
t220 = sin(qJ(3));
t222 = cos(qJ(3));
t241 = qJD(1) * qJD(3);
t238 = t222 * t241;
t201 = qJDD(1) * t220 + t238;
t239 = t220 * t241;
t202 = qJDD(1) * t222 - t239;
t141 = (-t201 - t238) * pkin(7) + (-t202 + t239) * pkin(3) + t151;
t171 = t217 * t197 + t218 * t199;
t152 = -pkin(2) * t225 + qJDD(1) * pkin(6) + t171;
t216 = -g(3) + qJDD(2);
t147 = t222 * t152 + t220 * t216;
t200 = (-pkin(3) * t222 - pkin(7) * t220) * qJD(1);
t224 = qJD(3) ^ 2;
t242 = t222 * qJD(1);
t144 = -pkin(3) * t224 + qJDD(3) * pkin(7) + t200 * t242 + t147;
t219 = sin(qJ(4));
t248 = cos(qJ(4));
t138 = t141 * t248 - t144 * t219;
t139 = t141 * t219 + t144 * t248;
t243 = qJD(1) * t220;
t195 = -qJD(3) * t248 + t219 * t243;
t196 = qJD(3) * t219 + t243 * t248;
t209 = qJD(4) - t242;
t153 = Ifges(6,5) * t196 + Ifges(6,6) * t209 + Ifges(6,3) * t195;
t156 = Ifges(5,4) * t196 - Ifges(5,2) * t195 + Ifges(5,6) * t209;
t158 = Ifges(5,1) * t196 - Ifges(5,4) * t195 + Ifges(5,5) * t209;
t167 = qJD(4) * t196 - qJDD(3) * t248 + t201 * t219;
t168 = -qJD(4) * t195 + qJDD(3) * t219 + t201 * t248;
t175 = mrSges(6,1) * t195 - mrSges(6,3) * t196;
t194 = qJDD(4) - t202;
t174 = pkin(4) * t195 - qJ(5) * t196;
t208 = t209 ^ 2;
t134 = -pkin(4) * t208 + qJ(5) * t194 + 0.2e1 * qJD(5) * t209 - t174 * t195 + t139;
t136 = -pkin(4) * t194 - qJ(5) * t208 + t174 * t196 + qJDD(5) - t138;
t157 = Ifges(6,1) * t196 + Ifges(6,4) * t209 + Ifges(6,5) * t195;
t232 = mrSges(6,1) * t136 - mrSges(6,3) * t134 - Ifges(6,4) * t168 - Ifges(6,2) * t194 - Ifges(6,6) * t167 - t157 * t195;
t180 = -mrSges(6,2) * t195 + mrSges(6,3) * t209;
t235 = -m(6) * t136 + t194 * mrSges(6,1) + t209 * t180;
t179 = -mrSges(6,1) * t209 + mrSges(6,2) * t196;
t240 = m(6) * t134 + t194 * mrSges(6,3) + t209 * t179;
t250 = -(-t156 + t153) * t196 + mrSges(5,1) * t138 - mrSges(5,2) * t139 + Ifges(5,5) * t168 - Ifges(5,6) * t167 + Ifges(5,3) * t194 + pkin(4) * (-t168 * mrSges(6,2) - t196 * t175 + t235) + qJ(5) * (-t167 * mrSges(6,2) - t195 * t175 + t240) + t195 * t158 - t232;
t146 = -t220 * t152 + t222 * t216;
t143 = -qJDD(3) * pkin(3) - t224 * pkin(7) + t200 * t243 - t146;
t137 = -0.2e1 * qJD(5) * t196 + (t195 * t209 - t168) * qJ(5) + (t196 * t209 + t167) * pkin(4) + t143;
t131 = m(6) * t137 + t167 * mrSges(6,1) - t168 * mrSges(6,3) - t196 * t179 + t180 * t195;
t234 = -mrSges(6,1) * t137 + mrSges(6,2) * t134;
t155 = Ifges(6,4) * t196 + Ifges(6,2) * t209 + Ifges(6,6) * t195;
t246 = -Ifges(5,5) * t196 + Ifges(5,6) * t195 - Ifges(5,3) * t209 - t155;
t115 = -mrSges(5,1) * t143 + mrSges(5,3) * t139 - pkin(4) * t131 + (t157 + t158) * t209 + t246 * t196 + (Ifges(5,6) - Ifges(6,6)) * t194 + (Ifges(5,4) - Ifges(6,5)) * t168 + (-Ifges(5,2) - Ifges(6,3)) * t167 + t234;
t231 = mrSges(6,2) * t136 - mrSges(6,3) * t137 + Ifges(6,1) * t168 + Ifges(6,4) * t194 + Ifges(6,5) * t167 + t209 * t153;
t116 = mrSges(5,2) * t143 - mrSges(5,3) * t138 + Ifges(5,1) * t168 - Ifges(5,4) * t167 + Ifges(5,5) * t194 - qJ(5) * t131 - t209 * t156 + t195 * t246 + t231;
t178 = mrSges(5,1) * t209 - mrSges(5,3) * t196;
t244 = -mrSges(5,1) * t195 - mrSges(5,2) * t196 - t175;
t247 = -mrSges(5,3) - mrSges(6,2);
t126 = m(5) * t139 - t194 * mrSges(5,2) + t167 * t247 - t209 * t178 + t195 * t244 + t240;
t177 = -mrSges(5,2) * t209 - mrSges(5,3) * t195;
t127 = m(5) * t138 + t194 * mrSges(5,1) + t168 * t247 + t209 * t177 + t196 * t244 + t235;
t122 = t126 * t248 - t127 * t219;
t128 = -m(5) * t143 - t167 * mrSges(5,1) - t168 * mrSges(5,2) - t195 * t177 - t178 * t196 - t131;
t186 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t220 + Ifges(4,2) * t222) * qJD(1);
t187 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t220 + Ifges(4,4) * t222) * qJD(1);
t249 = mrSges(4,1) * t146 - mrSges(4,2) * t147 + Ifges(4,5) * t201 + Ifges(4,6) * t202 + Ifges(4,3) * qJDD(3) + pkin(3) * t128 + pkin(7) * t122 + (t186 * t220 - t187 * t222) * qJD(1) + t248 * t115 + t219 * t116;
t198 = (-mrSges(4,1) * t222 + mrSges(4,2) * t220) * qJD(1);
t204 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t243;
t120 = m(4) * t147 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t202 - qJD(3) * t204 + t198 * t242 + t122;
t205 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t242;
t124 = m(4) * t146 + qJDD(3) * mrSges(4,1) - mrSges(4,3) * t201 + qJD(3) * t205 - t198 * t243 + t128;
t236 = t120 * t222 - t124 * t220;
t110 = m(3) * t171 - mrSges(3,1) * t225 - qJDD(1) * mrSges(3,2) + t236;
t121 = t126 * t219 + t127 * t248;
t228 = -m(4) * t151 + t202 * mrSges(4,1) - mrSges(4,2) * t201 - t204 * t243 + t205 * t242 - t121;
t114 = m(3) * t170 + qJDD(1) * mrSges(3,1) - mrSges(3,2) * t225 + t228;
t105 = t110 * t217 + t114 * t218;
t112 = t120 * t220 + t124 * t222;
t237 = t110 * t218 - t114 * t217;
t185 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t220 + Ifges(4,6) * t222) * qJD(1);
t101 = mrSges(4,2) * t151 - mrSges(4,3) * t146 + Ifges(4,1) * t201 + Ifges(4,4) * t202 + Ifges(4,5) * qJDD(3) - pkin(7) * t121 - qJD(3) * t186 - t115 * t219 + t116 * t248 + t185 * t242;
t107 = -mrSges(4,1) * t151 + mrSges(4,3) * t147 + Ifges(4,4) * t201 + Ifges(4,2) * t202 + Ifges(4,6) * qJDD(3) - pkin(3) * t121 + qJD(3) * t187 - t185 * t243 - t250;
t230 = mrSges(3,1) * t170 - mrSges(3,2) * t171 + Ifges(3,3) * qJDD(1) + pkin(2) * t228 + pkin(6) * t236 + t101 * t220 + t107 * t222;
t229 = mrSges(2,1) * t206 - mrSges(2,2) * t207 + Ifges(2,3) * qJDD(1) + pkin(1) * t105 + t230;
t103 = m(2) * t207 - mrSges(2,1) * t225 - qJDD(1) * mrSges(2,2) + t237;
t102 = m(2) * t206 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t225 + t105;
t99 = -mrSges(3,1) * t216 + mrSges(3,3) * t171 + t225 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t112 - t249;
t98 = mrSges(3,2) * t216 - mrSges(3,3) * t170 + Ifges(3,5) * qJDD(1) - Ifges(3,6) * t225 - pkin(6) * t112 + t101 * t222 - t107 * t220;
t97 = -mrSges(2,2) * g(3) - mrSges(2,3) * t206 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t225 - qJ(2) * t105 - t217 * t99 + t218 * t98;
t96 = Ifges(2,6) * qJDD(1) + t225 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t207 + t217 * t98 + t218 * t99 - pkin(1) * (m(3) * t216 + t112) + qJ(2) * t237;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t223 * t97 - t221 * t96 - pkin(5) * (t102 * t223 + t103 * t221), t97, t98, t101, t116, -t155 * t195 + t231; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t221 * t97 + t223 * t96 + pkin(5) * (-t102 * t221 + t103 * t223), t96, t99, t107, t115, -t196 * t153 - t232; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t229, t229, t230, t249, t250, Ifges(6,5) * t168 + Ifges(6,6) * t194 + Ifges(6,3) * t167 + t196 * t155 - t209 * t157 - t234;];
m_new = t1;
