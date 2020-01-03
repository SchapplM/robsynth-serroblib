% Calculate vector of cutting torques with Newton-Euler for
% S5RPRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPRRR7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:19
% EndTime: 2019-12-31 19:03:30
% DurationCPUTime: 6.18s
% Computational Cost: add. (85321->266), mult. (163822->336), div. (0->0), fcn. (102970->10), ass. (0->110)
t225 = sin(qJ(1));
t229 = cos(qJ(1));
t210 = t225 * g(1) - t229 * g(2);
t201 = qJDD(1) * pkin(1) + t210;
t211 = -t229 * g(1) - t225 * g(2);
t231 = qJD(1) ^ 2;
t203 = -t231 * pkin(1) + t211;
t220 = sin(pkin(9));
t221 = cos(pkin(9));
t180 = t221 * t201 - t220 * t203;
t168 = -qJDD(1) * pkin(2) - t231 * pkin(6) - t180;
t224 = sin(qJ(3));
t228 = cos(qJ(3));
t244 = qJD(1) * qJD(3);
t243 = t228 * t244;
t205 = t224 * qJDD(1) + t243;
t215 = t224 * t244;
t206 = t228 * qJDD(1) - t215;
t156 = (-t205 - t243) * pkin(7) + (-t206 + t215) * pkin(3) + t168;
t181 = t220 * t201 + t221 * t203;
t169 = -t231 * pkin(2) + qJDD(1) * pkin(6) + t181;
t219 = -g(3) + qJDD(2);
t163 = t228 * t169 + t224 * t219;
t204 = (-pkin(3) * t228 - pkin(7) * t224) * qJD(1);
t230 = qJD(3) ^ 2;
t245 = t228 * qJD(1);
t160 = -t230 * pkin(3) + qJDD(3) * pkin(7) + t204 * t245 + t163;
t223 = sin(qJ(4));
t227 = cos(qJ(4));
t143 = t227 * t156 - t223 * t160;
t246 = qJD(1) * t224;
t199 = t227 * qJD(3) - t223 * t246;
t176 = t199 * qJD(4) + t223 * qJDD(3) + t227 * t205;
t198 = qJDD(4) - t206;
t200 = t223 * qJD(3) + t227 * t246;
t213 = qJD(4) - t245;
t140 = (t199 * t213 - t176) * pkin(8) + (t199 * t200 + t198) * pkin(4) + t143;
t144 = t223 * t156 + t227 * t160;
t175 = -t200 * qJD(4) + t227 * qJDD(3) - t223 * t205;
t185 = t213 * pkin(4) - t200 * pkin(8);
t197 = t199 ^ 2;
t141 = -t197 * pkin(4) + t175 * pkin(8) - t213 * t185 + t144;
t222 = sin(qJ(5));
t226 = cos(qJ(5));
t139 = t222 * t140 + t226 * t141;
t162 = -t224 * t169 + t228 * t219;
t159 = -qJDD(3) * pkin(3) - t230 * pkin(7) + t204 * t246 - t162;
t142 = -t175 * pkin(4) - t197 * pkin(8) + t200 * t185 + t159;
t179 = t222 * t199 + t226 * t200;
t149 = -t179 * qJD(5) + t226 * t175 - t222 * t176;
t178 = t226 * t199 - t222 * t200;
t150 = t178 * qJD(5) + t222 * t175 + t226 * t176;
t212 = qJD(5) + t213;
t153 = Ifges(6,5) * t179 + Ifges(6,6) * t178 + Ifges(6,3) * t212;
t155 = Ifges(6,1) * t179 + Ifges(6,4) * t178 + Ifges(6,5) * t212;
t194 = qJDD(5) + t198;
t127 = -mrSges(6,1) * t142 + mrSges(6,3) * t139 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * t194 - t179 * t153 + t212 * t155;
t138 = t226 * t140 - t222 * t141;
t154 = Ifges(6,4) * t179 + Ifges(6,2) * t178 + Ifges(6,6) * t212;
t128 = mrSges(6,2) * t142 - mrSges(6,3) * t138 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * t194 + t178 * t153 - t212 * t154;
t170 = Ifges(5,5) * t200 + Ifges(5,6) * t199 + Ifges(5,3) * t213;
t172 = Ifges(5,1) * t200 + Ifges(5,4) * t199 + Ifges(5,5) * t213;
t164 = -t212 * mrSges(6,2) + t178 * mrSges(6,3);
t165 = t212 * mrSges(6,1) - t179 * mrSges(6,3);
t237 = m(6) * t142 - t149 * mrSges(6,1) + t150 * mrSges(6,2) - t178 * t164 + t179 * t165;
t161 = -t178 * mrSges(6,1) + t179 * mrSges(6,2);
t134 = m(6) * t138 + t194 * mrSges(6,1) - t150 * mrSges(6,3) - t179 * t161 + t212 * t164;
t135 = m(6) * t139 - t194 * mrSges(6,2) + t149 * mrSges(6,3) + t178 * t161 - t212 * t165;
t240 = -t222 * t134 + t226 * t135;
t108 = -mrSges(5,1) * t159 + mrSges(5,3) * t144 + Ifges(5,4) * t176 + Ifges(5,2) * t175 + Ifges(5,6) * t198 - pkin(4) * t237 + pkin(8) * t240 + t226 * t127 + t222 * t128 - t200 * t170 + t213 * t172;
t126 = t226 * t134 + t222 * t135;
t171 = Ifges(5,4) * t200 + Ifges(5,2) * t199 + Ifges(5,6) * t213;
t114 = mrSges(5,2) * t159 - mrSges(5,3) * t143 + Ifges(5,1) * t176 + Ifges(5,4) * t175 + Ifges(5,5) * t198 - pkin(8) * t126 - t222 * t127 + t226 * t128 + t199 * t170 - t213 * t171;
t182 = -t199 * mrSges(5,1) + t200 * mrSges(5,2);
t183 = -t213 * mrSges(5,2) + t199 * mrSges(5,3);
t124 = m(5) * t143 + t198 * mrSges(5,1) - t176 * mrSges(5,3) - t200 * t182 + t213 * t183 + t126;
t184 = t213 * mrSges(5,1) - t200 * mrSges(5,3);
t125 = m(5) * t144 - t198 * mrSges(5,2) + t175 * mrSges(5,3) + t199 * t182 - t213 * t184 + t240;
t122 = -t223 * t124 + t227 * t125;
t136 = -m(5) * t159 + t175 * mrSges(5,1) - t176 * mrSges(5,2) + t199 * t183 - t200 * t184 - t237;
t192 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t224 + Ifges(4,2) * t228) * qJD(1);
t193 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t224 + Ifges(4,4) * t228) * qJD(1);
t247 = mrSges(4,1) * t162 - mrSges(4,2) * t163 + Ifges(4,5) * t205 + Ifges(4,6) * t206 + Ifges(4,3) * qJDD(3) + pkin(3) * t136 + pkin(7) * t122 + t227 * t108 + t223 * t114 + (t224 * t192 - t228 * t193) * qJD(1);
t202 = (-mrSges(4,1) * t228 + mrSges(4,2) * t224) * qJD(1);
t208 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t246;
t120 = m(4) * t163 - qJDD(3) * mrSges(4,2) + t206 * mrSges(4,3) - qJD(3) * t208 + t202 * t245 + t122;
t209 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t245;
t130 = m(4) * t162 + qJDD(3) * mrSges(4,1) - t205 * mrSges(4,3) + qJD(3) * t209 - t202 * t246 + t136;
t241 = t228 * t120 - t224 * t130;
t111 = m(3) * t181 - t231 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t241;
t121 = t227 * t124 + t223 * t125;
t234 = -m(4) * t168 + t206 * mrSges(4,1) - t205 * mrSges(4,2) - t208 * t246 + t209 * t245 - t121;
t116 = m(3) * t180 + qJDD(1) * mrSges(3,1) - t231 * mrSges(3,2) + t234;
t105 = t220 * t111 + t221 * t116;
t113 = t224 * t120 + t228 * t130;
t242 = t221 * t111 - t220 * t116;
t191 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t224 + Ifges(4,6) * t228) * qJD(1);
t101 = mrSges(4,2) * t168 - mrSges(4,3) * t162 + Ifges(4,1) * t205 + Ifges(4,4) * t206 + Ifges(4,5) * qJDD(3) - pkin(7) * t121 - qJD(3) * t192 - t223 * t108 + t227 * t114 + t191 * t245;
t236 = -mrSges(6,1) * t138 + mrSges(6,2) * t139 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * t194 - t179 * t154 + t178 * t155;
t232 = mrSges(5,1) * t143 - mrSges(5,2) * t144 + Ifges(5,5) * t176 + Ifges(5,6) * t175 + Ifges(5,3) * t198 + pkin(4) * t126 + t200 * t171 - t199 * t172 - t236;
t107 = -mrSges(4,1) * t168 + mrSges(4,3) * t163 + Ifges(4,4) * t205 + Ifges(4,2) * t206 + Ifges(4,6) * qJDD(3) - pkin(3) * t121 + qJD(3) * t193 - t191 * t246 - t232;
t238 = mrSges(3,1) * t180 - mrSges(3,2) * t181 + Ifges(3,3) * qJDD(1) + pkin(2) * t234 + pkin(6) * t241 + t224 * t101 + t228 * t107;
t235 = mrSges(2,1) * t210 - mrSges(2,2) * t211 + Ifges(2,3) * qJDD(1) + pkin(1) * t105 + t238;
t103 = m(2) * t211 - t231 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t242;
t102 = m(2) * t210 + qJDD(1) * mrSges(2,1) - t231 * mrSges(2,2) + t105;
t99 = -mrSges(3,1) * t219 + mrSges(3,3) * t181 + t231 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t113 - t247;
t98 = mrSges(3,2) * t219 - mrSges(3,3) * t180 + Ifges(3,5) * qJDD(1) - t231 * Ifges(3,6) - pkin(6) * t113 + t228 * t101 - t224 * t107;
t97 = -mrSges(2,2) * g(3) - mrSges(2,3) * t210 + Ifges(2,5) * qJDD(1) - t231 * Ifges(2,6) - qJ(2) * t105 - t220 * t99 + t221 * t98;
t96 = Ifges(2,6) * qJDD(1) + t231 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t211 + t220 * t98 + t221 * t99 - pkin(1) * (m(3) * t219 + t113) + qJ(2) * t242;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t229 * t97 - t225 * t96 - pkin(5) * (t229 * t102 + t225 * t103), t97, t98, t101, t114, t128; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t225 * t97 + t229 * t96 + pkin(5) * (-t225 * t102 + t229 * t103), t96, t99, t107, t108, t127; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t235, t235, t238, t247, t232, -t236;];
m_new = t1;
