% Calculate vector of cutting torques with Newton-Euler for
% S5PRRPP2
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRPP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:08:49
% EndTime: 2019-12-05 16:08:56
% DurationCPUTime: 2.86s
% Computational Cost: add. (30758->254), mult. (65856->313), div. (0->0), fcn. (39890->8), ass. (0->97)
t219 = sin(pkin(7));
t248 = cos(pkin(7));
t204 = -t248 * g(1) - t219 * g(2);
t217 = -g(3) + qJDD(1);
t221 = sin(qJ(2));
t223 = cos(qJ(2));
t184 = t223 * t204 + t221 * t217;
t225 = qJD(2) ^ 2;
t176 = -t225 * pkin(2) + qJDD(2) * pkin(6) + t184;
t203 = t219 * g(1) - t248 * g(2);
t220 = sin(qJ(3));
t222 = cos(qJ(3));
t148 = -t220 * t176 - t222 * t203;
t242 = qJD(2) * qJD(3);
t240 = t222 * t242;
t200 = t220 * qJDD(2) + t240;
t143 = (-t200 + t240) * qJ(4) + (t220 * t222 * t225 + qJDD(3)) * pkin(3) + t148;
t149 = t222 * t176 - t220 * t203;
t201 = t222 * qJDD(2) - t220 * t242;
t244 = qJD(2) * t220;
t205 = qJD(3) * pkin(3) - qJ(4) * t244;
t216 = t222 ^ 2;
t144 = -t216 * t225 * pkin(3) + t201 * qJ(4) - qJD(3) * t205 + t149;
t218 = sin(pkin(8));
t243 = qJD(2) * t222;
t247 = cos(pkin(8));
t186 = t218 * t244 - t247 * t243;
t250 = -2 * qJD(4);
t140 = t218 * t143 + t247 * t144 + t186 * t250;
t172 = t218 * t200 - t247 * t201;
t187 = (t218 * t222 + t247 * t220) * qJD(2);
t180 = qJD(3) * mrSges(5,1) - t187 * mrSges(5,3);
t160 = t186 * pkin(4) - t187 * qJ(5);
t224 = qJD(3) ^ 2;
t133 = -t224 * pkin(4) + qJDD(3) * qJ(5) + 0.2e1 * qJD(5) * qJD(3) - t186 * t160 + t140;
t181 = -qJD(3) * mrSges(6,1) + t187 * mrSges(6,2);
t241 = m(6) * t133 + qJDD(3) * mrSges(6,3) + qJD(3) * t181;
t161 = t186 * mrSges(6,1) - t187 * mrSges(6,3);
t245 = -t186 * mrSges(5,1) - t187 * mrSges(5,2) - t161;
t249 = -mrSges(5,3) - mrSges(6,2);
t123 = m(5) * t140 - qJDD(3) * mrSges(5,2) - qJD(3) * t180 + t249 * t172 + t245 * t186 + t241;
t233 = t247 * t143 - t218 * t144;
t139 = t187 * t250 + t233;
t173 = t247 * t200 + t218 * t201;
t179 = -qJD(3) * mrSges(5,2) - t186 * mrSges(5,3);
t135 = -qJDD(3) * pkin(4) - t224 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t160) * t187 - t233;
t182 = -t186 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t237 = -m(6) * t135 + qJDD(3) * mrSges(6,1) + qJD(3) * t182;
t124 = m(5) * t139 + qJDD(3) * mrSges(5,1) + qJD(3) * t179 + t249 * t173 + t245 * t187 + t237;
t118 = t218 * t123 + t247 * t124;
t189 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t220 + Ifges(4,2) * t222) * qJD(2);
t190 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t220 + Ifges(4,4) * t222) * qJD(2);
t154 = Ifges(5,4) * t187 - Ifges(5,2) * t186 + Ifges(5,6) * qJD(3);
t156 = Ifges(5,1) * t187 - Ifges(5,4) * t186 + Ifges(5,5) * qJD(3);
t151 = Ifges(6,5) * t187 + Ifges(6,6) * qJD(3) + Ifges(6,3) * t186;
t155 = Ifges(6,1) * t187 + Ifges(6,4) * qJD(3) + Ifges(6,5) * t186;
t230 = mrSges(6,1) * t135 - mrSges(6,3) * t133 - Ifges(6,4) * t173 - Ifges(6,2) * qJDD(3) - Ifges(6,6) * t172 + t187 * t151 - t186 * t155;
t227 = mrSges(5,2) * t140 - t186 * t156 - qJ(5) * (-t172 * mrSges(6,2) - t186 * t161 + t241) - pkin(4) * (-t173 * mrSges(6,2) - t187 * t161 + t237) - mrSges(5,1) * t139 - t187 * t154 + Ifges(5,6) * t172 - Ifges(5,5) * t173 - Ifges(5,3) * qJDD(3) + t230;
t251 = mrSges(4,1) * t148 - mrSges(4,2) * t149 + Ifges(4,5) * t200 + Ifges(4,6) * t201 + Ifges(4,3) * qJDD(3) + pkin(3) * t118 + (t220 * t189 - t222 * t190) * qJD(2) - t227;
t153 = Ifges(6,4) * t187 + Ifges(6,2) * qJD(3) + Ifges(6,6) * t186;
t246 = -Ifges(5,5) * t187 + Ifges(5,6) * t186 - Ifges(5,3) * qJD(3) - t153;
t199 = (-mrSges(4,1) * t222 + mrSges(4,2) * t220) * qJD(2);
t207 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t243;
t114 = m(4) * t148 + qJDD(3) * mrSges(4,1) - t200 * mrSges(4,3) + qJD(3) * t207 - t199 * t244 + t118;
t206 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t244;
t238 = t247 * t123 - t218 * t124;
t115 = m(4) * t149 - qJDD(3) * mrSges(4,2) + t201 * mrSges(4,3) - qJD(3) * t206 + t199 * t243 + t238;
t112 = -t220 * t114 + t222 * t115;
t108 = m(3) * t184 - t225 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t112;
t183 = -t221 * t204 + t223 * t217;
t234 = -qJDD(2) * pkin(2) - t183;
t175 = -t225 * pkin(6) + t234;
t145 = -t201 * pkin(3) + qJDD(4) + t205 * t244 + (-qJ(4) * t216 - pkin(6)) * t225 + t234;
t137 = -0.2e1 * qJD(5) * t187 + (qJD(3) * t186 - t173) * qJ(5) + (qJD(3) * t187 + t172) * pkin(4) + t145;
t130 = m(6) * t137 + t172 * mrSges(6,1) - t173 * mrSges(6,3) - t187 * t181 + t186 * t182;
t229 = m(5) * t145 + t172 * mrSges(5,1) + t173 * mrSges(5,2) + t186 * t179 + t187 * t180 + t130;
t125 = -m(4) * t175 + t201 * mrSges(4,1) - t200 * mrSges(4,2) - t206 * t244 + t207 * t243 - t229;
t122 = m(3) * t183 + qJDD(2) * mrSges(3,1) - t225 * mrSges(3,2) + t125;
t239 = t223 * t108 - t221 * t122;
t236 = -mrSges(6,1) * t137 + mrSges(6,2) * t133;
t111 = t222 * t114 + t220 * t115;
t101 = mrSges(3,1) * t203 + mrSges(3,3) * t184 + t225 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t111 - t251;
t116 = -mrSges(5,1) * t145 + mrSges(5,3) * t140 - pkin(4) * t130 + t246 * t187 + (Ifges(5,4) - Ifges(6,5)) * t173 + (-Ifges(5,2) - Ifges(6,3)) * t172 + (Ifges(5,6) - Ifges(6,6)) * qJDD(3) + (t155 + t156) * qJD(3) + t236;
t231 = mrSges(6,2) * t135 - mrSges(6,3) * t137 + Ifges(6,1) * t173 + Ifges(6,4) * qJDD(3) + Ifges(6,5) * t172 + qJD(3) * t151;
t117 = mrSges(5,2) * t145 - mrSges(5,3) * t139 + Ifges(5,1) * t173 - Ifges(5,4) * t172 + Ifges(5,5) * qJDD(3) - qJ(5) * t130 - qJD(3) * t154 + t246 * t186 + t231;
t188 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t220 + Ifges(4,6) * t222) * qJD(2);
t102 = -mrSges(4,1) * t175 + mrSges(4,3) * t149 + Ifges(4,4) * t200 + Ifges(4,2) * t201 + Ifges(4,6) * qJDD(3) - pkin(3) * t229 + qJ(4) * t238 + qJD(3) * t190 + t247 * t116 + t218 * t117 - t188 * t244;
t103 = mrSges(4,2) * t175 - mrSges(4,3) * t148 + Ifges(4,1) * t200 + Ifges(4,4) * t201 + Ifges(4,5) * qJDD(3) - qJ(4) * t118 - qJD(3) * t189 - t218 * t116 + t247 * t117 + t188 * t243;
t99 = -mrSges(3,2) * t203 - mrSges(3,3) * t183 + Ifges(3,5) * qJDD(2) - t225 * Ifges(3,6) - pkin(6) * t111 - t220 * t102 + t222 * t103;
t232 = -mrSges(2,2) * t204 + pkin(5) * t239 + t223 * t101 + t221 * t99 + pkin(1) * (m(3) * t203 - t111) + mrSges(2,1) * t203;
t228 = mrSges(3,1) * t183 - mrSges(3,2) * t184 + Ifges(3,3) * qJDD(2) + pkin(2) * t125 + pkin(6) * t112 + t222 * t102 + t220 * t103;
t109 = (m(2) + m(3)) * t203 - t111;
t106 = t221 * t108 + t223 * t122;
t104 = m(2) * t204 + t239;
t97 = -mrSges(2,1) * t217 + mrSges(2,3) * t204 - pkin(1) * t106 - t228;
t96 = mrSges(2,2) * t217 - mrSges(2,3) * t203 - pkin(5) * t106 - t221 * t101 + t223 * t99;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t248 * t96 - t219 * t97 - qJ(1) * (t219 * t104 + t248 * t109), t96, t99, t103, t117, -t186 * t153 + t231; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t219 * t96 + t248 * t97 + qJ(1) * (t248 * t104 - t219 * t109), t97, t101, t102, t116, -t230; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t232, t232, t228, t251, -t227, Ifges(6,5) * t173 + Ifges(6,6) * qJDD(3) + Ifges(6,3) * t172 - qJD(3) * t155 + t187 * t153 - t236;];
m_new = t1;
