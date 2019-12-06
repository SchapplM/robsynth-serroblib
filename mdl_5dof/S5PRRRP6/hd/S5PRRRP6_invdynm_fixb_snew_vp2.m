% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:52
% EndTime: 2019-12-05 16:50:58
% DurationCPUTime: 2.91s
% Computational Cost: add. (34239->254), mult. (67976->311), div. (0->0), fcn. (41986->8), ass. (0->97)
t220 = sin(pkin(8));
t247 = cos(pkin(8));
t203 = -t247 * g(1) - t220 * g(2);
t219 = -g(3) + qJDD(1);
t223 = sin(qJ(2));
t225 = cos(qJ(2));
t184 = t225 * t203 + t223 * t219;
t226 = qJD(2) ^ 2;
t178 = -t226 * pkin(2) + qJDD(2) * pkin(6) + t184;
t202 = t220 * g(1) - t247 * g(2);
t222 = sin(qJ(3));
t224 = cos(qJ(3));
t160 = -t222 * t178 - t224 * t202;
t242 = qJD(2) * qJD(3);
t240 = t224 * t242;
t199 = t222 * qJDD(2) + t240;
t143 = (-t199 + t240) * pkin(7) + (t222 * t224 * t226 + qJDD(3)) * pkin(3) + t160;
t161 = t224 * t178 - t222 * t202;
t200 = t224 * qJDD(2) - t222 * t242;
t244 = qJD(2) * t222;
t206 = qJD(3) * pkin(3) - pkin(7) * t244;
t218 = t224 ^ 2;
t144 = -t218 * t226 * pkin(3) + t200 * pkin(7) - qJD(3) * t206 + t161;
t221 = sin(qJ(4));
t249 = cos(qJ(4));
t140 = t221 * t143 + t249 * t144;
t189 = (t221 * t224 + t249 * t222) * qJD(2);
t157 = t189 * qJD(4) + t221 * t199 - t249 * t200;
t217 = qJD(3) + qJD(4);
t180 = t217 * mrSges(5,1) - t189 * mrSges(5,3);
t243 = qJD(2) * t224;
t188 = t221 * t244 - t249 * t243;
t216 = qJDD(3) + qJDD(4);
t171 = t188 * pkin(4) - t189 * qJ(5);
t215 = t217 ^ 2;
t133 = -t215 * pkin(4) + t216 * qJ(5) + 0.2e1 * qJD(5) * t217 - t188 * t171 + t140;
t181 = -t217 * mrSges(6,1) + t189 * mrSges(6,2);
t241 = m(6) * t133 + t216 * mrSges(6,3) + t217 * t181;
t172 = t188 * mrSges(6,1) - t189 * mrSges(6,3);
t245 = -t188 * mrSges(5,1) - t189 * mrSges(5,2) - t172;
t248 = -mrSges(5,3) - mrSges(6,2);
t123 = m(5) * t140 - t216 * mrSges(5,2) + t248 * t157 - t217 * t180 + t245 * t188 + t241;
t139 = t249 * t143 - t221 * t144;
t158 = -t188 * qJD(4) + t249 * t199 + t221 * t200;
t179 = -t217 * mrSges(5,2) - t188 * mrSges(5,3);
t135 = -t216 * pkin(4) - t215 * qJ(5) + t189 * t171 + qJDD(5) - t139;
t182 = -t188 * mrSges(6,2) + t217 * mrSges(6,3);
t237 = -m(6) * t135 + t216 * mrSges(6,1) + t217 * t182;
t125 = m(5) * t139 + t216 * mrSges(5,1) + t248 * t158 + t217 * t179 + t245 * t189 + t237;
t118 = t221 * t123 + t249 * t125;
t186 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t222 + Ifges(4,2) * t224) * qJD(2);
t187 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t222 + Ifges(4,4) * t224) * qJD(2);
t165 = Ifges(5,4) * t189 - Ifges(5,2) * t188 + Ifges(5,6) * t217;
t167 = Ifges(5,1) * t189 - Ifges(5,4) * t188 + Ifges(5,5) * t217;
t162 = Ifges(6,5) * t189 + Ifges(6,6) * t217 + Ifges(6,3) * t188;
t166 = Ifges(6,1) * t189 + Ifges(6,4) * t217 + Ifges(6,5) * t188;
t231 = mrSges(6,1) * t135 - mrSges(6,3) * t133 - Ifges(6,4) * t158 - Ifges(6,2) * t216 - Ifges(6,6) * t157 + t189 * t162 - t188 * t166;
t228 = mrSges(5,2) * t140 - t188 * t167 - qJ(5) * (-t157 * mrSges(6,2) - t188 * t172 + t241) - pkin(4) * (-t158 * mrSges(6,2) - t189 * t172 + t237) - mrSges(5,1) * t139 - t189 * t165 + Ifges(5,6) * t157 - Ifges(5,5) * t158 - Ifges(5,3) * t216 + t231;
t250 = mrSges(4,1) * t160 - mrSges(4,2) * t161 + Ifges(4,5) * t199 + Ifges(4,6) * t200 + Ifges(4,3) * qJDD(3) + pkin(3) * t118 + (t222 * t186 - t224 * t187) * qJD(2) - t228;
t164 = Ifges(6,4) * t189 + Ifges(6,2) * t217 + Ifges(6,6) * t188;
t246 = -Ifges(5,5) * t189 + Ifges(5,6) * t188 - Ifges(5,3) * t217 - t164;
t198 = (-mrSges(4,1) * t224 + mrSges(4,2) * t222) * qJD(2);
t205 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t243;
t116 = m(4) * t160 + qJDD(3) * mrSges(4,1) - t199 * mrSges(4,3) + qJD(3) * t205 - t198 * t244 + t118;
t204 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t244;
t238 = t249 * t123 - t221 * t125;
t117 = m(4) * t161 - qJDD(3) * mrSges(4,2) + t200 * mrSges(4,3) - qJD(3) * t204 + t198 * t243 + t238;
t112 = -t222 * t116 + t224 * t117;
t108 = m(3) * t184 - t226 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t112;
t183 = -t223 * t203 + t225 * t219;
t234 = -qJDD(2) * pkin(2) - t183;
t177 = -t226 * pkin(6) + t234;
t145 = -t200 * pkin(3) + t206 * t244 + (-pkin(7) * t218 - pkin(6)) * t226 + t234;
t137 = -0.2e1 * qJD(5) * t189 + (t188 * t217 - t158) * qJ(5) + (t189 * t217 + t157) * pkin(4) + t145;
t126 = m(6) * t137 + t157 * mrSges(6,1) - t158 * mrSges(6,3) - t189 * t181 + t188 * t182;
t230 = m(5) * t145 + t157 * mrSges(5,1) + t158 * mrSges(5,2) + t188 * t179 + t189 * t180 + t126;
t120 = -m(4) * t177 + t200 * mrSges(4,1) - t199 * mrSges(4,2) - t204 * t244 + t205 * t243 - t230;
t119 = m(3) * t183 + qJDD(2) * mrSges(3,1) - t226 * mrSges(3,2) + t120;
t239 = t225 * t108 - t223 * t119;
t236 = -mrSges(6,1) * t137 + mrSges(6,2) * t133;
t111 = t224 * t116 + t222 * t117;
t101 = mrSges(3,1) * t202 + mrSges(3,3) * t184 + t226 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t111 - t250;
t113 = -mrSges(5,1) * t145 + mrSges(5,3) * t140 - pkin(4) * t126 + (t166 + t167) * t217 + (Ifges(5,6) - Ifges(6,6)) * t216 + t246 * t189 + (Ifges(5,4) - Ifges(6,5)) * t158 + (-Ifges(5,2) - Ifges(6,3)) * t157 + t236;
t232 = mrSges(6,2) * t135 - mrSges(6,3) * t137 + Ifges(6,1) * t158 + Ifges(6,4) * t216 + Ifges(6,5) * t157 + t217 * t162;
t114 = mrSges(5,2) * t145 - mrSges(5,3) * t139 + Ifges(5,1) * t158 - Ifges(5,4) * t157 + Ifges(5,5) * t216 - qJ(5) * t126 - t217 * t165 + t246 * t188 + t232;
t185 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t222 + Ifges(4,6) * t224) * qJD(2);
t102 = -mrSges(4,1) * t177 + mrSges(4,3) * t161 + Ifges(4,4) * t199 + Ifges(4,2) * t200 + Ifges(4,6) * qJDD(3) - pkin(3) * t230 + pkin(7) * t238 + qJD(3) * t187 + t249 * t113 + t221 * t114 - t185 * t244;
t103 = mrSges(4,2) * t177 - mrSges(4,3) * t160 + Ifges(4,1) * t199 + Ifges(4,4) * t200 + Ifges(4,5) * qJDD(3) - pkin(7) * t118 - qJD(3) * t186 - t221 * t113 + t249 * t114 + t185 * t243;
t99 = -mrSges(3,2) * t202 - mrSges(3,3) * t183 + Ifges(3,5) * qJDD(2) - t226 * Ifges(3,6) - pkin(6) * t111 - t222 * t102 + t224 * t103;
t233 = -mrSges(2,2) * t203 + pkin(5) * t239 + t225 * t101 + t223 * t99 + pkin(1) * (m(3) * t202 - t111) + mrSges(2,1) * t202;
t229 = mrSges(3,1) * t183 - mrSges(3,2) * t184 + Ifges(3,3) * qJDD(2) + pkin(2) * t120 + pkin(6) * t112 + t224 * t102 + t222 * t103;
t109 = (m(2) + m(3)) * t202 - t111;
t106 = t223 * t108 + t225 * t119;
t104 = m(2) * t203 + t239;
t97 = -mrSges(2,1) * t219 + mrSges(2,3) * t203 - pkin(1) * t106 - t229;
t96 = mrSges(2,2) * t219 - mrSges(2,3) * t202 - pkin(5) * t106 - t223 * t101 + t225 * t99;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t247 * t96 - t220 * t97 - qJ(1) * (t220 * t104 + t247 * t109), t96, t99, t103, t114, -t188 * t164 + t232; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t220 * t96 + t247 * t97 + qJ(1) * (t247 * t104 - t220 * t109), t97, t101, t102, t113, -t231; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t233, t233, t229, t250, -t228, Ifges(6,5) * t158 + Ifges(6,6) * t216 + Ifges(6,3) * t157 + t189 * t164 - t217 * t166 - t236;];
m_new = t1;
