% Calculate vector of cutting torques with Newton-Euler for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RPPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:09
% EndTime: 2019-12-31 17:45:11
% DurationCPUTime: 1.90s
% Computational Cost: add. (21394->203), mult. (41603->247), div. (0->0), fcn. (23232->8), ass. (0->94)
t200 = qJD(1) ^ 2;
t192 = sin(pkin(8));
t181 = t192 ^ 2;
t194 = cos(pkin(8));
t229 = t194 ^ 2 + t181;
t221 = t229 * mrSges(5,3);
t238 = t200 * t221;
t197 = sin(qJ(1));
t199 = cos(qJ(1));
t170 = t197 * g(1) - t199 * g(2);
t167 = qJDD(1) * pkin(1) + t170;
t171 = -t199 * g(1) - t197 * g(2);
t168 = -t200 * pkin(1) + t171;
t193 = sin(pkin(7));
t195 = cos(pkin(7));
t152 = t195 * t167 - t193 * t168;
t212 = -t200 * qJ(3) + qJDD(3) - t152;
t232 = -pkin(2) - qJ(4);
t237 = -(2 * qJD(1) * qJD(4)) + qJDD(1) * t232 + t212;
t153 = t193 * t167 + t195 * t168;
t236 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t153;
t235 = pkin(4) * t200;
t234 = mrSges(3,1) - mrSges(4,2);
t233 = -Ifges(3,6) + Ifges(4,5);
t231 = Ifges(5,6) * t192;
t189 = -g(3) + qJDD(2);
t225 = qJDD(1) * t194;
t230 = t237 * t194;
t125 = -pkin(6) * t225 + (-t194 * t235 - t189) * t192 + t230;
t131 = t194 * t189 + t237 * t192;
t226 = qJDD(1) * t192;
t126 = -pkin(6) * t226 - t181 * t235 + t131;
t196 = sin(qJ(5));
t198 = cos(qJ(5));
t123 = t198 * t125 - t196 * t126;
t216 = -t192 * t198 - t194 * t196;
t158 = t216 * qJD(1);
t215 = -t192 * t196 + t194 * t198;
t159 = t215 * qJD(1);
t145 = -t158 * mrSges(6,1) + t159 * mrSges(6,2);
t150 = t158 * qJD(5) + qJDD(1) * t215;
t154 = -qJD(5) * mrSges(6,2) + t158 * mrSges(6,3);
t120 = m(6) * t123 + qJDD(5) * mrSges(6,1) - t150 * mrSges(6,3) + qJD(5) * t154 - t159 * t145;
t124 = t196 * t125 + t198 * t126;
t149 = -t159 * qJD(5) + qJDD(1) * t216;
t155 = qJD(5) * mrSges(6,1) - t159 * mrSges(6,3);
t121 = m(6) * t124 - qJDD(5) * mrSges(6,2) + t149 * mrSges(6,3) - qJD(5) * t155 + t158 * t145;
t109 = t198 * t120 + t196 * t121;
t130 = -t192 * t189 + t230;
t213 = -qJDD(1) * mrSges(5,3) - t200 * (mrSges(5,1) * t192 + mrSges(5,2) * t194);
t106 = m(5) * t130 + t194 * t213 + t109;
t219 = -t196 * t120 + t198 * t121;
t107 = m(5) * t131 + t192 * t213 + t219;
t103 = t194 * t106 + t192 * t107;
t140 = -qJDD(1) * pkin(2) + t212;
t210 = -m(4) * t140 + t200 * mrSges(4,3) - t103;
t100 = m(3) * t152 - t200 * mrSges(3,2) + qJDD(1) * t234 + t210;
t138 = t200 * pkin(2) - t236;
t214 = qJDD(4) + t236;
t136 = t200 * t232 + t214;
t128 = pkin(4) * t226 + (-pkin(6) * t229 + t232) * t200 + t214;
t211 = -m(6) * t128 + t149 * mrSges(6,1) - t150 * mrSges(6,2) + t158 * t154 - t159 * t155;
t207 = -m(5) * t136 - mrSges(5,1) * t226 - mrSges(5,2) * t225 + t211;
t205 = -m(4) * t138 + t200 * mrSges(4,2) + qJDD(1) * mrSges(4,3) - t207;
t114 = t205 - qJDD(1) * mrSges(3,2) + (-mrSges(3,1) - t221) * t200 + m(3) * t153;
t95 = t195 * t100 + t193 * t114;
t228 = t200 * (Ifges(5,5) * t194 - t231);
t222 = Ifges(4,4) + t231;
t220 = -t193 * t100 + t195 * t114;
t104 = -t192 * t106 + t194 * t107;
t102 = m(4) * t189 + t104;
t218 = Ifges(5,1) * t194 - Ifges(5,4) * t192;
t217 = Ifges(5,4) * t194 - Ifges(5,2) * t192;
t142 = Ifges(6,4) * t159 + Ifges(6,2) * t158 + Ifges(6,6) * qJD(5);
t143 = Ifges(6,1) * t159 + Ifges(6,4) * t158 + Ifges(6,5) * qJD(5);
t209 = -mrSges(6,1) * t123 + mrSges(6,2) * t124 - Ifges(6,5) * t150 - Ifges(6,6) * t149 - Ifges(6,3) * qJDD(5) - t159 * t142 + t158 * t143;
t141 = Ifges(6,5) * t159 + Ifges(6,6) * t158 + Ifges(6,3) * qJD(5);
t110 = -mrSges(6,1) * t128 + mrSges(6,3) * t124 + Ifges(6,4) * t150 + Ifges(6,2) * t149 + Ifges(6,6) * qJDD(5) + qJD(5) * t143 - t159 * t141;
t111 = mrSges(6,2) * t128 - mrSges(6,3) * t123 + Ifges(6,1) * t150 + Ifges(6,4) * t149 + Ifges(6,5) * qJDD(5) - qJD(5) * t142 + t158 * t141;
t96 = -mrSges(5,1) * t136 + mrSges(5,3) * t131 + pkin(4) * t211 + pkin(6) * t219 + qJDD(1) * t217 + t198 * t110 + t196 * t111 - t194 * t228;
t98 = mrSges(5,2) * t136 - mrSges(5,3) * t130 - pkin(6) * t109 + qJDD(1) * t218 - t196 * t110 + t198 * t111 - t192 * t228;
t208 = mrSges(4,2) * t140 - mrSges(4,3) * t138 + Ifges(4,1) * qJDD(1) - qJ(4) * t103 - t192 * t96 + t194 * t98;
t206 = -mrSges(4,1) * t138 - pkin(3) * (t207 + t238) - qJ(4) * t104 - t192 * t98 - t194 * t96;
t204 = -mrSges(3,2) * t153 + pkin(2) * (-qJDD(1) * mrSges(4,2) + t210) + qJ(3) * (t205 - t238) + mrSges(3,1) * t152 + Ifges(3,3) * qJDD(1) + t208;
t203 = -mrSges(5,1) * t130 + mrSges(5,2) * t131 - Ifges(5,5) * t225 - pkin(4) * t109 + t209 + (-t192 * t218 - t194 * t217) * t200;
t202 = mrSges(2,1) * t170 - mrSges(2,2) * t171 + Ifges(2,3) * qJDD(1) + pkin(1) * t95 + t204;
t201 = -mrSges(4,1) * t140 - pkin(3) * t103 + t203;
t93 = m(2) * t171 - t200 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t220;
t92 = m(2) * t170 + qJDD(1) * mrSges(2,1) - t200 * mrSges(2,2) + t95;
t91 = -t201 + t233 * t200 + (mrSges(3,2) - mrSges(4,3)) * t189 + (Ifges(3,5) - t222) * qJDD(1) - mrSges(3,3) * t152 - qJ(3) * t102;
t90 = mrSges(3,3) * t153 - pkin(2) * t102 + (-Ifges(4,4) + Ifges(3,5)) * t200 - t234 * t189 - t233 * qJDD(1) + t206;
t89 = -mrSges(2,2) * g(3) - mrSges(2,3) * t170 + Ifges(2,5) * qJDD(1) - t200 * Ifges(2,6) - qJ(2) * t95 - t193 * t90 + t195 * t91;
t88 = Ifges(2,6) * qJDD(1) + t200 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t171 + t193 * t91 + t195 * t90 - pkin(1) * (m(3) * t189 + t102) + qJ(2) * t220;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t199 * t89 - t197 * t88 - pkin(5) * (t197 * t93 + t199 * t92), t89, t91, t208, t98, t111; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t197 * t89 + t199 * t88 + pkin(5) * (-t197 * t92 + t199 * t93), t88, t90, mrSges(4,3) * t189 - t200 * Ifges(4,5) + qJDD(1) * t222 + t201, t96, t110; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t202, t202, t204, -mrSges(4,2) * t189 + t200 * Ifges(4,4) + Ifges(4,5) * qJDD(1) - t206, -Ifges(5,6) * t226 - t203, -t209;];
m_new = t1;
