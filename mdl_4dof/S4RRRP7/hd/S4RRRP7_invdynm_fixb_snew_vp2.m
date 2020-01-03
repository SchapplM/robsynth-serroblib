% Calculate vector of cutting torques with Newton-Euler for
% S4RRRP7
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% m [3x5]
%   vector of cutting torques (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP7_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP7_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:20:11
% EndTime: 2019-12-31 17:20:15
% DurationCPUTime: 1.37s
% Computational Cost: add. (12963->231), mult. (25221->282), div. (0->0), fcn. (14391->6), ass. (0->85)
t185 = sin(qJ(1));
t187 = cos(qJ(1));
t174 = t185 * g(1) - t187 * g(2);
t189 = qJD(1) ^ 2;
t156 = -qJDD(1) * pkin(1) - t189 * pkin(5) - t174;
t184 = sin(qJ(2));
t186 = cos(qJ(2));
t203 = qJD(1) * qJD(2);
t200 = t186 * t203;
t169 = t184 * qJDD(1) + t200;
t201 = t184 * t203;
t170 = t186 * qJDD(1) - t201;
t117 = (-t169 - t200) * pkin(6) + (-t170 + t201) * pkin(2) + t156;
t175 = -t187 * g(1) - t185 * g(2);
t157 = -t189 * pkin(1) + qJDD(1) * pkin(5) + t175;
t149 = -t184 * g(3) + t186 * t157;
t168 = (-pkin(2) * t186 - pkin(6) * t184) * qJD(1);
t188 = qJD(2) ^ 2;
t204 = t186 * qJD(1);
t121 = -t188 * pkin(2) + qJDD(2) * pkin(6) + t168 * t204 + t149;
t183 = sin(qJ(3));
t210 = cos(qJ(3));
t114 = t210 * t117 - t183 * t121;
t115 = t183 * t117 + t210 * t121;
t205 = qJD(1) * t184;
t165 = -t210 * qJD(2) + t183 * t205;
t166 = t183 * qJD(2) + t210 * t205;
t177 = qJD(3) - t204;
t123 = Ifges(5,5) * t166 + Ifges(5,6) * t177 + Ifges(5,3) * t165;
t126 = Ifges(4,4) * t166 - Ifges(4,2) * t165 + Ifges(4,6) * t177;
t128 = Ifges(4,1) * t166 - Ifges(4,4) * t165 + Ifges(4,5) * t177;
t137 = t166 * qJD(3) - t210 * qJDD(2) + t183 * t169;
t138 = -t165 * qJD(3) + t183 * qJDD(2) + t210 * t169;
t142 = t165 * mrSges(5,1) - t166 * mrSges(5,3);
t164 = qJDD(3) - t170;
t141 = t165 * pkin(3) - t166 * qJ(4);
t176 = t177 ^ 2;
t110 = -t176 * pkin(3) + t164 * qJ(4) + 0.2e1 * qJD(4) * t177 - t165 * t141 + t115;
t113 = -t164 * pkin(3) - t176 * qJ(4) + t166 * t141 + qJDD(4) - t114;
t127 = Ifges(5,1) * t166 + Ifges(5,4) * t177 + Ifges(5,5) * t165;
t195 = mrSges(5,1) * t113 - mrSges(5,3) * t110 - Ifges(5,4) * t138 - Ifges(5,2) * t164 - Ifges(5,6) * t137 - t165 * t127;
t147 = -t165 * mrSges(5,2) + t177 * mrSges(5,3);
t198 = -m(5) * t113 + t164 * mrSges(5,1) + t177 * t147;
t146 = -t177 * mrSges(5,1) + t166 * mrSges(5,2);
t202 = m(5) * t110 + t164 * mrSges(5,3) + t177 * t146;
t212 = -(-t126 + t123) * t166 + mrSges(4,1) * t114 - mrSges(4,2) * t115 + Ifges(4,5) * t138 - Ifges(4,6) * t137 + Ifges(4,3) * t164 + pkin(3) * (-t138 * mrSges(5,2) - t166 * t142 + t198) + qJ(4) * (-t137 * mrSges(5,2) - t165 * t142 + t202) + t165 * t128 - t195;
t148 = -t186 * g(3) - t184 * t157;
t120 = -qJDD(2) * pkin(2) - t188 * pkin(6) + t168 * t205 - t148;
t112 = -0.2e1 * qJD(4) * t166 + (t165 * t177 - t138) * qJ(4) + (t166 * t177 + t137) * pkin(3) + t120;
t105 = m(5) * t112 + t137 * mrSges(5,1) - t138 * mrSges(5,3) - t166 * t146 + t165 * t147;
t144 = -t177 * mrSges(4,2) - t165 * mrSges(4,3);
t145 = t177 * mrSges(4,1) - t166 * mrSges(4,3);
t104 = -m(4) * t120 - t137 * mrSges(4,1) - t138 * mrSges(4,2) - t165 * t144 - t166 * t145 - t105;
t154 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t184 + Ifges(3,2) * t186) * qJD(1);
t155 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t184 + Ifges(3,4) * t186) * qJD(1);
t197 = -mrSges(5,1) * t112 + mrSges(5,2) * t110;
t125 = Ifges(5,4) * t166 + Ifges(5,2) * t177 + Ifges(5,6) * t165;
t208 = -Ifges(4,5) * t166 + Ifges(4,6) * t165 - Ifges(4,3) * t177 - t125;
t92 = -mrSges(4,1) * t120 + mrSges(4,3) * t115 - pkin(3) * t105 + (t127 + t128) * t177 + t208 * t166 + (Ifges(4,6) - Ifges(5,6)) * t164 + (Ifges(4,4) - Ifges(5,5)) * t138 + (-Ifges(4,2) - Ifges(5,3)) * t137 + t197;
t194 = mrSges(5,2) * t113 - mrSges(5,3) * t112 + Ifges(5,1) * t138 + Ifges(5,4) * t164 + Ifges(5,5) * t137 + t177 * t123;
t93 = mrSges(4,2) * t120 - mrSges(4,3) * t114 + Ifges(4,1) * t138 - Ifges(4,4) * t137 + Ifges(4,5) * t164 - qJ(4) * t105 - t177 * t126 + t208 * t165 + t194;
t206 = -t165 * mrSges(4,1) - t166 * mrSges(4,2) - t142;
t209 = -mrSges(4,3) - mrSges(5,2);
t102 = m(4) * t115 - t164 * mrSges(4,2) + t209 * t137 - t177 * t145 + t206 * t165 + t202;
t103 = m(4) * t114 + t164 * mrSges(4,1) + t209 * t138 + t177 * t144 + t206 * t166 + t198;
t99 = t210 * t102 - t183 * t103;
t211 = mrSges(3,1) * t148 - mrSges(3,2) * t149 + Ifges(3,5) * t169 + Ifges(3,6) * t170 + Ifges(3,3) * qJDD(2) + pkin(2) * t104 + pkin(6) * t99 + (t154 * t184 - t155 * t186) * qJD(1) + t183 * t93 + t210 * t92;
t167 = (-mrSges(3,1) * t186 + mrSges(3,2) * t184) * qJD(1);
t173 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t204;
t100 = m(3) * t148 + qJDD(2) * mrSges(3,1) - t169 * mrSges(3,3) + qJD(2) * t173 - t167 * t205 + t104;
t172 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t205;
t97 = m(3) * t149 - qJDD(2) * mrSges(3,2) + t170 * mrSges(3,3) - qJD(2) * t172 + t167 * t204 + t99;
t199 = -t184 * t100 + t186 * t97;
t98 = t183 * t102 + t210 * t103;
t192 = -m(3) * t156 + t170 * mrSges(3,1) - t169 * mrSges(3,2) - t172 * t205 + t173 * t204 - t98;
t153 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t184 + Ifges(3,6) * t186) * qJD(1);
t86 = mrSges(3,2) * t156 - mrSges(3,3) * t148 + Ifges(3,1) * t169 + Ifges(3,4) * t170 + Ifges(3,5) * qJDD(2) - pkin(6) * t98 - qJD(2) * t154 + t153 * t204 - t183 * t92 + t210 * t93;
t88 = -mrSges(3,1) * t156 + mrSges(3,3) * t149 + Ifges(3,4) * t169 + Ifges(3,2) * t170 + Ifges(3,6) * qJDD(2) - pkin(2) * t98 + qJD(2) * t155 - t153 * t205 - t212;
t193 = mrSges(2,1) * t174 - mrSges(2,2) * t175 + Ifges(2,3) * qJDD(1) + pkin(1) * t192 + pkin(5) * t199 + t184 * t86 + t186 * t88;
t94 = m(2) * t174 + qJDD(1) * mrSges(2,1) - t189 * mrSges(2,2) + t192;
t91 = t186 * t100 + t184 * t97;
t89 = m(2) * t175 - t189 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t199;
t84 = mrSges(2,1) * g(3) + mrSges(2,3) * t175 + t189 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t91 - t211;
t83 = -mrSges(2,2) * g(3) - mrSges(2,3) * t174 + Ifges(2,5) * qJDD(1) - t189 * Ifges(2,6) - pkin(5) * t91 - t184 * t88 + t186 * t86;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t187 * t83 - t185 * t84 - pkin(4) * (t185 * t89 + t187 * t94), t83, t86, t93, -t165 * t125 + t194; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t185 * t83 + t187 * t84 + pkin(4) * (-t185 * t94 + t187 * t89), t84, t88, t92, -t166 * t123 - t195; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t193, t193, t211, t212, Ifges(5,5) * t138 + Ifges(5,6) * t164 + Ifges(5,3) * t137 + t166 * t125 - t177 * t127 - t197;];
m_new = t1;
