% Calculate vector of cutting torques with Newton-Euler for
% S4RRRR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:27:26
% EndTime: 2019-12-31 17:27:32
% DurationCPUTime: 2.73s
% Computational Cost: add. (32959->235), mult. (65606->300), div. (0->0), fcn. (41828->8), ass. (0->96)
t189 = sin(qJ(1));
t193 = cos(qJ(1));
t179 = -t193 * g(1) - t189 * g(2);
t195 = qJD(1) ^ 2;
t163 = -t195 * pkin(1) + qJDD(1) * pkin(5) + t179;
t188 = sin(qJ(2));
t192 = cos(qJ(2));
t152 = -t192 * g(3) - t188 * t163;
t172 = (-pkin(2) * t192 - pkin(6) * t188) * qJD(1);
t194 = qJD(2) ^ 2;
t208 = qJD(1) * t188;
t136 = -qJDD(2) * pkin(2) - t194 * pkin(6) + t172 * t208 - t152;
t187 = sin(qJ(3));
t191 = cos(qJ(3));
t170 = t187 * qJD(2) + t191 * t208;
t206 = qJD(1) * qJD(2);
t205 = t192 * t206;
t173 = t188 * qJDD(1) + t205;
t145 = -t170 * qJD(3) + t191 * qJDD(2) - t187 * t173;
t169 = t191 * qJD(2) - t187 * t208;
t146 = t169 * qJD(3) + t187 * qJDD(2) + t191 * t173;
t207 = t192 * qJD(1);
t181 = qJD(3) - t207;
t150 = -t181 * mrSges(4,2) + t169 * mrSges(4,3);
t151 = t181 * mrSges(4,1) - t170 * mrSges(4,3);
t154 = t181 * pkin(3) - t170 * pkin(7);
t167 = t169 ^ 2;
t118 = -t145 * pkin(3) - pkin(7) * t167 + t154 * t170 + t136;
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t148 = t186 * t169 + t190 * t170;
t123 = -t148 * qJD(4) + t145 * t190 - t146 * t186;
t147 = t190 * t169 - t186 * t170;
t124 = t147 * qJD(4) + t145 * t186 + t146 * t190;
t180 = qJD(4) + t181;
t138 = -t180 * mrSges(5,2) + t147 * mrSges(5,3);
t139 = t180 * mrSges(5,1) - t148 * mrSges(5,3);
t200 = m(5) * t118 - t123 * mrSges(5,1) + t124 * mrSges(5,2) - t147 * t138 + t148 * t139;
t107 = -m(4) * t136 + t145 * mrSges(4,1) - t146 * mrSges(4,2) + t169 * t150 - t151 * t170 - t200;
t153 = -t188 * g(3) + t192 * t163;
t160 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t188 + Ifges(3,2) * t192) * qJD(1);
t161 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t188 + Ifges(3,4) * t192) * qJD(1);
t182 = t188 * t206;
t174 = t192 * qJDD(1) - t182;
t178 = t189 * g(1) - t193 * g(2);
t162 = -qJDD(1) * pkin(1) - t195 * pkin(5) - t178;
t134 = (-t173 - t205) * pkin(6) + (-t174 + t182) * pkin(2) + t162;
t137 = -t194 * pkin(2) + qJDD(2) * pkin(6) + t172 * t207 + t153;
t125 = t191 * t134 - t187 * t137;
t168 = qJDD(3) - t174;
t116 = (t169 * t181 - t146) * pkin(7) + (t169 * t170 + t168) * pkin(3) + t125;
t126 = t187 * t134 + t191 * t137;
t117 = -pkin(3) * t167 + t145 * pkin(7) - t154 * t181 + t126;
t115 = t116 * t186 + t117 * t190;
t127 = Ifges(5,5) * t148 + Ifges(5,6) * t147 + Ifges(5,3) * t180;
t129 = Ifges(5,1) * t148 + Ifges(5,4) * t147 + Ifges(5,5) * t180;
t164 = qJDD(4) + t168;
t104 = -mrSges(5,1) * t118 + mrSges(5,3) * t115 + Ifges(5,4) * t124 + Ifges(5,2) * t123 + Ifges(5,6) * t164 - t148 * t127 + t129 * t180;
t114 = t116 * t190 - t117 * t186;
t128 = Ifges(5,4) * t148 + Ifges(5,2) * t147 + Ifges(5,6) * t180;
t105 = mrSges(5,2) * t118 - mrSges(5,3) * t114 + Ifges(5,1) * t124 + Ifges(5,4) * t123 + Ifges(5,5) * t164 + t147 * t127 - t128 * t180;
t140 = Ifges(4,5) * t170 + Ifges(4,6) * t169 + Ifges(4,3) * t181;
t142 = Ifges(4,1) * t170 + Ifges(4,4) * t169 + Ifges(4,5) * t181;
t131 = -mrSges(5,1) * t147 + mrSges(5,2) * t148;
t111 = m(5) * t114 + mrSges(5,1) * t164 - t124 * mrSges(5,3) - t148 * t131 + t138 * t180;
t112 = m(5) * t115 - mrSges(5,2) * t164 + t123 * mrSges(5,3) + t147 * t131 - t139 * t180;
t203 = -t111 * t186 + t190 * t112;
t89 = -mrSges(4,1) * t136 + mrSges(4,3) * t126 + Ifges(4,4) * t146 + Ifges(4,2) * t145 + Ifges(4,6) * t168 - pkin(3) * t200 + pkin(7) * t203 + t190 * t104 + t186 * t105 - t170 * t140 + t181 * t142;
t103 = t190 * t111 + t186 * t112;
t141 = Ifges(4,4) * t170 + Ifges(4,2) * t169 + Ifges(4,6) * t181;
t93 = mrSges(4,2) * t136 - mrSges(4,3) * t125 + Ifges(4,1) * t146 + Ifges(4,4) * t145 + Ifges(4,5) * t168 - pkin(7) * t103 - t104 * t186 + t105 * t190 + t140 * t169 - t141 * t181;
t149 = -t169 * mrSges(4,1) + t170 * mrSges(4,2);
t101 = m(4) * t125 + mrSges(4,1) * t168 - t146 * mrSges(4,3) - t149 * t170 + t150 * t181 + t103;
t102 = m(4) * t126 - mrSges(4,2) * t168 + t145 * mrSges(4,3) + t149 * t169 - t151 * t181 + t203;
t99 = -t101 * t187 + t191 * t102;
t209 = mrSges(3,1) * t152 - mrSges(3,2) * t153 + Ifges(3,5) * t173 + Ifges(3,6) * t174 + Ifges(3,3) * qJDD(2) + pkin(2) * t107 + pkin(6) * t99 + t187 * t93 + t191 * t89 + (t160 * t188 - t192 * t161) * qJD(1);
t171 = (-mrSges(3,1) * t192 + mrSges(3,2) * t188) * qJD(1);
t177 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t207;
t106 = m(3) * t152 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t173 + qJD(2) * t177 - t171 * t208 + t107;
t176 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t208;
t97 = m(3) * t153 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t174 - qJD(2) * t176 + t171 * t207 + t99;
t204 = -t106 * t188 + t192 * t97;
t98 = t101 * t191 + t102 * t187;
t198 = -m(3) * t162 + t174 * mrSges(3,1) - mrSges(3,2) * t173 - t176 * t208 + t177 * t207 - t98;
t159 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t188 + Ifges(3,6) * t192) * qJD(1);
t86 = mrSges(3,2) * t162 - mrSges(3,3) * t152 + Ifges(3,1) * t173 + Ifges(3,4) * t174 + Ifges(3,5) * qJDD(2) - pkin(6) * t98 - qJD(2) * t160 + t159 * t207 - t187 * t89 + t191 * t93;
t199 = -mrSges(5,1) * t114 + mrSges(5,2) * t115 - Ifges(5,5) * t124 - Ifges(5,6) * t123 - Ifges(5,3) * t164 - t148 * t128 + t147 * t129;
t196 = mrSges(4,1) * t125 - mrSges(4,2) * t126 + Ifges(4,5) * t146 + Ifges(4,6) * t145 + Ifges(4,3) * t168 + pkin(3) * t103 + t170 * t141 - t169 * t142 - t199;
t88 = -mrSges(3,1) * t162 + mrSges(3,3) * t153 + Ifges(3,4) * t173 + Ifges(3,2) * t174 + Ifges(3,6) * qJDD(2) - pkin(2) * t98 + qJD(2) * t161 - t159 * t208 - t196;
t201 = mrSges(2,1) * t178 - mrSges(2,2) * t179 + Ifges(2,3) * qJDD(1) + pkin(1) * t198 + pkin(5) * t204 + t188 * t86 + t192 * t88;
t94 = m(2) * t178 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t195 + t198;
t92 = t106 * t192 + t188 * t97;
t90 = m(2) * t179 - mrSges(2,1) * t195 - qJDD(1) * mrSges(2,2) + t204;
t84 = mrSges(2,1) * g(3) + mrSges(2,3) * t179 + t195 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t92 - t209;
t83 = -mrSges(2,2) * g(3) - mrSges(2,3) * t178 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t195 - pkin(5) * t92 - t188 * t88 + t192 * t86;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t193 * t83 - t189 * t84 - pkin(4) * (t189 * t90 + t193 * t94), t83, t86, t93, t105; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t189 * t83 + t193 * t84 + pkin(4) * (-t189 * t94 + t193 * t90), t84, t88, t89, t104; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t201, t201, t209, t196, -t199;];
m_new = t1;
