% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR9
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR9_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR9_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:24
% EndTime: 2019-12-31 17:09:30
% DurationCPUTime: 2.60s
% Computational Cost: add. (29340->235), mult. (63010->302), div. (0->0), fcn. (39449->8), ass. (0->94)
t190 = sin(qJ(1));
t193 = cos(qJ(1));
t180 = -g(1) * t193 - g(2) * t190;
t195 = qJD(1) ^ 2;
t164 = -pkin(1) * t195 + qJDD(1) * pkin(5) + t180;
t189 = sin(qJ(2));
t192 = cos(qJ(2));
t147 = -t192 * g(3) - t189 * t164;
t172 = (-pkin(2) * t192 - qJ(3) * t189) * qJD(1);
t194 = qJD(2) ^ 2;
t208 = qJD(1) * t189;
t136 = -qJDD(2) * pkin(2) - qJ(3) * t194 + t172 * t208 + qJDD(3) - t147;
t186 = sin(pkin(7));
t187 = cos(pkin(7));
t168 = qJD(2) * t187 - t186 * t208;
t207 = t192 * qJD(1);
t150 = mrSges(4,2) * t207 + mrSges(4,3) * t168;
t169 = qJD(2) * t186 + t187 * t208;
t151 = -mrSges(4,1) * t207 - mrSges(4,3) * t169;
t206 = qJD(1) * qJD(2);
t205 = t192 * t206;
t174 = qJDD(1) * t189 + t205;
t152 = qJDD(2) * t187 - t174 * t186;
t153 = qJDD(2) * t186 + t174 * t187;
t154 = -pkin(3) * t207 - pkin(6) * t169;
t167 = t168 ^ 2;
t120 = -pkin(3) * t152 - pkin(6) * t167 + t154 * t169 + t136;
t188 = sin(qJ(4));
t191 = cos(qJ(4));
t145 = t168 * t188 + t169 * t191;
t125 = -qJD(4) * t145 + t152 * t191 - t153 * t188;
t144 = t168 * t191 - t169 * t188;
t126 = qJD(4) * t144 + t152 * t188 + t153 * t191;
t181 = qJD(4) - t207;
t138 = -mrSges(5,2) * t181 + mrSges(5,3) * t144;
t139 = mrSges(5,1) * t181 - mrSges(5,3) * t145;
t200 = m(5) * t120 - t125 * mrSges(5,1) + mrSges(5,2) * t126 - t144 * t138 + t139 * t145;
t112 = -m(4) * t136 + t152 * mrSges(4,1) - mrSges(4,2) * t153 + t168 * t150 - t151 * t169 - t200;
t148 = -g(3) * t189 + t192 * t164;
t161 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t189 + Ifges(3,2) * t192) * qJD(1);
t162 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t189 + Ifges(3,4) * t192) * qJD(1);
t182 = t189 * t206;
t175 = qJDD(1) * t192 - t182;
t179 = t190 * g(1) - t193 * g(2);
t163 = -qJDD(1) * pkin(1) - t195 * pkin(5) - t179;
t134 = (-t174 - t205) * qJ(3) + (-t175 + t182) * pkin(2) + t163;
t137 = -pkin(2) * t194 + qJDD(2) * qJ(3) + t172 * t207 + t148;
t118 = -0.2e1 * qJD(3) * t169 + t187 * t134 - t186 * t137;
t116 = (-t168 * t207 - t153) * pkin(6) + (t168 * t169 - t175) * pkin(3) + t118;
t119 = 0.2e1 * qJD(3) * t168 + t186 * t134 + t187 * t137;
t117 = -pkin(3) * t167 + pkin(6) * t152 + t154 * t207 + t119;
t115 = t116 * t188 + t117 * t191;
t127 = Ifges(5,5) * t145 + Ifges(5,6) * t144 + Ifges(5,3) * t181;
t129 = Ifges(5,1) * t145 + Ifges(5,4) * t144 + Ifges(5,5) * t181;
t171 = qJDD(4) - t175;
t104 = -mrSges(5,1) * t120 + mrSges(5,3) * t115 + Ifges(5,4) * t126 + Ifges(5,2) * t125 + Ifges(5,6) * t171 - t127 * t145 + t129 * t181;
t114 = t116 * t191 - t117 * t188;
t128 = Ifges(5,4) * t145 + Ifges(5,2) * t144 + Ifges(5,6) * t181;
t105 = mrSges(5,2) * t120 - mrSges(5,3) * t114 + Ifges(5,1) * t126 + Ifges(5,4) * t125 + Ifges(5,5) * t171 + t127 * t144 - t128 * t181;
t140 = Ifges(4,5) * t169 + Ifges(4,6) * t168 - Ifges(4,3) * t207;
t142 = Ifges(4,1) * t169 + Ifges(4,4) * t168 - Ifges(4,5) * t207;
t131 = -mrSges(5,1) * t144 + mrSges(5,2) * t145;
t110 = m(5) * t114 + mrSges(5,1) * t171 - mrSges(5,3) * t126 - t131 * t145 + t138 * t181;
t111 = m(5) * t115 - mrSges(5,2) * t171 + mrSges(5,3) * t125 + t131 * t144 - t139 * t181;
t203 = -t110 * t188 + t111 * t191;
t89 = -mrSges(4,1) * t136 + mrSges(4,3) * t119 + Ifges(4,4) * t153 + Ifges(4,2) * t152 - Ifges(4,6) * t175 - pkin(3) * t200 + pkin(6) * t203 + t191 * t104 + t188 * t105 - t169 * t140 - t142 * t207;
t103 = t110 * t191 + t111 * t188;
t141 = Ifges(4,4) * t169 + Ifges(4,2) * t168 - Ifges(4,6) * t207;
t93 = mrSges(4,2) * t136 - mrSges(4,3) * t118 + Ifges(4,1) * t153 + Ifges(4,4) * t152 - Ifges(4,5) * t175 - pkin(6) * t103 - t104 * t188 + t105 * t191 + t140 * t168 + t141 * t207;
t146 = -mrSges(4,1) * t168 + mrSges(4,2) * t169;
t101 = m(4) * t118 - mrSges(4,1) * t175 - mrSges(4,3) * t153 - t146 * t169 - t150 * t207 + t103;
t102 = m(4) * t119 + mrSges(4,2) * t175 + mrSges(4,3) * t152 + t146 * t168 + t151 * t207 + t203;
t99 = -t101 * t186 + t102 * t187;
t209 = mrSges(3,1) * t147 - mrSges(3,2) * t148 + Ifges(3,5) * t174 + Ifges(3,6) * t175 + Ifges(3,3) * qJDD(2) + pkin(2) * t112 + qJ(3) * t99 + t186 * t93 + t187 * t89 + (t161 * t189 - t162 * t192) * qJD(1);
t173 = (-mrSges(3,1) * t192 + mrSges(3,2) * t189) * qJD(1);
t178 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t207;
t106 = m(3) * t147 + qJDD(2) * mrSges(3,1) - mrSges(3,3) * t174 + qJD(2) * t178 - t173 * t208 + t112;
t177 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t208;
t97 = m(3) * t148 - qJDD(2) * mrSges(3,2) + mrSges(3,3) * t175 - qJD(2) * t177 + t173 * t207 + t99;
t204 = -t106 * t189 + t192 * t97;
t98 = t101 * t187 + t102 * t186;
t198 = -m(3) * t163 + t175 * mrSges(3,1) - mrSges(3,2) * t174 - t177 * t208 + t178 * t207 - t98;
t160 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t189 + Ifges(3,6) * t192) * qJD(1);
t86 = mrSges(3,2) * t163 - mrSges(3,3) * t147 + Ifges(3,1) * t174 + Ifges(3,4) * t175 + Ifges(3,5) * qJDD(2) - qJ(3) * t98 - qJD(2) * t161 + t160 * t207 - t186 * t89 + t187 * t93;
t199 = -mrSges(5,1) * t114 + mrSges(5,2) * t115 - Ifges(5,5) * t126 - Ifges(5,6) * t125 - Ifges(5,3) * t171 - t145 * t128 + t144 * t129;
t196 = -mrSges(4,1) * t118 + mrSges(4,2) * t119 - Ifges(4,5) * t153 - Ifges(4,6) * t152 - pkin(3) * t103 - t169 * t141 + t168 * t142 + t199;
t88 = t196 + (Ifges(3,2) + Ifges(4,3)) * t175 + Ifges(3,4) * t174 + qJD(2) * t162 - mrSges(3,1) * t163 + mrSges(3,3) * t148 + Ifges(3,6) * qJDD(2) - pkin(2) * t98 - t160 * t208;
t201 = mrSges(2,1) * t179 - mrSges(2,2) * t180 + Ifges(2,3) * qJDD(1) + pkin(1) * t198 + pkin(5) * t204 + t189 * t86 + t192 * t88;
t94 = m(2) * t179 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t195 + t198;
t92 = t106 * t192 + t189 * t97;
t90 = m(2) * t180 - mrSges(2,1) * t195 - qJDD(1) * mrSges(2,2) + t204;
t84 = mrSges(2,1) * g(3) + mrSges(2,3) * t180 + t195 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t92 - t209;
t83 = -mrSges(2,2) * g(3) - mrSges(2,3) * t179 + Ifges(2,5) * qJDD(1) - Ifges(2,6) * t195 - pkin(5) * t92 - t189 * t88 + t192 * t86;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t193 * t83 - t190 * t84 - pkin(4) * (t190 * t90 + t193 * t94), t83, t86, t93, t105; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t190 * t83 + t193 * t84 + pkin(4) * (-t190 * t94 + t193 * t90), t84, t88, t89, t104; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t201, t201, t209, -Ifges(4,3) * t175 - t196, -t199;];
m_new = t1;
