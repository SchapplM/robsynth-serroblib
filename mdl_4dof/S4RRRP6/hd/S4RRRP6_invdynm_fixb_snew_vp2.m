% Calculate vector of cutting torques with Newton-Euler for
% S4RRRP6
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
% Datum: 2019-12-31 17:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP6_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP6_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:18:12
% EndTime: 2019-12-31 17:18:15
% DurationCPUTime: 1.34s
% Computational Cost: add. (13387->230), mult. (26298->282), div. (0->0), fcn. (15108->6), ass. (0->86)
t185 = sin(qJ(3));
t188 = cos(qJ(3));
t186 = sin(qJ(2));
t209 = qJD(1) * t186;
t170 = t188 * qJD(2) - t185 * t209;
t189 = cos(qJ(2));
t207 = qJD(1) * qJD(2);
t203 = t189 * t207;
t174 = t186 * qJDD(1) + t203;
t143 = t170 * qJD(3) + t185 * qJDD(2) + t188 * t174;
t171 = t185 * qJD(2) + t188 * t209;
t145 = -t170 * mrSges(5,1) + t171 * mrSges(5,2);
t187 = sin(qJ(1));
t190 = cos(qJ(1));
t179 = t187 * g(1) - t190 * g(2);
t192 = qJD(1) ^ 2;
t160 = -qJDD(1) * pkin(1) - t192 * pkin(5) - t179;
t204 = t186 * t207;
t175 = t189 * qJDD(1) - t204;
t117 = (-t174 - t203) * pkin(6) + (-t175 + t204) * pkin(2) + t160;
t180 = -t190 * g(1) - t187 * g(2);
t161 = -t192 * pkin(1) + qJDD(1) * pkin(5) + t180;
t153 = -t186 * g(3) + t189 * t161;
t173 = (-pkin(2) * t189 - pkin(6) * t186) * qJD(1);
t191 = qJD(2) ^ 2;
t208 = t189 * qJD(1);
t122 = -t191 * pkin(2) + qJDD(2) * pkin(6) + t173 * t208 + t153;
t113 = t188 * t117 - t185 * t122;
t169 = qJDD(3) - t175;
t181 = qJD(3) - t208;
t107 = -0.2e1 * qJD(4) * t171 + (t170 * t181 - t143) * qJ(4) + (t170 * t171 + t169) * pkin(3) + t113;
t147 = -t181 * mrSges(5,2) + t170 * mrSges(5,3);
t206 = m(5) * t107 + t169 * mrSges(5,1) + t181 * t147;
t104 = -t143 * mrSges(5,3) - t171 * t145 + t206;
t114 = t185 * t117 + t188 * t122;
t128 = Ifges(4,4) * t171 + Ifges(4,2) * t170 + Ifges(4,6) * t181;
t129 = Ifges(5,1) * t171 + Ifges(5,4) * t170 + Ifges(5,5) * t181;
t130 = Ifges(4,1) * t171 + Ifges(4,4) * t170 + Ifges(4,5) * t181;
t142 = -t171 * qJD(3) + t188 * qJDD(2) - t185 * t174;
t149 = t181 * pkin(3) - t171 * qJ(4);
t168 = t170 ^ 2;
t110 = -t168 * pkin(3) + t142 * qJ(4) + 0.2e1 * qJD(4) * t170 - t181 * t149 + t114;
t127 = Ifges(5,4) * t171 + Ifges(5,2) * t170 + Ifges(5,6) * t181;
t198 = -mrSges(5,1) * t107 + mrSges(5,2) * t110 - Ifges(5,5) * t143 - Ifges(5,6) * t142 - Ifges(5,3) * t169 - t171 * t127;
t214 = mrSges(4,1) * t113 - mrSges(4,2) * t114 + Ifges(4,5) * t143 + Ifges(4,6) * t142 + Ifges(4,3) * t169 + pkin(3) * t104 + t171 * t128 - (t130 + t129) * t170 - t198;
t152 = -t189 * g(3) - t186 * t161;
t121 = -qJDD(2) * pkin(2) - t191 * pkin(6) + t173 * t209 - t152;
t148 = -t181 * mrSges(4,2) + t170 * mrSges(4,3);
t112 = -t142 * pkin(3) - t168 * qJ(4) + t171 * t149 + qJDD(4) + t121;
t201 = -m(5) * t112 + t142 * mrSges(5,1) + t170 * t147;
t150 = t181 * mrSges(5,1) - t171 * mrSges(5,3);
t210 = -t181 * mrSges(4,1) + t171 * mrSges(4,3) - t150;
t212 = -mrSges(4,2) - mrSges(5,2);
t103 = -m(4) * t121 + t142 * mrSges(4,1) + t212 * t143 + t170 * t148 + t210 * t171 + t201;
t158 = Ifges(3,6) * qJD(2) + (Ifges(3,4) * t186 + Ifges(3,2) * t189) * qJD(1);
t159 = Ifges(3,5) * qJD(2) + (Ifges(3,1) * t186 + Ifges(3,4) * t189) * qJD(1);
t125 = Ifges(5,5) * t171 + Ifges(5,6) * t170 + Ifges(5,3) * t181;
t126 = Ifges(4,5) * t171 + Ifges(4,6) * t170 + Ifges(4,3) * t181;
t199 = -mrSges(5,1) * t112 + mrSges(5,3) * t110 + Ifges(5,4) * t143 + Ifges(5,2) * t142 + Ifges(5,6) * t169 + t181 * t129;
t205 = m(5) * t110 + t142 * mrSges(5,3) + t170 * t145;
t91 = Ifges(4,4) * t143 + Ifges(4,2) * t142 + Ifges(4,6) * t169 + t181 * t130 - mrSges(4,1) * t121 + mrSges(4,3) * t114 - pkin(3) * (t143 * mrSges(5,2) - t201) + qJ(4) * (-t169 * mrSges(5,2) - t181 * t150 + t205) + (-pkin(3) * t150 - t125 - t126) * t171 + t199;
t197 = mrSges(5,2) * t112 - mrSges(5,3) * t107 + Ifges(5,1) * t143 + Ifges(5,4) * t142 + Ifges(5,5) * t169 + t170 * t125;
t93 = mrSges(4,2) * t121 - mrSges(4,3) * t113 + Ifges(4,1) * t143 + Ifges(4,4) * t142 + Ifges(4,5) * t169 - qJ(4) * t104 + t170 * t126 + (-t127 - t128) * t181 + t197;
t146 = -t170 * mrSges(4,1) + t171 * mrSges(4,2);
t101 = m(4) * t114 + t142 * mrSges(4,3) + t170 * t146 + t212 * t169 + t210 * t181 + t205;
t99 = m(4) * t113 + t169 * mrSges(4,1) + t181 * t148 + (-t145 - t146) * t171 + (-mrSges(4,3) - mrSges(5,3)) * t143 + t206;
t98 = t188 * t101 - t185 * t99;
t213 = mrSges(3,1) * t152 - mrSges(3,2) * t153 + Ifges(3,5) * t174 + Ifges(3,6) * t175 + Ifges(3,3) * qJDD(2) + pkin(2) * t103 + pkin(6) * t98 + t185 * t93 + t188 * t91 + (t186 * t158 - t189 * t159) * qJD(1);
t172 = (-mrSges(3,1) * t189 + mrSges(3,2) * t186) * qJD(1);
t178 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t208;
t102 = m(3) * t152 + qJDD(2) * mrSges(3,1) - t174 * mrSges(3,3) + qJD(2) * t178 - t172 * t209 + t103;
t177 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t209;
t96 = m(3) * t153 - qJDD(2) * mrSges(3,2) + t175 * mrSges(3,3) - qJD(2) * t177 + t172 * t208 + t98;
t202 = -t186 * t102 + t189 * t96;
t97 = t185 * t101 + t188 * t99;
t195 = -m(3) * t160 + t175 * mrSges(3,1) - t174 * mrSges(3,2) - t177 * t209 + t178 * t208 - t97;
t157 = Ifges(3,3) * qJD(2) + (Ifges(3,5) * t186 + Ifges(3,6) * t189) * qJD(1);
t85 = mrSges(3,2) * t160 - mrSges(3,3) * t152 + Ifges(3,1) * t174 + Ifges(3,4) * t175 + Ifges(3,5) * qJDD(2) - pkin(6) * t97 - qJD(2) * t158 + t157 * t208 - t185 * t91 + t188 * t93;
t87 = -mrSges(3,1) * t160 + mrSges(3,3) * t153 + Ifges(3,4) * t174 + Ifges(3,2) * t175 + Ifges(3,6) * qJDD(2) - pkin(2) * t97 + qJD(2) * t159 - t157 * t209 - t214;
t196 = mrSges(2,1) * t179 - mrSges(2,2) * t180 + Ifges(2,3) * qJDD(1) + pkin(1) * t195 + pkin(5) * t202 + t186 * t85 + t189 * t87;
t92 = m(2) * t179 + qJDD(1) * mrSges(2,1) - t192 * mrSges(2,2) + t195;
t90 = t189 * t102 + t186 * t96;
t88 = m(2) * t180 - t192 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t202;
t83 = mrSges(2,1) * g(3) + mrSges(2,3) * t180 + t192 * Ifges(2,5) + Ifges(2,6) * qJDD(1) - pkin(1) * t90 - t213;
t82 = -mrSges(2,2) * g(3) - mrSges(2,3) * t179 + Ifges(2,5) * qJDD(1) - t192 * Ifges(2,6) - pkin(5) * t90 - t186 * t87 + t189 * t85;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t190 * t82 - t187 * t83 - pkin(4) * (t187 * t88 + t190 * t92), t82, t85, t93, -t181 * t127 + t197; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t187 * t82 + t190 * t83 + pkin(4) * (-t187 * t92 + t190 * t88), t83, t87, t91, -t171 * t125 + t199; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t196, t196, t213, t214, -t170 * t129 - t198;];
m_new = t1;
