% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRP4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:28
% EndTime: 2019-12-31 17:34:29
% DurationCPUTime: 1.01s
% Computational Cost: add. (9104->196), mult. (16581->237), div. (0->0), fcn. (8150->6), ass. (0->79)
t165 = sin(pkin(7));
t166 = cos(pkin(7));
t150 = t165 * g(1) - t166 * g(2);
t148 = qJDD(2) - t150;
t151 = -t166 * g(1) - t165 * g(2);
t168 = sin(qJ(3));
t170 = cos(qJ(3));
t115 = t168 * t148 + t170 * t151;
t171 = qJD(3) ^ 2;
t112 = -t171 * pkin(3) + qJDD(3) * pkin(6) + t115;
t164 = g(3) - qJDD(1);
t169 = cos(qJ(4));
t158 = t169 * t164;
t167 = sin(qJ(4));
t108 = -t167 * t112 + t158;
t109 = t169 * t112 + t167 * t164;
t123 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t167 + Ifges(5,2) * t169) * qJD(3);
t124 = Ifges(6,5) * qJD(4) + (Ifges(6,1) * t167 + Ifges(6,4) * t169) * qJD(3);
t125 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t167 + Ifges(5,4) * t169) * qJD(3);
t189 = qJD(3) * qJD(4);
t185 = t169 * t189;
t143 = t167 * qJDD(3) + t185;
t144 = t169 * qJDD(3) - t167 * t189;
t188 = qJD(3) * qJD(5);
t196 = pkin(4) * t171;
t104 = qJDD(4) * pkin(4) + t158 + (-t143 + t185) * qJ(5) + (t169 * t196 - t112 - 0.2e1 * t188) * t167;
t191 = qJD(3) * t167;
t152 = qJD(4) * pkin(4) - qJ(5) * t191;
t163 = t169 ^ 2;
t105 = t144 * qJ(5) - qJD(4) * t152 - t163 * t196 + 0.2e1 * t169 * t188 + t109;
t122 = Ifges(6,6) * qJD(4) + (Ifges(6,4) * t167 + Ifges(6,2) * t169) * qJD(3);
t178 = -mrSges(6,1) * t104 + mrSges(6,2) * t105 - Ifges(6,5) * t143 - Ifges(6,6) * t144 - Ifges(6,3) * qJDD(4) - t122 * t191;
t141 = (-mrSges(6,1) * t169 + mrSges(6,2) * t167) * qJD(3);
t190 = qJD(3) * t169;
t155 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t190;
t187 = m(6) * t104 + qJDD(4) * mrSges(6,1) + qJD(4) * t155;
t99 = -t143 * mrSges(6,3) - t141 * t191 + t187;
t198 = mrSges(5,1) * t108 - mrSges(5,2) * t109 + Ifges(5,5) * t143 + Ifges(5,6) * t144 + Ifges(5,3) * qJDD(4) + pkin(4) * t99 - (-t167 * t123 + (t124 + t125) * t169) * qJD(3) - t178;
t195 = -mrSges(5,2) - mrSges(6,2);
t153 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t191;
t192 = -qJD(4) * mrSges(5,1) + mrSges(5,3) * t191 - t153;
t186 = m(6) * t105 + t144 * mrSges(6,3) + t141 * t190;
t142 = (-mrSges(5,1) * t169 + mrSges(5,2) * t167) * qJD(3);
t156 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t190;
t96 = m(5) * t108 + qJDD(4) * mrSges(5,1) + qJD(4) * t156 + (-mrSges(5,3) - mrSges(6,3)) * t143 + (-t141 - t142) * t191 + t187;
t97 = m(5) * t109 + t144 * mrSges(5,3) + t192 * qJD(4) + t195 * qJDD(4) + t142 * t190 + t186;
t93 = -t167 * t96 + t169 * t97;
t89 = m(4) * t115 - t171 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t93;
t114 = t170 * t148 - t168 * t151;
t182 = -qJDD(3) * pkin(3) - t114;
t111 = -t171 * pkin(6) + t182;
t107 = t152 * t191 - t144 * pkin(4) + qJDD(5) + (-qJ(5) * t163 - pkin(6)) * t171 + t182;
t184 = -m(6) * t107 + t144 * mrSges(6,1) + t155 * t190;
t98 = -m(5) * t111 + t144 * mrSges(5,1) + t195 * t143 + t156 * t190 + t192 * t191 + t184;
t94 = m(4) * t114 + qJDD(3) * mrSges(4,1) - t171 * mrSges(4,2) + t98;
t85 = -t168 * t94 + t170 * t89;
t183 = m(3) * t151 + t85;
t92 = t167 * t97 + t169 * t96;
t84 = t168 * t89 + t170 * t94;
t120 = Ifges(6,3) * qJD(4) + (Ifges(6,5) * t167 + Ifges(6,6) * t169) * qJD(3);
t121 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t167 + Ifges(5,6) * t169) * qJD(3);
t179 = -mrSges(6,1) * t107 + mrSges(6,3) * t105 + Ifges(6,4) * t143 + Ifges(6,2) * t144 + Ifges(6,6) * qJDD(4) + qJD(4) * t124;
t86 = Ifges(5,4) * t143 + Ifges(5,2) * t144 + Ifges(5,6) * qJDD(4) + qJD(4) * t125 - mrSges(5,1) * t111 + mrSges(5,3) * t109 - pkin(4) * (t143 * mrSges(6,2) - t184) + qJ(5) * (-qJDD(4) * mrSges(6,2) - qJD(4) * t153 + t186) + (-pkin(4) * t153 - t120 - t121) * t191 + t179;
t177 = mrSges(6,2) * t107 - mrSges(6,3) * t104 + Ifges(6,1) * t143 + Ifges(6,4) * t144 + Ifges(6,5) * qJDD(4) + t120 * t190;
t87 = t121 * t190 + mrSges(5,2) * t111 - mrSges(5,3) * t108 + Ifges(5,1) * t143 + Ifges(5,4) * t144 + Ifges(5,5) * qJDD(4) - qJ(5) * t99 + (-t122 - t123) * qJD(4) + t177;
t78 = mrSges(4,2) * t164 - mrSges(4,3) * t114 + Ifges(4,5) * qJDD(3) - t171 * Ifges(4,6) - pkin(6) * t92 - t167 * t86 + t169 * t87;
t79 = -mrSges(4,1) * t164 + mrSges(4,3) * t115 + t171 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t92 - t198;
t181 = mrSges(3,2) * t148 - pkin(5) * t84 - t168 * t79 + t170 * t78;
t180 = -m(3) * t148 - t84;
t176 = -pkin(2) * (-m(4) * t164 - t92) - pkin(5) * t85 - t168 * t78 - t170 * t79;
t175 = mrSges(4,1) * t114 - mrSges(4,2) * t115 + Ifges(4,3) * qJDD(3) + pkin(3) * t98 + pkin(6) * t93 + t167 * t87 + t169 * t86;
t174 = -mrSges(3,1) * t148 + mrSges(3,3) * t151 - pkin(2) * t84 - t175;
t172 = mrSges(2,1) * t150 - mrSges(2,2) * t151 + pkin(1) * t180 + qJ(2) * t183 + t174;
t90 = (-m(3) - m(4)) * t164 - t92;
t81 = m(2) * t151 + t183;
t80 = m(2) * t150 + t180;
t76 = -mrSges(2,3) * t150 - qJ(2) * t90 + (-mrSges(2,2) + mrSges(3,3)) * t164 + t181;
t75 = -pkin(1) * t90 + (mrSges(2,1) + mrSges(3,1)) * t164 + (mrSges(3,2) + mrSges(2,3)) * t151 + t176;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t166 * t76 - t165 * t75 - qJ(1) * (t165 * t81 + t166 * t80), t76, mrSges(3,3) * t164 + t181, t78, t87, -qJD(4) * t122 + t177; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t165 * t76 + t166 * t75 + qJ(1) * (-t165 * t80 + t166 * t81), t75, t174, t79, t86, -t120 * t191 + t179; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t172, t172, -mrSges(3,1) * t164 - mrSges(3,2) * t151 - t176, t175, t198, -t124 * t190 - t178;];
m_new = t1;
