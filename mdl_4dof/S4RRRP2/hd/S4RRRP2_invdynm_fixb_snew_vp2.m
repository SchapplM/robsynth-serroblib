% Calculate vector of cutting torques with Newton-Euler for
% S4RRRP2
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
% Datum: 2019-12-31 17:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRRP2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:12:55
% EndTime: 2019-12-31 17:12:56
% DurationCPUTime: 0.92s
% Computational Cost: add. (10007->190), mult. (12764->236), div. (0->0), fcn. (5708->6), ass. (0->76)
t167 = sin(qJ(1));
t170 = cos(qJ(1));
t152 = t167 * g(1) - t170 * g(2);
t144 = qJDD(1) * pkin(1) + t152;
t153 = -t170 * g(1) - t167 * g(2);
t171 = qJD(1) ^ 2;
t145 = -t171 * pkin(1) + t153;
t166 = sin(qJ(2));
t169 = cos(qJ(2));
t114 = t166 * t144 + t169 * t145;
t158 = qJD(1) + qJD(2);
t156 = t158 ^ 2;
t157 = qJDD(1) + qJDD(2);
t111 = -t156 * pkin(2) + t157 * pkin(6) + t114;
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t194 = t168 * g(3);
t107 = -t165 * t111 - t194;
t108 = -t165 * g(3) + t168 * t111;
t122 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t165 + Ifges(4,2) * t168) * t158;
t123 = Ifges(5,5) * qJD(3) + (Ifges(5,1) * t165 + Ifges(5,4) * t168) * t158;
t124 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t165 + Ifges(4,4) * t168) * t158;
t187 = qJD(3) * t158;
t183 = t168 * t187;
t138 = t165 * t157 + t183;
t139 = t168 * t157 - t165 * t187;
t186 = qJD(4) * t158;
t195 = pkin(3) * t156;
t103 = qJDD(3) * pkin(3) - t194 + (-t138 + t183) * qJ(4) + (t168 * t195 - t111 - 0.2e1 * t186) * t165;
t192 = t158 * t165;
t146 = qJD(3) * pkin(3) - qJ(4) * t192;
t164 = t168 ^ 2;
t104 = t139 * qJ(4) - qJD(3) * t146 - t164 * t195 + 0.2e1 * t168 * t186 + t108;
t121 = Ifges(5,6) * qJD(3) + (Ifges(5,4) * t165 + Ifges(5,2) * t168) * t158;
t178 = -mrSges(5,1) * t103 + mrSges(5,2) * t104 - Ifges(5,5) * t138 - Ifges(5,6) * t139 - Ifges(5,3) * qJDD(3) - t121 * t192;
t136 = (-mrSges(5,1) * t168 + mrSges(5,2) * t165) * t158;
t191 = t158 * t168;
t149 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t191;
t185 = m(5) * t103 + qJDD(3) * mrSges(5,1) + qJD(3) * t149;
t98 = -t138 * mrSges(5,3) - t136 * t192 + t185;
t197 = mrSges(4,1) * t107 - mrSges(4,2) * t108 + Ifges(4,5) * t138 + Ifges(4,6) * t139 + Ifges(4,3) * qJDD(3) + pkin(3) * t98 - (-t165 * t122 + (t123 + t124) * t168) * t158 - t178;
t193 = -mrSges(4,2) - mrSges(5,2);
t137 = (-mrSges(4,1) * t168 + mrSges(4,2) * t165) * t158;
t150 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t191;
t96 = m(4) * t107 + qJDD(3) * mrSges(4,1) + qJD(3) * t150 + (-t136 - t137) * t192 + (-mrSges(4,3) - mrSges(5,3)) * t138 + t185;
t184 = m(5) * t104 + t139 * mrSges(5,3) + t136 * t191;
t147 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t192;
t188 = -qJD(3) * mrSges(4,1) + mrSges(4,3) * t192 - t147;
t97 = m(4) * t108 + t139 * mrSges(4,3) + t188 * qJD(3) + t193 * qJDD(3) + t137 * t191 + t184;
t182 = -t165 * t96 + t168 * t97;
t87 = m(3) * t114 - t156 * mrSges(3,1) - t157 * mrSges(3,2) + t182;
t113 = t169 * t144 - t166 * t145;
t179 = -t157 * pkin(2) - t113;
t110 = -t156 * pkin(6) + t179;
t106 = t146 * t192 - t139 * pkin(3) + qJDD(4) + (-qJ(4) * t164 - pkin(6)) * t156 + t179;
t180 = -m(5) * t106 + t139 * mrSges(5,1) + t149 * t191;
t173 = -m(4) * t110 + t139 * mrSges(4,1) + t193 * t138 + t150 * t191 + t188 * t192 + t180;
t91 = m(3) * t113 + t157 * mrSges(3,1) - t156 * mrSges(3,2) + t173;
t80 = t166 * t87 + t169 * t91;
t89 = t165 * t97 + t168 * t96;
t181 = -t166 * t91 + t169 * t87;
t177 = -mrSges(5,1) * t106 + mrSges(5,3) * t104 + Ifges(5,4) * t138 + Ifges(5,2) * t139 + Ifges(5,6) * qJDD(3) + qJD(3) * t123;
t119 = Ifges(5,3) * qJD(3) + (Ifges(5,5) * t165 + Ifges(5,6) * t168) * t158;
t176 = mrSges(5,2) * t106 - mrSges(5,3) * t103 + Ifges(5,1) * t138 + Ifges(5,4) * t139 + Ifges(5,5) * qJDD(3) + t119 * t191;
t120 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t165 + Ifges(4,6) * t168) * t158;
t82 = Ifges(4,4) * t138 + Ifges(4,2) * t139 + Ifges(4,6) * qJDD(3) + qJD(3) * t124 - mrSges(4,1) * t110 + mrSges(4,3) * t108 - pkin(3) * (t138 * mrSges(5,2) - t180) + qJ(4) * (-qJDD(3) * mrSges(5,2) - qJD(3) * t147 + t184) + (-pkin(3) * t147 - t119 - t120) * t192 + t177;
t84 = t120 * t191 + mrSges(4,2) * t110 - mrSges(4,3) * t107 + Ifges(4,1) * t138 + Ifges(4,4) * t139 + Ifges(4,5) * qJDD(3) - qJ(4) * t98 + (-t121 - t122) * qJD(3) + t176;
t175 = mrSges(3,1) * t113 - mrSges(3,2) * t114 + Ifges(3,3) * t157 + pkin(2) * t173 + pkin(6) * t182 + t165 * t84 + t168 * t82;
t174 = mrSges(2,1) * t152 - mrSges(2,2) * t153 + Ifges(2,3) * qJDD(1) + pkin(1) * t80 + t175;
t78 = m(2) * t153 - t171 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t181;
t77 = m(2) * t152 + qJDD(1) * mrSges(2,1) - t171 * mrSges(2,2) + t80;
t76 = mrSges(3,1) * g(3) + mrSges(3,3) * t114 + t156 * Ifges(3,5) + Ifges(3,6) * t157 - pkin(2) * t89 - t197;
t75 = -mrSges(3,2) * g(3) - mrSges(3,3) * t113 + Ifges(3,5) * t157 - t156 * Ifges(3,6) - pkin(6) * t89 - t165 * t82 + t168 * t84;
t74 = -mrSges(2,2) * g(3) - mrSges(2,3) * t152 + Ifges(2,5) * qJDD(1) - t171 * Ifges(2,6) - pkin(5) * t80 - t166 * t76 + t169 * t75;
t73 = Ifges(2,6) * qJDD(1) + t171 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t153 + t166 * t75 + t169 * t76 - pkin(1) * (-m(3) * g(3) + t89) + pkin(5) * t181;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t170 * t74 - t167 * t73 - pkin(4) * (t167 * t78 + t170 * t77), t74, t75, t84, -qJD(3) * t121 + t176; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t167 * t74 + t170 * t73 + pkin(4) * (-t167 * t77 + t170 * t78), t73, t76, t82, -t119 * t192 + t177; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t174, t174, t175, t197, -t123 * t191 - t178;];
m_new = t1;
