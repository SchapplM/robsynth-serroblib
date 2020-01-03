% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR8
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
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR8_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:07
% EndTime: 2019-12-31 16:55:09
% DurationCPUTime: 1.02s
% Computational Cost: add. (9833->192), mult. (18988->240), div. (0->0), fcn. (10218->6), ass. (0->80)
t165 = sin(qJ(1));
t168 = cos(qJ(1));
t148 = -t168 * g(1) - t165 * g(2);
t180 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t148;
t191 = -pkin(1) - pkin(5);
t190 = mrSges(2,1) - mrSges(3,2);
t189 = Ifges(2,5) - Ifges(3,4);
t188 = (-Ifges(2,6) + Ifges(3,5));
t163 = sin(qJ(4));
t166 = cos(qJ(4));
t164 = sin(qJ(3));
t167 = cos(qJ(3));
t133 = (-t163 * t167 - t164 * t166) * qJD(1);
t185 = qJD(1) * qJD(3);
t141 = -t164 * qJDD(1) - t167 * t185;
t183 = t164 * t185;
t142 = t167 * qJDD(1) - t183;
t109 = t133 * qJD(4) + t163 * t141 + t166 * t142;
t134 = (-t163 * t164 + t166 * t167) * qJD(1);
t114 = -t133 * mrSges(5,1) + t134 * mrSges(5,2);
t154 = qJD(3) + qJD(4);
t121 = -t154 * mrSges(5,2) + t133 * mrSges(5,3);
t153 = qJDD(3) + qJDD(4);
t147 = t165 * g(1) - t168 * g(2);
t169 = qJD(1) ^ 2;
t179 = -t169 * qJ(2) + qJDD(2) - t147;
t124 = t191 * qJDD(1) + t179;
t116 = t164 * g(3) + t167 * t124;
t100 = (-t142 - t183) * pkin(6) + (-t164 * t167 * t169 + qJDD(3)) * pkin(3) + t116;
t117 = -t167 * g(3) + t164 * t124;
t186 = qJD(1) * t167;
t146 = qJD(3) * pkin(3) - pkin(6) * t186;
t160 = t164 ^ 2;
t101 = -t160 * t169 * pkin(3) + t141 * pkin(6) - qJD(3) * t146 + t117;
t98 = t166 * t100 - t163 * t101;
t95 = m(5) * t98 + t153 * mrSges(5,1) - t109 * mrSges(5,3) - t134 * t114 + t154 * t121;
t108 = -t134 * qJD(4) + t166 * t141 - t163 * t142;
t122 = t154 * mrSges(5,1) - t134 * mrSges(5,3);
t99 = t163 * t100 + t166 * t101;
t96 = m(5) * t99 - t153 * mrSges(5,2) + t108 * mrSges(5,3) + t133 * t114 - t154 * t122;
t86 = t163 * t96 + t166 * t95;
t187 = qJD(1) * t164;
t182 = -t163 * t95 + t166 * t96;
t140 = (mrSges(4,1) * t164 + mrSges(4,2) * t167) * qJD(1);
t144 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t187;
t83 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t142 * mrSges(4,3) + qJD(3) * t144 - t140 * t186 + t86;
t145 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t186;
t84 = m(4) * t117 - qJDD(3) * mrSges(4,2) + t141 * mrSges(4,3) - qJD(3) * t145 - t140 * t187 + t182;
t81 = -t164 * t83 + t167 * t84;
t80 = t164 * t84 + t167 * t83;
t104 = t146 * t186 - t141 * pkin(3) + (-pkin(6) * t160 + t191) * t169 + t180;
t178 = m(5) * t104 - t108 * mrSges(5,1) + t109 * mrSges(5,2) - t133 * t121 + t134 * t122;
t129 = -qJDD(1) * pkin(1) + t179;
t177 = -m(3) * t129 + (t169 * mrSges(3,3)) - t80;
t111 = Ifges(5,4) * t134 + Ifges(5,2) * t133 + Ifges(5,6) * t154;
t112 = Ifges(5,1) * t134 + Ifges(5,4) * t133 + Ifges(5,5) * t154;
t176 = mrSges(5,1) * t98 - mrSges(5,2) * t99 + Ifges(5,5) * t109 + Ifges(5,6) * t108 + Ifges(5,3) * t153 + t134 * t111 - t133 * t112;
t127 = t169 * pkin(1) - t180;
t123 = t191 * t169 + t180;
t130 = (Ifges(4,3) * qJD(3)) + (Ifges(4,5) * t167 - Ifges(4,6) * t164) * qJD(1);
t132 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t167 - Ifges(4,4) * t164) * qJD(1);
t110 = Ifges(5,5) * t134 + Ifges(5,6) * t133 + Ifges(5,3) * t154;
t87 = -mrSges(5,1) * t104 + mrSges(5,3) * t99 + Ifges(5,4) * t109 + Ifges(5,2) * t108 + Ifges(5,6) * t153 - t134 * t110 + t154 * t112;
t88 = mrSges(5,2) * t104 - mrSges(5,3) * t98 + Ifges(5,1) * t109 + Ifges(5,4) * t108 + Ifges(5,5) * t153 + t133 * t110 - t154 * t111;
t74 = -mrSges(4,1) * t123 + mrSges(4,3) * t117 + Ifges(4,4) * t142 + Ifges(4,2) * t141 + Ifges(4,6) * qJDD(3) - pkin(3) * t178 + pkin(6) * t182 + qJD(3) * t132 - t130 * t186 + t163 * t88 + t166 * t87;
t131 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t167 - Ifges(4,2) * t164) * qJD(1);
t76 = mrSges(4,2) * t123 - mrSges(4,3) * t116 + Ifges(4,1) * t142 + Ifges(4,4) * t141 + Ifges(4,5) * qJDD(3) - pkin(6) * t86 - qJD(3) * t131 - t130 * t187 - t163 * t87 + t166 * t88;
t175 = mrSges(3,2) * t129 - mrSges(3,3) * t127 + Ifges(3,1) * qJDD(1) - pkin(5) * t80 - t164 * t74 + t167 * t76;
t91 = -m(4) * t123 + t141 * mrSges(4,1) - t142 * mrSges(4,2) - t144 * t187 - t145 * t186 - t178;
t174 = -mrSges(3,1) * t127 - pkin(2) * t91 - pkin(5) * t81 - t164 * t76 - t167 * t74;
t172 = -m(3) * t127 + t169 * mrSges(3,2) + qJDD(1) * mrSges(3,3) - t91;
t173 = -mrSges(2,2) * t148 + mrSges(2,1) * t147 + Ifges(2,3) * qJDD(1) + t175 + pkin(1) * (-qJDD(1) * mrSges(3,2) + t177) + qJ(2) * t172;
t171 = mrSges(4,1) * t116 - mrSges(4,2) * t117 + Ifges(4,5) * t142 + Ifges(4,6) * t141 + Ifges(4,3) * qJDD(3) + pkin(3) * t86 + t131 * t186 + t132 * t187 + t176;
t170 = mrSges(3,1) * t129 + pkin(2) * t80 + t171;
t89 = m(2) * t148 - t169 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t172;
t79 = -m(3) * g(3) + t81;
t77 = m(2) * t147 - t169 * mrSges(2,2) + t190 * qJDD(1) + t177;
t73 = (t188 * t169) + t189 * qJDD(1) + (-mrSges(2,2) + mrSges(3,3)) * g(3) - mrSges(2,3) * t147 - qJ(2) * t79 + t170;
t72 = mrSges(2,3) * t148 - pkin(1) * t79 + t190 * g(3) - t188 * qJDD(1) + t189 * t169 + t174;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t168 * t73 - t165 * t72 - pkin(4) * (t165 * t89 + t168 * t77), t73, t175, t76, t88; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t165 * t73 + t168 * t72 + pkin(4) * (-t165 * t77 + t168 * t89), t72, -mrSges(3,3) * g(3) + Ifges(3,4) * qJDD(1) - (t169 * Ifges(3,5)) - t170, t74, t87; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t173, t173, mrSges(3,2) * g(3) + t169 * Ifges(3,4) + Ifges(3,5) * qJDD(1) - t174, t171, t176;];
m_new = t1;
