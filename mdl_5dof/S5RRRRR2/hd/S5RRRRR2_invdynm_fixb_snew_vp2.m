% Calculate vector of cutting torques with Newton-Euler for
% S5RRRRR2
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:53:01
% EndTime: 2019-12-05 18:53:06
% DurationCPUTime: 2.09s
% Computational Cost: add. (29993->224), mult. (40720->291), div. (0->0), fcn. (27881->10), ass. (0->89)
t165 = qJD(1) + qJD(2);
t168 = sin(qJ(4));
t169 = sin(qJ(3));
t173 = cos(qJ(4));
t174 = cos(qJ(3));
t146 = (t168 * t169 - t173 * t174) * t165;
t171 = sin(qJ(1));
t176 = cos(qJ(1));
t157 = -t176 * g(1) - t171 * g(2);
t177 = qJD(1) ^ 2;
t152 = -t177 * pkin(1) + t157;
t170 = sin(qJ(2));
t175 = cos(qJ(2));
t156 = t171 * g(1) - t176 * g(2);
t185 = qJDD(1) * pkin(1) + t156;
t137 = t175 * t152 + t170 * t185;
t132 = -t174 * g(3) - t169 * t137;
t133 = -t169 * g(3) + t174 * t137;
t144 = Ifges(4,6) * qJD(3) + (Ifges(4,4) * t169 + Ifges(4,2) * t174) * t165;
t145 = Ifges(4,5) * qJD(3) + (Ifges(4,1) * t169 + Ifges(4,4) * t174) * t165;
t163 = qJDD(1) + qJDD(2);
t189 = qJD(3) * t165;
t149 = t169 * t163 + t174 * t189;
t188 = t169 * t189;
t150 = t174 * t163 - t188;
t161 = t165 ^ 2;
t127 = (-t161 * t174 ^ 2 - qJD(3) ^ 2) * pkin(2) + t133;
t182 = (t161 * t169 * t174 + qJDD(3)) * pkin(2) + t132;
t110 = t173 * t127 + t168 * t182;
t136 = t170 * t152 - t175 * t185;
t123 = (-t150 + t188) * pkin(2) + t136;
t167 = sin(qJ(5));
t172 = cos(qJ(5));
t107 = t172 * t110 + t167 * t123;
t109 = t168 * t127 - t173 * t182;
t122 = -t146 * qJD(4) + t173 * t149 + t168 * t150;
t147 = (t168 * t174 + t169 * t173) * t165;
t164 = qJD(3) + qJD(4);
t139 = t172 * t147 + t167 * t164;
t162 = qJDD(3) + qJDD(4);
t111 = -t139 * qJD(5) - t167 * t122 + t172 * t162;
t138 = -t167 * t147 + t172 * t164;
t112 = t138 * qJD(5) + t172 * t122 + t167 * t162;
t142 = qJD(5) + t146;
t113 = Ifges(6,5) * t139 + Ifges(6,6) * t138 + Ifges(6,3) * t142;
t115 = Ifges(6,1) * t139 + Ifges(6,4) * t138 + Ifges(6,5) * t142;
t121 = -t147 * qJD(4) - t168 * t149 + t173 * t150;
t120 = qJDD(5) - t121;
t100 = -mrSges(6,1) * t109 + mrSges(6,3) * t107 + Ifges(6,4) * t112 + Ifges(6,2) * t111 + Ifges(6,6) * t120 - t139 * t113 + t142 * t115;
t106 = -t167 * t110 + t172 * t123;
t114 = Ifges(6,4) * t139 + Ifges(6,2) * t138 + Ifges(6,6) * t142;
t101 = mrSges(6,2) * t109 - mrSges(6,3) * t106 + Ifges(6,1) * t112 + Ifges(6,4) * t111 + Ifges(6,5) * t120 + t138 * t113 - t142 * t114;
t129 = Ifges(5,4) * t147 - Ifges(5,2) * t146 + Ifges(5,6) * t164;
t130 = Ifges(5,1) * t147 - Ifges(5,4) * t146 + Ifges(5,5) * t164;
t181 = mrSges(5,1) * t109 + mrSges(5,2) * t110 - Ifges(5,5) * t122 - Ifges(5,6) * t121 - Ifges(5,3) * t162 - t172 * t100 - t167 * t101 - t147 * t129 - t146 * t130;
t125 = -t142 * mrSges(6,2) + t138 * mrSges(6,3);
t126 = t142 * mrSges(6,1) - t139 * mrSges(6,3);
t131 = t146 * mrSges(5,1) + t147 * mrSges(5,2);
t140 = -t164 * mrSges(5,2) - t146 * mrSges(5,3);
t103 = t162 * mrSges(5,1) + t111 * mrSges(6,1) - t112 * mrSges(6,2) - t122 * mrSges(5,3) + t138 * t125 - t139 * t126 - t147 * t131 + t164 * t140 + (-m(5) - m(6)) * t109;
t117 = -t138 * mrSges(6,1) + t139 * mrSges(6,2);
t104 = m(6) * t106 + t120 * mrSges(6,1) - t112 * mrSges(6,3) - t139 * t117 + t142 * t125;
t105 = m(6) * t107 - t120 * mrSges(6,2) + t111 * mrSges(6,3) + t138 * t117 - t142 * t126;
t141 = t164 * mrSges(5,1) - t147 * mrSges(5,3);
t96 = m(5) * t110 - t162 * mrSges(5,2) + t121 * mrSges(5,3) - t167 * t104 + t172 * t105 - t146 * t131 - t164 * t141;
t93 = t173 * t103 + t168 * t96;
t192 = mrSges(4,1) * t132 - mrSges(4,2) * t133 + Ifges(4,5) * t149 + Ifges(4,6) * t150 + Ifges(4,3) * qJDD(3) + pkin(2) * t93 + (t169 * t144 - t174 * t145) * t165 - t181;
t191 = t165 * t169;
t190 = t165 * t174;
t143 = Ifges(4,3) * qJD(3) + (Ifges(4,5) * t169 + Ifges(4,6) * t174) * t165;
t179 = m(5) * t123 - t121 * mrSges(5,1) + t122 * mrSges(5,2) + t172 * t104 + t167 * t105 + t146 * t140 + t147 * t141;
t128 = Ifges(5,5) * t147 - Ifges(5,6) * t146 + Ifges(5,3) * t164;
t94 = mrSges(5,2) * t123 + mrSges(5,3) * t109 + Ifges(5,1) * t122 + Ifges(5,4) * t121 + Ifges(5,5) * t162 - t167 * t100 + t172 * t101 - t146 * t128 - t164 * t129;
t180 = mrSges(6,1) * t106 - mrSges(6,2) * t107 + Ifges(6,5) * t112 + Ifges(6,6) * t111 + Ifges(6,3) * t120 + t139 * t114 - t138 * t115;
t97 = -mrSges(5,1) * t123 + mrSges(5,3) * t110 + Ifges(5,4) * t122 + Ifges(5,2) * t121 + Ifges(5,6) * t162 - t147 * t128 + t164 * t130 - t180;
t87 = -mrSges(4,1) * t136 + mrSges(4,3) * t133 + Ifges(4,4) * t149 + Ifges(4,2) * t150 + Ifges(4,6) * qJDD(3) - pkin(2) * t179 + qJD(3) * t145 - t143 * t191 + t168 * t94 + t173 * t97;
t90 = mrSges(4,2) * t136 - mrSges(4,3) * t132 + Ifges(4,1) * t149 + Ifges(4,4) * t150 + Ifges(4,5) * qJDD(3) - qJD(3) * t144 + t143 * t190 - t168 * t97 + t173 * t94;
t184 = -mrSges(3,1) * t136 - mrSges(3,2) * t137 + Ifges(3,3) * t163 + t169 * t90 + t174 * t87;
t153 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t191;
t154 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t190;
t148 = (-mrSges(4,1) * t174 + mrSges(4,2) * t169) * t165;
t91 = m(4) * t132 + qJDD(3) * mrSges(4,1) - t149 * mrSges(4,3) + qJD(3) * t154 - t148 * t191 + t93;
t92 = m(4) * t133 - qJDD(3) * mrSges(4,2) + t150 * mrSges(4,3) - qJD(3) * t153 - t168 * t103 + t148 * t190 + t173 * t96;
t183 = -mrSges(2,2) * t157 + mrSges(2,1) * t156 + Ifges(2,3) * qJDD(1) + t184 + pkin(1) * (t170 * (m(3) * t137 - t161 * mrSges(3,1) - t163 * mrSges(3,2) - t169 * t91 + t174 * t92) + t175 * (t163 * mrSges(3,1) + t150 * mrSges(4,1) - t161 * mrSges(3,2) - t149 * mrSges(4,2) + (-t153 * t169 + t154 * t174) * t165 + (-m(3) - m(4)) * t136 - t179));
t88 = mrSges(3,1) * g(3) + mrSges(3,3) * t137 + t161 * Ifges(3,5) + Ifges(3,6) * t163 - t192;
t84 = -mrSges(3,2) * g(3) + mrSges(3,3) * t136 + Ifges(3,5) * t163 - t161 * Ifges(3,6) - t169 * t87 + t174 * t90;
t83 = -mrSges(2,2) * g(3) - mrSges(2,3) * t156 + Ifges(2,5) * qJDD(1) - t177 * Ifges(2,6) - t170 * t88 + t175 * t84;
t82 = Ifges(2,6) * qJDD(1) + t177 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t157 + t170 * t84 + t175 * t88 - pkin(1) * (-m(3) * g(3) + t169 * t92 + t174 * t91);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - t171 * t82 + t176 * t83, t83, t84, t90, t94, t101; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t171 * t83 + t176 * t82, t82, t88, t87, t97, t100; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t183, t183, t184, t192, -t181, t180;];
m_new = t1;
