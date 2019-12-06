% Calculate vector of cutting torques with Newton-Euler for
% S5PPRRP3
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRRP3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:10:36
% EndTime: 2019-12-05 15:10:39
% DurationCPUTime: 1.47s
% Computational Cost: add. (14931->200), mult. (25663->252), div. (0->0), fcn. (14425->8), ass. (0->80)
t166 = sin(pkin(7));
t168 = cos(pkin(7));
t154 = -t168 * g(1) - t166 * g(2);
t165 = sin(pkin(8));
t167 = cos(pkin(8));
t191 = g(3) - qJDD(1);
t124 = t167 * t154 - t165 * t191;
t153 = t166 * g(1) - t168 * g(2);
t152 = qJDD(2) - t153;
t170 = sin(qJ(3));
t172 = cos(qJ(3));
t119 = t172 * t124 + t170 * t152;
t174 = qJD(3) ^ 2;
t117 = -t174 * pkin(3) + qJDD(3) * pkin(6) + t119;
t169 = sin(qJ(4));
t123 = t165 * t154 + t167 * t191;
t171 = cos(qJ(4));
t192 = t171 * t123;
t113 = -t169 * t117 + t192;
t114 = t171 * t117 + t169 * t123;
t128 = Ifges(6,6) * qJD(4) + (Ifges(6,5) * t169 - Ifges(6,3) * t171) * qJD(3);
t131 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t169 + Ifges(5,2) * t171) * qJD(3);
t145 = (-mrSges(6,1) * t171 - mrSges(6,3) * t169) * qJD(3);
t186 = qJD(3) * qJD(4);
t147 = t169 * qJDD(3) + t171 * t186;
t148 = t171 * qJDD(3) - t169 * t186;
t144 = (-pkin(4) * t171 - qJ(5) * t169) * qJD(3);
t173 = qJD(4) ^ 2;
t187 = qJD(3) * t171;
t109 = -t173 * pkin(4) + qJDD(4) * qJ(5) + 0.2e1 * qJD(5) * qJD(4) + t144 * t187 + t114;
t111 = -qJDD(4) * pkin(4) - t173 * qJ(5) - t192 + qJDD(5) + (qJD(3) * t144 + t117) * t169;
t181 = -mrSges(6,1) * t111 + mrSges(6,3) * t109 + Ifges(6,4) * t147 + Ifges(6,2) * qJDD(4) - Ifges(6,6) * t148;
t158 = mrSges(6,2) * t187 + qJD(4) * mrSges(6,3);
t183 = -m(6) * t111 + qJDD(4) * mrSges(6,1) + qJD(4) * t158;
t188 = qJD(3) * t169;
t156 = -qJD(4) * mrSges(6,1) + mrSges(6,2) * t188;
t184 = m(6) * t109 + qJDD(4) * mrSges(6,3) + qJD(4) * t156 + t145 * t187;
t132 = Ifges(6,4) * qJD(4) + (Ifges(6,1) * t169 - Ifges(6,5) * t171) * qJD(3);
t189 = t132 + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t169 + Ifges(5,4) * t171) * qJD(3);
t195 = -(t189 * t171 + (t128 - t131) * t169) * qJD(3) + mrSges(5,1) * t113 - mrSges(5,2) * t114 + Ifges(5,5) * t147 + Ifges(5,6) * t148 + Ifges(5,3) * qJDD(4) + pkin(4) * (-t147 * mrSges(6,2) - t145 * t188 + t183) + qJ(5) * (t148 * mrSges(6,2) + t184) + t181;
t193 = mrSges(5,3) + mrSges(6,2);
t146 = (-mrSges(5,1) * t171 + mrSges(5,2) * t169) * qJD(3);
t155 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t188;
t101 = m(5) * t114 - qJDD(4) * mrSges(5,2) - qJD(4) * t155 + t146 * t187 + t193 * t148 + t184;
t157 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t187;
t102 = m(5) * t113 + qJDD(4) * mrSges(5,1) + qJD(4) * t157 - t193 * t147 + (-t145 - t146) * t188 + t183;
t98 = t171 * t101 - t169 * t102;
t94 = m(4) * t119 - t174 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t98;
t118 = -t170 * t124 + t172 * t152;
t116 = -qJDD(3) * pkin(3) - t174 * pkin(6) - t118;
t112 = -t148 * pkin(4) - t147 * qJ(5) + (-0.2e1 * qJD(5) * t169 + (pkin(4) * t169 - qJ(5) * t171) * qJD(4)) * qJD(3) + t116;
t104 = m(6) * t112 - t148 * mrSges(6,1) - t147 * mrSges(6,3) - t156 * t188 - t158 * t187;
t103 = -m(5) * t116 + t148 * mrSges(5,1) - t147 * mrSges(5,2) - t155 * t188 + t157 * t187 - t104;
t99 = m(4) * t118 + qJDD(3) * mrSges(4,1) - t174 * mrSges(4,2) + t103;
t90 = -t170 * t99 + t172 * t94;
t87 = m(3) * t124 + t90;
t97 = t169 * t101 + t171 * t102;
t95 = (-m(3) - m(4)) * t123 - t97;
t185 = -t165 * t95 + t167 * t87;
t182 = -mrSges(6,1) * t112 + mrSges(6,2) * t109;
t89 = t170 * t94 + t172 * t99;
t179 = -m(3) * t152 - t89;
t129 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t169 + Ifges(5,6) * t171) * qJD(3);
t130 = Ifges(6,2) * qJD(4) + (Ifges(6,4) * t169 - Ifges(6,6) * t171) * qJD(3);
t91 = -mrSges(5,1) * t116 + mrSges(5,3) * t114 - pkin(4) * t104 + (Ifges(5,2) + Ifges(6,3)) * t148 + (Ifges(5,4) - Ifges(6,5)) * t147 + (Ifges(5,6) - Ifges(6,6)) * qJDD(4) + t189 * qJD(4) + (-t129 - t130) * t188 + t182;
t178 = mrSges(6,2) * t111 - mrSges(6,3) * t112 + Ifges(6,1) * t147 + Ifges(6,4) * qJDD(4) - Ifges(6,5) * t148 + qJD(4) * t128 + t130 * t187;
t92 = mrSges(5,2) * t116 - mrSges(5,3) * t113 + Ifges(5,1) * t147 + Ifges(5,4) * t148 + Ifges(5,5) * qJDD(4) - qJ(5) * t104 - qJD(4) * t131 + t129 * t187 + t178;
t80 = mrSges(4,2) * t123 - mrSges(4,3) * t118 + Ifges(4,5) * qJDD(3) - t174 * Ifges(4,6) - pkin(6) * t97 - t169 * t91 + t171 * t92;
t84 = -mrSges(4,1) * t123 + mrSges(4,3) * t119 + t174 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t97 - t195;
t77 = mrSges(3,2) * t152 + mrSges(3,3) * t123 - pkin(5) * t89 - t170 * t84 + t172 * t80;
t176 = mrSges(4,1) * t118 - mrSges(4,2) * t119 + Ifges(4,3) * qJDD(3) + pkin(3) * t103 + pkin(6) * t98 + t169 * t92 + t171 * t91;
t79 = -mrSges(3,1) * t152 + mrSges(3,3) * t124 - pkin(2) * t89 - t176;
t180 = mrSges(2,1) * t153 - mrSges(2,2) * t154 + pkin(1) * t179 + qJ(2) * t185 + t165 * t77 + t167 * t79;
t177 = mrSges(3,1) * t123 + mrSges(3,2) * t124 - pkin(2) * (-m(4) * t123 - t97) - pkin(5) * t90 - t170 * t80 - t172 * t84;
t86 = m(2) * t153 + t179;
t83 = t165 * t87 + t167 * t95;
t81 = m(2) * t154 + t185;
t75 = mrSges(2,1) * t191 + mrSges(2,3) * t154 - pkin(1) * t83 + t177;
t74 = -mrSges(2,2) * t191 - mrSges(2,3) * t153 - qJ(2) * t83 - t165 * t79 + t167 * t77;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t168 * t74 - t166 * t75 - qJ(1) * (t166 * t81 + t168 * t86), t74, t77, t80, t92, t178; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t166 * t74 + t168 * t75 + qJ(1) * (-t166 * t86 + t168 * t81), t75, t79, t84, t91, (-t169 * t128 - t171 * t132) * qJD(3) + t181; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t180, t180, -t177, t176, t195, Ifges(6,5) * t147 + Ifges(6,6) * qJDD(4) - Ifges(6,3) * t148 - qJD(4) * t132 + t130 * t188 - t182;];
m_new = t1;
