% Calculate vector of cutting torques with Newton-Euler for
% S5PRRRR2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
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
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PRRRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:46
% EndTime: 2019-12-05 17:04:48
% DurationCPUTime: 1.35s
% Computational Cost: add. (22976->160), mult. (25726->203), div. (0->0), fcn. (13492->8), ass. (0->74)
t156 = sin(qJ(2));
t160 = cos(qJ(2));
t137 = t156 * g(1) - t160 * g(2);
t134 = qJDD(2) * pkin(2) + t137;
t138 = -t160 * g(1) - t156 * g(2);
t163 = qJD(2) ^ 2;
t135 = -t163 * pkin(2) + t138;
t155 = sin(qJ(3));
t159 = cos(qJ(3));
t121 = t159 * t134 - t155 * t135;
t149 = qJDD(2) + qJDD(3);
t118 = t149 * pkin(3) + t121;
t122 = t155 * t134 + t159 * t135;
t150 = qJD(2) + qJD(3);
t148 = t150 ^ 2;
t119 = -t148 * pkin(3) + t122;
t154 = sin(qJ(4));
t158 = cos(qJ(4));
t115 = t154 * t118 + t158 * t119;
t142 = qJDD(4) + t149;
t111 = t142 * pkin(6) + t115;
t152 = -g(3) + qJDD(1);
t153 = sin(qJ(5));
t157 = cos(qJ(5));
t109 = -t153 * t111 + t157 * t152;
t110 = t157 * t111 + t153 * t152;
t143 = qJD(4) + t150;
t124 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t153 + Ifges(6,2) * t157) * t143;
t125 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t153 + Ifges(6,4) * t157) * t143;
t177 = qJD(5) * t143;
t127 = t153 * t142 + t157 * t177;
t128 = t157 * t142 - t153 * t177;
t182 = mrSges(6,1) * t109 - mrSges(6,2) * t110 + Ifges(6,5) * t127 + Ifges(6,6) * t128 + Ifges(6,3) * qJDD(5) + (t124 * t153 - t125 * t157) * t143;
t114 = t158 * t118 - t154 * t119;
t141 = t143 ^ 2;
t112 = -t141 * pkin(6) - t114;
t180 = t143 * t153;
t132 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t180;
t179 = t143 * t157;
t133 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t179;
t103 = m(5) * t114 - m(6) * t112 + t142 * mrSges(5,1) + t128 * mrSges(6,1) - t141 * mrSges(5,2) - t127 * mrSges(6,2) + (-t132 * t153 + t133 * t157) * t143;
t126 = (-mrSges(6,1) * t157 + mrSges(6,2) * t153) * t143;
t107 = m(6) * t109 + qJDD(5) * mrSges(6,1) - t127 * mrSges(6,3) + qJD(5) * t133 - t126 * t180;
t108 = m(6) * t110 - qJDD(5) * mrSges(6,2) + t128 * mrSges(6,3) - qJD(5) * t132 + t126 * t179;
t173 = -t153 * t107 + t157 * t108;
t94 = m(5) * t115 - t141 * mrSges(5,1) - t142 * mrSges(5,2) + t173;
t91 = t158 * t103 + t154 * t94;
t88 = m(4) * t121 + t149 * mrSges(4,1) - t148 * mrSges(4,2) + t91;
t174 = -t154 * t103 + t158 * t94;
t89 = m(4) * t122 - t148 * mrSges(4,1) - t149 * mrSges(4,2) + t174;
t83 = t155 * t89 + t159 * t88;
t80 = m(3) * t137 + qJDD(2) * mrSges(3,1) - t163 * mrSges(3,2) + t83;
t175 = -t155 * t88 + t159 * t89;
t81 = m(3) * t138 - t163 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t175;
t181 = t156 * t81 + t160 * t80;
t178 = t157 * t107 + t153 * t108;
t176 = m(5) * t152 + t178;
t172 = m(4) * t152 + t176;
t123 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t153 + Ifges(6,6) * t157) * t143;
t100 = -mrSges(6,1) * t112 + mrSges(6,3) * t110 + Ifges(6,4) * t127 + Ifges(6,2) * t128 + Ifges(6,6) * qJDD(5) + qJD(5) * t125 - t123 * t180;
t101 = mrSges(6,2) * t112 - mrSges(6,3) * t109 + Ifges(6,1) * t127 + Ifges(6,4) * t128 + Ifges(6,5) * qJDD(5) - qJD(5) * t124 + t123 * t179;
t84 = mrSges(5,2) * t152 - mrSges(5,3) * t114 + Ifges(5,5) * t142 - t141 * Ifges(5,6) - pkin(6) * t178 - t153 * t100 + t157 * t101;
t95 = -mrSges(5,1) * t152 + mrSges(5,3) * t115 + t141 * Ifges(5,5) + Ifges(5,6) * t142 - t182;
t76 = -mrSges(4,1) * t152 + mrSges(4,3) * t122 + t148 * Ifges(4,5) + Ifges(4,6) * t149 - pkin(3) * t176 + pkin(5) * t174 + t154 * t84 + t158 * t95;
t77 = mrSges(4,2) * t152 - mrSges(4,3) * t121 + Ifges(4,5) * t149 - t148 * Ifges(4,6) - pkin(5) * t91 - t154 * t95 + t158 * t84;
t71 = -mrSges(3,1) * t152 + mrSges(3,3) * t138 + t163 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t172 + pkin(4) * t175 + t155 * t77 + t159 * t76;
t74 = mrSges(3,2) * t152 - mrSges(3,3) * t137 + Ifges(3,5) * qJDD(2) - t163 * Ifges(3,6) - pkin(4) * t83 - t155 * t76 + t159 * t77;
t170 = mrSges(2,2) * t152 + mrSges(2,3) * g(2) - t156 * t71 + t160 * t74;
t169 = -mrSges(2,1) * t152 - pkin(1) * (m(3) * t152 + t172) + t156 * t74 + t160 * t71;
t168 = mrSges(5,1) * t114 - mrSges(5,2) * t115 + Ifges(5,3) * t142 + pkin(6) * t173 + t157 * t100 + t153 * t101;
t167 = mrSges(4,1) * t121 - mrSges(4,2) * t122 + Ifges(4,3) * t149 + pkin(3) * t91 + t168;
t165 = mrSges(3,1) * t137 - mrSges(3,2) * t138 + Ifges(3,3) * qJDD(2) + pkin(2) * t83 + t167;
t164 = mrSges(2,2) * g(1) + pkin(1) * t181 + t165;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) - qJ(1) * (-m(2) * g(2) + t181) + t170, t170, t74, t77, t84, t101; mrSges(1,1) * g(3) + qJ(1) * (-t156 * t80 + t160 * t81) + (-qJ(1) * m(2) - mrSges(1,3) - mrSges(2,3)) * g(1) + t169, -mrSges(2,3) * g(1) + t169, t71, t76, t95, t100; t164 + mrSges(1,2) * g(1) + (-mrSges(1,1) - mrSges(2,1)) * g(2), -mrSges(2,1) * g(2) + t164, t165, t167, t168, t182;];
m_new = t1;
