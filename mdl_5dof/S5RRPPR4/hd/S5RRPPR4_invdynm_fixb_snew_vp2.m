% Calculate vector of cutting torques with Newton-Euler for
% S5RRPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5RRPPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR4_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:43
% EndTime: 2019-12-31 19:27:45
% DurationCPUTime: 1.53s
% Computational Cost: add. (24384->182), mult. (28059->222), div. (0->0), fcn. (11856->8), ass. (0->79)
t160 = sin(qJ(1));
t163 = cos(qJ(1));
t137 = t160 * g(1) - t163 * g(2);
t132 = qJDD(1) * pkin(1) + t137;
t138 = -t163 * g(1) - t160 * g(2);
t165 = qJD(1) ^ 2;
t133 = -t165 * pkin(1) + t138;
t159 = sin(qJ(2));
t162 = cos(qJ(2));
t120 = t159 * t132 + t162 * t133;
t151 = qJDD(1) + qJDD(2);
t152 = (qJD(1) + qJD(2));
t177 = t151 * qJ(3) + (2 * qJD(3) * t152) + t120;
t183 = (-pkin(2) - pkin(3));
t184 = t152 ^ 2;
t110 = (t183 * t184) + t177;
t119 = t162 * t132 - t159 * t133;
t174 = -qJ(3) * t184 + qJDD(3) - t119;
t114 = t183 * t151 + t174;
t156 = sin(pkin(8));
t157 = cos(pkin(8));
t108 = t157 * t110 + t156 * t114;
t105 = -(pkin(4) * t184) - t151 * pkin(7) + t108;
t155 = g(3) + qJDD(4);
t158 = sin(qJ(5));
t161 = cos(qJ(5));
t102 = -t158 * t105 + t161 * t155;
t103 = t161 * t105 + t158 * t155;
t122 = Ifges(6,6) * qJD(5) + (-Ifges(6,4) * t158 - Ifges(6,2) * t161) * t152;
t123 = Ifges(6,5) * qJD(5) + (-Ifges(6,1) * t158 - Ifges(6,4) * t161) * t152;
t179 = qJD(5) * t152;
t127 = -t158 * t151 - t161 * t179;
t128 = -t161 * t151 + t158 * t179;
t185 = mrSges(6,1) * t102 - mrSges(6,2) * t103 + Ifges(6,5) * t127 + Ifges(6,6) * t128 + Ifges(6,3) * qJDD(5) - (t122 * t158 - t123 * t161) * t152;
t182 = (mrSges(3,1) + mrSges(4,1));
t115 = -(pkin(2) * t184) + t177;
t126 = (mrSges(6,1) * t161 - mrSges(6,2) * t158) * t152;
t180 = t152 * t161;
t135 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t180;
t181 = t152 * t158;
t100 = m(6) * t102 + qJDD(5) * mrSges(6,1) - t127 * mrSges(6,3) + qJD(5) * t135 + t126 * t181;
t134 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t181;
t101 = m(6) * t103 - qJDD(5) * mrSges(6,2) + t128 * mrSges(6,3) - qJD(5) * t134 - t126 * t180;
t94 = -t158 * t100 + t161 * t101;
t90 = m(5) * t108 - (mrSges(5,1) * t184) + t151 * mrSges(5,2) + t94;
t107 = -t156 * t110 + t157 * t114;
t104 = t151 * pkin(4) - pkin(7) * t184 - t107;
t98 = -m(6) * t104 + t128 * mrSges(6,1) - t127 * mrSges(6,2) + t134 * t181 - t135 * t180;
t97 = m(5) * t107 - t151 * mrSges(5,1) - mrSges(5,2) * t184 + t98;
t88 = -t156 * t97 + t157 * t90;
t175 = m(4) * t115 + t151 * mrSges(4,3) + t88;
t81 = m(3) * t120 - t151 * mrSges(3,2) - (t182 * t184) + t175;
t117 = -t151 * pkin(2) + t174;
t87 = t156 * t90 + t157 * t97;
t173 = -m(4) * t117 + t151 * mrSges(4,1) + mrSges(4,3) * t184 - t87;
t83 = m(3) * t119 + t151 * mrSges(3,1) - mrSges(3,2) * t184 + t173;
t76 = t159 * t81 + t162 * t83;
t178 = -t159 * t83 + t162 * t81;
t93 = t161 * t100 + t158 * t101;
t92 = -m(5) * t155 - t93;
t121 = Ifges(6,3) * qJD(5) + (-Ifges(6,5) * t158 - Ifges(6,6) * t161) * t152;
t95 = -mrSges(6,1) * t104 + mrSges(6,3) * t103 + Ifges(6,4) * t127 + Ifges(6,2) * t128 + Ifges(6,6) * qJDD(5) + qJD(5) * t123 + t121 * t181;
t96 = mrSges(6,2) * t104 - mrSges(6,3) * t102 + Ifges(6,1) * t127 + Ifges(6,4) * t128 + Ifges(6,5) * qJDD(5) - qJD(5) * t122 - t121 * t180;
t78 = mrSges(5,2) * t155 - mrSges(5,3) * t107 - Ifges(5,5) * t151 - (Ifges(5,6) * t184) - pkin(7) * t93 - t158 * t95 + t161 * t96;
t86 = -mrSges(5,1) * t155 + mrSges(5,3) * t108 + Ifges(5,5) * t184 - Ifges(5,6) * t151 - pkin(4) * t93 - t185;
t172 = mrSges(4,2) * t117 + mrSges(4,3) * g(3) + Ifges(4,4) * t151 + (Ifges(4,6) * t184) - qJ(4) * t87 - t156 * t86 + t157 * t78;
t171 = mrSges(4,2) * t115 - pkin(3) * t92 - qJ(4) * t88 - t156 * t78 - t157 * t86;
t169 = mrSges(5,1) * t107 - mrSges(5,2) * t108 - Ifges(5,3) * t151 + pkin(4) * t98 + pkin(7) * t94 + t158 * t96 + t161 * t95;
t168 = -mrSges(4,1) * t117 + mrSges(4,3) * t115 + Ifges(4,2) * t151 - pkin(3) * t87 - t169;
t167 = -mrSges(3,2) * t120 + mrSges(3,1) * t119 + Ifges(3,3) * t151 + t168 + qJ(3) * (-mrSges(4,1) * t184 + t175) + pkin(2) * t173;
t166 = mrSges(2,1) * t137 - mrSges(2,2) * t138 + Ifges(2,3) * qJDD(1) + pkin(1) * t76 + t167;
t91 = -m(4) * g(3) + t92;
t74 = m(2) * t138 - t165 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t178;
t73 = m(2) * t137 + qJDD(1) * mrSges(2,1) - t165 * mrSges(2,2) + t76;
t72 = -mrSges(3,2) * g(3) - mrSges(3,3) * t119 + Ifges(3,5) * t151 - (Ifges(3,6) * t184) - qJ(3) * t91 + t172;
t71 = mrSges(3,3) * t120 - pkin(2) * t91 + (Ifges(3,6) - Ifges(4,6)) * t151 + (Ifges(4,4) + Ifges(3,5)) * t184 + t182 * g(3) + t171;
t70 = -mrSges(2,2) * g(3) - mrSges(2,3) * t137 + Ifges(2,5) * qJDD(1) - t165 * Ifges(2,6) - pkin(6) * t76 - t159 * t71 + t162 * t72;
t69 = Ifges(2,6) * qJDD(1) + t165 * Ifges(2,5) + mrSges(2,3) * t138 + t159 * t72 + t162 * t71 - pkin(1) * t92 + pkin(6) * t178 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t163 * t70 - t160 * t69 - pkin(5) * (t160 * t74 + t163 * t73), t70, t72, t172, t78, t96; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t160 * t70 + t163 * t69 + pkin(5) * (-t160 * t73 + t163 * t74), t69, t71, t168, t86, t95; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t166, t166, t167, -mrSges(4,1) * g(3) - Ifges(4,4) * t184 + Ifges(4,6) * t151 - t171, t169, t185;];
m_new = t1;
