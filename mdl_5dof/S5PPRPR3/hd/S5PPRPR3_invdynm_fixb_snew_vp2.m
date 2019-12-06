% Calculate vector of cutting torques with Newton-Euler for
% S5PPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:52
% EndTime: 2019-12-05 15:04:57
% DurationCPUTime: 1.81s
% Computational Cost: add. (23172->158), mult. (35445->207), div. (0->0), fcn. (23732->10), ass. (0->76)
t144 = sin(pkin(7));
t147 = cos(pkin(7));
t136 = -t147 * g(1) - t144 * g(2);
t141 = -g(3) + qJDD(1);
t143 = sin(pkin(8));
t146 = cos(pkin(8));
t122 = t146 * t136 + t143 * t141;
t135 = t144 * g(1) - t147 * g(2);
t134 = qJDD(2) - t135;
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t117 = -t149 * t122 + t151 * t134;
t115 = qJDD(3) * pkin(3) + t117;
t118 = t151 * t122 + t149 * t134;
t152 = qJD(3) ^ 2;
t116 = -t152 * pkin(3) + t118;
t142 = sin(pkin(9));
t145 = cos(pkin(9));
t112 = t142 * t115 + t145 * t116;
t109 = -t152 * pkin(4) + qJDD(3) * pkin(6) + t112;
t121 = t143 * t136 - t146 * t141;
t120 = qJDD(4) + t121;
t148 = sin(qJ(5));
t150 = cos(qJ(5));
t106 = -t148 * t109 + t150 * t120;
t107 = t150 * t109 + t148 * t120;
t124 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t148 + Ifges(6,2) * t150) * qJD(3);
t125 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t148 + Ifges(6,4) * t150) * qJD(3);
t165 = qJD(3) * qJD(5);
t131 = t148 * qJDD(3) + t150 * t165;
t132 = t150 * qJDD(3) - t148 * t165;
t168 = mrSges(6,1) * t106 - mrSges(6,2) * t107 + Ifges(6,5) * t131 + Ifges(6,6) * t132 + Ifges(6,3) * qJDD(5) + (t124 * t148 - t125 * t150) * qJD(3);
t130 = (-mrSges(6,1) * t150 + mrSges(6,2) * t148) * qJD(3);
t166 = qJD(3) * t150;
t138 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t166;
t167 = qJD(3) * t148;
t102 = m(6) * t106 + qJDD(5) * mrSges(6,1) - t131 * mrSges(6,3) + qJD(5) * t138 - t130 * t167;
t137 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t167;
t103 = m(6) * t107 - qJDD(5) * mrSges(6,2) + t132 * mrSges(6,3) - qJD(5) * t137 + t130 * t166;
t162 = -t148 * t102 + t150 * t103;
t87 = m(5) * t112 - t152 * mrSges(5,1) - qJDD(3) * mrSges(5,2) + t162;
t111 = t145 * t115 - t142 * t116;
t108 = -qJDD(3) * pkin(4) - t152 * pkin(6) - t111;
t156 = -m(6) * t108 + t132 * mrSges(6,1) - t131 * mrSges(6,2) - t137 * t167 + t138 * t166;
t98 = m(5) * t111 + qJDD(3) * mrSges(5,1) - t152 * mrSges(5,2) + t156;
t84 = t142 * t87 + t145 * t98;
t91 = t150 * t102 + t148 * t103;
t164 = -t142 * t98 + t145 * t87;
t82 = m(4) * t117 + qJDD(3) * mrSges(4,1) - t152 * mrSges(4,2) + t84;
t83 = m(4) * t118 - t152 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t164;
t78 = -t149 * t82 + t151 * t83;
t75 = m(3) * t122 + t78;
t161 = m(5) * t120 + t91;
t88 = (-m(3) - m(4)) * t121 - t161;
t163 = -t143 * t88 + t146 * t75;
t77 = t149 * t83 + t151 * t82;
t158 = -m(3) * t134 - t77;
t123 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t148 + Ifges(6,6) * t150) * qJD(3);
t95 = -mrSges(6,1) * t108 + mrSges(6,3) * t107 + Ifges(6,4) * t131 + Ifges(6,2) * t132 + Ifges(6,6) * qJDD(5) + qJD(5) * t125 - t123 * t167;
t96 = mrSges(6,2) * t108 - mrSges(6,3) * t106 + Ifges(6,1) * t131 + Ifges(6,4) * t132 + Ifges(6,5) * qJDD(5) - qJD(5) * t124 + t123 * t166;
t79 = mrSges(5,2) * t120 - mrSges(5,3) * t111 + Ifges(5,5) * qJDD(3) - t152 * Ifges(5,6) - pkin(6) * t91 - t148 * t95 + t150 * t96;
t80 = -mrSges(5,1) * t120 + mrSges(5,3) * t112 + t152 * Ifges(5,5) + Ifges(5,6) * qJDD(3) - pkin(4) * t91 - t168;
t68 = -mrSges(4,1) * t121 + mrSges(4,3) * t118 + t152 * Ifges(4,5) + Ifges(4,6) * qJDD(3) - pkin(3) * t161 + qJ(4) * t164 + t142 * t79 + t145 * t80;
t69 = mrSges(4,2) * t121 - mrSges(4,3) * t117 + Ifges(4,5) * qJDD(3) - t152 * Ifges(4,6) - qJ(4) * t84 - t142 * t80 + t145 * t79;
t65 = mrSges(3,2) * t134 + mrSges(3,3) * t121 - pkin(5) * t77 - t149 * t68 + t151 * t69;
t157 = mrSges(5,1) * t111 - mrSges(5,2) * t112 + Ifges(5,3) * qJDD(3) + pkin(4) * t156 + pkin(6) * t162 + t148 * t96 + t150 * t95;
t153 = mrSges(4,1) * t117 - mrSges(4,2) * t118 + Ifges(4,3) * qJDD(3) + pkin(3) * t84 + t157;
t67 = -mrSges(3,1) * t134 + mrSges(3,3) * t122 - pkin(2) * t77 - t153;
t159 = mrSges(2,1) * t135 - mrSges(2,2) * t136 + pkin(1) * t158 + qJ(2) * t163 + t143 * t65 + t146 * t67;
t154 = -mrSges(3,1) * t121 - mrSges(3,2) * t122 + pkin(2) * (-m(4) * t121 - t161) + pkin(5) * t78 + t149 * t69 + t151 * t68;
t74 = m(2) * t135 + t158;
t72 = t143 * t75 + t146 * t88;
t70 = m(2) * t136 + t163;
t63 = -mrSges(2,1) * t141 + mrSges(2,3) * t136 - pkin(1) * t72 - t154;
t62 = mrSges(2,2) * t141 - mrSges(2,3) * t135 - qJ(2) * t72 - t143 * t67 + t146 * t65;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t147 * t62 - t144 * t63 - qJ(1) * (t144 * t70 + t147 * t74), t62, t65, t69, t79, t96; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t144 * t62 + t147 * t63 + qJ(1) * (-t144 * t74 + t147 * t70), t63, t67, t68, t80, t95; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t159, t159, t154, t153, t157, t168;];
m_new = t1;
