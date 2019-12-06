% Calculate vector of cutting torques with Newton-Euler for
% S5PPRPR2
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:03
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPRPR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:03:02
% EndTime: 2019-12-05 15:03:04
% DurationCPUTime: 1.01s
% Computational Cost: add. (11364->158), mult. (17690->196), div. (0->0), fcn. (10502->8), ass. (0->75)
t173 = -pkin(3) - pkin(6);
t172 = mrSges(4,1) - mrSges(5,2);
t171 = -Ifges(5,4) + Ifges(4,5);
t170 = Ifges(5,5) - Ifges(4,6);
t149 = sin(qJ(3));
t151 = cos(qJ(3));
t145 = sin(pkin(7));
t147 = cos(pkin(7));
t134 = -t147 * g(1) - t145 * g(2);
t141 = -g(3) + qJDD(1);
t144 = sin(pkin(8));
t146 = cos(pkin(8));
t115 = -t144 * t134 + t146 * t141;
t116 = t146 * t134 + t144 * t141;
t110 = t151 * t115 - t149 * t116;
t152 = qJD(3) ^ 2;
t162 = -t152 * qJ(4) + qJDD(4) - t110;
t108 = -qJDD(3) * pkin(3) + t162;
t148 = sin(qJ(5));
t150 = cos(qJ(5));
t105 = t173 * qJDD(3) + t162;
t133 = t145 * g(1) - t147 * g(2);
t132 = qJDD(2) - t133;
t101 = t150 * t105 - t148 * t132;
t128 = (mrSges(6,1) * t148 + mrSges(6,2) * t150) * qJD(3);
t167 = qJD(3) * qJD(5);
t130 = t150 * qJDD(3) - t148 * t167;
t169 = qJD(3) * t148;
t135 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t169;
t168 = qJD(3) * t150;
t97 = m(6) * t101 + qJDD(5) * mrSges(6,1) - t130 * mrSges(6,3) + qJD(5) * t135 - t128 * t168;
t102 = t148 * t105 + t150 * t132;
t129 = -t148 * qJDD(3) - t150 * t167;
t136 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t168;
t98 = m(6) * t102 - qJDD(5) * mrSges(6,2) + t129 * mrSges(6,3) - qJD(5) * t136 - t128 * t169;
t85 = t148 * t98 + t150 * t97;
t160 = -m(5) * t108 + t152 * mrSges(5,3) - t85;
t79 = m(4) * t110 - t152 * mrSges(4,2) + t172 * qJDD(3) + t160;
t111 = t149 * t115 + t151 * t116;
t161 = qJDD(3) * qJ(4) + 0.2e1 * qJD(4) * qJD(3) + t111;
t106 = t152 * pkin(3) - t161;
t104 = t173 * t152 + t161;
t99 = -m(6) * t104 + t129 * mrSges(6,1) - t130 * mrSges(6,2) - t135 * t169 - t136 * t168;
t157 = -m(5) * t106 + t152 * mrSges(5,2) + qJDD(3) * mrSges(5,3) - t99;
t89 = m(4) * t111 - t152 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t157;
t77 = t149 * t89 + t151 * t79;
t86 = -t148 * t97 + t150 * t98;
t84 = m(5) * t132 + t86;
t75 = m(3) * t115 + t77;
t165 = -t149 * t79 + t151 * t89;
t76 = m(3) * t116 + t165;
t166 = -t144 * t75 + t146 * t76;
t164 = (-m(3) - m(4)) * t132 - t84;
t119 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t150 - Ifges(6,6) * t148) * qJD(3);
t121 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t150 - Ifges(6,4) * t148) * qJD(3);
t91 = -mrSges(6,1) * t104 + mrSges(6,3) * t102 + Ifges(6,4) * t130 + Ifges(6,2) * t129 + Ifges(6,6) * qJDD(5) + qJD(5) * t121 - t119 * t168;
t120 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t150 - Ifges(6,2) * t148) * qJD(3);
t92 = mrSges(6,2) * t104 - mrSges(6,3) * t101 + Ifges(6,1) * t130 + Ifges(6,4) * t129 + Ifges(6,5) * qJDD(5) - qJD(5) * t120 - t119 * t169;
t156 = -mrSges(5,1) * t106 - pkin(4) * t99 - pkin(6) * t86 - t148 * t92 - t150 * t91;
t72 = mrSges(4,3) * t111 - pkin(3) * t84 - t170 * qJDD(3) - t172 * t132 + t171 * t152 + t156;
t159 = mrSges(6,1) * t101 - mrSges(6,2) * t102 + Ifges(6,5) * t130 + Ifges(6,6) * t129 + Ifges(6,3) * qJDD(5) + t120 * t168 + t121 * t169;
t155 = mrSges(5,1) * t108 + pkin(4) * t85 + t159;
t73 = -mrSges(4,3) * t110 - qJ(4) * t84 + t170 * t152 + (mrSges(4,2) - mrSges(5,3)) * t132 + t171 * qJDD(3) + t155;
t66 = -mrSges(3,1) * t132 + mrSges(3,3) * t116 + t149 * t73 + t151 * t72 - pkin(2) * (m(4) * t132 + t84) + pkin(5) * t165;
t68 = mrSges(3,2) * t132 - mrSges(3,3) * t115 - pkin(5) * t77 - t149 * t72 + t151 * t73;
t163 = mrSges(2,1) * t133 - mrSges(2,2) * t134 + pkin(1) * t164 + qJ(2) * t166 + t144 * t68 + t146 * t66;
t158 = mrSges(5,2) * t108 - mrSges(5,3) * t106 + Ifges(5,1) * qJDD(3) - pkin(6) * t85 - t148 * t91 + t150 * t92;
t154 = -mrSges(4,2) * t111 + mrSges(4,1) * t110 + Ifges(4,3) * qJDD(3) + t158 + pkin(3) * (-qJDD(3) * mrSges(5,2) + t160) + qJ(4) * t157;
t153 = mrSges(3,1) * t115 - mrSges(3,2) * t116 + pkin(2) * t77 + t154;
t82 = m(2) * t133 + t164;
t71 = t144 * t76 + t146 * t75;
t69 = m(2) * t134 + t166;
t64 = -mrSges(2,1) * t141 + mrSges(2,3) * t134 - pkin(1) * t71 - t153;
t63 = mrSges(2,2) * t141 - mrSges(2,3) * t133 - qJ(2) * t71 - t144 * t66 + t146 * t68;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t147 * t63 - t145 * t64 - qJ(1) * (t145 * t69 + t147 * t82), t63, t68, t73, t158, t92; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t145 * t63 + t147 * t64 + qJ(1) * (-t145 * t82 + t147 * t69), t64, t66, t72, mrSges(5,3) * t132 + Ifges(5,4) * qJDD(3) - t152 * Ifges(5,5) - t155, t91; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t163, t163, t153, t154, -mrSges(5,2) * t132 + t152 * Ifges(5,4) + Ifges(5,5) * qJDD(3) - t156, t159;];
m_new = t1;
