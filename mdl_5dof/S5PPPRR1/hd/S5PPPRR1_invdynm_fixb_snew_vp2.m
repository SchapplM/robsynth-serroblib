% Calculate vector of cutting torques with Newton-Euler for
% S5PPPRR1
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
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPPRR1_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:56
% EndTime: 2019-12-05 14:57:59
% DurationCPUTime: 1.68s
% Computational Cost: add. (21339->147), mult. (32645->196), div. (0->0), fcn. (23732->10), ass. (0->74)
t134 = sin(pkin(7));
t137 = cos(pkin(7));
t127 = -t137 * g(1) - t134 * g(2);
t131 = -g(3) + qJDD(1);
t133 = sin(pkin(8));
t136 = cos(pkin(8));
t115 = t136 * t127 + t133 * t131;
t126 = t134 * g(1) - t137 * g(2);
t125 = qJDD(2) - t126;
t132 = sin(pkin(9));
t135 = cos(pkin(9));
t111 = -t132 * t115 + t135 * t125;
t112 = t135 * t115 + t132 * t125;
t139 = sin(qJ(4));
t141 = cos(qJ(4));
t108 = t139 * t111 + t141 * t112;
t142 = qJD(4) ^ 2;
t105 = -t142 * pkin(4) + qJDD(4) * pkin(6) + t108;
t114 = -t133 * t127 + t136 * t131;
t113 = qJDD(3) - t114;
t138 = sin(qJ(5));
t140 = cos(qJ(5));
t102 = -t138 * t105 + t140 * t113;
t103 = t140 * t105 + t138 * t113;
t117 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t138 + Ifges(6,2) * t140) * qJD(4);
t118 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t138 + Ifges(6,4) * t140) * qJD(4);
t154 = qJD(4) * qJD(5);
t122 = qJDD(4) * t138 + t140 * t154;
t123 = qJDD(4) * t140 - t138 * t154;
t157 = mrSges(6,1) * t102 - mrSges(6,2) * t103 + Ifges(6,5) * t122 + Ifges(6,6) * t123 + Ifges(6,3) * qJDD(5) + (t117 * t138 - t118 * t140) * qJD(4);
t121 = (-mrSges(6,1) * t140 + mrSges(6,2) * t138) * qJD(4);
t155 = qJD(4) * t140;
t129 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t155;
t156 = qJD(4) * t138;
t98 = m(6) * t102 + qJDD(5) * mrSges(6,1) - mrSges(6,3) * t122 + qJD(5) * t129 - t121 * t156;
t128 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t156;
t99 = m(6) * t103 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t123 - qJD(5) * t128 + t121 * t155;
t152 = -t138 * t98 + t140 * t99;
t83 = m(5) * t108 - mrSges(5,1) * t142 - qJDD(4) * mrSges(5,2) + t152;
t107 = t111 * t141 - t112 * t139;
t104 = -qJDD(4) * pkin(4) - pkin(6) * t142 - t107;
t146 = -m(6) * t104 + t123 * mrSges(6,1) - mrSges(6,2) * t122 - t128 * t156 + t129 * t155;
t94 = m(5) * t107 + qJDD(4) * mrSges(5,1) - mrSges(5,2) * t142 + t146;
t80 = t139 * t83 + t141 * t94;
t87 = t138 * t99 + t140 * t98;
t78 = m(4) * t111 + t80;
t151 = -t139 * t94 + t141 * t83;
t79 = m(4) * t112 + t151;
t74 = -t132 * t78 + t135 * t79;
t71 = m(3) * t115 + t74;
t85 = (-m(4) - m(5)) * t113 - t87;
t84 = m(3) * t114 + t85;
t153 = -t133 * t84 + t136 * t71;
t73 = t132 * t79 + t135 * t78;
t148 = -m(3) * t125 - t73;
t116 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t138 + Ifges(6,6) * t140) * qJD(4);
t91 = -mrSges(6,1) * t104 + mrSges(6,3) * t103 + Ifges(6,4) * t122 + Ifges(6,2) * t123 + Ifges(6,6) * qJDD(5) + qJD(5) * t118 - t116 * t156;
t92 = mrSges(6,2) * t104 - mrSges(6,3) * t102 + Ifges(6,1) * t122 + Ifges(6,4) * t123 + Ifges(6,5) * qJDD(5) - qJD(5) * t117 + t116 * t155;
t75 = mrSges(5,2) * t113 - mrSges(5,3) * t107 + Ifges(5,5) * qJDD(4) - Ifges(5,6) * t142 - pkin(6) * t87 - t138 * t91 + t140 * t92;
t76 = -mrSges(5,1) * t113 + mrSges(5,3) * t108 + t142 * Ifges(5,5) + Ifges(5,6) * qJDD(4) - pkin(4) * t87 - t157;
t64 = -mrSges(4,1) * t113 + mrSges(4,3) * t112 + t139 * t75 + t141 * t76 - pkin(3) * (m(5) * t113 + t87) + pkin(5) * t151;
t65 = mrSges(4,2) * t113 - mrSges(4,3) * t111 - pkin(5) * t80 - t139 * t76 + t141 * t75;
t61 = mrSges(3,2) * t125 - mrSges(3,3) * t114 - qJ(3) * t73 - t132 * t64 + t135 * t65;
t147 = mrSges(5,1) * t107 - mrSges(5,2) * t108 + Ifges(5,3) * qJDD(4) + pkin(4) * t146 + pkin(6) * t152 + t138 * t92 + t140 * t91;
t143 = mrSges(4,1) * t111 - mrSges(4,2) * t112 + pkin(3) * t80 + t147;
t63 = -mrSges(3,1) * t125 + mrSges(3,3) * t115 - pkin(2) * t73 - t143;
t149 = mrSges(2,1) * t126 - mrSges(2,2) * t127 + pkin(1) * t148 + qJ(2) * t153 + t133 * t61 + t136 * t63;
t144 = mrSges(3,1) * t114 - mrSges(3,2) * t115 + pkin(2) * t85 + qJ(3) * t74 + t132 * t65 + t135 * t64;
t70 = m(2) * t126 + t148;
t68 = t133 * t71 + t136 * t84;
t66 = m(2) * t127 + t153;
t59 = -mrSges(2,1) * t131 + mrSges(2,3) * t127 - pkin(1) * t68 - t144;
t58 = mrSges(2,2) * t131 - mrSges(2,3) * t126 - qJ(2) * t68 - t133 * t63 + t136 * t61;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t137 * t58 - t134 * t59 - qJ(1) * (t134 * t66 + t137 * t70), t58, t61, t65, t75, t92; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t134 * t58 + t137 * t59 + qJ(1) * (-t134 * t70 + t137 * t66), t59, t63, t64, t76, t91; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t149, t149, t144, t143, t147, t157;];
m_new = t1;
