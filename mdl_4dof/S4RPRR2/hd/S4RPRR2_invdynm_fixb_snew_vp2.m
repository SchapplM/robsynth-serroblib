% Calculate vector of cutting torques with Newton-Euler for
% S4RPRR2
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:09
% EndTime: 2019-12-31 16:48:10
% DurationCPUTime: 0.95s
% Computational Cost: add. (13076->149), mult. (18278->194), div. (0->0), fcn. (9296->8), ass. (0->67)
t128 = qJD(1) + qJD(3);
t134 = sin(qJ(4));
t137 = cos(qJ(4));
t106 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t134 + Ifges(5,2) * t137) * t128;
t107 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t134 + Ifges(5,4) * t137) * t128;
t127 = qJDD(1) + qJDD(3);
t151 = qJD(4) * t128;
t111 = t134 * t127 + t137 * t151;
t112 = t137 * t127 - t134 * t151;
t131 = -g(3) + qJDD(2);
t126 = t128 ^ 2;
t136 = sin(qJ(1));
t139 = cos(qJ(1));
t121 = t136 * g(1) - t139 * g(2);
t116 = qJDD(1) * pkin(1) + t121;
t122 = -t139 * g(1) - t136 * g(2);
t140 = qJD(1) ^ 2;
t117 = -t140 * pkin(1) + t122;
t132 = sin(pkin(7));
t133 = cos(pkin(7));
t103 = t133 * t116 - t132 * t117;
t100 = qJDD(1) * pkin(2) + t103;
t104 = t132 * t116 + t133 * t117;
t101 = -t140 * pkin(2) + t104;
t135 = sin(qJ(3));
t138 = cos(qJ(3));
t97 = t135 * t100 + t138 * t101;
t94 = -t126 * pkin(3) + t127 * pkin(6) + t97;
t91 = t137 * t131 - t134 * t94;
t92 = t134 * t131 + t137 * t94;
t154 = mrSges(5,1) * t91 - mrSges(5,2) * t92 + Ifges(5,5) * t111 + Ifges(5,6) * t112 + Ifges(5,3) * qJDD(4) + (t106 * t134 - t107 * t137) * t128;
t110 = (-mrSges(5,1) * t137 + mrSges(5,2) * t134) * t128;
t152 = t128 * t137;
t119 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t152;
t153 = t128 * t134;
t89 = m(5) * t91 + qJDD(4) * mrSges(5,1) - t111 * mrSges(5,3) + qJD(4) * t119 - t110 * t153;
t118 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t153;
t90 = m(5) * t92 - qJDD(4) * mrSges(5,2) + t112 * mrSges(5,3) - qJD(4) * t118 + t110 * t152;
t148 = -t134 * t89 + t137 * t90;
t76 = m(4) * t97 - t126 * mrSges(4,1) - t127 * mrSges(4,2) + t148;
t96 = t138 * t100 - t135 * t101;
t93 = -t127 * pkin(3) - t126 * pkin(6) - t96;
t144 = -m(5) * t93 + t112 * mrSges(5,1) - t111 * mrSges(5,2) - t118 * t153 + t119 * t152;
t84 = m(4) * t96 + t127 * mrSges(4,1) - t126 * mrSges(4,2) + t144;
t73 = t135 * t76 + t138 * t84;
t69 = m(3) * t103 + qJDD(1) * mrSges(3,1) - t140 * mrSges(3,2) + t73;
t147 = -t135 * t84 + t138 * t76;
t70 = m(3) * t104 - t140 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t147;
t64 = t132 * t70 + t133 * t69;
t78 = t134 * t90 + t137 * t89;
t150 = m(4) * t131 + t78;
t149 = -t132 * t69 + t133 * t70;
t105 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t134 + Ifges(5,6) * t137) * t128;
t81 = -mrSges(5,1) * t93 + mrSges(5,3) * t92 + Ifges(5,4) * t111 + Ifges(5,2) * t112 + Ifges(5,6) * qJDD(4) + qJD(4) * t107 - t105 * t153;
t82 = mrSges(5,2) * t93 - mrSges(5,3) * t91 + Ifges(5,1) * t111 + Ifges(5,4) * t112 + Ifges(5,5) * qJDD(4) - qJD(4) * t106 + t105 * t152;
t145 = mrSges(4,1) * t96 - mrSges(4,2) * t97 + Ifges(4,3) * t127 + pkin(3) * t144 + pkin(6) * t148 + t134 * t82 + t137 * t81;
t142 = mrSges(3,1) * t103 - mrSges(3,2) * t104 + Ifges(3,3) * qJDD(1) + pkin(2) * t73 + t145;
t141 = mrSges(2,1) * t121 - mrSges(2,2) * t122 + Ifges(2,3) * qJDD(1) + pkin(1) * t64 + t142;
t71 = -mrSges(4,1) * t131 + mrSges(4,3) * t97 + t126 * Ifges(4,5) + Ifges(4,6) * t127 - pkin(3) * t78 - t154;
t65 = mrSges(4,2) * t131 - mrSges(4,3) * t96 + Ifges(4,5) * t127 - t126 * Ifges(4,6) - pkin(6) * t78 - t134 * t81 + t137 * t82;
t62 = m(2) * t122 - t140 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t149;
t61 = m(2) * t121 + qJDD(1) * mrSges(2,1) - t140 * mrSges(2,2) + t64;
t60 = mrSges(3,2) * t131 - mrSges(3,3) * t103 + Ifges(3,5) * qJDD(1) - t140 * Ifges(3,6) - pkin(5) * t73 - t135 * t71 + t138 * t65;
t59 = -mrSges(3,1) * t131 + mrSges(3,3) * t104 + t140 * Ifges(3,5) + Ifges(3,6) * qJDD(1) - pkin(2) * t150 + pkin(5) * t147 + t135 * t65 + t138 * t71;
t58 = -mrSges(2,2) * g(3) - mrSges(2,3) * t121 + Ifges(2,5) * qJDD(1) - t140 * Ifges(2,6) - qJ(2) * t64 - t132 * t59 + t133 * t60;
t57 = Ifges(2,6) * qJDD(1) + t140 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t122 + t132 * t60 + t133 * t59 - pkin(1) * (m(3) * t131 + t150) + qJ(2) * t149;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t139 * t58 - t136 * t57 - pkin(4) * (t136 * t62 + t139 * t61), t58, t60, t65, t82; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t136 * t58 + t139 * t57 + pkin(4) * (-t136 * t61 + t139 * t62), t57, t59, t71, t81; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t141, t141, t142, t145, t154;];
m_new = t1;
