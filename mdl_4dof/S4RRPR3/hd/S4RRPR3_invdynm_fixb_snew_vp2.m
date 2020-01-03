% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR3
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR3_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:30
% EndTime: 2019-12-31 17:01:31
% DurationCPUTime: 0.98s
% Computational Cost: add. (14032->149), mult. (18278->194), div. (0->0), fcn. (9296->8), ass. (0->67)
t128 = qJD(1) + qJD(2);
t133 = sin(qJ(4));
t136 = cos(qJ(4));
t105 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t133 + Ifges(5,2) * t136) * t128;
t106 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t133 + Ifges(5,4) * t136) * t128;
t127 = qJDD(1) + qJDD(2);
t150 = qJD(4) * t128;
t110 = t133 * t127 + t136 * t150;
t111 = t136 * t127 - t133 * t150;
t130 = -g(3) + qJDD(3);
t126 = t128 ^ 2;
t135 = sin(qJ(1));
t138 = cos(qJ(1));
t120 = t135 * g(1) - t138 * g(2);
t115 = qJDD(1) * pkin(1) + t120;
t121 = -t138 * g(1) - t135 * g(2);
t139 = qJD(1) ^ 2;
t116 = -t139 * pkin(1) + t121;
t134 = sin(qJ(2));
t137 = cos(qJ(2));
t103 = t134 * t115 + t137 * t116;
t100 = -t126 * pkin(2) + t103;
t131 = sin(pkin(7));
t132 = cos(pkin(7));
t102 = t137 * t115 - t134 * t116;
t99 = t127 * pkin(2) + t102;
t96 = t132 * t100 + t131 * t99;
t93 = -t126 * pkin(3) + t127 * pkin(6) + t96;
t90 = t136 * t130 - t133 * t93;
t91 = t133 * t130 + t136 * t93;
t153 = mrSges(5,1) * t90 - mrSges(5,2) * t91 + Ifges(5,5) * t110 + Ifges(5,6) * t111 + Ifges(5,3) * qJDD(4) + (t105 * t133 - t106 * t136) * t128;
t109 = (-mrSges(5,1) * t136 + mrSges(5,2) * t133) * t128;
t151 = t128 * t136;
t118 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t151;
t152 = t128 * t133;
t88 = m(5) * t90 + qJDD(4) * mrSges(5,1) - t110 * mrSges(5,3) + qJD(4) * t118 - t109 * t152;
t117 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t152;
t89 = m(5) * t91 - qJDD(4) * mrSges(5,2) + t111 * mrSges(5,3) - qJD(4) * t117 + t109 * t151;
t147 = -t133 * t88 + t136 * t89;
t75 = m(4) * t96 - t126 * mrSges(4,1) - t127 * mrSges(4,2) + t147;
t95 = -t131 * t100 + t132 * t99;
t92 = -t127 * pkin(3) - t126 * pkin(6) - t95;
t143 = -m(5) * t92 + t111 * mrSges(5,1) - t110 * mrSges(5,2) - t117 * t152 + t118 * t151;
t83 = m(4) * t95 + t127 * mrSges(4,1) - t126 * mrSges(4,2) + t143;
t72 = t131 * t75 + t132 * t83;
t68 = m(3) * t102 + t127 * mrSges(3,1) - t126 * mrSges(3,2) + t72;
t148 = -t131 * t83 + t132 * t75;
t69 = m(3) * t103 - t126 * mrSges(3,1) - t127 * mrSges(3,2) + t148;
t63 = t134 * t69 + t137 * t68;
t77 = t133 * t89 + t136 * t88;
t149 = m(4) * t130 + t77;
t146 = -t134 * t68 + t137 * t69;
t104 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t133 + Ifges(5,6) * t136) * t128;
t80 = -mrSges(5,1) * t92 + mrSges(5,3) * t91 + Ifges(5,4) * t110 + Ifges(5,2) * t111 + Ifges(5,6) * qJDD(4) + qJD(4) * t106 - t104 * t152;
t81 = mrSges(5,2) * t92 - mrSges(5,3) * t90 + Ifges(5,1) * t110 + Ifges(5,4) * t111 + Ifges(5,5) * qJDD(4) - qJD(4) * t105 + t104 * t151;
t144 = mrSges(4,1) * t95 - mrSges(4,2) * t96 + Ifges(4,3) * t127 + pkin(3) * t143 + pkin(6) * t147 + t133 * t81 + t136 * t80;
t141 = mrSges(3,1) * t102 - mrSges(3,2) * t103 + Ifges(3,3) * t127 + pkin(2) * t72 + t144;
t140 = mrSges(2,1) * t120 - mrSges(2,2) * t121 + Ifges(2,3) * qJDD(1) + pkin(1) * t63 + t141;
t70 = -mrSges(4,1) * t130 + mrSges(4,3) * t96 + t126 * Ifges(4,5) + Ifges(4,6) * t127 - pkin(3) * t77 - t153;
t64 = mrSges(4,2) * t130 - mrSges(4,3) * t95 + Ifges(4,5) * t127 - t126 * Ifges(4,6) - pkin(6) * t77 - t133 * t80 + t136 * t81;
t61 = m(2) * t121 - t139 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t146;
t60 = m(2) * t120 + qJDD(1) * mrSges(2,1) - t139 * mrSges(2,2) + t63;
t59 = -mrSges(3,2) * g(3) - mrSges(3,3) * t102 + Ifges(3,5) * t127 - t126 * Ifges(3,6) - qJ(3) * t72 - t131 * t70 + t132 * t64;
t58 = mrSges(3,1) * g(3) + mrSges(3,3) * t103 + t126 * Ifges(3,5) + Ifges(3,6) * t127 - pkin(2) * t149 + qJ(3) * t148 + t131 * t64 + t132 * t70;
t57 = -mrSges(2,2) * g(3) - mrSges(2,3) * t120 + Ifges(2,5) * qJDD(1) - t139 * Ifges(2,6) - pkin(5) * t63 - t134 * t58 + t137 * t59;
t56 = Ifges(2,6) * qJDD(1) + t139 * Ifges(2,5) + mrSges(2,1) * g(3) + mrSges(2,3) * t121 + t134 * t59 + t137 * t58 - pkin(1) * (-m(3) * g(3) + t149) + pkin(5) * t146;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t138 * t57 - t135 * t56 - pkin(4) * (t135 * t61 + t138 * t60), t57, t59, t64, t81; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t135 * t57 + t138 * t56 + pkin(4) * (-t135 * t60 + t138 * t61), t56, t58, t70, t80; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t140, t140, t141, t144, t153;];
m_new = t1;
