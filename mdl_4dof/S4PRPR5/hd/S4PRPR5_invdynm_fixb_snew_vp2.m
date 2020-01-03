% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR5
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:23:00
% EndTime: 2019-12-31 16:23:01
% DurationCPUTime: 0.78s
% Computational Cost: add. (7906->138), mult. (12884->182), div. (0->0), fcn. (7300->8), ass. (0->63)
t128 = sin(qJ(4));
t130 = cos(qJ(4));
t105 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t128 + Ifges(5,2) * t130) * qJD(2);
t106 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t128 + Ifges(5,4) * t130) * qJD(2);
t143 = qJD(2) * qJD(4);
t111 = t128 * qJDD(2) + t130 * t143;
t112 = t130 * qJDD(2) - t128 * t143;
t125 = sin(pkin(6));
t127 = cos(pkin(6));
t115 = t125 * g(1) - t127 * g(2);
t114 = qJDD(3) - t115;
t132 = qJD(2) ^ 2;
t116 = -t127 * g(1) - t125 * g(2);
t123 = -g(3) + qJDD(1);
t129 = sin(qJ(2));
t131 = cos(qJ(2));
t102 = -t129 * t116 + t131 * t123;
t100 = qJDD(2) * pkin(2) + t102;
t103 = t131 * t116 + t129 * t123;
t101 = -t132 * pkin(2) + t103;
t124 = sin(pkin(7));
t126 = cos(pkin(7));
t97 = t124 * t100 + t126 * t101;
t94 = -t132 * pkin(3) + qJDD(2) * pkin(5) + t97;
t91 = t130 * t114 - t128 * t94;
t92 = t128 * t114 + t130 * t94;
t146 = mrSges(5,1) * t91 - mrSges(5,2) * t92 + Ifges(5,5) * t111 + Ifges(5,6) * t112 + Ifges(5,3) * qJDD(4) + (t105 * t128 - t106 * t130) * qJD(2);
t110 = (-mrSges(5,1) * t130 + mrSges(5,2) * t128) * qJD(2);
t144 = qJD(2) * t130;
t118 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t144;
t145 = qJD(2) * t128;
t87 = m(5) * t91 + qJDD(4) * mrSges(5,1) - t111 * mrSges(5,3) + qJD(4) * t118 - t110 * t145;
t117 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t145;
t88 = m(5) * t92 - qJDD(4) * mrSges(5,2) + t112 * mrSges(5,3) - qJD(4) * t117 + t110 * t144;
t141 = -t128 * t87 + t130 * t88;
t72 = m(4) * t97 - t132 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t141;
t96 = t126 * t100 - t124 * t101;
t93 = -qJDD(2) * pkin(3) - t132 * pkin(5) - t96;
t135 = -m(5) * t93 + t112 * mrSges(5,1) - t111 * mrSges(5,2) - t117 * t145 + t118 * t144;
t83 = m(4) * t96 + qJDD(2) * mrSges(4,1) - t132 * mrSges(4,2) + t135;
t69 = t124 * t72 + t126 * t83;
t76 = t128 * t88 + t130 * t87;
t142 = -t124 * t83 + t126 * t72;
t67 = m(3) * t102 + qJDD(2) * mrSges(3,1) - t132 * mrSges(3,2) + t69;
t68 = m(3) * t103 - t132 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t142;
t140 = -t129 * t67 + t131 * t68;
t139 = m(4) * t114 + t76;
t104 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t128 + Ifges(5,6) * t130) * qJD(2);
t80 = -mrSges(5,1) * t93 + mrSges(5,3) * t92 + Ifges(5,4) * t111 + Ifges(5,2) * t112 + Ifges(5,6) * qJDD(4) + qJD(4) * t106 - t104 * t145;
t81 = mrSges(5,2) * t93 - mrSges(5,3) * t91 + Ifges(5,1) * t111 + Ifges(5,4) * t112 + Ifges(5,5) * qJDD(4) - qJD(4) * t105 + t104 * t144;
t64 = mrSges(4,2) * t114 - mrSges(4,3) * t96 + Ifges(4,5) * qJDD(2) - t132 * Ifges(4,6) - pkin(5) * t76 - t128 * t80 + t130 * t81;
t65 = -mrSges(4,1) * t114 + mrSges(4,3) * t97 + t132 * Ifges(4,5) + Ifges(4,6) * qJDD(2) - pkin(3) * t76 - t146;
t58 = mrSges(3,1) * t115 + mrSges(3,3) * t103 + t132 * Ifges(3,5) + Ifges(3,6) * qJDD(2) - pkin(2) * t139 + qJ(3) * t142 + t124 * t64 + t126 * t65;
t60 = -mrSges(3,2) * t115 - mrSges(3,3) * t102 + Ifges(3,5) * qJDD(2) - t132 * Ifges(3,6) - qJ(3) * t69 - t124 * t65 + t126 * t64;
t137 = -mrSges(2,2) * t116 + pkin(4) * t140 + t129 * t60 + t131 * t58 + mrSges(2,1) * t115 + pkin(1) * (m(3) * t115 - t139);
t136 = mrSges(4,1) * t96 - mrSges(4,2) * t97 + Ifges(4,3) * qJDD(2) + pkin(3) * t135 + pkin(5) * t141 + t128 * t81 + t130 * t80;
t133 = mrSges(3,1) * t102 - mrSges(3,2) * t103 + Ifges(3,3) * qJDD(2) + pkin(2) * t69 + t136;
t73 = (m(2) + m(3)) * t115 - t139;
t63 = t129 * t68 + t131 * t67;
t61 = m(2) * t116 + t140;
t56 = -mrSges(2,1) * t123 + mrSges(2,3) * t116 - pkin(1) * t63 - t133;
t55 = mrSges(2,2) * t123 - mrSges(2,3) * t115 - pkin(4) * t63 - t129 * t58 + t131 * t60;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t127 * t55 - t125 * t56 - qJ(1) * (t125 * t61 + t127 * t73), t55, t60, t64, t81; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t125 * t55 + t127 * t56 + qJ(1) * (-t125 * t73 + t127 * t61), t56, t58, t65, t80; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t137, t137, t133, t136, t146;];
m_new = t1;
