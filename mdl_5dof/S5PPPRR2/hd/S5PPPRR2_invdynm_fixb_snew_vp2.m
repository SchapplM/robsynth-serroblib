% Calculate vector of cutting torques with Newton-Euler for
% S5PPPRR2
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
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S5PPPRR2_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR2_invdynm_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:29
% EndTime: 2019-12-05 14:59:31
% DurationCPUTime: 1.50s
% Computational Cost: add. (19323->147), mult. (29269->196), div. (0->0), fcn. (20964->10), ass. (0->74)
t128 = sin(qJ(5));
t130 = cos(qJ(5));
t107 = Ifges(6,6) * qJD(5) + (Ifges(6,4) * t128 + Ifges(6,2) * t130) * qJD(4);
t108 = Ifges(6,5) * qJD(5) + (Ifges(6,1) * t128 + Ifges(6,4) * t130) * qJD(4);
t141 = qJD(4) * qJD(5);
t113 = t128 * qJDD(4) + t130 * t141;
t114 = t130 * qJDD(4) - t128 * t141;
t124 = sin(pkin(7));
t127 = cos(pkin(7));
t118 = -t127 * g(1) - t124 * g(2);
t121 = -g(3) + qJDD(1);
t123 = sin(pkin(8));
t126 = cos(pkin(8));
t105 = t126 * t118 + t123 * t121;
t117 = t124 * g(1) - t127 * g(2);
t116 = qJDD(2) - t117;
t122 = sin(pkin(9));
t125 = cos(pkin(9));
t100 = t122 * t105 - t125 * t116;
t132 = qJD(4) ^ 2;
t101 = t125 * t105 + t122 * t116;
t104 = -t123 * t118 + t126 * t121;
t103 = qJDD(3) - t104;
t129 = sin(qJ(4));
t131 = cos(qJ(4));
t98 = t131 * t101 + t129 * t103;
t96 = -t132 * pkin(4) + qJDD(4) * pkin(6) + t98;
t93 = t130 * t100 - t128 * t96;
t94 = t128 * t100 + t130 * t96;
t144 = mrSges(6,1) * t93 - mrSges(6,2) * t94 + Ifges(6,5) * t113 + Ifges(6,6) * t114 + Ifges(6,3) * qJDD(5) + (t107 * t128 - t108 * t130) * qJD(4);
t143 = qJD(4) * t128;
t142 = qJD(4) * t130;
t112 = (-mrSges(6,1) * t130 + mrSges(6,2) * t128) * qJD(4);
t120 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t142;
t90 = m(6) * t93 + qJDD(5) * mrSges(6,1) - t113 * mrSges(6,3) + qJD(5) * t120 - t112 * t143;
t119 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t143;
t91 = m(6) * t94 - qJDD(5) * mrSges(6,2) + t114 * mrSges(6,3) - qJD(5) * t119 + t112 * t142;
t85 = -t128 * t90 + t130 * t91;
t82 = m(5) * t98 - t132 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t85;
t97 = -t129 * t101 + t131 * t103;
t95 = -qJDD(4) * pkin(4) - t132 * pkin(6) - t97;
t92 = -m(6) * t95 + t114 * mrSges(6,1) - t113 * mrSges(6,2) - t119 * t143 + t120 * t142;
t88 = m(5) * t97 + qJDD(4) * mrSges(5,1) - t132 * mrSges(5,2) + t92;
t79 = -t129 * t88 + t131 * t82;
t76 = m(4) * t101 + t79;
t84 = t128 * t91 + t130 * t90;
t80 = (-m(4) - m(5)) * t100 - t84;
t71 = -t122 * t80 + t125 * t76;
t68 = m(3) * t105 + t71;
t78 = t129 * t82 + t131 * t88;
t77 = -m(4) * t103 - t78;
t74 = m(3) * t104 + t77;
t140 = -t123 * t74 + t126 * t68;
t70 = t122 * t76 + t125 * t80;
t137 = -m(3) * t116 - t70;
t106 = Ifges(6,3) * qJD(5) + (Ifges(6,5) * t128 + Ifges(6,6) * t130) * qJD(4);
t86 = -mrSges(6,1) * t95 + mrSges(6,3) * t94 + Ifges(6,4) * t113 + Ifges(6,2) * t114 + Ifges(6,6) * qJDD(5) + qJD(5) * t108 - t106 * t143;
t87 = mrSges(6,2) * t95 - mrSges(6,3) * t93 + Ifges(6,1) * t113 + Ifges(6,4) * t114 + Ifges(6,5) * qJDD(5) - qJD(5) * t107 + t106 * t142;
t72 = mrSges(5,2) * t100 - mrSges(5,3) * t97 + Ifges(5,5) * qJDD(4) - t132 * Ifges(5,6) - pkin(6) * t84 - t128 * t86 + t130 * t87;
t73 = -mrSges(5,1) * t100 + mrSges(5,3) * t98 + t132 * Ifges(5,5) + Ifges(5,6) * qJDD(4) - pkin(4) * t84 - t144;
t61 = mrSges(4,2) * t103 + mrSges(4,3) * t100 - pkin(5) * t78 - t129 * t73 + t131 * t72;
t133 = mrSges(5,1) * t97 - mrSges(5,2) * t98 + Ifges(5,3) * qJDD(4) + pkin(4) * t92 + pkin(6) * t85 + t128 * t87 + t130 * t86;
t62 = -mrSges(4,1) * t103 + mrSges(4,3) * t101 - pkin(3) * t78 - t133;
t58 = mrSges(3,2) * t116 - mrSges(3,3) * t104 - qJ(3) * t70 - t122 * t62 + t125 * t61;
t135 = -mrSges(4,1) * t100 - mrSges(4,2) * t101 + pkin(3) * (-m(5) * t100 - t84) + pkin(5) * t79 + t129 * t72 + t131 * t73;
t60 = -mrSges(3,1) * t116 + mrSges(3,3) * t105 - pkin(2) * t70 - t135;
t138 = mrSges(2,1) * t117 - mrSges(2,2) * t118 + pkin(1) * t137 + qJ(2) * t140 + t123 * t58 + t126 * t60;
t134 = mrSges(3,1) * t104 - mrSges(3,2) * t105 + pkin(2) * t77 + qJ(3) * t71 + t122 * t61 + t125 * t62;
t67 = m(2) * t117 + t137;
t65 = t123 * t68 + t126 * t74;
t63 = m(2) * t118 + t140;
t56 = -mrSges(2,1) * t121 + mrSges(2,3) * t118 - pkin(1) * t65 - t134;
t55 = mrSges(2,2) * t121 - mrSges(2,3) * t117 - qJ(2) * t65 - t123 * t60 + t126 * t58;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t127 * t55 - t124 * t56 - qJ(1) * (t124 * t63 + t127 * t67), t55, t58, t61, t72, t87; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t124 * t55 + t127 * t56 + qJ(1) * (-t124 * t67 + t127 * t63), t56, t60, t62, t73, t86; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t138, t138, t134, t135, t133, t144;];
m_new = t1;
