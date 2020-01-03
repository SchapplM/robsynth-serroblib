% Calculate vector of cutting torques with Newton-Euler for
% S4PRPR4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4PRPR4_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:57
% EndTime: 2019-12-31 16:21:58
% DurationCPUTime: 0.47s
% Computational Cost: add. (4029->137), mult. (6861->171), div. (0->0), fcn. (3398->6), ass. (0->62)
t143 = -pkin(2) - pkin(5);
t142 = mrSges(3,1) - mrSges(4,2);
t141 = -Ifges(4,4) + Ifges(3,5);
t140 = Ifges(4,5) - Ifges(3,6);
t122 = sin(qJ(2));
t124 = cos(qJ(2));
t125 = qJD(2) ^ 2;
t121 = sin(qJ(4));
t123 = cos(qJ(4));
t103 = (mrSges(5,1) * t121 + mrSges(5,2) * t123) * qJD(2);
t137 = qJD(2) * qJD(4);
t105 = t123 * qJDD(2) - t121 * t137;
t139 = qJD(2) * t121;
t109 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t139;
t138 = qJD(2) * t123;
t116 = -g(3) + qJDD(1);
t119 = sin(pkin(6));
t120 = cos(pkin(6));
t107 = t119 * g(1) - t120 * g(2);
t108 = -t120 * g(1) - t119 * g(2);
t88 = t124 * t107 - t122 * t108;
t135 = -t125 * qJ(3) + qJDD(3) - t88;
t83 = t143 * qJDD(2) + t135;
t79 = -t121 * t116 + t123 * t83;
t76 = m(5) * t79 + qJDD(4) * mrSges(5,1) - t105 * mrSges(5,3) + qJD(4) * t109 - t103 * t138;
t104 = -t121 * qJDD(2) - t123 * t137;
t110 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t138;
t80 = t123 * t116 + t121 * t83;
t77 = m(5) * t80 - qJDD(4) * mrSges(5,2) + t104 * mrSges(5,3) - qJD(4) * t110 - t103 * t139;
t65 = t121 * t77 + t123 * t76;
t86 = -qJDD(2) * pkin(2) + t135;
t132 = -m(4) * t86 + t125 * mrSges(4,3) - t65;
t62 = m(3) * t88 - t125 * mrSges(3,2) + t142 * qJDD(2) + t132;
t89 = t122 * t107 + t124 * t108;
t133 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t89;
t82 = t143 * t125 + t133;
t75 = -m(5) * t82 + t104 * mrSges(5,1) - t105 * mrSges(5,2) - t109 * t139 - t110 * t138;
t84 = t125 * pkin(2) - t133;
t130 = -m(4) * t84 + t125 * mrSges(4,2) + qJDD(2) * mrSges(4,3) - t75;
t69 = m(3) * t89 - t125 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t130;
t60 = t122 * t69 + t124 * t62;
t66 = -t121 * t76 + t123 * t77;
t136 = -t122 * t62 + t124 * t69;
t64 = m(4) * t116 + t66;
t93 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t123 - Ifges(5,2) * t121) * qJD(2);
t94 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t123 - Ifges(5,4) * t121) * qJD(2);
t134 = mrSges(5,1) * t79 - mrSges(5,2) * t80 + Ifges(5,5) * t105 + Ifges(5,6) * t104 + Ifges(5,3) * qJDD(4) + t93 * t138 + t94 * t139;
t92 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t123 - Ifges(5,6) * t121) * qJD(2);
t71 = -mrSges(5,1) * t82 + mrSges(5,3) * t80 + Ifges(5,4) * t105 + Ifges(5,2) * t104 + Ifges(5,6) * qJDD(4) + qJD(4) * t94 - t92 * t138;
t72 = mrSges(5,2) * t82 - mrSges(5,3) * t79 + Ifges(5,1) * t105 + Ifges(5,4) * t104 + Ifges(5,5) * qJDD(4) - qJD(4) * t93 - t92 * t139;
t131 = mrSges(4,2) * t86 - mrSges(4,3) * t84 + Ifges(4,1) * qJDD(2) - pkin(5) * t65 - t121 * t71 + t123 * t72;
t129 = -mrSges(4,1) * t84 - pkin(3) * t75 - pkin(5) * t66 - t121 * t72 - t123 * t71;
t128 = mrSges(4,1) * t86 + pkin(3) * t65 + t134;
t127 = -mrSges(3,2) * t89 + Ifges(3,3) * qJDD(2) + t131 + pkin(2) * (-qJDD(2) * mrSges(4,2) + t132) + qJ(3) * t130 + mrSges(3,1) * t88;
t126 = mrSges(2,1) * t107 - mrSges(2,2) * t108 + pkin(1) * t60 + t127;
t58 = m(2) * t108 + t136;
t57 = m(2) * t107 + t60;
t56 = -mrSges(3,3) * t88 - qJ(3) * t64 + t140 * t125 + (mrSges(3,2) - mrSges(4,3)) * t116 + t141 * qJDD(2) + t128;
t55 = mrSges(3,3) * t89 - pkin(2) * t64 - t140 * qJDD(2) - t142 * t116 + t141 * t125 + t129;
t54 = mrSges(2,2) * t116 - mrSges(2,3) * t107 - pkin(4) * t60 - t122 * t55 + t124 * t56;
t53 = -mrSges(2,1) * t116 + mrSges(2,3) * t108 + t122 * t56 + t124 * t55 - pkin(1) * (m(3) * t116 + t64) + pkin(4) * t136;
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t120 * t54 - t119 * t53 - qJ(1) * (t119 * t58 + t120 * t57), t54, t56, t131, t72; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t119 * t54 + t120 * t53 + qJ(1) * (-t119 * t57 + t120 * t58), t53, t55, mrSges(4,3) * t116 + Ifges(4,4) * qJDD(2) - t125 * Ifges(4,5) - t128, t71; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t126, t126, t127, -mrSges(4,2) * t116 + t125 * Ifges(4,4) + Ifges(4,5) * qJDD(2) - t129, t134;];
m_new = t1;
