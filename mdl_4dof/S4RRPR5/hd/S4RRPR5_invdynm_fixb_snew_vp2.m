% Calculate vector of cutting torques with Newton-Euler for
% S4RRPR5
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function m_new = S4RRPR5_invdynm_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_invdynm_fixb_snew_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_m_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:29
% EndTime: 2019-12-31 17:03:30
% DurationCPUTime: 0.54s
% Computational Cost: add. (6519->150), mult. (7778->183), div. (0->0), fcn. (3398->6), ass. (0->66)
t150 = -pkin(2) - pkin(6);
t149 = mrSges(3,1) - mrSges(4,2);
t148 = -Ifges(4,4) + Ifges(3,5);
t147 = Ifges(4,5) - Ifges(3,6);
t127 = sin(qJ(2));
t130 = cos(qJ(2));
t123 = qJD(1) + qJD(2);
t121 = t123 ^ 2;
t122 = qJDD(1) + qJDD(2);
t126 = sin(qJ(4));
t129 = cos(qJ(4));
t102 = (mrSges(5,1) * t126 + mrSges(5,2) * t129) * t123;
t144 = qJD(4) * t123;
t104 = t129 * t122 - t126 * t144;
t146 = t123 * t126;
t110 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t146;
t145 = t123 * t129;
t128 = sin(qJ(1));
t131 = cos(qJ(1));
t113 = t128 * g(1) - t131 * g(2);
t108 = qJDD(1) * pkin(1) + t113;
t114 = -t131 * g(1) - t128 * g(2);
t132 = qJD(1) ^ 2;
t109 = -t132 * pkin(1) + t114;
t90 = t130 * t108 - t127 * t109;
t141 = -t121 * qJ(3) + qJDD(3) - t90;
t85 = t150 * t122 + t141;
t81 = t126 * g(3) + t129 * t85;
t78 = m(5) * t81 + qJDD(4) * mrSges(5,1) - t104 * mrSges(5,3) + qJD(4) * t110 - t102 * t145;
t103 = -t126 * t122 - t129 * t144;
t111 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t145;
t82 = -t129 * g(3) + t126 * t85;
t79 = m(5) * t82 - qJDD(4) * mrSges(5,2) + t103 * mrSges(5,3) - qJD(4) * t111 - t102 * t146;
t67 = t126 * t79 + t129 * t78;
t88 = -t122 * pkin(2) + t141;
t139 = -m(4) * t88 + t121 * mrSges(4,3) - t67;
t64 = m(3) * t90 - t121 * mrSges(3,2) + t149 * t122 + t139;
t91 = t127 * t108 + t130 * t109;
t142 = t122 * qJ(3) + 0.2e1 * qJD(3) * t123 + t91;
t84 = t150 * t121 + t142;
t76 = -m(5) * t84 + t103 * mrSges(5,1) - t104 * mrSges(5,2) - t110 * t146 - t111 * t145;
t86 = t121 * pkin(2) - t142;
t137 = -m(4) * t86 + t121 * mrSges(4,2) + t122 * mrSges(4,3) - t76;
t71 = m(3) * t91 - t121 * mrSges(3,1) - t122 * mrSges(3,2) + t137;
t62 = t127 * t71 + t130 * t64;
t68 = -t126 * t78 + t129 * t79;
t143 = -t127 * t64 + t130 * t71;
t95 = Ifges(5,6) * qJD(4) + (Ifges(5,4) * t129 - Ifges(5,2) * t126) * t123;
t96 = Ifges(5,5) * qJD(4) + (Ifges(5,1) * t129 - Ifges(5,4) * t126) * t123;
t140 = mrSges(5,1) * t81 - mrSges(5,2) * t82 + Ifges(5,5) * t104 + Ifges(5,6) * t103 + Ifges(5,3) * qJDD(4) + t95 * t145 + t96 * t146;
t94 = Ifges(5,3) * qJD(4) + (Ifges(5,5) * t129 - Ifges(5,6) * t126) * t123;
t73 = -mrSges(5,1) * t84 + mrSges(5,3) * t82 + Ifges(5,4) * t104 + Ifges(5,2) * t103 + Ifges(5,6) * qJDD(4) + qJD(4) * t96 - t94 * t145;
t74 = mrSges(5,2) * t84 - mrSges(5,3) * t81 + Ifges(5,1) * t104 + Ifges(5,4) * t103 + Ifges(5,5) * qJDD(4) - qJD(4) * t95 - t94 * t146;
t138 = mrSges(4,2) * t88 - mrSges(4,3) * t86 + Ifges(4,1) * t122 - pkin(6) * t67 - t126 * t73 + t129 * t74;
t136 = -mrSges(4,1) * t86 - pkin(3) * t76 - pkin(6) * t68 - t126 * t74 - t129 * t73;
t135 = mrSges(4,1) * t88 + pkin(3) * t67 + t140;
t134 = -mrSges(3,2) * t91 + Ifges(3,3) * t122 + t138 + pkin(2) * (-t122 * mrSges(4,2) + t139) + qJ(3) * t137 + mrSges(3,1) * t90;
t133 = mrSges(2,1) * t113 - mrSges(2,2) * t114 + Ifges(2,3) * qJDD(1) + pkin(1) * t62 + t134;
t66 = -m(4) * g(3) + t68;
t60 = m(2) * t114 - t132 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t143;
t59 = m(2) * t113 + qJDD(1) * mrSges(2,1) - t132 * mrSges(2,2) + t62;
t58 = -mrSges(3,3) * t90 - qJ(3) * t66 + t148 * t122 + t147 * t121 + (-mrSges(3,2) + mrSges(4,3)) * g(3) + t135;
t57 = mrSges(3,3) * t91 - pkin(2) * t66 + t149 * g(3) + t148 * t121 - t147 * t122 + t136;
t56 = -mrSges(2,2) * g(3) - mrSges(2,3) * t113 + Ifges(2,5) * qJDD(1) - t132 * Ifges(2,6) - pkin(5) * t62 - t127 * t57 + t130 * t58;
t55 = Ifges(2,6) * qJDD(1) + t132 * Ifges(2,5) + mrSges(2,3) * t114 + t127 * t58 + t130 * t57 - pkin(1) * t68 + pkin(5) * t143 + (mrSges(2,1) - pkin(1) * (-m(3) - m(4))) * g(3);
t1 = [-mrSges(1,2) * g(3) + mrSges(1,3) * g(2) + t131 * t56 - t128 * t55 - pkin(4) * (t128 * t60 + t131 * t59), t56, t58, t138, t74; mrSges(1,1) * g(3) - mrSges(1,3) * g(1) + t128 * t56 + t131 * t55 + pkin(4) * (-t128 * t59 + t131 * t60), t55, t57, -mrSges(4,3) * g(3) + Ifges(4,4) * t122 - t121 * Ifges(4,5) - t135, t73; -mrSges(1,1) * g(2) + mrSges(1,2) * g(1) + t133, t133, t134, mrSges(4,2) * g(3) + t121 * Ifges(4,4) + Ifges(4,5) * t122 - t136, t140;];
m_new = t1;
