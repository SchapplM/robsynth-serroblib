% Calculate vector of inverse dynamics joint torques for
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
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:28
% EndTime: 2019-12-31 17:03:30
% DurationCPUTime: 0.90s
% Computational Cost: add. (747->154), mult. (1028->204), div. (0->0), fcn. (405->8), ass. (0->87)
t64 = sin(qJ(4));
t148 = -t64 / 0.2e1;
t67 = cos(qJ(4));
t132 = t67 / 0.2e1;
t141 = mrSges(3,1) - mrSges(4,2);
t84 = mrSges(5,1) * t64 + mrSges(5,2) * t67;
t147 = mrSges(4,3) + t84;
t146 = -mrSges(5,3) - t141;
t145 = mrSges(3,2) - t147;
t122 = Ifges(5,4) * t67;
t123 = Ifges(5,4) * t64;
t113 = pkin(1) * qJD(1);
t65 = sin(qJ(2));
t103 = t65 * t113;
t60 = qJD(1) + qJD(2);
t33 = qJ(3) * t60 + t103;
t85 = mrSges(5,1) * t67 - mrSges(5,2) * t64;
t144 = ((-Ifges(5,1) * t64 - t122) * t132 + (-Ifges(5,2) * t67 - t123) * t148) * t60 + t33 * t85 + qJD(4) * (-Ifges(5,5) * t64 - Ifges(5,6) * t67) / 0.2e1;
t131 = pkin(1) * t65;
t119 = t64 * mrSges(5,3);
t32 = -qJD(4) * mrSges(5,2) - t119 * t60;
t117 = t67 * mrSges(5,3);
t34 = qJD(4) * mrSges(5,1) - t117 * t60;
t140 = (-t64 * t32 - t67 * t34) * t131;
t23 = t84 * t60;
t139 = -mrSges(4,3) * t60 - t23;
t59 = qJDD(1) + qJDD(2);
t138 = qJ(3) * t59 + qJD(3) * t60;
t107 = qJD(4) * t64;
t24 = -t107 * t60 + t59 * t67;
t106 = qJD(4) * t67;
t25 = -t106 * t60 - t59 * t64;
t137 = t67 * (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t24) + t64 * (-qJDD(4) * mrSges(5,2) + mrSges(5,3) * t25);
t63 = qJ(1) + qJ(2);
t56 = sin(t63);
t57 = cos(t63);
t136 = -g(1) * t56 + g(2) * t57;
t68 = cos(qJ(2));
t102 = t68 * t113;
t87 = qJD(3) - t102;
t135 = t145 * t56 + t146 * t57;
t70 = -pkin(2) - pkin(6);
t134 = t145 * t57 + (-m(5) * t70 - t146) * t56;
t22 = t60 * t70 + t87;
t112 = pkin(1) * qJD(2);
t101 = t65 * t112;
t129 = pkin(1) * t68;
t28 = -qJD(1) * t101 + qJDD(1) * t129;
t77 = qJDD(3) - t28;
t9 = t59 * t70 + t77;
t1 = -t107 * t22 + t67 * t9;
t2 = t106 * t22 + t64 * t9;
t88 = t1 * t67 + t2 * t64;
t133 = t88 * m(5) + t32 * t106 - t34 * t107 + t137;
t66 = sin(qJ(1));
t130 = pkin(1) * t66;
t69 = cos(qJ(1));
t58 = t69 * pkin(1);
t121 = t33 * t68;
t120 = t59 * mrSges(4,2);
t114 = t57 * pkin(2) + t56 * qJ(3);
t110 = qJD(2) * t68;
t104 = t58 + t114;
t100 = pkin(1) * t110;
t53 = -pkin(2) - t129;
t44 = t57 * qJ(3);
t97 = -pkin(2) * t56 + t44;
t91 = t60 * t102;
t89 = (t64 ^ 2 + t67 ^ 2) * t65 * t22;
t83 = Ifges(5,1) * t67 - t123;
t82 = -Ifges(5,2) * t64 + t122;
t29 = (qJD(1) * t110 + qJDD(1) * t65) * pkin(1);
t12 = t29 + t138;
t40 = qJD(3) + t100;
t50 = qJ(3) + t131;
t80 = t12 * t50 + t33 * t40;
t79 = t32 * t67 - t34 * t64;
t78 = t12 * qJ(3) + t33 * qJD(3);
t73 = -t33 * t60 + t136;
t17 = -pkin(2) * t59 + t77;
t20 = Ifges(5,6) * qJD(4) + t60 * t82;
t21 = Ifges(5,5) * qJD(4) + t60 * t83;
t71 = -t29 * mrSges(3,2) - t1 * t117 - t119 * t2 - t21 * t107 / 0.2e1 - t20 * t106 / 0.2e1 + t17 * mrSges(4,2) + t28 * mrSges(3,1) + (Ifges(5,4) * t24 + Ifges(5,2) * t25) * t148 + (Ifges(5,1) * t24 + Ifges(5,4) * t25) * t132 + t24 * t83 / 0.2e1 + t25 * t82 / 0.2e1 + (Ifges(3,3) + Ifges(4,1)) * t59 + t147 * t12 + (0.2e1 * Ifges(5,5) * t132 - Ifges(5,6) * t64) * qJDD(4) + t144 * qJD(4);
t51 = t57 * pkin(6);
t31 = -pkin(2) * t60 + t87;
t5 = -mrSges(5,1) * t25 + mrSges(5,2) * t24;
t3 = [m(5) * (t112 * t89 + t80) + m(4) * (t101 * t31 + t17 * t53 + t80) + t71 + t53 * t120 + m(3) * (t28 * t68 + t29 * t65) * pkin(1) + t50 * t5 + Ifges(2,3) * qJDD(1) - t139 * t40 - t140 * qJD(2) + (mrSges(3,1) * t129 - mrSges(3,2) * t131 + t50 * mrSges(4,3)) * t59 + t133 * (-pkin(6) + t53) + (-m(5) * (t51 + t104) - m(4) * t104 - m(3) * t58 - t69 * mrSges(2,1) + t66 * mrSges(2,2) + t135) * g(2) + (m(3) * t130 - m(4) * (t97 - t130) - m(5) * (t44 - t130) + t66 * mrSges(2,1) + t69 * mrSges(2,2) + t134) * g(1) + (-mrSges(3,2) * t100 - t141 * t101) * t60; mrSges(3,2) * t91 - pkin(2) * t120 + qJ(3) * t5 + t71 + t141 * t60 * t103 + t87 * t23 + t140 * qJD(1) + (-(t89 + t121) * t113 + t78) * m(5) + (-(t31 * t65 + t121) * t113 - pkin(2) * t17 + t78) * m(4) + (-m(5) * (t51 + t114) - m(4) * t114 + t135) * g(2) + (-m(4) * t97 - m(5) * t44 + t134) * g(1) + (-t91 + t138) * mrSges(4,3) + t133 * t70; t120 + t139 * t60 + t79 * qJD(4) + (t73 + t88) * m(5) + (t17 + t73) * m(4) + t137; t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t24 + Ifges(5,6) * t25 + Ifges(5,3) * qJDD(4) + g(3) * t84 - t79 * t22 + (t64 * t21 / 0.2e1 + t20 * t132 - t144) * t60 + t136 * t85;];
tau = t3;
