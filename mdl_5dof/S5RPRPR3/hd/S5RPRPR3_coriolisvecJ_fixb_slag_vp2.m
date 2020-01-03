% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:57
% EndTime: 2020-01-03 11:36:02
% DurationCPUTime: 1.00s
% Computational Cost: add. (1433->156), mult. (2925->242), div. (0->0), fcn. (1620->8), ass. (0->97)
t73 = sin(pkin(9));
t70 = t73 ^ 2;
t75 = cos(pkin(9));
t71 = t75 ^ 2;
t134 = (t70 + t71) * mrSges(5,3);
t72 = qJD(1) + qJD(3);
t126 = t72 * t73;
t77 = sin(qJ(5));
t79 = cos(qJ(5));
t94 = mrSges(6,1) * t77 + mrSges(6,2) * t79;
t48 = t94 * t126;
t133 = pkin(1) * sin(pkin(8));
t108 = qJD(1) * t133;
t67 = cos(pkin(8)) * pkin(1) + pkin(2);
t61 = t67 * qJD(1);
t78 = sin(qJ(3));
t80 = cos(qJ(3));
t50 = t80 * t108 + t78 * t61;
t40 = t72 * qJ(4) + t50;
t20 = -t75 * qJD(2) + t73 * t40;
t129 = t20 * t73;
t21 = t73 * qJD(2) + t75 * t40;
t92 = t21 * t75 + t129;
t138 = -m(5) * t92 - t72 * t134 - t73 * t48;
t116 = qJD(4) * t75;
t122 = t75 * t79;
t64 = t78 * t108;
t49 = t80 * t61 - t64;
t118 = qJ(4) * t75;
t60 = -t75 * pkin(4) - t73 * pkin(7) - pkin(3);
t51 = -t77 * t118 + t79 * t60;
t137 = t51 * qJD(5) + t79 * t116 - t49 * t122 - t77 * t50;
t123 = t75 * t77;
t52 = t79 * t118 + t77 * t60;
t136 = -t52 * qJD(5) - t77 * t116 + t49 * t123 - t79 * t50;
t97 = qJD(4) - t49;
t117 = qJD(3) * t80;
t44 = -qJD(3) * t64 + t61 * t117;
t34 = t72 * qJD(4) + t44;
t128 = t34 * t70;
t121 = t80 * t133 + t78 * t67;
t105 = t70 * qJD(5) * t72;
t132 = mrSges(6,3) * t73;
t131 = Ifges(6,4) * t77;
t130 = Ifges(6,4) * t79;
t127 = t34 * t71;
t124 = t75 * t72;
t119 = qJ(4) * t34;
t114 = qJD(5) * t73;
t113 = qJD(5) * t77;
t112 = qJD(5) * t79;
t65 = t78 * t133;
t110 = t77 * t132;
t109 = t79 * t132;
t107 = mrSges(6,3) * t114;
t106 = mrSges(6,3) * t112;
t103 = t80 * t67 - t65;
t102 = t77 * t107;
t101 = t73 * t106;
t17 = t60 * t72 + t97;
t8 = t79 * t17 - t77 * t21;
t9 = t77 * t17 + t79 * t21;
t98 = -t77 * t8 + t79 * t9;
t96 = -t75 * mrSges(5,1) + t73 * mrSges(5,2);
t95 = mrSges(6,1) * t79 - mrSges(6,2) * t77;
t93 = -Ifges(6,5) * t77 - Ifges(6,6) * t79;
t63 = qJD(5) - t124;
t46 = -t63 * mrSges(6,2) - t72 * t110;
t47 = t63 * mrSges(6,1) - t72 * t109;
t91 = t79 * t46 - t77 * t47;
t54 = -qJD(3) * t65 + t67 * t117;
t41 = -t103 + t60;
t56 = qJ(4) + t121;
t11 = t56 * t122 + t77 * t41;
t10 = -t56 * t123 + t79 * t41;
t90 = t77 * (-Ifges(6,2) * t79 - t131);
t89 = t79 * (-Ifges(6,1) * t77 - t130);
t88 = (t79 * Ifges(6,1) - t131) * t73;
t87 = (-t77 * Ifges(6,2) + t130) * t73;
t86 = t93 * qJD(5);
t85 = t95 * t114;
t26 = Ifges(6,6) * t63 + t72 * t87;
t27 = Ifges(6,5) * t63 + t72 * t88;
t45 = t50 * qJD(3);
t3 = qJD(5) * t8 + t34 * t122 + t77 * t45;
t4 = -qJD(5) * t9 - t34 * t123 + t79 * t45;
t81 = -t44 * mrSges(4,2) + t3 * (t75 * mrSges(6,2) - t110) + t20 * t85 + t94 * t128 + t4 * (-t75 * mrSges(6,1) - t109) + t8 * t102 + (-t124 / 0.2e1 + t63 / 0.2e1) * t73 * t86 + (t96 - mrSges(4,1)) * t45 + (-t90 + t89) * t105 + (t128 + t127) * mrSges(5,3) - ((t26 + t72 * (-Ifges(6,6) * t75 + t87)) * t79 + (t27 + t72 * (-Ifges(6,5) * t75 + t88)) * t77) * t114 / 0.2e1;
t57 = t96 * t72;
t55 = t121 * qJD(3);
t53 = qJD(4) + t54;
t43 = t72 * t85;
t39 = -t72 * pkin(3) + t97;
t23 = t70 * t119;
t15 = t56 * t128;
t6 = -t11 * qJD(5) - t53 * t123 + t79 * t55;
t5 = t10 * qJD(5) + t53 * t122 + t77 * t55;
t1 = [m(5) * (t56 * t127 + t15 + t45 * (-pkin(3) - t103) + t39 * t55) + m(4) * (-t45 * t103 + t44 * t121 - t49 * t55 + t50 * t54) + m(6) * (t4 * t10 + t3 * t11 + t9 * t5 + t8 * t6 + t15) + t73 * t56 * t43 + t5 * t46 + t6 * t47 + t55 * t57 + t81 - t9 * t101 + (-t55 * mrSges(4,1) - t54 * mrSges(4,2) + t10 * t102 - t11 * t101) * t72 + (m(6) * t129 - t138) * t53; -t75 * t43 + (-t77 ^ 2 - t79 ^ 2) * mrSges(6,3) * t105 + (-t46 * t113 - t47 * t112 + m(6) * (-t8 * t112 - t9 * t113 + t3 * t79 - t34 * t75 - t4 * t77)) * t73; (t50 * mrSges(4,1) + t49 * mrSges(4,2) + (t51 * t77 - t52 * t79) * t107 + t97 * t134) * t72 + (qJ(4) * t43 - t9 * t106 + t97 * t48) * t73 + t136 * t47 + t137 * t46 - t50 * t57 + t81 + (t97 * t129 + t136 * t8 + t137 * t9 + t3 * t52 + t4 * t51 + t23) * m(6) + (-t45 * pkin(3) + t71 * t119 - t39 * t50 + t97 * t92 + t23) * m(5); t91 * qJD(5) + m(5) * t45 + m(6) * (t98 * qJD(5) + t3 * t77 + t4 * t79) + (-t91 * t75 - m(6) * (t9 * t122 - t8 * t123 + t129) + t138) * t72; t4 * mrSges(6,1) - t3 * mrSges(6,2) - t8 * t46 + t9 * t47 + (-t20 * t95 + t77 * t27 / 0.2e1 + t79 * t26 / 0.2e1 - t63 * t93 / 0.2e1 + (-t89 / 0.2e1 + t90 / 0.2e1) * t126 + t86 + t98 * mrSges(6,3)) * t126;];
tauc = t1(:);
