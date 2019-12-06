% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-12-05 16:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:02:27
% EndTime: 2019-12-05 16:02:34
% DurationCPUTime: 2.18s
% Computational Cost: add. (1435->267), mult. (3432->392), div. (0->0), fcn. (2066->8), ass. (0->133)
t67 = sin(qJ(5));
t134 = Ifges(6,6) * t67;
t70 = cos(qJ(5));
t135 = Ifges(6,5) * t70;
t144 = t70 / 0.2e1;
t145 = -t67 / 0.2e1;
t114 = qJD(4) * t67;
t71 = cos(qJ(4));
t116 = qJD(2) * t71;
t53 = t70 * t116 + t114;
t148 = t53 / 0.2e1;
t138 = Ifges(6,4) * t53;
t113 = qJD(4) * t70;
t52 = -t67 * t116 + t113;
t68 = sin(qJ(4));
t118 = qJD(2) * t68;
t64 = qJD(5) + t118;
t17 = Ifges(6,2) * t52 + Ifges(6,6) * t64 + t138;
t48 = Ifges(6,4) * t52;
t18 = Ifges(6,1) * t53 + Ifges(6,5) * t64 + t48;
t119 = qJD(1) * t68;
t73 = -pkin(2) - pkin(7);
t65 = sin(pkin(5));
t120 = qJD(1) * t65;
t72 = cos(qJ(2));
t102 = t72 * t120;
t83 = qJD(3) - t102;
t44 = t73 * qJD(2) + t83;
t66 = cos(pkin(5));
t28 = -t66 * t119 + t44 * t71;
t21 = -qJD(4) * pkin(4) - t28;
t136 = Ifges(6,4) * t70;
t87 = -Ifges(6,2) * t67 + t136;
t137 = Ifges(6,4) * t67;
t89 = Ifges(6,1) * t70 - t137;
t90 = mrSges(6,1) * t67 + mrSges(6,2) * t70;
t127 = t66 * t71;
t29 = qJD(1) * t127 + t44 * t68;
t22 = qJD(4) * pkin(8) + t29;
t69 = sin(qJ(2));
t106 = t69 * t120;
t58 = pkin(4) * t68 - pkin(8) * t71 + qJ(3);
t39 = t58 * qJD(2) + t106;
t5 = -t22 * t67 + t39 * t70;
t6 = t22 * t70 + t39 * t67;
t92 = t5 * t70 + t6 * t67;
t166 = -t92 * mrSges(6,3) + t18 * t144 + t17 * t145 + (-t134 + t135) * t64 / 0.2e1 + t89 * t148 + t87 * t52 / 0.2e1 + t21 * t90;
t37 = -mrSges(6,2) * t64 + mrSges(6,3) * t52;
t38 = mrSges(6,1) * t64 - mrSges(6,3) * t53;
t158 = -m(6) * t92 - t67 * t37 - t70 * t38;
t109 = qJD(4) * qJD(2);
t100 = t71 * t109;
t108 = qJD(4) * qJD(5);
t111 = qJD(5) * t71;
t33 = t70 * t108 + (-t67 * t111 - t68 * t113) * qJD(2);
t19 = mrSges(6,1) * t100 - mrSges(6,3) * t33;
t34 = -t67 * t108 + (-t70 * t111 + t68 * t114) * qJD(2);
t20 = -mrSges(6,2) * t100 + mrSges(6,3) * t34;
t165 = qJD(5) * t158 - t67 * t19 + t70 * t20;
t112 = qJD(4) * t71;
t103 = t73 * t112;
t126 = t67 * t69;
t125 = t68 * t73;
t41 = t70 * t125 + t58 * t67;
t94 = pkin(4) * t71 + pkin(8) * t68;
t49 = t94 * qJD(4) + qJD(3);
t164 = -t41 * qJD(5) - t67 * t103 + t49 * t70 - (-t68 * t126 + t70 * t72) * t120;
t124 = t69 * t70;
t40 = -t67 * t125 + t58 * t70;
t163 = t40 * qJD(5) + t70 * t103 + t49 * t67 - (t68 * t124 + t67 * t72) * t120;
t122 = Ifges(5,5) * qJD(4);
t139 = Ifges(5,4) * t68;
t160 = qJD(2) / 0.2e1;
t56 = qJD(2) * qJ(3) + t106;
t162 = t56 * mrSges(5,2) + t122 / 0.2e1 + (t71 * Ifges(5,1) - t139) * t160 - t28 * mrSges(5,3) + t166;
t117 = qJD(2) * t69;
t105 = t65 * t117;
t14 = t44 * t112 + (-qJD(4) * t66 + t105) * t119;
t30 = (t49 + t102) * qJD(2);
t1 = t5 * qJD(5) + t14 * t70 + t30 * t67;
t2 = -t6 * qJD(5) - t14 * t67 + t30 * t70;
t93 = t1 * t70 - t2 * t67;
t159 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t33 + Ifges(6,6) * t34;
t123 = qJD(4) * mrSges(5,1) + mrSges(6,1) * t52 - mrSges(6,2) * t53 - mrSges(5,3) * t116;
t157 = m(6) * t21 - t123;
t155 = 2 * m(5);
t153 = t33 / 0.2e1;
t152 = t34 / 0.2e1;
t151 = -t52 / 0.2e1;
t149 = -t53 / 0.2e1;
t147 = -t64 / 0.2e1;
t110 = t29 * qJD(4);
t97 = qJD(1) * t105;
t15 = -t71 * t97 + t110;
t128 = t65 * t72;
t42 = t71 * t128 + t66 * t68;
t133 = t15 * t42;
t130 = t56 * t72;
t121 = Ifges(5,6) * qJD(4);
t115 = qJD(2) * t72;
t107 = t68 * t128;
t104 = t65 * t115;
t101 = -t15 * t73 / 0.2e1;
t99 = -t122 / 0.2e1;
t98 = -t121 / 0.2e1;
t96 = -t106 / 0.2e1;
t91 = t70 * mrSges(6,1) - t67 * mrSges(6,2);
t88 = Ifges(6,1) * t67 + t136;
t86 = Ifges(6,2) * t70 + t137;
t84 = Ifges(6,5) * t67 + Ifges(6,6) * t70;
t43 = -t107 + t127;
t26 = t65 * t124 - t43 * t67;
t27 = t65 * t126 + t43 * t70;
t50 = (qJD(3) + t102) * qJD(2);
t77 = t56 * t115 + t50 * t69;
t76 = t6 * mrSges(6,2) - t64 * Ifges(6,3) - t53 * Ifges(6,5) - t52 * Ifges(6,6) + t121 / 0.2e1 + (Ifges(5,4) * t71 - t68 * Ifges(5,2)) * t160 - t5 * mrSges(6,1) - t56 * mrSges(5,1);
t74 = qJD(2) ^ 2;
t63 = Ifges(6,3) * t100;
t61 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t118;
t55 = t94 * qJD(2);
t54 = (mrSges(5,1) * t68 + mrSges(5,2) * t71) * qJD(2);
t51 = -qJD(2) * pkin(2) + t83;
t47 = (mrSges(5,1) * t71 - mrSges(5,2) * t68) * t109;
t25 = -qJD(4) * t107 - t71 * t105 + t66 * t112;
t24 = -t42 * qJD(4) + t68 * t105;
t13 = t28 * t70 + t55 * t67;
t12 = -t28 * t67 + t55 * t70;
t9 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t8 = Ifges(6,1) * t33 + Ifges(6,4) * t34 + Ifges(6,5) * t100;
t7 = Ifges(6,4) * t33 + Ifges(6,2) * t34 + Ifges(6,6) * t100;
t4 = t26 * qJD(5) + t67 * t104 + t24 * t70;
t3 = -t27 * qJD(5) + t70 * t104 - t24 * t67;
t10 = [t26 * t19 + t27 * t20 + t24 * t61 + t3 * t38 + t4 * t37 + t42 * t9 - t123 * t25 + (-t42 * t68 - t43 * t71) * mrSges(5,3) * t109 + m(5) * (t14 * t43 + t24 * t29 - t25 * t28 + t133) + m(6) * (t1 * t27 + t2 * t26 + t21 * t25 + t3 * t5 + t4 * t6 + t133) + (t54 * t115 + t69 * t47 + m(5) * t77 + ((-mrSges(3,2) + mrSges(4,3)) * t72 + (-mrSges(3,1) + mrSges(4,2)) * t69) * t74 + (t51 * t117 - t72 * t97 + t77) * m(4)) * t65; qJ(3) * t47 + t40 * t19 + t41 * t20 + t83 * t54 + t164 * t38 + t163 * t37 + (t83 * qJD(2) + t50) * mrSges(4,3) + 0.2e1 * (-m(5) * t130 / 0.2e1 - (pkin(2) * t117 + t51 * t69 + t130) * m(4) / 0.2e1) * t120 + (-t61 * t106 + t63 / 0.2e1 + t50 * mrSges(5,1) - t14 * mrSges(5,3) + (t29 * t96 + t14 * t73 / 0.2e1) * t155 + (0.3e1 / 0.2e1 * Ifges(5,4) * t118 + t99 + (-m(5) * t28 + t157) * t73 - t162) * qJD(4) + t159) * t68 + (t8 * t144 + t7 * t145 + t89 * t153 + t87 * t152 + t50 * mrSges(5,2) - t73 * t9 + (mrSges(5,3) + t90) * t15 + (-t1 * t67 - t2 * t70) * mrSges(6,3) + (t28 * t96 + t101) * t155 + (t84 * t147 + t86 * t151 + t88 * t149 + t21 * t91 - t70 * t17 / 0.2e1 + t18 * t145 + (t5 * t67 - t6 * t70) * mrSges(6,3)) * qJD(5) + (t73 * t61 + t98 + (m(5) * t73 - mrSges(5,3)) * t29 + ((-0.3e1 / 0.2e1 * Ifges(5,4) + t135 / 0.2e1 - t134 / 0.2e1) * t71 + (0.3e1 / 0.2e1 * Ifges(5,2) + Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(5,1)) * t68) * qJD(2) - t76) * qJD(4) + t157 * t106) * t71 + (m(4) + m(5)) * (qJ(3) * t50 + qJD(3) * t56) + (t1 * t41 + 0.2e1 * t101 * t71 + t163 * t6 + t164 * t5 + t2 * t40) * m(6); -t74 * mrSges(4,3) + (-t9 + (t70 * t37 - t67 * t38 + t61) * qJD(4) + m(5) * (-t15 + t110) + m(6) * (t6 * t113 - t5 * t114 - t15)) * t71 + (-t123 * qJD(4) + m(5) * (-t28 * qJD(4) + t14) + m(6) * (qJD(4) * t21 + t93) + t165) * t68 + (-m(5) * t56 - t54 + (-t56 + t106) * m(4) + t158) * qJD(2); t88 * t153 + t86 * t152 + t7 * t144 + t67 * t8 / 0.2e1 - t28 * t61 - t13 * t37 - t12 * t38 - pkin(4) * t9 - t14 * mrSges(5,2) + t123 * t29 + (-mrSges(5,1) - t91) * t15 + t93 * mrSges(6,3) + t166 * qJD(5) + ((t98 + qJD(4) * t84 / 0.2e1 + t29 * mrSges(5,3) + Ifges(5,4) * t116 / 0.2e1 + t76) * t71 + (t99 + (-t139 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t71) * qJD(2) + t162) * t68) * qJD(2) + (-pkin(4) * t15 - t12 * t5 - t13 * t6 - t21 * t29) * m(6) + (m(6) * t93 + t165) * pkin(8); t63 - t21 * (mrSges(6,1) * t53 + mrSges(6,2) * t52) + (Ifges(6,1) * t52 - t138) * t149 + t17 * t148 + (Ifges(6,5) * t52 - Ifges(6,6) * t53) * t147 - t5 * t37 + t6 * t38 + (t5 * t52 + t53 * t6) * mrSges(6,3) + (-Ifges(6,2) * t53 + t18 + t48) * t151 + t159;];
tauc = t10(:);
