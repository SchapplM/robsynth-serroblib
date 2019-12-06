% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:38
% EndTime: 2019-12-05 15:46:43
% DurationCPUTime: 1.64s
% Computational Cost: add. (1581->222), mult. (3852->328), div. (0->0), fcn. (2610->8), ass. (0->118)
t86 = sin(pkin(9));
t80 = pkin(2) * t86 + pkin(6);
t132 = pkin(7) + t80;
t89 = sin(qJ(4));
t67 = t132 * t89;
t92 = cos(qJ(4));
t68 = t132 * t92;
t88 = sin(qJ(5));
t91 = cos(qJ(5));
t35 = -t67 * t88 + t68 * t91;
t104 = qJD(4) * t132;
t55 = t89 * t104;
t56 = t92 * t104;
t93 = cos(qJ(2));
t115 = qJD(1) * t93;
t90 = sin(qJ(2));
t116 = qJD(1) * t90;
t79 = t86 * t116;
t87 = cos(pkin(9));
t60 = t87 * t115 - t79;
t72 = t88 * t92 + t89 * t91;
t141 = -qJD(5) * t35 + t55 * t88 - t56 * t91 + t72 * t60;
t34 = -t67 * t91 - t68 * t88;
t98 = t88 * t89 - t91 * t92;
t140 = qJD(5) * t34 - t55 * t91 - t56 * t88 + t98 * t60;
t139 = -Ifges(5,1) / 0.2e1;
t113 = qJD(2) * t92;
t138 = -Ifges(5,4) * t113 / 0.2e1;
t70 = t86 * t93 + t87 * t90;
t29 = t98 * t70;
t85 = qJD(4) + qJD(5);
t137 = m(5) * t80 + mrSges(5,3);
t78 = qJD(2) * pkin(2) + t115;
t52 = t87 * t116 + t86 * t78;
t47 = qJD(2) * pkin(6) + t52;
t103 = pkin(7) * qJD(2) + t47;
t112 = qJD(3) * t89;
t31 = t103 * t92 + t112;
t64 = t98 * qJD(2);
t135 = -t64 / 0.2e1;
t65 = t72 * qJD(2);
t134 = t65 / 0.2e1;
t133 = pkin(2) * t87;
t131 = mrSges(6,3) * t64;
t130 = Ifges(5,4) * t89;
t69 = t86 * t90 - t87 * t93;
t61 = t69 * qJD(2);
t54 = qJD(1) * t61;
t124 = t54 * t89;
t41 = t47 * t92 + t112;
t19 = -qJD(4) * t41 + t124;
t129 = t19 * t89;
t128 = t31 * t88;
t127 = t31 * t91;
t84 = t92 * qJD(3);
t40 = -t47 * t89 + t84;
t126 = t40 * t89;
t59 = t70 * qJD(2);
t53 = qJD(1) * t59;
t125 = t53 * t69;
t123 = t65 * mrSges(6,3);
t122 = t65 * Ifges(6,4);
t42 = t85 * t98;
t36 = t42 * qJD(2);
t121 = t98 * t36;
t43 = t85 * t72;
t37 = t43 * qJD(2);
t120 = t72 * t37;
t119 = qJD(4) * t84 - t92 * t54;
t118 = Ifges(5,5) * qJD(4);
t117 = Ifges(5,6) * qJD(4);
t114 = qJD(2) * t89;
t111 = qJD(4) * t89;
t110 = qJD(5) * t88;
t109 = qJD(5) * t91;
t108 = pkin(4) * t111;
t107 = -pkin(4) * t92 - pkin(3);
t106 = t118 / 0.2e1;
t105 = -t117 / 0.2e1;
t51 = t78 * t87 - t79;
t46 = -qJD(2) * pkin(3) - t51;
t102 = t114 * t139 - t118 / 0.2e1 + t138 - t46 * mrSges(5,2);
t39 = mrSges(6,1) * t64 + mrSges(6,2) * t65;
t101 = -t39 + (mrSges(5,1) * t92 - mrSges(5,2) * t89 + mrSges(4,1)) * qJD(2);
t97 = t103 * t89;
t30 = t84 - t97;
t24 = qJD(4) * pkin(4) + t30;
t7 = t24 * t91 - t128;
t8 = t24 * t88 + t127;
t100 = t41 * t92 - t126;
t76 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t114;
t77 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t113;
t99 = -t89 * t76 + t92 * t77;
t96 = t41 * mrSges(5,3) + t117 / 0.2e1 + (Ifges(5,2) * t92 + t130) * qJD(2) / 0.2e1 - t46 * mrSges(5,1);
t58 = t70 * qJD(1);
t14 = -qJD(4) * t97 + t119;
t15 = -qJD(4) * t31 + t124;
t2 = qJD(5) * t7 + t14 * t91 + t15 * t88;
t26 = -t64 * Ifges(6,2) + t85 * Ifges(6,6) + t122;
t57 = Ifges(6,4) * t64;
t27 = t65 * Ifges(6,1) + t85 * Ifges(6,5) - t57;
t3 = -qJD(5) * t8 - t14 * t88 + t15 * t91;
t44 = t107 * qJD(2) - t51;
t95 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t7 * t131 + t26 * t134 - t44 * (mrSges(6,1) * t65 - mrSges(6,2) * t64) - t65 * (-Ifges(6,1) * t64 - t122) / 0.2e1 - t85 * (-Ifges(6,5) * t64 - Ifges(6,6) * t65) / 0.2e1 - Ifges(6,6) * t37 - Ifges(6,5) * t36 + (-Ifges(6,2) * t65 + t27 - t57) * t64 / 0.2e1;
t81 = -pkin(3) - t133;
t75 = t107 - t133;
t66 = (mrSges(5,1) * t89 + mrSges(5,2) * t92) * qJD(4) * qJD(2);
t49 = mrSges(6,1) * t85 - t123;
t48 = -mrSges(6,2) * t85 - t131;
t45 = (t58 + t108) * qJD(2);
t28 = t72 * t70;
t18 = -t47 * t111 + t119;
t13 = mrSges(6,1) * t37 - mrSges(6,2) * t36;
t10 = t30 * t91 - t128;
t9 = -t30 * t88 - t127;
t5 = t29 * t85 + t72 * t61;
t4 = -t43 * t70 + t61 * t98;
t1 = [t4 * t48 + t5 * t49 + (-mrSges(3,1) * t90 - mrSges(3,2) * t93) * qJD(2) ^ 2 + (t13 + t66) * t69 + (-t76 * t92 - t77 * t89) * t70 * qJD(4) + (-t28 * t36 + t29 * t37) * mrSges(6,3) - (-qJD(2) * mrSges(4,2) + t99) * t61 - t101 * t59 + m(4) * (-t51 * t59 - t52 * t61 - t54 * t70 + t125) + m(5) * (t46 * t59 + t125 - t100 * t61 + (t18 * t92 - t129 + (-t40 * t92 - t41 * t89) * qJD(4)) * t70) + m(6) * (-t2 * t29 - t28 * t3 + t4 * t8 + t44 * t59 + t45 * t69 + t5 * t7); t85 * (-Ifges(6,5) * t42 - Ifges(6,6) * t43) / 0.2e1 + t75 * t13 + t81 * t66 + t45 * (mrSges(6,1) * t98 + mrSges(6,2) * t72) + t44 * (mrSges(6,1) * t43 - mrSges(6,2) * t42) - t53 * mrSges(4,1) - t42 * t27 / 0.2e1 - t43 * t26 / 0.2e1 + t141 * t49 + t140 * t48 + (-t43 * t135 + t37 * t98) * Ifges(6,2) + (-t42 * t134 - t36 * t72) * Ifges(6,1) + (qJD(2) * t60 + t54) * mrSges(4,2) + t101 * t58 + (-t2 * t98 - t3 * t72 + t34 * t36 - t35 * t37 + t42 * t7 - t43 * t8) * mrSges(6,3) + (t53 * mrSges(5,2) - t19 * mrSges(5,3) + t60 * t76 + (-t80 * t77 + pkin(4) * t39 - 0.3e1 / 0.2e1 * Ifges(5,4) * t114 + t105 - t96) * qJD(4)) * t89 - m(5) * (-t60 * t126 + t46 * t58) + m(5) * (t53 * t81 + (-t111 * t41 - t129) * t80) + (-t43 * t134 - t42 * t135 - t120 + t121) * Ifges(6,4) + (-t53 * mrSges(5,1) + (-t80 * t76 + t106 - t137 * t40 + (0.3e1 / 0.2e1 * Ifges(5,4) * t92 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t89) * qJD(2) - t102) * qJD(4) + (-t41 * m(5) - t77) * t60 + t137 * t18) * t92 + (t2 * t35 + t3 * t34 + t45 * t75 + t140 * t8 + t141 * t7 + (t108 - t58) * t44) * m(6) + (t51 * t58 - t52 * t60 + (-t53 * t87 - t54 * t86) * pkin(2)) * m(4); -t42 * t48 - t43 * t49 + (-t120 - t121) * mrSges(6,3) + m(5) * (t18 * t89 + t19 * t92) + m(6) * (t2 * t72 - t3 * t98 - t42 * t8 - t43 * t7) + (m(5) * t100 + (-t89 ^ 2 - t92 ^ 2) * qJD(2) * mrSges(5,3) + t99) * qJD(4); -m(6) * (t10 * t8 + t7 * t9) + (m(6) * (t8 * t109 - t7 * t110 + t2 * t88 + t3 * t91) + t48 * t109 - t49 * t110 + (t36 * t91 - t37 * t88) * mrSges(6,3)) * pkin(4) + ((t40 * mrSges(5,3) + t102 + t106 + t138) * t92 + (t105 + (t130 / 0.2e1 + (Ifges(5,2) / 0.2e1 + t139) * t92) * qJD(2) + (-m(6) * t44 - t39) * pkin(4) + t96) * t89) * qJD(2) + t8 * t123 + t41 * t76 - t40 * t77 - t10 * t48 - t9 * t49 - t18 * mrSges(5,2) + t19 * mrSges(5,1) + t95; (t49 + t123) * t8 - t48 * t7 + t95;];
tauc = t1(:);
