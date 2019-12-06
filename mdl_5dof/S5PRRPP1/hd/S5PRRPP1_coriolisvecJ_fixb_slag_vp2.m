% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
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
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:06:12
% EndTime: 2019-12-05 16:06:20
% DurationCPUTime: 2.21s
% Computational Cost: add. (999->237), mult. (2725->321), div. (0->0), fcn. (1611->4), ass. (0->112)
t133 = mrSges(6,1) + mrSges(5,1);
t132 = Ifges(5,1) + Ifges(6,1);
t131 = Ifges(6,4) + Ifges(5,5);
t73 = sin(qJ(3));
t91 = qJD(2) * qJD(3);
t84 = t73 * t91;
t105 = mrSges(6,2) + mrSges(5,3);
t104 = -Ifges(5,4) + Ifges(6,5);
t130 = Ifges(6,6) - Ifges(5,6);
t72 = sin(pkin(8));
t74 = cos(qJ(3));
t99 = cos(pkin(8));
t81 = t99 * t74;
t95 = qJD(2) * t73;
t47 = -qJD(2) * t81 + t72 * t95;
t114 = Ifges(6,5) * t47;
t43 = Ifges(5,4) * t47;
t82 = t99 * t73;
t56 = t72 * t74 + t82;
t49 = t56 * qJD(2);
t129 = t131 * qJD(3) + t132 * t49 + t114 - t43;
t111 = t47 * mrSges(5,3);
t112 = t47 * mrSges(6,2);
t33 = qJD(3) * mrSges(6,3) - t112;
t102 = -qJD(3) * mrSges(5,2) - t111 + t33;
t109 = t49 * mrSges(5,3);
t110 = t49 * mrSges(6,2);
t101 = t133 * qJD(3) - t109 - t110;
t128 = -2 * pkin(2);
t126 = -t47 / 0.2e1;
t125 = t47 / 0.2e1;
t123 = t49 / 0.2e1;
t121 = pkin(3) * t72;
t103 = -qJ(4) - pkin(6);
t63 = t103 * t74;
t27 = -t103 * t82 - t63 * t72;
t80 = qJD(3) * t103;
t46 = qJD(4) * t74 + t73 * t80;
t92 = qJD(1) * qJD(3);
t68 = t74 * t92;
t26 = qJD(2) * t46 + t68;
t76 = -t73 * qJD(4) + t74 * t80;
t75 = qJD(2) * t76 - t73 * t92;
t3 = t26 * t72 - t75 * t99;
t120 = t27 * t3;
t77 = -t72 * t73 + t81;
t119 = t3 * t77;
t118 = t3 * t56;
t117 = Ifges(4,4) * t73;
t116 = Ifges(5,4) * t49;
t50 = t77 * qJD(3);
t41 = qJD(2) * t50;
t113 = t41 * mrSges(6,2);
t93 = t73 * qJD(1);
t45 = -qJD(2) * t63 + t93;
t108 = t72 * t45;
t107 = -qJD(3) / 0.2e1;
t4 = t99 * t26 + t72 * t75;
t34 = t99 * t45;
t100 = qJD(3) * pkin(3);
t71 = t74 * qJD(1);
t86 = t103 * t73;
t44 = qJD(2) * t86 + t71;
t39 = t44 + t100;
t10 = t72 * t39 + t34;
t98 = Ifges(4,5) * qJD(3);
t97 = Ifges(4,6) * qJD(3);
t69 = -pkin(3) * t74 - pkin(2);
t96 = qJD(2) * t69;
t94 = qJD(2) * t74;
t90 = pkin(3) * t95;
t89 = mrSges(4,3) * t95;
t88 = mrSges(4,3) * t94;
t87 = m(4) * pkin(6) + mrSges(4,3);
t85 = t99 * pkin(3);
t79 = pkin(3) * t84;
t59 = pkin(6) * t94 + t93;
t9 = t39 * t99 - t108;
t60 = qJD(4) + t96;
t48 = t56 * qJD(3);
t70 = Ifges(4,4) * t94;
t67 = -t85 - pkin(4);
t65 = qJ(5) + t121;
t62 = -qJD(3) * mrSges(4,2) + t88;
t61 = qJD(3) * mrSges(4,1) - t89;
t58 = -pkin(6) * t95 + t71;
t54 = t59 * qJD(3);
t53 = -pkin(6) * t84 + t68;
t52 = Ifges(4,1) * t95 + t70 + t98;
t51 = t97 + (t74 * Ifges(4,2) + t117) * qJD(2);
t42 = Ifges(6,5) * t49;
t40 = qJD(2) * t48;
t38 = t41 * mrSges(5,2);
t37 = t40 * mrSges(6,1);
t28 = -t63 * t99 + t72 * t86;
t22 = -pkin(4) * t77 - t56 * qJ(5) + t69;
t21 = mrSges(5,1) * t47 + mrSges(5,2) * t49;
t20 = mrSges(6,1) * t47 - mrSges(6,3) * t49;
t17 = -Ifges(5,2) * t47 + Ifges(5,6) * qJD(3) + t116;
t16 = Ifges(6,6) * qJD(3) + Ifges(6,3) * t47 + t42;
t15 = t46 * t99 + t72 * t76;
t14 = t44 * t99 - t108;
t13 = t46 * t72 - t76 * t99;
t12 = t44 * t72 + t34;
t11 = pkin(4) * t49 + qJ(5) * t47 + t90;
t8 = t47 * pkin(4) - t49 * qJ(5) + t60;
t7 = qJD(3) * qJ(5) + t10;
t6 = -qJD(3) * pkin(4) + qJD(5) - t9;
t5 = pkin(4) * t48 - qJ(5) * t50 - qJD(5) * t56 + t100 * t73;
t2 = qJD(3) * qJD(5) + t4;
t1 = pkin(4) * t40 - qJ(5) * t41 - qJD(5) * t49 + t79;
t18 = [t102 * t50 - t101 * t48 + (-t73 * t61 + t74 * t62 + (-t73 ^ 2 - t74 ^ 2) * qJD(2) * mrSges(4,3)) * qJD(3) + m(4) * (t53 * t73 - t54 * t74 + (-t58 * t73 + t59 * t74) * qJD(3)) + m(5) * (t10 * t50 + t4 * t56 - t48 * t9 - t119) + m(6) * (t2 * t56 + t48 * t6 + t50 * t7 - t119) + t105 * (-t40 * t56 - t41 * t77); (t87 * t53 + (-t58 * mrSges(4,3) + t52 / 0.2e1 + t98 / 0.2e1 + (mrSges(4,2) * t128 + 0.3e1 / 0.2e1 * Ifges(4,4) * t74) * qJD(2) + (-m(4) * t58 - t61) * pkin(6)) * qJD(3)) * t74 + (t87 * t54 + (-t59 * mrSges(4,3) - t51 / 0.2e1 - t97 / 0.2e1 + (-m(4) * t59 - t62) * pkin(6) + (mrSges(4,1) * t128 - 0.3e1 / 0.2e1 * t117 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t74) * qJD(2) + (t21 + qJD(2) * (-mrSges(5,1) * t77 + mrSges(5,2) * t56) + m(5) * (t60 + t96)) * pkin(3)) * qJD(3)) * t73 + (t4 * t77 + t118) * mrSges(5,3) + (t2 * t77 + t118) * mrSges(6,2) + m(6) * (t1 * t22 + t13 * t6 + t15 * t7 + t2 * t28 + t5 * t8 + t120) + m(5) * (t10 * t15 - t13 * t9 + t28 * t4 + t120) + (-t22 * mrSges(6,3) - t104 * t77 + t105 * t27 + t132 * t56) * t41 + (mrSges(5,1) * t69 + t104 * t56 - (Ifges(5,2) + Ifges(6,3)) * t77 - t105 * t28) * t40 - t101 * t13 + t102 * t15 + t69 * t38 + t1 * (-mrSges(6,1) * t77 - mrSges(6,3) * t56) + t22 * t37 + t5 * t20 + (t60 * mrSges(5,2) + mrSges(6,2) * t6 - mrSges(5,3) * t9 - t8 * mrSges(6,3) + Ifges(5,4) * t126 + Ifges(6,5) * t125) * t50 + (-t10 * mrSges(5,3) - t7 * mrSges(6,2) + t60 * mrSges(5,1) + t8 * mrSges(6,1) + t16 / 0.2e1 - t17 / 0.2e1 + Ifges(6,3) * t125 - Ifges(5,2) * t126) * t48 + (t104 * t48 + t132 * t50) * t123 + t129 * t50 / 0.2e1 + (t130 * t48 + t131 * t50) * qJD(3) / 0.2e1; (-t11 * t8 - t12 * t6 + t2 * t65 + t3 * t67 + (-t14 + qJD(5)) * t7) * m(6) - t9 * t111 + t51 * t95 / 0.2e1 - t21 * t90 - Ifges(4,6) * t84 / 0.2e1 - t133 * t3 + (-t10 * t14 + t12 * t9 - t60 * t90 + (-t3 * t99 + t4 * t72) * pkin(3)) * m(5) - t53 * mrSges(4,2) - t54 * mrSges(4,1) - t60 * (t49 * mrSges(5,1) - t47 * mrSges(5,2)) - t8 * (t49 * mrSges(6,1) + t47 * mrSges(6,3)) + qJD(5) * t33 - (-Ifges(4,2) * t95 + t52 + t70) * t94 / 0.2e1 + (-t73 * (Ifges(4,1) * t74 - t117) / 0.2e1 + pkin(2) * (mrSges(4,1) * t73 + mrSges(4,2) * t74)) * qJD(2) ^ 2 - t11 * t20 - t4 * mrSges(5,2) + t2 * mrSges(6,3) + (-t62 + t88) * t58 + (t61 + t89) * t59 + t10 * t109 + t7 * t110 + t6 * t112 + t67 * t113 + t17 * t123 + (Ifges(6,3) * t49 - t114) * t126 + t91 * Ifges(4,5) * t74 / 0.2e1 + t101 * t12 - t102 * t14 + (-Ifges(5,2) * t49 + t129 - t43) * t125 + (-mrSges(6,2) * t65 - mrSges(5,3) * t121 + t130) * t40 + (-mrSges(5,3) * t85 + t131) * t41 + (t130 * t49 - t131 * t47) * t107 - (-t132 * t47 - t116 + t16 + t42) * t49 / 0.2e1; t40 * mrSges(5,1) - t41 * mrSges(6,3) + t101 * t49 + t102 * t47 + t37 + t38 + (t47 * t7 - t49 * t6 + t1) * m(6) + (t10 * t47 + t49 * t9 + t79) * m(5); t113 - qJD(3) * t33 + t49 * t20 + 0.2e1 * (t3 / 0.2e1 + t7 * t107 + t8 * t123) * m(6);];
tauc = t18(:);
