% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
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
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:08:39
% EndTime: 2019-12-31 18:08:43
% DurationCPUTime: 2.11s
% Computational Cost: add. (1185->240), mult. (3054->322), div. (0->0), fcn. (1797->6), ass. (0->116)
t140 = mrSges(6,1) + mrSges(5,1);
t139 = Ifges(5,1) + Ifges(6,1);
t138 = Ifges(6,4) + Ifges(5,5);
t78 = sin(qJ(3));
t97 = qJD(1) * qJD(3);
t90 = t78 * t97;
t112 = mrSges(6,2) + mrSges(5,3);
t111 = -Ifges(5,4) + Ifges(6,5);
t137 = Ifges(6,6) - Ifges(5,6);
t103 = qJD(1) * t78;
t75 = sin(pkin(8));
t107 = cos(pkin(8));
t79 = cos(qJ(3));
t87 = t107 * t79;
t53 = -qJD(1) * t87 + t103 * t75;
t122 = Ifges(6,5) * t53;
t51 = Ifges(5,4) * t53;
t88 = t107 * t78;
t61 = t75 * t79 + t88;
t55 = t61 * qJD(1);
t136 = t138 * qJD(3) + t139 * t55 + t122 - t51;
t117 = t55 * mrSges(5,3);
t118 = t55 * mrSges(6,2);
t109 = qJD(3) * t140 - t117 - t118;
t119 = t53 * mrSges(5,3);
t120 = t53 * mrSges(6,2);
t40 = qJD(3) * mrSges(6,3) - t120;
t135 = qJD(3) * mrSges(5,2) + t119 - t40;
t70 = sin(pkin(7)) * pkin(1) + pkin(6);
t63 = t70 * qJD(1);
t85 = qJ(4) * qJD(1) + t63;
t99 = t78 * qJD(2);
t35 = t79 * t85 + t99;
t134 = -t53 / 0.2e1;
t133 = t53 / 0.2e1;
t131 = t55 / 0.2e1;
t129 = pkin(3) * t75;
t100 = qJD(4) * t79;
t101 = qJD(3) * t78;
t74 = t79 * qJD(2);
t41 = qJD(3) * t74 - t101 * t63;
t28 = (-qJ(4) * t101 + t100) * qJD(1) + t41;
t98 = t78 * qJD(4);
t80 = -qJD(1) * t98 - qJD(3) * t35;
t2 = -t107 * t80 + t28 * t75;
t108 = qJ(4) + t70;
t59 = t108 * t79;
t24 = t108 * t88 + t59 * t75;
t128 = t2 * t24;
t115 = t75 * t78;
t82 = t87 - t115;
t127 = t2 * t82;
t126 = t2 * t61;
t125 = Ifges(4,4) * t78;
t124 = Ifges(5,4) * t55;
t56 = t82 * qJD(3);
t46 = qJD(1) * t56;
t121 = t46 * mrSges(6,2);
t116 = t75 * t35;
t114 = -qJD(3) / 0.2e1;
t3 = t107 * t28 + t75 * t80;
t30 = t107 * t35;
t34 = -t78 * t85 + t74;
t32 = qJD(3) * pkin(3) + t34;
t8 = t75 * t32 + t30;
t106 = Ifges(4,5) * qJD(3);
t105 = Ifges(4,6) * qJD(3);
t93 = -cos(pkin(7)) * pkin(1) - pkin(2);
t62 = -pkin(3) * t79 + t93;
t104 = qJD(1) * t62;
t102 = qJD(1) * t79;
t96 = pkin(3) * t103;
t95 = mrSges(4,3) * t103;
t94 = mrSges(4,3) * t102;
t73 = Ifges(4,4) * t102;
t92 = m(4) * t70 + mrSges(4,3);
t91 = t107 * pkin(3);
t86 = qJD(3) * t108;
t84 = pkin(3) * t90;
t65 = t93 * qJD(1);
t48 = t63 * t79 + t99;
t7 = t107 * t32 - t116;
t81 = -t79 * t86 - t98;
t54 = t61 * qJD(3);
t52 = qJD(4) + t104;
t71 = -t91 - pkin(4);
t68 = qJ(5) + t129;
t66 = -qJD(3) * mrSges(4,2) + t94;
t64 = qJD(3) * mrSges(4,1) - t95;
t58 = Ifges(4,1) * t103 + t106 + t73;
t57 = t105 + (t79 * Ifges(4,2) + t125) * qJD(1);
t50 = Ifges(6,5) * t55;
t47 = -t63 * t78 + t74;
t45 = qJD(1) * t54;
t44 = t46 * mrSges(5,2);
t43 = t45 * mrSges(6,1);
t42 = t48 * qJD(3);
t36 = -t78 * t86 + t100;
t27 = mrSges(5,1) * t53 + mrSges(5,2) * t55;
t26 = mrSges(6,1) * t53 - mrSges(6,3) * t55;
t25 = t107 * t59 - t108 * t115;
t18 = -Ifges(5,2) * t53 + Ifges(5,6) * qJD(3) + t124;
t17 = Ifges(6,6) * qJD(3) + Ifges(6,3) * t53 + t50;
t16 = -pkin(4) * t82 - t61 * qJ(5) + t62;
t15 = pkin(4) * t55 + qJ(5) * t53 + t96;
t14 = t53 * pkin(4) - t55 * qJ(5) + t52;
t13 = t107 * t36 + t75 * t81;
t12 = -t107 * t81 + t36 * t75;
t11 = t107 * t34 - t116;
t10 = t34 * t75 + t30;
t9 = pkin(3) * t101 + pkin(4) * t54 - qJ(5) * t56 - qJD(5) * t61;
t6 = qJD(3) * qJ(5) + t8;
t5 = -qJD(3) * pkin(4) + qJD(5) - t7;
t4 = pkin(4) * t45 - qJ(5) * t46 - qJD(5) * t55 + t84;
t1 = qJD(3) * qJD(5) + t3;
t19 = [(t92 * t42 + (t65 * mrSges(4,1) - t57 / 0.2e1 - t70 * t66 - t105 / 0.2e1 - t92 * t48 + (t93 * mrSges(4,1) - 0.3e1 / 0.2e1 * t125 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t79) * qJD(1) + (qJD(1) * (-mrSges(5,1) * t82 + mrSges(5,2) * t61) + m(5) * (t52 + t104) + t27) * pkin(3)) * qJD(3)) * t78 + (t3 * t82 + t126) * mrSges(5,3) + (t1 * t82 + t126) * mrSges(6,2) + m(6) * (t1 * t25 + t12 * t5 + t13 * t6 + t14 * t9 + t16 * t4 + t128) + m(5) * (-t12 * t7 + t13 * t8 + t25 * t3 + t128) + t16 * t43 + t4 * (-mrSges(6,1) * t82 - mrSges(6,3) * t61) + (-t16 * mrSges(6,3) - t111 * t82 + t112 * t24 + t139 * t61) * t46 + (mrSges(5,1) * t62 + t111 * t61 - (Ifges(5,2) + Ifges(6,3)) * t82 - t112 * t25) * t45 + (t92 * t41 + (t58 / 0.2e1 - t70 * t64 + 0.3e1 / 0.2e1 * t73 + t106 / 0.2e1 - t92 * t47 + 0.2e1 * t65 * mrSges(4,2)) * qJD(3)) * t79 - t109 * t12 - t135 * t13 + t9 * t26 + t62 * t44 + (t52 * mrSges(5,2) + t5 * mrSges(6,2) - t7 * mrSges(5,3) - t14 * mrSges(6,3) + Ifges(5,4) * t134 + Ifges(6,5) * t133) * t56 + (-t8 * mrSges(5,3) - t6 * mrSges(6,2) + t14 * mrSges(6,1) + t52 * mrSges(5,1) + Ifges(6,3) * t133 - Ifges(5,2) * t134 + t17 / 0.2e1 - t18 / 0.2e1) * t54 + (t111 * t54 + t139 * t56) * t131 + t136 * t56 / 0.2e1 + (t137 * t54 + t138 * t56) * qJD(3) / 0.2e1; -t135 * t56 - t109 * t54 + (-t78 * t64 + t79 * t66 + (-t78 ^ 2 - t79 ^ 2) * qJD(1) * mrSges(4,3)) * qJD(3) + m(4) * (t41 * t78 - t42 * t79 + (-t47 * t78 + t48 * t79) * qJD(3)) + m(5) * (t3 * t61 - t54 * t7 + t56 * t8 - t127) + m(6) * (t1 * t61 + t5 * t54 + t56 * t6 - t127) + t112 * (-t45 * t61 - t46 * t82); -t140 * t2 + t57 * t103 / 0.2e1 - t27 * t96 + t109 * t10 + (-Ifges(5,2) * t55 + t136 - t51) * t133 + (-mrSges(6,2) * t68 - mrSges(5,3) * t129 + t137) * t45 + (-mrSges(5,3) * t91 + t138) * t46 + (t137 * t55 - t138 * t53) * t114 - (-t139 * t53 - t124 + t17 + t50) * t55 / 0.2e1 + (t10 * t7 - t11 * t8 - t52 * t96 + (-t107 * t2 + t3 * t75) * pkin(3)) * m(5) - (-Ifges(4,2) * t103 + t58 + t73) * t102 / 0.2e1 - t7 * t119 + t18 * t131 + (Ifges(6,3) * t55 - t122) * t134 + t8 * t117 + t6 * t118 + t5 * t120 + t71 * t121 + (t1 * t68 - t10 * t5 - t14 * t15 + t2 * t71 + (-t11 + qJD(5)) * t6) * m(6) - Ifges(4,6) * t90 / 0.2e1 - t14 * (t55 * mrSges(6,1) + t53 * mrSges(6,3)) - t52 * (t55 * mrSges(5,1) - t53 * mrSges(5,2)) + qJD(5) * t40 - t41 * mrSges(4,2) - t42 * mrSges(4,1) - t15 * t26 + t1 * mrSges(6,3) - t3 * mrSges(5,2) + (-t65 * (mrSges(4,1) * t78 + mrSges(4,2) * t79) - (Ifges(4,1) * t79 - t125) * t103 / 0.2e1) * qJD(1) + t135 * t11 + (-t66 + t94) * t47 + (t64 + t95) * t48 + t97 * Ifges(4,5) * t79 / 0.2e1; t45 * mrSges(5,1) - t46 * mrSges(6,3) + t109 * t55 - t135 * t53 + t43 + t44 + (-t5 * t55 + t53 * t6 + t4) * m(6) + (t53 * t8 + t55 * t7 + t84) * m(5); t121 - qJD(3) * t40 + t55 * t26 + 0.2e1 * (t2 / 0.2e1 + t6 * t114 + t14 * t131) * m(6);];
tauc = t19(:);
