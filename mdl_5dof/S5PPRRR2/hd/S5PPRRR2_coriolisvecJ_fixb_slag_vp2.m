% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:11
% EndTime: 2019-12-05 15:14:15
% DurationCPUTime: 1.40s
% Computational Cost: add. (1384->209), mult. (3635->309), div. (0->0), fcn. (2564->8), ass. (0->107)
t124 = -pkin(7) - pkin(6);
t83 = sin(qJ(4));
t72 = t124 * t83;
t86 = cos(qJ(4));
t73 = t124 * t86;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t44 = t72 * t82 - t73 * t85;
t80 = sin(pkin(9));
t81 = cos(pkin(9));
t84 = sin(qJ(3));
t87 = cos(qJ(3));
t63 = t84 * t80 - t81 * t87;
t54 = t63 * qJD(1);
t66 = t82 * t86 + t83 * t85;
t99 = qJD(4) * t124;
t68 = t83 * t99;
t69 = t86 * t99;
t131 = -qJD(5) * t44 - t66 * t54 - t68 * t82 + t69 * t85;
t43 = t72 * t85 + t73 * t82;
t92 = t82 * t83 - t85 * t86;
t130 = qJD(5) * t43 - t92 * t54 + t68 * t85 + t69 * t82;
t129 = -Ifges(5,1) / 0.2e1;
t104 = qJD(3) * t86;
t128 = -Ifges(5,4) * t104 / 0.2e1;
t64 = t80 * t87 + t84 * t81;
t28 = t92 * t64;
t79 = qJD(4) + qJD(5);
t106 = qJD(2) * t83;
t55 = t64 * qJD(1);
t50 = qJD(3) * pkin(6) + t55;
t96 = pkin(7) * qJD(3) + t50;
t37 = t86 * t96 + t106;
t60 = t92 * qJD(3);
t126 = -t60 / 0.2e1;
t61 = t66 * qJD(3);
t125 = t61 / 0.2e1;
t123 = mrSges(6,3) * t60;
t122 = Ifges(5,4) * t83;
t121 = Ifges(6,4) * t61;
t56 = t63 * qJD(3);
t51 = qJD(1) * t56;
t116 = t51 * t83;
t41 = t50 * t86 + t106;
t19 = -qJD(4) * t41 + t116;
t120 = t19 * t83;
t119 = t37 * t82;
t118 = t37 * t85;
t78 = t86 * qJD(2);
t40 = -t50 * t83 + t78;
t117 = t40 * t83;
t57 = t64 * qJD(3);
t52 = qJD(1) * t57;
t115 = t52 * t63;
t114 = t61 * mrSges(6,3);
t38 = t79 * t92;
t32 = t38 * qJD(3);
t113 = t92 * t32;
t39 = t79 * t66;
t33 = t39 * qJD(3);
t112 = t66 * t33;
t109 = qJD(4) * t78 - t86 * t51;
t108 = Ifges(5,5) * qJD(4);
t107 = Ifges(5,6) * qJD(4);
t105 = qJD(3) * t83;
t103 = qJD(4) * t83;
t102 = qJD(5) * t82;
t101 = qJD(5) * t85;
t100 = pkin(4) * t103;
t76 = -pkin(4) * t86 - pkin(3);
t98 = t108 / 0.2e1;
t97 = -t107 / 0.2e1;
t35 = mrSges(6,1) * t60 + mrSges(6,2) * t61;
t95 = -t35 + (mrSges(5,1) * t86 - mrSges(5,2) * t83 + mrSges(4,1)) * qJD(3);
t91 = t96 * t83;
t36 = t78 - t91;
t31 = qJD(4) * pkin(4) + t36;
t7 = t31 * t85 - t119;
t8 = t31 * t82 + t118;
t94 = t41 * t86 - t117;
t70 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t105;
t71 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t104;
t93 = -t83 * t70 + t86 * t71;
t49 = -qJD(3) * pkin(3) + t54;
t90 = t40 * mrSges(5,3) + t105 * t129 - t108 / 0.2e1 + t128 - t49 * mrSges(5,2);
t89 = t41 * mrSges(5,3) + t107 / 0.2e1 + (Ifges(5,2) * t86 + t122) * qJD(3) / 0.2e1 - t49 * mrSges(5,1);
t12 = -qJD(4) * t91 + t109;
t13 = -qJD(4) * t37 + t116;
t2 = qJD(5) * t7 + t12 * t85 + t13 * t82;
t25 = -Ifges(6,2) * t60 + Ifges(6,6) * t79 + t121;
t53 = Ifges(6,4) * t60;
t26 = Ifges(6,1) * t61 + Ifges(6,5) * t79 - t53;
t3 = -qJD(5) * t8 - t12 * t82 + t13 * t85;
t45 = qJD(3) * t76 + t54;
t88 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t7 * t123 + t25 * t125 - t45 * (mrSges(6,1) * t61 - mrSges(6,2) * t60) - t61 * (-Ifges(6,1) * t60 - t121) / 0.2e1 - t79 * (-Ifges(6,5) * t60 - Ifges(6,6) * t61) / 0.2e1 - Ifges(6,6) * t33 - Ifges(6,5) * t32 + (-Ifges(6,2) * t61 + t26 - t53) * t60 / 0.2e1;
t62 = (mrSges(5,1) * t83 + mrSges(5,2) * t86) * qJD(4) * qJD(3);
t47 = mrSges(6,1) * t79 - t114;
t46 = -mrSges(6,2) * t79 - t123;
t42 = (t55 + t100) * qJD(3);
t27 = t66 * t64;
t18 = -t103 * t50 + t109;
t11 = t36 * t85 - t119;
t10 = -t36 * t82 - t118;
t9 = mrSges(6,1) * t33 - mrSges(6,2) * t32;
t5 = t28 * t79 + t66 * t56;
t4 = -t39 * t64 + t56 * t92;
t1 = [t4 * t46 + t5 * t47 + (t62 + t9) * t63 + (-t70 * t86 - t71 * t83) * t64 * qJD(4) + (-t27 * t32 + t28 * t33) * mrSges(6,3) - t95 * t57 - (-qJD(3) * mrSges(4,2) + t93) * t56 + m(4) * (-t51 * t64 + t54 * t57 - t55 * t56 + t115) + m(5) * (t49 * t57 + t115 - t94 * t56 + (t18 * t86 - t120 + (-t40 * t86 - t41 * t83) * qJD(4)) * t64) + m(6) * (-t2 * t28 - t27 * t3 + t4 * t8 + t42 * t63 + t45 * t57 + t5 * t7); -t38 * t46 - t39 * t47 + (-t112 - t113) * mrSges(6,3) + m(5) * (t18 * t83 + t19 * t86) + m(6) * (t2 * t66 - t3 * t92 - t38 * t8 - t39 * t7) + (m(5) * t94 + (-t83 ^ 2 - t86 ^ 2) * qJD(3) * mrSges(5,3) + t93) * qJD(4); t76 * t9 + t79 * (-Ifges(6,5) * t38 - Ifges(6,6) * t39) / 0.2e1 + t42 * (mrSges(6,1) * t92 + t66 * mrSges(6,2)) - pkin(3) * t62 - t52 * mrSges(4,1) - t39 * t25 / 0.2e1 + t45 * (t39 * mrSges(6,1) - t38 * mrSges(6,2)) - t38 * t26 / 0.2e1 + t131 * t47 + t130 * t46 + (-t126 * t39 + t33 * t92) * Ifges(6,2) + (-t125 * t38 - t32 * t66) * Ifges(6,1) + (-qJD(3) * t54 + t51) * mrSges(4,2) + t95 * t55 + (t52 * mrSges(5,2) - t19 * mrSges(5,3) - t54 * t70 + (-pkin(6) * t71 + pkin(4) * t35 - 0.3e1 / 0.2e1 * Ifges(5,4) * t105 + t97 - t89) * qJD(4)) * t83 - m(5) * (t54 * t117 + t49 * t55) + m(5) * (-pkin(3) * t52 + (-t41 * t103 - t120) * pkin(6)) + (-t2 * t92 - t3 * t66 + t32 * t43 - t33 * t44 + t38 * t7 - t39 * t8) * mrSges(6,3) + (-t125 * t39 - t126 * t38 - t112 + t113) * Ifges(6,4) + (-t52 * mrSges(5,1) + (t98 + (-m(5) * t40 - t70) * pkin(6) + (0.3e1 / 0.2e1 * Ifges(5,4) * t86 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t83) * qJD(3) - t90) * qJD(4) - (-m(5) * t41 - t71) * t54 + (m(5) * pkin(6) + mrSges(5,3)) * t18) * t86 + (t2 * t44 + t3 * t43 + t42 * t76 + t130 * t8 + t131 * t7 + (t100 - t55) * t45) * m(6); t88 + t8 * t114 + t41 * t70 - t40 * t71 - t11 * t46 - t10 * t47 - t18 * mrSges(5,2) + t19 * mrSges(5,1) - m(6) * (t10 * t7 + t11 * t8) + (m(6) * (t101 * t8 - t102 * t7 + t2 * t82 + t3 * t85) + t46 * t101 - t47 * t102 + (t32 * t85 - t33 * t82) * mrSges(6,3)) * pkin(4) + ((t128 + t98 + t90) * t86 + (t97 + (t122 / 0.2e1 + (Ifges(5,2) / 0.2e1 + t129) * t86) * qJD(3) + (-m(6) * t45 - t35) * pkin(4) + t89) * t83) * qJD(3); t88 - t46 * t7 + (t47 + t114) * t8;];
tauc = t1(:);
