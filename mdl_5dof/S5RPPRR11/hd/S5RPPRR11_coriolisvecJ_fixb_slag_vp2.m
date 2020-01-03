% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:32
% EndTime: 2019-12-31 18:05:36
% DurationCPUTime: 1.49s
% Computational Cost: add. (1141->233), mult. (2425->332), div. (0->0), fcn. (1133->4), ass. (0->112)
t56 = cos(qJ(5));
t54 = sin(qJ(5));
t94 = qJD(4) * t54;
t57 = cos(qJ(4));
t99 = qJD(1) * t57;
t37 = t56 * t99 + t94;
t115 = Ifges(6,4) * t37;
t92 = qJD(4) * t56;
t36 = -t54 * t99 + t92;
t55 = sin(qJ(4));
t100 = qJD(1) * t55;
t46 = qJD(5) + t100;
t11 = Ifges(6,2) * t36 + Ifges(6,6) * t46 + t115;
t111 = Ifges(6,6) * t54;
t112 = Ifges(6,5) * t56;
t119 = t56 / 0.2e1;
t34 = Ifges(6,4) * t36;
t12 = Ifges(6,1) * t37 + Ifges(6,5) * t46 + t34;
t120 = -t54 / 0.2e1;
t123 = t37 / 0.2e1;
t48 = qJD(1) * qJ(2) + qJD(3);
t43 = -qJD(1) * pkin(6) + t48;
t106 = t43 * t57;
t33 = -qJD(4) * pkin(4) - t106;
t113 = Ifges(6,4) * t56;
t70 = -Ifges(6,2) * t54 + t113;
t114 = Ifges(6,4) * t54;
t72 = Ifges(6,1) * t56 - t114;
t73 = mrSges(6,1) * t54 + mrSges(6,2) * t56;
t53 = pkin(1) + qJ(3);
t40 = pkin(4) * t55 - pkin(7) * t57 + t53;
t28 = t40 * qJD(1) - qJD(2);
t107 = t43 * t55;
t32 = qJD(4) * pkin(7) + t107;
t8 = t28 * t56 - t32 * t54;
t9 = t28 * t54 + t32 * t56;
t76 = t54 * t9 + t56 * t8;
t136 = -t76 * mrSges(6,3) + t11 * t120 + t12 * t119 + (-t111 + t112) * t46 / 0.2e1 + t72 * t123 + t70 * t36 / 0.2e1 + t33 * t73;
t102 = Ifges(5,5) * qJD(4);
t116 = Ifges(5,4) * t55;
t134 = qJD(1) / 0.2e1;
t135 = t102 / 0.2e1 + (t57 * Ifges(5,1) - t116) * t134 + t136;
t133 = t53 * qJD(1);
t87 = qJD(1) * qJD(2);
t91 = qJD(4) * t57;
t27 = t43 * t91 + t55 * t87;
t80 = pkin(4) * t57 + pkin(7) * t55;
t35 = qJD(4) * t80 + qJD(3);
t29 = t35 * qJD(1);
t1 = qJD(5) * t8 + t27 * t56 + t29 * t54;
t2 = -qJD(5) * t9 - t27 * t54 + t29 * t56;
t78 = t1 * t56 - t2 * t54;
t85 = qJD(4) * qJD(5);
t88 = qJD(5) * t57;
t20 = t56 * t85 + (-t54 * t88 - t55 * t92) * qJD(1);
t93 = qJD(4) * t55;
t21 = -t54 * t85 + (t54 * t93 - t56 * t88) * qJD(1);
t132 = t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t20 + Ifges(6,6) * t21;
t24 = -mrSges(6,2) * t46 + mrSges(6,3) * t36;
t25 = mrSges(6,1) * t46 - mrSges(6,3) * t37;
t65 = -t54 * t24 - t56 * t25;
t131 = -m(6) * t76 + t65;
t130 = (t55 ^ 2 + t57 ^ 2) * t43;
t128 = t20 / 0.2e1;
t127 = t21 / 0.2e1;
t126 = -t36 / 0.2e1;
t124 = -t37 / 0.2e1;
t122 = -t46 / 0.2e1;
t26 = t43 * t93 - t57 * t87;
t52 = qJ(2) - pkin(6);
t110 = t26 * t52;
t104 = t52 * t55;
t96 = qJD(4) * mrSges(5,1);
t103 = -mrSges(6,1) * t36 + mrSges(6,2) * t37 + mrSges(5,3) * t99 - t96;
t101 = Ifges(5,6) * qJD(4);
t44 = -qJD(2) + t133;
t97 = qJD(3) * t44;
t95 = qJD(4) * mrSges(5,2);
t90 = qJD(5) * t54;
t89 = qJD(5) * t56;
t86 = qJD(1) * qJD(3);
t84 = qJD(1) * t91;
t83 = m(6) * t33 + t103;
t82 = m(3) * qJ(2) + mrSges(4,2) + mrSges(3,3);
t79 = -qJD(5) * t104 + t35;
t77 = -t1 * t54 - t2 * t56;
t75 = t8 * t54 - t9 * t56;
t74 = t56 * mrSges(6,1) - t54 * mrSges(6,2);
t71 = Ifges(6,1) * t54 + t113;
t69 = Ifges(6,2) * t56 + t114;
t67 = Ifges(6,5) * t54 + Ifges(6,6) * t56;
t13 = mrSges(6,1) * t84 - mrSges(6,3) * t20;
t14 = -mrSges(6,2) * t84 + mrSges(6,3) * t21;
t66 = -t54 * t13 + t56 * t14;
t41 = -mrSges(5,3) * t100 - t95;
t64 = t56 * t24 - t54 * t25 + t41;
t61 = qJD(2) * t55 + qJD(5) * t40 + t52 * t91;
t60 = t9 * mrSges(6,2) - t46 * Ifges(6,3) - t37 * Ifges(6,5) - t36 * Ifges(6,6) + t101 / 0.2e1 + (Ifges(5,4) * t57 - t55 * Ifges(5,2)) * t134 - t44 * mrSges(5,1) - t8 * mrSges(6,1);
t58 = qJD(1) ^ 2;
t45 = Ifges(6,3) * t84;
t39 = t80 * qJD(1);
t38 = qJD(1) * (mrSges(5,1) * t55 + mrSges(5,2) * t57);
t23 = t56 * t104 + t40 * t54;
t22 = -t54 * t104 + t40 * t56;
t16 = t56 * t106 + t39 * t54;
t15 = -t54 * t106 + t39 * t56;
t7 = -mrSges(6,1) * t21 + mrSges(6,2) * t20;
t6 = Ifges(6,1) * t20 + Ifges(6,4) * t21 + Ifges(6,5) * t84;
t5 = Ifges(6,4) * t20 + Ifges(6,2) * t21 + Ifges(6,6) * t84;
t4 = -t54 * t61 + t56 * t79;
t3 = t54 * t79 + t56 * t61;
t10 = [qJD(3) * t38 + t22 * t13 + t23 * t14 + t3 * t24 + t4 * t25 + 0.2e1 * (qJD(3) * mrSges(4,3) + qJD(2) * t82) * qJD(1) + m(6) * (t1 * t23 + t2 * t22 + t3 * t9 + t4 * t8) + m(5) * (qJD(2) * t130 + t53 * t86 + t97) + m(4) * (qJD(2) * t48 + t97 + (qJ(2) * qJD(2) + qJD(3) * t53) * qJD(1)) + (qJD(2) * t41 + mrSges(5,1) * t86 + t45 / 0.2e1 + (m(5) * t52 - mrSges(5,3)) * t27 + (0.3e1 / 0.2e1 * Ifges(5,4) * t100 - t102 / 0.2e1 + (-t44 - t133) * mrSges(5,2) + t83 * t52 - t135) * qJD(4) + t132) * t55 + (t72 * t128 + t70 * t127 - t52 * t7 + t6 * t119 + t5 * t120 + mrSges(5,2) * t86 + (mrSges(5,3) + t73) * t26 - t103 * qJD(2) + t77 * mrSges(6,3) + m(6) * (-qJD(2) * t33 - t110) - m(5) * t110 + (t67 * t122 + t69 * t126 + t71 * t124 + t33 * t74 - t56 * t11 / 0.2e1 + t12 * t120 + t75 * mrSges(6,3)) * qJD(5) + (t52 * t41 - t101 / 0.2e1 + ((-0.3e1 / 0.2e1 * Ifges(5,4) + t112 / 0.2e1 - t111 / 0.2e1) * t57 + t53 * mrSges(5,1) + (Ifges(6,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(5,2) - 0.3e1 / 0.2e1 * Ifges(5,1)) * t55) * qJD(1) - t60) * qJD(4)) * t57; -t24 * t89 + t25 * t90 - t54 * t14 - t56 * t13 + m(6) * (qJD(5) * t75 + t77) - t82 * t58 + ((-qJD(3) - t130) * m(5) + (-qJD(3) - t48) * m(4) + (t83 - t96) * t57 + (m(6) * t75 - t64 + t95) * t55) * qJD(1); -t58 * mrSges(4,3) + (-t7 + t64 * qJD(4) - m(5) * t26 + m(6) * (-t8 * t94 + t9 * t92 - t26)) * t57 + (t65 * qJD(5) + t103 * qJD(4) + m(5) * t27 + m(6) * (qJD(4) * t33 - t8 * t89 - t9 * t90 + t78) + t66) * t55 + (-m(5) * t44 - t38 + (qJD(2) - t44) * m(4) + t131) * qJD(1); t71 * t128 + t69 * t127 + t5 * t119 + t54 * t6 / 0.2e1 - t16 * t24 - t15 * t25 - t27 * mrSges(5,2) - pkin(4) * t7 + (-mrSges(5,1) - t74) * t26 + t66 * pkin(7) + t78 * mrSges(6,3) + (-t103 * t55 - t57 * t41) * t43 + (pkin(7) * t131 + t136) * qJD(5) + ((Ifges(5,4) * t99 / 0.2e1 + t60) * t57 + (t44 * mrSges(5,2) + (-t116 / 0.2e1 + (Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t57) * qJD(1) + t135) * t55 + (-Ifges(5,5) * t55 / 0.2e1 + (-Ifges(5,6) / 0.2e1 + t67 / 0.2e1) * t57) * qJD(4)) * qJD(1) + (-pkin(4) * t26 + pkin(7) * t78 - t33 * t107 - t15 * t8 - t16 * t9) * m(6); t45 - t33 * (mrSges(6,1) * t37 + mrSges(6,2) * t36) + (Ifges(6,1) * t36 - t115) * t124 + t11 * t123 + (Ifges(6,5) * t36 - Ifges(6,6) * t37) * t122 - t8 * t24 + t9 * t25 + (t36 * t8 + t37 * t9) * mrSges(6,3) + (-Ifges(6,2) * t37 + t12 + t34) * t126 + t132;];
tauc = t10(:);
