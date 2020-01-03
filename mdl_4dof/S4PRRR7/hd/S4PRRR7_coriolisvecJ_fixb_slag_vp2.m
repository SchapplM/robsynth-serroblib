% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:04
% EndTime: 2019-12-31 16:36:08
% DurationCPUTime: 1.57s
% Computational Cost: add. (997->229), mult. (2794->352), div. (0->0), fcn. (1769->8), ass. (0->120)
t146 = qJD(3) / 0.2e1;
t65 = sin(qJ(3));
t103 = qJD(2) * t65;
t145 = t103 / 0.2e1;
t87 = Ifges(4,5) * t146;
t63 = cos(pkin(4));
t104 = qJD(1) * t63;
t62 = sin(pkin(4));
t105 = qJD(1) * t62;
t66 = sin(qJ(2));
t92 = t66 * t105;
t53 = qJD(2) * pkin(6) + t92;
t68 = cos(qJ(3));
t37 = t68 * t104 - t53 * t65;
t108 = qJD(2) * pkin(2);
t69 = cos(qJ(2));
t91 = t69 * t105;
t54 = -t91 - t108;
t101 = qJD(2) * t68;
t61 = Ifges(4,4) * t101;
t64 = sin(qJ(4));
t117 = Ifges(5,6) * t64;
t67 = cos(qJ(4));
t119 = Ifges(5,5) * t67;
t129 = t67 / 0.2e1;
t130 = -t64 / 0.2e1;
t49 = qJD(3) * t64 + t67 * t103;
t132 = t49 / 0.2e1;
t123 = Ifges(5,4) * t49;
t95 = t67 * qJD(3);
t48 = -t64 * t103 + t95;
t59 = qJD(4) - t101;
t15 = Ifges(5,2) * t48 + Ifges(5,6) * t59 + t123;
t47 = Ifges(5,4) * t48;
t16 = Ifges(5,1) * t49 + Ifges(5,5) * t59 + t47;
t26 = -qJD(3) * pkin(3) - t37;
t121 = Ifges(5,4) * t67;
t78 = -Ifges(5,2) * t64 + t121;
t122 = Ifges(5,4) * t64;
t80 = Ifges(5,1) * t67 - t122;
t81 = mrSges(5,1) * t64 + mrSges(5,2) * t67;
t38 = t65 * t104 + t53 * t68;
t27 = qJD(3) * pkin(7) + t38;
t55 = -pkin(3) * t68 - pkin(7) * t65 - pkin(2);
t39 = t55 * qJD(2) - t91;
t7 = -t27 * t64 + t39 * t67;
t8 = t27 * t67 + t39 * t64;
t83 = t64 * t8 + t67 * t7;
t71 = -t83 * mrSges(5,3) + t59 * (-t117 + t119) / 0.2e1 + t48 * t78 / 0.2e1 + t80 * t132 + t26 * t81 + t16 * t129 + t15 * t130;
t144 = t54 * mrSges(4,2) - t37 * mrSges(4,3) + Ifges(4,1) * t145 + t87 + t61 / 0.2e1 + t71;
t143 = m(4) * pkin(6);
t86 = -Ifges(4,6) * qJD(3) / 0.2e1;
t100 = qJD(3) * t65;
t110 = t68 * t69;
t84 = pkin(3) * t65 - pkin(7) * t68;
t52 = t84 * qJD(3);
t96 = qJD(4) * t68;
t98 = qJD(4) * t64;
t142 = -(-t64 * t110 + t66 * t67) * t105 - t55 * t98 + t52 * t67 + (t64 * t100 - t67 * t96) * pkin(6);
t97 = qJD(4) * t67;
t141 = -(t67 * t110 + t64 * t66) * t105 + t55 * t97 + t52 * t64 + (-t64 * t96 - t65 * t95) * pkin(6);
t111 = t62 * t69;
t89 = qJD(2) * t111;
t74 = qJD(1) * (qJD(3) * t63 + t89);
t17 = -t53 * t100 + t68 * t74;
t34 = (t52 + t92) * qJD(2);
t1 = t7 * qJD(4) + t17 * t67 + t34 * t64;
t2 = -t8 * qJD(4) - t17 * t64 + t34 * t67;
t140 = t1 * t67 - t2 * t64;
t93 = qJD(3) * qJD(4);
t30 = t67 * t93 + (-t65 * t98 + t68 * t95) * qJD(2);
t99 = qJD(3) * t68;
t31 = -t64 * t93 + (-t64 * t99 - t65 * t97) * qJD(2);
t139 = -t2 * mrSges(5,1) + t1 * mrSges(5,2) - Ifges(5,5) * t30 - Ifges(5,6) * t31;
t109 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t48 - mrSges(5,2) * t49 - mrSges(4,3) * t103;
t138 = m(4) * t37 - m(5) * t26 + t109;
t136 = t30 / 0.2e1;
t135 = t31 / 0.2e1;
t134 = -t48 / 0.2e1;
t133 = -t49 / 0.2e1;
t131 = -t59 / 0.2e1;
t126 = pkin(6) * t68;
t18 = t53 * t99 + t65 * t74;
t112 = t62 * t66;
t42 = t65 * t112 - t63 * t68;
t115 = t18 * t42;
t102 = qJD(2) * t66;
t94 = qJD(2) * qJD(3);
t90 = t62 * t102;
t88 = t65 * t94;
t82 = t67 * mrSges(5,1) - t64 * mrSges(5,2);
t79 = Ifges(5,1) * t64 + t121;
t77 = Ifges(5,2) * t67 + t122;
t76 = Ifges(5,5) * t64 + Ifges(5,6) * t67;
t43 = t68 * t112 + t63 * t65;
t24 = -t67 * t111 - t43 * t64;
t75 = t64 * t111 - t43 * t67;
t73 = t54 * mrSges(4,1) + t7 * mrSges(5,1) + Ifges(5,3) * t59 + Ifges(5,6) * t48 + Ifges(5,5) * t49 + t86 - (Ifges(4,4) * t65 + Ifges(4,2) * t68) * qJD(2) / 0.2e1 - t38 * mrSges(4,3) - t8 * mrSges(5,2);
t70 = qJD(2) ^ 2;
t58 = Ifges(5,3) * t88;
t57 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t101;
t51 = t84 * qJD(2);
t50 = (-mrSges(4,1) * t68 + mrSges(4,2) * t65) * qJD(2);
t46 = (mrSges(4,1) * t65 + mrSges(4,2) * t68) * t94;
t41 = t67 * t126 + t55 * t64;
t40 = -t64 * t126 + t55 * t67;
t36 = mrSges(5,1) * t59 - mrSges(5,3) * t49;
t35 = -mrSges(5,2) * t59 + mrSges(5,3) * t48;
t23 = -t42 * qJD(3) + t68 * t89;
t22 = t43 * qJD(3) + t65 * t89;
t20 = -mrSges(5,2) * t88 + mrSges(5,3) * t31;
t19 = mrSges(5,1) * t88 - mrSges(5,3) * t30;
t13 = t37 * t67 + t51 * t64;
t12 = -t37 * t64 + t51 * t67;
t9 = -mrSges(5,1) * t31 + mrSges(5,2) * t30;
t6 = Ifges(5,1) * t30 + Ifges(5,4) * t31 + Ifges(5,5) * t88;
t5 = Ifges(5,4) * t30 + Ifges(5,2) * t31 + Ifges(5,6) * t88;
t4 = t24 * qJD(4) + t23 * t67 + t64 * t90;
t3 = t75 * qJD(4) - t23 * t64 + t67 * t90;
t10 = [t24 * t19 - t75 * t20 + t23 * t57 + t3 * t36 + t4 * t35 + t42 * t9 - t109 * t22 + (t42 * t68 - t43 * t65) * mrSges(4,3) * t94 + ((-mrSges(3,2) * t70 - t46) * t69 + (-mrSges(3,1) * t70 + qJD(2) * t50) * t66) * t62 + m(4) * (t17 * t43 + t115 - t22 * t37 + t23 * t38 + (t54 - t91) * t90) + m(5) * (-t1 * t75 + t2 * t24 + t22 * t26 + t3 * t7 + t4 * t8 + t115); -t50 * t92 - pkin(2) * t46 + t40 * t19 + t41 * t20 + t142 * t36 + t141 * t35 + 0.2e1 * (-t108 / 0.2e1 - t54 / 0.2e1) * m(4) * t92 + (-t58 / 0.2e1 + (mrSges(4,3) + t143) * t17 + (-mrSges(4,1) * t102 + (-m(4) * t38 - t57) * t69) * t105 + (0.3e1 / 0.2e1 * t61 + t87 - t138 * pkin(6) + t144) * qJD(3) + t139) * t68 + (t78 * t135 + t6 * t129 + t5 * t130 + pkin(6) * t9 + t80 * t136 + (mrSges(4,3) + t81) * t18 + (-t1 * t64 - t2 * t67) * mrSges(5,3) + (-pkin(6) * t57 + t73 + t86) * qJD(3) + (t76 * t131 + t77 * t134 + t79 * t133 + t26 * t82 - t67 * t15 / 0.2e1 + t16 * t130 + (t7 * t64 - t8 * t67) * mrSges(5,3)) * qJD(4) + (-t38 * qJD(3) + t18) * t143 + (mrSges(4,2) * t92 + ((t119 / 0.2e1 - t117 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t65 + (0.3e1 / 0.2e1 * Ifges(4,1) - Ifges(5,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,2)) * t68) * qJD(3)) * qJD(2)) * t65 + t138 * t65 * t91 + (pkin(6) * t18 * t65 + t1 * t41 + t141 * t8 + t142 * t7 + t2 * t40) * m(5); t79 * t136 + t77 * t135 + t5 * t129 + t64 * t6 / 0.2e1 - t37 * t57 - t12 * t36 - t13 * t35 - t17 * mrSges(4,2) - pkin(3) * t9 + t109 * t38 + (-mrSges(4,1) - t82) * t18 + t140 * mrSges(5,3) + t71 * qJD(4) + ((Ifges(4,4) * t145 + t76 * t146 - t73 + t86) * t65 + (t87 - t61 / 0.2e1 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t103 - t144) * t68) * qJD(2) + (-t18 * pkin(3) - t7 * t12 - t8 * t13 - t26 * t38) * m(5) + (-t64 * t19 + t67 * t20 + m(5) * t140 + (-m(5) * t83 - t64 * t35 - t67 * t36) * qJD(4)) * pkin(7); t58 - t26 * (mrSges(5,1) * t49 + mrSges(5,2) * t48) + (Ifges(5,1) * t48 - t123) * t133 + t15 * t132 + (Ifges(5,5) * t48 - Ifges(5,6) * t49) * t131 - t7 * t35 + t8 * t36 + (t48 * t7 + t49 * t8) * mrSges(5,3) + (-Ifges(5,2) * t49 + t16 + t47) * t134 - t139;];
tauc = t10(:);
