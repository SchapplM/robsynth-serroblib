% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:23
% EndTime: 2019-12-05 15:16:28
% DurationCPUTime: 1.61s
% Computational Cost: add. (1410->225), mult. (3697->341), div. (0->0), fcn. (2561->8), ass. (0->117)
t134 = -pkin(7) - pkin(6);
t87 = sin(qJ(4));
t76 = t134 * t87;
t90 = cos(qJ(4));
t77 = t134 * t90;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t47 = t76 * t86 - t77 * t89;
t84 = sin(pkin(9));
t120 = qJD(1) * t84;
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t67 = qJD(2) * t91 - t88 * t120;
t70 = t86 * t90 + t87 * t89;
t104 = qJD(4) * t134;
t72 = t87 * t104;
t73 = t90 * t104;
t144 = -t47 * qJD(5) + t70 * t67 - t72 * t86 + t73 * t89;
t46 = t76 * t89 + t77 * t86;
t97 = t86 * t87 - t89 * t90;
t143 = t46 * qJD(5) + t97 * t67 + t72 * t89 + t73 * t86;
t58 = qJD(3) * t67;
t141 = -Ifges(5,1) / 0.2e1;
t115 = qJD(3) * t90;
t140 = -Ifges(5,4) * t115 / 0.2e1;
t54 = t97 * t88;
t113 = qJD(4) * t87;
t85 = cos(pkin(9));
t119 = qJD(1) * t85;
t79 = t87 * t119;
t101 = qJD(4) * t79 - t58 * t87;
t112 = qJD(4) * t90;
t68 = qJD(2) * t88 + t91 * t120;
t57 = qJD(3) * pkin(6) + t68;
t23 = -t57 * t112 + t101;
t45 = t57 * t90 - t79;
t139 = -t45 * t113 - t23 * t87;
t83 = qJD(4) + qJD(5);
t64 = t97 * qJD(3);
t136 = -t64 / 0.2e1;
t65 = t70 * qJD(3);
t135 = t65 / 0.2e1;
t81 = -pkin(4) * t90 - pkin(3);
t51 = t81 * qJD(3) - t67;
t133 = m(6) * t51;
t132 = mrSges(6,3) * t64;
t131 = Ifges(5,4) * t87;
t100 = pkin(7) * qJD(3) + t57;
t39 = t100 * t90 - t79;
t129 = t39 * t86;
t128 = t39 * t89;
t59 = qJD(3) * t68;
t127 = t59 * t88;
t126 = t65 * mrSges(6,3);
t125 = t65 * Ifges(6,4);
t124 = t84 * t91;
t37 = mrSges(6,1) * t64 + mrSges(6,2) * t65;
t123 = t37 + (-mrSges(5,1) * t90 + mrSges(5,2) * t87) * qJD(3);
t122 = Ifges(5,5) * qJD(4);
t121 = Ifges(5,6) * qJD(4);
t117 = qJD(3) * t87;
t116 = qJD(3) * t88;
t114 = qJD(3) * t91;
t111 = qJD(5) * t86;
t110 = qJD(5) * t89;
t109 = qJD(3) * qJD(4);
t108 = pkin(4) * t113;
t107 = t90 * t119;
t106 = t84 * t116;
t103 = t122 / 0.2e1;
t102 = -t121 / 0.2e1;
t66 = (mrSges(5,1) * t87 + mrSges(5,2) * t90) * t109;
t40 = t83 * t97;
t34 = t40 * qJD(3);
t41 = t83 * t70;
t35 = t41 * qJD(3);
t9 = mrSges(6,1) * t35 - mrSges(6,2) * t34;
t92 = qJD(3) ^ 2;
t98 = t92 * mrSges(4,2) + t66 + t9;
t38 = -t100 * t87 - t107;
t33 = qJD(4) * pkin(4) + t38;
t7 = t33 * t89 - t129;
t8 = t33 * t86 + t128;
t60 = -t87 * t124 - t85 * t90;
t61 = t90 * t124 - t85 * t87;
t27 = t60 * t89 - t61 * t86;
t28 = t60 * t86 + t61 * t89;
t44 = -t57 * t87 - t107;
t56 = -qJD(3) * pkin(3) - t67;
t96 = t44 * mrSges(5,3) + t117 * t141 - t122 / 0.2e1 + t140 - t56 * mrSges(5,2);
t95 = t45 * mrSges(5,3) + t121 / 0.2e1 + (Ifges(5,2) * t90 + t131) * qJD(3) / 0.2e1 - t56 * mrSges(5,1);
t94 = -t92 * mrSges(4,1) + t123 * qJD(3);
t52 = t90 * t58;
t18 = t38 * qJD(4) + t52;
t19 = -t100 * t112 + t101;
t2 = qJD(5) * t7 + t18 * t89 + t19 * t86;
t25 = -t64 * Ifges(6,2) + t83 * Ifges(6,6) + t125;
t55 = Ifges(6,4) * t64;
t26 = t65 * Ifges(6,1) + t83 * Ifges(6,5) - t55;
t3 = -qJD(5) * t8 - t18 * t86 + t19 * t89;
t93 = t3 * mrSges(6,1) - t2 * mrSges(6,2) - t7 * t132 + t25 * t135 - t51 * (mrSges(6,1) * t65 - mrSges(6,2) * t64) - t65 * (-Ifges(6,1) * t64 - t125) / 0.2e1 - t83 * (-Ifges(6,5) * t64 - Ifges(6,6) * t65) / 0.2e1 - Ifges(6,6) * t35 - Ifges(6,5) * t34 + (-Ifges(6,2) * t65 + t26 - t55) * t64 / 0.2e1;
t75 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t115;
t74 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t117;
t53 = t70 * t88;
t50 = (t68 + t108) * qJD(3);
t49 = mrSges(6,1) * t83 - t126;
t48 = -mrSges(6,2) * t83 - t132;
t43 = t60 * qJD(4) - t90 * t106;
t42 = -t61 * qJD(4) + t87 * t106;
t22 = t44 * qJD(4) + t52;
t13 = -t70 * t114 + t83 * t54;
t12 = -t41 * t88 - t91 * t64;
t11 = t38 * t89 - t129;
t10 = -t38 * t86 - t128;
t5 = -t28 * qJD(5) + t42 * t89 - t43 * t86;
t4 = t27 * qJD(5) + t42 * t86 + t43 * t89;
t1 = [t4 * t48 + t42 * t74 + t43 * t75 + t5 * t49 + (t27 * t34 - t28 * t35) * mrSges(6,3) + (-t60 * t90 - t61 * t87) * mrSges(5,3) * t109 + m(5) * (t22 * t61 + t23 * t60 + t42 * t44 + t43 * t45) + m(6) * (t2 * t28 + t27 * t3 + t4 * t8 + t5 * t7) + (t94 * t91 + t98 * t88 + m(4) * (-t67 * t114 - t68 * t116 + t58 * t91 + t127) + m(5) * (t56 * t114 + t127) + m(6) * (t51 * t114 + t50 * t88)) * t84; t12 * t48 + t13 * t49 + (-t34 * t53 + t35 * t54) * mrSges(6,3) + m(6) * (t12 * t8 + t13 * t7 - t2 * t54 - t3 * t53) + ((-t87 * t74 + t90 * t75) * qJD(3) - t98 - m(6) * t50 + (t45 * t115 - t44 * t117 - t59) * m(5)) * t91 + ((-t74 * t90 - t75 * t87) * qJD(4) + t94 + qJD(3) * t133 + (qJD(3) * t56 - t44 * t112 + t22 * t90 + t139) * m(5)) * t88; t81 * t9 + t83 * (-Ifges(6,5) * t40 - Ifges(6,6) * t41) / 0.2e1 + t50 * (mrSges(6,1) * t97 + mrSges(6,2) * t70) - pkin(3) * t66 - t59 * mrSges(4,1) + t51 * (mrSges(6,1) * t41 - mrSges(6,2) * t40) - t40 * t26 / 0.2e1 - t41 * t25 / 0.2e1 + t144 * t49 + t143 * t48 + (-t41 * t136 + t35 * t97) * Ifges(6,2) + (-t40 * t135 - t34 * t70) * Ifges(6,1) + (qJD(3) * mrSges(4,1) - t123) * t68 + (t59 * mrSges(5,2) - t23 * mrSges(5,3) + t67 * t74 + (-pkin(6) * t75 + pkin(4) * t37 - 0.3e1 / 0.2e1 * Ifges(5,4) * t117 + t102 - t95) * qJD(4)) * t87 - m(5) * (-t44 * t67 * t87 + t56 * t68) + m(5) * (-pkin(3) * t59 + t139 * pkin(6)) + (-t2 * t97 - t3 * t70 + t34 * t46 - t35 * t47 + t40 * t7 - t41 * t8) * mrSges(6,3) + (-t41 * t135 - t40 * t136 + t34 * t97 - t35 * t70) * Ifges(6,4) + (-t59 * mrSges(5,1) + (t103 + (-m(5) * t44 - t74) * pkin(6) + (0.3e1 / 0.2e1 * Ifges(5,4) * t90 + (0.3e1 / 0.2e1 * Ifges(5,1) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t87) * qJD(3) - t96) * qJD(4) + (-t45 * m(5) - t75) * t67 + (pkin(6) * m(5) + mrSges(5,3)) * t22) * t90 + (t2 * t47 + t3 * t46 + t50 * t81 + t143 * t8 + t144 * t7 + (t108 - t68) * t51) * m(6); -m(6) * (t10 * t7 + t11 * t8) + t93 + (t48 * t110 - t49 * t111 + m(6) * (t8 * t110 - t7 * t111 + t2 * t86 + t3 * t89) + (t34 * t89 - t35 * t86) * mrSges(6,3)) * pkin(4) + ((t103 + t140 + t96) * t90 + (t102 + (t131 / 0.2e1 + (Ifges(5,2) / 0.2e1 + t141) * t90) * qJD(3) + (-t37 - t133) * pkin(4) + t95) * t87) * qJD(3) + t45 * t74 - t44 * t75 - t11 * t48 - t10 * t49 + t8 * t126 - t22 * mrSges(5,2) + t23 * mrSges(5,1); t93 - t48 * t7 + (t49 + t126) * t8;];
tauc = t1(:);
