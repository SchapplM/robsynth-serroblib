% Calculate time derivative of joint inertia matrix for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [6x6]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRRPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:20
% EndTime: 2019-03-08 21:04:23
% DurationCPUTime: 1.54s
% Computational Cost: add. (1938->270), mult. (4769->396), div. (0->0), fcn. (4413->10), ass. (0->131)
t154 = m(5) + m(6);
t84 = cos(qJ(6));
t119 = qJD(6) * t84;
t77 = sin(pkin(11));
t79 = cos(pkin(11));
t82 = sin(qJ(3));
t85 = cos(qJ(3));
t57 = t77 * t85 + t79 * t82;
t52 = t57 * qJD(3);
t56 = t77 * t82 - t79 * t85;
t81 = sin(qJ(6));
t95 = t56 * t119 + t81 * t52;
t128 = (mrSges(6,1) + mrSges(5,3));
t153 = 2 * t128;
t145 = pkin(4) + pkin(9);
t121 = qJD(3) * t85;
t122 = qJD(3) * t82;
t53 = t79 * t121 - t77 * t122;
t75 = pkin(3) * t122;
t93 = -qJ(5) * t53 - qJD(5) * t57 + t75;
t11 = t145 * t52 + t93;
t125 = -qJ(4) - pkin(8);
t112 = qJD(3) * t125;
t51 = qJD(4) * t85 + t82 * t112;
t91 = -qJD(4) * t82 + t85 * t112;
t24 = t51 * t77 - t79 * t91;
t14 = pkin(5) * t53 + t24;
t73 = -pkin(3) * t85 - pkin(2);
t97 = -qJ(5) * t57 + t73;
t23 = t145 * t56 + t97;
t63 = t125 * t82;
t65 = t125 * t85;
t41 = -t79 * t63 - t65 * t77;
t26 = pkin(5) * t57 + t41;
t6 = -t23 * t81 + t26 * t84;
t1 = t6 * qJD(6) + t11 * t84 + t14 * t81;
t7 = t23 * t84 + t26 * t81;
t2 = -t7 * qJD(6) - t11 * t81 + t14 * t84;
t151 = t1 * t81 + t2 * t84;
t78 = sin(pkin(6));
t86 = cos(qJ(2));
t135 = t78 * t86;
t83 = sin(qJ(2));
t136 = t78 * t83;
t80 = cos(pkin(6));
t54 = -t82 * t136 + t80 * t85;
t55 = t85 * t136 + t80 * t82;
t28 = -t79 * t54 + t55 * t77;
t19 = t81 * t135 + t84 * t28;
t124 = qJD(2) * t83;
t116 = t78 * t124;
t123 = qJD(2) * t86;
t115 = t78 * t123;
t39 = -t55 * qJD(3) - t82 * t115;
t40 = t54 * qJD(3) + t85 * t115;
t12 = -t79 * t39 + t40 * t77;
t96 = t84 * t135 - t81 * t28;
t3 = t96 * qJD(6) - t81 * t116 + t84 * t12;
t4 = t19 * qJD(6) + t84 * t116 + t81 * t12;
t150 = (t19 * t81 + t84 * t96) * qJD(6) - t84 * t3 - t81 * t4;
t149 = m(5) * pkin(3);
t148 = t56 / 0.2e1;
t147 = -t81 / 0.2e1;
t146 = t84 / 0.2e1;
t144 = pkin(3) * t77;
t143 = pkin(3) * t79;
t140 = mrSges(7,3) * t56;
t139 = Ifges(7,4) * t81;
t138 = Ifges(7,4) * t84;
t137 = Ifges(7,6) * t57;
t13 = t39 * t77 + t40 * t79;
t29 = t54 * t77 + t55 * t79;
t5 = t29 * t13;
t104 = Ifges(7,1) * t81 + t138;
t22 = t57 * Ifges(7,5) + t104 * t56;
t134 = t81 * t22;
t67 = Ifges(7,1) * t84 - t139;
t132 = t81 * t67;
t103 = Ifges(7,2) * t84 + t139;
t21 = t103 * t56 + t137;
t131 = t84 * t21;
t130 = t84 * t52;
t66 = -Ifges(7,2) * t81 + t138;
t129 = t84 * t66;
t127 = -mrSges(5,2) + mrSges(6,3);
t126 = -mrSges(6,2) + mrSges(5,1);
t120 = qJD(6) * t81;
t118 = 0.2e1 * t52;
t117 = 0.2e1 * t82;
t114 = t56 * t120;
t72 = -pkin(4) - t143;
t111 = t95 * Ifges(7,5) + Ifges(7,6) * t130 + Ifges(7,3) * t53;
t110 = t78 ^ 2 * t83 * t123;
t106 = t6 * t81 - t7 * t84;
t105 = -mrSges(4,1) * t85 + mrSges(4,2) * t82;
t64 = mrSges(7,1) * t81 + mrSges(7,2) * t84;
t102 = -Ifges(7,5) * t81 - Ifges(7,6) * t84;
t25 = t79 * t51 + t77 * t91;
t42 = t63 * t77 - t65 * t79;
t100 = t24 * t41 + t25 * t42;
t34 = mrSges(7,1) * t57 - t81 * t140;
t35 = -mrSges(7,2) * t57 + t84 * t140;
t99 = -t81 * t34 + t84 * t35;
t70 = qJ(5) + t144;
t98 = qJD(5) * t29 + t13 * t70;
t94 = t114 - t130;
t90 = t41 * t12 + t42 * t13 + t24 * t28 + t25 * t29;
t88 = -t39 * t82 + t40 * t85 + (-t54 * t85 - t55 * t82) * qJD(3);
t87 = m(7) * t150;
t69 = -pkin(9) + t72;
t62 = t104 * qJD(6);
t61 = t103 * qJD(6);
t60 = (mrSges(4,1) * t82 + mrSges(4,2) * t85) * qJD(3);
t59 = -mrSges(7,1) * t119 + mrSges(7,2) * t120;
t49 = t53 * mrSges(6,3);
t48 = t53 * mrSges(5,2);
t38 = -mrSges(6,2) * t56 - mrSges(6,3) * t57;
t37 = mrSges(5,1) * t56 + mrSges(5,2) * t57;
t33 = pkin(4) * t56 + t97;
t32 = (-mrSges(7,1) * t84 + mrSges(7,2) * t81) * t56;
t31 = -t52 * mrSges(6,2) - t49;
t30 = t52 * mrSges(5,1) + t48;
t27 = -pkin(5) * t56 + t42;
t18 = pkin(4) * t52 + t93;
t17 = mrSges(7,1) * t53 - t95 * mrSges(7,3);
t16 = -mrSges(7,2) * t53 - t94 * mrSges(7,3);
t15 = -t52 * pkin(5) + t25;
t10 = t94 * mrSges(7,1) + t95 * mrSges(7,2);
t9 = t95 * Ifges(7,1) - t94 * Ifges(7,4) + t53 * Ifges(7,5);
t8 = t95 * Ifges(7,4) - t94 * Ifges(7,2) + t53 * Ifges(7,6);
t20 = [0.2e1 * m(7) * (t19 * t3 - t4 * t96 + t5) + 0.2e1 * m(4) * (t54 * t39 + t55 * t40 - t110) + 0.2e1 * t154 * (t28 * t12 - t110 + t5); t29 * t10 + t13 * t32 - t96 * t16 + t19 * t17 + t3 * t34 + t4 * t35 + t88 * mrSges(4,3) + ((-t30 - t31 - t60) * t86 + (-t86 * mrSges(3,2) + (-mrSges(3,1) + t105 + t37 + t38) * t83) * qJD(2)) * t78 + m(5) * ((t73 * t124 - t86 * t75) * t78 + t90) + m(6) * ((t33 * t124 - t18 * t86) * t78 + t90) + m(7) * (-t1 * t96 + t13 * t27 + t15 * t29 + t19 * t2 + t3 * t6 + t4 * t7) + t128 * (t12 * t57 - t13 * t56 + t28 * t53 - t29 * t52) + (-pkin(2) * t116 + t88 * pkin(8)) * m(4); -0.2e1 * pkin(2) * t60 + 0.2e1 * t1 * t35 + 0.2e1 * t27 * t10 + 0.2e1 * t15 * t32 + 0.2e1 * t7 * t16 + 0.2e1 * t6 * t17 + 0.2e1 * t18 * t38 + 0.2e1 * t2 * t34 + 0.2e1 * t73 * t30 + 0.2e1 * t33 * t31 + t53 * t41 * t153 + (-Ifges(4,4) * t82 + pkin(3) * t37) * qJD(3) * t117 + 0.2e1 * m(5) * (t73 * t75 + t100) + 0.2e1 * m(6) * (t18 * t33 + t100) + 0.2e1 * m(7) * (t1 * t7 + t15 * t27 + t2 * t6) + (-t153 * t42 + t131 + t134) * t52 + (0.2e1 * Ifges(4,4) * t85 + (Ifges(4,1) - Ifges(4,2)) * t117) * t121 + ((-Ifges(5,4) - Ifges(6,6)) * t118 + t24 * t153 + ((2 * Ifges(5,1)) + (2 * Ifges(6,2)) + Ifges(7,3)) * t53 + t111) * t57 + (t84 * t8 + t81 * t9 + (Ifges(6,3) + Ifges(5,2)) * t118 - t25 * t153 + (-0.2e1 * Ifges(5,4) - 0.2e1 * Ifges(6,6) - t102) * t53 + (t84 * t22 + (-t21 - t137) * t81) * qJD(6)) * t56; t39 * mrSges(4,1) - t40 * mrSges(4,2) - t29 * t59 - t126 * t12 + (t64 + t127) * t13 + m(6) * (t12 * t72 + t98) + m(7) * t98 - t69 * t87 + (-t12 * t79 + t13 * t77) * t149 + t150 * mrSges(7,3); t70 * t10 + t15 * t64 - t27 * t59 + t127 * t25 - t126 * t24 + (-t56 * mrSges(6,1) + t32) * qJD(5) + (t69 * t17 - t61 * t148 - t2 * mrSges(7,3) + t9 / 0.2e1) * t84 + (t69 * t16 - t62 * t148 - t1 * mrSges(7,3) - t8 / 0.2e1) * t81 + m(6) * (qJD(5) * t42 + t24 * t72 + t25 * t70) + m(7) * (qJD(5) * t27 + t15 * t70 + t151 * t69) + (-t24 * t79 + t25 * t77) * t149 + (t72 * mrSges(6,1) - mrSges(5,3) * t143 + Ifges(7,5) * t146 + Ifges(7,6) * t147 - Ifges(6,4) + Ifges(5,5)) * t53 + (Ifges(4,5) * t85 - Ifges(4,6) * t82 + t105 * pkin(8)) * qJD(3) + (t129 / 0.2e1 + t132 / 0.2e1 - t70 * mrSges(6,1) + Ifges(6,5) - Ifges(5,6) - mrSges(5,3) * t144) * t52 + (-t131 / 0.2e1 - t134 / 0.2e1 + t57 * t102 / 0.2e1 + (t67 * t146 + t66 * t147) * t56 + t106 * mrSges(7,3) + (-m(7) * t106 + t99) * t69) * qJD(6); -0.2e1 * t59 * t70 + t61 * t81 - t62 * t84 + (-t129 - t132) * qJD(6) + 0.2e1 * (mrSges(6,3) + t64 + (m(6) + m(7)) * t70) * qJD(5); m(7) * (-t81 * t3 + t84 * t4 + (-t19 * t84 + t81 * t96) * qJD(6)) + t154 * t116; m(5) * t75 + t84 * t16 - t81 * t17 + t48 - t49 + t126 * t52 + (-t84 * t34 - t81 * t35) * qJD(6) + m(7) * (t1 * t84 - t2 * t81 + (-t6 * t84 - t7 * t81) * qJD(6)) + m(6) * t18; 0; 0; m(6) * t12 - t87; t53 * mrSges(6,1) + t81 * t16 + t84 * t17 + t99 * qJD(6) + m(7) * (-t106 * qJD(6) + t151) + m(6) * t24; 0; 0; 0; mrSges(7,1) * t3 - mrSges(7,2) * t4; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t114 + t111; (-t64 * t69 + t102) * qJD(6); t59; -t64 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t20(1) t20(2) t20(4) t20(7) t20(11) t20(16); t20(2) t20(3) t20(5) t20(8) t20(12) t20(17); t20(4) t20(5) t20(6) t20(9) t20(13) t20(18); t20(7) t20(8) t20(9) t20(10) t20(14) t20(19); t20(11) t20(12) t20(13) t20(14) t20(15) t20(20); t20(16) t20(17) t20(18) t20(19) t20(20) t20(21);];
Mq  = res;
