% Calculate time derivative of joint inertia matrix for
% S6RPPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP1_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:57:43
% EndTime: 2019-03-09 01:57:45
% DurationCPUTime: 1.67s
% Computational Cost: add. (2280->298), mult. (4996->417), div. (0->0), fcn. (4516->8), ass. (0->135)
t138 = Ifges(6,5) + Ifges(7,5);
t164 = Ifges(6,3) + Ifges(7,3);
t92 = cos(qJ(5));
t124 = qJD(5) * t92;
t87 = sin(pkin(10));
t88 = cos(pkin(10));
t91 = sin(qJ(4));
t93 = cos(qJ(4));
t61 = t87 * t93 + t91 * t88;
t119 = t61 * t124;
t60 = t87 * t91 - t93 * t88;
t56 = t60 * qJD(4);
t90 = sin(qJ(5));
t146 = t56 * t90;
t98 = t119 - t146;
t76 = sin(pkin(9)) * pkin(1) + qJ(3);
t153 = pkin(7) + t76;
t58 = t153 * t87;
t59 = t153 * t88;
t37 = -t91 * t58 + t59 * t93;
t32 = t92 * t37;
t99 = -cos(pkin(9)) * pkin(1) - pkin(3) * t88 - pkin(2);
t33 = pkin(4) * t60 - pkin(8) * t61 + t99;
t14 = t90 * t33 + t32;
t113 = (t90 ^ 2 + t92 ^ 2) * t56;
t163 = -t93 * t58 - t59 * t91;
t136 = -qJ(6) - pkin(8);
t68 = t136 * t90;
t71 = t136 * t92;
t162 = -t68 * t90 - t71 * t92;
t161 = 2 * m(5);
t160 = 2 * m(7);
t159 = -2 * mrSges(5,3);
t158 = -2 * mrSges(7,3);
t23 = t61 * qJD(3) + t37 * qJD(4);
t157 = 0.2e1 * t23;
t156 = -0.2e1 * t163;
t155 = m(7) * pkin(5);
t57 = t61 * qJD(4);
t154 = pkin(4) * t57;
t152 = mrSges(6,2) * t92;
t151 = Ifges(6,4) * t90;
t150 = Ifges(6,4) * t92;
t149 = Ifges(7,4) * t90;
t148 = Ifges(7,4) * t92;
t147 = t23 * t163;
t145 = t56 * t92;
t144 = t57 * t60;
t142 = t61 * t90;
t141 = t61 * t92;
t137 = -Ifges(6,6) - Ifges(7,6);
t125 = qJD(5) * t90;
t120 = t61 * t125;
t97 = t120 + t145;
t18 = mrSges(7,1) * t57 + t97 * mrSges(7,3);
t19 = mrSges(6,1) * t57 + t97 * mrSges(6,3);
t135 = t18 + t19;
t20 = -mrSges(7,2) * t57 - t98 * mrSges(7,3);
t21 = -mrSges(6,2) * t57 - t98 * mrSges(6,3);
t134 = t20 + t21;
t102 = -Ifges(7,2) * t90 + t148;
t25 = t60 * Ifges(7,6) + t102 * t61;
t103 = -Ifges(6,2) * t90 + t150;
t26 = t60 * Ifges(6,6) + t103 * t61;
t133 = -t25 - t26;
t104 = Ifges(7,1) * t92 - t149;
t27 = Ifges(7,5) * t60 + t104 * t61;
t105 = Ifges(6,1) * t92 - t151;
t28 = Ifges(6,5) * t60 + t105 * t61;
t132 = t27 + t28;
t41 = -mrSges(7,2) * t60 - mrSges(7,3) * t142;
t42 = -mrSges(6,2) * t60 - mrSges(6,3) * t142;
t131 = t41 + t42;
t43 = mrSges(7,1) * t60 - mrSges(7,3) * t141;
t44 = mrSges(6,1) * t60 - mrSges(6,3) * t141;
t130 = -t43 - t44;
t129 = -mrSges(6,1) * t92 + mrSges(6,2) * t90 - mrSges(5,1);
t62 = mrSges(7,1) * t125 + mrSges(7,2) * t124;
t127 = qJ(6) * t61;
t126 = qJD(5) * t61;
t22 = -t60 * qJD(3) + qJD(4) * t163;
t40 = pkin(8) * t56 + t154;
t123 = t33 * t124 + t92 * t22 + t90 * t40;
t122 = -mrSges(7,1) * t98 + mrSges(7,2) * t145;
t121 = pkin(5) * t125;
t118 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t72 = Ifges(7,2) * t92 + t149;
t73 = Ifges(6,2) * t92 + t151;
t117 = -t72 / 0.2e1 - t73 / 0.2e1;
t74 = Ifges(7,1) * t90 + t148;
t75 = Ifges(6,1) * t90 + t150;
t116 = t74 / 0.2e1 + t75 / 0.2e1;
t115 = -mrSges(6,1) - t155;
t114 = t137 * t90;
t51 = t56 * mrSges(5,2);
t112 = t57 * mrSges(5,1) - t51;
t111 = -t22 * t90 + t92 * t40;
t13 = t92 * t33 - t37 * t90;
t109 = qJD(5) * t136;
t108 = -t138 * t145 + t164 * t57;
t107 = -(2 * Ifges(5,4)) + t114;
t106 = mrSges(6,1) * t90 + t152;
t101 = -t163 * t57 + t23 * t60;
t100 = qJ(6) * t56 - qJD(6) * t61;
t15 = -mrSges(7,2) * t120 - t122;
t96 = t115 * t90 - t152;
t95 = t137 * t92 - t138 * t90;
t94 = t130 * t90 + t131 * t92;
t82 = Ifges(6,5) * t124;
t81 = Ifges(7,5) * t124;
t78 = -pkin(5) * t92 - pkin(4);
t69 = -mrSges(7,1) * t92 + mrSges(7,2) * t90;
t67 = t105 * qJD(5);
t66 = t104 * qJD(5);
t65 = t103 * qJD(5);
t64 = t102 * qJD(5);
t63 = t106 * qJD(5);
t55 = -qJD(6) * t90 + t92 * t109;
t54 = qJD(6) * t92 + t90 * t109;
t39 = t106 * t61;
t38 = (mrSges(7,1) * t90 + mrSges(7,2) * t92) * t61;
t24 = pkin(5) * t142 - t163;
t16 = t98 * mrSges(6,1) - t97 * mrSges(6,2);
t12 = t98 * pkin(5) + t23;
t10 = -t97 * Ifges(6,1) - t98 * Ifges(6,4) + t57 * Ifges(6,5);
t9 = -t97 * Ifges(7,1) - t98 * Ifges(7,4) + t57 * Ifges(7,5);
t8 = -t97 * Ifges(6,4) - t98 * Ifges(6,2) + t57 * Ifges(6,6);
t7 = -t97 * Ifges(7,4) - t98 * Ifges(7,2) + t57 * Ifges(7,6);
t6 = -t90 * t127 + t14;
t5 = pkin(5) * t60 - t92 * t127 + t13;
t4 = -qJD(5) * t14 + t111;
t3 = -t37 * t125 + t123;
t2 = -qJ(6) * t119 + (-qJD(5) * t37 + t100) * t90 + t123;
t1 = pkin(5) * t57 + t100 * t92 + (-t32 + (-t33 + t127) * t90) * qJD(5) + t111;
t11 = [0.2e1 * t99 * t112 + t16 * t156 + 0.2e1 * t12 * t38 + t39 * t157 + 0.2e1 * t2 * t41 + 0.2e1 * t3 * t42 + 0.2e1 * t1 * t43 + 0.2e1 * t4 * t44 + 0.2e1 * t5 * t18 + 0.2e1 * t13 * t19 + 0.2e1 * t6 * t20 + 0.2e1 * t14 * t21 + 0.2e1 * t24 * t15 + t37 * t57 * t159 - (mrSges(5,3) * t156 + t132 * t92 + t133 * t90) * t56 + 0.2e1 * m(6) * (t13 * t4 + t14 * t3 - t147) + (t22 * t37 - t147) * t161 + (t1 * t5 + t12 * t24 + t2 * t6) * t160 + (t22 * t159 + ((2 * Ifges(5,2)) + t164) * t57 - t107 * t56 + t108) * t60 + (mrSges(5,3) * t157 - 0.2e1 * Ifges(5,1) * t56 + (t10 + t9) * t92 + (-t7 - t8) * t90 + (t138 * t92 + t107) * t57 + (-t132 * t90 + t133 * t92 + t95 * t60) * qJD(5)) * t61 + 0.2e1 * (m(4) * t76 + mrSges(4,3)) * qJD(3) * (t87 ^ 2 + t88 ^ 2); (t15 + t16) * t60 + (t38 + t39) * t57 - t94 * t56 + m(7) * (t12 * t60 - t6 * t145 + t5 * t146 + t24 * t57) + m(6) * (t13 * t146 - t14 * t145 + t101) + m(5) * (-t37 * t56 + t101) + (t134 * t92 - t135 * t90 + (t130 * t92 - t131 * t90) * qJD(5) + m(7) * (-t1 * t90 - t5 * t124 - t6 * t125 + t2 * t92) + m(6) * (-t13 * t124 - t14 * t125 + t3 * t92 - t4 * t90) + m(5) * t22) * t61; (-t56 * t61 + t144) * t161 + 0.4e1 * (m(7) / 0.2e1 + m(6) / 0.2e1) * (-t61 * t113 + t144); t135 * t92 + t134 * t90 + t94 * qJD(5) + m(7) * (t1 * t92 + t2 * t90 + (-t5 * t90 + t6 * t92) * qJD(5)) + m(6) * (t3 * t90 + t4 * t92 + (-t13 * t90 + t14 * t92) * qJD(5)) + t112; 0; 0; t68 * t18 + t12 * t69 - t71 * t20 + t78 * t15 + t24 * t62 - t163 * t63 + t55 * t43 - Ifges(5,5) * t56 - Ifges(5,6) * t57 + t54 * t41 - t22 * mrSges(5,2) - pkin(4) * t16 + m(7) * (t1 * t68 + t12 * t78 - t2 * t71 + t5 * t55 + t54 * t6) + (t81 / 0.2e1 + t82 / 0.2e1) * t60 + (-m(6) * pkin(4) + t129) * t23 + (t9 / 0.2e1 + t10 / 0.2e1 - t1 * mrSges(7,3) - t4 * mrSges(6,3) + (-t64 / 0.2e1 - t65 / 0.2e1) * t61 + (Ifges(6,5) / 0.2e1 + Ifges(7,5) / 0.2e1) * t57 - t117 * t56 + (-t25 / 0.2e1 - t26 / 0.2e1 + pkin(5) * t38 - t6 * mrSges(7,3) - t14 * mrSges(6,3) - t116 * t61 + t118 * t60 + t24 * t155) * qJD(5) + (-m(6) * t4 - t19 + (-m(6) * t14 - t42) * qJD(5)) * pkin(8)) * t90 + (t7 / 0.2e1 + t8 / 0.2e1 + t2 * mrSges(7,3) + t3 * mrSges(6,3) + (t66 / 0.2e1 + t67 / 0.2e1) * t61 - t118 * t57 - t116 * t56 + (m(6) * t3 + t21) * pkin(8) + (t28 / 0.2e1 + t27 / 0.2e1 - t5 * mrSges(7,3) - t13 * mrSges(6,3) + t117 * t61 + (-m(6) * t13 - t44) * pkin(8)) * qJD(5)) * t92; t51 + (t62 + t63) * t60 + (t69 + t129) * t57 + m(6) * (-pkin(8) * t113 - t154) - (mrSges(6,3) + mrSges(7,3)) * t113 + (t60 * t121 + t57 * t78 + (t54 * t92 - t55 * t90 + (-t68 * t92 + t71 * t90) * qJD(5)) * t61 - t162 * t56) * m(7); m(7) * (qJD(5) * t162 + t54 * t90 + t55 * t92); (-t54 * t71 + t55 * t68) * t160 - 0.2e1 * pkin(4) * t63 + 0.2e1 * t78 * t62 + (t55 * t158 + t66 + t67 + (-t71 * t158 - t72 - t73 + 0.2e1 * (m(7) * t78 + t69) * pkin(5)) * qJD(5)) * t90 + (0.2e1 * t54 * mrSges(7,3) + t64 + t65 + (t68 * t158 + t74 + t75) * qJD(5)) * t92; mrSges(6,1) * t4 + mrSges(7,1) * t1 - mrSges(6,2) * t3 - mrSges(7,2) * t2 - t56 * t114 + (m(7) * t1 + t18) * pkin(5) + t95 * t126 + t108; -t96 * t56 + (t115 * t92 + (mrSges(6,2) + mrSges(7,2)) * t90) * t126 + t122; t96 * qJD(5) - t62; -mrSges(7,2) * t54 + t81 + t82 + (mrSges(7,1) + t155) * t55 + ((-mrSges(6,1) * pkin(8) - (mrSges(7,3) * pkin(5))) * t92 + (mrSges(6,2) * pkin(8) + t137) * t90) * qJD(5); 0; m(7) * t12 + t15; m(7) * t57; 0; m(7) * t121 + t62; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t11(1) t11(2) t11(4) t11(7) t11(11) t11(16); t11(2) t11(3) t11(5) t11(8) t11(12) t11(17); t11(4) t11(5) t11(6) t11(9) t11(13) t11(18); t11(7) t11(8) t11(9) t11(10) t11(14) t11(19); t11(11) t11(12) t11(13) t11(14) t11(15) t11(20); t11(16) t11(17) t11(18) t11(19) t11(20) t11(21);];
Mq  = res;
