% Calculate time derivative of joint inertia matrix for
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% MqD [5x5]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPRR12_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:04
% EndTime: 2019-12-31 20:29:09
% DurationCPUTime: 1.77s
% Computational Cost: add. (1800->252), mult. (3929->359), div. (0->0), fcn. (3172->6), ass. (0->116)
t64 = sin(qJ(5));
t67 = cos(qJ(5));
t102 = t64 ^ 2 + t67 ^ 2;
t66 = sin(qJ(2));
t100 = qJD(2) * t66;
t65 = sin(qJ(4));
t69 = cos(qJ(2));
t110 = t65 * t69;
t68 = cos(qJ(4));
t95 = qJD(4) * t68;
t99 = qJD(2) * t69;
t21 = -qJD(4) * t110 - t68 * t100 + t65 * t99 + t66 * t95;
t31 = t65 * t66 + t68 * t69;
t22 = (qJD(2) - qJD(4)) * t31;
t109 = t67 * t22;
t32 = t66 * t68 - t110;
t93 = qJD(5) * t64;
t87 = t32 * t93;
t73 = t87 - t109;
t7 = mrSges(6,1) * t21 + t73 * mrSges(6,3);
t122 = pkin(6) - pkin(7);
t40 = t122 * t100;
t48 = t122 * t69;
t81 = qJD(2) * t48;
t88 = t122 * t66;
t23 = t65 * t48 - t68 * t88;
t97 = qJD(4) * t23;
t11 = -t68 * t40 + t65 * t81 - t97;
t103 = qJ(3) * t99 + t66 * qJD(3);
t70 = -pkin(2) - pkin(3);
t25 = t70 * t100 + t103;
t6 = pkin(4) * t21 - pkin(8) * t22 + t25;
t139 = -t69 * pkin(2) - t66 * qJ(3);
t43 = -pkin(1) + t139;
t30 = t69 * pkin(3) - t43;
t15 = pkin(4) * t31 - pkin(8) * t32 + t30;
t24 = t68 * t48 + t65 * t88;
t10 = t15 * t64 + t24 * t67;
t91 = t10 * qJD(5);
t2 = -t11 * t64 + t6 * t67 - t91;
t82 = t2 + t91;
t133 = -m(6) * t82 - t7;
t121 = mrSges(6,3) * t32;
t18 = -mrSges(6,2) * t31 - t64 * t121;
t147 = -qJD(5) * t18 + t133;
t86 = t102 * mrSges(6,3);
t146 = (mrSges(5,2) - t86) * t68;
t145 = m(4) * pkin(6) + mrSges(4,2);
t44 = -t67 * mrSges(6,1) + mrSges(6,2) * t64;
t106 = mrSges(5,1) - t44;
t134 = -m(6) * pkin(4) - t106;
t143 = m(6) * pkin(8);
t41 = -t65 * qJ(3) + t68 * t70;
t26 = t68 * qJD(3) + t41 * qJD(4);
t141 = t26 * mrSges(5,2);
t42 = t68 * qJ(3) + t65 * t70;
t137 = t102 * t68;
t77 = mrSges(6,1) * t64 + mrSges(6,2) * t67;
t34 = t77 * qJD(5);
t45 = Ifges(6,5) * t64 + Ifges(6,6) * t67;
t52 = Ifges(6,6) * t93;
t92 = qJD(5) * t67;
t136 = (t45 / 0.2e1 - Ifges(5,6)) * t21 + Ifges(5,5) * t22 + t23 * t34 - t11 * mrSges(5,2) + t31 * (Ifges(6,5) * t92 - t52) / 0.2e1;
t19 = mrSges(6,1) * t31 - t67 * t121;
t111 = t64 * t22;
t74 = t32 * t92 + t111;
t8 = -mrSges(6,2) * t21 - t74 * mrSges(6,3);
t135 = -qJD(5) * t19 + t8;
t36 = Ifges(6,4) * t92 - Ifges(6,2) * t93;
t37 = Ifges(6,1) * t92 - Ifges(6,4) * t93;
t120 = Ifges(6,4) * t64;
t46 = Ifges(6,2) * t67 + t120;
t119 = Ifges(6,4) * t67;
t47 = Ifges(6,1) * t64 + t119;
t71 = -(t64 * t46 - t67 * t47) * qJD(5) + t67 * t36 + t64 * t37;
t132 = 2 * m(5);
t131 = 0.2e1 * m(6);
t130 = -0.2e1 * pkin(1);
t96 = qJD(4) * t24;
t12 = -t65 * t40 - t68 * t81 + t96;
t129 = 0.2e1 * t12;
t128 = 0.2e1 * t25;
t127 = -0.2e1 * t34;
t126 = 0.2e1 * t43;
t124 = -t32 / 0.2e1;
t123 = -t46 / 0.2e1;
t118 = Ifges(6,5) * t67;
t117 = t12 * t23;
t27 = t65 * qJD(3) + t42 * qJD(4);
t116 = t23 * t27;
t115 = t26 * t65;
t114 = t27 * t68;
t108 = t68 * t34;
t107 = qJD(5) / 0.2e1;
t105 = -Ifges(4,5) + Ifges(3,4);
t104 = Ifges(6,5) * t109 + Ifges(6,3) * t21;
t9 = t15 * t67 - t24 * t64;
t101 = t9 * qJD(5);
t39 = -pkin(8) + t42;
t94 = qJD(5) * t39;
t89 = 0.2e1 * t69;
t85 = t102 * t26;
t84 = -Ifges(6,6) * t64 - (2 * Ifges(5,4));
t1 = t11 * t67 + t6 * t64 + t101;
t83 = -t1 + t101;
t3 = -t73 * Ifges(6,4) - t74 * Ifges(6,2) + Ifges(6,6) * t21;
t80 = t3 / 0.2e1 + t22 * t47 / 0.2e1;
t4 = -t73 * Ifges(6,1) - t74 * Ifges(6,4) + Ifges(6,5) * t21;
t79 = t4 / 0.2e1 + t22 * t123;
t78 = -t69 * mrSges(4,1) - t66 * mrSges(4,3);
t38 = pkin(4) - t41;
t16 = t77 * t32;
t14 = Ifges(6,5) * t31 + (Ifges(6,1) * t67 - t120) * t32;
t13 = Ifges(6,6) * t31 + (-Ifges(6,2) * t64 + t119) * t32;
t5 = t74 * mrSges(6,1) - t73 * mrSges(6,2);
t17 = [t14 * t109 - t13 * t111 + 0.2e1 * t30 * (mrSges(5,1) * t21 + mrSges(5,2) * t22) + 0.2e1 * t23 * t5 + t16 * t129 + 0.2e1 * t1 * t18 + 0.2e1 * t2 * t19 + 0.2e1 * t9 * t7 + 0.2e1 * t10 * t8 + 0.2e1 * (-t21 * t24 + t22 * t23) * mrSges(5,3) + (t11 * t24 + t25 * t30 + t117) * t132 + (t1 * t10 + t2 * t9 + t117) * t131 + (mrSges(5,1) * t128 - 0.2e1 * t11 * mrSges(5,3) + t84 * t22 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t21 + t104) * t31 + (mrSges(5,2) * t128 + mrSges(5,3) * t129 + 0.2e1 * Ifges(5,1) * t22 - t64 * t3 + t67 * t4 + (t84 + t118) * t21 + (-t67 * t13 - t64 * t14 - t31 * t45) * qJD(5)) * t32 + ((mrSges(3,2) * t130 - 0.2e1 * t43 * mrSges(4,3) + t105 * t89) * t69 + (mrSges(3,1) * t130 + mrSges(4,1) * t126 - 0.2e1 * t105 * t66 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t89) * t66) * qJD(2) + (m(4) * t126 + 0.2e1 * t78) * (pkin(2) * t100 - t103); t38 * t5 + t27 * t16 + t106 * t12 + m(5) * (t11 * t42 - t12 * t41 + t24 * t26 + t116) + m(6) * (t12 * t38 + t116) + (-t21 * t42 - t22 * t41 - t26 * t31 + t27 * t32) * mrSges(5,3) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t69 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t66 + (m(4) * t139 - t69 * mrSges(3,1) + t66 * mrSges(3,2) + t78) * pkin(6)) * qJD(2) + (-qJD(5) * t14 / 0.2e1 + m(6) * (t1 * t39 + t10 * t26 - t9 * t94) + t39 * t8 + t26 * t18 - t19 * t94 + (-t37 / 0.2e1 + t46 * t107) * t32 + t83 * mrSges(6,3) - t80) * t67 + (t13 * t107 - t18 * t94 + (t36 / 0.2e1 + t47 * t107) * t32 + t82 * mrSges(6,3) - t79 + (-m(6) * t9 - t19) * t26 + t133 * t39) * t64 - t136 + t145 * qJD(3) * t69; t38 * t127 + 0.2e1 * t141 + (t27 * t38 + t39 * t85) * t131 + (t26 * t42 - t27 * t41) * t132 + t71 + 0.2e1 * t106 * t27 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) - 0.2e1 * t86 * t26; t145 * t99 + (-t22 * mrSges(5,3) - t5 - m(6) * t12 + m(5) * (-t12 + t96) + (-t31 * mrSges(5,3) + t67 * t18 - t64 * t19 + m(6) * (t10 * t67 - t64 * t9)) * qJD(4)) * t68 + (qJD(4) * t16 + (qJD(4) * t32 - t21) * mrSges(5,3) + m(6) * (-t9 * t92 + t97) + m(5) * (t11 + t97) + (m(6) * t1 + t135) * t67 + t147 * t64) * t65; t108 + m(6) * (t102 * t115 - t114) + m(5) * (-t114 + t115) + (t106 * t65 + t146 + m(6) * (t137 * t39 + t38 * t65) + m(5) * (-t41 * t65 + t42 * t68)) * qJD(4); (-0.1e1 + t102) * t65 * t95 * t131; -pkin(4) * t5 + t134 * t12 + (t32 * t37 / 0.2e1 + t1 * mrSges(6,3) + (t14 / 0.2e1 + t32 * t123 - t9 * mrSges(6,3)) * qJD(5) + (-m(6) * t83 + t135) * pkin(8) + t80) * t67 + (t36 * t124 - t2 * mrSges(6,3) + (-t13 / 0.2e1 + t47 * t124 - t10 * mrSges(6,3)) * qJD(5) + t147 * pkin(8) + t79) * t64 + t136; -t141 + (t38 + pkin(4)) * t34 + (mrSges(6,3) + t143) * t85 + t134 * t27 - t71; -t108 + (t134 * t65 + t137 * t143 - t146) * qJD(4); pkin(4) * t127 + t71; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,5) * t87 - t74 * Ifges(6,6) + t104; t52 - t77 * t26 + (t44 * t39 - t118) * qJD(5); (t65 * t93 - t67 * t95) * mrSges(6,2) + (-t64 * t95 - t65 * t92) * mrSges(6,1); -t52 + (t44 * pkin(8) + t118) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t17(1), t17(2), t17(4), t17(7), t17(11); t17(2), t17(3), t17(5), t17(8), t17(12); t17(4), t17(5), t17(6), t17(9), t17(13); t17(7), t17(8), t17(9), t17(10), t17(14); t17(11), t17(12), t17(13), t17(14), t17(15);];
Mq = res;
