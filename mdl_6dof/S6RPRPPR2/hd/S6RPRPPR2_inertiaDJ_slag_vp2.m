% Calculate time derivative of joint inertia matrix for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPRPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:17
% EndTime: 2019-03-09 02:41:19
% DurationCPUTime: 1.17s
% Computational Cost: add. (1740->220), mult. (3649->321), div. (0->0), fcn. (3211->8), ass. (0->108)
t129 = m(6) + m(5);
t60 = sin(pkin(10));
t61 = cos(pkin(10));
t64 = sin(qJ(3));
t66 = cos(qJ(3));
t42 = t60 * t66 + t61 * t64;
t37 = t42 * qJD(3);
t108 = t61 * t66;
t41 = t60 * t64 - t108;
t63 = sin(qJ(6));
t65 = cos(qJ(6));
t93 = qJD(6) * t65;
t74 = t63 * t37 + t41 * t93;
t128 = 2 * mrSges(6,1) + 2 * mrSges(5,3);
t118 = pkin(4) + pkin(8);
t89 = -cos(pkin(9)) * pkin(1) - pkin(2);
t46 = -pkin(3) * t66 + t89;
t71 = -t42 * qJ(5) + t46;
t14 = t118 * t41 + t71;
t54 = sin(pkin(9)) * pkin(1) + pkin(7);
t97 = qJ(4) + t54;
t39 = t97 * t64;
t40 = t97 * t66;
t22 = t61 * t39 + t40 * t60;
t17 = pkin(5) * t42 + t22;
t3 = -t14 * t63 + t17 * t65;
t96 = qJD(3) * t64;
t38 = qJD(3) * t108 - t60 * t96;
t57 = pkin(3) * t96;
t95 = qJD(5) * t42;
t72 = -qJ(5) * t38 + t57 - t95;
t8 = t118 * t37 + t72;
t86 = qJD(3) * t97;
t29 = qJD(4) * t66 - t64 * t86;
t70 = -t64 * qJD(4) - t66 * t86;
t15 = t29 * t60 - t61 * t70;
t9 = pkin(5) * t38 + t15;
t1 = qJD(6) * t3 + t63 * t9 + t65 * t8;
t4 = t14 * t65 + t17 * t63;
t2 = -qJD(6) * t4 - t63 * t8 + t65 * t9;
t126 = t1 * t63 + t2 * t65;
t13 = pkin(4) * t37 + t72;
t125 = -0.2e1 * t13;
t21 = t41 * pkin(4) + t71;
t124 = -0.2e1 * t21;
t123 = 0.2e1 * t46;
t122 = m(5) * pkin(3);
t121 = t41 / 0.2e1;
t120 = -t63 / 0.2e1;
t119 = t65 / 0.2e1;
t117 = pkin(3) * t60;
t116 = pkin(3) * t61;
t113 = mrSges(7,3) * t41;
t112 = Ifges(7,4) * t63;
t111 = Ifges(7,4) * t65;
t109 = t42 * Ifges(7,6);
t27 = t42 * t38;
t81 = Ifges(7,1) * t63 + t111;
t20 = t42 * Ifges(7,5) + t41 * t81;
t107 = t63 * t20;
t49 = Ifges(7,1) * t65 - t112;
t105 = t63 * t49;
t80 = Ifges(7,2) * t65 + t112;
t19 = t41 * t80 + t109;
t104 = t65 * t19;
t103 = t65 * t37;
t48 = -Ifges(7,2) * t63 + t111;
t102 = t65 * t48;
t100 = mrSges(6,2) - mrSges(5,1);
t34 = t38 * mrSges(5,2);
t35 = t38 * mrSges(6,3);
t99 = t35 - t34;
t98 = t63 ^ 2 + t65 ^ 2;
t94 = qJD(6) * t63;
t92 = 0.2e1 * t66;
t91 = t41 * t94;
t55 = -pkin(4) - t116;
t88 = t98 * t37;
t87 = t74 * Ifges(7,5) + Ifges(7,6) * t103 + Ifges(7,3) * t38;
t83 = t3 * t65 + t4 * t63;
t82 = t3 * t63 - t4 * t65;
t47 = mrSges(7,1) * t63 + mrSges(7,2) * t65;
t79 = -Ifges(7,5) * t63 - Ifges(7,6) * t65;
t16 = t61 * t29 + t60 * t70;
t23 = -t39 * t60 + t40 * t61;
t78 = t15 * t22 + t16 * t23;
t25 = mrSges(7,1) * t42 - t113 * t63;
t26 = -mrSges(7,2) * t42 + t113 * t65;
t77 = t25 * t65 + t26 * t63;
t76 = -t63 * t25 + t65 * t26;
t52 = qJ(5) + t117;
t75 = t38 * t52 + t95;
t73 = t91 - t103;
t68 = -t82 * qJD(6) + t126;
t11 = -mrSges(7,2) * t38 - mrSges(7,3) * t73;
t12 = mrSges(7,1) * t38 - mrSges(7,3) * t74;
t67 = qJD(6) * t76 + t63 * t11 + t65 * t12;
t51 = -pkin(8) + t55;
t45 = t81 * qJD(6);
t44 = t80 * qJD(6);
t43 = -mrSges(7,1) * t93 + mrSges(7,2) * t94;
t24 = (-mrSges(7,1) * t65 + mrSges(7,2) * t63) * t41;
t18 = -pkin(5) * t41 + t23;
t10 = -t37 * pkin(5) + t16;
t7 = mrSges(7,1) * t73 + mrSges(7,2) * t74;
t6 = Ifges(7,1) * t74 - Ifges(7,4) * t73 + t38 * Ifges(7,5);
t5 = Ifges(7,4) * t74 - Ifges(7,2) * t73 + t38 * Ifges(7,6);
t28 = [0.2e1 * t1 * t26 + 0.2e1 * t10 * t24 + 0.2e1 * t4 * t11 + 0.2e1 * t3 * t12 + 0.2e1 * t18 * t7 + 0.2e1 * t2 * t25 + t35 * t124 + t34 * t123 + 0.2e1 * m(5) * t78 + 0.2e1 * m(6) * (t13 * t21 + t78) + 0.2e1 * m(7) * (t1 * t4 + t10 * t18 + t2 * t3) + (mrSges(6,3) * t125 + t128 * t15 + t87) * t42 + (mrSges(6,2) * t125 + t65 * t5 + t63 * t6 - t16 * t128 + (t20 * t65 + (-t19 - t109) * t63) * qJD(6)) * t41 + ((t89 * mrSges(4,2) + Ifges(4,4) * t66) * t92 + (t122 * t123 + 0.2e1 * t89 * mrSges(4,1) + 0.2e1 * pkin(3) * (mrSges(5,1) * t41 + mrSges(5,2) * t42) - 0.2e1 * Ifges(4,4) * t64 + (Ifges(4,1) - Ifges(4,2)) * t92) * t64) * qJD(3) + (mrSges(5,1) * t123 + mrSges(6,2) * t124 - t23 * t128 + t104 + t107 + 0.2e1 * (-Ifges(5,4) - Ifges(6,6)) * t42 + 0.2e1 * (Ifges(5,2) + Ifges(6,3)) * t41) * t37 + (t22 * t128 + ((2 * Ifges(5,1)) + (2 * Ifges(6,2)) + Ifges(7,3)) * t42 + (-0.2e1 * Ifges(5,4) - 0.2e1 * Ifges(6,6) - t79) * t41) * t38; t38 * t24 + t42 * t7 + t77 * t37 + t67 * t41 + m(7) * (t10 * t42 + t18 * t38 + t83 * t37 + t68 * t41) + t129 * (t15 * t41 + t16 * t42 + t22 * t37 + t23 * t38); 0.2e1 * m(7) * (t41 * t88 + t27) + 0.2e1 * t129 * (t37 * t41 + t27); t10 * t47 - t18 * t43 + t52 * t7 + (-mrSges(5,2) + mrSges(6,3)) * t16 + t100 * t15 + (-t41 * mrSges(6,1) + t24) * qJD(5) + (t51 * t12 - t44 * t121 - t2 * mrSges(7,3) + t6 / 0.2e1) * t65 + (t51 * t11 - t45 * t121 - t1 * mrSges(7,3) - t5 / 0.2e1) * t63 + m(6) * (qJD(5) * t23 + t15 * t55 + t16 * t52) + m(7) * (qJD(5) * t18 + t10 * t52 + t126 * t51) + (-t15 * t61 + t16 * t60) * t122 + (t55 * mrSges(6,1) - mrSges(5,3) * t116 + Ifges(7,5) * t119 + Ifges(7,6) * t120 - Ifges(6,4) + Ifges(5,5)) * t38 + (Ifges(4,5) * t66 - Ifges(4,6) * t64 + (-mrSges(4,1) * t66 + mrSges(4,2) * t64) * t54) * qJD(3) + (t102 / 0.2e1 + t105 / 0.2e1 + Ifges(6,5) - Ifges(5,6) - t52 * mrSges(6,1) - mrSges(5,3) * t117) * t37 + (-t104 / 0.2e1 - t107 / 0.2e1 + t42 * t79 / 0.2e1 + (t119 * t49 + t120 * t48) * t41 + t82 * mrSges(7,3) + (-m(7) * t82 + t76) * t51) * qJD(6); t38 * t47 - t42 * t43 + (-mrSges(4,1) * t64 - mrSges(4,2) * t66) * qJD(3) + (-mrSges(7,3) * t98 + t100) * t37 + m(6) * (t37 * t55 + t75) + m(7) * (t51 * t88 + t75) + (-t37 * t61 + t38 * t60) * t122 + t99; -0.2e1 * t43 * t52 + t44 * t63 - t45 * t65 + (-t102 - t105) * qJD(6) + 0.2e1 * (mrSges(6,3) + t47 + (m(6) + m(7)) * t52) * qJD(5); m(5) * t57 + t65 * t11 - t63 * t12 - t100 * t37 - t77 * qJD(6) + m(7) * (-qJD(6) * t83 + t1 * t65 - t2 * t63) + m(6) * t13 - t99; 0; 0; 0; m(6) * t15 + m(7) * t68 + t38 * mrSges(6,1) + t67; 0.2e1 * (m(6) / 0.2e1 + m(7) * t98 / 0.2e1) * t37; 0; 0; 0; mrSges(7,1) * t2 - mrSges(7,2) * t1 - Ifges(7,6) * t91 + t87; -t7; (-t47 * t51 + t79) * qJD(6); t43; -t47 * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t28(1) t28(2) t28(4) t28(7) t28(11) t28(16); t28(2) t28(3) t28(5) t28(8) t28(12) t28(17); t28(4) t28(5) t28(6) t28(9) t28(13) t28(18); t28(7) t28(8) t28(9) t28(10) t28(14) t28(19); t28(11) t28(12) t28(13) t28(14) t28(15) t28(20); t28(16) t28(17) t28(18) t28(19) t28(20) t28(21);];
Mq  = res;
