% Calculate time derivative of joint inertia matrix for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6RPPRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:44:00
% EndTime: 2019-03-09 01:44:02
% DurationCPUTime: 1.45s
% Computational Cost: add. (1610->209), mult. (3163->318), div. (0->0), fcn. (2745->8), ass. (0->98)
t58 = sin(qJ(6));
t60 = cos(qJ(6));
t98 = t58 ^ 2 + t60 ^ 2;
t118 = cos(qJ(4));
t56 = sin(pkin(10));
t59 = sin(qJ(4));
t97 = cos(pkin(10));
t36 = t118 * t56 + t59 * t97;
t33 = t36 * qJD(4);
t75 = t97 * t118;
t35 = t56 * t59 - t75;
t94 = qJD(6) * t58;
t66 = -t60 * t33 + t35 * t94;
t22 = t35 * t33;
t134 = -Ifges(5,1) + Ifges(5,2);
t46 = sin(pkin(9)) * pkin(1) + qJ(3);
t40 = t59 * pkin(4) + t46;
t15 = pkin(5) * t36 + pkin(8) * t35 + t40;
t68 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t65 = qJ(5) - t68;
t34 = t65 * t59;
t64 = t118 * t68;
t63 = -qJ(5) * t118 + t64;
t17 = -t34 * t97 + t56 * t63;
t5 = t15 * t60 - t58 * t17;
t6 = t58 * t15 + t17 * t60;
t93 = qJD(6) * t60;
t133 = -t5 * t93 - t6 * t94;
t122 = m(6) * pkin(4);
t132 = t122 * t56;
t96 = qJD(4) * t59;
t32 = -qJD(4) * t75 + t56 * t96;
t131 = -t32 * t56 - t97 * t33;
t67 = t58 * t33 + t35 * t93;
t8 = -t32 * mrSges(7,1) - mrSges(7,3) * t66;
t9 = t32 * mrSges(7,2) + mrSges(7,3) * t67;
t130 = -t58 * t8 + t60 * t9;
t21 = t36 * t32;
t129 = m(6) * (-t21 + t22);
t41 = -mrSges(7,1) * t60 + t58 * mrSges(7,2);
t49 = -pkin(4) * t97 - pkin(5);
t128 = m(7) * t49 + t41;
t48 = pkin(4) * t56 + pkin(8);
t127 = (-m(7) * t48 - mrSges(7,3)) * t98;
t126 = t32 * t35 - t33 * t36;
t125 = t122 * t97 - t128;
t124 = -2 * mrSges(6,3);
t116 = Ifges(7,4) * t58;
t115 = Ifges(7,4) * t60;
t114 = Ifges(7,6) * t58;
t24 = qJD(4) * t63 - t59 * qJD(5);
t62 = -qJD(5) * t118 + t65 * t96;
t10 = t24 * t56 - t62 * t97;
t16 = -t34 * t56 - t63 * t97;
t113 = t10 * t16;
t111 = t32 * t58;
t110 = t32 * t60;
t109 = t33 * t16;
t108 = t35 * t58;
t107 = t35 * t60;
t104 = t58 * mrSges(7,3);
t42 = Ifges(7,2) * t60 + t116;
t102 = t58 * t42;
t101 = t60 * mrSges(7,3);
t43 = Ifges(7,1) * t58 + t115;
t99 = t60 * t43;
t95 = qJD(6) * t35;
t78 = qJD(4) * t118;
t44 = pkin(4) * t78 + qJD(3);
t92 = 0.2e1 * t32 * mrSges(6,3);
t91 = t35 * t124;
t86 = mrSges(5,1) * t96;
t79 = -t32 * mrSges(6,1) - t33 * mrSges(6,2);
t11 = t24 * t97 + t56 * t62;
t14 = -pkin(5) * t32 + pkin(8) * t33 + t44;
t1 = qJD(6) * t5 + t11 * t60 + t58 * t14;
t2 = -qJD(6) * t6 - t58 * t11 + t14 * t60;
t76 = t1 * t60 - t2 * t58;
t74 = mrSges(7,1) * t58 + mrSges(7,2) * t60;
t73 = Ifges(7,1) * t60 - t116;
t72 = -Ifges(7,2) * t58 + t115;
t71 = t10 * t36 - t16 * t32;
t70 = t35 * t10 + t109;
t19 = -mrSges(7,2) * t36 + t104 * t35;
t20 = t36 * mrSges(7,1) + t101 * t35;
t69 = t60 * t19 - t58 * t20;
t61 = t66 * Ifges(7,5) + Ifges(7,6) * t67 - Ifges(7,3) * t32;
t50 = Ifges(7,5) * t93;
t39 = t73 * qJD(6);
t38 = t72 * qJD(6);
t37 = t74 * qJD(6);
t18 = t74 * t35;
t13 = Ifges(7,5) * t36 - t35 * t73;
t12 = Ifges(7,6) * t36 - t35 * t72;
t7 = mrSges(7,1) * t67 - mrSges(7,2) * t66;
t4 = Ifges(7,1) * t66 + Ifges(7,4) * t67 - Ifges(7,5) * t32;
t3 = Ifges(7,4) * t66 + Ifges(7,2) * t67 - Ifges(7,6) * t32;
t23 = [t17 * t92 + 0.2e1 * Ifges(6,1) * t22 + (-0.2e1 * Ifges(5,4) * t118 + t134 * t59) * t78 + (0.2e1 * Ifges(5,4) * t59 + t118 * t134) * t96 + t66 * t13 + t67 * t12 + 0.2e1 * t44 * (mrSges(6,1) * t36 - mrSges(6,2) * t35) - 0.2e1 * t16 * t7 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t5 * t8 + 0.2e1 * t6 * t9 - t4 * t107 + 0.2e1 * t40 * t79 + 0.2e1 * t46 * (mrSges(5,1) * t118 - mrSges(5,2) * t59) * qJD(4) - 0.2e1 * t126 * Ifges(6,4) + (t11 * t36 + t109) * t124 + t36 * t61 + t3 * t108 + (t91 - 0.2e1 * t18) * t10 + (-(-Ifges(7,5) * t60 + t114) * t35 + (-(2 * Ifges(6,2)) - Ifges(7,3)) * t36) * t32 + 0.2e1 * m(6) * (t11 * t17 + t40 * t44 + t113) + 0.2e1 * m(7) * (t1 * t6 + t2 * t5 + t113) + 0.2e1 * (t59 * mrSges(5,1) + mrSges(5,2) * t118 + mrSges(4,3) + (m(5) + m(4)) * t46) * qJD(3); t32 * t18 - t36 * t7 + m(6) * (-t11 * t35 - t17 * t33 + t71) + m(7) * t71 + (-t33 * t19 - t35 * t9 + t20 * t95 + m(7) * (-t1 * t35 - t33 * t6 + t5 * t95)) * t60 + (t19 * t95 + t33 * t20 + t35 * t8 + m(7) * (t2 * t35 + t33 * t5 + t6 * t95)) * t58; 0.2e1 * t129 + 0.2e1 * m(7) * (t22 * t98 - t21); -t35 * t7 + (-t18 + t91) * t33 - t69 * t32 + m(7) * (-t6 * t110 + t5 * t111 + t70) + m(6) * (-t17 * t32 + t70) + (t92 + (-t58 * t19 - t60 * t20) * qJD(6) + m(7) * (t76 + t133) + m(6) * t11 + t130) * t36; m(7) * t126 * (-0.1e1 + t98); 0.2e1 * t129 + 0.2e1 * m(7) * (-t21 * t98 + t22); -t32 * (Ifges(7,5) * t58 + Ifges(7,6) * t60) / 0.2e1 + t60 * t3 / 0.2e1 + t58 * t4 / 0.2e1 + t16 * t37 - t49 * t7 + Ifges(6,6) * t32 - t39 * t107 / 0.2e1 - t2 * t104 - Ifges(5,5) * t96 - t12 * t94 / 0.2e1 + t36 * (-Ifges(7,6) * t94 + t50) / 0.2e1 - t68 * t86 - Ifges(5,6) * t78 - qJD(4) * mrSges(5,2) * t64 + t1 * t101 + (qJD(6) * t43 + t38) * t108 / 0.2e1 + (t35 * t42 + t13) * t93 / 0.2e1 - t131 * mrSges(6,3) * pkin(4) + (-mrSges(6,2) + t132) * t11 + t133 * mrSges(7,3) + (-Ifges(6,5) + t102 / 0.2e1 - t99 / 0.2e1) * t33 + (-mrSges(6,1) - t125) * t10 + (m(7) * ((-t5 * t60 - t58 * t6) * qJD(6) + t76) - t20 * t93 - t19 * t94 + t130) * t48; -mrSges(5,1) * t78 + mrSges(5,2) * t96 + t36 * t37 - t79 + t125 * t32 + (t127 - t132) * t33; -mrSges(5,2) * t78 - t86 + t131 * t122 + t35 * t37 + (-mrSges(6,1) + t128) * t33 + (mrSges(6,2) + t127) * t32; 0.2e1 * t49 * t37 + t38 * t60 + t58 * t39 + (t99 - t102) * qJD(6); t58 * t9 + t60 * t8 + t69 * qJD(6) + m(7) * (t58 * t1 + t2 * t60 + (-t5 * t58 + t6 * t60) * qJD(6)) + m(6) * t44 + t79; 0; 0; 0; 0; t2 * mrSges(7,1) - t1 * mrSges(7,2) + t61; t7; (t36 * t94 + t110) * mrSges(7,2) + (-t36 * t93 + t111) * mrSges(7,1); t50 + (t41 * t48 - t114) * qJD(6); -t37; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t23(1) t23(2) t23(4) t23(7) t23(11) t23(16); t23(2) t23(3) t23(5) t23(8) t23(12) t23(17); t23(4) t23(5) t23(6) t23(9) t23(13) t23(18); t23(7) t23(8) t23(9) t23(10) t23(14) t23(19); t23(11) t23(12) t23(13) t23(14) t23(15) t23(20); t23(16) t23(17) t23(18) t23(19) t23(20) t23(21);];
Mq  = res;
