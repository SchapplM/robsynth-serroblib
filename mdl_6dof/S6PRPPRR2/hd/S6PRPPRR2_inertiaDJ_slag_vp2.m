% Calculate time derivative of joint inertia matrix for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S6PRPPRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_inertiaDJ_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:44
% EndTime: 2019-03-08 19:17:47
% DurationCPUTime: 1.18s
% Computational Cost: add. (1044->211), mult. (2802->342), div. (0->0), fcn. (2457->10), ass. (0->112)
t56 = sin(qJ(6));
t59 = cos(qJ(6));
t97 = t56 ^ 2 + t59 ^ 2;
t132 = m(7) * (-0.1e1 + t97);
t52 = sin(pkin(11));
t54 = cos(pkin(11));
t58 = sin(qJ(2));
t61 = cos(qJ(2));
t125 = -t52 * t58 + t54 * t61;
t53 = sin(pkin(6));
t24 = t125 * t53;
t55 = cos(pkin(6));
t57 = sin(qJ(5));
t60 = cos(qJ(5));
t69 = -t24 * t60 - t55 * t57;
t94 = t69 * qJD(5);
t115 = pkin(9) * t57;
t116 = pkin(5) * t60;
t30 = qJD(4) + (t115 + t116) * qJD(5);
t80 = -pkin(2) * t54 - pkin(3);
t44 = -pkin(8) + t80;
t90 = qJD(5) * t60;
t83 = t59 * t90;
t103 = t56 * t57;
t114 = pkin(9) * t60;
t117 = pkin(5) * t57;
t45 = pkin(2) * t52 + qJ(4);
t28 = -t114 + t45 + t117;
t13 = -t44 * t103 + t28 * t59;
t86 = t13 * qJD(6);
t3 = t30 * t56 + t44 * t83 + t86;
t77 = t3 - t86;
t84 = t56 * t90;
t101 = t57 * t59;
t14 = t44 * t101 + t28 * t56;
t85 = t14 * qJD(6);
t4 = t30 * t59 - t44 * t84 - t85;
t130 = -t4 - t85;
t39 = t57 * mrSges(6,1) + t60 * mrSges(6,2);
t129 = mrSges(5,3) + t39;
t128 = t57 ^ 2 - t60 ^ 2;
t92 = qJD(5) * t57;
t81 = t56 * t92;
t87 = qJD(6) * t60;
t82 = t59 * t87;
t127 = t81 - t82;
t25 = (t52 * t61 + t54 * t58) * t53;
t22 = qJD(2) * t25;
t16 = -t24 * t57 + t55 * t60;
t93 = qJD(5) * t16;
t6 = -t22 * t60 + t93;
t112 = t6 * t60;
t5 = t22 * t57 + t94;
t124 = (t16 * t60 - t57 * t69) * qJD(5) + t5 * t57 - t112;
t38 = -mrSges(7,1) * t59 + mrSges(7,2) * t56;
t123 = -m(7) * pkin(5) - mrSges(6,1) + t38;
t122 = 0.2e1 * t44;
t121 = m(6) / 0.2e1;
t120 = m(7) / 0.2e1;
t119 = -t56 / 0.2e1;
t118 = t57 * t90 * t132;
t113 = t69 * t6;
t111 = Ifges(7,4) * t56;
t110 = Ifges(7,4) * t59;
t109 = Ifges(7,5) * t56;
t108 = Ifges(7,6) * t56;
t107 = Ifges(7,6) * t59;
t95 = qJD(2) * t53;
t23 = t125 * t95;
t9 = t25 * t23;
t102 = t57 * Ifges(7,6);
t99 = t60 * mrSges(7,3);
t35 = -mrSges(7,2) * t57 - t56 * t99;
t100 = t59 * t35;
t98 = Ifges(7,6) * t81 + Ifges(7,3) * t90;
t91 = qJD(5) * t59;
t89 = qJD(6) * t56;
t88 = qJD(6) * t59;
t79 = -Ifges(7,5) * t59 + (2 * Ifges(6,4));
t78 = -t6 + t93;
t8 = t16 * t59 + t25 * t56;
t1 = -t8 * qJD(6) + t23 * t59 - t5 * t56;
t7 = -t16 * t56 + t25 * t59;
t2 = t7 * qJD(6) + t23 * t56 + t5 * t59;
t76 = -t1 * t56 + t2 * t59;
t74 = t60 * mrSges(6,1) - t57 * mrSges(6,2);
t73 = mrSges(7,1) * t56 + mrSges(7,2) * t59;
t72 = Ifges(7,1) * t59 - t111;
t41 = Ifges(7,1) * t56 + t110;
t71 = -Ifges(7,2) * t56 + t110;
t40 = Ifges(7,2) * t59 + t111;
t68 = t25 * qJD(4) + t45 * t23;
t67 = t56 * t87 + t57 * t91;
t66 = qJD(5) * t56 * t7 - t8 * t91 + t6;
t64 = t130 * t56 + t77 * t59;
t19 = mrSges(7,1) * t90 + t67 * mrSges(7,3);
t20 = -mrSges(7,2) * t90 + t127 * mrSges(7,3);
t29 = t73 * t60;
t36 = mrSges(7,1) * t57 - t59 * t99;
t63 = qJD(5) * t29 - t56 * t19 + t59 * t20 - t35 * t89 - t36 * t88;
t62 = (-t7 * t88 - t8 * t89 + t76 - t94) * t120 + (t5 - t94) * t121;
t47 = Ifges(7,5) * t88;
t34 = t72 * qJD(6);
t33 = t71 * qJD(6);
t32 = t74 * qJD(5);
t31 = t73 * qJD(6);
t27 = Ifges(7,5) * t57 + t72 * t60;
t26 = t71 * t60 + t102;
t12 = t127 * mrSges(7,1) + t67 * mrSges(7,2);
t11 = -t41 * t87 + (Ifges(7,5) * t60 - t72 * t57) * qJD(5);
t10 = -t40 * t87 + (Ifges(7,6) * t60 - t71 * t57) * qJD(5);
t15 = [0.2e1 * m(7) * (t1 * t7 + t2 * t8 - t113) + 0.2e1 * m(6) * (t16 * t5 - t113 + t9) + 0.2e1 * (m(5) + m(4)) * (-t22 * t24 + t9); t1 * t36 + t69 * t12 + t7 * t19 + t2 * t35 + t8 * t20 + t25 * t32 + t6 * t29 + (-mrSges(4,1) + mrSges(5,2)) * t22 + (-mrSges(3,1) * t58 - mrSges(3,2) * t61) * t95 + (-mrSges(4,2) + t129) * t23 - t124 * mrSges(6,3) + m(7) * (t1 * t13 + t14 * t2 + t3 * t8 + t4 * t7) + m(6) * t68 + m(5) * (t22 * t80 + t68) + m(4) * (-t22 * t54 + t23 * t52) * pkin(2) + ((-t69 * t92 - t112) * t120 + t124 * t121) * t122; 0.2e1 * m(7) * (t13 * t4 + t14 * t3) + 0.2e1 * t3 * t35 + 0.2e1 * t14 * t20 + 0.2e1 * t4 * t36 + 0.2e1 * t13 * t19 + 0.2e1 * t45 * t32 + 0.2e1 * ((m(5) + m(6)) * t45 + t129) * qJD(4) + ((t122 * t29 + t56 * t26 - t59 * t27 + t57 * t79) * qJD(5) + t98) * t57 + (-t56 * t10 + t59 * t11 + t12 * t122 + (t57 * (-t107 - t109) - t56 * t27 - t59 * t26) * qJD(6) + ((-t79 - t108) * t60 + (-0.2e1 * m(7) * t44 ^ 2 - (2 * Ifges(6,1)) + (2 * Ifges(6,2)) + Ifges(7,3)) * t57) * qJD(5)) * t60; 0.2e1 * (t66 * t120 - m(6) * t78 / 0.2e1) * t57 + 0.2e1 * t62 * t60; -t57 * t12 + (m(7) * (-t14 * t101 + t13 * t103 + t128 * t44) - t57 * t100 + t36 * t103) * qJD(5) + (m(7) * t64 + t63) * t60; -0.2e1 * t118; m(5) * t22 + 0.2e1 * (-m(7) * t66 / 0.2e1 + t78 * t121) * t60 + 0.2e1 * t62 * t57; (t12 + (m(7) * (-t13 * t56 + t14 * t59) + t100 - t56 * t36) * qJD(5)) * t60 + (m(7) * (-0.2e1 * t44 * t90 + t64) + t63) * t57; -t128 * qJD(5) * t132; 0.2e1 * t118; -t5 * mrSges(6,2) - t69 * t31 + (m(7) * pkin(9) + mrSges(7,3)) * ((-t56 * t8 - t59 * t7) * qJD(6) + t76) + t123 * t6; pkin(5) * t12 + (t47 / 0.2e1 + (t123 * t44 - Ifges(6,5)) * qJD(5)) * t57 + (t11 / 0.2e1 - t4 * mrSges(7,3) + t40 * t92 / 0.2e1 + (-t102 / 0.2e1 - t26 / 0.2e1 - t14 * mrSges(7,3)) * qJD(6) + (m(7) * t130 - qJD(6) * t35 - t19) * pkin(9)) * t56 + (t10 / 0.2e1 + qJD(6) * t27 / 0.2e1 - t41 * t92 / 0.2e1 + t77 * mrSges(7,3) + (m(7) * t77 - qJD(6) * t36 + t20) * pkin(9)) * t59 + (t59 * t34 / 0.2e1 + t33 * t119 - t44 * t31 + (-t59 * t40 / 0.2e1 + t41 * t119) * qJD(6) + (t109 / 0.2e1 + t107 / 0.2e1 - Ifges(6,6) - t44 * mrSges(6,2)) * qJD(5)) * t60; t57 * t31 + (t60 * t38 + m(7) * (-t115 * t97 - t116) - t97 * t57 * mrSges(7,3) - t74) * qJD(5); -t60 * t31 + (t57 * t38 + m(7) * (t114 * t97 - t117) + t97 * t99 - t39) * qJD(5); -0.2e1 * pkin(5) * t31 + t33 * t59 + t34 * t56 + (-t40 * t56 + t41 * t59) * qJD(6); mrSges(7,1) * t1 - mrSges(7,2) * t2; mrSges(7,1) * t4 - mrSges(7,2) * t3 - Ifges(7,5) * t67 - Ifges(7,6) * t82 + t98; t12; (t57 * t89 - t83) * mrSges(7,2) + (-t57 * t88 - t84) * mrSges(7,1); t47 + (pkin(9) * t38 - t108) * qJD(6); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_6_matlab.m
res = [t15(1) t15(2) t15(4) t15(7) t15(11) t15(16); t15(2) t15(3) t15(5) t15(8) t15(12) t15(17); t15(4) t15(5) t15(6) t15(9) t15(13) t15(18); t15(7) t15(8) t15(9) t15(10) t15(14) t15(19); t15(11) t15(12) t15(13) t15(14) t15(15) t15(20); t15(16) t15(17) t15(18) t15(19) t15(20) t15(21);];
Mq  = res;
