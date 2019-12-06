% Calculate time derivative of joint inertia matrix for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:30:09
% EndTime: 2019-12-05 16:30:15
% DurationCPUTime: 1.43s
% Computational Cost: add. (1254->254), mult. (3433->414), div. (0->0), fcn. (3031->10), ass. (0->117)
t88 = sin(pkin(10));
t90 = cos(pkin(10));
t135 = mrSges(5,1) * t88 + mrSges(5,2) * t90;
t92 = sin(qJ(5));
t95 = cos(qJ(5));
t99 = t88 * t92 - t90 * t95;
t59 = t99 * qJD(5);
t134 = 2 * m(5);
t133 = 2 * m(6);
t132 = 2 * pkin(7);
t87 = t90 ^ 2;
t131 = m(5) / 0.2e1;
t130 = -t99 / 0.2e1;
t70 = t88 * t95 + t90 * t92;
t129 = t70 / 0.2e1;
t128 = t90 / 0.2e1;
t96 = cos(qJ(3));
t107 = qJD(3) * t96;
t55 = t135 * t107;
t60 = t70 * qJD(5);
t93 = sin(qJ(3));
t25 = -t107 * t99 - t60 * t93;
t26 = -t107 * t70 + t59 * t93;
t7 = -mrSges(6,1) * t26 + mrSges(6,2) * t25;
t127 = t7 + t55;
t124 = Ifges(5,4) * t88;
t123 = Ifges(5,4) * t90;
t122 = Ifges(5,2) * t88;
t97 = cos(qJ(2));
t109 = qJD(2) * t97;
t89 = sin(pkin(5));
t103 = t89 * t109;
t94 = sin(qJ(2));
t117 = t89 * t94;
t91 = cos(pkin(5));
t62 = t117 * t96 + t91 * t93;
t40 = qJD(3) * t62 + t103 * t93;
t121 = t40 * t93;
t61 = t117 * t93 - t91 * t96;
t41 = -qJD(3) * t61 + t103 * t96;
t120 = t41 * t96;
t18 = t61 * t40;
t119 = t88 * t93;
t118 = t88 * t96;
t116 = t89 * t97;
t115 = t90 * t93;
t114 = t90 * t96;
t113 = pkin(8) + qJ(4);
t112 = -Ifges(6,5) * t59 - Ifges(6,6) * t60;
t111 = -mrSges(5,1) * t90 + mrSges(5,2) * t88 - mrSges(4,1);
t108 = qJD(3) * t93;
t105 = pkin(7) * t108;
t58 = -qJD(4) * t93 + (pkin(3) * t93 - qJ(4) * t96) * qJD(3);
t36 = t88 * t105 + t58 * t90;
t75 = -pkin(3) * t96 - qJ(4) * t93 - pkin(2);
t50 = pkin(7) * t114 + t75 * t88;
t106 = -Ifges(6,5) * t25 - Ifges(6,6) * t26 - Ifges(6,3) * t108;
t104 = qJD(2) * t117;
t102 = pkin(4) * t88 + pkin(7);
t101 = -Ifges(5,5) * t90 + Ifges(5,6) * t88;
t20 = t104 * t90 - t41 * t88;
t21 = t104 * t88 + t41 * t90;
t100 = -t20 * t88 + t21 * t90;
t68 = t90 * t75;
t32 = -pkin(8) * t115 + t68 + (-pkin(7) * t88 - pkin(4)) * t96;
t42 = -pkin(8) * t119 + t50;
t8 = t32 * t95 - t42 * t92;
t9 = t32 * t92 + t42 * t95;
t38 = -t116 * t90 - t62 * t88;
t39 = -t116 * t88 + t62 * t90;
t10 = t38 * t95 - t39 * t92;
t11 = t38 * t92 + t39 * t95;
t76 = t113 * t88;
t78 = t113 * t90;
t43 = -t76 * t95 - t78 * t92;
t44 = -t76 * t92 + t78 * t95;
t98 = t107 * t61 + t121;
t84 = -pkin(4) * t90 - pkin(3);
t74 = t102 * t93;
t73 = (mrSges(4,1) * t93 + mrSges(4,2) * t96) * qJD(3);
t72 = -mrSges(5,1) * t96 - mrSges(5,3) * t115;
t71 = mrSges(5,2) * t96 - mrSges(5,3) * t119;
t66 = t102 * t107;
t65 = (mrSges(5,1) * t93 - mrSges(5,3) * t114) * qJD(3);
t64 = (-mrSges(5,2) * t93 - mrSges(5,3) * t118) * qJD(3);
t63 = t135 * t93;
t54 = t99 * t93;
t53 = t70 * t93;
t51 = t88 * t58;
t49 = -pkin(7) * t118 + t68;
t48 = (Ifges(5,5) * t93 + (Ifges(5,1) * t90 - t124) * t96) * qJD(3);
t47 = (Ifges(5,6) * t93 + (-t122 + t123) * t96) * qJD(3);
t46 = -mrSges(6,1) * t96 + mrSges(6,3) * t54;
t45 = mrSges(6,2) * t96 - mrSges(6,3) * t53;
t37 = -t105 * t90 + t51;
t35 = Ifges(6,1) * t70 - Ifges(6,4) * t99;
t34 = Ifges(6,4) * t70 - Ifges(6,2) * t99;
t33 = mrSges(6,1) * t99 + mrSges(6,2) * t70;
t31 = -Ifges(6,1) * t59 - Ifges(6,4) * t60;
t30 = -Ifges(6,4) * t59 - Ifges(6,2) * t60;
t29 = mrSges(6,1) * t60 - mrSges(6,2) * t59;
t28 = t51 + (-pkin(7) * t115 - pkin(8) * t118) * qJD(3);
t27 = mrSges(6,1) * t53 - mrSges(6,2) * t54;
t19 = (pkin(4) * t93 - pkin(8) * t114) * qJD(3) + t36;
t17 = -Ifges(6,1) * t54 - Ifges(6,4) * t53 - Ifges(6,5) * t96;
t16 = -Ifges(6,4) * t54 - Ifges(6,2) * t53 - Ifges(6,6) * t96;
t15 = -qJD(4) * t70 - qJD(5) * t44;
t14 = -qJD(4) * t99 + qJD(5) * t43;
t13 = -mrSges(6,2) * t108 + mrSges(6,3) * t26;
t12 = mrSges(6,1) * t108 - mrSges(6,3) * t25;
t6 = Ifges(6,1) * t25 + Ifges(6,4) * t26 + Ifges(6,5) * t108;
t5 = Ifges(6,4) * t25 + Ifges(6,2) * t26 + Ifges(6,6) * t108;
t4 = -qJD(5) * t9 + t19 * t95 - t28 * t92;
t3 = qJD(5) * t8 + t19 * t92 + t28 * t95;
t2 = -qJD(5) * t11 + t20 * t95 - t21 * t92;
t1 = qJD(5) * t10 + t20 * t92 + t21 * t95;
t22 = [0.2e1 * m(5) * (t20 * t38 + t21 * t39 + t18) + 0.2e1 * m(6) * (t1 * t11 + t10 * t2 + t18) + 0.2e1 * m(4) * (-t109 * t89 ^ 2 * t94 + t41 * t62 + t18); t1 * t45 + t10 * t12 + t11 * t13 + t2 * t46 + t20 * t72 + t21 * t71 + t38 * t65 + t39 * t64 + t127 * t61 + (t27 + t63) * t40 + (-t97 * t73 + (-t97 * mrSges(3,2) + (-t96 * mrSges(4,1) + t93 * mrSges(4,2) - mrSges(3,1)) * t94) * qJD(2)) * t89 + m(5) * (t20 * t49 + t21 * t50 + t36 * t38 + t37 * t39) + m(6) * (t1 * t9 + t10 * t4 + t11 * t3 + t2 * t8 + t40 * t74 + t61 * t66) - m(4) * pkin(2) * t104 + (t98 * t131 + m(4) * (-t62 * t108 + t120 + t98) / 0.2e1) * t132 + (t121 + t120 + (t61 * t96 - t62 * t93) * qJD(3)) * mrSges(4,3); -0.2e1 * pkin(2) * t73 + 0.2e1 * t8 * t12 + 0.2e1 * t9 * t13 + t26 * t16 + t25 * t17 + 0.2e1 * t66 * t27 + 0.2e1 * t3 * t45 + 0.2e1 * t36 * t72 + 0.2e1 * t37 * t71 + 0.2e1 * t4 * t46 + 0.2e1 * t49 * t65 - t53 * t5 + 0.2e1 * t50 * t64 - t54 * t6 + 0.2e1 * t74 * t7 + (t36 * t49 + t37 * t50) * t134 + (t3 * t9 + t4 * t8 + t66 * t74) * t133 + (t55 * t132 - t88 * t47 + t90 * t48 + (-Ifges(6,5) * t54 - Ifges(6,6) * t53 + (-(2 * Ifges(4,4)) - t101) * t93) * qJD(3)) * t93 + ((t63 * t132 + 0.2e1 * (Ifges(4,4) + t101) * t96 + (-(2 * Ifges(5,3)) + (2 * Ifges(4,1)) - (2 * Ifges(4,2)) - Ifges(6,3) + (pkin(7) ^ 2 * t134) + Ifges(5,1) * t87 + (t122 - 0.2e1 * t123) * t88) * t93) * qJD(3) + t106) * t96; -t41 * mrSges(4,2) + t61 * t29 + t100 * mrSges(5,3) + (t33 + t111) * t40 + m(5) * (-pkin(3) * t40 + (-t38 * t88 + t39 * t90) * qJD(4) + t100 * qJ(4)) + m(6) * (t1 * t44 + t10 * t15 + t11 * t14 + t2 * t43 + t40 * t84) + (-t1 * t99 + t10 * t59 - t11 * t60 - t2 * t70) * mrSges(6,3); -t96 * t112 / 0.2e1 + t84 * t7 + t6 * t129 + t74 * t29 - t59 * t17 / 0.2e1 - t60 * t16 / 0.2e1 + t66 * t33 + t5 * t130 - t53 * t30 / 0.2e1 - t54 * t31 / 0.2e1 - pkin(3) * t55 + t43 * t12 + t44 * t13 + t14 * t45 + t15 * t46 + t26 * t34 / 0.2e1 + t25 * t35 / 0.2e1 + m(6) * (t14 * t9 + t15 * t8 + t3 * t44 + t4 * t43 + t66 * t84) + (t47 / 0.2e1 + qJ(4) * t64 + qJD(4) * t71 + t37 * mrSges(5,3) + m(5) * (qJ(4) * t37 + qJD(4) * t50)) * t90 + (t48 / 0.2e1 - qJ(4) * t65 - qJD(4) * t72 - t36 * mrSges(5,3) + m(5) * (-qJ(4) * t36 - qJD(4) * t49)) * t88 + (-t3 * t99 - t4 * t70 + t59 * t8 - t60 * t9) * mrSges(6,3) + ((Ifges(5,5) * t88 / 0.2e1 + Ifges(5,6) * t128 + Ifges(6,5) * t129 + Ifges(6,6) * t130 - Ifges(4,6) + pkin(7) * mrSges(4,2)) * t93 + (Ifges(4,5) - t88 * (Ifges(5,2) * t90 + t124) / 0.2e1 + (Ifges(5,1) * t88 + t123) * t128 + (-m(5) * pkin(3) + t111) * pkin(7)) * t96) * qJD(3); (t14 * t44 + t15 * t43) * t133 - t59 * t35 + t70 * t31 + 0.2e1 * t84 * t29 - t60 * t34 - t99 * t30 + 0.2e1 * (-t14 * t99 - t15 * t70 + t43 * t59 - t44 * t60) * mrSges(6,3) + (qJ(4) * t134 + 0.2e1 * mrSges(5,3)) * qJD(4) * (t88 ^ 2 + t87); 0.2e1 * (t131 + m(6) / 0.2e1) * t40; m(5) * pkin(7) * t107 + m(6) * t66 + t127; t29; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1; mrSges(6,1) * t4 - mrSges(6,2) * t3 - t106; mrSges(6,1) * t15 - mrSges(6,2) * t14 + t112; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t22(1), t22(2), t22(4), t22(7), t22(11); t22(2), t22(3), t22(5), t22(8), t22(12); t22(4), t22(5), t22(6), t22(9), t22(13); t22(7), t22(8), t22(9), t22(10), t22(14); t22(11), t22(12), t22(13), t22(14), t22(15);];
Mq = res;
