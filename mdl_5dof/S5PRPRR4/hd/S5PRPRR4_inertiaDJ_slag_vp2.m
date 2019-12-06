% Calculate time derivative of joint inertia matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR4_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR4_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR4_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:49:43
% EndTime: 2019-12-05 15:49:48
% DurationCPUTime: 0.98s
% Computational Cost: add. (780->186), mult. (2233->303), div. (0->0), fcn. (1982->10), ass. (0->97)
t47 = sin(pkin(10));
t48 = sin(pkin(5));
t49 = cos(pkin(10));
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t21 = (t47 * t55 + t49 * t52) * t48;
t51 = sin(qJ(4));
t54 = cos(qJ(4));
t83 = cos(pkin(5));
t58 = -t21 * t51 + t54 * t83;
t81 = qJD(4) * t58;
t50 = sin(qJ(5));
t53 = cos(qJ(5));
t84 = t50 ^ 2 + t53 ^ 2;
t73 = qJD(5) * t53;
t78 = qJD(4) * t54;
t57 = -t50 * t78 - t51 * t73;
t74 = qJD(5) * t51;
t69 = t50 * t74;
t70 = t53 * t78;
t108 = t69 - t70;
t33 = -mrSges(6,1) * t53 + mrSges(6,2) * t50;
t107 = -m(6) * pkin(4) - mrSges(5,1) + t33;
t106 = 0.2e1 * m(6);
t41 = pkin(2) * t47 + pkin(7);
t105 = 0.2e1 * t41;
t104 = m(5) / 0.2e1;
t103 = m(6) / 0.2e1;
t94 = Ifges(6,4) * t50;
t34 = Ifges(6,2) * t53 + t94;
t102 = -t34 / 0.2e1;
t101 = -t50 / 0.2e1;
t100 = pkin(4) * t51;
t99 = pkin(8) * t54;
t60 = t47 * t52 - t49 * t55;
t82 = qJD(2) * t48;
t19 = t60 * t82;
t15 = t21 * t54 + t51 * t83;
t80 = qJD(4) * t15;
t5 = -t19 * t51 + t80;
t98 = t58 * t5;
t97 = t5 * t51;
t6 = -t19 * t54 + t81;
t96 = t6 * t54;
t95 = mrSges(6,3) * t51;
t93 = Ifges(6,4) * t53;
t92 = Ifges(6,5) * t50;
t91 = Ifges(6,6) * t50;
t90 = Ifges(6,6) * t53;
t89 = Ifges(6,6) * t54;
t18 = qJD(2) * t21;
t20 = t60 * t48;
t88 = t18 * t20;
t87 = t50 * t54;
t86 = t53 * t54;
t79 = qJD(4) * t51;
t85 = -Ifges(6,5) * t70 - Ifges(6,3) * t79;
t42 = -pkin(2) * t49 - pkin(3);
t25 = -pkin(4) * t54 - pkin(8) * t51 + t42;
t12 = t25 * t53 - t41 * t87;
t77 = qJD(5) * t12;
t13 = t25 * t50 + t41 * t86;
t76 = qJD(5) * t13;
t75 = qJD(5) * t50;
t72 = t41 * t79;
t67 = (2 * Ifges(5,4)) + t91;
t32 = (-t99 + t100) * qJD(4);
t3 = t32 * t50 - t53 * t72 + t77;
t66 = t3 - t77;
t8 = t15 * t53 + t20 * t50;
t1 = -qJD(5) * t8 + t18 * t53 - t50 * t6;
t7 = -t15 * t50 + t20 * t53;
t2 = qJD(5) * t7 + t18 * t50 + t53 * t6;
t65 = -t1 * t50 + t2 * t53;
t64 = t51 * mrSges(5,1) + t54 * mrSges(5,2);
t63 = mrSges(6,1) * t50 + mrSges(6,2) * t53;
t62 = Ifges(6,1) * t53 - t94;
t35 = Ifges(6,1) * t50 + t93;
t61 = -Ifges(6,2) * t50 + t93;
t59 = -t58 * t78 + t97;
t44 = Ifges(6,5) * t73;
t31 = -mrSges(6,1) * t54 - t53 * t95;
t30 = mrSges(6,2) * t54 - t50 * t95;
t29 = t62 * qJD(5);
t28 = t61 * qJD(5);
t27 = t64 * qJD(4);
t26 = t63 * qJD(5);
t24 = t63 * t51;
t23 = -Ifges(6,5) * t54 + t51 * t62;
t22 = t51 * t61 - t89;
t17 = -mrSges(6,2) * t79 + mrSges(6,3) * t57;
t16 = mrSges(6,1) * t79 + t108 * mrSges(6,3);
t11 = t57 * mrSges(6,1) + t108 * mrSges(6,2);
t10 = -t35 * t74 + (Ifges(6,5) * t51 + t54 * t62) * qJD(4);
t9 = -t34 * t74 + (Ifges(6,6) * t51 + t61 * t54) * qJD(4);
t4 = t32 * t53 + t50 * t72 - t76;
t14 = [0.2e1 * m(6) * (t1 * t7 + t2 * t8 - t98) + 0.2e1 * m(5) * (t15 * t6 + t88 - t98) + 0.2e1 * m(4) * (-t19 * t21 + t88); t19 * mrSges(4,2) + t1 * t31 + t58 * t11 + t7 * t16 + t8 * t17 + t2 * t30 + t20 * t27 + t5 * t24 + (-t54 * mrSges(5,1) + t51 * mrSges(5,2) - mrSges(4,1)) * t18 + (-mrSges(3,1) * t52 - mrSges(3,2) * t55) * t82 + m(6) * (t1 * t12 + t13 * t2 + t3 * t8 + t4 * t7) + m(5) * t42 * t18 + (t59 * t103 + (-t15 * t79 + t59 + t96) * t104) * t105 + m(4) * (-t18 * t49 - t19 * t47) * pkin(2) + (t97 + t96 + (-t15 * t51 - t54 * t58) * qJD(4)) * mrSges(5,3); (t12 * t4 + t13 * t3) * t106 + 0.2e1 * t3 * t30 + 0.2e1 * t13 * t17 + 0.2e1 * t4 * t31 + 0.2e1 * t12 * t16 + 0.2e1 * t42 * t27 + ((t105 * t24 - t50 * t22 + t53 * t23 + t54 * t67) * qJD(4) + t85) * t54 + (t53 * t10 - 0.2e1 * t41 * t11 - t50 * t9 + (-t50 * t23 - t53 * t22 - t54 * (-t90 - t92)) * qJD(5) + ((Ifges(6,5) * t53 - t67) * t51 + (t106 * t41 ^ 2 + (2 * Ifges(5,1)) - (2 * Ifges(5,2)) - Ifges(6,3)) * t54) * qJD(4)) * t51; 0.2e1 * ((-t5 + (-t50 * t7 + t53 * t8) * qJD(4)) * t103 + (-t5 + t80) * t104) * t54 + 0.2e1 * ((-t7 * t73 - t75 * t8 + t65 - t81) * t103 + (t6 - t81) * t104) * t51; t54 * t11 + (m(6) * (-t12 * t73 - t13 * t75 + t3 * t53 - t4 * t50) - t30 * t75 + t53 * t17 - t31 * t73 - t50 * t16) * t51 + (m(6) * (-t12 * t87 + t13 * t86 + (t51 ^ 2 - t54 ^ 2) * t41) + t30 * t86 - t31 * t87 + t51 * t24) * qJD(4); (-0.1e1 + t84) * t51 * t78 * t106; -t6 * mrSges(5,2) - t58 * t26 + (m(6) * pkin(8) + mrSges(6,3)) * ((-t50 * t8 - t53 * t7) * qJD(5) + t65) + t107 * t5; pkin(4) * t11 + (-t44 / 0.2e1 + (t107 * t41 + Ifges(5,5)) * qJD(4)) * t54 + (-t4 * mrSges(6,3) + t10 / 0.2e1 + t78 * t102 + (-t22 / 0.2e1 + t89 / 0.2e1 - t13 * mrSges(6,3)) * qJD(5) + (m(6) * (-t4 - t76) - qJD(5) * t30 - t16) * pkin(8)) * t50 + (qJD(5) * t23 / 0.2e1 + t9 / 0.2e1 + t35 * t78 / 0.2e1 + t66 * mrSges(6,3) + (m(6) * t66 - qJD(5) * t31 + t17) * pkin(8)) * t53 + (t28 * t101 + t41 * t26 + t53 * t29 / 0.2e1 + (t101 * t35 + t102 * t53) * qJD(5) + (-Ifges(5,6) + t92 / 0.2e1 + t90 / 0.2e1 + t41 * mrSges(5,2)) * qJD(4)) * t51; -t54 * t26 + (t51 * t33 + m(6) * (t84 * t99 - t100) + t84 * t54 * mrSges(6,3) - t64) * qJD(4); -0.2e1 * pkin(4) * t26 + t28 * t53 + t29 * t50 + (-t34 * t50 + t35 * t53) * qJD(5); mrSges(6,1) * t1 - mrSges(6,2) * t2; mrSges(6,1) * t4 - mrSges(6,2) * t3 - Ifges(6,5) * t69 + Ifges(6,6) * t57 - t85; t11; t44 + (pkin(8) * t33 - t91) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t14(1), t14(2), t14(4), t14(7), t14(11); t14(2), t14(3), t14(5), t14(8), t14(12); t14(4), t14(5), t14(6), t14(9), t14(13); t14(7), t14(8), t14(9), t14(10), t14(14); t14(11), t14(12), t14(13), t14(14), t14(15);];
Mq = res;
