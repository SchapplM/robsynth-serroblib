% Calculate time derivative of joint inertia matrix for
% S5PRPRR6
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
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:29
% EndTime: 2019-12-05 15:56:33
% DurationCPUTime: 1.08s
% Computational Cost: add. (1162->186), mult. (3010->292), div. (0->0), fcn. (2898->10), ass. (0->94)
t55 = sin(pkin(10));
t57 = cos(pkin(10));
t112 = (m(4) * qJ(3) + mrSges(4,3)) * (t55 ^ 2 + t57 ^ 2);
t58 = cos(pkin(5));
t56 = sin(pkin(5));
t61 = sin(qJ(2));
t91 = t56 * t61;
t34 = -t55 * t91 + t57 * t58;
t35 = t55 * t58 + t57 * t91;
t110 = -t34 * t55 + t35 * t57;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t46 = -mrSges(6,1) * t62 + mrSges(6,2) * t59;
t109 = -m(6) * pkin(4) - mrSges(5,1) + t46;
t108 = -2 * mrSges(5,3);
t88 = pkin(7) + qJ(3);
t44 = t88 * t55;
t45 = t88 * t57;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t29 = -t44 * t60 + t45 * t63;
t40 = t55 * t63 + t57 * t60;
t19 = qJD(3) * t40 + qJD(4) * t29;
t107 = 0.2e1 * t19;
t69 = -t44 * t63 - t45 * t60;
t106 = -0.2e1 * t69;
t105 = -t40 / 0.2e1;
t102 = Ifges(6,4) * t59;
t47 = Ifges(6,2) * t62 + t102;
t104 = -t47 / 0.2e1;
t103 = mrSges(6,3) * t40;
t101 = Ifges(6,4) * t62;
t37 = t40 * qJD(4);
t100 = Ifges(6,5) * t37;
t99 = Ifges(6,6) * t37;
t39 = t55 * t60 - t57 * t63;
t98 = Ifges(6,6) * t39;
t97 = Ifges(6,6) * t59;
t96 = t19 * t69;
t21 = t34 * t60 + t35 * t63;
t64 = cos(qJ(2));
t82 = qJD(2) * t64;
t78 = t56 * t82;
t11 = qJD(4) * t21 + t40 * t78;
t70 = t34 * t63 - t35 * t60;
t95 = t70 * t11;
t92 = t56 ^ 2 * t61;
t90 = t56 * t64;
t36 = t39 * qJD(4);
t89 = t62 * t36;
t87 = -Ifges(6,5) * t89 + Ifges(6,3) * t37;
t50 = -pkin(3) * t57 - pkin(2);
t25 = pkin(4) * t39 - pkin(8) * t40 + t50;
t8 = t25 * t62 - t29 * t59;
t84 = qJD(5) * t8;
t9 = t25 * t59 + t29 * t62;
t83 = qJD(5) * t9;
t81 = qJD(5) * t59;
t80 = qJD(5) * t62;
t79 = qJD(2) * t91;
t77 = t40 * t81;
t22 = t37 * mrSges(5,1) - mrSges(5,2) * t36;
t75 = -(2 * Ifges(5,4)) - t97;
t74 = mrSges(6,1) * t59 + mrSges(6,2) * t62;
t73 = Ifges(6,1) * t62 - t102;
t72 = -Ifges(6,2) * t59 + t101;
t71 = -t11 * t69 - t19 * t70;
t14 = -t21 * t59 - t62 * t90;
t68 = -t21 * t62 + t59 * t90;
t67 = -t59 * t36 + t40 * t80;
t66 = t77 + t89;
t51 = Ifges(6,5) * t80;
t48 = Ifges(6,1) * t59 + t101;
t43 = t73 * qJD(5);
t42 = t72 * qJD(5);
t41 = t74 * qJD(5);
t27 = mrSges(6,1) * t39 - t103 * t62;
t26 = -mrSges(6,2) * t39 - t103 * t59;
t24 = pkin(4) * t37 + pkin(8) * t36;
t23 = t74 * t40;
t18 = -qJD(3) * t39 + qJD(4) * t69;
t17 = t39 * Ifges(6,5) + t40 * t73;
t16 = t40 * t72 + t98;
t13 = -mrSges(6,2) * t37 - mrSges(6,3) * t67;
t12 = mrSges(6,1) * t37 + mrSges(6,3) * t66;
t10 = qJD(4) * t70 - t39 * t78;
t7 = mrSges(6,1) * t67 - mrSges(6,2) * t66;
t6 = -Ifges(6,1) * t66 - Ifges(6,4) * t67 + t100;
t5 = -Ifges(6,4) * t66 - Ifges(6,2) * t67 + t99;
t4 = qJD(5) * t68 - t59 * t10 + t62 * t79;
t3 = qJD(5) * t14 + t62 * t10 + t59 * t79;
t2 = -t18 * t59 + t24 * t62 - t83;
t1 = t18 * t62 + t24 * t59 + t84;
t15 = [0.2e1 * m(6) * (t14 * t4 - t3 * t68 - t95) + 0.2e1 * m(4) * (t110 * t56 - t92) * t82 + 0.2e1 * (t10 * t21 - t82 * t92 - t95) * m(5); t11 * t23 + t14 * t12 - t68 * t13 - t70 * t7 + t3 * t26 + t4 * t27 + m(6) * (-t1 * t68 + t14 * t2 + t3 * t9 + t4 * t8 + t71) + m(5) * (t10 * t29 + t18 * t21 + t71) + m(4) * t110 * qJD(3) + ((-m(4) * pkin(2) + m(5) * t50 - mrSges(4,1) * t57 + mrSges(5,1) * t39 + mrSges(4,2) * t55 + mrSges(5,2) * t40 - mrSges(3,1)) * qJD(2) * t61 + (-t22 + (-mrSges(3,2) + t112) * qJD(2)) * t64) * t56 + (-t10 * t39 + t11 * t40 - t21 * t37 + t36 * t70) * mrSges(5,3); t29 * t37 * t108 + 0.2e1 * t1 * t26 + 0.2e1 * t8 * t12 + 0.2e1 * t9 * t13 + t23 * t107 + 0.2e1 * t2 * t27 + 0.2e1 * t50 * t22 + t7 * t106 + 0.2e1 * m(6) * (t1 * t9 + t2 * t8 - t96) + 0.2e1 * m(5) * (t18 * t29 - t96) - (mrSges(5,3) * t106 - t16 * t59 + t17 * t62) * t36 + 0.2e1 * t112 * qJD(3) + (t18 * t108 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t37 - t75 * t36 + t87) * t39 + (mrSges(5,3) * t107 - 0.2e1 * Ifges(5,1) * t36 - t59 * t5 + t62 * t6 + (Ifges(6,5) * t62 + t75) * t37 + (-t62 * t16 - t59 * t17 + t39 * (-Ifges(6,5) * t59 - Ifges(6,6) * t62)) * qJD(5)) * t40; m(6) * (t59 * t3 + t62 * t4 + (-t14 * t59 - t62 * t68) * qJD(5)) + (m(4) + m(5)) * t79; m(6) * (t1 * t59 + t2 * t62 + (-t59 * t8 + t62 * t9) * qJD(5)) + t26 * t80 + t59 * t13 - t27 * t81 + t62 * t12 + t22; 0; -t10 * mrSges(5,2) - t70 * t41 + (m(6) * pkin(8) + mrSges(6,3)) * (t3 * t62 - t4 * t59 + (-t14 * t62 + t59 * t68) * qJD(5)) + t109 * t11; t39 * t51 / 0.2e1 - t69 * t41 - Ifges(5,5) * t36 - Ifges(5,6) * t37 - t18 * mrSges(5,2) - pkin(4) * t7 + t109 * t19 + (t42 * t105 - t36 * t104 - t2 * mrSges(6,3) + t6 / 0.2e1 + t100 / 0.2e1 + (-t16 / 0.2e1 + t48 * t105 - t9 * mrSges(6,3) - t98 / 0.2e1) * qJD(5) + (m(6) * (-t2 - t83) - t12 - qJD(5) * t26) * pkin(8)) * t59 + (t40 * t43 / 0.2e1 - t36 * t48 / 0.2e1 + t1 * mrSges(6,3) + t99 / 0.2e1 + t5 / 0.2e1 + (t17 / 0.2e1 + t40 * t104 - t8 * mrSges(6,3)) * qJD(5) + (m(6) * (t1 - t84) + t13 - qJD(5) * t27) * pkin(8)) * t62; 0; -0.2e1 * pkin(4) * t41 + t42 * t62 + t43 * t59 + (-t59 * t47 + t62 * t48) * qJD(5); mrSges(6,1) * t4 - mrSges(6,2) * t3; mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,5) * t77 - Ifges(6,6) * t67 + t87; -t41; t51 + (pkin(8) * t46 - t97) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
