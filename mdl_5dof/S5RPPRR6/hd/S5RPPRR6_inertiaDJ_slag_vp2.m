% Calculate time derivative of joint inertia matrix for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR6_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR6_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR6_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:57:41
% EndTime: 2019-12-31 17:57:43
% DurationCPUTime: 0.65s
% Computational Cost: add. (1053->153), mult. (2292->238), div. (0->0), fcn. (2062->8), ass. (0->74)
t90 = -2 * mrSges(5,3);
t47 = sin(pkin(9));
t48 = cos(pkin(9));
t51 = sin(qJ(4));
t53 = cos(qJ(4));
t32 = t47 * t51 - t53 * t48;
t28 = t32 * qJD(4);
t33 = t47 * t53 + t51 * t48;
t50 = sin(qJ(5));
t52 = cos(qJ(5));
t67 = qJD(5) * t52;
t55 = t50 * t28 - t33 * t67;
t68 = qJD(5) * t50;
t66 = t33 * t68;
t74 = t52 * t28;
t54 = t66 + t74;
t40 = sin(pkin(8)) * pkin(1) + qJ(3);
t85 = pkin(6) + t40;
t30 = t85 * t47;
t31 = t85 * t48;
t16 = -t51 * t30 + t31 * t53;
t11 = t33 * qJD(3) + t16 * qJD(4);
t89 = 0.2e1 * t11;
t88 = -t33 / 0.2e1;
t83 = Ifges(6,4) * t50;
t38 = Ifges(6,2) * t52 + t83;
t87 = -t38 / 0.2e1;
t29 = t33 * qJD(4);
t86 = pkin(4) * t29;
t84 = mrSges(6,3) * t33;
t82 = Ifges(6,4) * t52;
t81 = Ifges(6,5) * t29;
t80 = Ifges(6,6) * t29;
t79 = Ifges(6,6) * t32;
t78 = Ifges(6,6) * t50;
t57 = -t53 * t30 - t31 * t51;
t77 = t11 * t57;
t76 = t32 * t29;
t73 = -Ifges(6,5) * t74 + Ifges(6,3) * t29;
t56 = -cos(pkin(8)) * pkin(1) - pkin(3) * t48 - pkin(2);
t14 = pkin(4) * t32 - pkin(7) * t33 + t56;
t5 = t14 * t52 - t16 * t50;
t71 = qJD(5) * t5;
t6 = t14 * t50 + t16 * t52;
t70 = qJD(5) * t6;
t69 = qJD(5) * t33;
t64 = (t50 ^ 2 + t52 ^ 2) * t28;
t63 = t29 * mrSges(5,1) - t28 * mrSges(5,2);
t62 = -(2 * Ifges(5,4)) - t78;
t37 = -mrSges(6,1) * t52 + mrSges(6,2) * t50;
t61 = mrSges(6,1) * t50 + mrSges(6,2) * t52;
t60 = Ifges(6,1) * t52 - t83;
t59 = -Ifges(6,2) * t50 + t82;
t58 = t11 * t32 - t29 * t57;
t42 = Ifges(6,5) * t67;
t39 = Ifges(6,1) * t50 + t82;
t36 = t60 * qJD(5);
t35 = t59 * qJD(5);
t34 = t61 * qJD(5);
t20 = mrSges(6,1) * t32 - t52 * t84;
t19 = -mrSges(6,2) * t32 - t50 * t84;
t18 = pkin(7) * t28 + t86;
t17 = t61 * t33;
t13 = Ifges(6,5) * t32 + t60 * t33;
t12 = t59 * t33 + t79;
t10 = -t32 * qJD(3) + t57 * qJD(4);
t9 = -t29 * mrSges(6,2) + t55 * mrSges(6,3);
t8 = t29 * mrSges(6,1) + t54 * mrSges(6,3);
t7 = mrSges(6,1) * t55 + mrSges(6,2) * t54;
t4 = -t54 * Ifges(6,1) + t55 * Ifges(6,4) + t81;
t3 = -t54 * Ifges(6,4) + t55 * Ifges(6,2) + t80;
t2 = -t50 * t10 + t52 * t18 - t70;
t1 = t52 * t10 + t50 * t18 + t71;
t15 = [t16 * t29 * t90 + 0.2e1 * t56 * t63 + 0.2e1 * t1 * t19 + 0.2e1 * t2 * t20 + 0.2e1 * t5 * t8 + 0.2e1 * t6 * t9 + 0.2e1 * t57 * t7 + t17 * t89 + 0.2e1 * m(6) * (t1 * t6 + t2 * t5 - t77) + 0.2e1 * m(5) * (t10 * t16 - t77) - (-t50 * t12 + t13 * t52 + t57 * t90) * t28 + (t10 * t90 + ((2 * Ifges(5,2)) + Ifges(6,3)) * t29 - t62 * t28 + t73) * t32 + (mrSges(5,3) * t89 - 0.2e1 * Ifges(5,1) * t28 - t50 * t3 + t52 * t4 + (Ifges(6,5) * t52 + t62) * t29 + (t32 * (-Ifges(6,5) * t50 - Ifges(6,6) * t52) - t50 * t13 - t52 * t12) * qJD(5)) * t33 + 0.2e1 * (m(4) * t40 + mrSges(4,3)) * qJD(3) * (t47 ^ 2 + t48 ^ 2); t29 * t17 - t32 * t7 + m(5) * (t10 * t33 - t16 * t28 + t58) + m(6) * t58 + (-t28 * t19 + t33 * t9 - t20 * t69 + m(6) * (t1 * t33 - t28 * t6 - t5 * t69)) * t52 + (-t19 * t69 + t28 * t20 - t33 * t8 + m(6) * (-t2 * t33 + t28 * t5 - t6 * t69)) * t50; 0.2e1 * m(5) * (-t28 * t33 + t76) + 0.2e1 * m(6) * (-t33 * t64 + t76); m(6) * (t50 * t1 + t52 * t2 + (-t5 * t50 + t52 * t6) * qJD(5)) + t19 * t67 + t50 * t9 - t20 * t68 + t52 * t8 + t63; 0; 0; t32 * t42 / 0.2e1 - t57 * t34 - Ifges(5,6) * t29 - Ifges(5,5) * t28 - t10 * mrSges(5,2) + pkin(4) * t7 + (-m(6) * pkin(4) - mrSges(5,1) + t37) * t11 + (t35 * t88 - t28 * t87 - t2 * mrSges(6,3) + t4 / 0.2e1 + t81 / 0.2e1 + (-t12 / 0.2e1 - t79 / 0.2e1 + t39 * t88 - t6 * mrSges(6,3)) * qJD(5) + (-qJD(5) * t19 - t8 + m(6) * (-t2 - t70)) * pkin(7)) * t50 + (t33 * t36 / 0.2e1 - t28 * t39 / 0.2e1 + t1 * mrSges(6,3) + t80 / 0.2e1 + t3 / 0.2e1 + (t13 / 0.2e1 + t33 * t87 - t5 * mrSges(6,3)) * qJD(5) + (m(6) * (t1 - t71) + t9 - qJD(5) * t20) * pkin(7)) * t52; m(6) * (-pkin(7) * t64 - t86) + t29 * t37 + t32 * t34 - mrSges(6,3) * t64 - t63; 0; -0.2e1 * pkin(4) * t34 + t35 * t52 + t36 * t50 + (-t50 * t38 + t52 * t39) * qJD(5); mrSges(6,1) * t2 - mrSges(6,2) * t1 - Ifges(6,5) * t66 + t55 * Ifges(6,6) + t73; t7; -t34; t42 + (t37 * pkin(7) - t78) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
