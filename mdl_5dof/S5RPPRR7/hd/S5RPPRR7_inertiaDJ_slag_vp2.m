% Calculate time derivative of joint inertia matrix for
% S5RPPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 18:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR7_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR7_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR7_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR7_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:59:31
% EndTime: 2019-12-31 17:59:33
% DurationCPUTime: 0.71s
% Computational Cost: add. (520->144), mult. (1147->232), div. (0->0), fcn. (758->6), ass. (0->73)
t77 = 2 * qJD(3);
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t59 = t34 ^ 2 + t36 ^ 2;
t82 = m(6) * (-0.1e1 + t59);
t21 = -t36 * mrSges(6,1) + t34 * mrSges(6,2);
t80 = -mrSges(5,1) + t21;
t35 = sin(qJ(4));
t37 = cos(qJ(4));
t79 = t35 ^ 2 - t37 ^ 2;
t25 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(6);
t78 = -0.2e1 * t25;
t54 = qJD(4) * t37;
t76 = t35 * t54 * t82;
t75 = -t34 / 0.2e1;
t74 = pkin(4) * t37;
t73 = pkin(7) * t35;
t72 = t35 * pkin(4);
t71 = t37 * pkin(7);
t70 = Ifges(6,4) * t34;
t69 = Ifges(6,4) * t36;
t68 = Ifges(6,5) * t34;
t67 = Ifges(6,6) * t34;
t66 = Ifges(6,6) * t35;
t65 = Ifges(6,6) * t36;
t64 = t34 * t35;
t63 = t35 * t36;
t61 = t37 * mrSges(6,3);
t19 = -t35 * mrSges(6,2) - t34 * t61;
t62 = t36 * t19;
t55 = qJD(4) * t35;
t49 = t34 * t55;
t60 = Ifges(6,6) * t49 + Ifges(6,3) * t54;
t26 = sin(pkin(8)) * pkin(1) + qJ(3);
t13 = t26 - t71 + t72;
t6 = t36 * t13 - t25 * t64;
t57 = qJD(5) * t6;
t7 = t34 * t13 + t25 * t63;
t56 = qJD(5) * t7;
t53 = qJD(5) * t34;
t52 = qJD(5) * t36;
t51 = qJD(5) * t37;
t50 = t36 * t51;
t48 = t34 * t54;
t47 = t36 * t54;
t46 = -Ifges(6,5) * t36 + (2 * Ifges(5,4));
t15 = qJD(3) + (t73 + t74) * qJD(4);
t1 = t34 * t15 + t25 * t47 + t57;
t45 = t1 - t57;
t44 = mrSges(6,1) * t34 + mrSges(6,2) * t36;
t43 = Ifges(6,1) * t36 - t70;
t23 = Ifges(6,1) * t34 + t69;
t42 = -Ifges(6,2) * t34 + t69;
t22 = Ifges(6,2) * t36 + t70;
t41 = t34 * t51 + t36 * t55;
t40 = t49 - t50;
t2 = t36 * t15 - t25 * t48 - t56;
t39 = t1 * t36 - t2 * t34 - t6 * t52 - t7 * t53;
t10 = -mrSges(6,2) * t54 + t40 * mrSges(6,3);
t14 = t44 * t37;
t20 = t35 * mrSges(6,1) - t36 * t61;
t9 = mrSges(6,1) * t54 + t41 * mrSges(6,3);
t38 = qJD(4) * t14 + t36 * t10 - t19 * t53 - t20 * t52 - t34 * t9;
t28 = Ifges(6,5) * t52;
t18 = t43 * qJD(5);
t17 = t42 * qJD(5);
t16 = t44 * qJD(5);
t12 = Ifges(6,5) * t35 + t43 * t37;
t11 = t42 * t37 + t66;
t5 = -t40 * mrSges(6,1) - t41 * mrSges(6,2);
t4 = -t23 * t51 + (Ifges(6,5) * t37 - t43 * t35) * qJD(4);
t3 = -t22 * t51 + (Ifges(6,6) * t37 - t42 * t35) * qJD(4);
t8 = [0.2e1 * t1 * t19 + 0.2e1 * t7 * t10 + 0.2e1 * t2 * t20 + 0.2e1 * t6 * t9 + 0.2e1 * m(6) * (t7 * t1 + t6 * t2) + (mrSges(4,3) + (m(4) + m(5)) * t26) * t77 + (mrSges(5,1) * t77 + (-0.2e1 * t26 * mrSges(5,2) + t34 * t11 - t36 * t12 + 0.2e1 * t25 * t14 + t46 * t35) * qJD(4) + t60) * t35 + (mrSges(5,2) * t77 + t5 * t78 - t34 * t3 + t36 * t4 + (t35 * (-t65 - t68) - t34 * t12 - t36 * t11) * qJD(5) + (0.2e1 * t26 * mrSges(5,1) + (-t46 - t67) * t37 + (-0.2e1 * m(6) * t25 ^ 2 - (2 * Ifges(5,1)) + (2 * Ifges(5,2)) + Ifges(6,3)) * t35) * qJD(4)) * t37; t35 * t5 + (-t35 * t62 + t20 * t64 + m(6) * (t25 * t79 + t6 * t64 - t7 * t63)) * qJD(4) + (m(6) * t39 + t38) * t37; -0.2e1 * t76; (-t5 + (m(6) * (-t34 * t6 + t36 * t7) + t62 - t34 * t20) * qJD(4)) * t37 + (m(6) * (t54 * t78 + t39) + t38) * t35; -t79 * qJD(4) * t82; 0.2e1 * t76; -pkin(4) * t5 + (t28 / 0.2e1 + (-Ifges(5,5) + (-m(6) * pkin(4) + t80) * t25) * qJD(4)) * t35 + (t4 / 0.2e1 - t2 * mrSges(6,3) + t22 * t55 / 0.2e1 + (-t66 / 0.2e1 - t11 / 0.2e1 - t7 * mrSges(6,3)) * qJD(5) + (m(6) * (-t2 - t56) - qJD(5) * t19 - t9) * pkin(7)) * t34 + (t3 / 0.2e1 + qJD(5) * t12 / 0.2e1 - t23 * t55 / 0.2e1 + t45 * mrSges(6,3) + (m(6) * t45 - qJD(5) * t20 + t10) * pkin(7)) * t36 + (t36 * t18 / 0.2e1 + t17 * t75 - t25 * t16 + (-t36 * t22 / 0.2e1 + t23 * t75) * qJD(5) + (t68 / 0.2e1 + t65 / 0.2e1 - Ifges(5,6) - t25 * mrSges(5,2)) * qJD(4)) * t37; (m(6) * (-t59 * t73 - t74) + t80 * t37) * qJD(4) + (t16 + (-t59 * mrSges(6,3) + mrSges(5,2)) * qJD(4)) * t35; -t37 * t16 + (-t37 * mrSges(5,2) + m(6) * (t59 * t71 - t72) + t59 * t61 + t80 * t35) * qJD(4); -0.2e1 * pkin(4) * t16 + t36 * t17 + t34 * t18 + (-t22 * t34 + t23 * t36) * qJD(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) - t41 * Ifges(6,5) - Ifges(6,6) * t50 + t60; -t5; (t35 * t53 - t47) * mrSges(6,2) + (-t35 * t52 - t48) * mrSges(6,1); t28 + (t21 * pkin(7) - t67) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t8(1), t8(2), t8(4), t8(7), t8(11); t8(2), t8(3), t8(5), t8(8), t8(12); t8(4), t8(5), t8(6), t8(9), t8(13); t8(7), t8(8), t8(9), t8(10), t8(14); t8(11), t8(12), t8(13), t8(14), t8(15);];
Mq = res;
