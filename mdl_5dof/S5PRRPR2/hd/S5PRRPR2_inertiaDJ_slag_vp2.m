% Calculate time derivative of joint inertia matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:17:19
% EndTime: 2019-12-05 16:17:21
% DurationCPUTime: 0.49s
% Computational Cost: add. (444->117), mult. (1084->173), div. (0->0), fcn. (734->6), ass. (0->65)
t37 = cos(qJ(3));
t58 = pkin(2) * qJD(3);
t25 = t37 * t58 + qJD(4);
t35 = sin(qJ(3));
t27 = t35 * pkin(2) + qJ(4);
t79 = qJ(4) * t25 + qJD(4) * t27;
t32 = sin(pkin(9));
t28 = t32 ^ 2;
t33 = cos(pkin(9));
t29 = t33 ^ 2;
t78 = (t28 + t29) * mrSges(5,3);
t77 = qJD(4) + t25;
t76 = ((-t33 * mrSges(5,1) + t32 * mrSges(5,2) - mrSges(4,1)) * t35 - mrSges(4,2) * t37) * t58;
t75 = 2 * m(5);
t74 = 2 * m(6);
t34 = sin(qJ(5));
t36 = cos(qJ(5));
t53 = qJD(5) * t32;
t12 = (mrSges(6,1) * t36 - mrSges(6,2) * t34) * t53;
t73 = 0.2e1 * t12;
t42 = mrSges(6,1) * t34 + mrSges(6,2) * t36;
t17 = t42 * t32;
t72 = 0.2e1 * t17;
t67 = mrSges(6,3) * t32;
t19 = t33 * mrSges(6,2) - t34 * t67;
t71 = 0.2e1 * t19;
t20 = -t33 * mrSges(6,1) - t36 * t67;
t70 = 0.2e1 * t20;
t69 = t37 * pkin(2);
t66 = Ifges(6,4) * t34;
t65 = Ifges(6,4) * t36;
t64 = t27 * t25;
t63 = t33 * t34;
t62 = t33 * t36;
t61 = t34 * (-Ifges(6,2) * t36 - t66) * t53;
t60 = t79 * t28;
t56 = qJ(4) * t33;
t54 = qJD(4) * t33;
t52 = qJD(5) * t34;
t51 = qJD(5) * t36;
t50 = qJ(4) * qJD(4);
t49 = 0.2e1 * mrSges(6,3);
t48 = t35 * t58;
t13 = (-Ifges(6,5) * t34 - Ifges(6,6) * t36) * t53;
t47 = -t33 * t13 + t28 * (-Ifges(6,1) * t34 - t65) * t51;
t45 = t19 * t51 - t20 * t52;
t23 = -t33 * pkin(4) - t32 * pkin(7) - pkin(3);
t18 = t23 - t69;
t5 = t36 * t18 - t27 * t63;
t6 = t34 * t18 + t27 * t62;
t44 = t34 * t5 - t36 * t6;
t8 = t36 * t23 - t34 * t56;
t9 = t34 * t23 + t36 * t56;
t43 = t34 * t8 - t36 * t9;
t41 = -(-Ifges(6,6) * t33 + (-Ifges(6,2) * t34 + t65) * t32) * t36 - (-Ifges(6,5) * t33 + (Ifges(6,1) * t36 - t66) * t32) * t34;
t40 = 0.2e1 * t78;
t39 = -t19 * t52 - t20 * t51;
t38 = -t33 * t12 + (-t34 ^ 2 - t36 ^ 2) * t28 * qJD(5) * mrSges(6,3);
t26 = t28 * t50;
t15 = t28 * t64;
t4 = -qJD(5) * t9 - t34 * t54;
t3 = qJD(5) * t8 + t36 * t54;
t2 = -qJD(5) * t6 - t25 * t63 + t36 * t48;
t1 = qJD(5) * t5 + t25 * t62 + t34 * t48;
t7 = [0; (m(6) * (t1 * t36 - t2 * t34 - t25 * t33 - t5 * t51 - t6 * t52) + t39) * t32 + t38; t1 * t71 + t2 * t70 + t25 * t40 + 0.2e1 * t76 + (t6 * t1 + t5 * t2 + t15) * t74 + (t29 * t64 + t15 + (-pkin(3) - t69) * t48) * t75 + (t27 * t73 - t61 + t25 * t72 + (t44 * t49 + t41) * qJD(5)) * t32 + t47; (m(6) * (t3 * t36 - t34 * t4 - t8 * t51 - t9 * t52 - t54) + t39) * t32 + t38; (t4 + t2) * t20 + (t3 + t1) * t19 + t76 + m(5) * (-pkin(3) * t48 + t79 * t29 + t60) + m(6) * (t9 * t1 + t8 * t2 + t3 * t6 + t4 * t5 + t60) + (-t61 + t77 * t17 + (qJ(4) + t27) * t12 + (((-t6 - t9) * t36 + (t5 + t8) * t34) * mrSges(6,3) + t41) * qJD(5)) * t32 + t47 + t77 * t78; t3 * t71 + t4 * t70 + qJD(4) * t40 + (t9 * t3 + t8 * t4 + t26) * t74 + (t29 * t50 + t26) * t75 + (qJ(4) * t73 + qJD(4) * t72 - t61 + (t43 * t49 + t41) * qJD(5)) * t32 + t47; 0; m(6) * (-qJD(5) * t44 + t34 * t1 + t36 * t2) + m(5) * t48 + t45; m(6) * (-qJD(5) * t43 + t34 * t3 + t36 * t4) + t45; 0; -t12; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t13; t4 * mrSges(6,1) - t3 * mrSges(6,2) + t13; -t42 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t7(1), t7(2), t7(4), t7(7), t7(11); t7(2), t7(3), t7(5), t7(8), t7(12); t7(4), t7(5), t7(6), t7(9), t7(13); t7(7), t7(8), t7(9), t7(10), t7(14); t7(11), t7(12), t7(13), t7(14), t7(15);];
Mq = res;
