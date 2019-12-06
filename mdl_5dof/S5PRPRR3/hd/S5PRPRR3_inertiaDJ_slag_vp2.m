% Calculate time derivative of joint inertia matrix for
% S5PRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:46:39
% EndTime: 2019-12-05 15:46:40
% DurationCPUTime: 0.42s
% Computational Cost: add. (659->103), mult. (1640->180), div. (0->0), fcn. (1477->8), ass. (0->54)
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t30 = (t39 * mrSges(5,1) + t42 * mrSges(5,2)) * qJD(4);
t38 = sin(qJ(5));
t41 = cos(qJ(5));
t45 = t38 * t39 - t41 * t42;
t60 = qJD(4) + qJD(5);
t15 = t60 * t45;
t29 = t38 * t42 + t41 * t39;
t16 = t60 * t29;
t7 = t16 * mrSges(6,1) - t15 * mrSges(6,2);
t61 = t7 + t30;
t36 = sin(pkin(9));
t37 = cos(pkin(9));
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t27 = t36 * t43 + t37 * t40;
t9 = t45 * t27;
t59 = 2 * m(6);
t34 = t39 ^ 2;
t58 = m(6) * pkin(4);
t32 = t36 * pkin(2) + pkin(6);
t57 = pkin(7) + t32;
t21 = t27 * qJD(2);
t26 = t36 * t40 - t37 * t43;
t17 = t26 * t21;
t56 = t45 * t16;
t55 = t29 * t15;
t54 = t42 ^ 2 + t34;
t53 = qJD(4) * t39;
t52 = qJD(4) * t42;
t51 = pkin(4) * t53;
t22 = t26 * qJD(2);
t2 = -t16 * t27 + t22 * t45;
t3 = t29 * t22 + t60 * t9;
t50 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t33 = -t37 * pkin(2) - pkin(3);
t49 = t22 * t54;
t48 = qJD(4) * t57;
t47 = -t42 * mrSges(5,1) + t39 * mrSges(5,2);
t23 = t57 * t39;
t24 = t57 * t42;
t10 = -t41 * t23 - t38 * t24;
t11 = -t38 * t23 + t41 * t24;
t19 = t39 * t48;
t20 = t42 * t48;
t5 = t10 * qJD(5) - t41 * t19 - t38 * t20;
t6 = -t11 * qJD(5) + t38 * t19 - t41 * t20;
t44 = t6 * mrSges(6,1) - t5 * mrSges(6,2) - Ifges(6,5) * t15 - Ifges(6,6) * t16;
t31 = -t42 * pkin(4) + t33;
t25 = (-mrSges(6,1) * t38 - mrSges(6,2) * t41) * qJD(5) * pkin(4);
t18 = mrSges(6,1) * t45 + t29 * mrSges(6,2);
t8 = t29 * t27;
t1 = [0.2e1 * m(6) * (-t9 * t2 - t8 * t3 + t17) + 0.2e1 * m(5) * (-t27 * t49 + t17) + 0.2e1 * m(4) * (-t27 * t22 + t17); t61 * t26 + (-t40 * mrSges(3,1) - t43 * mrSges(3,2)) * qJD(2) - (mrSges(5,3) * t54 - mrSges(4,2)) * t22 + (-mrSges(4,1) + t18 + t47) * t21 + m(6) * (t10 * t3 + t11 * t2 + t31 * t21 + t26 * t51 - t5 * t9 - t6 * t8) + m(5) * (t33 * t21 - t32 * t49) + m(4) * (-t21 * t37 - t22 * t36) * pkin(2) + (-t8 * t15 + t9 * t16 - t2 * t45 - t3 * t29) * mrSges(6,3); -0.2e1 * Ifges(6,1) * t55 + 0.2e1 * Ifges(6,2) * t56 + (t10 * t6 + t11 * t5 + t31 * t51) * t59 + 0.2e1 * t18 * t51 + 0.2e1 * t31 * t7 + 0.2e1 * t33 * t30 + 0.2e1 * (-t34 * qJD(4) + t42 * t52) * Ifges(5,4) + 0.2e1 * (t15 * t45 - t29 * t16) * Ifges(6,4) + 0.2e1 * (t10 * t15 - t11 * t16 - t6 * t29 - t45 * t5) * mrSges(6,3) + 0.2e1 * (Ifges(5,1) - Ifges(5,2)) * t39 * t52; m(6) * (t15 * t9 + t16 * t8 + t29 * t2 - t3 * t45); m(6) * (-t16 * t10 - t15 * t11 + t29 * t5 - t45 * t6); (-t55 + t56) * t59; (t42 * t22 + t27 * t53) * mrSges(5,2) + (t39 * t22 - t27 * t52) * mrSges(5,1) + (t2 * t38 + t3 * t41 + (t38 * t8 - t41 * t9) * qJD(5)) * t58 + t50; (Ifges(5,5) * t42 - Ifges(5,6) * t39 + t32 * t47) * qJD(4) + (m(6) * (t38 * t5 + t41 * t6 + (-t10 * t38 + t11 * t41) * qJD(5)) + (t41 * t15 - t38 * t16 + (t29 * t38 - t41 * t45) * qJD(5)) * mrSges(6,3)) * pkin(4) + t44; (-t15 * t38 - t16 * t41 + (t29 * t41 + t38 * t45) * qJD(5)) * t58 - t61; 0.2e1 * t25; t50; t44; -t7; t25; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
