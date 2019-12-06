% Calculate time derivative of joint inertia matrix for
% S5PRRPR1
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:35
% EndTime: 2019-12-05 16:15:36
% DurationCPUTime: 0.33s
% Computational Cost: add. (459->88), mult. (1078->129), div. (0->0), fcn. (828->6), ass. (0->47)
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t22 = t37 * t33 + t35 * t34;
t18 = t22 * qJD(5);
t21 = -t35 * t33 + t37 * t34;
t60 = t21 * t18;
t17 = t21 * qJD(5);
t59 = t22 * t17;
t58 = 2 * mrSges(5,3);
t57 = 2 * mrSges(6,3);
t38 = cos(qJ(3));
t47 = pkin(2) * qJD(3);
t44 = t38 * t47;
t27 = qJD(4) + t44;
t48 = t33 ^ 2 + t34 ^ 2;
t42 = t48 * t27;
t56 = -t34 * mrSges(5,1) - t21 * mrSges(6,1) + t33 * mrSges(5,2) + t22 * mrSges(6,2);
t39 = t48 * qJD(4);
t55 = 2 * m(5);
t54 = 2 * m(6);
t11 = t18 * mrSges(6,1) + t17 * mrSges(6,2);
t51 = 0.2e1 * t11;
t50 = t38 * pkin(2);
t49 = Ifges(6,5) * t17 - Ifges(6,6) * t18;
t36 = sin(qJ(3));
t45 = t36 * t47;
t43 = 0.2e1 * Ifges(6,1) * t59 - 0.2e1 * Ifges(6,2) * t60 + 0.2e1 * (t17 * t21 - t18 * t22) * Ifges(6,4);
t29 = -t34 * pkin(4) - pkin(3);
t41 = t48 * qJ(4);
t28 = t36 * pkin(2) + qJ(4);
t19 = (-pkin(7) - t28) * t33;
t30 = t34 * pkin(7);
t20 = t34 * t28 + t30;
t9 = t37 * t19 - t35 * t20;
t10 = t35 * t19 + t37 * t20;
t24 = (-pkin(7) - qJ(4)) * t33;
t26 = t34 * qJ(4) + t30;
t13 = t37 * t24 - t35 * t26;
t14 = t35 * t24 + t37 * t26;
t23 = t29 - t50;
t8 = -t22 * qJD(4) - t14 * qJD(5);
t7 = t21 * qJD(4) + t13 * qJD(5);
t2 = -t10 * qJD(5) - t22 * t27;
t1 = t9 * qJD(5) + t21 * t27;
t3 = [(t59 - t60) * t54; m(6) * (t1 * t22 + t10 * t17 - t9 * t18 + t2 * t21); -0.2e1 * mrSges(4,2) * t44 + t23 * t51 + (t10 * t1 + t9 * t2) * t54 + t43 - 0.2e1 * (t9 * t17 + t2 * t22) * mrSges(6,3) + (t1 * t21 - t10 * t18) * t57 + (t28 * t55 + t58) * t42 + (-(2 * mrSges(4,1)) + t23 * t54 + (-pkin(3) - t50) * t55 + 0.2e1 * t56) * t45; m(6) * (-t13 * t18 + t14 * t17 + t8 * t21 + t7 * t22); (t23 + t29) * t11 + (-mrSges(4,2) * t38 + (-mrSges(4,1) + t56) * t36) * t47 + m(6) * (t14 * t1 + t7 * t10 + t13 * t2 + t29 * t45 + t8 * t9) + m(5) * (-pkin(3) * t45 + t27 * t41 + t28 * t39) + (t42 + t39) * mrSges(5,3) + ((-t2 - t8) * t22 + (t1 + t7) * t21 - (t10 + t14) * t18 + (-t13 - t9) * t17) * mrSges(6,3) + t43; t29 * t51 + (t13 * t8 + t14 * t7) * t54 + t41 * t55 * qJD(4) + t43 + t39 * t58 + (-t13 * t17 - t14 * t18 + t7 * t21 - t8 * t22) * t57; 0; (m(5) + m(6)) * t45 + t11; t11; 0; -t11; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t49; t8 * mrSges(6,1) - t7 * mrSges(6,2) + t49; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t3(1), t3(2), t3(4), t3(7), t3(11); t3(2), t3(3), t3(5), t3(8), t3(12); t3(4), t3(5), t3(6), t3(9), t3(13); t3(7), t3(8), t3(9), t3(10), t3(14); t3(11), t3(12), t3(13), t3(14), t3(15);];
Mq = res;
