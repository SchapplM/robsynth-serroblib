% Calculate time derivative of joint inertia matrix for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2020-01-03 11:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:20:06
% EndTime: 2020-01-03 11:20:09
% DurationCPUTime: 0.36s
% Computational Cost: add. (465->86), mult. (1084->153), div. (0->0), fcn. (971->8), ass. (0->45)
t32 = sin(pkin(8));
t29 = t32 ^ 2;
t31 = sin(pkin(9));
t33 = cos(pkin(9));
t35 = sin(qJ(5));
t36 = cos(qJ(5));
t37 = t35 * t31 - t33 * t36;
t21 = t37 * qJD(5);
t49 = 2 * m(6);
t25 = t31 * t36 + t35 * t33;
t14 = t25 * t32;
t48 = -0.2e1 * t14;
t15 = t37 * t32;
t47 = -0.2e1 * t15;
t28 = sin(pkin(7)) * pkin(1) + qJ(3);
t34 = cos(pkin(8));
t46 = t28 * t34;
t45 = t31 * t32;
t44 = t32 * t33;
t22 = t25 * qJD(5);
t12 = t32 * t22;
t13 = t32 * t21;
t43 = -Ifges(6,5) * t12 + Ifges(6,6) * t13;
t23 = -cos(pkin(7)) * pkin(1) - pkin(2) - t32 * qJ(4) - t34 * pkin(3);
t42 = t31 * t23 + t33 * t46;
t41 = qJD(3) * t28;
t40 = qJD(3) * t32;
t39 = qJD(3) * t34;
t38 = qJD(4) * t32;
t17 = t33 * t23;
t6 = -pkin(6) * t44 + t17 + (-t28 * t31 - pkin(4)) * t34;
t7 = -pkin(6) * t45 + t42;
t3 = -t35 * t7 + t36 * t6;
t4 = t35 * t6 + t36 * t7;
t5 = -t13 * mrSges(6,1) - t12 * mrSges(6,2);
t30 = t34 ^ 2;
t27 = t29 * t41;
t20 = -t31 * t38 + t33 * t39;
t19 = -t31 * t39 - t33 * t38;
t18 = (pkin(4) * t31 + t28) * t32;
t9 = -mrSges(6,1) * t34 + t15 * mrSges(6,3);
t8 = mrSges(6,2) * t34 - t14 * mrSges(6,3);
t2 = -t4 * qJD(5) + t36 * t19 - t35 * t20;
t1 = t3 * qJD(5) + t35 * t19 + t36 * t20;
t10 = [0.2e1 * t20 * (t34 * mrSges(5,2) - mrSges(5,3) * t45) + 0.2e1 * t19 * (-t34 * mrSges(5,1) - mrSges(5,3) * t44) + 0.2e1 * t18 * t5 - t34 * t43 + 0.2e1 * t1 * t8 + 0.2e1 * t2 * t9 + (0.2e1 * mrSges(6,3) * t4 + Ifges(6,4) * t47 + Ifges(6,2) * t48 - Ifges(6,6) * t34) * t13 - (-0.2e1 * mrSges(6,3) * t3 + Ifges(6,1) * t47 + Ifges(6,4) * t48 - Ifges(6,5) * t34) * t12 + ((0.2e1 * t14 * mrSges(6,1) + mrSges(6,2) * t47) * t32 + 0.2e1 * (mrSges(5,1) * t31 + mrSges(5,2) * t33) * t29 + 0.2e1 * (t30 + t29) * mrSges(4,3)) * qJD(3) + 0.2e1 * m(4) * (t30 * t41 + t27) + 0.2e1 * m(5) * (t42 * t20 + (-t31 * t46 + t17) * t19 + t27) + (t1 * t4 + t18 * t40 + t2 * t3) * t49; -t12 * t8 + t13 * t9 - t34 * t5 + (-t12 * t14 - t13 * t15) * mrSges(6,3) + m(5) * (-t19 * t31 + t20 * t33 - t39) * t32 + (-t1 * t15 - t12 * t4 + t13 * t3 - t14 * t2 - t39 * t32) * m(6); (t12 * t15 - t13 * t14) * t49; -t21 * t8 - t22 * t9 + (-t12 * t37 + t13 * t25) * mrSges(6,3) + m(6) * (t1 * t25 - t2 * t37 - t21 * t4 - t22 * t3) + m(5) * (t19 * t33 + t20 * t31); m(6) * (-t12 * t25 - t13 * t37 + t14 * t22 + t15 * t21); (-t21 * t25 + t22 * t37) * t49; (m(5) + m(6)) * t40 + t5; 0; 0; 0; mrSges(6,1) * t2 - mrSges(6,2) * t1 + t43; -t5; -mrSges(6,1) * t22 + mrSges(6,2) * t21; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
