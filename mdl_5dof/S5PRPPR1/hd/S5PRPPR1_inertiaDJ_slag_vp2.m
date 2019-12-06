% Calculate time derivative of joint inertia matrix for
% S5PRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:21:48
% EndTime: 2019-12-05 15:21:50
% DurationCPUTime: 0.35s
% Computational Cost: add. (384->84), mult. (1003->151), div. (0->0), fcn. (890->6), ass. (0->44)
t31 = sin(pkin(8));
t28 = t31 ^ 2;
t30 = sin(pkin(9));
t32 = cos(pkin(9));
t34 = sin(qJ(5));
t35 = cos(qJ(5));
t36 = t34 * t30 - t35 * t32;
t18 = t36 * qJD(5);
t48 = 2 * m(6);
t23 = t35 * t30 + t34 * t32;
t14 = t23 * t31;
t47 = -0.2e1 * t14;
t15 = t36 * t31;
t46 = -0.2e1 * t15;
t45 = t30 * t31;
t44 = t31 * t32;
t19 = t23 * qJD(5);
t12 = t31 * t19;
t13 = t31 * t18;
t43 = -Ifges(6,5) * t12 + Ifges(6,6) * t13;
t33 = cos(pkin(8));
t25 = -t33 * pkin(3) - t31 * qJ(4) - pkin(2);
t41 = qJ(3) * t33;
t42 = t30 * t25 + t32 * t41;
t40 = qJD(3) * t31;
t39 = qJD(3) * t33;
t38 = qJD(4) * t31;
t37 = qJ(3) * qJD(3);
t21 = t32 * t25;
t6 = -pkin(6) * t44 + t21 + (-qJ(3) * t30 - pkin(4)) * t33;
t7 = -pkin(6) * t45 + t42;
t3 = -t34 * t7 + t35 * t6;
t4 = t34 * t6 + t35 * t7;
t5 = -t13 * mrSges(6,1) - t12 * mrSges(6,2);
t29 = t33 ^ 2;
t27 = t28 * t37;
t24 = (pkin(4) * t30 + qJ(3)) * t31;
t17 = -t30 * t38 + t32 * t39;
t16 = -t30 * t39 - t32 * t38;
t9 = -t33 * mrSges(6,1) + t15 * mrSges(6,3);
t8 = t33 * mrSges(6,2) - t14 * mrSges(6,3);
t2 = -t4 * qJD(5) + t35 * t16 - t34 * t17;
t1 = t3 * qJD(5) + t34 * t16 + t35 * t17;
t10 = [(t15 * t12 - t14 * t13) * t48; -t12 * t8 + t13 * t9 - t33 * t5 + (-t14 * t12 - t15 * t13) * mrSges(6,3) + m(5) * (-t16 * t30 + t17 * t32 - t39) * t31 + (-t1 * t15 - t4 * t12 + t3 * t13 - t2 * t14 - t39 * t31) * m(6); -t33 * t43 + 0.2e1 * t1 * t8 + 0.2e1 * t2 * t9 + 0.2e1 * t24 * t5 + 0.2e1 * t17 * (t33 * mrSges(5,2) - mrSges(5,3) * t45) + 0.2e1 * t16 * (-t33 * mrSges(5,1) - mrSges(5,3) * t44) + (0.2e1 * t4 * mrSges(6,3) + Ifges(6,4) * t46 + Ifges(6,2) * t47 - Ifges(6,6) * t33) * t13 - (-0.2e1 * t3 * mrSges(6,3) + Ifges(6,1) * t46 + Ifges(6,4) * t47 - Ifges(6,5) * t33) * t12 + ((0.2e1 * t14 * mrSges(6,1) + mrSges(6,2) * t46) * t31 + 0.2e1 * (mrSges(5,1) * t30 + mrSges(5,2) * t32) * t28 + 0.2e1 * (t29 + t28) * mrSges(4,3)) * qJD(3) + (t4 * t1 + t3 * t2 + t24 * t40) * t48 + 0.2e1 * m(5) * (t42 * t17 + (-t30 * t41 + t21) * t16 + t27) + 0.2e1 * m(4) * (t29 * t37 + t27); m(6) * (-t23 * t12 - t13 * t36 + t19 * t14 + t18 * t15); -t18 * t8 - t19 * t9 + (-t12 * t36 + t23 * t13) * mrSges(6,3) + m(6) * (t23 * t1 - t18 * t4 - t19 * t3 - t2 * t36) + m(5) * (t32 * t16 + t30 * t17); (-t23 * t18 + t19 * t36) * t48; 0; (m(5) + m(6)) * t40 + t5; 0; 0; -t5; t2 * mrSges(6,1) - t1 * mrSges(6,2) + t43; -t19 * mrSges(6,1) + t18 * mrSges(6,2); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t10(1), t10(2), t10(4), t10(7), t10(11); t10(2), t10(3), t10(5), t10(8), t10(12); t10(4), t10(5), t10(6), t10(9), t10(13); t10(7), t10(8), t10(9), t10(10), t10(14); t10(11), t10(12), t10(13), t10(14), t10(15);];
Mq = res;
