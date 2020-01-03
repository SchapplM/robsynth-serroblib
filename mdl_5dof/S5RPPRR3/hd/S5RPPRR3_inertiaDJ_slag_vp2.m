% Calculate time derivative of joint inertia matrix for
% S5RPPRR3
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:30
% EndTime: 2020-01-03 11:27:32
% DurationCPUTime: 0.43s
% Computational Cost: add. (1003->93), mult. (2096->155), div. (0->0), fcn. (1988->8), ass. (0->47)
t39 = sin(pkin(9));
t36 = sin(pkin(8)) * pkin(1) + qJ(3);
t56 = pkin(6) + t36;
t30 = t56 * t39;
t40 = cos(pkin(9));
t31 = t56 * t40;
t43 = sin(qJ(4));
t45 = cos(qJ(4));
t19 = -t43 * t30 + t45 * t31;
t58 = m(6) * pkin(4);
t34 = t45 * t39 + t43 * t40;
t29 = t34 * qJD(4);
t57 = t29 * pkin(4);
t50 = t45 * t40;
t33 = -t43 * t39 + t50;
t42 = sin(qJ(5));
t44 = cos(qJ(5));
t21 = t42 * t33 + t44 * t34;
t28 = t33 * qJD(4);
t11 = -qJD(5) * t21 - t42 * t28 - t44 * t29;
t20 = t44 * t33 - t42 * t34;
t55 = t20 * t11;
t10 = qJD(5) * t20 + t44 * t28 - t42 * t29;
t54 = t21 * t10;
t53 = t29 * mrSges(5,1);
t52 = t33 * t29;
t51 = t34 * t28;
t25 = t45 * t30;
t4 = -t11 * mrSges(6,1) + t10 * mrSges(6,2);
t18 = -t43 * t31 - t25;
t27 = t28 * mrSges(5,2);
t48 = t27 + t4;
t14 = -qJD(4) * t25 + qJD(3) * t50 + (-qJD(3) * t39 - qJD(4) * t31) * t43;
t12 = -t29 * pkin(7) + t14;
t15 = -t34 * qJD(3) - t19 * qJD(4);
t13 = -t28 * pkin(7) + t15;
t16 = -t34 * pkin(7) + t18;
t17 = t33 * pkin(7) + t19;
t5 = t44 * t16 - t42 * t17;
t2 = qJD(5) * t5 + t44 * t12 + t42 * t13;
t6 = t42 * t16 + t44 * t17;
t3 = -qJD(5) * t6 - t42 * t12 + t44 * t13;
t47 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t11;
t46 = -cos(pkin(8)) * pkin(1) - t40 * pkin(3) - pkin(2);
t32 = (-mrSges(6,1) * t42 - mrSges(6,2) * t44) * qJD(5) * pkin(4);
t22 = -t33 * pkin(4) + t46;
t1 = [0.2e1 * Ifges(6,1) * t54 + 0.2e1 * Ifges(6,2) * t55 + 0.2e1 * (-t20 * mrSges(6,1) + t21 * mrSges(6,2)) * t57 + 0.2e1 * t22 * t4 + 0.2e1 * t46 * (t27 + t53) - 0.2e1 * Ifges(5,2) * t52 + 0.2e1 * Ifges(5,1) * t51 + 0.2e1 * m(6) * (t6 * t2 + t22 * t57 + t5 * t3) + 0.2e1 * m(5) * (t19 * t14 + t18 * t15) + 0.2e1 * (t10 * t20 + t21 * t11) * Ifges(6,4) + 0.2e1 * (t33 * t28 - t29 * t34) * Ifges(5,4) + 0.2e1 * (-t5 * t10 + t6 * t11 + t2 * t20 - t3 * t21) * mrSges(6,3) + 0.2e1 * (t14 * t33 - t15 * t34 - t18 * t28 - t19 * t29) * mrSges(5,3) + 0.2e1 * (m(4) * t36 + mrSges(4,3)) * qJD(3) * (t39 ^ 2 + t40 ^ 2); m(5) * (t14 * t34 + t15 * t33 - t18 * t29 + t19 * t28) + m(6) * (t6 * t10 + t5 * t11 + t2 * t21 + t3 * t20); 0.2e1 * m(5) * (t51 - t52) + 0.2e1 * m(6) * (t54 + t55); -(-mrSges(5,1) - t58) * t29 + t48; 0; 0; t15 * mrSges(5,1) - t14 * mrSges(5,2) + Ifges(5,5) * t28 - Ifges(5,6) * t29 + (m(6) * (t2 * t42 + t3 * t44 + (-t42 * t5 + t44 * t6) * qJD(5)) + (-t44 * t10 + t42 * t11 + (t20 * t44 + t21 * t42) * qJD(5)) * mrSges(6,3)) * pkin(4) + t47; -t53 + (t10 * t42 + t11 * t44 + (-t20 * t42 + t21 * t44) * qJD(5)) * t58 - t48; 0; 0.2e1 * t32; t47; -t4; 0; t32; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
