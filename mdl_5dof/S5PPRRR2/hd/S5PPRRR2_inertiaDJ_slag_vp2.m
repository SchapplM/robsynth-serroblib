% Calculate time derivative of joint inertia matrix for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:12
% EndTime: 2019-12-05 15:14:13
% DurationCPUTime: 0.40s
% Computational Cost: add. (568->97), mult. (1535->171), div. (0->0), fcn. (1380->8), ass. (0->52)
t37 = sin(qJ(4));
t40 = cos(qJ(4));
t26 = (t37 * mrSges(5,1) + t40 * mrSges(5,2)) * qJD(4);
t36 = sin(qJ(5));
t39 = cos(qJ(5));
t43 = t36 * t37 - t39 * t40;
t58 = qJD(4) + qJD(5);
t13 = t58 * t43;
t25 = t36 * t40 + t39 * t37;
t14 = t58 * t25;
t4 = t14 * mrSges(6,1) - t13 * mrSges(6,2);
t59 = t26 + t4;
t34 = sin(pkin(9));
t35 = cos(pkin(9));
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t23 = t41 * t34 + t38 * t35;
t9 = t43 * t23;
t57 = 2 * m(6);
t32 = t37 ^ 2;
t56 = m(6) * pkin(4);
t55 = -pkin(7) - pkin(6);
t20 = t23 * qJD(3);
t22 = t38 * t34 - t41 * t35;
t15 = t22 * t20;
t54 = t43 * t14;
t53 = t25 * t13;
t52 = t40 ^ 2 + t32;
t51 = qJD(4) * t37;
t50 = qJD(4) * t40;
t49 = pkin(4) * t51;
t19 = t22 * qJD(3);
t2 = -t14 * t23 + t43 * t19;
t3 = t25 * t19 + t58 * t9;
t48 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t47 = qJD(4) * t55;
t46 = t52 * t19;
t45 = -t40 * mrSges(5,1) + t37 * mrSges(5,2);
t29 = t55 * t37;
t30 = t55 * t40;
t17 = t39 * t29 + t36 * t30;
t18 = t36 * t29 - t39 * t30;
t27 = t37 * t47;
t28 = t40 * t47;
t6 = qJD(5) * t17 + t39 * t27 + t36 * t28;
t7 = -qJD(5) * t18 - t36 * t27 + t39 * t28;
t42 = t7 * mrSges(6,1) - t6 * mrSges(6,2) - Ifges(6,5) * t13 - Ifges(6,6) * t14;
t31 = -t40 * pkin(4) - pkin(3);
t21 = (-mrSges(6,1) * t36 - mrSges(6,2) * t39) * qJD(5) * pkin(4);
t16 = mrSges(6,1) * t43 + t25 * mrSges(6,2);
t8 = t25 * t23;
t1 = [0.2e1 * m(6) * (-t9 * t2 - t8 * t3 + t15) + 0.2e1 * m(5) * (-t23 * t46 + t15) + 0.2e1 * m(4) * (-t23 * t19 + t15); m(6) * (t13 * t9 + t14 * t8 + t25 * t2 - t3 * t43); (-t53 + t54) * t57; t59 * t22 + (-mrSges(4,1) + t16 + t45) * t20 - (t52 * mrSges(5,3) - mrSges(4,2)) * t19 + m(6) * (t17 * t3 + t18 * t2 + t31 * t20 + t22 * t49 - t6 * t9 - t7 * t8) + m(5) * (-pkin(3) * t20 - pkin(6) * t46) + (-t8 * t13 + t9 * t14 - t2 * t43 - t3 * t25) * mrSges(6,3); m(6) * (-t18 * t13 - t17 * t14 + t6 * t25 - t43 * t7); (t17 * t7 + t18 * t6 + t31 * t49) * t57 - 0.2e1 * Ifges(6,1) * t53 + 0.2e1 * Ifges(6,2) * t54 + 0.2e1 * t16 * t49 + 0.2e1 * t31 * t4 - 0.2e1 * pkin(3) * t26 + 0.2e1 * (-t32 * qJD(4) + t40 * t50) * Ifges(5,4) + 0.2e1 * (t13 * t43 - t25 * t14) * Ifges(6,4) + 0.2e1 * (t17 * t13 - t18 * t14 - t7 * t25 - t43 * t6) * mrSges(6,3) + 0.2e1 * (Ifges(5,1) - Ifges(5,2)) * t37 * t50; (t40 * t19 + t23 * t51) * mrSges(5,2) + (t37 * t19 - t23 * t50) * mrSges(5,1) + (t2 * t36 + t3 * t39 + (t36 * t8 - t39 * t9) * qJD(5)) * t56 + t48; (-t13 * t36 - t14 * t39 + (t25 * t39 + t36 * t43) * qJD(5)) * t56 - t59; (Ifges(5,5) * t40 - Ifges(5,6) * t37 + t45 * pkin(6)) * qJD(4) + (m(6) * (t36 * t6 + t39 * t7 + (-t17 * t36 + t18 * t39) * qJD(5)) + (t39 * t13 - t36 * t14 + (t25 * t36 - t39 * t43) * qJD(5)) * mrSges(6,3)) * pkin(4) + t42; 0.2e1 * t21; t48; -t4; t42; t21; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
