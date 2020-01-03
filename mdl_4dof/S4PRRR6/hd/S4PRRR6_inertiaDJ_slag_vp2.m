% Calculate time derivative of joint inertia matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:41
% DurationCPUTime: 0.36s
% Computational Cost: add. (342->81), mult. (961->146), div. (0->0), fcn. (759->6), ass. (0->41)
t30 = sin(qJ(2));
t28 = sin(qJ(4));
t29 = sin(qJ(3));
t31 = cos(qJ(4));
t32 = cos(qJ(3));
t34 = t28 * t29 - t31 * t32;
t16 = t34 * t30;
t26 = t29 ^ 2;
t43 = t32 ^ 2 + t26;
t47 = qJD(3) + qJD(4);
t46 = 2 * m(5);
t45 = -pkin(6) - pkin(5);
t33 = cos(qJ(2));
t42 = qJD(2) * t33;
t41 = qJD(3) * t29;
t40 = qJD(3) * t32;
t39 = pkin(3) * t41;
t19 = t28 * t32 + t31 * t29;
t11 = t47 * t19;
t3 = -t11 * t30 - t34 * t42;
t4 = t47 * t16 - t19 * t42;
t38 = t4 * mrSges(5,1) - t3 * mrSges(5,2);
t37 = qJD(3) * t45;
t10 = t47 * t34;
t23 = t45 * t29;
t24 = t45 * t32;
t13 = t31 * t23 + t28 * t24;
t21 = t29 * t37;
t22 = t32 * t37;
t6 = qJD(4) * t13 + t31 * t21 + t28 * t22;
t14 = t28 * t23 - t31 * t24;
t7 = -qJD(4) * t14 - t28 * t21 + t31 * t22;
t36 = t7 * mrSges(5,1) - t6 * mrSges(5,2) - Ifges(5,5) * t10 - Ifges(5,6) * t11;
t35 = -t32 * mrSges(4,1) + t29 * mrSges(4,2);
t25 = -t32 * pkin(3) - pkin(2);
t20 = (mrSges(4,1) * t29 + mrSges(4,2) * t32) * qJD(3);
t17 = (-mrSges(5,1) * t28 - mrSges(5,2) * t31) * qJD(4) * pkin(3);
t15 = t19 * t30;
t12 = mrSges(5,1) * t34 + t19 * mrSges(5,2);
t1 = t11 * mrSges(5,1) - t10 * mrSges(5,2);
t2 = [(-t15 * t4 - t16 * t3) * t46 + 0.4e1 * (m(4) * (-0.1e1 + t43) / 0.2e1 - m(5) / 0.2e1) * t30 * t42; m(5) * (t13 * t4 + t14 * t3 - t7 * t15 - t6 * t16) + (-t15 * t10 + t16 * t11 - t4 * t19 - t3 * t34) * mrSges(5,3) + (-m(4) * pkin(2) + m(5) * t25 - mrSges(3,1) + t12 + t35) * qJD(2) * t30 + (-t20 - t1 - m(5) * t39 + (-mrSges(3,2) + (m(4) * pkin(5) + mrSges(4,3)) * t43) * qJD(2)) * t33; 0.2e1 * t12 * t39 + 0.2e1 * t25 * t1 + 0.2e1 * t11 * Ifges(5,2) * t34 - 0.2e1 * t10 * t19 * Ifges(5,1) + (t13 * t7 + t14 * t6 + t25 * t39) * t46 - 0.2e1 * pkin(2) * t20 + 0.2e1 * (-t26 * qJD(3) + t32 * t40) * Ifges(4,4) + 0.2e1 * (t10 * t34 - t11 * t19) * Ifges(5,4) + 0.2e1 * (t13 * t10 - t14 * t11 - t7 * t19 - t34 * t6) * mrSges(5,3) + 0.2e1 * (Ifges(4,1) - Ifges(4,2)) * t29 * t40; (t30 * t41 - t32 * t42) * mrSges(4,2) + (-t29 * t42 - t30 * t40) * mrSges(4,1) + m(5) * (t28 * t3 + t31 * t4 + (t15 * t28 - t16 * t31) * qJD(4)) * pkin(3) + t38; (Ifges(4,5) * t32 - Ifges(4,6) * t29 + t35 * pkin(5)) * qJD(3) + (m(5) * (t28 * t6 + t31 * t7 + (-t13 * t28 + t14 * t31) * qJD(4)) + (t31 * t10 - t28 * t11 + (t19 * t28 - t31 * t34) * qJD(4)) * mrSges(5,3)) * pkin(3) + t36; 0.2e1 * t17; t38; t36; t17; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
