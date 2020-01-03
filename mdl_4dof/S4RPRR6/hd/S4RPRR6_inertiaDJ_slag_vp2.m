% Calculate time derivative of joint inertia matrix for
% S4RPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:52:29
% EndTime: 2019-12-31 16:52:30
% DurationCPUTime: 0.34s
% Computational Cost: add. (634->75), mult. (1441->133), div. (0->0), fcn. (1343->6), ass. (0->39)
t37 = sin(pkin(7));
t47 = pkin(5) + qJ(2);
t32 = t47 * t37;
t38 = cos(pkin(7));
t33 = t47 * t38;
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t20 = -t40 * t32 + t42 * t33;
t31 = t42 * t37 + t40 * t38;
t24 = t31 * qJD(3);
t49 = t24 * pkin(3);
t28 = t42 * t32;
t48 = t42 * t38;
t45 = -t38 * pkin(2) - pkin(1);
t30 = -t40 * t37 + t48;
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t18 = t39 * t30 + t41 * t31;
t23 = t30 * qJD(3);
t10 = -qJD(4) * t18 - t39 * t23 - t41 * t24;
t17 = t41 * t30 - t39 * t31;
t9 = qJD(4) * t17 + t41 * t23 - t39 * t24;
t44 = -t10 * mrSges(5,1) + t9 * mrSges(5,2);
t19 = -t40 * t33 - t28;
t13 = -qJD(3) * t28 + qJD(2) * t48 + (-qJD(2) * t37 - qJD(3) * t33) * t40;
t11 = -t24 * pkin(6) + t13;
t14 = -t31 * qJD(2) - t20 * qJD(3);
t12 = -t23 * pkin(6) + t14;
t15 = -t31 * pkin(6) + t19;
t16 = t30 * pkin(6) + t20;
t4 = t41 * t15 - t39 * t16;
t2 = qJD(4) * t4 + t41 * t11 + t39 * t12;
t5 = t39 * t15 + t41 * t16;
t3 = -qJD(4) * t5 - t39 * t11 + t41 * t12;
t43 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t10;
t26 = (-mrSges(5,1) * t39 - mrSges(5,2) * t41) * qJD(4) * pkin(3);
t22 = t23 * mrSges(4,2);
t21 = -t30 * pkin(3) + t45;
t1 = [0.2e1 * (-t17 * mrSges(5,1) + t18 * mrSges(5,2)) * t49 + 0.2e1 * t21 * t44 + 0.2e1 * t9 * t18 * Ifges(5,1) + 0.2e1 * t10 * Ifges(5,2) * t17 + 0.2e1 * t23 * t31 * Ifges(4,1) - 0.2e1 * t24 * Ifges(4,2) * t30 + 0.2e1 * t45 * (t24 * mrSges(4,1) + t22) + 0.2e1 * m(5) * (t5 * t2 + t21 * t49 + t4 * t3) + 0.2e1 * m(4) * (t20 * t13 + t19 * t14) + 0.2e1 * (t18 * t10 + t9 * t17) * Ifges(5,4) + 0.2e1 * (t23 * t30 - t31 * t24) * Ifges(4,4) + 0.2e1 * (t5 * t10 + t2 * t17 - t3 * t18 - t4 * t9) * mrSges(5,3) + 0.2e1 * (t13 * t30 - t14 * t31 - t19 * t23 - t20 * t24) * mrSges(4,3) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * qJD(2) * (t37 ^ 2 + t38 ^ 2); t22 - (-m(5) * pkin(3) - mrSges(4,1)) * t24 + t44; 0; t14 * mrSges(4,1) - t13 * mrSges(4,2) + Ifges(4,5) * t23 - Ifges(4,6) * t24 + (m(5) * (t2 * t39 + t3 * t41 + (-t39 * t4 + t41 * t5) * qJD(4)) + (t39 * t10 - t41 * t9 + (t17 * t41 + t18 * t39) * qJD(4)) * mrSges(5,3)) * pkin(3) + t43; 0; 0.2e1 * t26; t43; 0; t26; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
