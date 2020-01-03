% Calculate time derivative of joint inertia matrix for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:24
% DurationCPUTime: 0.38s
% Computational Cost: add. (704->95), mult. (1600->167), div. (0->0), fcn. (1400->6), ass. (0->49)
t49 = cos(qJ(2));
t43 = -t49 * pkin(2) - pkin(1);
t59 = 0.2e1 * t43;
t44 = sin(pkin(7));
t58 = pkin(2) * t44;
t57 = -qJ(3) - pkin(5);
t47 = sin(qJ(2));
t52 = qJD(2) * t57;
t31 = t49 * qJD(3) + t47 * t52;
t32 = -t47 * qJD(3) + t49 * t52;
t45 = cos(pkin(7));
t14 = t45 * t31 + t44 * t32;
t40 = t57 * t47;
t41 = t57 * t49;
t20 = t44 * t40 - t45 * t41;
t56 = 0.2e1 * t49;
t55 = qJD(2) * t47 * pkin(2);
t35 = -t44 * t47 + t45 * t49;
t36 = t44 * t49 + t45 * t47;
t46 = sin(qJ(4));
t48 = cos(qJ(4));
t18 = t46 * t35 + t48 * t36;
t33 = t36 * qJD(2);
t34 = t35 * qJD(2);
t10 = -t18 * qJD(4) - t48 * t33 - t46 * t34;
t17 = t48 * t35 - t46 * t36;
t9 = t17 * qJD(4) - t46 * t33 + t48 * t34;
t54 = -t10 * mrSges(5,1) + t9 * mrSges(5,2);
t53 = t33 * mrSges(4,1) + t34 * mrSges(4,2);
t13 = -t44 * t31 + t45 * t32;
t19 = t45 * t40 + t44 * t41;
t11 = -t34 * pkin(6) + t13;
t12 = -t33 * pkin(6) + t14;
t15 = -t36 * pkin(6) + t19;
t16 = t35 * pkin(6) + t20;
t4 = t48 * t15 - t46 * t16;
t2 = t4 * qJD(4) + t46 * t11 + t48 * t12;
t5 = t46 * t15 + t48 * t16;
t3 = -t5 * qJD(4) + t48 * t11 - t46 * t12;
t51 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t10;
t42 = t45 * pkin(2) + pkin(3);
t29 = t48 * t42 - t46 * t58;
t26 = t29 * qJD(4);
t30 = t46 * t42 + t48 * t58;
t27 = t30 * qJD(4);
t50 = -t27 * mrSges(5,1) - t26 * mrSges(5,2);
t22 = -t35 * pkin(3) + t43;
t21 = t33 * pkin(3) + t55;
t1 = [-0.2e1 * t35 * Ifges(4,2) * t33 + 0.2e1 * t36 * t34 * Ifges(4,1) + t53 * t59 + 0.2e1 * t17 * Ifges(5,2) * t10 + 0.2e1 * t9 * t18 * Ifges(5,1) + 0.2e1 * t21 * (-t17 * mrSges(5,1) + t18 * mrSges(5,2)) + 0.2e1 * t22 * t54 + 0.2e1 * m(4) * (t19 * t13 + t20 * t14) + 0.2e1 * m(5) * (t5 * t2 + t22 * t21 + t4 * t3) + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t49) * t56 + (-0.2e1 * pkin(1) * mrSges(3,1) + 0.2e1 * pkin(2) * (-t35 * mrSges(4,1) + t36 * mrSges(4,2)) + m(4) * pkin(2) * t59 - 0.2e1 * Ifges(3,4) * t47 + (Ifges(3,1) - Ifges(3,2)) * t56) * t47) * qJD(2) + 0.2e1 * (t10 * t18 + t17 * t9) * Ifges(5,4) + 0.2e1 * (-t36 * t33 + t35 * t34) * Ifges(4,4) + 0.2e1 * (t5 * t10 + t2 * t17 - t3 * t18 - t4 * t9) * mrSges(5,3) + 0.2e1 * (-t13 * t36 + t14 * t35 - t19 * t34 - t20 * t33) * mrSges(4,3); m(5) * (t30 * t2 + t26 * t5 - t27 * t4 + t29 * t3) + t13 * mrSges(4,1) - t14 * mrSges(4,2) - Ifges(4,6) * t33 + Ifges(4,5) * t34 + (Ifges(3,5) * t49 - Ifges(3,6) * t47 + (-mrSges(3,1) * t49 + mrSges(3,2) * t47) * pkin(5)) * qJD(2) + (m(4) * (t13 * t45 + t14 * t44) + (-t44 * t33 - t45 * t34) * mrSges(4,3)) * pkin(2) + (t30 * t10 + t26 * t17 + t27 * t18 - t29 * t9) * mrSges(5,3) + t51; 0.2e1 * m(5) * (t30 * t26 - t29 * t27) + 0.2e1 * t50; m(4) * t55 + m(5) * t21 + t53 + t54; 0; 0; t51; t50; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
