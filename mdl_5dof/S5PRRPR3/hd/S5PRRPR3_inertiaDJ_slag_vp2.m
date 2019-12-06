% Calculate time derivative of joint inertia matrix for
% S5PRRPR3
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:04
% EndTime: 2019-12-05 16:19:05
% DurationCPUTime: 0.47s
% Computational Cost: add. (926->114), mult. (2128->192), div. (0->0), fcn. (1896->6), ass. (0->55)
t50 = cos(qJ(3));
t44 = -t50 * pkin(3) - pkin(2);
t65 = 0.2e1 * t44;
t64 = m(5) * pkin(3);
t45 = sin(pkin(9));
t63 = pkin(3) * t45;
t46 = cos(pkin(9));
t48 = sin(qJ(3));
t36 = -t45 * t48 + t46 * t50;
t37 = t45 * t50 + t46 * t48;
t47 = sin(qJ(5));
t49 = cos(qJ(5));
t19 = t47 * t36 + t49 * t37;
t34 = t37 * qJD(3);
t35 = t36 * qJD(3);
t11 = -t19 * qJD(5) - t49 * t34 - t47 * t35;
t18 = t49 * t36 - t47 * t37;
t62 = t18 * t11;
t10 = t18 * qJD(5) - t47 * t34 + t49 * t35;
t61 = t19 * t10;
t60 = t36 * t34;
t59 = t37 * t35;
t58 = -qJ(4) - pkin(6);
t54 = qJD(3) * t58;
t32 = t50 * qJD(4) + t48 * t54;
t33 = -t48 * qJD(4) + t50 * t54;
t15 = t46 * t32 + t45 * t33;
t41 = t58 * t48;
t42 = t58 * t50;
t21 = t45 * t41 - t46 * t42;
t57 = 0.2e1 * t50;
t56 = qJD(3) * t48 * pkin(3);
t4 = -t11 * mrSges(6,1) + t10 * mrSges(6,2);
t55 = t34 * mrSges(5,1) + t35 * mrSges(5,2);
t14 = -t45 * t32 + t46 * t33;
t20 = t46 * t41 + t45 * t42;
t12 = -t35 * pkin(7) + t14;
t13 = -t34 * pkin(7) + t15;
t16 = -t37 * pkin(7) + t20;
t17 = t36 * pkin(7) + t21;
t5 = t49 * t16 - t47 * t17;
t2 = t5 * qJD(5) + t47 * t12 + t49 * t13;
t6 = t47 * t16 + t49 * t17;
t3 = -t6 * qJD(5) + t49 * t12 - t47 * t13;
t53 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t10 + Ifges(6,6) * t11;
t43 = t46 * pkin(3) + pkin(4);
t30 = t49 * t43 - t47 * t63;
t27 = t30 * qJD(5);
t31 = t47 * t43 + t49 * t63;
t28 = t31 * qJD(5);
t52 = -t28 * mrSges(6,1) - t27 * mrSges(6,2);
t51 = -t4 - t55;
t23 = -t36 * pkin(4) + t44;
t22 = t34 * pkin(4) + t56;
t1 = [0.2e1 * m(5) * (t59 - t60) + 0.2e1 * m(6) * (t61 + t62); m(5) * (t14 * t36 + t15 * t37 - t20 * t34 + t21 * t35) + m(6) * (t6 * t10 + t5 * t11 + t3 * t18 + t2 * t19); t55 * t65 - 0.2e1 * Ifges(5,2) * t60 + 0.2e1 * Ifges(5,1) * t59 + 0.2e1 * Ifges(6,2) * t62 + 0.2e1 * Ifges(6,1) * t61 + 0.2e1 * t22 * (-t18 * mrSges(6,1) + t19 * mrSges(6,2)) + 0.2e1 * t23 * t4 + 0.2e1 * m(5) * (t20 * t14 + t21 * t15) + 0.2e1 * m(6) * (t6 * t2 + t23 * t22 + t5 * t3) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) * t50) * t57 + (t64 * t65 - 0.2e1 * pkin(2) * mrSges(4,1) + 0.2e1 * pkin(3) * (-t36 * mrSges(5,1) + t37 * mrSges(5,2)) - 0.2e1 * Ifges(4,4) * t48 + (Ifges(4,1) - Ifges(4,2)) * t57) * t48) * qJD(3) + 0.2e1 * (t18 * t10 + t11 * t19) * Ifges(6,4) + 0.2e1 * (-t34 * t37 + t36 * t35) * Ifges(5,4) + 0.2e1 * (-t5 * t10 + t6 * t11 + t2 * t18 - t3 * t19) * mrSges(6,3) + 0.2e1 * (-t14 * t37 + t15 * t36 - t20 * t35 - t21 * t34) * mrSges(5,3); (-t48 * mrSges(4,1) - t50 * mrSges(4,2)) * qJD(3) + m(6) * (t31 * t10 + t30 * t11 - t28 * t18 + t27 * t19) + (-t34 * t46 + t35 * t45) * t64 + t51; Ifges(5,5) * t35 - Ifges(5,6) * t34 + m(6) * (t31 * t2 + t27 * t6 - t28 * t5 + t30 * t3) + t14 * mrSges(5,1) - t15 * mrSges(5,2) + (Ifges(4,5) * t50 - Ifges(4,6) * t48 + (-mrSges(4,1) * t50 + mrSges(4,2) * t48) * pkin(6)) * qJD(3) + (m(5) * (t14 * t46 + t15 * t45) + (-t45 * t34 - t46 * t35) * mrSges(5,3)) * pkin(3) + (-t30 * t10 + t31 * t11 + t27 * t18 + t28 * t19) * mrSges(6,3) + t53; 0.2e1 * m(6) * (t27 * t31 - t30 * t28) + 0.2e1 * t52; 0; m(5) * t56 + m(6) * t22 - t51; 0; 0; -t4; t53; t52; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
