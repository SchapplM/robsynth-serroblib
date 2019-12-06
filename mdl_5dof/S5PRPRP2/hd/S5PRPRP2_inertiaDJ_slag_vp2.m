% Calculate time derivative of joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:30:32
% EndTime: 2019-12-05 15:30:34
% DurationCPUTime: 0.41s
% Computational Cost: add. (284->100), mult. (748->157), div. (0->0), fcn. (524->4), ass. (0->55)
t61 = 2 * qJD(3);
t29 = sin(pkin(8));
t32 = cos(qJ(4));
t42 = qJD(4) * t32;
t12 = (pkin(4) * t42 + qJD(3)) * t29;
t21 = t29 * mrSges(6,1) * t42;
t31 = sin(qJ(4));
t43 = qJD(4) * t31;
t9 = -mrSges(6,2) * t29 * t43 + t21;
t60 = m(6) * t12 + t9;
t30 = cos(pkin(8));
t20 = -t30 * pkin(3) - t29 * pkin(6) - pkin(2);
t23 = t32 * t30 * qJ(3);
t8 = t31 * t20 + t23;
t59 = m(6) * pkin(4);
t57 = t31 * t29;
t56 = t32 * t29;
t55 = mrSges(5,2) + mrSges(6,2);
t54 = Ifges(6,4) + Ifges(5,4);
t53 = Ifges(6,5) + Ifges(5,5);
t52 = Ifges(5,6) + Ifges(6,6);
t45 = qJD(3) * t30;
t51 = t20 * t42 + t32 * t45;
t15 = t30 * mrSges(6,2) - mrSges(6,3) * t57;
t16 = t30 * mrSges(5,2) - mrSges(5,3) * t57;
t50 = t15 + t16;
t17 = -t30 * mrSges(6,1) - mrSges(6,3) * t56;
t18 = -t30 * mrSges(5,1) - mrSges(5,3) * t56;
t49 = -t17 - t18;
t47 = qJ(3) * t31;
t46 = qJ(5) * t29;
t44 = qJD(4) * t29;
t41 = qJD(5) * t29;
t40 = qJ(3) * qJD(3);
t39 = 0.2e1 * t29;
t38 = 0.2e1 * t30;
t37 = t30 * t47;
t36 = t32 * t46;
t35 = t31 * t45;
t34 = mrSges(6,1) + t59;
t33 = t32 * t39;
t26 = t30 ^ 2;
t25 = t29 ^ 2;
t24 = t25 * t40;
t19 = (pkin(4) * t31 + qJ(3)) * t29;
t14 = t32 * t20;
t10 = (mrSges(5,1) * t32 - mrSges(5,2) * t31) * t44;
t7 = t14 - t37;
t6 = -t31 * t46 + t8;
t5 = -t8 * qJD(4) - t35;
t4 = -qJD(4) * t37 + t51;
t3 = -t36 + t14 + (-pkin(4) - t47) * t30;
t2 = -t35 - t32 * t41 + (-t23 + (-t20 + t46) * t31) * qJD(4);
t1 = -t31 * t41 + (-t36 - t37) * qJD(4) + t51;
t11 = [0; (-t10 - t60) * t30 + ((-t31 * t50 + t32 * t49) * qJD(4) + m(5) * (-t31 * t5 + t32 * t4 - t42 * t7 - t43 * t8 - t45) + m(6) * (t1 * t32 - t2 * t31 - t3 * t42 - t43 * t6)) * t29 + (mrSges(5,3) + mrSges(6,3)) * t25 * qJD(4) * (-t31 ^ 2 - t32 ^ 2); 0.2e1 * t1 * t15 + 0.2e1 * t4 * t16 + 0.2e1 * t2 * t17 + 0.2e1 * t5 * t18 + 0.2e1 * t19 * t9 + (t25 + t26) * mrSges(4,3) * t61 + 0.2e1 * m(5) * (t8 * t4 + t7 * t5 + t24) + 0.2e1 * m(4) * (t26 * t40 + t24) + 0.2e1 * m(6) * (t6 * t1 + t19 * t12 + t3 * t2) + (0.2e1 * t12 * (mrSges(6,1) * t31 + mrSges(6,2) * t32) + 0.2e1 * qJ(3) * t10 + (mrSges(5,1) * t31 + mrSges(5,2) * t32) * t29 * t61 + ((-0.2e1 * t8 * mrSges(5,3) - 0.2e1 * t6 * mrSges(6,3) - t33 * t54 + t38 * t52) * t32 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * t7 * mrSges(5,3) + t54 * t31 * t39 + t53 * t38 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t33) * t31) * qJD(4)) * t29; 0; (t49 * t31 + t50 * t32) * qJD(4) + m(6) * (t31 * t1 + t32 * t2 + (-t3 * t31 + t32 * t6) * qJD(4)) + m(5) * (t31 * t4 + t32 * t5 + (-t31 * t7 + t32 * t8) * qJD(4)); 0; -t21 + ((-mrSges(5,1) - t59) * t32 + t55 * t31) * t44; t5 * mrSges(5,1) - t4 * mrSges(5,2) - t1 * mrSges(6,2) + t34 * t2 + (-t52 * t32 + (mrSges(6,3) * pkin(4) - t53) * t31) * t44; (-t55 * t32 + (-mrSges(5,1) - t34) * t31) * qJD(4); 0; 0; t60; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t11(1), t11(2), t11(4), t11(7), t11(11); t11(2), t11(3), t11(5), t11(8), t11(12); t11(4), t11(5), t11(6), t11(9), t11(13); t11(7), t11(8), t11(9), t11(10), t11(14); t11(11), t11(12), t11(13), t11(14), t11(15);];
Mq = res;
