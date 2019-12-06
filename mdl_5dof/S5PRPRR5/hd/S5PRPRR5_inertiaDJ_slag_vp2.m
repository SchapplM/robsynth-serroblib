% Calculate time derivative of joint inertia matrix for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR5_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR5_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:28
% EndTime: 2019-12-05 15:53:30
% DurationCPUTime: 0.60s
% Computational Cost: add. (1020->130), mult. (2577->211), div. (0->0), fcn. (2453->8), ass. (0->59)
t50 = sin(pkin(9));
t67 = pkin(6) + qJ(3);
t44 = t67 * t50;
t51 = cos(pkin(9));
t45 = t67 * t51;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t30 = -t53 * t44 + t56 * t45;
t68 = t51 * t56;
t59 = t50 * t53 - t68;
t35 = t59 * qJD(4);
t70 = m(6) * pkin(4);
t43 = t50 * t56 + t51 * t53;
t36 = t43 * qJD(4);
t69 = pkin(4) * t36;
t40 = t56 * t44;
t66 = t50 ^ 2 + t51 ^ 2;
t54 = sin(qJ(2));
t65 = qJD(2) * t54;
t57 = cos(qJ(2));
t64 = qJD(2) * t57;
t63 = t54 * t64;
t32 = t43 * t54;
t33 = t59 * t54;
t52 = sin(qJ(5));
t55 = cos(qJ(5));
t20 = -t32 * t55 + t33 * t52;
t24 = -t54 * t36 - t59 * t64;
t25 = t54 * t35 - t43 * t64;
t6 = t20 * qJD(5) + t24 * t55 + t25 * t52;
t21 = -t32 * t52 - t33 * t55;
t7 = -t21 * qJD(5) - t24 * t52 + t25 * t55;
t62 = t7 * mrSges(6,1) - t6 * mrSges(6,2);
t47 = -pkin(3) * t51 - pkin(2);
t61 = t66 * mrSges(4,3);
t27 = -t43 * t52 - t55 * t59;
t13 = t27 * qJD(5) - t35 * t55 - t36 * t52;
t28 = t43 * t55 - t52 * t59;
t14 = -t28 * qJD(5) + t35 * t52 - t36 * t55;
t4 = -mrSges(6,1) * t14 + t13 * mrSges(6,2);
t29 = -t45 * t53 - t40;
t60 = t66 * qJ(3);
t22 = -pkin(7) * t43 + t29;
t23 = -pkin(7) * t59 + t30;
t8 = t22 * t55 - t23 * t52;
t9 = t22 * t52 + t23 * t55;
t18 = -qJD(4) * t40 + qJD(3) * t68 + (-qJD(3) * t50 - qJD(4) * t45) * t53;
t16 = -pkin(7) * t36 + t18;
t19 = -t43 * qJD(3) - t30 * qJD(4);
t17 = pkin(7) * t35 + t19;
t2 = t8 * qJD(5) + t16 * t55 + t17 * t52;
t3 = -t9 * qJD(5) - t16 * t52 + t17 * t55;
t58 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t13 + Ifges(6,6) * t14;
t38 = (-mrSges(6,1) * t52 - mrSges(6,2) * t55) * qJD(5) * pkin(4);
t34 = t35 * mrSges(5,2);
t31 = pkin(4) * t59 + t47;
t26 = mrSges(5,1) * t36 - t34;
t15 = -mrSges(6,1) * t27 + mrSges(6,2) * t28;
t1 = [0.2e1 * m(6) * (t20 * t7 + t21 * t6 - t63) + 0.2e1 * m(5) * (-t33 * t24 - t32 * t25 - t63) + 0.2e1 * m(4) * (-0.1e1 + t66) * t63; (-t26 - t4) * t57 + ((-mrSges(3,2) + t61) * t57 + (-mrSges(4,1) * t51 + mrSges(5,1) * t59 + mrSges(4,2) * t50 + mrSges(5,2) * t43 - mrSges(3,1) + t15) * t54) * qJD(2) + m(6) * (t2 * t21 + t3 * t20 + t31 * t65 - t57 * t69 + t9 * t6 + t8 * t7) + m(5) * (-t18 * t33 - t19 * t32 + t24 * t30 + t25 * t29 + t47 * t65) + m(4) * (t66 * t54 * qJD(3) + (-pkin(2) * t54 + t57 * t60) * qJD(2)) + (-t13 * t20 + t14 * t21 + t27 * t6 - t28 * t7) * mrSges(6,3) + (-t24 * t59 - t25 * t43 - t32 * t35 + t33 * t36) * mrSges(5,3); -0.2e1 * t35 * t43 * Ifges(5,1) + 0.2e1 * t13 * t28 * Ifges(6,1) + 0.2e1 * t14 * Ifges(6,2) * t27 + 0.2e1 * t47 * t26 + 0.2e1 * t31 * t4 + 0.2e1 * m(6) * (t2 * t9 + t8 * t3 + t31 * t69) + 0.2e1 * m(5) * (t18 * t30 + t19 * t29) - 0.2e1 * (-Ifges(5,2) * t59 - pkin(4) * t15) * t36 + 0.2e1 * (t13 * t27 + t14 * t28) * Ifges(6,4) + 0.2e1 * (t35 * t59 - t36 * t43) * Ifges(5,4) + 0.2e1 * (-t13 * t8 + t14 * t9 + t2 * t27 - t28 * t3) * mrSges(6,3) + 0.2e1 * (-t18 * t59 - t19 * t43 + t29 * t35 - t30 * t36) * mrSges(5,3) + 0.2e1 * (m(4) * t60 + t61) * qJD(3); (m(4) + m(5) + m(6)) * t65; -t34 - (-mrSges(5,1) - t70) * t36 + t4; 0; t25 * mrSges(5,1) - t24 * mrSges(5,2) + (t52 * t6 + t55 * t7 + (-t20 * t52 + t21 * t55) * qJD(5)) * t70 + t62; t19 * mrSges(5,1) - t18 * mrSges(5,2) - Ifges(5,5) * t35 - Ifges(5,6) * t36 + (m(6) * (t2 * t52 + t3 * t55 + (-t52 * t8 + t55 * t9) * qJD(5)) + (-t55 * t13 + t52 * t14 + (t27 * t55 + t28 * t52) * qJD(5)) * mrSges(6,3)) * pkin(4) + t58; 0; 0.2e1 * t38; t62; t58; 0; t38; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
