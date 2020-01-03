% Calculate time derivative of joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
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
% Datum: 2020-01-03 11:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:25:09
% EndTime: 2020-01-03 11:25:12
% DurationCPUTime: 0.41s
% Computational Cost: add. (364->102), mult. (828->159), div. (0->0), fcn. (604->6), ass. (0->56)
t62 = 2 * qJD(3);
t30 = sin(pkin(8));
t33 = cos(qJ(4));
t42 = qJD(4) * t33;
t23 = t30 * mrSges(6,1) * t42;
t32 = sin(qJ(4));
t43 = qJD(4) * t32;
t12 = -mrSges(6,2) * t30 * t43 + t23;
t16 = (pkin(4) * t42 + qJD(3)) * t30;
t61 = m(6) * t16 + t12;
t31 = cos(pkin(8));
t15 = -cos(pkin(7)) * pkin(1) - pkin(2) - t30 * pkin(6) - t31 * pkin(3);
t25 = sin(pkin(7)) * pkin(1) + qJ(3);
t21 = t33 * t31 * t25;
t8 = t32 * t15 + t21;
t60 = m(6) * pkin(4);
t45 = qJD(3) * t31;
t58 = t15 * t42 + t33 * t45;
t57 = t25 * t32;
t56 = t32 * t30;
t55 = t33 * t30;
t54 = mrSges(5,2) + mrSges(6,2);
t53 = Ifges(6,4) + Ifges(5,4);
t52 = Ifges(6,5) + Ifges(5,5);
t51 = Ifges(5,6) + Ifges(6,6);
t17 = t31 * mrSges(6,2) - mrSges(6,3) * t56;
t18 = t31 * mrSges(5,2) - mrSges(5,3) * t56;
t50 = t17 + t18;
t19 = -t31 * mrSges(6,1) - mrSges(6,3) * t55;
t20 = -t31 * mrSges(5,1) - mrSges(5,3) * t55;
t49 = -t19 - t20;
t47 = qJ(5) * t30;
t46 = qJD(3) * t25;
t44 = qJD(4) * t30;
t41 = qJD(5) * t30;
t40 = 0.2e1 * t30;
t39 = 0.2e1 * t31;
t38 = t31 * t57;
t37 = t33 * t47;
t36 = t32 * t45;
t35 = mrSges(6,1) + t60;
t34 = t33 * t40;
t27 = t31 ^ 2;
t26 = t30 ^ 2;
t22 = t26 * t46;
t14 = (pkin(4) * t32 + t25) * t30;
t13 = (mrSges(5,1) * t33 - mrSges(5,2) * t32) * t44;
t11 = t33 * t15;
t7 = t11 - t38;
t6 = -t32 * t47 + t8;
t5 = -t8 * qJD(4) - t36;
t4 = -qJD(4) * t38 + t58;
t3 = -t37 + t11 + (-pkin(4) - t57) * t31;
t2 = -t36 - t33 * t41 + (-t21 + (-t15 + t47) * t32) * qJD(4);
t1 = -t32 * t41 + (-t37 - t38) * qJD(4) + t58;
t9 = [0.2e1 * t1 * t17 + 0.2e1 * t14 * t12 + 0.2e1 * t4 * t18 + 0.2e1 * t2 * t19 + 0.2e1 * t5 * t20 + (t26 + t27) * mrSges(4,3) * t62 + 0.2e1 * m(4) * (t27 * t46 + t22) + 0.2e1 * m(6) * (t6 * t1 + t14 * t16 + t3 * t2) + 0.2e1 * m(5) * (t8 * t4 + t7 * t5 + t22) + (0.2e1 * t16 * (mrSges(6,1) * t32 + mrSges(6,2) * t33) + 0.2e1 * t25 * t13 + (mrSges(5,1) * t32 + mrSges(5,2) * t33) * t30 * t62 + ((-0.2e1 * t8 * mrSges(5,3) - 0.2e1 * t6 * mrSges(6,3) - t34 * t53 + t39 * t51) * t33 + (0.2e1 * t3 * mrSges(6,3) + 0.2e1 * t7 * mrSges(5,3) + t53 * t32 * t40 + t52 * t39 + (-Ifges(5,1) - Ifges(6,1) + Ifges(5,2) + Ifges(6,2)) * t34) * t32) * qJD(4)) * t30; (-t13 - t61) * t31 + ((-t32 * t50 + t33 * t49) * qJD(4) + m(5) * (-t32 * t5 + t33 * t4 - t42 * t7 - t43 * t8 - t45) + m(6) * (t1 * t33 - t2 * t32 - t3 * t42 - t43 * t6)) * t30 + (mrSges(5,3) + mrSges(6,3)) * t26 * qJD(4) * (-t32 ^ 2 - t33 ^ 2); 0; (t49 * t32 + t50 * t33) * qJD(4) + m(6) * (t32 * t1 + t33 * t2 + (-t3 * t32 + t33 * t6) * qJD(4)) + m(5) * (t32 * t4 + t33 * t5 + (-t32 * t7 + t33 * t8) * qJD(4)); 0; 0; t5 * mrSges(5,1) - t4 * mrSges(5,2) - t1 * mrSges(6,2) + t35 * t2 + (-t51 * t33 + (mrSges(6,3) * pkin(4) - t52) * t32) * t44; -t23 + ((-mrSges(5,1) - t60) * t33 + t54 * t32) * t44; (-t54 * t33 + (-mrSges(5,1) - t35) * t32) * qJD(4); 0; t61; 0; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t9(1), t9(2), t9(4), t9(7), t9(11); t9(2), t9(3), t9(5), t9(8), t9(12); t9(4), t9(5), t9(6), t9(9), t9(13); t9(7), t9(8), t9(9), t9(10), t9(14); t9(11), t9(12), t9(13), t9(14), t9(15);];
Mq = res;
