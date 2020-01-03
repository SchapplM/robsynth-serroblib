% Calculate time derivative of joint inertia matrix for
% S5RPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:12:21
% EndTime: 2019-12-31 18:12:22
% DurationCPUTime: 0.42s
% Computational Cost: add. (451->100), mult. (1087->139), div. (0->0), fcn. (869->4), ass. (0->45)
t38 = 2 * mrSges(4,3) + 2 * mrSges(5,1);
t29 = sin(pkin(7));
t49 = pkin(6) + qJ(2);
t22 = t49 * t29;
t30 = cos(pkin(7));
t23 = t49 * t30;
t32 = sin(qJ(3));
t54 = cos(qJ(3));
t40 = qJD(3) * t54;
t41 = t54 * t30;
t6 = (qJD(2) * t29 + qJD(3) * t23) * t32 - qJD(2) * t41 + t22 * t40;
t53 = t32 * t29;
t20 = -t41 + t53;
t21 = t54 * t29 + t32 * t30;
t42 = -t30 * pkin(2) - pkin(1);
t34 = -t21 * qJ(4) + t42;
t50 = pkin(3) + qJ(5);
t5 = t50 * t20 + t34;
t56 = 0.2e1 * t5;
t55 = -2 * mrSges(6,1);
t52 = mrSges(4,1) - mrSges(5,2);
t51 = -mrSges(5,1) - mrSges(6,1);
t47 = qJD(4) * t20;
t46 = qJ(4) * qJD(4);
t45 = 0.2e1 * t20;
t44 = 0.2e1 * t21;
t43 = Ifges(6,6) - Ifges(4,4) - Ifges(5,6);
t11 = t54 * t22 + t32 * t23;
t12 = -t32 * t22 + t54 * t23;
t7 = t21 * qJD(2) + t12 * qJD(3);
t37 = t11 * t7 - t12 * t6;
t16 = qJD(3) * t53 - t30 * t40;
t36 = t16 * qJ(4) - t21 * qJD(4);
t17 = t21 * qJD(3);
t15 = t17 * mrSges(6,3);
t14 = t16 * mrSges(4,2);
t13 = t16 * mrSges(5,3);
t10 = t20 * pkin(3) + t34;
t9 = -t20 * pkin(4) + t12;
t8 = t21 * pkin(4) + t11;
t4 = t17 * pkin(3) + t36;
t3 = -t16 * pkin(4) + t7;
t2 = -t17 * pkin(4) - t6;
t1 = t20 * qJD(5) + t50 * t17 + t36;
t18 = [t15 * t56 - 0.2e1 * t42 * t14 + 0.2e1 * t1 * (-t21 * mrSges(6,2) + t20 * mrSges(6,3)) + 0.2e1 * t4 * (-t20 * mrSges(5,2) - t21 * mrSges(5,3)) + 0.2e1 * t10 * t13 + 0.2e1 * m(6) * (t5 * t1 + t9 * t2 + t8 * t3) + 0.2e1 * m(4) * t37 + 0.2e1 * m(5) * (t10 * t4 + t37) + (0.2e1 * t42 * mrSges(4,1) - 0.2e1 * t10 * mrSges(5,2) + t9 * t55 - t12 * t38 + t43 * t44 + (Ifges(6,2) + Ifges(5,3) + Ifges(4,2)) * t45) * t17 + (t8 * t55 + mrSges(6,2) * t56 - t11 * t38 + (-Ifges(4,1) - Ifges(5,2) - Ifges(6,3)) * t44 - t43 * t45) * t16 + (t6 * t20 + t7 * t21) * t38 + 0.2e1 * (-t2 * t20 + t3 * t21) * mrSges(6,1) + 0.2e1 * (m(3) * qJ(2) + mrSges(3,3)) * (t29 ^ 2 + t30 ^ 2) * qJD(2); m(5) * t4 + m(6) * t1 + t16 * mrSges(6,2) + t52 * t17 + t13 - t14 + t15; 0; -mrSges(5,1) * t47 + t2 * mrSges(6,2) - t3 * mrSges(6,3) - t52 * t7 + (mrSges(4,2) - mrSges(5,3)) * t6 + (-qJD(5) * t21 - t47) * mrSges(6,1) + m(6) * (qJ(4) * t2 + qJD(4) * t9 - qJD(5) * t8 - t3 * t50) + m(5) * (-pkin(3) * t7 - qJ(4) * t6 + qJD(4) * t12) + (t51 * qJ(4) + Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t17 + (pkin(3) * mrSges(5,1) + mrSges(6,1) * t50 + Ifges(5,4) - Ifges(4,5) - Ifges(6,5)) * t16; 0; 0.2e1 * m(5) * t46 + 0.2e1 * qJD(5) * mrSges(6,3) + 0.2e1 * m(6) * (qJD(5) * t50 + t46) + 0.2e1 * (mrSges(5,3) + mrSges(6,2)) * qJD(4); m(5) * t7 + m(6) * t3 + t51 * t16; 0; -m(6) * qJD(5); 0; m(6) * t2 - t17 * mrSges(6,1); 0; m(6) * qJD(4); 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;
