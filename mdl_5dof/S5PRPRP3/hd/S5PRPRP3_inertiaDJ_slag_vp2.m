% Calculate time derivative of joint inertia matrix for
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:40
% EndTime: 2019-12-05 15:32:41
% DurationCPUTime: 0.26s
% Computational Cost: add. (205->75), mult. (540->117), div. (0->0), fcn. (380->6), ass. (0->37)
t41 = -2 * mrSges(6,3);
t40 = m(6) * pkin(4);
t20 = sin(pkin(8));
t21 = cos(pkin(8));
t23 = sin(qJ(2));
t25 = cos(qJ(2));
t9 = t20 * t25 + t21 * t23;
t4 = t9 * qJD(2);
t8 = t20 * t23 - t21 * t25;
t39 = t8 * t4;
t24 = cos(qJ(4));
t14 = t20 * pkin(2) + pkin(6);
t33 = qJ(5) + t14;
t7 = t33 * t24;
t38 = t24 * t7;
t37 = mrSges(5,2) * t24;
t36 = mrSges(5,2) + mrSges(6,2);
t35 = Ifges(5,4) + Ifges(6,4);
t31 = qJD(4) * t24;
t22 = sin(qJ(4));
t32 = qJD(4) * t22;
t10 = mrSges(6,1) * t32 + mrSges(6,2) * t31;
t34 = t22 ^ 2 + t24 ^ 2;
t30 = 0.2e1 * t22;
t29 = mrSges(6,1) + t40;
t15 = -t21 * pkin(2) - pkin(3);
t5 = t8 * qJD(2);
t28 = t34 * t5;
t27 = qJD(4) * t33;
t26 = -mrSges(5,1) - t29;
t13 = -t24 * mrSges(6,1) + t22 * mrSges(6,2);
t12 = -t24 * pkin(4) + t15;
t11 = (mrSges(5,1) * t22 + t37) * qJD(4);
t6 = t33 * t22;
t3 = -t22 * qJD(5) - t24 * t27;
t2 = t24 * qJD(5) - t22 * t27;
t1 = [0.2e1 * m(4) * (-t9 * t5 + t39) + 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t28 * t9 + t39); (t11 + t10) * t8 + (-t23 * mrSges(3,1) - t25 * mrSges(3,2)) * qJD(2) + (-t24 * mrSges(5,1) + t22 * mrSges(5,2) - mrSges(4,1) + t13) * t4 - (-mrSges(4,2) + (mrSges(5,3) + mrSges(6,3)) * t34) * t5 + m(5) * (-t14 * t28 + t15 * t4) + m(4) * (-t20 * t5 - t21 * t4) * pkin(2) + (t12 * t4 - t5 * t38 + (pkin(4) * qJD(4) * t8 - t5 * t6) * t22 + (t2 * t24 + t31 * t6 + (-qJD(4) * t7 - t3) * t22) * t9) * m(6); 0.2e1 * t12 * t10 + 0.2e1 * t15 * t11 + 0.2e1 * m(6) * (t7 * t2 - t6 * t3) + (t3 * t41 + (t7 * t41 - t35 * t30 + 0.2e1 * (m(6) * t12 + t13) * pkin(4)) * qJD(4)) * t22 + (0.2e1 * t2 * mrSges(6,3) + (-t6 * t41 + 0.2e1 * t35 * t24 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2)) * t30) * qJD(4)) * t24; 0; m(6) * (t22 * t2 + t24 * t3 + (t22 * t6 + t38) * qJD(4)); 0; -(t22 * t26 - t24 * t36) * t5 + (t22 * t36 + t24 * t26) * t9 * qJD(4); -t2 * mrSges(6,2) + t29 * t3 + ((mrSges(5,2) * t14 - Ifges(5,6) - Ifges(6,6)) * t22 + (-mrSges(5,1) * t14 - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t24) * qJD(4); (-t37 + (-mrSges(5,1) - t40) * t22) * qJD(4) - t10; 0; m(6) * t4; t32 * t40 + t10; 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
