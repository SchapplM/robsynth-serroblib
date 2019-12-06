% Calculate time derivative of joint inertia matrix for
% S5PPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
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
% Datum: 2019-12-05 14:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPPRR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPPRR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPPRR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:57:51
% EndTime: 2019-12-05 14:57:52
% DurationCPUTime: 0.23s
% Computational Cost: add. (187->48), mult. (646->94), div. (0->0), fcn. (598->8), ass. (0->34)
t19 = sin(pkin(9));
t21 = cos(pkin(9));
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t40 = -t24 * t19 + t26 * t21;
t39 = m(6) * pkin(6) + mrSges(6,3);
t23 = sin(qJ(5));
t25 = cos(qJ(5));
t15 = -t25 * mrSges(6,1) + t23 * mrSges(6,2);
t38 = -m(6) * pkin(4) - mrSges(5,1) + t15;
t18 = t25 ^ 2;
t20 = sin(pkin(8));
t33 = qJD(4) * t20;
t6 = t40 * t33;
t12 = t26 * t19 + t24 * t21;
t7 = t12 * t20;
t37 = t7 * t6;
t10 = t12 * qJD(4);
t36 = t40 * t10;
t32 = qJD(5) * t23;
t9 = t40 * qJD(4);
t31 = (t23 ^ 2 + t18) * t9;
t30 = t10 * t7 - t40 * t6;
t22 = cos(pkin(8));
t8 = t40 * t20;
t28 = t22 * t23 - t25 * t8;
t3 = -t22 * t25 - t23 * t8;
t29 = -t23 * t3 - t25 * t28;
t5 = t12 * t33;
t1 = qJD(5) * t3 - t25 * t5;
t2 = t28 * qJD(5) + t23 * t5;
t27 = t1 * t25 - t2 * t23 + (t23 * t28 - t25 * t3) * qJD(5);
t13 = (mrSges(6,1) * t23 + mrSges(6,2) * t25) * qJD(5);
t4 = [0.2e1 * m(5) * (-t8 * t5 + t37) + 0.2e1 * m(6) * (-t1 * t28 + t3 * t2 + t37); m(5) * (-t12 * t5 + t9 * t8 + t30) + m(6) * (t27 * t12 + t29 * t9 + t30); 0.2e1 * m(5) * (t12 * t9 - t36) + 0.2e1 * m(6) * (t12 * t31 - t36); m(6) * (t29 * qJD(5) + t23 * t1 + t25 * t2); 0; 0; t5 * mrSges(5,2) + t7 * t13 + t39 * t27 + t38 * t6; -t9 * mrSges(5,2) + t38 * t10 - t13 * t40 + t39 * t31; 0; 0.2e1 * qJD(5) * t18 * Ifges(6,4) - 0.2e1 * pkin(4) * t13 + 0.2e1 * (-Ifges(6,4) * t23 + (Ifges(6,1) - Ifges(6,2)) * t25) * t32; t2 * mrSges(6,1) - t1 * mrSges(6,2); (t12 * t32 - t25 * t9) * mrSges(6,2) + (-qJD(5) * t12 * t25 - t23 * t9) * mrSges(6,1); -t13; (Ifges(6,5) * t25 - Ifges(6,6) * t23 + t15 * pkin(6)) * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
