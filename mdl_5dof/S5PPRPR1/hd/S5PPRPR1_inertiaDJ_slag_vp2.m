% Calculate time derivative of joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
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
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:00:54
% EndTime: 2019-12-05 15:00:54
% DurationCPUTime: 0.27s
% Computational Cost: add. (304->64), mult. (863->115), div. (0->0), fcn. (784->8), ass. (0->38)
t24 = sin(pkin(9));
t26 = cos(pkin(9));
t28 = sin(qJ(5));
t30 = cos(qJ(5));
t32 = t28 * t24 - t30 * t26;
t11 = t32 * qJD(5);
t40 = 2 * m(6);
t17 = t30 * t24 + t28 * t26;
t12 = t17 * qJD(5);
t39 = t32 * t12;
t25 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t18 = t31 * t25 + t29 * t27;
t14 = t18 * qJD(3);
t16 = t29 * t25 - t31 * t27;
t8 = t16 * t14;
t38 = t17 * t11;
t37 = pkin(6) + qJ(4);
t36 = t24 ^ 2 + t26 ^ 2;
t35 = t36 * mrSges(5,3);
t34 = t36 * t18;
t33 = t36 * qJ(4);
t19 = t37 * t24;
t20 = t37 * t26;
t9 = -t30 * t19 - t28 * t20;
t10 = -t28 * t19 + t30 * t20;
t21 = -t26 * pkin(4) - pkin(3);
t13 = t16 * qJD(3);
t7 = t12 * mrSges(6,1) - t11 * mrSges(6,2);
t6 = t32 * t18;
t5 = t17 * t18;
t4 = -t17 * qJD(4) - t10 * qJD(5);
t3 = -t32 * qJD(4) + t9 * qJD(5);
t2 = t18 * t11 + t17 * t13;
t1 = -t18 * t12 + t32 * t13;
t15 = [0.2e1 * m(4) * (-t18 * t13 + t8) + 0.2e1 * m(5) * (-t13 * t34 + t8) + 0.2e1 * m(6) * (-t6 * t1 - t5 * t2 + t8); m(6) * (t17 * t1 + t11 * t6 + t12 * t5 - t2 * t32); (-t38 + t39) * t40; t16 * t7 + (-t26 * mrSges(5,1) + mrSges(6,1) * t32 + t24 * mrSges(5,2) + t17 * mrSges(6,2) - mrSges(4,1)) * t14 - (-mrSges(4,2) + t35) * t13 + m(5) * (-pkin(3) * t14 + qJD(4) * t34 - t13 * t33) + m(6) * (t10 * t1 + t21 * t14 + t9 * t2 - t3 * t6 - t4 * t5) + (-t1 * t32 - t5 * t11 + t6 * t12 - t2 * t17) * mrSges(6,3); m(6) * (-t10 * t11 - t9 * t12 + t3 * t17 - t32 * t4); -0.2e1 * Ifges(6,1) * t38 + 0.2e1 * t21 * t7 + 0.2e1 * Ifges(6,2) * t39 + (t10 * t3 + t9 * t4) * t40 + 0.2e1 * (t11 * t32 - t17 * t12) * Ifges(6,4) + 0.2e1 * (-t10 * t12 + t9 * t11 - t4 * t17 - t3 * t32) * mrSges(6,3) + 0.2e1 * (m(5) * t33 + t35) * qJD(4); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t14; 0; t7; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2); -t7; t4 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,5) * t11 - Ifges(6,6) * t12; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t15(1), t15(2), t15(4), t15(7), t15(11); t15(2), t15(3), t15(5), t15(8), t15(12); t15(4), t15(5), t15(6), t15(9), t15(13); t15(7), t15(8), t15(9), t15(10), t15(14); t15(11), t15(12), t15(13), t15(14), t15(15);];
Mq = res;
