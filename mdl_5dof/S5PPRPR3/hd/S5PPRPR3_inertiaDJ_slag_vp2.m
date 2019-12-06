% Calculate time derivative of joint inertia matrix for
% S5PPRPR3
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
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:04:48
% EndTime: 2019-12-05 15:04:49
% DurationCPUTime: 0.24s
% Computational Cost: add. (210->60), mult. (703->113), div. (0->0), fcn. (639->8), ass. (0->37)
t21 = sin(pkin(9));
t23 = cos(pkin(9));
t26 = sin(qJ(3));
t28 = cos(qJ(3));
t11 = t21 * t26 - t23 * t28;
t27 = cos(qJ(5));
t20 = t27 ^ 2;
t42 = m(5) * pkin(3);
t22 = sin(pkin(8));
t35 = qJD(3) * t22;
t5 = t11 * t35;
t12 = t21 * t28 + t23 * t26;
t7 = t12 * t22;
t41 = t7 * t5;
t9 = t12 * qJD(3);
t40 = t11 * t9;
t25 = sin(qJ(5));
t15 = -t27 * mrSges(6,1) + t25 * mrSges(6,2);
t37 = mrSges(5,1) - t15;
t36 = t25 ^ 2 + t20;
t34 = qJD(5) * t25;
t10 = t11 * qJD(3);
t33 = t10 * t36;
t32 = -t11 * t5 + t9 * t7;
t24 = cos(pkin(8));
t8 = t11 * t22;
t3 = -t24 * t27 + t25 * t8;
t30 = t24 * t25 + t27 * t8;
t31 = -t25 * t3 - t27 * t30;
t6 = t12 * t35;
t1 = qJD(5) * t3 - t27 * t6;
t2 = t30 * qJD(5) + t25 * t6;
t29 = t1 * t27 - t2 * t25 + (t25 * t30 - t27 * t3) * qJD(5);
t18 = -t23 * pkin(3) - pkin(4);
t17 = t21 * pkin(3) + pkin(6);
t13 = (mrSges(6,1) * t25 + mrSges(6,2) * t27) * qJD(5);
t4 = [0.2e1 * m(5) * (t8 * t6 - t41) + 0.2e1 * m(6) * (-t1 * t30 + t3 * t2 - t41); m(5) * (t10 * t8 - t12 * t6 + t32) + m(6) * (-t31 * t10 + t29 * t12 + t32); 0.2e1 * m(5) * (-t12 * t10 + t40) + 0.2e1 * m(6) * (-t12 * t33 + t40); t6 * mrSges(5,2) + t7 * t13 + t37 * t5 + (-mrSges(4,1) * t28 + mrSges(4,2) * t26) * t35 + m(6) * (t29 * t17 - t18 * t5) + (-t21 * t6 + t23 * t5) * t42 + t29 * mrSges(6,3); t11 * t13 - t37 * t9 + (-t26 * mrSges(4,1) - t28 * mrSges(4,2)) * qJD(3) - (t36 * mrSges(6,3) - mrSges(5,2)) * t10 + m(6) * (-t17 * t33 + t18 * t9) + (-t10 * t21 - t23 * t9) * t42; 0.2e1 * qJD(5) * t20 * Ifges(6,4) + 0.2e1 * t18 * t13 + 0.2e1 * (-Ifges(6,4) * t25 + (Ifges(6,1) - Ifges(6,2)) * t27) * t34; m(6) * (t31 * qJD(5) + t25 * t1 + t27 * t2); 0; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2); (t10 * t27 + t12 * t34) * mrSges(6,2) + (-qJD(5) * t12 * t27 + t10 * t25) * mrSges(6,1); (Ifges(6,5) * t27 - Ifges(6,6) * t25 + t15 * t17) * qJD(5); -t13; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t4(1), t4(2), t4(4), t4(7), t4(11); t4(2), t4(3), t4(5), t4(8), t4(12); t4(4), t4(5), t4(6), t4(9), t4(13); t4(7), t4(8), t4(9), t4(10), t4(14); t4(11), t4(12), t4(13), t4(14), t4(15);];
Mq = res;
