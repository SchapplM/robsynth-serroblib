% Calculate time derivative of joint inertia matrix for
% S5PPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PPRRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:06:40
% EndTime: 2019-12-05 15:06:40
% DurationCPUTime: 0.23s
% Computational Cost: add. (169->69), mult. (490->108), div. (0->0), fcn. (338->6), ass. (0->35)
t39 = -2 * mrSges(6,3);
t38 = m(6) * pkin(4);
t18 = sin(pkin(8));
t19 = cos(pkin(8));
t21 = sin(qJ(3));
t23 = cos(qJ(3));
t7 = t23 * t18 + t21 * t19;
t5 = t7 * qJD(3);
t6 = t21 * t18 - t23 * t19;
t37 = t6 * t5;
t22 = cos(qJ(4));
t36 = mrSges(5,2) * t22;
t32 = -qJ(5) - pkin(6);
t12 = t32 * t22;
t35 = t12 * t22;
t34 = mrSges(5,2) + mrSges(6,2);
t33 = Ifges(5,4) + Ifges(6,4);
t29 = qJD(4) * t22;
t20 = sin(qJ(4));
t30 = qJD(4) * t20;
t8 = mrSges(6,1) * t30 + mrSges(6,2) * t29;
t31 = t20 ^ 2 + t22 ^ 2;
t28 = 0.2e1 * t20;
t27 = mrSges(6,1) + t38;
t4 = t6 * qJD(3);
t26 = t31 * t4;
t25 = qJD(4) * t32;
t24 = -mrSges(5,1) - t27;
t13 = -t22 * pkin(4) - pkin(3);
t11 = -t22 * mrSges(6,1) + t20 * mrSges(6,2);
t10 = t32 * t20;
t9 = (mrSges(5,1) * t20 + t36) * qJD(4);
t3 = -t20 * qJD(5) + t22 * t25;
t2 = t22 * qJD(5) + t20 * t25;
t1 = [0.2e1 * m(4) * (-t7 * t4 + t37) + 0.4e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * (-t7 * t26 + t37); 0; 0; (t9 + t8) * t6 + (-t22 * mrSges(5,1) + t20 * mrSges(5,2) - mrSges(4,1) + t11) * t5 - (-mrSges(4,2) + (mrSges(5,3) + mrSges(6,3)) * t31) * t4 + m(5) * (-pkin(3) * t5 - pkin(6) * t26) + (t13 * t5 + t4 * t35 + (pkin(4) * qJD(4) * t6 + t10 * t4) * t20 + (-t10 * t29 + t2 * t22 + (qJD(4) * t12 - t3) * t20) * t7) * m(6); m(6) * (t2 * t20 + t3 * t22 + (-t10 * t20 - t35) * qJD(4)); -0.2e1 * pkin(3) * t9 + 0.2e1 * t13 * t8 + 0.2e1 * m(6) * (t10 * t3 - t12 * t2) + (t3 * t39 + (-t12 * t39 - t33 * t28 + 0.2e1 * (m(6) * t13 + t11) * pkin(4)) * qJD(4)) * t20 + (0.2e1 * t2 * mrSges(6,3) + (t10 * t39 + 0.2e1 * t33 * t22 + (Ifges(5,1) + Ifges(6,1) - Ifges(5,2) - Ifges(6,2)) * t28) * qJD(4)) * t22; -(t24 * t20 - t34 * t22) * t4 + (t34 * t20 + t24 * t22) * t7 * qJD(4); (-t36 + (-mrSges(5,1) - t38) * t20) * qJD(4) - t8; -t2 * mrSges(6,2) + t27 * t3 + ((mrSges(5,2) * pkin(6) - Ifges(5,6) - Ifges(6,6)) * t20 + (-mrSges(5,1) * pkin(6) - mrSges(6,3) * pkin(4) + Ifges(5,5) + Ifges(6,5)) * t22) * qJD(4); 0; m(6) * t5; 0; t30 * t38 + t8; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
