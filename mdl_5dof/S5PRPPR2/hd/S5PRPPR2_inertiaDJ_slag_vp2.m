% Calculate time derivative of joint inertia matrix for
% S5PRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3,theta4]';
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
% Datum: 2019-12-05 15:25
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPPR2_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR2_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR2_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR2_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:23:59
% EndTime: 2019-12-05 15:24:01
% DurationCPUTime: 0.32s
% Computational Cost: add. (353->70), mult. (926->123), div. (0->0), fcn. (839->8), ass. (0->43)
t25 = sin(pkin(9));
t27 = cos(pkin(9));
t37 = t25 ^ 2 + t27 ^ 2;
t46 = mrSges(5,3) * t37;
t26 = sin(pkin(8));
t45 = t26 * pkin(2);
t28 = cos(pkin(8));
t44 = t28 * pkin(2);
t43 = m(5) * qJD(4);
t29 = sin(qJ(5));
t31 = cos(qJ(5));
t33 = t29 * t25 - t31 * t27;
t13 = t33 * qJD(5);
t42 = 2 * m(6);
t22 = qJ(4) + t45;
t41 = pkin(6) + t22;
t30 = sin(qJ(2));
t32 = cos(qJ(2));
t19 = t26 * t32 + t28 * t30;
t11 = t19 * qJD(2);
t17 = t26 * t30 - t28 * t32;
t10 = t17 * t11;
t20 = t31 * t25 + t29 * t27;
t14 = t20 * qJD(5);
t39 = t33 * t14;
t38 = t20 * t13;
t36 = -pkin(3) - t44;
t35 = t37 * t19;
t34 = t37 * t22;
t15 = t41 * t25;
t16 = t41 * t27;
t7 = -t31 * t15 - t29 * t16;
t8 = -t29 * t15 + t31 * t16;
t21 = -t27 * pkin(4) + t36;
t12 = t17 * qJD(2);
t9 = t14 * mrSges(6,1) - t13 * mrSges(6,2);
t6 = t33 * t19;
t5 = t20 * t19;
t4 = -qJD(4) * t20 - qJD(5) * t8;
t3 = -qJD(4) * t33 + qJD(5) * t7;
t2 = t12 * t20 + t19 * t13;
t1 = t12 * t33 - t14 * t19;
t18 = [0.2e1 * m(6) * (-t6 * t1 - t5 * t2 + t10) + 0.2e1 * m(4) * (-t19 * t12 + t10) + 0.2e1 * m(5) * (-t12 * t35 + t10); m(6) * (t8 * t1 + t7 * t2 - t3 * t6 - t4 * t5) + t17 * t9 + t35 * t43 + (-m(4) * t44 + m(5) * t36 + m(6) * t21 - t27 * mrSges(5,1) + mrSges(6,1) * t33 + t25 * mrSges(5,2) + t20 * mrSges(6,2) - mrSges(4,1)) * t11 + (-t30 * mrSges(3,1) - t32 * mrSges(3,2)) * qJD(2) + (-t1 * t33 - t5 * t13 + t6 * t14 - t2 * t20) * mrSges(6,3) + (-m(4) * t45 - m(5) * t34 + mrSges(4,2) - t46) * t12; -0.2e1 * Ifges(6,1) * t38 + 0.2e1 * Ifges(6,2) * t39 + 0.2e1 * t21 * t9 + (t8 * t3 + t7 * t4) * t42 + 0.2e1 * t34 * t43 + 0.2e1 * qJD(4) * t46 + 0.2e1 * (t13 * t33 - t20 * t14) * Ifges(6,4) + 0.2e1 * (t7 * t13 - t8 * t14 - t4 * t20 - t3 * t33) * mrSges(6,3); m(6) * (t20 * t1 + t13 * t6 + t14 * t5 - t2 * t33); m(6) * (-t13 * t8 - t14 * t7 + t20 * t3 - t33 * t4); (-t38 + t39) * t42; 0.2e1 * (m(6) / 0.2e1 + m(5) / 0.2e1) * t11; t9; 0; 0; t2 * mrSges(6,1) - t1 * mrSges(6,2); t4 * mrSges(6,1) - t3 * mrSges(6,2) - Ifges(6,5) * t13 - Ifges(6,6) * t14; -t9; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t18(1), t18(2), t18(4), t18(7), t18(11); t18(2), t18(3), t18(5), t18(8), t18(12); t18(4), t18(5), t18(6), t18(9), t18(13); t18(7), t18(8), t18(9), t18(10), t18(14); t18(11), t18(12), t18(13), t18(14), t18(15);];
Mq = res;
