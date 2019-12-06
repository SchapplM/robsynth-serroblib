% Calculate time derivative of joint inertia matrix for
% S5PRPRP1
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
% Datum: 2019-12-05 15:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5PRPRP1_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP1_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRP1_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRP1_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:28:19
% EndTime: 2019-12-05 15:28:21
% DurationCPUTime: 0.28s
% Computational Cost: add. (307->60), mult. (783->94), div. (0->0), fcn. (623->4), ass. (0->29)
t44 = (-mrSges(5,3) - mrSges(6,2));
t24 = cos(pkin(8));
t33 = -t24 * pkin(3) - pkin(2);
t43 = 0.2e1 * t33;
t42 = 2 * t44;
t41 = m(5) + m(6);
t23 = sin(pkin(8));
t25 = sin(qJ(4));
t26 = cos(qJ(4));
t16 = t25 * t23 - t26 * t24;
t13 = t16 * qJD(4);
t11 = t13 * mrSges(5,2);
t17 = t26 * t23 + t25 * t24;
t14 = t17 * qJD(4);
t12 = t14 * mrSges(6,1);
t3 = pkin(4) * t14 + t13 * qJ(5) - t17 * qJD(5);
t40 = -m(6) * t3 - t14 * mrSges(5,1) - t13 * mrSges(6,3) + t11 - t12;
t39 = m(6) * qJ(5) + mrSges(6,3);
t36 = pkin(6) + qJ(3);
t18 = t36 * t23;
t19 = t36 * t24;
t10 = -t25 * t18 + t26 * t19;
t30 = -t26 * t18 - t25 * t19;
t4 = -t16 * qJD(3) + t30 * qJD(4);
t5 = t17 * qJD(3) + t10 * qJD(4);
t34 = t10 * t4 - t30 * t5;
t31 = 2 * Ifges(6,5) - 2 * Ifges(5,4);
t7 = t16 * pkin(4) - t17 * qJ(5) + t33;
t1 = [0.2e1 * t41 * (-t17 * t13 + t16 * t14); t41 * (-t10 * t13 - t14 * t30 + t5 * t16 + t4 * t17); -t11 * t43 + 0.2e1 * t3 * (t16 * mrSges(6,1) - t17 * mrSges(6,3)) + 0.2e1 * t7 * t12 + 0.2e1 * m(6) * (t7 * t3 + t34) + 0.2e1 * m(5) * t34 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3) * (t23 ^ 2 + t24 ^ 2) + (mrSges(5,1) * t43 + t17 * t31 + 0.2e1 * (Ifges(6,3) + Ifges(5,2)) * t16 + t10 * t42) * t14 - (-0.2e1 * t7 * mrSges(6,3) + t16 * t31 + t30 * t42 + 0.2e1 * (Ifges(5,1) + Ifges(6,1)) * t17) * t13 - 0.2e1 * t44 * (-t4 * t16 + t5 * t17); 0; -t40; 0; t40; m(6) * qJD(5) * t10 + (-Ifges(5,6) + Ifges(6,6)) * t14 - (Ifges(5,5) + Ifges(6,4)) * t13 + (pkin(4) * t13 - qJ(5) * t14 - qJD(5) * t16) * mrSges(6,2) + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t5 + (-mrSges(5,2) + t39) * t4; 0; 0.2e1 * t39 * qJD(5); m(6) * t14; m(6) * t5 - t13 * mrSges(6,2); 0; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
