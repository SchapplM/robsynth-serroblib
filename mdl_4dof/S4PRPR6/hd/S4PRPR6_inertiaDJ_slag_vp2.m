% Calculate time derivative of joint inertia matrix for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MqD [4x4]
%   time derivative of inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRPR6_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:21
% DurationCPUTime: 0.21s
% Computational Cost: add. (173->52), mult. (505->96), div. (0->0), fcn. (407->6), ass. (0->30)
t19 = sin(pkin(7));
t20 = cos(pkin(7));
t30 = t19 ^ 2 + t20 ^ 2;
t27 = m(4) * t30;
t33 = qJD(3) * t27;
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t25 = t21 * t19 - t23 * t20;
t10 = t25 * qJD(4);
t32 = 2 * m(5);
t31 = pkin(5) + qJ(3);
t24 = cos(qJ(2));
t28 = qJD(2) * t24;
t26 = t30 * mrSges(4,3);
t14 = t31 * t19;
t15 = t31 * t20;
t6 = -t23 * t14 - t21 * t15;
t7 = -t21 * t14 + t23 * t15;
t13 = t23 * t19 + t21 * t20;
t11 = t13 * qJD(4);
t22 = sin(qJ(2));
t16 = -t20 * pkin(3) - pkin(2);
t9 = t25 * t22;
t8 = t13 * t22;
t5 = t11 * mrSges(5,1) - t10 * mrSges(5,2);
t4 = t22 * t10 - t13 * t28;
t3 = -t22 * t11 - t25 * t28;
t2 = -t13 * qJD(3) - t7 * qJD(4);
t1 = -t25 * qJD(3) + t6 * qJD(4);
t12 = [(-t9 * t3 - t8 * t4) * t32 + 0.4e1 * (m(4) * (-0.1e1 + t30) / 0.2e1 - m(5) / 0.2e1) * t22 * t28; -t24 * t5 + m(5) * (-t1 * t9 - t2 * t8 + t7 * t3 + t6 * t4) + t22 * t33 + ((-m(4) * pkin(2) + m(5) * t16 - t20 * mrSges(4,1) + mrSges(5,1) * t25 + t19 * mrSges(4,2) + t13 * mrSges(5,2) - mrSges(3,1)) * t22 + (qJ(3) * t27 - mrSges(3,2) + t26) * t24) * qJD(2) + (-t8 * t10 + t9 * t11 - t4 * t13 - t25 * t3) * mrSges(5,3); (t7 * t1 + t6 * t2) * t32 + 0.2e1 * t16 * t5 - 0.2e1 * t10 * t13 * Ifges(5,1) + 0.2e1 * t11 * Ifges(5,2) * t25 + 0.2e1 * qJ(3) * t33 + 0.2e1 * t26 * qJD(3) + 0.2e1 * (t10 * t25 - t13 * t11) * Ifges(5,4) + 0.2e1 * (-t1 * t25 + t6 * t10 - t7 * t11 - t2 * t13) * mrSges(5,3); (m(4) + m(5)) * t22 * qJD(2); t5; 0; t4 * mrSges(5,1) - t3 * mrSges(5,2); t2 * mrSges(5,1) - t1 * mrSges(5,2) - Ifges(5,5) * t10 - Ifges(5,6) * t11; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t12(1), t12(2), t12(4), t12(7); t12(2), t12(3), t12(5), t12(8); t12(4), t12(5), t12(6), t12(9); t12(7), t12(8), t12(9), t12(10);];
Mq = res;
