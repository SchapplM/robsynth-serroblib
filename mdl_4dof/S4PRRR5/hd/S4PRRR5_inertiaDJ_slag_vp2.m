% Calculate time derivative of joint inertia matrix for
% S4PRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:33:35
% EndTime: 2019-12-31 16:33:36
% DurationCPUTime: 0.29s
% Computational Cost: add. (190->51), mult. (545->92), div. (0->0), fcn. (383->6), ass. (0->33)
t22 = sin(qJ(4));
t25 = cos(qJ(4));
t36 = t22 ^ 2 + t25 ^ 2;
t23 = sin(qJ(3));
t24 = sin(qJ(2));
t26 = cos(qJ(3));
t27 = cos(qJ(2));
t11 = t23 * t24 - t26 * t27;
t51 = qJD(2) + qJD(3);
t5 = t51 * t11;
t32 = t36 * t5;
t14 = -t25 * mrSges(5,1) + t22 * mrSges(5,2);
t55 = -mrSges(4,1) + t14;
t54 = Ifges(5,1) - Ifges(5,2);
t52 = t36 * t26;
t50 = (mrSges(5,3) * t36 - mrSges(4,2)) * t26 + t55 * t23;
t49 = m(4) / 0.2e1;
t12 = t23 * t27 + t26 * t24;
t6 = t51 * t12;
t48 = t11 * t6;
t43 = Ifges(5,6) * t22;
t42 = t11 * t23;
t35 = pkin(2) * qJD(3);
t34 = qJD(4) * t22;
t33 = qJD(4) * t25;
t30 = mrSges(5,1) * t22 + mrSges(5,2) * t25;
t29 = (-0.2e1 * Ifges(5,4) * t22 + t54 * t25) * t34 + (0.2e1 * Ifges(5,4) * t25 + t54 * t22) * t33;
t13 = t30 * qJD(4);
t28 = t5 * mrSges(4,2) - mrSges(5,3) * t32 + t11 * t13 + t55 * t6;
t19 = Ifges(5,5) * t33;
t18 = -t26 * pkin(2) - pkin(3);
t17 = t23 * pkin(2) + pkin(6);
t1 = [0.2e1 * m(4) * (-t12 * t5 + t48) + 0.2e1 * m(5) * (-t12 * t32 + t48); m(5) * (-t17 * t32 + t18 * t6) + (-t24 * mrSges(3,1) - t27 * mrSges(3,2)) * qJD(2) + 0.2e1 * ((-t23 * t5 - t26 * t6) * t49 + ((t12 * t26 + t42) * t49 + m(5) * (t52 * t12 + t42) / 0.2e1) * qJD(3)) * pkin(2) + t28; 0.2e1 * t18 * t13 + 0.2e1 * (m(5) * (t52 * t17 + t18 * t23) + t50) * t35 + t29; m(5) * (-pkin(3) * t6 - pkin(6) * t32) + t28; (t18 - pkin(3)) * t13 + (m(5) * (-pkin(3) * t23 + t52 * pkin(6)) + t50) * t35 + t29; -0.2e1 * pkin(3) * t13 + t29; (t12 * t34 + t25 * t5) * mrSges(5,2) + (-t12 * t33 + t22 * t5) * mrSges(5,1); t19 - t30 * t26 * t35 + (t14 * t17 - t43) * qJD(4); t19 + (t14 * pkin(6) - t43) * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
