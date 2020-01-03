% Calculate time derivative of joint inertia matrix for
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4PRRR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:31
% DurationCPUTime: 0.23s
% Computational Cost: add. (266->58), mult. (680->107), div. (0->0), fcn. (512->4), ass. (0->29)
t33 = qJD(3) + qJD(4);
t32 = 2 * m(5);
t23 = cos(qJ(3));
t31 = -0.2e1 * t23 * pkin(3) - (2 * pkin(2));
t30 = m(5) * pkin(3);
t29 = -pkin(6) - pkin(5);
t20 = sin(qJ(4));
t21 = sin(qJ(3));
t22 = cos(qJ(4));
t13 = -t20 * t21 + t22 * t23;
t14 = t20 * t23 + t22 * t21;
t9 = t33 * t14;
t28 = t13 * t9;
t8 = t33 * t13;
t27 = t14 * t8;
t26 = 0.2e1 * t23;
t1 = t9 * mrSges(5,1) + t8 * mrSges(5,2);
t25 = qJD(3) * t29;
t17 = t29 * t21;
t18 = t29 * t23;
t10 = t22 * t17 + t20 * t18;
t15 = t21 * t25;
t16 = t23 * t25;
t3 = qJD(4) * t10 + t22 * t15 + t20 * t16;
t11 = t20 * t17 - t22 * t18;
t4 = -qJD(4) * t11 - t20 * t15 + t22 * t16;
t24 = t4 * mrSges(5,1) - t3 * mrSges(5,2) + Ifges(5,5) * t8 - Ifges(5,6) * t9;
t12 = (-mrSges(5,1) * t20 - mrSges(5,2) * t22) * qJD(4) * pkin(3);
t2 = [(t27 - t28) * t32; m(5) * (-t10 * t9 + t11 * t8 + t4 * t13 + t3 * t14); 0.2e1 * Ifges(5,1) * t27 - 0.2e1 * Ifges(5,2) * t28 + t1 * t31 + (t10 * t4 + t11 * t3) * t32 + ((-(pkin(2) * mrSges(4,2)) + Ifges(4,4) * t23) * t26 + (0.2e1 * pkin(3) * (-t13 * mrSges(5,1) + t14 * mrSges(5,2)) - (2 * pkin(2) * mrSges(4,1)) + t30 * t31 - 0.2e1 * Ifges(4,4) * t21 + (Ifges(4,1) - Ifges(4,2)) * t26) * t21) * qJD(3) + 0.2e1 * (t8 * t13 - t14 * t9) * Ifges(5,4) + 0.2e1 * (-t10 * t8 - t11 * t9 + t3 * t13 - t4 * t14) * mrSges(5,3); (-t21 * mrSges(4,1) - t23 * mrSges(4,2)) * qJD(3) + (t20 * t8 - t22 * t9 + (-t13 * t20 + t14 * t22) * qJD(4)) * t30 - t1; (Ifges(4,5) * t23 - Ifges(4,6) * t21 + (-mrSges(4,1) * t23 + mrSges(4,2) * t21) * pkin(5)) * qJD(3) + (m(5) * (t20 * t3 + t22 * t4 + (-t10 * t20 + t11 * t22) * qJD(4)) + (-t20 * t9 - t22 * t8 + (t13 * t22 + t14 * t20) * qJD(4)) * mrSges(5,3)) * pkin(3) + t24; 0.2e1 * t12; -t1; t24; t12; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t2(1), t2(2), t2(4), t2(7); t2(2), t2(3), t2(5), t2(8); t2(4), t2(5), t2(6), t2(9); t2(7), t2(8), t2(9), t2(10);];
Mq = res;
