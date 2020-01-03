% Calculate time derivative of joint inertia matrix for
% S4RRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR8_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR8_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR8_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:07:49
% EndTime: 2019-12-31 17:07:50
% DurationCPUTime: 0.40s
% Computational Cost: add. (338->93), mult. (772->147), div. (0->0), fcn. (530->4), ass. (0->42)
t22 = sin(qJ(2));
t41 = pkin(5) - pkin(6);
t19 = t41 * t22;
t21 = sin(qJ(4));
t23 = cos(qJ(4));
t49 = -t21 * mrSges(5,1) - t23 * mrSges(5,2);
t48 = qJD(2) - qJD(4);
t24 = cos(qJ(2));
t47 = (m(4) * pkin(5) + mrSges(4,2)) * t24;
t46 = 2 * m(5);
t45 = -2 * pkin(1);
t34 = t22 * qJ(3);
t29 = -t24 * pkin(2) - t34;
t18 = -pkin(1) + t29;
t44 = 0.2e1 * t18;
t42 = pkin(2) + pkin(3);
t16 = -t21 * qJ(3) - t23 * t42;
t8 = t23 * qJD(3) + qJD(4) * t16;
t40 = t8 * mrSges(5,2);
t17 = t23 * qJ(3) - t21 * t42;
t9 = -t21 * qJD(3) - qJD(4) * t17;
t39 = t9 * mrSges(5,1);
t36 = Ifges(4,5) - Ifges(3,4);
t35 = qJ(3) * t24;
t32 = t22 * qJD(3);
t31 = 0.2e1 * t24;
t20 = t41 * t24;
t30 = -t24 * mrSges(4,1) - t22 * mrSges(4,3);
t5 = t23 * t19 - t21 * t20;
t6 = t21 * t19 + t23 * t20;
t28 = t24 * t21 - t22 * t23;
t27 = t22 * t21 + t24 * t23;
t14 = qJD(2) * t19;
t15 = qJD(2) * t20;
t1 = qJD(4) * t5 - t23 * t14 + t21 * t15;
t2 = -qJD(4) * t6 + t21 * t14 + t23 * t15;
t3 = t48 * t28;
t4 = t48 * t27;
t26 = t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t4 - Ifges(5,6) * t3;
t11 = t42 * t24 + pkin(1) + t34;
t7 = t32 + (-t42 * t22 + t35) * qJD(2);
t10 = [0.2e1 * t3 * Ifges(5,2) * t27 - 0.2e1 * t4 * t28 * Ifges(5,1) + 0.2e1 * t7 * (mrSges(5,1) * t27 - mrSges(5,2) * t28) + 0.2e1 * t11 * (t3 * mrSges(5,1) + t4 * mrSges(5,2)) + (t6 * t1 + t11 * t7 + t5 * t2) * t46 + (((mrSges(3,2) * t45) - 0.2e1 * t18 * mrSges(4,3) - t36 * t31) * t24 + ((mrSges(3,1) * t45) + mrSges(4,1) * t44 + (Ifges(3,1) + Ifges(4,1) - Ifges(3,2) - Ifges(4,3)) * t31) * t22) * qJD(2) + 0.2e1 * (-t27 * t4 + t28 * t3) * Ifges(5,4) + 0.2e1 * (-t1 * t27 + t2 * t28 - t6 * t3 - t5 * t4) * mrSges(5,3) + 0.2e1 * t36 * t22 ^ 2 * qJD(2) + (m(4) * t44 + 0.2e1 * t30) * (-t32 + (pkin(2) * t22 - t35) * qJD(2)); m(5) * (t17 * t1 + t16 * t2 + t9 * t5 + t8 * t6) + (-t16 * t4 - t17 * t3 - t27 * t8 + t28 * t9) * mrSges(5,3) + ((-pkin(2) * mrSges(4,2) + Ifges(4,4) + Ifges(3,5)) * t24 + (-qJ(3) * mrSges(4,2) - Ifges(3,6) + Ifges(4,6)) * t22 + (m(4) * t29 - t24 * mrSges(3,1) + t22 * mrSges(3,2) + t30) * pkin(5)) * qJD(2) - t26 + qJD(3) * t47; (t16 * t9 + t17 * t8) * t46 + 0.2e1 * t40 - 0.2e1 * t39 + 0.2e1 * (m(4) * qJ(3) + mrSges(4,3)) * qJD(3); m(5) * (t21 * t1 + t23 * t2 + (-t21 * t5 + t23 * t6) * qJD(4)) + qJD(2) * t47 + (-t21 * t3 - t23 * t4 + (-t21 * t28 - t23 * t27) * qJD(4)) * mrSges(5,3); m(5) * (t21 * t8 + t23 * t9) + (m(5) * (-t16 * t21 + t17 * t23) - t49) * qJD(4); 0; t26; t39 - t40; t49 * qJD(4); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t10(1), t10(2), t10(4), t10(7); t10(2), t10(3), t10(5), t10(8); t10(4), t10(5), t10(6), t10(9); t10(7), t10(8), t10(9), t10(10);];
Mq = res;
