% Calculate time derivative of joint inertia matrix for
% S4RRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRPR4_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR4_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR4_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:27
% EndTime: 2019-12-31 17:02:27
% DurationCPUTime: 0.29s
% Computational Cost: add. (393->80), mult. (918->118), div. (0->0), fcn. (696->6), ass. (0->45)
t59 = 2 * mrSges(4,3);
t58 = 2 * mrSges(5,3);
t39 = cos(qJ(2));
t48 = pkin(1) * qJD(2);
t45 = t39 * t48;
t28 = qJD(3) + t45;
t34 = sin(pkin(7));
t35 = cos(pkin(7));
t49 = t34 ^ 2 + t35 ^ 2;
t43 = t49 * t28;
t36 = sin(qJ(4));
t38 = cos(qJ(4));
t22 = -t34 * t36 + t35 * t38;
t23 = t34 * t38 + t35 * t36;
t57 = -mrSges(4,1) * t35 - mrSges(5,1) * t22 + mrSges(4,2) * t34 + mrSges(5,2) * t23;
t40 = t49 * qJD(3);
t56 = 2 * m(4);
t55 = 2 * m(5);
t18 = t22 * qJD(4);
t19 = t23 * qJD(4);
t11 = t19 * mrSges(5,1) + t18 * mrSges(5,2);
t52 = 0.2e1 * t11;
t51 = t39 * pkin(1);
t50 = Ifges(5,5) * t18 - Ifges(5,6) * t19;
t37 = sin(qJ(2));
t46 = t37 * t48;
t44 = 0.2e1 * t23 * t18 * Ifges(5,1) - 0.2e1 * t22 * Ifges(5,2) * t19 + 0.2e1 * (t18 * t22 - t19 * t23) * Ifges(5,4);
t30 = -pkin(3) * t35 - pkin(2);
t42 = t49 * qJ(3);
t29 = pkin(1) * t37 + qJ(3);
t20 = (-pkin(6) - t29) * t34;
t31 = t35 * pkin(6);
t21 = t29 * t35 + t31;
t9 = t20 * t38 - t21 * t36;
t10 = t20 * t36 + t21 * t38;
t25 = (-pkin(6) - qJ(3)) * t34;
t27 = qJ(3) * t35 + t31;
t13 = t25 * t38 - t27 * t36;
t14 = t25 * t36 + t27 * t38;
t24 = t30 - t51;
t8 = -qJD(3) * t23 - qJD(4) * t14;
t7 = qJD(3) * t22 + qJD(4) * t13;
t2 = -qJD(4) * t10 - t23 * t28;
t1 = qJD(4) * t9 + t22 * t28;
t3 = [(t1 * t10 + t2 * t9) * t55 + t24 * t52 - 0.2e1 * mrSges(3,2) * t45 + t44 - 0.2e1 * (t9 * t18 + t2 * t23) * mrSges(5,3) + (t1 * t22 - t10 * t19) * t58 + (t29 * t56 + t59) * t43 + (t24 * t55 + (-pkin(2) - t51) * t56 - (2 * mrSges(3,1)) + 0.2e1 * t57) * t46; (t24 + t30) * t11 + (-mrSges(3,2) * t39 + (-mrSges(3,1) + t57) * t37) * t48 + m(5) * (t1 * t14 + t10 * t7 + t13 * t2 + t30 * t46 + t8 * t9) + m(4) * (-pkin(2) * t46 + t28 * t42 + t29 * t40) + (t43 + t40) * mrSges(4,3) + ((-t2 - t8) * t23 + (t1 + t7) * t22 - (t10 + t14) * t19 + (-t13 - t9) * t18) * mrSges(5,3) + t44; t30 * t52 + (t13 * t8 + t14 * t7) * t55 + t42 * t56 * qJD(3) + t44 + t40 * t59 + (-t13 * t18 - t14 * t19 + t7 * t22 - t8 * t23) * t58; (m(4) + m(5)) * t46 + t11; t11; 0; mrSges(5,1) * t2 - mrSges(5,2) * t1 + t50; mrSges(5,1) * t8 - mrSges(5,2) * t7 + t50; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t3(1), t3(2), t3(4), t3(7); t3(2), t3(3), t3(5), t3(8); t3(4), t3(5), t3(6), t3(9); t3(7), t3(8), t3(9), t3(10);];
Mq = res;
