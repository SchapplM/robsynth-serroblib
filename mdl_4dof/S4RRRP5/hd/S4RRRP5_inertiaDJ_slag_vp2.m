% Calculate time derivative of joint inertia matrix for
% S4RRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
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
% Datum: 2019-12-31 17:17
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RRRP5_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP5_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRP5_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRP5_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:16:38
% EndTime: 2019-12-31 17:16:40
% DurationCPUTime: 0.49s
% Computational Cost: add. (476->92), mult. (1161->149), div. (0->0), fcn. (848->4), ass. (0->41)
t47 = (mrSges(5,2) + mrSges(4,3));
t52 = 2 * t47;
t51 = -mrSges(4,1) - mrSges(5,1);
t50 = qJD(2) + qJD(3);
t32 = cos(qJ(2));
t27 = -t32 * pkin(2) - pkin(1);
t49 = 0.2e1 * t27;
t48 = -pkin(6) - pkin(5);
t29 = sin(qJ(3));
t46 = qJD(3) * t29;
t31 = cos(qJ(3));
t45 = qJD(3) * t31;
t44 = 0.2e1 * t32;
t43 = pkin(2) * t46;
t42 = pkin(2) * t45;
t30 = sin(qJ(2));
t21 = t48 * t30;
t22 = t48 * t32;
t36 = t31 * t21 + t29 * t22;
t41 = t36 * t46;
t15 = t29 * t21 - t31 * t22;
t39 = qJD(2) * t48;
t20 = t30 * t39;
t37 = t32 * t39;
t5 = qJD(3) * t36 + t31 * t20 + t29 * t37;
t6 = qJD(3) * t15 + t29 * t20 - t31 * t37;
t40 = t15 * t5 - t36 * t6;
t38 = -2 * Ifges(4,4) + 2 * Ifges(5,5);
t18 = t29 * t32 + t31 * t30;
t17 = t29 * t30 - t31 * t32;
t12 = t50 * t17;
t13 = t50 * t18;
t34 = t51 * t6 + (-mrSges(4,2) + mrSges(5,3)) * t5 + (-Ifges(4,6) + Ifges(5,6)) * t13 - (Ifges(5,4) + Ifges(4,5)) * t12;
t23 = qJD(4) + t42;
t33 = -mrSges(4,2) * t42 + t23 * mrSges(5,3) + t51 * t43;
t28 = qJD(4) * mrSges(5,3);
t26 = -t31 * pkin(2) - pkin(3);
t24 = t29 * pkin(2) + qJ(4);
t7 = t17 * pkin(3) - t18 * qJ(4) + t27;
t2 = qJD(2) * t30 * pkin(2) + t13 * pkin(3) + t12 * qJ(4) - t18 * qJD(4);
t1 = [(t13 * mrSges(4,1) - t12 * mrSges(4,2)) * t49 + 0.2e1 * t7 * (t13 * mrSges(5,1) + t12 * mrSges(5,3)) + 0.2e1 * m(4) * t40 + 0.2e1 * m(5) * (t7 * t2 + t40) + (-0.2e1 * t2 * mrSges(5,3) + t13 * t38 + t6 * t52 - 0.2e1 * (Ifges(4,1) + Ifges(5,1)) * t12) * t18 + (0.2e1 * t2 * mrSges(5,1) - t12 * t38 - 0.2e1 * t47 * t5 + 0.2e1 * (Ifges(4,2) + Ifges(5,3)) * t13) * t17 + ((-pkin(1) * mrSges(3,2) + Ifges(3,4) * t32) * t44 + (-0.2e1 * pkin(1) * mrSges(3,1) + m(4) * pkin(2) * t49 + 0.2e1 * pkin(2) * (t17 * mrSges(4,1) + t18 * mrSges(4,2)) - 0.2e1 * Ifges(3,4) * t30 + (Ifges(3,1) - Ifges(3,2)) * t44) * t30) * qJD(2) + (t12 * t36 - t15 * t13) * t52; m(5) * (t23 * t15 + t24 * t5 + t26 * t6) + (-t26 * t12 - t24 * t13 - t23 * t17) * mrSges(5,2) + (Ifges(3,5) * t32 - Ifges(3,6) * t30 + (-mrSges(3,1) * t32 + mrSges(3,2) * t30) * pkin(5)) * qJD(2) + (t18 * mrSges(5,2) * t46 - m(5) * t41 + m(4) * (t15 * t45 + t29 * t5 - t31 * t6 - t41) + (t31 * t12 - t29 * t13 + (-t17 * t31 + t18 * t29) * qJD(3)) * mrSges(4,3)) * pkin(2) + t34; 0.2e1 * m(5) * (t24 * t23 + t26 * t43) + 0.2e1 * t33; m(5) * (-pkin(3) * t6 + qJ(4) * t5 + qJD(4) * t15) + (pkin(3) * t12 - qJ(4) * t13 - qJD(4) * t17) * mrSges(5,2) + t34; t28 + m(5) * (-pkin(3) * t43 + qJ(4) * t23 + qJD(4) * t24) + t33; 0.2e1 * m(5) * qJ(4) * qJD(4) + 0.2e1 * t28; m(5) * t6 - t12 * mrSges(5,2); m(5) * t43; 0; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
