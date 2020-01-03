% Calculate time derivative of joint inertia matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRPPR3_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR3_inertiaDJ_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR3_inertiaDJ_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR3_inertiaDJ_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:26:28
% EndTime: 2019-12-31 19:26:29
% DurationCPUTime: 0.24s
% Computational Cost: add. (182->59), mult. (498->85), div. (0->0), fcn. (277->6), ass. (0->39)
t19 = sin(pkin(8));
t24 = cos(qJ(2));
t37 = pkin(1) * qJD(2);
t20 = cos(pkin(8));
t22 = sin(qJ(2));
t39 = t20 * t22;
t4 = (t19 * t24 + t39) * t37;
t48 = 0.2e1 * t4;
t23 = cos(qJ(5));
t47 = t23 ^ 2;
t21 = sin(qJ(5));
t9 = t21 * mrSges(6,1) + t23 * mrSges(6,2);
t46 = mrSges(5,3) + t9;
t45 = -mrSges(4,1) + mrSges(5,2);
t30 = mrSges(6,1) * t23 - mrSges(6,2) * t21;
t7 = t30 * qJD(5);
t44 = 0.2e1 * t7;
t13 = t19 * t22 * pkin(1);
t5 = t20 * t24 * t37 - qJD(2) * t13;
t1 = qJD(4) + t5;
t16 = t24 * pkin(1) + pkin(2);
t27 = pkin(1) * t39 + t19 * t16;
t3 = qJ(4) + t27;
t43 = t3 * t1;
t42 = t5 * mrSges(4,2);
t38 = t21 ^ 2 + t47;
t35 = -t20 * pkin(2) - pkin(3);
t34 = t38 * t4;
t33 = t38 * mrSges(6,3);
t32 = t20 * t16 - t13;
t31 = -pkin(3) - t32;
t29 = -Ifges(6,5) * t21 - Ifges(6,6) * t23;
t15 = t19 * pkin(2) + qJ(4);
t28 = qJD(4) * t3 + t15 * t1;
t26 = (-mrSges(3,1) * t22 - mrSges(3,2) * t24) * t37;
t25 = (-0.2e1 * t47 * Ifges(6,4) + (0.2e1 * Ifges(6,4) * t21 + 0.2e1 * (-Ifges(6,1) + Ifges(6,2)) * t23) * t21) * qJD(5);
t14 = -pkin(7) + t35;
t2 = -pkin(7) + t31;
t6 = [-0.2e1 * t42 + t3 * t44 + t45 * t48 + 0.2e1 * m(6) * (t2 * t34 + t43) + 0.2e1 * m(4) * (t27 * t5 - t32 * t4) + 0.2e1 * m(5) * (t31 * t4 + t43) + t25 + 0.2e1 * t46 * t1 + 0.2e1 * t26 - 0.2e1 * t33 * t4; -t42 + (t3 + t15) * t7 + t26 + (-t33 + t45) * t4 + m(6) * (t14 * t34 + t28) + m(4) * (t19 * t5 - t20 * t4) * pkin(2) + m(5) * (t35 * t4 + t28) + t25 + t46 * (t1 + qJD(4)); t15 * t44 + 0.2e1 * ((m(5) + m(6)) * t15 + t46) * qJD(4) + t25; 0; 0; 0; (m(6) * t38 / 0.2e1 + m(5) / 0.2e1) * t48; 0; 0; 0; t30 * t4 + (-t2 * t9 + t29) * qJD(5); (-t14 * t9 + t29) * qJD(5); -t7; -t9 * qJD(5); 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t6(1), t6(2), t6(4), t6(7), t6(11); t6(2), t6(3), t6(5), t6(8), t6(12); t6(4), t6(5), t6(6), t6(9), t6(13); t6(7), t6(8), t6(9), t6(10), t6(14); t6(11), t6(12), t6(13), t6(14), t6(15);];
Mq = res;
