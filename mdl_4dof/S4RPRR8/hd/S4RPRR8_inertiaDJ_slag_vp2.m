% Calculate time derivative of joint inertia matrix for
% S4RPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2019-12-31 16:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S4RPRR8_inertiaDJ_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR8_inertiaDJ_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR8_inertiaDJ_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR8_inertiaDJ_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_time_derivative_floatb_twist_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:55:06
% EndTime: 2019-12-31 16:55:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (354->61), mult. (729->107), div. (0->0), fcn. (550->4), ass. (0->36)
t27 = cos(qJ(3));
t48 = qJD(3) * t27;
t25 = sin(qJ(3));
t47 = -t25 * mrSges(4,1) - t27 * mrSges(4,2);
t46 = qJD(3) + qJD(4);
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t14 = -t24 * t27 - t26 * t25;
t38 = t26 * t27;
t41 = t24 * t25;
t15 = t38 - t41;
t36 = qJD(3) * t25;
t7 = -qJD(4) * t41 - t24 * t36 + t46 * t38;
t8 = t46 * t14;
t45 = -t24 * t7 - t26 * t8 + (t14 * t26 + t15 * t24) * qJD(4);
t44 = t14 * t7;
t43 = t15 * t8;
t28 = -pkin(1) - pkin(5);
t42 = pkin(6) - t28;
t35 = 2 * mrSges(5,3);
t34 = t8 * mrSges(5,1) - t7 * mrSges(5,2);
t17 = t42 * t27;
t33 = -t43 + t44;
t12 = t42 * t36;
t13 = qJD(3) * t17;
t16 = t42 * t25;
t9 = t24 * t16 - t26 * t17;
t2 = qJD(4) * t9 + t24 * t12 - t26 * t13;
t10 = -t26 * t16 - t24 * t17;
t3 = -qJD(4) * t10 + t26 * t12 + t24 * t13;
t31 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t8 - Ifges(5,6) * t7;
t29 = t7 * t10 - t14 * t2 + t15 * t3 + t8 * t9;
t21 = t25 * pkin(3) + qJ(2);
t18 = pkin(3) * t48 + qJD(2);
t11 = (-mrSges(5,1) * t24 - mrSges(5,2) * t26) * qJD(4) * pkin(3);
t1 = [0.2e1 * t18 * (-t14 * mrSges(5,1) + t15 * mrSges(5,2)) + 0.2e1 * t21 * (t7 * mrSges(5,1) + t8 * mrSges(5,2)) + 0.2e1 * Ifges(5,1) * t43 - 0.2e1 * Ifges(5,2) * t44 + 0.2e1 * m(5) * (t10 * t2 + t21 * t18 + t9 * t3) + 0.2e1 * (t8 * t14 - t15 * t7) * Ifges(5,4) - t29 * t35 + 0.2e1 * (qJ(2) * (mrSges(4,1) * t27 - mrSges(4,2) * t25) + (t25 ^ 2 - t27 ^ 2) * Ifges(4,4)) * qJD(3) + 0.2e1 * (mrSges(3,3) + (m(3) + m(4)) * qJ(2) - t47) * qJD(2) + 0.2e1 * (-Ifges(4,1) + Ifges(4,2)) * t25 * t48; m(5) * t29 + t33 * t35; -0.2e1 * m(5) * t33; ((-mrSges(4,2) * t28 - Ifges(4,6)) * t27 + (-mrSges(4,1) * t28 - Ifges(4,5)) * t25) * qJD(3) + (m(5) * (t2 * t24 + t26 * t3 + (t10 * t26 - t24 * t9) * qJD(4)) + t45 * mrSges(5,3)) * pkin(3) + t31; -m(5) * t45 * pkin(3) + t47 * qJD(3) + t34; 0.2e1 * t11; t31; t34; t11; 0;];
%% Postprocessing: Reshape Output
% From vec2symmat_4_matlab.m
res = [t1(1), t1(2), t1(4), t1(7); t1(2), t1(3), t1(5), t1(8); t1(4), t1(5), t1(6), t1(9); t1(7), t1(8), t1(9), t1(10);];
Mq = res;
